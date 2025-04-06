/************************************************************************/
/*                                                                      */
/*   utility to convert and extract images from 360 degree environments */
/*                                                                      */
/*            Copyright 2024 by Kay F. Jahnke                           */
/*                                                                      */
/*    The git repository for this software is at                        */
/*                                                                      */
/*    https://github.com/kfjahnke/envutil                               */
/*                                                                      */
/*    Please direct questions, bug reports, and contributions to        */
/*                                                                      */
/*    kfjahnke+envutil@gmail.com                                        */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

// A utility to extract an image from an environment. This program takes
// a 2:1 lat/lon environment or a 1:6 cubemap image as input and produces
// output in the specified orientation, projection, field of view and
// extent. For CL arguments, try 'envutil --help'. The program can also
// create environment images - just pass 'spherical' or 'cubemap' as
// output projection, 360 or 90 degrees hfov, respectively, and an
// appropriate output size. This ability can be used to convert from
// one environment to another, optionally with an arbitrary 3D rotation.
//
// The output projection can be one of "spherical", "cylindrical",
// "rectilinear", "stereographic", "fisheye" or "cubemap". The geometrical
// extent of the output is set up most conveniently by passing --hfov, the
// horizontal field of view of the output. The x0, x1, y0, and y1 parameters
// allow passing specific extent values (in model space units), which should
// rarely be necessary. To specify the orientation of the 'virtual camera',
// pass Euler angles yaw, pitch and roll - they default to zero: a view
// 'straight ahead' to the point corresponding to the center of the
// environment image with no camera roll. The size of the output is
// given by --width and --height. You must pass an output filename
// with --output; --input specifies the environment image.
//
// Per default, envutil uses 'twining' - inlined oversampling with
// subsequent weighted pixel binning. The default with this method is
// to use a simple box filter whose specific parameterization is set
// up automatically. Additional/ parameters can change the amount of
// oversampling and add gaussian weights to the filter parameters.
// To switch twining off and use direct interpolation from the source
// image data, pass --twine 0. 'twining' is quite fast (if the number
// of filter taps isn't very large. When down-scaling, the parameter
// 'twine' should be at least the same as the scaling factor to avoid
// aliasing. When upscaling, larger twining values will slighly soften
// the output and suppress the star-shaped artifacts typical for bilinear
// interpolation. Twining is new and this is a first approach. The method
// is intrinsically very flexible (it's based on a generalization of
// convolution), and the full flexibility isn't accessible in 'extract'
// with the parameterization as it stands now, but it's already quite
// useful with the few parameters I offer.
//
// The program uses zimt as it's 'strip-mining' and SIMD back-end, and
// sets up the pixel pipelines using zimt's functional composition tools.
// This allows for terse programming, and the use of a functional
// paradigm allows for many features to be freely combined - a property
// which is sometimes called 'orthogonality'. What you can't combine in
// 'extract' is twining and interpolation with OIIO - this is pointless,
// because OIIO offers all the anti-aliasing and quality interpolation
// one might want, and using twining on top would not improve the
// output. Currently, the build is set up to produce binary for AVX2-
// -capable CPUs - nowadays most 'better' CPUs support this SIMD ISA.
// When building for other (and non-i86) CPUs, suitable parameters should
// be passed to the compiler (you'll have to modify the CMakeLists.txt).
// I strongly suggest you install highway on your system - the build
// will detect and use it to good effect. This is a build-time dependency
// only. Next-best (when using i86 CPUs up to AVX2) is Vc, the fall-back
// is to use std::simd, and even that can be turned off if you want to
// rely on autovectorization; zimt structures the processing so that it's
// autovectorization-friendly and performance is still quite good that way.
// I have managed to build envutil on Linux, macOS (on intel CPUs) and
// Windows. The build adapts to the given system and expects a set of
// dependencies (OpenImageIO, Imath, ffmpeg), zimt code is provided in
// this repository. The macOS build fulfilled the dependencies with
// macPorts, the Windows build used msys2/mingw64.

// envutil_dispatch.cc provides the dispatch_base pointer to the
// ISA-specific rendering code via get_dispatch, which obtains the
// pointer via a call to HWY_DYNAMIC_DISPATCH. For single-ISA builds,
// there is no dispatching, and the payload code resides in namespace
// project::zsimd, following the zimt convention for the nested
// namespace's name if MULTI_SIMD_ISA is not defined. But we still
// have to define get_dispatch then - it simply delegates to
// _get_dispatch in the nested namespace.

#include "envutil_dispatch.h"
#include "pto.h"

#ifndef MULTI_SIMD_ISA

// for single-ISA builds, we directly code get_dispatch here:

namespace project
{
  namespace zsimd
  {
    const dispatch_base * const _get_dispatch() ;
  } ;

  const dispatch_base * const get_dispatch()
  {
    return zsimd::_get_dispatch() ;
  }
} ;

#endif

// a large part of the code in this file is dedicated to processing
// command line arguments. We use OpenImageIO's ArgParse object, which
// is similar to python's argparse. The result of gleaning the arguments
// is held in member variables of the 'args' object. If any command line
// arguments aren't acceptable, the program will terminate with an
// exception. Beyond simply parsing the arguments, the code in this
// object does some calculations to process the arguments further and
// provide more palatable values to the program. Note: angles are passed
// in degrees, but internally, only radians are used. The conversion is
// done right after the parameter acquisition.

#include <OpenImageIO/filesystem.h>
#include <OpenImageIO/argparse.h>

using OIIO::ArgParse ;
using OIIO::Filesystem::convert_native_arguments ;

bool facet_spec::init ( int argc , const char ** argv )
{
  convert_native_arguments(argc, (const char**)argv);
  ArgParse ap;

  ap.add_argument("--facet %s:IMAGE %s:PROJECTION %F:HFOV %F:YAW %F:PITCH %F:ROLL", &filename , &projection_str, &hfov, &yaw, &pitch, &roll)
    .help("load oriented source ('facet') image") ;

  if (ap.parse(argc, argv) < 0 ) {
      std::cerr << ap.geterror() << std::endl;
      ap.print_help();
      return false ;
  }

  if ( hfov <= 0.0 )
  {
    std::cerr << "facet hfov invalid: " << hfov << std::endl ;
    return false ;
  }

  // initially, the asset key is just the filename.

  asset_key = filename ;

  // determine facet image's projection
  int prj = 0 ;
  for ( const auto & p : projection_name )
  {
    if ( p == projection_str )
      break ;
    ++ prj ;
  }
  projection = projection_t ( prj ) ;
  if ( prj >= PRJ_NONE )
  {
    std::cerr << "facet projection invalid: " << projection_str << std::endl ;
    return false ;
  }

  // all seems well so far, let's open the image

  get_image_metrics() ;

  brighten = 1.0f ;

  return true ;
}

// to avoid having to pass the arguments around, we us a global
// 'args' object (declared extern in basic.h)

arguments args ;

void arguments::init ( int argc , const char ** argv )
{
  // we're using OIIO's argparse, since we're using OIIO anyway.
  // This is a convenient way to glean arguments on all supported
  // platforms - getopt isn't available everywhere.

  convert_native_arguments(argc, (const char**)argv);
  ArgParse ap;

  ap.intro("envutil: convert and create extracts from environment images\n")
    .usage("envutil [options] --input INPUT --output OUTPUT");

  ap.arg("-v", &verbose)
    .help("Verbose output");

  // the options are grouped thematically

  // mandatory options

  ap.separator("  mandatory options:");

  ap.arg("--output OUTPUT")
    .help("output file name (mandatory)")
    .metavar("OUTPUT");

  // important options which have defaults

  ap.separator("  important options which have defaults:");

  ap.arg("--projection PRJ")
    .help("projection used for the output image(s) (default: rectilinear)")
    .metavar("PRJ");

  ap.arg("--hfov ANGLE")
    .help("horiziontal field of view of the output (default: 90)")
    .metavar("ANGLE");

  ap.arg("--width EXTENT")
    .help("width of the output (default: 1024)")
    .metavar("EXTENT");

  ap.arg("--height EXTENT")
    .help("height of the output (default: same as width)")
    .metavar("EXTENT");

  ap.arg("--support_min EXTENT")
    .help("minimal additional support around the cube face proper")
    .metavar("EXTENT");

  ap.arg("--tile_size EXTENT")
    .help("tile size for the internal representation image")
    .metavar("EXTENT");

  // parameters for single-image output

  ap.separator("  additional parameters for single-image output:");

  // tentative:

  ap.arg("--single FACET")
    .help("render an image like facet FACET")
    .metavar("FACET");

  ap.arg("--split FORMAT_STRING")
    .help("create a 'single' facet for all facets in a PTO")
    .metavar("FORMAT_STRING");

  ap.arg("--yaw ANGLE")
    .help("yaw of the virtual camera")
    .metavar("ANGLE");

  ap.arg("--pitch ANGLE")
    .help("pitch of the virtual camera")
    .metavar("ANGLE");

  ap.arg("--roll ANGLE")
    .help("roll of the virtual camera")
    .metavar("ANGLE");

  ap.arg("--x0 EXTENT")
    .help("low end of the horizontal range")
    .metavar("EXTENT");

  ap.arg("--x1 EXTENT")
    .help("high end of the horizontal range")
    .metavar("EXTENT");

  ap.arg("--y0 EXTENT")
    .help("low end of the vertical range")
    .metavar("EXTENT");

  ap.arg("--y1 EXTENT")
    .help("high end of the vertical range")
    .metavar("EXTENT");

  // parameters for multi-image and video output

  ap.separator("  additional parameters for multi-image and video output:");

  ap.arg("--seqfile SEQFILE")
    .help("image sequence file name (optional)")
    .metavar("SEQFILE");

  ap.arg("--codec CODEC")
    .help("video codec for video sequence output (default: libx265)")
    .metavar("CODEC");

  ap.arg("--mbps MBPS")
    .help("output video with MBPS Mbit/sec (default: 8)")
    .metavar("MBPS");

  ap.arg("--fps FPS")
    .help("output video FPS frames/sec (default: 60)")
    .metavar("FPS");

  // interpolation options

  ap.separator("  interpolation options:");

  ap.arg("--prefilter DEG")
    .help("prefilter degree (>= 0) for b-spline-based interpolations")
    .metavar("DEG");

  ap.arg("--degree DEG")
    .help("degree of the spline (>= 0) for b-spline-based interpolations")
    .metavar("DEG");

  // parameters for twining

  ap.separator("  parameters for twining (--twine 0 switches twining off)");

  ap.arg("--twine TWINE")
    .help("use twine*twine oversampling - omit this arg for automatic twining")
    .metavar("TWINE");

  ap.arg("--twf_file TWF_FILE")
    .help("read twining filter kernel from TWF_FILE (switches twining on)")
    .metavar("TWF_FILE");

  ap.arg("--twine_normalize", &twine_normalize)
    .help("normalize twining filter weights gleaned from a file");

  ap.arg("--twine_precise", &twine_precise)
    .help("project twining basis vectors to tangent plane");

  ap.arg("--twine_width WIDTH")
    .help("widen the pick-up area of the twining filter")
    .metavar("WIDTH");

  ap.arg("--twine_density DENSITY")
    .help("increase tap count of an 'automatic' twining filter")
    .metavar("DENSITY");

  ap.arg("--twine_sigma SIGMA")
    .help("use a truncated gaussian for the twining filter (default: don't)")
    .metavar("SIGMA");

  ap.arg("--twine_threshold THR")
    .help("discard twining filter taps below this threshold")
    .metavar("THR");

  ap.separator("  parameters for mounted (facet) image input:");
  
  ap.arg("--pto PTOFILE")
    .help("panotools script in hugin PTO dialect (optional)")
    .metavar("PTOFILE");

  ap.add_argument("--facet %L:IMAGE %L:PROJECTION %L:HFOV %L:YAW %L:PITCH %L:ROLL",
                  &facet_name_v , &facet_projection_v, &facet_hfov_v, &facet_yaw_v, &facet_pitch_v, &facet_roll_v )
    .help("load oriented non-environment source image") ;

  ap.add_argument("--pto_line %L:LINE", &addenda )
    .help("add (trailing) line of PTO code") ;

  ap.add_argument("--solo FACET_INDEX")
    .help("show only this facet (indexes starting from zero)")
    .metavar("FACET_INDEX") ;

  ap.add_argument("--mask_for FACET_INDEX")
    .help("paint this facet white, all others black")
    .metavar("FACET_INDEX") ;

  ap.add_argument("--nchannels CHANNELS")
    .help("produce output with CHANNELS channels (1-4)")
    .metavar("CHANNELS") ;

  if (ap.parse(argc, argv) < 0 ) {
      std::cerr << ap.geterror() << std::endl;
      ap.print_help();
      assert ( false ) ;
  }

  if (!metamatch.empty()) {
      field_re.assign(metamatch, std::regex_constants::extended
                                  | std::regex_constants::icase);
  }

  // extract the arguments from the argparser, parse the projection

  output = ap["output"].as_string ( "" ) ;
  seqfile = ap["seqfile"].as_string ( "" ) ;
  pto_file = ap["pto"].as_string ( "" ) ;
  twf_file = ap["twf_file"].as_string ( "" ) ;
  split = ap["split"].as_string ( "" ) ;
  // fct_file = ap["fct_file"].as_string ( "" ) ;
  codec = ap["codec"].as_string ( "libx265" ) ; 
  mbps = ( 1000000.0 * ap["mbps"].get<float> ( 8.0 ) ) ;
  fps = ap["fps"].get<int>(60);
  prefilter_degree = ap["prefilter"].get<int>(-1);
  spline_degree = ap["degree"].get<int>(1);
  twine = ap["twine"].get<int>(-1);
  twine_width = ap["twine_width"].get<float>(1.0);
  twine_density = ap["twine_density"].get<float>(1.0);
  twine_sigma = ap["twine_sigma"].get<float>(0.0);
  twine_threshold = ap["twine_threshold"].get<float>(0.0);
  x0 = ap["x0"].get<float> ( 0.0 ) ;
  x1 = ap["x1"].get<float> ( 0.0 ) ;
  y0 = ap["y0"].get<float> ( 0.0 ) ;
  y1 = ap["y1"].get<float> ( 0.0 ) ;
  width = ap["width"].get<int> ( 0 ) ;
  height = ap["height"].get<int> ( 0 ) ;
  hfov = ap["hfov"].get<float>(90.0);
  tile_size = ap["tile_size"].get<int> ( 64 ) ;
  support_min = ap["support_min"].get<int> ( 8 ) ;
  if ( hfov != 0.0 )
    x0 = x1 = y0 = y1 = 0 ;
  yaw = ap["yaw"].get<float>(0.0);
  pitch = ap["pitch"].get<float>(0.0);
  roll = ap["roll"].get<float>(0.0);
  projection_str = ap["projection"].as_string ( "rectilinear" ) ;

  if ( prefilter_degree < 0 )
    prefilter_degree = spline_degree ;

  // determine output projection
  int prj = 0 ;
  for ( const auto & p : projection_name )
  {
    if ( p == projection_str )
      break ;
    ++ prj ;
  }
  projection = projection_t ( prj ) ;

  if ( pto_file == std::string() && addenda.size() == 0 )
    assert ( args.facet_name_v.size() > 0 ) ;
  assert ( output != std::string() || split != std::string() ) ;

  bool ignore_p_line = false ;

  if ( width == 0 )
  {
    width = 1024 ;
  }
  else
  {
    ignore_p_line = true ;
  }

  if ( projection == CUBEMAP || projection == BIATAN6 )
  {
    height = 6 * width ;
    assert ( hfov >= 90.0 ) ;
  }
  if ( projection == SPHERICAL && height == 0 )
  {
    if ( width & 1 )
      ++width ;
    height = width / 2 ;
  }
  if ( height == 0 )
    height = width ;

  bool p_line_present = false ;
  projection_t p_line_projection ;
  std::size_t p_line_width ;
  std::size_t p_line_height ;
  double p_line_hfov ;
  float eev_sum = 0.0f ;
  int eev_count = 0 ;

  if ( pto_file != std::string() || addenda.size() > 0 )
  {
    if ( verbose )
      std::cout << "processing PTO file " << pto_file << std::endl ;

    pto_parser_type parser ;
    auto success = parser.read_pto_file ( pto_file , addenda ) ;

    if ( verbose )
      std::cout << "PTO file parse "
                << ( success ? "succeeded" : "failed" )
                << std::endl ;

    if ( ! success )
      exit ( -1 ) ;

    // instead of using straight std::stod on the values in the
    // i-line, we use the lambda 'glean', which returns zero if
    // the value is not set at all. One example is Tpp/Tpy,
    // which isn't present in older panoramas.

    auto glean = [] ( const std::string & str ) -> double
    {
      if ( str == std::string() )
        return 0.0 ;
      return std::stod ( str ) ;
    } ;

    if ( ! ignore_p_line )
    {
      auto & p_line_list ( parser.line_group [ "p" ] ) ;
      if ( p_line_list.size() != 0 )
        p_line_present = true ;

      for ( auto & p_line : p_line_list )
      {
        auto & dir ( p_line.field_map ) ;
        int prj = std::stoi ( dir [ "f" ] ) ;
        if ( prj == 0 )
          p_line_projection = RECTILINEAR ;
        else if ( prj == 1 )
          p_line_projection = CYLINDRICAL ;
        else if ( prj == 2 )
          p_line_projection = SPHERICAL ;
        else if ( prj == 3 )
          p_line_projection = FISHEYE ;
        else if ( prj == 4 )
          p_line_projection = STEREOGRAPHIC ;
        else
        {
          // TODO: allow override, better error message, check
          // abort code (currently produces a zombie process)
          std::cerr << "can't handle PTO projection code "
                    << prj << " in p-line" << std::endl ;

          p_line_projection = PRJ_NONE ;
        }
        p_line_width = std::stoi ( dir [ "w" ] ) ;
        p_line_height = std::stoi ( dir [ "h" ] ) ;
        p_line_hfov = ( M_PI / 180.0 ) * std::stod ( dir [ "v" ] ) ;

        std::string crop_str = dir [ "S" ] ;
        if ( crop_str != std::string() )
        {
          have_crop = true ;
          std::regex crop_regex ( "([0-9]+),([0-9]+),([0-9]+),([0-9]+)" ) ;
          std::smatch parts ;
          std::regex_match ( crop_str , parts , crop_regex ) ;
          p_crop_x0 = std::stoi ( parts[1].str() ) ;
          p_crop_x1 = std::stoi ( parts[2].str() ) ;
          p_crop_y0 = std::stoi ( parts[3].str() ) ;
          p_crop_y1 = std::stoi ( parts[4].str() ) ;
        }
        break ; // additional p-lines are ignored for now.
      }
    }

    auto & i_line_list ( parser.line_group [ "i" ] ) ;
    for ( auto & i_line : i_line_list )
    {
      auto & dir ( i_line.field_map ) ;
      facet_spec f ;
      f.facet_no = nfacets++ ;
      f.filename = dir [ "n" ] ;
      f.have_crop = false ;
      f.have_pto_mask = false ;
      if ( f.filename [ 0 ] == '"' )
        f.filename = f.filename.substr ( 1 , f.filename.size() - 2 ) ;
      f.asset_key = f.filename ;
      int prj = std::stoi ( dir [ "f" ] ) ;
      if ( prj == 0 )
        f.projection = RECTILINEAR ;
      else if ( prj == 1 )
        f.projection = CYLINDRICAL ;
      else if ( prj == 2 || prj == 3 ) // TODO: elliptic crop f. 2
        f.projection = FISHEYE ;
      else if ( prj == 4 )
         f.projection = SPHERICAL ;
      else if ( prj == 10 )
        f.projection = STEREOGRAPHIC ;
      else
      {
        std::cerr << "can't handle PTO projection code "
                  << prj << " in i-line" << std::endl ;
        exit ( -1 ) ;
      }

      f.projection_str = projection_name [ f.projection ] ;
      f.hfov = ( M_PI / 180.0 ) * std::stod ( dir [ "v" ] ) ;
      f.get_image_metrics() ;
      f.yaw = ( M_PI / 180.0 ) * glean ( dir [ "y" ] ) ;
      f.pitch = ( M_PI / 180.0 ) * glean ( dir [ "p" ] ) ;
      f.roll = ( M_PI / 180.0 ) * glean ( dir [ "r" ] ) ;
      f.tr_x = glean ( dir [ "TrX" ] ) ;
      f.tr_y = glean ( dir [ "TrY" ] ) ;
      f.tr_z = - glean ( dir [ "TrZ" ] ) ;
      f.tp_y = ( M_PI / 180.0 ) * glean ( dir [ "Tpy" ] ) ;
      f.tp_p = ( M_PI / 180.0 ) * glean ( dir [ "Tpp" ] ) ;
      f.tp_r = 0.0 ;
      f.shear_g = glean ( dir [ "g" ] ) / f.height ;
      f.shear_t = glean ( dir [ "t" ] ) / f.width ;
      f.step = get_step ( f.projection , f.width ,
                          f.height , f.hfov ) ;
      (extent_type&)f = get_extent ( f.projection , f.width ,
                                     f.height , f.hfov ) ;
      f.a = glean ( dir [ "a" ] ) ;
      f.b = glean ( dir [ "b" ] ) ;
      f.c = glean ( dir [ "c" ] ) ;
      f.h = glean ( dir [ "d" ] ) ;
      f.v = glean ( dir [ "e" ] ) ;
      f.process_lc() ;
      f.brighten = glean ( dir [ "Eev" ] ) ;
      if ( f.brighten != 0.0f )
      {
        eev_sum += f.brighten ;
        eev_count++ ;
      }
      std::string crop_str = dir [ "S" ] ;
      if ( crop_str != std::string() )
      {
        f.have_crop = true ;
        std::regex crop_regex ( "([0-9]+),([0-9]+),([0-9]+),([0-9]+)" ) ;
        std::smatch parts ;
        std::regex_match ( crop_str , parts , crop_regex ) ;
        f.crop_x0 = std::stoi ( parts[1].str() ) ;
        f.crop_x1 = std::stoi ( parts[2].str() ) ;
        f.crop_y0 = std::stoi ( parts[3].str() ) ;
        f.crop_y1 = std::stoi ( parts[4].str() ) ;
      }
      facet_spec_v.push_back ( f ) ;
    }

    const std::regex mask_corner_regex
      ( "([+-]?[0-9.]+)\\s([+-]?[0-9.]+)" ) ;

    auto & k_line_list ( parser.line_group [ "k" ] ) ;

    if ( verbose && k_line_list.size() != 0 )
      std::cout << "processing " << k_line_list.size()
                << " masks" << std::endl ;

    int mask_no = 0 ;
    for ( auto & k_line : k_line_list )
    {
      auto & dir ( k_line.field_map ) ;
      pto_mask_type mask ;
      mask.image = std::stoi ( dir [ "i" ] ) ;
      facet_spec_v [ mask.image ] . have_pto_mask = true ;
      mask.variant = std::stoi ( dir [ "t" ] ) ;
      mask.vertex_list = dir [ "p" ] ;

      auto start = std::sregex_iterator ( mask.vertex_list.begin() ,
                                          mask.vertex_list.end() ,
                                          mask_corner_regex ) ;
      auto end = std::sregex_iterator() ;
      for ( auto i = start ; i != end ; ++i )
      {
        // iterate over the items and separate name and value

        auto item = i->str() ;
        std::smatch parts ;
        std::regex_match ( item , parts , mask_corner_regex ) ;

        auto xs = parts[1].str() ;
        auto ys = parts[2].str() ;

        mask.vx.push_back ( std::stod ( xs ) ) ;
        mask.vy.push_back ( std::stod ( ys ) ) ;

        // TODO: this works only for natively landscape images, for portrait
        // the coordinates have to be modified to refer to the rotated image
      }
      
      if ( mask.variant != 0 )
      {
        // TODO: we might also process exclude masks for 'all images with
        // the same lens', but so far, we have no concept of a 'lens' in
        // envutil. I'm also unsure about 'include masks'. The wiki is
        // a bit vague insofar as it states that PTO masks are more like
        // hints to the stitcher (especially include masks). Definitely,
        // an include mask is nothing which can be implemented as an
        // attribute of a facet - it affects all facets which are touched
        // by the mask.

        std::cerr << "warning: mask type not implemented: " << mask.variant
                  << " this mask will be ignored" << std::endl ;
      }

      // more for reference:

      pto_mask_v.push_back ( mask ) ;

      // the exclude masks are stored with the facet they refer to:

      auto & fct ( facet_spec_v [ mask.image ] ) ;
      fct.pto_mask_v.push_back ( mask ) ;

      // we suffix the sequential numbers of all used masks to the
      // asset key to avoid reusing facets with the same name but
      // different masking/cropping.

      std::string suffix ( "." ) ;
      if ( fct.filename == fct.asset_key )
      {
        suffix += args.pto_file + "." ;
        fct.have_pto_mask = true ;
      }
      suffix += std::to_string ( mask_no ) ;
      fct.asset_key += suffix ;
      ++ mask_no ;
    }
  }

  if ( verbose )
  {
    for ( const auto mask: pto_mask_v )
      std::cout << mask ;
    if ( have_crop )
      std::cout << "p-line crop: " << p_crop_x0 << " " << p_crop_x1
                << " " << p_crop_y0 << " " << p_crop_y1 << std::endl ;
  }

  // add facets given as single --facet arguments. Even if the argument
  // occurs before a --pto argument, the 'free' facets are processed
  // after the ones in the PTO file to simplify facet numbering - it's
  // easiest to keep the PTO-internal numbering as-is: there are back
  // references and other references to the facet number in PTO.

  int n_free_facets = facet_name_v.size() ;
  for ( int i = 0 ; i < n_free_facets ; i++ )
  {
    facet_spec fspec ;
  
    const char * spec[8] ;
    spec [ 0 ] = "facet_spec" ;
    spec [ 1 ] = "--facet" ;
    spec [ 2 ] = facet_name_v[i].c_str() ;
    spec [ 3 ] = facet_projection_v[i].c_str() ;
    spec [ 4 ] = facet_hfov_v[i].c_str() ;
    spec [ 5 ] = facet_yaw_v[i].c_str() ;
    spec [ 6 ] = facet_pitch_v[i].c_str() ;
    spec [ 7 ] = facet_roll_v[i].c_str() ;
    bool success = fspec.init ( 8 , spec ) ;
    if ( ! success )
    {
      std::cerr << "parse of 'facet' argument with index " << i
                << " failed" << std::endl ;
      exit ( -1 ) ;
    }
    fspec.facet_no = i + nfacets ;
    fspec.hfov *= M_PI / 180.0 ;
    fspec.yaw *= M_PI / 180.0 ;
    fspec.pitch *= M_PI / 180.0 ;
    fspec.roll *= M_PI / 180.0 ;
    fspec.step = get_step ( fspec.projection , fspec.width ,
                            fspec.height , fspec.hfov ) ;
    (extent_type&) fspec = get_extent ( fspec.projection , fspec.width ,
                                        fspec.height , fspec.hfov ) ;
    fspec.tr_x = fspec.tr_y = fspec.tr_z = 0.0 ;
    fspec.tp_y = fspec.tp_p = fspec.tp_r = 0.0 ;
    fspec.shear_g = fspec.shear_t = 0.0 ;
    fspec.a = fspec.b = fspec.c = fspec.h = fspec.v = 0.0 ;
    fspec.process_lc() ;
    fspec.shift_only = false ;
    fspec.have_crop = false ;
    fspec.have_pto_mask = false ;
    fspec.asset_key = fspec.filename ;
    fspec.brighten = 0.0 ;
    facet_spec_v.push_back ( fspec ) ;
  }
  nfacets += n_free_facets ;
  assert ( nfacets ) ;

  solo = ap["solo"].get<int> ( -1 ) ;
  single = ap["single"].get<int> ( -1 ) ;

  if ( solo != -1 )
    assert ( solo < nfacets ) ;

  if ( single != -1 )
    assert ( single < nfacets ) ;

  // if there is only one facet, we set 'solo' to zero, this will
  // also result in facet zero's 'active' field being set true

  if ( nfacets == 1 )
    solo = 0 ;

  mask_for = ap["mask_for"].get<int> ( -1 ) ;
  if ( mask_for != -1 )
    assert ( mask_for < nfacets ) ;

  nchannels = 1 ;
  bool alpha_seen = false ;

  if ( eev_count > 0 )
  {
    eev_sum /= eev_count ;
  }

  for ( auto & m : facet_spec_v )
  {
    // we calculate the facet's 'brighten' value from Eev values in
    // the PTO file (if present) - this is correct for linear light
    // only!

    if ( eev_count )
    {
      if ( m.brighten == 0.0 )
      {
        m.brighten = 1.0f ;
      }
      else
      {
        // m.brighten initially holds the Eev value. now we set it
        // to hold a multiplicative parameter to adapt brightness
        // when this facet is used for rendering. Note: higher Eev
        // means less light, one step Eev is doubling/halving the
        // brightness. hence:

        m.brighten = pow ( 2.0 , m.brighten - eev_sum ) ;
      }
    }
    else
    {
      // if there were no Eev values, brightnes won't be modified.

      m.brighten = 1.0f ;
    }

    if ( m.have_pto_mask || m.have_crop )
    {
      if ( m.nchannels == 1 || m.nchannels == 3 )
        m.nchannels ++ ;
    }
    if ( m.nchannels == 2 || m.nchannels == 4 )
      alpha_seen = true ;

    // if this facet has a higher channel count than what we've
    // seen so far, update nchannels

    if ( m.nchannels > nchannels )
      nchannels = m.nchannels ;

    // if there is a mask_for argument, set only this facet's
    // 'masked' field to 1 - all other facets' 'mask' field to
    // zero. If there is no 'mask_for' argument, the 'mask'
    // field is set to -1 in all facets.

    if ( mask_for == -1 )
    {
      m.masked = -1 ;
    }
    else
    {
      if ( m.facet_no == mask_for )
        m.masked = 1 ;
      else
        m.masked = 0 ;
    }

    if ( verbose )
    {
      std::cout << "facet " << m.facet_no
                << " '" << m.filename << "' "
                << projection_name[m.projection]
                << " " << m.width << "*" << m.height << "#" << m.nchannels
                << " hfov: " << m.hfov * 180.0 / M_PI
                << " step: " << m.step << std::endl
                << "orientation  y:" << m.yaw * 180.0 / M_PI 
                << " p:" << m.pitch * 180.0 / M_PI
                << " r:" << m.roll * 180.0 / M_PI << std::endl
                << "translation tr_x:" << m.tr_x 
                << " tr_y:" << m.tr_y
                << " tr_z:" << m.tr_z << std::endl
                << "reprojection plane tp_y:" << m.tp_y * 180.0 / M_PI 
                << " tp_p:" << m.tp_p * 180.0 / M_PI
                << " tp_r:" << m.tp_r * 180.0 / M_PI << std::endl
                << "shear g: " << m.shear_g
                << " shear t: " << m.shear_t << " (in texture units)"
                << std::endl
                << "brighten: " << m.brighten << std::endl ;
      if ( m.masked != -1 )
        std::cout << "  b/w mask only: " << ( m.masked == 0
                                              ? "black"
                                              : "white" ) << std::endl ;
      if ( m.have_crop )
        std::cout << "  cropping active: " << m.crop_x0 << " "
                  << m.crop_x1 << " " << m.crop_y0
                  << " " << m.crop_y1 << std::endl ;
    }
  }

  // special case: there was at least one 2-channel facet, which is
  // taken as greyscale with alpha, resulting in alpha_seen true. And
  // there were also three-channel facets, so that nchannels is now
  // three (the maximum seen). In this case we up nchannels to four
  // to enforce alpha processing.

  if ( alpha_seen && nchannels == 3 )
  {
    std::cout << "found at least one image with transparency"
              << std::endl ;
    nchannels = 4 ;
  }

  // if there is an 'nchannels' argument, it overrides the value
  // we have set up from looking at the facets (the maximum seen)
  // the facet's own nchannels value depends on the image it's
  // made from, and the pixels will be forced to the global
  // nchannels value during processing. So the default behaviour
  // is to render the target image with the highest channel count
  // seen in any of the facets, but with a global nchannels
  // override, the target image will be rendered with that number
  // of channels unconditionally.

  int nch = ap["nchannels"].get<int> ( 0 ) ;
  if ( nch > 0 )
  {
    std::cout << "global nchannels override in arguments" << std::endl ;
    nchannels = nch ;
  }

  if ( verbose )
    std::cout << "global nchannels set to: " << nchannels << std::endl ;

  // if ( seqfile == std::string() )
  {
    if ( single >= 0 )
    {
      // 'single' means: produce output with the same metrics as the
      // given facet. The idea is to even produce the same distortions;
      // for now, to test, just the base metrics:

      if ( verbose )
        std::cout << "using '--single' argument to set output metrics"
                  << std::endl ;

      const auto & fspec = facet_spec_v [ single ] ;

      // take over the facet's geometry to the target geometry in args

      (facet_base&) args = fspec ;
    }
    else if ( p_line_present )
    {
      hfov = p_line_hfov ;
      projection = p_line_projection ;
      projection_str = projection_name [ projection ] ;
      width = p_line_width ;
      height = p_line_height ;
    }
    else
    {
      // convert angles to radians
      hfov *= M_PI / 180.0 ;
      yaw *= M_PI / 180.0 ;
      pitch *= M_PI / 180.0 ;
      roll *= M_PI / 180.0 ;
    }

    // single-image output (no sequence file given)
    // if there is a sequence file, the variables which are set
    // here will be set for every frame specified in the sequence,
    // and single-image parameters from the CL have no effect.

    if ( verbose)
    {
      std::cout << "output:             " << output << std::endl ;
      std::cout << "output projection:  " << projection_str << std::endl ;
      std::cout << "output width:       " << width << std::endl ;
      std::cout << "output height:      " << height << std::endl ;

      if ( hfov > 0.0 )
        std::cout << "output hfov:      "
                  << hfov * 180.0 / M_PI << std::endl ;

      std::cout << "virtual camera yaw: "
                << yaw * 180.0 / M_PI
                << " pitch: "
                << pitch * 180.0 / M_PI
                << " roll: "
                << roll * 180.0 / M_PI << std::endl ;
    }

    // calculate extent - a non-zero hfov overrides x0, x1, y0, and y1

    step = 0.0 ;
    if ( hfov != 0.0 )
    {
      (extent_type&) args
        = get_extent ( projection , width , height , hfov  ) ;
    }
    assert ( x0 < x1 ) ;
    assert ( y0 < y1 ) ;

    step = ( x1 - x0 ) / width ;

    if ( verbose )
    {
      if ( hfov == 0.0 )
      {
        std::cout << "extent calculated from hfov:"
                  << std::endl ;
      }
      else
      {
        std::cout << "extent gleaned from command line arguments:"
                  << std::endl ;
      }
      std::cout << "x0: " << x0 << " x1: " << x1 << std::endl ;
      std::cout << "y0: " << y0 << " y1: " << y1 << std::endl ;
      std::cout << "step: " << step << std::endl ;
    }
  }
}

void make_spread ( std::vector < zimt::xel_t < float , 3 > > & trg ,
                   int w = 2 ,
                   int h = 0 ,
                   float d = 1.0f ,
                   float sigma = 0.0f ,
                   float threshold = 0.0f )
{
  if ( w <= 2 )
    w = 2 ;
  if ( h <= 0 )
    h = w ;
  float wgt = 1.0 / ( w * h ) ;
  double x0 = - ( w - 1.0 ) / ( 2.0 * w ) ;
  double dx = 1.0 / w ;
  double y0 = - ( h - 1.0 ) / ( 2.0 * h ) ;
  double dy = 1.0 / h ;
  trg.clear() ;
  sigma *= - x0 ;
  double sum = 0.0 ;

  for ( int y = 0 ; y < h ; y++ )
  {
    for ( int x = 0 ; x < w ; x++ )
    {
      float wf = 1.0 ;
      if ( sigma > 0.0 )
      {
        double wx = ( x0 + x * dx ) / sigma ;
        double wy = ( y0 + y * dy ) / sigma ;
        wf = exp ( - sqrt ( wx * wx + wy * wy ) ) ;
      }
      zimt::xel_t < float , 3 >
        v { float ( d * ( x0 + x * dx ) ) ,
            float ( d * ( y0 + y * dy ) ) ,
            wf * wgt } ;
      trg.push_back ( v ) ;
      sum += wf * wgt ;
    }
  }

  double th_sum = 0.0 ;
  bool renormalize = false ;

  if ( sigma != 0.0 )
  {
    for ( auto & v : trg )
    {
      v[2] /= sum ;
      if ( v[2] >= threshold )
      {        
        th_sum += v[2] ;
      }
      else
      {
        renormalize = true ;
        v[2] = 0.0f ;
      }
    }
    if ( renormalize )
    {
      for ( auto & v : trg )
      {
        v[2] /= th_sum ;
      }
    }
  }

  if ( args.verbose )
  {
    if ( sigma != 0.0 )
    {
      std::cout << "using this twining filter kernel:" << std::endl ;
      for ( int y = 0 ; y < h ; y++ )
      {
        for ( int x = 0 ; x < w ; x++ )
        {
          std::cout << '\t' << trg [ y * h + x ] [ 2 ] ;
        }
        std::cout << std::endl ;
      }
    }
    else
    {
      std::cout << "using box filter for twining" << std::endl ;
    }
  }

  if ( renormalize )
  {
    auto help = trg ;
    trg.clear() ;
    for ( auto v : help )
    {
      if ( v[2] > 0.0f )
        trg.push_back ( v ) ;
    }
    if ( args.verbose )
    {
      std::cout << "twining filter taps after after thresholding: "
                << trg.size() << std::endl ;
    }
  }
}

// read the twining filter from a twf-file.

void read_twf_file ( std::vector < zimt::xel_t < float , 3 > > & trg )
{
  zimt::xel_t < float , 3 > c ;

  std::ifstream ifs ( args.twf_file ) ;
  assert ( ifs.good() ) ;

  double sum = 0.0 ;
  trg.clear() ;

  while ( ifs.good() )
  {
    ifs >> c[0] >> c[1] >> c[2] ;
    if ( ifs.eof() )
      break ;
    trg.push_back ( c ) ;
    sum += c[2] ;
  }

  for ( auto & c : trg )
  {
    c[0] *= args.twine_width ;
    c[1] *= args.twine_width ;
    if ( args.twine_normalize )
      c[2] /= sum ;
  }

  if ( args.verbose )
  {
    std::cout << args.twf_file << " yields twining filter kernel:"
              << std::endl ;
    for ( const auto & c : trg )
    {
       std::cout << "x: " << c[0] << " y: " << c[1]
                 << " w: " << c[2] << std::endl ;
    }
    if ( args.twine_normalize )
    {
      std::cout << "twining filter weights sum: 1.0" << std::endl ;
    }
    else
    {
      std::cout << "twining filter weights sum: " << sum << std::endl ;
    }
  }
  ifs.close() ;
}

void arguments::twine_setup()
{
  // first, we initialize the twining parameters.

  if ( twine != -1 )
  {
    // the user has passed a --twine argument other than -1,
    // which is the default value if no --twine argument is
    // passed on the command line. If the value is
    // zero, we take this to mean 'no twining'. Other values
    // result in twining - even a value of --twine 1, which
    // still produces a result where the pickup location is
    // modified by the single twining coefficient. Negative
    // values are set to zero.

    if ( twine < 0 )
      twine = 0 ;

    if ( twine > 0 )
    {
      assert ( twine_width > 0.0f ) ;
    }
  }
  else
  {
    // the twine argument was not passed, or passed -1 explicitly
    // to trigger automatic twining, which is the default anyway.
    // so we set up twining parameters to fit the job at hand.

    // first, figure out the magnification in the image center.
    // To make sure that the twine value is sufficiently large
    // for all contributing facets, we look for the smallest
    // 'step' value in the facet population and calculate
    // 'twine' accordingly. A single hi-res facet can therefore
    // result in a high twining factor which is overkill for
    // most of the result, so this detection and method is
    // conservative but potentially computationally expensive.

    double smallest_step = std::numeric_limits < double > :: max() ;
    
    if ( nfacets == 1 || solo > 0 )
    {
      // in solo mode, we calculate the twining factor to be suitable
      // for the single facet we're processing. This is debatable:
      // the result will differ slightly from what this facet's image
      // will be in a synopsis. The difference should be sub-pixel
      // only, though, and we avoid overkill for facets which have
      // lower resolution.

      smallest_step = args.facet_spec_v[solo].step ;
    }
    else
    { 
      for ( const auto & fspec : args.facet_spec_v )
      {
        smallest_step = std::min ( fspec.step , smallest_step ) ;
      }
    }

    double mag = smallest_step / step ;

    if ( mag > 1.0 )
    {
      // if the transformation magnifies, we use a moderate twine
      // size and a twine_width equal to the magnification, to
      // avoid the star-shaped artifacts from the bilinear
      // interpolation used for the lookup of the contributing
      // rays. If mag is small, the star-shaped artifacts aren't
      // really an issue, and twine values beyond, say, five
      // do little to improve the filter response, so we cap
      // the twine value at five, but lower it when approaching
      // mag 1, down to two - where the next case down starts.

      if ( args.spline_degree > 1 )
      {
        // if the substrate is a spline with a degree greater than
        // one, we only use moderate twining if there are several
        // facets, to soften potential collisions. For single-facet
        // operation, a b-spline is already smooth, so we don't twine.
        // TODO: consider tapering - the abrupt change to twine 1 if
        // mag > 1.0 may be problematic.

        if ( args.nfacets > 1 )
          twine = 3 ;
        else if ( mag < 2.0 ) // << tapering
          twine = 2 ;
        else
          twine = 1 ;
      }
      else
      {
        // we're using binlinear interpolation. since this is a
        // magnifying view, there may be star-shaped artifacts in the
        // 'ground truth' rendition, so we use twining.

        twine = std::min ( 5 , int ( 1.0 + mag ) ) ;

        // the twine width is set equal to the magnification: We want
        // to form the weighted sum over a field of as many pixels in
        // the 'ground truth' signal as correspond to a single pixel
        // in the target. Note that we're looking at the 'ground truth'
        // signal here, not the 'raw' image data! In a magnifying view,
        // the ground truth signal is the magnified version of the
        // image data - we only use twining here to avoid shortcomings
        // of the bilinear interpolation.

        twine_width = mag ;
      }
    }
    else
    {
      // for down-scaling, we use a twine size which depends on the
      // downscaling factor (reciprocal of 'mag') and a twine_width
      // of 1.0: we only want anti-aliasing. picking a sufficiently
      // large twine value guarantees that we're not skipping any
      // pixels in the ground truth signal due to undersampling.

      twine = int ( 1.0 + 1.0 / mag ) ;
      twine_width = 1.0 ;
    }

    if ( verbose )
    {
      std::cout << "automatic twining for magnification " << mag
                << ":" << std::endl ;
      std::cout << "twine: " << twine
                << " twine_width: " << twine_width
                << std::endl ;
    }
  }

  // we now have values for twine and twine_width, either passed
  // explicitly by the user or set automatically in the code above.

  if ( twine_density != 1.0f )
  {
    // if the user has passed twine_density, we use it as a
    // multiplicative factor to change twine - typically
    // twine_density will be larger than one, so we'll get
    // more filter taps.
  
    twine = std::round ( twine * twine_density ) ;

    if ( verbose )
    {
      std::cout << "applied twine_density " << twine_density
                << ": twine is now " << twine << std::endl ;
    }
  }

  // now we can make the 'spread' by calculating the twining
  // coefficients or by reading thw twf file.

  if ( args.twf_file == std::string() )
  {
    make_spread ( args.twine_spread , args.twine , args.twine ,
                  args.twine_width , args.twine_sigma ,
                  args.twine_threshold ) ;
  }
  else
  {
    // why only read this file now? because automatic twining may
    // have recalculated twine_width, which modulates the values
    // taken from the twf file. One might consider buffering the
    // file's original content.

    read_twf_file ( args.twine_spread ) ;
  }


  if ( twine )
  {
    // if twining is active, we should have a 'spread' by now.

    assert ( twine_spread.size() ) ;
  }
}

// cumulated frame rendering time

long rt_cumulated = 0 ;

int main ( int argc , const char ** argv )
{
  // process command line arguments - the result is held in a bunch
  // of member variables in the global 'args' object

  args.init ( argc , argv ) ;

  // are we to process a sequence file? If so, open the file

  bool have_seq = ( args.seqfile != std::string() ) ;
  std::ifstream seqstream ;
  if ( have_seq )
  {
    seqstream.open ( args.seqfile ) ;
    assert ( seqstream.good() ) ;
  }

  // obtain a dispatch_base pointer to the ISA-specific code.
  // get_dispatch is in envutil_dispatch.cc

  auto dp = project::get_dispatch() ;
  if ( args.verbose )
  {
    std::cout << "using " << dp->hwy_target_name << " ISA" << std::endl ;
  }

  // This loop will iterate over the renditions of single images.
  // If we're not running a sequence, there will only be one
  // iteration.

  int twine_as_passed = args.twine ;

  std::size_t frames = 0 ;
  do
  {
    // if we're running a sequence, we'll overwrite a few variables
    // in 'args' with values derived from the current line of input
    // read from the sequence file.

    if ( have_seq )
    {
      double seq_hfov , seq_yaw , seq_pitch , seq_roll ;
      seqstream >> seq_hfov >> seq_yaw >> seq_pitch >> seq_roll ;
      if ( ! seqstream.good() )
        break ;
      std::cout << "from seqfile: hfov: " << seq_hfov
                << " yaw: " << seq_yaw << " pitch: " << seq_pitch
                << " roll: " << seq_roll << std::endl ;
      args.hfov = seq_hfov * M_PI / 180.0 ;
      args.yaw = seq_yaw * M_PI / 180.0 ;
      args.pitch = seq_pitch * M_PI / 180.0 ;
      args.roll = seq_roll * M_PI / 180.0 ;
    
      (extent_type &) args = get_extent ( args.projection , args.width ,
                                          args.height , args.hfov  ) ;
      assert ( args.x0 < args.x1 ) ;
      assert ( args.y0 < args.y1 ) ;
      args.step = ( args.x1 - args.x0 ) / args.width ;

      // (re-) set 'twine' to the value which was passed to envutil
      // initially. This is to trigger recalculation of the twining
      // parameters to adapt to possible changes in the course of
      // the sequence.

      args.twine = twine_as_passed ;
    }

    args.twine_setup() ;

    // find the parameters which are type-relevant to route to
    // the specialized code above. There are several stages of
    // 'roll_out' which move the parameterization given by run-time
    // arguments into type information.

    int nch = args.nchannels ;
    int ninp = ( args.twine == 0 ) ? 3 : 9 ;
    projection_t prj = args.projection ;

    // now call the ISA-specific rendering code via the dispatch_base
    // pointer received from get_dispatch.

    if ( args.split != std::string() )
    {
      image_series split_name ( args.split ) ;

      for ( int i = 0 ; i < args.nfacets ; i++ )
      {
        // if 'solo' is set (it would be -1 otherwise) we skip the
        // solo facet as a target - we already have it, so there's
        // little point in re-creating it.

        if ( i == args.solo )
          continue ;

        // for all other facets, we take the facet's geometry (encoded
        // in it's base class facet_base) over as target geometry and
        // then run a 'single' job with that target, which recreates
        // the facet from either the solo image (if that was given)
        // or from the synopsis formed from all source facets. The
        // latter variant (no solo facet) will show stitching artifacts
        // depending on the stitching process, and the 're-created'
        // facets will look as if they had been 're-created' from a
        // stitch done with envutil with the current settings. As of
        // this writing, this is quite rough-and-ready - it's
        // geometrically correct, but the facet's aren't blended,
        // so the seams may be visible.

        const auto & fspec = args.facet_spec_v [ i ] ;

        // take over the facet's geometry to the target geometry in args

        (facet_base&) args = fspec ;

        // request a 'single' job - TODO may become obsolete

        args.single = i ;

        // save the target with this filename:

        args.output = split_name[i] ;

        dp->payload ( nch , ninp , prj ) ;
      }
    }
    else
      dp->payload ( nch , ninp , prj ) ;
    ++frames ;
  }
  while ( have_seq ) ; // loop criterion for 'do' loop

  if ( args.verbose && have_seq )
  {
    double rt_avg = rt_cumulated ;
    rt_avg /= frames ;
    std::cout << "average frame rendering time: " << rt_avg
              << " msec" << std::endl ;
  }
}
