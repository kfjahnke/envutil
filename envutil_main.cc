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

// This program takes a 2:1 lat/lon environment, a 1:6 cubemap image or
// a set of several 'facet' images as input and produces output in the
// specified orientation, projection, field of view and extent. For CL
// arguments, try 'envutil --help'. Find documentation on envutil's
// project page on github:
// https://github.com/kfjahnke/envutil

// This file deals mainly with argument processing, it has the 'main'
// function and calls the 'payload' code which does the actual rendering.

// envutil_dispatch.cc provides the dispatch_base pointer to the
// ISA-specific rendering code via get_dispatch, which obtains the
// pointer via a call to HWY_DYNAMIC_DISPATCH. For single-ISA builds,
// there is no dispatching, and the payload code resides in namespace
// project::zsimd, following the zimt convention for the nested
// namespace's name if MULTI_SIMD_ISA is not defined. But we still
// have to define get_dispatch then - it simply delegates to
// _get_dispatch in the nested namespace.

#include <thread>
#include <chrono>

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
#include <OpenImageIO/color.h>
#include <OpenImageIO/imagebufalgo.h>

#include <iostream>
#include <string>

using namespace OIIO;


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

  bool read_hfov = false ;
  bool read_projection = false ;

  if ( hfov <= 0.0 )
  {
    if ( hfov == -1.0 )
    {
      // glean hfov from  image metadata
      read_hfov = true ;
    }
    else
    {
      std::cerr << "facet hfov invalid: " << hfov << std::endl ;
      return false ;
    }
  }

  // initially, the asset key is just the filename.

  asset_key = filename ;

  // determine facet image's projection

  int prj = 0 ;
  if ( projection_str == "metadata" )
  {
    read_projection = true ;
  }

  // all seems well so far, let's open the image

  get_image_metrics ( read_hfov , read_projection ) ;

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
    .usage("envutil [options...] --output OUTPUT");

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

  ap.arg("--synopsis MODE")
    .help("mode of composing several images (panorama or hdr_merge)")
    .metavar("MODE");

  ap.arg("--working_colour_space CSP")
    .help("colour space used for internal processing (default scene_linear)")
    .metavar("CSP");

  // parameters for single-image output

  ap.separator("  additional parameters for single-image output:");

  ap.arg("--output_colour_space CSP")
    .help("colour space used for output (default scene_linear)")
    .metavar("CSP");

  ap.arg("--single FACET")
    .help("render an image like facet FACET")
    .metavar("FACET");

  ap.arg("--split FORMAT_STRING")
    .help("create a 'single' image for all facets in a PTO")
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

  ap.arg("--brighten FACTOR")
    .help("multiplicative factor to darken/brighten output")
    .metavar("FACTOR");

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

  ap.arg("--twine_max MAX")
    .help("maximal twine factor when using automatic twining")
    .metavar("MAX");

  ap.separator("  parameters for mounted (facet) image input:");

  // tentative

  ap.add_argument("--photo %L:IMAGE", &photo_name_v )
    .help("load photographic image, interpreting metadata") ;

  ap.add_argument("--facet %L:IMAGE %L:PROJECTION %L:HFOV %L:YAW %L:PITCH %L:ROLL",
                  &facet_name_v , &facet_projection_v, &facet_hfov_v, &facet_yaw_v, &facet_pitch_v, &facet_roll_v )
    .help("load oriented non-environment source image") ;

  ap.arg( "--oiio %L:OPTION" , &oiio_option_v )
    .help("pass option to configure OIIO plugin (may be used repeatedly)") ;

  ap.arg("--input_colour_space CSP")
    .help("default colour space for input images (default: none)")
    .metavar("CSP");

  ap.arg("--pto PTOFILE")
    .help("panotools script in hugin PTO dialect (optional)")
    .metavar("PTOFILE");

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

  // colour spaces. The first one defines a blanket colour space for
  // all incoming facet images, which may be overridden with a 'Csp'
  // clause in an i-line. The second one is 'scene_linear' and used
  // for internal processing - this must use a linear representation
  // to be mathematically correct. The third is also scene_linear,
  // only JPG output is set to sRGB unconditionally. Note that both
  // 'sRGB' and 'scene_linear' are OIIO symbols, and if there is an
  // active OCIO config file (gleaned via the $OCIO environment
  // variable) these symbols are at times not processed as expected
  // or not understood. In such Situations you may need to use
  // colour space names from the OCIO config or delete the OCIO
  // environment variable. input_colour_space is, per default,
  // set to an empty string, expecting that OIIO will provide a
  // value fitting the data in the image file.

  std::string default_working_colour_space ;
  std::string default_output_colour_space ;

  // TODO: we need to find the currently used ColorConfig and query
  // it. this code does not work (gemini's idea)
  // if ( OpenColorIO::is_initialized() )
  // {
  //   if ( verbose )
  //     std::cout << "found OCIO configuration" << std::endl ;
  //   default_working_colour_space = "ACEScg" ;
  //   default_output_colour_space = "ACES2065-1" ;
  // }
  // else
  {
    // if ( verbose )
    //   std::cout << "OCIO configuration not found, using OIIO's built-in"
    //             << std::endl ;

    // I think that without an active OCIO config, scene_linear will
    // use Rec709. Will need some experimentation to get this right.
    // user can override the defaults and work with or without config.

    default_working_colour_space = "Linear" ;
    default_output_colour_space = "Linear" ;
  }

  // the input colour space should be gleaned from the image's metadata,
  // so we can't set a reasonable default here.
  input_colour_space
    = ap["input_colour_space"].as_string ( "" ) ;
  working_colour_space
    = ap["working_colour_space"].as_string ( default_working_colour_space ) ;
  colour_space
    = ap["output_colour_space"].as_string ( default_output_colour_space ) ;
  pto_file = ap["pto"].as_string ( "" ) ;
  twf_file = ap["twf_file"].as_string ( "" ) ;
  split = ap["split"].as_string ( "" ) ;
  synopsis = ap["synopsis"].as_string ( "panorama" ) ;
  prefilter_degree = ap["prefilter"].get<int>(-1);
  spline_degree = ap["degree"].get<int>(1);
  twine = ap["twine"].get<int>(-1);
  twine_width = ap["twine_width"].get<float>(1.0);
  twine_density = ap["twine_density"].get<float>(1.0);
  twine_sigma = ap["twine_sigma"].get<float>(0.0);
  twine_threshold = ap["twine_threshold"].get<float>(0.0);
  twine_max = ap["twine_max"].get<int>(8);
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
  brighten = ap["brighten"].get<float>(1.0);
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
    assert ( args.facet_name_v.size() > 0 || args.photo_name_v.size() > 0 ) ;
  assert ( output != std::string() || split != std::string() ) ;

  bool ignore_p_line = false ;

  // 'solo' may be set for 'unstitching' jobs, so we initialize it:

  solo = -1 ;

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
  double p_line_eev = 0.0 ;
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

    auto iglean = [] ( const std::string & str ) -> int
    {
      if ( str == std::string() )
        return 0 ;
      return std::stoi ( str ) ;
    } ;

    auto & c_line_list ( parser.line_group [ "c" ] ) ;

    for ( auto & c_line : c_line_list )
    {
      auto & dir ( c_line.field_map ) ;
      cp_t cp ;

      cp.t = iglean ( dir [ "t" ] ) ;
      cp.n = iglean ( dir [ "n" ] ) ;
      cp.N = iglean ( dir [ "N" ] ) ;

      cp.x = glean ( dir [ "x" ] ) ;
      cp.y = glean ( dir [ "y" ] ) ;
      cp.X = glean ( dir [ "X" ] ) ;
      cp.Y = glean ( dir [ "Y" ] ) ;

      cp_v.push_back ( cp ) ;
    }

    if ( verbose && cp_v.size() )
      std::cout << "PTO file contains " << cp_v.size()
                << " control points" << std::endl ;

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
        p_line_width = iglean ( dir [ "w" ] ) ;
        p_line_height = iglean ( dir [ "h" ] ) ;
        p_line_hfov = ( M_PI / 180.0 ) * glean ( dir [ "v" ] ) ;
        p_line_eev = glean ( dir [ "Eev" ] ) ;

        std::string crop_str = dir [ "S" ] ;
        if ( crop_str != std::string() )
        {
          store_cropped = true ;
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
      f.has_lens_crop = false ;
      f.has_pto_mask = false ;

      // envutil provides an extension to PT format: x 'Csp' clause.
      // It defines the colour space in which the facet is stored.
      // There seems to be soem uncertainty around which colour space
      // names OIIO/OCIO will recognize. As of this writing, my setup
      // requires to use names given in an OCIO config file while
      // 'generic' oiio names are not supported. Given here, in an
      // i-line, it decribes what's inside the image file. In a p-line,
      // it prescribes which colour space the date are to be converted
      // to when writing the data to disc. Internally, envutil uses
      // the scene_linear colour space. The Csp clause takes prefence
      // over the blanket input_colour_space argument, which is used
      // if no Csp clause is present.

      std::string csp = dir [ "Csp" ] ;
      if ( csp != std::string() )
      {
        if ( csp [ 0 ] == '"' )
          csp = csp.substr ( 1 , csp.size() - 2 ) ;
        if ( verbose )
          std::cout << "facet's native colour space (via Csp): "
                    << csp << std::endl ;
      }
      else
      {
        csp = args.input_colour_space ;
        if ( verbose )
          std::cout
            << "facet's colour space (via input_colour_space): "
            << csp << std::endl ;
      }
      f.colour_space = csp ;

      // envutil provides an extension to PT format: a 'Pano'
      // clause. This takes the data from a PTO file's p-line
      // over to the current facet. This is used for 'unstitching',
      // see the section on the '--split' argument in the
      // documentation.
      
      auto pano = dir [ "Pano" ] ;
      if ( pano != std::string() )
      {
        assert ( p_line_present == true ) ;
        f.filename = pano ;
        f.colour_space = colour_space ;
        f.asset_key = f.filename ;
        if ( f.filename [ 0 ] == '"' )
          f.filename = f.filename.substr ( 1 , f.filename.size() - 2 ) ;
        f.projection = p_line_projection ;
        f.hfov = p_line_hfov ;
        f.window_width = p_crop_x1 - p_crop_x0 ;
        f.window_height = p_crop_y1 - p_crop_y0 ;
        f.get_image_metrics() ;
        if ( store_cropped )
        {
          // make sure the image file has the crop window's size
          assert ( f.window_width == f.width ) ;
          assert ( f.window_height == f.height ) ;
          f.width = p_line_width ;
          f.height = p_line_height ;
          f.window_x_offset = p_crop_x0 ;
          f.window_y_offset = p_crop_y0 ;
          // make sure the metrics are realistic
          assert ( f.width >= f.window_x_offset + f.window_width ) ;
          assert ( f.height >= f.window_y_offset + f.window_height ) ;
        }
        else
        {
          f.window_width = f.width ;
          f.window_height = f.height ;
          f.window_x_offset = 0 ;
          f.window_y_offset = 0 ;
        }
        solo = f.facet_no ;
      }
      else
      {
        // 'regular' processing of the i-line

        f.filename = dir [ "n" ] ;
        if ( f.filename [ 0 ] == '"' )
          f.filename = f.filename.substr ( 1 , f.filename.size() - 2 ) ;
        f.asset_key = f.filename ;
        int prj = std::stoi ( dir [ "f" ] ) ;
        if ( prj == 0 )
          f.projection = RECTILINEAR ;
        else if ( prj == 1 )
          f.projection = CYLINDRICAL ;
        else if ( prj == 2 || prj == 3 )
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

        f.get_image_metrics() ;
        f.hfov = ( M_PI / 180.0 ) * std::stod ( dir [ "v" ] ) ;

        // The 'W' clause is another envutil extension to PTO format,
        // used to deal with cropped image input.

        std::string window_str = dir [ "W" ] ;

        if ( window_str != std::string() )
        {
          // as an extension to PTO format, we accept a parameter to
          // describe cropped input images. The four values passed after
          // the 'W' are the same which would be given with 'S' in a
          // p-line to describe output cropping. If the 'W' parameter
          // is given, the 'w' and 'h' paramters must also be present
          // and give the size of the 'total' or uncropped image, and
          // the hfov relates to this total size. The image size gleaned
          // by get_image_metrics is the size of the image on disk, so
          // it must be just the same as the size of the window given
          // with 'W'.

          std::regex crop_regex ( "([0-9]+),([0-9]+),([0-9]+),([0-9]+)" ) ;
          std::smatch parts ;
          std::regex_match ( window_str , parts , crop_regex ) ;
          int x0 = std::stoi ( parts[1].str() ) ;
          int x1 = std::stoi ( parts[2].str() ) ;
          int y0 = std::stoi ( parts[3].str() ) ;
          int y1 = std::stoi ( parts[4].str() ) ;
          f.window_x_offset = x0 ;
          f.window_y_offset = y0 ;
          f.window_width = x1 - x0 ;
          f.window_height = y1 - y0 ;
          assert ( f.window_width == f.width ) ;
          assert ( f.window_height == f.height ) ;
          f.width = iglean ( dir [ "w" ] ) ;
          f.height = iglean ( dir [ "h" ] ) ;
          assert ( f.width != 0 ) ;
          assert ( f.height != 0 ) ;
        }
      }

      f.projection_str = projection_name [ f.projection ] ;
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
      f.process_geometry() ;
      f.brighten = glean ( dir [ "Eev" ] ) ;
      if ( f.brighten != 0.0f )
      {
        eev_sum += f.brighten ;
        eev_count++ ;
      }
      std::string crop_str = dir [ "S" ] ;
      if ( crop_str != std::string() )
      {
        f.has_lens_crop = true ;
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
      facet_spec_v [ mask.image ] . has_pto_mask = true ;
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
        fct.has_pto_mask = true ;
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
    if ( store_cropped )
      std::cout << "p-line crop: " << p_crop_x0 << " " << p_crop_x1
                << " " << p_crop_y0 << " " << p_crop_y1 << std::endl ;
  }

  // TODO: process image metadata. For now all 'photos' are taken as
  // head-on 90-degree rectilinear.

  for ( const auto & filename : photo_name_v )
  {
    facet_name_v.push_back ( filename ) ;
    facet_projection_v.push_back ( "metadata" ) ; // use metadata
    facet_hfov_v.push_back ( "-1" ) ;             // use metadata
    facet_yaw_v.push_back ( "0" ) ;
    facet_pitch_v.push_back ( "0" ) ;
    facet_roll_v.push_back ( "0" ) ;
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
    fspec.process_geometry() ;
    fspec.has_lens_crop = false ;
    fspec.has_pto_mask = false ;
    fspec.asset_key = fspec.filename ;
    fspec.brighten = 0.0 ;
    facet_spec_v.push_back ( fspec ) ;
  }
  nfacets += n_free_facets ;
  assert ( nfacets ) ;

  // if 'solo' hasn't been set yet (e.g. by an 'unstitching' job)
  // we glean it now:

  if ( solo == -1 )
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
  if ( p_line_eev != 0.0 )
  {
    if ( verbose )
      std::cout << "p-line has Eev, hence Eev out = " << p_line_eev
                << std::endl ;
    eev_sum = p_line_eev ;
  }
  else
  {
    if ( verbose )
      std::cout << "no p-line Eev, hence Eev out = " << eev_sum
                << std::endl ;
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

    // finally, output brightness is multiplied on the facets' individual
    // brighten values - avoiding a final multiplication of the output
    // with that factor.

    if ( args.brighten != 1.0 )
    {
      m.brighten *= args.brighten ;
    }

    if ( m.has_pto_mask || m.has_lens_crop )
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
      if ( m.has_lens_crop )
        std::cout << " input cropping: " << m.crop_x0 << " "
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

      std::string save_csp = args.colour_space ;
      const auto & fspec = facet_spec_v [ single ] ;

      // take over the facet's geometry to the target geometry in args
      // but don't take over the colour space

      (facet_base&) args = fspec ;
      args.colour_space = save_csp ;
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
    assert ( x0 <= x1 ) ;
    assert ( y0 <= y1 ) ;

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
// TODO: improve file format: allow comments

void read_twf_file ( std::vector < zimt::xel_t < float , 3 > > & trg ,
                     double twine_width )
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
    std::cout << "processing twf file: " << args.twf_file << std::endl ;
    std::cout << "applying scaling factor twine_width: "
              << args.twine_width << std::endl ;
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
  if ( args.twf_file != std::string() )
  {
    // if there is a twf file, unconditionally switch twining on

    twine = 1 ;
  }

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

      // tentative: we clamp the twine value. without clamping,
      // renditions of e.g. rectiliear views with very large hfov
      // would take very long due to excessive twine values.
      // Clamping to a relatively generous value will still result
      // in a weighted sum of a good many 'ground truth' samples,
      // and since the sampling is evenly spread over the area
      // corresponding to a single output pixel, we get a population
      // of ground truth samples which is kind of 'representative',
      // and the desired antialiasing effect should still be achieved.
      // unwanted moiree-like artifacts may be further reduced by using
      // twining kernels with non-regular patterns, e.g. random point
      // clouds with gaussian weights - but these would have to be
      // introduced with w twf file, or one might consider dtithering.

      twine = std::min ( args.twine_max , twine ) ;
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

    read_twf_file ( args.twine_spread , args.twine_width ) ;
  }


  if ( twine )
  {
    // if twining is active, we should have a 'spread' by now.

    assert ( twine_spread.size() ) ;
  }

  if ( args.verbose )
  {
    std::cout << "final twining filter kernel:" << std::endl ;
    int ord = 0 ;
    for ( const auto & c : args.twine_spread )
    {
       std::cout << ord++ << "\tx:\t" << c[0] << "\ty:\t" << c[1]
                 << "\tw:\t" << c[2] << std::endl ;
    }
  }

}

// cumulated frame rendering time

long rt_cumulated = 0 ;

// 'core' is like a 'main' function, but only for a single job;
// typically the rendition of a single image or a split job. If
// the command line does not terminate with a single '-', 'core'
// is only called once with all arguments from the command line.
// If there is a trailing '-', envutil reads batches of arguments
// from cin, one line at a time. The batch of arguments is prepended
// with any arguments which may have occured before the '-' and the
// combined argument list is passed to 'core'. So the original CL
// arguments are repeated for each call of 'core'.
// This new modus operandi makes 'seqfiles' obsolete - the new
// method is much more flexible and powerful.

int core ( int argc , const char ** argv , bool tethered = false )
{
  // process command line arguments - the result is held in a bunch
  // of member variables in the global 'args' object

  args.facet_spec_v.clear() ;
  args.nfacets = 0 ;
  args.facet_name_v.clear() ;
  args.facet_projection_v.clear() ;
  args.facet_hfov_v.clear() ;
  args.facet_yaw_v.clear() ;
  args.facet_pitch_v.clear() ;
  args.facet_roll_v.clear() ;
  args.photo_name_v.clear() ;
  args.addenda.clear() ;
  args.init ( argc , argv ) ;
  args.tethered = tethered ;

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

      // request a 'single' job

      args.single = i ;
      args.store_cropped = false ;

      // save the target with this filename:

      args.output = split_name[i] ;

      dp->payload ( nch , ninp , prj ) ;
    }
  }
  else
  {
    if ( args.single != -1 )
      args.store_cropped = false ;
    dp->payload ( nch , ninp , prj ) ;
  }

  return 0 ;
}

// the simple tokenizer for arguments is in envutil_basic.cc

extern std::vector < std::string > tokenize ( const std::string input ) ;

// main handles command line and piped arguments and calls 'core',
// possibly several times, if a trailing '-' switches envutil into
// 'pipe' mode.

#include <chrono>

using namespace std::chrono_literals ;

#include "visor.h"

// handle_job receives a reference to a job_t object and processes
// the parameters in this object to create output. currently, this
// is merely a dummy. Later on, this will trigger a rendering job
// which attaches the produced data in the job_t object's 'data'
// member variable. for now, we simply move a pointer slowly through
// an array of pixel data. This avoids all allocation and data
// generation, and we get an idea about the raw speed of the process.

bool handle_job ( ipc_data_t & ipc , int job_pending )
{
  spec_t & spec ( ipc.spec_array [ job_pending ] ) ;

  if ( spec.serial_no == 0 )
    return false ;

  // for now we just produce data reflecting some of the spec's
  // data fields (yaw_cam, pitch_cam)

  // get the index of the frame buffer to use

  bool success = ipc.store.get ( spec.buffer_index ) ;
  assert ( success ) ;

  // and extract the corresponding buffer address

  auto * p = (std::uint8_t*)
    ipc.get_buffer_address ( spec.buffer_index ) . get() ;

  // guard against overly large windows

  if ( spec.snapshot == false )
  {
    // if it's not a snapshot job, we clamp the extent to the size
    // of the frame buffers in the shared memory. Currently, the clamp
    // is done silently.

    if (   spec.height_cam * spec.width_cam
         > ipc.desktop_width * ipc.desktop_height )
    {
      spec.height_cam = ipc.desktop_height ;
      spec.width_cam = ipc.desktop_width ;
    }
  }

  args.p_screen_data = (std::uint32_t*) p ;

  std::vector < const char * > visor_argv ;
  std::size_t visor_argc ;
  
  // for now, we set the argument vector up 'on foot':

  std::string sw = std::to_string ( spec.width_cam ) ;
  std::string sh = std::to_string ( spec.height_cam ) ;
  std::string sy = std::to_string ( spec.yaw_cam ) ;
  std::string sf = std::to_string ( spec.hfov_cam ) ;
  std::string sp = std::to_string ( spec.pitch_cam ) ;
  std::string sr = std::to_string ( spec.roll_cam ) ;
  std::string sb = std::to_string ( spec.brighten ) ;

  std::vector < const char *> nargv ;

  nargv.push_back ( "envutil" ) ;
  nargv.push_back ( "--output" ) ;
  if ( spec.snapshot )
  {
    nargv.push_back ( spec.filename.c_str() ) ;
  }
  else
  {
    nargv.push_back ( "none.jpg" ) ;
  }
  nargv.push_back ( "--twine" ) ;
  if ( spec.refine )
    nargv.push_back ( "-1" ) ;
  else
    nargv.push_back ( "0" ) ;
  nargv.push_back ( "--hfov" ) ;
  nargv.push_back ( "65" ) ;

  ipc.flat_args.extract ( visor_argc , visor_argv ) ;

  for ( std::size_t a = 1 ; a < visor_argc ; a++ )
    nargv.push_back ( visor_argv[a] ) ;
  
  nargv.push_back ( "--width" ) ;
  nargv.push_back ( sw.c_str() ) ;
  nargv.push_back ( "--height" ) ;
  nargv.push_back ( sh.c_str() ) ;
  nargv.push_back ( "--yaw" ) ;
  nargv.push_back ( sy.c_str() ) ;
  nargv.push_back ( "--pitch" ) ;
  nargv.push_back ( sp.c_str() ) ;
  nargv.push_back ( "--roll" ) ;
  nargv.push_back ( sr.c_str() ) ;
  nargv.push_back ( "--hfov" ) ;
  nargv.push_back ( sf.c_str() ) ;
  if ( spec.brighten != 1.0 )
  {
    nargv.push_back ( "--brighten" ) ;
    nargv.push_back ( sb.c_str() ) ;
  }

//   std::cout << "***** have " << nargv.size() << " args" << std::endl ;
//   
//   for ( std::size_t i = 1 ; i < nargv.size() ; i++ )
//   {
//     std::cout << "arg " << i << " "
//               << nargv [ i ] << std::endl ;
//   }

  int nargc = nargv.size() ;

  bool tethered = true ;
  if ( spec.snapshot )
  {
    tethered = false ;
    spec.snapshot = false ;
  }

  core ( nargc , nargv.data() , tethered ) ;

  return true ;
}

// code to test for availability of an OCIO config. this code and the code to
// invoke it (in the beginning of main() was generated by google's gemini.

bool is_ocio_config_active()
{
    // 1. Access the currently loaded global configuration (uses $OCIO or fallback).
    // Note: The accessor "default_colorconfig()" is the correct call for your version.
    const ColorConfig & config ( ColorConfig::default_colorconfig() ) ;
    
    std::string scene_linear_cs_name;

    // 2. Attempt to resolve the "scene_linear" role.
    // The role name will be empty or "linear" in the minimal fallback, 
    // but something like "ACEScg" in a full OCIO config.
    try {
        // This call is now correct for your OIIO version.
        scene_linear_cs_name = config.getColorSpaceNameByRole("scene_linear") ;
    } catch (...) {
        // If an exception occurs, treat it as inactive/fallback.
        return false;
    }

    // 3. Evaluation Logic

    // A. Check for number of colorspaces (The high-confidence indicator)
    // The fallback has 4-5 spaces (linear, sRGB, Rec709, raw). 
    // A full OCIO config has many more.
    if (config.getNumColorSpaces() > 10) { 
        return true; // Confident: Full OCIO config is active.
    }
    
    // B. Check the resolved name (The final check)
    // If the full OCIO config is active, the role usually maps to a name 
    // that is NOT the simple built-in 'linear' name.
    
    // We explicitly check if the resolved name is the simple fallback name.
    if (scene_linear_cs_name == "linear" || scene_linear_cs_name == "") {
        // If it maps to the built-in OIIO "linear" space, it's the fallback.
        return false;
    }

    // If the space count is low (<= 10) BUT the resolved role is neither "linear" nor empty,
    // we assume a custom, small OCIO config is loaded (i.e., OCIO is active).
    if (scene_linear_cs_name.length() > 0) {
        return true;
    }
    
    // If all else fails, it's the fallback.
    return false;
}

int main ( int argc , const char ** argv )
{
  if (is_ocio_config_active())
  {
    std::cout << "✅ OCIO is ACTIVE. Use 'scene_linear' (maps to " 
              << ColorConfig::default_colorconfig().getColorSpaceNameByRole("scene_linear") 
              << ") for the working space.\n";
  }
  else
  {
    std::cout << "⚠️ OCIO is INACTIVE (OIIO fallback). Use 'linear' for the working space.\n";
  }

  // test the last argument. if it's '+', we are running tethered
  // to a visor job, and if it's '-', we run in 'pipe mode', fetching
  // argument lines from cin.

  if ( std::string ( argv [ argc - 1 ] ) == "+" )
  {
    // tethered mode is handled by a commodity function in class
    // visor_protocol (in visor.h).

    visor_protocol::render_loop ( handle_job ) ;
    return 0 ;
  }

  if ( std::string ( argv [ argc - 1 ] ) != "-" )
    core ( argc , argv ) ;
  else
  {
    argc-- ;
    std::vector < std::string > sv ;
    std::vector < const char * > av ;
    while ( true )
    {
      std::string str ;
      if ( std::cin.eof() )
      {
        std::cout << "pipe has reached EOF" << std::endl ;
        break ;
      }
      else
      {
        if ( std::getline ( std::cin , str ) )
        {
          av.clear() ;
          args = arguments() ;
          sv = tokenize ( str ) ;
          for ( int i = 0 ; i < argc ; i++ )
            av.push_back ( argv[i] ) ;
          for ( const auto & t : sv )
          {
            std::cout << " <" << t << ">" ;
            av.push_back ( t.c_str() ) ;
          }
          std::cout << std::endl ;
          core ( av.size() , av.data() ) ;
        }
      }
    }
  }
  return 0 ;
}
