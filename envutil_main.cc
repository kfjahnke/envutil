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
// You can choose several different interpolation methods with the --itp
// cammand line argument. The default is --itp 1, which uses bilinear
// interpolation. This is fast and often good enough, especially if there
// are no great scale changes involved - so, if the output's resolution is
// similar to the input's. --itp -1 employs OpenImageIO (OIIO for short)
// for interpolation. Without further parameters, OIIO's default mode is
// used, which uses sophisticated, but slow methods to produce the output.
// All of OIIO's interpolation, mip-mapping and wrapping modes can be
// selected by using the relevant additional parameters. Finally, --itp -2
// uses 'twining' - inlined oversampling with subsequent weighted pixel
// binning. The default with this method is to use a simple box filter
// whose specific parameterization is set up automatically. Additional
// parameters can change the amount of oversampling and add gaussian
// weights to the filter parameters. 'twining' is quite fast (if the number
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

  ap.add_argument("--facet %s:IMAGE %s:PROJECTION %F:HFOV %F:YAW %F:PITCH %F:ROLL",
                  &filename , &projection_str, &hfov, &yaw, &pitch, &roll)
    .help("load oriented non-environment source image") ;

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

  auto inp = ImageInput::open ( filename ) ;

  if ( ! inp )
  {
    std::cerr << "failed to open facet image '"
              << filename << "'" << std::endl ;
    return false ;
  }

  const ImageSpec &spec = inp->spec() ;

  width = spec.width ;
  height = spec.height ;
  nchannels = spec.nchannels ;
  inp->close() ;

  return true ;
}

// to avoid having to pass the arguments around, we us a global
// 'args' object (declared extern in basic.h)

arguments args ;

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

  ap.arg("--input INPUT")
    .help("input file name (mandatory if no 'facet' is given)")
    .metavar("INPUT");

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

  // additional input parameters for cubemap input
  ap.separator("  additional input parameters for cubemap input:");

  ap.arg("--cbmfov ANGLE")
    .help("horiziontal field of view of cubemap input (default: 90)")
    .metavar("ANGLE");

  ap.arg("--cbm_prj PRJ")
    .help("projection for cubemaps (default: cubemap - alt. biatan6)")
    .metavar("PRJ");

  ap.arg("--support_min EXTENT")
    .help("minimal additional support around the cube face proper")
    .metavar("EXTENT");

  ap.arg("--tile_size EXTENT")
    .help("tile size for the internal representation image")
    .metavar("EXTENT");

  ap.arg("--ctc CTC")
    .help("pass '1' to interpret cbmfov as center-to-center (default 0)")
    .metavar("CTC");

  // parameters for single-image output
  ap.separator("  additional parameters for single-image output:");

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

  ap.arg("--itp ITP")
    .help("interpolator: 1 for b-spline, -1 for OIIO, -2 b-spline+twining")
    .metavar("ITP");

  ap.arg("--prefilter DEG")
    .help("prefilter degree (>= 0) for b-spline-based interpolations")
    .metavar("DEG");

  ap.arg("--degree DEG")
    .help("degree of the spline (>= 0) for b-spline-based interpolations")
    .metavar("DEG");

  // parameters for twining (with --itp -2)
  ap.separator("  parameters for twining (with --itp -2):");

  ap.arg("--twine TWINE")
    .help("use twine*twine oversampling - default: automatic settings")
    .metavar("TWINE");

  ap.arg("--twf_file TWF_FILE")
    .help("read twining filter kernel from TWF_FILE")
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

  // parameters for lookup with OpenImageIO (with --itp -1)
  ap.separator("  parameters for lookup with OpenImageIO (with --itp -1):");

  ap.arg("--tsoptions KVLIST")
    .help("OIIO TextureSystem Options: coma-separated key=value pairs")
    .metavar("KVLIST");

  ap.arg("--swrap WRAP")
    .help("OIIO Texture System swrap mode")
    .metavar("WRAP");

  ap.arg("--twrap WRAP")
    .help("OIIO Texture System twrap mode")
    .metavar("WRAP");

  ap.arg("--mip MIP")
    .help("OIIO Texture System mip mode")
    .metavar("MIP");

  ap.arg("--interp INTERP")
    .help("OIIO Texture System interp mode")
    .metavar("INTERP");

  ap.arg("--stwidth EXTENT")
    .help("swidth and twidth OIIO Texture Options")
    .metavar("EXTENT");

  ap.arg("--stblur EXTENT")
    .help("sblur and tblur OIIO Texture Options")
    .metavar("EXTENT");

  ap.arg("--conservative YESNO")
    .help("OIIO conservative_filter Texture Option - pass 0 or 1")
    .metavar("YESNO");

  ap.separator("  parameters for mounted (facet) image input:");
  
  ap.add_argument("--facet %L:IMAGE %L:PROJECTION %L:HFOV %L:YAW %L:PITCH %L:ROLL",
                  &facet_name_v , &facet_projection_v, &facet_hfov_v, &facet_yaw_v, &facet_pitch_v, &facet_roll_v)
    .help("load oriented non-environment source image") ;

  // TODO:
  // ap.arg("--fct_file FCT_FILE")
  //   .help("read multi-facet environment from FCT_FILE")
  //   .metavar("FCT_FILE");

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

  input = ap["input"].as_string ( "" ) ;
  output = ap["output"].as_string ( "" ) ;
  seqfile = ap["seqfile"].as_string ( "" ) ;
  twf_file = ap["twf_file"].as_string ( "" ) ;
  // fct_file = ap["fct_file"].as_string ( "" ) ;
  codec = ap["codec"].as_string ( "libx265" ) ; 
  mbps = ( 1000000.0 * ap["mbps"].get<float> ( 8.0 ) ) ;
  fps = ap["fps"].get<int>(60);
  itp = ap["itp"].get<int>(1);
  prefilter_degree = ap["prefilter"].get<int>(-1);
  spline_degree = ap["degree"].get<int>(1);
  twine = ap["twine"].get<int>(0);
  twine_width = ap["twine_width"].get<float>(1.0);
  twine_density = ap["twine_density"].get<float>(1.0);
  twine_sigma = ap["twine_sigma"].get<float>(0.0);
  twine_threshold = ap["twine_threshold"].get<float>(0.0);
  swrap = ap["swrap"].as_string ( "WrapDefault" ) ;
  twrap = ap["twrap"].as_string ( "WrapDefault" ) ;
  mip = ap["mip"].as_string ( "MipModeDefault" ) ;
  interp = ap["interp"].as_string ( "InterpSmartBicubic" ) ;
  tsoptions = ap["tsoptions"].as_string ( "automip=1" ) ;
  conservative_filter = ap["conservative"].get<int>(1) ;
  x0 = ap["x0"].get<float> ( 0.0 ) ;
  x1 = ap["x1"].get<float> ( 0.0 ) ;
  y0 = ap["y0"].get<float> ( 0.0 ) ;
  y1 = ap["y1"].get<float> ( 0.0 ) ;
  width = ap["width"].get<int> ( 0 ) ;
  stwidth = ap["stwidth"].get<float> ( 1 ) ;
  stblur = ap["stblur"].get<float> ( 0 ) ;
  height = ap["height"].get<int> ( 0 ) ;
  hfov = ap["hfov"].get<float>(90.0);
  cbmfov = ap["cbmfov"].get<float>(90.0);
  ctc = ap["ctc"].get<int>(0);
  tile_size = ap["tile_size"].get<int> ( 64 ) ;
  support_min = ap["support_min"].get<int> ( 8 ) ;
  if ( hfov != 0.0 )
    x0 = x1 = y0 = y1 = 0 ;
  yaw = ap["yaw"].get<float>(0.0);
  pitch = ap["pitch"].get<float>(0.0);
  roll = ap["roll"].get<float>(0.0);
  prj_str = ap["projection"].as_string ( "rectilinear" ) ;
  cbm_prj_str = ap["cbm_prj"].as_string ( "cubemap" ) ;

  if ( prefilter_degree < 0 )
    prefilter_degree = spline_degree ;

  // determine output projection
  int prj = 0 ;
  for ( const auto & p : projection_name )
  {
    if ( p == prj_str )
      break ;
    ++ prj ;
  }
  projection = projection_t ( prj ) ;

  // optionally, determine cubeface projection: cubemap = rectilinear,
  // biatan6 = biatan6

  prj = 0 ;
  for ( const auto & p : projection_name )
  {
    if ( p == cbm_prj_str )
      break ;
    ++ prj ;
  }
  cbm_prj = projection_t ( prj ) ;
  assert ( cbm_prj == CUBEMAP || cbm_prj == BIATAN6 ) ;

  if ( verbose )
    std::cout << "cbm_prj = " << projection_name [ cbm_prj ]
              << std::endl ;

  // TODO:
  // if ( fct_file != std::string() )
  // {
  //   // we set 'input' to the name of the mounted image; we use
  //   // 'input' as id for the environment 'asset', so we can't
  //   // just leave it blank.
  // 
  //   input = fct_file ;
  // }
  // else
  if ( args.facet_name_v.size() == 0 )
  {
    assert ( input != std::string() ) ;
  }
  assert ( output != std::string() ) ;
  if ( width == 0 )
    width = 1024 ;
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

  cbmfov *= M_PI / 180.0 ;

  // if ( mount_image != std::string() )
  // {
  //   mount_inp = ImageInput::open ( mount_image ) ;
  //   assert ( mount_inp ) ;
  // 
  //   const ImageSpec &spec = mount_inp->spec() ;
  // 
  //   mount_width = spec.width ;
  //   mount_height = spec.height ;
  //   nchannels = spec.nchannels ;
  //   mount_hfov *= M_PI / 180.0 ;
  //   switch ( mount_prj )
  //   {
  //     case SPHERICAL:
  //     case FISHEYE:
  //     case CYLINDRICAL:
  //       env_step = mount_hfov / mount_width ;
  //       break ;
  //     case RECTILINEAR:
  //       env_step = 2.0 * tan ( mount_hfov / 2.0 ) / mount_width ;
  //       break ;
  //     case STEREOGRAPHIC:
  //       env_step = 4.0 * tan ( mount_hfov / 4.0 ) / mount_width ;
  //       break ;
  //     default:
  //       break ;
  //   }
  //   mount_yaw *= M_PI / 180.0 ;
  //   mount_pitch *= M_PI / 180.0 ;
  //   mount_roll *= M_PI / 180.0 ;
  // }
  // else if ( fct_file != std::string() )
  // {
  //   std::cerr << "not implemented: input is a multi-facet file "
  //               << fct_file << std::endl ;
  // }
  // else
  if ( args.facet_name_v.size() == 0 )
  {
    // some member variables in the args object are gleaned from
    // the input image.

    // first we check for percent signs in the input filename. If
    // we find one, we assume that the input is a cubemap consisting
    // of six separate images following a naming scheme described by
    // the string in 'input' which is treated as a format string.

    multiple_input = false ;
    auto has_percent = input.find_first_of ( "%" ) ;
    if ( has_percent != std::string::npos )
    {
      // input must be a set of six cubeface images, that's the
      // only way how we accept a format string.

      cfs = cubeface_series ( input ) ;
      multiple_input = cfs.valid() ;
      if ( multiple_input )
      {
        // let's open the first cube face to extract the metrics.
        // the cubemap's 'load' routine will check all images in
        // turn, so we needn't do that here.

        inp = ImageInput::open ( cfs[0] ) ;
        assert ( inp ) ;

        const ImageSpec &spec = inp->spec() ;

        assert ( spec.width == spec.height ) ;
        env_width = spec.width ;
        env_height = spec.height * 6 ;
        nchannels = spec.nchannels ;
        env_projection = cbm_prj ;
        if ( ctc )
        {
          double half_md = tan ( cbmfov / 2.0 ) ;
          half_md *= ( ( env_width + 1.0 ) / env_width ) ;
          cbmfov = atan ( half_md ) * 2.0 ;
          if ( verbose )
            std::cout << "ctc is set, adjusted cbmfov to "
                      << ( cbmfov * 180.0 / M_PI ) << std::endl ;
        }
        env_step = cbmfov / env_width ;
      }
    }

    if ( ! multiple_input )
    {
      // we have a single image as input.

      inp = ImageInput::open ( input ) ;
      assert ( inp ) ;

      const ImageSpec &spec = inp->spec() ;

      env_width = spec.width ;
      env_height = spec.height ;
      nchannels = spec.nchannels ;

      assert (    env_width == env_height * 2
              || env_height == env_width * 6 ) ;

      if ( env_width == env_height * 2 )
      {
        env_projection = SPHERICAL ;
        env_step = 2.0 * M_PI / env_width ;
      }
      else if ( env_width * 6 == env_height )
      {
        if ( ctc )
        {
          double half_md = tan ( cbmfov / 2.0 ) ;
          half_md *= ( ( env_width + 1.0 ) / env_width ) ;
          cbmfov = atan ( half_md ) * 2.0 ;
          if ( verbose )
            std::cout << "ctc is set, adjusted cbmfov to "
                      << ( cbmfov * 180.0 / M_PI ) << std::endl ;
        }
        env_projection = cbm_prj ;
        // TODO: might be slightly different for biatan6, check!
        env_step = cbmfov / env_width ;
      }
      else
      {
        std::cerr << "input image must have 2:1 or 1:6 aspect ratio"
                  << std::endl ;
        exit ( -1 ) ;
      }
    }

    if ( verbose )
    {
      std::cout << "input: " << input << std::endl ;
      std::cout << "input width: " << env_width << std::endl ;
      std::cout << "input height: " << env_height << std::endl ;
      std::cout << "input has " << nchannels << " channels" << std::endl ;
      std::cout << "env_step: " << env_step << std::endl ;
      if ( itp == 1 )
        std::cout << "using direct bilinear interpolation"
        << std::endl ;
      else if ( itp == -1 )
        std::cout << "using OIIO for interpolation"
        << std::endl ;
      else 
        std::cout << "using twining over bilinear interpolation"
        << std::endl ;
      std::cout << "output width: " << width
                << " height: " << height << std::endl ;
    }
  }

  if ( seqfile == std::string() )
  {
    // single-image output (no sequence file given)
    // if there is a sequence file, the variables which are set
    // here will be set for every frame specified in the sequence,
    // and single-image parameters from the CL have no effect.

    if ( verbose)
    {
      std::cout << "output: " << output << std::endl ;
      std::cout << "output projection: " << prj_str << std::endl ;

      if ( hfov > 0.0 )
        std::cout << "output hfov: " << hfov << std::endl ;

      std::cout << "virtual camera yaw: " << yaw
                << " pitch: " << pitch
                << " roll: " << roll << std::endl ;
    }

    // convert angles to radians

    hfov *= M_PI / 180.0 ;
    yaw *= M_PI / 180.0 ;
    pitch *= M_PI / 180.0 ;
    roll *= M_PI / 180.0 ;

    if ( ( projection == CUBEMAP || projection == BIATAN6 ) && ctc )
    {
      double half_md = tan ( hfov / 2.0 ) ;
      half_md *= ( ( width + 1.0 ) / width ) ;
      hfov = atan ( half_md ) * 2.0 ;
      if ( verbose )
        std::cout << "cubemap output: ctc is set, adjusted hfov to "
                  << ( hfov * 180.0 / M_PI ) << std::endl ;
    }

    // calculate extent - a non-zero hfov overrides x0, x1, y0, and y1

    step = 0.0 ;
    if ( hfov != 0.0 )
    {
      auto extent = get_extent ( projection , width , height , hfov  ) ;
      x0 = extent.x0 ;
      x1 = extent.x1 ;
      y0 = extent.y0 ;
      y1 = extent.y1 ;
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

  if ( twf_file != std::string() )
  {
    // if user passes a twf-file, it's used to set up the twining
    // filter, and all the other twining-related arguments apart
    // from twine_width and twine_normalize are ignored.

    read_twf_file ( twine_spread ) ;
    assert ( twine_spread.size() ) ;
  }
  else if ( itp == -2 )
  {
    if ( twine == 0 )
    {
      // user has passed twine 0, or not passed anything - but itp
      // is -2. We set up the twining with automatically generated
      // parameters.

      // figure out the magnification in the image center as a
      // guideline

      double mag = env_step / step ;

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

        twine = std::min ( 5 , int ( 1.0 + mag ) ) ;
        twine_width = mag ;
      }
      else
      {
        // otherwise, we use a twine size which depends on the
        // downscaling factor (reciprocal of 'mag') and a twine_width
        // of 1.0: we only want anti-aliasing. picking a sufficiently
        // large twine value guarantees that we're not skipping any
        // pixels in the source (due to undersampling).

        twine = int ( 1.0 + 1.0 / mag ) ;
        twine_width = 1.0 ;
      }

      if ( twine_density != 1.0f )
      {
        // if the user has passed twine_density, we use it as a
        // multiplicative factor to change twine - typically
        // twine_density will be larger than one, so we'll get
        // more filter taps.

        double twine = twine * twine_density ;
        twine = std::round ( twine ) ;
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
    else
    {
      if ( verbose )
      {
        std::cout << "using fixed twine: " << twine
                  << " twine_width: " << twine_width
                  << std::endl ;
      }
    }
  }
}

// cumulated frame rendering time

long rt_cumulated = 0 ;

int main ( int argc , const char ** argv )
{
  // process command line arguments - the result is held in a bunch
  // of member variables in the global 'args' object

  args.init ( argc , argv ) ;

  // do we have facets?

  // TODO: this is a bit rough-and-ready - correct arguments will be
  // parsed correctly, but if the numeric parameters aren't receiving
  // a string representing a number, there is no error - the field
  // will contain zero.

  args.nfacets = args.facet_name_v.size() ;
  if ( args.nfacets )
  {
    facet_spec fspec ;
    for ( int i = 0 ; i < args.nfacets ; i++ )
    {
      const char * spec[8] ;
      spec [ 0 ] = "facet_spec" ;
      spec [ 1 ] = "--facet" ;
      spec [ 2 ] = args.facet_name_v[i].c_str() ;
      spec [ 3 ] = args.facet_projection_v[i].c_str() ;
      spec [ 4 ] = args.facet_hfov_v[i].c_str() ;
      spec [ 5 ] = args.facet_yaw_v[i].c_str() ;
      spec [ 6 ] = args.facet_pitch_v[i].c_str() ;
      spec [ 7 ] = args.facet_roll_v[i].c_str() ;
      bool success = fspec.init ( 8 , spec ) ;
      if ( ! success )
      {
        std::cerr << "parse of facet argument with index " << i
                  << " failed" << std::endl ;
        exit ( -1 ) ;
      }
      fspec.facet_no = i ;
      fspec.hfov *= M_PI / 180.0 ;
      fspec.yaw *= M_PI / 180.0 ;
      fspec.pitch *= M_PI / 180.0 ;
      fspec.roll *= M_PI / 180.0 ;
      args.facet_spec_v.push_back ( fspec ) ;
    }
    std::size_t nch = args.facet_spec_v[0].nchannels ;
    for ( auto & m : args.facet_spec_v )
    {
      assert ( m.nchannels == nch ) ;
      if ( args.verbose )
        std::cout << "facet " << m.facet_no
                  << " '" << m.filename << "' "
                  << projection_name[m.projection]
                  << " " << m.width << "*" << m.height << "#" << m.nchannels
                  << " hfov: " << m.hfov * 180.0 / M_PI
                  << " y:" << m.yaw * 180.0 / M_PI 
                  << " p:" << m.pitch * 180.0 / M_PI
                  << " r:" << m.roll * 180.0 / M_PI << std::endl ;
    }
    args.nchannels = nch ;
  }

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

      auto extent = get_extent ( args.projection , args.width ,
                                 args.height , args.hfov  ) ;
      args.x0 = extent.x0 ;
      args.x1 = extent.x1 ;
      args.y0 = extent.y0 ;
      args.y1 = extent.y1 ;
      assert ( args.x0 < args.x1 ) ;
      assert ( args.y0 < args.y1 ) ;
      args.step = ( args.x1 - args.x0 ) / args.width ;
      args.twine = 0 ;
    }

    // only now we calculate a twining 'spread' - if we run a sequence
    // file, the twining parameters need to be recalculated for every
    // output image to adapt to changes of the output hfov. A spread is
    // the generalized equivalent of a filter kernel. While a 'standard'
    // convolution kernel has pre-determined geometry (it's a matrix of
    // coefficients meant to be applied to an equally-shaped matrix of
    // data) - here we have three values for each coefficient: the
    // first two define the position of the look-up relative to the
    // 'central' position, and the third is the weight and corresponds
    // to a 'normal' convolution coefficient.

    make_spread ( args.twine_spread , args.twine , args.twine ,
                  args.twine_width , args.twine_sigma ,
                  args.twine_threshold ) ;

    // find the parameters which are type-relevant to route to
    // the specialized code above. There are several stages of
    // 'roll_out' which move the parameterization given by run-time
    // arguments into type information.

    int nch = args.nchannels ;
    int ninp = ( args.itp == 1 ) ? 3 : 9 ;
    projection_t prj = args.projection ;

    // now call the ISA-specific rendering code via the dispatch_base
    // pointer received from get_dispatch.

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
