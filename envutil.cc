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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C"
{
  #include <libavcodec/avcodec.h>
  #include <libavutil/opt.h>
  #include <libavutil/imgutils.h>
} ;

#include <fstream>

#include "stepper.h"

// To conveniently rotate with a rotational quaternion, we employ
// Imath's 'Quat' data type, packaged in a zimt::unary_functor.
// This is not factored out because it requires inclusion of
// some Imath headers, which I want to keep out of the other
// code, e.g. in geometry.h, where it would fit in nicely.

#include <Imath/ImathVec.h>
#include <Imath/ImathEuler.h>
#include <Imath/ImathQuat.h>
#include <Imath/ImathLine.h>

// rotate_3d uses a SIMDized Imath Quaternion to affect a 3D rotation
// of a 3D SIMDized coordinate. Imath::Quat<float> can't broadcast
// to handle SIMDized input, but if we use an Imath::Quat of the
// SIMDized type, we get the desired effect.

template < typename T , std::size_t L >
struct rotate_3d
: public zimt::unary_functor
    < zimt::xel_t < T , 3 > , zimt::xel_t < T , 3 > , L >
{
  typedef zimt::simdized_type < T , L > f_v ;
  typedef zimt::xel_t < T , 3 > crd3_t ;
  typedef zimt::simdized_type < crd3_t , L > crd3_v ;

  Imath::Quat < T > q ;

  rotate_3d ( T roll , T pitch , T yaw , bool inverse = false )
  {
    // set up the rotational quaternion. if 'inverse' is set, produce
    // the conjugate.

    if ( inverse )
    {
      Imath::Eulerf angles ( -yaw , -pitch , -roll , Imath::Eulerf::YXZ ) ;
      q = angles.toQuat() ;
    }
    else
    {
      Imath::Eulerf angles ( roll , pitch , yaw , Imath::Eulerf::ZXY ) ;
      q = angles.toQuat() ;
    }
  }

  // eval applies the quaternion, reinterpret-casting the zimt 'xel'
  // data we use throughoutthe program to Imath::Vec3 - both types
  // are compatible.

  template < typename U >
  void eval ( const zimt::xel_t < U , 3 > & in ,
              zimt::xel_t < U , 3 > & out ) const
  {
    auto const & in_e
      = reinterpret_cast < const Imath::Vec3 < U > & > ( in ) ;

    auto & out_e
      = reinterpret_cast < Imath::Vec3 < U > & > ( out ) ;

    out_e = in_e * Imath::Quat < U > ( q ) ;
  }

  // for convenience:

  template < typename U >
  zimt::xel_t < U , 3 > operator() ( const zimt::xel_t < U , 3 > & in )
  {
    zimt::xel_t < U , 3 > out ;
    eval ( in , out ) ;
    return out ;
  }
} ;

// for image I/O, we use OpenImageIO.

#include <OpenImageIO/imageio.h>

using OIIO::ImageInput ;
using OIIO::ImageOutput ;
using OIIO::TypeDesc ;
using OIIO::ImageSpec ;

// struct image_series holds a format string for a series of numbered
// images. It's a standard printf-type format string containing precisely
// one %-sequence accepting an integer - e.g. %04d
// It has operator[] to provide the filename fixed with the given index.

struct image_series
{
  std::string format_string ;
  std::size_t buffer_size ;

  image_series ( const std::string & _format_string )
  : format_string ( _format_string )
  {
    // scan the format string for percent signs. Only one of them
    // is allowed.

    int percent_count = 0 ;
    int percent_pos = 0 ;
    int pos = 0 ;

    for ( auto const & c : format_string )
    {
      if ( c == '%' )
      {
        percent_count++ ;
        percent_pos = pos ;
      }
      ++pos ;
    }

    if ( percent_count == 1 )
      buffer_size = pos + 16 ;
  }

  bool valid()
  {
    return ( buffer_size != 0 ) ;
  }

  std::string operator[] ( const std::size_t & index )
  {
    if ( buffer_size )
    {
      char buffer [ buffer_size ] ;
      snprintf ( buffer , buffer_size ,
                 format_string.c_str() , int(index) ) ;
      return buffer ;
    }
    else
    {
      return format_string ;
    }
  }
} ;

// similar class, used for six cube faces

struct cubeface_series
{
  std::string format_string ;
  std::vector < std::string > filename ;
  std::size_t buffer_size ;

  cubeface_series() = default ;

  cubeface_series ( const std::string & _format_string )
  : format_string ( _format_string )
  {
    // scan the format string for percent signs. Only one of them
    // is allowed.

    int percent_count = 0 ;
    int percent_pos = 0 ;
    int pos = 0 ;

    for ( auto const & c : format_string )
    {
      if ( c == '%' )
      {
        percent_count++ ;
        percent_pos = pos ;
      }
      ++pos ;
    }

    buffer_size = 0 ;

    if ( percent_count == 1 )
    {
      buffer_size = pos + 16 ;

      char buffer [ buffer_size ] ;
      snprintf ( buffer , buffer_size , format_string.c_str() , "left" ) ;
      filename.push_back ( buffer ) ;
      snprintf ( buffer , buffer_size , format_string.c_str() , "right" ) ;
      filename.push_back ( buffer ) ;
      snprintf ( buffer , buffer_size , format_string.c_str() , "top" ) ;
      filename.push_back ( buffer ) ;
      snprintf ( buffer , buffer_size , format_string.c_str() , "bottom" ) ;
      filename.push_back ( buffer ) ;
      snprintf ( buffer , buffer_size , format_string.c_str() , "front" ) ;
      filename.push_back ( buffer ) ;
      snprintf ( buffer , buffer_size , format_string.c_str() , "back" ) ;
      filename.push_back ( buffer ) ;
    }
  }

  bool valid()
  {
    return filename.size() == 6 ;
  }

  // we provide two operator[] overloads. The first one accesses the
  // cubeface filename by number, the second inserts a given orientation
  // string.

  std::string operator[] ( const std::size_t & index )
  {
    assert ( buffer_size != 0 && index < 6 ) ;
    return filename [ index ] ;
  }

  std::string operator[] ( const std::string & face )
  {
    char buffer [ buffer_size ] ;
    snprintf ( buffer , buffer_size ,
               format_string.c_str() ,
               face.c_str() ) ;
    return buffer ;
  }

  // return a reference to the set of six names

  const std::vector < std::string > & get_filenames()
  {
    return filename ;
  }
} ;

// geometry.h has the coordinate transformations needed when
// dealing with environment images, especially cubemaps.

#include "geometry.h"

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

#include <regex>
#include <OpenImageIO/filesystem.h>
#include <OpenImageIO/argparse.h>

using OIIO::ArgParse ;
using OIIO::Filesystem::convert_native_arguments ;

struct arguments
{
  bool verbose ;
  std::string input ;
  std::string output ;
  std::string mount_image ;
  double hfov ;
  float mount_hfov ;
  std::size_t width ;
  std::size_t height ;
  std::size_t mount_width ;
  std::size_t mount_height ;
  std::string prj_str ;
  std::string mount_prj_str ;
  projection_t projection ;
  projection_t mount_prj;

  double cbmfov ;
  std::size_t support_min ;
  std::size_t tile_size ;
  bool ctc ;

  double yaw , pitch , roll ;
  double x0 , x1 , y0 , y1 ;
 
  std::string seqfile ;
  std::string codec ;
  float mbps ;
  int fps ;

  int itp ;
  int twine  ;
  std::string twf_file ;
  bool twine_normalize ;
  bool twine_precise ;
  double twine_width , twine_density , twine_sigma , twine_threshold ;
  std::string swrap, twrap, mip, interp , tsoptions ;
  float stwidth , stblur ;
  bool conservative_filter ;

  // gleaned from other parameters or input images

  bool multiple_input ;
  cubeface_series cfs ;
  projection_t env_projection ;
  double step , env_step ;
  std::size_t env_width , env_height ;
  std::size_t nchannels ;
  std::unique_ptr<ImageInput> inp ;
  std::unique_ptr<ImageInput> mount_inp ;

  // technical variables for the argument parser

  std::string metamatch ;
  std::regex field_re ;

  // the 'arguments' object's 'init' takes the main program's argc
  // and argv.

  void init ( int argc , const char ** argv )
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
      .help("input file name (mandatory)")
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
      .help("interpolator: 1 for bilinear, -1 for OIIO, -2 bilinear+twining")
      .metavar("ITP");

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

    ap.separator("  parameters for mounted image input:");
    
    // std::string mount_image , mount_prj ;
    // float mount_hfov ;
    ap.add_argument("--mount %s:IMAGE %s:PROJECTION %f:HFOV",
                    &mount_image , &mount_prj_str, &mount_hfov)
      .help("load non-environment source image") ;

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
    codec = ap["codec"].as_string ( "libx265" ) ; 
    mbps = ( 1000000.0 * ap["mbps"].get<float> ( 8.0 ) ) ;
    fps = ap["fps"].get<int>(60);
    itp = ap["itp"].get<int>(1);
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
    int prj = 0 ;
    for ( const auto & p : projection_name )
    {
      if ( p == prj_str )
        break ;
      ++ prj ;
    }
    projection = projection_t ( prj ) ;
    if ( mount_image != std::string() )
    {
      prj = 0 ;
      for ( const auto & p : projection_name )
      {
        if ( p == mount_prj_str )
          break ;
        ++ prj ;
      }
      mount_prj = projection_t ( prj ) ;
      assert ( input == std::string() ) ;
    }
    else
    {
      assert ( input != std::string() ) ;
    }
    assert ( output != std::string() ) ;
    if ( width == 0 )
      width = 1024 ;
    if ( projection == CUBEMAP )
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

    if ( mount_image != std::string() )
    {
      std::cout << "create env from mount image " << mount_image << std::endl ;
      mount_inp = ImageInput::open ( mount_image ) ;
      assert ( mount_inp ) ;

      const ImageSpec &spec = mount_inp->spec() ;

      mount_width = spec.width ;
      mount_height = spec.height ;
      nchannels = spec.nchannels ;
      mount_hfov *= M_PI / 180.0 ;
      switch ( mount_prj )
      {
        case SPHERICAL:
        case FISHEYE:
        case CYLINDRICAL:
          env_step = mount_hfov / mount_width ;
          break ;
        case RECTILINEAR:
          env_step = 2.0 * tan ( mount_hfov / 2.0 ) / mount_width ;
          break ;
        case STEREOGRAPHIC:
          env_step = 4.0 * tan ( mount_hfov / 4.0 ) / mount_width ;
          break ;
        default:
          break ;
      }
    }
    else
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
          env_projection = CUBEMAP ;
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
          env_projection = CUBEMAP ;
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
        std::cout << "interpolation: "
                << ( itp == 1 ? "direct bilinear" : "uses OIIO" )
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

      if ( ( projection == CUBEMAP ) && ctc )
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
  }
} ;

// to avoid having to pass the arguments around, we us a global
// 'args' object.

arguments args ;

// helper function to save a zimt array of pixels to an image file, or
// to a set of six cube face images, if 'output' has a format string.

template < std::size_t nchannels >
void save_array ( const std::string & filename ,
                  const zimt::view_t
                    < 2 ,
                      zimt::xel_t < float , nchannels >
                    > & pixels ,
                  bool is_latlon = false )
{
  if ( args.projection == CUBEMAP )
  {
    // output is a cubemap, let's see if 'output' is a format
    // string for six separate cube faces

    assert ( pixels.shape[1] == 6 * pixels.shape[0] ) ;

    auto has_percent = filename.find_first_of ( "%" ) ;
    if ( has_percent != std::string::npos )
    {
      // input must be a set of six cubeface images, that's the
      // only way how we accept a format string.

      cubeface_series cfs ( filename ) ;
      if ( cfs.valid() )
      {
        // we'll call save_array recursively for the single images, hence:

        args.projection = RECTILINEAR ;

        // save six subarrays to individual images

        std::size_t w = pixels.shape[0] ;
        for ( std::size_t i = 0 ; i < 6 ; i++ )
        {
          save_array ( cfs[i] ,
                       pixels.window ( { 0L , long ( i * w ) } ,
                                       { long ( w ) , long ( ( i + 1 ) * w ) } ) ) ;
        }

        // restore the projection in 'args'

        args.projection = CUBEMAP ;

        // we're done.

        return ;
      }
    }
  }

  // if we land here, we're supposed to store an ordinary single image

  auto out = ImageOutput::create ( filename );
  assert ( out != nullptr ) ;
  ImageSpec ospec ( pixels.shape[0] , pixels.shape[1] ,
                    nchannels , TypeDesc::HALF ) ;
  out->open ( filename , ospec ) ;

  auto success = out->write_image ( TypeDesc::FLOAT , pixels.data() ) ;
  assert ( success ) ;
  out->close();
}

// quick shot at encoding video: I'm using an example file from the ffmpeg
// examples section, with some adaptations to bend the C code to C++,
// and additions to convert the float RGB data from the pixel pieline 
// to YUV for encoding with h264. C++ style comments (begining with //)
// are mine, to explain where I made alterations of the original C
// code or added stuff - and why.

// The code example's end is marked with 'end of copied ffmpeg example
// code'

// Begin copied example code:

/*
 * Copyright (c) 2001 Fabrice Bellard
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/**
 * @file libavcodec encoding video API usage example
 * @example encode_video.c
 *
 * Generate synthetic video data and encode it to an output file.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C"
{
#include <libavcodec/avcodec.h>
#include "libavformat/avformat.h"
#include <libavutil/opt.h>
#include <libavutil/imgutils.h>
#include <libswscale/swscale.h>
} ;

static void encode(AVCodecContext *enc_ctx, AVFrame *frame, AVPacket *pkt,
                   FILE *outfile)
{
    int ret;

    /* send the frame to the encoder */

    // the format strings used for the printf don't seem to work with
    // C++, so I commented them out:
  
    // if (frame)
    //     printf("Send frame %3"PRId64"\n", frame->pts);

    ret = avcodec_send_frame(enc_ctx, frame);
    if (ret < 0) {
        fprintf(stderr, "Error sending a frame for encoding\n");
        exit(1);
    }

    while (ret >= 0) {
        ret = avcodec_receive_packet(enc_ctx, pkt);
        if (ret == AVERROR(EAGAIN) || ret == AVERROR_EOF)
            return;
        else if (ret < 0) {
            fprintf(stderr, "Error during encoding\n");
            exit(1);
        }

        // printf("Write packet %3"PRId64" (size=%5d)\n", pkt->pts, pkt->size);
        fwrite(pkt->data, 1, pkt->size, outfile);
        av_packet_unref(pkt);
    }
}

// what was in the original example is now the body of a C++ object,
// and the C code's global variables are now class members. main() is
// replaced by the class' c'tor which sets up the object for sustained
// operation with the given parameterization, and the code to push a
// single frame to the codec is in a member function and will be called
// repeatedly by the process generating images from a sequence file.

struct frame_sink
{
  const char *filename, *codec_name;
  const AVCodec *codec;
  AVCodecContext *c= NULL;
  int i, ret, x, y;
  FILE *f;
  AVFrame *frame;
  AVPacket *pkt;
  uint8_t endcode[4] { 0, 0, 1, 0xb7 };

  frame_sink ( const char * _filename ,
               const char * _codec_name )
  : filename ( _filename ) ,
    codec_name ( _codec_name )
  {
    /* find the mpeg1video encoder */
    codec = avcodec_find_encoder_by_name(codec_name);
    if (!codec) {
        fprintf(stderr, "Codec '%s' not found\n", codec_name);
        exit(1);
    }

    c = avcodec_alloc_context3(codec);
    if (!c) {
        fprintf(stderr, "Could not allocate video codec context\n");
        exit(1);
    }

    pkt = av_packet_alloc();
    if (!pkt)
        exit(1);

    /* put sample parameters */
    c->bit_rate = args.mbps ;
    /* resolution must be a multiple of two */
    c->width = args.width ;
    c->height = args.height ;
    /* frames per second */
    c->time_base = (AVRational){1, args.fps};
    c->framerate = (AVRational){args.fps, 1};

    /* emit one intra frame every ten frames
     * check frame pict_type before passing frame
     * to encoder, if frame->pict_type is AV_PICTURE_TYPE_I
     * then gop_size is ignored and the output of encoder
     * will always be I frame irrespective to gop_size
     */
    c->gop_size = 10;
    c->max_b_frames = 1;

    // for now, there is no choice about the codec, it's hardcoded, and
    // this codec seems to require YUV420P pixels. H264 works with the
    // same code and is faster.

    c->pix_fmt = AV_PIX_FMT_YUV420P;

    if (codec->id == AV_CODEC_ID_H265)
        av_opt_set(c->priv_data, "preset", "slow", 0);

    /* open it */
    ret = avcodec_open2(c, codec, NULL);
    if (ret < 0) {
        fprintf(stderr, "Could not open codec: %s\n", av_err2str(ret));
        exit(1);
    }

    f = fopen(filename, "wb");
    if (!f) {
        fprintf(stderr, "Could not open %s\n", filename);
        exit(1);
    }

    frame = av_frame_alloc();
    if (!frame) {
        fprintf(stderr, "Could not allocate video frame\n");
        exit(1);
    }
    frame->format = c->pix_fmt;
    frame->width  = c->width;
    frame->height = c->height;

    ret = av_frame_get_buffer(frame, 0);
    if (ret < 0) {
        fprintf(stderr, "Could not allocate the video frame data\n");
        exit(1);
    }

  }

  template < std::size_t nchannels >
  void encode_frame ( const zimt::view_t
                        < 2 ,
                          zimt::xel_t < float , nchannels >
                        > & pixels )
  {
    static struct SwsContext *img_convert_ctx;
    static int sws_flags = SWS_BICUBIC;
    if (img_convert_ctx == NULL)
    {
      // here we set up the converter to move from RGB to YUV. This is an
      // adaptation from another ffmpeg example file by the same author,
      // licensed as the other ffmpeg example code:

      /*
       * Libavformat API example: Output a media file in any supported
       * libavformat format. The default codecs are used.
       *
       * Copyright (c) 2003 Fabrice Bellard
       * ... license follows, see above */

      img_convert_ctx = sws_getContext(c->width, c->height,
                        // we say we have incoming RGB24. 
                                       AV_PIX_FMT_RGB24,
                                       c->width, c->height,
                        // and we want outgoing YUV420P
                                       AV_PIX_FMT_YUV420P,
                                       sws_flags,
                                       NULL, NULL, NULL );
    }
    if (img_convert_ctx == NULL)
    {
      fprintf(stderr, "Cannot initialize the conversion context\n");
      exit(1);
    }

    // end of code adapted from the second ffmpeg example, back to the
    // first example

    static int i = 0;

    // instead of having this code in a loop, we now have it in a
    // separate function which is called once for each frame.

    fflush(stdout);

    /* Make sure the frame data is writable.
        On the first round, the frame is fresh from av_frame_get_buffer()
        and therefore we know it is writable.
        But on the next rounds, encode() will have called
        avcodec_send_frame(), and the codec may have kept a reference to
        the frame in its internal structures, that makes the frame
        unwritable.
        av_frame_make_writable() checks that and allocates a new buffer
        for the frame only if necessary.
      */
    ret = av_frame_make_writable(frame);
    if (ret < 0)
        exit(1);

    // internally, we have the frame in float RGB(A). We only output
    // three channels. Alpha should be associated and result in darkening
    // the pixels without additional processing. If it were unassociated,
    // we'd have to multiply with alpha now. TODO: check
    // the float values are in the range [0,1], just as OIIO initially
    // produced them from whichever input file was given, so now we
    // have to scale up to [0,255]. An alternative would be to hold
    // a zimt array of YUV values and chain the conversion to the pixel
    // pipeline in the 'act' functor. But this conversion isn't the most
    // time-consuming operation, it's running the H265 codec (at least
    // here - my GPU doesn't have it).
    // Note that the code seems to expect sRGB data. There is no provision
    // for processing linear RGB yet, this should be added later, plus,
    // eventually, handling of other colour spaces.

    unsigned char buffer [ args.width * args.height * 3 ] ;
    for ( int y = 0 ; y < args.height ; y++ )
    {
      for ( int x = 0 ; x < args.width ; x++ )
      {
        const auto & px ( pixels [ { x , y } ] ) ;
        buffer [ 3 * y * args.width + 3 * x ] = 255.0f * px[0] ;
        if constexpr ( nchannels > 1 )
          buffer [ 3 * y * args.width + 3 * x + 1 ] = 255.0f * px[1] ;
        else
          buffer [ 3 * y * args.width + 3 * x + 1 ] = 255.0f * px[0] ;
        if constexpr ( nchannels > 2 )
          buffer [ 3 * y * args.width + 3 * x + 2 ] = 255.0f * px[2] ;
        else
          buffer [ 3 * y * args.width + 3 * x + 2 ] = 255.0f * px[0] ;
      }
    }

    // we need to pass a pointer to a pointer, hence:

    unsigned char * p_buffer = buffer ;

    // w is in bytes, not pixels, hence:
  
    int w = args.width * 3 ;

    // now do the type conversion. we don't have a size change.

    sws_scale ( img_convert_ctx,
                &p_buffer, &w, 0, c->height,
                frame->data, frame->linesize); 

    frame->pts = i;

    /* encode the image */
    encode(c, frame, pkt, f);

    if ( args.verbose )
      std::cout << "encoded video frame " << i << std::endl ;

    ++i ;
  }

  // the d'tor has the code needed to shut everything down gracefully.

  ~frame_sink()
  {
    /* flush the encoder */
    encode(c, NULL, pkt, f);

    /* Add sequence end code to have a real MPEG file.
       It makes only sense because this tiny examples writes packets
       directly. This is called "elementary stream" and only works for some
       codecs. To create a valid file, you usually need to write packets
       into a proper file format or protocol; see mux.c.
     */
    if (codec->id == AV_CODEC_ID_MPEG1VIDEO || codec->id == AV_CODEC_ID_MPEG2VIDEO)
        fwrite(endcode, 1, sizeof(endcode), f);
    fclose(f);

    avcodec_free_context(&c);
    av_frame_free(&frame);
    av_packet_free(&pkt);
  }
} ;

// end of copied ffmpeg example code

// here's our interface to the modified ffmpeg example code:
// this function takes a frame's worth of data encoded as a
// zimt::view, sets up the frame_sink object when first called
// and uses it to encode one frame per call. Alternatively,
// if 'args.output' is parsed as a format string, individual
// frame images are produced.

template < std::size_t nchannels >
void push_video_frame ( const zimt::view_t
                        < 2 ,
                          zimt::xel_t < float , nchannels >
                        > & pixels )
{
  static image_series is ( args.output ) ;
  static int index = 0 ;
  if ( is.valid() )
  {
    save_array ( is [ index ] , pixels ) ;
    ++index ;
  }
  else
  {
    static frame_sink video ( args.output.c_str() , args.codec.c_str() ) ;
    video.encode_frame ( pixels ) ;
  }
}

// 'work' is the function where the state we've built up in the
// code below it culminates in the call to zimt::process, which
// obtains pick-up coordinates from the 'get' object, produces
// pixels for these coordinates with the 'act' object and stores
// these pixels with a zimt::storer, which is set up in this
// function. All the specificity of the code has now moved to
// the types get_t and act_t - the dispatch is complete and we
// can code a single 'work' template which is good for all
// variants.

template < typename get_t , typename act_t >
void work ( get_t & get , act_t & act )
{
  typedef typename act_t::out_type px_t ;
  static const size_t nchannels = px_t::size() ;

  // we have the functors set up. now we set up an array to
  // receive the output pixels. Again we use a static object:
  // the output size will remain the same, so the array can be
  // re-used every time and we don't have to deallocate and then
  // reallocate the memory.

  static zimt::array_t < 2 , px_t >
    trg ( { args.width , args.height } ) ;
  
  // set up a zimt::storer to populate the target array with
  // zimt::process. This is the third component needed for
  // zimt::process - we already have the get_t and act.
  
  zimt::storer < float , nchannels , 2 , 16 > cstor ( trg ) ;
  
  // use the get, act and put components with zimt::process
  // to produce the target images. This is the point where all
  // the state we have built up is finally put to use, running
  // a multithreaded pipeline which fills the target image.
  
  zimt::bill_t bill ;
  // bill.njobs = 1 ;
  zimt::process ( trg.shape , get , act , cstor , bill ) ;
  
  // store the result to disk - either as a single frame added
  // to the video file, or as a single image stored individually.

  if ( args.seqfile != std::string() )
  {
    push_video_frame ( trg ) ;
  }
  else
  {
    if ( args.verbose )
      std::cout << "saving output image: " << args.output << std::endl ;

    save_array < nchannels > ( args.output , trg ) ;
  }
}

// environment.h has most of the 'workhorse' code for this program.
// it provides functional constructs to yield pixels for coordinates.
// These functors are used with zimt::process to populate the output.

#include "environment.h"

// to call the appropriate instantiation of 'work' (above)
// we need to do some dispatching: picking types depending
// on run-time variables. We achieve this with a staged 'dispatch'
// routine: every stage processes one argument and dispatches to
// specialized code. The least specialized dispatch variant is
// the lowest one down, this here is the final stage where we have
// the number of channels and the stepper as template arguments.
// Here we proceed to set up more state which is common to all
// code paths, set up the 'act' functor which yields pixels, and
// finally call 'work' to run the pixel pipeline.

template < int NCH ,
           template < typename , std::size_t > class STP >
void dispatch ( int ninputs )
{
  // set up an orthonormal system of basis vectors for the view

  zimt::xel_t < double , 3 > xx { 1.0 , 0.0 , 0.0 } ;
  zimt::xel_t < double , 3 > yy { 0.0 , 1.0 , 0.0 } ;
  zimt::xel_t < double , 3 > zz { 0.0 , 0.0 , 1.0 } ;

  // the three vectors are rotated with the given yaw, pitch and
  // roll, and later passed on to the to 'steppers', the objects
  // which provide 3D 'ray' coordinates. They incorporate the
  // rotated basis in their ray generation, resulting in
  // appropriately oriented ray coordinates which can be formed
  // more efficiently in the steppers - first calculating the
  // rays and then rotating the rays in a second step takes
  // (many) more CPU cycles.

  rotate_3d < double , 16 > r3 ( args.roll , args.pitch , args.yaw ) ;

  xx = r3 ( xx ) ;
  yy = r3 ( yy ) ;
  zz = r3 ( zz ) ;

  // There are two code paths: one taking the getters yielding
  // simple single-ray 3D coordinates, and one taking the three-ray
  // variant needed to compute the derivatives. They use specific
  // types of 'environment' objects. The code for the environment
  // objects is in environment.h

  if ( ninputs == 3 )
  {
    // set up a simple single-coordinate stepper of the type
    // fixed by 'STP'. This route is taken with direct bilinear
    // interpolation (itp 1)

    STP < float , 16 > get_ray
      ( xx , yy , zz , args.width , args.height ,
        args.x0 , args.x1 , args.y0 , args.y1 ) ;

    // create a static 'environment' object. This will persist, so
    // if we're creating an image sequence, it will be reused for
    // each individual image.

    static environment < float , float , NCH , 16 > env ;

    // now we call the final 'work' template which uses the get_t
    // and act objects we've set up so far

    work ( get_ray , env ) ;
  }
  else // ninputs == 9
  {
    // do the same, but with a deriv_stepper and an 'environment9'
    // object. This code path is taken for lookup with 'ninepacks'
    // holding three rays - the additional two used to calculate
    // the derivatives.

    deriv_stepper < float , 16 , STP > get_ray
      ( xx , yy , zz , args.width , args.height ,
        args.x0 , args.x1 , args.y0 , args.y1 ) ;

    environment9 < NCH , 16 > env ;

    work ( get_ray , env ) ;
  }
}

// we have the number of channels as a template argument from the
// dispatch below, now we dispatch on the projection and instantiate
// the next dispatch level with a stepper type which fits the
// projection.

template < int NCH >
void dispatch ( int ninputs ,
                projection_t projection )
{
  switch ( projection )
  {
    case SPHERICAL:
      dispatch < NCH , spherical_stepper > ( ninputs ) ;
      break ;
    case CYLINDRICAL:
      dispatch < NCH , cylindrical_stepper > ( ninputs ) ;
      break ;
    case RECTILINEAR:
      dispatch < NCH , rectilinear_stepper > ( ninputs ) ;
      break ;
    case FISHEYE:
      dispatch < NCH , fisheye_stepper > ( ninputs ) ;
      break ;
    case STEREOGRAPHIC:
      dispatch < NCH , stereographic_stepper > ( ninputs ) ;
      break ;
    case CUBEMAP:
      dispatch < NCH , cubemap_stepper > ( ninputs ) ;
      break ;
    default:
      break ;
  }
}

// dispatch by the number of colour channels. We process one to four,
// where the usefulness of two channels isn't clear.

void dispatch ( int nchannels ,
                int ninputs ,
                projection_t projection )
{
  switch ( nchannels )
  {
    case 1:
      dispatch < 1 > ( ninputs , projection ) ;
      break ;
    case 2:
      dispatch < 2 > ( ninputs , projection ) ;
      break ;
    case 3:
      dispatch < 3 > ( ninputs , projection ) ;
      break ;
    case 4:
      dispatch < 4 > ( ninputs , projection ) ;
      break ;
  }
}

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

  // This loop will iterate over the renditions of single images.
  // If we're not running a sequence, there will only be one
  // iteration.

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

    // find the parameters which are type-relevant to dispatch to
    // the specialized code above. There are several stages of
    // 'dispatch' which move the parameterization given by run-time
    // arguments into type information.

    int nch = args.nchannels ;
    int ninp = ( args.itp == 1 ) ? 3 : 9 ;
    projection_t prj = args.projection ;

    dispatch ( nch , ninp , prj ) ;
  }
  while ( have_seq ) ; // loop criterion for 'do' loop
}
    
