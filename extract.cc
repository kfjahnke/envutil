/************************************************************************/
/*                                                                      */
/*    extract - extract a partial image from an environment             */
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
// extent. For CL arguments, try 'extract --help'.
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
// binning. The default with this method is to use a simple 2X2 box filter
// on a signal which is oversampled by a factor of four. Additional
// parameters can change the smaount of oversampling and add gaussian
// weights to the filter parameters. 'twining' is quite fast (if the number
// of filter taps isn't very large - when down-scaling, the parameter
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

  // eval applies the quaternion

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

// from the envutil project, we use geometry.h for the coordinate
// transformations needed when dealing with environment images,
// especially cubemaps.

#include "geometry.h"

// helper function to save a zimt array of pixels to an image file

template < std::size_t nchannels >
void save_array ( const std::string & filename ,
                  const zimt::view_t
                    < 2 ,
                      zimt::xel_t < float , nchannels >
                    > & pixels ,
                  bool is_latlon = false )
{
  auto out = ImageOutput::create ( filename );
  assert ( out != nullptr ) ;
  ImageSpec ospec ( pixels.shape[0] , pixels.shape[1] ,
                    nchannels , TypeDesc::HALF ) ;
  out->open ( filename , ospec ) ;

  auto success = out->write_image ( TypeDesc::FLOAT , pixels.data() ) ;
  assert ( success ) ;
  out->close();
}

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
  double yaw , pitch , roll ;
  projection_t projection , env_projection ;
  std::string prj_str ;
  double x0 , x1 , y0 , y1 ;
  double hfov ;
  double step , env_step ;
  std::size_t width , height ;
  std::size_t env_width , env_height ;
  std::string input ;
  std::string output ;
  std::string seqfile ;
  std::size_t nchannels ;
  int itp , twine  ;
  double twine_width , twine_sigma , twine_threshold ;
  std::string swrap, twrap, mip, interp , tsoptions ;
  float stwidth , stblur ;
  bool conservative_filter ;
  std::unique_ptr<ImageInput> inp ;

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
    ap.intro("extract: extract image from an environment\n")
      .usage("extract [options] --input INPUT --output OUTPUT");
    ap.arg("-v", &verbose)
      .help("Verbose output");
    ap.arg("--input INPUT")
      .help("input file name (mandatory)")
      .metavar("INPUT");
    ap.arg("--output OUTPUT")
      .help("output file name (mandatory)")
      .metavar("OUTPUT");
    ap.arg("--seqfile SEQFILE")
      .help("image sequence file name (optional)")
      .metavar("OUTPUT");
    ap.arg("--itp ITP")
      .help("interpolator: 1 for bilinear, -1 for OIIO, -2 bilinear+twining")
      .metavar("ITP");
    ap.arg("--width EXTENT")
      .help("width of the output")
      .metavar("EXTENT");
    ap.arg("--height EXTENT")
      .help("height of the output")
      .metavar("EXTENT");
    ap.arg("--projection PRJ")
      .help("target projection")
      .metavar("PRJ");
    ap.arg("--hfov ANGLE")
      .help("horiziontal field of view of the output (in degrees)")
      .metavar("ANGLE");
    ap.arg("--yaw ANGLE")
      .help("yaw of the virtual camera (in degrees)")
      .metavar("ANGLE");
    ap.arg("--pitch ANGLE")
      .help("pitch of the virtual camera (in degrees)")
      .metavar("ANGLE");
    ap.arg("--roll ANGLE")
      .help("roll of the virtual camera (in degrees)")
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
    ap.arg("--twine TWINE")
      .help("use twine*twine oversampling - use with itp -2")
      .metavar("TWINE");
    ap.arg("--twine_width TWINE_WIDTH")
      .help("widen the pick-up area of the twining filter")
      .metavar("TWINE_WIDTH");
    ap.arg("--twine_sigma TWINE_SIGMA")
      .help("use a truncated gaussian for the twining filter (default: don't)")
      .metavar("TWINE_SIGMA");
    ap.arg("--twine_threshold TWINE_THRESHOLD")
      .help("discard twining filter taps below this threshold")
      .metavar("TWINE_THRESHOLD");
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
    ap.arg("--conservative_filter YESNO")
      .help("OIIO conservative_filter Texture Option - pass 0 or 1")
      .metavar("YESNO");
    
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
    itp = ap["itp"].get<int>(1);
    twine = ap["twine"].get<int>(0);
    twine_width = ap["twine_width"].get<float>(1.0);
    twine_sigma = ap["twine_sigma"].get<float>(0.0);
    twine_threshold = ap["twine_threshold"].get<float>(0.0);
    swrap = ap["swrap"].as_string ( "WrapDefault" ) ;
    twrap = ap["twrap"].as_string ( "WrapDefault" ) ;
    mip = ap["mip"].as_string ( "MipModeDefault" ) ;
    interp = ap["interp"].as_string ( "InterpSmartBicubic" ) ;
    tsoptions = ap["tsoptions"].as_string ( "automip=1" ) ;
    conservative_filter = ap["conservative_filter"].get<int>(1) ;
    x0 = ap["x0"].get<float> ( 0.0 ) ;
    x1 = ap["x1"].get<float> ( 0.0 ) ;
    y0 = ap["y0"].get<float> ( 0.0 ) ;
    y1 = ap["y1"].get<float> ( 0.0 ) ;
    width = ap["width"].get<int> ( 0 ) ;
    stwidth = ap["stwidth"].get<float> ( 1 ) ;
    stblur = ap["stblur"].get<float> ( 0 ) ;
    height = ap["height"].get<int> ( 0 ) ;
    hfov = ap["hfov"].get<float>(0.0);
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

    assert ( input != std::string() ) ;
    assert ( output != std::string() || seqfile != std::string() ) ;
    assert ( projection != PRJ_NONE ) ;
    assert ( width > 0 ) ;
    assert ( height > 0 ) ;
    if ( projection == CUBEMAP )
      assert ( height = 6 * width ) ;

    // some member variables in the args object are gleaned from
    // the input image:

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
      env_projection = CUBEMAP ;
      env_step = M_PI_2 / env_width ;
    }
    else
    {
      std::cerr << "input image must have 2:1 or 1:6 aspect ratio"
                << std::endl ;
      exit ( -1 ) ;
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

// environment.h has most of the 'workhorse' code for this program.
// it provides functional constructs to yield pixels for coordinates.
// These functors are used with zimt::process to populate the output.

#include "environment.h"

// because we specialize for different channel counts with a template
// argument (for maximum efficiency), the 'workhorse' code is in a
// function template 'work'. There are two overloads:

// 'work' overload to produce source data from a lat/lon environment
// image using OIIO's 'environment' or 'texture' function, or bilinear
// interpolation with 'twining'. All these lookup methods use 'ninepacks'
// which are used to glean the derivatives of the coordinate
// transformation, on top of the actual pick-up coordinate.

template < std::size_t nchannels >
void work ( const zimt::grok_get_t < float , 9 , 2 , 16 > & get_ray )
{
  // the 'act' functor is fed 'ninepacks' - sets of three 3D ray
  // coordinates, holding information to allow the calculations
  // of the coordinate transformation's derivatives. This is
  // done inside the 'env' object and depends on it's specific
  // type - the 'act' functor yields pixels, so here we don't
  // have to deal with it's internal workings.

  typedef zimt::xel_t < float , 9 > crd9_t ;
  typedef zimt::xel_t < float , nchannels > px_t ;

  // the 'act' functor's 'inner' type is variable, so we use a
  // uniform 'outer' type by 'grokking' it (type erasure)

  zimt::grok_type < crd9_t , px_t , 16 > act ;

  if ( args.itp == -2 )
  {
    // use twining (--itp -2). first create an 'environment' object.
    // Note how this object will persist (we declare it static), so
    // it's only created right when control flow arrives here for
    // the very first time (parameterization is taken from 'args').
    // This is deliberate: subsequent invocations of 'work' are
    // supposed to use the same environment, so it would be
    // wasteful to set up a new one - it's an expensive asset.

    static environment < float , float , nchannels , 16 > env ;

    // set up the twining filter

    if ( args.twine == 0 )
    {
      // user has passed twine 0, or not passed anything - but itp
      // is -2. We set up the twining with automatically generated
      // parameters.

      // figure out the magnification in the image center as a
      // guideline

      double mag = args.env_step / args.step ;

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

        args.twine = std::min ( 5 , int ( 1.0 + mag ) ) ;
        args.twine_width = mag ;
      }
      else
      {
        // otherwise, we use a twine size which depends on the
        // downscaling factor (reciprocal of 'mag') and a twine_width
        // of 1.0: we only want anti-aliasing. picking a sufficiently
        // large twine value guarantees that we're not skipping any
        // pixels in the source (due to undersampling).

        args.twine = int ( 1.0 + 1.0 / mag ) ;
        args.twine_width = 1.0 ;
      }

      if ( verbose )
      {
        std::cout << "automatic twining for magnification " << mag
                  << ":" << std::endl ;
        std::cout << "twine: " << args.twine
                  << " twine_width: " << args.twine_width
                  << std::endl ;
      }
    }
    else
    {
      if ( verbose )
      {
        std::cout << "using fixed twine: " << args.twine
                  << " twine_width: " << args.twine_width
                  << std::endl ;
      }
    }

    // with the given twining parameters, we can now set up a 'spread':
    // the generalized equivalent of a filter kernel. While a 'standard'
    // convolution kernel has pre-determined geometry (it's a matrix of
    // coefficients meant to be applied to an equally-shaped matrix of
    // data) - here we have three values for each coefficient: the
    // first two define the position of the look-up relative to the
    // 'central' position, and the third is the weight and corresponds
    // to a 'normal' convolution coefficient.

    std::vector < zimt::xel_t < float , 3 > > spread ;
    make_spread ( spread , args.twine , args.twine ,
                  args.twine_width , args.twine_sigma ,
                  args.twine_threshold ) ;

    // wrap the 'environment' object in a twine_t object and assign
    // to 'act' - this 'groks' the twine_t object to act's type.
    // Note how this object is not declared static: the twining
    // parameters may change due to changing hfov from one invocation
    // of 'work' to the next.

    act = twine_t < nchannels , 16 > ( env , spread ) ;
  }
  else
  {
    // use OIIO's 'environment' or 'texture' lookup functions.
    // These code paths are coded in the 'latlon' and 'cubemap'
    // objects, which are created and also grokked to 'act'.
    // There, the 'ninepack' is used to calculate the derivatives
    // and then all the data needed for the look-up are passed to
    // the relevant (batched) OIIO look-up function.

    if ( args.env_width == args.env_height * 2 )
    {
      // set up an environment object picking up pixel values from
      // a lat/lon image using OIIO's 'environment' function. Again
      // we set up a static object, but it's assigned to 'act' in
      // every invocation of 'work'.

      static latlon < float , float , nchannels , 16 > ll ;
      act = ll ;
    }
    else
    {
      // set up an environment object picking up pixel values from
      // a texture representing the cubemap image, using OIIO's
      // 'texture' function

      static cubemap < float , float , nchannels , 16 > cbm ;
      act = cbm ;
    }
  }

  // we have the 'act' functor set up. now we set up an array to
  // receive the output pixels. Again we use a static object:
  // the output size will remain the same, so the array can be
  // re-used every time and we don't have to deallocate and then
  // reallocate the memory.

  static zimt::array_t < 2 , px_t > trg ( { args.width , args.height } ) ;
  
  // set up a zimt::storer to populate the target array with
  // zimt::process. This is the third component needed for
  // zimt::process - we already have the get_t and act.
  
  zimt::storer < float , nchannels , 2 , 16 > cstor ( trg ) ;
  
  // use the get, act and put components with zimt::process
  // to produce the target images. This is the point where all
  // the state we have built up is finally put to use, running
  // a multithreaded pipeline which fills the target image.
  
  zimt::process ( trg.shape , get_ray , act , cstor ) ;
  
  if ( args.verbose )
    std::cout << "saving output image: " << args.output << std::endl ;

  // store the result to disk

  save_array < nchannels > ( args.output , trg ) ;
}

// overload using an 'environment' object directly as data source.
// This uses bilinear interplation directly from the source image.
// Note how the 'get_ray' object has only three input channels,
// in contrast to the previous overload where it has nine. This
// is the 'fast lane'.

template < std::size_t nchannels >
void work ( const zimt::grok_get_t < float , 3 , 2 , 16 > & get_ray )
{
  typedef zimt::xel_t < float , 3 > crd3_t ;
  typedef zimt::xel_t < float , nchannels > px_t ;

  // this is to accomodate the 'act' functor

  zimt::grok_type < crd3_t , px_t , 16 > act ;

  // create the 'environment' object and 'grok' it to 'act'

  static environment < float , float , nchannels , 16 > env ;
  act = env ;

  // set up an array to receive the output pixels

  static zimt::array_t < 2 , px_t > trg ( { args.width , args.height } ) ;
  
  // set up a zimt::storer to populate the target array with
  // zimt::process
  
  zimt::storer < float , nchannels , 2 , 16 > cstor ( trg ) ;
  
  // use the get, act and put components with zimt::process
  // to produce the target images and store them to disk
  
  zimt::process ( trg.shape , get_ray , act , cstor ) ;

  // save the result

  if ( args.verbose )
    std::cout << "saving output image: " << args.output << std::endl ;
  
  save_array < nchannels > ( args.output , trg ) ;
}

// to call the apprpriate template instantiation and overload of
// 'work' we need to do some dispatching: picking types depending
// on run-time variables. We achieve this with a staged 'dispatch'
// routine: every stage processes one argument and dispatches to
// specialized code. The least specialized dispatch variant is
// the lowest one down, this here is the final stage where we have
// the number of channels and the stepper as template arguments.
// Here we proceed to set up more state which is common to all
// code paths and finally call 'work' to run the pixel pipeline.

template < int NCH , typename stepper_t >
void dispatch()
{
  // first do some processing here which is independent from the
  // number of channels in the pixels:

  bool have_seq = ( args.seqfile != std::string() ) ;
  std::ifstream seqstream ;
  if ( have_seq )
  {
    seqstream.open ( args.seqfile ) ;
    assert ( seqstream.good() ) ;
  }

  int seqno = 0 ;
  while ( true )
  {
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
      ++seqno ;
      char buffer [ 12 ] ;
      snprintf ( buffer , 12 , "%03d" , seqno ) ;
      args.output = "seq" + std::string(buffer) + ".jpg" ;
      args.twine = 0 ;
    }

    // orthonormal system of basis vectors for the view

    crd3_t xx { 1.0 , 0.0 , 0.0 } ;
    crd3_t yy { 0.0 , 1.0 , 0.0 } ;
    crd3_t zz { 0.0 , 0.0 , 1.0 } ;

    // the three vectors are rotated with the given yaw, pitch and
    // roll, and later passed on to the to 'steppers', the objects
    // which provide 3D 'ray' coordinates. They incorporate the
    // rotated basis in their ray generation, resulting in
    // appropriately oriented ray coordinates which can be formed
    // more efficiently in the steppers - first calculating the
    // rays and then rotating the rays in a second step takes
    // more CPU cycles.

    rotate_3d < float , 16 > r3 ( args.roll , args.pitch , args.yaw ) ;
    
    xx = r3 ( xx ) ;
    yy = r3 ( yy ) ;
    zz = r3 ( zz ) ;

    // set up the stepper. note the extents of the 2D manifold
    // given in model space uits. These objects will deliver 3D
    // 'ray' coordinates as input to the 'act' functor, or 9D
    // ray and neighbouring rays (deriv_stepper). Note how
    // each projection has it's distinct type of stepper, but
    // they are all assigned to a common type, a 'grok_get_t'.
    // This uses type erasure and captures the functionality
    // in std::functions, and the resulting object is only
    // characterized by it's input and output type and lane
    // count. The uniform type allows us to pass these objects
    // around with a common type - a convenient way of harnessing
    // groups of types with different implementation but equal
    // interface.

    stepper_t get_ray ( xx , yy , zz , args.width , args.height ,
                        args.x0 , args.x1 , args.y0 , args.y1 ) ;

    zimt::grok_get_t < float , stepper_t::size ,  2 , 16 >
          getter ( get_ray ) ;

    // now we can call the channel-specific 'work', with the
    // stepper-specific getter. There are two overloads:
    // one taking the getters yielding simple single-ray
    // coordinates and one taking the three-ray variant
    // needed to compute the derivatives.

    work < NCH > ( getter ) ;

    // if we're not working a sequence, we're done now.

    if ( ! have_seq )
      break ;
  }
}

// we already have the number of colour channels and the type of the
// stepper as template arguments, now we dispatch depending on
// 'ninputs': three means it's an ordinary lookup using the stepper
// directly, and nine means it's a lookup using ninepacks and
// working with the derivatives.

template < int NCH ,
           template < typename , std::size_t > class STP >
void dispatch ( int ninputs )
{
  switch ( ninputs )
  {
    case 3:
      dispatch < NCH , STP < float , 16 > >() ;
      break ;
    case 9:
      dispatch < NCH , deriv_stepper < float , 16 , STP > >() ;
      break ;
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

  // find the parameters which are type-relevant to dispatch to
  // the specialized code

  int nch = args.nchannels ;
  int ninp = ( args.itp == 1 ) ? 3 : 9 ;
  projection_t prj = args.projection ;

  dispatch ( nch , ninp , prj ) ;
}
    
