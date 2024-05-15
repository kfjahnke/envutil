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

// utility to extract an image from an environment. This program takes
// a 2:1 lat/lon environment or a 1:6 cubemap image as input and produces
// output in the specified orientation, projection, field of view and
// extent. this started out as a simple demo for 'steppers' (stepper.cc
// is still there with the initial code), but I thought that, with a bit
// of additional parameterization, it would make a useful tool.
// For the time being, it only uses bilinear interpolation, so the
// resolution of the output should be close to the input's. For CL
// arguments, try 'extract -?'. The output projection can be one of
// "spherical", "cylindrical", "rectilinear", "stereographic",
// "fisheye" or "cubemap". The geometrical extent of the output is
// set up most conveniently by passing --hfov, the horizontal field
// of view of the output. The x0, x1, y0, and y1 parameters allow
// passing specific extent values (in model space units), which should
// rarely be necessary. To specify a 3D rotation, pass Euler angles
// yaw, pitch and roll - they default to zero: no rotation. The size
// of the output is given by --width and --height. You must pass an
// output filename with --output; --input specifies the environment
// image.

#include "stepper.h"

// To conveniently rotate with a rotational quaternion, we employ
// Imath's 'Quat' data type, packaged in a zimt::unary_functor.
// This is not factored out because it requires inclusion of
// some Imath headers, which I want to keep out of the other
// code, e.g. in geometry.h, where it would fit in nicely.

#include <Imath/ImathVec.h>
#include <Imath/ImathEuler.h>
#include <Imath/ImathQuat.h>

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

// for image I/O, we use OpenImageIO. For this demo, this is all
// the OpenImageIO code we'll use:

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

  if ( is_latlon )
    ospec.attribute ( "textureformat" , "LatLong Environment" ) ;

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
// provide more palatable values to the program.

#include <regex>
#include <OpenImageIO/filesystem.h>
#include <OpenImageIO/argparse.h>

using OIIO::ArgParse ;
using OIIO::Filesystem::convert_native_arguments ;

struct arguments
{
  bool verbose ;
  double yaw , pitch , roll ;
  projection_t projection ;
  std::string prj_str ;
  double x0 , x1 , y0 , y1 ;
  double hfov ;
  std::size_t width , height ;
  std::string input ;
  std::string output ;
  std::string metamatch ;
  std::regex field_re ;
  std::size_t env_width , env_height ;
  std::size_t nchannels ;

  double get_vfov()
  {
    double vfov = 0.0 ;
    switch ( projection )
    {
      case RECTILINEAR:
      {
        // as a one-liner, this is probably clearer than the code below
         vfov = 2.0 * atan ( height * tan ( hfov / 2.0 ) / width ) ;
         break ;
      }
      case CYLINDRICAL:
      {
        double pixels_per_rad = width / hfov ;
        double h_rad = height / pixels_per_rad ;
        vfov = 2.0 * atan ( h_rad / 2.0 ) ;
        break ;
      }
      case STEREOGRAPHIC:
      {
        double w_rad = 2.0 * tan ( hfov / 4.0 ) ;
        double pixels_per_rad = width / w_rad ;
        double h_rad = height / pixels_per_rad ;
        vfov = 4.0 * atan ( h_rad / 2.0 ) ;
        break ;
      }
      case SPHERICAL:
      case FISHEYE:
      {
        vfov = hfov * height / width ;
        break ;
      }
      default:
      {
        vfov = hfov ; // debatable...
        break ;
      }
    }
    std::cout << "gleaned vfov = " << vfov << std::endl ;
    return vfov ;
  }

  void get_extent()
  {
    double alpha_x = - hfov / 2.0 ;
    double beta_x = hfov / 2.0 ;
    double beta_y = get_vfov() / 2.0 ;
    double alpha_y = - beta_y ;

    switch ( projection )
    {
      case SPHERICAL:
      case FISHEYE:
      {
        x0 = alpha_x ;
        x1 = beta_x ;

        y0 = alpha_y ;
        y1 = beta_y ;
        break ;
      }
      case CYLINDRICAL:
      {
        x0 = alpha_x ;
        x1 = beta_x ;

        y0 = tan ( alpha_y ) ;
        y1 = tan ( beta_y ) ;
        break ;
      }
      case RECTILINEAR:
      {
        x0 = tan ( alpha_x ) ;
        x1 = tan ( beta_x ) ;

        y0 = tan ( alpha_y ) ;
        y1 = tan ( beta_y ) ;
        break ;
      }
      case STEREOGRAPHIC:
      {
        x0 = 2.0 * tan ( alpha_x / 2.0 ) ;
        x1 = 2.0 * tan ( beta_x / 2.0 ) ;

        y0 = 2.0 * tan ( alpha_y / 2.0 ) ;
        y1 = 2.0 * tan ( beta_y / 2.0 ) ;
        break ;
      }
      case CUBEMAP:
      {
        x0 = tan ( alpha_x ) ;
        x1 = tan ( beta_x ) ;

        y0 = 6 * x0 ;
        y1 = 6 * x1 ;
        break ;
      }
      default:
      {
        break ;
      }
    }
  }

  arguments ( int argc , const char ** argv )
  {
    // we're using OIIO's argparse, since we're using OIIO anyway.
    // This is a convenient way to glean arguments on all supported
    // platforms - getopt isn't available everywhere.

    convert_native_arguments(argc, (const char**)argv);
    ArgParse ap;
    ap.intro("stepper: extract image from an environment\n")
      .usage("extract [options] --input INPUT --output OUTPUT");
    ap.arg("-v", &verbose)
      .help("Verbose output");
    ap.arg("--input INPUT")
      .help("input file name (mandatory)")
      .metavar("INPUT");
    ap.arg("--output OUTPUT")
      .help("output file name (mandatory)")
      .metavar("OUTPUT");
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
    
    if (ap.parse(argc, argv) < 0 ) {
        std::cerr << ap.geterror() << std::endl;
        ap.print_help();
        assert ( false ) ;
    }

    if (!metamatch.empty()) {
        field_re.assign(metamatch, std::regex_constants::extended
                                   | std::regex_constants::icase);
    }
  
    input = ap["input"].as_string ( "" ) ;
    output = ap["output"].as_string ( "" ) ;
    x0 = ap["x0"].get<float> ( 0.0 ) ;
    x1 = ap["x1"].get<float> ( 0.0 ) ;
    y0 = ap["y0"].get<float> ( 0.0 ) ;
    y1 = ap["y1"].get<float> ( 0.0 ) ;
    width = ap["width"].get<int> ( 0 ) ;
    height = ap["height"].get<int> ( 0 ) ;
    hfov = ap["hfov"].get<float>(0.0);
    if ( hfov != 0.0 )
    {
      hfov *= M_PI / 180.0 ;
      x0 = x1 = y0 = y1 = 0 ;
    }
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
    assert ( output != std::string() ) ;
    assert ( projection != PRJ_NONE ) ;
    assert ( width > 0 ) ;
    assert ( height > 0 ) ;
    if ( projection == CUBEMAP )
      assert ( height = 6 * width ) ;
    if ( hfov != 0.0 )
    {
      get_extent() ;
    }
    assert ( x0 < x1 ) ;
    assert ( y0 < y1 ) ;

    // some member variables in the args object are gleaned from
    // the input image:

    auto inp = ImageInput::open ( input ) ;
    assert ( inp ) ;

    const ImageSpec &spec = inp->spec() ;

    env_width = spec.width ;
    env_height = spec.height ;
    nchannels = spec.nchannels ;

    assert (    env_width == env_height * 2
             || env_height == env_width * 6 ) ;

    inp->close() ;

    if ( verbose )
    {
      std::cout << "input: " << input << std::endl ;
      std::cout << "input width: " << env_width << std::endl ;
      std::cout << "input height: " << env_height << std::endl ;
      std::cout << "input has " << nchannels << " channels" << std::endl ;
    
      std::cout << "output: " << output << std::endl ;
      std::cout << "output projection: " << prj_str << std::endl ;
      std::cout << "width: " << width
                << " height: " << height << std::endl ;

      if ( hfov > 0.0 )
      {
        std::cout << "hfov: " << hfov << std::endl ;
        std::cout << "extent gleaned from hfov:" << std::endl ;
      }
      else
      {
        std::cout << "extent gleaned from command line arguments:"
                  << std::endl ;
      }
      std::cout << "x0: " << x0 << " x1: " << x1 << std::endl ;
      std::cout << "y0: " << y0 << " y1: " << y1 << std::endl ;
    }
  }
} ;

// environment.h has most of the 'workhorse' code for this demo.

#include "environment.h"

template < std::size_t nchannels >
void work ( const arguments & args ,
            zimt::grok_get_t < float , 9 , 2 , 16 > & get_ray )
{
  // set up the environment object yielding content. This serves as
  // the 'act' functor for zimt::process

  environment < float , float , nchannels , 16 > env ( args.input ) ;

  typedef zimt::xel_t < float , nchannels > px_t ;
  
  zimt::array_t < 2 , px_t > trg ( { args.width , args.height } ) ;
  
  // set up zimt::storers to populate the target arrays with
  // zimt::process
  
  zimt::storer < float , nchannels , 2 , 16 > cstor ( trg ) ;
  
  // use the get, act and put components with zimt::process
  // to produce the target images and store them to disk
  
  if ( args.verbose )
    std::cout << "producing output" << std::endl ;
  
  zimt::process ( trg.shape , get_ray , env , cstor ) ;
  
  if ( args.verbose )
    std::cout << "saving output image: " << args.output << std::endl ;
  
  save_array < nchannels > ( args.output , trg ) ;
  
  if ( args.verbose )
    std::cout << "done." << std::endl ;
}

int main ( int argc , const char ** argv )
{
  // process command line arguments - the result is held in a bunch
  // of member variables in the 'args' object, to be passed around
  // conveniently.

  arguments args ( argc , argv ) ;

  // orthonormal system of basis vectors the view

  crd3_t xx { 1.0 , 0.0 , 0.0 } ;
  crd3_t yy { 0.0 , 1.0 , 0.0 } ;
  crd3_t zz { 0.0 , 0.0 , 1.0 } ;

  // the three vectors are rotated with the given yaw, pitch
  // and roll, and later passed on the to 'steppers', the objects
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

  // set up the steppers. note the extents of the 2D manifold
  // given in model space uits. These objects will deliver 3D
  // 'ray' coordinates as input to the 'act' functor. Not how
  // each projection has it's distinct type of stepper, but
  // they are all assigned to a common type, a 'grok_get_t'.
  // This uses type erasure and captures the functionality
  // in std::functions, and the resulting object is only
  // characterized by it's input and output type and lane
  // count.

  // zimt::grok_get_t < float , 3 , 2 , 16 > get_ray ;
  zimt::grok_get_t < float , 9 , 2 , 16 > get_ray ;

  switch ( args.projection )
  {
    case RECTILINEAR :
      get_ray = deriv_stepper < float , 16 , rectilinear_stepper >
                  ( xx , yy , zz , args.width , args.height ,
                    args.x0 , args.x1 , args.y0 , args.y1 ) ;
      // get_ray = rectilinear_stepper < float , 16 >
      //             ( xx , yy , zz , args.width , args.height ,
      //               args.x0 , args.x1 , args.y0 , args.y1 ) ;
      break ;
    case FISHEYE :
      get_ray = deriv_stepper < float , 16 , fisheye_stepper >
                  ( xx , yy , zz , args.width , args.height ,
                    args.x0 , args.x1 , args.y0 , args.y1 ) ;
      // get_ray = fisheye_stepper < float , 16 >
      //             ( xx , yy , zz , args.width , args.height ,
      //               args.x0 , args.x1 , args.y0 , args.y1 ) ;
      break ;
    case STEREOGRAPHIC :
      get_ray = deriv_stepper < float , 16 , stereographic_stepper >
                  ( xx , yy , zz , args.width , args.height ,
                    args.x0 , args.x1 , args.y0 , args.y1 ) ;
      // get_ray = stereographic_stepper < float , 16 >
      //             ( xx , yy , zz , args.width , args.height ,
      //               args.x0 , args.x1 , args.y0 , args.y1 ) ;
      break ;
    case SPHERICAL :
      get_ray = deriv_stepper < float , 16 , spherical_stepper >
                  ( xx , yy , zz , args.width , args.height ,
                    args.x0 , args.x1 , args.y0 , args.y1 ) ;
      // get_ray = spherical_stepper < float , 16 >
      //             ( xx , yy , zz , args.width , args.height ,
      //               args.x0 , args.x1 , args.y0 , args.y1 ) ;
      break ;
    case CYLINDRICAL :
      get_ray = deriv_stepper < float , 16 , cylindrical_stepper >
                  ( xx , yy , zz , args.width , args.height ,
                    args.x0 , args.x1 , args.y0 , args.y1 ) ;
      // get_ray = cylindrical_stepper < float , 16 >
      //             ( xx , yy , zz , args.width , args.height ,
      //               args.x0 , args.x1 , args.y0 , args.y1 ) ;
      break ;
    case CUBEMAP :
      get_ray = deriv_stepper < float , 16 , cubemap_stepper >
                  ( xx , yy , zz , args.width , args.height ,
                    args.x0 , args.x1 , args.y0 , args.y1 ) ;
      // get_ray = cubemap_stepper < float , 16 >
      //             ( xx , yy , zz , args.width , args.height ,
      //               args.x0 , args.x1 , args.y0 , args.y1 ) ;
    default:
      break ;
  }

  // the code to set up and extract data from the environment
  // depends on the number of channels. That's passed as a template
  // argument, so here we have a case switch over 'nchannels' which
  // dispatches to the appropriate instantiations.

  switch ( args.nchannels )
  {
    case 1:
      work < 1 > ( args , get_ray ) ;
      break ;
    case 2:
      work < 2 > ( args , get_ray ) ;
      break ;
    case 3:
      work < 3 > ( args , get_ray ) ;
      break ;
    case 4:
      work < 4 > ( args , get_ray ) ;
      break ;
  }
}
    
