/************************************************************************/
/*                                                                      */
/*    envutil - utility to convert between environment formats          */
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

// demo program to demonstrate the use of 'steppers', as defined in
// "stepper.h". This program takes a 2:1 lat/lon environment image as
// input (preferably in openEXR format) and produces several output
// images, named like their projection, which show part or all of
// the environment image. The environment can be given either as
// a 2:1 lat/lon image or a 1:6 cubemap image.
// Note how most supporting functionality has been factored out and
// placed in headers. The code to use these facilities is very terse.

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

// environment.h has most of the 'workhorse' code for this demo.

#include "environment.h"

int main ( int argc , char * argv[] )
{
  // obtain environment image data

  if ( argc <= 1 )
  {
    std::cerr << "pass a 2:1 lat/lon or 1:6 cubemap image" << std::endl ;
  }

  auto inp = ImageInput::open ( argv[1] ) ;
  assert (inp ) ;

  std::size_t w ;
  std::size_t h ;
  std::size_t nchannels ;

  const ImageSpec &spec = inp->spec() ;

  w = spec.width ;
  h = spec.height ;
  nchannels = spec.nchannels ;

  assert ( w == h * 2 || h == w * 6 ) ;

  inp->close() ;

  // for simplicity's sake, only produce RGB output

  typedef zimt::xel_t < float , 3 > px_t ;

  // we use three different output shapes: square for rectilinear,
  // fisheye and stereographic targets, 2:1 for spherical and
  // cylindrical ones, and finally 1:6 for the cubemap.

  zimt::array_t < 2 , px_t > trg ( { 3000 , 3000 } ) ;
  zimt::array_t < 2 , px_t > trg2 ( { 6000 , 3000 } ) ;
  zimt::array_t < 2 , px_t > trg3 ( { 1000 , 6000 } ) ;

  // set up zimt::storers to populate the target arrays with
  // zimt::process

  zimt::storer < float , 3 , 2 , 16 > cstor ( trg ) ;
  zimt::storer < float , 3 , 2 , 16 > cstor2 ( trg2 ) ;
  zimt::storer < float , 3 , 2 , 16 > cstor3 ( trg3 ) ;

  // orthonormal system for the view

  crd3_t xx { 1.0 , 0.0 , 0.0 } ;
  crd3_t yy { 0.0 , 1.0 , 0.0 } ;
  crd3_t zz { 0.0 , 0.0 , 1.0 } ;

  double m_pi_4 = M_PI / 4.0 ;

  // play with rotation: uncomment this bit to use a 45 degree
  // camera roll:

  // rotate_3d < float , 16 > r3 ( m_pi_4 , 0.0 , 0 ) ;
  // 
  // xx = r3 ( xx ) ;
  // yy = r3 ( yy ) ;
  // zz = r3 ( zz ) ;

  // set up the steppers. note the extents of the 2D manifold
  // given in model space uits.

  rectilinear_stepper < float , 16 >
    rst ( xx , yy , zz , 3000 , 3000 ,
          -m_pi_4 , m_pi_4 , -m_pi_4 , m_pi_4 ) ;

  fisheye_stepper < float , 16 >
    fst ( xx , yy , zz , 3000 , 3000 ,
          -M_PI , M_PI , -M_PI , M_PI ) ;
  
  // what's this unobvious extent in model space? This is owed to
  // the stereographic projection and corresponds to ninety degrees
  // field of view, which comes out slighly dilated from the
  // +/- pi/4 which we'd use for most other projections and this fov.

  stereographic_stepper < float , 16 >
    stst ( xx , yy , zz , 3000 , 3000 ,
            -.828427 , .828427 , -.828427 , .828427 ) ;
  
  spherical_stepper < float , 16 >
    sst ( xx , yy , zz , 6000 , 3000 ,
          -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;
  
  cylindrical_stepper < float , 16 >
    cst ( xx , yy , zz , 6000 , 3000 ,
          -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  cubemap_stepper < float , 16 >
    cbmst ( xx , yy , zz , 1000 , 6000 ,
          -1.0 , 1.0 , -6.0 , 6.0 ) ;

  // set up the environment object yielding content. This serves as
  // the 'act' functor for zimt::process

  environment < float , float , 3 , 16 > env ( argv[1] ) ;

  // use the get, act and put components with zimt::process
  // to produce the target images and store them to disk

  zimt::process ( trg.shape , rst , env , cstor ) ;
  save_array < 3 > ( "rectilinear.exr" , trg ) ;

  zimt::process ( trg.shape , fst , env , cstor ) ;
  save_array < 3 > ( "fisheye.exr" , trg ) ;
  
  zimt::process ( trg.shape , stst , env , cstor ) ;
  save_array < 3 > ( "stereographic.exr" , trg ) ;
  
  zimt::process ( trg2.shape , sst , env , cstor2 ) ;
  save_array < 3 > ( "spherical.exr" , trg2 ) ;
  
  zimt::process ( trg2.shape , cst , env , cstor2 ) ;
  save_array < 3 > ( "cylindrical.exr" , trg2 ) ;
  
  zimt::process ( trg3.shape , cbmst , env , cstor3 ) ;
  save_array < 3 > ( "cubemap.exr" , trg3 ) ;
}
    
