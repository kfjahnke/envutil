/************************************************************************/
/*                                                                      */
/*    zimt - abstraction layer for SIMD programming                     */
/*                                                                      */
/*            Copyright 2024 by Kay F. Jahnke                           */
/*                                                                      */
/*    The git repository for this software is at                        */
/*                                                                      */
/*    https://github.com/kfjahnke/zimt                                  */
/*                                                                      */
/*    Please direct questions, bug reports, and contributions to        */
/*                                                                      */
/*    kfjahnke+zimt@gmail.com                                           */
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

// This is a test program for the geometry code used in envutil.

// There are two different ways of compiling it. The first is to create 
// a single-ISA binary, fixing the SIMD ISA at compile time by passing
// appropriate flags to the compiler. This will work with all SIMD
// back-ends - here, I show the compiler invocations for an AVX2 version.
// Note that - especially when using highway - additional flags may be
// needed to get a program which actually uses the AVX2 instructions.
// When highway uses pragmas to set the stage for AVX2 code, it uses
// all these flags: -msse2 -mssse3 -msse4.1 -msse4.2 -mpclmul -maes
// -mavx -mavx2 -mbmi -mbmi2 -mfma -mf16c
// Note also the -I. directive to tell the compiler to find files to
// #include in the current folder as well.

//   clang++ -mavx2 geometry.cc -O3 -I. -DUSE_HWY -lhwy
//   clang++ -mavx2 geometry.cc -O3 -I. -DUSE_VC -lVc
//   clang++ -mavx2 geometry.cc -O3 -I. -DUSE_STDSIMD
//   clang++ -mavx2 geometry.cc -O3 -I.

// The second way is to use highway's automatic dispatch to embedded
// variants of the code running with different ISAs. This requires the
// definition of MULTI_SIMD_ISA and linkage to libhwy. Here, no
// architecture flags are passed to the compiler:

//   clang++ geometry.cc -O3 -I. -DMULTI_SIMD_ISA -DUSE_HWY -lhwy
//   clang++ geometry.cc -O3 -I. -DMULTI_SIMD_ISA -DUSE_STDSIMD -lhwy
//   clang++ geometry.cc -O3 -I. -DMULTI_SIMD_ISA -DUSE_VC -lhwy -lVc
//   clang++ geometry.cc -O3 -I. -DMULTI_SIMD_ISA -lhwy

// binaries made with the second method will dispatch to what is deemd
// the best SIMD ISA available on the CPU on which the binary is run.
// Because this is done meticulously by highway's CPU detection code,
// the binary variant which is picked is usually optimal and may
// out-perform single-ISA variants with 'manually' supplied ISA flags,
// if the set of flags isn't optimal as well. The disadvantage of the
// multi-SIMD-ISA variants is (much) longer compile time and code size.

// Due to the 'commodification' the source code itself doesn't have
// to be modified in any way to produce one variant or another.
// This suggests that during the implementation of a new program a
// fixed-ISA build can be used to evolve the code with fast turn-around
// times, adding dispatch capability later on by passing the relevant
// compiler flags.

// we define a dispatch base class. All the 'payload' code is called
// through virtual member functions of this class. In this example,
// we only have a single payload function. We have to enclose this
// base class definition in an include guard, because it must not
// be compiled repeatedly, which happens when highway's foreach_target
// mechansim is used. This definition might go to a header file, but
// it's better placed here, because it holds the declarations of the
// pure virtual member function(s) used as conduit to the ISA-specific
// code, and here, in an example, we want to see both these declarations
// and, later on, the implementations, in the same file.

#ifndef DISPATCH_BASE
#define DISPATCH_BASE

struct dispatch_base
{
  // in dispatch_base and derived classes, we keep two flags.
  // 'backend' holds a value indicating which of zimt's back-end
  // libraries is used. 'hwy_isa' is only set when the highway
  // backend is used and holds highway's HWY_TARGET value for
  // the given nested namespace.

  int backend = -1 ;
  unsigned long hwy_isa = 0 ;

  // next we have pure virtual member function definitions for
  // payload code. In this example, we only have one payload
  // function which calls what would be 'main' in a simple
  // program without multiple SIMD ISAs or SIMD back-ends

  virtual int payload ( int argc , char * argv[] ) const = 0 ;
} ;

#endif

#ifdef MULTI_SIMD_ISA

// if we're using MULTI_SIMD_ISA, we have to define HWY_TARGET_INCLUDE
// to tell the foreach_target mechanism which file should be repeatedly
// re-included and re-copmpiled with SIMD-ISA-specific flags

#undef HWY_TARGET_INCLUDE

/////////////// Tell highway which file to submit to foreach_target

#define HWY_TARGET_INCLUDE "geometry.cc"  // this very file

//--------------------------------------------------------------------

#include <hwy/foreach_target.h>  // must come before highway.h

#include <hwy/highway.h>

#endif // #ifdef MULTI_SIMD_ISA

//////////////// Put the #includes needed for your program here:

#include <random>
#include <chrono>
#include "geometry.h"
#include "stepper.h"

// to make highway's use of #pragma directives to the compiler
// effective, we surround the SIMD-ISA-specific code with
// HWY_BEFORE_NAMESPACE() and HWY_AFTER_NAMESPACE().

HWY_BEFORE_NAMESPACE() ;

// To conveniently rotate with a rotational quaternion, we employ
// Imath's 'Quat' data type, packaged in a zimt::unary_functor.
// This is not factored out because it requires inclusion of
// some Imath headers, which I want to keep out of the other
// code, e.g. in geometry.h, where it would fit in nicely.

#include <Imath/ImathVec.h>
#include <Imath/ImathEuler.h>
#include <Imath/ImathQuat.h>
#include <Imath/ImathLine.h>

//--------------------------------------------------------------------

// this macro puts us into a nested namespace inside namespace 'project'.
// For single-SIMD-ISA builds, this is conventionally project::zsimd,
// and for multi-SIMD-ISA builds it is project::HWY_NAMESPACE. The macro
// is defined in common.h. After the macro invocation, we can access
// all zimt names with a simple zimt:: prefix - both 'true' zimt names
// and SIMD-ISA-specific versions living in the nested namespace.

BEGIN_ZIMT_SIMD_NAMESPACE(project)

// rotate_3d uses a SIMDized Imath Quaternion to affect a 3D rotation
// of a 3D SIMDized coordinate. Imath::Quat<float> can't broadcast
// to handle SIMDized input, but if we use an Imath::Quat of the
// SIMDized type, we get the desired effect for simdized input -
// hence the broadcasting to Imath::Quat<U> in 'eval', which
// has no effect for scalars, but represents a broadcast of the
// rotation to all lanes of the simdized quaternion's components.
// We'll use this functor to compare the output of steppers with
// built-in rotation to unrotated steppers with a subsequent
// rotation of the resulting 3D ray.

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

  // for the actual rotation (the 'multiplication' with teh rotational
  // quaternion), we don't use Imath code. See comments in 'eval'.

  // calculate the cross product of two vectors

  template < typename U >
  U cross ( const U & x , const U & y ) const
  {
    U result ;
    result[0] = x[1] * y[2] - x[2] * y[1] ;
    result[1] = x[2] * y[0] - x[0] * y[2] ;
    result[2] = x[0] * y[1] - x[1] * y[0] ;
    return result ;
  }

  // perform the quaternion multiplication.

  template < typename U >
  U mulq ( const U & vec ) const
  {
    U result ;
    U qv { q.v[0] , q.v[1] , q.v[2] } ;
    auto a = cross ( qv , vec ) ;
    auto b = cross ( qv , a ) ;
    return vec + T(2) * (q.r * a + b);
  }

  // eval applies the quaternion. Note how we use a template of
  // typename U for the formulation. This way, we can handle both
  // scalar and simdized arguments.

  template < typename U >
  void eval ( const zimt::xel_t < U , 3 > & in ,
              zimt::xel_t < U , 3 > & out ) const
  {
    // using highway's foreach_target mechanism, coding the quaternion
    // multiplication by invoking Imath's operator* with an Imath::Quat
    // argument, produces slow code compared to single-ISA compiles.
    // I work around this problem by coding the operation 'manually'.
    // Code using Imath's operator* would do this:
    //
    auto const & in_e
      = reinterpret_cast < const Imath::Vec3 < U > & > ( in ) ;
    
    auto & out_e
      = reinterpret_cast < Imath::Vec3 < U > & > ( out ) ;
    
    out_e = in_e * Imath::Quat < U > ( q ) ;

    // instead, we calculate 'manually' like this, using the mulq
    // member function above.

    // out = mulq ( in ) ;
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

typedef zimt::xel_t < double , 2 > d2_t ;
typedef zimt::xel_t < double , 3 > d3_t ;
typedef zimt::simdized_type < d3_t , LANES > d3_v ;

template < typename pa , typename pb >
d3_t work ( d3_t _c3 )
{
  // set up a ray-to-ray functor containing the two functors
  // of type 'pa' and 'pb'. The ray-to-ray functor takes a
  // 3D ray coordinate, converts it to a 2D coordinate with
  // the projection codified in pa and then back to a ray
  // with the projection codified in pb. Here, we test pa
  // and pb using the same projection - in actual code, pa
  // and pb might use different projections.

  pa tf1 ;
  pb tf2 ;
  ray_to_ray < double , LANES > tf ( tf1 , tf2 ) ;

  // we want to make sure that the pa and pb functors do
  // their job with both scalar and simdized values, so we
  // run the test also with a simdized evaluation.

  d3_v cv ( _c3 ) ;
  d3_t c3 ;
  tf.eval ( _c3 , c3 ) ;
  tf.eval ( cv , cv ) ;

  // we test that the first lane of the vectorized output
  // is very close to the scalar output. With some back-ends,
  // the values are actually equal, but we can't rely on it.

  assert ( std::abs ( cv[0][0] - c3[0] ) < .0000000000001 ) ;
  assert ( std::abs ( cv[1][0] - c3[1] ) < .0000000000001 ) ;
  assert ( std::abs ( cv[2][0] - c3[2] ) < .0000000000001 ) ;

  // we compare input and output and make sure that they
  // differ if at all then only minimally. Note how we compare
  // the values after normalization - they may differ by a
  // factor, because the code does not necessarily produce
  // normalized output - nor is the input normalized.
  
  assert ( std::abs ( _c3[0] / norm ( _c3 ) - c3[0] / norm ( c3 ) )
           < .0000000000001 ) ;

  assert ( std::abs ( _c3[1] / norm ( _c3 ) - c3[1] / norm ( c3 ) )
           < .0000000000001 ) ;

  assert ( std::abs ( _c3[2] / norm ( _c3 ) - c3[2] / norm ( c3 ) )
           < .0000000000001 ) ;

  return c3 ;
}

// we use a two-level dispatch to get from run-time values for
// the projections to types of conversion functors

template < typename pa >
d3_t route2 ( projection_t pb ,
                 d3_t c3 )
{
  d3_t result ;

  switch ( pb )
  {
    case SPHERICAL:
      result = work < pa , ll_to_ray_t < double > > ( c3 ) ;
      break ;
    case CYLINDRICAL:
      result = work < pa , cyl_to_ray_t < double > > ( c3 ) ;
      break ;
    case RECTILINEAR:
      result = work < pa , rect_to_ray_t < double > > ( c3 ) ;
      break ;
    case FISHEYE:
      result = work < pa , fish_to_ray_t < double > > ( c3 ) ;
      break ;
    case STEREOGRAPHIC:
      result = work < pa , ster_to_ray_t < double > > ( c3 ) ;
      break ;
    case CUBEMAP:
      result = work < pa , ir_to_ray_t < double > > ( c3 ) ;
      break ;
    case BIATAN6:
      result = work < pa , ba6_to_ray_t < double > > ( c3 ) ;
      break ;
    default:
      break ;
  }

  return result ;
}

d3_t route ( projection_t pa ,
                projection_t pb ,
                d3_t c3 )
{
  d3_t result ;

  switch ( pa )
  {
    case SPHERICAL:
      result = route2 < ray_to_ll_t < double > > ( pb , c3 ) ;
      break ;
    case CYLINDRICAL:
      result = route2 < ray_to_cyl_t < double > > ( pb , c3 ) ;
      break ;
    case RECTILINEAR:
      // negative z axis doesn't work with rectilinear projection, hence:
      c3[2] = std::abs ( c3[2] ) ;

      result = route2 < ray_to_rect_t < double > > ( pb , c3 ) ;
      break ;
    case FISHEYE:
      result = route2 < ray_to_fish_t < double > > ( pb , c3 ) ;
      break ;
    case STEREOGRAPHIC:
      result = route2 < ray_to_ster_t < double > > ( pb , c3 ) ;
      break ;
    case CUBEMAP:
      result = route2 < ray_to_ir_t < double > > ( pb , c3 ) ;
      break ;
    case BIATAN6:
      result = route2 < ray_to_ba6_t < double > > ( pb , c3 ) ;
      break ;
    default:
      break ;
  }

  return result ;
}

void test_r2r ( d3_t d3 )
{
  d3_t result ;

  // for this test, we set up a ray_to_... functor as first
  // functor, then a ..._to_ray functor as second functor.
  // The projection is the same for the two functors, but
  // the direction of the transformation (3D->2D vs. 2D->3D)
  // is opposite.

  for ( int pa = SPHERICAL ; pa <= BIATAN6 ; pa++ )
  {
    int pb = pa ;
    {
      result = route ( projection_t(pa) ,
                       projection_t(pb) ,
                       d3 ) ;
    }
  }
}

int _payload ( int argc , char * argv[] )
{
  // first, we run a test where we project rays to a planar surface
  // and back to the ray. This should always succeed (except for
  // rectilinear projection and negative z axis, which is avoided)
  // because the resulting 2D coordinate can always regenerate the
  // ray precisely, whereas the opposite operation (start with a
  // planar coordinate, move to ray and back) will fail for certain
  // planar coordinates - e.g. with spherical projection and y==0
  // where the initial x coordinate can't be recovered.

  std::mt19937 gen ; // Standard mersenne_twister_engine
  std::uniform_real_distribution<> dis(-10.0, 10.0);

  for ( std::size_t i = 0 ; i < 10000 ; i++ )
  {
    test_r2r ( { dis(gen) , dis(gen) , dis(gen) } ) ;
  }

  // Next, we populate arrays of ray coordinates with 'steppers'
  // and test the result against arrays which are populated using
  // the ..._to_ray_t functors. We also test against values gleaned
  // by directly invoking the functors' eval function with input
  // generated by a lambda directly from discrete coordinates.

  // two equally-shaped target arrays take up the output of
  // the stepper and the output from the plain conversion
  // operator, which should be identical.

  zimt::array_t < 2 , d3_t > a3 ( { 1000 , 500 } ) ;
  zimt::array_t < 2 , d3_t > b3 ( { 1000 , 500 } ) ;

  // first set up the stepper.

  spherical_stepper < double , LANES >
    sphs ( { 1.0 , 0.0 , 0.0 } ,
           { 0.0 , 1.0 , 0.0 } ,
           { 0.0 , 0.0 , 1.0 } ,
           1000 ,
           500 ) ;

  // use a zimt::process run to fill the first target array with
  // output from the stepper. We use a pass_through act functor
  // because we're only interested in the stepper's output.

  zimt::process ( a3.shape , sphs ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  // now we set up parameters for the use of the conversion
  // operator. The 2D coordinates - in model space units - run
  // from (x0, y0) with a step of dx in the horizontal and dy
  // in the vertical. Note how we're working with centered
  // coordinates, as we would in a viewing context, where
  // the images are draped in model space so that their center
  // rides on the forward axis.

  double x0 = ( a3.shape[0] - 1 ) / -2.0 ;
  double y0 = ( a3.shape[1] - 1 ) / -2.0 ;
  double dx = 2.0 * M_PI / a3.shape[0] ;
  double dy = M_PI / a3.shape[1] ;
  x0 *= dx ;
  y0 *= dy ;

  // this is the conversion operator:

  ll_to_ray_t < double , LANES > ll_to_ray ;

  // we set up a zimt linspace generator as input

  zimt::linspace_t < double , 2 , 2 , LANES >
    ls ( { x0 , y0 } , { dx , dy } ) ;

  // and use zimt::process to populate the second target array

  zimt::process ( b3.shape , ls , ll_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  // for the 'manual' test we set up a lambda which provides
  // centered 2D model space coordinates for discrete coordinates

  auto to_md = [&] ( std::size_t sx , std::size_t sy )
  {
    zimt::xel_t < double , 2 > crd2 ;
    double & x ( crd2[0] ) ;
    double & y ( crd2[1] ) ;
    x = sx - ( a3.shape[0] - 1 ) / 2.0 ;
    y = sy - ( a3.shape[1] - 1 ) / 2.0 ;
    x *= 2.0 * M_PI / a3.shape[0] ;
    y *= M_PI / a3.shape[1] ;
    return crd2 ;
  } ;

  // now we feed the converter in a loop. We use the discrete
  // coordinates to test the result values in arrays a3 and b3
  // against the results we receive as we go.

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto crd2 = to_md ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      ll_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  // Now we set up a rotation functor, and a stepper with the
  // same built-in rotation. We set up the rotation with three
  // 'odd' values (no multiples of pi).

  rotate_3d < double , LANES > r3 ( 1.0 , 2.0 , 3.0 ) ;

  // Now we produce the three basis vectors for the rotated stepper

  d3_t ex { 1.0 , 0.0 , 0.0 } ;
  d3_t ey { 0.0 , 1.0 , 0.0 } ;
  d3_t ez { 0.0 , 0.0 , 1.0 } ;
  r3.eval ( ex , ex ) ;
  r3.eval ( ey , ey ) ;
  r3.eval ( ez , ez ) ;

  spherical_stepper < double , LANES >
    rsphs ( ex , ey , ez , 1000 , 500 ) ;

  // use a zimt::process run to fill the first target array with
  // output from the rotated stepper. We time this and the subsequent
  // run with the alternative computation scheme to see how performance
  // differs - on my system, using the rotated stepper is 30-50 % faster.
  // Of course, this may be due to the fact that the application of the
  // rotational quaternion with Imath uses a quaternion of simdized
  // components, even though all lanes of these components are equal.
  // So to really compare performance, we'd have to do the rotation
  // 'manually', applying the scalar components, but I assume that
  // the output would be roughly the same.
  // A note about the results of the speed measurements: with the code,
  // as it stands, compiling a single-SIMD-ISA binary with explicitly
  // stated ISA flags results in significantly faster code for the
  // second run - the one with the separate rotation after the coordinate
  // transformation with ll_to_ray - compared to the multi-SIMD-ISA
  // version of the program. Why is that? It's due to the use of Imath's
  // quaternion code for the rotation. Imath is not adapted to use
  // highway's foreach_target mechanism, so including the Imath headers
  // for the first time instantiates the Imath types finally, and this
  // happens with a low-grade ISA. Subsequent re-compilations with the
  // foreach_target mechanism can't reinstantiate the templates, because
  // the sentinels in the Imath headers blank out the code - Imath
  // 'thinks' they have already been dealt with. This results in
  // quaternion code which is hobbled to use only the instantiations
  // from the first compilation, and this significantly degrades
  // performance. For now, I see no way out of this dilemma, short of
  // falling back to a scheme with several separate TUs for the SIMD
  // ISAs. Note, though, that this problem affects only the code which
  // we use to test that the steppers perform as expected! The steppers
  // don't use Imath's quaternion code at all. This is why the first
  // run below comes out equally fast in both modes of compilation,
  // and it's another good reason to use steppers if possible.

  std::chrono::system_clock::time_point start
    = std::chrono::system_clock::now() ;
  
  for ( int times = 0 ; times < 100 ; times++ )
    zimt::process ( a3.shape , rsphs ,
                    zimt::pass_through < double , 3 , LANES > () ,
                    zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  std::chrono::system_clock::time_point end
    = std::chrono::system_clock::now() ;
  
  std::cout << "first run took:               "
            << std::chrono::duration_cast<std::chrono::milliseconds>
                ( end - start ) . count()
            << " ms" << std::endl ;

  // Now we fill b3 with the concatenation of the unrotated
  // stepper and the rotation functor

  start = std::chrono::system_clock::now() ;

  for ( int times = 0 ; times < 100 ; times++ )
    zimt::process ( b3.shape , ls , ll_to_ray + r3 ,
                    zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  end = std::chrono::system_clock::now() ;
  
  std::cout << "second run took:               "
            << std::chrono::duration_cast<std::chrono::milliseconds>
                ( end - start ) . count()
            << " ms" << std::endl ;

  // and look at the results

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto d = abs ( a3 [ { x , y } ] - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  // we repeat the process for all projections.

  cylindrical_stepper < double , LANES , false >
    cyls ( { 1.0 , 0.0 , 0.0 } ,
           { 0.0 , 1.0 , 0.0 } ,
           { 0.0 , 0.0 , 1.0 } ,
           1000 , 500 ,
           -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  zimt::process ( a3.shape , cyls ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  cyl_to_ray_t < double , LANES > cyl_to_ray ;

  zimt::process ( b3.shape , ls ,
                  cyl_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto crd2 = to_md ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      cyl_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  cylindrical_stepper < double , LANES , false >
    rcyls ( ex , ey , ez , 1000 , 500 ,
            -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  zimt::process ( a3.shape , rcyls ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  zimt::process ( b3.shape , ls , cyl_to_ray + r3 ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto d = abs ( a3 [ { x , y } ] - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  rectilinear_stepper < double , LANES , false >
    rects ( { 1.0 , 0.0 , 0.0 } ,
            { 0.0 , 1.0 , 0.0 } ,
            { 0.0 , 0.0 , 1.0 } ,
            1000 , 500 ,
            -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  zimt::process ( a3.shape , rects ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  rect_to_ray_t < double , LANES > rect_to_ray ;

  zimt::process ( b3.shape , ls ,
                  rect_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto crd2 = to_md ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      rect_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  rectilinear_stepper < double , LANES , false >
    rrects ( ex , ey , ez , 1000 , 500 ,
            -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  zimt::process ( a3.shape , rrects ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  zimt::process ( b3.shape , ls , rect_to_ray + r3 ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto d = abs ( a3 [ { x , y } ] - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  fisheye_stepper < double , LANES , false >
    fishs ( { 1.0 , 0.0 , 0.0 } ,
            { 0.0 , 1.0 , 0.0 } ,
            { 0.0 , 0.0 , 1.0 } ,
            1000 , 500 ,
            -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  zimt::process ( a3.shape , fishs ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  fish_to_ray_t < double , LANES > fish_to_ray ;

  zimt::process ( b3.shape , ls ,
                  fish_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto crd2 = to_md ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      fish_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  fisheye_stepper < double , LANES , false >
    rfishs ( ex , ey , ez , 1000 , 500 ,
            -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  zimt::process ( a3.shape , rfishs ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  zimt::process ( b3.shape , ls , fish_to_ray + r3 ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto d = abs ( a3 [ { x , y } ] - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  stereographic_stepper < double , LANES , false >
    sters ( { 1.0 , 0.0 , 0.0 } ,
            { 0.0 , 1.0 , 0.0 } ,
            { 0.0 , 0.0 , 1.0 } ,
            1000 , 500 ,
            -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  zimt::process ( a3.shape , sters ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  ster_to_ray_t < double , LANES > ster_to_ray ;

  zimt::process ( b3.shape , ls ,
                  ster_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto crd2 = to_md ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      ster_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  stereographic_stepper < double , LANES , false >
    rsters ( ex , ey , ez , 1000 , 500 ,
            -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  zimt::process ( a3.shape , rsters ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  zimt::process ( b3.shape , ls , ster_to_ray + r3 ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto d = abs ( a3 [ { x , y } ] - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  // for cubemap projection, we use differently-shaped
  // arrays - using a shpe which would be natural for a
  // 'standard' 1:6 vertically-stacked cubemap. Apart from
  // parameterization, the process is the same as for the
  // other projections.

  zimt::array_t < 2 , d3_t > a6 ( { 500 , 3000 } ) ;
  zimt::array_t < 2 , d3_t > b6 ( { 500 , 3000 } ) ;

  cubemap_stepper < double , LANES , false >
    cbms ( { 1.0 , 0.0 , 0.0 } ,
           { 0.0 , 1.0 , 0.0 } ,
           { 0.0 , 0.0 , 1.0 } ,
           500 , 3000 ,
           -1.0 , 1.0 , -6.0 , 6.0 ) ;

  zimt::process ( a6.shape , cbms ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a6 ) ) ;

  ir_to_ray_t < double , LANES > ir_to_ray ;

  x0 = ( a6.shape[0] - 1 ) / -2.0 ;
  y0 = ( a6.shape[1] - 1 ) / -2.0 ;
  dx = 2.0 / a6.shape[0] ;
  dy = 12.0 / a6.shape[1] ;
  x0 *= dx ;
  y0 *= dy ;

  zimt::linspace_t < double , 2 , 2 , LANES >
    ls2 ( { x0 , y0 } , { dx , dy } ) ;
           
  zimt::process ( b6.shape , ls2 ,
                  ir_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b6 ) ) ;

  auto to_md2 = [&] ( std::size_t sx , std::size_t sy )
  {
    zimt::xel_t < double , 2 > crd2 ;
    double & x ( crd2[0] ) ;
    double & y ( crd2[1] ) ;
    x = sx - ( a6.shape[0] - 1 ) / 2.0 ;
    y = sy - ( a6.shape[1] - 1 ) / 2.0 ;
    x *= 2.0 / a6.shape[0] ;
    y *= 12.0 / a6.shape[1] ;
    return crd2 ;
  } ;

  for ( std::size_t y = 0 ; y < a6.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a6.shape[0] ; x++ )
    {
      auto crd2 = to_md2 ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      ir_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a6 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b6 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  cubemap_stepper < double , LANES , false >
    rcbms ( ex , ey , ez , 500 , 3000 ,
           -1.0 , 1.0 , -6.0 , 6.0 ) ;

  zimt::process ( a6.shape , rcbms ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a6 ) ) ;

  zimt::process ( b6.shape , ls2 , ir_to_ray + r3 ,
                  zimt::storer < double , 3 , 2 , LANES > ( b6 ) ) ;

  for ( std::size_t y = 0 ; y < a6.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a6.shape[0] ; x++ )
    {
      auto d = abs ( a6 [ { x , y } ] - b6 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  // for biatan6 projection, we use the same data as for cubemap.

  biatan6_stepper < double , LANES , false >
    ba6s ( { 1.0 , 0.0 , 0.0 } ,
           { 0.0 , 1.0 , 0.0 } ,
           { 0.0 , 0.0 , 1.0 } ,
           500 , 3000 ,
           -1.0 , 1.0 , -6.0 , 6.0 ) ;

  zimt::process ( a6.shape , ba6s ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a6 ) ) ;

  ba6_to_ray_t < double , LANES > ba6_to_ray ;

  x0 = ( a6.shape[0] - 1 ) / -2.0 ;
  y0 = ( a6.shape[1] - 1 ) / -2.0 ;
  dx = 2.0 / a6.shape[0] ;
  dy = 12.0 / a6.shape[1] ;
  x0 *= dx ;
  y0 *= dy ;

  zimt::process ( b6.shape , ls2 ,
                  ba6_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b6 ) ) ;

  for ( std::size_t y = 0 ; y < a6.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a6.shape[0] ; x++ )
    {
      auto crd2 = to_md2 ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      ba6_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a6 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b6 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  biatan6_stepper < double , LANES , false >
    rba6s ( ex , ey , ez , 500 , 3000 ,
           -1.0 , 1.0 , -6.0 , 6.0 ) ;

  zimt::process ( a6.shape , rba6s ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a6 ) ) ;

  zimt::process ( b6.shape , ls2 , ba6_to_ray + r3 ,
                  zimt::storer < double , 3 , 2 , LANES > ( b6 ) ) ;

  for ( std::size_t y = 0 ; y < a6.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a6.shape[0] ; x++ )
    {
      auto d = abs ( a6 [ { x , y } ] - b6 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  // if none of the assertions failed and the program terminates
  // now, all is well and the functors work as expected.

  return 0 ;
}

// Here, we define the SIMD-ISA-specific derived 'dispatch' class:

struct dispatch
: public dispatch_base
{
  // We fit the derived dispatch class with a c'tor which fills in
  // information about the nested SIMD ISA we're currently in.

  dispatch()
  {
    backend = int ( zimt::simdized_type<int,4>::backend ) ;
    #if defined USE_HWY || defined MULTI_SIMD_ISA
      hwy_isa = HWY_TARGET ;
    #endif
  }

  // 'payload', the SIMD-ISA-specific overload of dispatch_base's
  // pure virtual member function, now has the code which was in
  // main() when this example was first coded without dispatch.
  // One might be more tight-fisted with which part of the former
  // 'main' should go here and which part should remain in the
  // new 'main', but the little extra code which wouldn't benefit
  // from vectorization doesn't make much of a difference here.
  // Larger projects would have both several payload-type functions
  // and a body of code which is independent of vectorization.

///////////////// write a payload function with a 'main' signature
  
  int payload ( int argc , char * argv[] ) const
  {
    // we can get information about the specific dispatch object:

    std::cout << "payload code is using back-end: "
              << zimt::backend_name [ backend ] << std::endl ;

    #if defined USE_HWY || defined MULTI_SIMD_ISA

    std::cout << "highway target: "
              << hwy::TargetName ( hwy_isa ) << std::endl ;

    #endif

    _payload ( argc , argv ) ;
    return 0 ;
  }
} ;

//--------------------------------------------------------------------

// we also code a local function _get_dispatch which returns a pointer
// to 'dispatch_base', which points to an object of the derived class
// 'dispatch'. This is used with highway's HWY_DYNAMIC_DISPATCH and
// returns the dispatch pointer for the SIMD ISA which highway deems
// most appropriate for the CPU on which the code is currently running.

const dispatch_base * const _get_dispatch()
{
  static dispatch d ;
  return &d ;
}

END_ZIMT_SIMD_NAMESPACE

HWY_AFTER_NAMESPACE() ;

// Now for code which isn't SIMD-ISA-specific. ZIMT_ONCE is defined
// as either HWY_ONCE (if MULTI_SIMD_ISA is #defined) or simply true
// otherwise - then, there is only one compilation anyway.

#if ZIMT_ONCE

namespace project {

#ifdef MULTI_SIMD_ISA

// we're using highway's foreach_target mechanism. To get access to the
// SIMD-ISA-specific variant of _get_dispatch (in project::HWY_NAMESPACE)
// we use the HWY_EXPORT macro:

HWY_EXPORT(_get_dispatch);

// now we can code get_dispatch: it simply uses HWY_DYNAMIC_DISPATCH
// to pick the SIMD-ISA-specific get_dispatch variant, which in turn
// yields the desired dispatch_base pointer.

const dispatch_base * const get_dispatch()
{
  return HWY_DYNAMIC_DISPATCH(_get_dispatch)() ;
}

#else // #ifdef MULTI_SIMD_ISA

// if we're not using highway's foreach_target mechanism, there is
// only a single _get_dispatch variant in namespace project::zsimd.
// So we call that one, to receive the desired dispatch_base pointer.

const dispatch_base * const get_dispatch()
{
  return zsimd::_get_dispatch() ;
}

#endif // #ifdef MULTI_SIMD_ISA

}  // namespace project

int main ( int argc , char * argv[] )
{
  // Here we use zimt's dispatch mechanism: first, we get a pointer
  // to the dispatcher, then we invoke a member function of the
  // dispatcher. What's the point? We can call a SIMD-ISA-specific
  // bit of code without having to concern ourselves with figuring
  // out which SIMD ISA to use on the current CPU: this happens via
  // highway's dispatch mechanism, or is fixed at compile time, but
  // in any case we receive a dispatch_base pointer routing to the
  // concrete variant. project::get_dispatch might even be coded
  // to provide pointers to dispatch objects in separate TUs, e.g.
  // when these TUs use different back-ends or compiler flags. Here,
  // we can remain unaware of how the concrete dispatch object is
  // set up and the pointer obtained.

  auto dp = project::get_dispatch() ;

  // we can get information about the specific dispatch object:

  std::cout << "calling payload with "
#ifdef MULTI_SIMD_ISA
            << "dynamic dispatch" << std::endl ;
#else
            << "static dispatch" << std::endl ;
#endif

  // now we call the payload via the dispatch_base pointer.

  int success = dp->payload ( argc , argv ) ;

  std::cout << std::endl ;

// to invoke other variants with explicit calls, uncomment this section
// (x86 only - other CPU types would need other N_... macros)

#ifdef MULTI_SIMD_ISA

  std::cout << "******** calling payloads with successively better ISAs:"
            << std::endl ;
  
  std::cout << "******** calling SSE2 payload" << std::endl ;
  project::N_SSE2::dispatch d3 ;
  d3.payload ( argc , argv ) ;
  
  std::cout << "******** calling SSSE3 payload" << std::endl ;
  project::N_SSSE3::dispatch d4 ;
  d4.payload ( argc , argv ) ;
  
  std::cout << "******** calling SSE4 payload" << std::endl ;
  project::N_SSE4::dispatch d5 ;
  d5.payload ( argc , argv ) ;
  
  std::cout << "******** calling AVX2 payload" << std::endl ;
  project::N_AVX2::dispatch d6 ;
  d6.payload ( argc , argv ) ;
  
  std::cout << "******** calling AVX3 payload" << std::endl ;
  project::N_AVX3::dispatch d7 ;
  d7.payload ( argc , argv ) ;

#endif
}

#endif  // ZIMT_ONCE

