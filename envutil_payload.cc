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

// this file has the dispatch to ISA-specific 'payload' code.

#include <iostream>

// This header defines all the macros having to do with targets:

#include <hwy/detect_targets.h>

// glean the target as 'TG_ISA' from outside - this file is intended
// to produce ISA-specific separate TUs containing only binary for
// one given ISA, but assuming that other files of similar structure,
// but for different ISAs will also be made and all linked together
// with more code which actually makes use of the single-ISA TUs.
// 'Slotting in' the target ISA from the build system is enough to
// produce a SIMD-ISA-specific TU - all the needed specifics are
// derived from this single information. detect_targets.h sets
// HWY_TARGET to HWY_STATIC_TARGET, so we #undef it and use the
// target specification from outside instead.

#undef HWY_TARGET
#define HWY_TARGET TG_ISA
#undef HWY_NAMESPACE
#define HWY_NAMESPACE NS_ISA

// now we #include highway.h - as we would do after foreach_target.h
// in a multi-ISA build. With foreach_target.h, the code is re-included
// several times, each time with a different ISA. Here we have set one
// specific ISA and there won't be any re-includes.

#include <hwy/highway.h>

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

// we #include xel.h in the ISA-specific namespace. It is guarded
// by a conventional sentinel, so further re-inclusions have no
// effect and the ISA environment is fixed for the entire TU.
// same for array.h.

#include "zimt/xel.h"
#include "zimt/array.h"

HWY_AFTER_NAMESPACE() ;

// envutil_dispatch.h has the definition of class dispatch_base and
// of the function get_dispatch.

#include "zimt/common.h"
#include "zimt/zimt.h"
#include "envutil_dispatch.h"

// environment.h has most of the 'workhorse' code for this program.
// it provides functional constructs to yield pixels for coordinates.
// These functors are used with zimt::process to populate the output.

#include "environment.h"
#include "stepper.h"

// now we invoke HWY_BEFORE_NAMESPACE() for this file, to make code
// down to HWY_AFTER_NAMESPACE() compile with settings for a specific
// target (via e.g. #pragmas to the compiler)

HWY_BEFORE_NAMESPACE() ;

BEGIN_ZIMT_SIMD_NAMESPACE(project)

// // we define 'rollout', which has the entry into actual payload code
// 
// int rollout ( int nchannels ,
//               int ninputs ,
//               projection_t projection )
// {
//   // just a stub for now
// 
//   std::cout << "call to rollout, target = "
//             << hwy::TargetName ( HWY_TARGET )
//             << std::endl ;
// 
//   return 0 ;
// }


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

// 'work' is the function where the state we've built up in the
// code below it culminates in the call to zimt::process, which
// obtains pick-up coordinates from the 'get' object, produces
// pixels for these coordinates with the 'act' object and stores
// these pixels with a zimt::storer, which is set up in this
// function. All the specificity of the code has now moved to
// the types get_t and act_t - the roll_out is complete and we
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

// to call the appropriate instantiation of 'work' (above)
// we need to do some dispatching: picking types depending
// on run-time variables. We achieve this with a staged 'roll_out'
// routine: every stage processes one argument and routes to
// specialized code. The least specialized roll_out variant is
// the lowest one down, this here is the final stage where we have
// the number of channels and the stepper as template arguments.
// Here we proceed to set up more state which is common to all
// code paths, set up the 'act' functor which yields pixels, and
// finally call 'work' to run the pixel pipeline.

template < int NCH ,
           template < typename , std::size_t > class STP >
void roll_out ( int ninputs )
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
// roll_out below, now we roll_out on the projection and instantiate
// the next roll_out level with a stepper type which fits the
// projection.

template < int NCH >
void roll_out ( int ninputs ,
                projection_t projection )
{
  switch ( projection )
  {
    case SPHERICAL:
      roll_out < NCH , spherical_stepper > ( ninputs ) ;
      break ;
    case CYLINDRICAL:
      roll_out < NCH , cylindrical_stepper > ( ninputs ) ;
      break ;
    case RECTILINEAR:
      roll_out < NCH , rectilinear_stepper > ( ninputs ) ;
      break ;
    case FISHEYE:
      roll_out < NCH , fisheye_stepper > ( ninputs ) ;
      break ;
    case STEREOGRAPHIC:
      roll_out < NCH , stereographic_stepper > ( ninputs ) ;
      break ;
    case CUBEMAP:
      roll_out < NCH , cubemap_stepper > ( ninputs ) ;
      break ;
    default:
      break ;
  }
}

// roll_out by the number of colour channels. We process one to four,
// where the usefulness of two channels isn't clear.

void roll_out ( int nchannels ,
                int ninputs ,
                projection_t projection )
{
  switch ( nchannels )
  {
    case 1:
      roll_out < 1 > ( ninputs , projection ) ;
      break ;
    case 2:
      roll_out < 2 > ( ninputs , projection ) ;
      break ;
    case 3:
      roll_out < 3 > ( ninputs , projection ) ;
      break ;
    case 4:
      roll_out < 4 > ( ninputs , projection ) ;
      break ;
  }
}

END_ZIMT_SIMD_NAMESPACE
HWY_AFTER_NAMESPACE() ;
