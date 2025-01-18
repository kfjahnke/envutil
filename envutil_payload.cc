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


#include "zimt/common.h"

#ifdef MULTI_SIMD_ISA

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

// now we #include highway.h - as we would do after foreach_target.h
// in a multi-ISA build. With foreach_target.h, the code is re-included
// several times, each time with a different ISA. Here we have set one
// specific ISA and there won't be any re-includes.

#include <hwy/highway.h>

#else // MULTI_SIMD_ISA

#define HWY_TARGET TG_ISA

#endif // MULTI_SIMD_ISA

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

// this asset_t object is a handle to the current 'enviromnet'
// providing image data

static zimt::asset_t current_env ;

// envutil_dispatch.h has the definition of class dispatch_base and
// of the function get_dispatch.

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

    Imath::Eulerf angles ( roll , pitch , yaw , Imath::Eulerf::ZXY ) ;
    q = angles.toQuat() ;

    if ( inverse )
      q.invert() ;
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

  std::chrono::system_clock::time_point start
    = std::chrono::system_clock::now() ;
  
  zimt::process ( trg.shape , get , act , cstor , bill ) ;
  
  std::chrono::system_clock::time_point end
    = std::chrono::system_clock::now() ;

  auto msec = std::chrono::duration_cast<std::chrono::milliseconds>
                ( end - start ) . count() ;

  if ( args.verbose )
  {
    std::cout << "frame rendering time: " << msec << " ms" << std::endl ;
  }

  rt_cumulated += msec ;

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
  typedef environment < float , float , NCH , 16 > env_t ;
  env_t * p_env = nullptr ;

  if constexpr ( NCH == 0 )
  {
    // just a stub:
    std::cout << "processing multi-facet file " << args.fct_file
              << " not yet implemented!" << std::endl ;
    return ;

    // TODO: extend class 'environment' with a c'tor taking a
    // facet file. The resulting object should behave like any
    // other environment object, but internally use several
    // sub-environments for the mounted facets and logic to
    // form a weighted sum of pixels from the contributing
    // facets.

    // just to show how the assets could be set up:

    typedef asset_t < std::string > named_asset_t ;
    auto p_image_v = new std::vector < std::shared_ptr < named_asset_t > >  ;
    
    std::string img_name1 = "dummy1" ;
    std::string img_name2 = "dummy2" ;

    // these pointers would be made to point to image data - right now we
    // just need something in dynamic memory:

    std::string * p_str1 = new std::string ( "cargo1" ) ;
    std::string * p_str2 = new std::string ( "cargo2" ) ;

    // TODO: alternatively hold e.g. pointers to class environment.

    auto p_image1 = std::make_shared < named_asset_t > ( img_name1 , p_str1 ) ;
    auto p_image2 = std::make_shared < named_asset_t > ( img_name2 , p_str2 ) ;

    // the shared_ptrs can be stored in a vector

    p_image_v->push_back ( p_image1 ) ;
    p_image_v->push_back ( p_image2 ) ;

    // after scanning the facet file, image_v should contain shared pointers
    // to all the assets corresponding to images files in the facet file.
    // Then, corresponding steppers can be set up, and the set of steppers
    // combined in a confluence_t object, with a synopsis function which
    // produces e.g. a voronoid by evaluating corresponding interpolators
    // and prioritizing them by distance-from-center.

    // image_v can be made the current_env by clearing current_env and
    // placing a pointer to a dynamically created asset vector in it, with
    // the facet file's name as ID:

    current_env.reset ( args.fct_file , p_image_v ) ;

    // now this multi-image environment will persist until a new environment
    // is set up or unti the program ends. environment-wise: there can only be
    // one!
  }
  else
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

    rotate_3d < double , 16 > r3 ( args.roll ,
                                   args.pitch ,
                                   args.yaw ,
                                   false ) ;

    xx = r3 ( xx ) ;
    yy = r3 ( yy ) ;
    zz = r3 ( zz ) ;

    // if we have a 'facet' - a mounted image with non-standard
    // orientation - we add a second rotation, now with the inverse
    // of the quaternion formed from the orientation Euler angles.
    // Why inverse? If the facet is oriented precisely like the
    // virtual camera, both rotations should cancel each other out
    // precisely.

    if (   args.mount_yaw != 0.0f
        || args.mount_pitch != 0.0f
        || args.mount_roll != 0.0f )
    {
      rotate_3d < double , 16 > rr3 ( args.mount_roll ,
                                      args.mount_pitch ,
                                      args.mount_yaw ,
                                      true ) ;
      xx = rr3 ( xx ) ;
      yy = rr3 ( yy ) ;
      zz = rr3 ( zz ) ;
    }

    if ( args.itp == 1 || args.itp == -2 )
    {
      // both direct spline interpolation and twining need an
      // 'environment' object.

      // create an 'environment' object. This will persist while the
      // input image remains the same, so if we're creating an image
      // sequence, it will be reused for each individual image.
      // we 'hold' the environment object in current_env, an 'asset'
      // object, which destructs and dealocates it's client object
      // only when itself destructs or is 'cleared'.

      p_env = (env_t*) current_env.has ( args.input ) ;

      if ( ! p_env )
      {
        current_env.clear() ;
        p_env = new env_t() ;
        current_env.reset ( args.input , p_env ) ;
      }
    }

    // There are two code paths: one taking the getters yielding
    // simple single-ray 3D coordinates, and one taking the three-ray
    // variant needed to compute the derivatives. They use specific
    // types of 'environment' objects. The code for the environment
    // objects is in environment.h

    if ( ninputs == 3 )
    {
      // set up a simple single-coordinate stepper of the type
      // fixed by 'STP'. This route is taken with direct b-spline
      // interpolation (itp 1)

      STP < float , 16 > get_ray
        ( xx , yy , zz , args.width , args.height ,
          args.x0 , args.x1 , args.y0 , args.y1 ) ;

      // now we call the final 'work' template which uses the get_t
      // and act objects we've set up so far

      work ( get_ray , *p_env ) ;
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

      // we pass p_env - a pointer to an environment object.
      // if itp is set to -1 (use OIIO) this pointer is nullptr.

      environment9 < NCH , 16 > env ( p_env ) ;

      work ( get_ray , env ) ;
    }
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
    case BIATAN6:
      roll_out < NCH , biatan6_stepper > ( ninputs ) ;
      break ;
    default:
      break ;
  }
}

// here we have the implementation of the derived class 'dispatch',
// and of _get_dispatch. 

typedef zimt::xel_t < float , 2 > crd_t ;
typedef zimt::xel_t < float , 3 > px_t ;

struct fake_eval
: public zimt::unary_functor < crd_t , px_t , 16 >
{
  template < typename in_t , typename out_t >
  void eval ( const in_t & in , out_t & out )
  {
    out[0] = in[0] ;
    out[1] = in[1] ;
    out[2] = in[0] + in[1] ;
  }
} ;

struct dispatch
: public dispatch_base
{
  // We fit the derived dispatch class with a c'tor which fills in
  // information about the nested SIMD ISA we're currently in.

  dispatch()
  {
    hwy_target = HWY_TARGET ;
    #ifdef HWY_TARGET_STR
      hwy_target_name = hwy::TargetName ( HWY_TARGET ) ;
      hwy_target_str = HWY_TARGET_STR ;
    #else
      hwy_target_name = "unknown" ;
      hwy_target_str = "unknown" ;
    #endif
  }

  int payload ( int nchannels ,
                int ninputs ,
                projection_t projection ) const
  {
    switch ( nchannels )
    {
      case 0:
        roll_out < 0 > ( ninputs , projection ) ;
        break ;
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

    /*
    typedef zimt::confluence_t < float , 3 , 2 , 16 , float , 2 > cf_t ;
    typedef zimt::loader < float , 2 , 2 , 16 > get_t ;
    typedef zimt::storer < float , 3 , 2 , 16 > put_t ;

    typedef zimt::pass_through < float , 3 , 16 > pass_t ;

    typedef zimt::simdized_type < crd_t , 16 > crd_v ;
    typedef zimt::simdized_type < px_t , 16 > px_v ;

    zimt::bill_t bill ;
    zimt::array_t < 2 , crd_t > src1 ( { 100 , 100 } ) ;
    zimt::array_t < 2 , crd_t > src2 ( { 100 , 100 } ) ;
    zimt::array_t < 2 , px_t > trg ( { 100 , 100 } ) ;

    for ( int y = 0 ; y < 100 ; y++ )
    {
      for ( int x = 0 ; x < 100 ; x++ )
      {
        src1 [ { x , y } ] = x + 100 * y * 3 ;
        src2 [ { x , y } ] = x + 100 * y ;
      }
    }
    src1 [ { 1 , 1 } ] = { .3f , .4f } ;

    typedef zimt::grok_get_t < float , 2 , 2 , 16 > gg_t ;
    std::vector < gg_t > get_v ;
    get_v.push_back ( get_t ( src1 , bill ) ) ;
    get_v.push_back ( get_t ( src2 , bill ) ) ;

    std::vector < fake_eval > evv ( 2 ) ;

    cf_t cf ( get_v ,
              [&evv] ( const cf_t & cfr ,
                       px_v & out ,
                       const std::size_t & cap )
              {
                // 'voronoi' function: of the set of coordinates produced
                // by the grok_get_t, pick the one with the smallest
                // squared norm, and evaluate the corresponding evaluator.
                auto peak = squared_norm ( cfr.src_v[0] ) ;
                zimt::simdized_type < int , 16 > winner ( 0 ) ;
                // evv[0].eval ( cfr.src_v[0] , out ) ;
                for ( std::size_t i = 1 ; i < cfr.sz ; i++ )
                {
                  auto sqn = squared_norm ( cfr.src_v[i] ) ;
                  auto mask = sqn < peak ;
                  if ( any_of ( mask ) )
                  {
                    peak ( mask ) = sqn ;
                    winner ( mask ) = i ;
                  }
                }
                for ( std::size_t i = 0 ; i < cfr.sz ; i++ )
                {
                  auto mask = ( winner == i ) ;
                  if ( any_of ( mask ) )
                  {
                    px_v help ;
                    evv[i].eval ( cfr.src_v[i] , help ) ;
                    out[0] ( mask ) = help[0] ;
                    out[1] ( mask ) = help[1] ;
                    out[2] ( mask ) = help[2] ;
                  }
                }
              } ) ;

    pass_t ps ;
    put_t p ( trg ) ;

    zimt::process ( src1.shape , cf , ps , p , bill ) ;

    std::cout << trg [ { 0 , 1 } ] << std::endl ;
    std::cout << trg [ { 1 , 1 } ] << std::endl ;
    */

    return 0 ;
  }
} ;

const dispatch_base * const _get_dispatch()
{
  static dispatch d ;
  return &d ;
}

END_ZIMT_SIMD_NAMESPACE
HWY_AFTER_NAMESPACE() ;
