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

  rotate_3d ( T roll = 0 , T pitch = 0 , T yaw = 0 , bool inverse = false )
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
  bill.njobs = 1 ;

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

// this class is an example for giving preference to one specific
// facet - as opposed to code mixing several values. There are two
// member functions. The first one processes a std::vector of simdized
// rays, selecting the rays which are most central (by calculating the
// angle with the 'forward' ray (0,0,1)) and evaluating there. The
// angle can be calculated in this way because the rays we receive
// here are in the facets' coordinate system - the orientations of
// the facets and the virtual camera are encoded in the 'steppers'
// which generate the rays we receive here. The second member function
// receives a std::vector of 'ninepacks' which are the basis for using
// 'twining'. It does in turn call the first member repeatedly - once
// for each of the twining kernel's coefficients. Note how this softens
// the hard discontinuity between 'colliding' facets: each of the
// component values of the weighted sum which is formed by the twining
// operations may stem from a *different* facet, so near the border,
// result pixels will contain a *blend* of components from both facets.
// While twining is computationally expensive, this is a strategy to
// avoid the ungainly staircase artifacts where facets collide. The
// final result is like a scaled-down version of a much larger image,
// and the scaling-down would have the same effect of blending pixels
// from colliding facets where pixels from two facets fall into the
// same bin in the scaled-down result.
// The result is a voronoi diagram - apart from the fact that where
// facets don't cover their full share of their tangent plane, other
// facets may be visible. If all facets are sufficiently large and
// there are no gaps, the result is much like a voronoi diagram. envutil
// also allows facets in other projections than rectilinear - with
// such facets, the image center is also the reference but the samples
// are - conceptually - no longer on a tangent plane but on the 'draped'
// 2D manifold - e.g. a sphere for spherically-projected facets. Still,
// the proximity criterion is calculated with 3D rays - but where
// rectilinear facets have straight intersections, non-rectilinear
// ones may have curved intersections. The ray-based proximity criterion
// is slower to compute than the 2D distance-to-center one might also
// use (on target coordinates) - but ray-based is universally applicable
// and produces values which compare across different projections,
// whereas the 2D target coordinate varies with the projection used.

// The code in voronoi_syn is only suitable for fully-opaque images:
// to handle data with transparency, figuring out a single 'winner'
// is not sufficient, because the winning pixel might be partially or
// fully transparent, in which case worse-ranking pixels should 'shine
// through'. To achieve this, z-buffering of the data is needed - like
// it is done in lux (see pv_rendering.cc, class 'multi_facet_type_4').
// class voronoi_syn_plus implements handling of transparent facets
// with alpha compositing in z-buffer sequence, see there!
// The lux code also shows that even with transparency, a lot of
// potential shortcuts can save time and not every pixel needs the
// 'full treatment'. lux also employs facet pre-selection, inspecting
// each facet's frustum and it's intersection with the target image's.
// This process of 'culling' can greatly reduce computational load if
// the current view only involves a small amount of facets. When
// stitching panoramas, though, it's often the case that all facets are
// involved, so the frustum inspection code does not help. In envutil,
// I intend to use a different method: coarse masking and tiling. The
// idea is to use a coarse mask of the intersection of each facet with
// the target image stored in a cubemap (preferably biatan6 type) with
// relatively low resolution which coincides with a tiling of the target
// image so that each pixel in the coarse mask is 'as large' as a tile
// of the target image. If, then, the target image is built tile by tile,
// the participating facets can be found by evaluating their corresponding
// coarse masks and only picking those which aren't zero at the given
// location. A biatan6-projected cubemap is quite close to ideal when it
// comes to evenness of represented spherical 'real estate' per pixel
// of the map. The operation would be more fine-grained compared to the
// frustum-based approach.

// This aside, here goes: this is the code for fully opaque facets:

template < typename ENV >
struct voronoi_syn
{
  typedef typename ENV::px_t px_t ;
  typedef simdized_type < px_t , 16 > px_v ;
  static const std::size_t nch = px_t::size() ;
  const int sz ;
  typedef zimt::xel_t < float , 3 > ray_t ;
  typedef zimt::xel_t < float , 9 > ninepack_t ;
  typedef simdized_type < float , 16 > f_v ;
  typedef simdized_type < ray_t , 16 > ray_v ;
  typedef simdized_type < ninepack_t , 16 > ninepack_v ;

  std::vector < ENV > & env_v ;

  // 'scratch' will hold all ray packets which will react with
  // a given twining coefficient 'cf'
  
  std::vector < ray_v > scratch ;


  voronoi_syn ( std::vector < ENV > & _env_v )
  : env_v ( _env_v ) ,
    sz ( _env_v.size() ) ,
    scratch ( _env_v.size() )
  { }

  // general calculation of the angles between rays:
  //
  // f_v angle ( const ray_v & a , const ray_v & b ) const
  // {
  //   auto dot = ( a * b ) . sum() ;
  //   auto costheta = dot / ( norm ( a ) * norm ( b ) ) ;
  //   return acos ( costheta ) ;
  // }

  // this helper function calculates the vector of angles between
  // a vector of rays and (0,0,1) - unit forward in our CS.

  // TODO: many steppers can be modified to produce normalized rays
  // without additional computational cost, and even if there is some
  // cost, it should not amount to more than the division by the norm
  // which we use here to obtain the angle. With incoming normalized
  // rays, it would be sufficient to look at the z component of the
  // incoming ray - due to the pre-rotation, z will fall from +1 with
  // true forward rays to -1 to rays pointing straigt back. So the
  // facet where the ray has the largest z coordinate is the winner.
  // ----> done.
  // .. but: normalizing the rays destroys some information, which
  // is needed e.g. when applying translation parameters (or so it
  // seems). Most steppers do produce normalized rays, so for these
  // the max_z or angle_against_001 criterion could be obtained
  // without calculation from the z coordinate. Only steppers where
  // the ray isn't automatically normalized would need a caclulation.
  // Maybe the value can be produced by the stepper?

  // f_v angle_against_001 ( const ray_v & a ) const
  // {
  //   return acos ( a[2] / norm ( a ) ) ;
  // }

  // this is the operator() overload for incoming ray data. The next
  // one down is for 'ninepacks'.

  void operator() (
        const std::vector < ray_v > & pv ,
        px_v & trg ,
        const std::size_t & cap = 16 ) const
  {
    // this vector will contain the index of the 'most suitable'
    // facet at the given location, according to the criterion.

    simdized_type < int , 16 > champion_v ( -1 ) ;

    // we initialize 'max_z' to a very small value - no
    // 'real' z value can ever be this small.

    simdized_type < float , 16 >
      max_z ( std::numeric_limits<float>::min() ) ;

    // next_best starts out with an invalid facet index, and if this
    // isn't updated during processing, it serves as an indicator
    // that no facet is hit by any ray - a special case which can
    // be dealt with very efficiently: produce 0000

    int next_best = -1 ;

    // test which of the rays in pv[0] pass through facet 0

    auto valid = env_v[0].get_mask ( pv[0] ) ;

    // do any of the rays hit the facet? then update 'next_best',
    // which, at this point, is still -1, indicating 'all misses'.
    // This can't be the case any more. Also update max_z:
    // overwrite the initial 'impossible' value with a real z
    // value where the rays hit the facet.

    if ( any_of ( valid ) )
    {
      next_best = 0 ;
      champion_v = 0 ;
      max_z ( valid ) = pv[0][2] ;
    }

    // now go through the facets in turn and replace the current
    // champions where the current facet provides better candidates

    for ( int i = 1 ; i < sz ; i++ )
    {
      // find the mask for facet i - if no rays hit the facet,
      // go to the next loop iteration straight away

      valid = env_v[i].get_mask ( pv[i] ) ;
      if ( none_of ( valid ) )
        continue ;

      // we initialize 'max_z' with an impossibly small value
      // and 'patch in' z values where rays hit the facet.
      // 'max_z' encodes our 'quality' criterion: the largest z
      // we find is the best.

      simdized_type < float , 16 >
        current_z ( std::numeric_limits<float>::min() ) ;

      current_z ( valid ) = pv[i][2] ;

      // now we test whether current z values are larger than the
      // largest we've seen so far

      auto mask = ( current_z > max_z ) ;

      // are any current z larger? then update 'next_best',
      // which, at this point, may still be -1, or already a value
      // set by a previous loop iteration. Also update max_z,
      // and, finally, update the vector of champions where new
      // maximal z values were found to the index of the facet we're
      // looking at now.

      if ( any_of ( mask ) )
      {
        next_best = i ;
        max_z ( mask ) = current_z ;
        champion_v ( mask ) = i ;
      }
    }

    // first check: if next_best is still -1, all facets were missed
    // and the result is 0000. This is the fastest outcome, and if
    // part of the output isn't covered by any facets, this will
    // save a good deal of cycles.

    if ( next_best == -1 )
    {
      trg = 0.0f ;
    }

    // now check: if all occupants of 'champion_v' are equal,
    // we can special-case. Since this is the most common case,
    // this is beneficial for performance. Note that we don't
    // check the contents of champion_v for equality - next_best
    // certainly holds a facet index to a facet which *is* hit,
    // (it isn't -1, so it is set to a 'hit' facet's index)
    // and if all champions are equal, they must be equal to this
    // one, which is the only one in-play in this special case.

    else if ( all_of ( champion_v == next_best ) )
    {
      env_v [ next_best ] . eval ( pv [ next_best ] , trg ) ;
    }
    else
    {
      // champion_v's occupants are not all equal, so we need to
      // evaluate every in-play facet and compose the result.
      // some rays may not have hit any facet, so we clear trg
      // first. We refrain from testing for occupants with value
      // -1, because forming a mask and performing a masked assignment
      // is more costly then clearing trg. this is also the 'worst
      // case' - preformance-wise - where the rays hit several facets. 

      trg = 0.0f ;

      for ( int i = 0 ; i < sz ; i++ )
      {
        // form a mask which is true where the champion is i.

        auto mask = ( champion_v == i ) ;

        // if this mask has occupants, we need to obtain pixel data
        // for this facet and transfer them to the corresponding lanes.

        if ( any_of ( mask ) )
        {
          px_v help ;
          env_v [ i ] . eval ( pv [ i ] , help ) ;
          trg ( mask ) = help ;
        }
      }
    }
  }

  // operator() for 'ninepacks'. This implements 'twining': evaluation
  // at several rays in the vicinity of the central ray and formation
  // of a weighted sum of these evaluations. The incoming 'ninepack'
  // holds three concatenated vectorized rays: the 'central' rays,
  // the rays which result from a discrete target coordinate set off
  // by one to the right, and the rays which result from a discrete
  // target coordinate set off one down. By differencing, we obtain
  // two vectors representing 'canonical' steps in the target image's
  // horizontal and vertical. We use these two vectors as a basis for
  // the 'twining' operation: each twining coefficient has two factors
  // describing it's spatial aspect, and multiplying them with the
  // corresponding basis vectors and summing up produces the offset
  // we need to add to the 'central' or 'pickup' ray to obtian a
  // 'deflected' ray. The substrate is evaluated at each deflected
  // ray, and the results are combined in a weighted sum, where the
  // weight comes from the third factor in the twining coefficient.

  void operator() ( const std::vector < ninepack_v > & pv ,
                    px_v & trg ,
                    const std::size_t & cap = 16 )
  {
    typedef simdized_type < float , 16 > f_v ;
    typedef simdized_type < zimt::xel_t < float , 3 > , 16 > ray_v ;
    typedef simdized_type < zimt::xel_t < float , 9 > , 16 > ray9_v;

    // trg is built up as a weighted sum of contributing values
    // resulting from the reaction with a given coefficient, so we
    // clear it befor we begin.
  
    trg = 0.0f ;

    // for each coefficient in the twining 'spread'

    for ( const auto & cf : args.twine_spread )
    {
      for ( std::size_t facet = 0 ; facet < sz ; facet++ )
      {
        const ray9_v & in ( pv[facet] ) ; // shorthand

        ray_v p0 { in[0] ,         in[1] ,         in[2] } ;
        ray_v du { in[3] - in[0] , in[4] - in[1] , in[5] - in[2] } ;
        ray_v dv { in[6] - in[0] , in[7] - in[1] , in[8] - in[2] } ;

        scratch [ facet ] = p0 + cf[0] * du + cf[1] * dv ;
      }

      // now we have sz entries in scratch, and we can call
      // the three-component form above to produce a synopsis
      // for the contributors reacting with the current coefficient

      px_v help ;
      operator() ( scratch , help , cap ) ;

      // 'help' is the pixel value, which is weighted with the
      // weighting value in the twining kernel and added to what's
      // in trg already.

      trg += cf[2] * help ;
    }  
    // and that's it - trg has the result.
  }
} ;

// variant with z-buffering for facets with alpha channel:

template < typename ENV , bool fix_duv = true >
struct voronoi_syn_plus
{
  typedef typename ENV::px_t px_t ;
  typedef simdized_type < px_t , 16 > px_v ;
  static const std::size_t nch = px_t::size() ;
  const int sz ;
  typedef zimt::xel_t < float , 3 > ray_t ;
  typedef zimt::xel_t < float , 9 > ninepack_t ;
  typedef simdized_type < float , 16 > f_v ;
  typedef simdized_type < ray_t , 16 > ray_v ;
  typedef simdized_type < ninepack_t , 16 > ninepack_v ;
  typedef simdized_type < int , 16 > index_v ;

  std::vector < ENV > & env_v ;
  std::vector < ray_v > scratch ;

  voronoi_syn_plus ( std::vector < ENV > & _env_v )
  : env_v ( _env_v ) ,
    sz ( _env_v.size() ) ,
    scratch ( _env_v.size() )
  { }

  // // this helper function calculates the vector of angles between
  // // a vector of rays and (0,0,1) - unit forward in our CS.
  // 
  // f_v angle_against_001 ( const ray_v & a ) const
  // {
  //   return acos ( a[2] / norm ( a ) ) ;
  // }

  // helper function to swap occupants in two vectorized objects
  // where a mask is true.

  template < typename VT , typename M >
  static void masked_swap ( VT & a , VT & b , const M & mask )
  {
    auto help = a ;
    a ( mask ) = b ;
    b ( mask ) = help ;
  }

  // this is the operator() overload for incoming ray data. The next
  // one down is for 'ninepacks'.

  void operator() (
        const std::vector < ray_v > & pv ,
        px_v & trg ,
        const std::size_t & cap = 16 ) const
  {
    // Here we have a set of sz index vectors to encode the
    // z-buffering. We'll use a 'trickle-up sort', and when
    // we're done, the first index vector should hold just the
    // same 'champions' as the single champion vector in class
    // voronoi_syn. The next one will hold next-best indices
    // where the 'maximal z' criterion came out smaller, and
    // so forth down to the 'worst' indices. We also need to
    // keep track of the z values: that's in max_z.

    std::vector < index_v > champion_v ( sz ) ;
    std::vector < f_v > max_z ( sz ) ;
    
    // next_best starts out with an invalid facet index, and if this
    // isn't updated during processing, it serves as an indicator
    // that no facet is hit by any ray - a special case which can
    // be dealt with very efficiently: produce 0000
    // 'layers' indicate how many layers we have to alpha-composit
    // to form the final result.

    int next_best = -1 ;
    int layers = 0 ;

    // we initialize max_z[0] to a very small value - no real
    // z value can ever be this small. This is so that when
    // this z value is compared to a real z, it will always
    // come out smaller.

    max_z[0] = std::numeric_limits<float>::min() ;

    // Initailly, we have no valid 'champions'

    champion_v[0] = -1 ;

    // test which of the rays in pv[0] pass through facet 0

    auto valid = env_v[0].get_mask ( pv[0] ) ;

    // do any of the rays hit the facet? then update 'next_best',
    // which, at this point, is still -1, indicating 'all misses'.
    // This can't be the case any more. We now have a first layer.
    // Also update max_z[0]: overwrite the initial 'impossible'
    // value with a true z value wherever the rays hit the facet.

    if ( any_of ( valid ) )
    {
      next_best = 0 ;
      layers = 1 ;
      champion_v[0] ( valid ) = 0 ;
      max_z[0] ( valid ) = pv[0][2] ;
    }

    // now go through the other facets in turn and perform the
    // 'trickle-up' with values for each current facet. This is
    // computationally intensive - the worst case it that we 'see'
    // the lowest angles last and they have to 'trickle up' all
    // the way. Schemes to pre-order the facets so that the most
    // likely candidate is processed first can reduce the need for
    // 'trickling up' - if the values come in in sorted order,
    // there's no need to sort them.

    for ( int i = 1 ; i < sz ; i++ )
    {
      // find the mask for facet i - if no rays hit the facet,
      // go to the next loop iteration straight away

      valid = env_v[i].get_mask ( pv[i] ) ;
      if ( none_of ( valid ) )
        continue ;

      // we let 'next_best' 'tag along.

      next_best = i ;

      // we initialize 'current_z' with an impossibly small value
      // and 'patch in' z values where rays hit the facet.
      // the z value encodes our 'quality' criterion: the
      // largest z we find is the best.

      auto & current_z ( max_z [ layers ] ) ;
      auto & current_champion ( champion_v [ layers ] ) ;

      current_z = std::numeric_limits<float>::min() ;
      current_z ( valid ) = pv[i][2] ;

      current_champion = -1 ;
      current_champion ( valid ) = i ;

      // now the trickle-up

      for ( std::size_t l = layers ; l > 0 ; --l )
      {
        // are any z values in layer l larger than in layer l-1?

        auto mask = ( max_z [ l ] > max_z [ l - 1 ] ) ;

        // if not, everything is in sorted order, we can leave
        // the loop prematurely, because all layers with lower
        // indices are already sorted.

        if ( none_of ( mask ) )
          break ;

        // if yes, swap these larger z values so that they move
        // to level l - 1, and do the same for the indices. This
        // is the 'trickle-up' - if a very large z came in
        // in some lane, the trickle-up might take it 'through'
        // all layers which were processed already, until it
        // 'hits the ceiling'.

        masked_swap ( max_z [ l ] , max_z [ l - 1 ] , mask ) ;
        masked_swap ( champion_v [ l ] , champion_v [ l - 1 ] , mask ) ;
      }

      // we need to update 'layers'

      ++layers ;
    }

    // first check: if 'layers'st is still 0, all facets were missed
    // and the result is 0000. This is the fastest outcome, and if
    // part of the output isn't covered by any facets, this will
    // save a good deal of cycles.

    trg = 0.0f ;

    if ( layers == 0 )
    {
      return ;
    }

    if ( all_of ( champion_v[0] == next_best ) )
    {
      // a very common scenario which we special-case

      px_v help ;
      env_v [ next_best ] . eval ( pv [ next_best ] , help ) ;
      const auto & alpha ( help [ nch - 1 ] ) ;
      if ( all_of ( alpha >= 1.0f ) )
      {
        // fully opaque homogeneous top layer. we're done:
        trg = help ;
        return ;
      }
    }

    // evaluate all facets. let's assume we have associated alpha.
    // We assume that the last channel in each pixel holds an alpha
    // value in the range of zero to one.
    // It's enough to evaluate the facets which have actually yielded
    // contributions, but we'd have to save and process that information.
    
    std::vector < px_v > lv ( sz ) ;
    for ( std::size_t i = 0 ; i < sz ; i++ )
    {
      env_v [ i ] . eval ( pv [ i ] , lv [ i ] ) ;
      // assert ( all_of ( lv [ i ] [ nch - 1 ] <= 1.0f ) ) ;
      // if we had unassociated alpha incoming:
      // for ( int ch = 0 ; ch < ( nch - 1 ) ; ch++ )
      //   lv [ i ] [ ch ] *= lv [ i ] [ nch - 1 ] ;
    }

    // now we have extracted the information in pv and the
    // corresponding pixel data from env_v, and we can work
    // our way through the layers to form a composite result.

    for ( std::size_t i = 0 ; i < layers ; i++ )
    {
      // at the current layer, which are the 'champions'?

      auto indexes = champion_v [ i ] ;

      // we're only interested in 'facet indices proper', not
      // in any lanes which hold -1 to indicate there are no
      // contributions in this layer.

      auto mask = ( champion_v[i] != -1 ) ;
      if ( none_of ( mask ) )
        continue ;
    
      // we convert the vector of 'champion' indices into a set
      // of indexes for gathering from the pixel values in lv.

      indexes *= ( 16 * nch ) ;
      indexes += index_v::iota() ;

      // we gather into 'help' - after the gathering is complete,
      // help will hold a composite, where each lane is set to the
      // corresponding value in that pixel vector which is the
      // champion for that lane.

      px_v help ;
      float * pf = (float*) ( & ( lv[0] ) ) ;
      for ( int ch = 0 ; ch < nch ; ch++ )
        help[ch].gather ( pf , indexes + ( ch * 16 ) ) ;

      // if we're processing the top layer, we simply copy 'help'
      // - omitting lanes which don't comntribute

      if ( i == 0 )
      {
        trg ( mask ) = help ;
      }

      // subsequent layers are alpha-composited onto the result
      // we have so far. Note that the alpha-compositing is done
      // with associated alpha, which makes it simple.

      else
      {
        for ( int ch = 0 ; ch < nch ; ch++ )
          trg[ch] ( mask ) += ( 1.0f - trg[nch-1] ) * help[ch] ;
      }

      // // if the result so far is already completely opaque, we
      // // can ignore the layers we haven't yet processed - they
      // // could not possibly contribute anything. (detrimental)
      // 
      // if ( all_of ( trg[nch-1] >= 1.0f ) )
      //   break ;
    }

    // for unassociated alpha, we'd now de-associate

    // for ( int ch = 0 ; ch < ( nch - 1 ) ; ch++ )
    // {
    //   trg[ch] /= trg[nch-1] ;
    //   trg[ch] ( trg[nch-1] == 0.0f ) = 0.0f ;
    // }
  }
    
  // operator() for 'ninepacks'. This is the same as in voronoi_syn,
  // see there for a commented version.

  void operator() ( const std::vector < ninepack_v > & pv ,
                    px_v & trg ,
                    const std::size_t & cap = 16 )
  {
    trg = 0.0f ;

    for ( const auto & cf : args.twine_spread )
    {
      for ( std::size_t facet = 0 ; facet < sz ; facet++ )
      {
        const auto & in ( pv[facet] ) ; // shorthand

        ray_v p0 { in[0] ,         in[1] ,         in[2] } ;
        ray_v du { in[3] - in[0] , in[4] - in[1] , in[5] - in[2] } ;
        ray_v dv { in[6] - in[0] , in[7] - in[1] , in[8] - in[2] } ;

        scratch [ facet ] = p0 + cf[0] * du + cf[1] * dv ;
      }

      px_v help ;
      operator() ( scratch , help , cap ) ;
      trg += cf[2] * help ;
    }  
  }
} ;

// helper class 'reproject_t' is needed for facets with translation
// parameters. It receives a ray in model space coordinates, produced
// by a stepper. The ray is projected to the reprojection plane and
// the translation parameters are subtracted, resulting in a ray
// in model space coordinates. This ray is finally rotated into the
// facet's coordinate system (by applying the basis, bs), and the
// result - a ray in the facet's coordinate system - is written to
// 'out'. The reproject_t object is concatenated with the stepper
// using a 'suffixed_t' object, and the resulting get_t produces
// the desired rays in the facet's coordinate system, which can be
// used to pick up pixel values or serve as a quality criterion.
// There are a few fine points to take into account: the tr_x etc.
// values can only be taken straight from the args object if the
// reprojection plane is not tilted (tp_* == 0) - if so, the tr_*
// values have to be rotated accordingly - this is done in the code
// setting up the reproject_t objects

typedef xel_t < float , 3 > ray_t ;

struct reproject_t 
: public unary_functor < ray_t , ray_t , 16 >
{
  typedef zimt::xel_t < zimt::xel_t < double , 3 > , 3 > basis_t ;
  typedef simdized_type < float , 16 > f_v ;
  typedef typename f_v::mask_type mask_t ;

  float tr_x , tr_y , tr_z ;
  basis_t bs ;

  reproject_t ( const ray_t & _trxyz ,
                const basis_t & _bs )
  : tr_x ( _trxyz[0] ) ,
    tr_y ( _trxyz[1] ) ,
    tr_z ( _trxyz[2] ) ,
    bs ( _bs )
  { }

  // reprojection may not be possible: if the incoming ray has
  // zero or negative z, it can't be projected onto the plane
  // one-forward. If that occurs, we return a straight-back ray,
  // which is unlikely to occur in the given context: the code
  // aims to produce an image mosaic of several images taken of
  // a flat surface with variying camera positions - typically
  // varying along the horizontal and with similar distance, and
  // also with little yaw or pitch. Incoming 'impossible' rays
  // occur only where the source images depict content which is
  // not on the plane of the flat surface (so, beyond it's
  // horizon). If the partial images don't provide content
  // 'all the way' to directly behind the camera's optical
  // axis, the result will be 0000 anyway.
  // To solve the problem in a truly general way, 'impossible'
  // rays would need to be tagged recognizably, or additional
  // information (like, a mask) would have to be passed alongside
  // the ray data - but since the reproject_t object is built into
  // a get_t object which produces rays (and no masks) this is not
  // possible with the current logic. I 'tag' the 'impossible' rays
  // with a z component of negative infinity, in case calling code
  // needs to recognize them.

  template < typename I , typename O >
  mask_t eval ( const I & _in , O & out )
  {
    auto mask = ( _in[2] <= 0.0f ) ;
    if ( all_of ( mask ) )
    {
      out[0] = 0.0f ;
      out[1] = 0.0f ;
      out[2] = - std::numeric_limits<float>::infinity() ;
      return mask_t ( false ) ;
    }
    I in ;
    in[0] = ( _in[0] / _in[2] ) - tr_x ;
    in[1] = ( _in[1] / _in[2] ) - tr_y ;
    in[2] = 1.0f - tr_z ;
    out = in[0] * bs[0] + in[1] * bs[1] + in[2] * bs[2] ;
    if ( any_of ( mask ) )
    {
      out[0] ( mask ) = 0.0f ;
      out[1] ( mask ) = 0.0f ;
      out[2] ( mask ) = - std::numeric_limits<float>::infinity() ;
    }
    return ! mask ;
  }
} ;

typedef xel_t < float , 9 > ninepack_t ;

struct reproject9_t 
: public unary_functor < ninepack_t , ninepack_t , 16 >
{
  typedef zimt::xel_t < zimt::xel_t < double , 3 > , 3 > basis_t ;

  float tr_x , tr_y , tr_z ;
  basis_t bs ;

  reproject9_t ( const ray_t & _trxyz ,
                 const basis_t & _bs )
  : tr_x ( _trxyz[0] ) ,
    tr_y ( _trxyz[1] ) ,
    tr_z ( _trxyz[2] ) ,
    bs ( _bs )
  { }

  // for twining, we have three rays to deal with: the current
  // central viewing ray and the two rays corresponding with it's
  // two 'canonical' neighbours. Here, we take the approach that
  // an 'impossible' central ray always is accompanied by two
  // other impossible rays, even if the canonical neighbours are
  // valid. The effect is that twining will go through the futile
  // exercise of evaluating the same impossible ray for each
  // twining coefficient and sum up the results, which should be
  // zero anyway. This seems wasteful, but should rarely occur,
  // due to the same resoning as above: we're trying to create
  // an image mosaic, and cases where 'impossible' rays occur
  // are contrary to the purpose.

  template < typename I , typename O >
  void eval ( const I & _in , O & out )
  {
    auto mask = ( _in[2] <= 0.0f ) ;
    if ( all_of ( mask ) )
    {
      out[0] = 0.0f ;
      out[1] = 0.0f ;
      out[2] =  - std::numeric_limits<float>::infinity() ;
      out[3] = 0.0f ;
      out[4] = 0.0f ;
      out[5] =  - std::numeric_limits<float>::infinity() ;
      out[6] = 0.0f ;
      out[7] = 0.0f ;
      out[8] =  - std::numeric_limits<float>::infinity() ;
      return ;
    }
    I in ;
    in[0] = ( _in[0] / _in[2] ) - tr_x ;
    in[1] = ( _in[1] / _in[2] ) - tr_y ;
    in[2] = 1.0f - tr_z ;
    in[3] = ( _in[3] / _in[5] ) - tr_x ;
    in[4] = ( _in[4] / _in[5] ) - tr_y ;
    in[5] = 1.0f - tr_z ;
    in[6] = ( _in[6] / _in[8] ) - tr_x ;
    in[7] = ( _in[7] / _in[8] ) - tr_y ;
    in[8] = 1.0f - tr_z ;
    auto help1 = in[0] * bs[0] + in[1] * bs[1] + in[2] * bs[2] ;
    auto help2 = in[3] * bs[0] + in[4] * bs[1] + in[5] * bs[2] ;
    auto help3 = in[6] * bs[0] + in[7] * bs[1] + in[8] * bs[2] ;
    out[0] = help1[0] ;
    out[1] = help1[1] ;
    out[2] = help1[2] ;
    out[3] = help2[0] ;
    out[4] = help2[1] ;
    out[5] = help2[2] ;
    out[6] = help3[0] ;
    out[7] = help3[1] ;
    out[8] = help3[2] ;
    if ( any_of ( mask ) )
    {
      out[0] ( mask ) = 0.0f ;
      out[1] ( mask ) = 0.0f ;
      out[2] ( mask ) =  - std::numeric_limits<float>::infinity() ;
      out[3] ( mask ) = 0.0f ;
      out[4] ( mask ) = 0.0f ;
      out[5] ( mask ) =  - std::numeric_limits<float>::infinity() ;
      out[6] ( mask ) = 0.0f ;
      out[7] ( mask ) = 0.0f ;
      out[8] ( mask ) =  - std::numeric_limits<float>::infinity() ;
    }
  }
} ;

template < int NCH ,
           template < typename , std::size_t , bool > class STP ,
           typename synopsis_t >
void fuse ( int ninputs )
{
  typedef zimt::xel_t < float , NCH > px_t ;
  typedef simdized_type < px_t , 16 > px_v ;

  typedef grok_get_t < float , 3 , 2 , 16 > gg_t ;
  typedef grok_get_t < float , 9 , 2 , 16 > gg9_t ;

  typedef environment < float , float , NCH , 16 > env_t ;

  typedef fusion_t < float , NCH , 2 , 16 , float , 3 , synopsis_t > fs_t ;
  typedef fusion_t < float , NCH , 2 , 16 , float , 9 , synopsis_t > fs9_t ;

  typedef zimt::xel_t < zimt::xel_t < double , 3 > , 3 > basis_t ;
  std::vector < basis_t > basis_v ;
  std::vector < basis_t > basis1_v ;
  std::vector < basis_t > basis2_v ;

  // we have a set of facet specs in args.facet_spec_v, which already
  // have the information about the facet images. Now we set up get_v,
  // a vector of grok_get_t containing pre-rotated steppers for the
  // facets. We also set up env_v, a vector of 'environment' objects
  // which provide access to the image data via 3D ray coordinates.
  // The pre-rotated steppers will provide 3D ray coordinates in the
  // facet images' coordinate system, the synopsis object will look
  // at the ray coordinates and evaluate some or all of the environment
  // objects at these coordinates and provide a single pixel result.
  // It's up to the synopsis object to decide how to compose the
  // result.
  // We don't automatically create an environment object for each
  // facet - special uses like caleidoscopic images may use the same
  // facet image several times. And once we implement switching to
  // a different image during a session, keeping images opened as
  // environments will stop the program from futilely re-opening an
  // image which was used before. For now, we don't use a mechanism
  // to flush images which are no longer needed, because we set up
  // the 'content map' as a static variable local to 'fuse', so when
  // the program teminates, the the environment objects held via the
  // shared_ptrs are destructed. Holding the map static keeps the
  // environments alive for image sequences as well.
  // Note that currently facets with the same filename but differing
  // hfov or projection will be re-loaded and don't share the same
  // environment object.

  static std::map < std::string , std::shared_ptr<env_t> > content_map ;
  std::vector < env_t > env_v ;

  std::vector < xel_t < double , 3 > > trxyzv ( args.nfacets ) ;

  typedef rotate_3d < double , 16 > rotate_t ;
  rotate_t r_camera ( args.roll , args.pitch , args.yaw , false ) ;

  std::vector < rotate_t > rot_plane ( args.nfacets ) ;
  std::vector < rotate_t > rot_plane_i ( args.nfacets ) ;
  std::vector < rotate_t > rot_facet ( args.nfacets ) ;

  for ( int i = 0 ; i < args.nfacets ; i++ )
  {
    auto const & fct ( args.facet_spec_v [ i ] ) ; // shorthand

    auto key =   fct.filename 
               + projection_name [ fct.projection ]
               + std::to_string ( fct.hfov )
               + "." + std::to_string ( fct.shear_g )
               + "." + std::to_string ( fct.shear_t ) ;

    auto it = content_map.find ( key ) ;
    if ( it == content_map.end() )
    {
      if ( args.verbose )
        std::cout << fct.filename << " is registered as environment "
        << key << std::endl ;
      content_map [ key ]
        = std::make_shared < env_t > ( fct ) ;
    }
    else
    {
      if ( args.verbose )
        std::cout << fct.filename << " already registered as "
        << key << " - reusing it"
        << std::endl ;
    }

    // we're now certain to have a suitable entry in the content map,
    // so we can push a copy to env_v. this is lightweight, because
    // the actual data are also held via a shared_ptr to a bspline
    // object.

    env_v.push_back ( * ( content_map [ key ] ) ) ;

    // we also need several rotations stemming from the
    // orientation of the virtual camera and the orientation
    // of the facet. We create one such 'basis' for each
    // of the facets and store them in a vector 'basis_v'

    zimt::xel_t < double , 3 > xx { 1.0 , 0.0 , 0.0 } ;
    zimt::xel_t < double , 3 > yy { 0.0 , 1.0 , 0.0 } ;
    zimt::xel_t < double , 3 > zz { 0.0 , 0.0 , 1.0 } ;

    basis_t neutral { xx , yy , zz } ;

    // all facets need to process two rotations: one due to
    // the viewer's orientation (global yaw, pitch and roll)
    // and one due to the facet's own orientation. The first
    // one is already set (it's the same throughout: rot_camera)
    // now we create and save the facet-specific one.

    rotate_3d < double , 16 > r_facet ( fct.roll ,
                                        fct.pitch ,
                                        fct.yaw ,
                                        true ) ;
    rot_facet[i] = r_facet ;

    // additionally, the reprojection plane for facets with
    // translation can be tilted. We need both the rotation and
    // it's inverse:

    rotate_3d < double , 16 > r_plane ( fct.tp_r ,
                                        fct.tp_p ,
                                        fct.tp_y ,
                                        true ) ;
    rot_plane[i] = r_plane ;

    rotate_3d < double , 16 > r_plane_i ( fct.tp_r ,
                                          fct.tp_p ,
                                          fct.tp_y ,
                                          false ) ;
    rot_plane_i[i] = r_plane_i ;

    if (    fct.tr_x == 0.0 && fct.tr_y == 0.0 && fct.tr_z == 0.0
         && fct.tp_r == 0.0 && fct.tp_p == 0.0 && fct.tp_y == 0.0 )
    {
      // there are no translation parameters. only the combined
      // rotations due to the orientation of the virtual camera
      // and of the facet's orientation are used and directly
      // fed to the stepper later on.

      xx = r_camera ( xx ) ;
      yy = r_camera ( yy ) ;
      zz = r_camera ( zz ) ;

      xx = r_facet ( xx ) ;
      yy = r_facet ( yy ) ;
      zz = r_facet ( zz ) ;

      basis_t bs { xx , yy , zz } ;
      basis_v.push_back ( bs ) ;

      // basis1_v and basis2_v are not used, but the facet's slot
      // has to be filled in with something, to keep all basis_v
      // equally sized.

      basis1_v.push_back ( neutral ) ;
      basis2_v.push_back ( neutral ) ;
    }
    else
    {
      // handle yaw and pitch of the reprojection plane. The plane
      // is always at unit distance from the master camera in the
      // origin, and the yaw and pitch both determine the point where
      // the unit sphere touches the plane and the plane's normal vector.

      // with translation parameters, the two rotations are kept
      // separate. the first brings us to the CS of the reprojection
      // plane:

      auto xx1 = r_camera ( xx ) ;
      auto yy1 = r_camera ( yy ) ;
      auto zz1 = r_camera ( zz ) ;

      xx1 = r_plane ( xx1 ) ;
      yy1 = r_plane ( yy1 ) ;
      zz1 = r_plane ( zz1 ) ;

      basis_t bs1 { xx1 , yy1 , zz1 } ;
      basis1_v.push_back ( bs1 ) ;

      // the second one brings us to the CS of the facet - rays which
      // have been through the chain of transformations will now be
      // usable with environment_t objects to produce pixel data.

      auto xx2 = r_plane_i ( xx ) ;
      auto yy2 = r_plane_i ( yy ) ;
      auto zz2 = r_plane_i ( zz ) ;

      xx2 = r_facet ( xx2 ) ;
      yy2 = r_facet ( yy2 ) ;
      zz2 = r_facet ( zz2 ) ;

      basis_t bs2 { xx2 , yy2 , zz2 } ;
      basis2_v.push_back ( bs2 ) ;

      basis_v.push_back ( neutral ) ;
    }
  }

  // now we have a set of cases to handle, depending on the type
  // of job we're running. We have to consider three points:

  // - ninputs is either 3 or 9 - 3 for direct interpolation
  //   and 9 for 'twining'

  // - we may have only a single facet, in which case we can
  //   handle the job with leaner code

  // - facets may come with or without translation parameters

  if ( ninputs == 3 )
  {
    // ninputs == 3 means that we're not using twining, and the
    // steppers produce 'ordinary' rays with three components.

    if ( args.nfacets == 1 )
    {
      // special case: there is only one facet. Using the multi-facet
      // code would work, but would produce futile overhead.
      // Note how we use a stepper which does not normalize it's result:
      // Since we don't compare the z component of the ray as quality
      // criterion (which we'd do for multiple facets) we can do without
      // the normalization.

      if ( args.verbose )
        std::cout << "using single-facet rendering" << std::endl ;

      const auto & fct ( args.facet_spec_v[0] ) ; // shorthand

      if (    fct.tr_x != 0.0 || fct.tr_y != 0.0 || fct.tr_z != 0.0
           || fct.tp_y != 0.0 || fct.tp_p != 0.0 || fct.tp_r != 0.0 )
      {
        // there are non-zero translation parameters for this facet,
        // so we have to add code to handle the reprojection of the
        // facet to the one-forward plane and the view of this plane
        // from the translated camera position. First the stepper:

        // note this setting VVV - we needn't normalize the ray here.

        STP < float , 16 , false > get_ray
          ( basis1_v[0][0] , basis1_v[0][1] , basis1_v[0][2] ,
            args.width , args.height ,
            args.x0 , args.x1 , args.y0 , args.y1 ) ;

        // the stepper, which only uses the rotation of the virtual
        // camera (basis1) is suffixed with a functor handling the
        // projection of the facet to the reprojection plane and the
        // application of the camera translation.

        ray_t trxyz { fct.tr_x , fct.tr_y , fct.tr_z } ;
        trxyz = rot_plane[0] ( trxyz ) ;

        reproject_t rprj ( trxyz , basis2_v[0] ) ;

        suffixed_t < float , 3 , 2 , 16 > sfg ( get_ray , rprj ) ;

        // the resulting object is a get_t in it's own right, and
        // all that's left to do now is to call work:

        work ( sfg , env_v[0] ) ;
      }
      else
      {
        // all translation parameters are zero. this is easy:

        STP < float , 16 , false > get_ray
          ( basis_v[0][0] , basis_v[0][1] , basis_v[0][2] ,
            args.width , args.height ,
            args.x0 , args.x1 , args.y0 , args.y1 ) ;

        work ( get_ray , env_v[0] ) ;
      }
    }
    else
    {
      // there are several facets.

      // we'll have several get_t objects, one for each facet,
      // producing rays. These get_t objects are 'grokked' to erase
      // ther specific type (STP) and stored in this vector:

      std::vector < gg_t > get_v ;

      for ( int i = 0 ; i < args.nfacets ; i++ )
      {
        const auto & fct ( args.facet_spec_v[i] ) ; // shorthand

        if (    fct.tr_x != 0.0 || fct.tr_y != 0.0 || fct.tr_z != 0.0
             || fct.tp_y != 0.0 || fct.tp_p != 0.0 || fct.tp_r != 0.0 )
        {
          // if there are translation parameters for this facet,
          // proceed as in the single-facet case: suffix the stepper
          // with a reproject_t and then push the resulting object
          // to get_v (it's 'grokked' in the process, because the
          // get_v vector is a vector of grok_get_t)

          STP < float , 16 , false > get_ray
            ( basis1_v[i][0] , basis1_v[i][1] , basis1_v[i][2] ,
              args.width , args.height ,
              args.x0 , args.x1 , args.y0 , args.y1 ) ;

          ray_t trxyz { fct.tr_x , fct.tr_y , fct.tr_z } ;
          trxyz = rot_plane[i] ( trxyz ) ;

          reproject_t rprj ( trxyz , basis2_v[i] ) ;

          suffixed_t < float , 3 , 2 , 16 > sfg ( get_ray , rprj ) ;
          get_v.push_back ( sfg ) ;
        }
        else
        {
          // no translation parameters.

          // set up a simple single-coordinate stepper of the type
          // fixed by 'STP'. This route is taken with direct b-spline
          // interpolation (itp 1)

          // note this setting VV - we need to normalize the ray here.

          STP < float , 16 , true > get_ray
            ( basis_v[i][0] , basis_v[i][1] , basis_v[i][2] ,
              args.width , args.height ,
              args.x0 , args.x1 , args.y0 , args.y1 ) ;

          get_v.push_back ( get_ray ) ;
        }
      }

      // for now, we use a hard-coded synopsis-forming object

      synopsis_t vs ( env_v ) ;

      // now we can create the fusion_t object

      fs_t fs ( get_v , vs ) ;

      // now we call the final 'work' template which uses the fusion_t
      // which directly produces a pixel value, so we use a pass_through
      // act functor, rather than having the interpolator step in the
      // act functor.

      zimt::pass_through < float , NCH , 16 > act ;

      work ( fs , act ) ;
    }
  }
  else // ninputs == 9
  {
    // this is the code we use for twining. instead of processing
    // simple 3D rays, we process 'ninepacks'.

    if ( args.nfacets == 1 )
    {
      const auto & fct ( args.facet_spec_v[0] ) ;

      // special case: there is only one facet. We have to wrap the
      // single evaluator in a twine_t object to process the 'ninepacks'
      // which the deriv_stepper in get_ray will produce. Using the
      // multi-facet code would work, but would produce futile overhead.
      // Note how we use a stepper which does not normalize it's result:
      // Since we don't compare the z component of the ray as quality
      // criterion (which we'd do for multiple facets) we can do without
      // the normalization.

      if ( args.verbose )
        std::cout << "using single-facet rendering" << std::endl ;

      if (    fct.tr_x != 0.0 || fct.tr_y != 0.0 || fct.tr_z != 0.0
           || fct.tp_y != 0.0 || fct.tp_p != 0.0 || fct.tp_r != 0.0 )
      {
        deriv_stepper < float , 16 , STP , true > get_ray
            ( basis1_v[0][0] , basis1_v[0][1] , basis1_v[0][2] ,
              args.width , args.height ,
              args.x0 , args.x1 , args.y0 , args.y1 ) ;

        twine_t < NCH , 16 > twenv ( env_v[0] , args.twine_spread ) ;
        ray_t trxyz { fct.tr_x , fct.tr_y , fct.tr_z } ;
        trxyz = rot_plane[0] ( trxyz ) ;
        reproject9_t rprj ( trxyz , basis2_v[0] ) ;
        suffixed_t < float , 9 , 2 , 16 > sfg ( get_ray , rprj ) ;
        work ( sfg , twenv ) ;
      }
      else
      {
        // for single-facet operation with twining, we also use steppers
        // with normalized rays:           VVVV

        deriv_stepper < float , 16 , STP , true > get_ray
            ( basis_v[0][0] , basis_v[0][1] , basis_v[0][2] ,
              args.width , args.height ,
              args.x0 , args.x1 , args.y0 , args.y1 ) ;

        twine_t < NCH , 16 > twenv ( env_v[0] , args.twine_spread ) ;

        work ( get_ray , twenv ) ;
      }
    }
    else
    {
      // again we set up a vector of grok_get_t, but this time they
      // are getters yielding 'ninepack' values

      std::vector < gg9_t > get_v ;

      for ( int i = 0 ; i < args.nfacets ; i++ )
      {
        const auto & fct ( args.facet_spec_v[i] ) ;

        if (    fct.tr_x != 0.0 || fct.tr_y != 0.0 || fct.tr_z != 0.0
             || fct.tp_y != 0.0 || fct.tp_p != 0.0 || fct.tp_r != 0.0 )
        {
          // we set up a deriv_stepper (specialized with the stepper
          // type given in template argument STP) which produces ninepacks

          // note this setting               VVV - we need normalized rays.

          deriv_stepper < float , 16 , STP , true > get_ray
              ( basis1_v[i][0] , basis1_v[i][1] , basis1_v[i][2] ,
                args.width , args.height ,
                args.x0 , args.x1 , args.y0 , args.y1 ) ;

          // and push it to the get_v vector

          ray_t trxyz { fct.tr_x , fct.tr_y , fct.tr_z } ;
          trxyz = rot_plane[i] ( trxyz ) ;

          reproject9_t rprj ( trxyz , basis2_v[i] ) ;

          suffixed_t < float , 9 , 2 , 16 > sfg ( get_ray , rprj ) ;

          get_v.push_back ( sfg ) ;
        }
        else
        {
          deriv_stepper < float , 16 , STP , true > get_ray
              ( basis_v[i][0] , basis_v[i][1] , basis_v[i][2] ,
                args.width , args.height ,
                args.x0 , args.x1 , args.y0 , args.y1 ) ;

          // and push it to the get_v vector

          get_v.push_back ( get_ray ) ;
        }
      }

      // we set up the synopsis-forming object and introduce it to
      // the fusion_t object - which is now used via it's ninepack-
      // processing member function.

      synopsis_t vs ( env_v ) ;

      fs9_t fs ( get_v , vs ) ;

      // now we call the final 'work' template which uses the fusion_t
      // which directly produces a pixel value, so we use a pass_through
      // act functor.

      zimt::pass_through < float , NCH , 16 > act ;

      work ( fs , act ) ;
    }
  }
}

// to call the appropriate instantiation of 'fuse' (above)
// we need to do some dispatching: picking types depending
// on run-time variables. We achieve this with a staged 'roll_out'
// routine: every stage processes one argument and routes to
// specialized code. The least specialized roll_out variant is
// the lowest one down, this here is the final stage where we have
// the number of channels and the stepper as template arguments.
// Finally, we dispatch to fuse, with the synopsis-forming
// object's type depending on whether we have data with alpha
// channel or not. 'fuse' special-cases if there is only a
// single facet in play.

template < int NCH ,
           template < typename , std::size_t > class STP >
void roll_out ( int ninputs )
{
  typedef environment < float , float , NCH , 16 > env_t ;

  // quick shot: assume one- and three-channel images are
  // without alpha channel, two- and four-channel images
  // are with alpha channel.

  if constexpr ( NCH == 1 || NCH == 3 )
    fuse < NCH , STP , voronoi_syn < env_t > > ( ninputs ) ;
  else
    fuse < NCH , STP , voronoi_syn_plus < env_t > > ( ninputs ) ;
}

// while evolving the code, I narrow the scope to three channels
// only and only three projections. This lowers turn-around time
// considerably.

#define NARROW_SCOPE

// we have the number of channels as a template argument from the
// roll_out below, now we roll_out on the projection and instantiate
// the next roll_out level with a stepper type which fits the
// projection. Note that a stepper steps through the *target* image;
// it's output are ray coordinates in the source image's coordinate
// system. Also note that roll_out delegates to 'fuse' if the input
// consists of facets.

template < int NCH >
void roll_out ( int ninputs ,
                projection_t projection )
{
  switch ( projection )
  {
    case SPHERICAL:
      roll_out < NCH , spherical_stepper > ( ninputs ) ;
      break ;
#ifndef NARROW_SCOPE
    case CYLINDRICAL:
      roll_out < NCH , cylindrical_stepper > ( ninputs ) ;
      break ;
#endif
    case RECTILINEAR:
      roll_out < NCH , rectilinear_stepper > ( ninputs ) ;
      break ;
    case FISHEYE:
      roll_out < NCH , fisheye_stepper > ( ninputs ) ;
      break ;
#ifndef NARROW_SCOPE
    case STEREOGRAPHIC:
      roll_out < NCH , stereographic_stepper > ( ninputs ) ;
      break ;
    case CUBEMAP:
      roll_out < NCH , cubemap_stepper > ( ninputs ) ;
      break ;
    case BIATAN6:
      roll_out < NCH , biatan6_stepper > ( ninputs ) ;
      break ;
#endif
    default:
      std::cerr << "unhandled projection # " << int(projection)
                << std::endl ;
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
    #ifdef HWY_TARGET_STR
      hwy_target = HWY_TARGET ;
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
#ifndef NARROW_SCOPE
      case 1:
        roll_out < 1 > ( ninputs , projection ) ;
        break ;
      case 2:
        roll_out < 2 > ( ninputs , projection ) ;
        break ;
#endif
      case 3:
        roll_out < 3 > ( ninputs , projection ) ;
        break ;
// #ifndef NARROW_SCOPE
      case 4:
        roll_out < 4 > ( ninputs , projection ) ;
        break ;
// #endif
      std::cerr << "unhandled channel count " << nchannels
                << std::endl ;
    }
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
