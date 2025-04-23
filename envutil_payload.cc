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
#include "lens_correction.h"

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

  template < typename U >
  void as_matrix ( r3_t < U > & m )
  {
    xel_t < double , 3 > ex { 1.0 , 0.0 , 0.0 } ;
    xel_t < double , 3 > ey { 0.0 , 1.0 , 0.0 } ;
    xel_t < double , 3 > ez { 0.0 , 0.0 , 1.0 } ;
    eval ( ex , ex ) ;
    eval ( ey , ey ) ;
    eval ( ez , ez ) ;
    m[0] = ex ;
    m[1] = ey ;
    m[2] = ez ;
  }
} ;

// directly produce an r3_t from Euler angles. This might also
// be done directly with a very large formula - this is easier:

template < typename T >
r3_t < T > make_r3_t ( T r , T p , T y ,
                       bool inverse = false )
{
  rotate_3d < T , 1 > rq3 ( r , p , y , inverse ) ;
  r3_t < T > r3m ;
  rq3.as_matrix ( r3m ) ;
  return r3m ;
}

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
  // reallocate the memory. If there is a cropping specification
  // for the output image (from a PTO file's p-line), the target
  // array is sized accordingly.

  std::size_t w , h ;
  if ( args.store_cropped )
  {
    w = args.p_crop_x1 - args.p_crop_x0 ;
    h = args.p_crop_y1 - args.p_crop_y0 ;
  }
  else
  {
    w = args.width ;
    h = args.height ;
  }

  zimt::array_t < 2 , px_t > trg ( { w , h } ) ;
  
  // set up a zimt::storer to populate the target array with
  // zimt::process. This is the third component needed for
  // zimt::process - we already have the get_t and act.
  
  zimt::storer < float , nchannels , 2 , LANES > cstor ( trg ) ;
  
  // use the get, act and put components with zimt::process
  // to produce the target images. This is the point where all
  // the state we have built up is finally put to use, running
  // a multithreaded pipeline which fills the target image.
  
  zimt::bill_t bill ;

  // If there is a cropping specification for the output image
  // (from a PTO file's p-line), the discrete coordinates fed into
  // the pixel pipeline have to be raised appropriately:

  if ( args.store_cropped )
  {
    bill.get_offset.push_back ( args.p_crop_x0 ) ;
    bill.get_offset.push_back ( args.p_crop_y0 ) ;
  }

  std::chrono::system_clock::time_point start
    = std::chrono::system_clock::now() ;

  // bill.njobs = 1 ; // to run the rendering single-threaded

  float unbrighten = 1.0f ;
  if ( args.single >= 0 )
  {
    const auto & fct ( args.facet_spec_v [ args.single ] ) ;
    if ( fct.brighten != 1.0 )
    {
      unbrighten = 1.0 / fct.brighten ;
    }
  }

  if ( unbrighten != 1.0f )
  {
    // for 'single' jobs, undo the target facet's 'brighten'
    // TODO tentative. linear vs. sRGB is still unresolved!
    // Where 'brighten' is used now - so, here, and 'on the
    // other end of the pipeline' - a whole complex of processing
    // needs to be 'slotted in' to handle e.g. colour space
    // conversions.

    px_t amplify = unbrighten ;
    if ( nchannels == 2 || nchannels == 4 )
      amplify [ nchannels - 1 ] = 1.0 ;
    amplify_type < px_t , px_t , px_t , LANES > amp ( amplify ) ;
    zimt::process ( trg.shape , get , act + amp , cstor , bill ) ;
  }
  else
  {
    zimt::process ( trg.shape , get , act , cstor , bill ) ;
  }
  
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

  // if ( args.seqfile != std::string() )
  // {
  //   push_video_frame ( trg ) ;
  // }
  // else
  {
    if ( args.verbose )
      std::cout << "saving output image: " << args.output << std::endl ;

    save_array < nchannels > ( args.output , trg ) ;
  }
}

// we'll use a common wrapper class for all synopsis-forming objects.
// This class wraps a specific synopsis-forming object handling only
// 'normal' ray coordinates as input and adds an operator() overload
// which handles 'ninepack' arguments and twining, using operator()
// of the 'wrappee' to evaluate individual rays in the twining spread.

template < std::size_t NCH , template < std::size_t > class SF >
struct synopsis_t
: public SF < NCH >
{
  typedef SF < NCH > base_t ;
  using base_t::sz ;
  using base_t::operator() ;
  using typename base_t::ray_v ;
  using typename base_t::px_v ;
  using typename base_t::env_t ;

  typedef zimt::xel_t < float , 9 > ninepack_t ;
  typedef simdized_type < ninepack_t , LANES > ninepack_v ;

  // 'scratch' will hold all ray packets which will react with
  // a given twining coefficient 'cf'
  
  std::vector < ray_v > scratch ;

  // we hold a copy of the spread here, which is modified in the c'tor
  // to incorporate the 'bias' values.

  std::vector < zimt::xel_t < float , 3 > > spread ;

  synopsis_t ( std::vector < env_t > & env_v ,
               const float & bias = 4.0f )
  : base_t ( env_v ) ,
    scratch ( env_v.size() ) ,
    spread ( args.twine_spread )
  {
    // we apply the 'bias' value to the twining coefficients, to avoid
    // the multiplication with the dx/dy value which would be needed
    // otherwise. This 'bias' is the reciprocal value of the 'bias_x'
    // and 'bias_y' values used in steppers to form slightly
    // offsetted planar coordinates for twining.

    for ( auto & cf : spread )
    {
      cf[0] *= bias ;
      cf[1] *= bias ;
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
                    const std::size_t & cap = LANES )
  {
    typedef simdized_type < float , LANES > f_v ;
    typedef simdized_type < zimt::xel_t < float , 3 > , LANES > ray_v ;
    typedef simdized_type < zimt::xel_t < float , 9 > , LANES > ray9_v;

    // trg is built up as a weighted sum of contributing values
    // resulting from the reaction with a given coefficient, so we
    // clear it befor we begin.
  
    trg = 0.0f ;

    // for each coefficient in the twining 'spread'

    for ( const auto & cf : spread )
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
      base_t::operator() ( scratch , help , cap ) ;

      // 'help' is the pixel value, which is weighted with the
      // weighting value in the twining kernel and added to what's
      // in trg already.

      trg += cf[2] * help ;
    }  
    // and that's it - trg has the result.
  }
} ;

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

template < std::size_t NCH >
struct _voronoi_syn
{
  typedef environment < float , float , NCH , LANES > env_t ;

  typedef typename env_t::px_t px_t ;
  typedef simdized_type < px_t , LANES > px_v ;
  static const std::size_t nch = px_t::size() ;
  const int sz ;
  typedef zimt::xel_t < float , 3 > ray_t ;
  typedef simdized_type < float , LANES > f_v ;
  typedef simdized_type < ray_t , LANES > ray_v ;

  std::vector < env_t > & env_v ;

  _voronoi_syn ( std::vector < env_t > & _env_v )
  : env_v ( _env_v ) ,
    sz ( _env_v.size() )
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

  // this is the operator() overload for incoming ray data.

  void operator() (
        const std::vector < ray_v > & pv ,
        px_v & trg ,
        const std::size_t & cap = LANES ) const
  {
    // this vector will contain the index of the 'most suitable'
    // facet at the given location, according to the criterion.

    simdized_type < int , LANES > champion_v ( -1 ) ;

    // we initialize 'max_z' to a very small value - no
    // 'real' z value can ever be this small.

    simdized_type < float , LANES >
      max_z ( std::numeric_limits<float>::lowest() ) ;

    // next_best starts out with an invalid facet index, and if this
    // isn't updated during processing, it serves as an indicator
    // that no facet is hit by any ray - a special case which can
    // be dealt with very efficiently: produce 0000

    int next_best = -1 ;

    // test which of the rays in pv[0] pass through facet 0

    auto valid = env_v[0].get_mask ( pv[0] ) ;

    // do any of the rays hit the facet? then update 'next_best',
    // which, at this point, is still -1, indicating 'all misses'.
    // This can't be the case any more. Where rays hit this first
    // facets (indicated by 'valid'), we now have a champion, which
    // we save in champion_v. We also update max_z: we overwrite
    // the initial 'impossibly small' value with a real z value
    // where the rays hit the facet.

    if ( any_of ( valid ) )
    {
      next_best = 0 ;
      champion_v ( valid ) = 0 ;
      max_z ( valid ) = pv[0][2] * env_v[0].recip_step ;
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

      simdized_type < float , LANES >
        current_z ( std::numeric_limits<float>::lowest() ) ;

      current_z ( valid ) = pv[i][2] * env_v[i].recip_step ;

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
      // Luckily, this is also usually quite rare - it only occurs
      // where facets collide or near facet edges.

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
} ;

template < std::size_t NCH >
using voronoi_syn = synopsis_t < NCH , _voronoi_syn > ;

// variant with z-buffering for facets with alpha channel:

template < std::size_t NCH , bool fix_duv = true >
struct _voronoi_syn_plus
{
  typedef environment < float , float , NCH , LANES > env_t ;

  typedef typename env_t::px_t px_t ;
  typedef simdized_type < px_t , LANES > px_v ;
  static const std::size_t nch = px_t::size() ;
  const int sz ;
  typedef zimt::xel_t < float , 3 > ray_t ;
  typedef simdized_type < float , LANES > f_v ;
  typedef simdized_type < ray_t , LANES > ray_v ;
  typedef simdized_type < int , LANES > index_v ;

  std::vector < env_t > & env_v ;

  _voronoi_syn_plus ( std::vector < env_t > & _env_v )
  : env_v ( _env_v ) ,
    sz ( _env_v.size() )
  { }

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
        const std::size_t & cap = LANES ) const
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

    max_z[0] = std::numeric_limits<float>::lowest() ;

    // Initially, we have no valid 'champions'

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
      max_z[0] ( valid ) = pv[0][2] * env_v[0].recip_step ;
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

      current_z = std::numeric_limits<float>::lowest() ;
      current_z ( valid ) = pv[i][2] * env_v[i].recip_step ;

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

      indexes *= int ( LANES * nch ) ;
      indexes += index_v::iota() ;

      // we gather into 'help' - after the gathering is complete,
      // help will hold a composite, where each lane is set to the
      // corresponding value in that pixel vector which is the
      // champion for that lane.

      px_v help ;
      float * pf = (float*) ( & ( lv[0] ) ) ;
      for ( int ch = 0 ; ch < nch ; ch++ )
        help[ch].gather ( pf , indexes + ( ch * LANES ) ) ;

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
} ;

template < std::size_t NCH >
using voronoi_syn_plus = synopsis_t < NCH , _voronoi_syn_plus > ;

// class hdr_merge_syn codes synopsis-forming objects which produce
// HDR-merged output of the facets in facet_spec_v. It expects input
// in scene-referred linear form, and if the data have transparency,
// it expects associated alpha (a.k.a premultiplied). There are also
// some implicit expectations about the content: it's assumed that
// the Eev values given for the facets are valid, and that the range
// of values found in the facets (as they appear after OIIO has provided
// them) are in the range of [0,1]. Normal facet processing looks at
// the Eev value and the desired Eev for the output (or, if that is
// missing the average Eev of all facets) and puts a 'brightening'
// factor into the processing chain, which 'harmonizes' facet
// brightness so that ideally facets showing identical content should
// hold matching intensity values. This pre-processing already presents
// the data as they are needed for HDR merging: here, the input consists
// of a 'bracket' of several exposures done with different exposure
// times (and, therefore, different Eev values), but apart from the
// varying brightness, the images should match very well - which may
// need a 'registration' step with tools like cpfind or align_image_stack
// from the hugin toolchain. There are two points in using an exposure
// bracket: short exposures will capture bright content which produces
// overexposure in longer exposures. longer exposures provide better
// quality data for dark content, which is more prone to noise and
// quantization artifacts. hdr_merge_syn takes both these points into
// account and composes a target image with extended dynamic range,
// picking only pixels from the image set which are not overexposed
// and giving more weight to the longer exposures in dark areas.
// The target image has it's grey point set to match with the target
// Eev (prescribed in the p-line of a PTO or the average of the input),
// so if, for example, an ordinary -2/0/+2 bracket is used as input,
// the target image will have the same overall brightness as the middle
// exposure - apart from very bright areas, where the middle exposure
// is already overexposed and data from the long exposure 'take over'.
// The target image - being HDR - will no longer 'fit' into the 'normal'
// dynamic range. You want to save it in a format which can handle an
// extended dynamic range, like OpenEXR. As of this writing, envutil
// does not do any colour space alterations, so you want an output
// format which can contain the scene-referred linear data, and OpenEXR
// is ideal for the purpose. What's good input? envutil can process RAW
// images, which is a good choice, but obtaining a PTO with registration
// pertaining to RAW images is tricky -
// (see my attempts here:https://bitbucket.org/kfj/oiio_fileio)
// another good option is to convert the images to linear TIFF with a
// raw converter like dcraw, which uses the same processing as OIIO's
// libraw plugin. The resulting TIFFs will have the same geometry as
// the RAW files, and if you have a PTO referring to these TIFFs, a
// handy trick is to do a search/replace on such a TIFF, replacing
// the TIFF extension by the original RAW file's extension. Such a
// 'dubbed' PTO file can be used with lux or envutil, which use OIIO
// for image IO, and can directly access the RAW file's content using
// OIIO's libraw plugin.

// Now for the technicalities. I mentioned that envutil does exclude
// overexposed pixels. envutil is strict and 'bans' pixels even if only
// one or two colour channels are overexposed. So it uses a stage of
// 'grey projection' where the *maximum of all colour channels* is taken
// as grey value, and this grey value is used to calculate the 'quality'
// criterion for each pixel. Another popular way to do the 'grey
// projection' is to use the average of the colour channels, but this
// can result in some overexposed content 'slipping into' the result,
// which, I think, is best avoided.
// The 'quality' is provided by comparing the grey value to a reference
// value - the 'optimum' value. This is the grey value right in the
// middle of the range of values found in the image data. Since we're
// working with normalized intensity values and pre-brightening, the
// 'optimum' is 0.5, times the brightening factor. The quality will be
// 'high' for grey values near the 'optimum' and lower the further away
// from the optimum it is. The 'quality' is also normalized, so that
// grey values coinciding with the optimum receive quality 1.0, and
// grey values of zero receive quality zero, and overexposed pixels
// receive a very small quality value which is just enough to produce
// a white output in places where there is only one facet and it's
// content is overexposed - if there are any facets providing better
// data for these areas, they will 'win' by a wide margin.
// The last step is to form a weighted sum of the pixels, weighting each
// with it's corresponding 'quality' value. I have outlined how the
// 'quality' value is calculated, and implied that all pixels which are
// not overexposed are qualified to participate in the output. This is
// debatable - one might exclude dark pixels from short exposures in
// favour of content from longer exposures, which should have better
// intensity resolution and less noise. There's room here for a lot of
// fine-tuning, and this is subject to reasearch, just as the grey
// projection step can be done differently, and the translation of
// 'distance to optimum' to a quality value can be done in different
// ways as well - common schemes use a simple triangular function
// (currently envutil does so) or a bell-shaped curve using a negative
// exponential of the distance-to-optimum.

template < std::size_t NCH >
struct _hdr_merge_syn
{
  typedef environment < float , float , NCH , LANES > env_t ;

  typedef typename env_t::px_t px_t ;
  typedef simdized_type < px_t , LANES > px_v ;
  static const std::size_t nch = px_t::size() ;
  const int sz ;
  typedef zimt::xel_t < float , 3 > ray_t ;
  typedef simdized_type < float , LANES > f_v ;
  typedef simdized_type < ray_t , LANES > ray_v ;

  std::vector < env_t > & env_v ;
  std::vector < float > optimum ;

  _hdr_merge_syn ( std::vector < env_t > & _env_v )
  : env_v ( _env_v ) ,
    sz ( _env_v.size() ) ,
    optimum ( _env_v.size() )
  {
    float lowest_brighten = 100000.0f ;
    float highest_brighten = -1.0f ;
    int lowest_brighten_fct ;
    int highest_brighten_fct ;

    for ( int i = 0 ; i < args.nfacets ; i++ )
    {
      const auto & fct ( args.facet_spec_v [ i ] ) ;
      if ( fct.brighten < lowest_brighten )
      {
        lowest_brighten = fct.brighten ;
        lowest_brighten_fct = i ;
      }
      if ( fct.brighten > highest_brighten )
      {
        highest_brighten = fct.brighten ;
        highest_brighten_fct = i ;
      }
    }
  }

  // we have a two-stage process for obtaining a 'quality' value.
  // These two functions calculate a 'quality' value for a grey value,
  // either without or with alpha present. Below is the function
  // processing pixel data - there, the grey value is formed and,
  // if present, the alpha value extracted, then the grey-processing
  // get_quality overload is called to produce the result.
  // TODO: the simple triangular function I use doesn't work so
  // well for the dark areas - I think content from long exposures
  // must be suppressed more for the dark areas, so I divide by
  // optimum squared. the extra division by optimum boosts long
  // exposures (they have low optimum).
  
  f_v get_quality ( f_v _m , float optimum ) const
  {
    f_v grey = _m ;
    grey -= optimum ;
    grey = optimum - abs ( grey ) ;
    grey /= ( optimum * optimum ) ;

    // incoming _m is either the single grey channel or the maximum
    // of R, G and B. The three operations above have yielded a
    // quality value which is 1.0 for 'optimum' _m and near zero
    // for _m at either end of the range. We want to produce a small
    // quality value only if we're at the upper end (_m overexposed)
    // to make the result come out at full intensity, hence the final
    // multiplication with _m. If _m is at the lower end (black pixel)
    // the final multiplication will result in zero quality.

    grey ( grey < 0.00001f ) = 0.00001f * _m ;
    return grey ;
  }

  f_v get_quality ( f_v _grey , f_v alpha , float optimum ) const
  {
    f_v grey = _grey ;
    auto ao = alpha * optimum ;
    grey -= ao ;
    grey = ao - abs ( grey ) ;
    grey /= ( ao * ao ) ;
    grey ( ! ( grey > 0.00001f ) ) = 0.00001f * _grey ;
    return grey ;
  }

  f_v get_quality ( const px_v & px , float optimum ) const
  {

    if constexpr ( px_v::size() == 1 )
    {
      f_v q = get_quality ( px[0] , optimum ) ;
      return q ;
    }
    else if constexpr ( px_v::size() == 3 )
    {
      // grey projection picks the maximum of R, G and B:
      f_v grey = max ( px[0] , px[1] ) ;
      grey = max ( grey , px[2] ) ;
      f_v q = get_quality ( grey , optimum ) ;
      return q ;
    }
    else if constexpr ( px_v::size() == 2 )
    {
      // alpha is in the last channel
      f_v alpha = px[1] ;
      f_v q = get_quality ( px[0] , alpha , optimum ) ;
      return q ;
    }
    else if constexpr ( px_v::size() == 4 )
    {
      // alpha is in the last channel
      f_v alpha = px[3] ;
      // grey projection picks the maximum of R, G and B:
      f_v grey = max ( px[0] , px[1] ) ;
      grey = max ( grey , px[2] ) ;
      f_v q = get_quality ( grey , alpha , optimum ) ;
      return q ;
    }
  }

  // this is the operator() overload for incoming ray data. The next
  // one down is for 'ninepacks'.

  void operator() (
        const std::vector < ray_v > & pv ,
        px_v & trg ,
        const std::size_t & cap = LANES ) const
  {
    // go through the facets in turn and evaluate the env, calculate
    // the quality criterion, form a weighted sum of the input pixels

    trg = 0.0f ;
    f_v qsum ( 0.0f ) ;

    for ( int i = 0 ; i < sz ; i++ )
    {
      // evaluate the environment object in env_v[i]

      px_v px ;
      env_v [ i ] . eval ( pv [ i ] , px ) ;

      // calculate the 'quality' value for the pixel data

      f_v quality = get_quality ( px , optimum[i] ) ;

      // add it to qsum, to normalize in the final step

      qsum += quality ;

      // we use 'if constexpr' to code for the variants in one
      // function:

      if constexpr ( px_v::size() == 1 || px_v::size() == 3 )
      {
        // no alpha present. this is simple:

        trg += ( px * quality ) ;
      }
      else if constexpr ( px_v::size() == 2 )
      {
        // two channels, interpreted as grey + alpha

        // deassociate. we're assuming alpha is associated!

        f_v grey = 0.0f ;
        grey ( px[1] > 0.000001f ) = px[0] / px[1] ;

        // apply weight to grey channel and sum up

        trg[0] += ( grey * quality ) ;

        // result alpha is max of alpha of all pixels

        trg[1] = max ( trg[1] , px[1] ) ;
      }
      else if constexpr ( px_v::size() == 4 )
      {
        // deassociate

        xel_t < f_v , 3 > rgb = f_v ( 0.0f ) ;

        rgb[0] ( px[3] > 0.000001f ) = px[0] / px[3] ;
        rgb[1] ( px[3] > 0.000001f ) = px[1] / px[3] ;
        rgb[2] ( px[3] > 0.000001f ) = px[2] / px[3] ;

        // apply weight to colour channels and sum up

        trg[0] += ( rgb[0] * quality ) ;
        trg[1] += ( rgb[1] * quality ) ;
        trg[2] += ( rgb[2] * quality ) ;

        // result alpha is max of alpha of all pixels

        trg[3] = max ( trg[3] , px[3] ) ;
      }
    }

    // final stage: normalization and re-association of alpha. Note
    // that the output will be set to zero where the sum of quality
    // values for all facets doesn't exceed zero. The formulation
    // with the negated mask is to mask out nan pixels as well.

    if constexpr ( px_v::size() == 1 || px_v::size() == 3 )
    {
      trg /= qsum ;      // normalize
      trg ( ! ( qsum > 0.0f ) ) = 0.0f ; // suppress bad q
    }
    else if constexpr ( px_v::size() == 2 )
    {
      trg[0] /= qsum ;   // normalize
      trg[0] ( ! ( qsum > 0.0f ) ) = 0.0f ; // suppress bad q
      trg[0] *= trg[1] ; // reassociate
    }
    else if constexpr ( px_v::size() == 4 )
    {
      // normalize

      trg[0] /= qsum ;
      trg[1] /= qsum ;
      trg[2] /= qsum ;

      // suppress low q values

      trg[0] ( ! ( qsum > 0.0f ) ) = 0.0f ;
      trg[1] ( ! ( qsum > 0.0f ) ) = 0.0f ;
      trg[2] ( ! ( qsum > 0.0f ) ) = 0.0f ;

      // reassociate alpha

      trg[0] *= trg[3] ;
      trg[1] *= trg[3] ;
      trg[2] *= trg[3] ;
    }
  }
} ;

template < std::size_t NCH >
using hdr_merge_syn = synopsis_t < NCH , _hdr_merge_syn > ;

template < typename T , std::size_t L >
struct generic_r3
: public unary_functor < xel_t<T,3> , xel_t<T,3> , L >
{
  grok_type < xel_t<T,3> , xel_t<T,3> , L > ev ;
  typedef r3_t < T > r_t ;

  generic_r3 ( const facet_base & ft ,
               const facet_base & fs )
  {
    // here, the facet passed in ft is in the 'camera' position
    // r_camera takes us from target coordinates to model space:

    r_t r_camera = make_r3_t ( ft.roll , ft.pitch , ft.yaw , false ) ;

    // this rotation takes us to the target's translation plane

    r_t rt_tp = make_r3_t ( ft.tp_r , ft.tp_p , ft.tp_y , true ) ;

    // this rotation takes us back to model space

    r_t rt_tpi = make_r3_t ( ft.tp_r , ft.tp_p , ft.tp_y , false ) ;

    // this rotation takes us from model space to the source facet's
    // translation plane

    r_t rs_tp = make_r3_t ( fs.tp_r , fs.tp_p , fs.tp_y , true ) ;

    // this rotation takes us back to model space

    r_t rs_tpi = make_r3_t ( fs.tp_r , fs.tp_p , fs.tp_y , false ) ;

    // this rotation takes us to the source facet's CS

    r_t r_facet = make_r3_t ( fs.roll , fs.pitch , fs.yaw , true ) ;

    bool have_ttp = ( ft.tr_x != 0 || ft.tr_y != 0 || ft.tr_z != 0 ) ;
    bool have_stp = ( fs.tr_x != 0 || fs.tr_y != 0 || fs.tr_z != 0 ) ;

    // the PTO 'translation' is implemented as a rotation to the
    // reprojection plane, the reprojection itself by 'casting the
    // rays down onto the plane' (division by z), a shift which moves
    // to the CS of the 'virtual camera' and another rotation back
    // to the next CS - e.g. model space. The position of the virtual
    // camera is given in model space coordinates, which I find a
    // debatable choice - I think it may make more sense to define
    // this position in the CS of to the reprojection plane, which
    // is how I've coded it here. Maybe using model space coordinates
    // works better when one uses plane definitions and rays in
    // model space, rather than doing the CS transformations I use
    // here.

    xel_t<T,3> shift_t { ft.tr_x , ft.tr_y , ft.tr_z } ;

    // the translation is given in the model space CS, but we want
    // to apply it in the CS of the translation plane. So we have
    // to rotate the shift accordingly.

    if ( ft.tp_y != 0 || ft.tp_p != 0 || ft.tp_r != 0 )
      shift_t = rotate ( xel_t<double,3> ( shift_t ) , rt_tp ) ;

    // dcp is a scaling factor, which is needed for re-cration of
    // single images. There, the reprojection plane, onto which the
    // rays are cast in the first step, is not at unit distance from
    // the origin, which, at that point, is where the 'virtual camera'
    // of the translation is located - instead it's closer or further
    // away, depending on the distance from the translation plane,
    // which we have in the z component of shift_t. Note that we can't
    // use tr_z here, because that is in model space units, but we need
    // the distance in the CS of the translation plane, to which we
    // have rotated shift_z, just above. hence:
  
    T dcp = ( 1.0 - shift_t[2] ) ;

    // this is the inverse transformation, hence:

    shift_t = - shift_t ;

    // same for the source facet, with slight modifications:

    xel_t<T,3> shift_s { fs.tr_x , fs.tr_y , fs.tr_z } ;

    // the translation is given in the model space CS, but we want
    // to apply it in the CS of the translation plane. So we have
    // to rotate the shift accordingly. Same as above for the
    // target image - but now, not inverting the sign.

    if ( fs.tp_y != 0 || fs.tp_p != 0 || fs.tp_r != 0 )
      shift_s = rotate ( xel_t<double,3> ( shift_s ) , rs_tp ) ;

    typedef r3_t < T > r_t ;

    if ( have_ttp )
    {
      if ( have_stp )
      {
        // std::cout << "case 1: ttp and stp" << std::endl ;
        r_t r_to_ttp = rotate ( r_camera , rt_tp ) ;
        tf3d_t < T , L > tf3d1 ( r_to_ttp , rt_tpi , shift_t , dcp ) ;
        r_t md_to_facet = rotate ( rs_tpi , r_facet ) ;
        tf3d_t < T , L > tf3d2 ( rs_tp , md_to_facet , shift_s ) ;
        ev = tf3d1 + tf3d2 ;
      }
      else
      {
        // std::cout << "case 2: ttp only" << std::endl ;
        r_t r_to_ttp = rotate ( r_camera , rt_tp ) ;
        r_t ttp_to_facet = rotate ( rt_tpi , r_facet ) ;
        tf3d_t < T , L > tf3d1 ( r_to_ttp , ttp_to_facet , shift_t , dcp ) ;
        ev = tf3d1 ;
      }
    }
    else
    {
      if ( have_stp )
      {
        // std::cout << "case 3: stp only" << std::endl ;
        r_t r_to_stp = rotate ( r_camera , rs_tp ) ;
        r_t stp_to_facet = rotate ( rs_tpi , r_facet ) ;
        tf3d_t < T , L > tf3d1 ( r_to_stp , stp_to_facet , shift_s ) ;
        ev = tf3d1 ;
      }
      else
      {
        // this is the simplest case: no translation in both facets
        // std::cout << "case 4: no translation" << std::endl ;
        r_t r_complete = rotate ( r_camera , r_facet ) ;
        ev = rotate_t < T , L > ( r_complete ) ;
      }
    }
  }

  // target is the regular output given in args, so there are no
  // transformations due to lens correction or translation on the
  // target side, but there may be translation on the source side,
  // which we process here - planar transformations are handled in
  // the 'environment' object.

  generic_r3 ( const facet_base & fs )
  {
    typedef r3_t < T > r_t ;
    
    // here, the facet passed in f is in the 'camera' position
    // r_camera takes us from target coordinates to model space:

    r_t r_camera = make_r3_t
      ( args.roll , args.pitch , args.yaw , false ) ;

    // this rotation takes us from model space to the source facet's
    // translation plane

    r_t rs_tp = make_r3_t ( fs.tp_r , fs.tp_p , fs.tp_y , true ) ;

    // this rotation takes us back to model space

    r_t rs_tpi = make_r3_t ( fs.tp_r , fs.tp_p , fs.tp_y , false ) ;

    // this rotation takes us to the source facet's CS

    r_t r_facet = make_r3_t ( fs.roll , fs.pitch , fs.yaw , true ) ;

    bool have_stp = ( fs.tr_x != 0 || fs.tr_y != 0 || fs.tr_z != 0 ) ;

    xel_t<T,3> shift_s { fs.tr_x , fs.tr_y , fs.tr_z } ;

    // the translation is given in the model space CS, but we want
    // to apply it in the CS of the translation plane. So we have
    // to rotate the shift accordingly.

    if ( fs.tp_y != 0 || fs.tp_p != 0 || fs.tp_r != 0 )
      shift_s = rotate ( xel_t<double,3> ( shift_s ) , rs_tp ) ;

    if ( have_stp )
    {
      // std::cout << "case 3: stp only" << std::endl ;
      r_t r_to_stp = rotate ( r_camera , rs_tp ) ;
      r_t stp_to_facet = rotate ( rs_tpi , r_facet ) ;
      tf3d_t < T , L > tf3d1 ( r_to_stp , stp_to_facet , shift_s ) ;
      ev = tf3d1 ;
    }
    else
    {
      // this is the simplest case: no translation in both facets
      // std::cout << "case 4: no translation" << std::endl ;
      r_t r_complete = rotate ( r_camera , r_facet ) ;
      ev = rotate_t < T , L > ( r_complete ) ;
    }
  }

  template < typename I , typename O >
  void eval ( const I & in , O & out )
  {
    ev.eval ( in , out ) ;
  }
} ;

// this functor handles the transformation of 2D model space coordinates
// pertaining to an output image to 3D ray coordinates. The speciality
// here is that the output is recreating one of the facets, so there
// is a '--single' argument, or we're looking at control points where
// we want to figure out ray coordinates for coordinates pertaining to
// a facet image - which is needed for control point inspection. The
// functor's argument receives a facet_spec object and gleans the
// transformation steps from it. So, again, the facet_spec we process
// here refers to the *target*, and not to the source, for which we
// normally use facet_spec objects. We want to implement the steps of
// the transformation which 'ordinary' steppers go through: that's
// the 'rise' from planar 2D to 3D ray coordinates and a rotation
// either to the source facet's CS or to the translation plane, if
// present. What we may have here, additionally, is the inverse of
// the 'normal' planar transformation due to lens correction, and the
// inverse of the translation in 3D.

template < typename T , std::size_t L >
struct tf_ex_facet
: public unary_functor < xel_t < T , 2 > ,
                         xel_t < T , 3 > ,
                         L >
{
  // 2D->2D transformation, may contain inverse of lens correction

  pto_planar < T , L , true > tf22 ;

  // 2D->3D transformation as per the facet's projection.

  grok_type < xel_t < T , 2 > , xel_t < T , 3 > , L > tf23 ;

  // 3D->3D transformation, may contain inverse translation

  generic_r3 < T , L > tf33 ;

  const bool has_2d_tf ;

  tf_ex_facet ( const facet_base & ft , const facet_base & fs )
  : tf22 ( ft ) ,
    tf23 ( roll_out_23 < T , L > ( ft.projection ) ) ,
    tf33 ( ft , fs ) ,
    has_2d_tf ( ft.has_2d_tf )
  { }

  // eval puts the three steps together and provides a ray coordinate
  // in the source facet's CS.

  template < typename I , typename O >
  void eval ( const I & _in , O & out )
  {
    I in ;
    if ( has_2d_tf )
      tf22.eval ( _in , in ) ;
    else
      in = _in ;

    tf23.eval ( in , out ) ;
    tf33.eval ( out , out ) ;
  }
} ;

template < int NCH ,
           template < typename , std::size_t , bool > class STP ,
           typename SYN >
void fuse ( int ninputs )
{
  typedef zimt::xel_t < float , NCH > px_t ;
  typedef simdized_type < px_t , LANES > px_v ;

  typedef grok_get_t < float , 3 , 2 , LANES > gg_t ;
  typedef grok_get_t < float , 9 , 2 , LANES > gg9_t ;

  typedef environment < float , float , NCH , LANES > env_t ;

  typedef fusion_t < float , NCH , 2 , LANES , float , 3 , SYN > fs_t ;
  typedef fusion_t < float , NCH , 2 , LANES , float , 9 , SYN > fs9_t ;

  typedef r3_t < double > r_t ;

  // we'll need some 3D rotations for the 'steppers':

  std::vector < r_t > basis_v ;

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

  std::vector < env_t > env_v ;

  // r_camera takes us from target coordinates to model space:

  r_t r_camera = make_r3_t
    ( args.roll , args.pitch , args.yaw , false ) ;

  for ( int i = 0 ; i < args.nfacets ; i++ )
  {
    auto const & fct ( args.facet_spec_v [ i ] ) ; // shorthand

    // first we create an environment object and push it to env_v.
    // the environment object takes care of redundancy: if the facet
    // image is already in RAM, it won't be loaded again.

    env_v.push_back ( env_t ( fct ) ) ;

    // for 'standard' steppers (encoded in template argument STP)
    // we need the rotation from the camera's to each facet's CS:

    // this rotation takes us to the source facet's CS

    r_t r_facet = make_r3_t ( fct.roll , fct.pitch , fct.yaw , true ) ;

    // and this is the combination of both rotations, taking us
    // from the camera's CS to the facet's

    r_t r_complete = rotate ( r_camera , r_facet ) ;

    basis_v.push_back ( r_complete ) ;
  }

  // now we have a set of cases to handle, depending on the type
  // of job we're running. We have to consider three points:

  // - ninputs is either 3 or 9 - 3 for direct interpolation
  //   and 9 for 'twining'

  // - we may have only a single source facet, in which case we
  //   can handle the job with leaner code

  // - facets may come with or without translation parameters
  //   or lens correction parameters

  // - the target may be an 'ordinary' image specified in 'args'
  //   or a 'single' image recreating a facet, which may make it
  //   necessary to add inverse translation and lens correction
  //   to recreate the facet as faithfully as possible.

  // we need to determine whether to use a generic stepper or
  // the stepper encoded in STP - the latter only handles a
  // single rotation directly to the target image, which can't
  // have attributes like lens correction or tanslation - this
  // is currently only possible when passing --single. If we do
  // have a --single argument, we use the generic stepper
  // if the facet has lens correction or translation active,
  // and set generic_target to this effect.

  bool generic_target = false ;
  bool generic_source = false ;

  // TODO: undo facet's 'brighten' when rendering a 'single'

  if ( args.single != -1 )
  {
    // let's look at the 'single' facet. If it has lens
    // correction or translation active, we set 'generic_target'
    // to trigger the use of a generic stepper.

    auto const & fct ( args.facet_spec_v [ args.single ] ) ;

    if ( fct.has_2d_tf || fct.has_translation )
      generic_target = true ;
  }

  if ( ninputs == 3 )
  {
    // ninputs == 3 means that we're not using twining, and the
    // steppers produce 'ordinary' rays with three components.

    if ( args.solo != -1 )
    {
      // special case: use only one facet. We don't need a synopsis-
      // -forming object, and we don't need normalized ray coordinates:
      // Since we don't compare the z component of the ray as quality
      // criterion (which we'd do for multiple facets) we can do without
      // the normalization.

      if ( args.verbose )
        std::cout << "using single-facet rendering" << std::endl ;

      int f = args.solo ;
      const auto & fct ( args.facet_spec_v [ f ] ) ;

      // if the source facet has translation parameters, we also use a
      // generic stepper. lens correction on the source side is handled
      // by the 'environment' object, so we don't need a generic stepper
      // for that.

      generic_source = fct.has_translation ;

      if ( generic_source || generic_target )
      {
        // source or target have lens correction or translation.
        // we use a 'generic stepper'. The target geometry in 'args'
        // is contained in it's base class, facet_base, so we can
        // simply pass 'args' as target facet geometry in the call
        // to tf_ex_facet.

        generic_stepper < float , LANES , false > get_ray
          ( args.width , args.height ,
            args.x0 , args.x1 , args.y0 , args.y1 ,
            0 , 0 , tf_ex_facet < float , LANES > ( args , fct ) ) ;

        work ( get_ray , env_v[f] ) ;
      }
      else
      {
        // neither source nor target have translation or lens control
        // set, so we can use the 'fast lane', using the type of stepper
        // encoded in 'STP' in the calling function.

        STP < float , LANES , false > get_ray
          ( basis_v[f][0] , basis_v[f][1] , basis_v[f][2] ,
            args.width , args.height ,
            args.x0 , args.x1 , args.y0 , args.y1 ) ;

        work ( get_ray , env_v[f] ) ;
      }
    }
    else
    {
      // there are several facets. The strategy is the same as above,
      // but now repeated for each source facet. We need normalized
      // ray coordinates, in order to have a simple criterion for the
      // synopsis-forming object (maximal z value) rather than having
      // to perform more involved calculations (like the angle vs. the
      // z axis).

      // we'll have several get_t objects, one for each facet,
      // producing rays. These get_t objects are 'grokked' to erase
      // ther specific type (STP) and stored in this vector:

      std::vector < gg_t > get_v ;

      for ( int i = 0 ; i < args.nfacets ; i++ )
      {
        const auto & fct ( args.facet_spec_v[i] ) ; // shorthand

        generic_source = fct.has_translation ;

        if ( generic_source || generic_target )
        {
          generic_stepper < float , LANES , true > get_ray
            ( args.width , args.height ,
              args.x0 , args.x1 , args.y0 , args.y1 ,
              0 , 0 , tf_ex_facet < float , LANES > ( args , fct ) ) ;
          get_v.push_back ( get_ray ) ;
        }
        else
        {
          STP < float , LANES , true > get_ray
            ( basis_v[i][0] , basis_v[i][1] , basis_v[i][2] ,
              args.width , args.height ,
              args.x0 , args.x1 , args.y0 , args.y1 ) ;
          get_v.push_back ( get_ray ) ;
        }
      }

      // for now, we use a hard-coded synopsis-forming object

      SYN vs ( env_v ) ;

      // now we can create the fusion_t object

      fs_t fs ( get_v , vs ) ;

      // now we call the final 'work' template which uses the fusion_t
      // which directly produces a pixel value, so we use a pass_through
      // act functor, rather than having the interpolator step in the
      // act functor.

      zimt::pass_through < float , NCH , LANES > act ;

      work ( fs , act ) ;
    }
  }
  else // ninputs == 9
  {
    // this is the code we use for twining. instead of processing
    // simple 3D rays, we process 'ninepacks'. The structure of the
    // code is very similar to the code without twining, but now we
    // 'wrap' the steppers in 'deriv_stepper' objects which contain
    // three separate steppers with slightly different parameters,
    // producing three slightly different rays which are combined in
    // the deriv_stepper's output - a 'ninepack'.

    if ( args.solo != -1 )
    {
      // special case: use only one facet.
      // Note how we use a stepper which does not normalize it's result:
      // Since we don't compare the z component of the ray as quality
      // criterion (which we'd do for multiple facets) we can do without
      // the normalization.

      if ( args.verbose )
        std::cout << "using single-facet rendering" << std::endl ;

      int f = args.solo ;
      const auto & fct ( args.facet_spec_v [ f ] ) ;

      generic_source = fct.has_translation ;

      if ( generic_source || generic_target )
      {
        deriv_stepper < float , LANES , generic_stepper , false > get_ray
          ( args.width , args.height ,
            args.x0 , args.x1 , args.y0 , args.y1 ,
            tf_ex_facet < float , LANES > ( args , fct ) ) ;
        twine_t < NCH , LANES > twenv ( env_v[f] , args.twine_spread ) ;
        work ( get_ray , twenv ) ;
      }
      else
      {
        deriv_stepper < float , LANES , STP , false > get_ray
            ( basis_v[f][0] , basis_v[f][1] , basis_v[f][2] ,
              args.width , args.height ,
              args.x0 , args.x1 , args.y0 , args.y1 ) ;
        twine_t < NCH , LANES > twenv ( env_v[f] , args.twine_spread ) ;
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

        generic_source = fct.has_translation ;

        if ( generic_source || generic_target )
        {
          deriv_stepper < float , LANES , generic_stepper , true > get_ray
            ( args.width , args.height ,
              args.x0 , args.x1 , args.y0 , args.y1 ,
              tf_ex_facet < float , LANES > ( args , fct ) ) ;
          get_v.push_back ( get_ray ) ;
        }
        else
        {
          deriv_stepper < float , LANES , STP , true > get_ray
              ( basis_v[i][0] , basis_v[i][1] , basis_v[i][2] ,
                args.width , args.height ,
                args.x0 , args.x1 , args.y0 , args.y1 ) ;
          get_v.push_back ( get_ray ) ;
        }
      }

      // we set up the synopsis-forming object and introduce it to
      // the fusion_t object - which is now used via it's ninepack-
      // processing member function. The synopsis-forming object
      // handles the twining as well.

      SYN vs ( env_v ) ;

      fs9_t fs ( get_v , vs ) ;

      // now we call the final 'work' template which uses the fusion_t
      // which directly produces a pixel value, so we use a pass_through
      // act functor.

      zimt::pass_through < float , NCH , LANES > act ;

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
  if ( args.synopsis == "hdr_merge" )
  {
    fuse < NCH , STP , hdr_merge_syn < NCH > > ( ninputs ) ;
  }
  else if ( args.synopsis == "panorama" )
  {
    // assume one- and three-channel images are without alpha
    // channel, two- and four-channel images are with alpha channel.

    if constexpr ( NCH == 1 || NCH == 3 )
      fuse < NCH , STP , voronoi_syn < NCH > > ( ninputs ) ;
    else
      fuse < NCH , STP , voronoi_syn_plus < NCH > > ( ninputs ) ;
  }
  else
  {
    assert ( false ) ;
  }
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
// #ifndef NARROW_SCOPE
    case CYLINDRICAL:
      roll_out < NCH , cylindrical_stepper > ( ninputs ) ;
      break ;
// #endif
    case RECTILINEAR:
      roll_out < NCH , rectilinear_stepper > ( ninputs ) ;
      break ;
#ifndef NARROW_SCOPE
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
: public zimt::unary_functor < crd_t , px_t , LANES >
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
    conclude_cycle() ;
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
