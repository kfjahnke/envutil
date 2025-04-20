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

#include "zimt/zimt.h"

#if defined(ENVUTIL_STEPPER_H) == defined(HWY_TARGET_TOGGLE)
  #ifdef ENVUTIL_STEPPER_H
    #undef ENVUTIL_STEPPER_H
  #else
    #define ENVUTIL_STEPPER_H
  #endif

HWY_BEFORE_NAMESPACE() ;
BEGIN_ZIMT_SIMD_NAMESPACE(project)

// This header provides a set of 'steppers': concrete implementations
// of objects which can serve as data source for zimt::process. The
// 'steppers' are quite specific: the provide 3D 'ray' coordinates
// to sample an environment so that these coordinates will pick pixels
// to populate a given 2D image in a specific projection.
// The steppers only produce the coordinates - providing pixel values
// to populate the image has to be done in a second step, with the
// 'act' functor of zimt::process - or by 'suffixing' the stepper with
// further processing steps. The latter method 'pulls' functionality
// into the stepper. There's ready-made code for the purpose in
// zimt/get.h, class suffixed_t: this class suffixes a get_t type
// object with an act type functor and produces a get_t type object
// as result. Steppers are regular get_t type objects, so they can
// be 'suffixed' like any other get_t type.

// The steppers have the added feature of arbitrary orientation. Their
// c'tors take three unit vectors for the three axes of the stepper's
// orientation - if you want no rotation, pass (1,0,0), (0,1,0) and
// (0,0,1), and if you want rotation, rotate all three vectors with
// e.g. the same rotational quaternion. The resulting steppers will
// have the rotation 'built-in', so if you do e.g. a cut-out from
// an environment, the rotation can be used to set up the orientation
// of the 'virtual camera' for the cut-out.
// Adding this feature comes at little additional cost compared to
// what an implementation without rotation would use, but increases
// the usefulness. The alternative - using steppers without built-in
// rotation and starting the 'act' functor with a rotation - e.g. with
// a rotational quaternion - is slower, because the rotation has to be
// applied to every ray coordinate, and this takes quite a few CPU
// cycles. The built-in rotation used here can also exploit certain
// invariants and mathematical tricks which make it more efficient.
// The intended use is for extraction of partial images from
// environments, reprojection, panorama viewing and the likes.
// A concrete program using this code is 'extract.cc' in this repo.

// I have finally grasped the principle how to avoid 'rotating'
// coordinates at all: if one first calculates the 3D ray for the
// 'archetypal' unrotated 2D manifold, it's three components are
// used as factors for the three rotated vectors representing the
// re-oriented, rotated manifold. If no invariants apply, this
// needs nine multiplications and six additions, but many steppers
// have invariant parts throughout the processing of a segment.
// This hints at further possibilities: instead of using an
// orthonormal set of vectors, any set of three vecors might be
// used for all kinds of effects. The steppers take the vectors
// as arguments, so there is ample room for experiments.

// My realization is simply understanding the application of a
// rotational matrix. The key to understanding is that there 'is
// no rotation' - what's happening is that a given 3D vector is
// put into a different coordinate system with the same origin.
// This can be interpreted as the *result of a rotation*, but it
// does not mean that anything *actually does rotate*. If something
// actually does rotate, then rotational matrices and quaternions
// can be used to calculate the result of the rotation(s), but the
// essential fact is that the original vector is simply put into
// a different coordinate system. The typical modus operandi of using
// a stepper is to use a target which is 'draped' on the 'archetypal'
// 2D manifold representing it's projection. So the target is
// thought to be upright and it's center at unit distance forward;
// the target's x axis is parallel to the 'right' axis (x axis in
// lux convention) and the target's y axis is parallel to the 'down'
// axis (y axis in lux convention). Rather than thinking along the
// lines of a rotated camera, we can operate in the target image's
// coordinate system and 'pull data' from a rotated environment.
// This is the reverse transformation typically used for geometrical
// transformations in a panorama stitching and viewing context:
// we start out with a planar discrete coordinate pertaining to
// a pixel in the target image and ask: "which source image pixel(s)
// should affect this target pixel?" the location of the source
// pixel(s) is obtained by appropriate geometrical transformation,
// and how the source pixels in the vicinity of this location are
// combined is a matter of the chosen interpolation method. When
// stitching several images together, an additional criterion is
// the preference - or weight - which is given to what can be gleaned
// from each contributing source image in this fashion.

// Steppers also take limits for the planar coordinates, to define
// the area which will actually be sampled. The limits are given
// in model space units, which requires knowledge about the extent
// of the 2D archetype in question. An alternative would be to use
// texture coordinates, but to work with texture coordinates, we'd
// have to move to model space units for the internal caclculations,
// so 'feeding' model space units is more efficient. envutil
// can accept a horizontal field of view on the command line and
// calculates the extent's limits from it. This should be the most
// common way to go - the extent values produced in this way are
// symmetrical and isotropic. Here, in the implementation, we
// process the extent values, because they allow more flexibility.

// Steppers as such do not have a 'notion' of an image which they
// might be used for - beyond the limits, which can be passed to
// their c'tor. Within these limits, they will provide a uniform
// sampling over the 2D manifold, but emit 3D rays which coincide
// with the 2D sample points. Note that the sampling is uniform
// in relation to a 'flattened-out' 2D manifold. If we're using
// a spherical_stepper, the 'archetypal form' is a spherical
// surface which can't be sampled uniformly. But the lat/lon
// representation of this spherical surface can be sampled
// uniformly, and this is what the spherical_stepper does.
// cylindrical_stepper has a cylindrical surface as it's
// archetypal form, but this can be 'rolled out' without
// deforming it, and the other singe-image steppers also have
// 'flat' archetypes because they are, conceptually, projections
// to the plane anyway. cubemap_stepper is slightly different.
// The corresponding archetype is the surface of a cube, but
// we use a 1X6 stripe representation. If the cubemap_stepper
// is set up for a minimal planar representation (section_md
// == 2 and refc_md == 1), the stripe has model space extent
// of (2,12), and the sampling of the environment will be as
// 'compact' as possible. With extra space around the cube faces
// proper, some areas near the cube's edges will be sampled
// more often (because they are 'seen' from several of the
// 'padded' cube faces. The set of rays generated from a 'tight'
// cubemap_stepper is a 'quite' even sampling of the environment:
// The distance between neighbouring rays does vary because of
// the rectilinear projection used for the cube faces, but the
// variation is nowhere near as extreme as what happens near the
// poles with a spherical_stepper, or at the 'back pole' of a
// fisheye_stepper. So to sample an entire environment I think
// a cubemap_stepper is quite a good choice, and combined with
// an interpolator which takes the variations of ray distance into
// account (e.g. twining), the results should be good. Other schemes
// to 'cover' an entire 360X180 environment as uniformly as possible
// require more calculation (e.g. triangulations) and/or do not map
// to a uniform grid. A disadvantage of the 1X6 stripe representation
// are the discontinuities between the six partial images, which
// can't be avoided: there is no way to form a stripe of six partial
// images without unrelated edges 'colliding'. This is the reason
// for using 'padding': if the padding is sufficiently wide (so that
// the interpolator's support remains within a single padded image)
// the discontinuity has no effect. To merely sample an environment
// with roughly equally spaced rays, a 'tight' cubemap_stepper is a
// good choice, though - this fact can be exploited e.g. for 'maps'
// of the environment which indicate properties of parts of the
// environment - often scaled-down versions of the environment used
// as masks to avoid unnecessary calculations e.g. for 'void' areas.
// With the underlying rectangular structure, such schemes can be
// used together with tiling schemes, marking which tiles will be
// used in an intended manipulation. There is a variant of class
// cubemap_stepper which uses a planar trasformation on the in-face
// coordinates to mitigate the shortcomings of the rectilinear
// projection which is normally used for cubemaps, using atan/tan
// and a linear factor. This is quite close to an optimal compromise
// between resolution, computability and capacity: biatan6_stepper.
// See there for more comments.

// we set up a base class for 'steppers' which encodes the '2D part'
// of setting up a sampling: We have a horizontal and a vertical
// range of 2D coordinates pertaining to the 2D manifold 'draped'
// in 'model space'. I refer to this 2D manifold as the 'archetypal'
// representation, and it is clearly defined for several projections:
// e.g. for spherical projection, it's a spherical surface with unit
// radius, with spherical coordinates centered on (0,0,1) - this
// is where the 'planar' comes from. The 'draping' basically fixes
// the scale: the archetypal representation has unit radius, and
// the width and height parameters, together with the 'extent' of
// the draped manifold, determine the distance between the sampling
// points on the 2D manifold taken as a plane.
// stepper_base is very similar to zimt's linspace_t, which might
// be used instead. The focus here is on precision and parameterization
// which fits the intended use - envutil has the 'extent' of 2D manifolds
// readily at hand and uses 'edge-to-edge' semantics for the extent.

template < typename T ,     // fundamental type
           std::size_t L >  // lane count
struct stepper_base
{
  typedef T value_t ;
  static const std::size_t size = 3 ;
  typedef zimt::xel_t < T , 3 > crd3_t ;
  typedef zimt::simdized_type < crd3_t , L > crd3_v ;
  typedef zimt::xel_t < T , 2 > crd2_t ;
  typedef zimt::simdized_type < crd2_t , L > crd2_v ;
  typedef zimt::simdized_type < T , L > f_v ;
  typedef zimt::xel_t < int , 2 > crd_t ;
  typedef zimt::simdized_type < crd_t , L > crd_v ;
  typedef typename crd_v::value_type crd_ele_v ;

  // delta is the increase from one (SIMDized) batch of data to
  // the next, in the horizontal. This equals the increase from
  // one draped pixel to the next times the SIMD lane count.

  const T delta ;

  // 'planar' is the current SIMDized planar coordinate. It's vertical
  // component remains constant throughout a line, the horizontal
  // component varies. The invariant vertical is what we aim to
  // exploit with the stepper: oftentimes the output of the stepper,
  // a 3D ray coordinate, can be calculated more efficiently if
  // the vertical component remains constant.

  crd2_v planar ;

  // some helper variables, used for efficiency. The f... values
  // are basically scaling factors - internally, we use doubled
  // discrete coordinates, but the results are in model space
  // scale and in floating point. Why doubled discrete coordinates?
  // Because they can be used to represent half unit steps
  // precisely, and this is useful, because we're often dealing
  // with half unit steps: with edge-to-edge semantics, the extent
  // is half a sampling step further out than the marginal pixels,
  // and the center of even-shaped images is at a half unit as well.

  const T fx0 , fx1 , fy0 , fy1 , bias_x , bias_y ;

  // this is a SIMD constant holding values 0, 2, 4, ...
  // so, twice what the ordinary 'iota' would provide, hence the name.

  const crd_ele_v iota2 ;

  // width and height give the number of sample steps over the entire
  // surface covered by the image's extent. These two values define
  // the stepper's - and it's generated sampling's - resolution.
  // The width and height are independent of the extent and of
  // each other, so anisotropic sampling is allowed and a given
  // extent can be sampled with any horizontal and vertical
  // resolution.

  const int width , height ;

  // stepper_base's standard c'tor receives the width and height,
  // and the image's extent, given as the extrema first in
  // horizontal, then in vertical direction. Note that we use
  // edge-to-edge semantics to place the planar coordinates
  // inside the extent.
  // Note the two recently added parameters '_bias_x' and '_bias_y':
  // These are used for twining. They are small offsets in the
  // horizontal/vertical which are added to the planar coordinate
  // to produce slightly offsetted rays. In the first implementation
  // I used an offset of one 'canonical' step, but this can 'land'
  // the planar coordinate in the domain of another pixel, which,
  // for some steppers, may correspond with a completely different
  // ray. The bias values I use now are chosen to keep the offsetted
  // steppers in the same pixel's domain (so the step is not 1, but
  // the bias value instead, which is chosen smaller than .5 - I
  // use .25. When the resulting rays are used to form the derivatives,
  // The 'raw' value obtained from differencing needs to be pulled up
  // to correspond with a canonical step, so after the differencing,
  // the reciprocal of the bias used here is applied. Note that this
  // does not produce additional computations on this side apart from
  // in the initialization of the planar coordinate in 'init'.

  stepper_base ( int _width , int _height ,
                 T _a0 , T _a1 , T _b0 , T _b1 ,
                 T _bias_x = 0 , T _bias_y = 0 )
  : width ( _width ) ,
    height ( _height) ,
    fx1 ( _a1 / ( 2.0 * _width ) ) ,
    fx0 ( _a0 / ( 2.0 * _width ) ) ,
    fy1 ( _b1 / ( 2.0 * _height ) ) ,
    fy0 ( _b0 / ( 2.0 * _height ) ) ,
    bias_x ( _bias_x * ( _a1 - _a0 ) / T ( _width ) ) ,
    bias_y ( _bias_y * ( _b1 - _b0 ) / T ( _height ) ) ,
    delta ( L * ( _a1 - _a0 ) / T ( _width ) ) ,
    iota2 ( 2 * crd_ele_v::iota() )
  { }

  // init sets up the internal state ('planar') resulting from
  // the initial discrete coordinate for the current segment.
  // This is called by the 'init' member function in the
  // concrete stepper implementations. Note how the incoming
  // discrete coordinate is doubled in processing, and then
  // mapped to model space scale by application of the scaling
  // factors. Note also the way how the planar coordinate is
  // calculated: by pulling both extrema into the calculation
  // so that the corresponding opposite term disappears when
  // the calculation 'hits the margin', we make certain that
  // the resulting values are precisely in-range. This takes
  // a few more cycles than the 'starting point + x * delta'
  // approach, but with the latter, with large x, we might
  // well end up a good way off the upper end of the range.

  void init ( const crd_t & crd )
  {
    auto ll0 = iota2 + T ( crd[0] * 2 + 1 ) ;

    planar[0] = bias_x + ll0 * fx1 + ( T ( 2 * width ) - ll0 ) * fx0 ;

    auto ll1 = crd[1] * 2 + 1 ;

    planar[1] = bias_y + ll1 * fy1 + T ( 2 * height - ll1 ) * fy0 ;
  }

  // increase modifies it's argument to contain the next value
  // here we add 'delta' - an alternative would be to keep state
  // representing the current discrete coordinate and to derive
  // the current planar[0] value from it. This would be more
  // precise, but since we calculate a new starting point for
  // planar[0] with every segment (with the call to 'init') the
  // drift should be negligible.
  // increase is called by the concrete stepper implementation's
  // 'increase' member function. Note that increase is 'dumb'
  // and will produce values outside the extent eventually;
  // the conrete stepper implementation will take care of that.

  void increase()
  {
    planar[0] += delta ;
  }
} ;

template < typename T ,     // fundamental type
           std::size_t L ,  // lane count
           bool normalize = true >
struct generic_stepper
: public stepper_base < T , L >
{
  typedef grok_type < xel_t < float , 2 > ,
                      xel_t < float , 3 > ,
                      L > tf_t ;

  typedef stepper_base < T , L > base_t ;
  using typename base_t::value_t ;
  using base_t::size ;
  using typename base_t::crd_t ;
  using typename base_t::crd3_t ;
  using typename base_t::crd3_v ;
  using typename base_t::f_v ;
  using base_t::planar ;
  using base_t::init ;
  using base_t::increase ;

  // c'tor arguments:

  // _width and _height are the extent, in pixels, of the field
  // of sample values which this get_t will produce - the
  // 'notional' shape for zimt::process. The width and height
  // don't need to result in an isotropic sampling; whichever
  // the values are, the spread of sampling locations will be
  // uniform in the horizontal and in the vertical. But note
  // that the sampling does not start at (_a0,_b0) but half a
  // sampling step 'inwards' - and it also ends half a sampling
  // step inwards from (_a1,_b1): we're using 'edge-to-edge
  // semantics'

  // The next four arguments are the 'extent' of the 2D manifold
  // which will be sampled, in model space units. The default
  // values will spread the sampling locations evenly over an
  // entire generic surface.

  // the 2D part of the work - producing the 'planar' coordinate
  // is handled by the base type.

  // this functor will handle the transformation from planar
  // target coordinates all the way to 3D ray coordinates in the
  // source facet's CS.

  tf_t tf ;

  generic_stepper ( // crd3_t _xx , crd3_t _yy , crd3_t _zz ,
                    int _width ,
                    int _height ,
                    T _a0 ,
                    T _a1 ,
                    T _b0 ,
                    T _b1 ,
                    T _bias_x ,
                    T _bias_y ,
                    const tf_t & _tf
                  )
  : base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 ,
             _bias_x , _bias_y ) ,
    tf ( _tf )
  { }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and this
  // function infers the start value from it.

  void init ( crd3_v & trg , const crd_t & crd )
  {
    // the base type's 'init' initializes 'planar'

    init ( crd ) ;

    // 'tf' transforms the planar coordinate to a 3D ray

    tf.eval ( planar , trg ) ;
    if constexpr ( normalize )
    {
      trg /= norm ( trg ) ;
    }
  }

  // 'capped' variant. This is only needed if the current segment is
  // so short that no vectors can be formed at all. We fill up the
  // target value with the last valid datum. The code for this
  // function is the same for all steppers, but since it uses
  // 'init', which is specific to a given stepper, and since I
  // don't want to use a virtual init function, I repeat it
  // in all steppers.

  void init ( crd3_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    init ( trg , crd ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }

  // increase modifies it's argument to contain the next values:
  // the next 3D 'ray' coordinates in model space units, which are
  // used as input to the 'act' functor invoked by zimt::process.

  void increase ( crd3_v & trg )
  {
    // the base class increase increases planar[0]

    increase() ;

    // 'tf' transforms the planar coordinate to a 3D ray

    // crd3_v help ;
    tf.eval ( planar , trg ) ;
    if constexpr ( normalize )
    {
      trg /= norm ( trg ) ;
    }
    // trg = help[0] * xx + help[1] * yy + help[2] * zz ;
  }

  // 'capped' variant. This is called after all vectors in the current
  // segment have been processed. If we'd 'go over the lanes', we might
  // only set entries below cap and rely on the higher lanes holding
  // viable values from processing the last full vector, but it's
  // probably more efficient to fill the target with a vector operation
  // and stuff unconditionally.

  void increase ( crd3_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    increase ( trg ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }
} ;

// struct spherical_stepper is an object to be used as a zimt get_t.
// It provides 3D ray coordinates for a sampling of the surface of a
// sphere with unit radius. What's special here is that the sampling
// works in an arbitrarily rotated coordinate system - or, if you
// think of the sphere as having it's own coordinate system, it works
// with an arbitrarily rotated sphere. Because the code can exploit
// certain invariants, this is more efficient than producing the
// sampling from a 'normal' sphere and subsequently rotating the
// 3D coordinates - e.g. with a rotational quaternion. The
// mathematics behind the code were gleaned from here:
// https://math.stackexchange.com/questions/73237/parametric-equation-of-a-circle-in-3d-space
// The key point is that, following a latitude line on an arbitrarily
// rotated sphere, the latitude line is still a circle which can be
// produced by a parametric function of the angle. The parameters of
// the equation are established and the fact that a zimt::get_t only
// varies it's output along the 'hot axis' while it processes a segment
// is exploited to good effect, reducing the calculatory expense to
// an invocation of sincos, nine additions and six multiplications
// for each batch of produced coordinates (counting SIMD operations;
// the vector width may result in the nominal vectors being represented
// by more than one register's worth).
// This object is intended for use with lux: lux already uses a similar
// scheme for rectilinear targets, but I did not realize that the
// method can be extended to other target projections.

template < typename T ,     // fundamental type
           std::size_t L ,  // lane count
           bool normalize = true >
struct spherical_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
  using typename base_t::value_t ;
  using base_t::size ;
  using typename base_t::crd_t ;
  using typename base_t::crd3_t ;
  using typename base_t::crd3_v ;
  using typename base_t::f_v ;
  using base_t::planar ;
  using base_t::init ;
  using base_t::increase ;

  // we have three 3D vectors of unit length. They are orthogonal and
  // can be produced e.g. by applying a rotational quaternion to the
  // three cardinal vectors (1,0,0), (0,1,0), (0,0,1)

  const crd3_t xx ; 
  const crd3_t yy ;
  const crd3_t zz ;

  // some helper variables, used for efficiency.

  f_v x , r , y , z ;
  crd3_v xxx , yyy , zzz ;

public:

  // c'tor arguments:

  // _xx, _yy and _zz are 3D directional vectors corresponding to
  // (1,0,0), (0,1,0) and (0,0,1) with the rotation applied.
  // If there is no rotation, we could use a simplified object,
  // but this implementation is perfectly general. An alternative
  // way to introduce these three vectors would be to pass in
  // a rotational quaternion and rotate the three base vectors
  // by it, but since this creates a dependency on the relevant
  // quaternion code, I prefer to pass in the three vectors.
  // The vectors are expected to be normalized.
  // To put it differently: the three vectors _xx, _yy and _zz
  // form the orthonormal basis for the rotated coordinate
  // system in which we operate.

  // _width and _height are the extent, in pixels, of the field
  // of sample values which this get_t will produce - the
  // 'notional' shape for zimt::process. The width and height
  // don't need to result in an isotropic sampling; whichever
  // the values are, the spread of sampling locations will be
  // uniform in the horizontal and in the vertical. But note
  // that the sampling does not start at (_a0,_b0) but half a
  // sampling step 'inwards' - and it also ends half a sampling
  // step inwards from (_a1,_b1): we're using 'edge-to-edge
  // semantics'

  // The last four arguments are the 'extent' of the 2D manifold
  // which will be sampled, in model space units. The default
  // values will spread the sampling locations evenly over an
  // entire spherical surface.

  // the 2D part of the work - producing the 'planar' coordinate
  // is handled by the base type.

  spherical_stepper ( crd3_t _xx , crd3_t _yy , crd3_t _zz ,
                      int _width ,
                      int _height ,
                      T _a0 ,
                      T _a1 ,
                      T _b0 ,
                      T _b1 ,
                      T _bias_x = 0 ,
                      T _bias_y = 0
                     )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 ,
             _bias_x , _bias_y )
  { }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and this
  // function infers the start value from it.

  void init ( crd3_v & trg , const crd_t & crd )
  {
    // the base type's 'init' initializes 'planar'

    init ( crd ) ;

    // these values pertain to the vertical:

    sincos ( planar[1] , y , r ) ;

    // these values pertain to the horizontal:

    sincos ( planar[0] , x , z ) ;

    // for efficiency, we form a few products now - they remain
    // constant for the entire segment (so, until 'init' is called
    // again by zimt::process)

    xxx = xx * r ;
    yyy = yy * y ;
    zzz = zz * r ;

    // now we can initialize the target value

    trg = xxx * x + zzz * z + yyy ;
  }

  // 'capped' variant. This is only needed if the current segment is
  // so short that no vectors can be formed at all. We fill up the
  // target value with the last valid datum. The code for this
  // function is the same for all steppers, but since it uses
  // 'init', which is specific to a given stepper, and since I
  // don't want to use a virtual init function, I repeat it
  // in all steppers.

  void init ( crd3_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    init ( trg , crd ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }

  // increase modifies it's argument to contain the next values:
  // the next 3D 'ray' coordinates in model space units, which are
  // used as input to the 'act' functor invoked by zimt::process.

  void increase ( crd3_v & trg )
  {
    // the base class increase increases planar[0]

    increase() ;

    // update x and z only - y and r are already fixed in the
    // internal xxx, yyy and zzz variables.

    sincos ( planar[0] , x , z ) ;

    // finally, calculate the new target value

    trg = xxx * x + zzz * z + yyy ;
  }

  // 'capped' variant. This is called after all vectors in the current
  // segment have been processed. If we'd 'go over the lanes', we might
  // only set entries below cap and rely on the higher lanes holding
  // viable values from processing the last full vector, but it's
  // probably more efficient to fill the target with a vector operation
  // and stuff unconditionally.

  void increase ( crd3_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    increase ( trg ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }
} ;

// cylindrical_stepper is very similar to spherical_stepper. The
// difference is the different treatment of the vertical, and that
// the x and z components aren't scaled with the radius, because
// it's always 1.0
// TODO: all rays along a horizontal line have the same length, so
// if we calculate this length in 'init' and form it's reciprocal
// value, we obtain a factor which will normalize the result in
// 'target' without needing the costly division by the norm. The
// factor can be applied to the y component in 'init' because it
// remains static throughout the segment, x and z would need to be
// scaled after the 'sincos', but an alternative would be to scale
// xx and zz, which should produce the same result, so that we might
// obtain a normalized result without any additional operations in
// the course of processing a segment.
// another alternative: use spherical_stepper, but not with planar[1]
// but with it's atan. Since planar[1] is also only used in 'init',
// the cost would only by due once per segment.

template < typename T ,     // fundamental type
           std::size_t L ,  // lane count
           bool normalize = true >
struct cylindrical_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
  using typename base_t::value_t ;
  using base_t::size ;
  using typename base_t::crd_t ;
  using typename base_t::crd3_t ;
  using typename base_t::crd3_v ;
  using typename base_t::f_v ;
  using base_t::planar ;
  using base_t::init ;
  using base_t::increase ;

  // we have three 3D vectors of unit length. The are orthogonal and
  // can be produced e.g. by applying a rotational quaternion to the
  // three cardinal vectors (1,0,0), (0,1,0), (0,0,1)

  const crd3_t xx ; 
  const crd3_t yy ;
  const crd3_t zz ;

  // some helper variables, used for efficiency.

  f_v x , r , y , z ;
  crd3_v yyy ;
  f_v rcp_length ;

public:

  cylindrical_stepper ( crd3_t _xx , crd3_t _yy , crd3_t _zz ,
                        int _width ,
                        int _height ,
                        T _a0 ,
                        T _a1 ,
                        T _b0 ,
                        T _b1 ,
                        T _bias_x = 0 ,
                        T _bias_y = 0
                       )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 ,
             _bias_x , _bias_y )

  { }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and this
  // function infers the start value from it.

  void init ( crd3_v & trg , const crd_t & crd )
  {
    // the base type's init initializes 'planar'

    init ( crd ) ;

    // this value pertains to the vertical:

    y = planar[1] ; // constant while 'latitude' is constant

    // these values pertain to the horizontal:

    sincos ( planar[0] , x , z ) ;

    // for efficiency, we form a few products now - they remain constant
    // for the entire segment

    yyy = yy * y ;

    // now we can initialize the target value

    trg = xx * x + zz * z + yyy ;

    if constexpr ( normalize )
    {
      rcp_length = 1.0f / norm ( trg ) ;
      trg *= rcp_length ;
    }
  }

  // increase modifies it's argument to contain the next value

  void increase ( crd3_v & trg )
  {
    // the base class increase increases planar[0]

    increase() ;

    // update x and z

    sincos ( planar[0] , x , z ) ;

    // and calculate the new target value

    trg = xx * x + zz * z + yyy ;

    if constexpr ( normalize )
    {
      trg *= rcp_length ;
    }
  }

  void init ( crd3_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    init ( trg , crd ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }

  void increase ( crd3_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    increase ( trg ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }
} ;

// rectilinear_stepper is the simplest of the family, because
// the planar coordinate can be mapped more or less directly
// to the output.
// To create normalized output, it's probably most efficient to
// divide by the norm - an alternative would be to form the atan
// of the progressing x coordinate to obtain an angle, which
// might be processed as in spherical_stepper (sincos, multiply
// with basis vectors) - but that requires two transcendental
// functions, vs. the srqt, which should be faster (test!)

template < typename T ,     // fundamental type
           std::size_t L ,  // lane count
           bool normalize = true >
struct rectilinear_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
  using typename base_t::value_t ;
  using base_t::size ;
  using typename base_t::crd_t ;
  using typename base_t::crd3_t ;
  using typename base_t::crd3_v ;
  using typename base_t::f_v ;
  using base_t::planar ;
  using base_t::init ;
  using base_t::increase ;

  // we have three 3D vectors of unit length. The are orthogonal and
  // can be produced e.g. by applying a rotational quaternion to the
  // three cardinal vectors (1,0,0), (0,1,0), (0,0,1)

  const crd3_t xx ; 
  const crd3_t yy ;
  const crd3_t zz ;

  // some helper variables, used for efficiency.

  f_v x , r , y , z ;
  crd3_v ddd ; 

public:

  rectilinear_stepper ( crd3_t _xx , crd3_t _yy , crd3_t _zz ,
                        int _width ,
                        int _height ,
                        T _a0 ,
                        T _a1 ,
                        T _b0 ,
                        T _b1 ,
                        T _bias_x = 0 ,
                        T _bias_y = 0 
                       )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 ,
             _bias_x , _bias_y )
  { }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and this
  // function infers the start value from it.

  void init ( crd3_v & trg , const crd_t & crd )
  {
    // the base type's init initializes 'planar'

    init ( crd ) ;

    // for efficiency, we form 'ddd' now - it remains constant
    // for the entire segment

    ddd = yy * planar[1] + zz ;

    // now we can initialize the target value

    trg = xx * planar[0] + ddd ;

    if constexpr ( normalize )
    {
      trg /= norm ( trg ) ;
    }
  }

  // increase modifies it's argument to contain the next value

  void increase ( crd3_v & trg )
  {
    // the base class increase increases planar[0]

    increase() ;

    // calculate the new target value

    trg = xx * planar[0] + ddd ;

    if constexpr ( normalize )
    {
      trg /= norm ( trg ) ;
    }
  }

  // the capped variants:

  void init ( crd3_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    init ( trg , crd ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }

  void increase ( crd3_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    increase ( trg ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }
} ;

template < typename T >
bool near_normal ( const T & x )
{
  auto n = abs ( norm ( x ) - 1.0f ) ;
  return ( all_of ( n < .000001f ) ) ;
}

// fisheye_stepper is a bit more involved mathematically, and
// there are no handy invariances, because the fisheye projection
// is a radial function, whereas the traversal of the 2D manifold
// is a tensor, so there is no common mathematical ground and
// each sample point (or, rather, each batch of sample points)
// has to be calculated 'all the way' from the planar coordinate.

template < typename T ,     // fundamental type
           std::size_t L ,  // lane count
           bool normalize = true >
struct fisheye_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
  using typename base_t::value_t ;
  using base_t::size ;
  using typename base_t::crd_t ;
  using typename base_t::crd3_t ;
  using typename base_t::crd3_v ;
  using typename base_t::f_v ;
  using base_t::planar ;
  using base_t::init ;
  using base_t::increase ;

  // we have three 3D vectors of unit length. The are orthogonal and
  // can be produced e.g. by applying a rotational quaternion to the
  // three cardinal vectors (1,0,0), (0,1,0), (0,0,1)

  const crd3_t xx ; 
  const crd3_t yy ;
  const crd3_t zz ;

  // some helper variables, used for efficiency.

  f_v x , r , y , z ;
  crd3_v yyy ; 
  fisheye_stepper ( crd3_t _xx , crd3_t _yy , crd3_t _zz ,
                    int _width ,
                    int _height ,
                    T _a0 ,
                    T _a1 ,
                    T _b0 ,
                    T _b1 ,
                    T _bias_x = 0 ,
                    T _bias_y = 0 
                  )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 ,
             _bias_x , _bias_y )
  { }

  // for fisheye images, we can't use any handy invariants - each
  // sample point has to be calculated 'all the way', and we use
  // this function for the calculation:

  void work ( crd3_v & trg )
  {
    // 2D distance from center, equals the angle from the central
    // axis zz (the rotated 'forward' axis (0,0,1) and the same
    // as the 'latitude'

    auto a =  M_PI_2 - norm ( planar ) ;

    // angle from the horizontal, same as the 'longitude'

    auto b = atan2 ( planar[0] , planar[1] ) ;

    // these values pertain to the vertical:

    sincos ( a , z , r ) ;

    // these values pertain to the horizontal:

    sincos ( b , x , y ) ;

    // now we can initialize the target value

    trg = xx * r * x + zz * z + yy * r * y ;
  }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and this
  // function infers the start value from it.

  void init ( crd3_v & trg , const crd_t & crd )
  {
    // the base type's init initializes 'planar'

    init ( crd ) ;
    work ( trg ) ;
  }

  // increase modifies it's argument to contain the next value

  void increase ( crd3_v & trg )
  {
    // the base class increase increases planar[0]

    increase() ;
    work ( trg ) ;
  }

  void init ( crd3_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    init ( trg , crd ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }

  void increase ( crd3_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    increase ( trg ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }
} ;

// even more involved than fisheye_stepper, this stepper for
// stereographic targets has to deal with the peculiarities of
// the stereographic projection, resulting in a more complex
// temp for 'a', the equivalent of the latitude: with fisheye
// images, this is simply the distance-to-center of the planar,
// but for stereographic images, we need to calculate our way
// by first looking at the line from the unit sphere's pole
// to the target point, then intersect that with the unit
// sphere to get the ray.

template < typename T ,     // fundamental type
           std::size_t L ,  // lane count
           bool normalize = true >
struct stereographic_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
  using typename base_t::value_t ;
  using base_t::size ;
  using typename base_t::crd_t ;
  using typename base_t::crd3_t ;
  using typename base_t::crd3_v ;
  using typename base_t::f_v ;
  using base_t::planar ;
  using base_t::init ;
  using base_t::increase ;

  // we have three 3D vectors of unit length. The are orthogonal and
  // can be produced e.g. by applying a rotational quaternion to the
  // three cardinal vectors (1,0,0), (0,1,0), (0,0,1)

  const crd3_t xx ; 
  const crd3_t yy ;
  const crd3_t zz ;

  // some helper variables, used for efficiency.

  f_v x , r , y , z ;
  crd3_v yyy ; 
  stereographic_stepper ( crd3_t _xx , crd3_t _yy , crd3_t _zz ,
                    int _width ,
                    int _height ,
                    T _a0 ,
                    T _a1 ,
                    T _b0 ,
                    T _b1 ,
                    T _bias_x = 0 ,
                    T _bias_y = 0 
                  )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 ,
             _bias_x , _bias_y )
  { }

  // for stereographic images, we can't use any handy invariants
  // - each sample point has to be calculated 'all the way', and
  // we use this function for the calculation:

  void work ( crd3_v & trg )
  {
    auto a = M_PI_2 - 2.0 * atan ( norm ( planar ) / 2.0 ) ;

    // angle from the horizontal, same as the 'longitude'

    auto b = atan2 ( planar[0] , planar[1] ) ;

    // these values pertain to the vertical:

    sincos ( a , z , r ) ;

    // these values pertain to the horizontal:

    sincos ( b , x , y ) ;

    // now we can initialize the target value

    trg = xx * r * x + zz * z + yy * r * y ;
  }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and this
  // function infers the start value from it.

  void init ( crd3_v & trg , const crd_t & crd )
  {
    // the base type's init initializes 'planar'

    init ( crd ) ;
    work ( trg ) ;
  }

  // increase modifies it's argument to contain the next value

  void increase ( crd3_v & trg )
  {
    // the base class increase increases planar[0]

    increase() ;
    work ( trg ) ;
  }

  void init ( crd3_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    init ( trg , crd ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }

  void increase ( crd3_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    increase ( trg ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }
} ;

// with the new notion of the IR image of a cubemap, which gives
// a direct conversion between IR and ray coordinates both ways,
// we can now set up a cubemap stepper, which will 'cast rays'
// to populate a cubemap image. With 'tight' extent of (-1,-6),
// (1,6) we get a minimal cubemap with precisely 90 degrees
// fov - larger extent poduces images with larger fov.

template < typename T ,     // fundamental type
           std::size_t L ,  // lane count
           bool normalize = true >
struct cubemap_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
  using typename base_t::value_t ;
  using base_t::size ;
  using typename base_t::crd_t ;
  using typename base_t::crd3_t ;
  using typename base_t::crd3_v ;
  using typename base_t::f_v ;
  using base_t::planar ;
  using base_t::init ;
  using base_t::increase ;
  using base_t::width ;

  // we have three 3D vectors of unit length. The are orthogonal and
  // can be produced e.g. by applying a rotational quaternion to the
  // three cardinal vectors (1,0,0), (0,1,0), (0,0,1)

  const crd3_t xx ; 
  const crd3_t yy ;
  const crd3_t zz ;
  int face ;

  // some helper variables, used for efficiency.

  crd3_v ccc , vvv ;
  T section_md , refc_md ;

public:

  cubemap_stepper ( crd3_t _xx , crd3_t _yy , crd3_t _zz ,
                    int _width ,
                    int _height ,
                    T _a0 ,
                    T _a1 ,
                    T _b0 ,
                    T _b1 ,
                    T _bias_x = 0 ,
                    T _bias_y = 0
                  )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 ,
             _bias_x , _bias_y ) ,
    section_md ( _a1 - _a0 ) ,
    refc_md ( ( _a1 - _a0 ) / 2.0 )
  { }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and
  // this function infers the start value from it.

  void init ( crd3_v & trg , const crd_t & crd )
  {
    // the base type's init initializes 'planar'

    init ( crd ) ;

    // when traversing a line of a cubemap image, the cube face
    // remains the same. The variation is along the planar
    // horizontal, which, for the given cube face, is collinear
    // with one of the three vectors xx, yy or zz. The remaining
    // two vectors can be contracted into a constant expression.
    // Here we get maximal benefits from the invariants: once set
    // up, the code to produce the rays is as simple as for the
    // rectilinear case. The set-up is quite complex, though.

    face = crd[1] / width ;

    // find p1, the in-face coordinate derived from the planar
    // y coordinate in planar[1]

    auto p1 = planar[1] + ( 3 - face ) * section_md - refc_md ;

    // now we set up two vectors: ccc for the part which remains
    // constant throughout the segment, and vvv, which is in the
    // direction of the variable component given by planar[0]
    // the code is built along the case switch in ir_to_ray
    // (in geometry.h) - we could use this functor here, but
    // that would not exploit the invariants, which is the
    // whole point of setting up the stepper.

    switch ( face )
    {
      case CM_LEFT :
        ccc = -1.0 * xx + p1 * yy ;
        vvv = zz ;
        break ;
      case CM_RIGHT :
        ccc =  1.0 * xx + p1 * yy ;
        vvv = - zz ;
        break ;
      case CM_TOP :
        ccc = -1.0 * yy - p1 * zz ;
        vvv = - xx ;
        break ;
      case CM_BOTTOM :
        ccc =  1.0 * yy + p1 * zz ;
        vvv = - xx ;
        break ;
      case CM_FRONT :
        ccc =  p1 * yy + 1.0 * zz ;
        vvv = xx ;
        break ;
      case CM_BACK :
      default :
        ccc =  p1 * yy - 1.0 * zz ;
        vvv = - xx ;
        break ;
    }

    // now we can produce the target value with just three
    // multiplications and three additions:

    trg = ccc + ( planar[0] ) * vvv ;

    if constexpr ( normalize )
    {
      trg /= norm ( trg ) ;
    }
  }

  // increase modifies it's argument to contain the next value

  void increase ( crd3_v & trg )
  {
    // the base class increase increases planar[0]

    increase() ;

    trg = ccc + ( planar[0] ) * vvv ;

    if constexpr ( normalize )
    {
      trg /= norm ( trg ) ;
    }
  }

  // the capped variants:

  void init ( crd3_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    init ( trg , crd ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }

  void increase ( crd3_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    increase ( trg ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }
} ;

// This is a variation on cubemap_stepper, where the planar coordinate
// is modified to produce a more even sampling over the spherical
// surface. cubemaps use rectilinear projection, which produces
// denser sampling near the edges. Here, the sampling proceeds
// along equal angular steps in the horizontal and vertical,
// resulting in steps which become larger towards the edges in the
// plane (due to the tangens growing faster than the angle).
// Since this stepper uses angular steps, there should be efficient
// ways to produce a normalized output.

template < typename T ,     // fundamental type
           std::size_t L ,  // lane count
           bool normalize = true >
struct biatan6_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
  using typename base_t::value_t ;
  using base_t::size ;
  using typename base_t::crd_t ;
  using typename base_t::crd3_t ;
  using typename base_t::crd3_v ;
  using typename base_t::f_v ;
  using base_t::planar ;
  using base_t::init ;
  using base_t::increase ;
  using base_t::width ;

  // we have three 3D vectors of unit length. The are orthogonal and
  // can be produced e.g. by applying a rotational quaternion to the
  // three cardinal vectors (1,0,0), (0,1,0), (0,0,1)

  const crd3_t xx ; 
  const crd3_t yy ;
  const crd3_t zz ;
  int face ;

  // some helper variables, used for efficiency.

  crd3_v ccc , vvv ;
  T section_md , refc_md ;

public:

  biatan6_stepper ( crd3_t _xx , crd3_t _yy , crd3_t _zz ,
                       int _width ,
                       int _height ,
                       T _a0 ,
                       T _a1 ,
                       T _b0 ,
                       T _b1 ,
                       T _bias_x = 0 ,
                       T _bias_y = 0
                     )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 ,
             _bias_x , _bias_y ) ,
    section_md ( _a1 - _a0 ) ,
    refc_md ( ( _a1 - _a0 ) / 2.0 )
  { }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and
  // this function infers the start value from it.

  void init ( crd3_v & trg , const crd_t & crd )
  {
    // the base type's init initializes 'planar'

    init ( crd ) ;

    // when traversing a line of a cubemap image, the cube face
    // remains the same. The variation is along the planar
    // horizontal, which, for the given cube face, is collinear
    // with one of the three vectors xx, yy or zz. The remaining
    // two vectors can be contracted into a constant expression.
    // Here we get maximal benefits from the invariants: once set
    // up, the code to produce the rays is as simple as for the
    // rectilinear case. The set-up is quite complex, though.

    face = crd[1] / width ;

    // find p1, the in-face coordinate derived from the planar
    // y coordinate in planar[1]

    auto p1 = planar[1] + ( 3 - face ) * section_md - refc_md ;
    auto p0 = planar[0] ;

    // instead of using p0 and p1 directly (which would be right
    // for rectilinear cube faces) we scale from +/- 1 to +/- pi/4
    // and calculate the tangens. The result is again in +/- 1.
    // This is the inverse in-plane function to the one used on
    // 'biatan6'-type cubemap input, namely: atan ( p ) * 4 / pi
    // The result is a compression of the content towards the edges
    // of the cube faces (where it is normally unduly stretched
    // because of the rectilinear projection) and a widening in
    // the center, where there is now more space due to the
    // compression near the edges. Where the sample steps along
    // a horizon line (central horizontal of one of the surrounding
    // four cube faces) correspond with rays whose angular difference
    // decreases towards the edges with rectilinear projection, with
    // this in-plane transformation the corresponding rays are all
    // separated by the same angular step. There is still a distortion
    // which can't be helped (we're modelling a curved 2D manifold
    // on a plane), but it amounts to a maximum of 4/pi in the
    // center of a vertical near a horizontal edge (or vice versa)
    // due to the scaling factor used both ways.

    p1 = tan ( p1 * T ( M_PI / 4.0 ) ) ;
    p0 = tan ( p0 * T ( M_PI / 4.0 ) ) ;

    // now we set up two vectors: ccc for the part which remains
    // constant throughout the segment, and vvv, which is in the
    // direction of the variable component given by planar[0]
    // the code is built along the case switch in ir_to_ray
    // (in geometry.h) - we could use this functor here, but
    // that would not exploit the invariants, which is the
    // whole point of setting up the stepper.

    switch ( face )
    {
      case CM_LEFT :
        ccc = -1.0 * xx + p1 * yy ;
        vvv = zz ;
        break ;
      case CM_RIGHT :
        ccc =  1.0 * xx + p1 * yy ;
        vvv = - zz ;
        break ;
      case CM_TOP :
        ccc = -1.0 * yy - p1 * zz ;
        vvv = - xx ;
        break ;
      case CM_BOTTOM :
        ccc =  1.0 * yy + p1 * zz ;
        vvv = - xx ;
        break ;
      case CM_FRONT :
        ccc =  p1 * yy + 1.0 * zz ;
        vvv = xx ;
        break ;
      case CM_BACK :
      default :
        ccc =  p1 * yy - 1.0 * zz ;
        vvv = - xx ;
        break ;
    }

    // now we can produce the target value with just three
    // multiplications and three additions:

    trg = ccc + p0 * vvv ;

    if constexpr ( normalize )
    {
      trg /= norm ( trg ) ;
    }
  }

  // increase modifies it's argument to contain the next value

  void increase ( crd3_v & trg )
  {
    // the base class increase increases planar[0]

    increase() ;

    trg = ccc + ( tan ( planar[0] * T ( M_PI / 4.0 ) ) ) * vvv ;

    if constexpr ( normalize )
    {
      trg /= norm ( trg ) ;
    }
  }

  // the capped variants:

  void init ( crd3_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    init ( trg , crd ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }

  void increase ( crd3_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    increase ( trg ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }
} ;

// variant of stepper which yields the ray coordinate for
// the pick-up location itself and two more locations, which
// are one 'bias' step in x- or y-direction away. The caller
// can use these additional rays to approximate the derivatives,
// which is needed for 'twining'. The default for 'bias' is
// .25 - it should be less than .5

template < typename T ,     // fundamental type
           std::size_t L ,  // lane count
           template < typename , size_t , bool > class S ,
           bool normalize = true >
struct deriv_stepper
{
  typedef T value_t ;
  static const std::size_t size = 9 ;
  typedef zimt::xel_t < T , 3 > crd3_t ;
  typedef zimt::simdized_type < crd3_t , L > crd3_v ;
  typedef zimt::xel_t < T , 9 > crd9_t ;
  typedef zimt::simdized_type < crd9_t , L > crd9_v ;
  typedef zimt::xel_t < T , 2 > crd2_t ;
  typedef zimt::simdized_type < crd2_t , L > crd2_v ;
  typedef zimt::simdized_type < T , L > f_v ;
  typedef zimt::xel_t < int , 2 > crd_t ;
  typedef zimt::simdized_type < crd_t , L > crd_v ;
  typedef typename crd_v::value_type crd_ele_v ;

  S < T , L , normalize > r00 ; // yields current ray coordinate
  S < T , L , normalize > r10 ; // ray coordinate with x += bias
  S < T , L , normalize > r01 ; // ray coordinate with y += bias

  deriv_stepper ( crd3_t _xx , crd3_t _yy , crd3_t _zz ,
                  int _width ,
                  int _height ,
                  T _a0 ,
                  T _a1 ,
                  T _b0 ,
                  T _b1 ,
                  T bias = .25
                )
  : r00 ( _xx , _yy , _zz , _width , _height ,
          _a0 , _a1 , _b0 , _b1 , 0 , 0 ) ,
    r10 ( _xx , _yy , _zz , _width , _height ,
          _a0 , _a1 , _b0 , _b1 , bias , 0 ) ,
    r01 ( _xx , _yy , _zz , _width , _height ,
          _a0 , _a1 , _b0 , _b1 , 0 , bias )
    { }

  // c'tor overload for generic steppers. Here, we don't have the
  // three basis vectors - the rotation is incorporated in the
  // tf_t functor.

  typedef grok_type < xel_t < float , 2 > ,
                      xel_t < float , 3 > ,
                      L > tf_t ;

  deriv_stepper ( int _width ,
                  int _height ,
                  T _a0 ,
                  T _a1 ,
                  T _b0 ,
                  T _b1 ,
                  const tf_t & tf ,
                  T bias = .25
                )
  : r00 ( _width , _height ,
          _a0 , _a1 , _b0 , _b1 , 0 , 0 , tf ) ,
    r10 ( _width , _height ,
          _a0 , _a1 , _b0 , _b1 , bias , 0 , tf ) ,
    r01 ( _width , _height ,
          _a0 , _a1 , _b0 , _b1 , 0 , bias , tf )
    { }

  void init ( crd9_v & trg , const crd_t & crd )
  {
    crd3_v trg00 , trg10 , trg01 ;

    // we initialize the three sub-steppers to the current discrete
    // coordinate. Note that r10 and r01 have been set up with bias
    // values to produce slightly offsetted planar coordinates.

    r00.init ( trg00 , crd ) ;
    r10.init ( trg10 , crd ) ;
    r01.init ( trg01 , crd ) ;

    // now we move the resulting rays to the 'ninepack':

    trg[0] = trg00[0] ;
    trg[1] = trg00[1] ;
    trg[2] = trg00[2] ;
    trg[3] = trg10[0] ;
    trg[4] = trg10[1] ;
    trg[5] = trg10[2] ;
    trg[6] = trg01[0] ;
    trg[7] = trg01[1] ;
    trg[8] = trg01[2] ;
  }

  void increase ( crd9_v & trg )
  {
    crd3_v trg00 , trg10 , trg01 ;

    r00.increase ( trg00 ) ;
    r10.increase ( trg10 ) ;
    r01.increase ( trg01 ) ;

    trg[0] = trg00[0] ;
    trg[1] = trg00[1] ;
    trg[2] = trg00[2] ;
    trg[3] = trg10[0] ;
    trg[4] = trg10[1] ;
    trg[5] = trg10[2] ;
    trg[6] = trg01[0] ;
    trg[7] = trg01[1] ;
    trg[8] = trg01[2] ;
  }

  // the capped variants:

  void init ( crd9_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    init ( trg , crd ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }

  void increase ( crd9_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    increase ( trg ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }
} ;

// this 2D 'stepper' class does not actually produce 3D ray
// coordinates, but rather 2D planar coordinates. This is used if
// planar reprojections need to be applied before the 3D ray
// coordinate is formed. This class is pretty much 'just the base
// type', with the complete set of four member functions needed
// to produce a get_t type object usable by zimt::process. Note
// again that this class provides 2D coordinates!
// We can build the functionality of the other steppers by starting
// out with a planar_stepper object and 'suffixing' the formation
// of ray coordinates. This is 'just as good' arithmetically, but
// the code resulting from the use of the 3D steppers above is more
// efficient, because it is more specific and can use some invariants.

template < typename T ,     // fundamental type
           std::size_t L >  // lane count
struct planar_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
  using typename base_t::value_t ;
  using base_t::size ;
  using typename base_t::crd_t ;
  using typename base_t::crd2_t ;
  using typename base_t::crd2_v ;
  using typename base_t::f_v ;
  using base_t::planar ;
  using base_t::init ;
  using base_t::increase ;

  // without the 3D part, we only initialize the base type.

  planar_stepper ( int _width ,
                   int _height ,
                   T _a0 ,
                   T _a1 ,
                   T _b0 ,
                   T _b1 ,
                   T _bias_x = 0 ,
                   T _bias_y = 0 )
  : base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 ,
             _bias_x , _bias_y )
  { }

  void init ( crd2_v & trg , const crd_t & crd )
  {
    init ( crd ) ;
    trg = planar ;
  }

  void init ( crd2_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    init ( trg , crd ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }

  void increase ( crd2_v & trg )
  {
    increase() ;
    trg = planar ;
  }

  void increase ( crd2_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    increase ( trg ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }
} ;

END_ZIMT_SIMD_NAMESPACE
HWY_AFTER_NAMESPACE() ;

#endif // sentinel
