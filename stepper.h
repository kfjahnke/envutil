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

// This header provides a set of 'steppers': concrete implementations
// of objects which can serve as data source for zimt::process. The
// 'steppers' are quite specific: the provide 3D 'ray' coordinates
// to sample an environment so that these coordinates will pick pixels
// to populate a given 2D image in a specific projection.
// The steppers only produce the coordinates - providing pixel values
// to populate the image has to be done in a second step, with the
// 'act' functor of zimt::process.
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
// cycles. The built-in rotation used here can exploit certain
// invariants and mathematical tricks which make it more efficient.
// The intended use is for extraction of partial images from
// environments, reprojection, panorama viewing and the likes.
// A concrete program using this code is 'extract.cc' in this repo.

// I have finally grasped the principle how to avoid rotating
// coordinates at all: if one first calculates the 3D ray for the
// 'archetypal' unrotated 2D manifold, it's three components are
// used as factors for the three rotated vectors representing the
// re-oriented, rotated manifold. If no invariants apply, this
// needs nine multiplications and six additions, but many steppers
// have invariant parts throughout the processing of a segment.
// This hints at further possibilities: instead of using an
// orthonormal set of vectors, any set of three vecors might be
// used for all kinds of effects. The steppers take the vectors
// as arguments, so there is ample room for experiments. The
// steppers also take limits for the planar coordinates, to define
// the area which will actually be sampled. The limits are given
// in model space units, which requires knowledge about the extent
// of the 2D archetype in question. An alternative would be to use
// texture coordinates, but to work with texture coordinates, we'd
// have to move to model space units for the internal caclculations,
// so 'feeding' model space units is more efficient. extract.cc
// can accept a horizontal field of view on the command line and
// calculates the extent's limits from t. This should be the most
// common way to go - the extent values produced in this way are
// symmetrical and isotropic. Here, in the implementation, we
// process the extent values, because they allow more flexibility.

#include "zimt/zimt.h"
#include "geometry.h"

typedef zimt::xel_t < double , 2 > crd2_t ;
typedef zimt::xel_t < double , 3 > crd3_t ;

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
// which fits the intended use - lux has the 'extent' of 2D manifolds
// readily at hand and uses 'edge-to-edge' semantics for the extent.

template < typename T ,     // fundamental type
           std::size_t L >  // lane count
struct stepper_base
{
  typedef zimt::xel_t < T , 3 > crd3_t ;
  typedef zimt::simdized_type < crd3_t , L > crd3_v ;
  typedef zimt::xel_t < T , 2 > crd2_t ;
  typedef zimt::simdized_type < crd2_t , L > crd2_v ;
  typedef zimt::simdized_type < T , L > f_v ;
  typedef zimt::xel_t < long , 2 > crd_t ;
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

  const T fx0 , fx1 , fy0 , fy1 ;

  // this is a SIMD constant holding values 0, 2, 4, ...

  const crd_ele_v iota2 ;

  // width and height give the number of sample steps over the entire
  // surface covered by the image's extent. These two values define
  // the stepper's - and it's generated sampling's - resolution.
  // The width and height are independent of the extent and of
  // each other, so anisotropic sampling is allowed and a given
  // extent can be sampled with any horizontal and vertical
  // resolution.

  const std::size_t width , height ;

  // stepper_base's standard c'tor receives the width and height,
  // and the image's extent, given as the extrema first in
  // horizontal, then in vertical direction. Note that we use
  // edge-to-edge semantics to place the planar coordinates
  // inside the extent.

  stepper_base ( std::size_t _width , std::size_t _height ,
                 T _a0 , T _a1 , T _b0 , T _b1 )
  : width ( _width ) ,
    height ( _height) ,
    fx1 ( _a1 / ( 2.0 * _width ) ) ,
    fx0 ( _a0 / ( 2.0 * _width ) ) ,
    fy1 ( _b1 / ( 2.0 * _height ) ) ,
    fy0 ( _b0 / ( 2.0 * _height ) ) ,
    delta ( L * ( _a1 - _a0 ) / long ( _width ) ) ,
    iota2 ( 2L * crd_ele_v::iota() )
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
    auto ll0 = iota2 + ( crd[0] * 2L + 1L ) ;

    planar[0] = ll0 * fx1 + ( 2L * width - ll0 ) * fx0 ;

    auto ll1 = crd[1] * 2L + 1L ;

    planar[1] = ll1 * fy1 + ( 2L * height - ll1 ) * fy0 ;
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
           std::size_t L >  // lane count
class spherical_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
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
                      std::size_t _width ,
                      std::size_t _height ,
                      T _a0 = - M_PI ,
                      T _a1 = M_PI ,
                      T _b0 = - M_PI_2 ,
                      T _b1 = M_PI_2 ,
                      int _dy = 0 )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 )
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

template < typename T ,     // fundamental type
           std::size_t L >  // lane count
class cylindrical_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
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

public:

  cylindrical_stepper ( crd3_t _xx , crd3_t _yy , crd3_t _zz ,
                        std::size_t _width ,
                        std::size_t _height ,
                        T _a0 ,
                        T _a1 ,
                        T _b0 ,
                        T _b1 ,
                        int _dy = 0 
                       )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 )
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
// the planar coordinate can be mapped more or less diretcly
// to the output.

template < typename T ,     // fundamental type
           std::size_t L >  // lane count
class rectilinear_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
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
                        std::size_t _width ,
                        std::size_t _height ,
                        T _a0 ,
                        T _a1 ,
                        T _b0 ,
                        T _b1 ,
                        int _dy = 0 
                       )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 )
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
  }

  // increase modifies it's argument to contain the next value

  void increase ( crd3_v & trg )
  {
    // the base class increase increases planar[0]

    increase() ;

    // calculate the new target value

    trg = xx * planar[0] + ddd ;
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

// fisheye_stepper is a bit more involved mathematically, and
// there are no handy invariances, because the fisheye projection
// is a radial function, whereas the traversal of the 2D manifold
// is a tensor, so there is no common mathematical ground and
// each sample point (or, rather, each batch of sample points)
// has to be calculated 'all the way' from the planar coordinate.

template < typename T ,     // fundamental type
           std::size_t L >  // lane count
struct fisheye_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
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
                    std::size_t _width ,
                    std::size_t _height ,
                    T _a0 ,
                    T _a1 ,
                    T _b0 ,
                    T _b1 ,
                    int _dy = 0 
                  )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 )
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
           std::size_t L >  // lane count
struct stereographic_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
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
                    std::size_t _width ,
                    std::size_t _height ,
                    T _a0 ,
                    T _a1 ,
                    T _b0 ,
                    T _b1 ,
                    int _dy = 0 
                  )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 )
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
           std::size_t L >  // lane count
class cubemap_stepper
: public stepper_base < T , L >
{
  typedef stepper_base < T , L > base_t ;
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
  const int dy ;
  int face ;

  // some helper variables, used for efficiency.

  crd3_v ccc , vvv ;
  T section_md , refc_md ;

public:

  cubemap_stepper ( crd3_t _xx , crd3_t _yy , crd3_t _zz ,
                    std::size_t _width ,
                    std::size_t _height ,
                    T _a0 ,
                    T _a1 ,
                    T _b0 ,
                    T _b1 ,
                    int _dy = 0
                  )
  : xx ( _xx ) ,
    yy ( _yy ) ,
    zz ( _zz ) ,
    base_t ( _width , _height , _a0 , _a1 , _b0 , _b1 ) ,
    section_md ( _a1 - _a0 ) ,
    refc_md ( ( _a1 - _a0 ) / 2.0 ) ,
    dy ( _dy )
  { }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and
  // this function infers the start value from it.

  void init ( crd3_v & trg , const crd_t & crd )
  {
    // the base type's init initializes 'planar'

    init ( crd ) ;

    // when using derivative_stepper, two additional steppers are
    // employed: one where the coordinate is shifted by one to the
    // right and one where it's shifted by one downwards. In the
    // cubemap stepper, this would result in the latter value being
    // calculated with a different face index if
    // crd[1] % ( width - 1 ) is zero. The 'dy' term is to counteract
    // this problem. It's passed to the c'tor and the face index is
    // calculated with the correct crd[1].

    face = ( crd[1] + dy ) / width ;
    assert ( face >= 0 && face <= 5 ) ;

    // when traversing a line of a cubemap image, the cube face
    // remains the same. The variation is along the planar
    // horizontal, which, for the given cube face, is collinear
    // with one of the three vectors xx, yy or zz. The remaining
    // two vectors can be contracted into a constant expression.
    // Here we get maximal benefits from the invariants: once set
    // up, the code to produce the rays is as simple as for the
    // rectilinear case. The set-up is quite complex, though.

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
  }

  // increase modifies it's argument to contain the next value

  void increase ( crd3_v & trg )
  {
    // the base class increase increases planar[0]

    increase() ;

    trg = ccc + ( planar[0] ) * vvv ;
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

// variant of stepper which yields the coordiate and the
// derivatives for a a canonical step in x and y direction

template < typename T ,     // fundamental type
           std::size_t L ,  // lane count
           template < typename , size_t > class S >
struct deriv_stepper
{
  typedef zimt::xel_t < T , 3 > crd3_t ;
  typedef zimt::simdized_type < crd3_t , L > crd3_v ;
  typedef zimt::xel_t < T , 9 > crd9_t ;
  typedef zimt::simdized_type < crd9_t , L > crd9_v ;
  typedef zimt::xel_t < T , 2 > crd2_t ;
  typedef zimt::simdized_type < crd2_t , L > crd2_v ;
  typedef zimt::simdized_type < T , L > f_v ;
  typedef zimt::xel_t < long , 2 > crd_t ;
  typedef zimt::simdized_type < crd_t , L > crd_v ;
  typedef typename crd_v::value_type crd_ele_v ;

  S < T , L > r00 ; // yields current ray coordinate
  S < T , L > r10 ; // yields ray coordinates with x += 1
  S < T , L > r01 ; // yields ray coordinates with y -= 1

  deriv_stepper ( crd3_t _xx , crd3_t _yy , crd3_t _zz ,
                  std::size_t _width ,
                  std::size_t _height ,
                  T _a0 ,
                  T _a1 ,
                  T _b0 ,
                  T _b1
                )
  : r00 ( _xx , _yy , _zz , _width , _height ,
          _a0 , _a1 , _b0 , _b1 ) ,
    r10 ( _xx , _yy , _zz , _width , _height ,
          _a0 , _a1 , _b0 , _b1 ) ,

    // note that we pass a delta of -1. This is needed for the cubemap
    // stepper only and has no effect for other steppers, the parameter
    // is in all stepper signatures to make them compatible. The delta
    // gives the stepper a way to figure out that it's calculating
    // rays for coordinates which were offsetted by one canonical step
    // in the vertical, which is essential for the correct function of
    // the cubemap stepper (to use the same face index for the offsetted
    // calculation). Other steppers which derive subimages from the
    // coordinate would follow this pattern (e.g. dual fisheye images
    // stacked vertically)

    r01 ( _xx , _yy , _zz , _width , _height ,
          _a0 , _a1 , _b0 , _b1 , -1 )
    { }

  void init ( crd9_v & trg , const crd_t & crd )
  {
    crd3_v trg00 , trg10 , trg01 ;
    r00.init ( trg00 , crd ) ;
    r10.init ( trg10 , { crd[0] + 1 , crd[1] } ) ;
    r01.init ( trg01 , { crd[0] , crd[1] + 1 } ) ;
    trg[0] = trg00[0] ;
    trg[1] = -trg00[1] ;
    trg[2] = -trg00[2] ;
    trg[3] = trg10[0] - trg00[0] ;
    trg[4] = - ( trg10[1] - trg00[1] ) ;
    trg[5] = - ( trg10[2] - trg00[2] ) ;
    trg[6] = trg01[0] - trg00[0] ;
    trg[7] = - ( trg01[1] - trg00[1] ) ;
    trg[8] = - ( trg01[2] - trg00[2] ) ;
  }

  void increase ( crd9_v & trg )
  {
    crd3_v trg00 , trg10 , trg01 ;

    r00.increase ( trg00 ) ;
    r10.increase ( trg10 ) ;
    r01.increase ( trg01 ) ;

    trg[0] = trg00[0] ;
    trg[1] = -trg00[1] ;
    trg[2] = -trg00[2] ;
    trg[3] = trg10[0] - trg00[0] ;
    trg[4] = - ( trg10[1] - trg00[1] ) ;
    trg[5] = - ( trg10[2] - trg00[2] ) ;
    trg[6] = trg01[0] - trg00[0] ;
    trg[7] = - ( trg01[1] - trg00[1] ) ;
    trg[8] = - ( trg01[2] - trg00[2] ) ;
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
