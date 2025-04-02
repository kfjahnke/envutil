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

// This header has all the geometrical transformations used to convert
// 2D coordinates in several projections to 3D 'ray' coordinates and
// back. To make this code as versatile as possible, the functors for
// the transformations also provide a scalar 'eval' variant - when used
// with zimt::process, this is not used, but it's nice to have.
// To test the code, there is a simple test program 'geometry.cc' which
// uses all the functors in both scalar and simdized form. The actual
// code for the transformations is largely from lux, except for the
// cubemap code which is largely from the envutil project itself.

#include "envutil_basic.h"
#include "zimt/zimt.h"

#if defined(ENVUTIL_GEOMETRY_H) == defined(HWY_TARGET_TOGGLE)
  #ifdef ENVUTIL_GEOMETRY_H
    #undef ENVUTIL_GEOMETRY_H
  #else
    #define ENVUTIL_GEOMETRY_H
  #endif

HWY_BEFORE_NAMESPACE() ;
BEGIN_ZIMT_SIMD_NAMESPACE(project)

// zimt types for 2D and 3D coordinates and pixels

typedef zimt::xel_t < int , 2 > v2i_t ;
typedef zimt::xel_t < int , 3 > v3i_t ;
typedef zimt::xel_t < int , 2 > index_type ;
typedef zimt::xel_t < int , 2 > shape_type ;

typedef zimt::xel_t < float , 2 > v2_t ;
typedef zimt::xel_t < float , 3 > v3_t ;

// our 'own brand' of rotation matrix

template < typename T = double >
using r3_t = zimt::xel_t < xel_t < T , 3 > , 3 > ;

// rotate a 3D vector 'lhs' with the rotation matrix 'rhs'

template < typename T = double , typename U = double >
xel_t<U,3> rotate ( const xel_t<U,3> & lhs ,
                    const r3_t<T> & rhs )
{
  return { lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2] } ;
}

// rotate a rotation matrix 'lhs' with the rotation matrix 'rhs'
// this is a concatenation of the two rotations.

template < typename T = double >
r3_t<T> rotate ( const r3_t<T> & lhs , const r3_t<T> & rhs )
{
  return { rotate ( lhs[0] , rhs ) ,
           rotate ( lhs[1] , rhs ) ,
           rotate ( lhs[2] , rhs ) } ;
}

// transpose the rotation matrix. rotating with this matrix undoes
// a previous rotation with r. Note that transposition swaps rows
// and columns.

template < typename T = float >
r3_t<T> transpose ( const r3_t<T> & r )
{
  return { xel_t<T,3> ( { r[0][0] , r[1][0] , r[2][0] } ) ,
           xel_t<T,3> ( { r[0][1] , r[1][1] , r[2][1] } ) ,
           xel_t<T,3> ( { r[0][2] , r[1][2] , r[2][2] } ) } ;
}

template < typename T , std::size_t L >
struct rotate_t
: public unary_functor < xel_t<T,3> , xel_t<T,3> ,  L >
{
  const r3_t < T > r3 ;

  rotate_t ( const r3_t < T > & _r3 )
  : r3 ( _r3 )
  { }

  template < typename I , typename O >
  void eval ( const I & in , O & out )
  {
    out = rotate ( in , r3 ) ;
  }
} ;

// some SIMDized types we'll use. I use 16 SIMD lanes for now,
// which is also the lane count currently supported by OIIO.

#define LANES 16

typedef zimt::simdized_type < float , LANES > f_v ;
typedef zimt::simdized_type < int , LANES > i_v ;
typedef zimt::simdized_type < v2_t ,  LANES > crd2_v ;
typedef zimt::simdized_type < v2i_t , LANES > v2i_v ;
typedef zimt::simdized_type < v3i_t , LANES > v3i_v ;
typedef zimt::simdized_type < v3_t ,  LANES > crd3_v ;
typedef zimt::simdized_type < index_type , LANES > index_v ;

// coordinate transformations, coded as templates in zimt 'act'
// functor style, returning the result via a reference argument

// code to convert lat/lon to 3D ray. Note that I am using lux
// coordinate convention ('book order': x is right, y down and
// z forward). The lat/lon values coming in are angles in radians,
// and the resulting 'ray' coordinates are in 'model space' units.
// We're using an implicit radius of 1.0 for the conversion from
// spherical to cartesian coordinates, so the output has unit length.

template < typename T = float , std::size_t L = LANES >
struct ll_to_ray_t
: public zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 3 > , L >
{
  typedef zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 3 > , L > base_t ;

  using typename base_t::in_type ;
  using typename base_t::in_ele_v ;
  using typename base_t::in_v ;

  using typename base_t::out_type ;
  using typename base_t::out_v ;

  void eval ( const in_v & in ,
              out_v & out ) const
  {
    // incoming, we have lon/lat coordinates.

    auto const & lon ( in[0] ) ;
    auto const & lat ( in[1] ) ;

    // outgoing, we have a directional vector.

    auto & right   ( out [ RIGHT   ] ) ;
    auto & down    ( out [ DOWN    ] ) ;
    auto & forward ( out [ FORWARD ] ) ;

    // we measure angles so that a view to directly ahead (0,0,1)
    // corresponds to latitude and longitude zero. longitude increases
    // into positive values when the view moves toward the right and
    // latitude increases into positive values when the view moves
    // downwards from the view straight ahead.
    // The code benefits from using sincos, where available. Back-ends
    // lacking direct sincos support will map to use of sin and cos.

    in_ele_v sinlat , coslat , sinlon , coslon ;
    sincos ( lat , sinlat , coslat ) ;
    sincos ( lon , sinlon , coslon ) ;

    // the x component, pointing to the right in lux, is zero at
    // longitude zero, which is affected by the sine term. The
    // cosine term affects a scaling of this 'raw' value which
    // is one for latitude zero and decreases both ways.

    right = sinlon * coslat ;

    // The z component, pointing forward, is one at longitude and
    // latitude zero, and decreases with both increasing and decreasing
    // longitude and latitude.

    forward = coslon * coslat ;

    // The y component, pointing down, is zero for the view straight
    // ahead and increases with the latitude. For latitudes above
    // the equator, we'll see negative values, and positive values
    // for views into the 'southern hemisphere'. This component is not
    // affected by the longitude.

    down = sinlat ;
  }

  void eval ( const in_type & in ,
              out_type & out ) const
  {
    // incoming, we have lon/lat coordinates.

    auto const & lon ( in[0] ) ;
    auto const & lat ( in[1] ) ;

    // outgoing, we have a directional vector.

    auto & right   ( out [ RIGHT   ] ) ;
    auto & down    ( out [ DOWN    ] ) ;
    auto & forward ( out [ FORWARD ] ) ;

    // we measure angles so that a view to directly ahead (0,0,1)
    // corresponds to latitude and longitude zero. longitude increases
    // into positive values when the view moves toward the right and
    // latitude increases into positive values when the view moves
    // downwards from the view straight ahead.
    // The code benefits from using sincos, where available. Back-ends
    // lacking direct sincos support will map to use of sin and cos.

    T sinlat , coslat , sinlon , coslon ;
    sincos ( lat , &sinlat , &coslat ) ;
    sincos ( lon , &sinlon , &coslon ) ;

    // the x component, pointing to the right in lux, is zero at
    // longitude zero, which is affected by the sine term. The
    // cosine term affects a scaling of this 'raw' value which
    // is one for latitude zero and decreases both ways.

    right = sinlon * coslat ;

    // The z component, pointing forward, is one at longitude and
    // latitude zero, and decreases with both increasing and decreasing
    // longitude and latitude.

    forward = coslon * coslat ;

    // The y component, pointing down, is zero for the view straight
    // ahead and increases with the latitude. For latitudes above
    // the equator, we'll see negative values, and positive values
    // for views into the 'southern hemisphere'. This component is not
    // affected by the longitude.

    down = sinlat ;
  }
} ;

// code to move from 3D ray coordinates to lat/lon. This is the
// reverse operation to ll_to_ray above and follows the same
// conventions. Incoming, we have 3D ray coordinates, and
// outgoing, 2D lon/lat coordinates. Note how the use of atan2
// allow us to take rays in any scale.
// The output longitude is in [-pi, pi], and it's zero for the view
// 'straight ahead' in lux convention, which coincides with the
// center of the full spherical image. Values increase as the ray
// proceeds from left to right; the wrap-around point is on the
// 'back' axis.
// The output latitude is in [-pi/2,pi/2]. Negative values pertain
// to the upper hemisphere, the values increase as the ray proceeds
// from zenith to nadir.

template < typename T = float , std::size_t L = LANES >
struct ray_to_ll_t
: public zimt::unary_functor
    < zimt::xel_t < T , 3 > , zimt::xel_t < T , 2 > , L >
{
  template < typename in_type , typename out_type >
  static void eval ( const in_type & in ,
                     out_type & out )
  {
    // incoming, we have a 3D directional vector
    
    auto const & right ( in[RIGHT] ) ;
    auto const & down ( in[DOWN] ) ;
    auto const & forward ( in[FORWARD] ) ;

    // outgoing, we have a 2D lat/lon coordinate.

    auto & lon ( out[0] ) ;
    auto & lat ( out[1] ) ;

    auto s = sqrt ( right * right + forward * forward ) ;
    lat = atan2 ( down , s ) ;
    lon = atan2 ( right , forward ) ;
  }
} ;

template < typename T = float , std::size_t L = LANES >
struct ll_to_px_t
: public zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 2 > , L >
{
  T scale ;
  ll_to_px_t ( std::size_t _height )
  : scale (  T ( _height ) / M_PI )
  { }

  template < typename in_type , typename out_type >
  void eval ( const in_type & in ,
              out_type & out ) const
  {
    const zimt::xel_t < T , 2 > shift { M_PI , M_PI_2 } ;
    out = ( in + shift ) * scale ;
    out -= .5f ;
  }
} ;

// more transformations from 3D ray coordinates to various
// 2D projections: rectilinear, cylindrical, stereographic
// and fisheye

template < typename T = float , std::size_t L = LANES >
struct ray_to_rect_t
: public zimt::unary_functor
    < zimt::xel_t < T , 3 > , zimt::xel_t < T , 2 > , L >
{
  template < typename in_type , typename out_type >
  static void eval ( const in_type & in ,
                     out_type & out )
  {
    // incoming, we have a 3D directional vector
    
    auto const & right ( in[RIGHT] ) ;
    auto const & down ( in[DOWN] ) ;
    auto const & forward ( in[FORWARD] ) ;

    // outgoing, we have a 2D rectilinear coordinate.

    auto & horizontal ( out[0] ) ;
    auto & vertical ( out[1] ) ;

    // dividing by the z coordinate projects the rays to a plane
    // at unit distance. Note that we don't handle unsuitable input:
    // forward == 0 will result in Inf values, and rays pointing
    // to the 'back' hemisphere (negative z coordinate) will
    // also produce output.

    horizontal = right / forward ;
    vertical = down / forward ;
  }
} ;

// reverse operation, mapping a point on a plane to a 3D ray
// coordinate. The plane is taken to be at unit distance forward.
// The resulting 3D ray coordinate is not normalized.

template < typename T = float , std::size_t L = LANES >
struct rect_to_ray_t
: public zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 3 > , L >
{
  template < typename in_type , typename out_type >
  static void eval ( const in_type & in ,
                     out_type & out )
  {
    // incoming, we have a 2D planar coordinate
    
    auto const & horizontal ( in[0] ) ;
    auto const & vertical ( in[1] ) ;

    // outgoing, we have a 3D ray coordinate.

    auto & right ( out[RIGHT] ) ;
    auto & down ( out[DOWN] ) ;
    auto & forward ( out[FORWARD] ) ;

    right = horizontal ;
    down = vertical ;
    forward = 1 ;
  }
} ;

template < typename T = float , std::size_t L = LANES >
struct ray_to_cyl_t
: public zimt::unary_functor
    < zimt::xel_t < T , 3 > , zimt::xel_t < T , 2 > , L >
{
  template < typename in_type , typename out_type >
  static void eval ( const in_type & in ,
                     out_type & out )
  {
    // incoming, we have a 3D directional vector
    
    auto const & right ( in[RIGHT] ) ;
    auto const & down ( in[DOWN] ) ;
    auto const & forward ( in[FORWARD] ) ;

    // outgoing, we have a 2D lat/lon coordinate.

    auto & lon ( out[0] ) ;
    auto & lat ( out[1] ) ;

    auto s = sqrt ( right * right + forward * forward ) ;
    lat = down / s ;
    lon = atan2 ( right , forward ) ;
  }
} ;

// reverse operation

template < typename T = float , std::size_t L = LANES >
struct cyl_to_ray_t
: public zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 3 > , L >
{
  template < typename in_type , typename out_type >
  static void eval ( const in_type & in ,
                     out_type & out )
  {
    // incoming, we have a 2D planar coordinate
    
    auto const & horizontal ( in[0] ) ;
    auto const & vertical ( in[1] ) ;

    // outgoing, we have a 3D ray coordinate.

    auto & right ( out[RIGHT] ) ;
    auto & down ( out[DOWN] ) ;
    auto & forward ( out[FORWARD] ) ;

    // TODO: test

    forward = cos ( horizontal ) ;
    right = sin ( horizontal ) ;
    down = vertical ;
  }
} ;

template < typename T = float , std::size_t L = LANES >
struct ray_to_ster_t
: public zimt::unary_functor
    < zimt::xel_t < T , 3 > , zimt::xel_t < T , 2 > , L >
{
  template < typename in_type , typename out_type >
  static void eval ( const in_type & in ,
                     out_type & out )
  {
    auto reciprocal_norm = T ( 1 ) / sqrt (   in[0] * in[0]
                                            + in[1] * in[1]
                                            + in[2] * in[2] ) ;

    // project 3D view ray to unit sphere surface by applying the norm

    auto right = in[RIGHT] * reciprocal_norm ;
    auto down = in[DOWN] * reciprocal_norm ;
    auto forward = in[FORWARD] * reciprocal_norm ;

    // 'factor' projects x and y to stereographic: x+1 puts us to the point
    // on the sphere opposite the center of the view, and the 2.0 accounts for
    // the fact that the plane is now 2u distant instead of just 1 when seen
    // from the origin. If x gets very close to -1, we produce FLT_MAX as the
    // result, which shoud be outside the valid range

    auto factor = T ( 2 ) / ( forward + T ( 1 ) ) ;

    out[0] = right * factor ;
    // TODO: out[0] ( forward <= T ( -1 ) + FLT_EPSILON ) = FLT_MAX ;
    out[1] = down * factor ;
    // TODO: out[1] ( forward <= T ( -1 ) + FLT_EPSILON ) = FLT_MAX ;
  }
} ;

// reverse operation

template < typename T = float , std::size_t L = LANES >
struct ster_to_ray_t
: public zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 3 > , L >
{
  template < typename in_type , typename out_type >
  static void eval ( const in_type & in ,
                     out_type & out )
  {
    // incoming, we have a 2D planar coordinate
    
    auto const & horizontal ( in[0] ) ;
    auto const & vertical ( in[1] ) ;

    // outgoing, we have a 3D ray coordinate.

    auto & right ( out[RIGHT] ) ;
    auto & down ( out[DOWN] ) ;
    auto & forward ( out[FORWARD] ) ;

    // TODO: test

    auto r = sqrt ( horizontal * horizontal + vertical * vertical ) ;
    auto theta = atan ( r / 2.0 ) * 2.0f ;
    auto phi = atan2 ( horizontal , - vertical ) ;

    forward = cos ( theta ) ;
    down = - sin ( theta ) * cos ( phi ) ;
    right = sin ( theta ) * sin ( phi ) ;
  }
} ;

template < typename T = float , std::size_t L = LANES >
struct ray_to_fish_t
: public zimt::unary_functor
    < zimt::xel_t < T , 3 > , zimt::xel_t < T , 2 > , L >
{
  template < typename in_type , typename out_type >
  static void eval ( const in_type & in ,
                     out_type & out )
  {
    auto const & right ( in[RIGHT] ) ;
    auto const & down ( in[DOWN] ) ;
    auto const & forward ( in[FORWARD] ) ;

    auto s = sqrt ( right * right + down * down ) ;

    auto r = T ( M_PI_2 ) - atan2 ( forward , s ) ;

    auto phi = atan2 ( down , right ) ;

    out[0] = r * cos ( phi ) ;
    out[1] = r * sin ( phi ) ;
  }
} ;

// reverse operation

template < typename T = float , std::size_t L = LANES >
struct fish_to_ray_t
: public zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 3 > , L >
{
  template < typename in_type , typename out_type >
  static void eval ( const in_type & in ,
                     out_type & out )
  {
    // incoming, we have a 2D planar coordinate
    
    auto const & horizontal ( in[0] ) ;
    auto const & vertical ( in[1] ) ;

    // outgoing, we have a 3D ray coordinate.

    auto & right ( out[RIGHT] ) ;
    auto & down ( out[DOWN] ) ;
    auto & forward ( out[FORWARD] ) ;

    // TODO: test

    auto r = sqrt ( horizontal * horizontal + vertical * vertical ) ;
    auto phi = atan2 ( horizontal , - vertical ) ;

    forward = cos ( r ) ;
    down = - sin ( r ) * cos ( phi ) ;
    right = sin ( r ) * sin ( phi ) ;
  }
} ;

// this functor template converts incoming in-face coordinates
// to ray coordinates for a given face index, which is passed
// as a template argument - so the sixfold 'if constexpr ...' is
// not a conditional, it's just a handy way of putting the code
// into a single function without having to write partial template
// specializations for the six possible face indices.
// currently unused.

template < face_index_t F  , typename T = float , std::size_t L = LANES >
struct in_face_to_ray
: public zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 3 > , L >
{
  // incoming, we have a 2D in-face coordinate in model space
  // units. For a cubemap with no additional support, the in-face
  // coordinates are in the range [-1,1]. Maybe the naming is
  // slightly misleading, because the incoming coordinates may
  // well come from outside the 'cube face proper', but they
  // are on the given cube face's plane. The template
  // argument F fixes the face for which the calculation is
  // performed. To create the ray, a third component is needed,
  // which is added with untit magnitude. The resulting 3D ray
  // coordinate is not normalized.

  template < typename I , typename O >
  static void eval ( const I & crd2 , O & crd3 )
  {
    if constexpr ( F == CM_FRONT )
    {
      crd3[RIGHT]   =   crd2[RIGHT] ;
      crd3[DOWN]    =   crd2[DOWN]  ;
      crd3[FORWARD] =   1.0f ;
    }
    else if constexpr ( F == CM_BACK )
    {
      crd3[RIGHT]   = - crd2[RIGHT] ;
      crd3[DOWN]    =   crd2[DOWN]  ;
      crd3[FORWARD] = - 1.0f ;
    }
    else if constexpr ( F == CM_RIGHT )
    {
      crd3[RIGHT] =     1.0f ;
      crd3[DOWN] =      crd2[DOWN]  ;
      crd3[FORWARD] = - crd2[RIGHT] ;
    }
    else if constexpr ( F == CM_LEFT )
    {
      crd3[RIGHT] =   - 1.0f ;
      crd3[DOWN] =      crd2[DOWN]  ;
      crd3[FORWARD] =   crd2[RIGHT] ;
    }

    // for bottom and top, note that we're using openEXR convention.
    // to use lux convention, invert the signs.

    else if constexpr ( F == CM_BOTTOM )
    {
      crd3[RIGHT] =   - crd2[RIGHT] ;
      crd3[DOWN] =      1.0f ;
      crd3[FORWARD] =   crd2[DOWN]  ;
    }
    else if constexpr ( F == CM_TOP )
    {
      crd3[RIGHT] =   - crd2[RIGHT] ;
      crd3[DOWN] =    - 1.0f ;
      crd3[FORWARD] = - crd2[DOWN]  ;
    }
  }
} ;

// this functor converts incoming 2D coordinates pertaining
// to the entire IR image to 3D ray coordinates in lux convention.
// This functor can serve to populate the IR image: set up a
// functor yielding model space coordinates pertaining to pixels
// in the IR image, pass these model space coordinates to this
// functor, receive ray coordinates, then glean pixel values
// for the given ray by evaluating some functor taking ray
// coordinates and yielding pixels.
// The c'tor takes two arguments: first the 'section size':
// this is equal to the IR image's width, expressed in model
// space units. If the IR image does not have additional support
// and holds cube face images of precisely ninety degrees fov,
// the value would be 2.0 precisely. With added support, it's
// slightly larger. The second argument is the distance, in model
// space units, from the left margin of a section to the cube
// face image's center. If the cube face image has even width,
// this is precisely half the section size, but with odd width,
// this isn't possible, hence the extra argument.
// To accept UL-base coordinates, instantiate with
// use_centered_coordinates false.

template < typename T = float ,
           std::size_t L = LANES ,
           bool use_centered_coordinates = true >
struct ir_to_ray_t
: public zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 3 > , L >
{
  typedef zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 3 > , L > base_t ;

  using typename base_t::in_type ;
  using typename base_t::in_ele_v ;
  using typename base_t::in_v ;

  using typename base_t::out_type ;
  using typename base_t::out_v ;

  const double section_md ;
  const double refc_md ;
  const zimt::xel_t < T , 2 > ul2c ;

  ir_to_ray_t ( double _section_md = 2.0 , double _refc_md = 1.0 )
  : section_md ( _section_md ) ,
    refc_md ( _refc_md ) ,
    ul2c { _refc_md , 3 * _section_md }
  { }

  // incoming, we have 2D model space coordinates.
  // This functor takes centered planar cordinates - note the
  // addition of ul2c at the beginning of the eval functions, which
  // converts the coordinates from to center-based to ul-based for
  // processing.

  void eval ( const in_v & _crd2 , out_v & crd3 ) const
  {
    in_v crd2 ( _crd2 ) ;

    if constexpr ( use_centered_coordinates )
      crd2 += ul2c ;

    // The numerical constants for the cube faces/sections are set
    // up so that a simple division of the y coordinate yields the
    // corresponding section index.

    i_v section ( crd2[1] / section_md ) ;

    // The incoming coordinates are relative to the upper left
    // corner of the IR image. Now we move to in-face coordinates,
    // which are centered on the cube face we're dealing with.

    crd2[1] -= section * section_md ;
    crd2    -= refc_md ;

    // the section number can also yield the 'dominant' axis
    // by dividing the value by two (another property which is
    // deliberate):

    i_v dom ( section >> 1 ) ;

    // again we use a conditional to avoid lengthy calculations
    // when there aren't any populated lanes for the given predicate

    if ( any_of ( dom == 0 ) )
    {
      auto m = ( section == CM_RIGHT ) ;
      if ( any_of ( m ) )
      {
        crd3[RIGHT](m) =     1.0f ;
        crd3[DOWN](m) =      crd2[DOWN] ;
        crd3[FORWARD](m) = - crd2[RIGHT] ;
      }
      m = ( section == CM_LEFT ) ;
      if ( any_of ( m ) )
      {
        crd3[RIGHT](m) =   - 1.0f ;
        crd3[DOWN](m) =      crd2[DOWN] ;
        crd3[FORWARD](m) =   crd2[RIGHT] ;
      }
    }
    if ( any_of ( dom == 1 ) )
    {
      auto m = ( section == CM_BOTTOM ) ;
      if ( any_of ( m ) )
      {
        crd3[RIGHT](m) =   - crd2[RIGHT] ;
        crd3[DOWN](m) =      1.0f ;
        crd3[FORWARD](m) =   crd2[DOWN] ;
      }
      m = ( section == CM_TOP ) ;
      if ( any_of ( m ) )
      {
        crd3[RIGHT](m) =   - crd2[RIGHT] ;
        crd3[DOWN](m) =    - 1.0f ;
        crd3[FORWARD](m) = - crd2[DOWN] ;
      }
    }
    if ( any_of ( dom == 2 ) )
    {
      auto m = ( section == CM_FRONT ) ;
      if ( any_of ( m ) )
      {
        crd3[RIGHT](m)   =   crd2[RIGHT] ;
        crd3[DOWN](m)    =   crd2[DOWN] ;
        crd3[FORWARD](m) =   1.0f ;
      }
      m = ( section == CM_BACK ) ;
      if ( any_of ( m ) )
      {
        crd3[RIGHT](m)   = - crd2[RIGHT] ;
        crd3[DOWN](m)    =   crd2[DOWN] ;
        crd3[FORWARD](m) = - 1.0f ;
      }
    }
  }

  // for completeness' sake, the scalar eval

  void eval ( const in_type & _crd2 , out_type & crd3 ) const
  {
    in_type crd2 ( _crd2 ) ;

    if constexpr ( use_centered_coordinates )
      crd2 += ul2c ;

    // The numerical constants for the cube faces/sections are set
    // up so that a simple division of the y coordinate yields the
    // corresponding section index.

    int section ( crd2[1] / section_md ) ;

    // The incoming coordinates are relative to the upper left
    // corner of the IR image. Now we move to in-face coordinates,
    // which are centered on the cube face we're dealing with.

    crd2[1] -= section * section_md ;
    crd2    -= refc_md ;

    // the section number can also yield the 'dominant' axis
    // by dividing the value by two (another property which is
    // deliberate):

    int dom ( section >> 1 ) ;

    // again we use a conditional to avoid lengthy calculations
    // when there aren't any populated lanes for the given predicate

    if ( dom == 0 )
    {
      if ( section == CM_RIGHT )
      {
        crd3[RIGHT] =     1.0f ;
        crd3[DOWN] =      crd2[DOWN] ;
        crd3[FORWARD] = - crd2[RIGHT] ;
      }
      else
      {
        crd3[RIGHT] =   - 1.0f ;
        crd3[DOWN] =      crd2[DOWN] ;
        crd3[FORWARD] =   crd2[RIGHT] ;
      }
    }
    else if ( dom == 1 )
    {
      if ( section == CM_BOTTOM )
      {
        crd3[RIGHT] =   - crd2[RIGHT] ;
        crd3[DOWN] =      1.0f ;
        crd3[FORWARD] =   crd2[DOWN] ;
      }
      else
      {
        crd3[RIGHT] =   - crd2[RIGHT] ;
        crd3[DOWN] =    - 1.0f ;
        crd3[FORWARD] = - crd2[DOWN] ;
      }
    }
    else
    {
      if ( section == CM_FRONT )
      {
        crd3[RIGHT]   =   crd2[RIGHT] ;
        crd3[DOWN]    =   crd2[DOWN] ;
        crd3[FORWARD] =   1.0f ;
      }
      else
      {
        crd3[RIGHT]   = - crd2[RIGHT] ;
        crd3[DOWN]    =   crd2[DOWN] ;
        crd3[FORWARD] = - 1.0f ;
      }
    }
  }
} ;

// same for biatan6 - done by copy-and-paste. the only difference
// to the class above is the in-plane transformation.

template < typename T = float ,
           std::size_t L = LANES ,
           bool use_centered_coordinates = true >
struct ba6_to_ray_t
: public zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 3 > , L >
{
  typedef zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 3 > , L > base_t ;

  using typename base_t::in_type ;
  using typename base_t::in_ele_v ;
  using typename base_t::in_v ;

  using typename base_t::out_type ;
  using typename base_t::out_v ;

  const double section_md ;
  const double refc_md ;
  const zimt::xel_t < T , 2 > ul2c ;

  ba6_to_ray_t ( double _section_md = 2.0 , double _refc_md = 1.0 )
  : section_md ( _section_md ) ,
    refc_md ( _refc_md ) ,
    ul2c { _refc_md , 3 * _section_md }
  { }

  // incoming, we have 2D model space coordinates.
  // This functor takes centered planar cordinates - note the
  // addition of ul2c at the beginning of the eval functions, which
  // converts the coordinates from to center-based to ul-based for
  // processing.

  void eval ( const in_v & _crd2 , out_v & crd3 ) const
  {
    in_v crd2 ( _crd2 ) ;

    if constexpr ( use_centered_coordinates )
      crd2 += ul2c ;

    // The numerical constants for the cube faces/sections are set
    // up so that a simple division of the y coordinate yields the
    // corresponding section index.

    i_v section ( crd2[1] / section_md ) ;

    // The incoming coordinates are relative to the upper left
    // corner of the IR image. Now we move to in-face coordinates,
    // which are centered on the cube face we're dealing with.

    crd2[1] -= section * section_md ;
    crd2    -= refc_md ;

    // apply the in-plane biatan6 transformation

    crd2 = tan ( crd2 * T ( M_PI / 4.0 ) ) ;

    // the section number can also yield the 'dominant' axis
    // by dividing the value by two (another property which is
    // deliberate):

    i_v dom ( section >> 1 ) ;

    // again we use a conditional to avoid lengthy calculations
    // when there aren't any populated lanes for the given predicate

    if ( any_of ( dom == 0 ) )
    {
      auto m = ( section == CM_RIGHT ) ;
      if ( any_of ( m ) )
      {
        crd3[RIGHT](m) =     1.0f ;
        crd3[DOWN](m) =      crd2[DOWN] ;
        crd3[FORWARD](m) = - crd2[RIGHT] ;
      }
      m = ( section == CM_LEFT ) ;
      if ( any_of ( m ) )
      {
        crd3[RIGHT](m) =   - 1.0f ;
        crd3[DOWN](m) =      crd2[DOWN] ;
        crd3[FORWARD](m) =   crd2[RIGHT] ;
      }
    }
    if ( any_of ( dom == 1 ) )
    {
      auto m = ( section == CM_BOTTOM ) ;
      if ( any_of ( m ) )
      {
        crd3[RIGHT](m) =   - crd2[RIGHT] ;
        crd3[DOWN](m) =      1.0f ;
        crd3[FORWARD](m) =   crd2[DOWN] ;
      }
      m = ( section == CM_TOP ) ;
      if ( any_of ( m ) )
      {
        crd3[RIGHT](m) =   - crd2[RIGHT] ;
        crd3[DOWN](m) =    - 1.0f ;
        crd3[FORWARD](m) = - crd2[DOWN] ;
      }
    }
    if ( any_of ( dom == 2 ) )
    {
      auto m = ( section == CM_FRONT ) ;
      if ( any_of ( m ) )
      {
        crd3[RIGHT](m)   =   crd2[RIGHT] ;
        crd3[DOWN](m)    =   crd2[DOWN] ;
        crd3[FORWARD](m) =   1.0f ;
      }
      m = ( section == CM_BACK ) ;
      if ( any_of ( m ) )
      {
        crd3[RIGHT](m)   = - crd2[RIGHT] ;
        crd3[DOWN](m)    =   crd2[DOWN] ;
        crd3[FORWARD](m) = - 1.0f ;
      }
    }
  }

  // for completeness' sake, the scalar eval

  void eval ( const in_type & _crd2 , out_type & crd3 ) const
  {
    in_type crd2 ( _crd2 ) ;

    if constexpr ( use_centered_coordinates )
      crd2 += ul2c ;

    // The numerical constants for the cube faces/sections are set
    // up so that a simple division of the y coordinate yields the
    // corresponding section index.

    int section ( crd2[1] / section_md ) ;

    // The incoming coordinates are relative to the upper left
    // corner of the IR image. Now we move to in-face coordinates,
    // which are centered on the cube face we're dealing with.

    crd2[1] -= section * section_md ;
    crd2    -= refc_md ;

    // apply the in-plane biatan6 transformation

    crd2 = tan ( crd2 * T ( M_PI / 4.0 ) ) ;

    // the section number can also yield the 'dominant' axis
    // by dividing the value by two (another property which is
    // deliberate):

    int dom ( section >> 1 ) ;

    // again we use a conditional to avoid lengthy calculations
    // when there aren't any populated lanes for the given predicate

    if ( dom == 0 )
    {
      if ( section == CM_RIGHT )
      {
        crd3[RIGHT] =     1.0f ;
        crd3[DOWN] =      crd2[DOWN] ;
        crd3[FORWARD] = - crd2[RIGHT] ;
      }
      else
      {
        crd3[RIGHT] =   - 1.0f ;
        crd3[DOWN] =      crd2[DOWN] ;
        crd3[FORWARD] =   crd2[RIGHT] ;
      }
    }
    else if ( dom == 1 )
    {
      if ( section == CM_BOTTOM )
      {
        crd3[RIGHT] =   - crd2[RIGHT] ;
        crd3[DOWN] =      1.0f ;
        crd3[FORWARD] =   crd2[DOWN] ;
      }
      else
      {
        crd3[RIGHT] =   - crd2[RIGHT] ;
        crd3[DOWN] =    - 1.0f ;
        crd3[FORWARD] = - crd2[DOWN] ;
      }
    }
    else
    {
      if ( section == CM_FRONT )
      {
        crd3[RIGHT]   =   crd2[RIGHT] ;
        crd3[DOWN]    =   crd2[DOWN] ;
        crd3[FORWARD] =   1.0f ;
      }
      else
      {
        crd3[RIGHT]   = - crd2[RIGHT] ;
        crd3[DOWN]    =   crd2[DOWN] ;
        crd3[FORWARD] = - 1.0f ;
      }
    }
  }
} ;

// to make the conversion efficient and transparent, I refrain from
// using ir_to_ray and a subsequent coordinate transformation in
// favour of this dedicated functor, which is basically a copy of
// the one above, but with different component indexes and signs
// inverted where necessary (namely for the left and up direction,
// which are the negative of zimt's right and down).

// currently unused.
/*
template < typename T = float , std::size_t L = LANES >
struct ir_to_exr
: public zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 3 > , L >
{
  const double section_md ;
  const double refc_md ;

  ir_to_exr ( double _section_md , double _refc_md )
  : section_md ( _section_md ) ,
    refc_md ( _refc_md )
  { }

  // incoming, we have 2D model space coordinates, with the origin
  // at the (total!) IR image's upper left corner.

  template < typename I , typename O >
  void eval ( const I & _crd2 , O & crd3 ) const
  {
    I crd2 ( _crd2 ) ;

    // The numerical constants for the cube faces/sections are set
    // up so that a simple division of the y coordinate yields the
    // corresponding section index.

    i_v section ( crd2[1] / section_md ) ;

    // The incoming coordinates are relative to the upper left
    // corner of the IR image. Now we move to in-face coordinates,
    // which are centered on the cube face we're dealing with.

    crd2[1] -= section * section_md ;
    crd2    -= refc_md ;

    // the section number can also yield the 'dominant' axis
    // by dividing the value by two (another property which is
    // deliberate):

    i_v dom ( section >> 1 ) ;

    // again we use a conditional to avoid lengthy calculations
    // when there aren't any populated lanes for the given predicate

    if ( any_of ( dom == 0 ) )
    {
      auto m = ( section == CM_RIGHT ) ;
      if ( any_of ( m ) )
      {
        crd3[EXR_LEFT](m) =    - 1.0f ;
        crd3[EXR_UP](m) =      - crd2[DOWN] ;
        crd3[EXR_FORWARD](m) = - crd2[RIGHT] ;
      }
      m = ( section == CM_LEFT ) ;
      if ( any_of ( m ) )
      {
        crd3[EXR_LEFT](m) =      1.0f ;
        crd3[EXR_UP](m) =      - crd2[DOWN] ;
        crd3[EXR_FORWARD](m) =   crd2[RIGHT] ;
      }
    }
    if ( any_of ( dom == 1 ) )
    {
      auto m = ( section == CM_BOTTOM ) ;
      if ( any_of ( m ) )
      {
        crd3[EXR_LEFT](m) =      crd2[RIGHT] ;
        crd3[EXR_UP](m) =      - 1.0f ;
        crd3[EXR_FORWARD](m) =   crd2[DOWN] ;
      }
      m = ( section == CM_TOP ) ;
      if ( any_of ( m ) )
      {
        crd3[EXR_LEFT](m) =      crd2[RIGHT] ;
        crd3[EXR_UP](m) =        1.0f ;
        crd3[EXR_FORWARD](m) = - crd2[DOWN] ;
      }
    }
    if ( any_of ( dom == 2 ) )
    {
      auto m = ( section == CM_FRONT ) ;
      if ( any_of ( m ) )
      {
        crd3[EXR_LEFT](m)   = - crd2[RIGHT] ;
        crd3[EXR_UP](m)    =  - crd2[DOWN] ;
        crd3[EXR_FORWARD](m) =  1.0f ;
      }
      m = ( section == CM_BACK ) ;
      if ( any_of ( m ) )
      {
        crd3[EXR_LEFT](m)   =    crd2[RIGHT] ;
        crd3[EXR_UP](m)    =   - crd2[DOWN] ;
        crd3[EXR_FORWARD](m) = - 1.0f ;
      }
    }
  }
} ;
*/

// given a 3D 'ray' coordinate, find the corresponding cube face
// and the in-face coordinate - note the two references which take
// the result values. The incoming 'ray' coordinate does not have
// to be normalized. The resulting in-face coordinates are in the
// range of [-1,1] - in 'model space units' pertaining to planes
// 'draped' at unit distance from the origin and perpendicular
// to one of the axes.
// Note how the results of this functor can be fed to one of the 
// get_pickup_coordinate_... functions provided in metrics.h to
// get 2D coordinates pertaining to a cubemap's IR representation.
// first, the simdized version:

template < typename T , std::size_t L >
void ray_to_cubeface
  ( const zimt::xel_t < zimt::simdized_type < T , L > , 3 > & c ,
          zimt::simdized_type < int , L > & face ,
          zimt::xel_t < zimt::simdized_type < T , L > , 2 > & in_face )
{
  // form three masks with relations of the numerical values of
  // the 'ray' coordinate. These are sufficient to find out which
  // component has the largest absolute value (the 'dominant' one,
  // along the 'dominant' axis)

  auto m1 = ( abs ( c[RIGHT] ) >= abs ( c[DOWN] ) ) ;
  auto m2 = ( abs ( c[RIGHT] ) >= abs ( c[FORWARD] ) ) ;
  auto m3 = ( abs ( c[DOWN] )  >= abs ( c[FORWARD] ) ) ;

  // form a mask which is true where a specific axis is 'dominant'.
  // We start out looking at the x axis: where the numerical value
  // along the x axis (pointing right) is larger than the other two,
  // 'dom' will be true.

  auto dom = m1 & m2 ;

  // Now we can assign face indexes, stored in f. If the coordinate
  // value is negative along the dominant axis, we're looking at
  // the face opposite and assign a face index one higher. Note
  // how this SIMD code does the job for LANES coordinates in
  // parallel, avoiding conditionals and using masking instead.
  // we also find in-face coordinates:
  // we divide the two non-dominant coordinate values by the
  // dominant one. One of the axes comes out just right when
  // dividing by the absolute value - e.g. the vertical axis
  // points downwards for all the four cube faces around
  // the center. The other axis is divided by the 'major'
  // coordinate value as-is; the resulting coordinate runs
  // one way for positive major values and backwards for
  // negative ones. Note that we might capture the absolute
  // values (which we've used before) in variables, but the
  // compiler will recognize the common subexpressions and
  // do it for us. While it's generally preferable to avoid
  // conditionals in inner-loop code, I use conditionals here
  // because most of the time all coordinates will 'land' in
  // the same cube face, so for two cases, the rather expensive
  // code to calculate the face index and in-face coordinate
  // can be omitted. TODO: One might test whether omitting the
  // conditionals is actually slower.

  if ( any_of ( dom ) )
  {
    // extract in-face coordinates for the right and left cube
    // face. the derivation of the x coordinate uses opposites for
    // the two faces, the direction of the y coordinate is equal.
    // Note that some lanes in c[RIGHT] may be zero and result in
    // an Inf result, but these will never be the ones which end
    // up in the result, because only those where c[RIGHT] is
    // 'dominant' will 'make it through', and where c[RIGHT] is
    // dominant, it's certainly not zero. But we rely on the
    // system not to throw a division-by-zero exception, which
    // would spoil our scheme.

    face = CM_RIGHT ;
    face ( dom & ( c[RIGHT] < 0 ) ) = CM_LEFT ;

    in_face[0] ( dom ) = - c[FORWARD] / c[RIGHT] ;
    in_face[1] ( dom ) = c[DOWN] / abs ( c[RIGHT] ) ;
  }

  // set dom true where the z axis (pointing forward) has the
  // largest numerical value
  
  dom = ( ! m2 ) & ( ! m3 ) ;

  if ( any_of ( dom ) )
  {
    // test for front and back faces

    face ( dom ) = CM_FRONT ;
    face ( dom & ( c[FORWARD] < 0 ) ) = CM_BACK ;

    in_face[0] ( dom ) = c[RIGHT] / c[FORWARD] ;
    in_face[1] ( dom ) = c[DOWN] / abs ( c[FORWARD] ) ;
  }

  // now set dom true where the y axis (pointing down) has the
  // largest numerical value

  dom = ( ! m1 ) & m3 ; 

  if ( any_of ( dom ) )
  {
    // same for the top and bottom cube faces - here the x coordinate
    // corresponds to the right 3D axis, the y coordinate depends on
    // which of the faces we're looking at (hence no abs)
    // the top and bottom images could each be oriented in four
    // different ways. The orientation I expect in this program
    // is openEXR's cubemap format, where the top and bottom
    // image align with the 'back' image (the last one down).
    // For lux conventions (the top and bottom image aligning
    // with the 'front' cube face) swap the signs in both expressions.

    face ( dom ) = CM_BOTTOM ;
    face ( dom & ( c[DOWN] < 0 ) ) = CM_TOP ;

    // lux convention:
    // in_face[0] ( dom ) =   c[RIGHT] / abs ( c[DOWN] ) ;
    // in_face[1] ( dom ) = - c[FORWARD] / c[DOWN] ;

    in_face[0] ( dom ) = - c[RIGHT] / abs ( c[DOWN] ) ;
    in_face[1] ( dom ) =   c[FORWARD] / c[DOWN] ;

  }

}

// the scalar version is much simpler:

template < typename T >
void ray_to_cubeface ( const zimt::xel_t < T , 3 > & c ,
                       int & face ,
                       zimt::xel_t < T , 2 > & in_face )
{
  auto m1 = ( std::abs ( c[RIGHT] ) >= std::abs ( c[DOWN] ) ) ;
  auto m2 = ( std::abs ( c[RIGHT] ) >= std::abs ( c[FORWARD] ) ) ;

  if ( m1 & m2 )
  {
    face = ( ( c[RIGHT] < 0 ) ? CM_LEFT : CM_RIGHT ) ;

    in_face[0] = - c[FORWARD] / c[RIGHT] ;
    in_face[1] = c[DOWN] / std::abs ( c[RIGHT] ) ;
  }
  else
  {
    auto m3 = ( std::abs ( c[DOWN] )  >= std::abs ( c[FORWARD] ) ) ;

    if ( ( ! m2 ) & ( ! m3 ) )
    {
      face = ( ( c[FORWARD] < 0 ) ? CM_BACK : CM_FRONT ) ;

      in_face[0] = c[RIGHT] / c[FORWARD] ;
      in_face[1] = c[DOWN] / std::abs ( c[FORWARD] ) ;
    }
    else
    {
      face = ( ( c[DOWN] < 0 ) ? CM_TOP : CM_BOTTOM ) ;

      in_face[0] = - c[RIGHT] / std::abs ( c[DOWN] ) ;
      in_face[1] =   c[FORWARD] / c[DOWN] ;

    }
  }

  // alternative:

  // if ( std::abs ( c[RIGHT] ) >= std::abs ( c[DOWN] ) )
  // {
  //   if ( std::abs ( c[RIGHT] ) >= std::abs ( c[FORWARD] ) )
  //   {
  //     face = ( ( c[RIGHT] < 0 ) ? CM_LEFT : CM_RIGHT ) ;
  // 
  //     in_face[0] = - c[FORWARD] / c[RIGHT] ;
  //     in_face[1] = c[DOWN] / std::abs ( c[RIGHT] ) ;
  //     return ;
  //   }
  // }
  // else
  // {
  //   if ( std::abs ( c[DOWN] )  >= std::abs ( c[FORWARD] ) )
  //   {
  //     face = ( ( c[DOWN] < 0 ) ? CM_TOP : CM_BOTTOM ) ;
  // 
  //     in_face[0] = - c[RIGHT] / std::abs ( c[DOWN] ) ;
  //     in_face[1] =   c[FORWARD] / c[DOWN] ;
  //     return ;
  //   }
  // }
  // face = ( ( c[FORWARD] < 0 ) ? CM_BACK : CM_FRONT ) ;
  // 
  // in_face[0] = c[RIGHT] / c[FORWARD] ;
  // in_face[1] = c[DOWN] / std::abs ( c[FORWARD] ) ;
}

// variant which takes a given face vector. This is used to
// approximate the first derivative, which is done by subtracting
// the result of the coordinate transformation of incoming
// coordinates which were offset by one sampling step in either
// canonical direction. We have to look at the same cube face,
// to avoid the possibility that the cube face changes between
// the calculation of the result for the actual cordinate and
// it's offsetted neighbours, which would likely result in a
// large discontinuity in the in-face coordinate, spoiling the
// result.
// So, incoming, we have a 3D ray coordinate and a set of cube
// face indices, and we'll get in-face coordinates as output.
// Note how the different use case results in the second argument
// being a const&, in contrast to the previous function, where
// it is a non-const reference - a result value.

template < typename T , std::size_t L >
void ray_to_cubeface_fixed
  ( const zimt::xel_t < zimt::simdized_type < T , L > , 3 > & c ,
    const zimt::simdized_type < int , L > & face ,
          zimt::xel_t < zimt::simdized_type < T , L > , 2 > & in_plane )
{
  // form a mask which is true where a specific axis is 'dominant'.
  // since we have the face indices already, this is simple: it's
  // the face indices shifted to the right to remove their least
  // significant bit, which codes for the sign along the dominant
  // axis. We benefit from working with a strictly codified setup;
  // the nubering of the face indices is designed specifically to
  // make such calculations as efficient as possible. The reaminder
  // of the function proceeds as above, but note how due to the
  // fact that we prescribe given target planes, the resulting
  // in-plane coordinates are no longer necessarily in (-1, 1) as
  // the in-face coordinates yielded by the previous functor.
  // We silently assume that the incoming rays are 'not very far
  // off' the face prescribed for them - e.g. not oriented so
  // that they can't intersect with the prescribed plane.

  auto dom_v = face >> 1 ;

  auto dom = ( dom_v == 0 ) ;

  if ( any_of ( dom ) )
  {
    in_plane[0] ( dom ) = - c[FORWARD] / c[RIGHT] ;
    in_plane[1] ( dom ) = c[DOWN] / abs ( c[RIGHT] ) ;
  }

  dom = ( dom_v == 1 ) ; 

  if ( any_of ( dom ) )
  {
    in_plane[0] ( dom ) = - c[RIGHT] / abs ( c[DOWN] ) ;
    in_plane[1] ( dom ) =   c[FORWARD] / c[DOWN] ;
  }

  dom = ( dom_v == 2 ) ;

  if ( any_of ( dom ) )
  {
    in_plane[0] ( dom ) = c[RIGHT] / c[FORWARD] ;
    in_plane[1] ( dom ) = c[DOWN] / abs ( c[FORWARD] ) ;
  }
}

template < typename T >
void ray_to_cubeface_fixed ( const zimt::xel_t < T , 3 > & c ,
                             const int & face ,
                             zimt::xel_t < T , 2 > & in_plane )
{
  auto dom = face >> 1 ;
  if ( dom == 0 )
  {
    in_plane[0] = - c[FORWARD] / c[RIGHT] ;
    in_plane[1] = c[DOWN] / std::abs ( c[RIGHT] ) ;
  }
  else if ( dom == 1 )
  {
    in_plane[0] = - c[RIGHT] / std::abs ( c[DOWN] ) ;
    in_plane[1] =   c[FORWARD] / c[DOWN] ;
  }
  else
  {
    in_plane[0] = c[RIGHT] / c[FORWARD] ;
    in_plane[1] = c[DOWN] / std::abs ( c[FORWARD] ) ;
  }
}

// this functor converts incoming 3D ray coordinates to 2D
// coordinates pertaining to the IR image (in model space units)
// Similar code is used in sixfold_t - here, we have a stand-alone
// version which only needs two values from the metrics_t object.
// This functor is the inversion of ir_to_ray_t, which needs the
// same two values only.
// This functor produces centered planar cordinates - note the
// subtraction of ul2c at the end of the eval functions, which
// convert the coordinates from ul-based to center-based.
// To produce UL-base coordinates, instantiate with
// use_centered_coordinates false.

template < typename T = float , std::size_t L = LANES ,
           bool use_centered_coordinates = true >
struct ray_to_ir_t
: public zimt::unary_functor
    < zimt::xel_t < T , 3 > , zimt::xel_t < T , 2 > , L >
{
  typedef zimt::unary_functor
    < zimt::xel_t < T , 3 > , zimt::xel_t < T , 2 > , L > base_t ;

  using typename base_t::in_type ;
  using typename base_t::in_ele_v ;
  using typename base_t::in_v ;

  using typename base_t::out_type ;
  using typename base_t::out_v ;

  using typename base_t::ic_v ;

  // to obtain IR coordinates, we need the width of a 'section'
  // of the IR image in model space units and the distance to
  // the center of the cube face images from the left/top margin:

  const T section_md ;
  const T refc_md ;
  const zimt::xel_t < T , 2 > ul2c ;

  ray_to_ir_t ( const T & _section_md = 2.0 ,
                const T & _refc_md = 1.0 )
  : section_md ( _section_md ) ,
    refc_md ( _refc_md ) ,
    ul2c { _refc_md , 3 * _section_md }
  { }

  void eval ( const in_type & c , out_type & ir_c ) const
  {
    // first, we glean the face indices and in-face coordinates

    int face ;

    ray_to_cubeface < T > ( c , face , ir_c ) ;

    ir_c += refc_md ;

    // for the vertical, we add the face index times the section
    // size (in model space units)

    ir_c[1] += face * section_md ;

    if constexpr ( use_centered_coordinates )
      ir_c -= ul2c ;
  }

  void eval ( const in_v & c , out_v & ir_c ) const
  {
    // first, we glean the face indices and in-face coordinates

    ic_v face ;

    ray_to_cubeface < T , L > ( c , face , ir_c ) ;

    ir_c += refc_md ;

    // for the vertical, we add the face index times the section
    // size (in model space units)

    ir_c[1] += face * section_md ;

    if constexpr ( use_centered_coordinates )
      ir_c -= ul2c ;
  }
} ;

template < typename T = float , std::size_t L = LANES ,
           bool use_centered_coordinates = true >
struct ray_to_ba6_t
: public zimt::unary_functor
    < zimt::xel_t < T , 3 > , zimt::xel_t < T , 2 > , L >
{
  typedef zimt::unary_functor
    < zimt::xel_t < T , 3 > , zimt::xel_t < T , 2 > , L > base_t ;

  using typename base_t::in_type ;
  using typename base_t::in_ele_v ;
  using typename base_t::in_v ;

  using typename base_t::out_type ;
  using typename base_t::out_v ;

  using typename base_t::ic_v ;

  // to obtain IR coordinates, we need the width of a 'section'
  // of the IR image in model space units and the distance to
  // the center of the cube face images from the left/top margin:

  const T section_md ;
  const T refc_md ;
  const zimt::xel_t < T , 2 > ul2c ;

  ray_to_ba6_t ( const T & _section_md = 2.0 ,
                 const T & _refc_md = 1.0 )
  : section_md ( _section_md ) ,
    refc_md ( _refc_md ) ,
    ul2c { _refc_md , 3 * _section_md }
  { }

  void eval ( const in_type & c , out_type & ir_c ) const
  {
    // first, we glean the face indices and in-face coordinates

    int face ;

    ray_to_cubeface < T > ( c , face , ir_c ) ;

    // perform the biatan6 in-plane transformation

    ir_c = T( 4.0 / M_PI ) * atan ( ir_c ) ;

    ir_c += refc_md ;

    // for the vertical, we add the face index times the section
    // size (in model space units)

    ir_c[1] += face * section_md ;

    if constexpr ( use_centered_coordinates )
      ir_c -= ul2c ;
  }

  void eval ( const in_v & c , out_v & ir_c ) const
  {
    // first, we glean the face indices and in-face coordinates

    ic_v face ;

    ray_to_cubeface < T , L > ( c , face , ir_c ) ;

    // perform the biatan6 in-plane transformation

    ir_c = T ( 4.0 / M_PI ) * atan ( ir_c ) ;

    ir_c += refc_md ;

    // for the vertical, we add the face index times the section
    // size (in model space units)

    ir_c[1] += face * section_md ;

    if constexpr ( use_centered_coordinates )
      ir_c -= ul2c ;
  }
} ;

// when doing geometrical transformations of images, the coordinate
// transformations can often be expressed as the chaining of two
// functors: one transforming the incoming coordinate to a 3D ray,
// and a second one transforming the 3D ray to another 2D coordinate.
// struct plane_to_plane offers just such a construct, optionally
// with a 3D-to-3D operation inserted in between, e.g. a 3D rotation.
// This is pure SIMDized code, there is no scalar overload.
// currently unused.
// This is more of a pointer to additional code: while many of the
// 2d->2D mappings have to be coded using an intermediate 3D ray
// coordinates, some mappings can do without - specifically mapping
// from one rectilinear image to another can be done without the
// 3D intermediate. The code in this repo does now rely on 'stepper'
// objects to generate 3D ray coordinates directly, which are then
// projected to a plane, so 2D->2D mappings aren't used. steppers
// have the advantage of having rotation built-in, so the resulting
// calculations are very efficient, avoiding the rotation step in
// the 3D domain, which would be necessary when using 2D->2D mappings.
// It's good, though, to have several routes to a desired end, to
// be able to double-check.
// TODO: in order to do just that (doublecheck), we should have
// the transformations from rectlinear, stereographic, cylindrical
// and fisheye coordinates to 3D ray coordinates as well. Until now,
// they were missing, because steppers take the part. I've copied
// code from lux to the effect, but it needs to be tested.

template < typename T , std::size_t L >
struct plane_to_plane
: public zimt::unary_functor
    < zimt::xel_t < T , 2 > , zimt::xel_t < T , 2 > , L >
{
  typedef zimt::simdized_type < zimt::xel_t < T , 2 > , L > crd2_v ;
  typedef zimt::simdized_type < zimt::xel_t < T , 3 > , L > crd3_v ;
  std::function < void ( const crd2_v & , crd2_v & ) > tf ;

  template < typename FA , typename FB >
  plane_to_plane ( FA fa , FB fb )
  {
    tf = [=] ( const crd2_v & in , crd2_v & out )
    {
      crd3_v intermediate ;
      fa.eval ( in , intermediate ) ;
      fb.eval ( intermediate , out ) ;
    } ;
  }

  template < typename FA , typename FB , typename FC >
  plane_to_plane ( FA fa , FB fb , FC fc )
  {
    tf = [=] ( const crd2_v & in , crd2_v & out )
    {
      crd3_v intermediate ;
      fa.eval ( in , intermediate ) ;
      fb.eval ( intermediate , intermediate ) ;
      fc.eval ( intermediate , out ) ;
    } ;
  }

  void eval ( const crd2_v & in , crd2_v & out )
  {
    tf ( in , out ) ;
  }
} ;

// this functor projects a ray to a plane, then produces another
// ray. With the first constructor, there is no in-between
// processing, so if FA and FB are plain conversions, we should
// get the same ray out - provided the projection can handle it.
// 'full' projections like spherical. fisheye and cubemap should
// handle all rays correctly (minus inevitable small errors),
// other projections should do so in their defined range.
// This can be used to test the projections above. Note that
// a test starting out on and returning to the plane via a ray
// intermediate will only work for cubemaps - for sphericals
// and fisheyes, sone locations on the plane have no defined
// associated ray (both poles for sphericals, the 'back pole'
// for fisheyes) - and near those locations, results will be
// imprecise.

template < typename T , std::size_t L >
struct ray_to_ray
: public zimt::unary_functor
    < zimt::xel_t < T , 3 > , zimt::xel_t < T , 3 > , L >
{
  typedef zimt::xel_t < T , 3 > crd3_t ;
  typedef zimt::xel_t < T , 2 > crd2_t ;
  typedef zimt::simdized_type < zimt::xel_t < T , 3 > , L > crd3_v ;
  typedef zimt::simdized_type < zimt::xel_t < T , 2 > , L > crd2_v ;
  std::function < void ( const crd3_v & , crd3_v & ) > tfv ;
  std::function < void ( const crd3_t & , crd3_t & ) > tf ;

  template < typename FA , typename FB >
  ray_to_ray ( FA fa , FB fb )
  {
    tfv = [=] ( const crd3_v & in , crd3_v & out )
    {
      crd2_v intermediate ;
      fa.eval ( in , intermediate ) ;
      fb.eval ( intermediate , out ) ;
    } ;
    tf = [=] ( const crd3_t & in , crd3_t & out )
    {
      crd2_t intermediate ;
      fa.eval ( in , intermediate ) ;
      fb.eval ( intermediate , out ) ;
    } ;
  }

  template < typename FA , typename FB , typename FC >
  ray_to_ray ( FA fa , FB  fb , FC fc )
  {
    tfv = [=] ( const crd3_v & in , crd3_v & out )
    {
      crd2_v intermediate ;
      fa.eval ( in , intermediate ) ;
      fb.eval ( intermediate , intermediate ) ;
      fc.eval ( intermediate , out ) ;
    } ;
    tf = [=] ( const crd3_t & in , crd3_t & out )
    {
      crd2_t intermediate ;
      fa.eval ( in , intermediate ) ;
      fb.eval ( intermediate , intermediate ) ;
      fc.eval ( intermediate , out ) ;
    } ;
  }

  void eval ( const crd3_t & in , crd3_t & out )
  {
    tf ( in , out ) ;
  }

  void eval ( const crd3_v & in , crd3_v & out )
  {
    tfv ( in , out ) ;
  }
} ;

// helper functions to obtain a functor for a plane-to-ray or
// ray-to-plane transformation, given the projection of the
// planar representation. The specific functor is 'grokked'
// to return a uniform type. This is handy to quickly set up
// 'transform stacks' in a type-safe manner - e.g. chaining
// a planar stepper with a to_ray_t, then an operation in ray
// space and finally a to_plane_t for a complete plane-to-plane
// transformation. TODO: what about out-of-bounds access?

// typedef xel_t < float , 2 > crd2_t ;
// typedef xel_t < float , 3 > crd3_t ;
// 
// typedef zimt::grok_type < crd2_t , crd3_t , LANES > to_ray_t ;
// typedef zimt::grok_type < crd3_t , crd2_t , LANES > to_plane_t ;

template < typename T , std::size_t L >
grok_type < xel_t < T , 3 > , xel_t < T , 2 > , L >
roll_out_32 ( projection_t projection )
{
  grok_type < xel_t < T , 3 > , xel_t < T , 2 > , L > result ;

  switch ( projection )
  {
    case SPHERICAL:
      result = ray_to_ll_t < T , L > () ;
      break ;
    case CYLINDRICAL:
      result = ray_to_cyl_t < T , L > () ;
      break ;
    case RECTILINEAR:
      result = ray_to_rect_t < T , L > () ;
      break ;
    case FISHEYE:
      result = ray_to_fish_t < T , L > () ;
      break ;
    case STEREOGRAPHIC:
      result = ray_to_ster_t < T , L > () ;
      break ;
    case CUBEMAP:
      result = ray_to_ir_t < T , L > () ;
      break ;
    case BIATAN6:
      result = ray_to_ba6_t < T , L > () ;
      break ;
    default:
      std::cerr << "unhandled projection # " << int(projection)
                << std::endl ;
      break ;
  }
  return result ;
}

template < typename T , std::size_t L >
grok_type < xel_t < T , 2 > , xel_t < T , 3 > , L >
roll_out_23 ( projection_t projection )
{
  grok_type < xel_t < T , 2 > , xel_t < T , 3 > , L > result ;
  switch ( projection )
  {
    case SPHERICAL:
      result = ll_to_ray_t < T , L > () ;
      break ;
    case CYLINDRICAL:
      result = cyl_to_ray_t < T , L > () ;
      break ;
    case RECTILINEAR:
      result = rect_to_ray_t < T , L > () ;
      break ;
    case FISHEYE:
      result = fish_to_ray_t < T , L > () ;
      break ;
    case STEREOGRAPHIC:
      result = ster_to_ray_t < T , L > () ;
      break ;
    case CUBEMAP:
      result = ir_to_ray_t < T , L > () ;
      break ;
    case BIATAN6:
      result = ba6_to_ray_t < T , L > () ;
      break ;
    default:
      std::cerr << "unhandled projection # " << int(projection)
                << std::endl ;
      break ;
  }
  return result ;
}

// ray-to-ray transformation moving from target image coordinates
// to source image coordinates, and it's inverse. This is a two-step
// rotation, first rotating to model space, then from model space
// onwards. The idea is to be able to 'slot in' more code working
// on the vector in model space - like a shift, or more rotations.
// It doesn't have to be model space - some other intermediate will
// work just as well - e.g. the orientation where the reprojection
// plane of a facet with translation parameters is at (0,0,1). When
// the CS is changed to that CS in a first step, the shift due to the
// positioning of the virtual cameras can be applied, then the second
// rotation moves on to source coordinates. And the inversion simply
// proceeds through the series backwards, using the transposed
// rotations, and the shift with inverted sign.

template < typename T , std::size_t L > // , bool invert = false >
struct tf3d_t
: public zimt::unary_functor < xel_t<T,3> , xel_t<T,3> , L >
{
  const r3_t<T> trg_to_md ;
  const r3_t<T> md_to_trg ;
  const r3_t<T> md_to_src ;
  const r3_t<T> src_to_md ;
  const r3_t<T> trg_to_src ;
  const r3_t<T> src_to_trg ;
  const xel_t<T,3> shift ;
  bool has_shift ;
  T dcp ;

  // the standard c'tor takes the two rotations and an optional shift.

  tf3d_t ( const r3_t<T> & _trg_to_md ,
           const r3_t<T> & _md_to_src ,
           const xel_t<T,3> & _shift = xel_t<T,3> ( 0 ) ,
           T _dcp = T(1) )
  : trg_to_md ( _trg_to_md ) ,
    md_to_trg ( transpose ( _trg_to_md ) ) ,
    md_to_src ( _md_to_src ) ,
    src_to_md ( transpose ( _md_to_src ) ) ,
    trg_to_src ( rotate ( _trg_to_md , _md_to_src ) ) ,
    src_to_trg ( transpose ( rotate ( _trg_to_md , _md_to_src ) ) ) ,
    shift ( _shift ) ,
    has_shift ( _shift[0] != 0 || _shift[1] != 0 || _shift[2] != 0 ) ,
    dcp ( _dcp )
  { }

  // // this c'tor is for convenience. the first two rotations are to
  // // model space and onwards, the third is the rotation from model
  // // space to the reprojection plane. Here, passing the shift is
  // // mandatory - it would be pointless without.
  // 
  // tf3d_t ( const r3_t<T> & _trg_to_md ,
  //          const r3_t<T> & _md_to_src ,
  //          const r3_t<T> & _md_to_tp ,
  //          const xel_t<T,3> & _shift ,
  //          T _dcp = T(1) )
  // : tf3d_t ( rotate ( _trg_to_md , _md_to_tp ) ,
  //            rotate ( transpose ( _md_to_tp ) , _md_to_src ) ,
  //            _shift , _dcp )
  // { }

  template < typename I , typename O >
  void eval ( const I & in , O & out )
  {
    // the transformation is from target image coordinates
    // to source image coordinates - both in model space units.
    // the first case uses two separate rotations and 'slots in'
    // a shift. the second case does it in one go, which is faster.

    if ( has_shift )
    {
      // std::cout << "in b4 rot.       " << in << std::endl ;
      out = rotate ( in , trg_to_md ) ;
      // std::cout << "in after rot.    " << out << std::endl ;
      auto mask = ( out[2] <= 0.0f ) ;
      if ( all_of ( mask ) )
      {
        out[0] = 0.0f ;
        out[1] = 0.0f ;
        out[2] = - std::numeric_limits<float>::infinity() ;
      }
      else
      {
        out[0] /= out[2] ;
        out[1] /= out[2] ;
        out[2] = 1.0f ;
        // std::cout << "after prj to p " << out << std::endl ;

        // the scaling with dcp is only needed when recreating a
        // single source image with '--single' - for the 'normal'
        // direction of transformation, dcp is 1.0

        out *= dcp ;

        out -= shift ;
        // std::cout << "after shift     " << out << std::endl ;
        out = rotate ( out , md_to_src ) ;
        // std::cout << "after final rot " << out << std::endl ;

        if ( any_of ( mask ) )
        {
          out[0] ( mask ) = 0.0f ;
          out[1] ( mask ) = 0.0f ;
          out[2] ( mask ) = - std::numeric_limits<float>::infinity() ;
        }
      }
    }
    else
    {
      out = rotate ( in , trg_to_src ) ;
    }
  }
} ;

END_ZIMT_SIMD_NAMESPACE
HWY_AFTER_NAMESPACE() ;

#endif // sentinel
