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

#ifndef GEOMETRY_H
#define GEOMETRY_H

// This header has all the geometrical transformations used to convert
// between lat/lon and cubemap coordinates. Then it proceeds to give
// conversions from coordinates pertaining to the IR image to 3D 'ray'
// coordinates, both in lux and openEXR convention.

#include "common.h"

// we start out with conventions used for directions and how they
// are encoded in 3D coordinates stored as zimt::xel_t of three
// fundamentals.

// enum encoding the sequence of cube face images in the cubemap
// This is the sequence used for openEXR cubmap layout. The top
// and bottom squares are oriented so as to align with the back
// image. Of course, the labels are debatable: my understanding
// of 'front' is 'aligned with the image center'. If one were to
// associate 'front' with the wrap-around point of the full
// spherical, the labels would be different.

typedef enum
{
  CM_LEFT ,
  CM_RIGHT ,
  CM_TOP ,
  CM_BOTTOM ,
  CM_FRONT ,
  CM_BACK
} face_index_t ;

// we use lux coordinate system convention. I call it 'latin book
// order': if you have a stack of prints in front of you and read
// them, your eyes move first left to right inside the line, then
// top to bottom from line to line, then, moving to the next pages,
// forward in the stack. Using this order also makes the first two
// components agree with normal image indexing conventions, namely
// x is to the right and y down. Note that I put the fastest-moving
// index first, which is 'fortran' style, whereas C/C++ use the
// opposite order for nD arrays.

enum { RIGHT , DOWN , FORWARD } ;

// openEXR uses different 3D axis semantics, and if we want to use
// OIIO's environment lookup function, we need openEXR 3D coordinates.

// Here's what the openEXR documentation sys about their axis
// order (next to a drawing which says differently, see this issue:
// https://github.com/AcademySoftwareFoundation/openexr/issues/1687)

// quote:
// We assume that a camera is located at the origin, O, of a 3D
// camera coordinate system. The camera looks along the positive z
// axis. The positive x and y axes correspond to the camera’s left
// and up directions.
// end quote

// so we'd get this axis order, assuming they store x,y,z:

enum { EXR_LEFT , EXR_UP , EXR_FORWARD } ;

// the cubemap comes out right this way, so I assume that their text
// is correct and the drawing is wrong.

typedef enum
{
  SPHERICAL ,
  CYLINDRICAL ,
  RECTILINEAR ,
  STEREOGRAPHIC ,
  FISHEYE ,
  CUBEMAP ,
  PRJ_NONE
}  projection_t ;

const char * const projection_name[]
{
  "spherical" ,
  "cylindrical" ,
  "rectilinear" ,
  "stereographic" ,
  "fisheye" ,
  "cubemap"
} ;

struct extent_type
{
  double x0 , x1 , y0 , y1 ;
} ;

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
    < zimt::xel_t < T , 2 > ,
      zimt::xel_t < T , 3 > ,
      L
    >
{
  typedef zimt::simdized_type < T , L > f_v ;

  template < typename in_type , typename out_type >
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
    // The code benefits from using sincos, where available.

    // #if defined USE_HWY or defined USE_VC

      f_v sinlat , coslat , sinlon , coslon ;
      sincos ( lat , sinlat , coslat ) ;
      sincos ( lon , sinlon , coslon ) ;

    // #else
    // 
    //   f_v sinlat = sin ( lat ) ;
    //   f_v coslat = cos ( lat ) ;
    //   f_v sinlon = sin ( lon ) ;
    //   f_v coslon = cos ( lon ) ;
    // 
    // #endif

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

struct ray_to_ll_t
: public zimt::unary_functor < v3_t , v2_t , LANES >
{
  template < typename in_type , typename out_type >
  void eval ( const in_type & in ,
              out_type & out ) const
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

// this functor template converts incoming in-face coordinates
// to ray coordinates for a given face index, which is passed
// as a template argument - so the sixfold 'if constexpr ...' is
// not a conditional, it's just a handy way of putting the code
// into a single function without having to write partial template
// specializations for the six possible face indices.
// currently unused.

template < face_index_t F >
struct in_face_to_ray
: public zimt::unary_functor < v2_t , v3_t , LANES >
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
  void eval ( const I & crd2 , O & crd3 ) const
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
// to the entire IR array to 3D ray coordinates in lux convention.
// This functor can serve to populate the IR image: set up a
// functor yielding model space coordinates pertaining to pixels
// in the IR image, pass these model space coordinates to this
// functor, receive ray coordinates, then glean pixel values
// for the given ray by evaluating some functor taking ray
// coordinates and yielding pixels.
// The functor takes two arguments: first the 'section size':
// this is equal to the IR image's width, expressed in model
// space units. If the IR image does not have additional support
// and holds cube face images of precisely ninety degrees fov,
// the value would be 2.0 precisely. With added support, it's
// slightly larger. The second argument is the distance, in model
// space units, from the upper left corner of a section to the
// cube face image's center. If the cube face image has even width,
// this is precisely half the section size, but with odd width,
// this isn't possible, hence the extra argument.

struct ir_to_ray
: public zimt::unary_functor < v2_t , v3_t , LANES >
{
  const double section_md ;
  const double refc_md ;

  ir_to_ray ( double _section_md , double _refc_md )
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
} ;

// to make the conversion efficient and transparent, I refrain from
// using ir_to_ray_gen and a subsequent coordinate transformation
// in favour of this dedicated functor, which is basically a copy of
// the one above, but with different component indexes and signs
// inverted where necessary (namely for the left and up direction,
// which are the negative of zimt's right and down).

struct ir_to_exr
: public zimt::unary_functor < v2_t , v3_t , LANES >
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

// given a 3D 'ray' coordinate, find the corresponding cube face
// and the in-face coordinate - note the two references which take
// the result values. The incoming 'ray' coordinate does not have
// to be normalized. The resulting in-face coordinates are in the
// range of [-1,1] - in 'model space units' pertaining to planes
// 'draped' at unit distance from the origin and perpendicular
// to one of the axes. This function is coded as pure SIMD code,
// we don't need it for scalars.

void ray_to_cubeface ( const crd3_v & c ,
                       i_v & face ,
                       crd2_v & in_face )
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

  // set dom true where the z axis (pointing forward) has the
  // largest numerical value
  
  dom = ( ! m2 ) & ( ! m3 ) ;

  if ( any_of ( dom ) )
  {
    // finally the front and back faces

    face ( dom ) = CM_FRONT ;
    face ( dom & ( c[FORWARD] < 0 ) ) = CM_BACK ;

    in_face[0] ( dom ) = c[RIGHT] / c[FORWARD] ;
    in_face[1] ( dom ) = c[DOWN] / abs ( c[FORWARD] ) ;
  }
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

void ray_to_cubeface_fixed ( const crd3_v & c ,
                             const i_v & face ,
                             crd2_v & in_face )
{
  // form a mask which is true where a specific axis is 'dominant'.
  // since we have the face indices already, this is simple: it's
  // the face indices shifted to the right to remove their least
  // significant bit, which codes for the sign along the dominant
  // axis.

  auto dom_v = face >> 1 ;

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

  auto dom = ( dom_v == 0 ) ;
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
    // Since we know the face value already, all we need here is
    // the deivision by the dominant component.

    in_face[0] ( dom ) = - c[FORWARD] / c[RIGHT] ;
    in_face[1] ( dom ) = c[DOWN] / abs ( c[RIGHT] ) ;
  }

  // now set dom true where the y axis (pointing down) has the
  // largest numerical value

  dom = ( dom_v == 1 ) ; 

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

    // lux convention:
    // in_face[0] ( dom ) =   c[RIGHT] / abs ( c[DOWN] ) ;
    // in_face[1] ( dom ) = - c[FORWARD] / c[DOWN] ;

    in_face[0] ( dom ) = - c[RIGHT] / abs ( c[DOWN] ) ;
    in_face[1] ( dom ) =   c[FORWARD] / c[DOWN] ;

  }

  // set dom true where the z axis (pointing forward) has the
  // largest numerical value
  
  dom = ( dom_v == 2 ) ;

  if ( any_of ( dom ) )
  {
    // finally the front and back faces

    in_face[0] ( dom ) = c[RIGHT] / c[FORWARD] ;
    in_face[1] ( dom ) = c[DOWN] / abs ( c[FORWARD] ) ;
  }
}

// this functor converts incoming 3D ray coordinates to 2D
// coordinates pertaining to the IR image (in model space units)
// Similar code is used in sixfold_t - here, we have a stand-alone
// version which only needs two values from the metrics_t object.
// This functor is the inversion of ir_to_ray, which needs the
// same two values only. With this pair, we can slot the conversions
// into lux' rendering code.

struct ray_to_ir
: public zimt::unary_functor < v3_t , v2_t , LANES >
{
  // to obtain IR coordinates, we need the width of a 'section'
  // of the IR image in model space units and the distance to
  // the center of the cube face images from the left/top margin:

  float section_md ;
  float refc_md ;

  ray_to_ir ( const float & _section_md ,
              const float & _refc_md )
  : section_md ( _section_md ) ,
    refc_md ( _refc_md )
  { }

  template < typename I , typename O >
  void eval ( const I & c , O & ir_c ) const
  {
    // first, we glean the face indices and in-face coordinates

    i_v face ;

    ray_to_cubeface ( c , face , ir_c ) ;

    ir_c += refc_md ;

    // for the vertical, we add the face index times the section
    // size (in model space units)

    ir_c[1] += face * section_md ;
  }
} ;

// when doing geometrical transformations of images, the coordinate
// transformations can often be expressed as the chaining of two
// functors: one transforming the incoming coordinate to a 3D ray,
// and a second one transforming the 3D ray to another 2D coordinate.
// struct plane_to_plane odders just such a construct, optionally
// with a 3D-to-3D operation inserted in between, e.g. a 3D rotation.
// This is pure SIMDized code, there is no scalar overload.

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

  template < typename FA , typename F3 , typename FB >
  plane_to_plane ( FA fa , F3  f3 , FB fb )
  {
    tf = [=] ( const crd2_v & in , crd2_v & out )
    {
      crd3_v intermediate ;
      fa.eval ( in , intermediate ) ;
      f3.eval ( intermediate , intermediate ) ;
      fb.eval ( intermediate , out ) ;
    } ;
  }

  void eval ( const crd2_v & in , crd2_v & out )
  {
    tf ( in , out ) ;
  }
} ;

#endif // #ifndef GEOMETRY_H
