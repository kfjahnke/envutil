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

// coordinate transformations, coded as templates in zimt 'act'
// functor style, returning the result via a reference argument

// code to convert lat/lon to 3D ray. Note that I am using lux
// coordinate convention ('book order': x is right, y down and
// z forward). The lat/lon values coming in are angles in radians,
// and the resulting 'ray' coordinates are in 'model space' units.
// We're using an implicit radius of 1.0 for the conversion from
// spherical to cartesian coordinates, so the output has unit length.

struct ll_to_ray_t
: public zimt::unary_functor < v2_t , v3_t , LANES >
{
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

    #if defined USE_HWY or defined USE_VC

      f_v sinlat , coslat , sinlon , coslon ;
      sincos ( lat , sinlat , coslat ) ;
      sincos ( lon , sinlon , coslon ) ;

    #else

      f_v sinlat = sin ( lat ) ;
      f_v coslat = cos ( lat ) ;
      f_v sinlon = sin ( lon ) ;
      f_v coslon = cos ( lon ) ;

    #endif

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

/*

template < face_index_t F , int nchannels >
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

*/

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

template < int nchannels >
struct ir_to_ray
: public zimt::unary_functor < v2_t , v3_t , LANES >
{
  const double section_md ;
  const double refc ;

  ir_to_ray ( double _section_md , double _refc )
  : section_md ( _section_md ) ,
    refc ( _refc )
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
    crd2    -= refc ;

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

template < int nchannels >
struct ir_to_exr
: public zimt::unary_functor < v2_t , v3_t , LANES >
{
  const double section_md ;
  const double refc ;

  ir_to_exr ( double _section_md , double _refc )
  : section_md ( _section_md ) ,
    refc ( _refc )
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
    crd2    -= refc ;

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

