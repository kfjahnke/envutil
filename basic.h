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

// this header has basic enums and types which do not depend on zimt,
// and some helper functions which don't use zimt.

#ifndef ENVUTIL_BASIC_H
#define ENVUTIL_BASIC_H

#include <cmath>

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

// assuming isotropic sampling (same sampling resolution in the horizontal
// and vertical), calculate the vertical field of view from the horizontal
// field of view, under the given projection.
// Note that this function is for centered images only
// (x0 == -x1), (y0 == -y1)

double get_vfov ( projection_t projection ,
                  int width ,
                  int height ,
                  double hfov ) ;

// the 'step' of an image is the angle - in radians - which
// corresponds to the width of one pixel in the image center.
// for some projections and in certain directions, this value
// will be usable at non-central points (e.g. for spherical
// images along the horizon). In any case it can be used as a
// 'rule of thumb' indicator of the image's resolution.
// If we have the 'extent' of a spherical or cylindrical image
// already, we can calculate the step as hfov / width;
// for other projections this simple formula doesn't apply.
// Note that this function is for centered images only
// (x0 == -x1), (y0 == -y1)
// currently unused.

double get_step ( projection_t projection ,
                  int width ,
                  int height ,
                  double hfov ) ;

// extent_type is a handy struct to contain the horizontal and
// vertical extent of a rectangular section of the plane.

struct extent_type
{
  double x0 , x1 , y0 , y1 ;
} ;

// extract internally uses the notion of an image's 'extent' in 'model
// space'. The image is thought to be 'draped' to an 'archetypal 2D
// manifold' - the surface of a sphere or cylinder with unit radius
// or a plane at unit distance forward - where the sample points are
// placed on the 2D manifold so that rays from the origin to the
// scene point which corresponds with the sample point intersect there.
// To put it differently: the sample point cloud is scaled and shifted
// to come to lie on the 'archetypal' 2D manifolds. This makes for
// efficient calculations. The image is taken to be centered on the
// 'forward' ray.

extent_type get_extent ( projection_t projection ,
                         int width ,
                         int height ,
                         double hfov ) ;

#endif // #ifndef ENVUTIL_BASIC_H
