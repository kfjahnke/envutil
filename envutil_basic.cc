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

// this file has some helper functions which don't use zimt and which
// aren't performance-critical. The code is largely from lux.

#include "envutil_basic.h"

// assuming isotropic sampling (same sampling resolution in the horizontal
// and vertical), calculate the vertical field of view from the horizontal
// field of view, under the given projection.
// Note that this function is for centered images only
// (x0 == -x1), (y0 == -y1)

double get_vfov ( projection_t projection ,
                  int width ,
                  int height ,
                  double hfov )
{
  double vfov = 0.0 ;
  switch ( projection )
  {
    case RECTILINEAR:
    {
      // as a one-liner, this is probably clearer than the code below
      vfov = 2.0 * atan ( height * tan ( hfov / 2.0 ) / width ) ;
      break ;
    }
    case CYLINDRICAL:
    {
      double pixels_per_rad = width / hfov ;
      double h_rad = height / pixels_per_rad ;
      vfov = 2.0 * atan ( h_rad / 2.0 ) ;
      break ;
    }
    case STEREOGRAPHIC:
    {
      double w_rad = 2.0 * tan ( hfov / 4.0 ) ;
      double pixels_per_rad = width / w_rad ;
      double h_rad = height / pixels_per_rad ;
      vfov = 4.0 * atan ( h_rad / 2.0 ) ;
      break ;
    }
    case SPHERICAL:
    case FISHEYE:
    {
      vfov = hfov * height / width ;
      break ;
    }
    case CUBEMAP:
    case BIATAN6:
    {
      vfov = 2.0 * M_PI ;
    }
    default:
    {
      vfov = hfov ; // debatable...
      break ;
    }
  }
  return vfov ;
}

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
                  double hfov )
{
  double step = 0.0 ;
  switch ( projection )
  {
    case RECTILINEAR:
    case CUBEMAP:
    {
      step = atan ( 2.0 * tan ( hfov / 2.0 ) / width ) ;
      break ;
    }
    case BIATAN6:
    case SPHERICAL:
    case CYLINDRICAL:
    case FISHEYE:
    {
      step = hfov / width ;
      break ;
    }
    case STEREOGRAPHIC:
    {
      step = atan ( 4.0 * tan ( hfov / 4.0 ) / width ) ;
      break ;
    }
    default:
    {
      break ;
    }
  }
  return step ;
}

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
                         double hfov )
{
  double x0 , x1 , y0 , y1 ;

  double alpha_x = - hfov / 2.0 ;
  double beta_x = hfov / 2.0 ;
  double beta_y = get_vfov ( projection , width , height , hfov ) / 2.0 ;
  double alpha_y = - beta_y ;

  switch ( projection )
  {
    case SPHERICAL:
    case FISHEYE:
    {
      x0 = alpha_x ;
      x1 = beta_x ;

      y0 = alpha_y ;
      y1 = beta_y ;
      break ;
    }
    case CYLINDRICAL:
    {
      x0 = alpha_x ;
      x1 = beta_x ;

      y0 = tan ( alpha_y ) ;
      y1 = tan ( beta_y ) ;
      break ;
    }
    case RECTILINEAR:
    {
      x0 = tan ( alpha_x ) ;
      x1 = tan ( beta_x ) ;

      y0 = tan ( alpha_y ) ;
      y1 = tan ( beta_y ) ;
      break ;
    }
    case STEREOGRAPHIC:
    {
      x0 = 2.0 * tan ( alpha_x / 2.0 ) ;
      x1 = 2.0 * tan ( beta_x / 2.0 ) ;

      y0 = 2.0 * tan ( alpha_y / 2.0 ) ;
      y1 = 2.0 * tan ( beta_y / 2.0 ) ;
      break ;
    }
    case CUBEMAP:
    case BIATAN6:
    {
      x0 = tan ( alpha_x ) ;
      x1 = tan ( beta_x ) ;

      y0 = 6 * x0 ;
      y1 = 6 * x1 ;
      break ;
    }
    default:
    {
      x0 = x1 = y0 = y1 = 0.0 ;
      break ;
    }
  }
  return { x0 , x1 , y0 , y1 } ;
}
