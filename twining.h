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

// This header provides code for 'twining': inlined oversampling with
// subsequent weighted averaging.

#include <fstream>
#include "envutil_basic.h"
#include "zimt/unary_functor.h"

#if defined(ENVUTIL_TWINING_H) == defined(HWY_TARGET_TOGGLE)
  #ifdef ENVUTIL_TWINING_H
    #undef ENVUTIL_TWINING_H
  #else
    #define ENVUTIL_TWINING_H
  #endif

HWY_BEFORE_NAMESPACE() ;
BEGIN_ZIMT_SIMD_NAMESPACE(project)

// this might go to environment.h

// class twine_t takes data from a deriv_stepper which yields
// 3D ray coordinates of the pick-up point and two more points
// corresponding to the right and lower neighbours in the target
// image - a fractional step away. It wraps an 'act' type
// functor which can produce pixels (with nchannel channels) from
// 3D ray coordinates. The twine_t object receives the three ray
// coordinates (the 'ninepack'), forms the two differences, and adds
// weighted linear combinations of the differences to the first ray
// (the one representing the pick-up location). The modified rays
// are fed to the wrapped act-type functor in turn, yielding pixel
// values for (slightly) varying pick-up locations in the vicinity
// of the 'central' pick-up location. These pixel values are
// combined in a weighted sum which constitutes the final output.

template < std::size_t nchannels ,
           std::size_t L ,
           bool deriv_tangential = false >
struct twine_t
: public zimt::unary_functor < zimt::xel_t < float , 9 > ,
                               zimt::xel_t < float , nchannels > ,
                               L
                             >
{
  typedef zimt::unary_functor < zimt::xel_t < float , 9 > ,
                                zimt::xel_t < float , nchannels > ,
                                L
                              > base_t ;

  typedef zimt::grok_type < zimt::xel_t < float , 3 > ,
                            zimt::xel_t < float , nchannels > ,
                            L
                          > act_t ;

  act_t inner ;
  std::vector < zimt::xel_t < float , 3 > > spread ;

  using typename base_t::in_ele_v ;

  // new addition: the 'bias' value is a factor which is applied to
  // the 'derivative' values dx and dy once they have been gleaned
  // from the ninepack. The ninepack can be calculated with a
  // distance between the central pickup point and it's neighbours
  // which is proportionally smaller. This is desirable to keep all
  // three sub-pickups confined to the area of a single pixel, to
  // avoid effects due to large ray differences occuring where
  // neighbouring pixels don't correspond to 'neighbouring' rays,
  // e.g. when processing cubemap IR images.

  twine_t ( const act_t & _inner ,
            const std::vector < zimt::xel_t < float , 3 > > & _spread ,
            const float & bias = 4.0 )
  : inner ( _inner ) ,
    spread ( _spread )
  {
    // we apply the 'bias' value to the twining coefficients, to avoid
    // the multiplication with the dx/dy value which would be needed
    // otherwise.

    for ( auto & cf : spread )
    {
      cf[0] *= bias ;
      cf[1] *= bias ;
    }
  }

  // eval function. incoming, we have 'ninepacks' with the ray data.
  // we form the two differences and evaluate 'inner' in a loop,
  // building up the weighted sum.

  template < typename in_type , typename out_type >
  void eval ( const in_type & in ,
              out_type & out )
  {
    typedef typename act_t::in_v crd_v ;

    // ray coordinate of the 'central', unmodified pick-up location

    crd_v pickup { in[0] , in[1] , in[2] } ;

    // we want the 'derivatives' in x and y directions. What we have
    // are two additional rays pointing to the two 'nearby'
    // neighbours of the 'central' pickup ray. If we're sloppy,
    // we can form the 'derivatives' by simple differencing, but
    // if the rays aren't very close together (they are, normally,
    // which is the reason differencing works pretty much
    // always) the vectors received from differencing are not
    // orthogonal to the pickup ray and the plane they define
    // is not the tangential plane of 'pickup', but instead it's
    // tilted a bit. A reasonably inexpensive way of obtaining
    // two vectors in the tangent plane which can serve as a base
    // to calculate the actual pickup locations from the 'spread'
    // is using an orthogonal projection of the neighbouring points
    // onto the tangential plane. We only have to do this once per
    // 'twined' lookup, so the cost is tolerable. To deactivate this
    // code, set deriv_tangential false - plain differencing works
    // most of the time, for this reason: actual pick-up points have
    // to be so closely spaced that they don't 'skip' pixels in the
    // source image - they have to obey the sampling theorem, or
    // otherwise we simply get visibly overlaid image copies (try
    // unduly large twine_width values to see the effect). With
    // such close spacing, the difference between having the actual
    // pick-up points on the tangent plane or on a very slightly
    // tilted plane is negligible. Nevertheless, projecting the
    // neighbours to the tangent plane is mathematically sounder
    // than using simple differencing.
    
    // dx and dy are the two 3D vectors which are multiplied with
    // successive x and y values from each filter coefficient,
    // yielding ray coordinates for the sub-pickups (filter taps)
    // which are weighted with the corresponding w components
    // and then added to the sum.

    crd_v dx , dy ;

    if constexpr ( deriv_tangential )
    {  
      // alternative approach: draw a line through origin and
      // pickup, then, on this line, find the closest point to
      // lp10 or lp01 respectively, and use the diferences as
      // basis vectors. This should compute faster, but have
      // slightly less magnitude than projecting to the
      // tangential plane.

      // the second and third point in the ninepack represent
      // rays - they have unit distance from the origin, but
      // they aren't on 'pickup's' tangent plane, which is where
      // we'd want two vectors as a basis (dx, dy) to produce
      // the additional pick-up points (the in_k below).
      // we start out by manifesting these two points:

      crd_v p10 { in[3] , in[4] , in[5] } ;
      crd_v p01 { in[6] , in[7] , in[8] } ;

      // we use two lines parallel to the 'pickup' ray and passing
      // through p10 and p01, respectively. Imath comes to play:

      Imath::Line3 < in_ele_v > lp10 , lp01 ;

      lp10.pos = p10 ;
      lp10.dir = pickup ;

      lp01.pos = p01 ;
      lp01.dir = pickup ;

      // now we calculate the closest points on these lines to
      // 'pickup': these points are on the tangential plane.
      // they are, to put it differently, the orthogonal projection
      // of p01 and p10 to the tangential plane, and we'll use them
      // to form our basis.

      // zimt and Imath are compatible, but Imath doesn't know that,
      // so we reinterpret_cast. Note how we can form an Imath::Vec3
      // of 'in_ele_v' - a SIMD vector of floats - and use Imath on
      // this data type: the SIMD data type (provided by zimt) has all
      // necessary operators and functions defined to be used by Imath,
      // and the code involved to calculate 'closestPointTo' does not
      // use conditionals, so we're go, and we'll received SIMD results
      // in dx and dy without further ado.

      auto const & pi
        = reinterpret_cast < Imath::Vec3 < in_ele_v > const & > ( pickup ) ; 
      auto & dxi = reinterpret_cast < Imath::Vec3 < in_ele_v > & > ( dx ) ; 
      auto & dyi = reinterpret_cast < Imath::Vec3 < in_ele_v > & > ( dy ) ; 

      dxi = lp10.closestPointTo ( pi ) ;
      dyi = lp01.closestPointTo ( pi ) ;

      // subtracting 'pickup' yields the desired vectors coplanar to
      // the tangential plane

      dx -= pickup ;
      dy -= pickup ;
    }
    else
    {
      // this is the alternative code using simple differencing, which
      // is fine for the 'normal' scenario: very close neighbours.

      dx = crd_v ( { in[3] - in[0] , in[4] - in[1] , in[5] - in[2] } ) ;
      dy = crd_v ( { in[6] - in[0] , in[7] - in[1] , in[8] - in[2] } ) ;
    }

    // with the dx and dy vectors set up, we can now calculate the
    // filter's response:

    // 'out' will receive the weighted sum, we start out by zeroing it
    // px_k is for the value at the sub-pickup as gleaned from 'inner'

    out = 0 ;
    out_type px_k ;

    for ( auto const & coefficient : spread )
    {
      // form the slightly offsetted sub-pickup ray coordinate

      auto in_k = pickup + coefficient[0] * dx + coefficient[1] * dy ;

      // evaluate 'inner' to yield the partial result

      inner.eval ( in_k , px_k ) ;

      // weight it and add it to the final result

      out += coefficient[2] * px_k ;
    }
  }
} ;

END_ZIMT_SIMD_NAMESPACE
HWY_AFTER_NAMESPACE() ;

#endif // sentinel
