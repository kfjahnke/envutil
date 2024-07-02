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

#ifndef TWINING_H
#define TWINING_H

// This header provides code for 'twining': inlined oversampling with
// subsequent weighted averaging.

#include <fstream>
#include "common.h"

// this function sets up a simple box filter to use with the
// twine_t functor above. The given w and h values determine
// the number of pick-up points in the horizontal and vertical
// direction - typically, you'd use the same value for both.
// The deltas are set up so that, over all pick-ups, they produce
// a uniform sampling.
// additional parameters allow to apply gaussian weights and to
// apply a threshold to suppress small weighting factors; in a
// final step the weights are normalized.

void make_spread ( std::vector < zimt::xel_t < float , 3 > > & trg ,
                   int w = 2 ,
                   int h = 0 ,
                   float d = 1.0f ,
                   float sigma = 0.0f ,
                   float threshold = 0.0f )
{
  if ( w <= 2 )
    w = 2 ;
  if ( h <= 0 )
    h = w ;
  float wgt = 1.0 / ( w * h ) ;
  double x0 = - ( w - 1.0 ) / ( 2.0 * w ) ;
  double dx = 1.0 / w ;
  double y0 = - ( h - 1.0 ) / ( 2.0 * h ) ;
  double dy = 1.0 / h ;
  trg.clear() ;
  sigma *= - x0 ;
  double sum = 0.0 ;

  for ( int y = 0 ; y < h ; y++ )
  {
    for ( int x = 0 ; x < w ; x++ )
    {
      float wf = 1.0 ;
      if ( sigma > 0.0 )
      {
        double wx = ( x0 + x * dx ) / sigma ;
        double wy = ( y0 + y * dy ) / sigma ;
        wf = exp ( - sqrt ( wx * wx + wy * wy ) ) ;
      }
      zimt::xel_t < float , 3 >
        v { float ( d * ( x0 + x * dx ) ) ,
            float ( d * ( y0 + y * dy ) ) ,
            wf * wgt } ;
      trg.push_back ( v ) ;
      sum += wf * wgt ;
    }
  }

  double th_sum = 0.0 ;
  bool renormalize = false ;

  if ( sigma != 0.0 )
  {
    for ( auto & v : trg )
    {
      v[2] /= sum ;
      if ( v[2] >= threshold )
      {        
        th_sum += v[2] ;
      }
      else
      {
        renormalize = true ;
        v[2] = 0.0f ;
      }
    }
    if ( renormalize )
    {
      for ( auto & v : trg )
      {
        v[2] /= th_sum ;
      }
    }
  }

  if ( args.verbose )
  {
    if ( sigma != 0.0 )
    {
      std::cout << "using this twining filter kernel:" << std::endl ;
      for ( int y = 0 ; y < h ; y++ )
      {
        for ( int x = 0 ; x < w ; x++ )
        {
          std::cout << '\t' << trg [ y * h + x ] [ 2 ] ;
        }
        std::cout << std::endl ;
      }
    }
    else
    {
      std::cout << "using box filter for twining" << std::endl ;
    }
  }

  if ( renormalize )
  {
    auto help = trg ;
    trg.clear() ;
    for ( auto v : help )
    {
      if ( v[2] > 0.0f )
        trg.push_back ( v ) ;
    }
    if ( args.verbose )
    {
      std::cout << "twining filter taps after after thresholding: "
                << trg.size() << std::endl ;
    }
  }
}

// read the twining filter from a twf-file.

void read_twf_file ( std::vector < zimt::xel_t < float , 3 > > & trg )
{
  zimt::xel_t < float , 3 > c ;

  std::ifstream ifs ( args.twf_file ) ;
  assert ( ifs.good() ) ;

  double sum = 0.0 ;

  while ( ifs.good() )
  {
    ifs >> c[0] >> c[1] >> c[2] ;
    if ( ifs.eof() )
      break ;
    trg.push_back ( c ) ;
    sum += c[2] ;
  }

  for ( auto & c : trg )
  {
    c[0] *= args.twine_width ;
    c[1] *= args.twine_width ;
    if ( args.twine_normalize )
      c[2] /= sum ;
  }

  if ( verbose )
  {
    std::cout << args.twf_file << " yields twining filter kernel:"
              << std::endl ;
    for ( const auto & c : trg )
    {
       std::cout << "x: " << c[0] << " y: " << c[1]
                 << " w: " << c[2] << std::endl ;
    }
    if ( args.twine_normalize )
    {
      std::cout << "twining filter weights sum: 1.0" << std::endl ;
    }
    else
    {
      std::cout << "twining filter weights sum: " << sum << std::endl ;
    }
  }
  ifs.close() ;
}

// class twine_t takes data from a deriv_stepper which yields
// 3D ray coordinates of the pick-up point and two more points
// corresponding to the right and lower neighbours in the target
// image - so, one 'canonical' step away. It wraps an 'act' type
// functor which can produce pixels (with nchannel channels) from
// 3D ray coordinates. The twine_t object receives the three ray
// coordinates (the 'ninepack'), forms the two differences, and adds
// weighted linear combinations of the differences to the first ray
// (the one representing the pick-up location). The modified rays
// are fed to the wrapped act-type functor in turn, yielding pixel
// values for (slightly) varying pick-up locations in the vicinity
// of the 'central' pick-up location. These pixel values are
// combined in a weighted sum which constitutes the final output.

template < std::size_t nchannels , std::size_t L >
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
  const std::vector < zimt::xel_t < float , 3 > > spread ;

  using typename base_t::in_ele_v ;

  twine_t ( const act_t & _inner ,
            const std::vector < zimt::xel_t < float , 3 > > & _spread )
  : inner ( _inner ) ,
    spread ( _spread )
  { }

  // eval function. incoming, we have 'ninepacks' with the ray data.
  // we form the two differences and evaluate 'inner' in a loop,
  // building up the weighted sum.

  template < typename in_type , typename out_type >
  void eval ( const in_type & in ,
              out_type & out )
  {
    out = 0 ;
    out_type px_k ;

    typedef typename act_t::in_v crd_v ;

    // ray coordinate of the 'central', unmodified pick-up location

    crd_v pickup { in[0] , in[1] , in[2] } ;

    // we want the 'derivatives' in x and y directions. What we have
    // are two additional rays pointing to the two next 'canonial'
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
    // 'twined' lookup, so the cost is tolerable. To activate this
    // code, #define DERIV_TANGENTIAL - I use plain differencing for
    // now, for this reason: if the twining is to work like a filter,
    // the actual pick-up points have to be so closely spaced that
    // they don't 'skip' pixels in the source image - they have to
    // obey the sampling theorem, or otherwise we simply get visibly
    // overlaid image copies (try unduly large twine_width values
    // to see the effect). With such close spacing, the difference
    // between having the actual pick-up points on the tangent plane
    // or on a very slightly tiled plane is negligible.

// KFJ 2024-07-02 now using the projection to the tangent plane.
    
#define DERIV_TANGENTIAL

#ifdef DERIV_TANGENTIAL
  
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

    crd_v dx , dy ;

    // zimt and Imath are compatible, but Imath doesn't know that,
    // so we reinterpret_cast:

    auto const & pi = reinterpret_cast < Imath::Vec3 < in_ele_v > const & >
                        ( pickup ) ; 
    auto & dxi = reinterpret_cast < Imath::Vec3 < in_ele_v > & > ( dx ) ; 
    auto & dyi = reinterpret_cast < Imath::Vec3 < in_ele_v > & > ( dy ) ; 

    dxi = lp10.closestPointTo ( pi ) ;
    dyi = lp01.closestPointTo ( pi ) ;

    // subtracting 'pickup' yields the desired vectors coplanar to
    // the tangential plane

    dx -= pickup ;
    dy -= pickup ;

#else

    // this is the alternative code using simple differencing, which
    // is fine for the 'normal' scenario: very close neighbours.

    crd_v dx { in[3] - in[0] , in[4] - in[1] , in[5] - in[2] } ;
    crd_v dy { in[6] - in[0] , in[7] - in[1] , in[8] - in[2] } ;

#endif

    for ( auto const & contrib : spread )
    {
      // form the slightly offsetted pick-up ray coordinate

      auto in_k = pickup + contrib[0] * dx + contrib[1] * dy ;

      // evaluate 'inner' to yield the partial result

      inner.eval ( in_k , px_k ) ;

      // weight it and add it to the final result

      out += contrib[2] * px_k ;
    }
  }
} ;

#endif // TWINING_H
