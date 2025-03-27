/************************************************************************/
/*                                                                      */
/*    envutil - utility to convert between environment formats          */
/*                                                                      */
/*            Copyright 2025 by Kay F. Jahnke                           */
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

// this header implements use of a panotools-compatible lens correction
// polynomial, using the same semantics as the original. The polynomial
// itself is a trivial thing, what's trickier to implement is the inverse
// of the polynomial, which is desirable to make to geometrical
// transformations in the rendering process fully invertible. Instead
// of even trying to provide a formula for the inverse polynomial, I
// only sample it, using Newton's method, and then set up a b-spline
// over the samples. The result is a functor which performs well enough
// for mass coordinate transformations needed in rendering images
// in contrast to the iteration to find the inverse in a given point which
// is quite slow. The process should be applicable for all 'naturally
// occuring' lens correction polynomials - these should be monotonously
// rising over the range covered, anything else wouldn't make sense
// in the given context.

#if defined(ENVUTIL_LENS_CORRECTION_H) == defined(HWY_TARGET_TOGGLE)
  #ifdef ENVUTIL_LENS_CORRECTION_H
    #undef ENVUTIL_LENS_CORRECTION_H
  #else
    #define ENVUTIL_LENS_CORRECTION_H
  #endif

#include "zimt/eval.h"

HWY_BEFORE_NAMESPACE() ;
BEGIN_ZIMT_SIMD_NAMESPACE(project)

// class eu_polynomial encodes a polynomial of degree DEGREE and
// provides member functions to evaluate the polynomial, it's
// first derivative and it's inverse. The calculation of the
// inverse is done with an iterative process, which isn't guaranteed
// to succeed - but in the context here where we're working with
// quite 'tame' polynomials used for lens correction, failures
// should not occur - we have an assertion in place which will
// throw an exception if this assumption fails. The iteration is
// is done with Newton's method, which succeeds without fail if
// the polynomial's derivative does not change sign, and if we
// take the radius transformation as the polynomial, this should
// certainly be the case for the interval we're looking at.
// With the given three member functions, we add two more with
// 'standard' eval semantics, providing the evaluation of the
// polynomial and the inverse.

template < typename value_type ,
           std::size_t DEGREE ,
           std::size_t L >
struct eu_polynomial
{
  value_type coefficients [ DEGREE + 1 ] ;
  value_type derivative_coefficients [ DEGREE + 1 ] ;
  eu_polynomial() = default ;

  // evaluation of the poynomial and it's derivative is trivial:

  template < typename T >
  T function ( const T & x )
  {
    T sum = 0.0 ;
    T power = 1.0 ;
    for ( std::size_t i = 0 ; i <= DEGREE ; i++ )
    {
      std::size_t j = DEGREE - i ;
      sum += coefficients [ j ] * power ;
      power *= x ;
    }
    return sum ;
  }

  template < typename T >
  T derivative ( const T & x )
  {
    T sum = 0.0 ;
    T power = 1.0 ;
    for ( std::size_t i = 0 ; i < DEGREE ; i++ )
    {
      std::size_t j = DEGREE - i - 1 ;
      sum += derivative_coefficients [ j ] * power ;
      power *= x ;
    }
    return sum ;
  }

  // we use Newton's method to find the inverse. We trust that
  // within the interval we're looking at, the polynomial is
  // rising monotonously. If the result is within the required
  // tolerance, we return true, false otherwise. We're aiming
  // at very high precision; the result we get is usually
  // precise to ca. 2*std::numeric_limits<T>::epsilon(). The
  // iteration tends to find the solution quite quickly. It's
  // essantial to provide a good guess for the start of the
  // iteration (passed in as required_input - see the code in
  // inverse_lcp, where 'inverse' is invoked repeatedly to find
  // knot point values for a b-spline. With a good starting point,
  // we avoid the problem that the iteration might 'drift off'
  // to some other point on the polynomial with the same y value,
  // in case that the entire polynomial is not monotonously rising.

  template < typename T >
  bool inverse ( const T & desired_output ,
                 T & required_input ,
                 const T & tolerance
                   = 100 * std::numeric_limits<T>::epsilon() )
  {
    T current = required_input ;
    T result , difference , last_difference ;

    last_difference = std::numeric_limits<T>::max() ;

    for ( int count = 0 ; count < 16 ; count++ )
    {
      result = function ( current ) ;
      difference = desired_output - result ;

      // std::cout << "*** " << count
      //           << " current: " << current
      //           << " result " << result
      //           << " difference " << difference << std::endl ;

      // we iterate until the difference does not change anymore. This
      // usually takes only a few iterations.

      if ( last_difference == difference )
        break ;

      if ( fabs ( difference ) <= tolerance )
        break ;

      last_difference = difference ;

      current = current + difference / derivative ( current ) ;
    }

    if ( fabs ( difference ) < tolerance )
    {
      required_input = current ;
      return true ;
    }

    return false ;
  }

  // standard eval-semantic member functions for evaluation and
  // inverse

  template < typename T >
  void eval ( const T & x , T & y )
  {
    y = function ( x ) ;
  }

  template < typename T >
  void reval ( const T & y , T & x )
  {
    bool success = inverse ( y , x ) ;
    assert ( success ) ;
  }

  eu_polynomial ( const value_type * p_coefficients )
  {
    std::size_t power = DEGREE ;
    for ( std::size_t i = 0 ; i <= DEGREE ; i++ )
    {
      coefficients [ i ] = p_coefficients [ i ] ;
      derivative_coefficients [ i ] = coefficients [ i ] * power ;
      --power ;
    }
  }

  eu_polynomial ( const std::vector < value_type > & vcf )
  : eu_polynomial ( vcf.data() )
  {
    assert ( vcf.size() == DEGREE + 1 ) ;
  }
} ;

// to calculate the *scaling factor* for the radius which is needed
// to apply the lens correction polynomial, a degree-3 polynomial
// is enough. We provide this as a zimt::unary_functor to be used
// in rendering pipelines. Note that pixel pipelines will need to
// calculate the radius as argument to this functor, but the resulting
// factor will be used typically to scale centered 2D coordinates
// rather than to scale the radius, which would only be useful when
// using polar coordinates.

template < typename T , std::size_t L >
struct lcp
: public zimt::unary_functor < T , T , L > ,
  public eu_polynomial < T , 3 , L >
{
  typedef eu_polynomial < T , 3 , L > pl_t ;

  lcp ( T _a , T _b , T _c , T r_max )
  : pl_t ( { _a , _b , _c , T(1) - ( _a + _b + _c ) } )
  {
    // test whether rr(r) rises over the intended knot point interval
    // not entirely safe - there might be a very small wiggle we don't
    // catch.
    double scale = r_max / 31.0 ;
    double y0 = 0.0 ;
    for ( int i = 1 ; i < 36 ; i++ )
    {
      double x = i * scale ;
      double y ;
      eval ( x , y ) ;
      y *= x ;
      // std::cout << "sample: " << x << " -> " << y << std::endl ;
      assert ( y > y0 ) ;
      y0 = y ;
    }
  }

  using pl_t::eval ; 
} ;

// struct inverse_lcp calculates the radial scaling factor to scale
// a radius value (or a centered 2D coordinate) with lens control
// polynimial applied to a radius value without lcp applied. It's
// used to build the inverse coordinate transformation - inverse to
// the 'normal' direction where coordinate values progress from
// discrete target coordinates to discrete source coordinates.

// The 'normal' direction is, technically speaking, an inverse
// transformation, because it proceeds from target coordinates
// 'back' to source coordinates, but this is the 'normal' direction
// for the type of rendering used in envutil, lux, etc., so we think
// of inverse_lcp as the inverse of the 'normal' operation.

// With this functor, producing images exhibiting lens distortion
// becomes feasible. The inversion of the actual polynomial is
// tricky - here we use an iterative process to arrive at very
// exact double precision values. To make processing fast and also
// allow for SIMDized evaluation, only a set of precise values are
// calculated which are used to create a b-spline. The b-spline
// is evaluated to yield the desired scaling factor. With this
// modelling of the inverse lcp as a b-spline, the functor becomes
// fast enough to employ it in mass data processing - so it's not
// only suitable to work with single coordinates (like when processing
// control point data) but also for rendering jobs with millions of
// pixels. Note the fact that the output is *not the radius* but a
// *scaling factor* for the radius. This functor is used to modify
// centered 2D coordinates, which need to be scaled with the factor.
// If the radius itself were returned by the eval function, the
// scaling factor would have to be calculated after receiving it by
// forming the ratio between it and the current radius, requiring
// a division. By performing the division in double precision during
// the set-up of the inverse_lcp object, the values are more precise
// and the division is no longer needed on the receiving end - the
// 2D coordinate can simply be multiplied with the interpolated factor.

template < typename T , std::size_t L >
struct inverse_lcp
: public zimt::unary_functor < T , T , L >
{
  // even though there are only three coefficients a, b, and c for
  // the lens correction polynomial, the function to yield an
  // appropriately scaled radius for an unscaled one is a fourth-degree
  // polynomial. The fourth coefficient is one minus the sum of the
  // other coefficients, to make the polynomial evaluate to 1.0 at
  // x = 1.0. The last coefficient has to be zero, So that f(0) == 0

  double coefficients[5] ;

  // both radii are chosen slightly higher than the incoming argument
  // to have a few extra knot points towards the end of the spline.
  // if the lens correction polynomial yealing the scaled radius is
  // p(x) = a xxxx + b xxx + c xx + ( 1 - a - b - c ) x
  // we have
  // rr_max = p ( r_max )

  double rr_max ;  // maximal incoming radius to transform back
  double r_max ;   // maximal original radius

  // b-spline and evaluator use type T - usually float - because they
  // are employed in rendering jobs. When working CPs, use double as T.
  // We take note of the knot point count as 'nk'

  int nk ;
  zimt::bspline < T , 1 > inv_model ;
  zimt::grok_type < T , T , L > fev ;

  // We use a cubic spline, which combines good precision and reasonable
  // processing times. I was unsure initially whether using only control
  // points starting at zero would be okay, but since the lcp is a
  // polynomial through the origin without DC and predominant linear
  // component, the curve is pretty much a straight line near zero,
  // so using NATURAL BCs is fine - putting a zero into the second
  // derivative at origin is close to the actual curve, which has next
  // to no curvature at the origin anyway. At the right end of the spline,
  // there is curvature, and putting a zero into the second derivative
  // there is not safe - that's why it makes sense to add a few
  // extra control points at the tail end - we're silently assuming
  // that for the extra CPs here the polynomial is still rising.
  // I found that near the origin the resolution which results from
  // linear steps in the radius as knot point positions for the spline
  // is often insufficient unless the number of knot points is increased
  // drastically. So I add a nonlinear transformation to the stepping
  // function, resulting in ever smaller steps near the origin, at the
  // cost of larger ones near the upper end of the spline.

  inverse_lcp ( double _a , double _b , double _c ,
                double _r_max , int sz = 32 )
  : coefficients { _a , _b , _c , 1.0 - ( _a + _b + _c ) , 0.0 } ,
    inv_model ( sz + 4 , 3 , zimt::NATURAL ) ,
    r_max ( _r_max * ( ( sz + 3.0 ) / sz ) )
  {
    eu_polynomial < double , 4 , L > p ( coefficients ) ;
    rr_max = p.function ( r_max ) ;
    nk = inv_model.core.shape[0] ; // == sz + 4

    for ( std::size_t i = 0 ; i < nk ; i++ )
    {
      // now we set the spline's control points, proceeding
      // some steps beyond the coefficient corresponding to the
      // radius at the corner of the image.

      // first we calculate 'notches' from zero to 1.0

      double notch = double(i) / ( nk - 1 ) ;

      // now we use a nonlinear scaling function. The resulting
      // range of notches is still [0...1]

      notch *= notch ;

      // now we extend the range to [0...rr_max]

      notch *= rr_max ;

      // we want a good starting point for the iteration. A reasonable
      // approximation is to form a straight line from the origin to
      // end of the covered interval and intersect it with the desired
      // y value

      double out = i * r_max / sz ;

      p.reval ( notch , out ) ; // will throw if inverse isn't found

      // as the b-spline's knot points, we don't store the radius
      // but rather the factor we need to apply to the incoming radius
      // to obtain the outgoing radius.

      inv_model.core [ i ] = (    notch == 0.0
                                ? 1.0 / p.derivative ( 0.0 )
                                : out / notch ) ;
    }

    // new we prefilter and set up an evaluator for the spline.

    inv_model.prefilter() ;

    fev = zimt::make_safe_evaluator < decltype ( inv_model ) , T , L >
      ( inv_model ) ;
  }

  // the eval function first scales the incoming radius value (for which
  // we seek the scaling factor to transform it back to the 'original'
  // radius value) to the interval [0...1], then applies the inverse of
  // the nonlinear sampling function we used to set up the spline's knot
  // points, and finally moves from [0...1] to the spline's range to
  // evaluate. output is the desired scaling factor, which can be used
  // directly to scale a radius or centered 2D coordinates.
  // The square root in the evaluation is annoying (costly to compute)
  // but the nonlinear sampling reduces the spline's knot point count
  // by a fair deal. An alternative is to omit the nonlinear sampling
  // and increase the knot point count considerably.

  template < typename I , typename O >
  void eval ( const I & _in , O & out )
  {
    auto in = _in / rr_max ; // move to [0...1]
    in = sqrt ( in ) ;       // apply inverse of nonlinear stepping
    in *= ( nk - 1 ) ;       // move to spline coordinates

    fev.eval ( in , out ) ;  // evaluate the spline
  }
} ;

// a class for the entire planar transformations occuring in PTO.
// with 'invert' false, we have the 'normal' operation from target
// to source coordinates. 'invert' true yields the inverse of that
// transformation - it should be needed less frequently, e.g. for
// the inspection of control points or creation of synthetic images
// with the same geometrical flaws as given input facets from an
// already-stitched panorama. Note the terminological ambiguity:
// 'normal' operation is actually an inverse transformation - namely
// from target to source coordinates. Here the template parameter
// 'inverse' refers to the inverse of the 'normal' operation.

template < typename T , std::size_t L , bool invert = false >
struct pto_planar
: public zimt::unary_functor < xel_t<T,2> , xel_t<T,2> , L >
{
  lcp < T , L > polynomial ;
  inverse_lcp < T , L > inv_polynomial ;
  const T h , v , g , t ;

  pto_planar ( double _a , double _b , double _c , double r_max ,
               double _d = 0.0 , double _e = 0.0 ,
               double _g = 0.0 , double _t = 0.0 )
  : polynomial ( _a , _b , _c , r_max ) ,
    inv_polynomial ( _a , _b , _c , r_max , 57 ) ,
    h ( _d ) ,
    v ( _e ) ,
    g ( _g ) ,
    t ( _t )
  { }

  template < typename I , typename O >
  void eval ( const I & _in , O & out )
  {
    if constexpr ( invert == false )
    {
      // the transformation is from target image coordinates
      // to source image coordinates - both in model space units.

      T factor ;
      polynomial.eval ( norm(_in) , factor ) ;
      out = _in * factor ;
      if ( h != 0.0 || v != 0.0 )
      {
        out[0] += h ;
        out[1] += v ;
      }
      if ( g != 0.0 || t != 0.0 )
      {
        out = { out[0] + ( out[1] * g ) ,
                out[1] + ( out[0] * t ) } ;
      }
    }
    else
    {
      // the transformation is from source image coordinates
      // to target image coordinates - both in model space units.

      out = _in ;
      if ( g != 0.0 || t != 0.0 )
      {
        // I adapted this formula from panotools' math.c (see shearInv):
        out[1]  = (_in[1] - t * _in[0]) / (1 - t * g);
        out[0]  = (_in[0] - g * out[1]);
      }
      if ( h != 0.0 || v != 0.0 )
      {
        out[0] -= h ;
        out[1] -= v ;
      }
      T factor ;
      inv_polynomial.eval ( norm(out) , factor ) ;
      out *= factor ;
    }
  }
} ;

END_ZIMT_SIMD_NAMESPACE
HWY_AFTER_NAMESPACE() ;

#endif // sentinel
