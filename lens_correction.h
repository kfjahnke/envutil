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
// throw an exception if this assumption fails.
// With the given three member functions, we add two more with
// 'standard' eval semantics, providing the evaluation of the
// polynomial and the inverse.

template < typename value_type ,
           std::size_t DEGREE ,
           std::size_t LANES >
struct eu_polynomial
{
  value_type coefficients [ DEGREE + 1 ] ;
  value_type derivative_coefficients [ DEGREE + 1 ] ;
  eu_polynomial() = default ;

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
  // rising monotonously.

  template < typename T >
  bool inverse ( const T & desired_output ,
                 T & required_input ,
                 const T & tolerance
                   = 100 * std::numeric_limits<T>::epsilon() )
  {
    T current = required_input ;
    T result , difference ;
    // we'll try to get really close:
    T aim = 3 * std::numeric_limits<T>::epsilon() ;

    for ( int count = 0 ; count < 10 ; count++ )
    {
      result = function ( current ) ;
      difference = desired_output - result ;

      // std::cout << "current: " << current
      //           << " deriv: " << derivative ( current )
      //           << " result: " << result
      //           << " diff: " << difference
      //           << " tol " << tolerance
      //           << std::endl ;

      if ( fabs ( difference ) <= aim )
      {
        required_input = current ;
        return true ;
      }

      current = current + difference / derivative ( current ) ;
    }
    if ( fabs ( difference ) < tolerance )
      return true ;

    return false ;
  }

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

template < typename T , std::size_t LANES >
struct lcp
: public zimt::unary_functor < T , T , LANES > ,
  public eu_polynomial < T , 3 , LANES >
{
  typedef eu_polynomial < T , 3 , LANES > pl_t ;

  lcp ( T _a , T _b , T _c )
  : pl_t ( { _a , _b , _c , T(1) - ( _a + _b + _c ) } )
  { }

  using pl_t::eval ; 
} ;

// struct inverse_lcp calculates the radial scaling factor to convert
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

template < typename T , std::size_t LANES >
struct inverse_lcp
: public zimt::unary_functor < T , T , LANES >
{
  // even though there are only three coefficients a, b, and c for
  // the mens correction polynomial, the function to yield an
  // appropriately scaled radius for an unscaled one is a fourth-degree
  // polynomial. The fourth coefficient is one minus the sum of the
  // other coefficients, to make the polynomial evaluate to 1.0 at
  // x = 1.0. The last coefficient has to be zero, So that f(0) == 0

  double coefficients[5] ;

  // The scaling factor scales incoming radii in [0...r_max] to
  // [0...16], to be used to evaluate the b-spline

  double r_max ;
  double scale ;

  // b-spline and evaluator use type T - usually float - because they
  // are employed in rendering jobs. When working CPs, use double as T.

  zimt::bspline < T , 1 > inv_model ;
  zimt::grok_type < T , T , LANES > fev ;

  // note how we create a b-spline which has several control points
  // outside of the area in which we'll evaluate it: this is to avoid
  // margin effects. We use a cubic spline, which combines good precision
  // and reasonable processing times.

  inverse_lcp ( double _a , double _b , double _c , double _r_max )
  : coefficients { _a , _b , _c , 1.0 - ( _a + _b + _c ) , 0.0 } ,
    inv_model ( 16 , 3 , zimt::NATURAL ) ,
    r_max ( _r_max ) ,
    scale ( 16.0 / _r_max )
  {
    eu_polynomial < double , 4 , LANES > p ( coefficients ) ;
    for ( std::size_t i = 0 ; i < inv_model.core.shape[0] ; i++ )
    {
      // now we set the spline's control points, starting further
      // to the left to avoid margin effects, and proceeding
      // some steps beyond the coefficient corresponding to the
      // radius at the corner of the image

      double notch = i / scale ;
      double out = notch ; // reasonable starting point of iteration
      p.reval ( notch , out ) ; // will throw if inverse isn't found

      // as the b-spline's knot points, we don't store the radius
      // but rather the factor we need to apply to the incoming radius
      // to obtain the outgoing radius.

      inv_model.core [ i ] = (    notch == 0.0
                                ? 1.0 / p.derivative ( 0.0 )
                                : out / notch ) ;

      std::cout << "notch " << i << " rr " << notch
                << " r " << out << " f " << inv_model.core [ i ]
                << std::endl ;
    }
    // new we set up an evaluator for the spline.

    inv_model.prefilter() ;
    fev = zimt::make_safe_evaluator ( inv_model ) ;
  }

  // the eval function scales and shifts incoming radius values, which
  // sould be in [0,r_max], to evaluate in the central half of the spline.
  // output is the desired scaling factor, which can be used directly
  // to scale centered 2D coordinates.

  template < typename I , typename O >
  void eval ( const I & in , O & out )
  {
    fev.eval ( in * scale , out ) ;
  }
} ;

END_ZIMT_SIMD_NAMESPACE
HWY_AFTER_NAMESPACE() ;

#endif // sentinel
