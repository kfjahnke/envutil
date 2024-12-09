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

// test program for the projections in geometry.h and the
// 'steppers' in stepper.h
// TODO: test rotated steppers agains 'plain' conversions
// followed by a rotation.

#include <random>
#include "stepper.h"

typedef zimt::xel_t < double , 2 > d2_t ;
typedef zimt::xel_t < double , 3 > d3_t ;
typedef zimt::simdized_type < d3_t , LANES > d3_v ;

template < typename pa , typename pb >
d3_t work ( d3_t _c3 )
{
  // set up a ray-to-ray functor containing the two functors
  // of type 'pa' and 'pb'. The ray-to-ray functor takes a
  // 3D ray coordinate, converts it to a 2D coordinate with
  // the projection codified in pa and then back to a ray
  // with the projection codified in pb. Here, we test pa
  // and pb using the same projection - in actual code, pa
  // and pb might use different projections.

  pa tf1 ;
  pb tf2 ;
  ray_to_ray < double , LANES > tf ( tf1 , tf2 ) ;

  // we want to make sure that the pa and pb functors do
  // their job with both scalar and simdized values, so we
  // run the test also with a simdized evaluation.

  d3_v cv ( _c3 ) ;
  d3_t c3 ;
  tf.eval ( _c3 , c3 ) ;
  tf.eval ( cv , cv ) ;

  // we test that the first lane of the vectorized output
  // is very close to the scalar output. With some back-ends,
  // the values are actually equal, but we can't rely on it.

  assert ( std::abs ( cv[0][0] - c3[0] ) < .0000000000001 ) ;
  assert ( std::abs ( cv[1][0] - c3[1] ) < .0000000000001 ) ;
  assert ( std::abs ( cv[2][0] - c3[2] ) < .0000000000001 ) ;

  // we compare input and output and make sure that they
  // differ if at all then only minimally. Note how we compare
  // the values after normalization - they may differ by a
  // factor, because the code does not necessarily produce
  // normalized output - nor is the input normalized.
  
  assert ( std::abs ( _c3[0] / norm ( _c3 ) - c3[0] / norm ( c3 ) )
           < .0000000000001 ) ;

  assert ( std::abs ( _c3[1] / norm ( _c3 ) - c3[1] / norm ( c3 ) )
           < .0000000000001 ) ;

  assert ( std::abs ( _c3[2] / norm ( _c3 ) - c3[2] / norm ( c3 ) )
           < .0000000000001 ) ;

  return c3 ;
}

// we use a two-level dispatch to get from run-time values for
// the projections to types of conversion functors

template < typename pa >
d3_t dispatch2 ( projection_t pb ,
                 d3_t c3 )
{
  d3_t result ;

  switch ( pb )
  {
    case SPHERICAL:
      result = work < pa , ll_to_ray_t < double > > ( c3 ) ;
      break ;
    case CYLINDRICAL:
      result = work < pa , cyl_to_ray_t < double > > ( c3 ) ;
      break ;
    case RECTILINEAR:
      result = work < pa , rect_to_ray_t < double > > ( c3 ) ;
      break ;
    case FISHEYE:
      result = work < pa , fish_to_ray_t < double > > ( c3 ) ;
      break ;
    case STEREOGRAPHIC:
      result = work < pa , ster_to_ray_t < double > > ( c3 ) ;
      break ;
    case CUBEMAP:
      result = work < pa , ir_to_ray_t < double > > ( c3 ) ;
      break ;
    default:
      break ;
  }

  return result ;
}

d3_t dispatch ( projection_t pa ,
                projection_t pb ,
                d3_t c3 )
{
  d3_t result ;

  switch ( pa )
  {
    case SPHERICAL:
      result = dispatch2 < ray_to_ll_t < double > > ( pb , c3 ) ;
      break ;
    case CYLINDRICAL:
      result = dispatch2 < ray_to_cyl_t < double > > ( pb , c3 ) ;
      break ;
    case RECTILINEAR:
      // negative z axis doesn't work with rectilinear projection, hence:
      c3[2] = std::abs ( c3[2] ) ;

      result = dispatch2 < ray_to_rect_t < double > > ( pb , c3 ) ;
      break ;
    case FISHEYE:
      result = dispatch2 < ray_to_fish_t < double > > ( pb , c3 ) ;
      break ;
    case STEREOGRAPHIC:
      result = dispatch2 < ray_to_ster_t < double > > ( pb , c3 ) ;
      break ;
    case CUBEMAP:
      result = dispatch2 < ray_to_ir_t < double > > ( pb , c3 ) ;
      break ;
    default:
      break ;
  }

  return result ;
}

void test_r2r ( d3_t d3 )
{
  d3_t result ;

  // for this test, we set up a ray_to_... functor as first
  // functor, then a ..._to_ray functor as second functor.
  // The projection is the same for the two functors, but
  // the direction of the transformation (3D->2D vs. 2D->3D)
  // is opposite.

  for ( int pa = SPHERICAL ; pa <= CUBEMAP ; pa++ )
  {
    int pb = pa ;
    {
      result = dispatch ( projection_t(pa) ,
                          projection_t(pb) , d3 ) ;
    }
  }
}

int main ( int argc , char * argv[] )
{
  // first, we run a test where we project rays to a planar surface
  // and back to the ray. This should always succeed (except for
  // rectilinear projection and negative z axis, which is avoided)
  // because the resulting 2D coordinate can always regenerate the
  // ray precisely, whereas the opposite operation (start with a
  // planar coordinate, move to ray and back) will fail for certain
  // planar coordinates - e.g. with spherical projection and y==0
  // where the initial x coordinate can't be recovered.

  std::mt19937 gen ; // Standard mersenne_twister_engine
  std::uniform_real_distribution<> dis(-10.0, 10.0);

  for ( std::size_t i = 0 ; i < 10000 ; i++ )
  {
    test_r2r ( { dis(gen) , dis(gen) , dis(gen) } ) ;
  }

  // Next, we populate arrays of ray coordinates with 'steppers'
  // and test the result against arrays which are populated using
  // the ..._to_ray_t functors. We also test against values gleaned
  // by directly invoking the functors' eval function with input
  // generated by a lambda directly from discrete coordinates.

  // two equally-shaped target arrays take up the output of
  // the stepper and the output from the plain conversion
  // operator, which should be identical.

  zimt::array_t < 2 , d3_t > a3 ( { 1000 , 500 } ) ;
  zimt::array_t < 2 , d3_t > b3 ( { 1000 , 500 } ) ;

  // first set up the stepper.

  spherical_stepper < double , LANES >
    sphs ( { 1.0 , 0.0 , 0.0 } ,
           { 0.0 , 1.0 , 0.0 } ,
           { 0.0 , 0.0 , 1.0 } ,
           1000 ,
           500 ) ;

  // use a zimt::process run to fill the first target array with
  // output from the stepper. We use a pass_through act functor
  // because we're only interested in the stepper's output.

  zimt::process ( a3.shape , sphs ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  // now we set up parameters for the use of the conversion
  // operator. The 2D coordinates - in model space units - run
  // from (x0, y0) with a step of dx in the horizontal and dy
  // in the vertical. Note how we're working with centered
  // coordinates, as we would in a viewing context, where
  // the images are draped in model space so that their center
  // rides on the forward axis.

  double x0 = ( a3.shape[0] - 1 ) / -2.0 ;
  double y0 = ( a3.shape[1] - 1 ) / -2.0 ;
  double dx = 2.0 * M_PI / a3.shape[0] ;
  double dy = M_PI / a3.shape[1] ;
  x0 *= dx ;
  y0 *= dy ;

  // this is the conversion operator:

  ll_to_ray_t < double , LANES > ll_to_ray ;

  // we set up a zimt linspace generator as input

  zimt::linspace_t < double , 2 , 2 , LANES >
    ls ( { x0 , y0 } , { dx , dy } ) ;

  // and use zimt::process to populate the second target array

  zimt::process ( b3.shape , ls , ll_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  // for the 'manual' test we set up a lambda which provides
  // centered 2D model space coordinates for discrete coordinates

  auto to_md = [&] ( std::size_t sx , std::size_t sy )
  {
    zimt::xel_t < double , 2 > crd2 ;
    double & x ( crd2[0] ) ;
    double & y ( crd2[1] ) ;
    x = sx - ( a3.shape[0] - 1 ) / 2.0 ;
    y = sy - ( a3.shape[1] - 1 ) / 2.0 ;
    x *= 2.0 * M_PI / a3.shape[0] ;
    y *= M_PI / a3.shape[1] ;
    return crd2 ;
  } ;

  // now we feed the converter in a loop. We use the discrete
  // coordinates to test the result values in arrays a3 and b3
  // against the results we receive as we go.

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto crd2 = to_md ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      ll_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  // we repeat the process for all projections.

  cylindrical_stepper < double , LANES >
    cyls ( { 1.0 , 0.0 , 0.0 } ,
           { 0.0 , 1.0 , 0.0 } ,
           { 0.0 , 0.0 , 1.0 } ,
           1000 , 500 ,
           -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  zimt::process ( a3.shape , cyls ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  cyl_to_ray_t < double , LANES > cyl_to_ray ;

  zimt::process ( b3.shape , ls ,
                  cyl_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto crd2 = to_md ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      cyl_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  rectilinear_stepper < double , LANES >
    rects ( { 1.0 , 0.0 , 0.0 } ,
            { 0.0 , 1.0 , 0.0 } ,
            { 0.0 , 0.0 , 1.0 } ,
            1000 , 500 ,
            -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  zimt::process ( a3.shape , rects ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  rect_to_ray_t < double , LANES > rect_to_ray ;

  zimt::process ( b3.shape , ls ,
                  rect_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto crd2 = to_md ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      rect_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  fisheye_stepper < double , LANES >
    fishs ( { 1.0 , 0.0 , 0.0 } ,
            { 0.0 , 1.0 , 0.0 } ,
            { 0.0 , 0.0 , 1.0 } ,
            1000 , 500 ,
            -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  zimt::process ( a3.shape , fishs ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  fish_to_ray_t < double , LANES > fish_to_ray ;

  zimt::process ( b3.shape , ls ,
                  fish_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto crd2 = to_md ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      fish_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  stereographic_stepper < double , LANES >
    sters ( { 1.0 , 0.0 , 0.0 } ,
            { 0.0 , 1.0 , 0.0 } ,
            { 0.0 , 0.0 , 1.0 } ,
            1000 , 500 ,
            -M_PI , M_PI , -M_PI_2 , M_PI_2 ) ;

  zimt::process ( a3.shape , sters ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a3 ) ) ;

  ster_to_ray_t < double , LANES > ster_to_ray ;

  zimt::process ( b3.shape , ls ,
                  ster_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b3 ) ) ;

  for ( std::size_t y = 0 ; y < a3.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a3.shape[0] ; x++ )
    {
      auto crd2 = to_md ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      ster_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b3 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  // for cubemap projection, we use differently-shaped
  // arrays - using a shpe which would be natural for a
  // 'standard' 1:6 vertically-stacked cubemap. Apart from
  // parameterization, the process is the same as for the
  // other projections.

  zimt::array_t < 2 , d3_t > a6 ( { 500 , 3000 } ) ;
  zimt::array_t < 2 , d3_t > b6 ( { 500 , 3000 } ) ;

  cubemap_stepper < double , LANES >
    cbms ( { 1.0 , 0.0 , 0.0 } ,
           { 0.0 , 1.0 , 0.0 } ,
           { 0.0 , 0.0 , 1.0 } ,
           500 , 3000 ,
           -1.0 , 1.0 , -6.0 , 6.0 ) ;

  zimt::process ( a6.shape , cbms ,
                  zimt::pass_through < double , 3 , LANES > () ,
                  zimt::storer < double , 3 , 2 , LANES > ( a6 ) ) ;

  ir_to_ray_t < double , LANES > ir_to_ray ;

  x0 = ( a6.shape[0] - 1 ) / -2.0 ;
  y0 = ( a6.shape[1] - 1 ) / -2.0 ;
  dx = 2.0 / a6.shape[0] ;
  dy = 12.0 / a6.shape[1] ;
  x0 *= dx ;
  y0 *= dy ;

  zimt::linspace_t < double , 2 , 2 , LANES >
    ls2 ( { x0 , y0 } , { dx , dy } ) ;
           
  zimt::process ( b6.shape , ls2 ,
                  ir_to_ray ,
                  zimt::storer < double , 3 , 2 , LANES > ( b6 ) ) ;

  auto to_md2 = [&] ( std::size_t sx , std::size_t sy )
  {
    zimt::xel_t < double , 2 > crd2 ;
    double & x ( crd2[0] ) ;
    double & y ( crd2[1] ) ;
    x = sx - ( a6.shape[0] - 1 ) / 2.0 ;
    y = sy - ( a6.shape[1] - 1 ) / 2.0 ;
    x *= 2.0 / a6.shape[0] ;
    y *= 12.0 / a6.shape[1] ;
    return crd2 ;
  } ;

  for ( std::size_t y = 0 ; y < a6.shape[1] ; y++ )
  {
    for ( std::size_t x = 0 ; x < a6.shape[0] ; x++ )
    {
      auto crd2 = to_md2 ( x , y ) ;
      zimt::xel_t < double , 3 > crd3 ;
      ir_to_ray.eval ( crd2 , crd3 ) ;
      auto d = abs ( crd3 - a6 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
      d = abs ( crd3 - b6 [ { x , y } ] ) . sum() ;
      assert ( d < .0000000000001 ) ;
    }
  }

  // if none of the assertions failed and the program terminates
  // now, all is well and the functors work as expected.
}
