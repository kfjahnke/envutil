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

// the classes in this header are used to create 'mask' images. They are
// made by producing pixels in a given set colour rather than the greyscale
// or RGB component contained in an image. For single-image processing this
// isn't very useful, it's used for multi-facet jobs: one facet is set to
// produce 'white' pixels, the other facets are set to produce 'black'
// pixels. When subsequent processing produces output, it will be 'white'
// only where greyscale/RGB content from the singled-out facet would be
// visible in an 'ordinary' job. This mask can be used as spatial 'quality'
// indicator in the Burt&Adelson image splining algorithm, which is why I
// have implemented the feature - the B&A algorithm is not yet available
// in envutil, but the masking jobs are here already, so the masks plus
// 'solo' images can be fed to other programs which can do the job.
// images which have an alpha channel will result in masks with an alpha
// channel, where transparent sections of the image will result in
// transparent sections in the mask image. This may be unwanted - to
// create a single-channel mask, the result can be 'cast down' to one
// channel.

#if defined(ENVUTIL_MASKING_H) == defined(HWY_TARGET_TOGGLE)
  #ifdef ENVUTIL_MASKING_H
    #undef ENVUTIL_MASKING_H
  #else
    #define ENVUTIL_MASKING_H
  #endif

#include "zimt/bspline.h"
#include "zimt/eval.h"

HWY_BEFORE_NAMESPACE() ;
BEGIN_ZIMT_SIMD_NAMESPACE(project)

// the masking_t functor yields a pixel with C channels painted
// 'paint' unconditionally.

template < std::size_t D , std::size_t C , std::size_t L >
struct masking_t
: public zimt::unary_functor < zimt::xel_t < float , D > ,
                               zimt::xel_t < float , C > ,
                               L >
{
  const float paint ;

  masking_t ( const float & _paint )
  : paint ( _paint )
  { }

  template < typename I , typename O >
  void eval ( const I & in , O & out )
  {
    out = paint ;
  }
} ;

template < std::size_t C , std::size_t L >
struct alpha_masking_t
: public zimt::unary_functor < zimt::xel_t < float , 2 > ,
                               zimt::xel_t < float , C > ,
                               L >
{
  const float paint ;
  typedef zimt::xel_t < float , C > px_t ;
  typedef zimt::bspline < px_t , 2 > spl_t ;
  typedef zimt::bspline < float , 2 > spl1_t ;
  typedef simdized_type < float , L > f_v ;

  std::shared_ptr < spl_t > p_bspl ;
  grok_type < zimt::xel_t < float , 2 > ,
              zimt::xel_t < float , C > ,
              L > aev ;

  alpha_masking_t ( const float & _paint ,
                    std::shared_ptr < spl_t > _p_bspl )
  : paint ( _paint ) ,
    p_bspl ( _p_bspl ) ,
    aev ( make_safe_evaluator < spl_t , float , L >
           ( * _p_bspl ) )
  { }

  template < typename I , typename O >
  void eval ( const I & in , O & out )
  {
    // TODO: a bit wasteful (unless the optimizer 'gets it')
    // - I'd rather use a channel_view_type, but I can't get
    // the code to compile...

    aev.eval ( in , out ) ;
    out [ 0 ] = paint * out [ C - 1 ] ;
    if constexpr ( C == 4 )
    {
      out [ 1 ] = out [ 2 ] = out[0] ;
    }
  }
} ;

END_ZIMT_SIMD_NAMESPACE
HWY_AFTER_NAMESPACE() ;

#endif // sentinel

