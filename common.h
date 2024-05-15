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

// this header has code common to the entire program, like common types.

#include "zimt/zimt.h"

// zimt types for 2D and 3D coordinates and pixels

typedef zimt::xel_t < int , 2 > v2i_t ;
typedef zimt::xel_t < int , 3 > v3i_t ;
typedef zimt::xel_t < long , 2 > index_type ;
typedef zimt::xel_t < std::size_t , 2 > shape_type ;

typedef zimt::xel_t < float , 2 > v2_t ;
typedef zimt::xel_t < float , 3 > v3_t ;

// some SIMDized types we'll use. I use 16 SIMD lanes for now,
// which is also the lane count currently supported by OIIO.

// #ifdef OIIO_TEXTURE_SIMD_BATCH_WIDTH
// #define LANES OIIO_TEXTURE_SIMD_BATCH_WIDTH
// #else
#define LANES 16
// #endif

typedef zimt::simdized_type < float , LANES > f_v ;
typedef zimt::simdized_type < int , LANES > i_v ;
typedef zimt::simdized_type < v2_t ,  LANES > crd2_v ;
typedef zimt::simdized_type < v2i_t , LANES > v2i_v ;
typedef zimt::simdized_type < v3i_t , LANES > v3i_v ;
typedef zimt::simdized_type < v3_t ,  LANES > crd3_v ;
typedef zimt::simdized_type < index_type , LANES > index_v ;

