/************************************************************************/
/*                                                                      */
/*   utility to convert and extract images from 360 degree environments */
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

// this file handles the dispatch to ISA-specific code. This only makes
// sense for code which benefits from SIMD, the remainder of the program,
// which deals with I/O and arguments, can be compiled with a single
// common-denominator ISA.

#include "envutil_dispatch.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "envutil_dispatch.cc" 

#include <hwy/foreach_target.h> 
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE() ;

namespace project
{
  namespace HWY_NAMESPACE
  {
    // there isn't any payload code in this file - the implementations
    // of the ISA-specific derived 'dispatch' classes is coded in
    // envutil_payload.cc. We want to keep definitions of ISA-specific
    // things in the payload TU, so that we can invoke this code
    // for TUs which are made without highway: such TUs are best
    // 'pulled into' the main program if they offer the same handles.

    struct dispatch
    : public dispatch_base
    {
      // We fit the derived dispatch class with a c'tor which fills in
      // information about the nested SIMD ISA we're currently in.

      dispatch() ;

      // This is the declaration of the single 'payload' function:

      int payload ( int nchannels ,
                    int ninputs ,
                    projection_t projection ) const ;
    } ;

    // _get_dispatch is also declared here, because we'll use HWY_EXPORT
    // and HWY_DYNAMIC_DISPATCH in this file

    const dispatch_base * const _get_dispatch() ;

  } ;
} ;

HWY_AFTER_NAMESPACE() ;

#if HWY_ONCE

namespace project
{
  // we're using highway's foreach_target mechanism. To get access to the
  // SIMD-ISA-specific variant of _get_dispatch (in project::HWY_NAMESPACE)
  // we use the HWY_EXPORT macro:

  HWY_EXPORT ( _get_dispatch ) ;

  // now we can code get_dispatch: it simply uses HWY_DYNAMIC_DISPATCH
  // to pick the SIMD-ISA-specific get_dispatch variant, which in turn
  // yields the desired dispatch_base pointer.
  // The main program calls 'get_dispatch' and receives a dispatch_base
  // pointer, which can be used to route to ISA-specific code best suited
  // for the CPU running the code:
  // auto dp = get_dispatch() ;
  // dp->payload ( ... ) ;

  const dispatch_base * const get_dispatch()
  {
    return HWY_DYNAMIC_DISPATCH ( _get_dispatch ) () ;
  }

} ;

#endif
