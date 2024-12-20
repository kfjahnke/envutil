/************************************************************************/
/*                                                                      */
/*    zimt - abstraction layer for SIMD programming                     */
/*                                                                      */
/*            Copyright 2024 by Kay F. Jahnke                           */
/*                                                                      */
/*    The git repository for this software is at                        */
/*                                                                      */
/*    https://github.com/kfjahnke/zimt                                  */
/*                                                                      */
/*    Please direct questions, bug reports, and contributions to        */
/*                                                                      */
/*    kfjahnke+zimt@gmail.com                                           */
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

// This is a test program for zimt's recently acquired b-spline
// processing capabilites and also serves to measure performance of the
// b-spline evaluation code with splines of various degrees and boundary
// conditions and varying SIMD back-ends/ISAs.

// This variant (of geometry.cc) uses separate SIMD_ISA-sepecific TUs
// for the 'payload' code which are linked in. This allows for more
// freedom in compiling those payload TUs. See 'inset.cc' for the code
// for the ISA-specific TUs and for compile instructions. This variant
// is for use with MULTI_SIMD_ISA - the Vc and std::simd back-ends
// can't be used with it to good effect, unless the SIMD-ISA-specific
// TUs are built with SIMD-ISA-specific compiler flags. Here, we
// use highway's 'internal' dispatch mechanism (using #pragmas inside
// the code) to affect the ISA-specific compilation. Of course, the
// general method of obtaining a dispatch_base pointer to some
// derived dispatch object somewhere else is perfectly general, but
// here we want to exploit the highway code base and namespace
// structure as much as we can to make the job as easy as possible
// and to demonstrate how little is needed to get the multi-ISA
// binary even with this more involved method.

// we define a dispatch base class. All the 'payload' code is called
// through virtual member functions of this class. In this example,
// we only have a single payload function. This definition might go
// to a header file, but it's better placed here, because it holds
// the declarations of the pure virtual member function(s) used as
// conduit to the ISA-specific code.
// Because this file is re-included several times, we need a sentinel
// to protect this class definition.

#ifndef DISPATCH_BASE_DEFINED
#define DISPATCH_BASE_DEFINED

struct dispatch_base
{
  // in dispatch_base and derived classes, we keep two flags.
  // 'backend' holds a value indicating which of zimt's back-end
  // libraries is used. 'hwy_isa' is only set when the highway
  // backend is used and holds highway's HWY_TARGET value for
  // the given nested namespace.

  int backend = -1 ;
  unsigned long hwy_isa = 0 ;

  // next we have pure virtual member function definitions for
  // payload code. In this example, we only have one payload
  // function which calls what would be 'main' in a simple
  // program without multiple SIMD ISAs or SIMD back-ends

  virtual int payload ( int argc , char * argv[] ) const = 0 ;
} ;

#endif

// we're using MULTI_SIMD_ISA, so we have to define HWY_TARGET_INCLUDE
// to tell the foreach_target mechanism which file should be repeatedly
// re-included and re-copmpiled with SIMD-ISA-specific flags

#undef HWY_TARGET_INCLUDE

/////////////// Tell highway which file to submit to foreach_target

#define HWY_TARGET_INCLUDE "disp_to_tu.cc"  // this very file

//--------------------------------------------------------------------

#include <hwy/foreach_target.h>  // must come before highway.h

#include <hwy/highway.h>

//////////////// Put the #includes needed for your program here:

#include <iostream>
#include "zimt/zimt.h"

//--------------------------------------------------------------------

// to make highway's use of #pragma directives to the compiler
// effective, we surround the SIMD-ISA-specific code with
// HWY_BEFORE_NAMESPACE() and HWY_AFTER_NAMESPACE().

HWY_BEFORE_NAMESPACE() ;

// this macro puts us into a nested namespace inside namespace 'project'.
// For multi-SIMD-ISA builds it is project::HWY_NAMESPACE. The macro
// is defined in common.h. After the macro invocation, we can access
// all zimt names with a simple zimt:: prefix - both 'true' zimt names
// and SIMD-ISA-specific versions living in the nested namespace.

BEGIN_ZIMT_SIMD_NAMESPACE(project)

// we also use a local function _get_dispatch which returns a pointer
// to 'dispatch_base', which points to an object of the derived class
// 'dispatch'. This is used with highway's HWY_DYNAMIC_DISPATCH and
// returns the dispatch pointer for the SIMD ISA which highway deems
// most appropriate for the CPU on which the binary is executing.
// Note how this function is actually defined in the ISA-specific TUs,
// it's enough that we provide an external declaration here.

extern const dispatch_base * const _get_dispatch() ;

// that's us done with the ISA-specific nested namespace and ISA-specific
// compilation.

END_ZIMT_SIMD_NAMESPACE

HWY_AFTER_NAMESPACE() ;

// Now for code which isn't SIMD-ISA-specific.

#if HWY_ONCE

namespace project {

// we're using highway's foreach_target mechanism. To get access to the
// SIMD-ISA-specific variant of _get_dispatch (in project::HWY_NAMESPACE)
// we use the HWY_EXPORT macro:

HWY_EXPORT(_get_dispatch);

// now we can code get_dispatch: it simply uses HWY_DYNAMIC_DISPATCH
// to pick the SIMD-ISA-specific get_dispatch variant, which in turn
// yields the desired dispatch_base pointer.

const dispatch_base * const get_dispatch()
{
  return HWY_DYNAMIC_DISPATCH(_get_dispatch)() ;
}

}  // namespace project

int main ( int argc , char * argv[] )
{
  // Here we use zimt's dispatch mechanism: first, we get a pointer
  // to the dispatcher, then we invoke a member function of the
  // dispatcher. What's the point? We can call a SIMD-ISA-specific
  // bit of code without having to concern ourselves with figuring
  // out which SIMD ISA to use on the current CPU: this happens via
  // highway's dispatch mechanism, or is fixed at compile time, but
  // in any case we receive a dispatch_base pointer routing to the
  // concrete variant. project::get_dispatch might even be coded
  // to provide pointers to dispatch objects in separate TUs, e.g.
  // when these TUs use different back-ends or compiler flags. Here,
  // we can remain unaware of how the concrete dispatch object is
  // set up and the pointer obtained.

  auto dp = project::get_dispatch() ;

  // we can get information about the specific dispatch object:

  std::cout << "calling payload with "
#ifdef MULTI_SIMD_ISA
            << "dynamic dispatch" << std::endl ;
#else
            << "static dispatch" << std::endl ;
#endif

  // now we call the payload via the dispatch_base pointer.

  int success = dp->payload ( argc , argv ) ;

  std::cout << "payload returned " << success << std::endl ;

  //to invoke other variants with explicit calls, uncomment
  // this section (x86 only - other CPU types would need other
  // N_... macros)

#ifdef MULTI_SIMD_ISA

  std::cout << "******** calling payloads with successively better ISAs:"
            << std::endl ;
  
  std::cout << "******** calling SSE2 payload" << std::endl ;
  dp = project::N_SSE2::_get_dispatch() ;
  dp->payload ( argc , argv ) ;
  
  std::cout << "******** calling SSSE3 payload" << std::endl ;
  dp = project::N_SSSE3::_get_dispatch() ;
  dp->payload ( argc , argv ) ;
  
  std::cout << "******** calling SSE4 payload" << std::endl ;
  dp = project::N_SSE4::_get_dispatch() ;
  dp->payload ( argc , argv ) ;
  
  std::cout << "******** calling AVX2 payload" << std::endl ;
  dp = project::N_AVX2::_get_dispatch() ;
  dp->payload ( argc , argv ) ;
  
  std::cout << "******** calling AVX3 payload" << std::endl ;
  dp = project::N_AVX3::_get_dispatch() ;
  dp->payload ( argc , argv ) ;

#endif
}

#endif  // ZIMT_ONCE

