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

#if defined(ENVUTIL_ENVIRONMENT_H) == defined(HWY_TARGET_TOGGLE)
  #ifdef ENVUTIL_ENVIRONMENT_H
    #undef ENVUTIL_ENVIRONMENT_H
  #else
    #define ENVUTIL_ENVIRONMENT_H
  #endif

#include <filesystem>
#include <OpenImageIO/imageio.h>

#include "zimt/prefilter.h"
#include "zimt/bspline.h"
#include "zimt/eval.h"
#include "zimt/convolve.h"

#include "twining.h"
#include "geometry.h"
#include "cubemap.h"
#include "masking.h"
#include "lens_correction.h"

HWY_BEFORE_NAMESPACE() ;
BEGIN_ZIMT_SIMD_NAMESPACE(project)

using OIIO::ImageInput ;
using OIIO::TypeDesc ;
using OIIO::ImageSpec ;

// source_t provides pixel values from a mounted 2D manifold. The
// incoming 2D coordinates are pixel coordinates guaranteed to be
// in the range of the spline:
// ( -0.5 ... width - 0.5 ) , ( -0.5 ... height - 0.5 )

template < std::size_t nchannels , std::size_t ncrd , std::size_t L >
struct source_t
: public zimt::unary_functor < zimt::xel_t < float , ncrd > ,
                               zimt::xel_t < float , nchannels > ,
                               L
                             >
{
  typedef zimt::unary_functor < zimt::xel_t < float , ncrd > ,
                                zimt::xel_t < float , nchannels > ,
                                L
                              > base_t ;

  typedef zimt::xel_t < float , ncrd > pkg_t ;
  typedef zimt::simdized_type < pkg_t , L > pkg_v ;
  typedef zimt::xel_t < float , 2 > crd_t ;
  typedef zimt::simdized_type < crd_t , L > crd_v ;
  typedef zimt::xel_t < float , nchannels > px_t ;
  typedef zimt::simdized_type < px_t , L > px_v ;

  static_assert ( ncrd == 2 || ncrd == 6 ) ;

  int width ;
  int height ;
  zimt::xel_t < int , 2 > strides ;

  typedef zimt::bspline < px_t , 2 > spl_t ;
  std::shared_ptr < spl_t > p_bspl ;
  zimt::grok_type < crd_t , px_t , L > bsp_ev ;

  // we may receive a nullptr in p_bspl. this can occur if there is
  // a 'solo' argument, in which case loading image data for the other
  // facets is futile. This will leave bsp_ev in it's default-c'ted
  // state, which is okay, since it won't be called.

  source_t ( std::shared_ptr < spl_t > _p_bspl , int masked )
  : p_bspl ( _p_bspl ) ,
    width ( _p_bspl ? _p_bspl->core.shape[0] : 1 ) ,
    height ( _p_bspl ? _p_bspl->core.shape[1] : 1 )
  {
    if ( masked == -1 )
    {
      if ( p_bspl )
        bsp_ev = make_safe_evaluator < spl_t , float , L > ( *p_bspl ) ;
    }
    else if ( nchannels == 1 || nchannels == 3 )
    {
      bsp_ev = masking_t < 2 , nchannels , L > ( masked ) ;
    }
    else
    {
      if ( p_bspl )
        bsp_ev = alpha_masking_t < nchannels , L >
                 ( masked , p_bspl ) ;
    }
  }

  void eval ( const pkg_v & _crd , px_v & px )
  {
    // incoming is const &, hence:

    pkg_v crd ( _crd ) ;

    if constexpr ( ncrd == 2 )
    {
      // crd[0] *= width ;
      // crd[1] *= height ;
      // crd -= .5f ;
      bsp_ev.eval ( crd , px ) ;
    }
    else
    {
      assert ( false ) ;
    }
  }
} ;

// struct mount_t provides data from a rectangular 2D manifold holding
// pixel data in a given projection which may not cover the entire
// 360X180 degree environment. Typical candidates would be rectilinear
// patches or cropped images. This class handles channel counts up to
// four and paints pixels outside the covered range black. For RGBA
// pixels, 'outside' pixels are painted 0000, assuming associated alpha.
// The 'ncrd' template argument is either three for lookups without
// or nine for lookups with two next neighbour's coordinates which
// can be used to find derivatives. This class relies on an 'inner'
// functor of class source_t to actually provide pixel data, this here
// class deals with the geometrical transformations needed for the
// different types of projections which the mounted images can have,
// and with masking out parts where the mounted image has no data.

template < std::size_t nchannels ,
           projection_t P ,
           std::size_t L >
struct mount_t
: public zimt::unary_functor < zimt::xel_t < float , 3 > ,
                               zimt::xel_t < float , nchannels > ,
                               L
                             >
{
  typedef zimt::unary_functor < zimt::xel_t < float , 3 > ,
                                zimt::xel_t < float , nchannels > ,
                                L
                              > base_t ;

  typedef zimt::xel_t < float , 3 > ray_t ;
  typedef zimt::simdized_type < ray_t , L > ray_v ;
  typedef zimt::xel_t < float , 2 > crd_t ;
  typedef zimt::xel_t < float , 3 > crd3_t ;
  typedef zimt::simdized_type < crd_t , L > crd_v ;
  typedef zimt::simdized_type < crd3_t , L > crd3_v ;
  typedef zimt::xel_t < float , nchannels > px_t ;
  typedef zimt::simdized_type < px_t , L > px_v ;

  using typename base_t::in_v ;
  using typename base_t::in_ele_v ;
  typedef typename in_ele_v::mask_type mask_t ;
  typedef std::function < void ( crd_v & ) > planar_f ;

  extent_type extent ;
  crd_t center ;
  crd_t rgirth ;
  source_t < nchannels , 2 , L > inner ;
  planar_f pf ;

  static_assert (    P == SPHERICAL
                  || P == CYLINDRICAL
                  || P == RECTILINEAR
                  || P == STEREOGRAPHIC
                  || P == FISHEYE ) ;

  static crd_t get_center
    ( const source_t < nchannels , 2 , L > & src ,
      const extent_type & ext )
  {
    // first get the center in model space coordinates

    crd_t center { ( ext.x0 + ext.x1 ) * .5 ,
                   ( ext.y0 + ext.y1 ) * .5 } ;

    // then the distance from the left/upper margin of the extent

    crd_t dc { center[0] - ext.x0 ,
               center[1] - ext.y0 } ;

    // now move to image coordinates

    dc[0] *= src.width ;
    dc[1] *= src.height ;
    dc[0] -= .5 ;
    dc[1] -= .5 ;

    return dc ;
  }
                  
  // calculate the scaling factor from model space to pixel units

  static crd_t get_rgirth
    ( const source_t < nchannels , 2 , L > & src ,
      const extent_type & ext )
  {
    crd_t rgirth { src.width / ( ext.x1 - ext.x0 ) ,
                   src.height / ( ext.y1 - ext.y0 ) } ;
    return rgirth ;
  }
                            
  mount_t ( extent_type _extent ,
            source_t < nchannels , 2 , L > & _inner ,
            const planar_f & _pf = [] ( crd_v & ) { } )
  : extent ( _extent ) ,
    center ( get_center ( _inner , _extent ) ) ,
    rgirth ( get_rgirth ( _inner , _extent ) ) ,
    inner ( _inner ) ,
    pf ( _pf )
  { }

  // shade provides a pixel value for a 2D coordinate inside the
  // 2D manifold's 'extent' by delegating to the 'inner' functor
  // of class source_t

  px_v shade ( crd_v crd )
  {
    // move to image coordinates. The eval code uses clamping,
    // so absolute precision isn't required.

    crd[0] = ( crd[0] - extent.x0 ) * rgirth[0] - .5f ;
    crd[1] = ( crd[1] - extent.y0 ) * rgirth[1] - .5f ;

    // use 'inner' to provide the pixel

    px_v result ;

    inner.eval ( crd , result ) ;
    return result ;
  }

  // plain coordinate transformation without masking. The mask is
  // calculated one function down. The geometrical transformations
  // are coded in 'geometry.h'.
  
  void get_coordinate_nomask ( const crd3_v & crd3 , crd_v & crd ) const
  {
    // depending on the source image's geometry, we convert the
    // incoming 3D ray coordinate to a planar coordinate pertaining
    // to the source image plane

    if constexpr ( P == RECTILINEAR )
    {
      ray_to_rect_t<float,L>::eval ( crd3 , crd ) ;
    }
    else if constexpr ( P == SPHERICAL )
    {
      ray_to_ll_t<float,L>::eval ( crd3 , crd ) ;
    }
    else if constexpr ( P == CYLINDRICAL )
    {
      ray_to_cyl_t<float,L>::eval ( crd3 , crd ) ;
    }
    else if constexpr ( P == STEREOGRAPHIC )
    {
      ray_to_ster_t<float,L>::eval ( crd3 , crd ) ;
    }
    else if constexpr ( P == FISHEYE )
    {
      ray_to_fish_t<float,L>::eval ( crd3 , crd ) ;
    }

    // we now have the 'raw' 2D coordinate in 'crd'. The last step
    // applies an optional in-plane transformation, which is used
    // for stuff like lens correction and shear.

    pf ( crd ) ;
  }

  // get_coordinate yields a coordinate into the mounted 2D manifold
  // and a true mask value if the ray passes through the draped 2D
  // manifold, or a valid coordinate (the center of the 'extent')
  // and a false mask value if the ray does not 'hit' the 2D manifold.

  mask_t get_coordinate ( const crd3_v & crd3 , crd_v & crd ) const
  {
    // first, obtain the 'raw' 2D source coordinate

    get_coordinate_nomask ( crd3 , crd ) ;

    // test this coordinate against the boundaries of the image
    // encoded in 'extent'

    auto mask =    ( crd[0] >= extent.x0 )
                && ( crd[0] <= extent.x1 )
                && ( crd[1] >= extent.y0 )
                && ( crd[1] <= extent.y1 ) ;

    // only for rectilinear images: clear the mask where the z
    // coordinate of the ray is zero or negative

    if constexpr ( P == RECTILINEAR )
    {
      mask &= ( crd3[2] > 0.0f ) ;
    }

    // now, where the mask is not set, assign 'center' to the
    // 2D coordinate - this is a safe value avoiding problems
    // if this value is used e.g. to produce a pixel value.
    // We trust that the mask will prevent such pixels from
    // actually being used. So this is merely a precaution.

    crd ( ! mask ) = center ;

    // finally, we return the *mask* - the resulting 2D
    // coordinate is passed back via the non-const reference
    // we received as 'crd' argument

    return mask ;
  }

  // get_mask does everything which get_coordinate does, minus
  // the actual coordinate transformation. This seems to produce
  // redundant work, but I trust the compiler will recognize the
  // common subexpression.

  mask_t get_mask ( const crd3_v & crd3 ) const
  {
    crd_v crd ;
    get_coordinate_nomask ( crd3 , crd ) ;

    auto mask =    ( crd[0] >= extent.x0 )
                && ( crd[0] <= extent.x1 )
                && ( crd[1] >= extent.y0 )
                && ( crd[1] <= extent.y1 ) ;

    if constexpr ( P == RECTILINEAR )
      mask &= ( crd3[2] > 0.0f ) ;

    return mask ;
  }

  // eval puts it all together and yields pixel sets with 'misses' masked
  // out to zero.

  mask_t eval ( const in_v & ray , px_v & px )
  {
    // the first three components have the 3D pickup coordinate itself

    crd3_v crd3 { ray[0] , ray[1] , ray[2] } ;
    crd_v crd ;

    auto mask = get_coordinate ( crd3 , crd ) ;
    if ( none_of ( mask ) )
    {
      px = 0.0f ;
      return mask ;
    }

    px = shade ( crd ) ;

    // mask out 'misses' to all-zero

    if ( ! all_of ( mask ) )
    {
      px ( ! mask ) = 0.0f ;
    }

    return mask ;
  }
} ;

/// for full spherical images, we have to perform the prefiltering and
/// bracing of the b-splines 'manually', because these images don't fit
/// any of vspline's standard modes. The problem is that while the images
/// are horizontally periodic, the vertical periodicity is not manifest:
/// the vertical periodicity is along great circles passing through the
/// poles, and of these great circles, one half is in the left half and
/// one in the right half of the standard 2:1 spherical format.
/// The prefiltering for the vertical can be done with vspline, but the
/// image halves have to be put into separate views and 'stacked'.
/// This is done first. Next comes the bracing, which also needs to be
/// done manually, since the vertical continuation for one half is always
/// in the other half, going in the opposite direction. This requires
/// more index and view artistry. Coding for full sphericals has become
/// possible only with the introduction of stacked array prefiltering.
/// While the coding presents a certain degree of effort, the reward is
/// mathematically sound treatment of full spherical data, which had been
/// processed with REFLECT BCs after the demise of BC mode
/// zimt::SPHERICAL. So now the specialized code is where it should
/// reside: here in pv_rendering_common.cc, rather than in vspline,
/// where it was misplaced, because it's so specific. And with the new
/// code, we can handle even the smallest possible spherical pamorama,
/// consisting of only two pixels, correctly :D

// TODO: while readymade full sphericals do indeed normally have 2:1
// aspect ratio, the width can be halved and the code functions as
// intended. But odd width can occur, and because it's not meaningless,
// but simply unusual, it should be handled correctly, and it isn't.
// With odd width, we can't split the image into a left and a right
// half-image to stack them on top of each other. It would be necessary
// to generate a stripe on top and a stripe on the bottom which are
// interpolated to yield values just *between* the normal grid columns,
// which would approximate the continuation 'beyond the poles'. This
// could be prefiltered as a stack. To use the fully correct calculation
// which needs periodicity, one would use only one additional stripe
// with full height.

template < typename dtype >
void spherical_prefilter
  ( zimt::bspline < dtype , 2 > & bspl ,
    zimt::view_t < 2 , dtype > & input ,
    int njobs )
{
  auto output = bspl.core ;
  typedef zimt::view_t < 2 , dtype > view_type ;

  // the type of filter we want to use is a b-spline prefilter operating
  // in single-precision float, using Vc::SimdArrays with default vector
  // width for aggregation. Note how this code will use vspline's SIMD
  // emulation code if Vc is not used: this code is *not* part of the pixel
  // pipeline, which has to be reduced to VSIZE 1 if Vc is not present.

  typedef zimt::recursive_filter
          < zimt::simdized_type , float >
    stripe_handler_type ;

  // we set up the specifications for the filter. it is periodic in all
  // cases (the stacked half-images together are periodic vertically)

  int degree = bspl.spline_degree ; // spline degree
  int npoles = degree / 2 ;         // number of filter poles

  zimt::iir_filter_specs specs
    ( zimt::PERIODIC ,
      degree / 2 ,
      zimt_constants::precomputed_poles [ degree ] ,
      0.0001 ) ;

  // we call separable_filter's overload taking MultiArrayViews and an axis
  // to do the horizontal prefiltering step on the unmodified data

  if ( degree > 1 || input.data() != output.data() )
    zimt::detail::separable_filter
      < view_type , view_type , stripe_handler_type >()
        ( input , output , 0 , specs ) ;

  // now we stack the left half and the verticall flipped right half
  // to obtain a stack which is periodic vertically

  std::vector < view_type > vertical_view ;

  // note that for scaled-down spline there is a possibility that the
  // horizontal extent is not even, as it is for the original image.
  // we silently ignore this for now - it may result in a pixel fault.
  // KFJ 2022-01-20 this became an issue. Here's a quick fix for odd
  // width, assuming shape[0], the width, is 'large'. The correct way
  // would be to interpolate at N/2 distance.

  auto shape = output.shape ;
  if ( shape[0] & 1 )
  {
    if ( shape[0] > 1 )
    {
      std::cerr << "warning: input is a full spherical with odd width "
                << shape[0] << std::endl ;
      std::cerr
        << "you may notice slight errors near poles and image boundaries"
        << std::endl ;
    }
    shape [ 0 ] = shape [ 0 ] / 2 + 1 ;
  }
  else
  {
    shape [ 0 ] /= 2 ;
  }

  auto strides = output.strides ;
  strides [ 1 ] = - strides [ 1 ] ;

  view_type upper ( shape , output.strides , output.data() ) ;

  view_type lower ( shape , strides ,
                    output.data()
                    + ( shape[1] - 1 ) * output.strides[1]
                    + shape[0] * output.strides[0] ) ;

  vertical_view.push_back ( upper ) ;
  vertical_view.push_back ( lower ) ;

  // the stack is fed to the overload taking stacks and an axis,
  // but we only need to actually apply the filter if degree is
  // greater than one, otherwise bracing is all that's required.

  if ( degree > 1 )
  {
    zimt::detail::separable_filter
      < view_type , view_type , stripe_handler_type >()
      ( vertical_view , vertical_view , 1 , specs , njobs ) ;
  }

  // now we apply the special bracing for spherical images
  // first get views to the left and right half-image

  view_type left ( shape , output.strides , output.data() ) ;
  view_type right ( shape , output.strides , output.data() + shape[0] ) ;

  // we initialize the indices for the lines we copy from source to target
  // and the limits for these indices to make sure we don't produce memory
  // faults by overshooting. Note how we use views to what amounts to the
  // spline's core, and write to lines outside these views, since we can be
  // sure that this is safe: these writes go to the spline's frame, and we
  // can obtain the frame's size from the bspline object.

  std::ptrdiff_t upper_margin = bspl.left_frame[1] ;
  std::ptrdiff_t y0 = - upper_margin ;
  std::ptrdiff_t lower_margin = bspl.right_frame[1] ;
  std::ptrdiff_t y1 = shape[1] + lower_margin ;

  std::ptrdiff_t upper_source = 0 ;
  std::ptrdiff_t upper_target = -1 ;

  std::ptrdiff_t lower_source = shape[1] - 1 ;
  std::ptrdiff_t lower_target = shape[1] ;

  while ( true )
  {
    // check all indices
    bool check1 = upper_source >= y0 && upper_source < y1 ;
    bool check2 = upper_target >= y0 && upper_target < y1 ;
    bool check3 = lower_source >= y0 && lower_source < y1 ;
    bool check4 = lower_target >= y0 && lower_target < y1 ;

    if ( ( ! check2 ) && ( ! check4 ) )
    {
      // both target indices are out of range, break, we're done
      break ;
    }

    if ( check2 )
    {
      // upper_target is in range
      if ( check1 )
      {
        // so is upper_source, do the copy
        left.slice ( 1 , upper_target )
          . copy_data ( right.slice ( 1 , upper_source ) ) ;
        right.slice ( 1 , upper_target )
          . copy_data ( left.slice ( 1 , upper_source ) ) ;
        upper_source ++ ;
      }
      upper_target -- ;
    }

    if ( check4 )
    {
      // lower_target is in range
      if ( check3 )
      {
        // so is lower_source, do the copy
        left.slice ( 1 , lower_target )
          . copy_data ( right.slice ( 1 , lower_source ) ) ;
        right.slice ( 1 , lower_target )
          . copy_data ( left.slice ( 1 , lower_source ) ) ;
        lower_source -- ;
      }
      lower_target ++ ;
    }
  }

  // vertical bracing is done, apply horizontal bracing

  bspl.brace ( 0 ) ;
  bspl.prefiltered = true ;
}

// this class translates from pixels with N channels to pixels with
// M channels. It's silently assumed that incoming two- and four-
// -channel data carry transparency in the last channel, whereas
// one- and three-channel data are opaque. It's also assumed that
// transparency is given as associated alpha.

template < typename T , std::size_t IN , std::size_t OUT , std::size_t L >
struct repix_t
: public unary_functor < xel_t < T , IN > , xel_t < T , OUT > , L >
{
  typedef unary_functor < xel_t < T , IN > ,
                          xel_t < T , OUT > , L > base_t ;

  using typename base_t::in_v ;
  using typename base_t::out_v ;

  void eval ( const in_v & in , out_v & out ,
              const std::size_t & cap = L ) const
  {
    static_assert ( IN > 0 && OUT > 0 && IN < 5 && OUT < 5 ) ;

    if constexpr ( IN == OUT )
    {
      out = in ;
    }
    else if constexpr ( IN == 1 )
    {
      if constexpr ( OUT == 3 )
      {
        out = in[0] ;
      }
      else
      {
        if constexpr ( OUT == 2 )
        {
          out[0] = in[0] ;
          out[1] = 1.0f ;
        }
        else if constexpr ( OUT == 4 )
        {
          out[0] = in[0] ;
          out[1] = in[0] ;
          out[2] = in[0] ;
          out[3] = 1.0f ;
        }
      }
    }
    else if constexpr ( IN == 2 )
    {
      if constexpr ( OUT == 1 )
      {
        out[0] = in[0] / in[1] ;
        out[0] ( in[1] == 0.0f ) = 0.0f ;
      }
      else if constexpr ( OUT == 3 )
      {
        out[0] = in[0] / in[1] ;
        out[0] ( in[1] == 0.0f ) = 0.0f ;
        out[1] = out[0] ;
        out[2] = out[0] ;
      }
      else if constexpr ( OUT == 4 )
      {
        out[0] = in[0] ;
        out[1] = in[0] ;
        out[2] = in[0] ;
        out[3] = in[1] ;
      }
    }
    else if constexpr ( IN == 3 )
    {
      if constexpr ( OUT == 1 )
      {
        out[0] = in.sum() / 3.0f ;
      }
      else if constexpr ( OUT == 2 )
      {
        out[0] = in.sum() / 3.0f ;
        out[1] = 1.0f ;
      }
      else if constexpr ( OUT == 4 )
      {
        out[0] = in[0] ;
        out[1] = in[1] ;
        out[2] = in[2] ;
        out[3] = 1.0f ;
      }
    }
    else if constexpr ( IN == 4 )
    {
      if constexpr ( OUT == 1 )
      {
        out[0] = ( in[0] + in[1] + in[2] ) / 3.0f ;
        out[0] /= in[3] ;
        out[0] ( in[3] == 0.0f ) = 0.0f ;
      }
      else if constexpr ( OUT == 2 )
      {
        out[0] = ( in[0] + in[1] + in[2] ) / 3.0f ;
        out[1] = in[3] ;
      }
      else if constexpr ( OUT == 3 )
      {
        out[0] = in[0] / in[3] ;
        out[1] = in[1] / in[3] ;
        out[2] = in[2] / in[3] ;
        out ( in[3] == 0.0f ) = 0.0f ;
      }
    }
  }
} ;

// for masking jobs, the R, G and B values are irrelevant. The 'colour'
// is either 1.0 or 0.0, and if there are three 'colour' channels, they
// are equal. So where the class above would form averages, we can simply
// pick any of the three channels. Alpha is preserved and kept associated.
// Only one- and two-chanel output is allowed, that's the point after all.
// Note that we have to propagate the transparency of the facet: here we
// produce the input to the synopsis-forming object, and for jobs with
// alpha channel, this object has to do alpha compositing to create correct
// output. The final product after the synopsis will carry two (equal)
// channels - the 'colour' can only be 0 or 1, and since we're working with
// associated alpha, the first channel is 1 * alpha, so == alpha. Only
// after the synopsis, the alpha channel can be discarded to form a plain
// grey-scale mask.

template < typename T , std::size_t IN , std::size_t OUT , std::size_t L >
struct mono_t
: public unary_functor < xel_t < T , IN > , xel_t < T , OUT > , L >
{
  typedef unary_functor < xel_t < T , IN > ,
                          xel_t < T , OUT > , L > base_t ;

  using typename base_t::in_v ;
  using typename base_t::out_v ;

  void eval ( const in_v & in , out_v & out ,
              const std::size_t & cap = L ) const
  {
    static_assert ( IN > 0 && OUT > 0 && IN < 5 && OUT < 5 ) ;
    assert ( OUT == 1 || OUT == 2 ) ;

    if constexpr ( IN == OUT )
    {
      out = in ;
    }
    else if constexpr ( IN == 1 )
    {
      // OUT must be 2
      out[0] = in[0] ;
      out[1] = 1.0f ;
    }
    else if constexpr ( IN == 2 )
    {
      // OUT must be 1
      out[0] = in[0] / in[1] ;
      out[0] ( in[1] == 0.0f ) = 0.0f ;
    }
    else if constexpr ( IN == 3 )
    {
      if constexpr ( OUT == 1 )
      {
        out[0] = in[0] ;
      }
      else if constexpr ( OUT == 2 )
      {
        out[0] = in[0] ;
        out[1] = 1.0f ;
      }
    }
    else if constexpr ( IN == 4 )
    {
      if constexpr ( OUT == 1 )
      {
        out[0] = in[0] ;
        out[0] /= in[3] ;
        out[0] ( in[3] == 0.0f ) = 0.0f ;
      }
      else if constexpr ( OUT == 2 )
      {
        out[0] = in[0] ;
        out[1] = in[3] ;
      }
    }
  }
} ;

// cubemap_view_t is a stripped-down version of cubemap_t which
// is strictly for viewing. It depends on a ready-made b-spline
// holding the IR image of a cubemap (conveniently made by using
// cubemap_t originally) and a minimal set of parameters and
// code to produce pixel values for ray coordinates. It can also
// produce masks. I prefer to have this class here, because it
// is quite specific and I want to keep details like masking
// (which requires classes masking_t or alpha_masking_t) out of
// the header cubemap.h.

template < std::size_t nchannels , projection_t cbm_prj >
struct cubemap_view_t
: public zimt::unary_functor
   < zimt::xel_t < float , 3 > ,
     zimt::xel_t < float , nchannels > ,
     LANES >
{
  typedef zimt::xel_t < float , 2 > crd2_t ;
  typedef zimt::xel_t < float , 3 > ray_t ;
  typedef zimt::xel_t < float , nchannels > px_t ;

  typedef zimt::simdized_type < int , LANES > i_v ;
  typedef zimt::simdized_type < crd2_t , LANES > crd2_v ;
  typedef zimt::simdized_type < ray_t , LANES > ray_v ;
  typedef zimt::simdized_type < px_t , LANES > px_v ;

  typedef zimt::bspline < px_t , 2 > spl_t ;

  std::shared_ptr < spl_t > p_bsp ;

  // we hold the functor producing pixels as a grok_type, which
  // gives us more flexibility because we can override it later on,
  // e.g. to produce masks instead of image data.

  zimt::grok_type < crd2_t , px_t , LANES > ev ;

  // these are the variables which are needed from cubemap_t to
  // evaluate the b-spline holding the IR image

  const float refc_md , model_to_px ;
  const int section_px ;

  cubemap_view_t ( float _refc_md ,
                   float _model_to_px ,
                   std::shared_ptr < spl_t > _p_bsp ,
                   int masked )
  : refc_md ( _refc_md ) ,
    model_to_px ( _model_to_px ) ,
    p_bsp ( _p_bsp ) ,
    section_px ( _p_bsp->core.shape [ 0 ] )
  {
    if ( masked != -1 )
    {
      assert ( nchannels == 2 || nchannels == 4 ) ;
      ev = alpha_masking_t < nchannels , LANES > ( masked , p_bsp ) ;
    }
    else
    {
      // std::cout << "making an evaluator from " << p_bsp << std::endl ;
      ev = make_safe_evaluator < spl_t , float , LANES > ( * p_bsp ) ;
    }
  }

  // this is the same code as in metrics.h - it's all we need to
  // obtain the pick-up coordinate in spline units

  void get_pickup_coordinate_px ( const i_v & face ,
                                  const crd2_v & in_face ,
                                  crd2_v & pickup ) const
  {
    pickup = in_face + refc_md ;
    pickup *= model_to_px ;
    pickup[1] += ( face * section_px ) ;
    pickup -= .5f ;
  }

  // and this is the same code as in cubemap.h

  void cubemap_to_pixel ( const i_v & face ,
                          const crd2_v & in_face ,
                          px_v & px )
  {
    crd2_v pickup ;
    get_pickup_coordinate_px ( face , in_face , pickup ) ;
    ev.eval ( pickup , px ) ;
  }

  void eval ( const ray_v & ray , px_v & px )
  {
    i_v face ;
    crd2_v in_face ;

    ray_to_cubeface < float , LANES > ( ray , face , in_face ) ;

    if constexpr ( cbm_prj == BIATAN6 )
      in_face = float ( 4.0 / M_PI ) * atan ( in_face ) ;

    // std::cout << "face: " << face << " in_face " << in_face << std::endl ;
    cubemap_to_pixel ( face , in_face , px ) ;
    // std::cout << "px: " << px << std::endl ;
  }

} ;

// the 'environment' template codes objects which can serve as 'act'
// functor in zimt::process. It's coded as a zimt::unary_functor
// taking 3D 'ray' coordinates and producing pixels with C channels.
// 'environment' objects also provide the substrate for 'twining'.
// All types of environment objects are based on cardinal b-splines
// over the image data. The default is to use degree-1 b-splines,
// which is also known as 'bilinear interpolation'. Other spline
// degrees can be specified by passing 'spline_degree'. If necessary,
// the image data are prefiltered, so as to provide an interpolating
// spline.

template < typename T , typename U , std::size_t C , std::size_t L >
struct _environment
: public zimt::unary_functor < zimt::xel_t < T , 3 > ,
                               zimt::xel_t < U , C > ,
                               L >
{
  typedef zimt::unary_functor < zimt::xel_t < T , 3 > ,
                               zimt::xel_t < U , C > ,
                               L > base_t ;

  using typename base_t::in_v ;
  using typename base_t::in_ele_v ;
  using typename base_t::out_v ;
  typedef typename in_ele_v::mask_type mask_t ;

  typedef zimt::xel_t < T , 2 > crd2_t ;
  typedef simdized_type < crd2_t , L > crd_v ;
  typedef zimt::xel_t < T , 3 > ray_t ;
  typedef simdized_type < ray_t , L > ray_v ;
  typedef zimt::xel_t < U , C > px_t ;
  typedef simdized_type < px_t , L > px_v ;
  typedef zimt::xel_t < std::size_t , 2 > shape_type ;
  typedef zimt::xel_t < float , C > in_px_t ;

  // this member variable will hold the type-erased functor encoding
  // the various code paths we handle with this object. The 'eval'
  // member function simply calls this functor.

  zimt::grok_type < ray_t , px_t , L > env ;

  // we'll hold on to image data via a std::shared_ptr to a zimt::bspline.

  typedef zimt::bspline < in_px_t , 2 > spl_t ;
  std::shared_ptr < spl_t > p_bspl = nullptr ;

  typedef std::function < mask_t ( const in_v & ) > mask_f ;
  mask_f get_mask ;

  // create an 'environment' object from a facet specification
  // Some facets may cover the full 360X180 degrees, but most will
  // only cover a part. The std::function get_mask, above, will
  // produce a mask which indicates which of the gleaned pixels
  // originate from the facet image - for facets with 'full cover',
  // it will return an all-true mask unconditionally, and we hope
  // that the optimizer can use this information as well to
  // streamline processing for such facets.

  _environment ( const facet_spec & fct )
  {
    // sneaky: this might as well be a member variable, but having
    // it static to a member function saves the need for external
    // instantiation

    #define PSPL(NCH) std::shared_ptr \
      < bspline < xel_t < float , NCH > , 2 > >

    static std::map < std::string , PSPL(C) > spl_map ;

    #undef PSPL

    // if the facet is in a cubemap format, we special-case here
    // and directly form a cubemap_t object with the facet_spec
    // as argument. get_mask in this case unconditionally returns
    // an all-true mask.

    if ( fct.projection == CUBEMAP || fct.projection == BIATAN6 )
    {
      metrics_t cbm_metrics ( fct.width , fct.hfov ,
                              args.support_min ,
                              args.tile_size ) ;

      // collision of terminology: get_mask gets a SIMD mask which
      // is true for rays which hit the facet. obviously, this must
      // be an all-tru mask for cubemaps, because they cover the
      // entire 360X180. The masks we refer to next are masks
      // which are 1 - or white - where a facet produces visible
      // contribution to the final image. To form such masks,
      // what we're doing conceptually is to replace one facet
      // image's pixels with entirely white pixels and all other
      // facet images' pixels with black pixels. If we form an
      // output image from such modified data, we receive the
      // desired mask.

      get_mask = []( const ray_v & ray ) { return mask_t ( true ) ; } ;

      // if we're doing a masking job for a facet without alpha
      // channel (one or three channels, greyscale and RGB), we
      // create a special kind of environment object where we
      // don't even have to open the image: all we need is it's
      // metrics, which can be gleaned from the facet_spec in fct.

      if ( ( fct.masked != -1 ) && ( C == 1 || C == 3 ) )
      {
        env = masking_t < 3 , C , L > ( fct.masked ) ;
      }
      else 
      {
        // for 'normal' operation, and also for masking jobs with
        // facets with alpha channel, we actually need image data.
        // image data are valuable assets, so we keep track of
        // images we load in 'spl_map' and reuse the data in RAM
        // if we already have them.
        // Note: for cubemaps, masking and cropping are not available,
        // we assume that cubemaps always have complete, opaque data.
        // So we reuse cubemap image data if the name matches. Normal
        // facets arealways reloaded if they have cropping or masking. 

        auto it = spl_map.find ( fct.asset_key ) ;

        if ( ( it == spl_map.end() ) )
        {
          if ( args.verbose )
            std::cout << "cubemap " << fct.asset_key
                      << " is now loaded from disk" << std::endl ;

          // The b-spline used inside a cubemap_t object is not simply
          // a spline over the image data, but the result of quite
          // elaborate processing, and possibly a composite of six
          // discrete cube face images. So we set up a cubemap_t
          // object and then extract it's b-spline. Later on, we'll
          // us a cubemap_view_t object to obtain pixel data from
          // that spline.

          if ( fct.projection == CUBEMAP )
          {
            cubemap_t < C , CUBEMAP > cbm ( fct.width , fct.hfov ,
                                            args.support_min ,
                                            args.tile_size ) ;
            cbm.load ( fct.filename ) ;
            p_bspl = cbm.p_bsp ;
            spl_map [ fct.asset_key ] = p_bspl ;
          }
          else
          {
            cubemap_t < C , BIATAN6 > cbm ( fct.width , fct.hfov ,
                                            args.support_min ,
                                            args.tile_size ) ;
            cbm.load ( fct.filename ) ;
            p_bspl = cbm.p_bsp ;
            spl_map [ fct.asset_key ] = p_bspl ;
          }
        }
        else
        {
          // we're in luck - this cubemap was already loaded
          // into a cubemap_t object previously and the resulting
          // bspline object was extracted.

          if ( args.verbose )
            std::cout << "cubemap " << fct.filename
                      << " is already present in RAM" << std::endl ;

          p_bspl = it->second ;
        }

        // now that we have the b-spline, we can set up a cubemap_view_t.
        // This only needs a few parameters characterizing the meaning
        // of the spline (refc_md and model_to_px), the spline itself
        // and a variable indicating whether we want pixel or masking
        // data.

        if ( fct.projection == CUBEMAP )
        {
          env = cubemap_view_t < C , CUBEMAP >
                 ( cbm_metrics.refc_md ,
                   cbm_metrics.model_to_px ,
                   p_bspl , fct.masked ) ;
        }
        else
        {
          env = cubemap_view_t < C , BIATAN6 >
                 ( cbm_metrics.refc_md ,
                   cbm_metrics.model_to_px ,
                   p_bspl , fct.masked ) ;
        }
      }
      // we're done with the case that the facet is a cubemap, so we
      // return early - the remainder of this function deals with other
      // types of facets.

      return ;
    }

    // the facet isn't a cubemap. First we check if the facet provides
    // full 360X180 degree coverage, in which case we can use slightly
    // more performant code (we need no SIMD masks to indicate where
    // the facet provides data and where it doesn't).

    shape_type shape { fct.width , fct.height } ;
    bool full_environment = false ;

    // most facets are rendered with REFLECT boundary conditions,
    // but if the facet covers the full 360 degrees in the horiontal,
    // we'll use PERIODIC boundary conditions instead.

    zimt::bc_code bc0 = zimt::REFLECT ;
    if (    fct.projection == SPHERICAL
         || fct.projection == CYLINDRICAL )
    {
      if ( fabs ( fct.hfov - 2.0 * M_PI ) < .000001 )
        bc0 = PERIODIC ;
    }

    // fct.masked == -1 means 'normal operation': we set up a
    // b-spline over the image data. The image data are also
    // needed for a masking job for facets with alpha channel.
    // masking jobs for images without an alpha channel don't
    // need image data, so we don't set up a b-spline at all
    // and leave p_bspl at it's default value of nullptr.
    // for 'solo' jobs, it's also futile to load image data
    // for all but the 'solo' facet, so if we have a 'solo'
    // argument and the currently handled facet is not the solo
    // facet, we also skip over the image-loading code and leave
    // p_bspl nullptr.

    bool load_facet = ( fct.masked == -1 || ( C == 2 || C == 4 ) ) ;

    if ( args.solo >= 0 && fct.facet_no != args.solo )
      load_facet = false ;

    if ( load_facet )
    {
      auto it = spl_map.find ( fct.asset_key ) ;

      if ( it != spl_map.end() ) // && alpha_modified == false )
      {
        // we're in luck - this image file was already loaded
        // into a b-spline previously.

        if ( args.verbose )
          std::cout << "asset " << fct.asset_key
                    << " is already present in RAM" << std::endl ;

        p_bspl = it->second ;
      }
      else
      {
        // no luck - we need to load the image data from disk

        if ( args.verbose )
          std::cout << "asset " << fct.asset_key
                    << " is now loaded from disk" << std::endl ;

        // currently building with raw::user_flip set to zero, to load
        // raw images in memory order without EXIF rotation. This only
        // affects raw images.

        ImageSpec config;
        config [ "raw:user_flip" ] = 0 ;
        config [ "raw:ColorSpace" ] = "sRGB-linear" ;

        auto inp = ImageInput::open ( fct.filename , &config ) ;

        const ImageSpec &spec = inp->spec() ;

        p_bspl.reset ( new spl_t ( shape , args.spline_degree ,
                                    { bc0 , zimt::REFLECT } ) ) ;

        bool success = inp->read_image (
          0 , 0 , 0 , C ,
          TypeDesc::FLOAT ,
          p_bspl->core.data() ,
          sizeof ( in_px_t ) ,
          p_bspl->core.strides[1] * sizeof ( in_px_t ) ) ;
        assert ( success ) ;

        inp->close() ;

        int native_nchannels = spec.nchannels ;

        if ( native_nchannels != fct.nchannels )
        {
          // this occurs only when cropping or masking affect the facet,
          // in which case the facet's 'nchannels' value is raised from
          // one to two or from three to four, to make room for an alpha
          // channel if that isn't already present.

          assert ( fct.have_crop || fct.have_pto_mask ) ;

          // we initialize the alpha channel to 1.0f

          auto p_alpha_data = (float*) ( p_bspl->container.data() ) ;
          p_alpha_data += ( C - 1 ) ;

          view_t < 2 , float > alpha_view
                        ( p_alpha_data ,
                          p_bspl->container.strides * C ,
                          p_bspl->container.shape ) ;

          alpha_view.set_data ( 1.0f ) ;
        }

        // now we process masking and cropping information

        if ( fct.have_crop || fct.have_pto_mask )
        {
          assert ( C == 2 || C == 4 ) ;

          array_t < 2 , float > alpha ( p_bspl->core.shape ) ;
          alpha.set_data ( 1.0f ) ;

          int w = alpha.shape[0] ;
          int h = alpha.shape[1] ;

          if ( fct.have_pto_mask )
          {
            auto clear = [&] ( int x , int y )
            {
              alpha [ { x , y } ] = 0.0f ;
            } ;

            for ( const auto & mask : fct.pto_mask_v )
            {
              if ( mask.variant == 0 )
              {
                fill_polygon ( mask.vx , mask.vy ,
                               0 , 0 , w , h , clear ) ;
              }
            }
          }
  
          if ( fct.have_crop )
          {
            float a = fabs ( fct.crop_x1 - fct.crop_x0 ) / 2.0 ;
            float b = fabs ( fct.crop_y1 - fct.crop_y0 ) / 2.0 ;

            int w = p_bspl->core.shape [ 0 ] ;
            int h = p_bspl->core.shape [ 1 ] ;

            if ( fct.projection == FISHEYE )
            {
              std::cout << "elliptic crop" << std::endl ;
              float mx = ( fct.crop_x0 + fct.crop_x1 ) / 2.0 ;
              float my = ( fct.crop_y0 + fct.crop_y1 ) / 2.0 ;
              for ( int y = 0 ; y < h ; y++ )
              {
                auto dy = fabs ( y - my ) ;
                // if the y coordinate is outside the ellipse, mask out
                // the entire line
                if ( dy > b )
                {
                  for ( int x = 0 ; x < w ; x++ )  
                  {
                    alpha [ { x , y } ] = 0 ;
                  }
                  continue ;
                }
                // else, find the half width of the ellipse at given y and
                // mask out points with x outside
                // for the ellipse, we have x*x / a*a + y*y / b*b = 1, hence

                float xmargin = sqrt ( ( a * a ) * ( 1.0 - ( dy * dy ) / ( b * b ) ) ) ;
                for ( int x = 0 ; x < w ; x++ )
          
                {
                  auto dx = fabs ( x - mx ) ;
                  if ( dx > xmargin )
                    alpha [ { x , y } ] = 0 ;
                }
              }
            }
            else
            {
              std::cout << "rectangular crop" << std::endl ;
              for ( int y = 0 ; y < h ; y++ )
              {
                for ( int x = 0 ; x < w ; x++ )
                {
                  if (   ( x < fct.crop_x0 )
                      || ( x >= fct.crop_x1 )
                      || ( y < fct.crop_y0 )
                      || ( y >= fct.crop_y1 ) )
                    alpha [ { x , y } ] = 0 ;
                }
              }
            }
          }

          // finally, apply a low pass filter to the masking alpha channel
          // this mitigates the issue of the hard mask boundaries, but it
          // will 'pull in' a bit of masked-out content, so if the mask is
          // cut very close to unwanted features, they may bleed in. The
          // large-ish binomial is tentative, but seems to work quite well.

          convolve
          ( alpha ,
            alpha ,
            { REFLECT , REFLECT } ,
            { 1.0 / 16.0 ,
              4.0 / 16.0 ,
              6.0 / 16.0 ,
              4.0 / 16.0 ,
              1.0 / 16.0
            } ,
            2 ) ;

          // we might use a loop to apply the alpha mask, like this:

          // for ( std::size_t y = 0 ; y < h ; y++ )
          // {
          //   for ( std::size_t x = 0 ; x < w ; x++ )
          //   {
          //     p_bspl->core [ { x , y } ] *= alpha [ { x , y } ] ;
          //   }
          // }

          // but we can do better and do the job in multithreaded
          // SIMD code using zimt:

          // we set up two loaders, one for the image data in the
          // b-spline's core, one for the alpha mask we've just made.
          
          loader < float , C , 2 , LANES > ldpx ( p_bspl->core ) ;
          loader < float , 1 , 2 , LANES > lda ( alpha ) ;
          
          // we use a synopsis-forming lambda which produces the
          // product of it's two inputs
          
          auto syn = [] ( const px_v & v1 , const f_v & v2 ,
                          px_v & v3 , std::size_t cap = LANES )
          {
            v3 = v1 * v2 ;
          } ;
          
          // and set up a zip_t wiring the two loaders and the
          // synopsis object to produce the masked image data
          
          zip_t < float , C , 2 , LANES , decltype ( syn ) ,
                  float , 1 , float , C > zip ( ldpx , lda , syn ) ;
          
          // the act functor is not used
          
          pass_through < float , C , LANES > pass ;
          
          // the masked image data go back into the b-spline's core
          
          storer < float , C , 2 , LANES > store ( p_bspl->core ) ;
          
          // showtime
          
          process ( alpha.shape , zip , pass , store ) ;
        }

        // we save the shared_ptr to the newly made b-spline in
        // spl_map, to make it available for reuse if needed.

        spl_map [ fct.asset_key ] = p_bspl ;

        // the spline needs to be prefiltered appropriately.

        if (    fct.projection == SPHERICAL
            && fabs ( fct.hfov - 2.0 * M_PI ) < .000001
            && fct.width == 2 * fct.height )
        {
          // assume the mounted image is a full spherical

          p_bspl->spline_degree = args.prefilter_degree ;
          spherical_prefilter ( *p_bspl , p_bspl->core , zimt::default_njobs ) ;
          p_bspl->spline_degree = args.spline_degree ;
          full_environment = true ;
        }
        else
        {
          // assume it's a partial spherical, use ordinary
          // prefilter. Note that we may prefilter with a
          // different degree (prefilter_degree) if that was
          // specified in args.

          p_bspl->spline_degree = args.prefilter_degree ;
          p_bspl->prefilter() ;
          p_bspl->spline_degree = args.spline_degree ;
        }
      }
    }

    // if we're making a mask for an image without alpha channel,
    // p_bspl is nullptr (we don't need image data), and fct.masked
    // is 1 for the facet for which we make the mask and zero for
    // all other facets. source_t's c'tor will result in a suitable
    // functor: one producing only zero, or one, respectively.
    // For 'normal' operation and masking jobs for images with
    // alpha channel, p_bspl points to a b-spline.

    source_t < C , 2 , L > src ( p_bspl , fct.masked ) ;

    // for now, we mount images to the center; other types of cropping
    // might be added by providing suitable parameterization.

    auto extent = get_extent ( fct.projection , fct.width ,
                               fct.height , fct.hfov  ) ;

    // Some facets require additional processing of the planar (image)
    // coordinates which the source_t object receives. This is needed
    // if the PTO file specifies shear and/or lens correction parameters.
    // This can't currently be affected by command line parameters.
    // We encode such processing in a std::function modifying the
    // SIMD vector of 2D coordinates:

    std::function < void ( crd_v & ) > pf = [] ( crd_v & ) { } ;

    // with lens control active, we can use the handy 'pto_planar'
    // class which combines lens correction polynomial, lens shift
    // and shear:

    if ( fct.lens_correction_active )
    {
      pf = pto_planar < float , LANES >
             ( fct.a , fct.b , fct.c , fct.s , fct.r_max ,
               fct.h , fct.v , fct.shear_g , fct.shear_t ) ;
    }
    else if ( fct.shear_g != 0.0 || fct.shear_t != 0.0 )
    {
      pf = [=] ( crd_v & crd )
      {
        crd = {  crd[0] + ( crd[1] * fct.shear_g ) ,
                 crd[1] + ( crd[0] * fct.shear_t ) } ;
      } ;
    }

    // we fix the projection as a template argument to class mount_t,
    // passing the planar function as well. Then we also set up the
    // get_mask function - this is specific to the mount_t, and since
    // we 'grok' the mount_t by assigning to 'env' we can't perform
    // the operation later on.

    switch ( fct.projection )
    {
      case RECTILINEAR:
      {
        mount_t < C , RECTILINEAR , L > mnt ( extent , src , pf ) ;
        env = mnt ;
        get_mask = [=] ( const in_v & crd3 )
          { return mnt.get_mask ( crd3 ) ; } ;
        break ;
      }
      case SPHERICAL:
      {
        mount_t < C , SPHERICAL , L > mnt ( extent , src , pf ) ;
        env = mnt ;
        if ( full_environment )
          get_mask = []( const ray_v & ray )
            { return mask_t ( true ) ; } ;
        else
          get_mask = [=] ( const in_v & crd3 )
            { return mnt.get_mask ( crd3 ) ; } ;
        break ;
      }
      case CYLINDRICAL:
      {
        mount_t < C , CYLINDRICAL , L > mnt ( extent , src , pf ) ;
        env = mnt ;
        get_mask = [=] ( const in_v & crd3 )
          { return mnt.get_mask ( crd3 ) ; } ;
        break ;
      }
      case STEREOGRAPHIC:
      {
        mount_t < C , STEREOGRAPHIC , L > mnt ( extent , src , pf ) ;
        env = mnt ;
        get_mask = [=] ( const in_v & crd3 )
          { return mnt.get_mask ( crd3 ) ; } ;
        break ;
      }
      case FISHEYE:
      {
        mount_t < C , FISHEYE , L > mnt ( extent , src , pf ) ;
        env = mnt ;
        if ( fct.hfov >= M_PI * 2.0 )
          get_mask = []( const ray_v & ray )
            { return mask_t ( true ) ; } ;
        else
          get_mask = [=] ( const in_v & crd3 )
            { return mnt.get_mask ( crd3 ) ; } ;
        break ;
      }
      default:
      {
        std::cerr << "unknown projection: "
                  << fct.projection_str << std::endl ;
        assert ( false ) ;
        break ;
      }
    }
  }

  // eval simply delegates to 'env'

  void eval ( const in_v & in , out_v & out )
  {
    env.eval ( in , out ) ;
  }

} ;

// class environment is a wrapper around class _environment,
// which has te workhorse code. the wrapper takes care of
// converting pixels with the number of channels specific
// to the facet to pixels with C channels, the number of
// channels requested by the given template argument.
// the conversion is handled by a repix_t object.

template < typename T , typename U , std::size_t C ,
           std::size_t L >
struct environment
: public zimt::unary_functor < zimt::xel_t < T , 3 > ,
                               zimt::xel_t < U , C > ,
                               L >
{
  typedef zimt::unary_functor < zimt::xel_t < T , 3 > ,
                                zimt::xel_t < U , C > ,
                                L > base_t ;

  using typename base_t::in_v ;
  using typename base_t::in_ele_v ;
  using typename base_t::out_v ;
  typedef typename in_ele_v::mask_type mask_t ;

  typedef zimt::xel_t < T , 2 > crd2_t ;
  typedef simdized_type < crd2_t , L > crd_v ;
  typedef zimt::xel_t < T , 3 > ray_t ;
  typedef simdized_type < ray_t , L > ray_v ;
  typedef zimt::xel_t < U , C > px_t ;
  typedef zimt::xel_t < std::size_t , 2 > shape_type ;
  typedef zimt::xel_t < float , C > in_px_t ;

  typedef std::function < mask_t ( const in_v & ) > mask_f ;
  mask_f get_mask ;
  const float recip_step ;
  const float brighten ;

  grok_type < zimt::xel_t < T , 3 > ,
              zimt::xel_t < U , C > ,
              L > env ;

  void eval ( const in_v & in , out_v & out )
  {
    env.eval ( in , out ) ;

    if ( brighten != 1.0f )
    {
      if constexpr ( C == 1 || C == 3 )
      {
        out *= brighten ;
      }
      else if constexpr ( C == 2 )
      {
        out[0] *= brighten ;
      }
      else if constexpr ( C == 4 )
      {
        out[0] *= brighten ;
        out[1] *= brighten ;
        out[2] *= brighten ;
      }
    }
  }

  environment ( const facet_spec & fct )
  : recip_step ( 1.0 / fct.step ) ,
    brighten ( fct.brighten )
  {
    if ( fct.nchannels == C )
    {
      // the facet has the required channel count

      _environment < T , U , C , L > e ( fct ) ;
      env = e.env ;
      get_mask = e.get_mask ;
      return ;
    }

    // the required channel count is different, use a repix_t

    if ( fct.masked == -1 )
    {
      switch ( fct.nchannels )
      {
        case 1 :
        {
          repix_t < U , 1 , C , L > repix ;
          _environment < T , U , 1 , L > e ( fct ) ;
          env = e.env + repix ;
          get_mask = e.get_mask ;
          break ;
        }
        case 2 :
        {
          repix_t < U , 2 , C , L > repix ;
          _environment < T , U , 2 , L > e ( fct ) ;
          env = e.env + repix ;
          get_mask = e.get_mask ;
          break ;
        }
        case 3 :
        {
          repix_t < U , 3 , C , L > repix ;
          _environment < T , U , 3 , L > e ( fct ) ;
          env = e.env + repix ;
          get_mask = e.get_mask ;
          break ;
        }
        case 4 :
        {
          repix_t < U , 4 , C , L > repix ;
          _environment < T , U , 4 , L > e ( fct ) ;
          env = e.env + repix ;
          get_mask = e.get_mask ;
          break ;
        }
        default:
        {
          std::cerr << "fct.nchannels invalid: " << fct.nchannels
                    << std::endl ;
          exit ( -1 ) ;
        }
      } ;
    }
    else
    {
      // for masking jobs, we use mono_t instead of repix_t. This will
      // only work if C == 1 or C == 2.

      switch ( fct.nchannels )
      {
        case 1 :
        {
          mono_t < U , 1 , C , L > mono ;
          _environment < T , U , 1 , L > e ( fct ) ;
          env = e.env + mono ;
          get_mask = e.get_mask ;
          break ;
        }
        case 2 :
        {
          mono_t < U , 2 , C , L > mono ;
          _environment < T , U , 2 , L > e ( fct ) ;
          env = e.env + mono ;
          get_mask = e.get_mask ;
          break ;
        }
        case 3 :
        {
          mono_t < U , 3 , C , L > mono ;
          _environment < T , U , 3 , L > e ( fct ) ;
          env = e.env + mono ;
          get_mask = e.get_mask ;
          break ;
        }
        case 4 :
        {
          mono_t < U , 4 , C , L > mono ;
          _environment < T , U , 4 , L > e ( fct ) ;
          env = e.env + mono ;
          get_mask = e.get_mask ;
          break ;
        }
        default:
        {
          std::cerr << "fct.nchannels invalid: " << fct.nchannels
                    << std::endl ;
          exit ( -1 ) ;
        }
      } ;
    }
  }
} ;

// environment9 objects mediate lookup with derivatives. Their eval
// member function takes 'ninepacks' and provides pixels with
// nchannels channels. This is used with 'twining'.

template < std::size_t nchannels , std::size_t L >
struct environment9
: public zimt::unary_functor < zimt::xel_t < float , 9 > ,
                               zimt::xel_t < float , nchannels > ,
                               L
                             >
{
  typedef zimt::unary_functor < zimt::xel_t < float , 9 > ,
                                zimt::xel_t < float , nchannels > ,
                                L
                              > base_t ;

  typedef zimt::xel_t < float , 9 > crd9_t ;
  typedef zimt::xel_t < float , nchannels > px_t ;

  // the 'act' functor's 'inner' type is variable, so we use a
  // uniform 'outer' type by 'grokking' it (type erasure)

  zimt::grok_type < crd9_t , px_t , L > act ;

  typedef environment < float , float , nchannels , L > env_t ;

  environment9 ( env_t * p_env )
  {
    // wrap the 'environment' object in a twine_t object and assign
    // to 'act' - this 'groks' the twine_t object to act's type.
    // Note how this object is re-created for every run: the twining
    // parameters may change due to changing hfov from one invocation
    // of 'work' to the next.

    if ( args.twine_precise )
    {
      // TODO: I think this is futile.
      act = twine_t < nchannels , L , true >
        ( *p_env , args.twine_spread ) ;
    }
    else
    {
      act = twine_t < nchannels , L , false >
        ( *p_env , args.twine_spread ) ;
    }
  }

  // 'eval' simply delegates to the grokked specific environment code

  template < typename I , typename O >
  void eval ( const I & in , O & out )
  {
    act.eval ( in , out ) ;
  }
} ;

END_ZIMT_SIMD_NAMESPACE
HWY_AFTER_NAMESPACE() ;

#endif // sentinel
