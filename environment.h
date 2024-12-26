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


#include <filesystem>
#include <OpenImageIO/imageio.h>
#include <OpenImageIO/texture.h>

// some versions of OIIO yield a plain pointer, some a shared_ptr.
typedef decltype ( OIIO::TextureSystem::create() ) ts_ptr_t ;

#include "zimt/prefilter.h"
#include "zimt/bspline.h"
#include "zimt/eval.h"

#include "metrics.h"
#include "twining.h"
#include "geometry.h"

#if defined(ENVUTIL_ENVIRONMENT_H) == defined(HWY_TARGET_TOGGLE)
  #ifdef ENVUTIL_ENVIRONMENT_H
    #undef ENVUTIL_ENVIRONMENT_H
  #else
    #define ENVUTIL_ENVIRONMENT_H
  #endif

HWY_BEFORE_NAMESPACE() ;
BEGIN_ZIMT_SIMD_NAMESPACE(project)

using OIIO::ImageInput ;
using OIIO::TypeDesc ;
using OIIO::ImageSpec ;

bool verbose = true ;

// This functor 'looks at' a full spherical image. It receives
// 2D lat/lon coordinates and yields pixel values. Pick-up is done
// with bilinear interpolation. For this demo program, we keep it
// simple, so this is a stripped-down version of what's used in
// envutil - we omit using OIIO's 'environment' function and
// also the use of 'twining'.

template < std::size_t nchannels >
struct eval_latlon
: public zimt::unary_functor
           < v2_t , zimt::xel_t < float , nchannels > , LANES >
{
  typedef zimt::xel_t < float , nchannels > px_t ;
  typedef zimt::simdized_type < px_t , LANES > px_v ;

  // latlon contains a view to an image in 2:1 aspect ratio in
  // spherical projection.

  const zimt::view_t < 2 , px_t > & latlon ;

  // scaling factor to move from model space coordinates to
  // image coordinates

  const float scale ;

  eval_latlon ( const zimt::view_t < 2 , px_t > & _latlon )
  : latlon ( _latlon ) ,
    scale ( _latlon.shape[1] / M_PI )
  { }

  // direct lookup from the lat/lon image with bilinear interpolation.
  // we expect incoming coordinates which are outside the range of
  // the lat/lon image to be 'close' to the range and avoid the
  // calculations which would be required to accept completely
  // arbitrary coordinates (this would require a divison)
  // The same applies to neighbouring pixels needed for bilinear
  // interpolation. The special geometry of a full spherical is
  // honoured: periodicity in the horizontal and periodicity with
  // a longitute jump of pi over the poles.

  template < typename I , typename O >
  void eval ( const I & _in , O & out )
  {
    // we have incoming model space coordinates, in[0] is the
    // longitude in the range of [-pi,pi] and in[1] is the
    // latitude in the range of [-pi/2,pi/2]. First we move
    // to image coordinates.

    const zimt::xel_t < float , 2 > shift { M_PI , M_PI_2 } ;

    auto in = ( _in + shift ) * scale ;

    // if the coordinate is now precisely zero, this corresponds
    // to a pick-up point at the top or left margin of the UL
    // pixel, which is 0.5 away from it's center

    in -= .5f ;

    // now we move to discrete values for the index calculations
    // the nearest discrete coordinate with components smaller
    // than or equal to 'in' is found by applying 'floor'. We
    // store the result in 'uli' for 'upper left index'.

    v2i_v uli { floor ( in[0] ) , floor ( in[1] ) } ;

    // the distance of this coordinate, both in horizontal and
    // vertical direction, is stored in 'wr' ('weight right'),
    // and it will be used to weight the right-hand part of the
    // constituents (ur and lr)

    auto wr = in - uli ;

    // 'wl' ('weight left') is used for the opposite constituents,
    // namely ul and ll.

    auto wl = 1.0f - wr ;

    // from the discrete upper left value, we derive it's three
    // neighbours downwards and to the right.

    v2i_v uri ( uli ) ;
    uri[0] += 1 ;
    
    v2i_v lli ( uli ) ;
    lli[1] += 1 ;
    
    v2i_v lri ( lli ) ;
    lri[0] += 1 ;
    
    // shorthand for width and height

    int w = latlon.shape[0] ;
    int h = latlon.shape[1] ;

    // first we look at the vertical axis and map excessive values
    // back into the range. Note how this isn't as straightforward
    // as handling the horizontal: The continuation is on the
    // opposite hemisphere (add w/2 to the horizontal component)
    // then mirror on the pole. Note also that this won't work with
    // completely arbitrary coordinates, which would be more expensive
    // mathematically and isn't needed in this context.

    uli[0] ( uli[1] < 0 ) += w / 2 ;
    uli[1] ( uli[1] < 0 ) = 1 - uli[1] ; // e.g. -1 -> 0, opposite

    uri[0] ( uri[1] < 0 ) += w / 2 ;
    uri[1] ( uri[1] < 0 ) = 1 - uri[1] ; // e.g. -1 -> 0, opposite

    lli[0] ( lli[1] >= h ) += w / 2 ;
    lli[1] ( lli[1] >= h ) = ( h - 1 ) - ( lli[1] - h ) ;

    lri[0] ( lri[1] >= h ) += w / 2 ;
    lri[1] ( lri[1] >= h ) = ( h - 1 ) - ( lri[1] - h ) ;

    // now we look at the horizontal axis - the longitude axis -
    // and map any coordinates which are outside the range back in,
    // exploiting the periodicity. Note that the code for the vertical
    // may have put values way outside the range (by adding w/2), but
    // not so far as that they wouldn't be mapped back into the range
    // now. We don't expect uri and lri to ever have negative values.

    uli[0] ( uli[0] < 0 ) += w ;
    lli[0] ( lli[0] < 0 ) += w ;

    uli[0] ( uli[0] >= w ) -= w ;
    uri[0] ( uri[0] >= w ) -= w ;
    lli[0] ( lli[0] >= w ) -= w ;
    lri[0] ( lri[0] >= w ) -= w ;

    // base pointer for the gather operation

    const auto * const p = (const float* const) ( latlon.data() ) ;

    // problems on my mac: vectorized long doesn't work
    zimt::xel_t < int , 2 > llstrides ( latlon.strides ) ;
    int nch = nchannels ;

    // obtain the four constituents by first truncating their
    // coordinate to int and then gathering from p.

    index_v idsdxl { uli[0] , uli[1] } ;
    auto ofs = ( idsdxl * llstrides ) . sum() * nch ;
    px_v pxul ;
    pxul.gather ( p , ofs ) ;

    index_v idsdxr { uri[0] , uri[1] } ;
    ofs = ( idsdxr * llstrides ) . sum() * nch ;
    px_v pxur ;
    pxur.gather ( p , ofs ) ;

    index_v idxll { lli[0] , lli[1] } ;
    ofs = ( idxll * llstrides ) . sum() * nch ;
    px_v pxll ;
    pxll.gather ( p , ofs ) ;

    index_v idxlr { lri[0] , lri[1] } ;
    ofs = ( idxlr * llstrides ) . sum() * nch ;
    px_v pxlr ;
    pxlr.gather ( p , ofs ) ;

    // apply the bilinear formula with the weights gleaned above

    out  = wl[1] * ( wl[0] * pxul + wr[0] * pxur ) ;
    out += wr[1] * ( wl[0] * pxll + wr[0] * pxlr ) ;
  }
} ;

// some helper code to pass OIIO texture option parameters as
// strings from the command line

std::map < std::string , OIIO::TextureOpt::Wrap >
wrap_map
{
  { "WrapDefault" , OIIO::TextureOpt::WrapDefault } ,
  { "WrapBlack" , OIIO::TextureOpt::WrapBlack } ,
  { "WrapClamp" , OIIO::TextureOpt::WrapClamp } ,
  { "WrapPeriodic" , OIIO::TextureOpt::WrapPeriodic } ,
  { "WrapMirror" , OIIO::TextureOpt::WrapMirror } ,
  { "WrapPeriodicPow2" , OIIO::TextureOpt::WrapPeriodicPow2 } ,
  { "WrapPeriodicSharedBorder" ,
      OIIO::TextureOpt::WrapPeriodicSharedBorder }
} ;

std::map < std::string , OIIO::TextureOpt::MipMode >
mipmode_map
{
  { "MipModeDefault" , OIIO::TextureOpt::MipModeDefault } ,
  { "MipModeNoMIP" , OIIO::TextureOpt::MipModeNoMIP } ,
  { "MipModeOneLevel" , OIIO::TextureOpt::MipModeOneLevel } ,
  { "MipModeTrilinear" , OIIO::TextureOpt::MipModeTrilinear } ,
  { "MipModeAniso" , OIIO::TextureOpt::MipModeAniso } ,
  // { "MipModeStochasticTrilinear" , OIIO::TextureOpt::MipModeStochasticTrilinear } ,
  // { "MipModeStochasticAniso" , OIIO::TextureOpt::MipModeStochasticAniso }
} ;

std::map < std::string , OIIO::TextureOpt::InterpMode >
interpmode_map
{
  { "InterpClosest" , OIIO::TextureOpt::InterpClosest } ,
  { "InterpBilinear" , OIIO::TextureOpt::InterpBilinear } ,
  { "InterpBicubic" , OIIO::TextureOpt::InterpBicubic } ,
  { "InterpSmartBicubic" , OIIO::TextureOpt::InterpSmartBicubic }
} ;

template < std::size_t nchannels >
struct eval_env
: public zimt::unary_functor
   < zimt::xel_t < float , 9 > ,
     zimt::xel_t < float , nchannels > ,
     LANES >
{
  ts_ptr_t ts ;
  OIIO::TextureOptBatch batch_options ;
  OIIO::TextureSystem::TextureHandle * th ;

  // pull in the c'tor arguments

  eval_env()
  {
    ts = OIIO::TextureSystem::create() ;
    ts->attribute ( "options" , args.tsoptions ) ;

    for ( int i = 0 ; i < 16 ; i++ )
      batch_options.swidth[i] = batch_options.twidth[i]
        = args.stwidth ;
    
    for ( int i = 0 ; i < 16 ; i++ )
      batch_options.sblur[i] = batch_options.tblur[i]
        = args.stblur ;

    batch_options.conservative_filter = args.conservative_filter ;

    auto wrap_it = wrap_map.find ( args.swrap ) ;
    assert ( wrap_it != wrap_map.end() ) ;
    batch_options.swrap = OIIO::Tex::Wrap ( wrap_it->second ) ;
  
    wrap_it = wrap_map.find ( args.twrap ) ;
    assert ( wrap_it != wrap_map.end() ) ;
    batch_options.twrap = OIIO::Tex::Wrap ( wrap_it->second ) ;
  
    auto mip_it = mipmode_map.find ( args.mip ) ;
    assert ( mip_it != mipmode_map.end() ) ;
    batch_options.mipmode = OIIO::Tex::MipMode ( mip_it->second ) ;
  
    auto interp_it = interpmode_map.find ( args.interp ) ;
    assert ( interp_it != interpmode_map.end() ) ;
    batch_options.interpmode = OIIO::Tex::InterpMode ( interp_it->second ) ;
  
    OIIO::ustring uenvironment ( args.input , 0 ) ;
    th = ts->get_texture_handle ( uenvironment ) ;
  }

  // set up the eval function.

  template < typename I , typename O >
  void eval ( const I & crd9 , O & px )
  {
    // we convert to OIIO-compatible coordinates as we go, this
    // requires changing the sign of the x and y axis (hence the
    // seemingly wrong order of the differences for ds[0], ds[1],
    // dt[0] and dt[1]). This is the only place where we need to
    // do this in 'extract', because only here we are feeding 3D
    // ray coordinates directly to OIIO code. OIIO's axis convention
    // is the same as Imath's. They have the positive x axis pointing
    // left and the positive y axis pointing down, whereas lux
    // convention, which I use in the remainder of the program, has
    // the positive x axis pointing right and the positive y axis
    // pointing down ('latin book order'). Both systems have the
    // positive z axis pointing forward.

    crd3_v c3 { - crd9[0] ,  // note the sign change to move to
                - crd9[1] ,  // OIIO-compatible 3D ray coordinates
                  crd9[2] } ;

    crd3_v ds { crd9[0] - crd9[3] ,    // sic
                crd9[1] - crd9[4] ,    // sic
                crd9[5] - crd9[2] } ;

    crd3_v dt { crd9[0] - crd9[6] ,    // sic
                crd9[1] - crd9[7] ,    // sic
                crd9[8] - crd9[2] } ;

    // now we can call 'environment', but depending on the SIMD
    // back-end, we provide the pointers which 'environment' needs
    // in different ways. The first form would in fact work for all
    // back-ends, but the second form is more concise. Note that,
    // as of this writing, the OIIO code accepts the batched arguments
    // but then proceeds to loop over them with single-point lookups,
    // which kind of defeats the purpose of a SIMDized pipeline. But
    // 'on this side' we're doing 'proper SIMD', and hope that OIIO
    // will also eventually provide it.

#if defined USE_VC or defined USE_STDSIMD

    // to interface with zimt's Vc and std::simd backends, we need to
    // extract the data from the SIMDized objects and re-package the
    // ouptut as a SIMDized object. The compiler will likely optimize
    // this away and work the entire operation in registers, so let's
    // call this a 'semantic manoevre'.

    float scratch [ 4 * nchannels * LANES ] ;

    c3.store ( scratch ) ;
    ds.store ( scratch + nchannels * LANES ) ;
    dt.store ( scratch + nchannels * LANES ) ;

    ts->environment ( th , nullptr, batch_options ,
                      OIIO::Tex::RunMaskOn ,
                      scratch ,
                      scratch + nchannels * LANES ,
                      scratch + 2 * nchannels * LANES ,
                      nchannels ,
                      scratch + 3 * nchannels * LANES ) ;

    px.load ( scratch + 3 * nchannels * LANES ) ;

#else

    // the highway and zimt's own backend have an internal representation
    // as a C vector of fundamentals, so we van use data() on them, making
    // the code even simpler - though the code above would work just the
    // same.

    ts->environment ( th , nullptr, batch_options ,
                      OIIO::Tex::RunMaskOn ,
                      c3[0].data() ,
                      ds[0].data() ,
                      dt[0].data() ,
                      nchannels ,
                      px[0].data() ) ;

#endif

    // and that's us done! the 'environment' function has returned pixel
    // values, which have been passed as result to the zimt::transform
    // process invoking this functor.
  }
} ;

// functor to convert between pixels with different channel count.
// currently the only case occuring in this program is the very first
// one: the output will be the same as the input and the code is
// optimized away. This object is not well-tested. currently unused.

/*
template < typename T , std::size_t A , std::size_t B , std::size_t L >
struct repix
: public zimt::unary_functor < zimt::xel_t < T , A > ,
                               zimt::xel_t < T , B > ,
                               L >
{
  typedef zimt::unary_functor < zimt::xel_t < T , A > ,
                                zimt::xel_t < T , B > ,
                                L >
          base_t ;

  using typename base_t::in_v ;
  using typename base_t::out_v ;

  void eval ( const in_v & in , out_v & out )
  {
    if constexpr ( A == B )
    {
      out = in ;
    }
    else if constexpr ( A == 1 )
    {
      out = in[0] ;
      if constexpr ( B == 4 )
        out[3] = 1.0 ;
    }
    else if constexpr ( B == 1 )
    {
      if constexpr ( A < 4 )
        out = in.sum() / float ( A ) ;
      else
      {
        out[0] = in[0] / in[4] ;
        out[0] += in[1] / in[4] ;
        out[0] += in[2] / in[4] ;
        out[0] /= 3.0f ;
        out[0] ( in[4] == 0 ) = 0 ;
      }
    }
    else if constexpr ( A == 4 && B == 3 )
    {
      out[0] = in[0] ;
      out[1] = in[1] ;
      out[2] = in[2] ;
      // TODO: should be (associated alpha)
      // out[0] = in[0] / in[4] ;
      // out[1] = in[1] / in[4] ;
      // out[2] = in[2] / in[4] ;
      // out[0] ( in[4] == 0 ) = 0 ;
      // out[1] ( in[4] == 0 ) = 0 ;
      // out[2] ( in[4] == 0 ) = 0 ;
    }
    else if constexpr ( A == 3 && B == 4 )
    {
      out[0] = in[0] ;
      out[1] = in[1] ;
      out[2] = in[2] ;
      out[4] = 1 ;
    }
  }
} ;
*/

// sixfold_t provides an internal representation of a cubemap
// with widened support, to provide for easy mip-mapping and
// interpolation.

template < std::size_t nchannels >
struct sixfold_t
: public metrics_t
{
  // shorthand for pixels and SIMDized pixels

  typedef zimt::xel_t < float , 2 > crd2_t ;
  typedef zimt::xel_t < float , nchannels > px_t ;
  typedef zimt::simdized_type < px_t , LANES > px_v ;

  // pointers to an OIIO texture system and an OIIO texture handle.
  // These are only used if the pick-up is done using OIIO's
  // 'texture' function.

  std::filesystem::path texture_file ;
  ts_ptr_t ts ;
  OIIO::TextureSystem::TextureHandle * th ;

  // This array holds all the image data. I'll refer to this
  // array as the 'IR image': the internal representation.
  // The array's size is gleaned from the base class metrics_t,
  // which is initialized first, so we can already use the
  // members we inherit from it (here: section_px, the width
  // of the IR image).

  zimt::array_t < 2 , px_t > store ;

  // pointer to the upper left corner of the topmost cubeface
  // image inside the 'store' array

  px_t * p_ul ;

  typedef zimt::bspline < px_t , 2 > spl_t ;
  std::shared_ptr < spl_t > p_bsp ;
  zimt::grok_type < crd2_t , px_t , LANES > bsp_ev ;

  sixfold_t ( std::size_t _face_px ,
              double _face_fov = M_PI_2 ,
              std::size_t _support_min_px = 4UL ,
              std::size_t _tile_px = 64UL )
  : metrics_t ( _face_px ,
                _face_fov ,
                _support_min_px ,
                _tile_px ) ,
    store ( { section_px , 6 * section_px } )
  {
    // let the user know the field of view of IR sections. This
    // is only possible if the cube face images have even size,
    // otherwise the cube face image's center can't coincide
    // with a section's center - unless the section is also
    // odd-sized, which only occurs when the tile size is one.

    if ( verbose )
    {
      if ( left_frame_px != right_frame_px )
        std::cout << "cube face is not centered in IR" << std::endl ;
      else
        std::cout << "IR sections have "
                  << atan ( section_md / 2.0 ) * 360.0 / M_PI
                  << " degrees fov" << std::endl ;
    }

    // find the location of the first cube face's upper left corner
    // in the 'store' array

    p_ul = store.data() ;
    p_ul += ( left_frame_px * store.strides ) . sum() ;
    
    auto * ps = new spl_t ( store , store , args.spline_degree ,
                            { zimt::REFLECT , zimt::REFLECT } ,
                            -1 , left_frame_px ) ;
    p_bsp.reset ( ps ) ;
    bsp_ev = zimt::make_evaluator < spl_t , float , LANES > ( *p_bsp ) ;
  }

  // function to read the cube face image data from disk (via
  // an OIIO-provided inp) - this version reads a single image
  // with 1:6 aspect ratio containing six square cube face
  // images, concatenated vertically, in openEXR order (left,
  // right, top, bottom, front, back) and orientation (The top
  // and bottom images are oriented so that they align vertically
  // with the 'back' image)
  // we can load cube maps with arbitrary field of view (as long
  // as the field of view is at least ninety degrees). 'Normal'
  // cubemaps will hold cube face images with precisely ninety
  // degrees field of view, measured either center-to-center
  // or edge-to-edge (see 'ctc'), but the scope of the 'load'
  // functions is wider - the face_fov argument determines the
  // precise extent.
  // The six cube face images are separated and each one is
  // stored in the corresponding section - right in the center
  // for even cube face sizes, and as central as possible for
  // odd cube face sizes - so the data are unmodified, but
  // shifted to their place in the IR image. After the load,
  // the frame of support is filled with interpolated values.

  void load ( const std::unique_ptr<ImageInput> & inp )
  {
    assert ( inp != nullptr ) ;

    const ImageSpec &spec = inp->spec() ;
    int xres = spec.width ;
    int yres = spec.height ;

    // the calling code should already have looked at the image and
    // gleaned the face width, but here we check again:

    assert ( xres == face_px ) ;
    assert ( yres == 6 * face_px ) ;

    // read the six cube face images from the 1:6 stripe into the
    // appropriate slots in the sixfold_t object's 'store' array.

    if ( inp->supports ( "scanlines" ) )
    {
      if ( verbose )
      {
        std::cout << "input supports scanline-based access"
        << std::endl ;
      }

      // if the input is scanline-based, we copy batches of
      // scanlines into the cube face slots in the store

      for ( int face = 0 ; face < 6 ; face++ )
      {
        auto * p_trg = p_ul + face * offset_px ;
      
        // for reference: OIIO's read_scanlines' signature
      
        // virtual bool read_scanlines ( int subimage, int miplevel,
        //                               int ybegin, int yend,
        //                               int z, int chbegin, int chend,
        //                               TypeDesc format, void *data,
        //                               stride_t xstride = AutoStride,
        //                               stride_t ystride = AutoStride)
      
        // note how we read face_px scanlines in one go, using
        // appropriate strides to place the image data inside the
        // larger 'store' array, converting to float as we go along.
        // The channels are capped at nchannels. We ask OIIO to
        // provide float data.
      
        // TODO: read_scanlines doesn't work with tile-based files
      
        auto success =
        inp->read_scanlines ( 0 , 0 ,
                              face * face_px , (face+1) * face_px ,
                              0 , 0 , nchannels ,
                              TypeDesc::FLOAT , p_trg ,
                              nchannels * 4 ,
                              nchannels * 4 * store.strides[1] ) ;
        assert ( success ) ;
      }
    }
    else
    {
      if ( verbose )
      {
        std::cout << "input is tiled" << std::endl ;
      }

      // if the input is not scanline-based, we read the entire image
      // into a buffer, then copy the cube faces to the store. We can
      // use a 3D array as a target, with the six cube faces 'on
      // top of each other' populating the third spatial dimension.
      // We need a large-ish amount of memory temporarily, but it's
      // easier that way - if memory were an issue we'd have to work
      // through the tiles and split them up to separate cube faces
      // if the cube face images' size is not a multiple of the
      // tile size.

      zimt::array_t < 3 , px_t >
        buffer ( { face_px , face_px , std::size_t ( 6 ) } ) ;
      
      // and a view to the store with the same shape, but strides to
      // match the metrics of the target memory area in the store
        
      zimt::view_t < 3 , px_t >
        target ( p_ul ,
                { 1L , long(section_px) , long(offset_px) } ,
                { face_px , face_px , std::size_t ( 6 ) } ) ;

      // read the image data into the buffer

      inp->read_image ( 0 , 0 , 0 , nchannels ,
                        TypeDesc::FLOAT ,
                        buffer.data() ) ;

      // zimt handles the data transfer from the buffer to the view

      target.copy_data ( buffer ) ;
    }
  }

  // load six separate cube face images. It's assumed that all
  // images exist and have the same geometry. This is essentially
  // the same as the function above, but the code is simplified
  // to always read an entire cube face image in one go.

  void load ( const std::vector < std::string > & filename6 )
  {
    // make sure there are six file names

    assert ( filename6.size() == 6 ) ;

    // set up a buffer for a single cube face image

    zimt::array_t < 2 , px_t >
      buffer ( { face_px , face_px } ) ;

    // read the six cube face images into the appropriate slots in
    // the sixfold_t object's 'store' array.

    for ( std::size_t face = 0 ; face < 6 ; face++ )
    {
      auto inp = ImageInput::open ( filename6[face] ) ;
      assert ( inp != nullptr ) ;

      if ( verbose )
        std::cout << "load cube face image " << filename6[face]
                  << std::endl ;

      const ImageSpec &spec = inp->spec() ;
      std::size_t w = spec.width ;
      std::size_t h = spec.height ;

      assert ( w == h ) ;
      assert ( w == face_px ) ;

      inp->read_image ( 0 , 0 , 0 , nchannels ,
                        TypeDesc::FLOAT ,
                        p_ul + face * offset_px ,
                        store.strides[0] * nchannels * 4 ,
                        store.strides[1] * nchannels * 4 ) ;

      inp->close() ;
    }
  }

  // We can use the fall-back bilinear interpolation to pick up
  // pixel values from the IR image, but to use OIIO's 'texture'
  // function, we need access to the texture system and - for
  // speed - the texture handle. But we don't have these data
  // when the sixfold_t is created - the texture has to be
  // generated first, then stored to disk and then fed to the
  // texture system. So we can only call this function later:

  void gen_texture ( const std::string & ts_options )
  {
    // to load the texture with OIIO's texture system code, it
    // has to be in a file, so we store the IR image to a temporary
    // file:

    auto temp_path = std::filesystem::temp_directory_path() ;
    texture_file = temp_path / "temp_texture.exr" ;

    if ( verbose )
      std::cout << "saving generated texture to "
                << texture_file.string() << std::endl ;

    save_array ( texture_file.string() , store ) ;

    // now we can create to the texture system and receive
    // a texture handle for the temporary file for fast access

    ts = OIIO::TextureSystem::create() ;

    if ( ts_options != std::string() )
    {
      if ( verbose )
        std::cout << "adding texture system options: " << ts_options
                  << std::endl ;
    
      ts->attribute ( "options" , ts_options ) ;
    }

    OIIO::ustring uenvironment ( texture_file.string() , 0 ) ;
    th = ts->get_texture_handle ( uenvironment ) ;
  }

  // given the cube face and the in-face coordinate, extract the
  // corresponding pixel values from the internal representation
  // of the cubemap held in the sixfold_t object. We have two
  // variants here, the first one using nearest neighbour
  // interpolation, the second bilinear. The first one is currently
  // unused. The pick-up is not guaranteed to look up pixel data
  // strictly inside the 90-degree cube face 'proper' but may
  // glean some information from the support frame, so this has
  // to be present. Initially we provide a one-pixel-wide
  // support frame of mirrored pixels (if necessary - if the
  // incoming partial images have 'inherent support' because
  // they span more than ninety degrees, this is not necessary)
  // which is enough for bilinear interpolation. Once we have
  // filled in the support frame, we can use interpolators with
  // wider support.

  // TODO: really, this member function could also be labeled
  // 'const', because it does not modify anything. But the
  // contained functor bsp_ev does not have a const eval.

  void cubemap_to_pixel ( const i_v & face ,
                          const crd2_v & in_face ,
                          px_v & px )
  {
    crd2_v pickup ;
    get_pickup_coordinate_px ( face , in_face , pickup ) ;
    bsp_ev.eval ( pickup , px ) ;
  }

  // this is the equivalent function to 'cubemap_to_pixel', above,
  // using OIIO's 'texture' function to perform the gleaning of
  // pixel data from the IR image. The IR image is not accessed
  // directly, but provided via OIIO's texture system.
  // On top of the pick-up coordinate (coming in in texture
  // units in [0,1], rather than pixel units which are used in
  // 'cubemap_to_pixel') we have the four derivatives and the
  // OIIO texture batch options.

  void get_filtered_px ( const crd2_v & pickup ,
                         px_v & px ,
                         const f_v & dsdx ,
                         const f_v & dtdx ,
                         const f_v & dsdy ,
                         const f_v & dtdy ,
                         OIIO::TextureOptBatch batch_options ) const
  {
    assert ( all_of ( pickup[0] >= 0.0f ) ) ;
    assert ( all_of ( pickup[0] <= 1.0f ) ) ;
    assert ( all_of ( pickup[1] >= 0.0f ) ) ;
    assert ( all_of ( pickup[1] <= 1.0f ) ) ;

    // code to truncate the pickup coordinate to int and gather,
    // after restoring the pickup to pixel units. Uncomment this
    // code to directly pick up pixels from the IR image with NN
    // interpolation - this is to verify that the pick-up coordinate
    // is indeed correct.

//     auto x = pickup[0] * float ( store.shape[0] ) ;
//     auto y = pickup[1] * float ( store.shape[1] ) ;
//     index_v idx { x , y } ;
//     const auto ofs = ( idx * store.strides ) . sum() * nchannels ;
//     const auto * p = (float*) ( store.data() ) ;
//     
//     px.gather ( p , ofs ) ;
//     return ;

    // OIIO's 'texture' signature is quite a mouthful:

    // virtual bool texture (    TextureHandle *texture_handle,
                              // Perthread *thread_info,
                              // TextureOptBatch &options,
                              // Tex::RunMask mask,
                              // const float *s, const float *t,
                              // const float *dsdx, const float *dtdx,
                              // const float *dsdy, const float *dtdy,
                              // std::size_t nchannels,
                              // float *result,
                              // float *dresultds = nullptr,
                              // float *dresultdt = nullptr)

    #if defined USE_VC or defined USE_STDSIMD

    // to interface with zimt's Vc and std::simd backends, we need to
    // extract the data from the SIMDized objects and re-package the
    // ouptut as a SIMDized object. The compiler will likely optimize
    // this away and work the entire operation in registers, so let's
    // call this a 'semantic manoevre'.

    float scratch [ 6 * LANES + nchannels * LANES ] ;

    pickup.store ( scratch ) ; // stores 2 * LANES
    dsdx.store ( scratch + 2 * LANES ) ;
    dtdx.store ( scratch + 3 * LANES ) ;
    dsdy.store ( scratch + 4 * LANES ) ;
    dtdy.store ( scratch + 5 * LANES ) ;

    bool result =
    ts->texture ( th , nullptr , batch_options , OIIO::Tex::RunMaskOn ,
                  scratch , scratch + LANES ,
                  scratch + 2 * LANES , scratch + 3 * LANES ,
                  scratch + 4 * LANES , scratch + 5 * LANES ,
                  nchannels , scratch + 6 * LANES ) ;

    assert ( result ) ;
    px.load ( scratch + 6 * LANES ) ;

    #else

    // zimt's own and the highway backend have a representation as
    // a C vector of fundamentals and provide a 'data' function
    // to yield it's address. This simplifies matters, we can pass
    // these pointers to OIIO directly.

    bool result =
    ts->texture ( th , nullptr , batch_options , OIIO::Tex::RunMaskOn ,
                  pickup[0].data() , pickup[1].data() ,
                  dsdx.data() , dtdx.data() ,
                  dsdy.data() , dtdy.data() ,
                  nchannels , (float*) ( px[0].data() ) ) ;

    #endif
  }

  // variant taking derivatives of the in-face coordinate, which
  // are approximated by calculating the difference to a canonical
  // (target image) coordinate one sample step to the right (x)
  // or below (y), respectively. The derivatives are in texture
  // units aready, and we also convert the pickup coordinate to
  // texture units.

  void cubemap_to_pixel ( const i_v & face ,
                          crd2_v in_face ,
                          px_v & px ,
                          const f_v & dsdx ,
                          const f_v & dtdx ,
                          const f_v & dsdy ,
                          const f_v & dtdy ,
                          OIIO::TextureOptBatch & bo ) const
  {
    crd2_v pickup_tx ;

    // obtain the pickup coordinate in texture units (in [0,1])

    get_pickup_coordinate_tx ( face , in_face , pickup_tx ) ;

    // use OIIO to get the pixel value

    get_filtered_px ( pickup_tx , px , dsdx , dtdx , dsdy , dtdy , bo ) ;
  }

  // After the cube faces have been read from disk, they are surrounded
  // by black (or even undefined) pixels. We want to provide minimal
  // support, this support's quality is not crucial, but it should not
  // be black, but rather like mirroring on the edge, which this function
  // does - it produces a one-pixel-wide frame with mirrored pixels
  // around each of the cube faces. In the next step, we want to fill
  // in the support frame around the cube faces proper, and we may have
  // to access image data close to the margin. Rather than implementing
  // the mirroring on the edge by manipulating the coordinate of the
  // pick-up (e.g. clamping it) having the one-pixel-wide minimal
  // support allows us to do the next stage without looking at the
  // coordinates: we can be sure that the pick-up will not exceed
  // the support area.

  void mirror_around()
  {
    auto * p_base = store.data() ;

    for ( int face = 0 ; face < 6 ; face++ )
    {
      // get a pointer to the upper left of the cube face 'proper'

      auto * p_frame = p_base + face * section_px * store.strides[1]
                              + left_frame_px * store.strides[1]
                              + left_frame_px * store.strides[0] ;

      // get a zimt view to the current cube face

      zimt::view_t < 2 , px_t > cubeface
        ( p_frame , store.strides , { face_px , face_px } ) ;

      // we use 2D discrete coordinates

      typedef zimt::xel_t < int , 2 > ix_t ;

      // mirror the horizontal edges

      for ( int x = -1 ; x <= int ( face_px ) ; x++ )
      {
        ix_t src { x ,  0 } ;
        ix_t trg { x , -1 } ;
        cubeface [ trg ] = cubeface [ src ] ;
        src [ 1 ] = face_px - 1 ;
        trg [ 1 ] = face_px ;
        cubeface [ trg ] = cubeface [ src ] ;
      }

      // and the vertical edges

      for ( int y = -1 ; y <= int ( face_px ) ; y++ )
      {
        ix_t src {  0 , y } ;
        ix_t trg { -1 , y } ;
        cubeface [ trg ] = cubeface [ src ] ;
        src [ 0 ] = face_px - 1 ;
        trg [ 0 ] = face_px ;
        cubeface [ trg ] = cubeface [ src ] ;
      }
    }
  }

  // We set up an internal functor which we'll use to fill the frame
  // of support pixels. This functor is local to struct sixfold_t
  // because it has no practical use outside of this context and
  // is tailor-made for the purpose (the member function fill_support
  // just below it)

  // this functor is used to fill the frame of support pixels in the
  // array in the sixfold_t object. incoming, we have 2D coordinates,
  // which, for the purpose at hand, will lie outside the cube face.
  // But we can still convert these image coordinates to planar
  // coordinates in 'model space' and then further to 'pickup'
  // coordinates into the array in the sixfold_t object. When we
  // pick up data from the neighbouring cube faces, we produce
  // support around the cube face proper, 'regenerating' what would
  // have been there in the first place, if we had had partial images
  // with larger field of view. Due to the geometric transformation
  // and interpolation, the regenerated data in the frame will not
  // be 'as good' as genuine image data, but we'll never actually
  // 'look at' the regenerated data: they only serve as support for
  // filtering and mip-mapping, and for that purpose, they are certainly
  // 'good enough'.
  // The eval functor could of course also produce pixel values for
  // 2D coordinates inside the cube face image. Note that we have 
  // discrete incoming coordinates - the very coordinates of pixels
  // inside the frame areas which need filling.in

  struct fill_frame_t
  : public zimt::unary_functor
      < v2i_t , zimt::xel_t < float , nchannels > , LANES >
  {
    const sixfold_t < nchannels > & sf ;
    const int face ;
    const int ithird ;

    typedef zimt::xel_t < float , nchannels > px_t ;
    typedef zimt::simdized_type < px_t , LANES > px_v ;

    // we use a degree-1 b-spline (bilinear interpolation) to
    // fill the support frame
    // TODO: when specializing with 1, invocation with
    // --spline_degree 0 crashes.

    zimt::evaluator < crd2_t , px_t , LANES , 0 > ev ;

    // note the factor of two in the initialization of 'ithird':
    // incoming coordinates are doubled (!) for the purpose at hand.

    fill_frame_t ( const sixfold_t < nchannels > & _cubemap ,
                   const int & _face )
    : sf ( _cubemap ) ,
      face ( _face ) ,
      ithird ( _cubemap.model_to_px * 2 ) ,
      ev ( * _cubemap.p_bsp )
    { }

    void eval ( const v2i_v & crd2 , px_v & px ) const
    {
      crd3_v crd3 ;

      // since 'face' is const, this case switch should be optimized away.
      // here, we move from discrete image coordinates to float 3D
      // coordinates, but we don't scale to model coordinates and do the
      // arithmetic in integer until the final conversion to float:
      // the scale does not matter, because the next processing step
      // is insensitive to scale.
      // Note that we're obtaining crd2 values from a linspace_t which
      // provides readily shifted and doubled (!) coordinates, hence the
      // factor of two in the initialization of ithird.

      switch ( face )
      {
        case CM_FRONT :
          crd3[RIGHT]   =   crd2[RIGHT] ;
          crd3[DOWN]    =   crd2[DOWN] ;
          crd3[FORWARD] =   ithird ;
          break ;
        case CM_BACK :
          crd3[RIGHT]   = - crd2[RIGHT] ;
          crd3[DOWN]    =   crd2[DOWN] ;
          crd3[FORWARD] = - ithird ;
          break ;
        case CM_RIGHT :
          crd3[RIGHT] =     ithird ;
          crd3[DOWN] =      crd2[DOWN] ;
          crd3[FORWARD] = - crd2[RIGHT] ;
          break ;
        case CM_LEFT :
          crd3[RIGHT] =   - ithird ;
          crd3[DOWN] =      crd2[DOWN] ;
          crd3[FORWARD] =   crd2[RIGHT] ;
          break ;
        // for bottom and top, note that we're using openEXR convention.
        // to use lux convention, invert the signs.
        case CM_BOTTOM :
          crd3[RIGHT] =   - crd2[RIGHT] ;
          crd3[DOWN] =      ithird ;
          crd3[FORWARD] =   crd2[DOWN] ;
          break ;
        case CM_TOP :
          crd3[RIGHT] =   - crd2[RIGHT] ;
          crd3[DOWN] =    - ithird ;
          crd3[FORWARD] = - crd2[DOWN] ;
          break ;
      } ;

      // next we use the 3D coordinates we have just obtained to
      // find a source cube face and in-face coordinates - this will
      // be a different cube face to 'face', and it will be one
      // which actually can provide data for location we're filling
      // in. So we re-project the content from adjoining cube faces
      // into the support area around the cube face we're currently
      // surrounding with a support frame.
      // Note how we use sixfold_t's member functions for most of
      // the work - gleaning pixel values from the sixfold_t object
      // works the same for both purposes. The difference here is
      // the source of the 3D ray coordinates: here they originate
      // from target locations on the support frame, whereas the
      // eval member function in ll_to_px_t produces them from
      // incoming lat/lon coordinates.

      i_v fv ;
      crd2_v in_face ;
      ray_to_cubeface < float , LANES > ( crd3 , fv , in_face ) ;
      crd2_v pickup ;
      sf.get_pickup_coordinate_px ( fv , in_face , pickup ) ;

      // finally we use this information to obtain pixel values,
      // which are written to the target location.

      ev.eval ( pickup , px ) ;
    }
  } ;

  // fill_support uses the fill_frame_t functor to populate the
  // frame of support. The structure of the code is similar to
  // 'mirror_around', iterating over the six sections of the
  // array and manipulating each in turn. But here we fill in
  // the entire surrounding frame, not just a pixel-wide line,
  // and we pick up data from neighbouring cube faces.

  void fill_support()
  {
    if ( left_frame_px == 0 && right_frame_px == 0 )
      return ;

    auto * p_base = store.data() ;

    for ( int face = 0 ; face < 6 ; face++ )
    {
      // set up the 'gleaning' functor

      fill_frame_t fill_frame ( *this , face ) ;
    
      // get a pointer to the upper left of the section of the IR

      auto * p_frame = p_base + face * section_px * store.strides[1] ;
      
      // we form a view to the current section of the array in the
      // sixfold_t object. The shape of this view is the 'notional
      // shape' zimt::process will work with.

      zimt::view_t < 2 , px_t >
        section ( p_frame , store.strides ,
                { section_px , section_px } ) ;

      auto & shp = section.shape ;

      // we'll use a 'loading bill' to narrow the filling-in down
      // to the areas which are outside the cube face.

      zimt::bill_t bill ;
      
      // now we fill in the lower and upper limits. This is a good
      // demonstration of how these parameters can be put to use.
      // We use a 'notional' shape which encompasses the entire
      // section, but the iteration will only visit those coordinates
      // which are in the range given by the limits.
      // our linspace_t will generate doubled coordinates -
      // scaling up by a factor of two is irrelevant for the next
      // stage of the processing, and it allows us to code the
      // sampling mathematics entirely in int until the stage
      // where we do the unavoidable division to project the
      // ray onto a plane. Without the doubling, we'd have to
      // start out at ( section_px -1 ) / 2, which can't be
      // expressed precisely in int.

      auto ishift = section_px - 1 ;
      zimt::linspace_t < int , 2 , 2 , LANES > ls ( -ishift , 2 ) ;

      zimt::storer < float , nchannels , 2 , LANES > st ( section ) ;
      
      // fill in the stripe above the cube face

      if ( left_frame_px > 0 )
      {
        bill.lower_limit = { 0 , 0 } ;
        bill.upper_limit = { long(section_px) , long(left_frame_px) } ;
        zimt::process ( shp , ls , fill_frame , st , bill ) ;
      }

      // fill in the stripe below the cube face

      if ( right_frame_px > 0 )
      {
        bill.lower_limit = { 0 , long(section_px - right_frame_px) } ;
        bill.upper_limit = { long(section_px) , long(section_px) } ;
        zimt::process ( shp , ls , fill_frame , st , bill ) ;
      }

      // fill in the stripe to the left of the cube face

      if ( left_frame_px > 0 )
      {
        bill.lower_limit = { 0 , long(left_frame_px) } ;
        bill.upper_limit = { long(left_frame_px) ,
                            long(section_px) - long(right_frame_px) } ;
        zimt::process ( shp , ls , fill_frame , st , bill ) ;
      }

      // fill in the stripe to the right of the cube face

      if ( right_frame_px > 0 )
      {
        bill.lower_limit = { long(left_frame_px + face_px) ,
                             long(left_frame_px) } ;
        bill.upper_limit = { long(section_px) ,
                             long(section_px) - long(right_frame_px) } ;
      }

      zimt::process ( shp , ls , fill_frame , st , bill ) ;
    }
  }

    // if we're using b-spline interpolation for the cubemap,
    // we need to prefilter the sections appropriately, TODO:
    // may be better in a spearate function

  void prefilter ( int spline_degree )
  {
    auto * p_base = store.data() ;

    if ( spline_degree > 1 )
    {
      for ( int face = 0 ; face < 6 ; face++ )
      {
        // get a pointer to the upper left of the section of the IR

        auto * p_frame = p_base + face * section_px * store.strides[1] ;
        
        // we form a view to the current section of the array in the
        // sixfold_t object. The shape of this view is the 'notional
        // shape' zimt::process will work with.

        zimt::view_t < 2 , px_t >
          section ( p_frame , store.strides ,
                  { section_px , section_px } ) ;

        zimt::prefilter ( section , section ,
                          { zimt::NATURAL , zimt::NATURAL } ,
                          spline_degree ) ;
      }
    }
  }

} ; // end of struct sixfold_t

// functor to obtain pixel values for ray coordinates from a cubemap

template < std::size_t nchannels >
struct cbm_to_px_t
: public zimt::unary_functor
   < zimt::xel_t < float , 3 > ,
     zimt::xel_t < float , nchannels > ,
     LANES >
{
  // some types

  typedef zimt::xel_t < float , 3 > crd3_t ;
  typedef zimt::simdized_type < crd3_t , LANES > crd3_v ;
  typedef zimt::xel_t < float , nchannels > px_t ;
  typedef zimt::simdized_type < px_t , LANES > px_v ;

  sixfold_t < nchannels > cubemap ;

  // cbm_to_px_t's c'tor obtains a const reference to the sixfold_t
  // object holding pixel data. It receives 3D ray coordinates and
  // produces pixel values gleaned with bilinear interpolation.

  cbm_to_px_t ( const sixfold_t < nchannels > & _cubemap )
  : cubemap ( _cubemap )
  { }

  void eval ( const crd3_v & crd3 , px_v & px )
  {
    // find the cube face and in-face coordinate for 'c'

    i_v face ;
    crd2_v in_face ;
    ray_to_cubeface < float , LANES > ( crd3 , face , in_face ) ;

    cubemap.cubemap_to_pixel ( face , in_face , px ) ;
  }

} ;

template < std::size_t nchannels >
struct cbm_to_px_t2
: public zimt::unary_functor
   < zimt::xel_t < float , 9 > ,
     zimt::xel_t < float , nchannels > ,
     LANES >
{
  // some types

  typedef zimt::xel_t < float , 2 > crd2_t ;
  typedef zimt::simdized_type < crd2_t , LANES > crd2_v ;
  typedef zimt::xel_t < float , 3 > crd3_t ;
  typedef zimt::simdized_type < crd3_t , LANES > crd3_v ;
  typedef zimt::xel_t < float , nchannels > px_t ;
  typedef zimt::simdized_type < px_t , LANES > px_v ;

  sixfold_t < nchannels > cubemap ;

  // for lookup with OIIO's texture system, we need batch options.
  // I'd keep them in the code which actually uses them, but the
  // OIIO 'texture' function expects an lvalue.

  OIIO::TextureOptBatch batch_options ;

  // scaling factor to move from model space units to texture units
  // (separate for the s and t direction - these factors are applied
  // to coordinates pertaining to the IR image of the cubemap, which
  // has 1:6 aspect ratio)

  const float scale_s ;
  const float scale_t ;

  // cbm_to_px_t's c'tor obtains a const reference to the sixfold_t
  // object holding pixel data. It receives 3D ray coordinates and
  // produces pixel values gleaned with bilinear interpolation.

  cbm_to_px_t2 ( const sixfold_t < nchannels > & _cubemap ,
                 const OIIO::TextureOptBatch & _batch_options )
  : cubemap ( _cubemap ) ,
    scale_s ( _cubemap.model_to_px / _cubemap.store.shape[0] ) ,
    scale_t ( _cubemap.model_to_px / _cubemap.store.shape[1] ) ,
    batch_options ( _batch_options )
  { }

  template < typename I , typename O >
  void eval ( const I & crd9 , O & px )
  {
    zimt::xel_t < float , 2 > scale { scale_s , scale_t } ;

    // extract the three 3D ray coordinates from the nine-pack

    crd3_v c3 { crd9[0] , crd9[1] , crd9[2] } ;
    crd3_v dx { crd9[3] , crd9[4] , crd9[5] } ;
    crd3_v dy { crd9[6] , crd9[7] , crd9[8] } ;

    // find the cube face and in-face coordinate for 'c3'

    i_v face ;
    crd2_v in_face_00 , in_face_10 , in_face_01 ;
    ray_to_cubeface < float , LANES > ( c3 , face , in_face_00 ) ;

    // do the same for the two offsetted rays
  
    ray_to_cubeface_fixed < float , LANES >
      ( dx , face , in_face_10 ) ;

    ray_to_cubeface_fixed < float , LANES >
      ( dy , face , in_face_01 ) ;

    // obtain the pickup coordinate in texture units (in [0,1])

    crd2_v c3_tx ;
    cubemap.get_pickup_coordinate_tx ( face , in_face_00 , c3_tx ) ;

    // convert the neighbouring rays to texture coordinates
    
    in_face_00 *= scale ;
    in_face_10 *= scale ;
    in_face_01 *= scale ;

    // form the differences

    crd2_v dx_tx = in_face_10 - in_face_00 ;
    crd2_v dy_tx = in_face_01 - in_face_00 ;

    // use OIIO to get the pixel value

    // OIIO wants the derivatives in this order:
    // float *dsdx, const float *dtdx,
    // const float *dsdy, const float *dtdy

    cubemap.get_filtered_px ( c3_tx , px ,
                              dx_tx[0] , dx_tx[1] ,
                              dy_tx[0] , dy_tx[1] ,
                              batch_options ) ;
  }
} ;

// the 'latlon' template codes objects which can serve as 'act'
// functor in zimt::process. It's coded as a zimt::unary_functor
// taking 'ninepack' coordinates and returning pixels with C
// channels

template < typename T , typename U , std::size_t C , std::size_t L >
class latlon
: public zimt::unary_functor < zimt::xel_t < T , 9 > ,
                               zimt::xel_t < U , C > ,
                               L >
{
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , C > > > pa ;

public:
  
  typedef zimt::xel_t < T , 9 > ray_t ;
  typedef zimt::xel_t < U , C > px_t ;
  typedef zimt::xel_t < std::size_t , 2 > shape_type ;

  zimt::grok_type < ray_t , px_t , L > env ;

  latlon()
  {
    auto & inp ( args.inp ) ;
    auto & w ( args.env_width ) ;
    auto & h ( args.env_height ) ;

    shape_type shape { w , h } ;

    assert ( w == 2 * h ) ;

    if ( verbose )
      std::cout << "input has 2:1 aspect ratio, assuming latlon"
                << std::endl ;

    typedef zimt::xel_t < float , C > in_px_t ;
    pa = std::make_shared < zimt::array_t < 2 , in_px_t > >
      ( shape ) ;

    bool success = inp->read_image ( 0 , 0 , 0 , C ,
                                    TypeDesc::FLOAT , pa->data() ) ;
    assert ( success ) ;

    env = zimt::grok_type < ray_t , px_t , L > ( eval_env<C>() ) ;
  }

  void eval ( const typename zimt::grok_type < ray_t , px_t , L >::in_v & in ,
              typename zimt::grok_type < ray_t , px_t , L >::out_v & out )
  {
    env.eval ( in , out ) ;
  }

} ;

template < typename T , typename U , std::size_t C , std::size_t L >
class cubemap
: public zimt::unary_functor < zimt::xel_t < T , 9 > ,
                               zimt::xel_t < U , C > ,
                               L >
{
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , C > > > pa4 ;

public:
  
  typedef zimt::xel_t < T , 9 > ray_t ;
  typedef zimt::xel_t < U , C > px_t ;
  typedef zimt::xel_t < std::size_t , 2 > shape_type ;

  zimt::grok_type < ray_t , px_t , L > env ;

  cubemap()
  {
    auto & inp ( args.inp ) ;
    auto & w ( args.env_width ) ;
    auto & h ( args.env_height ) ;

    shape_type shape { w , h } ;

    assert ( h == 6 * w ) ;

    OIIO::TextureOptBatch batch_options ;

    for ( int i = 0 ; i < 16 ; i++ )
      batch_options.swidth[i] = batch_options.twidth[i]
        = args.stwidth ;
    
    for ( int i = 0 ; i < 16 ; i++ )
      batch_options.sblur[i] = batch_options.tblur[i]
        = args.stblur ;

    batch_options.conservative_filter = args.conservative_filter ;

    auto wrap_it = wrap_map.find ( args.swrap ) ;
    assert ( wrap_it != wrap_map.end() ) ;
    batch_options.swrap = OIIO::Tex::Wrap ( wrap_it->second ) ;
  
    wrap_it = wrap_map.find ( args.twrap ) ;
    assert ( wrap_it != wrap_map.end() ) ;
    batch_options.twrap = OIIO::Tex::Wrap ( wrap_it->second ) ;
  
    auto mip_it = mipmode_map.find ( args.mip ) ;
    assert ( mip_it != mipmode_map.end() ) ;
    batch_options.mipmode = OIIO::Tex::MipMode ( mip_it->second ) ;
  
    auto interp_it = interpmode_map.find ( args.interp ) ;
    assert ( interp_it != interpmode_map.end() ) ;
    batch_options.interpmode = OIIO::Tex::InterpMode ( interp_it->second ) ;
  
    if ( verbose )
    {
      std::cout << "input has 1:6 aspect ratio, assuming cubemap"
                << std::endl ;
    }

    sixfold_t<C> sf ( args.env_width , args.cbmfov ,
                      args.support_min , args.tile_size ) ;
    if ( args.multiple_input )
      sf.load ( args.cfs.get_filenames() ) ;
    else
      sf.load ( inp ) ;
    inp->close() ;
    sf.mirror_around() ;
    sf.fill_support() ;
    sf.gen_texture ( args.tsoptions ) ;

    env =   cbm_to_px_t2 < C > ( sf , batch_options ) ;
  }

  void eval ( const typename zimt::grok_type < ray_t , px_t , L >::in_v & in ,
              typename zimt::grok_type < ray_t , px_t , L >::out_v & out )
  {
    env.eval ( in , out ) ;
  }

} ;

// source_t provides pixel values from a mounted 2D manifold. The
// incoming 2D coordinates are texture coordinates guaranteed to be
// in the range of [0,1].

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
  const float * const p = nullptr ;
  ts_ptr_t ts ;
  OIIO::TextureSystem::TextureHandle * th ;
  OIIO::TextureOptBatch batch_options ;

  // the c'tor with atguments receives a view to image data and extracts
  // the relevant values to access these data directly. This is the route
  // taken for direct bilinear interpolation of the source data. Note how
  // we cast data() to float to 'shed' the channel count from the type

  source_t ( const zimt::view_t < 2 , px_t > & src )
  : width ( src.shape [ 0 ] ) ,
    height ( src.shape [ 1 ] ) ,
    strides ( src.strides ) ,
    p ( (const float* const) src.data() )
  { }

  // The c'tor without arguments routes to code using OIIO's 'texture'
  // function to handle the lookup. Here, we rely on OIIO's TextureSystem
  // to provide the data - we don't need to set up image storage.

  source_t()
  {
    // itp == -1 means 'use OIIO's 'texture' function

    assert ( args.itp == -1 ) ;

    // transfer all OIIO-related options from the arguments. This might
    // be factored out (the same code is also used for cubemap processing
    // via an OIIO texture) - but one might consider using different
    // options for the two processes, so for now I just copy and paste.

    for ( int i = 0 ; i < 16 ; i++ )
      batch_options.swidth[i] = batch_options.twidth[i]
        = args.stwidth ;
    
    for ( int i = 0 ; i < 16 ; i++ )
      batch_options.sblur[i] = batch_options.tblur[i]
        = args.stblur ;

    batch_options.conservative_filter = args.conservative_filter ;

    auto wrap_it = wrap_map.find ( args.swrap ) ;
    assert ( wrap_it != wrap_map.end() ) ;
    batch_options.swrap = OIIO::Tex::Wrap ( wrap_it->second ) ;
  
    wrap_it = wrap_map.find ( args.twrap ) ;
    assert ( wrap_it != wrap_map.end() ) ;
    batch_options.twrap = OIIO::Tex::Wrap ( wrap_it->second ) ;
  
    auto mip_it = mipmode_map.find ( args.mip ) ;
    assert ( mip_it != mipmode_map.end() ) ;
    batch_options.mipmode = OIIO::Tex::MipMode ( mip_it->second ) ;
  
    auto interp_it = interpmode_map.find ( args.interp ) ;
    assert ( interp_it != interpmode_map.end() ) ;
    batch_options.interpmode = OIIO::Tex::InterpMode ( interp_it->second ) ;
  
    ts = OIIO::TextureSystem::create() ;

    if ( args.tsoptions != std::string() )
    {
      if ( verbose )
        std::cout << "adding texture system options: " << args.tsoptions
                  << std::endl ;
    
      ts->attribute ( "options" , args.tsoptions ) ;
    }

    OIIO::ustring uenvironment ( args.mount_image , 0 ) ;
    th = ts->get_texture_handle ( uenvironment ) ;
  }

  void eval ( const pkg_v & _crd , px_v & px ) const
  {
    // incoming is const &, hence:

    pkg_v crd ( _crd ) ;

    if constexpr ( ncrd == 2 )
    {
      // bilinear interpolation with clamping at the edges.
      // We rely on crd being in range [0,1], but the code at hand
      // can tolerate small (up to less than .5) overshoots.

      crd[0] *= width ;
      crd[1] *= height ;
      crd -= .5f ;

      // get the upper left neighbouring discrete coordinate

      v2i_v uli { floor ( crd[0] ) , floor ( crd[1] ) } ;

      // and derive weights for the bilinear interpolation

      auto wr = crd - uli ;
      auto wl = 1.0f - wr ;

      // we also want the other three neighbouring dicrete coordinates

      v2i_v uri ( uli ) ;
      uri[0] += 1 ;
      
      v2i_v lli ( uli ) ;
      lli[1] += 1 ;
      
      v2i_v lri ( lli ) ;
      lri[0] += 1 ;

      // clamp the values to [0,width-1] and [0,height-1]

      uli[0] = uli[0].at_least ( 0 ) ;
      lli[0] = lli[0].at_least ( 0 ) ;
      
      uli[1] = uli[1].at_least ( 0 ) ;
      uri[1] = uri[1].at_least ( 0 ) ;
      
      uri[0] = uri[0].at_most ( width - 1 ) ;
      lri[0] = lri[0].at_most ( width - 1 ) ;
      
      lli[1] = lli[1].at_most ( height - 1 ) ;
      lri[1] = lri[1].at_most ( height - 1 ) ;

      // now gather, apply weights, sum up

      // obtain the four constituents by first truncating their
      // coordinate to int and then gathering from p.

      zimt::xel_t < int , 2 > istrides ( strides ) ;
      int inchannels ( nchannels ) ;

      index_v idsdxl { uli[0] , uli[1] } ;
      auto ofs = ( idsdxl * istrides ) . sum() * inchannels ;
      px_v pxul ;
      pxul.gather ( p , ofs ) ;

      index_v idsdxr { uri[0] , uri[1] } ;
      ofs = ( idsdxr * istrides ) . sum() * inchannels ;
      px_v pxur ;
      pxur.gather ( p , ofs ) ;

      index_v idxll { lli[0] , lli[1] } ;
      ofs = ( idxll * istrides ) . sum() * inchannels ;
      px_v pxll ;
      pxll.gather ( p , ofs ) ;

      index_v idxlr { lri[0] , lri[1] } ;
      ofs = ( idxlr * istrides ) . sum() * inchannels ;
      px_v pxlr ;
      pxlr.gather ( p , ofs ) ;

      // apply the bilinear formula with the weights gleaned above

      px  = wl[1] * ( wl[0] * pxul + wr[0] * pxur ) ;
      px += wr[1] * ( wl[0] * pxll + wr[0] * pxlr ) ;
    }
    else // to ( ncrd == 2 )
    {
      // use OIIO's 'texture' function to extract the pixel value

      // we have six incoming coordinate values, namely the pick-up
      // coordinate itself and two more derived from the target
      // coordinate's right and lower neighbour. We form an approximation
      // of the derivative by differencing. This is not strictly correct
      // (see the remarks about forming two vectors coplanar to the
      // tangent plane at 'crd' elsewhere) but it's 'good enough'.

      crd_v dx_tx = { abs ( crd[2] - crd[0] ) , abs ( crd[3] - crd[1] ) } ;
      crd_v dy_tx = { abs ( crd[4] - crd[0] ) , abs ( crd[5] - crd[1] ) } ;

      // OIIO's 'texture' functions requires a non-const batch_options
      // object, the compiler should take care of that with copy elision

      OIIO::TextureOptBatch _batch_options ( batch_options ) ;

      // for correct wrap-around, overly large values are assumed to
      // originate from opposite edges of the texture. (note: set twrap,
      // swrap for 360X180 to avoid black line at the wrap-around point)

      // TODO: think about this some more, re. different projections
      // note: with 360 degree fisheyes we get flawed output near the
      // 'back' pole with standard OIIO pickup.

      dx_tx[0] ( dx_tx[0] > .9f ) = 1.0f - dx_tx[0] ;
      dx_tx[1] ( dx_tx[1] > .9f ) = 1.0f - dx_tx[1] ;
      dy_tx[0] ( dy_tx[0] > .9f ) = 1.0f - dy_tx[0] ;
      dy_tx[1] ( dy_tx[1] > .9f ) = 1.0f - dy_tx[1] ;

      #if defined USE_VC or defined USE_STDSIMD

      // to interface with zimt's Vc and std::simd backends, we need to
      // extract the data from the SIMDized objects and re-package the
      // ouptut as a SIMDized object. The compiler will likely optimize
      // this away and work the entire operation in registers, so let's
      // call this a 'semantic manoevre'.

      float scratch [ 6 * LANES + nchannels * LANES ] ;

      crd[0].store ( scratch ) ;
      crd[1].store ( scratch + LANES ) ;
      dx_tx[0].store ( scratch + 2 * LANES ) ;
      dx_tx[1].store ( scratch + 3 * LANES ) ;
      dy_tx[0].store ( scratch + 4 * LANES ) ;
      dy_tx[1].store ( scratch + 5 * LANES ) ;

      bool result =
      ts->texture ( th , nullptr , _batch_options , OIIO::Tex::RunMaskOn ,
                    scratch , scratch + LANES ,
                    scratch + 2 * LANES , scratch + 3 * LANES ,
                    scratch + 4 * LANES , scratch + 5 * LANES ,
                    nchannels , scratch + 6 * LANES ) ;

      assert ( result ) ;
      px.load ( scratch + 6 * LANES ) ;

      #else

      // zimt's own and the highway backend have a representation as
      // a C vector of fundamentals and provide a 'data' function
      // to yield it's address. This simplifies matters, we can pass
      // these pointers to OIIO directly.

      bool result =
      ts->texture ( th , nullptr , _batch_options , OIIO::Tex::RunMaskOn ,
                    crd[0].data() , crd[1].data() ,
                    dx_tx[0].data() , dx_tx[1].data() ,
                    dy_tx[0].data() , dy_tx[1].data() ,
                    nchannels , (float*) ( px[0].data() ) ) ;

      assert ( result ) ;
      #endif
    }
  }
} ;

// struct mount_t provides data from a rectanglular 2D manifold holding
// pixel data in a given projection which do not cover the entire
// 360X180 degree environment. Typical candidates would be rectilinear
// patches or cropped images. This class handles channel counts up to
// four and paints pixels outside the covered range black. Four RGBA
// pixels, 'outside' pixels are painted 0000, assuming associated alpha.
// The 'ncrd' template argument is either three for lookups without
// or nine for lookups with two next neighbour's coordinates which
// can be used to find derivatives. This class relies on an 'inner'
// functor of class source_t to actually provide pixel data, this here
// class deals with the geometrical transformations needed for the
// different types of projections which the mounted images can have,
// and with masking out parts where the mounted image has no data.

template < std::size_t nchannels , std::size_t ncrd ,
           projection_t P , std::size_t L >
struct mount_t
: public zimt::unary_functor < zimt::xel_t < float , ncrd > ,
                               zimt::xel_t < float , nchannels > ,
                               L
                             >
{
  typedef zimt::unary_functor < zimt::xel_t < float , ncrd > ,
                                zimt::xel_t < float , nchannels > ,
                                L
                              > base_t ;

  typedef zimt::xel_t < float , ncrd > ray_t ;
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

  extent_type extent ;
  crd_t center ;
  crd_t rgirth ;
  source_t < nchannels , 2 * ncrd / 3 , L > & inner ;

  static_assert ( ncrd == 3 || ncrd == 9 ) ;

  static_assert (    P == SPHERICAL
                  || P == CYLINDRICAL
                  || P == RECTILINEAR
                  || P == STEREOGRAPHIC
                  || P == FISHEYE ) ;

  mount_t ( extent_type _extent ,
            source_t < nchannels , 2 * ncrd / 3 , L > & _inner )
  : extent ( _extent ) ,
    center { ( _extent.x0 + _extent.x1 ) * .5 ,
             ( _extent.y0 + _extent.y1 ) * .5 } ,
    rgirth { 1.0  / ( _extent.x1 - _extent.x0 ) ,
             1.0 / ( _extent.y1 - _extent.y0 ) } ,
    inner ( _inner )
  { }

  // shade provides a pixel value for a 2D coordinate inside the
  // 2D manifold's 'extent' by delegating to the 'inner' functor
  // of class source_t

  px_v shade ( crd_v crd ) const
  {
    // move to texture coordinates. The eval code uses clamping,
    // so absolute precision isn't required and multiplying with
    // the reciprocal girth should be sufficient.

    crd[0] = ( crd[0] - extent.x0 ) * rgirth[0] ;
    crd[1] = ( crd[1] - extent.y0 ) * rgirth[1] ;

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
    if constexpr ( P == RECTILINEAR )
    {
      crd[0] = crd3[0] / crd3[2] ;
      crd[1] = crd3[1] / crd3[2] ;
    }
    else if constexpr ( P == SPHERICAL )
    {
      static ray_to_ll_t ray_to_ll ;
      ray_to_ll.eval ( crd3 , crd ) ;
    }
    else if constexpr ( P == CYLINDRICAL )
    {
      static ray_to_cyl_t ray_to_cyl ;
      ray_to_cyl.eval ( crd3 , crd ) ;
    }
    else if constexpr ( P == STEREOGRAPHIC )
    {
      static ray_to_ster_t ray_to_ster ;
      ray_to_ster.eval ( crd3 , crd ) ;
    }
    else if constexpr ( P == FISHEYE )
    {
      static ray_to_fish_t ray_to_fish ;
      ray_to_fish.eval ( crd3 , crd ) ;
    }
  }

  // get_coordinate yields a coordinate into the mounted 2D manifold
  // and a true mask value if the ray passes through the draped 2D
  // manifold, or a valid coordinate (the center of the 'extent')
  // and a false mask value if the ray does not 'hit' the 2D manifold.

  mask_t get_coordinate ( const crd3_v & crd3 , crd_v & crd ) const
  {
    get_coordinate_nomask ( crd3 , crd ) ;

    auto mask =    ( crd[0] < extent.x0 )
                || ( crd[0] > extent.x1 )
                || ( crd[1] < extent.y0 )
                || ( crd[1] > extent.y1 ) ;

    if constexpr ( P == RECTILINEAR )
      mask |= ( crd3[2] <= 0.0f ) ;

    crd[0] ( mask ) = center[0] ;
    crd[1] ( mask ) = center[1] ;
    return mask ;
  }

  // eval puts it all together and yields pixel sets with 'misses' masked
  // out to zero.

  void eval ( const in_v & ray , px_v & px ) const
  {
    // the first three components have the 3D pickup coordinate itself

    crd3_v crd3 { ray[0] , ray[1] , ray[2] } ;
    crd_v crd ;

    auto mask = get_coordinate ( crd3 , crd ) ;

    if constexpr ( ncrd == 3 )
    {
      // thsi is easy: just call shade

      px = shade ( crd ) ;
    }
    else
    {
      // extract the 3D ray coordinates for the two adjacent locations

      crd3_v dx3 { ray[3] , ray[4] , ray[5] } ;
      crd3_v dy3 { ray[6] , ray[7] , ray[8] } ;

      // obtain their corresponding 2D source image coordinates

      crd_v dx2 , dy2 ;
      get_coordinate_nomask ( dx3 , dx2 ) ;
      get_coordinate_nomask ( dy3 , dy2 ) ;

      // scale to texture coordinates

      crd[0] = ( crd[0] - extent.x0 ) * rgirth[0] ;
      dx2[0] = ( dx2[0] - extent.x0 ) * rgirth[0] ;
      dy2[0] = ( dy2[0] - extent.x0 ) * rgirth[0] ;
      crd[1] = ( crd[1] - extent.y0 ) * rgirth[1] ;
      dx2[1] = ( dx2[1] - extent.y0 ) * rgirth[1] ;
      dy2[1] = ( dy2[1] - extent.y0 ) * rgirth[1] ;

      // delegate to 'inner' to glean the pixel value

      inner.eval ( { crd[0] , crd[1] ,
                     dx2[0] , dx2[1] ,
                     dy2[0] , dy2[1] } , px ) ; 
    }
    // mask out 'misses' to all-zero

    if ( any_of ( mask ) )
    {
      for ( int i = 0 ; i < nchannels ; i++ )
        px[i] ( mask ) = 0.0f ;
    }
  }
} ;

// the 'environment' template codes objects which can serve as 'act'
// functor in zimt::process. It's coded as a zimt::unary_functor
// taking 3D 'ray' coordinates and producing pixels with C channels.

template < typename T , typename U , std::size_t C , std::size_t L >
class environment
: public zimt::unary_functor < zimt::xel_t < T , 3 > ,
                               zimt::xel_t < U , C > ,
                               L >
{
  typedef zimt::unary_functor < zimt::xel_t < T , 3 > ,
                               zimt::xel_t < U , C > ,
                               L > base_t ;

  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , C > > > pa ;

public:
  
  using typename base_t::in_v ;
  using typename base_t::out_v ;

  typedef zimt::xel_t < T , 3 > ray_t ;
  typedef zimt::xel_t < U , C > px_t ;
  typedef zimt::xel_t < std::size_t , 2 > shape_type ;
  typedef zimt::xel_t < float , C > in_px_t ;

  // this member variable will hold the type-erased functor encoding
  // the various code paths we handle with this object. The 'eval'
  // member function simply calls this functor.

  zimt::grok_type < ray_t , px_t , L > env ;

  environment()
  {
    // we have two different code paths to deal with. The first one
    // is for 'mounted' images in various projections. They require
    // specific code to deal with the projection and mask out pixels
    // which aren't covered by the source image data.

    if ( args.mount_image != std::string() )
    {
      shape_type shape { args.mount_width , args.mount_height } ;

      if ( verbose )
        std::cout << "processing mounted image " << args.mount_image
                  << " shape " << shape << std::endl ;

      pa = std::make_shared < zimt::array_t < 2 , in_px_t > >
        ( shape ) ;

      bool success = args.mount_inp->read_image
              ( 0 , 0 , 0 , C , TypeDesc::FLOAT , pa->data() ) ;

      assert ( success ) ;

      static source_t < C , 2 , 16 > src ( *pa ) ;

      // for now, we mount images to the center; other types of cropping
      // might be added by providing suitable parameterization.
  
      auto extent = get_extent ( args.mount_prj , args.mount_width ,
                                 args.mount_height , args.mount_hfov  ) ;

      // we fix the projection as a template argument to class mount_t.

      switch ( args.mount_prj )
      {
        case RECTILINEAR:
        {
          mount_t < C , 3 , RECTILINEAR , 16 > mnt ( extent , src ) ;
          env = mnt ;
          break ;
        }
        case SPHERICAL:
        {
          mount_t < C , 3 , SPHERICAL , 16 > mnt ( extent , src ) ;
          env = mnt ;
          break ;
        }
        case CYLINDRICAL:
        {
          mount_t < C , 3 , CYLINDRICAL , 16 > mnt ( extent , src ) ;
          env = mnt ;
          break ;
        }
        case STEREOGRAPHIC:
        {
          mount_t < C , 3 , STEREOGRAPHIC , 16 > mnt ( extent , src ) ;
          env = mnt ;
          break ;
        }
        case FISHEYE:
        {
          mount_t < C , 3 , FISHEYE , 16 > mnt ( extent , src ) ;
          env = mnt ;
          break ;
        }
        default:
        {
          assert ( false ) ;
          break ;
        }
      }
    }
    else
    {
      // this code path is for 'true' environments which provide pixel
      // data for the entire 360X180 degrees. For now, we accepts either
      // a 2:1 full equirect or a 1:6 vertically stacked cubemap.

      auto & inp ( args.inp ) ;
      auto & w ( args.env_width ) ;
      auto & h ( args.env_height ) ;

      shape_type shape { w , h } ;

      if ( w == 2 * h )
      {
        if ( verbose )
          std::cout << "environment has 2:1 aspect ratio, assuming latlon"
                    << std::endl ;

        pa = std::make_shared < zimt::array_t < 2 , in_px_t > >
          ( shape ) ;

        bool success = inp->read_image ( 0 , 0 , 0 , C ,
                                          TypeDesc::FLOAT , pa->data() ) ;
        assert ( success ) ;

        env =   ray_to_ll_t()
              + eval_latlon < C > ( *pa ) ;
      }
      else if ( h == 6 * w )
      {
        if ( verbose )
        {
          std::cout << "environment has 1:6 aspect ratio, assuming cubemap"
                    << std::endl ;
        }
        sixfold_t<C> sf ( args.env_width , args.cbmfov ,
                          args.support_min , args.tile_size ) ;
        if ( args.multiple_input )
          sf.load ( args.cfs.get_filenames() ) ;
        else
          sf.load ( inp ) ;
        inp->close() ;
        sf.mirror_around() ;
        sf.fill_support( ) ;
        sf.prefilter ( args.prefilter_degree ) ;

        env =   cbm_to_px_t < C > ( sf ) ;
      }
      else
      {
        std::cerr << "environment doesn't have 2:1 or 1:6 aspect ratio"
                  << std::endl ;

        assert ( false ) ;
      }
    }
  }

  // eval simply delegates to 'env'

  void eval ( const in_v & in , out_v & out )
  {
    env.eval ( in , out ) ;
  }

} ;

// environment9 objects mediate lookup with derivatives. Their eval
// member function takes 'ninepacks' and provides pixels with
// nchannels channels.

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

  zimt::grok_type < crd9_t , px_t , 16 > act ;

  environment9()
  {
    if ( args.itp == -2 )
    {
      // use twining (--itp -2). first create an 'environment' object.
      // This is an object holding pixel data and using direct
      // bilinear interpolation to get pixel values. But this
      // object won't directly yield output data - the 'twining'
      // will use it to obtain several pixel values in close
      // vicinity and form a weighted sum from them.
  
      // Note how this object will persist (we declare it static), so
      // it's only created right when control flow arrives here for
      // the very first time (parameterization is taken from 'args').
      // This is deliberate: subsequent invocations of 'work' are
      // supposed to use the same environment, so it would be
      // wasteful to set up a new one - it's an expensive asset.

      // static environment < float , float , nchannels , 16 > env ;

      typedef environment < float , float , nchannels , 16 > env_t ;
      env_t * p_env = (env_t*) current_env.has ( args.input ) ;

      if ( ! p_env )
      {
        current_env.clear() ;
        p_env = new env_t() ;
        current_env.reset ( args.input , p_env ) ;
      }

      // set up the twining filter

      std::vector < zimt::xel_t < float , 3 > > spread ;

      if ( args.twf_file != std::string() )
      {
        // if user passes a twf-file, it's used to set up the twining
        // filter, and all the other twining-related arguments apart
        // from twine_width and twine_normalize are ignored.

        read_twf_file ( spread ) ;
        assert ( spread.size() ) ;
      }
      else
      {
        if ( args.twine == 0 )
        {
          // user has passed twine 0, or not passed anything - but itp
          // is -2. We set up the twining with automatically generated
          // parameters.

          // figure out the magnification in the image center as a
          // guideline

          double mag = args.env_step / args.step ;

          if ( mag > 1.0 )
          {
            // if the transformation magnifies, we use a moderate twine
            // size and a twine_width equal to the magnification, to
            // avoid the star-shaped artifacts from the bilinear
            // interpolation used for the lookup of the contributing
            // rays. If mag is small, the star-shaped artifacts aren't
            // really an issue, and twine values beyond, say, five
            // do little to improve the filter response, so we cap
            // the twine value at five, but lower it when approaching
            // mag 1, down to two - where the next case down starts.

            args.twine = std::min ( 5 , int ( 1.0 + mag ) ) ;
            args.twine_width = mag ;
          }
          else
          {
            // otherwise, we use a twine size which depends on the
            // downscaling factor (reciprocal of 'mag') and a twine_width
            // of 1.0: we only want anti-aliasing. picking a sufficiently
            // large twine value guarantees that we're not skipping any
            // pixels in the source (due to undersampling).

            args.twine = int ( 1.0 + 1.0 / mag ) ;
            args.twine_width = 1.0 ;
          }

          if ( args.twine_density != 1.0f )
          {
            // if the user has passed twine_density, we use it as a
            // multiplicative factor to change args.twine - typically
            // twine_density will be larger than one, so we'll get
            // more filter taps.

            double twine = args.twine * args.twine_density ;
            args.twine = std::round ( twine ) ;
          }

          if ( verbose )
          {
            std::cout << "automatic twining for magnification " << mag
                      << ":" << std::endl ;
            std::cout << "twine: " << args.twine
                      << " twine_width: " << args.twine_width
                      << std::endl ;
          }
        }
        else
        {
          if ( verbose )
          {
            std::cout << "using fixed twine: " << args.twine
                      << " twine_width: " << args.twine_width
                      << std::endl ;
          }
        }

        // with the given twining parameters, we can now set up a 'spread':
        // the generalized equivalent of a filter kernel. While a 'standard'
        // convolution kernel has pre-determined geometry (it's a matrix of
        // coefficients meant to be applied to an equally-shaped matrix of
        // data) - here we have three values for each coefficient: the
        // first two define the position of the look-up relative to the
        // 'central' position, and the third is the weight and corresponds
        // to a 'normal' convolution coefficient.

          make_spread ( spread , args.twine , args.twine ,
                        args.twine_width , args.twine_sigma ,
                        args.twine_threshold ) ;
      }

      // wrap the 'environment' object in a twine_t object and assign
      // to 'act' - this 'groks' the twine_t object to act's type.
      // Note how this object is not declared static: the twining
      // parameters may change due to changing hfov from one invocation
      // of 'work' to the next.

      if ( args.twine_precise )
      {
        act = twine_t < nchannels , 16 , true > ( *p_env , spread ) ;
      }
      else
      {
        act = twine_t < nchannels , 16 , false > ( *p_env , spread ) ;
      }
    }
    else
    {
      // use OIIO's 'environment' or 'texture' lookup functions.
      // These code paths are coded in the 'latlon' and 'cubemap'
      // objects, which are created and also grokked to 'act'.
      // There, the 'ninepack' is used to calculate the derivatives
      // and then all the data needed for the look-up are passed to
      // the relevant (batched) OIIO look-up function. mounted
      // images use an object of type 'mount_t' which, in turn,
      // uses an 'inner functor' of type source_t. So we first
      // set up the source_t object, then create the projection-
      // -specific mount_t object.

      if ( args.mount_image != std::string() )
      {
        typedef source_t < nchannels , 6 , 16 > env_t ;
        env_t * p_env = (env_t*) current_env.has ( args.mount_image ) ;

        if ( ! p_env )
        {
          current_env.clear() ;
          p_env = new env_t() ;
          current_env.reset ( args.mount_image , p_env ) ;
        }
        // static source_t < nchannels , 6 , 16 > src ;
        auto extent = get_extent ( args.mount_prj , args.mount_width ,
                                   args.mount_height , args.mount_hfov  ) ;

        // same case switch as for the single-coordinate code path;
        // might be factored out - for now we just copy and paste

        switch ( args.mount_prj )
        {
          case RECTILINEAR:
          {
            mount_t < nchannels , 9 , RECTILINEAR , 16 >
              mnt ( extent , *p_env ) ;
            act = mnt ;
            break ;
          }
          case SPHERICAL:
          {
            mount_t < nchannels , 9 , SPHERICAL , 16 >
              mnt ( extent , *p_env ) ;
            act = mnt ;
            break ;
          }
          case CYLINDRICAL:
          {
            mount_t < nchannels , 9 , CYLINDRICAL , 16 >
              mnt ( extent , *p_env ) ;
            act = mnt ;
            break ;
          }
          case STEREOGRAPHIC:
          {
            mount_t < nchannels , 9 , STEREOGRAPHIC , 16 >
              mnt ( extent , *p_env ) ;
            act = mnt ;
            break ;
          }
          case FISHEYE:
          {
            mount_t < nchannels , 9 , FISHEYE , 16 >
              mnt ( extent , *p_env ) ;
            act = mnt ;
            break ;
          }
          default:
          {
            assert ( false ) ;
            break ;
          }
        }
      }
      else if ( args.env_width == args.env_height * 2 )
      {
        // set up an environment object picking up pixel values from
        // a lat/lon image using OIIO's 'environment' function. Again
        // we set up a static object, but it's assigned to 'act' in
        // every invocation of 'work'.

        // static latlon < float , float , nchannels , 16 > ll ;
        typedef latlon < float , float , nchannels , 16 > env_t ;
        env_t * p_env = (env_t*) current_env.has ( args.input ) ;

        if ( ! p_env )
        {
          current_env.clear() ;
          p_env = new env_t() ;
          current_env.reset ( args.input , p_env ) ;
        }
        act = *p_env ;
      }
      else
      {
        // set up an environment object picking up pixel values from
        // a texture representing the cubemap image, using OIIO's
        // 'texture' function

        // static cubemap < float , float , nchannels , 16 > cbm ;
        typedef cubemap < float , float , nchannels , 16 > env_t ;
        env_t * p_env = (env_t*) current_env.has ( args.input ) ;

        if ( ! p_env )
        {
          current_env.clear() ;
          p_env = new env_t() ;
          current_env.reset ( args.input , p_env ) ;
        }
        act = *p_env ;
      }
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
