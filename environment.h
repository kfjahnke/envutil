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

#include <OpenImageIO/imageio.h>
#include <OpenImageIO/texture.h>

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

  const double scale ;

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

    const zimt::xel_t < double , 2 > shift { M_PI , M_PI_2 } ;

    auto in = ( _in + shift ) * scale ;

    // if the coordinate is now precisely zero, this corresponds
    // to a pick-up point at the top or left margin of the UL
    // pixel, which is 0.5 away from it's center

    in -= .5 ;

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

    auto wl = 1.0 - wr ;

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

    const auto * p = (float*) ( latlon.data() ) ;

    // obtain the four constituents by first truncating their
    // coordinate to int and then gathering from p.

    index_v idsdxl { uli[0] , uli[1] } ;
    auto ofs = ( idsdxl * latlon.strides ) . sum() * nchannels ;
    px_v pxul ;
    pxul.gather ( p , ofs ) ;

    index_v idsdxr { uri[0] , uri[1] } ;
    ofs = ( idsdxr * latlon.strides ) . sum() * nchannels ;
    px_v pxur ;
    pxur.gather ( p , ofs ) ;

    index_v idxll { lli[0] , lli[1] } ;
    ofs = ( idxll * latlon.strides ) . sum() * nchannels ;
    px_v pxll ;
    pxll.gather ( p , ofs ) ;

    index_v idxlr { lri[0] , lri[1] } ;
    ofs = ( idxlr * latlon.strides ) . sum() * nchannels ;
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
  { "MipModeStochasticTrilinear" , OIIO::TextureOpt::MipModeStochasticTrilinear } ,
  { "MipModeStochasticAniso" , OIIO::TextureOpt::MipModeStochasticAniso }
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
  OIIO::TextureSystem * ts ;
  OIIO::TextureOptBatch batch_options ;
  OIIO::TextureSystem::TextureHandle * th ;

  // pull in the c'tor arguments

  eval_env ( const arguments & args )
  : ts ( OIIO::TextureSystem::create() )
  {
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
  
    OIIO::ustring uenvironment ( args.input.c_str() ) ;
    th = ts->get_texture_handle ( uenvironment ) ;
  }

  // set up the eval function.

  template < typename I , typename O >
  void eval ( const I & crd9 , O & px )
  {
    // we convert to OIIO-compatible coordinates as we go, this
    // requires changing the sign of the y and z axis (hance the
    // seemingly wrong order of the differences for ds and dt)

    crd3_v c3 {   crd9[0] ,
                - crd9[1] ,
                - crd9[2] } ;

    crd3_v ds { crd9[3] - crd9[0] ,
                crd9[1] - crd9[4] ,
                crd9[2] - crd9[5] } ;

    crd3_v dt { crd9[6] - crd9[0] ,
                crd9[1] - crd9[7] ,
                crd9[2] - crd9[8] } ;

    // now we can call 'environment', but depending on the SIMD
    // back-end, we provide the pointers which 'environemnt' needs
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
// optimized away. This object is not well-tested.

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

#include "metrics.h"
#include <filesystem>

// stripped-down version of sixfold_t, omitting the use of OIIO's
// 'texture' function and twining, using bilinear interpolation
// only. The focus here is on testing the steppers in stepper.h.

template < std::size_t nchannels >
struct sixfold_t
: public metrics_t
{
  // shorthand for pixels and SIMDized pixels

  typedef zimt::xel_t < float , nchannels > px_t ;
  typedef zimt::simdized_type < px_t , LANES > px_v ;

  // pointers to an OIIO texture system and an OIIO texture handle.
  // These are only used if the pick-up is done using OIIO's
  // 'texture' function.

  std::filesystem::path texture_file ;
  OIIO::TextureSystem * ts ;
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

  sixfold_t ( std::size_t _face_px ,
              std::size_t _support_min_px = 4UL ,
              std::size_t _tile_px = 64UL ,
              double _face_fov = M_PI_2 )
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
        buffer ( { face_px , face_px , 6UL } ) ;
      
      // and a view to the store with the same shape, but strides to
      // match the metrics of the target memory area in the store
        
      zimt::view_t < 3 , px_t >
        target ( p_ul ,
                { 1L , long(section_px) , long(offset_px) } ,
                { face_px , face_px , 6UL } ) ;

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
                << texture_file.c_str() << std::endl ;

    save_array ( texture_file.c_str() , store ) ;

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

    OIIO::ustring uenvironment ( texture_file.c_str() ) ;
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

  void cubemap_to_pixel ( const i_v & face ,
                          const crd2_v & in_face ,
                          px_v & px ) const
  {
    // get the pick-up coordinate, in pixel units

    crd2_v pickup ;
    get_pickup_coordinate_px ( face , in_face , pickup ) ;

    // bilinear interpolation. Since we have incoming coordinates
    // in the range of (-0.5, width-0,5) relative to the 'ninety
    // degrees proper', this interpolator may 'look at' pixels
    // outside the 'ninety degrees proper'. So we need sufficient
    // support here: pickup coordinates with negative values will,
    // for example, look at pixels just outside the 'ninety degree
    // zone'. Note that we want correct support, rather than using
    // the next-best thing (like mirroring on the edge), and we'll
    // generate this support. Then we can be sure that the output
    // is optimal and doesn't carry in any artifacts from the edges.
    // So, to be quite clear: this bit of code does a bilinear
    // interpolation without bounds checking, assuming that there
    // is sufficient support for every expected pick-up coordinate.
    // Fuzzing this with arbitrary coordinates would fail.

    const auto * p = (float*) ( store.data() ) ;
    px_v px2, help ;

    // find the floor of the pickup coordinate by truncation

    v2i_v low { pickup[0] , pickup[1] } ;

    // how far is the pick-up coordinate from the floor value?

    auto diff = pickup - crd2_v ( low ) ;

    // gather data pertaining to the pixels at the 'low' coordinate

    v2i_v idx { low[0] , low[1] } ;
    auto ofs = ( idx * store.strides ) . sum() * nchannels ;
    px.gather ( p , ofs ) ;

    // and weight them according to distance from the 'low' value

    auto one = f_v::One() ;

    px *= ( one - diff[0] ) ;

    // repeat the process for the low coordinate's neighbours

    idx[0] += 1 ; ;
    ofs = ( idx * store.strides ) . sum() * nchannels ;
    help.gather ( p , ofs ) ;
    px += help * diff[0] ;

    // the first partial sum is also weighted, now according to
    // vertical distance

    px *= ( one - diff[1] ) ;

    idx[0] -= 1 ; ;
    idx[1] += 1 ; ;
    ofs = ( idx * store.strides ) . sum() * nchannels ;
    px2.gather ( p , ofs ) ;
    px2 *= ( one - diff[0] ) ;

    idx[0] += 1 ; ;
    ofs = ( idx * store.strides ) . sum() * nchannels ;
    help.gather ( p , ofs ) ;
    px2 += help * diff[0] ;

    px += px2 * diff[1] ;
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
    const int degree ;
    const int ithird ;

    typedef zimt::xel_t < float , nchannels > px_t ;
    typedef zimt::simdized_type < px_t , LANES > px_v ;

    // note the factor of two in the initialization of 'ithird':
    // incoming coordinates are doubled (!) for the purpose at hand.

    fill_frame_t ( const sixfold_t < nchannels > & _cubemap ,
                  const int & _face ,
                  const int & _degree )
    : sf ( _cubemap ) ,
      face ( _face ) ,
      degree ( _degree ) ,
      ithird ( _cubemap.model_to_px * 2 )
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
      ray_to_cubeface ( crd3 , fv , in_face ) ;

      // finally we use this information to obtain pixel values,
      // which are written to the target location.

      sf.cubemap_to_pixel ( fv , in_face , px ) ;
    }
  } ;

  // fill_support uses the fill_frame_t functor to populate the
  // frame of support. The structure of the code is similar to
  // 'mirror_around', iterating over the six sections of the
  // array and manipulating each in turn. But here we fill in
  // the entire surrounding frame, not just a pixel-wide line,
  // and we pick up data from neighbouring cube faces.

  void fill_support ( int degree )

  {
    if ( left_frame_px == 0 && right_frame_px == 0 )
      return ;

    auto * p_base = store.data() ;

    for ( int face = 0 ; face < 6 ; face++ )
    {
      // set up the 'gleaning' functor

      fill_frame_t fill_frame ( *this , face , degree ) ;
    
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
    ray_to_cubeface ( crd3 , face , in_face ) ;

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
    ray_to_cubeface ( c3 , face , in_face_00 ) ;

    // do the same for the two offsetted rays
  
    ray_to_cubeface_fixed ( dx , face , in_face_10 ) ;
    ray_to_cubeface_fixed ( dy , face , in_face_01 ) ;

    // obtain the pickup coordinate in texture units (in [0,1])

    crd2_v c3_tx ;
    cubemap.get_pickup_coordinate_tx ( face , in_face_00 , c3_tx ) ;

    // convert the neighbouring rays to texture coordinates
    
    // in_face_00 += 1.0f ;
    in_face_00 *= scale ;
    // in_face_10 += 1.0f ;
    in_face_10 *= scale ;
    // in_face_01 += 1.0f ;
    in_face_01 *= scale ;

    // form the differences

    crd2_v dx_tx = in_face_10 - in_face_00 ;
    crd2_v dy_tx = in_face_01 - in_face_00 ;

    // cubemap.get_pickup_coordinate_tx ( face , in_face_10 , dx_tx ) ;
    // cubemap.get_pickup_coordinate_tx ( face , in_face_01 , dy_tx ) ;

//     // form the differences to get the derivatives
//   
//     dx_tx -= c3_tx ;
//     dy_tx -= c3_tx ;

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
// taking 3D 'ray' coordinates and producing pixels with C channels.
// struct repix is used to convert the output to the desired number
// of channels.

template < typename T , typename U , std::size_t C , std::size_t L >
class latlon
: public zimt::unary_functor < zimt::xel_t < T , 9 > ,
                               zimt::xel_t < U , C > ,
                               L >
{
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , 1 > > > pa1 ;
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , 2 > > > pa2 ;
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , 3 > > > pa3 ;
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , 4 > > > pa4 ;

public:
  
  typedef zimt::xel_t < T , 9 > ray_t ;
  typedef zimt::xel_t < U , C > px_t ;
  typedef zimt::xel_t < std::size_t , 2 > shape_type ;

  zimt::grok_type < ray_t , px_t , L > env ;

  latlon ( const std::unique_ptr<ImageInput> & inp ,
           const std::size_t & w ,
           const std::size_t & h ,
           const std::size_t & nchannels ,
           const arguments & args )
  {
    shape_type shape { w , h } ;

    assert ( w == 2 * h ) ;

    if ( verbose )
      std::cout << "input has 2:1 aspect ratio, assuming latlon"
                << std::endl ;

    if ( nchannels >= 4 )
    {
      typedef zimt::xel_t < float , 4 > in_px_t ;
      pa4 = std::make_shared < zimt::array_t < 2 , in_px_t > >
        ( shape ) ;

      bool success = inp->read_image ( 0 , 0 , 0 , nchannels ,
                                      TypeDesc::FLOAT , pa4->data() ) ;
      assert ( success ) ;

      env = zimt::grok_type < ray_t , px_t , L >
        ( eval_env < C > ( args ) ) ;
    }
    else if ( nchannels == 3 )
    {
      typedef zimt::xel_t < float , 3 > in_px_t ;
      pa3 = std::make_shared < zimt::array_t < 2 , in_px_t > >
        ( shape ) ;

      bool success = inp->read_image ( 0 , 0 , 0 , nchannels ,
                                      TypeDesc::FLOAT , pa3->data() ) ;
      assert ( success ) ;

      env = zimt::grok_type < ray_t , px_t , L >
        ( eval_env < C > ( args ) ) ;
    }
    else if ( nchannels == 2 )
    {
      typedef zimt::xel_t < float , 2 > in_px_t ;
      pa2 = std::make_shared < zimt::array_t < 2 , in_px_t > >
        ( shape ) ;
    
      bool success = inp->read_image ( 0 , 0 , 0 , nchannels ,
                                      TypeDesc::FLOAT , pa2->data() ) ;
      assert ( success ) ;
    
      env = zimt::grok_type < ray_t , px_t , L >
        ( eval_env < C > ( args ) ) ;
    }
    else if ( nchannels == 1 )
    {
      typedef zimt::xel_t < float , 1 > in_px_t ;
      pa1 = std::make_shared < zimt::array_t < 2 , in_px_t > >
        ( shape ) ;
    
      bool success = inp->read_image ( 0 , 0 , 0 , nchannels ,
                                      TypeDesc::FLOAT , pa2->data() ) ;
      assert ( success ) ;
    
      env = zimt::grok_type < ray_t , px_t , L >
        ( eval_env < C > ( args ) ) ;
    }
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
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , 1 > > > pa1 ;
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , 2 > > > pa2 ;
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , 3 > > > pa3 ;
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , 4 > > > pa4 ;

public:
  
  typedef zimt::xel_t < T , 9 > ray_t ;
  typedef zimt::xel_t < U , C > px_t ;
  typedef zimt::xel_t < std::size_t , 2 > shape_type ;

  zimt::grok_type < ray_t , px_t , L > env ;

  cubemap ( const std::unique_ptr<ImageInput> & inp ,
            const std::size_t & w ,
            const std::size_t & h ,
            const std::size_t & nchannels ,
            const arguments & args )
  {
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
    if ( nchannels >= 4 )
    {
      sixfold_t<4> sf ( w ) ;
      sf.load ( inp ) ;
      inp->close() ;
      sf.mirror_around() ;
      sf.fill_support ( 1 ) ;
      sf.gen_texture ( args.tsoptions ) ;
      std::cout << "cbm set env 4" << std::endl ;

      env =   cbm_to_px_t2 < 4 > ( sf , batch_options )
            + repix < U , 4 , C , L > () ;
    }
    else if ( nchannels == 3 )
    {
      sixfold_t<3> sf ( w ) ;
      sf.load ( inp ) ;
      inp->close() ;
      sf.mirror_around() ;
      sf.fill_support ( 1 ) ;
      sf.gen_texture ( args.tsoptions ) ;
      std::cout << "cbm set env 3" << std::endl ;

      env =   cbm_to_px_t2 < 3 > ( sf , batch_options )
            + repix < U , 3 , C , L > () ;
    }
    else if ( nchannels == 2 )
    {
      sixfold_t<2> sf ( w ) ;
      sf.load ( inp ) ;
      inp->close() ;
      sf.mirror_around() ;
      sf.fill_support ( 1 ) ;
      sf.gen_texture ( args.tsoptions ) ;
      std::cout << "cbm set env 2" << std::endl ;

      env =   cbm_to_px_t2 < 2 > ( sf , batch_options )
            + repix < U , 2 , C , L > () ;
    }
    else
    {
      sixfold_t<1> sf ( w ) ;
      sf.load ( inp ) ;
      inp->close() ;
      sf.mirror_around() ;
      sf.fill_support ( 1 ) ;
      sf.gen_texture ( args.tsoptions ) ;
      std::cout << "cbm set env 1" << std::endl ;

      env =   cbm_to_px_t2 < 1 > ( sf , batch_options )
            + repix < U , 1 , C , L > () ;
    }
  }

  void eval ( const typename zimt::grok_type < ray_t , px_t , L >::in_v & in ,
              typename zimt::grok_type < ray_t , px_t , L >::out_v & out )
  {
    env.eval ( in , out ) ;
  }

} ;

// the 'environment' template codes objects which can serve as 'act'
// functor in zimt::process. It's coded as a zimt::unary_functor
// taking 3D 'ray' coordinates and producing pixels with C channels.
// struct repix is used to convert the output to the desired number
// of channels.

template < typename T , typename U , std::size_t C , std::size_t L >
class environment
: public zimt::unary_functor < zimt::xel_t < T , 3 > ,
                               zimt::xel_t < U , C > ,
                               L >
{
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , 1 > > > pa1 ;
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , 2 > > > pa2 ;
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , 3 > > > pa3 ;
  std::shared_ptr < zimt::array_t < 2 , zimt::xel_t < T , 4 > > > pa4 ;

public:
  
  typedef zimt::xel_t < T , 3 > ray_t ;
  typedef zimt::xel_t < U , C > px_t ;
  typedef zimt::xel_t < std::size_t , 2 > shape_type ;

  zimt::grok_type < ray_t , px_t , L > env ;

  environment ( const std::unique_ptr<ImageInput> & inp ,
                const std::size_t & w ,
                const std::size_t & h ,
                const std::size_t & nchannels ,
                const arguments & args )
  {
    assert ( inp ) ;
    shape_type shape { w , h } ;

    if ( w == 2 * h )
    {
      if ( verbose )
        std::cout << "input has 2:1 aspect ratio, assuming latlon"
                  << std::endl ;

      if ( nchannels >= 4 )
      {
        typedef zimt::xel_t < float , 4 > in_px_t ;
        pa4 = std::make_shared < zimt::array_t < 2 , in_px_t > >
          ( shape ) ;

        bool success = inp->read_image ( 0 , 0 , 0 , nchannels ,
                                         TypeDesc::FLOAT , pa4->data() ) ;
        assert ( success ) ;

        env =   ray_to_ll_t()
              + eval_latlon < 4 > ( *pa4 )
              + repix < U , 4 , C , L > () ;
      }
      else if ( nchannels == 3 )
      {
        typedef zimt::xel_t < float , 3 > in_px_t ;
        pa3 = std::make_shared < zimt::array_t < 2 , in_px_t > >
          ( shape ) ;

        bool success = inp->read_image ( 0 , 0 , 0 , nchannels ,
                                         TypeDesc::FLOAT , pa3->data() ) ;
        assert ( success ) ;

        env =   ray_to_ll_t()
              + eval_latlon < 3 > ( *pa3 )
              + repix < U , 3 , C , L > () ;
      }
      else if ( nchannels == 2 )
      {
        typedef zimt::xel_t < float , 2 > in_px_t ;
        pa2 = std::make_shared < zimt::array_t < 2 , in_px_t > >
          ( shape ) ;
     
        bool success = inp->read_image ( 0 , 0 , 0 , nchannels ,
                                         TypeDesc::FLOAT , pa2->data() ) ;
        assert ( success ) ;
      
        env =   ray_to_ll_t()
              + eval_latlon < 2 > ( *pa2 )
              + repix < U , 2 , C , L > () ;
      }
      else
      {
        typedef zimt::xel_t < float , 1 > in_px_t ;
        pa1 = std::make_shared < zimt::array_t < 2 , in_px_t > >
          ( shape ) ;
     
        bool success = inp->read_image ( 0 , 0 , 0 , nchannels ,
                                         TypeDesc::FLOAT , pa1->data() ) ;
        assert ( success ) ;
      
        env =   ray_to_ll_t()
              + eval_latlon < 1 > ( *pa1 )
              + repix < U , 1 , C , L > () ;
      }
    }
    else if ( h == 6 * w )
    {
      if ( verbose )
      {
        std::cout << "input has 1:6 aspect ratio, assuming cubemap"
                  << std::endl ;
      }
      if ( nchannels >= 4 )
      {
        sixfold_t<4> sf ( w ) ;
        sf.load ( inp ) ;
        inp->close() ;
        sf.mirror_around() ;
        sf.fill_support ( 1 ) ;
        std::cout << "cbm set env 4" << std::endl ;

        env =   cbm_to_px_t < 4 > ( sf )
              + repix < U , 4 , C , L > () ;
      }
      else if ( nchannels == 3 )
      {
        sixfold_t<3> sf ( w ) ;
        sf.load ( inp ) ;
        inp->close() ;
        sf.mirror_around() ;
        sf.fill_support ( 1 ) ;
        std::cout << "cbm set env 3" << std::endl ;

        env =   cbm_to_px_t < 3 > ( sf )
              + repix < U , 3 , C , L > () ;
      }
      else if ( nchannels == 2 )
      {
        sixfold_t<2> sf ( w ) ;
        sf.load ( inp ) ;
        inp->close() ;
        sf.mirror_around() ;
        sf.fill_support ( 1 ) ;
        std::cout << "cbm set env 2" << std::endl ;

        env =   cbm_to_px_t < 2 > ( sf )
              + repix < U , 2 , C , L > () ;
      }
      else
      {
        sixfold_t<1> sf ( w ) ;
        sf.load ( inp ) ;
        inp->close() ;
        sf.mirror_around() ;
        sf.fill_support ( 1 ) ;
        std::cout << "cbm set env 1" << std::endl ;

        env =   cbm_to_px_t < 1 > ( sf )
              + repix < U , 1 , C , L > () ;
      }
    }
  }

  void eval ( const typename zimt::grok_type < ray_t , px_t , L >::in_v & in ,
              typename zimt::grok_type < ray_t , px_t , L >::out_v & out )
  {
    env.eval ( in , out ) ;
  }

} ;
