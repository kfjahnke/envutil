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

// sixfold_t provides an internal representation of a cubemap
// with widened support, to provide for easy mip-mapping and
// interpolation. It's also a zimt::unary_functor providing
// pixel values for 3D ray coordinates. Since a cubemap can
// provide pixel values for all directions, there is no code
// to provide masks.

template < std::size_t nchannels , projection_t cbm_prj >
struct sixfold_t
: public metrics_t ,
  public zimt::unary_functor
   < zimt::xel_t < float , 3 > ,
     zimt::xel_t < float , nchannels > ,
     LANES >
{
  // shorthand for pixels and SIMDized pixels

  typedef zimt::xel_t < float , 2 > crd2_t ;
  typedef zimt::xel_t < float , nchannels > px_t ;
  typedef zimt::simdized_type < px_t , LANES > px_v ;

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

    if ( args.verbose )
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

private:

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
  // the support area. For clarity: the next stage will not pick up
  // pixel values from the frame of mirrored pixels - this frame only
  // serves to avoid having black pixels in the support needed for
  // the actual pick-ups, which will all happen at coordinates inside
  // the cube faces proper. This function is only used internally by
  // 'fill_support'.

  void mirror_around()
  {
    auto * p_base = store.data() ;

    for ( int face = 0 ; face < 6 ; face++ )
    {
      // get a pointer to the upper left of the cube face 'proper'

      auto * p_frame =   p_base 
                       + face * section_px * store.strides[1]
                       + ( left_frame_px * store.strides ) .sum() ;

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
    const sixfold_t < nchannels , cbm_prj > & sf ;
    const int face ;
    const int ithird ;

    typedef zimt::xel_t < float , nchannels > px_t ;
    typedef zimt::simdized_type < px_t , LANES > px_v ;

    // we use a degree-1 b-spline (bilinear interpolation) to
    // fill the support frame
    // TODO: when specializing with 1, invocation with
    // --spline_degree 0 crashes.

    zimt::evaluator < crd2_t , px_t , LANES , 1 > ev ;

    // note the factor of two in the initialization of 'ithird':
    // incoming coordinates are doubled (!) for the purpose at hand.

    fill_frame_t ( const sixfold_t < nchannels , cbm_prj > & _cubemap ,
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

      // initially I coded to use the in-plane transformation here:

      // if constexpr ( cbm_prj == BIATAN6 )
      //   in_face = float ( 4.0 / M_PI ) * atan ( in_face ) ;

      // but that's an error: the incoming image is already in
      // the given in-plane representation, and we're only copying
      // content from a different cube face.

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
    mirror_around() ;

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
        zimt::process ( shp , ls , fill_frame , st , bill ) ;
      }
    }
  }

  // if we're using b-spline interpolation for the cubemap,
  // we need to prefilter the sections appropriately. To avoid
  // the ringing artifacts we'd produce when running the prefilter
  // over the entire IR image, we prefilter the six sections
  // separately with 'NATURAL' boundary conditions. If the spline
  // degree given for prefiltering is less than two, the prefilter
  // is not needed.

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

public:

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
      if ( args.verbose )
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
      if ( args.verbose )
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
    fill_support() ;
    prefilter ( args.prefilter_degree ) ;
  }

  void load ( const std::string & filename )
  {
    auto inp = ImageInput::open ( filename ) ;
    load ( inp ) ;
    inp->close() ;
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

      if ( args.verbose )
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
    fill_support() ;
    prefilter ( args.prefilter_degree ) ;
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

  void eval ( const crd3_v & crd3 , px_v & px )
  {
    // find the cube face and in-face coordinate for 'c'

    i_v face ;
    crd2_v in_face ;
    ray_to_cubeface < float , LANES > ( crd3 , face , in_face ) ;

    if constexpr ( cbm_prj == BIATAN6 )
      in_face = float ( 4.0 / M_PI ) * atan ( in_face ) ;

    cubemap_to_pixel ( face , in_face , px ) ;
  }

} ; // end of struct sixfold_t

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

  typedef zimt::bspline < px_t , 2 > spl_t ;
  std::shared_ptr < spl_t > p_bspl ;
  zimt::grok_type < crd_t , px_t , L > bsp_ev ;

  // the c'tor with arguments receives a view to image data and extracts
  // the relevant values to access these data directly. This is the route
  // taken for direct bilinear interpolation of the source data. Note how
  // we cast data() to float to 'shed' the channel count from the type

  source_t ( std::shared_ptr < spl_t > _p_bspl )
  : p_bspl ( _p_bspl ) ,
    width ( _p_bspl->core.shape[0] ) ,
    height ( _p_bspl->core.shape[1] )
  {
    bsp_ev = zimt::make_safe_evaluator
                     < spl_t , float , L > ( *p_bspl ) ;
  }

  void eval ( const pkg_v & _crd , px_v & px )
  {
    // incoming is const &, hence:

    pkg_v crd ( _crd ) ;

    if constexpr ( ncrd == 2 )
    {
      crd[0] *= width ;
      crd[1] *= height ;
      crd -= .5f ;
      bsp_ev.eval ( crd , px ) ;
    }
    else
    {
      assert ( false ) ;
    }
  }
} ;

// struct mount_t provides data from a rectanglular 2D manifold holding
// pixel data in a given projection which do not cover the entire
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
  source_t < nchannels , 2 * ncrd / 3 , L > inner ;

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
    rgirth { 1.0 / ( _extent.x1 - _extent.x0 ) ,
             1.0 / ( _extent.y1 - _extent.y0 ) } ,
    inner ( _inner )
  { }

  // shade provides a pixel value for a 2D coordinate inside the
  // 2D manifold's 'extent' by delegating to the 'inner' functor
  // of class source_t

  px_v shade ( crd_v crd )
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
  }

  // get_coordinate yields a coordinate into the mounted 2D manifold
  // and a true mask value if the ray passes through the draped 2D
  // manifold, or a valid coordinate (the center of the 'extent')
  // and a false mask value if the ray does not 'hit' the 2D manifold.

  mask_t get_coordinate ( const crd3_v & crd3 , crd_v & crd ) const
  {
    get_coordinate_nomask ( crd3 , crd ) ;

    auto mask =    ( crd[0] >= extent.x0 )
                && ( crd[0] <= extent.x1 )
                && ( crd[1] >= extent.y0 )
                && ( crd[1] <= extent.y1 ) ;

    if constexpr ( P == RECTILINEAR )
      mask &= ( crd3[2] > 0.0f ) ;

    crd[0] ( ! mask ) = center[0] ;
    crd[1] ( ! mask ) = center[1] ;

    return mask ;
  }

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

  void eval ( const in_v & ray , px_v & px )
  {
    // the first three components have the 3D pickup coordinate itself

    crd3_v crd3 { ray[0] , ray[1] , ray[2] } ;
    crd_v crd ;

    auto mask = get_coordinate ( crd3 , crd ) ;
    if ( none_of ( mask ) )
    {
      px = 0.0f ;
      return ;
    }

    if constexpr ( ncrd == 3 )
    {
      // this is easy: just call shade

      px = shade ( crd ) ;
    }
    else
    {
      assert ( false ) ;
    }
    // mask out 'misses' to all-zero

    if ( ! all_of ( mask ) )
    {
      for ( int i = 0 ; i < nchannels ; i++ )
        px[i] ( ! mask ) = 0.0f ;
    }
  }

  // variant which returns a mask and does not paint misses 0000

  mask_t mask_eval ( const in_v & ray , px_v & px )
  {
    // the first three components have the 3D pickup coordinate itself

    crd3_v crd3 { ray[0] , ray[1] , ray[2] } ;
    crd_v crd ;

    auto mask = get_coordinate ( crd3 , crd ) ;

    if ( any_of ( mask ) )
    {
      if constexpr ( ncrd == 3 )
      {
        // this is easy: just call shade

        px = shade ( crd ) ;
      }
      else
      {
        assert ( false ) ;
      }
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

// the 'environment' template codes objects which can serve as 'act'
// functor in zimt::process. It's coded as a zimt::unary_functor
// taking 3D 'ray' coordinates and producing pixels with C channels.
// 'environment' objects also provide the substrate for 'twining'.
// All types of environment objects are based on cardinal b-splines
// over the image data. The default is to use degree-1 b-splines,
// which is also known as 'bilinear interpolation'. Other spline
// degrees can be specified by passing 'spline_degree'. If necessary,
// the image data are prefilered, so as to provide an interpolating
// spline.

template < typename T , typename U , std::size_t C , std::size_t L >
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
  typedef zimt::xel_t < T , 3 > ray_t ;
  typedef zimt::xel_t < U , C > px_t ;
  typedef zimt::xel_t < std::size_t , 2 > shape_type ;
  typedef zimt::xel_t < float , C > in_px_t ;

  zimt::grok_type < ray_t , px_t , L > env ;

  // we'll hold on to image data via a std::shared_ptr to a zimt::bspline.

  typedef zimt::bspline < in_px_t , 2 > spl_t ;
  std::shared_ptr < spl_t > p_bspl ;

  // this member variable will hold the type-erased functor encoding
  // the various code paths we handle with this object. The 'eval'
  // member function simply calls this functor.

  typedef std::function < mask_t ( const in_v & ) > mask_f ;
  mask_f get_mask ;

  // this c'tor routes to OIIO code

  // environment ( const std::size_t facet_no )
  // {
    // source_t < C , 2 , 16 > src ( facet_no ) ;
    // 
    // // for now, we mount images to the center; other types of cropping
    // // might be added by providing suitable parameterization.
    // 
    // auto const & fct ( args.facet_spec_v [ facet_no ] ) ;
    // auto extent = get_extent ( fct.projection , fct.width ,
    //                            fct.height , fct.hfov  ) ;
    // 
    // // we fix the projection as a template argument to class mount_t.
    // 
    // switch ( fct.projection )
    // {
    //   case RECTILINEAR:
    //   {
    //     mount_t < C , 3 , RECTILINEAR , 16 > mnt ( extent , src ) ;
    //     env = mnt ;
    //     get_mask = [=] ( const in_v & crd3 )
    //       { return mnt.get_mask ( crd3 ) ; } ;
    //     break ;
    //   }
    //   case SPHERICAL:
    //   {
    //     mount_t < C , 3 , SPHERICAL , 16 > mnt ( extent , src ) ;
    //     env = mnt ;
    //     get_mask = [=] ( const in_v & crd3 )
    //       { return mnt.get_mask ( crd3 ) ; } ;
    //     break ;
    //   }
    //   case CYLINDRICAL:
    //   {
    //     mount_t < C , 3 , CYLINDRICAL , 16 > mnt ( extent , src ) ;
    //     env = mnt ;
    //     get_mask = [=] ( const in_v & crd3 )
    //       { return mnt.get_mask ( crd3 ) ; } ;
    //     break ;
    //   }
    //   case STEREOGRAPHIC:
    //   {
    //     mount_t < C , 3 , STEREOGRAPHIC , 16 > mnt ( extent , src ) ;
    //     env = mnt ;
    //     get_mask = [=] ( const in_v & crd3 )
    //       { return mnt.get_mask ( crd3 ) ; } ;
    //     break ;
    //   }
    //   case FISHEYE:
    //   {
    //     mount_t < C , 3 , FISHEYE , 16 > mnt ( extent , src ) ;
    //     env = mnt ;
    //     get_mask = [=] ( const in_v & crd3 )
    //       { return mnt.get_mask ( crd3 ) ; } ;
    //     break ;
    //   }
    //   default:
    //   {
    //     std::cerr << "unknown projection: "
    //               << fct.projection_str << std::endl ;
    //     assert ( false ) ;
    //     break ;
    //   }
    // }
  // }

  environment ( const facet_spec & fct )
  {
    shape_type shape { fct.width , fct.height } ;

    zimt::bc_code bc0 = zimt::REFLECT ;
    if (    fct.projection == SPHERICAL
         || fct.projection == CYLINDRICAL )
    {
      if ( fabs ( fct.hfov - 2.0 * M_PI ) < .000001 )
        bc0 = PERIODIC ;
    }
    p_bspl.reset ( new spl_t ( shape , args.spline_degree ,
                                { bc0 , zimt::REFLECT } ) ) ;

    auto inp = ImageInput::open ( fct.filename ) ;

    bool success = inp->read_image (
      0 , 0 , 0 , C ,
      TypeDesc::FLOAT ,
      p_bspl->core.data() ,
      sizeof ( in_px_t ) ,
      p_bspl->core.strides[1] * sizeof ( in_px_t ) ) ;
    assert ( success ) ;
    inp->close() ;

    if (    fct.projection == SPHERICAL
          && fabs ( fct.hfov - 2.0 * M_PI ) < .000001
          && fct.width == 2 * fct.height )
    {
      // assume the mounted image is a full spherical
      p_bspl->spline_degree = args.prefilter_degree ;
      spherical_prefilter ( *p_bspl , p_bspl->core , zimt::default_njobs ) ;
      p_bspl->spline_degree = args.spline_degree ;
    }
    else
    {
      // assume it's a partial spherical, use ordinary prefilter
      p_bspl->spline_degree = args.prefilter_degree ;
      p_bspl->prefilter() ;
      p_bspl->spline_degree = args.spline_degree ;
    }

    source_t < C , 2 , 16 > src ( p_bspl ) ;

    // for now, we mount images to the center; other types of cropping
    // might be added by providing suitable parameterization.

    auto extent = get_extent ( fct.projection , fct.width ,
                                fct.height , fct.hfov  ) ;

    // we fix the projection as a template argument to class mount_t.

    switch ( fct.projection )
    {
      case RECTILINEAR:
      {
        mount_t < C , 3 , RECTILINEAR , 16 > mnt ( extent , src ) ;
        env = mnt ;
        get_mask = [=] ( const in_v & crd3 )
          { return mnt.get_mask ( crd3 ) ; } ;
        break ;
      }
      case SPHERICAL:
      {
        mount_t < C , 3 , SPHERICAL , 16 > mnt ( extent , src ) ;
        env = mnt ;
        get_mask = [=] ( const in_v & crd3 )
          { return mnt.get_mask ( crd3 ) ; } ;
        break ;
      }
      case CYLINDRICAL:
      {
        mount_t < C , 3 , CYLINDRICAL , 16 > mnt ( extent , src ) ;
        env = mnt ;
        get_mask = [=] ( const in_v & crd3 )
          { return mnt.get_mask ( crd3 ) ; } ;
        break ;
      }
      case STEREOGRAPHIC:
      {
        mount_t < C , 3 , STEREOGRAPHIC , 16 > mnt ( extent , src ) ;
        env = mnt ;
        get_mask = [=] ( const in_v & crd3 )
          { return mnt.get_mask ( crd3 ) ; } ;
        break ;
      }
      case FISHEYE:
      {
        mount_t < C , 3 , FISHEYE , 16 > mnt ( extent , src ) ;
        env = mnt ;
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

  environment()
  {
    // we have two different code paths to deal with. The first one
    // is for 'mounted' images in various projections. They require
    // specific code to deal with the projection and mask out pixels
    // which aren't covered by the source image data.

//     if ( args.mount_image != std::string() )
//     {
//       shape_type shape { args.mount_width , args.mount_height } ;
// 
//       if ( args.verbose )
//         std::cout << "processing mounted image " << args.mount_image
//                   << " shape " << shape << std::endl ;
// 
//       zimt::bc_code bc0 = zimt::REFLECT ;
//       if ( args.mount_prj == SPHERICAL || args.mount_prj == CYLINDRICAL )
//       {
//         if ( fabs ( args.mount_hfov - 2.0 * M_PI ) < .000001 )
//           bc0 = PERIODIC ;
//       }
//       p_bspl.reset ( new spl_t ( shape , args.spline_degree ,
//                                   { bc0 , zimt::REFLECT } ) ) ;
// 
//       bool success = args.mount_inp->read_image (
//         0 , 0 , 0 , C ,
//         TypeDesc::FLOAT ,
//         p_bspl->core.data() ,
//         sizeof ( in_px_t ) ,
//         p_bspl->core.strides[1] * sizeof ( in_px_t ) ) ;
//       assert ( success ) ;
// 
//       if (    args.mount_prj == SPHERICAL
//            && fabs ( args.mount_hfov - 2.0 * M_PI ) < .000001
//            && args.mount_width == 2 * args.mount_height )
//       {
//         // assume the mounted image is a full spherical
//         p_bspl->spline_degree = args.prefilter_degree ;
//         spherical_prefilter ( *p_bspl , p_bspl->core , zimt::default_njobs ) ;
//         p_bspl->spline_degree = args.spline_degree ;
//       }
//       else
//       {
//         // assume it's a partial spherical, use ordinary prefilter
//         p_bspl->spline_degree = args.prefilter_degree ;
//         p_bspl->prefilter() ;
//         p_bspl->spline_degree = args.spline_degree ;
//       }
// 
//       source_t < C , 2 , 16 > src ( p_bspl ) ;
// 
//       // for now, we mount images to the center; other types of cropping
//       // might be added by providing suitable parameterization.
//   
//       auto extent = get_extent ( args.mount_prj , args.mount_width ,
//                                  args.mount_height , args.mount_hfov  ) ;
// 
//       // we fix the projection as a template argument to class mount_t.
// 
//       switch ( args.mount_prj )
//       {
//         case RECTILINEAR:
//         {
//           mount_t < C , 3 , RECTILINEAR , 16 > mnt ( extent , src ) ;
//           env = mnt ;
//           get_mask = [=] ( const in_v & crd3 )
//             { return mnt.get_mask ( crd3 ) ; } ;
//           break ;
//         }
//         case SPHERICAL:
//         {
//           mount_t < C , 3 , SPHERICAL , 16 > mnt ( extent , src ) ;
//           env = mnt ;
//           get_mask = [=] ( const in_v & crd3 )
//             { return mnt.get_mask ( crd3 ) ; } ;
//           break ;
//         }
//         case CYLINDRICAL:
//         {
//           mount_t < C , 3 , CYLINDRICAL , 16 > mnt ( extent , src ) ;
//           env = mnt ;
//           get_mask = [=] ( const in_v & crd3 )
//             { return mnt.get_mask ( crd3 ) ; } ;
//           break ;
//         }
//         case STEREOGRAPHIC:
//         {
//           mount_t < C , 3 , STEREOGRAPHIC , 16 > mnt ( extent , src ) ;
//           env = mnt ;
//           get_mask = [=] ( const in_v & crd3 )
//             { return mnt.get_mask ( crd3 ) ; } ;
//           break ;
//         }
//         case FISHEYE:
//         {
//           mount_t < C , 3 , FISHEYE , 16 > mnt ( extent , src ) ;
//           env = mnt ;
//           get_mask = [=] ( const in_v & crd3 )
//             { return mnt.get_mask ( crd3 ) ; } ;
//           break ;
//         }
//         default:
//         {
//           std::cerr << "unknown projection: "
//                     << args.mount_prj_str << std::endl ;
//           assert ( false ) ;
//           break ;
//         }
//       }
//     }
//     else
    {
      // this code path is for 'true' environments which provide pixel
      // data for the entire 360X180 degrees. For now, we accept either
      // a 2:1 full equirect or a 1:6 vertically stacked cubemap.

      auto & inp ( args.inp ) ;
      auto & w ( args.env_width ) ;
      auto & h ( args.env_height ) ;

      shape_type shape { w , h } ;

      if ( w == 2 * h )
      {
        if ( args.verbose )
          std::cout << "environment has 2:1 aspect ratio, assuming latlon"
                    << std::endl ;

        p_bspl.reset ( new spl_t ( shape , args.spline_degree ,
                                   { zimt::PERIODIC , zimt::REFLECT } ) ) ;
        bool success = inp->read_image (
          0 , 0 , 0 , C ,
          TypeDesc::FLOAT ,
          p_bspl->core.data() ,
          sizeof ( in_px_t ) ,
          p_bspl->core.strides[1] * sizeof ( in_px_t ) ) ;
        assert ( success ) ;

        p_bspl->spline_degree = args.prefilter_degree ;
        spherical_prefilter ( *p_bspl , p_bspl->core , zimt::default_njobs ) ;
        p_bspl->spline_degree = args.spline_degree ;
        auto bsp_ev = zimt::make_evaluator < spl_t , float , LANES > ( *p_bspl ) ;

        env = ray_to_ll_t() + ll_to_px_t ( h ) + bsp_ev ;
      }
      else if ( h == 6 * w )
      {
        if ( args.verbose )
        {
          std::cout << "environment has 1:6 aspect ratio, assuming cubemap"
                    << std::endl ;
        }
        if ( args.cbm_prj == CUBEMAP )
        {
          sixfold_t < C , CUBEMAP > sf ( args.env_width , args.cbmfov ,
                                         args.support_min , args.tile_size ) ;
          if ( args.multiple_input )
            sf.load ( args.cfs.get_filenames() ) ;
          else
            sf.load ( inp ) ;
          inp->close() ;

          // env =   cbm_to_px_t < C , CUBEMAP > ( sf ) ;
          env = sf ;
        }
        else
        {
          sixfold_t < C , BIATAN6 > sf ( args.env_width , args.cbmfov ,
                                         args.support_min , args.tile_size ) ;
          if ( args.multiple_input )
            sf.load ( args.cfs.get_filenames() ) ;
          else
            sf.load ( inp ) ;
          inp->close() ;

          // env =   cbm_to_px_t < C , BIATAN6 > ( sf ) ;
          env = sf ;
        }
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

  mask_t masked_eval ( const in_v & in , out_v & out )
  {
    mask_t mask = get_mask() ;
    env.eval ( in , out ) ;
    return mask ;
  }

} ;

// environment9 objects mediate lookup with derivatives. Their eval
// member function takes 'ninepacks' and provides pixels with
// nchannels channels. We have two variants: The first one,
// which is created when p_env is passed a pointer to an
// 'environment' object, is used with 'twining'. The second
// one, receiving nullptr, is used with OIIO's texture system
// and does not use an environment object.

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

  typedef environment < float , float , nchannels , 16 > env_t ;

  environment9 ( env_t * p_env = nullptr )
  {
    if ( args.itp == -2 )
    {
      // with itp -2, we expect a non-nullptr p_env argument

      assert ( p_env != nullptr ) ;

      // set up the twining filter

      // std::vector < zimt::xel_t < float , 3 > > spread ;
      // 
      // if ( args.twf_file != std::string() )
      // {
      //   // if user passes a twf-file, it's used to set up the twining
      //   // filter, and all the other twining-related arguments apart
      //   // from twine_width and twine_normalize are ignored.
      // 
      //   read_twf_file ( spread ) ;
      //   assert ( spread.size() ) ;
      // }
      // else
      // {
      //   if ( args.twine == 0 )
      //   {
      //     // user has passed twine 0, or not passed anything - but itp
      //     // is -2. We set up the twining with automatically generated
      //     // parameters.
      // 
      //     // figure out the magnification in the image center as a
      //     // guideline
      // 
      //     double mag = args.env_step / args.step ;
      // 
      //     if ( mag > 1.0 )
      //     {
      //       // if the transformation magnifies, we use a moderate twine
      //       // size and a twine_width equal to the magnification, to
      //       // avoid the star-shaped artifacts from the bilinear
      //       // interpolation used for the lookup of the contributing
      //       // rays. If mag is small, the star-shaped artifacts aren't
      //       // really an issue, and twine values beyond, say, five
      //       // do little to improve the filter response, so we cap
      //       // the twine value at five, but lower it when approaching
      //       // mag 1, down to two - where the next case down starts.
      // 
      //       args.twine = std::min ( 5 , int ( 1.0 + mag ) ) ;
      //       args.twine_width = mag ;
      //     }
      //     else
      //     {
      //       // otherwise, we use a twine size which depends on the
      //       // downscaling factor (reciprocal of 'mag') and a twine_width
      //       // of 1.0: we only want anti-aliasing. picking a sufficiently
      //       // large twine value guarantees that we're not skipping any
      //       // pixels in the source (due to undersampling).
      // 
      //       args.twine = int ( 1.0 + 1.0 / mag ) ;
      //       args.twine_width = 1.0 ;
      //     }
      // 
      //     if ( args.twine_density != 1.0f )
      //     {
      //       // if the user has passed twine_density, we use it as a
      //       // multiplicative factor to change args.twine - typically
      //       // twine_density will be larger than one, so we'll get
      //       // more filter taps.
      // 
      //       double twine = args.twine * args.twine_density ;
      //       args.twine = std::round ( twine ) ;
      //     }
      // 
      //     if ( args.verbose )
      //     {
      //       std::cout << "automatic twining for magnification " << mag
      //                 << ":" << std::endl ;
      //       std::cout << "twine: " << args.twine
      //                 << " twine_width: " << args.twine_width
      //                 << std::endl ;
      //     }
      //   }
      //   else
      //   {
      //     if ( args.verbose )
      //     {
      //       std::cout << "using fixed twine: " << args.twine
      //                 << " twine_width: " << args.twine_width
      //                 << std::endl ;
      //     }
      //   }
      // 
      //   // with the given twining parameters, we can now set up a 'spread':
      //   // the generalized equivalent of a filter kernel. While a 'standard'
      //   // convolution kernel has pre-determined geometry (it's a matrix of
      //   // coefficients meant to be applied to an equally-shaped matrix of
      //   // data) - here we have three values for each coefficient: the
      //   // first two define the position of the look-up relative to the
      //   // 'central' position, and the third is the weight and corresponds
      //   // to a 'normal' convolution coefficient.
      // 
      //     make_spread ( spread , args.twine , args.twine ,
      //                   args.twine_width , args.twine_sigma ,
      //                   args.twine_threshold ) ;
      // }

      // wrap the 'environment' object in a twine_t object and assign
      // to 'act' - this 'groks' the twine_t object to act's type.
      // Note how this object is re-created for every run: the twining
      // parameters may change due to changing hfov from one invocation
      // of 'work' to the next.

      if ( args.twine_precise )
      {
        act = twine_t < nchannels , 16 , true >
          ( *p_env , args.twine_spread ) ;
      }
      else
      {
        act = twine_t < nchannels , 16 , false >
          ( *p_env , args.twine_spread ) ;
      }
    }
    // else
    // {
    //   assert ( args.itp == -1 ) ;
    // 
    //   // with itp -1, we expect a nullptr p_env argument
    // 
    //   assert ( p_env == nullptr ) ;
    // 
    //   // use OIIO's 'environment' or 'texture' lookup functions.
    //   // These code paths are coded in the 'latlon' and 'cubemap'
    //   // objects, which are created and also grokked to 'act'.
    //   // There, the 'ninepack' is used to calculate the derivatives
    //   // and then all the data needed for the look-up are passed to
    //   // the relevant (batched) OIIO look-up function. mounted
    //   // images use an object of type 'mount_t' which, in turn,
    //   // uses an 'inner functor' of type source_t. So we first
    //   // set up the source_t object, then create the projection-
    //   // -specific mount_t object.
    // 
    //   // if ( args.mount_image != std::string() )
    //   // {
    //   //   // note how we create a source_t with no arguments to the c'tor.
    //   //   // this routes to the code using OIIO. The parameters to set up
    //   //   // the source_t object are taken from 'args'.
    //   // 
    //   //   source_t < nchannels , 6 , 16 > src ;
    //   // 
    //   //   auto extent = get_extent ( args.mount_prj , args.mount_width ,
    //   //                              args.mount_height , args.mount_hfov  ) ;
    //   // 
    //   //   // same case switch as for the single-coordinate code path;
    //   //   // might be factored out - for now we just copy and paste.
    //   //   // note the second template argument, '9'. This produces
    //   //   // an object evaluating 'ninepacks' - three sets of 3D
    //   //   // coordinates.
    //   // 
    //   //   switch ( args.mount_prj )
    //   //   {
    //   //     case RECTILINEAR:
    //   //     {
    //   //       mount_t < nchannels , 9 , RECTILINEAR , 16 >
    //   //         mnt ( extent , src ) ;
    //   //       act = mnt ;
    //   //       break ;
    //   //     }
    //   //     case SPHERICAL:
    //   //     {
    //   //       mount_t < nchannels , 9 , SPHERICAL , 16 >
    //   //         mnt ( extent , src ) ;
    //   //       act = mnt ;
    //   //       break ;
    //   //     }
    //   //     case CYLINDRICAL:
    //   //     {
    //   //       mount_t < nchannels , 9 , CYLINDRICAL , 16 >
    //   //         mnt ( extent , src ) ;
    //   //       act = mnt ;
    //   //       break ;
    //   //     }
    //   //     case STEREOGRAPHIC:
    //   //     {
    //   //       mount_t < nchannels , 9 , STEREOGRAPHIC , 16 >
    //   //         mnt ( extent , src ) ;
    //   //       act = mnt ;
    //   //       break ;
    //   //     }
    //   //     case FISHEYE:
    //   //     {
    //   //       mount_t < nchannels , 9 , FISHEYE , 16 >
    //   //         mnt ( extent , src ) ;
    //   //       act = mnt ;
    //   //       break ;
    //   //     }
    //   //     default:
    //   //     {
    //   //       assert ( false ) ;
    //   //       break ;
    //   //     }
    //   //   }
    //   // }
    //   // else
    //   if ( args.env_width == args.env_height * 2 )
    //   {
    //     // set up an environment object picking up pixel values from
    //     // a lat/lon image using OIIO's 'environment' function.
    // 
    //     typedef latlon < float , float , nchannels , 16 > env_t ;
    //     env_t * p_env = (env_t*) current_env.has ( args.input ) ;
    // 
    //     if ( ! p_env )
    //     {
    //       current_env.clear() ;
    //       p_env = new env_t() ;
    //       current_env.reset ( args.input , p_env ) ;
    //     }
    //     act = *p_env ;
    //   }
    //   else
    //   {
    //     // set up an environment object picking up pixel values from
    //     // a texture representing the cubemap image, using OIIO's
    //     // 'texture' function
    // 
    //     typedef cubemap < float , float , nchannels , 16 > env_t ;
    //     env_t * p_env = (env_t*) current_env.has ( args.input ) ;
    // 
    //     if ( ! p_env )
    //     {
    //       current_env.clear() ;
    //       p_env = new env_t() ;
    //       current_env.reset ( args.input , p_env ) ;
    //     }
    //     act = *p_env ;
    //   }
    // }
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
