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

// This header provides class cubemap_t, an internal representation
// of a cubemap with cube faces in rectilinear or 'biatan6' projection.
// The latter projection adds an in-plane reprojection to the cube faces
// to minimize distrortion near the edges. The cube face images are
// arranged vertically in a single stripe image, with each individual
// cube face surrounded with additional 'support', which is simply an
// extension of the cube face's plane, filled in with data from
// adjoining cube faces. The extra support is needed for interpolators
// beyond nearest-neighbour (so, spline degree greater than zero) to
// allow for interpolation near the edges without special-casing. The
// internal representation uses class bspline, and interpolation of
// pixel values will be done with b-spline interpolation.

#if defined(ENVUTIL_CUBEMAP_H) == defined(HWY_TARGET_TOGGLE)
  #ifdef ENVUTIL_CUBEMAP_H
    #undef ENVUTIL_CUBEMAP_H
  #else
    #define ENVUTIL_CUBEMAP_H
  #endif

#include <OpenImageIO/imageio.h>

#include "zimt/prefilter.h"
#include "zimt/bspline.h"
#include "zimt/eval.h"

#include "geometry.h"

HWY_BEFORE_NAMESPACE() ;
BEGIN_ZIMT_SIMD_NAMESPACE(project)

using OIIO::ImageInput ;
using OIIO::TypeDesc ;
using OIIO::ImageSpec ;

// this structure is used to calculate the metrics of a sixfold_t.
// These values depend on four input values: the tile size, the
// size of a cube face image, as found in the input image(s), the
// horizontal field of view of a cube face image, and a minimum
// value for the size of the 'support'. The first two are in units of
// pixels, the third is a floating point value. The support size
// (in pixel units) is the minimal width of the surrounding frame
// which we need to set up good-quality interpolators.

struct metrics_t
{
  // the metrics_t object describes the geometry of the internal
  // representation (the 'IR image') of a cubemap as we use it
  // in the zimt examples dealing with cubemaps, and this structure
  // is also intended for use in lux.
  // we state a few preliminary rules: the square shapes which
  // constitute the geometrical object we describe are symmetric
  // both in the horizontal and the vertical, and each of the
  // six square 'sections' which are 'stacked' vertically in
  // the whole IR image are also symmetrical. The sections have
  // three distinct features: their total extent, the extent of
  // the image found in the input, and the extent of the central
  // section representing a 'cube face proper', which we refer to
  // as the 'core'. The first two are discrete values, the third
  // one may be discrete, so it's not guaranteed that the core has
  // discrete extent, but if it does, we set a flag (discrete90).
  // We add one more requirement: the cube face images must have
  // even width. This is to enforce perfect symmetry: the extent
  // of the entire IR image is a multiple of the tile size which,
  // in turn, must be a power of two. If we were to allow odd
  // sizes for the cube face images, their center would not be at
  // the same distance from the edges of the IR image they are
  // embedded in. We may be able to relax this requirement later
  // on, but as it stands, we rely on it.

  // we use variable names ending in _px to indicate discrete
  // values in pixel units, and variables ending in _md to
  // indicate floating point values in model space units.

  // the face width is the width (and height, both must agree)
  // of all the six square images in the input. This is a given
  // value which is received by the c'tor

  const std::size_t face_px ;

  // the face field of view is the field of view which the square
  // images correspond to. This is also a given, passed to the c'tor.
  // the field of view is understood as the angle from the left edge
  // of the leftmost pixel to the right edge of the rightmost pixel,
  // understanding pixels as small square areas of uniform colour
  // which have a size of one pixel unit. This indicates that if
  // the incoming images have their outermost pixels coinciding
  // precisely with the virtual cube's edges, they are slightly
  // larger than ninety degrees if we stick to the notion of the
  // field of view given above.

  const double face_fov ;

  // the width of a section must be a multiple of the tile width,
  // so we have section_px == n_tiles * tile_px. The tile
  // width is also a given, passed to the c'tor, but n_tiles and
  // section_px are calculated during set-up.

  const std::size_t tile_px ;

  std::size_t n_tiles ;
  std::size_t section_px ;
  std::size_t offset_px ;

  // size of the frame of pixels which is put around the cube face
  // image. left_frame_px is also the width of the frame on top
  // of the cube face image, right_frame_px also for the bottom.
  // With even face_px, both values are equal, with odd face_px,
  // they differ by one.

  std::size_t left_frame_px ;
  std::size_t right_frame_px ;

  // the minimal support, in pixel units, is the amount of pixels
  // we require as additional 'headroom' around the central part which
  // covers the 'cube face proper' covering precisely ninety degrees.
  // To make a clear definition: It's the number of pixels from the
  // margin inwards which are wholly outside the 'cube face proper'.
  // The support is the difference between the section width and the
  // core size, truncated to integer. support_min_px is also a given,
  // passed to the c'tor

  const std::size_t support_min_px ;

  // this value gives the number of pixels in the cube face image,
  // going inwards from the edge, which are entirely outside the
  // central ninety degree wide 'cube face proper'.

  std::size_t inherent_support_px ;

  // in the metrics_t object, we hold two scaling factors, which
  // can be used to translate between pixel units and 'model space
  // units'. 'model space units' refer to the IR image 'draped'
  // on a plane at unit distance, where the sample points
  // are laid out so that the 3D viewing rays they concide with
  // intersect with the plane. The center of the entire IR image,
  // draped in this way, is between the third and fourth section
  // vertically, and conincides with the line through the centers
  // of all sections horizontally. A point on a cube face located
  // on a corner of a 'core' has unit coordinate values (tan (pi/4))
  // the factors are multiplicative factors, and the pair given
  // is reciprocal, so, model_to_px == 1.0 / px_to_model, within
  // floating point precision. This is to facilitate conversions,
  // where we wish to avoid divisions which are harder to compute.

  double model_to_px ;
  double px_to_model ;

  // for conversion between model space units and texture units, we
  // need separate factors for x and y axis due to the 1:6 aspect ratio

  zimt::xel_t < double , 2 > model_to_tx ;
  zimt::xel_t < double , 2 > tx_to_model ;

  // section_md is the extent of a section in model space units,
  // and we have: section_md == section_px * px_to_model.
  // Note that the 'size' of the core in model space units is
  // precisely 2.0: it's 2 * tan ( pi/4 ) - so there is no variable
  // 'core_md' because it's implicit. We don't give a size
  // corresponding to the cube face images in the input.

  double section_md ;

  // refc_md is the distance, in model space units, from the edge of
  // a section to the cube face's center.

  double refc_md ;

  // ref90_md is the distance, in model space units, from the edge of
  // a section to the core, the 'cube face proper'. This is redundant;
  // it's refc_md - 1.

  double ref90_md ;

  // discrete90 is a flag which indicates whether ref90_md is discrete
  // (in pixel units) - within floating point precision.

  bool discrete90 ;

  double overscan_md ;
  double radius_md ;
  double diameter_md ;

  // the c'tor only needs one argument: the size of an individual
  // cube face image. This is the size covering the entire image,
  // whose field of view may be ninety degrees or more. The default
  // is ninety degrees, as we would expect with 'standard'
  // cubemaps, but we want to cater for cube face images which
  // already contain some support themselves. A minimal support of
  // four is ample for most direct interpolators, and the tile width
  // of 64 is a common choice, allowing for a good many levels of
  // mip-mapping with simple 4:1 pixel binning.

  metrics_t ( std::size_t _face_px ,
              double _face_fov = M_PI_2 ,
              std::size_t _support_min_px = 4UL ,
              std::size_t _tile_px = 64UL
            )
  : face_px ( _face_px ) ,
    face_fov ( _face_fov ) ,
    support_min_px ( _support_min_px ) ,
    tile_px ( _tile_px )
  {
    // first make sure that certain minimal requirements are met

    // the cube face images must have at least 90 degrees fov

    assert ( face_fov >= M_PI_2 ) ;

    // the tile width must be at least one

    assert ( tile_px > 0 ) ;

    // the tile width must be a power of two

    assert ( ( tile_px & ( tile_px - 1 ) ) == 0 ) ;

    // given the face image's field of view, how much support does
    // the face image already contain? We start out by calculating
    // The cube face image's diameter (in model space units) and the
    // 'overscan' - by how much the diameter exceeds the 2.0 diameter
    // which occurs with a cube face image of precisely ninety degrees
    // field of view. If the cube face images cover precisely ninety
    // degrees, we leave the values as they are initialized here:

    overscan_md = 0.0 ;
    radius_md = 1.0 ;
    diameter_md = 2.0 ;
    inherent_support_px = 0 ;

    // calculate the diameter in model space units. The diameter
    // for a ninety degree face image would be precisely 2, but
    // if the partial image has larger field of view, it will be
    // larger.

    if ( face_fov > M_PI_2 )
    {
      radius_md = tan ( face_fov / 2.0 ) ;
      diameter_md = 2.0 * radius_md ;
      overscan_md = radius_md - 1.0 ;
    }

    // calculate scaling factors from model space units to pixel
    // units and from pixel units to model space units. This factor
    // is crucial: it depends on the width of the cube face image
    // and it's field of view. 'draping' the image to model space
    // maps it to a plane at unit distance from the origin, so
    // that vectors to the 'draped' pixels correspond with
    // directional vectors to the scene points the pixels represent.
    // An alternative approach would be to hold the sample step
    // constant - e.g. 1.0 - and 'drape' the image so far from
    // the origin that the same correspondence would hold true.
    // A plane at unit distance has the advantage of one component
    // of a 3D coordinate on it being precisely 1.0.

    model_to_px = double ( face_px ) / diameter_md ;
    px_to_model = diameter_md / double ( face_px ) ;

    // The diameter, expressed in pixels - is just the cube face
    // image's width, and it's a discrete value. The 'cuber face
    // proper' - the section covering precisely ninety degrees,
    // may correspond to a discrete value or not, but what's
    // more important, is whether it's edges map to discrete
    // value, so that we can simply copy this section out when
    // we store the cubemap to a file. We can find out about
    // this property by looking at the overscan: obviously,
    // if the overscan corresponds to a discrete number of
    // pixels, so does the central section.
    // if we truncate the overscan - expressed in pixel units
    // - to an integer, we get the number of pixels which the
    // cube face image covers but which are wholly outside of
    // the central ninety degree 'cube face proper', which is
    // our definition of the 'inherent support'.

    double px_overscan = model_to_px * overscan_md ;
    inherent_support_px = std::trunc ( px_overscan ) ;

    if ( px_overscan - std::trunc ( px_overscan ) < 0.0000001 )
    {
      discrete90 = true ;
    }
    else
    {
      discrete90 = false ;
    }

    // if there is at least as much inherent support as the required
    // minimal support, we don't need to provide additional support.
    // If the inherent support is too small, we do need additional
    // support.
  
    std::size_t additional_support_px = 0 ;

    if ( inherent_support_px < support_min_px )
    {
      additional_support_px = support_min_px - inherent_support_px ;
    }

    // given the additional support we need - if any - how many tiles
    // are needed to contain a cube face image and it's support?

    std::size_t px_min = face_px + 2 * additional_support_px ;

    // if the user passes a tile size of one, n_tiles will end up
    // equal to px_min.

    n_tiles = px_min / tile_px ;
    if ( n_tiles * tile_px < px_min )
      n_tiles++ ;

    // this gives us the 'outer width': the size, in pixels, which
    // the required number of tiles will occupy, and the frame size,
    // the number of pixels from the edge of the IR section to
    // the first pixel originating from the cube face image.
    // Note that left_frame_px and right_frame_px will differ by
    // one for odd cube face sizes - for even sizes, they will be
    // equal. If they differ, which of the two we select for the
    // left or right value is arbitrary.

    section_px = n_tiles * tile_px ;

    // what's the offset to the same location in the next section?

    offset_px = section_px * section_px ;

    std::size_t frame_total = section_px - face_px ;

    left_frame_px = frame_total / 2 ;
    right_frame_px = frame_total - left_frame_px ;

    // paranoid

    assert (    ( left_frame_px + right_frame_px + face_px )
             == section_px ) ;

    // we also want the section's extent in model space units:

    section_md = px_to_model * section_px ;

    // the central part of each partial image, which covers
    // precisely ninety degrees - the cube face proper - is at
    // a specific distance from the edge of the total section.
    // we know that it's precisely 1.0 from the cube face image's
    // center in model space units, so we calculate that first:

    double refc_px =   double ( left_frame_px )
                     + double ( face_px ) / 2.0 ;

    refc_md = px_to_model * refc_px ;
    ref90_md = refc_md - 1.0 ;

    // we have the section's extent in model space units, so now we
    // can also provide the conversion factor to texture units, which
    // is anisotropic because the texture is not square.

    model_to_tx[0] = 1.0 / section_md ;
    model_to_tx[1] = 1.0 / ( 6.0 * section_md ) ;

    tx_to_model[0] = section_md ;
    tx_to_model[1] = 6.0 * section_md ;
  }

  // get_pickup_coordinate receives the in-face coordinate and
  // the face index and produces the corresponding coordinate,
  // in pixel units, which is equivalent in the total internal
  // representation image. So this coordinate can be used to
  // 'pick up' the pixel value from the IR image by using the
  // desired interpolation method, and it can be picked up
  // with a single interpolator invocation without any case
  // switching or conditionals to choose the correct cube face
  // image - all the images are combined in the IR.

  // The function is coded as a template to allow for scalar
  // and SIMDized pickups alike. Incoming we have face index(es)
  // (in the range (0, 5)) and in-face coordinate(s) in model
  // space units (in the range (-1, 1)), and outgoing we have
  // a coordinate in pixel units pertaining to the entire
  // IR image. Note that the outgoing value is in pixel units,
  // but it's not discrete.

  template < typename face_index_t , typename crd_t >
  void get_pickup_coordinate_px ( const face_index_t & face_index ,
                                  const crd_t & in_face_coordinate ,
                                  crd_t & target ) const
  {
    target =    in_face_coordinate // in (-1,1)
              + refc_md ;          // center of the cube face

    // we could add the per-face offset here, but see below:

    // target[1] += face_index * section_md ;

    // move from model space units to pixel units. This yields us
    // a coordinate in pixel units pertaining to the section.

    target *= model_to_px ;

    // add the per-section offset - Doing this after the move to
    // pixel units makes the calculation more precise, because
    // we can derive the offset from integer values. Here we reap
    // the benefit of using a fixed layout: the face index(es)
    // translate neatly into an offset/offsets in a single SIMD
    // multiplication, without any more case-switching.

    target[1] += ( face_index * int(section_px) ) ;

    // Subtract 0.5 - we look at pixels as small squares with an
    // extent of 1 pixel unit, and an incoming coordinate which
    // is precisely on the margin of the cube face (a value of
    // +/- 1) has to be mapped to the outermost pixel's margin as
    // well. The output of this function produces coordinate
    // values which coincide with the discrete pixel indices,
    // So a value of (0,0) refers to the upper left pixel's
    // center, and direct interpolation would yield this pixel's
    // value precisely.
    // Note that the output is now in the range (-0.5, width-0.5)
    // and using this with interpolators needing support may
    // access pixels outside the 'ninety degrees proper'.
    // We usually provide support to cater for that, but its'
    // good to keep the fact in mind. Even nearest-neighbour
    // pickup might fall outside the 'ninety degrees proper'
    // due to small imprecisions and subsequent rounding.

    target -= .5f ;
  }

  // variant to get the pick-up coordinate in model space units.
  // This isn't likely to be useful, but I leave it in for now.

  template < typename face_index_t , typename crd_t >
  void get_pickup_coordinate_md ( const face_index_t & face_index ,
                                  const crd_t & in_face_coordinate ,
                                  crd_t & target ) const
  {
    target =    in_face_coordinate // in (-1,1)
              + refc_md ;          // center of the cube face

    // we add the per-face offset and we're done.

    // problem on my mac: can't multiply int vector and double

    target[1] += ( face_index * float ( section_md ) ) ;
  }

  // variant of get_pickup_coordinate which yields the pickup
  // coordinate in texture units with extent in [0,1]. Same as
  // above, with an added multiplication to move from model space
  // coordinates to texture coordinates, taking into account the
  // non-square shape of the IR texture.

  template < typename face_index_t , typename crd_t >
  void get_pickup_coordinate_tx ( const face_index_t & face_index ,
                                  const crd_t & in_face_coordinate ,
                                  crd_t & target ) const
  {
    target =    in_face_coordinate // in (-1,1)
              + refc_md ;          // center of the cube face

    // problem on my mac: can't multiply int vector and double

    target[1] += ( face_index * float ( section_md ) ) ;

    // move from model space units to pixel units. This yields us
    // a coordinate in pixel units pertaining to the section.

    target *= model_to_tx ;
  }
} ;

// cubemap_t provides an internal representation of a cubemap
// with widened support, to provide for easy mip-mapping and
// interpolation. It's also a zimt::unary_functor providing
// pixel values for 3D ray coordinates. Since a cubemap can
// provide pixel values for all directions, there is no code
// to provide masks.

template < std::size_t nchannels , projection_t cbm_prj >
struct cubemap_t
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

  zimt::view_t < 2 , px_t > store ;

  // pointer to the upper left corner of the topmost cubeface
  // image inside the 'store' array

  px_t * p_ul ;

  typedef zimt::bspline < px_t , 2 > spl_t ;
  std::shared_ptr < spl_t > p_bsp ;
  zimt::grok_type < crd2_t , px_t , LANES > bsp_ev ;

  cubemap_t ( std::size_t _face_px ,
              double _face_fov = M_PI_2 ,
              std::size_t _support_min_px = 4UL ,
              std::size_t _tile_px = 64UL )
  : metrics_t ( _face_px ,
                _face_fov ,
                _support_min_px ,
                _tile_px )
    // p_bsp ( std::make_shared < spl_t >
    //           ( { section_px , 6 * section_px } , true ,
    //             args.spline_degree ,
    //             { zimt::REFLECT , zimt::REFLECT } ) )
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

/*    
    auto * ps = new spl_t ( store , store , args.spline_degree ,
                            { zimt::REFLECT , zimt::REFLECT } ,
                            -1 , left_frame_px ) ;
    p_bsp.reset ( ps ) ;*/
    auto * ps = new spl_t ( { section_px , 6 * section_px } , true ,
                            args.spline_degree ,
                            { zimt::REFLECT , zimt::REFLECT } ) ;
    p_bsp.reset ( ps ) ;
    store = p_bsp->container ;
    p_ul = store.data() ;
    p_ul += ( left_frame_px * store.strides ) . sum() ;
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

      int cmin = left_frame_px > 0 ? -1 : 0 ;
      int cmax = right_frame_px > 0 ? face_px : face_px - 1 ;

      // we use 2D discrete coordinates

      typedef zimt::xel_t < int , 2 > ix_t ;

      // mirror the horizontal edges

      for ( int x = cmin ; x <= cmax ; x++ )
      {
        ix_t src { x ,  0 } ;
        ix_t trg { x , -1 } ;
        if ( left_frame_px )
          cubeface [ trg ] = cubeface [ src ] ;
        src [ 1 ] = face_px - 1 ;
        trg [ 1 ] = face_px ;
        if ( right_frame_px )
          cubeface [ trg ] = cubeface [ src ] ;
      }

      // and the vertical edges

      for ( int y = cmin ; y <= cmax ; y++ )
      {
        ix_t src {  0 , y } ;
        ix_t trg { -1 , y } ;
        if ( left_frame_px )
          cubeface [ trg ] = cubeface [ src ] ;
        src [ 0 ] = face_px - 1 ;
        trg [ 0 ] = face_px ;
        if ( right_frame_px )
          cubeface [ trg ] = cubeface [ src ] ;
      }
    }
  }

  // We set up an internal functor which we'll use to fill the frame
  // of support pixels. This functor is local to struct cubemap_t
  // because it has no practical use outside of this context and
  // is tailor-made for the purpose (the member function fill_support
  // just below it)

  // this functor is used to fill the frame of support pixels in the
  // array in the cubemap_t object. incoming, we have 2D coordinates,
  // which, for the purpose at hand, will lie outside the cube face.
  // But we can still convert these image coordinates to planar
  // coordinates in 'model space' and then further to 'pickup'
  // coordinates into the array in the cubemap_t object. When we
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
    const cubemap_t < nchannels , cbm_prj > & sf ;
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

    fill_frame_t ( const cubemap_t < nchannels , cbm_prj > & _cubemap ,
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
      // Note how we use cubemap_t's member functions for most of
      // the work - gleaning pixel values from the cubemap_t object
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
    if ( left_frame_px == 0 && right_frame_px == 0 )
      return ;

    mirror_around() ;

    auto * p_base = store.data() ;

    for ( int face = 0 ; face < 6 ; face++ )
    {
      // set up the 'gleaning' functor

      fill_frame_t fill_frame ( *this , face ) ;
    
      // get a pointer to the upper left of the section of the IR

      auto * p_frame = p_base + face * section_px * store.strides[1] ;
      
      // we form a view to the current section of the array in the
      // cubemap_t object. The shape of this view is the 'notional
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
        // cubemap_t object. The shape of this view is the 'notional
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
    // appropriate slots in the cubemap_t object's 'store' array.

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
    // the cubemap_t object's 'store' array.

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
  // of the cubemap held in the cubemap_t object. We have two
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
} ; // end of struct cubemap_t

END_ZIMT_SIMD_NAMESPACE
HWY_AFTER_NAMESPACE() ;

#endif // sentinel
