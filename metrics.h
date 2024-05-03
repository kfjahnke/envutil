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

#include "common.h"

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
  // and in-face coordinate(s) in model space units, and outgoing
  // we have a coordinate in pixel units pertaining to the entire
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
    // we can derive the offset from integer values.

    target[1] += f_v ( face_index * int(section_px) ) ;

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

  // variant to get the pick-up coordinate in model space units

  template < typename face_index_t , typename crd_t >
  void get_pickup_coordinate_md ( const face_index_t & face_index ,
                                  const crd_t & in_face_coordinate ,
                                  crd_t & target ) const
  {
    target =    in_face_coordinate // in (-1,1)
              + refc_md ;          // center of the cube face

    // we add the per-face offset and we're done.

    target[1] += face_index * section_md ;
  }

  // variant of get_pickup_coordinate which yields the pickup
  // coordinate in texture units with extent in [0,1].

  template < typename face_index_t , typename crd_t >
  void get_pickup_coordinate_tx ( const face_index_t & face_index ,
                                  const crd_t & in_face_coordinate ,
                                  crd_t & target ) const
  {
    target =    in_face_coordinate // in (-1,1)
              + refc_md ;          // center of the cube face

    target[1] += face_index * section_md ;

    // move from model space units to pixel units. This yields us
    // a coordinate in pixel units pertaining to the section.

    target *= model_to_tx ;
  }
} ;
