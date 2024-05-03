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

// utility program producing a lat/lon 'environment map' from a cubemap.
// and the reverse.

// AFAICT there are two formats used to represent a complete 360X180
// degree environment. The more common one is a lat/lon environment,
// which captures the environment in a single image in spherical
// projection. Using this format is quite straightforward, but it
// requires using transcendental functions to move between the
// lat/lon spherical coordinates and 3D 'ray' geometry. The second
// format, which - I have been told - is less common, is the
// 'cubemap' or 'sky box' format. It captures the environment in
// six square images representing the faces of a virtual cube
// surrounding the origin. For viewing purposes, this format has
// some advantages and some disadvantages. On the plus side is the
// fact that, since the cube faces are stored in rectilinear projection,
// reprojection to a rectilinear view can be done without transcendental
// functions. On the negative side, handling the six discrete images,
// which requires picking the right one to pick up image information
// to go to a specific target location in the output, requires a fair
// amount of logic, and the sequence and (for the top and bottom square)
// orientation of the cube faces isn't obvious and leaves room for error.
// There are standards, though - openEXR defines a cubemap format with
// specific sequence and orientations, and this is the convention I use
// in this program. openEXR offers a tool to convert between the two
// types of environment, but I found that it doesn't seem to work
// correctly - see this issue:
// https://github.com/AcademySoftwareFoundation/openexr/issues/1675
// I now think that this may be partly due to a difference in conception:
// Some documentation I've looked at seems to suggest that the cube
// face images openEXR expects are what I would consider slightly
// wider than ninety degrees: Their outermost pixels coincide with
// the cube's edges, whereas the cube faces I use have their outermost
// pixels half a sample step inwards from the cube's edges. The openEXR
// way has some merits - e.g. bilinear interpolation is immediately
// possible over the entire cube face - but I won't cater for that
// type of cube face in this program, and the images I render may
// not be 100% compatible to the openEXR standard.
// Since I am processing both types of environment representation in
// lux, I have a special interest in them and writing some image
// processing code in zimt+OpenImageIO offers a good opportunity to
// deepen my understanding and come up with efficient ways of dealing
// with both formats. I started out with 'cubemap.cc', which converts
// lat/lon format to openEXR-compatible cubemap format. To pick up
// data from a lat/lon environment map, there is ready-made code in
// OpenImageIO. I use this code to good effect, producing output
// with a very 'proper' anisotropic anti-aliasing filter. The
// reverse transformation - from a cubemap to lat/lon format - is
// what I've added in this program. Here, employing OpenImageIO for
// the task of picking up data from the environment is not available
// out-of-the-box - OpenImageIO does not support this format. This
// may well be because it is less used and harder to handle. In my
// opinion, the greatest problem with this format is the fact that the
// six cube face images each cover precisely ninety degrees. This is
// sufficient to regenerate the environment, but it's awkward, due
// to the fact that near the edges there is not enough correct
// support in the individual images to use good interpolators.
// Instead, to use such interpolators, this support has to be
// gleaned from adjoining cube faces, with proper reprojection.
// The second stumbling stone - when trying to use such cubemaps
// with OpenImageIO - is the fact that to use mip-mapping on them
// needs even larger support around the cube faces (unless they happen
// to have a size which is a multiple of the tile size) so as to
// avoid mixing data from cube faces which are located next to each
// other in the 1:6 stripe but have no other relation - there is
// a hard discontiuity from one image's bottom to the next image's
// top. The third stumbling stone is the fact that OpenImageIO's
// texture system code is file-based, and the data I produce from
// the 'raw' cubemap input to deal with the first two issues can't
// simply be fed back to OIIO's texture system, because they are
// in memory, whereas the texture system wants them on disk. I'd
// like to find a way to feed them directly, but for now I use
// an intermediate image on disk for the task.
//
// I have mentioned that I have dealt with the two first issues,
// so I'll give a quick outline here - the code is amply commented
// to explain the details. I start out by setting up a single array
// of pixel data which has enough 'headroom' to accomodate the cube
// faces plus a surrounding frame of support pixels which is large
// enough to allow for good interpolators and mip-mapping, so each
// cube face is embedded in a square section of the array which has
// a multiple of the tile size as it's extent. I proceed to import
// the cube face data, which I surround initially with a single-pixel
// frame of mirrored data to give enough support for bilinear
// interpolation even near the edges. Then I fill in the support frames
// using bilinear interpolation. The array is now filled with six square
// images which have more than ninety degrees field of view - the part
// beyond ninety degrees generated from adjoining cube faces. This is
// the texture from which I pick output data. Due to the ample support,
// it could be subjected to filters with large-ish support, and the
// array of six 'widened' cube faces might even stand as a useful
// format by itself, which would be quite easy to describe formally
// because it's a derivative of the openEXR cubemap layout.
// You can store this intermediate in a file called internal.exr
// by setting the global boolean 'save_ir' to true.
// The program has grown since it's inception to provide more code
// on the topic, which I used to verify that the conversion is
// correct and to be able to look at intermediate images. The
// central object, the 'sixfold_t', which holds the internal
// representation of the data, might be a good candidate to factor
// out into a separate TU. I also coded the conversion from a lat/lon
// environment into a cubemap - there is 'cubemap.cc' which also
// does that, but it's more a brief zimt demo program; the conversion
// I offer in this program is more flexible and comprehensive.

// I use some unconventional terminology: 'model space units' are
// coordinates pertaining to 'archetypal' 2D manifolds 'draped'
// in space. The image plane is draped at unit distance forward
// from the origin, and the image points are distributed on it
// so that their other two coordinates coincide with intersections
// of rays pointing at them. Using this scheme makes conversions
// easier and provides a common frame of reference. For spherical
// data, I use the surface of a sphere with unit radius as the
// 'archetypal' 2D manifold, and the image points are located
// where the rays pointing towards them intersect with this sphere.
// This results in the 'draped' image points residing at unit
// distance from the origin. A full spherical image is draped so
// that it's center coincides with the center of the image plane,
// and a rectilinear image is draped in the same way.
// I use lux coordinate system convention: x axis points right,
// y axis points down, z axis points forward. Where needed, I
// switch to openEXR axis convention - OIIO's lat/lon environment
// lookup expects 3D ray coordinates using opneEXR convention.
// 'simdized' values are in SoA format, so a simdized pixel is
// made up from three vectors with LANES elements each.

#include <array>
#include <filesystem>
#include "zimt/zimt.h"
#include <OpenImageIO/texture.h>
#include <OpenImageIO/filesystem.h>
#include "metrics.h"
#include "geometry.h"

using namespace OIIO ;

// some globals, which will be set via argparse in main

#include <regex>
#include <OpenImageIO/argparse.h>

static bool verbose = false;
static bool ctc = false ;
static bool help    = false ;
static bool store6 = false ;
static bool store6_lux = false ;
static std::string metamatch;
static std::regex field_re;
std::string input , output , save_ir , ts_options ;
int extent , support_min_px , tile_px ;
int itp , twine , twine_px ;
double face_fov , twine_sigma , twine_threshold ;

// helper function to save a zimt array to an image file. I am
// working with openEXR images here, hence the HALF data type.
// If the output format can't produce HALF, it will use the
// 'next-best' thing. Note that input and output should agree
// on the colour space. If you stay within one format, that's
// not an issue.

template < std::size_t nchannels >
void save_array ( const std::string & filename ,
                  const zimt::view_t
                    < 2 ,
                      zimt::xel_t < float , nchannels >
                    > & pixels ,
                  bool is_latlon = false )
{
  auto out = ImageOutput::create ( filename );
  assert ( out != nullptr ) ;
  ImageSpec ospec ( pixels.shape[0] , pixels.shape[1] ,
                    nchannels , TypeDesc::HALF ) ;
  out->open ( filename , ospec ) ;

  if ( is_latlon )
    ospec.attribute ( "textureformat" , "LatLong Environment" ) ;

  auto success = out->write_image ( TypeDesc::FLOAT , pixels.data() ) ;
  assert ( success ) ;
  out->close();
}

// helper function to create a set of cube face image names, by
// infixing an underscore, followed by the orientation, between
// the base name and the extension. While lux does not yet
// conform to openEXR standard re. the orientation there is a
// special case if store6_lux is set which names with the notion
// of the horizontal direction rotated by 180 degrees.

void six_names ( const std::string & name ,
                 std::vector < std::string > & name6 )
{
  auto point_pos = name.find_last_of ( "." ) ;
  assert ( point_pos != std::string::npos ) ;
  std::string base = name.substr ( 0 , point_pos ) ;
  std::string ext = name.substr ( point_pos ) ;
  name6.clear() ;
  if ( store6_lux )
  {
    name6.push_back ( base + "_right" + ext ) ;
    name6.push_back ( base + "_left" + ext ) ;
    name6.push_back ( base + "_top" + ext ) ;
    name6.push_back ( base + "_bottom" + ext ) ;
    name6.push_back ( base + "_back" + ext ) ;
    name6.push_back ( base + "_front" + ext ) ;
  }
  else
  {
    name6.push_back ( base + "_left" + ext ) ;
    name6.push_back ( base + "_right" + ext ) ;
    name6.push_back ( base + "_top" + ext ) ;
    name6.push_back ( base + "_bottom" + ext ) ;
    name6.push_back ( base + "_front" + ext ) ;
    name6.push_back ( base + "_back" + ext ) ;
  }
}

// as an augmentation which also can reduce alisasing
// if certain criteria are met, we introduce class twine_t.
// The signal is evaluated several times in the close vicinity
// of the pick-up point, the values are weighted and summed up.
// this functor wraps an inner functor, to which single-location
// lookups are routed. The inner functor doesn't have to be of
// very high quality - bilinear is perfectly sufficient, but
// twining is independent of the 'inner' interpolator. Using it
// with an interpolator which already does good antialiasing
// (like OIIO's default) is pointless, though, and will only
// result in very slight additional blur, so it's used with
// itp == 1, bilinear interpolation.
// With a 'spread' calculated with 'make_spread' (below), the
// overall result is precisely the same as oversampling the
// signal resulting from bilinear interpolation and subsequently
// applying a small box filter.
// Note that the parameterization allows not only for 'conventional'
// filters with evenly sampled kernel values, but for a generalized
// form of filter, which can locate the contributing pick-ups at
// arbitrary distance from the central ('un-twined') target
// coordinate and with arbitrary weights. This is not exploited
// here. Likely candidates would be e.g. circular patterns.

template < std::size_t nchannels >
struct twine_t
: public zimt::unary_functor
           < v2_t , zimt::xel_t < float , nchannels > , LANES >
{
  zimt::grok_type 
    < v2_t , zimt::xel_t < float , nchannels > , LANES > inner ;

  const std::vector < zimt::xel_t < float , 3 > > & spread ;
  
  twine_t ( const zimt::grok_type
               < v2_t , zimt::xel_t < float , nchannels > , LANES >
                 & _inner ,
             const std::vector < zimt::xel_t < float , 3 > > & _spread )
  : inner ( _inner ) ,
    spread ( _spread )
  { }

  template < typename in_type , typename out_type >
  void eval ( const in_type & in ,
              out_type & out )
  {
    out = 0 ;
    out_type px_k ;

    for ( auto const & contrib : spread )
    {
      auto in_k ( in ) ;
      in_k[0] += contrib[0] ;
      in_k[1] += contrib[1] ;
      inner.eval ( in_k , px_k ) ;
      out += contrib[2] * px_k ;
    }
  }
} ;

// this function sets up a simple box filter to use with the
// twine_t functor above. The output is the average of the
// contributing partial values. The given w and h values
// determine the number of pick-up points in the horizontal
// and vertical direction - typically, you'd use the same value
// for both. The deltas are set up so that, over all pick-ups,
// they produce a uniform sampling.
// additional parameters alow to apply gaussian weights and
// apply a threshold to suppress very small weighting factors;
// in a final step the weights are normalized.

void make_spread ( int w , int h , float d ,
                   float sigma , float threshold ,
                   std::vector < zimt::xel_t < float , 3 > > & trg )
{
  float wgt = 1.0 / ( w * h ) ;
  double x0 = - ( w - 1.0 ) / ( 2.0 * w ) ;
  double dx = 1.0 / w ;
  double y0 = - ( h - 1.0 ) / ( 2.0 * h ) ;
  double dy = 1.0 / h ;
  trg.clear() ;
  sigma *= - x0 ;
  double sum = 0.0 ;

  for ( int y = 0 ; y < h ; y++ )
  {
    for ( int x = 0 ; x < w ; x++ )
    {
      float wf = 1.0 ;
      if ( sigma > 0.0 )
      {
        double wx = ( x0 + x * dx ) / sigma ;
        double wy = ( y0 + y * dy ) / sigma ;
        wf = exp ( - sqrt ( wx * wx + wy * wy ) ) ;
      }
      zimt::xel_t < float , 3 >
        v { float ( d * ( x0 + x * dx ) ) ,
            float ( d * ( y0 + y * dy ) ) ,
            wf * wgt } ;
      trg.push_back ( v ) ;
      sum += wf * wgt ;
    }
  }

  double th_sum = 0.0 ;
  bool renormalize = false ;

  if ( sigma != 0.0 )
  {
    for ( auto & v : trg )
    {
      v[2] /= sum ;
      if ( v[2] >= threshold )
      {        
        th_sum += v[2] ;
      }
      else
      {
        renormalize = true ;
        v[2] = 0.0f ;
      }
    }
    if ( renormalize )
    {
      for ( auto & v : trg )
      {
        v[2] /= th_sum ;
      }
    }
  }

  if ( verbose )
  {
    if ( sigma != 0.0 )
    {
      std::cout << "using this twining filter kernel:" << std::endl ;
      for ( int y = 0 ; y < h ; y++ )
      {
        for ( int x = 0 ; x < w ; x++ )
        {
          std::cout << '\t' << trg [ y * h + x ] [ 2 ] ;
        }
        std::cout << std::endl ;
      }
    }
    else
    {
      std::cout << "using box filter for twining" << std::endl ;
    }
  }

  if ( renormalize )
  {
    auto help = trg ;
    trg.clear() ;
    for ( auto v : help )
    {
      if ( v[2] > 0.0f )
        trg.push_back ( v ) ;
    }
    if ( verbose )
    {
      std::cout << "twining filter taps after after thresholding: "
                << trg.size() << std::endl ;
    }
  }
}

// sixfold_t contains data for a six-square sky box with support
// around the cube faces to allow for proper filtering and
// interpolation. The content is set up following the layout
// of openEXR environment maps. cubemaps in openEXR contain six
// square images concatenated vertically. The sequence of the
// images, from top to bottom, is: left, right, top, bottom,
// front, back. The top and bottom images are oriented so that
// they align vertically with the 'back' image (lux aligns with
// the front image) - but note that labels like 'front' are
// somewhat arbitrary - I associate 'front' with the center of
// a full spherical image and 'left' with the first image in the
// cubemap. Note, again, that openEXR's 'own' cubemap format
// may use slightly larger cube faces than what I use here.
// The sixfold_t object also combines the six images in one
// array, but it adds a frame of additional 'support' pixels
// around each square image to allow for interpolation with
// interpolators needing support around the interpolation locus.
// The supporting frame around each square image is chosen to
// yield a total per-cube-face size which is a multiple of the
// tile size, to make mip-mapping easy. With a bit of indexing
// magic, we'll be able to form 'pick-up' coordinates pertaining
// directly to the entire array held in the sixfold_t, which is
// more efficient than combining pixel values from per-cube-face
// lookups - and we can code the entire process in SIMD code. The
// frames around the cube faces make it possible to do this with
// mip-mapping and correct anti-aliasing, using OIIO's planar
// texture lookup. As a fall-back, there's an implementation of
// simple bilinear interpolation directly from the IR image.
// The template argument 'nchannels' is the number of channels
// in the image. We accept up to four in main - RGBA should
// come in with associated alpha, so that we can formulate
// arithmetic with pixel values generically.

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

  TextureSystem * ts ;
  TextureSystem::TextureHandle * th ;

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
    store ( { section_px , 6 * section_px } ) ,
    ts ( nullptr ) ,
    th ( nullptr )
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

  // We can use the fall-back bilinear interpolation to pick up
  // pixel values from the IR image, but to use OIIO's 'texture'
  // function, we need access to the texture system and - for
  // speed - the texture handle. But we don't have these data
  // when the sixfold_t is created - the texture has to be
  // generated first, then stored to disk and then fed to the
  // texture system. So we can only call this function later:

  void set_ts ( TextureSystem * _ts ,
                TextureSystem::TextureHandle * _th )
  {
    ts = _ts ;
    th = _th ;
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

  void load ( std::unique_ptr<ImageInput> & inp )
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

  // store_cubemap stores a standard cubemap with cube faces
  // with a field of view of face_fov. This always corresponds
  // with a discrete set of pixels held in the store - that's
  // how the IR image is set up - so we can simply copy out
  // the data as they are.

  void store_cubemap ( const std::string & filename ) const
  {
    const int xres = face_px ;
    const int yres = 6 * face_px ;

    if ( verbose )
    {
      std::cout << "storing " << xres
                << " X " << yres << " cubemap" << std::endl ;
    }

    std::unique_ptr<ImageOutput> out
      = ImageOutput::create ( filename.c_str() ) ;
    assert ( out != nullptr ) ;

    // we'll store a few metadata along with the cubemap.
    // passing textureformat = "CubeFace Environment" may
    // be misleading if the cubemap isn't precisely to the
    // 'standard' specifications, so we also store two more
    // values, which indicate the field of view, and whether
    // the field of view was measured center-to-center or
    // edge-to-edge.

    ImageSpec spec ( xres , yres , nchannels , TypeDesc::HALF ) ;
    spec.attribute ( "textureformat" , "CubeFace Environment" ) ;
    double nominal_fov = face_fov ;
    if ( ctc )
    {
      double half_md = tan ( face_fov / 2.0 ) ;
      half_md *= ( xres / ( xres + 1.0 ) ) ;
      nominal_fov = atan ( half_md ) * 2.0 ;
    }
    spec.attribute ( "field-of-view" ,
                     float ( nominal_fov * 180.0 / M_PI ) ) ;
    spec.attribute ( "fov-measured" ,
                     ctc ? "center-to-center" : "edge-to-edge" ) ;
    out->open ( filename.c_str() , spec ) ;

    auto p_base = store.data() ;
    p_base += ( left_frame_px * store.strides ) . sum() ;

    for ( int face = 0 ; face < 6 ; face++ )
    {
      // virtual bool write_scanlines ( int ybegin , int yend , int z ,
      //                                TypeDesc format ,
      //                                const void * data ,
      //                                stride_t xstride = AutoStride ,
      //                                stride_t ystride = AutoStride )

      auto success =
      out->write_scanlines (   face * face_px ,
                             ( face + 1 ) * face_px ,
                               0 ,
                               TypeDesc::FLOAT ,
                               p_ul + face * offset_px ,
                               store.strides[0] * nchannels * 4 ,
                               store.strides[1] * nchannels * 4 ) ;

      assert ( success ) ;
    }

    out->close() ;
  }

  // overload of store_cubemap storing six separate cube face
  // images

  void store_cubemap ( const std::vector < std::string > & filename6 ) const
  {
    const int xres = face_px ;
    const int yres = 6 * face_px ;

    if ( verbose )
    {
      std::cout << "storing six " << xres << " X "
                << yres << " cube face images" << std::endl ;
    }

    auto p_base = store.data() ;
    p_base += ( left_frame_px * store.strides ) . sum() ;

    for ( int face = 0 ; face < 6 ; face++ )
    {
      std::unique_ptr<ImageOutput> out
        = ImageOutput::create ( filename6[face].c_str() ) ;
      assert ( out != nullptr ) ;

      ImageSpec spec ( xres , yres , nchannels , TypeDesc::HALF ) ;
      out->open ( filename6[face].c_str() , spec ) ;

      // virtual bool write_scanlines ( int ybegin , int yend , int z ,
      //                                TypeDesc format ,
      //                                const void * data ,
      //                                stride_t xstride = AutoStride ,
      //                                stride_t ystride = AutoStride )

      auto success =
      out->write_scanlines (   0 ,
                               face_px ,
                               0 ,
                               TypeDesc::FLOAT ,
                               p_base + face * offset_px ,
                               store.strides[0] * nchannels * 4 ,
                               store.strides[1] * nchannels * 4 ) ;

      assert ( success ) ;
      out->close() ;
    }
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
                          px_v & px ,
                          const int & degree ) const
  {
    // get the pick-up coordinate, in pixel units

    crd2_v pickup ;
    get_pickup_coordinate_px ( face , in_face , pickup ) ;

    if ( degree == 0 )
    {
      // simple nearest-neighbour lookup. This is not currently used,
      // but it's instructive.

      // convert the in-face coordinates to integer. If the in-face
      // coordinate is right on the edge, the pick-up may fall to
      // pixels outside the 'ninety degrees poper'.

      index_v idx { round ( pickup[0] ) , round ( pickup[1] ) } ;

      // calculate corresponding offsets, using the strides, and
      // obtain pixel data by gathering with this set of offsets.
      // Note how the first multiplication broadcasts the pair of
      // strides to the pair of int vectors

      const auto ofs = ( idx * store.strides ) . sum() * nchannels ;
      const auto * p = (float*) ( store.data() ) ;

      px.gather ( p , ofs ) ;
    }
    else if ( degree == 1 )
    {
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
                         TextureOptBatch & batch_options ) const
  {
    // code to truncate the pickup coordinate to int and gather,
    // after restoring the pickup to pixel units

    // pickup[0] *= float ( store.shape[0] ) ;
    // pickup[1] *= float ( store.shape[1] ) ;
    // index_v idx { pickup[0] , pickup[1] } ;
    // const auto ofs = ( idx * store.strides ) . sum() * nchannels ;
    // const auto * p = (float*) ( store.data() ) ;
    // 
    // px.gather ( p , ofs ) ;
    // return ;

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
    ts->texture ( th , nullptr , batch_options , Tex::RunMaskOn ,
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
    ts->texture ( th , nullptr , batch_options , Tex::RunMaskOn ,
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
                          TextureOptBatch & bo ) const
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

      sf.cubemap_to_pixel ( fv , in_face , px , degree ) ;
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

// ll_to_px_t is the functor used as 'act' functor for zimt::process
// to produce pixel values for lat/lon coordinates.
// This functor accesses the data in the sixfold_t object, yielding
// pixel data for incoming 2D lat/lon coordinates. It implements a
// typical processing pipeline, using three distinct stages: first,
// the incoming 2D lat/lon coordinates are converted to 3D 'ray'
// coordinates, then these are used to figure out the corresponding
// cube face and the 2D in-face coordinates, and finally these values
// are used to obtain pixel data from the sixfold_t object.
// Note how we code the entire process as a SIMD operation: we're
// not handling individual coordinates or pixels, but their 'SIMDized'
// equivalents. With this functor, we can glean pixel values for
// arbitrary discrete 2D coordinates located on any of the six
// planes which surround the origin at 1.0 units distance in model
// space.

template < std::size_t nchannels >
struct ll_to_px_t
: public zimt::unary_functor
   < v2_t , zimt::xel_t < float , nchannels > , LANES >
{
  // some types

  typedef zimt::xel_t < float , nchannels > px_t ;
  typedef zimt::simdized_type < px_t , LANES > px_v ;

  // infrastructure
  
  const ll_to_ray_t ll_to_ray ;

  // for lookup with OIIO's texture system, we need batch options.
  // I'd keep them in the code which actually uses them, but the
  // OIIO 'texture' function expects an lvalue.

  TextureOptBatch & batch_options ;

  // source of pixel data

  const sixfold_t < nchannels > & cubemap ;

  // parameter to choose the interpolator

  const int degree ;

  // sampling step

  v2_t delta ;

  // scaling factor to move from model space units to texture units
  // (separate for the s and t direction - these factors are applied
  // to coordinates pertaining to the IR image of the cubemap, which
  // has 1:6 aspect ratio)

  const float scale_s ;
  const float scale_t ;


  // ll_to_px_t's c'tor obtains a const reference to the sixfold_t
  // object holding pixel data, the degree of the interpolator and
  // the delta of the sampling step, to form derivatives.

  ll_to_px_t ( const sixfold_t < nchannels > & _cubemap ,
               const int & _degree ,
               const v2_t & _delta ,
               TextureOptBatch & _batch_options )
  : cubemap ( _cubemap ) ,
    degree ( _degree ) ,
    delta ( _delta ) ,
    batch_options ( _batch_options ) ,
    scale_s ( _cubemap.model_to_px / _cubemap.store.shape[0] ) ,
    scale_t ( _cubemap.model_to_px / -cubemap.store.shape[1] ) ,
    ll_to_ray() // g++ is picky.
  { }

  // 'eval' function which will be called by zimt::process.
  // incoming: 2D coordinate as lat/lon in radians. Internally, we
  // first calculate a 3D directional vector, then the index of the
  // cube face and then the 2D in-face coordinate, where we pick up
  // the pixel value from the cube face image. We form 2D pick-up
  // coordinates into the array held in the sixfold_t object, which
  // contains the images 'embedded' in additional support space,
  // which offers 'headroom' for interpolators which need support
  // and also provides each cube face image with just so much frame
  // that the sixth of the total array it inhabits can be tiled
  // exactly with the given tile size. This is to enable correct
  // mip-mapping for filtered texture lookup.

  void eval ( const crd2_v & lat_lon , px_v & px )
  {
    // convert 2D spherical coordinate to 3D ray coordinate 'crd3'
  
    crd3_v crd3 ;
    ll_to_ray.eval ( lat_lon , crd3 ) ;

    // find the cube face and in-face coordinate for 'c'

    i_v face ;
    crd2_v in_face ;
    ray_to_cubeface ( crd3 , face , in_face ) ;

    if ( degree == -1 )
    {
      // degree == -1 stands for "use OIIO's 'texture' function".
      // for this interpolation, we need the derivatives, in
      // pixel coordinates pertaining to the IR image. We use
      // the same approach as in cubemap.cc: we obtain source
      // image coordinates for two points, each one step in
      // one of the canonical directions away. Then we subtract
      // the unmodified pickup coordinate.

      // get a copy of the incoming 2D lat/lon coordinates and
      // add the sample step width to the horizontal component

      crd2_v dx1 = lat_lon ;
      dx1[0] += delta[0] ;

      // convert the results to 3D ray coordinates

      crd3_v dx1_3 ;
      ll_to_ray.eval ( dx1 , dx1_3 ) ;

      // and then to in-face coordinates, using the same face
      // indices we gleaned above. This is important: if we were
      // to compare in-face coordinates from different cube faces,
      // we'd get totally wrong results.

      crd2_v dx1_if ;
      ray_to_cubeface_fixed ( dx1_3 , face , dx1_if ) ;

      // now we can calculate the approximation of the derivative
      // by forming the difference from the in-face coordinate of
      // the pick-up location. We get two components, which we
      // handle separately, and we scale to texture coordinates.

      auto dsdx = ( dx1_if[0] - in_face[0] ) * scale_s ;
      auto dtdx = ( dx1_if[1] - in_face[1] ) * scale_t ;

      // we repeat the process for a coordinate one sample step away
      // along the vertical axis

      crd2_v dy1 = lat_lon ;
      dy1[1] += delta[1] ;
      crd3_v dy1_3 ;
      ll_to_ray.eval ( dy1 , dy1_3 ) ;
      crd2_v dy1_if ;
      ray_to_cubeface_fixed ( dy1_3 , face , dy1_if ) ;

      auto dsdy = ( dy1_if[0] - in_face[0] ) * scale_s ;
      auto dtdy = ( dy1_if[1] - in_face[1] ) * scale_t ;

      // now we can call the cubemap_to_pixel variant which is based
      // on OIIO's texture lookup and takes derivatives

      cubemap.cubemap_to_pixel ( face , in_face , px ,
                                 dsdx , dtdx , dsdy , dtdy ,
                                 batch_options ) ;
    }
    else
    {
      // use 'face' and 'in_face' to obtain pixel values directly
      // from the IR image with NN (degree 0) or bilinear (degree 1)
      // interpolation.

      cubemap.cubemap_to_pixel ( face , in_face , px , degree ) ;
    }
  }

} ;

// This functor 'looks at' a full spherical image. It receives
// 2D lat/lon coordinates and yields pixel values. Pick-up is done
// with bilinear interpolation. This is needed when we generate a
// cubemap from a lat/lon image. It's the equivalent of what OIIO
// does with it's 'environment' function. While OIIO's 'environemnt'
// receives derivatives and can cover a wide range of scale changes,
// this functor is only suitable for moderate up-scaling, and when
// down-scaling, 'twining' should be used to avoid aliasing.

template < std::size_t nchannels >
struct eval_latlon
: public zimt::unary_functor
           < v2_t , zimt::xel_t < float , nchannels > , LANES >
{
  typedef zimt::xel_t < float , nchannels > px_t ;
  typedef zimt::simdized_type < px_t , 16 > px_v ;

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

    // this can go later:

    // assert ( all_of ( uli[0] >= 0 ) ) ;
    // assert ( all_of ( uri[0] >= 0 ) ) ;
    // assert ( all_of ( lli[0] >= 0 ) ) ;
    // assert ( all_of ( lri[0] >= 0 ) ) ;
    // 
    // assert ( all_of ( uli[0] < w ) ) ;
    // assert ( all_of ( uri[0] < w ) ) ;
    // assert ( all_of ( lli[0] < w ) ) ;
    // assert ( all_of ( lri[0] < w ) ) ;
    // 
    // assert ( all_of ( uli[1] >= 0 ) ) ;
    // assert ( all_of ( uri[1] >= 0 ) ) ;
    // assert ( all_of ( lli[1] >= 0 ) ) ;
    // assert ( all_of ( lri[1] >= 0 ) ) ;
    // 
    // assert ( all_of ( uli[1] < h ) ) ;
    // assert ( all_of ( uri[1] < h ) ) ;
    // assert ( all_of ( lli[1] < h ) ) ;
    // assert ( all_of ( lri[1] < h ) ) ;

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

// next we have a functor converting discrete cubemap coordinates
// (so, pixel units, starting at (0,0) for the upper left corner)
// to pixel values gleaned from an openEXR lat/lon environment
// map, using OIIO's 'environment' function. The environment is
// introduced via it's texture handle th referring to it's internal
// rfpresentation in the texture system ts. The calling code can
// set the batch options to influence the rendition - with the
// default options, OIIO uses a quite intense antialiasing filter
// which removes high frequency content, giving the result a
// slightly blurred appearance.

template < std::size_t nchannels >
struct eval_env
: public zimt::unary_functor
   < v2i_t , zimt::xel_t < float , nchannels > , LANES >
{
  TextureSystem * ts ;
  TextureOptBatch & batch_options ;
  TextureSystem::TextureHandle * th ;
  int width ;
  double px2_to_model ;
  const sixfold_t < nchannels > & sf ;
  ir_to_exr < nchannels > to_exr ;

  // pull in the c'tor arguments

  eval_env ( TextureSystem * _ts ,
             TextureOptBatch & _batch_options ,
             TextureSystem::TextureHandle * _th ,
             const sixfold_t < nchannels > & _sf )
  : ts ( _ts ) ,
    batch_options ( _batch_options ) ,
    th ( _th ) ,
    width ( int ( _sf.section_px ) ) ,
    px2_to_model ( _sf.px_to_model * 0.5 ) ,
    sf ( _sf ) ,
    to_exr ( _sf.section_md , _sf.refc_md )
   { }

  // set up the eval function.

  template < typename I , typename O >
  void eval ( const I & crd2 , O & px ) const
  {
    // Incoming, we have 2D discrete (!) coordinates pertaining
    // to the IR image. We convert them to 3D discrete coordinates
    // with doubled (!) value at appropriate distance. This way,
    // we can process incoming discrete coordinates as integer
    // values and only move to float when we need to. This may
    // be considered 'playful' - it's more to demonstrate how
    // zimt::transform without a source argument feeds discrete
    // target coodinates.

    v3i_v crdi3 { 2 * crd2[0] ,
                  2 * crd2[1] ,
                  i_v ( 2 * width ) } ;

    // to get the right and lower neighbour, we add two (!) to the
    // appropriate component (we're working in doubled coordinates!)

    auto crdi3_x1 = crdi3 ;
    crdi3_x1[0] += 2 ;

    auto crdi3_y1 = crdi3 ;
    crdi3_y1[1] += 2 ;

    // now we move to floating point and model space units, note
    // the factor px2_to_model, which also takes care of halving
    // the doubled coordinates

    crd3_v p00 ( crdi3 ) ;
    p00 *= px2_to_model ;

    crd3_v p10 ( crdi3_x1 ) ;
    p10 *= px2_to_model ;

    crd3_v p01 ( crdi3_y1 ) ;
    p01 *= px2_to_model ;

    // now we obtain ray coordinates from model space coordinates.
    // We have a ready-made functor for the purpose already set up:
    
    crd3_v p00r , p10r , p01r ;
    
    to_exr.eval ( p00 , p00r ) ;
    to_exr.eval ( p10 , p10r ) ;
    to_exr.eval ( p01 , p01r ) ;

    // with the ray coordinates for the current coordinate and it's
    // two neighbours, we can obtain a reasonable approximation of
    // the derivatives in canonical x and y direction by forming the
    // difference

    auto ds = p10r - p00r ;
    auto dt = p01r - p00r ;

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

    p00r.store ( scratch ) ;
    ds.store ( scratch + nchannels * LANES ) ;
    dt.store ( scratch + nchannels * LANES ) ;

    ts->environment ( th , nullptr, batch_options ,
                      Tex::RunMaskOn ,
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
                      Tex::RunMaskOn ,
                      p00r[0].data() ,
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

// cubemap_to_latlon converts a cubemap to a lat/lon environment
// map, a.k.a full spherical panorama.

template < std::size_t nchannels >
void cubemap_to_latlon ( const std::string & input ,
                         std::size_t height ,
                         const std::string & latlon ,
                         int degree )
{
  std::size_t w ;
  std::size_t h ;
  int xres ;
  std::vector < std::string > filename6 ;

  auto inp = ImageInput::open ( input ) ;

  if ( inp == nullptr && store6 )
  {
    // no such image. if store6 is set, we may have a cubemap made
    // up of separate images

    six_names ( input , filename6 ) ;
    inp = ImageInput::open ( filename6[0] ) ;
  }

  assert ( inp ) ;
  const ImageSpec &spec = inp->spec() ;
  xres = spec.width ;
  w = spec.width ;
  h = spec.height ;

  if ( verbose )
    std::cout << "input cube face width: " << w << std::endl ;

  // we set up the sixfold_t object, 'preparing the ground'
  // to pull in the image data

  sixfold_t<nchannels> sf ( w ,
                            support_min_px ,
                            tile_px ,
                            face_fov ) ;

  if ( filename6.size() == 6 )
  {
    assert ( w == h ) ;

    // load the cube faces into slots in the array in the
    // sixfold_t object

    sf.load ( filename6 ) ;
  }
  else
  {
    assert ( w * 6 == h ) ;

    if ( verbose )
      std::cout << "load cubemap image " << input << std::endl ;

    sf.load ( inp ) ;
    inp->close() ;
  }

  // now we build up the 'support'. The square images for the
  // cube faces are cut off at the ninety-degree point, and to
  // interpolate in the vicinity of a cube face's margin, we
  // need some support in good quality. Further out from the
  // margins, the quality of the support data isn't so critical -
  // these pixels only come to play when the cube face is
  // mip-mapped. But going over these pixels again and filling
  // them in with values from bilinear interpolation over the
  // data in the neighbouring cube faces does no harm - we might
  // reduce the target area to some thinner stripe near the
  // edge, but since this is only working over a small-ish part
  // of the data, for now, we'll just redo the lot a few times.

  // we start out mirroring out the square's 1-pixel edge if
  // there is no 'inherent support', which we have only if the
  // cube face image spans more than ninety degrees. Note that
  // the possible contribution of the pixels provided by
  // mirroring is very small, and even quite negligible with
  // 'normal' face widths. But if the value at that position
  // weren't initialized and it would enter the interpolation,
  // there's no telling what might happen - if it were very
  // large, even a small contribution would spoil the result.
  // Hence the mirroring-around. With it, we can rest assured
  // that the result will be near-perfect, and with the
  // subsequent support-gleaning runs we're well on the safe
  // side. A sloppy approach would be to initialize the area
  // with, say, medium grey.

  if ( sf.inherent_support_px == 0 )
    sf.mirror_around() ;

  // to refine the result, we generate support with bilinear
  // interpolation. We assume that the face width will not be
  // 'very small' and the errors from using the mirrored pixels
  // will be negligible.

  sf.fill_support ( 1 ) ;

  // to load the texture with OIIO's texture system code, it
  // has to be in a file, so we store it to a temporary file:

  auto temp_path = std::filesystem::temp_directory_path() ;
  auto temp_filename = temp_path / "temp_texture.exr" ;

  if ( degree == -1 )
  {
    if ( verbose )
      std::cout << "saving generated texture to " << temp_filename.c_str()
                << std::endl ;

    save_array ( temp_filename.c_str() , sf.store ) ;

    // now we can introduce it to the texture system and receive
    // a texture handle for fast access

    TextureSystem * ts = TextureSystem::create() ;
    ustring uenvironment ( temp_filename.c_str() ) ;
    auto th = ts->get_texture_handle ( uenvironment ) ;

    sf.set_ts ( ts , th ) ;

    if ( ts_options != std::string() )
    {
      if ( verbose )
        std::cout << "adding texture system options: " << ts_options
                  << std::endl ;

      ts->attribute ( "options" , ts_options ) ;
    }
  }

  // with proper support near the edges, we can now run the
  // actual payload code - the conversion to lon/lat - with
  // bilinear interpolation, which is a step up from the
  // nearest-neighbour interpolation we had until now and
  // is good enough if the resolution of input and output
  // don't differ too much. But the substrate we have
  // generated is also good for mip-mapping, so we have the
  // option to switch to interpolation methods with larger
  // support, like OIIO's texture access.
  // The internal representation - as it is now - is also
  // useful for other interpolation schemes which need
  // support - if the frame is wide enough, running a b-spline
  // prefilter is quite feasible, because the disturbances due
  // to the imperfections further out from the areas of interest
  // will be negligible.
  // Note also that with the conversion of an incoming 2D or 3D
  // coordinate to a single 2D coordinate pertaining to the
  // internal representation (rather than the pair of cube face
  // index and in-face coordinate) picking up the result is also
  // more efficient and may outperform the current method used
  // in lux, making cubemaps even better suited as environments.
  // It might be a good idea to store the entire internal
  // representation - cubefaces plus support - to an image file
  // for faster access, avoiding the production of the support
  // area from outher cube faces, at the cost of saving a few
  // extra pixels. The process as it stands now seems to produce
  // support which is just as good as support which is created
  // by rendering cube face images with slightly more than ninety
  // degrees (as it can be done in lux), so we can process the
  // 'orthodox' format and yet avoid it's shortcomings. For now
  // I'll stick to using the the cubemap format with cube face
  // images covering precisely ninety degrees.

  // time to do the remaining work.

  // The target array will receive the output pixel data

  typedef zimt::xel_t < float , nchannels > px_t ;

  zimt::array_t < 2, px_t > trg ( { 2 * height , height } ) ;

  // we directly provide input in 'model space units' with a
  // zimt linspace_t object, rather than feeding discrete coordinates
  // and scaling and shifting them, which would work just the same.

  // set up a linspace_t over the lon/lat sample points as get_t
  // (a.k.a input generator). d is the step width from one sample
  // to the next:

  double d = M_PI / height ;

  // note that the first sample point is not at -pi, -pi/2 but
  // half a sample step in the horizonal and vertical toward the
  // center. This way, the sample points are distributed evenly
  // around the image center - the last sample point will be
  // at pi-d/2, pi/2-d/2.
  
  v2_t start { - M_PI + d / 2.0 , - M_PI_2 + d / 2.0 } ;
  v2_t step { d , d } ;

  zimt::linspace_t < float , 2 , 2 , LANES > linspace ( start , step ) ;

  // set up a zimt::storer writing to to the target array

  zimt::storer < float , nchannels , 2 , LANES > st ( trg ) ;

  TextureOptBatch batch_options ;

  // TextureOptBatch's c'tor does not initialize these members, hence:

  for ( int i = 0 ; i < 16 ; i++ )
    batch_options.swidth[i] = batch_options.twidth[i] = 1 ;
  
  for ( int i = 0 ; i < 16 ; i++ )
    batch_options.sblur[i] = batch_options.tblur[i] = 0 ;

  // set up ll_to_px_t using bilinear interpolation (argument 1)
  // or argument -1 to use pick-up with OIIO's 'texture' function.
  // The value is taken from the global 'itp' which is set with
  // the CL parameter of the same name.
  // This is the functor which takes lon/lat coordinates and produces
  // pixel data from the internal representation of the cubemap held
  // in the sixfold_t object.

  ll_to_px_t<nchannels> act ( sf , degree , step , batch_options ) ;

  if ( twine > 1 )
  {
    // with twine > 1, we use a 'twine_t' object wrapping the
    // act functor which we have assembled so far. This has the
    // effect of picking up several values in the close vicinity
    // of the given pick-up point and averaging them. This is the
    // same, mathematically, as oversampling and applying a box
    // filter to the resulting oversampled signal.

    if ( verbose )
      std::cout << "applying a twine of " << twine << std::endl ;

    std::vector < zimt::xel_t < float , 3 > > twine_v ;
    make_spread ( twine , twine ,
                  sf.px_to_model * twine_px ,
                    twine_sigma , twine_threshold , twine_v ) ;

    twine_t < nchannels > twined_act ( act , twine_v ) ;

    zimt::process ( trg.shape , linspace , twined_act , st ) ;
  }
  else
  {
    zimt::process ( trg.shape , linspace , act , st ) ;
  }

  // we don't need the temporary texture any more. If save_ir is
  // set, we save it under the given name.
  if ( save_ir != std::string() )
  {
    if ( verbose )
      std::cout << "saving internal representation to '" << save_ir
                << "'" << std::endl ;
    if ( degree == -1 )
      std::filesystem::rename ( temp_filename , save_ir ) ;
    else
      save_array ( save_ir , sf.store ) ;
  }
  else if ( degree == -1 )
  {
    if ( verbose )
      std::cout << "removing temporary texture file "
                << temp_filename.c_str() << std::endl ;
    std::filesystem::remove ( temp_filename ) ;
  }

  // finally we store the data to an image file - note how we have
  // float data in 'trg', and OIIO will convert these on-the-fly to
  // HALF, as specified in the write_image invocation.
  // Note that the target will receive image data in the same colour
  // space as the input. If you feed, e.g. openEXR, and store to JPEG,
  // the image will look too dark, because the linear RGB data are
  // stored as if they were sRGB.
  // the output is a lat/lon environment with openEXR image order and
  // orientation, so we set the textureformat tag, but TODO: there may
  // be slight differences in the precise format of the individual
  // cube faces - the cube faces we store here are precisely ninety
  // degrees from the left edge of the leftmost pixel to the right
  // edge of the rightmost pixel (and alike for top/bottom), whereas
  // the openEXR spec may measure the ninety degrees from the pixel
  // centers, so the output would not conform to their spec precisely

  if ( verbose )
    std::cout << "saving lat/lon environment map to '"
              << latlon << std::endl ;

  save_array<nchannels> ( latlon , trg , true ) ;
}

template < std::size_t nchannels >
using pix_t = zimt::xel_t < float , nchannels > ;

// this function takes a lat/lon image as it's input and transforms
// it into the IR image of a sixfold_t. The source image is passed
// in as a zimt::view_t, the target as a reference to sixfold_t.
// This function is only used if the cubemap is 'regenerated' from
// the lat/lon image we have just created, or if the cubemap is
// made from a lat/lon image passed in as a file.

template < std::size_t nchannels >
void latlon_to_ir ( const zimt::view_t < 2 , pix_t<nchannels> > & latlon ,
                    sixfold_t < nchannels > & sf )
{
  // we set up a linspace_t object to step through the sample points.
  // in image coordinates, we'd use these start and step values:

  v2_t start { 0.5 , 0.5 } ;

  v2_t step { 1.0 , 1.0 } ;

  // but we want to feed model space coordinates to the 'act'
  // functor, hence we scale:

  start *= sf.px_to_model ;
  step *= sf.px_to_model ;

  zimt::linspace_t < float , 2 , 2 , LANES > ls ( start , step ) ;

  // the data are to be stored to the IR image, held in sf.store

  zimt::storer < float , nchannels , 2 , LANES > st ( sf.store ) ;

  // the act functor is a chain of three separate functors, which we
  // set up first:

  // this one converts coordinates pertaining to the IR image to
  // 3D ray coordinates.

  ir_to_ray < nchannels > itr
    ( sf.section_md , sf.refc_md ) ;

  // this one converts 3D ray coordinates to lat/lon values

  ray_to_ll_t rtl ;

  // and this one does the pick-up from the lat/lon image

  eval_latlon < nchannels > ltp ( latlon ) ;

  // now we form the act functor by chaining these three functors

  auto act = itr + rtl + ltp ;

  if ( twine > 1 )
  {
    // with twine > 1, we use a 'twine_t' object wrapping the
    // act functor which we have assembled so far. This has the
    // effect of picking up several values in the close vicinity
    // of the given pick-up point and averaging them. This is the
    // same, mathematically, as oversampling and applying a box
    // filter to the resulting oversampled signal.

    if ( verbose )
      std::cout << "applying a twine of " << twine << std::endl ;

    std::vector < zimt::xel_t < float , 3 > > twine_v ;
    make_spread ( twine , twine ,
                  sf.px_to_model * twine_px ,
                    twine_sigma , twine_threshold , twine_v ) ;

    twine_t < nchannels > twined_act ( act , twine_v ) ;

    zimt::process ( sf.store.shape , ls , twined_act , st ) ;
  }
  else
  {
    // no twine given, use act as it is

    zimt::process ( sf.store.shape , ls , act , st ) ;
  }
}

// convert a lat/lon environment map into a cubemap. The 'width'
// parameter sets the width of the cubemap's cube face images, it's
// height is six times that. The 'degree' parameter selects which
// interpolator is used: degree 1 uses direct bilinear interpolation,
// and degree -1 uses OIIO's anisotropic antialiasing filter.

template < std::size_t nchannels >
void latlon_to_cubemap ( const std::string & latlon ,
                         std::size_t width ,
                         const std::string & cubemap ,
                         int degree )
{
  sixfold_t<nchannels> sf ( width ,
                            support_min_px ,
                            tile_px ,
                            face_fov ) ;

  if ( degree == 1 )
  {
    typedef zimt::xel_t < float , nchannels > px_t ;

    auto inp = ImageInput::open ( latlon ) ;
    assert ( inp != nullptr ) ;

    const ImageSpec &spec = inp->spec() ;
    std::size_t w = spec.width ;
    std::size_t h = spec.height ;
    assert ( w == h * 2 ) ;

    zimt::array_t < 2 , px_t > src ( { w , h } ) ;
    
    bool success = inp->read_image ( 0 , 0 , 0 , nchannels ,
                                    TypeDesc::FLOAT , src.data() ) ;

    assert ( success ) ;

    latlon_to_ir < nchannels > ( src , sf ) ;
  }
  else
  {
    TextureSystem * ts = TextureSystem::create() ;
    if ( ts_options != std::string() )
    {
      if ( verbose )
        std::cout << "adding texture system options: " << ts_options
                  << std::endl ;

      ts->attribute ( "options" , ts_options ) ;
    }
    ustring uenvironment ( latlon.c_str() ) ;
    auto th = ts->get_texture_handle ( uenvironment ) ;
    TextureOptBatch batch_options ;

    // TextureOptBatch's c'tor does not initialize these members, hence:

    for ( int i = 0 ; i < 16 ; i++ )
      batch_options.swidth[i] = batch_options.twidth[i] = 1 ;
    
    for ( int i = 0 ; i < 16 ; i++ )
      batch_options.sblur[i] = batch_options.tblur[i] = 0 ;

    // we have a dedicated functor going all the way from discrete
    // cubemap image coordinates to pixel values:

    assert ( ts != nullptr ) ;
    assert ( th != nullptr ) ;

    eval_env < nchannels > act ( ts , batch_options , th , sf ) ;

    // showtime! notice the call signature: we pass no source, because
    // we want discrete coordinates to start from. Omitting the source
    // parameter does just that: the input to the act functor will be
    // discrete coordinates of the target location for which the functor
    // is supposed to calculate content.

    if ( twine > 1 )
    {
      if ( verbose )
        std::cout << "applying a twine of " << twine << std::endl ;

      std::vector < zimt::xel_t < float , 3 > > twine_v ;
      make_spread ( twine , twine , twine_px ,
                    twine_sigma , twine_threshold , twine_v ) ;

      twine_t < nchannels > twined_act ( act , twine_v ) ;

      zimt::transform ( twined_act , sf.store ) ;
    }
    else
    {
      zimt::transform ( act , sf.store ) ;
    }
  }

  // the result of the zimt::transform is slightly more than we need,
  // namely the entire sixfold_t object's IR image. We only store the
  // 'cubemap proper':

  if ( verbose )
    std::cout << "saving cubemap to '" << cubemap << "'" << std::endl ;

  if ( store6 )
  {
    std::vector < std::string > filename6 ;
    six_names ( cubemap , filename6 ) ;
    sf.store_cubemap ( filename6 ) ;
    if ( store6_lux )
    {
      // write a .lux file
      auto point_pos = cubemap.find_last_of ( "." ) ;
      assert ( point_pos != std::string::npos ) ;
      std::string base = cubemap.substr ( 0 , point_pos ) ;
      std::string ext ( ".lux" ) ;
      std::ofstream ofs ( base + ext ) ;
      ofs << "# lux cubemap script file made with envutil" << std::endl ;
      ofs << "projection=cubemap" << std::endl ;
      ofs << "cubeface_fov=" << 180.0 / M_PI * sf.face_fov << std::endl ;
      ofs << "cube_right=" << filename6[0] << std::endl ;
      ofs << "cube_left=" << filename6[1] << std::endl ;
      ofs << "cube_top=" << filename6[2] << std::endl ;
      ofs << "cube_bottom=" << filename6[3] << std::endl ;
      ofs << "cube_back=" << filename6[4] << std::endl ;
      ofs << "cube_front=" << filename6[5] << std::endl ;
      ofs.close() ;
    }
  }
  else
  {
    sf.store_cubemap ( cubemap ) ;
  }

  if ( save_ir != std::string() )
  {
    if ( verbose )
      std::cout << "saving internal representation to '" << save_ir
                << "'" << std::endl ;

    save_array ( save_ir , sf.store ) ;
  }
}

int main ( int argc , const char ** argv )
{
  // we're using OIIO's argparse, since we're using OIIO anyway.
  // This is a convenient way to glean arguments on all supported
  // platforms - getopt isn't available everywhere.

  Filesystem::convert_native_arguments(argc, (const char**)argv);
  ArgParse ap;
  ap.intro("envutil -- convert between lat/lon and cubemap format\n")
    .usage("envutil [options] --input INPUT --output OUTPUT");
  ap.arg("-v", &verbose)
    .help("Verbose output");
  ap.arg("--input INPUT")
    .help("input file name (mandatory)")
    .metavar("INPUT");
  ap.arg("--output OUTPUT")
    .help("output file name (mandatory)")
    .metavar("OUTPUT");
  ap.arg("--save_ir INTERNAL")
    .help("save IR image to this file")
    .metavar("INTERNAL");
  ap.arg("--ts_options OPTIONS")
    .help("pass comma-separates k=v list of options to OIIO's texture system")
    .metavar("OPTIONS");
  ap.arg("--extent EXTENT")
    .help("width of the cubemap / height of the envmap")
    .metavar("EXTENT");
  ap.arg("--itp ITP")
    .help("interpolator: 1 for direct bilinear, -1 for OIIO's anisotropic")
    .metavar("ITP");
  ap.arg("--twine TWINE")
    .help("use twine*twine oversampling and box filter - best with itp1")
    .metavar("TWINE");
  ap.arg("--twine_px TWINE_WIDTH")
    .help("widen the pick-up area of the twining filter")
    .metavar("TWINE_WIDTH");
  ap.arg("--twine_sigma TWINE_SIGMA")
    .help("use a truncated gaussian for the twining filter (default: don't)")
    .metavar("TWINE_SIGMA");
  ap.arg("--twine_threshold TWINE_THRESHOLD")
    .help("discard twining filter taps below this threshold")
    .metavar("TWINE_THRESHOLD");
  ap.arg("--face_fov FOV")
    .help("field of view of the cube faces of a cubemap input (in degrees)")
    .metavar("FOV");
  ap.arg("--support_min_px EXTENT")
    .help("minimal additional support around the cube face proper")
    .metavar("EXTENT");
  ap.arg("--tile_px EXTENT")
    .help("tile width for the internal representation image")
    .metavar("EXTENT");
  ap.arg("--ctc", &ctc)
    .help("flag indicating fov is measured between marginal pixel centers");
  ap.arg("--6", &store6)
    .help("use six separate cube face images");
  ap.arg("--lux", &store6_lux)
    .help("use cube face images named in lux convention");
    
  if (ap.parse(argc, argv) < 0 ) {
      std::cerr << ap.geterror() << std::endl;
      ap.print_help();
      return help ? EXIT_SUCCESS : EXIT_FAILURE ;
  }

  if (!metamatch.empty()) {
      field_re.assign(metamatch, std::regex_constants::extended
                                      | std::regex_constants::icase);
  }
  
  // extract the CL arguments from the argument parser

  input = ap["input"].as_string("");
  output = ap["output"].as_string("");
  ts_options = ap["ts_options"].as_string("");
  save_ir = ap["save_ir"].as_string("");
  extent = ap["extent"].get<int>(0);
  itp = ap["itp"].get<int>(-1);
  twine = ap["twine"].get<int>(1);
  support_min_px = ap["support_min_px"].get<int>(4);
  tile_px = ap["tile_px"].get<int>(64);
  twine_px = ap["twine_px"].get<float>(1.0);
  twine_sigma = ap["twine_sigma"].get<float>(0.0);
  twine_threshold = ap["twine_threshold"].get<float>(0.0);
  face_fov = ap["face_fov"].get<float>(90.0);
  face_fov *= M_PI / 180.0 ;

  assert ( input != std::string() ) ;
  assert ( output != std::string() ) ;

  if ( verbose )
  {    
    std::cout << "input: " << input << std::endl ;
    std::cout << "output: " << output << std::endl ;
    std::cout << "interpolation: "
              << ( itp == 1 ? "direct bilinear" : "OIIO anisotropic" )
              << std::endl ;
    std::cout << "cube face fov: " << face_fov << std::endl ;
  }

  auto inp = ImageInput::open ( input ) ;

  std::size_t w ;
  std::size_t h ;
  std::size_t nchannels ;
  std::vector < std::string > filename6 ;

  if ( inp == nullptr && store6 )
  {
    // no such image. if store6 is set, we may have a cubemap made
    // up of separate images

    six_names ( input , filename6 ) ;
    inp = ImageInput::open ( filename6[0] ) ;
  }
  assert ( inp ) ;

  const ImageSpec &spec = inp->spec() ;
  w = spec.width ;
  h = spec.height ;
  nchannels = spec.nchannels ;

  inp->close() ;

  if ( verbose )
    std::cout << "input has " << nchannels << " channels" << std::endl ;

  if ( w == 2 * h )
  {
    if ( verbose )
      std::cout << "input has 2:1 aspect ratio, assuming latlon"
                << std::endl ;

    if ( extent == 0 )
    {
      double e = h * 2.0 / M_PI ;
      extent = e ;
      if ( extent % 64 )
        extent = ( ( extent / 64 ) + 1 ) * 64 ;
      if ( verbose )
        std::cout << "no extent given, using " << extent << std::endl ;
    }

    // the 'ctc' flag indicates that the field of view of the
    // cube face images is given referring to the angle
    // between the leftmost and rightmost pixel's center, as
    // opposed to edge-to-edge, which is standard. Some cubemaps
    // are said to be made this way - if the face fov is given
    // as ninety degrees, we'll calculate with the adjusted
    // slightly larger value here to come out right.

    if ( ctc )
    {
      double half_md = tan ( face_fov / 2.0 ) ;
      half_md *= ( ( extent + 1.0 ) / extent ) ;
      face_fov = atan ( half_md ) * 2.0 ;
      if ( verbose )
        std::cout << "ctc is set, adjusted face_fov to "
                  << face_fov << std::endl ;
    }

    if ( nchannels >= 4 )
    {
      latlon_to_cubemap<4> ( input , extent , output , itp ) ;
    }
    else if ( nchannels == 3 )
    {
      latlon_to_cubemap<3> ( input , extent , output , itp ) ;
    }
    else if ( nchannels == 1 )
    {
      latlon_to_cubemap<1> ( input , extent , output , itp ) ;
    }
    else
    {
      std::cerr << "input format error: need 1,3 or >=4 channels"
                << std::endl ;
      exit ( EXIT_FAILURE ) ;
    }
  }
  else if ( h == 6 * w || filename6.size() == 6 )
  {
    if ( verbose )
    {
      if ( h == 6 * w )
        std::cout << "input has 1:6 aspect ratio, assuming cubemap"
                << std::endl ;
      else
        std::cout << "expecting six separate cube face images"
                  << std::endl ;
    }

    // the 'ctc' flag indicates that the field of view of the
    // cube face images is given referring to the angle
    // between the leftmost and rightmost pixel's center, as
    // opposed to edge-to-edge, which is standard. Some cubemaps
    // are said to be made this way - if the face fov is given
    // as ninety degrees, we'll calculate with the adjusted
    // slightly larger value here to come out right.

    if ( ctc )
    {
      double half_md = tan ( face_fov / 2.0 ) ;
      half_md *= ( ( w + 1.0 ) / w ) ;
      face_fov = atan ( half_md ) * 2.0 ;
      if ( verbose )
        std::cout << "ctc is set, adjusted face_fov to "
                  << face_fov << std::endl ;
    }

    if ( extent == 0 )
    {
      extent = 4 * w ;
      if ( extent % 64 )
        extent = ( ( extent / 64 ) + 1 ) * 64 ;
      if ( verbose )
        std::cout << "no extent given, using " << extent << std::endl ;
    }

    if ( nchannels >= 4 )
    {
      cubemap_to_latlon<4> ( input , extent , output , itp ) ;
    }
    else if ( nchannels == 3 )
    {
      cubemap_to_latlon<3> ( input , extent , output , itp ) ;
    }
    else if ( nchannels == 1 )
    {
      cubemap_to_latlon<1> ( input , extent , output , itp ) ;
    }
    else
    {
      std::cerr << "input format error: need 1,3 or >=4 channels"
                << std::endl ;
      exit ( EXIT_FAILURE ) ;
    }
  }
  else
  {
    std::cerr << "input format error: need lat/lon or cubemap input" << std::endl ;
    exit ( EXIT_FAILURE ) ;
  }
  if ( verbose )
    std::cout << "conversion complete. exiting." << std::endl ;
  OIIO::geterror() ;
  OIIO::shutdown() ;
  exit ( EXIT_SUCCESS ) ;
}

