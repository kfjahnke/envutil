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

// this header has basic enums and types which do not depend on zimt,
// and some helper functions which don't use zimt.

#ifndef ENVUTIL_BASIC_H
#define ENVUTIL_BASIC_H

#include <cmath>
#include <iostream>

// enum encoding the sequence of cube face images in the cubemap
// This is the sequence used for openEXR cubmap layout. The top
// and bottom squares are oriented so as to align with the back
// image. Of course, the labels are debatable: my understanding
// of 'front' is 'aligned with the image center'. If one were to
// associate 'front' with the wrap-around point of the full
// spherical, the labels would be different.

typedef enum
{
  CM_LEFT ,
  CM_RIGHT ,
  CM_TOP ,
  CM_BOTTOM ,
  CM_FRONT ,
  CM_BACK
} face_index_t ;

// we use lux coordinate system convention. I call it 'latin book
// order': if you have a stack of prints in front of you and read
// them, your eyes move first left to right inside the line, then
// top to bottom from line to line, then, moving to the next pages,
// forward in the stack. Using this order also makes the first two
// components agree with normal image indexing conventions, namely
// x is to the right and y down. Note that I put the fastest-moving
// index first, which is 'fortran' style, whereas C/C++ use the
// opposite order for nD arrays.

enum { RIGHT , DOWN , FORWARD } ;

// openEXR uses different 3D axis semantics, and if we want to use
// OIIO's environment lookup function, we need openEXR 3D coordinates.

// Here's what the openEXR documentation sys about their axis
// order (next to a drawing which says differently, see this issue:
// https://github.com/AcademySoftwareFoundation/openexr/issues/1687)

// quote:
// We assume that a camera is located at the origin, O, of a 3D
// camera coordinate system. The camera looks along the positive z
// axis. The positive x and y axes correspond to the camera’s left
// and up directions.
// end quote

// so we'd get this axis order, assuming they store x,y,z:

enum { EXR_LEFT , EXR_UP , EXR_FORWARD } ;

// the cubemap comes out right this way, so I assume that their text
// is correct and the drawing is wrong.

typedef enum
{
  SPHERICAL ,
  CYLINDRICAL ,
  RECTILINEAR ,
  STEREOGRAPHIC ,
  FISHEYE ,
  CUBEMAP ,
  BIATAN6 ,
  PRJ_NONE
}  projection_t ;

const char * const projection_name[]
{
  "spherical" ,
  "cylindrical" ,
  "rectilinear" ,
  "stereographic" ,
  "fisheye" ,
  "cubemap" ,
  "biatan6" ,
  "unsupported"
} ;

extern long rt_cumulated ;

// assuming isotropic sampling (same sampling resolution in the horizontal
// and vertical), calculate the vertical field of view from the horizontal
// field of view, under the given projection.
// Note that this function is for centered images only
// (x0 == -x1), (y0 == -y1)

double get_vfov ( projection_t projection ,
                  int width ,
                  int height ,
                  double hfov ) ;

// the 'step' of an image is the angle - in radians - which
// corresponds to the width of one pixel in the image center.
// for some projections and in certain directions, this value
// will be usable at non-central points (e.g. for spherical
// images along the horizon). In any case it can be used as a
// 'rule of thumb' indicator of the image's resolution.
// If we have the 'extent' of a spherical or cylindrical image
// already, we can calculate the step as hfov / width;
// for other projections this simple formula doesn't apply.
// Note that this function is for centered images only
// (x0 == -x1), (y0 == -y1)
// currently unused.

double get_step ( projection_t projection ,
                  int width ,
                  int height ,
                  double hfov ) ;

// extent_type is a handy struct to contain the horizontal and
// vertical extent of a rectangular section of the plane.

struct extent_type
{
  double x0 , x1 , y0 , y1 ;

  friend std::ostream & operator<<
    ( std::ostream & osr , const extent_type & e )
  {
    osr << "ext { " << e.x0 << " " << e.x1 << " "
        << e.y0 << " " << e.y1 << " }" ;
    return osr ;
  }
} ;

// extract internally uses the notion of an image's 'extent' in 'model
// space'. The image is thought to be 'draped' to an 'archetypal 2D
// manifold' - the surface of a sphere or cylinder with unit radius
// or a plane at unit distance forward - where the sample points are
// placed on the 2D manifold so that rays from the origin to the
// scene point which corresponds with the sample point intersect there.
// To put it differently: the sample point cloud is scaled and shifted
// to come to lie on the 'archetypal' 2D manifolds. This makes for
// efficient calculations. The image is taken to be centered on the
// 'forward' ray.

extent_type get_extent ( projection_t projection ,
                         int width ,
                         int height ,
                         double hfov ) ;

struct cp_t
{
  int t ;
  int n ;
  int N ;
  double x ;
  double y ;
  double X ;
  double Y ;
} ;

#include <regex>
#include <OpenImageIO/imageio.h>
#include <OpenImageIO/imagebuf.h>
#include <OpenImageIO/imagebufalgo.h>

using OIIO::ImageInput ;
using OIIO::ImageOutput ;
using OIIO::TypeDesc ;
using OIIO::ImageSpec ;

// struct image_series holds a format string for a series of numbered
// images. It's a standard printf-type format string containing precisely
// one %-sequence accepting an integer - e.g. %04d
// It has operator[] to provide the filename fixed with the given index.

struct image_series
{
  std::string format_string ;
  std::size_t buffer_size ;

  image_series ( const std::string & _format_string )
  : format_string ( _format_string )
  {
    // scan the format string for percent signs. Only one of them
    // is allowed.

    int percent_count = 0 ;
    int percent_pos = 0 ;
    int pos = 0 ;

    for ( auto const & c : format_string )
    {
      if ( c == '%' )
      {
        percent_count++ ;
        percent_pos = pos ;
      }
      ++pos ;
    }

    if ( percent_count == 1 )
      buffer_size = pos + 16 ;
  }

  bool valid()
  {
    return ( buffer_size != 0 ) ;
  }

  std::string operator[] ( const std::size_t & index )
  {
    if ( buffer_size )
    {
      char buffer [ buffer_size ] ;
      snprintf ( buffer , buffer_size ,
                 format_string.c_str() , int(index) ) ;
      return buffer ;
    }
    else
    {
      return format_string ;
    }
  }
} ;

// similar class, used for six cube faces

struct cubeface_series
{
  std::string format_string ;
  std::vector < std::string > filename ;
  std::size_t buffer_size ;

  cubeface_series() = default ;

  cubeface_series ( const std::string & _format_string )
  : format_string ( _format_string )
  {
    // scan the format string for percent signs. Only one of them
    // is allowed.

    int percent_count = 0 ;
    int percent_pos = 0 ;
    int pos = 0 ;

    for ( auto const & c : format_string )
    {
      if ( c == '%' )
      {
        percent_count++ ;
        percent_pos = pos ;
      }
      ++pos ;
    }

    buffer_size = 0 ;

    if ( percent_count == 1 )
    {
      buffer_size = pos + 16 ;

      char buffer [ buffer_size ] ;
      snprintf ( buffer , buffer_size , format_string.c_str() , "left" ) ;
      filename.push_back ( buffer ) ;
      snprintf ( buffer , buffer_size , format_string.c_str() , "right" ) ;
      filename.push_back ( buffer ) ;
      snprintf ( buffer , buffer_size , format_string.c_str() , "top" ) ;
      filename.push_back ( buffer ) ;
      snprintf ( buffer , buffer_size , format_string.c_str() , "bottom" ) ;
      filename.push_back ( buffer ) ;
      snprintf ( buffer , buffer_size , format_string.c_str() , "front" ) ;
      filename.push_back ( buffer ) ;
      snprintf ( buffer , buffer_size , format_string.c_str() , "back" ) ;
      filename.push_back ( buffer ) ;
    }
  }

  bool valid()
  {
    return filename.size() == 6 ;
  }

  // we provide two operator[] overloads. The first one accesses the
  // cubeface filename by number, the second inserts a given orientation
  // string.

  std::string operator[] ( const std::size_t & index )
  {
    assert ( buffer_size != 0 && index < 6 ) ;
    return filename [ index ] ;
  }

  std::string operator[] ( const std::string & face )
  {
    char buffer [ buffer_size ] ;
    snprintf ( buffer , buffer_size ,
               format_string.c_str() ,
               face.c_str() ) ;
    return buffer ;
  }

  // return a reference to the set of six names

  const std::vector < std::string > & get_filenames()
  {
    return filename ;
  }
} ;

// calculate the cross or dot product of two vectors

template < typename U >
U cross ( const U & x , const U & y )
{
  U result ;
  result[0] = x[1] * y[2] - x[2] * y[1] ;
  result[1] = x[2] * y[0] - x[0] * y[2] ;
  result[2] = x[0] * y[1] - x[1] * y[0] ;
  return result ;
}

template < typename U >
typename U::value_type dot ( const U & x , const U & y )
{
  return ( x * y ) . sum() ;
}

// general calculation of the angles between rays:

template < typename U >
typename U::value_type angle ( const U & a , const U & b )
{
  auto costheta = dot ( a , b ) / ( norm ( a ) * norm ( b ) ) ;
  return acos ( costheta ) ;
}

#include "zimt/array.h"

struct pto_mask_type
{
  int image ;
  int variant ;
  std::string vertex_list ;
  std::vector < float > vx ;
  std::vector < float > vy ;

  pto_mask_type()
  : image ( -1 ) ,
    variant ( -1 )
    { }

  friend std::ostream & operator<< ( std::ostream & osr ,
                                     const pto_mask_type & mask )
  {
    std::cout << "type " << mask.variant
              << " mask for facet " << mask.image
              << std::endl ;
    for ( int i = 0 ; i < mask.vx.size() ; i++ )
      std::cout << "  ( " << mask.vx[i] << ", "
                << mask.vy[i] << ")" << std::endl ;
    return osr ;
  }
} ;

void fill_polygon ( const std::vector<float> & px ,
                    const std::vector<float> & py ,
                    int IMAGE_LEFT , int IMAGE_TOP ,
                    int IMAGE_RIGHT , int IMAGE_BOT ,
                    std::function < void ( int , int ) > fillPixel ) ;

// class facet_base is a common base type for images with PTO attributes.
// It's used both for 'facet' images - single oriented source images
// gleaned directly from the command line or a PTO script - and for the
// target image - currently, class 'args' directly inherits from
// facet_base. The factoring-out to this base class makes handling
// the geometries easier - earlier I started out with discrete member
// variables, but I think this approach is better. The common base
// class is especially handy when it comes to 'single' and 'split'
// jobs, where a facet's geometry is taken over as the target's
// geometry to re-create single images from already-stitched ones.
// facet_base itself inherits from extent_type - a standard way of
// expressing the dimensions of a 2D raster in 'model space'.

struct facet_base
: public extent_type
{
  std::string colour_space ;
  std::string projection_str ;
  projection_t projection ;
  double hfov ;
  double step ;
  double yaw , pitch , roll ;
  int width ;
  int height ;
  int window_width ;
  int window_height ;
  int window_x_offset ;
  int window_y_offset ;

  double tr_x , tr_y , tr_z ;
  double tp_y , tp_p , tp_r ;
  double shear_g , shear_t ;
  double s, a, b, c, d, h, v, cap_radius , r_max ;

  bool has_shift ;
  bool has_lcp ;
  bool has_shear ;
  bool has_2d_tf ;
  bool has_translation ;
} ;

struct facet_spec
: public facet_base
{
  int facet_no ;
  int nchannels ;

  std::string filename ;
  std::string asset_key ;

  // PTO can specify that certein marginal parts of a facet image should
  // be 'blacked out' (like a passepartout - rectangular or circular). In
  // envutil, we realize this by an alpha channel manipulation, and if the
  // source image comes without an alpha channel, we add one. The unwanted
  // parts are set to transparent black.

  bool has_lens_crop ;
  int crop_x0 , crop_x1 , crop_y0 , crop_y1 ;

  // PTO also allows the application of polygonal masks to 'black out'
  // unwanted content. Like the feature above, we use alpha and create the
  // A channel if needed.

  bool has_pto_mask ;
  std::vector < pto_mask_type > pto_mask_v ;

  bool init ( int argc , const char ** argv ) ;

  // this member has nothing to do with the PTO masks above. It's used
  // when producing mask images (see masking.h).

  int masked ;

  float brighten ;

  // this member function inspects the geometry-related parameters
  // and sets the internal state accordingly. It also sets a few flags
  // to make it easier for other code to decide which components are
  // present and which can be left out in processing.

  void process_geometry()
  {
    has_shift = ( h != 0.0 || v != 0.0 ) ;
    has_lcp = ( a != 0.0 || b != 0.0 || c != 0.0 ) ;
    has_shear = ( shear_g != 0.0 || shear_t != 0.0 ) ;
    has_2d_tf = ( has_shift || has_lcp || has_shear ) ;

    has_translation = ( tr_x != 0 || tr_y != 0 || tr_z != 0 ) ;

    // reference radius in PTO is half the extent of the smaller edge

    double dv = fabs ( y1 - y0 ) / 2.0 ;
    double dh = fabs ( x1 - x0 ) / 2.0 ;

    s = ( dh < dv ) ? dh : dv ;

    // so the larger of the two is larger by this factor:

    double aspect = ( dh >= dv ) ? dh / dv : dv / dh ;

    // which gives us r_max (expressed in units of s):
    
    r_max = sqrt ( 1 + aspect * aspect ) ;

    // set d so that the image is not scaled

    d = 1.0 - ( a + b + c ) ;

    // the PTO d and e parameters are in pixels, h and v in unit radii

    double factor = fabs ( x1 - x0 ) / width ;
    h *= factor ;
    v *= factor ;

    auto d1 = x0 * x0 + y0 + y0 ;
    auto d2 = x1 * x1 + y0 + y0 ;
    auto d3 = x0 * x0 + y1 + y1 ;
    auto d4 = x1 * x1 + y1 + y1 ;

    d1 = std::max ( d1 , d2 ) ;
    d1 = std::max ( d1 , d3 ) ;
    d1 = std::max ( d1 , d4 ) ;

    cap_radius = sqrt ( d1 ) ;
  }

  void get_image_metrics ( bool read_hfov = false ,
                           bool read_projection = false )
  {
    // currently building with raw::user_flip set to zero, to load
    // raw images in memory order without EXIF rotation. This only
    // affects raw images.

    ImageSpec config;
    config [ "raw:user_flip" ] = 0 ;
    auto inp = ImageInput::open ( filename , &config ) ;

    if ( ! inp )
    {
      // TODO: refactor

      auto has_percent = filename.find_first_of ( "%" ) ;
      if ( has_percent != std::string::npos )
      {
        OIIO::geterror() ;

        // input must be a set of six cubeface images, that's the
        // only way how we accept a format string.

        cubeface_series cfs ( filename ) ;
        if ( cfs.valid() )
        {
          inp = ImageInput::open ( cfs[0] , &config ) ;
        }
      }
      else
      {
        std::cerr << "failed to open facet image '"
                  << filename << "'" << std::endl ;
        exit ( -1 ) ;
      }
    }

    const ImageSpec &spec = inp->spec() ;

    width = window_width = spec.width ;
    height = window_height = spec.height ;
    window_x_offset = window_y_offset = 0 ;
    nchannels = spec.nchannels ;

    if ( read_hfov )
    {
      std::cout << "***** try read hfov metadatum" << std::endl ;
      // glean projection and hfov from metadata
      float fv ;
      bool ok = spec.getattribute ( "Hfov", OIIO::TypeFloat, &fv ) ;
      if ( ok )
      {
        std::cout << "found hfov in metadata: " << fv << std::endl ;
        hfov = fv ;
      }
      else
      {
        std::cout << "no 'Hfov' metadatum found; assuming 65 degrees"
                  << std::endl ;
        hfov = 65.0 ;
      }
    }
    if ( read_projection )
    {
      std::cout << "***** try read projection metadatum" << std::endl ;
      // glean projection and hfov from metadata
      std::string metadatum ;
      metadatum = spec [ "Projection" ] ;
      if ( metadatum != std::string() )
      {
        projection_str = metadatum ;
        std::cout << "found projection in metadata: "
                  << projection_str << std::endl ;
      }
      else
      {
        std::cout << "no 'Projection' metadatum found; assuming 'rectilinear'"
                  << std::endl ;

        projection_str = "rectilinear" ;
        projection = RECTILINEAR ;
      }
    }

    inp->close() ;
  }
} ;

struct arguments
: public facet_base
{
  bool verbose ;
  std::string working_colour_space ;
  std::string input_colour_space ;
  std::string output ;
  std::string pano ;
  std::string split ;
  std::vector < std::string > oiio_option_v ;
  std::size_t support_min ;
  std::size_t tile_size ;
  std::string pto_file ;
  std::string synopsis ;
  int prefilter_degree ;
  int spline_degree ;
  int twine  ;
  std::string twf_file ;
  bool twine_normalize ;
  bool twine_precise ;
  double twine_width , twine_density , twine_sigma , twine_threshold ;
  std::vector < zimt::xel_t < float , 3 > > twine_spread ;
  std::vector < cp_t > cp_v ;

  // gleaned from other parameters or input images

  std::size_t nchannels ;

  // technical variables for the argument parser

  std::string metamatch ;
  std::regex field_re ;

  std::size_t nfacets ;
  std::vector < std::string > facet_name_v ;
  std::vector < std::string > facet_projection_v ;
  std::vector < std::string > facet_hfov_v ;
  std::vector < std::string > facet_yaw_v ;
  std::vector < std::string > facet_pitch_v ;
  std::vector < std::string > facet_roll_v ;
  std::vector < std::string > photo_name_v ;
  std::vector < std::string > addenda ;

  std::vector < facet_spec > facet_spec_v ;
  std::vector < pto_mask_type > pto_mask_v ;

  // if output cropping is specified in the p-line, this flag is
  // set, and the four values defining the cropping window are in
  // p_crop_...

  bool store_cropped ;
  int p_crop_x0 , p_crop_x1 , p_crop_y0 , p_crop_y1 ;

  int solo ;
  int single ;
  int mask_for ;

  // the 'arguments' object's 'init' takes the main program's argc
  // and argv.

  void init ( int argc , const char ** argv ) ;

  // twine_setup initializes the twining coefficients

  void twine_setup() ;
} ;

extern arguments args ;

template < class VT , class NT = float >
VT RGB2sRGB ( VT value )
{
  VT result = NT(1.055) * pow ( value , NT(0.41666666666666667) ) - 0.055 ;

  // using vspline::assign_if to code uniformly for scalar and vectorized VT:

  if ( value <= NT(0.0031308) )
    result = NT(12.92) * value ;

  return result ;
}

// helper function to save a zimt array of pixels to an image file, or
// to a set of six cube face images, if 'output' has a format string.

template < std::size_t nchannels >
void save_array ( const std::string & filename ,
                  zimt::view_t
                    < 2 ,
                      zimt::xel_t < float , nchannels >
                    > pixels ,
                  bool is_latlon = false )
{
  if ( args.projection == CUBEMAP || args.projection == BIATAN6 )
  {
    auto original_prj = args.projection ;

    // output is a cubemap, let's see if 'output' is a format
    // string for six separate cube faces

    assert ( pixels.shape[1] == 6 * pixels.shape[0] ) ;

    auto has_percent = filename.find_first_of ( "%" ) ;
    if ( has_percent != std::string::npos )
    {
      // input must be a set of six cubeface images, that's the
      // only way how we accept a format string.

      cubeface_series cfs ( filename ) ;
      if ( cfs.valid() )
      {
        // we'll call save_array recursively for the single images, hence:

        args.projection = RECTILINEAR ;

        // save six subarrays to individual images

        std::size_t w = pixels.shape[0] ;
        for ( std::size_t i = 0 ; i < 6 ; i++ )
        {
          save_array ( cfs[i] ,
                       pixels.window ( { 0L , long ( i * w ) } ,
                                       { long ( w ) , long ( ( i + 1 ) * w ) } ) ) ;
        }

        // restore the projection in 'args'

        args.projection = original_prj ;

        // we're done.

        return ;
      }
    }
  }

  // if we land here, we're supposed to store an ordinary single image

  ImageSpec ospec ( pixels.shape[0] , pixels.shape[1] ,
                    nchannels , TypeDesc::FLOAT ) ;

  ospec["ImageDescription"] = "image processed by envutil";
  ospec["Projection"] = projection_name [ args.projection ] ;
  ospec["Hfov"] = ( 180.0 / M_PI ) * args.hfov ;

  // we wrap the pixel data in on OIIO ImageBuf. OIIO uses byte
  // strides, so we need to scale up the pixel strides.

  static const size_t px_bytes = sizeof ( float ) * nchannels ;

  OIIO::ImageBuf out_buf ( ospec ,
                           pixels.data() ,
                           pixels.strides[0] * px_bytes ,
                           pixels.strides[1] * px_bytes ) ;

  auto sz = filename.size() ;
  assert ( sz > 4 ) ;
  auto extension = filename.substr ( sz - 4 , sz ) ;

  // special treatment for jpeg output: convert to sRGB

  std::string target_csp = args.colour_space ;

  if ( extension == ".JPG" || extension == ".jpg" )
  {
    if ( args.verbose )
      std::cout << "enforcing sRGB for JPEG output" << std::endl ;

    // internally, we're working in the scene_linear colour space, but
    // jpeg files must be in sRGB. The default output colour space is
    // also scene_linear, but for jpeg output we prescribe sRGB. This
    // might be extended for other output formats - for now the
    // automatic choice of sRGB output is limited to JPG.

    target_csp = "sRGB" ;
  }

  if ( args.working_colour_space != target_csp )
  {
    if ( args.verbose )
      std::cout << "converting from internal csp "
                << args.working_colour_space
                << " to " << target_csp
                << std::endl ;

    OIIO::ImageBufAlgo::colorconvert
      ( out_buf , out_buf ,
        args.working_colour_space ,
        target_csp ) ;
  }

  auto success = out_buf.write ( filename ) ;
  assert ( success ) ;
}

#include "zimt/bspline.h"

template < std::size_t NCH >
using px_t = zimt::xel_t < float , NCH > ;

// template < std::size_t NCH >
// using spl_t = zimt::bspline < px_t < NCH > , 2 > ;

template < std::size_t NCH >
bool read_image_data ( zimt::view_t < 2 , px_t < NCH > > & trg ,
                       // const facet_spec & fct ,
                       const std::string & filename ,
                       const std::string & colour_space ,
                       int & native_nchannels )
{
  typedef px_t < NCH > in_px_t ;

  if ( args.verbose )
    std::cout << "file " << filename
              << " is now loaded from disk" << std::endl ;

  // set up a config telling OIIO we want float data with the
  // given width and height.

  ImageSpec config ;
  config.format = TypeDesc::FLOAT ;
  config.width = trg.shape[0] ;
  config.height = trg.shape[1] ;

  // add further attributes gleaned from the CL to the config

  for ( const auto & attr : args.oiio_option_v )
  {
    std::string oiio_arg , oiio_type , oiio_val ;

    auto pos = attr.find_first_of ( "=" ) ;
    if ( pos != attr.npos )
    {
      auto lhs = attr.substr ( 0 , pos ) ;
      auto at_pos = lhs.find_first_of ( "@" ) ;
      if ( at_pos != lhs.npos )
      {
        // this argument is suffixed with an OIIO typestring
        oiio_arg = lhs.substr ( 0 , at_pos ) ;
        oiio_type = lhs.substr ( at_pos + 1 ) ;
      }
      else
      {
        // no typestring
        oiio_arg = lhs ;
        oiio_type = std::string() ;
      }
      // take the remainder after the '=' as the attribute's value
      oiio_val = attr.substr ( pos + 1 ) ;
    }
    else
    {
      oiio_arg = attr ;
      oiio_type = std::string() ;
      oiio_val = std::string() ;
    }
    if ( oiio_type.size() )
    {
      if ( args.verbose )
        std::cout << "processing typed oiio argument: " << oiio_arg
                  << " type: " << oiio_type
                  << " value: " << oiio_val
                  << std::endl ;

      // typed argument. OIIO recognizes it's own brand of
      // typestring.

      auto typedesc = TypeDesc ( oiio_type ) ;

      // with a type descriptor, we can process the value
      // in string form. The user should separate individual
      // values of multi-value rhs with space or tab.

      config.attribute ( oiio_arg , typedesc , oiio_val ) ;
    }
    else
    {
      if ( args.verbose )
        std::cout << "processing untyped oiio argument: "
                  << oiio_arg
                  << " value: " << oiio_val
                  << std::endl ;

      // untyped argument

      config [ oiio_arg ] = oiio_val ;
    }
  }

  // we set up an OIIO ImageBuf to pull in the image data, using
  // the config we've set up above

  // OIIO::ImageBuf in_buf ( config ,
  //                         (float*) (p_bspl->core.data()) ,
  //                         sizeof ( in_px_t ) ,
  //                         p_bspl->core.strides[1] * sizeof ( in_px_t )
  //                       ) ;

  OIIO::ImageBuf in_buf ( config ,
                          (float*) (trg.data()) ,
                          sizeof ( in_px_t ) ,
                          trg.strides[1] * sizeof ( in_px_t )
                        ) ;

  // To read the buffer with acknowledgement of the config, which
  // may contain additional config attributes, this syntax works.
  // but a format request (setting config.format) seems not to have
  // an effect and the buffer still contains the native data format.

  OIIO::ImageBuf read_buf ( filename , 0 , 0 , nullptr , &config ) ;

  // To carry on with float data, I next set up an ImageBuf 'occupying'
  // the b-spline's 'core' area and issue a 'copy' command with an
  // explicit type conversion (second parameter). directly reading into
  // such a buffer does not work reliably, because it does not honour
  // any additional config attributes (e.g raw:ColorSpace=ACES)

  in_buf.init_spec ( filename , 0 , 0 ) ;
  in_buf.copy ( read_buf , TypeDesc::FLOAT ) ;
  
  // re-read the spec, doublecheck data format

  const ImageSpec &spec = in_buf.spec() ;
  assert ( in_buf.spec().format == TypeDesc::FLOAT ) ;
  native_nchannels = in_buf.nchannels() ;

  // the facet may have inherent colour space information via
  // a 'Csp' parameter in a PTO i-line. This takes precedence:

  auto csp = colour_space ;

  // if there is no such spec, we rely on the OIIO attribute

  if ( csp == std::string() )
    csp = spec.get_string_attribute ( "oiio:ColorSpace" ) ;

  if ( args.verbose )
    std::cout << "facet's colour space: "
              << csp << std::endl ;

  // if the two coour spaces differ, we need to convert:

  if ( csp != args.working_colour_space )
  {
    if ( args.verbose )
      std::cout << "converting from facet's csp " << csp
                << " to internal csp "
                << args.working_colour_space << std::endl ;

    bool success = OIIO::ImageBufAlgo::colorconvert
      ( in_buf , in_buf ,
        csp ,
        args.working_colour_space ) ;

    assert ( success ) ;
  }

  // obtain the data via the ImageBuf, which will apply the
  // colour space transformation if it was specified above
  // by calling 'colorconvert' on the buffer
  
  bool success = in_buf.get_pixels ( OIIO::ROI() ,
                                     OIIO::TypeFloat ,
              trg.data() ,
              sizeof ( in_px_t ) ,
              trg.strides[1] * sizeof ( in_px_t ) ) ;
  
  return ( success ) ;
}

#endif // #ifndef ENVUTIL_BASIC_H
