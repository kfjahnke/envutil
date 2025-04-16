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

#include <regex>
#include <OpenImageIO/imageio.h>

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

  void get_image_metrics()
  {
    // currently building with raw::user_flip set to zero, to load
    // raw images in memory order without EXIF rotation. This only
    // affects raw images.

    ImageSpec config;
    config [ "raw:user_flip" ] = 0 ;
    auto inp = ImageInput::open ( filename , &config ) ;

    if ( ! inp )
    {
      std::cerr << "failed to open facet image '"
                << filename << "'" << std::endl ;
      exit ( -1 ) ;
    }

    const ImageSpec &spec = inp->spec() ;

    width = window_width = spec.width ;
    height = window_height = spec.height ;
    window_x_offset = window_y_offset = 0 ;
    nchannels = spec.nchannels ;
    inp->close() ;
  }
} ;

struct arguments
: public facet_base
{
  bool verbose ;
  std::string output ;
  std::string pano ;
  std::string split ;
  std::size_t support_min ;
  std::size_t tile_size ;
  std::string pto_file ;
  // std::string seqfile ;
  // std::string codec ;
  // float mbps ;
  // int fps ;
  int prefilter_degree ;
  int spline_degree ;
  int twine  ;
  std::string twf_file ;
  bool twine_normalize ;
  bool twine_precise ;
  double twine_width , twine_density , twine_sigma , twine_threshold ;
  std::vector < zimt::xel_t < float , 3 > > twine_spread ;

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

// helper function to save a zimt array of pixels to an image file, or
// to a set of six cube face images, if 'output' has a format string.

template < std::size_t nchannels >
void save_array ( const std::string & filename ,
                  const zimt::view_t
                    < 2 ,
                      zimt::xel_t < float , nchannels >
                    > & pixels ,
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

  auto out = ImageOutput::create ( filename );
  assert ( out != nullptr ) ;
  ImageSpec ospec ( pixels.shape[0] , pixels.shape[1] ,
                    nchannels , TypeDesc::HALF ) ;
  out->open ( filename , ospec ) ;

  auto success = out->write_image ( TypeDesc::FLOAT , pixels.data() ) ;
  assert ( success ) ;
  out->close();
}

#endif // #ifndef ENVUTIL_BASIC_H
