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
  "biatan6"
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

struct facet_spec
{
  int facet_no ;
  std::string filename ;
  std::string projection_str ;
  projection_t projection ;
  double hfov ;
  double step ;
  double yaw , pitch , roll ;
  std::size_t width ;
  std::size_t height ;
  std::size_t nchannels ;
  double tr_x , tr_y , tr_z ;
  double tp_y , tp_p , tp_r ;

  bool init ( int argc , const char ** argv ) ;
} ;

struct arguments
{
  bool verbose ;
  std::string output ;
  double hfov ;
  std::size_t width ;
  std::size_t height ;
  std::string prj_str ;
  projection_t projection ;

  std::size_t support_min ;
  std::size_t tile_size ;

  double yaw , pitch , roll ;
  double x0 , x1 , y0 , y1 ;
 
  std::string seqfile ;
  std::string codec ;
  float mbps ;
  int fps ;

  int itp ;
  int prefilter_degree ;
  int spline_degree ;
  int twine  ;
  std::string twf_file ;
  bool twine_normalize ;
  bool twine_precise ;
  double twine_width , twine_density , twine_sigma , twine_threshold ;
  std::vector < zimt::xel_t < float , 3 > > twine_spread ;

  // gleaned from other parameters or input images

  double step ;
  std::size_t nchannels ;

  // technical variables for the argument parser

  std::string metamatch ;
  std::regex field_re ;

  // the 'arguments' object's 'init' takes the main program's argc
  // and argv.

  void init ( int argc , const char ** argv ) ;

  std::size_t nfacets ;
  std::vector < std::string > facet_name_v ;
  std::vector < std::string > facet_projection_v ;
  std::vector < std::string> facet_hfov_v ;
  std::vector < std::string > facet_yaw_v ;
  std::vector < std::string > facet_pitch_v ;
  std::vector < std::string > facet_roll_v ;
  std::vector < std::string > facet_trx_v ;
  std::vector < std::string > facet_try_v ;
  std::vector < std::string > facet_trz_v ;
  std::vector < std::string > facet_tpy_v ;
  std::vector < std::string > facet_tpp_v ;
  std::vector < std::string > facet_tpr_v ;
  std::vector < facet_spec > facet_spec_v ;
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

// quick shot at encoding video: I'm using an example file from the ffmpeg
// examples section, with some adaptations to bend the C code to C++,
// and additions to convert the float RGB data from the pixel pieline 
// to YUV for encoding with h264. C++ style comments (begining with //)
// are mine, to explain where I made alterations of the original C
// code or added stuff - and why.

extern "C"
{
#include <libavcodec/avcodec.h>
#include "libavformat/avformat.h"
#include <libavutil/opt.h>
#include <libavutil/imgutils.h>
#include <libswscale/swscale.h>
} ;

// The code example's end is marked with 'end of copied ffmpeg example
// code'

// Begin copied example code:

/*
 * Copyright (c) 2001 Fabrice Bellard
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/**
 * @file libavcodec encoding video API usage example
 * @example encode_video.c
 *
 * Generate synthetic video data and encode it to an output file.
 */

static void encode(AVCodecContext *enc_ctx, AVFrame *frame, AVPacket *pkt,
                   FILE *outfile)
{
    int ret;

    /* send the frame to the encoder */

    // the format strings used for the printf don't seem to work with
    // C++, so I commented them out:
  
    // if (frame)
    //     printf("Send frame %3"PRId64"\n", frame->pts);

    ret = avcodec_send_frame(enc_ctx, frame);
    if (ret < 0) {
        fprintf(stderr, "Error sending a frame for encoding\n");
        exit(1);
    }

    while (ret >= 0) {
        ret = avcodec_receive_packet(enc_ctx, pkt);
        if (ret == AVERROR(EAGAIN) || ret == AVERROR_EOF)
            return;
        else if (ret < 0) {
            fprintf(stderr, "Error during encoding\n");
            exit(1);
        }

        // printf("Write packet %3"PRId64" (size=%5d)\n", pkt->pts, pkt->size);
        fwrite(pkt->data, 1, pkt->size, outfile);
        av_packet_unref(pkt);
    }
}

// what was in the original example is now the body of a C++ object,
// and the C code's global variables are now class members. main() is
// replaced by the class' c'tor which sets up the object for sustained
// operation with the given parameterization, and the code to push a
// single frame to the codec is in a member function and will be called
// repeatedly by the process generating images from a sequence file.

struct frame_sink
{
  const char *filename, *codec_name;
  const AVCodec *codec;
  AVCodecContext *c= NULL;
  int i, ret, x, y;
  FILE *f;
  AVFrame *frame;
  AVPacket *pkt;
  uint8_t endcode[4] { 0, 0, 1, 0xb7 };

  frame_sink ( const char * _filename ,
               const char * _codec_name )
  : filename ( _filename ) ,
    codec_name ( _codec_name )
  {
    /* find the mpeg1video encoder */
    codec = avcodec_find_encoder_by_name(codec_name);
    if (!codec) {
        fprintf(stderr, "Codec '%s' not found\n", codec_name);
        exit(1);
    }

    c = avcodec_alloc_context3(codec);
    if (!c) {
        fprintf(stderr, "Could not allocate video codec context\n");
        exit(1);
    }

    pkt = av_packet_alloc();
    if (!pkt)
        exit(1);

    /* put sample parameters */
    c->bit_rate = args.mbps ;
    /* resolution must be a multiple of two */
    c->width = args.width ;
    c->height = args.height ;
    /* frames per second */
    c->time_base = (AVRational){1, args.fps};
    c->framerate = (AVRational){args.fps, 1};

    /* emit one intra frame every ten frames
     * check frame pict_type before passing frame
     * to encoder, if frame->pict_type is AV_PICTURE_TYPE_I
     * then gop_size is ignored and the output of encoder
     * will always be I frame irrespective to gop_size
     */
    c->gop_size = 10;
    c->max_b_frames = 1;

    // for now, there is no choice about the codec, it's hardcoded, and
    // this codec seems to require YUV420P pixels. H264 works with the
    // same code and is faster.

    c->pix_fmt = AV_PIX_FMT_YUV420P;

    if (codec->id == AV_CODEC_ID_H265)
        av_opt_set(c->priv_data, "preset", "slow", 0);

    /* open it */
    ret = avcodec_open2(c, codec, NULL);
    if (ret < 0) {
        fprintf(stderr, "Could not open codec: %s\n", av_err2str(ret));
        exit(1);
    }

    f = fopen(filename, "wb");
    if (!f) {
        fprintf(stderr, "Could not open %s\n", filename);
        exit(1);
    }

    frame = av_frame_alloc();
    if (!frame) {
        fprintf(stderr, "Could not allocate video frame\n");
        exit(1);
    }
    frame->format = c->pix_fmt;
    frame->width  = c->width;
    frame->height = c->height;

    ret = av_frame_get_buffer(frame, 0);
    if (ret < 0) {
        fprintf(stderr, "Could not allocate the video frame data\n");
        exit(1);
    }

  }

  template < std::size_t nchannels >
  void encode_frame ( const zimt::view_t
                        < 2 ,
                          zimt::xel_t < float , nchannels >
                        > & pixels )
  {
    static struct SwsContext *img_convert_ctx;
    static int sws_flags = SWS_BICUBIC;
    if (img_convert_ctx == NULL)
    {
      // here we set up the converter to move from RGB to YUV. This is an
      // adaptation from another ffmpeg example file by the same author,
      // licensed as the other ffmpeg example code:

      /*
       * Libavformat API example: Output a media file in any supported
       * libavformat format. The default codecs are used.
       *
       * Copyright (c) 2003 Fabrice Bellard
       * ... license follows, see above */

      img_convert_ctx = sws_getContext(c->width, c->height,
                        // we say we have incoming RGB24. 
                                       AV_PIX_FMT_RGB24,
                                       c->width, c->height,
                        // and we want outgoing YUV420P
                                       AV_PIX_FMT_YUV420P,
                                       sws_flags,
                                       NULL, NULL, NULL );
    }
    if (img_convert_ctx == NULL)
    {
      fprintf(stderr, "Cannot initialize the conversion context\n");
      exit(1);
    }

    // end of code adapted from the second ffmpeg example, back to the
    // first example

    static int i = 0;

    // instead of having this code in a loop, we now have it in a
    // separate function which is called once for each frame.

    fflush(stdout);

    /* Make sure the frame data is writable.
        On the first round, the frame is fresh from av_frame_get_buffer()
        and therefore we know it is writable.
        But on the next rounds, encode() will have called
        avcodec_send_frame(), and the codec may have kept a reference to
        the frame in its internal structures, that makes the frame
        unwritable.
        av_frame_make_writable() checks that and allocates a new buffer
        for the frame only if necessary.
      */
    ret = av_frame_make_writable(frame);
    if (ret < 0)
        exit(1);

    // internally, we have the frame in float RGB(A). We only output
    // three channels. Alpha should be associated and result in darkening
    // the pixels without additional processing. If it were unassociated,
    // we'd have to multiply with alpha now. TODO: check
    // the float values are in the range [0,1], just as OIIO initially
    // produced them from whichever input file was given, so now we
    // have to scale up to [0,255]. An alternative would be to hold
    // a zimt array of YUV values and chain the conversion to the pixel
    // pipeline in the 'act' functor. But this conversion isn't the most
    // time-consuming operation, it's running the H265 codec (at least
    // here - my GPU doesn't have it).
    // Note that the code seems to expect sRGB data. There is no provision
    // for processing linear RGB yet, this should be added later, plus,
    // eventually, handling of other colour spaces.

    unsigned char buffer [ args.width * args.height * 3 ] ;
    for ( int y = 0 ; y < args.height ; y++ )
    {
      for ( int x = 0 ; x < args.width ; x++ )
      {
        const auto & px ( pixels [ { x , y } ] ) ;
        buffer [ 3 * y * args.width + 3 * x ] = 255.0f * px[0] ;
        if constexpr ( nchannels > 1 )
          buffer [ 3 * y * args.width + 3 * x + 1 ] = 255.0f * px[1] ;
        else
          buffer [ 3 * y * args.width + 3 * x + 1 ] = 255.0f * px[0] ;
        if constexpr ( nchannels > 2 )
          buffer [ 3 * y * args.width + 3 * x + 2 ] = 255.0f * px[2] ;
        else
          buffer [ 3 * y * args.width + 3 * x + 2 ] = 255.0f * px[0] ;
      }
    }

    // we need to pass a pointer to a pointer, hence:

    unsigned char * p_buffer = buffer ;

    // w is in bytes, not pixels, hence:
  
    int w = args.width * 3 ;

    // now do the type conversion. we don't have a size change.

    sws_scale ( img_convert_ctx,
                &p_buffer, &w, 0, c->height,
                frame->data, frame->linesize); 

    frame->pts = i;

    /* encode the image */
    encode(c, frame, pkt, f);

    if ( args.verbose )
      std::cout << "encoded video frame " << i << std::endl ;

    ++i ;
  }

  // the d'tor has the code needed to shut everything down gracefully.

  ~frame_sink()
  {
    /* flush the encoder */
    encode(c, NULL, pkt, f);

    /* Add sequence end code to have a real MPEG file.
       It makes only sense because this tiny examples writes packets
       directly. This is called "elementary stream" and only works for some
       codecs. To create a valid file, you usually need to write packets
       into a proper file format or protocol; see mux.c.
     */
    if (codec->id == AV_CODEC_ID_MPEG1VIDEO || codec->id == AV_CODEC_ID_MPEG2VIDEO)
        fwrite(endcode, 1, sizeof(endcode), f);
    fclose(f);

    avcodec_free_context(&c);
    av_frame_free(&frame);
    av_packet_free(&pkt);
  }
} ;

// end of copied ffmpeg example code

// here's our interface to the modified ffmpeg example code:
// this function takes a frame's worth of data encoded as a
// zimt::view, sets up the frame_sink object when first called
// and uses it to encode one frame per call. Alternatively,
// if 'args.output' is parsed as a format string, individual
// frame images are produced.

template < std::size_t nchannels >
void push_video_frame ( const zimt::view_t
                        < 2 ,
                          zimt::xel_t < float , nchannels >
                        > & pixels )
{
  static image_series is ( args.output ) ;
  static int index = 0 ;
  if ( is.valid() )
  {
    save_array ( is [ index ] , pixels ) ;
    ++index ;
  }
  else
  {
    static frame_sink video ( args.output.c_str() , args.codec.c_str() ) ;
    video.encode_frame ( pixels ) ;
  }
}

#endif // #ifndef ENVUTIL_BASIC_H
