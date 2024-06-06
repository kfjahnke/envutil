# envutil - a utility program to convert between lat/lon and cubemap environment maps and extract partial images from them

This is a stand-alone repository for envutil, which started out as a demo program
for my library [zimt](https://github.com/kfjahnke/zimt). The program has grown
beyond the limits of what I think is sensible for a demo program, and I also
think it's a useful tool.

This program takes a 2:1 lat/lon environment or a 1:6 cubemap image as input
and produces output in the specified orientation, projection, field of view
and extent. For CL arguments, try 'envutil --help'. Panorama photographers
may not be familiar with the term 'lat/lon environment' - to them, this is
a 'full spherical panorama' or 'full equirect'. The difference is merely
terminology, both are the same. Cubemaps are rarely used in panorama
photography, but I aim at reviving them by building a code base to
process them efficiently, hone the standard, and integrate processing of
single-image cubemap format (as it is used by envutil) into [lux](https://bitbucket.org/kfj/pv), my
image and panorama viewer - currently lux cubemaps need six separate images
and a special litle script and the code is less optimized than the code I
offer here. The single-image cubemap format I process in envutil leans
on the openEXR standard as far as the order and orientation of the cube
face images is concerned: the images are expected stacked vertically (in
this sequence: left, right, top, bottom, front, back; top and bottom align
with the 'back' image), but the precise geometry of the images may differ
from openEXR's standard - I measure cube face image field of view
'edge-to-edge', meaning that I consider pixels small square areas and
X degrees field of view extend over just as many of these small areas as
the image is wide and high. A different notion is 'center-to-center'
measuring, where the X degrees refer to the area enclosed between the
marginal pixels seen as points in space (the convex hull). Cubemaps with
this type of filed of view mesurement can be processed by passing --ctc 1
(for 'center-to-center'). envutil can also process cubemaps with cube face
images with larger field of view than ninety degrees (pass --cmbfov). 
Cubemaps can also be passed as six single cube face images. The images
have to follow a naming scheme which the caller passes via a format
string. You pass the format string as input, and it's used to generate six
filenames, where the format string is expanded with "left", "right"...
for the six cube faces. All six image file names generated in this fashion
must resolve to square images of equal size and field of view.

The output projection can be one of "spherical", "cylindrical",
"rectilinear", "stereographic", "fisheye" or "cubemap". The geometrical
extent of the output is set up most conveniently by passing --hfov, the
horizontal field of view of the output. The x0, x1, y0, and y1 parameters
allow passing specific extent values (in model space units), which should
rarely be necessary. To specify the orientation of the 'virtual camera',
pass Euler angles yaw, pitch and roll - they default to zero: a view
'straight ahead' to the point corresponding to the center of the
environment image with no camera roll. The size of the output is
given by --width and --height. You must pass an output filename
with --output; --input specifies the environment image. Note that you
can use 'envutil' to reproject environment images, e.g convert a cubemap
to a lat/lon image and the reverse. To get a full 360 degree lat/lon
image from a cubemap, you must pass --hfov 360 and --projection spherical,
and for the reverse operation, pass --hfov 90 and --projection cubemap.
Export of a cubemap as six separate images can be achieved by passing a
format string - the same mechanism as for input is used.

envutil can also produce image series. The individual images are all
created from the same environment, but the horizontal field of view
and yaw, pitch and roll of the virtual camera for each image are taken
from a 'sequence file', passed as --seqfile. The four values are given
four in a line, each line provides settings for one image. If you
pass a format string for the output, the images will be named accordingly;
use format specs like %03d (to make the images lexically sortable).
If output is not a format string, it's interpreted as a video file name,
and the consecutive images are used to create a video with ffmpeg.
This is a handy way to create video sequences of zooms and pans,
but you'll want a script to produce the seqfile - it's a bit like a
slicer file for a 3D printer - you don't want to write that by hand,
either. Default is H265 with 60 fps and 8 Mbit/sec.

You can choose several different interpolation methods with the --itp
command line argument. The default is --itp 1, which uses bilinear
interpolation. This is fast and often good enough, especially if there
are no great scale changes involved - so, if the output's resolution is
similar to the input's. --itp -1 employs OpenImageIO (OIIO for short)
for interpolation. Without further parameters, OIIO's default mode is
used, which uses sophisticated, but slow methods to produce the output.
All of OIIO's interpolation, mip-mapping and wrapping modes can be
selected by using the relevant additional parameters. Finally, --itp -2
uses 'twining' - inlined oversampling with subsequent weighted pixel
binning. The default with this method is to use a simple box filter
on a signal which is interpolated with bilinear interpolation; the
number of taps and the footprint of the pick-up are set automatically.
Additional parameters can change the fottprint and the amount of
oversampling and add gaussian weights to the filter parameters.
Twining is quite fast (if the number of filter taps isn't
very large); when down-scaling, the parameter 'twine' should be at
least the same as the scaling factor to avoid aliasing. When upscaling,
larger twining values will slighly soften the output and suppress the
star-shaped artifacts typical for bilinear interpolation. Twining is
new and this is a first approach. The method is intrinsically very
flexible (it's based on a generalization of convolution), and the full
flexibility isn't accessible in 'envutil' with the parameterization
as it stands now, but it's already quite useful with the few parameters
I offer. Automatic twining is also used for video output, where fixed
twining parameters would fail to adapt to changing hfov. Use of twining
or OIIO's lookup (which also adapts to the field of view) is strongly
recommended for video output to avoid aliasing.

The program uses [zimt](https://github.com/kfjahnke/zimt) as it's 'strip-mining' and SIMD back-end, and
sets up the pixel pipelines using zimt's functional composition tools.
This allows for terse programming, and the use of the functional
paradigm allows for many features to be freely combined - a property
which is sometimes called 'orthogonality'. What you can't combine in
'envutil' is twining and interpolation with OIIO - this is pointless,
because OIIO offers all the anti-aliasing and quality interpolation
one might want, and using twining on top would not improve the
output - it's meant as an alternative to OIIO's lookup code, but it's
a new development and hasn't been tested much. I'd welcome enterprising
users to compare it with other anti-aliasing and interpolation methods
and share the results!

Currently, the build is set up to produce binary for AVX2-capable
CPUs - nowadays most 'better' CPUs support this SIMD ISA. When
building for other (and non-i86) CPUs, suitable parameters should
be passed to the compiler (you'll have to modify the CMakeLists.txt).
I strongly suggest you install highway on your system - the build
will detect and use it to good effect. This is a build-time dependency
only. Next-best (when using i86 CPUs up to AVX2) is Vc, the fall-back
is to use std::simd, and even that can be turned off if you want to
rely on autovectorization; zimt structures the processing so that it's
autovectorization-friendly and performance is still quite good that way.

The program is built with CMake.

The only mandatory dependency is [OpenImageIO](https://github.com/AcademySoftwareFoundation/OpenImageIO) - OIIO in short. envutil also
links to some ffmpeg libraries, but these should come with an installation
of OpenImageIO - OIIO can use ffmpeg to open videos, but at the time of
this writing it could not write to videos, so I had to do it 'by hand',
using ffmpeg code directly.

If you have highway installed, the build will use it, if not it tries for Vc,
and if that isn't present either it uses std::simd. If you dont' want any of
the three, set the option USE_GOADING=ON. The functionality and results should
be the same with either of the back-ends, but processing time will vary,
especially if you are using 'twining'. My own library code (from the zimt
library) is provided here as a stripped-down copy - it's header-only C++
template metacode and doesn't link to anything, so it's a compile-time
dependency only.

It's recommended to build with clang++

To build, try this:

    mkdir build
    cd build
    cmake [options] ..
    make

recommended options for the build:

    -DCMAKE_CXX_COMPILER=clang++
    -DCMAKE_C_COMPILER=clang
    -DCPACK_GENERATOR=DEB

the last option is only useful if you want to build debian packages with
'make package'. Compilation with g++/gcc may or may not work - I don't check
regularly.

With version 0.1.1, on top of building on febian12,  I have also managed to
build envutil on on an intel mac running macOS 12.7.5 (using macPorts for the
dependencies) and on windows 11 using mingw64. I haven't yet manages to produce
video files with the mac build which played on my mac, but otherwise the builds
seem viable.
    
'make' should produce a binary named 'envutil' or 'envutil.exe''.

envutil --help gives a summary of command line options:

    --help                 Print help message
    -v                     Verbose output
    --input INPUT          input file name (mandatory)
    --output OUTPUT        output file name (mandatory)
    --width EXTENT         width of the output
    --height EXTENT        height of the output
    --projection PRJ       projection used for the output image(s)
    --hfov ANGLE           horiziontal field of view of the output
    --cbmfov ANGLE         horiziontal field of view of cubemap input (default: 90)
    --support_min EXTENT   minimal additional support around the cube face proper
    --tile_size EXTENT     tile size for the internal representation image
    --ctc CTC              pass '1' to interpret cbmfov as center-to-center (default 0)
    --yaw ANGLE            yaw of the virtual camera
    --pitch ANGLE          pitch of the virtual camera
    --roll ANGLE           roll of the virtual camera
    --x0 EXTENT            low end of the horizontal range
    --x1 EXTENT            high end of the horizontal range
    --y0 EXTENT            low end of the vertical range
    --y1 EXTENT            high end of the vertical range
    --seqfile SEQFILE      image sequence file name (optional)
    --codec CODEC          video codec for video sequence output (default: libx265)
    --mbps MBPS            output video with MBPS Mbit/sec (default: 8)
    --fps FPS              output video FPS frames/sec (default: 60)
    --itp ITP              interpolator: 1 for bilinear, -1 for OIIO, -2 bilinear+twining
    --twine TWINE          use twine*twine oversampling - default: automatic settings
    --twine_width WIDTH    widen the pick-up area of the twining filter
    --twine_sigma SIGMA    use a truncated gaussian for the twining filter (default: don't)
    --twine_threshold THR  discard twining filter taps below this threshold
    --tsoptions KVLIST     OIIO TextureSystem Options: coma-separated key=value pairs
    --swrap WRAP           OIIO Texture System swrap mode
    --twrap WRAP           OIIO Texture System twrap mode
    --mip MIP              OIIO Texture System mip mode
    --interp INTERP        OIIO Texture System interp mode
    --stwidth EXTENT       swidth and twidth OIIO Texture Options
    --stblur EXTENT        sblur and tblur OIIO Texture Options
    --conservative YESNO   OIIO conservative_filter Texture Option - pass 0 or 1
    
The input can be either a lat/lon environment image (a.k.a. 'full spherical' or 'full equirect' or '360X180 degree panorama') - or a 'cubemap' - a set of six square images in rectilinear projection showing the view to the six cardinal directions (left, right, up, down, front, back). The cubemap can be provided as a single image with the images concatenated vertically, or as six separate images with 'left', 'right' etc. in their - otherwise identical - filenames, which are introduced via a format string.

envutil only processes sRGB and linear RGB data, the output will be in the same
colour space as the input. If you use the same format for input and output, this
will automatically be the case, if not, you may get faulty output if the default
coulour spaces of the formats don't match - your output will look too bright or
too dark. The number of colour channels will also be the same in the output as
in the input - it's recommended you use the same file format for both, e.g.
produce JPEG output from JPEG input, but if the formats are compatible, you're
free to move from one to the other. To get mathematically correct results, note
that your input should be in linear RGB (e.g. from an openEXR file), because
internally all calculations are done *as if* the image data were linear RGB.
This is not mathematically correct for sRGB data. Future versions of envutil
may add a stage of processing to deal with non-linear-RGB input. Linear RGBA
is fine; the import with OIIO should provide associated alpha which can be
processed correctly by envutil as four-channel data, producing correct output
for images with alpha channel. Single-channel data are also fine, and even
two-channel data should work, though it's hard to find them 'out there' -
maybe anaglyphs? Channels beyond four are ignored.

# envutil Command Line Options

The options are given with a headline made from the argument parser's help
text. The capitalized word following the parameter is a placeholder for the
actual value you pass.

## --input INPUT         input file name (mandatory)

Any image file which has 2:1 aspect ratio will be accepted as lat/lon environment
map, and any image file with 1:6 aspect ratio will be accepted as a cubemap.
If you pass a string containing a percent sign, envutil will consider this as
a format string which can be used to specify six separate cube face images.
The string will be used like a 'normal' C format string to generate six filenames,
passing the strings "left", "right", "top", "bottom", "front" and "back", in turn,
to replace the format sequence in the string - you'd typically just use '%s' here.
Only one percent sign is allowed, and the images are all treated alike. And they
must, of course, all exist and have 1:1 aspect ratio. Once they've been loaded
into the internal representation, processing continues as if the input had been
a single cubemap image. When invoked with an input containing a percent sign,
envutil first tries to find the left cube face image. If that fails, the attempt
to load six images is aborted and instead envutil will assume you have actually
passed a filename containing a verbatim percent sign and proceed accordingly -
and likely fail.

## --output OUTPUT       output file name (mandatory)

The output will be stored under this name. If you are generating cubemaps
(pass --projection cubemap) you may pass a format string, just as for input.
You'll get six separate cube face images instead of the single-image cubemap.
When generating image sequences and you want single-image output, you also
have to pass a format string, but with a format element specifying an integer.
This is replaced with the output image's sequential number. Use a format
element like %03d to get alphabetically sortable images - in preference to
e.g. %d which does not produce leading zeroes. If you are generating image
sequences and your output parameter is *not* a format string, envutil will
produce video output to the given filename.

## --projection PRJ  target projection

Pass one of the supported output projections: "spherical", "cylindrical",
"rectilinear", "stereographic", "fisheye" or "cubemap". The default is
"rectilinear".

## --hfov ANGLE      horiziontal field of view of the output (in degrees)

The default here is ninety degrees, which all projections can handle.
Spherical, cylindrical and fisheye output will automatically periodize
the image if the hfov exceeds 360 degrees. Rectilinear images can't handle
more than 180 degrees hfov, and at this hfov, they won't produce usable
output. Stereographic images can theoretically accomodate 360 degrees fov,
but just as rectilinear images aren't usable near their limit of 180 degrees,
they aren't usable when their limit is approached: most of the content
becomes concentrated in the center and around that a lot of space is wasted
on a bit of content which is stretched extremely in radial direction. For
cubemaps, you should specify ninety degrees hfov, but you can produce
cubemaps with different hfov - they just won't conform to any standards
and won't be usable with other software unless that software offers
suitable options. envutil supports wider-angle cubemaps, see --cbmfov

## --width EXTENT    width of the output

in pixel units. For lat/lon environment images, this should be precisely
twice the height, so this value should be an even number. I recommeend that
you pick a multiple of a small-ish power of two (e.g. 64) to make it easier
for software wanting to mip-map the data. When producing cubemaps from
full sphericals, I recommend using a width which is ca. 1/pi times the width
of the input. For the reverse operation, just use four times the cubemap's
width. These factors preserve the resolution.

## --height EXTENT   height of the output

in pixel units. For cubemaps, this is automatically set to six times the width.
For spherical output, if height is not passed, it is set to half the width,
increasing 'width' to the next even value. For other projections, if you don't
pass 'height', the default is to use the same as the width, so to render a
square image.

# Additional Parameters for Cubemap Input

## --cbmfov ANGLE         horiziontal field of view of cubemap input (default: 90)

If the environment given as input is a cubemap, you can specify the horizontal
field of view of the cube face images with this parameter. The default is
precisely ninety degrees, but values greater than ninety are allowed as well.
You can create such wider-angle cubemaps by passing --hfov greater than ninety
when creating a cubemap with envutil.

## --support_min EXTENT  minimal additional support around the cube face proper

This is a technical value, probably best to leave it at it's default of 8.
when a cubemap is converted to it's internal representation, the ninety degree
'cube face proper' is surrounded with a frame of additional 'support' pixels
which help interpolation. You can specify here how wide this frame should be
at the least.

## --tile_width EXTENT   tile width for the internal representation image

Also best left at the default of 64.
This value takes care of widening the support frame further - if necessary - to
make the size of the square images in the internal representation a multiple
of this size. This does not magnify the image but adds pixels reprojected from
other cube faces, so it doesn't affect image quality. Using a small-ish power
of two here is especially useful when using OIIO for lookup, to help it
mip-map the texture generated from the internal representation - OIIO can't
natively process cubemap environments, so envutil generates a texture file
and feeds that to OIIO's texture system for the look-up.

## --ctc  flag indicating cube face fov is measured between marginal pixel centers

The standard way of measuring the field of view of cube face images in envutil is
to consider pixels as small square areas of constant colour with an extent of one
pixel unit. If an image is W pixels wide, a field of view of D degrees is taken to
coincide with the angle between rays to the left margin of the leftmost pixel and
the right margin of the rightmost pixel (same for top and botttom). If you pass
--ctc 1, D will instead coincide with the angle between rays to the centers of the
marginal pixels.
So usually (--ctc 0), we have D = atan ( f * W / 2 ), but with --ctc 1, we have
D = atan ( f * ( W - 1 ) / 2 ) This is hard to see, but some cubemaps seem to use this convention, and using them without --ctc will lead to subtle errors. Internally, envutil uses the first notion, and simply recalculates the field of view to be used internally to the slightly larger value which results form the edge-to-edge notion;
You could do the same 'manually' and pass a slightly higer value for cbmfov - using
--ctc is merely a convenience saving you the calculation. 

### A Side Note on lat/lon Images

Note that --ctc does not affect the processing of lat/lon environment maps - I may
add a separate option for that. lat/lon environment maps are always processed with
edge-to-edge semantics, assuming the image is periodic in the horizontal - the first
column follows again after the last. The other assumption about lat/lon images is
that they also follow edge-to-edge semantics for the vertical: the image lines at
the very top and bottom of the image represent a (very small) circle around the
pole with half a pixel width radius, opposed to some lat/lon formats which use
'center-to-center' semantics and repeat a single value (namely that for the pole)
for all pixels in the top, and another for the bottom row. From this it shoud be
clear why envutil expects lat/lon images to have precisely 2:1 aspect ratio.
Envutil honours the peculiarities of the spherical (a.k.a. equirectangular)
projection and interpolations near the poles will 'look at' pixels which
are nearby on the spherical surface, even if they are on opposite sides of the
pole, so you can e.g. safely extract nadir caps.

# Parameters for Single-Image Output

## --yaw ANGLE       yaw of the virtual camera (in degrees)
## --pitch ANGLE     pitch of the virtual camera (in degrees)
## --roll ANGLE      roll of the virtual camera (in degrees)

These three angles are applied to the 'virtual camera' taking the envutiled
view. They default to zero. It's okay to pass none or just one or two. yaw
is taken as moving the camera to the right, pitch is taken as upward movement,
and roll as a clockwise rotation. Note that the orientation of the *virtual
camera* is modified; when looking at the resulting images, objects seen on
them seem to move the opposite way. Negative values have the opposite effect.
Panorama photographers: to extract nadir patches, pass --pitch -90
These angles are known as the 'Euler Angles' and are easy to understand, as
opposed to the quaternions which envutil uses internally to represent rotations.

## --x0 EXTENT       low end of the horizontal range
## --x1 EXTENT       high end of the horizontal range
## --y0 EXTENT       low end of the vertical range
## --y1 EXTENT       high end of the vertical range

These are special values which can be used to specify the extent, in model
space units, of the output. This requires some understanding of the inner
workings of this program - if you use -v, the verbose output will tell you
for each rendering which extent values are generated from a field of view
parameter, given a specific projection. This can help you figure out specific
values you may want to pass, e.g. to produce anisotropic output or cropped
images.

# Parameters for Multi-Image and Video Output

## --seq SEQFILE         image sequence file name

Alternative output of an image sequence with hfov, yaw, pitch and roll given
as four consecutive numbers per line, one line per image. This is a new
experimental feature, output files are numbered consecutively when 'output'
can be interpreted as a format string (use something like --output img%03d.jpg).
If 'output' is not a format string, it's interpreted as the name of a video
file, which is put together from the image sequence using ffmpeg and the
codec given with --codec (default is libx265).

Single images can be combined into a video with ffmpeg, like so:

ffmpeg -f image2 -pattern_type glob -framerate 60 -i 'seq*.jpg' -s 960x540
       -c:v libx264 foo.mp4

This may be preferable, because all ffmpeg options can be exploited that way,
and the code in envutil is derived from a very old example program which may
not be at the height of time.

## --codec CODEC         video codec for video sequence output (default: libx265)

This is a string passed to ffmpeg to specify the video codec. Internally, frames
for video are encoded in YUV, which may not work with other codecs than the
hevc family (H264, H265). If direct video output fails, pass a format string
as output and then combine the separate images with ffmpeg (see --seqfile, above)
for H264 output, pass libx264. This is a bit enigmatic - I haven't yet figured
out which strings to pass here for all the codecs which ffmpeg supports, and
which of the codecs work with the pixel type I supply.

## --mbps MBPS           output video with MBPS Mbit/sec (default: 8.0)

Megabits/second value for the codec. This determines how much the data will be
compressed by the video codec. Note that this only affects video output - if you
generate separate images, they aren't affected by this parameter.

## --fps FPS             output video FPS frames/sec (default: 60)

This does not change the image sequence, but it changes the speed of the video
when it's displayed.

# Interpolation Options

## --itp ITP

This is an integer value determining the interpolation method. There are
currently three modes of interpolation:

    1 - use simple bilinear interpolation directly on the source image
        this is the fastest option, and unless there is a significant scale
        change involved, the output should be 'good enough' for most still
        image renditions. This is the default.

    -1 - use OpenImageIO's 'environemnt' or 'texture' function for lookup.
         without additional arguments, this will use a sophisticated
         interpolator with good antialiasing and bicubic interpolation.

    -2 - use 'twining' - this is a method which first super-samples and then
         combines several pixels to one output pixel ('binning'). This is my
         own invention. It's quite fast and produces good quality output.
         This method should see community review to compare it with other
         methods. If you only pass --itp -2 and not --twine, envutil will
         set up twining parameters calculated to fit well with the
         transformation at hand, and that's also done for image sequence
         output, where the parameters have to adapt to the changing geometry.

Why use negative values for ITP for the second and third mode? This is similar
to the values used in lux for 'decimators' - in lux, positive values are
reserved for degrees of a b-spline reconstruction filter used as low-pass
filter. bilinear interpolation is the same as a degree-1 b-spline, hence the
value 1 for bilinear interpolation. I may add other spline degrees to envutil,
so I reserve the positive numbers for future use.

# Twining-specific options

These options control the 'twining' filter. These options only have an effect if
you activate twining with --itp -2. envutil will use bilinear interpolation on
the source image for single-point lookups, but it will perform more lookups and
then combine several neighbouring pixels from the oversampled result into each
target pixel.

## --twine TWINE         use twine*twine oversampling and box filter

The effect of 'twining' is the same as oversampling and subsequent application
of a box filter. The filter is sized so that the oversampling is uniform over
the data, but the direct result of the oversampling is never saved - all
samples falling into a common output pixel are pooled and only the average
is stored. This keeps the pipeline 'afloat' in SIMD registers, which is fast
(as is the arithmetic) - especially when highway or Vc are used, which
increase SIMD performance particularly well.

If you activate twining by selecting --itp -2, but don't pass --twine
(or pass --twine 0 explicitly), envutil will set up twining parameters
automatically, so that they fit the relation of input and output. If the
output magnifies (in it's center), the twine width will be widened to
avoid star-shaped artifacts in the output. Otherwise, the twine factor
will be raised to avoid aliasing. You can see the values which envutil
calculates if you pass -v. If the effect isn't to your liking, you can
take the automatically calculated values as a starting point, but even if
you pass --twine, the values will adapt in an image sequence. --twine
only has an effect for single-image output

## --twine_width TWINE_WIDTH  widen the pick-up area of the twining filter

A second parameter affecting 'twining'. If the source image has smaller
resolution than the target image, the output reflects the interpolator's
shortcomings, so with e.g. bilinear interpolation and large scale change
(magnification) the output may show star-shaped and staircase artifacts.
A 'standard' twine with twine_width 1.0 will pool look-ups which all
correspond to target locations inside the *target* pixel's boundaries.
For magnifying views, this becomes ever more pointless with increasing
magnification - the look-up locations will all fall into a small area
of the source image - so, they'll be very much one like the other, and
pooling several of them is futile. So for magnifying views, you want to
widen the area in which look-ups are done to an area which is in the
same order of magnitude as a *source* image pixel. To get this effect,
try and pass a twine_width up to the magnitude of the scale change,
or rely on automatic twining, which calculates a good value for you.

Input with low resolution is often insufficiently band-limited which will
result in artifacts in the output or become very blurred when you try to
counteract the artifacts with excessive blurring. There's little to be
gained from scaling up anyway - the lost detail can't be regained. If
you don't get satisfactory results with twining and adequate twine_width,
you may be better off with one of the better OIIO interpolators, e.g.
bicubic.

The twine_width parameter also only affects single-image output - for image
sequences, it's set automatically to fit the relation of input and output.
If --twine is not passed (or passed 0) this parameter will be calculated
automatically, and any value passed here is then overridden by the
automatics.

## --twine_sigma TWINE_SIGMA  use a truncated gaussian for the twining filter (default: don't)

If you don't pass --twine_sigma, envutil will use a simple box filter to combine the result of supersampling into single output pixels values. If you pass twine_sigma, the kernel will be derived from a gaussian with a sigma equivalent to twine_sigma times the half kernel width. This gives more weight to supersamples near the center of the pick-up.  If you pass -v, the filter is echoed to std::cout for inspection.

You can combine this parameter with automatic twining - the twine factor and
the twine width will be calculated automatically, then the gaussian is applied
to the initially equal-weighted box filter. Keep in mind that applying gaussian
weights will 'narrow' the filter, so you may need to pass a larger twine_width
to counteract that effect. This may be a bit counterintuitive, because gaussians
are commonly associated with blurring - but a gaussian kernel produces 'sharper'
output than a box filter. Also consider the next parameter which eliminates weights
below a given threshold to save CPU time.

## --twine_threshold TWINE_THRESHOLD  discard twining filter taps below this threshold

If you pass twine_sigma, marginal twining kernel values may become quite small and using them as filter taps makes no sense. Pass a threshold here to suppress kernel values below the threshold. This is mainly to reduce processing time. Use -v to
display the kernel and see which kernel values 'survive' the thresholding.
This parameter makes no sense without --twine_sigma (see above): if all weights
are equal, they'd either all be above or below the threshold.

You can combine this parameter with automatic twining - the twine factor and
the twine width will be calculated automatically, then the twine_sigma is applied,
and finally the thresholding eliminates small weights.

After thresholding, the weights are 'normalized' to produce a filter with 'unit
gain' - you can see that all the weights add up to 1.0 precisely. Without the
normlization, just eliminating the sub-threshold taps would darken the output.

# OpenImageIO-specific options

These options can be used to configure the look-up with OpenImageIO's
'environemnt' and 'texture' functions. You can read up the options in the
[OpenImageIO documentation](https://openimageio.readthedocs.io/en/v2.5.11.0/texturesys.html#)

## --tsoptions KVLIST  OIIO TextureSystem Options: coma-separated key=value pairs

This argument can be used to pass a set of comma-separated key=value pairs
to the 'catch-all' OIIO texture system attribute 'options'. It's a handy way
to handle the transfer, because OIIO has a parser for such lists.

## --swrap WRAP        OIIO Texture System swrap mode
## --twrap WRAP        OIIO Texture System twrap mode

Separate wrapping modes, determining how texture access outside the texture
area proper is handled.

## --mip MIP           OIIO Texture System mip mode

mip-mapping mode. I use a default of 'automip=1' for tsoptions

## --interp INTERP     OIIO Texture System interp mode

interpolator. This is the interpolator OIIO uses for it's lookup, not the
interpolator specified with the --itp option. To use OIIO's lookup, you pass
--itp -1. Then you can use --interp to specify OIIO's interpolator, like
--interp InterpSmartBicubic

## --stwidth EXTENT    swidth and twidth OIIO Texture Options

I lump together the swidth and twidth options - OIIO allows to pass them
separately, but I don't see a need to do so in envutil. The most useful
value here is to pass zero, which makes OIIO ignore the pickup-point's
derivatives. This can speed up the calculations.

## --stblur EXTENT      sblur and tblur OIIO Texture Options

I also lump together pre-blur along the s and t axis. This is to add more blur
to the processing.

## --conservative YESNO  OIIO conservative_filter Texture Option

pass this to switch on OIIO's 'conservative filter' on or off (pass 0 or 1);
the default is to have it on.

# Additional Technical Notes

One problem with cubemaps is that they are normally stored as concatenations of
six square images with precisely ninety degrees fov (field of view). This makes
access to pixels near the edges tricky - if the interpolator used for the purpose
needs support, marginal pixels don't have it all around, and it has to be gleaned
from adjoining cube faces, reprojecting content to provide the support. envutil
uses an internal representation (IR for short) of the cubemap which holds six
square images with a fov slightly larger than ninety degrees, where the part
exceeding the 'ninety degrees proper' cubeface image is augmented with pixels
reprojected from adjoining cube faces. With this additional 'frame', interpolators
needing support can operate without need for special-casing access to marginal
pixels. The formation and use of the IR image is transparent, it's used automatically.
There are two parameters which influence the size of the 'support margin', namely
--support_min and --tile_size - usually it's best to leave them at their defaults.
To access pixels in lat/lon environment maps which are marginal, envutil exploits the
inherent periodicity of the lat/lon image - simple periodicity in the horizontal and
'over-the-pole' periodicity in the vertical (that's reflecting in the vertical plus
an offset of half the image's width in the horizontal). One can in fact generate an
image from a full spherical which is periodic vertically by cutting of the right
half, rotating it by 180 degrees and pasting it to the bottom of the first half.
lux does just that for it's 'spherical filter', so that the image pyramids used
internally can be created with geometrically correct down-scaling. envutil does
not use mip-mapping or other pyramid-like code itself (even though OIIO's texture
system code does so), so this is merely a technical hint.

As an alternative to the antialiasing and interpolation provided by OIIO, envutil
offers processing with bilinear interpolation and it's own antialiasing filter, using
a method which I call 'twining'. Both for OIIO-based lookup and for twining, the
beginning of the pixel pipeline receives not just one 3D directional 'ray' coordinate
pointing 'into' the environment, but three: the second and third ray are calculated
to coincide with a target coordinate one discrete step toward the right or downwards,
respectively, in the *target* image. The pixel pipeline which receives these sets
of three ray coordinates can glean the difference between the first ray and the
two other two rays and adapt the lookup based on this additional information. For
lookup, simple differencing produces an approximation of the derivatives, or,
to interpret the difference differently, a pair of 3D vectors which can be used
to construct a plane and place several lookup points on that plane, to combine
their results to form the output as a weighted sum.

From visual inspection of the results, envutil's use of OIIO's lookup seems to produce
quite soft images. OIIO's lookup is - with the settings envutil uses - conservative
and makes an effort to avoid aliasing, erring on the side of caution. Of course I
can't rule out that my use of OIIO's lookup is not coded correctly, but I'm quite
confident that I've got it right. As of this writing, OIIO's lookup offers
functions with SIMD parameter set, processing batches of coordinates. Internally,
though, the coordinates are processed one after the other (albeit with use of
vertical SIMDization for the single-coordinate lookups). This is one of the reasons
why this process isn't very fast. OIIO's method of inspecting the derivatives has
the advantage of being perfectly general, and the lookup can decide for each pixel
whether it needs to be interpolated or antialiased - and how much.

In contrast, using bilinear interpolation (--itp 1) is very fast, but it's not
adequate for all situations - especially not if there are large differences in
resolution between input and output. If the scale of the output is much smaller
than the input's, the use of 'twining' should provide adequate antialiasing.
When scaling up, bilinear interpolation only produces visible artifacts if the
scale change is very large (the typical 'star-shaped artifacts'), which makes
little sense, because there's no gain to be had from scaling up by large
factors - the content won't improve. Up-scaling with OIIO's lookup uses bicubic
interpolation, which may be preferable. Ultimately it's up to the user to find
a suitable process by inspecting the output. envutil's 'twining' can cover a
certain range of input/output relations with a given parameter set, but for
extreme distortions, it may be suboptimal. This should be analyzed in depth,
for now, it's just a promising new method which seems to work 'well enough'.

I think that twining offers a good compromise beween speed and quality. I'm sure
I'm not the first one to think of this method, but I haven't done research to
see if I can find similar code 'out there'. What I am sure of is that my
implementation is fast due to the use of multithreaded horizontal SIMD code
through the entire processing chain, so I think it's an attractive offer. I'd
welcome external evaluation of the results and a discussion of the methods;
please don't hesitate to open issues on the issue tracker if you'd like to
discuss or report back! In my tests, images generated with twining seemed to
come out somewhat 'crisper' than renditions with OIIO's default mode, so they
may be preferable for photographic work, whereas OIIO's renditions make extra
sure there is no aliasing, which is good for video output.

There seems to be ambiguity of what constitutes a 'correct' cube face image with
ninety degrees field of view. In envutil, I code so that each pixel is taken to
represent a small square section of the image with constant colour. So the
'ninety degrees proper' extends form the leftmost pixel's left margin to the
rightmost pixel's right margin. Some cubemap formats provide images where the
centers of the marginal pixels coincide with the virtual cube's edges, repeating
each edge in the image it joins up with in the cube. If you process such cubemaps
with envutil, pass --ctc 1 (which stands for center-to-center). Otherwise, there
will be subtle errors along the cube face edges which can easily go unnoticed.
Make sure you figure out which 'flavour' your cubemaps are.
