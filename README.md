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
process them efficiently, honing the standard, and integrating processing of
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
this type of field of view mesurement can be processed by passing --ctc 1
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

envutil can also 'mount' images in the supported five projections as
input (use --mount ...) - currently, the images have to be cropped
symmetrically, meaning that the optical axis is taken to pass through
the image center. For mounted input, you can specify the horizontal
field of view and the projection together with the image filename.
If the view 'looks at' areas not covered by the mounted image, it's
painted black or transparent black for RGBA.

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
command line argument. The default is --itp 1, which uses b-spline
interpolation. This is fast and often good enough, especially if there
are no great scale changes involved - so, if the output's resolution is
similar to the input's. --itp -1 employs OpenImageIO (OIIO for short)
for interpolation. Without further parameters, OIIO's default mode is
used, which uses sophisticated, but slow methods to produce the output.
All of OIIO's interpolation, mip-mapping and wrapping modes can be
selected by using the relevant additional parameters. Finally, --itp -2
uses 'twining' - inlined oversampling with subsequent weighted pixel
binning. The default with this method is to use a simple box filter
on a signal which is interpolated with b-spline interpolation; the
number of taps and the footprint of the pick-up are set automatically.
Additional parameters can change the footprint and the amount of
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

There is now code to read the twining filter kernel from a file
(use --twf_file); this filter can be scalled by additionally passing
--twine_width. This allows for arbitrary filters.

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

Currently, single-ISA build are set up to produce binary for AVX2-capable
CPUs - nowadays most 'better' x86 CPUs support this SIMD ISA. When
building for other (and non-i86) CPUs, suitable parameters should
be passed to the compiler (you'll have to modify the CMakeLists.txt).
I strongly suggest you install highway on your system - the build
will detect and use it to good effect. This is a build-time dependency
only. Next-best (when using i86 CPUs up to AVX2) is Vc, the fall-back
is to use std::simd, and even that can be turned off if you want to
rely on autovectorization; zimt structures the processing so that it's
autovectorization-friendly and performance is still quite good that way.

With version 0.1.2 I have changed the code to cooperate with highway's
foreach_target mechanism, which is now the default (but requires highway).
The resulting binary will contain specialized machine code for each ISA
which is common on a given CPU line (e.g. SSE*, SSSE3, AVX and AVX2 on
x86 CPUs) and dispatch internally to the best ISA which the current CPU
can process. This gives the resulting binary a much wider 'feeding
spectrum' - it should run on any x86 CPU with near-optimal ISA - for
x86 highway produces machine code up to AVX512 - and never produce
illegal instruction errors, because the CPU detection makes sure no
'too good' ISA instructions will be executed. If you prefer single-ISA
builds - e.g. if you're producing a binary only for a specific machine
or if you're modifying the code (recompilation is much quicker with
only a single ISA) - pass -DMULTI_SIMD_ISA=OFF to cmake. multi-ISA
builds require highway, but they can also be used with zimt's 'goading'
back-end. Pass -DUSE_GOADING=ON to make the build ignore explicit SIMD
libraries (highway, Vc or std::simd) which would otherwise be used in
this order of preference.

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

With version 0.1.1, on top of building on debian12,  I have also managed to
build envutil on on an intel mac running macOS 12.7.5 (using macPorts for the
dependencies) and on windows 11 using mingw64. I haven't yet managed to produce
video files with the mac build which played on my mac, but otherwise the builds
seem viable.
    
'make' should produce a binary named 'envutil' or 'envutil.exe''.

envutil --help gives a summary of command line options:

    --help                         Print help message
    -v                             Verbose output
    mandatory options:
      --input INPUT                  input file name (mandatory)
      --output OUTPUT                output file name (mandatory)
    important options which have defaults:
      --projection PRJ               projection used for the output image(s) (default: rectilinear)
      --hfov ANGLE                   horiziontal field of view of the output (default: 90)
      --width EXTENT                 width of the output (default: 1024)
      --height EXTENT                height of the output (default: same as width)
    additional input parameters for cubemap input:
      --cbmfov ANGLE                 horiziontal field of view of cubemap input (default: 90)
      --cbm_prj PRJ                  in-plane transformation for cube faces (default: "cubemap")
      --support_min EXTENT           minimal additional support around the cube face proper
      --tile_size EXTENT             tile size for the internal representation image
      --ctc CTC                      pass '1' to interpret cbmfov as center-to-center (default 0)
    additional parameters for single-image output:
      --yaw ANGLE                    yaw of the virtual camera
      --pitch ANGLE                  pitch of the virtual camera
      --roll ANGLE                   roll of the virtual camera
      --x0 EXTENT                    low end of the horizontal range
      --x1 EXTENT                    high end of the horizontal range
      --y0 EXTENT                    low end of the vertical range
      --y1 EXTENT                    high end of the vertical range
    additional parameters for multi-image and video output:
      --seqfile SEQFILE              image sequence file name (optional)
      --codec CODEC                  video codec for video sequence output (default: libx265)
      --mbps MBPS                    output video with MBPS Mbit/sec (default: 8)
      --fps FPS                      output video FPS frames/sec (default: 60)
    interpolation options:
      --itp ITP                      interpolator: 1 for spline, -1 for OIIO, -2 spline+twining
      --prefilter DEG                prefilter degree (>= 0) for the spline used with --ipt 1
      --degree DEG                   degree of the spline (>= 0) used with --ipt 1
    parameters for twining (with --itp -2):
      --twine TWINE                  use twine*twine oversampling - default: automatic settings
      --twf_file TWF_FILE            read twining filter kernel from TWF_FILE
      --twine_normalize              normalize twining filter weights gleaned from a file
      --twine_precise                project twining basis vectors to tangent plane
      --twine_width WIDTH            widen the pick-up area of the twining filter
      --twine_density DENSITY        increase tap count of an 'automatic' twining filter
      --twine_sigma SIGMA            use a truncated gaussian for the twining filter (default: don't)
      --twine_threshold THR          discard twining filter taps below this threshold
    parameters for lookup with OpenImageIO (with --itp -1):
      --tsoptions KVLIST             OIIO TextureSystem Options: coma-separated key=value pairs
      --swrap WRAP                   OIIO Texture System swrap mode
      --twrap WRAP                   OIIO Texture System twrap mode
      --mip MIP                      OIIO Texture System mip mode
      --interp INTERP                OIIO Texture System interp mode
      --stwidth EXTENT               swidth and twidth OIIO Texture Options
      --stblur EXTENT                sblur and tblur OIIO Texture Options
      --conservative YESNO           OIIO conservative_filter Texture Option - pass 0 or 1
    parameters for mounted image input:
      --mount IMAGE PROJECTION HFOV  load non-environment source image
      --facet IMAGE PROJECTION HFOV YAW PITCH ROLL
                                    load oriented non-environment source image
    
The input can be either a lat/lon environment image (a.k.a. 'full spherical' or 'full equirect' or '360X180 degree panorama') - or a 'cubemap' - a set of six square images in rectilinear projection showing the view to the six cardinal directions (left, right, up, down, front, back). The cubemap can be provided as a single image with the images concatenated vertically, or as six separate images with 'left', 'right' etc. in their - otherwise identical - filenames, which are introduced via a format string.

On top of 'full' environment images, envutil can also process 'mounted' images:
images with a given projection and hfov which may or may not cover the entire
360X180 degrees. If the 'full environment' isn't covered, lookup in areas outside
the mounted image are returned as black or transparent black. A variation of
mounted images are 'facets' - they aren't mounted 'straight ahead' but assigned
an orientation with three Euler angles.

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
"rectilinear", "stereographic", "fisheye", "cubemap" or "biatan6". The
default is "rectilinear". "biatan6" is a recent addition: it's a cubemap
with an additional in-plane transformation to sample the sphere more
evenly than can be done with rectilinear cube faces. See the remarks on
biatan6 in the "--cbm_prj PRJ" chapter.

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

## --cbm_prj PRJ          in-plane transformation for cube faces (default: "cubemap")

On top of the default, "cubemap", which does not use an in-plane transformation
(the cube faces are in rectilinear transformation and used just so), envutil
now supports "biatan6" in-plane transformation. This is a transformation which,
when used to create a cubemap, The result is a compression of the content
towards the edges of the cube faces (where it is normally unduly stretched
because of the rectilinear projection) and a widening in the center, where
there is now more space due to the compression near the edges. Where the
sample steps along a horizon line (central horizontal of one of the surrounding
four cube faces) correspond with rays whose angular difference decreases
towards the edges with rectilinear projection, with this in-plane 
transformation the corresponding rays are all separated by the same angular
step. There is still a distortion which can't be helped (we're modelling
a curved 2D manifold on a plane), but it amounts to a maximum of 4/pi in
the center of a vertical near a horizontal edge (or vice versa) due to the
scaling factor used both ways. Why 'biatan6'? because it uses the arcus tangens
(atan) in both the vertical and horizontal of all six cube faces. The tangens
and it's inverse, the arcus tangens, translate between an angle and the length
of a tangent - the latter is what's used by the rectlinear projection, which
projects to the tangent plane, and is the default for cube faces. Using
equal angular steps is closer to an even sampling of the sphere. The in-plane
transformation functions are used on the in-plane coordinates of the cube faces,
which, in envutil, range from -1 to 1 in 'model space' (the cube is modelled
to have unit center-to-plane distance, so as to just enclose a unit sphere).
There is a pair of them, which are inverses of each other, meaning that
applying them to a 2D coordinate one after the other will reproduce the
initial value. These are the two functions:

    float ( 4.0 / M_PI ) * atan ( in_face ) ;

    tan ( in_face * float ( M_PI / 4.0 ) ) ;

Using transcendental functions (tan, atan) is costly in terms of CPU
cycles, but both functions have SIMD implementations which make the cost
acceptable, because they can still process several values at once. An
alternative would be to use a pair of two functions roughly modelling
the two curves above, but with the same property of being inverses of
each other. I leave this option for further development - another idea
would be to introduce the in-plane transformation as a functional paramter,
so that user code can 'slot in' such a transformation.
The advantage of the 'biatan' transformation is that it transforms each
2X2 square to another 2X2 square - if one were to use e.g. spherical
transformation, there would be redundant parts in several images. So
with biatan transformation, each point in the cubemap has precisely one
correspondence on the sphere, just as with rectilinear projection. 

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

This may be preferable, because all ffmpeg's options can be exploited that way,
and the code in envutil is derived from a very old example program which may
not be at the height of time, and my understanding of ffmpeg is not very good.
On my system, the native video player displays the video, but vlc and lux
don't - obviously I am missing something, and help with ffmpeg would be welcome!
This feature might be modified slightly to accept parameter sets from a stream.
Here's an example of the first few lines of a sequfile for a plain pan with
camera hfov of 90 degrees in 1-degree steps (second column):

    90 0 0 0
    90 1 0 0
    90 2 0 0
    90 3 0 0
    90 4 0 0
    90 5 0 0
    90 6 0 0
    90 7 0 0
    90 8 0 0
    90 9 0 0
    90 10 0 0
    90 11 0 0
    90 12 0 0
  
Since every frame is represented by a single line, this format is quite verbose,
requiring 60 lines per second for 60 fps video - you obviously don't want to
hand-code such files, but rather generate them with software. It's a bit like
slicer format for 3D printers, which you also don't write manually.
The format might be beefed up (e.g. 'progressing' variables, switch of source
image and other parameterization) - the current status quo is just to show
that it can be done and to check that anti-aliasing works as expected - to see
the effect of aliasing, still images aren't usually sufficient.

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

    1 - use b-spline interpolation directly on the source image. this
        is the fastest option, and unless there is a significant scale
        change involved, the output should be 'good enough' for most still
        image renditions. This is the default, with a spline degree of 1,
        a.k.a. bilinear interpolation. Other spline degrees can be chosen
        by passing --degree D, where D is the degree of the spline. Per
        default, splines with degree > 1 are 'prefiltered'. You may pass
        --prefilter D to apply a prefilter for a different degree than the
        one used for evaluation. This can be used e.g. to blur the output
        (use a smaller value for the prefilter degree than for the spline
        degree). Note that b-splines may 'overshoot', unless you omit
        the prefilter. This depends on the signal - if it's sufficiently
        band-limited (nothing above half Nyquist frequency), the spline
        won't overshoot. Using a degree-2 b-spline without prefilter
        (--degree 2 --prefilter 0) introduces only slight blur and won't
        overshoot - it's often a good compromise and also avoids the
        star-shaped artifacts typical of degree-1 b-splines when the
        signal is magnified a lot.

    -1 - use OpenImageIO's 'environemnt' or 'texture' function for lookup.
         without additional arguments, this will use a sophisticated
         interpolator with good antialiasing and bicubic interpolation,
         but without prefiltering - the output will look slightly blurred.

    -2 - use 'twining' - this is a method which first super-samples and then
         combines several pixels to one output pixel ('binning'). This is my
         own invention. It's quite fast and produces good quality output.
         This method should see community review to compare it with other
         methods. If you only pass --itp -2 and not --twine, envutil will
         set up twining parameters calculated to fit well with the
         transformation at hand, and that's also done for image sequence
         output, where the parameters have to adapt to the changing geometry.
         The 'twining' interpolator is 'grafted' onto the 'substrate'
         interpolator - that is the b-spline from the source image, which
         you parameterize just as for --itp 1. So if you pass --degree 3,
         the 'substrate' of the twining operator will be a cubic b-spline,
         rather than a degree-1 b-spline (a.k.a bilinear interpolation).

In general, producing visible output is often a two-step process. The first
step is to provide some sort of 'ground truth' - an internal representation
of the image data which will provide a specific value for a specific pick-up
location. This step tends to aim for speed and precision, without taking
into account considerations like aliasing or artifacts introduced by the
interpolation. The signal which is provided by the first step is also
continuous due to interpolation. Sometimes, the first stage will not use
true interpolation, meaning that the value of the first-stage signal is
not necessarily equal to the image data at discrete coordinates. If so,
the signal is typically blurred - e.g. by using a b-spline kernel on the
raw data without adequate prefiltering.

The second step - if present - operates on the first-stage signal. This step
is often added to avoid problems from using the first-stage signal directly.
The most important effect which the second stage tries to produce is
anti-aliasing. If the output is a scaled-down version of the input, direct
interpolation at the pick-up points will produce aliasing where the input
signal has high-frequency content. A typical strategy to avoid this is to
consider a section of the first-stage signal corresponding to the 'footprint'
of each target pixel and form a weighted sum over source pixels in that area.
The precise way of how this is done varies - OIIO can use an elliptic shape
placed over the first-stage signal and produces an average over this area,
whereas 'twining' gathers several point-samples in this area and forms a
weighted sum. Both methods produce similar results.

A third strategy adds scaled-down versions of the image data, which are then
used to provide 'ground truth' for rendering jobs which would down-scale
from the original resolution. The archetypal construct for this strategy is
an 'image pyramid' - a set of images where each has half the resolution
of the one 'below' it. OIIO can use this strategy (look for mip-mapping),
and it's advantage is that pickup kernels can remain relatively small,
because there are fewer pixels to form a sum over to avoid aliasing.
Other rendering schemes in envutil do not currently use image pyramids;
twining deals with the aliasing problem by picking adequately large
kernels directly on the first-stage signal. It would be feasible, though,
to add pyramid schemes. One problem with image pyramids is the fact that
they have to be produced at all: they take up memory and generating them
costs computational resources. envutil often does 'one-shot' jobs, and
rather than producing an image pyramid with *all* resolutions, producing
output with just one resolution - even if that requires a large kernel -
is usually more efficient. If an image pyramid is present, it's an option
to mix 'ground truth' data from two adjacent pyramid levels for output
whose resolution is between that of the two pyramid levels (look e.g.
for 'trilinear interpolation'). twining, which uses a scalable filter,
can adapt the filter size to the change in resolution, so it doesn't use
this method.

Why use negative values for ITP for the second and third mode? This is similar
to the values used in lux for 'decimators' - in lux, positive values are
reserved for degrees of a b-spline reconstruction filter used as low-pass
filter. bilinear interpolation is the same as a degree-1 b-spline, hence the
value 1 for bilinear interpolation. In envutil, I have decided against using
positive 'itp' values from 'meaning' b-splines -here, you pass --itp 1 for
all splines, and --degree and, optionally, --prefilter additionally condition
the spline:

##  --degree DEG      degree of the b-spline (>= 0)
##  --prefilter DEG   prefilter degree (>= 0) for the b-spline

Images created with itp set to 1 or -2 both use b-splines as their substrate.
If you don't pass --degree or --prefilter, a degree-1 b-spline is used - this
is also known as 'bilinear interpolation' and already 'quite good'. But higher
spline degrees can produce even better output, especially if the view magnifies
the source image, and the star-shaped artifacts of the bilinear interpolation
become visible. If you only pass 'degree', the spline will be set up as an
interpolating spline, meaning it will yield precisely the same pixel values
as the input when evaluated at discrete coordinates. This requires prefiltering
of the spline coefficients with a prefilter of the same degree as the spline.
If you pass a different prefilter degree, the coefficients are prefiltered
*as if* the spline degree were so, whereas the evaluation is done with the
given degree. You can use this to produce smoother output (lower prefilter
degree than spline degree) or to sharpen them (higher prefilter degree than
spline degree). A disadvantage of interpolating splines is that they will
produce ringing artifacts if the input signal isn't band-limited to half the
Nyquist frequency. With raw images straight from the camera this is usually
the case, but processed or generated images often have high frequency content
which will result in these artifacts, which become annoyingly visible at high
magnifications. One way to deal with this problem is to accept a certain amount
of smoothing by omitting the prefilter - so, passing --prefilter 0 and --degree
greater than one. With a spline degree of two and no prefiltering, there is
mild suppression of high frequencies, but there are no ringing artifacts, and
this is often a good compromise.

# Twining-specific options

These options control the 'twining' filter. These options only have an effect if
you activate twining with --itp -2. envutil will use b-spline interpolation on
the source image for single-point lookups, and it will perform more lookups and
then combine several neighbouring pixels from the oversampled result into each
target pixel. See the documentation on 'itp' for more details about the
parameterization of the b-spline - the default is to use a degree-1 b-spline
(a.k.a bilinear interpolation).

The operation of the twining filter differs conceptually from OIIO's pick-up
with derivatives: OIIO's filter (as I understand it) looks at the difference
of the 2D coordinates into the source texture and produces a filter with a
footprint proportional to that. twining instead inspects the difference in
3D coordinates and places it's sub-pickups to coincide with rays which are
produced by processing this difference. This is a subtle difference, but it
has an important consequence: OIIO's pick-up will encompass a sufficiently
large area in the source texture (or one of it's mip levels) to generate an
output pixel, so even if this area is very large, the pick-up will be
adequately filtered. twining will spread out it's sub-pickups according to
the number of kernel coefficients and kernel size, but the number of kernel
coefficients does not vary with the pick-up location, only the area over
which the sub-pickups are spread will vary. So to filter correctly with a
twining filter (obeying the sampling theorem), the twining filter needs to
have sufficiently many coefficients spread 'quite evenly' over the pickup
area to avoid sampling at lower rates than the source texture's sampling
rate. While OIIO's filter forms a weighted sum over pixels in one of the
texture's mip levels, twining relies on a b-spline interpolation of the
originally-sized texture as it's substrate. If the sub-pickup locations are
'too close' to each other, the shortcomings of a bilinear interpolation
may 'shine through' if the transformation magnifies the image (you'll see
the typical star-shaped artifacts). If the sub-pick-ups are 'too far apart',
you may notice aliasing - of course depending on the input's spectrum as well.

Keeping this in mind, if you want to set up a twining filter yourself, either
by passing twining-related parameters setting the number and 'spread' of
the filter coefficients or by passing a file with filter coefficients, you
must be careful to operate within these constraints - using automatic
twining, on the other hand (just pass --itp -2 and no other twining-related
arguments) will figure out a parameter set which should keep the filter
within the constraints, so you neither get aliasing for down-scaling views
nor bilinear artifacts in up-scaling views. When creating a twining filter
externally, be generous: processing is fast even with many sub-pickups,
so using a larger number than strictly necessary won't do much harm.

## --twf_file TWF_FILE   read twining filter from a file

Passing a twf file with this option reads the twining filter from this file,
and ignores most other twining-related options. The file is a simple text file
with three float values per line, let's call them x, y and w.

The x and y values are offsets from the pick-up location and w is a weight.
The x and y values are offsets in the horizontal/vertical of the target
image (!): an offset of (1,0) will pick up from the same location as the
unmodified source coordinate of the target pixel one step to the right.
How far away - in source texture coordinates - this is, depends on the
relative geometry - projection and size of the source and target image.

twining is based on an approximation of the derivatives of the coordinate
transformation at the locus of interpolation: if you calculate the pick-up
coordinate for a pixel one step to the right or one step down, respectively,
you can form the differences to the actual pick-up coordinate and obtainin
two vectors xv and yv which approximate the derivative. x and y are used as
multipliers for these vectors, resulting in offsets which are added to the
actual pick-up coordinate for each set of x, y and w. The offsetted coordinate
is used to obtain a pixel value, this value is multiplied with the weight, w,
and all such products are summed up to form the final result. The vectors xv
and yv are 3D vectors and represent the difference of two 3D rays: the ray to
the central pickup location and a ray to another sub-pickup location in it's
vicinity.

A twining filter read from a file will apply the given weights as they are,
there is no automatic normalization - pass 'twine_normalize' to have the
filter values normalized after they have been read. Beware: using a filter
with gain greater than 1 may produce invalid output! The only other parameter
affecting a twining filter read from a file is 'twine_width': if it is not
1.0 (the default), it will be applied as a multiplicative factor to the 'x'
and 'y' values, modifying the filter's footprint.

With an externally-generated twining filter, there is no fixed geometry,
and any kind of non-recursive filter can be realized. envutil applies this
filter to 3D ray coordinates rather than to the 2D coordinates which occur
after projecting the 3D ray coordinate to a source image coordinate, so
the effect is similar to what you get from OIIO's 'environment' lookup
with the elliptic filter.

Because the two vectors, xv and yv, are recalculated for every contributing
coordinate, the filer adapts to the location, and handles varying local
geometry within it's capability: the sampling theorem still has to be
obeyed insofar as generated sub-pickups in a cluster mustn't spread out
too far and when they are too close, shortcomings of the underlying
interpolation may show. And there have to be sufficiently many
sub-pickups. Automatic twining produces a 'reasonable' filter which
tries to address all these issues - play with automatic twining and
various magnifications to see what envutil comes up with (pass -v
to see the filter coefficients).

## --twine_normalize   normalize twining filter weights gleaned from a file

If you pass a file with twining kernel values, the default is to use
them as they are: you are free to apply any possible twining filter,
including, e.g. differentiators which may yield negative output. Most of
the time, though, your kernel's weights will all be positive and the
filter will have unit gain, meaning that all weights add up to one.
At times it's easier to just calculate filter weights without making
sure that the sum of the weights is one. Pass --twine_normalize to let
envutil do the normalization for you.

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

## --twine_width TWINE_WIDTH  alter the size of the pick-up area of the twining filter

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
On the other hand, passing a twine_width smaller than 1.0 will make the
kernel proportionally smaller, producing 'sharper' output, at the cost of
potentially producing more aliasing.

twine_width is the only one of the twining-related options which also
affect a twining filter read from a file - apart from --twf_file itself,
obviously. The twine_width is applied as a multiplicative factor to the
first two parameters of each coefficient, so the default of 1.0 results
in the filter geometry being unchanged.

Input with low resolution is often insufficiently band-limited which will
result in artifacts in the output or become very blurred when you try to
counteract the artifacts with excessive blurring. There's little to be
gained from scaling up anyway - the lost detail can't be regained.
Scaling up with a b-spline of degree two or higher is a good idea, but
to see why this is so, a bit of explanation is needed. Most image viewers
produce up-scaled views by rendering the source image's pixels as small
rectangles, producing a 'pixelated' view. A b-spline rendition is smooth,
though. This leads to the common misconception that the b-spline is
somehow 'blurred' compared to the 'pixelated' rendition. If the b-spline
is made from properly prefiltered coefficients, there is no blurring:
the smoothness is just so much as to adequately represent the uncertainty
created by sampling the ambient information. The 'pixelated' look shows
you the size of the pixels, but a lot of it's visual quality consists
of artifacts: namely the distinct rectangular shapes of the 'pixels'
with discontinuities where one pixel borders on the next one. These
artifacts are often the most prominent part of a 'pixelated' image,
and users are often surprised that the seemingly 'blurred' image from
a b-spline rendition is easier to decipher than the 'pixelated' one.
In fact, if 'pixelation' is used to make faces unrecognizable and the
size of the 'pixels' is not chosen large enough, *smoothing* the
'pixelated' image may reveal a recognizable face, because the information
is no longer 'spoiled' by the high-frequency artifacts resulting from
the artifical rectangles in the pixelated area. Only if prefiltering
is omitted blurring does indeed happen. envutil decouples the prefilter
and spline degree, so you can observe the effects of picking each of
them individually.

If you don't get good results with twining and adequate twine_width,
you may be better off with one of the better OIIO interpolators, e.g.
bicubic. If the signal is sufficiently band-limited, using a b-spline with
a degree of two or higher is a good idea: there are no more star-shaped
artifacts, and for magnifying views, evaluating the spline yields a good
result without further processing. If the signal is not band-limited, but
you can live with slight blurring, a degree-2 b-spline with no prefiltering
(--prefilter 0) is often a good compromise, especially if the resolution
of the input is 'excessive', like the sensor data from very small sensors
and very high resolution - and this very high resolution is not actually
exploited (or usable) to provide fine detail. AFAICT, OIIO's spline-based
interpolators (e.g. bicubic) don't use prefiltering and therefore produce
noticeable blurring. Try a b-spline with proper prefiltering and no twining
(--itp 1 -degree 3) to see the difference - you'll see no blurring, but if
the view magnifies and the signal isn't band-limited, you may see ringing
artifacts.

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

## --twine_precise      project twining basis vectors to tangent plane

This parameter affects all twining filters, but it's effect is rarely
noticeable. To explain what this parameter does, a bit of explanation
is needed. When envutil does a lookup from an environment, it starts out
with the target coordinate - the discrete coordinate where an output
pixel will appear. Next, it figures out a 3D 'ray' coordinate which
encodes the direction where the lookup should 'look for' content. This
obviously depends on the target projection and the orientation of the
virtual camera. With 'simple' interpolation - like the bilinear
interpolation you get with --itp 1, the next step is to find a 2D
coordinate into the source image which has the content corresponding
to the given ray and interpolate a pixel value from the source image
at that coordinate.

The more sophisticated lookup methods in envutil use 'derivatives':
they don't merely look at the single lookup point, but take into
account what happens when the lookup progresses to the next neighbours
of the target point, both in the horizontal and the vertical: the
neighbours have different corresponding rays and therefore different
corresponding 2D source image coordinates. A simple way to obtain the
derivative is to subtract the current pick-up coordinate from those
corresponding to it's neighbours, yielding two 2D differences which
encode an approximation of the derivative. But there is also a
different way: the difference can be formed between the 3D ray
coordinates, and this yields two 3D differences. These differences
are (usually small) vectors from one point on the sphere to another.
envutil's twining filter uses these two vectors to place the pick-up
coordinates for each kernel coefficient: it multiplies the first one
with the 'x' value of the kernel coefficient, the second one with the
'y' value', and adds the result to the current pickup ray coordinate,
yielding a slightly altered ray which is used as sub-pick-up. But there
is a problem here: The two small difference vectors are not perpendicular
to the current pick-up ray, because they are not on it's tangent plane,
but instead connect two points on a spherical surface. If the difference
is small, they are still almost parallel to the tangent plane and this
distinction is irrelevant, but to be very precise, and for larger
differences, it's better to project the difference vectors onto the
tangent plane and use the resulting vectors to calculate the sub-pick-up
rays. Without this projection to the tangent plane, the pattern of
3D ray coordinates is tilted slightly off the tangent plane, and a
sub-pick-up with negative kernel x and y values is not at precisely
the same distance from the current pick-up location as one with positive
x and y of equal magnitude. --twine_precise takes care of this problem
and does the projection to the tangent sphere. The calculation needs to
be done once per target pixel, so with increasing kernel size it does
not require additional computations. Most of the time the change which
results from this parameter is so minimal that it's not worth the extra
effort, that's why it's not the default.

So what's the point of calculating the 'derivative' as a 3D vector, vs.
a 2D difference in source image coordinates? Suppose you have a sub-pick
whose ray is substantially different from the current pick-up location.
What if it's projection to the source image plane lands outside the source
image, or on the opposite edge (which can happen e.g. with full
spherical source images)? To deal with this situation, you need to add
code which recognizes and mends this problem, and you need code for each
type of projection to handle such issues properly. Because this is all
inner-loop code and introduces conditionals, it's detrimental to performance,
where you want conditional-free stencil code which is the same for each
coordinate. Calculating the 'derivatives' in 3D avoids this problem:
only the 3D ray coordinates of each sub-pick-up are projected to source
image coordinates and the current pick-up coordinate is never 'taken
all the way' to the source image - the result pixel is a weighted sum
of the sub-pick-ups. The sub-pick-ups either yield a pixel value or they
don't, the 2D difference never occurs and therefore doesn't produce any
problems.

## --twine_density DENSITY increase tap count of an 'automatic' twining filter

'automatic' twining uses a number of twining kernel coefficients which
is deemed sufficient to produce an artifact-free result (or to keep the
unavoidable artifacts so small that they won't be noticeable). You may
want to be 'more generous' and increase the number of kernel coefficients,
and twine_density does that: if your twining kernel is generated
automatically, it multiplies the 'twine' value with this factor, rounds,
and assigns the result to 'twine'.

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

# Parameters for mounted image input

When --mount or --facet are given, this overrides any --input parameter
you may have passed - you wouldn't, usually, but leave 'input' unset if
you use a mounted image as input.

## --mount IMAGE PROJECTION HFOV  load non-environment source image

envutil can 'mount' images in various projections and hfov which may only cover
a part of the full 360X180 degree environment. This routes to different code,
because the parts of the output which don't receive input have to be masked
out, and because projections apart from spherical (which is present for
lat/lon environments anyway) have to be dealt with. Supported projections
for mounted images are 'rectilinear', 'spherical', 'cylindrical',
'stereographic'  and 'fisheye'. hfov is in degrees. All three values
(image filename, projection and hfov) must be passed after --mount, separated
by space. Currently, there is an issue with itp -1 (OIIO pick-up) with the
default settings when the mounted image is a 360 degree fisheye - some artifacts
show near the 'back pole'. Use other interpolators, or disable the code
relying on derivatives (use --stwidth 0).

## --facet IMAGE PROJECTION HFOV YAW PITCH ROLL
##     load oriented non-environment source image

This is an extension to 'mounted images', taking three more parameters: these
are Euler angles specifying the orientation of the mounted image - 'plain'
mounted images are mounted 'straight ahead'. If you specify three angles
y, p and r, then an orientation of the virtual camera with the same Euler
angles will cancel the orientation precisely. It's intended to allow several
fcets, but this isn't implemented yet.

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
system code does so), so this is merely a technical hint. If you use b-spline
interpolation with full-spherical source images, a specialized prefilter function
is used which is correctly set up to handle both dimensions as periodic signals,
producing artifact-free evaluation near the poles.

As an alternative to the antialiasing and interpolation provided by OIIO, envutil
offers processing with b-spline interpolation and it's own antialiasing filter, using
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
When using b-splines 'better' than degree-1 (a.k.a. bilinear interpolation),
upscaling should be artifact-free - apart from possible ringing artifacts if
the signal isn't band-limited to half the Nyquist frequency. The ringing
artifacts can be avoided by omitting the prefilter, which will produce a
smoothed and bounded signal, but remove some high frequency content (if there
is any).

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

# Style

I use a coding style which avoids forward declarations, so at times it's
better to read the code from bottom to top to follow the flow of control.
My code formatting is quite 'old-school', with white-space-separated
tokens and curly braces in separate lines - and generally a lot of white
space. I think this style makes the code easier to grasp and also helps
with reading it.

The code itself is 'optimistic' - I don't make efforts to prevent errors
from happening - if they happen, most of the time an exception will terminate
the program. I tend not to analyze parameters and trust the user to pass
sensible ones.

The code is complex - it's template metacode, and control flow may be hard
to grasp, due to the use of functional programming. Functional components
follow the pattern I have established in vspline and zimt: the functors
have an 'eval' method which takes a const reference for input and writes
output to another reference. The eval method itself is void. Usually
I code it as a template to avoid a verbose signature - the types are
implicit from the functor's template arguments, and oftentimes the template
can be used for scalar and SIMDized arguments alike.

# An Aside - Image Pyramids Reloaded

Image pyramids are a popular method of providing multi-resolution
representations of images. The traditional method to generate an image
pyramid is to apply a low-pass filter to the image, then decimate it
to obtain the next level. This level is again low-pass-filtered and
decimated, and this is repeated until only one or a few pixels are left.

With a scalable filter (like twining) there is an alternative approach,
which I find promising and which has the potential to produce better
pyramids: You start out with a twining kernel with many coefficients
(say, 16X16) and directly produce the first few pyramid levels by
producing a scaled-down output from the original (level-0) image,
halving the size with each rendition and at the same time doubling the
kernel width (by passing a doubled 'twine_width' parameter). Depending
on the kernel, you may want to start with a twine_width less than one
for the first rendition and start doubling from that - e.g .125, .25 ...
Finally, you reach a point where another doubling of the kernel width
would violate the sampling theorem (kernel coefficients further apart
than the source image's sampling step), and *only then* you start using
the next level of the pyramid as input, now keeping twine_width constant
for successive renditions, rather than carrying on with the doubling
of the kernel width. By using a large kernel on a level 'a few steps
down' you obtain better quality pyramid levels, because they suffer less
from the inevitable degradation due to decimation, which necessarily
introduces some aliasing, because the low-pass never successfully
eliminates the upper half of the frequency band completely. Especially
if one limits the choice of filter to all-positive, there is a good deal
of high frequency left in the result which will alias when the signal
is decimated. Using the scheme described above, by the time you 'reach
down' to a level 'a few steps down', no intermediate decimations have
taken place, so the only effect is what the twining filter does to the
signal - and you can shape that kernel any way you like.
