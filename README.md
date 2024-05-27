# This repository has the code for [envutil](#envutil) and [extract](#extract)

# envutil
# utility to convert between lat/lon and cubemap environment maps

This is a stand-alone repository for envutil, which started out as a demo program
for my library [zimt](https://github.com/kfjahnke/zimt). The program has grown
beyond the limits of what I think is sensible for a demo program, and I also
think it's a useful tool. There is now a second useful program in this repository,
called 'extract'. It will also be built and installed by the cmake script. See further down for documentation.

The program is built with CMake.

The only mandatory dependency is [OpenImageIO](https://github.com/AcademySoftwareFoundation/OpenImageIO) - OIIO in short.
If you have highway installed, the build will use it, if not it tries for Vc,
and if that isn't present either it uses std::simd. If you dont' want any of
the three, set the option USE_GOADING=ON. The functionality and results should
be the same with either of the back-ends, but processing time will vary,
especially if you are using 'twining'.

It's recommended to build with clang++ (call cmake with -DCMAKE_CXX_COMPILER=clang++).

To build, try this:

    mkdir build
    cd build
    cmake [options] ..
    make

this should produce a binary named 'envutil' or 'envutil.exe', same for 'extract'.

envutil --help gives a summary of command line options:

    $ envutil --help
    
    envutil -- convert between lat/lon and cubemap format

    Usage: envutil [options] --input INPUT --output OUTPUT

    ## --help                Print help message
        -v                    Verbose output
    ## --input INPUT         input file name (mandatory)
    ## --output OUTPUT       output file name (mandatory)
    ## --save_ir INTERNAL    save IR image to this file
    ## --ts_options OPTIONS  pass comma-separates k=v list of options to
                                OIIO's texture system
    ## --extent EXTENT       width of the cubemap / height of the envmap
    ## --itp ITP             interpolator: 1 for direct bilinear, -1 for
                                OIIO's anisotropic
    ## --twine TWINE         use twine*twine oversampling and box filter -
                                best with itp1
    ## --twine_width TWINE_WIDTH  widen the pick-up area of the twining filter
    ## --twine_sigma TWINE_SIGMA  use a truncated gaussian for the twining
                                         filter (default: don't)
    ## --twine_threshold TWINE_THRESHOLD  discard twining filter taps below this threshold
    ## --face_fov FOV        field of view of the cube faces of a cubemap
                                input (in degrees)
    ## --support_min EXTENT  minimal additional support around the cube face
                                proper
    ## --tile_width EXTENT   tile width for the internal representation image
    ## --ctc                 flag indicating fov is measured between marginal
                                pixel centers
    ## --6                   use six separate cube face images
    ## --lux                 use cube face images named in lux convention
    
The specific conversion will depend on the input (given with --input ...).
If you pass a lat/lon environment map, it will be converted to a cubemap, and
if you pass a cubemap, it will be converted to a lat/lon environment map. The
result of the conversion is stored under the file name given by --output ...
and the input is not affected.

envutil only processes sRGB and linear RGB data, the output will be in the same
colour space as the input. If you use the same format for input and output, this
will automatically be the case, if not, you may get faulty output if the default
coulour spaces of the formats don't match - your output will look too bright or
too dark.

## --input INPUT         input file name (mandatory)

Any image file which has 2:1 aspect ratio will be accepted as lat/lon environment
map, and any image file with 1:6 aspect ratio will be accepted as a cubemap.
See --6 for a way to load a cubemap from six separate cube face image files. 

## --output OUTPUT       output file name (mandatory)

The output will be stored under this name. See --6 for a way to store cubemaps
to six separate cube face image files.

## --save_ir INTERNAL    save IR image to this file

envutil holds an internal representation (IR) of the cubemap, which augments each
cube face image with a frame of additional pixels, producing a vertical concatenation
of six square images with as large or larger field of view than the cubemap's
representation on disk. The IR is also sized to be a multiple of a given tile
width, to make it easy to mipmap it. To save the IR to an image file, use this
option. The IR image can serve as input to envutil, but you need to tell envutil
the field of view, which is typically larger than ninety degrees. Since you
already have the frame of support and a size which is a multiple of the tile
size, you'd pass --min_support 0 and --tile_size 1 - unless you want to change
the values. Using an IR file as input will only work for even cube face image sizes.
For those, when you run envutil with -v, it will give you the fov for the IR
image.

## --ts_options OPTIONS  pass comma-separates k=v list of options to OIIO's texture system

You can do this, but since we're just running a single image to single image operation,
we needn't really pass any texture system argumets.

## --extent EXTENT       width of the cubemap / height of the envmap

This value is for the output. The extent of the input is gleaned from the file.
If you don't pass an extent, envutil will pick a 'sensible' value: if the source
is a lat/lon environment, it will pick a multiple of 64 which is just larger than
1/pi times the lat/lon input's width, to make up for the distortions of the
rectilinear projection towards the edges - the extent is chosen so that the
resolution in the center of the cube faces is the same as the lat/lon image's
at the equator. For the reverse conversion, it will use twice the cubemap's
width, which also results in at least equal resolution.

## --itp ITP             interpolator: 1 for direct bilinear, -1 for OIIO's anisotropic

So far, there are these two interpolators. the 'direct blinear' is very fast and it's
often good enough if the resolution of source and target don't differ much. If aliasing
is an issue and good quality is required, OIIO's anisotropic is recommended. This is
slow and tends to come out quite soft, losing some detail. I offer 'twining' with the
next option, which can improve the quality of the result when the direct bilinear
interpolator is used - it can be used with OIIO's ansotropic, but it won't have much
of an effect there.

## --twine TWINE         use twine*twine oversampling and box filter

The effect of 'twining' is the same as oversampling and subsequent application of a box
filter. The filter is sized so that the oversampling is uniform over the data, but the
direct result of the oversampling is never saved - all samples falling into a common
output pixel are pooled and only the average is stored. This keeps the pipeline afloat
in SIMD registers, which is fast (as is the arithmetiic) - especially when highway or
Vc are used, which increase SIMD performance particularly well.

## --twine_width TWINE_WIDTH  widen the pick-up area of the twining filter

A second parameter affecting 'twining'. If the source image has smaller resolution
than the target image, the output reflects the interpolator's shortcomings, so
with e.g. bilinear interpolation and large scale change (magnification) the output
may show star-shaped and staircase artifacts. To counteract this problem, try and pass
a twine_width up to roughly half the magnitude of the scale change.
Input with low resolution is often insufficiently band-limited which will result
in artifacts in the output or become very blurred when you try to counteract the
artifacts with excessive blurring. There's little to be gained from scaling up
anyway - the lost detail can't be regained.

## --twine_sigma TWINE_SIGMA  use a truncated gaussian for the twining filter (default: don't)

If you don't pass --twine_sigma, emvutil will use a simple box filter to combine the result of supersampling into single output pixels values. If you pass twine_sigma, the kernel will be derived from a gaussian with a sigma equivalent to twine_sigma times the half kernel width. This gives more weight to supersamples near the center of the pick-up.

## --twine_threshold TWINE_THRESHOLD  discard twining filter taps below this threshold

If you pass twine_sigma, marginal twining kernel values may become quite small and using them as filter taps makes no sense. Pass a threshold here to suppress kernel values below the threshold. This is mainly to reduce processing time. Use -v to display the kernel and see which kernel values 'survive' the thresholding.

## --face_fov FOV        field of view of the cube faces (in degrees)

Normally, the cube face images in a cubemap have ninety degrees field of view.
But you can pass a larger value here. If the cubemap is the input, passing
face_fov specifies that this is not a 'standard' cubemap with ninety degrees
field of view, but that it contains images with the given field of view.
If the cubemap is made (like, from the lat/lon environment map) and you
pass face_fov, the output cubemap will have this nonstandard field of view.
Writing cubemaps, envutil will try and set appropriate metadata. Please
refer to the 'ctc' argument as well!

## --support_min EXTENT  minimal additional support around the cube face proper

This is a technical value, probably best to leave it at it's default of 4.

## --tile_width EXTENT   tile width for the internal representation image

Also best left at the default of 64.

## --ctc  flag indicating cube face fov is measured between marginal pixel centers

The standard way of measuring the field of view of cube face images in envutil is
to consider pixels as small square areas of constant colour with an extent of one
pixel unit. If an image is W pixels wide, a field of view of D degrees is taken to
coincide with the angle between rays to the left margin of the leftmost pixel and
the right margin of the rightmost pixel (same for top and botttom). If you pass
--ctc, D will instead coincide with the angle between rays to the centers of the
marginal pixels.
So usually, we have D = atan ( f * W / 2 ), with --ctc D = atan ( f * ( W - 1 ) / 2 ) This is hard to see, but some cubemaps seem to use this convention, and using them without --ctc will lead to subtle errors. Internally, envutil uses the first notion, and simply recalculates the field of view to be used internally to the slightly larger value which results form the edge-to-edge notion. Note that --ctc does not affect
the processing of lat/lon environment maps - I may add a separate option for that.
lat/lon environment maps are processed with edge-to-edge semantics.

## --6                   use six separate cube face images

This flag pertains both to input and output and tells envutil to expect or produce
six separate cube face images. The images have to follow a naming scheme: all six
file names are derived from a common name template. This template is split into
base name and extension, the the base name is suffixed with an underscore followed
by the view direction (left, right, top, bottom, front, back). Note that no file
by the name of the template name is used for output - if such a file exists, this
does not matter to envutil. For input, if such a file exists, envutil will try
and open it as input, because envutil looks at the input to figure out which
conversion is wanted, so the file's presence will prevent the processing of
the six separate cube face images. Also see the next option.

## --lux                 use cube face images named in lux convention

In lux (as of this writing) front and back, and left and right are swapped, so
lux effectively 'looks' to a view rotated by 180 degrees around the vertical.
To create a lux cubemap script file, this has to be taken into account, and
to process cube face images made for processing with lux it's also essential:
otherwise, the top and bottom images won't agree with the rest. Note that
lux currently does not produce entirely flawless views with cube face
images of precisely ninety degrees. Just a little bit of extra foc is enough;
you might use --ctc to the effect.
The comments about the directions holds true for lux up to 1.2.2 - I intend
to change lux to follow openEXR convention, and eventually to use code from
this program for better (and faster) cubemap handling.

# Technical Notes

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
--support_min and --tile_width - usually it's best to leave them at their defaults.
To access pixels in lat/lon environment maps which are marginal, envutil exploits the
inherent periodicity of the lat/lon image - simple periodicity in the horizontal and
'over-the-pole' periodicity in the vertical (that's mirroring plus an offset of half
the image's width).

As an alternative to the antialiasing and interpolation provided by OIIO, envutil
offers processing with bilinear interpolation and it's own antialiasing filter, using
a method which I call 'twining'. Twining exploits the fact that the processing builds
up a pixel pipeline as a functional construct, and parts of the chain of functors can
be 'wrapped' in other functors which modify their input or output. So twining picks
the part of the processing chain which handles the conversion from a 2D pick-up
coordinate to a pixel value and wraps it in an outer functor which calls that bit
of processing repeatedly with slightly altered coordinates, gathers the results and
forms a weighted sum from them. The weighted sum is then produced as output of the
outer functor. The remainder of the pixel pipeline remains unaffected, so the
mechanism is independent of other precessing steps and can be used in combination
with all of the interpolators which envutil offers. So using 'twining' together
with OIIO's environment lookup is quite possible, though it's not sensible, because
OIIO provides it's own filtering.

From visual inspection of the results, envutil's use of OIIO's lookup seems to produce
quite soft images. OIIO's lookup is - with the settings envutil uses - conservative
and makes an effort to avoid aliasing, erring on the side of caution. Of course I
can't rule out that my use of OIIO's lookup is not coded correctly, but I'm quite
confident that I've got it right. Using OIIO's lookup is quite involved, because
it needs the derivatives of the pick-up coordinate relative to the progress of
the canonical target image coordinate. As of this writing, OIIO's lookup offers
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
a suitable process by inspecting the output. I use OIIO's lookup as default to
honour it's being well-established and well-thought-out, whereas the use of
bilinear interpolation, optionally with 'twining', still has to prove it's
suitability, and 'twining' is my invention and not yet tested much - I'm sure
I'm not the first one to think of this method, but I haven't done research to
see if I can find similar code 'out there'. What I am sure of is that my
implementation is fast due to the use of multithreaded horizontal SIMD code
through the entire processing chain, so I think it's an attractive offer. I'd
welcome external evaluation of the results and a discussion of the methods;
please don't hesitate to open issues on the issue tracker if you'd like to
discuss or report back!

There seems to be ambiguity of what constitutes a 'correct' cube face image with
ninety degrees field of view. In envutil, I code so that each pixel is taken to
represent a small square section of the image with constant colour. So the
'ninety degrees proper' extends form the leftmost pixel's left margin to the
rightmost pixel's right margin. Some cubemap formats provide images where the
centers of the marginal pixels coincide with the virtual cube's edges, repeating
each edge in the image it joins up with in the cube. If you process such cubemaps
with envutil, pass --ctc (which stands for center-to-center). Otherwise, there will
be subtle errors along the cube face edges which can easily go unnoticed. Make sure
you figure out which 'flavour' your cubemaps are.

# extract
# utility to extract a reprojected image from an environment.

This program takes a 2:1 lat/lon environment or a 1:6 cubemap image as input
and produces output in the specified orientation, projection, field of view
and extent. For CL arguments, try 'extract --help'. Panorama photographers
may not be familiar with the term 'lat/lon environment' - to them, this is
a 'full spherical panorama' or 'full equirect'. The difference is merely
terminology, both are the same. Cubemaps are rarely used in panorama
photography, but I aim at reviving them by building a code base to
process them efficiently, hone the standard, and integrate processing of
single-image cubemap format (as it is used by extract) into [lux](https://bitbucket.org/kfj/pv), my
image and panorama viewer - currently lux cubemaps need six separate images
and a special litle script and the code is less optimized than the code I
offer here. The single-image cubemap format I process in extract leans
on the openEXR standard as far as the order and orientation of the cube
face images is concerned: the images are expected stacked vertically (in
this sequence: left, right, top, bottom, front, back; top and bottom align
with the 'back' image), but the precise geometry of the images may differ
from openEXR's standard - I measure cube face image field of view
'edge-to-edge', meaning that I consider pixels small square areas and
X degrees field of view extend over just as many of these small areas as
the image is wide and high. A different notion is 'center-to-center'
measuring, where the X degrees refer to the area enclosed between the
marginal pixels seen as points in space (the convex hull) - this notion
is not currently supported in extract.

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
can use 'extract' to reproject environment images, e.g convert a cubemap
to a lat/lon image and the reverse. To get a full 360 degree lat/lon
image from a cubemap, you must pass --hfov 360 and --projection spherical,
and for the reverse operation, pass --hfov 90 and --projection cubemap.
The parameterization for cubemaps is not quite as flexible as what
envutil offers (e.g. the --ctc option is missing) but otherwise
extract can do roughly the same job, and I may opt to discontinue
envutil in it's favour. Another feature which extract doesn't currently
offer is processing of six separate single cube face images.

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
binning. The default with this method is to use a simple 2X2 box filter
on a signal which is interpolated with bilinear interpolation, and oversampled by a factor of four. Additional parameters can change the
amount of oversampling and add gaussian weights to the filter
parameters. Twining is quite fast (if the number of filter taps isn't
very large); when down-scaling, the parameter 'twine' should be at
least the same as the scaling factor to avoid aliasing. When upscaling,
larger twining values will slighly soften the output and suppress the
star-shaped artifacts typical for bilinear interpolation. Twining is
new and this is a first approach. The method is intrinsically very
flexible (it's based on a generalization of convolution), and the full
flexibility isn't accessible in 'extract' with the parameterization
as it stands now, but it's already quite useful with the few parameters
I offer.

The program uses [zimt](https://github.com/kfjahnke/zimt) as it's 'strip-mining' and SIMD back-end, and
sets up the pixel pipelines using zimt's functional composition tools.
This allows for terse programming, and the use of the functional
paradigm allows for many features to be freely combined - a property
which is sometimes called 'orthogonality'. What you can't combine in
'extract' is twining and interpolation with OIIO - this is pointless,
because OIIO offers all the anti-aliasing and quality interpolation
one might want, and using twining on top would not improve the
output.

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

If you pass -v, you'll get informations about what is going on internally,
and passing --help will produce this command line argument synopsis:

    --help                             Print help message
    -v                                 Verbose output
    --input INPUT                      input file name (mandatory)
    --output OUTPUT                    output file name (mandatory)
    --itp ITP                          interpolator: 1 for bilinear, -1 for OIIO, -2 bilinear+twining
    --width EXTENT                     width of the output
    --height EXTENT                    height of the output
    --projection PRJ                   target projection
    --hfov ANGLE                       horiziontal field of view of the output (in degrees)
    --yaw ANGLE                        yaw of the virtual camera (in degrees)
    --pitch ANGLE                      pitch of the virtual camera (in degrees)
    --roll ANGLE                       roll of the virtual camera (in degrees)
    --x0 EXTENT                        low end of the horizontal range
    --x1 EXTENT                        high end of the horizontal range
    --y0 EXTENT                        low end of the vertical range
    --y1 EXTENT                        high end of the vertical range
    --twine TWINE                      use twine*twine oversampling - use with itp -2
    --twine_width TWINE_WIDTH             widen the pick-up area of the twining filter
    --twine_sigma TWINE_SIGMA          use a truncated gaussian for the twining filter (default: don't)
    --twine_threshold TWINE_THRESHOLD  discard twining filter taps below this threshold
    --tsoptions KVLIST                 OIIO TextureSystem Options: coma-separated key=value pairs
    --swrap WRAP                       OIIO Texture System swrap mode
    --twrap WRAP                       OIIO Texture System twrap mode
    --mip MIP                          OIIO Texture System mip mode
    --interp INTERP                    OIIO Texture System interp mode
    --stwidth EXTENT                   swidth and twidth OIIO Texture Options
    --stblur EXTENT                    sblur and tblur OIIO Texture Options
    --conservative_filter              OIIO conservative_filter Texture Option


## --input INPUT         input file name (mandatory)

Any image file which has 2:1 aspect ratio will be accepted as lat/lon environment
map, and any image file with 1:6 aspect ratio will be accepted as a cubemap.
See --6 for a way to load a cubemap from six separate cube face image files. 

## --output OUTPUT       output file name (mandatory)

The output will be stored under this name. See --6 for a way to store cubemaps
to six separate cube face image files.

## --itp ITP

This is an integer value determining the interpolation method. There are
currently three modes of interpolation:

    1 - use simple bilinear interpolation directly on the source image
        this is the fastest option, and unless there is a significant scale
        change involved, the output should be 'good enough' for most purposes.
        This is the default.

    -1 - use OpenImageIO's 'environemnt' or 'texture' function for lookup.
         without additional arguments, this will use a sophisticated
         interpolator with good antialiasing and bicubic interpolation.

    -2 - use 'twining' - this is a method which first super-samples and then
         combines several pixels to one output pixel ('binning'). This is my
         own invention. It's quite fast and produces good quality output.
         This method should see community review to compare it with other
         methods. If you only pass --itp -2 and not --twine, extract will
         set up twining parameters calculated to fit well with the
         transformation at hand.

## --width EXTENT    width of the output

in pixel units. For lat/lon environment images, this should be precisely
twice the height.

## --height EXTENT   height of the output

in pixel units. For cubemaps, this must be precisely six times the width.

## --projection PRJ  target projection

Pass on of the supported output projections: "spherical", "cylindrical",
"rectilinear", "stereographic", "fisheye" or "cubemap". The default is
"rectilinear".

## --hfov ANGLE      horiziontal field of view of the output (in degrees)

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
suitable options.

## --yaw ANGLE       yaw of the virtual camera (in degrees)
## --pitch ANGLE     pitch of the virtual camera (in degrees)
## --roll ANGLE      roll of the virtual camera (in degrees)

These three angles are applied to the 'virtual camera' taking the extracted
view. They default to zero. It's okay to pass none or just one or two. yaw
is taken as moving the camera to the right, pitch is taken as upward movement,
and roll as a clockwise rotation. Note that the orientation of the *virtual
camera* is modified; when looking at the resulting images, objects seen on
them seem to move the opposite way. Negative values hove the opposite effect.
Panorama photographers: to extract nadir patches, pass --pitch -90

## --x0 EXTENT       low end of the horizontal range
## --x1 EXTENT       high end of the horizontal range
## --y0 EXTENT       low end of the vertical range
## --y1 EXTENT       high end of the vertical range

These are special values which can be used to specify the extent, in model
space units, of the output. This requires some understanding of the inner
workings of this program - if you use -v, the verbose output will tell you
for each extraction which extent values are generated from a field of view
parameter, given a specific projection. This can help you figure out specific
values you may want to pass, e.g. to produce anisotropic output or cropped
images.

# OIIO-specific options

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
separately, but I don't see a need to do so in envutil and extract. The most
useful value here is to pass zero, which makes OIIO ignore the pickup-point's
derivatives. This can speed up the calculations.

# --stblur EXTENT      sblur and tblur OIIO Texture Options

I also lump together pre-blur along the s and t axis. This is to add more blur
to the processing.

# --conservative_filter YESNO     OIIO conservative_filter Texture Option

pass this to switch on OIIO's 'conservative filter' on or off (pass 0 or 1);
the default is to have it on.

# Twining-specific options

These options control the 'twining' filter. These options only hve an effect if
you activate twining with --itp -2. extract will use bilinear interpolation on
the source image for single-point lookups, but it will perform more lookups and
then combine several neighbouring pixels from the oversampled result into each
target pixel. The meaning of the parameters is the same as in envutil:

## --twine TWINE         use twine*twine oversampling and box filter

The effect of 'twining' is the same as oversampling and subsequent application
of a box filter. The filter is sized so that the oversampling is uniform over
the data, but the direct result of the oversampling is never saved - all
samples falling into a common output pixel are pooled and only the average
is stored. This keeps the pipeline afloat in SIMD registers, which is fast
(as is the arithmetiic) - especially when highway or Vc are used, which
increase SIMD performance particularly well.
If you activate twining by selecting --itp -2, but don't pass --twine (or if
you pass --twine 0 explicitly), extract will set up twining parameters
automatically, so that they fit the relation of input and output. If the
output magnifies (in it's center), the twine width will be widened to
avoid star-shaped artifacts in the output. Otherwise, the twine factor
will be raised to avoid aliasing.

## --twine_width TWINE_WIDTH  widen the pick-up area of the twining filter

A second parameter affecting 'twining'. If the source image has smaller
resolution than the target image, the output reflects the interpolator's
shortcomings, so with e.g. bilinear interpolation and large scale change
(magnification) the output may show star-shaped and staircase artifacts.
To counteract this problem, try and pass a twine_width up to the
magnitude of the scale change. Input with low resolution is often
insufficiently band-limited which will result in artifacts in the output or
become very blurred when you try to counteract the artifacts with excessive
blurring. There's little to be gained from scaling up anyway - the lost
detail can't be regained. If you don't get satisfactory results with
twining and adequate twine_width, you may be better off with one of the
better OIIO interpolators, e.g. bicubic.

## --twine_sigma TWINE_SIGMA  use a truncated gaussian for the twining filter (default: don't)

If you don't pass --twine_sigma, extract will use a simple box filter to combine the result of supersampling into single output pixels values. If you pass twine_sigma, the kernel will be derived from a gaussian with a sigma equivalent to twine_sigma times the half kernel width. This gives more weight to supersamples near the center of the pick-up.

## --twine_threshold TWINE_THRESHOLD  discard twining filter taps below this threshold

If you pass twine_sigma, marginal twining kernel values may become quite small and using them as filter taps makes no sense. Pass a threshold here to suppress kernel values below the threshold. This is mainly to reduce processing time. Use -v to display the kernel and see which kernel values 'survive' the thresholding.

