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

This program takes a 2:1 lat/lon environment or a 1:6 cubemap image as input and
produces output in the specified orientation, projection, field of view and extent. The program is new and will still need some tweaking, but it's
already functional.
This started out as a simple demo for 'steppers' (stepper.cc is still there with
the initial code), but I thought that, with a bit of additional parameterization,
it would make a useful tool. For the time being, it only uses bilinear
interpolation, so the resolution of the output should be close to the input's.
The output projection can be one of "spherical", "cylindrical", "rectilinear",
"stereographic", "fisheye" or "cubemap". The geometrical extent of the output is
set up most conveniently by passing --hfov, the horizontal field of view of the
output. The x0, x1, y0, and y1 parameters allow passing specific extent values
(in model space units), which should rarely be necessary. To specify a 3D
rotation, pass Euler angles yaw, pitch and roll - they default to zero: no
rotation. The size of the output is given by --width and --height. You must pass
an output filename with --output; --input specifies the environment image.

## --input INPUT         input file name (mandatory)

Any image file which has 2:1 aspect ratio will be accepted as lat/lon environment
map, and any image file with 1:6 aspect ratio will be accepted as a cubemap.
See --6 for a way to load a cubemap from six separate cube face image files. 

## --output OUTPUT       output file name (mandatory)

The output will be stored under this name. See --6 for a way to store cubemaps
to six separate cube face image files.

## --width EXTENT    width of the output

in pixel units. For cubemaps, this should be precisely six times the height.

## --height EXTENT   height of the output

in pixel units.

## --projection PRJ  target projection

Pass on of the supported output projections: "spherical", "cylindrical",
"rectilinear", "stereographic", "fisheye" or "cubemap". The default is
"rectilinear".

## --hfov ANGLE      horiziontal field of view of the output (in degrees)

For cubemap output, you want precisely 90 degrees - or a bit more, if your
cubemap processing software can handle such cubemaps. The 'natural' limit for
rectilinear images is 180 degrees, but this does not produce 'sensible' output.
Same for stereographic images. spherical, cylindrical and fisheye output can
have any field of view - the environment image will be periodized for hfov
larger than 360 degrees.

## --yaw ANGLE       yaw of the virtual camera (in degrees)
## --pitch ANGLE     pitch of the virtual camera (in degrees)
## --roll ANGLE      roll of the virtual camera (in degrees)

These three angles are applied to the 'virtual camera' taking the extracted
view. They default to zero. It's okay to pass just one or two.

## --x0 EXTENT       low end of the horizontal range
## --x1 EXTENT       high end of the horizontal range
## --y0 EXTENT       low end of the vertical range
## --y1 EXTENT       high end of the vertical range

These are special values which can be used to specify the extent, in model space
units, of the output. This requires some understanding of the inner workings of
this program - if you use -v, the verbose output will tell you for each extraction
which extent values are generated from a field of view parameter, given a specific
projection. This can help you figure out specific values you may want to pass,
e.g. to produce anisotropic output.


