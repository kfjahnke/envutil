# envutil

utility to convert between lat/lon and cubemap environment maps

This is a stnd-alone repository for envutil, which started out as a demo program
for my library [zimt](https://github.com/kfjahnke/zimt). The program has grown
beyond the limits of what I think is sensible for a demo program, and I also
think it's a useful tool.

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
    cmake [options] ...
    make

this should produce a binary named 'envutil' or 'envutil.exe'.

envutil --help gives a summary of command line options:

    $ envutil --help
    
    envutil -- convert between lat/lon and cubemap format

    Usage: envutil [options] --input INPUT --output OUTPUT

        --help                Print help message
        -v                    Verbose output
        --input INPUT         input file name (mandatory)
        --output OUTPUT       output file name (mandatory)
        --save_ir INTERNAL    save IR image to this file
        --ts_options OPTIONS  pass comma-separates k=v list of options to
                                OIIO's texture system
        --extent EXTENT       width of the cubemap / height of the envmap
        --itp ITP             interpolator: 1 for direct bilinear, -1 for
                                OIIO's anisotropic
        --twine TWINE         use twine*twine oversampling and box filter -
                                best with itp1
        --twine_width TWINE_WIDTH  widen the pick-up area of the twining filter
        --twine_sigma TWINE_SIGMA  use a truncated gaussian for the twining
                                         filter (default: don't)
        --twine_threshold TWINE_THRESHOLD  discard twining filter taps below this threshold
        --face_fov FOV        field of view of the cube faces of a cubemap
                                input (in degrees)
        --support_min EXTENT  minimal additional support around the cube face
                                proper
        --tile_width EXTENT   tile width for the internal representation image
        --ctc                 flag indicating fov is measured between marginal
                                pixel centers
    
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
map, and any image file with 1:6 aspect ratio will be accepted as a cubemap. The
only requirement for cubemaps is that their width should be even.

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
the values.

## --ts_options OPTIONS  pass comma-separates k=v list of options to OIIO's texture system

You can do this, but since we're just running a single image to single image operation,
we needn't really pass any texture system argumets.

## --extent EXTENT       width of the cubemap / height of the envmap

This value is for the output. The extent of the input is gleaned from the file.
If you don't pass an extent, envutil will pick a 'sensible' value: if the source
is a lat/lon environment, it will pick a multiple of 64 which is just larger than
4/pi times the lat/lon input's width, to make up for the distortions of the
rectilinear projection towards the edges. For the reverse conversion, it will
use four times the cubemap's width, which results in at least equal resolution.

## --itp ITP             interpolator: 1 for direct bilinear, -1 for OIIO's anisotropic

So far, there are these two interpolators. the 'direct blinear' is very fast and it's often good enough if the resolution of source and target don't differ much. If aliasing
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

## --face_fov FOV        field of view of the cube faces of a cubemap input (in degrees)

If an input cubemap has other than ninety degrees field of view per cubeface image, pass
this value here.

## --support_min EXTENT  minimal additional support around the cube face proper

This is a technical value, probably best to leave it at it's default of 4.

## --tile_width EXTENT   tile width for the internal representation image

Also best left at the default of 64.

## --ctc                 flag indicating fov is measured between marginal pixel centers

Some cubemaps have marginal pixels which coincide geometrically with the virtual cube's
edges. For such cubemaps, set this flag. The default is to consider the pixels evenly
distributed so that the marginal samples are half a sampling step from the virtual cube's
(and the cube face image's) edge.

# Technical Notes

One problem with cubemaps is that they are normally stored as concatenations of
six square images with precisely ninety degrees fov (field of view). This makes
access to pixels near the edges tricky - if the interpolator used for the purpose
needs support, marginal pixels don't have them all around, and it has to be gleaned
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
whether it needs to be interpolated or antialiased.

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
