# envutil

utility to convert between lat/lon and cubemap environment maps

This is a stnd-alone repository for envutil, which started out as a demo program
for my library [zimt](https://github.com/kfjahnke/zimt). The program has grown
beyond the limits of what I think is sensible for a demo program, and I also
think it's a useful tool.

The program is built with CMake.

The only mandatory dependency is [OpenImageIO](https://github.com/AcademySoftwareFoundation/OpenImageIO). If you have highway installed,
the build will use it, if not it tries for Vc, and if that isn't
present either it uses std::simd. If you dont' want any of the three,
set the option USE_GOADING=ON.

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
option.

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
