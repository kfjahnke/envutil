# envutil

utility to convert between lat/lon and cubemap environment maps

This is a stnd-alone repository for envutil, which started out as a demo program
for my library [zimt](https://github.com/kfjahnke/zimt). The program has grown
beyond the limits of what I think is sensible for a demo program, and I also
think it's a useful tool.

This initial version will, per default, not use any of the SIMD back-ends which
might be employed for better SIMD code but instead rely on small-loop autovectorization
using 'zimt's own' SIMD emulation back-end. Since performance is largely I/O-bound,
this is not a big drawback - eventually I'll add support for the other back-ends
as well. The only dependency is [OpenImageIO](https://github.com/AcademySoftwareFoundation/OpenImageIO).

The program is built with CMake. To build, try this:

    mkdir build
    cd build
    cmake ..
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
        --twine ITP           use twine*twine oversampling and box filter -
                                best with itp1
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

Any image file which has 2:1 aspect ratio will be accepted as lat/lon environment
map, and any image file with 1:6 aspect ratio will be accepted as a cebemap. The
only requirement for cubemaps is that their width should be even.

envutil only processes sRGB and linear RGB data, the output will be in the same
colour space as the input. If you use the same format for input and output, this
will automatically be the case, if not, you may get faulty output if the default
coulour spaces of the formats don't match - your output will look too bright or
too dark.
