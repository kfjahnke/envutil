#          cmake file for envutil
#          
#             Kay F. Jahnke <kfjahnke+envutil@gmail.com> 2024
#
#  The git repository for this software is at
#
#  https://github.com/kfjahnke/envutil
#
#  Please direct questions, bug reports, and contributions to
#
#  kfjahnke+envutil@gmail.com
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# You need some dependencies: OpenImageIO and Imath.

# To build, try this for starters (ninja instead of make should also work):

# mkdir build
# cd build
# cmake ..
# make

cmake_minimum_required(VERSION 3.10)

get_filename_component(_varsource "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component(_varbinary "${CMAKE_BINARY_DIR}" REALPATH)

# prevent in-tree building

if("${_varsource}" STREQUAL "${_varbinary}")
	message(FATAL_ERROR "In-source builds are not allowed")
endif()

# # envutil should be built with clang++
# 
# set(ENV{CXX} clang++)
# set(ENV{CC} clang)

set(envutil_major 0)
set(envutil_minor 1)
set(envutil_patch 2)

message(STATUS "VERSION = ${envutil_major}.${envutil_minor}.${envutil_patch}")

# Project name 'envutil'

set(_project_name envutil)
project (${_project_name})

# specify some directories for installation of components. Note that
# some builds fail if CMAKE_INSTALL_PREFIX points to folders with
# white space or non-alphanumeric characters, like 'Program Files (x86)'
# on windows. If that happens, use a 'sane' path like '/usr/local'.

message(STATUS "***** CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}")

# work out the target processor architecture by checking for all values
# of CMAKE_SYSTEM_PROCESSOR which, on any platform, indicates an intel
# processor.

message(STATUS "***** SYSTEM PROCESSOR ${CMAKE_SYSTEM_PROCESSOR}")

if (     ${CMAKE_SYSTEM_PROCESSOR} STREQUAL x86_64
      OR ${CMAKE_SYSTEM_PROCESSOR} STREQUAL AMD64
      OR ${CMAKE_SYSTEM_PROCESSOR} STREQUAL i686 )
  set(i86 TRUE)
  message(STATUS "***** setting i86 TRUE for an intel/AMD target")
else()
  set(i86 FALSE)
  message(STATUS "***** setting i86 FALSE; not an intel/AMD target")
endif()

# rules to create the object files

set(CMAKE_VERBOSE_MAKEFILE ON)

## OIIO anticipates moving to C++17 in 2024

set(CMAKE_CXX_STANDARD 17)
message(STATUS "***** using C++17 standard")

find_package ( OpenImageIO CONFIG REQUIRED )
if ( OpenImageIO_FOUND )
  message ( STATUS "find_package: OpenImageIO = true" )
endif()

find_package ( Imath CONFIG REQUIRED )
if ( Imath_FOUND )
  message ( STATUS "find_package: Imath = true" )
endif()

# highway has added code to find it with cmake config mode with
# the release of 1.0.4. If we can't find highway, we try Vc, and
# if that doesn't work we use std::simd.

option(MULTI_SIMD_ISA "AUTOMATIC SIMD ISA DISPATCH (REQUIRES HIGHWAY)" ON)
option(USE_GOADING "USE SIMD EMULATION WITH SMALL LOOPS" OFF)

if ( USE_GOADING )
  # goading can work together with MULTI_SIMD_ISA, but then we need
  # to have and link to libhwy.
  if ( MULTI_SIMD_ISA )
    find_package ( HWY CONFIG )
    if ( HWY_FOUND )
      message ( STATUS "find_package: HWY_FOUND = true" )
      set ( required_libraries ${required_libraries} hwy )
    endif()
  endif()
else()
  # if we're not using goading, look for SIMD libraries which zimt
  # can use. prefer hwy over Vc and Vc over std::simd.
  find_package ( HWY CONFIG )
  if ( HWY_FOUND )
    set ( SIMD_LIBRARY "-D USE_HWY" )
    message ( STATUS "find_package: HWY_FOUND = true" )
    if ( MULTI_SIMD_ISA )
      set ( required_libraries ${required_libraries} hwy )
    endif()
  else()
    # with Vc and std::simd we can't use highway's foreach_target
    set ( MULTI_SIMD_ISA OFF )
    find_package ( Vc )
    if ( Vc_FOUND )
      message ( STATUS "find_package: Vc_FOUND = true" )
      set ( required_libraries ${required_libraries} Vc )
      set ( SIMD_LIBRARY "-D USE_VC" )
    else()
      if ( NOT APPLE )
        # the std::simd build doesn't work on my mac
        set ( SIMD_LIBRARY "-D USE_STDSIMD" )
      endif()
    endif()
  endif()
endif()

# For now, the only windows build is using msys2/mingw. We look at
# a few CMake variables and set "WINDOWS" accordingly.

if(WIN32)
  message(STATUS "***** found CMake::WIN32, setting WINDOWS TRUE")
  set(WINDOWS TRUE)
elseif(MINGW)
  message(STATUS "***** found CMake::MINGW, setting WINDOWS TRUE")
  set(WINDOWS TRUE)
elseif(CYGWIN)
  message(STATUS "***** found CMake::CYGWIN, setting WINDOWS TRUE")
  set(WINDOWS TRUE)
else()
  set(WINDOWS FALSE)
endif()

# on the mac, macPorts installs to /opt/local, hence:

if ( APPLE )
  link_directories ( /usr/local/lib /opt/local/lib )
  include_directories ( /usr/local/include /opt/local/include )
elseif ( WINDOWS )
  include_directories ( /usr/local/include /mingw64/include )
  link_directories ( /usr/local/lib /mingw64/lib )
endif()

set( required_libraries
      ${required_libraries} OpenImageIO OpenImageIO_Util )

add_executable(${_project_name} envutil.cc basic.cc)

target_link_libraries(${_project_name} pthread ${required_libraries} avutil avcodec swscale hwy)

# TODO: currently, optimization with -Ofast does not work

# note: compiler options for AVX2, now used by default:

#      "-O3 -mavx2 -march=haswell -mpclmul -maes ${SIMD_LIBRARY}")

if ( MULTI_SIMD_ISA )
  set_source_files_properties ( envutil.cc PROPERTIES COMPILE_FLAGS
        "-O3 -Wno-enum-compare -g -DMULTI_SIMD_ISA -I.. ${SIMD_LIBRARY}")
else()
set_source_files_properties ( envutil.cc PROPERTIES COMPILE_FLAGS
      "-O3 -Wno-enum-compare -mavx2 -march=haswell -mpclmul -maes ${SIMD_LIBRARY}")
endif()

# install the 'useful' binaries (stepper is only a demo)

install(TARGETS envutil DESTINATION bin)

# Some more options for building debian packages. The first one - with

set(CPACK_DEBIAN_PACKAGE_NAME "envutil")
set(CPACK_DEBIAN_PACKAGE_VERSION "0.1.1")

set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CMAKE_INSTALL_PREFIX})
set(CPACK_SET_DESTDIR "ON")

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Kay F. Jahnke")

# insert license

file(READ LICENSE CPACK_PACKAGE_LICENSE)

# insert description

set(CPACK_PACKAGE_DESCRIPTION_FILE
  "${CMAKE_SOURCE_DIR}/brief_description.txt")

set(CPACK_DEBIAN_PACKAGE_DESCRIPTION
  "envutil: utility program for 360X180 degree images")

set(CPACK_DEBIAN_PACKAGE_SECTION "universe/graphics")
set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://github.com/kfjahnke/envutil")

# my debian system did not find the OIIO dependency automatically, hence

set(CPACK_DEBIAN_PACKAGE_DEPENDS "libopenimageio2.4")

# This option generates dependencies automatically - it requires
# the dpkg-dev package to function

set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

# is 'VENDOR' the right term here? This is FOSS

set(CPACK_PACKAGE_VENDOR "Kay F. Jahnke")

set(CPACK_STRIP_FILES 1)
set(CPACK_PACKAGE_VERSION_MAJOR ${envutil_major})
set(CPACK_PACKAGE_VERSION_MINOR ${envutil_minor})
set(CPACK_PACKAGE_VERSION_PATCH ${envutil_patch})

set(CPACK_PACKAGE_CONTACT "Kay F. Jahnke <kfjahnke+envutil@gmail.com>")

set(CPACK_PACKAGE_ICON scripts/PV_512x512x32.png)

install(FILES scripts/256X256/lux.png
        DESTINATION share/icons/hicolor/256x256/apps)

install(FILES scripts/512X512/lux.png 
        DESTINATION share/icons/hicolor/512x512/apps)

install(FILES scripts/1024X1024/lux.png
        DESTINATION share/icons/hicolor/1024x1024/apps)

install(FILES scripts/envutil.desktop 
        DESTINATION share/applications)

# envutil is MIT-licensed

set(CPACK_DEBIAN_PACKAGE_LICENSE "MIT License")
set(CPACK_PACKAGE_LICENSE "MIT License")
set(CPACK_DEBIAN_PACKAGE_DEBUG ON)

message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
include(CPack)

# TODO: figuring out the architecture is only tested for x64_86

if (     ${CMAKE_SYSTEM_PROCESSOR} STREQUAL x86_64
      OR ${CMAKE_SYSTEM_PROCESSOR} STREQUAL AMD64 )
  set(x86_64 TRUE)
  message(STATUS "***** setting i86 TRUE for an intel/AMD target")
else()
  set(x86_64 FALSE)
  message(STATUS "***** setting i86 FALSE; not an intel/AMD target")

  # tentative: I don't know what and ARM system would produce

  if ( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL arm64 )
    message(STATUS "***** setting arm64 TRUE for an ARM target")
    set ( arm64 TRUE )
  else()
    message(STATUS "***** setting arm64 FALSE, not an ARM target")
    set ( arm64 FALSE )
    if ( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL aarch64 )
      message(STATUS "***** setting aarch64 TRUE for an SVE target")
      set ( aarch64 TRUE )
    else()
      message(STATUS "***** setting aarch64 FALSE, not an SVE target")
      set ( aarch64 FALSE )
    endif()
  endif()
endif()

# assuming we have detected the architecture, we list ISAs

if ( x86_64 )
  list ( APPEND isa_l SSE2 SSSE3 SSE4 AVX2 AVX3 AVX3_ZEN4 AVX3_SPR )
# the NEON_BF16 ISA seems not to be built by default. We've used this
# code to build envutil on a macBook pro with M1 processor:
elseif ( arm64 )
  list ( APPEND isa_l NEON_WITHOUT_AES NEON ) # NEON_BF16
# tentative, we have no SVE systems
elseif ( aarch64 )
  list ( APPEND isa_l ALL_SVE SVE SVE2 SVE_256 SVE2_128 )
else()
  # tentative. on my system, this doesn't work.
  list ( APPEND isa_l EMU128 )
endif()

# for the test program composed of several ISA-specific TUs,
# we have the main program disp_to_tu.cc, which does the
# dispatching, and basic.cc, which has ISA-independent code

add_executable(disp_to_tu disp_to_tu.cc basic.cc)

# we need to link with libhwy

target_link_libraries(disp_to_tu hwy)

# the main program needs specific compile options:

set_source_files_properties ( disp_to_tu.cc PROPERTIES COMPILE_FLAGS
                "-DUSE_HWY -O3 -std=gnu++17 -g -DMULTI_SIMD_ISA -I.." )

# for the ISA-specific object files holding 'payload' code,
# we use cmake 'object libraries'. This places the ISA-specific
# object files in separate directories, for which we use the
# same name as the ISA. For this program, each of the object
# libraries will only contain a single object file made from
# inset.cc with ISA-specific compilation instructions. Since
# we're already running a loop over the ISAs, we add a line
# to tell cmake to link the object file in.

foreach ( isa IN LISTS isa_l )

    add_library ( ${isa} OBJECT inset.cc )

    target_compile_options ( ${isa} PUBLIC -DTG_ISA=HWY_${isa}
                 -DMULTI_SIMD_ISA -DUSE_HWY -O3 -std=gnu++17 )

    target_link_libraries(disp_to_tu $<TARGET_OBJECTS:${isa}>)

endforeach()

# we set up the same sort of build for 'envutil_hwy', a highway-only
# multi-SIMD-ISA build of envutil using 'sharded' payload code.

add_executable(envutil_hwy envutil_main.cc envutil_dispatch.cc envutil_basic.cc)

# we need to link with libhwy

target_link_libraries(envutil_hwy pthread ${required_libraries} avutil avcodec swscale hwy)

# the main program needs specific compile options:

set_source_files_properties ( envutil_main.cc envutil_dispatch.cc envutil_basic.cc PROPERTIES COMPILE_FLAGS
                "-DUSE_HWY -O3 -std=gnu++17 -DMULTI_SIMD_ISA -I.." )

foreach ( isa IN LISTS isa_l )

    add_library ( eu_${isa} OBJECT envutil_payload.cc )

    target_compile_options ( eu_${isa} PUBLIC -DTG_ISA=HWY_${isa}
               -g -DMULTI_SIMD_ISA -DUSE_HWY -O3 -std=gnu++17 )

    target_link_libraries(envutil_hwy $<TARGET_OBJECTS:eu_${isa}>)

endforeach()




