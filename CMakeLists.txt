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
set(envutil_minor 2)
set(envutil_patch 0)

message(STATUS "VERSION = ${envutil_major}.${envutil_minor}.${envutil_patch}")

# Project name 'envutil'

set(project_name envutil)
project (envutil)

# specify some directories for installation of components. Note that
# some builds fail if CMAKE_INSTALL_PREFIX points to folders with
# white space or non-alphanumeric characters, like 'Program Files (x86)'
# on windows. If that happens, use a 'sane' path like '/usr/local'.

message(STATUS "***** CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}")

# work out the target processor architecture by checking for all values
# of CMAKE_SYSTEM_PROCESSOR which, on any platform, indicates an intel
# processor.

message(STATUS "***** SYSTEM PROCESSOR ${CMAKE_SYSTEM_PROCESSOR}")

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
  list ( APPEND isa_l SSE2 SSSE3 SSE4 AVX2 AVX3 AVX3_ZEN4 AVX3_SPR AVX3_DL AVX10_2)
# the NEON_BF16 ISA seems not to be built by default. We've used this
# code to build envutil on a macBook pro with M1 processor:
elseif ( arm64 )
  list ( APPEND isa_l NEON_WITHOUT_AES NEON NEON_BF16 )
# tentative, we have no SVE systems
elseif ( aarch64 )
  list ( APPEND isa_l ALL_SVE SVE SVE2 SVE_256 SVE2_128 )
else()
  # tentative. on my system, this doesn't work.
  list ( APPEND isa_l EMU128 )
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

option(MULTI_SIMD_ISA "AUTOMATIC SIMD ISA DISPATCH (REQUIRES HIGHWAY)" ON)

# if the user hasn't #defined 'BACKEND', we try for highway. if that's
# available, we use it as back-end, and std::simd if it's not.

if ( NOT DEFINED BACKEND )
  message ( STATUS "BACKEND variable is not set, looking for highway" )
  find_package ( HWY CONFIG )
  if ( HWY_FOUND )
    SET ( BACKEND "highway" )
    message ( STATUS "found highway, using highway back-end" )
  else()
    SET ( BACKEND "std::simd" )
    message ( STATUS "highway not found, using std::simd back-end" )
  endif()
endif()

# now we set 'SIMD_BACKEND' according to the 'BACKEND' variable

if ( "${BACKEND}" STREQUAL "highway" )
  set ( SIMD_BACKEND "-DUSE_HWY" )
elseif ( "${BACKEND}" STREQUAL "std::simd" )
  set ( SIMD_BACKEND "-DUSE_STDSIMD" )
elseif ( "${BACKEND}" STREQUAL "Vc" )
  # if the user wants to use the Vc back-end, we must have libVc, and
  # we add Vc to the required libraries to link with.
  set ( SIMD_BACKEND "-DUSE_VC" )
  find_package ( Vc REQUIRED )
  if ( Vc_FOUND )
    message ( STATUS "find_package: Vc_FOUND = true" )
    set ( required_libraries ${required_libraries} Vc )
  endif()
else()
  # -DUSE_GOADING doesn't actually have an effect, but it's clearer
  # than simply clearing 'SIMD_BACKEND'
  set ( SIMD_BACKEND "-DUSE_GOADING" )
endif()

if ( MULTI_SIMD_ISA OR ( "${BACKEND}" STREQUAL "highway" ) )

  message ( STATUS "build is using highway; backend: ${BACKEND}" )

  # for multi-ISA builds, and for use of the highway SIMD back-end,
  # we positively need highway. If it wasn't already found earlier
  # we now do the find_package with REQUIRED flag.

  if ( NOT HWY_FOUND )
    find_package ( HWY CONFIG REQUIRED )
    if ( HWY_FOUND )
      message ( STATUS "find_package: HWY_FOUND = true" )
    endif()
  endif()

  # if we arrive here, the find_package call must have succeeded.

  set ( required_libraries ${required_libraries} hwy )

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

# this section was provided by Karl Krieger:

# on the mac, macPorts installs to /opt/local, hence:
# 2025-01-06 Comment Karl: for homebrew we have to use /opt/homebrew
# Code added to handle link and include directories on Apple architecture
# for either MacPorts, Homebrew or none.

if (APPLE)
  if (EXISTS "/opt/local")
    message(STATUS "Detected MacPorts installation. Adding /opt/local/lib and /opt/local/include.")
    link_directories(/usr/local/lib /opt/local/lib)
    include_directories(/usr/local/include /opt/local/include)
  elseif (EXISTS "/opt/homebrew")
    message(STATUS "Detected Homebrew installation. Adding /opt/homebrew/lib and /opt/homebrew/include.")
    link_directories(/usr/local/lib /opt/homebrew/lib)
    include_directories(/usr/local/include /opt/homebrew/include)
  else()
    message(WARNING "Neither /opt/local (MacPorts) nor /opt/homebrew (Homebrew) directories exist. Defaulting to /usr/local.")
    link_directories(/usr/local/lib)
    include_directories(/usr/local/include)
  endif()
elseif (WINDOWS)
  include_directories(/usr/local/include /mingw64/include)
  link_directories(/usr/local/lib /mingw64/lib)
endif()

# replacing my original code here:

# # on the mac, macPorts installs to /opt/local, hence:
# 
# if ( APPLE )
#   link_directories ( /usr/local/lib /opt/local/lib )
#   include_directories ( /usr/local/include /opt/local/include )
# elseif ( WINDOWS )
#   include_directories ( /usr/local/include /mingw64/include )
#   link_directories ( /usr/local/lib /mingw64/lib )
# endif()

set( required_libraries
      ${required_libraries} OpenImageIO OpenImageIO_Util )

set ( COMMON_ARGS "-O3" )

if ( MULTI_SIMD_ISA )

  set_source_files_properties (
    envutil_main.cc envutil_dispatch.cc envutil_basic.cc
    PROPERTIES COMPILE_FLAGS "${COMMON_ARGS} -DMULTI_SIMD_ISA -I.." )

  add_executable ( envutil
                   envutil_main.cc envutil_dispatch.cc envutil_basic.cc)

  foreach ( isa IN LISTS isa_l )

      add_library ( eu_${isa} OBJECT envutil_payload.cc )

      # Here, highway sets the compiler flags for each ISA internally,
      # the TG_ISA variable 'tells' the build which ISA to use for the
      # object file at hand.

      # TODO: would use '${COMMON_ARGS} here, but build fails,
      # hence stating explicitly: -O3

      target_compile_options (
        eu_${isa} PUBLIC -DTG_ISA=HWY_${isa} -DMULTI_SIMD_ISA -O3 ${SIMD_BACKEND} )

      target_link_libraries ( envutil $<TARGET_OBJECTS:eu_${isa}>)

  endforeach()

else()

  set_source_files_properties (
    envutil_main.cc envutil_dispatch.cc envutil_basic.cc
    PROPERTIES COMPILE_FLAGS ${COMMON_ARGS} )

  # if you're building envutil for a specific target ISA, you can add
  # ISA-specific compiler flags via 'ISA_SPECIFIC_ARGS'

  set_source_files_properties ( envutil_payload.cc
                  PROPERTIES COMPILE_FLAGS
                  "${COMMON_ARGS} ${SIMD_BACKEND} ${ISA_SPECIFIC_ARGS}" )

  add_executable ( envutil
                   envutil_main.cc envutil_basic.cc envutil_payload.cc)

endif()

target_link_libraries( envutil ${required_libraries} pthread)

# install the 'useful' binaries (stepper is only a demo)

install(TARGETS envutil DESTINATION bin)
install(FILES copyright LICENSE DESTINATION share/doc/envutil)

# Some more options for building debian packages. The first one - with

set(CPACK_DEBIAN_PACKAGE_NAME "envutil")
set(CPACK_DEBIAN_PACKAGE_VERSION "0.2.0")

set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CMAKE_INSTALL_PREFIX})
set(CPACK_SET_DESTDIR "ON")

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Kay F. Jahnke")

# insert license

file(READ LICENSE CPACK_PACKAGE_LICENSE)

# insert description

set(CPACK_PACKAGE_DESCRIPTION_FILE
  "${CMAKE_SOURCE_DIR}/brief_description.txt")

set(CPACK_DEBIAN_PACKAGE_DESCRIPTION
  "envutil: utility to process oriented images and environments")

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

# envutil is MIT-licensed. The license info doesn't seem to make it
# into the packet, though.

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")

set(CPACK_DEBIAN_PACKAGE_LICENSE mit)
set(CPACK_PACKAGE_LICENSE mit)

set(CPACK_DEBIAN_PACKAGE_DEBUG ON)

install(FILES scripts/256X256/lux.png
        DESTINATION share/icons/hicolor/256x256/apps)

install(FILES scripts/512X512/lux.png 
        DESTINATION share/icons/hicolor/512x512/apps)

install(FILES scripts/1024X1024/lux.png
        DESTINATION share/icons/hicolor/1024x1024/apps)

install(FILES scripts/envutil.desktop 
        DESTINATION share/applications)

message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
include(CPack)
