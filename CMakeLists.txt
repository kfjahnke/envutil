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
set(envutil_minor 0)
set(envutil_patch 0)

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
if ( OpenImageIO_FOUND )
  message ( STATUS "find_package: Imath = true" )
endif()

# highway has added code to find it with cmake config mode with
# the release of 1.0.4. If we can't find highway, we try Vc, and
# if that doesn't work we use std::simd.

option(USE_GOADING "USE SIMD EMULATION WITH SMALL LOOPS" OFF)

if ( NOT USE_GOADING )
  find_package ( HWY CONFIG )
  if ( HWY_FOUND )
    message ( STATUS "find_package: HWY_FOUND = true" )
    set ( required_libraries ${required_libraries} hwy )
    set ( SIMD_LIBRARY "-D USE_HWY" )
  else()
    find_package ( Vc )
    if ( VC_FOUND )
      message ( STATUS "find_package: VC_FOUND = true" )
      set ( required_libraries ${required_libraries} Vc )
      set ( SIMD_LIBRARY "-D USE_VC" )
    else()
      set ( SIMD_LIBRARY "-D USE_STDSIMD" )
    endif()
  endif()
endif()

set( required_libraries
      ${required_libraries} OpenImageIO OpenImageIO_Util )

add_executable(${_project_name} envutil.cc)
# add_executable(stepper stepper.cc)
add_executable(extract extract.cc)

target_link_libraries(${_project_name} pthread ${required_libraries})
# target_link_libraries(stepper pthread ${required_libraries})
target_link_libraries(extract pthread ${required_libraries})

# TODO: currently, optimization with -Ofast does not work

# note: compiler options for AVX2:
#      "-O3 -mavx2 -march=haswell -mpclmul -maes ${SIMD_LIBRARY}")

set_source_files_properties ( envutil.cc extract.cc PROPERTIES COMPILE_FLAGS
      "-O3 -Wno-enum-compare -mavx2 -march=haswell -mpclmul -maes ${SIMD_LIBRARY}")
      # "-O3 -Wno-enum-compare ${SIMD_LIBRARY}")

# install the 'useful' binaries (stepper is only a demo)

install(TARGETS extract DESTINATION bin)
