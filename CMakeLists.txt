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

# envutil should be built with clang++

set(ENV{CXX} clang++)
set(ENV{CC} clang)

# get pv-version string from envutil_version_first envutil_version_second
# and envutil_version_third

macro(getversion file out)
	set(_vers_lines "")
	file(STRINGS ${file} _vers_lines)
	set(number 0)
	foreach(_l ${_vers_lines} )
		if(_l MATCHES "^\([0-9]+\)$")
			# the last one takes over
			set(number ${CMAKE_MATCH_1})
		endif()
	endforeach()

	set(${out} ${number})
endmacro()

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

# We switch of some warnings cluttering the build, assuming it's safe
# to ignore them.

set(RENDERING_FLAGS "-fno-math-errno -Wno-unused-value")

# handle use of SIMD libraries. for now, we're not using any of the
# SIMD back-ends, but rely on small loops and autovectorization.

option(USE_GOADING "USE SIMD EMULATION WITH SMALL LOOPS" ON)
option(USE_HWY_LIBRARY "ACTIVATES USE OF HIGHWAY SIMD LIBRARY" OFF)
option(USE_VC_LIBRARY "ACTIVATES USE OF VC SIMD LIBRARY" OFF)
option(USE_STDSIMD "ACTIVATES USE OF STD::SIMD" OFF)

# OIIO anticipates moving to C++17 in 2024, and I may want to access
# std::filesystem coming with C++17. So I raise the build to using C++17
# as overall standard.

set(CMAKE_CXX_STANDARD 17)
message(STATUS "***** using C++17 standard")

set_source_files_properties ( envutil.cc PROPERTIES COMPILE_FLAGS "-O3")

add_executable(${_project_name} envutil.cc)

if ( NOT APPLE )
  find_package ( OpenImageIO CONFIG REQUIRED )
  if ( OpenImageIO_FOUND )
    message ( STATUS "find_package: OpenImageIO = true" )
  endif()

  if ( NOT "${OpenImageIO_INCLUDE_DIR}" STREQUAL "" )
    message ( STATUS "OpenImageIO_INCLUDE_DIR = ${OpenImageIO_INCLUDE_DIR}" )
    include_directories ( ${OpenImageIO_INCLUDE_DIR} )
  endif()
endif()

set( required_libraries
      ${required_libraries} OpenImageIO OpenImageIO_Util )

target_link_libraries(${_project_name} pthread ${required_libraries})

# install the binary

install(TARGETS ${_project_name} DESTINATION bin)
