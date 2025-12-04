#             cmake file for envutil
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

# You need some dependencies: OpenImageIO and Imath. highway is recommended.

# To build, try this for starters (ninja instead of make should also work):

# mkdir build
# cd build
# cmake ..
# make

# I usually use clang, and to make a fully optimized build, I use:

# cmake -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release ..

# I recently recreated the cmake file with the help of chatGPT.

cmake_minimum_required(VERSION 3.16)
project(envutil VERSION 0.3.0 LANGUAGES CXX)

# -----------------------
# Basic project settings
# -----------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(MULTI_SIMD_ISA "Build per-ISA object files and use custom dispatching" ON)
option(ENABLE_STRICT_WARNINGS "Enable -Wall -Wextra on supported compilers" OFF)

message(STATUS "Project: ${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "MULTI_SIMD_ISA = ${MULTI_SIMD_ISA}")

# -----------------------
# Detect platform (simple)
# TODO: handle e.g. RiscV properly
# -----------------------

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)$")
  set(x86_64 TRUE)
else()
  set(x86_64 FALSE)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)$")
  set(arm64 TRUE)
else()
  set(arm64 FALSE)
endif()

# --------
# ISA list
# TODO: handle e.g. RiscV properly
# -----------------------
if (x86_64)
  set(isa_l SSE2 SSSE3 SSE4 AVX2 AVX3 AVX3_ZEN4 AVX3_SPR AVX3_DL AVX10_2)
elseif (arm64)
  set(isa_l NEON_WITHOUT_AES NEON NEON_BF16)
else()
  set(isa_l EMU128)
endif()

# -----------------------
# Find dependencies
# -----------------------

find_package(OpenImageIO CONFIG REQUIRED)
find_package(Imath CONFIG REQUIRED)
find_package(Threads REQUIRED)

# Highway discovery (optional).

find_package(HWY CONFIG)
if(HWY_FOUND)
  message(STATUS "Highway found: HWY_FOUND=TRUE")
else()
  message(STATUS "Highway not found")
endif()

# -----------------------
# Common compile flags (target-local)
# -----------------------

# Keep aggressive optimization in Release configuration only

set(COMMON_RELEASE_FLAGS "-O3")

# SIMD backend selection: user may set BACKEND via -DBACKEND=highway or std.simd

if(NOT DEFINED BACKEND)
  if(HWY_FOUND)
    set(BACKEND "highway")
  else()
    set(BACKEND "std::simd")
  endif()
endif()
message(STATUS "SIMD BACKEND = ${BACKEND}")

if("${BACKEND}" STREQUAL "highway")
  set(SIMD_BACKEND_FLAG "-DUSE_HWY")
elseif("${BACKEND}" STREQUAL "std::simd")
  set(SIMD_BACKEND_FLAG "-DUSE_STDSIMD")
elseif("${BACKEND}" STREQUAL "Vc")
  set(SIMD_BACKEND_FLAG "-DUSE_VC")
  find_package(Vc QUIET)
  if(TARGET Vc::Vc)
    list(APPEND EXTRA_LINK_LIBS Vc::Vc)
  else()
    message(FATAL_ERROR "Vc backend requested but Vc not found")
  endif()
else()
  set(SIMD_BACKEND_FLAG "-DUSE_GOADING")
endif()

# -----------------------
# Source files
# -----------------------

set(SRC_MAIN
  envutil_main.cc
  envutil_dispatch.cc
  envutil_basic.cc
)

set(SRC_PAYLOAD envutil_payload.cc)

# -----------------------
# Build path: multi-ISA vs single
# -----------------------

if(MULTI_SIMD_ISA)
  if(NOT HWY_FOUND)
    message(WARNING "MULTI_SIMD_ISA requested but HWY not found. Continuing but you may need highway for TG_ISA macros to make sense.")
  endif()

  # Create object libraries, one per ISA

  set(OBJECT_LIBS "")
  foreach(isa IN LISTS isa_l)
    string(TOUPPER "${isa}" isa_u)
    set(obj_name "eu_${isa}")

    add_library(${obj_name} OBJECT ${SRC_PAYLOAD})

    # Per-object compile options: define TG_ISA macro and enable multi-simd
    # define. We choose to set the macro via compile options here so token
    # pasting works

    # TODO code is a bit dense

    target_compile_options(${obj_name} PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:${COMMON_RELEASE_FLAGS} -DTG_ISA=HWY_${isa_u} -DMULTI_SIMD_ISA ${SIMD_BACKEND_FLAG}>
      $<$<CXX_COMPILER_ID:MSVC>:${COMMON_RELEASE_FLAGS} -DTG_ISA=HWY_${isa_u} -DMULTI_SIMD_ISA ${SIMD_BACKEND_FLAG}>
    )

    # Optionally add include path pointing to project root for payload includes

    target_include_directories(${obj_name} PRIVATE ${PROJECT_SOURCE_DIR})

    list(APPEND OBJECT_LIBS $<TARGET_OBJECTS:${obj_name}>)
  endforeach()

  # Main executable links to object files

  add_executable(envutil ${SRC_MAIN})
  target_include_directories(envutil PRIVATE ${PROJECT_SOURCE_DIR})

  if(OBJECT_LIBS)
    target_link_libraries(envutil PRIVATE ${OBJECT_LIBS})
  endif()

else()
  # Non multi-ISA: compile payload into the main target

  add_executable(envutil ${SRC_MAIN} ${SRC_PAYLOAD})
  target_include_directories(envutil PRIVATE ${PROJECT_SOURCE_DIR})

  # compile flags for payload are inherited from target settings below

endif()

# -----------------------
# Link libraries
# -----------------------
# Prefer imported targets if available

if(TARGET OpenImageIO::OpenImageIO)
  target_link_libraries(envutil PRIVATE OpenImageIO::OpenImageIO OpenImageIO::OpenImageIO_Util)
else()
  target_link_libraries(envutil PRIVATE OpenImageIO OpenImageIO_Util)
endif()

if(TARGET Imath::Imath)
  target_link_libraries(envutil PRIVATE Imath::Imath)
else()
  target_link_libraries(envutil PRIVATE Imath)
endif()

if(EXTRA_LINK_LIBS)
  target_link_libraries(envutil PRIVATE ${EXTRA_LINK_LIBS})
endif()

target_link_libraries(envutil PRIVATE Threads::Threads)

target_link_libraries(envutil PRIVATE hwy)

# -----------------------
# Compiler options for the main target
# -----------------------
# Release-only optimization

target_compile_options(envutil PRIVATE $<$<CONFIG:Release>:${COMMON_RELEASE_FLAGS}>)

# Expose SIMD backend define to the top-level too (useful for headers)

target_compile_definitions(envutil PRIVATE $<$<BOOL:${MULTI_SIMD_ISA}>:MULTI_SIMD_ISA>)
if("${BACKEND}" STREQUAL "highway")
  target_compile_definitions(envutil PRIVATE USE_HWY)
elseif("${BACKEND}" STREQUAL "std::simd")
  target_compile_definitions(envutil PRIVATE USE_STDSIMD)
endif()

# Strict warnings if requested

if(ENABLE_STRICT_WARNINGS)
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(envutil PRIVATE -Wall -Wextra -Wpedantic)
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(envutil PRIVATE /W4)
  endif()
endif()

# chatGPT proposed adding IPO optimization, but this maxed out my system
# and I stopped it after some time, so I took out this section:

# -----------------------
# IPO / LTO if available
# -----------------------
# include(CheckIPOSupported)
# check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
# if(ipo_supported)
  # set_property(TARGET envutil PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  # # also enable for object-libs if needed
  # foreach(isa IN LISTS isa_l)
    # if(TARGET eu_${isa})
      # set_property(TARGET eu_${isa} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    # endif()
  # endforeach()
# endif()

# -----------------------
# Install rules & packaging
# -----------------------

install(TARGETS envutil RUNTIME DESTINATION bin)

install(FILES copyright LICENSE DESTINATION share/doc/envutil)
install(FILES scripts/envutil.desktop DESTINATION share/applications)
install(FILES scripts/256X256/lux.png DESTINATION share/icons/hicolor/256x256/apps)
install(FILES scripts/512X512/lux.png DESTINATION share/icons/hicolor/512x512/apps)
install(FILES scripts/1024X1024/lux.png DESTINATION share/icons/hicolor/1024x1024/apps)

# CPack minimal config

include(CPack)
set(CPACK_PACKAGE_NAME "envutil")
set(CPACK_PACKAGE_VENDOR "Kay F. Jahnke")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})

# -----------------------
# Diagnostics
# -----------------------
message(STATUS "Build summary:")
message(STATUS "  MULTI_SIMD_ISA = ${MULTI_SIMD_ISA}")
message(STATUS "  BACKEND = ${BACKEND}")
message(STATUS "  HWY_FOUND = ${HWY_FOUND}")
message(STATUS "  ISA list = ${isa_l}")

