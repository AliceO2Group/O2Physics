# Copyright 2019-2025 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

include_guard()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# Enabled warnings supported by Clang and GCC, not treated as errors
set(O2PHYSICS_WARNINGS_COMMON_NO_ERROR "")

# Enabled warnings supported by Clang only, not treated as errors
set(O2PHYSICS_WARNINGS_CLANG_NO_ERROR "")

# Enabled warnings supported by GCC only, not treated as errors
set(O2PHYSICS_WARNINGS_GCC_NO_ERROR "")

# Function to build a list of warning flags from their names
function(o2_build_warning_flags)
  cmake_parse_arguments(PARSE_ARGV 0 A "" "PREFIX;OUTPUTVARNAME" "WARNINGS")
  if(A_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Unexpected unparsed arguments: ${A_UNPARSED_ARGUMENTS}")
  endif()
  list(TRANSFORM A_WARNINGS STRIP)
  list(TRANSFORM A_WARNINGS PREPEND ${A_PREFIX})
  string(JOIN " " OUTPUT ${A_WARNINGS})
  set(${A_OUTPUTVARNAME} ${OUTPUT} PARENT_SCOPE)
endfunction()

message(STATUS "O2PHYSICS_WARNINGS_AS_ERRORS: ${O2PHYSICS_WARNINGS_AS_ERRORS}")

# Treat warnings as errors.
if(O2PHYSICS_WARNINGS_AS_ERRORS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  # Set warning flags for all platforms.
  o2_build_warning_flags(PREFIX "-Wno-error="
    OUTPUTVARNAME O2PHYSICS_CXX_WARNINGS_COMMON_NO_ERROR
    WARNINGS ${O2PHYSICS_WARNINGS_COMMON_NO_ERROR})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${O2PHYSICS_CXX_WARNINGS_COMMON_NO_ERROR}")
  if(APPLE)
    # Set warning flags for macOS only.
    o2_build_warning_flags(PREFIX "-Wno-error="
      OUTPUTVARNAME O2PHYSICS_CXX_WARNINGS_APPLE_NO_ERROR
      WARNINGS ${O2PHYSICS_WARNINGS_CLANG_NO_ERROR})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${O2PHYSICS_CXX_WARNINGS_APPLE_NO_ERROR}")
  elseif(UNIX)
    # Set warning flags for Linux only.
    o2_build_warning_flags(PREFIX "-Wno-error="
      OUTPUTVARNAME O2PHYSICS_CXX_WARNINGS_UNIX_NO_ERROR
      WARNINGS ${O2PHYSICS_WARNINGS_GCC_NO_ERROR})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${O2PHYSICS_CXX_WARNINGS_UNIX_NO_ERROR}")
  endif()
endif()

if(ENABLE_TIMETRACE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftime-trace")
endif()

if(APPLE)
else()
set(SPLIT_DWARF "-gsplit-dwarf=single")
endif()

set(CMAKE_CXX_FLAGS_COVERAGE "${SPLIT_DWARF} -O2 -fprofile-arcs -ftest-coverage")
set(CMAKE_C_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_COVERAGE}")
set(CMAKE_Fortran_FLAGS_COVERAGE "${SPLIT_DWARF} -O2 -fprofile-arcs -ftest-coverage")
set(CMAKE_LINK_FLAGS_COVERAGE "--coverage -fprofile-arcs  -fPIC")

mark_as_advanced(
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_Fortran_FLAGS_COVERAGE
    CMAKE_LINK_FLAGS_COVERAGE)

#Check the compiler and set the compile and link flags
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Set BuildType to DEBUG")
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()

if(ENABLE_CASSERT) #For the CI, we want to have <cassert> assertions enabled
  set(CMAKE_CXX_FLAGS_RELEASE "-O2 -pipe")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 ${SPLIT_DWARF} -pipe")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG -pipe")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 ${SPLIT_DWARF} -DNDEBUG -pipe")
  if(CMAKE_BUILD_TYPE STREQUAL "RELEASE"
      OR CMAKE_BUILD_TYPE STREQUAL "RELWITHDEBINFO")
    set(FAIR_MIN_SEVERITY "info")
  endif()
endif()
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 ${SPLIT_DWARF}")
# make sure Debug build not optimized (does not seem to work without CACHE + FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "${SPLIT_DWARF} -O0" CACHE STRING "Debug mode build flags" FORCE)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "Debug mode build flags" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG "${SPLIT_DWARF} -O0" CACHE STRING "Debug mode build flags" FORCE)

if(APPLE)
elseif(UNIX)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined") # avoid undefined in our libs
endif()

message(STATUS "Using build type: ${CMAKE_BUILD_TYPE} - CXXFLAGS: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
