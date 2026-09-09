# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

# 3.27 for $<COMPILE_ONLY:>, which O2PhysicsAddWorkflow.cmake uses to give a
# workflow the precompiled header carrier's compile usage requirements without
# also putting it on the link line. Builds use CMake 4.1 (alidist cmake.sh).
cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

project(O2Physics
        VERSION 0.0.1
        LANGUAGES C CXX
        DESCRIPTION "Physics Analysis for O2")

cmake_policy(SET CMP0144 NEW)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set_property(GLOBAL PROPERTY REPORT_UNDEFINED_PROPERTIES)

# How much memory this build may actually use, in MiB.
#
# cmake_host_system_information reads the HOST, which is wrong under a container
# memory limit: there is no lxcfs in the CI images, so on a 251 GiB builder
# capped at 50 GiB it reported 251 GiB and the pool below never bound. The
# result was an OOM-killed build that looked like a compiler crash:
#
#   c++: fatal error: Killed signal terminated program cc1plus
#
# So prefer the cgroup limit when there is one. cgroup v2 first, then v1; both
# report bytes, and both say "no limit" with a sentinel ("max", or a number so
# large it is really 2^63 rounded down) rather than by being absent.
cmake_host_system_information(RESULT _totalmem QUERY TOTAL_PHYSICAL_MEMORY)
foreach(_cg /sys/fs/cgroup/memory.max /sys/fs/cgroup/memory/memory.limit_in_bytes)
  if(EXISTS ${_cg})
    file(READ ${_cg} _cglimit)
    string(STRIP "${_cglimit}" _cglimit)
    if(_cglimit MATCHES "^[0-9]+$")
      math(EXPR _cgmib "${_cglimit}/1048576")
      # Anything at or above the host's own memory is "unlimited" expressed as a
      # number, not a smaller budget.
      if(_cgmib GREATER 0 AND _cgmib LESS ${_totalmem})
        set(_totalmem ${_cgmib})
        message(STATUS "Memory limit from ${_cg}: ${_totalmem} MiB")
        break()
      endif()
    endif()
  endif()
endforeach()

# 8 GiB per job, not 2. Measured on the CI builders with GCC 14.2: the heaviest
# analysis producers peak at 7.4-10.1 GiB each (derivedDataWriter 10.09,
# femtoProducer 8.27, candidateCreatorXic0Omegac0 7.36), and -fmem-report shows
# ~6.4 GB of that is frontend AST held live -- template instantiation over the
# distinct soa::Join types each task subscribes to. At 2048 the pool permitted
# four times what the memory allows, and ninja builds in declaration order, so
# the heavy producers arrive as a cluster rather than spread out.
set(ANALYSIS_COMPILE_MEMORY_MB 8192 CACHE STRING
    "Assumed peak memory of one analysis compilation, in MiB")
math(EXPR _total_analysis_jobs "(${_totalmem}-2048)/${ANALYSIS_COMPILE_MEMORY_MB}")
if(_total_analysis_jobs LESS_EQUAL 0)
  set(_total_analysis_jobs 1)
endif()
set(ANALYSIS_COMPILE_POOL ${_total_analysis_jobs} CACHE STRING "How many parallel analysis compilation jobs")
set_property(GLOBAL PROPERTY JOB_POOLS analysis=${ANALYSIS_COMPILE_POOL})
message(STATUS "Limiting workflow compilations in parallel to: ${ANALYSIS_COMPILE_POOL}")

# C++ standard
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
endif()

include(O2PhysicsDefineOptions)
o2physics_define_options()

include(O2PhysicsDefineOutputPaths)
o2physics_define_output_paths()

include(O2PhysicsDefineRPATH)
o2physics_define_rpath()

# External dependencies
include(dependencies/CMakeLists.txt)

# Include macros
include(O2PhysicsNameTarget)
include(O2PhysicsAddExecutable)
include(O2PhysicsAddWorkflow)
include(O2PhysicsAddLibrary)
include(O2PhysicsAddHeaderOnlyLibrary)
include(O2PhysicsTargetRootDictionary)
include(O2PhysicsSetROOTPCMDependencies)
include(O2PhysicsDataFile)
include(AddRootDictionary)

# Main targets of the project in various subdirectories. Order matters.
add_subdirectory(Common)
add_subdirectory(ALICE3)
add_subdirectory(DPG)

add_subdirectory(PWGCF)
add_subdirectory(PWGDQ)
add_subdirectory(PWGEM)
add_subdirectory(PWGHF)
add_subdirectory(PWGJE)
add_subdirectory(PWGLF)
add_subdirectory(PWGMM)
add_subdirectory(PWGUD)

add_subdirectory(Tools)
add_subdirectory(Tutorials)
add_subdirectory(EventFiltering)

add_subdirectory(Scripts)

# Testing and packaging only needed if we are the top level directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  add_subdirectory(packaging)
endif()

set_root_pcm_dependencies()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/ClangBuildAnalyzer.ini.in
               ${CMAKE_CURRENT_BINARY_DIR}/ClangBuildAnalyzer.ini @ONLY)
