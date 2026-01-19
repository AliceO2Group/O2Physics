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

cmake_minimum_required(VERSION 3.23 FATAL_ERROR)

project(O2Physics
        VERSION 0.0.1
        LANGUAGES C CXX
        DESCRIPTION "Physics Analysis for O2")

cmake_policy(SET CMP0144 NEW)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set_property(GLOBAL PROPERTY REPORT_UNDEFINED_PROPERTIES)

cmake_host_system_information(RESULT _totalmem QUERY TOTAL_PHYSICAL_MEMORY)
math(EXPR _total_analysis_jobs "(${_totalmem}-2048)/2048")
if(_total_analysis_jobs LESS_EQUAL 0)
  set(_total_analysis_jobs 1)
endif()
set(ANALYSIS_COMPILE_POOL ${_total_analysis_jobs} CACHE STRING "How many parallel analysis compilation jobs")
set(HIGH_MEMORY_COMPILE_POOL 1 CACHE STRING "Pool for tasks which use a lot of memory")
set_property(GLOBAL PROPERTY JOB_POOLS analysis=${ANALYSIS_COMPILE_POOL} serial=${HIGH_MEMORY_COMPILE_POOL})
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
