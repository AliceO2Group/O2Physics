<<<<<<< HEAD
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

cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

project(O2Physics
        VERSION 0.0.1
        LANGUAGES CXX
        DESCRIPTION "Physics Analysis for O2")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set_property(GLOBAL PROPERTY REPORT_UNDEFINED_PROPERTIES)

include(GNUInstallDirs)

cmake_host_system_information(RESULT _totalmem QUERY TOTAL_PHYSICAL_MEMORY)
math(EXPR _total_analysis_jobs "(${_totalmem}-4096)/10240")
if(_total_analysis_jobs LESS_EQUAL 0)
  set(_total_analysis_jobs 1)
endif()
set(ANALYSIS_COMPILE_POOL ${_total_analysis_jobs} CACHE STRING "How many parallel analysis compilation jobs")
set_property(GLOBAL PROPERTY JOB_POOLS analysis=${ANALYSIS_COMPILE_POOL})

# C++ standard
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
endif()

# Dependencies
find_package(O2 CONFIG REQUIRED)

feature_summary(WHAT ALL FATAL_ON_MISSING_REQUIRED_PACKAGES)

include(O2PhysicsDefineOptions)
o2physics_define_options()

# External dependencies
include(dependencies/CMakeLists.txt)

# Include macros
include(O2PhysicsNameTarget)
include(O2PhysicsAddExecutable)
include(O2PhysicsAddWorkflow)

# Main targets of the project in various subdirectories. Order matters.
add_subdirectory(Examples)
=======
# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

add_subdirectory(Core)
add_subdirectory(DataModel)
add_subdirectory(Tasks)
add_subdirectory(Tutorials)
>>>>>>> Analysis: update directory structure (#2868)
