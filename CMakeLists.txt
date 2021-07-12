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
        DESCRIPTION "Analysis framework for O2")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set_property(GLOBAL PROPERTY REPORT_UNDEFINED_PROPERTIES)

include(GNUInstallDirs)

# Set the default build type to "RelWithDebInfo"
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo"
      CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel Coverage." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# C++ standard
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
endif()

# Dependencies
find_package(O2 CONFIG REQUIRED)

feature_summary(WHAT ALL FATAL_ON_MISSING_REQUIRED_PACKAGES)

# Include macros
include(O2PhysicsNameTarget)
include(O2PhysicsAddExecutable)
include(O2PhysicsAddWorkflow)

# Main targets of the project in various subdirectories. Order matters.
add_subdirectory(Examples)
