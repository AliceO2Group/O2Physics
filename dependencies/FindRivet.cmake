# Copyright 2024-2099 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.
#
# - Try to find Rivet
#
# This searches for the configuration script "rivet-config".  If this
# script is not found, then the package is not considered found.
#
# The targets
#
#  HepMC3
#  FastJet::Contrib
#  Yoda
#
# _must_ previously have been found.
#
# Defines the target "Rivet" as well as the variables
#
#  Rivet_FOUND
#  RIVET_VERSION
#  RIVET_PREFIX
#  RIVET_BIN_DIR
#  RIVET_CPPFLAGS
#  RIVET_CXXFLAGS
#  RIVET_LDFLAGS
#  RIVET_INCLUDE_DIR
#  RIVET_LIBRARY
#  RIVET_LIB_DIR


if (NOT HepMC3_FOUND)
  # HepMC3 was not previously found, give up
  message(HepMC3 not found, disabling Rivet)
  set(Rivet_FOUND FALSE)
  return()
endif()

if (NOT fjcontrib_FOUND)
  # HepMC3 was not previously found, give up
  message(FastJet::Contrib not found, disabling Rivet)
  set(Rivet_FOUND FALSE)
  return()
endif()

if (NOT Yoda_FOUND)
  # Yoda was not previously found, give up
  message(Yoda not found, disabling Rivet)
  set(Rivet_FOUND FALSE)
  return()
endif()

if (TARGET Rivet)
   # nothing to do, target is already there
   set(Rivet_FOUND TRUE)
   return()
endif()

# the following disables all default paths (either from cmake, from environment)
FIND_PATH (RIVET_BIN_DIR rivet-config
  HINTS ${RIVET_ROOT} ENV RIVET_ROOT
)

if (RIVET_BIN_DIR)
  # MESSAGE("Executing ${RIVET_BIN_DIR}/rivet-config")
  EXECUTE_PROCESS (COMMAND ${RIVET_BIN_DIR}/rivet-config --version
    OUTPUT_VARIABLE RIVET_VERSION
  )
  # Remove white-space
  STRING (REGEX REPLACE "[ \t\r\n]+" "" RIVET_VERSION "${RIVET_VERSION}")
  # Remove trailing patch
  STRING (REGEX REPLACE "-.*" "" RIVET_VERSION "${RIVET_VERSION}")
  # message(STATUS "Found Rivet version: ${RIVET_VERSION}")

  EXECUTE_PROCESS (COMMAND ${RIVET_BIN_DIR}/rivet-config --prefix
    OUTPUT_VARIABLE RIVET_PREFIX
  )
  # message(STATUS "Found Rivet prefix: ${RIVET_PREFIX}")

  EXECUTE_PROCESS (COMMAND ${RIVET_BIN_DIR}/rivet-config --includedir
    OUTPUT_VARIABLE RIVET_INCLUDE_DIR)
  string (REGEX REPLACE "\n$" "" RIVET_INCLUDE_DIR "${RIVET_INCLUDE_DIR}")
  # message(STATUS "RIVET_INCLUDE_DIR=${RIVET_INCLUDE_DIR}")

  EXECUTE_PROCESS (COMMAND ${RIVET_BIN_DIR}/rivet-config --libdir
    OUTPUT_VARIABLE RIVET_LIB_DIR)
  string (REGEX REPLACE "\n$" "" RIVET_LIB_DIR "${RIVET_LIB_DIR}")
  # message(STATUS "RIVET_LIB_DIR=${RIVET_LIB_DIR}")

  EXECUTE_PROCESS (COMMAND ${RIVET_BIN_DIR}/rivet-config --cppflags
    OUTPUT_VARIABLE RIVET_CPPFLAGS)
  string (REGEX REPLACE "\n$" "" RIVET_CPPFLAGS "${RIVET_CPPFLAGS}")
  # message(STATUS "RIVET_CPPFLAGS=${RIVET_CPPFLAGS}")

  EXECUTE_PROCESS (COMMAND ${RIVET_BIN_DIR}/rivet-config --cxxflags
    OUTPUT_VARIABLE RIVET_CXXFLAGS)
  string (REGEX REPLACE "\n$" "" RIVET_CXXFLAGS "${RIVET_CXXFLAGS}")
  # message(STATUS "RIVET_CXXFLAGS=${RIVET_CXXFLAGS}")

  EXECUTE_PROCESS (COMMAND ${RIVET_BIN_DIR}/rivet-config --ldflags
    OUTPUT_VARIABLE RIVET_LDFLAGS)
  string (REGEX REPLACE "\n$" "" RIVET_LDFLAGS "${RIVET_LDFLAGS}")
  # message(STATUS "RIVET_LDFLAGS=${RIVET_LDFLAGS}")

  EXECUTE_PROCESS (COMMAND ${RIVET_BIN_DIR}/rivet-config --libs
    OUTPUT_VARIABLE RIVET_LIBS)
  string (REGEX REPLACE "\n$" "" RIVET_LIBS "${RIVET_LIBS}")
  # message(STATUS "RIVET_LIBS=${RIVET_LIBS}")

  #set(RIVET_LIBRARIES ${RIVET_LIBRARY})
  #  mark_as_advanced(Rivet_FOUND)

  string (REGEX REPLACE "\n$" "" RIVET_CXXFLAGS "${RIVET_CXXFLAGS}")
  string (REPLACE " " ";" RIVET_CXXFLAGS "${RIVET_CXXFLAGS}")

  find_library(RIVET_LIBRARY
    NAMES Rivet
    HINTS ${RIVET_LIB_DIR} ENV ${RIVET_ROOT}
  )

  add_library(Rivet SHARED IMPORTED)
  add_library(Rivet::Rivet ALIAS Rivet)
  set_target_properties(Rivet PROPERTIES IMPORTED_LOCATION
    ${RIVET_LIBRARY}
    INTERFACE_LINK_DIRECTORIES ${RIVET_LIB_DIR}
    INTERFACE_INCLUDE_DIRECTORIES ${RIVET_INCLUDE_DIR})
  set_target_properties(Rivet PROPERTIES IMPORTED_GLOBAL TRUE)
endif(RIVET_BIN_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Rivet
  REQUIRED_VARS RIVET_LIBRARY
  VERSION_VAR RIVET_VERSION
)

if (Rivet_FOUND)
  MESSAGE (STATUS
    "Compatible version of Rivet found: " ${RIVET_VERSION} " at " ${RIVET_PREFIX})
else(Rivet_FOUND)
  MESSAGE (STATUS
    "No compatible version of Rivet found")
endif(Rivet_FOUND)
#
# EOF
#


