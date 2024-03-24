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
# - Try to find Yoda
#
# This searches for the configuration script "yoda-config".  If this
# script is not found, then the package is not considered found.
#
# The targets
#
#  HepMC3
#  FastJet::Contrib
#
# _must_ previously have been found.
#
# Defines the target "Yoda" as well as the variables
#
#  Yoda_FOUND
#  YODA_VERSION
#  YODA_PREFIX
#  YODA_BIN_DIR
#  YODA_CPPFLAGS
#  YODA_CXXFLAGS
#  YODA_LDFLAGS
#  YODA_INCLUDE_DIR
#  YODA_LIBRARY
#  YODA_LIB_DIR

if (TARGET Yoda)
   # nothing to do, target is already there
   set(Yoda_FOUND TRUE)
   return()
endif()

# the following disables all default paths (either from cmake, from environment)
FIND_PATH (YODA_BIN_DIR yoda-config
  HINTS ${YODA_ROOT} ENV YODA_ROOT
)

if (YODA_BIN_DIR)
  MESSAGE("Executing ${YODA_BIN_DIR}/yoda-config")
  EXECUTE_PROCESS (COMMAND ${YODA_BIN_DIR}/yoda-config --version
    OUTPUT_VARIABLE YODA_VERSION
  )
  # Remove white-space
  STRING (REGEX REPLACE "[ \t\r\n]+" "" YODA_VERSION "${YODA_VERSION}")
  # Remove trailing patch
  STRING (REGEX REPLACE "-.*" "" YODA_VERSION "${YODA_VERSION}")
  message(STATUS "Found Yoda version: ${YODA_VERSION}")

  EXECUTE_PROCESS (COMMAND ${YODA_BIN_DIR}/yoda-config --prefix
    OUTPUT_VARIABLE YODA_PREFIX
  )
  # message(STATUS "Found Yoda prefix: ${YODA_PREFIX}")

  EXECUTE_PROCESS (COMMAND ${YODA_BIN_DIR}/yoda-config --includedir
    OUTPUT_VARIABLE YODA_INCLUDE_DIR)
  string (REGEX REPLACE "\n$" "" YODA_INCLUDE_DIR "${YODA_INCLUDE_DIR}")
  # message(STATUS "YODA_INCLUDE_DIR=${YODA_INCLUDE_DIR}")

  EXECUTE_PROCESS (COMMAND ${YODA_BIN_DIR}/yoda-config --libdir
    OUTPUT_VARIABLE YODA_LIB_DIR)
  string (REGEX REPLACE "\n$" "" YODA_LIB_DIR "${YODA_LIB_DIR}")
  # message(STATUS "YODA_LIB_DIR=${YODA_LIB_DIR}")

  EXECUTE_PROCESS (COMMAND ${YODA_BIN_DIR}/yoda-config --cppflags
    OUTPUT_VARIABLE YODA_CPPFLAGS)
  string (REGEX REPLACE "\n$" "" YODA_CPPFLAGS "${YODA_CPPFLAGS}")
  # message(STATUS "YODA_CPPFLAGS=${YODA_CPPFLAGS}")

  EXECUTE_PROCESS (COMMAND ${YODA_BIN_DIR}/yoda-config --cxxflags
    OUTPUT_VARIABLE YODA_CXXFLAGS)
  string (REGEX REPLACE "\n$" "" YODA_CXXFLAGS "${YODA_CXXFLAGS}")
  # message(STATUS "YODA_CXXFLAGS=${YODA_CXXFLAGS}")

  EXECUTE_PROCESS (COMMAND ${YODA_BIN_DIR}/yoda-config --ldflags
    OUTPUT_VARIABLE YODA_LDFLAGS)
  string (REGEX REPLACE "\n$" "" YODA_LDFLAGS "${YODA_LDFLAGS}")
  # message(STATUS "YODA_LDFLAGS=${YODA_LDFLAGS}")

  EXECUTE_PROCESS (COMMAND ${YODA_BIN_DIR}/yoda-config --libs
    OUTPUT_VARIABLE YODA_LIBS)
  string (REGEX REPLACE "\n$" "" YODA_LIBS "${YODA_LIBS}")
  # message(STATUS "YODA_LIBS=${YODA_LIBS}")

  #set(YODA_LIBRARIES ${YODA_LIBRARY})
  #  mark_as_advanced(Yoda_FOUND)

  string (REGEX REPLACE "\n$" "" YODA_CXXFLAGS "${YODA_CXXFLAGS}")
  string (REPLACE " " ";" YODA_CXXFLAGS "${YODA_CXXFLAGS}")

  find_library(YODA_LIBRARY
    NAMES YODA
    HINTS ${YODA_LIB_DIR} ENV ${YODA_ROOT}
  )

  add_library(Yoda SHARED IMPORTED)
  add_library(Yoda::Yoda ALIAS Yoda)
  set_target_properties(Yoda PROPERTIES IMPORTED_LOCATION
    ${YODA_LIBRARY}
    INTERFACE_LINK_DIRECTORIES ${YODA_LIB_DIR}
    INTERFACE_INCLUDE_DIRECTORIES ${YODA_INCLUDE_DIR})
  set_target_properties(Yoda PROPERTIES IMPORTED_GLOBAL TRUE)
endif(YODA_BIN_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Yoda
  REQUIRED_VARS YODA_LIBRARY
  VERSION_VAR YODA_VERSION
)

if (Yoda_FOUND)
  MESSAGE (STATUS
    "Compatible version of Yoda found: " ${YODA_VERSION} " at " ${YODA_PREFIX})
else(Yoda_FOUND)
  MESSAGE (STATUS
    "No compatible version of Yoda found")
endif(Yoda_FOUND)
#
# EOF
#


