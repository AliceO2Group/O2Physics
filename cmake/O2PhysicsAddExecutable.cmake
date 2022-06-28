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

include_guard()

#
# o2physics_add_executable(basename SOURCES ...) add an executable with the
# given sources.
#
# * SOURCES (required) gives the list of source files to compile into the
#   executable
# * PUBLIC_LINK_LIBRARIES (needed in most cases) indicates the list of targets
#   this executable depends on.
#
# The installed executable will be named o2physics[-exeType][-component_name]-basename
#
# where :
#
# * exeType is `test` if IS_TEST is set, `bench` if IS_BENCHMARK is set or void
#   otherwise
# * COMPONENT_NAME (optional) is typically used to indicate a subsystem name for
#   regular executables (e.g. o2physics-tpc-... or o2physics-mch-...) or the
#   origin target for tests (to help locate the source file in the source
#   directory hierarchy, e.g. o2physics-test-DataFormats-...)
#
# Note that the _target_ corresponding to the executable will be named
# O2Physicsexe[-exeType][component]-basename and can be retrieved with the
# TARGETVARNAME parameter if needed.
#
# For instance after a call to:
#
# o2physics_add_executable(toto SOURCES ... TARGETVARNAME titi IS_TEST)
#
# ${titi} will contain something like `O2Physicsexe-test-toto` (for the exact
# naming see the o2physics_name_target function) and an executable named
# o2-test-toto will be created upon build)

function(o2physics_add_executable baseTargetName)

  cmake_parse_arguments(PARSE_ARGV
                        1
                        A
                        "IS_TEST;NO_INSTALL;IS_BENCHMARK"
                        "COMPONENT_NAME;TARGETVARNAME"
                        "SOURCES;PUBLIC_LINK_LIBRARIES;JOB_POOL")

  if(A_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Got trailing arguments ${A_UNPARSED_ARGUMENTS}")
  endif()

  # check naming conventions for executable
  if (NOT ${baseTargetName} MATCHES "^[a-z0-9\-]*$")
    message(FATAL_ERROR "Executable name can only contain lower case letters, numbers and -. Violated by ${baseTargetName}")
  endif()

  # set the executable name following our coding convention
  if(A_IS_TEST)
    set(exeType -test)
  elseif(A_IS_BENCHMARK)
    set(exeType -bench)
  endif()

  if(A_COMPONENT_NAME)
    string(TOLOWER ${A_COMPONENT_NAME} component)
    set(comp -${component})
  endif()

  # Extract PWG name from folder path (if exists)
  # First get a relative source directory
  string(REPLACE ${PROJECT_SOURCE_DIR} "" RELATIVE_LIST_DIR ${CMAKE_CURRENT_LIST_DIR})
  # Match PWG
  string(REGEX MATCH "PWG[A-Z][A-Z]" PWG ${RELATIVE_LIST_DIR})
  if(PWG)
    string(REPLACE "PWG" "" PWGSUFFIX ${PWG})
    string(TOLOWER ${PWGSUFFIX} pwg_lower)
    set(pwg -${pwg_lower})
  endif()

  set(exeName o2${exeType}${comp}${pwg}-${baseTargetName})

  if(A_IS_TEST)
    set(isTest "IS_TEST")
  endif()
  if(A_IS_BENCH)
    set(isBench "IS_BENCH")
  endif()

  # get the target name. the convention might be different from the executable
  # convention.
  o2physics_name_target(${baseTargetName}
                        NAME
                        targetName
                        IS_EXE
                        ${isTest}
                        ${isBench}
                        PWG ${PWGSUFFIX})

  set(target ${targetName})

  if(A_TARGETVARNAME)
    set(${A_TARGETVARNAME} ${target} PARENT_SCOPE)
  endif()

  # add the executable with its sources
  add_executable(${target} ${A_SOURCES})

  # set the executable output name
  set_property(TARGET ${target} PROPERTY OUTPUT_NAME ${exeName})
  if(A_JOB_POOL)
    set_property(TARGET ${target} PROPERTY JOB_POOL_COMPILE ${A_JOB_POOL})
    set_property(TARGET ${target} PROPERTY JOB_POOL_LINK ${A_JOB_POOL})
  endif()

  if(A_IS_TEST)
    # tests go in a separate directory
    get_filename_component(outdir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/../tests
                           ABSOLUTE)
    set_property(TARGET ${target} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${outdir})
  endif()

  # use its dependencies
  foreach(lib IN LISTS A_PUBLIC_LINK_LIBRARIES)
    string(FIND ${lib} "::" NS)
    if(${NS} EQUAL -1)
      message(FATAL_ERROR "Trying to use a non-namespaced target ${lib}")
    endif()
    target_link_libraries(${target} PUBLIC ${lib})
  endforeach()

  if(NOT A_NO_INSTALL)
    # install the executable

    if(A_IS_TEST)
      install(TARGETS ${target}
              RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/${testsDir})
    else()
      install(TARGETS ${target} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()
  endif()

endfunction()
