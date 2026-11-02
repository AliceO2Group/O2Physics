# Copyright 2019-2026 CERN and copyright holders of ALICE O2.
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
# o2physics_add_dpl_workflow(basename SOURCES ...) add a new dpl workflow executable.
# Besides building the executable, this will also generate a configuration json
# file via <executable-name> --dump so that it can be easily registered and
# deployed e.g. in the train infrastructure.
#
# For most of the options please see how o2physics_add_executable works.
#
# The installed executable will be named as regular executables (see
# o2physics_add_executable for details)
#
# The installed JSON file will be named <executable-name>.json and installed
# under share/dpl directory
#

function(o2physics_add_dpl_workflow baseTargetName)
  cmake_parse_arguments(PARSE_ARGV 1 A "" "COMPONENT_NAME;TARGETVARNAME;REUSE_FROM"
                        "SOURCES;PUBLIC_LINK_LIBRARIES")

  if(A_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Got trailing arguments ${A_UNPARSED_ARGUMENTS}")
  endif()

  o2physics_add_executable(${baseTargetName}
    COMPONENT_NAME ${A_COMPONENT_NAME} TARGETVARNAME targetExeName
    SOURCES ${A_SOURCES}
    PUBLIC_LINK_LIBRARIES O2::Framework ${A_PUBLIC_LINK_LIBRARIES})

  if(A_TARGETVARNAME)
    set(${A_TARGETVARNAME}
        ${targetExeName}
        PARENT_SCOPE)
  endif()
  set_property(TARGET ${targetExeName} PROPERTY JOB_POOL_COMPILE analysis)
  set_property(TARGET ${targetExeName} PROPERTY JOB_POOL_LINK analysis)

  # Reuse a precompiled header. Without an explicit REUSE_FROM, fall back to the
  # shared AnalysisPCH (Common/Core): a workflow translation unit spends most of
  # its time parsing the framework headers, and there are ~1500 of them, so the
  # default is worth more than the handful of targets that name their own.
  #
  # Set O2PHYSICS_DEFAULT_PCH to an empty string to opt out globally, e.g. when
  # bisecting a PCH-related build failure.
  if(NOT DEFINED O2PHYSICS_DEFAULT_PCH)
    set(O2PHYSICS_DEFAULT_PCH AnalysisPCH)
  endif()
  set(_pch "${A_REUSE_FROM}")
  if(NOT _pch)
    set(_pch "${O2PHYSICS_DEFAULT_PCH}")
  endif()
  # A target cannot reuse its own PCH, and the carrier is not built when recc
  # is caching compilations remotely instead.
  if(_pch AND NOT _pch STREQUAL targetExeName AND NOT DEFINED ENV{USE_RECC})
    # GCC refuses a PCH built with a different preprocessor state than the
    # consumer's, and -Werror turns that refusal into a build failure:
    #   cmake_pch.hxx.gch: not used because `RANS_ENABLE_JSON' not defined
    # The carrier links O2Physics::AnalysisCore, which reaches O2::rANS and its
    # INTERFACE -DRANS_ENABLE_JSON, while a workflow that links only
    # O2::Framework (the converters, the tutorials) never sees it.
    #
    # Copying COMPILE_DEFINITIONS fixed that case and revealed another: CI then
    # failed on `_REENTRANT' not defined, on a compile line that DID carry
    # -DRANS_ENABLE_JSON -- so the copy was working and simply does not reach
    # far enough. Whatever supplies _REENTRANT arrives at the carrier as
    # something other than a compile definition (-pthread, which travels in
    # INTERFACE_COMPILE_OPTIONS, is the likely route), so adding options alone
    # would only move the goalposts to whichever kind of usage requirement goes
    # missing next.
    #
    # $<COMPILE_ONLY:> applies a target's *compile* usage requirements --
    # definitions, options, include directories, features -- without placing it
    # on the link line or creating a link dependency, which is the constraint
    # that ruled out simply linking the carrier. It transfers the whole
    # preprocessor state the carrier compiled with, and that state is exactly
    # what GCC compares, rather than one property of it at a time.
    # Requires CMake >= 3.27.
    target_link_libraries(${targetExeName} PRIVATE $<COMPILE_ONLY:${_pch}>)
    target_precompile_headers(${targetExeName} REUSE_FROM ${_pch})
  endif()

  set(jsonFile $<TARGET_FILE_BASE_NAME:${targetExeName}>.json)

  add_custom_command(
    TARGET ${targetExeName} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E env "LD_LIBRARY_PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}:$$LD_LIBRARY_PATH" $<TARGET_FILE:${targetExeName}> -b --dump-workflow --dump-workflow-file ${jsonFile})
  add_dependencies(${targetExeName} O2::FrameworkAnalysisSupport)

  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${jsonFile}
    DESTINATION ${CMAKE_INSTALL_DATADIR}/dpl)

endfunction()
