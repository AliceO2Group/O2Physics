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

# o2physics_add_header_only_library creates a header-only target.
#
# * INCLUDE_DIRECTORIES the relative path(s) to the headers of that library if
#   not specified will be set as "include" simply (which should work just fine
#   in most cases)
#
function(o2physics_add_header_only_library baseTargetName)

  cmake_parse_arguments(PARSE_ARGV
                        1
                        A
                        ""
                        ""
                        "INCLUDE_DIRECTORIES;INTERFACE_LINK_LIBRARIES;HEADERS")

  if(A_UNPARSED_ARGUMENTS)
    message(
      FATAL_ERROR "Unexpected unparsed arguments: ${A_UNPARSED_ARGUMENTS}")
  endif()

  o2physics_name_target(${baseTargetName} NAME target)

  # define the target and its O2:: alias
  add_library(${target} INTERFACE)
  add_library(O2Physics::${baseTargetName} ALIAS ${target})

  # set the export name so that packages using O2Physics can reference the target as
  # O2Physics::${baseTargetName} as well (assuming the export is installed with
  # namespace O2Physics::)
  set_property(TARGET ${target} PROPERTY EXPORT_NAME ${baseTargetName})

  if(NOT A_INCLUDE_DIRECTORIES)
    get_filename_component(dir include ABSOLUTE)
    if(EXISTS ${dir})
      set(A_INCLUDE_DIRECTORIES $<BUILD_INTERFACE:${dir}>)
    else()
      set(A_INCLUDE_DIRECTORIES $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>)
    endif()
  endif()

  # specify only the BUILD_INTERFACE directories here.
  # the INSTALL_INTERFACE is taken care of by the
  # install(TARGETS ... EXPORT ... INCLUDES DESTINATION) below
  target_include_directories(
    ${target}
    INTERFACE ${A_INCLUDE_DIRECTORIES})

  if(A_INTERFACE_LINK_LIBRARIES)
    target_link_libraries(${target} INTERFACE ${A_INTERFACE_LINK_LIBRARIES})
  endif()

  if(A_HEADERS)
    target_sources(${target} INTERFACE FILE_SET HEADERS FILES ${A_HEADERS})
  endif()

  # All the directories given after
  # INCLUDES DESTINATION are added to the INTERFACE_INCLUDE_DIRECTORIES
  # property of each installed target listed after TARGETS
  install(TARGETS ${target}
          EXPORT O2PhysicsTargets
          FILE_SET HEADERS
          INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  # Link subdirectories
  install(
    SCRIPT ${CMAKE_SOURCE_DIR}/cmake/O2PhysicsLinkAllSubDirs.cmake
    CODE " o2physics_link_all_subdirs(${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}) "
  )

endfunction()
