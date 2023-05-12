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

include(O2PhysicsNameTarget)

#
# o2physics_add_library(baseTargetName SOURCES c1.cxx c2.cxx .....) defines a new
# target of type "library" composed of the given sources. It also defines an
# alias named O2Physics::baseTargetName. The generated library will be called
# libO2Physics[baseTargetName].(dylib|so|.a) (for exact naming see the
# o2physics_name_target function). For each source c1.cxx a header c1.h is installed
# if it exists in the same directory.
#
# The library will be static or shared depending on the BUILD_SHARED_LIBS option
# (which is normally ON for O2 project)
#
# Parameters:
#
# * SOURCES (required) : the list of source files to compile into this library
#
# * PUBLIC_LINK_LIBRARIES (needed in most cases) : the list of targets this
#   library depends on (e.g. ROOT::Hist, O2::CommonConstants). It is mandatory
#   to use the fully qualified target name (i.e. including the namespace part)
#   even for internal (O2) targets.
#
# * INSTALL_HEADERS (not needed in most cases): the list of additional headers
#   which should be installed with the library. Not needed for each source
#   c1.cxx where the header c1.h is found in the same folder. Those are installed
#   automatically.
#
# * PUBLIC_INCLUDE_DIRECTORIES (not needed in most cases) : the list of include
#   directories where to find the include files needed to compile this library
#   and that will be needed as well by the consumers of that library. By default
#   the current source directory is taken into
#   account, which should cover most of the use cases. Use this parameter only
#   for special cases then. Note that if you do specify this parameter it
#   replaces the default, it does not add to them.
#
# * PRIVATE_LINK_LIBRARIES (not needed in most cases) : the list of targets this
#   library needs at compile time (i.e. those dependencies won't be propagated
#   to the targets depending on this one). It is mandatory to use the fully
#   qualified target name (i.e. including the namespace part) even for internal
#   (O2) targets.
#
# * PRIVATE_INCLUDE_DIRECTORIES (not needed in most cases) : the list of include
#   directories where to find the include files needed to compile this library,
#   but that will _not_ be needed by its consumers. But default we add the
#   ${CMAKE_CURRENT_BINARY_DIR} here to cover use case of generated headers
#   (e.g. by protobuf). Note that if you do specify this parameter it replaces
#   the default, it does not add to them.
#
function(o2physics_add_library baseTargetName)

  cmake_parse_arguments(
    PARSE_ARGV
    1
    A
    ""
    "TARGETVARNAME"
    "SOURCES;PUBLIC_INCLUDE_DIRECTORIES;PUBLIC_LINK_LIBRARIES;PRIVATE_INCLUDE_DIRECTORIES;PRIVATE_LINK_LIBRARIES;INSTALL_HEADERS"
    )

  if(A_UNPARSED_ARGUMENTS)
    message(
      FATAL_ERROR "Unexpected unparsed arguments: ${A_UNPARSED_ARGUMENTS}")
  endif()

  o2physics_name_target(${baseTargetName} NAME targetName)
  set(target ${targetName})

  # If -DSTANDALONE_EXECUTABLES is passed to cmake,
  # we only build object libraries so that each executable
  # is fully selfcontained for what concerns O2Physics
  # code and can be deployed with a simple cp, reusing
  # a compatible O2 environment
  if (STANDALONE_EXECUTABLES)
    add_library(${target} OBJECT ${A_SOURCES})
  else()
    add_library(${target} ${A_SOURCES})
  endif()
  # define the target and its O2Physics:: alias
  add_library(O2Physics::${baseTargetName} ALIAS ${target})

  # set the export name so that packages using O2Physics can reference the target as
  # O2Physics::${baseTargetName} as well (assuming the export is installed with
  # namespace O2Physics::)
  set_property(TARGET ${target} PROPERTY EXPORT_NAME ${baseTargetName})

  # output name of the lib will be libO2Physics[baseTargetName].(so|dylib|a)
  set_property(TARGET ${target} PROPERTY OUTPUT_NAME O2Physics${baseTargetName})

  if(A_TARGETVARNAME)
    set(${A_TARGETVARNAME} ${target} PARENT_SCOPE)
  endif()

  # Start by adding the public dependencies to other targets
  if(A_PUBLIC_LINK_LIBRARIES)
    foreach(L IN LISTS A_PUBLIC_LINK_LIBRARIES)
      string(FIND ${L} "::" NS)
      if(${NS} EQUAL -1)
        message(FATAL_ERROR "Trying to use a non-namespaced target ${L}")
      endif()
      target_link_libraries(${target} PUBLIC ${L})
    endforeach()
  endif()

  # Then add the private dependencies to other targets
  if(A_PRIVATE_LINK_LIBRARIES)
    foreach(L IN LISTS A_PRIVATE_LINK_LIBRARIES)
      string(FIND ${L} "::" NS)
      if(${NS} EQUAL -1)
        message(FATAL_ERROR "Trying to use a non-namespaced target ${L}")
      endif()
      target_link_libraries(${target} PRIVATE ${L})
    endforeach()
  endif()
  # set the public include directories if available
  if(A_PUBLIC_INCLUDE_DIRECTORIES)
    foreach(d IN LISTS A_PUBLIC_INCLUDE_DIRECTORIES)
      get_filename_component(adir ${d} ABSOLUTE)
      if(NOT IS_DIRECTORY ${adir})
        message(
          FATAL_ERROR "Trying to append non existing include directory ${d}")
      endif()
      target_include_directories(${target} PUBLIC $<BUILD_INTERFACE:${adir}>)
    endforeach()
  else()
    # default
    target_include_directories(
      ${target}
      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    # add top level directory so e.g. the #include "Common/Core/xxx.h" will work
    target_include_directories(
      ${target}
      PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>)
  endif()

  # set the private include directories if available
  if(A_PRIVATE_INCLUDE_DIRECTORIES)
    foreach(d IN LISTS A_PRIVATE_INCLUDE_DIRECTORIES)
      get_filename_component(adir ${d} ABSOLUTE)
      if(NOT IS_DIRECTORY ${adir})
        message(
          FATAL_ERROR "Trying to append non existing include directory ${d}")
      endif()
      target_include_directories(${target} PRIVATE $<BUILD_INTERFACE:${d}>)
    endforeach()
  else()
    # use sane(?) default
    target_include_directories(
      ${target}
      PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
    get_filename_component(adir ${CMAKE_CURRENT_LIST_DIR}/src ABSOLUTE)
    if(EXISTS ${adir})
      target_include_directories(
        ${target}
        PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src>)
    endif()
  endif()

  # The INCLUDES DESTINATION adds ${CMAKE_INSTALL_INCLUDEDIR} to the
  # INTERFACE_INCLUDE_DIRECTORIES property
  #
  # The EXPORT must come first in the list of parameters
  #
  install(TARGETS ${target}
          EXPORT O2PhysicsTargets
          INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

  # Install header for each source (if exists)
  if(A_SOURCES)
    foreach(d IN LISTS A_SOURCES)
      # Replace .cxx -> .h
      string(REGEX REPLACE "[.]cxx$" ".h" OUTPUT_HEADER ${d})
      get_filename_component(OUTPUT_HEADER_ABS ${OUTPUT_HEADER} ABSOLUTE)
      if(EXISTS ${OUTPUT_HEADER_ABS})
        install(FILES ${OUTPUT_HEADER}
                DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
      endif()
    endforeach()
  endif()

  # install additional headers
  if(A_INSTALL_HEADERS)
    install(FILES ${CMAKE_CURRENT_LIST_DIR}/${A_INSTALL_HEADERS}
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  endif()

  # Link subdirectories
  install(
    SCRIPT ${CMAKE_SOURCE_DIR}/cmake/O2PhysicsLinkAllSubDirs.cmake
    CODE " o2physics_link_all_subdirs(${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}) "
  )

endfunction()
