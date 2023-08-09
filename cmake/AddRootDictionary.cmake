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

configure_file(${CMAKE_CURRENT_LIST_DIR}/rootcling_wrapper.sh.in
               ${CMAKE_BINARY_DIR}/rootcling_wrapper.sh @ONLY)

#
# add_root_dictionary generates one dictionary to be added to a target.
#
# Besides the dictionary source itself a pcm file is also created.
# It will be installed alongside the target's library file
#
# arguments :
#
# * 1st parameter (required) is the target the dictionary should be added to
#
# * HEADERS (required) is a list of relative filepaths needed for the dictionary
#   definition
#
# * LINKDEF (required) is a single relative filepath to the LINKDEF file needed
#   by rootcling.
#
# LINKDEF and HEADERS must contain relative paths only (relative to the
# CMakeLists.txt that calls this add_root_dictionary function).
#
# The target must be of course defined _before_ calling this function (i.e.
# add_library(target ...) has been called).
#
# In addition :
#
# * target_include_directories _must_ have be called as well, in order to be
#   able to compute the list of include directories needed to _compile_ the
#   dictionary
#
# Note also that the generated dictionary is added to PRIVATE SOURCES list of
# the target.
#
function(add_root_dictionary target)
  cmake_parse_arguments(PARSE_ARGV
                        1
                        A
                        ""
                        "LINKDEF"
                        "HEADERS;BASENAME")
  if(A_UNPARSED_ARGUMENTS)
    message(
      FATAL_ERROR "Unexpected unparsed arguments: ${A_UNPARSED_ARGUMENTS}")
  endif()

  if(A_BASENAME)
    message(STATUS "BASENAME parameter is deprecated. Will be ignored")
  endif()

  set(required_args "LINKDEF;HEADERS")
  foreach(required_arg IN LISTS required_args)
    if(NOT A_${required_arg})
      message(FATAL_ERROR "Missing required argument: ${required_arg}")
    endif()
  endforeach()

  # check all given filepaths are relative ones
  foreach(h ${A_HEADERS} ${A_LINKDEF})
    if(IS_ABSOLUTE ${h})
      message(
        FATAL_ERROR "add_root_dictionary only accepts relative paths, but the"
                    "following path is absolute : ${h}")
    endif()
  endforeach()

  # convert all relative paths to absolute ones. LINKDEF must be the last one.
  foreach(h ${A_HEADERS} ${A_LINKDEF})
    get_filename_component(habs ${CMAKE_CURRENT_LIST_DIR}/${h} ABSOLUTE)
    list(APPEND headers ${habs})
  endforeach()

  # check all given filepaths actually exist
  foreach(h ${headers})
    get_filename_component(f ${h} ABSOLUTE)
    if(NOT EXISTS ${f})
      message(
        FATAL_ERROR
          "add_root_dictionary was given an inexistant input include ${f}")
    endif()
  endforeach()

  # Generate the pcm files alongside the library
  get_property(lib_output_dir
               TARGET ${target}
               PROPERTY LIBRARY_OUTPUT_DIRECTORY)
  if(NOT lib_output_dir)
    set(lib_output_dir ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  # Define the names of generated files
  get_property(basename TARGET ${target} PROPERTY OUTPUT_NAME)
  if(NOT basename)
    set(basename ${target})
  endif()
  set(dictionary G__${basename})
  set(dictionaryFile ${CMAKE_CURRENT_BINARY_DIR}/${dictionary}.cxx)
  set(pcmBase ${dictionary}_rdict.pcm)
  set(pcmFile ${lib_output_dir}/${pcmBase})
  set(o2physics_targetPCMMAP_TARGET "${o2physics_targetPCMMAP_TARGET};${target}" CACHE INTERNAL "target/PCM map (target)")
  set(o2physics_targetPCMMAP_PCM "${o2physics_targetPCMMAP_PCM};${pcmFile}" CACHE INTERNAL "target/PCM map (pcm)")

  # get the list of compile_definitions
  set(prop $<TARGET_PROPERTY:${target},COMPILE_DEFINITIONS>)

  # Build the LD_LIBRARY_PATH required to get rootcling running fine
  #
  # Need at least root core library
  get_filename_component(LD_LIBRARY_PATH ${ROOT_Core_LIBRARY} DIRECTORY)
  # and possibly toolchain libs if we are using a toolchain
  if(DEFINED ENV{GCC_TOOLCHAIN_ROOT})
    set(LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:$ENV{GCC_TOOLCHAIN_ROOT}/lib")
    set(LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:$ENV{GCC_TOOLCHAIN_ROOT}/lib64")
  endif()

  set(includeDirs $<TARGET_PROPERTY:${target},INCLUDE_DIRECTORIES>)

  # add a custom command to generate the dictionary using rootcling
  # cmake-format: off
  add_custom_command(
    OUTPUT ${dictionaryFile} ${pcmFile}
    VERBATIM
    COMMAND
    ${CMAKE_BINARY_DIR}/rootcling_wrapper.sh
      --dictionary_file ${dictionaryFile}
      --ld_library_path ${LD_LIBRARY_PATH}
      --include_dirs -I$<JOIN:${includeDirs},$<SEMICOLON>-I>
      $<$<BOOL:${prop}>:--compile_defs>
      $<$<BOOL:${prop}>:-D$<JOIN:${prop},$<SEMICOLON>-D>>
      --pcmdeps "$<REMOVE_DUPLICATES:${list_pcm_deps_${target}}>"
      --headers "${headers}"
    COMMAND
    ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/${pcmBase} ${pcmFile}
    DEPENDS ${headers} "$<REMOVE_DUPLICATES:${list_pcm_deps_${target}}>")
  # cmake-format: on

  # add dictionary source to the target sources
  target_sources(${target} PRIVATE ${dictionaryFile})

  get_property(libs TARGET ${target} PROPERTY INTERFACE_LINK_LIBRARIES)
  if(NOT ROOT::RIO IN_LIST libs)
    # add ROOT::IO if not already there as a target that has a Root dictionary
    # has to depend on ... Root
    target_link_libraries(${target} PUBLIC ROOT::RIO)
  endif()

  # Get the list of include directories that will be required to compile the
  # dictionary itself and add them as private include directories
  foreach(h IN LISTS A_HEADERS)
    if(IS_ABSOLUTE ${h})
      message(FATAL_ERROR "Path ${h} should be relative, not absolute")
    endif()
    get_filename_component(a ${h} ABSOLUTE)
    string(REPLACE "${h}" "" d "${a}")
    list(APPEND dirs ${d})
  endforeach()
  list(REMOVE_DUPLICATES dirs)
  target_include_directories(${target} PRIVATE ${dirs})

  # will install the pcm files alongside the target's lib
  get_filename_component(dict ${dictionaryFile} NAME_WE)
  install(FILES ${pcmFile} DESTINATION ${CMAKE_INSTALL_LIBDIR})

endfunction()
