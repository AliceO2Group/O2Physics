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

function(o2physics_link_all_subdirs mainDir currentDir installDir)
  # Find the relative subdirectory of currentDir wrt mainDir
  # and create for each of the subdirectory levels a symlink in installDir

  # Find the relative subdirectory
  string(REPLACE "${mainDir}" "" relDir "${currentDir}")
  # Splite the relative directory to a list of subdirectories
  string(REPLACE "/" ";" listDirs "${relDir}")
  foreach(dir ${listDirs})
    # Link each subdirectory inside include
    # This helps to use #include "A/B/C/x.h" from the flat installation of header files
    message("Linking ${dir} in ${installDir}")
    execute_process(
      COMMAND ln -sfn . ${dir}
      WORKING_DIRECTORY ${installDir}
      )
  endforeach()
endfunction()
