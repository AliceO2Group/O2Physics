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

# add KFParticle::KFParticle as library to targets depending on KFParticle


#if (NOT DEFINED KFParticle_DIR)
#        set(KFParticle_DIR "$ENV{KFPARTICLE_ROOT}")
#endif()

find_path(KFPARTICLE_INCLUDE_DIR KFParticle.h
        PATH_SUFFIXES "include" "include/KFParticle"
        HINTS "$ENV{KFPARTICLE_ROOT}")
find_library(KFPARTICLE_LIBPATH "KFParticle"
        PATH_SUFFIXES "lib"
        HINTS "$ENV{KFPARTICLE_ROOT}")

mark_as_advanced(KFPARTICLE_LIBPATH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(KFParticle DEFAULT_MSG
                                  KFPARTICLE_LIBPATH KFPARTICLE_INCLUDE_DIR)

if(KFPARTICLE_FOUND)
   set(KFPARTICLE_LIBRARIES ${KFPARTICLE_LIBPATH})
   set(KFPARTICLE_INCLUDE_DIRS ${KFPARTICLE_INCLUDE_DIR})
   # add target
   if(NOT TARGET KFParticle::KFParticle)
      add_library(KFParticle::KFParticle IMPORTED INTERFACE)
      set_target_properties(KFParticle::KFParticle PROPERTIES
        INTERFACE_LINK_LIBRARIES "${KFPARTICLE_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${KFPARTICLE_INCLUDE_DIR}")
   endif()
endif()
