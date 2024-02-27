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

if (RIVET_FOUND)
  # o2physics_add_header_only_library(RivetOutput
  #   HEADERS RivetAOs.h)

  o2physics_add_library(RivetOutput
    SOURCES RivetAOs.cxx
    PUBLIC_LINK_LIBRARIES
    ROOT::Physics
    Rivet::Rivet
    Yoda::Yoda
    # O2::Mergers
  )

  o2physics_target_root_dictionary(RivetOutput
    HEADERS RivetAOs.h
    LINKDEF RivetOutputLinkDef.h)

  # Install header
  o2physics_name_target(RivetOutput NAME target)
  install(TARGETS ${target}
          # EXPORT O2PhysicsTargets
          FILE_SET RivetAOs.h
          INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PWGMM)

endif(RIVET_FOUND)

#
# EOF
#
