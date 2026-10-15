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

if(FastJet_FOUND)
o2physics_add_library(PWGJECore
               SOURCES  FastJetUtilities.cxx
                        JetFinder.cxx
                        JetBkgSubUtils.cxx
                        emcalCrossTalkEmulation.cxx
               PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore FastJet::FastJet FastJet::Contrib ONNXRuntime::ONNXRuntime O2::EMCALBase O2::EMCALReconstruction)

o2physics_target_root_dictionary(PWGJECore
              HEADERS JetFinder.h
                      JetUtilities.h
                      FastJetUtilities.h
                      JetTaggingUtilities.h
                      JetBkgSubUtils.h
                      emcalCrossTalkEmulation.h
                      utilsTrackMatchingEMC.h
              LINKDEF PWGJECoreLinkDef.h)
endif()

# Precompiled header shared across PWGJE. Percentages are the share of the
# group's workflow translation units that include the header, measured over the tree.
# Jet.h reaches 90% of the group; the existing substructure carriers already
# precompile it, so this one serves everything they do not cover.
#
# runDataProcessing.h is deliberately NOT precompiled: it requires each workflow's
# customize() to be declared before it is included, so precompiling it silently
# unregisters every option (the binary aborts at --dump-workflow, but links fine).
#
# Skipped under recc, which caches compilations remotely instead.
add_library(JetDataModelPCH OBJECT jetDataModelPCH.cxx)
target_link_libraries(JetDataModelPCH PUBLIC O2::Framework O2Physics::AnalysisCore O2Physics::PWGJECore)
if(NOT DEFINED ENV{USE_RECC})
  target_precompile_headers(JetDataModelPCH PRIVATE
    [["PWGJE/DataModel/Jet.h"]]                # 90%
  )
endif()
