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

o2physics_add_dpl_workflow(qa-event-track
        SOURCES qaEventTrack.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::ReconstructionDataFormats
        COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-event-track-lite
        SOURCES qaEventTrackLite.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::ReconstructionDataFormats
        COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-efficiency
        SOURCES qaEfficiency.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
        COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-match-eff
        SOURCES qaMatchEff.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
        COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-impact-parameter
        SOURCES qaImpPar.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework
        O2Physics::AnalysisCore
        O2::ReconstructionDataFormats
        O2::DetectorsBase
        O2::DetectorsCommonDataFormats
        O2::DetectorsVertexing
        COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-k0s-tracking-efficiency
                    SOURCES qaK0sTrackingEfficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# PID QA

## TOF
o2physics_add_dpl_workflow(pid-tof-qa
        SOURCES qaPIDTOF.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
        COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-tof-qa-mc
        SOURCES qaPIDTOFMC.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
        COMPONENT_NAME Analysis)

## TPC
o2physics_add_dpl_workflow(pid-tpc-qa
        SOURCES qaPIDTPC.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
        COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-tpc-qa-mc
        SOURCES qaPIDTPCMC.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
        COMPONENT_NAME Analysis)
##FDD
  o2physics_add_dpl_workflow(fdd-qa
          SOURCES qaFDD.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase COMPONENT_NAME Analysis)


add_subdirectory(TPC)