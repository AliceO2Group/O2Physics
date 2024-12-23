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

o2physics_add_dpl_workflow(femtouniverse-pair-track-track
          SOURCES femtoUniversePairTaskTrackTrack.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-pair-track-track-mc-truth
          SOURCES femtoUniversePairTaskTrackTrackMcTruth.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-pair-track-track-extended
          SOURCES femtoUniversePairTaskTrackTrackExtended.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-pair-track-track-mult-kt-extended
          SOURCES femtoUniversePairTaskTrackTrackMultKtExtended.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-pair-track-nucleus
          SOURCES femtoUniversePairTaskTrackNucleus.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-pair-track-track-threedrelmom-mult-kt-extended
          SOURCES femtoUniversePairTaskTrackTrack3DMultKtExtended.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(o2-analysis-cf-femtouniverse-pair-track-track-spherhar-mult-kt-extended
          SOURCES femtoUniversePairTaskTrackTrackSpherHarMultKtExtended.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-debug-track
          SOURCES femtoUniverseDebugTrack.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-pair-track-v0-extended
          SOURCES femtoUniversePairTaskTrackV0Extended.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-pair-track-cascade-extended
          SOURCES femtoUniversePairTaskTrackCascadeExtended.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-pair-track-d0
          SOURCES femtoUniversePairTaskTrackD0.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-pair-track-phi
          SOURCES femtoUniversePairTaskTrackPhi.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-debug-v0
          SOURCES femtoUniverseDebugV0.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-hash
          SOURCES femtoUniverseHashTask.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-pair-track-track-m-c
          SOURCES femtoUniversePairTaskTrackTrackMC.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-efficiency-task
          SOURCES femtoUniverseEfficiencyTask.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtouniverse-efficiency-base
          SOURCES femtoUniverseEfficiencyBase.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_executable(femtouniverse-cutculator
          SOURCES femtoUniverseCutCulator.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)