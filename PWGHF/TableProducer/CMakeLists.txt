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

o2physics_add_dpl_workflow(hf-track-index-skims-creator
                    SOURCES HFTrackIndexSkimsCreator.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-candidate-creator-2prong
                    SOURCES HFCandidateCreator2Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-tree-creator-d0-tokpi
                    SOURCES HFTreeCreatorD0ToKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-candidate-creator-cascade
                    SOURCES HFCandidateCreatorCascade.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-candidate-creator-3prong
                    SOURCES HFCandidateCreator3Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-candidate-creator-bplus
                    SOURCES HFCandidateCreatorBPlus.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-candidate-creator-xicc
                    SOURCES HFCandidateCreatorXicc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-candidate-creator-x
                    SOURCES HFCandidateCreatorX.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-tree-creator-lc-topkpi
                    SOURCES HFTreeCreatorLcToPKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-tree-creator-x-tojpsipipi
                    SOURCES HFTreeCreatorXToJpsiPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-tree-creator-xicc-topkpipi
                    SOURCES HFTreeCreatorXiccToPKPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-d0-candidate-selector
                    SOURCES HFD0CandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-dplus-topikpi-candidate-selector
                    SOURCES HFDplusToPiKPiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-lc-candidate-selector
                    SOURCES HFLcCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-jpsi-candidate-selector
                    SOURCES HFJpsiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-xic-topkpi-candidate-selector
                    SOURCES HFXicToPKPiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-xicc-topkpipi-candidate-selector
                    SOURCES HFXiccToPKPiPiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-lc-tok0sp-candidate-selector
                    SOURCES HFLcK0sPCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-bplus-tod0pi-candidate-selector
                    SOURCES HFBPlusToD0PiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-x-tojpsipipi-candidate-selector
                    SOURCES HFXToJpsiPiPiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)
