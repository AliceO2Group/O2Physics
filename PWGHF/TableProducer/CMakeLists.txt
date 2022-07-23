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

o2physics_add_dpl_workflow(track-index-skims-creator
                    SOURCES HFTrackIndexSkimsCreator.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pvrefit-skim-converted
                    SOURCES HFPvRefitSkimConverted.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-2prong
                    SOURCES HFCandidateCreator2Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-dstar
                    SOURCES candidateCreatorDstar.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-d0-tokpi
                    SOURCES HFTreeCreatorD0ToKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-cascade
                    SOURCES HFCandidateCreatorCascade.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-3prong
                    SOURCES HFCandidateCreator3Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-bplus
                    SOURCES HFCandidateCreatorBPlus.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-xicc
                    SOURCES HFCandidateCreatorXicc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-chic
                    SOURCES HFCandidateCreatorChic.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-lb
                    SOURCES HFCandidateCreatorLb.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-x
                    SOURCES HFCandidateCreatorX.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-lc-topkpi
                    SOURCES HFTreeCreatorLcToPKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-lb-tolcpi
                    SOURCES HFTreeCreatorLbToLcPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-x-tojpsipipi
                    SOURCES HFTreeCreatorXToJpsiPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-xicc-topkpipi
                    SOURCES HFTreeCreatorXiccToPKPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-bplus-tod0pi
                    SOURCES HFTreeCreatorBplusToD0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-chic-tojpsigamma
                   SOURCES HFTreeCreatorChicToJpsiGamma.cxx
                   PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                   COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(d0-candidate-selector
                    SOURCES HFD0CandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(d0-candidate-selector-alice3-barrel
                    SOURCES HFD0CandidateSelectorALICE3Barrel.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(d0-candidate-selector-parametrized-pid
                    SOURCES HFD0CandidateSelectorparametrizedPID.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lc-candidate-selector-parametrized-pid
                    SOURCES HFLcCandidateSelectorparametrizedPID.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lc-candidate-selector-alice3
                    SOURCES HFLcCandidateSelectorALICE3.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(d0-candidate-selector-alice3-forward
                    SOURCES HFD0CandidateSelectorALICE3Forward.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dplus-topikpi-candidate-selector
                    SOURCES HFDplusToPiKPiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(ds-tokkpi-candidate-selector
                    SOURCES HFDsToKKPiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lc-candidate-selector
                    SOURCES HFLcCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jpsi-candidate-selector
                    SOURCES HFJpsiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(xic-topkpi-candidate-selector
                    SOURCES HFXicToPKPiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(xicc-topkpipi-candidate-selector
                    SOURCES HFXiccToPKPiPiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lc-tok0sp-candidate-selector
                    SOURCES HFLcK0sPCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(bplus-tod0pi-candidate-selector
                    SOURCES HFBPlusToD0PiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lb-tolcpi-candidate-selector
                    SOURCES HFLbToLcPiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(x-tojpsipipi-candidate-selector
                    SOURCES HFXToJpsiPiPiCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(chic-tojpsigamma-candidate-selector
                    SOURCES HFChicToJpsiGammaCandidateSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)
