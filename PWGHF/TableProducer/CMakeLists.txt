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

# Skimming

o2physics_add_dpl_workflow(track-index-skim-creator
                    SOURCES trackIndexSkimCreator.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing O2::DCAFitter O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(refit-pv-dummy
                    SOURCES refitPvDummy.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# Candidate creators

o2physics_add_dpl_workflow(candidate-creator-2prong
                    SOURCES candidateCreator2Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter KFParticle::KFParticle
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-3prong
                    SOURCES candidateCreator3Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-b0
                    SOURCES candidateCreatorB0.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-bplus
                    SOURCES candidateCreatorBplus.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-bs
                    SOURCES candidateCreatorBs.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-cascade
                    SOURCES candidateCreatorCascade.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-dstar
                    SOURCES candidateCreatorDstar.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-lb
                    SOURCES candidateCreatorLb.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-sigmac0plusplus
                    SOURCES candidateCreatorSigmac0plusplus.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-to-xi-pi
                    SOURCES candidateCreatorToXiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-xicc
                    SOURCES candidateCreatorXicc.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

# Candidate selectors

o2physics_add_dpl_workflow(candidate-selector-b0-to-d-pi
                    SOURCES candidateSelectorB0ToDPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-bplus-to-d0-pi
                    SOURCES candidateSelectorBplusToD0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-bs-to-ds-pi
                    SOURCES candidateSelectorBsToDsPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-d0
                    SOURCES candidateSelectorD0.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-dplus-to-pi-k-pi
                    SOURCES candidateSelectorDplusToPiKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-ds-to-k-k-pi
                    SOURCES candidateSelectorDsToKKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lb-to-lc-pi
                    SOURCES candidateSelectorLbToLcPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc
                    SOURCES candidateSelectorLc.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc-ml
                    SOURCES candidateSelectorLcMl.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc-to-k0s-p
                    SOURCES candidateSelectorLcToK0sP.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-to-xi-pi
                    SOURCES candidateSelectorToXiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-xic-to-p-k-pi
                    SOURCES candidateSelectorXicToPKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-xicc-to-p-k-pi-pi
                    SOURCES candidateSelectorXiccToPKPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# Tree creators

o2physics_add_dpl_workflow(tree-creator-b0-to-d-pi
                    SOURCES treeCreatorB0ToDPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-bplus-to-d0-pi
                    SOURCES treeCreatorBplusToD0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-bs-to-ds-pi
                    SOURCES treeCreatorBsToDsPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-d0-to-k-pi
                    SOURCES treeCreatorD0ToKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-dplus-to-pi-k-pi
                    SOURCES treeCreatorDplusToPiKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-ds-to-k-k-pi
                    SOURCES treeCreatorDsToKKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-lb-to-lc-pi
                    SOURCES treeCreatorLbToLcPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-lc-to-k0s-p
                    SOURCES treeCreatorLcToK0sP.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-lc-to-p-k-pi
                    SOURCES treeCreatorLcToPKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-to-xi-pi
                    SOURCES treeCreatorToXiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-xic-to-p-k-pi
                    SOURCES treeCreatorXicToPKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-xicc-to-p-k-pi-pi
                    SOURCES treeCreatorXiccToPKPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)


