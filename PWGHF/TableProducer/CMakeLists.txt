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

o2physics_add_dpl_workflow(track-index-skim-creator
                    SOURCES trackIndexSkimCreator.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsVertexing O2::DCAFitter O2Physics::AnalysisCore O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(refit-pv-dummy
                    SOURCES refitPvDummy.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-2prong
                    SOURCES candidateCreator2Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-dstar
                    SOURCES candidateCreatorDstar.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-d0-to-k-pi
                    SOURCES treeCreatorD0ToKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-cascade
                    SOURCES candidateCreatorCascade.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-3prong
                    SOURCES candidateCreator3Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-b0
                    SOURCES candidateCreatorB0.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-bplus
                    SOURCES candidateCreatorBplus.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-xicc
                    SOURCES candidateCreatorXicc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-to-xi-pi
                    SOURCES candidateCreatorToXiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-sigmac0plusplus
                    SOURCES candidateCreatorSigmac0plusplus.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-chic
                    SOURCES candidateCreatorChic.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-lb
                    SOURCES candidateCreatorLb.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-x
                    SOURCES candidateCreatorX.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-lc-to-p-k-pi
                    SOURCES treeCreatorLcToPKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-lb-to-lc-pi
                    SOURCES treeCreatorLbToLcPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-x-to-jpsi-pi-pi
                    SOURCES treeCreatorXToJpsiPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-xicc-to-p-k-pi-pi
                    SOURCES treeCreatorXiccToPKPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-bplus-to-d0-pi
                    SOURCES treeCreatorBplusToD0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-chic-to-jpsi-gamma
                   SOURCES treeCreatorChicToJpsiGamma.cxx
                   PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                   COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-lc-to-k0s-p
                    SOURCES treeCreatorLcToK0sP.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-d0
                    SOURCES candidateSelectorD0.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-d0-alice3-barrel
                    SOURCES candidateSelectorD0Alice3Barrel.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-d0-parametrized-pid
                    SOURCES candidateSelectorD0ParametrizedPid.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc-parametrized-pid
                    SOURCES candidateSelectorLcParametrizedPid.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc-alice3
                    SOURCES candidateSelectorLcAlice3.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-d0-alice3-forward
                    SOURCES candidateSelectorD0Alice3Forward.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-dplus-to-pi-k-pi
                    SOURCES candidateSelectorDplusToPiKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-ds-to-k-k-pi
                    SOURCES candidateSelectorDsToKKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc
                    SOURCES candidateSelectorLc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc-ml
                    SOURCES candidateSelectorLcMl.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::CCDB ONNXRuntime::ONNXRuntime O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-jpsi
                    SOURCES candidateSelectorJpsi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-xic-to-p-k-pi
                    SOURCES candidateSelectorXicToPKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-xicc-to-p-k-pi-pi
                    SOURCES candidateSelectorXiccToPKPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-to-xi-pi
                    SOURCES candidateSelectorToXiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc-to-k0s-p
                    SOURCES candidateSelectorLcToK0sP.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-b0-to-d-pi
                    SOURCES candidateSelectorB0ToDPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-bplus-to-d0-pi
                    SOURCES candidateSelectorBplusToD0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lb-to-lc-pi
                    SOURCES candidateSelectorLbToLcPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-x-to-jpsi-pi-pi
                    SOURCES candidateSelectorXToJpsiPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-chic-to-jpsi-gamma
                    SOURCES candidateSelectorChicToJpsiGamma.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
