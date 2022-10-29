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
                    SOURCES trackIndexSkimCreator.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pvrefit-skim-converted
                    SOURCES refitPVDummy.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-2prong
                    SOURCES candidateCreator2Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-dstar
                    SOURCES candidateCreatorDstar.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-d0-tokpi
                    SOURCES treeCreatorD0ToKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-cascade
                    SOURCES candidateCreatorCascade.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-3prong
                    SOURCES candidateCreator3Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-b0
                    SOURCES candidateCreatorB0.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-bplus
                    SOURCES candidateCreatorBplus.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-xicc
                    SOURCES candidateCreatorXicc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-chic
                    SOURCES candidateCreatorChic.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-lb
                    SOURCES candidateCreatorLb.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-x
                    SOURCES candidateCreatorX.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-lc-topkpi
                    SOURCES treeCreatorLcToPKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-lb-tolcpi
                    SOURCES treeCreatorLbToLcPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-x-tojpsipipi
                    SOURCES treeCreatorXToJpsiPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-xicc-topkpipi
                    SOURCES treeCreatorXiccToPKPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-bplus-tod0pi
                    SOURCES treeCreatorBplusToD0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-chic-tojpsigamma
                   SOURCES treeCreatorChicToJpsiGamma.cxx
                   PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing ROOT::EG
                   COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(d0-candidate-selector
                    SOURCES candidateSelectorD0.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(d0-candidate-selector-alice3-barrel
                    SOURCES candidateSelectorD0Alice3Barrel.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(d0-candidate-selector-parametrized-pid
                    SOURCES candidateSelectorD0ParametrizedPid.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lc-candidate-selector-parametrized-pid
                    SOURCES candidateSelectorLcParametrizedPid.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lc-candidate-selector-alice3
                    SOURCES candidateSelectorLcAlice3.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(d0-candidate-selector-alice3-forward
                    SOURCES candidateSelectorD0Alice3Forward.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dplus-topikpi-candidate-selector
                    SOURCES candidateSelectorDplusToPiKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(ds-tokkpi-candidate-selector
                    SOURCES candidateSelectorDsToKKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lc-candidate-selector
                    SOURCES candidateSelectorLc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jpsi-candidate-selector
                    SOURCES candidateSelectorJpsi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(xic-topkpi-candidate-selector
                    SOURCES candidateSelectorXicToPKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(xicc-topkpipi-candidate-selector
                    SOURCES candidateSelectorXiccToPKPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lc-tok0sp-candidate-selector
                    SOURCES candidateSelectorLcK0sP.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(b0-todpi-candidate-selector
                    SOURCES candidateSelectorB0ToDPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(bplus-tod0pi-candidate-selector
                    SOURCES candidateSelectorBplusToD0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lb-tolcpi-candidate-selector
                    SOURCES candidateSelectorLbToLcPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(x-tojpsipipi-candidate-selector
                    SOURCES candidateSelectorXToJpsiPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(chic-tojpsigamma-candidate-selector
                    SOURCES candidateSelectorChicToJpsiGamma.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)
