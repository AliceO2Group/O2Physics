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
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing O2::DCAFitter O2Physics::AnalysisCCDB O2Physics::MLCore O2Physics::SGCutParHolder
                    COMPONENT_NAME Analysis)

# Helpers

o2physics_add_dpl_workflow(refit-pv-dummy
                    SOURCES refitPvDummy.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-creator
                    SOURCES pidCreator.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mc-pid-tof
                    SOURCES mcPidTof.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::TOFWorkflowUtils
                    COMPONENT_NAME Analysis)

# Candidate creators

o2physics_add_dpl_workflow(candidate-creator-2prong
                    SOURCES candidateCreator2Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter KFParticle::KFParticle O2Physics::SGCutParHolder O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-3prong
                    SOURCES candidateCreator3Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter KFParticle::KFParticle O2Physics::SGCutParHolder O2Physics::AnalysisCCDB
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
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-dstar
                    SOURCES candidateCreatorDstar.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-lb
                    SOURCES candidateCreatorLb.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-sigmac0plusplus
                    SOURCES candidateCreatorSigmac0plusplus.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-sigmac0plusplus-cascade
                    SOURCES candidateCreatorSigmac0plusplusCascade.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-xic0-omegac0
                    SOURCES candidateCreatorXic0Omegac0.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter KFParticle::KFParticle O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-xic-to-xi-pi-pi
                    SOURCES candidateCreatorXicToXiPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter KFParticle::KFParticle O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-xicc
                    SOURCES candidateCreatorXicc.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-mc-gen
                    SOURCES candidateCreatorMcGen.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# Candidate selectors

o2physics_add_dpl_workflow(candidate-selector-b0-to-d-pi
                    SOURCES candidateSelectorB0ToDPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-bplus-to-d0-pi
                    SOURCES candidateSelectorBplusToD0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
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

o2physics_add_dpl_workflow(candidate-selector-dstar-to-d0-pi
                    SOURCES candidateSelectorDstarToD0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lb-to-lc-pi
                    SOURCES candidateSelectorLbToLcPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc
                    SOURCES candidateSelectorLc.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc-pid-ml
                    SOURCES candidateSelectorLcPidMl.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc-to-k0s-p
                    SOURCES candidateSelectorLcToK0sP.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-omegac0-to-omega-ka
                    SOURCES candidateSelectorOmegac0ToOmegaKa.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-omegac0-to-omega-pi
                    SOURCES candidateSelectorOmegac0ToOmegaPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-omegac0-xic0-to-omega-ka
                    SOURCES candidateSelectorOmegac0Xic0ToOmegaKa.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-xic0-to-xi-pi-kf
                    SOURCES candidateSelectorXic0ToXiPiKf.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-to-xi-pi
                    SOURCES candidateSelectorToXiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-xic-to-p-k-pi
                    SOURCES candidateSelectorXicToPKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-xic-to-xi-pi-pi
                    SOURCES candidateSelectorXicToXiPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-xicc-to-p-k-pi-pi
                    SOURCES candidateSelectorXiccToPKPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

 o2physics_add_dpl_workflow(candidate-selector-cd
                     SOURCES candidateSelectorCd.cxx
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

o2physics_add_dpl_workflow(tree-creator-sigmac-corr-bkg
                    SOURCES treeCreatorSigmacCorrBkg.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-omegac-st
                    SOURCES treeCreatorOmegacSt.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-omegac0-to-omega-ka
                    SOURCES treeCreatorOmegac0ToOmegaKa.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-omegac0-to-omega-pi
                    SOURCES treeCreatorOmegac0ToOmegaPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-to-xi-pi
                    SOURCES treeCreatorToXiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-xic0-to-xi-pi-kf
                    SOURCES treeCreatorXic0ToXiPiKf.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-xic-to-p-k-pi
                    SOURCES treeCreatorXicToPKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-xic-to-xi-pi-pi
                    SOURCES treeCreatorXicToXiPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-xicc-to-p-k-pi-pi
                    SOURCES treeCreatorXiccToPKPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-dstar-to-d0-pi
                    SOURCES treeCreatorDstarToD0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-tcc-to-d0-d0-pi
                    SOURCES treeCreatorTccToD0D0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
# Derived-data creators

o2physics_add_dpl_workflow(derived-data-creator-b0-to-d-pi
                    SOURCES derivedDataCreatorB0ToDPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derived-data-creator-bplus-to-d0-pi
                    SOURCES derivedDataCreatorBplusToD0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derived-data-creator-d0-to-k-pi
                    SOURCES derivedDataCreatorD0ToKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derived-data-creator-dplus-to-pi-k-pi
                    SOURCES derivedDataCreatorDplusToPiKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derived-data-creator-ds-to-k-k-pi
                    SOURCES derivedDataCreatorDsToKKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derived-data-creator-dstar-to-d0-pi
                    SOURCES derivedDataCreatorDstarToD0Pi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derived-data-creator-lc-to-p-k-pi
                    SOURCES derivedDataCreatorLcToPKPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derived-data-creator-xic-to-xi-pi-pi
                    SOURCES derivedDataCreatorXicToXiPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# Converters

o2physics_add_dpl_workflow(converter-dstar-indices
                    SOURCES converterDstarIndices.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
