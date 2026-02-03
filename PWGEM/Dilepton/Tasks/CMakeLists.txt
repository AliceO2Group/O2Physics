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


o2physics_add_dpl_workflow(efficiency-ee
                    SOURCES emEfficiencyEE.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lmee-lf-cocktail
                  SOURCES lmeeLFCocktail.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lmee-hf-cocktail
                  SOURCES lmeeHFCocktail.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mc-templates
                    SOURCES MCtemplates.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(smearing
                    SOURCES smearing.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(create-resolution-map
                    SOURCES createResolutionMap.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::GlobalTracking
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(test-bremsstrahlung
                    SOURCES testBremsstrahlung.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(table-reader-barrel
                    SOURCES tableReaderBarrel.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-qc
                    SOURCES eventQC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(single-electron-qc
                    SOURCES singleElectronQC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2Physics::PWGEMDileptonCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(single-electron-qc-mc
                    SOURCES singleElectronQCMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(single-muon-qc
                    SOURCES singleMuonQC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGEMDileptonCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(single-muon-qc-mc
                    SOURCES singleMuonQCMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dielectron
                    SOURCES dielectron.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2Physics::PWGEMDileptonCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dielectron-mc
                    SOURCES dielectronMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dimuon
                    SOURCES dimuon.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGEMDileptonCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dimuon-mc
                    SOURCES dimuonMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(prefilter-dielectron
                    SOURCES prefilterDielectron.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::MLCore O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(associate-mccollision-to-collision
                    SOURCES associateMCcollision.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(study-mc-truth
                    SOURCES studyMCTruth.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(matching-mft
                    SOURCES matchingMFT.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::GlobalTracking
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tagging-hfe
                    SOURCES taggingHFE.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dielectron-hadron-mpc
                    SOURCES dielectronHadronMPC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2Physics::PWGEMDileptonCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dimuon-hadron-mpc
                    SOURCES dimuonHadronMPC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2Physics::PWGEMDileptonCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mc-particle-predictions-otf
                    SOURCES mcParticlePredictionsOTF.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(study-dcafitter
                    SOURCES studyDCAFitter.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(evaluate-acceptance
                    SOURCES evaluateAcceptance.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dilepton-polarization
                    SOURCES dileptonPolarization.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(check-mc-template
                    SOURCES checkMCTemplate.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::GlobalTracking
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(check-mc-pair-template
                    SOURCES checkMCPairTemplate.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2Physics::PWGEMDileptonCore O2Physics::PWGEMDileptonUtils
                    COMPONENT_NAME Analysis)
