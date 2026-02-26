# Copyright 2019-2026 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.
add_library(JetSubstructureHFPCH OBJECT jetSubstructureHFPCH.cxx)
target_link_libraries(JetSubstructureHFPCH
  PUBLIC O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore)
target_precompile_headers(JetSubstructureHFPCH PRIVATE
  [["PWGJE/DataModel/Jet.h"]]
  [["PWGJE/DataModel/JetSubstructure.h"]]
  [["PWGJE/DataModel/JetSubtraction.h"]]
  [["PWGJE/Tasks/jetSubstructureHF.h"]]
  <Framework/ASoA.h>
  <Framework/AnalysisTask.h>
  <Framework/ConfigContext.h>
  <Framework/DataProcessorSpec.h>
  <Framework/runDataProcessing.h>

  <vector>
  <memory>
)

o2physics_add_dpl_workflow(emc-cellmonitor
                    SOURCES emcCellMonitor.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(emc-clustermonitor
                    SOURCES emcClusterMonitor.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(emc-eventselection-qa
                    SOURCES emcEventSelectionQA.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(emc-vertexselection-qa
                    SOURCES emcVertexSelectionQA.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(emcal-gamma-gamma-bc-wise
                    SOURCES emcalGammaGammaBcWise.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(emc-pi0-energyscale-calib
                    SOURCES emcalPi0EnergyScaleCalib.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(emc-tmmonitor
                    SOURCES emcTmMonitor.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(mc-generator-studies
                    SOURCES mcGeneratorStudies.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(photon-isolation-qa
                    SOURCES photonIsolationQA.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(photon-charged-trigger-correlation
                    SOURCES photonChargedTriggerCorrelation.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(task-emc-extensive-mc-qa
                    SOURCES taskEmcExtensiveMcQa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(photon-charged-trigger-producer
                    SOURCES photonChargedTriggerProducer.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(hadron-photon-correlation
										SOURCES hadronPhotonCorrelation.cxx
										PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
										COMPONENT_NAME Analysis)

if(FastJet_FOUND)
    o2physics_add_dpl_workflow(jet-background-analysis
                        SOURCES jetBackgroundAnalysis.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure
                        SOURCES jetSubstructure.cxx
                        REUSE_FROM JetSubstructureHFPCH
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-output
                        SOURCES jetSubstructureOutput.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-d0
                        SOURCES jetSubstructureD0.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        REUSE_FROM JetSubstructureHFPCH
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-dplus
                        SOURCES jetSubstructureDplus.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        REUSE_FROM JetSubstructureHFPCH
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-ds
                        SOURCES jetSubstructureDs.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        REUSE_FROM JetSubstructureHFPCH
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-dstar
                        SOURCES jetSubstructureDstar.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        REUSE_FROM JetSubstructureHFPCH
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-lc
                        SOURCES jetSubstructureLc.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        REUSE_FROM JetSubstructureHFPCH
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-b0
                        SOURCES jetSubstructureB0.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        REUSE_FROM JetSubstructureHFPCH
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-bplus
                        SOURCES jetSubstructureBplus.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        REUSE_FROM JetSubstructureHFPCH
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-xictoxipipi
                        SOURCES jetSubstructureXicToXiPiPi.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        REUSE_FROM JetSubstructureHFPCH
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-dielectron
                        SOURCES jetSubstructureDielectron.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        REUSE_FROM JetSubstructureHFPCH
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-d0-output
                        SOURCES jetSubstructureD0Output.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-dplus-output
                        SOURCES jetSubstructureDplusOutput.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-ds-output
                        SOURCES jetSubstructureDsOutput.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-dstar-output
                        SOURCES jetSubstructureDstarOutput.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-lc-output
                        SOURCES jetSubstructureLcOutput.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-b0-output
                        SOURCES jetSubstructureB0Output.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-bplus-output
                        SOURCES jetSubstructureBplusOutput.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-xictoxipipi-output
                        SOURCES jetSubstructureXicToXiPiPiOutput.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-dielectron-output
                        SOURCES jetSubstructureDielectronOutput.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-fragmentation
                        SOURCES jetFragmentation.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-v0-spectra
                        SOURCES v0JetSpectra.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-v0qa
                        SOURCES v0QA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-charged-qa
                        SOURCES jetFinderQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-outlier-qa
                        SOURCES jetOutlierQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-charged-v2
                        SOURCES jetChargedV2.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-d0-qa
                        SOURCES jetFinderD0QA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-dplus-qa
                        SOURCES jetFinderDplusQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-ds-qa
                        SOURCES jetFinderDsQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-dstar-qa
                        SOURCES jetFinderDstarQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-lc-qa
                        SOURCES jetFinderLcQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-b0-qa
                        SOURCES jetFinderB0QA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-bplus-qa
                        SOURCES jetFinderBplusQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-xictoxipipi-qa
                        SOURCES jetFinderXicToXiPiPiQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-dielectron-qa
                        SOURCES jetFinderDielectronQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-full-qa
                        SOURCES jetFinderFullQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-v0-qa
                        SOURCES jetFinderV0QA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-spectra-charged
                        SOURCES jetSpectraCharged.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(charged-jet-hadron
                        SOURCES chargedJetHadron.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(trigger-correlations
                        SOURCES triggerCorrelations.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-trigger-charged-qa
                        SOURCES jetTriggerChargedQa.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-full-trigger-qa
                        SOURCES fullJetTriggerQATask.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-matching-qa
                        SOURCES jetMatchingQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-validation-qa
                        SOURCES jetValidationQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-tutorial
                        SOURCES jetTutorial.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-tutorial-skeleton
                        SOURCES jetTutorialSkeleton.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(track-jet-qa
                        SOURCES trackJetQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(track-efficiency
                        SOURCES trackEfficiency.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-hadron-recoil
                        SOURCES jetHadronRecoil.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(recoil-jets
                        SOURCES recoilJets.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-nsubjettiness
                        SOURCES nsubjettiness.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(phi-in-jets
                        SOURCES phiInJets.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(nuclei-in-jets
                        SOURCES nucleiInJets.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                        O2Physics::AnalysisCCDB
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-taggerhf-qa
                        SOURCES jetTaggerHFQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-lund-reclustering
                        SOURCES jetLundReclustering.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore FastJet::FastJet FastJet::Contrib
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(hf-fragmentation-function
                        SOURCES hfFragmentationFunction.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-planarflow
                        SOURCES jetPlanarFlow.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-ch-corr
                        SOURCES jetChCorr.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(bjet-tree-creator
                        SOURCES bjetTreeCreator.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(stat-prompt-photon
                        SOURCES statPromptPhoton.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(full-jet-spectra
                        SOURCES fullJetSpectra.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                        O2Physics::AnalysisCCDB
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(bjet-tagging-ml
                        SOURCES bjetTaggingML.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2Physics::MLCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-spectra-ese
                        SOURCES jetSpectraEseTask.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-cross-section-efficiency
                        SOURCES jetCrossSectionEfficiency.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(gamma-jet-tree-producer
                        SOURCES gammaJetTreeProducer.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(dijet-finder-charged-qa
                        SOURCES dijetFinderQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(bjet-tagging-gnn
                        SOURCES bjetTaggingGnn.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                        O2Physics::AnalysisCCDB
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-shape
                        SOURCES jetShape.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-formationtimereclustering
                        SOURCES jetFormationTimeReclustering.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(hf-debug
                        SOURCES hfDebug.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-debug
                        SOURCES jetDebug.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(substructure-debug
                        SOURCES substructureDebug.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-ds-spectrum-subs
                        SOURCES jetDsSpectrumAndSubstructure.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(bjet-cent-mult
                        SOURCES bjetCentMult.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
endif()
