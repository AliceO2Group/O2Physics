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


o2physics_add_dpl_workflow(emc-cellmonitor
                    SOURCES emccellmonitor.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                     COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(emc-clustermonitor
                    SOURCES emcclustermonitor.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(emc-eventselection-qa
                    SOURCES emceventselectionqa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(emc-vertexselection-qa
                    SOURCES emcvertexselectionqa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(emc-pi0-energyscale-calib
                    SOURCES emcalPi0EnergyScaleCalib.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(emc-tmmonitor
                    SOURCES emctmmonitor.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(mc-generator-studies
                    SOURCES mcgeneratorstudies.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

if(FastJet_FOUND)
    o2physics_add_dpl_workflow(jet-substructure
                        SOURCES jetsubstructure.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-output
                        SOURCES jetsubstructureoutput.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-hf
                        SOURCES jetsubstructurehf.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-substructure-hf-output
                        SOURCES jetsubstructurehfoutput.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-fragmentation
                        SOURCES jetfragmentation.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-v0-spectra
                        SOURCES v0jetspectra.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-charged-qa
                        SOURCES jetfinderQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-hf-qa
                        SOURCES jetfinderhfQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-full-qa
                        SOURCES jetfinderfullQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-finder-v0-qa
                        SOURCES jetfinderv0QA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(trigger-correlations
                        SOURCES triggerCorrelations.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-charged-trigger-qa
                        SOURCES ChJetTriggerQATask.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
     o2physics_add_dpl_workflow(jet-full-trigger-qa
                         SOURCES FullJetTriggerQATask.cxx
                         PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                         COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-matching-qa
                        SOURCES jetmatchingqa.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-validation-qa
                        SOURCES jetvalidationqa.cxx
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
                        SOURCES trackJetqa.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-hadron-recoil
                        SOURCES jetHadronRecoil.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-nsubjettiness
                        SOURCES nSubjettiness.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(phi-in-jets
                        SOURCES phiInJets.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-taggerhf-qa
                        SOURCES jettaggerhfQA.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-lund-reclustering
                        SOURCES jetLundReclustering.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore FastJet::FastJet FastJet::Contrib
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-hf-fragmentation
                        SOURCES hffragmentationfunction.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-planarflow
                        SOURCES jetplanarflow.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
    o2physics_add_dpl_workflow(jet-ch-corr
                        SOURCES jetChCorr.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                        COMPONENT_NAME Analysis)
endif()
