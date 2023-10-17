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

add_subdirectory(QC)

o2physics_add_dpl_workflow(hyperon-reco-test
                    SOURCES hyperon-reco-test.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nuclei-batask
                    SOURCES LFNucleiBATask.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mc-spectra-efficiency
                    SOURCES mcspectraefficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(spectra-tof
                    SOURCES spectraTOF.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(spectra-tof-run2
                    SOURCES spectraTOFRun2.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(spectra-tpc
                    SOURCES spectraTPC.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(spectra-tpc-pikapr
                    SOURCES spectraTPCPiKaPr.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(spectra-tpc-tiny
                    SOURCES spectraTPCtiny.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(spectra-tpc-tiny-pikapr
                    SOURCES spectraTPCtinyPiKaPr.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzeroanalysis
                    SOURCES lambdakzeroanalysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzeroanalysis-mc
                    SOURCES lambdakzeroanalysisMC.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadeanalysis
                    SOURCES cascadeanalysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadeanalysismc
                    SOURCES cascadeanalysisMC.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypertritonanalysis
                    SOURCES hypertritonAnalysis.cxx
                    PUBLIC_LINK_LIBRARIES O2::DetectorsBase O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(spectra-charged
                    SOURCES spectraCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(id-raa
                    SOURCES identifiedraa.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(rsnanalysis
                    SOURCES rsnanalysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phianalysis
                    SOURCES phianalysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(k892analysis
                    SOURCES k892analysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambda1520analysis
                    SOURCES lambda1520analysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nuclei-hist
                    SOURCES NucleiHistTask.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-hist
                    SOURCES QAHistTask.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascpostprocessing
                    SOURCES cascpostprocessing.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hstrangecorrelation
                    SOURCES hStrangeCorrelation.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stqa
                    SOURCES stqa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::ReconstructionDataFormats O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(non-prompt-cascade
                    SOURCES nonPromptCascade.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::ReconstructionDataFormats O2Physics::AnalysisCore O2::DetectorsBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hyhefour-analysis
                    SOURCES hyhe4analysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(k1analysis
                    SOURCES k1analysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phianalysisrun3
                    SOURCES phianalysisrun3.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(f0980analysis
                    SOURCES f0980analysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0postprocessing
                    SOURCES v0postprocessing.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(antimatter-abs-hmpid
                    SOURCES AntimatterAbsorptionHMPID.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2::ReconstructionDataFormats O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadecorrelations
                    SOURCES cascadecorrelations.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(kinkanalysis
                    SOURCES kinkAnalysis.cxx
                    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambda1520spherocityanalysis
                    SOURCES lambda1520SpherocityAnalysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypertriton3bodyanalysis
                   SOURCES hypertriton3bodyanalysis.cxx
                   PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
                   COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypertriton3bodymcqa
                   SOURCES hypertriton3bodyMCQA.cxx
                   PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
                   COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(vzero-cascade-absorption
                    SOURCES vzero_cascade_absorption.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(delta-analysis
                    SOURCES deltaanalysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

 o2physics_add_dpl_workflow(rhoanalysis
                    SOURCES rhoanalysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(rsn-analysis-thnsparse
                    SOURCES rsnanalysisTHnSparse.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
