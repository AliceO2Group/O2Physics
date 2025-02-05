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

o2physics_add_dpl_workflow(nuclei-batask
    SOURCES LFNucleiBATask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypertritonanalysis
    SOURCES hypertritonAnalysis.cxx
    PUBLIC_LINK_LIBRARIES O2::DetectorsBase O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nuclei-hist
    SOURCES NucleiHistTask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypertriton3bodyanalysis
    SOURCES hypertriton3bodyanalysis.cxx
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypertriton3bodymcqa
    SOURCES hypertriton3bodyMcqa.cxx
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore O2::TOFBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nuclei-in-jets
    SOURCES nucleiInJets.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(helium-flow
    SOURCES helium_flow.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(antimatter-abs-hmpid
    SOURCES AntimatterAbsorptionHMPID.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2::ReconstructionDataFormats O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hyhefour-analysis
    SOURCES hyhe4analysis.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(antid-lambda-ebye
    SOURCES antidLambdaEbye.cxx
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

o2physics_add_dpl_workflow(spectra-charged
    SOURCES spectraCharged.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(id-raa
    SOURCES identifiedraa.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-hist
    SOURCES QAHistTask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qc-spectra-tpc
    SOURCES QCspectraTPC.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(deut-rt-task
    SOURCES deutRtTask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hadronnucleicorrelation
    SOURCES hadronnucleicorrelation.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nuclei-efficiency-hist
    SOURCES NucleiEfficiencyTask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nuclei-ebye
    SOURCES nucleiEbye.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nuclei-toward-transv
    SOURCES nuclei_in_toward_transv_regions.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(ebye-mult
    SOURCES ebyeMult.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nuclei-from-hypertriton-map
    SOURCES nucleiFromHypertritonMap.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

if(FastJet_FOUND)
o2physics_add_dpl_workflow(angular-correlations-in-jets
    SOURCES AngularCorrelationsInJets.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib
    COMPONENT_NAME Analysis)
endif()
