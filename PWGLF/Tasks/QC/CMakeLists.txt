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

o2physics_add_dpl_workflow(v0cascades-qa
    SOURCES v0cascadesqa.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(track-checks
    SOURCES trackchecks.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(resonanceqa
    SOURCES resonanceqa.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-qa
    SOURCES lfpidqa.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strarecostudy
    SOURCES straRecoStudy.cxx
    PUBLIC_LINK_LIBRARIES O2::DetectorsBase O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(kfperformancestudy
    SOURCES kfPerformanceStudy.cxx
    PUBLIC_LINK_LIBRARIES O2::DetectorsBase O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(its-tpc-matching-qa
    SOURCES lfITSTPCMatchingQA.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangeness-tracking-qc
    SOURCES strangenessTrackingQC.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2::ReconstructionDataFormats O2Physics::AnalysisCore O2::DetectorsBase
    COMPONENT_NAME Analysis)


o2physics_add_dpl_workflow(mc-signal-loss
    SOURCES mcSignalLoss.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lfpropstudy
    SOURCES lfpropStudy.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tpc-dedx-postcalibration
    SOURCES tpc_dEdx_postcalibration.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tpc-dedx-qa
    SOURCES tpc_dedx_qa.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(vertexqa
    SOURCES vertexQA.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(efficiencyqa
    SOURCES efficiencyQA.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangenessqcpp
    SOURCES strangenessQCPP.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangepidqa
    SOURCES strangepidqa.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stqa
    SOURCES stqa.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2::ReconstructionDataFormats O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(kfstrangenessstudy
    SOURCES kfStrangenessStudy.cxx
    PUBLIC_LINK_LIBRARIES O2::DetectorsBase O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(findable-study
    SOURCES findableStudy.cxx
    PUBLIC_LINK_LIBRARIES O2::DetectorsBase O2Physics::AnalysisCore O2Physics::v0SelectionGroup
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mcinelgt0
    SOURCES mcinelgt0.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tracked-cascade-properties
    SOURCES trackedCascadeProperties.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mc-particle-predictions
    SOURCES mcParticlePrediction.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strange-derived-qa
    SOURCES strangederivedqa.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(str-derived-genqa
    SOURCES strderivedGenQA.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(its-tpc-matching-vzeros
    SOURCES lfITSTPCMatchingSecondaryTracksQA.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0assoqa
    SOURCES v0assoqa.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter KFParticle::KFParticle O2Physics::TPCDriftManager
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(systematics-mapping
    SOURCES systematicsMapping.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)
