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

add_subdirectory(PID)
add_subdirectory(V0Cascades)

o2physics_add_dpl_workflow(qa-trackselection
    SOURCES qaTrackSelection.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-event-track
    SOURCES qaEventTrack.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::ReconstructionDataFormats
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-event-track-lite
    SOURCES qaEventTrackLite.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::ReconstructionDataFormats
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-event-track-lite-producer
    SOURCES qaEventTrackLiteProducer.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::ReconstructionDataFormats
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-fakehits
    SOURCES qaFakeHits.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-efficiency
    SOURCES qaEfficiency.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-dca-mc
    SOURCES qaDcaMC.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-sign-charge-mc
    SOURCES qaSignChargeMC.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-match-eff
    SOURCES qaMatchEff.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-impact-parameter
    SOURCES qaImpPar.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    O2::ReconstructionDataFormats
    O2::DetectorsCommonDataFormats
    O2::DetectorsVertexing
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(monitorfilterbit
    SOURCES MonitorFilterBit.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-prim-vtx-vs-time
    SOURCES qaPrimVtxVsTime.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-tracksplitting
    SOURCES qaTrackSplitting.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(unit-test-for-reconstruction
    SOURCES unitTestForReconstruction.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tag-and-probe-dmesons
    SOURCES tagAndProbeDmesons.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing O2Physics::MLCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derived-data-creator-d0-calibration
    SOURCES derivedDataCreatorD0Calibration.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCore O2::DCAFitter O2Physics::MLCore
    COMPONENT_NAME Analysis)

