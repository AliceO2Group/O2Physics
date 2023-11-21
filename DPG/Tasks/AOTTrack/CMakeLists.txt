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

o2physics_add_dpl_workflow(qa-efficiency
                           SOURCES qaEfficiency.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-efficiency-v0s
                           SOURCES qaEfficiencyV0s.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
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

o2physics_add_dpl_workflow(qa-k0s-tracking-efficiency
                           SOURCES qaK0sTrackingEfficiency.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(perf-k0s-resolution
                           SOURCES perfK0sResolution.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-cascades
                           SOURCES qaCascades.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(monitorfilterbit
                           SOURCES MonitorFilterBit.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-prim-vtx-vs-time
                           SOURCES qaPrimVtxVsTime.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)
