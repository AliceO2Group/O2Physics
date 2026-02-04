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

add_subdirectory(Converters)

o2physics_add_dpl_workflow(dgcand-producer
                           SOURCES DGCandProducer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGCutparHolder  O2Physics::AnalysisCCDB O2Physics::EventFilteringUtils
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sgcand-producer
                           SOURCES SGCandProducer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::SGCutParHolder O2Physics::AnalysisCCDB
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dgbccand-producer
                           SOURCES DGBCCandProducer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGCutparHolder
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upccand-producer
                           SOURCES UPCCandidateProducer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::UPCCutparHolder
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tau-event-table-producer
                           SOURCES tauEventTableProducer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(two-tracks-event-table-producer
                           SOURCES twoTracksEventTableProducer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(fwd-track-propagation
                           SOURCES fwdTrackPropagation.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::GlobalTracking
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(udmcparticles-to-udtracks
                           SOURCES udMcParticles2udTracks.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(udmccollisions-to-udcollisions
                           SOURCES udMcCollisions2udCollisions.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tau-three-prong-event-table-producer
                           SOURCES tauThreeProngEventTableProducer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-cand-producer-muon
                           SOURCES upcCandProducerMuon.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::UPCCutparHolder O2::GlobalTracking
                           COMPONENT_NAME Analysis)
