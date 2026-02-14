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
add_subdirectory(PID)

o2physics_add_dpl_workflow(trackextension
                    SOURCES trackextension.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(trackselection
                    SOURCES trackselection.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-selection-service
                    SOURCES eventSelectionService.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCCDB
                                          O2::DataFormatsITSMFT
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-selection-service-test-ccdb
                    SOURCES eventSelectionServiceTestCCDB.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCCDB
                                          O2::DataFormatsITSMFT
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-selection-service-run2
                    SOURCES eventSelectionServiceRun2.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCCDB
                                          O2::DataFormatsITSMFT
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(multiplicity-table
                    SOURCES multiplicityTable.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(multcenttable
                    SOURCES multCentTable.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(multiplicity-extra-table
                    SOURCES multiplicityExtraTable.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(centrality-table
                    SOURCES centralityTable.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(timestamp
                    SOURCES timestamp.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(weak-decay-indices
                    SOURCES weakDecayIndices.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(ft0-corrected-table
                    SOURCES ft0CorrectedTable.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(track-propagation
                    SOURCES trackPropagation.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(propagationservice
                    SOURCES propagationService.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter KFParticle::KFParticle O2Physics::TPCDriftManager
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(propagationservice-run2
                    SOURCES propagationServiceRun2.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter KFParticle::KFParticle O2Physics::TPCDriftManager
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(track-dca-cov-filler-run2
                    SOURCES trackDcaCovFillerRun2.cxx
                    PUBLIC_LINK_LIBRARIES O2::DetectorsBase O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(calo-clusters
                    SOURCES caloClusterProducer.cxx
                    PUBLIC_LINK_LIBRARIES O2::DataFormatsPHOS O2::PHOSBase O2::PHOSReconstruction O2Physics::DataModel
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(fwdtrackextension
                    SOURCES fwdtrackextension.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                                          O2::DetectorsBase
                                          O2::DetectorsCommonDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(track-to-collision-associator
                    SOURCES trackToCollisionAssociator.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(fwdtrack-to-collision-associator
                    SOURCES fwdtrackToCollisionAssociator.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mccollisionextra
                    SOURCES mcCollsExtra.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qvector-table
                    SOURCES qVectorsTable.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework
                                          O2Physics::AnalysisCore
                                          O2::DetectorsCommonDataFormats
                                          O2Physics::DataModel
                                          O2::CCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mftmchmatchingml
                    SOURCES mftmchMatchingML.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework
                                          O2Physics::AnalysisCore
                                          ONNXRuntime::ONNXRuntime
                                          O2::CCDB O2Physics::MLCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(match-mft-ft0
                    SOURCES match-mft-ft0.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                                          O2::DetectorsBase O2::DetectorsCommonDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(zdc-extra-table-producer
                    SOURCES zdcExtraTableProducer.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(ese-table-producer
                    SOURCES eseTableProducer.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mftmch-matching-data
                    SOURCES match-mft-mch-data.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::GlobalTracking
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(occ-table-producer
                    SOURCES occupancyTableProducer.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mftmch-matching-data-mc
                    SOURCES match-mft-mch-data-mc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::GlobalTracking
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(muon-realignment
                    SOURCES muonRealignment.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2::DetectorsCommonDataFormats O2::MathUtils O2::MCHTracking O2::DataFormatsMCH O2::GlobalTracking O2::MCHBase O2::MCHGeometryTransformer O2::CommonUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(selectionstudytable
                    SOURCES selectionStudyTable.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(zdc-task-light-ions
                    SOURCES zdcTaskLightIons.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
