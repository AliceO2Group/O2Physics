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

# tableReader_withAssoc.cxx is 5000 lines of task definitions, and seven
# workflow variants below #include it verbatim -- differing only in which
# adaptAnalysisTask<> entries they list, with identical preprocessor state and
# identical link libraries. Without this it is parsed and instantiated seven
# times, and it is the heaviest thing in the PWG: a cc1plus compiling one of
# these was OOM-killed at ~8 GB RSS.
#
# Skipped under recc, which caches compilations remotely instead.
add_library(TableReaderWithAssocPCH OBJECT tableReaderWithAssocPCH.cxx)
target_link_libraries(TableReaderWithAssocPCH PUBLIC O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::PWGDQCore O2Physics::MLCore)
if(NOT DEFINED ENV{USE_RECC})
  target_precompile_headers(TableReaderWithAssocPCH PRIVATE
    [["PWGDQ/Tasks/tableReader_withAssoc.cxx"]])
endif()

o2physics_add_dpl_workflow(table-reader
                    SOURCES tableReader.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(table-reader-with-assoc
                    SOURCES tableReader_withAssoc_workflowSpec.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::PWGDQCore O2Physics::MLCore
                    REUSE_FROM TableReaderWithAssocPCH
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(table-reader-with-assoc-same-event-pairing-barrel-only
                    SOURCES tableReader_withAssoc_SameEventPairingBarrelOnly_workflowSpec.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::PWGDQCore O2Physics::MLCore
                    REUSE_FROM TableReaderWithAssocPCH
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(table-reader-with-assoc-same-event-pairing-muon-only
                    SOURCES tableReader_withAssoc_SameEventPairingMuonOnly_workflowSpec.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::PWGDQCore O2Physics::MLCore
                    REUSE_FROM TableReaderWithAssocPCH
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(table-reader-with-assoc-same-event-pairing
                    SOURCES tableReader_withAssoc_SameEventPairing_workflowSpec.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::PWGDQCore O2Physics::MLCore
                    REUSE_FROM TableReaderWithAssocPCH
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(table-reader-with-assoc-asymmetric-pairing
                    SOURCES tableReader_withAssoc_AsymmetricPairing_workflowSpec.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::PWGDQCore O2Physics::MLCore
                    REUSE_FROM TableReaderWithAssocPCH
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(table-reader-with-assoc-dilepton-track-pairing
                    SOURCES tableReader_withAssoc_DileptonTrackPairing_workflowSpec.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::PWGDQCore O2Physics::MLCore
                    REUSE_FROM TableReaderWithAssocPCH
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(table-reader-with-assoc-dilepton-track-track-pairing
                    SOURCES tableReader_withAssoc_DileptonTrackTrackPairing_workflowSpec.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::PWGDQCore O2Physics::MLCore
                    REUSE_FROM TableReaderWithAssocPCH
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(table-reader-with-assoc-direct
                    SOURCES tableReader_withAssoc_direct.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::PWGDQCore O2Physics::MLCore O2::ReconstructionDataFormats O2::DetectorsCommonDataFormats O2::DetectorsVertexing O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(efficiency
                    SOURCES dqEfficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(efficiency-with-assoc
                    SOURCES dqEfficiency_withAssoc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(efficiency-with-assoc-direct
                    SOURCES dqEfficiency_withAssoc_direct.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore O2::ReconstructionDataFormats O2::DetectorsCommonDataFormats O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(energy-correlator-direct
                    SOURCES dqEnergyCorrelator_direct.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(filter-pp
                    SOURCES filterPP.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(filter-pp-with-association
                    SOURCES filterPPwithAssociation.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(filter-pb-pb
                    SOURCES filterPbPb.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore O2Physics::SGCutParHolder
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0-selector
                    SOURCES v0selector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DCAFitter O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dalitz-selection
                    SOURCES DalitzSelection.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow
                    SOURCES dqFlow.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-muon-mch-trk-eff
                    SOURCES taskMuonMchTrkEfficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-j-psi-hf
                    SOURCES taskJpsiHf.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-muon-dca
                    SOURCES muonDCA.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlation
                    SOURCES dqCorrelation.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-mch-align-record
                    SOURCES mchAlignRecord.cxx
                    PUBLIC_LINK_LIBRARIES
                        O2::Framework
                        O2Physics::AnalysisCore
                        O2Physics::PWGDQCore
                        O2::CommonUtils
                        O2::MCHClustering
                        O2::DPLUtils
                        O2::CCDB
                        O2::DataFormatsParameters
                        O2::MCHBase
                        O2::MCHTracking
                        O2::DataFormatsMCH
                        O2::DetectorsBase
                        O2::MCHGeometryTransformer
                        O2::MathUtils
                        O2::MCHAlign
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-muon-mid-eff
                    SOURCES MIDefficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::MIDBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-fwd-track-pid
                    SOURCES taskFwdTrackPid.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(quarkonia-to-hyperons
                    SOURCES quarkoniaToHyperons.cxx
                    PUBLIC_LINK_LIBRARIES O2::DetectorsBase O2::Framework O2::DCAFitter KFParticle::KFParticle O2Physics::AnalysisCore O2Physics::MLCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(model-converter-mult-pv
                    SOURCES ModelConverterMultPv.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(model-converter-event-extended
                    SOURCES ModelConverterEventExtended.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(model-converter-mc-reduced-event
                    SOURCES ModelConverterReducedMCEvents.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tag-and-probe
                    SOURCES TagAndProbe.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-matching
                    SOURCES qaMatching.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mft-mch-matcher
                    SOURCES mftMchMatcher.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(muon-global-alignment
                    SOURCES muonGlobalAlignment.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore O2::MCHGeometryTransformer
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(j-psi-muon-correlations
                    SOURCES jPsiMuonCorrelations.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore O2::MCHGeometryTransformer
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(global-muon-matcher
                    SOURCES global-muon-matcher.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore O2::MCHGeometryTransformer
                    COMPONENT_NAME Analysis)
