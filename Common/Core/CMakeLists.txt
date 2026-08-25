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

o2physics_add_library(AnalysisCore
        SOURCES TrackSelection.cxx
        OrbitRange.cxx
        PID/ParamBase.cxx
        PID/PIDTOF.cxx
        PID/PIDTOFParamService.cxx
        CollisionAssociation.cxx
        TrackSelectionDefaults.cxx
        EventPlaneHelper.cxx
        TableHelper.cxx
        MetadataHelper.cxx
        CollisionTypeHelper.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2::DataFormatsParameters ROOT::EG O2::CCDB ROOT::Physics O2::FT0Base O2::FV0Base O2::DataFormatsParamTOF)

# Precompiled header shared by DPL analysis workflows. Parsing the framework
# headers dominates the cost of a workflow translation unit, and there are
# ~1500 of them; o2physics_add_dpl_workflow reuses this by default.
#
# Skipped under recc, which caches compilations remotely instead -- the two do
# not combine, as a PCH is a local artefact of one compiler invocation.
#
# The reusing targets must compile with the same preprocessor state as this
# carrier, which does not follow from the link list: AnalysisCore reaches
# O2::rANS and its INTERFACE -DRANS_ENABLE_JSON, whereas a bare workflow links
# only O2::Framework. o2physics_add_dpl_workflow therefore copies this target's
# compile definitions onto every consumer; keep that in mind when changing what
# is linked here.
add_library(AnalysisPCH OBJECT analysisPCH.cxx)
target_link_libraries(AnalysisPCH PUBLIC O2::Framework O2Physics::AnalysisCore)
if(NOT DEFINED ENV{USE_RECC})
  target_precompile_headers(AnalysisPCH PRIVATE
    <Framework/ASoA.h>
    <Framework/AnalysisTask.h>
    <Framework/ConfigContext.h>
    <Framework/DataProcessorSpec.h>
    # NOT runDataProcessing.h. It states the contract itself: "we need to
    # declare the customize method before include this file". Precompiling it
    # makes it the first include in every workflow, so a source's customize() is
    # never seen, its options are never registered, and the binary aborts at
    # --dump-workflow with "missing option: ...". The link succeeds, so this
    # surfaces as a runtime failure with no hint that a PCH caused it.

    <memory>
    <string>
    <vector>
  )
endif()

o2physics_target_root_dictionary(AnalysisCore
        HEADERS TrackSelection.h
        TrackSelectionDefaults.h
        OrbitRange.h
        TableHelper.h
        MetadataHelper.h
        EventPlaneHelper.h
        PID/ParamBase.h
        PID/DetectorResponse.h
        PID/PIDTOF.h
        PID/TPCPIDResponse.h
        PID/PIDTOFParamService.h
        CollisionTypeHelper.h
        LINKDEF AnalysisCoreLinkDef.h)

o2physics_add_library(EventFilteringUtils
    SOURCES Zorro.cxx ZorroSummary.cxx
    INSTALL_HEADERS ZorroHelper.h ZorroSummary.h
    PUBLIC_LINK_LIBRARIES O2::Framework O2::CCDB ROOT::EG O2::CCDB ROOT::Physics Arrow::arrow_shared)

o2physics_target_root_dictionary(EventFilteringUtils
   HEADERS ZorroHelper.h ZorroSummary.h
   LINKDEF EventFilteringUtilsLinkDef.h)

o2physics_add_header_only_library(TPCDriftManager
        HEADERS TPCVDriftManager.h
        INTERFACE_LINK_LIBRARIES)

o2physics_add_library(FFitWeights
        SOURCES FFitWeights.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2::CCDB ROOT::Physics)

o2physics_target_root_dictionary(FFitWeights
        HEADERS FFitWeights.h
        LINKDEF FFitWeightsLinkDef.h)

