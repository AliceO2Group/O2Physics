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

o2physics_add_library(UDGoodRunSelector
                      SOURCES UDGoodRunSelector.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore ROOT::EG RapidJSON::RapidJSON)

o2physics_target_root_dictionary(UDGoodRunSelector
                                 HEADERS UDGoodRunSelector.h
                                 LINKDEF UDGoodRunSelectorLinkDef.h)

o2physics_add_library(UDFSParser
                      SOURCES UDFSParser.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore)

o2physics_target_root_dictionary(UDFSParser
                                 HEADERS UDFSParser.h
                                 LINKDEF UDFSParserLinkDef.h)

o2physics_add_library(DGCutparHolder
                      SOURCES DGCutparHolder.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore)

o2physics_target_root_dictionary(DGCutparHolder
                                 HEADERS DGCutparHolder.h
                                 LINKDEF DGCutparHolderLinkDef.h)

o2physics_add_library(SGCutParHolder
                      SOURCES SGCutParHolder.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore)

o2physics_target_root_dictionary(SGCutParHolder
                      HEADERS SGCutParHolder.h
                      LINKDEF SGCutParHolderLinkDef.h)

o2physics_add_library(DGPIDSelector
                      SOURCES DGPIDSelector.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore ROOT::EG)

o2physics_target_root_dictionary(DGPIDSelector
                                 HEADERS DGPIDSelector.h
                                 LINKDEF DGPIDSelectorLinkDef.h)

o2physics_add_library(UPCCutparHolder
                      SOURCES UPCCutparHolder.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore)

o2physics_target_root_dictionary(UPCCutparHolder
                                 HEADERS UPCCutparHolder.h
                                 LINKDEF UPCCutparHolderLinkDef.h)

o2physics_add_library(decayTree
                      SOURCES decayTree.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore ROOT::EG RapidJSON::RapidJSON)

o2physics_target_root_dictionary(decayTree
                                 HEADERS decayTree.h
                                 LINKDEF decayTreeLinkDef.h)

o2physics_add_library(FITCutParHolder
                      SOURCES FITCutParHolder.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore)

o2physics_target_root_dictionary(FITCutParHolder
                                 HEADERS FITCutParHolder.h
                                 LINKDEF FITCutParHolderLinkDef.h)
# Precompiled header shared across PWGUD. Percentages are the share of the
# group's workflow translation units that include the header, measured over the tree.
# UDTables.h is the group's single shared data model.
#
# runDataProcessing.h is deliberately NOT precompiled: it requires each workflow's
# customize() to be declared before it is included, so precompiling it silently
# unregisters every option (the binary aborts at --dump-workflow, but links fine).
#
# Skipped under recc, which caches compilations remotely instead.
add_library(UDTablesPCH OBJECT udTablesPCH.cxx)
# The PCH compile is one of the heaviest translation units in the tree -- it
# instantiates the whole data model up front -- but add_library() does not put
# it in the 'analysis' job pool that o2physics_add_dpl_workflow() uses for the
# workflows. Left out, it runs unthrottled alongside the pooled compilations
# and defeats the per-job memory budget the pool exists to enforce.
set_property(TARGET UDTablesPCH PROPERTY JOB_POOL_COMPILE analysis)
target_link_libraries(UDTablesPCH PUBLIC O2::Framework O2Physics::AnalysisCore)
if(NOT DEFINED ENV{USE_RECC})
  target_precompile_headers(UDTablesPCH PRIVATE
    # Framework baseline — same set AnalysisPCH precompiles. A target with
    # REUSE_FROM uses ONLY this PCH, so anything left out here is parsed from
    # scratch by every consumer: without this block these workflows were
    # precompiling LESS than the default carrier they opted out of.
    <Framework/ASoA.h>
    <Framework/AnalysisDataModel.h>
    <Framework/AnalysisHelpers.h>
    <Framework/AnalysisTask.h>
    <Framework/ConfigContext.h>
    <Framework/Configurable.h>
    <Framework/DataProcessorSpec.h>
    <Framework/HistogramRegistry.h>
    <Framework/HistogramSpec.h>
    <Framework/InitContext.h>
    [["PWGUD/DataModel/UDTables.h"]]           # 79%, ~650 lines
  )
endif()
