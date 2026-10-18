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

o2physics_add_library(v0SelectionGroup
    SOURCES v0SelectionGroup.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore)

o2physics_target_root_dictionary(v0SelectionGroup
    HEADERS v0SelectionGroup.h
    LINKDEF v0SelectionGroupLinkDef.h)

o2physics_add_library(collisionCutsGroup
    SOURCES collisionCutsGroup.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore)

o2physics_target_root_dictionary(collisionCutsGroup
    HEADERS collisionCutsGroup.h
    LINKDEF collisionCutsGroupLinkDef.h)

o2physics_add_library(LfDecay3bodyBuilderHelper
    SOURCES decay3bodyBuilderHelper.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore KFParticle::KFParticle)

o2physics_add_library(LfStrangenessBuilderHelper
    SOURCES strangenessBuilderHelper.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore KFParticle::KFParticle)

o2physics_add_library(LfStrangenessBuilderModule
    SOURCES strangenessBuilderModule.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore KFParticle::KFParticle)

o2physics_add_library(LfCascadeMlResponse
    SOURCES cascadeMlResponse.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::MLCore)

# Precompiled header shared across the PWGLF strangeness workflows. Percentages
# are the share of the group's workflow translation units that include the header,
# measured over the tree. LFStrangenessPIDTables.h is deliberately left out: it
# reaches only 26% of them, and it includes LFStrangenessTables.h anyway.
#
# runDataProcessing.h is deliberately NOT precompiled: it requires each workflow's
# customize() to be declared before it is included, so precompiling it silently
# unregisters every option (the binary aborts at --dump-workflow, but links fine).
#
# Skipped under recc, which caches compilations remotely instead.
add_library(LfStrangenessPCH OBJECT strangenessPCH.cxx)
# The PCH compile is one of the heaviest translation units in the tree -- it
# instantiates the whole data model up front -- but add_library() does not put
# it in the 'analysis' job pool that o2physics_add_dpl_workflow() uses for the
# workflows. Left out, it runs unthrottled alongside the pooled compilations
# and defeats the per-job memory budget the pool exists to enforce.
set_property(TARGET LfStrangenessPCH PROPERTY JOB_POOL_COMPILE analysis)
target_link_libraries(LfStrangenessPCH PUBLIC O2::Framework O2Physics::AnalysisCore)
if(NOT DEFINED ENV{USE_RECC})
  target_precompile_headers(LfStrangenessPCH PRIVATE
    # Framework baseline -- same set AnalysisPCH precompiles. A target with
    # REUSE_FROM uses ONLY this PCH, so anything left out here is parsed from
    # scratch by every consumer: without this block these workflows precompile
    # LESS than the default carrier they opted out of.
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
    [["PWGLF/DataModel/LFStrangenessTables.h"]]  # 100%, 2045 lines
  )
endif()
