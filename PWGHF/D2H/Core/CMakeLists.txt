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

o2physics_add_library(HfDataCreationCharmReso
  SOURCES DataCreationCharmReso.cxx
  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore KFParticle::KFParticle)

# Precompiled header shared by the D2H reduced-data-model workflows. 28 of them
# include ReducedDataModel.h, which alone reaches CandidateReconstructionTables.h,
# RecoDecay.h, the PID/EventSelection/Centrality models and AnalysisDataModel.h,
# so precompiling it covers most of what the family parses.
#
# Defined in Core because it is added before TableProducer and Tasks, and a
# REUSE_FROM target must exist before the targets that name it.
#
# Skipped under recc, which caches compilations remotely instead -- the two do
# not combine, as a PCH is a local artefact of one compiler invocation.
add_library(D2HReducedPCH OBJECT d2hReducedPCH.cxx)
target_link_libraries(D2HReducedPCH PUBLIC O2::Framework O2Physics::AnalysisCore)
if(NOT DEFINED ENV{USE_RECC})
  target_precompile_headers(D2HReducedPCH PRIVATE
    [["PWGHF/D2H/DataModel/ReducedDataModel.h"]]
    [["PWGHF/DataModel/CandidateSelectionTables.h"]]
    <Framework/ASoA.h>
    <Framework/AnalysisHelpers.h>
    <Framework/AnalysisTask.h>
    <Framework/ConfigContext.h>
    <Framework/Configurable.h>
    <Framework/DataProcessorSpec.h>
    <Framework/HistogramRegistry.h>
    <Framework/HistogramSpec.h>
    <Framework/InitContext.h>
    <Framework/Logger.h>
    # NOT runDataProcessing.h, although all 28 consumers include it and
    # PWGJE's carriers do precompile it. It states the contract itself: the
    # customize() of a workflow must be declared before it is included.
    # Precompiling it makes it the first include in every consumer, so a
    # source's customize() is never seen, its options are never registered,
    # and the binary aborts at --dump-workflow with "missing option: ...".
    # The link succeeds, so this surfaces as a runtime failure with no hint
    # that a PCH caused it.

    <array>
    <cmath>
    <cstdint>
    <memory>
    <numeric>
    <string>
    <vector>
  )
endif()
