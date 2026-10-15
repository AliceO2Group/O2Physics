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

o2physics_add_library(HfMlCore
  SOURCES HfMlResponse.cxx
  PUBLIC_LINK_LIBRARIES O2Physics::MLCore)

# Precompiled header shared across PWGHF. Percentages are the share of the
# group's workflow translation units that include the header, measured over the tree.
# D2H already has its own carrier (D2HReducedPCH); those targets keep it.
#
# runDataProcessing.h is deliberately NOT precompiled: it requires each workflow's
# customize() to be declared before it is included, so precompiling it silently
# unregisters every option (the binary aborts at --dump-workflow, but links fine).
#
# Skipped under recc, which caches compilations remotely instead.
add_library(HfCandidatesPCH OBJECT hfCandidatesPCH.cxx)
target_link_libraries(HfCandidatesPCH PUBLIC O2::Framework O2Physics::AnalysisCore)
if(NOT DEFINED ENV{USE_RECC})
  target_precompile_headers(HfCandidatesPCH PRIVATE
    [["PWGHF/DataModel/CandidateReconstructionTables.h"]]  # 78%, 2479 lines
    [["PWGHF/DataModel/CandidateSelectionTables.h"]]       # 72%
    [["Common/Core/RecoDecay.h"]]                          # 65%, 1523 lines
  )
endif()
