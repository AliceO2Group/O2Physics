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

o2physics_add_library(PWGDQCore
               SOURCES  VarManager.cxx
                        HistogramManager.cxx
                        HistogramsLibrary.cxx
                        CutsLibrary.cxx
                        MixingLibrary.cxx
                        MCSignalLibrary.cxx
                        MixingHandler.cxx
                        AnalysisCut.cxx
                        AnalysisCompositeCut.cxx
                        MCProng.cxx
                        MCSignal.cxx
               PUBLIC_LINK_LIBRARIES O2::Framework O2::DCAFitter O2::GlobalTracking O2Physics::AnalysisCore KFParticle::KFParticle O2Physics::MLCore)

# Precompiled header shared across PWGDQ. Percentages are the share of the
# group's workflow translation units that include the header, measured over the tree.
# VarManager.h alone is ~8000 lines and is parsed by two thirds of the group.
#
# runDataProcessing.h is deliberately NOT precompiled: it requires each workflow's
# customize() to be declared before it is included, so precompiling it silently
# unregisters every option (the binary aborts at --dump-workflow, but links fine).
#
# Skipped under recc, which caches compilations remotely instead.
add_library(DQAnalysisPCH OBJECT dqAnalysisPCH.cxx)
target_link_libraries(DQAnalysisPCH PUBLIC O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore)
if(NOT DEFINED ENV{USE_RECC})
  target_precompile_headers(DQAnalysisPCH PRIVATE
    [["PWGDQ/DataModel/ReducedInfoTables.h"]]  # 65%
    [["PWGDQ/Core/VarManager.h"]]              # 65%, ~8000 lines
  )
endif()
