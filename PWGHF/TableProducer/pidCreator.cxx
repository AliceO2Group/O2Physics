// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file pidCreator.cxx
/// \brief Workflow to produce tables with TPC+TOF combined n sigma
///
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#include "PWGHF/DataModel/AliasTables.h" // IWYU pragma: keep
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsPid.h"

#include "Common/Core/TableHelper.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::pid_tpc_tof_utils;

struct HfPidCreator {
  Produces<aod::PidTpcTofFullEl> trackPidFullEl;
  Produces<aod::PidTpcTofTinyEl> trackPidTinyEl;
  Produces<aod::PidTpcTofFullMu> trackPidFullMu;
  Produces<aod::PidTpcTofTinyMu> trackPidTinyMu;
  Produces<aod::PidTpcTofFullPi> trackPidFullPi;
  Produces<aod::PidTpcTofTinyPi> trackPidTinyPi;
  Produces<aod::PidTpcTofFullKa> trackPidFullKa;
  Produces<aod::PidTpcTofTinyKa> trackPidTinyKa;
  Produces<aod::PidTpcTofFullPr> trackPidFullPr;
  Produces<aod::PidTpcTofTinyPr> trackPidTinyPr;
  Produces<aod::PidTpcTofFullDe> trackPidFullDe;
  Produces<aod::PidTpcTofTinyDe> trackPidTinyDe;
  Produces<aod::PidTpcTofFullTr> trackPidFullTr;
  Produces<aod::PidTpcTofTinyTr> trackPidTinyTr;
  Produces<aod::PidTpcTofFullHe> trackPidFullHe;
  Produces<aod::PidTpcTofTinyHe> trackPidTinyHe;
  Produces<aod::PidTpcTofFullAl> trackPidFullAl;
  Produces<aod::PidTpcTofTinyAl> trackPidTinyAl;
  /// Function to check whether the process function flag matches the need for filling the table
  /// \param initContext  workflow context (argument of the init function)
  /// \param table  name of the table
  /// \param doprocess  process function flag
  template <typename TFlag>
  void checkTableSwitch(InitContext& initContext, const std::string& table, TFlag& doprocess)
  {
    auto isNeeded = isTableRequiredInWorkflow(initContext, table);
    if (isNeeded && !doprocess.value) {
      LOGF(fatal, "Table %s is needed but not requested. Enable the corresponding process function!", table);
    }
    if (!isNeeded && doprocess.value) {
      LOGF(warn, "Table %s is requested but not needed. Disabling the corresponding process function!", table);
      // NOTE: This does not remove the input table subscription from the context! The input table is still considered consumed.
      doprocess.value = false;
    }
  }

  void init(InitContext& initContext)
  {
    // Check whether the right process functions are enabled.
    checkTableSwitch(initContext, "PidTpcTofFullEl", doprocessFullEl);
    checkTableSwitch(initContext, "PidTpcTofTinyEl", doprocessTinyEl);
    checkTableSwitch(initContext, "PidTpcTofFullMu", doprocessFullMu);
    checkTableSwitch(initContext, "PidTpcTofTinyMu", doprocessTinyMu);
    checkTableSwitch(initContext, "PidTpcTofFullPi", doprocessFullPi);
    checkTableSwitch(initContext, "PidTpcTofTinyPi", doprocessTinyPi);
    checkTableSwitch(initContext, "PidTpcTofFullKa", doprocessFullKa);
    checkTableSwitch(initContext, "PidTpcTofTinyKa", doprocessTinyKa);
    checkTableSwitch(initContext, "PidTpcTofFullPr", doprocessFullPr);
    checkTableSwitch(initContext, "PidTpcTofTinyPr", doprocessTinyPr);
    checkTableSwitch(initContext, "PidTpcTofFullDe", doprocessFullDe);
    checkTableSwitch(initContext, "PidTpcTofTinyDe", doprocessTinyDe);
    checkTableSwitch(initContext, "PidTpcTofFullTr", doprocessFullTr);
    checkTableSwitch(initContext, "PidTpcTofTinyTr", doprocessTinyTr);
    checkTableSwitch(initContext, "PidTpcTofFullHe", doprocessFullHe);
    checkTableSwitch(initContext, "PidTpcTofTinyHe", doprocessTinyHe);
    checkTableSwitch(initContext, "PidTpcTofFullAl", doprocessFullAl);
    checkTableSwitch(initContext, "PidTpcTofTinyAl", doprocessTinyAl);
  }

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(HfPidCreator, processDummy, "Process nothing", true);

// Macro for declaring process functions per species
#define PROCESS_PID(_Species_)                                                                                            \
  void processFull##_Species_(aod::TracksPid##_Species_ const& tracks)                                                    \
  {                                                                                                                       \
    for (const auto& track : tracks) {                                                                                    \
      trackPidFull##_Species_(combineNSigma<false>(track.tpcNSigma##_Species_(), track.tofNSigma##_Species_()));          \
    }                                                                                                                     \
  }                                                                                                                       \
  PROCESS_SWITCH(HfPidCreator, processFull##_Species_, "Process full " #_Species_, false);                                \
                                                                                                                          \
  void processTiny##_Species_(aod::TracksPidTiny##_Species_ const& tracks)                                                \
  {                                                                                                                       \
    for (const auto& track : tracks) {                                                                                    \
      trackPidTiny##_Species_(combineNSigma<true>(track.tpcNSigmaStore##_Species_(), track.tofNSigmaStore##_Species_())); \
    }                                                                                                                     \
  }                                                                                                                       \
  PROCESS_SWITCH(HfPidCreator, processTiny##_Species_, "Process tiny " #_Species_, false);

  // Declare process functions for all species.
  PROCESS_PID(El)
  PROCESS_PID(Mu)
  PROCESS_PID(Pi)
  PROCESS_PID(Ka)
  PROCESS_PID(Pr)
  PROCESS_PID(De)
  PROCESS_PID(Tr)
  PROCESS_PID(He)
  PROCESS_PID(Al)

#undef PROCESS_PID
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfPidCreator>(cfgc)};
}
