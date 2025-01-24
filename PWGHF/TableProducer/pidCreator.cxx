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

#include <string>

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TableHelper.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"

using namespace o2;
using namespace o2::framework;

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

  static constexpr float defaultNSigmaTolerance = .1f;
  static constexpr float defaultNSigma = -999.f + defaultNSigmaTolerance; // -999.f is the default value set in TPCPIDResponse.h and PIDTOF.h

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
  }

  /// Function to combine TPC and TOF NSigma
  /// \param tiny switch between full and tiny (binned) PID tables
  /// \param tpcNSigma is the (binned) NSigma separation in TPC (if tiny = true)
  /// \param tofNSigma is the (binned) NSigma separation in TOF (if tiny = true)
  /// \return combined NSigma of TPC and TOF
  template <bool tiny, typename T1>
  T1 combineNSigma(T1 tpcNSigma, T1 tofNSigma)
  {
    if constexpr (tiny) {
      tpcNSigma *= aod::pidtpc_tiny::binning::bin_width;
      tofNSigma *= aod::pidtof_tiny::binning::bin_width;
    }
    if ((tpcNSigma > defaultNSigma) && (tofNSigma > defaultNSigma)) { // TPC and TOF
      return std::sqrt(.5f * (tpcNSigma * tpcNSigma + tofNSigma * tofNSigma));
    }
    if (tpcNSigma > defaultNSigma) { // only TPC
      return std::abs(tpcNSigma);
    }
    if (tofNSigma > defaultNSigma) { // only TOF
      return std::abs(tofNSigma);
    }
    return tofNSigma; // no TPC nor TOF
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

#undef PROCESS_PID
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfPidCreator>(cfgc)};
}
