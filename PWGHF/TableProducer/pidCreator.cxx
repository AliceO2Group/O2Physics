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

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TableHelper.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"

using namespace o2;
using namespace o2::framework;

struct HfPidCreator {
  Produces<aod::TracksPidFullElS> tracksPidFullElS;
  Produces<aod::TracksPidTinyElS> tracksPidTinyElS;

  static constexpr float defaultNSigmaTolerance = .1f;
  static constexpr float defaultNSigma = -999.f + defaultNSigmaTolerance; // -999.f is the default value set in TPCPIDResponse.h and PIDTOF.h

  template <typename TSwitch>
  void checkTableSwitch(InitContext& initContext, const std::string& table, const TSwitch& doprocess)
  {
    auto isNeeded = isTableRequiredInWorkflow(initContext, table);
    if (isNeeded && !doprocess.value) {
      LOGF(fatal, "Table %s is needed but not requested. Enable the corresponding process function!", table);
    }
    if (!isNeeded && doprocess.value) {
      LOGF(warn, "Table %s is requested but not needed. Disable the corresponding process function!", table);
    }
  }

  void init(InitContext& initContext)
  {
    checkTableSwitch(initContext, "TracksPidFullElS", doprocessFullEl);
    checkTableSwitch(initContext, "TracksPidTinyElS", doprocessTinyEl);
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

  void processFullEl(aod::TracksPidEl const& tracks)
  {
    for (const auto& track : tracks) {
      tracksPidFullElS(combineNSigma<false>(track.tpcNSigmaEl(), track.tofNSigmaEl()));
    }
  }
  PROCESS_SWITCH(HfPidCreator, processFullEl, "Process full El ", false);

  void processTinyEl(aod::TracksPidTinyEl const& tracks)
  {
    for (const auto& track : tracks) {
      tracksPidTinyElS(combineNSigma<true>(track.tpcNSigmaStoreEl(), track.tofNSigmaStoreEl()));
    }
  }
  PROCESS_SWITCH(HfPidCreator, processTinyEl, "Process tiny El", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfPidCreator>(cfgc)};
}
