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

///
/// \file   pidTOFbeta.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to produce TOF beta tables
///         QA histograms for the TOF PID can be produced by adding `--add-qa 1` to the workflow
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "pidTOFBase.h"
#include "DPG/Tasks/qaPIDTOF.h"

using namespace o2;
using namespace o2::pid;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Produce TOF PID QA histograms"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

struct tofPidBeta {
  Produces<aod::pidTOFbeta> tablePIDBeta;
  Configurable<float> expreso{"tof-expreso", 80, "Expected resolution for the computation of the expected beta"};

  void init(o2::framework::InitContext&)
  {
    if (doprocessNoEvTime == true && doprocessEvTime == true) {
      LOGF(fatal, "Cannot enable processNoEvTime and processEvTime at the same time. Please choose one.");
    }
    responseBeta.mExpectedResolution = expreso.value;
    responseBetaEvTime.mExpectedResolution = expreso.value;
  }

  using TrksEvTime = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime>;
  tof::Beta<TrksEvTime::iterator> responseBetaEvTime;
  template <o2::track::PID::ID pid>
  using ResponseImplementationEvTime = o2::pid::tof::ExpTimes<TrksEvTime::iterator, pid>;
  void processEvTime(TrksEvTime const& tracks, aod::Collisions const&)
  {
    tablePIDBeta.reserve(tracks.size());
    for (auto const& trk : tracks) {
      tablePIDBeta(responseBetaEvTime.GetBeta(trk, trk.tofEvTime()),
                   responseBetaEvTime.GetExpectedSigma(trk),
                   responseBetaEvTime.GetExpectedSignal<o2::track::PID::Electron>(trk),
                   responseBetaEvTime.GetExpectedSigma(trk),
                   responseBetaEvTime.GetSeparation<o2::track::PID::Electron>(trk));
    }
  }
  PROCESS_SWITCH(tofPidBeta, processEvTime, "Produce TOF beta with TOF event time", false);

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal>;
  tof::Beta<Trks::iterator> responseBeta;
  template <o2::track::PID::ID pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<Trks::iterator, pid>;
  void processNoEvTime(Trks const& tracks, aod::Collisions const&)
  {
    tablePIDBeta.reserve(tracks.size());
    for (auto const& trk : tracks) {
      tablePIDBeta(responseBeta.GetBeta(trk),
                   responseBeta.GetExpectedSigma(trk),
                   responseBeta.GetExpectedSignal<o2::track::PID::Electron>(trk),
                   responseBeta.GetExpectedSigma(trk),
                   responseBeta.GetSeparation<o2::track::PID::Electron>(trk));
    }
  }
  PROCESS_SWITCH(tofPidBeta, processNoEvTime, "Produce TOF beta without TOF event time, standard for Run 2", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<tofPidBeta>(cfgc)};
  if (cfgc.options().get<int>("add-qa")) {
    workflow.push_back(adaptAnalysisTask<tofPidBetaQa>(cfgc));
  }
  return workflow;
}
