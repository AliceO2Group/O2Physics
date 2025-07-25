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

/// \file simpleApplyPidOnnxInterface.cxx
/// \brief A simple example for using PID obtained from the PID ML ONNX Interface. See README.md for more detailed instructions.
///
/// \author Maja Kabus <mkabus@cern.ch>

#include "Tools/PIDML/pidOnnxInterface.h"
//
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/CcdbApi.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace mlpidresult
{
DECLARE_SOA_INDEX_COLUMN(Track, track);       //! Track index
DECLARE_SOA_COLUMN(Pid, pid, int);            //! Pid to be tested by the model
DECLARE_SOA_COLUMN(Accepted, accepted, bool); //! Whether the model accepted particle to be of given kind
} // namespace mlpidresult
DECLARE_SOA_TABLE(MlPidResults, "AOD", "MLPIDRESULTS", o2::soa::Index<>, mlpidresult::TrackId, mlpidresult::Pid, mlpidresult::Accepted);
} // namespace o2::aod

struct SimpleApplyPidOnnxInterface {
  Configurable<LabeledArray<double>> ptCuts{"ptCuts", {pidml_pt_cuts::Cuts[0], pidml_pt_cuts::NPids, pidml_pt_cuts::NCutVars, pidml_pt_cuts::pidLabels, pidml_pt_cuts::cutVarLabels}, "pT cuts for each output pid and each detector configuration"};
  Configurable<std::vector<int>> pdgPids{"pdgPids", std::vector<int>{pidml_pt_cuts::pidsV}, "PIDs to predict"};
  Configurable<std::vector<double>> mlIdentCertaintyThresholds{"mlIdentCertaintyThresholds", std::vector<double>{pidml_pt_cuts::certaintiesV}, "Min certainties of the models to accept given particle to be of given kind"};
  Configurable<bool> autoMode{"autoMode", true, "Use automatic model matching: default pT cuts and min certainties"};

  Configurable<std::string> ccdbPath{"ccdbPath", "Users/m/mkabus/PIDML", "base path to the CCDB directory with ONNX models"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<bool> useCcdb{"useCcdb", true, "Whether to autofetch ML model from CCDB. If false, local file will be used."};
  Configurable<std::string> localPath{"localPath", "/home/mkabus/PIDML/", "base path to the local directory with ONNX models"};

  Configurable<bool> useFixedTimestamp{"useFixedTimestamp", false, "Whether to use fixed timestamp from configurable instead of timestamp calculated from the data"};
  Configurable<uint64_t> fixedTimestamp{"fixedTimestamp", 1524176895000, "Hardcoded timestamp for tests"};

  o2::ccdb::CcdbApi ccdbApi;
  int currentRunNumber = -1;

  Produces<o2::aod::MlPidResults> pidMLResults;

  Filter trackFilter = requireGlobalTrackInFilter();
  // Minimum table requirements for sample model:
  // TPC signal (FullTracks), TOF signal (TOFSignal), TOF beta (pidTOFbeta), dcaXY and dcaZ (TracksDCA)
  // Filter on isGlobalTrack (TracksSelection)
  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal>>;

  PidONNXInterface<BigTracks> pidInterface; // One instance to manage all needed ONNX models

  void init(InitContext const&)
  {
    if (useCcdb) {
      ccdbApi.init(ccdbUrl);
    } else {
      pidInterface = PidONNXInterface<BigTracks>(localPath.value, ccdbPath.value, useCcdb.value, ccdbApi, -1, pdgPids.value, ptCuts.value, mlIdentCertaintyThresholds.value, autoMode.value);
    }
  }

  void processCollisions(aod::Collisions const& collisions, BigTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
    if (useCcdb && bc.runNumber() != currentRunNumber) {
      uint64_t timestamp = useFixedTimestamp ? fixedTimestamp.value : bc.timestamp();
      pidInterface = PidONNXInterface<BigTracks>(localPath.value, ccdbPath.value, useCcdb.value, ccdbApi, timestamp, pdgPids.value, ptCuts.value, mlIdentCertaintyThresholds.value, autoMode.value);
    }

    for (const auto& track : tracks) {
      for (const int& pid : pdgPids.value) {
        bool accepted = pidInterface.applyModelBoolean(track, pid);
        LOGF(info, "collision id: %d track id: %d pid: %d accepted: %d p: %.3f; x: %.3f, y: %.3f, z: %.3f",
             track.collisionId(), track.index(), pid, accepted, track.p(), track.x(), track.y(), track.z());
        pidMLResults(track.index(), pid, accepted);
      }
    }
  }
  PROCESS_SWITCH(SimpleApplyPidOnnxInterface, processCollisions, "Process with collisions and bcs for CCDB", true);

  void processTracksOnly(BigTracks const& tracks)
  {
    for (const auto& track : tracks) {
      for (const int& pid : pdgPids.value) {
        bool accepted = pidInterface.applyModelBoolean(track, pid);
        LOGF(info, "collision id: %d track id: %d pid: %d accepted: %d p: %.3f; x: %.3f, y: %.3f, z: %.3f",
             track.collisionId(), track.index(), pid, accepted, track.p(), track.x(), track.y(), track.z());
        pidMLResults(track.index(), pid, accepted);
      }
    }
  }
  PROCESS_SWITCH(SimpleApplyPidOnnxInterface, processTracksOnly, "Process with tracks only -- faster but no CCDB", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SimpleApplyPidOnnxInterface>(cfgc)};
}
