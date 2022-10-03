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

/// \file simpleApplyPidInterface
/// \brief A simple example for using PID obtained from the PID ML ONNX Interface. See https://github.com/saganatt/PID_ML_in_O2 for more detailed instructions.
///
/// \author Maja Kabus <mkabus@cern.ch>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CCDB/CcdbApi.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Tools/PIDML/pidOnnxInterface.h"

#include <string>

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

struct SimpleApplyOnnxInterface {
  PidONNXInterface pidInterface; // One instance to manage all needed ONNX models

  Configurable<LabeledArray<double>> cfgPTCuts{"pT_cuts", {pidml_pt_cuts::cuts[0], pidml_pt_cuts::nPids, pidml_pt_cuts::nCutVars, pidml_pt_cuts::pidLabels, pidml_pt_cuts::cutVarLabels}, "pT cuts for each output pid and each detector configuration"};
  Configurable<std::vector<int>> cfgPids{"pids", std::vector<int>{pidml_pt_cuts::pids_v}, "PIDs to predict"};
  Configurable<std::vector<double>> cfgCertainties{"certainties", std::vector<double>{pidml_pt_cuts::certainties_v}, "Min certainties of the models to accept given particle to be of given kind"};
  Configurable<bool> cfgAutoMode{"autoMode", true, "Use automatic model matching: default pT cuts and min certainties"};

  Configurable<std::string> cfgPathCCDB{"ccdb-path", "Users/m/mkabus/PIDML", "base path to the CCDB directory with ONNX models"};
  Configurable<std::string> cfgCCDBURL{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<bool> cfgUseCCDB{"useCCDB", true, "Whether to autofetch ML model from CCDB. If false, local file will be used."};
  Configurable<std::string> cfgPathLocal{"local-path", "/home/mkabus/PIDML/", "base path to the local directory with ONNX models"};

  o2::ccdb::CcdbApi ccdbApi;
  int currentRunNumber = -1;

  Produces<o2::aod::MlPidResults> pidMLResults;

  Filter trackFilter = requireGlobalTrackInFilter();
  // Minimum table requirements for sample model:
  // TPC signal (FullTracks), TOF signal (TOFSignal), TOF beta (pidTOFbeta), dcaXY and dcaZ (TracksDCA)
  // Filter on isGlobalTrack (TracksSelection)
  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal>>;

  void init(InitContext const&)
  {
    if (cfgUseCCDB) {
      ccdbApi.init(cfgCCDBURL);
    } else {
      pidInterface = PidONNXInterface(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, -1, cfgPids.value, cfgPTCuts.value, cfgCertainties.value, cfgAutoMode.value);
    }
  }

  void processCollisions(aod::Collisions const& collisions, BigTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
    if (cfgUseCCDB && bc.runNumber() != currentRunNumber) {
      pidInterface = PidONNXInterface(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, bc.timestamp(), cfgPids.value, cfgPTCuts.value, cfgCertainties.value, cfgAutoMode.value);
    }

    for (auto& track : tracks) {
      for (int pid : cfgPids.value) {
        bool accepted = pidInterface.applyModelBoolean(track, pid);
        LOGF(info, "collision id: %d track id: %d pid: %d accepted: %d p: %.3f; x: %.3f, y: %.3f, z: %.3f",
             track.collisionId(), track.index(), pid, accepted, track.p(), track.x(), track.y(), track.z());
        pidMLResults(track.index(), pid, accepted);
      }
    }
  }
  PROCESS_SWITCH(SimpleApplyOnnxInterface, processCollisions, "Process with collisions and bcs for CCDB", true);

  void processTracksOnly(BigTracks const& tracks)
  {
    for (auto& track : tracks) {
      for (int pid : cfgPids.value) {
        bool accepted = pidInterface.applyModelBoolean(track, pid);
        LOGF(info, "collision id: %d track id: %d pid: %d accepted: %d p: %.3f; x: %.3f, y: %.3f, z: %.3f",
             track.collisionId(), track.index(), pid, accepted, track.p(), track.x(), track.y(), track.z());
        pidMLResults(track.index(), pid, accepted);
      }
    }
  }
  PROCESS_SWITCH(SimpleApplyOnnxInterface, processTracksOnly, "Process with tracks only -- faster but no CCDB", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SimpleApplyOnnxInterface>(cfgc)};
}
