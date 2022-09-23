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

/// \file simpleApplyPidModel
/// \brief A simple example for using PID obtained from the PID ML ONNX Model. See https://github.com/saganatt/PID_ML_in_O2 for more detailed instructions.
///
/// \author Maja Kabus <mkabus@cern.ch>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CCDB/CcdbApi.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Tools/PIDML/pidOnnxModel.h"

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

struct SimpleApplyOnnxModel {
  PidONNXModel pidModel; // One instance per model, e.g., one per each pid to predict
  Configurable<uint32_t> cfgDetector{"detector", kTPCTOFTRD, "What detectors to use: 0: TPC only, 1: TPC + TOF, 2: TPC + TOF + TRD"};
  Configurable<int> cfgPid{"pid", 211, "PID to predict"};
  Configurable<float> cfgCertainty{"certainty", 0.5f, "Min certainty of the model to accept given particle to be of given kind"};

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
      pidModel = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, -1, cfgPid.value, static_cast<PidMLDetector>(cfgDetector.value), cfgCertainty.value);
    }
  }

  void processCollisions(aod::Collisions const& collisions, BigTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
    if (cfgUseCCDB && bc.runNumber() != currentRunNumber) {
      pidModel = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, bc.timestamp(), cfgPid.value, static_cast<PidMLDetector>(cfgDetector.value), cfgCertainty.value);
    }

    for (auto& track : tracks) {
      bool accepted = pidModel.applyModelBoolean(track);
      LOGF(info, "collision id: %d track id: %d accepted: %d p: %.3f; x: %.3f, y: %.3f, z: %.3f",
           track.collisionId(), track.index(), accepted, track.p(), track.x(), track.y(), track.z());
      pidMLResults(track.index(), cfgPid.value, accepted);
    }
  }
  PROCESS_SWITCH(SimpleApplyOnnxModel, processCollisions, "Process with collisions and bcs for CCDB", true);

  void processTracksOnly(BigTracks const& tracks)
  {
    for (auto& track : tracks) {
      bool accepted = pidModel.applyModelBoolean(track);
      LOGF(info, "collision id: %d track id: %d accepted: %d p: %.3f; x: %.3f, y: %.3f, z: %.3f",
           track.collisionId(), track.index(), accepted, track.p(), track.x(), track.y(), track.z());
      pidMLResults(track.index(), cfgPid.value, accepted);
    }
  }
  PROCESS_SWITCH(SimpleApplyOnnxModel, processTracksOnly, "Process with tracks only -- faster but no CCDB", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SimpleApplyOnnxModel>(cfgc)};
}
