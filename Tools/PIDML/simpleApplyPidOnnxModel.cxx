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

/// \file simpleApplyPidModelTask
/// \brief A simple example for using PID obtained from the PID ML ONNX Model. See https://github.com/saganatt/PID_ML_in_O2 for more detailed instructions.
///
/// \author Maja Kabus <mkabus@cern.ch>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
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
DECLARE_SOA_INDEX_COLUMN(Track, track); //! Track index
DECLARE_SOA_COLUMN(MlPid, mlPid, int);  //! PID predicted by the model
} // namespace mlpidresult
DECLARE_SOA_TABLE(MlPidResults, "AOD", "MLPIDRESULTS", o2::soa::Index<>, mlpidresult::TrackId, mlpidresult::MlPid);
} // namespace o2::aod

struct SimpleApplyOnnxModelTask {
  PidONNXModel pidModel; // One instance per model, e.g., one per each pid to predict
  Configurable<bool> cfgUseTOF{"useTOF", true, "Use ML model with TOF signal"};
  Configurable<int> cfgPid{"pid", 211, "PID to predict"};

  // TODO: CCDB does not support ROOT files yet. Temporary solution: git repository with ML models
  // Paths given here are relative to the main models directory $MLMODELS_ROOT
  Configurable<std::string> cfgScalingParamsFile{"scaling-params", "train_208_mc_with_beta_and_sigmas_scaling_params.json", "JSON file with scaling parameters from training"};

  // Configurable<std::string> cfgModelDir{"model-dir", "http://alice-ccdb.cern.ch/Users/m/mkabus/pidml/onnx_models", "base path to the directory with ONNX models"};
  // Configurable<std::string> cfgScalingParamsFile{"scaling-params", "http://alice-ccdb.cern.ch/Users/m/mkabus/pidml/onnx_models/train_208_mc_with_beta_and_sigmas_scaling_params.json", "base path to the ccdb JSON file with scaling parameters from training"};
  //  Paths to local files
  //  Configurable<std::string> cfgModelDir{"model-dir", "/home/maja/CERN_part/CERN/PIDML/onnx_models", "base path to the directory with ONNX models"};
  //  Configurable<std::string> cfgScalingParamsFile{"scaling-params", "/home/maja/CERN_part/CERN/PIDML/onnx_models/train_208_mc_with_beta_and_sigmas_scaling_params.json", "JSON file with scaling parameters from training"};

  Produces<o2::aod::MlPidResults> pidMLResults;

  Filter trackFilter = requireGlobalTrackInFilter();
  // Minimum table requirements for sample model:
  // TPC signal (FullTracks), TOF signal (TOFSignal), TOF beta (pidTOFbeta), dcaXY and dcaZ (TracksDCA)
  // Filter on isGlobalTrack (TracksSelection)
  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal>>;

  void init(InitContext const&)
  {
    pidModel = PidONNXModel(cfgScalingParamsFile.value, cfgPid.value, cfgUseTOF.value);
  }

  void process(BigTracks const& tracks)
  {
    for (auto& track : tracks) {
      float pid = pidModel.applyModel(track);
      // pid > 0 --> track is predicted to be of this kind; pid < 0 --> rejected
      LOGF(info, "collision id: %d track id: %d pid: %.3f p: %.3f; x: %.3f, y: %.3f, z: %.3f",
           track.collisionId(), track.index(), pid, track.p(), track.x(), track.y(), track.z());
      pidMLResults(track.index(), pid);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SimpleApplyOnnxModelTask>(cfgc)};
}
