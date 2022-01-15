// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Tools/PIDML/pidONNXInferer.h"

#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// See https://github.com/saganatt/PID_ML_in_O2 for instructions

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
  PidONNXInferer pidInferer; // We cannot use copy constructor, so it needs to be a pointer
  Configurable<bool> useGPU{"useGPU", false, "Use GPU for ML inference if true"};

  // TODO: This is a quick solution that uses CCDB to store files for Hyperloop.
  // In the future, they should be on cvmfs or somewhere else, where files different from ROOT can be stored.
  // modelFile should be in ONNX format, and trainColumnsFile and scalingParamsFile in JSON format.
  // For now ONNX and JSONs are converted to ROOT TStrings
  Configurable<std::string> modelFile{"ccdb-onnx", "Users/m/mkabus/pidml/onnx_models/All/simple_model_211.root", "base path to the ccdb ONNX model file"};
  Configurable<std::string> trainColumnsFile{"train-columns", "Users/m/mkabus/pidml/onnx_models/All/columns_for_training.root", "base path to the ccdb JSON file with track properties used for training"};
  Configurable<std::string> scalingParamsFile{"scaling-params", "Users/m/mkabus/pidml/onnx_models/train_208_mc_with_beta_and_sigmas_scaling_params.root", "base path to the ccdb JSON file with scaling parameters from training"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<long> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  Produces<o2::aod::MlPidResults> pidMLResults;

  Filter trackFilter = aod::track::isGlobalTrack == (uint8_t) true;
  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksExtended, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal>>;

  void init(InitContext const&)
  {
    pidInferer = PidONNXInferer(modelFile.value, trainColumnsFile.value, scalingParamsFile.value, url.value, nolaterthan.value, useGPU.value);
  }

  void process(BigTracks const& tracks)
  {
    for (auto& track : tracks) {
      float pid = pidInferer.applyModel(track);
      LOGF(info, "collision id: %d track id: %d pid: %d eta: %.3f; p: %.3f; x: %.3f, y: %.3f, z: %.3f",
           track.collisionId(), track.index(), pid, track.eta(), track.p(), track.x(), track.y(), track.z());
      pidMLResults(track.index(), pid);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SimpleApplyOnnxModelTask>(cfgc)};
}
