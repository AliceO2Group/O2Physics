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
/// \file   pidTPCFull.cxx
///
/// \author Christian Sonnabend christian.sonnabend@cern.ch
///
/// \brief  Showcase application of an ONNX model in O2Physics
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Logger.h"
#include "CCDB/CcdbApi.h"
#include "Tools/ML/model.h"
#include "TRandom3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::ml;

struct applyModel {

  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<bool> useLocalFile{"useLocalFile", true, "If true: Uses network from modelPathLocally; else: Load from CCDB to modelPathLocally"};
  Configurable<std::string> modelPathLocally{"modelPathLocally", "network.onnx", "(std::string) Path to the local .onnx file"};
  Configurable<std::string> modelPathCCDB{"modelPathCCDB", "Analysis/PID/TPC/ML", "Model-path on CCDB"};
  Configurable<bool> enableOptimizations{"enableOptimizations", true, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};

  OnnxModel network;
  std::map<std::string, std::string> metadata;
  o2::ccdb::CcdbApi ccdbApi;

  void init(InitContext const&)
  {
    if (useLocalFile) {
      network.initModel(modelPathLocally.value, enableOptimizations.value);
    } else {
      ccdbApi.init(url);
      bool retrieveSuccess = ccdbApi.retrieveBlob(modelPathCCDB.value, ".", metadata, 1, false, modelPathLocally.value); // see the ccdb api header "O2/CCDB/include/CCDB/CcdbApi.h"; Fetching an arbitrary network for showcasing now
      if (retrieveSuccess) {
        network.initModel(modelPathLocally.value, enableOptimizations.value); // initializes the model, prints out shape of the input and output already
      }
    }
  }

  void process(aod::Tracks const&)
  {
    std::vector<float> modelInput = std::vector<float>(network.getNumInputNodes() * 10, 0.); // Use an input of size 100 * (number of model inputs)
    TRandom3* fRndm = new TRandom3(0);
    for (int i = 0; i < network.getNumInputNodes() * 10; i++) {
      modelInput[i] = fRndm->Rndm();
    }
    float* modelOutput = network.evalModel(modelInput); // evaluate
    for (int i = 0; i < 10; i++) {
      LOG(info) << "Input: [" << modelInput[i * 6] << "," << modelInput[i * 6 + 1] << "," << modelInput[i * 6 + 2] << "," << modelInput[i * 6 + 3] << "," << modelInput[i * 6 + 4] << "," << modelInput[i * 6 + 5] << "] -> Output: " << modelOutput[i];
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<applyModel>(cfgc),
  };
}
