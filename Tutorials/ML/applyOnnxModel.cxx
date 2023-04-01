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

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Logger.h"
#include "CCDB/CcdbApi.h"
#include "Tools/ML/model.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::ml;

struct applyModel {

  Configurable<std::string> modelPathLocally{"modelPath", "test_net.onnx", "(std::string) Path to the local .onnx file"};
  Configurable<bool> enableOptimizations{"enableOptimizations", true, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};

  OnnxModel network;

  void init(InitContext const&)
  {
    network.initModel(modelPathLocally.value, enableOptimizations.value);
  }

  void process()
  {
    // This is our input data from the tutorial
    std::vector<float> modelInput = std::vector<float>{0.5, 3.1415926536, 0., 0.001, -3.};

    // Here we evaluate the model
    float* modelOutput = network.evalModel(modelInput);

    // And now we print the output
    for (int i = 0; i < 5; i++) {
      LOG(info) << "Input: " << modelInput[i] << ", Output: " << modelOutput[i];
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<applyModel>(cfgc),
  };
}
