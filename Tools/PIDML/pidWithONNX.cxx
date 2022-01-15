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
#include "Framework/AnalysisDataModel.h"

#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#include <string>

using namespace o2;
using namespace o2::framework;

// See https://github.com/saganatt/PID_ML_in_O2 for instructions

struct ApplyOnnxModelTask {
  Configurable<std::string> onnxFileConf{"ONNX file", "/home/maja/CERN_part/CERN/PID_ML_in_O2/models/Simple_example.onnx", "ONNX file"};

  std::vector<std::string> mInputNames;
  std::vector<std::vector<int64_t>> mInputShapes;
  std::vector<std::string> mOutputNames;
  std::vector<std::vector<int64_t>> mOutputShapes;

  Ort::Env mEnv;
  Ort::SessionOptions mSessionOptions;
  Ort::Experimental::Session mSession;
  // TODO: What is this allocator for? How is it initialized?
  //Ort::AllocatorWithDefaultOptions mAllocator;
  //Ort::MemoryInfo mMemoryInfo;

  // pretty prints a shape dimension vector
  std::string print_shape(const std::vector<int64_t>& v)
  {
    std::stringstream ss("");
    for (size_t i = 0; i < v.size() - 1; i++)
      ss << v[i] << "x";
    ss << v[v.size() - 1];
    return ss.str();
  }

  void init(InitContext const&)
  {
    mEnv = Ort::Env{ORT_LOGGING_LEVEL_WARNING, "pid-ml-model-inference"};
    mSession = Ort::Experimental::Session{mEnv, (std::string)onnxFileConf, mSessionOptions};
    // TODO: Ask what these are needed for
    //OrtAllocatorType allocatorType;
    //OrtMemType memoryType;
    //mMemoryInfo(Ort::MemoryInfo::CreateCpu(allocatorType, memoryType)) {}

    // TODO: In their code is more manual. Why?
    mInputNames = session.GetInputNames();
    mInputShapes = session.GetInputShapes();
    mOutputNames = session.GetOutputNames();
    mOutputShapes = session.GetOutputShapes();

    // print name/shape of inputs
    LOG(info) << "Input Node Name/Shape (" << mInputNames.size() << "):";
    for (size_t i = 0; i < mInputNames.size(); i++) {
      LOG(info) << "\t" << mInputNames[i] << " : " << print_shape(mInputShapes[i]);
    }

    // print name/shape of outputs
    LOG(info) << "Output Node Name/Shape (" << mOutputNames.size() << "):";
    for (size_t i = 0; i < mOutputNames.size(); i++) {
      LOG(info) << "\t" << mOutputNames[i] << " : " << print_shape(mOutputShapes[i]);
    }

    // Assume model has 1 input node and 1 output node.
    //assert(mInputNames.size() == 1 && mOutputNames.size() == 1);
  }

  void process(aod::Tracks const& tracks)
  {
    auto input_shape = mInputShapes[0];
    for (auto& track : tracks) {
      LOGF(info, "collision id: %d; eta: %.3f; p: %.3f; x: %.3f, y: %.3f, z: %.3f",
           track.collisionId(), track.eta(), track.p(), track.x(), track.y(), track.z());

      std::vector<float> inputTensorValues{track.eta(), track.p(), track.x(), track.y(), track.z()};
      std::vector<Ort::Value> inputTensors;
      inputTensors.push_back(Ort::Experimental::Value::CreateTensor<float>(inputTensorValues.data(), inputTensorValues.size(), input_shape));

      // double-check the dimensions of the input tensor
      assert(inputTensors[0].IsTensor()
               inputTensors[0]
                 .GetTensorTypeAndShapeInfo()
                 .GetShape() == input_shape);
      LOG(info) << "input_tensor shape: " << print_shape(inputTensors[0].GetTensorTypeAndShapeInfo().GetShape());

      // pass data through model
      LOG(info) << "Running model...";
      try {
        auto output_tensors = session.Run(mInputNames, inputTensors, mOutputNames);
        LOG(info) << "done";
        //LOG(info) << "Number of output tensors: " << output_tensors.size();

        // double-check the dimensions of the output tensors
        // NOTE: the number of output tensors is equal to the number of output nodes specifed in the Run() call
        assert(output_tensors.size() == mOutputNames.size() &&
               output_tensors[0].IsTensor());
        LOG(info) << "output_tensor_shape: " << print_shape(output_tensors[0].GetTensorTypeAndShapeInfo().GetShape());

        for (auto& output : output_tensors) {
          const float* output_value = output.GetTensorData<float>();
          LOG(info) << "output: " << *output_value;
          onnxResults.get<TH1>(HIST("results"))->Fill(*output_value);
        }

      } catch (const Ort::Exception& exception) {
        LOG(error) << "Error running model inference: " << exception.what();
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ApplyOnnxModelTask>(cfgc)};
}
