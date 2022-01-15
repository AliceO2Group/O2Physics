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
  std::string onnxFile = (std::string)onnxFileConf;

  std::vector<std::string> input_names;
  std::vector<std::vector<int64_t>> input_shapes;
  std::vector<std::string> output_names;
  std::vector<std::vector<int64_t>> output_shapes;

  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "example-model-explorer"};
  Ort::SessionOptions session_options;
  Ort::Experimental::Session session{env, onnxFile, session_options};

  HistogramRegistry onnxResults{"onnxResults", {{"results", "results", {HistType::kTH1F, {{2, 0., 1.}}}}}};

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
    input_names = session.GetInputNames();
    input_shapes = session.GetInputShapes();
    output_names = session.GetOutputNames();
    output_shapes = session.GetOutputShapes();

    // print name/shape of inputs
    LOG(info) << "Input Node Name/Shape (" << input_names.size() << "):";
    for (size_t i = 0; i < input_names.size(); i++) {
      LOG(info) << "\t" << input_names[i] << " : " << print_shape(input_shapes[i]);
    }

    // print name/shape of outputs
    LOG(info) << "Output Node Name/Shape (" << output_names.size() << "):";
    for (size_t i = 0; i < output_names.size(); i++) {
      LOG(info) << "\t" << output_names[i] << " : " << print_shape(output_shapes[i]);
    }

    // Assume model has 1 input node and 1 output node.
    assert(input_names.size() == 1 && output_names.size() == 1);
  }

  void process(aod::Tracks const& tracks)
  {
    auto input_shape = input_shapes[0];
    for (auto& track : tracks) {
      LOGF(info, "collision id: %d; eta: %.3f; p: %.3f; x: %.3f, y: %.3f, z: %.3f",
           track.collisionId(), track.eta(), track.p(), track.x(), track.y(), track.z());

      std::vector<float> input_tensor_values{track.eta(), track.p(), track.x(), track.y(), track.z()};
      std::vector<Ort::Value> input_tensors;
      input_tensors.push_back(Ort::Experimental::Value::CreateTensor<float>(input_tensor_values.data(), input_tensor_values.size(), input_shape));

      // double-check the dimensions of the input tensor
      assert(input_tensors[0].IsTensor()
               input_tensors[0]
                 .GetTensorTypeAndShapeInfo()
                 .GetShape() == input_shape);
      LOG(info) << "input_tensor shape: " << print_shape(input_tensors[0].GetTensorTypeAndShapeInfo().GetShape());

      // pass data through model
      LOG(info) << "Running model...";
      try {
        auto output_tensors = session.Run(input_names, input_tensors, output_names);
        LOG(info) << "done";
        //LOG(info) << "Number of output tensors: " << output_tensors.size();

        // double-check the dimensions of the output tensors
        // NOTE: the number of output tensors is equal to the number of output nodes specifed in the Run() call
        assert(output_tensors.size() == output_names.size() &&
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
