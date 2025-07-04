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
/// \file     model.cxx
///
/// \author   Christian Sonnabend <christian.sonnabend@cern.ch>
///
/// \brief    A general-purpose class with functions for ONNX model applications
///

#include "Tools/ML/model.h"

#include <Framework/Logger.h>

#include <TSystem.h>

#include <onnxruntime_cxx_api.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace o2
{

namespace ml
{

std::string OnnxModel::printShape(const std::vector<int64_t>& v)
{
  std::stringstream ss("");
  for (size_t i = 0; i < v.size() - 1; i++)
    ss << v[i] << "x";
  ss << v[v.size() - 1];
  return ss.str();
}

bool OnnxModel::checkHyperloop(bool verbose)
{
  /// Testing hyperloop core settings
  const char* alienCores = gSystem->Getenv("ALIEN_JDL_CPUCORES");
  bool alienCoresFound = (alienCores != NULL);
  if (alienCoresFound) {
    if (verbose) {
      LOGP(info, "Hyperloop test/Grid job detected! Number of cores = {}. Setting threads anyway to 1.", alienCores);
    }
    activeThreads = 1;
    sessionOptions.SetIntraOpNumThreads(activeThreads);
  } else {
    if (verbose) {
      LOGP(info, "Not running on Hyperloop.");
    }
  }

  return alienCoresFound;
}

void OnnxModel::initModel(std::string localPath, bool enableOptimizations, int threads, uint64_t from, uint64_t until)
{

  assert(from <= until);

  LOG(info) << "--- ONNX-ML model ---";
  modelPath = localPath;
  activeThreads = threads;

  /// Running on Hyperloop
  if (!checkHyperloop(true)) {
    sessionOptions.SetIntraOpNumThreads(activeThreads);
  }

  /// Enableing optimizations
  if (enableOptimizations) {
    sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);
  }

  mEnv = std::make_shared<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "onnx-model");
  mSession = std::make_shared<Ort::Session>(*mEnv, modelPath.c_str(), sessionOptions);

  Ort::AllocatorWithDefaultOptions tmpAllocator;
  for (size_t i = 0; i < mSession->GetInputCount(); ++i) {
    mInputNames.push_back(mSession->GetInputNameAllocated(i, tmpAllocator).get());
  }
  for (size_t i = 0; i < mSession->GetInputCount(); ++i) {
    mInputShapes.emplace_back(mSession->GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape());
  }
  for (size_t i = 0; i < mSession->GetOutputCount(); ++i) {
    mOutputNames.push_back(mSession->GetOutputNameAllocated(i, tmpAllocator).get());
  }
  for (size_t i = 0; i < mSession->GetOutputCount(); ++i) {
    mOutputShapes.emplace_back(mSession->GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape());
  }
  LOG(info) << "Input Nodes:";
  for (size_t i = 0; i < mInputNames.size(); i++) {
    LOG(info) << "\t" << mInputNames[i] << " : " << printShape(mInputShapes[i]);
  }

  LOG(info) << "Output Nodes:";
  for (size_t i = 0; i < mOutputNames.size(); i++) {
    LOG(info) << "\t" << mOutputNames[i] << " : " << printShape(mOutputShapes[i]);
  }

  validFrom = from;
  validUntil = until;

  LOG(info) << "Model validity - From: " << validFrom << ", Until: " << validUntil;

  LOG(info) << "--- Model initialized! ---";
}

void OnnxModel::setActiveThreads(int threads)
{
  activeThreads = threads;
  if (!checkHyperloop(false)) {
    sessionOptions.SetIntraOpNumThreads(activeThreads);
  }
}

// ONNX protobuf model surgery
void OnnxModel::modelSurgery(std::string inpath, std::string outpath, std::unordered_map<std::string, float> specs)
{
  // Example for specs with model input of size Nx7:
  // {
  //   {"1", 1},
  //   {"2", 1},
  //   {"3", 1},
  //   {"4", 1},
  //   {"5", 1},
  //   {"6", 1},
  //   {"7", 1},
  // }
  std::ostringstream buffer;

  onnx::ModelProto model;
  std::ifstream input(inpath.c_str(), std::ios::in | std::ios::binary);

  if (!input) {
    throw std::runtime_error("Failed to open model file: " + inpath);
  }
  if (!model.ParseFromIstream(&input)) {
    throw std::runtime_error("Failed to parse ONNX model from stream.");
  }

  model.set_ir_version(onnx::IR_VERSION);
  model.set_producer_name("example_linear");

  onnx::GraphProto* graph = model.mutable_graph();
  graph->set_name("LinearGraph");
  add_concat_to_input(model, specs);

  std::cout << "=== ONNX Model Summary ===\n";
  std::cout << "Producer: " << model.producer_name() << "\n";
  std::cout << "IR version: " << model.ir_version() << "\n";
  std::cout << "Graph name: " << model.graph().name() << "\n";

  // Inputs
  {
    onnx::GraphProto const& graph = model.graph();
    std::cout << "\nInputs:\n";
    for (const auto& input : graph.input()) {
      std::cout << "  " << input.name() << " : ";
      if (input.has_type() && input.type().has_tensor_type()) {
        const auto& shape = input.type().tensor_type().shape();
        print_shape(shape);
      }
      std::cout << "\n";
    }

    // Outputs
    std::cout << "\nOutputs:\n";
    for (const auto& output : graph.output()) {
      std::cout << "  " << output.name() << " : ";
      if (output.has_type() && output.type().has_tensor_type()) {
        const auto& shape = output.type().tensor_type().shape();
        print_shape(shape);
      }
      std::cout << "\n";
    }

    // Nodes
    std::cout << "\nNodes:\n";
    for (const auto& node : graph.node()) {
      std::cout << "  [" << node.op_type() << "] ";
      for (const auto& input : node.input()) {
        std::cout << input << " ";
      }
      std::cout << "-> ";
      for (const auto& output : node.output()) {
        std::cout << output << " ";
      }
      std::cout << "\n";
    }
  }

  // Export the modified model to ./model.onnx
  std::ofstream outModel(outpath.c_str(), std::ios::out | std::ios::binary);
  if (!outModel) {
    std::cerr << "Failed to open output file for writing ONNX model.\n";
    return;
  }
  if (!model.SerializeToOstream(&outModel)) {
    std::cerr << "Failed to serialize ONNX model to file.\n";
    return;
  }
  std::cout << "Modified ONNX model exported to ./model.onnx\n";

  if (!model.SerializeToOstream(&buffer)) {
    throw std::runtime_error("Failed to serialize modified model.");
  }
}

void OnnxModel::add_concat_to_input(onnx::ModelProto& model, std::unordered_map<std::string, float> originals)
{
  auto* graph = model.mutable_graph();

  // 1. Save original input name and remove it
  std::string old_input_name = graph->input(0).name();
  graph->mutable_input()->DeleteSubrange(0, 1);

  // 2. Create N new input tensors: input_0, ..., input_N-1
  std::vector<std::string> input_names;
  int i = 0;
  for (const auto& kv : originals) {
    std::string name = "input_" + kv.first;
    input_names.push_back(name);

    auto* input = graph->add_input();
    input->set_name(name);
    auto* tensor_type = input->mutable_type()->mutable_tensor_type();
    tensor_type->set_elem_type(onnx::TensorProto_DataType_FLOAT);
    auto* shape = tensor_type->mutable_shape();
    shape->add_dim()->set_dim_param("N");
    shape->add_dim()->set_dim_value(static_cast<int64_t>(kv.second));
    ++i;
  }

  // 3. Add Concat node
  auto* concat = graph->add_node();
  concat->set_op_type("Concat");
  concat->set_name("concat_inputs");
  for (const auto& name : input_names) {
    concat->add_input(name);
  }
  concat->add_output(old_input_name); // Output replaces the original input
  auto* attr = concat->add_attribute();
  attr->set_name("axis");
  attr->set_type(onnx::AttributeProto_AttributeType_INT);
  attr->set_i(1);
}

void OnnxModel::print_shape(const onnx::TensorShapeProto& shape)
{
  std::cout << "[";
  for (int i = 0; i < shape.dim_size(); ++i) {
    if (i > 0) {
      std::cout << ", ";
    }
    if (shape.dim(i).has_dim_value()) {
      std::cout << shape.dim(i).dim_value();
    } else if (shape.dim(i).has_dim_param()) {
      std::cout << shape.dim(i).dim_param();
    } else {
      std::cout << "?";
    }
  }
  std::cout << "]";
}

} // namespace ml

} // namespace o2
