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
/// \file     pidTPCML.h
///
/// \author   Christian Sonnabend <christian.sonnabend@cern.ch>
/// \author   Nicolò Jacazio <nicolo.jacazio@cern.ch>
///
/// \brief    A class for loading an ONNX neural network and evaluating it for the TPC PID response
///

// C++ and system includes
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#include <vector>

// ROOT includes
#include "TSystem.h"

// O2 includes
#include "Framework/Logger.h"
#include "Common/TableProducer/PID/pidTPCML.h"
#include "ReconstructionDataFormats/PID.h"

namespace o2::pid::tpc
{

std::string Network::printShape(const std::vector<int64_t>& v)
{
  std::stringstream ss("");
  for (size_t i = 0; i < v.size() - 1; i++)
    ss << v[i] << "x";
  ss << v[v.size() - 1];
  return ss.str();
}

Network::Network(std::string path,
                 bool enableOptimization = true,
                 int numThreads = 0)
{

  /*
  Constructor: Creating a class instance from a file and enabling optimizations with the boolean option.
  - Input:
    -- pathLocally:         std::string   ; Path to the model file;
    -- loadFromAlien:       bool          ; Download network from AliEn directory (true) or use local file (false)
    -- pathAlien:           std::string   ; if loadFromAlien is true, then the network will be downloaded from pathAlien
    -- enableOptimization:  bool          ; enabling optimizations for the loaded model in the session options;
  */

  LOG(info) << "--- Neural Network for the TPC PID response correction ---";

  mEnv = std::make_shared<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "pid-neural-network");
  if (enableOptimization) {
    sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);
  }
  if (numThreads > 0) {
    // Setting the number of threads, potentially try: unsigned int nThreads = std::thread::hardware_concurrency();
    sessionOptions.SetIntraOpNumThreads(numThreads);
  }

  mSession.reset(new Ort::Experimental::Session{*mEnv, path, sessionOptions});

  mInputNames = mSession->GetInputNames();
  mInputShapes = mSession->GetInputShapes();
  mOutputNames = mSession->GetOutputNames();
  mOutputShapes = mSession->GetOutputShapes();

  LOG(info) << "Input Nodes:";
  for (size_t i = 0; i < mInputNames.size(); i++) {
    LOG(info) << "\t" << mInputNames[i] << " : " << printShape(mInputShapes[i]);
  }

  LOG(info) << "Output Nodes:";
  for (size_t i = 0; i < mOutputNames.size(); i++) {
    LOG(info) << "\t" << mOutputNames[i] << " : " << printShape(mOutputShapes[i]);
  }

  LOG(info) << "--- Network initialized! ---";

} // Network::Network(std::string, bool)

Network::Network(std::string path,
                 uint64_t start,
                 uint64_t end,
                 bool enableOptimization = true,
                 int numThreads = 0)
{

  /*
  Constructor: Creating a class instance from a file and enabling optimizations with the boolean option.
  - Input:
    -- path:                std::string   ; Local path to the model file;
    -- start:               uint64_t ; Timestamp validity of model (start)
    -- pathAlien:           uint64_t ; Timestamp validity of model (end)
    -- enableOptimization:  bool          ; enabling optimizations for the loaded model in the session options;
  */

  LOG(info) << "--- Neural Network for the TPC PID response correction ---";

  mEnv = std::make_shared<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "pid-neural-network");
  if (enableOptimization) {
    sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);
  }
  if (numThreads > 0) {
    // Setting the number of threads, potentially try: unsigned int nThreads = std::thread::hardware_concurrency();
    sessionOptions.SetIntraOpNumThreads(numThreads);
  }

  mSession.reset(new Ort::Experimental::Session{*mEnv, path, sessionOptions});

  mInputNames = mSession->GetInputNames();
  mInputShapes = mSession->GetInputShapes();
  mOutputNames = mSession->GetOutputNames();
  mOutputShapes = mSession->GetOutputShapes();

  LOG(info) << "Input Nodes:";
  for (size_t i = 0; i < mInputNames.size(); i++) {
    LOG(info) << "\t" << mInputNames[i] << " : " << printShape(mInputShapes[i]);
  }

  LOG(info) << "Output Nodes:";
  for (size_t i = 0; i < mOutputNames.size(); i++) {
    LOG(info) << "\t" << mOutputNames[i] << " : " << printShape(mOutputShapes[i]);
  }

  valid_from = start;
  valid_until = end;
  LOG(info) << "Range of validity: Valid-From: " << valid_from << " Valid-Until: " << valid_until;

  LOG(info) << "--- Network initialized! ---";

} // Network::Network(std::string, uint64_t, uint64_t, bool)

Network& Network::operator=(Network& inst)
{

  /*
  Operator: Setting instances of one class to the instances of an input class
  - Input:
    -- inst:    const Network&  ; An instance of the network class;
  - Output:
    -- *this:   Network&        ; An instance with the private properties of the "inst" (input) instance;
  */

  mEnv = inst.mEnv;
  mSession = inst.mSession;
  // sessionOptions = inst.sessionOptions; // Comment (Christian): Somehow the session options throw an error when trying to copy
  mInputNames = inst.mInputNames;
  mInputShapes = inst.mInputShapes;
  mOutputNames = inst.mOutputNames;
  mOutputShapes = inst.mOutputShapes;

  valid_from = inst.valid_from;
  valid_until = inst.valid_until;

  LOG(debug) << "Network copied!";

  return *this;

} // Network& Network::operator=(const Network &)

std::vector<Ort::Value> Network::createTensor(std::array<float, 6> input) const
{

  /*
  Function: Creating a std::vector<Ort::Value> from a std::vector<float>
  - Input:
    -- input:        std::vector<float>       ; The vector from which the tensor should be created;
  - Output:
    -- inputValues:  std::vector<Ort::Value>  ; An ONNX tensor;
  */

  int64_t size = input.size();
  std::vector<int64_t> input_shape{size / mInputShapes[0][1], mInputShapes[0][1]};
  std::vector<Ort::Value> inputTensors;
  inputTensors.emplace_back(Ort::Experimental::Value::CreateTensor<float>(input.data(), size, input_shape));

  return inputTensors;

} // function std::vector<Ort::Value> Network::createTensor(std::vector<float>)

float* Network::evalNetwork(std::vector<Ort::Value> input)
{

  /*
  Function: Evaluating the network for a std::vector<Ort::Value>
  - Input:
    -- input:         std::vector<Ort::Value>     ; The tensor which should be evaluated by the network;
  - Output:
    -- output_values: const float*                ; A float array which can be indexed;
  */

  try {
    LOG(debug) << "Shape of input (tensor): " << printShape(input[0].GetTensorTypeAndShapeInfo().GetShape());

    auto outputTensors = mSession->Run(mInputNames, input, mOutputNames);
    float* output_values = outputTensors[0].GetTensorMutableData<float>();
    LOG(debug) << "Shape of output (tensor): " << printShape(outputTensors[0].GetTensorTypeAndShapeInfo().GetShape());

    // std::vector<float> output_vals(std::begin(output_values), std::end(output_values));

    return output_values;

  } catch (const Ort::Exception& exception) {
    LOG(error) << "Error running model inference: " << exception.what();

    return 0;
  }

} // function Network::evalNetwork(std::vector<Ort::Value>)

float* Network::evalNetwork(std::vector<float> input)
{

  /*
 Function: Evaluating the network for a std::vector<float>
 - Input:
   -- input:         std::vector<float>    ; The vector which should be evaluated by the network.
                                             The network will evaluate n inputs where n = input.size()/input_nodes with input_nodes = #(input neurons);
 - Output:
   -- output_values: const float*          ; A float array which can be indexed;
 */

  int64_t size = input.size();
  std::vector<int64_t> input_shape{size / mInputShapes[0][1], mInputShapes[0][1]};
  std::vector<Ort::Value> inputTensors;
  inputTensors.emplace_back(Ort::Experimental::Value::CreateTensor<float>(input.data(), size, input_shape));

  try {

    LOG(debug) << "Shape of input (vector): " << printShape(input_shape);
    auto outputTensors = mSession->Run(mInputNames, inputTensors, mOutputNames);
    LOG(debug) << "Shape of output (tensor): " << printShape(outputTensors[0].GetTensorTypeAndShapeInfo().GetShape());
    float* output_values = outputTensors[0].GetTensorMutableData<float>();

    return output_values;

  } catch (const Ort::Exception& exception) {

    LOG(error) << "Error running model inference: " << exception.what();

    return 0;
  }

} // function Network::evalNetwork(std::vector<float>)

} // namespace o2::pid::tpc
