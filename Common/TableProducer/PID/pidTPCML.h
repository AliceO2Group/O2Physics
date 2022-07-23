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

#ifndef O2_PID_TPC_ML_H_
#define O2_PID_TPC_ML_H_

#include <vector>

#include "TSystem.h"

// O2 includes
#include "Framework/Logger.h"
#include "ReconstructionDataFormats/PID.h"
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>

namespace o2::pid::tpc
{
class Network
{

 public:
  // Constructor, destructor and copy-constructor
  Network() = default;
  Network(std::string, bool);
  ~Network() = default;

  // Operators
  Network& operator=(Network&); // copy a network into an existing class instance

  // Functions
  template <typename C, typename T>
  std::vector<float> createInputFromTrack(const C&, const T&, const uint8_t) const; // create a std::vector<float> with all the inputs for the network
  std::vector<Ort::Value> createTensor(std::vector<float>) const;                   // create a std::vector<Ort::Value> (= ONNX tensor) for model input
  float* evalNetwork(std::vector<Ort::Value>);                                      // evaluate the network on a std::vector<Ort::Value> (= ONNX tensor)
  float* evalNetwork(std::vector<float>);                                           // evaluate the network on a std::vector<float>

  // Getters & Setters
  int getInputDimensions() const { return mInputShapes[0][1]; };
  int getOutputDimensions() const { return mOutputShapes[0][1]; };
  void SetNClNormalization(const float nclnorm) { nClNorm = nclnorm; }
  const float GetNClNormalization() const { return nClNorm; }

 private:
  // Environment variables for the ONNX runtime
  std::shared_ptr<Ort::Env> mEnv = nullptr;
  std::shared_ptr<Ort::Experimental::Session> mSession = nullptr;
  Ort::SessionOptions sessionOptions;

  // Input & Output specifications of the loaded network
  std::vector<std::string> mInputNames;
  std::vector<std::vector<int64_t>> mInputShapes;
  std::vector<std::string> mOutputNames;
  std::vector<std::vector<int64_t>> mOutputShapes;

  float nClNorm = 152.f;

  // Internal function for printing the shape of tensors: See https://github.com/saganatt/PID_ML_in_O2 or O2Physics/Tools/PIDML/simpleApplyPidOnnxModel.cxx
  std::string printShape(const std::vector<int64_t>& v)
  {
    std::stringstream ss("");
    for (size_t i = 0; i < v.size() - 1; i++)
      ss << v[i] << "x";
    ss << v[v.size() - 1];
    return ss.str();
  }

  // Class version
  ClassDefNV(Network, 2);

}; // class Network

Network::Network(std::string path,
                 bool enableOptimization = true)
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

  LOG(debug) << "Network copied!";

  return *this;

} // Network& Network::operator=(const Network &)

template <typename C, typename T>
std::vector<float> Network::createInputFromTrack(const C& collision_it, const T& track, const uint8_t id) const
{

  /*
  Function: Creating a std::vector<float> from a track with the variables that the network has been trained on
  - Input:
    -- collisions_it    const C&            ;   An iterator of a collisions table of the form: soa::Join<aod::Collisions, aod::Mults>::iterator const& collision
    -- track:           const T&            ;   A track, typically from soa::Join<...> tables or of their iterators;
    -- id:              uint8_t             ;   The id of a particle used for the mass assignment with o2::track::pid_constants::sMasses[id];
  - Output:
    -- inputValues:     std::vector<float>  ;   A std::vector<float> with the input variables for the network;
  */

  const float p = track.tpcInnerParam();
  const float tgl = track.tgl();
  const float signed1Pt = track.signed1Pt();
  const float mass = o2::track::pid_constants::sMasses[id];
  const float multTPC = collision_it.multTPC() / 11000.;
  const float ncl = std::sqrt(nClNorm / track.tpcNClsFound());

  std::vector<float> inputValues{p, tgl, signed1Pt, mass, multTPC, ncl};

  return inputValues;

} // std::vector<float> Network::createInputFromTrack(const C&, const T&, uint8_t)

std::vector<Ort::Value> Network::createTensor(std::vector<float> input) const
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

#endif // O2_PID_TPC_ML_H_