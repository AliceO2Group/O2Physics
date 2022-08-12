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

// O2 includes
#include "ReconstructionDataFormats/PID.h"
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#include <vector>
#include <array>
#include <string>

namespace o2::pid::tpc
{
class Network
{

 public:
  // Constructor, destructor and copy-constructor
  Network() = default;
  Network(std::string, bool);
  Network(std::string, unsigned long, unsigned long, bool); // initialization with timestamps
  ~Network() = default;

  // Operators
  Network& operator=(Network&); // copy a network into an existing class instance

  // Functions
  template <typename C, typename T>
  std::array<float, 6> createInputFromTrack(const C&, const T&, const uint8_t, const float) const; // create a std::vector<float> with all the inputs for the network
  std::vector<Ort::Value> createTensor(std::array<float, 6>) const;                                // create a std::vector<Ort::Value> (= ONNX tensor) for model input
  float* evalNetwork(std::vector<Ort::Value>);                                                     // evaluate the network on a std::vector<Ort::Value> (= ONNX tensor)
  float* evalNetwork(std::vector<float>);                                                          // evaluate the network on a std::vector<float>

  // Getters & Setters
  int getInputDimensions() const { return mInputShapes[0][1]; }
  int getOutputDimensions() const { return mOutputShapes[0][1]; }
  unsigned long getValidityFrom() const { return valid_from; }
  unsigned long getValidityUntil() const { return valid_until; }
  void setValidityFrom(unsigned long t) { valid_from = t; }
  void setValidityUntil(unsigned long t) { valid_until = t; }

 private:
  // Range of validity in timestamps
  unsigned long valid_from = 0;
  unsigned long valid_until = 0;

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
  std::string printShape(const std::vector<int64_t>& v);

  // Class version
  ClassDefNV(Network, 3);

}; // class Network

} // namespace o2::pid::tpc

#endif // O2_PID_TPC_ML_H_