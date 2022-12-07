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

#ifndef COMMON_TABLEPRODUCER_PID_PIDTPCML_H_
#define COMMON_TABLEPRODUCER_PID_PIDTPCML_H_

// C++ and system includes
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#include <vector>
#include <array>
#include <string>
#include <memory>

// O2 includes
#include "ReconstructionDataFormats/PID.h"

namespace o2::pid::tpc
{
class Network
{

 public:
  // Constructor, destructor and copy-constructor
  Network() = default;
  Network(std::string, bool, int);
  Network(std::string, uint64_t, uint64_t, bool, int); // initialization with timestamps
  ~Network() = default;

  // Operators
  Network& operator=(Network&); // copy a network into an existing class instance

  // Functions
  std::vector<Ort::Value> createTensor(std::array<float, 6>) const; // create a std::vector<Ort::Value> (= ONNX tensor) for model input
  float* evalNetwork(std::vector<Ort::Value>);                      // evaluate the network on a std::vector<Ort::Value> (= ONNX tensor)
  float* evalNetwork(std::vector<float>);                           // evaluate the network on a std::vector<float>

  // Getters & Setters
  int getInputDimensions() const { return mInputShapes[0][1]; }
  int getOutputDimensions() const { return mOutputShapes[0][1]; }
  uint64_t getValidityFrom() const { return valid_from; }
  uint64_t getValidityUntil() const { return valid_until; }
  void setValidityFrom(uint64_t t) { valid_from = t; }
  void setValidityUntil(uint64_t t) { valid_until = t; }

 private:
  // Range of validity in timestamps
  uint64_t valid_from = 0;
  uint64_t valid_until = 0;

  // Environment variables for the ONNX runtime
  std::shared_ptr<Ort::Env> mEnv = nullptr;
  std::shared_ptr<Ort::Experimental::Session> mSession = nullptr;
  Ort::SessionOptions sessionOptions;

  // Input & Output specifications of the loaded network
  std::vector<std::string> mInputNames;
  std::vector<std::vector<int64_t>> mInputShapes;
  std::vector<std::string> mOutputNames;
  std::vector<std::vector<int64_t>> mOutputShapes;

  // Internal function for printing the shape of tensors: See https://github.com/saganatt/PID_ML_in_O2 or O2Physics/Tools/PIDML/simpleApplyPidOnnxModel.cxx
  std::string printShape(const std::vector<int64_t>& v);

  // Class version
  ClassDefNV(Network, 4);

}; // class Network

} // namespace o2::pid::tpc

#endif // COMMON_TABLEPRODUCER_PID_PIDTPCML_H_
