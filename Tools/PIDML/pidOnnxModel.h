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

/// \file pidONNXModel.h
/// \brief A class that wraps PID ML ONNX model. See README.md for more detailed instructions.
///
/// \author Maja Kabus <mkabus@cern.ch>

#ifndef TOOLS_PIDML_PIDONNXMODEL_H_
#define TOOLS_PIDML_PIDONNXMODEL_H_

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <algorithm>
#include <map>
#include <utility>
#include <memory>
#include <vector>
#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#else
#include <onnxruntime_cxx_api.h>
#endif

#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include "CCDB/CcdbApi.h"
#include "Tools/PIDML/pidUtils.h"

using namespace pidml::pidutils;

enum PidMLDetector {
  kTPCOnly = 0,
  kTPCTOF,
  kTPCTOFTRD,
  kNDetectors ///< number of available detectors configurations
};

namespace pidml_pt_cuts
{
// TODO: for now first limit wouldn't be used,
// network needs TPC, so we can either do not cut it by p or return 0.0f as prediction
constexpr std::array<double, kNDetectors> defaultModelPLimits({0.0, 0.5, 0.8});
} // namespace pidml_pt_cuts

// TODO: Copied from cefpTask, shall we put it in some common utils code?
namespace
{
bool readJsonFile(const std::string& config, rapidjson::Document& d)
{
  FILE* fp = fopen(config.data(), "rb");
  if (!fp) {
    LOG(error) << "Missing configuration json file: " << config;
    return false;
  }

  char readBuffer[65536];
  rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  d.ParseStream(is);
  fclose(fp);
  return true;
}
} // namespace

struct PidONNXModel {
 public:
  PidONNXModel(std::string& localPath, std::string& ccdbPath, bool useCCDB, o2::ccdb::CcdbApi& ccdbApi, uint64_t timestamp,
               int pid, double minCertainty, const double* pLimits = &pidml_pt_cuts::defaultModelPLimits[0])
    : mPid(pid), mMinCertainty(minCertainty), mPLimits(pLimits, pLimits + kNDetectors)
  {
    assert(mPLimits.size() == kNDetectors);

    std::string modelFile;
    loadInputFiles(localPath, ccdbPath, useCCDB, ccdbApi, timestamp, pid, modelFile);

    Ort::SessionOptions sessionOptions;
    mEnv = std::make_shared<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "pid-onnx-inferer");
    LOG(info) << "Loading ONNX model from file: " << modelFile;
#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
    mSession.reset(new Ort::Experimental::Session{*mEnv, modelFile, sessionOptions});
#else
    mSession.reset(new Ort::Session{*mEnv, modelFile.c_str(), sessionOptions});
#endif
    LOG(info) << "ONNX model loaded";

#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
    mInputNames = mSession->GetInputNames();
    mInputShapes = mSession->GetInputShapes();
    mOutputNames = mSession->GetOutputNames();
    mOutputShapes = mSession->GetOutputShapes();
#else
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
#endif

    LOG(debug) << "Input Node Name/Shape (" << mInputNames.size() << "):";
    for (size_t i = 0; i < mInputNames.size(); i++) {
      LOG(debug) << "\t" << mInputNames[i] << " : " << printShape(mInputShapes[i]);
    }

    LOG(debug) << "Output Node Name/Shape (" << mOutputNames.size() << "):";
    for (size_t i = 0; i < mOutputNames.size(); i++) {
      LOG(debug) << "\t" << mOutputNames[i] << " : " << printShape(mOutputShapes[i]);
    }

    // Assume model has 1 input node and 1 output node.
    assert(mInputNames.size() == 1 && mOutputNames.size() == 1);
  }
  PidONNXModel() = default;
  PidONNXModel(PidONNXModel&&) = default;
  PidONNXModel& operator=(PidONNXModel&&) = default;
  PidONNXModel(const PidONNXModel&) = delete;
  PidONNXModel& operator=(const PidONNXModel&) = delete;
  ~PidONNXModel() = default;

  template <typename T>
  float applyModel(const T& track)
  {
    return getModelOutput(track);
  }

  template <typename T>
  bool applyModelBoolean(const T& track)
  {
    return getModelOutput(track) >= mMinCertainty;
  }

  int mPid;
  double mMinCertainty;

 private:
  void getModelPaths(std::string const& path, std::string& modelDir, std::string& modelFile, std::string& modelPath, int pid, std::string const& ext)
  {
    modelDir = path;
    modelFile = "attention_model_";

    if (pid < 0) {
      modelFile += "0" + std::to_string(-pid);
    } else {
      modelFile += std::to_string(pid);
    }

    modelFile += ext;
    modelPath = modelDir + "/" + modelFile;
  }

  void downloadFromCCDB(o2::ccdb::CcdbApi& ccdbApi, std::string const& ccdbFile, uint64_t timestamp, std::string const& localDir, std::string const& localFile)
  {
    std::map<std::string, std::string> metadata;
    bool retrieveSuccess = ccdbApi.retrieveBlob(ccdbFile, localDir, metadata, timestamp, false, localFile);
    if (retrieveSuccess) {
      std::map<std::string, std::string> headers = ccdbApi.retrieveHeaders(ccdbFile, metadata, timestamp);
      LOG(info) << "Network file downloaded from: " << ccdbFile << " to: " << localDir << "/" << localFile;
    } else {
      LOG(fatal) << "Error encountered while fetching/loading the network from CCDB! Maybe the network doesn't exist yet for this run number/timestamp?";
    }
  }

  void loadInputFiles(std::string const& localPath, std::string const& ccdbPath, bool useCCDB, o2::ccdb::CcdbApi& ccdbApi, uint64_t timestamp, int pid, std::string& modelPath)
  {
    rapidjson::Document trainColumnsDoc;
    rapidjson::Document scalingParamsDoc;

    std::string localDir, localModelFile;
    std::string trainColumnsFile = "columns_for_training";
    std::string scalingParamsFile = "scaling_params";
    getModelPaths(localPath, localDir, localModelFile, modelPath, pid, ".onnx");
    std::string localTrainColumnsPath = localDir + "/" + trainColumnsFile + ".json";
    std::string localScalingParamsPath = localDir + "/" + scalingParamsFile + ".json";

    if (useCCDB) {
      std::string ccdbDir, ccdbModelFile, ccdbModelPath;
      getModelPaths(ccdbPath, ccdbDir, ccdbModelFile, ccdbModelPath, pid, "");
      std::string ccdbTrainColumnsPath = ccdbDir + "/" + trainColumnsFile;
      std::string ccdbScalingParamsPath = ccdbDir + "/" + scalingParamsFile;
      downloadFromCCDB(ccdbApi, ccdbModelPath, timestamp, localDir, localModelFile);
      downloadFromCCDB(ccdbApi, ccdbTrainColumnsPath, timestamp, localDir, "columns_for_training.json");
      downloadFromCCDB(ccdbApi, ccdbScalingParamsPath, timestamp, localDir, "scaling_params.json");
    }

    LOG(info) << "Using configuration files: " << localTrainColumnsPath << ", " << localScalingParamsPath;
    if (readJsonFile(localTrainColumnsPath, trainColumnsDoc)) {
      for (auto& param : trainColumnsDoc["columns_for_training"].GetArray()) {
        mTrainColumns.emplace_back(param.GetString());
      }
    }
    if (readJsonFile(localScalingParamsPath, scalingParamsDoc)) {
      for (auto& param : scalingParamsDoc["data"].GetArray()) {
        mScalingParams[param[0].GetString()] = std::make_pair(param[1].GetFloat(), param[2].GetFloat());
      }
    }
  }

  template <typename T>
  std::vector<float> createInputsSingle(const T& track)
  {
    // TODO: Hardcoded for now. Planning to implement RowView extension to get runtime access to selected columns
    // sign is short, trackType and tpcNClsShared uint8_t

    float scaledTPCSignal = (track.tpcSignal() - mScalingParams.at("fTPCSignal").first) / mScalingParams.at("fTPCSignal").second;

    std::vector<float> inputValues{scaledTPCSignal};

    // When TRD Signal shouldn't be used we pass quiet_NaNs to the network
    if (!inPLimit(track, mPLimits[kTPCTOFTRD]) || trdMissing(track)) {
      inputValues.push_back(std::numeric_limits<float>::quiet_NaN());
      inputValues.push_back(std::numeric_limits<float>::quiet_NaN());
    } else {
      float scaledTRDSignal = (track.trdSignal() - mScalingParams.at("fTRDSignal").first) / mScalingParams.at("fTRDSignal").second;
      inputValues.push_back(scaledTRDSignal);
      inputValues.push_back(track.trdPattern());
    }

    // When TOF Signal shouldn't be used we pass quiet_NaNs to the network
    if (!inPLimit(track, mPLimits[kTPCTOF]) || tofMissing(track)) {
      inputValues.push_back(std::numeric_limits<float>::quiet_NaN());
      inputValues.push_back(std::numeric_limits<float>::quiet_NaN());
    } else {
      float scaledTOFSignal = (track.tofSignal() - mScalingParams.at("fTOFSignal").first) / mScalingParams.at("fTOFSignal").second;
      float scaledBeta = (track.beta() - mScalingParams.at("fBeta").first) / mScalingParams.at("fBeta").second;
      inputValues.push_back(scaledTOFSignal);
      inputValues.push_back(scaledBeta);
    }

    float scaledX = (track.x() - mScalingParams.at("fX").first) / mScalingParams.at("fX").second;
    float scaledY = (track.y() - mScalingParams.at("fY").first) / mScalingParams.at("fY").second;
    float scaledZ = (track.z() - mScalingParams.at("fZ").first) / mScalingParams.at("fZ").second;
    float scaledAlpha = (track.alpha() - mScalingParams.at("fAlpha").first) / mScalingParams.at("fAlpha").second;
    float scaledTPCNClsShared = (static_cast<float>(track.tpcNClsShared()) - mScalingParams.at("fTPCNClsShared").first) / mScalingParams.at("fTPCNClsShared").second;
    float scaledDcaXY = (track.dcaXY() - mScalingParams.at("fDcaXY").first) / mScalingParams.at("fDcaXY").second;
    float scaledDcaZ = (track.dcaZ() - mScalingParams.at("fDcaZ").first) / mScalingParams.at("fDcaZ").second;

    inputValues.insert(inputValues.end(), {track.p(), track.pt(), track.px(), track.py(), track.pz(), static_cast<float>(track.sign()), scaledX, scaledY, scaledZ, scaledAlpha, static_cast<float>(track.trackType()), scaledTPCNClsShared, scaledDcaXY, scaledDcaZ});

    return inputValues;
  }

  template <typename T>
  float getModelOutput(const T& track)
  {
    // First rank of the expected model input is -1 which means that it is dynamic axis.
    // Axis is exported as dynamic to make it possible to run model inference with the batch of
    // tracks at once in the future (batch would need to have the same amount of quiet_NaNs in each row).
    // For now we hardcode 1.
    static constexpr int64_t batch_size = 1;
    auto input_shape = mInputShapes[0];
    input_shape[0] = batch_size;

    std::vector<float> inputTensorValues = createInputsSingle(track);
    std::vector<Ort::Value> inputTensors;

#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
    inputTensors.emplace_back(Ort::Experimental::Value::CreateTensor<float>(inputTensorValues.data(), inputTensorValues.size(), input_shape));
#else
    Ort::MemoryInfo mem_info = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);
    inputTensors.emplace_back(Ort::Value::CreateTensor<float>(mem_info, inputTensorValues.data(), inputTensorValues.size(), input_shape.data(), input_shape.size()));
#endif

    // Double-check the dimensions of the input tensor
    assert(inputTensors[0].IsTensor() &&
           inputTensors[0].GetTensorTypeAndShapeInfo().GetShape() == input_shape);
    LOG(debug) << "input tensor shape: " << printShape(inputTensors[0].GetTensorTypeAndShapeInfo().GetShape());

    try {
#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
      auto outputTensors = mSession->Run(mInputNames, inputTensors, mOutputNames);
#else
      Ort::RunOptions runOptions;
      std::vector<const char*> inputNamesChar(mInputNames.size(), nullptr);
      std::transform(std::begin(mInputNames), std::end(mInputNames), std::begin(inputNamesChar),
                     [&](const std::string& str) { return str.c_str(); });

      std::vector<const char*> outputNamesChar(mOutputNames.size(), nullptr);
      std::transform(std::begin(mOutputNames), std::end(mOutputNames), std::begin(outputNamesChar),
                     [&](const std::string& str) { return str.c_str(); });
      auto outputTensors = mSession->Run(runOptions, inputNamesChar.data(), inputTensors.data(), inputTensors.size(), outputNamesChar.data(), outputNamesChar.size());
#endif

      // Double-check the dimensions of the output tensors
      // The number of output tensors is equal to the number of output nodes specified in the Run() call
      assert(outputTensors.size() == mOutputNames.size() && outputTensors[0].IsTensor());
      LOG(debug) << "output tensor shape: " << printShape(outputTensors[0].GetTensorTypeAndShapeInfo().GetShape());

      const float* output_value = outputTensors[0].GetTensorData<float>();
      float certainty = *output_value;
      return certainty;
    } catch (const Ort::Exception& exception) {
      LOG(error) << "Error running model inference: " << exception.what();
    }
    return false; // unreachable code
  }

  // Pretty prints a shape dimension vector
  std::string printShape(const std::vector<int64_t>& v)
  {
    std::stringstream ss("");
    for (size_t i = 0; i < v.size() - 1; i++)
      ss << v[i] << "x";
    ss << v[v.size() - 1];
    return ss.str();
  }

  std::vector<std::string> mTrainColumns;
  std::map<std::string, std::pair<float, float>> mScalingParams;

  std::shared_ptr<Ort::Env> mEnv = nullptr;
  // No empty constructors for Session, we need a pointer
#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
  std::shared_ptr<Ort::Experimental::Session> mSession = nullptr;
#else
  std::shared_ptr<Ort::Session> mSession = nullptr;
#endif

  std::vector<double> mPLimits;
  std::vector<std::string> mInputNames;
  std::vector<std::vector<int64_t>> mInputShapes;
  std::vector<std::string> mOutputNames;
  std::vector<std::vector<int64_t>> mOutputShapes;
};

#endif // TOOLS_PIDML_PIDONNXMODEL_H_
