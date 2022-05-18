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
/// \brief A class that wraps PID ML ONNX model. See https://github.com/saganatt/PID_ML_in_O2 for more detailed instructions.
///
/// \author Maja Kabus <mkabus@cern.ch>

#ifndef O2_ANALYSIS_PIDONNXMODEL_H_
#define O2_ANALYSIS_PIDONNXMODEL_H_

#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <string>

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
  PidONNXModel(const std::string& scalingParamsFile, int pid = 211, bool useTOF = false)
  {
    std::string modelFile;
    loadInputFiles(scalingParamsFile, useTOF, pid, modelFile);

    Ort::SessionOptions sessionOptions;
    mEnv = std::make_shared<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "pid-onnx-inferer");
    mSession.reset(new Ort::Experimental::Session{*mEnv, modelFile, sessionOptions});

    mInputNames = mSession->GetInputNames();
    mInputShapes = mSession->GetInputShapes();
    mOutputNames = mSession->GetOutputNames();
    mOutputShapes = mSession->GetOutputShapes();

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
  PidONNXModel(PidONNXModel& other) = default;
  ~PidONNXModel() = default;

  template <typename T>
  float applyModel(const T& track)
  {
    auto input_shape = mInputShapes[0];
    std::vector<float> inputTensorValues = createInputsSingle(track);
    std::vector<Ort::Value> inputTensors;
    inputTensors.emplace_back(Ort::Experimental::Value::CreateTensor<float>(inputTensorValues.data(), inputTensorValues.size(), input_shape));

    // Double-check the dimensions of the input tensor
    assert(inputTensors[0].IsTensor() &&
           inputTensors[0].GetTensorTypeAndShapeInfo().GetShape() == input_shape);
    LOG(debug) << "input tensor shape: " << printShape(inputTensors[0].GetTensorTypeAndShapeInfo().GetShape());

    try {
      auto outputTensors = mSession->Run(mInputNames, inputTensors, mOutputNames);

      // Double-check the dimensions of the output tensors
      // The number of output tensors is equal to the number of output nodes specifed in the Run() call
      assert(outputTensors.size() == mOutputNames.size() && outputTensors[0].IsTensor());
      LOG(debug) << "output tensor shape: " << printShape(outputTensors[0].GetTensorTypeAndShapeInfo().GetShape());

      const float* output_value = outputTensors[0].GetTensorData<float>();
      return *output_value;
    } catch (const Ort::Exception& exception) {
      LOG(error) << "Error running model inference: " << exception.what();
    }
    return -1.0f; // unreachable code
  }

 private:
  void loadInputFiles(const std::string& scalingParamsFile, bool useTOF, int pid, std::string& modelFile)
  {
    rapidjson::Document trainColumnsDoc;
    rapidjson::Document scalingParamsDoc;

    char* mlmodelsDir = getenv("MLMODELS_ROOT");
    if (mlmodelsDir == NULL) {
      LOG(fatal) << "Path to ML models undefined, did you load MLModels environment?";
    }

    std::string modelSubdir = useTOF ? "All" : "TPC";
    std::ostringstream tmp;
    tmp << mlmodelsDir << "/models/PID_ML/";
    tmp << modelSubdir << "/simple_model_" << pid << ".onnx";
    modelFile = tmp.str();
    tmp.str("");
    tmp.clear();
    tmp << mlmodelsDir << "/models/PID_ML/";
    tmp << modelSubdir << "/columns_for_training.json";
    std::string trainColumnsFile = tmp.str();
    tmp.str("");
    tmp.clear();
    tmp << mlmodelsDir << "/models/PID_ML/";
    tmp << scalingParamsFile;
    std::string scalingParamsFilePath = tmp.str();

    if (readJsonFile(trainColumnsFile, trainColumnsDoc)) {
      for (auto& param : trainColumnsDoc["columns_for_training"].GetArray()) {
        mTrainColumns.emplace_back(param.GetString());
      }
    }
    if (readJsonFile(scalingParamsFilePath, scalingParamsDoc)) {
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
    float scaledTOFSignal = (track.tofSignal() - mScalingParams.at("fTOFSignal").first) / mScalingParams.at("fTOFSignal").second;
    float scaledX = (track.x() - mScalingParams.at("fX").first) / mScalingParams.at("fX").second;
    float scaledY = (track.y() - mScalingParams.at("fY").first) / mScalingParams.at("fY").second;
    float scaledZ = (track.z() - mScalingParams.at("fZ").first) / mScalingParams.at("fZ").second;
    float scaledAlpha = (track.alpha() - mScalingParams.at("fAlpha").first) / mScalingParams.at("fAlpha").second;
    float scaledTPCNClsShared = ((float)track.tpcNClsShared() - mScalingParams.at("fTPCNClsShared").first) / mScalingParams.at("fTPCNClsShared").second;
    float scaledDcaXY = (track.dcaXY() - mScalingParams.at("fDcaXY").first) / mScalingParams.at("fDcaXY").second;
    float scaledDcaZ = (track.dcaZ() - mScalingParams.at("fDcaZ").first) / mScalingParams.at("fDcaZ").second;

    std::vector<float> inputValues{scaledTPCSignal, scaledTOFSignal, track.beta(), track.px(), track.py(), track.pz(), (float)track.sign(), scaledX, scaledY, scaledZ, scaledAlpha, (float)track.trackType(), scaledTPCNClsShared, scaledDcaXY, scaledDcaZ, track.p()};

    return inputValues;
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

  std::string mModelDir;
  std::vector<std::string> mTrainColumns;
  std::map<std::string, std::pair<float, float>> mScalingParams;

  std::shared_ptr<Ort::Env> mEnv = nullptr;
  // No empty constructors for Session, we need a pointer
  std::shared_ptr<Ort::Experimental::Session> mSession = nullptr;

  std::vector<std::string> mInputNames;
  std::vector<std::vector<int64_t>> mInputShapes;
  std::vector<std::string> mOutputNames;
  std::vector<std::vector<int64_t>> mOutputShapes;
};

#endif // O2_ANALYSIS_PIDONNXMODEL_H_
