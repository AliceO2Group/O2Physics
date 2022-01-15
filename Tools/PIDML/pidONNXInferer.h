// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file pidONNXInferer.h
/// \brief A class that wraps PID ML ONNX model inference
///
/// \author Maja Kabus <mkabus@cern.ch>

#ifndef O2_ANALYSIS_PIDONNXINFERER_H_
#define O2_ANALYSIS_PIDONNXINFERER_H_

#include "CCDB/BasicCCDBManager.h"

#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <string>
#include <chrono>

// TODO: Copied from cefpTask, shall we put it in some common utils code?
namespace
{
bool readJsonFile(std::string& config, rapidjson::Document& d)
{
  FILE* fp = fopen(config.data(), "rb");
  if (!fp) {
    LOG(warning) << "Missing configuration json file: " << config;
    return false;
  }

  char readBuffer[65536];
  rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  d.ParseStream(is);
  fclose(fp);
  return true;
}
} // namespace

struct PidONNXInferer {
 public:
  PidONNXInferer(const std::string& modelFile, const std::string& trainColumnsFile, const std::string& scalingParamsFile, const std::string& url, long nolaterthan, bool useGPU)
  {
    TString* onnxModel = nullptr;
    loadInputFiles(modelFile, trainColumnsFile, scalingParamsFile, url, nolaterthan, onnxModel);

    Ort::SessionOptions sessionOptions;
    if (useGPU) {
      LOG(fatal) << "GPU not supported yet in PID ONNX Inferer" << std::endl;
      // FIXME: This can be used only if ONNXRuntime is build with CUDA support
      // Do we need GPU for inference?
      //Ort::ThrowOnError(OrtSessionOptionsAppendExecutionProvider_CUDA(sessionOptions, 0));
    } else {
      //mMemoryInfo = std::make_shared<Ort::MemoryInfo>(Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator, OrtMemTypeCPU));
    }
    mEnv = std::make_shared<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "pid-onnx-inferer");
    void* onnxData = (void*)onnxModel->Data(); // to get rid of const in TString::Data
    mSession = std::make_shared<Ort::Experimental::Session>(*mEnv, onnxData, onnxModel->Length(), sessionOptions);
    // TODO: This will be used when we will have direct access to the ONNX model file
    //mSession = Ort::Session{mEnv, modelFile, sessionOptions};
    //
    mInputNames = mSession->GetInputNames();
    mInputShapes = mSession->GetInputShapes();
    mOutputNames = mSession->GetOutputNames();
    mOutputShapes = mSession->GetOutputShapes();

    // print name/shape of inputs
    LOG(info) << "Input Node Name/Shape (" << mInputNames.size() << "):";
    for (size_t i = 0; i < mInputNames.size(); i++) {
      LOG(info) << "\t" << mInputNames[i] << " : " << printShape(mInputShapes[i]);
    }

    // print name/shape of outputs
    LOG(info) << "Output Node Name/Shape (" << mOutputNames.size() << "):";
    for (size_t i = 0; i < mOutputNames.size(); i++) {
      LOG(info) << "\t" << mOutputNames[i] << " : " << printShape(mOutputShapes[i]);
    }

    // Assume model has 1 input node and 1 output node.
    assert(mInputNames.size() == 1 && mOutputNames.size() == 1);
  }
  PidONNXInferer() = default;
  PidONNXInferer(PidONNXInferer& other) = default; // Ort::Env and unique pointers aren't copyable
  ~PidONNXInferer() = default;

  template <typename T>
  float applyModel(const T& track)
  {
    auto input_shape = mInputShapes[0];
    // TODO: We have different data types in the input but vector/tensor must have a single type. Narrowing conversion warnings.
    std::vector<float> inputTensorValues = createInputsSingle(track);
    std::vector<Ort::Value> inputTensors;
    if (mMemoryInfo != nullptr) {
      inputTensors.emplace_back(Ort::Value::CreateTensor<float>(*mMemoryInfo, inputTensorValues.data(), inputTensorValues.size(), input_shape.data(), input_shape.size()));
    } else {
      inputTensors.emplace_back(Ort::Experimental::Value::CreateTensor<float>(inputTensorValues.data(), inputTensorValues.size(), input_shape));
    }

    // double-check the dimensions of the input tensor
    assert(inputTensors[0].IsTensor() &&
           inputTensors[0].GetTensorTypeAndShapeInfo().GetShape() == input_shape);
    LOG(info) << "input tensor shape: " << printShape(inputTensors[0].GetTensorTypeAndShapeInfo().GetShape());

    // pass data through model
    LOG(info) << "Running model...";
    try {
      auto outputTensors = mSession->Run(mInputNames, inputTensors, mOutputNames);
      LOG(info) << "done";
      LOG(info) << "Number of output tensors: " << outputTensors.size();

      // double-check the dimensions of the output tensors
      // NOTE: the number of output tensors is equal to the number of output nodes specifed in the Run() call
      assert(outputTensors.size() == mOutputNames.size() && outputTensors[0].IsTensor());
      LOG(info) << "output tensor shape: " << printShape(outputTensors[0].GetTensorTypeAndShapeInfo().GetShape());

      for (auto& output : outputTensors) {
        const float* output_value = output.GetTensorData<float>();
        LOG(info) << "output: " << *output_value;
      }
      const float* output_value = outputTensors[0].GetTensorData<float>();
      return *output_value;

    } catch (const Ort::Exception& exception) {
      LOG(error) << "Error running model inference: " << exception.what();
    }
    return -1.0f; // unreachable code
  }

 private:
  // TODO: Temporary quick solution with CCDB and TStrings in ROOT files
  void loadInputFiles(const std::string& modelFile, const std::string& trainColumnsFile, const std::string& scalingParamsFile, const std::string& url, long nolaterthan, TString* onnxModel)
  {
    o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
    // Set CCDB url
    ccdb->setURL(url);
    // Enabling object caching, otherwise each call goes to the CCDB server
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    // Not later than now, will be replaced by the value of the train creation
    // This avoids that users can replace objects **while** a train is running
    ccdb->setCreatedNotAfter(nolaterthan);

    LOGF(info, "Getting object %s", modelFile.data());
    onnxModel = ccdb->getForTimeStamp<TString>(modelFile, nolaterthan);
    if (!onnxModel) {
      LOGF(fatal, "ONNX model %s not found!", modelFile);
    }
    LOGF(info, "Getting object %s", trainColumnsFile.data());
    TString* trainColumnsJson = ccdb->getForTimeStamp<TString>(trainColumnsFile, nolaterthan);
    if (!trainColumnsJson) {
      LOGF(fatal, "Train columns file %s not found!", trainColumnsFile);
    }
    LOGF(info, "Getting object %s", scalingParamsFile.data());
    TString* scalingParamsJson = ccdb->getForTimeStamp<TString>(scalingParamsFile, nolaterthan);
    if (!scalingParamsJson) {
      LOGF(fatal, "Scaling parameters file %s not found!", scalingParamsFile);
    }

    rapidjson::Document trainColumnsDoc;
    rapidjson::Document scalingParamsDoc;
    trainColumnsDoc.Parse(trainColumnsJson->Data());
    for (auto& param : trainColumnsDoc["columns_for_training"].GetArray()) {
      mTrainColumns.emplace_back(param.GetString());
    }
    scalingParamsDoc.Parse(scalingParamsJson->Data());
    for (auto& param : scalingParamsDoc["data"].GetArray()) {
      mScalingParams[param[0].GetString()] = std::make_pair(param[1].GetFloat(), param[2].GetFloat());
    }

    // TODO: This will be used when we will have direct access to JSON files
    //if (readJsonFile(trainColumnsFile, trainColumnsDoc)) {
    //  mTrainColumns = trainColumnsDoc["columns_for_training"].GetArray();
    //}
    //if (readJsonFile(scalingParamsFile, scalingParamsDoc)) {
    //  for (auto& param : scalingParamsDoc["data"].GetArray()) {
    //    mScalingParams[param[0]] = std::make_pair(param[1], param[2]);
    //  }
    //}
  }

  // TODO: Any more elegant way to select columns? This doesn't compile.
#define GET_VALUE_FOR_COLUMN(_TableIt_, _ColumnName_, _Val_) \
  _Val_ = static_cast<o2::aod::track::#_ColumnName_>(_TableIt_).getIterator().mCurrentPos;

  template <typename T>
  std::vector<float> createInputsSingle(const T& track)
  {
    std::vector<float> inputValues{track.tpcSignal(), track.beta(), track.px(), track.py(), track.pz(), track.sign(), track.x(), track.y(), track.z(), track.alpha(), track.trackType(), track.tpcNClsShared(), track.dcaXY(), track.dcaZ(), track.p()};
    //for (auto& columnName : mTrainColumns) {
    //  float* val;
    //  GET_VALUE_FOR_COLUMN(track, columnName, val);
    //  inputValues.push_back(*val);
    //}
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

  std::vector<std::string> mTrainColumns;
  std::map<std::string, std::pair<float, float>> mScalingParams;

  std::shared_ptr<Ort::Env> mEnv = nullptr;
  // No empty constructors for Session and MemoryInfo, we need pointers
  std::shared_ptr<Ort::Experimental::Session> mSession = nullptr;
  // Optional - manage memory (CPU vs GPU)
  std::shared_ptr<Ort::MemoryInfo> mMemoryInfo = nullptr;

  std::vector<std::string> mInputNames;
  std::vector<std::vector<int64_t>> mInputShapes;
  std::vector<std::string> mOutputNames;
  std::vector<std::vector<int64_t>> mOutputShapes;
};

#endif // O2_ANALYSIS_PIDONNXINFERER_H_
