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
/// \brief A class that wraps PID ML ONNX model inference. See https://github.com/saganatt/PID_ML_in_O2 for more detailed instructions.
///
/// \author Maja Kabus <mkabus@cern.ch>

#ifndef O2_ANALYSIS_PIDONNXINFERER_H_
#define O2_ANALYSIS_PIDONNXINFERER_H_

#include "CCDB/BasicCCDBManager.h"

#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <string>
#include <cctype>
#include <chrono>

// Temporary solution for Hyperloop tests - model files placed on CCDB. In the future, they will be taken probably from cvmfs (to be discussed).
#ifndef USE_CCDB
#define USE_CCDB 0
#endif

#if USE_CCDB == 0
// TODO: Copied from cefpTask, shall we put it in some common utils code?
namespace
{
bool readJsonFile(const std::string& config, rapidjson::Document& d)
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
#endif

struct PidONNXInferer {
 public:
  PidONNXInferer(std::string& modelFile, const std::string& trainColumnsFile, const std::string& scalingParamsFile, const std::string& url, long nolaterthan, bool useGPU) : mUseGPU(useGPU)
  {
    TString* onnxModel = nullptr;
    loadInputFiles(modelFile, trainColumnsFile, scalingParamsFile, url, nolaterthan, onnxModel);
#if USE_CCDB == 1
    //if (onnxModel) {
    //  LOG(info) << "Loaded ONNX model: " << onnxModel->Data();
    //} else {
    //  LOG(error) << "ONNX model is null!";
    //}
#endif

    Ort::SessionOptions sessionOptions;
    if (useGPU) {
      LOG(fatal) << "GPU not supported yet in PID ONNX Inferer" << std::endl;
      // FIXME: This can be used only if ONNXRuntime is build with CUDA support
      // Do we need GPU for inference?
      //Ort::ThrowOnError(OrtSessionOptionsAppendExecutionProvider_CUDA(sessionOptions, 0));
    } else {
      LOG(info) << "Using CPU";
      // Memory info is set each time in Ort::Experimental::Value::CreateTensor()
    }
    LOG(info) << "Creating ONNX env";
    mEnv = std::make_shared<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "pid-onnx-inferer");
#if USE_CCDB == 1
    char* onnxData = (char*)onnxModel->Data(); // to get rid of const in TString::Data
    mSession = std::make_shared<Ort::Experimental::Session>(*mEnv, onnxData, onnxModel->Length() * sizeof(char), sessionOptions);
#else
    mSession.reset(new Ort::Experimental::Session{*mEnv, modelFile, sessionOptions});
#endif

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
  PidONNXInferer() = default;
  PidONNXInferer(PidONNXInferer& other) = default;
  ~PidONNXInferer() = default;

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
  void loadInputFiles(std::string& modelFile, const std::string& trainColumnsFile, const std::string& scalingParamsFile, const std::string& url, long nolaterthan, TString* onnxModel)
  {
    rapidjson::Document trainColumnsDoc;
    rapidjson::Document scalingParamsDoc;

#if USE_CCDB == 1
    // TODO: Temporary quick solution with CCDB and TStrings in ROOT files
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

    trainColumnsDoc.Parse(trainColumnsJson->Data());
    for (auto& param : trainColumnsDoc["columns_for_training"].GetArray()) {
      mTrainColumns.emplace_back(param.GetString());
    }
    scalingParamsDoc.Parse(scalingParamsJson->Data());
    for (auto& param : scalingParamsDoc["data"].GetArray()) {
      mScalingParams[param[0].GetString()] = std::make_pair(param[1].GetFloat(), param[2].GetFloat());
    }

#else
    // TODO: This will be used when we will have direct access to JSON files
    if (readJsonFile(trainColumnsFile, trainColumnsDoc)) {
      for (auto& param : trainColumnsDoc["columns_for_training"].GetArray()) {
        mTrainColumns.emplace_back(param.GetString());
      }
    }
    if (readJsonFile(scalingParamsFile, scalingParamsDoc)) {
      for (auto& param : scalingParamsDoc["data"].GetArray()) {
        mScalingParams[param[0].GetString()] = std::make_pair(param[1].GetFloat(), param[2].GetFloat());
      }
    }
#endif
  }

#define GET_VALUE_FOR_COLUMN(_TableIt_, _ColumnGetter_) \
  static_assert(false, _ColumnGetter_);                 \
  float INPUT = static_cast<o2::aod::track::#_ColumnGetter_>(_TableIt_).getIterator().mCurrentPos;

  template <typename T>
  std::vector<float> createInputsSingle(const T& track)
  {
    //TODO: Hardcoded for now, needs to be modifiable
    //sign is short, trackType and tpcNClsShared uint8_t
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

    // TODO: Any more elegant way to select columns? This doesn't compile, a macro cannot interpret a string.
    //for (auto& columnLabel : mTrainColumns) {
    //  std::string columnGetter = columnLabel.substr(1);
    //  columnGetter[0] = std::tolower((unsigned char)columnGetter[0]);
    //  GET_VALUE_FOR_COLUMN(track, columnGetter);
    //  auto scaleIt = std::find(mScalingParams.begin(), mScalingParams.end(), columnLabel);
    //  if (scaleIt != mScalingParams.end()) {
    //    INPUT = (INPUT - scaleIt->second.first) / scaleIt->second.second;
    //  }
    //  inputValues.push_back(INPUT);
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

  bool mUseGPU;
  std::shared_ptr<Ort::Env> mEnv = nullptr;
  // No empty constructors for Session, we need a pointer
  std::shared_ptr<Ort::Experimental::Session> mSession = nullptr;

  std::vector<std::string> mInputNames;
  std::vector<std::vector<int64_t>> mInputShapes;
  std::vector<std::string> mOutputNames;
  std::vector<std::vector<int64_t>> mOutputShapes;
};

#endif // O2_ANALYSIS_PIDONNXINFERER_H_
