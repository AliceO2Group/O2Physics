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

/// \file HFMLResponse.h
/// \brief Class to compute the ML response for HF-analysis selections
/// \author Fabio Catalano <fabio.catalano@cern.ch>, Universita' and INFN Torino

#ifndef O2_ANALYSIS_HFMLRESPONSE_H_
#define O2_ANALYSIS_HFMLRESPONSE_H_

#include <vector>
#include <string>

// ML application
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>

#include "Framework/Array2D.h"
#include "PWGHF/Core/SelectorCuts.h"
#include <CCDB/BasicCCDBManager.h>
#include "Tools/ML/model.h"



namespace o2::analysis
{

template <typename T = float>
class HFMLResponse
{
  public:
    /// Default constructor
    HFMLResponse() = default;
    /// Constructor initializing cuts and cut directions
    HFMLResponse(const std::vector<double>& binsLimits, const o2::framework::LabeledArray<double>& cuts, const std::vector<int>& cutDir) {
      mBinsLimits = binsLimits;
      mCuts = cuts;
      mCutDir = cutDir;
 
      mNetworks = std::vector<o2::ml::OnnxModel>(mNModels);
      mPaths = std::vector<std::string>(mNModels);
      mNModels = binsLimits.size() - 1;
      mNClasses = (cutDir.size() >= 3) ? cutDir.size() : 2;
    }
    /// Default destructor
    virtual ~HFMLResponse() = default;


    /// Access model from CCDB
    /// \param onnxFiles is a vector of onnx file names, one for each bin
    /// \param ccdbApi is the CCDB API
    /// \param pathCCDB is the model path in CCDB
    /// \param timestampCCDB is the CCDB timestamp
    void accessModelFromCCDB(const std::vector<std::string>& onnxFiles, o2::ccdb::CcdbApi& ccdbApi, std::string pathCCDB, long timestampCCDB) {
      uint8_t counterModel{0};
      for (const auto& onnxFile : onnxFiles) {
        std::map<std::string, std::string> metadata;
        bool retrieveSuccess = ccdbApi.retrieveBlob(pathCCDB, ".", metadata, timestampCCDB, false, onnxFile);
        if (retrieveSuccess) {
          mPaths[counterModel] = onnxFile;
        } else {
          LOG(fatal) << "Error encountered while accessing the ML model from CCDB! Maybe the ML model doesn't exist yet for this runnumber/timestamp?";
        }
        ++counterModel;
      }
    }

    /// Initialize OnnxModels
    /// \param paths is a vector of onnx model paths
    /// \param enableOptimizations is a switch no enable optimizations
    /// \param threads
    void init(const std::vector<std::string>& paths, bool enableOptimizations, int threads = 0) {
      uint8_t counterModel{0};
      for (const auto& path : paths) {
        mPaths[counterModel] = path;
        mNetworks[counterModel].initModel(path, enableOptimizations, threads);
        ++counterModel;
      }
    }

    // Initialize OnnxModels when paths already defined
    /// \param paths is a vector of onnx model paths
    /// \param enableOptimizations is a switch no enable optimizations
    /// \param threads
    void init(bool enableOptimizations, int threads = 0) {
      return init(mPaths, enableOptimizations, threads);
    }

    /// Get vector with model predictions
    /// \param input a vector containing the values of features used in BDT model
    /// \param nModel is the model index
    /// \return model prediction for each class and each model
    template <typename T1, typename T2>
    std::vector<T> getModelOutput(T1& input, const T2& nModel)
    {
        T* outputPtr = mNetworks[nModel].evalModel(input);
        std::vector<T> output(outputPtr, outputPtr + mNClasses);
        return output;
    }

    /// ML selections 
    /// \param input is the input features
    /// \param nModel is the model index
    /// \return boolean telling if model predictions pass the cuts
    template <typename T1, typename T2>
    bool isSelectedML(T1& input, const T2& nModel) {
      auto output = getModelOutput(input, nModel);
      uint8_t iClass{0};
      for (const auto& outputValue : output) {
        uint8_t dir = mCutDir.at(iClass);
        if (dir != hf_cuts_ml::CutDirection::CutNot) {
          if (dir == hf_cuts_ml::CutDirection::CutGreater && outputValue > mCuts.get(nModel, iClass)) {
            return false;
          }
          else if (dir == hf_cuts_ml::CutDirection::CutSmaller && outputValue < mCuts.get(nModel, iClass)) {
            return false;
          }
        }
        ++iClass;
      }
      return true;
    }

    /// ML selections 
    /// \param input is the input features
    /// \param nModel is the model index
    /// \param output is a container to be filled with model output
    /// \return boolean telling if model predictions pass the cuts
    template <typename T1, typename T2>
    bool isSelectedML(T1& input, const T2& nModel, std::vector<T>& output) {
      output = getModelOutput(input, nModel);
      uint8_t iClass{0};
      for (const auto& outputValue : output) {
        uint8_t dir = mCutDir.at(iClass);
        if (dir != hf_cuts_ml::CutDirection::CutNot) {
          if (dir == hf_cuts_ml::CutDirection::CutGreater && outputValue > mCuts.get(nModel, iClass)) {
            return false;
          }
          else if (dir == hf_cuts_ml::CutDirection::CutSmaller && outputValue < mCuts.get(nModel, iClass)) {
            return false;
          }
        }
        ++iClass;
      }
      return true;
    }

    /// Tag BDT predicted class
    /// \param scores is a vector with BDT out scores
    /// \return 0 if rejected, otherwise bitmap with BIT(RecoDecay::OriginType::Prompt) and/or BIT(RecoDecay::OriginType::NonPrompt) on
    template <typename T1, typename T2>
    uint8_t tagBDT(const T1& scores, const T2& nModel)
    {
      uint8_t tag{0};

      if (mNClasses == 2) { // binary classification
        if (scores[0] > mCuts.get(nModel, RecoDecay::OriginType::Prompt - 1)) {
          SETBIT(tag, RecoDecay::OriginType::Prompt); // here Prompt stands for Signal
        }
      } else if (mNClasses == 3) { // 3-class classification
        if (scores[0] > mCuts.get(nModel, RecoDecay::OriginType::None)) {
          return tag;
        }
        if (scores[1] > mCuts.get(nModel, RecoDecay::OriginType::Prompt)) {
          SETBIT(tag, RecoDecay::OriginType::Prompt);
        }
        if (scores[2] > mCuts.get(nModel, RecoDecay::OriginType::NonPrompt)) {
          SETBIT(tag, RecoDecay::OriginType::NonPrompt);
        }
      }
      return tag;
    }

    // FIXME : getter for debugging
    /// Get pointer to model of index 0
    o2::ml::OnnxModel* getModel() {
      return &mNetworks[0];
    }

  private:
    std::vector<o2::ml::OnnxModel> mNetworks; // OnnxModel objects, one for each bin
    uint8_t mNModels = 1; // number of bins
    uint8_t mNClasses = 3; // number of model classes
    std::vector<double> mBinsLimits = {}; // bin limits of the variable (e.g. pT) used to select which model to use
    std::vector<std::string> mPaths = {""}; // paths to the models, one for each bin
    std::vector<int> mCutDir = {}; // direction of the cuts on the model scores (no cut is also supported)
    o2::framework::LabeledArray<double> mCuts = {}; // array of cut values to apply on the model scores
};

}  // namespace o2::analysis

#endif // O2_ANALYSIS_HFMLRESPONSE_H_
