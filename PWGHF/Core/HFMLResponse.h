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

template <std::size_t nModels = 1, typename T = float>
class HFMLResponse
{
  public:
    HFMLResponse() = default;
    virtual ~HFMLResponse() = default;


    /// Access model from CCDB
    /// \param onnxFiles is an array of onnx file names, one for each bin
    /// \param ccdbApi is the CCDB API
    /// \param pathCCDB is the model path in CCDB
    /// \param timestampCCDB is the CCDB timestamp
    void accessModelFromCCDB(std::array<std::string, nModels> onnxFiles, o2::ccdb::CcdbApi& ccdbApi, std::string pathCCDB, long timestampCCDB) {
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

// TODO add cuts and bin configs in init method
    void init(std::array<std::string, nModels> paths, bool enableOptimizations, int threads = 0) {
      mPaths = paths;
      uint8_t counterModel{0};
      for (const auto& path : paths) {
        mNetworks[counterModel].initModel(path, enableOptimizations, threads);
        ++counterModel;
      }
    }

    void init(bool enableOptimizations, int threads = 0) {
      return init(mPaths, enableOptimizations, threads);
    }

    /// Get array with model prediction
    /// \param input a vector containing the values of features used in BDT model
    /// \return model prediction for each class and each model
    std::vector<T> getModelOutputValues(int nModel, std::vector<T>& input)
    {
        T* outputValues = mNetworks[nModel].evalModel(input);
        mOutputValues = (std::vector<T>)*outputValues;
        return mOutputValues;
    }

    template <typename T1, typename T2, typename T3>
    bool isSelectedML(const T1& nModel, const T2& cutDir, const T3& cutsML) {
      uint8_t iClass{0};
      for (const auto& outputValue : mOutputValues) {
        uint8_t dir = cutDir->at(iClass);
        if (dir != hf_cuts_ml::CutDirection::CutNot) {
          if (dir == hf_cuts_ml::CutDirection::CutGreater && outputValue > cutsML->get(nModel, iClass)) {
            return false;
          }
          else if (dir == hf_cuts_ml::CutDirection::CutSmaller && outputValue < cutsML->get(nModel, iClass)) {
            return false;
          }
        }
        ++iClass;
      }
      return true;
    }

  private:
    std::array<o2::ml::OnnxModel, nModels> mNetworks; // OnnxModel objects, one for each bin

    std::vector<T> mOutputValues; 

    // std::vector<double> mBinsLimits = {}; // bin limits of the variable (e.g. pT) used to select which model to use
    std::array<std::string, nModels> mPaths = {}; // paths to the models, one for each bin
    std::vector<o2::analysis::hf_cuts_ml::CutDirection> mCutDir = {}; // direction of the cuts on the model scores (no cut is also supported)
    // o2::framework::LabeledArray<double> mCuts = {}; // array of cut values to apply on the model scores
};

}  // namespace o2::analysis

#endif // O2_ANALYSIS_HFMLRESPONSE_H_