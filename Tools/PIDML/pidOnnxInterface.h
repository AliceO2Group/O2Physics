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

/// \file pidONNXInterface.h
/// \brief A class that wraps PID ML ONNX model. See https://github.com/saganatt/PID_ML_in_O2 for more detailed instructions.
///
/// \author Maja Kabus <mkabus@cern.ch>

#ifndef O2_ANALYSIS_PIDONNXINTERFACE_H_
#define O2_ANALYSIS_PIDONNXINTERFACE_H_

#include "Tools/PIDML/pidOnnxModel.h"

#include <string>
#include <array>

struct PidConfig {
  int pid;
  std::array<float, kNDetectors - 1> pTLimits; // What detectors use for a track: switch on subsequent detectors when each next pT limit is passed
  float minCertainty;                          // value in [0, 1] describing min certainty of the model required to accept a particle
};

struct PidONNXInterface {
  PidONNXInterface(std::string& localPath, std::string& ccdbPath, bool useCCDB, o2::ccdb::CcdbApi& ccdbApi, uint64_t timestamp, std::vector<int> const& pids, std::vector<std::array<float, kNDetectors - 1>> const& pTLimits, std::vector<float> const& minCertainties, bool autoMode) : mAutoMode(autoMode), mPids{pids}, mPTLimits{pTLimits}, mMinCertainties{minCertainties}
  {
    if (pids.size() == 0) {
      LOG(fatal) << "PID ML Interface needs at least 1 output pid to predict";
    }

    mModels.clear();
    if (autoMode) {
      fillDefaultConfiguration();
    }
    mModels.resize(kNDetectors * pids.size());
    for (int i = 0; i < pids.size(); i++) {
      for (int j = 0; j < kNDetectors; j++) {
        mModels[i * kNDetectors + j] = PidONNXModel(localPath, ccdbPath, useCCDB, ccdbApi, timestamp, pids[i], kTPCOnly + j, minCertainties[i]);
      }
    }
  }

  template <typename T>
  float applyModel(const T& track, int pid)
  {
    for (int i = 0; i < mModels.size(); i += kNDetectors) {
      if (mModels[i].mPid == pid) {
        int j = 0;
        while (j < kNDetectors && track.pt() < mPTLimits[i / kNDetectors][j]) {
          j++;
        }
        return mModels[i + j - 1].applyModel(track);
      }
    }
  }

  template <typename T>
  bool applyModelBoolean(const T& track, int pid)
  {
    for (int i = 0; i < mModels.size(); i += kNDetectors) {
      if (mModels[i].mPid == pid) {
        int j = 0;
        while (j < kNDetectors && track.pt() < mPTLimits[i / kNDetectors][j]) {
          j++;
        }
        return mModels[i + j - 1].applyModelBoolean(track);
      }
    }
  }

 private:
  void fillDefaultConfiguration()
  {
    // FIXME: A more sophisticated strategy should be based on pid values as well
    for (int i = 0; i < mConfigs.size(); i++) {
      mPTLimits[i] = {0.5, 0.8};
      mMinCertainties[i] = 0.5;
    }
  }

  // static constexpr int mNumKinds = 10;
  //  pions, protons, kaons, electrons, muons, and their antiparticles
  // static constexpr int mPdgs[mNumKinds] = {211, -211, 2212, -2212, 321, -321, 11, -11, 13, -13};
  // static constexpr int mPdgs[mNumKinds] = {11, 13, 211, 321, 2212, -11, -13, -211, -321, -2212};
  std::vector<PidONNXModel> mModels;
  std::vector<int> mPids;
  std::vector<std::array<float, kNDetectors - 1>> mPTLimits;
  std::vector<float> mMinCertainties;
  bool mAutoMode;
};
#endif // O2_ANALYSIS_PIDONNXINTERFACE_H_
