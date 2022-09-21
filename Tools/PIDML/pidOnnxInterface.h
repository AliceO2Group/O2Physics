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

// What detectors use for a track:
// always baseDetector, then switch on subsequent detectors when each next limit is passed
// It's possible to use just one detector configuration if you set limits to very big values
struct PidConfig {
  int pid;
  PidMLDetector baseDetector;
  std::array<float, kNDetectors - 1> pTLimits;
  float minCertainty; // value in [0, 1] describing min certainty of the model required to accept a particle
};

struct PidONNXInterface {
  PidONNXInterface(std::string& localPath, std::string& ccdbPath, bool useCCDB, o2::ccdb::CcdbApi& ccdbApi, uint64_t timestamp, std::vector<int> const& pids, std::vector<PidConfig> const& configs, bool autoMode) : mAutoMode(autoMode), mConfigs{configs}
  {
    mModels.clear();
    if (autoMode) {
      fillDefaultConfiguration();
    }
    mModels.resize(kNDetectors * pids.size());
    for (int i = 0; i < pids.size(); i++) {
      for (int j = 0; j < kNDetectors; j++) {
        mModels[i * kNDetectors + j] = PidONNXModel(localPath, ccdbPath, useCCDB, ccdbApi, timestamp, configs[i].pid, kTPCOnly + j, configs[i].minCertainty);
      }
    }
  }

  template <typename T>
  float applyModel(const T& track, int pid)
  {
    for (int i = 0; i < mModels.size(); i += kNDetectors) {
      if (mModels[i].mPid == pid) {
        int j = mConfigs[i / kNDetectors].baseDetector;
        while (j < kNDetectors && track.pt() < mConfigs[i / kNDetectors].pTLimits[j]) {
          j++;
        }
        return mModels[i + j].applyModel(track);
      }
    }
  }

  template <typename T>
  bool applyModelBoolean(const T& track, int pid)
  {
    for (int i = 0; i < mModels.size(); i += kNDetectors) {
      if (mModels[i].mPid == pid) {
        int j = mConfigs[i / kNDetectors].baseDetector;
        while (j < kNDetectors && track.pt() < mConfigs[i / kNDetectors].pTLimits[j]) {
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
      mConfigs[i].baseDetector = kTPCOnly;
      mConfigs[i].pTLimits = {0.5, 0.8};
      mConfigs[i].minCertainty = 0.5;
    }
  }

  // static constexpr int mNumKinds = 10;
  //  pions, protons, kaons, electrons, muons, and their antiparticles
  // static constexpr int mPdgs[mNumKinds] = {211, -211, 2212, -2212, 321, -321, 11, -11, 13, -13};
  // static constexpr int mPdgs[mNumKinds] = {11, 13, 211, 321, 2212, -11, -13, -211, -321, -2212};
  std::vector<PidONNXModel> mModels;
  std::vector<PidConfig> mConfigs;
  bool mAutoMode;
};
#endif // O2_ANALYSIS_PIDONNXINTERFACE_H_
