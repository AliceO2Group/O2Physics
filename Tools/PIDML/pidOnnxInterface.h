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
/// \brief A class that wraps PID ML ONNX model. See README.md for more detailed instructions.
///
/// \author Maja Kabus <mkabus@cern.ch>

#ifndef TOOLS_PIDML_PIDONNXINTERFACE_H_
#define TOOLS_PIDML_PIDONNXINTERFACE_H_

#include <string>
#include <array>
#include <set>
#include <vector>

#include "Framework/Array2D.h"
#include "Tools/PIDML/pidOnnxModel.h"

namespace pidml_pt_cuts
{
static constexpr int nPids = 6;
static constexpr int nCutVars = kNDetectors;
constexpr int pids[nPids] = {211, 321, 2212, -211, -321, -2212};
auto pids_v = std::vector<int>{pids, pids + nPids};
constexpr double certainties[nPids] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
auto certainties_v = std::vector<double>{certainties, certainties + nPids};

// default values for the cuts
constexpr double cuts[nPids][nCutVars] = {{0.0, 0.5, 0.8}, {0.0, 0.5, 0.8}, {0.0, 0.5, 0.8}, {0.0, 0.5, 0.8}, {0.0, 0.5, 0.8}, {0.0, 0.5, 0.8}};

// row labels
static const std::vector<std::string> pidLabels = {
  "211", "321", "2212", "0211", "0321", "02212"};
// column labels
static const std::vector<std::string> cutVarLabels = {
  "TPC", "TPC + TOF", "TPC + TOF + TRD"};
} // namespace pidml_pt_cuts

struct PidONNXInterface {
  PidONNXInterface(std::string& localPath, std::string& ccdbPath, bool useCCDB, o2::ccdb::CcdbApi& ccdbApi, uint64_t timestamp, std::vector<int> const& pids, o2::framework::LabeledArray<double> const& pTLimits, std::vector<double> const& minCertainties, bool autoMode) : mNPids{pids.size()}, mPTLimits{pTLimits}
  {
    if (pids.size() == 0) {
      LOG(fatal) << "PID ML Interface needs at least 1 output pid to predict";
    }
    std::set<int> tmp;
    for (auto& pid : pids) {
      if (!tmp.insert(pid).second) {
        LOG(fatal) << "PID M Interface: output pids cannot repeat!";
      }
    }

    std::vector<double> minCertaintiesFilled;
    if (autoMode) {
      fillDefaultConfiguration(minCertaintiesFilled);
    } else {
      if (minCertainties.size() != mNPids) {
        LOG(fatal) << "PID ML Interface: min certainties vector must be of the same size as the output pids vector";
      }
      minCertaintiesFilled = minCertainties;
    }
    for (std::size_t i = 0; i < mNPids; i++) {
      for (uint32_t j = 0; j < kNDetectors; j++) {
        mModels.emplace_back(localPath, ccdbPath, useCCDB, ccdbApi, timestamp, pids[i], (PidMLDetector)(kTPCOnly + j), minCertaintiesFilled[i]);
      }
    }
  }
  PidONNXInterface() = default;
  PidONNXInterface(PidONNXInterface&&) = default;
  PidONNXInterface& operator=(PidONNXInterface&&) = default;
  PidONNXInterface(const PidONNXInterface&) = delete;
  PidONNXInterface& operator=(const PidONNXInterface&) = delete;
  ~PidONNXInterface() = default;

  template <typename T>
  float applyModel(const T& track, int pid)
  {
    for (std::size_t i = 0; i < mNPids; i++) {
      if (mModels[i * kNDetectors].mPid == pid) {
        for (uint32_t j = 0; j < kNDetectors; j++) {
          if (track.pt() >= mPTLimits[i][j] && (j == kNDetectors - 1 || track.pt() < mPTLimits[i][j + 1])) {
            return mModels[i * kNDetectors + j].applyModel(track);
          }
        }
      }
    }
    LOG(error) << "No suitable PID ML model found for track: " << track.globalIndex() << " from collision: " << track.collision().globalIndex() << " and expected pid: " << pid;
    return -1.0f;
  }

  template <typename T>
  bool applyModelBoolean(const T& track, int pid)
  {
    for (std::size_t i = 0; i < mNPids; i++) {
      if (mModels[i * kNDetectors].mPid == pid) {
        for (uint32_t j = 0; j < kNDetectors; j++) {
          if (track.pt() >= mPTLimits[i][j] && (j == kNDetectors - 1 || track.pt() < mPTLimits[i][j + 1])) {
            return mModels[i * kNDetectors + j].applyModelBoolean(track);
          }
        }
      }
    }
    LOG(error) << "No suitable PID ML model found for track: " << track.globalIndex() << " from collision: " << track.collision().globalIndex() << " and expected pid: " << pid;
    return false;
  }

 private:
  void fillDefaultConfiguration(std::vector<double>& minCertainties)
  {
    // FIXME: A more sophisticated strategy should be based on pid values as well
    mPTLimits = o2::framework::LabeledArray{pidml_pt_cuts::cuts[0], pidml_pt_cuts::nPids, pidml_pt_cuts::nCutVars, pidml_pt_cuts::pidLabels, pidml_pt_cuts::cutVarLabels};
    minCertainties = std::vector<double>(mNPids, 0.5);
  }

  std::vector<PidONNXModel> mModels;
  std::size_t mNPids;
  o2::framework::LabeledArray<double> mPTLimits;
};
#endif // TOOLS_PIDML_PIDONNXINTERFACE_H_
