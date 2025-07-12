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

/// \file pidOnnxInterface.h
/// \brief A class that wraps PID ML ONNX model. See README.md for more detailed instructions.
///
/// \author Maja Kabus <mkabus@cern.ch>

#ifndef TOOLS_PIDML_PIDONNXINTERFACE_H_
#define TOOLS_PIDML_PIDONNXINTERFACE_H_

#include "Tools/PIDML/pidOnnxModel.h"
//
#include <CCDB/CcdbApi.h>
#include <Framework/Array2D.h>
#include <Framework/Logger.h>

#include <cstddef>
#include <cstdint>
#include <set>
#include <string>
#include <vector>

namespace pidml_pt_cuts
{
static constexpr int NPids = 6;
static constexpr int NCutVars = kNDetectors;
constexpr int Pids[NPids] = {211, 321, 2212, -211, -321, -2212};
auto pidsV = std::vector<int>{Pids, Pids + NPids};
constexpr double Certainties[NPids] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
auto certaintiesV = std::vector<double>{Certainties, Certainties + NPids};

// default values for the cuts
constexpr double Cuts[NPids][NCutVars] = {{0.0, 0.5, 0.8}, {0.0, 0.5, 0.8}, {0.0, 0.5, 0.8}, {0.0, 0.5, 0.8}, {0.0, 0.5, 0.8}, {0.0, 0.5, 0.8}};
// row labels
static const std::vector<std::string> pidLabels = {
  "211", "321", "2212", "0211", "0321", "02212"};
// column labels
static const std::vector<std::string> cutVarLabels = {
  "TPC", "TPC + TOF", "TPC + TOF + TRD"};

} // namespace pidml_pt_cuts

template <typename T>
struct PidONNXInterface {
  PidONNXInterface(std::string& localPath, std::string& ccdbPath, bool useCCDB, o2::ccdb::CcdbApi& ccdbApi, uint64_t timestamp, std::vector<int> const& Pids, o2::framework::LabeledArray<double> const& pLimits, std::vector<double> const& minCertainties, bool autoMode) : mNPids{Pids.size()}, mPLimits{pLimits}
  {
    if (Pids.size() == 0) {
      LOG(fatal) << "PID ML Interface needs at least 1 output pid to predict";
    }
    std::set<int> tmp;
    for (const auto& pid : Pids) {
      if (!tmp.insert(pid).second) {
        LOG(fatal) << "PID ML Interface: output Pids cannot repeat!";
      }
    }

    std::vector<double> minCertaintiesFilled;
    if (autoMode) {
      fillDefaultConfiguration(minCertaintiesFilled);
    } else {
      if (minCertainties.size() != mNPids) {
        LOG(fatal) << "PID ML Interface: min Certainties vector must be of the same size as the output Pids vector";
      }
      minCertaintiesFilled = minCertainties;
    }
    for (std::size_t i = 0; i < mNPids; i++) {
      mModels.emplace_back(localPath, ccdbPath, useCCDB, ccdbApi, timestamp, Pids[i], minCertaintiesFilled[i], mPLimits[i]);
    }
  }
  PidONNXInterface() = default;
  PidONNXInterface(PidONNXInterface&&) = default;
  PidONNXInterface& operator=(PidONNXInterface&&) = default;
  PidONNXInterface(const PidONNXInterface&) = delete;
  PidONNXInterface& operator=(const PidONNXInterface&) = delete;
  ~PidONNXInterface() = default;

  float applyModel(const T::iterator& track, int pid)
  {
    for (std::size_t i = 0; i < mNPids; i++) {
      if (mModels[i].mPid == pid) {
        return mModels[i].applyModel(track);
      }
    }
    LOG(error) << "No suitable PID ML model found for track: " << track.globalIndex() << " from collision: " << track.collision().globalIndex() << " and expected pid: " << pid;
    return -1.0f;
  }

  bool applyModelBoolean(const T::iterator& track, int pid)
  {
    for (std::size_t i = 0; i < mNPids; i++) {
      if (mModels[i].mPid == pid) {
        return mModels[i].applyModelBoolean(track);
      }
    }
    LOG(error) << "No suitable PID ML model found for track: " << track.globalIndex() << " from collision: " << track.collision().globalIndex() << " and expected pid: " << pid;
    return false;
  }

 private:
  void fillDefaultConfiguration(std::vector<double>& minCertainties)
  {
    // FIXME: A more sophisticated strategy should be based on pid values as well
    mPLimits = o2::framework::LabeledArray{pidml_pt_cuts::Cuts[0], pidml_pt_cuts::NPids, pidml_pt_cuts::NCutVars, pidml_pt_cuts::pidLabels, pidml_pt_cuts::cutVarLabels};
    minCertainties = std::vector<double>(mNPids, 0.5);
  }

  std::vector<PidONNXModel<T>> mModels;
  std::size_t mNPids{0};
  o2::framework::LabeledArray<double> mPLimits;
};
#endif // TOOLS_PIDML_PIDONNXINTERFACE_H_
