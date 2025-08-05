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
//
// Contact: jseo@cern.ch
//
// Class to compute the ML response for DQ-analysis selections
//

#ifndef PWGDQ_CORE_DQMLRESPONSE_H_
#define PWGDQ_CORE_DQMLRESPONSE_H_

#include "Tools/ML/MlResponse.h"

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis
{

enum class InputFeatures : uint8_t { // refer to DielectronsAll, TODO: add more features if needed
  kMass = 0,
  kPt,
  kEta,
  kPhi,
  kPt1,
  kITSChi2NCl1,
  kTPCNClsCR1,
  kTPCNClsFound1,
  kTPCChi2NCl1,
  kDcaXY1,
  kDcaZ1,
  kTPCNSigmaEl1,
  kTPCNSigmaPi1,
  kTPCNSigmaPr1,
  kTOFNSigmaEl1,
  kTOFNSigmaPi1,
  kTOFNSigmaPr1,
  kPt2,
  kITSChi2NCl2,
  kTPCNClsCR2,
  kTPCNClsFound2,
  kTPCChi2NCl2,
  kDcaXY2,
  kDcaZ2,
  kTPCNSigmaEl2,
  kTPCNSigmaPi2,
  kTPCNSigmaPr2,
  kTOFNSigmaEl2,
  kTOFNSigmaPi2,
  kTOFNSigmaPr2
};

static const std::map<InputFeatures, std::string> gFeatureNameMap = {
  {InputFeatures::kMass, "kMass"},
  {InputFeatures::kPt, "kPt"},
  {InputFeatures::kEta, "kEta"},
  {InputFeatures::kPhi, "kPhi"},
  {InputFeatures::kPt1, "kPt1"},
  {InputFeatures::kITSChi2NCl1, "kITSChi2NCl1"},
  {InputFeatures::kTPCNClsCR1, "kTPCNClsCR1"},
  {InputFeatures::kTPCNClsFound1, "kTPCNClsFound1"},
  {InputFeatures::kTPCChi2NCl1, "kTPCChi2NCl1"},
  {InputFeatures::kDcaXY1, "kDcaXY1"},
  {InputFeatures::kDcaZ1, "kDcaZ1"},
  {InputFeatures::kTPCNSigmaEl1, "kTPCNSigmaEl1"},
  {InputFeatures::kTPCNSigmaPi1, "kTPCNSigmaPi1"},
  {InputFeatures::kTPCNSigmaPr1, "kTPCNSigmaPr1"},
  {InputFeatures::kTOFNSigmaEl1, "kTOFNSigmaEl1"},
  {InputFeatures::kTOFNSigmaPi1, "kTOFNSigmaPi1"},
  {InputFeatures::kTOFNSigmaPr1, "kTOFNSigmaPr1"},
  {InputFeatures::kPt2, "kPt2"},
  {InputFeatures::kITSChi2NCl2, "kITSChi2NCl2"},
  {InputFeatures::kTPCNClsCR2, "kTPCNClsCR2"},
  {InputFeatures::kTPCNClsFound2, "kTPCNClsFound2"},
  {InputFeatures::kTPCChi2NCl2, "kTPCChi2NCl2"},
  {InputFeatures::kDcaXY2, "kDcaXY2"},
  {InputFeatures::kDcaZ2, "kDcaZ2"},
  {InputFeatures::kTPCNSigmaEl2, "kTPCNSigmaEl2"},
  {InputFeatures::kTPCNSigmaPi2, "kTPCNSigmaPi2"},
  {InputFeatures::kTPCNSigmaPr2, "kTPCNSigmaPr2"},
  {InputFeatures::kTOFNSigmaEl2, "kTOFNSigmaEl2"},
  {InputFeatures::kTOFNSigmaPi2, "kTOFNSigmaPi2"},
  {InputFeatures::kTOFNSigmaPr2, "kTOFNSigmaPr2"}};

template <typename TypeOutputScore = float>
class DQMlResponse : public MlResponse<TypeOutputScore>
{
 public:
  DQMlResponse() = default;
  virtual ~DQMlResponse() = default;

  void setBinsCent(const std::vector<std::pair<double, double>>& bins) { binsCent = bins; }
  void setBinsPt(const std::vector<std::pair<double, double>>& bins) { binsPt = bins; }
  void setCentType(std::string& type) { centType = type; }

  const std::vector<std::pair<double, double>>& getBinsCent() const { return binsCent; }
  const std::vector<std::pair<double, double>>& getBinsPt() const { return binsPt; }
  const std::string& getCentType() const { return centType; }

  /// Method to get the input features vector needed for ML inference
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename TValues>
  std::vector<float> getInputFeatures(const T1& t1,
                                      const T2& t2,
                                      const TValues& fg) const
  {
    std::vector<float> dqInputFeatures;
    dqInputFeatures.reserve(MlResponse<TypeOutputScore>::mCachedIndices.size());

    for (auto idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      auto enumIdx = static_cast<InputFeatures>(idx);
      auto mapIdx = gFeatureNameMap.find(enumIdx);
      if (mapIdx == gFeatureNameMap.end()) {
        LOG(fatal) << "Unknown InputFeatures index: " << static_cast<int>(enumIdx);
      }

      const auto& name = mapIdx->second;
      if (name == "kMass") {
        dqInputFeatures.push_back(fg[VarManager::fgVarNamesMap["kMass"]]);
      } else if (name == "kPt") {
        dqInputFeatures.push_back(fg[VarManager::fgVarNamesMap["kPt"]]);
      } else if (name == "kEta") {
        dqInputFeatures.push_back(fg[VarManager::fgVarNamesMap["kEta"]]);
      } else if (name == "kPhi") {
        dqInputFeatures.push_back(fg[VarManager::fgVarNamesMap["kPhi"]]);
      } else if (name == "kPt1") {
        dqInputFeatures.push_back(t1.pt());
      } else if (name == "kITSChi2NCl1") {
        dqInputFeatures.push_back(t1.itsChi2NCl());
      } else if (name == "kTPCNClsCR1") {
        dqInputFeatures.push_back(t1.tpcNClsCrossedRows());
      } else if (name == "kTPCNClsFound1") {
        dqInputFeatures.push_back(t1.tpcNClsFound());
      } else if (name == "kTPCChi2NCl1") {
        dqInputFeatures.push_back(t1.tpcChi2NCl());
      } else if (name == "kDcaXY1") {
        dqInputFeatures.push_back(t1.dcaXY());
      } else if (name == "kDcaZ1") {
        dqInputFeatures.push_back(t1.dcaZ());
      } else if (name == "kTPCNSigmaEl1") {
        dqInputFeatures.push_back(t1.tpcNSigmaEl());
      } else if (name == "kTPCNSigmaPi1") {
        dqInputFeatures.push_back(t1.tpcNSigmaPi());
      } else if (name == "kTPCNSigmaPr1") {
        dqInputFeatures.push_back(t1.tpcNSigmaPr());
      } else if (name == "kTOFNSigmaEl1") {
        dqInputFeatures.push_back(t1.tofNSigmaEl());
      } else if (name == "kTOFNSigmaPi1") {
        dqInputFeatures.push_back(t1.tofNSigmaPi());
      } else if (name == "kTOFNSigmaPr1") {
        dqInputFeatures.push_back(t1.tofNSigmaPr());
      } else if (name == "kPt2") {
        dqInputFeatures.push_back(t2.pt());
      } else if (name == "kITSChi2NCl2") {
        dqInputFeatures.push_back(t2.itsChi2NCl());
      } else if (name == "kTPCNClsCR2") {
        dqInputFeatures.push_back(t2.tpcNClsCrossedRows());
      } else if (name == "kTPCNClsFound2") {
        dqInputFeatures.push_back(t2.tpcNClsFound());
      } else if (name == "kTPCChi2NCl2") {
        dqInputFeatures.push_back(t2.tpcChi2NCl());
      } else if (name == "kDcaXY2") {
        dqInputFeatures.push_back(t2.dcaXY());
      } else if (name == "kDcaZ2") {
        dqInputFeatures.push_back(t2.dcaZ());
      } else if (name == "kTPCNSigmaEl2") {
        dqInputFeatures.push_back(t2.tpcNSigmaEl());
      } else if (name == "kTPCNSigmaPi2") {
        dqInputFeatures.push_back(t2.tpcNSigmaPi());
      } else if (name == "kTPCNSigmaPr2") {
        dqInputFeatures.push_back(t2.tpcNSigmaPr());
      } else if (name == "kTOFNSigmaEl2") {
        dqInputFeatures.push_back(t2.tofNSigmaEl());
      } else if (name == "kTOFNSigmaPi2") {
        dqInputFeatures.push_back(t2.tofNSigmaPi());
      } else if (name == "kTOFNSigmaPr2") {
        dqInputFeatures.push_back(t2.tofNSigmaPr());
      } else {
        LOG(fatal) << "Missing accessor for feature: " << name;
      }
    }
    LOG(debug) << "Total features collected: " << dqInputFeatures.size();
    return dqInputFeatures;
  }

 protected:
  std::vector<std::pair<double, double>> binsCent;
  std::vector<std::pair<double, double>> binsPt;
  std::string centType;

  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      {"kMass", static_cast<uint8_t>(InputFeatures::kMass)},
      {"kPt", static_cast<uint8_t>(InputFeatures::kPt)},
      {"kEta", static_cast<uint8_t>(InputFeatures::kEta)},
      {"kPhi", static_cast<uint8_t>(InputFeatures::kPhi)},
      {"kPt1", static_cast<uint8_t>(InputFeatures::kPt1)},
      {"kITSChi2NCl1", static_cast<uint8_t>(InputFeatures::kITSChi2NCl1)},
      {"kTPCNClsCR1", static_cast<uint8_t>(InputFeatures::kTPCNClsCR1)},
      {"kTPCNClsFound1", static_cast<uint8_t>(InputFeatures::kTPCNClsFound1)},
      {"kTPCChi2NCl1", static_cast<uint8_t>(InputFeatures::kTPCChi2NCl1)},
      {"kDcaXY1", static_cast<uint8_t>(InputFeatures::kDcaXY1)},
      {"kDcaZ1", static_cast<uint8_t>(InputFeatures::kDcaZ1)},
      {"kTPCNSigmaEl1", static_cast<uint8_t>(InputFeatures::kTPCNSigmaEl1)},
      {"kTPCNSigmaPi1", static_cast<uint8_t>(InputFeatures::kTPCNSigmaPi1)},
      {"kTPCNSigmaPr1", static_cast<uint8_t>(InputFeatures::kTPCNSigmaPr1)},
      {"kTOFNSigmaEl1", static_cast<uint8_t>(InputFeatures::kTOFNSigmaEl1)},
      {"kTOFNSigmaPi1", static_cast<uint8_t>(InputFeatures::kTOFNSigmaPi1)},
      {"kTOFNSigmaPr1", static_cast<uint8_t>(InputFeatures::kTOFNSigmaPr1)},
      {"kPt2", static_cast<uint8_t>(InputFeatures::kPt2)},
      {"kITSChi2NCl2", static_cast<uint8_t>(InputFeatures::kITSChi2NCl2)},
      {"kTPCNClsCR2", static_cast<uint8_t>(InputFeatures::kTPCNClsCR2)},
      {"kTPCNClsFound2", static_cast<uint8_t>(InputFeatures::kTPCNClsFound2)},
      {"kTPCChi2NCl2", static_cast<uint8_t>(InputFeatures::kTPCChi2NCl2)},
      {"kDcaXY2", static_cast<uint8_t>(InputFeatures::kDcaXY2)},
      {"kDcaZ2", static_cast<uint8_t>(InputFeatures::kDcaZ2)},
      {"kTPCNSigmaEl2", static_cast<uint8_t>(InputFeatures::kTPCNSigmaEl2)},
      {"kTPCNSigmaPi2", static_cast<uint8_t>(InputFeatures::kTPCNSigmaPi2)},
      {"kTPCNSigmaPr2", static_cast<uint8_t>(InputFeatures::kTPCNSigmaPr2)},
      {"kTOFNSigmaEl2", static_cast<uint8_t>(InputFeatures::kTOFNSigmaEl2)},
      {"kTOFNSigmaPi2", static_cast<uint8_t>(InputFeatures::kTOFNSigmaPi2)},
      {"kTOFNSigmaPr2", static_cast<uint8_t>(InputFeatures::kTOFNSigmaPr2)}};
  }
};

} // namespace o2::analysis

#endif // PWGDQ_CORE_DQMLRESPONSE_H_
