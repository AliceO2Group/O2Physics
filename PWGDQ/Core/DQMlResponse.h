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

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP(FEATURE) \
  {                       \
    #FEATURE, static_cast<uint8_t>(InputFeatures::FEATURE)}

namespace o2::analysis
{

enum class InputFeatures : uint8_t { // refer to DielectronsAll
  fMass = 0,
  fPt,
  fEta,
  fPhi,
  fPt1,
  fITSChi2NCl1,
  fTPCNClsCR1,
  fTPCNClsFound1,
  fTPCChi2NCl1,
  fDcaXY1,
  fDcaZ1,
  fTPCNSigmaEl1,
  fTPCNSigmaPi1,
  fTPCNSigmaPr1,
  fTOFNSigmaEl1,
  fTOFNSigmaPi1,
  fTOFNSigmaPr1,
  fPt2,
  fITSChi2NCl2,
  fTPCNClsCR2,
  fTPCNClsFound2,
  fTPCChi2NCl2,
  fDcaXY2,
  fDcaZ2,
  fTPCNSigmaEl2,
  fTPCNSigmaPi2,
  fTPCNSigmaPr2,
  fTOFNSigmaEl2,
  fTOFNSigmaPi2,
  fTOFNSigmaPr2,
};

static const std::map<InputFeatures, std::string> gFeatureNameMap = {
  {InputFeatures::fMass, "fMass"},
  {InputFeatures::fPt, "fPt"},
  {InputFeatures::fEta, "fEta"},
  {InputFeatures::fPhi, "fPhi"},
  {InputFeatures::fPt1, "fPt1"},
  {InputFeatures::fITSChi2NCl1, "fITSChi2NCl1"},
  {InputFeatures::fTPCNClsCR1, "fTPCNClsCR1"},
  {InputFeatures::fTPCNClsFound1, "fTPCNClsFound1"},
  {InputFeatures::fTPCChi2NCl1, "fTPCChi2NCl1"},
  {InputFeatures::fDcaXY1, "fDcaXY1"},
  {InputFeatures::fDcaZ1, "fDcaZ1"},
  {InputFeatures::fTPCNSigmaEl1, "fTPCNSigmaEl1"},
  {InputFeatures::fTPCNSigmaPi1, "fTPCNSigmaPi1"},
  {InputFeatures::fTPCNSigmaPr1, "fTPCNSigmaPr1"},
  {InputFeatures::fTOFNSigmaEl1, "fTOFNSigmaEl1"},
  {InputFeatures::fTOFNSigmaPi1, "fTOFNSigmaPi1"},
  {InputFeatures::fTOFNSigmaPr1, "fTOFNSigmaPr1"},
  {InputFeatures::fPt2, "fPt2"},
  {InputFeatures::fITSChi2NCl2, "fITSChi2NCl2"},
  {InputFeatures::fTPCNClsCR2, "fTPCNClsCR2"},
  {InputFeatures::fTPCNClsFound2, "fTPCNClsFound2"},
  {InputFeatures::fTPCChi2NCl2, "fTPCChi2NCl2"},
  {InputFeatures::fDcaXY2, "fDcaXY2"},
  {InputFeatures::fDcaZ2, "fDcaZ2"},
  {InputFeatures::fTPCNSigmaEl2, "fTPCNSigmaEl2"},
  {InputFeatures::fTPCNSigmaPi2, "fTPCNSigmaPi2"},
  {InputFeatures::fTPCNSigmaPr2, "fTPCNSigmaPr2"},
  {InputFeatures::fTOFNSigmaEl2, "fTOFNSigmaEl2"},
  {InputFeatures::fTOFNSigmaPi2, "fTOFNSigmaPi2"},
  {InputFeatures::fTOFNSigmaPr2, "fTOFNSigmaPr2"}};

template <typename TypeOutputScore = float>
class DQMlResponse : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  DQMlResponse() = default;
  /// Default destructor
  virtual ~DQMlResponse() = default;

  /// Method to get the input features vector needed for ML inference
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename TValues>
  std::vector<float> getInputFeatures(const T1& t1,
                                      const T2& t2,
                                      const TValues& fg) const
  {
    using Accessor = std::function<float(const T1&, const T2&, const TValues&)>;
    static const std::unordered_map<std::string, Accessor> featureMap{
      {"fMass", [](auto const&, auto const&, auto const& v) { return v[VarManager::kMass]; }},
      {"fPt", [](auto const&, auto const&, auto const& v) { return v[VarManager::kPt]; }},
      {"fEta", [](auto const&, auto const&, auto const& v) { return v[VarManager::kEta]; }},
      {"fPhi", [](auto const&, auto const&, auto const& v) { return v[VarManager::kPhi]; }},

      {"fPt1", [](auto const& t1, auto const&, auto const&) { return t1.pt(); }},
      {"fITSChi2NCl1", [](auto const& t1, auto const&, auto const&) { return t1.itsChi2NCl(); }},
      {"fTPCNClsCR1", [](auto const& t1, auto const&, auto const&) { return t1.tpcNClsCrossedRows(); }},
      {"fTPCNClsFound1", [](auto const& t1, auto const&, auto const&) { return t1.tpcNClsFound(); }},
      {"fTPCChi2NCl1", [](auto const& t1, auto const&, auto const&) { return t1.tpcChi2NCl(); }},
      {"fDcaXY1", [](auto const& t1, auto const&, auto const&) { return t1.dcaXY(); }},
      {"fDcaZ1", [](auto const& t1, auto const&, auto const&) { return t1.dcaZ(); }},
      {"fTPCNSigmaEl1", [](auto const& t1, auto const&, auto const&) { return t1.tpcNSigmaEl(); }},
      {"fTPCNSigmaPi1", [](auto const& t1, auto const&, auto const&) { return t1.tpcNSigmaPi(); }},
      {"fTPCNSigmaPr1", [](auto const& t1, auto const&, auto const&) { return t1.tpcNSigmaPr(); }},
      {"fTOFNSigmaEl1", [](auto const& t1, auto const&, auto const&) { return t1.tofNSigmaEl(); }},
      {"fTOFNSigmaPi1", [](auto const& t1, auto const&, auto const&) { return t1.tofNSigmaPi(); }},
      {"fTOFNSigmaPr1", [](auto const& t1, auto const&, auto const&) { return t1.tofNSigmaPr(); }},

      {"fPt2", [](auto const&, auto const& t2, auto const&) { return t2.pt(); }},
      {"fITSChi2NCl2", [](auto const&, auto const& t2, auto const&) { return t2.itsChi2NCl(); }},
      {"fTPCNClsCR2", [](auto const&, auto const& t2, auto const&) { return t2.tpcNClsCrossedRows(); }},
      {"fTPCNClsFound2", [](auto const&, auto const& t2, auto const&) { return t2.tpcNClsFound(); }},
      {"fTPCChi2NCl2", [](auto const&, auto const& t2, auto const&) { return t2.tpcChi2NCl(); }},
      {"fDcaXY2", [](auto const&, auto const& t2, auto const&) { return t2.dcaXY(); }},
      {"fDcaZ2", [](auto const&, auto const& t2, auto const&) { return t2.dcaZ(); }},
      {"fTPCNSigmaEl2", [](auto const&, auto const& t2, auto const&) { return t2.tpcNSigmaEl(); }},
      {"fTPCNSigmaPi2", [](auto const&, auto const& t2, auto const&) { return t2.tpcNSigmaPi(); }},
      {"fTPCNSigmaPr2", [](auto const&, auto const& t2, auto const&) { return t2.tpcNSigmaPr(); }},
      {"fTOFNSigmaEl2", [](auto const&, auto const& t2, auto const&) { return t2.tofNSigmaEl(); }},
      {"fTOFNSigmaPi2", [](auto const&, auto const& t2, auto const&) { return t2.tofNSigmaPi(); }},
      {"fTOFNSigmaPr2", [](auto const&, auto const& t2, auto const&) { return t2.tofNSigmaPr(); }}};

    std::vector<float> dqInputFeatures;
    dqInputFeatures.reserve(MlResponse<TypeOutputScore>::mCachedIndices.size());

    for (auto idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      auto enumIdx = static_cast<InputFeatures>(idx);
      const auto& name = gFeatureNameMap.at(enumIdx);

      auto acc = featureMap.find(name);
      if (acc == featureMap.end()) {
        LOG(error) << "Missing accessor for " << name;
        continue;
      } else {
        dqInputFeatures.push_back(acc->second(t1, t2, fg));
      }
    }
    return dqInputFeatures;
  }

 protected:
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP(fMass),
      FILL_MAP(fPt),
      FILL_MAP(fEta),
      FILL_MAP(fPhi),
      FILL_MAP(fPt1),
      FILL_MAP(fITSChi2NCl1),
      FILL_MAP(fTPCNClsCR1),
      FILL_MAP(fTPCNClsFound1),
      FILL_MAP(fTPCChi2NCl1),
      FILL_MAP(fDcaXY1),
      FILL_MAP(fDcaZ1),
      FILL_MAP(fTPCNSigmaEl1),
      FILL_MAP(fTPCNSigmaPi1),
      FILL_MAP(fTPCNSigmaPr1),
      FILL_MAP(fTOFNSigmaEl1),
      FILL_MAP(fTOFNSigmaPi1),
      FILL_MAP(fTOFNSigmaPr1),
      FILL_MAP(fPt2),
      FILL_MAP(fITSChi2NCl2),
      FILL_MAP(fTPCNClsCR2),
      FILL_MAP(fTPCNClsFound2),
      FILL_MAP(fTPCChi2NCl2),
      FILL_MAP(fDcaXY2),
      FILL_MAP(fDcaZ2),
      FILL_MAP(fTPCNSigmaEl2),
      FILL_MAP(fTPCNSigmaPi2),
      FILL_MAP(fTPCNSigmaPr2),
      FILL_MAP(fTOFNSigmaEl2),
      FILL_MAP(fTOFNSigmaPi2),
      FILL_MAP(fTOFNSigmaPr2)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP

#endif // PWGDQ_CORE_DQMLRESPONSE_H_
