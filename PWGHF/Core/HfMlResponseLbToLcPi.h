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

/// \file HfMlResponseLbToLcPi.h
/// \brief Class to compute the ML response for Lb → Lc∓ π± analysis selections
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg University

#ifndef PWGHF_CORE_HfMlResponseLBTOLCPI_H_
#define PWGHF_CORE_HfMlResponseLBTOLCPI_H_

#include <map>
#include <string>
#include <vector>

#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_Lb(FEATURE)                                       \
  {                                                                \
    #FEATURE, static_cast<uint8_t>(InputFeaturesLbToLcPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_Lb_FULL(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesLbToLcPi::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());               \
    break;                                                     \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the GETTER function taking OBJECT in argument
#define CHECK_AND_FILL_VEC_Lb_FUNC(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesLbToLcPi::FEATURE): { \
    inputFeatures.emplace_back(GETTER(OBJECT));                \
    break;                                                     \
  }

// Specific case of CHECK_AND_FILL_VEC_Lb_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_Lb(GETTER)                         \
  case static_cast<uint8_t>(InputFeaturesLbToLcPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());           \
    break;                                                    \
  }

namespace o2::analysis
{

enum class InputFeaturesLbToLcPi : uint8_t {
  ptProng0 = 0,
  ptProng1,
  impactParameter0,
  impactParameter1,
  impactParameterProduct,
  chi2PCA,
  decayLength,
  decayLengthXY,
  decayLengthNormalised,
  decayLengthXYNormalised,
  cpa,
  cpaXY,
  maxNormalisedDeltaIP,
  prong0MlScoreBkg,
  prong0MlScorePrompt,
  prong0MlScoreNonprompt,
  tpcNSigmaPi1,
  tofNSigmaPi1,
  tpcTofNSigmaPi1
};

template <typename TypeOutputScore = float>
class HfMlResponseLbToLcPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseLbToLcPi() = default;
  /// Default destructor
  virtual ~HfMlResponseLbToLcPi() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Lb candidate
  /// \param prong1 is the candidate's prong1
  /// \return inputFeatures vector
  template <bool withDmesMl, typename T1, typename T2>
  std::vector<float> getInputFeatures(T1 const& candidate,
                                      T2 const& prong1)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      if constexpr (withDmesMl) {
        switch (idx) {
          CHECK_AND_FILL_VEC_Lb(ptProng0);
          CHECK_AND_FILL_VEC_Lb(ptProng1);
          CHECK_AND_FILL_VEC_Lb(impactParameter0);
          CHECK_AND_FILL_VEC_Lb(impactParameter1);
          CHECK_AND_FILL_VEC_Lb(impactParameterProduct);
          CHECK_AND_FILL_VEC_Lb(chi2PCA);
          CHECK_AND_FILL_VEC_Lb(decayLength);
          CHECK_AND_FILL_VEC_Lb(decayLengthXY);
          CHECK_AND_FILL_VEC_Lb(decayLengthNormalised);
          CHECK_AND_FILL_VEC_Lb(decayLengthXYNormalised);
          CHECK_AND_FILL_VEC_Lb(cpa);
          CHECK_AND_FILL_VEC_Lb(cpaXY);
          CHECK_AND_FILL_VEC_Lb(maxNormalisedDeltaIP);
          CHECK_AND_FILL_VEC_Lb(prong0MlScoreBkg);
          CHECK_AND_FILL_VEC_Lb(prong0MlScorePrompt);
          CHECK_AND_FILL_VEC_Lb(prong0MlScoreNonprompt);
          // TPC PID variable
          CHECK_AND_FILL_VEC_Lb_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
          // TOF PID variable
          CHECK_AND_FILL_VEC_Lb_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
          // Combined PID variables
          CHECK_AND_FILL_VEC_Lb_FUNC(prong1, tpcTofNSigmaPi1, o2::pid_tpc_tof_utils::getTpcTofNSigmaPi1);
        }
      } else {
        switch (idx) {
          CHECK_AND_FILL_VEC_Lb(ptProng0);
          CHECK_AND_FILL_VEC_Lb(ptProng1);
          CHECK_AND_FILL_VEC_Lb(impactParameter0);
          CHECK_AND_FILL_VEC_Lb(impactParameter1);
          CHECK_AND_FILL_VEC_Lb(impactParameterProduct);
          CHECK_AND_FILL_VEC_Lb(chi2PCA);
          CHECK_AND_FILL_VEC_Lb(decayLength);
          CHECK_AND_FILL_VEC_Lb(decayLengthXY);
          CHECK_AND_FILL_VEC_Lb(decayLengthNormalised);
          CHECK_AND_FILL_VEC_Lb(decayLengthXYNormalised);
          CHECK_AND_FILL_VEC_Lb(cpa);
          CHECK_AND_FILL_VEC_Lb(cpaXY);
          CHECK_AND_FILL_VEC_Lb(maxNormalisedDeltaIP);
          // TPC PID variable
          CHECK_AND_FILL_VEC_Lb_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
          // TOF PID variable
          CHECK_AND_FILL_VEC_Lb_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
          // Combined PID variables
          CHECK_AND_FILL_VEC_Lb_FUNC(prong1, tpcTofNSigmaPi1, o2::pid_tpc_tof_utils::getTpcTofNSigmaPi1);
        }
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_Lb(ptProng0),
      FILL_MAP_Lb(ptProng1),
      FILL_MAP_Lb(impactParameter0),
      FILL_MAP_Lb(impactParameter1),
      FILL_MAP_Lb(impactParameterProduct),
      FILL_MAP_Lb(chi2PCA),
      FILL_MAP_Lb(decayLength),
      FILL_MAP_Lb(decayLengthXY),
      FILL_MAP_Lb(decayLengthNormalised),
      FILL_MAP_Lb(decayLengthXYNormalised),
      FILL_MAP_Lb(cpa),
      FILL_MAP_Lb(cpaXY),
      FILL_MAP_Lb(maxNormalisedDeltaIP),
      FILL_MAP_Lb(prong0MlScoreBkg),
      FILL_MAP_Lb(prong0MlScorePrompt),
      FILL_MAP_Lb(prong0MlScoreNonprompt),
      // TPC PID variable
      FILL_MAP_Lb(tpcNSigmaPi1),
      // TOF PID variable
      FILL_MAP_Lb(tofNSigmaPi1),
      // Combined PID variable
      FILL_MAP_Lb(tpcTofNSigmaPi1)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_Lb
#undef CHECK_AND_FILL_VEC_Lb_FULL
#undef CHECK_AND_FILL_VEC_Lb_FUNC
#undef CHECK_AND_FILL_VEC_Lb

#endif // PWGHF_CORE_HfMlResponseLBTOLCPI_H_
