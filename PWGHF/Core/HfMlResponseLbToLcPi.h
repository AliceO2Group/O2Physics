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

#ifndef PWGHF_CORE_HFMLRESPONSELBTOLCPI_H_
#define PWGHF_CORE_HFMLRESPONSELBTOLCPI_H_

#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"

#include "Tools/ML/MlResponse.h"

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_LB(FEATURE)                                       \
  {                                                                \
    #FEATURE, static_cast<uint8_t>(InputFeaturesLbToLcPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_LB_FULL(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesLbToLcPi::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());               \
    break;                                                     \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the GETTER function taking OBJECT in argument
#define CHECK_AND_FILL_VEC_LB_FUNC(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesLbToLcPi::FEATURE): { \
    inputFeatures.emplace_back(GETTER(OBJECT));                \
    break;                                                     \
  }

// Specific case of (OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_LB(GETTER)                         \
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
          CHECK_AND_FILL_VEC_LB(ptProng0);
          CHECK_AND_FILL_VEC_LB(ptProng1);
          CHECK_AND_FILL_VEC_LB(impactParameter0);
          CHECK_AND_FILL_VEC_LB(impactParameter1);
          CHECK_AND_FILL_VEC_LB(impactParameterProduct);
          CHECK_AND_FILL_VEC_LB(chi2PCA);
          CHECK_AND_FILL_VEC_LB(decayLength);
          CHECK_AND_FILL_VEC_LB(decayLengthXY);
          CHECK_AND_FILL_VEC_LB(decayLengthNormalised);
          CHECK_AND_FILL_VEC_LB(decayLengthXYNormalised);
          CHECK_AND_FILL_VEC_LB(cpa);
          CHECK_AND_FILL_VEC_LB(cpaXY);
          CHECK_AND_FILL_VEC_LB(maxNormalisedDeltaIP);
          CHECK_AND_FILL_VEC_LB(prong0MlScoreBkg);
          CHECK_AND_FILL_VEC_LB(prong0MlScorePrompt);
          CHECK_AND_FILL_VEC_LB(prong0MlScoreNonprompt);
          // TPC PID variable
          CHECK_AND_FILL_VEC_LB_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
          // TOF PID variable
          CHECK_AND_FILL_VEC_LB_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
          // Combined PID variables
          CHECK_AND_FILL_VEC_LB_FUNC(prong1, tpcTofNSigmaPi1, o2::pid_tpc_tof_utils::getTpcTofNSigmaPi1);
        }
      } else {
        switch (idx) {
          CHECK_AND_FILL_VEC_LB(ptProng0);
          CHECK_AND_FILL_VEC_LB(ptProng1);
          CHECK_AND_FILL_VEC_LB(impactParameter0);
          CHECK_AND_FILL_VEC_LB(impactParameter1);
          CHECK_AND_FILL_VEC_LB(impactParameterProduct);
          CHECK_AND_FILL_VEC_LB(chi2PCA);
          CHECK_AND_FILL_VEC_LB(decayLength);
          CHECK_AND_FILL_VEC_LB(decayLengthXY);
          CHECK_AND_FILL_VEC_LB(decayLengthNormalised);
          CHECK_AND_FILL_VEC_LB(decayLengthXYNormalised);
          CHECK_AND_FILL_VEC_LB(cpa);
          CHECK_AND_FILL_VEC_LB(cpaXY);
          CHECK_AND_FILL_VEC_LB(maxNormalisedDeltaIP);
          // TPC PID variable
          CHECK_AND_FILL_VEC_LB_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
          // TOF PID variable
          CHECK_AND_FILL_VEC_LB_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
          // Combined PID variables
          CHECK_AND_FILL_VEC_LB_FUNC(prong1, tpcTofNSigmaPi1, o2::pid_tpc_tof_utils::getTpcTofNSigmaPi1);
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
      FILL_MAP_LB(ptProng0),
      FILL_MAP_LB(ptProng1),
      FILL_MAP_LB(impactParameter0),
      FILL_MAP_LB(impactParameter1),
      FILL_MAP_LB(impactParameterProduct),
      FILL_MAP_LB(chi2PCA),
      FILL_MAP_LB(decayLength),
      FILL_MAP_LB(decayLengthXY),
      FILL_MAP_LB(decayLengthNormalised),
      FILL_MAP_LB(decayLengthXYNormalised),
      FILL_MAP_LB(cpa),
      FILL_MAP_LB(cpaXY),
      FILL_MAP_LB(maxNormalisedDeltaIP),
      FILL_MAP_LB(prong0MlScoreBkg),
      FILL_MAP_LB(prong0MlScorePrompt),
      FILL_MAP_LB(prong0MlScoreNonprompt),
      // TPC PID variable
      FILL_MAP_LB(tpcNSigmaPi1),
      // TOF PID variable
      FILL_MAP_LB(tofNSigmaPi1),
      // Combined PID variable
      FILL_MAP_LB(tpcTofNSigmaPi1)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_LB
#undef CHECK_AND_FILL_VEC_LB_FULL
#undef CHECK_AND_FILL_VEC_LB_FUNC
#undef CHECK_AND_FILL_VEC_LB

#endif // PWGHF_CORE_HFMLRESPONSELBTOLCPI_H_
