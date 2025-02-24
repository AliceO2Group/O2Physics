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

/// \file HfMlResponseB0ToDPi.h
/// \brief Class to compute the ML response for B0 → D∓ π± analysis selections
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#ifndef PWGHF_CORE_HFMLRESPONSEB0TODPI_H_
#define PWGHF_CORE_HFMLRESPONSEB0TODPI_H_

#include <map>
#include <string>
#include <vector>

#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_B0(FEATURE)                                      \
  {                                                               \
    #FEATURE, static_cast<uint8_t>(InputFeaturesB0ToDPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_B0_FULL(OBJECT, FEATURE, GETTER)   \
  case static_cast<uint8_t>(InputFeaturesB0ToDPi::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());              \
    break;                                                    \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the GETTER function taking OBJECT in argument
#define CHECK_AND_FILL_VEC_B0_FUNC(OBJECT, FEATURE, GETTER)   \
  case static_cast<uint8_t>(InputFeaturesB0ToDPi::FEATURE): { \
    inputFeatures.emplace_back(GETTER(OBJECT));               \
    break;                                                    \
  }

// Specific case of CHECK_AND_FILL_VEC_B0_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_B0(GETTER)                        \
  case static_cast<uint8_t>(InputFeaturesB0ToDPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());          \
    break;                                                   \
  }

namespace o2::analysis
{

enum class InputFeaturesB0ToDPi : uint8_t {
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
class HfMlResponseB0ToDPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseB0ToDPi() = default;
  /// Default destructor
  virtual ~HfMlResponseB0ToDPi() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the B0 candidate
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
          CHECK_AND_FILL_VEC_B0(ptProng0);
          CHECK_AND_FILL_VEC_B0(ptProng1);
          CHECK_AND_FILL_VEC_B0(impactParameter0);
          CHECK_AND_FILL_VEC_B0(impactParameter1);
          CHECK_AND_FILL_VEC_B0(impactParameterProduct);
          CHECK_AND_FILL_VEC_B0(chi2PCA);
          CHECK_AND_FILL_VEC_B0(decayLength);
          CHECK_AND_FILL_VEC_B0(decayLengthXY);
          CHECK_AND_FILL_VEC_B0(decayLengthNormalised);
          CHECK_AND_FILL_VEC_B0(decayLengthXYNormalised);
          CHECK_AND_FILL_VEC_B0(cpa);
          CHECK_AND_FILL_VEC_B0(cpaXY);
          CHECK_AND_FILL_VEC_B0(maxNormalisedDeltaIP);
          CHECK_AND_FILL_VEC_B0(prong0MlScoreBkg);
          CHECK_AND_FILL_VEC_B0(prong0MlScorePrompt);
          CHECK_AND_FILL_VEC_B0(prong0MlScoreNonprompt);
          // TPC PID variable
          CHECK_AND_FILL_VEC_B0_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
          // TOF PID variable
          CHECK_AND_FILL_VEC_B0_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
          // Combined PID variables
          CHECK_AND_FILL_VEC_B0_FUNC(prong1, tpcTofNSigmaPi1, o2::pid_tpc_tof_utils::getTpcTofNSigmaPi1);
        }
      } else {
        switch (idx) {
          CHECK_AND_FILL_VEC_B0(ptProng0);
          CHECK_AND_FILL_VEC_B0(ptProng1);
          CHECK_AND_FILL_VEC_B0(impactParameter0);
          CHECK_AND_FILL_VEC_B0(impactParameter1);
          CHECK_AND_FILL_VEC_B0(impactParameterProduct);
          CHECK_AND_FILL_VEC_B0(chi2PCA);
          CHECK_AND_FILL_VEC_B0(decayLength);
          CHECK_AND_FILL_VEC_B0(decayLengthXY);
          CHECK_AND_FILL_VEC_B0(decayLengthNormalised);
          CHECK_AND_FILL_VEC_B0(decayLengthXYNormalised);
          CHECK_AND_FILL_VEC_B0(cpa);
          CHECK_AND_FILL_VEC_B0(cpaXY);
          CHECK_AND_FILL_VEC_B0(maxNormalisedDeltaIP);
          // TPC PID variable
          CHECK_AND_FILL_VEC_B0_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
          // TOF PID variable
          CHECK_AND_FILL_VEC_B0_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
          // Combined PID variables
          CHECK_AND_FILL_VEC_B0_FUNC(prong1, tpcTofNSigmaPi1, o2::pid_tpc_tof_utils::getTpcTofNSigmaPi1);
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
      FILL_MAP_B0(ptProng0),
      FILL_MAP_B0(ptProng1),
      FILL_MAP_B0(impactParameter0),
      FILL_MAP_B0(impactParameter1),
      FILL_MAP_B0(impactParameterProduct),
      FILL_MAP_B0(chi2PCA),
      FILL_MAP_B0(decayLength),
      FILL_MAP_B0(decayLengthXY),
      FILL_MAP_B0(decayLengthNormalised),
      FILL_MAP_B0(decayLengthXYNormalised),
      FILL_MAP_B0(cpa),
      FILL_MAP_B0(cpaXY),
      FILL_MAP_B0(maxNormalisedDeltaIP),
      FILL_MAP_B0(prong0MlScoreBkg),
      FILL_MAP_B0(prong0MlScorePrompt),
      FILL_MAP_B0(prong0MlScoreNonprompt),
      // TPC PID variable
      FILL_MAP_B0(tpcNSigmaPi1),
      // TOF PID variable
      FILL_MAP_B0(tofNSigmaPi1),
      // Combined PID variable
      FILL_MAP_B0(tpcTofNSigmaPi1)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_B0
#undef CHECK_AND_FILL_VEC_B0_FULL
#undef CHECK_AND_FILL_VEC_B0_FUNC
#undef CHECK_AND_FILL_VEC_B0

#endif // PWGHF_CORE_HFMLRESPONSEB0TODPI_H_
