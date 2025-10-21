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

/// \file HfMlResponseBsToDsPi.h
/// \brief Class to compute the ML response for Bs → Ds∓ π± analysis selections
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN

#ifndef PWGHF_CORE_HFMLRESPONSEBSTODSPI_H_
#define PWGHF_CORE_HFMLRESPONSEBSTODSPI_H_

#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"

#include "Tools/ML/MlResponse.h"

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_BS(FEATURE)                                       \
  {                                                                \
    #FEATURE, static_cast<uint8_t>(InputFeaturesBsToDsPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_BS_FULL(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesBsToDsPi::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());               \
    break;                                                     \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the GETTER function taking OBJECT in argument
#define CHECK_AND_FILL_VEC_BS_FUNC(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesBsToDsPi::FEATURE): { \
    inputFeatures.emplace_back(GETTER(OBJECT));                \
    break;                                                     \
  }

// Specific case of CHECK_AND_FILL_VEC_BS_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_BS(GETTER)                         \
  case static_cast<uint8_t>(InputFeaturesBsToDsPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());           \
    break;                                                    \
  }

namespace o2::analysis
{

enum class InputFeaturesBsToDsPi : uint8_t {
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
class HfMlResponseBsToDsPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseBsToDsPi() = default;
  /// Default destructor
  virtual ~HfMlResponseBsToDsPi() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Bs candidate
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
          CHECK_AND_FILL_VEC_BS(ptProng0);
          CHECK_AND_FILL_VEC_BS(ptProng1);
          CHECK_AND_FILL_VEC_BS(impactParameter0);
          CHECK_AND_FILL_VEC_BS(impactParameter1);
          CHECK_AND_FILL_VEC_BS(impactParameterProduct);
          CHECK_AND_FILL_VEC_BS(chi2PCA);
          CHECK_AND_FILL_VEC_BS(decayLength);
          CHECK_AND_FILL_VEC_BS(decayLengthXY);
          CHECK_AND_FILL_VEC_BS(decayLengthNormalised);
          CHECK_AND_FILL_VEC_BS(decayLengthXYNormalised);
          CHECK_AND_FILL_VEC_BS(cpa);
          CHECK_AND_FILL_VEC_BS(cpaXY);
          CHECK_AND_FILL_VEC_BS(maxNormalisedDeltaIP);
          CHECK_AND_FILL_VEC_BS(prong0MlScoreBkg);
          CHECK_AND_FILL_VEC_BS(prong0MlScorePrompt);
          CHECK_AND_FILL_VEC_BS(prong0MlScoreNonprompt);
          //  Pion PID variables
          CHECK_AND_FILL_VEC_BS_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
          CHECK_AND_FILL_VEC_BS_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
          CHECK_AND_FILL_VEC_BS_FUNC(prong1, tpcTofNSigmaPi1, o2::pid_tpc_tof_utils::getTpcTofNSigmaPi1);
        }
      } else {
        switch (idx) {
          CHECK_AND_FILL_VEC_BS(ptProng0);
          CHECK_AND_FILL_VEC_BS(ptProng1);
          CHECK_AND_FILL_VEC_BS(impactParameter0);
          CHECK_AND_FILL_VEC_BS(impactParameter1);
          CHECK_AND_FILL_VEC_BS(impactParameterProduct);
          CHECK_AND_FILL_VEC_BS(chi2PCA);
          CHECK_AND_FILL_VEC_BS(decayLength);
          CHECK_AND_FILL_VEC_BS(decayLengthXY);
          CHECK_AND_FILL_VEC_BS(decayLengthNormalised);
          CHECK_AND_FILL_VEC_BS(decayLengthXYNormalised);
          CHECK_AND_FILL_VEC_BS(cpa);
          CHECK_AND_FILL_VEC_BS(cpaXY);
          CHECK_AND_FILL_VEC_BS(maxNormalisedDeltaIP);
          // Pion PID variables
          CHECK_AND_FILL_VEC_BS_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
          CHECK_AND_FILL_VEC_BS_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
          CHECK_AND_FILL_VEC_BS_FUNC(prong1, tpcTofNSigmaPi1, o2::pid_tpc_tof_utils::getTpcTofNSigmaPi1);
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
      FILL_MAP_BS(ptProng0),
      FILL_MAP_BS(ptProng1),
      FILL_MAP_BS(impactParameter0),
      FILL_MAP_BS(impactParameter1),
      FILL_MAP_BS(impactParameterProduct),
      FILL_MAP_BS(chi2PCA),
      FILL_MAP_BS(decayLength),
      FILL_MAP_BS(decayLengthXY),
      FILL_MAP_BS(decayLengthNormalised),
      FILL_MAP_BS(decayLengthXYNormalised),
      FILL_MAP_BS(cpa),
      FILL_MAP_BS(cpaXY),
      FILL_MAP_BS(maxNormalisedDeltaIP),
      FILL_MAP_BS(prong0MlScoreBkg),
      FILL_MAP_BS(prong0MlScorePrompt),
      FILL_MAP_BS(prong0MlScoreNonprompt),
      // Pion PID variables
      FILL_MAP_BS(tpcNSigmaPi1),
      FILL_MAP_BS(tofNSigmaPi1),
      FILL_MAP_BS(tpcTofNSigmaPi1)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_BS
#undef CHECK_AND_FILL_VEC_BS_FULL
#undef CHECK_AND_FILL_VEC_BS_FUNC
#undef CHECK_AND_FILL_VEC_BS

#endif // PWGHF_CORE_HFMLRESPONSEBSTODSPI_H_
