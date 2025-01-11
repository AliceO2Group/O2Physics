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

/// \file HfMlResponseBplusToD0Pi.h
/// \brief Class to compute the ML response for B± → D0(bar) π± analysis selections
/// \author Antonio Palasciano <antonio.palasciano@ba.infn.it>, INFN Bari

#ifndef PWGHF_CORE_HFMLRESPONSEBPLUSTOD0PI_H_
#define PWGHF_CORE_HFMLRESPONSEBPLUSTOD0PI_H_

#include <map>
#include <string>
#include <vector>

#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_BPLUS(FEATURE)                                       \
  {                                                                   \
    #FEATURE, static_cast<uint8_t>(InputFeaturesBplusToD0Pi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_BPLUS_FULL(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesBplusToD0Pi::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());                  \
    break;                                                        \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the GETTER function taking OBJECT in argument
#define CHECK_AND_FILL_VEC_BPLUS_FUNC(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesBplusToD0Pi::FEATURE): { \
    inputFeatures.emplace_back(GETTER(OBJECT));                   \
    break;                                                        \
  }

// Specific case of CHECK_AND_FILL_VEC_BPLUS_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_BPLUS(GETTER)                         \
  case static_cast<uint8_t>(InputFeaturesBplusToD0Pi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());              \
    break;                                                       \
  }

// where OBJECT is named candidateD , FEATURE = GETTER and INDEX is the index of the vector
#define CHECK_AND_FILL_VEC_D0_INDEX(FEATURE, GETTER1, GETTER2, INDEX) \
  case static_cast<uint8_t>(InputFeaturesBplusToD0Pi::FEATURE): {     \
    if (pdgCode == o2::constants::physics::kD0) {                     \
      inputFeatures.emplace_back((candidateD0.GETTER1())[INDEX]);     \
    } else {                                                          \
      inputFeatures.emplace_back((candidateD0.GETTER2())[INDEX]);     \
    }                                                                 \
    break;                                                            \
  }

namespace o2::analysis
{

enum class InputFeaturesBplusToD0Pi : uint8_t {
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
  prong0MlProbBkg,
  prong0MlProbPrompt,
  prong0MlProbNonPrompt,
  tpcNSigmaPi1,
  tofNSigmaPi1,
  tpcTofNSigmaPi1
};

template <typename TypeOutputScore = float>
class HfMlResponseBplusToD0Pi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseBplusToD0Pi() = default;
  /// Default destructor
  virtual ~HfMlResponseBplusToD0Pi() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the B+ candidate
  /// \param prong1 is the candidate's prong1
  /// \return inputFeatures vector
  template <bool withDmesMl, typename T1, typename T2, typename T3>
  std::vector<float> getInputFeatures(T1 const& candidate,
                                      T2 const& candidateD0,
                                      int const& pdgCode,
                                      T3 const& prong1)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      if constexpr (withDmesMl) {
        switch (idx) {
          CHECK_AND_FILL_VEC_BPLUS(ptProng0);
          CHECK_AND_FILL_VEC_BPLUS(ptProng1);
          CHECK_AND_FILL_VEC_BPLUS(impactParameter0);
          CHECK_AND_FILL_VEC_BPLUS(impactParameter1);
          CHECK_AND_FILL_VEC_BPLUS(impactParameterProduct);
          CHECK_AND_FILL_VEC_BPLUS(chi2PCA);
          CHECK_AND_FILL_VEC_BPLUS(decayLength);
          CHECK_AND_FILL_VEC_BPLUS(decayLengthXY);
          CHECK_AND_FILL_VEC_BPLUS(decayLengthNormalised);
          CHECK_AND_FILL_VEC_BPLUS(decayLengthXYNormalised);
          CHECK_AND_FILL_VEC_BPLUS(cpa);
          CHECK_AND_FILL_VEC_BPLUS(cpaXY);
          CHECK_AND_FILL_VEC_BPLUS(maxNormalisedDeltaIP);
          CHECK_AND_FILL_VEC_D0_INDEX(prong0MlProbBkg, mlProbD0, mlProbD0bar, 0);
          CHECK_AND_FILL_VEC_D0_INDEX(prong0MlProbPrompt, mlProbD0, mlProbD0bar, 1);
          CHECK_AND_FILL_VEC_D0_INDEX(prong0MlProbNonPrompt, mlProbD0, mlProbD0bar, 2);
          // TPC PID variable
          CHECK_AND_FILL_VEC_BPLUS_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
          // TOF PID variable
          CHECK_AND_FILL_VEC_BPLUS_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
          // Combined PID variables
          CHECK_AND_FILL_VEC_BPLUS_FUNC(prong1, tpcTofNSigmaPi1, o2::pid_tpc_tof_utils::getTpcTofNSigmaPi1);
        }
      } else {
        switch (idx) {
          CHECK_AND_FILL_VEC_BPLUS(ptProng0);
          CHECK_AND_FILL_VEC_BPLUS(ptProng1);
          CHECK_AND_FILL_VEC_BPLUS(impactParameter0);
          CHECK_AND_FILL_VEC_BPLUS(impactParameter1);
          CHECK_AND_FILL_VEC_BPLUS(impactParameterProduct);
          CHECK_AND_FILL_VEC_BPLUS(chi2PCA);
          CHECK_AND_FILL_VEC_BPLUS(decayLength);
          CHECK_AND_FILL_VEC_BPLUS(decayLengthXY);
          CHECK_AND_FILL_VEC_BPLUS(decayLengthNormalised);
          CHECK_AND_FILL_VEC_BPLUS(decayLengthXYNormalised);
          CHECK_AND_FILL_VEC_BPLUS(cpa);
          CHECK_AND_FILL_VEC_BPLUS(cpaXY);
          CHECK_AND_FILL_VEC_BPLUS(maxNormalisedDeltaIP);
          // TPC PID variable
          CHECK_AND_FILL_VEC_BPLUS_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
          // TOF PID variable
          CHECK_AND_FILL_VEC_BPLUS_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
          // Combined PID variables
          CHECK_AND_FILL_VEC_BPLUS_FUNC(prong1, tpcTofNSigmaPi1, o2::pid_tpc_tof_utils::getTpcTofNSigmaPi1);
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
      FILL_MAP_BPLUS(ptProng0),
      FILL_MAP_BPLUS(ptProng1),
      FILL_MAP_BPLUS(impactParameter0),
      FILL_MAP_BPLUS(impactParameter1),
      FILL_MAP_BPLUS(impactParameterProduct),
      FILL_MAP_BPLUS(chi2PCA),
      FILL_MAP_BPLUS(decayLength),
      FILL_MAP_BPLUS(decayLengthXY),
      FILL_MAP_BPLUS(decayLengthNormalised),
      FILL_MAP_BPLUS(decayLengthXYNormalised),
      FILL_MAP_BPLUS(cpa),
      FILL_MAP_BPLUS(cpaXY),
      FILL_MAP_BPLUS(maxNormalisedDeltaIP),
      FILL_MAP_BPLUS(prong0MlProbBkg),
      FILL_MAP_BPLUS(prong0MlProbPrompt),
      FILL_MAP_BPLUS(prong0MlProbNonPrompt),
      // TPC PID variable
      FILL_MAP_BPLUS(tpcNSigmaPi1),
      // TOF PID variable
      FILL_MAP_BPLUS(tofNSigmaPi1),
      // Combined PID variable
      FILL_MAP_BPLUS(tpcTofNSigmaPi1)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_BPLUS
#undef CHECK_AND_FILL_VEC_BPLUS_FULL
#undef CHECK_AND_FILL_VEC_BPLUS_FUNC
#undef CHECK_AND_FILL_VEC_BPLUS

#endif // PWGHF_CORE_HFMLRESPONSEBPLUSTOD0PI_H_
