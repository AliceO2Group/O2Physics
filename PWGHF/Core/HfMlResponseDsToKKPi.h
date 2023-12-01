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

/// \file HfMlResponseDsToKKPi.h
/// \brief Class to compute the ML response for Ds∓ → K∓ K∓ π± analysis selections
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>

#ifndef PWGHF_CORE_HFMLRESPONSEDSTOKKPI_H_
#define PWGHF_CORE_HFMLRESPONSEDSTOKKPI_H_

#include <vector>

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponse.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_DS(FEATURE)                                         \
  {                                                                  \
#FEATURE, static_cast < uint8_t>(InputFeaturesDsToKKPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_DS_FULL(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesDsToKKPi::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());               \
    break;                                                     \
  }

// Specific case of CHECK_AND_FILL_VEC_DS_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_DS(GETTER)                         \
  case static_cast<uint8_t>(InputFeaturesDsToKKPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());           \
    break;                                                    \
  }

// Variation of CHECK_AND_FILL_VEC_DS_FULL(OBJECT, FEATURE, GETTER)
// where GETTER is a method of hfHelper
#define CHECK_AND_FILL_VEC_DS_HFHELPER(OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesDsToKKPi::FEATURE): {  \
    inputFeatures.emplace_back(hfHelper.GETTER(OBJECT));        \
    break;                                                      \
  }

// Variation of CHECK_AND_FILL_VEC_DS_HFHELPER(OBJECT, FEATURE, GETTER)
// where GETTER1 and GETTER2 are methods of hfHelper, and the variable
// is filled depending on whether it is a DsToKKPi or a DsToPiKK
#define CHECK_AND_FILL_VEC_DS_HFHELPER_SIGNED(OBJECT, FEATURE, GETTER1, GETTER2) \
  case static_cast<uint8_t>(InputFeaturesDsToKKPi::FEATURE): {                   \
    if (caseDsToKKPi) {                                                          \
      inputFeatures.emplace_back(hfHelper.GETTER1(OBJECT));                      \
    } else {                                                                     \
      inputFeatures.emplace_back(hfHelper.GETTER2(OBJECT));                      \
    }                                                                            \
    break;                                                                       \
  }

namespace o2::analysis
{
enum class InputFeaturesDsToKKPi : uint8_t {
  chi2PCA = 0,
  decayLength,
  decayLengthXY,
  decayLengthNormalised,
  decayLengthXYNormalised,
  maxNormalisedDeltaIP,
  cpa,
  cpaXY,
  ptProng0,
  ptProng1,
  ptProng2,
  impactParameterXY0,
  impactParameterXY1,
  impactParameterXY2,
  nSigTpcTofPi0,
  nSigTpcTofPi1,
  nSigTpcTofPi2,
  nSigTpcTofKa0,
  nSigTpcTofKa1,
  nSigTpcTofKa2,
  cos3PiK,
  deltaMassPhi
};

template <typename TypeOutputScore = float>
class HfMlResponseDsToKKPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseDsToKKPi() = default;
  /// Default destructor
  virtual ~HfMlResponseDsToKKPi() = default;

  HfHelper hfHelper;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Ds candidate
  /// \param prong0 is the candidate's prong0
  /// \param prong1 is the candidate's prong1
  /// \param prong2 is the candidate's prong2
  /// \param caseDsToKKPi used to divide the case DsToKKPi from DsToPiKK
  /// \return inputFeatures vector
  template <typename T1, typename T2>
  std::vector<float> getInputFeatures(T1 const& candidate,
                                      T2 const& prong0, T2 const& prong1, T2 const& prong2, bool const& caseDsToKKPi)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_DS(chi2PCA);
        CHECK_AND_FILL_VEC_DS(decayLength);
        CHECK_AND_FILL_VEC_DS(decayLengthXY);
        CHECK_AND_FILL_VEC_DS(decayLengthNormalised);
        CHECK_AND_FILL_VEC_DS(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_DS(maxNormalisedDeltaIP);
        CHECK_AND_FILL_VEC_DS(cpa);
        CHECK_AND_FILL_VEC_DS(cpaXY);
        CHECK_AND_FILL_VEC_DS(ptProng0);
        CHECK_AND_FILL_VEC_DS(ptProng1);
        CHECK_AND_FILL_VEC_DS(ptProng2);
        CHECK_AND_FILL_VEC_DS_FULL(candidate, impactParameterXY0, impactParameter0);
        CHECK_AND_FILL_VEC_DS_FULL(candidate, impactParameterXY1, impactParameter1);
        CHECK_AND_FILL_VEC_DS_FULL(candidate, impactParameterXY2, impactParameter2);
        // Combined PID variables
        CHECK_AND_FILL_VEC_DS_FULL(prong0, nSigTpcTofPi0, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_DS_FULL(prong1, nSigTpcTofPi1, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_DS_FULL(prong2, nSigTpcTofPi2, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_DS_FULL(prong0, nSigTpcTofKa0, tpcTofNSigmaKa);
        CHECK_AND_FILL_VEC_DS_FULL(prong1, nSigTpcTofKa1, tpcTofNSigmaKa);
        CHECK_AND_FILL_VEC_DS_FULL(prong2, nSigTpcTofKa2, tpcTofNSigmaKa);

        // Ds specific variables
        CHECK_AND_FILL_VEC_DS_HFHELPER_SIGNED(candidate, cos3PiK, cos3PiKDsToKKPi, cos3PiKDsToPiKK);
        CHECK_AND_FILL_VEC_DS_HFHELPER_SIGNED(candidate, deltaMassPhi, deltaMassPhiDsToKKPi, deltaMassPhiDsToPiKK);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_DS(chi2PCA),
      FILL_MAP_DS(decayLength),
      FILL_MAP_DS(decayLengthXY),
      FILL_MAP_DS(decayLengthNormalised),
      FILL_MAP_DS(decayLengthXYNormalised),
      FILL_MAP_DS(maxNormalisedDeltaIP),
      FILL_MAP_DS(cpa),
      FILL_MAP_DS(cpaXY),
      FILL_MAP_DS(ptProng0),
      FILL_MAP_DS(ptProng1),
      FILL_MAP_DS(ptProng2),
      FILL_MAP_DS(impactParameterXY0),
      FILL_MAP_DS(impactParameterXY1),
      FILL_MAP_DS(impactParameterXY2),
      // Combined PID variables
      FILL_MAP_DS(nSigTpcTofPi0),
      FILL_MAP_DS(nSigTpcTofPi1),
      FILL_MAP_DS(nSigTpcTofPi2),
      FILL_MAP_DS(nSigTpcTofKa0),
      FILL_MAP_DS(nSigTpcTofKa1),
      FILL_MAP_DS(nSigTpcTofKa2),

      // Ds specific variables
      FILL_MAP_DS(cos3PiK),
      FILL_MAP_DS(deltaMassPhi)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_DS
#undef CHECK_AND_FILL_VEC_DS_FULL
#undef CHECK_AND_FILL_VEC_DS
#undef CHECK_AND_FILL_VEC_DS_HFHELPER
#undef CHECK_AND_FILL_VEC_DS_HFHELPER_SIGNED

#endif // PWGHF_CORE_HFMLRESPONSEDSTOKKPI_H_
