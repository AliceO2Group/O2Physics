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

/// \file HfMlResponseD0ToKPi.h
/// \brief Class to compute the ML response for D0 → K∓ π± analysis selections
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>, IJCLab Orsay

#ifndef PWGHF_CORE_HFMLRESPONSED0TOKPI_H_
#define PWGHF_CORE_HFMLRESPONSED0TOKPI_H_

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponse.h"

#include "Tools/ML/MlResponse.h"

#include <CommonConstants/PhysicsConstants.h>

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_D0(FEATURE)                                      \
  {                                                               \
    #FEATURE, static_cast<uint8_t>(InputFeaturesD0ToKPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_D0_FULL(OBJECT, FEATURE, GETTER)   \
  case static_cast<uint8_t>(InputFeaturesD0ToKPi::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());              \
    break;                                                    \
  }

// Specific case of CHECK_AND_FILL_VEC_D0_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_D0(GETTER)                        \
  case static_cast<uint8_t>(InputFeaturesD0ToKPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());          \
    break;                                                   \
  }

// Variation of CHECK_AND_FILL_VEC_D0_FULL(OBJECT, FEATURE, GETTER)
// where GETTER is a method of HfHelper
#define CHECK_AND_FILL_VEC_D0_HFHELPER(OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesD0ToKPi::FEATURE): {   \
    inputFeatures.emplace_back(HfHelper::GETTER(OBJECT));       \
    break;                                                      \
  }

// Variation of CHECK_AND_FILL_VEC_D0_HFHELPER(OBJECT, FEATURE, GETTER)
// where GETTER1 and GETTER2 are methods of HfHelper, and the variable
// is filled depending on whether it is a D0 or a D0bar
#define CHECK_AND_FILL_VEC_D0_HFHELPER_SIGNED(OBJECT, FEATURE, GETTER1, GETTER2) \
  case static_cast<uint8_t>(InputFeaturesD0ToKPi::FEATURE): {                    \
    if (pdgCode == o2::constants::physics::kD0) {                                \
      inputFeatures.emplace_back(HfHelper::GETTER1(OBJECT));                     \
    } else {                                                                     \
      inputFeatures.emplace_back(HfHelper::GETTER2(OBJECT));                     \
    }                                                                            \
    break;                                                                       \
  }

// Variation of CHECK_AND_FILL_VEC_D0_HFHELPER(OBJECT, FEATURE, GETTER)
// where GETTER1 and GETTER2 are methods of HfHelper, and the variable
// is filled depending on whether it is a D0 or a D0bar
#define CHECK_AND_FILL_VEC_D0_OBJECT_HFHELPER_SIGNED(OBJECT1, OBJECT2, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesD0ToKPi::FEATURE): {                           \
    if (pdgCode == o2::constants::physics::kD0) {                                       \
      inputFeatures.emplace_back(OBJECT1.GETTER());                                     \
    } else {                                                                            \
      inputFeatures.emplace_back(OBJECT2.GETTER());                                     \
    }                                                                                   \
    break;                                                                              \
  }

// Variation of CHECK_AND_FILL_VEC_D0_HFHELPER_SIGNED(OBJECT, FEATURE, GETTER1, GETTER2)
// where GETTER1 and GETTER2 are methods of the OBJECT, and the variable
// is filled depending on whether it is a D0 or a D0bar
#define CHECK_AND_FILL_VEC_D0_SIGNED(OBJECT, FEATURE, GETTER1, GETTER2) \
  case static_cast<uint8_t>(InputFeaturesD0ToKPi::FEATURE): {           \
    if (pdgCode == o2::constants::physics::kD0) {                       \
      inputFeatures.emplace_back(OBJECT.GETTER1());                     \
    } else {                                                            \
      inputFeatures.emplace_back(OBJECT.GETTER2());                     \
    }                                                                   \
    break;                                                              \
  }

// Variation of CHECK_AND_FILL_VEC_D0_SIGNED(OBJECT, FEATURE, GETTER1, GETTER2)
// where GETTER1 and GETTER2 are methods of the OBJECT, the variable
// is filled depending on whether it is a D0 or a D0bar
// and INDEX is the index of the vector
#define CHECK_AND_FILL_VEC_D0_ML(OBJECT, FEATURE, GETTER1, GETTER2, INDEX) \
  case static_cast<uint8_t>(InputFeaturesD0ToKPi::FEATURE): {              \
    if constexpr (usingMl) {                                               \
      if (pdgCode == o2::constants::physics::kD0) {                        \
        inputFeatures.emplace_back(OBJECT.GETTER1()[INDEX]);               \
      } else {                                                             \
        inputFeatures.emplace_back(OBJECT.GETTER2()[INDEX]);               \
      }                                                                    \
    }                                                                      \
    break;                                                                 \
  }

namespace o2::analysis
{
enum class InputFeaturesD0ToKPi : uint8_t {
  chi2PCA = 0,
  decayLength,
  decayLengthXY,
  decayLengthNormalised,
  decayLengthXYNormalised,
  impactParameterXYNormalised0,
  ptProng0,
  ptProng1,
  impactParameterXY0,
  impactParameterXY1,
  impactParameterZ0,
  impactParameterZ1,
  nSigTpcPi0,
  nSigTpcKa0,
  nSigTofPi0,
  nSigTofKa0,
  nSigTpcTofPi0,
  nSigTpcTofKa0,
  nSigTpcPi1,
  nSigTpcKa1,
  nSigTofPi1,
  nSigTofKa1,
  nSigTpcTofPi1,
  nSigTpcTofKa1,
  nSigTpcPiExpPi,
  nSigTpcKaExpPi,
  nSigTpcPiExpKa,
  nSigTpcKaExpKa,
  nSigTofPiExpPi,
  nSigTofKaExpPi,
  nSigTofPiExpKa,
  nSigTofKaExpKa,
  nSigTpcTofPiExpPi,
  nSigTpcTofKaExpPi,
  nSigTpcTofPiExpKa,
  nSigTpcTofKaExpKa,
  maxNormalisedDeltaIP,
  impactParameterProduct,
  bdtOutputBkg,
  bdtOutputNonPrompt,
  bdtOutputPrompt,
  cosThetaStar,
  cpa,
  cpaXY,
  ct
};

template <typename TypeOutputScore = float>
class HfMlResponseD0ToKPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseD0ToKPi() = default;
  /// Default destructor
  virtual ~HfMlResponseD0ToKPi() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the D0 candidate
  /// \return inputFeatures vector
  template <bool usingMl = false, typename T1>
  std::vector<float> getInputFeatures(T1 const& candidate, int const& pdgCode)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_D0(chi2PCA);
        CHECK_AND_FILL_VEC_D0(decayLength);
        CHECK_AND_FILL_VEC_D0(decayLengthXY);
        CHECK_AND_FILL_VEC_D0(decayLengthNormalised);
        CHECK_AND_FILL_VEC_D0(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_D0(ptProng0);
        CHECK_AND_FILL_VEC_D0(ptProng1);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, impactParameterXY0, impactParameter0);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, impactParameterXY1, impactParameter1);
        CHECK_AND_FILL_VEC_D0(impactParameterZ0);
        CHECK_AND_FILL_VEC_D0(impactParameterZ1);
        // TPC PID variables
        CHECK_AND_FILL_VEC_D0_FULL(candidate, nSigTpcPi0, /*getter*/ nSigTpcPi0);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, nSigTpcKa0, /*getter*/ nSigTpcKa0);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, nSigTpcPi1, /*getter*/ nSigTpcPi1);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, nSigTpcKa1, /*getter*/ nSigTpcKa1);
        CHECK_AND_FILL_VEC_D0_SIGNED(candidate, nSigTpcPiExpPi, nSigTpcPi0, nSigTpcPi1);
        CHECK_AND_FILL_VEC_D0_SIGNED(candidate, nSigTpcKaExpPi, nSigTpcKa0, nSigTpcKa1);
        CHECK_AND_FILL_VEC_D0_SIGNED(candidate, nSigTpcPiExpKa, nSigTpcPi1, nSigTpcPi0);
        CHECK_AND_FILL_VEC_D0_SIGNED(candidate, nSigTpcKaExpKa, nSigTpcKa1, nSigTpcKa0);
        // TOF PID variables
        CHECK_AND_FILL_VEC_D0_FULL(candidate, nSigTofPi0, /*getter*/ nSigTofPi0);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, nSigTofKa0, /*getter*/ nSigTofKa0);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, nSigTofPi1, /*getter*/ nSigTofPi1);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, nSigTofKa1, /*getter*/ nSigTofKa1);
        CHECK_AND_FILL_VEC_D0_SIGNED(candidate, nSigTofPiExpPi, nSigTofPi0, nSigTofPi1);
        CHECK_AND_FILL_VEC_D0_SIGNED(candidate, nSigTofKaExpPi, nSigTofKa0, nSigTofKa1);
        CHECK_AND_FILL_VEC_D0_SIGNED(candidate, nSigTofPiExpKa, nSigTofPi1, nSigTofPi0);
        CHECK_AND_FILL_VEC_D0_SIGNED(candidate, nSigTofKaExpKa, nSigTofKa1, nSigTofKa0);
        // Combined PID variables
        CHECK_AND_FILL_VEC_D0_FULL(candidate, nSigTpcTofPi0, tpcTofNSigmaPi0);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, nSigTpcTofKa0, tpcTofNSigmaKa0);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, nSigTpcTofPi1, tpcTofNSigmaPi1);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, nSigTpcTofKa1, tpcTofNSigmaKa1);
        CHECK_AND_FILL_VEC_D0_SIGNED(candidate, nSigTpcTofPiExpPi, tpcTofNSigmaPi0, tpcTofNSigmaPi1);
        CHECK_AND_FILL_VEC_D0_SIGNED(candidate, nSigTpcTofKaExpPi, tpcTofNSigmaKa0, tpcTofNSigmaKa1);
        CHECK_AND_FILL_VEC_D0_SIGNED(candidate, nSigTpcTofPiExpKa, tpcTofNSigmaPi1, tpcTofNSigmaPi0);
        CHECK_AND_FILL_VEC_D0_SIGNED(candidate, nSigTpcTofKaExpKa, tpcTofNSigmaKa1, tpcTofNSigmaKa0);

        CHECK_AND_FILL_VEC_D0_ML(candidate, bdtOutputBkg, mlProbD0, mlProbD0bar, 0);
        CHECK_AND_FILL_VEC_D0_ML(candidate, bdtOutputNonPrompt, mlProbD0, mlProbD0bar, 1);
        CHECK_AND_FILL_VEC_D0_ML(candidate, bdtOutputPrompt, mlProbD0, mlProbD0bar, 2);

        CHECK_AND_FILL_VEC_D0(maxNormalisedDeltaIP);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, impactParameterProduct, impactParameterProduct);
        CHECK_AND_FILL_VEC_D0_HFHELPER_SIGNED(candidate, cosThetaStar, cosThetaStarD0, cosThetaStarD0bar);
        CHECK_AND_FILL_VEC_D0(cpa);
        CHECK_AND_FILL_VEC_D0(cpaXY);
        CHECK_AND_FILL_VEC_D0_HFHELPER(candidate, ct, ctD0);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_D0(chi2PCA),
      FILL_MAP_D0(decayLength),
      FILL_MAP_D0(decayLengthXY),
      FILL_MAP_D0(decayLengthNormalised),
      FILL_MAP_D0(decayLengthXYNormalised),
      FILL_MAP_D0(ptProng0),
      FILL_MAP_D0(ptProng1),
      FILL_MAP_D0(impactParameterXY0),
      FILL_MAP_D0(impactParameterXY1),
      FILL_MAP_D0(impactParameterZ0),
      FILL_MAP_D0(impactParameterZ1),
      // TPC PID variables
      FILL_MAP_D0(nSigTpcPi0),
      FILL_MAP_D0(nSigTpcKa0),
      FILL_MAP_D0(nSigTpcPi1),
      FILL_MAP_D0(nSigTpcKa1),
      FILL_MAP_D0(nSigTpcPiExpPi),
      FILL_MAP_D0(nSigTpcKaExpPi),
      FILL_MAP_D0(nSigTpcPiExpKa),
      FILL_MAP_D0(nSigTpcKaExpKa),
      // TOF PID variables
      FILL_MAP_D0(nSigTofPi0),
      FILL_MAP_D0(nSigTofKa0),
      FILL_MAP_D0(nSigTofPi1),
      FILL_MAP_D0(nSigTofKa1),
      FILL_MAP_D0(nSigTofPiExpPi),
      FILL_MAP_D0(nSigTofKaExpPi),
      FILL_MAP_D0(nSigTofPiExpKa),
      FILL_MAP_D0(nSigTofKaExpKa),
      // Combined PID variables
      FILL_MAP_D0(nSigTpcTofPi0),
      FILL_MAP_D0(nSigTpcTofKa0),
      FILL_MAP_D0(nSigTpcTofPi1),
      FILL_MAP_D0(nSigTpcTofKa1),
      FILL_MAP_D0(nSigTpcTofPiExpPi),
      FILL_MAP_D0(nSigTpcTofKaExpPi),
      FILL_MAP_D0(nSigTpcTofPiExpKa),
      FILL_MAP_D0(nSigTpcTofKaExpKa),
      // ML variables
      FILL_MAP_D0(bdtOutputBkg),
      FILL_MAP_D0(bdtOutputNonPrompt),
      FILL_MAP_D0(bdtOutputPrompt),

      FILL_MAP_D0(maxNormalisedDeltaIP),
      FILL_MAP_D0(impactParameterProduct),
      FILL_MAP_D0(cosThetaStar),
      FILL_MAP_D0(cpa),
      FILL_MAP_D0(cpaXY),
      FILL_MAP_D0(ct)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_D0
#undef CHECK_AND_FILL_VEC_D0_FULL
#undef CHECK_AND_FILL_VEC_D0
#undef CHECK_AND_FILL_VEC_D0_HFHELPER
#undef CHECK_AND_FILL_VEC_D0_HFHELPER_SIGNED
#undef CHECK_AND_FILL_VEC_D0_OBJECT_HFHELPER_SIGNED
#undef CHECK_AND_FILL_VEC_D0_ML

#endif // PWGHF_CORE_HFMLRESPONSED0TOKPI_H_
