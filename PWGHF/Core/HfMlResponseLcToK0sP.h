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

/// \file HfMlResponseLcToK0sP.h
/// \brief Class to compute the ML response for Lc± → K0s p analysis selections
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg
/// \author Elisa Meninno <elisa.meninno@cern.ch>, SMI Vienna

#ifndef PWGHF_CORE_HFMLRESPONSELCTOK0SP_H_
#define PWGHF_CORE_HFMLRESPONSELCTOK0SP_H_

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponse.h"

#include "Tools/ML/MlResponse.h"

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_LC(FEATURE)                                       \
  {                                                                \
    #FEATURE, static_cast<uint8_t>(InputFeaturesLcToK0sP::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_LC_FULL(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesLcToK0sP::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());               \
    break;                                                     \
  }

// Specific case of CHECK_AND_FILL_VEC_LC_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_LC(GETTER)                         \
  case static_cast<uint8_t>(InputFeaturesLcToK0sP::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());           \
    break;                                                    \
  }

// Variation of CHECK_AND_FILL_VEC_LC_FULL(OBJECT, FEATURE, GETTER)
// where GETTER is a method of HfHelper
#define CHECK_AND_FILL_VEC_LC_HFHELPER(OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesLcToK0sP::FEATURE): {  \
    inputFeatures.emplace_back(HfHelper::GETTER(OBJECT));       \
    break;                                                      \
  }

namespace o2::analysis
{
// possible input features for ML
enum class InputFeaturesLcToK0sP : uint8_t {
  chi2PCA,
  rSecondaryVertex,
  decayLength,
  decayLengthXY,
  decayLengthNormalised,
  decayLengthXYNormalised,
  impactParameterNormalised0,
  ptProng0,
  impactParameterNormalised1,
  ptProng1,
  impactParameter0,
  impactParameter1,
  v0Radius,
  v0cosPA,
  v0MLambda,
  v0MAntiLambda,
  v0MK0Short,
  v0MGamma,
  ctV0,
  dcaV0daughters,
  ptV0Pos,
  dcaPosToPV,
  ptV0Neg,
  dcaNegToPV,
  nSigmaTpcPr0,
  nSigmaTofPr0,
  nSigmaTpcTofPr0,
  cpa,
  cpaXY,
  ct
};

template <typename TypeOutputScore = float>
class HfMlResponseLcToK0sP : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseLcToK0sP() = default;
  /// Default destructor
  virtual ~HfMlResponseLcToK0sP() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Lc candidate
  /// \param bach is the bachelor candidate (proton)
  /// \return inputFeatures vector
  template <typename T1, typename T2>
  std::vector<float> getInputFeatures(T1 const& candidate,
                                      T2 const& bach)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_LC(chi2PCA);
        CHECK_AND_FILL_VEC_LC(rSecondaryVertex);
        CHECK_AND_FILL_VEC_LC(decayLength);
        CHECK_AND_FILL_VEC_LC(decayLengthXY);
        CHECK_AND_FILL_VEC_LC(decayLengthNormalised);
        CHECK_AND_FILL_VEC_LC(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_LC(impactParameterNormalised0);
        CHECK_AND_FILL_VEC_LC(ptProng0);
        CHECK_AND_FILL_VEC_LC(impactParameterNormalised1);
        CHECK_AND_FILL_VEC_LC(ptProng1);
        CHECK_AND_FILL_VEC_LC(impactParameter0);
        CHECK_AND_FILL_VEC_LC(impactParameter1);
        CHECK_AND_FILL_VEC_LC_FULL(candidate, v0Radius, v0radius);
        CHECK_AND_FILL_VEC_LC(v0cosPA);
        CHECK_AND_FILL_VEC_LC_FULL(candidate, v0MLambda, mLambda);
        CHECK_AND_FILL_VEC_LC_FULL(candidate, v0MAntiLambda, mAntiLambda);
        CHECK_AND_FILL_VEC_LC_FULL(candidate, v0MK0Short, mK0Short);
        CHECK_AND_FILL_VEC_LC_FULL(candidate, v0MGamma, mGamma);
        CHECK_AND_FILL_VEC_LC_HFHELPER(candidate, ctV0, ctV0K0s);
        // CHECK_AND_FILL_VEC_LC_HFHELPER(candidate, ctV0, ctV0Lambda);
        CHECK_AND_FILL_VEC_LC(dcaV0daughters);
        CHECK_AND_FILL_VEC_LC(ptV0Pos);
        CHECK_AND_FILL_VEC_LC_FULL(candidate, dcaPosToPV, dcapostopv);
        CHECK_AND_FILL_VEC_LC(ptV0Neg);
        CHECK_AND_FILL_VEC_LC_FULL(candidate, dcaNegToPV, dcanegtopv);
        CHECK_AND_FILL_VEC_LC(cpa);
        CHECK_AND_FILL_VEC_LC(cpaXY);
        CHECK_AND_FILL_VEC_LC_HFHELPER(candidate, ct, ctLc);
        // TPC PID variables
        CHECK_AND_FILL_VEC_LC_FULL(bach, nSigmaTpcPr0, tpcNSigmaPr);
        // TOF PID variables
        CHECK_AND_FILL_VEC_LC_FULL(bach, nSigmaTofPr0, tofNSigmaPr);
        // Combined nSigma variable
        CHECK_AND_FILL_VEC_LC_FULL(bach, nSigmaTpcTofPr0, tpcTofNSigmaPr);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_LC(chi2PCA),
      FILL_MAP_LC(rSecondaryVertex),
      FILL_MAP_LC(decayLength),
      FILL_MAP_LC(decayLengthXY),
      FILL_MAP_LC(decayLengthNormalised),
      FILL_MAP_LC(decayLengthXYNormalised),
      FILL_MAP_LC(impactParameterNormalised0),
      FILL_MAP_LC(ptProng0),
      FILL_MAP_LC(impactParameterNormalised1),
      FILL_MAP_LC(ptProng1),
      FILL_MAP_LC(impactParameter0),
      FILL_MAP_LC(impactParameter1),
      FILL_MAP_LC(v0Radius),
      FILL_MAP_LC(v0cosPA),
      FILL_MAP_LC(v0MLambda),
      FILL_MAP_LC(v0MAntiLambda),
      FILL_MAP_LC(v0MK0Short),
      FILL_MAP_LC(v0MGamma),
      FILL_MAP_LC(ctV0),
      FILL_MAP_LC(dcaV0daughters),
      FILL_MAP_LC(ptV0Pos),
      FILL_MAP_LC(dcaPosToPV),
      FILL_MAP_LC(ptV0Neg),
      FILL_MAP_LC(dcaNegToPV),
      FILL_MAP_LC(cpa),
      FILL_MAP_LC(cpaXY),
      FILL_MAP_LC(ct),
      // TPC PID variables
      FILL_MAP_LC(nSigmaTpcPr0),
      // TOF PID variables
      FILL_MAP_LC(nSigmaTofPr0),
      // Combined nSigma variable
      FILL_MAP_LC(nSigmaTpcTofPr0)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_LC
#undef CHECK_AND_FILL_VEC_LC_FULL
#undef CHECK_AND_FILL_VEC_LC
#undef CHECK_AND_FILL_VEC_LC_HFHELPER

#endif // PWGHF_CORE_HFMLRESPONSELCTOK0SP_H_
