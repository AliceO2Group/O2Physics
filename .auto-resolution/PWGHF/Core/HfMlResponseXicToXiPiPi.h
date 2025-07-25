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

/// \file HfMlResponseXicToXiPiPi.h
/// \brief Class to compute the ML response for Ξc± → Ξ∓ π± π± analysis selections
/// \author Jaeyoon Cho <jaeyoon.cho@cern.ch>, Inha University

#ifndef PWGHF_CORE_HFMLRESPONSEXICTOXIPIPI_H_
#define PWGHF_CORE_HFMLRESPONSEXICTOXIPIPI_H_

#include "PWGHF/Core/HfMlResponse.h"

#include "Tools/ML/MlResponse.h"

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_XICTOXIPIPI(FEATURE)                                 \
  {                                                                   \
    #FEATURE, static_cast<uint8_t>(InputFeaturesXicToXiPiPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_XICTOXIPIPI_FULL(OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesXicToXiPiPi::FEATURE): {    \
    inputFeatures.emplace_back(OBJECT.GETTER());                     \
    break;                                                           \
  }

// Specific case of CHECK_AND_FILL_VEC_XICTOXIPIPI_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_XICTOXIPIPI(GETTER)                   \
  case static_cast<uint8_t>(InputFeaturesXicToXiPiPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());              \
    break;                                                       \
  }

namespace o2::analysis
{

enum class InputFeaturesXicToXiPiPi : uint8_t {
  ptProng0 = 0,
  ptProng1,
  ptProng2,
  chi2PCA,
  decayLength,
  decayLengthNormalised,
  decayLengthXY,
  decayLengthXYNormalised,
  cpa,
  cpaXY,
  cpaXi,
  cpaXYXi,
  cpaLambda,
  cpaXYLambda,
  impactParameterXi,
  impactParameterPi0,
  impactParameterPi1,
  invMassXi,
  invMassLambda,
  dcaXiDaughters,
  dcaV0Daughters,
  dcaPosToPV,
  dcaNegToPV,
  dcaBachelorToPV,
  dcaXYCascToPV,
  dcaZCascToPV,
  nSigTpcPiFromXicPlus0,
  nSigTpcPiFromXicPlus1,
  nSigTpcBachelorPi,
  nSigTpcPiFromLambda,
  nSigTpcPrFromLambda,
  nSigTofPiFromXicPlus0,
  nSigTofPiFromXicPlus1,
  nSigTofBachelorPi,
  nSigTofPiFromLambda,
  nSigTofPrFromLambda
};

template <typename TypeOutputScore = float>
class HfMlResponseXicToXiPiPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseXicToXiPiPi() = default;
  /// Default destructor
  virtual ~HfMlResponseXicToXiPiPi() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Xic candidate
  /// \return inputFeatures vector
  template <typename T1>
  std::vector<float> getInputFeatures(T1 const& candidate)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_XICTOXIPIPI(ptProng0);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(ptProng1);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(ptProng2);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(chi2PCA);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(decayLength);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(decayLengthNormalised);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(decayLengthXY);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(cpa);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(cpaXY);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(cpaXi);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(cpaXYXi);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(cpaLambda);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(cpaXYLambda);
        CHECK_AND_FILL_VEC_XICTOXIPIPI_FULL(candidate, impactParameterXi, impactParameter0);
        CHECK_AND_FILL_VEC_XICTOXIPIPI_FULL(candidate, impactParameterPi0, impactParameter1);
        CHECK_AND_FILL_VEC_XICTOXIPIPI_FULL(candidate, impactParameterPi1, impactParameter2);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(invMassXi);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(invMassLambda);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(dcaXiDaughters);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(dcaV0Daughters);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(dcaPosToPV);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(dcaNegToPV);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(dcaBachelorToPV);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(dcaXYCascToPV);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(dcaZCascToPV);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(nSigTpcPiFromXicPlus0);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(nSigTpcPiFromXicPlus1);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(nSigTpcBachelorPi);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(nSigTpcPiFromLambda);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(nSigTpcPrFromLambda);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(nSigTofPiFromXicPlus0);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(nSigTofPiFromXicPlus1);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(nSigTofBachelorPi);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(nSigTofPiFromLambda);
        CHECK_AND_FILL_VEC_XICTOXIPIPI(nSigTofPrFromLambda);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_XICTOXIPIPI(ptProng0),
      FILL_MAP_XICTOXIPIPI(ptProng1),
      FILL_MAP_XICTOXIPIPI(ptProng2),
      FILL_MAP_XICTOXIPIPI(chi2PCA),
      FILL_MAP_XICTOXIPIPI(decayLength),
      FILL_MAP_XICTOXIPIPI(decayLengthNormalised),
      FILL_MAP_XICTOXIPIPI(decayLengthXY),
      FILL_MAP_XICTOXIPIPI(decayLengthXYNormalised),
      FILL_MAP_XICTOXIPIPI(cpa),
      FILL_MAP_XICTOXIPIPI(cpaXY),
      FILL_MAP_XICTOXIPIPI(cpaXi),
      FILL_MAP_XICTOXIPIPI(cpaXYXi),
      FILL_MAP_XICTOXIPIPI(cpaLambda),
      FILL_MAP_XICTOXIPIPI(cpaXYLambda),
      FILL_MAP_XICTOXIPIPI(impactParameterXi),
      FILL_MAP_XICTOXIPIPI(impactParameterPi0),
      FILL_MAP_XICTOXIPIPI(impactParameterPi1),
      FILL_MAP_XICTOXIPIPI(invMassXi),
      FILL_MAP_XICTOXIPIPI(invMassLambda),
      FILL_MAP_XICTOXIPIPI(dcaXiDaughters),
      FILL_MAP_XICTOXIPIPI(dcaV0Daughters),
      FILL_MAP_XICTOXIPIPI(dcaPosToPV),
      FILL_MAP_XICTOXIPIPI(dcaNegToPV),
      FILL_MAP_XICTOXIPIPI(dcaBachelorToPV),
      FILL_MAP_XICTOXIPIPI(dcaXYCascToPV),
      FILL_MAP_XICTOXIPIPI(dcaZCascToPV),
      FILL_MAP_XICTOXIPIPI(nSigTpcPiFromXicPlus0),
      FILL_MAP_XICTOXIPIPI(nSigTpcPiFromXicPlus1),
      FILL_MAP_XICTOXIPIPI(nSigTpcBachelorPi),
      FILL_MAP_XICTOXIPIPI(nSigTpcPiFromLambda),
      FILL_MAP_XICTOXIPIPI(nSigTpcPrFromLambda),
      FILL_MAP_XICTOXIPIPI(nSigTofPiFromXicPlus0),
      FILL_MAP_XICTOXIPIPI(nSigTofPiFromXicPlus1),
      FILL_MAP_XICTOXIPIPI(nSigTofBachelorPi),
      FILL_MAP_XICTOXIPIPI(nSigTofPiFromLambda),
      FILL_MAP_XICTOXIPIPI(nSigTofPrFromLambda)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_XICTOXIPIPI
#undef CHECK_AND_FILL_VEC_XICTOXIPIPI_FULL
#undef CHECK_AND_FILL_VEC_XICTOXIPIPI

#endif // PWGHF_CORE_HFMLRESPONSEXICTOXIPIPI_H_
