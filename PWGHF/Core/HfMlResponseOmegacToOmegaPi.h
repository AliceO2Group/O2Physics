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

/// \file HfMlResponseOmegacToOmegaPi.h
/// \brief Class to compute the ML response for Ωc0 → Ω∓ π± analysis selections
/// \author Yunfan Liu <yunfan.liu@cern.ch>

#ifndef PWGHF_CORE_HFMLRESPONSEOMEGACTOOMEGAPI_H_
#define PWGHF_CORE_HFMLRESPONSEOMEGACTOOMEGAPI_H_

#include <map>
#include <string>
#include <vector>

#include "CommonConstants/PhysicsConstants.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponse.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_OMEGAC0(FEATURE)                                           \
  {                                                                         \
#FEATURE, static_cast < uint8_t>(InputFeaturesOmegacToOmegaPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_OMEGAC0_FULL(OBJECT, FEATURE, GETTER)      \
  case static_cast<uint8_t>(InputFeaturesOmegacToOmegaPi::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());                      \
    break;                                                            \
  }

// Specific case of CHECK_AND_FILL_VEC_OMEGAC0_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_OMEGAC0(GETTER)                           \
  case static_cast<uint8_t>(InputFeaturesOmegacToOmegaPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());                  \
    break;                                                           \
  }

// Variation of CHECK_AND_FILL_VEC_OMEGAC0_FULL(OBJECT, FEATURE, GETTER)
// where GETTER is a method of hfHelper
#define CHECK_AND_FILL_VEC_OMEGAC0_HFHELPER(OBJECT, FEATURE, GETTER)  \
  case static_cast<uint8_t>(InputFeaturesOmegacToOmegaPi::FEATURE): { \
    inputFeatures.emplace_back(hfHelper.GETTER(OBJECT));              \
    break;                                                            \
  }

// Variation of CHECK_AND_FILL_VEC_OMEGAC0_HFHELPER(OBJECT, FEATURE, GETTER)
// where GETTER1 and GETTER2 are methods of hfHelper, and the variable
// is filled depending on whether it is a OMEGAC0 or a OMEGAC0bar
#define CHECK_AND_FILL_VEC_OMEGAC0_HFHELPER_SIGNED(OBJECT, FEATURE, GETTER1, GETTER2) \
  case static_cast<uint8_t>(InputFeaturesOmegacToOmegaPi::FEATURE): {                 \
    if (pdgCode == o2::constants::physics::kOmegaC0) {                                \
      inputFeatures.emplace_back(hfHelper.GETTER1(OBJECT));                           \
    } else {                                                                          \
      inputFeatures.emplace_back(hfHelper.GETTER2(OBJECT));                           \
    }                                                                                 \
    break;                                                                            \
  }

namespace o2::analysis
{
enum class InputFeaturesOmegacToOmegaPi : uint8_t {
  xPv = 0,
  yPv,
  zPv,
  xDecayVtxCharmBaryon,
  yDecayVtxCharmBaryon,
  zDecayVtxCharmBaryon,
  xDecayVtxCascade,
  yDecayVtxCascade,
  zDecayVtxCascade,
  xDecayVtxV0,
  yDecayVtxV0,
  zDecayVtxV0,
  signDecay,
  covVtxCharmBaryon0,
  covVtxCharmBaryon1,
  covVtxCharmBaryon2,
  covVtxCharmBaryon3,
  covVtxCharmBaryon4,
  covVtxCharmBaryon5,
  pxCharmBaryon,
  pyCharmBaryon,
  pzCharmBaryon,
  pxCasc,
  pyCasc,
  pzCasc,
  pxPiFromCharmBaryon,
  pyPiFromCharmBaryon,
  pzPiFromCharmBaryon,
  pxLambda,
  pyLambda,
  pzLambda,
  pxKaFromCasc,
  pyKaFromCasc,
  pzKaFromCasc,
  pxPosV0Dau,
  pyPosV0Dau,
  pzPosV0Dau,
  pxNegV0Dau,
  pyNegV0Dau,
  pzNegV0Dau,
  impactParCascXY,
  impactParPiFromCharmBaryonXY,
  errImpactParCascXY,
  errImpactParPiFromCharmBaryonXY,
  invMassLambda,
  invMassCascade,
  invMassCharmBaryon,
  cosPAV0,
  cosPACharmBaryon,
  cosPACasc,
  cosPAXYV0,
  cosPAXYCharmBaryon,
  cosPAXYCasc,
  ctauOmegac,
  ctauCascade,
  ctauV0,
  etaV0PosDau,
  etaV0NegDau,
  etaKaFromCasc,
  etaPiFromCharmBaryon,
  etaCharmBaryon,
  etaCascade,
  etaV0,
  dcaXYToPvV0Dau0,
  dcaXYToPvV0Dau1,
  dcaXYToPvCascDau,
  dcaCascDau,
  dcaV0Dau,
  dcaCharmBaryonDau,
  decLenCharmBaryon,
  decLenCascade,
  decLenV0,
  errorDecayLengthCharmBaryon,
  errorDecayLengthXYCharmBaryon,
  ptCharmBaryon,
  ptCasc,
  ptPiFromCharmBaryon,
  nSigTpcPi0,
  nSigTpcPr0,
  nSigTofPi0,
  nSigTofPr0,
  nSigTpcPi1,
  nSigTpcPr1,
  nSigTofPi1,
  nSigTofPr1,
  nSigTpcKaCasc,
  nSigTofKaCasc,
  nSigTpcPiCharmBaryon,
  nSigTofPiCharmBaryon
};

template <typename TypeOutputScore = float>
class HfMlResponseOmegaCToOmegaPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseOmegaCToOmegaPi() = default;
  /// Default destructor
  virtual ~HfMlResponseOmegaCToOmegaPi() = default;

  HfHelper hfHelper;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the OMEGAC0 candidate
  /// \param lamProng0 is the candidate's lamProng0
  /// \param lamProng1 is the candidate's lamProng1
  /// \return inputFeatures vector
  template <typename T1, typename T2>
  std::vector<float> getInputFeatures(T1 const& candidate, T2 const& lamProng0, T2 const& lamProng1, T2 const& cascProng, T2 const& charmBaryonProng, int const& pdgCode)
  // std::vector<float> getInputFeatures(T1 const& candidate,T2 const& lamProng0, T2 const& lamProng1, int const& pdgCode)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_OMEGAC0(xPv);
        CHECK_AND_FILL_VEC_OMEGAC0(yPv);
        CHECK_AND_FILL_VEC_OMEGAC0(zPv);
        CHECK_AND_FILL_VEC_OMEGAC0(xDecayVtxCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(yDecayVtxCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(zDecayVtxCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(xDecayVtxCascade);
        CHECK_AND_FILL_VEC_OMEGAC0(yDecayVtxCascade);
        CHECK_AND_FILL_VEC_OMEGAC0(zDecayVtxCascade);
        CHECK_AND_FILL_VEC_OMEGAC0(xDecayVtxV0);
        CHECK_AND_FILL_VEC_OMEGAC0(yDecayVtxV0);
        CHECK_AND_FILL_VEC_OMEGAC0(zDecayVtxV0);
        CHECK_AND_FILL_VEC_OMEGAC0(signDecay);
        CHECK_AND_FILL_VEC_OMEGAC0(covVtxCharmBaryon0);
        CHECK_AND_FILL_VEC_OMEGAC0(covVtxCharmBaryon1);
        CHECK_AND_FILL_VEC_OMEGAC0(covVtxCharmBaryon2);
        CHECK_AND_FILL_VEC_OMEGAC0(covVtxCharmBaryon3);
        CHECK_AND_FILL_VEC_OMEGAC0(covVtxCharmBaryon4);
        CHECK_AND_FILL_VEC_OMEGAC0(covVtxCharmBaryon5);
        CHECK_AND_FILL_VEC_OMEGAC0(pxCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(pyCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(pzCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(pxCasc);
        CHECK_AND_FILL_VEC_OMEGAC0(pyCasc);
        CHECK_AND_FILL_VEC_OMEGAC0(pzCasc);
        CHECK_AND_FILL_VEC_OMEGAC0(pxPiFromCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(pyPiFromCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(pzPiFromCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(pxLambda);
        CHECK_AND_FILL_VEC_OMEGAC0(pyLambda);
        CHECK_AND_FILL_VEC_OMEGAC0(pzLambda);
        CHECK_AND_FILL_VEC_OMEGAC0(pxKaFromCasc);
        CHECK_AND_FILL_VEC_OMEGAC0(pyKaFromCasc);
        CHECK_AND_FILL_VEC_OMEGAC0(pzKaFromCasc);
        CHECK_AND_FILL_VEC_OMEGAC0(pxPosV0Dau);
        CHECK_AND_FILL_VEC_OMEGAC0(pyPosV0Dau);
        CHECK_AND_FILL_VEC_OMEGAC0(pzPosV0Dau);
        CHECK_AND_FILL_VEC_OMEGAC0(pxNegV0Dau);
        CHECK_AND_FILL_VEC_OMEGAC0(pyNegV0Dau);
        CHECK_AND_FILL_VEC_OMEGAC0(pzNegV0Dau);
        CHECK_AND_FILL_VEC_OMEGAC0(impactParCascXY);
        CHECK_AND_FILL_VEC_OMEGAC0(impactParPiFromCharmBaryonXY);
        CHECK_AND_FILL_VEC_OMEGAC0(errImpactParCascXY);
        CHECK_AND_FILL_VEC_OMEGAC0(errImpactParPiFromCharmBaryonXY);
        CHECK_AND_FILL_VEC_OMEGAC0(invMassLambda);
        CHECK_AND_FILL_VEC_OMEGAC0(invMassCascade);
        CHECK_AND_FILL_VEC_OMEGAC0(invMassCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(cosPAV0);
        CHECK_AND_FILL_VEC_OMEGAC0(cosPACharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(cosPACasc);
        CHECK_AND_FILL_VEC_OMEGAC0(cosPAXYV0);
        CHECK_AND_FILL_VEC_OMEGAC0(cosPAXYCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(cosPAXYCasc);
        CHECK_AND_FILL_VEC_OMEGAC0(ctauOmegac);
        CHECK_AND_FILL_VEC_OMEGAC0(ctauCascade);
        CHECK_AND_FILL_VEC_OMEGAC0(ctauV0);
        CHECK_AND_FILL_VEC_OMEGAC0(etaV0PosDau);
        CHECK_AND_FILL_VEC_OMEGAC0(etaV0NegDau);
        CHECK_AND_FILL_VEC_OMEGAC0(etaKaFromCasc);
        CHECK_AND_FILL_VEC_OMEGAC0(etaPiFromCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(etaCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(etaCascade);
        CHECK_AND_FILL_VEC_OMEGAC0(etaV0);
        CHECK_AND_FILL_VEC_OMEGAC0(dcaXYToPvV0Dau0);
        CHECK_AND_FILL_VEC_OMEGAC0(dcaXYToPvV0Dau1);
        CHECK_AND_FILL_VEC_OMEGAC0(dcaXYToPvCascDau);
        CHECK_AND_FILL_VEC_OMEGAC0(dcaCascDau);
        CHECK_AND_FILL_VEC_OMEGAC0(dcaV0Dau);
        CHECK_AND_FILL_VEC_OMEGAC0(dcaCharmBaryonDau);
        CHECK_AND_FILL_VEC_OMEGAC0(decLenCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(decLenCascade);
        CHECK_AND_FILL_VEC_OMEGAC0(decLenV0);
        CHECK_AND_FILL_VEC_OMEGAC0(errorDecayLengthCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(errorDecayLengthXYCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(ptCharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0(ptCasc);
        CHECK_AND_FILL_VEC_OMEGAC0(ptPiFromCharmBaryon);
        // CHECK_AND_FILL_VEC_OMEGAC0(decayLength);
        // CHECK_AND_FILL_VEC_OMEGAC0(decayLengthXY);
        // CHECK_AND_FILL_VEC_OMEGAC0(decayLengthNormalised);
        // CHECK_AND_FILL_VEC_OMEGAC0(decayLengthXYNormalised);
        // CHECK_AND_FILL_VEC_OMEGAC0(ptProng0);
        // CHECK_AND_FILL_VEC_OMEGAC0(ptProng1);
        // CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, impactParameter0, impactParameter0);
        // CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, impactParameter1, impactParameter1);
        //  TPC PID variables
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProng0, nSigTpcPi0, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProng0, nSigTpcPr0, tpcNSigmaPr);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProng1, nSigTpcPi1, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProng1, nSigTpcPr1, tpcNSigmaPr);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(cascProng, nSigTpcKaCasc, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(charmBaryonProng, nSigTpcPiCharmBaryon, tpcNSigmaPi);
        // TOF PID variables
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProng0, nSigTofPi0, tofNSigmaPi);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProng0, nSigTofPr0, tofNSigmaPr);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProng1, nSigTofPi1, tofNSigmaPi);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProng1, nSigTofPr1, tofNSigmaPr);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(cascProng, nSigTofKaCasc, tofNSigmaKa);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(charmBaryonProng, nSigTofPiCharmBaryon, tofNSigmaPi);
        // Combined PID variables
        // CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProng0, nSigTpcTofPi0, tpcTofNSigmaPi);
        // CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProng0, nSigTpcTofPr0, tpcTofNSigmaPr);
        // CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProng1, nSigTpcTofPi1, tpcTofNSigmaPi);
        // CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProng1, nSigTpcTofPr1, tpcTofNSigmaPr);
        // CHECK_AND_FILL_VEC_OMEGAC0_FULL(cascProng, nSigTpcTofKaCasc, tpcTofNSigmaKa);
        // CHECK_AND_FILL_VEC_OMEGAC0_FULL(charmBaryonProng, nSigTpcTofPiCharmBaryon, tpcTofNSigmaPi);
        //
        // CHECK_AND_FILL_VEC_OMEGAC0(maxNormalisedDeltaIP);
        // CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, impactParameterProduct, impactParameterProduct);
        // CHECK_AND_FILL_VEC_OMEGAC0_HFHELPER_SIGNED(candidate, cosThetaStar, cosThetaStarOMEGAC0, cosThetaStarOMEGAC0bar);
        // CHECK_AND_FILL_VEC_OMEGAC0(cpa);
        // CHECK_AND_FILL_VEC_OMEGAC0(cpaXY);
        // CHECK_AND_FILL_VEC_OMEGAC0_HFHELPER(candidate, ct, ctOMEGAC0);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_OMEGAC0(xPv),
      FILL_MAP_OMEGAC0(yPv),
      FILL_MAP_OMEGAC0(zPv),
      FILL_MAP_OMEGAC0(xDecayVtxCharmBaryon),
      FILL_MAP_OMEGAC0(yDecayVtxCharmBaryon),
      FILL_MAP_OMEGAC0(zDecayVtxCharmBaryon),
      FILL_MAP_OMEGAC0(xDecayVtxCascade),
      FILL_MAP_OMEGAC0(yDecayVtxCascade),
      FILL_MAP_OMEGAC0(zDecayVtxCascade),
      FILL_MAP_OMEGAC0(xDecayVtxV0),
      FILL_MAP_OMEGAC0(yDecayVtxV0),
      FILL_MAP_OMEGAC0(zDecayVtxV0),
      FILL_MAP_OMEGAC0(signDecay),
      FILL_MAP_OMEGAC0(covVtxCharmBaryon0),
      FILL_MAP_OMEGAC0(covVtxCharmBaryon1),
      FILL_MAP_OMEGAC0(covVtxCharmBaryon2),
      FILL_MAP_OMEGAC0(covVtxCharmBaryon3),
      FILL_MAP_OMEGAC0(covVtxCharmBaryon4),
      FILL_MAP_OMEGAC0(covVtxCharmBaryon5),
      FILL_MAP_OMEGAC0(pxCharmBaryon),
      FILL_MAP_OMEGAC0(pyCharmBaryon),
      FILL_MAP_OMEGAC0(pzCharmBaryon),
      FILL_MAP_OMEGAC0(pxCasc),
      FILL_MAP_OMEGAC0(pyCasc),
      FILL_MAP_OMEGAC0(pzCasc),
      FILL_MAP_OMEGAC0(pxPiFromCharmBaryon),
      FILL_MAP_OMEGAC0(pyPiFromCharmBaryon),
      FILL_MAP_OMEGAC0(pzPiFromCharmBaryon),
      FILL_MAP_OMEGAC0(pxLambda),
      FILL_MAP_OMEGAC0(pyLambda),
      FILL_MAP_OMEGAC0(pzLambda),
      FILL_MAP_OMEGAC0(pxKaFromCasc),
      FILL_MAP_OMEGAC0(pyKaFromCasc),
      FILL_MAP_OMEGAC0(pzKaFromCasc),
      FILL_MAP_OMEGAC0(pxPosV0Dau),
      FILL_MAP_OMEGAC0(pyPosV0Dau),
      FILL_MAP_OMEGAC0(pzPosV0Dau),
      FILL_MAP_OMEGAC0(pxNegV0Dau),
      FILL_MAP_OMEGAC0(pyNegV0Dau),
      FILL_MAP_OMEGAC0(pzNegV0Dau),
      FILL_MAP_OMEGAC0(impactParCascXY),
      FILL_MAP_OMEGAC0(impactParPiFromCharmBaryonXY),
      FILL_MAP_OMEGAC0(errImpactParCascXY),
      FILL_MAP_OMEGAC0(errImpactParPiFromCharmBaryonXY),
      FILL_MAP_OMEGAC0(invMassLambda),
      FILL_MAP_OMEGAC0(invMassCascade),
      FILL_MAP_OMEGAC0(invMassCharmBaryon),
      FILL_MAP_OMEGAC0(cosPAV0),
      FILL_MAP_OMEGAC0(cosPACharmBaryon),
      FILL_MAP_OMEGAC0(cosPACasc),
      FILL_MAP_OMEGAC0(cosPAXYV0),
      FILL_MAP_OMEGAC0(cosPAXYCharmBaryon),
      FILL_MAP_OMEGAC0(cosPAXYCasc),
      FILL_MAP_OMEGAC0(ctauOmegac),
      FILL_MAP_OMEGAC0(ctauCascade),
      FILL_MAP_OMEGAC0(ctauV0),
      FILL_MAP_OMEGAC0(etaV0PosDau),
      FILL_MAP_OMEGAC0(etaV0NegDau),
      FILL_MAP_OMEGAC0(etaKaFromCasc),
      FILL_MAP_OMEGAC0(etaPiFromCharmBaryon),
      FILL_MAP_OMEGAC0(etaCharmBaryon),
      FILL_MAP_OMEGAC0(etaCascade),
      FILL_MAP_OMEGAC0(etaV0),
      FILL_MAP_OMEGAC0(dcaXYToPvV0Dau0),
      FILL_MAP_OMEGAC0(dcaXYToPvV0Dau1),
      FILL_MAP_OMEGAC0(dcaXYToPvCascDau),
      FILL_MAP_OMEGAC0(dcaCascDau),
      FILL_MAP_OMEGAC0(dcaV0Dau),
      FILL_MAP_OMEGAC0(dcaCharmBaryonDau),
      FILL_MAP_OMEGAC0(decLenCharmBaryon),
      FILL_MAP_OMEGAC0(decLenCascade),
      FILL_MAP_OMEGAC0(decLenV0),
      FILL_MAP_OMEGAC0(errorDecayLengthCharmBaryon),
      FILL_MAP_OMEGAC0(errorDecayLengthXYCharmBaryon),
      FILL_MAP_OMEGAC0(ptCharmBaryon),
      FILL_MAP_OMEGAC0(ptCasc),
      FILL_MAP_OMEGAC0(ptPiFromCharmBaryon),
      // FILL_MAP_OMEGAC0(decayLength),
      // FILL_MAP_OMEGAC0(decayLengthXY),
      // FILL_MAP_OMEGAC0(decayLengthNormalised),
      // FILL_MAP_OMEGAC0(decayLengthXYNormalised),
      // FILL_MAP_OMEGAC0(ptProng0),
      // FILL_MAP_OMEGAC0(ptProng1),
      // FILL_MAP_OMEGAC0(impactParameter0),
      // FILL_MAP_OMEGAC0(impactParameter1),
      //  TPC PID variables
      FILL_MAP_OMEGAC0(nSigTpcPi0),
      FILL_MAP_OMEGAC0(nSigTpcPr0),
      FILL_MAP_OMEGAC0(nSigTpcPi1),
      FILL_MAP_OMEGAC0(nSigTpcPr1),
      FILL_MAP_OMEGAC0(nSigTpcKaCasc),
      FILL_MAP_OMEGAC0(nSigTpcPiCharmBaryon),
      // TOF PID variables
      FILL_MAP_OMEGAC0(nSigTofPi0),
      FILL_MAP_OMEGAC0(nSigTofPr0),
      FILL_MAP_OMEGAC0(nSigTofPi1),
      FILL_MAP_OMEGAC0(nSigTofPr1),
      FILL_MAP_OMEGAC0(nSigTofKaCasc),
      FILL_MAP_OMEGAC0(nSigTofPiCharmBaryon)
      // Combined PID variables
      // FILL_MAP_OMEGAC0(nSigTpcTofPi0),
      // FILL_MAP_OMEGAC0(nSigTpcTofPr0),
      // FILL_MAP_OMEGAC0(nSigTpcTofPi1),
      // FILL_MAP_OMEGAC0(nSigTpcTofPr1),
      // FILL_MAP_OMEGAC0(nSigTpcTofKaCasc),
      // FILL_MAP_OMEGAC0(nSigTpcTofPiCharmBaryon)
      //
      // FILL_MAP_OMEGAC0(maxNormalisedDeltaIP),
      // FILL_MAP_OMEGAC0(impactParameterProduct),
      // FILL_MAP_OMEGAC0(cosThetaStar),
      // FILL_MAP_OMEGAC0(cpa),
      // FILL_MAP_OMEGAC0(cpaXY),
      // FILL_MAP_OMEGAC0(ct)
    };
  }
};

} // namespace o2::analysis

#undef FILL_MAP_OMEGAC0
#undef CHECK_AND_FILL_VEC_OMEGAC0_FULL
#undef CHECK_AND_FILL_VEC_OMEGAC0
#undef CHECK_AND_FILL_VEC_OMEGAC0_HFHELPER
#undef CHECK_AND_FILL_VEC_OMEGAC0_HFHELPER_SIGNED

#endif // PWGHF_CORE_HFMLRESPONSEOMEGACTOOMEGAPI_H_
