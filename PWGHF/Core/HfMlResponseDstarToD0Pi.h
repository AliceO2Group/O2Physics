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

/// \file HfMlResponseDstarToD0Pi.h
/// \brief Class to compute the ML response for D*± → D0(bar) π± analysis selections
/// \author Mingze Li <Mingze.li@cern.ch>

#ifndef PWGHF_CORE_HFMLRESPONSEDSTARTOD0PI_H_
#define PWGHF_CORE_HFMLRESPONSEDSTARTOD0PI_H_

#include <map>
#include <string>
#include <vector>

#include "PWGHF/Core/HfMlResponse.h"
#include "CommonConstants/PhysicsConstants.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_DSTAR(FEATURE)                                         \
  {                                                                     \
#FEATURE, static_cast < uint8_t>(InputFeaturesDstarToD0Pi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_DSTAR_FULL(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesDstarToD0Pi::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());                  \
    break;                                                        \
  }

// Specific case of CHECK_AND_FILL_VEC_DSTAR_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_DSTAR(GETTER)                         \
  case static_cast<uint8_t>(InputFeaturesDstarToD0Pi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());              \
    break;                                                       \
  }

// Very specific case of CHECK_AND_FILL_VEC_DSTAR_FULL(OBJECT, FEATURE, GETTER)
// Use for push back different value for D*+ or D*- candidate
#define CHECK_AND_FILL_VEC_DSTAR_CHARGEBASE(POSGETTER, NEGGETTER, FEATURENAME) \
  case static_cast<uint8_t>(InputFeaturesDstarToD0Pi::FEATURENAME): {          \
    if (candidate.signSoftPi() > 0) {                                          \
      inputFeatures.emplace_back(candidate.POSGETTER());                       \
    } else {                                                                   \
      inputFeatures.emplace_back(candidate.NEGGETTER());                       \
    }                                                                          \
    break;                                                                     \
  }

// Very specific case of CHECK_AND_FILL_VEC_DSTAR_FULL(OBJECT, FEATURE, GETTER)
// Use for push back deltaMassD0 for D*+ or D*- candidate
#define CHECK_AND_FILL_VEC_DSTAR_DELTA_MASS_D0(FEATURENAME)                                  \
  case static_cast<uint8_t>(InputFeaturesDstarToD0Pi::FEATURENAME): {                        \
    if (candidate.signSoftPi() > 0) {                                                        \
      inputFeatures.emplace_back(candidate.invMassD0() - o2::constants::physics::MassD0);    \
    } else {                                                                                 \
      inputFeatures.emplace_back(candidate.invMassD0Bar() - o2::constants::physics::MassD0); \
    }                                                                                        \
    break;                                                                                   \
  }

namespace o2::analysis
{

enum class InputFeaturesDstarToD0Pi : uint8_t {
  chi2PCAD0 = 0,
  decayLengthD0,
  decayLengthXYD0,
  decayLengthNormalisedD0,
  decayLengthXYNormalisedD0,
  cpaD0,
  cpaXYD0,
  deltaIPNormalisedMaxD0,
  impactParameterProductD0,
  ptProng0,
  ptProng1,
  ptSoftPi,
  impactParameter0,
  impactParameter1,
  impParamSoftPi,
  impactParameterNormalised0,
  impactParameterNormalised1,
  normalisedImpParamSoftPi,
  cosThetaStarD0,
  massD0,
  deltaMassD0,
  nSigmaTPCPiPr0,
  nSigmaTPCKaPr0,
  nSigmaTOFPiPr0,
  nSigmaTOFKaPr0,
  nSigmaTPCTOFPiPr0,
  nSigmaTPCTOFKaPr0,
  nSigmaTPCPiPr1,
  nSigmaTPCKaPr1,
  nSigmaTOFPiPr1,
  nSigmaTOFKaPr1,
  nSigmaTPCTOFPiPr1,
  nSigmaTPCTOFKaPr1,
  nSigmaTPCPiPrSoftPi,
  nSigmaTPCKaPrSoftPi,
  nSigmaTOFPiPrSoftPi,
  nSigmaTOFKaPrSoftPi,
  nSigmaTPCTOFPiPrSoftPi,
  nSigmaTPCTOFKaPrSoftPi,
};

template <typename TypeOutputScore = float>
class HfMlResponseDstarToD0Pi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseDstarToD0Pi() = default;
  /// Default destructor
  virtual ~HfMlResponseDstarToD0Pi() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Dstar candidate
  /// \param prong0 is the candidate's prong0
  /// \param prong1 is the candidate's prong1
  /// \param prongSoftPi is the candidate's prongSoftPi
  /// \return inputFeatures vector
  template <typename T1, typename T2>
  std::vector<float> getInputFeatures(T1 const& candidate,
                                      T2 const& prong0, T2 const& prong1, T2 const& prongSoftPi)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_DSTAR(chi2PCAD0);
        CHECK_AND_FILL_VEC_DSTAR(decayLengthD0);
        CHECK_AND_FILL_VEC_DSTAR(decayLengthXYD0);
        CHECK_AND_FILL_VEC_DSTAR(decayLengthNormalisedD0);
        CHECK_AND_FILL_VEC_DSTAR(decayLengthXYNormalisedD0);
        CHECK_AND_FILL_VEC_DSTAR(cpaD0);
        CHECK_AND_FILL_VEC_DSTAR(cpaXYD0);
        CHECK_AND_FILL_VEC_DSTAR(deltaIPNormalisedMaxD0);
        CHECK_AND_FILL_VEC_DSTAR(impactParameterProductD0);
        CHECK_AND_FILL_VEC_DSTAR(ptProng0);
        CHECK_AND_FILL_VEC_DSTAR(ptProng1);
        CHECK_AND_FILL_VEC_DSTAR(ptSoftPi);
        CHECK_AND_FILL_VEC_DSTAR(impactParameter0);
        CHECK_AND_FILL_VEC_DSTAR(impactParameter1);
        CHECK_AND_FILL_VEC_DSTAR(impParamSoftPi);
        CHECK_AND_FILL_VEC_DSTAR(impactParameterNormalised0);
        CHECK_AND_FILL_VEC_DSTAR(impactParameterNormalised1);
        CHECK_AND_FILL_VEC_DSTAR(normalisedImpParamSoftPi);
        CHECK_AND_FILL_VEC_DSTAR_CHARGEBASE(cosThetaStarD0, cosThetaStarD0Bar, cosThetaStarD0);
        CHECK_AND_FILL_VEC_DSTAR_CHARGEBASE(invMassD0, invMassD0Bar, massD0);
        CHECK_AND_FILL_VEC_DSTAR_DELTA_MASS_D0(deltaMassD0);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prong0, nSigmaTPCPiPr0, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prong0, nSigmaTPCKaPr0, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prong0, nSigmaTOFPiPr0, tofNSigmaPi);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prong0, nSigmaTOFKaPr0, tofNSigmaKa);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prong0, nSigmaTPCTOFPiPr0, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prong0, nSigmaTPCTOFKaPr0, tpcTofNSigmaKa);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prong1, nSigmaTPCPiPr1, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prong1, nSigmaTPCKaPr1, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prong1, nSigmaTOFPiPr1, tofNSigmaPi);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prong1, nSigmaTOFKaPr1, tofNSigmaKa);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prong1, nSigmaTPCTOFPiPr1, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prong1, nSigmaTPCTOFKaPr1, tpcTofNSigmaKa);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prongSoftPi, nSigmaTPCPiPrSoftPi, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prongSoftPi, nSigmaTPCKaPrSoftPi, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prongSoftPi, nSigmaTOFPiPrSoftPi, tofNSigmaPi);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prongSoftPi, nSigmaTOFKaPrSoftPi, tofNSigmaKa);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prongSoftPi, nSigmaTPCTOFPiPrSoftPi, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_DSTAR_FULL(prongSoftPi, nSigmaTPCTOFKaPrSoftPi, tpcTofNSigmaKa);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_DSTAR(chi2PCAD0),
      FILL_MAP_DSTAR(decayLengthD0),
      FILL_MAP_DSTAR(decayLengthXYD0),
      FILL_MAP_DSTAR(decayLengthNormalisedD0),
      FILL_MAP_DSTAR(decayLengthXYNormalisedD0),
      FILL_MAP_DSTAR(cpaD0),
      FILL_MAP_DSTAR(cpaXYD0),
      FILL_MAP_DSTAR(deltaIPNormalisedMaxD0),
      FILL_MAP_DSTAR(impactParameterProductD0),
      FILL_MAP_DSTAR(ptProng0),
      FILL_MAP_DSTAR(ptProng1),
      FILL_MAP_DSTAR(ptSoftPi),
      FILL_MAP_DSTAR(impactParameter0),
      FILL_MAP_DSTAR(impactParameter1),
      FILL_MAP_DSTAR(impParamSoftPi),
      FILL_MAP_DSTAR(impactParameterNormalised0),
      FILL_MAP_DSTAR(impactParameterNormalised1),
      FILL_MAP_DSTAR(normalisedImpParamSoftPi),
      FILL_MAP_DSTAR(cosThetaStarD0),
      FILL_MAP_DSTAR(massD0),
      FILL_MAP_DSTAR(deltaMassD0),
      FILL_MAP_DSTAR(nSigmaTPCPiPr0),
      FILL_MAP_DSTAR(nSigmaTPCKaPr0),
      FILL_MAP_DSTAR(nSigmaTOFPiPr0),
      FILL_MAP_DSTAR(nSigmaTOFKaPr0),
      FILL_MAP_DSTAR(nSigmaTPCTOFPiPr0),
      FILL_MAP_DSTAR(nSigmaTPCTOFKaPr0),
      FILL_MAP_DSTAR(nSigmaTPCPiPr1),
      FILL_MAP_DSTAR(nSigmaTPCKaPr1),
      FILL_MAP_DSTAR(nSigmaTOFPiPr1),
      FILL_MAP_DSTAR(nSigmaTOFKaPr1),
      FILL_MAP_DSTAR(nSigmaTPCTOFPiPr1),
      FILL_MAP_DSTAR(nSigmaTPCTOFKaPr1),
      FILL_MAP_DSTAR(nSigmaTPCPiPrSoftPi),
      FILL_MAP_DSTAR(nSigmaTPCKaPrSoftPi),
      FILL_MAP_DSTAR(nSigmaTOFPiPrSoftPi),
      FILL_MAP_DSTAR(nSigmaTOFKaPrSoftPi),
      FILL_MAP_DSTAR(nSigmaTPCTOFPiPrSoftPi),
      FILL_MAP_DSTAR(nSigmaTPCTOFKaPrSoftPi)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_DSTAR
#undef CHECK_AND_FILL_VEC_DSTAR_FULL
#undef CHECK_AND_FILL_VEC_DSTAR
#undef CHECK_AND_FILL_VEC_DSTAR_CHARGEBASE
#undef CHECK_AND_FILL_VEC_DSTAR_DELTA_MASS_D0

#endif // PWGHF_CORE_HFMLRESPONSEDSTARTOD0PI_H_
