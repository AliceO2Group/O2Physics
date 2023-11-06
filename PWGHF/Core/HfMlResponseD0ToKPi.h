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

#include <map>
#include <string>
#include <vector>

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponse.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_D0(FEATURE)                                        \
  {                                                                 \
#FEATURE, static_cast < uint8_t>(InputFeaturesD0ToKPi::FEATURE) \
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
// where GETTER is a method of hfHelper
#define CHECK_AND_FILL_VEC_D0_HFHELPER(OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesD0ToKPi::FEATURE): {   \
    inputFeatures.emplace_back(hfHelper.GETTER(OBJECT));        \
    break;                                                      \
  }

namespace o2::analysis
{
enum class InputFeaturesD0ToKPi : uint8_t {
  collisionId = 0,
  posX,
  posY,
  posZ,
  xSecondaryVertex,
  ySecondaryVertex,
  zSecondaryVertex,
  errorDecayLength,
  errorDecayLengthXY,
  chi2PCA,
  rSecondaryVertex,
  decayLength,
  decayLengthXY,
  decayLengthNormalised,
  decayLengthXYNormalised,
  impactParameterNormalised0,
  ptProng0,
  pt2Prong0,
  ptProng1,
  pt2Prong1,
  pxProng0,
  pyProng0,
  pzProng0,
  pxProng1,
  pyProng1,
  pzProng1,
  impactParameter0,
  impactParameter1,
  errorImpactParameter0,
  errorImpactParameter1,
  nSigTpcPi0,
  nSigTpcKa0,
  nSigTofPi0,
  nSigTofKa0,
  nSigTpcPi1,
  nSigTpcKa1,
  nSigTofPi1,
  nSigTofKa1,
  m,
  maxNormalisedDeltaIP,
  impactParameterProduct,
  cosThetaStar,
  pt,
  p,
  cpa,
  cpaXY,
  ct,
  eta,
  phi,
  y,
  e
};

template <typename TypeOutputScore = float>
class HfMlResponseD0ToKPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseD0ToKPi() = default;
  /// Default destructor
  virtual ~HfMlResponseD0ToKPi() = default;

  HfHelper hfHelper;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Dplus candidate
  /// \param prong0 is the candidate's prong0
  /// \param prong1 is the candidate's prong1
  /// \return inputFeatures vector
  template <typename T1, typename T2>
  std::vector<float> getInputFeatures(T1 const& candidate,
                                      T2 const& prong0, T2 const& prong1)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_D0(collisionId);
        CHECK_AND_FILL_VEC_D0(posX);
        CHECK_AND_FILL_VEC_D0(posY);
        CHECK_AND_FILL_VEC_D0(posZ);
        CHECK_AND_FILL_VEC_D0(xSecondaryVertex);
        CHECK_AND_FILL_VEC_D0(ySecondaryVertex);
        CHECK_AND_FILL_VEC_D0(zSecondaryVertex);
        CHECK_AND_FILL_VEC_D0(errorDecayLength);
        CHECK_AND_FILL_VEC_D0(errorDecayLengthXY);
        CHECK_AND_FILL_VEC_D0(chi2PCA);
        CHECK_AND_FILL_VEC_D0(rSecondaryVertex);
        CHECK_AND_FILL_VEC_D0(decayLength);
        CHECK_AND_FILL_VEC_D0(decayLengthXY);
        CHECK_AND_FILL_VEC_D0(decayLengthNormalised);
        CHECK_AND_FILL_VEC_D0(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_D0(ptProng0);
        CHECK_AND_FILL_VEC_D0(ptProng1);
        CHECK_AND_FILL_VEC_D0(pxProng0);
        CHECK_AND_FILL_VEC_D0(pyProng0);
        CHECK_AND_FILL_VEC_D0(pzProng0);
        CHECK_AND_FILL_VEC_D0(pxProng1);
        CHECK_AND_FILL_VEC_D0(pyProng1);
        CHECK_AND_FILL_VEC_D0(pzProng1);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, impactParameter0, impactParameter0);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, impactParameter1, impactParameter1);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, errorImpactParameter0, errorImpactParameter0);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, errorImpactParameter1, errorImpactParameter1);
        // TPC PID variables
        CHECK_AND_FILL_VEC_D0_FULL(prong0, nSigTpcPi0, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_D0_FULL(prong0, nSigTpcKa0, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_D0_FULL(prong1, nSigTpcPi1, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_D0_FULL(prong1, nSigTpcKa1, tpcNSigmaKa);
        // TOF PID variables
        CHECK_AND_FILL_VEC_D0_FULL(prong0, nSigTofPi0, tofNSigmaPi);
        CHECK_AND_FILL_VEC_D0_FULL(prong0, nSigTofKa0, tofNSigmaKa);
        CHECK_AND_FILL_VEC_D0_FULL(prong1, nSigTofPi1, tofNSigmaPi);
        CHECK_AND_FILL_VEC_D0_FULL(prong1, nSigTofKa1, tofNSigmaKa);

        CHECK_AND_FILL_VEC_D0_HFHELPER(candidate, m, invMassD0ToPiK);
        CHECK_AND_FILL_VEC_D0(maxNormalisedDeltaIP);
        CHECK_AND_FILL_VEC_D0_FULL(candidate, impactParameterProduct, impactParameterProduct);
        CHECK_AND_FILL_VEC_D0_HFHELPER(candidate, cosThetaStar, cosThetaStarD0);
        CHECK_AND_FILL_VEC_D0(pt);
        CHECK_AND_FILL_VEC_D0(p);
        CHECK_AND_FILL_VEC_D0(cpa);
        CHECK_AND_FILL_VEC_D0(cpaXY);
        CHECK_AND_FILL_VEC_D0_HFHELPER(candidate, ct, ctD0);
        CHECK_AND_FILL_VEC_D0(eta);
        CHECK_AND_FILL_VEC_D0(phi);
        CHECK_AND_FILL_VEC_D0_HFHELPER(candidate, y, yD0);
        CHECK_AND_FILL_VEC_D0_HFHELPER(candidate, e, eD0);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_D0(collisionId),
      FILL_MAP_D0(posX),
      FILL_MAP_D0(posY),
      FILL_MAP_D0(posZ),
      FILL_MAP_D0(xSecondaryVertex),
      FILL_MAP_D0(ySecondaryVertex),
      FILL_MAP_D0(zSecondaryVertex),
      FILL_MAP_D0(errorDecayLength),
      FILL_MAP_D0(errorDecayLengthXY),
      FILL_MAP_D0(chi2PCA),
      FILL_MAP_D0(rSecondaryVertex),
      FILL_MAP_D0(decayLength),
      FILL_MAP_D0(decayLengthXY),
      FILL_MAP_D0(decayLengthNormalised),
      FILL_MAP_D0(decayLengthXYNormalised),
      FILL_MAP_D0(ptProng0),
      FILL_MAP_D0(ptProng1),
      FILL_MAP_D0(pxProng0),
      FILL_MAP_D0(pyProng0),
      FILL_MAP_D0(pzProng0),
      FILL_MAP_D0(pxProng1),
      FILL_MAP_D0(pyProng1),
      FILL_MAP_D0(pzProng1),
      FILL_MAP_D0(impactParameter0),
      FILL_MAP_D0(impactParameter1),
      FILL_MAP_D0(errorImpactParameter0),
      FILL_MAP_D0(errorImpactParameter1),
      // TPC PID variables
      FILL_MAP_D0(nSigTpcPi0),
      FILL_MAP_D0(nSigTpcKa0),
      FILL_MAP_D0(nSigTpcPi1),
      FILL_MAP_D0(nSigTpcKa1),
      // TOF PID variables
      FILL_MAP_D0(nSigTofPi0),
      FILL_MAP_D0(nSigTofKa0),
      FILL_MAP_D0(nSigTofPi1),
      FILL_MAP_D0(nSigTofKa1),

      FILL_MAP_D0(m),
      FILL_MAP_D0(maxNormalisedDeltaIP),
      FILL_MAP_D0(impactParameterProduct),
      FILL_MAP_D0(cosThetaStar),
      FILL_MAP_D0(pt),
      FILL_MAP_D0(p),
      FILL_MAP_D0(cpa),
      FILL_MAP_D0(cpaXY),
      FILL_MAP_D0(ct),
      FILL_MAP_D0(eta),
      FILL_MAP_D0(phi),
      FILL_MAP_D0(y),
      FILL_MAP_D0(e)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_D0
#undef CHECK_AND_FILL_VEC_D0_FULL
#undef CHECK_AND_FILL_VEC_D0
#undef CHECK_AND_FILL_VEC_D0_HFHELPER

#endif // PWGHF_CORE_HFMLRESPONSED0TOKPI_H_