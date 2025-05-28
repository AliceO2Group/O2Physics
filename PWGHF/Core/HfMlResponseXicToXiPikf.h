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

/// \file HfMlResponseXicToXiPikf.h
/// \brief Class to compute the ML response for Ξc^0 → Ξ∓ π± analysis selections
/// \author Tao Fang <tao.fang@cern.ch>, Central China Normal University

#ifndef PWGHF_CORE_HFMLRESPONSEXICTOXIPIKF_H_
#define PWGHF_CORE_HFMLRESPONSEXICTOXIPIKF_H_

#include <map>
#include <string>
#include <vector>

#include "PWGHF/Core/HfMlResponse.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_XICTOXIPI(FEATURE)                                 \
  {                                                                 \
    #FEATURE, static_cast<uint8_t>(InputFeaturesXicToXiPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_XICTOXIPI_FULL(OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesXicToXiPi::FEATURE): {    \
    inputFeatures.emplace_back(OBJECT.GETTER());                   \
    break;                                                         \
  }

// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_XICTOXIPI(GETTER)                   \
  case static_cast<uint8_t>(InputFeaturesXicToXiPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());            \
    break;                                                     \
  }

namespace o2::analysis
{

enum class InputFeaturesXicToXiPi : uint8_t {
  tpcNSigmaPiFromLambda,
  tpcNSigmaPiFromCasc,
  tpcNSigmaPiFromCharmBaryon,
  dcaCascDau,
  dcaCharmBaryonDau,
  kfDcaXYPiFromXic,
  kfDcaXYCascToPv,
  cascChi2OverNdf,
  xicChi2OverNdf,
  cascldl,
  chi2TopoCascToPv,
  chi2TopoCascToXic,
  cosPaCascToXic,
  decayLenXYCasc
};

template <typename TypeOutputScore = float>
class HfMlResponseXicToXiPikf : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseXicToXiPikf() = default;
  /// Default destructor
  virtual ~HfMlResponseXicToXiPikf() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Xic candidate
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename T3>
  // std::vector<float> getInputFeatures(T1 const& candidate)
  std::vector<float> getInputFeatures(T1 const& candidate, T2 const& lamProngPi, T2 const& cascProngPi, T3 const& charmBaryonProngPi)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        // PID variables
        CHECK_AND_FILL_VEC_XICTOXIPI_FULL(lamProngPi, tpcNSigmaPiFromLambda, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_XICTOXIPI_FULL(cascProngPi, tpcNSigmaPiFromCasc, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_XICTOXIPI_FULL(charmBaryonProngPi, tpcNSigmaPiFromCharmBaryon, tpcNSigmaPi);
        // DCA
        CHECK_AND_FILL_VEC_XICTOXIPI(dcaCascDau)
        CHECK_AND_FILL_VEC_XICTOXIPI(dcaCharmBaryonDau);
        CHECK_AND_FILL_VEC_XICTOXIPI(kfDcaXYPiFromXic);
        CHECK_AND_FILL_VEC_XICTOXIPI(kfDcaXYCascToPv);
        // Chi2Geo
        CHECK_AND_FILL_VEC_XICTOXIPI(cascChi2OverNdf);
        CHECK_AND_FILL_VEC_XICTOXIPI(xicChi2OverNdf);
        // ldl
        CHECK_AND_FILL_VEC_XICTOXIPI(cascldl);
        // Chi2Topo
        CHECK_AND_FILL_VEC_XICTOXIPI(chi2TopoCascToPv);
        CHECK_AND_FILL_VEC_XICTOXIPI(chi2TopoCascToXic);
        // CosPa
        CHECK_AND_FILL_VEC_XICTOXIPI(cosPaCascToXic);
        // Decay length
        CHECK_AND_FILL_VEC_XICTOXIPI(decayLenXYCasc);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_XICTOXIPI(tpcNSigmaPiFromLambda),
      FILL_MAP_XICTOXIPI(tpcNSigmaPiFromCasc),
      FILL_MAP_XICTOXIPI(tpcNSigmaPiFromCharmBaryon),
      FILL_MAP_XICTOXIPI(dcaCascDau),
      FILL_MAP_XICTOXIPI(dcaCharmBaryonDau),
      FILL_MAP_XICTOXIPI(kfDcaXYPiFromXic),
      FILL_MAP_XICTOXIPI(kfDcaXYCascToPv),
      FILL_MAP_XICTOXIPI(cascChi2OverNdf),
      FILL_MAP_XICTOXIPI(xicChi2OverNdf),
      FILL_MAP_XICTOXIPI(cascldl),
      FILL_MAP_XICTOXIPI(chi2TopoCascToPv),
      FILL_MAP_XICTOXIPI(chi2TopoCascToXic),
      FILL_MAP_XICTOXIPI(cosPaCascToXic),
      FILL_MAP_XICTOXIPI(decayLenXYCasc),
    };
  }
};

} // namespace o2::analysis

#undef FILL_MAP_XICTOXIPI
#undef CHECK_AND_FILL_VEC_XICTOXIPI_FULL
#undef CHECK_AND_FILL_VEC_XICTOXIPI

#endif // PWGHF_CORE_HFMLRESPONSEXICTOXIPIKF_H_
