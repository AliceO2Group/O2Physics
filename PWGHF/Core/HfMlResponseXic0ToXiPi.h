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

/// \file HfMlResponseXic0ToXiPi.h
/// \brief Class to compute the ML response for Ξc^0 → Ξ∓ π± analysis selections
/// \author Tao Fang <tao.fang@cern.ch>, Central China Normal University

#ifndef PWGHF_CORE_HFMLRESPONSEXIC0TOXIPI_H_
#define PWGHF_CORE_HFMLRESPONSEXIC0TOXIPI_H_

#include "PWGHF/Core/HfMlResponse.h"

#include "Tools/ML/MlResponse.h"

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_XIC0TOXIPI(FEATURE)                                 \
  {                                                                  \
    #FEATURE, static_cast<uint8_t>(InputFeaturesXic0ToXiPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_XIC0TOXIPI_FULL(OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesXic0ToXiPi::FEATURE): {    \
    inputFeatures.emplace_back(OBJECT.GETTER());                    \
    break;                                                          \
  }

// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_XIC0TOXIPI(GETTER)                   \
  case static_cast<uint8_t>(InputFeaturesXic0ToXiPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());             \
    break;                                                      \
  }

namespace o2::analysis
{

enum class InputFeaturesXic0ToXiPi : uint8_t {
  tpcNSigmaPiFromLambda,
  tpcNSigmaPiFromCasc,
  tpcNSigmaPiFromCharmBaryon,
  dcaCascDau,
  dcaCharmBaryonDau,
  cosPACharmBaryon,
  cosPACasc,
  cosPAV0,
  impactParBachFromCharmBaryonXY,
  impactParCascXY
};

template <typename TypeOutputScore = float>
class HfMlResponseXic0ToXiPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseXic0ToXiPi() = default;
  /// Default destructor
  virtual ~HfMlResponseXic0ToXiPi() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Xic0 candidate
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename T3>
  // std::vector<float> getInputFeatures(T1 const& candidate)
  std::vector<float> getInputFeatures(T1 const& candidate, T2 const& lamProngPi, T2 const& cascProngPi, T3 const& charmBaryonProngPi)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        // PID variables
        CHECK_AND_FILL_VEC_XIC0TOXIPI_FULL(lamProngPi, tpcNSigmaPiFromLambda, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_XIC0TOXIPI_FULL(cascProngPi, tpcNSigmaPiFromCasc, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_XIC0TOXIPI_FULL(charmBaryonProngPi, tpcNSigmaPiFromCharmBaryon, tpcNSigmaPi);
        // DCA
        CHECK_AND_FILL_VEC_XIC0TOXIPI(dcaCascDau);
        CHECK_AND_FILL_VEC_XIC0TOXIPI(dcaCharmBaryonDau);
        // CosPA
        CHECK_AND_FILL_VEC_XIC0TOXIPI(cosPACharmBaryon);
        CHECK_AND_FILL_VEC_XIC0TOXIPI(cosPACasc);
        CHECK_AND_FILL_VEC_XIC0TOXIPI(cosPAV0);
        // ImpactPar
        CHECK_AND_FILL_VEC_XIC0TOXIPI(impactParBachFromCharmBaryonXY);
        CHECK_AND_FILL_VEC_XIC0TOXIPI(impactParCascXY);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_XIC0TOXIPI(tpcNSigmaPiFromLambda),
      FILL_MAP_XIC0TOXIPI(tpcNSigmaPiFromCasc),
      FILL_MAP_XIC0TOXIPI(tpcNSigmaPiFromCharmBaryon),
      FILL_MAP_XIC0TOXIPI(dcaCascDau),
      FILL_MAP_XIC0TOXIPI(dcaCharmBaryonDau),
      FILL_MAP_XIC0TOXIPI(cosPACharmBaryon),
      FILL_MAP_XIC0TOXIPI(cosPACasc),
      FILL_MAP_XIC0TOXIPI(cosPAV0),
      FILL_MAP_XIC0TOXIPI(impactParBachFromCharmBaryonXY),
      FILL_MAP_XIC0TOXIPI(impactParCascXY)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_XIC0TOXIPI
#undef CHECK_AND_FILL_VEC_XIC0TOXIPI_FULL
#undef CHECK_AND_FILL_VEC_XIC0TOXIPI

#endif // PWGHF_CORE_HFMLRESPONSEXIC0TOXIPI_H_
