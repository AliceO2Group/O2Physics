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

/// \file HfMlResponseXic0ToXiPiKf.h
/// \brief Class to compute the ML response for Ξc^0 → Ξ∓ π± kf analysis selections
/// \author Tao Fang <tao.fang@cern.ch>, Central China Normal University

#ifndef PWGHF_CORE_HFMLRESPONSEXIC0TOXIPIKF_H_
#define PWGHF_CORE_HFMLRESPONSEXIC0TOXIPIKF_H_

#include "PWGHF/Core/HfMlResponse.h"

#include "Tools/ML/MlResponse.h"

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_XIC0TOXIPIKF(FEATURE)                                 \
  {                                                                    \
    #FEATURE, static_cast<uint8_t>(InputFeaturesXic0ToXiPiKf::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_XIC0TOXIPIKF_FULL(OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesXic0ToXiPiKf::FEATURE): {    \
    inputFeatures.emplace_back(OBJECT.GETTER());                      \
    break;                                                            \
  }

// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_XIC0TOXIPIKF(GETTER)                   \
  case static_cast<uint8_t>(InputFeaturesXic0ToXiPiKf::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());               \
    break;                                                        \
  }

namespace o2::analysis
{

enum class InputFeaturesXic0ToXiPiKf : uint8_t {
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
class HfMlResponseXic0ToXiPiKf : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseXic0ToXiPiKf() = default;
  /// Default destructor
  virtual ~HfMlResponseXic0ToXiPiKf() = default;

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
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF_FULL(lamProngPi, tpcNSigmaPiFromLambda, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF_FULL(cascProngPi, tpcNSigmaPiFromCasc, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF_FULL(charmBaryonProngPi, tpcNSigmaPiFromCharmBaryon, tpcNSigmaPi);
        // DCA
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(dcaCascDau);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(dcaCharmBaryonDau);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(kfDcaXYPiFromXic);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(kfDcaXYCascToPv);
        // Chi2Geo
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(cascChi2OverNdf);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(xicChi2OverNdf);
        // ldl
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(cascldl);
        // Chi2Topo
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(chi2TopoCascToPv);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(chi2TopoCascToXic);
        // CosPa
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(cosPaCascToXic);
        // Decay length
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(decayLenXYCasc);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_XIC0TOXIPIKF(tpcNSigmaPiFromLambda),
      FILL_MAP_XIC0TOXIPIKF(tpcNSigmaPiFromCasc),
      FILL_MAP_XIC0TOXIPIKF(tpcNSigmaPiFromCharmBaryon),
      FILL_MAP_XIC0TOXIPIKF(dcaCascDau),
      FILL_MAP_XIC0TOXIPIKF(dcaCharmBaryonDau),
      FILL_MAP_XIC0TOXIPIKF(kfDcaXYPiFromXic),
      FILL_MAP_XIC0TOXIPIKF(kfDcaXYCascToPv),
      FILL_MAP_XIC0TOXIPIKF(cascChi2OverNdf),
      FILL_MAP_XIC0TOXIPIKF(xicChi2OverNdf),
      FILL_MAP_XIC0TOXIPIKF(cascldl),
      FILL_MAP_XIC0TOXIPIKF(chi2TopoCascToPv),
      FILL_MAP_XIC0TOXIPIKF(chi2TopoCascToXic),
      FILL_MAP_XIC0TOXIPIKF(cosPaCascToXic),
      FILL_MAP_XIC0TOXIPIKF(decayLenXYCasc),
    };
  }
};

} // namespace o2::analysis

#undef FILL_MAP_XIC0TOXIPIKF
#undef CHECK_AND_FILL_VEC_XIC0TOXIPIKF_FULL
#undef CHECK_AND_FILL_VEC_XIC0TOXIPIKF

#endif // PWGHF_CORE_HFMLRESPONSEXIC0TOXIPIKF_H_
