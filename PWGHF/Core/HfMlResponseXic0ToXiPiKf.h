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

#include "Common/Core/RecoDecay.h"
#include "Tools/ML/MlResponse.h"

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_XIC0TOXIPIKF(FEATURE) \
  {                                    \
    #FEATURE, static_cast<uint8_t>(InputFeaturesXic0ToXiPiKf::FEATURE)}

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

#define CHECK_AND_FILL_VEC_XIC0TOXIPIKF_VALUE(FEATURE, VALUE)      \
  case static_cast<uint8_t>(InputFeaturesXic0ToXiPiKf::FEATURE): { \
    inputFeatures.emplace_back(VALUE);                             \
    break;                                                         \
  }

#define FILL_MAP_XIC0TOXIPIKF_AS(NAME, FEATURE) \
  {NAME, static_cast<uint8_t>(InputFeaturesXic0ToXiPiKf::FEATURE)}

namespace o2::analysis
{

enum class InputFeaturesXic0ToXiPiKf : uint8_t {
  // Existing feature keys retained for backwards compatibility.
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
  decayLenXYCasc,
  // Remaining features used by the local PbPb Transformer.
  decayLenXYXic,
  chi2TopoXicToPv,
  chi2TopoPiFromXicToPv,
  kfptPiFromXic,
  cosPaCascToPv,
  dcaXYToPvCascDau,
  decayLenXYLambda,
  v0Chi2OverNdf,
  chi2TopoV0ToPv,
  chi2TopoV0ToCasc,
  cosPaV0ToPv,
  cosPaV0ToCasc,
  dcaV0Dau,
  dcaXYToPvV0Dau0,
  dcaXYToPvV0Dau1,
  tpcNSigmaPrFromLambda
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
  std::vector<float> getInputFeatures(T1 const& candidate, T2 const& lamProngPi, T2 const& lamProngPr,
                                      T2 const& cascProngPi, T3 const& charmBaryonProngPi)
  {
    std::vector<float> inputFeatures;
    inputFeatures.reserve(MlResponse<TypeOutputScore>::mCachedIndices.size());

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        // Xic candidate topology.
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(decayLenXYXic);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(xicChi2OverNdf);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(chi2TopoXicToPv);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(dcaCharmBaryonDau);
        // Bachelor pion from Xic.
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(kfDcaXYPiFromXic);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(chi2TopoPiFromXicToPv);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF_VALUE(kfptPiFromXic, RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCharmBaryon(), candidate.pyBachFromCharmBaryon()));
        // Cascade.
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(decayLenXYCasc);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(cascChi2OverNdf);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(chi2TopoCascToPv);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(chi2TopoCascToXic);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF_FULL(candidate, cosPaCascToPv, cosPACasc);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(cosPaCascToXic);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(dcaCascDau);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(kfDcaXYCascToPv);
        // Pion from cascade.
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF_FULL(cascProngPi, tpcNSigmaPiFromCasc, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(dcaXYToPvCascDau);
        // Lambda and its charged daughters. The selector passes sign-corrected prongs.
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(decayLenXYLambda);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(v0Chi2OverNdf);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(chi2TopoV0ToPv);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(chi2TopoV0ToCasc);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF_FULL(candidate, cosPaV0ToPv, cosPAV0);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(cosPaV0ToCasc);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(dcaV0Dau);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(dcaXYToPvV0Dau0);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(dcaXYToPvV0Dau1);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF_FULL(lamProngPr, tpcNSigmaPrFromLambda, tpcNSigmaPr);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF_FULL(lamProngPi, tpcNSigmaPiFromLambda, tpcNSigmaPi);
        // Legacy features that are not part of the local 28-feature model.
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF_FULL(charmBaryonProngPi, tpcNSigmaPiFromCharmBaryon, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_XIC0TOXIPIKF(cascldl);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_XIC0TOXIPIKF_AS("fDecayLenXYXic", decayLenXYXic),
      FILL_MAP_XIC0TOXIPIKF(xicChi2OverNdf),
      FILL_MAP_XIC0TOXIPIKF_AS("fXicChi2OverNdf", xicChi2OverNdf),
      FILL_MAP_XIC0TOXIPIKF_AS("fChi2TopoXicToPv", chi2TopoXicToPv),
      FILL_MAP_XIC0TOXIPIKF(dcaCharmBaryonDau),
      FILL_MAP_XIC0TOXIPIKF_AS("fDcaCharmBaryonDau", dcaCharmBaryonDau),
      FILL_MAP_XIC0TOXIPIKF(kfDcaXYPiFromXic),
      FILL_MAP_XIC0TOXIPIKF_AS("fKfDcaXYPiFromXic", kfDcaXYPiFromXic),
      FILL_MAP_XIC0TOXIPIKF_AS("fChi2TopoPiFromXicToPv", chi2TopoPiFromXicToPv),
      FILL_MAP_XIC0TOXIPIKF_AS("fKfptPiFromXic", kfptPiFromXic),
      FILL_MAP_XIC0TOXIPIKF(decayLenXYCasc),
      FILL_MAP_XIC0TOXIPIKF_AS("fDecayLenXYCasc", decayLenXYCasc),
      FILL_MAP_XIC0TOXIPIKF(cascChi2OverNdf),
      FILL_MAP_XIC0TOXIPIKF_AS("fCascChi2OverNdf", cascChi2OverNdf),
      FILL_MAP_XIC0TOXIPIKF(chi2TopoCascToPv),
      FILL_MAP_XIC0TOXIPIKF_AS("fChi2TopoCascToPv", chi2TopoCascToPv),
      FILL_MAP_XIC0TOXIPIKF(chi2TopoCascToXic),
      FILL_MAP_XIC0TOXIPIKF_AS("fChi2TopoCascToXic", chi2TopoCascToXic),
      FILL_MAP_XIC0TOXIPIKF_AS("fCosPaCascToPv", cosPaCascToPv),
      FILL_MAP_XIC0TOXIPIKF(cosPaCascToXic),
      FILL_MAP_XIC0TOXIPIKF_AS("fCosPaCascToXic", cosPaCascToXic),
      FILL_MAP_XIC0TOXIPIKF(dcaCascDau),
      FILL_MAP_XIC0TOXIPIKF_AS("fDcaCascDau", dcaCascDau),
      FILL_MAP_XIC0TOXIPIKF(kfDcaXYCascToPv),
      FILL_MAP_XIC0TOXIPIKF_AS("fKfDcaXYCascToPv", kfDcaXYCascToPv),
      FILL_MAP_XIC0TOXIPIKF(tpcNSigmaPiFromCasc),
      FILL_MAP_XIC0TOXIPIKF_AS("fTpcNSigmaPiFromCasc", tpcNSigmaPiFromCasc),
      FILL_MAP_XIC0TOXIPIKF_AS("fDcaXYToPvCascDau", dcaXYToPvCascDau),
      FILL_MAP_XIC0TOXIPIKF_AS("fDecayLenXYLambda", decayLenXYLambda),
      FILL_MAP_XIC0TOXIPIKF_AS("fV0Chi2OverNdf", v0Chi2OverNdf),
      FILL_MAP_XIC0TOXIPIKF_AS("fChi2TopoV0ToPv", chi2TopoV0ToPv),
      FILL_MAP_XIC0TOXIPIKF_AS("fChi2TopoV0ToCasc", chi2TopoV0ToCasc),
      FILL_MAP_XIC0TOXIPIKF_AS("fCosPaV0ToPv", cosPaV0ToPv),
      FILL_MAP_XIC0TOXIPIKF_AS("fCosPaV0ToCasc", cosPaV0ToCasc),
      FILL_MAP_XIC0TOXIPIKF_AS("fDcaV0Dau", dcaV0Dau),
      FILL_MAP_XIC0TOXIPIKF_AS("fDcaXYToPvV0Dau0", dcaXYToPvV0Dau0),
      FILL_MAP_XIC0TOXIPIKF_AS("fDcaXYToPvV0Dau1", dcaXYToPvV0Dau1),
      FILL_MAP_XIC0TOXIPIKF_AS("fTpcNSigmaPrFromLambda", tpcNSigmaPrFromLambda),
      FILL_MAP_XIC0TOXIPIKF(tpcNSigmaPiFromLambda),
      FILL_MAP_XIC0TOXIPIKF_AS("fTpcNSigmaPiFromLambda", tpcNSigmaPiFromLambda),
      FILL_MAP_XIC0TOXIPIKF(tpcNSigmaPiFromCharmBaryon),
      FILL_MAP_XIC0TOXIPIKF(cascldl),
    };
  }
};

} // namespace o2::analysis

#undef FILL_MAP_XIC0TOXIPIKF
#undef FILL_MAP_XIC0TOXIPIKF_AS
#undef CHECK_AND_FILL_VEC_XIC0TOXIPIKF_FULL
#undef CHECK_AND_FILL_VEC_XIC0TOXIPIKF
#undef CHECK_AND_FILL_VEC_XIC0TOXIPIKF_VALUE

#endif // PWGHF_CORE_HFMLRESPONSEXIC0TOXIPIKF_H_
