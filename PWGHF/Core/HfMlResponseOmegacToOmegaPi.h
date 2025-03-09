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
/// \brief Class to compute the ML response for Ωc± → Ω∓ π±  analysis selections
/// \author Yunfan Liu <yunfan.liu@cern.ch>, China University of Geosciences

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
#define FILL_MAP_OMEGAC0(FEATURE)                                         \
  {                                                                       \
    #FEATURE, static_cast<uint8_t>(InputFeaturesOmegacToOmegaPi::FEATURE) \
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
namespace o2::analysis
{
enum class InputFeaturesOmegacToOmegaPi : uint8_t {

  cosPaOmegacToPv = 0,
  kfDcaXYPiFromOmegac,
  cosThetaStarPiFromOmegac,
  chi2TopoPiFromOmegacToPv,
  dcaCharmBaryonDau,
  invMassCascade,
  massCascChi2OverNdf,
  cosPaCascToPv,
  kfDcaXYCascToPv,
  nSigmaTPCPiFromV0,
  nSigmaTPCPiFromOmegac,
  nSigmaTPCKaFromCasc

};

template <typename TypeOutputScore = float>
class HfMlResponseOmegacToOmegaPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseOmegacToOmegaPi() = default;
  /// Default destructor
  virtual ~HfMlResponseOmegacToOmegaPi() = default;

  HfHelper hfHelper;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the OMEGAC0 candidate
  /// \param lamProngPi is the candidate's lamProngPi
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename T3>
  std::vector<float> getInputFeatures(T1 const& candidate, T2 const& lamProngPi, T2 const& cascProng, T3 const& charmBaryonProng)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {

        CHECK_AND_FILL_VEC_OMEGAC0(kfDcaXYPiFromOmegac);
        CHECK_AND_FILL_VEC_OMEGAC0(cosThetaStarPiFromOmegac);
        CHECK_AND_FILL_VEC_OMEGAC0(chi2TopoPiFromOmegacToPv);
        CHECK_AND_FILL_VEC_OMEGAC0(dcaCharmBaryonDau);
        CHECK_AND_FILL_VEC_OMEGAC0(invMassCascade);
        CHECK_AND_FILL_VEC_OMEGAC0(massCascChi2OverNdf);
        CHECK_AND_FILL_VEC_OMEGAC0(kfDcaXYCascToPv);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, cosPaOmegacToPv, cosPACharmBaryon);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, cosPaCascToPv, cosPACasc);
        // TPC PID variables
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProngPi, nSigmaTPCPiFromV0, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(cascProng, nSigmaTPCKaFromCasc, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(charmBaryonProng, nSigmaTPCPiFromOmegac, tpcNSigmaPi);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {

      FILL_MAP_OMEGAC0(invMassCascade),
      FILL_MAP_OMEGAC0(cosPaOmegacToPv),
      FILL_MAP_OMEGAC0(dcaCharmBaryonDau),
      FILL_MAP_OMEGAC0(kfDcaXYPiFromOmegac),
      FILL_MAP_OMEGAC0(cosThetaStarPiFromOmegac),
      FILL_MAP_OMEGAC0(chi2TopoPiFromOmegacToPv),
      FILL_MAP_OMEGAC0(massCascChi2OverNdf),
      FILL_MAP_OMEGAC0(cosPaCascToPv),
      FILL_MAP_OMEGAC0(kfDcaXYCascToPv),
      // TPC PID variables
      FILL_MAP_OMEGAC0(nSigmaTPCPiFromV0),
      FILL_MAP_OMEGAC0(nSigmaTPCKaFromCasc),
      FILL_MAP_OMEGAC0(nSigmaTPCPiFromOmegac),

    };
  }
};

} // namespace o2::analysis

#undef FILL_MAP_OMEGAC0
#undef CHECK_AND_FILL_VEC_OMEGAC0_FULL
#undef CHECK_AND_FILL_VEC_OMEGAC0
#undef CHECK_AND_FILL_VEC_OMEGAC0_HFHELPER
#endif // PWGHF_CORE_HFMLRESPONSEOMEGACTOOMEGAPI_H_
