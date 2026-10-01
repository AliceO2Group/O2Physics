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

/// \file HfMlResponseOmegacToOmegaPiQa.h
/// \brief Class to compute the ML response for Ωc0 → Ω- π+ (and charge-conjugate) analysis selections
/// \author Yunfan Liu <yunfan.liu@cern.ch>, China University of Geosciences
/// \author Maria Fernanda Torres Cabrera <maria.fernanda.torres.cabrera@cern.ch>, University of Houston

#ifndef PWGHF_CORE_HFMLRESPONSEOMEGACTOOMEGAPIQA_H_
#define PWGHF_CORE_HFMLRESPONSEOMEGACTOOMEGAPIQA_H_

#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

#include "Tools/ML/MlResponse.h"

#include <cstdint>
#include <map>
#include <string>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_OMEGAC0(FEATURE) \
  {                               \
    #FEATURE, static_cast<uint8_t>(InputFeaturesOmegacToOmegaPi::FEATURE)}

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
namespace o2::analysis
{
enum class InputFeaturesOmegacToOmegaPi : uint8_t {

  // Common variables
  dcaCharmBaryonDau = 0,
  dcaCascDau,
  nSigmaTPCPiFromV0,
  nSigmaTPCPiFromOmegac,
  nSigmaTPCKaFromCasc,

  // DCAFitter variables
  impactParBachFromCharmBaryonXY,
  impactParBachFromCharmBaryonZ,
  impactParCascXY,
  impactParCascZ,
  cosPACharmBaryon,
  cosPACasc,

  // KFParticle variables
  cosPaOmegacToPv,
  kfDcaXYPiFromOmegac,
  chi2TopoPiFromOmegacToPv,
  invMassCascade,
  massCascChi2OverNdf,
  kfDcaXYCascToPv,
  cosPaCascToPv,
  cosThetaStarPiFromOmegac,
  chi2NdfTopoOmegacToPv,
  ldlCasc,
  cosPaCascToOmegac,
  decayLenXYCasc,
  decayLenXYOmegac,
  ldlOmegac,
  chi2NdfTopoCascToOmegac,
  chi2NdfTopoCascToPv,
  chi2GeoOmegac,
  chi2GeoCasc,
};

template <typename TypeOutputScore = float, aod::hf_cand_casc_lf::ConstructMethod reconstructionType = aod::hf_cand_casc_lf::ConstructMethod::DcaFitter>
class HfMlResponseOmegacToOmegaPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseOmegacToOmegaPi() = default;
  /// Default destructor
  virtual ~HfMlResponseOmegacToOmegaPi() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the OMEGAC0 candidate
  /// \param lamProngPi is the candidate's lamProngPi
  /// \param cascProng is the bachelor track from the Omega decay
  /// \param charmBaryonProng is the bachelor pion from the OmegaC decay
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename T3>
  std::vector<float> getInputFeatures(T1 const& candidate, T2 const& lamProngPi, T2 const& cascProng, T3 const& charmBaryonProng)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {

      // Variables common to DCAFitter and KFParticle
      switch (idx) {

        CHECK_AND_FILL_VEC_OMEGAC0(dcaCharmBaryonDau);
        CHECK_AND_FILL_VEC_OMEGAC0(dcaCascDau);
        // TPC PID variables
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProngPi, nSigmaTPCPiFromV0, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(cascProng, nSigmaTPCKaFromCasc, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(charmBaryonProng, nSigmaTPCPiFromOmegac, tpcNSigmaPi);
      }

      // DCAFitter variables
      if constexpr (reconstructionType == aod::hf_cand_casc_lf::ConstructMethod::DcaFitter) {

        switch (idx) {

          CHECK_AND_FILL_VEC_OMEGAC0(impactParBachFromCharmBaryonXY);
          CHECK_AND_FILL_VEC_OMEGAC0(impactParBachFromCharmBaryonZ);
          CHECK_AND_FILL_VEC_OMEGAC0(impactParCascXY);
          CHECK_AND_FILL_VEC_OMEGAC0(impactParCascZ);
          CHECK_AND_FILL_VEC_OMEGAC0(cosPACharmBaryon);
          CHECK_AND_FILL_VEC_OMEGAC0(cosPACasc);
        }
      }

      // KFParticle variables
      if constexpr (reconstructionType == aod::hf_cand_casc_lf::ConstructMethod::KfParticle) {

        switch (idx) {

          CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, cosPaOmegacToPv, cosPACharmBaryon);
          CHECK_AND_FILL_VEC_OMEGAC0(kfDcaXYPiFromOmegac);
          CHECK_AND_FILL_VEC_OMEGAC0(chi2TopoPiFromOmegacToPv);
          CHECK_AND_FILL_VEC_OMEGAC0(invMassCascade);
          CHECK_AND_FILL_VEC_OMEGAC0(massCascChi2OverNdf);
          CHECK_AND_FILL_VEC_OMEGAC0(kfDcaXYCascToPv);
          CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, cosPaCascToPv, cosPACasc);
          CHECK_AND_FILL_VEC_OMEGAC0(cosThetaStarPiFromOmegac);
          CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, chi2NdfTopoOmegacToPv, chi2TopoOmegacToPv);
          CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, ldlCasc, cascldl);
          CHECK_AND_FILL_VEC_OMEGAC0(cosPaCascToOmegac);
          CHECK_AND_FILL_VEC_OMEGAC0(decayLenXYCasc);
          CHECK_AND_FILL_VEC_OMEGAC0(decayLenXYOmegac);
          CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, ldlOmegac, omegacldl);
          CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, chi2NdfTopoCascToOmegac, chi2TopoCascToOmegac);
          CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, chi2NdfTopoCascToPv, chi2TopoCascToPv);
          CHECK_AND_FILL_VEC_OMEGAC0(chi2GeoOmegac);
          CHECK_AND_FILL_VEC_OMEGAC0(chi2GeoCasc);
        }
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    // Common variables
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {

      FILL_MAP_OMEGAC0(dcaCharmBaryonDau),
      FILL_MAP_OMEGAC0(dcaCascDau),
      FILL_MAP_OMEGAC0(nSigmaTPCPiFromV0),
      FILL_MAP_OMEGAC0(nSigmaTPCPiFromOmegac),
      FILL_MAP_OMEGAC0(nSigmaTPCKaFromCasc)};

    // DCAFitter variables
    if constexpr (reconstructionType == aod::hf_cand_casc_lf::ConstructMethod::DcaFitter) {

      std::map<std::string, uint8_t> mapDcaFeatures{

        FILL_MAP_OMEGAC0(impactParBachFromCharmBaryonXY),
        FILL_MAP_OMEGAC0(impactParBachFromCharmBaryonZ),
        FILL_MAP_OMEGAC0(impactParCascXY),
        FILL_MAP_OMEGAC0(impactParCascZ),
        FILL_MAP_OMEGAC0(cosPACharmBaryon),
        FILL_MAP_OMEGAC0(cosPACasc)};

      MlResponse<TypeOutputScore>::mAvailableInputFeatures.insert(mapDcaFeatures.begin(), mapDcaFeatures.end());
    }

    // KFParticle variables
    if constexpr (reconstructionType == aod::hf_cand_casc_lf::ConstructMethod::KfParticle) {

      std::map<std::string, uint8_t> mapKfFeatures{

        FILL_MAP_OMEGAC0(cosPaOmegacToPv),
        FILL_MAP_OMEGAC0(kfDcaXYPiFromOmegac),
        FILL_MAP_OMEGAC0(chi2TopoPiFromOmegacToPv),
        FILL_MAP_OMEGAC0(invMassCascade),
        FILL_MAP_OMEGAC0(massCascChi2OverNdf),
        FILL_MAP_OMEGAC0(kfDcaXYCascToPv),
        FILL_MAP_OMEGAC0(cosPaCascToPv),
        FILL_MAP_OMEGAC0(cosThetaStarPiFromOmegac),
        FILL_MAP_OMEGAC0(chi2NdfTopoOmegacToPv),
        FILL_MAP_OMEGAC0(ldlCasc),
        FILL_MAP_OMEGAC0(cosPaCascToOmegac),
        FILL_MAP_OMEGAC0(decayLenXYCasc),
        FILL_MAP_OMEGAC0(decayLenXYOmegac),
        FILL_MAP_OMEGAC0(ldlOmegac),
        FILL_MAP_OMEGAC0(chi2NdfTopoCascToOmegac),
        FILL_MAP_OMEGAC0(chi2NdfTopoCascToPv),
        FILL_MAP_OMEGAC0(chi2GeoOmegac),
        FILL_MAP_OMEGAC0(chi2GeoCasc)};

      MlResponse<TypeOutputScore>::mAvailableInputFeatures.insert(mapKfFeatures.begin(), mapKfFeatures.end());
    }
  }
};

} // namespace o2::analysis

#undef FILL_MAP_OMEGAC0
#undef CHECK_AND_FILL_VEC_OMEGAC0_FULL
#undef CHECK_AND_FILL_VEC_OMEGAC0
#endif // PWGHF_CORE_HFMLRESPONSEOMEGACTOOMEGAPIQA_H_
