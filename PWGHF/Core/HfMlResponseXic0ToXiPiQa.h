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

/// \file HfMlResponseXic0ToXiPiQa.h
/// \brief Class to compute the ML response for XiC0 to Xi Pi selections
/// \author Fabio Catalano <fabio.catalano@cern.ch>, University of Houston

#ifndef PWGHF_CORE_HFMLRESPONSEXIC0TOXIPIQA_H_
#define PWGHF_CORE_HFMLRESPONSEXIC0TOXIPIQA_H_

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
#define FILL_MAP_XIC0(FEATURE) \
  {                            \
    #FEATURE, static_cast<uint8_t>(InputFeaturesXic0ToXiPiQa::FEATURE)}

// Specific case of CHECK_AND_FILL_VEC_XIC0_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_XIC0(GETTER)                           \
  case static_cast<uint8_t>(InputFeaturesXic0ToXiPiQa::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());               \
    break;                                                        \
  }

// Check whether the cached feature index matches FEATURE, then append
// OBJECT.GETTER() to the input feature vector
#define CHECK_AND_FILL_VEC_XIC0_FULL(OBJECT, FEATURE, GETTER)      \
  case static_cast<uint8_t>(InputFeaturesXic0ToXiPiQa::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());                   \
    break;                                                         \
  }

// Variation of CHECK_AND_FILL_VEC_XIC0_FULL for HfHelper methods
#define CHECK_AND_FILL_VEC_XIC0_HFHELPER(OBJECT, FEATURE, GETTER)  \
  case static_cast<uint8_t>(InputFeaturesXic0ToXiPiQa::FEATURE): { \
    inputFeatures.emplace_back(HfHelper::GETTER(OBJECT));          \
    break;                                                         \
  }

namespace o2::analysis
{
enum class InputFeaturesXic0ToXiPiQa : uint8_t {
  // Common variables
  impactParameterXYCharmBachelor = 0,
  impactParameterZCharmBachelor,
  impactParameterXYCascade,
  impactParameterZCascade,
  impactParameterXYCascadeBach,
  impactParameterXYV0PosDau,
  impactParameterXYV0NegDau,
  dcaDaughters,
  dcaCascadeDaughters,
  dcaV0Daughters,
  cpa,
  cpaCascade,
  cpaV0,
  cpaCascadeToCharmHadron,
  cpaV0ToCascade,
  cosThetaStarPiFromXic0,
  decayLengthXY,
  decayLengthXYCascade,
  decayLengthXYV0,
  nSigmaTPCPiFromV0,
  nSigmaTPCPrFromV0,
  nSigmaTPCPiFromXic0,
  nSigmaTPCPiFromCasc,

  // KFParticle variables
  decayLengthNormKf,
  decayLengthNormKfCascade,
  chi2GeoOverNdf,
  chi2TopoOverNdfToPv,
  chi2GeoOverNdfCascade,
  chi2TopoOverNdfCascadeToCharmHadron,
  chi2TopoOverNdfCascadeToPv,
  chi2GeoOverNdfV0,
  chi2TopoOverNdfV0ToPv,
  chi2TopoOverNdfV0ToCascade

  // Legacy KF inputs not stored in the shared QA candidate tables
  // chi2TopoPiFromXicToPv,
  // kfptPiFromXic
};

template <typename TypeOutputScore = float, aod::hf_cand_casc_lf::ConstructMethod reconstructionType = aod::hf_cand_casc_lf::ConstructMethod::DcaFitter>
class HfMlResponseXic0ToXiPiQa : public HfMlResponse<TypeOutputScore>
{
 public:
  HfMlResponseXic0ToXiPiQa() = default;
  virtual ~HfMlResponseXic0ToXiPiQa() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the XiC0 candidate
  /// \param lamProngPi is the candidate's pion from the Lambda decay
  /// \param lamProngPr is the candidate's proton from the Lambda decay
  /// \param cascProng is the bachelor track from the Xi decay
  /// \param charmBaryonProng is the bachelor pion from the XiC decay
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename T3>
  std::vector<float> getInputFeatures(T1 const& candidate, T2 const& lamProngPi, T2 const& lamProngPr,
                                      T2 const& cascProng, T3 const& charmBaryonProng)
  {
    std::vector<float> inputFeatures;
    inputFeatures.reserve(MlResponse<TypeOutputScore>::mCachedIndices.size());

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {

      // Variables common to DCAFitter and KFParticle
      switch (idx) {
        CHECK_AND_FILL_VEC_XIC0_FULL(candidate, impactParameterXYCharmBachelor, impactParameter1);
        CHECK_AND_FILL_VEC_XIC0_FULL(candidate, impactParameterZCharmBachelor, impactParameterZ1);
        CHECK_AND_FILL_VEC_XIC0_FULL(candidate, impactParameterXYCascade, impactParameter0);
        CHECK_AND_FILL_VEC_XIC0_FULL(candidate, impactParameterZCascade, impactParameterZ0);
        CHECK_AND_FILL_VEC_XIC0(impactParameterXYCascadeBach);
        CHECK_AND_FILL_VEC_XIC0(impactParameterXYV0PosDau);
        CHECK_AND_FILL_VEC_XIC0(impactParameterXYV0NegDau);
        CHECK_AND_FILL_VEC_XIC0(dcaDaughters);
        CHECK_AND_FILL_VEC_XIC0(dcaCascadeDaughters);
        CHECK_AND_FILL_VEC_XIC0(dcaV0Daughters);
        CHECK_AND_FILL_VEC_XIC0(cpa);
        CHECK_AND_FILL_VEC_XIC0(cpaCascade);
        CHECK_AND_FILL_VEC_XIC0(cpaV0);
        CHECK_AND_FILL_VEC_XIC0(cpaCascadeToCharmHadron);
        CHECK_AND_FILL_VEC_XIC0(cpaV0ToCascade);
        CHECK_AND_FILL_VEC_XIC0_HFHELPER(candidate, cosThetaStarPiFromXic0, cosThetaStarXiczeroToXiPi);
        CHECK_AND_FILL_VEC_XIC0(decayLengthXY);
        CHECK_AND_FILL_VEC_XIC0(decayLengthXYCascade);
        CHECK_AND_FILL_VEC_XIC0(decayLengthXYV0);
        CHECK_AND_FILL_VEC_XIC0_FULL(lamProngPi, nSigmaTPCPiFromV0, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_XIC0_FULL(lamProngPr, nSigmaTPCPrFromV0, tpcNSigmaPr);
        CHECK_AND_FILL_VEC_XIC0_FULL(charmBaryonProng, nSigmaTPCPiFromXic0, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_XIC0_FULL(cascProng, nSigmaTPCPiFromCasc, tpcNSigmaPi);
      }

      // KFParticle variables
      if constexpr (reconstructionType == aod::hf_cand_casc_lf::ConstructMethod::KfParticle) {
        switch (idx) {
          CHECK_AND_FILL_VEC_XIC0(decayLengthNormKf);
          CHECK_AND_FILL_VEC_XIC0(decayLengthNormKfCascade);
          CHECK_AND_FILL_VEC_XIC0(chi2GeoOverNdf);
          CHECK_AND_FILL_VEC_XIC0(chi2TopoOverNdfToPv);
          CHECK_AND_FILL_VEC_XIC0(chi2GeoOverNdfCascade);
          CHECK_AND_FILL_VEC_XIC0(chi2TopoOverNdfCascadeToCharmHadron);
          CHECK_AND_FILL_VEC_XIC0(chi2TopoOverNdfCascadeToPv);
          CHECK_AND_FILL_VEC_XIC0(chi2GeoOverNdfV0);
          CHECK_AND_FILL_VEC_XIC0(chi2TopoOverNdfV0ToPv);
          CHECK_AND_FILL_VEC_XIC0(chi2TopoOverNdfV0ToCascade);
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
      FILL_MAP_XIC0(impactParameterXYCharmBachelor),
      FILL_MAP_XIC0(impactParameterZCharmBachelor),
      FILL_MAP_XIC0(impactParameterXYCascade),
      FILL_MAP_XIC0(impactParameterZCascade),
      FILL_MAP_XIC0(impactParameterXYCascadeBach),
      FILL_MAP_XIC0(impactParameterXYV0PosDau),
      FILL_MAP_XIC0(impactParameterXYV0NegDau),
      FILL_MAP_XIC0(dcaDaughters),
      FILL_MAP_XIC0(dcaCascadeDaughters),
      FILL_MAP_XIC0(dcaV0Daughters),
      FILL_MAP_XIC0(cpa),
      FILL_MAP_XIC0(cpaCascade),
      FILL_MAP_XIC0(cpaV0),
      FILL_MAP_XIC0(cpaCascadeToCharmHadron),
      FILL_MAP_XIC0(cpaV0ToCascade),
      FILL_MAP_XIC0(cosThetaStarPiFromXic0),
      FILL_MAP_XIC0(decayLengthXY),
      FILL_MAP_XIC0(decayLengthXYCascade),
      FILL_MAP_XIC0(decayLengthXYV0),
      FILL_MAP_XIC0(nSigmaTPCPiFromV0),
      FILL_MAP_XIC0(nSigmaTPCPrFromV0),
      FILL_MAP_XIC0(nSigmaTPCPiFromXic0),
      FILL_MAP_XIC0(nSigmaTPCPiFromCasc)};

    // KFParticle variables
    if constexpr (reconstructionType == aod::hf_cand_casc_lf::ConstructMethod::KfParticle) {
      std::map<std::string, uint8_t> mapKfFeatures{
        FILL_MAP_XIC0(decayLengthNormKf),
        FILL_MAP_XIC0(decayLengthNormKfCascade),
        FILL_MAP_XIC0(chi2GeoOverNdf),
        FILL_MAP_XIC0(chi2TopoOverNdfToPv),
        FILL_MAP_XIC0(chi2GeoOverNdfCascade),
        FILL_MAP_XIC0(chi2TopoOverNdfCascadeToCharmHadron),
        FILL_MAP_XIC0(chi2TopoOverNdfCascadeToPv),
        FILL_MAP_XIC0(chi2GeoOverNdfV0),
        FILL_MAP_XIC0(chi2TopoOverNdfV0ToPv),
        FILL_MAP_XIC0(chi2TopoOverNdfV0ToCascade)};

      MlResponse<TypeOutputScore>::mAvailableInputFeatures.insert(mapKfFeatures.begin(), mapKfFeatures.end());
    }
  }
};

} // namespace o2::analysis

#undef FILL_MAP_XIC0
#undef CHECK_AND_FILL_VEC_XIC0
#undef CHECK_AND_FILL_VEC_XIC0_FULL
#undef CHECK_AND_FILL_VEC_XIC0_HFHELPER

#endif // PWGHF_CORE_HFMLRESPONSEXIC0TOXIPIQA_H_
