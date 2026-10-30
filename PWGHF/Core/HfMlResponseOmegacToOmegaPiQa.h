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
#include "PWGHF/DataModel/CandidateReconstructionTablesQa.h"

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

// Specific case of CHECK_AND_FILL_VEC_OMEGAC0_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_OMEGAC0(GETTER)                           \
  case static_cast<uint8_t>(InputFeaturesOmegacToOmegaPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());                  \
    break;                                                           \
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

// Variation of CHECK_AND_FILL_VEC_OMEGAC0_FULL(OBJECT, FEATURE, GETTER)
// where GETTER is a method of HfHelper
#define CHECK_AND_FILL_VEC_OMEGAC0_HFHELPER(OBJECT, FEATURE, GETTER)  \
  case static_cast<uint8_t>(InputFeaturesOmegacToOmegaPi::FEATURE): { \
    inputFeatures.emplace_back(HfHelperQa::GETTER(OBJECT));           \
    break;                                                            \
  }

namespace o2::analysis
{
enum class InputFeaturesOmegacToOmegaPi : uint8_t { // TODO: check if all the features needed are convered
  // Common variables
  impactParameterXYCharmBachelor = 0,
  impactParameterZCharmBachelor,
  impactParameterXYCascade,
  impactParameterZCascade,
  dcaDaughters,
  dcaCascadeDaughters,
  cpa,
  cpaCascade,
  cpaCascadeToCharmHadron,
  cosThetaStarPiFromOmegac,
  decayLengthXY,
  decayLengthXYNormalised,
  decayLengthXYCascade,
  massCascade,
  nSigmaTPCPiFromV0,
  nSigmaTPCPiFromOmegac,
  nSigmaTPCKaFromCasc,

  // KFParticle variables
  decayLengthNormKf,
  decayLengthNormKfCascade,
  chi2Geo,
  chi2TopoOverNdfToPv,
  chi2GeoCascade,
  chi2TopoOverNdfCascadeToCharmHadron,
  chi2TopoOverNdfCascadeToPv,

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
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, impactParameterXYCharmBachelor, impactParameter1);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, impactParameterZCharmBachelor, impactParameterZ1);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, impactParameterXYCascade, impactParameter0);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(candidate, impactParameterZCascade, impactParameterZ0);
        CHECK_AND_FILL_VEC_OMEGAC0(dcaDaughters);
        CHECK_AND_FILL_VEC_OMEGAC0(dcaCascadeDaughters);
        CHECK_AND_FILL_VEC_OMEGAC0(cpa);
        CHECK_AND_FILL_VEC_OMEGAC0(cpaCascade);
        CHECK_AND_FILL_VEC_OMEGAC0(cpaCascadeToCharmHadron);
        CHECK_AND_FILL_VEC_OMEGAC0_HFHELPER(candidate, cosThetaStarPiFromOmegac, cosThetaStarOmegaczeroToOmegaPi);
        CHECK_AND_FILL_VEC_OMEGAC0(decayLengthXY);
        CHECK_AND_FILL_VEC_OMEGAC0(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_OMEGAC0(decayLengthXYCascade);
        CHECK_AND_FILL_VEC_OMEGAC0(massCascade);
        // TPC PID variables
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(lamProngPi, nSigmaTPCPiFromV0, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(cascProng, nSigmaTPCKaFromCasc, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_OMEGAC0_FULL(charmBaryonProng, nSigmaTPCPiFromOmegac, tpcNSigmaPi);
      }

      // KFParticle variables
      if constexpr (reconstructionType == aod::hf_cand_casc_lf::ConstructMethod::KfParticle) {
        switch (idx) {
          CHECK_AND_FILL_VEC_OMEGAC0(decayLengthNormKf);
          CHECK_AND_FILL_VEC_OMEGAC0(decayLengthNormKfCascade);
          CHECK_AND_FILL_VEC_OMEGAC0(chi2Geo);
          CHECK_AND_FILL_VEC_OMEGAC0(chi2TopoOverNdfCascadeToPv);
          CHECK_AND_FILL_VEC_OMEGAC0(chi2GeoCascade);
          CHECK_AND_FILL_VEC_OMEGAC0(chi2TopoOverNdfCascadeToCharmHadron);
          CHECK_AND_FILL_VEC_OMEGAC0(chi2TopoOverNdfToPv);
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
      FILL_MAP_OMEGAC0(impactParameterXYCharmBachelor),
      FILL_MAP_OMEGAC0(impactParameterZCharmBachelor),
      FILL_MAP_OMEGAC0(impactParameterXYCascade),
      FILL_MAP_OMEGAC0(impactParameterZCascade),
      FILL_MAP_OMEGAC0(dcaDaughters),
      FILL_MAP_OMEGAC0(dcaCascadeDaughters),
      FILL_MAP_OMEGAC0(cpa),
      FILL_MAP_OMEGAC0(cpaCascade),
      FILL_MAP_OMEGAC0(cpaCascadeToCharmHadron),
      FILL_MAP_OMEGAC0(cosThetaStarPiFromOmegac),
      FILL_MAP_OMEGAC0(decayLengthXY),
      FILL_MAP_OMEGAC0(decayLengthXYNormalised),
      FILL_MAP_OMEGAC0(decayLengthXYCascade),
      FILL_MAP_OMEGAC0(massCascade),
      FILL_MAP_OMEGAC0(nSigmaTPCPiFromV0),
      FILL_MAP_OMEGAC0(nSigmaTPCPiFromOmegac),
      FILL_MAP_OMEGAC0(nSigmaTPCKaFromCasc)};

    // KFParticle variables
    if constexpr (reconstructionType == aod::hf_cand_casc_lf::ConstructMethod::KfParticle) {

      std::map<std::string, uint8_t> mapKfFeatures{
        FILL_MAP_OMEGAC0(decayLengthNormKf),
        FILL_MAP_OMEGAC0(decayLengthNormKfCascade),
        FILL_MAP_OMEGAC0(chi2Geo),
        FILL_MAP_OMEGAC0(chi2TopoOverNdfToPv),
        FILL_MAP_OMEGAC0(chi2GeoCascade),
        FILL_MAP_OMEGAC0(chi2TopoOverNdfCascadeToCharmHadron),
        FILL_MAP_OMEGAC0(chi2TopoOverNdfCascadeToPv)};

      MlResponse<TypeOutputScore>::mAvailableInputFeatures.insert(mapKfFeatures.begin(), mapKfFeatures.end());
    }
  }
};

} // namespace o2::analysis

#undef FILL_MAP_OMEGAC0
#undef CHECK_AND_FILL_VEC_OMEGAC0
#undef CHECK_AND_FILL_VEC_OMEGAC0_FULL
#undef CHECK_AND_FILL_VEC_OMEGAC0_HFHELPER
#endif // PWGHF_CORE_HFMLRESPONSEOMEGACTOOMEGAPIQA_H_
