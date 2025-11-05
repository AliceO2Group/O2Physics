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

/// \file HfMlResponseLcToPKPi.h
/// \brief Class to compute the ML response for Lc+ → p K- π+ analysis selections
/// \author Grazia Luparello <grazia.luparello@cern.ch>

#ifndef PWGHF_CORE_HFMLRESPONSELCTOPKPI_H_
#define PWGHF_CORE_HFMLRESPONSELCTOPKPI_H_

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
#define FILL_MAP_LCTOPKPI(FEATURE)                                 \
  {                                                                \
    #FEATURE, static_cast<uint8_t>(InputFeaturesLcToPKPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_LCTOPKPI_FULL(OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesLcToPKPi::FEATURE): {    \
    inputFeatures.emplace_back(OBJECT.GETTER());                  \
    break;                                                        \
  }

// Specific case of CHECK_AND_FILL_VEC_LCTOPKPI_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_LCTOPKPI(GETTER)                   \
  case static_cast<uint8_t>(InputFeaturesLcToPKPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());           \
    break;                                                    \
  }

// Variation of CHECK_AND_FILL_VEC_LCTOPKPI_FULL(OBJECT, FEATURE, GETTER)
// where GETTER is a method of HfHelper
#define CHECK_AND_FILL_VEC_LCTOPKPI_HFHELPER(OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesLcToPKPi::FEATURE): {        \
    inputFeatures.emplace_back(HfHelper::GETTER(OBJECT));             \
    break;                                                            \
  }

// Variation of CHECK_AND_FILL_VEC_LCTOPKPI_OBJECT_SIGNED(OBJECT1, OBJECT2, FEATURE, GETTER)
// where OBJECT1 and OBJECT2 are the objects from which we call the GETTER method, and the variable
// is filled depending on whether it is a LcToPKPi or a LcToPiKP
#define CHECK_AND_FILL_VEC_LCTOPKPI_OBJECT_SIGNED(OBJECT1, OBJECT2, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesLcToPKPi::FEATURE): {                       \
    if (caseLcToPKPi) {                                                              \
      inputFeatures.emplace_back(OBJECT1.GETTER());                                  \
    } else {                                                                         \
      inputFeatures.emplace_back(OBJECT2.GETTER());                                  \
    }                                                                                \
    break;                                                                           \
  }

// Variation of CHECK_AND_FILL_VEC_LCTOPKPI_HFHELPER_SIGNED(OBJECT, FEATURE, GETTER1, GETTER2)
// where GETTER1 and GETTER2 are methods of the OBJECT
#define CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(OBJECT, FEATURE, GETTER1, GETTER2) \
  case static_cast<uint8_t>(InputFeaturesLcToPKPi::FEATURE): {                \
    if (caseLcToPKPi) {                                                       \
      inputFeatures.emplace_back(OBJECT.GETTER1());                           \
    } else {                                                                  \
      inputFeatures.emplace_back(OBJECT.GETTER2());                           \
    }                                                                         \
    break;                                                                    \
  }

namespace o2::analysis
{
enum class InputFeaturesLcToPKPi : uint8_t {
  ptProng0 = 0,
  ptProng1,
  ptProng2,
  impactParameterXY0,
  impactParameterXY1,
  impactParameterXY2,
  impactParameterZ0,
  impactParameterZ1,
  impactParameterZ2,
  decayLength,
  decayLengthXY,
  decayLengthXYNormalised,
  cpa,
  cpaXY,
  chi2PCA,
  tpcNSigmaPr0, // 0
  tpcNSigmaKa0, // 0
  tpcNSigmaPi0, // 0
  tpcNSigmaPr1, // 1
  tpcNSigmaKa1, // 1
  tpcNSigmaPi1, // 1
  tpcNSigmaPr2, // 2
  tpcNSigmaKa2, // 2
  tpcNSigmaPi2, // 2
  tofNSigmaPr0, //
  tofNSigmaKa0, //
  tofNSigmaPi0, //
  tofNSigmaPr1,
  tofNSigmaKa1,
  tofNSigmaPi1,
  tofNSigmaPr2,
  tofNSigmaKa2,
  tofNSigmaPi2,
  tpcTofNSigmaPi0,
  tpcTofNSigmaPi1,
  tpcTofNSigmaPi2,
  tpcTofNSigmaKa0,
  tpcTofNSigmaKa1,
  tpcTofNSigmaKa2,
  tpcTofNSigmaPr0,
  tpcTofNSigmaPr1,
  tpcTofNSigmaPr2,
  tpcNSigmaPrExpPr0,
  tpcNSigmaPiExpPi2,
  tofNSigmaPrExpPr0,
  tofNSigmaPiExpPi2,
  tpcTofNSigmaPrExpPr0,
  tpcTofNSigmaPiExpPi2,
  kfChi2PrimProton,
  kfChi2PrimKaon,
  kfChi2PrimPion,
  kfChi2GeoKaonPion,
  kfChi2GeoProtonPion,
  kfChi2GeoProtonKaon,
  kfDcaKaonPion,
  kfDcaProtonPion,
  kfDcaProtonKaon,
  kfChi2Geo,
  kfChi2Topo,
  kfDecayLengthNormalised
};

template <typename TypeOutputScore = float, aod::hf_cand::VertexerType reconstructionType = aod::hf_cand::VertexerType::DCAFitter>
class HfMlResponseLcToPKPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseLcToPKPi() = default;
  /// Default destructor
  virtual ~HfMlResponseLcToPKPi() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Lc candidate
  /// \param prong0 is the candidate's prong0
  /// \param prong1 is the candidate's prong1
  /// \param prong2 is the candidate's prong2
  /// \return inputFeatures vector
  template <typename T1>
  std::vector<float> getInputFeatures(T1 const& candidate, bool const caseLcToPKPi)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_LCTOPKPI(ptProng0);
        CHECK_AND_FILL_VEC_LCTOPKPI(ptProng1);
        CHECK_AND_FILL_VEC_LCTOPKPI(ptProng2);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, impactParameterXY0, impactParameter0);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, impactParameterXY1, impactParameter1);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, impactParameterXY2, impactParameter2);
        CHECK_AND_FILL_VEC_LCTOPKPI(impactParameterZ0);
        CHECK_AND_FILL_VEC_LCTOPKPI(impactParameterZ1);
        CHECK_AND_FILL_VEC_LCTOPKPI(impactParameterZ2);
        CHECK_AND_FILL_VEC_LCTOPKPI(decayLength);
        CHECK_AND_FILL_VEC_LCTOPKPI(decayLengthXY);
        CHECK_AND_FILL_VEC_LCTOPKPI(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_LCTOPKPI(cpa);
        CHECK_AND_FILL_VEC_LCTOPKPI(cpaXY);
        CHECK_AND_FILL_VEC_LCTOPKPI(chi2PCA);
        // TPC PID variables
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcNSigmaPr0, nSigTpcPr0);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcNSigmaKa0, nSigTpcKa0);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcNSigmaPi0, nSigTpcPi0);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcNSigmaPr1, nSigTpcPr1);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcNSigmaKa1, nSigTpcKa1);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcNSigmaPi1, nSigTpcPi1);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcNSigmaPr2, nSigTpcPr2);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcNSigmaKa2, nSigTpcKa2);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcNSigmaPi2, nSigTpcPi2);
        CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(candidate, tpcNSigmaPrExpPr0, nSigTpcPr0, nSigTpcPr2);
        CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(candidate, tpcNSigmaPiExpPi2, nSigTpcPi2, nSigTpcPi0);
        // TOF PID variables
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tofNSigmaPr0, nSigTofPr0);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tofNSigmaKa0, nSigTofKa0);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tofNSigmaPi0, nSigTofPi0);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tofNSigmaPr1, nSigTofPr1);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tofNSigmaKa1, nSigTofKa1);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tofNSigmaPi1, nSigTofPi1);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tofNSigmaPr2, nSigTofPr2);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tofNSigmaKa2, nSigTofKa2);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tofNSigmaPi2, nSigTofPi2);
        CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(candidate, tofNSigmaPrExpPr0, nSigTofPr0, nSigTofPr2);
        CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(candidate, tofNSigmaPiExpPi2, nSigTofPi2, nSigTofPi0);
        // Combined PID variables
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcTofNSigmaPi0, tpcTofNSigmaPi0);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcTofNSigmaPi1, tpcTofNSigmaPi1);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcTofNSigmaPi2, tpcTofNSigmaPi2);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcTofNSigmaKa0, tpcTofNSigmaKa0);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcTofNSigmaKa1, tpcTofNSigmaKa1);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcTofNSigmaKa2, tpcTofNSigmaKa2);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcTofNSigmaPr0, tpcTofNSigmaPr0);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcTofNSigmaPr1, tpcTofNSigmaPr1);
        CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, tpcTofNSigmaPr2, tpcTofNSigmaPr2);
        CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(candidate, tpcTofNSigmaPrExpPr0, tpcTofNSigmaPr0, tpcTofNSigmaPr2);
        CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(candidate, tpcTofNSigmaPiExpPi2, tpcTofNSigmaPi2, tpcTofNSigmaPi0);
      }
      if constexpr (reconstructionType == aod::hf_cand::VertexerType::KfParticle) {
        switch (idx) {
          CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(candidate, kfChi2PrimProton, kfChi2PrimProng0, kfChi2PrimProng2);
          CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, kfChi2PrimKaon, kfChi2PrimProng1);
          CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(candidate, kfChi2PrimPion, kfChi2PrimProng2, kfChi2PrimProng0);
          CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(candidate, kfChi2GeoKaonPion, kfChi2GeoProng1Prong2, kfChi2GeoProng0Prong1);
          CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, kfChi2GeoProtonPion, kfChi2GeoProng0Prong2);
          CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(candidate, kfChi2GeoProtonKaon, kfChi2GeoProng0Prong1, kfChi2GeoProng1Prong2);
          CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(candidate, kfDcaKaonPion, kfDcaProng1Prong2, kfDcaProng0Prong1);
          CHECK_AND_FILL_VEC_LCTOPKPI_FULL(candidate, kfDcaProtonPion, kfDcaProng0Prong2);
          CHECK_AND_FILL_VEC_LCTOPKPI_SIGNED(candidate, kfDcaProtonKaon, kfDcaProng0Prong1, kfDcaProng1Prong2);
          CHECK_AND_FILL_VEC_LCTOPKPI(kfChi2Geo);
          CHECK_AND_FILL_VEC_LCTOPKPI(kfChi2Topo);
          case static_cast<uint8_t>(InputFeaturesLcToPKPi::kfDecayLengthNormalised): {
            inputFeatures.emplace_back(candidate.kfDecayLength() / candidate.kfDecayLengthError());
            break;
          }
        }
      }
    }
    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_LCTOPKPI(ptProng0),
      FILL_MAP_LCTOPKPI(ptProng1),
      FILL_MAP_LCTOPKPI(ptProng2),
      FILL_MAP_LCTOPKPI(impactParameterXY0),
      FILL_MAP_LCTOPKPI(impactParameterXY1),
      FILL_MAP_LCTOPKPI(impactParameterXY2),
      FILL_MAP_LCTOPKPI(impactParameterZ0),
      FILL_MAP_LCTOPKPI(impactParameterZ1),
      FILL_MAP_LCTOPKPI(impactParameterZ2),
      FILL_MAP_LCTOPKPI(decayLength),
      FILL_MAP_LCTOPKPI(decayLengthXY),
      FILL_MAP_LCTOPKPI(decayLengthXYNormalised),
      FILL_MAP_LCTOPKPI(cpa),
      FILL_MAP_LCTOPKPI(cpaXY),
      FILL_MAP_LCTOPKPI(chi2PCA),
      // TPC PID variables
      FILL_MAP_LCTOPKPI(tpcNSigmaPr0),
      FILL_MAP_LCTOPKPI(tpcNSigmaKa0),
      FILL_MAP_LCTOPKPI(tpcNSigmaPi0),
      FILL_MAP_LCTOPKPI(tpcNSigmaPr1),
      FILL_MAP_LCTOPKPI(tpcNSigmaKa1),
      FILL_MAP_LCTOPKPI(tpcNSigmaPi1),
      FILL_MAP_LCTOPKPI(tpcNSigmaPr2),
      FILL_MAP_LCTOPKPI(tpcNSigmaKa2),
      FILL_MAP_LCTOPKPI(tpcNSigmaPi2),
      FILL_MAP_LCTOPKPI(tpcNSigmaPrExpPr0),
      FILL_MAP_LCTOPKPI(tpcNSigmaPiExpPi2),
      // TOF PID variables
      FILL_MAP_LCTOPKPI(tofNSigmaPr0),
      FILL_MAP_LCTOPKPI(tofNSigmaKa0),
      FILL_MAP_LCTOPKPI(tofNSigmaPi0),
      FILL_MAP_LCTOPKPI(tofNSigmaPr1),
      FILL_MAP_LCTOPKPI(tofNSigmaKa1),
      FILL_MAP_LCTOPKPI(tofNSigmaPi1),
      FILL_MAP_LCTOPKPI(tofNSigmaPr2),
      FILL_MAP_LCTOPKPI(tofNSigmaKa2),
      FILL_MAP_LCTOPKPI(tofNSigmaPi2),
      FILL_MAP_LCTOPKPI(tofNSigmaPrExpPr0),
      FILL_MAP_LCTOPKPI(tofNSigmaPiExpPi2),
      // Combined PID variables
      FILL_MAP_LCTOPKPI(tpcTofNSigmaPi0),
      FILL_MAP_LCTOPKPI(tpcTofNSigmaPi1),
      FILL_MAP_LCTOPKPI(tpcTofNSigmaPi2),
      FILL_MAP_LCTOPKPI(tpcTofNSigmaKa0),
      FILL_MAP_LCTOPKPI(tpcTofNSigmaKa1),
      FILL_MAP_LCTOPKPI(tpcTofNSigmaKa2),
      FILL_MAP_LCTOPKPI(tpcTofNSigmaPr0),
      FILL_MAP_LCTOPKPI(tpcTofNSigmaPr1),
      FILL_MAP_LCTOPKPI(tpcTofNSigmaPr2),
      FILL_MAP_LCTOPKPI(tpcTofNSigmaPrExpPr0),
      FILL_MAP_LCTOPKPI(tpcTofNSigmaPiExpPi2)};
    if constexpr (reconstructionType == aod::hf_cand::VertexerType::KfParticle) {
      std::map<std::string, uint8_t> mapKfFeatures{
        // KFParticle variables
        FILL_MAP_LCTOPKPI(kfChi2PrimProton),
        FILL_MAP_LCTOPKPI(kfChi2PrimKaon),
        FILL_MAP_LCTOPKPI(kfChi2PrimPion),
        FILL_MAP_LCTOPKPI(kfChi2GeoKaonPion),
        FILL_MAP_LCTOPKPI(kfChi2GeoProtonPion),
        FILL_MAP_LCTOPKPI(kfChi2GeoProtonKaon),
        FILL_MAP_LCTOPKPI(kfDcaKaonPion),
        FILL_MAP_LCTOPKPI(kfDcaProtonPion),
        FILL_MAP_LCTOPKPI(kfDcaProtonKaon),
        FILL_MAP_LCTOPKPI(kfChi2Geo),
        FILL_MAP_LCTOPKPI(kfChi2Topo),
        FILL_MAP_LCTOPKPI(kfDecayLengthNormalised)};
      MlResponse<TypeOutputScore>::mAvailableInputFeatures.insert(mapKfFeatures.begin(), mapKfFeatures.end());
    }
  }
};

} // namespace o2::analysis

#undef FILL_MAP_LCTOPKPI
#undef CHECK_AND_FILL_VEC_LCTOPKPI_FULL
#undef CHECK_AND_FILL_VEC_LCTOPKPI
#undef CHECK_AND_FILL_VEC_LCTOPKPI_HFHELPER
#undef CHECK_AND_FILL_VEC_LCTOPKPI_OBJECT_SIGNED

#endif // PWGHF_CORE_HFMLRESPONSELCTOPKPI_H_
