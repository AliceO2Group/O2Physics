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

/// \file HfMlResponseB0ToJpsiK0StarReduced.h
/// \brief Class to compute the ML response for B0 → J/Psi K*0 analysis selections in the reduced format
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università degli Studi and INFN Torino

#ifndef PWGHF_CORE_HFMLRESPONSEB0TOJPSIK0STARREDUCED_H_
#define PWGHF_CORE_HFMLRESPONSEB0TOJPSIK0STARREDUCED_H_

#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"

#include "Tools/ML/MlResponse.h"

#include <CommonConstants/PhysicsConstants.h>

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_B0(FEATURE) \
  {                          \
    #FEATURE, static_cast<uint8_t>(InputFeaturesB0ToJpsiK0StarReduced::FEATURE)}

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_B0_FULL(OBJECT, FEATURE, GETTER)                 \
  case static_cast<uint8_t>(InputFeaturesB0ToJpsiK0StarReduced::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());                            \
    break;                                                                  \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the GETTER function taking OBJECT in argument
#define CHECK_AND_FILL_VEC_B0_FUNC(OBJECT, FEATURE, GETTER)                 \
  case static_cast<uint8_t>(InputFeaturesB0ToJpsiK0StarReduced::FEATURE): { \
    inputFeatures.emplace_back(GETTER(OBJECT));                             \
    break;                                                                  \
  }

// Specific case of CHECK_AND_FILL_VEC_B0_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_B0(GETTER)                                      \
  case static_cast<uint8_t>(InputFeaturesB0ToJpsiK0StarReduced::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());                        \
    break;                                                                 \
  }

// Specific case of CHECK_AND_FILL_VEC_B0_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER and mass hypotheses are needed
#define CHECK_AND_FILL_VEC_B0_WITH_MASS_HYPO(GETTER)                       \
  case static_cast<uint8_t>(InputFeaturesB0ToJpsiK0StarReduced::GETTER): { \
    std::array massesKPi{                                                  \
      o2::constants::physics::MassMuon,                                    \
      o2::constants::physics::MassMuon,                                    \
      o2::constants::physics::MassKPlus,                                   \
      o2::constants::physics::MassPiPlus};                                 \
    std::array massesPiK{                                                  \
      o2::constants::physics::MassMuon,                                    \
      o2::constants::physics::MassMuon,                                    \
      o2::constants::physics::MassPiPlus,                                  \
      o2::constants::physics::MassKPlus};                                  \
    if (caseK0StarToKPi) {                                                 \
      inputFeatures.emplace_back(candidate.GETTER(massesKPi));             \
    } else {                                                               \
      inputFeatures.emplace_back(candidate.GETTER(massesPiK));             \
    }                                                                      \
    break;                                                                 \
  }

namespace o2::analysis
{

enum class InputFeaturesB0ToJpsiK0StarReduced : uint8_t {
  ptProng0 = 0,
  ptProng1,
  impactParameter0,
  impactParameter1,
  impactParameter2,
  impactParameter3,
  impactParameterProduct,
  impactParameterProductJpsi,
  impactParameterProductK0Star,
  chi2PCA,
  decayLength,
  decayLengthXY,
  decayLengthNormalised,
  decayLengthXYNormalised,
  cpa,
  cpaXY,
  maxNormalisedDeltaIP,
  ctXY,
  tpcNSigmaKa0,
  tofNSigmaKa0,
  tpcTofNSigmaKa0,
  tpcNSigmaKa1,
  tofNSigmaKa1,
  tpcTofNSigmaKa1
};

template <typename TypeOutputScore = float>
class HfMlResponseB0ToJpsiK0StarReduced : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseB0ToJpsiK0StarReduced() = default;
  /// Default destructor
  virtual ~HfMlResponseB0ToJpsiK0StarReduced() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the B0 candidate
  /// \param prong1 is the candidate's prong1
  /// \param prong2 is the candidate's prong2
  /// \param caseK0StarToKPi whether we have K*0 → K+ Pi- or K*0bar → K- Pi+
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename T3>
  std::vector<float> getInputFeatures(T1 const& candidate,
                                      T2 const& prong1,
                                      T3 const& prong2,
                                      bool const caseK0StarToKPi)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_B0(ptProng0);
        CHECK_AND_FILL_VEC_B0(ptProng1);
        CHECK_AND_FILL_VEC_B0(impactParameter0);
        CHECK_AND_FILL_VEC_B0(impactParameter1);
        CHECK_AND_FILL_VEC_B0(impactParameter2);
        CHECK_AND_FILL_VEC_B0(impactParameter3);
        CHECK_AND_FILL_VEC_B0(impactParameterProduct);
        CHECK_AND_FILL_VEC_B0(impactParameterProductJpsi);
        CHECK_AND_FILL_VEC_B0(impactParameterProductK0Star);
        CHECK_AND_FILL_VEC_B0(chi2PCA);
        CHECK_AND_FILL_VEC_B0(decayLength);
        CHECK_AND_FILL_VEC_B0(decayLengthXY);
        CHECK_AND_FILL_VEC_B0(decayLengthNormalised);
        CHECK_AND_FILL_VEC_B0(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_B0(cpa);
        CHECK_AND_FILL_VEC_B0(cpaXY);
        CHECK_AND_FILL_VEC_B0(maxNormalisedDeltaIP);
        CHECK_AND_FILL_VEC_B0_WITH_MASS_HYPO(ctXY);
        // TPC PID variable
        CHECK_AND_FILL_VEC_B0_FULL(prong1, tpcNSigmaKa0, tpcNSigmaKa);
        // TOF PID variable
        CHECK_AND_FILL_VEC_B0_FULL(prong1, tofNSigmaKa0, tofNSigmaKa);
        // Combined PID variables
        CHECK_AND_FILL_VEC_B0_FUNC(prong1, tpcTofNSigmaKa0, o2::pid_tpc_tof_utils::getTpcTofNSigmaKa1);
        // TPC PID variable
        CHECK_AND_FILL_VEC_B0_FULL(prong2, tpcNSigmaKa1, tpcNSigmaKa);
        // TOF PID variable
        CHECK_AND_FILL_VEC_B0_FULL(prong2, tofNSigmaKa1, tofNSigmaKa);
        // Combined PID variables
        CHECK_AND_FILL_VEC_B0_FUNC(prong2, tpcTofNSigmaKa1, o2::pid_tpc_tof_utils::getTpcTofNSigmaKa1);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_B0(ptProng0),
      FILL_MAP_B0(ptProng1),
      FILL_MAP_B0(impactParameter0),
      FILL_MAP_B0(impactParameter1),
      FILL_MAP_B0(impactParameter2),
      FILL_MAP_B0(impactParameter3),
      FILL_MAP_B0(impactParameterProduct),
      FILL_MAP_B0(impactParameterProductJpsi),
      FILL_MAP_B0(impactParameterProductK0Star),
      FILL_MAP_B0(chi2PCA),
      FILL_MAP_B0(decayLength),
      FILL_MAP_B0(decayLengthXY),
      FILL_MAP_B0(decayLengthNormalised),
      FILL_MAP_B0(decayLengthXYNormalised),
      FILL_MAP_B0(cpa),
      FILL_MAP_B0(cpaXY),
      FILL_MAP_B0(maxNormalisedDeltaIP),
      FILL_MAP_B0(ctXY),
      // TPC PID variable
      FILL_MAP_B0(tpcNSigmaKa0),
      // TOF PID variable
      FILL_MAP_B0(tofNSigmaKa0),
      // Combined PID variable
      FILL_MAP_B0(tpcTofNSigmaKa0),
      // TPC PID variable
      FILL_MAP_B0(tpcNSigmaKa1),
      // TOF PID variable
      FILL_MAP_B0(tofNSigmaKa1),
      // Combined PID variable
      FILL_MAP_B0(tpcTofNSigmaKa1)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_B0
#undef CHECK_AND_FILL_VEC_B0_FULL
#undef CHECK_AND_FILL_VEC_B0_FUNC
#undef CHECK_AND_FILL_VEC_B0
#undef CHECK_AND_FILL_VEC_B0_WITH_MASS_HYPO

#endif // PWGHF_CORE_HFMLRESPONSEB0TOJPSIK0STARREDUCED_H_
