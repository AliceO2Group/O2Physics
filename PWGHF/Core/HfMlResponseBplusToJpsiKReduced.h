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

/// \file HfMlResponseBplusToJpsiKReduced.h
/// \brief Class to compute the ML response for B± → J/Psi K± analysis selections in the reduced format
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università degli Studi and INFN Torino

#ifndef PWGHF_CORE_HFMLRESPONSEBPLUSTOJPSIKREDUCED_H_
#define PWGHF_CORE_HFMLRESPONSEBPLUSTOJPSIKREDUCED_H_

#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"

#include "Tools/ML/MlResponse.h"

#include <CommonConstants/PhysicsConstants.h>

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_BPLUS(FEATURE)                                               \
  {                                                                           \
    #FEATURE, static_cast<uint8_t>(InputFeaturesBplusToJpsiKReduced::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_BPLUS_FULL(OBJECT, FEATURE, GETTER)            \
  case static_cast<uint8_t>(InputFeaturesBplusToJpsiKReduced::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());                          \
    break;                                                                \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the GETTER function taking OBJECT in argument
#define CHECK_AND_FILL_VEC_BPLUS_FUNC(OBJECT, FEATURE, GETTER)            \
  case static_cast<uint8_t>(InputFeaturesBplusToJpsiKReduced::FEATURE): { \
    inputFeatures.emplace_back(GETTER(OBJECT));                           \
    break;                                                                \
  }

// Specific case of CHECK_AND_FILL_VEC_BPLUS_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_BPLUS(GETTER)                                 \
  case static_cast<uint8_t>(InputFeaturesBplusToJpsiKReduced::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());                      \
    break;                                                               \
  }

// Specific case of CHECK_AND_FILL_VEC_BPLUS_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate, FEATURE = GETTER, and args are needed
#define CHECK_AND_FILL_VEC_BPLUS_WITH_ARGS(GETTER, ARGS...)              \
  case static_cast<uint8_t>(InputFeaturesBplusToJpsiKReduced::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER(ARGS));                  \
    break;                                                               \
  }

namespace o2::analysis
{

enum class InputFeaturesBplusToJpsiKReduced : uint8_t {
  ptProng0 = 0,
  ptProng1,
  impactParameter0,
  impactParameter1,
  impactParameter2,
  impactParameterProduct,
  impactParameterProductJpsi,
  chi2PCA,
  decayLength,
  decayLengthXY,
  decayLengthNormalised,
  decayLengthXYNormalised,
  cpa,
  cpaXY,
  maxNormalisedDeltaIP,
  ctXY,
  tpcNSigmaKa1,
  tofNSigmaKa1,
  tpcTofNSigmaKa1
};

template <typename TypeOutputScore = float>
class HfMlResponseBplusToJpsiKReduced : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseBplusToJpsiKReduced() = default;
  /// Default destructor
  virtual ~HfMlResponseBplusToJpsiKReduced() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the B+ candidate
  /// \param prong1 is the candidate's prong1
  /// \return inputFeatures vector
  template <typename T1, typename T2>
  std::vector<float> getInputFeatures(T1 const& candidate,
                                      T2 const& prong1)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_BPLUS(ptProng0);
        CHECK_AND_FILL_VEC_BPLUS(ptProng1);
        CHECK_AND_FILL_VEC_BPLUS(impactParameter0);
        CHECK_AND_FILL_VEC_BPLUS(impactParameter1);
        CHECK_AND_FILL_VEC_BPLUS(impactParameter2);
        CHECK_AND_FILL_VEC_BPLUS(impactParameterProduct);
        CHECK_AND_FILL_VEC_BPLUS(impactParameterProductJpsi);
        CHECK_AND_FILL_VEC_BPLUS(chi2PCA);
        CHECK_AND_FILL_VEC_BPLUS(decayLength);
        CHECK_AND_FILL_VEC_BPLUS(decayLengthXY);
        CHECK_AND_FILL_VEC_BPLUS(decayLengthNormalised);
        CHECK_AND_FILL_VEC_BPLUS(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_BPLUS(cpa);
        CHECK_AND_FILL_VEC_BPLUS(cpaXY);
        CHECK_AND_FILL_VEC_BPLUS(maxNormalisedDeltaIP);
        CHECK_AND_FILL_VEC_BPLUS_WITH_ARGS(ctXY, std::array{o2::constants::physics::MassMuon, o2::constants::physics::MassMuon, o2::constants::physics::MassKPlus});
        // TPC PID variable
        CHECK_AND_FILL_VEC_BPLUS_FULL(prong1, tpcNSigmaKa1, tpcNSigmaKa);
        // TOF PID variable
        CHECK_AND_FILL_VEC_BPLUS_FULL(prong1, tofNSigmaKa1, tofNSigmaKa);
        // Combined PID variables
        CHECK_AND_FILL_VEC_BPLUS_FUNC(prong1, tpcTofNSigmaKa1, o2::pid_tpc_tof_utils::getTpcTofNSigmaKa1);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_BPLUS(ptProng0),
      FILL_MAP_BPLUS(ptProng1),
      FILL_MAP_BPLUS(impactParameter0),
      FILL_MAP_BPLUS(impactParameter1),
      FILL_MAP_BPLUS(impactParameter2),
      FILL_MAP_BPLUS(impactParameterProduct),
      FILL_MAP_BPLUS(impactParameterProductJpsi),
      FILL_MAP_BPLUS(chi2PCA),
      FILL_MAP_BPLUS(decayLength),
      FILL_MAP_BPLUS(decayLengthXY),
      FILL_MAP_BPLUS(decayLengthNormalised),
      FILL_MAP_BPLUS(decayLengthXYNormalised),
      FILL_MAP_BPLUS(cpa),
      FILL_MAP_BPLUS(cpaXY),
      FILL_MAP_BPLUS(maxNormalisedDeltaIP),
      FILL_MAP_BPLUS(ctXY),
      // TPC PID variable
      FILL_MAP_BPLUS(tpcNSigmaKa1),
      // TOF PID variable
      FILL_MAP_BPLUS(tofNSigmaKa1),
      // Combined PID variable
      FILL_MAP_BPLUS(tpcTofNSigmaKa1)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_BPLUS
#undef CHECK_AND_FILL_VEC_BPLUS_FULL
#undef CHECK_AND_FILL_VEC_BPLUS_FUNC
#undef CHECK_AND_FILL_VEC_BPLUS

#endif // PWGHF_CORE_HFMLRESPONSEBPLUSTOJPSIKREDUCED_H_
