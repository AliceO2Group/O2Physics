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

/// \file HfMlResponseDplusToPiKPi.h
/// \brief Class to compute the ML response for D± → π± K∓ π± analysis selections
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#ifndef PWGHF_CORE_HFMLRESPONSEDPLUSTOPIKPI_H_
#define PWGHF_CORE_HFMLRESPONSEDPLUSTOPIKPI_H_

#include <map>
#include <string>
#include <vector>

#include "PWGHF/Core/HfMlResponse.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_DPLUS(FEATURE)                                          \
  {                                                                      \
#FEATURE, static_cast < uint8_t>(InputFeaturesDplusToPiKPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_DPLUS_FULL(OBJECT, FEATURE, GETTER)     \
  case static_cast<uint8_t>(InputFeaturesDplusToPiKPi::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());                   \
    break;                                                         \
  }

// Specific case of CHECK_AND_FILL_VEC_DPLUS_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_DPLUS(GETTER)                          \
  case static_cast<uint8_t>(InputFeaturesDplusToPiKPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());               \
    break;                                                        \
  }

namespace o2::analysis
{

enum class InputFeaturesDplusToPiKPi : uint8_t {
  ptProng0 = 0,
  ptProng1,
  ptProng2,
  impactParameterXY0,
  impactParameterXY1,
  impactParameterXY2,
  decayLength,
  decayLengthXYNormalised,
  cpa,
  cpaXY,
  maxNormalisedDeltaIP,
  tpcNSigmaPi0,
  tpcNSigmaKa0,
  tpcNSigmaPi1,
  tpcNSigmaKa1,
  tpcNSigmaPi2,
  tpcNSigmaKa2,
  tofNSigmaPi0,
  tofNSigmaKa0,
  tofNSigmaPi1,
  tofNSigmaKa1,
  tofNSigmaPi2,
  tofNSigmaKa2,
  tpcTofNSigmaPi0,
  tpcTofNSigmaPi1,
  tpcTofNSigmaPi2,
  tpcTofNSigmaKa0,
  tpcTofNSigmaKa1,
  tpcTofNSigmaKa2
};

template <typename TypeOutputScore = float>
class HfMlResponseDplusToPiKPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseDplusToPiKPi() = default;
  /// Default destructor
  virtual ~HfMlResponseDplusToPiKPi() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Dplus candidate
  /// \param prong0 is the candidate's prong0
  /// \param prong1 is the candidate's prong1
  /// \param prong2 is the candidate's prong2
  /// \return inputFeatures vector
  template <typename T1, typename T2>
  std::vector<float> getInputFeatures(T1 const& candidate,
                                      T2 const& prong0, T2 const& prong1, T2 const& prong2)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_DPLUS(ptProng0);
        CHECK_AND_FILL_VEC_DPLUS(ptProng1);
        CHECK_AND_FILL_VEC_DPLUS(ptProng2);
        CHECK_AND_FILL_VEC_DPLUS_FULL(candidate, impactParameterXY0, impactParameter0);
        CHECK_AND_FILL_VEC_DPLUS_FULL(candidate, impactParameterXY1, impactParameter1);
        CHECK_AND_FILL_VEC_DPLUS_FULL(candidate, impactParameterXY2, impactParameter2);
        CHECK_AND_FILL_VEC_DPLUS(decayLength);
        CHECK_AND_FILL_VEC_DPLUS(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_DPLUS(cpa);
        CHECK_AND_FILL_VEC_DPLUS(cpaXY);
        CHECK_AND_FILL_VEC_DPLUS(maxNormalisedDeltaIP);
        // TPC PID variables
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong0, tpcNSigmaPi0, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong0, tpcNSigmaKa0, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong1, tpcNSigmaKa1, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong2, tpcNSigmaPi2, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong2, tpcNSigmaKa2, tpcNSigmaKa);
        // TOF PID variables
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong0, tofNSigmaPi0, tofNSigmaPi);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong0, tofNSigmaKa0, tofNSigmaKa);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong1, tofNSigmaKa1, tofNSigmaKa);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong2, tofNSigmaPi2, tofNSigmaPi);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong2, tofNSigmaKa2, tofNSigmaKa);
        // Combined PID variables
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong0, tpcTofNSigmaPi0, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong1, tpcTofNSigmaPi1, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong2, tpcTofNSigmaPi2, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong0, tpcTofNSigmaKa0, tpcTofNSigmaKa);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong1, tpcTofNSigmaKa1, tpcTofNSigmaKa);
        CHECK_AND_FILL_VEC_DPLUS_FULL(prong2, tpcTofNSigmaKa2, tpcTofNSigmaKa);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_DPLUS(ptProng0),
      FILL_MAP_DPLUS(ptProng1),
      FILL_MAP_DPLUS(ptProng2),
      FILL_MAP_DPLUS(impactParameterXY0),
      FILL_MAP_DPLUS(impactParameterXY1),
      FILL_MAP_DPLUS(impactParameterXY2),
      FILL_MAP_DPLUS(decayLength),
      FILL_MAP_DPLUS(decayLengthXYNormalised),
      FILL_MAP_DPLUS(cpa),
      FILL_MAP_DPLUS(cpaXY),
      FILL_MAP_DPLUS(maxNormalisedDeltaIP),
      // TPC PID variables
      FILL_MAP_DPLUS(tpcNSigmaPi0),
      FILL_MAP_DPLUS(tpcNSigmaKa0),
      FILL_MAP_DPLUS(tpcNSigmaPi1),
      FILL_MAP_DPLUS(tpcNSigmaKa1),
      FILL_MAP_DPLUS(tpcNSigmaPi2),
      FILL_MAP_DPLUS(tpcNSigmaKa2),
      // TOF PID variables
      FILL_MAP_DPLUS(tofNSigmaPi0),
      FILL_MAP_DPLUS(tofNSigmaKa0),
      FILL_MAP_DPLUS(tofNSigmaPi1),
      FILL_MAP_DPLUS(tofNSigmaKa1),
      FILL_MAP_DPLUS(tofNSigmaPi2),
      FILL_MAP_DPLUS(tofNSigmaKa2),
      // Combined PID variables
      FILL_MAP_DPLUS(tpcTofNSigmaPi0),
      FILL_MAP_DPLUS(tpcTofNSigmaPi1),
      FILL_MAP_DPLUS(tpcTofNSigmaPi2),
      FILL_MAP_DPLUS(tpcTofNSigmaKa0),
      FILL_MAP_DPLUS(tpcTofNSigmaKa1),
      FILL_MAP_DPLUS(tpcTofNSigmaKa2)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_DPLUS
#undef CHECK_AND_FILL_VEC_DPLUS_FULL
#undef CHECK_AND_FILL_VEC_DPLUS

#endif // PWGHF_CORE_HFMLRESPONSEDPLUSTOPIKPI_H_
