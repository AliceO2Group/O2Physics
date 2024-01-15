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

#include <vector>

#include "PWGHF/Core/HfMlResponse.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_LC(FEATURE)                                         \
  {                                                                  \
#FEATURE, static_cast < uint8_t>(InputFeaturesLcToPKPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_LC_FULL(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesLcToPKPi::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());               \
    break;                                                     \
  }

// Specific case of CHECK_AND_FILL_VEC_LC_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_LC(GETTER)                         \
  case static_cast<uint8_t>(InputFeaturesLcToPKPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());           \
    break;                                                    \
  }

// Variation of CHECK_AND_FILL_VEC_LC_FULL(OBJECT, FEATURE, GETTER)
// where GETTER is a method of hfHelper
#define CHECK_AND_FILL_VEC_LC_HFHELPER(OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesLcToPKPi::FEATURE): {  \
    inputFeatures.emplace_back(hfHelper.GETTER(OBJECT));        \
    break;                                                      \
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
  decayLength,
  decayLengthXY,
  decayLengthXYNormalised,
  cpa,
  cpaXY,
  tpcNSigmaP0,  // 0
  tpcNSigmaKa0, // 0
  tpcNSigmaPi0, // 0
  tpcNSigmaP1,  // 1
  tpcNSigmaKa1, // 1
  tpcNSigmaPi1, // 1
  tpcNSigmaP2,  // 2
  tpcNSigmaKa2, // 2
  tpcNSigmaPi2, // 2
  tofNSigmaP0,  //
  tofNSigmaKa0, //
  tofNSigmaPi0, //
  tofNSigmaP1,
  tofNSigmaKa1,
  tofNSigmaPi1,
  tofNSigmaP2,
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
  tpcTofNSigmaPr2
};

template <typename TypeOutputScore = float>
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
  template <typename T1, typename T2>
  std::vector<float> getInputFeatures(T1 const& candidate,
                                      T2 const& prong0, T2 const& prong1, T2 const& prong2)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_LC(ptProng0);
        CHECK_AND_FILL_VEC_LC(ptProng1);
        CHECK_AND_FILL_VEC_LC(ptProng2);
        CHECK_AND_FILL_VEC_LC_FULL(candidate, impactParameterXY0, impactParameter0);
        CHECK_AND_FILL_VEC_LC_FULL(candidate, impactParameterXY1, impactParameter1);
        CHECK_AND_FILL_VEC_LC_FULL(candidate, impactParameterXY2, impactParameter2);
        CHECK_AND_FILL_VEC_LC(decayLength);
        CHECK_AND_FILL_VEC_LC(decayLengthXY);
        CHECK_AND_FILL_VEC_LC(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_LC(cpa);
        CHECK_AND_FILL_VEC_LC(cpaXY);
        // TPC PID variables
        CHECK_AND_FILL_VEC_LC_FULL(prong0, tpcNSigmaP0, tpcNSigmaPr);
        CHECK_AND_FILL_VEC_LC_FULL(prong0, tpcNSigmaKa0, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_LC_FULL(prong0, tpcNSigmaPi0, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_LC_FULL(prong1, tpcNSigmaP1, tpcNSigmaPr);
        CHECK_AND_FILL_VEC_LC_FULL(prong1, tpcNSigmaKa1, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_LC_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_LC_FULL(prong2, tpcNSigmaP2, tpcNSigmaPr);
        CHECK_AND_FILL_VEC_LC_FULL(prong2, tpcNSigmaKa2, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_LC_FULL(prong2, tpcNSigmaPi2, tpcNSigmaPi);
        // TOF PID variables
        CHECK_AND_FILL_VEC_LC_FULL(prong0, tofNSigmaP0, tofNSigmaPr);
        CHECK_AND_FILL_VEC_LC_FULL(prong0, tofNSigmaKa0, tofNSigmaKa);
        CHECK_AND_FILL_VEC_LC_FULL(prong0, tofNSigmaPi0, tofNSigmaPi);
        CHECK_AND_FILL_VEC_LC_FULL(prong1, tofNSigmaP1, tofNSigmaPr);
        CHECK_AND_FILL_VEC_LC_FULL(prong1, tofNSigmaKa1, tofNSigmaKa);
        CHECK_AND_FILL_VEC_LC_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
        CHECK_AND_FILL_VEC_LC_FULL(prong2, tofNSigmaP2, tofNSigmaPr);
        CHECK_AND_FILL_VEC_LC_FULL(prong2, tofNSigmaKa2, tofNSigmaKa);
        CHECK_AND_FILL_VEC_LC_FULL(prong2, tofNSigmaPi2, tofNSigmaPi);
        // Combined PID variables
        CHECK_AND_FILL_VEC_LC_FULL(prong0, tpcTofNSigmaPi0, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_LC_FULL(prong1, tpcTofNSigmaPi1, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_LC_FULL(prong2, tpcTofNSigmaPi2, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_LC_FULL(prong0, tpcTofNSigmaKa0, tpcTofNSigmaKa);
        CHECK_AND_FILL_VEC_LC_FULL(prong1, tpcTofNSigmaKa1, tpcTofNSigmaKa);
        CHECK_AND_FILL_VEC_LC_FULL(prong2, tpcTofNSigmaKa2, tpcTofNSigmaKa);
        CHECK_AND_FILL_VEC_LC_FULL(prong0, tpcTofNSigmaPr0, tpcTofNSigmaPr);
        CHECK_AND_FILL_VEC_LC_FULL(prong1, tpcTofNSigmaPr1, tpcTofNSigmaPr);
        CHECK_AND_FILL_VEC_LC_FULL(prong2, tpcTofNSigmaPr2, tpcTofNSigmaPr);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_LC(ptProng0),
      FILL_MAP_LC(ptProng1),
      FILL_MAP_LC(ptProng2),
      FILL_MAP_LC(impactParameterXY0),
      FILL_MAP_LC(impactParameterXY1),
      FILL_MAP_LC(impactParameterXY2),
      FILL_MAP_LC(decayLength),
      FILL_MAP_LC(decayLengthXY),
      FILL_MAP_LC(decayLengthXYNormalised),
      FILL_MAP_LC(cpa),
      FILL_MAP_LC(cpaXY),
      // TPC PID variables
      FILL_MAP_LC(tpcNSigmaP0),
      FILL_MAP_LC(tpcNSigmaKa0),
      FILL_MAP_LC(tpcNSigmaPi0),
      FILL_MAP_LC(tpcNSigmaP1),
      FILL_MAP_LC(tpcNSigmaKa1),
      FILL_MAP_LC(tpcNSigmaPi1),
      FILL_MAP_LC(tpcNSigmaP2),
      FILL_MAP_LC(tpcNSigmaKa2),
      FILL_MAP_LC(tpcNSigmaPi2),
      // TOF PID variables
      FILL_MAP_LC(tofNSigmaP0),
      FILL_MAP_LC(tofNSigmaKa0),
      FILL_MAP_LC(tofNSigmaPi0),
      FILL_MAP_LC(tofNSigmaP1),
      FILL_MAP_LC(tofNSigmaKa1),
      FILL_MAP_LC(tofNSigmaPi1),
      FILL_MAP_LC(tofNSigmaP2),
      FILL_MAP_LC(tofNSigmaKa2),
      FILL_MAP_LC(tofNSigmaPi2),
      // Combined PID variables
      FILL_MAP_LC(tpcTofNSigmaPi0),
      FILL_MAP_LC(tpcTofNSigmaPi1),
      FILL_MAP_LC(tpcTofNSigmaPi2),
      FILL_MAP_LC(tpcTofNSigmaKa0),
      FILL_MAP_LC(tpcTofNSigmaKa1),
      FILL_MAP_LC(tpcTofNSigmaKa2),
      FILL_MAP_LC(tpcTofNSigmaPr0),
      FILL_MAP_LC(tpcTofNSigmaPr1),
      FILL_MAP_LC(tpcTofNSigmaPr2)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_LC
#undef CHECK_AND_FILL_VEC_LC_FULL
#undef CHECK_AND_FILL_VEC_LC
#undef CHECK_AND_FILL_VEC_LC_HFHELPER

#endif // PWGHF_CORE_HFMLRESPONSELCTOPKPI_H_
