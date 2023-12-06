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

/// \file HfMlResponseXicToPKPi.h
/// \brief Class to compute the ML response for Xic+ → p K- π+ (based on task for  D± → π± K∓ π±) analysis selections
/// \author Cristina Terrevoli

#ifndef PWGHF_CORE_HFMLRESPONSEXICTOPKPI_H_
#define PWGHF_CORE_HFMLRESPONSEXICTOPKPI_H_

#include <map>
#include <string>
#include <vector>

#include "PWGHF/Core/HfMlResponse.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_XIC(FEATURE)                                         \
  {                                                                   \
#FEATURE, static_cast < uint8_t>(InputFeaturesXicToPKPi::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_XIC_FULL(OBJECT, FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesXicToPKPi::FEATURE): { \
    inputFeatures.emplace_back(OBJECT.GETTER());                \
    break;                                                      \
  }

// Specific case of CHECK_AND_FILL_VEC_XIC_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_XIC(GETTER)                         \
  case static_cast<uint8_t>(InputFeaturesXicToPKPi::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());            \
    break;                                                     \
  }

namespace o2::analysis
{
enum class InputFeaturesXicToPKPi : uint8_t {
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
  tofNSigmaPi2
};

template <typename TypeOutputScore = float>
class HfMlResponseXicToPKPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseXicToPKPi() = default;
  /// Default destructor
  virtual ~HfMlResponseXicToPKPi() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Xic candidate
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
        CHECK_AND_FILL_VEC_XIC(ptProng0);
        CHECK_AND_FILL_VEC_XIC(ptProng1);
        CHECK_AND_FILL_VEC_XIC(ptProng2);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, impactParameterXY0, impactParameter0);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, impactParameterXY1, impactParameter1);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, impactParameterXY2, impactParameter2);
        CHECK_AND_FILL_VEC_XIC(decayLength);
        CHECK_AND_FILL_VEC_XIC(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_XIC(cpa);
        CHECK_AND_FILL_VEC_XIC(cpaXY);
        // TPC PID variables
        CHECK_AND_FILL_VEC_XIC_FULL(prong0, tpcNSigmaP0, tpcNSigmaPr);
        CHECK_AND_FILL_VEC_XIC_FULL(prong0, tpcNSigmaKa0, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_XIC_FULL(prong0, tpcNSigmaPi0, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_XIC_FULL(prong1, tpcNSigmaP1, tpcNSigmaPr);
        CHECK_AND_FILL_VEC_XIC_FULL(prong1, tpcNSigmaKa1, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_XIC_FULL(prong1, tpcNSigmaPi1, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_XIC_FULL(prong2, tpcNSigmaP2, tpcNSigmaPr);
        CHECK_AND_FILL_VEC_XIC_FULL(prong2, tpcNSigmaKa2, tpcNSigmaKa);
        CHECK_AND_FILL_VEC_XIC_FULL(prong2, tpcNSigmaPi2, tpcNSigmaPi);
        // TOF PID variables
        CHECK_AND_FILL_VEC_XIC_FULL(prong0, tofNSigmaP0, tofNSigmaPr);
        CHECK_AND_FILL_VEC_XIC_FULL(prong0, tofNSigmaKa0, tofNSigmaKa);
        CHECK_AND_FILL_VEC_XIC_FULL(prong0, tofNSigmaPi0, tofNSigmaPi);
        CHECK_AND_FILL_VEC_XIC_FULL(prong1, tofNSigmaP1, tofNSigmaPr);
        CHECK_AND_FILL_VEC_XIC_FULL(prong1, tofNSigmaKa1, tofNSigmaKa);
        CHECK_AND_FILL_VEC_XIC_FULL(prong1, tofNSigmaPi1, tofNSigmaPi);
        CHECK_AND_FILL_VEC_XIC_FULL(prong2, tofNSigmaP2, tofNSigmaPr);
        CHECK_AND_FILL_VEC_XIC_FULL(prong2, tofNSigmaKa2, tofNSigmaKa);
        CHECK_AND_FILL_VEC_XIC_FULL(prong2, tofNSigmaPi2, tofNSigmaPi);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_XIC(ptProng0),
      FILL_MAP_XIC(ptProng1),
      FILL_MAP_XIC(ptProng2),
      FILL_MAP_XIC(impactParameterXY0),
      FILL_MAP_XIC(impactParameterXY1),
      FILL_MAP_XIC(impactParameterXY2),
      FILL_MAP_XIC(decayLength),
      FILL_MAP_XIC(decayLengthXYNormalised),
      FILL_MAP_XIC(cpa),
      FILL_MAP_XIC(cpaXY),
      // TPC PID variables
      FILL_MAP_XIC(tpcNSigmaP0),
      FILL_MAP_XIC(tpcNSigmaKa0),
      FILL_MAP_XIC(tpcNSigmaPi0),
      FILL_MAP_XIC(tpcNSigmaP1),
      FILL_MAP_XIC(tpcNSigmaKa1),
      FILL_MAP_XIC(tpcNSigmaPi1),
      FILL_MAP_XIC(tpcNSigmaP2),
      FILL_MAP_XIC(tpcNSigmaKa2),
      FILL_MAP_XIC(tpcNSigmaPi2),
      // TOF PID variables
      FILL_MAP_XIC(tofNSigmaP0),
      FILL_MAP_XIC(tofNSigmaKa0),
      FILL_MAP_XIC(tofNSigmaPi0),
      FILL_MAP_XIC(tofNSigmaP1),
      FILL_MAP_XIC(tofNSigmaKa1),
      FILL_MAP_XIC(tofNSigmaPi1),
      FILL_MAP_XIC(tofNSigmaP2),
      FILL_MAP_XIC(tofNSigmaKa2),
      FILL_MAP_XIC(tofNSigmaPi2)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_XIC
#undef CHECK_AND_FILL_VEC_XIC_FULL
#undef CHECK_AND_FILL_VEC_XIC

#endif // PWGHF_CORE_HFMLRESPONSEXICTOPKPI_H_
