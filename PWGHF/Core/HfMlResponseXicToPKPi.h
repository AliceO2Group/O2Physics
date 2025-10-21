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

#include "PWGHF/Core/HfMlResponse.h"

#include "Tools/ML/MlResponse.h"

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_XIC(FEATURE)                                       \
  {                                                                 \
    #FEATURE, static_cast<uint8_t>(InputFeaturesXicToPKPi::FEATURE) \
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

// Variation of CHECK_AND_FILL_VEC_XIC_OBJECT_SIGNED(OBJECT1, OBJECT2, FEATURE, GETTER)
// where OBJECT1 and OBJECT2 are the objects from which we call the GETTER method, and the variable
// is filled depending on whether it is a XicToPKPi or a XicToPiKP
#define CHECK_AND_FILL_VEC_XIC_OBJECT_SIGNED(OBJECT1, OBJECT2, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesXicToPKPi::FEATURE): {                 \
    if (caseXicToPKPi) {                                                        \
      inputFeatures.emplace_back(OBJECT1.GETTER());                             \
    } else {                                                                    \
      inputFeatures.emplace_back(OBJECT2.GETTER());                             \
    }                                                                           \
    break;                                                                      \
  }

// Variation of CHECK_AND_FILL_VEC_XIC_OBJECT_SIGNED(OBJECT, FEATURE, GETTER1, GETTER2)
// where GETTER1 and GETTER2 are methods of the OBJECT, and used
// depending on whether the candidate is a XicToPKPi or a XicToPiKP
#define CHECK_AND_FILL_VEC_XIC_SIGNED(OBJECT, FEATURE, GETTER1, GETTER2) \
  case static_cast<uint8_t>(InputFeaturesXicToPKPi::FEATURE): {          \
    if (caseXicToPKPi) {                                                 \
      inputFeatures.emplace_back(OBJECT.GETTER1());                      \
    } else {                                                             \
      inputFeatures.emplace_back(OBJECT.GETTER2());                      \
    }                                                                    \
    break;                                                               \
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
  tpcTofNSigmaPiExpPi2
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
  template <typename T1>
  std::vector<float> getInputFeatures(T1 const& candidate, bool const caseXicToPKPi)
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
        CHECK_AND_FILL_VEC_XIC(impactParameterZ0);
        CHECK_AND_FILL_VEC_XIC(impactParameterZ1);
        CHECK_AND_FILL_VEC_XIC(impactParameterZ2);
        CHECK_AND_FILL_VEC_XIC(decayLength);
        CHECK_AND_FILL_VEC_XIC(decayLengthXY);
        CHECK_AND_FILL_VEC_XIC(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC_XIC(cpa);
        CHECK_AND_FILL_VEC_XIC(cpaXY);
        CHECK_AND_FILL_VEC_XIC(chi2PCA);
        // TPC PID variables
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcNSigmaPr0, nSigTpcPr0);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcNSigmaKa0, nSigTpcKa0);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcNSigmaPi0, nSigTpcPi0);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcNSigmaPr1, nSigTpcPr1);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcNSigmaKa1, nSigTpcKa1);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcNSigmaPi1, nSigTpcPi1);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcNSigmaPr2, nSigTpcPr2);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcNSigmaKa2, nSigTpcKa2);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcNSigmaPi2, nSigTpcPi2);
        CHECK_AND_FILL_VEC_XIC_SIGNED(candidate, tpcNSigmaPrExpPr0, nSigTpcPr0, nSigTpcPr2);
        CHECK_AND_FILL_VEC_XIC_SIGNED(candidate, tpcNSigmaPiExpPi2, nSigTpcPi2, nSigTpcPi0);
        // TOF PID variables
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tofNSigmaPr0, nSigTofPr0);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tofNSigmaKa0, nSigTofKa0);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tofNSigmaPi0, nSigTofPi0);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tofNSigmaPr1, nSigTofPr1);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tofNSigmaKa1, nSigTofKa1);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tofNSigmaPi1, nSigTofPi1);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tofNSigmaPr2, nSigTofPr2);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tofNSigmaKa2, nSigTofKa2);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tofNSigmaPi2, nSigTofPi2);
        CHECK_AND_FILL_VEC_XIC_SIGNED(candidate, tofNSigmaPrExpPr0, nSigTofPr0, nSigTofPr2);
        CHECK_AND_FILL_VEC_XIC_SIGNED(candidate, tofNSigmaPiExpPi2, nSigTofPi2, nSigTofPi0);
        // Combined PID variables
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcTofNSigmaPi0, tpcTofNSigmaPi0);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcTofNSigmaPi1, tpcTofNSigmaPi1);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcTofNSigmaPi2, tpcTofNSigmaPi2);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcTofNSigmaKa0, tpcTofNSigmaKa0);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcTofNSigmaKa1, tpcTofNSigmaKa1);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcTofNSigmaKa2, tpcTofNSigmaKa2);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcTofNSigmaPr0, tpcTofNSigmaPr0);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcTofNSigmaPr1, tpcTofNSigmaPr1);
        CHECK_AND_FILL_VEC_XIC_FULL(candidate, tpcTofNSigmaPr2, tpcTofNSigmaPr2);
        CHECK_AND_FILL_VEC_XIC_SIGNED(candidate, tpcTofNSigmaPrExpPr0, tpcTofNSigmaPr0, tpcTofNSigmaPr2);
        CHECK_AND_FILL_VEC_XIC_SIGNED(candidate, tpcTofNSigmaPiExpPi2, tpcTofNSigmaPi2, tpcTofNSigmaPi0);
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
      FILL_MAP_XIC(impactParameterZ0),
      FILL_MAP_XIC(impactParameterZ1),
      FILL_MAP_XIC(impactParameterZ2),
      FILL_MAP_XIC(decayLength),
      FILL_MAP_XIC(decayLengthXY),
      FILL_MAP_XIC(decayLengthXYNormalised),
      FILL_MAP_XIC(cpa),
      FILL_MAP_XIC(cpaXY),
      FILL_MAP_XIC(chi2PCA),
      // TPC PID variables
      FILL_MAP_XIC(tpcNSigmaPr0),
      FILL_MAP_XIC(tpcNSigmaKa0),
      FILL_MAP_XIC(tpcNSigmaPi0),
      FILL_MAP_XIC(tpcNSigmaPr1),
      FILL_MAP_XIC(tpcNSigmaKa1),
      FILL_MAP_XIC(tpcNSigmaPi1),
      FILL_MAP_XIC(tpcNSigmaPr2),
      FILL_MAP_XIC(tpcNSigmaKa2),
      FILL_MAP_XIC(tpcNSigmaPi2),
      FILL_MAP_XIC(tpcNSigmaPrExpPr0),
      FILL_MAP_XIC(tpcNSigmaPiExpPi2),
      // TOF PID variables
      FILL_MAP_XIC(tofNSigmaPr0),
      FILL_MAP_XIC(tofNSigmaKa0),
      FILL_MAP_XIC(tofNSigmaPi0),
      FILL_MAP_XIC(tofNSigmaPr1),
      FILL_MAP_XIC(tofNSigmaKa1),
      FILL_MAP_XIC(tofNSigmaPi1),
      FILL_MAP_XIC(tofNSigmaPr2),
      FILL_MAP_XIC(tofNSigmaKa2),
      FILL_MAP_XIC(tofNSigmaPi2),
      FILL_MAP_XIC(tofNSigmaPrExpPr0),
      FILL_MAP_XIC(tofNSigmaPiExpPi2),
      // Combined PID variables
      FILL_MAP_XIC(tpcTofNSigmaPi0),
      FILL_MAP_XIC(tpcTofNSigmaPi1),
      FILL_MAP_XIC(tpcTofNSigmaPi2),
      FILL_MAP_XIC(tpcTofNSigmaKa0),
      FILL_MAP_XIC(tpcTofNSigmaKa1),
      FILL_MAP_XIC(tpcTofNSigmaKa2),
      FILL_MAP_XIC(tpcTofNSigmaPr0),
      FILL_MAP_XIC(tpcTofNSigmaPr1),
      FILL_MAP_XIC(tpcTofNSigmaPr2),
      FILL_MAP_XIC(tpcTofNSigmaPrExpPr0),
      FILL_MAP_XIC(tpcTofNSigmaPiExpPi2)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_XIC
#undef CHECK_AND_FILL_VEC_XIC_FULL
#undef CHECK_AND_FILL_VEC_XIC
#undef CHECK_AND_FILL_VEC_XIC_OBJECT_SIGNED

#endif // PWGHF_CORE_HFMLRESPONSEXICTOPKPI_H_
