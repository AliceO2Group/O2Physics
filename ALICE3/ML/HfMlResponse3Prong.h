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

/// \file Alice3MlRResponse3Prong.h
/// \brief Class to compute the ML response for HF 3-prong candidates
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Polytechnic University of Turin and INFN Turin

#ifndef ALICE3_ML_HFMLRESPONSE3PRONG_H_
#define ALICE3_ML_HFMLRESPONSE3PRONG_H_

#include "Tools/ML/MlResponse.h"

#include <cstdint>
#include <map>
#include <string>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_3PRONG(FEATURE)                                 \
  {                                                              \
    #FEATURE, static_cast<uint8_t>(InputFeatures3Prong::FEATURE) \
  }

// Specific case of CHECK_AND_FILL_ML_ALICE3_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_ML_ALICE3(GETTER)                    \
  case static_cast<uint8_t>(InputFeatures3Prong::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());         \
    break;                                                  \
  }

namespace o2::analysis
{
enum class InputFeatures3Prong : uint8_t {
  ptProng0 = 0,
  ptProng1,
  ptProng2,
  impactParameterY0,
  impactParameterY1,
  impactParameterY2,
  impactParameterZ0,
  impactParameterZ1,
  impactParameterZ2,
  decayLength,
  decayLengthXY,
  decayLengthXYNormalised,
  cpa,
  cpaXY,
  chi2PCA,
  nSigRichPr0,   // 0
  nSigRichKa0,   // 0
  nSigRichPi0,   // 0
  nSigRichPr1,   // 1
  nSigRichKa1,   // 1
  nSigRichPi1,   // 1
  nSigRichPr2,   // 2
  nSigRichKa2,   // 2
  nSigRichPi2,   // 2
  nSigInnTofPr0, // 0
  nSigInnTofKa0, // 0
  nSigInnTofPi0, // 0
  nSigInnTofPr1, // 1
  nSigInnTofKa1, // 1
  nSigInnTofPi1, // 1
  nSigInnTofPr2, // 2
  nSigInnTofKa2, // 2
  nSigInnTofPi2, // 2
  nSigOutTofPr0, // 0
  nSigOutTofKa0, // 0
  nSigOutTofPi0, // 0
  nSigOutTofPr1, // 1
  nSigOutTofKa1, // 1
  nSigOutTofPi1, // 1
  nSigOutTofPr2, // 2
  nSigOutTofKa2, // 2
  nSigOutTofPi2, // 2
  nSigTrkPr0,    // 0
  nSigTrkKa0,    // 0
  nSigTrkPi0,    // 0
  nSigTrkPr1,    // 1
  nSigTrkKa1,    // 1
  nSigTrkPi1,    // 1
  nSigTrkPr2,    // 2
  nSigTrkKa2,    // 2
  nSigTrkPi2     // 2
};

template <typename TypeOutputScore = float>
class HfMlResponse3Prong : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponse3Prong() = default;
  /// Default destructor
  virtual ~HfMlResponse3Prong() = default;

  /// Method to get the input features vector needed for ML inference
  /// \tparam T1 type of the 3-prong candidate
  /// \param candidate is the 3-prong candidate
  /// \return inputFeatures vector
  template <typename T1>
  std::vector<float> getInputFeatures(T1 const& candidate)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_ML_ALICE3(ptProng0);
        CHECK_AND_FILL_ML_ALICE3(ptProng1);
        CHECK_AND_FILL_ML_ALICE3(ptProng2);
        CHECK_AND_FILL_ML_ALICE3(impactParameterY0);
        CHECK_AND_FILL_ML_ALICE3(impactParameterY1);
        CHECK_AND_FILL_ML_ALICE3(impactParameterY2);
        CHECK_AND_FILL_ML_ALICE3(impactParameterZ0);
        CHECK_AND_FILL_ML_ALICE3(impactParameterZ1);
        CHECK_AND_FILL_ML_ALICE3(impactParameterZ2);
        CHECK_AND_FILL_ML_ALICE3(decayLength);
        CHECK_AND_FILL_ML_ALICE3(decayLengthXY);
        CHECK_AND_FILL_ML_ALICE3(decayLengthXYNormalised);
        CHECK_AND_FILL_ML_ALICE3(cpa);
        CHECK_AND_FILL_ML_ALICE3(cpaXY);
        CHECK_AND_FILL_ML_ALICE3(chi2PCA);
        // TRACKER PID variables
        CHECK_AND_FILL_ML_ALICE3(nSigTrkPr0);
        CHECK_AND_FILL_ML_ALICE3(nSigTrkKa1);
        CHECK_AND_FILL_ML_ALICE3(nSigTrkPi2);
        // RICH PID variables
        CHECK_AND_FILL_ML_ALICE3(nSigRichPr0);
        CHECK_AND_FILL_ML_ALICE3(nSigRichKa1);
        CHECK_AND_FILL_ML_ALICE3(nSigRichPi2);
        // INNER TOF PID variables
        CHECK_AND_FILL_ML_ALICE3(nSigInnTofPr0);
        CHECK_AND_FILL_ML_ALICE3(nSigInnTofKa1);
        CHECK_AND_FILL_ML_ALICE3(nSigInnTofPi2);
        // OUTER TOF PID variables
        CHECK_AND_FILL_ML_ALICE3(nSigOutTofPr0);
        CHECK_AND_FILL_ML_ALICE3(nSigOutTofKa1);
        CHECK_AND_FILL_ML_ALICE3(nSigOutTofPi2);
      }
    }
    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_3PRONG(ptProng0),
      FILL_MAP_3PRONG(ptProng1),
      FILL_MAP_3PRONG(ptProng2),
      FILL_MAP_3PRONG(impactParameterY0),
      FILL_MAP_3PRONG(impactParameterY1),
      FILL_MAP_3PRONG(impactParameterY2),
      FILL_MAP_3PRONG(impactParameterZ0),
      FILL_MAP_3PRONG(impactParameterZ1),
      FILL_MAP_3PRONG(impactParameterZ2),
      FILL_MAP_3PRONG(decayLength),
      FILL_MAP_3PRONG(decayLengthXY),
      FILL_MAP_3PRONG(decayLengthXYNormalised),
      FILL_MAP_3PRONG(cpa),
      FILL_MAP_3PRONG(cpaXY),
      FILL_MAP_3PRONG(chi2PCA),
      // RICH PID variables
      FILL_MAP_3PRONG(nSigRichPr0),
      FILL_MAP_3PRONG(nSigRichKa0),
      FILL_MAP_3PRONG(nSigRichPi0),
      FILL_MAP_3PRONG(nSigRichPr1),
      FILL_MAP_3PRONG(nSigRichKa1),
      FILL_MAP_3PRONG(nSigRichPi1),
      FILL_MAP_3PRONG(nSigRichPr2),
      FILL_MAP_3PRONG(nSigRichKa2),
      FILL_MAP_3PRONG(nSigRichPi2),
      // INNER TOF PID variables
      FILL_MAP_3PRONG(nSigInnTofPr0),
      FILL_MAP_3PRONG(nSigInnTofKa0),
      FILL_MAP_3PRONG(nSigInnTofPi0),
      FILL_MAP_3PRONG(nSigInnTofPr1),
      FILL_MAP_3PRONG(nSigInnTofKa1),
      FILL_MAP_3PRONG(nSigInnTofPi1),
      FILL_MAP_3PRONG(nSigInnTofPr2),
      FILL_MAP_3PRONG(nSigInnTofKa2),
      FILL_MAP_3PRONG(nSigInnTofPi2),
      // OUTER TOF PID variables
      FILL_MAP_3PRONG(nSigOutTofPr0),
      FILL_MAP_3PRONG(nSigOutTofKa0),
      FILL_MAP_3PRONG(nSigOutTofPi0),
      FILL_MAP_3PRONG(nSigOutTofPr1),
      FILL_MAP_3PRONG(nSigOutTofKa1),
      FILL_MAP_3PRONG(nSigOutTofPi1),
      FILL_MAP_3PRONG(nSigOutTofPr2),
      FILL_MAP_3PRONG(nSigOutTofKa2),
      FILL_MAP_3PRONG(nSigOutTofPi2),
      // TRACKER PID variables
      FILL_MAP_3PRONG(nSigTrkPr0),
      FILL_MAP_3PRONG(nSigTrkKa0),
      FILL_MAP_3PRONG(nSigTrkPi0),
      FILL_MAP_3PRONG(nSigTrkPr1),
      FILL_MAP_3PRONG(nSigTrkKa1),
      FILL_MAP_3PRONG(nSigTrkPi1),
      FILL_MAP_3PRONG(nSigTrkPr2),
      FILL_MAP_3PRONG(nSigTrkKa2),
      FILL_MAP_3PRONG(nSigTrkPi2)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_3PRONG
#undef CHECK_AND_FILL_ML_ALICE3

#endif // ALICE3_ML_HFMLRESPONSE3PRONG_H_
