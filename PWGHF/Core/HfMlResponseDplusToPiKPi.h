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

namespace o2::analysis
{

enum InputFeaturesDplusToPiKPi : int8_t {
  ptProng0 = 0,
  ptProng1,
  ptProng2,
  impactParameter0,
  impactParameter1,
  impactParameter2,
  decayLength,
  decayLengthXYNormalised,
  cpa,
  cpaXY,
  maxNormalisedDeltaIP,
  tpcTofNSigmaPi0,
  tpcTofNSigmaPi1,
  tpcTofNSigmaPi2,
  tpcTofNSigmaKa0,
  tpcTofNSigmaKa1,
  tpcTofNSigmaKa2
};

template <typename T = float, typename U = InputFeaturesDplusToPiKPi>
class HfMlResponseDplusToPiKPi : public HfMlResponse<T, U>
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
  std::vector<T> getInputFeatures(T1 const& candidate, T2 const& prong0, T2 const& prong1, T2 const& prong2) const
  {
    std::vector<T> inputFeatures;

    for (const auto& idx : MlResponse<T, U>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC(ptProng0);
        CHECK_AND_FILL_VEC(ptProng1);
        CHECK_AND_FILL_VEC(ptProng2);
        CHECK_AND_FILL_VEC(impactParameter0);
        CHECK_AND_FILL_VEC(impactParameter1);
        CHECK_AND_FILL_VEC(impactParameter2);
        CHECK_AND_FILL_VEC(decayLength);
        CHECK_AND_FILL_VEC(decayLengthXYNormalised);
        CHECK_AND_FILL_VEC(cpa);
        CHECK_AND_FILL_VEC(cpaXY);
        CHECK_AND_FILL_VEC(maxNormalisedDeltaIP);
        CHECK_AND_FILL_VEC_FULL(prong0, tpcTofNSigmaPi0, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_FULL(prong1, tpcTofNSigmaPi1, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_FULL(prong2, tpcTofNSigmaPi2, tpcTofNSigmaPi);
        CHECK_AND_FILL_VEC_FULL(prong0, tpcTofNSigmaKa0, tpcTofNSigmaKa);
        CHECK_AND_FILL_VEC_FULL(prong1, tpcTofNSigmaKa1, tpcTofNSigmaKa);
        CHECK_AND_FILL_VEC_FULL(prong2, tpcTofNSigmaKa2, tpcTofNSigmaKa);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<T, U>::mAvailableInputFeatures = {
      FILL_MAP(ptProng0),
      FILL_MAP(ptProng1),
      FILL_MAP(ptProng2),
      FILL_MAP(impactParameter0),
      FILL_MAP(impactParameter1),
      FILL_MAP(impactParameter2),
      FILL_MAP(decayLength),
      FILL_MAP(decayLengthXYNormalised),
      FILL_MAP(cpa),
      FILL_MAP(cpaXY),
      FILL_MAP(maxNormalisedDeltaIP),
      FILL_MAP(tpcTofNSigmaPi0),
      FILL_MAP(tpcTofNSigmaPi1),
      FILL_MAP(tpcTofNSigmaPi2),
      FILL_MAP(tpcTofNSigmaKa0),
      FILL_MAP(tpcTofNSigmaKa1),
      FILL_MAP(tpcTofNSigmaKa2)};
  }
};

} // namespace o2::analysis

#endif // PWGHF_CORE_HFMLRESPONSEDPLUSTOPIKPI_H_
