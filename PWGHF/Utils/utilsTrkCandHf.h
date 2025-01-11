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

/// \file utilsTrkCandHf.h
/// \brief Utility for track and candidate selections for HF analyses
/// \author Mattia Faggin <mfaggin@cern.ch>, CERN

#ifndef PWGHF_UTILS_UTILSTRKCANDHF_H_
#define PWGHF_UTILS_UTILSTRKCANDHF_H_

#include "Framework/HistogramSpec.h"

namespace o2::hf_trkcandsel
{
// secondary vertex fitting monitoring
enum SVFitting {
  BeforeFit = 0,
  FitOk,
  Fail,
  NCases
};

o2::framework::AxisSpec axisCands = {SVFitting::NCases, -0.5f, static_cast<float>(SVFitting::NCases) - 0.5f, ""};

/// @brief Function to put labels on candidate monitoring histogram
/// \param hCandidates is the histogram
template <typename Histo>
void setLabelHistoCands(Histo& hCandidates)
{
  hCandidates->SetTitle("HF candidate counter;;candidates");
  hCandidates->GetXaxis()->SetBinLabel(SVFitting::BeforeFit + 1, "Before secondary vertexing");
  hCandidates->GetXaxis()->SetBinLabel(SVFitting::FitOk + 1, "With secondary vertex");
  hCandidates->GetXaxis()->SetBinLabel(SVFitting::Fail + 1, "Run-time error in secondary vertexing");
}

/// @brief Function to evaluate number of ones in a binary representation of the argument
/// \param num is the input argument
int countOnesInBinary(uint8_t num)
{
  int count = 0;

  for (int iBit = 0; iBit < 8; iBit++) {
    count += TESTBIT(num, iBit);
  }

  return count;
}

} // namespace o2::hf_trkcandsel

#endif // PWGHF_UTILS_UTILSTRKCANDHF_H_
