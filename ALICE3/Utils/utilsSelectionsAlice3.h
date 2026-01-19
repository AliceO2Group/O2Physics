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

/// \file utilsSelectionsAlice3s.h
/// \brief Default pT bins and cut arrays for selections in ALICE3 performance analysis tasks
///
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Polytechnic University of Turin and INFN Turin

#ifndef ALICE3_UTILS_UTILSSELECTIONSALICE3_H_
#define ALICE3_UTILS_UTILSSELECTIONSALICE3_H_

#include <string> // std::string
#include <vector> // std::vector

namespace o2::analysis
{
namespace hf_cuts_3prongs_alice3
{
static constexpr int NBinsPt = 10;
static constexpr int NCutVars = 10;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0.,
  1.,
  2.,
  3.,
  4.,
  5.,
  6.,
  8.,
  12.,
  24.,
  36.};
const auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts                m,  ptP, ptK, ptPi, chi2PCA, cosp, dL, dLXY, NdLXY, ImpParXY
constexpr double Cuts[NBinsPt][NCutVars] = {{0.4, 0.4, 0.4, 0.4, 10000., 0.005, 0., 0., 0., 1e+10},  /* 0  < pT < 1  */
                                            {0.4, 0.4, 0.4, 0.4, 10000., 0.005, 0., 0., 0., 1e+10},  /* 1  < pT < 2  */
                                            {0.4, 0.4, 0.4, 0.4, 10000., 0.005, 0., 0., 0., 1e+10},  /* 2  < pT < 3  */
                                            {0.4, 0.4, 0.4, 0.4, 10000., 0.005, 0., 0., 0., 1e+10},  /* 3  < pT < 4  */
                                            {0.4, 0.4, 0.4, 0.4, 10000., 0.005, 0., 0., 0., 1e+10},  /* 4  < pT < 5  */
                                            {0.4, 0.4, 0.4, 0.4, 10000., 0.005, 0., 0., 0., 1e+10},  /* 5  < pT < 6  */
                                            {0.4, 0.4, 0.4, 0.4, 10000., 0.005, 0., 0., 0., 1e+10},  /* 6  < pT < 8  */
                                            {0.4, 0.4, 0.4, 0.4, 10000., 0.005, 0., 0., 0., 1e+10},  /* 8  < pT < 12 */
                                            {0.4, 0.4, 0.4, 0.4, 10000., 0.005, 0., 0., 0., 1e+10},  /* 12 < pT < 24 */
                                            {0.4, 0.4, 0.4, 0.4, 10000., 0.005, 0., 0., 0., 1e+10}}; /* 24 < pT < 36 */

// row labels
static const std::vector<std::string> labelsPt = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9"};

// column labels
static const std::vector<std::string> labelsCutVar = {"m", "pT prong 0", "pT prong 1", "pT prong 2", "Chi2PCA", "cos pointing angle", "decay length", "decLengthXY", "normDecLXY", "impParXY"};
} // namespace hf_cuts_3prongs_alice3

} // namespace o2::analysis

#endif // ALICE3_UTILS_UTILSSELECTIONSALICE3_H_
