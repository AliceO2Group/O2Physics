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
//
/// \file SelectorCutsRedDataFormat.h
/// \brief Default pT bins and cut arrays for heavy-flavour selectors and analysis tasks
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, Università degli Studi di Torino

// namespace with D selections for reduced charmed-resonances analysis

#ifndef PWGHF_D2H_CORE_SELECTORCUTSREDDATAFORMAT_H_
#define PWGHF_D2H_CORE_SELECTORCUTSREDDATAFORMAT_H_

#include <string> // std::string
#include <vector> // std::vector

namespace o2::analysis
{
namespace hf_cuts_d_daughter
{
static constexpr int NBinsPt = 7;
static constexpr int NCutVars = 6;
constexpr double BinsPt[NBinsPt + 1] = {
  1.,
  2.,
  4.,
  6.,
  8.,
  12.,
  24.,
  1000.};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};
// default values for the cuts
constexpr double Cuts[NBinsPt][NCutVars] = {{1.84, 1.89, 1.77, 1.81, 1.92, 1.96},  /* 1   < pt < 2 */
                                            {1.84, 1.89, 1.77, 1.81, 1.92, 1.96},  /* 2 < pt < 4 */
                                            {1.84, 1.89, 1.77, 1.81, 1.92, 1.96},  /* 4   < pt < 6 */
                                            {1.84, 1.89, 1.77, 1.81, 1.92, 1.96},  /* 6 < pt < 8 */
                                            {1.84, 1.89, 1.77, 1.81, 1.92, 1.96},  /* 8 < pt < 12 */
                                            {1.84, 1.89, 1.77, 1.81, 1.92, 1.96},  /* 12   < pt < 24 */
                                            {1.84, 1.89, 1.77, 1.81, 1.92, 1.96}}; /* 24   < pt < 1000 */
// row labels
static const std::vector<std::string> labelsPt{
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6"};
// column labels
static const std::vector<std::string> labelsCutVar = {"invMassSignalLow", "invMassSignalHigh", "invMassLeftSBLow", "invMassLeftSBHigh", "invMassRightSBLow", "invMassRightSBHigh"};
} // namespace hf_cuts_d_daughter

// namespace with v0 selections for reduced charmed-resonances analysis
namespace hf_cuts_v0_daughter
{
static constexpr int NBinsPt = 7;
static constexpr int NCutVars = 5;
constexpr double BinsPt[NBinsPt + 1] = {
  0.,
  1.,
  2.,
  4.,
  8.,
  12.,
  24.,
  1000.};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};
// default values for the cuts
constexpr double Cuts[NBinsPt][NCutVars] = {{0.48, 0.52, 0.99, 1., 0.9},  /* 1   < pt < 2 */
                                            {0.48, 0.52, 0.99, 1., 0.9},  /* 2 < pt < 4 */
                                            {0.48, 0.52, 0.99, 1., 0.9},  /* 4   < pt < 6 */
                                            {0.48, 0.52, 0.99, 1., 0.9},  /* 6 < pt < 8 */
                                            {0.48, 0.52, 0.99, 1., 0.9},  /* 8 < pt < 12 */
                                            {0.48, 0.52, 0.99, 1., 0.9},  /* 12   < pt < 24 */
                                            {0.48, 0.52, 0.99, 1., 0.9}}; /* 24   < pt < 1000 */
// row labels
static const std::vector<std::string> labelsPt{};
// column labels
static const std::vector<std::string> labelsCutVar = {"invMassLow", "invMassHigh", "cpaMin", "dcaMax", "radiusMin"};
} // namespace hf_cuts_v0_daughter
} // namespace o2::analysis
#endif // PWGHF_D2H_CORE_SELECTORCUTSREDDATAFORMAT_H_
