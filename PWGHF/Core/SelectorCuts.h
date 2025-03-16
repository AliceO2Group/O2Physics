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

/// \file SelectorCuts.h
/// \brief Default pT bins and cut arrays for heavy-flavour selectors and analysis tasks
///
/// \author Anton Alkin <anton.alkin@cern.ch>, CERN

#ifndef PWGHF_CORE_SELECTORCUTS_H_
#define PWGHF_CORE_SELECTORCUTS_H_

#include <string> // std::string
#include <vector> // std::vector

namespace o2::analysis
{

// namespace per channel

namespace hf_cuts_single_track
{
static constexpr int NBinsPtTrack = 6;
static constexpr int NCutVarsTrack = 4;
// default values for the pT bin edges (can be used to configure histogram axis)
// common for any candidate type (2-prong, 3-prong)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPtTrack[NBinsPtTrack + 1] = {
  0,
  0.5,
  1.0,
  1.5,
  2.0,
  3.0,
  1000.0};
auto vecBinsPtTrack = std::vector<double>{BinsPtTrack, BinsPtTrack + NBinsPtTrack + 1};

// default values for the dca_xy and dca_z cuts of displaced tracks
constexpr double CutsTrack[NBinsPtTrack][NCutVarsTrack] = {{0.0000, 10., 0.0000, 100.},  /* 0   < pt < 0.5 */
                                                           {0.0000, 10., 0.0000, 100.},  /* 0.5 < pt < 1 */
                                                           {0.0000, 10., 0.0000, 100.},  /* 1   < pt < 1.5 */
                                                           {0.0000, 10., 0.0000, 100.},  /* 1.5 < pt < 2 */
                                                           {0.0000, 10., 0.0000, 100.},  /* 2   < pt < 3 */
                                                           {0.0000, 10., 0.0000, 100.}}; /* 3   < pt < 1000 */

// default values for the dca_xy and dca_z cuts of primary tracks (e.g. D* soft pions)
constexpr double CutsTrackPrimary[NBinsPtTrack][NCutVarsTrack] = {{0.0000, 2., 0.0000, 100.},  /* 0   < pt < 0.5 */
                                                                  {0.0000, 2., 0.0000, 100.},  /* 0.5 < pt < 1 */
                                                                  {0.0000, 2., 0.0000, 100.},  /* 1   < pt < 1.5 */
                                                                  {0.0000, 2., 0.0000, 100.},  /* 1.5 < pt < 2 */
                                                                  {0.0000, 2., 0.0000, 100.},  /* 2   < pt < 3 */
                                                                  {0.0000, 2., 0.0000, 100.}}; /* 3   < pt < 1000 */

// row labels
static const std::vector<std::string> labelsPtTrack{};

// column labels
static const std::vector<std::string> labelsCutVarTrack = {"min_dcaxytoprimary", "max_dcaxytoprimary", "min_dcaztoprimary", "max_dcaztoprimary"};
} // namespace hf_cuts_single_track

namespace hf_presel_pid
{
// default values for the PID cuts for protons in the track-index-skim-creator
constexpr float CutsPid[4][6] = {{0.f, 1000.f, 5.f, 0.f, 1000.f, 5.f},
                                 {0.f, 1000.f, 5.f, 0.f, 1000.f, 5.f},
                                 {0.f, 1000.f, 5.f, 0.f, 1000.f, 5.f},
                                 {0.f, 1000.f, 5.f, 0.f, 1000.f, 5.f}};
static const std::vector<std::string> labelsCutsPid = {"minPtTpc", "maxPtTpc", "nSigmaMaxTpc", "minPtTof", "maxPtTof", "nSigmaMaxTof"};
static const std::vector<std::string> labelsRowsPid = {"ProtonInLcToPKPi", "ProtonInXicToPKPi", "ProtonInLcToPK0S", "KaonIn3Prongs"};
} // namespace hf_presel_pid

namespace hf_cuts_bdt_multiclass
{
static constexpr int NBinsPt = 1;
static constexpr int NCutBdtScores = 3;
// default values for the pT bin edges (can be used to configure histogram axis)
// common for any charm candidate
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0.,
  1000.0};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
constexpr double Cuts[NBinsPt][NCutBdtScores] = {{0.1, 0.5, 0.5}};

// row labels
static const std::vector<std::string> labelsPt{};

// column labels
static const std::vector<std::string> labelsCutBdt = {"BDTbkg", "BDTprompt", "BDTnonprompt"};

// ONNX file names
static const std::string onnxFileNameSpecies[5][1] = {{"ModelHandler_onnx_D0ToKPi.onnx"}, {"ModelHandler_onnx_DplusToPiKPi.onnx"}, {"ModelHandler_onnx_LcToPKPi.onnx"}, {"ModelHandler_onnx_DsToKKPi.onnx"}, {""}}; // NOLINT: runtime/string

// row labels for ONNX file names
static const std::vector<std::string> labelsSpecies = {"D0", "Dplus", "Lc", "Ds", "Xic"};

// column label for ONNX file names
static const std::vector<std::string> labelsModels = {"Model path"};

} // namespace hf_cuts_bdt_multiclass

namespace hf_cuts_ml
{
// direction of the cut
enum CutDirection {
  CutGreater = 0, // require score < cut value
  CutSmaller,     // require score > cut value
  CutNot          // do not cut on score
};

static constexpr int NBinsPt = 12;
static constexpr int NCutScores = 3;
// default values for the pT bin edges, offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0.,
  1.,
  2.,
  3.,
  4.,
  5.,
  6.,
  8.,
  10.,
  12.,
  16.,
  24.,
  50.};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the ML model paths, one model per pT bin
static const std::vector<std::string> modelPaths = {
  ""};

// default values for the cut directions
constexpr int CutDir[NCutScores] = {CutGreater, CutSmaller, CutSmaller};
auto vecCutDir = std::vector<int>{CutDir, CutDir + NCutScores};

// default values for the cuts
constexpr double Cuts[NBinsPt][NCutScores] = {
  {0.5, 0.5, 0.5},
  {0.5, 0.5, 0.5},
  {0.5, 0.5, 0.5},
  {0.5, 0.5, 0.5},
  {0.5, 0.5, 0.5},
  {0.5, 0.5, 0.5},
  {0.5, 0.5, 0.5},
  {0.5, 0.5, 0.5},
  {0.5, 0.5, 0.5},
  {0.5, 0.5, 0.5},
  {0.5, 0.5, 0.5},
  {0.5, 0.5, 0.5}};

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
  "pT bin 9",
  "pT bin 10",
  "pT bin 11"};

// column labels
static const std::vector<std::string> labelsCutScore = {"score class 1", "score class 2", "score class 3"};
static const std::vector<std::string> labelsDmesCutScore = {"ML score charm bkg", "ML score charm prompt", "ML score charm nonprompt"};
} // namespace hf_cuts_ml

namespace hf_cuts_presel_2prong
{
static constexpr int NBinsPt = 2;
static constexpr int NCutVars = 4;
// default values for the pT bin edges (can be used to configure histogram axis)
// common for any 2-prong candidate
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  1.,
  5.,
  1000.0};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
constexpr double Cuts[NBinsPt][NCutVars] = {{1.65, 2.15, 0.5, 100.},  /* 1 < pt < 5 */
                                            {1.65, 2.15, 0.5, 100.}}; /* 5 < pt < 1000 */

// row labels
static const std::vector<std::string> labelsPt{};

// column labels
static const std::vector<std::string> labelsCutVar = {"massMin", "massMax", "cosp", "d0d0"};
} // namespace hf_cuts_presel_2prong

namespace hf_cuts_presel_3prong
{
static constexpr int NBinsPt = 2;
static constexpr int NCutVars = 4;
// default values for the pT bin edges (can be used to configure histogram axis)
// common for any 3-prong candidate
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  1.,
  5.,
  1000.0};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
constexpr double Cuts[NBinsPt][NCutVars] = {{1.75, 2.05, 0.7, 0.02},  /* 1 < pt < 5 */
                                            {1.75, 2.05, 0.5, 0.02}}; /* 5 < pt < 1000 */

// row labels
static const std::vector<std::string> labelsPt{};

// column labels
static const std::vector<std::string> labelsCutVar = {"massMin", "massMax", "cosp", "decL"};
} // namespace hf_cuts_presel_3prong

namespace hf_cuts_presel_ds
{
static constexpr int NBinsPt = 2;
static constexpr int NCutVars = 5;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  1.,
  5.,
  1000.0};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
constexpr double Cuts[NBinsPt][NCutVars] = {{1.70, 2.15, 0.7, 0.02, 0.02},  /* 1 < pt < 5 */
                                            {1.70, 2.15, 0.5, 0.02, 0.02}}; /* 5 < pt < 1000 */

// row labels
static const std::vector<std::string> labelsPt{};

// column labels
static const std::vector<std::string> labelsCutVar = {"massMin", "massMax", "cosp", "decL", "deltaMassKK"};
} // namespace hf_cuts_presel_ds

namespace hf_cuts_presel_dstar
{
static constexpr int NBinsPt = 2;
static constexpr int NCutVars = 2;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  1.,
  5.,
  1000.0};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
constexpr double Cuts[NBinsPt][NCutVars] = {{0.17, 0.05},  /* 1 < pt < 5 */
                                            {0.17, 0.08}}; /* 5 < pt < 1000 */

// row labels
static const std::vector<std::string> labelsPt{};

// column labels
static const std::vector<std::string> labelsCutVar = {"deltaMassMax", "deltaMassD0"};
} // namespace hf_cuts_presel_dstar

namespace hf_cuts_d0_to_pi_k
{
static constexpr int NBinsPt = 25;
static constexpr int NCutVars = 15;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0,
  0.5,
  1.0,
  1.5,
  2.0,
  2.5,
  3.0,
  3.5,
  4.0,
  4.5,
  5.0,
  5.5,
  6.0,
  6.5,
  7.0,
  7.5,
  8.0,
  9.0,
  10.0,
  12.0,
  16.0,
  20.0,
  24.0,
  36.0,
  50.0,
  100.0};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
constexpr double Cuts[NBinsPt][NCutVars] = {{0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5},   /* 0   < pT < 0.5 */
                                            {0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5},   /* 0.5 < pT < 1   */
                                            {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5},  /* 1   < pT < 1.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5},  /* 1.5 < pT < 2   */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0., 10., 10., 0.06, 0.5},  /* 2   < pT < 2.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0., 10., 10., 0.06, 0.5},  /* 2.5 < pT < 3   */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  /* 3   < pT < 3.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  /* 3.5 < pT < 4   */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 4   < pT < 4.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 4.5 < pT < 5   */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 5   < pT < 5.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 5.5 < pT < 6   */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 6   < pT < 6.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 6.5 < pT < 7   */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 7   < pT < 7.5 */
                                            {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 7.5 < pT < 8   */
                                            {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 8   < pT < 9   */
                                            {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 9   < pT < 10  */
                                            {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 10  < pT < 12  */
                                            {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 10000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 12  < pT < 16  */
                                            {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  /* 16  < pT < 20  */
                                            {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  /* 20  < pT < 24  */
                                            {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  /* 24  < pT < 36  */
                                            {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  /* 36  < pT < 50  */
                                            {0.400, 300. * 1E-4, 1.0, 0.6, 0.6, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5}}; /* 50  < pT < 100 */

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
  "pT bin 9",
  "pT bin 10",
  "pT bin 11",
  "pT bin 12",
  "pT bin 13",
  "pT bin 14",
  "pT bin 15",
  "pT bin 16",
  "pT bin 17",
  "pT bin 18",
  "pT bin 19",
  "pT bin 20",
  "pT bin 21",
  "pT bin 22",
  "pT bin 23",
  "pT bin 24"};

// column labels
static const std::vector<std::string> labelsCutVar = {"m", "DCA", "cos theta*", "pT K", "pT Pi", "d0K", "d0pi", "d0d0", "cos pointing angle", "cos pointing angle xy", "min norm decay length XY", "max decay length", "max decay length XY", "min decay length", "norm dauImpPar XY"};
} // namespace hf_cuts_d0_to_pi_k

namespace hf_cuts_dstar_to_d0_pi
{
static constexpr int NBinsPt = 25;
static constexpr int NCutVars = 8;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0.,
  0.5,
  1.0,
  1.5,
  2.0,
  2.5,
  3.0,
  3.5,
  4.0,
  4.5,
  5.0,
  5.5,
  6.0,
  6.5,
  7.0,
  7.5,
  8.0,
  9.0,
  10.0,
  12.0,
  16.0,
  20.0,
  24.0,
  36.0,
  50.0,
  100.0};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

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
  "pT bin 9",
  "pT bin 10",
  "pT bin 11",
  "pT bin 12",
  "pT bin 13",
  "pT bin 14",
  "pT bin 15",
  "pT bin 16",
  "pT bin 17",
  "pT bin 18",
  "pT bin 19",
  "pT bin 20",
  "pT bin 21",
  "pT bin 22",
  "pT bin 23",
  "pT bin 24"};

// column label
static const std::vector<std::string> labelsCutVar = {"ptSoftPiMin", "ptSoftPiMax", "d0SoftPi", "d0SoftPiNormalised", "deltaMInvDstar", "chi2PCA", "d0Prong0Normalised", "d0Prong1Normalised"};

// default values for the cuts
constexpr double Cuts[NBinsPt][NCutVars] = {{0.05, 0.2, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 0.2, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 0.3, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 0.3, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 0.4, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 0.4, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 0.6, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 0.6, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 0.6, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0},
                                            {0.05, 100, 0.1, 1000.0, 0.2, 300.0, 0.0, 0.0}};
} // namespace hf_cuts_dstar_to_d0_pi

namespace hf_cuts_lc_to_p_k_pi
{
static constexpr int NBinsPt = 10;
static constexpr int NCutVars = 11;
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
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts                m,  ptP, ptK, ptPi, chi2PCA, dL, cosp, dLXY, NdLXY, ImpParXY, mass(Kpi)
constexpr double Cuts[NBinsPt][NCutVars] = {{0.4, 0.4, 0.4, 0.4, 0., 0.005, 0., 0., 0., 1e+10, -1.},  /* 0  < pT < 1  */
                                            {0.4, 0.4, 0.4, 0.4, 0., 0.005, 0., 0., 0., 1e+10, -1.},  /* 1  < pT < 2  */
                                            {0.4, 0.4, 0.4, 0.4, 0., 0.005, 0., 0., 0., 1e+10, -1.},  /* 2  < pT < 3  */
                                            {0.4, 0.4, 0.4, 0.4, 0., 0.005, 0., 0., 0., 1e+10, -1.},  /* 3  < pT < 4  */
                                            {0.4, 0.4, 0.4, 0.4, 0., 0.005, 0., 0., 0., 1e+10, -1.},  /* 4  < pT < 5  */
                                            {0.4, 0.4, 0.4, 0.4, 0., 0.005, 0., 0., 0., 1e+10, -1.},  /* 5  < pT < 6  */
                                            {0.4, 0.4, 0.4, 0.4, 0., 0.005, 0., 0., 0., 1e+10, -1.},  /* 6  < pT < 8  */
                                            {0.4, 0.4, 0.4, 0.4, 0., 0.005, 0., 0., 0., 1e+10, -1.},  /* 8  < pT < 12 */
                                            {0.4, 0.4, 0.4, 0.4, 0., 0.005, 0., 0., 0., 1e+10, -1.},  /* 12 < pT < 24 */
                                            {0.4, 0.4, 0.4, 0.4, 0., 0.005, 0., 0., 0., 1e+10, -1.}}; /* 24 < pT < 36 */

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
static const std::vector<std::string> labelsCutVar = {"m", "pT p", "pT K", "pT Pi", "Chi2PCA", "decay length", "cos pointing angle", "decLengthXY", "normDecLXY", "impParXY", "mass (Kpi)"};
} // namespace hf_cuts_lc_to_p_k_pi

namespace hf_cuts_lc_to_k0s_p
{
static constexpr int NBinsPt = 8;
static constexpr int NCutVars = 8;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  1.,
  2.,
  3.,
  4.,
  5.,
  6.,
  8.,
  12.,
  24.};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
// mK0s(GeV)     mLambdas(GeV)    mGammas(GeV)    ptp     ptK0sdau     ptK0s     d0p     d0K0
constexpr double Cuts[NBinsPt][NCutVars] = {{0.008, 0.005, 0.1, 0.5, 0.3, 0.6, 0.05, 999999.},  // 1 < pt < 2
                                            {0.008, 0.005, 0.1, 0.5, 0.4, 1.3, 0.05, 999999.},  // 2 < pt < 3
                                            {0.009, 0.005, 0.1, 0.6, 0.4, 1.3, 0.05, 999999.},  // 3 < pt < 4
                                            {0.011, 0.005, 0.1, 0.6, 0.4, 1.4, 0.05, 999999.},  // 4 < pt < 5
                                            {0.013, 0.005, 0.1, 0.6, 0.4, 1.4, 0.06, 999999.},  // 5 < pt < 6
                                            {0.013, 0.005, 0.1, 0.9, 0.4, 1.6, 0.09, 999999.},  // 6 < pt < 8
                                            {0.016, 0.005, 0.1, 0.9, 0.4, 1.7, 0.10, 999999.},  // 8 < pt < 12
                                            {0.019, 0.005, 0.1, 1.0, 0.4, 1.9, 0.20, 999999.}}; // 12 < pt < 24

// row labels
static const std::vector<std::string> labelsPt = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7"};

// column labels
static const std::vector<std::string> labelsCutVar = {"mK0s", "mLambda", "mGamma", "ptBach", "ptV0Dau", "ptV0", "d0Bach", "d0V0"};
} // namespace hf_cuts_lc_to_k0s_p

namespace hf_cuts_dplus_to_pi_k_pi
{
static constexpr int NBinsPt = 12;
static constexpr int NCutVars = 8;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  1.,
  2.,
  3.,
  4.,
  5.,
  6.,
  7.,
  8.,
  10.,
  12.,
  16.,
  24.,
  36.};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
// selections from pp at 5 TeV 2017 analysis https://alice-notes.web.cern.ch/node/808
constexpr double Cuts[NBinsPt][NCutVars] = {{0.2, 0.3, 0.3, 0.07, 6., 0.96, 0.985, 2.5},  /* 1  < pT < 2  */
                                            {0.2, 0.3, 0.3, 0.07, 5., 0.96, 0.985, 2.5},  /* 2  < pT < 3  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.980, 2.5},  /* 3  < pT < 4  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 4  < pT < 5  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 5  < pT < 6  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 6  < pT < 7  */
                                            {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 7  < pT < 8  */
                                            {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 8  < pT < 10 */
                                            {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 10 < pT < 12 */
                                            {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 12 < pT < 16 */
                                            {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 16 < pT < 24 */
                                            {0.2, 0.3, 0.3, 0.20, 5., 0.94, 0.000, 2.5}}; /* 24 < pT < 36 */

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
  "pT bin 9",
  "pT bin 10",
  "pT bin 11"};

// column labels
static const std::vector<std::string> labelsCutVar = {"deltaM", "pT Pi", "pT K", "decay length", "normalized decay length XY", "cos pointing angle", "cos pointing angle XY", "max normalized deltaIP"};
} // namespace hf_cuts_dplus_to_pi_k_pi

namespace hf_cuts_ds_to_k_k_pi
{
static constexpr int NBinsPt = 8;
static constexpr int NCutVars = 11;
// momentary cuts
constexpr double BinsPt[NBinsPt + 1] = {
  2.,
  3.,
  4.,
  5.,
  6.,
  8.,
  12.,
  16.,
  24.};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
// selections from pp at 5 TeV 2017 analysis https://alice-notes.web.cern.ch/node/808
constexpr double Cuts[NBinsPt][NCutVars] = {{0.2, 0.3, 0.3, 0.02, 4., 0.92, 0.92, 0.014, 0.010, 0.10, 5},  /* 2  < pT < 3  */
                                            {0.2, 0.3, 0.3, 0.02, 4., 0.92, 0.92, 0.014, 0.010, 0.10, 5},  /* 3  < pT < 4  */
                                            {0.2, 0.3, 0.3, 0.03, 4., 0.90, 0.90, 0.012, 0.010, 0.05, 5},  /* 4  < pT < 5  */
                                            {0.2, 0.3, 0.3, 0.03, 4., 0.90, 0.90, 0.012, 0.010, 0.05, 5},  /* 5  < pT < 6  */
                                            {0.2, 0.3, 0.3, 0.03, 4., 0.90, 0.90, 0.012, 0.010, 0.05, 5},  /* 6  < pT < 8  */
                                            {0.2, 0.3, 0.3, 0.03, 4., 0.90, 0.90, 0.012, 0.010, 0.00, 5},  /* 8  < pT < 12 */
                                            {0.2, 0.3, 0.3, 0.05, 4., 0.85, 0.85, 0.012, 0.015, 0.00, 5},  /* 12 < pT < 16 */
                                            {0.2, 0.3, 0.3, 0.05, 4., 0.85, 0.85, 0.012, 0.015, 0.00, 5}}; /* 16 < pT < 24 */

// row labels
static const std::vector<std::string> labelsPt = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7"};

// column labels
static const std::vector<std::string> labelsCutVar = {"deltaM", "pT Pi", "pT K", "decay length", "normalized decay length XY", "cos pointing angle", "cos pointing angle XY", "impact parameter XY", "deltaM Phi", "cos^3 theta_PiK", "chi2PCA"};
} // namespace hf_cuts_ds_to_k_k_pi

namespace hf_cuts_omegac_to_omega_pi
{
static constexpr int NBinsPt = 4;
static constexpr int NCutVars = 1;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {

  1.0,
  2.0,
  4.0,
  6.0,
  12.0};

auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
// pi_pT
constexpr double Cuts[NBinsPt][NCutVars] = {{0.2},  /* 1 < pt < 2 */
                                            {0.2},  /* 2 < pt < 4 */
                                            {0.6},  /* 4 < pt < 6 */
                                            {0.8}}; /* 6 < pt < 12 */

// row labels
static const std::vector<std::string> labelsPt = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3"};

// column labels
static const std::vector<std::string> labelsCutVar = {"pT pi from Omegac"};
} // namespace hf_cuts_omegac_to_omega_pi

namespace hf_cuts_xic_to_p_k_pi
{
static constexpr int NBinsPt = 10;
static constexpr int NCutVars = 11;
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
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts                m    ptP  ptK  ptPi chi2PCA dL   cosp, dLXY, NdL, ct, ImpParXY
constexpr double Cuts[NBinsPt][NCutVars] = {{0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8, 0.005, 4., 2., 0.},  /* 0  < pT < 1  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8, 0.005, 4., 2., 0.},  /* 1  < pT < 2  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8, 0.005, 4., 2., 0.},  /* 2  < pT < 3  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8, 0.005, 4., 2., 0.},  /* 3  < pT < 4  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8, 0.005, 4., 2., 0.},  /* 4  < pT < 5  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8, 0.005, 4., 2., 0.},  /* 5  < pT < 6  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8, 0.005, 4., 2., 0.},  /* 6  < pT < 8  */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8, 0.005, 4., 2., 0.},  /* 8  < pT < 12 */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8, 0.005, 4., 2., 0.},  /* 12 < pT < 24 */
                                            {0.400, 0.4, 0.4, 0.4, 1e-5, 0.005, 0.8, 0.005, 4., 2., 0.}}; /* 24 < pT < 36 */

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
static const std::vector<std::string> labelsCutVar = {"m", "pT p", "pT K", "pT Pi", "chi2PCA", "decay length", "cos pointing angle", "decLengthXY", "normDecLXY", "ct", "impParXY"};
} // namespace hf_cuts_xic_to_p_k_pi

namespace hf_cuts_xic_to_xi_pi_pi
{
static constexpr int NBinsPt = 10;
static constexpr int NCutVars = 12;
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
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts               m   ptXi ptPi0 ptPi1 chi2PCA dL dLXY cosp cospXY impParXY Xi Pi0 Pi1
constexpr double Cuts[NBinsPt][NCutVars] = {{0.4, 0.4, 0.4, 0.4, 1e-5, 0.5, 0.5, 0.9, 0.9, 0.1, 0.1, 0.1},  /* 0  < pT < 1  */
                                            {0.4, 0.4, 0.4, 0.4, 1e-5, 0.5, 0.5, 0.9, 0.9, 0.1, 0.1, 0.1},  /* 1  < pT < 2  */
                                            {0.4, 0.4, 0.4, 0.4, 1e-5, 0.5, 0.5, 0.9, 0.9, 0.1, 0.1, 0.1},  /* 2  < pT < 3  */
                                            {0.4, 0.4, 0.4, 0.4, 1e-5, 0.5, 0.5, 0.9, 0.9, 0.1, 0.1, 0.1},  /* 3  < pT < 4  */
                                            {0.4, 0.4, 0.4, 0.4, 1e-5, 0.5, 0.5, 0.9, 0.9, 0.1, 0.1, 0.1},  /* 4  < pT < 5  */
                                            {0.4, 0.4, 0.4, 0.4, 1e-5, 0.5, 0.5, 0.9, 0.9, 0.1, 0.1, 0.1},  /* 5  < pT < 6  */
                                            {0.4, 0.4, 0.4, 0.4, 1e-5, 0.5, 0.5, 0.9, 0.9, 0.1, 0.1, 0.1},  /* 6  < pT < 8  */
                                            {0.4, 0.4, 0.4, 0.4, 1e-5, 0.5, 0.5, 0.9, 0.9, 0.1, 0.1, 0.1},  /* 8  < pT < 10 */
                                            {0.4, 0.4, 0.4, 0.4, 1e-5, 0.5, 0.5, 0.9, 0.9, 0.1, 0.1, 0.1},  /* 12 < pT < 24 */
                                            {0.4, 0.4, 0.4, 0.4, 1e-5, 0.5, 0.5, 0.9, 0.9, 0.1, 0.1, 0.1}}; /* 24 < pT < 36 */

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
static const std::vector<std::string> labelsCutVar = {"m", "pT Xi", "pT Pi0", "pT Pi1", "chi2PCA", "max decay length", "max decay length XY", "cos pointing angle", "cos pointing angle XY", "max impParXY Xi", "max impParXY Pi0", "max impParXY Pi1"};
} // namespace hf_cuts_xic_to_xi_pi_pi

namespace hf_cuts_xicc_to_p_k_pi_pi
{
static constexpr int NBinsPt = 10;
static constexpr int NCutVars = 14;
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
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
constexpr double Cuts[NBinsPt][NCutVars] = {{0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 0  < pT < 1  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 1  < pT < 2  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 2  < pT < 3  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 3  < pT < 4  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 4  < pT < 5  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 5  < pT < 6  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 6  < pT < 8  */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 8  < pT < 12 */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8},  /* 12 < pT < 24 */
                                            {0.400, 0.5, 0.2, 1.e-3, 10.0, 1.e-3, 10.0, 9999., 1.e-3, 0.0, 50.0, 50.0, 0.8, 0.8}}; /* 24 < pT < 36 */

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
static const std::vector<std::string> labelsCutVar = {"m", "pT Xic", "pT Pi", "min d0 Xic", "max d0 Xic", "min d0 Pi", "max d0 Pi", "d0d0", "chi2PCA", "min decay length", "max decay length", "max decay length XY", "cos pointing angle", "cos pointing angle XY"};
} // namespace hf_cuts_xicc_to_p_k_pi_pi

namespace hf_cuts_jpsi_to_e_e
{
static constexpr int NBinsPt = 9;
static constexpr int NCutVars = 5;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0,
  0.5,
  1.0,
  2.0,
  3.0,
  4.0,
  5.0,
  7.0,
  10.0,
  15.0,
};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
constexpr double Cuts[NBinsPt][NCutVars] = {{0.5, 0.2, 0.4, 1, 1.},  /* 0   < pT < 0.5 */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 0.5 < pT < 1   */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 1   < pT < 2   */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 2   < pT < 3   */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 3   < pT < 4   */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 4   < pT < 5   */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 5   < pT < 7   */
                                            {0.5, 0.2, 0.4, 1, 1.},  /* 7   < pT < 10  */
                                            {0.5, 0.2, 0.4, 1, 1.}}; /* 10  < pT < 15  */

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
  "pT bin 8"};

// column labels
static const std::vector<std::string> labelsCutVar = {"m", "DCA_xy", "DCA_z", "pT El", "chi2PCA"};
} // namespace hf_cuts_jpsi_to_e_e

namespace hf_cuts_b0_to_d_pi
{
static constexpr int NBinsPt = 12;
static constexpr int NCutVars = 12;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0,
  0.5,
  1.0,
  2.0,
  3.0,
  4.0,
  5.0,
  7.0,
  10.0,
  13.0,
  16.0,
  20.0,
  24.0};

auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
// DeltaM CPA chi2PCA d0D d0Pi pTD pTPi B0DecayLength B0DecayLengthXY IPProd DeltaMD CthetaStr
constexpr double Cuts[NBinsPt][NCutVars] = {{1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 0 < pt < 0.5 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 0.5 < pt < 1 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 1 < pt < 2 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 2 < pt < 3 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 3 < pt < 4 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 4 < pt < 5 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 5 < pt < 7 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 7 < pt < 10 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 10 < pt < 13 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 13 < pt < 16 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 16 < pt < 20 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8}}; /* 20 < pt < 24 */
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
  "pT bin 9",
  "pT bin 10",
  "pT bin 11"};

// column labels
static const std::vector<std::string> labelsCutVar = {"m", "CPA", "Chi2PCA", "d0 D", "d0 Pi", "pT D", "pT Pi", "B0 decLen", "B0 decLenXY", "Imp. Par. Product", "DeltaMD", "Cos ThetaStar"};
} // namespace hf_cuts_b0_to_d_pi

namespace hf_cuts_bs_to_ds_pi
{
static constexpr int NBinsPt = 10;
static constexpr int NCutVars = 10;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0,
  1.0,
  2.0,
  3.0,
  4.0,
  5.0,
  8.0,
  10.0,
  12.0,
  16.0,
  24.0};

auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
// DeltaM CPA chi2PCA d0Ds d0Pi pTDs pTPi BsDecayLength BsDecayLengthXY IPProd
constexpr double Cuts[NBinsPt][NCutVars] = {{1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0.},  /* 0 < pt < 1 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0.},  /* 1 < pt < 2 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0.},  /* 2 < pt < 3 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0.},  /* 3 < pt < 4 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0.},  /* 4 < pt < 5 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0.},  /* 5 < pt < 6 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0.},  /* 8 < pt < 10 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0.},  /* 10 < pt < 12 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0.},  /* 12 < pt < 16 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0.}}; /* 16 < pt < 24 */

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
static const std::vector<std::string> labelsCutVar = {"m", "CPA", "Chi2PCA", "d0 Ds", "d0 Pi", "pT Ds", "pT Pi", "Bs decLen", "Bs decLenXY", "Imp. Par. Product"};
} // namespace hf_cuts_bs_to_ds_pi

namespace hf_cuts_bplus_to_d0_pi
{
static constexpr int NBinsPt = 12;
static constexpr int NCutVars = 11;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0,
  0.5,
  1.0,
  2.0,
  3.0,
  4.0,
  5.0,
  7.0,
  10.0,
  13.0,
  16.0,
  20.0,
  24.0};

auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
// DeltaM CPA d0D0 d0Pi pTD0 pTPi BDecayLength BDecayLengthXY IPProd DeltaMD0 CthetaStr
constexpr double Cuts[NBinsPt][NCutVars] = {{1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 0 < pt < 0.5 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 0.5 < pt < 1 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 1 < pt < 2 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 2 < pt < 3 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 3 < pt < 4 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 4 < pt < 5 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 5 < pt < 7 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 7 < pt < 10 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 10 < pt < 13 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 13 < pt < 16 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 16 < pt < 20 */
                                            {1., 0.8, 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8}}; /* 20 < pt < 24 */
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
  "pT bin 9",
  "pT bin 10",
  "pT bin 11"};

// column labels
static const std::vector<std::string> labelsCutVar = {"m", "CPA", "d0 D0", "d0 Pi", "pT D0", "pT Pi", "B decLen", "B decLenXY", "Imp. Par. Product", "DeltaMD0", "Cos ThetaStar"};
} // namespace hf_cuts_bplus_to_d0_pi

namespace hf_cuts_lb_to_lc_pi
{
static constexpr int NBinsPt = 12;
static constexpr int NCutVars = 12;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0,
  0.5,
  1.0,
  2.0,
  3.0,
  4.0,
  5.0,
  7.0,
  10.0,
  13.0,
  16.0,
  20.0,
  24.0};

auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
// DeltaM CPA chi2PCA d0Lc d0Pi pTLc pTPi LbDecayLength LbDecayLengthXY IPProd DeltaMLc CthetaStr
constexpr double Cuts[NBinsPt][NCutVars] = {{1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 0 < pt < 0.5 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 0.5 < pt < 1 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 1 < pt < 2 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 2 < pt < 3 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 3 < pt < 4 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 4 < pt < 5 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 5 < pt < 7 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 7 < pt < 10 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 10 < pt < 13 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 13 < pt < 16 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8},  /* 16 < pt < 20 */
                                            {1., 0.8, 1., 0.01, 0.01, 1.0, 0.15, 0.05, 0.05, 0., 0.1, 0.8}}; /* 20 < pt < 24 */
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
  "pT bin 9",
  "pT bin 10",
  "pT bin 11"};

// column labels
static const std::vector<std::string> labelsCutVar = {"m", "CPA", "Chi2PCA", "d0 Lc+", "d0 Pi", "pT Lc+", "pT Pi", "Lb decLen", "Lb decLenXY", "Imp. Par. Product", "DeltaMLc", "Cos ThetaStar"};
} // namespace hf_cuts_lb_to_lc_pi

namespace hf_cuts_x_to_jpsi_pi_pi
{
static constexpr int NBinsPt = 9;
static constexpr int NCutVars = 7;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0,
  0.5,
  1.0,
  2.0,
  3.0,
  4.0,
  5.0,
  7.0,
  10.0,
  15.0,
};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
//                                            m   CPA  d0Jpsi  d0Pi pTJpsi pTPi chi2PCA
constexpr double Cuts[NBinsPt][NCutVars] = {{0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 0<pt<0.5 */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 0.5<pt<1 */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 1<pt<2   */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 2<pt<3   */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 3<pt<4   */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 4<pt<5   */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 5<pt<7   */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.},  /* 7<pt<10  */
                                            {0.5, 0.80, 0.001, 0.001, 3.0, 0.15, 1.}}; /* 10<pt<15 */
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
  "pT bin 8"};
// column labels
static const std::vector<std::string> labelsCutVar = {"m", "CPA", "d0 Jpsi", "d0 Pi", "pT Jpsi", "pT Pi", "chi2PCA"};
} // namespace hf_cuts_x_to_jpsi_pi_pi

namespace hf_cuts_chic_to_jpsi_gamma
{
// dummy selections for chic --> TO BE IMPLEMENTED
static constexpr int NBinsPt = 9;
static constexpr int NCutVars = 7;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0,
  0.5,
  1.0,
  2.0,
  3.0,
  4.0,
  5.0,
  7.0,
  10.0,
  15.0,
};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
//                                            m   CPA  d0Jpsi  d0gamma pTJpsi pTgamma chi2PCA
constexpr double Cuts[NBinsPt][NCutVars] = {{3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 0<pt<0.5 */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 0.5<pt<1 */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 1<pt<2   */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 2<pt<3   */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 3<pt<4   */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 4<pt<5   */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 5<pt<7   */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.},  /* 7<pt<10  */
                                            {3.0, -1., 0.001, 0.001, 0.5, 0.15, 1.}}; /* 10<pt<15 */
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
  "pT bin 8"};
// column labels
static const std::vector<std::string> labelsCutVar = {"m", "CPA", "d0 Jpsi", "d0 Gamma", "pT Jpsi", "pT Gamma", "chi2PCA"};
} // namespace hf_cuts_chic_to_jpsi_gamma

namespace hf_cuts_sigmac_to_p_k_pi
{
static constexpr int NBinsPt = 10;
static constexpr int NCutVars = 2;
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
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};

// default values for the cuts
constexpr double Cuts[NBinsPt][NCutVars] = {{0.03, 0.03},  /* 0  < pT < 1  */
                                            {0.03, 0.03},  /* 1  < pT < 2  */
                                            {0.03, 0.03},  /* 2  < pT < 3  */
                                            {0.03, 0.03},  /* 3  < pT < 4  */
                                            {0.03, 0.03},  /* 4  < pT < 5  */
                                            {0.03, 0.03},  /* 5  < pT < 6  */
                                            {0.03, 0.03},  /* 6  < pT < 8  */
                                            {0.03, 0.03},  /* 8  < pT < 12 */
                                            {0.03, 0.03},  /* 12 < pT < 24 */
                                            {0.03, 0.03}}; /* 24 < pT < 36 */

// row labels
static const std::vector<std::string> labelsPt = {
  "pT bin 0 Lc",
  "pT bin 1 Lc",
  "pT bin 2 Lc",
  "pT bin 3 Lc",
  "pT bin 4 Lc",
  "pT bin 5 Lc",
  "pT bin 6 Lc",
  "pT bin 7 Lc",
  "pT bin 8 Lc",
  "pT bin 9 Lc"};

// column labels
static const std::vector<std::string> labelsCutVar = {"max pKpi mass Lc", "max piKp mass Lc"};
} // namespace hf_cuts_sigmac_to_p_k_pi

} // namespace o2::analysis

#endif // PWGHF_CORE_SELECTORCUTS_H_
