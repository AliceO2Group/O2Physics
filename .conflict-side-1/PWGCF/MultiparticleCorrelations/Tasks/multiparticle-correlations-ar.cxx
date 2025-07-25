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

/// \file multiparticle-correlations-ar.cxx
/// \brief multiparticle-correlations-ar - Task belonging to Anton Riedel for computing multiparticle correlations
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include <TMath.h>
#include "fairlogger/Logger.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <string>
#include <vector>
#include <array>
#include <numeric>
#include "TComplex.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Expressions.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Define useful constants and enums in a separate name space
namespace MultiParticleCorrelationsARTaskGlobalVariables
{
// for before and after applying cuts
enum BeforeAfterEnum {
  kBEFORE,
  kAFTER,
  kLAST_BA
};
static constexpr std::string_view BeforeAfter[kLAST_BA] = {"before/", "after/"};

// for reconstructed and simulated data
enum RecoSimEnum {
  kRECO,
  kSIM,
  kLAST_RS
};
static constexpr std::string_view RecoSim[kLAST_RS] = {"reco/", "sim/"};

// for applying cuts
enum CutEnum {
  kCUTLOW,    // low cut
  kCUTHIGH,   // high cut
  kCUTOPTION, // option whether to apply a cut at all
  kLAST_CUT
};

// for differential computation of correlators with event variables
enum CorEventDep {
  kINTEGRATED,
  kMULDEP,
  kCENDEP,
  kLAST_CorEventDep
};
const std::string CorEventDepNames[kLAST_CorEventDep] = {
  "[kINTEGRATED]",
  "[kMULDEP]",
  "[kCENDEP]"};
std::vector<std::vector<double>> CorEventDepDefaultBins = {
  {1., 0., 1.},      // kINTEGRATED
  {300., 0., 3000.}, // kMULDEP
  {120., 0., 120.}}; // kCENDEP

// for differential computation of correlators with track variables
enum CorTrackDep {
  kPTDEP,
  kETADEP,
  kLAST_CorTrackDep
};
const std::string CorTrackDepNames[kLAST_CorTrackDep] = {
  "[kPTDEP]",
  "[kETADEP]"};
std::vector<std::vector<double>> CorTrackDepDefaultBins = {
  {VARIABLE_WIDTH, 0.2, 0.34, 0.5, 0.7, 1., 2., 5.}, // kPT
  {VARIABLE_WIDTH, -0.8, -0.4, 0., 0.4, 0.8}};       // kETA

// event variables
enum EventVariable {
  kVX,
  kVY,
  kVZ,
  kVABS,
  kCEN,
  kMULQ,
  kMULW,
  kMULNC,
  kMULTPC,
  kLAST_EventVariable
};
static constexpr std::string_view EventVariableNames[kLAST_EventVariable] = {
  "EventVariable_VertexX",
  "EventVariable_VertexY",
  "EventVariable_VertexZ",
  "EventVariable_VertexAbs",
  "EventVariable_Centrality",
  "EventVariable_MultiplicityQvector",
  "EventVariable_MultiplicityWeights",
  "EventVariable_MultiplicityNumContrib",
  "EventVariable_MultiplicityTPC"};
std::vector<std::vector<double>> EventVariableDefaultBins = {
  {1000, -5., 5.},     // kVX
  {1000., -5., 5.},    // kVY
  {1000., -15., 15.},  // kVZ
  {1000, 0., 15.},     // kVABS
  {120., 0., 120.},    // kCEN
  {3000., 0., 3000.},  // kMULQ
  {3000., 0., 3000.},  // kMULW
  {3000., 0., 3000.},  // kMULNC
  {3000., 0., 3000.}}; // kMULTPC
const std::vector<std::array<float, kLAST_CUT>> EventVariableDefaultCuts = {
  {-2., 2., 1.},    // kVX
  {-2., 2., 1.},    // kVY
  {-10., 10., 1.},  // kVZ
  {1.e-5, 15., 1.}, // kVABS
  {0., 80., 1.},    // kCEN
  {2., 3000., 1.},  // kMULQ
  {2., 3000., 1.},  // kMULW
  {2., 3000., 1.},  // kMULNC
  {2., 3000., 1.}}; // kMULTPC

// track variables
enum TrackVariable {
  kPT,
  kPHI,
  kETA,
  kCHARGE,
  kDCAZ,
  kDCAXY,
  kTPCCLUSTERS,
  kTPCCROSSEDROWS,
  kTPCCHI2,
  kITSCLUSTERS,
  kLAST_TrackVariable
};
static constexpr std::string_view TrackVariableNames[kLAST_TrackVariable] = {
  "TrackVariable_Pt",
  "TrackVariable_Phi",
  "TrackVariable_Eta",
  "TrackVariable_Charge",
  "TrackVariable_DCAZ",
  "TrackVariable_DCAXY",
  "TrackVariable_TPCClusters",
  "TrackVariable_TPCCrossedRows",
  "TrackVariable_TPCChi2",
  "TrackVariable_ITSClusters"};
std::vector<std::vector<double>> TrackVariableDefaultBins = {
  {1000, 0., 10.},           // kPT
  {720., 0, TMath::TwoPi()}, // kPHI
  {1000., -1., 1.},          // kETA
  {5., -2.5, 2.5},           // kCHARGE
  {1000., -5., 5.},          // kDCAZ
  {1000., -5., 5.},          // kDCAXY
  {160., 0., 160.},          // kTPCCLUSTERS
  {160., 0., 160.},          // kTPCCROSSEDROWS
  {1000., 0., 5.},           // kTPCCHI2
  {6., 0., 6.}};             // kITSCLUSTERS
std::vector<std::array<float, kLAST_CUT>> TrackVariableDefaultCuts = {
  {0.2, 5., 1.},           // kPT
  {0., TMath::TwoPi(), 0}, // kPHI
  {-0.8, 0.8, 1.},         // kETA
  {-1.5, 1.5, 1},          // kCHARGE
  {-2.4, 2.4, 1.},         // kDCAZ
  {-3.2, 3.2, 1.},         // kDCAXY
  {80., 161., 1.},         // kTPCCLUSTERS
  {80., 161., 1.},         // kTPCCROSSEDROWS
  {0.4, 5., 1.},           // kTPCCHI2
  {0., 6., 1.}};           // kITSCLUSTERS

// common info string for all configurables
const std::string BinningSuffix = std::string("-Binning");
const std::string CutSuffix = std::string("-Cut");
const std::string CorrelatorHistPrefix = std::string("CorDepBinning_");
const std::vector<std::string> CutInfo = {std::string("Low"), std::string("High"), std::string("Option")};

const int MaxHarmonic = 10;
const int MaxPower = 10;

// function for computing the absolute distance of the primary vertex form the origin
inline float abs(float vx, float vy, float vz)
{
  return std::sqrt(vx * vx + vy * vy * vz * vz);
}

// generice function for checking if the value of a variable passes a cut
inline bool SurviveCut(LabeledArray<float> ConfigValue, float Value)
{
  uint32_t row = 0;
  bool flag = true;
  // check if the cut is configured to be use in the first place
  if (ConfigValue.get(row, kCUTOPTION) > 0.) {
    // check if the value of the variable is not lower than the lower bound and
    // not not larger than the upper bound
    if (!(Value > ConfigValue.get(row, kCUTLOW) && Value < ConfigValue.get(row, kCUTHIGH))) {
      flag = false;
    }
  }
  return flag;
}
// split an integer into its digit
// used for parsing symmetric cumulants as input like
// 23 -> {2,3} <- SC(2,3)
inline std::vector<int> SplitNumber(int n)
{
  std::vector<int> digits;
  while (n >= 10) {
    digits.push_back(n % 10);
    n /= 10;
  }
  digits.push_back(n);
  return digits;
}
}; // namespace MultiParticleCorrelationsARTaskGlobalVariables

// use an alias for our namespace
namespace AR = MultiParticleCorrelationsARTaskGlobalVariables;

struct MultiParticleCorrelationsARTask {

  Configurable<bool> cfgVerbosity = {"verbose", false, "Set to false to silence all info and warning messages"};

  // event control histogram binning and cuts
  ConfigurableAxis cfgEventAxisVX{
    std::string(AR::EventVariableNames[AR::kVX]) + AR::BinningSuffix,
    std::move(AR::EventVariableDefaultBins.at(AR::kVX)),
    ""};
  Configurable<LabeledArray<float>> cfgEventCutVX{
    std::string(AR::EventVariableNames[AR::kVX]) + AR::CutSuffix,
    {AR::EventVariableDefaultCuts.at(AR::kVX).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgEventAxisVY{
    std::string(AR::EventVariableNames[AR::kVY]) + AR::BinningSuffix,
    std::move(AR::EventVariableDefaultBins.at(AR::kVY)),
    ""};
  Configurable<LabeledArray<float>> cfgEventCutVY{
    std::string(AR::EventVariableNames[AR::kVY]) + AR::CutSuffix,
    {AR::EventVariableDefaultCuts.at(AR::kVY).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgEventAxisVZ{
    std::string(AR::EventVariableNames[AR::kVZ]) + AR::BinningSuffix,
    std::move(AR::EventVariableDefaultBins.at(AR::kVZ)),
    ""};
  Configurable<LabeledArray<float>> cfgEventCutVZ{
    std::string(AR::EventVariableNames[AR::kVZ]) + AR::CutSuffix,
    {AR::EventVariableDefaultCuts.at(AR::kVZ).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgEventAxisVABS{
    std::string(AR::EventVariableNames[AR::kVABS]) + AR::BinningSuffix,
    std::move(AR::EventVariableDefaultBins.at(AR::kVABS)),
    ""};
  Configurable<LabeledArray<float>> cfgEventCutVABS{
    std::string(AR::EventVariableNames[AR::kVABS]) + AR::CutSuffix,
    {AR::EventVariableDefaultCuts.at(AR::kVABS).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgEventAxisCEN{
    std::string(AR::EventVariableNames[AR::kCEN]) + AR::BinningSuffix,
    std::move(AR::EventVariableDefaultBins.at(AR::kCEN)),
    ""};
  Configurable<LabeledArray<float>> cfgEventCutCEN{
    std::string(AR::EventVariableNames[AR::kCEN]) + AR::CutSuffix,
    {AR::EventVariableDefaultCuts.at(AR::kCEN).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgEventAxisMULQ{
    std::string(AR::EventVariableNames[AR::kMULQ]) + AR::BinningSuffix,
    std::move(AR::EventVariableDefaultBins.at(AR::kMULQ)),
    ""};
  Configurable<LabeledArray<float>> cfgEventCutMULQ{
    std::string(AR::EventVariableNames[AR::kMULQ]) + AR::CutSuffix,
    {AR::EventVariableDefaultCuts.at(AR::kMULQ).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgEventAxisMULW{
    std::string(AR::EventVariableNames[AR::kMULW]) + AR::BinningSuffix,
    std::move(AR::EventVariableDefaultBins.at(AR::kMULW)),
    ""};
  Configurable<LabeledArray<float>> cfgEventCutMULW{
    std::string(AR::EventVariableNames[AR::kMULW]) + AR::CutSuffix,
    {AR::EventVariableDefaultCuts.at(AR::kMULW).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgEventAxisMULNC{
    std::string(AR::EventVariableNames[AR::kMULNC]) + AR::BinningSuffix,
    std::move(AR::EventVariableDefaultBins.at(AR::kMULNC)),
    ""};
  Configurable<LabeledArray<float>> cfgEventCutMULNC{
    std::string(AR::EventVariableNames[AR::kMULNC]) + AR::CutSuffix,
    {AR::EventVariableDefaultCuts.at(AR::kMULNC).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgEventAxisMULTPC{
    std::string(AR::EventVariableNames[AR::kMULTPC]) + AR::BinningSuffix,
    std::move(AR::EventVariableDefaultBins.at(AR::kMULTPC)),
    ""};
  Configurable<LabeledArray<float>> cfgEventCutMULTPC{
    std::string(AR::EventVariableNames[AR::kMULTPC]) + AR::CutSuffix,
    {AR::EventVariableDefaultCuts.at(AR::kMULTPC).data(), 3, AR::CutInfo},
    ""};
  std::vector<ConfigurableAxis> cfgEventBinning = {
    cfgEventAxisVX,
    cfgEventAxisVY,
    cfgEventAxisVZ,
    cfgEventAxisVABS,
    cfgEventAxisCEN,
    cfgEventAxisMULQ,
    cfgEventAxisMULW,
    cfgEventAxisMULNC,
    cfgEventAxisMULTPC};

  // track control histogram binning and cuts
  ConfigurableAxis cfgTrackAxisPT{
    std::string(AR::TrackVariableNames[AR::kPT]) + AR::BinningSuffix,
    std::move(AR::TrackVariableDefaultBins.at(AR::kPT)),
    ""};
  Configurable<LabeledArray<float>> cfgTrackCutPT{
    std::string(AR::TrackVariableNames[AR::kPT]) + AR::CutSuffix,
    {AR::TrackVariableDefaultCuts.at(AR::kPT).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgTrackAxisPHI{
    std::string(AR::TrackVariableNames[AR::kPHI]) + AR::BinningSuffix,
    std::move(AR::TrackVariableDefaultBins.at(AR::kPHI)),
    ""};
  Configurable<LabeledArray<float>> cfgTrackCutPHI{
    std::string(AR::TrackVariableNames[AR::kPHI]) + AR::CutSuffix,
    {AR::TrackVariableDefaultCuts.at(AR::kPHI).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgTrackAxisETA{
    std::string(AR::TrackVariableNames[AR::kETA]) + AR::BinningSuffix,
    std::move(AR::TrackVariableDefaultBins.at(AR::kETA)),
    ""};
  Configurable<LabeledArray<float>> cfgTrackCutETA{
    std::string(AR::TrackVariableNames[AR::kETA]) + AR::CutSuffix,
    {AR::TrackVariableDefaultCuts.at(AR::kETA).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgTrackAxisCHARGE{
    std::string(AR::TrackVariableNames[AR::kCHARGE]) + AR::BinningSuffix,
    std::move(AR::TrackVariableDefaultBins.at(AR::kCHARGE)),
    ""};
  Configurable<LabeledArray<float>> cfgTrackCutCHARGE{
    std::string(AR::TrackVariableNames[AR::kCHARGE]) + AR::CutSuffix,
    {AR::TrackVariableDefaultCuts.at(AR::kCHARGE).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgTrackAxisDCAZ{
    std::string(AR::TrackVariableNames[AR::kDCAZ]) + AR::BinningSuffix,
    std::move(AR::TrackVariableDefaultBins.at(AR::kDCAZ)),
    ""};
  Configurable<LabeledArray<float>> cfgTrackCutDCAZ{
    std::string(AR::TrackVariableNames[AR::kDCAZ]) + AR::CutSuffix,
    {AR::TrackVariableDefaultCuts.at(AR::kDCAZ).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgTrackAxisDCAXY{
    std::string(AR::TrackVariableNames[AR::kDCAXY]) + AR::BinningSuffix,
    std::move(AR::TrackVariableDefaultBins.at(AR::kDCAXY)),
    ""};
  Configurable<LabeledArray<float>> cfgTrackCutDCAXY{
    std::string(AR::TrackVariableNames[AR::kDCAXY]) + AR::CutSuffix,
    {AR::TrackVariableDefaultCuts.at(AR::kDCAXY).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgTrackAxisTPCCLUSTERS{
    std::string(AR::TrackVariableNames[AR::kTPCCLUSTERS]) + AR::BinningSuffix,
    std::move(AR::TrackVariableDefaultBins.at(AR::kTPCCLUSTERS)),
    ""};
  Configurable<LabeledArray<float>> cfgTrackCutTPCCLUSTERS{
    std::string(AR::TrackVariableNames[AR::kTPCCLUSTERS]) + AR::CutSuffix,
    {AR::TrackVariableDefaultCuts.at(AR::kTPCCLUSTERS).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgTrackAxisTPCCROSSEDROWS{
    std::string(AR::TrackVariableNames[AR::kTPCCROSSEDROWS]) + AR::BinningSuffix,
    std::move(AR::TrackVariableDefaultBins.at(AR::kTPCCROSSEDROWS)),
    ""};
  Configurable<LabeledArray<float>> cfgTrackCutTPCCROSSEDROWS{
    std::string(AR::TrackVariableNames[AR::kTPCCROSSEDROWS]) + AR::CutSuffix,
    {AR::TrackVariableDefaultCuts.at(AR::kTPCCROSSEDROWS).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgTrackAxisTPCCHI2{
    std::string(AR::TrackVariableNames[AR::kTPCCHI2]) + AR::BinningSuffix,
    std::move(AR::TrackVariableDefaultBins.at(AR::kTPCCHI2)),
    ""};
  Configurable<LabeledArray<float>> cfgTrackCutTPCCHI2{
    std::string(AR::TrackVariableNames[AR::kTPCCHI2]) + AR::CutSuffix,
    {AR::TrackVariableDefaultCuts.at(AR::kTPCCHI2).data(), 3, AR::CutInfo},
    ""};
  ConfigurableAxis cfgTrackAxisITSCLUSTERS{
    std::string(AR::TrackVariableNames[AR::kITSCLUSTERS]) + AR::BinningSuffix,
    std::move(AR::TrackVariableDefaultBins.at(AR::kITSCLUSTERS)),
    ""};
  Configurable<LabeledArray<float>> cfgTrackCutITSCLUSTERS{
    std::string(AR::TrackVariableNames[AR::kITSCLUSTERS]) + AR::CutSuffix,
    {AR::TrackVariableDefaultCuts.at(AR::kITSCLUSTERS).data(), 3, AR::CutInfo},
    ""};
  std::vector<ConfigurableAxis> cfgTrackBinning = {
    cfgTrackAxisPT,
    cfgTrackAxisPHI,
    cfgTrackAxisETA,
    cfgTrackAxisCHARGE,
    cfgTrackAxisDCAZ,
    cfgTrackAxisDCAXY,
    cfgTrackAxisTPCCLUSTERS,
    cfgTrackAxisTPCCROSSEDROWS,
    cfgTrackAxisTPCCHI2,
    cfgTrackAxisITSCLUSTERS};

  // configurables for differential analysis of correlators using event variables
  ConfigurableAxis cfgCorrelatorAxisINTEGRATED{
    AR::CorrelatorHistPrefix + std::string(AR::CorEventDepNames[AR::kINTEGRATED]) + AR::BinningSuffix,
    std::move(AR::CorEventDepDefaultBins.at(AR::kINTEGRATED)),
    ""};
  ConfigurableAxis cfgCorrelatorAxisMULDEP{
    AR::CorrelatorHistPrefix + std::string(AR::CorEventDepNames[AR::kMULDEP]) + AR::BinningSuffix,
    std::move(AR::CorEventDepDefaultBins.at(AR::kMULDEP)),
    ""};
  ConfigurableAxis cfgCorrelatorAxisCENDEP{
    AR::CorrelatorHistPrefix + std::string(AR::CorEventDepNames[AR::kCENDEP]) + AR::BinningSuffix,
    std::move(AR::CorEventDepDefaultBins.at(AR::kCENDEP)),
    ""};
  std::vector<ConfigurableAxis> cfgCorEventDep = {
    cfgCorrelatorAxisINTEGRATED,
    cfgCorrelatorAxisMULDEP,
    cfgCorrelatorAxisCENDEP};

  // configurables for differential analysis of correlators using track variables
  ConfigurableAxis cfgCorrelatorAxisPTDEP{
    AR::CorrelatorHistPrefix + std::string(AR::CorTrackDepNames[AR::kPTDEP]) + AR::BinningSuffix,
    std::move(AR::CorTrackDepDefaultBins.at(AR::kPTDEP)),
    ""};
  ConfigurableAxis cfgCorrelatorAxisETADEP{
    AR::CorrelatorHistPrefix + std::string(AR::CorTrackDepNames[AR::kETADEP]) + AR::BinningSuffix,
    std::move(AR::CorTrackDepDefaultBins.at(AR::kETADEP)),
    ""};
  std::vector<ConfigurableAxis> cfgCorTrackDep = {
    cfgCorrelatorAxisPTDEP,
    cfgCorrelatorAxisETADEP};

  // configurable for specifying which symmetric cumulants should be computed
  // the symmetric cumulant is only specified so the appropriate correlators are computed and saved in the output
  // the actual value of the symmetric cumulant is computed in the post processing by parsing the output
  Configurable<std::vector<int>> cfgSC = {"SymmetricCumulants", {23, 24, 34}, "Symmetric Cumulants to be computed"};

  // declare histogram registry
  HistogramRegistry fRegistry{
    "MultiParticleCorrelationsARTask",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    false};

  // declare objects for computing qvectors
  // global object holding the qvectors
  std::array<std::array<TComplex, AR::MaxHarmonic>, AR::MaxPower> fQvectors;
  // global object holding all azimuthal angles and weights, used for computing correlators as a function of event variables
  std::vector<double> fAzimuthalAnglesAll;
  std::vector<double> fWeightsAll;
  // histgram for computing the bins of the event variables
  std::array<TH1F*, AR::kLAST_CorEventDep> fEventDepHists;
  // global object holding all azimuthal angles and weights in bins of the track variables, used for computing correlators as a function of track variables
  std::array<std::vector<std::vector<double>>, AR::kLAST_CorTrackDep> fAzimuthalAnglesTrackDep;
  std::array<std::vector<std::vector<double>>, AR::kLAST_CorTrackDep> fWeightsTrackDep;
  // histgram for computing the bins of the track variables
  std::array<TH1F*, AR::kLAST_CorTrackDep> fTrackDepHists;

  // mapping symmetric cumulants to correlators
  std::map<int, std::vector<std::vector<int>>> fMapScToCor;
  // mapping correlators to index in the list of all correlators
  std::map<std::vector<int>, int> fMapCorToIndex;

  // vector of all correlators
  std::vector<std::vector<int>> fCorrelators;

  // list holding nested lists for each correlator
  // each list contains the correlators as a function of event and track variables
  OutputObj<TList> fOutput{"fOutput", OutputObjHandlingPolicy::AnalysisObject};
  TList* fCorrelatorList;

  void init(InitContext&)
  {

    // add control histograms to registry
    for (int rs = 0; rs < AR::kLAST_RS; rs++) {
      for (int ba = 0; ba < AR::kLAST_BA; ba++) {

        // iterate over event configurables
        for (int e = 0; e < AR::kLAST_EventVariable; e++) {
          fRegistry.add((std::string(AR::RecoSim[rs]) +
                         std::string("EventControl/") +
                         std::string(AR::BeforeAfter[ba]) +
                         std::string(AR::EventVariableNames[e]))
                          .c_str(),
                        "",
                        HistType::kTH1D,
                        {cfgEventBinning.at(e)});
        }
        // iterate over track configurables
        for (int t = 0; t < AR::kLAST_TrackVariable; t++) {
          fRegistry.add((std::string(AR::RecoSim[rs]) +
                         std::string("TrackControl/") +
                         std::string(AR::BeforeAfter[ba]) +
                         std::string(AR::TrackVariableNames[t]))
                          .c_str(),
                        "",
                        HistType::kTH1D,
                        {cfgTrackBinning.at(t)});
        }
      }
    }

    // create histograms for binning the correlators with respect to event variables
    for (int eventDep = 0; eventDep < AR::kLAST_CorEventDep; eventDep++) {
      fEventDepHists[eventDep] = new TH1F(AR::CorEventDepNames[eventDep].c_str(),
                                          AR::CorEventDepNames[eventDep].c_str(),
                                          cfgCorEventDep.at(eventDep).value.at(0),
                                          cfgCorEventDep.at(eventDep).value.at(1),
                                          cfgCorEventDep.at(eventDep).value.at(2));
    }
    std::vector<std::vector<double>> tmp1 = {};
    std::vector<float> tmp2 = {};
    // create histograms for binning the correlators with respect to event variables
    for (int trackDep = 0; trackDep < AR::kLAST_CorTrackDep; trackDep++) {
      // remove the first entry, which is just a 0 indication the variable bin width
      tmp2.insert(tmp2.begin(), ++cfgCorTrackDep.at(trackDep).value.begin(), cfgCorTrackDep.at(trackDep).value.end());
      fTrackDepHists[trackDep] = new TH1F(AR::CorTrackDepNames[trackDep].c_str(),
                                          AR::CorTrackDepNames[trackDep].c_str(),
                                          tmp2.size() - 1,
                                          tmp2.data());

      // push back empty vector, each on corresponding to a different bin of a track variable
      fAzimuthalAnglesTrackDep[trackDep].clear();
      tmp1.clear();
      for (std::size_t j = 0; j < tmp2.size() - 1; j++) {
        tmp1.push_back(std::vector<double>{});
      }
      fAzimuthalAnglesTrackDep[trackDep] = tmp1;
      fWeightsTrackDep[trackDep] = tmp1;

      tmp2.clear();
    }

    // create master list holding the results of all the correlators
    fCorrelatorList = new TList();
    fCorrelatorList->SetOwner(true);
    fOutput.setObject(fCorrelatorList);

    // convert SC, given as input, into a unique list of all needed correlators
    GetCorreltors();
    // book a list for each correlator, containing a TProfile for each track and event variable used as an dependency
    BookCorrelators();
  }

  std::vector<std::vector<int>>
    MapSCToCor(int SC)
  {
    // map symmetric cumulant to the correlators needed for its computation
    // the sc is given as an integer, i.e. 23 -> SC(2,3) -> { {-3,-2,2,3}, {-3,3}, {-2,2} }

    std::vector<int> sc = AR::SplitNumber(SC);
    std::sort(sc.begin(), sc.end());
    std::vector<std::vector<int>> correlators;

    switch (sc.size()) {
      case 2:
        correlators = {
          {-sc.at(0), -sc.at(1), sc.at(1), sc.at(0)}, // <4>_{-l,-k,k,l}
          {-sc.at(0), sc.at(0)},                      // <2>_{-k, k}
          {-sc.at(1), sc.at(1)}                       // <2>_{-l, l}
        };
        break;
      case 3:
        correlators = {
          {-sc.at(0), -sc.at(1), -sc.at(2), sc.at(2), sc.at(1), sc.at(0)}, // <6>_{-k,-l,-n,n,l,k}
          {-sc.at(0), -sc.at(1), sc.at(1), sc.at(0)},                      // <4>_{-k,-l,l,k}
          {-sc.at(0), -sc.at(2), sc.at(2), sc.at(0)},                      // <4>_{-k,-n,n,k}
          {-sc.at(1), -sc.at(2), sc.at(2), sc.at(1)},                      // <4>_{-l,-n,n,l}
          {-sc.at(0), sc.at(0)},                                           // <2>_{-k, k}
          {-sc.at(1), sc.at(1)},                                           // <2>_{-l, l}
          {-sc.at(2), sc.at(2)}};                                          // <2>_{ -n, n }
        break;
      default:
        LOG(fatal) << "Symmetric Cumulants of order " << sc.size() << " are not implemented yet";
    }
    return correlators;
  }

  void GetCorreltors()
  {
    // convert symmetric cumulants, given as input, to a unique list of correlators
    int Index = 0;
    std::vector<std::vector<int>> Correlators;
    for (auto SC : cfgSC.value) {
      Correlators = MapSCToCor(SC);
      fMapScToCor.insert({SC, Correlators});
      for (auto cor : Correlators) {
        if (std::find(fCorrelators.begin(), fCorrelators.end(), cor) !=
            fCorrelators.end()) {
          continue;
        } else {
          fCorrelators.push_back(cor);
          // use a map, so we can later figure out at which list index the correlators resides at
          fMapCorToIndex.insert({cor, Index});
          Index++;
        }
      }
    }
  }

  void BookCorrelators()
  {
    // Book profiles holding correlator as a function of event and track variables

    TList* corList;
    TProfile* profileEventDep[AR::kLAST_CorEventDep];
    TProfile* profileTrackDep[AR::kLAST_CorTrackDep];
    std::string corListName;
    std::string corName;

    // create a uniqe name for each correlator
    for (std::size_t i = 0; i < fCorrelators.size(); i++) {
      corListName = std::string("v_{");
      for (std::size_t j = 0; j < fCorrelators.at(i).size(); j++) {
        if (j != fCorrelators.at(i).size() - 1) {
          corListName += std::to_string(fCorrelators.at(i).at(j));
        } else {
          corListName += std::to_string(fCorrelators.at(i).at(j));
        }
      }

      corListName += std::string("}");
      // create a new list for each correlator
      corList = new TList();
      corList->SetName(corListName.c_str());

      // create profiles for event variables
      for (int eventDep = 0; eventDep < AR::kLAST_CorEventDep; eventDep++) {
        profileEventDep[eventDep] = new TProfile((corListName + AR::CorEventDepNames[eventDep]).c_str(),
                                                 (corListName + AR::CorEventDepNames[eventDep]).c_str(),
                                                 cfgCorEventDep.at(eventDep).value.at(0),
                                                 cfgCorEventDep.at(eventDep).value.at(1),
                                                 cfgCorEventDep.at(eventDep).value.at(2));
        corList->Add(profileEventDep[eventDep]);
      }

      // create profiles for track variables
      std::vector<float> tmp;
      for (int trackDep = 0; trackDep < AR::kLAST_CorTrackDep; trackDep++) {
        // remove the first entry, which is just a 0 indication the variable bin width
        tmp.insert(tmp.begin(), ++cfgCorTrackDep.at(trackDep).value.begin(), cfgCorTrackDep.at(trackDep).value.end());
        profileTrackDep[trackDep] = new TProfile((corListName + AR::CorTrackDepNames[trackDep]).c_str(),
                                                 (corListName + AR::CorTrackDepNames[trackDep]).c_str(),
                                                 tmp.size() - 1,
                                                 tmp.data());
        tmp.clear();
        corList->Add(profileTrackDep[trackDep]);
      }
      // add the list to the master list
      // to access the profiles in list be aware that
      // for the event profiles just use the corresponding enum
      // and for track profiles use the corresponding enum with AR::kLAST_CorEventDep as offset
      fCorrelatorList->Add(corList);
    }
  }

  // function for filling event control histograms
  // exclude MultiplicityQ, i.e. number of tracks in the QVector, and
  // MultiplicityW, i.e. the weighted number of tracks in the QVector, and fill them separably
  template <AR::RecoSimEnum rs, AR::BeforeAfterEnum ba, typename CollisionObject>
  void FillEventControlHist(CollisionObject const& collision, HistogramRegistry& registry)
  {
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kVX]),
                  collision.posX());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kVY]),
                  collision.posY());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kVZ]),
                  collision.posZ());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kVABS]),
                  std::sqrt(std::pow(collision.posX(), 2) + std::pow(collision.posY(), 2) + std::pow(collision.posZ(), 2)));
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kCEN]),
                  collision.centRun2V0M());
    // registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULQ]), collision.size()); fill separately
    // registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULW]), collision.size()); fill separately
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULNC]),
                  collision.numContrib());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULTPC]),
                  collision.multTPC());
  }

  // function for filling event control histograms for Multiplicity{Q,W},
  // since they have to be computed separately
  template <AR::RecoSimEnum rs, AR::BeforeAfterEnum ba>
  void FillEventControlHistMul(HistogramRegistry& registry, double mulQvector, double mulWeights)
  {
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULQ]),
                  mulQvector);
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULW]),
                  mulWeights);
  }

  // function for filling track control histograms
  template <AR::RecoSimEnum rs, AR::BeforeAfterEnum ba, typename TrackObject>
  void FillTrackControlHist(TrackObject const& track, HistogramRegistry& registry)
  {
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kPT]),
                  track.pt());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kPHI]),
                  track.phi());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kETA]),
                  track.eta());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kCHARGE]),
                  track.sign());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kDCAZ]),
                  track.dcaZ());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kDCAXY]),
                  track.dcaXY());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kTPCCLUSTERS]),
                  track.tpcNClsFound());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kTPCCROSSEDROWS]),
                  track.tpcNClsCrossedRows());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kTPCCHI2]),
                  track.tpcChi2NCl());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kITSCLUSTERS]),
                  track.itsNCls());
  }

  // function for checking if collision survives event cuts
  template <typename CollisionObject, typename TrackObject>
  bool SurviveEventCuts(CollisionObject collision, TrackObject tracks)
  {

    // Check if event survives event cuts, where we can get the values for the variables immediately
    if (!(AR::SurviveCut(cfgEventCutVX.value, collision.posX()) &&
          AR::SurviveCut(cfgEventCutVY.value, collision.posY()) &&
          AR::SurviveCut(cfgEventCutVZ.value, collision.posZ()) &&
          AR::SurviveCut(cfgEventCutVABS.value, AR::abs(collision.posX(), collision.posY(), collision.posZ())) &&
          AR::SurviveCut(cfgEventCutCEN.value, collision.centRun2V0M()) &&
          AR::SurviveCut(cfgEventCutMULNC.value, collision.numContrib()) &&
          AR::SurviveCut(cfgEventCutMULTPC.value, collision.multTPC()))) {
      return false;
    }

    int MultiplicityQ = 0.;
    float MultiplicityW = 0.;

    // only compute Multiplicity{Q,W} if the other checks pass, i.e. our flag is still set to true
    for (auto& track : tracks) {
      if (SurviveTrackCuts(track) == true) {
        MultiplicityQ += 1.;
        MultiplicityW += 1.;
      }
    }

    // at last, check if event also passes multiplicity cuts
    return AR::SurviveCut(cfgEventCutMULQ.value, MultiplicityQ) && AR::SurviveCut(cfgEventCutMULW.value, MultiplicityW);
  }

  // function for checking if track survices trach cuts
  template <typename TrackObject>
  bool SurviveTrackCuts(TrackObject track)
  {
    // if all SurviveCut return true, the function will return true
    // if at least one fails, it will return false
    return AR::SurviveCut(cfgTrackCutPT.value, track.pt()) &&
           AR::SurviveCut(cfgTrackCutPHI.value, track.phi()) &&
           AR::SurviveCut(cfgTrackCutETA.value, track.eta()) &&
           AR::SurviveCut(cfgTrackCutCHARGE.value, track.sign()) &&
           AR::SurviveCut(cfgTrackCutDCAZ.value, track.dcaZ()) &&
           AR::SurviveCut(cfgTrackCutDCAXY.value, track.dcaXY()) &&
           AR::SurviveCut(cfgTrackCutTPCCLUSTERS.value, track.tpcNClsFound()) &&
           AR::SurviveCut(cfgTrackCutTPCCROSSEDROWS.value, track.tpcNClsCrossedRows()) &&
           AR::SurviveCut(cfgTrackCutTPCCHI2.value, track.tpcChi2NCl()) &&
           AR::SurviveCut(cfgTrackCutITSCLUSTERS.value, track.itsNCls());
  }

  template <typename TrackObject>
  void FillAzimuthalAngle(TrackObject track)
  {
    double angle = track.phi();
    double weight = GetWeight<TrackObject>(track);

    // all event angles, for computing correlators as a function of event variable
    fAzimuthalAnglesAll.push_back(angle);
    fWeightsAll.push_back(weight);

    // fill angles into bins of a track variable, for computing correlators as a function of track variable
    std::array<double, AR::kLAST_CorTrackDep> TrackDep = {track.pt(), track.eta()};
    int bin;
    for (int trackDep = 0; trackDep < AR::kLAST_CorTrackDep; trackDep++) {
      bin = fTrackDepHists[trackDep]->FindBin(TrackDep[trackDep]) - 1; // in root, bins start at 1
      fAzimuthalAnglesTrackDep[trackDep].at(bin).push_back(angle);
      fWeightsTrackDep[trackDep].at(bin).push_back(weight);
    }
  };

  template <typename TrackObject>
  double GetWeight(TrackObject /*track*/)
  {
    // for efficiency corrections, tbi
    return 1.;
  }

  // Calculate all Q-vectors
  void CalculateQvectors(std::vector<double> AzimuthalAngles, std::vector<double> Weights)
  {
    // Make sure all Q-vectors are initially zero
    for (int h = 0; h < AR::MaxHarmonic; h++) {
      for (int p = 0; p < AR::MaxPower; p++) {
        fQvectors[h][p] = TComplex(0., 0.);
      }
    }
    // Calculate Q-vectors for available angles and weights
    double dPhi = 0.;
    double wPhi = 1.;         // particle weight
    double wPhiToPowerP = 1.; // particle weight raised to power p
    for (std::size_t i = 0; i < AzimuthalAngles.size(); i++) {
      dPhi = AzimuthalAngles.at(i);
      wPhi = Weights.at(i);
      for (int h = 0; h < AR::MaxHarmonic; h++) {
        for (int p = 0; p < AR::MaxPower; p++) {
          wPhiToPowerP = TMath::Power(wPhi, p);
          fQvectors[h][p] += TComplex(wPhiToPowerP * TMath::Cos(h * dPhi),
                                      wPhiToPowerP * TMath::Sin(h * dPhi));
        }
      }
    }
  }

  void FillCorrelators(double Centrality)
  {

    double corr = 0.0;
    double weight = 1.0;

    // compute q vectors for all angles in the event
    CalculateQvectors(fAzimuthalAnglesAll, fWeightsAll);
    int Index;
    for (int eventDep = 0; eventDep < AR::kLAST_CorEventDep; eventDep++) {
      for (auto correlator : fCorrelators) {
        if (fAzimuthalAnglesAll.size() <= correlator.size()) {
          if (cfgVerbosity.value) {
            LOG(warning) << "BEGIN WARNING";
            LOG(warning) << "Not enough tracks in the event to compute the correlator v_{";
            std::for_each(correlator.begin(), correlator.end(), [](const auto& e) { LOG(warning) << e << ","; });
            LOG(warning) << "}!";
            LOG(warning) << "END WARNING";
          }
          continue;
        }
        // get index of the correlator in fCorrelator list
        Index = fMapCorToIndex[correlator];
        // compute the correlator
        ComputeCorrelator(correlator, &corr, &weight);
        // fill the correlator depending on the event variable
        dynamic_cast<TProfile*>(dynamic_cast<TList*>(fCorrelatorList->At(Index))->At(AR::kINTEGRATED))->Fill(0.5, corr, weight);
        dynamic_cast<TProfile*>(dynamic_cast<TList*>(fCorrelatorList->At(Index))->At(AR::kMULDEP))->Fill(fEventDepHists[AR::kMULDEP]->FindBin(fAzimuthalAnglesAll.size()), corr, weight);
        dynamic_cast<TProfile*>(dynamic_cast<TList*>(fCorrelatorList->At(Index))->At(AR::kCENDEP))->Fill(fEventDepHists[AR::kCENDEP]->FindBin(Centrality), corr, weight);
      }
    }

    // loop over all track variables
    for (int trackDep = 0; trackDep < AR::kLAST_CorTrackDep; trackDep++) {
      // loop over all bins of the track variable
      // note that in ROOT the bins of a histogram start at 1
      for (std::size_t bin = 0; bin < fAzimuthalAnglesTrackDep[trackDep].size(); bin++) {
        // compute the qvectors in each bin
        CalculateQvectors(fAzimuthalAnglesTrackDep[trackDep].at(bin),
                          fWeightsTrackDep[trackDep].at(bin));
        corr = 0.;
        weight = 1.;
        // loop over all correlators
        for (auto correlator : fCorrelators) {
          // check if there are enough tracks for computing the correlator
          if (fAzimuthalAnglesTrackDep[trackDep].at(bin).size() <= correlator.size()) {
            if (cfgVerbosity.value) {
              LOG(warning) << "BEGIN WARNING";
              LOG(warning) << "Not enough tracks in track bin to compute the correlator v_{";
              std::for_each(correlator.begin(), correlator.end(), [](const auto& e) { LOG(warning) << e << ","; });
              LOG(warning) << "}!";
              LOG(warning) << "Track variable: " << AR::CorTrackDepNames[trackDep];
              LOG(warning) << "Bin: " << fTrackDepHists[trackDep]->GetBinLowEdge(bin + 1) << " - " << fTrackDepHists[trackDep]->GetBinLowEdge(bin + 2);
              LOG(warning) << "END WARNING";
            }
            continue;
          }
          Index = fMapCorToIndex[correlator];
          // compute the correlator in this bin of the track variable
          ComputeCorrelator(correlator, &corr, &weight);
          // fill it into the corresponding profile
          dynamic_cast<TProfile*>(dynamic_cast<TList*>(fCorrelatorList->At(Index))->At(AR::kLAST_CorEventDep + trackDep))->Fill(fTrackDepHists[trackDep]->GetBinCenter(bin + 1), corr, weight);
        }
      }
    }
  };

  void ComputeCorrelator(std::vector<int> Correlator, double* value, double* weight)
  {
    // compute a correlator
    // write back its value and weight into the passed pointers
    double Value = 0., Weight = 1.;
    switch (static_cast<int>(Correlator.size())) {
      case 2:
        Value = Two(Correlator.at(0), Correlator.at(1)).Re();
        Weight = Two(0, 0).Re();
        break;
      case 3:
        Value = Three(Correlator.at(0), Correlator.at(1), Correlator.at(2)).Re();
        Weight = Three(0, 0, 0).Re();
        break;
      case 4:
        Value = Four(Correlator.at(0), Correlator.at(1), Correlator.at(2), Correlator.at(3)).Re();
        Weight = Four(0, 0, 0, 0).Re();
        break;
      case 5:
        Value = Five(Correlator.at(0), Correlator.at(1), Correlator.at(2), Correlator.at(3), Correlator.at(4)).Re();
        Weight = Five(0, 0, 0, 0, 0).Re();
        break;
      case 6:
        Value = Six(Correlator.at(0), Correlator.at(1), Correlator.at(2), Correlator.at(3), Correlator.at(4), Correlator.at(5)).Re();
        Weight = Six(0, 0, 0, 0, 0, 0).Re();
        break;
      default:
        Value = Recursion(Correlator.size(), Correlator.data()).Re();
        Weight =
          Recursion(Correlator.size(), std::vector<int>(Correlator.size(), 0).data()).Re();
    }
    // correlators are not normalized yet
    Value /= Weight;
    // write back the values
    *value = Value;
    *weight = Weight;
  }

  TComplex Q(int n, int p)
  {
    // return Qvector from fQvectors array
    if (n > AR::MaxHarmonic || p > AR::MaxPower) {
      LOG(fatal) << "Harmonic " << n << ">" << AR::MaxHarmonic << "(MaxHarmonic) or the power " << p << ">" << AR::MaxPower << "(MaxPower)";
    }
    if (n >= 0) {
      return fQvectors[n][p];
    }
    return TComplex::Conjugate(fQvectors[-n][p]);
  }

  TComplex Two(int n1, int n2)
  {
    // Generic two-particle correlation <exp[i(n1*phi1+n2*phi2)]>.
    TComplex two = Q(n1, 1) * Q(n2, 1) - Q(n1 + n2, 2);
    return two;
  }

  TComplex Three(int n1, int n2, int n3)
  {
    // Generic three-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3)]>.
    TComplex three = Q(n1, 1) * Q(n2, 1) * Q(n3, 1) - Q(n1 + n2, 2) * Q(n3, 1) -
                     Q(n2, 1) * Q(n1 + n3, 2) - Q(n1, 1) * Q(n2 + n3, 2) +
                     2. * Q(n1 + n2 + n3, 3);
    return three;
  }

  TComplex Four(int n1, int n2, int n3, int n4)
  {
    // Generic four-particle correlation
    // <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4)]>.
    TComplex four =
      Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) -
      Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) -
      Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) -
      Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) + 2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) -
      Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) + Q(n2 + n3, 2) * Q(n1 + n4, 2) -
      Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) + Q(n1 + n3, 2) * Q(n2 + n4, 2) +
      2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) - Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) +
      Q(n1 + n2, 2) * Q(n3 + n4, 2) + 2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) +
      2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) - 6. * Q(n1 + n2 + n3 + n4, 4);

    return four;
  }

  TComplex Five(int n1, int n2, int n3, int n4, int n5)
  {
    // Generic five-particle correlation
    // <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5)]>.
    TComplex five = Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) -
                    Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) -
                    Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) * Q(n5, 1) -
                    Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) * Q(n5, 1) +
                    2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) * Q(n5, 1) -
                    Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) * Q(n5, 1) +
                    Q(n2 + n3, 2) * Q(n1 + n4, 2) * Q(n5, 1) -
                    Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) * Q(n5, 1) +
                    Q(n1 + n3, 2) * Q(n2 + n4, 2) * Q(n5, 1) +
                    2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) * Q(n5, 1) -
                    Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) * Q(n5, 1) +
                    Q(n1 + n2, 2) * Q(n3 + n4, 2) * Q(n5, 1) +
                    2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) * Q(n5, 1) +
                    2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) * Q(n5, 1) -
                    6. * Q(n1 + n2 + n3 + n4, 4) * Q(n5, 1) -
                    Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n1 + n5, 2) +
                    Q(n2 + n3, 2) * Q(n4, 1) * Q(n1 + n5, 2) +
                    Q(n3, 1) * Q(n2 + n4, 2) * Q(n1 + n5, 2) +
                    Q(n2, 1) * Q(n3 + n4, 2) * Q(n1 + n5, 2) -
                    2. * Q(n2 + n3 + n4, 3) * Q(n1 + n5, 2) -
                    Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n2 + n5, 2) +
                    Q(n1 + n3, 2) * Q(n4, 1) * Q(n2 + n5, 2) +
                    Q(n3, 1) * Q(n1 + n4, 2) * Q(n2 + n5, 2) +
                    Q(n1, 1) * Q(n3 + n4, 2) * Q(n2 + n5, 2) -
                    2. * Q(n1 + n3 + n4, 3) * Q(n2 + n5, 2) +
                    2. * Q(n3, 1) * Q(n4, 1) * Q(n1 + n2 + n5, 3) -
                    2. * Q(n3 + n4, 2) * Q(n1 + n2 + n5, 3) -
                    Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n3 + n5, 2) +
                    Q(n1 + n2, 2) * Q(n4, 1) * Q(n3 + n5, 2) +
                    Q(n2, 1) * Q(n1 + n4, 2) * Q(n3 + n5, 2) +
                    Q(n1, 1) * Q(n2 + n4, 2) * Q(n3 + n5, 2) -
                    2. * Q(n1 + n2 + n4, 3) * Q(n3 + n5, 2) +
                    2. * Q(n2, 1) * Q(n4, 1) * Q(n1 + n3 + n5, 3) -
                    2. * Q(n2 + n4, 2) * Q(n1 + n3 + n5, 3) +
                    2. * Q(n1, 1) * Q(n4, 1) * Q(n2 + n3 + n5, 3) -
                    2. * Q(n1 + n4, 2) * Q(n2 + n3 + n5, 3) -
                    6. * Q(n4, 1) * Q(n1 + n2 + n3 + n5, 4) -
                    Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4 + n5, 2) +
                    Q(n1 + n2, 2) * Q(n3, 1) * Q(n4 + n5, 2) +
                    Q(n2, 1) * Q(n1 + n3, 2) * Q(n4 + n5, 2) +
                    Q(n1, 1) * Q(n2 + n3, 2) * Q(n4 + n5, 2) -
                    2. * Q(n1 + n2 + n3, 3) * Q(n4 + n5, 2) +
                    2. * Q(n2, 1) * Q(n3, 1) * Q(n1 + n4 + n5, 3) -
                    2. * Q(n2 + n3, 2) * Q(n1 + n4 + n5, 3) +
                    2. * Q(n1, 1) * Q(n3, 1) * Q(n2 + n4 + n5, 3) -
                    2. * Q(n1 + n3, 2) * Q(n2 + n4 + n5, 3) -
                    6. * Q(n3, 1) * Q(n1 + n2 + n4 + n5, 4) +
                    2. * Q(n1, 1) * Q(n2, 1) * Q(n3 + n4 + n5, 3) -
                    2. * Q(n1 + n2, 2) * Q(n3 + n4 + n5, 3) -
                    6. * Q(n2, 1) * Q(n1 + n3 + n4 + n5, 4) -
                    6. * Q(n1, 1) * Q(n2 + n3 + n4 + n5, 4) +
                    24. * Q(n1 + n2 + n3 + n4 + n5, 5);
    return five;
  }

  TComplex Six(int n1, int n2, int n3, int n4, int n5, int n6)
  {
    // Generic six-particle correlation
    // <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6)]>.
    TComplex six =
      Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) -
      Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) -
      Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) -
      Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) +
      2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) -
      Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n6, 1) +
      Q(n2 + n3, 2) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n6, 1) -
      Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n6, 1) +
      Q(n1 + n3, 2) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n6, 1) +
      2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) * Q(n5, 1) * Q(n6, 1) -
      Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n6, 1) +
      Q(n1 + n2, 2) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n6, 1) +
      2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) * Q(n5, 1) * Q(n6, 1) +
      2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) * Q(n5, 1) * Q(n6, 1) -
      6. * Q(n1 + n2 + n3 + n4, 4) * Q(n5, 1) * Q(n6, 1) -
      Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n6, 1) +
      Q(n2 + n3, 2) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n6, 1) +
      Q(n3, 1) * Q(n2 + n4, 2) * Q(n1 + n5, 2) * Q(n6, 1) +
      Q(n2, 1) * Q(n3 + n4, 2) * Q(n1 + n5, 2) * Q(n6, 1) -
      2. * Q(n2 + n3 + n4, 3) * Q(n1 + n5, 2) * Q(n6, 1) -
      Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n6, 1) +
      Q(n1 + n3, 2) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n6, 1) +
      Q(n3, 1) * Q(n1 + n4, 2) * Q(n2 + n5, 2) * Q(n6, 1) +
      Q(n1, 1) * Q(n3 + n4, 2) * Q(n2 + n5, 2) * Q(n6, 1) -
      2. * Q(n1 + n3 + n4, 3) * Q(n2 + n5, 2) * Q(n6, 1) +
      2. * Q(n3, 1) * Q(n4, 1) * Q(n1 + n2 + n5, 3) * Q(n6, 1) -
      2. * Q(n3 + n4, 2) * Q(n1 + n2 + n5, 3) * Q(n6, 1) -
      Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n6, 1) +
      Q(n1 + n2, 2) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n6, 1) +
      Q(n2, 1) * Q(n1 + n4, 2) * Q(n3 + n5, 2) * Q(n6, 1) +
      Q(n1, 1) * Q(n2 + n4, 2) * Q(n3 + n5, 2) * Q(n6, 1) -
      2. * Q(n1 + n2 + n4, 3) * Q(n3 + n5, 2) * Q(n6, 1) +
      2. * Q(n2, 1) * Q(n4, 1) * Q(n1 + n3 + n5, 3) * Q(n6, 1) -
      2. * Q(n2 + n4, 2) * Q(n1 + n3 + n5, 3) * Q(n6, 1) +
      2. * Q(n1, 1) * Q(n4, 1) * Q(n2 + n3 + n5, 3) * Q(n6, 1) -
      2. * Q(n1 + n4, 2) * Q(n2 + n3 + n5, 3) * Q(n6, 1) -
      6. * Q(n4, 1) * Q(n1 + n2 + n3 + n5, 4) * Q(n6, 1) -
      Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n6, 1) +
      Q(n1 + n2, 2) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n6, 1) +
      Q(n2, 1) * Q(n1 + n3, 2) * Q(n4 + n5, 2) * Q(n6, 1) +
      Q(n1, 1) * Q(n2 + n3, 2) * Q(n4 + n5, 2) * Q(n6, 1) -
      2. * Q(n1 + n2 + n3, 3) * Q(n4 + n5, 2) * Q(n6, 1) +
      2. * Q(n2, 1) * Q(n3, 1) * Q(n1 + n4 + n5, 3) * Q(n6, 1) -
      2. * Q(n2 + n3, 2) * Q(n1 + n4 + n5, 3) * Q(n6, 1) +
      2. * Q(n1, 1) * Q(n3, 1) * Q(n2 + n4 + n5, 3) * Q(n6, 1) -
      2. * Q(n1 + n3, 2) * Q(n2 + n4 + n5, 3) * Q(n6, 1) -
      6. * Q(n3, 1) * Q(n1 + n2 + n4 + n5, 4) * Q(n6, 1) +
      2. * Q(n1, 1) * Q(n2, 1) * Q(n3 + n4 + n5, 3) * Q(n6, 1) -
      2. * Q(n1 + n2, 2) * Q(n3 + n4 + n5, 3) * Q(n6, 1) -
      6. * Q(n2, 1) * Q(n1 + n3 + n4 + n5, 4) * Q(n6, 1) -
      6. * Q(n1, 1) * Q(n2 + n3 + n4 + n5, 4) * Q(n6, 1) +
      24. * Q(n1 + n2 + n3 + n4 + n5, 5) * Q(n6, 1) -
      Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n6, 2) +
      Q(n2 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n6, 2) +
      Q(n3, 1) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n1 + n6, 2) +
      Q(n2, 1) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n1 + n6, 2) -
      2. * Q(n2 + n3 + n4, 3) * Q(n5, 1) * Q(n1 + n6, 2) +
      Q(n3, 1) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n1 + n6, 2) -
      Q(n3 + n4, 2) * Q(n2 + n5, 2) * Q(n1 + n6, 2) +
      Q(n2, 1) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n1 + n6, 2) -
      Q(n2 + n4, 2) * Q(n3 + n5, 2) * Q(n1 + n6, 2) -
      2. * Q(n4, 1) * Q(n2 + n3 + n5, 3) * Q(n1 + n6, 2) +
      Q(n2, 1) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n1 + n6, 2) -
      Q(n2 + n3, 2) * Q(n4 + n5, 2) * Q(n1 + n6, 2) -
      2. * Q(n3, 1) * Q(n2 + n4 + n5, 3) * Q(n1 + n6, 2) -
      2. * Q(n2, 1) * Q(n3 + n4 + n5, 3) * Q(n1 + n6, 2) +
      6. * Q(n2 + n3 + n4 + n5, 4) * Q(n1 + n6, 2) -
      Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n2 + n6, 2) +
      Q(n1 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n2 + n6, 2) +
      Q(n3, 1) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n2 + n6, 2) +
      Q(n1, 1) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n2 + n6, 2) -
      2. * Q(n1 + n3 + n4, 3) * Q(n5, 1) * Q(n2 + n6, 2) +
      Q(n3, 1) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n2 + n6, 2) -
      Q(n3 + n4, 2) * Q(n1 + n5, 2) * Q(n2 + n6, 2) +
      Q(n1, 1) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n2 + n6, 2) -
      Q(n1 + n4, 2) * Q(n3 + n5, 2) * Q(n2 + n6, 2) -
      2. * Q(n4, 1) * Q(n1 + n3 + n5, 3) * Q(n2 + n6, 2) +
      Q(n1, 1) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n2 + n6, 2) -
      Q(n1 + n3, 2) * Q(n4 + n5, 2) * Q(n2 + n6, 2) -
      2. * Q(n3, 1) * Q(n1 + n4 + n5, 3) * Q(n2 + n6, 2) -
      2. * Q(n1, 1) * Q(n3 + n4 + n5, 3) * Q(n2 + n6, 2) +
      6. * Q(n1 + n3 + n4 + n5, 4) * Q(n2 + n6, 2) +
      2. * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n2 + n6, 3) -
      2. * Q(n3 + n4, 2) * Q(n5, 1) * Q(n1 + n2 + n6, 3) -
      2. * Q(n4, 1) * Q(n3 + n5, 2) * Q(n1 + n2 + n6, 3) -
      2. * Q(n3, 1) * Q(n4 + n5, 2) * Q(n1 + n2 + n6, 3) +
      4. * Q(n3 + n4 + n5, 3) * Q(n1 + n2 + n6, 3) -
      Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n5, 1) * Q(n3 + n6, 2) +
      Q(n1 + n2, 2) * Q(n4, 1) * Q(n5, 1) * Q(n3 + n6, 2) +
      Q(n2, 1) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n3 + n6, 2) +
      Q(n1, 1) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n3 + n6, 2) -
      2. * Q(n1 + n2 + n4, 3) * Q(n5, 1) * Q(n3 + n6, 2) +
      Q(n2, 1) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n3 + n6, 2) -
      Q(n2 + n4, 2) * Q(n1 + n5, 2) * Q(n3 + n6, 2) +
      Q(n1, 1) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n3 + n6, 2) -
      Q(n1 + n4, 2) * Q(n2 + n5, 2) * Q(n3 + n6, 2) -
      2. * Q(n4, 1) * Q(n1 + n2 + n5, 3) * Q(n3 + n6, 2) +
      Q(n1, 1) * Q(n2, 1) * Q(n4 + n5, 2) * Q(n3 + n6, 2) -
      Q(n1 + n2, 2) * Q(n4 + n5, 2) * Q(n3 + n6, 2) -
      2. * Q(n2, 1) * Q(n1 + n4 + n5, 3) * Q(n3 + n6, 2) -
      2. * Q(n1, 1) * Q(n2 + n4 + n5, 3) * Q(n3 + n6, 2) +
      6. * Q(n1 + n2 + n4 + n5, 4) * Q(n3 + n6, 2) +
      2. * Q(n2, 1) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n3 + n6, 3) -
      2. * Q(n2 + n4, 2) * Q(n5, 1) * Q(n1 + n3 + n6, 3) -
      2. * Q(n4, 1) * Q(n2 + n5, 2) * Q(n1 + n3 + n6, 3) -
      2. * Q(n2, 1) * Q(n4 + n5, 2) * Q(n1 + n3 + n6, 3) +
      4. * Q(n2 + n4 + n5, 3) * Q(n1 + n3 + n6, 3) +
      2. * Q(n1, 1) * Q(n4, 1) * Q(n5, 1) * Q(n2 + n3 + n6, 3) -
      2. * Q(n1 + n4, 2) * Q(n5, 1) * Q(n2 + n3 + n6, 3) -
      2. * Q(n4, 1) * Q(n1 + n5, 2) * Q(n2 + n3 + n6, 3) -
      2. * Q(n1, 1) * Q(n4 + n5, 2) * Q(n2 + n3 + n6, 3) +
      4. * Q(n1 + n4 + n5, 3) * Q(n2 + n3 + n6, 3) -
      6. * Q(n4, 1) * Q(n5, 1) * Q(n1 + n2 + n3 + n6, 4) +
      6. * Q(n4 + n5, 2) * Q(n1 + n2 + n3 + n6, 4) -
      Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n5, 1) * Q(n4 + n6, 2) +
      Q(n1 + n2, 2) * Q(n3, 1) * Q(n5, 1) * Q(n4 + n6, 2) +
      Q(n2, 1) * Q(n1 + n3, 2) * Q(n5, 1) * Q(n4 + n6, 2) +
      Q(n1, 1) * Q(n2 + n3, 2) * Q(n5, 1) * Q(n4 + n6, 2) -
      2. * Q(n1 + n2 + n3, 3) * Q(n5, 1) * Q(n4 + n6, 2) +
      Q(n2, 1) * Q(n3, 1) * Q(n1 + n5, 2) * Q(n4 + n6, 2) -
      Q(n2 + n3, 2) * Q(n1 + n5, 2) * Q(n4 + n6, 2) +
      Q(n1, 1) * Q(n3, 1) * Q(n2 + n5, 2) * Q(n4 + n6, 2) -
      Q(n1 + n3, 2) * Q(n2 + n5, 2) * Q(n4 + n6, 2) -
      2. * Q(n3, 1) * Q(n1 + n2 + n5, 3) * Q(n4 + n6, 2) +
      Q(n1, 1) * Q(n2, 1) * Q(n3 + n5, 2) * Q(n4 + n6, 2) -
      Q(n1 + n2, 2) * Q(n3 + n5, 2) * Q(n4 + n6, 2) -
      2. * Q(n2, 1) * Q(n1 + n3 + n5, 3) * Q(n4 + n6, 2) -
      2. * Q(n1, 1) * Q(n2 + n3 + n5, 3) * Q(n4 + n6, 2) +
      6. * Q(n1 + n2 + n3 + n5, 4) * Q(n4 + n6, 2) +
      2. * Q(n2, 1) * Q(n3, 1) * Q(n5, 1) * Q(n1 + n4 + n6, 3) -
      2. * Q(n2 + n3, 2) * Q(n5, 1) * Q(n1 + n4 + n6, 3) -
      2. * Q(n3, 1) * Q(n2 + n5, 2) * Q(n1 + n4 + n6, 3) -
      2. * Q(n2, 1) * Q(n3 + n5, 2) * Q(n1 + n4 + n6, 3) +
      4. * Q(n2 + n3 + n5, 3) * Q(n1 + n4 + n6, 3) +
      2. * Q(n1, 1) * Q(n3, 1) * Q(n5, 1) * Q(n2 + n4 + n6, 3) -
      2. * Q(n1 + n3, 2) * Q(n5, 1) * Q(n2 + n4 + n6, 3) -
      2. * Q(n3, 1) * Q(n1 + n5, 2) * Q(n2 + n4 + n6, 3) -
      2. * Q(n1, 1) * Q(n3 + n5, 2) * Q(n2 + n4 + n6, 3) +
      4. * Q(n1 + n3 + n5, 3) * Q(n2 + n4 + n6, 3) -
      6. * Q(n3, 1) * Q(n5, 1) * Q(n1 + n2 + n4 + n6, 4) +
      6. * Q(n3 + n5, 2) * Q(n1 + n2 + n4 + n6, 4) +
      2. * Q(n1, 1) * Q(n2, 1) * Q(n5, 1) * Q(n3 + n4 + n6, 3) -
      2. * Q(n1 + n2, 2) * Q(n5, 1) * Q(n3 + n4 + n6, 3) -
      2. * Q(n2, 1) * Q(n1 + n5, 2) * Q(n3 + n4 + n6, 3) -
      2. * Q(n1, 1) * Q(n2 + n5, 2) * Q(n3 + n4 + n6, 3) +
      4. * Q(n1 + n2 + n5, 3) * Q(n3 + n4 + n6, 3) -
      6. * Q(n2, 1) * Q(n5, 1) * Q(n1 + n3 + n4 + n6, 4) +
      6. * Q(n2 + n5, 2) * Q(n1 + n3 + n4 + n6, 4) -
      6. * Q(n1, 1) * Q(n5, 1) * Q(n2 + n3 + n4 + n6, 4) +
      6. * Q(n1 + n5, 2) * Q(n2 + n3 + n4 + n6, 4) +
      24. * Q(n5, 1) * Q(n1 + n2 + n3 + n4 + n6, 5) -
      Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5 + n6, 2) +
      Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) * Q(n5 + n6, 2) +
      Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) * Q(n5 + n6, 2) +
      Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) * Q(n5 + n6, 2) -
      2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) * Q(n5 + n6, 2) +
      Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) * Q(n5 + n6, 2) -
      Q(n2 + n3, 2) * Q(n1 + n4, 2) * Q(n5 + n6, 2) +
      Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) * Q(n5 + n6, 2) -
      Q(n1 + n3, 2) * Q(n2 + n4, 2) * Q(n5 + n6, 2) -
      2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) * Q(n5 + n6, 2) +
      Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) * Q(n5 + n6, 2) -
      Q(n1 + n2, 2) * Q(n3 + n4, 2) * Q(n5 + n6, 2) -
      2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) * Q(n5 + n6, 2) -
      2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) * Q(n5 + n6, 2) +
      6. * Q(n1 + n2 + n3 + n4, 4) * Q(n5 + n6, 2) +
      2. * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n1 + n5 + n6, 3) -
      2. * Q(n2 + n3, 2) * Q(n4, 1) * Q(n1 + n5 + n6, 3) -
      2. * Q(n3, 1) * Q(n2 + n4, 2) * Q(n1 + n5 + n6, 3) -
      2. * Q(n2, 1) * Q(n3 + n4, 2) * Q(n1 + n5 + n6, 3) +
      4. * Q(n2 + n3 + n4, 3) * Q(n1 + n5 + n6, 3) +
      2. * Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n2 + n5 + n6, 3) -
      2. * Q(n1 + n3, 2) * Q(n4, 1) * Q(n2 + n5 + n6, 3) -
      2. * Q(n3, 1) * Q(n1 + n4, 2) * Q(n2 + n5 + n6, 3) -
      2. * Q(n1, 1) * Q(n3 + n4, 2) * Q(n2 + n5 + n6, 3) +
      4. * Q(n1 + n3 + n4, 3) * Q(n2 + n5 + n6, 3) -
      6. * Q(n3, 1) * Q(n4, 1) * Q(n1 + n2 + n5 + n6, 4) +
      6. * Q(n3 + n4, 2) * Q(n1 + n2 + n5 + n6, 4) +
      2. * Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n3 + n5 + n6, 3) -
      2. * Q(n1 + n2, 2) * Q(n4, 1) * Q(n3 + n5 + n6, 3) -
      2. * Q(n2, 1) * Q(n1 + n4, 2) * Q(n3 + n5 + n6, 3) -
      2. * Q(n1, 1) * Q(n2 + n4, 2) * Q(n3 + n5 + n6, 3) +
      4. * Q(n1 + n2 + n4, 3) * Q(n3 + n5 + n6, 3) -
      6. * Q(n2, 1) * Q(n4, 1) * Q(n1 + n3 + n5 + n6, 4) +
      6. * Q(n2 + n4, 2) * Q(n1 + n3 + n5 + n6, 4) -
      6. * Q(n1, 1) * Q(n4, 1) * Q(n2 + n3 + n5 + n6, 4) +
      6. * Q(n1 + n4, 2) * Q(n2 + n3 + n5 + n6, 4) +
      24. * Q(n4, 1) * Q(n1 + n2 + n3 + n5 + n6, 5) +
      2. * Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4 + n5 + n6, 3) -
      2. * Q(n1 + n2, 2) * Q(n3, 1) * Q(n4 + n5 + n6, 3) -
      2. * Q(n2, 1) * Q(n1 + n3, 2) * Q(n4 + n5 + n6, 3) -
      2. * Q(n1, 1) * Q(n2 + n3, 2) * Q(n4 + n5 + n6, 3) +
      4. * Q(n1 + n2 + n3, 3) * Q(n4 + n5 + n6, 3) -
      6. * Q(n2, 1) * Q(n3, 1) * Q(n1 + n4 + n5 + n6, 4) +
      6. * Q(n2 + n3, 2) * Q(n1 + n4 + n5 + n6, 4) -
      6. * Q(n1, 1) * Q(n3, 1) * Q(n2 + n4 + n5 + n6, 4) +
      6. * Q(n1 + n3, 2) * Q(n2 + n4 + n5 + n6, 4) +
      24. * Q(n3, 1) * Q(n1 + n2 + n4 + n5 + n6, 5) -
      6. * Q(n1, 1) * Q(n2, 1) * Q(n3 + n4 + n5 + n6, 4) +
      6. * Q(n1 + n2, 2) * Q(n3 + n4 + n5 + n6, 4) +
      24. * Q(n2, 1) * Q(n1 + n3 + n4 + n5 + n6, 5) +
      24. * Q(n1, 1) * Q(n2 + n3 + n4 + n5 + n6, 5) -
      120. * Q(n1 + n2 + n3 + n4 + n5 + n6, 6);
    return six;
  }

  TComplex Recursion(int n, int* harmonic, int mult = 1, int skip = 0)
  {
    // Calculate multi-particle correlators by using recursion (an improved
    // faster version) originally developed by Kristjan Gulbrandsen
    // (gulbrand@nbi.dk).
    int nm1 = n - 1;
    TComplex c(Q(harmonic[nm1], mult));
    if (nm1 == 0)
      return c;
    c *= Recursion(nm1, harmonic);
    if (nm1 == skip)
      return c;
    int multp1 = mult + 1;
    int nm2 = n - 2;
    int counter1 = 0;
    int hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
    int counter2 = n - 3;
    while (counter2 >= skip) {
      harmonic[nm2] = harmonic[counter1];
      harmonic[counter1] = hhold;
      ++counter1;
      hhold = harmonic[counter1];
      harmonic[counter1] = harmonic[nm2];
      harmonic[nm2] = hhold + harmonic[nm1];
      c2 += Recursion(nm1, harmonic, multp1, counter2);
      --counter2;
    }
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    if (mult == 1)
      return c - c2;
    return c - double(mult) * c2;
  }

  using CollisionsInstance = soa::Join<aod::Collisions, aod::CentRun2V0Ms, aod::Mults>;
  using TracksInstance = soa::Join<aod::Tracks, aod::TracksDCA, aod::TracksExtra>;

  using CollisionsInstanceIterator = CollisionsInstance::iterator;
  // using TracksInstanceIterator = TracksInstance::iterator;

  void process(CollisionsInstanceIterator const& collision, TracksInstance const& tracks)
  {
    // clear angles and weights
    fAzimuthalAnglesAll.clear();
    fWeightsAll.clear();
    for (int trackDep = 0; trackDep < AR::kLAST_CorTrackDep; trackDep++) {
      for (std::size_t bin = 0; bin < fAzimuthalAnglesTrackDep[trackDep].size(); bin++) {
        fAzimuthalAnglesTrackDep[trackDep].at(bin).clear();
        fWeightsTrackDep[trackDep].at(bin).clear();
      }
    }

    if (cfgVerbosity.value) {
      LOG(info) << "Process event: " << collision.index();
    }

    FillEventControlHist<AR::kRECO, AR::kBEFORE, CollisionsInstanceIterator>(collision, fRegistry);
    FillEventControlHistMul<AR::kRECO, AR::kBEFORE>(fRegistry, collision.size(), collision.size());

    if (!SurviveEventCuts(collision, tracks)) {
      if (cfgVerbosity.value) {
        LOG(info) << "Event was CUT";
      }
      return;
    }

    FillEventControlHist<AR::kRECO, AR::kAFTER, CollisionsInstanceIterator>(collision, fRegistry);

    // loop over all tracks in the event
    for (auto const& track : tracks) {
      // fill track control histograms before track cut
      FillTrackControlHist<AR::kRECO, AR::kBEFORE, TracksInstance::iterator>(track, fRegistry);

      if (!SurviveTrackCuts(track)) {
        continue;
      }
      // fill track control histograms after cut
      FillTrackControlHist<AR::kRECO, AR::kAFTER, TracksInstance::iterator>(track, fRegistry);

      // fill angles and weights into vectors for processing
      FillAzimuthalAngle<TracksInstance::iterator>(track);
    }

    FillEventControlHistMul<AR::kRECO, AR::kAFTER>(fRegistry, fAzimuthalAnglesAll.size(), std::accumulate(fWeightsAll.begin(), fWeightsAll.end(), 0.));

    // compute correlators
    FillCorrelators(collision.centRun2V0M());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiParticleCorrelationsARTask>(cfgc)};
}
