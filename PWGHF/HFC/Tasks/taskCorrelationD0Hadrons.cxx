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

/// \file taskCorrelationD0Hadrons.cxx
/// \brief D0-Hadron correlator task - data-like, MC-reco and MC-kine analyses.
/// \note Extended from taskCorrelationDDbar
///
/// \author Samrangy Sadhu <samrangy.sadhu@cern.ch>, INFN Bari
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>, IIT Indore

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_correlation_d0_hadron;

namespace o2::aod
{
using DHadronPairFull = soa::Join<aod::DHadronPair, aod::DHadronRecoInfo>;
} // namespace o2::aod

// string definitions, used for histogram axis labels
const TString stringPtD = "#it{p}_{T}^{D} (GeV/#it{c});";
const TString stringPtHadron = "#it{p}_{T}^{Hadron} (GeV/#it{c});";
const TString stringDeltaEta = "#it{#eta}^{Hadron}-#it{#eta}^{D};";
const TString stringDeltaPhi = "#it{#varphi}^{Hadron}-#it{#varphi}^{D} (rad);";
const TString stringDHadron = "D,Hadron candidates ";
const TString stringSignal = "signal region;";
const TString stringSideband = "sidebands;";
const TString stringMcParticles = "MC gen - D,Hadron particles;";
const TString stringMcReco = "MC reco - D,Hadron candidates ";

// histogram axes definition
AxisSpec axisDeltaEta = {100, -2., 2., ""};
AxisSpec axisDeltaPhi = {64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf, ""};
AxisSpec axisPtD = {10, 0., 10., ""};
AxisSpec axisPtHadron = {11, 0., 11., ""};
AxisSpec axisPoolBin = {9, 0., 9., ""};

// definition of vectors for standard ptbin and invariant mass configurables
const int nPtBinsCorrelations = 8;
const double pTBinsCorrelations[nPtBinsCorrelations + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
auto vecPtBinsCorrelations = std::vector<double>{pTBinsCorrelations, pTBinsCorrelations + nPtBinsCorrelations + 1};
const double signalRegionLeftDefault[nPtBinsCorrelations] = {1.810, 1.810, 1.810, 1.810, 1.810, 1.810, 1.810, 1.810};
const double signalRegionRightDefault[nPtBinsCorrelations] = {1.922, 1.922, 1.922, 1.922, 1.922, 1.922, 1.922, 1.922};
const double sidebandLeftInnerDefault[nPtBinsCorrelations] = {1.754, 1.754, 1.754, 1.754, 1.754, 1.754, 1.754, 1.754};
const double sidebandLeftOuterDefault[nPtBinsCorrelations] = {1.642, 1.642, 1.642, 1.642, 1.642, 1.642, 1.642, 1.642};
const double sidebandRightInnerDefault[nPtBinsCorrelations] = {1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978};
const double sidebandRightOuterDefault[nPtBinsCorrelations] = {2.090, 2.090, 2.090, 2.090, 2.090, 2.090, 2.090, 2.090};
auto vecsignalRegionLeft = std::vector<double>{signalRegionLeftDefault, signalRegionLeftDefault + nPtBinsCorrelations};
auto vecsignalRegionRight = std::vector<double>{signalRegionRightDefault, signalRegionRightDefault + nPtBinsCorrelations};
auto vecSidebandLeftInner = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + nPtBinsCorrelations};
auto vecSidebandLeftOuter = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + nPtBinsCorrelations};
auto vecSidebandRightInner = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + nPtBinsCorrelations};
auto vecSidebandRightOuter = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + nPtBinsCorrelations};
const int nPtBinsEfficiency = o2::analysis::hf_cuts_d0_to_pi_k::nBinsPt;
const double efficiencyDmesonDefault[nPtBinsEfficiency] = {};
auto vecEfficiencyDmeson = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + nPtBinsEfficiency};
const double ptHadronMax = 10.0;

struct HfTaskCorrelationD0Hadrons {

  // pT ranges for correlation plots: the default values are those embedded in hf_cuts_d0_to_pi_k (i.e. the mass pT bins), but can be redefined via json files
  Configurable<std::vector<double>> binsCorrelations{"ptBinsForCorrelations", std::vector<double>{vecPtBinsCorrelations}, "pT bin limits for correlation plots"};
  // pT bins for effiencies: same as above
  Configurable<std::vector<double>> binsEfficiency{"ptBinsForEfficiency", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for efficiency"};
  // signal and sideband region edges, to be defined via json file (initialised to empty)
  Configurable<std::vector<double>> signalRegionLeft{"signalRegionLeft", std::vector<double>{vecsignalRegionLeft}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionRight{"signalRegionRight", std::vector<double>{vecsignalRegionRight}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{vecSidebandLeftInner}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{vecSidebandLeftOuter}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{vecSidebandRightInner}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{vecSidebandRightOuter}, "Outer values of right sideband vs pT"};
  Configurable<std::vector<double>> efficiencyDmeson{"efficiencyDmeson", std::vector<double>{vecEfficiencyDmeson}, "Efficiency values for D meson specie under study"};
  Configurable<int> applyEfficiency{"efficiencyFlagD", 1, "Flag for applying efficiency weights"};

  HistogramRegistry registry{
    "registry",
    {{"hDeltaEtaPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hCorrel2DPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hCorrel2DVsPtSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hDeltaEtaPtIntSignalRegionSoftPi", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSignalRegionSoftPi", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hCorrel2DPtIntSignalRegionSoftPi", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hCorrel2DVsPtSignalRegionSoftPi", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hDeltaEtaPtIntSidebands", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hCorrel2DPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hCorrel2DVsPtSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hDeltaEtaPtIntSidebandsSoftPi", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSidebandsSoftPi", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hCorrel2DPtIntSidebandsSoftPi", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hCorrel2DVsPtSidebandsSoftPi", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtRecSig", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DVsPtRecBg", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     // correlation histograms for MCRec for signal only
     {"hCorrel2DVsPtSignalRegionRecSig", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DVsPtSignalRegionSoftPiRecSig", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DPtIntSignalRegionRecSig", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hCorrel2DPtIntSignalRegionSoftPiRecSig", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hDeltaEtaPtIntSignalRegionRecSig", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaEtaPtIntSignalRegionSoftPiRecSig", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSignalRegionRecSig", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hDeltaPhiPtIntSignalRegionSoftPiRecSig", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hCorrel2DVsPtSidebandsRecSig", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DVsPtSidebandsSoftPiRecSig", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DPtIntSidebandsRecSig", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hCorrel2DPtIntSidebandsSoftPiRecSig", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hDeltaEtaPtIntSidebandsRecSig", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaEtaPtIntSidebandsSoftPiRecSig", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSidebandsRecSig", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hDeltaPhiPtIntSidebandsSoftPiRecSig", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     // correlation histograms for MCRec for reflection candidates only
     {"hCorrel2DVsPtSignalRegionRecRef", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DPtIntSignalRegionRecRef", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hDeltaEtaPtIntSignalRegionRecRef", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSignalRegionRecRef", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hCorrel2DVsPtSidebandsRecRef", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DPtIntSidebandsRecRef", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hDeltaEtaPtIntSidebandsRecRef", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSidebandsRecRef", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hCorrel2DVsPtSignalRegionSoftPiRecRef", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DPtIntSignalRegionSoftPiRecRef", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hDeltaEtaPtIntSignalRegionSoftPiRecRef", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSignalRegionSoftPiRecRef", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hCorrel2DVsPtSidebandsSoftPiRecRef", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DPtIntSidebandsSoftPiRecRef", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hDeltaEtaPtIntSidebandsSoftPiRecRef", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSidebandsSoftPiRecRef", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     // correlation histograms for MCRec for background candidates only
     {"hCorrel2DVsPtSignalRegionRecBg", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DPtIntSignalRegionRecBg", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hDeltaEtaPtIntSignalRegionRecBg", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSignalRegionRecBg", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hCorrel2DVsPtSidebandsRecBg", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DPtIntSidebandsRecBg", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hDeltaEtaPtIntSidebandsRecBg", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSidebandsRecBg", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hCorrel2DVsPtSignalRegionSoftPiRecBg", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DPtIntSignalRegionSoftPiRecBg", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hDeltaEtaPtIntSignalRegionSoftPiRecBg", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSignalRegionSoftPiRecBg", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     {"hCorrel2DVsPtSidebandsSoftPiRecBg", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DPtIntSidebandsSoftPiRecBg", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hDeltaEtaPtIntSidebandsSoftPiRecBg", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntSidebandsSoftPiRecBg", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
     // correlation histograms for MCGen
     {"hCorrel2DVsPtGen", stringMcParticles + stringDeltaPhi + stringDeltaEta + stringPtD + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}}, // note: axes 3 and 4 (the pT) are updated in the init(),
     {"hCorrel2DPtIntGen", stringMcParticles + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
     {"hDeltaEtaPtIntGen", stringMcParticles + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
     {"hDeltaPhiPtIntGen", stringMcParticles + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}}}};

  void init(InitContext&)
  {
    int nBinsPtAxis = binsCorrelations->size() - 1;
    const double* valuesPtAxis = binsCorrelations->data();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionSoftPi"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsSoftPi"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionSoftPi"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsSoftPi"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionRecSig"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionSoftPiRecSig"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsRecSig"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsSoftPiRecSig"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionRecSig"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionSoftPiRecSig"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsRecSig"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsSoftPiRecSig"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionRecRef"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsRecRef"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionRecRef"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsRecRef"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionSoftPiRecRef"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsSoftPiRecRef"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionSoftPiRecRef"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsSoftPiRecRef"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionRecBg"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsRecBg"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionRecBg"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsRecBg"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionSoftPiRecBg"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsSoftPiRecBg"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionSoftPiRecBg"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsSoftPiRecBg"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtRecSig"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtRecSig"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtRecBg"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtRecBg"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtGen"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtGen"))->Sumw2();
  }

  /// D-h correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  /// Works on both USL and LS analyses pair tables
  void processData(aod::DHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      double massD = pairEntry.mD();
      double massDbar = pairEntry.mDbar();
      int signalStatus = pairEntry.signalStatus();
      int effBinD = o2::analysis::findBin(binsEfficiency, ptD);
      int ptBinD = o2::analysis::findBin(binsCorrelations, ptD);
      int poolBin = pairEntry.poolBin();

      // reject entries outside pT ranges of interest
      if (ptBinD < 0 || effBinD < 0) {
        continue;
      }
      if (ptHadron > ptHadronMax) {
        ptHadron = ptHadronMax + 0.5;
      }

      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyDmeson->at(o2::analysis::findBin(binsEfficiency, ptD))); // ***** track efficiency to be implemented *****
      }

      // reject entries outside pT ranges of interest
      if (ptBinD == -1) { // at least one particle outside accepted pT range
        continue;
      }
      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if ((massD > signalRegionLeft->at(ptBinD) && massD < signalRegionRight->at(ptBinD)) && ((signalStatus == ParticleTypeData::D0Only) || (signalStatus == ParticleTypeData::D0D0barBoth))) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }

      if ((massD > signalRegionLeft->at(ptBinD) && massD < signalRegionRight->at(ptBinD)) && ((signalStatus == ParticleTypeData::D0OnlySoftPi) || (signalStatus >= ParticleTypeData::D0D0barBothSoftPi))) {
        // in signal region, fills for soft pion only in ME
        registry.fill(HIST("hCorrel2DVsPtSignalRegionSoftPi"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionSoftPi"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionSoftPi"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionSoftPi"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && ((signalStatus == ParticleTypeData::D0barOnly) || (signalStatus == ParticleTypeData::D0D0barBoth))) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && (signalStatus >= ParticleTypeData::D0barOnlySoftPi)) {
        // in signal region, fills for soft pion only in ME
        registry.fill(HIST("hCorrel2DVsPtSignalRegionSoftPi"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionSoftPi"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionSoftPi"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionSoftPi"), deltaPhi, efficiencyWeight);
      }

      if (((massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) ||
           (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD))) &&
          ((signalStatus == ParticleTypeData::D0Only) || (signalStatus == ParticleTypeData::D0D0barBoth))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }

      if (((massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) ||
           (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD))) &&
          ((signalStatus == ParticleTypeData::D0OnlySoftPi) || (signalStatus >= ParticleTypeData::D0D0barBothSoftPi))) {
        // in sideband region, fills for soft pion only in ME
        registry.fill(HIST("hCorrel2DVsPtSidebandsSoftPi"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsSoftPi"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsSoftPi"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsSoftPi"), deltaPhi, efficiencyWeight);
      }

      if (((massDbar > sidebandLeftOuter->at(ptBinD) && massDbar < sidebandLeftInner->at(ptBinD)) ||
           (massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD))) &&
          ((signalStatus == ParticleTypeData::D0barOnly) || (signalStatus == ParticleTypeData::D0D0barBoth))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }

      if (((massDbar > sidebandLeftOuter->at(ptBinD) && massDbar < sidebandLeftInner->at(ptBinD)) ||
           (massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD))) &&
          (signalStatus >= ParticleTypeData::D0barOnlySoftPi)) {
        // in sideband region, fills for soft pion only in ME
        registry.fill(HIST("hCorrel2DVsPtSidebandsSoftPi"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsSoftPi"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsSoftPi"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsSoftPi"), deltaPhi, efficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationD0Hadrons, processData, "Process data", false);

  void processMcRec(aod::DHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      double massD = pairEntry.mD();
      double massDbar = pairEntry.mDbar();
      int signalStatus = pairEntry.signalStatus();
      int ptBinD = o2::analysis::findBin(binsCorrelations, ptD);
      int poolBin = pairEntry.poolBin();

      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyDmeson->at(o2::analysis::findBin(binsEfficiency, ptD)));
      }

      // fill correlation plots for signal/bagkground correlations
      if (pairEntry.signalStatus()) {
        registry.fill(HIST("hCorrel2DVsPtRecSig"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);

      } else {
        registry.fill(HIST("hCorrel2DVsPtRecBg"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
      }

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots

      // ---------------------- Fill plots for signal case, D0 ->1, D0bar ->8 ---------------------------------------------
      if ((massD > signalRegionLeft->at(ptBinD) && massD < signalRegionRight->at(ptBinD)) && (TESTBIT(signalStatus, ParticleTypeMcRec::D0Sig))) {
        // in signal region, tests bit ParticleTypeMcRec::D0Sig, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSignalRegionRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionRecSig"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionRecSig"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionRecSig"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && (TESTBIT(signalStatus, ParticleTypeMcRec::D0barSig))) {
        // in signal region, tests bit ParticleTypeMcRec::D0barSig, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSignalRegionRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionRecSig"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionRecSig"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionRecSig"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && (signalStatus >= BIT(ParticleTypeMcRec::SoftPi))) {
        // in signal region, fills for soft pion only for event mixing, tests bit at least ParticleTypeMcRec::SoftPi
        registry.fill(HIST("hCorrel2DVsPtSignalRegionSoftPiRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionSoftPiRecSig"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionSoftPiRecSig"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionSoftPiRecSig"), deltaPhi, efficiencyWeight);
      }

      if ((((massD > sidebandLeftOuter->at(ptBinD)) && (massD < sidebandLeftInner->at(ptBinD))) ||
           ((massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)))) &&
          ((TESTBIT(signalStatus, ParticleTypeMcRec::D0Sig)))) {
        // in sideband region, tests bit ParticleTypeMcRec::D0Sig, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSidebandsRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsRecSig"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsRecSig"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsRecSig"), deltaPhi, efficiencyWeight);
      }

      if ((((massDbar > sidebandLeftOuter->at(ptBinD)) && (massDbar < sidebandLeftInner->at(ptBinD))) ||
           ((massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD)))) &&
          (TESTBIT(signalStatus, ParticleTypeMcRec::D0barSig))) {
        // in sideband region, tests bit ParticleTypeMcRec::D0barSig, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSidebandsRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsRecSig"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsRecSig"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsRecSig"), deltaPhi, efficiencyWeight);
      }

      if ((((massDbar > sidebandLeftOuter->at(ptBinD)) && (massDbar < sidebandLeftInner->at(ptBinD))) ||
           ((massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD)))) &&
          (signalStatus >= BIT(ParticleTypeMcRec::SoftPi))) {
        // in sideband region, fills for soft pion only for event mixing, tests bit at least ParticleTypeMcRec::SoftPi
        registry.fill(HIST("hCorrel2DVsPtSidebandsSoftPiRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsSoftPiRecSig"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsSoftPiRecSig"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsSoftPiRecSig"), deltaPhi, efficiencyWeight);
      }

      // ---------------------- Fill plots for reflection case, D0 ->2, D0bar ->16 ---------------------------------------------
      if ((massD > signalRegionLeft->at(ptBinD) && massD < signalRegionRight->at(ptBinD)) && TESTBIT(signalStatus, ParticleTypeMcRec::D0Ref)) {
        // in signal region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSignalRegionRecRef"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionRecRef"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionRecRef"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionRecRef"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && TESTBIT(signalStatus, ParticleTypeMcRec::D0barRef)) {
        // in signal region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSignalRegionRecRef"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionRecRef"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionRecRef"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionRecRef"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && (signalStatus >= BIT(ParticleTypeMcRec::SoftPi))) {
        // in signal region, fills for soft pion only for event mixing, tests bit at least ParticleTypeMcRec::SoftPi
        registry.fill(HIST("hCorrel2DVsPtSignalRegionSoftPiRecRef"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionSoftPiRecRef"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionSoftPiRecRef"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionSoftPiRecRef"), deltaPhi, efficiencyWeight);
      }

      if ((((massD > sidebandLeftOuter->at(ptBinD)) && (massD < sidebandLeftInner->at(ptBinD))) ||
           ((massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)))) &&
          TESTBIT(signalStatus, ParticleTypeMcRec::D0Ref)) {
        // in sideband region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSidebandsRecRef"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsRecRef"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsRecRef"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsRecRef"), deltaPhi, efficiencyWeight);
      }

      if ((((massDbar > sidebandLeftOuter->at(ptBinD)) && (massDbar < sidebandLeftInner->at(ptBinD))) ||
           ((massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD)))) &&
          TESTBIT(signalStatus, ParticleTypeMcRec::D0barRef)) {
        // in sideband region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSidebandsRecRef"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsRecRef"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsRecRef"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsRecRef"), deltaPhi, efficiencyWeight);
      }

      if ((((massDbar > sidebandLeftOuter->at(ptBinD)) && (massDbar < sidebandLeftInner->at(ptBinD))) ||
           ((massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD)))) &&
          (signalStatus >= BIT(ParticleTypeMcRec::SoftPi))) {
        // in sideband region, fills for soft pion only for event mixing, tests bit at least ParticleTypeMcRec::SoftPi
        registry.fill(HIST("hCorrel2DVsPtSidebandsSoftPiRecRef"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsSoftPiRecRef"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsSoftPiRecRef"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsSoftPiRecRef"), deltaPhi, efficiencyWeight);
      }

      // ---------------------- Fill plots for background case, D0 ->4, D0bar ->32 ---------------------------------------------
      if ((massD > signalRegionLeft->at(ptBinD) && massD < signalRegionRight->at(ptBinD)) && TESTBIT(signalStatus, ParticleTypeMcRec::D0Bg)) {
        // in signal region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSignalRegionRecBg"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionRecBg"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionRecBg"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionRecBg"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && TESTBIT(signalStatus, ParticleTypeMcRec::D0barBg)) {
        // in signal region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSignalRegionRecBg"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionRecBg"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionRecBg"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionRecBg"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && (signalStatus >= BIT(ParticleTypeMcRec::SoftPi))) {
        // in signal region, fills for soft pion only for event mixing, tests bit at least ParticleTypeMcRec::SoftPi
        registry.fill(HIST("hCorrel2DVsPtSignalRegionSoftPiRecBg"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionSoftPiRecBg"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionSoftPiRecBg"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionSoftPiRecBg"), deltaPhi, efficiencyWeight);
      }

      if ((((massD > sidebandLeftOuter->at(ptBinD)) && (massD < sidebandLeftInner->at(ptBinD))) ||
           ((massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)))) &&
          TESTBIT(signalStatus, ParticleTypeMcRec::D0Bg)) {
        // in sideband region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSidebandsRecBg"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsRecBg"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsRecBg"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsRecBg"), deltaPhi, efficiencyWeight);
      }

      if ((((massDbar > sidebandLeftOuter->at(ptBinD)) && (massDbar < sidebandLeftInner->at(ptBinD))) ||
           ((massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD)))) &&
          TESTBIT(signalStatus, ParticleTypeMcRec::D0barBg)) {
        // in sideband region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSidebandsRecBg"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsRecBg"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsRecBg"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsRecBg"), deltaPhi, efficiencyWeight);
      }

      if ((((massDbar > sidebandLeftOuter->at(ptBinD)) && (massDbar < sidebandLeftInner->at(ptBinD))) ||
           ((massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD)))) &&
          (signalStatus >= BIT(ParticleTypeMcRec::SoftPi))) {
        // in sideband region, fills for soft pion only for event mixing, tests bit at least ParticleTypeMcRec::SoftPi
        registry.fill(HIST("hCorrel2DVsPtSidebandsSoftPiRecBg"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsSoftPiRecBg"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsSoftPiRecBg"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsSoftPiRecBg"), deltaPhi, efficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationD0Hadrons, processMcRec, "Process MC Reco mode", true);

  /// D-Hadron correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::DHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      // reject entries outside pT ranges of interest
      if (o2::analysis::findBin(binsCorrelations, ptD) < 0) {
        continue;
      }
      if (ptHadron > ptHadronMax) {
        ptHadron = ptHadronMax + 0.5;
      }

      registry.fill(HIST("hCorrel2DVsPtGen"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
      registry.fill(HIST("hCorrel2DPtIntGen"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntGen"), deltaPhi);
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationD0Hadrons, processMcGen, "Process MC Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationD0Hadrons>(cfgc)};
}
