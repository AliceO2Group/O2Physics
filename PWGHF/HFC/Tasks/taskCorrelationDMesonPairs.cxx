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

/// \file taskCorrelationDMesonPairs.cxx
/// \brief D-Dbar analysis task - data-like, MC-reco and MC-kine analyses for D0 and DPlus pairs.
///
/// \author Fabio Colamaria <fabio.colamaria@ba.infn.it>, INFN Bari
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>, IJCLab Orsay

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/HFC/DataModel/DMesonPairsTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
using D0PairFull = soa::Join<aod::D0Pair, aod::D0PairRecoInfo>;
using DplusPairFull = soa::Join<aod::DplusPair, aod::DplusPairRecoInfo>;
} // namespace o2::aod

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiD, double phiDbar)
{
  return RecoDecay::constrainAngle(phiDbar - phiD, -o2::constants::math::PIHalf);
}

namespace
{
enum CandidateTypeSel {
  SelectedD = 0, // This particle is selected as a D
  SelectedDbar,  // This particle is selected as a Dbar
  TrueD,         // This particle is a true D
  TrueDbar       // This particle is a true Dbar
};

enum DMesonType {
  Default = 0, // Default value
  Signal,      // This particle is a signal D meson
  Reflected,   // This particle is a reflected D meson
  Bkg          // This particle is background of D meson
};

enum PairTypeSel {
  DD = 0,  // Analyse D0-D0 or DPlus-DPlus correlations
  DDbar,   // Analyse D0-D0bar or DPlus-DMinus correlations
  DbarDbar // Analyse D0bar-D0bar or DMinus-DMinus correlations
};
} // namespace

///
/// Returns phi of candidate/particle evaluated from x and y components of segment connecting primary and secondary vertices
///
double evaluatePhiByVertex(double xVertex1, double xVertex2, double yVertex1, double yVertex2)
{
  return RecoDecay::phi(xVertex2 - xVertex1, yVertex2 - yVertex1);
}

// string definitions, used for histogram axis labels
const TString stringPtD = "#it{p}_{T}^{D_{1}} (GeV/#it{c});";
const TString stringPtDbar = "#it{p}_{T}^{D_{2}} (GeV/#it{c});";
const TString stringDeltaPt = "#it{p}_{T}^{D_{2}}-#it{p}_{T}^{D_{1}} (GeV/#it{c});";
const TString stringDeltaPtMaxMin = "#it{p}_{T}^{max}-#it{p}_{T}^{min} (GeV/#it{c});";
const TString stringDeltaEta = "#it{#eta}^{D_{2}}-#it{#eta}^{D_{1}};";
const TString stringDeltaY = "#it{y}^{D_{2}}-#it{y}^{D_{1}};";
const TString stringDeltaPhi = "#it{#varphi}^{D_{2}}-#it{#varphi}^{D_{1}} (rad);";
const TString stringDDbar = "D meson pair candidates ";
const TString stringSignal = "signal region;";
const TString stringSideband = "sidebands;";
const TString stringMCParticles = "MC gen - D meson pair particles;";
const TString stringMCReco = "MC reco - D meson pair candidates ";
const TString stringMassD = "inv. mass D_{1} (GeV/#it{c}^{2});";
const TString stringMassDbar = "inv. mass D_{2} (GeV/#it{c}^{2});";

// definition of vectors for standard ptbin and invariant mass configurables
const int nPtBinsCorrelations = 8;
const double ptBinsCorrelations[nPtBinsCorrelations + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
auto vecPtBinsCorrelations = std::vector<double>{ptBinsCorrelations, ptBinsCorrelations + nPtBinsCorrelations + 1};
const double signalRegionInnerDefault[nPtBinsCorrelations] = {1.810, 1.810, 1.810, 1.810, 1.810, 1.810, 1.810, 1.810};
const double signalRegionOuterDefault[nPtBinsCorrelations] = {1.922, 1.922, 1.922, 1.922, 1.922, 1.922, 1.922, 1.922};
const double sidebandLeftInnerDefault[nPtBinsCorrelations] = {1.642, 1.642, 1.642, 1.642, 1.642, 1.642, 1.642, 1.642};
const double sidebandLeftOuterDefault[nPtBinsCorrelations] = {1.754, 1.754, 1.754, 1.754, 1.754, 1.754, 1.754, 1.754};
const double sidebandRightInnerDefault[nPtBinsCorrelations] = {1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978};
const double sidebandRightOuterDefault[nPtBinsCorrelations] = {2.090, 2.090, 2.090, 2.090, 2.090, 2.090, 2.090, 2.090};
auto vecSignalRegionInner = std::vector<double>{signalRegionInnerDefault, signalRegionInnerDefault + nPtBinsCorrelations};
auto vecSignalRegionOuter = std::vector<double>{signalRegionOuterDefault, signalRegionOuterDefault + nPtBinsCorrelations};
auto vecSidebandLeftInner = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + nPtBinsCorrelations};
auto vecSidebandLeftOuter = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + nPtBinsCorrelations};
auto vecSidebandRightInner = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + nPtBinsCorrelations};
auto vecSidebandRightOuter = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + nPtBinsCorrelations};

struct HfTaskCorrelationDMesonPairs {
  // Enable histograms with finer pT and y binning
  Configurable<bool> enableFinerBinning{"enableFinerBinning", false, "Enable histograms with finer pT and y binning"};
  // pT ranges for correlation plots: the default values are those embedded in hf_cuts_d0_to_pi_k (i.e. the mass pT bins), but can be redefined via json files
  Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{vecPtBinsCorrelations}, "pT bin limits for correlation plots"};
  // signal and sideband region edges, to be defined via json file (initialised to empty)
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", std::vector<double>{vecSignalRegionInner}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", std::vector<double>{vecSignalRegionOuter}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{vecSidebandLeftInner}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{vecSidebandLeftOuter}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{vecSidebandRightInner}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{vecSidebandRightOuter}, "Outer values of right sideband vs pT"};
  Configurable<int> pairType{"pairType", 0, "Pair type: 0 = DD, 1=DDbar, 2 = DbarDbar"};

  // HistoTypes
  HistogramConfigSpec hTHnMass2DCorrPairs{HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {10, 0., 10.}, {10, 0., 10.}, {10, -1, 1}, {10, -1, 1}}};                        // note: axes 3 and 4 (the pT) are updated in the init();
  HistogramConfigSpec hTHnCorrel2DVsPt{HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}; // note: axes 3 and 4 (the pT) are updated in the init()
  HistogramConfigSpec hTH1Y{HistType::kTH1F, {{200, -10., 10.}}};
  HistogramConfigSpec hTH1DeltaPtDDbar{HistType::kTH1F, {{144, -36., 36.}}};
  HistogramConfigSpec hTH1DeltaPtMaxMin{HistType::kTH1F, {{72, 0., 36.}}};
  HistogramConfigSpec hTH1Phi{HistType::kTH1F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}}};
  HistogramConfigSpec hTH2CorrelPt{HistType::kTH2F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {200, -10., 10.}}};
  HistogramConfigSpec hTHnMass2DCorrPairsFinerBinning{HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {60, 1., 6.}, {60, 1., 6.}, {160, -0.8, 0.8}, {160, -0.8, 0.8}}};

  HistogramRegistry registry{
    "registry",
    // NOTE: use hMassD0 (from correlator task) for normalisation, and hMass2DCorrelationPairs for 2D-sideband-subtraction purposes
    {{"hMass2DCorrelationPairs", stringDDbar + "2D;inv. mass D (GeV/#it{c}^{2});inv. mass Dbar (GeV/#it{c}^{2});" + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairs},
     {"hDeltaEtaPtIntSignalRegion", stringDDbar + stringSignal + stringDeltaEta + "entries", hTH1Y},
     {"hDeltaYPtIntSignalRegion", stringDDbar + stringSignal + stringDeltaY + "entries", hTH1Y},
     {"hDeltaPhiPtIntSignalRegion", stringDDbar + stringSignal + stringDeltaPhi + "entries", hTH1Phi},
     {"hCorrel2DPtIntSignalRegion", stringDDbar + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", hTH2CorrelPt},
     {"hCorrel2DVsPtSignalRegion", stringDDbar + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hDeltaPtDDbarSignalRegion", stringDDbar + stringSignal + stringDeltaPt + "entries", hTH1DeltaPtDDbar},
     {"hDeltaPtMaxMinSignalRegion", stringDDbar + stringSignal + stringDeltaPtMaxMin + "entries", hTH1DeltaPtMaxMin},
     {"hDeltaEtaPtIntSidebands", stringDDbar + stringSideband + stringDeltaEta + "entries", hTH1Y},
     {"hDeltaYPtIntSidebands", stringDDbar + stringSideband + stringDeltaY + "entries", hTH1Y},
     {"hDeltaPhiPtIntSidebands", stringDDbar + stringSideband + stringDeltaPhi + "entries", hTH1Phi},
     {"hCorrel2DPtIntSidebands", stringDDbar + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", hTH2CorrelPt},
     {"hCorrel2DVsPtSidebands", stringDDbar + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hDeltaPtDDbarSidebands", stringDDbar + stringSideband + stringDeltaPt + "entries", hTH1DeltaPtDDbar},
     {"hDeltaPtMaxMinSidebands", stringDDbar + stringSideband + stringDeltaPtMaxMin + "entries", hTH1DeltaPtMaxMin},
     {"hMass2DCorrelationPairsMCRecBkgBkg", stringDDbar + "2D BkgBkg - MC reco;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairs},
     {"hMass2DCorrelationPairsMCRecBkgRef", stringDDbar + "2D BkgRef - MC reco;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairs},
     {"hMass2DCorrelationPairsMCRecBkgSig", stringDDbar + "2D BkgSig - MC reco;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairs},
     {"hMass2DCorrelationPairsMCRecRefBkg", stringDDbar + "2D RefBkg - MC reco;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairs},
     {"hMass2DCorrelationPairsMCRecRefRef", stringDDbar + "2D RefRef - MC reco;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairs},
     {"hMass2DCorrelationPairsMCRecRefSig", stringDDbar + "2D RefSig - MC reco;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairs},
     {"hMass2DCorrelationPairsMCRecSigBkg", stringDDbar + "2D SigBkg - MC reco;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairs},
     {"hMass2DCorrelationPairsMCRecSigRef", stringDDbar + "2D SigRef - MC reco;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairs},
     {"hMass2DCorrelationPairsMCRecSigSig", stringDDbar + "2D SigSig - MC reco;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairs},
     {"hDeltaEtaPtIntSignalRegionMCRec", stringMCReco + stringSignal + stringDeltaEta + "entries", hTH1Y},
     {"hDeltaPhiPtIntSignalRegionMCRec", stringMCReco + stringSignal + stringDeltaPhi + "entries", hTH1Phi},
     {"hDeltaYPtIntSignalRegionMCRec", stringMCReco + stringSignal + stringDeltaY + "entries", hTH1Y},
     {"hCorrel2DPtIntSignalRegionMCRec", stringMCReco + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", hTH2CorrelPt},
     {"hDeltaPtDDbarSignalRegionMCRec", stringMCReco + stringSignal + stringDeltaPt + "entries", hTH1DeltaPtDDbar},
     {"hDeltaPtMaxMinSignalRegionMCRec", stringMCReco + stringSignal + stringDeltaPtMaxMin + "entries", hTH1DeltaPtMaxMin},
     {"hCorrel2DVsPtSignalRegionMCRecBkgBkg", stringMCReco + "BkgBkg" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSignalRegionMCRecBkgRef", stringMCReco + "BkgRef" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSignalRegionMCRecBkgSig", stringMCReco + "BkgSig" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSignalRegionMCRecRefBkg", stringMCReco + "RefBkg" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSignalRegionMCRecRefRef", stringMCReco + "RefRef" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSignalRegionMCRecRefSig", stringMCReco + "RefSig" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSignalRegionMCRecSigBkg", stringMCReco + "SigBkg" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSignalRegionMCRecSigRef", stringMCReco + "SigRef" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSignalRegionMCRecSigSig", stringMCReco + "SigSig" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hDeltaEtaPtIntSidebandsMCRec", stringMCReco + stringSideband + stringDeltaEta + "entries", hTH1Y},
     {"hDeltaPhiPtIntSidebandsMCRec", stringMCReco + stringSideband + stringDeltaPhi + "entries", hTH1Phi},
     {"hCorrel2DPtIntSidebandsMCRec", stringMCReco + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", hTH2CorrelPt},
     {"hDeltaPtDDbarSidebandsMCRec", stringMCReco + stringSideband + stringDeltaPt + "entries", hTH1DeltaPtDDbar},
     {"hDeltaPtMaxMinSidebandsMCRec", stringMCReco + stringSideband + stringDeltaPtMaxMin + "entries", hTH1DeltaPtMaxMin},
     {"hCorrel2DVsPtSidebandsMCRecBkgBkg", stringMCReco + "BkgBkg" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSidebandsMCRecBkgRef", stringMCReco + "BkgRef" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSidebandsMCRecBkgSig", stringMCReco + "BkgSig" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSidebandsMCRecRefBkg", stringMCReco + "RefBkg" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSidebandsMCRecRefRef", stringMCReco + "RefRef" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSidebandsMCRecRefSig", stringMCReco + "RefSig" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSidebandsMCRecSigBkg", stringMCReco + "SigBkg" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSidebandsMCRecSigRef", stringMCReco + "SigRef" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hCorrel2DVsPtSidebandsMCRecSigSig", stringMCReco + "SigSig" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hDeltaEtaPtIntMCGen", stringMCParticles + stringDeltaEta + "entries", hTH1Y},
     {"hDeltaYPtIntMCGen", stringMCParticles + stringDeltaY + "entries", hTH1Y},
     {"hDeltaPhiPtIntMCGen", stringMCParticles + stringDeltaPhi + "entries", hTH1Phi},
     {"hCorrel2DPtIntMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + "entries", hTH2CorrelPt},
     {"hCorrel2DVsPtMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", hTHnCorrel2DVsPt},
     {"hDeltaPtDDbarMCGen", stringMCParticles + stringDeltaPt + "entries", hTH1DeltaPtDDbar},
     {"hDeltaPtMaxMinMCGen", stringMCParticles + stringDeltaPtMaxMin + "entries", hTH1DeltaPtMaxMin}}};

  void init(InitContext&)
  {
    // redefinition of pT axes for THnSparse holding correlation entries
    int nBinspTaxis = binsPtCorrelations->size() - 1;
    const double* valuespTaxis = binsPtCorrelations->data();

    if (enableFinerBinning) {
      registry.add("hMass2DCorrelationPairsFinerBinning", stringDDbar + "2D Finer Binning;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairsFinerBinning);
      registry.add("hMass2DCorrelationPairsMCRecBkgBkgFinerBinning", stringDDbar + "2D BkgBkg - MC reco Finer Binning;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairsFinerBinning);
      registry.add("hMass2DCorrelationPairsMCRecBkgRefFinerBinning", stringDDbar + "2D BkgBkg - MC reco Finer Binning;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairsFinerBinning);
      registry.add("hMass2DCorrelationPairsMCRecBkgSigFinerBinning", stringDDbar + "2D BkgBkg - MC reco Finer Binning;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairsFinerBinning);
      registry.add("hMass2DCorrelationPairsMCRecRefBkgFinerBinning", stringDDbar + "2D BkgBkg - MC reco Finer Binning;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairsFinerBinning);
      registry.add("hMass2DCorrelationPairsMCRecRefRefFinerBinning", stringDDbar + "2D BkgBkg - MC reco Finer Binning;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairsFinerBinning);
      registry.add("hMass2DCorrelationPairsMCRecRefSigFinerBinning", stringDDbar + "2D BkgBkg - MC reco Finer Binning;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairsFinerBinning);
      registry.add("hMass2DCorrelationPairsMCRecSigBkgFinerBinning", stringDDbar + "2D BkgBkg - MC reco Finer Binning;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairsFinerBinning);
      registry.add("hMass2DCorrelationPairsMCRecSigRefFinerBinning", stringDDbar + "2D BkgBkg - MC reco Finer Binning;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairsFinerBinning);
      registry.add("hMass2DCorrelationPairsMCRecSigSigFinerBinning", stringDDbar + "2D BkgBkg - MC reco Finer Binning;" + stringMassD + stringMassDbar + stringPtD + stringPtDbar + "entries", hTHnMass2DCorrPairsFinerBinning);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsFinerBinning"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgBkgFinerBinning"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgRefFinerBinning"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgSigFinerBinning"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefBkgFinerBinning"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefRefFinerBinning"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefSigFinerBinning"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigBkgFinerBinning"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigRefFinerBinning"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigSigFinerBinning"))->Sumw2();
    }

    for (int i = 2; i <= 3; i++) {
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairs"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairs"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgBkg"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgRef"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgSig"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefBkg"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefRef"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefSig"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigBkg"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigRef"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigSig"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecBkgBkg"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecBkgRef"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecBkgSig"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecRefBkg"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecRefRef"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecRefSig"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecSigBkg"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecSigRef"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecSigSig"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecBkgBkg"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecBkgRef"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecBkgSig"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecRefBkg"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecRefRef"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecRefSig"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecSigBkg"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecSigRef"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecSigSig"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgBkg"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgRef"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgSig"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefBkg"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefRef"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefSig"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigBkg"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigRef"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigSig"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecBkgBkg"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecBkgRef"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecBkgSig"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecRefBkg"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecRefRef"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecRefSig"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecSigBkg"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecSigRef"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecSigSig"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecBkgBkg"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecBkgRef"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecBkgSig"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecRefBkg"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecRefRef"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecRefSig"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecSigBkg"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecSigRef"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecSigSig"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtMCGen"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtMCGen"))->Sumw2();
    }
  }
  ///
  /// FUNCTIONS
  ///

  // Register whether our D meson candidate is Sig, Ref or Bkg
  uint getDMesonType(uint const& candidateType)
  {
    if (TESTBIT(candidateType, SelectedD) && TESTBIT(candidateType, TrueD) && !TESTBIT(candidateType, TrueDbar)) { // Signal
      return Signal;
    } else if (TESTBIT(candidateType, SelectedD) && TESTBIT(candidateType, TrueDbar)) { // Reflected
      return Reflected;
    } else if (TESTBIT(candidateType, SelectedD) && !(TESTBIT(candidateType, TrueD) && TESTBIT(candidateType, TrueDbar))) { // Background
      return Bkg;
    }
    return Default;
  }

  // Register whether our Dbar meson candidate is Sig, Ref or Bkg
  uint getDMesonBarType(uint const& candidateType)
  {
    if (TESTBIT(candidateType, SelectedDbar) && TESTBIT(candidateType, TrueDbar) && !TESTBIT(candidateType, TrueD)) { // Signal
      return Signal;
    } else if (TESTBIT(candidateType, SelectedDbar) && TESTBIT(candidateType, TrueD)) { // Reflected
      return Reflected;
    } else if (TESTBIT(candidateType, SelectedDbar) && !(TESTBIT(candidateType, TrueD) && TESTBIT(candidateType, TrueDbar))) { // Background
      return Bkg;
    }
    return Default;
  }

  // Fill Mass correlation histograms
  void fillMassCorrHists(std::shared_ptr<THnSparse> hMassCorrArray[3][3], uint const& candLabel1, uint const& candLabel2, double const& massCand1, double const& massCand2, double const& ptCand1, double const& ptCand2, double const& yCand1, double const& yCand2)
  {
    if (candLabel1 != 0 && candLabel2 != 0) {
      hMassCorrArray[candLabel1 - 1][candLabel2 - 1]->Fill(massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
    }
  }

  // Fill angular correlation histograms
  void fillAngularCorrelHists(std::shared_ptr<THnSparse> hCorrelArray[3][3], uint const& candLabel1, uint const& candLabel2, double const& deltaPhi, double const& deltaEta, double const& ptCand1, double const& ptCand2)
  {
    if (candLabel1 != 0 && candLabel2 != 0) {
      hCorrelArray[candLabel1 - 1][candLabel2 - 1]->Fill(deltaPhi, deltaEta, ptCand1, ptCand2);
    }
  }

  // Fill kinematic pair info histograms in signal region
  void fillKinematicSignalHists(int const& dataType, double const& yCand1, double const& yCand2, double const& deltaPhi, double const& deltaEta, double const& ptCand1, double const& ptCand2)
  {
    if (dataType == 0) { // Data
      registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptCand1, ptCand2);
      registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi);
      registry.fill(HIST("hDeltaYPtIntSignalRegion"), yCand2 - yCand1);
      registry.fill(HIST("hDeltaPtDDbarSignalRegion"), ptCand2 - ptCand1);
      registry.fill(HIST("hDeltaPtMaxMinSignalRegion"), std::abs(ptCand2 - ptCand1));
    } else if (dataType == 1) { // MC Reco
      registry.fill(HIST("hCorrel2DPtIntSignalRegionMCRec"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntSignalRegionMCRec"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntSignalRegionMCRec"), deltaPhi);
      registry.fill(HIST("hDeltaYPtIntSignalRegionMCRec"), yCand2 - yCand1);
      registry.fill(HIST("hDeltaPtDDbarSignalRegionMCRec"), ptCand2 - ptCand1);
      registry.fill(HIST("hDeltaPtMaxMinSignalRegionMCRec"), std::abs(ptCand2 - ptCand1));
    } else { // MC gen
      registry.fill(HIST("hCorrel2DVsPtMCGen"), deltaPhi, deltaEta, ptCand1, ptCand2);
      registry.fill(HIST("hCorrel2DPtIntMCGen"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntMCGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMCGen"), deltaPhi);
      registry.fill(HIST("hDeltaYPtIntMCGen"), yCand2 - yCand1);
      registry.fill(HIST("hDeltaPtDDbarMCGen"), ptCand2 - ptCand1);
      registry.fill(HIST("hDeltaPtMaxMinMCGen"), std::abs(ptCand2 - ptCand2));
    }
  }

  // Fill kinematic pair info histograms in sideband region
  void fillKinematicSidebandHists(int const& dataType, double const& yCand1, double const& yCand2, double const& deltaPhi, double const& deltaEta, double const& ptCand1, double const& ptCand2)
  {
    if (dataType == 0) { // Data
      registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptCand1, ptCand2);
      registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi);
      registry.fill(HIST("hDeltaPtDDbarSidebands"), ptCand2 - ptCand1);
      registry.fill(HIST("hDeltaYPtIntSidebands"), yCand2 - yCand1);
      registry.fill(HIST("hDeltaPtMaxMinSidebands"), std::abs(ptCand2 - ptCand1));
    } else { // MC Reco
      registry.fill(HIST("hCorrel2DPtIntSidebandsMCRec"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntSidebandsMCRec"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntSidebandsMCRec"), deltaPhi);
      registry.fill(HIST("hDeltaPtDDbarSidebandsMCRec"), ptCand2 - ptCand1);
      registry.fill(HIST("hDeltaPtMaxMinSidebandsMCRec"), std::abs(ptCand2 - ptCand1));
    }
  }

  // Common code to analyse correlations at data level
  template <typename T>
  void analyseData(const T& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      if (pairEntry.dataType() != 0) {
        continue;
      }
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptCand1 = pairEntry.ptCand1();
      double ptCand2 = pairEntry.ptCand2();
      double massCand1 = pairEntry.mCand1();
      double massCand2 = pairEntry.mCand2();
      double yCand1 = pairEntry.yCand1();
      double yCand2 = pairEntry.yCand2();

      int ptBinCand1 = o2::analysis::findBin(binsPtCorrelations, ptCand1);
      int ptBinCand2 = o2::analysis::findBin(binsPtCorrelations, ptCand2);

      // fill 2D invariant mass plots
      if (pairType == DD && (TESTBIT(pairEntry.candidateType1(), SelectedD) && TESTBIT(pairEntry.candidateType2(), SelectedD))) {
        registry.fill(HIST("hMass2DCorrelationPairs"), massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
        if (enableFinerBinning) {
          registry.fill(HIST("hMass2DCorrelationPairsFinerBinning"), massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
        }
      } else if (pairType == DDbar && ((TESTBIT(pairEntry.candidateType1(), SelectedD) && TESTBIT(pairEntry.candidateType2(), SelectedDbar)) || (TESTBIT(pairEntry.candidateType1(), SelectedDbar) && TESTBIT(pairEntry.candidateType2(), SelectedD)))) {
        registry.fill(HIST("hMass2DCorrelationPairs"), massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
        if (enableFinerBinning) {
          registry.fill(HIST("hMass2DCorrelationPairsFinerBinning"), massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
        }
      } else if (pairType == DbarDbar && (TESTBIT(pairEntry.candidateType1(), SelectedDbar) && TESTBIT(pairEntry.candidateType2(), SelectedDbar))) {
        registry.fill(HIST("hMass2DCorrelationPairs"), massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
        if (enableFinerBinning) {
          registry.fill(HIST("hMass2DCorrelationPairsFinerBinning"), massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
        }
      }

      // reject entries outside pT ranges of interest
      if (ptBinCand1 == -1 || ptBinCand2 == -1) { // at least one particle outside accepted pT range
        continue;
      }

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massCand1 > signalRegionInner->at(ptBinCand1) && massCand1 < signalRegionOuter->at(ptBinCand1) && massCand2 > signalRegionInner->at(ptBinCand2) && massCand2 < signalRegionOuter->at(ptBinCand2)) {
        // in signal region
        if (pairType == DD && (TESTBIT(pairEntry.candidateType1(), SelectedD) && TESTBIT(pairEntry.candidateType2(), SelectedD))) {
          fillKinematicSignalHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
        } else if (pairType == DDbar && ((TESTBIT(pairEntry.candidateType1(), SelectedD) && TESTBIT(pairEntry.candidateType2(), SelectedDbar)) || (TESTBIT(pairEntry.candidateType1(), SelectedDbar) && TESTBIT(pairEntry.candidateType2(), SelectedD)))) {
          fillKinematicSignalHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
        } else if (pairType == DbarDbar && (TESTBIT(pairEntry.candidateType1(), SelectedDbar) && TESTBIT(pairEntry.candidateType2(), SelectedDbar))) {
          fillKinematicSignalHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
        }
      }

      bool leftSidebandCand1 = massCand1 > sidebandLeftInner->at(ptBinCand1) && massCand1 < sidebandLeftOuter->at(ptBinCand1);
      bool rightSidebandCand1 = massCand1 > sidebandRightInner->at(ptBinCand1) && massCand1 < sidebandRightOuter->at(ptBinCand1);
      bool leftSidebandCand2 = massCand2 > sidebandLeftInner->at(ptBinCand2) && massCand2 < sidebandLeftOuter->at(ptBinCand2);
      bool rightSidebandCand2 = massCand2 > sidebandRightInner->at(ptBinCand2) && massCand2 < sidebandRightOuter->at(ptBinCand2);

      if ((leftSidebandCand1 || rightSidebandCand1) && (leftSidebandCand2 || rightSidebandCand2)) {
        // in sideband region
        if (pairType == DD && (TESTBIT(pairEntry.candidateType1(), SelectedD) && TESTBIT(pairEntry.candidateType2(), SelectedD))) {
          fillKinematicSidebandHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
        } else if (pairType == DDbar && ((TESTBIT(pairEntry.candidateType1(), SelectedD) && TESTBIT(pairEntry.candidateType2(), SelectedDbar)) || (TESTBIT(pairEntry.candidateType1(), SelectedDbar) && TESTBIT(pairEntry.candidateType2(), SelectedD)))) {
          fillKinematicSidebandHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
        } else if (pairType == DbarDbar && (TESTBIT(pairEntry.candidateType1(), SelectedDbar) && TESTBIT(pairEntry.candidateType2(), SelectedDbar))) {
          fillKinematicSidebandHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
        }
      }
    }
  }

  // Common code to analyse correlations at Mc reco level
  template <typename T>
  void analyseMcRec(const T& pairEntries)
  {
    // Array definitions to later be used to fill histograms
    std::shared_ptr<THnSparse> hMassCorrArray[3][3] = {{registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigSig")), registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigRef")), registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigBkg"))},
                                                       {registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefSig")), registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefRef")), registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefBkg"))},
                                                       {registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgSig")), registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgRef")), registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgBkg"))}};

    std::shared_ptr<THnSparse> hCorrelSignalArray[3][3] = {{registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecSigSig")), registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecSigRef")), registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecSigBkg"))},
                                                           {registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecRefSig")), registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecRefRef")), registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecRefBkg"))},
                                                           {registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecBkgSig")), registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecBkgRef")), registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRecBkgBkg"))}};

    std::shared_ptr<THnSparse> hCorrelSidebandsArray[3][3] = {{registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecSigSig")), registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecSigRef")), registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecSigBkg"))},
                                                              {registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecRefSig")), registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecRefRef")), registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecRefBkg"))},
                                                              {registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecBkgSig")), registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecBkgRef")), registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRecBkgBkg"))}};

    std::shared_ptr<THnSparse> hMassCorrArrayFinerBinning[3][3] = {{registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigSigFinerBinning")), registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigRefFinerBinning")), registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecSigBkgFinerBinning"))},
                                                                   {registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefSigFinerBinning")), registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefRefFinerBinning")), registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecRefBkgFinerBinning"))},
                                                                   {registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgSigFinerBinning")), registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgRefFinerBinning")), registry.get<THnSparse>(HIST("hMass2DCorrelationPairsMCRecBkgBkgFinerBinning"))}};

    for (const auto& pairEntry : pairEntries) {
      if (pairEntry.dataType() != 1) { // Assure that we only analyse Mc reco elements
        continue;
      }
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptCand1 = pairEntry.ptCand1();
      double ptCand2 = pairEntry.ptCand2();
      double massCand1 = pairEntry.mCand1();
      double massCand2 = pairEntry.mCand2();
      double yCand1 = pairEntry.yCand1();
      double yCand2 = pairEntry.yCand2();

      int ptBinCand1 = o2::analysis::findBin(binsPtCorrelations, ptCand1);
      int ptBinCand2 = o2::analysis::findBin(binsPtCorrelations, ptCand2);

      uint dMesonCand1 = 0, dMesonCand2 = 0;
      uint dMesonBarCand1 = 0, dMesonBarCand2 = 0;
      // 0: default; 1: Signal; 2: Reflected; 3: Background
      dMesonCand1 = getDMesonType(pairEntry.candidateType1());
      dMesonBarCand1 = getDMesonBarType(pairEntry.candidateType1());
      dMesonCand2 = getDMesonType(pairEntry.candidateType2());
      dMesonBarCand2 = getDMesonBarType(pairEntry.candidateType2());

      // fill 2D invariant mass plots
      switch (pairType) {
        case DD: // D0 D0
          fillMassCorrHists(hMassCorrArray, dMesonCand1, dMesonCand2, massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
          if (enableFinerBinning) {
            fillMassCorrHists(hMassCorrArrayFinerBinning, dMesonCand1, dMesonCand2, massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
          }
          break;
        case DDbar: // D0 D0bar
          fillMassCorrHists(hMassCorrArray, dMesonCand1, dMesonBarCand2, massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
          fillMassCorrHists(hMassCorrArray, dMesonBarCand1, dMesonCand2, massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
          if (enableFinerBinning) {
            fillMassCorrHists(hMassCorrArrayFinerBinning, dMesonCand1, dMesonBarCand2, massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
            fillMassCorrHists(hMassCorrArrayFinerBinning, dMesonBarCand1, dMesonCand2, massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
          }
          break;
        case DbarDbar: // D0bar D0bar
          fillMassCorrHists(hMassCorrArray, dMesonBarCand1, dMesonBarCand2, massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
          if (enableFinerBinning) {
            fillMassCorrHists(hMassCorrArrayFinerBinning, dMesonBarCand1, dMesonBarCand2, massCand1, massCand2, ptCand1, ptCand2, yCand1, yCand2);
          }
          break;
      }

      // reject entries outside pT ranges of interest
      if (ptBinCand1 == -1 || ptBinCand2 == -1) { // at least one particle outside accepted pT range
        continue;
      }

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massCand1 > signalRegionInner->at(ptBinCand1) && massCand1 < signalRegionOuter->at(ptBinCand1) && massCand2 > signalRegionInner->at(ptBinCand2) && massCand2 < signalRegionOuter->at(ptBinCand2)) {
        // in signal region
        // Fill histograms depending on the type of pair we are analysing
        if (pairType == DD && (dMesonCand1 != 0 && dMesonCand2 != 0)) {
          fillKinematicSignalHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
        } else if (pairType == DDbar && ((dMesonCand1 != 0 && dMesonBarCand2 != 0) || (dMesonBarCand1 != 0 && dMesonCand2 != 0))) {
          fillKinematicSignalHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
        } else if (pairType == DbarDbar && (dMesonBarCand1 != 0 && dMesonBarCand2 != 0)) {
          fillKinematicSignalHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
        }

        // fill 2D angular correlations plots
        switch (pairType) {
          case DD: // D0 D0
            fillAngularCorrelHists(hCorrelSignalArray, dMesonCand1, dMesonCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
            break;
          case DDbar: // D0 D0bar
            fillAngularCorrelHists(hCorrelSignalArray, dMesonCand1, dMesonBarCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
            fillAngularCorrelHists(hCorrelSignalArray, dMesonBarCand1, dMesonCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
            break;
          case DbarDbar: // D0bar D0bar
            fillAngularCorrelHists(hCorrelSignalArray, dMesonBarCand1, dMesonBarCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
            break;
        }
      }

      bool leftSidebandCand1 = massCand1 > sidebandLeftInner->at(ptBinCand1) && massCand1 < sidebandLeftOuter->at(ptBinCand1);
      bool rightSidebandCand1 = massCand1 > sidebandRightInner->at(ptBinCand1) && massCand1 < sidebandRightOuter->at(ptBinCand1);
      bool leftSidebandCand2 = massCand2 > sidebandLeftInner->at(ptBinCand2) && massCand2 < sidebandLeftOuter->at(ptBinCand2);
      bool rightSidebandCand2 = massCand2 > sidebandRightInner->at(ptBinCand2) && massCand2 < sidebandRightOuter->at(ptBinCand2);

      if ((leftSidebandCand1 || rightSidebandCand1) && (leftSidebandCand2 || rightSidebandCand2)) {
        // in sideband region
        // Fill histograms depending on the type of pair we are analysing
        if (pairType == DD && (dMesonCand1 != 0 && dMesonCand2 != 0)) {
          fillKinematicSidebandHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
        } else if (pairType == DDbar && ((dMesonCand1 != 0 && dMesonBarCand2 != 0) || (dMesonBarCand1 != 0 && dMesonCand2 != 0))) {
          fillKinematicSidebandHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
        } else if (pairType == DbarDbar && (dMesonBarCand1 != 0 && dMesonBarCand2 != 0)) {
          fillKinematicSidebandHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
        }

        // fill 2D angular correlations plots
        switch (pairType) {
          case DD: // D0 D0
            fillAngularCorrelHists(hCorrelSidebandsArray, dMesonCand1, dMesonCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
            break;
          case DDbar: // D0 D0bar
            fillAngularCorrelHists(hCorrelSidebandsArray, dMesonCand1, dMesonBarCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
            fillAngularCorrelHists(hCorrelSidebandsArray, dMesonBarCand1, dMesonCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
            break;
          case DbarDbar: // D0bar D0bar
            fillAngularCorrelHists(hCorrelSidebandsArray, dMesonBarCand1, dMesonBarCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
            break;
        }
      }
    } // end loop
  }

  // Common code to analyse correlations at Mc gen level
  template <typename T>
  void analyseMcGen(const T& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      if (pairEntry.dataType() != 2) { // Assure that we only have Mc gen elements
        continue;
      }
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptCand1 = pairEntry.ptCand1();
      double ptCand2 = pairEntry.ptCand2();
      double yCand1 = pairEntry.yCand1();
      double yCand2 = pairEntry.yCand2();

      uint dMesonCand1 = 0, dMesonCand2 = 0;
      uint dMesonBarCand1 = 0, dMesonBarCand2 = 0;
      // 0: default; 1: Signal; 2: Reflected; 3: Background
      dMesonCand1 = getDMesonType(pairEntry.candidateType1());
      dMesonBarCand1 = getDMesonBarType(pairEntry.candidateType1());
      dMesonCand2 = getDMesonType(pairEntry.candidateType2());
      dMesonBarCand2 = getDMesonBarType(pairEntry.candidateType2());

      // reject entries outside pT ranges of interest
      if (o2::analysis::findBin(binsPtCorrelations, ptCand1) == -1 || o2::analysis::findBin(binsPtCorrelations, ptCand2) == -1) {
        continue;
      }

      // Fill histograms depending on the type of pair we are analysing
      if (pairType == DD && (dMesonCand1 != 0 && dMesonCand2 != 0)) {
        fillKinematicSignalHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
      } else if (pairType == DDbar && ((dMesonCand1 != 0 && dMesonBarCand2 != 0) || (dMesonBarCand1 != 0 && dMesonCand2 != 0))) {
        fillKinematicSignalHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
      } else if (pairType == DbarDbar && (dMesonBarCand1 != 0 && dMesonBarCand2 != 0)) {
        fillKinematicSignalHists(pairEntry.dataType(), yCand1, yCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
      }
    }
  } // end loop

  ///
  /// PROCESS FUNCTIONS
  ///
  /// D0-D0bar correlation pair filling task, from pair tables - for data-level analysis (no filter/selection, only true signal)
  void processDataD0(aod::D0Pair const& pairEntries)
  {
    analyseData(pairEntries);
  }
  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processDataD0, "Process data D0", true);

  /// DPlus-DMinus correlation pair filling task, from pair tables - for data-level analysis (no filter/selection, only true signal)
  void processDataDPlus(aod::DplusPair const& pairEntries)
  {
    analyseData(pairEntries);
  }

  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processDataDPlus, "Process data Dplus", false);

  /// D0-D0bar correlation pair filling task, from pair tables - for MC reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRecD0(aod::D0PairFull const& pairEntries)
  {
    analyseMcRec(pairEntries);
  }

  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processMcRecD0, "Process MC D0 Rec mode", false);

  /// DPlus-DMinus correlation pair filling task, from pair tables - for MC reco-level analysis (no filter/selection, only true signal)
  void processMcRecDPlus(aod::DplusPairFull const& pairEntries)
  {
    analyseMcRec(pairEntries);
  }

  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processMcRecDPlus, "Process MC Dplus Reco mode", false);

  /// D0-D0bar correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGenD0(aod::D0PairFull const& pairEntries)
  {
    analyseMcGen(pairEntries);
  }

  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processMcGenD0, "Process MC D0 Gen mode", false);

  /// DPlus-DMinus correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGenDPlus(aod::DplusPairFull const& pairEntries)
  {
    analyseMcGen(pairEntries);
  }

  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processMcGenDPlus, "Process MC DPlus Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationDMesonPairs>(cfgc)};
}
