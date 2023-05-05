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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_d0_to_pi_k;
using namespace o2::constants::math;

namespace o2::aod
{
using D0PairFull = soa::Join<aod::D0Pair, aod::D0PairRecoInfo>;
using DPlusPairFull = soa::Join<aod::DPlusPair, aod::DPlusPairRecoInfo>;
} // namespace o2::aod

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiD, double phiDbar)
{
  return RecoDecay::constrainAngle(phiDbar - phiD, -o2::constants::math::PIHalf);
}

///
/// Returns phi of candidate/particle evaluated from x and y components of segment connecting primary and secondary vertices
///
double evaluatePhiByVertex(double xVertex1, double xVertex2, double yVertex1, double yVertex2)
{
  return RecoDecay::phi(xVertex2 - xVertex1, yVertex2 - yVertex1);
}

// string definitions, used for histogram axis labels
const TString stringPtD = "#it{p}_{T}^{D} (GeV/#it{c});";
const TString stringPtDbar = "#it{p}_{T}^{Dbar} (GeV/#it{c});";
const TString stringDeltaPt = "#it{p}_{T}^{Dbar}-#it{p}_{T}^{D} (GeV/#it{c});";
const TString stringDeltaPtMaxMin = "#it{p}_{T}^{max}-#it{p}_{T}^{min} (GeV/#it{c});";
const TString stringDeltaEta = "#it{#eta}^{Dbar}-#it{#eta}^{D};";
const TString stringDeltaY = "#it{y}^{Dbar}-#it{y}^{D};";
const TString stringDeltaPhi = "#it{#varphi}^{Dbar}-#it{#varphi}^{D} (rad);";
const TString stringDDbar = "D,Dbar candidates ";
const TString stringSignal = "signal region;";
const TString stringSideband = "sidebands;";
const TString stringMCParticles = "MC gen - D,Dbar particles;";
const TString stringMCReco = "MC reco - D,Dbar candidates ";

// definition of vectors for standard ptbin and invariant mass configurables
const int npTBinsCorrelations = 8;
const double pTBinsCorrelations[npTBinsCorrelations + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
auto pTBinsCorrelations_v = std::vector<double>{pTBinsCorrelations, pTBinsCorrelations + npTBinsCorrelations + 1};
const double signalRegionInnerDefault[npTBinsCorrelations] = {1.810, 1.810, 1.810, 1.810, 1.810, 1.810, 1.810, 1.810};
const double signalRegionOuterDefault[npTBinsCorrelations] = {1.922, 1.922, 1.922, 1.922, 1.922, 1.922, 1.922, 1.922};
const double sidebandLeftInnerDefault[npTBinsCorrelations] = {1.642, 1.642, 1.642, 1.642, 1.642, 1.642, 1.642, 1.642};
const double sidebandLeftOuterDefault[npTBinsCorrelations] = {1.754, 1.754, 1.754, 1.754, 1.754, 1.754, 1.754, 1.754};
const double sidebandRightInnerDefault[npTBinsCorrelations] = {1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978};
const double sidebandRightOuterDefault[npTBinsCorrelations] = {2.090, 2.090, 2.090, 2.090, 2.090, 2.090, 2.090, 2.090};
auto signalRegionInner_v = std::vector<double>{signalRegionInnerDefault, signalRegionInnerDefault + npTBinsCorrelations};
auto signalRegionOuter_v = std::vector<double>{signalRegionOuterDefault, signalRegionOuterDefault + npTBinsCorrelations};
auto sidebandLeftInner_v = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + npTBinsCorrelations};
auto sidebandLeftOuter_v = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + npTBinsCorrelations};
auto sidebandRightInner_v = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + npTBinsCorrelations};
auto sidebandRightOuter_v = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + npTBinsCorrelations};

struct HfTaskCorrelationDMesonPairs {
  // pT ranges for correlation plots: the default values are those embedded in hf_cuts_d0_to_pi_k (i.e. the mass pT bins), but can be redefined via json files
  Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{pTBinsCorrelations_v}, "pT bin limits for correlation plots"};
  // signal and sideband region edges, to be defined via json file (initialised to empty)
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", std::vector<double>{signalRegionInner_v}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", std::vector<double>{signalRegionOuter_v}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{sidebandLeftInner_v}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{sidebandLeftOuter_v}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{sidebandRightInner_v}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{sidebandRightOuter_v}, "Outer values of right sideband vs pT"};
  Configurable<int> pairType{"pairType", 0, "Pair type: 0 = D0D0, 1=D0D0bar, 2 = D0barD0bar"};

  HistogramRegistry registry{
    "registry",
    // NOTE: use hMassD0 (from correlator task) for normalisation, and hMass2DCorrelationPairs for 2D-sideband-subtraction purposes
    {{"hMass2DCorrelationPairs", stringDDbar + "2D;inv. mass D (GeV/#it{c}^{2});inv. mass Dbar (GeV/#it{c}^{2});" + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hDeltaEtaPtIntSignalRegion", stringDDbar + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
     {"hDeltaYPtIntSignalRegion", stringDDbar + stringSignal + stringDeltaY + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
     {"hDeltaPhiPtIntSignalRegion", stringDDbar + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
     {"hCorrel2DPtIntSignalRegion", stringDDbar + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
     {"hCorrel2DVsPtSignalRegion", stringDDbar + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hDeltaPtDDbarSignalRegion", stringDDbar + stringSignal + stringDeltaPt + "entries", {HistType::kTH1F, {{144, -36., 36.}}}},
     {"hDeltaPtMaxMinSignalRegion", stringDDbar + stringSignal + stringDeltaPtMaxMin + "entries", {HistType::kTH1F, {{72, 0., 36.}}}},
     {"hDeltaEtaPtIntSidebands", stringDDbar + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
     {"hDeltaYPtIntSidebands", stringDDbar + stringSideband + stringDeltaY + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
     {"hDeltaPhiPtIntSidebands", stringDDbar + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
     {"hCorrel2DPtIntSidebands", stringDDbar + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
     {"hCorrel2DVsPtSidebands", stringDDbar + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hDeltaPtDDbarSidebands", stringDDbar + stringSideband + stringDeltaPt + "entries", {HistType::kTH1F, {{144, -36., 36.}}}},
     {"hDeltaPtMaxMinSidebands", stringDDbar + stringSideband + stringDeltaPtMaxMin + "entries", {HistType::kTH1F, {{72, 0., 36.}}}},
     {"hMass2DCorrelationPairsMCRecBkgBkg", stringDDbar + "2D BkgBkg - MC reco;inv. mass D (GeV/#it{c}^{2});inv. mass Dbar (GeV/#it{c}^{2});" + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hMass2DCorrelationPairsMCRecBkgRef", stringDDbar + "2D BkgRef - MC reco;inv. mass D (GeV/#it{c}^{2});inv. mass Dbar (GeV/#it{c}^{2});" + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hMass2DCorrelationPairsMCRecBkgSig", stringDDbar + "2D BkgSig - MC reco;inv. mass D (GeV/#it{c}^{2});inv. mass Dbar (GeV/#it{c}^{2});" + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hMass2DCorrelationPairsMCRecRefBkg", stringDDbar + "2D RefBkg - MC reco;inv. mass D (GeV/#it{c}^{2});inv. mass Dbar (GeV/#it{c}^{2});" + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hMass2DCorrelationPairsMCRecRefRef", stringDDbar + "2D RefRef - MC reco;inv. mass D (GeV/#it{c}^{2});inv. mass Dbar (GeV/#it{c}^{2});" + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hMass2DCorrelationPairsMCRecRefSig", stringDDbar + "2D RefSig - MC reco;inv. mass D (GeV/#it{c}^{2});inv. mass Dbar (GeV/#it{c}^{2});" + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hMass2DCorrelationPairsMCRecSigBkg", stringDDbar + "2D SigBkg - MC reco;inv. mass D (GeV/#it{c}^{2});inv. mass Dbar (GeV/#it{c}^{2});" + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hMass2DCorrelationPairsMCRecSigRef", stringDDbar + "2D SigRef - MC reco;inv. mass D (GeV/#it{c}^{2});inv. mass Dbar (GeV/#it{c}^{2});" + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hMass2DCorrelationPairsMCRecSigSig", stringDDbar + "2D SigSig - MC reco;inv. mass D (GeV/#it{c}^{2});inv. mass Dbar (GeV/#it{c}^{2});" + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hDeltaEtaPtIntSignalRegionMCRec", stringMCReco + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
     {"hDeltaPhiPtIntSignalRegionMCRec", stringMCReco + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
     {"hDeltaYPtIntSignalRegionMCRec", stringMCReco + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
     {"hCorrel2DPtIntSignalRegionMCRec", stringMCReco + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
     {"hDeltaPtDDbarSignalRegionMCRec", stringMCReco + stringSignal + stringDeltaPt + "entries", {HistType::kTH1F, {{144, -36., 36.}}}},
     {"hDeltaPtMaxMinSignalRegionMCRec", stringMCReco + stringSignal + stringDeltaPtMaxMin + "entries", {HistType::kTH1F, {{72, 0., 36.}}}},
     {"hCorrel2DVsPtSignalRegionMCRecBkgBkg", stringMCReco + "BkgBkg" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSignalRegionMCRecBkgRef", stringMCReco + "BkgRef" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSignalRegionMCRecBkgSig", stringMCReco + "BkgSig" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSignalRegionMCRecRefBkg", stringMCReco + "RefBkg" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSignalRegionMCRecRefRef", stringMCReco + "RefRef" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSignalRegionMCRecRefSig", stringMCReco + "RefSig" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSignalRegionMCRecSigBkg", stringMCReco + "SigBkg" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSignalRegionMCRecSigRef", stringMCReco + "SigRef" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSignalRegionMCRecSigSig", stringMCReco + "SigSig" + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hDeltaEtaPtIntSidebandsMCRec", stringMCReco + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
     {"hDeltaPhiPtIntSidebandsMCRec", stringMCReco + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
     {"hCorrel2DPtIntSidebandsMCRec", stringMCReco + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
     {"hDeltaPtDDbarSidebandsMCRec", stringMCReco + stringSideband + stringDeltaPt + "entries", {HistType::kTH1F, {{144, -36., 36.}}}},
     {"hDeltaPtMaxMinSidebandsMCRec", stringMCReco + stringSideband + stringDeltaPtMaxMin + "entries", {HistType::kTH1F, {{72, 0., 36.}}}},
     {"hCorrel2DVsPtSidebandsMCRecBkgBkg", stringMCReco + "BkgBkg" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSidebandsMCRecBkgRef", stringMCReco + "BkgRef" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSidebandsMCRecBkgSig", stringMCReco + "BkgSig" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSidebandsMCRecRefBkg", stringMCReco + "RefBkg" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSidebandsMCRecRefRef", stringMCReco + "RefRef" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSidebandsMCRecRefSig", stringMCReco + "RefSig" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSidebandsMCRecSigBkg", stringMCReco + "SigBkg" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSidebandsMCRecSigRef", stringMCReco + "SigRef" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtSidebandsMCRecSigSig", stringMCReco + "SigSig" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hDeltaEtaPtIntMCGen", stringMCParticles + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
     {"hDeltaYPtIntMCGen", stringMCParticles + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
     {"hDeltaPhiPtIntMCGen", stringMCParticles + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
     {"hCorrel2DPtIntMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
     {"hCorrel2DVsPtMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtDbar + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {10, 0., 10.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hDeltaPtDDbarMCGen", stringMCParticles + stringDeltaPt + "entries", {HistType::kTH1F, {{144, -36., 36.}}}},
     {"hDeltaPtMaxMinMCGen", stringMCParticles + stringDeltaPtMaxMin + "entries", {HistType::kTH1F, {{72, 0., 36.}}}}}};

  void init(o2::framework::InitContext&)
  {
    // redefinition of pT axes for THnSparse holding correlation entries
    int nBinspTaxis = binsPtCorrelations->size() - 1;
    const double* valuespTaxis = binsPtCorrelations->data();

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

  // Possible values of candidateType:
  // |1: Sig D0              |10: Ref D0              |21: Sig D0 + Ref D0bar     |102: Bkg D0 + Sig D0bar      |201: Sig D0 + Bkg D0bar
  // |2: Sig D0bar           |12: Ref D0 + Sig D0bar  |30: Ref D0 + Ref D0bar     |120: Bkg D0 + Ref D0bar      |210: Ref D0 + Bkg D0bar
  // |3: Sig D0 + Sig D0bar  |20: Ref D0bar           |100: Bkg D0                |200: Bkg D0bar               |200: Bkg D0 + Bkg D0bar

  // Register whether our D0 candidate is Sig, Ref or Bkg
  uint getDMesonType(uint const& candidateType)
  {
    if (candidateType == 1 || candidateType == 3 || candidateType == 21 || candidateType == 201) { // Signal
      return 1;
    } else if (candidateType == 10 || candidateType == 12 || candidateType == 30 || candidateType == 210) { // Reflected
      return 2;
    } else if (candidateType == 100 || candidateType == 102 || candidateType == 120 || candidateType == 300) { // Background
      return 3;
    } else {
      return 0;
    }
  }

  // Register whether our D0bar candidate is Sig, Ref or Bkg
  uint getDMesonBarType(uint const& candidateType)
  {
    if (candidateType == 2 || candidateType == 3 || candidateType == 12 || candidateType == 102) { // Signal
      return 1;
    } else if (candidateType == 20 || candidateType == 21 || candidateType == 30 || candidateType == 120) { // Reflected
      return 2;
    } else if (candidateType == 200 || candidateType == 201 || candidateType == 210 || candidateType == 300) { // Background
      return 3;
    } else {
      return 0;
    }
  }

  // Plot Mass correlations
  void drawMassCorrPlots(uint const& candLabel1, uint const& candLabel2, double const& massCand1, double const& massCand2, double const& ptCand1, double const& ptCand2)
  {
    if (candLabel1 == 1 && candLabel2 == 1) {
      registry.fill(HIST("hMass2DCorrelationPairsMCRecSigSig"), massCand1, massCand2, ptCand1, ptCand2);
    } else if (candLabel1 == 1 && candLabel2 == 2) {
      registry.fill(HIST("hMass2DCorrelationPairsMCRecSigRef"), massCand1, massCand2, ptCand1, ptCand2);
    } else if (candLabel1 == 1 && candLabel2 == 3) {
      registry.fill(HIST("hMass2DCorrelationPairsMCRecSigBkg"), massCand1, massCand2, ptCand1, ptCand2);
    } else if (candLabel1 == 2 && candLabel2 == 1) {
      registry.fill(HIST("hMass2DCorrelationPairsMCRecRefSig"), massCand1, massCand2, ptCand1, ptCand2);
    } else if (candLabel1 == 2 && candLabel2 == 2) {
      registry.fill(HIST("hMass2DCorrelationPairsMCRecRefRef"), massCand1, massCand2, ptCand1, ptCand2);
    } else if (candLabel1 == 2 && candLabel2 == 3) {
      registry.fill(HIST("hMass2DCorrelationPairsMCRecRefBkg"), massCand1, massCand2, ptCand1, ptCand2);
    } else if (candLabel1 == 3 && candLabel2 == 1) {
      registry.fill(HIST("hMass2DCorrelationPairsMCRecBkgSig"), massCand1, massCand2, ptCand1, ptCand2);
    } else if (candLabel1 == 3 && candLabel2 == 2) {
      registry.fill(HIST("hMass2DCorrelationPairsMCRecBkgRef"), massCand1, massCand2, ptCand1, ptCand2);
    } else if (candLabel1 == 3 && candLabel2 == 3) {
      registry.fill(HIST("hMass2DCorrelationPairsMCRecBkgBkg"), massCand1, massCand2, ptCand1, ptCand2);
    }
  }

  // Plot angular correlations in signal region
  void drawCorrelSignalPlots(uint const& candLabel1, uint const& candLabel2, double const& deltaPhi, double const& deltaEta, double const& ptCand1, double const& ptCand2)
  {
    if (candLabel1 == 1 && candLabel2 == 1) {
      registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRecSigSig"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 1 && candLabel2 == 2) {
      registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRecSigRef"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 1 && candLabel2 == 3) {
      registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRecSigBkg"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 2 && candLabel2 == 1) {
      registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRecRefSig"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 2 && candLabel2 == 2) {
      registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRecRefRef"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 2 && candLabel2 == 3) {
      registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRecRefBkg"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 3 && candLabel2 == 1) {
      registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRecBkgSig"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 3 && candLabel2 == 2) {
      registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRecBkgRef"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 3 && candLabel2 == 3) {
      registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRecBkgBkg"), deltaPhi, deltaEta, ptCand1, ptCand2);
    }
  }

  // Plot angular correlations in sideband region
  void drawCorrelSidebandsPlots(uint const& candLabel1, uint const& candLabel2, double const& deltaPhi, double const& deltaEta, double const& ptCand1, double const& ptCand2)
  {
    if (candLabel1 == 1 && candLabel2 == 1) {
      registry.fill(HIST("hCorrel2DVsPtSidebandsMCRecSigSig"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 1 && candLabel2 == 2) {
      registry.fill(HIST("hCorrel2DVsPtSidebandsMCRecSigRef"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 1 && candLabel2 == 3) {
      registry.fill(HIST("hCorrel2DVsPtSidebandsMCRecSigBkg"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 2 && candLabel2 == 1) {
      registry.fill(HIST("hCorrel2DVsPtSidebandsMCRecRefSig"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 2 && candLabel2 == 2) {
      registry.fill(HIST("hCorrel2DVsPtSidebandsMCRecRefRef"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 2 && candLabel2 == 3) {
      registry.fill(HIST("hCorrel2DVsPtSidebandsMCRecRefBkg"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 3 && candLabel2 == 1) {
      registry.fill(HIST("hCorrel2DVsPtSidebandsMCRecBkgSig"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 3 && candLabel2 == 2) {
      registry.fill(HIST("hCorrel2DVsPtSidebandsMCRecBkgRef"), deltaPhi, deltaEta, ptCand1, ptCand2);
    } else if (candLabel1 == 3 && candLabel2 == 3) {
      registry.fill(HIST("hCorrel2DVsPtSidebandsMCRecBkgBkg"), deltaPhi, deltaEta, ptCand1, ptCand2);
    }
  }

  template <typename T>
  void analyseData(const T& pairEntries)
  {
    for (auto& pairEntry : pairEntries) {
      if (pairEntry.isData() == 1) {
        // define variables for widely used quantities
        double deltaPhi = pairEntry.deltaPhi();
        double deltaEta = pairEntry.deltaEta();
        double ptCand1 = pairEntry.ptCand1();
        double ptCand2 = pairEntry.ptCand2();
        double massCand1 = pairEntry.mCand1();
        double massCand2 = pairEntry.mCand2();
        double yCand1 = pairEntry.yCand1();
        double yCand2 = pairEntry.yCand2();

        int pTBinCand1 = o2::analysis::findBin(binsPtCorrelations, ptCand1);
        int pTBinCand2 = o2::analysis::findBin(binsPtCorrelations, ptCand2);

        // fill 2D invariant mass plots
        registry.fill(HIST("hMass2DCorrelationPairs"), massCand1, massCand2, ptCand1, ptCand2);

        // reject entries outside pT ranges of interest
        if (pTBinCand1 == -1 || pTBinCand2 == -1) { // at least one particle outside accepted pT range
          continue;
        }

        // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
        if (massCand1 > signalRegionInner->at(pTBinCand1) && massCand1 < signalRegionOuter->at(pTBinCand1) && massCand2 > signalRegionInner->at(pTBinCand2) && massCand2 < signalRegionOuter->at(pTBinCand2)) {
          // in signal region
          registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptCand1, ptCand2);
          registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta);
          registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta);
          registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi);
          registry.fill(HIST("hDeltaYPtIntSignalRegion"), yCand2 - yCand1);
          registry.fill(HIST("hDeltaPtDDbarSignalRegion"), ptCand2 - ptCand1);
          registry.fill(HIST("hDeltaPtMaxMinSignalRegion"), std::abs(ptCand2 - ptCand1));
        }

        if ((massCand1 > sidebandLeftInner->at(pTBinCand1) && massCand1 < sidebandLeftOuter->at(pTBinCand1) && massCand2 > sidebandLeftInner->at(pTBinCand2) && massCand2 < sidebandRightOuter->at(pTBinCand2)) ||
            (massCand1 > sidebandRightInner->at(pTBinCand1) && massCand1 < sidebandRightOuter->at(pTBinCand1) && massCand2 > sidebandLeftInner->at(pTBinCand2) && massCand2 < sidebandRightOuter->at(pTBinCand2)) ||
            (massCand1 > sidebandLeftInner->at(pTBinCand1) && massCand1 < sidebandRightOuter->at(pTBinCand1) && massCand2 > sidebandLeftInner->at(pTBinCand2) && massCand2 < sidebandLeftOuter->at(pTBinCand2)) ||
            (massCand1 > sidebandLeftInner->at(pTBinCand1) && massCand1 < sidebandRightOuter->at(pTBinCand1) && massCand2 > sidebandRightInner->at(pTBinCand2) && massCand2 < sidebandRightOuter->at(pTBinCand2))) {
          // in sideband region
          registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptCand1, ptCand2);
          registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta);
          registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta);
          registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi);
          registry.fill(HIST("hDeltaPtDDbarSidebands"), ptCand2 - ptCand1);
          registry.fill(HIST("hDeltaYPtIntSidebands"), yCand2 - yCand1);
          registry.fill(HIST("hDeltaPtMaxMinSidebands"), std::abs(ptCand2 - ptCand1));
        }
      }
    }
  }

  template <typename T>
  void analyseMcRec(const T& pairEntries)
  {
    for (auto& pairEntry : pairEntries) {
      if (pairEntry.isData() == 0 && pairEntry.isReco() == 1) {
        // define variables for widely used quantities
        double deltaPhi = pairEntry.deltaPhi();
        double deltaEta = pairEntry.deltaEta();
        double ptCand1 = pairEntry.ptCand1();
        double ptCand2 = pairEntry.ptCand2();
        double massCand1 = pairEntry.mCand1();
        double massCand2 = pairEntry.mCand2();
        double yCand1 = pairEntry.yCand1();
        double yCand2 = pairEntry.yCand2();

        int pTBinCand1 = o2::analysis::findBin(binsPtCorrelations, ptCand1);
        int pTBinCand2 = o2::analysis::findBin(binsPtCorrelations, ptCand2);

        // fill 2D invariant mass plots
        uint dMesonCand1 = 0, dMesonCand2 = 0;
        uint dMesonBarCand1 = 0, dMesonBarCand2 = 0;
        // 0: default; 1: Signal; 2: Reflected; 3: Background
        dMesonCand1 = getDMesonType(pairEntry.candidateType1());
        dMesonBarCand1 = getDMesonBarType(pairEntry.candidateType1());
        dMesonCand2 = getDMesonType(pairEntry.candidateType2());
        dMesonBarCand2 = getDMesonBarType(pairEntry.candidateType2());

        switch (pairType) {
          case 0: // D0 D0
            drawMassCorrPlots(dMesonCand1, dMesonCand2, massCand1, massCand2, ptCand1, ptCand2);
            break;
          case 1: // D0 D0bar
            drawMassCorrPlots(dMesonCand1, dMesonBarCand2, massCand1, massCand2, ptCand1, ptCand2);
            break;
          case 2: // D0bar D0bar
            drawMassCorrPlots(dMesonBarCand1, dMesonBarCand2, massCand1, massCand2, ptCand1, ptCand2);
            break;
        }

        // reject entries outside pT ranges of interest
        if (pTBinCand1 == -1 || pTBinCand2 == -1) { // at least one particle outside accepted pT range
          continue;
        }

        // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
        if (massCand1 > signalRegionInner->at(pTBinCand1) && massCand1 < signalRegionOuter->at(pTBinCand1) && massCand2 > signalRegionInner->at(pTBinCand2) && massCand2 < signalRegionOuter->at(pTBinCand2)) {
          // in signal region
          registry.fill(HIST("hCorrel2DPtIntSignalRegionMCRec"), deltaPhi, deltaEta);
          registry.fill(HIST("hDeltaEtaPtIntSignalRegionMCRec"), deltaEta);
          registry.fill(HIST("hDeltaPhiPtIntSignalRegionMCRec"), deltaPhi);
          registry.fill(HIST("hDeltaYPtIntSignalRegionMCRec"), yCand2 - yCand1);
          registry.fill(HIST("hDeltaPtDDbarSignalRegionMCRec"), ptCand2 - ptCand1);
          registry.fill(HIST("hDeltaPtMaxMinSignalRegionMCRec"), std::abs(ptCand2 - ptCand1));

          switch (pairType) {
            case 0: // D0 D0
              drawCorrelSignalPlots(dMesonCand1, dMesonCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
              break;
            case 1: // D0 D0bar
              drawCorrelSignalPlots(dMesonCand1, dMesonBarCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
              break;
            case 2: // D0bar D0bar
              drawCorrelSignalPlots(dMesonBarCand1, dMesonBarCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
              break;
          }
        }

        if ((massCand1 > sidebandLeftInner->at(pTBinCand1) && massCand1 < sidebandLeftOuter->at(pTBinCand1) && massCand2 > sidebandLeftInner->at(pTBinCand2) && massCand2 < sidebandRightOuter->at(pTBinCand2)) ||
            (massCand1 > sidebandRightInner->at(pTBinCand1) && massCand1 < sidebandRightOuter->at(pTBinCand1) && massCand2 > sidebandLeftInner->at(pTBinCand2) && massCand2 < sidebandRightOuter->at(pTBinCand2)) ||
            (massCand1 > sidebandLeftInner->at(pTBinCand1) && massCand1 < sidebandRightOuter->at(pTBinCand1) && massCand2 > sidebandLeftInner->at(pTBinCand2) && massCand2 < sidebandLeftOuter->at(pTBinCand2)) ||
            (massCand1 > sidebandLeftInner->at(pTBinCand1) && massCand1 < sidebandRightOuter->at(pTBinCand1) && massCand2 > sidebandRightInner->at(pTBinCand2) && massCand2 < sidebandRightOuter->at(pTBinCand2))) {
          // in sideband region
          registry.fill(HIST("hCorrel2DPtIntSidebandsMCRec"), deltaPhi, deltaEta);
          registry.fill(HIST("hDeltaEtaPtIntSidebandsMCRec"), deltaEta);
          registry.fill(HIST("hDeltaPhiPtIntSidebandsMCRec"), deltaPhi);
          registry.fill(HIST("hDeltaPtDDbarSidebandsMCRec"), ptCand2 - ptCand1);
          registry.fill(HIST("hDeltaPtMaxMinSidebandsMCRec"), std::abs(ptCand2 - ptCand1));

          switch (pairType) {
            case 0: // D0 D0
              drawCorrelSidebandsPlots(dMesonCand1, dMesonCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
              break;
            case 1: // D0 D0bar
              drawCorrelSidebandsPlots(dMesonCand1, dMesonBarCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
              break;
            case 2: // D0bar D0bar
              drawCorrelSidebandsPlots(dMesonBarCand1, dMesonBarCand2, deltaPhi, deltaEta, ptCand1, ptCand2);
              break;
          }
        }
      }
    } // end loop
  }

  template <typename T>
  void analyseMcGen(const T& pairEntries)
  {
    for (auto& pairEntry : pairEntries) {
      if (pairEntry.isData() == 0 && pairEntry.isReco() == 0) {
        // define variables for widely used quantities
        double deltaPhi = pairEntry.deltaPhi();
        double deltaEta = pairEntry.deltaEta();
        double ptCand1 = pairEntry.ptCand1();
        double ptCand2 = pairEntry.ptCand2();
        double yCand1 = pairEntry.yCand1();
        double yCand2 = pairEntry.yCand2();

        // reject entries outside pT ranges of interest
        if (o2::analysis::findBin(binsPtCorrelations, ptCand1) == -1 || o2::analysis::findBin(binsPtCorrelations, ptCand2) == -1) {
          continue;
        }

        registry.fill(HIST("hCorrel2DVsPtMCGen"), deltaPhi, deltaEta, ptCand1, ptCand2);
        registry.fill(HIST("hCorrel2DPtIntMCGen"), deltaPhi, deltaEta);
        registry.fill(HIST("hDeltaEtaPtIntMCGen"), deltaEta);
        registry.fill(HIST("hDeltaPhiPtIntMCGen"), deltaPhi);
        registry.fill(HIST("hDeltaYPtIntMCGen"), yCand2 - yCand1);
        registry.fill(HIST("hDeltaPtDDbarMCGen"), ptCand2 - ptCand1);
        registry.fill(HIST("hDeltaPtMaxMinMCGen"), std::abs(ptCand2 - ptCand2));
      }
    }
  } // end loop

  ///
  /// PROCESS FUNCTIONS
  ///
  void processDataD0(aod::D0Pair const& pairEntries)
  {
    analyseData(pairEntries);
  }
  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processDataD0, "Process data D0", true);

  void processDataDPlus(aod::DPlusPair const& pairEntries)
  {
    analyseData(pairEntries);
  }

  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processDataDPlus, "Process data Dplus", false);

  /// D-Dbar correlation pair filling task, from pair tables - for MC reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRecD0(aod::D0PairFull const& pairEntries)
  {
    analyseMcRec(pairEntries);
  }

  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processMcRecD0, "Process MC D0 Rec mode", true);

  void processMcRecDPlus(aod::DPlusPairFull const& pairEntries)
  {
    analyseMcRec(pairEntries);
  }

  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processMcRecDPlus, "Process MC Dplus Reco mode", false);

  /// D-Dbar correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal) - Ok for both USL and LS analyses
  void processMcGenD0(aod::D0PairFull const& pairEntries)
  {
    analyseMcGen(pairEntries);
  }

  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processMcGenD0, "Process MC D0 Gen mode", false);

  void processMcGenDPlus(aod::DPlusPairFull const& pairEntries)
  {
    analyseMcGen(pairEntries);
  }

  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processMcGenDPlus, "Process MC DPlus Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationDMesonPairs>(cfgc)};
}
