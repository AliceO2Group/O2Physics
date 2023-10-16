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

/// \file taskCorrelationLcHadrons.cxx
/// \author Marianna Mazzilli <marianna.mazzilli@cern.ch>
/// \author Zhen Zhang <zhenz@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
using LcHadronPairFull = soa::Join<aod::LcHadronPair, aod::LcHadronRecoInfo>;
} // namespace o2::aod

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiLc, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiLc, -o2::constants::math::PIHalf);
}

///
/// Returns phi of candidate/particle evaluated from x and y components of segment connecting primary and secondary vertices
///
double evaluatePhiByVertex(double xVertex1, double xVertex2, double yVertex1, double yVertex2)
{
  return RecoDecay::phi(xVertex2 - xVertex1, yVertex2 - yVertex1);
}

// string definitions, used for histogram axis labels
const TString stringPtLc = "#it{p}_{T}^{#Lambda_c} (GeV/#it{c});";
const TString stringPtHadron = "#it{p}_{T}^{Hadron} (GeV/#it{c});";
const TString stringDeltaEta = "#it{#eta}^{Hadron}-#it{#eta}^{#Lambda_c};";
const TString stringDeltaPhi = "#it{#varphi}^{Hadron}-#it{#varphi}^{#Lambda_c} (rad);";
const TString stringDHadron = "#Lambda_c,Hadron candidates ";
const TString stringSignal = "signal region;";
const TString stringSideband = "sidebands;";
const TString stringMcParticles = "Mc gen - #Lambda_c,Hadron particles;";
const TString stringMcReco = "Mc reco - #Lambda_c,Hadron candidates ";

const int npTBinsCorrelations = 8;
const double pTBinsCorrelations[npTBinsCorrelations + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
auto pTBinsCorrelations_v = std::vector<double>{pTBinsCorrelations, pTBinsCorrelations + npTBinsCorrelations + 1};
const double signalRegionInnerDefault[npTBinsCorrelations] = {2.269, 2.269, 2.269, 2.269, 2.269, 2.269, 2.269, 2.269};
const double signalRegionOuterDefault[npTBinsCorrelations] = {2.309, 2.309, 2.309, 2.309, 2.309, 2.309, 2.309, 2.309};
const double sidebandLeftOuterDefault[npTBinsCorrelations] = {2.209, 2.209, 2.209, 2.209, 2.209, 2.209, 2.209, 2.209};
const double sidebandLeftInnerDefault[npTBinsCorrelations] = {2.249, 2.249, 2.249, 2.249, 2.249, 2.249, 2.249, 2.249};
const double sidebandRightInnerDefault[npTBinsCorrelations] = {2.329, 2.329, 2.329, 2.329, 2.329, 2.329, 2.329, 2.329};
const double sidebandRightOuterDefault[npTBinsCorrelations] = {2.369, 2.369, 2.369, 2.369, 2.369, 2.369, 2.369, 2.369};
auto signalRegionInner_v = std::vector<double>{signalRegionInnerDefault, signalRegionInnerDefault + npTBinsCorrelations};
auto signalRegionOuter_v = std::vector<double>{signalRegionOuterDefault, signalRegionOuterDefault + npTBinsCorrelations};
auto sidebandLeftInner_v = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + npTBinsCorrelations};
auto sidebandLeftOuter_v = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + npTBinsCorrelations};
auto sidebandRightInner_v = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + npTBinsCorrelations};
auto sidebandRightOuter_v = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + npTBinsCorrelations};
const int npTBinsEfficiency = o2::analysis::hf_cuts_lc_to_p_k_pi::nBinsPt;
const double efficiencyLcDefault[npTBinsEfficiency] = {};
auto efficiencyLc_v = std::vector<double>{efficiencyLcDefault, efficiencyLcDefault + npTBinsEfficiency};

/// Lc-Hadron correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via Mc truth)
struct HfTaskCorrelationLcHadrons {
  // pT ranges for correlation plots: the default values are those embedded in hf_cuts_lc_to_p_k_pi (i.e. the mass pT bins), but can be redefined via json files
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying efficiency weights"};
  Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{pTBinsCorrelations_v}, "pT bin limits for correlation plots"};
  Configurable<std::vector<double>> binsPtEfficiency{"binsPtEfficiency", std::vector<double>{o2::analysis::hf_cuts_lc_to_p_k_pi::vecBinsPt}, "pT bin limits for efficiency"};
  // signal and sideband region edges, to be defined via json file (initialised to empty)
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", std::vector<double>{signalRegionInner_v}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", std::vector<double>{signalRegionOuter_v}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{sidebandLeftInner_v}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{sidebandLeftOuter_v}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{sidebandRightInner_v}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{sidebandRightOuter_v}, "Outer values of right sideband vs pT"};
  Configurable<std::vector<double>> efficiencyLc{"efficiencyLc", std::vector<double>{efficiencyLc_v}, "Efficiency values for Lc "};

  HistogramRegistry registry{
    "registry",
    {
      {"hDeltaEtaPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -4., 4.}}}},
      {"hDeltaPhiPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -4., 4.}}}},
      {"hCorrel2DVsPtSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -4., 4.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
      {"hDeltaEtaPtIntSidebands", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -4., 4.}}}},
      {"hDeltaPhiPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -4., 4.}}}},
      {"hCorrel2DVsPtSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -4., 4.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
      {"hDeltaEtaPtIntSignalRegionMcRec", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -4., 4.}}}},
      {"hDeltaPhiPtIntSignalRegionMcRec", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DPtIntSignalRegionMcRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -4., 4.}}}},
      {"hCorrel2DVsPtSignalRegionMcRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -4., 4.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}},
      {"hCorrel2DVsPtSignalMcRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -4., 4.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}},
      {"hCorrel2DVsPtBkgMcRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -4., 4.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}},
      {"hDeltaEtaPtIntSidebandsMcRec", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -4., 4.}}}},
      {"hDeltaPhiPtIntSidebandsMcRec", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DPtIntSidebandsMcRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -4., 4.}}}},
      {"hCorrel2DVsPtSidebandsMcRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -4., 4.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}},
      {"hDeltaEtaPtIntMcGen", stringMcParticles + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -4., 4.}}}},
      {"hDeltaPhiPtIntMcGen", stringMcParticles + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DPtIntMcGen", stringMcParticles + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -4., 4.}}}},
      {"hCorrel2DVsPtMcGen", stringMcParticles + stringDeltaPhi + stringDeltaEta + stringPtLc + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -4., 4.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
    }};

  void init(o2::framework::InitContext&)
  {
    // redefinition of pT axes for THnSparse holding correlation entries
    int nBinspTaxis = binsPtCorrelations->size() - 1;
    const double* valuespTaxis = binsPtCorrelations->data();

    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMcRec"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMcRec"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMcRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMcRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalMcRec"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalMcRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtBkgMcRec"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtBkgMcRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGen"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGen"))->Sumw2();
  }

  void processData(aod::LcHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptLc = pairEntry.ptLc();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      double massLc = pairEntry.mLc();
      int effBinLc = o2::analysis::findBin(binsPtEfficiency, ptLc);
      int pTBinLc = o2::analysis::findBin(binsPtCorrelations, ptLc);
      // reject entries outside pT ranges of interest
      if (pTBinLc < 0 || effBinLc < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }
      double efficiencyWeight = 1.;
      double efficiencyHadron = 1.; // Note: To be implemented later on
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyLc->at(effBinLc) * efficiencyHadron);
      }
      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massLc > signalRegionInner->at(pTBinLc) && massLc < signalRegionOuter->at(pTBinLc)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }

      if ((massLc > sidebandLeftOuter->at(pTBinLc) && massLc < sidebandLeftInner->at(pTBinLc)) ||
          (massLc > sidebandRightInner->at(pTBinLc) && massLc < sidebandRightOuter->at(pTBinLc))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationLcHadrons, processData, "Process data", false);

  /// Lc-Hadron correlation pair filling task, from pair tables - for Mc reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRec(aod::LcHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptLc = pairEntry.ptLc();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      double massLc = pairEntry.mLc();
      int effBinLc = o2::analysis::findBin(binsPtEfficiency, ptLc);
      int pTBinLc = o2::analysis::findBin(binsPtCorrelations, ptLc);
      if (pTBinLc < 0 || effBinLc < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }
      double efficiencyWeight = 1.;
      double efficiencyHadron = 1.; // Note: To be implemented later on
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyLc->at(effBinLc) * efficiencyHadron);
      }
      // fill correlation plots for signal/bagkground correlations
      if (pairEntry.signalStatus()) {
        registry.fill(HIST("hCorrel2DVsPtSignalMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
      } else {
        registry.fill(HIST("hCorrel2DVsPtBkgMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
      }
      // reject entries outside pT ranges of interest

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massLc > signalRegionInner->at(pTBinLc) && massLc < signalRegionOuter->at(pTBinLc)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegionMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionMcRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionMcRec"), deltaPhi, efficiencyWeight);
      }

      if (((massLc > sidebandLeftOuter->at(pTBinLc)) && (massLc < sidebandLeftInner->at(pTBinLc))) ||
          ((massLc > sidebandRightInner->at(pTBinLc) && massLc < sidebandRightOuter->at(pTBinLc)))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebandsMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsMcRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsMcRec"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationLcHadrons, processMcRec, "Process Mc Reco mode", true);

  /// Lc-Hadron correlation pair filling task, from pair tables - for Mc gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::LcHadronPair const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptLc = pairEntry.ptLc();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      // reject entries outside pT ranges of interest
      if (o2::analysis::findBin(binsPtCorrelations, ptLc) < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }

      registry.fill(HIST("hCorrel2DVsPtMcGen"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin);
      registry.fill(HIST("hCorrel2DPtIntMcGen"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntMcGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMcGen"), deltaPhi);
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationLcHadrons, processMcGen, "Process Mc Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationLcHadrons>(cfgc)};
}
