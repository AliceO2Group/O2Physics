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

/// \file taskCorrelationDplusHadrons.cxx
/// \author Shyam Kumar <shyam.kumar@cern.ch>

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

namespace o2::aod
{
using DplusHadronPairFull = soa::Join<aod::DplusHadronPair, aod::DplusHadronRecoInfo>;
} // namespace o2::aod

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiD, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PIHalf);
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
const TString stringPtHadron = "#it{p}_{T}^{Hadron} (GeV/#it{c});";
const TString stringDeltaEta = "#it{#eta}^{Hadron}-#it{#eta}^{D};";
const TString stringDeltaPhi = "#it{#varphi}^{Hadron}-#it{#varphi}^{D} (rad);";
const TString stringDHadron = "D,Hadron candidates ";
const TString stringSignal = "signal region;";
const TString stringSideband = "sidebands;";
const TString stringMCParticles = "MC gen - D,Hadron particles;";
const TString stringMCReco = "MC reco - D,Hadron candidates ";

const int npTBinsCorrelations = 8;
const double pTBinsCorrelations[npTBinsCorrelations + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
auto pTBinsCorrelations_v = std::vector<double>{pTBinsCorrelations, pTBinsCorrelations + npTBinsCorrelations + 1};
const double signalRegionInnerDefault[npTBinsCorrelations] = {1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490};
const double signalRegionOuterDefault[npTBinsCorrelations] = {1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890};
const double sidebandLeftOuterDefault[npTBinsCorrelations] = {1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690};
const double sidebandLeftInnerDefault[npTBinsCorrelations] = {1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250};
const double sidebandRightInnerDefault[npTBinsCorrelations] = {1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130};
const double sidebandRightOuterDefault[npTBinsCorrelations] = {1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690};
auto signalRegionInner_v = std::vector<double>{signalRegionInnerDefault, signalRegionInnerDefault + npTBinsCorrelations};
auto signalRegionOuter_v = std::vector<double>{signalRegionOuterDefault, signalRegionOuterDefault + npTBinsCorrelations};
auto sidebandLeftInner_v = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + npTBinsCorrelations};
auto sidebandLeftOuter_v = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + npTBinsCorrelations};
auto sidebandRightInner_v = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + npTBinsCorrelations};
auto sidebandRightOuter_v = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + npTBinsCorrelations};
const int npTBinsEfficiency = o2::analysis::hf_cuts_dplus_to_pi_k_pi::nBinsPt;
std::vector<double> efficiencyDmeson(npTBinsEfficiency + 1);

/// Dplus-Hadron correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfTaskCorrelationDplusHadrons {
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying efficiency weights"};
  // pT ranges for correlation plots: the default values are those embedded in hf_cuts_dplus_to_pi_k_pi (i.e. the mass pT bins), but can be redefined via json files
  Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{pTBinsCorrelations_v}, "pT bin limits for correlation plots"};
  Configurable<std::vector<double>> binsPtEfficiency{"binsPtEfficiency", std::vector<double>{o2::analysis::hf_cuts_dplus_to_pi_k_pi::vecBinsPt}, "pT bin limits for efficiency"};
  // signal and sideband region edges, to be defined via json file (initialised to empty)
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", std::vector<double>{signalRegionInner_v}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", std::vector<double>{signalRegionOuter_v}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{sidebandLeftInner_v}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{sidebandLeftOuter_v}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{sidebandRightInner_v}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{sidebandRightOuter_v}, "Outer values of right sideband vs pT"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", std::vector<double>{efficiencyDmeson}, "Efficiency values for D meson specie under study"};

  HistogramRegistry registry{
    "registry",
    {
      {"hDeltaEtaPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}}}},
      {"hCorrel2DPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {200, -10., 10.}}}},
      {"hCorrel2DVsPtSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
      {"hDeltaEtaPtIntSidebands", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}}}},
      {"hCorrel2DPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {200, -10., 10.}}}},
      {"hCorrel2DVsPtSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
      {"hDeltaEtaPtIntSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}}}},
      {"hDeltaEtaPtIntSidebandsMCRec", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hCorrel2DPtIntSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {200, -10., 10.}}}},
      {"hCorrel2DVsPtSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}},
      {"hCorrel2DVsPtSignalMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}},
      {"hCorrel2DVsPtBkgMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}},
      {"hDeltaPhiPtIntSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}}}},
      {"hCorrel2DPtIntSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {200, -10., 10.}}}},
      {"hCorrel2DVsPtSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}},
      {"hDeltaEtaPtIntMCGen", stringMCParticles + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntMCGen", stringMCParticles + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}}}},
      {"hCorrel2DPtIntMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {200, -10., 10.}}}},
      {"hCorrel2DVsPtMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + stringPtD + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}, {9, 0., 9.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
    }};

  void init(InitContext&)
  {
    // redefinition of pT axes for THnSparse holding correlation entries
    int nBinspTaxis = binsPtCorrelations->size() - 1;
    const double* valuespTaxis = binsPtCorrelations->data();

    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRec"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRec"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalMCRec"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalMCRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtBkgMCRec"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtBkgMCRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMCGen"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMCGen"))->Sumw2();
  }

  void processData(aod::DplusHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      double massD = pairEntry.mD();
      int effBinD = o2::analysis::findBin(binsPtEfficiency, ptD);
      int pTBinD = o2::analysis::findBin(binsPtCorrelations, ptD);

      // reject entries outside pT ranges of interest
      if (pTBinD < 0 || effBinD < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }
      double efficiencyWeight = 1.;
      double efficiencyHadron = 1.; // Note: To be implemented later on
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(effBinD) * efficiencyHadron);
      }
      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }

      if ((massD > sidebandLeftOuter->at(pTBinD) && massD < sidebandLeftInner->at(pTBinD)) ||
          (massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusHadrons, processData, "Process data", false);

  /// D-Hadron correlation pair filling task, from pair tables - for MC reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRec(aod::DplusHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      double massD = pairEntry.mD();
      int effBinD = o2::analysis::findBin(binsPtEfficiency, ptD);
      int pTBinD = o2::analysis::findBin(binsPtCorrelations, ptD);
      if (pTBinD < 0 || effBinD < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }
      double efficiencyWeight = 1.;
      double efficiencyHadron = 1.; // Note: To be implemented later on
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(effBinD) * efficiencyHadron);
      }
      // fill correlation plots for signal/bagkground correlations
      if (pairEntry.signalStatus()) {
        registry.fill(HIST("hCorrel2DVsPtSignalMCRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
      } else {
        registry.fill(HIST("hCorrel2DVsPtBkgMCRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
      }
      // reject entries outside pT ranges of interest

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionMCRec"), deltaPhi, efficiencyWeight);
      }

      if (((massD > sidebandLeftOuter->at(pTBinD)) && (massD < sidebandLeftInner->at(pTBinD))) ||
          ((massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD)))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebandsMCRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsMCRec"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusHadrons, processMcRec, "Process MC Reco mode", true);

  /// D-Hadron correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::DplusHadronPair const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      // reject entries outside pT ranges of interest
      if (o2::analysis::findBin(binsPtCorrelations, ptD) < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }

      registry.fill(HIST("hCorrel2DVsPtMCGen"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
      registry.fill(HIST("hCorrel2DPtIntMCGen"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntMCGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMCGen"), deltaPhi);
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusHadrons, processMcGen, "Process MC Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationDplusHadrons>(cfgc)};
}
