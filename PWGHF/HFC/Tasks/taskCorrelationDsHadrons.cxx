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

/// \file taskCorrelationDsHadrons.cxx
/// \author Grazia Luparello <Grazia.Luparello@cern.ch>
/// \author Samuele Cattaruzzi <Samuele.Cattaruzzi@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::aod::hf_correlation_ds_hadron;
using namespace o2::analysis::hf_cuts_ds_to_k_k_pi;
using namespace o2::constants::math;

namespace o2::aod
{
using DsHadronPairFull = soa::Join<aod::DsHadronPair, aod::DsHadronRecoInfo>;
}

/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
double getDeltaPhi(double phiD, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PIHalf);
}

/// Returns phi of candidate/particle evaluated from x and y components of segment connecting primary and secondary vertices
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

const int nBinsPtCorrelations = 8;
const double pTBinsCorrelations[nBinsPtCorrelations + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
auto vecBinsPtCorrelations = std::vector<double>{pTBinsCorrelations, pTBinsCorrelations + nBinsPtCorrelations + 1};
const double signalRegionInnerDefault[nBinsPtCorrelations] = {1.9440, 1.9440, 1.9440, 1.9440, 1.9440, 1.9440, 1.9440, 1.9440};
const double signalRegionOuterDefault[nBinsPtCorrelations] = {1.9920, 1.9920, 1.9920, 1.9920, 1.9920, 1.9920, 1.9920, 1.9920};
const double sidebandLeftOuterDefault[nBinsPtCorrelations] = {1.9040, 1.9040, 1.9040, 1.9040, 1.9040, 1.9040, 1.9040, 1.9040};
const double sidebandLeftInnerDefault[nBinsPtCorrelations] = {1.9360, 1.9360, 1.9360, 1.9360, 1.9360, 1.9360, 1.9360, 1.9360};
const double sidebandRightInnerDefault[nBinsPtCorrelations] = {2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000};
const double sidebandRightOuterDefault[nBinsPtCorrelations] = {2.0320, 2.0320, 2.0320, 2.0320, 2.0320, 2.0320, 2.0320, 2.0320};
auto vecSignalRegionInner = std::vector<double>{signalRegionInnerDefault, signalRegionInnerDefault + nBinsPtCorrelations};
auto vecSignalRegionOuter = std::vector<double>{signalRegionOuterDefault, signalRegionOuterDefault + nBinsPtCorrelations};
auto vecSidebandLeftInner = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + nBinsPtCorrelations};
auto vecSidebandLeftOuter = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + nBinsPtCorrelations};
auto vecSidebandRightInner = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + nBinsPtCorrelations};
auto vecSidebandRightOuter = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + nBinsPtCorrelations};
const int nBinsPtEfficiency = o2::analysis::hf_cuts_ds_to_k_k_pi::nBinsPt;
const double efficiencyDmesonDefault[nBinsPtEfficiency] = {};
auto vecEfficiencyDmeson = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + nBinsPtEfficiency};

/// Ds-Hadron correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfTaskCorrelationDsHadrons {
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying efficiency weights"};
  // pT ranges for correlation plots: the default values are those embedded in hf_cuts_ds_to_k_k_pi (i.e. the mass pT bins), but can be redefined via json files
  Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{vecBinsPtCorrelations}, "pT bin limits for correlation plots"};
  Configurable<std::vector<double>> binsPtEfficiency{"binsPtEfficiency", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for efficiency"};
  // signal and sideband region edges, to be defined via json file (initialised to empty)
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", std::vector<double>{vecSignalRegionInner}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", std::vector<double>{vecSignalRegionOuter}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{vecSidebandLeftInner}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{vecSidebandLeftOuter}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{vecSidebandRightInner}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{vecSidebandRightOuter}, "Outer values of right sideband vs pT"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", std::vector<double>{vecEfficiencyDmeson}, "Efficiency values for D meson specie under study"};

  HistogramRegistry registry{
    "registry",
    {
      {"hDeltaEtaPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -4., 4.}}}},
      {"hDeltaPhiPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}}}},
      {"hCorrel2DPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {200, -4., 4.}}}},
      {"hCorrel2DVsPtSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
      {"hDeltaEtaPtIntSidebands", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -4., 4.}}}},
      {"hDeltaPhiPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}}}},
      {"hCorrel2DPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {200, -4., 4.}}}},
      {"hCorrel2DVsPtSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
      {"hDeltaEtaPtIntSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -4., 4.}}}},
      {"hDeltaPhiPtIntSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}}}},
      {"hDeltaEtaPtIntSidebandsMCRec", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -4., 4.}}}},
      {"hCorrel2DPtIntSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {200, -4., 4.}}}},
      {"hCorrel2DVsPtSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}},
      {"hCorrel2DVsPtSignalMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}},
      {"hCorrel2DVsPtBkgMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}},
      {"hDeltaPhiPtIntSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}}}},
      {"hCorrel2DPtIntSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {200, -4., 4.}}}},
      {"hCorrel2DVsPtSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}},
      {"hDeltaEtaPtIntMCGen", stringMCParticles + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -4., 4.}}}},
      {"hDeltaPhiPtIntMCGen", stringMCParticles + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}}}},
      {"hCorrel2DPtIntMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {200, -4., 4.}}}},
      {"hCorrel2DVsPtMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + stringPtD + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
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

  void processData(aod::DsHadronPairFull const& pairEntries)
  {
    for (auto const& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      double massD = pairEntry.mD();
      int effBinD = o2::analysis::findBin(binsPtEfficiency, ptD);
      int pTBinD = o2::analysis::findBin(binsPtCorrelations, ptD);
      // reject entries outside pT ranges of interest
      if (pTBinD < 0 || effBinD < 0) {
        continue;
      }
      if (ptHadron > 10.0) { // all Hadrons with pT > 10 are put in the 11th bin of the axis of the histograms filled with them
        ptHadron = 10.5;
      }
      double efficiencyWeight = 1.;
      double efficiencyHadron = 1.; // Note: To be implemented later on ??
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(effBinD) * efficiencyHadron);
      }

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }

      if ((massD > sidebandLeftOuter->at(pTBinD) && massD < sidebandLeftInner->at(pTBinD)) ||
          (massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processData, "Process data", true);

  /// D-Hadron correlation pair filling task, from pair tables - for MC reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRec(aod::DsHadronPairFull const& pairEntries)
  {
    for (auto const& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      double massD = pairEntry.mD();
      int effBinD = o2::analysis::findBin(binsPtEfficiency, ptD);
      int pTBinD = o2::analysis::findBin(binsPtCorrelations, ptD);
      if (pTBinD < 0 || effBinD < 0) {
        continue;
      }
      if (ptHadron > 10.0) { // all Hadrons with pT > 10 are put in the 11th bin of the axis of the histograms filled with them
        ptHadron = 10.5;
      }
      double efficiencyWeight = 1.;
      double efficiencyHadron = 1.; // Note: To be implemented later on
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(effBinD) * efficiencyHadron);
      }
      // fill correlation plots for signal/bagkground correlations
      if (pairEntry.signalStatus()) {
        registry.fill(HIST("hCorrel2DVsPtSignalMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
      } else {
        registry.fill(HIST("hCorrel2DVsPtBkgMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
      }
      // reject entries outside pT ranges of interest

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionMCRec"), deltaPhi, efficiencyWeight);
      }

      if (((massD > sidebandLeftOuter->at(pTBinD)) && (massD < sidebandLeftInner->at(pTBinD))) ||
          ((massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD)))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebandsMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsMCRec"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processMcRec, "Process MC Reco mode", false);

  /// D-Hadron correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::DsHadronPair const& pairEntries)
  {
    for (auto const& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      // reject entries outside pT ranges of interest
      if (o2::analysis::findBin(binsPtCorrelations, ptD) < 0) {
        continue;
      }
      if (ptHadron > 10.0) { // all Hadrons with pT > 10 are put in the 11th bin of the axis of the histograms filled with them
        ptHadron = 10.5;
      }

      registry.fill(HIST("hCorrel2DVsPtMCGen"), deltaPhi, deltaEta, ptD, ptHadron);
      registry.fill(HIST("hCorrel2DPtIntMCGen"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntMCGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMCGen"), deltaPhi);
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processMcGen, "Process MC Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationDsHadrons>(cfgc)};
}
