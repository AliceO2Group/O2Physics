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
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::aod::hf_correlation_dplushadron;
using namespace o2::analysis::hf_cuts_dplus_topikpi;
using namespace o2::constants::math;

namespace o2::aod
{
using DplusHadronPairFull = soa::Join<aod::DplusHadronPair, aod::DplusHadronRecoInfo>;
} // namespace o2::aod

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiD, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PI / 2.);
}

///
/// Returns phi of candidate/particle evaluated from x and y components of segment connecting primary and secondary vertices
///
double evaluatePhiByVertex(double xVertex1, double xVertex2, double yVertex1, double yVertex2)
{
  return RecoDecay::Phi(xVertex2 - xVertex1, yVertex2 - yVertex1);
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

// definition of vectors for standard ptbin and invariant mass configurables
//   int arr[] = { 10, 20, 30 };
//   vector<int> vect(arr, arr + n);  for (int x : vect) cout << x << " "; 10,20, 30
const int npTBinsCorrelations = 8;
const double pTBinsCorrelations[npTBinsCorrelations + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
auto pTBinsCorrelations_v = std::vector<double>{pTBinsCorrelations, pTBinsCorrelations + npTBinsCorrelations + 1};
const double signalRegionInnerDefault[npTBinsCorrelations] = {1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490};
const double signalRegionOuterDefault[npTBinsCorrelations] = {1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890};
const double sidebandLeftInnerDefault[npTBinsCorrelations] = {1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690};
const double sidebandLeftOuterDefault[npTBinsCorrelations] = {1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250};
const double sidebandRightInnerDefault[npTBinsCorrelations] = {1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130};
const double sidebandRightOuterDefault[npTBinsCorrelations] = {1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690};
auto signalRegionInner_v = std::vector<double>{signalRegionInnerDefault, signalRegionInnerDefault + npTBinsCorrelations};
auto signalRegionOuter_v = std::vector<double>{signalRegionOuterDefault, signalRegionOuterDefault + npTBinsCorrelations};
auto sidebandLeftInner_v = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + npTBinsCorrelations};
auto sidebandLeftOuter_v = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + npTBinsCorrelations};
auto sidebandRightInner_v = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + npTBinsCorrelations};
auto sidebandRightOuter_v = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + npTBinsCorrelations};
const int npTBinsEfficiency = o2::analysis::hf_cuts_dplus_topikpi::npTBins;
const double efficiencyDmesonDefault[npTBinsEfficiency] = {};
const double efficiencyHadronDefault[npTBinsEfficiency] = {};
auto efficiencyDmeson_v = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + npTBinsEfficiency};
auto efficiencyHadron_v = std::vector<double>{efficiencyHadronDefault, efficiencyHadronDefault + npTBinsEfficiency};

/// Dplus-Hadron correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfTaskCorrelationDplusHadrons {

  HistogramRegistry registry{
    "registry",
    // NOTE: use hMassDplus (from correlator task) for normalisation, and hMass2DCorrelationPairs for 2D-sideband-subtraction purposes
    {
      {"hDeltaEtaPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hCorrel2DVsPtSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {1000, 0., 100.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
      {"hDeltaEtaPtIntSidebands", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hCorrel2DVsPtSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {1000, 0., 100.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
      {"hDeltaEtaPtIntSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hDeltaEtaPtIntSidebandsMCRec", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hCorrel2DPtIntSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hCorrel2DVsPtSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {1000, 0., 100.}}}},
      {"hDeltaPhiPtIntSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DPtIntSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hCorrel2DVsPtSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {1000, 0., 100.}}}},
      {"hDeltaEtaPtIntMCGen", stringMCParticles + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntMCGen", stringMCParticles + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DPtIntMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hCorrel2DVsPtMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + stringPtD + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {120, -6., 6.}, {10, 0., 10.}, {1000, 0., 100.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
    }};
  // pT ranges for correlation plots: the default values are those embedded in hf_cuts_dplus_topikpi (i.e. the mass pT bins), but can be redefined via json files
  Configurable<std::vector<double>> binsCorrelations{"ptBinsForCorrelations", std::vector<double>{pTBinsCorrelations_v}, "pT bin limits for correlation plots"};
  // pT bins for effiencies: same as above
  Configurable<std::vector<double>> binsEfficiency{"ptBinsForEfficiency", std::vector<double>{o2::analysis::hf_cuts_dplus_topikpi::pTBins_v}, "pT bin limits for efficiency"};
  // signal and sideband region edges, to be defined via json file (initialised to empty)
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", std::vector<double>{signalRegionInner_v}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", std::vector<double>{signalRegionOuter_v}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{sidebandLeftInner_v}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{sidebandLeftOuter_v}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{sidebandRightInner_v}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{sidebandRightOuter_v}, "Outer values of right sideband vs pT"};
  Configurable<std::vector<double>> efficiencyDmeson{"efficiencyDmeson", std::vector<double>{efficiencyDmeson_v}, "Efficiency values for D meson specie under study"};
  Configurable<std::vector<double>> efficiencyHadron{"efficiencyHadron", std::vector<double>{efficiencyHadron_v}, "Efficiency values for Hadron specie under study"};
  Configurable<int> flagApplyEfficiency{"efficiencyFlagDplus", 1, "Flag for applying efficiency weights"};

  void init(o2::framework::InitContext&)
  {
    // redefinition of pT axes for THnSparse holding correlation entries
    int nBinspTaxis = binsCorrelations->size() - 1;
    const double* valuespTaxis = binsCorrelations->data();

    for (int i = 1; i < 2; i++) {
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->GetAxis(i + 1)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->GetAxis(i + 1)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRec"))->GetAxis(i + 1)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRec"))->GetAxis(i + 1)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtMCGen"))->GetAxis(i + 1)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hCorrel2DVsPtMCGen"))->Sumw2();
    }
  }

  void processData(aod::DplusHadronPairFull const& pairEntries)
  {
    for (auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      double massD = pairEntry.mD();
      int pTBinD = o2::analysis::findBin(binsCorrelations, ptD);
      double efficiencyWeight = 1.;
      // if (flagApplyEfficiency) {
      //  efficiencyWeight = 1. / (efficiencyDmeson->at(o2::analysis::findBin(binsEfficiency, ptD)) * efficiencyHadron->at(o2::analysis::findBin(binsEfficiency, ptHadron)));
      // }

      // reject entries outside pT ranges of interest
      if (pTBinD == -1) { // at least one particle outside accepted pT range
        continue;
      }

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }

      if ((massD > sidebandLeftInner->at(pTBinD) && massD < sidebandLeftOuter->at(pTBinD)) ||
          (massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
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
    for (auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      double massD = pairEntry.mD();
      int pTBinD = o2::analysis::findBin(binsCorrelations, ptD);
      double efficiencyWeight = 1.;
      // if (flagApplyEfficiency) {
      //  efficiencyWeight = 1. / (efficiencyDmeson->at(o2::analysis::findBin(binsEfficiency, ptD)) * efficiencyHadron->at(o2::analysis::findBin(binsEfficiency, ptHadron)));
      // }

      // reject entries outside pT ranges of interest
      if (pTBinD == -1) { // at least one particle outside accepted pT range
        continue;
      }

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionMCRec"), deltaPhi, efficiencyWeight);
      }

      if (((massD > sidebandLeftInner->at(pTBinD)) && (massD < sidebandLeftOuter->at(pTBinD))) ||
          ((massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD)))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebandsMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsMCRec"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusHadrons, processMcRec, "Process MC Reco mode", true);

  /// D-Hadron correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal) - Ok for both USL and LS analyses
  void processMcGen(aod::DplusHadronPair const& pairEntries)
  {
    for (auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      // reject entries outside pT ranges of interest
      if (o2::analysis::findBin(binsCorrelations, ptD) == -1) {
        continue;
      }

      registry.fill(HIST("hCorrel2DVsPtMCGen"), deltaPhi, deltaEta, ptD, ptHadron);
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
