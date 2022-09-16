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
///
/// \author Samrangy Sadhu <samrangy.sadhu@cern.ch>, INFN Bari

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_correlation_d0hadron;
using namespace o2::analysis::hf_cuts_d0_topik;
using namespace o2::constants::math;

namespace o2::aod
{
using DHadronPairFull = soa::Join<aod::DHadronPair, aod::DHadronRecoInfo>;
} // namespace o2::aod


///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiD, double phiDbar)
{
  return RecoDecay::constrainAngle(phiDbar - phiD, -o2::constants::math::PI / 2.);
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

//definition of vectors for standard ptbin and invariant mass configurables
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
const int npTBinsEfficiency = o2::analysis::hf_cuts_d0_topik::npTBins;
const double efficiencyDmesonDefault[npTBinsEfficiency] = {};
auto efficiencyDmeson_v = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + npTBinsEfficiency};

struct taskCorrelationD0Hadrons {

  HistogramRegistry registry{
    "registry",
    {
         {"hDeltaEtaPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
         {"hDeltaPhiPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hCorrel2DVsPtSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
      {"hDeltaEtaPtIntSidebands", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hCorrel2DVsPtSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
      {"hCorrel2DVsPtSignalMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}},
      {"hCorrel2DVsPtBkgMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}},
      {"hCorrel2DVsPtSigSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}},
      {"hCorrel2DPtIntSigSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hDeltaEtaPtIntSigSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntSigSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DVsPtSigSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}},
      {"hCorrel2DPtIntSigSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hDeltaEtaPtIntSigSidebandsMCRec", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntSigSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      
      {"hCorrel2DVsPtRefSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}},
      {"hCorrel2DPtIntRefSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hDeltaEtaPtIntRefSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntRefSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DVsPtRefSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}},
      {"hCorrel2DPtIntRefSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hDeltaEtaPtIntRefSidebandsMCRec", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntRefSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      
      {"hCorrel2DVsPtBkgSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}},
      {"hCorrel2DPtIntBkgSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hDeltaEtaPtIntBkgSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntBkgSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DVsPtBkgSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}},
      {"hCorrel2DPtIntBkgSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hDeltaEtaPtIntBkgSidebandsMCRec", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntBkgSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}},
      {"hCorrel2DVsPtMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + stringPtD + "entries", {HistType::kTHnSparseD, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {40, -2., 2.}, {10, 0., 10.}, {11, 0., 11.}}}}, // note: axes 3 and 4 (the pT) are updated in the init(),
      {"hCorrel2DPtIntMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}, {200, -10., 10.}}}},
      {"hDeltaEtaPtIntMCGen", stringMCParticles + stringDeltaEta + "entries", {HistType::kTH1F, {{200, -10., 10.}}}},
      {"hDeltaPhiPtIntMCGen", stringMCParticles + stringDeltaPhi + "entries", {HistType::kTH1F, {{64, -o2::constants::math::PI / 2., 3. * o2::constants::math::PI / 2.}}}}

     }};
 
  //pT ranges for correlation plots: the default values are those embedded in hf_cuts_d0_topik (i.e. the mass pT bins), but can be redefined via json files
  Configurable<std::vector<double>> binsCorrelations{"ptBinsForCorrelations", std::vector<double>{pTBinsCorrelations_v}, "pT bin limits for correlation plots"};
  //pT bins for effiencies: same as above
  Configurable<std::vector<double>> binsEfficiency{"ptBinsForEfficiency", std::vector<double>{o2::analysis::hf_cuts_d0_topik::pTBins_v}, "pT bin limits for efficiency"};
  //signal and sideband region edges, to be defined via json file (initialised to empty)
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", std::vector<double>{signalRegionInner_v}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", std::vector<double>{signalRegionOuter_v}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{sidebandLeftInner_v}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{sidebandLeftOuter_v}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{sidebandRightInner_v}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{sidebandRightOuter_v}, "Outer values of right sideband vs pT"};
  Configurable<std::vector<double>> efficiencyDmeson{"efficiencyDmeson", std::vector<double>{efficiencyDmeson_v}, "Efficiency values for D meson specie under study"};
  Configurable<int> flagApplyEfficiency{"efficiencyFlagD", 1, "Flag for applying efficiency weights"};

  void init(o2::framework::InitContext&)
  {

   int nBinspTaxis = binsCorrelations->size() - 1;
    const double* valuespTaxis = binsCorrelations->data(); 
     registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSigSignalRegionMCRec"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSigSidebandsMCRec"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSigSignalRegionMCRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSigSidebandsMCRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalMCRec"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalMCRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtBkgMCRec"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtBkgMCRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMCGen"))->GetAxis(2)->Set(nBinspTaxis, valuespTaxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMCGen"))->Sumw2();
  }

  /// D-h correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  /// Works on both USL and LS analyses pair tables
  void processData(aod::DHadronPairFull const& pairEntries)
  {
   for (auto& pairEntry : pairEntries) {
   //define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      double massD = pairEntry.mD();
      double massDbar = pairEntry.mDbar();
      int SignalStatus = pairEntry.signalStatus();
      int effBinD = o2::analysis::findBin(binsEfficiency, ptD);
      int pTBinD = o2::analysis::findBin(binsCorrelations, ptD);

// reject entries outside pT ranges of interest
      if (pTBinD < 0 || effBinD < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }

      double efficiencyWeight = 1.;
      if (flagApplyEfficiency) {
        efficiencyWeight = 1. / (efficiencyDmeson->at(o2::analysis::findBin(binsEfficiency, ptD))); // ***** track efficiency to be implemented *****
      }

      

      //reject entries outside pT ranges of interest
      if (pTBinD == -1) { //at least one particle outside accepted pT range
        continue;
      }
// check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if ((massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) && ((SignalStatus =1) || (SignalStatus =3))) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }
      
      if ((massDbar > signalRegionInner->at(pTBinD) && massDbar < signalRegionOuter->at(pTBinD)) && (SignalStatus >=2)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }
      
      if (((massD > sidebandLeftOuter->at(pTBinD) && massD < sidebandLeftInner->at(pTBinD)) ||
          (massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD))) && ((SignalStatus =1) || (SignalStatus =3))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
      
      if (((massDbar > sidebandLeftOuter->at(pTBinD) && massDbar < sidebandLeftInner->at(pTBinD)) ||
          (massDbar > sidebandRightInner->at(pTBinD) && massDbar < sidebandRightOuter->at(pTBinD))) && (SignalStatus >=2)) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
      
   } 
  }
  PROCESS_SWITCH(taskCorrelationD0Hadrons, processData, "Process data", false);
  
  void processMcRec(aod::DHadronPairFull const& pairEntries)
  {
    for (auto& pairEntry : pairEntries) {
      //define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      double massD = pairEntry.mD();
      double massDbar = pairEntry.mDbar();
      int SignalStatus = pairEntry.signalStatus();

      int pTBinD = o2::analysis::findBin(binsCorrelations, ptD);


      double efficiencyWeight = 1.;
      if (flagApplyEfficiency) {
        efficiencyWeight = 1. / (efficiencyDmeson->at(o2::analysis::findBin(binsEfficiency, ptD)));
      }
      
      // fill correlation plots for signal/bagkground correlations
      if (pairEntry.signalStatus()) {
        registry.fill(HIST("hCorrel2DVsPtSignalMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        
      } else {
        registry.fill(HIST("hCorrel2DVsPtBkgMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
      }
      
      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      
      // ---------------------- Fill plots for signal case, D0 ->1, D0bar ->8 ---------------------------------------------
      if ((massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) && (SignalStatus =1)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSigSignalRegionMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSigSignalRegionMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSigSignalRegionMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSigSignalRegionMCRec"), deltaPhi, efficiencyWeight);
      }
      
      if ((massDbar > signalRegionInner->at(pTBinD) && massDbar < signalRegionOuter->at(pTBinD)) && (SignalStatus =8)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSigSignalRegionMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSigSignalRegionMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSigSignalRegionMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSigSignalRegionMCRec"), deltaPhi, efficiencyWeight);
      }

      if ((((massD > sidebandLeftOuter->at(pTBinD)) && (massD < sidebandLeftInner->at(pTBinD))) ||
          ((massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD)))) && (SignalStatus =1)) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSigSidebandsMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSigSidebandsMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSigSidebandsMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSigSidebandsMCRec"), deltaPhi, efficiencyWeight);
      }
      
      if ((((massDbar > sidebandLeftOuter->at(pTBinD)) && (massDbar < sidebandLeftInner->at(pTBinD))) ||
          ((massDbar > sidebandRightInner->at(pTBinD) && massDbar < sidebandRightOuter->at(pTBinD)))) && (SignalStatus =8)) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSigSidebandsMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSigSidebandsMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSigSidebandsMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSigSidebandsMCRec"), deltaPhi, efficiencyWeight);
      }
      
      // ---------------------- Fill plots for reflection case, D0 ->2, D0bar ->16 ---------------------------------------------
      if ((massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) && (SignalStatus =2)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtRefSignalRegionMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntRefSignalRegionMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntRefSignalRegionMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntRefSignalRegionMCRec"), deltaPhi, efficiencyWeight);
      }
      
      if ((massDbar > signalRegionInner->at(pTBinD) && massDbar < signalRegionOuter->at(pTBinD)) && (SignalStatus =16)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtRefSignalRegionMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntRefSignalRegionMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntRefSignalRegionMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntRefSignalRegionMCRec"), deltaPhi, efficiencyWeight);
      }

      if ((((massD > sidebandLeftOuter->at(pTBinD)) && (massD < sidebandLeftInner->at(pTBinD))) ||
          ((massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD)))) && (SignalStatus =2)) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtRefSidebandsMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntRefSidebandsMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntRefSidebandsMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntRefSidebandsMCRec"), deltaPhi, efficiencyWeight);
      }
      
      if ((((massDbar > sidebandLeftOuter->at(pTBinD)) && (massDbar < sidebandLeftInner->at(pTBinD))) ||
          ((massDbar > sidebandRightInner->at(pTBinD) && massDbar < sidebandRightOuter->at(pTBinD)))) && (SignalStatus =16)) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtRefSidebandsMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntRefSidebandsMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntRefSidebandsMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntRefSidebandsMCRec"), deltaPhi, efficiencyWeight);
      }
      
      
      // ---------------------- Fill plots for background case, D0 ->4, D0bar ->32 ---------------------------------------------
      if ((massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) && (SignalStatus =4)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtBkgSignalRegionMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntBkgSignalRegionMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntBkgSignalRegionMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntBkgSignalRegionMCRec"), deltaPhi, efficiencyWeight);
      }
      
      if ((massDbar > signalRegionInner->at(pTBinD) && massDbar < signalRegionOuter->at(pTBinD)) && (SignalStatus =32)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtBkgSignalRegionMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntBkgSignalRegionMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntBkgSignalRegionMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntBkgSignalRegionMCRec"), deltaPhi, efficiencyWeight);
      }

      if ((((massD > sidebandLeftOuter->at(pTBinD)) && (massD < sidebandLeftInner->at(pTBinD))) ||
          ((massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD)))) && (SignalStatus =4)) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtBkgSidebandsMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntBkgSidebandsMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntBkgSidebandsMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntBkgSidebandsMCRec"), deltaPhi, efficiencyWeight);
      }
      
      if ((((massDbar > sidebandLeftOuter->at(pTBinD)) && (massDbar < sidebandLeftInner->at(pTBinD))) ||
          ((massDbar > sidebandRightInner->at(pTBinD) && massDbar < sidebandRightOuter->at(pTBinD)))) && (SignalStatus =32)) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtBkgSidebandsMCRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntBkgSidebandsMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntBkgSidebandsMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntBkgSidebandsMCRec"), deltaPhi, efficiencyWeight);
      }
      
     }
    }
   PROCESS_SWITCH(taskCorrelationD0Hadrons, processMcRec, "Process MC Reco mode", true); 
   
   /// D-Hadron correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::DHadronPairFull const& pairEntries)
  {
    for (auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      // reject entries outside pT ranges of interest
      if (o2::analysis::findBin(binsCorrelations, ptD) < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }

      registry.fill(HIST("hCorrel2DVsPtMCGen"), deltaPhi, deltaEta, ptD, ptHadron);
      registry.fill(HIST("hCorrel2DPtIntMCGen"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntMCGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMCGen"), deltaPhi);
    } // end loop
  }
  PROCESS_SWITCH(taskCorrelationD0Hadrons, processMcGen, "Process MC Gen mode", false);
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<taskCorrelationD0Hadrons>(cfgc)};
}
