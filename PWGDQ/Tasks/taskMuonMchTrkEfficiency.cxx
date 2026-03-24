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

/// \file taskMuonMchTrkEfficiency.cxx
/// \brief task to read the table created by muonMchEfficiency.cxx and build the plots for evaluating the muon tracking efficiency
///
/// @param muon MCH tracking efficiency table
/// Struct for reading the table created by muonMchEfficiency.cxx and filling the histos needed to compute the
///    muon tracking efficiency in the MCH detector
///
/// \author Zaida Conesa del Valle <zaida.conesa.del.valle@cern.ch>
///

#include "PWGDQ/DataModel/MchTrkEffTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// default histogram output binning
namespace muon_trk_eff_bins
{
static constexpr int nBinsPt = 24;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double binsPt[nBinsPt + 1] = {
  0., 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0,
  2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0,
  15.0, 18.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0};
auto vecBinsPt = std::vector<double>{binsPt, binsPt + nBinsPt + 1};

// row labels
static const std::vector<std::string> labelsPtTrack{
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

} // namespace muon_trk_eff_bins

struct taskMuonMchTrkEfficiency {

  Configurable<int> muonSelType{"muonSelType", 4, "Selected muon track type"}; /// Muon track type to be selected if value >=0 (no selection by default)

  Configurable<double> ptMuonMin{"ptMin", 0., "Lower bound of pT"};       /// Muon minimum pt to be studied
  Configurable<double> etaMuonMin{"etaMin", 2.5, "Lower bound of |eta|"}; /// Muon minimum |eta| to be studied
  Configurable<double> etaMuonMax{"etaMax", 4.0, "Upper bound of |eta|"}; /// Muon maximum |eta| to be studied

  Configurable<std::vector<double>> binsMuonPt{"binsPt", std::vector<double>{muon_trk_eff_bins::vecBinsPt}, "pT bin limits"}; /// Pt intervals for the histograms
  Configurable<int> nEtaBins{"nEtaBins", 12, "Number of Eta bins"};                                                           /// Number of eta bins for output histograms
  Configurable<int> nPhiBins{"nPhiBins", 6, "Number of Phi bins"};                                                            /// Number of phi bins for output histograms

  using muonFull = soa::Join<aod::MchTrkEffBase, aod::MchTrkEffGen>;

  /// Histogram registry: an object to hold your histograms
  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  ///  Initialize: configure, create specifics
  void init(o2::framework::InitContext&)
  {
    LOGF(debug, "Initialization");

    // define axes to be used
    // and the correspondent histograms
    auto vbins = (std::vector<double>)binsMuonPt;
    const AxisSpec axisEta{nEtaBins, etaMuonMin, etaMuonMax, "|#eta|"};
    const AxisSpec axisPt{vbins, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPhi{nPhiBins, -3.14, 3.14, "#varphi"};

    const AxisSpec axisEtaGen{nEtaBins, etaMuonMin, etaMuonMax, "|#eta| Gen"};
    const AxisSpec axisPtGen{vbins, "#it{p}_{T} (GeV/#it{c}) Gen"};
    const AxisSpec axisPhiGen{nPhiBins, -3.14, 3.14, "#varphi Gen"};
    const AxisSpec axismchBitmap{1031, -0.5, 1030.5, "mchBitmap"};

    const AxisSpec axisNhits{16, -0.5, 15.5, ""};
    // Labels for the chambers hit per station, numbering starting from 1
    // i.e. (Nij, N0j, Ni0) correspond to hit on i-j, hit on j not on i, hit on i not on j
    const char* elabels[15] = {"N12", "N10", "N02", "N34", "N30", "N04", "N56", "N50", "N06", "N78", "N70", "N08", "N910", "N90", "N010"};

    HistogramConfigSpec defaultNhitsEtaPtPhi({HistType::kTHnF, {{axisNhits}, {axisEta}, {axisPt}, {axisPhi}}});

    // define histograms to be added
    LOGF(debug, " Creating histograms");

    registry.add("hmchBitmap", "hmchBitmap", {HistType::kTH1F, {axismchBitmap}});

    registry.add("hPtRecPtGen", "hPtRecPtGen", {HistType::kTH2F, {{axisPt}, {axisPtGen}}});
    registry.add("hEtaRecEtaGen", "hEtaRecEtaGen", {HistType::kTH2F, {{axisEta}, {axisEtaGen}}});
    registry.add("hPhiRecPhiGen", "hPhiRecPhiGen", {HistType::kTH2F, {{axisPhi}, {axisPhiGen}}});

    registry.add("hHitsEtaPtPhi", "hHitsEtaPtPhi", defaultNhitsEtaPtPhi, false);
    auto hHitsEtaPtPhi = registry.get<THn>(HIST("hHitsEtaPtPhi"));
    for (int i = 0; i < 15; i++)
      hHitsEtaPtPhi->GetAxis(0)->SetBinLabel(i + 1, elabels[i]);

  } //! end of Initialize: configure, create specifics

  /// check whether a given chamber has hits
  bool ischamberhit(uint16_t map, int ich)
  { // i = 0..9
    LOGF(debug, " map %i --> %i", map, (map >> ich) & 1);
    return (map >> ich) & 1;
  }

  /// Kinematic selection
  bool IsInKinematics(double eta, double pt)
  {
    bool isSelected = true;

    if (pt < ptMuonMin) {
      isSelected = false;
    }
    if ((eta < etaMuonMin) || (eta > etaMuonMax)) {
      isSelected = false;
    }
    return isSelected;
  }

  /// Method to define and apply the weight required to optimise the generated distributions
  ///  TO BE UPDATED with realistic distributions and probably using input TF1 as weights
  void FillHistosWeight(double eta, double pt, double phi, uint16_t map, double /*etaGen*/, double /*ptGen*/, double /*phiGen*/)
  {
    double weighteta = 1, weightpt = 1, weightphi = 1; // default weight set to unity for now: no effect
    double etaw = eta * weighteta;
    double ptw = pt * weightpt;
    double phiw = phi * weightphi;
    FillHistos(etaw, ptw, phiw, map);
  }

  /// Filling histograms from generated & reconstructed information
  void FillHistosMC(double eta, double pt, double phi, uint16_t map, double etaGen, double ptGen, double phiGen)
  {
    registry.fill(HIST("hPtRecPtGen"), pt, ptGen);
    registry.fill(HIST("hEtaRecEtaGen"), eta, etaGen);
    registry.fill(HIST("hPhiRecPhiGen"), phi, phiGen);
    FillHistosWeight(eta, pt, phi, map, etaGen, ptGen, phiGen);
  }

  /// Filling histograms from reconstructed quantities
  void FillHistos(double eta, double pt, double phi, uint16_t map)
  {
    registry.fill(HIST("hmchBitmap"), map);

    bool iN[10];
    for (int i = 0; i < 10; i++) {
      bool ishit = ischamberhit(map, i);
      iN[i] = false;
      if (ishit)
        iN[i] = true;
    }
    if (iN[0] && iN[1])
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(0., eta, pt, phi);
    if (iN[0] && (!iN[1]))
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(1., eta, pt, phi);
    if ((!iN[0]) && iN[1])
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(2., eta, pt, phi);
    if (iN[2] && iN[3])
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(3., eta, pt, phi);
    if (iN[2] && (!iN[3]))
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(4., eta, pt, phi);
    if ((!iN[2]) && iN[3])
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(5., eta, pt, phi);
    if (iN[4] && iN[5])
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(6., eta, pt, phi);
    if (iN[4] && (!iN[5]))
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(7., eta, pt, phi);
    if ((!iN[4]) && iN[5])
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(8., eta, pt, phi);
    if (iN[6] && iN[7])
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(9., eta, pt, phi);
    if (iN[6] && (!iN[7]))
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(10., eta, pt, phi);
    if ((!iN[6]) && iN[7])
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(11., eta, pt, phi);
    if (iN[8] && iN[9])
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(12., eta, pt, phi);
    if (iN[8] && (!iN[9]))
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(13., eta, pt, phi);
    if ((!iN[8]) && iN[9])
      registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(14., eta, pt, phi);
  }

  //! process function for reconstructed muon information
  void processReco(aod::MchTrkEffBase const& mchtrkeffbases)
  {
    for (auto& mchtrkeffbase : mchtrkeffbases) {
      bool isSel = true;
      if (!IsInKinematics(mchtrkeffbase.eta(), mchtrkeffbase.pt()))
        isSel = false;
      if (muonSelType >= 0) {
        if (mchtrkeffbase.trackType() != muonSelType)
          isSel = false;
      }
      if (isSel)
        FillHistos(mchtrkeffbase.eta(), mchtrkeffbase.pt(), mchtrkeffbase.phi(), mchtrkeffbase.mchBitMap());
    }
  }
  PROCESS_SWITCH(taskMuonMchTrkEfficiency, processReco, "process reconstructed information", true);

  //! process function for simulated muon information
  void processSim(muonFull const& mchtrkeffbases)
  {
    for (auto& mchtrkeffbase : mchtrkeffbases) {
      if (IsInKinematics(mchtrkeffbase.eta(), mchtrkeffbase.pt()))
        FillHistosMC(mchtrkeffbase.eta(), mchtrkeffbase.pt(), mchtrkeffbase.phi(), mchtrkeffbase.mchBitMap(), mchtrkeffbase.etaGen(), mchtrkeffbase.ptGen(), mchtrkeffbase.phiGen());
    }
  }
  PROCESS_SWITCH(taskMuonMchTrkEfficiency, processSim, "process simulated information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskMuonMchTrkEfficiency>(cfgc)};
}
