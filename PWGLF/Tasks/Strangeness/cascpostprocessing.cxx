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
///
/// \brief post processing for Cascade analysis
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \modified by Chiara De Martin (chiara.de.martin@cern.ch)
/// \since March 20, 2023
/// \modified by Roman Nepeivoda (roman.nepeivoda@cern.ch)
/// \since June 1, 2023

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/cascqaanalysis.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

// constants
const float ctauxiPDG = 4.91;     // from PDG
const float ctauomegaPDG = 2.461; // from PDG

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct cascpostprocessing {
  // Xi or Omega
  Configurable<bool> isXi{"isXi", 1, "Apply cuts for Xi identification"};
  Configurable<bool> hastof{"hastof", 0, "Apply nsigmaTOF to daughter tracks of cascade"};
  Configurable<float> ptthrtof{"ptthrtof", 2, "NsigmaTOF selection is only applied if the track has pt larger than ptthr"};

  // Selection criteria
  Configurable<float> minpt{"minpt", 0.0, "Min p_{T} (GeV/c)"};
  Configurable<float> rap{"rap", 0.5, "Rapidity"};
  Configurable<float> masswin{"masswin", 0.075, "Mass window limit"};
  Configurable<float> masswintpc{"masswintpc", 0.075, "Mass window limit for Nsigma TPC daughters"};
  Configurable<float> rejcomp{"rejcomp", 0.008, "Competing Cascade rejection"};
  Configurable<float> lambdamasswin{"lambdamasswin", 0.008, "V0 Mass window limit"};
  Configurable<float> v0radiusMin{"v0radiusMin", 1.2, "V0 Minimum Radius"};
  Configurable<float> v0radiusMax{"v0radiusMax", 1000, "V0 Maximum Radius"};
  Configurable<float> cascradiusMin{"cascradiusMin", 0.6, "Casc Minimum Radius"};
  Configurable<float> cascradiusMax{"cascradiusMax", 1000, "Casc Maximum Radius"};
  Configurable<double> casccospa{"casccospa", 0.97, "Casc CosPA"};
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.03, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.03, "DCA Pos To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", 0.03, "DCA Bach To PV"};
  Configurable<float> dcacascdau{"dcacascdau", 1.3, "DCA Casc Daughters"};
  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};
  Configurable<float> dcav0topv{"dcav0topv", 0.06, "DCA V0 To PV"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};
  Configurable<float> proplifetime{"proplifetime", 6, "ctau/<ctau>"};
  Configurable<float> nsigmatpcPi{"nsigmatpcPi", 6, "N sigma TPC Pion"};
  Configurable<float> nsigmatpcPr{"nsigmatpcPr", 6, "N sigma TPC Proton"};
  Configurable<float> nsigmatpcKa{"nsigmatpcKa", 6, "N sigma TPC Kaon"};
  Configurable<float> nsigmatofPi{"nsigmatofPi", 6, "N sigma TOF Pion"};
  Configurable<float> nsigmatofPr{"nsigmatofPr", 6, "N sigma TOF Proton"};
  Configurable<float> nsigmatofKa{"nsigmatofKa", 6, "N sigma TOF Kaon"};
  Configurable<int> mintpccrrows{"mintpccrrows", 50, "min N TPC crossed rows"};
  Configurable<int> dooobrej{"dooobrej", 0, "OOB rejection: 0 no selection, 1 = ITS||TOF, 2 = TOF only for pT > ptthrtof"};

  Configurable<bool> isSelectBachBaryon{"isSelectBachBaryon", 0, "Bachelor-baryon cascade selection"};
  Configurable<float> bachBaryonCosPA{"bachBaryonCosPA", 0.9999, "Bachelor baryon CosPA"};
  Configurable<float> bachBaryonDCAxyToPV{"bachBaryonDCAxyToPV", 0.05, "DCA bachelor baryon to PV"};

  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};

  Configurable<int> evSelFlag{"evSelFlag", 2, "1 - INEL; 2 - INEL>0; 3 - INEL>1"};

  HistogramRegistry registry{"registryts"};

  // Necessary for particle charges
  Service<o2::framework::O2DatabasePDG> pdgDB;

  void init(InitContext const&)
  {

    AxisSpec ximassAxis = {200, 1.28f, 1.36f};
    AxisSpec omegamassAxis = {200, 1.59f, 1.75f};
    AxisSpec massAxis = ximassAxis;
    if (!isXi)
      massAxis = omegamassAxis;
    AxisSpec ptAxis = {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisTopoVar = {50, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisPID = {50, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis etaAxis{"etaAxis", {40, -2.0f, 2.0f}, "#eta"};

    ConfigurableAxis centFT0MAxis{"FT0M",
                                  {VARIABLE_WIDTH, 0., 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 101, 105.5},
                                  "FT0M (%)"};
    ConfigurableAxis centFV0AAxis{"FV0A",
                                  {VARIABLE_WIDTH, 0., 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 101, 105.5},
                                  "FV0A (%)"};

    ConfigurableAxis nChargedFT0MGenAxis{"nChargedFT0MGenAxis", {300, 0, 300}, "N_{FT0M, gen.}"};

    AxisSpec rapidityAxis = {200, -2.0f, 2.0f, "y"};
    AxisSpec phiAxis = {100, -TMath::Pi() / 2, 3. * TMath::Pi() / 2, "#varphi"};

    TString CutLabel[26] = {"All", "MassWin", "y", "EtaDau", "DCADauToPV", "CascCosPA", "V0CosPA", "DCACascDau", "DCAV0Dau", "rCasc", "rCascMax", "rV0", "rV0Max", "DCAV0ToPV", "LambdaMass", "TPCPr", "TPCPi", "TOFPr", "TOFPi", "TPCBach", "TOFBach", "ctau", "CompDecayMass", "Bach-baryon", "NTPCrows", "OOBRej"};
    TString CutLabelSummary[29] = {"MassWin", "y", "EtaDau", "dcapostopv", "dcanegtopv", "dcabachtopv", "CascCosPA", "V0CosPA", "DCACascDau", "DCAV0Dau", "rCasc", "rV0", "DCAV0ToPV", "LambdaMass", "TPCPr", "TPCPi", "TOFPr", "TOFPi", "TPCBach", "TOFBach", "proplifetime", "rejcomp", "ptthrtof", "bachBaryonCosPA", "bachBaryonDCAxyToPV", "NTPCrows", "OOBRej", "rCascMax", "rV0Max"};

    registry.add("hCandidate", "hCandidate", HistType::kTH1F, {{26, -0.5, 25.5}});
    for (Int_t n = 1; n <= registry.get<TH1>(HIST("hCandidate"))->GetNbinsX(); n++) {
      registry.get<TH1>(HIST("hCandidate"))->GetXaxis()->SetBinLabel(n, CutLabel[n - 1]);
    }

    registry.add("CascadeSelectionSummary", "CascadeSelectionSummary", HistType::kTH1F, {{29, -0.5, 28.5}});
    for (Int_t n = 1; n <= registry.get<TH1>(HIST("CascadeSelectionSummary"))->GetNbinsX(); n++) {
      registry.get<TH1>(HIST("CascadeSelectionSummary"))->GetXaxis()->SetBinLabel(n, CutLabelSummary[n - 1]);
    }
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(1, masswin);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(2, rap);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(3, etadau);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(4, dcapostopv);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(5, dcanegtopv);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(6, dcabachtopv);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(7, casccospa);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(8, v0cospa);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(9, dcacascdau);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(10, dcav0dau);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(11, cascradiusMin);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(12, v0radiusMin);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(13, dcav0topv);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(14, lambdamasswin);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(15, nsigmatpcPr);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(16, nsigmatpcPi);
    if (hastof)
      registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(17, nsigmatofPr);
    if (hastof)
      registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(18, nsigmatofPi);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(19, nsigmatpcKa);
    if (hastof)
      registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(20, nsigmatofKa);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(21, proplifetime);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(22, rejcomp);
    if (hastof)
      registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(23, ptthrtof);
    if (isSelectBachBaryon) {
      registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(24, bachBaryonCosPA);
      registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(25, bachBaryonDCAxyToPV);
    }
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(26, mintpccrrows);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(27, dooobrej);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(28, cascradiusMax);
    registry.get<TH1>(HIST("CascadeSelectionSummary"))->SetBinContent(29, v0radiusMax);

    registry.add("hPt", "hPt", {HistType::kTH1F, {ptAxis}});
    registry.add("hCascMinusInvMassvsPt", "hCascMinusInvMassvsPt", HistType::kTH2F, {ptAxis, massAxis});
    registry.add("hCascPlusInvMassvsPt", "hCascPlusInvMassvsPt", HistType::kTH2F, {ptAxis, massAxis});
    registry.add("hCascMinusInvMassvsPt_FT0M", "hCascMinusInvMassvsPt_FT0M", HistType::kTH3F, {centFT0MAxis, ptAxis, massAxis});
    registry.add("hCascPlusInvMassvsPt_FT0M", "hCascPlusInvMassvsPt_FT0M", HistType::kTH3F, {centFT0MAxis, ptAxis, massAxis});
    registry.add("hCascMinusInvMassvsPt_FV0A", "hCascMinusInvMassvsPt_FV0A", HistType::kTH3F, {centFV0AAxis, ptAxis, massAxis});
    registry.add("hCascPlusInvMassvsPt_FV0A", "hCascPlusInvMassvsPt_FV0A", HistType::kTH3F, {centFV0AAxis, ptAxis, massAxis});
    registry.add("hXiMinusInvMassvsPt_BefSels", "hXiMinusInvMassvsPt_BefSels", HistType::kTH2F, {ptAxis, ximassAxis});
    registry.add("hOmegaMinusInvMassvsPt_BefSels", "hOmegaMinusInvMassvsPt_BefSels", {HistType::kTH2F, {ptAxis, omegamassAxis}});
    registry.add("hXiPlusInvMassvsPt_BefSels", "hXiPlusInvMassvsPt_BefSels", {HistType::kTH2F, {ptAxis, ximassAxis}});
    registry.add("hOmegaPlusInvMassvsPt_BefSels", "hOmegaPlusInvMassvsPt_BefSels", {HistType::kTH2F, {ptAxis, omegamassAxis}});

    // topo
    registry.add("hDCANegToPV", "hDCANegToPV", {HistType::kTH2F, {ptAxisTopoVar, {200, -1.0f, 1.0f}}});
    registry.add("hDCAPosToPV", "hDCAPosToPV", {HistType::kTH2F, {ptAxisTopoVar, {200, -1.0f, 1.0f}}});
    registry.add("hDCABachToPV", "hDCABachToPV", {HistType::kTH2F, {ptAxisTopoVar, {200, -1.0f, 1.0f}}});
    registry.add("hCascCosPA", "hCascCosPA", {HistType::kTH2F, {ptAxisTopoVar, {100, 0.9f, 1.0f}}});
    registry.add("hV0CosPA", "hV0CosPA", {HistType::kTH2F, {ptAxisTopoVar, {100, 0.9f, 1.0f}}});
    registry.add("hCascRadius", "hCascRadius", {HistType::kTH2D, {ptAxisTopoVar, {500, 0.0f, 50.0f}}});
    registry.add("hV0Radius", "hV0Radius", {HistType::kTH2D, {ptAxisTopoVar, {500, 0.0f, 50.0f}}});
    registry.add("hDCACascDaughters", "hDCACascDaughters", {HistType::kTH2F, {ptAxisTopoVar, {55, 0.0f, 2.20f}}});
    registry.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH2F, {ptAxisTopoVar, {55, 0.0f, 2.20f}}});
    registry.add("hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH2F, {ptAxisTopoVar, {55, 0.0f, 2.20f}}});
    registry.add("hMassLambdaDau", "hMassLambdaDau", {HistType::kTH2F, {ptAxis, {60, 1.1f, 1.13f}}});

    registry.add("hBachBaryonCosPA", "hBachBaryonCosPA", {HistType::kTH2F, {ptAxisTopoVar, {100, 0.0f, 1.0f}}});
    registry.add("hBachBaryonDCAxyToPV", "hBachBaryonDCAxyToPV", {HistType::kTH2F, {ptAxisTopoVar, {300, -3.0f, 3.0f}}});

    // kine
    registry.add("hEtaMinus", "hEtaMinus", {HistType::kTH2F, {ptAxis, etaAxis}});
    registry.add("hEtaPlus", "hEtaPlus", {HistType::kTH2F, {ptAxis, etaAxis}});
    registry.add("hRapMinus", "hRapMinus", {HistType::kTH2F, {ptAxis, rapidityAxis}});
    registry.add("hRapMinus1D", "hRapMinus1D", {HistType::kTH1F, {rapidityAxis}});
    registry.add("hRapPlus", "hRapPlus", {HistType::kTH2F, {ptAxis, rapidityAxis}});
    registry.add("hRapPlus1D", "hRapPlus1D", {HistType::kTH1F, {rapidityAxis}});
    registry.add("hPhiMinus", "hPhiMinus", {HistType::kTH2F, {ptAxis, phiAxis}});
    registry.add("hPhiPlus", "hPhiPlus", {HistType::kTH2F, {ptAxis, phiAxis}});
    registry.add("hCtauMinus", "hCtauMinus", {HistType::kTH1F, {{100, 0.0f, 40.0f}}});
    registry.add("hCtauPlus", "hCtauPlus", {HistType::kTH1F, {{100, 0.0f, 40.0f}}});

    // daughter tracks
    registry.add("hTPCNSigmaPosPi", "hTPCNSigmaPosPi", {HistType::kTH2F, {ptAxisPID, {120, -6.0f, 6.0f}}});
    registry.add("hTPCNSigmaNegPi", "hTPCNSigmaNegPi", {HistType::kTH2F, {ptAxisPID, {120, -6.0f, 6.0f}}});
    registry.add("hTPCNSigmaPosPr", "hTPCNSigmaPosPr", {HistType::kTH2F, {ptAxisPID, {120, -6.0f, 6.0f}}});
    registry.add("hTPCNSigmaNegPr", "hTPCNSigmaNegPr", {HistType::kTH2F, {ptAxisPID, {120, -6.0f, 6.0f}}});
    registry.add("hTPCNSigmaBachPi", "hTPCNSigmaBachPi", {HistType::kTH2F, {ptAxisPID, {120, -6.0f, 6.0f}}});
    registry.add("hTPCNSigmaBachKa", "hTPCNSigmaBachKa", {HistType::kTH2F, {ptAxisPID, {120, -6.0f, 6.0f}}});
    registry.add("hTOFNSigmaPosPi", "hTOFNSigmaPosPi", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTOFNSigmaNegPi", "hTOFNSigmaNegPi", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTOFNSigmaPosPr", "hTOFNSigmaPosPr", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTOFNSigmaNegPr", "hTOFNSigmaNegPr", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTOFNSigmaBachPi", "hTOFNSigmaBachPi", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTOFNSigmaBachKa", "hTOFNSigmaBachKa", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hCascMinusEtaPos", "hCascMinusEtaPos", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hCascMinusEtaNeg", "hCascMinusEtaNeg", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hCascMinusEtaBach", "hCascMinusEtaBach", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});

    // Info for eff x acc from MC
    registry.add("hPtCascPlusTrueRec", "hPtCascPlusTrueRec", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});
    registry.add("hPtCascMinusTrueRec", "hPtCascMinusTrueRec", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});

    registry.add("hCascMinusMassvsPtTrueRec", "hCascMinusMassvsPtTrueRec", {HistType::kTH3F, {ptAxis, massAxis, centFT0MAxis}});
    registry.add("hCascPlusMassvsPtTrueRec", "hCascPlusMassvsPtTrueRec", {HistType::kTH3F, {ptAxis, massAxis, centFT0MAxis}});
    registry.add("hCascMinusMassvsPtBG", "hCascMinusMassvsPtBG", {HistType::kTH2F, {ptAxis, massAxis}});
    registry.add("hCascPlusMassvsPtBG", "hCascPlusMassvsPtBG", {HistType::kTH2F, {ptAxis, massAxis}});
    if (isMC) {
      registry.add("hPtXiPlusTrue", "hPtXiPlusTrue", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});
      registry.add("hPtXiMinusTrue", "hPtXiMinusTrue", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});
      registry.add("hPtOmegaPlusTrue", "hPtOmegaPlusTrue", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});
      registry.add("hPtOmegaMinusTrue", "hPtOmegaMinusTrue", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});
      registry.add("hPtXiPlusTrueAssocWithSelColl", "hPtXiPlusTrueAssocWithSelColl", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});
      registry.add("hPtXiMinusTrueAssocWithSelColl", "hPtXiMinusTrueAssocWithSelColl", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});
      registry.add("hPtOmegaPlusTrueAssocWithSelColl", "hPtOmegaPlusTrueAssocWithSelColl", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});
      registry.add("hPtOmegaMinusTrueAssocWithSelColl", "hPtOmegaMinusTrueAssocWithSelColl", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});
    }
  }

  void processRec(aod::MyCascades const& mycascades)
  {
    float invmass = 0;
    float ctau = 0;
    float rapidity = 0;
    bool isCandidate = 0;
    int counter = -1;
    bool isCorrectlyRec = 0;

    for (auto& candidate : mycascades) {

      switch (evSelFlag) {
        case 1: {
          if (!candidate.isINEL())
            continue;
          break;
        }
        case 2: {
          if (!candidate.isINELgt0())
            continue;
          break;
        }
        case 3: {
          if (!candidate.isINELgt1())
            continue;
          break;
        }
        default:
          LOGF(fatal, "incorrect evSelFlag in cascpostprocessing task");
          break;
      }

      counter = -1;
      registry.fill(HIST("hCandidate"), ++counter);

      // To have trace of how it was before selections
      if (candidate.sign() < 0) {
        registry.fill(HIST("hXiMinusInvMassvsPt_BefSels"), candidate.pt(), candidate.massxi());
        registry.fill(HIST("hOmegaMinusInvMassvsPt_BefSels"), candidate.pt(), candidate.massomega());
      }
      if (candidate.sign() > 0) {
        registry.fill(HIST("hXiPlusInvMassvsPt_BefSels"), candidate.pt(), candidate.massxi());
        registry.fill(HIST("hOmegaPlusInvMassvsPt_BefSels"), candidate.pt(), candidate.massomega());
      }

      if (isXi) {
        if (TMath::Abs(candidate.massxi() - pdgDB->Mass(3312)) > masswin)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
        if (TMath::Abs(candidate.rapxi()) > rap)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
      } else {
        if (TMath::Abs(candidate.massomega() - pdgDB->Mass(3334)) > masswin)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
        if (TMath::Abs(candidate.rapomega()) > rap)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
      }

      // Apply selections
      if (TMath::Abs(candidate.poseta()) > etadau || TMath::Abs(candidate.negeta()) > etadau || TMath::Abs(candidate.bacheta()) > etadau)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      if (TMath::Abs(candidate.dcanegtopv()) < dcanegtopv || TMath::Abs(candidate.dcapostopv()) < dcapostopv || TMath::Abs(candidate.dcabachtopv()) < dcabachtopv)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      if (candidate.casccospa() < casccospa)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      if (candidate.v0cospa() < v0cospa)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      if (candidate.dcacascdaughters() > dcacascdau)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      if (candidate.dcav0daughters() > dcav0dau)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      if (candidate.cascradius() < cascradiusMin)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      if (candidate.cascradius() > cascradiusMax)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      if (candidate.v0radius() < v0radiusMin)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      if (candidate.v0radius() > v0radiusMax)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      if (TMath::Abs(candidate.dcav0topv()) < dcav0topv)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      if (TMath::Abs(candidate.masslambdadau() - pdgDB->Mass(3122)) > lambdamasswin)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      if (candidate.sign() < 0) {
        if (TMath::Abs(candidate.ntpcsigmapospr()) > nsigmatpcPr)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
        if (TMath::Abs(candidate.ntpcsigmanegpi()) > nsigmatpcPi)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
      } else if (candidate.sign() > 0) {
        if (TMath::Abs(candidate.ntpcsigmanegpr()) > nsigmatpcPr)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
        if (TMath::Abs(candidate.ntpcsigmapospi()) > nsigmatpcPi)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
      }
      // TOF required only if hastof is set to one. In this case, daughter tracks with pt > threshold should have a hit in the tof
      if (hastof) {
        if (candidate.sign() < 0) {
          if (TMath::Abs(candidate.ntofsigmapospr()) > nsigmatofPr && candidate.pospt() > ptthrtof && candidate.poshastof())
            continue;
          registry.fill(HIST("hCandidate"), ++counter);
          if (TMath::Abs(candidate.ntofsigmanegpi()) > nsigmatofPi && candidate.negpt() > ptthrtof && candidate.neghastof())
            continue;
          registry.fill(HIST("hCandidate"), ++counter);
        } else if (candidate.sign() > 0) {
          if (TMath::Abs(candidate.ntofsigmanegpr()) > nsigmatofPr && candidate.negpt() > ptthrtof && candidate.neghastof())
            continue;
          registry.fill(HIST("hCandidate"), ++counter);
          if (TMath::Abs(candidate.ntofsigmapospi()) > nsigmatofPi && candidate.pospt() > ptthrtof && candidate.poshastof())
            continue;
          registry.fill(HIST("hCandidate"), ++counter);
        }
      } else {
        counter += 2;
      }
      if (isXi) {
        if (TMath::Abs(candidate.ntpcsigmabachpi()) > nsigmatpcPi)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
        if (hastof && TMath::Abs(candidate.ntofsigmabachpi()) > nsigmatofPi && candidate.bachpt() > ptthrtof && candidate.bachhastof())
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
        if (candidate.ctauxi() > proplifetime * ctauxiPDG)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
        if (TMath::Abs(candidate.massomega() - pdgDB->Mass(3334)) < rejcomp)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
        rapidity = candidate.rapxi();
        ctau = candidate.ctauxi();
        invmass = candidate.massxi();
      } else {
        if (TMath::Abs(candidate.ntpcsigmabachka()) > nsigmatpcKa)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
        if (hastof && TMath::Abs(candidate.ntofsigmabachka()) > nsigmatofKa && candidate.bachpt() > ptthrtof && candidate.bachhastof())
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
        if (candidate.ctauomega() > proplifetime * ctauomegaPDG)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
        if (TMath::Abs(candidate.massxi() - pdgDB->Mass(3312)) < rejcomp)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
        rapidity = candidate.rapomega();
        ctau = candidate.ctauomega();
        invmass = candidate.massomega();
      }
      if (isSelectBachBaryon && (candidate.bachBaryonCosPA() > bachBaryonCosPA || fabs(candidate.bachBaryonDCAxyToPV()) < bachBaryonDCAxyToPV)) { // Bach-baryon selection if required
        continue;
      }
      registry.fill(HIST("hCandidate"), ++counter);
      if (candidate.posntpccrrows() < mintpccrrows || candidate.negntpccrrows() < mintpccrrows || candidate.bachntpccrrows() < mintpccrrows)
        continue;
      registry.fill(HIST("hCandidate"), ++counter);
      bool kHasTOF = (candidate.poshastof() || candidate.neghastof() || candidate.bachhastof());
      bool kHasITS = ((candidate.positshits() > 1) || (candidate.negitshits() > 1) || (candidate.bachitshits() > 1));
      if (dooobrej == 1) {
        if (!kHasTOF && !kHasITS)
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
      } else if (dooobrej == 2) {
        if (!kHasTOF && (candidate.pt() > ptthrtof))
          continue;
        registry.fill(HIST("hCandidate"), ++counter);
      } else {
        ++counter;
      }

      registry.fill(HIST("hPt"), candidate.pt());
      registry.fill(HIST("hDCANegToPV"), candidate.pt(), candidate.dcanegtopv());
      registry.fill(HIST("hDCAPosToPV"), candidate.pt(), candidate.dcapostopv());
      registry.fill(HIST("hDCABachToPV"), candidate.pt(), candidate.dcabachtopv());
      registry.fill(HIST("hCascCosPA"), candidate.pt(), candidate.casccospa());
      registry.fill(HIST("hV0CosPA"), candidate.pt(), candidate.v0cospa());
      registry.fill(HIST("hCascRadius"), candidate.pt(), candidate.cascradius());
      registry.fill(HIST("hV0Radius"), candidate.pt(), candidate.v0radius());
      registry.fill(HIST("hDCACascDaughters"), candidate.pt(), candidate.dcacascdaughters());
      registry.fill(HIST("hDCAV0Daughters"), candidate.pt(), candidate.dcav0daughters());
      registry.fill(HIST("hDCAV0ToPV"), candidate.pt(), candidate.dcav0topv());
      registry.fill(HIST("hMassLambdaDau"), candidate.pt(), candidate.masslambdadau());

      registry.fill(HIST("hBachBaryonCosPA"), candidate.pt(), candidate.bachBaryonCosPA());
      registry.fill(HIST("hBachBaryonDCAxyToPV"), candidate.pt(), candidate.bachBaryonDCAxyToPV());
      if (candidate.sign() > 0) {
        registry.fill(HIST("hCtauPlus"), ctau);
        registry.fill(HIST("hEtaPlus"), candidate.pt(), candidate.eta());
        registry.fill(HIST("hRapPlus"), candidate.pt(), rapidity);
        registry.fill(HIST("hRapPlus1D"), rapidity);
        // registry.fill(HIST("hPhiPlus"), candidate.pt(), candidate.phi());
      } else {
        registry.fill(HIST("hCtauMinus"), ctau);
        registry.fill(HIST("hEtaMinus"), candidate.pt(), candidate.eta());
        registry.fill(HIST("hRapMinus"), candidate.pt(), rapidity);
        registry.fill(HIST("hRapMinus1D"), rapidity);
        // registry.fill(HIST("hPhiMinus"), candidate.pt(), candidate.phi());
      }

      if (isXi) {
        isCorrectlyRec = ((TMath::Abs(candidate.mcPdgCode()) == 3312) && (candidate.isPrimary() == 1)) ? 1 : 0;
        if (TMath::Abs(candidate.massxi() - pdgDB->Mass(3312)) < masswintpc) {
          isCandidate = 1;
        }
      } else if (!isXi) {
        isCorrectlyRec = ((TMath::Abs(candidate.mcPdgCode()) == 3334) && (candidate.isPrimary() == 1)) ? 1 : 0;
        if (TMath::Abs(candidate.massomega() - pdgDB->Mass(3334)) < masswintpc) {
          isCandidate = 1;
        }
      }
      if (isCandidate) {
        if (candidate.sign() < 0) {
          registry.fill(HIST("hTPCNSigmaPosPr"), candidate.pt(), candidate.ntpcsigmapospr());
          registry.fill(HIST("hTPCNSigmaNegPi"), candidate.pt(), candidate.ntpcsigmanegpi());
          registry.fill(HIST("hTOFNSigmaPosPr"), candidate.ntofsigmapospr());
          registry.fill(HIST("hTOFNSigmaNegPi"), candidate.ntofsigmanegpi());
          registry.fill(HIST("hCascMinusEtaPos"), candidate.poseta());
          registry.fill(HIST("hCascMinusEtaNeg"), candidate.negeta());
          registry.fill(HIST("hCascMinusEtaBach"), candidate.bacheta());
        } else {
          registry.fill(HIST("hTPCNSigmaPosPi"), candidate.pt(), candidate.ntpcsigmapospi());
          registry.fill(HIST("hTPCNSigmaNegPr"), candidate.pt(), candidate.ntpcsigmanegpr());
          registry.fill(HIST("hTOFNSigmaPosPi"), candidate.ntofsigmapospi());
          registry.fill(HIST("hTOFNSigmaNegPr"), candidate.ntofsigmanegpr());
        }
        if (isXi) {
          registry.fill(HIST("hTPCNSigmaBachPi"), candidate.pt(), candidate.ntpcsigmabachpi());
          registry.fill(HIST("hTOFNSigmaBachPi"), candidate.ntofsigmabachpi());
        } else {
          registry.fill(HIST("hTPCNSigmaBachKa"), candidate.pt(), candidate.ntpcsigmabachka());
          registry.fill(HIST("hTOFNSigmaBachKa"), candidate.ntofsigmabachka());
        }
      }
      // registry.fill(HIST("hPosITSHits"), candidate.positshits());
      // registry.fill(HIST("hNegITSHits"), candidate.negitshits());
      // registry.fill(HIST("hBachITSHits"), candidate.bachitshits());

      if (candidate.sign() < 0) {
        if (isCorrectlyRec) {
          registry.fill(HIST("hPtCascMinusTrueRec"), candidate.pt(), rapidity, candidate.centFT0M()); // 3rd axis is from MC calibration
          registry.fill(HIST("hCascMinusMassvsPtTrueRec"), candidate.pt(), invmass, candidate.centFT0M());
        } else {
          registry.fill(HIST("hCascMinusMassvsPtBG"), candidate.pt(), invmass);
        }
        registry.fill(HIST("hCascMinusInvMassvsPt"), candidate.pt(), invmass);
        registry.fill(HIST("hCascMinusInvMassvsPt_FT0M"), candidate.centFT0M(), candidate.pt(), invmass);
        registry.fill(HIST("hCascMinusInvMassvsPt_FV0A"), candidate.centFV0A(), candidate.pt(), invmass);
      }
      if (candidate.sign() > 0) {
        if (isCorrectlyRec) {
          registry.fill(HIST("hPtCascPlusTrueRec"), candidate.pt(), rapidity, candidate.centFT0M()); // 3rd axis is from MC calibration
          registry.fill(HIST("hCascPlusMassvsPtTrueRec"), candidate.pt(), invmass, candidate.centFT0M());
        } else {
          registry.fill(HIST("hCascPlusMassvsPtBG"), candidate.pt(), invmass);
        }
        registry.fill(HIST("hCascPlusInvMassvsPt"), candidate.pt(), invmass);
        registry.fill(HIST("hCascPlusInvMassvsPt_FT0M"), candidate.centFT0M(), candidate.pt(), invmass);
        registry.fill(HIST("hCascPlusInvMassvsPt_FV0A"), candidate.centFV0A(), candidate.pt(), invmass);
      }
    }
  }

  PROCESS_SWITCH(cascpostprocessing, processRec, "Process Run 3 reconstructed data", true);

  void processGen(aod::MyMCCascades const& myMCcascades)
  {
    for (auto& genCascade : myMCcascades) {
      if (genCascade.isPrimary() == 0)
        continue; // Consider only primaries

      switch (evSelFlag) {
        case 1: {
          if (!genCascade.isINEL())
            continue;
          break;
        }
        case 2: {
          if (!genCascade.isINELgt0())
            continue;
          break;
        }
        case 3: {
          if (!genCascade.isINELgt1())
            continue;
          break;
        }
        default:
          LOGF(fatal, "incorrect evSelFlag in cascpostprocessing task");
          break;
      }

      if (TMath::Abs(genCascade.y()) > rap)
        continue;

      // Histos of generated cascades from generated events with accepted z vrtx + chosen event type (evSelFlag) (for signal loss correction)
      if (genCascade.pdgCode() == -3312) {
        registry.fill(HIST("hPtXiPlusTrue"), genCascade.pt(), genCascade.y(), genCascade.centFT0M());
      }
      if (genCascade.pdgCode() == 3312) {
        registry.fill(HIST("hPtXiMinusTrue"), genCascade.pt(), genCascade.y(), genCascade.centFT0M());
      }
      if (genCascade.pdgCode() == -3334) {
        registry.fill(HIST("hPtOmegaPlusTrue"), genCascade.pt(), genCascade.y(), genCascade.centFT0M());
      }
      if (genCascade.pdgCode() == 3334) {
        registry.fill(HIST("hPtOmegaMinusTrue"), genCascade.pt(), genCascade.y(), genCascade.centFT0M());
      }

      // Histos of generated cascades from generated events with good z vrtx + chosen event type (evSelFlag) + associated to the accepted reconstructed event of the same type (for signal loss + efficiency x acceptance correction)
      switch (evSelFlag) {
        case 1: {
          if (!genCascade.isINELassoc())
            continue;
          break;
        }
        case 2: {
          if (!genCascade.isINELgt0assoc())
            continue;
          break;
        }
        case 3: {
          if (!genCascade.isINELgt1assoc())
            continue;
          break;
        }
        default:
          LOGF(fatal, "incorrect evSelFlag in cascpostprocessing task");
          break;
      }

      if (genCascade.pdgCode() == -3312) {
        registry.fill(HIST("hPtXiPlusTrueAssocWithSelColl"), genCascade.pt(), genCascade.y(), genCascade.centFT0M());
      }
      if (genCascade.pdgCode() == 3312) {
        registry.fill(HIST("hPtXiMinusTrueAssocWithSelColl"), genCascade.pt(), genCascade.y(), genCascade.centFT0M());
      }
      if (genCascade.pdgCode() == -3334) {
        registry.fill(HIST("hPtOmegaPlusTrueAssocWithSelColl"), genCascade.pt(), genCascade.y(), genCascade.centFT0M());
      }
      if (genCascade.pdgCode() == 3334) {
        registry.fill(HIST("hPtOmegaMinusTrueAssocWithSelColl"), genCascade.pt(), genCascade.y(), genCascade.centFT0M());
      }
    }
  }
  PROCESS_SWITCH(cascpostprocessing, processGen, "Process Run 3 MC generated data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascpostprocessing>(cfgc, TaskName{"lf-cascpostprocessing"})};
}
