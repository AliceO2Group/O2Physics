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
/// \brief this is a starting point for the third session of the tutorial
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/v0qaanalysis.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct v0postprocessing {

  Configurable<float> radius{"radius", 0.5, "Radius"};
  Configurable<float> maxradius{"maxradius", 100000, "Radius"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<double> cospaK0s{"cospaK0s", 0.97, "K0s CosPA"};
  Configurable<double> cospaLambda{"cospaLambda", 0.995, "Lambda CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> rap{"rap", 0.5, "Rapidity"};
  Configurable<float> ctauK0s{"ctauK0s", 20, "C tau K0s(cm)"};
  Configurable<float> ctauLambda{"ctauLambda", 30, "C tau Lambda (cm)"};
  Configurable<float> v0rejK0s{"v0rejK0s", 0.005, "V0 rej K0s"};
  Configurable<float> v0rejLambda{"v0rejLambda", 0.01, "V0 rej K0s"};
  Configurable<float> ntpcsigma{"ntpcsigma", 5, "N sigma TPC"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};
  Configurable<float> minITShits{"minITShits", 2, "min ITS hits"};
  Configurable<float> min_TPC_nClusters{"min_TPC_nClusters", 80, "min_TPC_nClusters"};
  Configurable<float> max_tpcSharedCls{"max_tpcSharedCls", 100, "max_tpcSharedCls"};
  Configurable<float> max_chi2_ITS{"max_chi2_ITS", 36, "max_chi2_ITS"};
  Configurable<float> max_chi2_TPC{"max_chi2_TPC", 4, "max_chi2_TPC"};
  Configurable<bool> isMC{"isMC", 1, "isMC"};
  Configurable<bool> evSel{"evSel", 1, "evSel"};
  Configurable<bool> hasTOF2Leg{"hasTOF2Leg", 0, "hasTOF2Leg"};
  Configurable<bool> hasTOF1Leg{"hasTOF1Leg", 0, "hasTOF1Leg"};
  Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2, "parameter Armenteros Cut"};
  Configurable<bool> doArmenterosCut{"doArmenterosCut", 1, "do Armenteros Cut for K0s"};
  Configurable<bool> doArmenterosCutLam{"doArmenterosCutLam", 1, "do Armenteros Cut for Lam"};
  Configurable<bool> doQA{"doQA", 1, "fill QA histograms"};

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    registry.add("hV0Cuts", ";Sel", {HistType::kTH1D, {{22, 0., 22.}}});
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(1, "all");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(2, "Event selection");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(3, "Radius");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(4, "Eta Daughters");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(5, "Dau DCA to PV");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(6, "DCA Daughters");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(7, "min ITS hits");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(8, "has TOF 1 Leg");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(9, "has TOF 2 Legs");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(10, "TPC NCl");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(11, "TPC Cls Shared");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(12, "ITS Chi2");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(13, "TPC Chi2");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(14, "cosPA K0s");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(15, "cosPA Lambda");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(16, "rapidity");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(17, "ctau K0s");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(18, "ctau Lambda");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(19, "v0 rej K0s");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(20, "v0 rej Lambda");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(21, "TPC nsigma Dau");
    registry.get<TH1>(HIST("hV0Cuts"))->GetXaxis()->SetBinLabel(22, "Armenteros-Podolansky");

    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(1, 1);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(2, evSel);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(3, radius);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(4, etadau);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(5, dcanegtopv);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(6, dcav0dau);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(7, minITShits);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(8, hasTOF1Leg);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(9, hasTOF2Leg);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(10, min_TPC_nClusters);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(11, max_tpcSharedCls);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(12, max_chi2_ITS);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(13, max_chi2_TPC);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(14, cospaK0s);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(15, cospaLambda);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(16, rap);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(17, ctauK0s);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(18, ctauLambda);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(19, v0rejK0s);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(20, v0rejLambda);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(21, ntpcsigma);
    registry.get<TH1>(HIST("hV0Cuts"))->SetBinContent(22, paramArmenterosCut * doArmenterosCut);

    registry.add("hMassK0Short", ";M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH1F, {{200, 0.4f, 0.6f}}});
    registry.add("hMassVsPtK0Short", ";p_{T} [GeV/c];M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH2F, {{250, 0.0f, 25.0f}, {200, 0.4f, 0.6f}}});
    registry.add("hMassVsPtK0ShortVsCentFT0M", ";p_{T} [GeV/c]; CentFT0M; M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 0.4f, 0.6f}}});
    registry.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
    registry.add("hMassVsPtLambda", "hMassVsPtLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});
    registry.add("hMassVsPtLambdaVsCentFT0M", ";p_{T} [GeV/c]; CentFT0M; M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 1.016f, 1.216f}}});
    registry.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
    registry.add("hMassVsPtAntiLambda", "hMassVsPtAntiLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});
    registry.add("hMassVsPtAntiLambdaVsCentFT0M", ";p_{T} [GeV/c]; CentFT0M; M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 1.016f, 1.216f}}});

    if (isMC) {
      registry.add("hMassK0Short_MC", ";M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH1F, {{200, 0.4f, 0.6f}}});
      registry.add("hMassVsPtK0ShortVsCentFT0M_MC", ";p_{T} [GeV/c];M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 0.4f, 0.6f}}});
      registry.add("hMassLambda_MC", "hMassLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
      registry.add("hMassVsPtLambdaVsCentFT0M_MC", ";p_{T} [GeV/c];M_{p^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 1.016f, 1.216f}}});
      registry.add("hMassAntiLambda_MC", "hMassAntiLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
      registry.add("hFDVsPtLambdaVsMotherPt_DoubleCharged_MC", ";p_{T} [GeV/c] (V0);p_{T}^{gen} [GeV/c] (#Xi^{-}); percentile", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {250, 0.0f, 25.0f}, {100, 0.f, 100.f}}});
      registry.add("hFDVsPtLambdaVsMotherPt_MCRatio_MC", ";p_{T} [GeV/c] (V0);p_{T}^{gen} [GeV/c] (#Xi^{-/0}); percentile", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {250, 0.0f, 25.0f}, {100, 0.f, 100.f}}});
      registry.add("hMassVsPtAntiLambdaVsCentFT0M_MC", ";p_{T} [GeV/c];M_{p^{-}#pi^{+}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 1.016f, 1.216f}}});
      registry.add("hFDVsPtAntiLambdaVsMotherPt_DoubleCharged_MC", ";p_{T} [GeV/c] (V0);p_{T}^{gen} [GeV/c] (#bar{#Xi}^{+});percentile", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {250, 0.0f, 25.0f}, {100, 0.f, 100.f}}});
      registry.add("hFDVsPtAntiLambdaVsMotherPt_MCRatio_MC", ";p_{T} [GeV/c] (V0);p_{T}^{gen} [GeV/c] (#bar{#Xi}^{+/0});percentile", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {250, 0.0f, 25.0f}, {100, 0.f, 100.f}}});
    }

    if (doQA) {
      registry.add("QA/hK0sSelection", ";Sel", {HistType::kTH1D, {{22, 0., 22.}}});
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(1, "all");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(2, "Event selection");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(3, "Radius");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(4, "Eta Daughters");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(5, "Dau DCA to PV");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(6, "DCA Daughters");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(7, "min ITS hits");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(8, "has TOF 1 Leg");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(9, "has TOF 2 Legs");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(10, "TPC NCl");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(11, "TPC Cls Shared");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(12, "ITS Chi2");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(13, "TPC Chi2");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(14, "cosPA");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(15, "rapidity");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(16, "ctau");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(17, "v0 rej");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(18, "TPC nsigma Neg Dau");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(19, "TPC nsigma Pos Dau");
      registry.get<TH1>(HIST("QA/hK0sSelection"))->GetXaxis()->SetBinLabel(20, "Armenteros-Podolansky");

      // common
      registry.add("QA/hV0_EvFlag", "hV0_EvFlag", {HistType::kTH1D, {{2, 0.0f, 2.0f}}});
      registry.add("QA/hV0_Radius", "hV0_Radius", {HistType::kTH1D, {{1000, 0.0f, 100.0f}}});
      registry.add("QA/hV0_DCADauToPV", "hV0_DCADauToPV", {HistType::kTH1F, {{200, -1.0f, 1.0f}}});
      registry.add("QA/hV0_DCADaughters", "hV0_DCADaughters", {HistType::kTH1F, {{200, 0.0f, 2.0f}}});
      registry.add("QA/hV0_EtaDau", "hV0_EtaDau", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      registry.add("QA/hV0_ITShits", "hV0_ITShits", {HistType::kTH1F, {{10, .0f, 10.0f}}});
      registry.add("QA/hV0_TPCNCls", "hV0_TPCNCls", {HistType::kTH1F, {{200, .0f, 200.0f}}});
      registry.add("QA/hV0_TPCNClsShared", "hV0_TPCNClsShared", {HistType::kTH1F, {{150, .0f, 1.5f}}});
      registry.add("QA/hV0_ITSChi2", "hV0_ITSChi2", {HistType::kTH1F, {{10, .0f, 10.0f}}});
      registry.add("QA/hV0_TPCChi2", "hV0_TPCChi2", {HistType::kTH1F, {{100, .0f, 100.0f}}});
      // K0s
      registry.add("QA/hK0s_ArmenterosPodolanski", "QA/hK0s_ArmenterosPodolanski", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}});
      registry.add("QA/hK0s_Rap", "hK0s_Rap", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      registry.add("QA/hK0s_CosPA", "hK0s_CosPA", {HistType::kTH1F, {{100, 0.95f, 1.0f}}});
      registry.add("QA/hK0s_Ctau", "hK0s_Ctau", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      registry.add("QA/hK0s_TPCNSigmaPosPi", "hK0s_TPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("QA/hK0s_TPCNSigmaNegPi", "hK0s_TPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      // Lambda
      registry.add("QA/hLambda_ArmenterosPodolanski", "QA/hLambda_ArmenterosPodolanski", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}});
      registry.add("QA/hLambda_Rap", "hLambda_Rap", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      registry.add("QA/hLambda_CosPA", "hLambda_CosPA", {HistType::kTH1F, {{100, 0.95f, 1.0f}}});
      registry.add("QA/hLambda_Ctau", "hLambda_Ctau", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      registry.add("QA/hLambda_TPCNSigmaPosPi", "hLambda_TPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("QA/hLambda_TPCNSigmaNegPi", "hLambda_TPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      // AntiLambda
      registry.add("QA/hAntiLambda_ArmenterosPodolanski", "QA/hAntiLambda_ArmenterosPodolanski", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}});
      registry.add("QA/hAntiLambda_Rap", "hAntiLambda_Rap", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      registry.add("QA/hAntiLambda_CosPA", "hAntiLambda_CosPA", {HistType::kTH1F, {{100, 0.95f, 1.0f}}});
      registry.add("QA/hAntiLambda_Ctau", "hAntiLambda_Ctau", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      registry.add("QA/hAntiLambda_TPCNSigmaPosPi", "hAntiLambda_TPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("QA/hAntiLambda_TPCNSigmaNegPi", "hAntiLambda_TPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});

      // common
      registry.add("QA/hV0_Sel_EvFlag", "hV0_Sel_EvFlag", {HistType::kTH1D, {{2, 0.0f, 2.0f}}});
      registry.add("QA/hV0_Sel_Radius", "hV0_Sel_Radius", {HistType::kTH1D, {{1000, 0.0f, 100.0f}}});
      registry.add("QA/hV0_Sel_DCADauToPV", "hV0_Sel_DCADauToPV", {HistType::kTH1F, {{200, -1.0f, 1.0f}}});
      registry.add("QA/hV0_Sel_DCADaughters", "hV0_Sel_DCADaughters", {HistType::kTH1F, {{200, 0.0f, 2.0f}}});
      registry.add("QA/hV0_Sel_EtaDau", "hV0_Sel_EtaDau", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      registry.add("QA/hV0_Sel_ITShits", "hV0_Sel_ITShits", {HistType::kTH1F, {{10, .0f, 10.0f}}});
      registry.add("QA/hV0_Sel_TPCNCls", "hV0_Sel_TPCNCls", {HistType::kTH1F, {{200, .0f, 200.0f}}});
      registry.add("QA/hV0_Sel_TPCNClsShared", "hV0_Sel_TPCNClsShared", {HistType::kTH1F, {{150, .0f, 1.5f}}});
      registry.add("QA/hV0_Sel_ITSChi2", "hV0_Sel_ITSChi2", {HistType::kTH1F, {{10, .0f, 10.0f}}});
      registry.add("QA/hV0_Sel_TPCChi2", "hV0_Sel_TPCChi2", {HistType::kTH1F, {{100, .0f, 100.0f}}});
      // K0s
      registry.add("QA/hK0s_Sel_ArmenterosPodolanski", "QA/hK0s_ArmenterosPodolanski", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}});
      registry.add("QA/hK0s_Sel_Rap", "hK0s_Sel_Rap", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      registry.add("QA/hK0s_Sel_CosPA", "hK0s_Sel_CosPA", {HistType::kTH1F, {{100, 0.95f, 1.0f}}});
      registry.add("QA/hK0s_Sel_Ctau", "hK0s_Sel_Ctau", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      registry.add("QA/hK0s_Sel_TPCNSigmaPosPi", "hK0s_Sel_TPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("QA/hK0s_Sel_TPCNSigmaNegPi", "hK0s_Sel_TPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      // Lambda
      registry.add("QA/hLambda_Sel_ArmenterosPodolanski", "QA/hLambda_Sel_ArmenterosPodolanski", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}});
      registry.add("QA/hLambda_Sel_Rap", "hLambda_Sel_Rap", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      registry.add("QA/hLambda_Sel_CosPA", "hLambda_Sel_CosPA", {HistType::kTH1F, {{100, 0.95f, 1.0f}}});
      registry.add("QA/hLambda_Sel_Ctau", "hLambda_Sel_Ctau", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      registry.add("QA/hLambda_Sel_TPCNSigmaPosPr", "hLambda_Sel_TPCNSigmaPosPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("QA/hLambda_Sel_TPCNSigmaNegPi", "hLambda_Sel_TPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      // AntiLambda
      registry.add("QA/hAntiLambda_Sel_ArmenterosPodolanski", "QA/hAntiLambda_Sel_ArmenterosPodolanski", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}});
      registry.add("QA/hAntiLambda_Sel_Rap", "hAntiLambda_Sel_Rap", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      registry.add("QA/hAntiLambda_Sel_CosPA", "hAntiLambda_Sel_CosPA", {HistType::kTH1F, {{100, 0.95f, 1.0f}}});
      registry.add("QA/hAntiLambda_Sel_Ctau", "hAntiLambda_Sel_Ctau", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      registry.add("QA/hAntiLambda_Sel_TPCNSigmaPosPi", "hAntiLambda_Sel_TPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("QA/hAntiLambda_Sel_TPCNSigmaNegPr", "hAntiLambda_Sel_TPCNSigmaNegPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    }
  }

  // V0 selection
  template <typename TV0Type>
  bool QAK0s(TV0Type const& candidate)
  {
    if (candidate.v0cospa() <= cospaK0s)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 13.5);

    if (candidate.rapk0short() >= rap)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 14.5);

    if (candidate.ctauk0short() >= ctauK0s)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 15.5);

    if (std::abs(candidate.masslambda() - o2::constants::physics::MassLambda0) <= v0rejK0s)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 16.5);

    if (std::abs(candidate.ntpcsigmanegpi()) > ntpcsigma)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 17.5);

    if (std::abs(candidate.ntpcsigmapospi()) > ntpcsigma)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 18.5);

    if (doArmenterosCut && candidate.qtarm() < (paramArmenterosCut * std::abs(candidate.alpha())))
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 19.5);

    return true;
  }

  // V0 selection
  template <typename TV0Type>
  bool AcceptV0(TV0Type const& candidate)
  {
    if (evSel && candidate.evflag() < 1)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 1.5);

    if (candidate.v0radius() < radius && candidate.v0radius() > maxradius)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 2.5);

    if (std::abs(candidate.v0poseta()) > etadau)
      return false;
    if (std::abs(candidate.v0negeta()) > etadau)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 3.5);

    if (std::abs(candidate.v0dcanegtopv()) < dcanegtopv)
      return false;
    if (std::abs(candidate.v0dcapostopv()) < dcapostopv)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 4.5);

    if (candidate.v0dcav0daughters() > dcav0dau)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 5.5);

    if (candidate.v0positshits() < minITShits)
      return false;
    if (candidate.v0negitshits() < minITShits)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 6.5);

    if (hasTOF1Leg && !candidate.poshastof() && !candidate.neghastof())
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 7.5);

    if (hasTOF2Leg && (!candidate.poshastof() || !candidate.neghastof()))
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 8.5);

    if (candidate.v0postpcCrossedRows() < min_TPC_nClusters)
      return false;
    if (candidate.v0negtpcCrossedRows() < min_TPC_nClusters)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 9.5);

    if (candidate.v0postpcNClsShared() > max_tpcSharedCls)
      return false;
    if (candidate.v0negtpcNClsShared() > max_tpcSharedCls)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 10.5);

    if (candidate.v0positsChi2NCl() > max_chi2_ITS)
      return false;
    if (candidate.v0negitsChi2NCl() > max_chi2_ITS)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 11.5);

    if (candidate.v0postpcChi2NCl() > max_chi2_TPC)
      return false;
    if (candidate.v0negtpcChi2NCl() > max_chi2_TPC)
      return false;
    registry.fill(HIST("QA/hK0sSelection"), 12.5);

    return true;
  }

  void process(aod::MyV0Candidates const& myv0s)
  {
    for (auto& candidate : myv0s) {

      if (doQA) {
        registry.fill(HIST("QA/hK0sSelection"), 0.5);
        registry.fill(HIST("QA/hV0_EvFlag"), candidate.evflag());
        registry.fill(HIST("QA/hV0_Radius"), candidate.v0radius());
        registry.fill(HIST("QA/hV0_DCADauToPV"), candidate.v0dcanegtopv());
        registry.fill(HIST("QA/hV0_DCADaughters"), candidate.v0dcav0daughters());
        registry.fill(HIST("QA/hV0_EtaDau"), candidate.v0poseta());
        registry.fill(HIST("QA/hV0_EtaDau"), candidate.v0negeta());
        registry.fill(HIST("QA/hV0_ITShits"), candidate.v0negitshits());
        registry.fill(HIST("QA/hV0_TPCNCls"), candidate.v0postpcCrossedRows());
        registry.fill(HIST("QA/hV0_TPCNCls"), candidate.v0negtpcCrossedRows());
        registry.fill(HIST("QA/hV0_TPCNClsShared"), candidate.v0postpcNClsShared());
        registry.fill(HIST("QA/hV0_TPCNClsShared"), candidate.v0negtpcNClsShared());
        registry.fill(HIST("QA/hV0_ITSChi2"), candidate.v0positsChi2NCl());
        registry.fill(HIST("QA/hV0_ITSChi2"), candidate.v0negitsChi2NCl());
        registry.fill(HIST("QA/hV0_TPCChi2"), candidate.v0postpcChi2NCl());
        registry.fill(HIST("QA/hV0_TPCChi2"), candidate.v0negtpcChi2NCl());
        registry.fill(HIST("QA/hK0s_ArmenterosPodolanski"), candidate.alpha(), candidate.qtarm());
        registry.fill(HIST("QA/hK0s_CosPA"), candidate.v0cospa());
        registry.fill(HIST("QA/hK0s_Rap"), candidate.rapk0short());
        registry.fill(HIST("QA/hK0s_Ctau"), candidate.ctauk0short());
        registry.fill(HIST("QA/hK0s_TPCNSigmaPosPi"), candidate.ntpcsigmapospi());
        registry.fill(HIST("QA/hK0s_TPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
        registry.fill(HIST("QA/hLambda_ArmenterosPodolanski"), candidate.alpha(), candidate.qtarm());
        registry.fill(HIST("QA/hLambda_CosPA"), candidate.v0cospa());
        registry.fill(HIST("QA/hLambda_Rap"), candidate.rapk0short());
        registry.fill(HIST("QA/hLambda_Ctau"), candidate.ctauk0short());
        registry.fill(HIST("QA/hLambda_TPCNSigmaPosPi"), candidate.ntpcsigmapospi());
        registry.fill(HIST("QA/hLambda_TPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
        registry.fill(HIST("QA/hAntiLambda_ArmenterosPodolanski"), candidate.alpha(), candidate.qtarm());
        registry.fill(HIST("QA/hAntiLambda_CosPA"), candidate.v0cospa());
        registry.fill(HIST("QA/hAntiLambda_Rap"), candidate.rapk0short());
        registry.fill(HIST("QA/hAntiLambda_Ctau"), candidate.ctauk0short());
        registry.fill(HIST("QA/hAntiLambda_TPCNSigmaPosPi"), candidate.ntpcsigmapospi());
        registry.fill(HIST("QA/hAntiLambda_TPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
      }

      // Apply common V0 selection
      if (!AcceptV0(candidate)) {
        continue;
      }

      QAK0s(candidate);

      if (doQA) {
        registry.fill(HIST("QA/hV0_Sel_EvFlag"), candidate.evflag());
        registry.fill(HIST("QA/hV0_Sel_Radius"), candidate.v0radius());
        registry.fill(HIST("QA/hV0_Sel_DCADauToPV"), candidate.v0dcanegtopv());
        registry.fill(HIST("QA/hV0_Sel_DCADaughters"), candidate.v0dcav0daughters());
        registry.fill(HIST("QA/hV0_Sel_EtaDau"), candidate.v0poseta());
        registry.fill(HIST("QA/hV0_Sel_EtaDau"), candidate.v0negeta());
        registry.fill(HIST("QA/hV0_Sel_ITShits"), candidate.v0negitshits());
        registry.fill(HIST("QA/hV0_Sel_ITShits"), candidate.v0positshits());
        registry.fill(HIST("QA/hV0_Sel_TPCNCls"), candidate.v0postpcCrossedRows());
        registry.fill(HIST("QA/hV0_Sel_TPCNCls"), candidate.v0negtpcCrossedRows());
        registry.fill(HIST("QA/hV0_Sel_TPCNClsShared"), candidate.v0postpcNClsShared());
        registry.fill(HIST("QA/hV0_Sel_TPCNClsShared"), candidate.v0negtpcNClsShared());
        registry.fill(HIST("QA/hV0_Sel_ITSChi2"), candidate.v0positsChi2NCl());
        registry.fill(HIST("QA/hV0_Sel_ITSChi2"), candidate.v0negitsChi2NCl());
        registry.fill(HIST("QA/hV0_Sel_TPCChi2"), candidate.v0postpcChi2NCl());
        registry.fill(HIST("QA/hV0_Sel_TPCChi2"), candidate.v0negtpcChi2NCl());
      }

      //////////////////////////////////
      //////////// K0Short /////////////
      //////////////////////////////////

      if (candidate.v0cospa() > cospaK0s &&
          std::abs(candidate.rapk0short()) < rap &&
          candidate.ctauk0short() < ctauK0s &&
          std::abs(candidate.massk0short() - o2::constants::physics::MassK0Short) < 0.1 &&
          std::abs(candidate.masslambda() - o2::constants::physics::MassLambda0) > v0rejK0s &&
          std::abs(candidate.ntpcsigmanegpi()) <= ntpcsigma &&
          std::abs(candidate.ntpcsigmapospi()) <= ntpcsigma &&
          (!doArmenterosCut || candidate.qtarm() > (paramArmenterosCut * std::abs(candidate.alpha())))) {

        registry.fill(HIST("hMassK0Short"), candidate.massk0short());
        registry.fill(HIST("hMassVsPtK0Short"), candidate.v0pt(), candidate.massk0short());
        registry.fill(HIST("hMassVsPtK0ShortVsCentFT0M"), candidate.v0pt(), candidate.multft0m(), candidate.massk0short());

        if (isMC &&
            candidate.pdgcode() == 310 &&
            candidate.isdauk0short() &&
            candidate.isphysprimary() == 1) {

          registry.fill(HIST("hMassK0Short_MC"), candidate.massk0short());
          registry.fill(HIST("hMassVsPtK0ShortVsCentFT0M_MC"), candidate.v0pt(), candidate.multft0m(), candidate.massk0short());
        }

        if (doQA) {
          registry.fill(HIST("QA/hK0s_Sel_ArmenterosPodolanski"), candidate.alpha(), candidate.qtarm());
          registry.fill(HIST("QA/hK0s_Sel_Rap"), candidate.rapk0short());
          registry.fill(HIST("QA/hK0s_Sel_CosPA"), candidate.v0cospa());
          registry.fill(HIST("QA/hK0s_Sel_Ctau"), candidate.ctauk0short());
          registry.fill(HIST("QA/hK0s_Sel_TPCNSigmaPosPi"), candidate.ntpcsigmapospi());
          registry.fill(HIST("QA/hK0s_Sel_TPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
        }
      }

      //////////////////////////////////
      ////// Lambda / AntiLambda ///////
      //////////////////////////////////

      if (candidate.v0cospa() > cospaLambda &&
          std::abs(candidate.raplambda()) < rap &&
          std::abs(candidate.massk0short() - o2::constants::physics::MassK0Short) > v0rejLambda) {

        //////////////////////////////////
        ///////////// Lambda /////////////
        //////////////////////////////////

        if (std::abs(candidate.ntpcsigmanegpi()) <= ntpcsigma &&
            std::abs(candidate.ntpcsigmapospr()) <= ntpcsigma &&
            candidate.ctaulambda() < ctauLambda &&
            std::abs(candidate.masslambda() - o2::constants::physics::MassLambda0) < 0.075 &&
            (!doArmenterosCutLam || candidate.qtarm() < (paramArmenterosCut * std::abs(candidate.alpha())))) {

          registry.fill(HIST("hMassLambda"), candidate.masslambda());
          registry.fill(HIST("hMassVsPtLambda"), candidate.v0pt(), candidate.masslambda());
          registry.fill(HIST("hMassVsPtLambdaVsCentFT0M"), candidate.v0pt(), candidate.multft0m(), candidate.masslambda());

          if (isMC && candidate.pdgcode() == 3122 && candidate.isdaulambda()) {

            if (candidate.isphysprimary() == 1) {
              registry.fill(HIST("hMassLambda_MC"), candidate.masslambda());
              registry.fill(HIST("hMassVsPtLambdaVsCentFT0M_MC"), candidate.v0pt(), candidate.multft0m(), candidate.masslambda());
            } else if (std::abs(candidate.masslambda() - o2::constants::physics::MassLambda0) < 0.01) {
              if (candidate.pdgcodemother() == 3312) {
                registry.fill(HIST("hFDVsPtLambdaVsMotherPt_DoubleCharged_MC"), candidate.v0pt(), candidate.v0motherpt(), candidate.multft0m());
              }
              if (candidate.pdgcodemother() == 3312 || candidate.pdgcodemother() == 3322) {
                registry.fill(HIST("hFDVsPtLambdaVsMotherPt_MCRatio_MC"), candidate.v0pt(), candidate.v0motherpt(), candidate.multft0m());
              }
            }
          }

          if (doQA) {
            registry.fill(HIST("QA/hLambda_Sel_ArmenterosPodolanski"), candidate.alpha(), candidate.qtarm());
            registry.fill(HIST("QA/hLambda_Sel_Rap"), candidate.rapk0short());
            registry.fill(HIST("QA/hLambda_Sel_CosPA"), candidate.v0cospa());
            registry.fill(HIST("QA/hLambda_Sel_Ctau"), candidate.ctauk0short());
            registry.fill(HIST("QA/hLambda_Sel_TPCNSigmaPosPr"), candidate.ntpcsigmapospr());
            registry.fill(HIST("QA/hLambda_Sel_TPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
          }
        }

        //////////////////////////////////
        /////////// AntiLambda ///////////
        //////////////////////////////////

        if (std::abs(candidate.ntpcsigmanegpr()) <= ntpcsigma &&
            std::abs(candidate.ntpcsigmapospi()) <= ntpcsigma &&
            candidate.ctauantilambda() < ctauLambda &&
            std::abs(candidate.massantilambda() - o2::constants::physics::MassLambda0) < 0.075 &&
            (!doArmenterosCutLam || candidate.qtarm() < (paramArmenterosCut * std::abs(candidate.alpha())))) {

          registry.fill(HIST("hMassAntiLambda"), candidate.massantilambda());
          registry.fill(HIST("hMassVsPtAntiLambda"), candidate.v0pt(), candidate.massantilambda());
          registry.fill(HIST("hMassVsPtAntiLambdaVsCentFT0M"), candidate.v0pt(), candidate.multft0m(), candidate.massantilambda());

          if (isMC && candidate.pdgcode() == -3122 && candidate.isdauantilambda()) {

            if (candidate.isphysprimary() == 1) {
              registry.fill(HIST("hMassAntiLambda_MC"), candidate.massantilambda());
              registry.fill(HIST("hMassVsPtAntiLambdaVsCentFT0M_MC"), candidate.v0pt(), candidate.v0motherpt(), candidate.massantilambda());
            } else if (std::abs(candidate.massantilambda() - o2::constants::physics::MassLambda0) < 0.01) {
              if (candidate.pdgcodemother() == -3312) {
                registry.fill(HIST("hFDVsPtAntiLambdaVsMotherPt_DoubleCharged_MC"), candidate.v0pt(), candidate.v0motherpt(), candidate.multft0m());
              }
              if (candidate.pdgcodemother() == -3312 || candidate.pdgcodemother() == -3322) {
                registry.fill(HIST("hFDVsPtAntiLambdaVsMotherPt_MCRatio_MC"), candidate.v0pt(), candidate.v0motherpt(), candidate.multft0m());
              }
            }
          }

          if (doQA) {
            registry.fill(HIST("QA/hAntiLambda_Sel_ArmenterosPodolanski"), candidate.alpha(), candidate.qtarm());
            registry.fill(HIST("QA/hAntiLambda_Sel_Rap"), candidate.rapk0short());
            registry.fill(HIST("QA/hAntiLambda_Sel_CosPA"), candidate.v0cospa());
            registry.fill(HIST("QA/hAntiLambda_Sel_Ctau"), candidate.ctauk0short());
            registry.fill(HIST("QA/hAntiLambda_Sel_TPCNSigmaPosPi"), candidate.ntpcsigmapospi());
            registry.fill(HIST("QA/hAntiLambda_Sel_TPCNSigmaNegPr"), candidate.ntpcsigmanegpr());
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0postprocessing>(cfgc, TaskName{"lf-v0postprocessing"})};
}
