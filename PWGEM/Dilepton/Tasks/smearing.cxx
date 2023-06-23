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
//
//
// Analysis task to produce smeared pt,eta,phi for electrons/muons in dilepton analysis
//    Please write to: daiki.sekihata@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGEM/Dilepton/Utils/MomentumSmearer.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

struct ApplySmearing {
  Produces<aod::SmearedTracks> smearedtrack;

  // Run for electrons or muons (For the moment the task is not designed for both at the same time)
  Configurable<int> fPdgCode{"cfgPdgCode", 11, "Set the type of particle to be smeared"};
  // Maps
  Configurable<std::string> fConfigResFileName{"cfgResFileName", "", "name of resolution file"};
  Configurable<std::string> fConfigResPtHistName{"cfgResPtHistName", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
  Configurable<std::string> fConfigResEtaHistName{"cfgResEtaHistName", "EtaResArr", "histogram name for eta in resolution file"};
  Configurable<std::string> fConfigResPhiPosHistName{"cfgResPhiPosHistName", "PhiPosResArr", "histogram name for phi pos in resolution file"};
  Configurable<std::string> fConfigResPhiNegHistName{"cfgResPhiNegHistName", "PhiEleResArr", "hisogram for phi neg in resolution file"};

  MomentumSmearer smearer;

  void init(InitContext& context)
  {
    smearer.setResFileName(TString(fConfigResFileName));
    smearer.setResPtHistName(TString(fConfigResPtHistName));
    smearer.setResEtaHistName(TString(fConfigResEtaHistName));
    smearer.setResPhiPosHistName(TString(fConfigResPhiPosHistName));
    smearer.setResPhiNegHistName(TString(fConfigResPhiNegHistName));
    smearer.init();
  }

  template <typename TTracksMC>
  void applySmearing(TTracksMC const& tracksMC)
  {
    for (auto& mctrack : tracksMC) {
      float ptgen = mctrack.pt();
      float etagen = mctrack.eta();
      float phigen = mctrack.phi();

      int pdgCode = mctrack.pdgCode();
      if (abs(pdgCode) == fPdgCode) {
        int ch = -1;
        if (pdgCode < 0) {
          ch = 1;
        }
        // apply smearing for electrons or muons.
        float ptsmeared, etasmeared, phismeared;
        smearer.applySmearing(ch, ptgen, etagen, phigen, ptsmeared, etasmeared, phismeared);
        smearedtrack(ptsmeared, etasmeared, phismeared);
      } else {
        // don't apply smearing
        smearedtrack(ptgen, etagen, phigen);
      }
    }
  }

  void processMCanalysis(ReducedMCTracks const& tracksMC)
  {
    applySmearing(tracksMC);
  }

  void processCocktail(aod::McParticles const& tracksMC)
  {
    applySmearing(tracksMC);
  }

  void processDummyCocktail(aod::McParticles const& tracksMC) {}

  void processDummyMCanalysis(ReducedMCTracks const& tracksMC) {}

  PROCESS_SWITCH(ApplySmearing, processMCanalysis, "Run for MC analysis", false);
  PROCESS_SWITCH(ApplySmearing, processCocktail, "Run for cocktail analysis", false);
  PROCESS_SWITCH(ApplySmearing, processDummyMCanalysis, "Dummy process function", false);
  PROCESS_SWITCH(ApplySmearing, processDummyCocktail, "Dummy process function", true);
};

struct CheckSmearing {
  using MyReducedTracks = soa::Join<ReducedMCTracks, SmearedTracks>;
  using MyCocktailTracks = soa::Join<aod::McParticles, SmearedTracks>;

  // Run for electrons or muons
  Configurable<int> fPdgCode{"cfgPdgCode", 11, "Set the type of particle to be checked"};

  // Resolution histos as cross check
  Configurable<bool> fConfigUsePtVecRes{"cfgUsePtVecRes", true, "If true, non-linear pt bins predefined in res histos"};
  ConfigurableAxis ptResBinsVec{"ptResBinsVec", {0., 0., 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.20, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285, 0.29, 0.295, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80, 2.00, 2.40, 2.80, 3.20, 3.70, 4.50, 6.00, 8.00, 10., 12.0, 14., 16., 18., 20.}, "Pt binning vector for resolution"};
  ConfigurableAxis ptResBins{"ptResBins", {20, 0.f, 20.f}, "Pt binning for resolution"};
  ConfigurableAxis deltaptResBins{"deltaptResBins", {500, -1.f, 1.f}, "DeltaPt binning for resolution"};
  ConfigurableAxis deltaetaResBins{"deltaetaResBins", {500, -0.5f, 0.5f}, "DeltaEta binning for resolution"};
  ConfigurableAxis deltaphiResBins{"deltaphiResBins", {500, -0.5f, 0.5f}, "DeltaPhi binning for resolution"};

  HistogramRegistry registry{"HistoAnalysisTrackSelection", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    registry.add<TH2>("hCorrelation_Pt", "pT correlation", {HistType::kTH2F, {{1000, 0.0f, 10.0f}, {1000, 0.0f, 10.0f}}});
    registry.add<TH2>("hCorrelation_Eta", "eta correlation", {HistType::kTH2F, {{200, -1.0f, +1.0f}, {200, -1.0f, +1.0f}}});
    registry.add<TH2>("hCorrelation_Phi", "phi correlation", {HistType::kTH2F, {{100, 0.0f, TMath::TwoPi()}, {100, 0.0f, TMath::TwoPi()}}});

    // Binning for resolution
    AxisSpec axisPtRes{ptResBins, "#it{p}^{gen}_{T,e} (GeV/#it{c})"};
    AxisSpec axisDeltaptRes{deltaptResBins, "(p^{gen}_{T} - p^{rec}_{T}) / p^{gen}_{T}"};
    AxisSpec axisDeltaetaRes{deltaetaResBins, "#eta^{gen} - #eta^{rec}"};
    AxisSpec axisDeltaphiRes{deltaphiResBins, "#varphi^{gen} - #varphi^{rec} (rad)"};

    if (!fConfigUsePtVecRes) {
      registry.add<TH2>("PtGen_DeltaPtOverPtGen", "", HistType::kTH2D, {axisPtRes, axisDeltaptRes}, true);
      registry.add<TH2>("PtGen_DeltaEta", "", HistType::kTH2D, {axisPtRes, axisDeltaetaRes}, true);
      registry.add<TH2>("PtGen_DeltaPhi_Ele", "", HistType::kTH2D, {axisPtRes, axisDeltaphiRes}, true);
      registry.add<TH2>("PtGen_DeltaPhi_Pos", "", HistType::kTH2D, {axisPtRes, axisDeltaphiRes}, true);
    } else {
      registry.add<TH2>("PtGen_DeltaPtOverPtGen", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,e} (GeV/#it{c})"}, axisDeltaptRes}, true);
      registry.add<TH2>("PtGen_DeltaEta", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,e} (GeV/#it{c})"}, axisDeltaetaRes}, true);
      registry.add<TH2>("PtGen_DeltaPhi_Ele", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,e} (GeV/#it{c})"}, axisDeltaphiRes}, true);
      registry.add<TH2>("PtGen_DeltaPhi_Pos", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,e} (GeV/#it{c})"}, axisDeltaphiRes}, true);
    }
  }

  template <typename TTracksMC>
  void Check(TTracksMC const& tracksMC)
  {
    for (auto& mctrack : tracksMC) {
      if (abs(mctrack.pdgCode()) != fPdgCode) {
        continue;
      }

      float deltaptoverpt = -1000.;
      if (mctrack.pt() > 0.)
        deltaptoverpt = (mctrack.pt() - mctrack.ptSmeared()) / mctrack.pt();
      float deltaeta = mctrack.eta() - mctrack.etaSmeared();
      float deltaphi = mctrack.phi() - mctrack.phiSmeared();
      registry.fill(HIST("PtGen_DeltaPtOverPtGen"), mctrack.pt(), deltaptoverpt);
      registry.fill(HIST("PtGen_DeltaEta"), mctrack.pt(), deltaeta);
      if (mctrack.pdgCode() < 0) {
        registry.fill(HIST("PtGen_DeltaPhi_Ele"), mctrack.pt(), deltaphi);
      } else {
        registry.fill(HIST("PtGen_DeltaPhi_Pos"), mctrack.pt(), deltaphi);
      }
      registry.fill(HIST("hCorrelation_Pt"), mctrack.pt(), mctrack.ptSmeared());
      registry.fill(HIST("hCorrelation_Eta"), mctrack.eta(), mctrack.etaSmeared());
      registry.fill(HIST("hCorrelation_Phi"), mctrack.phi(), mctrack.phiSmeared());
    } // end of mctrack loop
  }

  void processCheckMCanalysis(MyReducedTracks const& tracksMC)
  {
    Check(tracksMC);
  }

  void processCheckCocktail(MyCocktailTracks const& tracksMC)
  {
    Check(tracksMC);
  }

  void processDummyMCanalysis(ReducedMCTracks const& tracksMC) {}
  void processDummyCocktail(aod::McParticles const& tracksMC) {}

  PROCESS_SWITCH(CheckSmearing, processCheckMCanalysis, "Run for MC analysis", false);
  PROCESS_SWITCH(CheckSmearing, processCheckCocktail, "Run for cocktail analysis", false);
  PROCESS_SWITCH(CheckSmearing, processDummyMCanalysis, "Dummy process function", false);
  PROCESS_SWITCH(CheckSmearing, processDummyCocktail, "Dummy process function", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ApplySmearing>(cfgc, TaskName{"apply-smearing"}),
    adaptAnalysisTask<CheckSmearing>(cfgc, TaskName{"check-smearing"})};
}
