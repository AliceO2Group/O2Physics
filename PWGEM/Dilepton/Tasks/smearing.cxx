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
// Analysis task to produce smeared pt, eta, phi for electrons/muons in dilepton analysis
//    Please write to: daiki.sekihata@cern.ch

#include <CCDB/BasicCCDBManager.h>
#include <chrono>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h" // remove this later, because 2 data tables (covariant matrix) in this header confilict against EM tables.
#include "PWGEM/Dilepton/Utils/MomentumSmearer.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

struct ApplySmearing {
  Produces<aod::SmearedElectrons> smearedelectron;
  Produces<aod::SmearedMuons> smearedmuon;

  // Maps
  Configurable<std::string> fConfigResFileName_Electron{"cfgResFileName_Electron", "", "name of resolution file"};
  Configurable<std::string> fConfigResPtHistName_Electron{"cfgResPtHistName_Electron", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
  Configurable<std::string> fConfigResEtaHistName_Electron{"cfgResEtaHistName_Electron", "EtaResArr", "histogram name for eta in resolution file"};
  Configurable<std::string> fConfigResPhiPosHistName_Electron{"cfgResPhiPosHistName_Electron", "PhiPosResArr", "histogram name for phi pos in resolution file"};
  Configurable<std::string> fConfigResPhiNegHistName_Electron{"cfgResPhiNegHistName_Electron", "PhiEleResArr", "hisogram for phi neg in resolution file"};
  Configurable<std::string> fConfigEffFileName_Electron{"cfgEffFileName_Electron", "", "name of efficiency file"};
  Configurable<std::string> fConfigEffHistName_Electron{"cfgEffHistName_Electron", "fhwEffpT", "name of efficiency histogram"};
  Configurable<std::string> fConfigDCAFileName_Electron{"cfgDCAFileName_Electron", "", "name of DCA template file"};
  Configurable<std::string> fConfigDCAHistName_Electron{"cfgDCAHistName_Electron", "fh_DCAtempaltes", "histogram name of the DCA templates"};

  Configurable<std::string> fConfigResFileName_StandaloneMuon{"cfgResFileName_StandaloneMuon", "", "name of resolution file"};
  Configurable<std::string> fConfigResPtHistName_StandaloneMuon{"cfgResPtHistName_StandaloneMuon", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
  Configurable<std::string> fConfigResEtaHistName_StandaloneMuon{"cfgResEtaHistName_StandaloneMuon", "EtaResArr", "histogram name for eta in resolution file"};
  Configurable<std::string> fConfigResPhiPosHistName_StandaloneMuon{"cfgResPhiPosHistName_StandaloneMuon", "PhiPosResArr", "histogram name for phi pos in resolution file"};
  Configurable<std::string> fConfigResPhiNegHistName_StandaloneMuon{"cfgResPhiNegHistName_StandaloneMuon", "PhiEleResArr", "hisogram for phi neg in resolution file"};
  Configurable<std::string> fConfigEffFileName_StandaloneMuon{"cfgEffFileName_StandaloneMuon", "", "name of efficiency file"};
  Configurable<std::string> fConfigEffHistName_StandaloneMuon{"cfgEffHistName_StandaloneMuon", "fhwEffpT", "name of efficiency histogram"};
  Configurable<std::string> fConfigDCAFileName_StandaloneMuon{"cfgDCAFileName_StandaloneMuon", "", "name of DCA template file"};
  Configurable<std::string> fConfigDCAHistName_StandaloneMuon{"cfgDCAHistName_StandaloneMuon", "fh_DCAtempaltes", "histogram name of the DCA templates"};

  Configurable<std::string> fConfigResFileName_GlobalMuon{"cfgResFileName_GlobalMuon", "", "name of resolution file"};
  Configurable<std::string> fConfigResPtHistName_GlobalMuon{"cfgResPtHistName_GlobalMuon", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
  Configurable<std::string> fConfigResEtaHistName_GlobalMuon{"cfgResEtaHistName_GlobalMuon", "EtaResArr", "histogram name for eta in resolution file"};
  Configurable<std::string> fConfigResPhiPosHistName_GlobalMuon{"cfgResPhiPosHistName_GlobalMuon", "PhiPosResArr", "histogram name for phi pos in resolution file"};
  Configurable<std::string> fConfigResPhiNegHistName_GlobalMuon{"cfgResPhiNegHistName_GlobalMuon", "PhiEleResArr", "hisogram for phi neg in resolution file"};
  Configurable<std::string> fConfigEffFileName_GlobalMuon{"cfgEffFileName_GlobalMuon", "", "name of efficiency file"};
  Configurable<std::string> fConfigEffHistName_GlobalMuon{"cfgEffHistName_GlobalMuon", "fhwEffpT", "name of efficiency histogram"};
  Configurable<std::string> fConfigDCAFileName_GlobalMuon{"cfgDCAFileName_GlobalMuon", "", "name of DCA template file"};
  Configurable<std::string> fConfigDCAHistName_GlobalMuon{"cfgDCAHistName_GlobalMuon", "fh_DCAtempaltes", "histogram name of the DCA templates"};

  Configurable<bool> fFromCcdb{"cfgFromCcdb", false, "get resolution and efficiency histos from CCDB"};
  Configurable<std::string> fConfigCcdbPathRes{"cfgCcdbPathRes", "", "path to the ccdb object for resolution"};
  Configurable<std::string> fConfigCcdbPathEff{"cfgCcdbPahtEff", "", "path to the ccdb object for efficiency"};
  Configurable<std::string> fConfigCcdbPathDCA{"cfgCcdbPahtDCA", "", "path to the ccdb object for dca"};
  Configurable<std::string> fConfigCcdbUrl{"cfgCcdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> fConfigCcdbNoLaterThan{"cfgCcdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  MomentumSmearer smearer_Electron;
  MomentumSmearer smearer_StandaloneMuon;
  MomentumSmearer smearer_GlobalMuon;
  Service<ccdb::BasicCCDBManager> ccdb;

  void init(InitContext&)
  {
    smearer_Electron.setResFileName(TString(fConfigResFileName_Electron));
    smearer_Electron.setResPtHistName(TString(fConfigResPtHistName_Electron));
    smearer_Electron.setResEtaHistName(TString(fConfigResEtaHistName_Electron));
    smearer_Electron.setResPhiPosHistName(TString(fConfigResPhiPosHistName_Electron));
    smearer_Electron.setResPhiNegHistName(TString(fConfigResPhiNegHistName_Electron));
    smearer_Electron.setEffFileName(TString(fConfigEffFileName_Electron));
    smearer_Electron.setEffHistName(TString(fConfigEffHistName_Electron));
    smearer_Electron.setDCAFileName(TString(fConfigDCAFileName_Electron));
    smearer_Electron.setDCAHistName(TString(fConfigDCAHistName_Electron));

    smearer_StandaloneMuon.setResFileName(TString(fConfigResFileName_StandaloneMuon));
    smearer_StandaloneMuon.setResPtHistName(TString(fConfigResPtHistName_StandaloneMuon));
    smearer_StandaloneMuon.setResEtaHistName(TString(fConfigResEtaHistName_StandaloneMuon));
    smearer_StandaloneMuon.setResPhiPosHistName(TString(fConfigResPhiPosHistName_StandaloneMuon));
    smearer_StandaloneMuon.setResPhiNegHistName(TString(fConfigResPhiNegHistName_StandaloneMuon));
    smearer_StandaloneMuon.setEffFileName(TString(fConfigEffFileName_StandaloneMuon));
    smearer_StandaloneMuon.setEffHistName(TString(fConfigEffHistName_StandaloneMuon));
    smearer_StandaloneMuon.setDCAFileName(TString(fConfigDCAFileName_StandaloneMuon));
    smearer_StandaloneMuon.setDCAHistName(TString(fConfigDCAHistName_StandaloneMuon));

    smearer_GlobalMuon.setResFileName(TString(fConfigResFileName_GlobalMuon));
    smearer_GlobalMuon.setResPtHistName(TString(fConfigResPtHistName_GlobalMuon));
    smearer_GlobalMuon.setResEtaHistName(TString(fConfigResEtaHistName_GlobalMuon));
    smearer_GlobalMuon.setResPhiPosHistName(TString(fConfigResPhiPosHistName_GlobalMuon));
    smearer_GlobalMuon.setResPhiNegHistName(TString(fConfigResPhiNegHistName_GlobalMuon));
    smearer_GlobalMuon.setEffFileName(TString(fConfigEffFileName_GlobalMuon));
    smearer_GlobalMuon.setEffHistName(TString(fConfigEffHistName_GlobalMuon));
    smearer_GlobalMuon.setDCAFileName(TString(fConfigDCAFileName_GlobalMuon));
    smearer_GlobalMuon.setDCAHistName(TString(fConfigDCAHistName_GlobalMuon));

    if (fFromCcdb) {
      ccdb->setURL(fConfigCcdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(fConfigCcdbNoLaterThan);

      smearer_Electron.setCcdbPathRes(TString(fConfigCcdbPathRes));
      smearer_Electron.setCcdbPathEff(TString(fConfigCcdbPathEff));
      smearer_Electron.setCcdbPathDCA(TString(fConfigCcdbPathDCA));
      smearer_Electron.setTimestamp(fConfigCcdbNoLaterThan);
      smearer_Electron.setCcdb(ccdb);

      smearer_StandaloneMuon.setCcdbPathRes(TString(fConfigCcdbPathRes));
      smearer_StandaloneMuon.setCcdbPathEff(TString(fConfigCcdbPathEff));
      smearer_StandaloneMuon.setCcdbPathDCA(TString(fConfigCcdbPathDCA));
      smearer_StandaloneMuon.setTimestamp(fConfigCcdbNoLaterThan);
      smearer_StandaloneMuon.setCcdb(ccdb);

      smearer_GlobalMuon.setCcdbPathRes(TString(fConfigCcdbPathRes));
      smearer_GlobalMuon.setCcdbPathEff(TString(fConfigCcdbPathEff));
      smearer_GlobalMuon.setCcdbPathDCA(TString(fConfigCcdbPathDCA));
      smearer_GlobalMuon.setTimestamp(fConfigCcdbNoLaterThan);
      smearer_GlobalMuon.setCcdb(ccdb);
    }
    smearer_Electron.init();
    smearer_StandaloneMuon.init();
    smearer_GlobalMuon.init();
  }

  template <typename TTracksMC>
  void applySmearing(TTracksMC const& tracksMC)
  {
    for (auto& mctrack : tracksMC) {
      float ptgen = mctrack.pt();
      float etagen = mctrack.eta();
      float phigen = mctrack.phi();
      float efficiency = 1.;
      float dca = 0.;

      int pdgCode = mctrack.pdgCode();
      if (abs(pdgCode) == 11) {
        int ch = -1;
        if (pdgCode < 0) {
          ch = 1;
        }
        // apply smearing for electrons or muons.
        float ptsmeared, etasmeared, phismeared;
        smearer_Electron.applySmearing(ch, ptgen, etagen, phigen, ptsmeared, etasmeared, phismeared);
        // get the efficiency
        efficiency = smearer_Electron.getEfficiency(ptgen, etagen, phigen);
        // get DCA
        dca = smearer_Electron.getDCA(ptsmeared);
        // fill the table
        smearedelectron(ptsmeared, etasmeared, phismeared, efficiency, dca);
        smearedmuon(ptgen, etagen, phigen, 1.f, 0.f, ptgen, etagen, phigen, 1.f, 0.f);
      } else if (abs(pdgCode) == 13) {
        int ch = -1;
        if (pdgCode < 0) {
          ch = 1;
        }
        // apply smearing for muons based on resolution map of standalone muons
        float ptsmeared_sa = 0.f, etasmeared_sa = 0.f, phismeared_sa = 0.f, efficiency_sa = 1.f, dca_sa = 0.f;
        smearer_StandaloneMuon.applySmearing(ch, ptgen, etagen, phigen, ptsmeared_sa, etasmeared_sa, phismeared_sa);
        efficiency_sa = smearer_StandaloneMuon.getEfficiency(ptgen, etagen, phigen);
        dca_sa = smearer_StandaloneMuon.getDCA(ptsmeared_sa);

        float ptsmeared_gl = 0.f, etasmeared_gl = 0.f, phismeared_gl = 0.f, efficiency_gl = 1.f, dca_gl = 0.f;
        // apply smearing for muons based on resolution map of global muons
        smearer_GlobalMuon.applySmearing(ch, ptgen, etagen, phigen, ptsmeared_gl, etasmeared_gl, phismeared_gl);
        efficiency_gl = smearer_GlobalMuon.getEfficiency(ptgen, etagen, phigen);
        dca_gl = smearer_GlobalMuon.getDCA(ptsmeared_gl);
        smearedmuon(ptsmeared_sa, etasmeared_sa, phismeared_sa, efficiency_sa, dca_sa, ptsmeared_gl, etasmeared_gl, phismeared_gl, efficiency_gl, dca_gl);

        smearedelectron(ptgen, etagen, phigen, 1.f, 0.f);
      } else {
        // don't apply smearing
        smearedelectron(ptgen, etagen, phigen, efficiency, dca);
        smearedmuon(ptgen, etagen, phigen, efficiency, dca, ptgen, etagen, phigen, efficiency, dca);
      }
    }
  }

  void processMCanalysisEM(aod::EMMCParticles const& tracksMC)
  {
    applySmearing(tracksMC);
  }

  void processMCanalysisDQ(ReducedMCTracks const& tracksMC)
  {
    applySmearing(tracksMC);
  }

  void processCocktail(aod::McParticles const& tracksMC)
  {
    applySmearing(tracksMC);
  }

  void processDummyCocktail(aod::McParticles const&) {}

  void processDummyMCanalysis(ReducedMCTracks const&) {}

  PROCESS_SWITCH(ApplySmearing, processMCanalysisEM, "Run for MC analysis", false);
  PROCESS_SWITCH(ApplySmearing, processMCanalysisDQ, "Run for MC analysis", false);
  PROCESS_SWITCH(ApplySmearing, processCocktail, "Run for cocktail analysis", false);
  PROCESS_SWITCH(ApplySmearing, processDummyMCanalysis, "Dummy process function", false);
  PROCESS_SWITCH(ApplySmearing, processDummyCocktail, "Dummy process function", true);
};

struct CheckSmearing {
  using EMMCParticlesWithSmearing = soa::Join<aod::EMMCParticles, aod::SmearedElectrons>; // this is only for electrons
  using MyReducedTracks = soa::Join<ReducedMCTracks, aod::SmearedElectrons>;              // this is only for electrons
  using MyCocktailTracks = soa::Join<aod::McParticles, aod::SmearedElectrons>;            // this is only for electrons

  // Run for electrons or muons
  Configurable<int> fPdgCode{"cfgPdgCode", 11, "Set the type of particle to be checked"};

  // Resolution histos as cross check
  Configurable<bool> fConfigUsePtVecRes{"cfgUsePtVecRes", true, "If true, non-linear pt bins predefined in res histos"};
  ConfigurableAxis ptResBinsVec{"ptResBinsVec", {0., 0., 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.20, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285, 0.29, 0.295, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80, 2.00, 2.40, 2.80, 3.20, 3.70, 4.50, 6.00, 8.00, 10., 12.0, 14., 16., 18., 20.}, "Pt binning vector for resolution"};
  ConfigurableAxis ptResBins{"ptResBins", {20, 0.f, 20.f}, "Pt binning for resolution"};
  ConfigurableAxis deltaptResBins{"deltaptResBins", {500, -1.f, 1.f}, "DeltaPt binning for resolution"};
  ConfigurableAxis deltaetaResBins{"deltaetaResBins", {500, -0.5f, 0.5f}, "DeltaEta binning for resolution"};
  ConfigurableAxis deltaphiResBins{"deltaphiResBins", {500, -0.5f, 0.5f}, "DeltaPhi binning for resolution"};

  HistogramRegistry registry{"output", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    registry.add<TH2>("hCorrelation_Pt", "pT correlation", {HistType::kTH2F, {{1000, 0.0f, 10.0f}, {1000, 0.0f, 10.0f}}});
    registry.add<TH2>("hCorrelation_Eta", "eta correlation", {HistType::kTH2F, {{200, -1.0f, +1.0f}, {200, -1.0f, +1.0f}}});
    registry.add<TH2>("hCorrelation_Phi", "phi correlation", {HistType::kTH2F, {{100, 0.0f, TMath::TwoPi()}, {100, 0.0f, TMath::TwoPi()}}});

    // Binning for resolution
    AxisSpec axisPtRes{ptResBins, "#it{p}^{gen}_{T,l} (GeV/#it{c})"};
    AxisSpec axisDeltaptRes{deltaptResBins, "(p^{gen}_{T} - p^{rec}_{T}) / p^{gen}_{T}"};
    AxisSpec axisDeltaetaRes{deltaetaResBins, "#eta^{gen} - #eta^{rec}"};
    AxisSpec axisDeltaphiRes{deltaphiResBins, "#varphi^{gen} - #varphi^{rec} (rad.)"};

    if (!fConfigUsePtVecRes) {
      registry.add<TH2>("PtGen_DeltaPtOverPtGen", "", HistType::kTH2D, {axisPtRes, axisDeltaptRes}, true);
      registry.add<TH2>("PtGen_DeltaEta", "", HistType::kTH2D, {axisPtRes, axisDeltaetaRes}, true);
      registry.add<TH2>("PtGen_DeltaPhi_Neg", "", HistType::kTH2D, {axisPtRes, axisDeltaphiRes}, true);
      registry.add<TH2>("PtGen_DeltaPhi_Pos", "", HistType::kTH2D, {axisPtRes, axisDeltaphiRes}, true);
    } else {
      registry.add<TH2>("PtGen_DeltaPtOverPtGen", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,l} (GeV/#it{c})"}, axisDeltaptRes}, true);
      registry.add<TH2>("PtGen_DeltaEta", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,l} (GeV/#it{c})"}, axisDeltaetaRes}, true);
      registry.add<TH2>("PtGen_DeltaPhi_Neg", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,l} (GeV/#it{c})"}, axisDeltaphiRes}, true);
      registry.add<TH2>("PtGen_DeltaPhi_Pos", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,l} (GeV/#it{c})"}, axisDeltaphiRes}, true);
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
        registry.fill(HIST("PtGen_DeltaPhi_Neg"), mctrack.pt(), deltaphi);
      } else {
        registry.fill(HIST("PtGen_DeltaPhi_Pos"), mctrack.pt(), deltaphi);
      }
      registry.fill(HIST("hCorrelation_Pt"), mctrack.pt(), mctrack.ptSmeared());
      registry.fill(HIST("hCorrelation_Eta"), mctrack.eta(), mctrack.etaSmeared());
      registry.fill(HIST("hCorrelation_Phi"), mctrack.phi(), mctrack.phiSmeared());
    } // end of mctrack loop
  }

  void processCheckMCanalysisEM(EMMCParticlesWithSmearing const& tracksMC)
  {
    Check(tracksMC);
  }

  void processCheckMCanalysisDQ(MyReducedTracks const& tracksMC)
  {
    Check(tracksMC);
  }

  void processCheckCocktail(MyCocktailTracks const& tracksMC)
  {
    Check(tracksMC);
  }

  void processDummyMCanalysis(ReducedMCTracks const&) {}
  void processDummyCocktail(aod::McParticles const&) {}

  PROCESS_SWITCH(CheckSmearing, processCheckMCanalysisEM, "Run for MC analysis", false);
  PROCESS_SWITCH(CheckSmearing, processCheckMCanalysisDQ, "Run for MC analysis", false);
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
