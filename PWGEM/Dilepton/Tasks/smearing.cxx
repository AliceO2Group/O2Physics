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

#include <array>
#include <string>
#include <chrono>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
// #include "PWGDQ/DataModel/ReducedInfoTables.h" // remove this later, because 2 data tables (covariant matrix) in this header confilict against EM tables.
#include "PWGEM/Dilepton/Utils/MomentumSmearer.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

namespace o2::aod::pwgem::dilepton::smearing
{
enum class EMAnaType : int {
  kEfficiency = 0,
  kCocktail = 1,
};
} // namespace o2::aod::pwgem::dilepton::smearing

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsCent, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyMCCollisions = soa::Join<aod::EMMCEvents, aod::MostProbableEMEventIdsInMC>;
using MyMCCollision = MyMCCollisions::iterator;

struct ApplySmearing {

  Produces<aod::SmearedElectrons> smearedelectron;
  Produces<aod::SmearedMuons> smearedmuon;

  Configurable<bool> fFromCcdb{"cfgFromCcdb", false, "get resolution and efficiency histos from CCDB"};
  Configurable<std::string> fConfigCcdbUrl{"cfgCcdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> fTimestamp{"cfgCcdbTimestamp", 10, "valid timestamp of CCDB object"};
  Configurable<float> fCentralityForCocktail{"cfgCentralityForCocktail", 5, "average centrality for cocktail"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};

  struct : ConfigurableGroup {
    std::string prefix = "electron_filename_group";
    Configurable<bool> fConfigNDSmearing{"cfgNDSmearing", false, "apply ND-correlated smearing"};
    Configurable<std::string> fConfigResFileName{"cfgResFileName", "", "name of resolution file"};
    Configurable<std::string> fConfigResNDHistName{"cfgResNDHistName", "hs_reso", "name of ND resolution file"};
    Configurable<std::string> fConfigResPtHistName{"cfgResPtHistName", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
    Configurable<std::string> fConfigResEtaHistName{"cfgResEtaHistName", "EtaResArr", "histogram name for eta in resolution file"};
    Configurable<std::string> fConfigResPhiPosHistName{"cfgResPhiPosHistName", "PhiPosResArr", "histogram name for phi pos in resolution file"};
    Configurable<std::string> fConfigResPhiNegHistName{"cfgResPhiNegHistName", "PhiEleResArr", "hisogram for phi neg in resolution file"};
    Configurable<std::string> fConfigEffFileName{"cfgEffFileName", "", "name of efficiency file"};
    Configurable<std::string> fConfigEffHistName{"cfgEffHistName", "fhwEffpT", "name of efficiency histogram"};
    Configurable<std::string> fConfigDCAFileName{"cfgDCAFileName", "", "name of DCA template file"};
    Configurable<std::string> fConfigDCAHistName{"cfgDCAHistName", "fh_DCAtempaltes", "histogram name of the DCA templates"};
    Configurable<std::string> fConfigCcdbPathRes{"cfgCcdbPathRes", "", "path to the ccdb object for resolution"};
    Configurable<std::string> fConfigCcdbPathEff{"cfgCcdbPahtEff", "", "path to the ccdb object for efficiency"};
    Configurable<std::string> fConfigCcdbPathDCA{"cfgCcdbPahtDCA", "", "path to the ccdb object for dca"};
    Configurable<float> fConfigMinPt{"cfgMinPt", -1, "if ptgen is smaller than this threshold, this value is used as input for ptgen."};
  } electron_filenames;

  struct : ConfigurableGroup {
    std::string prefix = "sa_muon_filename_group";
    Configurable<bool> fConfigNDSmearing{"cfgNDSmearing", false, "apply ND-correlated smearing"};
    Configurable<std::string> fConfigResFileName{"cfgResFileName", "", "name of resolution file"};
    Configurable<std::string> fConfigResNDHistName{"cfgResNDHistName", "hs_reso", "name of ND resolution file"};
    Configurable<std::string> fConfigResPtHistName{"cfgResPtHistName", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
    Configurable<std::string> fConfigResEtaHistName{"cfgResEtaHistName", "EtaResArr", "histogram name for eta in resolution file"};
    Configurable<std::string> fConfigResPhiPosHistName{"cfgResPhiPosHistName", "PhiPosResArr", "histogram name for phi pos in resolution file"};
    Configurable<std::string> fConfigResPhiNegHistName{"cfgResPhiNegHistName", "PhiEleResArr", "hisogram for phi neg in resolution file"};
    Configurable<std::string> fConfigEffFileName{"cfgEffFileName", "", "name of efficiency file"};
    Configurable<std::string> fConfigEffHistName{"cfgEffHistName", "fhwEffpT", "name of efficiency histogram"};
    Configurable<std::string> fConfigDCAFileName{"cfgDCAFileName", "", "name of DCA template file"};
    Configurable<std::string> fConfigDCAHistName{"cfgDCAHistName", "fh_DCAtempaltes", "histogram name of the DCA templates"};
    Configurable<std::string> fConfigCcdbPathRes{"cfgCcdbPathRes", "", "path to the ccdb object for resolution"};
    Configurable<std::string> fConfigCcdbPathEff{"cfgCcdbPahtEff", "", "path to the ccdb object for efficiency"};
    Configurable<std::string> fConfigCcdbPathDCA{"cfgCcdbPahtDCA", "", "path to the ccdb object for dca"};
    Configurable<float> fConfigMinPt{"cfgMinPt", -1, "if ptgen is smaller than this threshold, this value is used as input for ptgen."};
  } sa_muon_filenames;

  struct : ConfigurableGroup {
    std::string prefix = "gl_muon_filename_group";
    Configurable<bool> fConfigNDSmearing{"cfgNDSmearing", false, "apply ND-correlated smearing"};
    Configurable<std::string> fConfigResFileName{"cfgResFileName", "", "name of resolution file"};
    Configurable<std::string> fConfigResNDHistName{"cfgResNDHistName", "hs_reso", "name of ND resolution file"};
    Configurable<std::string> fConfigResPtHistName{"cfgResPtHistName", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
    Configurable<std::string> fConfigResEtaHistName{"cfgResEtaHistName", "EtaResArr", "histogram name for eta in resolution file"};
    Configurable<std::string> fConfigResPhiPosHistName{"cfgResPhiPosHistName", "PhiPosResArr", "histogram name for phi pos in resolution file"};
    Configurable<std::string> fConfigResPhiNegHistName{"cfgResPhiNegHistName", "PhiEleResArr", "hisogram for phi neg in resolution file"};
    Configurable<std::string> fConfigEffFileName{"cfgEffFileName", "", "name of efficiency file"};
    Configurable<std::string> fConfigEffHistName{"cfgEffHistName", "fhwEffpT", "name of efficiency histogram"};
    Configurable<std::string> fConfigDCAFileName{"cfgDCAFileName", "", "name of DCA template file"};
    Configurable<std::string> fConfigDCAHistName{"cfgDCAHistName", "fh_DCAtempaltes", "histogram name of the DCA templates"};
    Configurable<std::string> fConfigCcdbPathRes{"cfgCcdbPathRes", "", "path to the ccdb object for resolution"};
    Configurable<std::string> fConfigCcdbPathEff{"cfgCcdbPahtEff", "", "path to the ccdb object for efficiency"};
    Configurable<std::string> fConfigCcdbPathDCA{"cfgCcdbPahtDCA", "", "path to the ccdb object for dca"};
    Configurable<float> fConfigMinPt{"cfgMinPt", -1, "if ptgen is smaller than this threshold, this value is used as input for ptgen."};
  } gl_muon_filenames;

  MomentumSmearer smearer_Electron;
  MomentumSmearer smearer_StandaloneMuon;
  MomentumSmearer smearer_GlobalMuon;
  Service<ccdb::BasicCCDBManager> ccdb;

  void init(InitContext&)
  {
    smearer_Electron.setNDSmearing(electron_filenames.fConfigNDSmearing.value);
    smearer_Electron.setResFileName(TString(electron_filenames.fConfigResFileName));
    smearer_Electron.setResNDHistName(TString(electron_filenames.fConfigResNDHistName));
    smearer_Electron.setResPtHistName(TString(electron_filenames.fConfigResPtHistName));
    smearer_Electron.setResEtaHistName(TString(electron_filenames.fConfigResEtaHistName));
    smearer_Electron.setResPhiPosHistName(TString(electron_filenames.fConfigResPhiPosHistName));
    smearer_Electron.setResPhiNegHistName(TString(electron_filenames.fConfigResPhiNegHistName));
    smearer_Electron.setEffFileName(TString(electron_filenames.fConfigEffFileName));
    smearer_Electron.setEffHistName(TString(electron_filenames.fConfigEffHistName));
    smearer_Electron.setDCAFileName(TString(electron_filenames.fConfigDCAFileName));
    smearer_Electron.setDCAHistName(TString(electron_filenames.fConfigDCAHistName));
    smearer_Electron.setMinPt(electron_filenames.fConfigMinPt);

    smearer_StandaloneMuon.setNDSmearing(sa_muon_filenames.fConfigNDSmearing.value);
    smearer_StandaloneMuon.setResFileName(TString(sa_muon_filenames.fConfigResFileName));
    smearer_StandaloneMuon.setResNDHistName(TString(sa_muon_filenames.fConfigResNDHistName));
    smearer_StandaloneMuon.setResPtHistName(TString(sa_muon_filenames.fConfigResPtHistName));
    smearer_StandaloneMuon.setResEtaHistName(TString(sa_muon_filenames.fConfigResEtaHistName));
    smearer_StandaloneMuon.setResPhiPosHistName(TString(sa_muon_filenames.fConfigResPhiPosHistName));
    smearer_StandaloneMuon.setResPhiNegHistName(TString(sa_muon_filenames.fConfigResPhiNegHistName));
    smearer_StandaloneMuon.setEffFileName(TString(sa_muon_filenames.fConfigEffFileName));
    smearer_StandaloneMuon.setEffHistName(TString(sa_muon_filenames.fConfigEffHistName));
    smearer_StandaloneMuon.setDCAFileName(TString(sa_muon_filenames.fConfigDCAFileName));
    smearer_StandaloneMuon.setDCAHistName(TString(sa_muon_filenames.fConfigDCAHistName));
    smearer_StandaloneMuon.setMinPt(sa_muon_filenames.fConfigMinPt);

    smearer_GlobalMuon.setNDSmearing(gl_muon_filenames.fConfigNDSmearing.value);
    smearer_GlobalMuon.setResFileName(TString(gl_muon_filenames.fConfigResFileName));
    smearer_GlobalMuon.setResNDHistName(TString(gl_muon_filenames.fConfigResNDHistName));
    smearer_GlobalMuon.setResPtHistName(TString(gl_muon_filenames.fConfigResPtHistName));
    smearer_GlobalMuon.setResEtaHistName(TString(gl_muon_filenames.fConfigResEtaHistName));
    smearer_GlobalMuon.setResPhiPosHistName(TString(gl_muon_filenames.fConfigResPhiPosHistName));
    smearer_GlobalMuon.setResPhiNegHistName(TString(gl_muon_filenames.fConfigResPhiNegHistName));
    smearer_GlobalMuon.setEffFileName(TString(gl_muon_filenames.fConfigEffFileName));
    smearer_GlobalMuon.setEffHistName(TString(gl_muon_filenames.fConfigEffHistName));
    smearer_GlobalMuon.setDCAFileName(TString(gl_muon_filenames.fConfigDCAFileName));
    smearer_GlobalMuon.setDCAHistName(TString(gl_muon_filenames.fConfigDCAHistName));
    smearer_GlobalMuon.setMinPt(gl_muon_filenames.fConfigMinPt);

    if (fFromCcdb) {
      ccdb->setURL(fConfigCcdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()); // now

      smearer_Electron.setCcdbPathRes(TString(electron_filenames.fConfigCcdbPathRes));
      smearer_Electron.setCcdbPathEff(TString(electron_filenames.fConfigCcdbPathEff));
      smearer_Electron.setCcdbPathDCA(TString(electron_filenames.fConfigCcdbPathDCA));
      smearer_Electron.setTimestamp(fTimestamp);
      smearer_Electron.setCcdb(ccdb);

      smearer_StandaloneMuon.setCcdbPathRes(TString(sa_muon_filenames.fConfigCcdbPathRes));
      smearer_StandaloneMuon.setCcdbPathEff(TString(sa_muon_filenames.fConfigCcdbPathEff));
      smearer_StandaloneMuon.setCcdbPathDCA(TString(sa_muon_filenames.fConfigCcdbPathDCA));
      smearer_StandaloneMuon.setTimestamp(fTimestamp);
      smearer_StandaloneMuon.setCcdb(ccdb);

      smearer_GlobalMuon.setCcdbPathRes(TString(gl_muon_filenames.fConfigCcdbPathRes));
      smearer_GlobalMuon.setCcdbPathEff(TString(gl_muon_filenames.fConfigCcdbPathEff));
      smearer_GlobalMuon.setCcdbPathDCA(TString(gl_muon_filenames.fConfigCcdbPathDCA));
      smearer_GlobalMuon.setTimestamp(fTimestamp);
      smearer_GlobalMuon.setCcdb(ccdb);
    }
    smearer_Electron.init();
    smearer_StandaloneMuon.init();
    smearer_GlobalMuon.init();
  }

  template <o2::aod::pwgem::dilepton::smearing::EMAnaType type, typename TTracksMC, typename TCollisions, typename TMCCollisions>
  void applySmearing(TTracksMC const& tracksMC, TCollisions const& collisions, TMCCollisions const&)
  {
    for (auto& mctrack : tracksMC) {
      float ptgen = mctrack.pt();
      float etagen = mctrack.eta();
      float phigen = mctrack.phi();
      float efficiency = 1.;
      float dca = 0.;

      float ptsmeared = 0, etasmeared = 0, phismeared = 0;
      float centrality = -1.f;
      if constexpr (type == o2::aod::pwgem::dilepton::smearing::EMAnaType::kEfficiency) {
        auto mccollision = mctrack.template emmcevent_as<TMCCollisions>();
        if (mccollision.mpemeventId() > 0) { // if mc collisions are not reconstructed, such mc collisions should not enter efficiency calculation.
          auto collision = collisions.rawIteratorAt(mccollision.mpemeventId());
          centrality = std::array{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()}[cfgCentEstimator];
        } else {
          ptsmeared = ptgen;
          etasmeared = etagen;
          phismeared = phigen;
        }
      } else {
        centrality = fCentralityForCocktail;
      }

      int pdgCode = mctrack.pdgCode();
      if (abs(pdgCode) == 11) {
        int ch = -1;
        if (pdgCode < 0) {
          ch = 1;
        }
        // apply smearing for electrons or muons.
        smearer_Electron.applySmearing(centrality, ch, ptgen, etagen, phigen, ptsmeared, etasmeared, phismeared);
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
        smearer_StandaloneMuon.applySmearing(centrality, ch, ptgen, etagen, phigen, ptsmeared_sa, etasmeared_sa, phismeared_sa);
        efficiency_sa = smearer_StandaloneMuon.getEfficiency(ptgen, etagen, phigen);
        dca_sa = smearer_StandaloneMuon.getDCA(ptsmeared_sa);

        float ptsmeared_gl = 0.f, etasmeared_gl = 0.f, phismeared_gl = 0.f, efficiency_gl = 1.f, dca_gl = 0.f;
        // apply smearing for muons based on resolution map of global muons
        smearer_GlobalMuon.applySmearing(centrality, ch, ptgen, etagen, phigen, ptsmeared_gl, etasmeared_gl, phismeared_gl);
        efficiency_gl = smearer_GlobalMuon.getEfficiency(ptgen, etagen, phigen);
        dca_gl = smearer_GlobalMuon.getDCA(ptsmeared_gl);
        smearedmuon(ptsmeared_sa, etasmeared_sa, phismeared_sa, efficiency_sa, dca_sa, ptsmeared_gl, etasmeared_gl, phismeared_gl, efficiency_gl, dca_gl);

        smearedelectron(ptgen, etagen, phigen, 1.f, 0.f);
      } else {
        // don't apply smearing
        smearedelectron(ptgen, etagen, phigen, efficiency, dca);
        smearedmuon(ptgen, etagen, phigen, efficiency, dca, ptgen, etagen, phigen, efficiency, dca);
      }
    } // end of mc track loop
  }

  void processMCanalysisEM(aod::EMMCParticles const& tracksMC, MyCollisions const& collisions, MyMCCollisions const& mccollisions)
  {
    applySmearing<o2::aod::pwgem::dilepton::smearing::EMAnaType::kEfficiency>(tracksMC, collisions, mccollisions);
  }

  // void processMCanalysisDQ(ReducedMCTracks const& tracksMC)
  // {
  //   applySmearing<EMAnaType::kEfficiency>(tracksMC);
  // }

  void processCocktail(aod::McParticles const& tracksMC)
  {
    applySmearing<o2::aod::pwgem::dilepton::smearing::EMAnaType::kCocktail>(tracksMC, nullptr, nullptr);
  }

  void processDummyCocktail(aod::McParticles const& tracksMC)
  {
    // don't apply smearing
    for (auto& mctrack : tracksMC) {
      int pdgCode = mctrack.pdgCode();
      if (abs(pdgCode) == 11) {
        smearedelectron(mctrack.pt(), mctrack.eta(), mctrack.phi(), 1.0, 0.0);
      } else if (abs(pdgCode) == 13) {
        smearedmuon(mctrack.pt(), mctrack.eta(), mctrack.phi(), 1.0, 0.0, mctrack.pt(), mctrack.eta(), mctrack.phi(), 1.0, 0.0);
      } else {
        smearedelectron(mctrack.pt(), mctrack.eta(), mctrack.eta(), 1.0, 0.0);
        smearedmuon(mctrack.pt(), mctrack.eta(), mctrack.phi(), 1.0, 0.0, mctrack.pt(), mctrack.eta(), mctrack.phi(), 1.0, 0.0);
      }
    }
  }

  void processDummyMCanalysisEM(aod::EMMCParticles const&) {}
  // void processDummyMCanalysisDQ(ReducedMCTracks const&) {}

  PROCESS_SWITCH(ApplySmearing, processMCanalysisEM, "Run for MC analysis which uses skimmed EM data format", false);
  // PROCESS_SWITCH(ApplySmearing, processMCanalysisDQ, "Run for MC analysis which uses skimmed DQ data format", false);
  PROCESS_SWITCH(ApplySmearing, processCocktail, "Run for cocktail analysis", false);
  PROCESS_SWITCH(ApplySmearing, processDummyMCanalysisEM, "Dummy process function", false);
  // PROCESS_SWITCH(ApplySmearing, processDummyMCanalysisDQ, "Dummy process function", false);
  PROCESS_SWITCH(ApplySmearing, processDummyCocktail, "Dummy process function", true);
};

struct CheckSmearing {
  using EMMCParticlesWithSmearing = soa::Join<aod::EMMCParticles, aod::SmearedElectrons>; // this is only for electrons
  // using MyReducedTracks = soa::Join<ReducedMCTracks, aod::SmearedElectrons>;              // this is only for electrons
  using MyCocktailTracks = soa::Join<aod::McParticles, aod::SmearedElectrons>; // this is only for electrons

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
    registry.add<TH2>("hCorrelation_Pt", "pT correlation;p_{T,l}^{gen} (GeV/c);p_{T,l}^{smeared} (GeV/c)", {HistType::kTH2F, {{1000, 0.0f, 10.0f}, {1000, 0.0f, 10.0f}}});
    registry.add<TH2>("hCorrelation_Eta", "eta correlation;#eta_{l}^{gen};#eta_{l}^{smeared}", {HistType::kTH2F, {{200, -1.0f, +1.0f}, {200, -1.0f, +1.0f}}});
    registry.add<TH2>("hCorrelation_Phi", "phi correlation;#varphi_{l}^{gen} (rad.);#varphi_{l}^{smeared} (rad.)", {HistType::kTH2F, {{100, 0.0f, TMath::TwoPi()}, {100, 0.0f, TMath::TwoPi()}}});

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

  template <o2::aod::pwgem::dilepton::smearing::EMAnaType type, typename TTracksMC, typename TMCCollisions>
  void Check(TTracksMC const& tracksMC, TMCCollisions const&)
  {
    for (auto& mctrack : tracksMC) {
      if (abs(mctrack.pdgCode()) != fPdgCode) {
        continue;
      }

      if constexpr (type == o2::aod::pwgem::dilepton::smearing::EMAnaType::kEfficiency) {
        auto mccollision = mctrack.template emmcevent_as<TMCCollisions>();
        if (mccollision.mpemeventId() < 0) { // if mc collisions are not reconstructed, such mc collisions should not enter efficiency calculation.
          continue;
        }
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

  void processCheckMCanalysisEM(EMMCParticlesWithSmearing const& tracksMC, MyMCCollisions const& mccollisions)
  {
    Check<o2::aod::pwgem::dilepton::smearing::EMAnaType::kEfficiency>(tracksMC, mccollisions);
  }

  // void processCheckMCanalysisDQ(MyReducedTracks const& tracksMC)
  // {
  //   Check(tracksMC);
  // }

  void processCheckCocktail(MyCocktailTracks const& tracksMC)
  {
    Check<o2::aod::pwgem::dilepton::smearing::EMAnaType::kCocktail>(tracksMC, nullptr);
  }

  void processDummyMCanalysisEM(aod::EMMCParticles const&) {}
  // void processDummyMCanalysisDQ(ReducedMCTracks const&) {}
  void processDummyCocktail(aod::McParticles const&) {}

  PROCESS_SWITCH(CheckSmearing, processCheckMCanalysisEM, "Run for MC analysis", false);
  // PROCESS_SWITCH(CheckSmearing, processCheckMCanalysisDQ, "Run for MC analysis", false);
  PROCESS_SWITCH(CheckSmearing, processCheckCocktail, "Run for cocktail analysis", false);
  PROCESS_SWITCH(CheckSmearing, processDummyMCanalysisEM, "Dummy process function", false);
  // PROCESS_SWITCH(CheckSmearing, processDummyMCanalysisDQ, "Dummy process function", false);
  PROCESS_SWITCH(CheckSmearing, processDummyCocktail, "Dummy process function", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ApplySmearing>(cfgc, TaskName{"apply-smearing"}),
    adaptAnalysisTask<CheckSmearing>(cfgc, TaskName{"check-smearing"})};
}
