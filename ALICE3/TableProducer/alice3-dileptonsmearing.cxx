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

#include "PWGEM/Dilepton/Utils/MomentumSmearer.h"

#include "ALICE3/DataModel/prefilterDilepton.h"
#include "ALICE3/DataModel/tracksAlice3.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TString.h>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyTracks = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksAlice3>;
using MyTracksWithSmearing = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksAlice3, aod::SmearedAlice3Dilepton>;

struct alice3dileptonsmearer {

  Produces<aod::SmearedAlice3Dilepton> smearedelectron;

  Configurable<bool> cfgFromCcdb{"cfgFromCcdb", false, "get resolution and efficiency histos from CCDB"};
  Configurable<std::string> cfgCcdbUrl{"cfgCcdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> cfgCcdbTimestamp{"cfgCcdbTimestamp", 10, "valid timestamp of CCDB object"};

  struct : ConfigurableGroup {
    std::string prefix = "electron_filename_group";
    Configurable<bool> cfgNDSmearing{"cfgNDSmearing", false, "apply ND-correlated smearing"};
    Configurable<std::string> cfgResFileName{"cfgResFileName", "", "name of resolution file"};
    Configurable<std::string> cfgResNDHistName{"cfgResNDHistName", "hs_reso", "name of ND resolution file"};
    Configurable<std::string> cfgResPtHistName{"cfgResPtHistName", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
    Configurable<std::string> cfgResEtaHistName{"cfgResEtaHistName", "EtaResArr", "histogram name for eta in resolution file"};
    Configurable<std::string> cfgResPhiPosHistName{"cfgResPhiPosHistName", "PhiPosResArr", "histogram name for phi pos in resolution file"};
    Configurable<std::string> cfgResPhiNegHistName{"cfgResPhiNegHistName", "PhiEleResArr", "hisogram for phi neg in resolution file"};
    Configurable<std::string> cfgEffFileName{"cfgEffFileName", "", "name of efficiency file"};
    Configurable<std::string> cfgEffHistName{"cfgEffHistName", "fhwEffpT", "name of efficiency histogram"};
    Configurable<std::string> cfgCcdbPathRes{"cfgCcdbPathRes", "", "path to the ccdb object for resolution"};
    Configurable<std::string> cfgCcdbPathEff{"cfgCcdbPathEff", "", "path to the ccdb object for efficiency"};
    Configurable<float> cfgMinPt{"cfgMinPt", -1, "if ptgen is smaller than this threshold, this value is used as input for ptgen."};
  } electron_filenames;

  struct : ConfigurableGroup {
    std::string prefix = "pion_filename_group";
    Configurable<bool> cfgNDSmearing{"cfgNDSmearing", false, "apply ND-correlated smearing"};
    Configurable<std::string> cfgResFileName{"cfgResFileName", "", "name of resolution file"};
    Configurable<std::string> cfgResNDHistName{"cfgResNDHistName", "hs_reso", "name of ND resolution file"};
    Configurable<std::string> cfgResPtHistName{"cfgResPtHistName", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
    Configurable<std::string> cfgResEtaHistName{"cfgResEtaHistName", "EtaResArr", "histogram name for eta in resolution file"};
    Configurable<std::string> cfgResPhiPosHistName{"cfgResPhiPosHistName", "PhiPosResArr", "histogram name for phi pos in resolution file"};
    Configurable<std::string> cfgResPhiNegHistName{"cfgResPhiNegHistName", "PhiEleResArr", "hisogram for phi neg in resolution file"};
    Configurable<std::string> cfgEffFileName{"cfgEffFileName", "", "name of efficiency file"};
    Configurable<std::string> cfgEffHistName{"cfgEffHistName", "fhwEffpT", "name of efficiency histogram"};
    Configurable<std::string> cfgCcdbPathRes{"cfgCcdbPathRes", "", "path to the ccdb object for resolution"};
    Configurable<std::string> cfgCcdbPathEff{"cfgCcdbPathEff", "", "path to the ccdb object for efficiency"};
    Configurable<float> cfgMinPt{"cfgMinPt", -1, "if ptgen is smaller than this threshold, this value is used as input for ptgen."};
  } pion_filenames;

  struct : ConfigurableGroup {
    std::string prefix = "kaon_filename_group";
    Configurable<bool> cfgNDSmearing{"cfgNDSmearing", false, "apply ND-correlated smearing"};
    Configurable<std::string> cfgResFileName{"cfgResFileName", "", "name of resolution file"};
    Configurable<std::string> cfgResNDHistName{"cfgResNDHistName", "hs_reso", "name of ND resolution file"};
    Configurable<std::string> cfgResPtHistName{"cfgResPtHistName", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
    Configurable<std::string> cfgResEtaHistName{"cfgResEtaHistName", "EtaResArr", "histogram name for eta in resolution file"};
    Configurable<std::string> cfgResPhiPosHistName{"cfgResPhiPosHistName", "PhiPosResArr", "histogram name for phi pos in resolution file"};
    Configurable<std::string> cfgResPhiNegHistName{"cfgResPhiNegHistName", "PhiEleResArr", "hisogram for phi neg in resolution file"};
    Configurable<std::string> cfgEffFileName{"cfgEffFileName", "", "name of efficiency file"};
    Configurable<std::string> cfgEffHistName{"cfgEffHistName", "fhwEffpT", "name of efficiency histogram"};
    Configurable<std::string> cfgCcdbPathRes{"cfgCcdbPathRes", "", "path to the ccdb object for resolution"};
    Configurable<std::string> cfgCcdbPathEff{"cfgCcdbPathEff", "", "path to the ccdb object for efficiency"};
    Configurable<float> cfgMinPt{"cfgMinPt", -1, "if ptgen is smaller than this threshold, this value is used as input for ptgen."};
  } kaon_filenames;

  struct : ConfigurableGroup {
    std::string prefix = "proton_filename_group";
    Configurable<bool> cfgNDSmearing{"cfgNDSmearing", false, "apply ND-correlated smearing"};
    Configurable<std::string> cfgResFileName{"cfgResFileName", "", "name of resolution file"};
    Configurable<std::string> cfgResNDHistName{"cfgResNDHistName", "hs_reso", "name of ND resolution file"};
    Configurable<std::string> cfgResPtHistName{"cfgResPtHistName", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
    Configurable<std::string> cfgResEtaHistName{"cfgResEtaHistName", "EtaResArr", "histogram name for eta in resolution file"};
    Configurable<std::string> cfgResPhiPosHistName{"cfgResPhiPosHistName", "PhiPosResArr", "histogram name for phi pos in resolution file"};
    Configurable<std::string> cfgResPhiNegHistName{"cfgResPhiNegHistName", "PhiEleResArr", "hisogram for phi neg in resolution file"};
    Configurable<std::string> cfgEffFileName{"cfgEffFileName", "", "name of efficiency file"};
    Configurable<std::string> cfgEffHistName{"cfgEffHistName", "fhwEffpT", "name of efficiency histogram"};
    Configurable<std::string> cfgCcdbPathRes{"cfgCcdbPathRes", "", "path to the ccdb object for resolution"};
    Configurable<std::string> cfgCcdbPathEff{"cfgCcdbPathEff", "", "path to the ccdb object for efficiency"};
    Configurable<float> cfgMinPt{"cfgMinPt", -1, "if ptgen is smaller than this threshold, this value is used as input for ptgen."};
  } proton_filenames;

  MomentumSmearer smearer_Electron;
  MomentumSmearer smearer_Pion;
  MomentumSmearer smearer_Kaon;
  MomentumSmearer smearer_Proton;
  Service<ccdb::BasicCCDBManager> ccdb;

  void init(InitContext&)
  {
    if (cfgCcdbTimestamp < 0) {
      LOG(fatal) << "Don't use time stamp = -1";
    }
    initResolutionMap(cfgCcdbTimestamp);
  }

  void initResolutionMap(const int64_t timestamp)
  {
    smearer_Electron.setNDSmearing(electron_filenames.cfgNDSmearing.value);
    smearer_Electron.setResFileName(TString(electron_filenames.cfgResFileName));
    smearer_Electron.setResNDHistName(TString(electron_filenames.cfgResNDHistName));
    smearer_Electron.setResPtHistName(TString(electron_filenames.cfgResPtHistName));
    smearer_Electron.setResEtaHistName(TString(electron_filenames.cfgResEtaHistName));
    smearer_Electron.setResPhiPosHistName(TString(electron_filenames.cfgResPhiPosHistName));
    smearer_Electron.setResPhiNegHistName(TString(electron_filenames.cfgResPhiNegHistName));
    smearer_Electron.setEffFileName(TString(electron_filenames.cfgEffFileName));
    smearer_Electron.setEffHistName(TString(electron_filenames.cfgEffHistName));
    smearer_Electron.setDCAFileName("");
    smearer_Electron.setDCAHistName("");
    smearer_Electron.setMinPt(electron_filenames.cfgMinPt);

    smearer_Pion.setNDSmearing(pion_filenames.cfgNDSmearing.value);
    smearer_Pion.setResFileName(TString(pion_filenames.cfgResFileName));
    smearer_Pion.setResNDHistName(TString(pion_filenames.cfgResNDHistName));
    smearer_Pion.setResPtHistName(TString(pion_filenames.cfgResPtHistName));
    smearer_Pion.setResEtaHistName(TString(pion_filenames.cfgResEtaHistName));
    smearer_Pion.setResPhiPosHistName(TString(pion_filenames.cfgResPhiPosHistName));
    smearer_Pion.setResPhiNegHistName(TString(pion_filenames.cfgResPhiNegHistName));
    smearer_Pion.setEffFileName(TString(pion_filenames.cfgEffFileName));
    smearer_Pion.setEffHistName(TString(pion_filenames.cfgEffHistName));
    smearer_Pion.setDCAFileName("");
    smearer_Pion.setDCAHistName("");
    smearer_Pion.setMinPt(pion_filenames.cfgMinPt);

    smearer_Kaon.setNDSmearing(kaon_filenames.cfgNDSmearing.value);
    smearer_Kaon.setResFileName(TString(kaon_filenames.cfgResFileName));
    smearer_Kaon.setResNDHistName(TString(kaon_filenames.cfgResNDHistName));
    smearer_Kaon.setResPtHistName(TString(kaon_filenames.cfgResPtHistName));
    smearer_Kaon.setResEtaHistName(TString(kaon_filenames.cfgResEtaHistName));
    smearer_Kaon.setResPhiPosHistName(TString(kaon_filenames.cfgResPhiPosHistName));
    smearer_Kaon.setResPhiNegHistName(TString(kaon_filenames.cfgResPhiNegHistName));
    smearer_Kaon.setEffFileName(TString(kaon_filenames.cfgEffFileName));
    smearer_Kaon.setEffHistName(TString(kaon_filenames.cfgEffHistName));
    smearer_Kaon.setDCAFileName("");
    smearer_Kaon.setDCAHistName("");
    smearer_Kaon.setMinPt(kaon_filenames.cfgMinPt);

    smearer_Proton.setNDSmearing(proton_filenames.cfgNDSmearing.value);
    smearer_Proton.setResFileName(TString(proton_filenames.cfgResFileName));
    smearer_Proton.setResNDHistName(TString(proton_filenames.cfgResNDHistName));
    smearer_Proton.setResPtHistName(TString(proton_filenames.cfgResPtHistName));
    smearer_Proton.setResEtaHistName(TString(proton_filenames.cfgResEtaHistName));
    smearer_Proton.setResPhiPosHistName(TString(proton_filenames.cfgResPhiPosHistName));
    smearer_Proton.setResPhiNegHistName(TString(proton_filenames.cfgResPhiNegHistName));
    smearer_Proton.setEffFileName(TString(proton_filenames.cfgEffFileName));
    smearer_Proton.setEffHistName(TString(proton_filenames.cfgEffHistName));
    smearer_Proton.setDCAFileName("");
    smearer_Proton.setDCAHistName("");
    smearer_Proton.setMinPt(proton_filenames.cfgMinPt);

    if (cfgFromCcdb) {
      ccdb->setURL(cfgCcdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()); // now

      smearer_Electron.setCcdbPathRes(TString(electron_filenames.cfgCcdbPathRes));
      smearer_Electron.setCcdbPathEff(TString(electron_filenames.cfgCcdbPathEff));
      // smearer_Electron.setCcdbPathDCA(TString(electron_filenames.fConfigCcdbPathDCA));
      smearer_Electron.setCcdbPathDCA("");
      smearer_Electron.setTimestamp(timestamp);
      smearer_Electron.setCcdb(ccdb);

      smearer_Pion.setCcdbPathRes(TString(pion_filenames.cfgCcdbPathRes));
      smearer_Pion.setCcdbPathEff(TString(pion_filenames.cfgCcdbPathEff));
      // smearer_Pion.setCcdbPathDCA(TString(pion_filenames.fConfigCcdbPathDCA));
      smearer_Pion.setCcdbPathDCA("");
      smearer_Pion.setTimestamp(timestamp);
      smearer_Pion.setCcdb(ccdb);

      smearer_Kaon.setCcdbPathRes(TString(kaon_filenames.cfgCcdbPathRes));
      smearer_Kaon.setCcdbPathEff(TString(kaon_filenames.cfgCcdbPathEff));
      // smearer_Kaon.setCcdbPathDCA(TString(kaon_filenames.fConfigCcdbPathDCA));
      smearer_Kaon.setCcdbPathDCA("");
      smearer_Kaon.setTimestamp(timestamp);
      smearer_Kaon.setCcdb(ccdb);

      smearer_Proton.setCcdbPathRes(TString(proton_filenames.cfgCcdbPathRes));
      smearer_Proton.setCcdbPathEff(TString(proton_filenames.cfgCcdbPathEff));
      // smearer_Proton.setCcdbPathDCA(TString(proton_filenames.fConfigCcdbPathDCA));
      smearer_Proton.setCcdbPathDCA("");
      smearer_Proton.setTimestamp(timestamp);
      smearer_Proton.setCcdb(ccdb);
    }

    smearer_Electron.init();
    smearer_Pion.init();
    smearer_Kaon.init();
    smearer_Proton.init();
  }

  void processACTSHybrid(MyTracks const& tracks, const aod::McParticles& /*mcParticles*/)
  {
    for (const auto& track : tracks) {
      float ptgen = track.pt();
      float etagen = track.eta();
      float phigen = track.phi();
      bool selected = true;
      float centrality = -1.f;
      if (track.has_mcParticle() && track.isReconstructed()) {
        const auto mcParticle = track.mcParticle_as<aod::McParticles>();
        ptgen = mcParticle.pt();
        etagen = mcParticle.eta();
        phigen = mcParticle.phi();
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kElectron) {
          int ch = -1;
          if (mcParticle.pdgCode() < 0) {
            ch = 1;
          }
          float ptsmeared = 0, etasmeared = 0, phismeared = 0;
          // apply smearing for electrons or muons.
          smearer_Electron.applySmearing(centrality, ch, ptgen, etagen, phigen, ptsmeared, etasmeared, phismeared);
          // get the efficiency
          float efficiency = smearer_Electron.getEfficiency(ptgen, etagen, phigen);
          // Generate a random double between 0 and 1
          double myRandom = gRandom->Uniform(0, 1);
          // Select
          if (myRandom < efficiency) {
            selected = true;
          } else {
            selected = false;
          }
          // fill the table
          smearedelectron(ptsmeared, etasmeared, phismeared, selected);
        } else if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus) {
          int ch = -1;
          if (mcParticle.pdgCode() < 0) {
            ch = 1;
          }
          float ptsmeared = 0, etasmeared = 0, phismeared = 0;
          // apply smearing for electrons or muons.
          smearer_Pion.applySmearing(centrality, ch, ptgen, etagen, phigen, ptsmeared, etasmeared, phismeared);
          // get the efficiency if there, otherwise applied zero
          if ((TString(pion_filenames.cfgCcdbPathEff).CompareTo("") != 0) || (TString(pion_filenames.cfgEffFileName).CompareTo("") != 0)) {
            float efficiency = smearer_Pion.getEfficiency(ptgen, etagen, phigen);
            // Generate a random double between 0 and 1
            double myRandom = gRandom->Uniform(0, 1);
            // Select
            if (myRandom < efficiency) {
              selected = true;
            } else {
              selected = false;
            }
          } else {
            selected = false;
          }
          // fill the table
          smearedelectron(ptsmeared, etasmeared, phismeared, selected);
        } else if (std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus) {
          int ch = -1;
          if (mcParticle.pdgCode() < 0) {
            ch = 1;
          }
          float ptsmeared = 0, etasmeared = 0, phismeared = 0;
          // apply smearing for electrons or muons.
          smearer_Kaon.applySmearing(centrality, ch, ptgen, etagen, phigen, ptsmeared, etasmeared, phismeared);
          // get the efficiency if there, otherwise applied zero
          if ((TString(kaon_filenames.cfgCcdbPathEff).CompareTo("") != 0) || (TString(kaon_filenames.cfgEffFileName).CompareTo("") != 0)) {
            float efficiency = smearer_Kaon.getEfficiency(ptgen, etagen, phigen);
            // Generate a random double between 0 and 1
            double myRandom = gRandom->Uniform(0, 1);
            // Select
            if (myRandom < efficiency) {
              selected = true;
            } else {
              selected = false;
            }
          } else {
            selected = false;
          }
          // fill the table
          smearedelectron(ptsmeared, etasmeared, phismeared, selected);
        } else if (std::abs(mcParticle.pdgCode()) == PDG_t::kProton) {
          int ch = -1;
          if (mcParticle.pdgCode() < 0) {
            ch = 1;
          }
          float ptsmeared = 0, etasmeared = 0, phismeared = 0;
          // apply smearing for electrons or muons.
          smearer_Proton.applySmearing(centrality, ch, ptgen, etagen, phigen, ptsmeared, etasmeared, phismeared);
          // get the efficiency if there, otherwise applied zero
          if ((TString(proton_filenames.cfgCcdbPathEff).CompareTo("") != 0) || (TString(proton_filenames.cfgEffFileName).CompareTo("") != 0)) {
            float efficiency = smearer_Proton.getEfficiency(ptgen, etagen, phigen);
            // Generate a random double between 0 and 1
            double myRandom = gRandom->Uniform(0, 1);
            // Select
            if (myRandom < efficiency) {
              selected = true;
            } else {
              selected = false;
            }
          } else {
            selected = false;
          }
          // fill the table
          smearedelectron(ptsmeared, etasmeared, phismeared, selected);
        } else {
          // don't apply smearing and reject completely
          smearedelectron(ptgen, etagen, phigen, false);
        }
      } else {
        // don't apply smearing
        smearedelectron(ptgen, etagen, phigen, false);
      }
    } // end of mc track loop
  }

  void processDummyACTSHybrid(MyTracks const& tracks)
  {
    for (const auto& track : tracks) {
      smearedelectron(track.pt(), track.eta(), track.phi(), true);
    }
  }

  PROCESS_SWITCH(alice3dileptonsmearer, processACTSHybrid, "Run for alice 3 ACTS hybrid", false);
  PROCESS_SWITCH(alice3dileptonsmearer, processDummyACTSHybrid, "Dummy process function", true);
};

struct alice3dileptonchecksmearer {
  // Resolution histos as cross check
  Configurable<bool> cfgUsePtVecRes{"cfgUsePtVecRes", true, "If true, non-linear pt bins predefined in res histos"};
  ConfigurableAxis ptResBinsVec{"ptResBinsVec", {0., 0., 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.20, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285, 0.29, 0.295, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80, 2.00, 2.40, 2.80, 3.20, 3.70, 4.50, 6.00, 8.00, 10., 12.0, 14., 16., 18., 20.}, "Pt binning vector for resolution"};
  ConfigurableAxis ptResBins{"ptResBins", {20, 0.f, 20.f}, "Pt binning for resolution"};
  ConfigurableAxis deltaptResBins{"deltaptResBins", {500, -1.f, 1.f}, "DeltaPt binning for resolution"};
  ConfigurableAxis deltaetaResBins{"deltaetaResBins", {500, -0.5f, 0.5f}, "DeltaEta binning for resolution"};
  ConfigurableAxis deltaphiResBins{"deltaphiResBins", {500, -0.5f, 0.5f}, "DeltaPhi binning for resolution"};
  ConfigurableAxis etaBins{"etaBins", {200, -1.f, 1.f}, "Eta binning for efficiency"};

  HistogramRegistry registry{"output", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    if (doprocessCheckACTSHybrid) {
      registry.add<TH2>("Electron/hCorrelation_Pt", "pT correlation;p_{T,l}^{gen} (GeV/c);p_{T,l}^{smeared} (GeV/c)", {HistType::kTH2F, {{1000, 0.0f, 10.0f}, {1000, 0.0f, 10.0f}}});
      registry.add<TH2>("Electron/hCorrelation_Eta", "eta correlation;#eta_{l}^{gen};#eta_{l}^{smeared}", {HistType::kTH2F, {{200, -1.0f, +1.0f}, {200, -1.0f, +1.0f}}});
      registry.add<TH2>("Electron/hCorrelation_Phi", "phi correlation;#varphi_{l}^{gen} (rad.);#varphi_{l}^{smeared} (rad.)", {HistType::kTH2F, {{100, 0.0f, TMath::TwoPi()}, {100, 0.0f, TMath::TwoPi()}}});

      // Binning for resolution
      AxisSpec axisPtRes{ptResBins, "#it{p}^{gen}_{T,l} (GeV/#it{c})"};
      AxisSpec axisDeltaptRes{deltaptResBins, "(p^{gen}_{T} - p^{rec}_{T}) / p^{gen}_{T}"};
      AxisSpec axisDeltaetaRes{deltaetaResBins, "#eta^{gen} - #eta^{rec}"};
      AxisSpec axisDeltaphiRes{deltaphiResBins, "#varphi^{gen} - #varphi^{rec} (rad.)"};
      // Binning for efficiency
      AxisSpec axiseta{etaBins, "#eta"};

      if (!cfgUsePtVecRes) {
        registry.add<TH2>("Electron/PtGen_DeltaPtOverPtGen", "", HistType::kTH2D, {axisPtRes, axisDeltaptRes}, true);
        registry.add<TH2>("Electron/PtGen_DeltaEta", "", HistType::kTH2D, {axisPtRes, axisDeltaetaRes}, true);
        registry.add<TH2>("Electron/PtGen_DeltaPhi_Neg", "", HistType::kTH2D, {axisPtRes, axisDeltaphiRes}, true);
        registry.add<TH2>("Electron/PtGen_DeltaPhi_Pos", "", HistType::kTH2D, {axisPtRes, axisDeltaphiRes}, true);
        registry.add<TH2>("Electron/PtEtaGen", "", HistType::kTH2D, {axisPtRes, axiseta}, true);
        registry.add<TH2>("Electron/PtEtaRec", "", HistType::kTH2D, {axisPtRes, axiseta}, true);
      } else {
        registry.add<TH2>("Electron/PtGen_DeltaPtOverPtGen", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,l} (GeV/#it{c})"}, axisDeltaptRes}, true);
        registry.add<TH2>("Electron/PtGen_DeltaEta", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,l} (GeV/#it{c})"}, axisDeltaetaRes}, true);
        registry.add<TH2>("Electron/PtGen_DeltaPhi_Neg", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,l} (GeV/#it{c})"}, axisDeltaphiRes}, true);
        registry.add<TH2>("Electron/PtGen_DeltaPhi_Pos", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,l} (GeV/#it{c})"}, axisDeltaphiRes}, true);
        registry.add<TH2>("Electron/PtEtaGen", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,l} (GeV/#it{c})"}, axiseta}, true);
        registry.add<TH2>("Electron/PtEtaRec", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,l} (GeV/#it{c})"}, axiseta}, true);
      }

      registry.addClone("Electron/", "Other/");
      registry.addClone("Electron/", "Pion/");
      registry.addClone("Electron/", "Kaon/");
      registry.addClone("Electron/", "Proton/");
    }
  }

  void processCheckACTSHybrid(MyTracksWithSmearing const& tracks, const aod::McParticles& /*mcParticles*/)
  {
    for (const auto& track : tracks) {
      if (!track.has_mcParticle() || !track.isReconstructed()) {
        continue;
      }
      const auto mcParticle = track.mcParticle_as<aod::McParticles>();
      float deltaptoverpt = -1000.f;
      if (mcParticle.pt() > 0.f) {
        deltaptoverpt = (mcParticle.pt() - track.ptSmeared()) / mcParticle.pt();
      }
      float deltaeta = mcParticle.eta() - track.etaSmeared();
      float deltaphi = mcParticle.phi() - track.phiSmeared();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kElectron) {
        registry.fill(HIST("Electron/PtGen_DeltaPtOverPtGen"), mcParticle.pt(), deltaptoverpt);
        registry.fill(HIST("Electron/PtGen_DeltaEta"), mcParticle.pt(), deltaeta);
        if (mcParticle.pdgCode() < 0) { // e+
          registry.fill(HIST("Electron/PtGen_DeltaPhi_Pos"), mcParticle.pt(), deltaphi);
        } else { // e-
          registry.fill(HIST("Electron/PtGen_DeltaPhi_Neg"), mcParticle.pt(), deltaphi);
        }
        registry.fill(HIST("Electron/hCorrelation_Pt"), mcParticle.pt(), track.ptSmeared());
        registry.fill(HIST("Electron/hCorrelation_Eta"), mcParticle.eta(), track.etaSmeared());
        registry.fill(HIST("Electron/hCorrelation_Phi"), mcParticle.phi(), track.phiSmeared());
        // efficiency
        registry.fill(HIST("Electron/PtEtaGen"), mcParticle.pt(), mcParticle.eta());
        if (track.selected()) {
          registry.fill(HIST("Electron/PtEtaRec"), mcParticle.pt(), mcParticle.eta());
        }
      } else if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus) {
        registry.fill(HIST("Pion/PtGen_DeltaPtOverPtGen"), mcParticle.pt(), deltaptoverpt);
        registry.fill(HIST("Pion/PtGen_DeltaEta"), mcParticle.pt(), deltaeta);
        if (mcParticle.pdgCode() < 0) { // e+
          registry.fill(HIST("Pion/PtGen_DeltaPhi_Pos"), mcParticle.pt(), deltaphi);
        } else { // e-
          registry.fill(HIST("Pion/PtGen_DeltaPhi_Neg"), mcParticle.pt(), deltaphi);
        }
        registry.fill(HIST("Pion/hCorrelation_Pt"), mcParticle.pt(), track.ptSmeared());
        registry.fill(HIST("Pion/hCorrelation_Eta"), mcParticle.eta(), track.etaSmeared());
        registry.fill(HIST("Pion/hCorrelation_Phi"), mcParticle.phi(), track.phiSmeared());
        // efficiency
        registry.fill(HIST("Pion/PtEtaGen"), mcParticle.pt(), mcParticle.eta());
        if (track.selected()) {
          registry.fill(HIST("Pion/PtEtaRec"), mcParticle.pt(), mcParticle.eta());
        }
      } else if (std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus) {
        registry.fill(HIST("Kaon/PtGen_DeltaPtOverPtGen"), mcParticle.pt(), deltaptoverpt);
        registry.fill(HIST("Kaon/PtGen_DeltaEta"), mcParticle.pt(), deltaeta);
        if (mcParticle.pdgCode() < 0) { // e+
          registry.fill(HIST("Kaon/PtGen_DeltaPhi_Pos"), mcParticle.pt(), deltaphi);
        } else { // e-
          registry.fill(HIST("Kaon/PtGen_DeltaPhi_Neg"), mcParticle.pt(), deltaphi);
        }
        registry.fill(HIST("Kaon/hCorrelation_Pt"), mcParticle.pt(), track.ptSmeared());
        registry.fill(HIST("Kaon/hCorrelation_Eta"), mcParticle.eta(), track.etaSmeared());
        registry.fill(HIST("Kaon/hCorrelation_Phi"), mcParticle.phi(), track.phiSmeared());
        // efficiency
        registry.fill(HIST("Kaon/PtEtaGen"), mcParticle.pt(), mcParticle.eta());
        if (track.selected()) {
          registry.fill(HIST("Kaon/PtEtaRec"), mcParticle.pt(), mcParticle.eta());
        }
      } else if (std::abs(mcParticle.pdgCode()) == PDG_t::kProton) {
        registry.fill(HIST("Proton/PtGen_DeltaPtOverPtGen"), mcParticle.pt(), deltaptoverpt);
        registry.fill(HIST("Proton/PtGen_DeltaEta"), mcParticle.pt(), deltaeta);
        if (mcParticle.pdgCode() < 0) { // e+
          registry.fill(HIST("Proton/PtGen_DeltaPhi_Pos"), mcParticle.pt(), deltaphi);
        } else { // e-
          registry.fill(HIST("Proton/PtGen_DeltaPhi_Neg"), mcParticle.pt(), deltaphi);
        }
        registry.fill(HIST("Proton/hCorrelation_Pt"), mcParticle.pt(), track.ptSmeared());
        registry.fill(HIST("Proton/hCorrelation_Eta"), mcParticle.eta(), track.etaSmeared());
        registry.fill(HIST("Proton/hCorrelation_Phi"), mcParticle.phi(), track.phiSmeared());
        // efficiency
        registry.fill(HIST("Proton/PtEtaGen"), mcParticle.pt(), mcParticle.eta());
        if (track.selected()) {
          registry.fill(HIST("Proton/PtEtaRec"), mcParticle.pt(), mcParticle.eta());
        }
      } else {
        registry.fill(HIST("Other/PtGen_DeltaPtOverPtGen"), mcParticle.pt(), deltaptoverpt);
        registry.fill(HIST("Other/PtGen_DeltaEta"), mcParticle.pt(), deltaeta);
        if (mcParticle.pdgCode() < 0) { // e+
          registry.fill(HIST("Other/PtGen_DeltaPhi_Pos"), mcParticle.pt(), deltaphi);
        } else { // e-
          registry.fill(HIST("Other/PtGen_DeltaPhi_Neg"), mcParticle.pt(), deltaphi);
        }
        registry.fill(HIST("Other/hCorrelation_Pt"), mcParticle.pt(), track.ptSmeared());
        registry.fill(HIST("Other/hCorrelation_Eta"), mcParticle.eta(), track.etaSmeared());
        registry.fill(HIST("Other/hCorrelation_Phi"), mcParticle.phi(), track.phiSmeared());
        // efficiency
        registry.fill(HIST("Other/PtEtaGen"), mcParticle.pt(), mcParticle.eta());
        if (track.selected()) {
          registry.fill(HIST("Other/PtEtaRec"), mcParticle.pt(), mcParticle.eta());
        }
      }
    } // end of loop
  }

  void processDummyACTSHybrid(aod::McParticles const&) {}

  PROCESS_SWITCH(alice3dileptonchecksmearer, processCheckACTSHybrid, "Run for alice3 analysis", false);
  PROCESS_SWITCH(alice3dileptonchecksmearer, processDummyACTSHybrid, "Dummy process function", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<alice3dileptonsmearer>(cfgc),
    adaptAnalysisTask<alice3dileptonchecksmearer>(cfgc)};
}
