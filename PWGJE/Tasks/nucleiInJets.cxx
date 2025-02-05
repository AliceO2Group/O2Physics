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

// author: Arvind Khuntia (arvind.khuntia@cern.ch) INFN Bologna, Italy

#include <string>
#include <vector>
#include <TLorentzVector.h>
#include <TVector2.h>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"
#include "ReconstructionDataFormats/Track.h"

#include "PWGLF/DataModel/LFParticleIdentification.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
struct nucleiInJets {

  enum Particle {
    kPion = 1,     // π+
    kKaon = 2,     // K+
    kProton = 3,   // p
    kDeuteron = 4, // d
    kTriton = 5,   // Tr
    kHelium = 6    // He
  };

  int mapPDGToValue(int pdgCode)
  {
    switch (pdgCode) {
      case 211: // π+
        return Particle::kPion;
      case -211: // π-
        return -Particle::kPion;

      case 321: // k+
        return Particle::kKaon;
      case -321: // k-
        return -Particle::kKaon;

      case 2212: // p
        return Particle::kProton;
      case -2212: // antip
        return -Particle::kProton;

      case 1000010020: // Deuteron
        return Particle::kDeuteron;
      case -1000010020: // AntiDeuteron
        return -Particle::kDeuteron;

      case 1000010030: // Triton
        return Particle::kTriton;
      case -1000010030: // AntiTriton
        return -Particle::kTriton;

      case 1000020030: // Helium
        return Particle::kHelium;
      case -1000020030: // AntiHelium
        return -Particle::kHelium;
      default:
        return 0; // Default case for unknown or unmapped PDG codes
    }
  }

  Configurable<std::string> cfgtrackSelections{"cfgtrackSelections", "globalTracks", "set track selections"};

  Configurable<bool> isMC{"isMC", false, "flag for the MC"};
  Configurable<bool> isWithJetEvents{"isWithJetEvents", true, "Events with at least one jet"};
  Configurable<bool> isWithLeadingJet{"isWithLeadingJet", true, "Events with leading jet"};
  Configurable<bool> useLfTpcPid{"useLfTpcPid", true, "Events with custom TPC parameters"};

  Configurable<double> cfgtrkMinPt{"cfgtrkMinPt", 0.15, "set track min pT"};
  Configurable<double> cfgtrkMaxEta{"cfgtrkMaxEta", 0.8, "set track max Eta"};
  Configurable<double> cfgtrkMaxRap{"cfgtrkMaxRap", 0.5, "set track max y"};
  Configurable<double> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.12, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 1.0, "Track DCAz cut to PV Maximum"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};
  Configurable<bool> cfgConnectedToPV{"cfgConnectedToPV", true, "PV contributor track selection"};           // PV Contriuibutor
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<double> cfgnFindableTPCClusters{"cfgnFindableTPCClusters", 120, "nFindable TPC Clusters"};
  Configurable<double> cfgnTPCCrossedRows{"cfgnTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgnTPCChi2{"cfgnTPChi2", 4.0, "nTPC Chi2 per Cluster"};
  Configurable<double> cfgnITSChi2{"cfgnITShi2", 36.0, "nITS Chi2 per Cluster"};

  Configurable<float> cfgjetPtMin{"cfgjetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> cfgjetR{"cfgjetR", 0.4, "jet resolution parameter"};
  Configurable<int> cDebugLevel{"cDebugLevel", 0, "print debug msg"};
  Configurable<int> cMaxPt{"cMaxPt", 10, "max pt for Hist"};

  Configurable<int> cfgnTPCPIDPr{"cfgnTPCPIDPr", 3, "nTPC PID Pr"};
  Configurable<int> cfgnTPCPIDDe{"cfgnTPCPIDDe", 3, "nTPC PID De"};
  Configurable<int> cfgnTPCPIDHe{"cfgnTPCPIDHe", 3, "nTPC PID He"};
  Configurable<int> cfgnTPCPIDTr{"cfgnTPCPIDTr", 3, "nTPC PID Tr"};

  Configurable<bool> cEnableProtonQA{"cEnableProtonQA", true, "nTPC PID Pr"};
  Configurable<bool> cEnableDeuteronQA{"cEnableDeuteronQA", true, "nTPC PID De"};
  Configurable<bool> cEnableHeliumQA{"cEnableHeliumQA", true, "nTPC PID He"};
  Configurable<bool> cEnableTritonQA{"cEnableTritonQA", true, "nTPC PID Tr"};
  Configurable<bool> addTOFplots{"addTOFplots", true, "add TOF plots"};
  Configurable<int> useTPCpreSel{"useTPCpreSel", 3, "add TPC nsgma preselection for TOF: (0) no selection (!0) selction on TPC"};

  Configurable<bool> addpik{"addpik", true, "add pion and kaon hist"};

  ConfigurableAxis binsDCA{"binsDCA", {400, -1.f, 1.f}, ""};
  ConfigurableAxis binsdEdx{"binsdEdx", {1000, 0.f, 1000.f}, ""};
  ConfigurableAxis binsBeta{"binsBeta", {120, 0.0, 1.2}, ""};

  ConfigurableAxis binsMassPr{"binsMassPr", {100, -1., 1.f}, ""};
  ConfigurableAxis binsMassDe{"binsMassDe", {180, -1.8, 1.8f}, ""};
  ConfigurableAxis binsMassTr{"binsMassTr", {250, -2.5, 2.5f}, ""};
  ConfigurableAxis binsMassHe{"binsMassHe", {300, -3., 3.f}, ""};

  ConfigurableAxis binsPtZHe{"binsPtZHe", {VARIABLE_WIDTH, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0}, ""};

  static constexpr float gMassProton = 0.93827208f;
  static constexpr float gMassDeuteron = 1.87561f;
  static constexpr float gMassTriton = 2.80892f;
  static constexpr float gMassHelium = 2.80839f;
  static constexpr int PDGProton = 2212;
  static constexpr int PDGDeuteron = 1000010020;
  static constexpr int PDGTriton = 1000010030;
  static constexpr int PDGHelium = 1000020030;

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs>; // , aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa,
                                    aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTPCFullTr, aod::pidTOFFullPi, aod::pidTOFFullKa,
                                    aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFbeta, aod::TOFSignal>;
  using TrackCandidatesLfPid = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa,
                                         aod::pidTPCLfFullPr, aod::pidTPCLfFullDe, aod::pidTPCLfFullHe, aod::pidTPCLfFullTr, aod::pidTOFFullPi, aod::pidTOFFullKa,
                                         aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFbeta, aod::TOFSignal>;
  using TrackCandidatesMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa,
                                      aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTPCFullTr, aod::pidTOFFullPi,
                                      aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe,
                                      aod::TOFSignal /*, aod::McTrackLabels*/>;

  Filter jetCuts = aod::jet::pt > cfgjetPtMin&& aod::jet::r == nround(cfgjetR.node() * 100.0f);

  using chargedJetstrack = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>;
  using JetMCPartTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  using JetMCDetTable = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;

  SliceCache cache;
  HistogramRegistry jetHist{"jetHist", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec PtAxis = {100, 0, 10.0};
    const AxisSpec PtJetAxis = {300, 0, 30.0};
    const AxisSpec MultAxis = {100, 0, 100};
    const AxisSpec dRAxis = {100, 0, 3.6};
    const AxisSpec dcaxyAxis{binsDCA, "DCAxy (cm)"};
    const AxisSpec dcazAxis{binsDCA, "DCAz (cm)"};
    const AxisSpec dedxAxis{binsdEdx, "d#it{E}/d#it{x} A.U."};

    const AxisSpec betaAxis{binsBeta, "TOF #beta"};
    const AxisSpec ptZHeAxis{binsPtZHe, "#it{p}_{T}"};

    const AxisSpec massPrAxis{binsMassPr, ""};
    const AxisSpec massDeAxis{binsMassDe, ""};
    const AxisSpec massTrAxis{binsMassTr, ""};
    const AxisSpec massHeAxis{binsMassHe, ""};

    // jet property
    jetHist.add("jet/h1JetPt", "jet_{p_{T}}", kTH1F, {PtJetAxis});
    jetHist.add("jet/h1JetEvents", "NumbeOfJetEvents", kTH1F, {{1, 0, 1}});
    jetHist.add("jet/h1JetEta", "jet_{#eta}", kTH1F, {{100, -1.0, 1.0}});
    jetHist.add("jet/h1JetPhi", "jet_{#phi}", kTH1F, {{80, -1.0, 7.}});
    jetHist.add("jet/nJetsPerEvent", "nJetsPerEvent", kTH1F, {{15, .0, 15.}});
    jetHist.add("mcpJet/nJetsPerEvent", "nJetsPerEvent", kTH1F, {{15, .0, 15.}});
    jetHist.add("mcdJet/nJetsPerEvent", "nJetsPerEvent", kTH1F, {{15, .0, 15.}});
    jetHist.add("jet/vertexZ", "vertexZ (Jet flag)", kTH1F, {{100, -15.0, 15.0}});
    jetHist.add("vertexZ", "vertexZ (all)", kTH1F, {{100, -15.0, 15.0}});
    jetHist.add("jetOut/vertexZ", "vertexZ (without z-flag)", kTH1F, {{100, -15.0, 15.0}});
    ////////////////////////////
    //           MC
    ////////////////////////////
    jetHist.add("mcpJet/eventStat", "vertexZ (All)", kTH1F, {{5, .0, 5.0}});
    auto h = jetHist.get<TH1>(HIST("mcpJet/eventStat"));
    h->GetXaxis()->SetBinLabel(1, "All");
    h->GetXaxis()->SetBinLabel(2, "Sel8-goodRecJet");
    h->GetXaxis()->SetBinLabel(3, "vz < 10");
    h->GetXaxis()->SetBinLabel(4, "ingt0");

    jetHist.add("mcpJet/vertexZ", "vertexZ (All)", kTH1F, {{100, -15.0, 15.0}});
    jetHist.add("mcdJet/vertexZ", "vertexZ (All)", kTH1F, {{100, -15.0, 15.0}});
    jetHist.add("mcdJet/eventStat", "vertexZ (All)", kTH1F, {{5, .0, 5.0}});
    auto h1 = jetHist.get<TH1>(HIST("mcdJet/eventStat"));
    h1->GetXaxis()->SetBinLabel(1, "All");
    h1->GetXaxis()->SetBinLabel(2, "Sel8-goodRecJet");
    h1->GetXaxis()->SetBinLabel(3, "vz< 10");
    h1->GetXaxis()->SetBinLabel(4, "ingt0");

    jetHist.add("recmatched/vertexZ", "vertexZ (All)", kTH1F, {{100, -15.0, 15.0}});
    jetHist.add("genmatched/vertexZ", "vertexZ (All)", kTH1F, {{100, -15.0, 15.0}});

    //////////////////////////////////////////////
    //              inside jet
    //////////////////////////////////////////////
    if (addpik) {
      jetHist.add<TH3>("tracks/pion/h3PtVsPionNSigmaTPCVsPtJet_jet", "pT(p) vs NSigmaTPC (p) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
      jetHist.add<TH3>("tracks/antiPion/h3PtVsPionNSigmaTPCVsPtJet_jet", "pT(p) vs NSigmaTPC (p) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
      jetHist.add<TH3>("tracks/kaon/h3PtVsKaonNSigmaTPCVsPtJet_jet", "pT(p) vs NSigmaTPC (p) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
      jetHist.add<TH3>("tracks/antiKaon/h3PtVsKaonNSigmaTPCVsPtJet_jet", "pT(p) vs NSigmaTPC (p) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    }
    jetHist.add<TH3>("tracks/proton/h3PtVsProtonNSigmaTPCVsPtJet_jet", "pT(p) vs NSigmaTPC (p) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/antiProton/h3PtVsantiProtonNSigmaTPCVsPtJet_jet", "pT(#bar{p}) vs NSigmaTPC (#bar{p}) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/deuteron/h3PtVsDeuteronNSigmaTPCVsPtJet_jet", "pT(d) vs NSigmaTPC (d) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/antiDeuteron/h3PtVsantiDeuteronNSigmaTPCVsPtJet_jet", "pT(#bar{d}) vs NSigmaTPC (#bar{d}) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/helium/h3PtVsHeliumNSigmaTPCVsPtJet_jet", "pT(He) vs NSigmaTPC (He) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/antiHelium/h3PtVsantiHeliumNSigmaTPCVsPtJet_jet", "pT(#bar{He}) vs NSigmaTPC (#bar{He}) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/triton/h3PtVsTritonNSigmaTPCVsPtJet_jet", "pT(Tr) vs NSigmaTPC (Tr) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/antiTriton/h3PtVsantiTritonNSigmaTPCVsPtJet_jet", "pT(#bar{Tr}) vs NSigmaTPC (#bar{Tr}) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});

    if (cEnableProtonQA) {
      jetHist.add<TH2>("tracks/proton/dca/after/hDCAxyVsPtProton_jet", "DCAxy vs Pt (p)", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/antiProton/dca/after/hDCAxyVsPtantiProton_jet", "DCAxy vs Pt (#bar{p})", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/proton/dca/after/hDCAzVsPtProton_jet", "DCAz vs Pt (p)", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/antiProton/dca/after/hDCAzVsPtantiProton_jet", "DCAz vs Pt (#bar{p})", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
    }

    if (cEnableDeuteronQA) {
      jetHist.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtDeuteron_jet", "DCAxy vs Pt (d)", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/antiDeuteron/dca/after/hDCAxyVsPtantiDeuteron_jet", "DCAxy vs Pt (#bar{d})", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtDeuteron_jet", "DCAz vs Pt (d)", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/antiDeuteron/dca/after/hDCAzVsPtantiDeuteron_jet", "DCAz vs Pt (#bar{d})", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
    }

    if (cEnableTritonQA) {
      jetHist.add<TH2>("tracks/triton/dca/after/hDCAxyVsPtTriton_jet", "DCAxy vs Pt (t)", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/antiTriton/dca/after/hDCAxyVsPtantiTriton_jet", "DCAxy vs Pt (#bar{t})", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/triton/dca/after/hDCAzVsPtTriton_jet", "DCAz vs Pt (t)", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/antiTriton/dca/after/hDCAzVsPtantiTriton_jet", "DCAz vs Pt (#bar{t})", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
    }
    if (cEnableHeliumQA) {
      jetHist.add<TH2>("tracks/helium/dca/after/hDCAxyVsPtHelium_jet", "DCAxy vs Pt (He)", HistType::kTH2F, {{dcaxyAxis}, {450, 0.5f, 5.f}});
      jetHist.add<TH2>("tracks/antiHelium/dca/after/hDCAxyVsPtantiHelium_jet", "DCAxy vs Pt (#bar{He})", HistType::kTH2F, {{dcaxyAxis}, {450, 0.5f, 5.f}});
      jetHist.add<TH2>("tracks/helium/dca/after/hDCAzVsPtHelium_jet", "DCAz vs Pt (He)", HistType::kTH2F, {{dcazAxis}, {450, 0.5f, 5.f}});
      jetHist.add<TH2>("tracks/antiHelium/dca/after/hDCAzVsPtantiHelium_jet", "DCAz vs Pt (#bar{He})", HistType::kTH2F, {{dcazAxis}, {450, 0.5f, 5.f}});
    }

    jetHist.add<TH2>("tracks/h2TPCsignVsTPCmomentum", "TPC <-dE/dX> vs #it{p}/Z;#it{p}/Z (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{500, -5.f, 5.f}, {dedxAxis}});
    jetHist.add<TH2>("tracks/h2TPCsignVsTPCmomentum_Jet", "TPC <-dE/dX> vs #it{p}/Z;#it{p}/Z (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{500, -5.f, 5.f}, {dedxAxis}});
    jetHist.add<TH2>("tracks/h2TPCsignVsTPCmomentum_OutJet", "TPC <-dE/dX> vs #it{p}/Z;#it{p}/Z (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{500, -5.f, 5.f}, {dedxAxis}});
    jetHist.add<TH2>("tracks/perpCone/h2TPCsignVsTPCmomentum", "TPC <-dE/dX> vs #it{p}/Z;#it{p}/Z (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{500, -5.f, 5.f}, {dedxAxis}});

    jetHist.add<TH2>("tracks/h2TOFbetaVsP_Jet", "TOF #beta vs #it{p}/Z; #it{p}/Z (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{250, -5.f, 5.f}, {betaAxis}});
    jetHist.add<TH2>("tracks/h2TOFbetaVsP_OutJet", "TOF #beta vs #it{p}/Z; #it{p}/Z (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{250, -5.f, 5.f}, {betaAxis}});
    jetHist.add<TH2>("tracks/h2TOFbetaVsP", "TOF #beta vs #it{p}/Z; #it{p}/Z (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{250, -5.f, 5.f}, {betaAxis}});

    // TOF hist
    jetHist.add<TH2>("tracks/proton/h2TOFmassProtonVsPt_jet", "h2TOFmassProtonVsPt_jet; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 4.}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiProton/h2TOFmassantiProtonVsPt_jet", "h2TOFmassantiProtonVsPt_jet; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 4.}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/proton/h2TOFmass2ProtonVsPt_jet", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T} (GeV)", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
    jetHist.add<TH2>("tracks/antiProton/h2TOFmass2antiProtonVsPt_jet", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T} (GeV)", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});

    jetHist.add<TH2>("tracks/deuteron/h2TOFmassDeuteronVsPt_jet", "h2TOFmassDeuteronVsPt_jet; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 4.}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiDeuteron/h2TOFmassantiDeuteronVsPt_jet", "h2TOFmassantiDeuteronVsPt_jet; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 4.}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/deuteron/h2TOFmass2DeuteronVsPt_jet", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T} (GeV)", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
    jetHist.add<TH2>("tracks/antiDeuteron/h2TOFmass2antiDeuteronVsPt_jet", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T} (GeV)", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});

    jetHist.add<TH2>("tracks/triton/h2TOFmassTritonVsPt_jet", "h2TOFmassTritonVsPt_jet; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 4.}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiTriton/h2TOFmassantiTritonVsPt_jet", "h2TOFmassantiTritonVsPt_jet; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 4.}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/triton/h2TOFmass2TritonVsPt_jet", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T} (GeV)", HistType::kTH2F, {{massTrAxis}, {250, 0., 5.}});
    jetHist.add<TH2>("tracks/antiTriton/h2TOFmass2antiTritonVsPt_jet", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T} (GeV)", HistType::kTH2F, {{massTrAxis}, {250, 0., 5.}});

    jetHist.add<TH2>("tracks/helium/h2TOFmassHeliumVsPt_jet", "h2TOFmassHeliumVsPt_jet; TOFmass; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {ptZHeAxis}});
    jetHist.add<TH2>("tracks/antiHelium/h2TOFmassantiHeliumVsPt_jet", "h2TOFmassantiHeliumVsPt_jet; TOFmass; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {ptZHeAxis}});
    jetHist.add<TH2>("tracks/helium/h2TOFmass2HeliumVsPt_jet", "#Delta M^{2} (t) vs #it{p}_{T}t; TOFmass2; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{massHeAxis}, {ptZHeAxis}});
    jetHist.add<TH2>("tracks/antiHelium/h2TOFmass2antiHeliumVsPt_jet", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{massHeAxis}, {ptZHeAxis}});

    // TOF hist nSigma
    if (addpik) {
      jetHist.add<TH2>("tracks/pion/h2TofNsigmaPionVsPt_jet", "h2TofNsigmaPionVsPt_jet; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
      jetHist.add<TH2>("tracks/antiPion/h2TofNsigmaantiPionVsPt_jet", "h2TofNsigmaantiPionVsPt_jet; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
      jetHist.add<TH2>("tracks/kaon/h2TofNsigmaKaonVsPt_jet", "h2TofNsigmaKaonVsPt_jet; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
      jetHist.add<TH2>("tracks/antiKaon/h2TofNsigmaantiKaonVsPt_jet", "h2TofNsigmaantiKaonVsPt_jet; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    }
    jetHist.add<TH2>("tracks/proton/h2TofNsigmaProtonVsPt_jet", "h2TofNsigmaProtonVsPt_jet; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiProton/h2TofNsigmaantiProtonVsPt_jet", "h2TofNsigmaantiProtonVsPt_jet; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/deuteron/h2TofNsigmaDeuteronVsPt_jet", "h2TofNsigmaDeuteronVsPt_jet; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiDeuteron/h2TofNsigmaantiDeuteronVsPt_jet", "h2TofNsigmaantiDeuteronVsPt_jet; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/triton/h2TofNsigmaTritonVsPt_jet", "h2TofNsigmaTritonVsPt_jet; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiTriton/h2TofNsigmaantiTritonVsPt_jet", "h2TofNsigmaantiTritonVsPt_jet; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/helium/h2TofNsigmaHeliumVsPt_jet", "h2TofNsigmaHeliumVsPt_jet; TofNsigma; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiHelium/h2TofNsigmaantiHeliumVsPt_jet", "h2TofNsigmaantiHeliumVsPt_jet; TofNsigma; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    /////////////
    // perp cone
    /////////////
    jetHist.add<TH2>("tracks/perpCone/proton/h2TofNsigmaProtonVsPt", "h2TofNsigmaProtonVsPt; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/perpCone/antiProton/h2TofNsigmaantiProtonVsPt", "h2TofNsigmaantiProtonVsPt; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/perpCone/deuteron/h2TofNsigmaDeuteronVsPt", "h2TofNsigmaDeuteronVsPt; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/perpCone/antiDeuteron/h2TofNsigmaantiDeuteronVsPt", "h2TofNsigmaantiDeuteronVsPt; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/perpCone/triton/h2TofNsigmaTritonVsPt", "h2TofNsigmaTritonVsPt; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/perpCone/antiTriton/h2TofNsigmaantiTritonVsPt", "h2TofNsigmaantiTritonVsPt; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/perpCone/helium/h2TofNsigmaHeliumVsPt", "h2TofNsigmaHeliumVsPt; TofNsigma; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/perpCone/antiHelium/h2TofNsigmaantiHeliumVsPt", "h2TofNsigmaantiHeliumVsPt; TofNsigma; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    //////////////////////////////////////////////
    //               outside jet
    //////////////////////////////////////////////
    jetHist.add<TH2>("tracks/proton/h3PtVsProtonNSigmaTPC", "pT(p) vs NSigmaTPC (p);  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;", HistType::kTH2F, {{PtAxis}, {200, -10, 10}});
    jetHist.add<TH2>("tracks/antiProton/h3PtVsantiProtonNSigmaTPC", "pT(#bar{p}) vs NSigmaTPC (#bar{p});  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;", HistType::kTH2F, {{PtAxis}, {200, -10, 10}});
    jetHist.add<TH2>("tracks/deuteron/h3PtVsDeuteronNSigmaTPC", "pT(d) vs NSigmaTPC (d);  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;", HistType::kTH2F, {{PtAxis}, {200, -10, 10}});
    jetHist.add<TH2>("tracks/antiDeuteron/h3PtVsantiDeuteronNSigmaTPC", "pT(#bar{d}) vs NSigmaTPC (#bar{d});  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;", HistType::kTH2F, {{PtAxis}, {200, -10, 10}});
    jetHist.add<TH2>("tracks/helium/h3PtVsHeliumNSigmaTPC", "pT(He) vs NSigmaTPC (He);  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;", HistType::kTH2F, {{PtAxis}, {200, -10, 10}});
    jetHist.add<TH2>("tracks/antiHelium/h3PtVsantiHeliumNSigmaTPC", "pT(#bar{He}) vs NSigmaTPC (#bar{He});  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;", HistType::kTH2F, {{PtAxis}, {200, -10, 10}});
    jetHist.add<TH2>("tracks/triton/h3PtVsTritonNSigmaTPC", "pT(Tr) vs NSigmaTPC(Tr);  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;", HistType::kTH2F, {{PtAxis}, {200, -10, 10}});
    jetHist.add<TH2>("tracks/antiTriton/h3PtVsantiTritonNSigmaTPC", "pT(#barTr}) vs NSigmaTPC (#bar{Tr});  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;", HistType::kTH2F, {{PtAxis}, {200, -10, 10}});

    jetHist.add<TH3>("tracks/perpCone/proton/h3PtVsProtonNSigmaTPCVsPtJet", "pT(p) vs NSigmaTPC (p) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/perpCone/antiProton/h3PtVsantiProtonNSigmaTPCVsPtJet", "pT(#bar{p}) vs NSigmaTPC (#bar{p}) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/perpCone/deuteron/h3PtVsDeuteronNSigmaTPCVsPtJet", "pT(d) vs NSigmaTPC (d) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/perpCone/antiDeuteron/h3PtVsantiDeuteronNSigmaTPCVsPtJet", "pT(#bar{d}) vs NSigmaTPC (#bar{d}) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/perpCone/helium/h3PtVsHeliumNSigmaTPCVsPtJet", "pT(He) vs NSigmaTPC (He) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/perpCone/antiHelium/h3PtVsantiHeliumNSigmaTPCVsPtJet", "pT(#bar{He}) vs NSigmaTPC (#bar{He}) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/perpCone/triton/h3PtVsTritonNSigmaTPCVsPtJet", "pT(Tr) vs NSigmaTPC (Tr) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
    jetHist.add<TH3>("tracks/perpCone/antiTriton/h3PtVsantiTritonNSigmaTPCVsPtJet", "pT(#bar{Tr}) vs NSigmaTPC (#bar{Tr}) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});

    if (cEnableProtonQA) {
      jetHist.add<TH2>("tracks/proton/dca/after/hDCAxyVsPtProton", "DCAxy vs Pt (p)", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/antiProton/dca/after/hDCAxyVsPtantiProton", "DCAxy vs Pt (#bar{p})", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/proton/dca/after/hDCAzVsPtProton", "DCAz vs Pt (p)", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/antiProton/dca/after/hDCAzVsPtantiProton", "DCAz vs Pt (#bar{p})", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
    }
    if (cEnableDeuteronQA) {
      jetHist.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtDeuteron", "DCAxy vs Pt (d)", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/antiDeuteron/dca/after/hDCAxyVsPtantiDeuteron", "DCAxy vs Pt (#bar{d})", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtDeuteron", "DCAz vs Pt (d)", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/antiDeuteron/dca/after/hDCAzVsPtantiDeuteron", "DCAz vs Pt (#bar{d})", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
    }
    if (cEnableTritonQA) {
      jetHist.add<TH2>("tracks/triton/dca/after/hDCAxyVsPtTriton", "DCAxy vs Pt (t)", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/antiTriton/dca/after/hDCAxyVsPtantiTriton", "DCAxy vs Pt (#bar{t})", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/triton/dca/after/hDCAzVsPtTriton", "DCAz vs Pt (t)", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
      jetHist.add<TH2>("tracks/antiTriton/dca/after/hDCAzVsPtantiTriton", "DCAz vs Pt (#bar{t})", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
    }
    if (cEnableHeliumQA) {
      jetHist.add<TH2>("tracks/helium/dca/after/hDCAxyVsPtHelium", "DCAxy vs Pt (He)", HistType::kTH2F, {{dcaxyAxis}, {450, 0.5f, 5.f}});
      jetHist.add<TH2>("tracks/antiHelium/dca/after/hDCAxyVsPtantiHelium", "DCAxy vs Pt (#bar{He})", HistType::kTH2F, {{dcaxyAxis}, {450, 0.5f, 5.f}});
      jetHist.add<TH2>("tracks/helium/dca/after/hDCAzVsPtHelium", "DCAz vs Pt (He)", HistType::kTH2F, {{dcazAxis}, {450, 0.5f, 5.f}});
      jetHist.add<TH2>("tracks/antiHelium/dca/after/hDCAzVsPtantiHelium", "DCAz vs Pt (#bar{He})", HistType::kTH2F, {{dcazAxis}, {450, 0.5f, 5.f}});
    }

    // TOF hist #DeltaMass2
    jetHist.add<TH2>("tracks/proton/h2TOFmassProtonVsPt", "h2TOFmassProtonVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 4.}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiProton/h2TOFmassantiProtonVsPt", "h2TOFmassantiProtonVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 4.}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T} (GeV)", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
    jetHist.add<TH2>("tracks/antiProton/h2TOFmass2antiProtonVsPt", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T} (GeV)", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});

    jetHist.add<TH2>("tracks/deuteron/h2TOFmassDeuteronVsPt", "h2TOFmassDeuteronVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 4.}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiDeuteron/h2TOFmassantiDeuteronVsPt", "h2TOFmassantiDeuteronVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 4.}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T} (GeV)", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
    jetHist.add<TH2>("tracks/antiDeuteron/h2TOFmass2antiDeuteronVsPt", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T} (GeV)", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});

    jetHist.add<TH2>("tracks/triton/h2TOFmassTritonVsPt", "h2TOFmassTritonVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 4.}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiTriton/h2TOFmassantiTritonVsPt", "h2TOFmassantiTritonVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 4.}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/triton/h2TOFmass2TritonVsPt", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T} (GeV)", HistType::kTH2F, {{massTrAxis}, {250, 0., 5.}});
    jetHist.add<TH2>("tracks/antiTriton/h2TOFmass2antiTritonVsPt", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T} (GeV)", HistType::kTH2F, {{massTrAxis}, {250, 0., 5.}});

    jetHist.add<TH2>("tracks/helium/h2TOFmassHeliumVsPt", "h2TOFmassHeliumVsPt; TOFmass; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {ptZHeAxis}});
    jetHist.add<TH2>("tracks/antiHelium/h2TOFmassantiHeliumVsPt", "h2TOFmassantiHeliumVsPt; TOFmass; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {ptZHeAxis}});
    jetHist.add<TH2>("tracks/helium/h2TOFmass2HeliumVsPt", "#Delta M^{2} (t) vs #it{p}_{T}t; TOFmass2; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{massHeAxis}, {ptZHeAxis}});
    jetHist.add<TH2>("tracks/antiHelium/h2TOFmass2antiHeliumVsPt", "#Delta M^{2} (t) vs #it{p}_{T}; TOFmass2; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{massHeAxis}, {ptZHeAxis}});

    // TOF hist nSigma
    jetHist.add<TH2>("tracks/proton/h2TofNsigmaProtonVsPt", "h2TofNsigmaProtonVsPt; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiProton/h2TofNsigmaantiProtonVsPt", "h2TofNsigmaantiProtonVsPt; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/deuteron/h2TofNsigmaDeuteronVsPt", "h2TofNsigmaDeuteronVsPt; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiDeuteron/h2TofNsigmaantiDeuteronVsPt", "h2TofNsigmaantiDeuteronVsPt; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/triton/h2TofNsigmaTritonVsPt", "h2TofNsigmaTritonVsPt; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiTriton/h2TofNsigmaantiTritonVsPt", "h2TofNsigmaantiTritonVsPt; TofNsigma; #it{p}_{T} (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/helium/h2TofNsigmaHeliumVsPt", "h2TofNsigmaHeliumVsPt; TofNsigma; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});
    jetHist.add<TH2>("tracks/antiHelium/h2TofNsigmaantiHeliumVsPt", "h2TofNsigmaantiHeliumVsPt; TofNsigma; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{100, -5, 5}, {50, 0., 5.}});

    if (isMC) {
      // inside jet
      jetHist.add<TH3>("tracks/mc/proton/h3PtVsProtonNSigmaTPCVsPtJet_jet", "pT(p) vs NSigmaTPC (p) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
      jetHist.add<TH3>("tracks/mc/antiProton/h3PtVsantiProtonNSigmaTPCVsPtJet_jet", "pT(#bar{p}) vs NSigmaTPC (#bar{p}) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
      jetHist.add<TH3>("tracks/mc/deuteron/h3PtVsDeuteronNSigmaTPCVsPtJet_jet", "pT(d) vs NSigmaTPC (d) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
      jetHist.add<TH3>("tracks/mc/antiDeuteron/h3PtVsantiDeuteronNSigmaTPCVsPtJet_jet", "pT(#bar{d}) vs NSigmaTPC (#bar{d}) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
      jetHist.add<TH3>("tracks/mc/helium/h3PtVsHeliumNSigmaTPCVsPtJet_jet", "pT(He) vs NSigmaTPC (He) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
      jetHist.add<TH3>("tracks/mc/antiHelium/h3PtVsantiHeliumNSigmaTPCVsPtJet_jet", "pT(#bar{He}) vs NSigmaTPC (#bar{He}) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
      jetHist.add<TH3>("tracks/mc/triton/h3PtVsTritonNSigmaTPCVsPtJet_jet", "pT(Tr) vs NSigmaTPC (Tr) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});
      jetHist.add<TH3>("tracks/mc/antiTriton/h3PtVsantiTritonNSigmaTPCVsPtJet_jet", "pT(#bar{Tr}) vs NSigmaTPC (#bar{Tr}) vs jet pT;  #it{p}_{T} (GeV/#it{c}; NSigmaTPC;  p^{jet}_{T}", HistType::kTH3F, {{PtAxis}, {200, -10, 10}, {PtJetAxis}});

      if (cEnableProtonQA) {
        jetHist.add<TH2>("tracks/proton/dca/before/hDCAxyVsPtProton_jet", "DCAxy vs Pt (p)", HistType::kTH2F, {{PtAxis}, {dcaxyAxis}});
        jetHist.add<TH2>("tracks/antiProton/dca/before/hDCAxyVsPtantiProton_jet", "DCAxy vs Pt (#bar{p})", HistType::kTH2F, {{PtAxis}, {dcaxyAxis}});
        jetHist.add<TH2>("tracks/proton/dca/before/hDCAzVsPtProton_jet", "DCAz vs Pt (p)", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
        jetHist.add<TH2>("tracks/antiProton/dca/before/hDCAzVsPtantiProton_jet", "DCAz vs Pt (#bar{p})", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
      }

      if (cEnableDeuteronQA) {
        jetHist.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtDeuteron_jet", "DCAxy vs Pt (d)", HistType::kTH2F, {{PtAxis}, {dcaxyAxis}});
        jetHist.add<TH2>("tracks/antiDeuteron/dca/before/hDCAxyVsPtantiDeuteron_jet", "DCAxy vs Pt (#bar{d})", HistType::kTH2F, {{PtAxis}, {dcaxyAxis}});
        jetHist.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtDeuteron_jet", "DCAz vs Pt (d)", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
        jetHist.add<TH2>("tracks/antiDeuteron/dca/before/hDCAzVsPtantiDeuteron_jet", "DCAz vs Pt (#bar{d})", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
      }

      if (cEnableTritonQA) {
        jetHist.add<TH2>("tracks/triton/dca/before/hDCAxyVsPtTriton_jet", "DCAxy vs Pt (t)", HistType::kTH2F, {{PtAxis}, {dcaxyAxis}});
        jetHist.add<TH2>("tracks/antiTriton/dca/before/hDCAxyVsPtantiTriton_jet", "DCAxy vs Pt (#bar{t})", HistType::kTH2F, {{PtAxis}, {dcaxyAxis}});
        jetHist.add<TH2>("tracks/triton/dca/before/hDCAzVsPtTriton_jet", "DCAz vs Pt (t)", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
        jetHist.add<TH2>("tracks/antiTriton/dca/before/hDCAzVsPtantiTriton_jet", "DCAz vs Pt (#bar{t})", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
      }

      if (cEnableHeliumQA) {
        jetHist.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtHelium_jet", "DCAxy vs Pt (He)", HistType::kTH2F, {{dcaxyAxis}, {450, 0.5f, 5.f}});
        jetHist.add<TH2>("tracks/antiHelium/dca/before/hDCAxyVsPtantiHelium_jet", "DCAxy vs Pt (#bar{He})", HistType::kTH2F, {{dcaxyAxis}, {450, 0.5f, 5.f}});
        jetHist.add<TH2>("tracks/helium/dca/before/hDCAzVsPtHelium_jet", "DCAz vs Pt (He)", HistType::kTH2F, {{dcazAxis}, {450, 0.5f, 5.f}});
        jetHist.add<TH2>("tracks/antiHelium/dca/before/hDCAzVsPtantiHelium_jet", "DCAz vs Pt (#bar{He})", HistType::kTH2F, {{dcazAxis}, {450, 0.5f, 5.f}});
      }

      // outside jet
      if (cEnableProtonQA) {
        jetHist.add<TH2>("tracks/proton/dca/before/hDCAxyVsPtProton", "DCAxy vs Pt (p)", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
        jetHist.add<TH2>("tracks/antiProton/dca/before/hDCAxyVsPtantiProton", "DCAxy vs Pt (#bar{p})", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
        jetHist.add<TH2>("tracks/proton/dca/before/hDCAzVsPtProton", "DCAz vs Pt (p)", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
        jetHist.add<TH2>("tracks/antiProton/dca/before/hDCAzVsPtantiProton", "DCAz vs Pt (#bar{p})", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
      }

      if (cEnableDeuteronQA) {
        jetHist.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtDeuteron", "DCAxy vs Pt (d)", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
        jetHist.add<TH2>("tracks/antiDeuteron/dca/before/hDCAxyVsPtantiDeuteron", "DCAxy vs Pt (#bar{d})", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
        jetHist.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtDeuteron", "DCAz vs Pt (d)", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
        jetHist.add<TH2>("tracks/antiDeuteron/dca/before/hDCAzVsPtantiDeuteron", "DCAz vs Pt (#bar{d})", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
      }

      if (cEnableTritonQA) {
        jetHist.add<TH2>("tracks/triton/dca/before/hDCAxyVsPtTriton", "DCAxy vs Pt (t)", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
        jetHist.add<TH2>("tracks/antiTriton/dca/before/hDCAxyVsPtantiTriton", "DCAxy vs Pt (#bar{t})", HistType::kTH2F, {{dcaxyAxis}, {PtAxis}});
        jetHist.add<TH2>("tracks/triton/dca/before/hDCAzVsPtTriton", "DCAz vs Pt (t)", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
        jetHist.add<TH2>("tracks/antiTriton/dca/before/hDCAzVsPtantiTriton", "DCAz vs Pt (#bar{t})", HistType::kTH2F, {{dcazAxis}, {PtAxis}});
      }

      if (cEnableHeliumQA) {
        jetHist.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtHelium", "DCAxy vs Pt (He)", HistType::kTH2F, {{450, 0.5f, 5.f}, {dcaxyAxis}});
        jetHist.add<TH2>("tracks/antiHelium/dca/before/hDCAxyVsPtantiHelium", "DCAxy vs Pt (#bar{He})", HistType::kTH2F, {{450, 0.5f, 5.f}, {dcaxyAxis}});
        jetHist.add<TH2>("tracks/helium/dca/before/hDCAzVsPtHelium", "DCAz vs Pt (He)", HistType::kTH2F, {{dcazAxis}, {450, 0.5f, 5.f}});
        jetHist.add<TH2>("tracks/antiHelium/dca/before/hDCAzVsPtantiHelium", "DCAz vs Pt (#bar{He})", HistType::kTH2F, {{dcazAxis}, {450, 0.5f, 5.f}});
      }

      // PartilceJet-constituents
      jetHist.add<TH3>("mcpJet/pt/PtParticleType", "Pt (p) vs jetflag vs particletype", HistType::kTH3D, {{100, 0.f, 10.f}, {2, 0, 2}, {14, -7, 7}});
      // detectorJet-constituents
      jetHist.add<TH3>("mcdJet/pt/PtParticleType", "Pt (p) vs jetflag vs particletype", HistType::kTH3D, {{100, 0.f, 10.f}, {2, 0, 2}, {14, -7, 7}});
      jetHist.add<TH2>("mcdJet/pt/perpCone/PtParticleType", "Pt (p) vs particletype", HistType::kTH2D, {{100, 0.f, 10.f}, {14, -7, 7}});

      jetHist.add<TH1>("mcpJet/hJetPt", "Pt (jet)", HistType::kTH1F, {{100, 0.f, 50.f}});
      jetHist.add<TH1>("mcpJet/hJetEta", "Eta (jet)", HistType::kTH1F, {{100, 1.5, 1.5}});
      jetHist.add<TH1>("mcpJet/hJetPhi", "Phi (jet)", HistType::kTH1F, {{70, 0.f, 7.f}});

      jetHist.add<TH1>("mcdJet/hJetPt", "Pt (jet)", HistType::kTH1F, {{100, 0.f, 50.f}});
      jetHist.add<TH1>("mcdJet/hJetEta", "Eta (jet)", HistType::kTH1F, {{100, 1.5, 1.5}});
      jetHist.add<TH1>("mcdJet/hJetPhi", "Phi (jet)", HistType::kTH1F, {{70, 0.f, 7.f}});

      // rec matched
      jetHist.add<TH2>("recmatched/hRecMatchedJetPt", "matched jet pT (Rec level);#it{p}_{T,jet part} (GeV/#it{c}); #it{p}_{T,jet part} - #it{p}_{T,jet det}", HistType::kTH2F, {{100, 0., 100.}, {400, -20., 20.}});
      jetHist.add<TH2>("recmatched/hRecMatchedJetPhi", "matched jet #varphi (Rec level);#varphi_{T,jet part}; #varphi_{jet part}-#varphi_{jet det}", HistType::kTH2F, {{700, 0., 7.}, {200, -5., 5.}});
      jetHist.add<TH2>("recmatched/hRecMatchedJetEta", "matched jet #eta (Rec level);#eta_{T,jet part}; #eta_{jet part}-#eta_{jet det} ", HistType::kTH2F, {{200, -1., 1.}, {500, -2.5, 2.5}});

      jetHist.add<TH2>("recmatched/h2ResponseMatrix", "matched jet pT;#it{p}_{T} (true); #it{p}_{T} (measured)", HistType::kTH2F, {{40, 0., 100.}, {40, 0., 100.}});
      /////////
      jetHist.add<TH1>("recmatched/hRecJetPt", "matched jet pT (Rec level);#it{p}_{T,jet part} (GeV/#it{c}); #it{p}_{T,jet part} - #it{p}_{T,jet det}", HistType::kTH1F, {{100, 0., 100.}});
      jetHist.add<TH1>("recmatched/hGenJetPt", "matched jet pT (Rec level);#it{p}_{T,jet part} (GeV/#it{c}); #it{p}_{T,jet part} - #it{p}_{T,jet det}", HistType::kTH1F, {{100, 0., 100.}});

      jetHist.add<TH3>("recmatched/pt/PtParticleType", "Pt (p) vs jetflag vs particletype", HistType::kTH3D, {{100, 0.f, 10.f}, {2, 0, 2}, {14, -7, 7}});

      // gen matched
      jetHist.add<TH2>("genmatched/hRecMatchedJetPt", "matched jet pT (Rec level);#it{p}_{T,jet part} (GeV/#it{c}); #it{p}_{T,jet part} - #it{p}_{T,jet det}", HistType::kTH2F, {{100, 0., 100.}, {400, -20., 20.}});
      jetHist.add<TH2>("genmatched/hRecMatchedJetPhi", "matched jet #varphi (Rec level);#varphi_{T,jet part}; #varphi_{jet part}-#varphi_{jet det}", HistType::kTH2F, {{700, 0., 7.}, {200, -5., 5.}});
      jetHist.add<TH2>("genmatched/hRecMatchedJetEta", "matched jet #eta (Rec level);#eta_{T,jet part}; #eta_{jet part}-#eta_{jet det} ", HistType::kTH2F, {{200, -1., 1.}, {500, -2.5, 2.5}});

      /////////
      jetHist.add<TH1>("genmatched/hRecJetPt", "matched jet pT (Rec level);#it{p}_{T,jet part} (GeV/#it{c}); #it{p}_{T,jet part} - #it{p}_{T,jet det}", HistType::kTH1F, {{100, 0., 100.}});
      jetHist.add<TH1>("genmatched/hGenJetPt", "matched jet pT (Rec level);#it{p}_{T,jet part} (GeV/#it{c}); #it{p}_{T,jet part} - #it{p}_{T,jet det}", HistType::kTH1F, {{100, 0., 100.}});
      jetHist.add<TH1>("genmatched/hRecJetWithGenPt", "matched jet pT (Rec level);#it{p}_{T,jet part} (GeV/#it{c}); #it{p}_{T,jet part} - #it{p}_{T,jet det}", HistType::kTH1F, {{100, 0., 100.}});
      jetHist.add<TH1>("genmatched/hGenJetPtMatched", "matched jet pT (Rec level);#it{p}_{T,jet part} (GeV/#it{c}); #it{p}_{T,jet part} - #it{p}_{T,jet det}", HistType::kTH1F, {{100, 0., 100.}});

      jetHist.add<TH3>("genmatched/pt/PtParticleType", "Pt (p) vs jetflag vs particletype", HistType::kTH3D, {{100, 0.f, 10.f}, {2, 0, 2}, {14, -7, 7}});
    }
  }

  std::array<float, 2> getPerpendicuarPhi(float jetPhi)
  {
    std::array<float, 2> PerpendicularConeAxisPhi = {-999.0f, -999.0f};
    // build 2 perp cones in phi around the leading jet (right and left of the jet)
    PerpendicularConeAxisPhi[0] = RecoDecay::constrainAngle<float, float>(jetPhi + (M_PI / 2.)); // This will contrain the angel between 0-2Pi
    PerpendicularConeAxisPhi[1] = RecoDecay::constrainAngle<float, float>(jetPhi - (M_PI / 2.)); // This will contrain the angel between 0-2Pi
    return PerpendicularConeAxisPhi;
  }

  template <typename TrackType>
  bool isTrackSelected(const TrackType track)
  {
    // standard track selection
    if (track.pt() < cfgtrkMinPt)
      return false;
    if (std::abs(track.eta()) > cfgtrkMaxEta)
      return false;
    if (std::abs(track.dcaXY()) > cfgMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cfgMaxDCAzToPVcut)
      return false;
    if (track.tpcNClsFindable() < cfgnFindableTPCClusters)
      return false;
    if (track.tpcNClsCrossedRows() < cfgnTPCCrossedRows)
      return false;
    if (track.tpcChi2NCl() > cfgnTPCChi2)
      return false;
    if (track.itsChi2NCl() > cfgnITSChi2)
      return false;
    if (cfgConnectedToPV && !track.isPVContributor())
      return false;
    return true;
  }

  int nEvents = 0;
  template <bool IsMC, typename TracksType, typename JetType>
  void fillTrackInfo(const TracksType& trk, const JetType& jets, std::vector<float>& leadingJetPtEtaPhi)
  {
    if (!isTrackSelected(trk))
      return;
    if (trk.pt() > cMaxPt)
      return;
    jetHist.fill(HIST("tracks/h2TPCsignVsTPCmomentum"), trk.tpcInnerParam() / (1.f * trk.sign()), trk.tpcSignal());
    bool jetFlag = false;
    bool jetFlagPerpCone = false;
    float jetPt = -999.;

    if (isWithLeadingJet) {
      double delPhi = TVector2::Phi_mpi_pi(leadingJetPtEtaPhi[2] - trk.phi());
      double delEta = leadingJetPtEtaPhi[1] - trk.eta();
      double R = TMath::Sqrt((delEta * delEta) + (delPhi * delPhi));
      if (R < cfgjetR)
        jetFlag = true;
      jetPt = leadingJetPtEtaPhi[0];
      // Get perpCone
      std::array<float, 2> perpConePhiJet = getPerpendicuarPhi(leadingJetPtEtaPhi[2]);
      double delPhiPerpCone1 = TVector2::Phi_mpi_pi(perpConePhiJet[0] - trk.phi());
      double delPhiPerpCone2 = TVector2::Phi_mpi_pi(perpConePhiJet[1] - trk.phi());
      double RPerpCone1 = TMath::Sqrt((delEta * delEta) + (delPhiPerpCone1 * delPhiPerpCone1));
      double RPerpCone2 = TMath::Sqrt((delEta * delEta) + (delPhiPerpCone2 * delPhiPerpCone2));
      if (RPerpCone1 < cfgjetR || RPerpCone2 < cfgjetR)
        jetFlagPerpCone = true;
    } else {
      for (auto const& jet : jets) {
        double delPhi = TVector2::Phi_mpi_pi(jet.phi() - trk.phi());
        double delEta = jet.eta() - trk.eta();
        double R = TMath::Sqrt((delEta * delEta) + (delPhi * delPhi));
        if (R < cfgjetR)
          jetFlag = true;
        jetPt = jet.pt();
        break;
      }
    }
    // tof
    // float gamma =-999;
    float massTOF = -999;
    if (trk.hasTOF()) {
      // gamma = 1.f / TMath::Sqrt(1.f - (trk.beta() * trk.beta()));
      massTOF = trk.p() * TMath::Sqrt(1.f / (trk.beta() * trk.beta()) - 1.f);
    }

    if (addTOFplots && trk.hasTOF()) {
      jetHist.fill(HIST("tracks/h2TOFbetaVsP"), trk.p() / (1.f * trk.sign()), trk.beta());
    }
    if (jetFlag) {
      jetHist.fill(HIST("tracks/h2TPCsignVsTPCmomentum_Jet"), trk.tpcInnerParam() / (1.f * trk.sign()), trk.tpcSignal());
      if (addTOFplots && trk.hasTOF()) {
        jetHist.fill(HIST("tracks/h2TOFbetaVsP_Jet"), trk.p() / (1.f * trk.sign()), trk.beta());
      }
      if (trk.sign() > 0) { // particle info
        if (addpik) {
          jetHist.fill(HIST("tracks/pion/h3PtVsPionNSigmaTPCVsPtJet_jet"), trk.pt(), trk.tpcNSigmaPi(), jetPt);
          jetHist.fill(HIST("tracks/kaon/h3PtVsKaonNSigmaTPCVsPtJet_jet"), trk.pt(), trk.tpcNSigmaKa(), jetPt);
        }
        jetHist.fill(HIST("tracks/proton/h3PtVsProtonNSigmaTPCVsPtJet_jet"), trk.pt(), trk.tpcNSigmaPr(), jetPt);
        jetHist.fill(HIST("tracks/deuteron/h3PtVsDeuteronNSigmaTPCVsPtJet_jet"), trk.pt(), trk.tpcNSigmaDe(), jetPt);
        jetHist.fill(HIST("tracks/helium/h3PtVsHeliumNSigmaTPCVsPtJet_jet"), trk.pt(), trk.tpcNSigmaHe(), jetPt);
        jetHist.fill(HIST("tracks/triton/h3PtVsTritonNSigmaTPCVsPtJet_jet"), trk.pt(), trk.tpcNSigmaTr(), jetPt);

        if (cEnableProtonQA && std::abs(trk.tpcNSigmaPr()) < cfgnTPCPIDPr) {
          jetHist.fill(HIST("tracks/proton/dca/after/hDCAxyVsPtProton_jet"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/proton/dca/after/hDCAzVsPtProton_jet"), trk.dcaZ(), trk.pt());
        }
        if (cEnableDeuteronQA && std::abs(trk.tpcNSigmaDe()) < cfgnTPCPIDDe) {
          jetHist.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtDeuteron_jet"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtDeuteron_jet"), trk.dcaZ(), trk.pt());
        }
        if (cEnableTritonQA && std::abs(trk.tpcNSigmaHe()) < cfgnTPCPIDTr) {
          jetHist.fill(HIST("tracks/triton/dca/after/hDCAxyVsPtTriton_jet"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/triton/dca/after/hDCAzVsPtTriton_jet"), trk.dcaZ(), trk.pt());
        }
        if (cEnableHeliumQA && std::abs(trk.tpcNSigmaHe()) < cfgnTPCPIDHe) {
          jetHist.fill(HIST("tracks/helium/dca/after/hDCAxyVsPtHelium_jet"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/helium/dca/after/hDCAzVsPtHelium_jet"), trk.dcaZ(), trk.pt());
        }

        if (addTOFplots && trk.hasTOF()) {
          if (!useTPCpreSel) {
            jetHist.fill(HIST("tracks/proton/h2TOFmassProtonVsPt_jet"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/proton/h2TOFmass2ProtonVsPt_jet"), massTOF * massTOF - gMassProton * gMassProton, trk.pt());

            jetHist.fill(HIST("tracks/deuteron/h2TOFmassDeuteronVsPt_jet"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/deuteron/h2TOFmass2DeuteronVsPt_jet"), massTOF * massTOF - gMassDeuteron * gMassDeuteron, trk.pt());

            jetHist.fill(HIST("tracks/triton/h2TOFmassTritonVsPt_jet"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/triton/h2TOFmass2TritonVsPt_jet"), massTOF * massTOF - gMassTriton * gMassTriton, trk.pt());

            jetHist.fill(HIST("tracks/helium/h2TOFmassHeliumVsPt_jet"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/helium/h2TOFmass2HeliumVsPt_jet"), massTOF * massTOF - gMassHelium * gMassHelium, trk.pt());

            jetHist.fill(HIST("tracks/proton/h2TofNsigmaProtonVsPt_jet"), trk.tofNSigmaPr(), trk.pt());
            jetHist.fill(HIST("tracks/deuteron/h2TofNsigmaDeuteronVsPt_jet"), trk.tofNSigmaDe(), trk.pt());
            jetHist.fill(HIST("tracks/helium/h2TofNsigmaHeliumVsPt_jet"), trk.tofNSigmaHe(), trk.pt());
            jetHist.fill(HIST("tracks/triton/h2TofNsigmaTritonVsPt_jet"), trk.tofNSigmaTr(), trk.pt());

          } else {
            if (trk.tpcNSigmaPr() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/proton/h2TOFmassProtonVsPt_jet"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/proton/h2TOFmass2ProtonVsPt_jet"), massTOF * massTOF - gMassProton * gMassProton, trk.pt());
              jetHist.fill(HIST("tracks/proton/h2TofNsigmaProtonVsPt_jet"), trk.tofNSigmaPr(), trk.pt());
            }

            if (trk.tpcNSigmaDe() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/deuteron/h2TOFmassDeuteronVsPt_jet"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/deuteron/h2TOFmass2DeuteronVsPt_jet"), massTOF * massTOF - gMassDeuteron * gMassDeuteron, trk.pt());
              jetHist.fill(HIST("tracks/deuteron/h2TofNsigmaDeuteronVsPt_jet"), trk.tofNSigmaDe(), trk.pt());
            }
            if (trk.tpcNSigmaTr() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/triton/h2TOFmassTritonVsPt_jet"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/triton/h2TOFmass2TritonVsPt_jet"), massTOF * massTOF - gMassTriton * gMassTriton, trk.pt());
              jetHist.fill(HIST("tracks/helium/h2TofNsigmaHeliumVsPt_jet"), trk.tofNSigmaHe(), trk.pt());
            }
            if (trk.tpcNSigmaHe() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/helium/h2TOFmassHeliumVsPt_jet"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/helium/h2TOFmass2HeliumVsPt_jet"), massTOF * massTOF - gMassHelium * gMassHelium, trk.pt());
              jetHist.fill(HIST("tracks/triton/h2TofNsigmaTritonVsPt_jet"), trk.tofNSigmaTr(), trk.pt());
            }
          }
          // nSigma
          if (addpik) {
            jetHist.fill(HIST("tracks/pion/h2TofNsigmaPionVsPt_jet"), trk.tofNSigmaPi(), trk.pt());
            jetHist.fill(HIST("tracks/kaon/h2TofNsigmaKaonVsPt_jet"), trk.tofNSigmaKa(), trk.pt());
          }

        } // tof info

      } else { // anti-particle info
        if (addpik) {
          jetHist.fill(HIST("tracks/antiPion/h3PtVsPionNSigmaTPCVsPtJet_jet"), trk.pt(), trk.tpcNSigmaPi(), jetPt);
          jetHist.fill(HIST("tracks/antiKaon/h3PtVsKaonNSigmaTPCVsPtJet_jet"), trk.pt(), trk.tpcNSigmaKa(), jetPt);
        }
        jetHist.fill(HIST("tracks/antiProton/h3PtVsantiProtonNSigmaTPCVsPtJet_jet"), trk.pt(), trk.tpcNSigmaPr(), jetPt);
        jetHist.fill(HIST("tracks/antiDeuteron/h3PtVsantiDeuteronNSigmaTPCVsPtJet_jet"), trk.pt(), trk.tpcNSigmaDe(), jetPt);
        jetHist.fill(HIST("tracks/antiHelium/h3PtVsantiHeliumNSigmaTPCVsPtJet_jet"), trk.pt(), trk.tpcNSigmaHe(), jetPt);
        jetHist.fill(HIST("tracks/antiTriton/h3PtVsantiTritonNSigmaTPCVsPtJet_jet"), trk.pt(), trk.tpcNSigmaTr(), jetPt);

        if (cEnableProtonQA && std::abs(trk.tpcNSigmaPr()) < cfgnTPCPIDPr) {
          jetHist.fill(HIST("tracks/antiProton/dca/after/hDCAxyVsPtantiProton_jet"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/antiProton/dca/after/hDCAzVsPtantiProton_jet"), trk.dcaZ(), trk.pt());
        }
        if (cEnableDeuteronQA && std::abs(trk.tpcNSigmaDe()) < cfgnTPCPIDDe) {
          jetHist.fill(HIST("tracks/antiDeuteron/dca/after/hDCAxyVsPtantiDeuteron_jet"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/antiDeuteron/dca/after/hDCAzVsPtantiDeuteron_jet"), trk.dcaZ(), trk.pt());
        }
        if (cEnableHeliumQA && std::abs(trk.tpcNSigmaHe()) < cfgnTPCPIDHe) {
          jetHist.fill(HIST("tracks/antiTriton/dca/after/hDCAxyVsPtantiTriton_jet"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/antiTriton/dca/after/hDCAzVsPtantiTriton_jet"), trk.dcaZ(), trk.pt());
        }
        if (cEnableTritonQA && std::abs(trk.tpcNSigmaHe()) < cfgnTPCPIDHe) {
          jetHist.fill(HIST("tracks/antiHelium/dca/after/hDCAxyVsPtantiHelium_jet"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/antiHelium/dca/after/hDCAzVsPtantiHelium_jet"), trk.dcaZ(), trk.pt());
        }

        if (addTOFplots && trk.hasTOF()) {
          if (!useTPCpreSel) {
            jetHist.fill(HIST("tracks/antiProton/h2TOFmassantiProtonVsPt_jet"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/antiProton/h2TOFmass2antiProtonVsPt_jet"), massTOF * massTOF - gMassProton * gMassProton, trk.pt());
            jetHist.fill(HIST("tracks/antiProton/h2TofNsigmaantiProtonVsPt_jet"), trk.tofNSigmaPr(), trk.pt());

            jetHist.fill(HIST("tracks/antiDeuteron/h2TOFmassantiDeuteronVsPt_jet"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/antiDeuteron/h2TOFmass2antiDeuteronVsPt_jet"), massTOF * massTOF - gMassDeuteron * gMassDeuteron, trk.pt());
            jetHist.fill(HIST("tracks/antiDeuteron/h2TofNsigmaantiDeuteronVsPt_jet"), trk.tofNSigmaDe(), trk.pt());

            jetHist.fill(HIST("tracks/antiTriton/h2TOFmassantiTritonVsPt_jet"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/antiTriton/h2TOFmass2antiTritonVsPt_jet"), massTOF * massTOF - gMassTriton * gMassTriton, trk.pt());
            jetHist.fill(HIST("tracks/antiTriton/h2TofNsigmaantiTritonVsPt_jet"), trk.tofNSigmaTr(), trk.pt());

            jetHist.fill(HIST("tracks/antiHelium/h2TOFmassantiHeliumVsPt_jet"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/antiHelium/h2TOFmass2antiHeliumVsPt_jet"), massTOF * massTOF - gMassHelium * gMassHelium, trk.pt());
            jetHist.fill(HIST("tracks/antiHelium/h2TofNsigmaantiHeliumVsPt_jet"), trk.tofNSigmaHe(), trk.pt());
          } else {
            if (trk.tpcNSigmaPr() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/antiProton/h2TOFmassantiProtonVsPt_jet"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/antiProton/h2TOFmass2antiProtonVsPt_jet"), massTOF * massTOF - gMassProton * gMassProton, trk.pt());
              jetHist.fill(HIST("tracks/antiProton/h2TofNsigmaantiProtonVsPt_jet"), trk.tofNSigmaPr(), trk.pt());
            }
            if (trk.tpcNSigmaDe() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/antiDeuteron/h2TOFmassantiDeuteronVsPt_jet"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/antiDeuteron/h2TOFmass2antiDeuteronVsPt_jet"), massTOF * massTOF - gMassDeuteron * gMassDeuteron, trk.pt());
              jetHist.fill(HIST("tracks/antiDeuteron/h2TofNsigmaantiDeuteronVsPt_jet"), trk.tofNSigmaDe(), trk.pt());
            }
            if (trk.tpcNSigmaTr() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/antiTriton/h2TOFmassantiTritonVsPt_jet"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/antiTriton/h2TOFmass2antiTritonVsPt_jet"), massTOF * massTOF - gMassTriton * gMassTriton, trk.pt());
              jetHist.fill(HIST("tracks/antiTriton/h2TofNsigmaantiTritonVsPt_jet"), trk.tofNSigmaTr(), trk.pt());
            }
            if (trk.tpcNSigmaHe() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/antiHelium/h2TOFmassantiHeliumVsPt_jet"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/antiHelium/h2TOFmass2antiHeliumVsPt_jet"), massTOF * massTOF - gMassHelium * gMassHelium, trk.pt());
              jetHist.fill(HIST("tracks/antiHelium/h2TofNsigmaantiHeliumVsPt_jet"), trk.tofNSigmaHe(), trk.pt());
            }
          }

          if (addpik) {
            if (!useTPCpreSel) {
              jetHist.fill(HIST("tracks/antiPion/h2TofNsigmaantiPionVsPt_jet"), trk.tofNSigmaPi(), trk.pt());
              jetHist.fill(HIST("tracks/antiKaon/h2TofNsigmaantiKaonVsPt_jet"), trk.tofNSigmaKa(), trk.pt());
            } else {
              if (trk.tpcNSigmaPi() < useTPCpreSel) {
                jetHist.fill(HIST("tracks/antiPion/h2TofNsigmaantiPionVsPt_jet"), trk.tofNSigmaPi(), trk.pt());
              }
              if (trk.tpcNSigmaKa() < useTPCpreSel) {
                jetHist.fill(HIST("tracks/antiKaon/h2TofNsigmaantiKaonVsPt_jet"), trk.tofNSigmaKa(), trk.pt());
              }
            }
          } // pikEnd
        }
      } // anti-particle
        ////////////////////////////////////////
      // within jet end
      //////////////////////////////////////////
    } else {
      jetHist.fill(HIST("tracks/h2TPCsignVsTPCmomentum_OutJet"), trk.tpcInnerParam() / (1.f * trk.sign()), trk.tpcSignal());
      if (jetFlagPerpCone && isWithLeadingJet) {
        jetHist.fill(HIST("tracks/perpCone/h2TPCsignVsTPCmomentum"), trk.tpcInnerParam() / (1.f * trk.sign()), trk.tpcSignal());
      }
      if (addTOFplots && trk.hasTOF()) {
        jetHist.fill(HIST("tracks/h2TOFbetaVsP_OutJet"), trk.p() / (1.f * trk.sign()), trk.beta());
      }
      if (trk.sign() > 0) {
        jetHist.fill(HIST("tracks/proton/h3PtVsProtonNSigmaTPC"), trk.pt(), trk.tpcNSigmaPr());     // Pr
        jetHist.fill(HIST("tracks/deuteron/h3PtVsDeuteronNSigmaTPC"), trk.pt(), trk.tpcNSigmaDe()); // De
        jetHist.fill(HIST("tracks/helium/h3PtVsHeliumNSigmaTPC"), trk.pt(), trk.tpcNSigmaHe());     // He
        jetHist.fill(HIST("tracks/triton/h3PtVsTritonNSigmaTPC"), trk.pt(), trk.tpcNSigmaTr());     // Tr
        // perpCone
        if (jetFlagPerpCone && isWithLeadingJet) {
          jetHist.fill(HIST("tracks/perpCone/proton/h3PtVsProtonNSigmaTPCVsPtJet"), trk.pt(), trk.tpcNSigmaPr(), jetPt);     // Pr
          jetHist.fill(HIST("tracks/perpCone/deuteron/h3PtVsDeuteronNSigmaTPCVsPtJet"), trk.pt(), trk.tpcNSigmaDe(), jetPt); // De
          jetHist.fill(HIST("tracks/perpCone/helium/h3PtVsHeliumNSigmaTPCVsPtJet"), trk.pt(), trk.tpcNSigmaHe(), jetPt);     // He
          jetHist.fill(HIST("tracks/perpCone/triton/h3PtVsTritonNSigmaTPCVsPtJet"), trk.pt(), trk.tpcNSigmaTr(), jetPt);     // Tr
        }

        if (cEnableProtonQA && std::abs(trk.tpcNSigmaPr()) < cfgnTPCPIDPr) {
          jetHist.fill(HIST("tracks/proton/dca/after/hDCAxyVsPtProton"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/proton/dca/after/hDCAzVsPtProton"), trk.dcaZ(), trk.pt());
        }
        if (cEnableDeuteronQA && std::abs(trk.tpcNSigmaDe()) < cfgnTPCPIDDe) {
          jetHist.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtDeuteron"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtDeuteron"), trk.dcaZ(), trk.pt());
        }
        if (cEnableTritonQA && std::abs(trk.tpcNSigmaHe()) < cfgnTPCPIDTr) {
          jetHist.fill(HIST("tracks/triton/dca/after/hDCAxyVsPtTriton"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/triton/dca/after/hDCAzVsPtTriton"), trk.dcaZ(), trk.pt());
        }
        if (cEnableHeliumQA && std::abs(trk.tpcNSigmaHe()) < cfgnTPCPIDHe) {
          jetHist.fill(HIST("tracks/helium/dca/after/hDCAxyVsPtHelium"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/helium/dca/after/hDCAzVsPtHelium"), trk.dcaZ(), trk.pt());
        }
        if (addTOFplots && trk.hasTOF()) {
          if (!useTPCpreSel) {
            jetHist.fill(HIST("tracks/proton/h2TOFmassProtonVsPt"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - gMassProton * gMassProton, trk.pt());
            jetHist.fill(HIST("tracks/proton/h2TofNsigmaProtonVsPt"), trk.tofNSigmaPr(), trk.pt());

            jetHist.fill(HIST("tracks/deuteron/h2TOFmassDeuteronVsPt"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - gMassDeuteron * gMassDeuteron, trk.pt());
            jetHist.fill(HIST("tracks/deuteron/h2TofNsigmaDeuteronVsPt"), trk.tofNSigmaDe(), trk.pt());

            jetHist.fill(HIST("tracks/triton/h2TOFmassTritonVsPt"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/triton/h2TOFmass2TritonVsPt"), massTOF * massTOF - gMassTriton * gMassTriton, trk.pt());
            jetHist.fill(HIST("tracks/helium/h2TofNsigmaHeliumVsPt"), trk.tofNSigmaHe(), trk.pt());

            jetHist.fill(HIST("tracks/helium/h2TOFmassHeliumVsPt"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/helium/h2TOFmass2HeliumVsPt"), massTOF * massTOF - gMassHelium * gMassHelium, trk.pt());
            jetHist.fill(HIST("tracks/triton/h2TofNsigmaTritonVsPt"), trk.tofNSigmaTr(), trk.pt());
          } else {
            if (trk.tpcNSigmaPr() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/proton/h2TOFmassProtonVsPt"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - gMassProton * gMassProton, trk.pt());
              jetHist.fill(HIST("tracks/proton/h2TofNsigmaProtonVsPt"), trk.tofNSigmaPr(), trk.pt());
              if (jetFlagPerpCone && isWithLeadingJet)
                jetHist.fill(HIST("tracks/perpCone/proton/h2TofNsigmaProtonVsPt"), trk.tofNSigmaPr(), trk.pt());
            }
            if (trk.tpcNSigmaDe() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/deuteron/h2TOFmassDeuteronVsPt"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - gMassDeuteron * gMassDeuteron, trk.pt());
              jetHist.fill(HIST("tracks/deuteron/h2TofNsigmaDeuteronVsPt"), trk.tofNSigmaDe(), trk.pt());
              if (jetFlagPerpCone && isWithLeadingJet)
                jetHist.fill(HIST("tracks/perpCone/deuteron/h2TofNsigmaDeuteronVsPt"), trk.tofNSigmaDe(), trk.pt());
            }
            if (trk.tpcNSigmaTr() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/triton/h2TOFmassTritonVsPt"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/triton/h2TOFmass2TritonVsPt"), massTOF * massTOF - gMassTriton * gMassTriton, trk.pt());
              jetHist.fill(HIST("tracks/triton/h2TofNsigmaTritonVsPt"), trk.tofNSigmaTr(), trk.pt());
              if (jetFlagPerpCone && isWithLeadingJet)
                jetHist.fill(HIST("tracks/perpCone/triton/h2TofNsigmaTritonVsPt"), trk.tofNSigmaTr(), trk.pt());
            }
            if (trk.tpcNSigmaHe() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/helium/h2TOFmassHeliumVsPt"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/helium/h2TOFmass2HeliumVsPt"), massTOF * massTOF - gMassHelium * gMassHelium, trk.pt());
              jetHist.fill(HIST("tracks/helium/h2TofNsigmaHeliumVsPt"), trk.tofNSigmaHe(), trk.pt());
              if (jetFlagPerpCone && isWithLeadingJet)
                jetHist.fill(HIST("tracks/perpCone/helium/h2TofNsigmaHeliumVsPt"), trk.tofNSigmaHe(), trk.pt());
            }
          }

        } // tof info
      } else {
        jetHist.fill(HIST("tracks/antiProton/h3PtVsantiProtonNSigmaTPC"), trk.pt(), trk.tpcNSigmaPr());     // Pr
        jetHist.fill(HIST("tracks/antiDeuteron/h3PtVsantiDeuteronNSigmaTPC"), trk.pt(), trk.tpcNSigmaDe()); // De
        jetHist.fill(HIST("tracks/antiHelium/h3PtVsantiHeliumNSigmaTPC"), trk.pt(), trk.tpcNSigmaHe());     // He
        jetHist.fill(HIST("tracks/antiTriton/h3PtVsantiTritonNSigmaTPC"), trk.pt(), trk.tpcNSigmaTr());     // Tr

        // perpCone
        if (jetFlagPerpCone && isWithLeadingJet) {
          // antiparticle info
          jetHist.fill(HIST("tracks/perpCone/antiProton/h3PtVsantiProtonNSigmaTPCVsPtJet"), trk.pt(), trk.tpcNSigmaPr(), jetPt);     // Pr
          jetHist.fill(HIST("tracks/perpCone/antiDeuteron/h3PtVsantiDeuteronNSigmaTPCVsPtJet"), trk.pt(), trk.tpcNSigmaDe(), jetPt); // De
          jetHist.fill(HIST("tracks/perpCone/antiHelium/h3PtVsantiHeliumNSigmaTPCVsPtJet"), trk.pt(), trk.tpcNSigmaHe(), jetPt);     // He
          jetHist.fill(HIST("tracks/perpCone/antiTriton/h3PtVsantiTritonNSigmaTPCVsPtJet"), trk.pt(), trk.tpcNSigmaTr(), jetPt);     // Tr
        }

        if (cEnableProtonQA && std::abs(trk.tpcNSigmaPr()) < cfgnTPCPIDPr) {
          jetHist.fill(HIST("tracks/antiProton/dca/after/hDCAxyVsPtantiProton"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/antiProton/dca/after/hDCAzVsPtantiProton"), trk.dcaZ(), trk.pt());
        }
        if (cEnableDeuteronQA && std::abs(trk.tpcNSigmaDe()) < cfgnTPCPIDDe) {
          jetHist.fill(HIST("tracks/antiDeuteron/dca/after/hDCAxyVsPtantiDeuteron"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/antiDeuteron/dca/after/hDCAzVsPtantiDeuteron"), trk.dcaZ(), trk.pt());
        }
        if (cEnableHeliumQA && std::abs(trk.tpcNSigmaHe()) < cfgnTPCPIDHe) {
          jetHist.fill(HIST("tracks/antiTriton/dca/after/hDCAxyVsPtantiTriton"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/antiTriton/dca/after/hDCAzVsPtantiTriton"), trk.dcaZ(), trk.pt());
        }
        if (cEnableTritonQA && std::abs(trk.tpcNSigmaHe()) < cfgnTPCPIDHe) {
          jetHist.fill(HIST("tracks/antiHelium/dca/after/hDCAxyVsPtantiHelium"), trk.dcaXY(), trk.pt());
          jetHist.fill(HIST("tracks/antiHelium/dca/after/hDCAzVsPtantiHelium"), trk.dcaZ(), trk.pt());
        }

        if (addTOFplots && trk.hasTOF()) {
          if (!useTPCpreSel) {
            jetHist.fill(HIST("tracks/antiProton/h2TOFmassantiProtonVsPt"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/antiProton/h2TOFmass2antiProtonVsPt"), massTOF * massTOF - gMassProton * gMassProton, trk.pt());
            jetHist.fill(HIST("tracks/antiProton/h2TofNsigmaantiProtonVsPt"), trk.tofNSigmaPr(), trk.pt());

            jetHist.fill(HIST("tracks/antiDeuteron/h2TOFmassantiDeuteronVsPt"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/antiDeuteron/h2TOFmass2antiDeuteronVsPt"), massTOF * massTOF - gMassDeuteron * gMassDeuteron, trk.pt());
            jetHist.fill(HIST("tracks/antiDeuteron/h2TofNsigmaantiDeuteronVsPt"), trk.tofNSigmaDe(), trk.pt());

            jetHist.fill(HIST("tracks/antiTriton/h2TOFmassantiTritonVsPt"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/antiTriton/h2TOFmass2antiTritonVsPt"), massTOF * massTOF - gMassTriton * gMassTriton, trk.pt());
            jetHist.fill(HIST("tracks/antiHelium/h2TofNsigmaantiHeliumVsPt"), trk.tofNSigmaHe(), trk.pt());

            jetHist.fill(HIST("tracks/antiHelium/h2TOFmassantiHeliumVsPt"), massTOF, trk.pt());
            jetHist.fill(HIST("tracks/antiHelium/h2TOFmass2antiHeliumVsPt"), massTOF * massTOF - gMassHelium * gMassHelium, trk.pt());
            jetHist.fill(HIST("tracks/antiTriton/h2TofNsigmaantiTritonVsPt"), trk.tofNSigmaTr(), trk.pt());
          } else {
            if (trk.tpcNSigmaPr() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/antiProton/h2TOFmassantiProtonVsPt"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/antiProton/h2TOFmass2antiProtonVsPt"), massTOF * massTOF - gMassProton * gMassProton, trk.pt());
              jetHist.fill(HIST("tracks/antiProton/h2TofNsigmaantiProtonVsPt"), trk.tofNSigmaPr(), trk.pt());
              if (jetFlagPerpCone && isWithLeadingJet)
                jetHist.fill(HIST("tracks/perpCone/antiProton/h2TofNsigmaantiProtonVsPt"), trk.tofNSigmaPr(), trk.pt());
            }
            if (trk.tpcNSigmaDe() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/antiDeuteron/h2TOFmassantiDeuteronVsPt"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/antiDeuteron/h2TOFmass2antiDeuteronVsPt"), massTOF * massTOF - gMassDeuteron * gMassDeuteron, trk.pt());
              jetHist.fill(HIST("tracks/antiDeuteron/h2TofNsigmaantiDeuteronVsPt"), trk.tofNSigmaDe(), trk.pt());
              if (jetFlagPerpCone && isWithLeadingJet)
                jetHist.fill(HIST("tracks/perpCone/antiDeuteron/h2TofNsigmaantiDeuteronVsPt"), trk.tofNSigmaDe(), trk.pt());
            }
            if (trk.tpcNSigmaTr() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/antiTriton/h2TOFmassantiTritonVsPt"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/antiTriton/h2TOFmass2antiTritonVsPt"), massTOF * massTOF - gMassTriton * gMassTriton, trk.pt());
              jetHist.fill(HIST("tracks/antiTriton/h2TofNsigmaantiTritonVsPt"), trk.tofNSigmaTr(), trk.pt());
              if (jetFlagPerpCone && isWithLeadingJet)
                jetHist.fill(HIST("tracks/perpCone/antiTriton/h2TofNsigmaantiTritonVsPt"), trk.tofNSigmaTr(), trk.pt());
            }
            if (trk.tpcNSigmaHe() < useTPCpreSel) {
              jetHist.fill(HIST("tracks/antiHelium/h2TOFmassantiHeliumVsPt"), massTOF, trk.pt());
              jetHist.fill(HIST("tracks/antiHelium/h2TOFmass2antiHeliumVsPt"), massTOF * massTOF - gMassHelium * gMassHelium, trk.pt());
              jetHist.fill(HIST("tracks/antiHelium/h2TofNsigmaantiHeliumVsPt"), trk.tofNSigmaHe(), trk.pt());
              if (jetFlagPerpCone && isWithLeadingJet)
                jetHist.fill(HIST("tracks/perpCone/antiHelium/h2TofNsigmaantiHeliumVsPt"), trk.tofNSigmaHe(), trk.pt());
            }
          }
        }
      }

    } ////////////////////////////////////////
      // outside jet end
    ////////////////////////////////////////
  }

  void processJetTracksData(aod::JetCollision const& collision, chargedJetstrack const& chargedjets, soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks, TrackCandidates const&)
  {

    if (fabs(collision.posZ()) > 10)
      return;
    if (!jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel8")))
      return;

    int nJets = 0;
    std::vector<float> leadingJetWithPtEtaPhi(3);
    float leadingJetPt = -1.0f;
    for (const auto& chargedjet : chargedjets) {
      jetHist.fill(HIST("jet/h1JetPt"), chargedjet.pt());
      jetHist.fill(HIST("jet/h1JetEta"), chargedjet.eta());
      jetHist.fill(HIST("jet/h1JetPhi"), chargedjet.phi());

      if (chargedjet.pt() > leadingJetPt) {
        leadingJetWithPtEtaPhi[0] = chargedjet.pt();
        leadingJetWithPtEtaPhi[1] = chargedjet.eta();
        leadingJetWithPtEtaPhi[2] = chargedjet.phi();
      }
      nJets++;
    }

    jetHist.fill(HIST("jet/nJetsPerEvent"), nJets);
    jetHist.fill(HIST("vertexZ"), collision.posZ());

    if (nJets > 0)
      jetHist.fill(HIST("jet/vertexZ"), collision.posZ());
    else
      jetHist.fill(HIST("jetOut/vertexZ"), collision.posZ());

    if (isWithJetEvents && nJets == 0)
      return;
    jetHist.fill(HIST("jet/h1JetEvents"), 0.5);

    for (auto& track : tracks) {
      if (useLfTpcPid) {
        auto trk = track.track_as<TrackCandidatesLfPid>();
        fillTrackInfo<false>(trk, chargedjets, leadingJetWithPtEtaPhi);
      } else {
        auto trk = track.track_as<TrackCandidates>();
        fillTrackInfo<false>(trk, chargedjets, leadingJetWithPtEtaPhi);
      }
    }
  }

  void processMCGen(o2::aod::JetMcCollision const& collision, /*soa::SmallGroups<soa::Join<aod::JMcCollisionLbs, aod::JCollisions>> const& recoColls,*/ aod::JetParticles const& mcParticles, soa::Filtered<aod::ChargedMCParticleLevelJets> const& mcpjets)
  {
    jetHist.fill(HIST("mcpJet/eventStat"), 0.5);
    jetHist.fill(HIST("mcpJet/eventStat"), 1.5);

    if (fabs(collision.posZ()) > 10) // bad vertex
      return;

    jetHist.fill(HIST("mcpJet/eventStat"), 2.5);

    jetHist.fill(HIST("mcpJet/vertexZ"), collision.posZ());

    bool INELgt0 = false;
    for (const auto& mcParticle : mcParticles) {
      if (fabs(mcParticle.eta()) < cfgtrkMaxEta) {
        INELgt0 = true;
        break;
      }
    }
    if (!INELgt0) // not true INEL
      return;

    jetHist.fill(HIST("mcpJet/eventStat"), 3.5);

    int nJets = 0;
    for (auto& mcpjet : mcpjets) {
      jetHist.fill(HIST("mcpJet/hJetPt"), mcpjet.pt());
      jetHist.fill(HIST("mcpJet/hJetEta"), mcpjet.eta());
      jetHist.fill(HIST("mcpJet/hJetPhi"), mcpjet.phi());
      nJets++;
    }
    jetHist.fill(HIST("mcpJet/nJetsPerEvent"), nJets);

    for (const auto& mcParticle : mcParticles) {

      if (!mcParticle.isPhysicalPrimary())
        continue;
      if (fabs(mcParticle.eta()) > cfgtrkMaxEta)
        continue;
      if (fabs(mcParticle.y()) > cfgtrkMaxRap)
        continue;

      bool jetFlag = false;
      // float jetPt = -999.;
      for (auto& mcpjet : mcpjets) {
        double delPhi = TVector2::Phi_mpi_pi(mcpjet.phi() - mcParticle.phi());
        double delEta = mcpjet.eta() - mcParticle.eta();
        double R = TMath::Sqrt((delEta * delEta) + (delPhi * delPhi));
        if (R < cfgjetR)
          jetFlag = true;
        // jetPt = mcpjet.pt();
        break;
      } // jet
      if (mapPDGToValue(mcParticle.pdgCode()) != 0) {
        jetHist.fill(HIST("mcpJet/pt/PtParticleType"), mcParticle.pt(), jetFlag, mapPDGToValue(mcParticle.pdgCode()));
      }

    } // track
  } // process mc

  void processMCRec(o2::aod::JetCollision const& collisionJet, soa::Join<aod::JetTracks, aod::JTrackPIs, aod::JMcTrackLbs> const& tracks,
                    soa::Filtered<aod::ChargedMCDetectorLevelJets> const& mcdjets, TrackCandidatesMC const&, aod::JetParticles const&)
  {
    jetHist.fill(HIST("mcdJet/eventStat"), 0.5);
    // JEhistos.fill(HIST("nEvents_MCRec"), 0.5);

    if (!jetderiveddatautilities::selectCollision(collisionJet, jetderiveddatautilities::initialiseEventSelectionBits("sel8")))
      return;
    // bool jetFlag = kFALSE;

    jetHist.fill(HIST("mcdJet/eventStat"), 1.5);

    if (fabs(collisionJet.posZ()) > 10)
      return;

    jetHist.fill(HIST("mcdJet/eventStat"), 2.5);
    bool INELgt0 = false;
    for (const auto& track : tracks) {
      if (fabs(track.eta()) < cfgtrkMaxEta) {
        INELgt0 = true;
        break;
      }
    }
    if (!INELgt0)
      return;
    jetHist.fill(HIST("mcdJet/eventStat"), 3.5);
    int nJets = 0;
    std::vector<float> leadingJetWithPtEtaPhi(3);
    float leadingJetPt = -1.0f;
    for (auto& mcdjet : mcdjets) {
      jetHist.fill(HIST("mcdJet/hJetPt"), mcdjet.pt());
      jetHist.fill(HIST("mcdJet/hJetEta"), mcdjet.eta());
      jetHist.fill(HIST("mcdJet/hJetPhi"), mcdjet.phi());
      if (mcdjet.pt() > leadingJetPt) {
        leadingJetWithPtEtaPhi[0] = mcdjet.pt();
        leadingJetWithPtEtaPhi[1] = mcdjet.eta();
        leadingJetWithPtEtaPhi[2] = mcdjet.phi();
      }
      nJets++;
    }

    jetHist.fill(HIST("mcdJet/vertexZ"), collisionJet.posZ());
    jetHist.fill(HIST("mcdJet/nJetsPerEvent"), nJets);

    if (isWithJetEvents && nJets == 0)
      return;

    for (auto& track : tracks) {
      auto fullTrack = track.track_as<TrackCandidatesMC>();
      if (!isTrackSelected(fullTrack))
        continue;
      if (!track.has_mcParticle())
        continue;
      auto mcTrack = track.mcParticle_as<aod::JetParticles>();
      if (fabs(mcTrack.eta()) > cfgtrkMaxEta)
        continue;
      if (!mcTrack.isPhysicalPrimary())
        continue;

      bool jetFlag = false;
      bool jetFlagPerpCone = false;
      // float jetPt = -999.;
      if (isWithLeadingJet) {
        double delPhi = TVector2::Phi_mpi_pi(leadingJetWithPtEtaPhi[2] - track.phi());
        double delEta = leadingJetWithPtEtaPhi[1] - track.eta();
        double R = TMath::Sqrt((delEta * delEta) + (delPhi * delPhi));
        if (R < cfgjetR)
          jetFlag = true;
        std::array<float, 2> perpConePhiJet = getPerpendicuarPhi(leadingJetWithPtEtaPhi[2]);
        double delPhiPerpCone1 = TVector2::Phi_mpi_pi(perpConePhiJet[0] - track.phi());
        double delPhiPerpCone2 = TVector2::Phi_mpi_pi(perpConePhiJet[1] - track.phi());
        double RPerpCone1 = TMath::Sqrt((delEta * delEta) + (delPhiPerpCone1 * delPhiPerpCone1));
        double RPerpCone2 = TMath::Sqrt((delEta * delEta) + (delPhiPerpCone2 * delPhiPerpCone2));
        if (RPerpCone1 < cfgjetR || RPerpCone2 < cfgjetR)
          jetFlagPerpCone = true;
      } else {
        for (auto& mcdjet : mcdjets) {
          double delPhi = TVector2::Phi_mpi_pi(mcdjet.phi() - track.phi());
          double delEta = mcdjet.eta() - track.eta();
          double R = TMath::Sqrt((delEta * delEta) + (delPhi * delPhi));
          if (R < cfgjetR)
            jetFlag = true;
          //  jetPt = mcdjet.pt();
          break;
        } // jet
      }

      if (mapPDGToValue(mcTrack.pdgCode()) != 0) {
        jetHist.fill(HIST("mcdJet/pt/PtParticleType"), mcTrack.pt(), jetFlag, mapPDGToValue(mcTrack.pdgCode()));
        if (jetFlagPerpCone)
          jetHist.fill(HIST("mcdJet/pt/perpCone/PtParticleType"), mcTrack.pt(), mapPDGToValue(mcTrack.pdgCode()));
      }

    } // tracks
  }

  void processRecMatched(aod::JetCollision const& collision, JetMCDetTable const& mcdjets,
                         soa::Join<aod::JetTracks, aod::JTrackPIs, aod::JMcTrackLbs> const& tracks,
                         JetMCPartTable const&, TrackCandidatesMC const&, aod::JetParticles const&)
  {
    if (fabs(collision.posZ()) > 10)
      return;
    if (!jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel8")))
      return;

    jetHist.fill(HIST("recmatched/vertexZ"), collision.posZ());
    bool INELgt0 = false;
    for (const auto& track : tracks) {
      if (fabs(track.eta()) < cfgtrkMaxEta) {
        INELgt0 = true;
        break;
      }
    }
    if (!INELgt0)
      return;

    std::vector<double> mcdJetPt{};
    std::vector<double> mcdJetPhi{};
    std::vector<double> mcdJetEta{};
    std::vector<double> mcpJetPt{};
    std::vector<double> mcpJetPhi{};
    std::vector<double> mcpJetEta{};

    for (auto& mcdjet : mcdjets) {

      if (!mcdjet.has_matchedJetGeo())
        continue;
      for (auto& mcpjet : mcdjet.template matchedJetGeo_as<JetMCPartTable>()) {
        if (!mcpjet.has_matchedJetGeo())
          continue;

        mcdJetPt.push_back(mcdjet.pt());
        mcdJetPhi.push_back(mcdjet.phi());
        mcdJetEta.push_back(mcdjet.eta());
        mcpJetPt.push_back(mcpjet.pt());
        mcpJetPhi.push_back(mcpjet.phi());
        mcpJetEta.push_back(mcpjet.eta());

        jetHist.fill(HIST("recmatched/hRecMatchedJetPt"), mcpjet.pt(), mcpjet.pt() - mcdjet.pt());
        jetHist.fill(HIST("recmatched/hRecMatchedJetPhi"), mcpjet.phi(), mcpjet.phi() - mcdjet.phi());
        jetHist.fill(HIST("recmatched/hRecMatchedJetEta"), mcpjet.eta(), mcpjet.eta() - mcdjet.eta());

        jetHist.fill(HIST("recmatched/hRecJetPt"), mcdjet.pt());
        jetHist.fill(HIST("recmatched/hGenJetPt"), mcpjet.pt());
        jetHist.fill(HIST("recmatched/h2ResponseMatrix"), mcpjet.pt(), mcdjet.pt());

      } // mcpJet

    } // mcdJet

    for (const auto& track : tracks) {
      auto completeTrack = track.track_as<TrackCandidatesMC>();
      if (fabs(completeTrack.eta()) > cfgtrkMaxEta)
        continue;
      if (!isTrackSelected(completeTrack))
        continue;
      if (!track.has_mcParticle())
        continue;
      auto mcTrack = track.mcParticle_as<aod::JetParticles>();
      // add pid later

      bool jetFlag = false;
      for (std::size_t iDJet = 0; iDJet < mcdJetPt.size(); iDJet++) {
        double delPhi = TVector2::Phi_mpi_pi(mcdJetPhi[iDJet] - track.phi());
        double delEta = mcdJetEta[iDJet] - track.eta();
        double R = TMath::Sqrt((delEta * delEta) + (delPhi * delPhi));

        if (R < cfgjetR) {
          jetFlag = true;
          break;
        }
      }

      if (mapPDGToValue(mcTrack.pdgCode()) != 0) {
        jetHist.fill(HIST("recmatched/pt/PtParticleType"), mcTrack.pt(), jetFlag, mapPDGToValue(mcTrack.pdgCode()));
      }
    } // tracks
  } // process

  int nprocessSimJEEvents = 0;
  void processGenMatched(aod::JetMcCollision const& collision,
                         /*soa::SmallGroups<soa::Join<aod::JMcCollisionLbs, aod::JCollisions>> const& recocolls,*/
                         JetMCDetTable const&, JetMCPartTable const& mcpjets, aod::JetParticles const& mcParticles)
  {

    if (cDebugLevel > 0) {
      nprocessSimJEEvents++;
      if ((nprocessSimJEEvents + 1) % 100000 == 0)
        LOG(debug) << "Jet Events: " << nprocessSimJEEvents;
    }
    if (fabs(collision.posZ()) > 10)
      return;

    jetHist.fill(HIST("genmatched/vertexZ"), collision.posZ());

    bool INELgt0 = false;
    for (const auto& mcParticle : mcParticles) {
      if (fabs(mcParticle.eta()) < cfgtrkMaxEta) {
        INELgt0 = true;
        break;
      }
    }
    if (!INELgt0)
      return;

    std::vector<double> mcdJetPt{};
    std::vector<double> mcdJetPhi{};
    std::vector<double> mcdJetEta{};
    std::vector<double> mcpJetPt{};
    std::vector<double> mcpJetPhi{};
    std::vector<double> mcpJetEta{};

    for (auto& mcpjet : mcpjets) {
      jetHist.fill(HIST("genmatched/hGenJetPt"), mcpjet.pt());
      if (!mcpjet.has_matchedJetGeo())
        continue;
      jetHist.fill(HIST("genmatched/hGenJetPtMatched"), mcpjet.pt());
      for (auto& mcdjet : mcpjet.template matchedJetGeo_as<JetMCDetTable>()) {
        if (!mcdjet.has_matchedJetGeo())
          continue;
        mcdJetPt.push_back(mcdjet.pt());
        mcdJetPhi.push_back(mcdjet.phi());
        mcdJetEta.push_back(mcdjet.eta());
        mcpJetPt.push_back(mcpjet.pt());
        mcpJetPhi.push_back(mcpjet.phi());
        mcpJetEta.push_back(mcpjet.eta());

        jetHist.fill(HIST("genmatched/hRecJetPt"), mcpjet.pt());
        jetHist.fill(HIST("genmatched/hRecJetWithGenPt"), mcdjet.pt());
        jetHist.fill(HIST("genmatched/hRecMatchedJetPt"), mcpjet.pt(), mcpjet.pt() - mcdjet.pt());
        jetHist.fill(HIST("genmatched/hRecMatchedJetPhi"), mcpjet.phi(), mcpjet.phi() - mcdjet.phi());
        jetHist.fill(HIST("genmatched/hRecMatchedJetEta"), mcpjet.eta(), mcpjet.eta() - mcdjet.eta());

      } // mcdJet
    } // mcpJet

    for (const auto& mcParticle : mcParticles) {
      if (fabs(mcParticle.eta()) > cfgtrkMaxEta)
        continue;
      // add pid later

      bool jetFlag = false;
      for (std::size_t iDJet = 0; iDJet < mcpJetPt.size(); iDJet++) {
        double delPhi = TVector2::Phi_mpi_pi(mcpJetPhi[iDJet] - mcParticle.phi());
        double delEta = mcpJetEta[iDJet] - mcParticle.eta();
        double R = TMath::Sqrt((delEta * delEta) + (delPhi * delPhi));

        if (R < cfgjetR) {
          jetFlag = true;
          break;
        }
      } // DetJet
      if (mapPDGToValue(mcParticle.pdgCode()) != 0) {
        jetHist.fill(HIST("genmatched/pt/PtParticleType"), mcParticle.pt(), jetFlag, mapPDGToValue(mcParticle.pdgCode()));
      }
    } // jet constituents

  } // process

  PROCESS_SWITCH(nucleiInJets, processJetTracksData, "nuclei in Jets data", true);
  PROCESS_SWITCH(nucleiInJets, processMCRec, "nuclei in Jets for detectorlevel Jets", true);
  PROCESS_SWITCH(nucleiInJets, processMCGen, "nuclei in Jets MC particlelevel Jets", false);
  PROCESS_SWITCH(nucleiInJets, processRecMatched, "nuclei in Jets rec matched", false);
  PROCESS_SWITCH(nucleiInJets, processGenMatched, "nuclei in Jets gen matched", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<nucleiInJets>(cfgc, TaskName{"nuclei-in-jets"})};
};
