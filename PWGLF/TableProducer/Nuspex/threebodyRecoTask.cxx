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
/// \file threebodyRecoTask.cxx
/// \brief Analysis task for 3-body decay process (now mainly for hypertriton)
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>

#include <cmath>
#include <array>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <TLorentzVector.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "PWGLF/DataModel/Reduced3BodyTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"
#include "CCDB/BasicCCDBManager.h"

#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using ReducedCols = soa::Join<aod::RedCollisions, aod::RedCentFT0Cs>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe>;
using MCLabeledTracksIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

struct Candidate3body {
  // Index
  int mcmotherId;
  int track0Id;
  int track1Id;
  int track2Id;
  // Collision
  float colCentFT0C;
  // sv and candidate
  bool isMatter;
  float invmass;
  float ct;
  float cosPA;
  float dcadaughters;
  float dcacandtopv;
  float vtxradius;
  // daughter tracks
  TLorentzVector lcand;
  TLorentzVector lproton;
  TLorentzVector lpion;
  TLorentzVector lbachelor;
  uint8_t dautpcNclusters[3]; // 0 - proton, 1 - pion, 2 - bachelor
  uint32_t dauitsclussize[3]; // 0 - proton, 1 - pion, 2 - bachelor
  float daudcaxytopv[3];      // 0 - proton, 1 - pion, 2 - bachelor
  float daudcatopv[3];        // 0 - proton, 1 - pion, 2 - bachelor
  float dautpcNsigma[3];      // 0 - proton, 1 - pion, 2 - bachelor
  float dauinnermostR[3];     // 0 - proton, 1 - pion, 2 - bachelor  !!! TracksIU required !!!
  float bachelortofNsigma;
  bool isBachPrimary = false;
  // MC infomartion
  TLorentzVector lgencand = {0, 0, 0, 0};
  float genct = -1;
  float genrapidity = -999;
  bool isSignal = false;
  bool isReco = false;
  int pdgCode = -1;
  bool survivedEventSelection = false;
};

struct ThreebodyRecoTask {

  Produces<aod::Hyp3BodyCands> outputDataTable;
  Produces<aod::MCHyp3BodyCands> outputMCTable;
  std::vector<Candidate3body> candidates3body;
  std::vector<unsigned int> filledMothers;
  std::vector<bool> isGoodCollision;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  //------------------------------------------------------------------
  Preslice<aod::Vtx3BodyDatas> perCollisionVtx3BodyDatas = o2::aod::vtx3body::collisionId;

  // Configuration to enable like-sign analysis
  Configurable<bool> cfgLikeSignAnalysis{"cfgLikeSignAnalysis", false, "Enable like-sign analysis"};
  // Selection criteria
  Configurable<double> vtxcospa{"vtxcospa", 0.99, "Vtx CosPA"};         // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcavtxdau{"dcavtxdau", 1.0, "DCA Vtx Daughters"}; // loose cut
  Configurable<float> dcapiontopv{"dcapiontopv", .05, "DCA Pion To PV"};
  Configurable<float> etacut{"etacut", 0.9, "etacut"};
  Configurable<float> rapiditycut{"rapiditycut", 1, "rapiditycut"};
  Configurable<float> tofPIDNSigmaMin{"tofPIDNSigmaMin", -5, "tofPIDNSigmaMin"};
  Configurable<float> tofPIDNSigmaMax{"tofPIDNSigmaMax", 5, "tofPIDNSigmaMax"};
  Configurable<float> tpcPIDNSigmaCut{"tpcPIDNSigmaCut", 5, "tpcPIDNSigmaCut"};
  Configurable<bool> eventSel8Cut{"eventSel8Cut", true, "flag to enable event sel8 selection"};
  Configurable<bool> mcEventCut{"mcEventCut", true, "flag to enable mc event selection: kIsTriggerTVX and kNoTimeFrameBorder"};
  Configurable<bool> eventPosZCut{"eventPosZCut", true, "flag to enable event posZ selection"};
  Configurable<float> lifetimecut{"lifetimecut", 40., "lifetimecut"}; // ct
  Configurable<float> minProtonPt{"minProtonPt", 0.3, "minProtonPt"};
  Configurable<float> maxProtonPt{"maxProtonPt", 5, "maxProtonPt"};
  Configurable<float> minPionPt{"minPionPt", 0.1, "minPionPt"};
  Configurable<float> maxPionPt{"maxPionPt", 1.2, "maxPionPt"};
  Configurable<float> minDeuteronPt{"minDeuteronPt", 0.6, "minDeuteronPt"};
  Configurable<float> maxDeuteronPt{"maxDeuteronPt", 10, "maxDeuteronPt"};
  Configurable<float> minDeuteronPUseTOF{"minDeuteronPUseTOF", 1, "minDeuteronPt Enable TOF PID"};
  Configurable<float> h3LMassLowerlimit{"h3LMassLowerlimit", 2.96, "Hypertriton mass lower limit"};
  Configurable<float> h3LMassUpperlimit{"h3LMassUpperlimit", 3.04, "Hypertriton mass upper limit"};
  Configurable<int> mintpcNClsproton{"mintpcNClsproton", 90, "min tpc Nclusters for proton"};
  Configurable<int> mintpcNClspion{"mintpcNClspion", 70, "min tpc Nclusters for pion"};
  Configurable<int> mintpcNClsdeuteron{"mintpcNClsdeuteron", 100, "min tpc Nclusters for deuteron"};

  Configurable<float> mcsigma{"mcsigma", 0.0015, "sigma of mc invariant mass fit"}; // obtained from MC
  Configurable<int> bachelorPdgCode{"bachelorPdgCode", 1000010020, "pdgCode of bachelor daughter"};
  Configurable<int> motherPdgCode{"motherPdgCode", 1010010030, "pdgCode of mother track"};

  // 3sigma region for Dalitz plot
  float lowersignallimit = o2::constants::physics::MassHyperTriton - 3 * mcsigma;
  float uppersignallimit = o2::constants::physics::MassHyperTriton + 3 * mcsigma;

  // CCDB options
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "", "Path to the PID response object"};

  // Zorro counting
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};

  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.0f, 5.0f}}}},
      {"hCentFT0C", "hCentFT0C", {HistType::kTH1F, {{100, 0.0f, 100.0f, "FT0C Centrality"}}}},
      {"hCandidatesCounter", "hCandidatesCounter", {HistType::kTH1F, {{12, 0.0f, 12.0f}}}},
      {"hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {{80, 2.96f, 3.04f}}}},
      {"hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {{80, 2.96f, 3.04f}}}},
      {"hMassHypertritonTotal", "hMassHypertritonTotal", {HistType::kTH1F, {{300, 2.9f, 3.2f}}}},
      {"hTOFPIDDeuteron", "hTOFPIDDeuteron", {HistType::kTH1F, {{2000, -100.0f, 100.0f}}}},
      {"hProtonTPCBB", "hProtonTPCBB", {HistType::kTH2F, {{160, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hPionTPCBB", "hPionTPCBB", {HistType::kTH2F, {{160, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hDeuteronTPCBB", "hDeuteronTPCBB", {HistType::kTH2F, {{160, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hProtonTPCVsPt", "hProtonTPCVsPt", {HistType::kTH2F, {{50, 0.0f, 5.0f, "#it{p}_{T} (GeV/c)"}, {240, -6.0f, 6.0f, "TPC n#sigma"}}}},
      {"hPionTPCVsPt", "hPionTPCVsPt", {HistType::kTH2F, {{20, 0.0f, 2.0f, "#it{p}_{T} (GeV/c)"}, {240, -6.0f, 6.0f, "TPC n#sigma"}}}},
      {"hDeuteronTPCVsPt", "hDeuteronTPCVsPt", {HistType::kTH2F, {{80, 0.0f, 8.0f, "#it{p}_{T} (GeV/c)"}, {240, -6.0f, 6.0f, "TPC n#sigma"}}}},
      {"hDeuteronTOFVsPBeforeTOFCut", "hDeuteronTOFVsPBeforeTOFCut", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},
      {"hDeuteronTOFVsPAtferTOFCut", "hDeuteronTOFVsPAtferTOFCut", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},

      {"hDalitz", "hDalitz", {HistType::kTH2F, {{120, 7.85, 8.45, "M^{2}(dp) (GeV^{2}/c^{4})"}, {60, 1.1, 1.4, "M^{2}(p#pi) (GeV^{2}/c^{4})"}}}},
    },
  };

  //------------------------------------------------------------------
  // Fill stats histograms
  enum Vtxstep { kCandAll = 0,
                 kCandDauEta,
                 kCandDauPt,
                 kCandTPCNcls,
                 kCandTPCPID,
                 kCandTOFPID,
                 kCandDcaToPV,
                 kCandRapidity,
                 kCandct,
                 kCandCosPA,
                 kCandDcaDau,
                 kCandInvMass,
                 kNCandSteps };

  struct {
    std::array<int32_t, kNCandSteps> candstats;
    std::array<int32_t, kNCandSteps> truecandstats;
  } statisticsRegistry;

  void resetHistos()
  {
    for (int ii = 0; ii < kNCandSteps; ii++) {
      statisticsRegistry.candstats[ii] = 0;
      statisticsRegistry.truecandstats[ii] = 0;
    }
  }
  void fillCandCounter(int kn, bool istrue = false)
  {
    statisticsRegistry.candstats[kn]++;
    if (istrue) {
      statisticsRegistry.truecandstats[kn]++;
    }
  }
  void fillHistos()
  {
    for (int ii = 0; ii < kNCandSteps; ii++) {
      registry.fill(HIST("hCandidatesCounter"), ii, statisticsRegistry.candstats[ii]);
      if (doprocessMC == true) {
        registry.fill(HIST("hTrueHypertritonCounter"), ii, statisticsRegistry.truecandstats[ii]);
      }
    }
  }

  int mRunNumber;

  void init(InitContext const&)
  {
    zorroSummary.setObject(zorro.getZorroSummary());
    mRunNumber = 0;

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(1, "total");
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(2, "sel8");
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(3, "vertexZ");
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(4, "Zorro H3L 3body event");
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(5, "has Candidate");

    // Check for selection criteria  !!! TracksIU required !!!
    registry.add("hDiffRVtxProton", "hDiffRVtxProton", HistType::kTH1F, {{100, -10, 10}});     // difference between the radius of decay vertex and minR of proton
    registry.add("hDiffRVtxPion", "hDiffRVtxPion", HistType::kTH1F, {{100, -10, 10}});         // difference between the radius of decay vertex and minR of pion
    registry.add("hDiffRVtxDeuteron", "hDiffRVtxDeuteron", HistType::kTH1F, {{100, -10, 10}}); // difference between the radius of decay vertex and minR of deuteron
    registry.add("hDiffDaughterR", "hDiffDaughterR", HistType::kTH1F, {{10000, -100, 100}});   // difference between minR of pion&proton and R of deuteron(bachelor)

    if (cfgLikeSignAnalysis) {
      registry.add("hInvMassCorrectSign", "hInvMassCorrectSign", HistType::kTH1F, {{80, 2.96f, 3.04f}}); // check if there are contamination of possible signals which are caused by unexpected PID
    }

    if (doprocessMC == true) {
      registry.add("hTrueHypertritonCounter", "hTrueHypertritonCounter", HistType::kTH1F, {{12, 0.0f, 12.0f}});
      auto hGeneratedHypertritonCounter = registry.add<TH1>("hGeneratedHypertritonCounter", "hGeneratedHypertritonCounter", HistType::kTH1F, {{2, 0.0f, 2.0f}});
      hGeneratedHypertritonCounter->GetXaxis()->SetBinLabel(1, "Total");
      hGeneratedHypertritonCounter->GetXaxis()->SetBinLabel(2, "3-body decay");
      registry.add("hPtGeneratedHypertriton", "hPtGeneratedHypertriton", HistType::kTH1F, {{200, 0.0f, 10.0f}});
      registry.add("hctGeneratedHypertriton", "hctGeneratedHypertriton", HistType::kTH1F, {{50, 0, 50, "ct(cm)"}});
      registry.add("hEtaGeneratedHypertriton", "hEtaGeneratedHypertriton", HistType::kTH1F, {{40, -2.0f, 2.0f}});
      registry.add("hRapidityGeneratedHypertriton", "hRapidityGeneratedHypertriton", HistType::kTH1F, {{40, -2.0f, 2.0f}});
      registry.add("hPtGeneratedAntiHypertriton", "hPtGeneratedAntiHypertriton", HistType::kTH1F, {{200, 0.0f, 10.0f}});
      registry.add("hctGeneratedAntiHypertriton", "hctGeneratedAntiHypertriton", HistType::kTH1F, {{50, 0, 50, "ct(cm)"}});
      registry.add("hEtaGeneratedAntiHypertriton", "hEtaGeneratedAntiHypertriton", HistType::kTH1F, {{40, -2.0f, 2.0f}});
      registry.add("hRapidityGeneratedAntiHypertriton", "hRapidityGeneratedAntiHypertriton", HistType::kTH1F, {{40, -2.0f, 2.0f}});
    }

    TString candCounterbinLabel[kNCandSteps] = {"Total", "TrackEta", "DauPt", "TPCNcls", "TPCPID", "d TOFPID", "PionDcatoPV", "MomRapidity", "Lifetime", "VtxCosPA", "VtxDcaDau", "InvMass"};
    for (int i{0}; i < kNCandSteps; i++) {
      registry.get<TH1>(HIST("hCandidatesCounter"))->GetXaxis()->SetBinLabel(i + 1, candCounterbinLabel[i]);
      if (doprocessMC == true) {
        registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(i + 1, candCounterbinLabel[i]);
      }
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fH3L3Body");
      zorro.populateHistRegistry(registry, bc.runNumber());
    }

    mRunNumber = bc.runNumber();
  }

  //------------------------------------------------------------------
  // Check if the mcparticle is hypertriton which decays into 3 daughters
  template <class TMCTrackTo, typename TMCParticle>
  bool is3bodyDecayed(TMCParticle const& particle)
  {
    if (std::abs(particle.pdgCode()) != motherPdgCode) {
      return false;
    }
    bool haveProton = false, havePion = false, haveBachelor = false;
    bool haveAntiProton = false, haveAntiPion = false, haveAntiBachelor = false;
    for (const auto& mcparticleDaughter : particle.template daughters_as<TMCTrackTo>()) {
      if (mcparticleDaughter.pdgCode() == 2212)
        haveProton = true;
      if (mcparticleDaughter.pdgCode() == -2212)
        haveAntiProton = true;
      if (mcparticleDaughter.pdgCode() == 211)
        havePion = true;
      if (mcparticleDaughter.pdgCode() == -211)
        haveAntiPion = true;
      if (mcparticleDaughter.pdgCode() == bachelorPdgCode)
        haveBachelor = true;
      if (mcparticleDaughter.pdgCode() == -bachelorPdgCode)
        haveAntiBachelor = true;
    }
    if (haveProton && haveAntiPion && haveBachelor && particle.pdgCode() > 0) {
      return true;
    } else if (haveAntiProton && havePion && haveAntiBachelor && particle.pdgCode() < 0) {
      return true;
    }
    return false;
  }

  //------------------------------------------------------------------
  // Event Selection
  template <typename TCollision>
  bool eventSelection(TCollision const& collision)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    registry.fill(HIST("hEventCounter"), 0.5);

    if (eventSel8Cut && !collision.sel8()) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), 1.5);
    if (eventPosZCut && std::abs(collision.posZ()) > 10.f) { // 10cm
      return false;
    }
    registry.fill(HIST("hEventCounter"), 2.5);
    registry.fill(HIST("hCentFT0C"), collision.centFT0C());

    if (cfgSkimmedProcessing) {
      bool zorroSelected = zorro.isSelected(bc.globalBC()); /// Just let Zorro do the accounting
      if (zorroSelected) {
        registry.fill(HIST("hEventCounter"), 3.5);
      }
    }

    return true;
  }

  //------------------------------------------------------------------
  // Fill candidate table
  template <typename TCollisionTable, typename TTrackTable, typename TCandTable>
  void fillCand(TCollisionTable const& collision, TCandTable const& candData, TTrackTable const& trackProton, TTrackTable const& trackPion, TTrackTable const& trackDeuteron, bool isMatter, bool isTrueCand = false, int lLabel = -1, TLorentzVector lmother = {0, 0, 0, 0}, double mcLifetime = -1, bool isBachPrimary = false)
  {
    double cospa = candData.vtxcosPA(collision.posX(), collision.posY(), collision.posZ());
    double ct = candData.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassHyperTriton;

    Candidate3body cand3body;
    cand3body.isMatter = isMatter;
    if (isMatter == true) {
      cand3body.lproton.SetXYZM(candData.pxtrack0(), candData.pytrack0(), candData.pztrack0(), o2::constants::physics::MassProton);
      cand3body.lpion.SetXYZM(candData.pxtrack1(), candData.pytrack1(), candData.pztrack1(), o2::constants::physics::MassPionCharged);
    } else {
      cand3body.lproton.SetXYZM(candData.pxtrack1(), candData.pytrack1(), candData.pztrack1(), o2::constants::physics::MassPionCharged);
      cand3body.lpion.SetXYZM(candData.pxtrack0(), candData.pytrack0(), candData.pztrack0(), o2::constants::physics::MassProton);
    }

    cand3body.mcmotherId = lLabel;
    cand3body.track0Id = candData.track0Id();
    cand3body.track1Id = candData.track1Id();
    cand3body.track2Id = candData.track2Id();
    cand3body.colCentFT0C = collision.centFT0C();
    cand3body.invmass = cand3body.isMatter ? candData.mHypertriton() : candData.mAntiHypertriton();
    cand3body.lcand.SetXYZM(candData.px(), candData.py(), candData.pz(), o2::constants::physics::MassHyperTriton);
    cand3body.ct = ct;
    cand3body.cosPA = cospa;
    cand3body.dcadaughters = candData.dcaVtxdaughters();
    cand3body.dcacandtopv = candData.dcavtxtopv(collision.posX(), collision.posY(), collision.posZ());
    cand3body.vtxradius = candData.vtxradius();
    cand3body.lbachelor.SetXYZM(candData.pxtrack2(), candData.pytrack2(), candData.pztrack2(), o2::constants::physics::MassDeuteron);
    cand3body.dautpcNclusters[0] = trackProton.tpcNClsFound();
    cand3body.dautpcNclusters[1] = trackPion.tpcNClsFound();
    cand3body.dautpcNclusters[2] = trackDeuteron.tpcNClsFound();
    cand3body.dauitsclussize[0] = trackProton.itsClusterSizes();
    cand3body.dauitsclussize[1] = trackPion.itsClusterSizes();
    cand3body.dauitsclussize[2] = trackDeuteron.itsClusterSizes();
    cand3body.dautpcNsigma[0] = trackProton.tpcNSigmaPr();
    cand3body.dautpcNsigma[1] = trackPion.tpcNSigmaPi();
    cand3body.dautpcNsigma[2] = trackDeuteron.tpcNSigmaDe();
    cand3body.daudcaxytopv[0] = cand3body.isMatter ? candData.dcaXYtrack0topv() : candData.dcaXYtrack1topv();
    cand3body.daudcaxytopv[1] = cand3body.isMatter ? candData.dcaXYtrack1topv() : candData.dcaXYtrack0topv();
    cand3body.daudcaxytopv[2] = candData.dcaXYtrack2topv();
    cand3body.daudcatopv[0] = cand3body.isMatter ? candData.dcatrack0topv() : candData.dcatrack1topv();
    cand3body.daudcatopv[1] = cand3body.isMatter ? candData.dcatrack1topv() : candData.dcatrack0topv();
    cand3body.daudcatopv[2] = candData.dcatrack2topv();
    cand3body.dauinnermostR[0] = trackProton.x();
    cand3body.dauinnermostR[1] = trackPion.x();
    cand3body.dauinnermostR[2] = trackDeuteron.x();

    cand3body.bachelortofNsigma = candData.tofNSigmaBachDe();
    cand3body.isBachPrimary = isBachPrimary;
    if (isTrueCand) {
      cand3body.mcmotherId = lLabel;
      cand3body.lgencand = lmother;
      cand3body.genct = mcLifetime;
      cand3body.genrapidity = lmother.Rapidity();
      cand3body.isSignal = true;
      cand3body.isReco = true;
      cand3body.pdgCode = cand3body.isMatter ? motherPdgCode : -motherPdgCode;
      cand3body.survivedEventSelection = true;
      filledMothers.push_back(lLabel);
    }

    candidates3body.push_back(cand3body);

    registry.fill(HIST("hProtonTPCBB"), trackProton.sign() * trackProton.p(), trackProton.tpcSignal());
    registry.fill(HIST("hPionTPCBB"), trackPion.sign() * trackPion.p(), trackPion.tpcSignal());
    registry.fill(HIST("hDeuteronTPCBB"), trackDeuteron.sign() * trackDeuteron.p(), trackDeuteron.tpcSignal());
    registry.fill(HIST("hProtonTPCVsPt"), trackProton.pt(), trackProton.tpcNSigmaPr());
    registry.fill(HIST("hPionTPCVsPt"), trackProton.pt(), trackPion.tpcNSigmaPi());
    registry.fill(HIST("hDeuteronTPCVsPt"), trackDeuteron.pt(), trackDeuteron.tpcNSigmaDe());
    registry.fill(HIST("hTOFPIDDeuteron"), candData.tofNSigmaBachDe());
    registry.fill(HIST("hDiffRVtxProton"), trackProton.x() - candData.vtxradius());
    registry.fill(HIST("hDiffRVtxPion"), trackPion.x() - candData.vtxradius());
    registry.fill(HIST("hDiffRVtxDeuteron"), trackDeuteron.x() - candData.vtxradius());
    float diffTrackR = trackDeuteron.x() - std::min(trackProton.x(), trackPion.x());
    registry.fill(HIST("hDiffDaughterR"), diffTrackR);
  }

  //------------------------------------------------------------------
  // Selections for candidates
  template <typename TCollisionTable, typename TTrackTable, typename TCandTable>
  bool selectCand(TCollisionTable const& collision, TCandTable const& candData, TTrackTable const& trackProton, TTrackTable const& trackPion, TTrackTable const& trackDeuteron, bool isMatter, bool isTrueCand = false)
  {
    fillCandCounter(kCandAll, isTrueCand);

    // Selection on daughters
    if (std::abs(trackProton.eta()) > etacut || std::abs(trackPion.eta()) > etacut || std::abs(trackDeuteron.eta()) > etacut) {
      return false;
    }
    fillCandCounter(kCandDauEta, isTrueCand);

    if (trackProton.pt() < minProtonPt || trackProton.pt() > maxProtonPt || trackPion.pt() < minPionPt || trackPion.pt() > maxPionPt || trackDeuteron.pt() < minDeuteronPt || trackDeuteron.pt() > maxDeuteronPt) {
      return false;
    }
    fillCandCounter(kCandDauPt, isTrueCand);

    if (trackProton.tpcNClsFound() < mintpcNClsproton || trackPion.tpcNClsFound() < mintpcNClspion || trackDeuteron.tpcNClsFound() < mintpcNClsdeuteron) {
      return false;
    }
    fillCandCounter(kCandTPCNcls, isTrueCand);

    if (std::abs(trackProton.tpcNSigmaPr()) > tpcPIDNSigmaCut || std::abs(trackPion.tpcNSigmaPi()) > tpcPIDNSigmaCut || std::abs(trackDeuteron.tpcNSigmaDe()) > tpcPIDNSigmaCut) {
      return false;
    }
    fillCandCounter(kCandTPCPID, isTrueCand);

    registry.fill(HIST("hDeuteronTOFVsPBeforeTOFCut"), trackDeuteron.sign() * trackDeuteron.p(), candData.tofNSigmaBachDe());
    if ((candData.tofNSigmaBachDe() < tofPIDNSigmaMin || candData.tofNSigmaBachDe() > tofPIDNSigmaMax) && trackDeuteron.p() > minDeuteronPUseTOF) {
      return false;
    }
    fillCandCounter(kCandTOFPID, isTrueCand);
    registry.fill(HIST("hDeuteronTOFVsPAtferTOFCut"), trackDeuteron.sign() * trackDeuteron.p(), candData.tofNSigmaBachDe());

    double dcapion = isMatter ? candData.dcatrack1topv() : candData.dcatrack0topv();
    if (std::abs(dcapion) < dcapiontopv) {
      return false;
    }
    fillCandCounter(kCandDcaToPV, isTrueCand);

    // Selection on candidate hypertriton
    if (std::abs(candData.yHypertriton()) > rapiditycut) {
      return false;
    }
    fillCandCounter(kCandRapidity, isTrueCand);

    double ct = candData.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassHyperTriton;
    if (ct > lifetimecut) {
      return false;
    }
    fillCandCounter(kCandct, isTrueCand);

    double cospa = candData.vtxcosPA(collision.posX(), collision.posY(), collision.posZ());
    if (cospa < vtxcospa) {
      return false;
    }
    fillCandCounter(kCandCosPA, isTrueCand);

    if (candData.dcaVtxdaughters() > dcavtxdau) {
      return false;
    }
    fillCandCounter(kCandDcaDau, isTrueCand);

    if ((isMatter && candData.mHypertriton() > h3LMassLowerlimit && candData.mHypertriton() < h3LMassUpperlimit)) {
      // Hypertriton
      registry.fill(HIST("hMassHypertriton"), candData.mHypertriton());
      registry.fill(HIST("hMassHypertritonTotal"), candData.mHypertriton());
      if (candData.mHypertriton() > lowersignallimit && candData.mHypertriton() < uppersignallimit) {
        registry.fill(HIST("hDalitz"), RecoDecay::m2(std::array{std::array{candData.pxtrack0(), candData.pytrack0(), candData.pztrack0()}, std::array{candData.pxtrack2(), candData.pytrack2(), candData.pztrack2()}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}), RecoDecay::m2(std::array{std::array{candData.pxtrack0(), candData.pytrack0(), candData.pztrack0()}, std::array{candData.pxtrack1(), candData.pytrack1(), candData.pztrack1()}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}));
      }
    } else if ((!isMatter && candData.mAntiHypertriton() > h3LMassLowerlimit && candData.mAntiHypertriton() < h3LMassUpperlimit)) {
      // AntiHypertriton
      registry.fill(HIST("hMassAntiHypertriton"), candData.mAntiHypertriton());
      registry.fill(HIST("hMassHypertritonTotal"), candData.mAntiHypertriton());
      if (candData.mAntiHypertriton() > lowersignallimit && candData.mAntiHypertriton() < uppersignallimit) {
        registry.fill(HIST("hDalitz"), RecoDecay::m2(std::array{std::array{candData.pxtrack1(), candData.pytrack1(), candData.pztrack1()}, std::array{candData.pxtrack2(), candData.pytrack2(), candData.pztrack2()}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}), RecoDecay::m2(std::array{std::array{candData.pxtrack1(), candData.pytrack1(), candData.pztrack1()}, std::array{candData.pxtrack0(), candData.pytrack0(), candData.pztrack0()}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}));
      }
    } else {
      return false;
    }
    fillCandCounter(kCandInvMass, isTrueCand);

    return true;
  }

  //------------------------------------------------------------------
  // Analysis process for a single candidate
  template <class TTrackClass, typename TCollisionTable, typename TCandTable>
  void candidateAnalysis(TCollisionTable const& collision, TCandTable const& candData, bool& ifHasCandidate, bool isTrueCand = false, int lLabel = -1, TLorentzVector lmother = {0, 0, 0, 0}, double mcLifetime = -1, double isBachPrimary = false)
  {
    auto track0 = candData.template track0_as<TTrackClass>();
    auto track1 = candData.template track1_as<TTrackClass>();
    auto track2 = candData.template track2_as<TTrackClass>();

    bool isMatter = track2.sign() > 0; // true if the candidate is hypertriton (p pi- d)

    auto& trackProton = isMatter ? track0 : track1;
    auto& trackPion = isMatter ? track1 : track0;
    auto& trackDeuteron = track2;

    if (selectCand(collision, candData, trackProton, trackPion, trackDeuteron, isMatter, isTrueCand)) {
      ifHasCandidate = true;
      fillCand(collision, candData, trackProton, trackPion, trackDeuteron, isMatter, isTrueCand, lLabel, lmother, mcLifetime, isBachPrimary);
    }
  }

  //------------------------------------------------------------------
  // Analysis process for like-sign candidates : (p pi- anti-d) or (anti-p pi+ d)
  template <class TTrackClass, typename TCollisionTable, typename TCandTable>
  void likeSignAnalysis(TCollisionTable const& collision, TCandTable const& candData, bool& ifHasCandidate, bool isTrueCand = false, int lLabel = -1, TLorentzVector lmother = {0, 0, 0, 0}, double mcLifetime = -1, double isBachPrimary = false)
  {
    auto track0 = candData.template track0_as<TTrackClass>();
    auto track1 = candData.template track1_as<TTrackClass>();
    auto track2 = candData.template track2_as<TTrackClass>();

    bool isMatter = track2.sign() < 0; // true if seach for background consists of (p pi- anti-d)

    // Assume proton has an oppisite charge with deuteron
    auto& trackProton = isMatter ? track0 : track1;
    auto& trackPion = isMatter ? track1 : track0;
    auto& trackDeuteron = track2;

    if (selectCand(collision, candData, trackProton, trackPion, trackDeuteron, isMatter, isTrueCand)) {
      ifHasCandidate = true;
      fillCand(collision, candData, trackProton, trackPion, trackDeuteron, isMatter, isTrueCand, lLabel, lmother, mcLifetime, isBachPrimary);
      // QA for if signals have the possibility to be reconginzed as a like-sign background
      if (isMatter) {
        registry.fill(HIST("hInvMassCorrectSign"), candData.mHypertriton());
      } else {
        registry.fill(HIST("hInvMassCorrectSign"), candData.mAntiHypertriton());
      }
    }
  }

  //------------------------------------------------------------------
  // Analysis process for reduced data
  template <class TTrackClass, typename TCollisionTable, typename TCandTable, typename TTracks>
  void reducedAnalysis(TCollisionTable const& collision, TCandTable const& candData, TTracks tracks, bool isTrueCand = false, int lLabel = -1, TLorentzVector lmother = {0, 0, 0, 0}, double mcLifetime = -1, double isBachPrimary = false)
  {
    auto track0 = tracks.rawIteratorAt(candData.track0Id());
    auto track1 = tracks.rawIteratorAt(candData.track1Id());
    auto track2 = tracks.rawIteratorAt(candData.track2Id());

    bool isMatter = track2.sign() > 0; // true if the candidate is hypertriton (p pi- d)

    auto& trackProton = isMatter ? track0 : track1;
    auto& trackPion = isMatter ? track1 : track0;
    auto& trackDeuteron = track2;

    if (selectCand(collision, candData, trackProton, trackPion, trackDeuteron, isMatter, isTrueCand)) {
      fillCand(collision, candData, trackProton, trackPion, trackDeuteron, isMatter, isTrueCand, lLabel, lmother, mcLifetime, isBachPrimary);
    }
  }

  //------------------------------------------------------------------
  // Analysis process for reduced data with like-sign candidates
  template <class TTrackClass, typename TCollisionTable, typename TCandTable, typename TTracks>
  void reducedLikeSignAnalysis(TCollisionTable const& collision, TCandTable const& candData, TTracks tracks, bool isTrueCand = false, int lLabel = -1, TLorentzVector lmother = {0, 0, 0, 0}, double mcLifetime = -1, bool isBachPrimary = false)
  {
    auto track0 = tracks.rawIteratorAt(candData.track0Id());
    auto track1 = tracks.rawIteratorAt(candData.track1Id());
    auto track2 = tracks.rawIteratorAt(candData.track2Id());

    bool isMatter = track2.sign() < 0; // true if seach for background consists of (p pi- anti-d)

    // Assume proton has an oppisite charge with deuteron
    auto& trackProton = isMatter ? track0 : track1;
    auto& trackPion = isMatter ? track1 : track0;
    auto& trackDeuteron = track2;

    if (selectCand(collision, candData, trackProton, trackPion, trackDeuteron, isMatter, isTrueCand)) {
      fillCand(collision, candData, trackProton, trackPion, trackDeuteron, isMatter, isTrueCand, lLabel, lmother, mcLifetime, isBachPrimary);
      // QA for if signals have the possibility to be reconginzed as a like-sign background
      if (isMatter) {
        registry.fill(HIST("hInvMassCorrectSign"), candData.mHypertriton());
      } else {
        registry.fill(HIST("hInvMassCorrectSign"), candData.mAntiHypertriton());
      }
    }
  }

  //------------------------------------------------------------------
  void fillOutputDataTable(Candidate3body const& cand3body)
  {
    outputDataTable(cand3body.colCentFT0C,
                    cand3body.isMatter, cand3body.invmass, cand3body.lcand.P(), cand3body.lcand.Pt(), cand3body.ct,
                    cand3body.cosPA, cand3body.dcadaughters, cand3body.dcacandtopv, cand3body.vtxradius,
                    cand3body.lproton.Pt(), cand3body.lproton.Eta(), cand3body.lproton.Phi(), cand3body.dauinnermostR[0],
                    cand3body.lpion.Pt(), cand3body.lpion.Eta(), cand3body.lpion.Phi(), cand3body.dauinnermostR[1],
                    cand3body.lbachelor.Pt(), cand3body.lbachelor.Eta(), cand3body.lbachelor.Phi(), cand3body.dauinnermostR[2],
                    cand3body.dautpcNclusters[0], cand3body.dautpcNclusters[1], cand3body.dautpcNclusters[2],
                    cand3body.dauitsclussize[0], cand3body.dauitsclussize[1], cand3body.dauitsclussize[2],
                    cand3body.dautpcNsigma[0], cand3body.dautpcNsigma[1], cand3body.dautpcNsigma[2], cand3body.bachelortofNsigma,
                    cand3body.daudcaxytopv[0], cand3body.daudcaxytopv[1], cand3body.daudcaxytopv[2],
                    cand3body.daudcatopv[0], cand3body.daudcatopv[1], cand3body.daudcatopv[2]);
  }

  //------------------------------------------------------------------
  // collect information for generated hypertriton (should be called after event selection)
  void getGeneratedH3LInfo(aod::McParticles const& particlesMC)
  {
    for (const auto& mcparticle : particlesMC) {
      if (std::abs(mcparticle.pdgCode()) != motherPdgCode) {
        continue;
      }
      registry.fill(HIST("hGeneratedHypertritonCounter"), 0.5);

      bool haveProton = false, havePionPlus = false, haveDeuteron = false;
      bool haveAntiProton = false, havePionMinus = false, haveAntiDeuteron = false;
      double mcLifetime = -1;
      for (const auto& mcparticleDaughter : mcparticle.template daughters_as<aod::McParticles>()) {
        if (mcparticleDaughter.pdgCode() == 2212)
          haveProton = true;
        if (mcparticleDaughter.pdgCode() == -2212)
          haveAntiProton = true;
        if (mcparticleDaughter.pdgCode() == 211)
          havePionPlus = true;
        if (mcparticleDaughter.pdgCode() == -211)
          havePionMinus = true;
        if (mcparticleDaughter.pdgCode() == bachelorPdgCode) {
          haveDeuteron = true;
          mcLifetime = RecoDecay::sqrtSumOfSquares(mcparticleDaughter.vx() - mcparticle.vx(), mcparticleDaughter.vy() - mcparticle.vy(), mcparticleDaughter.vz() - mcparticle.vz()) * o2::constants::physics::MassHyperTriton / mcparticle.p();
        }
        if (mcparticleDaughter.pdgCode() == -bachelorPdgCode) {
          haveAntiDeuteron = true;
          mcLifetime = RecoDecay::sqrtSumOfSquares(mcparticleDaughter.vx() - mcparticle.vx(), mcparticleDaughter.vy() - mcparticle.vy(), mcparticleDaughter.vz() - mcparticle.vz()) * o2::constants::physics::MassHyperTriton / mcparticle.p();
        }
      }
      if (haveProton && havePionMinus && haveDeuteron && mcparticle.pdgCode() == motherPdgCode) {
        registry.fill(HIST("hGeneratedHypertritonCounter"), 1.5);
        registry.fill(HIST("hPtGeneratedHypertriton"), mcparticle.pt());
        registry.fill(HIST("hctGeneratedHypertriton"), mcLifetime);
        registry.fill(HIST("hEtaGeneratedHypertriton"), mcparticle.eta());
        registry.fill(HIST("hRapidityGeneratedHypertriton"), mcparticle.y());
      } else if (haveAntiProton && havePionPlus && haveAntiDeuteron && mcparticle.pdgCode() == -motherPdgCode) {
        registry.fill(HIST("hGeneratedHypertritonCounter"), 1.5);
        registry.fill(HIST("hPtGeneratedAntiHypertriton"), mcparticle.pt());
        registry.fill(HIST("hctGeneratedAntiHypertriton"), mcLifetime);
        registry.fill(HIST("hEtaGeneratedAntiHypertriton"), mcparticle.eta());
        registry.fill(HIST("hRapidityGeneratedAntiHypertriton"), mcparticle.y());
      }
    }
  }

  //------------------------------------------------------------------
  // process real data analysis
  void processData(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions, aod::Vtx3BodyDatas const& vtx3bodydatas, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      candidates3body.clear();

      if (!eventSelection(collision)) {
        continue;
      }

      bool ifHasCandidate = false;
      auto d3BodyCands = vtx3bodydatas.sliceBy(perCollisionVtx3BodyDatas, collision.globalIndex());
      for (const auto& vtx : d3BodyCands) {
        if (cfgLikeSignAnalysis) {
          likeSignAnalysis<FullTracksExtIU>(collision, vtx, ifHasCandidate);
        } else {
          candidateAnalysis<FullTracksExtIU>(collision, vtx, ifHasCandidate);
        }
      }
      if (ifHasCandidate)
        registry.fill(HIST("hEventCounter"), 4.5);
      fillHistos();
      resetHistos();

      for (const auto& cand3body : candidates3body) {
        fillOutputDataTable(cand3body);
      }
    }
  }
  PROCESS_SWITCH(ThreebodyRecoTask, processData, "Real data reconstruction", true);

  //------------------------------------------------------------------
  // process reduced data analysis
  void processReducedData(ReducedCols const& collisions, aod::Vtx3BodyDatas const& vtx3bodydatas, aod::RedIUTracks const& tracks)
  {
    candidates3body.clear();

    for (const auto& vtx : vtx3bodydatas) {
      const auto& collision = collisions.iteratorAt(vtx.collisionId());
      if (cfgLikeSignAnalysis) {
        reducedLikeSignAnalysis<aod::RedIUTracks>(collision, vtx, tracks);
      } else {
        reducedAnalysis<aod::RedIUTracks>(collision, vtx, tracks);
      }
      for (const auto& cand3body : candidates3body) {
        fillOutputDataTable(cand3body);
      }
      candidates3body.clear();
    }
    fillHistos();
    resetHistos();
  }
  PROCESS_SWITCH(ThreebodyRecoTask, processReducedData, "Reduced data reconstruction", false);

  //------------------------------------------------------------------
  // process mc analysis
  void processMC(soa::Join<aod::Collisions, o2::aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs> const& collisions, aod::Vtx3BodyDatas const& vtx3bodydatas, aod::McParticles const& particlesMC, MCLabeledTracksIU const& /*tracks*/, aod::McCollisions const& mcCollisions)
  {
    filledMothers.clear();
    getGeneratedH3LInfo(particlesMC);
    isGoodCollision.resize(mcCollisions.size(), false);

    for (const auto& collision : collisions) {
      candidates3body.clear();
      registry.fill(HIST("hEventCounter"), 0.5);
      if (mcEventCut && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1.5);
      if (eventPosZCut && std::abs(collision.posZ()) > 10.f) { // 10cm
        continue;
      }
      registry.fill(HIST("hEventCounter"), 2.5);
      registry.fill(HIST("hCentFT0C"), collision.centFT0C());

      if (collision.mcCollisionId() >= 0) {
        isGoodCollision[collision.mcCollisionId()] = true;
      }

      bool ifHasCandidate = false;
      auto vtxsThisCol = vtx3bodydatas.sliceBy(perCollisionVtx3BodyDatas, collision.globalIndex());

      for (const auto& vtx : vtxsThisCol) {
        bool isBachPrimary = false;
        int lLabel = -1;
        int lPDG = -1;
        double mcLifetime = -1;
        TLorentzVector lmother;
        bool isTrueCand = false;
        auto track0 = vtx.track0_as<MCLabeledTracksIU>();
        auto track1 = vtx.track1_as<MCLabeledTracksIU>();
        auto track2 = vtx.track2_as<MCLabeledTracksIU>();
        if (track0.has_mcParticle() && track1.has_mcParticle() && track2.has_mcParticle()) {
          auto lMCTrack0 = track0.mcParticle_as<aod::McParticles>();
          auto lMCTrack1 = track1.mcParticle_as<aod::McParticles>();
          auto lMCTrack2 = track2.mcParticle_as<aod::McParticles>();

          if (lMCTrack2.isPhysicalPrimary()) {
            isBachPrimary = true;
          }

          if (lMCTrack0.has_mothers() && lMCTrack1.has_mothers() && lMCTrack2.has_mothers()) {
            for (const auto& lMother0 : lMCTrack0.mothers_as<aod::McParticles>()) {
              for (const auto& lMother1 : lMCTrack1.mothers_as<aod::McParticles>()) {
                for (const auto& lMother2 : lMCTrack2.mothers_as<aod::McParticles>()) {
                  if (lMother0.globalIndex() == lMother1.globalIndex() && lMother0.globalIndex() == lMother2.globalIndex()) {
                    lLabel = lMother0.globalIndex();
                    lPDG = lMother0.pdgCode();
                    if ((lPDG == motherPdgCode && lMCTrack0.pdgCode() == 2212 && lMCTrack1.pdgCode() == -211 && lMCTrack2.pdgCode() == bachelorPdgCode) ||
                        (lPDG == -motherPdgCode && lMCTrack0.pdgCode() == 211 && lMCTrack1.pdgCode() == -2212 && lMCTrack2.pdgCode() == -bachelorPdgCode)) {
                      isTrueCand = true;
                      mcLifetime = RecoDecay::sqrtSumOfSquares(lMCTrack2.vx() - lMother2.vx(), lMCTrack2.vy() - lMother2.vy(), lMCTrack2.vz() - lMother2.vz()) * o2::constants::physics::MassHyperTriton / lMother2.p();
                      lmother.SetXYZM(lMother0.px(), lMother0.py(), lMother0.pz(), o2::constants::physics::MassHyperTriton);
                    }
                  }
                }
              }
            }
          }
        }

        candidateAnalysis<MCLabeledTracksIU>(collision, vtx, ifHasCandidate, isTrueCand, lLabel, lmother, mcLifetime, isBachPrimary);
      }

      if (ifHasCandidate)
        registry.fill(HIST("hEventCounter"), 4.5);
      fillHistos();
      resetHistos();

      for (const auto& cand3body : candidates3body) {
        outputMCTable(cand3body.colCentFT0C,
                      cand3body.isMatter, cand3body.invmass, cand3body.lcand.P(), cand3body.lcand.Pt(), cand3body.ct,
                      cand3body.cosPA, cand3body.dcadaughters, cand3body.dcacandtopv, cand3body.vtxradius,
                      cand3body.lproton.Pt(), cand3body.lproton.Eta(), cand3body.lproton.Phi(), cand3body.dauinnermostR[0],
                      cand3body.lpion.Pt(), cand3body.lpion.Eta(), cand3body.lpion.Phi(), cand3body.dauinnermostR[1],
                      cand3body.lbachelor.Pt(), cand3body.lbachelor.Eta(), cand3body.lbachelor.Phi(), cand3body.dauinnermostR[2],
                      cand3body.dautpcNclusters[0], cand3body.dautpcNclusters[1], cand3body.dautpcNclusters[2],
                      cand3body.dauitsclussize[0], cand3body.dauitsclussize[1], cand3body.dauitsclussize[2],
                      cand3body.dautpcNsigma[0], cand3body.dautpcNsigma[1], cand3body.dautpcNsigma[2], cand3body.bachelortofNsigma,
                      cand3body.daudcaxytopv[0], cand3body.daudcaxytopv[1], cand3body.daudcaxytopv[2],
                      cand3body.daudcatopv[0], cand3body.daudcatopv[1], cand3body.daudcatopv[2],
                      cand3body.isBachPrimary,
                      cand3body.lgencand.P(), cand3body.lgencand.Pt(), cand3body.genct, cand3body.lgencand.Phi(), cand3body.lgencand.Eta(), cand3body.lgencand.Rapidity(),
                      cand3body.isSignal, cand3body.isReco, cand3body.pdgCode, cand3body.survivedEventSelection);
      }
    }

    // now we fill only the signal candidates that were not reconstructed
    for (const auto& mcparticle : particlesMC) {
      if (!is3bodyDecayed<aod::McParticles>(mcparticle)) {
        continue;
      }
      if (std::find(filledMothers.begin(), filledMothers.end(), mcparticle.globalIndex()) != std::end(filledMothers)) {
        continue;
      }
      bool isSurEvSelection = isGoodCollision[mcparticle.mcCollisionId()];
      std::array<float, 3> posSV{0.f};
      for (const auto& mcDaughter : mcparticle.daughters_as<aod::McParticles>()) {
        if (std::abs(mcDaughter.pdgCode()) == bachelorPdgCode) {
          posSV = {mcDaughter.vx(), mcDaughter.vy(), mcDaughter.vz()};
        }
      }
      double mcLifetime = RecoDecay::sqrtSumOfSquares(posSV[0] - mcparticle.vx(), posSV[1] - mcparticle.vy(), posSV[2] - mcparticle.vz()) * o2::constants::physics::MassHyperTriton / mcparticle.p();
      outputMCTable(-1,
                    -1, -1, -1, -1, -1,
                    -1, -1, -1, -1,
                    -1, -1, -1, -1,
                    -1, -1, -1, -1,
                    -1, -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    false,
                    mcparticle.p(), mcparticle.pt(), mcLifetime, mcparticle.phi(), mcparticle.eta(), mcparticle.y(),
                    true, false, mcparticle.pdgCode(), isSurEvSelection);
    }
  }
  PROCESS_SWITCH(ThreebodyRecoTask, processMC, "MC reconstruction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ThreebodyRecoTask>(cfgc),
  };
}
