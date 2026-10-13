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

/// \file upcTrackVertexingQA.cxx
/// \brief task to study the performance of vertexing for low-multiplicity UPC collisions
/// \author Andrea Tavira Garcia a.tavira@cern.ch
/// \author Andrea Giovanni Riffero andrea.giovanni.riffero@cern.ch

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>
#include <TH1.h>
#include <TH2.h>
#include <TPDGCode.h>

#include <array>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum class CandSpecies {
  kRho = 0,
  kJpsi = 1
};

struct UpcTrackVertexingQA {
  // Configurables
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 0., "min. cand. pT (GeV/c)"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track pT (GeV/c)"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.9, "max. |eta| of tracks"};
  Configurable<int> nMinTpcClusters{"nMinTpcClusters", 60, "min. number of TPC clusters"};

  // Axes (the mass axis matches the default rho window; override it for the J/psi)
  ConfigurableAxis axisMass{"axisMass", {200, 0.5, 1.3}, "m_{#pi#pi} (GeV/#it{c}^{2})"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.175, 0.20, 0.25, 0.30, 0.40, 0.50}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisY{"axisY", {100, -1., 1.}, "Rapidity"};
  ConfigurableAxis axisPtTrk{"axisPtTrk", {80, 0.0, 0.8}, "#it{p}_{T} (GeV/#it{c})"};

  HistogramConfigSpec hVtxZ{HistType::kTH1F, {{80, -20., 20.}}};
  HistogramConfigSpec hVtxX{HistType::kTH1F, {{50, -0.04, -0.015}}};
  HistogramConfigSpec hVtxY{HistType::kTH1F, {{60, -0.015, 0.015}}};
  HistogramConfigSpec hEtaTrk{HistType::kTH1F, {{100, -1.0, 1.0}}};

  HistogramRegistry registry{"registry"};

  // services
  Service<o2::ccdb::BasicCCDBManager> ccdb; // access to database

  // variables to store run info
  int runNumberBc = 0;     // run number used to process BCs
  int runNumberCol = 0;    // run number used to process collisions
  int64_t sor = 0;         // best known timestamp for the start of run
  int64_t orbitsPerTF = 0; // number of orbits per TF
  int64_t bcSOR = 0;       // first bc of the first orbit
  int64_t nBCsPerTF = 0;   // duration of TF in bcs
  int64_t currentTF = -1;  // current time frame being looked at
  int64_t nTF = 0;         // number of time frames in run

  // constants
  const unsigned int twoTracks = 2;
  const int unassignedDCA = 999;

  static constexpr float MassPion = o2::constants::physics::MassPionCharged;
  static constexpr float MassMuon = o2::constants::physics::MassMuon;

  // folders that contains plots for all tracks and associated tracks
  static constexpr std::array<const char*, 2> TrackFolders = {"Trk/", "TrkColl/"};

  // Name shortenings
  // passed* columns are in TrackSelectionExtension; isGlobalTrack* and trackCutFlag in TrackSelection
  using TracksExtraSels = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                    aod::TrackSelection, aod::TrackSelectionExtension>;
  using TracksExtraWPidPi = soa::Join<TracksExtraSels, aod::pidTPCFullPi>;
  using TracksExtraWPidMu = soa::Join<TracksExtraSels, aod::pidTPCFullMu>;
  using TracksExtraWPidPiMc = soa::Join<TracksExtraWPidPi, aod::McTrackLabels>;
  using TracksExtraWPidMuMc = soa::Join<TracksExtraWPidMu, aod::McTrackLabels>;
  using CollisionsWithMc = soa::Join<aod::Collisions, aod::McCollisionLabels>;
  using TracksExtraIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA,
                                  aod::TrackSelection, aod::TrackSelectionExtension>;

  void init(InitContext&)
  {
    // collision level
    registry.add("Coll/hNContrib", ";N_{PV contributors};entries", {HistType::kTH1F, {{10, -0.5, 9.5}}});
    registry.add("Coll/hBCid", ";BCid;entries", {HistType::kTH1F, {{1000, 0., 100000.}}});
    registry.add("Coll/hVtxZ", ";#it{z}_{vtx} (cm);entries", hVtxZ);
    registry.add("Coll/hnTF", ";TF;entries", {HistType::kTH1F, {{2000, 0., 2000}}});
    registry.add("Coll/hVtxX", ";#it{x}_{vtx} (cm);entries", hVtxX);
    registry.add("Coll/hVtxY", ";#it{y}_{vtx} (cm);entries", hVtxY);
    registry.add("Coll/hVtxChi2", ";#chi^{2} vtx;entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // track level (all good tracks, before grouping by collision)
    registry.add("Trk/hPt", ";#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtTrk}});
    registry.add("Trk/hEta", ";#eta;entries", hEtaTrk);
    registry.add("Trk/hChi2NCl", ";#chi^{2}/N_{cls} TPC;entries", {HistType::kTH1F, {{50, 0., 5.}}});
    registry.add("Trk/hTpcSignalVsP", ";#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x}", {HistType::kTH2F, {{200, 0., 2.}, {300, 0., 300.}}});
    registry.add("Trk/hNSigmaVsP", ";#it{p} (GeV/#it{c});n#sigma^{TPC}", {HistType::kTH2F, {{200, 0., 2.}, {100, -10., 10.}}});
    registry.add("Trk/hHasIts", ";has ITS;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}});
    registry.add("Trk/hIsPvContrib", ";is PV contributor;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}});
    registry.add("Trk/hTpcNClsFound", ";N_{cls} TPC;entries", {HistType::kTH1F, {{90, 60, 150.}}});
    registry.add("Trk/hItsChi2NCl", ";#chi^{2}/N_{cls} ITS;entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("Trk/hItsNCls", ";N_{cls} ITS;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}});
    registry.add("Trk/hItsNClsInnerBarrel", ";N_{cls} ITS Inner Barrel;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}});
    registry.add("Trk/hDcaXY", ";DCA_{xy} (cm);entries", {HistType::kTH1F, {{140, -0.035, 0.035}}});
    registry.add("Trk/hDcaZ", ";DCA_{z} (cm);entries", {HistType::kTH1F, {{200, -0.05, 0.05}}});

    // track level (prongs of the candidate, after collision matching)
    registry.addClone("Trk/", "TrkColl/");
    registry.add("TrkColl/hTrackSizePerEvent", ";#it{N}_{tracks} per event;entries", {HistType::kTH1F, {{10, -0.5, 9.5}}});

    // collision Id counter
    registry.add("Coll/hCollIdCounter", ";collision ID;entries", {HistType::kTH1F, {{3, -1.5, 1.5}}});

    // candidate level plots
    registry.add("Cand/hMass", ";m_{#pi#pi} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMass}});
    registry.add("Cand/hPt", ";#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Cand/hRapidity", ";#it{y};entries", {HistType::kTH1F, {axisY}});
    registry.add("Cand/hPtVsRapidity", ";#it{p}_{T} (GeV/#it{c});#it{y}", {HistType::kTH2F, {axisPt, axisY}});

    // MC histograms are only booked if at least one MC process function is enabled
    const bool needMc = doprocessRhoTracksBeforeGroupingMcInfo || doprocessJpsiTracksBeforeGroupingMcInfo ||
                        doprocessRhoCandMcInfo || doprocessJpsiCandMcInfo || doprocessMcGenExclusive || doprocessMcGenEmbedded;
    if (!needMc) {
      return;
    }

    // MC reco matched plots
    registry.addClone("Cand/", "Matched/Cand/");

    registry.addClone("Coll/hNContrib", "Matched/Coll/hNContrib");
    registry.addClone("Coll/hVtxZ", "Matched/Coll/hVtxZ");
    registry.addClone("Coll/hVtxX", "Matched/Coll/hVtxX");
    registry.addClone("Coll/hVtxY", "Matched/Coll/hVtxY");

    registry.add("Matched/Trk/hPt", ";#it{p}_{T}^{MC matched} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtTrk}});
    registry.add("Matched/Trk/hEta", ";#eta^{MC};entries", hEtaTrk);
    registry.add("Matched/Trk/hMcIsPhysicalPrimary", ";is physical primary;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}});
    registry.add("Matched/Trk/hMotherPdg", ";mother PDG code;entries", {HistType::kTH1F, {{1001, 0.5, 1000.5}}});

    registry.addClone("Matched/Trk/", "Matched/TrkColl/");
  
    // MC generated plots
    registry.addClone("Matched/Coll/", "McGen/Coll/");

    registry.addClone("Matched/Trk/", "McGen/Part/");

    registry.add("McGen/MotherPart/hMotherPdg", ";mother PDG code;entries", {HistType::kTH1F, {{1001, 0.5, 1000.5}}});
    registry.add("McGen/MotherPart/hPt", ";#it{p}_{T}^{MC} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("McGen/MotherPart/hRapidity", ";#it{y}^{MC};entries", {HistType::kTH1F, {axisY}});
    registry.add("McGen/MotherPart/hPtVsRapidity", ";#it{p}_{T}^{MC} (GeV/#it{c});#it{y}^{MC}", {HistType::kTH2F, {axisPt, axisY}});

    registry.addClone("TrkColl/hTrackSizePerEvent", "McGen/Part/hTrackSizePerEvent");

  } // end init()

  // FUNCTIONS: all functions defined here

  //--------------------------------------------------------------------------------
  // get info on the run, given the run number. Fills info needed to use getTimeFrame()
  void getRunInfo(int run)
  {
    auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(ccdb->instance(), run);
    sor = runInfo.sor; // in ms
    auto orbitSOR = runInfo.orbitSOR;
    auto orbitEOR = runInfo.orbitEOR;
    orbitsPerTF = runInfo.orbitsPerTF;
    bcSOR = orbitSOR * o2::constants::lhc::LHCMaxBunches;        // first bc of the first orbit
    nBCsPerTF = orbitsPerTF * o2::constants::lhc::LHCMaxBunches; // duration of TF in bcs
    nTF = std::ceil((orbitEOR - orbitSOR) / orbitsPerTF);
  }

  //--------------------------------------------------------------------------------
  // compute and return the TF for the given BC
  int64_t getTimeFrame(int64_t globalBC)
  {
    return (globalBC - bcSOR) / nBCsPerTF;
  }

  //--------------------------------------------------------------------------------
  // Basic single-track selection
  template <typename TTrack>
  bool isGoodTrack(TTrack const& track)
  {
    if (!track.hasTPC()) {
      return false;
    }
    if (track.tpcNClsFound() < nMinTpcClusters) {
      return false;
    }
    if (track.pt() < ptTrackMin || std::abs(track.eta()) > etaTrackMax) {
      return false;
    }
    if (std::abs(track.dcaZ()) >= unassignedDCA || std::abs(track.dcaXY()) >= unassignedDCA) {
      return false;
    }
    return true;
  }

  //--------------------------------------------------------------------------------
  // Select the PID column at compile time
  template <CandSpecies species, typename TTrack>
  float getNSigma(TTrack const& track)
  {
    if constexpr (species == CandSpecies::kRho) {
      return track.tpcNSigmaPi();
    } else {
      return track.tpcNSigmaMu();
    }
  }

  //--------------------------------------------------------------------------------
  // fills the histograms with info on track properties
  // works for both all and associated tracks, and if MC is available or not
  template <CandSpecies species, int folder, bool isMc, typename TTrack>
  void fillTrackHistograms(TTrack const& track)
  {
    // fill the histos only for good tracks
    if (!isGoodTrack(track)) {
      return;
    }

    registry.fill(HIST(TrackFolders[folder]) + HIST("hPt"), track.pt());
    registry.fill(HIST(TrackFolders[folder]) + HIST("hEta"), track.eta());
    registry.fill(HIST(TrackFolders[folder]) + HIST("hChi2NCl"), track.tpcChi2NCl());
    registry.fill(HIST(TrackFolders[folder]) + HIST("hTpcSignalVsP"), track.p(), track.tpcSignal());
    registry.fill(HIST(TrackFolders[folder]) + HIST("hNSigmaVsP"), track.p(), getNSigma<species>(track));
    registry.fill(HIST(TrackFolders[folder]) + HIST("hHasIts"), static_cast<int>(track.hasITS()));
    registry.fill(HIST(TrackFolders[folder]) + HIST("hIsPvContrib"), static_cast<int>(track.isPVContributor()));
    registry.fill(HIST(TrackFolders[folder]) + HIST("hTpcNClsFound"), track.tpcNClsFound());
    registry.fill(HIST(TrackFolders[folder]) + HIST("hItsChi2NCl"), track.itsChi2NCl());
    registry.fill(HIST(TrackFolders[folder]) + HIST("hItsNCls"), track.itsNCls());
    registry.fill(HIST(TrackFolders[folder]) + HIST("hItsNClsInnerBarrel"), track.itsNClsInnerBarrel());
    registry.fill(HIST(TrackFolders[folder]) + HIST("hDcaXY"), track.dcaXY());
    registry.fill(HIST(TrackFolders[folder]) + HIST("hDcaZ"), track.dcaZ());

    // if the MC info are available, fill the corresponding histograms
    if constexpr (isMc) {
      if (!track.has_mcParticle()) {
        return;
      }
      auto mcPart = track.template mcParticle_as<aod::McParticles>();
      registry.fill(HIST("Matched/") + HIST(TrackFolders[folder]) + HIST("hPt"), mcPart.pt());
      registry.fill(HIST("Matched/") + HIST(TrackFolders[folder]) + HIST("hEta"), mcPart.eta());
      registry.fill(HIST("Matched/") + HIST(TrackFolders[folder]) + HIST("hMcIsPhysicalPrimary"), static_cast<int>(mcPart.isPhysicalPrimary()));
      if (mcPart.has_mothers()) {
        for (auto const& mcMother : mcPart.template mothers_as<aod::McParticles>()) {
          registry.fill(HIST("Matched/") + HIST(TrackFolders[folder]) + HIST("hMotherPdg"), std::abs(mcMother.pdgCode()));
        }
      }
    }
  }

  //--------------------------------------------------------------------------------
  // loop over tracks and fill the histograms with info on track properties
  template <CandSpecies species, int folder, bool isMc, typename TTrack>
  void loopOverTracksBeforeGrouping(TTrack const& tracks){
    for (auto const& track : tracks) {
      fillTrackHistograms<species, folder, isMc>(track);
    }
  }

  //--------------------------------------------------------------------------------
  // fills the histograms with info on candidate properties
  template <CandSpecies species, bool isMc, typename TTrack, typename TCollision>
  void fillCandidateHistograms(TCollision const& collision, TTrack const& tracks, float candMass)
  {
    registry.fill(HIST("Coll/hNContrib"), collision.numContrib());
    registry.fill(HIST("Coll/hBCid"), collision.bcId());
    registry.fill(HIST("Coll/hVtxZ"), collision.posZ());
    registry.fill(HIST("Coll/hVtxX"), collision.posX());
    registry.fill(HIST("Coll/hVtxY"), collision.posY());
    registry.fill(HIST("Coll/hVtxChi2"), collision.chi2());
    if constexpr (isMc) {
      if (collision.has_mcCollision()) {
        auto mcColl = collision.template mcCollision_as<aod::McCollisions>();
        registry.fill(HIST("Matched/Coll/hNContrib"), collision.numContrib());
        registry.fill(HIST("Matched/Coll/hVtxX"), mcColl.posX());
        registry.fill(HIST("Matched/Coll/hVtxY"), mcColl.posY());
        registry.fill(HIST("Matched/Coll/hVtxZ"), mcColl.posZ());
      }
    }

    // Select tracks for the candidate
    std::vector<decltype(tracks.begin())> goodTracks;
    for (auto const& track : tracks) {
      if (isGoodTrack(track)) {
        goodTracks.push_back(track);
        // use TrkColl directory
        fillTrackHistograms<species, 1, isMc>(track);
      }
    }

    registry.fill(HIST("TrkColl/hTrackSizePerEvent"), goodTracks.size());

    // build UPC candidates using events with exactly two opposite sign good tracks
    if (goodTracks.size() != twoTracks) {
      return;
    }
    auto const& track0 = goodTracks[0];
    auto const& track1 = goodTracks[1];
    if (track0.sign() * track1.sign() >= 0) {
      return;
    }

    ROOT::Math::PxPyPzMVector p0(track0.px(), track0.py(), track0.pz(), candMass);
    ROOT::Math::PxPyPzMVector p1(track1.px(), track1.py(), track1.pz(), candMass);
    auto candidate = p0 + p1;

    registry.fill(HIST("Cand/hMass"), candidate.M());
    registry.fill(HIST("Cand/hPt"), candidate.Pt());
    registry.fill(HIST("Cand/hRapidity"), candidate.Rapidity());
    registry.fill(HIST("Cand/hPtVsRapidity"), candidate.Pt(), candidate.Rapidity());

    if constexpr (isMc) {
      if (track0.has_mcParticle() && track1.has_mcParticle()) {
        auto mcPart0 = track0.template mcParticle_as<aod::McParticles>();
        auto mcPart1 = track1.template mcParticle_as<aod::McParticles>();
        ROOT::Math::PxPyPzMVector mcP0(mcPart0.px(), mcPart0.py(), mcPart0.pz(), candMass);
        ROOT::Math::PxPyPzMVector mcP1(mcPart1.px(), mcPart1.py(), mcPart1.pz(), candMass);
        auto mcCandidate = mcP0 + mcP1;
        registry.fill(HIST("Matched/Cand/hMass"), mcCandidate.M());
        registry.fill(HIST("Matched/Cand/hPt"), mcCandidate.Pt());
        registry.fill(HIST("Matched/Cand/hRapidity"), mcCandidate.Rapidity());
        registry.fill(HIST("Matched/Cand/hPtVsRapidity"), mcCandidate.Pt(), mcCandidate.Rapidity());
      }
    }
  }

  //--------------------------------------------------------------------------------
  // generated (MC truth) level: daughters and their mothers
  template <bool isEmbedded>
  void checkMcGen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles){

    registry.fill(HIST("McGen/Coll/hVtxX"), mcCollision.posX());
    registry.fill(HIST("McGen/Coll/hVtxY"), mcCollision.posY());
    registry.fill(HIST("McGen/Coll/hVtxZ"), mcCollision.posZ());

    std::vector<aod::McParticle> goodMcGenParticles;
    for (auto const& mcPart : mcParticles) {
      
      if (!mcPart.isPhysicalPrimary() || mcPart.pt() < ptTrackMin || std::abs(mcPart.eta()) > etaTrackMax) {
        continue;      
      }
      registry.fill(HIST("McGen/Part/hPt"), mcPart.pt());
      registry.fill(HIST("McGen/Part/hEta"), mcPart.eta());
      
      goodMcGenParticles.push_back(mcPart);
    }
    registry.fill(HIST("McGen/Part/hTrackSizePerEvent"), goodMcGenParticles.size());

    // consider only j/psi->mumu and rho->pipi events
    auto genId = mcCollision.getGeneratorId();
    LOGF(info, "Gen Id: {%d}", genId);
    auto sourceId = mcCollision.getSourceId();
    LOGF(info, "Source Id: {%d}", sourceId);

    if constexpr (isEmbedded) {
      if (!sourceId) {
        return;
      }
    } else {
      if (sourceId) {
        LOGF(fatal, "Source Id: {%d} is not zero for non-embedded MC", sourceId);
        return;
      }
    }
    
    if (goodMcGenParticles.size() != twoTracks) {
      return;
    }

    // we have 2 good tracks, with the same mother, taking it from the first one    
    auto mothers = goodMcGenParticles[0].template mothers_as<aod::McParticles>();

    if (mothers.empty()) {
      LOGF(fatal, "Selected particle has no mother: pdgCode=%d", goodMcGenParticles[0].pdgCode());
      return;
    }
    auto motherPart = mothers[0];

    registry.fill(HIST("McGen/MotherPart/hMotherPdg"), std::abs(motherPart.pdgCode()));
    registry.fill(HIST("McGen/MotherPart/hPt"), motherPart.pt());
    registry.fill(HIST("McGen/MotherPart/hRapidity"), motherPart.y());
    registry.fill(HIST("McGen/MotherPart/hPtVsRapidity"), motherPart.pt(), motherPart.y());
    
  }

  //--------------------------------------------------------------------------------
  // PROCESS FUNCTIONS
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // check collisionId and TF
  void processTracksIU(TracksExtraIU const& tracks, aod::Collisions const&, aod::BCs const&) {
    // track loop
    for (const auto& track : tracks) {

      if (!isGoodTrack(track)) {
        continue;
      }

      int32_t collId = track.collisionId();

      if (collId < 0) {
        registry.fill(HIST("Coll/hCollIdCounter"), -1);
      }
      else if (collId > 0) {
        registry.fill(HIST("Coll/hCollIdCounter"), 1);
      } else {
        registry.fill(HIST("Coll/hCollIdCounter"), 0);
      }

      if (collId < 0) {
        continue;
      }

      auto collision = track.collision();
      auto bc = collision.bc();
      
      getRunInfo(bc.runNumber());

      registry.fill(HIST("Coll/hnTF"), getTimeFrame(bc.globalBC()));
    }
  }

  // look at all tracks
  void processRhoTracksBeforeGrouping(TracksExtraWPidPi const& tracks)
  {
    loopOverTracksBeforeGrouping<CandSpecies::kRho, 0, false>(tracks);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processRhoTracksBeforeGrouping, "Process rho tracks before asking for collisions", true);

  void processJpsiTracksBeforeGrouping(TracksExtraWPidMu const& tracks)
  {
    loopOverTracksBeforeGrouping<CandSpecies::kJpsi, 0, false>(tracks);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processJpsiTracksBeforeGrouping, "Process J/Psi tracks before asking for collisions", false);

  void processRhoTracksBeforeGroupingMcInfo(TracksExtraWPidPiMc const& tracks, aod::McParticles const&)
  {
    loopOverTracksBeforeGrouping<CandSpecies::kRho, 0, true>(tracks);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processRhoTracksBeforeGroupingMcInfo, "Process rho tracks before asking for collisions", false);

  void processJpsiTracksBeforeGroupingMcInfo(TracksExtraWPidMuMc const& tracks, aod::McParticles const&)
  {
    loopOverTracksBeforeGrouping<CandSpecies::kJpsi, 0, true>(tracks);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processJpsiTracksBeforeGroupingMcInfo, "Process J/Psi tracks before asking for collisions", false);

  // look at matched tracks
  void processRhoCand(aod::Collision const& collision, TracksExtraWPidPi const& tracks)
  {
    fillCandidateHistograms<CandSpecies::kRho, false>(collision, tracks, MassPion);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processRhoCand, "Rho -> pi pi candidates and track QA", true);

  void processJpsiCand(aod::Collision const& collision, TracksExtraWPidMu const& tracks)
  {
    fillCandidateHistograms<CandSpecies::kJpsi, false>(collision, tracks, MassMuon);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processJpsiCand, "J/Psi -> mu mu candidates and track QA", false);

  void processRhoCandMcInfo(CollisionsWithMc::iterator const& collision, TracksExtraWPidPiMc const& tracks,
                            aod::McParticles const&, aod::McCollisions const&)
  {
    fillCandidateHistograms<CandSpecies::kRho, true>(collision, tracks, MassPion);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processRhoCandMcInfo, "Rho -> pi pi candidates and track QA", false);

  void processJpsiCandMcInfo(CollisionsWithMc::iterator const& collision, TracksExtraWPidMuMc const& tracks,
                             aod::McParticles const&, aod::McCollisions const&)
  {
    fillCandidateHistograms<CandSpecies::kJpsi, true>(collision, tracks, MassMuon);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processJpsiCandMcInfo, "J/Psi -> mu mu candidates and track QA", false);

  // look at gen MC
  void processMcGenEmbedded(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    checkMcGen<true>(mcCollision, mcParticles);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processMcGenEmbedded, "Generated embedded events", false);

  void processMcGenExclusive(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    checkMcGen<false>(mcCollision, mcParticles);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processMcGenExclusive, "Generated exclusive events", false);

  PROCESS_SWITCH(UpcTrackVertexingQA, processTracksIU, "Track QA for IU tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UpcTrackVertexingQA>(cfgc)};
}