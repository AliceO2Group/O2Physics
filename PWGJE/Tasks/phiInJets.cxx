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

/// \file phiInJets.cxx
/// \brief Reconstruction of Phi yield through track-track Minv correlations for resonance hadrochemistry analysis.
///
///
/// \author Adrian Fereydon Nassirpour <adrian.fereydon.nassirpour@cern.ch>

#include <TLorentzVector.h>
#include <TVector2.h>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
//#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGLF/DataModel/LFResonanceTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phiInJets {
  SliceCache cache;
  // Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry JEhistos{"JEhistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry LFhistos{"LFhistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::string> cfgeventSelections{"cfgeventSelections", "sel8", "choose event selection"};
  Configurable<std::string> cfgtrackSelections{"cfgtrackSelections", "globalTracks", "set track selections"};

  Configurable<double> cfgtrkMinPt{"cfgtrkMinPt", 0.15, "set track min pT"};
  Configurable<double> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgConnectedToPV{"cfgConnectedToPV", true, "PV contributor track selection"};           // PV Contriuibutor
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<double> cfgnFindableTPCClusters{"cfgnFindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfgnTPCCrossedRows{"cfgnTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgnRowsOverFindable{"cfgnRowsOverFindable", 1.2, "nRowsOverFindable TPC CLusters"};
  Configurable<double> cfgnTPCChi2{"cfgnTPChi2", 4.0, "nTPC Chi2 per Cluster"};
  Configurable<double> cfgnITSChi2{"cfgnITShi2", 36.0, "nITS Chi2 per Cluster"};
  Configurable<int> cfgnTPCPID{"cfgnTPCPID", 4, "nTPC PID"};
  Configurable<int> cfgnTOFPID{"cfgnTOFPID", 4, "nTOF PID"};
  Configurable<float> cfgjetPtMin{"cfgjetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> cfgjetR{"cfgjetR", 0.4, "jet resolution parameter"};
  // CONFIG DONE
  /////////////////////////////////////////  //INIT

  int eventSelection = -1;

  void init(o2::framework::InitContext&)
  {
    // HISTOGRAMS
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPhi{200, -1, +7, "#phi"};
    const AxisSpec axisPt{200, 0, +200, "#pt"};
    const AxisSpec MinvAxis = {400, 0.95, 1.35};
    const AxisSpec PtAxis = {200, 0, 20.0};
    const AxisSpec MultAxis = {100, 0, 100};

    JEhistos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    JEhistos.add("phiHistogram", "phiHistogram", kTH1F, {axisPhi});
    JEhistos.add("FJetaHistogram", "FJetaHistogram", kTH1F, {axisEta});
    JEhistos.add("FJphiHistogram", "FJphiHistogram", kTH1F, {axisPhi});
    JEhistos.add("FJptHistogram", "FJptHistogram", kTH1F, {axisPt});

    JEhistos.add("nEvents", "nEvents", kTH1F, {{4, 0.0, 4.0}});
    JEhistos.add("hDCArToPv", "DCArToPv", kTH1F, {{300, 0.0, 3.0}});
    JEhistos.add("hDCAzToPv", "DCAzToPv", kTH1F, {{300, 0.0, 3.0}});
    JEhistos.add("rawpT", "rawpT", kTH1F, {{1000, 0.0, 10.0}});
    JEhistos.add("rawDpT", "rawDpT", kTH2F, {{1000, 0.0, 10.0}, {300, -1.5, 1.5}});
    JEhistos.add("hIsPrim", "hIsPrim", kTH1F, {{2, -0.5, +1.5}});
    JEhistos.add("hIsGood", "hIsGood", kTH1F, {{2, -0.5, +1.5}});
    JEhistos.add("hIsPrimCont", "hIsPrimCont", kTH1F, {{2, -0.5, +1.5}});
    JEhistos.add("hFindableTPCClusters", "hFindableTPCClusters", kTH1F, {{200, 0, 200}});
    JEhistos.add("hFindableTPCRows", "hFindableTPCRows", kTH1F, {{200, 0, 200}});
    JEhistos.add("hClustersVsRows", "hClustersVsRows", kTH1F, {{200, 0, 2}});
    JEhistos.add("hTPCChi2", "hTPCChi2", kTH1F, {{200, 0, 100}});
    JEhistos.add("hITSChi2", "hITSChi2", kTH1F, {{200, 0, 100}});

    JEhistos.add("hUSS", "hUSS", kTH3F, {MultAxis, PtAxis, MinvAxis});
    JEhistos.add("hUSS_1D", "hUSS_1D", kTH1F, {MinvAxis});
    JEhistos.add("hUSS_1D_2_3", "hUSS_1D_2_3", kTH1F, {MinvAxis});

    JEhistos.add("hLSS", "hLSS", kTH3F, {MultAxis, PtAxis, MinvAxis});
    JEhistos.add("hLSS_1D", "hLSS_1D", kTH1F, {MinvAxis});
    JEhistos.add("hLSS_1D_2_3", "hLSS_1D_2_3", kTH1F, {MinvAxis});

    JEhistos.add("hUSS_INSIDE", "hUSS_INSIDE", kTH3F, {MultAxis, PtAxis, MinvAxis});
    JEhistos.add("hUSS_INSIDE_1D", "hUSS_INSIDE_1D", kTH1F, {MinvAxis});
    JEhistos.add("hUSS_INSIDE_1D_2_3", "hUSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});

    JEhistos.add("hLSS_INSIDE", "hLSS_INSIDE", kTH3F, {MultAxis, PtAxis, MinvAxis});
    JEhistos.add("hLSS_INSIDE_1D", "hLSS_INSIDE_1D", kTH1F, {MinvAxis});
    JEhistos.add("hLSS_INSIDE_1D_2_3", "hLSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});

    JEhistos.add("hUSS_OUTSIDE", "hUSS_OUTSIDE", kTH3F, {MultAxis, PtAxis, MinvAxis});
    JEhistos.add("hUSS_OUTSIDE_1D", "hUSS_OUTSIDE_1D", kTH1F, {MinvAxis});
    JEhistos.add("hUSS_OUTSIDE_1D_2_3", "hUSS_OUTSIDE_1D_2_3", kTH1F, {MinvAxis});

    JEhistos.add("hLSS_OUTSIDE", "hLSS_OUTSIDE", kTH3F, {MultAxis, PtAxis, MinvAxis});
    JEhistos.add("hLSS_OUTSIDE_1D", "hLSS_OUTSIDE_1D", kTH1F, {MinvAxis});
    JEhistos.add("hLSS_OUTSIDE_1D_2_3", "hLSS_OUTSIDE_1D_2_3", kTH1F, {MinvAxis});

    JEhistos.add("hMultFT0M", "hMultFT0M", kTH1F, {MultAxis});

    LFhistos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    LFhistos.add("phiHistogram", "phiHistogram", kTH1F, {axisPhi});
    LFhistos.add("nEvents", "nEvents", kTH1F, {{4, 0.0, 4.0}});
    LFhistos.add("hDCArToPv", "DCArToPv", kTH1F, {{300, 0.0, 3.0}});
    LFhistos.add("hDCAzToPv", "DCAzToPv", kTH1F, {{300, 0.0, 3.0}});
    LFhistos.add("rawpT", "rawpT", kTH1F, {{1000, 0.0, 10.0}});
    LFhistos.add("rawDpT", "rawDpT", kTH2F, {{1000, 0.0, 10.0}, {300, -1.5, 1.5}});
    LFhistos.add("hIsPrim", "hIsPrim", kTH1F, {{2, -0.5, +1.5}});
    LFhistos.add("hIsGood", "hIsGood", kTH1F, {{2, -0.5, +1.5}});
    LFhistos.add("hIsPrimCont", "hIsPrimCont", kTH1F, {{2, -0.5, +1.5}});
    LFhistos.add("hFindableTPCClusters", "hFindableTPCClusters", kTH1F, {{200, 0, 200}});
    LFhistos.add("hFindableTPCRows", "hFindableTPCRows", kTH1F, {{200, 0, 200}});
    LFhistos.add("hClustersVsRows", "hClustersVsRows", kTH1F, {{200, 0, 2}});
    LFhistos.add("hTPCChi2", "hTPCChi2", kTH1F, {{200, 0, 100}});
    LFhistos.add("hITSChi2", "hITSChi2", kTH1F, {{200, 0, 100}});

    LFhistos.add("hUSS", "hUSS", kTH3F, {MultAxis, PtAxis, MinvAxis});
    LFhistos.add("hUSS_1D", "hUSS_1D", kTH1F, {MinvAxis});
    LFhistos.add("hUSS_1D_2_3", "hUSS_1D_2_3", kTH1F, {MinvAxis});

    LFhistos.add("hLSS", "hLSS", kTH3F, {MultAxis, PtAxis, MinvAxis});
    LFhistos.add("hLSS_1D", "hLSS_1D", kTH1F, {MinvAxis});
    LFhistos.add("hLSS_1D_2_3", "hLSS_1D_2_3", kTH1F, {MinvAxis});

    LFhistos.add("hMultFT0M", "hMultFT0M", kTH1F, {MultAxis});

    // EVENT SELECTION
    eventSelection = JetDerivedDataUtilities::initialiseEventSelection(static_cast<std::string>(cfgeventSelections));

  } // end of init

  double massKa = o2::constants::physics::MassKPlus;
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs>; // , aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullKa, aod::pidTOFFullKa>;
  Filter jetCuts = aod::jet::pt > cfgjetPtMin&& aod::jet::r == nround(cfgjetR.node() * 100.0f);

  // Function for track quality cuts
  template <typename TrackType>
  bool trackSelection(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cfgtrkMinPt)
      return false;

    if (std::abs(track.dcaXY()) > cfgMaxDCArToPVcut)
      return false;

    if (std::abs(track.dcaZ()) > cfgMaxDCAzToPVcut)
      return false;

    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;

    if (std::abs(track.tpcNClsFindable()) < cfgnFindableTPCClusters)
      return false;

    if (std::abs(track.tpcNClsCrossedRows()) < cfgnTPCCrossedRows)
      return false;

    if (std::abs(track.tpcCrossedRowsOverFindableCls()) > cfgnRowsOverFindable)
      return false;

    if (std::abs(track.tpcChi2NCl() > cfgnTPCChi2))
      return false;

    if (std::abs(track.itsChi2NCl() > cfgnITSChi2))
      return false;

    if (cfgConnectedToPV && !track.isPVContributor())
      return false;
    return true;
  };

  template <typename T>
  bool trackPID(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaKa()) < cfgnTPCPID)
      tpcPIDPassed = true;

    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaKa()) < cfgnTOFPID) {
        tofPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

  template <bool IsMC, bool IsMix, typename TracksType, typename JetType>
  void minvReconstruction(double mult, const TracksType& trk1, const TracksType& trk2, const JetType& jets)
  {
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;

    if (!trackSelection(trk1) || !trackSelection(trk2))
      return;

    if (!trackPID(trk1) || !trackPID(trk2))
      return;

    if (trk1.index() == trk2.index())
      return; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

    lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
    lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
    lResonance = lDecayDaughter1 + lDecayDaughter2;

    if (std::abs(lResonance.Rapidity()) > 0.5)
      return;

    if (trk1.sign() * trk2.sign() < 0) {
      JEhistos.fill(HIST("hUSS_1D"), lResonance.M());

      if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
        JEhistos.fill(HIST("hUSS_1D_2_3"), lResonance.M());

      JEhistos.fill(HIST("hUSS"), mult, lResonance.Pt(), lResonance.M());
    } else if (trk1.sign() * trk2.sign() > 0) {

      JEhistos.fill(HIST("hLSS_1D"), lResonance.M());
      if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
        JEhistos.fill(HIST("hLSS_1D_2_3"), lResonance.M());

      JEhistos.fill(HIST("hLSS"), mult, lResonance.Pt(), lResonance.M());
    }

    bool jetFlag = false;
    for (auto const& jet : jets) {

      double phidiff = TVector2::Phi_mpi_pi(jet.phi() - lResonance.Phi());
      double etadiff = jet.eta() - lResonance.Eta();
      double R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));
      if (R < cfgjetR)
        jetFlag = true;
    }
    if (jetFlag) {
      if (trk1.sign() * trk2.sign() < 0) {
        JEhistos.fill(HIST("hUSS_INSIDE_1D"), lResonance.M());

        if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
          JEhistos.fill(HIST("hUSS_INSIDE_1D_2_3"), lResonance.M());

        JEhistos.fill(HIST("hUSS_INSIDE"), mult, lResonance.Pt(), lResonance.M());
      } else if (trk1.sign() * trk2.sign() > 0) {

        JEhistos.fill(HIST("hLSS_INSIDE_1D"), lResonance.M());
        if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
          JEhistos.fill(HIST("hLSS_INSIDE_1D_2_3"), lResonance.M());

        JEhistos.fill(HIST("hLSS_INSIDE"), mult, lResonance.Pt(), lResonance.M());
      }
    } // jetflag
    if (!jetFlag) {
      if (trk1.sign() * trk2.sign() < 0) {
        JEhistos.fill(HIST("hUSS_OUTSIDE_1D"), lResonance.M());

        if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
          JEhistos.fill(HIST("hUSS_OUTSIDE_1D_2_3"), lResonance.M());

        JEhistos.fill(HIST("hUSS_OUTSIDE"), mult, lResonance.Pt(), lResonance.M());
      } else if (trk1.sign() * trk2.sign() > 0) {

        JEhistos.fill(HIST("hLSS_OUTSIDE_1D"), lResonance.M());
        if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
          JEhistos.fill(HIST("hLSS_OUTSIDE_1D_2_3"), lResonance.M());

        JEhistos.fill(HIST("hLSS_OUTSIDE"), mult, lResonance.Pt(), lResonance.M());
      }
    } //! jetflag
  }   // MinvReconstruction

  int nEvents = 0;
  void processJetTracks(aod::JCollision const& collision, soa::Filtered<aod::FullJets> const& fulljets, soa::Join<aod::JTracks, aod::JTrackPIs> const& tracks, TrackCandidates const&)
  {
    JEhistos.fill(HIST("nEvents"), 0.5);

    if (!JetDerivedDataUtilities::selectCollision(collision, eventSelection))
      return;

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      auto trk1 = track1.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>();
      auto trk2 = track2.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>();
      minvReconstruction<false, false>(1.0, trk1, trk2, fulljets);
    }
    for (auto fulljet : fulljets) {
      JEhistos.fill(HIST("FJetaHistogram"), fulljet.eta());
      JEhistos.fill(HIST("FJphiHistogram"), fulljet.phi());
      JEhistos.fill(HIST("FJptHistogram"), fulljet.phi());
    }

    JEhistos.fill(HIST("nEvents"), 1.5);
    for (auto const& track : tracks) {
      auto originalTrack = track.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>();
      JEhistos.fill(HIST("hDCArToPv"), originalTrack.dcaXY());
      JEhistos.fill(HIST("hDCAzToPv"), originalTrack.dcaZ());
      JEhistos.fill(HIST("rawpT"), originalTrack.pt());
      JEhistos.fill(HIST("rawDpT"), track.pt(), track.pt() - originalTrack.pt());
      JEhistos.fill(HIST("hIsPrim"), originalTrack.isPrimaryTrack());
      JEhistos.fill(HIST("hIsGood"), originalTrack.isGlobalTrackWoDCA());
      JEhistos.fill(HIST("hIsPrimCont"), originalTrack.isPVContributor());
      JEhistos.fill(HIST("hFindableTPCClusters"), originalTrack.tpcNClsFindable());
      JEhistos.fill(HIST("hFindableTPCRows"), originalTrack.tpcNClsCrossedRows());
      JEhistos.fill(HIST("hClustersVsRows"), originalTrack.tpcCrossedRowsOverFindableCls());
      JEhistos.fill(HIST("hTPCChi2"), originalTrack.tpcChi2NCl());
      JEhistos.fill(HIST("hITSChi2"), originalTrack.itsChi2NCl());

      if (!trackSelection(originalTrack))
        continue;

      JEhistos.fill(HIST("etaHistogram"), track.eta());
      JEhistos.fill(HIST("phiHistogram"), track.phi());
    } // JTrack Loop

  }; // Process Switch
  PROCESS_SWITCH(phiInJets, processJetTracks, "process JE Framework", true);

  void processLFTracks(EventCandidates::iterator const& collision, TrackCandidates const& tracks)
  {

    LFhistos.fill(HIST("nEvents"), 0.5);
    if (!collision.sel8())
      return;

    LFhistos.fill(HIST("hMultFT0M"), collision.multFT0M());

    LFhistos.fill(HIST("nEvents"), 1.5);
    for (auto const& track : tracks) {
      LFhistos.fill(HIST("hDCArToPv"), track.dcaXY());
      LFhistos.fill(HIST("hDCAzToPv"), track.dcaZ());
      LFhistos.fill(HIST("rawpT"), track.pt());
      LFhistos.fill(HIST("hIsPrim"), track.isPrimaryTrack());
      LFhistos.fill(HIST("hIsGood"), track.isGlobalTrackWoDCA());
      LFhistos.fill(HIST("hIsPrimCont"), track.isPVContributor());
      LFhistos.fill(HIST("hFindableTPCClusters"), track.tpcNClsFindable());
      LFhistos.fill(HIST("hFindableTPCRows"), track.tpcNClsCrossedRows());
      LFhistos.fill(HIST("hClustersVsRows"), track.tpcCrossedRowsOverFindableCls());
      LFhistos.fill(HIST("hTPCChi2"), track.tpcChi2NCl());
      LFhistos.fill(HIST("hITSChi2"), track.itsChi2NCl());

      if (!trackSelection(track))
        continue;

      LFhistos.fill(HIST("etaHistogram"), track.eta());
      LFhistos.fill(HIST("phiHistogram"), track.phi());
    } // LFTrack Loop
  };
  PROCESS_SWITCH(phiInJets, processLFTracks, "Process tracks in resonance framework", true);

}; // end of main struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phiInJets>(cfgc)};
};
