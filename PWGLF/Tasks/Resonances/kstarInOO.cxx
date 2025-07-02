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

/// \file kstarInOO.cxx
/// \author Jimun Lee <jimun.lee@cern.ch>

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TF1.h"
#include <TLorentzVector.h>
#include <TVector2.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct kstarInOO {
  SliceCache cache;
  HistogramRegistry OOhistos{"OOhistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::string> cfgeventSelections{"cfgeventSelections", "sel8", "choose event selection"};
  Configurable<std::string> cfgtrackSelections{"cfgtrackSelections", "globalTracks", "set track selections"};

  Configurable<double> cfgtrkMinPt{"cfgtrkMinPt", 0.15, "set track min pT"};
  Configurable<double> cfgtrkMaxEta{"cfgtrkMaxEta", 0.9, "set track max Eta"};
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
  Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0, "V_z cut selection"};
  Configurable<int> cDebugLevel{"cDebugLevel", 0, "Resolution of Debug"};

  void init(o2::framework::InitContext&)
  {
    // HISTOGRAMS
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPhi{200, -1, +7, "#phi"};
    const AxisSpec PtAxis = {200, 0, 20.0};
    const AxisSpec PIDAxis = {120, -6, 6};

    OOhistos.add("nEvents", "nEvents", kTH1F, {{4, 0.0, 4.0}});
    OOhistos.add("h_rawpT", "h_rawpT", kTH1F, {{1000, 0.0, 10.0}});
    OOhistos.add("h_rawpT_Kaon", "h_rawpT_Kaon", kTH1F, {{1000, 0.0, 10.0}});
    OOhistos.add("h_eta", "h_eta", kTH1F, {axisEta});
    OOhistos.add("h_phi", "h_phi", kTH1F, {axisPhi});

    OOhistos.add("QA_nSigma_pion_TPC", "QA_nSigma_pion_TPC", {HistType::kTH2F, {PtAxis, PIDAxis}});
    OOhistos.add("QA_nSigma_pion_TOF", "QA_nSigma_pion_TOF", {HistType::kTH2F, {PtAxis, PIDAxis}});
    OOhistos.add("QA_pion_TPC_TOF", "QA_pion_TPC_TOF", {HistType::kTH2F, {PIDAxis, PIDAxis}});
    OOhistos.add("QA_nSigma_kaon_TPC", "QA_nSigma_kaon_TPC", {HistType::kTH2F, {PtAxis, PIDAxis}});
    OOhistos.add("QA_nSigma_kaon_TOF", "QA_nSigma_kaon_TOF", {HistType::kTH2F, {PtAxis, PIDAxis}});
    OOhistos.add("QA_kaon_TPC_TOF", "QA_kaon_TPC_TOF", {HistType::kTH2F, {PIDAxis, PIDAxis}});

  } // end of init

  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiMinus;

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs>; // , aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPi, aod::pidTOFPi>;

  //==================================
  // 0. Track quality cuts
  //==================================
  // for PID QA TrackType
  template <typename TrackType>
  bool trackSelection(const TrackType track)
  {

    if (track.pt() < cfgtrkMinPt)
      return false;

    if (std::abs(track.eta()) > cfgtrkMaxEta)
      return false;

    if (std::abs(track.dcaXY()) > cfgMaxDCArToPVcut)
      return false;

    if (std::abs(track.dcaZ()) > cfgMaxDCAzToPVcut)
      return false;

    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;

    if (track.tpcNClsFindable() < cfgnFindableTPCClusters)
      return false;

    if (track.tpcNClsCrossedRows() < cfgnTPCCrossedRows)
      return false;

    if (track.tpcCrossedRowsOverFindableCls() > cfgnRowsOverFindable)
      return false;

    if (track.tpcChi2NCl() > cfgnTPCChi2)
      return false;

    if (track.itsChi2NCl() > cfgnITSChi2)
      return false;

    if (cfgConnectedToPV && !track.isPVContributor())
      return false;

    return true;
  };

  //---------------------------------------
  // 1-2. Check whether it passes tpc&tof
  //---------------------------------------
  // Kaon
  template <typename T>
  bool trackPIDKaon(const T& candidate, bool QA = false)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    // TPC
    if (QA) {
      OOhistos.fill(HIST("QA_nSigma_kaon_TPC"), candidate.pt(), candidate.tpcNSigmaKa());
      OOhistos.fill(HIST("QA_nSigma_kaon_TOF"), candidate.pt(), candidate.tofNSigmaKa());
    }
    if (std::abs(candidate.tpcNSigmaKa()) < cfgnTPCPID)
      tpcPIDPassed = true;

    // TOF
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaKa()) < cfgnTOFPID)
        tofPIDPassed = true;
      else
        tofPIDPassed = true;
    }

    // TPC & TOF
    if (tpcPIDPassed && tofPIDPassed)
      return true;

    return false;
  }

  // Pion
  template <typename T>
  bool trackPIDPion(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaPi()) < cfgnTPCPID)
      tpcPIDPassed = true;

    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPi()) < cfgnTOFPID)
        tofPIDPassed = true;
      else
        tofPIDPassed = true;
    }

    if (tpcPIDPassed && tofPIDPassed)
      return true;

    return false;
  }

  //================================
  // 3. Basic PID QA (Pion, Kaon)
  //================================
  // template <typename TrackType>
  // void fillHistograms(TrackType const& dTracks1, TrackType const& dTracks2)
  // {
  // for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2)))
  //   {
  //   // Full index policy is needed to consider all possible combinations
  //   if (trk1.index() == trk2.index())
  //     continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

  //   //// Initialize variables
  //   // trk1: Pion, trk2: Kaon
  //   // apply the track cut
  //   if (!trackSelection(trk1) || !trackSelection(trk2))
  //     continue;

  //   auto isTrk1hasTOF = trk1.hasTOF();
  //   auto isTrk2hasTOF = trk2.hasTOF();
  //   auto trk1ptPi = trk1.pt();
  //   auto trk1NSigmaPiTPC = trk1.tpcNSigmaPi();
  //   auto trk1NSigmaPiTOF = (isTrk1hasTOF) ? trk1.tofNSigmaPi() : -999.;
  //   auto trk2ptKa = trk2.pt();
  //   auto trk2NSigmaKaTPC = trk2.tpcNSigmaKa();
  //   auto trk2NSigmaKaTOF = (isTrk2hasTOF) ? trk2.tofNSigmaKa() : -999.;

  //   if (!trackPIDPion(trk1) || !trackPIDKaon(trk2))
  //     continue;

  //   // PID QA Pion
  //   OOhistos.fill(HIST("QA_nSigma_pion_TPC"), trk1ptPi, trk1NSigmaPiTPC);
  //   OOhistos.fill(HIST("QA_nSigma_pion_TOF"), trk1ptPi, trk1NSigmaPiTOF);
  //   OOhistos.fill(HIST("QA_pion_TPC_TOF"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);

  //   // PID QA Kaon
  //   OOhistos.fill(HIST("QA_nSigma_kaon_TPC"), trk2ptKa, trk2NSigmaKaTPC);
  //   OOhistos.fill(HIST("QA_nSigma_kaon_TOF"), trk2ptKa, trk2NSigmaKaTOF);
  //   OOhistos.fill(HIST("QA_kaon_TPC_TOF"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
  //   }
  // }

  //=================================
  // 1. nEvents Selection
  //=================================
  int nprocessEvents = 0;
  void processEvents(EventCandidates::iterator const& collision, TrackCandidates const& tracks)
  {
    // 1. All events
    if (cDebugLevel > 0) {
      nprocessEvents++;
      if ((nprocessEvents + 1) % 10000 == 0) {
        std::cout << "Processed Events: " << nprocessEvents << std::endl;
      }
    }
    OOhistos.fill(HIST("nEvents"), 0.5);
    if (std::fabs(collision.posZ()) > cfgVtxCut)
      return;

    // 2. The events passed a condition
    bool INELgt0 = false;
    for (const auto& track : tracks) {
      if (std::fabs(track.eta()) < cfgtrkMaxEta) {
        INELgt0 = true;
        break;
      }
    }
    if (!INELgt0) // not INEL
      return;

    OOhistos.fill(HIST("nEvents"), 1.5);

    //=====================================
    // 2. Basic track QA ( pt, phi, eta )
    //=====================================
    for (auto& track : tracks) {
      // auto originalTrack = track_as<TrackCandidates>();

      if (!trackSelection(track))
        continue;

      OOhistos.fill(HIST("h_rawpT"), track.pt());
      OOhistos.fill(HIST("h_eta"), track.eta());
      OOhistos.fill(HIST("h_phi"), track.phi());

      if (!trackPIDKaon(track, true)) // Once it sets the value is true, but later, should be change to false
        continue;

      OOhistos.fill(HIST("h_rawpT_Kaon"), track.pt());
    }
  }
  PROCESS_SWITCH(kstarInOO, processEvents, "Jimun Code Go!", true);

  void processEventsDummy(EventCandidates::iterator const&, TrackCandidates const&)
  {
    return;
  }
  PROCESS_SWITCH(kstarInOO, processEventsDummy, "dummy", false);

}; // kstarInOO

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<kstarInOO>(cfgc)};
};
