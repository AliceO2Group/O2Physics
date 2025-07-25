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
/// \brief the pT spectra of k*0(892) resonance analysis in OO collisions
/// \author Jimun Lee <jimun.lee@cern.ch>

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"
#include <Framework/ASoAHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include "TRandom.h"
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMathBase.h>
#include <TVector2.h>

#include <RtypesCore.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <stdlib.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct kstarInOO {
  SliceCache cache;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry OOhistos{"OOhistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  //==================================
  //||
  //||         Selection
  //||
  //==================================

  // Event Selection
  Configurable<std::string> cfg_Event_Selections{"cfg_Event_Selections", "sel8", "choose event selection"};
  Configurable<float> cfg_Event_VtxCut{"cfg_Event_VtxCut", 10.0, "V_z cut selection"};

  ConfigurableAxis cfg_CentAxis{"cfg_CentAxis", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};

  // Track Selection
  // General
  Configurable<std::string> cfg_Track_Selections{"cfg_Track_Selections", "globalTracks", "set track selections"};
  Configurable<double> cfg_Track_MinPt{"cfg_Track_MinPt", 0.15, "set track min pT"};
  Configurable<double> cfg_Track_MaxEta{"cfg_Track_MaxEta", 0.9, "set track max Eta"};
  Configurable<double> cfg_Track_MaxDCArToPVcut{"cfg_Track_MaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfg_Track_MaxDCAzToPVcut{"cfg_Track_MaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<bool> cfg_Track_PrimaryTrack{"cfg_Track_PrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfg_Track_ConnectedToPV{"cfg_Track_ConnectedToPV", true, "PV contributor track selection"};           // PV Contriuibutor
  Configurable<bool> cfg_Track_GlobalWoDCATrack{"cfg_Track_GlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  // TPC
  Configurable<double> cfg_Track_nFindableTPCClusters{"cfg_Track_FindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfg_Track_nTPCCrossedRows{"cfg_Track_TPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfg_Track_nRowsOverFindable{"cfg_Track_RowsOverFindable", 1.2, "nRowsOverFindable TPC CLusters"};
  Configurable<double> cfg_Track_nTPCChi2{"cfg_Track_TPCChi2", 4.0, "nTPC Chi2 per Cluster"};

  // ITS
  Configurable<double> cfg_Track_nITSChi2{"cfg_Track_ITSChi2", 36.0, "nITS Chi2 per Cluster"};

  // PID
  Configurable<bool> cfg_Track_TPCPID{"cfg_Track_TPCPID", true, "Enables TPC PID"};
  Configurable<bool> cfg_Track_TOFPID{"cfg_Track_TOFPID", true, "Enables TOF PID"};
  Configurable<float> cfg_Track_TPCPID_nSig{"cfg_Track_TPCPID_nSig", 4.0, "nTPC PID sigma"};
  Configurable<float> cfg_Track_TOFPID_nSig{"cfg_Track_TOFPID_nSig", 4.0, "nTOF PID sigma"};
  Configurable<int> cDebugLevel{"cDebugLevel", 0, "Resolution of Debug"};

  // Mixing
  ConfigurableAxis cfg_bins_MixVtx{"cfg_bins_MixVtx", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfg_bins_MixMult{"cfg_bins_MixMult", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f}, "Mixing bins - z-vertex"};
  Configurable<int> cfg_Mix_NMixedEvents{"cfg_Mix_NMixedEvents", 5, "Number of mixed events per event"};

  // Pair
  Configurable<int> cfg_MinvNBins{"cfg_MinvNBins", 300, "Number of bins for Minv axis"};
  Configurable<float> cfg_MinvMin{"cfg_MinvMin", 0.60, "Minimum Minv value"};
  Configurable<float> cfg_MinvMax{"cfg_MinvMax", 1.20, "Maximum Minv value"};

  // Histogram
  Configurable<bool> cfg_Track_CutQA{"cfg_Track_CutQA", false, "Enable Track QA Hists"};

  // std::vector<int> eventSelectionBits;

  void init(o2::framework::InitContext&)
  {
    // HISTOGRAMS
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPhi{200, -1, +7, "#phi"};
    const AxisSpec PtAxis = {200, 0, 20.0};
    const AxisSpec PIDAxis = {120, -6, 6};
    const AxisSpec MinvAxis = {cfg_MinvNBins, cfg_MinvMin, cfg_MinvMax};

    if (cfg_Track_CutQA) {
      OOhistos.add("h_rawpT", "h_rawpT", kTH1F, {{1000, 0.0, 10.0}});
      OOhistos.add("h_rawpT_Kaon", "h_rawpT_Kaon", kTH1F, {{1000, 0.0, 10.0}});
      OOhistos.add("h_rawpT_Pion", "h_rawpT_Pion", kTH1F, {{1000, 0.0, 10.0}});
      OOhistos.add("h_eta", "h_eta", kTH1F, {axisEta});
      OOhistos.add("h_phi", "h_phi", kTH1F, {axisPhi});

      OOhistos.add("QA_nSigma_pion_TPC", "QA_nSigma_pion_TPC", {HistType::kTH2F, {PtAxis, PIDAxis}});
      OOhistos.add("QA_nSigma_pion_TOF", "QA_nSigma_pion_TOF", {HistType::kTH2F, {PtAxis, PIDAxis}});
      OOhistos.add("QA_pion_TPC_TOF", "QA_pion_TPC_TOF", {HistType::kTH2F, {PIDAxis, PIDAxis}});
      OOhistos.add("QA_nSigma_kaon_TPC", "QA_nSigma_kaon_TPC", {HistType::kTH2F, {PtAxis, PIDAxis}});
      OOhistos.add("QA_nSigma_kaon_TOF", "QA_nSigma_kaon_TOF", {HistType::kTH2F, {PtAxis, PIDAxis}});
      OOhistos.add("QA_kaon_TPC_TOF", "QA_kaon_TPC_TOF", {HistType::kTH2F, {PIDAxis, PIDAxis}});
    }

    // MC histos
    OOhistos.add("hMC_USS", "hMC_USS", kTHnSparseF, {cfg_CentAxis, PtAxis, MinvAxis});
    OOhistos.add("hMC_LSS", "hMC_LSS", kTHnSparseF, {cfg_CentAxis, PtAxis, MinvAxis});
    OOhistos.add("hMC_USS_Mix", "hMC_USS_Mix", kTHnSparseF, {cfg_CentAxis, PtAxis, MinvAxis});
    OOhistos.add("hMC_LSS_Mix", "hMC_LSS_Mix", kTHnSparseF, {cfg_CentAxis, PtAxis, MinvAxis});

    // OOhistos.add("hMC_pt_Pion", "hMC_pt_Pion", kTH1F, {PtAxis});
    // OOhistos.add("hMC_pt_Kaon", "hMC_pt_Kaon", kTH1F, {PtAxis});
    // OOhistos.add("hMC_pt_Proton", "hMC_pt_Proton", kTH1F, {PtAxis});

    // Event Histograms
    OOhistos.add("nEvents_MC", "nEvents_MC", kTH1F, {{4, 0.0, 4.0}});
    OOhistos.add("nEvents_MC_Mix", "nEvents_MC_Mix", kTH1F, {{4, 0.0, 4.0}});

  } // end of init

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Cs>; //, aod::CentFT0Ms, aod::CentFT0As
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using TrackCandidates_MC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection,
                                       aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>;

  // For Mixed Event
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  Partition<TrackCandidates_MC> Kaon_MC = (!cfg_Track_TPCPID || (nabs(aod::pidtpc::tpcNSigmaKa) <= cfg_Track_TPCPID_nSig));
  Partition<TrackCandidates_MC> Pion_MC = (!cfg_Track_TPCPID || (nabs(aod::pidtpc::tpcNSigmaPi) <= cfg_Track_TPCPID_nSig));

  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiMinus;

  //==================================
  //||
  //||       Helper Templates
  //||
  //==================================
  template <typename EventType>
  bool eventSelection(const EventType event)
  {
    if (!event.sel8())
      return false;

    if (!event.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if (!event.selection_bit(aod::evsel::kNoSameBunchPileup))
      return false;

    if (!event.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;
    if (!event.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;

    if (!event.selection_bit(aod::evsel::kNoCollInTimeRangeStandard))
      return false;

    return true;
  };

  template <typename TracksType>
  bool trackSelection(const TracksType track)
  {
    if (track.pt() < cfg_Track_MinPt)
      return false;

    if (std::abs(track.eta()) > cfg_Track_MaxEta)
      return false;

    if (std::abs(track.dcaXY()) > cfg_Track_MaxDCArToPVcut)
      return false;

    if (std::abs(track.dcaZ()) > cfg_Track_MaxDCAzToPVcut)
      return false;

    if (cfg_Track_PrimaryTrack && !track.isPrimaryTrack())
      return false;

    if (cfg_Track_GlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;

    if (track.tpcNClsFindable() < cfg_Track_nFindableTPCClusters)
      return false;

    if (track.tpcNClsCrossedRows() < cfg_Track_nTPCCrossedRows)
      return false;

    if (track.tpcCrossedRowsOverFindableCls() > cfg_Track_nRowsOverFindable)
      return false;

    if (track.tpcChi2NCl() > cfg_Track_nTPCChi2)
      return false;

    if (track.itsChi2NCl() > cfg_Track_nITSChi2)
      return false;

    if (cfg_Track_ConnectedToPV && !track.isPVContributor())
      return false;

    return true;
  };

  template <typename TrackPID>
  bool trackPIDKaon(const TrackPID& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    // TPC
    if (cfg_Track_CutQA) {
      OOhistos.fill(HIST("QA_nSigma_kaon_TPC"), candidate.pt(), candidate.tpcNSigmaKa());
      OOhistos.fill(HIST("QA_nSigma_kaon_TOF"), candidate.pt(), candidate.tofNSigmaKa());
      OOhistos.fill(HIST("QA_kaon_TPC_TOF"), candidate.tpcNSigmaKa(), candidate.tofNSigmaKa());
    }
    if (std::abs(candidate.tpcNSigmaKa()) < cfg_Track_TPCPID_nSig)
      tpcPIDPassed = true;

    // TOF
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaKa()) < cfg_Track_TOFPID_nSig)
        tofPIDPassed = true;
      else
        tofPIDPassed = true;
    }

    // TPC & TOF
    if (tpcPIDPassed && tofPIDPassed)
      return true;

    return false;
  }

  template <typename TrackPID>
  bool trackPIDPion(const TrackPID& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    // TPC
    if (cfg_Track_CutQA) {
      OOhistos.fill(HIST("QA_nSigma_pion_TPC"), candidate.pt(), candidate.tpcNSigmaPi());
      OOhistos.fill(HIST("QA_nSigma_pion_TOF"), candidate.pt(), candidate.tofNSigmaPi());
      OOhistos.fill(HIST("QA_pion_TPC_TOF"), candidate.tpcNSigmaPi(), candidate.tofNSigmaPi());
    }

    if (std::abs(candidate.tpcNSigmaPi()) < cfg_Track_TPCPID_nSig)
      tpcPIDPassed = true;

    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPi()) < cfg_Track_TOFPID_nSig)
        tofPIDPassed = true;
      else
        tofPIDPassed = true;
    }

    // TPC & TOF
    if (tpcPIDPassed && tofPIDPassed)
      return true;

    return false;
  }

  template <typename CollisionType, typename TracksType>
  void TrackSlicing_MC(const CollisionType& collision1, const TracksType&, const CollisionType& collision2, const TracksType&, const bool IsMix)
  {
    auto tracks1 = Kaon_MC->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
    auto tracks2 = Pion_MC->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
    auto centrality = collision1.centFT0C();

    for (auto& [trk1, trk2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
      auto [KstarPt, Minv] = minvReconstruction(trk1, trk2);
      if (Minv < 0)
        continue;

      double conjugate = trk1.sign() * trk2.sign();
      if (!IsMix) {
        if (conjugate < 0) {
          OOhistos.fill(HIST("hMC_USS"), centrality, KstarPt, Minv);
        } else if (conjugate > 0) {
          OOhistos.fill(HIST("hMC_LSS"), centrality, KstarPt, Minv);
        }
      } else {
        if (conjugate < 0) {
          OOhistos.fill(HIST("hMC_USS_Mix"), centrality, KstarPt, Minv);
        } else if (conjugate > 0) {
          OOhistos.fill(HIST("hMC_LSS_Mix"), centrality, KstarPt, Minv);
        }
      }
    }
  }

  template <typename TracksType>
  std::pair<double, double> minvReconstruction(const TracksType& trk1, const TracksType& trk2)
  {
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;

    if (!trackSelection(trk1) || !trackSelection(trk2))
      return {-1.0, -1.0};
    if (!trackPIDKaon(trk1) || !trackPIDPion(trk2))
      return {-1.0, -1.0};

    if (trk1.globalIndex() == trk2.globalIndex())
      return {-1.0, -1.0}; // For Kstar, we need to run (0,1), (1,0) pairs as well. but same id pairs are not neede.

    lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
    lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massPi);
    lResonance = lDecayDaughter1 + lDecayDaughter2;

    if (std::abs(lResonance.Eta()) > cfg_Track_MaxEta)
      return {-1.0, -1.0};

    return {lResonance.Pt(), lResonance.M()};
  }

  //=======================================================
  //|
  //|                  MC STUFF (SE)
  //|
  //=======================================================

  int nEvents_MC = 0;
  void processSameEvent_MC(EventCandidates::iterator const& collision, TrackCandidates_MC const& tracks, aod::McParticles const&)
  {
    if (cDebugLevel > 0) {
      nEvents_MC++;
      if ((nEvents_MC + 1) % 10000 == 0) {
        double histmem = OOhistos.getSize();
        std::cout << histmem << std::endl;
        std::cout << "process_SameEvent_MC: " << nEvents_MC << std::endl;
      }
    }

    auto goodEv = eventSelection(collision);
    OOhistos.fill(HIST("nEvents_MC"), 0.5);
    if (!goodEv)
      return;

    if (std::fabs(collision.posZ()) > cfg_Event_VtxCut)
      return;

    bool INELgt0 = false;
    for (const auto& track : tracks) {
      if (std::fabs(track.eta()) < cfg_Track_MaxEta) {
        INELgt0 = true;
        break;
      }
    }
    if (!INELgt0)
      return;

    OOhistos.fill(HIST("nEvents_MC"), 1.5);

    TrackSlicing_MC(collision, tracks, collision, tracks, false);

  } // processSameEvents_MC
  PROCESS_SWITCH(kstarInOO, processSameEvent_MC, "process Same Event MC", true);

  //=======================================================
  //|
  //|                  MC STUFF (ME)
  //|
  //=======================================================

  int nEvents_MC_Mix = 0;
  void processMixedEvent_MC(EventCandidates const& collisions, TrackCandidates_MC const& tracks, aod::McParticles const&)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningType colBinning{{cfg_bins_MixVtx, cfg_bins_MixMult}, true}; // true is for 'ignore overflows' (true by default)
    SameKindPair<EventCandidates, TrackCandidates_MC, BinningType> pairs{colBinning, cfg_Mix_NMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (cDebugLevel > 0) {
        nEvents_MC_Mix++;
        if ((nEvents_MC_Mix + 1) % 10000 == 0) {
          std::cout << "Processed Mixed Events: " << nEvents_MC_Mix << std::endl;
        }
      }

      auto goodEv1 = eventSelection(collision1);
      auto goodEv2 = eventSelection(collision2);
      OOhistos.fill(HIST("nEvents_MC_Mix"), 0.5);

      if (!goodEv1 || !goodEv2)
        continue;

      OOhistos.fill(HIST("nEvents_MC_Mix"), 1.5);

      TrackSlicing_MC(collision1, tracks1, collision2, tracks2, true);
    } // mixing
  } // processMixedEvent_MC
  PROCESS_SWITCH(kstarInOO, processMixedEvent_MC, "process Mixed Event MC", false);

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
