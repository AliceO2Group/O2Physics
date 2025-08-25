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
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  //==================================
  //||
  //||         Selection
  //||
  //==================================

  // Event Selection
  Configurable<float> cfgEventVtxCut{"cfgEventVtxCut", 10.0, "V_z cut selection"};
  ConfigurableAxis cfgCentAxis{"cfgCentAxis", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};

  // Track Selection
  // General
  Configurable<double> cfgTrackMinPt{"cfgTrackMinPt", 0.15, "set track min pT"};
  Configurable<double> cfgTrackMaxEta{"cfgTrackMaxEta", 0.9, "set track max Eta"};
  Configurable<double> cfgTrackMaxDCArToPVcut{"cfgTrackMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfgTrackMaxDCAzToPVcut{"cfgTrackMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<bool> cfgTrackPrimaryTrack{"cfgTrackPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgTrackConnectedToPV{"cfgTrackConnectedToPV", true, "PV contributor track selection"};           // PV Contriuibutor
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", true, "Global track selection"};                                 // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgTrackGlobalWoDCATrack{"cfgTrackGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  // TPC
  Configurable<double> cfgTracknFindableTPCClusters{"cfgTrackFindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfgTracknTPCCrossedRows{"cfgTrackTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgTracknRowsOverFindable{"cfgTrackRowsOverFindable", 1.2, "nRowsOverFindable TPC CLusters"};
  Configurable<double> cfgTracknTPCChi2{"cfgTrackTPCChi2", 4.0, "nTPC Chi2 per Cluster"};

  // ITS
  Configurable<double> cfgTracknITSChi2{"cfgTrackITSChi2", 36.0, "nITS Chi2 per Cluster"};

  // PID
  Configurable<bool> cfgTrackTPCPID{"cfgTrackTPCPID", true, "Enables TPC PID"};
  Configurable<bool> cfgTrackTOFPID{"cfgTrackTOFPID", true, "Enables TOF PID"};
  Configurable<float> cfgTrackTPCPIDnSig{"cfgTrackTPCPIDnSig", 4.0, "nTPC PID sigma"};
  Configurable<float> cfgTrackTOFPIDnSig{"cfgTrackTOFPIDnSig", 4.0, "nTOF PID sigma"};
  Configurable<int> cDebugLevel{"cDebugLevel", 0, "Resolution of Debug"};

  // Mixing
  ConfigurableAxis cfgBinsMixMult{"cfgBinsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
  ConfigurableAxis cfgBinsMixVtx{"cfgBinsMixVtx", {VARIABLE_WIDTH, -10.0f, -5.f, 0.f, 5.f, 10.f}, "Mixing bins - z-vertex"};
  Configurable<int> cfgMixNMixedEvents{"cfgMixNMixedEvents", 10, "Number of mixed events per event"};

  // MCGen
  Configurable<bool> cfgForceGenReco{"cfgForceGenReco", false, "Only consider events which are reconstructed (neglect event-loss)"};
  // Configurable<bool> cfgForceBR{"cfgForceBR", false, "Only consider K*0->K(pm)pi(mp)"};
  // Configurable<bool> cfgForceKaonAcceptence{"cfgForceKaonAcceptence", false, "Only consider K*0's whose daughters decay inside acceptence (no signal loss)"};
  // Configurable<bool> cfgForcePionAcceptence{"cfgForcePionAcceptence", false, "Only consider K*0's whose daughters decay inside acceptence (no signal loss)"};

  // Pair
  Configurable<int> cfgMinvNBins{"cfgMinvNBins", 300, "Number of bins for Minv axis"};
  Configurable<float> cfgMinvMin{"cfgMinvMin", 0.60, "Minimum Minv value"};
  Configurable<float> cfgMinvMax{"cfgMinvMax", 1.20, "Maximum Minv value"};

  ConfigurableAxis binsDCAz{"binsDCAz", {40, -0.2, 0.2}, ""};
  ConfigurableAxis binsDCAxy{"binsDCAxy", {40, -0.2, 0.2}, ""};

  // Histogram
  Configurable<bool> cfgEventCutQA{"cfgEventCutsQA", false, "Enable Event QA Hists"};
  Configurable<bool> cfgTrackCutQA{"cfgTrackCutQA", false, "Enable Track QA Hists"};
  Configurable<bool> cfgDataHistos{"cfgDataHistos", false, "Enable Data Hists"};
  Configurable<bool> cfgMcHistos{"cfgMcHistos", false, "Enable MC Hists"};

  void init(o2::framework::InitContext&)
  {
    // HISTOGRAMS
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPhi{200, -1, +7, "#phi"};
    const AxisSpec ptAxis = {200, 0, 20.0};
    const AxisSpec pidAxis = {120, -6, 6};
    const AxisSpec minvAxis = {cfgMinvNBins, cfgMinvMin, cfgMinvMax};
    const AxisSpec axisDCAz{binsDCAz, "DCA_{z}"};
    const AxisSpec axisDCAxy{binsDCAxy, "DCA_{XY}"};

    if (cfgEventCutQA) {
      histos.add("hPosZ_BC", "hPosZ_Bc", kTH1F, {{300, -15.0, 15.0}});
      histos.add("hPosZ_AC", "hPosZ_AC", kTH1F, {{300, -15.0, 15.0}});
      histos.add("hcentFT0C_BC", "centFT0C_BC", kTH1F, {{110, 0.0, 110.0}});
      histos.add("hcentFT0C_AC", "centFT0C_AC", kTH1F, {{110, 0.0, 110.0}});
    }

    if (cfgTrackCutQA) {
      histos.add("hDCArToPv_BC", "DCArToPv_BC", kTH1F, {axisDCAxy});
      histos.add("hDCAzToPv_BC", "DCAzToPv_BC", kTH1F, {axisDCAz});
      histos.add("hIsPrim_BC", "hIsPrim_BC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hIsGood_BC", "hIsGood_BC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hIsPrimCont_BC", "hIsPrimCont_BC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hFindableTPCClusters_BC", "hFindableTPCClusters_BC", kTH1F, {{200, 0, 200}});
      histos.add("hFindableTPCRows_BC", "hFindableTPCRows_BC", kTH1F, {{200, 0, 200}});
      histos.add("hClustersVsRows_BC", "hClustersVsRows_BC", kTH1F, {{200, 0, 2}});
      histos.add("hTPCChi2_BC", "hTPCChi2_BC", kTH1F, {{200, 0, 100}});
      histos.add("QA_nSigma_pion_TPC_BC", "QA_nSigma_pion_TPC_BC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_nSigma_pion_TOF_BC", "QA_nSigma_pion_TOF_BC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_pion_TPC_TOF_BC", "QA_pion_TPC_TOF_BC", {HistType::kTH2F, {pidAxis, pidAxis}});
      histos.add("QA_nSigma_kaon_TPC_BC", "QA_nSigma_kaon_TPC_BC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_nSigma_kaon_TOF_BC", "QA_nSigma_kaon_TOF_BC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_kaon_TPC_TOF_BC", "QA_kaon_TPC_TOF_BC", {HistType::kTH2F, {pidAxis, pidAxis}});
      histos.add("QA_track_pT_BC", "QA_track_pT_BC", kTH1F, {{13, 0.0, 13.0}});

      histos.add("hDCArToPv_AC", "DCArToPv_AC", kTH1F, {axisDCAxy});
      histos.add("hDCAzToPv_AC", "DCAzToPv_AC", kTH1F, {axisDCAz});
      histos.add("hIsPrim_AC", "hIsPrim_AC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hIsGood_AC", "hIsGood_AC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hIsPrimCont_AC", "hIsPrimCont_AC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hFindableTPCClusters_AC", "hFindableTPCClusters_AC", kTH1F, {{200, 0, 200}});
      histos.add("hFindableTPCRows_AC", "hFindableTPCRows_AC", kTH1F, {{200, 0, 200}});
      histos.add("hClustersVsRows_AC", "hClustersVsRows_AC", kTH1F, {{200, 0, 2}});
      histos.add("hTPCChi2_AC", "hTPCChi2_AC", kTH1F, {{200, 0, 100}});
      histos.add("QA_nSigma_pion_TPC_AC", "QA_nSigma_pion_TPC_AC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_nSigma_pion_TOF_AC", "QA_nSigma_pion_TOF_AC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_pion_TPC_TOF_AC", "QA_pion_TPC_TOF_AC", {HistType::kTH2F, {pidAxis, pidAxis}});
      histos.add("QA_nSigma_kaon_TPC_AC", "QA_nSigma_kaon_TPC_AC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_nSigma_kaon_TOF_AC", "QA_nSigma_kaon_TOF_AC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_kaon_TPC_TOF_AC", "QA_kaon_TPC_TOF_AC", {HistType::kTH2F, {pidAxis, pidAxis}});
      histos.add("QA_track_pT_AC", "QA_track_pT_AC", kTH1F, {{13, 0.0, 13.0}});
    }

    if (cfgDataHistos) {
      histos.add("nEvents", "nEvents", kTH1F, {{4, 0.0, 4.0}});
      histos.add("hUSS", "hUSS", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hLSS", "hLSS", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hUSS_Mix", "hUSS_Mix", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hLSS_Mix", "hLSS_Mix", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    }

    if (cfgMcHistos) {
      histos.add("nEvents_MC", "nEvents_MC", kTH1F, {{4, 0.0, 4.0}});
      histos.add("nEvents_MC_True", "nEvents_MC_True", kTH1F, {{4, 0.0, 4.0}});

      histos.add("hMC_USS", "hMC_USS", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMC_LSS", "hMC_LSS", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMC_USS_Mix", "hMC_USS_Mix", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMC_LSS_Mix", "hMC_LSS_Mix", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMC_USS_True", "hMC_USS_True", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMC_kstar_True", "hMC_kstar_True", kTHnSparseF, {cfgCentAxis, ptAxis});
    }
  } // end of init

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Cs>; //, aod::CentFT0Ms, aod::CentFT0As
  using EventCandidatesTrue = aod::McCollisions;

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using TrackCandidatesMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection,
                                      aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>;

  // For Mixed Event
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  Partition<TrackCandidates> kaon = !cfgTrackTPCPID || (nabs(aod::pidtpc::tpcNSigmaKa) <= cfgTrackTPCPIDnSig);
  Partition<TrackCandidates> pion = !cfgTrackTPCPID || (nabs(aod::pidtpc::tpcNSigmaPi) <= cfgTrackTPCPIDnSig);

  Partition<TrackCandidatesMC> kaonMC = !cfgTrackTPCPID || (nabs(aod::pidtpc::tpcNSigmaKa) <= cfgTrackTPCPIDnSig);
  Partition<TrackCandidatesMC> pionMC = !cfgTrackTPCPID || (nabs(aod::pidtpc::tpcNSigmaPi) <= cfgTrackTPCPIDnSig);

  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiMinus;

  //==================================
  //||
  //||       Helper Templates
  //||
  //==================================
  template <typename EventType>
  bool eventSelection(const EventType event, const bool QA)
  {
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hPosZ_BC"), event.posZ());
      histos.fill(HIST("hcentFT0C_BC"), event.centFT0C());
    }
    if (!event.sel8())
      return false;
    if (std::abs(event.posZ()) > cfgEventVtxCut)
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

    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hPosZ_AC"), event.posZ());
      histos.fill(HIST("hcentFT0C_AC"), event.centFT0C());
    }
    return true;
  };

  template <typename TracksType>
  bool trackSelection(const TracksType track, bool QA)
  {
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("hDCArToPv_BC"), track.dcaXY());
      histos.fill(HIST("hDCAzToPv_BC"), track.dcaZ());
      histos.fill(HIST("hIsPrim_BC"), track.isPrimaryTrack());
      histos.fill(HIST("hIsGood_BC"), track.isGlobalTrackWoDCA());
      histos.fill(HIST("hIsPrimCont_BC"), track.isPVContributor());
      histos.fill(HIST("hFindableTPCClusters_BC"), track.tpcNClsFindable());
      histos.fill(HIST("hFindableTPCRows_BC"), track.tpcNClsCrossedRows());
      histos.fill(HIST("hClustersVsRows_BC"), track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("hTPCChi2_BC"), track.tpcChi2NCl());
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (track.pt() < cfgTrackMinPt)
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (std::abs(track.eta()) > cfgTrackMaxEta)
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (!cfgGlobalTrack && !track.isGlobalTrack())
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (std::abs(track.dcaXY()) > cfgTrackMaxDCArToPVcut)
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (std::abs(track.dcaZ()) > cfgTrackMaxDCAzToPVcut)
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (cfgTrackPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (cfgTrackGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (track.tpcNClsFindable() < cfgTracknFindableTPCClusters)
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (track.tpcNClsCrossedRows() < cfgTracknTPCCrossedRows)
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (track.tpcCrossedRowsOverFindableCls() > cfgTracknRowsOverFindable)
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (track.tpcChi2NCl() > cfgTracknTPCChi2)
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (track.itsChi2NCl() > cfgTracknITSChi2)
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (cfgTrackConnectedToPV && !track.isPVContributor())
      return false;
    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("QA_track_pT_BC"), track.pt());
    }

    if (cfgTrackCutQA && QA) {
      histos.fill(HIST("hDCArToPv_AC"), track.dcaXY());
      histos.fill(HIST("hDCAzToPv_AC"), track.dcaZ());
      histos.fill(HIST("hIsPrim_AC"), track.isPrimaryTrack());
      histos.fill(HIST("hIsGood_AC"), track.isGlobalTrackWoDCA());
      histos.fill(HIST("hIsPrimCont_AC"), track.isPVContributor());
      histos.fill(HIST("hFindableTPCClusters_AC"), track.tpcNClsFindable());
      histos.fill(HIST("hFindableTPCRows_AC"), track.tpcNClsCrossedRows());
      histos.fill(HIST("hClustersVsRows_AC"), track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("hTPCChi2_AC"), track.tpcChi2NCl());
      histos.fill(HIST("QA_track_pT_AC"), track.pt());
    }
    return true;
  };

  template <typename TrackPID>
  bool trackPIDKaon(const TrackPID& candidate, const bool QA)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (cfgTrackCutQA && QA) {
      // kaon
      histos.fill(HIST("QA_nSigma_kaon_TPC_BC"), candidate.pt(), candidate.tpcNSigmaKa());
      histos.fill(HIST("QA_nSigma_kaon_TOF_BC"), candidate.pt(), candidate.tofNSigmaKa());
      histos.fill(HIST("QA_kaon_TPC_TOF_BC"), candidate.tpcNSigmaKa(), candidate.tofNSigmaKa());
    }
    // TPC
    if (std::abs(candidate.tpcNSigmaKa()) < cfgTrackTPCPIDnSig)
      tpcPIDPassed = true;

    // TOF
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaKa()) < cfgTrackTOFPIDnSig) {
        tofPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
    }

    // TPC & TOF
    if (tpcPIDPassed && tofPIDPassed) {
      if (cfgTrackCutQA && QA) {
        // kaon
        histos.fill(HIST("QA_nSigma_kaon_TPC_AC"), candidate.pt(), candidate.tpcNSigmaKa());
        histos.fill(HIST("QA_nSigma_kaon_TOF_AC"), candidate.pt(), candidate.tofNSigmaKa());
        histos.fill(HIST("QA_kaon_TPC_TOF_AC"), candidate.tpcNSigmaKa(), candidate.tofNSigmaKa());
      }
      return true;
    }
    return false;
  }

  template <typename TrackPID>
  bool trackPIDPion(const TrackPID& candidate, const bool QA)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (cfgTrackCutQA && QA) {
      // pion
      histos.fill(HIST("QA_nSigma_pion_TPC_BC"), candidate.pt(), candidate.tpcNSigmaPi());
      histos.fill(HIST("QA_nSigma_pion_TOF_BC"), candidate.pt(), candidate.tofNSigmaPi());
      histos.fill(HIST("QA_pion_TPC_TOF_BC"), candidate.tpcNSigmaPi(), candidate.tofNSigmaPi());
    }

    // TPC
    if (std::abs(candidate.tpcNSigmaPi()) < cfgTrackTPCPIDnSig)
      tpcPIDPassed = true;
    // TOF
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPi()) < cfgTrackTOFPIDnSig) {
        tofPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
    }

    // TPC & TOF
    if (tpcPIDPassed && tofPIDPassed) {
      if (cfgTrackCutQA && QA) {
        // pion
        histos.fill(HIST("QA_nSigma_pion_TPC_AC"), candidate.pt(), candidate.tpcNSigmaPi());
        histos.fill(HIST("QA_nSigma_pion_TOF_AC"), candidate.pt(), candidate.tofNSigmaPi());
        histos.fill(HIST("QA_pion_TPC_TOF_AC"), candidate.tpcNSigmaPi(), candidate.tofNSigmaPi());
      }
      return true;
    }
    return false;
  }

  template <typename CollisionType, typename TracksType>
  void TrackSlicing(const CollisionType& collision1, const TracksType&, const CollisionType& collision2, const TracksType&, const bool QA, const bool IsMix)
  {
    auto tracks1 = kaon->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
    auto tracks2 = pion->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
    auto centrality = collision1.centFT0C();

    for (const auto& [trk1, trk2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
      auto [KstarPt, Minv] = minvReconstruction(trk1, trk2, QA);
      if (Minv < 0)
        continue;

      double conjugate = trk1.sign() * trk2.sign();
      if (cfgDataHistos) {
        if (!IsMix) {
          if (conjugate < 0) {
            histos.fill(HIST("hUSS"), centrality, KstarPt, Minv);
          } else if (conjugate > 0) {
            histos.fill(HIST("hLSS"), centrality, KstarPt, Minv);
          }
        } else {
          if (conjugate < 0) {
            histos.fill(HIST("hUSS_Mix"), centrality, KstarPt, Minv);
          } else if (conjugate > 0) {
            histos.fill(HIST("hLSS_Mix"), centrality, KstarPt, Minv);
          }
        }
      } // cfgDataHistos
    } // for
  } // TrackSlicing

  template <typename CollisionType, typename TracksType>
  void TrackSlicingMC(const CollisionType& collision1, const TracksType&, const CollisionType& collision2, const TracksType&, const bool QA, const bool IsMix)
  {
    auto tracks1 = kaonMC->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
    auto tracks2 = pionMC->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
    auto centrality = collision1.centFT0C();

    for (const auto& [trk1, trk2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
      auto [KstarPt, Minv] = minvReconstruction(trk1, trk2, QA);
      if (Minv < 0)
        continue;

      double conjugate = trk1.sign() * trk2.sign();
      if (cfgMcHistos) {
        if (!IsMix) {
          if (conjugate < 0) {
            histos.fill(HIST("hMC_USS"), centrality, KstarPt, Minv);
          } else if (conjugate > 0) {
            histos.fill(HIST("hMC_LSS"), centrality, KstarPt, Minv);
          }
        } else {
          if (conjugate < 0) {
            histos.fill(HIST("hMC_USS_Mix"), centrality, KstarPt, Minv);
          } else if (conjugate > 0) {
            histos.fill(HIST("hMC_LSS_Mix"), centrality, KstarPt, Minv);
          }
        }
      } // cfgMcHistos

      //======================
      // Gen MC
      if (!trk1.has_mcParticle() || !trk2.has_mcParticle())
        continue;

      auto particle1 = trk1.mcParticle();
      auto particle2 = trk2.mcParticle();
      if (std::fabs(particle1.pdgCode()) != 321)
        continue; // Not Kaon
      if (std::fabs(particle2.pdgCode()) != 211)
        continue; // Not Pion

      if (!particle1.has_mothers()) {
        continue;
      }
      if (!particle2.has_mothers())
        continue;

      std::vector<int> mothers1{};
      std::vector<int> mothers1PDG{};
      for (auto& particle1_mom : particle1.template mothers_as<aod::McParticles>()) {
        mothers1.push_back(particle1_mom.globalIndex());
        mothers1PDG.push_back(particle1_mom.pdgCode());
      }

      std::vector<int> mothers2{};
      std::vector<int> mothers2PDG{};
      for (auto& particle2_mom : particle2.template mothers_as<aod::McParticles>()) {
        mothers2.push_back(particle2_mom.globalIndex());
        mothers2PDG.push_back(particle2_mom.pdgCode());
      }

      if (mothers1PDG[0] != 313)
        continue; // mother not K*0
      if (mothers2PDG[0] != 313)
        continue; // mothers not K*0

      if (mothers1[0] != mothers2[0])
        continue; // Kaon and pion not from the same K*0

      if (cfgMcHistos) {
        histos.fill(HIST("hMC_USS_True"), centrality, KstarPt, Minv);
      }
    } // for
  } // TrackSlicingMC

  template <typename TracksType>
  std::pair<double, double> minvReconstruction(const TracksType& trk1, const TracksType& trk2, const bool QA)
  {
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;

    if (!trackSelection(trk1, QA) || !trackSelection(trk2, QA))
      return {-1.0, -1.0};
    if (!trackPIDKaon(trk1, QA) || !trackPIDPion(trk2, QA))
      return {-1.0, -1.0};

    if (trk1.globalIndex() == trk2.globalIndex()) {
      return {-1.0, -1.0}; // For Kstar, we need to run (0,1), (1,0) pairs as well. but same id pairs are not need.
    }

    lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
    lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massPi);
    lResonance = lDecayDaughter1 + lDecayDaughter2;

    if (std::abs(lResonance.Eta()) > cfgTrackMaxEta)
      return {-1.0, -1.0};

    return {lResonance.Pt(), lResonance.M()};
  }

  //=======================================================
  //|
  //|                  DATA STUFF (SE)
  //|
  //=======================================================

  int nEvents = 0;
  void processDataSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks)
  {
    if (cDebugLevel > 0) {
      nEvents++;
      if ((nEvents + 1) % 10000 == 0) {
        std::cout << "Processed Data Events: " << nEvents << std::endl;
      }
    }

    auto goodEv = eventSelection(collision, true);
    if (cfgDataHistos) {
      histos.fill(HIST("nEvents"), 0.5);
    }
    if (!goodEv)
      return;

    bool INELgt0 = false;
    for (const auto& track : tracks) {
      if (!trackSelection(track, true))
        continue;
      if (std::fabs(track.eta()) < cfgTrackMaxEta) {
        INELgt0 = true;
      }
    }
    if (!INELgt0)
      return;

    if (cfgDataHistos) {
      histos.fill(HIST("nEvents"), 1.5);
    }
    TrackSlicing(collision, tracks, collision, tracks, true, false);

  } // processSameEvents
  PROCESS_SWITCH(kstarInOO, processDataSameEvent, "process Data Same Event", false);

  //=======================================================
  //|
  //|                  DATA STUFF (ME)
  //|
  //=======================================================

  int nEventsMix = 0;
  void processDataMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningType colBinning{{cfgBinsMixVtx, cfgBinsMixMult}, true}; // true is for 'ignore overflows' (true by default)
    SameKindPair<EventCandidates, TrackCandidates, BinningType> pairs{colBinning, cfgMixNMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (cDebugLevel > 0) {
        nEventsMix++;
        if ((nEventsMix + 1) % 10000 == 0) {
          std::cout << "Processed DATA Mixed Events : " << nEventsMix << std::endl;
        }
      }
      auto goodEv1 = eventSelection(collision1, false);
      auto goodEv2 = eventSelection(collision2, false);

      if (!goodEv1 || !goodEv2)
        continue;

      TrackSlicing(collision1, tracks1, collision2, tracks2, false, true);
    }
  }
  PROCESS_SWITCH(kstarInOO, processDataMixedEvent, "process DATA Mixed Event", false);

  //=======================================================
  //|
  //|                  MC STUFF (SE)
  //|
  //=======================================================

  int nEventsMC = 0;
  void processMCSameEvent(EventCandidates::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&)
  {
    if (cDebugLevel > 0) {
      nEventsMC++;
      if ((nEventsMC + 1) % 10000 == 0) {
        double histmem = histos.getSize();
        std::cout << histmem << std::endl;
        std::cout << "process MC Same Event : " << nEventsMC << std::endl;
      }
    }

    auto goodEv = eventSelection(collision, true);
    if (cfgMcHistos) {
      histos.fill(HIST("nEvents_MC"), 0.5);
    }
    if (!goodEv)
      return;

    bool INELgt0 = false;
    for (const auto& track : tracks) {
      if (std::fabs(track.eta()) < cfgTrackMaxEta) {
        INELgt0 = true;
      }
    }
    if (!INELgt0)
      return;

    if (cfgMcHistos) {
      histos.fill(HIST("nEvents_MC"), 1.5);
    }
    TrackSlicingMC(collision, tracks, collision, tracks, true, false);

  } // processSameEvents_MC
  PROCESS_SWITCH(kstarInOO, processMCSameEvent, "process MC Same Event", true);

  //=======================================================
  //|
  //|                  MC STUFF (ME)
  //|
  //=======================================================

  int nEventsMCMix = 0;
  void processMCMixedEvent(EventCandidates const& collisions, TrackCandidatesMC const& tracks, aod::McParticles const&)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningType colBinning{{cfgBinsMixVtx, cfgBinsMixMult}, true}; // true is for 'ignore overflows' (true by default)
    SameKindPair<EventCandidates, TrackCandidatesMC, BinningType> pairs{colBinning, cfgMixNMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (cDebugLevel > 0) {
        nEventsMCMix++;
        if ((nEventsMCMix + 1) % 10000 == 0) {
          std::cout << "Processed MC Mixed Events : " << nEventsMCMix << std::endl;
        }
      }
      auto goodEv1 = eventSelection(collision1, false);
      auto goodEv2 = eventSelection(collision2, false);
      if (!goodEv1 || !goodEv2)
        continue;
      TrackSlicingMC(collision1, tracks1, collision2, tracks2, false, true);
    } // mixing
  } // processMixedEvent_MC
  PROCESS_SWITCH(kstarInOO, processMCMixedEvent, "process MC Mixed Event", false);

  //=======================================================
  //|
  //|             GENERATED MC STUFF (TRUE)
  //|
  //=======================================================

  int nEventsTrue = 0;
  void processMCTrue(EventCandidatesTrue::iterator const& collision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, EventCandidates>> const& recocolls, aod::McParticles const& particles)
  {
    if (cDebugLevel > 0) {
      ++nEventsTrue;
    }

    if (fabs(collision.posZ()) > cfgEventVtxCut)
      return;
    if (recocolls.size() <= 0) { // not reconstructed
      if (cfgForceGenReco) {
        return;
      }
    }

    double centrality = -1;
    for (auto& recocoll : recocolls) {
      centrality = recocoll.centFT0C();
      auto goodEv = eventSelection(recocoll, false);

      if (cfgMcHistos) {
        histos.fill(HIST("nEvents_MC_True"), 0.5);
      }
      if (!goodEv)
        continue;
    } // for

    for (auto& particle : particles) {
      if (particle.pdgCode() != 313)
        continue; // Not K*0
      if (std::fabs(particle.eta()) > cfgTrackMaxEta)
        continue;

      if (cfgMcHistos) {
        histos.fill(HIST("hMC_kstar_True"), centrality, particle.pt());
      }

    } // loop over particles

  } // processMCTrue
  PROCESS_SWITCH(kstarInOO, processMCTrue, "process MC True", false);

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
