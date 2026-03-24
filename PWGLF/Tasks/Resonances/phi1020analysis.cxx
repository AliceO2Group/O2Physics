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

/// \file phi1020analysis.cxx
/// \brief Reconstruction of Phi yield through track-track Minv correlations for resonance OO analysis
///
///
/// \author Adrian Fereydon Nassirpour <adrian.fereydon.nassirpour@cern.ch
/// \author Hirak Kumar Koley <hirak.kumar.koley@cern.ch

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

struct phi1020analysis {

  SliceCache cache;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Event configurables
  Configurable<std::string> cfg_Event_Sel{"cfg_Event_Sel", "sel8", "choose event selection"};
  Configurable<float> cfg_Event_VtxCut{"cfg_Event_VtxCut", 10.0, "V_z cut selection"};
  Configurable<bool> cfg_Event_Timeframe{"cfg_Event_Timeframe", true, "Timeframe border cut"};
  Configurable<bool> cfg_Event_Timerange{"cfg_Event_Timerange", true, "Timerange border cut"};
  Configurable<bool> cfg_Event_Centrality{"cfg_Event_Centrality", true, "Centrality cut"};
  Configurable<float> cfg_Event_CentralityMax{"cfg_Event_CentralityMax", 100, "CentralityMax cut"};
  Configurable<bool> cfg_Event_Pileup{"cfg_Event_Pileup", true, "Pileup border cut"};
  Configurable<bool> cfg_Event_OccupancyCut{"cfg_Event_OccupancyCut", true, "Occupancy border cut"};
  Configurable<float> cfg_Event_MaxOccupancy{"cfg_Event_MaxOccupancy", 1, "Max TPC Occupancy"};
  Configurable<int> centEstimator{"centEstimator", 0, "Select centrality estimator: 0 - FT0M, 1 - FT0A, 2 - FT0C"};

  // Track configurables
  Configurable<std::string> cfg_Track_Sel{"cfg_Track_Sel", "globalTracks", "set track selections"};
  Configurable<float> cfg_Track_MinPt{"cfg_Track_MinPt", 0.15, "set track min pT"};
  Configurable<float> cfg_Track_MaxEta{"cfg_Track_MaxEta", 0.9, "set track max Eta"};
  Configurable<double> cfg_Track_MaxDCArToPVcut{"cfg_Track_MaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfg_Track_MaxDCAzToPVcut{"cfg_Track_MaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<bool> cfg_Track_PrimaryTrack{"cfg_Track_PrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfg_Track_ConnectedToPV{"cfg_Track_ConnectedToPV", true, "PV contributor track selection"};           // PV Contributor
  Configurable<bool> cfg_Track_GlobalWoDCATrack{"cfg_Track_GlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<double> cfg_Track_nFindableTPCClusters{"cfg_Track_FindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfg_Track_nTPCCrossedRows{"cfg_Track_TPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfg_Track_nRowsOverFindable{"cfg_Track_RowsOverFindable", 1.2, "nRowsOverFindable TPC CLusters"};
  Configurable<double> cfg_Track_nTPCChi2{"cfg_Track_TPCChi2", 4.0, "nTPC Chi2 per Cluster"};
  Configurable<double> cfg_Track_nITSChi2{"cfg_Track_ITSChi2", 36.0, "nITS Chi2 per Cluster"};
  Configurable<bool> cfg_Track_TPCPID{"cfg_Track_TPCPID", true, "Enables TPC PID"};
  Configurable<bool> cfg_Track_TOFPID{"cfg_Track_TOFPID", true, "Enables TOF PID"};
  Configurable<bool> cfg_Track_Hard_TOFPID{"cfg_Track_Hard_TOFPID", true, "Enables STRICT TOF Reqruirement"};
  Configurable<float> cfg_Track_TPCPID_nSig{"cfg_Track_TPCPID_nSig", 4, "nTPC PID sigma"};
  Configurable<float> cfg_Track_TOFPID_nSig{"cfg_Track_TOFPID_nSig", 4, "nTOF PID sigma"};
  Configurable<bool> cfg_Track_Explicit_PID{"cfg_Track_Explicit_PID", true, "Enables explicit pid cehck"};
  Configurable<int> cDebugLevel{"cDebugLevel", 0, "Resolution of Debug"};

  // Pair configurables
  Configurable<int> cfg_Pair_MinvBins{"cfg_Pair_MinvBins", 300, "Number of bins for Minv axis"};
  Configurable<float> cfg_Pair_MinvMin{"cfg_Pair_MinvMin", 0.90, "Minimum Minv value"};
  Configurable<float> cfg_Pair_MinvMax{"cfg_Pair_MinvMax", 1.50, "Maximum Minv value"};

  Configurable<int> cfg_Mix_NMixedEvents{"cfg_Mix_NMixedEvents", 5, "Number of mixed events per event"};

  // MCGen configurables
  Configurable<bool> cfg_Force_GenReco{"cfg_Force_GenReco", false, "Only consider events which are reconstructed (neglect event-loss)"};
  Configurable<bool> cfg_Force_BR{"cfg_Force_BR", false, "Only consider phi->K+K-"};
  Configurable<bool> cfg_Force_Kaon_Acceptence{"cfg_Force_Kaon_Acceptence", false, "Only consider phi's whose daughters decay inside acceptence (no signal loss)"};

  // Histogram Configurables
  Configurable<bool> cfg_Event_CutQA{"cfg_Event_CutsQA", true, "Enables Track QA plots"};
  Configurable<bool> cfg_Track_CutQA{"cfg_Track_CutsQA", true, "Enables Track QA plots"};

  // Configurables for axis
  ConfigurableAxis binsDCAz{"binsDCAz", {40, -0.2, 0.2}, ""};
  ConfigurableAxis binsDCAxy{"binsDCAxy", {40, -0.2, 0.2}, ""};
  ConfigurableAxis cfg_bins_Cent{"cfg_bins_Cent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
  ConfigurableAxis cfg_bins_MixVtx{"cfg_bins_MixVtx", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfg_bins_MixMult{"cfg_bins_MixMult", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f}, "Mixing bins - z-vertex"};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec MinvAxis = {cfg_Pair_MinvBins, cfg_Pair_MinvMin, cfg_Pair_MinvMax};
    const AxisSpec PtAxis = {200, 0, 20.0};
    const AxisSpec MultAxis = {100, 0, 100};
    const AxisSpec dRAxis = {100, 0, 100};
    const AxisSpec pidAxis = {100, -5, 5};
    const AxisSpec axisDCAz{binsDCAz, "DCA_{z}"};
    const AxisSpec axisDCAxy{binsDCAxy, "DCA_{XY}"};

    // Event QA
    if (cfg_Event_CutQA) {
      histos.add("hPosZ_BC", "PosZ_BC", kTH1F, {{240, -12.0, 12.0}});
      histos.add("hcentFT0_BC", "centFT0_BC", kTH1F, {{110, 0.0, 110.0}});
      histos.add("hOccupancy_BC", "Occupancy_BC", kTH1F, {{100, 0.0, 20000}});
      //
      histos.add("hcentFT0_AC", "centFT0_AC", kTH1F, {{110, 0.0, 110.0}});
      histos.add("hPosZ_AC", "PosZ_AC", kTH1F, {{240, -12.0, 12.0}});
      histos.add("hOccupancy_AC", "Occupancy_AC", kTH1F, {{100, 0.0, 20000}});
    }
    // Track QA
    if (cfg_Track_CutQA) {
      histos.add("hDCArToPv_BC", "DCArToPv_BC", kTH1F, {axisDCAxy});
      histos.add("hDCAzToPv_BC", "DCAzToPv_BC", kTH1F, {axisDCAz});
      histos.add("hIsPrim_BC", "hIsPrim_BC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hIsGood_BC", "hIsGood_BC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hIsPrimCont_BC", "hIsPrimCont_BC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hFindableTPCClusters_BC", "hFindableTPCClusters_BC", kTH1F, {{200, 0, 200}});
      histos.add("hFindableTPCRows_BC", "hFindableTPCRows_BC", kTH1F, {{200, 0, 200}});
      histos.add("hClustersVsRows_BC", "hClustersVsRows_BC", kTH1F, {{200, 0, 2}});
      histos.add("hTPCChi2_BC", "hTPCChi2_BC", kTH1F, {{200, 0, 100}});
      histos.add("hITSChi2_BC", "hITSChi2_BC", kTH1F, {{200, 0, 100}});
      histos.add("hTPC_nSigma_BC", "hTPC_nSigma_BC", kTH1F, {pidAxis});
      histos.add("hTOF_nSigma_BC", "hTOF_nSigma_BC", kTH1F, {pidAxis});
      histos.add("hTPC_nSigma_v_pt_BC", "hTPC_nSigma_v_pt_BC", HistType::kTHnSparseD, {pidAxis, PtAxis});
      histos.add("hTOF_nSigma_v_pt_BC", "hTOF_nSigma_v_pt_BC", HistType::kTHnSparseD, {pidAxis, PtAxis});
      //
      histos.add("hDCArToPv_AC", "DCArToPv_AC", kTH1F, {axisDCAxy});
      histos.add("hDCAzToPv_AC", "DCAzToPv_AC", kTH1F, {axisDCAz});
      histos.add("hIsPrim_AC", "hIsPrim_AC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hIsGood_AC", "hIsGood_AC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hIsPrimCont_AC", "hIsPrimCont_AC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hFindableTPCClusters_AC", "hFindableTPCClusters_AC", kTH1F, {{200, 0, 200}});
      histos.add("hFindableTPCRows_AC", "hFindableTPCRows_AC", kTH1F, {{200, 0, 200}});
      histos.add("hClustersVsRows_AC", "hClustersVsRows_AC", kTH1F, {{200, 0, 2}});
      histos.add("hTPCChi2_AC", "hTPCChi2_AC", kTH1F, {{200, 0, 100}});
      histos.add("hITSChi2_AC", "hITSChi2_AC", kTH1F, {{200, 0, 100}});
      histos.add("hTPC_nSigma_AC", "hTPC_nSigma_AC", kTH1F, {pidAxis});
      histos.add("hTOF_nSigma_AC", "hTOF_nSigma_AC", kTH1F, {pidAxis});
      histos.add("hTPC_nSigma_v_pt_AC", "hTPC_nSigma_v_pt_AC", HistType::kTHnSparseD, {pidAxis, PtAxis});
      histos.add("hTOF_nSigma_v_pt_AC", "hTOF_nSigma_v_pt_AC", HistType::kTHnSparseD, {pidAxis, PtAxis});
    }
    // Data histos
    histos.add("hUSS", "hUSS", HistType::kTHnSparseD, {cfg_bins_Cent, MinvAxis, PtAxis});
    histos.add("hLSS", "hLSS", HistType::kTHnSparseD, {cfg_bins_Cent, MinvAxis, PtAxis});
    histos.add("hUSS_Mix", "hUSS_Mix", HistType::kTHnSparseD, {cfg_bins_Cent, MinvAxis, PtAxis});
    histos.add("hLSS_Mix", "hLSS_Mix", HistType::kTHnSparseD, {cfg_bins_Cent, MinvAxis, PtAxis});

    // MC histos
    histos.add("hMC_USS", "hMC_USS", HistType::kTHnSparseD, {cfg_bins_Cent, MinvAxis, PtAxis});
    histos.add("hMC_LSS", "hMC_LSS", HistType::kTHnSparseD, {cfg_bins_Cent, MinvAxis, PtAxis});
    histos.add("hMC_USS_Mix", "hMC_USS_Mix", HistType::kTHnSparseD, {cfg_bins_Cent, MinvAxis, PtAxis});
    histos.add("hMC_LSS_Mix", "hMC_LSS_Mix", HistType::kTHnSparseD, {cfg_bins_Cent, MinvAxis, PtAxis});
    histos.add("hMC_USS_True", "hMC_USS_True", HistType::kTHnSparseD, {cfg_bins_Cent, MinvAxis, PtAxis});
    histos.add("hMC_Phi_True", "hMC_Phi_True", HistType::kTHnSparseD, {cfg_bins_Cent, PtAxis});

    // Event Histograms
    histos.add("hnEvents", "Event selection decision", kTH1I, {{10, -0.5, 9.5}});
    histos.add("hnEvents_MC", "Event selection decision", kTH1I, {{10, -0.5, 9.5}});
    histos.add("hnEvents_MC_True", "Event selection decision", kTH1I, {{10, -0.5, 9.5}});

  } // end of init

  Filter collisionFilter = nabs(aod::collision::posZ) <= cfg_Event_VtxCut;
  Filter collisionFilter_MC = nabs(aod::mccollision::posZ) <= cfg_Event_VtxCut;
  // Filter centralityFilter = nabs(aod::cent::centFT0C) <= cfg_Event_CentralityMax;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfg_Track_MaxEta && nabs(aod::track::pt) >= cfg_Track_MinPt);
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using EventCandidates_True = soa::Filtered<aod::McCollisions>;

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                  aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTOFbeta>>;
  using TrackCandidates_MC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                     aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTOFbeta, aod::McTrackLabels>>;

  using BinningTypeVtxCent = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;

  Partition<TrackCandidates_MC> PosKaon_MC =
    (aod::track::signed1Pt > static_cast<float>(0)) &&
    (!cfg_Track_TPCPID || (nabs(aod::pidtpc::tpcNSigmaKa) <= cfg_Track_TPCPID_nSig));
  Partition<TrackCandidates_MC> NegKaon_MC =
    (aod::track::signed1Pt < static_cast<float>(0)) &&
    (!cfg_Track_TPCPID || (nabs(aod::pidtpc::tpcNSigmaKa) <= cfg_Track_TPCPID_nSig));
  Partition<TrackCandidates> PosKaon =
    (aod::track::signed1Pt > static_cast<float>(0)) &&
    (!cfg_Track_TPCPID || (nabs(aod::pidtpc::tpcNSigmaKa) <= cfg_Track_TPCPID_nSig));
  Partition<TrackCandidates> NegKaon =
    (aod::track::signed1Pt < static_cast<float>(0)) &&
    (!cfg_Track_TPCPID || (nabs(aod::pidtpc::tpcNSigmaKa) <= cfg_Track_TPCPID_nSig));

  double massKa = o2::constants::physics::MassKPlus;

  // Centralicity estimator selection
  template <typename Coll>
  float centEst(Coll collisions)
  {
    float returnValue = -999.0f;
    switch (centEstimator) {
      case 0:
        returnValue = collisions.centFT0M();
        break;
      case 1:
        returnValue = collisions.centFT0A();
        break;
      case 2:
        returnValue = collisions.centFT0C();
        break;
      default:
        returnValue = collisions.centFT0M();
        break;
    }
    return returnValue;
  }

  //***********************************//
  // First, we declare some helper functions
  template <typename objType>
  void fillQA(const bool pass, const objType& obj, const int objecttype = 0)
  {

    if (objecttype == 1) {
      if constexpr (requires { obj.posZ(); }) {
        if (!pass) {
          histos.fill(HIST("hPosZ_BC"), obj.posZ());
          histos.fill(HIST("hcentFT0_BC"), centEst(obj));
        } else {
          histos.fill(HIST("hPosZ_AC"), obj.posZ());
          histos.fill(HIST("hcentFT0_AC"), centEst(obj));
        }
      }
    }
    if constexpr (requires { obj.tpcCrossedRowsOverFindableCls(); }) {
      if (objecttype == 2) {
        if (!pass) {
          histos.fill(HIST("hDCArToPv_BC"), obj.dcaXY());
          histos.fill(HIST("hDCAzToPv_BC"), obj.dcaZ());
          histos.fill(HIST("hIsPrim_BC"), obj.isPrimaryTrack());
          histos.fill(HIST("hIsGood_BC"), obj.isGlobalTrackWoDCA());
          histos.fill(HIST("hIsPrimCont_BC"), obj.isPVContributor());
          histos.fill(HIST("hFindableTPCClusters_BC"), obj.tpcNClsFindable());
          histos.fill(HIST("hFindableTPCRows_BC"), obj.tpcNClsCrossedRows());
          histos.fill(HIST("hClustersVsRows_BC"), obj.tpcCrossedRowsOverFindableCls());
          histos.fill(HIST("hTPCChi2_BC"), obj.tpcChi2NCl());
        } else {
          histos.fill(HIST("hDCArToPv_AC"), obj.dcaXY());
          histos.fill(HIST("hDCAzToPv_AC"), obj.dcaZ());
          histos.fill(HIST("hIsPrim_AC"), obj.isPrimaryTrack());
          histos.fill(HIST("hIsGood_AC"), obj.isGlobalTrackWoDCA());
          histos.fill(HIST("hIsPrimCont_AC"), obj.isPVContributor());
          histos.fill(HIST("hFindableTPCClusters_AC"), obj.tpcNClsFindable());
          histos.fill(HIST("hFindableTPCRows_AC"), obj.tpcNClsCrossedRows());
          histos.fill(HIST("hClustersVsRows_AC"), obj.tpcCrossedRowsOverFindableCls());
          histos.fill(HIST("hTPCChi2_AC"), obj.tpcChi2NCl());
        }
      }
      if (objecttype == 3) {
        if (!pass) {
          histos.fill(HIST("hTPC_nSigma_BC"), obj.tpcNSigmaKa());
          histos.fill(HIST("hTOF_nSigma_BC"), obj.tofNSigmaKa());
          histos.fill(HIST("hTPC_nSigma_v_pt_BC"), obj.tpcNSigmaKa(), obj.pt());
          histos.fill(HIST("hTOF_nSigma_v_pt_BC"), obj.tofNSigmaKa(), obj.pt());
        } else {
          histos.fill(HIST("hTPC_nSigma_AC"), obj.tpcNSigmaKa());
          histos.fill(HIST("hTOF_nSigma_AC"), obj.tofNSigmaKa());
          histos.fill(HIST("hTPC_nSigma_v_pt_AC"), obj.tpcNSigmaKa(), obj.pt());
          histos.fill(HIST("hTOF_nSigma_v_pt_AC"), obj.tofNSigmaKa(), obj.pt());
        }
      }
    }
  };
  //***********************************//

  // evsel
  template <typename EventType>
  std::pair<bool, int> eventSelection(const EventType event, const bool QA)
  {

    if (cfg_Track_CutQA && QA)
      fillQA(false, event, 1);

    if (!event.sel8())
      return {false, 1};
    if (std::abs(event.posZ()) > cfg_Event_VtxCut)
      return {false, 2};
    if (cfg_Event_Timeframe && (!event.selection_bit(aod::evsel::kNoTimeFrameBorder) || !event.selection_bit(aod::evsel::kNoITSROFrameBorder)))
      return {false, 3};
    if (cfg_Event_Timerange && (!event.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)))
      return {false, 4};
    if (cfg_Event_Pileup && (!event.selection_bit(aod::evsel::kNoSameBunchPileup) || !event.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)))
      return {false, 5};
    if (cfg_Event_Centrality && (centEst(event) > cfg_Event_CentralityMax))
      return {false, 6};
    if (cfg_Event_OccupancyCut && (event.trackOccupancyInTimeRange() > cfg_Event_MaxOccupancy))
      return {false, 7};

    if (cfg_Track_CutQA && QA)
      fillQA(true, event, 1);

    return {true, 8};
  };

  // tracksel
  template <typename TrackType>
  bool trackSelection(const TrackType track, const bool QA)
  {
    if (cfg_Track_CutQA && QA)
      fillQA(false, track, 2);

    // basic track cuts
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

    if (cfg_Track_CutQA && QA)
      fillQA(true, track, 2);
    return true;
  };

  // trackpid

  template <typename TrackPID>
  bool trackPIDKaon(const TrackPID& candidate, const bool QA)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};

    if (cfg_Track_CutQA && QA)
      fillQA(false, candidate, 3);

    if (!cfg_Track_TPCPID) {
      tpcPIDPassed = true;
    } else {
      if (std::abs(candidate.tpcNSigmaKa()) < cfg_Track_TPCPID_nSig)
        tpcPIDPassed = true;
    }

    if (!cfg_Track_TOFPID) {
      tofPIDPassed = true;
    } else {
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaKa()) < cfg_Track_TOFPID_nSig) {
          tofPIDPassed = true;
        }
      } else if (!cfg_Track_Hard_TOFPID) {
        tofPIDPassed = true;
      }
      if (!candidate.hasTOF()) {
        // std::cout << candidate.tofNSigmaKa() << std::endl;
      }
    }
    if (tpcPIDPassed && tofPIDPassed) {
      if (cfg_Track_CutQA && QA) {
        fillQA(true, candidate, 3);
      }
      return true;
    }
    return false;
  }

  template <typename CollisionType, typename TracksType>
  void TrackSlicing(const CollisionType& collision1, const TracksType&, const CollisionType& collision2, const TracksType&, const bool QA, const bool IsMix)
  {
    auto slicedtracks1 = PosKaon->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
    auto slicedtracks2 = NegKaon->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
    auto centrality = centEst(collision1);
    for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(slicedtracks1, slicedtracks2))) {
      auto [Minv, PhiPt] = minvReconstruction(track1, track2, QA);
      if (Minv < 0)
        continue;
      double conjugate = track1.sign() * track2.sign();
      if (!IsMix) {
        if (conjugate < 0) {
          histos.fill(HIST("hUSS"), centrality, Minv, PhiPt);
        } else if (conjugate > 0) {
          histos.fill(HIST("hLSS"), centrality, Minv, PhiPt);
        }
      } else {
        if (conjugate < 0) {
          histos.fill(HIST("hUSS_Mix"), centrality, Minv, PhiPt);
        } else if (conjugate > 0) {
          histos.fill(HIST("hLSS_Mix"), centrality, Minv, PhiPt);
        }
      }
    }
  } // TrackSlicing

  template <typename CollisionType, typename TracksType>
  void TrackSlicing_MC(const CollisionType& collision1, const TracksType&, const CollisionType& collision2, const TracksType&, const bool QA, const bool IsMix)
  {
    auto slicedtracks1 = PosKaon_MC->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
    auto slicedtracks2 = NegKaon_MC->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
    auto centrality = centEst(collision1);
    for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(slicedtracks1, slicedtracks2))) {
      auto [Minv, PhiPt] = minvReconstruction(track1, track2, QA);
      if (Minv < 0)
        continue;
      double conjugate = track1.sign() * track2.sign();
      if (!IsMix) {
        if (conjugate < 0) {
          histos.fill(HIST("hMC_USS"), centrality, Minv, PhiPt);
        } else if (conjugate > 0) {
          histos.fill(HIST("hMC_LSS"), centrality, Minv, PhiPt);
        }
      } else {
        if (conjugate < 0) {
          histos.fill(HIST("hMC_USS_Mix"), centrality, Minv, PhiPt);
        } else if (conjugate > 0) {
          histos.fill(HIST("hMC_LSS_Mix"), centrality, Minv, PhiPt);
        }
      }
      // now we do mc true
      if (!track1.has_mcParticle() || !track2.has_mcParticle())
        continue;
      auto part1 = track1.mcParticle();
      auto part2 = track2.mcParticle();
      if (std::fabs(part1.pdgCode()) != 321)
        continue; // Not Kaon
      if (std::fabs(part2.pdgCode()) != 321)
        continue; // Not Kaon

      if (!part1.has_mothers())
        continue; // Not decaying Kaon
      if (!part2.has_mothers())
        continue; // Not decaying Kaon

      std::vector<int> mothers1{};
      std::vector<int> mothers1PDG{};
      for (auto& part1_mom : part1.template mothers_as<aod::McParticles>()) {
        mothers1.push_back(part1_mom.globalIndex());
        mothers1PDG.push_back(part1_mom.pdgCode());
      }

      std::vector<int> mothers2{};
      std::vector<int> mothers2PDG{};
      for (auto& part2_mom : part2.template mothers_as<aod::McParticles>()) {
        mothers2.push_back(part2_mom.globalIndex());
        mothers2PDG.push_back(part2_mom.pdgCode());
      }

      if (mothers1PDG[0] != 333)
        continue; // mother not phi
      if (mothers2PDG[0] != 333)
        continue; // mother not phi

      if (mothers1[0] != mothers2[0])
        continue; // Kaons not from the same phi

      histos.fill(HIST("hMC_USS_True"), centrality, Minv, PhiPt);
    }
  } // TrackSlicing

  // Invariant mass
  template <typename TracksType>
  std::pair<double, double> minvReconstruction(const TracksType& trk1, const TracksType& trk2, const bool QA)
  {
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    //====================================================

    if (!trackSelection(trk1, QA) || !trackSelection(trk2, false))
      return {-1.0, -1.0};

    if (cfg_Track_Explicit_PID) {
      if (!trackPIDKaon(trk1, QA) || !trackPIDKaon(trk2, false))
        return {-1.0, -1.0};
    }

    if (trk1.globalIndex() >= trk2.globalIndex())
      return {-1.0, -1.0};

    lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
    lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
    lResonance = lDecayDaughter1 + lDecayDaughter2;

    return {lResonance.M(), lResonance.Pt()};
  } // MinvReconstruction

  //***************//
  // DATA
  //***************//

  int nEvents = 0;
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks)
  {
    if (cDebugLevel > 0) {
      ++nEvents;
      if (nEvents % 10000 == 0) {
        std::cout << "Processed Data Events: " << nEvents << std::endl;
      }
    }

    auto [goodEv, code] = eventSelection(collision, true);
    histos.fill(HIST("hnEvents"), code);
    if (!goodEv)
      return;
    TrackSlicing(collision, tracks, collision, tracks, true, false);

  } // end of process

  PROCESS_SWITCH(phi1020analysis, processSameEvent, "Process Same events", true);

  //***************//
  // DATA (MIX)
  //***************//

  int nEvents_Mix = 0;
  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVtxCent colBinning{{cfg_bins_MixVtx, cfg_bins_MixMult}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVtxCent> pairs{colBinning, cfg_Mix_NMixedEvents, -1, collisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (cDebugLevel > 0) {
        ++nEvents_Mix;
        if (nEvents_Mix % 10000 == 0) {
          std::cout << "Processed Mixed Events: " << nEvents_Mix << std::endl;
        }
      }
      auto [goodEv1, code1] = eventSelection(collision1, false);
      auto [goodEv2, code2] = eventSelection(collision2, false);
      if (!goodEv1 || !goodEv2)
        continue;
      TrackSlicing(collision1, tracks1, collision2, tracks2, false, true);
    } // mixing
  } // end of process
  PROCESS_SWITCH(phi1020analysis, processMixedEvent, "Process Mixed events", false);

  //***************//
  // RECONSTRUCTED MC
  //***************//

  int nEvents_MC = 0;
  void processSameEvent_MC(EventCandidates::iterator const& collision, TrackCandidates_MC const& tracks, aod::McParticles const&)
  {
    if (cDebugLevel > 0) {
      ++nEvents_MC;
      if (nEvents_MC % 10000 == 0) {
        std::cout << "Processed MC (REC) Events: " << nEvents_MC << std::endl;
      }
    }

    auto [goodEv, code] = eventSelection(collision, true);
    histos.fill(HIST("hnEvents_MC"), code);
    if (!goodEv)
      return;
    TrackSlicing_MC(collision, tracks, collision, tracks, true, false);

  } // end of process
  PROCESS_SWITCH(phi1020analysis, processSameEvent_MC, "Process Same events (MC)", true);

  //***************//
  // RECONSTRUCTED MC (MIX)
  //***************//

  int nEvents_MC_Mix = 0;
  void processMixedEvent_MC(EventCandidates const& collisions, TrackCandidates_MC const& tracks, aod::McParticles const&)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVtxCent colBinning{{cfg_bins_MixVtx, cfg_bins_MixMult}, true};
    SameKindPair<EventCandidates, TrackCandidates_MC, BinningTypeVtxCent> pairs{colBinning, cfg_Mix_NMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (cDebugLevel > 0) {
        ++nEvents_MC_Mix;
        if (nEvents_MC_Mix % 10000 == 0) {
          std::cout << "Processed Mixed Events: " << nEvents_MC_Mix << std::endl;
        }
      }
      auto [goodEv1, code1] = eventSelection(collision1, false);
      auto [goodEv2, code2] = eventSelection(collision2, false);
      if (!goodEv1 || !goodEv2)
        continue;
      TrackSlicing_MC(collision1, tracks1, collision2, tracks2, false, true);
    } // mixing
  } // end of process
  PROCESS_SWITCH(phi1020analysis, processMixedEvent_MC, "Process Mixed events (MC)", false);

  //***************//
  // GENERATED MC
  //***************//

  int nEvents_True = 0;
  void processParticles(EventCandidates_True::iterator const& collision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, EventCandidates>> const& recocolls, aod::McParticles const& particles)
  {
    if (cDebugLevel > 0) {
      ++nEvents_True;
      if (nEvents_True % 10000 == 0) {
        std::cout << "Processed MC (GEN) Events: " << nEvents_True << std::endl;
      }
    }

    if (fabs(collision.posZ()) > cfg_Event_VtxCut)
      return;

    if (recocolls.size() <= 0) { // not reconstructed
      if (cfg_Force_GenReco) {
        return;
      }
    }

    double centrality = -1;
    for (auto& recocoll : recocolls) { // poorly reconstructed
      centrality = centEst(recocoll);
      auto [goodEv, code] = eventSelection(recocoll, false);
      histos.fill(HIST("hnEvents_MC_True"), code);
      if (!goodEv)
        return;
    }

    for (auto& particle : particles) {
      if (particle.pdgCode() != 333)
        continue;
      if (std::fabs(particle.eta()) > cfg_Track_MaxEta)
        continue;

      if (cfg_Force_BR) {
        bool baddecay = false;
        for (auto& phidaughter : particle.daughters_as<aod::McParticles>()) {
          if (std::fabs(phidaughter.pdgCode()) != 321) {
            baddecay = true;
            break;
          }
          if (cfg_Force_Kaon_Acceptence) {
            if (std::fabs(phidaughter.eta()) > cfg_Track_MaxEta) {
              baddecay = true;
              break;
            }
          }
        } // loop over daughters

        if (baddecay)
          continue;
      } // enforce BR restriction

      histos.fill(HIST("hMC_Phi_True"), centrality, particle.pt());
    } // loop over particles

  } // end of process
  PROCESS_SWITCH(phi1020analysis, processParticles, "Process Particles", false);

}; // end of main struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phi1020analysis>(cfgc)};
};
