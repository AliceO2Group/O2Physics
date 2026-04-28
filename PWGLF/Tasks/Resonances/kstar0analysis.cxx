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
/// \file kstar0analysis.cxx
/// \brief the reconstruction of k*0(892) resonance analysis using K\pi decay channel
/// \author Hirak Kumar Koley <hirak.koley@cern.ch>

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
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <stdlib.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct kstar0analysis {
  SliceCache cache;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  //==================================
  //||
  //||         Selection
  //||
  //==================================

  // Event Selection
  Configurable<std::string> cfgEventSelections{"cfgEventSelections", "sel8", "Set event selection"};
  Configurable<float> cfgEventVtxCut{"cfgEventVtxCut", 10.0, "V_z cut selection"};
  ConfigurableAxis cfgCentAxis{"cfgCentAxis", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0}, "Binning of the centrality axis"};
  Configurable<bool> cfgOccupancySel{"cfgOccupancySel", false, "Occupancy selection"};
  Configurable<float> cfgOccupancyMax{"cfgOccupancyMax", 999999., "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<float> cfgOccupancyMin{"cfgOccupancyMin", -100., "minimum occupancy of tracks in neighbouring collisions in a given time range"};

  // Track Selection
  // General
  Configurable<double> cfgTrackMinPt{"cfgTrackMinPt", 0.15, "set track min pT"};
  Configurable<double> cfgTrackMaxEta{"cfgTrackMaxEta", 0.8, "set track max Eta"};
  Configurable<double> cfgTrackMaxDCArToPVcut{"cfgTrackMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfgTrackMaxDCAzToPVcut{"cfgTrackMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<bool> cfgTrackGlobalSel{"cfgTrackGlobalSel", true, "Global track selection"};
  Configurable<bool> cfgTrackPrimaryTrack{"cfgTrackPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgTrackConnectedToPV{"cfgTrackConnectedToPV", true, "PV contributor track selection"};           // PV Contriuibutor
  Configurable<bool> cfgTrackGlobalWoDCATrack{"cfgTrackGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)

  // TPC
  Configurable<double> cfgTrackFindableTPCClusters{"cfgTrackFindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfgTrackTPCCrossedRows{"cfgTrackTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgTrackRowsOverFindable{"cfgTrackRowsOverFindable", 1.2, "nRowsOverFindable TPC CLusters"};
  Configurable<double> cfgTrackTPCChi2{"cfgTrackTPCChi2", 4.0, "nTPC Chi2 per Cluster"};

  // ITS
  Configurable<double> cfgTrackITSChi2{"cfgTrackITSChi2", 36.0, "nITS Chi2 per Cluster"};

  // PID
  Configurable<bool> cfgTrackTPCPID{"cfgTrackTPCPID", true, "Enables TPC PID"};
  Configurable<bool> cfgTrackTOFPID{"cfgTrackTOFPID", true, "Enables TOF PID"};
  Configurable<bool> cfgTrackSquarePIDCut{"cfgTrackSqurePIDCut", true, "Enables PID cut shape square switch"};
  Configurable<bool> cfgTrackCirclePIDCut{"cfgTrackCirclePIDCut", true, "Enables PID cut shape circle switch"};
  Configurable<int> cfgTrackCircleValue{"cfgTrackCircleValue", 2, "Enables TOF TPC PID circle cut value"};
  Configurable<bool> cfgTrackTOFHard{"cfgTrackTOFHard", false, "Enables TOF Hard"};

  Configurable<float> cfgTrackTPCPIDnSig{"cfgTrackTPCPIDnSig", 4.0, "nTPC PID sigma"};
  Configurable<float> cfgTrackTOFPIDnSig{"cfgTrackTOFPIDnSig", 4.0, "nTOF PID sigma"};
  Configurable<int> cDebugLevel{"cDebugLevel", 0, "Resolution of Debug"};

  // Mixing
  ConfigurableAxis cfgBinsMixMult{"cfgBinsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
  ConfigurableAxis cfgBinsMixVtx{"cfgBinsMixVtx", {VARIABLE_WIDTH, -10.0f, -5.f, 0.f, 5.f, 10.f}, "Mixing bins - z-vertex"};
  Configurable<int> cfgMixNMixedEvents{"cfgMixNMixedEvents", 10, "Number of mixed events per event"};
  Configurable<int> cfgVtxMixCut{"cfgVtxMixCut", 10, "Vertex Mix Cut"};

  // MCGen
  Configurable<bool> cfgForceGenReco{"cfgForceGenReco", false, "Only consider events which are reconstructed (neglect event-loss)"};
  Configurable<bool> cfgReqMcEffPID{"cfgReqMcEffPID", false, "Request McEfficiency PID"};

  // Pair
  Configurable<int> cfgMinvNBins{"cfgMinvNBins", 300, "Number of bins for Minv axis"};
  Configurable<float> cfgMinvMin{"cfgMinvMin", 0.60, "Minimum Minv value"};
  Configurable<float> cfgMinvMax{"cfgMinvMax", 1.20, "Maximum Minv value"};

  // Histogram
  ConfigurableAxis binsDCAz{"binsDCAz", {40, -0.2, 0.2}, ""};
  ConfigurableAxis binsDCAxy{"binsDCAxy", {40, -0.2, 0.2}, ""};
  Configurable<bool> cfgEventCutQA{"cfgEventCutsQA", false, "Enable Event QA Hists"};
  Configurable<bool> cfgTrackCutQA{"cfgTrackCutQA", false, "Enable Track QA Hists"};

  Configurable<bool> cfgMCHistos{"cfgMCHistos", false, "Enable MC Hists"};
  Configurable<bool> cfgMixedHistos{"cfgMixedHistos", false, "Enable Mixed Histos"};
  Configurable<bool> cfgManualEvSel{"cfgManualEvSel", false, "Enable Manual EvSel"};

  // Main
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
      histos.add("hEvent_Cut", "Number of event after cuts", kTH1D, {{13, -0.5, 12.5}});
      histos.add("hPosZ_BC", "hPosZ_BC", kTH1F, {{300, -15.0, 15.0}});
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
      histos.add("hITSChi2_BC", "hITSChi2_BC", kTH1F, {{200, 0, 100}});
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
      histos.add("hITSChi2_AC", "hITSChi2_AC", kTH1F, {{200, 0, 100}});
      histos.add("QA_nSigma_pion_TPC_AC", "QA_nSigma_pion_TPC_AC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_nSigma_pion_TOF_AC", "QA_nSigma_pion_TOF_AC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_pion_TPC_TOF_AC", "QA_pion_TPC_TOF_AC", {HistType::kTH2F, {pidAxis, pidAxis}});
      histos.add("QA_nSigma_kaon_TPC_AC", "QA_nSigma_kaon_TPC_AC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_nSigma_kaon_TOF_AC", "QA_nSigma_kaon_TOF_AC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_kaon_TPC_TOF_AC", "QA_kaon_TPC_TOF_AC", {HistType::kTH2F, {pidAxis, pidAxis}});
      histos.add("QA_track_pT_AC", "QA_track_pT_AC", kTH1F, {{13, 0.0, 13.0}});
    }

    ////////////////////////////////////
    histos.add("nEvents", "nEvents", kTH1F, {{7, 0.0, 7.0}});
    histos.add("hUSS_KPi", "hUSS_KPi", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    histos.add("hUSS_PiK", "hUSS_PiK", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    histos.add("hLSS_KPi", "hLSS_KPi", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    histos.add("hLSS_PiK", "hLSS_PiK", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});

    if (cfgMixedHistos) {
      histos.add("hUSS_KPi_Mix", "hUSS_KPi_Mix", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hUSS_PiK_Mix", "hUSS_PiK_Mix", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    }
    if (cfgMCHistos) {
      histos.add("nEvents_Gen", "nEvents_Gen", kTH1F, {{4, 0.0, 4.0}});
      histos.add("hUSS_TrueRec", "hUSS_TrueRec", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hGen_pT_Raw", "Gen_pT_Raw (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hGen_pT_GoodEv", "hGen_pT_GoodEv", kTHnSparseF, {cfgCentAxis, ptAxis});
    }

    if (cfgEventCutQA) {
      std::shared_ptr<TH1> hCutFlow = histos.get<TH1>(HIST("hEvent_Cut"));
      std::vector<std::string> eventCutLabels = {
        "All Events",
        "sel8",
        Form("|Vz| < %.1f", cfgEventVtxCut.value),
        "kIsGoodZvtxFT0vsPV",
        "kNoSameBunchPileup",
        "kNoTimeFrameBorder",
        "kNoITSROFrameBorder",
        "kNoCollInTimeRangeStandard",
        "kIsGoodITSLayersAll",
        Form("Occupancy < %.0f", cfgOccupancyMax.value),
        "All passed events"};

      for (size_t i = 0; i < eventCutLabels.size(); ++i) {
        hCutFlow->GetXaxis()->SetBinLabel(i + 1, eventCutLabels[i].c_str());
      }
    }
  } // end of init

  //======================
  //|| For LF Analysis
  //======================
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
  template <typename objType>
  void fillQA(const bool pass, const objType& obj, const int objecttype = 0)
  {
    if (objecttype == 1) {
      if constexpr (requires { obj.posZ(); }) {
        if (!pass) {
          histos.fill(HIST("hPosZ_BC"), obj.posZ());
          histos.fill(HIST("hcentFT0C_BC"), obj.centFT0C());
        } else {
          histos.fill(HIST("hPosZ_AC"), obj.posZ());
          histos.fill(HIST("hcentFT0C_AC"), obj.centFT0C());
        }
      }
    } // eventSelection histogram

    if constexpr (requires { obj.tpcCrossedRowsOverFindableCls(); }) {
      if (objecttype == 3) {
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
          histos.fill(HIST("hITSChi2_BC"), obj.itsChi2NCl());
          histos.fill(HIST("QA_track_pT_BC"), obj.pt());
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
          histos.fill(HIST("hITSChi2_AC"), obj.itsChi2NCl());
          histos.fill(HIST("QA_track_pT_AC"), obj.pt());
        }
      }
    } // trackSelection
    if (objecttype == 4) {
      if constexpr (requires { obj.pt(); }) {
        if (!pass) {
          histos.fill(HIST("QA_nSigma_kaon_TPC_BC"), obj.pt(), obj.tpcNSigmaKa());
          histos.fill(HIST("QA_nSigma_kaon_TOF_BC"), obj.pt(), obj.tofNSigmaKa());
          histos.fill(HIST("QA_kaon_TPC_TOF_BC"), obj.tpcNSigmaKa(), obj.tofNSigmaKa());
        } else {
          histos.fill(HIST("QA_nSigma_kaon_TPC_AC"), obj.pt(), obj.tpcNSigmaKa());
          histos.fill(HIST("QA_nSigma_kaon_TOF_AC"), obj.pt(), obj.tofNSigmaKa());
          histos.fill(HIST("QA_kaon_TPC_TOF_AC"), obj.tpcNSigmaKa(), obj.tofNSigmaKa());
        }
      }
    } // kaon pid Selection
    if (objecttype == 5) {
      if constexpr (requires { obj.pt(); }) {
        if (!pass) {
          histos.fill(HIST("QA_nSigma_pion_TPC_BC"), obj.pt(), obj.tpcNSigmaPi());
          histos.fill(HIST("QA_nSigma_pion_TOF_BC"), obj.pt(), obj.tofNSigmaPi());
          histos.fill(HIST("QA_pion_TPC_TOF_BC"), obj.tpcNSigmaPi(), obj.tofNSigmaPi());
        } else {
          histos.fill(HIST("QA_nSigma_pion_TPC_AC"), obj.pt(), obj.tpcNSigmaPi());
          histos.fill(HIST("QA_nSigma_pion_TOF_AC"), obj.pt(), obj.tofNSigmaPi());
          histos.fill(HIST("QA_pion_TPC_TOF_AC"), obj.tpcNSigmaPi(), obj.tofNSigmaPi());
        }
      }
    } // pion pid Selection
  } // fill QA

  enum class objectType { MB,
                          MBRecParticle };

  template <typename TrackType>
  void fillMinv(objectType type, const TrackType& trk1, const TrackType& trk2, const ROOT::Math::PxPyPzMVector& lReso, double centrality, bool IsMix, bool flip)
  {
    double conjugate = trk1.sign() * trk2.sign();
    switch (type) {
      case objectType::MB:
        if (IsMix && cfgMixedHistos) {
          if (conjugate < 0) {
            if (!flip)
              histos.fill(HIST("hUSS_KPi_Mix"), centrality, lReso.Pt(), lReso.M());
            else
              histos.fill(HIST("hUSS_PiK_Mix"), centrality, lReso.Pt(), lReso.M());
          }
        } else {
          if (conjugate < 0) {
            if (!flip)
              histos.fill(HIST("hUSS_KPi"), centrality, lReso.Pt(), lReso.M());
            else
              histos.fill(HIST("hUSS_PiK"), centrality, lReso.Pt(), lReso.M());
          } else if (conjugate > 0) {
            if (!flip)
              histos.fill(HIST("hLSS_KPi"), centrality, lReso.Pt(), lReso.M());
            else
              histos.fill(HIST("hLSS_PiK"), centrality, lReso.Pt(), lReso.M());
          }
        }
        break;

      case objectType::MBRecParticle:
        if (cfgMCHistos) {
          if (conjugate < 0) {
            histos.fill(HIST("hUSS_TrueRec"), centrality, lReso.Pt(), lReso.M());
          }
        }
        break;
    } // switch
  } // fillMinv
  //======================================================================

  template <typename EventType>
  std::pair<bool, int> eventSelection(const EventType event, const bool QA)
  {
    if (cfgEventCutQA && QA) {
      fillQA(false, event, 1);
      histos.fill(HIST("hEvent_Cut"), 0);
    }

    if (!event.sel8())
      return {false, 1};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 1);
    }
    if (std::abs(event.posZ()) > cfgEventVtxCut)
      return {false, 2};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 2);
    }
    if (!event.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return {false, 3};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 3);
    }
    if (!event.selection_bit(aod::evsel::kNoSameBunchPileup))
      return {false, 4};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 4);
    }
    if (!event.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return {false, 5};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 5);
    }
    if (!event.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return {false, 6};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 6);
    }
    if (!event.selection_bit(aod::evsel::kNoCollInTimeRangeStandard))
      return {false, 7};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 7);
    }
    if (!event.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      return {false, 8};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 8);
    }
    if (cfgOccupancySel && (event.trackOccupancyInTimeRange() > cfgOccupancyMax || event.trackOccupancyInTimeRange() < cfgOccupancyMin))
      return {false, 9};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 9);
    }

    if (cfgEventCutQA && QA) {
      fillQA(true, event, 1);
      histos.fill(HIST("hEvent_Cut"), 10);
    }
    return {true, 11};
  };

  template <typename TracksType>
  bool trackSelection(const TracksType track, const bool QA)
  {
    if (cfgTrackCutQA && QA) {
      fillQA(false, track, 3);
    }

    if (cfgTrackGlobalSel && !track.isGlobalTrack())
      return false;
    if (track.pt() < cfgTrackMinPt)
      return false;
    if (std::abs(track.eta()) > cfgTrackMaxEta)
      return false;
    if (std::abs(track.dcaXY()) > cfgTrackMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cfgTrackMaxDCAzToPVcut)
      return false;
    if (cfgTrackPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgTrackGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgTrackFindableTPCClusters > 0 && track.tpcNClsFindable() < cfgTrackFindableTPCClusters)
      return false;
    if (track.tpcNClsCrossedRows() < cfgTrackTPCCrossedRows)
      return false;
    if (cfgTrackRowsOverFindable > 0 && track.tpcCrossedRowsOverFindableCls() > cfgTrackRowsOverFindable)
      return false;
    if (track.tpcChi2NCl() > cfgTrackTPCChi2)
      return false;
    if (track.itsChi2NCl() > cfgTrackITSChi2)
      return false;
    if (cfgTrackConnectedToPV && !track.isPVContributor())
      return false;

    if (cfgTrackCutQA && QA) {
      fillQA(true, track, 3);
    }
    return true;
  };

  template <typename TrackPID>
  bool trackPIDKaon(const TrackPID& candidate, const bool QA)
  {
    double tpcpid = 0;
    double tofpid = 0;
    bool tpcPIDPassed{false}, tofPIDPassed{false};

    if (cfgTrackCutQA && QA) {
      fillQA(false, candidate, 4);
    }

    // TPC
    if (cfgTrackSquarePIDCut) {
      if (std::abs(candidate.tpcNSigmaKa()) < cfgTrackTPCPIDnSig)
        tpcPIDPassed = true;
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaKa()) < cfgTrackTOFPIDnSig) {
          tofPIDPassed = true;
        }
      } else {
        if (!cfgTrackTOFHard) {
          tofPIDPassed = true;
        } else {
          tofPIDPassed = false;
        }
      }
    } // end of square cut
    if (cfgTrackCirclePIDCut) {
      if (std::abs(candidate.tpcNSigmaKa()) < cfgTrackTPCPIDnSig)
        tpcpid = std::abs(candidate.tpcNSigmaKa());
      tofpid = 0;

      if (candidate.hasTOF()) {
        tofpid = std::abs(candidate.tofNSigmaKa());
      } else {
        if (cfgTrackTOFHard) {
          tofpid = 999;
        }
      }
      if (tpcpid * tpcpid + tofpid * tofpid < cfgTrackCircleValue) {
        tpcPIDPassed = true;
        tofPIDPassed = true;
      }
    } // circular cut

    // TPC & TOF
    if (tpcPIDPassed && tofPIDPassed) {
      if (cfgTrackCutQA && QA) {
        fillQA(true, candidate, 4);
      }
      return true;
    }
    return false;
  }

  template <typename TrackPID>
  bool trackPIDPion(const TrackPID& candidate, const bool QA)
  {
    double tpcpid = 0;
    double tofpid = 0;
    bool tpcPIDPassed{false}, tofPIDPassed{false};

    if (cfgTrackCutQA && QA) {
      fillQA(false, candidate, 5);
    }

    // TPC
    if (cfgTrackSquarePIDCut) {
      if (std::abs(candidate.tpcNSigmaPi()) < cfgTrackTPCPIDnSig)
        tpcPIDPassed = true;
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaPi()) < cfgTrackTOFPIDnSig) {
          tofPIDPassed = true;
        }
      } else {
        if (!cfgTrackTOFHard) {
          tofPIDPassed = true;
        } else {
          tofPIDPassed = false;
        }
      }
    } // end of square cut
    if (cfgTrackCirclePIDCut) {
      if (std::abs(candidate.tpcNSigmaPi()) < cfgTrackTPCPIDnSig)
        tpcpid = std::abs(candidate.tpcNSigmaPi());
      tofpid = 0;

      if (candidate.hasTOF()) {
        tofpid = std::abs(candidate.tofNSigmaPi());
      } else {
        if (cfgTrackTOFHard) {
          tofpid = 999;
        }
      }
      if (tpcpid * tpcpid + tofpid * tofpid < cfgTrackCircleValue) {
        tpcPIDPassed = true;
        tofPIDPassed = true;
      }
    } // circular cut

    // TPC & TOF
    if (tpcPIDPassed && tofPIDPassed) {
      if (cfgTrackCutQA && QA) {
        fillQA(true, candidate, 5);
      }
      return true;
    }
    return false;
  }

  template <typename TracksType>
  ROOT::Math::PxPyPzMVector minvReconstruction(const TracksType& trk1, const TracksType& trk2, const bool QA, const bool flip)
  {
    if (!trackSelection(trk1, false) || !trackSelection(trk2, false))
      return {};

    if (!trackPIDKaon(trk1, QA) || !trackPIDPion(trk2, QA))
      return {};

    if (trk1.globalIndex() >= trk2.globalIndex())
      return {};

    ROOT::Math::PxPyPzMVector lDecayDaughter1, lDecayDaughter2, lResonance;
    if (!flip) {
      lDecayDaughter1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPi);
    } else {
      lDecayDaughter1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPi);
      lDecayDaughter2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massKa);
    }
    lResonance = lDecayDaughter1 + lDecayDaughter2;

    if (std::abs(lResonance.Eta()) > cfgTrackMaxEta)
      return {};
    return {lResonance};
  }

  template <typename TracksType>
  ROOT::Math::PxPyPzMVector TrueReconstruction(const TracksType& trk1, const TracksType& trk2)
  {
    double conjugate = trk1.sign() * trk2.sign();
    if (conjugate > 0)
      return {};

    auto particle1 = trk1.mcParticle();
    auto particle2 = trk2.mcParticle();

    if (!particle1.has_mothers() || !particle2.has_mothers()) {
      return {};
    }

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
      return {}; // mother not K*0
    if (mothers2PDG[0] != 313)
      return {}; // mothers not K*0

    if (mothers1[0] != mothers2[0])
      return {}; // Kaon and pion not from the same K*0

    if (std::fabs(particle1.pdgCode()) != 211 && std::fabs(particle1.pdgCode()) != 321)
      return {};
    if (std::fabs(particle2.pdgCode()) != 211 && std::fabs(particle2.pdgCode()) != 321)
      return {};

    double track1_mass, track2_mass;
    if (std::fabs(particle1.pdgCode()) == 211) {
      track1_mass = massPi;
    } else {
      track1_mass = massKa;
    }

    if (std::fabs(particle2.pdgCode()) == 211) {
      track2_mass = massPi;
    } else {
      track2_mass = massKa;
    }

    if (track1_mass == track2_mass) {
      return {};
    }

    ROOT::Math::PxPyPzMVector lTrueDaughter1, lTrueDaughter2, lTrueReso;
    lTrueDaughter1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), track1_mass);
    lTrueDaughter2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), track2_mass);
    lTrueReso = lTrueDaughter1 + lTrueDaughter2;

    if (lTrueReso.M() < 0)
      return {};

    return {lTrueReso};
  }

  template <typename CollisionType, typename TracksType>
  void TrackSlicing(const CollisionType& collision1, const TracksType&, const CollisionType& collision2, const TracksType&, const bool IsMix, const bool QA)
  {
    auto tracks1 = kaon->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
    auto tracks2 = pion->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
    auto centrality = collision1.centFT0C();

    for (const auto& [trk1, trk2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
      for (bool flip : {false, true}) {
        auto lReso = minvReconstruction(trk1, trk2, QA, flip);
        if (lReso.M() < 0)
          continue;

        fillMinv(objectType::MB, trk1, trk2, lReso, centrality, IsMix, flip);
      } // flip
    } // for
  } // TrackSlicing

  template <typename CollisionType, typename TracksType>
  void TrackSlicingMC(const CollisionType& collision1, const TracksType&, const CollisionType& collision2, const TracksType&, const bool IsMix, const bool QA)
  {
    auto tracks1 = kaonMC->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
    auto tracks2 = pionMC->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
    auto centrality = collision1.centFT0C();

    for (const auto& [trk1, trk2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
      if (!trk1.has_mcParticle() || !trk2.has_mcParticle())
        continue;

      for (bool flip : {false, true}) {
        auto lReso = minvReconstruction(trk1, trk2, QA, flip);
        if (lReso.M() < 0)
          continue;

        fillMinv(objectType::MB, trk1, trk2, lReso, centrality, IsMix, flip);
      } // flip

      //============================
      //|  True Reconstruction
      //============================
      auto lTrueReso = TrueReconstruction(trk1, trk2);
      fillMinv(objectType::MBRecParticle, trk1, trk2, lTrueReso, centrality, IsMix, false);
    } // tracks
  } // TrackSlicingMC

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

    auto [goodEv, code] = eventSelection(collision, true);
    histos.fill(HIST("nEvents"), 0.5);

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
    histos.fill(HIST("nEvents"), 1.5);

    TrackSlicing(collision, tracks, collision, tracks, false, true);

  } // processSameEvents
  PROCESS_SWITCH(kstar0analysis, processDataSameEvent, "process Data Same Event", false);

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
      auto [goodEv1, code1] = eventSelection(collision1, false);
      auto [goodEv2, code2] = eventSelection(collision2, false);
      bool VtxMixFlag = false;
      bool CentMixFlag = false;
      // bool OccupanacyMixFlag = false;
      if (std::fabs(collision1.posZ() - collision2.posZ()) <= cfgVtxMixCut) // set default to maybe 10
        VtxMixFlag = true;
      if (std::fabs(collision1.centFT0C() - collision2.centFT0C()) <= cfgVtxMixCut) // set default to maybe 10
        CentMixFlag = true;

      if (!goodEv1 || !goodEv2)
        continue;
      if (!CentMixFlag)
        continue;
      if (!VtxMixFlag)
        continue;

      TrackSlicing(collision1, tracks1, collision2, tracks2, true, false);
    }
  }
  PROCESS_SWITCH(kstar0analysis, processDataMixedEvent, "process DATA Mixed Event", false);

  //=======================================================
  //|
  //|                  MC STUFF (SE)
  //|
  //=========================================================
  int nEventsMC = 0;
  void processSameEventMC(EventCandidates::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&)
  {
    if (cDebugLevel > 0) {
      nEventsMC++;
      if ((nEventsMC + 1) % 10000 == 0) {
        double histmem = histos.getSize();
        std::cout << histmem << std::endl;
        std::cout << "process_SameEvent_MC: " << nEventsMC << std::endl;
      }
    }
    auto [goodEv, code] = eventSelection(collision, true);

    histos.fill(HIST("nEvents"), 0.5);

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

    histos.fill(HIST("nEvents"), 1.5);

    TrackSlicingMC(collision, tracks, collision, tracks, false, true);
  } // processSameEvents_MC
  PROCESS_SWITCH(kstar0analysis, processSameEventMC, "process Same Event MC", false);

  //=======================================================
  //|
  //|                  MC STUFF (ME)
  //|
  //=======================================================
  int nEventsMCMix = 0;
  void processMixedEventMC(EventCandidates const& collisions, TrackCandidatesMC const& tracks, aod::McParticles const&)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningType colBinning{{cfgBinsMixVtx, cfgBinsMixMult}, true}; // true is for 'ignore overflows' (true by default)
    SameKindPair<EventCandidates, TrackCandidatesMC, BinningType> pairs{colBinning, cfgMixNMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (cDebugLevel > 0) {
        nEventsMCMix++;
        if ((nEventsMCMix + 1) % 10000 == 0) {
          std::cout << "Processed Mixed Events: " << nEventsMCMix << std::endl;
        }
      }
      auto [goodEv1, code1] = eventSelection(collision1, false);
      auto [goodEv2, code2] = eventSelection(collision2, false);
      if (!goodEv1 || !goodEv2) {
        continue;
      }

      TrackSlicingMC(collision1, tracks1, collision2, tracks2, true, false);
    } // mixing
  } // processMixedEvent_MC
  PROCESS_SWITCH(kstar0analysis, processMixedEventMC, "process Mixed Event MC", false);

  //======================================
  //|
  //|           GENERATED STUFF
  //|
  //======================================
  int nEventsGen = 0;
  void processGen(EventCandidatesTrue::iterator const& collision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, EventCandidates>> const& recocolls, aod::McParticles const& mcParticles)
  {
    if (cDebugLevel > 0) {
      ++nEventsGen;
      if (nEventsGen % 10000 == 0) {
        std::cout << "Processed MC (GEN) Events: " << nEventsGen << std::endl;
      }
    }
    //=======================
    //|| Event & Signal loss
    //=======================
    if (cfgMCHistos) {
      histos.fill(HIST("nEvents_Gen"), 0.5);
    }

    for (auto& particle : mcParticles) {
      if (particle.pdgCode() != 313)
        continue;
      if (std::fabs(particle.eta()) > cfgTrackMaxEta)
        continue;
      if (fabs(collision.posZ()) > cfgEventVtxCut)
        break;

      if (cfgMCHistos) {
        histos.fill(HIST("hGen_pT_Raw"), particle.pt());
      }
    } // Unreconstructed collisions(=Raw coll) for correction

    if (recocolls.size() <= 0) { // not reconstructed
      if (cfgForceGenReco) {
        return;
      }
    }
    double centrality = -1;
    for (auto& recocoll : recocolls) {
      centrality = recocoll.centFT0C();
      auto [goodEv, code] = eventSelection(recocoll, true);

      if (cfgMCHistos) {
        histos.fill(HIST("nEvents_Gen"), 1.5);
      }
      if (!goodEv)
        continue;
    } // recocolls (=reconstructed collisions)

    //=================
    //|| Efficiency
    //=================
    for (auto& particle : mcParticles) {
      if (particle.pdgCode() != 313)
        continue; // Not K*0
      if (std::fabs(particle.eta()) > cfgTrackMaxEta)
        continue;

      if (cfgMCHistos) {
        histos.fill(HIST("nEvents_Gen"), 2.5);
        histos.fill(HIST("hGen_pT_GoodEv"), centrality, particle.pt());
      } // cfgMCHistos
    } // loop over particles
  } // processMCTrue
  PROCESS_SWITCH(kstar0analysis, processGen, "process Generated Particles", false);

  void processEventsDummy(EventCandidates::iterator const&, TrackCandidates const&)
  {
    return;
  }
  PROCESS_SWITCH(kstar0analysis, processEventsDummy, "dummy", false);
}; // kstar0analysis
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<kstar0analysis>(cfgc)};
};
