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

/// \file f0980pbpbanalysis.cxx
/// \brief f0980 resonance analysis in PbPb collisions
/// \author Junlee Kim (jikim1290@gmail.com), Sangwoo Park (sangwoo.park@cern.ch)

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
// #include <iostream>
#include <memory>
#include <string>
#include <vector>

// #include "TLorentzVector.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include <Framework/SliceCache.h>

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVector2.h"
#include <TMath.h>

// from phi
#include "Common/DataModel/PIDResponseITS.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct F0980pbpbanalysis {
  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};

  // Evnet Selection Configurables
  Configurable<float> cfgEventCutVertex{"cfgEventCutVertex", 10.0, "PV selection"};
  Configurable<bool> cfgEventQvecSel{"cfgEventQvecSel", true, "Reject events when no QVector"};
  Configurable<bool> cfgEventOccupancySel{"cfgEventOccupancySel", false, "Occupancy selection"};
  Configurable<int> cfgEventOccupancyMax{"cfgEventOccupancyMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgEventOccupancyMin{"cfgEventOccupancyMin", -100, "minimum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<bool> cfgEventGoodZvtxSel{"cfgEventGoodZvtxSel", true, "kIsGoodZvtxFT0vsPV selection"};
  Configurable<bool> cfgEventNSamePileupSel{"cfgEventNSamePileupSel", true, "kNoSameBunchPileup selection"};
  Configurable<bool> cfgEventNCollinTRSel{"cfgEventNCollinTRSel", true, "kNoCollInTimeRangeStandard selection"};
  Configurable<bool> cfgEventPVSel{"cfgEventPVSel", false, "Additional PV selection flag for syst"};
  Configurable<float> cfgEventPV{"cfgEventPV", 8.0, "Additional PV selection range for syst"};

  Configurable<float> cfgEventCentMax{"cfgEventCentMax", 80., "CentralityMax cut"};
  Configurable<int> cfgEventCentEst{"cfgEventCentEst", 1, "Centrality estimator, 1: FT0C, 2: FT0M"};

  // Track Selection Configurables
  Configurable<float> cfgTrackPtMin{"cfgTrackPtMin", 0.15, "Minimum transverse momentum for charged track"};
  Configurable<float> cfgTrackEtaMax{"cfgTrackEtaMax", 0.8, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgTrackDCArToPVcutMax{"cfgTrackDCArToPVcutMax", 0.5, "Maximum transverse DCA"};
  Configurable<float> cfgTrackDCAzToPVcutMax{"cfgTrackDCAzToPVcutMax", 2.0, "Maximum longitudinal DCA"};
  Configurable<float> cfgTrackRapMin{"cfgTrackRapMin", -0.5, "Minimum rapidity for pair"};
  Configurable<float> cfgTrackRapMax{"cfgTrackRapMax", 0.5, "Maximum rapidity for pair"};

  Configurable<bool> cfgTrackIsPrimaryTrack{"cfgTrackIsPrimaryTrack", true, "Primary track selection"};
  Configurable<bool> cfgTrackIsGlobalWoDCATrack{"cfgTrackIsGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<bool> cfgTrackIsPVContributor{"cfgTrackIsPVContributor", true, "PV contributor track selection"};

  Configurable<double> cfgTrackNTPCCrossedRows{"cfgTrackNTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgTrackNFindableTPCClusters{"cfgTrackNFindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfgTrackNRowsOverFindable{"cfgTrackNRowsOverFindable", 1.2, "nRowsOverFindable TPC CLusters"};
  Configurable<double> cfgTrackNTPCChi2{"cfgTrackNTPCChi2", 4.0, "nTPC Chi2 per Cluster"};

  Configurable<double> cfgTrackNITSChi2{"cfgTrackNITSChi2", 36.0, "nITS Chi2 per Cluster"};

  // PID Configurables
  Configurable<bool> cfgPIDUSETOF{"cfgPIDUSETOF", true, "TOF usage"};

  Configurable<double> cfgPIDMaxTOFnSigmaPion{"cfgPIDMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"}; // TOF
  Configurable<double> cfgPIDMaxTPCnSigmaPion{"cfgPIDMaxTPCnSigmaPion", 5.0, "TPC nSigma cut for Pion"}; // TPC
  Configurable<double> cfgPIDMaxTPCnSigmaPionS{"cfgPIDMaxTPCnSigmaPionS", 3.0, "TPC nSigma cut for Pion as a standalone"};
  Configurable<double> cfgPIDMaxTiednSigmaPion{"cfgPIDMaxTiednSigmaPion", 3.0, "Combined nSigma cut for Pion"};

  // Flow Configurables
  Configurable<int> cfgQvecNMods{"cfgQvecNMods", 1, "The number of modulations of interest starting from 2"};
  Configurable<int> cfgQvecNum{"cfgQvecNum", 7, "The number of total Qvectors for looping over the task"};

  Configurable<std::string> cfgQvecDetName{"cfgQvecDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgQvecRefAName{"cfgQvecRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgQvecRefBName{"cfgQvecRefBName", "TPCneg", "The name of detector for reference B"};

  // Rotational Background Configurables
  Configurable<bool> cfgBkgRotSel{"cfgBkgRotSel", true, "flag to construct rotational backgrounds"};
  Configurable<int> cfgBkgRotNum{"cfgBkgRotNum", 10, "the number of rotational backgrounds"};

  // Mixed Event Background Configurables
  SliceCache cache;
  Configurable<int> cfgBkgMixedNum{"cfgBkgMixedNum", 10, "Number of mixed events per event"};
  ConfigurableAxis mixAxisVertex{"mixAxisVertex", {10, -10, 10}, "Vertex axis for mixing bin"};
  ConfigurableAxis mixAxisCent{"mixAxisCent", {VARIABLE_WIDTH, 0, 10, 20, 50, 100}, "multiplicity percentile for mixing bin"};
  // ConfigurableAxis mixingAxisMultiplicity{"mixingAxisMultiplicity", {2000, 0, 10000}, "TPC multiplicity for bin"};

  // List Configurables
  Configurable<int> cfgListPID{"cfgListPID", 0, "PID selection type"};
  Configurable<int> cfgListPtl{"cfgListPtl", 0, "Particle selection type"};
  Configurable<int> cfgListPair{"cfgListPair", 2, "Pair selection type"};

  // Histogram QA Configurables
  Configurable<bool> cfgQAEventCut{"cfgQAEventCut", true, "Enable Event QA Hists"};
  Configurable<bool> cfgQATrackCut{"cfgQATrackCut", true, "Enable Track QA Hists"};
  Configurable<bool> cfgQAPIDCut{"cfgQAPIDCut", true, "Enable PID QA Hists"};
  Configurable<bool> cfgQAEPCut{"cfgQAEPCut", true, "Enable Event Plane QA Hists"};

  ConfigurableAxis histAxisDCAz{"histAxisDCAz", {40, -0.2, 0.2}, "DCAz axis"};
  ConfigurableAxis histAxisDCAr{"histAxisDCAr", {40, -0.2, 0.2}, "DCAxy axis"};
  ConfigurableAxis histAxisOccupancy{"histAxisOccupancy", {100, 0.0, 20000}, "Occupancy axis"};

  // Configurable for axis
  ConfigurableAxis axisMass{"axisMass", {400, 0.2, 2.2}, "Invariant mass axis"};
  ConfigurableAxis axisPT{"axisPT", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 13.0, 20.0}, "Transverse momentum Binning"};
  ConfigurableAxis axisCent{"axisCent", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100}, "Centrality interval"};
  ConfigurableAxis axisEp{"axisEp", {6, 0.0, o2::constants::math::TwoPI}, "EP axis"};

  // for phi test
  Configurable<bool> cfgPhiITSClsSel{"cfgPhiITSClsSel", false, "ITS cluster selection flag"};
  Configurable<int> cfgPhiITScluster{"cfgPhiITScluster", 0, "Number of ITS cluster"};
  Configurable<bool> cfgPhiTOFBetaSel{"cfgPhiTOFBetaSel", false, "TOF beta cut selection flag"};
  Configurable<float> cfgPhiTOFBetaCut{"cfgPhiTOFBetaCut", 0.0, "cut TOF beta"};
  Configurable<bool> cfgPhiDeepAngleSel{"cfgPhiDeepAngleSel", false, "Deep Angle cut"};
  Configurable<double> cfgPhiDeepAngle{"cfgPhiDeepAngle", 0.04, "Deep Angle cut value"};

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;

  int nmode = 2;
  static constexpr double QvecAmpMin = 1e-4;

  int detId;
  int refAId;
  int refBId;

  int qVecDetInd;
  int qVecRefAInd;
  int qVecRefBInd;

  double eventPlaneDet;
  double eventPlaneRefA;
  double eventPlaneRefB;

  float centrality;

  double angle;
  double relPhi;
  double relPhiRot;
  double relPhiMix;

  // double massPi = o2::constants::physics::MassPionCharged;
  double massPtl;

  int nTotalEvents = 0;

  enum CentEstList {
    FT0C = 0,
    FT0M = 1,
  };

  enum PIDList {
    PIDRun3 = 0,
    PIDRun2 = 1,
    PIDTest = 2,
  };

  enum PtlList {
    PtlPion = 0,
    PtlKaon = 1,
  };

  enum IndexSelList {
    None = 0,
    woSame = 1,
    leq = 2
  };

  enum QAList {
    QAEvent = 1,
    QAEP = 2,
    QATrack = 3,
    QAPID = 4,
    QAPIDS = 5
  };

  TRandom* rn = new TRandom();

  using EventCandidatesOrigin = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults, aod::Qvectors>;
  using TrackCandidatesOrigin = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFbeta>;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgEventCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgTrackEtaMax && nabs(aod::track::pt) > cfgTrackPtMin);
  Filter cutDCAFilter = (nabs(aod::track::dcaXY) < cfgTrackDCArToPVcutMax) && (nabs(aod::track::dcaZ) < cfgTrackDCAzToPVcutMax);

  using EventCandidates = soa::Filtered<EventCandidatesOrigin>;
  using TrackCandidates = soa::Filtered<TrackCandidatesOrigin>;

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  template <typename T>
  int getDetId(const T& name)
  {
    if (name.value == "FT0C") {
      return 0;
    } else if (name.value == "FT0A") {
      return 1;
    } else if (name.value == "FT0M") {
      return 2;
    } else if (name.value == "FV0A") {
      return 3;
    } else if (name.value == "TPCpos") {
      return 4;
    } else if (name.value == "TPCneg") {
      return 5;
    } else {
      return 0;
    }
  }

  template <typename objType>
  void fillQA(const bool pass, const objType& obj, const int objecttype = 0)
  {
    if constexpr (requires { obj.posZ(); }) {
      if (objecttype == QAEvent) {
        if (!pass) {
          histos.fill(HIST("EventQA/Vz_BC"), obj.posZ(), 1.0);
          histos.fill(HIST("EventQA/CentDist_BC"), centrality, 1.0);
          histos.fill(HIST("EventQA/Occupancy_BC"), obj.trackOccupancyInTimeRange(), 1.0);
        } else {
          histos.fill(HIST("EventQA/Vz_AC"), obj.posZ(), 1.0);
          histos.fill(HIST("EventQA/CentDist_AC"), centrality, 1.0);
          histos.fill(HIST("EventQA/Occupancy_AC"), obj.trackOccupancyInTimeRange(), 1.0);
        }
      }
      if (objecttype == QAEP) {
        if (!pass) {
          histos.fill(HIST("EventQA/EPhist_BC"), centrality, eventPlaneDet);
          histos.fill(HIST("EventQA/EPhistAB_BC"), centrality, std::cos(static_cast<float>(nmode) * (eventPlaneDet - eventPlaneRefA)));
          histos.fill(HIST("EventQA/EPhistAC_BC"), centrality, std::cos(static_cast<float>(nmode) * (eventPlaneDet - eventPlaneRefB)));
          histos.fill(HIST("EventQA/EPhistBC_BC"), centrality, std::cos(static_cast<float>(nmode) * (eventPlaneRefA - eventPlaneRefB)));
        } else {
          histos.fill(HIST("EventQA/EPhist_AC"), centrality, eventPlaneDet);
          histos.fill(HIST("EventQA/EPhistAB_AC"), centrality, std::cos(static_cast<float>(nmode) * (eventPlaneDet - eventPlaneRefA)));
          histos.fill(HIST("EventQA/EPhistAC_AC"), centrality, std::cos(static_cast<float>(nmode) * (eventPlaneDet - eventPlaneRefB)));
          histos.fill(HIST("EventQA/EPhistBC_AC"), centrality, std::cos(static_cast<float>(nmode) * (eventPlaneRefA - eventPlaneRefB)));
        }
      }
    }
    if constexpr (requires { obj.tpcCrossedRowsOverFindableCls(); }) {
      if (objecttype == QATrack) {
        if (!pass) {
          histos.fill(HIST("TrackQA/DCArToPv_BC"), obj.dcaXY());
          histos.fill(HIST("TrackQA/DCAzToPv_BC"), obj.dcaZ());
          histos.fill(HIST("TrackQA/IsPrim_BC"), obj.isPrimaryTrack());
          histos.fill(HIST("TrackQA/IsGood_BC"), obj.isGlobalTrackWoDCA());
          histos.fill(HIST("TrackQA/IsPrimCont_BC"), obj.isPVContributor());
          histos.fill(HIST("TrackQA/FindableTPCClusters_BC"), obj.tpcNClsFindable());
          histos.fill(HIST("TrackQA/FindableTPCRows_BC"), obj.tpcNClsCrossedRows());
          histos.fill(HIST("TrackQA/ClustersVsRows_BC"), obj.tpcCrossedRowsOverFindableCls());
          histos.fill(HIST("TrackQA/TPCChi2_BC"), obj.tpcChi2NCl());
        } else {
          histos.fill(HIST("TrackQA/DCArToPv_AC"), obj.dcaXY());
          histos.fill(HIST("TrackQA/DCAzToPv_AC"), obj.dcaZ());
          histos.fill(HIST("TrackQA/IsPrim_AC"), obj.isPrimaryTrack());
          histos.fill(HIST("TrackQA/IsGood_AC"), obj.isGlobalTrackWoDCA());
          histos.fill(HIST("TrackQA/IsPrimCont_AC"), obj.isPVContributor());
          histos.fill(HIST("TrackQA/FindableTPCClusters_AC"), obj.tpcNClsFindable());
          histos.fill(HIST("TrackQA/FindableTPCRows_AC"), obj.tpcNClsCrossedRows());
          histos.fill(HIST("TrackQA/ClustersVsRows_AC"), obj.tpcCrossedRowsOverFindableCls());
          histos.fill(HIST("TrackQA/TPCChi2_AC"), obj.tpcChi2NCl());
        }
      }
      if (objecttype == QAPID) {
        if (!pass) {
          histos.fill(HIST("PIDQA/Nsigma_TPC_BC"), obj.pt(), getTpcNSigma(obj));
          histos.fill(HIST("PIDQA/Nsigma_TOF_BC"), obj.pt(), getTofNSigma(obj));
          histos.fill(HIST("PIDQA/TPC_TOF_BC"), getTpcNSigma(obj), getTofNSigma(obj));
        } else {
          histos.fill(HIST("PIDQA/Nsigma_TPC_AC"), obj.pt(), getTpcNSigma(obj));
          histos.fill(HIST("PIDQA/Nsigma_TOF_AC"), obj.pt(), getTofNSigma(obj));
          histos.fill(HIST("PIDQA/TPC_TOF_AC"), getTpcNSigma(obj), getTofNSigma(obj));
        }
      }
      if (objecttype == QAPIDS) {
        if (pass) {
          histos.fill(HIST("PIDQA/Nsigma_TPC_selected"), obj.pt(), getTpcNSigma(obj));
          histos.fill(HIST("PIDQA/Nsigma_TOF_selected"), obj.pt(), getTofNSigma(obj));
          histos.fill(HIST("PIDQA/TPC_TOF_selected"), getTpcNSigma(obj), getTofNSigma(obj));
        }
      }
    }
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const bool QA)
  {
    if (cfgQAEventCut && QA)
      fillQA(false, collision, 1);
    if (cfgQAEPCut && QA)
      fillQA(false, collision, 2);
    //
    // histos.fill(HIST("EventQA/hnEvents"), 0);
    if (std::abs(collision.posZ()) > cfgEventCutVertex) {
      return 0;
    }
    histos.fill(HIST("EventQA/hnEvents"), 1);
    if (!collision.sel8()) {
      return 0;
    }
    histos.fill(HIST("EventQA/hnEvents"), 2);
    if (cfgEventGoodZvtxSel && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    histos.fill(HIST("EventQA/hnEvents"), 3);
    if (cfgEventNSamePileupSel && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    histos.fill(HIST("EventQA/hnEvents"), 4);
    if (cfgEventNCollinTRSel && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
    }
    histos.fill(HIST("EventQA/hnEvents"), 5);
    if (cfgEventQvecSel && (collision.qvecAmp()[detId] < QvecAmpMin || collision.qvecAmp()[refAId] < QvecAmpMin || collision.qvecAmp()[refBId] < QvecAmpMin)) {
      return 0;
    }
    histos.fill(HIST("EventQA/hnEvents"), 6);
    if (cfgEventOccupancySel && (collision.trackOccupancyInTimeRange() > cfgEventOccupancyMax || collision.trackOccupancyInTimeRange() < cfgEventOccupancyMin)) {
      return 0;
    }
    histos.fill(HIST("EventQA/hnEvents"), 7);
    if (cfgEventCentMax < centrality) {
      return 0;
    }
    /*
        auto multNTracksPV = collision.multNTracksPV();
        if (multNTracksPV < fMultPVCutLow->Eval(centrality)) {
          return 0;
        }
        if (multNTracksPV > fMultPVCutHigh->Eval(centrality)) {
          return 0;
        }
    */
    histos.fill(HIST("EventQA/hnEvents"), 8);
    if (cfgEventPVSel && std::abs(collision.posZ()) > cfgEventPV) {
      return 0;
    }
    histos.fill(HIST("EventQA/hnEvents"), 9);
    histos.fill(HIST("EventQA/hnEvents"), 10);
    return 1;
  } // event selection

  template <typename TrackType>
  bool trackSelected(const TrackType track, const bool QA)
  {
    if (cfgQATrackCut && QA)
      fillQA(false, track, 3);
    //
    if (std::abs(track.pt()) < cfgTrackPtMin) {
      return 0;
    }
    if (std::abs(track.eta()) > cfgTrackEtaMax) {
      return 0;
    }
    if (std::abs(track.dcaXY()) > cfgTrackDCArToPVcutMax) {
      return 0;
    }
    if (std::abs(track.dcaZ()) > cfgTrackDCAzToPVcutMax) {
      return 0;
    }
    if (cfgTrackIsPVContributor && !track.isPVContributor()) {
      return 0;
    }
    if (cfgTrackIsPrimaryTrack && !track.isPrimaryTrack()) {
      return 0;
    }
    if (cfgTrackIsGlobalWoDCATrack && !track.isGlobalTrackWoDCA()) {
      return 0;
    }
    if (cfgTrackNTPCCrossedRows > 0 && track.tpcNClsCrossedRows() < cfgTrackNTPCCrossedRows) {
      return 0;
    }
    if (cfgTrackNFindableTPCClusters > 0 && track.tpcNClsFindable() < cfgTrackNFindableTPCClusters) {
      return 0;
    }
    if (cfgTrackNRowsOverFindable > 0 && track.tpcCrossedRowsOverFindableCls() > cfgTrackNRowsOverFindable) {
      return 0;
    }
    if (cfgTrackNTPCChi2 > 0 && track.tpcChi2NCl() > cfgTrackNTPCChi2) {
      return 0;
    }
    if (cfgTrackNITSChi2 > 0 && track.itsChi2NCl() > cfgTrackNITSChi2) {
      return 0;
    }
    return 1;
  }

  template <typename TrackType>
  bool selectionPID(const TrackType track, const bool QA)
  {
    if (cfgQAPIDCut && QA)
      fillQA(false, track, 4);
    //
    if (cfgListPID == PIDList::PIDRun3) {
      if (cfgPIDUSETOF) {
        if (std::abs(track.tofNSigmaPi()) > cfgPIDMaxTOFnSigmaPion) {
          return 0;
        }
        if (std::abs(track.tpcNSigmaPi()) > cfgPIDMaxTPCnSigmaPion) {
          return 0;
        }
      }
      if (std::abs(track.tpcNSigmaPi()) > cfgPIDMaxTPCnSigmaPionS) {
        return 0;
      }
    } else if (cfgListPID == PIDList::PIDRun2) {
      if (cfgPIDUSETOF) {
        if (track.hasTOF()) {
          if (std::abs(track.tofNSigmaPi()) > cfgPIDMaxTOFnSigmaPion) {
            return 0;
          }
          if (std::abs(track.tpcNSigmaPi()) > cfgPIDMaxTPCnSigmaPion) {
            return 0;
          }
        } else {
          if (std::abs(track.tpcNSigmaPi()) > cfgPIDMaxTPCnSigmaPionS) {
            return 0;
          }
        }
      } else {
        if (std::abs(track.tpcNSigmaPi()) > cfgPIDMaxTPCnSigmaPionS) {
          return 0;
        }
      }
    } else if (cfgListPID == PIDList::PIDTest) {
      if (cfgPIDUSETOF) {
        if (track.hasTOF()) {
          if ((getTpcNSigma(track) * getTpcNSigma(track) + getTofNSigma(track) * getTofNSigma(track)) > (cfgPIDMaxTiednSigmaPion * cfgPIDMaxTiednSigmaPion)) {
            return 0;
          }
        } else {
          if (std::abs(getTpcNSigma(track)) > cfgPIDMaxTPCnSigmaPionS) {
            return 0;
          }
        }
      } else {
        if (std::abs(getTpcNSigma(track)) > cfgPIDMaxTPCnSigmaPionS) {
          return 0;
        }
      }
    }
    return 1;
  }

  template <typename TrackType1, typename TrackType2>
  bool pairIndexSelection(const TrackType1 track1, const TrackType2 track2)
  {
    if (cfgListPair == IndexSelList::woSame) {
      if (track2.globalIndex() == track1.globalIndex()) {
        return 0;
      }
    } else if (cfgListPair == IndexSelList::leq) {
      if (track2.globalIndex() <= track1.globalIndex()) {
        return 0;
      }
    }
    return 1;
  }

  template <typename TrackType1, typename TrackType2>
  bool pairAngleSelection(const TrackType1 track1, const TrackType2 track2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = track1.pt();
    pt2 = track2.pt();
    pz1 = track1.pz();
    pz2 = track2.pz();
    p1 = track1.p();
    p2 = track2.p();
    angle = std::acos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (cfgPhiDeepAngleSel && angle < cfgPhiDeepAngle) {
      return 0;
    }
    return 1;
  }

  template <typename TrackType>
  float getTpcNSigma(const TrackType track)
  {
    if (cfgListPtl == PtlList::PtlPion) {
      return track.tpcNSigmaPi();
    } else {
      return track.tpcNSigmaKa();
    }
  }

  template <typename TrackType>
  float getTofNSigma(const TrackType track)
  {
    if (cfgListPtl == PtlList::PtlPion) {
      return track.tofNSigmaPi();
    } else {
      return track.tofNSigmaKa();
    }
  }

  template <bool IsMC, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks)
  {
    qVecDetInd = detId * 4 + 3 + (nmode - 2) * cfgQvecNum * 4;
    qVecRefAInd = refAId * 4 + 3 + (nmode - 2) * cfgQvecNum * 4;
    qVecRefBInd = refBId * 4 + 3 + (nmode - 2) * cfgQvecNum * 4;

    eventPlaneDet = std::atan2(collision.qvecIm()[qVecDetInd], collision.qvecRe()[qVecDetInd]) / static_cast<float>(nmode);
    eventPlaneRefA = std::atan2(collision.qvecIm()[qVecRefAInd], collision.qvecRe()[qVecRefAInd]) / static_cast<float>(nmode);
    eventPlaneRefB = std::atan2(collision.qvecIm()[qVecRefBInd], collision.qvecRe()[qVecRefBInd]) / static_cast<float>(nmode);

    fillQA(true, collision, 2); // EP QA

    ROOT::Math::PxPyPzMVector pion1, pion2, pion2Rot, reco, recoRot;
    for (const auto& trk1 : dTracks) {
      if (!trackSelected(trk1, true)) {
        continue;
      }
      fillQA(true, trk1, 3);

      if (!selectionPID(trk1, true)) {
        continue;
      }
      fillQA(true, trk1, 4);

      for (const auto& trk2 : dTracks) {
        if (!trackSelected(trk2, false)) {
          continue;
        }

        // PID
        if (!selectionPID(trk2, false)) {
          continue;
        }

        if (trk1.index() == trk2.index()) {
          fillQA(true, trk1, 5);
        }

        if (!pairIndexSelection(trk1, trk2)) {
          continue;
        }

        if (cfgPhiDeepAngleSel && !pairAngleSelection(trk1, trk2)) {
          continue;
        }

        pion1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPtl);
        pion2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPtl);
        reco = pion1 + pion2;

        if (reco.Rapidity() > cfgTrackRapMax || reco.Rapidity() < cfgTrackRapMin) {
          continue;
        }

        relPhi = TVector2::Phi_0_2pi((reco.Phi() - eventPlaneDet) * static_cast<float>(nmode));

        if (trk1.sign() * trk2.sign() < 0) {
          histos.fill(HIST("hInvMass_f0980_US_EPA"), reco.M(), reco.Pt(), centrality, relPhi);
        } else if (trk1.sign() > 0 && trk2.sign() > 0) {
          histos.fill(HIST("hInvMass_f0980_LSpp_EPA"), reco.M(), reco.Pt(), centrality, relPhi);
        } else if (trk1.sign() < 0 && trk2.sign() < 0) {
          histos.fill(HIST("hInvMass_f0980_LSmm_EPA"), reco.M(), reco.Pt(), centrality, relPhi);
        }

        if (cfgBkgRotSel && trk1.sign() * trk2.sign() < 0) {
          for (int nr = 0; nr < cfgBkgRotNum; nr++) {
            auto randomPhi = rn->Uniform(o2::constants::math::PI * 5.0 / 6.0, o2::constants::math::PI * 7.0 / 6.0);
            randomPhi += pion2.Phi();
            pion2Rot = ROOT::Math::PxPyPzMVector(pion2.Pt() * std::cos(randomPhi), pion2.Pt() * std::sin(randomPhi), trk2.pz(), massPtl);
            recoRot = pion1 + pion2Rot;
            relPhiRot = TVector2::Phi_0_2pi((recoRot.Phi() - eventPlaneDet) * static_cast<float>(nmode));
            histos.fill(HIST("hInvMass_f0980_USRot_EPA"), recoRot.M(), recoRot.Pt(), centrality, relPhiRot);
          }
        }
      }
    }
  }

  void processEventMixing(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    // nmode = 2; // second order
    qVecDetInd = detId * 4 + 3 + (nmode - 2) * cfgQvecNum * 4;

    auto trackTuple = std::make_tuple(tracks);
    BinningTypeVertexContributor binningOnPositions{{mixAxisVertex, mixAxisCent}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgBkgMixedNum, -1, collisions, trackTuple, &cache};
    ROOT::Math::PxPyPzMVector ptl1, ptl2, recoPtl;
    for (const auto& [c1, t1, c2, t2] : pair) {
      if (cfgEventCentEst == CentEstList::FT0C) {
        centrality = c1.centFT0C();
      } else if (cfgEventCentEst == CentEstList::FT0M) {
        centrality = c1.centFT0M();
      }
      if (!eventSelected(c1, false) || !eventSelected(c2, false)) {
        continue;
      }
      if (c1.bcId() == c2.bcId()) {
        continue;
      }
      double eventPlaneDet = std::atan2(c1.qvecIm()[qVecDetInd], c1.qvecRe()[qVecDetInd]) / static_cast<float>(nmode);

      for (const auto& trk1 : t1) {
        if (!trackSelected(trk1, false)) {
          continue;
        }
        if (!selectionPID(trk1, false)) {
          continue;
        }

        for (const auto& trk2 : t2) {
          if (!trackSelected(trk2, false)) {
            continue;
          }
          if (!selectionPID(trk2, false)) {
            continue;
          }
          // if (!pairIndexSelection(trk1, trk2)) {
          //   continue;
          // }
          if (cfgPhiDeepAngleSel && !pairAngleSelection(trk1, trk2)) {
            continue;
          }
          ptl1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPtl);
          ptl2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPtl);
          recoPtl = ptl1 + ptl2;
          if (recoPtl.Rapidity() > cfgTrackRapMax || recoPtl.Rapidity() < cfgTrackRapMin) {
            continue;
          }

          relPhiMix = TVector2::Phi_0_2pi((recoPtl.Phi() - eventPlaneDet) * static_cast<float>(nmode));

          if (trk1.sign() * trk2.sign() < 0) {
            histos.fill(HIST("hInvMass_f0980_MixedUS_EPA"), recoPtl.M(), recoPtl.Pt(), centrality, relPhiMix);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(F0980pbpbanalysis, processEventMixing, "Process Event mixing", true);

  void processOnce(EventCandidatesOrigin const& events)
  {
    if (cfgQAEventCut) {
      nTotalEvents += events.size();
      auto hTotalEvents = histos.get<TH1>(HIST("EventQA/hnEvents"));
      if (hTotalEvents)
        hTotalEvents->SetBinContent(1, static_cast<double>(nTotalEvents));
      // std::cout << "Total number of events processed: " << nTotalEvents << std::endl;
    }
  }
  PROCESS_SWITCH(F0980pbpbanalysis, processOnce, "fill Total nEvents once", true);

  void init(o2::framework::InitContext&)
  {
    AxisSpec qaCentAxis = {110, 0, 110};
    AxisSpec qaVzAxis = {100, -20, 20};
    AxisSpec qaPIDAxis = {100, -10, 10};
    AxisSpec qaPtAxis = {200, 0, 20};
    AxisSpec qaEpAxis = {100, -1.0 * o2::constants::math::PI, o2::constants::math::PI};
    AxisSpec epresAxis = {102, -1.02, 1.02};

    // Event QA
    if (cfgQAEventCut) {
      histos.add("EventQA/CentDist_BC", "", {HistType::kTH1F, {qaCentAxis}});
      histos.add("EventQA/Vz_BC", "", {HistType::kTH1F, {qaVzAxis}});
      histos.add("EventQA/Occupancy_BC", "", kTH1F, {{histAxisOccupancy}});
    }
    histos.add("EventQA/CentDist_AC", "", {HistType::kTH1F, {qaCentAxis}});
    histos.add("EventQA/Vz_AC", "", {HistType::kTH1F, {qaVzAxis}});
    histos.add("EventQA/Occupancy_AC", "", kTH1F, {{histAxisOccupancy}});

    // Track QA
    if (cfgQATrackCut) {
      histos.add("TrackQA/DCArToPv_BC", "", {HistType::kTH1F, {histAxisDCAz}});
      histos.add("TrackQA/DCAzToPv_BC", "", {HistType::kTH1F, {histAxisDCAz}});
      histos.add("TrackQA/IsPrim_BC", "", kTH1F, {{2, -0.5, 1.5}});
      histos.add("TrackQA/IsGood_BC", "", kTH1F, {{2, -0.5, 1.5}});
      histos.add("TrackQA/IsPrimCont_BC", "", kTH1F, {{2, -0.5, 1.5}});
      histos.add("TrackQA/FindableTPCClusters_BC", "", kTH1F, {{200, 0, 200}});
      histos.add("TrackQA/FindableTPCRows_BC", "", kTH1F, {{200, 0, 200}});
      histos.add("TrackQA/ClustersVsRows_BC", "", kTH1F, {{200, 0, 2}});
      histos.add("TrackQA/TPCChi2_BC", "", kTH1F, {{200, 0, 100}});
      histos.add("TrackQA/ITSChi2_BC", "", kTH1F, {{200, 0, 100}});
      //
      histos.add("TrackQA/DCArToPv_AC", "", {HistType::kTH1F, {histAxisDCAz}});
      histos.add("TrackQA/DCAzToPv_AC", "", {HistType::kTH1F, {histAxisDCAz}});
      histos.add("TrackQA/IsPrim_AC", "", kTH1F, {{2, -0.5, 1.5}});
      histos.add("TrackQA/IsGood_AC", "", kTH1F, {{2, -0.5, 1.5}});
      histos.add("TrackQA/IsPrimCont_AC", "", kTH1F, {{2, -0.5, 1.5}});
      histos.add("TrackQA/FindableTPCClusters_AC", "", kTH1F, {{200, 0, 200}});
      histos.add("TrackQA/FindableTPCRows_AC", "", kTH1F, {{200, 0, 200}});
      histos.add("TrackQA/ClustersVsRows_AC", "", kTH1F, {{200, 0, 2}});
      histos.add("TrackQA/TPCChi2_AC", "", kTH1F, {{200, 0, 100}});
      histos.add("TrackQA/ITSChi2_AC", "", kTH1F, {{200, 0, 100}});
    }

    // PID QA
    if (cfgQAPIDCut) {
      histos.add("PIDQA/Nsigma_TPC_BC", "", {HistType::kTH2F, {qaPtAxis, qaPIDAxis}});
      histos.add("PIDQA/Nsigma_TOF_BC", "", {HistType::kTH2F, {qaPtAxis, qaPIDAxis}});
      histos.add("PIDQA/TPC_TOF_BC", "", {HistType::kTH2F, {qaPIDAxis, qaPIDAxis}});
    }
    histos.add("PIDQA/Nsigma_TPC_AC", "", {HistType::kTH2F, {qaPtAxis, qaPIDAxis}});
    histos.add("PIDQA/Nsigma_TOF_AC", "", {HistType::kTH2F, {qaPtAxis, qaPIDAxis}});
    histos.add("PIDQA/TPC_TOF_AC", "", {HistType::kTH2F, {qaPIDAxis, qaPIDAxis}});
    //
    histos.add("PIDQA/Nsigma_TPC_selected", "", {HistType::kTH2F, {qaPtAxis, qaPIDAxis}});
    histos.add("PIDQA/Nsigma_TOF_selected", "", {HistType::kTH2F, {qaPtAxis, qaPIDAxis}});
    histos.add("PIDQA/TPC_TOF_selected", "", {HistType::kTH2F, {qaPIDAxis, qaPIDAxis}});

    // Event Plane QA
    if (cfgQAEPCut) {
      histos.add("EventQA/EPhist_BC", "", {HistType::kTH2F, {qaCentAxis, qaEpAxis}});
      histos.add("EventQA/EPhistAB_BC", "", {HistType::kTH2F, {qaCentAxis, epresAxis}});
      histos.add("EventQA/EPhistAC_BC", "", {HistType::kTH2F, {qaCentAxis, epresAxis}});
      histos.add("EventQA/EPhistBC_BC", "", {HistType::kTH2F, {qaCentAxis, epresAxis}});
    }
    //
    histos.add("EventQA/EPhist_AC", "", {HistType::kTH2F, {qaCentAxis, qaEpAxis}});
    histos.add("EventQA/EPhistAB_AC", "", {HistType::kTH2F, {qaCentAxis, epresAxis}});
    histos.add("EventQA/EPhistAC_AC", "", {HistType::kTH2F, {qaCentAxis, epresAxis}});
    histos.add("EventQA/EPhistBC_AC", "", {HistType::kTH2F, {qaCentAxis, epresAxis}});

    histos.add("hInvMass_f0980_US_EPA", "unlike invariant mass",
               {HistType::kTHnSparseF, {axisMass, axisPT, axisCent, axisEp}});
    histos.add("hInvMass_f0980_LSpp_EPA", "++ invariant mass",
               {HistType::kTHnSparseF, {axisMass, axisPT, axisCent, axisEp}});
    histos.add("hInvMass_f0980_LSmm_EPA", "-- invariant mass",
               {HistType::kTHnSparseF, {axisMass, axisPT, axisCent, axisEp}});
    histos.add("hInvMass_f0980_USRot_EPA", "unlike invariant mass Rotation",
               {HistType::kTHnSparseF, {axisMass, axisPT, axisCent, axisEp}});
    histos.add("hInvMass_f0980_MixedUS_EPA", "unlike invariant mass EventMixing",
               {HistType::kTHnSparseF, {axisMass, axisPT, axisCent, axisEp}});
    //    if (doprocessMCLight) {
    //      histos.add("MCL/hpT_f0980_GEN", "generated f0 signals", HistType::kTH1F, {qaPtAxis});
    //      histos.add("MCL/hpT_f0980_REC", "reconstructed f0 signals", HistType::kTH3F, {axisMass, qaPtAxis, axisCent});
    //    }

    // Event Histograms
    histos.add("EventQA/hnEvents", "Event selection steps", {HistType::kTH1F, {{11, -0.5, 10.5}}});
    std::shared_ptr<TH1> hEventsCutFlow = histos.get<TH1>(HIST("EventQA/hnEvents"));
    std::vector<std::string> eventCutLabels = {
      "All Events",
      "Zvtx",
      "sel8",
      "GoodZvtxFT0vsPV",
      "NoSameBunchPileup",
      "NoCollInTimeRangeStandard",
      "Qvec Amplitude",
      "Occupancy",
      "Centrality",
      "Additional PV cut",
      "Passed Events"};
    for (size_t i = 0; i < eventCutLabels.size(); ++i) {
      hEventsCutFlow->GetXaxis()->SetBinLabel(i + 1, eventCutLabels[i].c_str());
    }

    detId = getDetId(cfgQvecDetName);
    refAId = getDetId(cfgQvecRefAName);
    refBId = getDetId(cfgQvecRefBName);

    if (detId == refAId || detId == refBId || refAId == refBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      detId = 0;
      refAId = 4;
      refBId = 5;
    }

    if (cfgListPtl == PtlList::PtlPion) {
      massPtl = o2::constants::physics::MassPionCharged;
    } else if (cfgListPtl == PtlList::PtlKaon) {
      massPtl = o2::constants::physics::MassKaonCharged;
    }

    fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
    fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);

    ccdb->setURL(cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  void processData(EventCandidates::iterator const& collision,
                   TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    if (cfgEventCentEst == CentEstList::FT0C) {
      centrality = collision.centFT0C();
    } else if (cfgEventCentEst == CentEstList::FT0M) {
      centrality = collision.centFT0M();
    }
    if (!eventSelected(collision, true)) {
      return;
    }
    fillQA(true, collision, 1); // Event QA

    fillHistograms<false>(collision, tracks);
  };
  PROCESS_SWITCH(F0980pbpbanalysis, processData, "Process Event for data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<F0980pbpbanalysis>(cfgc)};
}
