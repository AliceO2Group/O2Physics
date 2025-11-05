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
/// \author Junlee Kim (jikim1290@gmail.com)

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
// #include <iostream>
#include <string>

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

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0, "PV selection"};
  Configurable<bool> cfgQvecSel{"cfgQvecSel", true, "Reject events when no QVector"};
  Configurable<bool> cfgOccupancySel{"cfgOccupancySel", false, "Occupancy selection"};
  Configurable<int> cfgOccupancyMax{"cfgOccupancyMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgOccupancyMin{"cfgOccupancyMin", -100, "minimum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<bool> cfgNCollinTR{"cfgNCollinTR", false, "Additional selection for the number of coll in time range"};
  Configurable<bool> cfgPVSel{"cfgPVSel", false, "Additional PV selection flag for syst"};
  Configurable<float> cfgPV{"cfgPV", 8.0, "Additional PV selection range for syst"};

  Configurable<float> cfgCentSel{"cfgCentSel", 80., "Centrality selection"};
  Configurable<int> cfgCentEst{"cfgCentEst", 1, "Centrality estimator, 1: FT0C, 2: FT0M"};

  Configurable<float> cfgMinPt{"cfgMinPt", 0.15, "Minimum transverse momentum for charged track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5, "Maximum transverse DCA"};
  Configurable<float> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0, "Maximum longitudinal DCA"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<float> cfgRatioTPCRowsOverFindableCls{"cfgRatioTPCRowsOverFindableCls", 0.8, "TPC Crossed Rows to Findable Clusters"};

  Configurable<float> cfgRapMin{"cfgRapMin", -0.5, "Minimum rapidity for pair"};
  Configurable<float> cfgRapMax{"cfgRapMax", 0.5, "Maximum rapidity for pair"};

  Configurable<bool> cfgIsPrimaryTrack{"cfgIsPrimaryTrack", true, "Primary track selection"};
  Configurable<bool> cfgIsGlobalWoDCATrack{"cfgIsGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<bool> cfgIsPVContributor{"cfgIsPVContributor", true, "PV contributor track selection"};

  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"}; // TOF
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 5.0, "TPC nSigma cut for Pion"}; // TPC
  Configurable<double> cMaxTPCnSigmaPionS{"cMaxTPCnSigmaPionS", 3.0, "TPC nSigma cut for Pion as a standalone"};
  Configurable<bool> cfgUSETOF{"cfgUSETOF", false, "TPC usage"};
  Configurable<int> cfgSelectPID{"cfgSelectPID", 0, "PID selection type"};
  Configurable<int> cfgSelectPtl{"cfgSelectPtl", 0, "Particle selection type"};

  Configurable<int> cfgNMods{"cfgNMods", 1, "The number of modulations of interest starting from 2"};
  Configurable<int> cfgNQvec{"cfgNQvec", 7, "The number of total Qvectors for looping over the task"};

  Configurable<std::string> cfgQvecDetName{"cfgQvecDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgQvecRefAName{"cfgQvecRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgQvecRefBName{"cfgQvecRefBName", "TPCneg", "The name of detector for reference B"};

  Configurable<bool> cfgRotBkgSel{"cfgRotBkgSel", true, "flag to construct rotational backgrounds"};
  Configurable<int> cfgRotBkgNum{"cfgRotBkgNum", 10, "the number of rotational backgrounds"};

  // for phi test
  Configurable<bool> cfgRatioTPCRowsFinableClsSel{"cfgRatioTPCRowsFinableClsSel", true, "TPC Crossed Rows to Findable Clusters selection flag"};
  Configurable<bool> cfgITSClsSel{"cfgITSClsSel", false, "ITS cluster selection flag"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<bool> cfgTOFBetaSel{"cfgTOFBetaSel", false, "TOF beta cut selection flag"};
  Configurable<float> cfgTOFBetaCut{"cfgTOFBetaCut", 0.0, "cut TOF beta"};
  Configurable<bool> cfgDeepAngleSel{"cfgDeepAngleSel", true, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<int> cfgTrackIndexSelType{"cfgTrackIndexSelType", 1, "Index selection type"};
  Configurable<double> cMaxTiednSigmaPion{"cMaxTiednSigmaPion", 3.0, "Combined nSigma cut for Pion"};

  ConfigurableAxis massAxis{"massAxis", {400, 0.2, 2.2}, "Invariant mass axis"};
  ConfigurableAxis ptAxis{"ptAxis", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 13.0, 20.0}, "Transverse momentum Binning"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100}, "Centrality interval"};
  ConfigurableAxis epAxis{"epAxis", {6, 0.0, o2::constants::math::TwoPI}, "EP axis"};

  // for event mixing
  SliceCache cache;
  Configurable<int> cfgNMixedEvents{"cfgNMixedEvents", 10, "Number of mixed events per event"};
  ConfigurableAxis mixingAxisVertex{"mixingAxisVertex", {10, -10, 10}, "Vertex axis for mixing bin"};
  ConfigurableAxis mixingAxisMultiplicity{"mixingAxisMultiplicity", {VARIABLE_WIDTH, 0, 10, 20, 50, 100}, "multiplicity percentile for mixing bin"};
  // ConfigurableAxis mixingAxisMultiplicity{"mixingAxisMultiplicity", {2000, 0, 10000}, "TPC multiplicity for bin"};

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;

  int detId;
  int refAId;
  int refBId;

  int qVecDetInd;
  int qVecRefAInd;
  int qVecRefBInd;

  float centrality;

  double angle;
  double relPhi;
  double relPhiRot;
  double relPhiMix;

  // double massPi = o2::constants::physics::MassPionCharged;
  double massPtl;

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

  TRandom* rn = new TRandom();
  // float theta2;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgMaxEta && nabs(aod::track::pt) > cfgMinPt);
  Filter cutDCAFilter = (nabs(aod::track::dcaXY) < cfgMaxDCArToPVcut) && (nabs(aod::track::dcaZ) < cfgMaxDCAzToPVcut);
  // from phi
  //  Filter centralityFilter = nabs(aod::cent::centFT0C) < cfgCentSel;
  //  Filter PIDcutFilter = nabs(aod::pidtpc::tpcNSigmaKa) < cMaxTPCnSigmaPion;
  //  Filter PIDcutFilter = nabs(aod::pidTPCFullKa::tpcNSigmaKa) < cMaxTPCnSigmaPion;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults, aod::Qvectors>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFbeta>>;
  // aod::pidTOFbeta 추가됨

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

  template <typename TCollision>
  bool eventSelected(TCollision collision)
  {
    constexpr double QvecAmpMin = 1e-4;
    if (!collision.sel8()) {
      return 0;
    }

    if (cfgCentSel < centrality) {
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
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    if (cfgQvecSel && (collision.qvecAmp()[detId] < QvecAmpMin || collision.qvecAmp()[refAId] < QvecAmpMin || collision.qvecAmp()[refBId] < QvecAmpMin)) {
      return 0;
    }
    if (cfgOccupancySel && (collision.trackOccupancyInTimeRange() > cfgOccupancyMax || collision.trackOccupancyInTimeRange() < cfgOccupancyMin)) {
      return 0;
    }
    if (cfgNCollinTR && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
    }
    if (cfgPVSel && std::abs(collision.posZ()) > cfgPV) {
      return 0;
    }
    return 1;
  } // event selection

  template <typename TrackType>
  bool trackSelected(const TrackType track)
  {
    if (std::abs(track.pt()) < cfgMinPt) {
      return 0;
    }
    if (std::fabs(track.eta()) > cfgMaxEta) {
      return 0;
    }
    if (std::fabs(track.dcaXY()) > cfgMaxDCArToPVcut) {
      return 0;
    }
    if (std::fabs(track.dcaZ()) > cfgMaxDCAzToPVcut) {
      return 0;
    }
    if (cfgIsPVContributor && !track.isPVContributor()) {
      return 0;
    }
    if (cfgIsPrimaryTrack && !track.isPrimaryTrack()) {
      return 0;
    }
    if (cfgIsGlobalWoDCATrack && !track.isGlobalTrackWoDCA()) {
      return 0;
    }
    if (track.tpcNClsFound() < cfgTPCcluster) {
      return 0;
    }
    if (cfgRatioTPCRowsFinableClsSel && track.tpcCrossedRowsOverFindableCls() < cfgRatioTPCRowsOverFindableCls) {
      return 0;
    }
    if (cfgITSClsSel && track.itsNCls() < cfgITScluster) {
      return 0;
    }
    return 1;
  }

  template <typename TrackType>
  bool selectionPID(const TrackType track)
  {
    if (cfgSelectPID == PIDList::PIDRun3) {
      if (cfgUSETOF) {
        if (std::fabs(track.tofNSigmaPi()) > cMaxTOFnSigmaPion) {
          return 0;
        }
        if (std::fabs(track.tpcNSigmaPi()) > cMaxTPCnSigmaPion) {
          return 0;
        }
      }
      if (std::fabs(track.tpcNSigmaPi()) > cMaxTPCnSigmaPionS) {
        return 0;
      }
    } else if (cfgSelectPID == PIDList::PIDRun2) {
      if (cfgUSETOF) {
        if (track.hasTOF()) {
          if (std::fabs(track.tofNSigmaPi()) > cMaxTOFnSigmaPion) {
            return 0;
          }
          if (std::fabs(track.tpcNSigmaPi()) > cMaxTPCnSigmaPion) {
            return 0;
          }
        } else {
          if (std::fabs(track.tpcNSigmaPi()) > cMaxTPCnSigmaPionS) {
            return 0;
          }
        }
      } else {
        if (std::fabs(track.tpcNSigmaPi()) > cMaxTPCnSigmaPionS) {
          return 0;
        }
      }
    } else if (cfgSelectPID == PIDList::PIDTest) {
      if (cfgUSETOF) {
        if (track.hasTOF()) {
          if ((getTpcNSigma(track) * getTpcNSigma(track) + getTofNSigma(track) * getTofNSigma(track)) > (cMaxTiednSigmaPion * cMaxTiednSigmaPion)) {
            return 0;
          }
        } else {
          if (std::fabs(getTpcNSigma(track)) > cMaxTPCnSigmaPionS) {
            return 0;
          }
        }
      } else {
        if (std::fabs(getTpcNSigma(track)) > cMaxTPCnSigmaPionS) {
          return 0;
        }
      }
    }
    return 1;
  }

  template <typename TrackType1, typename TrackType2>
  bool indexSelection(const TrackType1 track1, const TrackType2 track2)
  {
    if (cfgTrackIndexSelType == IndexSelList::woSame) {
      if (track2.globalIndex() == track1.globalIndex()) {
        return 0;
      }
    } else if (cfgTrackIndexSelType == IndexSelList::leq) {
      if (track2.globalIndex() <= track1.globalIndex()) {
        return 0;
      }
    }
    return 1;
  }

  template <typename TrackType1, typename TrackType2>
  bool selectionPair(const TrackType1 track1, const TrackType2 track2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = track1.pt();
    pt2 = track2.pt();
    pz1 = track1.pz();
    pz2 = track2.pz();
    p1 = track1.p();
    p2 = track2.p();
    angle = std::acos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (cfgDeepAngleSel && angle < cfgDeepAngle) {
      return 0;
    }
    return 1;
  }

  template <typename TrackType>
  float getTpcNSigma(const TrackType track)
  {
    if (cfgSelectPtl == PtlList::PtlPion) {
      return track.tpcNSigmaPi();
    } else {
      return track.tpcNSigmaKa();
    }
  }

  template <typename TrackType>
  float getTofNSigma(const TrackType track)
  {
    if (cfgSelectPtl == PtlList::PtlPion) {
      return track.tofNSigmaPi();
    } else {
      return track.tofNSigmaKa();
    }
  }

  template <bool IsMC, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision,
                      const TracksType& dTracks, int nmode)
  {
    qVecDetInd = detId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    qVecRefAInd = refAId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    qVecRefBInd = refBId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;

    double eventPlaneDet = std::atan2(collision.qvecIm()[qVecDetInd], collision.qvecRe()[qVecDetInd]) / static_cast<float>(nmode);
    double eventPlaneRefA = std::atan2(collision.qvecIm()[qVecRefAInd], collision.qvecRe()[qVecRefAInd]) / static_cast<float>(nmode);
    double eventPlaneRefB = std::atan2(collision.qvecIm()[qVecRefBInd], collision.qvecRe()[qVecRefBInd]) / static_cast<float>(nmode);

    histos.fill(HIST("QA/EPhist"), centrality, eventPlaneDet);
    histos.fill(HIST("QA/EPResAB"), centrality, std::cos(static_cast<float>(nmode) * (eventPlaneDet - eventPlaneRefA)));
    histos.fill(HIST("QA/EPResAC"), centrality, std::cos(static_cast<float>(nmode) * (eventPlaneDet - eventPlaneRefB)));
    histos.fill(HIST("QA/EPResBC"), centrality, std::cos(static_cast<float>(nmode) * (eventPlaneRefA - eventPlaneRefB)));

    ROOT::Math::PxPyPzMVector pion1, pion2, pion2Rot, reco, recoRot;
    for (const auto& trk1 : dTracks) {
      if (!trackSelected(trk1)) {
        continue;
      }

      if (!selectionPID(trk1)) {
        continue;
      }

      histos.fill(HIST("QA/Nsigma_TPC"), trk1.pt(), getTpcNSigma(trk1));
      histos.fill(HIST("QA/Nsigma_TOF"), trk1.pt(), getTofNSigma(trk1));
      histos.fill(HIST("QA/TPC_TOF"), getTpcNSigma(trk1), getTofNSigma(trk1));

      for (const auto& trk2 : dTracks) {
        if (!trackSelected(trk2)) {
          continue;
        }

        // PID
        if (!selectionPID(trk2)) {
          continue;
        }

        if (trk1.index() == trk2.index()) {
          histos.fill(HIST("QA/Nsigma_TPC_selected"), trk1.pt(), getTpcNSigma(trk2));
          histos.fill(HIST("QA/Nsigma_TOF_selected"), trk1.pt(), getTofNSigma(trk2));
          histos.fill(HIST("QA/TPC_TOF_selected"), getTpcNSigma(trk2), getTofNSigma(trk2));
        }

        if (!indexSelection(trk1, trk2)) {
          continue;
        }

        if (!selectionPair(trk1, trk2)) {
          continue;
        }

        pion1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPtl);
        pion2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPtl);
        reco = pion1 + pion2;

        if (reco.Rapidity() > cfgRapMax || reco.Rapidity() < cfgRapMin) {
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

        if (cfgRotBkgSel && trk1.sign() * trk2.sign() < 0) {
          for (int nr = 0; nr < cfgRotBkgNum; nr++) {
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
    int nmode = 2; // second order
    qVecDetInd = detId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;

    auto trackTuple = std::make_tuple(tracks);
    BinningTypeVertexContributor binningOnPositions{{mixingAxisVertex, mixingAxisMultiplicity}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgNMixedEvents, -1, collisions, trackTuple, &cache};
    ROOT::Math::PxPyPzMVector ptl1, ptl2, recoPtl;
    for (const auto& [c1, t1, c2, t2] : pair) {
      if (cfgCentEst == CentEstList::FT0C) {
        centrality = c1.centFT0C();
      } else if (cfgCentEst == CentEstList::FT0M) {
        centrality = c1.centFT0M();
      }
      if (!eventSelected(c1) || !eventSelected(c2)) {
        continue;
      }
      if (c1.bcId() == c2.bcId()) {
        continue;
      }
      double eventPlaneDet = std::atan2(c1.qvecIm()[qVecDetInd], c1.qvecRe()[qVecDetInd]) / static_cast<float>(nmode);

      for (const auto& trk1 : t1) {
        if (!trackSelected(trk1)) {
          continue;
        }
        if (!selectionPID(trk1)) {
          continue;
        }

        for (const auto& trk2 : t2) {
          if (!trackSelected(trk2)) {
            continue;
          }
          if (!selectionPID(trk2)) {
            continue;
          }
          if (!indexSelection(trk1, trk2)) {
            continue;
          }
          if (!selectionPair(trk1, trk2)) {
            continue;
          }
          ptl1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPtl);
          ptl2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPtl);
          recoPtl = ptl1 + ptl2;
          if (recoPtl.Rapidity() > cfgRapMax || recoPtl.Rapidity() < cfgRapMin) {
            continue;
          }

          relPhiMix = TVector2::Phi_0_2pi((recoPtl.Phi() - eventPlaneDet) * static_cast<float>(nmode));

          if (trk1.sign() * trk2.sign() < 0) {
            histos.fill(HIST("hInvMass_f0980_MixedUS_EPA"), recoPtl.M(), recoPtl.Pt(), centrality, relPhiMix);
          } // else if (trk1.sign() > 0 && trk2.sign() > 0) {
          //   histos.fill(HIST("hInvMass_f0980_MixedLSpp_EPA"), recoPtl.M(), recoPtl.Pt(), centrality, relPhiMix);
          // } else if (trk1.sign() < 0 && trk2.sign() < 0) {
          //   histos.fill(HIST("hInvMass_f0980_MixedLSmm_EPA"), recoPtl.M(), recoPtl.Pt(), centrality, relPhiMix);
          // }
        }
      }
      // for (auto& [trk1, trk2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(t1, t2))) {
      //   if (!trackSelected(trk1) || !trackSelected(trk2)) {
      //     continue;
      //   }
      //   if (!selectionPID(trk1) || !selectionPID(trk2)) {
      //     continue;
      //   }
      //   if (!indexSelection(trk1, trk2)) {
      //     continue;
      //   }
      //   if (!selectionPair(trk1, trk2)) {
      //     continue;
      //   }
      //   ptl1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPtl);
      //   ptl2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPtl);
      //   recoPtl = ptl1 + ptl2;
      //   if (recoPtl.Rapidity() > cfgRapMax || recoPtl.Rapidity() < cfgRapMin) {
      //     continue;
      //   }

      //   relPhiMix = TVector2::Phi_0_2pi((recoPtl.Phi() - eventPlaneDet) * static_cast<float>(nmode));

      //   if (trk1.sign() * trk2.sign() < 0) {
      //     histos.fill(HIST("hInvMass_f0980_MixedUS_EPA"), recoPtl.M(), recoPtl.Pt(), centrality, relPhiMix);
      //   } else if (trk1.sign() > 0 && trk2.sign() > 0) {
      //     histos.fill(HIST("hInvMass_f0980_MixedLSpp_EPA"), recoPtl.M(), recoPtl.Pt(), centrality, relPhiMix);
      //   } else if (trk1.sign() < 0 && trk2.sign() < 0) {
      //     histos.fill(HIST("hInvMass_f0980_MixedLSmm_EPA"), recoPtl.M(), recoPtl.Pt(), centrality, relPhiMix);
      //   }
      // }
    }
  }
  PROCESS_SWITCH(F0980pbpbanalysis, processEventMixing, "Process Event mixing", true);

  void init(o2::framework::InitContext&)
  {
    AxisSpec qaCentAxis = {110, 0, 110};
    AxisSpec qaVzAxis = {100, -20, 20};
    AxisSpec qaPIDAxis = {100, -10, 10};
    AxisSpec qaPtAxis = {200, 0, 20};
    AxisSpec qaEpAxis = {100, -1.0 * o2::constants::math::PI, o2::constants::math::PI};
    AxisSpec epresAxis = {102, -1.02, 1.02};

    histos.add("QA/CentDist", "", {HistType::kTH1F, {qaCentAxis}});
    histos.add("QA/Vz", "", {HistType::kTH1F, {qaVzAxis}});

    histos.add("QA/Nsigma_TPC", "", {HistType::kTH2F, {qaPtAxis, qaPIDAxis}});
    histos.add("QA/Nsigma_TOF", "", {HistType::kTH2F, {qaPtAxis, qaPIDAxis}});
    histos.add("QA/TPC_TOF", "", {HistType::kTH2F, {qaPIDAxis, qaPIDAxis}});

    histos.add("QA/Nsigma_TPC_selected", "", {HistType::kTH2F, {qaPtAxis, qaPIDAxis}});
    histos.add("QA/Nsigma_TOF_selected", "", {HistType::kTH2F, {qaPtAxis, qaPIDAxis}});
    histos.add("QA/TPC_TOF_selected", "", {HistType::kTH2F, {qaPIDAxis, qaPIDAxis}});

    histos.add("QA/EPhist", "", {HistType::kTH2F, {qaCentAxis, qaEpAxis}});
    histos.add("QA/EPResAB", "", {HistType::kTH2F, {qaCentAxis, epresAxis}});
    histos.add("QA/EPResAC", "", {HistType::kTH2F, {qaCentAxis, epresAxis}});
    histos.add("QA/EPResBC", "", {HistType::kTH2F, {qaCentAxis, epresAxis}});

    histos.add("hInvMass_f0980_US_EPA", "unlike invariant mass",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, epAxis}});
    histos.add("hInvMass_f0980_LSpp_EPA", "++ invariant mass",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, epAxis}});
    histos.add("hInvMass_f0980_LSmm_EPA", "-- invariant mass",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, epAxis}});
    histos.add("hInvMass_f0980_USRot_EPA", "unlike invariant mass Rotation",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, epAxis}});
    histos.add("hInvMass_f0980_MixedUS_EPA", "unlike invariant mass EventMixing",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, epAxis}});
    //    if (doprocessMCLight) {
    //      histos.add("MCL/hpT_f0980_GEN", "generated f0 signals", HistType::kTH1F, {qaPtAxis});
    //      histos.add("MCL/hpT_f0980_REC", "reconstructed f0 signals", HistType::kTH3F, {massAxis, qaPtAxis, centAxis});
    //    }

    detId = getDetId(cfgQvecDetName);
    refAId = getDetId(cfgQvecRefAName);
    refBId = getDetId(cfgQvecRefBName);

    if (detId == refAId || detId == refBId || refAId == refBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      detId = 0;
      refAId = 4;
      refBId = 5;
    }

    if (cfgSelectPtl == PtlList::PtlPion) {
      massPtl = o2::constants::physics::MassPionCharged;
    } else if (cfgSelectPtl == PtlList::PtlKaon) {
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
    if (cfgCentEst == CentEstList::FT0C) {
      centrality = collision.centFT0C();
    } else if (cfgCentEst == CentEstList::FT0M) {
      centrality = collision.centFT0M();
    }
    if (!eventSelected(collision)) {
      return;
    }
    histos.fill(HIST("QA/CentDist"), centrality, 1.0);
    histos.fill(HIST("QA/Vz"), collision.posZ(), 1.0);

    fillHistograms<false>(collision, tracks, 2); // second order
  };
  PROCESS_SWITCH(F0980pbpbanalysis, processData, "Process Event for data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<F0980pbpbanalysis>(cfgc)};
}
