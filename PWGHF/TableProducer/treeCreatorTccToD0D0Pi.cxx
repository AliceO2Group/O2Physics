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

/// \file treeCreatorTccToD0D0Pi.cxx
/// \brief tree creator for studying the charm exotic state Tcc to D0D0pi
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg University
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include <set>
#include <vector>
#include <algorithm>
#include <utility>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_trkcandsel;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(PtD1, ptD1, float);
DECLARE_SOA_COLUMN(PtD2, ptD2, float);
DECLARE_SOA_COLUMN(PtPi, ptPi, float);
DECLARE_SOA_COLUMN(PxProng0D1, pxProng0D1, float);
DECLARE_SOA_COLUMN(PxProng1D1, pxProng1D1, float);
DECLARE_SOA_COLUMN(PyProng0D1, pyProng0D1, float);
DECLARE_SOA_COLUMN(PyProng1D1, pyProng1D1, float);
DECLARE_SOA_COLUMN(PzProng0D1, pzProng0D1, float);
DECLARE_SOA_COLUMN(PzProng1D1, pzProng1D1, float);
DECLARE_SOA_COLUMN(PxProng0D2, pxProng0D2, float);
DECLARE_SOA_COLUMN(PxProng1D2, pxProng1D2, float);
DECLARE_SOA_COLUMN(PyProng0D2, pyProng0D2, float);
DECLARE_SOA_COLUMN(PyProng1D2, pyProng1D2, float);
DECLARE_SOA_COLUMN(PzProng0D2, pzProng0D2, float);
DECLARE_SOA_COLUMN(PzProng1D2, pzProng1D2, float);
DECLARE_SOA_COLUMN(PxSoftPi, pxSoftPi, float);
DECLARE_SOA_COLUMN(PySoftPi, pySoftPi, float);
DECLARE_SOA_COLUMN(PzSoftPi, pzSoftPi, float);
DECLARE_SOA_COLUMN(SelFlagD1, selFlagD1, int8_t);
DECLARE_SOA_COLUMN(SelFlagD2, selFlagD2, int8_t);
DECLARE_SOA_COLUMN(MD1, mD1, float);
DECLARE_SOA_COLUMN(MD2, mD2, float);
DECLARE_SOA_COLUMN(DeltaMD1, deltaMD1, float);
DECLARE_SOA_COLUMN(DeltaMD2, deltaMD2, float);
DECLARE_SOA_COLUMN(MDD, mDD, float);
DECLARE_SOA_COLUMN(MDPi1, mDPi1, float);
DECLARE_SOA_COLUMN(MDPi2, mDPi2, float);
DECLARE_SOA_COLUMN(MDDPi, mDDPi, float);
DECLARE_SOA_COLUMN(DeltaMDDPi, deltaMDDPi, float);
DECLARE_SOA_COLUMN(EtaD1, etaD1, float);
DECLARE_SOA_COLUMN(EtaD2, etaD2, float);
DECLARE_SOA_COLUMN(EtaSoftPi, etaSoftPi, float);
DECLARE_SOA_COLUMN(PhiD1, phiD1, float);
DECLARE_SOA_COLUMN(PhiD2, phiD2, float);
DECLARE_SOA_COLUMN(PhiSoftPi, phiSoftPi, float);
DECLARE_SOA_COLUMN(YD1, yD1, float);
DECLARE_SOA_COLUMN(YD2, yD2, float);
DECLARE_SOA_COLUMN(YSoftPi, ySoftPi, float);
DECLARE_SOA_COLUMN(NSigTpcSoftPi, nSigTpcSoftPi, float);
DECLARE_SOA_COLUMN(NSigTofSoftPi, nSigTofSoftPi, float);
DECLARE_SOA_COLUMN(MlScoreD1, mlScoreD1, float);
DECLARE_SOA_COLUMN(MlScoreD2, mlScoreD2, float);
DECLARE_SOA_COLUMN(DecayLengthD1, decayLengthD1, float);
DECLARE_SOA_COLUMN(DecayLengthD2, decayLengthD2, float);
DECLARE_SOA_COLUMN(CpaD1, cpaD1, float);
DECLARE_SOA_COLUMN(CpaD2, cpaD2, float);
DECLARE_SOA_COLUMN(SignSoftPi, signSoftPi, float);
DECLARE_SOA_COLUMN(DcaXYSoftPi, dcaXYSoftPi, float);
DECLARE_SOA_COLUMN(DcaZSoftPi, dcaZSoftPi, float);
DECLARE_SOA_COLUMN(NITSClsSoftPi, nITSClsSoftPi, float);
DECLARE_SOA_COLUMN(NTPCClsCrossedRowsSoftPi, nTPCClsCrossedRowsSoftPi, float);
DECLARE_SOA_COLUMN(NTPCChi2NClSoftPi, nTPCChi2NClSoftPi, float);
DECLARE_SOA_COLUMN(CentOfCand, centOfCand, float);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

DECLARE_SOA_TABLE(HfCandTccLites, "AOD", "HFCANDTCCLITE",
                  full::PtD1,
                  full::PtD2,
                  full::PtPi,
                  full::PxProng0D1,
                  full::PxProng1D1,
                  full::PyProng0D1,
                  full::PyProng1D1,
                  full::PzProng0D1,
                  full::PzProng1D1,
                  full::PxProng0D2,
                  full::PxProng1D2,
                  full::PyProng0D2,
                  full::PyProng1D2,
                  full::PzProng0D2,
                  full::PzProng1D2,
                  full::PxSoftPi,
                  full::PySoftPi,
                  full::PzSoftPi,
                  full::SelFlagD1,
                  full::SelFlagD2,
                  full::MD1,
                  full::MD2,
                  full::DeltaMD1,
                  full::DeltaMD2,
                  full::MDD,
                  full::MDPi1,
                  full::MDPi2,
                  full::MDDPi,
                  full::DeltaMDDPi,
                  full::EtaD1,
                  full::EtaD2,
                  full::EtaSoftPi,
                  full::PhiD1,
                  full::PhiD2,
                  full::PhiSoftPi,
                  full::YD1,
                  full::YD2,
                  full::YSoftPi,
                  full::NSigTpcSoftPi,
                  full::NSigTofSoftPi,
                  full::MlScoreD1,
                  full::MlScoreD2,
                  full::DecayLengthD1,
                  full::DecayLengthD2,
                  full::CpaD1,
                  full::CpaD2,
                  full::SignSoftPi,
                  full::DcaXYSoftPi,
                  full::DcaZSoftPi,
                  full::NITSClsSoftPi,
                  full::NTPCClsCrossedRowsSoftPi,
                  full::NTPCChi2NClSoftPi,
                  full::CentOfCand);

DECLARE_SOA_TABLE(HfCandTccFullEvs, "AOD", "HFCANDTCCFULLEV",
                  full::CollisionId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorTccToD0D0Pi {
  Produces<o2::aod::HfCandTccLites> rowCandidateLite;
  Produces<o2::aod::HfCandTccFullEvs> rowCandidateFullEvents;

  Configurable<float> ptMinSoftPion{"ptMinSoftPion", 0.0, "Min pt for the soft pion"};
  Configurable<bool> usePionIsGlobalTrackWoDCA{"usePionIsGlobalTrackWoDCA", true, "check isGlobalTrackWoDCA status for pions"};

  // Configurable<float> softPiEtaMax{"softPiEtaMax", 0.9f, "Soft pion max value for pseudorapidity (abs vale)"};
  // Configurable<float> softPiChi2Max{"softPiChi2Max", 36.f, "Soft pion max value for chi2 ITS"};
  // Configurable<int> softPiItsHitMap{"softPiItsHitMap", 127, "Soft pion ITS hitmap"};
  // Configurable<int> softPiItsHitsMin{"softPiItsHitsMin", 1, "Minimum number of ITS layers crossed by the soft pion among those in \"softPiItsHitMap\""};
  Configurable<float> softPiDcaXYMax{"softPiDcaXYMax", 0.065, "Soft pion max dcaXY (cm)"};
  Configurable<float> softPiDcaZMax{"softPiDcaZMax", 0.065, "Soft pion max dcaZ (cm)"};
  Configurable<float> deltaMassCanMax{"deltaMassCanMax", 2, "delta candidate max mass (DDPi-D0D0) ((GeV/c2)"};
  Configurable<float> massCanMax{"massCanMax", 4.0, "candidate max mass (DDPi) ((GeV/c2)"};

  HfHelper hfHelper;

  using TracksPid = soa::Join<aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using TracksWPid = soa::Join<aod::TracksWCovDcaExtra, TracksPid, aod::TrackSelection>;

  using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;

  using SelectedCandidatesMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0, aod::HfMlD0>>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_d0::isSelD0 >= 1 || aod::hf_sel_candidate_d0::isSelD0bar >= 1;

  Preslice<SelectedCandidatesMl> candsD0PerCollisionWithMl = aod::track_association::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  // Partition<SelectedCandidatesMl> candidatesMlAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;

  void init(InitContext const&)
  {

    std::array<bool, 3> doprocess{doprocessDataWithMl, doprocessDataWithMlWithFT0C, doprocessDataWithMlWithFT0M};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
  }

  template <typename T>
  void fillEvent(const T& collision, int isEventReject, int runNumber)
  {
    rowCandidateFullEvents(
      collision.globalIndex(),
      collision.numContrib(),
      collision.posX(),
      collision.posY(),
      collision.posZ(),
      isEventReject,
      runNumber);
  }

  /// Evaluate centrality/multiplicity percentile (centrality estimator is automatically selected based on the used table)
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision
  template <typename Coll>
  float evaluateCentralityColl(const Coll& collision)
  {
    return o2::hf_centrality::getCentralityColl<Coll>(collision);
  }

  template <typename CollType, typename CandType, typename TrkType>
  void runCandCreatorData(CollType const& collision,
                          CandType const& candidates,
                          aod::TrackAssoc const& trackIndices,
                          TrkType const& track, aod::BCs const&)
  {
    for (const auto& candidateD1 : candidates) {
      for (auto candidateD2 = candidateD1 + 1; candidateD2 != candidates.end(); ++candidateD2) {
        for (const auto& trackId : trackIndices) {
          auto trackPion = trackId.template track_as<TrkType>();
          if (usePionIsGlobalTrackWoDCA && !trackPion.isGlobalTrackWoDCA()) {
            continue;
          }
          // minimum pT selection
          if (trackPion.pt() < ptMinSoftPion) {
            continue;
          }
          if (std::abs(trackPion.dcaXY()) > softPiDcaXYMax || std::abs(trackPion.dcaZ()) > softPiDcaZMax) {
            continue;
          }
          // avoid shared tracks
          if (
            (candidateD1.prong0Id() == candidateD2.prong0Id()) ||
            (candidateD1.prong0Id() == candidateD2.prong1Id()) ||
            (candidateD1.prong1Id() == candidateD2.prong0Id()) ||
            (candidateD1.prong1Id() == candidateD2.prong1Id()) ||
            (candidateD1.prong0Id() == trackPion.globalIndex()) ||
            (candidateD1.prong1Id() == trackPion.globalIndex()) ||
            (candidateD2.prong0Id() == trackPion.globalIndex()) ||
            (candidateD2.prong1Id() == trackPion.globalIndex())) {
            continue;
          }
          // Retrieve properties of the two D0 candidates
          float yD1 = hfHelper.yD0(candidateD1);
          float yD2 = hfHelper.yD0(candidateD2);
          float massD01 = -999;
          float massD02 = -999;
          float deltaMassD01 = -999;
          float deltaMassD02 = -999;
          int candFlagD1 = -999;
          int candFlagD2 = -999;
          float cent = evaluateCentralityColl(collision);

          std::vector<float> mlScoresD1;
          std::vector<float> mlScoresD2;

          if (candidateD1.isSelD0()) {
            candFlagD1 = (candidateD1.isSelD0bar()) ? 3 : 1;
            std::copy(candidateD1.mlProbD0().begin(), candidateD1.mlProbD0().end(), std::back_inserter(mlScoresD1));
            massD01 = hfHelper.invMassD0ToPiK(candidateD1);
          } else if (candidateD1.isSelD0bar()) {
            candFlagD1 = 2;
            std::copy(candidateD1.mlProbD0bar().begin(), candidateD1.mlProbD0bar().end(), std::back_inserter(mlScoresD1));
            massD01 = hfHelper.invMassD0barToKPi(candidateD1);
          }

          if (candidateD2.isSelD0()) {
            candFlagD2 = (candidateD2.isSelD0bar()) ? 3 : 1;
            std::copy(candidateD2.mlProbD0().begin(), candidateD2.mlProbD0().end(), std::back_inserter(mlScoresD2));
            massD02 = hfHelper.invMassD0ToPiK(candidateD2);

          } else if (candidateD2.isSelD0bar()) {
            candFlagD2 = 2;
            std::copy(candidateD2.mlProbD0bar().begin(), candidateD2.mlProbD0bar().end(), std::back_inserter(mlScoresD2));
            massD02 = hfHelper.invMassD0barToKPi(candidateD2);
          }

          // LOG(info) << " candidateD1.collisionId() " << candidateD1.collisionId()<<" massD01 "<<massD01<<" massD02 "<<massD02 <<"  trackPion.pt() "<< trackPion.pt();

          auto trackPosD1Dau = track.rawIteratorAt(candidateD1.prong0Id()); // positive daughter D01
          auto trackNegD1Dau = track.rawIteratorAt(candidateD1.prong1Id()); // negative daughter D01

          auto trackPosD2Dau = track.rawIteratorAt(candidateD2.prong0Id()); // positive daughter D02
          auto trackNegD2Dau = track.rawIteratorAt(candidateD2.prong1Id()); // negative daughter D02

          std::array<float, 3> pVecPosD1Dau{trackPosD1Dau.pVector()};
          std::array<float, 3> pVecNegD1Dau{trackNegD1Dau.pVector()};
          std::array<float, 3> pVecPosD2Dau{trackPosD2Dau.pVector()};
          std::array<float, 3> pVecNegD2Dau{trackNegD2Dau.pVector()};
          std::array<float, 3> pVecSoftPion = {trackPion.pVector()};
          std::array<double, 2> massD1Daus{MassPiPlus, MassKPlus};
          std::array<double, 2> massD2Daus{MassPiPlus, MassKPlus};

          if (candidateD1.isSelD0bar()) {

            massD1Daus[0] = MassKPlus;
            massD1Daus[1] = MassPiPlus;
          }
          if (candidateD2.isSelD0bar()) {
            massD2Daus[0] = MassKPlus;
            massD2Daus[1] = MassPiPlus;
          }

          auto massKpipi1 = RecoDecay::m(std::array{pVecPosD1Dau, pVecNegD1Dau, pVecSoftPion}, std::array{massD1Daus[0], massD1Daus[1], MassPiPlus});
          auto massKpipi2 = RecoDecay::m(std::array{pVecPosD2Dau, pVecNegD2Dau, pVecSoftPion}, std::array{massD2Daus[0], massD2Daus[1], MassPiPlus});

          deltaMassD01 = massKpipi1 - massD01;
          deltaMassD02 = massKpipi2 - massD02;

          std::array<float, 3> pVecD1{candidateD1.px(), candidateD1.py(), candidateD1.pz()};
          std::array<float, 3> pVecD2{candidateD2.px(), candidateD2.py(), candidateD2.pz()};
          auto arrayMomentaDDpi = std::array{pVecD1, pVecD2, pVecSoftPion};
          const auto massD0D0Pi = RecoDecay::m(std::move(arrayMomentaDDpi), std::array{MassD0, MassD0, MassPiPlus});
          const auto deltaMassD0D0Pi = massD0D0Pi - (massD01 + massD02);
          const auto massD0D0Pair = RecoDecay::m(std::array{pVecD1, pVecD2}, std::array{MassD0, MassD0});

          if (deltaMassD0D0Pi > deltaMassCanMax || massD0D0Pi > massCanMax) {
            continue;
          }

          rowCandidateLite(
            candidateD1.pt(),
            candidateD2.pt(),
            trackPion.pt(),
            candidateD1.pxProng0(),
            candidateD1.pxProng1(),
            candidateD1.pyProng0(),
            candidateD1.pyProng1(),
            candidateD1.pzProng0(),
            candidateD1.pzProng1(),
            candidateD2.pxProng0(),
            candidateD2.pxProng1(),
            candidateD2.pyProng0(),
            candidateD2.pyProng1(),
            candidateD2.pzProng0(),
            candidateD2.pzProng1(),
            trackPion.px(),
            trackPion.py(),
            trackPion.pz(),
            candFlagD1,
            candFlagD2,
            massD01,
            massD02,
            deltaMassD01,
            deltaMassD02,
            massD0D0Pair,
            massKpipi1,
            massKpipi2,
            massD0D0Pi,
            deltaMassD0D0Pi,
            candidateD1.eta(),
            candidateD2.eta(),
            trackPion.eta(),
            candidateD1.phi(),
            candidateD2.phi(),
            trackPion.phi(),
            yD1,
            yD2,
            trackPion.y(),
            trackPion.tpcNSigmaPi(),
            trackPion.tofNSigmaPi(),
            mlScoresD1[0],
            mlScoresD2[0],
            candidateD1.decayLength(),
            candidateD2.decayLength(),
            candidateD1.cpa(),
            candidateD2.cpa(),
            trackPion.sign(),
            trackPion.dcaXY(),
            trackPion.dcaZ(),
            trackPion.itsNCls(),
            trackPion.tpcNClsCrossedRows(),
            trackPion.tpcChi2NCl(),
            cent);
        } // end of loop track
      } // end of loop second D0
    } // end of loop first D0
  }

  void processDataWithMl(Collisions const& collisions,
                         SelectedCandidatesMl const& candidates,
                         aod::TrackAssoc const& trackIndices,
                         TracksWPid const& tracks,
                         aod::BCs const& bcs)
  {
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, collision.bc().runNumber());
      auto thisCollId = collision.globalIndex();
      auto candwD0ThisColl = candidates.sliceBy(candsD0PerCollisionWithMl, thisCollId);
      if (candwD0ThisColl.size() <= 1)
        continue; // only loop the collision that include at least 2 D candidates
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runCandCreatorData(collision, candwD0ThisColl, trackIdsThisCollision, tracks, bcs);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorTccToD0D0Pi, processDataWithMl, "Process data with DCAFitterN with the ML method and without centrality", false);

  void processDataWithMlWithFT0C(CollisionsWithFT0C const& collisions,
                                 SelectedCandidatesMl const& candidates,
                                 aod::TrackAssoc const& trackIndices,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, collision.bc().runNumber());
      auto thisCollId = collision.globalIndex();
      auto candwD0ThisColl = candidates.sliceBy(candsD0PerCollisionWithMl, thisCollId);
      if (candwD0ThisColl.size() <= 1)
        continue; // only loop the collision that include at least 2 D candidates
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runCandCreatorData(collision, candwD0ThisColl, trackIdsThisCollision, tracks, bcs);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorTccToD0D0Pi, processDataWithMlWithFT0C, "Process data with DCAFitterN with the ML method and with FT0C centrality", true);

  void processDataWithMlWithFT0M(CollisionsWithFT0M const& collisions,
                                 SelectedCandidatesMl const& candidates,
                                 aod::TrackAssoc const& trackIndices,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, collision.bc().runNumber());
      auto thisCollId = collision.globalIndex();
      auto candwD0ThisColl = candidates.sliceBy(candsD0PerCollisionWithMl, thisCollId);
      if (candwD0ThisColl.size() <= 1)
        continue; // only loop the collision that include at least 2 D candidates
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runCandCreatorData(collision, candwD0ThisColl, trackIdsThisCollision, tracks, bcs);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorTccToD0D0Pi, processDataWithMlWithFT0M, "Process data with DCAFitterN with the ML method and with FT0M centrality", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorTccToD0D0Pi>(cfgc)};
}
