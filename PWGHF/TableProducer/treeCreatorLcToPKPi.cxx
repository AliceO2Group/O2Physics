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

/// \file treeCreatorLcToPKPi.cxx
/// \brief Writer of the 3 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
///
/// \author Nicolo' Jacazio <nicolo.jacazio@cern.ch>, CERN
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch>, CERN

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::constants::physics;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(PtProng2, ptProng2, float);
DECLARE_SOA_COLUMN(PProng2, pProng2, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised2, impactParameterNormalised2, float);
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(MassKPi, massKPi, float); // invariant mass of the candidate Kpi daughters
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(NSigTpcPi0, nSigTpcPi0, float);
DECLARE_SOA_COLUMN(NSigTpcPr0, nSigTpcPr0, float);
DECLARE_SOA_COLUMN(NSigTofPi0, nSigTofPi0, float);
DECLARE_SOA_COLUMN(NSigTofPr0, nSigTofPr0, float);
DECLARE_SOA_COLUMN(NSigTpcKa1, nSigTpcKa1, float);
DECLARE_SOA_COLUMN(NSigTofKa1, nSigTofKa1, float);
DECLARE_SOA_COLUMN(NSigTpcPi2, nSigTpcPi2, float);
DECLARE_SOA_COLUMN(NSigTpcPr2, nSigTpcPr2, float);
DECLARE_SOA_COLUMN(NSigTofPi2, nSigTofPi2, float);
DECLARE_SOA_COLUMN(NSigTofPr2, nSigTofPr2, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr0, nSigTpcTofPi0, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi0, nSigTpcTofPr0, float);
DECLARE_SOA_COLUMN(NSigTpcTofKa1, nSigTpcTofKa1, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr2, nSigTpcTofPi2, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi2, nSigTpcTofPr2, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(Cpa, cpa, float);
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(FlagMc, flagMc, int8_t);
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);
DECLARE_SOA_COLUMN(IsCandidateSwapped, isCandidateSwapped, int8_t);
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfCand3Prong, "_0");
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);
DECLARE_SOA_COLUMN(Channel, channel, int8_t); // direct or resonant
// Events
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(CentFV0A, centFV0A, float);
DECLARE_SOA_COLUMN(CentFDDM, centFDDM, float);
DECLARE_SOA_COLUMN(MultZeqNTracksPV, multZeqNTracksPV, float);
} // namespace

namespace kf
{
DECLARE_SOA_COLUMN(Chi2PrimProton, chi2PrimProton, float);
DECLARE_SOA_COLUMN(Chi2PrimKaon, chi2PrimKaon, float);
DECLARE_SOA_COLUMN(Chi2PrimPion, chi2PrimPion, float);
DECLARE_SOA_COLUMN(DCAProtonKaon, dcaProtonKaon, float);
DECLARE_SOA_COLUMN(DCAProtonPion, dcaProtonPion, float);
DECLARE_SOA_COLUMN(DCAPionKaon, dcaPionKaon, float);
DECLARE_SOA_COLUMN(Chi2geoProtonKaon, chi2geoProtonKaon, float);
DECLARE_SOA_COLUMN(Chi2geoProtonPion, chi2geoProtonPion, float);
DECLARE_SOA_COLUMN(Chi2geoPionKaon, chi2geoPionKaon, float);
DECLARE_SOA_COLUMN(Chi2geo, chi2geo, float);     //! chi2 geo of the full candidate
DECLARE_SOA_COLUMN(Chi2topo, chi2topo, float);     //! chi2 topo of the full candidate (chi2prim of candidate to PV)
DECLARE_SOA_COLUMN(L, l, float);     //! decay length
DECLARE_SOA_COLUMN(DeltaL, deltaL, float);     //! decay length error
DECLARE_SOA_COLUMN(T, t, float);     //! lifetime
DECLARE_SOA_COLUMN(DeltaT, deltat, float);     //! lifetime error
DECLARE_SOA_COLUMN(MassInv, massInv, float);     //! invariant mass
DECLARE_SOA_COLUMN(P, p, float);     //! momentum
DECLARE_SOA_COLUMN(Pt, pt, float);     //! transverse momentum
DECLARE_SOA_COLUMN(IsSelected, isSelected, int);     //! transverse momentum
DECLARE_SOA_COLUMN(SigBgStatus, sigBgStatus, int);     //! transverse momentum
}

DECLARE_SOA_TABLE(HfCandLcKFs, "AOD", "HFCANDLCKF",
                  kf::Chi2PrimProton, kf::Chi2PrimKaon, kf::Chi2PrimPion,
                  kf::DCAProtonKaon, kf::DCAProtonPion, kf::DCAPionKaon,
                  kf::Chi2geoProtonKaon, kf::Chi2geoProtonPion, kf::Chi2geoPionKaon,
                  kf::Chi2geo, kf::Chi2topo, kf::L, kf::DeltaL, kf::T, kf::DeltaT,
                  kf::MassInv, kf::P, kf::Pt, kf::IsSelected, kf::SigBgStatus
);

DECLARE_SOA_TABLE(HfCandLcLites, "AOD", "HFCANDLCLITE",
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_cand::NProngsContributorsPV,
                  hf_cand::BitmapProngsContributorsPV,
                  hf_cand::Chi2PCA,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::PtProng0,
                  full::PtProng1,
                  full::PtProng2,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand::ImpactParameter2,
                  full::NSigTpcPi0,
                  full::NSigTpcPr0,
                  full::NSigTofPi0,
                  full::NSigTofPr0,
                  full::NSigTpcKa1,
                  full::NSigTofKa1,
                  full::NSigTpcPi2,
                  full::NSigTpcPr2,
                  full::NSigTofPi2,
                  full::NSigTofPr2,
                  full::NSigTpcTofPi0,
                  full::NSigTpcTofPr0,
                  full::NSigTpcTofKa1,
                  full::NSigTpcTofPi2,
                  full::NSigTpcTofPr2,
                  full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::Cpa,
                  full::CpaXY,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::FlagMc,
                  full::OriginMcRec,
                  full::IsCandidateSwapped,
                  full::Channel,
                  full::MassKPi);

DECLARE_SOA_TABLE(HfCollIdLCLite, "AOD", "HFCOLLIDLCLITE",
                  full::CollisionId);

DECLARE_SOA_TABLE(HfCandLcFulls, "AOD", "HFCANDLCFULL",
                  full::CollisionId,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_cand::NProngsContributorsPV,
                  hf_cand::BitmapProngsContributorsPV,
                  hf_cand::XSecondaryVertex,
                  hf_cand::YSecondaryVertex,
                  hf_cand::ZSecondaryVertex,
                  hf_cand::ErrorDecayLength,
                  hf_cand::ErrorDecayLengthXY,
                  hf_cand::Chi2PCA,
                  full::RSecondaryVertex,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::ImpactParameterNormalised0,
                  full::PtProng0,
                  full::PProng0,
                  full::ImpactParameterNormalised1,
                  full::PtProng1,
                  full::PProng1,
                  full::ImpactParameterNormalised2,
                  full::PtProng2,
                  full::PProng2,
                  hf_cand::PxProng0,
                  hf_cand::PyProng0,
                  hf_cand::PzProng0,
                  hf_cand::PxProng1,
                  hf_cand::PyProng1,
                  hf_cand::PzProng1,
                  hf_cand::PxProng2,
                  hf_cand::PyProng2,
                  hf_cand::PzProng2,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand::ImpactParameter2,
                  hf_cand::ErrorImpactParameter0,
                  hf_cand::ErrorImpactParameter1,
                  hf_cand::ErrorImpactParameter2,
                  full::NSigTpcPi0,
                  full::NSigTpcPr0,
                  full::NSigTofPi0,
                  full::NSigTofPr0,
                  full::NSigTpcKa1,
                  full::NSigTofKa1,
                  full::NSigTpcPi2,
                  full::NSigTpcPr2,
                  full::NSigTofPi2,
                  full::NSigTofPr2,
                  full::NSigTpcTofPi0,
                  full::NSigTpcTofPr0,
                  full::NSigTpcTofKa1,
                  full::NSigTpcTofPi2,
                  full::NSigTpcTofPr2,
                  full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::P,
                  full::Cpa,
                  full::CpaXY,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::E,
                  full::FlagMc,
                  full::OriginMcRec,
                  full::IsCandidateSwapped,
                  full::CandidateId,
                  full::Channel,
                  full::MassKPi);

DECLARE_SOA_TABLE(HfCandLcFullEvs, "AOD", "HFCANDLCFULLEV",
                  full::CollisionId,
                  full::McCollisionId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber,
                  full::CentFT0A,
                  full::CentFT0C,
                  full::CentFT0M,
                  full::CentFV0A,
                  full::CentFDDM,
                  full::MultZeqNTracksPV);

DECLARE_SOA_TABLE(HfCandLcFullPs, "AOD", "HFCANDLCFULLP",
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::FlagMc,
                  full::OriginMcGen);



} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorLcToPKPi {
  Produces<o2::aod::HfCandLcFulls> rowCandidateFull;
  Produces<o2::aod::HfCandLcLites> rowCandidateLite;
  Produces<o2::aod::HfCandLcKFs> rowCandidateKF;
  Produces<o2::aod::HfCollIdLCLite> rowCollisionId;
  Produces<o2::aod::HfCandLcFullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandLcFullPs> rowCandidateFullParticles;

  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  Configurable<bool> fillCollIdTable{"fillCollIdTable", false, "Fill a single-column table with collision index"};
  Configurable<bool> keepOnlySignalMc{"keepOnlySignalMc", false, "Fill MC tree only with signal candidates"};
  Configurable<bool> keepOnlyBkg{"keepOnlyBkg", false, "Fill MC tree only with background candidates"};
  Configurable<double> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of candidates to store in the tree"};
  Configurable<float> downSampleBkgPtMax{"downSampleBkgPtMax", 100.f, "Max. pt for background downsampling"};

  HfHelper hfHelper;

  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  using Cents = soa::Join<aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs>;

  template <typename CandType>
  int DetermineSignalBgStatus(const CandType& candidate, int CandFlag) {
    const int flag = candidate.flagMcMatchRec();
    const int origin = candidate.originMcRec();
    const int swapped = candidate.isCandidateSwapped();
    int status{-1}; // 0 bg, 1 prompt, 2 non-prompt, 3 wrong order of prongs, -1 default value (illegal, should not be the case)

    if(std::abs(flag) == (1 << o2::aod::hf_cand_3prong::DecayType::LcToPKPi)) {
      if(swapped == 0) {
        if(CandFlag == 0) {
          if(origin == RecoDecay::OriginType::Prompt) status = 1;
          else if(origin == RecoDecay::OriginType::NonPrompt) status = 2;
        } else {
          status = 3;
        }
      } else {
        if(CandFlag == 1) {
          if(origin == RecoDecay::OriginType::Prompt) status = 1;
          else if(origin == RecoDecay::OriginType::NonPrompt) status = 2;
        } else {
          status = 3;
        }
      }
    } else {
      status = 0;
    }

    if(status == -1) {
      throw std::runtime_error("DetermineSignalBgStatus(): status == -1");
    }

    return status;
  }

  void init(InitContext const&)
  {
    std::array<bool, 6> processes = {doprocessDataNoCentrality, doprocessDataWithCentrality, doprocessMcNoCentralityWithDCAFitterN, doprocessMcWithCentralityWithDCAFitterN, doprocessMcNoCentralityWithKFParticle, doprocessMcWithCentralityWithKFParticle};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }
  }

  /// \brief core function to fill tables in MC
  /// \param collisions Collision table
  /// \param mcCollisions MC collision table
  /// \param candidates Lc->pKpi candidate table
  /// \param particles Generated particle table
  template <bool useCentrality, int reconstructionType, typename Colls, typename CandType>
  void fillTablesMc(Colls const& collisions,
                    aod::McCollisions const&,
                    CandType const& candidates,
                    soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                    TracksWPid const&, aod::BCs const&)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {

      float centFT0A = -1.f;
      float centFT0C = -1.f;
      float centFT0M = -1.f;
      float centFV0A = -1.f;
      float centFDDM = -1.f;
      if constexpr (useCentrality) {
        centFT0A = collision.centFT0A();
        centFT0C = collision.centFT0C();
        centFT0M = collision.centFT0M();
        centFV0A = collision.centFV0A();
        centFDDM = collision.centFDDM();
      }

      rowCandidateFullEvents(
        collision.globalIndex(),
        collision.mcCollisionId(),
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        0,
        collision.bc().runNumber(),
        centFT0A,
        centFT0C,
        centFT0M,
        centFV0A,
        centFDDM,
        collision.multZeqNTracksPV());
    }

    // Filling candidate properties
    if constexpr (reconstructionType == aod::hf_cand::VertexerType::DCAFitter) {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(candidates.size()*2);
      } else {
        rowCandidateFull.reserve(candidates.size()*2);
      }
    } else {
      rowCandidateKF.reserve(candidates.size()*2);
    }
    if (fillCollIdTable) {
      /// save also candidate collision indices
      rowCollisionId.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      auto trackPos1 = candidate.template prong0_as<TracksWPid>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.template prong1_as<TracksWPid>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.template prong2_as<TracksWPid>(); // positive daughter (negative for the antiparticles)
      bool isMcCandidateSignal = std::abs(candidate.flagMcMatchRec()) == (1 << o2::aod::hf_cand_3prong::DecayType::LcToPKPi);
      auto fillTable = [&](int CandFlag,
                           int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY,
                           float FunctionE,
                           float FunctionInvMassKPi) {
        double pseudoRndm = trackPos1.pt() * 1000. - (int64_t)(trackPos1.pt() * 1000);
        if (FunctionSelection >= selectionFlagLc && (/*keep all*/ (!keepOnlySignalMc && !keepOnlyBkg) || /*keep only signal*/ (keepOnlySignalMc && isMcCandidateSignal) || /*keep only background and downsample it*/ (keepOnlyBkg && !isMcCandidateSignal && (candidate.pt() > downSampleBkgPtMax || (pseudoRndm < downSampleBkgFactor && candidate.pt() < downSampleBkgPtMax))))) {
          if constexpr (reconstructionType == aod::hf_cand::VertexerType::DCAFitter) {
            if (fillCandidateLiteTable) {
              rowCandidateLite(
                candidate.posX(),
                candidate.posY(),
                candidate.posZ(),
                candidate.nProngsContributorsPV(),
                candidate.bitmapProngsContributorsPV(),
                candidate.chi2PCA(),
                candidate.decayLength(),
                candidate.decayLengthXY(),
                candidate.ptProng0(),
                candidate.ptProng1(),
                candidate.ptProng2(),
                candidate.impactParameter0(),
                candidate.impactParameter1(),
                candidate.impactParameter2(),
                trackPos1.tpcNSigmaPi(),
                trackPos1.tpcNSigmaPr(),
                trackPos1.tofNSigmaPi(),
                trackPos1.tofNSigmaPr(),
                trackNeg.tpcNSigmaKa(),
                trackNeg.tofNSigmaKa(),
                trackPos2.tpcNSigmaPi(),
                trackPos2.tpcNSigmaPr(),
                trackPos2.tofNSigmaPi(),
                trackPos2.tofNSigmaPr(),
                trackPos1.tpcTofNSigmaPi(),
                trackPos1.tpcTofNSigmaPr(),
                trackNeg.tpcTofNSigmaKa(),
                trackPos2.tpcTofNSigmaPi(),
                trackPos2.tpcTofNSigmaPr(),
                1 << CandFlag,
                FunctionInvMass,
                candidate.pt(),
                candidate.cpa(),
                candidate.cpaXY(),
                FunctionCt,
                candidate.eta(),
                candidate.phi(),
                FunctionY,
                candidate.flagMcMatchRec(),
                candidate.originMcRec(),
                candidate.isCandidateSwapped(),
                candidate.flagMcDecayChanRec(),
                FunctionInvMassKPi);

              if (fillCollIdTable) {
                /// save also candidate collision indices
                rowCollisionId(candidate.collisionId());
              }

            } else {
              rowCandidateFull(
                candidate.collisionId(),
                candidate.posX(),
                candidate.posY(),
                candidate.posZ(),
                candidate.nProngsContributorsPV(),
                candidate.bitmapProngsContributorsPV(),
                candidate.xSecondaryVertex(),
                candidate.ySecondaryVertex(),
                candidate.zSecondaryVertex(),
                candidate.errorDecayLength(),
                candidate.errorDecayLengthXY(),
                candidate.chi2PCA(),
                candidate.rSecondaryVertex(),
                candidate.decayLength(),
                candidate.decayLengthXY(),
                candidate.decayLengthNormalised(),
                candidate.decayLengthXYNormalised(),
                candidate.impactParameterNormalised0(),
                candidate.ptProng0(),
                RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
                candidate.impactParameterNormalised1(),
                candidate.ptProng1(),
                RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
                candidate.impactParameterNormalised2(),
                candidate.ptProng2(),
                RecoDecay::p(candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()),
                candidate.pxProng0(),
                candidate.pyProng0(),
                candidate.pzProng0(),
                candidate.pxProng1(),
                candidate.pyProng1(),
                candidate.pzProng1(),
                candidate.pxProng2(),
                candidate.pyProng2(),
                candidate.pzProng2(),
                candidate.impactParameter0(),
                candidate.impactParameter1(),
                candidate.impactParameter2(),
                candidate.errorImpactParameter0(),
                candidate.errorImpactParameter1(),
                candidate.errorImpactParameter2(),
                trackPos1.tpcNSigmaPi(),
                trackPos1.tpcNSigmaPr(),
                trackPos1.tofNSigmaPi(),
                trackPos1.tofNSigmaPr(),
                trackNeg.tpcNSigmaKa(),
                trackNeg.tofNSigmaKa(),
                trackPos2.tpcNSigmaPi(),
                trackPos2.tpcNSigmaPr(),
                trackPos2.tofNSigmaPi(),
                trackPos2.tofNSigmaPr(),
                trackPos1.tpcTofNSigmaPi(),
                trackPos1.tpcTofNSigmaPr(),
                trackNeg.tpcTofNSigmaKa(),
                trackPos2.tpcTofNSigmaPi(),
                trackPos2.tpcTofNSigmaPr(),
                1 << CandFlag,
                FunctionInvMass,
                candidate.pt(),
                candidate.p(),
                candidate.cpa(),
                candidate.cpaXY(),
                FunctionCt,
                candidate.eta(),
                candidate.phi(),
                FunctionY,
                FunctionE,
                candidate.flagMcMatchRec(),
                candidate.originMcRec(),
                candidate.isCandidateSwapped(),
                candidate.globalIndex(),
                candidate.flagMcDecayChanRec(),
                FunctionInvMassKPi);
            }
          } else {
            const float chi2prim_proton = CandFlag == 0 ? candidate.kfChi2PrimProng0() : candidate.kfChi2PrimProng2();
            const float chi2prim_kaon = candidate.kfChi2PrimProng1();
            const float chi2prim_pion = CandFlag == 0 ? candidate.kfChi2PrimProng2() : candidate.kfChi2PrimProng0();
            const float dca_proton_kaon = CandFlag == 0 ? candidate.kfDCAProng0Prong1() : candidate.kfDCAProng1Prong2();
            const float dca_proton_pion = candidate.kfDCAProng0Prong2();
            const float dca_pion_kaon = CandFlag == 0 ? candidate.kfDCAProng1Prong2() : candidate.kfDCAProng0Prong1();
            const float chi2geo_proton_kaon = CandFlag == 0 ? candidate.kfChi2geoProng0Prong1() : candidate.kfChi2geoProng1Prong2();
            const float chi2geo_proton_pion = candidate.kfChi2geoProng0Prong2();
            const float chi2geo_pion_kaon = CandFlag == 0 ? candidate.kfChi2geoProng1Prong2() : candidate.kfChi2geoProng0Prong1();
            const float chi2geo = candidate.kfChi2geo();
            const float chi2topo = candidate.kfChi2topo();
            const float l = candidate.kfL();
            const float dl = candidate.kfDeltaL();
            const float pt = std::sqrt(candidate.kfPx()*candidate.kfPx() + candidate.kfPy()*candidate.kfPy());
            const float p = std::sqrt(pt*pt + candidate.kfPz()*candidate.kfPz());
            const float T = l * MassLambdaCPlus / LightSpeedCm2PS / p;
            const float deltaT = dl * MassLambdaCPlus / LightSpeedCm2PS / p;
            const float mass = CandFlag == 0 ? candidate.kfMassPKPi() : candidate.kfMassPiKP();
            const int selectedStatus = CandFlag == 0 ? candidate.isSelLcToPKPi() : candidate.isSelLcToPiKP();
            const int sigbgstatus = DetermineSignalBgStatus(candidate, CandFlag);
            rowCandidateKF(
              chi2prim_proton, chi2prim_kaon, chi2prim_pion,
              dca_proton_kaon, dca_proton_pion, dca_pion_kaon,
              chi2geo_proton_kaon, chi2geo_proton_pion, chi2geo_pion_kaon,
              chi2geo, chi2topo, l, dl, T, deltaT,
              mass, p, pt, selectedStatus, sigbgstatus
            );

            if (fillCollIdTable) {
              /// save also candidate collision indices
              rowCollisionId(candidate.collisionId());
            }
          }
        }
      };

      fillTable(0, candidate.isSelLcToPKPi(), hfHelper.invMassLcToPKPi(candidate), hfHelper.ctLc(candidate), hfHelper.yLc(candidate), hfHelper.eLc(candidate), hfHelper.invMassKPiPairLcToPKPi(candidate));
      fillTable(1, candidate.isSelLcToPiKP(), hfHelper.invMassLcToPiKP(candidate), hfHelper.ctLc(candidate), hfHelper.yLc(candidate), hfHelper.eLc(candidate), hfHelper.invMassKPiPairLcToPiKP(candidate));
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (const auto& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi) {
        rowCandidateFullParticles(
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus),
          particle.flagMcMatchGen(),
          particle.originMcGen());
      }
    }
  }

  /// \brief process function for MC w/o centrality
  /// \param collisions Collision table w/o join of the centrality table
  /// \param mcCollisions MC collision table
  /// \param candidates Lc->pKpi candidate table
  /// \param particles Generated particle table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processMcNoCentralityWithDCAFitterN(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::PVMultZeqs> const& collisions,
                             aod::McCollisions const& mcCollisions,
                             soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc> const& candidates,
                             soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                             TracksWPid const& tracks, aod::BCs const& bcs)
  {
    fillTablesMc<false, aod::hf_cand::VertexerType::DCAFitter>(collisions, mcCollisions, candidates, particles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processMcNoCentralityWithDCAFitterN, "Process MC tree writer w/o centrality with DCAFitterN", false);

  /// \brief process function for MC with centrality
  /// \param collisions Collision table with join of the centrality table
  /// \param mcCollisions MC collision table
  /// \param candidates Lc->pKpi candidate table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processMcWithCentralityWithDCAFitterN(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::PVMultZeqs, Cents> const& collisions,
                               aod::McCollisions const& mcCollisions,
                               soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc> const& candidates,
                               soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                               TracksWPid const& tracks, aod::BCs const& bcs)
  {
    fillTablesMc<true, aod::hf_cand::VertexerType::DCAFitter>(collisions, mcCollisions, candidates, particles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processMcWithCentralityWithDCAFitterN, "Process MC tree writer with centrality with DCAFitterN", false);

  /// \brief process function for MC w/o centrality
  /// \param collisions Collision table w/o join of the centrality table
  /// \param mcCollisions MC collision table
  /// \param candidates Lc->pKpi candidate table
  /// \param particles Generated particle table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processMcNoCentralityWithKFParticle(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::PVMultZeqs> const& collisions,
                             aod::McCollisions const& mcCollisions,
                             soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc, aod::HfCand3ProngKF> const& candidates,
                             soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                             TracksWPid const& tracks, aod::BCs const& bcs)
  {
    fillTablesMc<false, aod::hf_cand::VertexerType::KfParticle>(collisions, mcCollisions, candidates, particles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processMcNoCentralityWithKFParticle, "Process MC tree writer w/o centrality with KFParticle", false);

  /// \brief process function for MC with centrality
  /// \param collisions Collision table with join of the centrality table
  /// \param mcCollisions MC collision table
  /// \param candidates Lc->pKpi candidate table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processMcWithCentralityWithKFParticle(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::PVMultZeqs, Cents> const& collisions,
                               aod::McCollisions const& mcCollisions,
                               soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc, aod::HfCand3ProngKF> const& candidates,
                               soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                               TracksWPid const& tracks, aod::BCs const& bcs)
  {
    fillTablesMc<true, aod::hf_cand::VertexerType::KfParticle>(collisions, mcCollisions, candidates, particles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processMcWithCentralityWithKFParticle, "Process MC tree writer with centrality with KFParticle", false);

  /// \brief core function to fill tables in data
  /// \param collisions Collision table
  /// \param candidates Lc->pKpi candidate table
  template <bool useCentrality, typename Colls>
  void fillTablesData(Colls const& collisions,
                      soa::Join<aod::HfCand3Prong, aod::HfSelLc> const& candidates,
                      TracksWPid const&, aod::BCs const&)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {

      float centFT0A = -1.f;
      float centFT0C = -1.f;
      float centFT0M = -1.f;
      float centFV0A = -1.f;
      float centFDDM = -1.f;
      if constexpr (useCentrality) {
        centFT0A = collision.centFT0A();
        centFT0C = collision.centFT0C();
        centFT0M = collision.centFT0M();
        centFV0A = collision.centFV0A();
        centFDDM = collision.centFDDM();
      }

      rowCandidateFullEvents(
        collision.globalIndex(),
        -1,
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        0,
        collision.bc().runNumber(),
        centFT0A,
        centFT0C,
        centFT0M,
        centFV0A,
        centFDDM,
        collision.multZeqNTracksPV());
    }

    // Filling candidate properties
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidates.size());
    } else {
      rowCandidateFull.reserve(candidates.size());
    }
    if (fillCollIdTable) {
      /// save also candidate collision indices
      rowCollisionId.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      auto trackPos1 = candidate.prong0_as<TracksWPid>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<TracksWPid>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<TracksWPid>(); // positive daughter (negative for the antiparticles)
      auto fillTable = [&](int CandFlag,
                           int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY,
                           float FunctionE,
                           float FunctionInvMassKPi) {
        double pseudoRndm = trackPos1.pt() * 1000. - (int64_t)(trackPos1.pt() * 1000);
        if (FunctionSelection >= 1 && (candidate.pt() > downSampleBkgPtMax || (pseudoRndm < downSampleBkgFactor && candidate.pt() < downSampleBkgPtMax))) {
          if (fillCandidateLiteTable) {
            rowCandidateLite(
              candidate.posX(),
              candidate.posY(),
              candidate.posZ(),
              candidate.nProngsContributorsPV(),
              candidate.bitmapProngsContributorsPV(),
              // candidate.errorDecayLength(),
              // candidate.errorDecayLengthXY(),
              candidate.chi2PCA(),
              candidate.decayLength(),
              candidate.decayLengthXY(),
              // candidate.decayLengthNormalised(),
              // candidate.decayLengthXYNormalised(),
              // candidate.impactParameterNormalised0(),
              candidate.ptProng0(),
              // candidate.impactParameterNormalised1(),
              candidate.ptProng1(),
              // candidate.impactParameterNormalised2(),
              candidate.ptProng2(),
              candidate.impactParameter0(),
              candidate.impactParameter1(),
              candidate.impactParameter2(),
              // candidate.errorImpactParameter0(),
              // candidate.errorImpactParameter1(),
              // candidate.errorImpactParameter2(),
              trackPos1.tpcNSigmaPi(),
              trackPos1.tpcNSigmaPr(),
              trackPos1.tofNSigmaPi(),
              trackPos1.tofNSigmaPr(),
              trackNeg.tpcNSigmaKa(),
              trackNeg.tofNSigmaKa(),
              trackPos2.tpcNSigmaPi(),
              trackPos2.tpcNSigmaPr(),
              trackPos2.tofNSigmaPi(),
              trackPos2.tofNSigmaPr(),
              trackPos1.tpcTofNSigmaPi(),
              trackPos1.tpcTofNSigmaPr(),
              trackNeg.tpcTofNSigmaKa(),
              trackPos2.tpcTofNSigmaPi(),
              trackPos2.tpcTofNSigmaPr(),
              1 << CandFlag,
              FunctionInvMass,
              candidate.pt(),
              candidate.cpa(),
              candidate.cpaXY(),
              FunctionCt,
              candidate.eta(),
              candidate.phi(),
              FunctionY,
              0.,
              0.,
              0.,
              -1,
              FunctionInvMassKPi);
            // candidate.globalIndex());

            if (fillCollIdTable) {
              /// save also candidate collision indices
              rowCollisionId(candidate.collisionId());
            }

          } else {
            rowCandidateFull(
              candidate.collisionId(),
              candidate.posX(),
              candidate.posY(),
              candidate.posZ(),
              candidate.nProngsContributorsPV(),
              candidate.bitmapProngsContributorsPV(),
              candidate.xSecondaryVertex(),
              candidate.ySecondaryVertex(),
              candidate.zSecondaryVertex(),
              candidate.errorDecayLength(),
              candidate.errorDecayLengthXY(),
              candidate.chi2PCA(),
              candidate.rSecondaryVertex(),
              candidate.decayLength(),
              candidate.decayLengthXY(),
              candidate.decayLengthNormalised(),
              candidate.decayLengthXYNormalised(),
              candidate.impactParameterNormalised0(),
              candidate.ptProng0(),
              RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
              candidate.impactParameterNormalised1(),
              candidate.ptProng1(),
              RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
              candidate.impactParameterNormalised2(),
              candidate.ptProng2(),
              RecoDecay::p(candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()),
              candidate.pxProng0(),
              candidate.pyProng0(),
              candidate.pzProng0(),
              candidate.pxProng1(),
              candidate.pyProng1(),
              candidate.pzProng1(),
              candidate.pxProng2(),
              candidate.pyProng2(),
              candidate.pzProng2(),
              candidate.impactParameter0(),
              candidate.impactParameter1(),
              candidate.impactParameter2(),
              candidate.errorImpactParameter0(),
              candidate.errorImpactParameter1(),
              candidate.errorImpactParameter2(),
              trackPos1.tpcNSigmaPi(),
              trackPos1.tpcNSigmaPr(),
              trackPos1.tofNSigmaPi(),
              trackPos1.tofNSigmaPr(),
              trackNeg.tpcNSigmaKa(),
              trackNeg.tofNSigmaKa(),
              trackPos2.tpcNSigmaPi(),
              trackPos2.tpcNSigmaPr(),
              trackPos2.tofNSigmaPi(),
              trackPos2.tofNSigmaPr(),
              trackPos1.tpcTofNSigmaPi(),
              trackPos1.tpcTofNSigmaPr(),
              trackNeg.tpcTofNSigmaKa(),
              trackPos2.tpcTofNSigmaPi(),
              trackPos2.tpcTofNSigmaPr(),
              1 << CandFlag,
              FunctionInvMass,
              candidate.pt(),
              candidate.p(),
              candidate.cpa(),
              candidate.cpaXY(),
              FunctionCt,
              candidate.eta(),
              candidate.phi(),
              FunctionY,
              FunctionE,
              0.,
              0.,
              0.,
              candidate.globalIndex(),
              -1,
              FunctionInvMassKPi);
          }
        }
      };

      fillTable(0, candidate.isSelLcToPKPi(), hfHelper.invMassLcToPKPi(candidate), hfHelper.ctLc(candidate), hfHelper.yLc(candidate), hfHelper.eLc(candidate), hfHelper.invMassKPiPairLcToPKPi(candidate));
      fillTable(1, candidate.isSelLcToPiKP(), hfHelper.invMassLcToPiKP(candidate), hfHelper.ctLc(candidate), hfHelper.yLc(candidate), hfHelper.eLc(candidate), hfHelper.invMassKPiPairLcToPiKP(candidate));
    }
  }

  /// \brief process function for data w/o centrality
  /// \param collisions Collision table w/o join of the centrality table
  /// \param candidates Lc->pKpi candidate table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processDataNoCentrality(soa::Join<aod::Collisions, aod::PVMultZeqs> const& collisions,
                               soa::Join<aod::HfCand3Prong, aod::HfSelLc> const& candidates,
                               TracksWPid const& tracks, aod::BCs const& bcs)
  {
    fillTablesData<false>(collisions, candidates, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processDataNoCentrality, "Process data tree writer w/o centrality", false);

  /// \brief process function for data with centrality
  /// \param collisions Collision table with join of the centrality table
  /// \param candidates Lc->pKpi candidate table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processDataWithCentrality(soa::Join<aod::Collisions, aod::PVMultZeqs, Cents> const& collisions,
                                 soa::Join<aod::HfCand3Prong, aod::HfSelLc> const& candidates,
                                 TracksWPid const& tracks, aod::BCs const& bcs)
  {
    fillTablesData<true>(collisions, candidates, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processDataWithCentrality, "Process data tree writer with centrality", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorLcToPKPi>(cfgc));
  return workflow;
}
