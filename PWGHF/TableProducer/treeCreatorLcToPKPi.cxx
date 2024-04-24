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
} // namespace full

DECLARE_SOA_TABLE(HfCandLcLites, "AOD", "HFCANDLCLITE",
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_cand::NProngsContributorsPV,
                  // hf_cand::ErrorDecayLength,
                  // hf_cand::ErrorDecayLengthXY,
                  hf_cand::Chi2PCA,
                  full::DecayLength,
                  full::DecayLengthXY,
                  // full::DecayLengthNormalised,
                  // full::DecayLengthXYNormalised,
                  // full::ImpactParameterNormalised0,
                  full::PtProng0,
                  // full::ImpactParameterNormalised1,
                  full::PtProng1,
                  // full::ImpactParameterNormalised2,
                  full::PtProng2,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand::ImpactParameter2,
                  // hf_cand::ErrorImpactParameter0,
                  // hf_cand::ErrorImpactParameter1,
                  // hf_cand::ErrorImpactParameter2,
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
                  full::IsCandidateSwapped);

DECLARE_SOA_TABLE(HfCollIdLCLite, "AOD", "HFCOLLIDLCLITE",
                  full::CollisionId);

DECLARE_SOA_TABLE(HfCandLcFulls, "AOD", "HFCANDLCFULL",
                  full::CollisionId,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_cand::NProngsContributorsPV,
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
                  full::CandidateId);

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
                  full::McCollisionId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::FlagMc,
                  full::OriginMcGen,
                  full::McParticleId);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorLcToPKPi {
  Produces<o2::aod::HfCandLcFulls> rowCandidateFull;
  Produces<o2::aod::HfCandLcLites> rowCandidateLite;
  Produces<o2::aod::HfCollIdLCLite> rowCollisionId;
  Produces<o2::aod::HfCandLcFullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandLcFullPs> rowCandidateFullParticles;

  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  Configurable<bool> fillCollIdTable{"fillCollIdTable", false, "Fill a single-column table with collision index"};
  Configurable<bool> keepOnlySignalMc{"keepOnlySignalMc", false, "Fill MC tree only with signal candidates"};
  Configurable<bool> keepOnlyBkg{"keepOnlyBkg", false, "Fill MC tree only with background candidates"};
  Configurable<double> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of candidates to store in the tree"};
  Configurable<float> downSampleBkgPtMax{"downSampleBkgPtMax", 100.f, "Max. pt for background downsampling"};

  HfHelper hfHelper;

  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  using Cents = soa::Join<aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs>;

  void init(InitContext const&)
  {
  }

  void processMc(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::PVMultZeqs, Cents> const& collisions,
                 aod::McCollisions const&,
                 soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                 TracksWPid const&, aod::BCs const&)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      rowCandidateFullEvents(
        collision.globalIndex(),
        collision.mcCollisionId(),
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        0,
        collision.bc().runNumber(),
        collision.centFT0A(),
        collision.centFT0C(),
        collision.centFT0M(),
        collision.centFV0A(),
        collision.centFDDM(),
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
      bool isMcCandidateSignal = std::abs(candidate.flagMcMatchRec()) == (1 << o2::aod::hf_cand_3prong::DecayType::LcToPKPi);
      auto fillTable = [&](int CandFlag,
                           int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY,
                           float FunctionE) {
        double pseudoRndm = trackPos1.pt() * 1000. - (int64_t)(trackPos1.pt() * 1000);
        if (FunctionSelection >= 1 && (/*keep all*/ (!keepOnlySignalMc && !keepOnlyBkg) || /*keep only signal*/ (keepOnlySignalMc && isMcCandidateSignal) || /*keep only background and downsample it*/ (keepOnlyBkg && !isMcCandidateSignal && (candidate.pt() > downSampleBkgPtMax || (pseudoRndm < downSampleBkgFactor && candidate.pt() < downSampleBkgPtMax))))) {
          if (fillCandidateLiteTable) {
            rowCandidateLite(
              candidate.posX(),
              candidate.posY(),
              candidate.posZ(),
              candidate.nProngsContributorsPV(),
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
              candidate.flagMcMatchRec(),
              candidate.originMcRec(),
              candidate.isCandidateSwapped());
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
              candidate.globalIndex());
          }
        }
      };

      fillTable(0, candidate.isSelLcToPKPi(), hfHelper.invMassLcToPKPi(candidate), hfHelper.ctLc(candidate), hfHelper.yLc(candidate), hfHelper.eLc(candidate));
      fillTable(1, candidate.isSelLcToPiKP(), hfHelper.invMassLcToPiKP(candidate), hfHelper.ctLc(candidate), hfHelper.yLc(candidate), hfHelper.eLc(candidate));
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (const auto& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi) {
        rowCandidateFullParticles(
          particle.mcCollisionId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus),
          particle.flagMcMatchGen(),
          particle.originMcGen(),
          particle.globalIndex());
      }
    }
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processMc, "Process MC tree writer", true);

  void processData(soa::Join<aod::Collisions, aod::PVMultZeqs, Cents> const& collisions,
                   soa::Join<aod::HfCand3Prong, aod::HfSelLc> const& candidates,
                   TracksWPid const&, aod::BCs const&)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      rowCandidateFullEvents(
        collision.globalIndex(),
        -1,
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        0,
        collision.bc().runNumber(),
        collision.centFT0A(),
        collision.centFT0C(),
        collision.centFT0M(),
        collision.centFV0A(),
        collision.centFDDM(),
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
                           float FunctionE) {
        double pseudoRndm = trackPos1.pt() * 1000. - (int64_t)(trackPos1.pt() * 1000);
        if (FunctionSelection >= 1 && (candidate.pt() > downSampleBkgPtMax || (pseudoRndm < downSampleBkgFactor && candidate.pt() < downSampleBkgPtMax))) {
          if (fillCandidateLiteTable) {
            rowCandidateLite(
              candidate.posX(),
              candidate.posY(),
              candidate.posZ(),
              candidate.nProngsContributorsPV(),
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
              0.);
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
              candidate.globalIndex());
          }
        }
      };

      fillTable(0, candidate.isSelLcToPKPi(), hfHelper.invMassLcToPKPi(candidate), hfHelper.ctLc(candidate), hfHelper.yLc(candidate), hfHelper.eLc(candidate));
      fillTable(1, candidate.isSelLcToPiKP(), hfHelper.invMassLcToPiKP(candidate), hfHelper.ctLc(candidate), hfHelper.yLc(candidate), hfHelper.eLc(candidate));
    }
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processData, "Process data tree writer", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorLcToPKPi>(cfgc));
  return workflow;
}
