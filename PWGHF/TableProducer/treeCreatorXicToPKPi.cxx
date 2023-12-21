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

/// \file treeCreatorXicToPKPi.cxx
/// \brief Writer of XiC -> pKpi candidates in the form of flat tables to be stored in TTrees.
/// inspired from file treeCreatorLcToPKPi.cxx and to treeCreatorDplusToPiKPi.cxx

/// \author Himanshu Sharma <himanshu.sharma@cern.ch>, INFN Padova
/// \author Cristina Terrevoli <cristina.terrevoli@cern.ch>, INFN Bari

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

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
DECLARE_SOA_COLUMN(NSigTpcKa0, nSigTpcKa0, float);
DECLARE_SOA_COLUMN(NSigTpcPr0, nSigTpcPr0, float);
DECLARE_SOA_COLUMN(NSigTofPi0, nSigTofPi0, float);
DECLARE_SOA_COLUMN(NSigTofKa0, nSigTofKa0, float);
DECLARE_SOA_COLUMN(NSigTofPr0, nSigTofPr0, float);
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);
DECLARE_SOA_COLUMN(NSigTpcKa1, nSigTpcKa1, float);
DECLARE_SOA_COLUMN(NSigTpcPr1, nSigTpcPr1, float);
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);
DECLARE_SOA_COLUMN(NSigTofKa1, nSigTofKa1, float);
DECLARE_SOA_COLUMN(NSigTofPr1, nSigTofPr1, float);
DECLARE_SOA_COLUMN(NSigTpcPi2, nSigTpcPi2, float);
DECLARE_SOA_COLUMN(NSigTpcKa2, nSigTpcKa2, float);
DECLARE_SOA_COLUMN(NSigTpcPr2, nSigTpcPr2, float);
DECLARE_SOA_COLUMN(NSigTofPi2, nSigTofPi2, float);
DECLARE_SOA_COLUMN(NSigTofKa2, nSigTofKa2, float);
DECLARE_SOA_COLUMN(NSigTofPr2, nSigTofPr2, float);
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
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

DECLARE_SOA_TABLE(HfCandXicLites, "AOD", "HFCANDXICLITE",
                  hf_cand::Chi2PCA,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::PtProng0,
                  full::PtProng1,
                  full::PtProng2,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand::ImpactParameter2,
                  full::NSigTpcPi0,
                  full::NSigTpcKa0,
                  full::NSigTpcPr0,
                  full::NSigTofPi0,
                  full::NSigTofKa0,
                  full::NSigTofPr0,
                  full::NSigTpcPi1,
                  full::NSigTpcKa1,
                  full::NSigTpcPr1,
                  full::NSigTofPi1,
                  full::NSigTofKa1,
                  full::NSigTofPr1,
                  full::NSigTpcPi2,
                  full::NSigTpcKa2,
                  full::NSigTpcPr2,
                  full::NSigTofPi2,
                  full::NSigTofKa2,
                  full::NSigTofPr2,
                  full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::Cpa,
                  full::CpaXY,
                  full::Eta,
                  full::Phi,
                  full::FlagMc,
                  full::OriginMcRec,
                  full::IsCandidateSwapped)

DECLARE_SOA_TABLE(HfCandXicFulls, "AOD", "HFCANDXICFULL",
                  full::CollisionId,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
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
                  full::NSigTpcKa0,
                  full::NSigTpcPr0,
                  full::NSigTofPi0,
                  full::NSigTofKa0,
                  full::NSigTofPr0,
                  full::NSigTpcPi1,
                  full::NSigTpcKa1,
                  full::NSigTpcPr1,
                  full::NSigTofPi1,
                  full::NSigTofKa1,
                  full::NSigTofPr1,
                  full::NSigTpcPi2,
                  full::NSigTpcKa2,
                  full::NSigTpcPr2,
                  full::NSigTofPi2,
                  full::NSigTofKa2,
                  full::NSigTofPr2,
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

DECLARE_SOA_TABLE(HfCandXicFullEvs, "AOD", "HFCANDXICFULLEV",
                  full::CollisionId,
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandXicFullPs, "AOD", "HFCANDXICFULLP",
                  full::CollisionId,
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::FlagMc,
                  full::OriginMcGen,
                  full::CandidateId);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorXicToPKPi {
  Produces<o2::aod::HfCandXicFulls> rowCandidateFull;
  Produces<o2::aod::HfCandXicFullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandXicFullPs> rowCandidateFullParticles;
  Produces<o2::aod::HfCandXicLites> rowCandidateLite;

  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  // parameters for production of training samples
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};
  Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background for ML trainings"};

  HfHelper hfHelper;

  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::TracksPidKa, aod::TracksPidPr>;

  void init(InitContext const&)
  {
  }

  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const& mcCollisions,
                 soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelXicToPKPi> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                 TracksWPid const& tracks)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      rowCandidateFullEvents(
        collision.globalIndex(),
        collision.bcId(),
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        0,
        1);
    }

    // Filling candidate properties
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidates.size());
    } else {
      rowCandidateFull.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      auto trackPos1 = candidate.prong0_as<TracksWPid>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<TracksWPid>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<TracksWPid>(); // positive daughter (negative for the antiparticles)
      if (fillOnlyBackground && downSampleBkgFactor < 1.) {
        double pseudoRndm = trackPos1.pt() * 1000. - (int64_t)(trackPos1.pt() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
      }
      auto fillTable = [&](int CandFlag,
                           int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY,
                           float FunctionE) {
        if (FunctionSelection >= 1) {
          if (fillCandidateLiteTable) {

            rowCandidateLite(
              candidate.chi2PCA(),
              candidate.decayLength(),
              candidate.decayLengthXY(),
              candidate.decayLengthNormalised(),
              candidate.decayLengthXYNormalised(),
              candidate.ptProng0(),
              candidate.ptProng1(),
              candidate.ptProng2(),
              candidate.impactParameter0(),
              candidate.impactParameter1(),
              candidate.impactParameter2(),
              trackPos1.tpcNSigmaPi(),
              trackPos1.tpcNSigmaKa(),
              trackPos1.tpcNSigmaPr(),
              trackPos1.tofNSigmaPi(),
              trackPos1.tofNSigmaKa(),
              trackPos1.tofNSigmaPr(),
              trackNeg.tpcNSigmaPi(),
              trackNeg.tpcNSigmaKa(),
              trackNeg.tpcNSigmaPr(),
              trackNeg.tofNSigmaPi(),
              trackNeg.tofNSigmaKa(),
              trackNeg.tofNSigmaPr(),
              trackPos2.tpcNSigmaPi(),
              trackPos2.tpcNSigmaKa(),
              trackPos2.tpcNSigmaPr(),
              trackPos2.tofNSigmaPi(),
              trackPos2.tofNSigmaKa(),
              trackPos2.tofNSigmaPr(),
              1 << CandFlag,
              FunctionInvMass,
              candidate.pt(),
              candidate.cpa(),
              candidate.cpaXY(),
              candidate.eta(),
              candidate.phi(),
              candidate.flagMcMatchRec(),
              candidate.originMcRec(),
              candidate.isCandidateSwapped());
          } else {
            rowCandidateFull(
              candidate.collisionId(),
              candidate.posX(),
              candidate.posY(),
              candidate.posZ(),
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
              trackPos1.tpcNSigmaKa(),
              trackPos1.tpcNSigmaPr(),
              trackPos1.tofNSigmaPi(),
              trackPos1.tofNSigmaKa(),
              trackPos1.tofNSigmaPr(),
              trackNeg.tpcNSigmaPi(),
              trackNeg.tpcNSigmaKa(),
              trackNeg.tpcNSigmaPr(),
              trackNeg.tofNSigmaPi(),
              trackNeg.tofNSigmaKa(),
              trackNeg.tofNSigmaPr(),
              trackPos2.tpcNSigmaPi(),
              trackPos2.tpcNSigmaKa(),
              trackPos2.tpcNSigmaPr(),
              trackPos2.tofNSigmaPi(),
              trackPos2.tofNSigmaKa(),
              trackPos2.tofNSigmaPr(),
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

      fillTable(0, candidate.isSelXicToPKPi(), hfHelper.invMassXicToPKPi(candidate), hfHelper.ctXic(candidate), hfHelper.yXic(candidate), hfHelper.eXic(candidate));
      fillTable(1, candidate.isSelXicToPiKP(), hfHelper.invMassXicToPiKP(candidate), hfHelper.ctXic(candidate), hfHelper.yXic(candidate), hfHelper.eXic(candidate));
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (const auto& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::XicToPKPi) {
        rowCandidateFullParticles(
          particle.mcCollision().globalIndex(),
          particle.mcCollision().bcId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::constants::physics::MassXiCPlus),
          particle.flagMcMatchGen(),
          particle.originMcGen(),
          particle.globalIndex());
      }
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToPKPi, processMc, "Process MC tree writer", true);

  void processData(aod::Collisions const& collisions,
                   soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi> const& candidates,
                   TracksWPid const& tracks)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      rowCandidateFullEvents(
        collision.globalIndex(),
        collision.bcId(),
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        0,
        1);
    }

    // Filling candidate properties
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidates.size());
    } else {
      rowCandidateFull.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      auto trackPos1 = candidate.prong0_as<TracksWPid>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<TracksWPid>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<TracksWPid>(); // positive daughter (negative for the antiparticles)
      if (downSampleBkgFactor < 1.) {
        double pseudoRndm = trackPos1.pt() * 1000. - (int64_t)(trackPos1.pt() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
      }
      auto fillTable = [&](int CandFlag,
                           int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY,
                           float FunctionE) {
        if (FunctionSelection >= 1) {
          if (fillCandidateLiteTable) {
            rowCandidateLite(
              candidate.chi2PCA(),
              candidate.decayLength(),
              candidate.decayLengthXY(),
              candidate.decayLengthNormalised(),
              candidate.decayLengthXYNormalised(),
              candidate.ptProng0(),
              candidate.ptProng1(),
              candidate.ptProng2(),
              candidate.impactParameter0(),
              candidate.impactParameter1(),
              candidate.impactParameter2(),
              trackPos1.tpcNSigmaPi(),
              trackPos1.tpcNSigmaKa(),
              trackPos1.tpcNSigmaPr(),
              trackPos1.tofNSigmaPi(),
              trackPos1.tofNSigmaKa(),
              trackPos1.tofNSigmaPr(),
              trackNeg.tpcNSigmaPi(),
              trackNeg.tpcNSigmaKa(),
              trackNeg.tpcNSigmaPr(),
              trackNeg.tofNSigmaPi(),
              trackNeg.tofNSigmaKa(),
              trackNeg.tofNSigmaPr(),
              trackPos2.tpcNSigmaPi(),
              trackPos2.tpcNSigmaKa(),
              trackPos2.tpcNSigmaPr(),
              trackPos2.tofNSigmaPi(),
              trackPos2.tofNSigmaKa(),
              trackPos2.tofNSigmaPr(),
              1 << CandFlag,
              FunctionInvMass,
              candidate.pt(),
              candidate.cpa(),
              candidate.cpaXY(),
              candidate.eta(),
              candidate.phi(),
              0.,
              0.,
              0.);
          } else {
            rowCandidateFull(
              candidate.collisionId(),
              candidate.posX(),
              candidate.posY(),
              candidate.posZ(),
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
              trackPos1.tpcNSigmaKa(),
              trackPos1.tpcNSigmaPr(),
              trackPos1.tofNSigmaPi(),
              trackPos1.tofNSigmaKa(),
              trackPos1.tofNSigmaPr(),
              trackNeg.tpcNSigmaPi(),
              trackNeg.tpcNSigmaKa(),
              trackNeg.tpcNSigmaPr(),
              trackNeg.tofNSigmaPi(),
              trackNeg.tofNSigmaKa(),
              trackNeg.tofNSigmaPr(),
              trackPos2.tpcNSigmaPi(),
              trackPos2.tpcNSigmaKa(),
              trackPos2.tpcNSigmaPr(),
              trackPos2.tofNSigmaPi(),
              trackPos2.tofNSigmaKa(),
              trackPos2.tofNSigmaPr(),
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

      fillTable(0, candidate.isSelXicToPKPi(), hfHelper.invMassXicToPKPi(candidate), hfHelper.ctXic(candidate), hfHelper.yXic(candidate), hfHelper.eXic(candidate));
      fillTable(1, candidate.isSelXicToPiKP(), hfHelper.invMassXicToPiKP(candidate), hfHelper.ctXic(candidate), hfHelper.yXic(candidate), hfHelper.eXic(candidate));
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToPKPi, processData, "Process data tree writer", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorXicToPKPi>(cfgc));
  return workflow;
}
