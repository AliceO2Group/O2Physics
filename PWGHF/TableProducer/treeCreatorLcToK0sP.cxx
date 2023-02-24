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

/// \file treeCreatorLcToK0sP.cxx
/// \brief Writer of the cascade candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
///        Modified version of treeCreatorLcToPKPi.cxx
///
/// \author Daniel Samitz <daniel.samitz@cern.ch>

#include "Common/Core/trackUtilities.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(ImpactParameterXY0, impactParameterXY0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(ImpactParameterXY1, impactParameterXY1, float);
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(NSigmaTPCPr0, nSigmaTPCPr0, float);
DECLARE_SOA_COLUMN(NSigmaTOFPr0, nSigmaTOFPr0, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(CPA, cpa, float);
DECLARE_SOA_COLUMN(CPAXY, cpaXY, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(PtV0Pos, ptV0Pos, float);
DECLARE_SOA_COLUMN(PtV0Neg, ptV0Neg, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(V0CosPA, v0CosPA, float);
DECLARE_SOA_COLUMN(V0MLambda, v0MLambda, float);
DECLARE_SOA_COLUMN(V0MAntiLambda, v0MAntiLambda, float);
DECLARE_SOA_COLUMN(V0MK0Short, v0MK0Short, float);
DECLARE_SOA_COLUMN(V0MGamma, v0MGamma, float);
DECLARE_SOA_COLUMN(FlagMc, flagMc, int8_t);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

DECLARE_SOA_TABLE(HfCandCascFull, "AOD", "HFCANDCASCFull",
                  collision::BCId,
                  collision::NumContrib,
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
                  hf_cand::PxProng0,
                  hf_cand::PyProng0,
                  hf_cand::PzProng0,
                  hf_cand::PxProng1,
                  hf_cand::PyProng1,
                  hf_cand::PzProng1,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0,
                  hf_cand::ErrorImpactParameter1,
                  hf_cand_casc::V0X,
                  hf_cand_casc::V0Y,
                  hf_cand_casc::V0Z,
                  full::V0Radius,
                  full::V0CosPA,
                  full::V0MLambda,
                  full::V0MAntiLambda,
                  full::V0MK0Short,
                  full::V0MGamma,
                  v0data::DCAV0Daughters,
                  v0data::PxPos,
                  v0data::PyPos,
                  v0data::PzPos,
                  full::PtV0Pos,
                  v0data::DCAPosToPV,
                  v0data::PxNeg,
                  v0data::PyNeg,
                  v0data::PzNeg,
                  full::PtV0Neg,
                  v0data::DCANegToPV,
                  full::NSigmaTPCPr0,
                  full::NSigmaTOFPr0,
                  full::M,
                  full::Pt,
                  full::P,
                  full::CPA,
                  full::CPAXY,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::E,
                  full::FlagMc);

DECLARE_SOA_TABLE(HfCandCascFullEvents, "AOD", "HFCANDCASCFullE",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ);

DECLARE_SOA_TABLE(HfCandCascFullParticles, "AOD", "HFCANDCASCFullP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::FlagMc);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorLcToK0sP {
  Produces<o2::aod::HfCandCascFull> rowCandidateFull;
  Produces<o2::aod::HfCandCascFullEvents> rowCandidateFullEvents;
  Produces<o2::aod::HfCandCascFullParticles> rowCandidateFullParticles;

  Configurable<double> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of candidates to store in the tree"};

  void init(InitContext const&)
  {
  }

  template <typename T, typename U>
  void fillCandidate(const T& candidate, const U& bach, int8_t flagMc)
  {
    rowCandidateFull(
      bach.collision().bcId(),
      bach.collision().numContrib(),
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
      candidate.pxProng0(),
      candidate.pyProng0(),
      candidate.pzProng0(),
      candidate.pxProng1(),
      candidate.pyProng1(),
      candidate.pzProng1(),
      candidate.impactParameter0(),
      candidate.impactParameter1(),
      candidate.errorImpactParameter0(),
      candidate.errorImpactParameter1(),
      candidate.v0x(),
      candidate.v0y(),
      candidate.v0z(),
      candidate.v0radius(),
      candidate.v0cosPA(),
      candidate.mLambda(),
      candidate.mAntiLambda(),
      candidate.mK0Short(),
      candidate.mGamma(),
      candidate.dcaV0daughters(),
      candidate.pxpos(),
      candidate.pypos(),
      candidate.pzpos(),
      candidate.ptV0Pos(),
      candidate.dcapostopv(),
      candidate.pxneg(),
      candidate.pyneg(),
      candidate.pzneg(),
      candidate.ptV0Neg(),
      candidate.dcanegtopv(),
      bach.tpcNSigmaPr(),
      bach.tofNSigmaPr(),
      o2::aod::hf_cand_casc::invMassLcToK0sP(candidate),
      candidate.pt(),
      candidate.p(),
      candidate.cpa(),
      candidate.cpaXY(),
      o2::aod::hf_cand_3prong::ctLc(candidate),
      candidate.eta(),
      candidate.phi(),
      o2::aod::hf_cand_3prong::yLc(candidate),
      o2::aod::hf_cand_3prong::eLc(candidate),
      flagMc);
  }

  template <typename T>
  void fillEvent(const T& collision)
  {
    rowCandidateFullEvents(
      collision.bcId(),
      collision.numContrib(),
      collision.posX(),
      collision.posY(),
      collision.posZ());
  }

  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const& mccollisions,
                 soa::Join<aod::HfCandCascade, aod::HfCandCascadeMcRec, aod::HfSelLcToK0sP> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandCascadeMcGen> const& particles,
                 aod::BigTracksPID const& tracks)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (auto const& collision : collisions) {
      fillEvent(collision);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (auto const& candidate : candidates) {
      auto bach = candidate.prong0_as<aod::BigTracksPID>(); // bachelor
      double pseudoRndm = bach.pt() * 1000. - (int16_t)(bach.pt() * 1000);
      if (candidate.isSelLcToK0sP() >= 1 && pseudoRndm < downSampleBkgFactor) {
        fillCandidate(candidate, bach, candidate.flagMcMatchRec());
      }
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (auto const& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == 1) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode())),
          particle.flagMcMatchGen());
      }
    }
  }
  PROCESS_SWITCH(HfTreeCreatorLcToK0sP, processMc, "Process MC tree writer", true);

  void processData(aod::Collisions const& collisions,
                   soa::Join<aod::HfCandCascade, aod::HfSelLcToK0sP> const& candidates,
                   aod::BigTracksPID const& tracks)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (auto const& collision : collisions) {
      fillEvent(collision);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (auto const& candidate : candidates) {
      auto bach = candidate.prong0_as<aod::BigTracksPID>(); // bachelor
      double pseudoRndm = bach.pt() * 1000. - (int16_t)(bach.pt() * 1000);
      if (candidate.isSelLcToK0sP() >= 1 && pseudoRndm < downSampleBkgFactor) {
        fillCandidate(candidate, bach, 0);
      }
    }
  }
  PROCESS_SWITCH(HfTreeCreatorLcToK0sP, processData, "Process data tree writer", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTreeCreatorLcToK0sP>(cfgc),
  };
}
