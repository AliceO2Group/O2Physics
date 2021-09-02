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

/// \file HFTreeCreator3Prong.cxx
/// \brief Writer of the 3 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
///
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN
/// \note based on O2Physics/Tasks/PWGHF/HFTreeCreatorLcToPKPi.cxx

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_x;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
DECLARE_SOA_COLUMN(PtProng2, ptProng2, float);
DECLARE_SOA_COLUMN(PProng2, pProng2, float);
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(CPA, cpa, float);
DECLARE_SOA_COLUMN(CPAXY, cpaXY, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(MCflag, mcflag, int8_t);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

DECLARE_SOA_TABLE(HfCandXFull, "AOD", "HFCANDXFull",
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::PtProng0,
                  full::PProng0,
                  full::PtProng1,
                  full::PProng1,
                  full::PtProng2,
                  full::PProng2,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand::ImpactParameter2,
                  full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::P,
                  full::CPA,
                  full::CPAXY,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::MCflag);

DECLARE_SOA_TABLE(HfCandXFullEvents, "AOD", "HFCANDXFullE",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandXFullParticles, "AOD", "HFCANDXFullP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::MCflag);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorXTojpsipipi {
  Produces<o2::aod::HfCandXFull> rowCandidateFull;
  Produces<o2::aod::HfCandXFullEvents> rowCandidateFullEvents;
  Produces<o2::aod::HfCandXFullParticles> rowCandidateFullParticles;

  void init(InitContext const&)
  {
  }

  void process(aod::Collisions const& collisions,
               aod::McCollisions const& mccollisions,
               soa::Join<aod::HfCandX, aod::HfCandXMCRec, aod::HFSelXToJpsiPiPiCandidate> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandXMCGen> const& particles,
               aod::BigTracksPID const& tracks)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (auto& collision : collisions) {
      rowCandidateFullEvents(
        collision.bcId(),
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        0,
        1);
    }

    // Filling candidate properties
    int indexCand = 0;
    rowCandidateFull.reserve(candidates.size());
    for (auto& candidate : candidates) {
      indexCand++;
      auto fillTable = [&](int CandFlag,
                           int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY) {
        if (FunctionSelection >= 1) {
          rowCandidateFull(
            candidate.posX(),
            candidate.posY(),
            candidate.posZ(),
            candidate.ptProng0(),
            RecoDecay::P(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
            candidate.ptProng1(),
            RecoDecay::P(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
            candidate.ptProng2(),
            RecoDecay::P(candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()),
            candidate.impactParameter0(),
            candidate.impactParameter1(),
            candidate.impactParameter2(),
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
            candidate.flagMCMatchRec());
        }
      };

      fillTable(0, candidate.isSelXToJpsiPiPi(), InvMassXToJpsiPiPi(candidate), CtX(candidate), YX(candidate));
    }

    // Filling particle properties
    float massX = 3.872;
    rowCandidateFullParticles.reserve(particles.size());
    for (auto& particle : particles) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << DecayType::XToJpsiPiPi) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, massX),
          particle.flagMCMatchGen());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorXTojpsipipi>(cfgc));
  return workflow;
}
