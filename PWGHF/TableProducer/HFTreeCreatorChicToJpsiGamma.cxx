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
/// \author Alessandro De Falco <alessandro.de.falco@ca.infn.it>, Università/INFN Cagliari
/// \note based on O2Physics/Tasks/PWGHF/HFTreeCreatorXToJPsiPiPi.cxx
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand_chic;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
DECLARE_SOA_COLUMN(Alpha, alpha, float);
DECLARE_SOA_COLUMN(Qt, qt, float);
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(CPA, cpa, float);
DECLARE_SOA_COLUMN(CPAXY, cpaXY, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(MCflag, mcflag, int8_t);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

DECLARE_SOA_TABLE(HfCandChicFull, "AOD", "HFCANDChicFull",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::ImpactParameterNormalised0,
                  full::RSecondaryVertex,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::PtProng0,
                  full::PProng0,
                  full::PtProng1,
                  full::PProng1,
                  full::Alpha,
                  full::Qt,
                  hf_cand::Chi2PCA,
                  hf_cand::ImpactParameter0,
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
                  hf_cand_chic::JpsiToMuMuMass,
                  full::MCflag);

DECLARE_SOA_TABLE(HfCandChicFullEvents, "AOD", "HFCANDChicFullE",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandChicFullParticles, "AOD", "HFCANDChicFullP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  hf_cand_chic::JpsiToMuMuMass,
                  full::MCflag);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorChicToJpsiGamma {
  Produces<o2::aod::HfCandChicFull> rowCandidateFull;
  Produces<o2::aod::HfCandChicFullEvents> rowCandidateFullEvents;
  Produces<o2::aod::HfCandChicFullParticles> rowCandidateFullParticles;

  void init(InitContext const&)
  {
  }

  void process(aod::Collisions const& collisions,
               aod::McCollisions const& mccollisions,
               soa::Join<aod::HfCandChic, aod::HfCandChicMCRec, aod::HFSelChicToJpsiGammaCandidate> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandChicMCGen> const& particles,
               aod::BigTracksPID const& tracks,
               aod::HfCandProng2 const& jpsiCands)
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
      array<float, 3> pvecChic = {candidate.px(), candidate.py(), candidate.pz()};
      array<float, 3> pvecJpsi = {candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()};
      array<float, 3> pvecGamma = {candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()};
      auto pchic = RecoDecay::p(pvecChic);
      auto pjpsi = RecoDecay::p(pvecJpsi);
      auto pl1 = std::abs(RecoDecay::dotProd(pvecChic, pvecJpsi)) / pchic;
      auto pl2 = std::abs(RecoDecay::dotProd(pvecChic, pvecGamma)) / pchic;
      auto alpha = (pl1 - pl2) / (pl1 + pl2);
      auto qt = std::sqrt(pjpsi * pjpsi - pl1 * pl1);

      indexCand++;
      auto fillTable = [&](int CandFlag,
                           int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY) {
        if (FunctionSelection >= 1) {
          rowCandidateFull(
            candidate.index0().index0_as<aod::BigTracksPID>().collision().bcId(),
            candidate.index0().index0_as<aod::BigTracksPID>().collision().numContrib(),
            candidate.posX(),
            candidate.posY(),
            candidate.posZ(),
            candidate.impactParameterNormalised0(),
            candidate.rSecondaryVertex(),
            candidate.decayLength(),
            candidate.decayLengthXY(),
            candidate.ptProng0(),
            RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
            candidate.ptProng1(),
            RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
            alpha,
            qt,
            candidate.chi2PCA(),
            candidate.impactParameter0(),
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
            candidate.jpsiToMuMuMass(),
            candidate.flagMCMatchRec());
        }
      };
      fillTable(0, candidate.isSelChicToJpsiToMuMuGamma(), InvMassChicToJpsiGamma(candidate), CtChic(candidate), YChic(candidate));
      //      fillTable(1, candidate.isSelChicToJpsiToEEGamma(), InvMassChicToJpsiGamma(candidate), CtChic(candidate), YChic(candidate));
    }

    // Filling particle properties
    float massChic = RecoDecay::getMassPDG(pdg::Code::kChic1);
    rowCandidateFullParticles.reserve(particles.size());
    for (auto& particle : particles) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << DecayType::ChicToJpsiToEEGamma || std::abs(particle.flagMCMatchGen()) == 1 << DecayType::ChicToJpsiToMuMuGamma) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, massChic),
          0., // put here the jpsi mass
          particle.flagMCMatchGen());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorChicToJpsiGamma>(cfgc));
  return workflow;
}
