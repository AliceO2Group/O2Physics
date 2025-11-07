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

/// \file treeCreatorXToJpsiPiPi.cxx
/// \brief Writer of the 3 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
/// \note based on treeCreatorLcToPKPi.cxx
///
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
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
DECLARE_SOA_COLUMN(Q, q, float);
DECLARE_SOA_COLUMN(DR1, dr1, float);
DECLARE_SOA_COLUMN(DR2, dr2, float);
DECLARE_SOA_COLUMN(PiBalance, piBalance, float);
DECLARE_SOA_COLUMN(NSigmaTOFPi1, nSigmaTOFPi1, float);
DECLARE_SOA_COLUMN(NSigmaTOFKa1, nSigmaTOFKa1, float);
DECLARE_SOA_COLUMN(NSigmaTOFPr1, nSigmaTOFPr1, float);
DECLARE_SOA_COLUMN(NSigmaTOFPi2, nSigmaTOFPi2, float);
DECLARE_SOA_COLUMN(NSigmaTOFKa2, nSigmaTOFKa2, float);
DECLARE_SOA_COLUMN(NSigmaTOFPr2, nSigmaTOFPr2, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(CPA, cpa, float);
DECLARE_SOA_COLUMN(CPAXY, cpaXY, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(MCflag, mcflag, int8_t);
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

DECLARE_SOA_TABLE(HfCandXFulls, "AOD", "HFCANDXFULL",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::ImpactParameterNormalised0,
                  full::ImpactParameterNormalised1,
                  full::RSecondaryVertex,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::PtProng0,
                  full::PProng0,
                  full::PtProng1,
                  full::PProng1,
                  full::PtProng2,
                  full::PProng2,
                  full::NSigmaTOFPi1,
                  full::NSigmaTOFKa1,
                  full::NSigmaTOFPr1,
                  full::NSigmaTOFPi2,
                  full::NSigmaTOFKa2,
                  full::NSigmaTOFPr2,
                  hf_cand::Chi2PCA,
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
                  full::Q,
                  full::DR1,
                  full::DR2,
                  full::PiBalance,
                  full::MCflag,
                  full::OriginMcRec);

DECLARE_SOA_TABLE(HfCandXFullEvs, "AOD", "HFCANDXFULLEV",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandXFullPs, "AOD", "HFCANDXFULLP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::MCflag,
                  full::OriginMcGen);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorXToJpsiPiPi {
  Produces<o2::aod::HfCandXFulls> rowCandidateFull;
  Produces<o2::aod::HfCandXFullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandXFullPs> rowCandidateFullParticles;

  using TracksWPid = soa::Join<aod::Tracks, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

  void init(InitContext const&)
  {
  }

  void process(aod::Collisions const& collisions,
               aod::McCollisions const&,
               soa::Join<aod::HfCandX, aod::HfCandXMcRec, aod::HfSelXToJpsiPiPi> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandXMcGen> const& particles,
               TracksWPid const&)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
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
    // int indexCand = 0;
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      if (!candidate.isSelXToJpsiToMuMuPiPi()) {
        continue;
      }
      // indexCand++;
      auto fillTable = [&](int CandFlag,
                           int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY,
                           float FunctionQ,
                           float FunctionDR1,
                           float FunctionDR2,
                           float FunctionPiBalance) {
        if (FunctionSelection >= 1) {
          rowCandidateFull(
            candidate.prong1_as<TracksWPid>().collision().bcId(),
            candidate.prong1_as<TracksWPid>().collision().numContrib(),
            candidate.posX(),
            candidate.posY(),
            candidate.posZ(),
            candidate.chi2PCA(),
            candidate.rSecondaryVertex(),
            candidate.decayLength(),
            candidate.decayLengthXY(),
            candidate.impactParameterNormalised0(),
            candidate.ptProng0(),
            RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
            candidate.impactParameterNormalised1(),
            candidate.ptProng1(),
            RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
            candidate.ptProng2(),
            RecoDecay::p(candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()),
            candidate.impactParameter0(),
            candidate.impactParameter1(),
            candidate.impactParameter2(),
            candidate.prong1_as<TracksWPid>().tofNSigmaPi(),
            candidate.prong1_as<TracksWPid>().tofNSigmaKa(),
            candidate.prong1_as<TracksWPid>().tofNSigmaPr(),
            candidate.prong2_as<TracksWPid>().tofNSigmaPi(),
            candidate.prong2_as<TracksWPid>().tofNSigmaKa(),
            candidate.prong2_as<TracksWPid>().tofNSigmaPr(),
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
            FunctionQ,
            FunctionDR1,
            FunctionDR2,
            FunctionPiBalance,
            candidate.flagMcMatchRec(),
            candidate.originMcRec());
        }
      };

      fillTable(0, candidate.isSelXToJpsiToMuMuPiPi(), HfHelper::invMassXToJpsiPiPi(candidate), HfHelper::ctX(candidate), HfHelper::yX(candidate), HfHelper::qX(candidate), HfHelper::dRX(candidate, 1), HfHelper::dRX(candidate, 2), HfHelper::balancePtPionsX(candidate));
    }

    // Filling particle properties
    float massX = 3.872;
    rowCandidateFullParticles.reserve(particles.size());
    for (const auto& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << hf_cand_x::DecayType::XToJpsiToEEPiPi || std::abs(particle.flagMcMatchGen()) == 1 << hf_cand_x::DecayType::XToJpsiToMuMuPiPi) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(particle.pVector(), massX),
          particle.flagMcMatchGen(),
          particle.originMcGen());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorXToJpsiPiPi>(cfgc));
  return workflow;
}
