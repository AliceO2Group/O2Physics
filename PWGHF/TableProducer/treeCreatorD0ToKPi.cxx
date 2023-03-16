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

/// \file treeCreatorD0ToKPi.cxx
/// \brief Writer of the 2 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
///
/// \author Nicolo' Jacazio <nicolo.jacazio@cern.ch>, CERN

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_2prong;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(NSigTPCPi0, nsigTPCPi0, float);
DECLARE_SOA_COLUMN(NSigTPCKa0, nsigTPCKa0, float);
DECLARE_SOA_COLUMN(NSigTOFPi0, nsigTOFPi0, float);
DECLARE_SOA_COLUMN(NSigTOFKa0, nsigTOFKa0, float);
DECLARE_SOA_COLUMN(NSigTPCPi1, nsigTPCPi1, float);
DECLARE_SOA_COLUMN(NSigTPCKa1, nsigTPCKa1, float);
DECLARE_SOA_COLUMN(NSigTOFPi1, nsigTOFPi1, float);
DECLARE_SOA_COLUMN(NSigTOFKa1, nsigTOFKa1, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(CPA, cpa, float);
DECLARE_SOA_COLUMN(CPAXY, cpaXY, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);
DECLARE_SOA_COLUMN(CosThetaStar, cosThetaStar, float);
DECLARE_SOA_COLUMN(Mcflag, mcflag, int8_t);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t); // is prompt or non-prompt, reco level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t); // is prompt or non-prompt, Gen level
} // namespace full

DECLARE_SOA_TABLE(HfCand2ProngFull, "AOD", "HFCAND2PFull",
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
                  full::NSigTPCPi0,
                  full::NSigTPCKa0,
                  full::NSigTOFPi0,
                  full::NSigTOFKa0,
                  full::NSigTPCPi1,
                  full::NSigTPCKa1,
                  full::NSigTOFPi1,
                  full::NSigTOFKa1,
                  full::CandidateSelFlag,
                  full::M,
                  full::ImpactParameterProduct,
                  full::CosThetaStar,
                  full::Pt,
                  full::P,
                  full::CPA,
                  full::CPAXY,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::E,
                  full::Mcflag,
                  full::OriginMcRec);

DECLARE_SOA_TABLE(HfCand2ProngFullEvents, "AOD", "HFCAND2PFullE",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCand2ProngFullParticles, "AOD", "HFCAND2PFullP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::Mcflag,
                  full::OriginMcGen);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorD0ToKPi {
  Produces<o2::aod::HfCand2ProngFull> rowCandidateFull;
  Produces<o2::aod::HfCand2ProngFullEvents> rowCandidateFullEvents;
  Produces<o2::aod::HfCand2ProngFullParticles> rowCandidateFullParticles;

  void init(InitContext const&)
  {
  }

  template <typename T>
  void fillEvent(const T& collision, int8_t isEventReject, int8_t runNumber)
  {
    rowCandidateFullEvents(
        collision.bcId(),
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        isEventReject,
        runNumber);
  }

template <typename T, typename U, typename V>
auto fillTable(const T& candidate, const U& prong0, const V& prong1, int CandFlag, int FunctionSelection, double FunctionInvMass, double FunctionCosThetaStar,
            double FunctionCt, double FunctionY, double FunctionE, int8_t flagMc, int8_t origin) {
        if (FunctionSelection >= 1) {
          rowCandidateFull(
            prong0.collision().bcId(),
            prong0.collision().numContrib(),
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
            prong0.tpcNSigmaPi(),
            prong0.tpcNSigmaKa(),
            prong0.tofNSigmaPi(),
            prong0.tofNSigmaKa(),
            prong1.tpcNSigmaPi(),
            prong1.tpcNSigmaKa(),
            prong1.tofNSigmaPi(),
            prong1.tofNSigmaKa(),
            1 << CandFlag,
            FunctionInvMass,
            candidate.impactParameterProduct(),
            FunctionCosThetaStar,
            candidate.pt(),
            candidate.p(),
            candidate.cpa(),
            candidate.cpaXY(),
            FunctionCt,
            candidate.eta(),
            candidate.phi(),
            FunctionY,
            FunctionE,
            flagMc,
            origin);
        }
  }


  void processData(aod::Collisions const& collisions,
               soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates,
               aod::BigTracksPID const& tracks)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (auto& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (auto& candidate : candidates) {
      auto prong0 = candidate.prong0_as<aod::BigTracksPID>();
      auto prong1 = candidate.prong1_as<aod::BigTracksPID>();
      fillTable(candidate, prong0, prong1, 0, candidate.isSelD0(), invMassD0ToPiK(candidate), cosThetaStarD0(candidate), ctD0(candidate), yD0(candidate), eD0(candidate), 0, 0);
      fillTable(candidate, prong0, prong1, 1, candidate.isSelD0bar(), invMassD0barToKPi(candidate), cosThetaStarD0bar(candidate), ctD0(candidate), yD0(candidate), eD0(candidate), 0, 0);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processData, "Process data", true);

  void processMc(aod::Collisions const& collisions,
               aod::McCollisions const& mccollisions,
               soa::Join<aod::HfCand2Prong, aod::HfCand2ProngMcRec, aod::HfSelD0> const& candidates,
               soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& particles,
               aod::BigTracksPID const& tracks)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (auto& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (auto& candidate : candidates) {
      auto prong0 = candidate.prong0_as<aod::BigTracksPID>();
      auto prong1 = candidate.prong0_as<aod::BigTracksPID>();
      fillTable(candidate, prong0, prong1, 0, candidate.isSelD0(), invMassD0ToPiK(candidate), cosThetaStarD0(candidate), ctD0(candidate), yD0(candidate), eD0(candidate), candidate.flagMcMatchRec(), candidate.originMcRec());
      fillTable(candidate, prong0, prong1, 1, candidate.isSelD0bar(), invMassD0barToKPi(candidate), cosThetaStarD0bar(candidate), ctD0(candidate), yD0(candidate), eD0(candidate), candidate.flagMcMatchRec(), candidate.originMcRec());
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (auto& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << DecayType::D0ToPiK) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode())),
          particle.flagMcMatchGen(),
          particle.originMcGen());
      }
    }
  }

  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMc, "Process Mc", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorD0ToKPi>(cfgc));
  return workflow;
}