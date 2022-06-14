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

/// \file HFTreeCreatorBplustoD0Pi.cxx
/// \brief Writer of the 2 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
/// \note Extended from HFTreeCreatorD0ToKPi
///
/// \author Antonio Palasciano <antonio.palasciano@ba.infn.it>, Università & INFN, Bari

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
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_bplus;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
//DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(CPA, cpa, float);
DECLARE_SOA_COLUMN(CPAXY, cpaXY, float);
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);
DECLARE_SOA_COLUMN(ImpactParameter0, impactParameter0, float);
DECLARE_SOA_COLUMN(ImpactParameter1, impactParameter1, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(Chi2PCA, chi2PCA, float);
DECLARE_SOA_COLUMN(NSigmaTOFBachPi, nSigmaTOFBachPi, float);
DECLARE_SOA_COLUMN(NSigmaTOFBachKa, nSigmaTOFBachKa, float);
DECLARE_SOA_COLUMN(NSigmaTPCBachPi, nSigmaTPCBachPi, float);
DECLARE_SOA_COLUMN(NSigmaTPCBachKa, nSigmaTPCBachKa, float);
DECLARE_SOA_COLUMN(MCflag, mcflag, int8_t);
// D0 (Prong0) selection variable
DECLARE_SOA_COLUMN(D0M, d0M, float);
DECLARE_SOA_COLUMN(D0Ct, d0Ct, float);
DECLARE_SOA_COLUMN(D0PtProng0, d0ptProng0, float);
DECLARE_SOA_COLUMN(D0PtProng1, d0ptProng1, float);
DECLARE_SOA_COLUMN(D0Y, d0Y, float);
DECLARE_SOA_COLUMN(D0Eta, d0Eta, float);
DECLARE_SOA_COLUMN(D0CPA, d0CPA, float);
DECLARE_SOA_COLUMN(D0CPAXY, d0CPAXY, float);
DECLARE_SOA_COLUMN(D0Chi2PCA, d0Chi2PCA, float);
DECLARE_SOA_COLUMN(D0DecayLength, d0DecayLength, float);
DECLARE_SOA_COLUMN(D0DecayLengthXY, d0DecayLengthXY, float);
DECLARE_SOA_COLUMN(D0DecayLengthNormalised, d0DecayLengthNormalised, float);
DECLARE_SOA_COLUMN(D0DecayLengthXYNormalised, d0decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(D0ImpactParameterProduct, d0impactParameterProduct, float);
DECLARE_SOA_COLUMN(D0ImpactParameter0, d0impactParameter0, float);
DECLARE_SOA_COLUMN(D0ImpactParameter1, d0impactParameter1, float);
DECLARE_SOA_COLUMN(D0ImpactParameterNormalised0, d0impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(D0ImpactParameterNormalised1, d0impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk0Ka, nSigmaTOFTrk0Ka, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk0Pi, nSigmaTOFTrk0Pi, float);
DECLARE_SOA_COLUMN(NSigmaTPCTrk0Ka, nSigmaTPCTrk0Ka, float);
DECLARE_SOA_COLUMN(NSigmaTPCTrk0Pi, nSigmaTPCTrk0Pi, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk1Ka, nSigmaTOFTrk1Ka, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk1Pi, nSigmaTOFTrk1Pi, float);
DECLARE_SOA_COLUMN(NSigmaTPCTrk1Ka, nSigmaTPCTrk1Ka, float);
DECLARE_SOA_COLUMN(NSigmaTPCTrk1Pi, nSigmaTPCTrk1Pi, float);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

// put the arguments into the table
DECLARE_SOA_TABLE(HfCandBplusFull, "AOD", "HFCANDBPFull",
                  full::RSecondaryVertex,
                  full::PtProng0,
                  full::PProng0,
                  full::PtProng1,
                  full::PProng1,
                  //full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::P,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::CPA,
                  full::CPAXY,
                  full::ImpactParameterProduct,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  full::ImpactParameterNormalised0,
                  full::ImpactParameterNormalised1,
                  hf_cand::PxProng0,
                  hf_cand::PyProng0,
                  hf_cand::PzProng0,
                  hf_cand::PxProng1,
                  hf_cand::PyProng1,
                  hf_cand::PzProng1,
                  hf_cand::Chi2PCA,
                  full::NSigmaTOFBachPi,
                  full::NSigmaTOFBachKa,
                  full::NSigmaTPCBachPi,
                  full::NSigmaTPCBachKa,
                  full::MCflag,
                  full::D0M,
                  full::D0PtProng0,
                  full::D0PtProng1,
                  full::D0Y,
                  full::D0Eta,
                  full::D0CPA,
                  full::D0CPAXY,
                  full::D0Chi2PCA,
                  full::D0DecayLength,
                  full::D0DecayLengthXY,
                  full::D0DecayLengthNormalised,
                  full::D0DecayLengthXYNormalised,
                  full::D0ImpactParameterProduct,
                  full::D0ImpactParameter0,
                  full::D0ImpactParameter1,
                  full::D0ImpactParameterNormalised0,
                  full::D0ImpactParameterNormalised1,
                  full::NSigmaTOFTrk0Pi,
                  full::NSigmaTOFTrk0Ka,
                  full::NSigmaTPCTrk0Pi,
                  full::NSigmaTPCTrk0Ka,
                  full::NSigmaTOFTrk1Pi,
                  full::NSigmaTOFTrk1Ka,
                  full::NSigmaTPCTrk1Pi,
                  full::NSigmaTPCTrk1Ka);

DECLARE_SOA_TABLE(HfCandBplusFullEvents, "AOD", "HFCANDBPFullE",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandBplusFullParticles, "AOD", "HFCANDBPFullP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::MCflag);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorBplusToD0Pi {
  Produces<o2::aod::HfCandBplusFull> rowCandidateFull;
  Produces<o2::aod::HfCandBplusFullEvents> rowCandidateFullEvents;
  Produces<o2::aod::HfCandBplusFullParticles> rowCandidateFullParticles;

  void init(InitContext const&)
  {
  }

  Configurable<int> isSignal{"isSignal", 1, "save only MC matched candidates"};

  void process(aod::Collisions const& collisions,
               aod::McCollisions const& mccollisions,
               soa::Join<aod::HfCandBPlus, aod::HfCandBPMCRec, aod::HFSelBPlusToD0PiCandidate> const& candidates,
               soa::Join<aod::McParticles_000, aod::HfCandBPMCGen> const& particles,
               aod::BigTracksPID const& tracks,
               aod::HfCandProng2 const&)
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
    rowCandidateFull.reserve(candidates.size());
    for (auto& candidate : candidates) {
      auto fillTable = [&](int CandFlag,
                           //int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY) {
        auto d0Cand = candidate.index0();
        auto piCand = candidate.index1_as<aod::BigTracksPID>();
        //adding D0 daughters to the table
        auto d0Daughter0 = d0Cand.index0_as<aod::BigTracksPID>();
        auto d0Daughter1 = d0Cand.index1_as<aod::BigTracksPID>();

        auto invMassD0 = 0.;
        if (piCand.sign() > 0) {
          invMassD0 = o2::aod::hf_cand_prong2::InvMassD0bar(d0Cand);
        } else if (piCand.sign() < 0) {
          invMassD0 = o2::aod::hf_cand_prong2::InvMassD0(d0Cand);
        }

        //if (FunctionSelection >= 1) {
        if (std::abs(candidate.flagMCMatchRec()) >= isSignal) {

          rowCandidateFull(
            candidate.rSecondaryVertex(),
            candidate.ptProng0(),
            RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
            candidate.ptProng1(),
            RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
            //1 << CandFlag,
            FunctionInvMass,
            candidate.pt(),
            candidate.p(),
            FunctionCt,
            candidate.eta(),
            candidate.phi(),
            FunctionY,
            candidate.decayLength(),
            candidate.decayLengthXY(),
            candidate.decayLengthNormalised(),
            candidate.decayLengthXYNormalised(),
            candidate.cpa(),
            candidate.cpaXY(),
            candidate.impactParameterProduct(),
            candidate.impactParameter0(),
            candidate.impactParameter1(),
            candidate.impactParameterNormalised0(),
            candidate.impactParameterNormalised1(),
            candidate.pxProng0(),
            candidate.pyProng0(),
            candidate.pzProng0(),
            candidate.pxProng1(),
            candidate.pyProng1(),
            candidate.pzProng1(),
            candidate.chi2PCA(),
            piCand.tofNSigmaPi(),
            piCand.tofNSigmaKa(),
            piCand.tpcNSigmaPi(),
            piCand.tpcNSigmaKa(),
            candidate.flagMCMatchRec(),
            invMassD0,
            d0Cand.ptProng0(),
            d0Cand.ptProng1(),
            o2::aod::hf_cand_prong2::YD0(d0Cand),
            d0Cand.eta(),
            d0Cand.cpa(),
            d0Cand.cpaXY(),
            d0Cand.chi2PCA(),
            d0Cand.decayLength(),
            d0Cand.decayLengthXY(),
            d0Cand.decayLengthNormalised(),
            d0Cand.decayLengthXYNormalised(),
            d0Cand.impactParameterProduct(),
            d0Cand.impactParameter0(),
            d0Cand.impactParameter1(),
            d0Cand.impactParameterNormalised0(),
            d0Cand.impactParameterNormalised1(),
            d0Daughter0.tofNSigmaPi(),
            d0Daughter0.tofNSigmaKa(),
            d0Daughter0.tpcNSigmaPi(),
            d0Daughter0.tpcNSigmaKa(),
            d0Daughter1.tofNSigmaPi(),
            d0Daughter1.tofNSigmaKa(),
            d0Daughter1.tpcNSigmaPi(),
            d0Daughter1.tpcNSigmaKa());
        }
      };

      //fillTable(0, candidate.isSelBPlusToD0Pi(), InvMassBPlus(candidate), CtBPlus(candidate), YBPlus(candidate));
      fillTable(0, InvMassBPlus(candidate), CtBPlus(candidate), YBPlus(candidate));
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (auto& particle : particles) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << DecayType::BPlusToD0Pi) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode())),
          particle.flagMCMatchGen());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorBplusToD0Pi>(cfgc));
  return workflow;
}
