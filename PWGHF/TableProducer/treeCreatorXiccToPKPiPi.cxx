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

/// \file treeCreatorXiccToPKPiPi.cxx
/// \brief Writer of the 2 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
/// \note Extended from treeCreatorD0ToKPi.cxx, treeCreatorLcToPKPi.cxx, treeCreatorXToJpsiPiPi.cxx
///
/// \author Jinjoo Seo <jin.joo.seo@cern.ch>, Inha University

#include "Common/Core/trackUtilities.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_xicc;

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
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(NSigmaTOFBachPi, nSigmaTOFBachPi, float);
DECLARE_SOA_COLUMN(Chi2PCA, chi2PCA, float);
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
// Xic selection variable
DECLARE_SOA_COLUMN(XicM, xicM, float);
DECLARE_SOA_COLUMN(XicCt, xicCt, float);
DECLARE_SOA_COLUMN(XicY, xicY, float);
DECLARE_SOA_COLUMN(XicE, xicE, float);
DECLARE_SOA_COLUMN(XicEta, xicEta, float);
DECLARE_SOA_COLUMN(XicCPA, xicCPA, float);
DECLARE_SOA_COLUMN(XicCPAXY, xicCPAXY, float);
DECLARE_SOA_COLUMN(XicChi2PCA, xicChi2PCA, float);
DECLARE_SOA_COLUMN(XicDecayLength, xicDecayLength, float);
DECLARE_SOA_COLUMN(XicDecayLengthXY, xicDecayLengthXY, float);
DECLARE_SOA_COLUMN(XicDecayLengthNormalised, xicDecayLengthNormalised, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk1Pr, nSigmaTOFTrk1Pr, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk1Pi, nSigmaTOFTrk1Pi, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk2Ka, nSigmaTOFTrk2Ka, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk3Pr, nSigmaTOFTrk3Pr, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk3Pi, nSigmaTOFTrk3Pi, float);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

// put the arguments into the table
DECLARE_SOA_TABLE(HfCandXiccFull, "AOD", "HFCANDXiccFull",
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
                  hf_cand::Chi2PCA,
                  full::NSigmaTOFBachPi,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0,
                  hf_cand::ErrorImpactParameter1,
                  full::ImpactParameterProduct,
                  full::XicM,
                  full::XicCt,
                  full::XicY,
                  full::XicE,
                  full::XicEta,
                  full::XicCPA,
                  full::XicCPAXY,
                  full::XicChi2PCA,
                  full::XicDecayLength,
                  full::XicDecayLengthXY,
                  full::XicDecayLengthNormalised,
                  full::NSigmaTOFTrk1Pr,
                  full::NSigmaTOFTrk1Pi,
                  full::NSigmaTOFTrk2Ka,
                  full::NSigmaTOFTrk3Pr,
                  full::NSigmaTOFTrk3Pi,
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
                  full::MCflag,
                  full::OriginMcRec);

DECLARE_SOA_TABLE(HfCandXiccFullEvents, "AOD", "HFCANDXiccFullE",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandXiccFullParticles, "AOD", "HFCANDXiccFullP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::MCflag,
                  full::OriginMcGen);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorXiccToPKPiPi {
  Produces<o2::aod::HfCandXiccFull> rowCandidateFull;
  Produces<o2::aod::HfCandXiccFullEvents> rowCandidateFullEvents;
  Produces<o2::aod::HfCandXiccFullParticles> rowCandidateFullParticles;

  void init(InitContext const&)
  {
  }

  void process(aod::Collisions const& collisions,
               aod::McCollisions const& mccollisions,
               soa::Join<aod::HfCandXicc, aod::HfCandXiccMcRec, aod::HfSelXiccToPKPiPi> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandXiccMcGen> const& particles,
               aod::BigTracksPID const& tracks,
               aod::HfCand3Prong const&)
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
                           int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY) {
        if (FunctionSelection >= 1) {
          auto xicCand = candidate.prong0();

          rowCandidateFull(
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
            candidate.chi2PCA(),
            candidate.prong1_as<aod::BigTracksPID>().tofNSigmaPi(),
            candidate.impactParameter0(),
            candidate.impactParameter1(),
            candidate.errorImpactParameter0(),
            candidate.errorImpactParameter1(),
            candidate.impactParameterProduct(),
            o2::aod::hf_cand_3prong::invMassXicToPKPi(xicCand),
            o2::aod::hf_cand_3prong::ctXic(xicCand),
            o2::aod::hf_cand_3prong::yXic(xicCand),
            o2::aod::hf_cand_3prong::eXic(xicCand),
            xicCand.eta(),
            xicCand.cpa(),
            xicCand.cpaXY(),
            xicCand.chi2PCA(),
            xicCand.decayLength(),
            xicCand.decayLengthXY(),
            xicCand.decayLengthXYNormalised(),
            xicCand.prong0_as<aod::BigTracksPID>().tofNSigmaPr(),
            xicCand.prong0_as<aod::BigTracksPID>().tofNSigmaPi(),
            xicCand.prong1_as<aod::BigTracksPID>().tofNSigmaKa(),
            xicCand.prong2_as<aod::BigTracksPID>().tofNSigmaPr(),
            xicCand.prong2_as<aod::BigTracksPID>().tofNSigmaPi(),
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
            candidate.flagMcMatchRec(),
            candidate.originMcRec());
        }
      };

      fillTable(0, candidate.isSelXiccToPKPiPi(), invMassXiccToXicPi(candidate), ctXicc(candidate), yXicc(candidate));
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (auto& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << DecayType::XiccToXicPi) {
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorXiccToPKPiPi>(cfgc));
  return workflow;
}
