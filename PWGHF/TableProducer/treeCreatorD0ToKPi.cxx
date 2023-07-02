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
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>, IJCLab

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

enum SelectionType {
  SelD0 = BIT(0),
  SelD0bar = BIT(1),
  MatchedRec = BIT(2),
  ReflectedD0 = BIT(3),
  ReflectedD0bar = BIT(4),
  Prompt = BIT(5),
  NonPrompt = BIT(6)
};

DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(CandidateType, candidateType, int8_t);
DECLARE_SOA_COLUMN(MassD0, massD0, float);
DECLARE_SOA_COLUMN(MassD0bar, massD0bar, float);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(NSigTpcPi0, nSigTpcPi0, float);
DECLARE_SOA_COLUMN(NSigTpcKa0, nSigTpcKa0, float);
DECLARE_SOA_COLUMN(NSigTofPi0, nSigTofPi0, float);
DECLARE_SOA_COLUMN(NSigTofKa0, nSigTofKa0, float);
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);
DECLARE_SOA_COLUMN(NSigTpcKa1, nSigTpcKa1, float);
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);
DECLARE_SOA_COLUMN(NSigTofKa1, nSigTofKa1, float);
DECLARE_SOA_COLUMN(NSigCombTpcTofPi0, nSigCombTpcTofPi0, float);
DECLARE_SOA_COLUMN(NSigCombTpcTofPi1, nSigCombTpcTofPi1, float);
DECLARE_SOA_COLUMN(NSigCombTpcTofKa0, nSigCombTpcTofKa0, float);
DECLARE_SOA_COLUMN(NSigCombTpcTofKa1, nSigCombTpcTofKa1, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(ImpactParameterXY, impactParameterXY, float);
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);
DECLARE_SOA_COLUMN(Cpa, cpa, float);
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);
DECLARE_SOA_COLUMN(CosThetaStar, cosThetaStar, float);
DECLARE_SOA_COLUMN(CosThetaStarD0, cosThetaStarD0, float);
DECLARE_SOA_COLUMN(CosThetaStarD0bar, cosThetaStarD0bar, float);
DECLARE_SOA_COLUMN(PtB, ptB, float);
DECLARE_SOA_COLUMN(FlagMc, flagMc, int8_t);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t); // is prompt or non-prompt, reco level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t); // is prompt or non-prompt, Gen level
} // namespace full

DECLARE_SOA_TABLE(HfCand2ProngLite, "AOD", "HFCAND2PLite",
                  full::CandidateType,
                  full::MassD0,
                  full::MassD0bar,
                  full::Pt,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::ImpactParameterXY,
                  full::ImpactParameterProduct,
                  full::Cpa,
                  full::CpaXY,
                  full::CosThetaStarD0,
                  full::CosThetaStarD0bar,
                  full::MaxNormalisedDeltaIP,
                  full::PtProng0,
                  full::PtProng1,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  full::NSigTpcPi0,
                  full::NSigTpcKa0,
                  full::NSigTofPi0,
                  full::NSigTofKa0,
                  full::NSigTpcPi1,
                  full::NSigTpcKa1,
                  full::NSigTofPi1,
                  full::NSigTofKa1,
                  full::NSigCombTpcTofPi0,
                  full::NSigCombTpcTofPi1,
                  full::NSigCombTpcTofKa0,
                  full::NSigCombTpcTofKa1,
                  full::FlagMc,
                  full::OriginMcRec);

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
                  full::NSigTpcPi0,
                  full::NSigTpcKa0,
                  full::NSigTofPi0,
                  full::NSigTofKa0,
                  full::NSigTpcPi1,
                  full::NSigTpcKa1,
                  full::NSigTofPi1,
                  full::NSigTofKa1,
                  full::CandidateSelFlag,
                  full::M,
                  full::ImpactParameterProduct,
                  full::CosThetaStar,
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
                  full::FlagMc,
                  full::OriginMcGen);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorD0ToKPi {

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<int> selectionFlagHf{"selectionFlagHf", 1, "Selection Flag for HF flagged candidates"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};

  Produces<o2::aod::HfCand2ProngLite> rowCandidateLite;
  Produces<o2::aod::HfCand2ProngFull> rowCandidateFull;
  Produces<o2::aod::HfCand2ProngFullEvents> rowCandidateFullEvents;
  Produces<o2::aod::HfCand2ProngFullParticles> rowCandidateFullParticles;

  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedCandidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<soa::Join<aod::HfCand2Prong, aod::HfCand2ProngMcRec, aod::HfSelD0>> recoFlag2Prongs = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;

  void init(InitContext const&)
  {
  }

  double combineNsigmaTPCTOF(double nsigmaTPC, double nsigmaTOF)
  {
    if (nsigmaTPC > -998. && nsigmaTOF > -998.) {
      return TMath::Sqrt((nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) / 2);
    } else if (nsigmaTPC > -998. && nsigmaTOF < -998.) {
      return TMath::Abs(nsigmaTPC);
    } else if (nsigmaTPC < -998. && nsigmaTOF > -998.) {
      return TMath::Abs(nsigmaTOF);
    } else {
      return -999.;
    }
  }

  template <typename T>
  void fillEvent(const T& collision, int isEventReject, int runNumber)
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

  template <typename T, typename U>
  auto fillLiteTable(const T& candidate, const U& prong0, const U& prong1, int candFlag, double invMassD0, double invMassD0bar, double ctsD0, double ctsD0bar,
                     double nsigCombPi0, double nsigCombPi1, double nsigCombKa0, double nsigCombKa1, int8_t flagMc, int8_t origin)
  {
    rowCandidateLite(
      candFlag,
      invMassD0,
      invMassD0bar,
      candidate.pt(),
      candidate.decayLength(),
      candidate.decayLengthXY(),
      candidate.decayLengthNormalised(),
      candidate.decayLengthXYNormalised(),
      candidate.impactParameterXY(),
      candidate.impactParameterProduct(),
      candidate.cpa(),
      candidate.cpaXY(),
      ctsD0,
      ctsD0bar,
      candidate.maxNormalisedDeltaIP(),
      candidate.ptProng0(),
      candidate.ptProng1(),
      candidate.impactParameter0(),
      candidate.impactParameter1(),
      prong0.tpcNSigmaPi(),
      prong0.tpcNSigmaKa(),
      prong0.tofNSigmaPi(),
      prong0.tofNSigmaKa(),
      prong1.tpcNSigmaPi(),
      prong1.tpcNSigmaKa(),
      prong1.tofNSigmaPi(),
      prong1.tofNSigmaKa(),
      nsigCombPi0,
      nsigCombPi1,
      nsigCombKa0,
      nsigCombKa1,
      flagMc,
      origin);
  }

  template <typename T, typename U>
  auto fillTable(const T& candidate, const U& prong0, const U& prong1, int candFlag, int selection, double invMass, double cosThetaStar,
                 double ct, double y, double e, int8_t flagMc, int8_t origin)
  {
    if (selection >= 1) {
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
        1 << candFlag,
        invMass,
        candidate.impactParameterProduct(),
        cosThetaStar,
        candidate.pt(),
        candidate.p(),
        candidate.cpa(),
        candidate.cpaXY(),
        ct,
        candidate.eta(),
        candidate.phi(),
        y,
        e,
        flagMc,
        origin);
    }
  }

  void processData(aod::Collisions const& collisions,
                   soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates,
                   aod::BigTracksPID const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (auto const& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(selectedCandidates.size());
      for (auto const& candidate : selectedCandidates) {
        if (!(candidate.hfflag() & 1 << DecayType::D0ToPiK)) {
          continue;
        }
        if (yCandMax >= 0. && std::abs(yD0(candidate)) > yCandMax) {
          continue;
        }
        int candType = 0;
        auto prong0 = candidate.prong0_as<aod::BigTracksPID>();
        auto prong1 = candidate.prong1_as<aod::BigTracksPID>();
        double massD0Cand = invMassD0ToPiK(candidate);
        double massD0barCand = invMassD0barToKPi(candidate);
        double ctsD0 = cosThetaStarD0(candidate);
        double ctsD0bar = cosThetaStarD0bar(candidate);
        double nSigCombPi0 = combineNsigmaTPCTOF(prong0.tpcNSigmaPi(), prong0.tofNSigmaPi());
        double nSigCombPi1 = combineNsigmaTPCTOF(prong1.tpcNSigmaPi(), prong1.tofNSigmaPi());
        double nSigCombKa0 = combineNsigmaTPCTOF(prong0.tpcNSigmaKa(), prong0.tofNSigmaKa());
        double nSigCombKa1 = combineNsigmaTPCTOF(prong1.tpcNSigmaKa(), prong1.tofNSigmaKa());
        if (candidate.isSelD0()) {
          candType |= o2::aod::full::SelD0;
        }
        if (candidate.isSelD0bar()) {
          candType |= o2::aod::full::SelD0bar;
        }
        fillLiteTable(candidate, prong0, prong1, candType, massD0Cand, massD0barCand, ctsD0, ctsD0bar, nSigCombPi0, nSigCombPi1, nSigCombKa0, nSigCombKa1, 0, 0);
      }
    } else {
      rowCandidateFull.reserve(candidates.size());
      for (auto const& candidate : candidates) {
        auto prong0 = candidate.prong0_as<aod::BigTracksPID>();
        auto prong1 = candidate.prong1_as<aod::BigTracksPID>();
        double yD = yD0(candidate);
        double eD = eD0(candidate);
        double ctD = ctD0(candidate);
        fillTable(candidate, prong0, prong1, 0, candidate.isSelD0(), invMassD0ToPiK(candidate), cosThetaStarD0(candidate), ctD, yD, eD, 0, 0);
        fillTable(candidate, prong0, prong1, 1, candidate.isSelD0bar(), invMassD0barToKPi(candidate), cosThetaStarD0bar(candidate), ctD, yD, eD, 0, 0);
      }
    }
  }

  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processData, "Process data", true);

  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const&,
                 soa::Join<aod::HfCand2Prong, aod::HfCand2ProngMcRec, aod::HfSelD0> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& particles,
                 aod::BigTracksPID const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (auto const& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    if (fillCandidateLiteTable) {
      rowCandidateFull.reserve(recoFlag2Prongs.size());
      for (auto const& candidate : recoFlag2Prongs) {
        if (!(candidate.hfflag() & 1 << DecayType::D0ToPiK)) {
          continue;
        }
        if (yCandMax >= 0. && std::abs(yD0(candidate)) > yCandMax) {
          continue;
        }
        int candType = 0;
        auto prong0 = candidate.prong0_as<aod::BigTracksPID>();
        auto prong1 = candidate.prong1_as<aod::BigTracksPID>();
        double massD0Cand = invMassD0ToPiK(candidate);
        double massD0barCand = invMassD0barToKPi(candidate);
        double ctsD0 = cosThetaStarD0(candidate);
        double ctsD0bar = cosThetaStarD0bar(candidate);
        double nSigCombPi0 = combineNsigmaTPCTOF(prong0.tpcNSigmaPi(), prong0.tofNSigmaPi());
        double nSigCombPi1 = combineNsigmaTPCTOF(prong1.tpcNSigmaPi(), prong1.tofNSigmaPi());
        double nSigCombKa0 = combineNsigmaTPCTOF(prong0.tpcNSigmaKa(), prong0.tofNSigmaKa());
        double nSigCombKa1 = combineNsigmaTPCTOF(prong1.tpcNSigmaKa(), prong1.tofNSigmaKa());
        if (candidate.isSelD0() >= selectionFlagD0) {
          candType |= o2::aod::full::SelD0;
          if (candidate.flagMcMatchRec() == (1 << DecayType::D0ToPiK) || candidate.flagMcMatchRec() == -(1 << DecayType::D0ToPiK)) {
            candType |= o2::aod::full::MatchedRec;
            if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
              candType |= o2::aod::full::Prompt;
            } else {
              candType |= o2::aod::full::NonPrompt;
            }
            if (candidate.flagMcMatchRec() == -(1 << DecayType::D0ToPiK)) {
              candType |= o2::aod::full::ReflectedD0;
            }
          }
        }
        if (candidate.isSelD0bar() >= selectionFlagD0bar) {
          candType |= o2::aod::full::SelD0bar;
          if (candidate.flagMcMatchRec() == (1 << DecayType::D0ToPiK) || candidate.flagMcMatchRec() == -(1 << DecayType::D0ToPiK)) {
            candType |= o2::aod::full::MatchedRec;
            if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
              candType |= o2::aod::full::Prompt;
            } else {
              candType |= o2::aod::full::NonPrompt;
            }
            if (candidate.flagMcMatchRec() == (1 << DecayType::D0ToPiK)) {
              candType |= o2::aod::full::ReflectedD0bar;
            }
          }
        }
        fillLiteTable(candidate, prong0, prong1, candType, massD0Cand, massD0barCand, ctsD0, ctsD0bar, nSigCombPi0, nSigCombPi1, nSigCombKa0, nSigCombKa1, candidate.flagMcMatchRec(), candidate.originMcRec());
      }
    } else {
      rowCandidateFull.reserve(candidates.size());
      for (auto const& candidate : candidates) {
        auto prong0 = candidate.prong0_as<aod::BigTracksPID>();
        auto prong1 = candidate.prong0_as<aod::BigTracksPID>();
        double yD = yD0(candidate);
        double eD = eD0(candidate);
        double ctD = ctD0(candidate);
        fillTable(candidate, prong0, prong1, 0, candidate.isSelD0(), invMassD0ToPiK(candidate), cosThetaStarD0(candidate), ctD, yD, eD, candidate.flagMcMatchRec(), candidate.originMcRec());
        fillTable(candidate, prong0, prong1, 1, candidate.isSelD0bar(), invMassD0barToKPi(candidate), cosThetaStarD0bar(candidate), ctD, yD, eD, candidate.flagMcMatchRec(), candidate.originMcRec());
      }
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (auto const& particle : particles) {
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

  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorD0ToKPi>(cfgc));
  return workflow;
}
