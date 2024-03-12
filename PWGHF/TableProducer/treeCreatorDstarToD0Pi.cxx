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

/// \file treeCreatorDstarToD0Pi.cxx
/// \brief Writer of D*+ → D0 ( → π+ K-) π+ candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug, local optimization of analysis on small samples or ML training.
///        In this file are defined and filled the output tables
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
// D0 related variables
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(MD0, mD0, float);
DECLARE_SOA_COLUMN(PtD0, ptD0, float);
DECLARE_SOA_COLUMN(PD0, pD0, float);
DECLARE_SOA_COLUMN(EtaD0, etaD0, float);
DECLARE_SOA_COLUMN(PhiD0, phiD0, float);
DECLARE_SOA_COLUMN(YD0, yD0, float);
DECLARE_SOA_COLUMN(NSigTpcPi0, nSigTpcPi0, float);
DECLARE_SOA_COLUMN(NSigTpcKa0, nSigTpcKa0, float);
DECLARE_SOA_COLUMN(NSigTofPi0, nSigTofPi0, float);
DECLARE_SOA_COLUMN(NSigTofKa0, nSigTofKa0, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi0, nSigTpcTofPi0, float);
DECLARE_SOA_COLUMN(NSigTpcTofKa0, nSigTpcTofKa0, float);
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);
DECLARE_SOA_COLUMN(NSigTpcKa1, nSigTpcKa1, float);
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);
DECLARE_SOA_COLUMN(NSigTofKa1, nSigTofKa1, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi1, nSigTpcTofPi1, float);
DECLARE_SOA_COLUMN(NSigTpcTofKa1, nSigTpcTofKa1, float);
DECLARE_SOA_COLUMN(DecayLengthD0, decayLengthD0, float);
DECLARE_SOA_COLUMN(DecayLengthXYD0, decayLengthXYD0, float);
DECLARE_SOA_COLUMN(DecayLengthNormalisedD0, decayLengthNormalisedD0, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalisedD0, decayLengthXYNormalisedD0, float);
DECLARE_SOA_COLUMN(CpaD0, cpaD0, float);
DECLARE_SOA_COLUMN(CpaXYD0, cpaXYD0, float);
DECLARE_SOA_COLUMN(CosThetaStarD0, cosThetaStarD0, float);
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIPD0, maxNormalisedDeltaIPD0, float);
DECLARE_SOA_COLUMN(CtD0, ctD0, float);
DECLARE_SOA_COLUMN(ImpactParameterProductD0, impactParameterProductD0, float);
// soft pion related variables
DECLARE_SOA_COLUMN(PtSoftPi, ptSoftPi, float);
DECLARE_SOA_COLUMN(PSoftPi, pSoftPi, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalisedSoftPi, impactParameterNormalisedSoftPi, float);
DECLARE_SOA_COLUMN(NSigTpcPiSoftPi, nSigTpcPiSoftPi, float);
DECLARE_SOA_COLUMN(NSigTpcKaSoftPi, nSigTpcKaSoftPi, float);
DECLARE_SOA_COLUMN(NSigTofPiSoftPi, nSigTofPiSoftPi, float);
DECLARE_SOA_COLUMN(NSigTofKaSoftPi, nSigTofKaSoftPi, float);
DECLARE_SOA_COLUMN(NSigTpcTofPiSoftPi, nSigTpcTofPiSoftPi, float);
DECLARE_SOA_COLUMN(NSigTpcTofKaSoftPi, nSigTpcTofKaSoftPi, float);
// Dstar related variables
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);

// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

DECLARE_SOA_TABLE(HfCandDstLites, "AOD", "HFCANDDSTLITE",
                  hf_cand_dstar::Chi2PCAD0,
                  full::DecayLengthD0,
                  full::DecayLengthXYD0,
                  full::DecayLengthNormalisedD0,
                  full::DecayLengthXYNormalisedD0,
                  full::CpaD0,
                  full::CpaXYD0,
                  full::MaxNormalisedDeltaIPD0,
                  full::ImpactParameterProductD0,
                  full::CosThetaStarD0,
                  full::PtProng0,
                  full::PtProng1,
                  full::PtSoftPi,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand_dstar::ImpParamSoftPi,
                  full::ImpactParameterNormalised0,
                  full::ImpactParameterNormalised1,
                  full::ImpactParameterNormalisedSoftPi,
                  full::NSigTpcPi0,
                  full::NSigTpcKa0,
                  full::NSigTofPi0,
                  full::NSigTofKa0,
                  full::NSigTpcTofPi0,
                  full::NSigTpcTofKa0,
                  full::NSigTpcPi1,
                  full::NSigTpcKa1,
                  full::NSigTofPi1,
                  full::NSigTofKa1,
                  full::NSigTpcTofPi1,
                  full::NSigTpcTofKa1,
                  full::NSigTpcPiSoftPi,
                  full::NSigTpcKaSoftPi,
                  full::NSigTofPiSoftPi,
                  full::NSigTofKaSoftPi,
                  full::NSigTpcTofPiSoftPi,
                  full::NSigTpcTofKaSoftPi,
                  full::MD0,
                  full::PtD0,
                  full::EtaD0,
                  full::PhiD0,
                  full::YD0,
                  full::M,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::CandidateSelFlag,
                  hf_cand_dstar::FlagMcMatchRec,
                  hf_cand_dstar::OriginMcRec)

DECLARE_SOA_TABLE(HfCandDstFulls, "AOD", "HFCANDDSTFULL",
                  full::CollisionId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_cand_dstar::XSecondaryVertexD0,
                  hf_cand_dstar::YSecondaryVertexD0,
                  hf_cand_dstar::ZSecondaryVertexD0,
                  hf_cand_dstar::ErrorDecayLengthD0,
                  hf_cand_dstar::ErrorDecayLengthXYD0,
                  hf_cand_dstar::Chi2PCAD0,
                  full::DecayLengthD0,
                  full::DecayLengthXYD0,
                  full::DecayLengthNormalisedD0,
                  full::DecayLengthXYNormalisedD0,
                  full::CpaD0,
                  full::CpaXYD0,
                  full::MaxNormalisedDeltaIPD0,
                  full::ImpactParameterProductD0,
                  full::CosThetaStarD0,
                  full::PProng0,
                  full::PProng1,
                  full::PSoftPi,
                  full::PtProng0,
                  full::PtProng1,
                  full::PtSoftPi,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand_dstar::ImpParamSoftPi,
                  full::ImpactParameterNormalised0,
                  full::ImpactParameterNormalised1,
                  full::ImpactParameterNormalisedSoftPi,
                  hf_cand::ErrorImpactParameter0,
                  hf_cand::ErrorImpactParameter1,
                  hf_cand_dstar::ErrorImpParamSoftPi,
                  full::NSigTpcPi0,
                  full::NSigTpcKa0,
                  full::NSigTofPi0,
                  full::NSigTofKa0,
                  full::NSigTpcTofPi0,
                  full::NSigTpcTofKa0,
                  full::NSigTpcPi1,
                  full::NSigTpcKa1,
                  full::NSigTofPi1,
                  full::NSigTofKa1,
                  full::NSigTpcTofPi1,
                  full::NSigTpcTofKa1,
                  full::NSigTpcPiSoftPi,
                  full::NSigTpcKaSoftPi,
                  full::NSigTofPiSoftPi,
                  full::NSigTofKaSoftPi,
                  full::NSigTpcTofPiSoftPi,
                  full::NSigTpcTofKaSoftPi,
                  full::MD0,
                  full::PtD0,
                  full::PD0,
                  full::CtD0,
                  full::EtaD0,
                  full::PhiD0,
                  full::YD0,
                  full::M,
                  full::Pt,
                  full::P,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::E,
                  full::CandidateSelFlag,
                  hf_cand_dstar::FlagMcMatchRec,
                  hf_cand_dstar::OriginMcRec);

DECLARE_SOA_TABLE(HfCandDstFullEvs, "AOD", "HFCANDDSTFULLEV",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandDstFullPs, "AOD", "HFCANDDSTFULLP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  hf_cand_dstar::FlagMcMatchGen,
                  hf_cand_dstar::OriginMcGen);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorDstarToD0Pi {
  Produces<o2::aod::HfCandDstFulls> rowCandidateFull;
  Produces<o2::aod::HfCandDstFullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandDstFullPs> rowCandidateFullParticles;
  Produces<o2::aod::HfCandDstLites> rowCandidateLite;

  Configurable<bool> selectionFlagDstarToD0Pi{"selectionFlagDstarToD0Pi", true, "Selection Flag for D* decay to D0 & Pi"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  // parameters for production of training samples
  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to fill derived tables with signal for ML trainings"};
  Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background for ML trainings"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  using CandDstarWSelFlag = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstar, aod::HfSelDstarToD0Pi>>;
  using CandDstarWSelFlagMcRec = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstar, aod::HfSelDstarToD0Pi, aod::HfCandDstarMcRec>>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;
  using CandDstarMcGen = soa::Filtered<soa::Join<aod::McParticles, aod::HfCandDstarMcGen>>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstarToD0Pi;
  Filter filterMcGenMatching = nabs(aod::hf_cand_dstar::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_dstar::DecayType::DstarToD0Pi));

  Partition<CandDstarWSelFlagMcRec> reconstructedCandSig = nabs(aod::hf_cand_dstar::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_dstar::DecayType::DstarToD0Pi));
  Partition<CandDstarWSelFlagMcRec> reconstructedCandBkg = nabs(aod::hf_cand_dstar::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_dstar::DecayType::DstarToD0Pi));

  void init(InitContext const&)
  {
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

  template <bool doMc = false, typename T>
  void fillCandidateTable(const T& candidate)
  {
    int8_t flagMc{0};
    int8_t originMc{0};
    if constexpr (doMc) {
      flagMc = candidate.flagMcMatchRec();
      originMc = candidate.originMcRec();
    }

    auto prong0 = candidate.template prong0_as<TracksWPid>();
    auto prong1 = candidate.template prong1_as<TracksWPid>();
    auto prongSoftPi = candidate.template prongPi_as<TracksWPid>();

    float massD0{-1.f};
    float massDStar{-1.f};
    float cosThetaD0{-1.f};
    if (candidate.signSoftPi() > 0) {
      massD0 = candidate.invMassD0();
      massDStar = candidate.invMassDstar();
      cosThetaD0 = candidate.cosThetaStarD0();
    } else {
      massD0 = candidate.invMassD0Bar();
      massDStar = candidate.invMassAntiDstar();
      cosThetaD0 = candidate.cosThetaStarD0Bar();
    }

    if (fillCandidateLiteTable) {
      rowCandidateLite(
        candidate.chi2PCAD0(),
        candidate.decayLengthD0(),
        candidate.decayLengthXYD0(),
        candidate.decayLengthNormalisedD0(),
        candidate.decayLengthXYNormalisedD0(),
        candidate.cpaD0(),
        candidate.cpaXYD0(),
        candidate.deltaIPNormalisedMaxD0(),
        candidate.impactParameterProductD0(),
        cosThetaD0,
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.ptSoftPi(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impParamSoftPi(),
        candidate.impactParameterNormalised0(),
        candidate.impactParameterNormalised1(),
        candidate.normalisedImpParamSoftPi(),
        prong0.tpcNSigmaPi(),
        prong0.tpcNSigmaKa(),
        prong0.tofNSigmaPi(),
        prong0.tofNSigmaKa(),
        prong0.tpcTofNSigmaPi(),
        prong0.tpcTofNSigmaKa(),
        prong1.tpcNSigmaPi(),
        prong1.tpcNSigmaKa(),
        prong1.tofNSigmaPi(),
        prong1.tofNSigmaKa(),
        prong1.tpcTofNSigmaPi(),
        prong1.tpcTofNSigmaKa(),
        prongSoftPi.tpcNSigmaPi(),
        prongSoftPi.tpcNSigmaKa(),
        prongSoftPi.tofNSigmaPi(),
        prongSoftPi.tofNSigmaKa(),
        prongSoftPi.tpcTofNSigmaPi(),
        prongSoftPi.tpcTofNSigmaKa(),
        massD0,
        candidate.ptD0(),
        candidate.etaD0(),
        candidate.phiD0(),
        candidate.yD0(),
        massDStar,
        candidate.pt(),
        candidate.eta(),
        candidate.phi(),
        candidate.y(constants::physics::MassDStar),
        candidate.isSelDstarToD0Pi(),
        flagMc,
        originMc);
    } else {
      rowCandidateFull(
        candidate.collision().bcId(),
        candidate.collision().numContrib(),
        candidate.posX(),
        candidate.posY(),
        candidate.posZ(),
        candidate.xSecondaryVertexD0(),
        candidate.ySecondaryVertexD0(),
        candidate.zSecondaryVertexD0(),
        candidate.errorDecayLengthD0(),
        candidate.errorDecayLengthXYD0(),
        candidate.chi2PCAD0(),
        candidate.decayLengthD0(),
        candidate.decayLengthXYD0(),
        candidate.decayLengthNormalisedD0(),
        candidate.decayLengthXYNormalisedD0(),
        candidate.cpaD0(),
        candidate.cpaXYD0(),
        candidate.deltaIPNormalisedMaxD0(),
        candidate.impactParameterProductD0(),
        cosThetaD0,
        RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
        RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
        RecoDecay::p(candidate.pxSoftPi(), candidate.pySoftPi(), candidate.pzSoftPi()),
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.ptSoftPi(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impParamSoftPi(),
        candidate.impactParameterNormalised0(),
        candidate.impactParameterNormalised1(),
        candidate.normalisedImpParamSoftPi(),
        candidate.errorImpactParameter0(),
        candidate.errorImpactParameter1(),
        candidate.errorImpParamSoftPi(),
        prong0.tpcNSigmaPi(),
        prong0.tpcNSigmaKa(),
        prong0.tofNSigmaPi(),
        prong0.tofNSigmaKa(),
        prong0.tpcTofNSigmaPi(),
        prong0.tpcTofNSigmaKa(),
        prong1.tpcNSigmaPi(),
        prong1.tpcNSigmaKa(),
        prong1.tofNSigmaPi(),
        prong1.tofNSigmaKa(),
        prong1.tpcTofNSigmaPi(),
        prong1.tpcTofNSigmaKa(),
        prongSoftPi.tpcNSigmaPi(),
        prongSoftPi.tpcNSigmaKa(),
        prongSoftPi.tofNSigmaPi(),
        prongSoftPi.tofNSigmaKa(),
        prongSoftPi.tpcTofNSigmaPi(),
        prongSoftPi.tpcTofNSigmaKa(),
        massD0,
        candidate.ptD0(),
        candidate.pD0(),
        candidate.ctD0(),
        candidate.etaD0(),
        candidate.phiD0(),
        candidate.yD0(),
        massDStar,
        candidate.pt(),
        candidate.p(),
        candidate.eta(),
        candidate.phi(),
        candidate.y(constants::physics::MassDStar),
        candidate.e(constants::physics::MassDStar),
        candidate.isSelDstarToD0Pi(),
        flagMc,
        originMc);
    }
  }

  void processData(aod::Collisions const& collisions,
                   CandDstarWSelFlag const& candidates,
                   TracksWPid const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidates.size());
    } else {
      rowCandidateFull.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      if (downSampleBkgFactor < 1.) {
        float pseudoRndm = candidate.ptProng0() * 1000. - (int64_t)(candidate.ptProng0() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
      }
      fillCandidateTable(candidate);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorDstarToD0Pi, processData, "Process data", true);

  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const&,
                 CandDstarWSelFlagMcRec const& candidates,
                 CandDstarMcGen const& particles,
                 TracksWPid const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    if (fillOnlySignal) {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(reconstructedCandSig.size());
      } else {
        rowCandidateFull.reserve(reconstructedCandSig.size());
      }
      for (const auto& candidate : reconstructedCandSig) {
        fillCandidateTable<true>(candidate);
      }
    } else if (fillOnlyBackground) {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(reconstructedCandBkg.size());
      } else {
        rowCandidateFull.reserve(reconstructedCandBkg.size());
      }
      for (const auto& candidate : reconstructedCandBkg) {
        if (downSampleBkgFactor < 1.) {
          float pseudoRndm = candidate.ptProng0() * 1000. - (int64_t)(candidate.ptProng0() * 1000);
          if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
            continue;
          }
        }
        fillCandidateTable<true>(candidate);
      }
    } else {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(candidates.size());
      } else {
        rowCandidateFull.reserve(candidates.size());
      }
      for (const auto& candidate : candidates) {
        fillCandidateTable<true>(candidate);
      }
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (const auto& particle : particles) {
      rowCandidateFullParticles(
        particle.mcCollision().bcId(),
        particle.pt(),
        particle.eta(),
        particle.phi(),
        RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::constants::physics::MassDStar),
        particle.flagMcMatchGen(),
        particle.originMcGen());
    }
  }

  PROCESS_SWITCH(HfTreeCreatorDstarToD0Pi, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorDstarToD0Pi>(cfgc)};
}
