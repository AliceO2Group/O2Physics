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

/// \file treeCreatorB0ToDPi.cxx
/// \brief Writer of B0 → D- π+ candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug, local optimization of analysis on small samples or ML training.
///        In this file are defined and filled the output tables
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);                     //! Radius of secondary vertex (cm)
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);                                     //! Transverse momentum of prong0 (GeV/c)
DECLARE_SOA_COLUMN(PProng0, pProng0, float);                                       //! Momentum of prong0 (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float); //! Normalised impact parameter of prong0
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);                                     //! Transverse momentum of prong1 (GeV/c)
DECLARE_SOA_COLUMN(PProng1, pProng1, float);                                       //! Momentum of prong1 (in GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float); //! Normalised impact parameter of prong1
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int);                       //! Selection flag of candidate (output of candidateSelector)
DECLARE_SOA_COLUMN(M, m, float);                                                   //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                                                 //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(P, p, float);                                                   //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                                   //! Rapidity of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);                                               //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                                               //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(E, e, float);                                                   //! Energy of candidate (GeV)
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);                                 //! TPC Nsigma separation for prong1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);                                 //! TOF Nsigma separation for prong1 with pion mass hypothesis
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                               //! Decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                           //! Transverse decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);           //! Normalised decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);       //! Normalised transverse decay length of candidate
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                               //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                           //! Cosine pointing angle of candidate in transverse plane
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);             //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(Ct, ct, float);                                                 //! Proper lifetime times c of candidate (cm)
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int); //! Event rejection flag
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);         //! Run number
} // namespace full

DECLARE_SOA_TABLE(HfCandB0Lites, "AOD", "HFCANDB0LITE",
                  hf_cand::Chi2PCA,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::PtProng0,
                  full::PtProng1,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  full::NSigTpcPi1,
                  full::NSigTofPi1,
                  full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::Cpa,
                  full::CpaXY,
                  full::MaxNormalisedDeltaIP,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcRec);

DECLARE_SOA_TABLE(HfCandB0Fulls, "AOD", "HFCANDB0FULL",
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
                  full::NSigTpcPi1,
                  full::NSigTofPi1,
                  full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::P,
                  full::Cpa,
                  full::CpaXY,
                  full::MaxNormalisedDeltaIP,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::E,
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcRec);

DECLARE_SOA_TABLE(HfCandB0FullEvs, "AOD", "HFCANDB0FULLEV",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandB0FullPs, "AOD", "HFCANDB0FULLP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcGen);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorB0ToDPi {
  Produces<o2::aod::HfCandB0Fulls> rowCandidateFull;
  Produces<o2::aod::HfCandB0FullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandB0FullPs> rowCandidateFullParticles;
  Produces<o2::aod::HfCandB0Lites> rowCandidateLite;

  Configurable<int> selectionFlagB0{"selectionB0", 1, "Selection Flag for B0"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  // parameters for production of training samples
  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to fill derived tables with signal for ML trainings"};
  Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background for ML trainings"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  HfHelper hfHelper;

  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCandB0, aod::HfCandB0McRec, aod::HfSelB0ToDPi>>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_b0::isSelB0ToDPi >= selectionFlagB0;

  Partition<SelectedCandidatesMc> recSig = nabs(aod::hf_cand_b0::flagMcMatchRec) == (int8_t)BIT(aod::hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi);
  Partition<SelectedCandidatesMc> recBg = nabs(aod::hf_cand_b0::flagMcMatchRec) != (int8_t)BIT(aod::hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi);

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

  template <bool doMc = false, typename T, typename U>
  void fillCandidateTable(const T& candidate, const U& prong1)
  {
    int8_t flagMc = 0;
    int8_t originMc = 0;
    if constexpr (doMc) {
      flagMc = candidate.flagMcMatchRec();
      originMc = candidate.originMcRec();
    }
    if (fillCandidateLiteTable) {
      rowCandidateLite(
        candidate.chi2PCA(),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXYNormalised(),
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        prong1.tpcNSigmaPi(),
        prong1.tofNSigmaPi(),
        candidate.isSelB0ToDPi(),
        hfHelper.invMassB0ToDPi(candidate),
        candidate.pt(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        candidate.eta(),
        candidate.phi(),
        hfHelper.yB0(candidate),
        flagMc,
        originMc);
    } else {
      rowCandidateFull(
        prong1.collision().bcId(),
        prong1.collision().numContrib(),
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
        prong1.tpcNSigmaPi(),
        prong1.tofNSigmaPi(),
        candidate.isSelB0ToDPi(),
        hfHelper.invMassB0ToDPi(candidate),
        candidate.pt(),
        candidate.p(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        hfHelper.ctB0(candidate),
        candidate.eta(),
        candidate.phi(),
        hfHelper.yB0(candidate),
        hfHelper.eB0(candidate),
        flagMc,
        originMc);
    }
  }

  void processData(aod::Collisions const& collisions,
                   soa::Filtered<soa::Join<aod::HfCandB0, aod::HfSelB0ToDPi>> const& candidates,
                   TracksWPid const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      if (fillOnlyBackground) {
        float pseudoRndm = candidate.ptProng1() * 1000. - (int64_t)(candidate.ptProng1() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
      }
      auto prong1 = candidate.prong1_as<TracksWPid>();
      fillCandidateTable(candidate, prong1);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorB0ToDPi, processData, "Process data", true);

  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const&,
                 SelectedCandidatesMc const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandB0McGen> const& particles,
                 TracksWPid const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    if (fillOnlySignal) {
      rowCandidateFull.reserve(recSig.size());
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(recSig.size());
      }
      for (const auto& candidate : recSig) {
        auto prong1 = candidate.prong1_as<TracksWPid>();
        fillCandidateTable<true>(candidate, prong1);
      }
    } else if (fillOnlyBackground) {
      rowCandidateFull.reserve(recBg.size());
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(recBg.size());
      }
      for (const auto& candidate : recBg) {
        float pseudoRndm = candidate.ptProng1() * 1000. - (int64_t)(candidate.ptProng1() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
        auto prong1 = candidate.prong1_as<TracksWPid>();
        fillCandidateTable<true>(candidate, prong1);
      }
    } else {
      rowCandidateFull.reserve(candidates.size());
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(candidates.size());
      }
      for (const auto& candidate : candidates) {
        auto prong1 = candidate.prong1_as<TracksWPid>();
        fillCandidateTable<true>(candidate, prong1);
      }
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (const auto& particle : particles) {
      if (TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi)) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::constants::physics::MassB0),
          particle.flagMcMatchGen(),
          particle.originMcGen());
      }
    }
  }

  PROCESS_SWITCH(HfTreeCreatorB0ToDPi, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorB0ToDPi>(cfgc)};
}
