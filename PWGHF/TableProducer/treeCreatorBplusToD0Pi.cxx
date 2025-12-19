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

/// \file treeCreatorBplusToD0Pi.cxx
/// \brief Writer of the 2 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
/// \note Extended from treeCreatorD0ToKPi.cxx
///
/// \author Antonio Palasciano <antonio.palasciano@ba.infn.it>, Università & INFN, Bari

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <array>
#include <cstdint>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_decay::hf_cand_beauty;

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
DECLARE_SOA_COLUMN(Ct, ct, float);
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
DECLARE_SOA_TABLE(HfCandBpFulls, "AOD", "HFCANDBPFULL",
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
                  hf_cand_2prong::FlagMcMatchRec,
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

DECLARE_SOA_TABLE(HfCandBpFullEvs, "AOD", "HFCANDBPFULLEV",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandBpFullPs, "AOD", "HFCANDBPFULLP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  hf_cand_2prong::FlagMcMatchRec,
                  hf_cand_2prong::OriginMcGen);

DECLARE_SOA_TABLE(HfCandBpLites, "AOD", "HFCANDBPLITE",
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
                  hf_cand_2prong::FlagMcMatchRec,
                  hf_cand_2prong::OriginMcRec);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorBplusToD0Pi {
  Produces<o2::aod::HfCandBpFulls> rowCandidateFull;
  Produces<o2::aod::HfCandBpFullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandBpFullPs> rowCandidateFullParticles;
  Produces<o2::aod::HfCandBpLites> rowCandidateLite;

  Configurable<int> selectionFlagBplus{"selectionBplus", 1, "Selection Flag for B+"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  // parameters for production of training samples
  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to fill derived tables with signal for ML trainings"};
  Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background for ML trainings"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCandBplus, aod::HfCandBplusMcRec, aod::HfSelBplusToD0Pi>>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::TracksPidKa>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= selectionFlagBplus;

  Partition<SelectedCandidatesMc> recSig = nabs(aod::hf_cand_bplus::flagMcMatchRec) == static_cast<int8_t>(DecayChannelMain::BplusToD0Pi);
  Partition<SelectedCandidatesMc> recBg = nabs(aod::hf_cand_bplus::flagMcMatchRec) != static_cast<int8_t>(DecayChannelMain::BplusToD0Pi);

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

  template <bool DoMc = false, typename T, typename U>
  void fillCandidateTable(const T& candidate, const U& prong1)
  {
    int8_t flagMc = 0;
    int8_t originMc = 0;
    if constexpr (DoMc) {
      flagMc = candidate.flagMcMatchRec();
      originMc = candidate.originMcRec();
    }
    auto d0Cand = candidate.prong0();
    // adding D0 daughters to the table
    auto d0Daughter0 = d0Cand.template prong0_as<TracksWPid>();
    auto d0Daughter1 = d0Cand.template prong1_as<TracksWPid>();
    auto invMassD0 = 0.;
    if (prong1.signed1Pt() > 0) {
      invMassD0 = HfHelper::invMassD0barToKPi(d0Cand);
    } else if (prong1.signed1Pt() < 0) {
      invMassD0 = HfHelper::invMassD0ToPiK(d0Cand);
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
        candidate.isSelBplusToD0Pi(),
        HfHelper::invMassBplusToD0Pi(candidate),
        candidate.pt(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        candidate.eta(),
        candidate.phi(),
        HfHelper::yBplus(candidate),
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
        candidate.isSelBplusToD0Pi(),
        HfHelper::invMassBplusToD0Pi(candidate),
        candidate.pt(),
        candidate.p(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        HfHelper::ctBplus(candidate),
        candidate.eta(),
        candidate.phi(),
        HfHelper::yBplus(candidate),
        HfHelper::eBplus(candidate),
        flagMc,
        invMassD0,
        d0Cand.ptProng0(),
        d0Cand.ptProng1(),
        HfHelper::yD0(d0Cand),
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
  }

  void processData(aod::Collisions const& collisions,
                   soa::Filtered<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>> const& candidates,
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
        float const pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
      }
      auto prong1 = candidate.prong1_as<TracksWPid>();
      fillCandidateTable(candidate, prong1);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorBplusToD0Pi, processData, "Process data", true);

  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const&,
                 SelectedCandidatesMc const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandBplusMcGen> const& particles,
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
        float const pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
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
      if (std::abs(particle.flagMcMatchGen()) == DecayChannelMain::BplusToD0Pi) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.mcCollisionId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::constants::physics::MassBPlus),
          particle.flagMcMatchGen());
      }
    }
  }

  PROCESS_SWITCH(HfTreeCreatorBplusToD0Pi, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorBplusToD0Pi>(cfgc)};
}
