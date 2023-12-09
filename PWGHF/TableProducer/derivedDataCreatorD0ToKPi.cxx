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

/// \file derivedDataCreatorD0ToKPi.cxx
/// \brief Producer of derived tables of D0 candidates, collisions and MC particles
///
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/DerivedTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
/// Table with basic candidate properties used in the analyses
// candidates for removal:
// E
DECLARE_SOA_TABLE(HfD0Bases, "AOD", "HFD0BASE",
                  hf_cand_base::Pt,
                  hf_cand_base::Eta,
                  hf_cand_base::Phi,
                  hf_cand_base::M,
                  hf_cand_base::P,
                  hf_cand_base::E,
                  hf_cand_base::Y);

/// Table with candidate properties used for selection
// candidates for removal:
// PxProng0, PyProng0, PzProng0,... (same for 1, 2), we can keep P, Pt, Phi instead
// XY: CpaXY, DecayLengthXY, ErrorDecayLengthXY
// normalised: DecayLengthNormalised, DecayLengthXYNormalised, ImpactParameterNormalised0
DECLARE_SOA_TABLE(HfD0Pars, "AOD", "HFD0PAR",
                  hf_cand::Chi2PCA,
                  hf_cand_par::DecayLength,
                  hf_cand_par::DecayLengthXY,
                  hf_cand_par::DecayLengthNormalised,
                  hf_cand_par::DecayLengthXYNormalised,
                  hf_cand_par::PtProng0,
                  hf_cand_par::PtProng1,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand_par::ImpactParameterNormalised0,
                  hf_cand_par::ImpactParameterNormalised1,
                  hf_cand_par::NSigTpcPi0,
                  hf_cand_par::NSigTpcKa0,
                  hf_cand_par::NSigTofPi0,
                  hf_cand_par::NSigTofKa0,
                  hf_cand_par::NSigTpcTofPi0,
                  hf_cand_par::NSigTpcTofKa0,
                  hf_cand_par::NSigTpcPi1,
                  hf_cand_par::NSigTpcKa1,
                  hf_cand_par::NSigTofPi1,
                  hf_cand_par::NSigTofKa1,
                  hf_cand_par::NSigTpcTofPi1,
                  hf_cand_par::NSigTpcTofKa1,
                  hf_cand_par::Cpa,
                  hf_cand_par::CpaXY,
                  hf_cand_par::MaxNormalisedDeltaIP,
                  hf_cand_par::ImpactParameterProduct);

DECLARE_SOA_TABLE(HfD0ParEs, "AOD", "HFD0PARE",
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_cand::XSecondaryVertex,
                  hf_cand::YSecondaryVertex,
                  hf_cand::ZSecondaryVertex,
                  hf_cand::ErrorDecayLength,
                  hf_cand::ErrorDecayLengthXY,
                  hf_cand::KfTopolChi2OverNdf,
                  hf_cand_par::RSecondaryVertex,
                  hf_cand_par::PProng0,
                  hf_cand_par::PProng1,
                  hf_cand::PxProng0,
                  hf_cand::PyProng0,
                  hf_cand::PzProng0,
                  hf_cand::PxProng1,
                  hf_cand::PyProng1,
                  hf_cand::PzProng1,
                  hf_cand::ErrorImpactParameter0,
                  hf_cand::ErrorImpactParameter1,
                  hf_cand_par::CosThetaStar,
                  hf_cand_par::Ct);

// Table with candidate selection flags
DECLARE_SOA_TABLE(HfD0Sels, "AOD", "HFD0SEL",
                  hf_cand_sel::CandidateSelFlag);

// Table with global indices for candidates
DECLARE_SOA_TABLE(HfD0Ids, "AOD", "HFD0ID",
                  hf_track_index::CollisionId,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id,
                  hf_cand_index::Candidate2PId);

// Table with candidate MC info
DECLARE_SOA_TABLE(HfD0Mcs, "AOD", "HFD0MC",
                  hf_cand_mc::FlagMc,
                  hf_cand_mc::OriginMcRec);

// Table with basic collision info
DECLARE_SOA_TABLE(HfD0CollBases, "AOD", "HFD0COLLBASE",
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_coll_base::IsEventReject,
                  hf_coll_base::RunNumber);

// Table with global indices for collisions
DECLARE_SOA_TABLE(HfD0CollIds, "AOD", "HFD0COLLID",
                  hf_cand_index::CollisionId);

// Table with MC particle info
DECLARE_SOA_TABLE(HfD0Ps, "AOD", "HFD0P",
                  hf_cand_base::Pt,
                  hf_cand_base::Eta,
                  hf_cand_base::Phi,
                  hf_cand_base::Y,
                  hf_cand_mc::FlagMc,
                  hf_cand_mc::OriginMcGen);

// Table with global indices for MC particles
DECLARE_SOA_TABLE(HfD0PIds, "AOD", "HFD0PID",
                  hf_cand_index::McCollisionId,
                  hf_cand_index::McParticleId);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfDerivedDataCreatorD0ToKPi {
  Produces<o2::aod::HfD0Bases> rowCandidateBase;
  Produces<o2::aod::HfD0Pars> rowCandidatePar;
  Produces<o2::aod::HfD0ParEs> rowCandidateParE;
  Produces<o2::aod::HfD0Sels> rowCandidateSel;
  Produces<o2::aod::HfD0Ids> rowCandidateId;
  Produces<o2::aod::HfD0Mcs> rowCandidateMc;
  Produces<o2::aod::HfD0CollBases> rowCollBase;
  Produces<o2::aod::HfD0CollIds> rowCollId;
  Produces<o2::aod::HfD0Ps> rowParticleBase;
  Produces<o2::aod::HfD0PIds> rowParticleId;

  // Switches for filling tables
  Configurable<bool> fillCandidateBase{"fillCandidateBase", true, "Fill candidate base properties"};
  Configurable<bool> fillCandidatePar{"fillCandidatePar", true, "Fill candidate parameters"};
  Configurable<bool> fillCandidateParE{"fillCandidateParE", true, "Fill candidate extended parameters"};
  Configurable<bool> fillCandidateSel{"fillCandidateSel", true, "Fill candidate selection flags"};
  Configurable<bool> fillCandidateId{"fillCandidateId", true, "Fill candidate indices"};
  Configurable<bool> fillCandidateMc{"fillCandidateMc", true, "Fill candidate MC info"};
  Configurable<bool> fillCollBase{"fillCollBase", true, "Fill collision base properties"};
  Configurable<bool> fillCollId{"fillCollId", true, "Fill collision indices"};
  Configurable<bool> fillParticleBase{"fillParticleBase", true, "Fill particle properties"};
  Configurable<bool> fillParticleId{"fillParticleId", true, "Fill particle indices"};
  // Parameters for production of training samples
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  HfHelper hfHelper;

  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPiExt, aod::TracksPidKaExt>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfCand2ProngMcRec, aod::HfSelD0>>;
  using SelectedCandidatesMcKf = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfCand2ProngKF, aod::HfCand2ProngMcRec, aod::HfSelD0>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_d0::isSelD0 >= 1 || aod::hf_sel_candidate_d0::isSelD0bar >= 1;
  Filter filterMcGenMatching = nabs(aod::hf_cand_2prong::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));

  Partition<SelectedCandidatesMc> reconstructedCandSig = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMc> reconstructedCandBkg = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMcKf> reconstructedCandSigKF = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));
  Partition<SelectedCandidatesMcKf> reconstructedCandBkgKF = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK));

  void init(InitContext const&)
  {
    std::array<bool, 8> doprocess{doprocessDataWithDCAFitterN, doprocessDataWithKFParticle, doprocessMcWithDCAFitterOnlySig, doprocessMcWithDCAFitterOnlyBkg, doprocessMcWithDCAFitterAll, doprocessMcWithKFParticleOnlySig, doprocessMcWithKFParticleOnlyBkg, doprocessMcWithKFParticleAll};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
  }

  template <typename T>
  void reserveTable(T& table, const Configurable<bool>& enabled, const uint64_t size)
  {
    if (enabled.value) {
      table.reserve(size);
    }
  };

  template <typename T>
  void fillTablesCollision(const T& collision, int isEventReject, int runNumber)
  {
    if (fillCollBase) {
      rowCollBase(
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        isEventReject,
        runNumber);
    }
    if (fillCollId) {
      rowCollId(
        collision.globalIndex());
    }
  }

  template <typename T, typename U>
  auto fillTablesCandidate(const T& candidate, const U& prong0, const U& prong1, int candFlag, double invMass, double cosThetaStar, double topoChi2,
                           double ct, double y, double e, int8_t flagMc, int8_t origin)
  {
    if (fillCandidateBase) {
      rowCandidateBase(
        candidate.pt(),
        candidate.eta(),
        candidate.phi(),
        invMass,
        candidate.p(),
        e,
        y);
    }
    if (fillCandidatePar) {
      rowCandidatePar(
        candidate.chi2PCA(),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXYNormalised(),
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impactParameterNormalised0(),
        candidate.impactParameterNormalised1(),
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
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        candidate.impactParameterProduct());
    }
    if (fillCandidateParE) {
      rowCandidateParE(
        candidate.posX(),
        candidate.posY(),
        candidate.posZ(),
        candidate.xSecondaryVertex(),
        candidate.ySecondaryVertex(),
        candidate.zSecondaryVertex(),
        candidate.errorDecayLength(),
        candidate.errorDecayLengthXY(),
        topoChi2,
        candidate.rSecondaryVertex(),
        RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
        RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
        candidate.pxProng0(),
        candidate.pyProng0(),
        candidate.pzProng0(),
        candidate.pxProng1(),
        candidate.pyProng1(),
        candidate.pzProng1(),
        candidate.errorImpactParameter0(),
        candidate.errorImpactParameter1(),
        cosThetaStar,
        ct);
    }
    if (fillCandidateSel) {
      rowCandidateSel(
        BIT(candFlag));
    }
    if (fillCandidateId) {
      rowCandidateId(
        candidate.collisionId(),
        candidate.prong0Id(),
        candidate.prong1Id(),
        candidate.globalIndex());
    }
    if (fillCandidateMc) {
      rowCandidateMc(
        flagMc,
        origin);
    }
  }

  template <typename T, typename U>
  void fillTablesParticle(const T& particle, U mass)
  {
    if (fillParticleBase) {
      rowParticleBase(
        particle.pt(),
        particle.eta(),
        particle.phi(),
        RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, mass),
        particle.flagMcMatchGen(),
        particle.originMcGen());
    }
    if (fillParticleId) {
      rowParticleId(
        particle.mcCollisionId(),
        particle.globalIndex());
    }
  }

  template <int reconstructionType, bool isMc, bool onlyBkg, bool onlySig, typename CandType>
  void processCandidates(aod::Collisions const& collisions,
                         CandType const& candidates,
                         TracksWPid const&,
                         aod::BCs const&)
  {
    // Fill collision properties
    auto sizeTableColl = collisions.size();
    reserveTable(rowCollBase, fillCollBase, sizeTableColl);
    reserveTable(rowCollId, fillCollId, sizeTableColl);
    for (const auto& collision : collisions) {
      fillTablesCollision(collision, 0, collision.bc().runNumber());
    }

    // Fill candidate properties
    auto sizeTableCand = candidates.size();
    reserveTable(rowCandidateBase, fillCandidateBase, sizeTableCand);
    reserveTable(rowCandidatePar, fillCandidatePar, sizeTableCand);
    reserveTable(rowCandidateParE, fillCandidateParE, sizeTableCand);
    reserveTable(rowCandidateSel, fillCandidateSel, sizeTableCand);
    reserveTable(rowCandidateId, fillCandidateId, sizeTableCand);
    if constexpr (isMc) {
      reserveTable(rowCandidateMc, fillCandidateMc, sizeTableCand);
    }
    int8_t flagMcRec, origin;
    for (const auto& candidate : candidates) {
      if constexpr (isMc) {
        flagMcRec = candidate.flagMcMatchRec();
        origin = candidate.originMcRec();
        if constexpr (onlyBkg) {
          if (TESTBIT(std::abs(flagMcRec), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
            continue;
          }
          if (downSampleBkgFactor < 1.) {
            float pseudoRndm = candidate.ptProng0() * 1000. - (int64_t)(candidate.ptProng0() * 1000);
            if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
              continue;
            }
          }
        }
        if constexpr (onlySig) {
          if (!TESTBIT(std::abs(flagMcRec), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
            continue;
          }
        }
      } else {
        flagMcRec = 0;
        origin = 0;
      }
      auto prong0 = candidate.template prong0_as<TracksWPid>();
      auto prong1 = candidate.template prong1_as<TracksWPid>();
      double yD = hfHelper.yD0(candidate);
      double eD = hfHelper.eD0(candidate);
      double ctD = hfHelper.ctD0(candidate);
      float massD0, massD0bar;
      float topolChi2PerNdf = -999.;
      if constexpr (reconstructionType == aod::hf_cand::VertexerType::KfParticle) {
        massD0 = candidate.kfGeoMassD0();
        massD0bar = candidate.kfGeoMassD0bar();
        topolChi2PerNdf = candidate.kfTopolChi2OverNdf();
      } else {
        massD0 = hfHelper.invMassD0ToPiK(candidate);
        massD0bar = hfHelper.invMassD0barToKPi(candidate);
      }
      if (candidate.isSelD0()) {
        fillTablesCandidate(candidate, prong0, prong1, 0, massD0, hfHelper.cosThetaStarD0(candidate), topolChi2PerNdf, ctD, yD, eD, flagMcRec, origin);
      }
      if (candidate.isSelD0bar()) {
        fillTablesCandidate(candidate, prong0, prong1, 1, massD0bar, hfHelper.cosThetaStarD0bar(candidate), topolChi2PerNdf, ctD, yD, eD, flagMcRec, origin);
      }
    }
  }

  template <typename ParticleType>
  void processMcParticles(ParticleType const& mcParticles)
  {
    // Fill MC particle properties
    auto sizeTablePart = mcParticles.size();
    reserveTable(rowParticleBase, fillParticleBase, sizeTablePart);
    reserveTable(rowParticleId, fillParticleId, sizeTablePart);
    for (const auto& particle : mcParticles) {
      if (!TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      fillTablesParticle(particle, o2::constants::physics::MassD0);
    }
  }

  void processDataWithDCAFitterN(aod::Collisions const& collisions,
                                 soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> const& candidates,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, false, false, false>(collisions, candidates, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processDataWithDCAFitterN, "Process data with DCAFitterN", true);

  void processDataWithKFParticle(aod::Collisions const& collisions,
                                 soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfCand2ProngKF, aod::HfSelD0>> const& candidates,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::KfParticle, false, false, false>(collisions, candidates, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processDataWithKFParticle, "Process data with KFParticle", false);

  void processMcWithDCAFitterOnlySig(aod::Collisions const& collisions,
                                     SelectedCandidatesMc const&,
                                     MatchedGenCandidatesMc const& mcParticles,
                                     TracksWPid const& tracks,
                                     aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, false, true>(collisions, reconstructedCandSig, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterOnlySig, "Process MC with DCAFitterN only for signals", false);

  void processMcWithDCAFitterOnlyBkg(aod::Collisions const& collisions,
                                     SelectedCandidatesMc const&,
                                     MatchedGenCandidatesMc const& mcParticles,
                                     TracksWPid const& tracks,
                                     aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, true, false>(collisions, reconstructedCandBkg, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterOnlyBkg, "Process MC with DCAFitterN only for background", false);

  void processMcWithDCAFitterAll(aod::Collisions const& collisions,
                                 SelectedCandidatesMc const& candidates,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, false, false>(collisions, candidates, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterAll, "Process MC with DCAFitterN", false);

  void processMcWithKFParticleOnlySig(aod::Collisions const& collisions,
                                      SelectedCandidatesMcKf const&,
                                      MatchedGenCandidatesMc const& mcParticles,
                                      TracksWPid const& tracks,
                                      aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, false, true>(collisions, reconstructedCandSigKF, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleOnlySig, "Process MC with KFParticle only for signals", false);

  void processMcWithKFParticleOnlyBkg(aod::Collisions const& collisions,
                                      SelectedCandidatesMcKf const&,
                                      MatchedGenCandidatesMc const& mcParticles,
                                      TracksWPid const& tracks,
                                      aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, true, false>(collisions, reconstructedCandBkgKF, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleOnlyBkg, "Process MC with KFParticle only for background", false);

  void processMcWithKFParticleAll(aod::Collisions const& collisions,
                                  SelectedCandidatesMcKf const& candidates,
                                  MatchedGenCandidatesMc const& mcParticles,
                                  TracksWPid const& tracks,
                                  aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, false, false>(collisions, candidates, tracks, bcs);
    processMcParticles(mcParticles);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleAll, "Process MC with KFParticle", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorD0ToKPi>(cfgc)};
}
