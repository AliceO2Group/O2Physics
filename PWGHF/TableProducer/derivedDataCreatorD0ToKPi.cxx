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
/// \brief Writer of the 2 prong candidates in the form of flat tables to be stored in TTrees.
///
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/DerivedTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
DECLARE_SOA_TABLE(HfCandD0Lites, "AOD", "HFCANDD0LITE",
                  hf_cand::Chi2PCA,
                  hf_cand_param::DecayLength,
                  hf_cand_param::DecayLengthXY,
                  hf_cand_param::DecayLengthNormalised,
                  hf_cand_param::DecayLengthXYNormalised,
                  hf_cand_param::PtProng0,
                  hf_cand_param::PtProng1,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand_param::ImpactParameterNormalised0,
                  hf_cand_param::ImpactParameterNormalised1,
                  hf_cand_param::NSigTpcPi0,
                  hf_cand_param::NSigTpcKa0,
                  hf_cand_param::NSigTofPi0,
                  hf_cand_param::NSigTofKa0,
                  hf_cand_param::NSigTpcTofPi0,
                  hf_cand_param::NSigTpcTofKa0,
                  hf_cand_param::NSigTpcPi1,
                  hf_cand_param::NSigTpcKa1,
                  hf_cand_param::NSigTofPi1,
                  hf_cand_param::NSigTofKa1,
                  hf_cand_param::NSigTpcTofPi1,
                  hf_cand_param::NSigTpcTofKa1,
                  hf_cand_flag::CandidateSelFlag,
                  hf_cand_analysis::M,
                  hf_cand_analysis::Pt,
                  hf_cand_param::Cpa,
                  hf_cand_param::CpaXY,
                  hf_cand_param::MaxNormalisedDeltaIP,
                  hf_cand_param::ImpactParameterProduct,
                  hf_cand_analysis::Eta,
                  hf_cand_analysis::Phi,
                  hf_cand_analysis::Y,
                  hf_cand_flag::FlagMc,
                  hf_cand_flag::OriginMcRec)

DECLARE_SOA_TABLE(HfCandD0Fulls, "AOD", "HFCANDD0FULL",
                  hf_index::CollisionId,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_cand::XSecondaryVertex,
                  hf_cand::YSecondaryVertex,
                  hf_cand::ZSecondaryVertex,
                  hf_cand::ErrorDecayLength,
                  hf_cand::ErrorDecayLengthXY,
                  hf_cand::Chi2PCA,
                  hf_cand::KfTopolChi2OverNdf,
                  hf_cand_param::RSecondaryVertex,
                  hf_cand_param::DecayLength,
                  hf_cand_param::DecayLengthXY,
                  hf_cand_param::DecayLengthNormalised,
                  hf_cand_param::DecayLengthXYNormalised,
                  hf_cand_param::ImpactParameterNormalised0,
                  hf_cand_param::PtProng0,
                  hf_cand_param::PProng0,
                  hf_cand_param::ImpactParameterNormalised1,
                  hf_cand_param::PtProng1,
                  hf_cand_param::PProng1,
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
                  hf_cand_param::NSigTpcPi0,
                  hf_cand_param::NSigTpcKa0,
                  hf_cand_param::NSigTofPi0,
                  hf_cand_param::NSigTofKa0,
                  hf_cand_param::NSigTpcTofPi0,
                  hf_cand_param::NSigTpcTofKa0,
                  hf_cand_param::NSigTpcPi1,
                  hf_cand_param::NSigTpcKa1,
                  hf_cand_param::NSigTofPi1,
                  hf_cand_param::NSigTofKa1,
                  hf_cand_param::NSigTpcTofPi1,
                  hf_cand_param::NSigTpcTofKa1,
                  hf_cand_flag::CandidateSelFlag,
                  hf_cand_analysis::M,
                  hf_cand_param::MaxNormalisedDeltaIP,
                  hf_cand_param::ImpactParameterProduct,
                  hf_cand_param::CosThetaStar,
                  hf_cand_analysis::Pt,
                  hf_cand_analysis::P,
                  hf_cand_param::Cpa,
                  hf_cand_param::CpaXY,
                  hf_cand_param::Ct,
                  hf_cand_analysis::Eta,
                  hf_cand_analysis::Phi,
                  hf_cand_analysis::Y,
                  hf_cand_analysis::E,
                  hf_cand_flag::FlagMc,
                  hf_cand_flag::OriginMcRec,
                  hf_index::Candidate2PId);

DECLARE_SOA_TABLE(HfCandD0FullEvs, "AOD", "HFCANDD0FULLEV",
                  hf_index::CollisionId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_collision::IsEventReject,
                  hf_collision::RunNumber);

DECLARE_SOA_TABLE(HfCandD0FullPs, "AOD", "HFCANDD0FULLP",
                  hf_index::McCollisionId,
                  hf_cand_analysis::Pt,
                  hf_cand_analysis::Eta,
                  hf_cand_analysis::Phi,
                  hf_cand_analysis::Y,
                  hf_cand_flag::FlagMc,
                  hf_cand_flag::OriginMcGen,
                  hf_index::McParticleId);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorD0ToKPi {
  Produces<o2::aod::HfCandD0Fulls> rowCandidateFull;
  Produces<o2::aod::HfCandD0FullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandD0FullPs> rowCandidateFullParticles;
  Produces<o2::aod::HfCandD0Lites> rowCandidateLite;

  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  // parameters for production of training samples
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
  void fillEvent(const T& collision, int isEventReject, int runNumber)
  {
    rowCandidateFullEvents(
      collision.globalIndex(),
      collision.numContrib(),
      collision.posX(),
      collision.posY(),
      collision.posZ(),
      isEventReject,
      runNumber);
  }

  template <typename T, typename U>
  auto fillTable(const T& candidate, const U& prong0, const U& prong1, int candFlag, double invMass, double cosThetaStar, double topoChi2,
                 double ct, double y, double e, int8_t flagMc, int8_t origin)
  {
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
        1 << candFlag,
        invMass,
        candidate.pt(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        candidate.impactParameterProduct(),
        candidate.eta(),
        candidate.phi(),
        y,
        flagMc,
        origin);
    } else {
      rowCandidateFull(
        candidate.collisionId(),
        candidate.posX(),
        candidate.posY(),
        candidate.posZ(),
        candidate.xSecondaryVertex(),
        candidate.ySecondaryVertex(),
        candidate.zSecondaryVertex(),
        candidate.errorDecayLength(),
        candidate.errorDecayLengthXY(),
        candidate.chi2PCA(),
        topoChi2,
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
        prong0.tpcTofNSigmaPi(),
        prong0.tpcTofNSigmaKa(),
        prong1.tpcNSigmaPi(),
        prong1.tpcNSigmaKa(),
        prong1.tofNSigmaPi(),
        prong1.tofNSigmaKa(),
        prong1.tpcTofNSigmaPi(),
        prong1.tpcTofNSigmaKa(),
        1 << candFlag,
        invMass,
        candidate.maxNormalisedDeltaIP(),
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
        origin,
        candidate.globalIndex());
    }
  }

  template <int reconstructionType, typename CandType>
  void processData(aod::Collisions const& collisions,
                   CandType const& candidates,
                   TracksWPid const&, aod::BCs const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, collision.bc().runNumber());
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
        fillTable(candidate, prong0, prong1, 0, massD0, hfHelper.cosThetaStarD0(candidate), topolChi2PerNdf, ctD, yD, eD, 0, 0);
      }
      if (candidate.isSelD0bar()) {
        fillTable(candidate, prong0, prong1, 1, massD0bar, hfHelper.cosThetaStarD0bar(candidate), topolChi2PerNdf, ctD, yD, eD, 0, 0);
      }
    }
  }

  void processDataWithDCAFitterN(aod::Collisions const& collisions,
                                 soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> const& candidates,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processData<aod::hf_cand::VertexerType::DCAFitter>(collisions, candidates, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processDataWithDCAFitterN, "Process data with DCAFitterN", true);

  void processDataWithKFParticle(aod::Collisions const& collisions,
                                 soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfCand2ProngKF, aod::HfSelD0>> const& candidates,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processData<aod::hf_cand::VertexerType::KfParticle>(collisions, candidates, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processDataWithKFParticle, "Process data with KFParticle", false);

  template <int reconstructionType, bool onlyBkg, bool onlySig, typename CandType>
  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const&,
                 CandType const& candidates,
                 MatchedGenCandidatesMc const& mcParticles,
                 TracksWPid const&,
                 aod::BCs const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, collision.bc().runNumber());
    }

    // Filling candidate properties
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidates.size());
    } else {
      rowCandidateFull.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      if constexpr (onlyBkg) {
        if (TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
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
        if (!TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          continue;
        }
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
        fillTable(candidate, prong0, prong1, 0, massD0, hfHelper.cosThetaStarD0(candidate), topolChi2PerNdf, ctD, yD, eD, candidate.flagMcMatchRec(), candidate.originMcRec());
      }
      if (candidate.isSelD0bar()) {
        fillTable(candidate, prong0, prong1, 1, massD0bar, hfHelper.cosThetaStarD0bar(candidate), topolChi2PerNdf, ctD, yD, eD, candidate.flagMcMatchRec(), candidate.originMcRec());
      }
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(mcParticles.size());
    for (const auto& particle : mcParticles) {
      if (TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        rowCandidateFullParticles(
          particle.mcCollisionId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::constants::physics::MassD0),
          particle.flagMcMatchGen(),
          particle.originMcGen(),
          particle.globalIndex());
      }
    }
  }

  void processMcWithDCAFitterOnlySig(aod::Collisions const& collisions,
                                     aod::McCollisions const& mcCollisions,
                                     SelectedCandidatesMc const&,
                                     MatchedGenCandidatesMc const& mcParticles,
                                     TracksWPid const& tracks,
                                     aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, false, true>(collisions, mcCollisions, reconstructedCandSig, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithDCAFitterOnlySig, "Process MC with DCAFitterN only for signals", false);

  void processMcWithDCAFitterOnlyBkg(aod::Collisions const& collisions,
                                     aod::McCollisions const& mcCollisions,
                                     SelectedCandidatesMc const&,
                                     MatchedGenCandidatesMc const& mcParticles,
                                     TracksWPid const& tracks,
                                     aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, true, false>(collisions, mcCollisions, reconstructedCandBkg, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithDCAFitterOnlyBkg, "Process MC with DCAFitterN only for background", false);

  void processMcWithDCAFitterAll(aod::Collisions const& collisions,
                                 aod::McCollisions const& mcCollisions,
                                 SelectedCandidatesMc const& candidates,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 TracksWPid const& tracks,
                                 aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, false, false>(collisions, mcCollisions, candidates, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithDCAFitterAll, "Process MC with DCAFitterN", false);

  void processMcWithKFParticleOnlySig(aod::Collisions const& collisions,
                                      aod::McCollisions const& mcCollisions,
                                      SelectedCandidatesMcKf const&,
                                      MatchedGenCandidatesMc const& mcParticles,
                                      TracksWPid const& tracks,
                                      aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, false, true>(collisions, mcCollisions, reconstructedCandSigKF, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithKFParticleOnlySig, "Process MC with KFParticle only for signals", false);

  void processMcWithKFParticleOnlyBkg(aod::Collisions const& collisions,
                                      aod::McCollisions const& mcCollisions,
                                      SelectedCandidatesMcKf const&,
                                      MatchedGenCandidatesMc const& mcParticles,
                                      TracksWPid const& tracks,
                                      aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, true, false>(collisions, mcCollisions, reconstructedCandBkgKF, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithKFParticleOnlyBkg, "Process MC with KFParticle only for background", false);

  void processMcWithKFParticleAll(aod::Collisions const& collisions,
                                  aod::McCollisions const& mcCollisions,
                                  SelectedCandidatesMcKf const& candidates,
                                  MatchedGenCandidatesMc const& mcParticles,
                                  TracksWPid const& tracks,
                                  aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, false, false>(collisions, mcCollisions, candidates, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithKFParticleAll, "Process MC with KFParticle", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorD0ToKPi>(cfgc));
  return workflow;
}
