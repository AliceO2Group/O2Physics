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

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
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
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <array>
#include <cstdint>
#include <numeric>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
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
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(Cpa, cpa, float);
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);
DECLARE_SOA_COLUMN(CosThetaStar, cosThetaStar, float);
DECLARE_SOA_COLUMN(FlagMc, flagMc, int8_t);
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t);
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t);
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t); // is prompt or non-prompt, reco level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t); // is prompt or non-prompt, Gen level
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfCand2Prong, "_0");
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);
} // namespace full
namespace ml
{
DECLARE_SOA_COLUMN(BdtOutputBkg, bdtOutputBkg, float);
DECLARE_SOA_COLUMN(BdtOutputPrompt, bdtOutputPrompt, float);
DECLARE_SOA_COLUMN(BdtOutputNonPrompt, bdtOutputNonPrompt, float);
} // namespace ml

DECLARE_SOA_TABLE(HfCandD0Lites, "AOD", "HFCANDD0LITE",
                  hf_cand::Chi2PCA,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::PtProng0,
                  full::PtProng1,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  full::ImpactParameterNormalised0,
                  full::ImpactParameterNormalised1,
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
                  full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::Cpa,
                  full::CpaXY,
                  full::MaxNormalisedDeltaIP,
                  full::ImpactParameterProduct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::FlagMc,
                  full::FlagMcDecayChanRec,
                  full::OriginMcRec)

DECLARE_SOA_TABLE(HfCandD0Fulls, "AOD", "HFCANDD0FULL",
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
                  full::NSigTpcTofPi0,
                  full::NSigTpcTofKa0,
                  full::NSigTpcPi1,
                  full::NSigTpcKa1,
                  full::NSigTofPi1,
                  full::NSigTofKa1,
                  full::NSigTpcTofPi1,
                  full::NSigTpcTofKa1,
                  full::CandidateSelFlag,
                  full::M,
                  full::MaxNormalisedDeltaIP,
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
                  full::FlagMcDecayChanRec,
                  full::OriginMcRec);

DECLARE_SOA_TABLE(HfCandD0FullEvs, "AOD", "HFCANDD0FULLEV",
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandD0FullPs, "AOD", "HFCANDD0FULLP",
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::FlagMc,
                  full::FlagMcDecayChanGen,
                  full::OriginMcGen);

DECLARE_SOA_TABLE(HfCandD0Mls, "AOD", "HFCANDD0ML",
                  ml::BdtOutputBkg,
                  ml::BdtOutputNonPrompt,
                  ml::BdtOutputPrompt);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorD0ToKPi {
  Produces<o2::aod::HfCandD0Fulls> rowCandidateFull;
  Produces<o2::aod::HfCandD0FullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandD0FullPs> rowCandidateFullParticles;
  Produces<o2::aod::HfCandD0Lites> rowCandidateLite;
  Produces<o2::aod::HfCandD0Mls> rowCandidateMl;

  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  // parameters for production of training samples
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};
  Configurable<bool> fillCorrBkgs{"fillCorrBkgs", false, "Flag to fill derived tables with correlated background candidates"};

  // using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngMcRec, aod::HfSelD0>>;
  using SelectedCandidatesMcMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngMcRec, aod::HfSelD0, aod::HfMlD0>>;
  using SelectedCandidatesMcKf = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF, aod::HfCand2ProngMcRec, aod::HfSelD0>>;
  using SelectedCandidatesMcKfMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF, aod::HfCand2ProngMcRec, aod::HfSelD0, aod::HfMlD0>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_d0::isSelD0 >= 1 || aod::hf_sel_candidate_d0::isSelD0bar >= 1;
  Filter filterMcGenMatching = nabs(aod::hf_cand_2prong::flagMcMatchGen) == static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) || (fillCorrBkgs && (nabs(aod::hf_cand_2prong::flagMcMatchGen) != 0));

  Partition<SelectedCandidatesMc> reconstructedCandSig = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) || (fillCorrBkgs && nabs(aod::hf_cand_2prong::flagMcMatchRec) != 0);
  Partition<SelectedCandidatesMc> reconstructedCandBkg = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);
  Partition<SelectedCandidatesMcKf> reconstructedCandSigKF = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) || (fillCorrBkgs && nabs(aod::hf_cand_2prong::flagMcMatchRec) != 0);
  Partition<SelectedCandidatesMcKf> reconstructedCandBkgKF = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);

  Partition<SelectedCandidatesMcMl> reconstructedCandSigMl = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) || (fillCorrBkgs && nabs(aod::hf_cand_2prong::flagMcMatchRec) != 0);
  Partition<SelectedCandidatesMcMl> reconstructedCandBkgMl = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);
  Partition<SelectedCandidatesMcKfMl> reconstructedCandSigKFMl = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) || (fillCorrBkgs && nabs(aod::hf_cand_2prong::flagMcMatchRec) != 0);
  Partition<SelectedCandidatesMcKfMl> reconstructedCandBkgKFMl = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);

  void init(InitContext const&)
  {
    std::array<bool, 16> doprocess{doprocessDataWithDCAFitterN, doprocessDataWithKFParticle, doprocessMcWithDCAFitterOnlySig, doprocessMcWithDCAFitterOnlyBkg,
                                   doprocessMcWithDCAFitterAll, doprocessMcWithKFParticleOnlySig, doprocessMcWithKFParticleOnlyBkg, doprocessMcWithKFParticleAll,
                                   doprocessDataWithDCAFitterNMl, doprocessDataWithKFParticleMl, doprocessMcWithDCAFitterOnlySigMl, doprocessMcWithDCAFitterOnlyBkgMl,
                                   doprocessMcWithDCAFitterAllMl, doprocessMcWithKFParticleOnlySigMl, doprocessMcWithKFParticleOnlyBkgMl, doprocessMcWithKFParticleAllMl};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
  }

  template <typename T>
  void fillEvent(const T& collision, int isEventReject, int runNumber)
  {
    rowCandidateFullEvents(
      collision.numContrib(),
      collision.posX(),
      collision.posY(),
      collision.posZ(),
      isEventReject,
      runNumber);
  }

  template <bool ApplyMl, typename T>
  auto fillTable(const T& candidate, int candFlag, double invMass, double topoChi2,
                 double ct, double y, double e, int8_t flagMc, int8_t flagMcDecay, int8_t origin)
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
        candidate.nSigTpcPi0(),
        candidate.nSigTpcKa0(),
        candidate.nSigTofPi0(),
        candidate.nSigTofKa0(),
        candidate.tpcTofNSigmaPi0(),
        candidate.tpcTofNSigmaKa0(),
        candidate.nSigTpcPi1(),
        candidate.nSigTpcKa1(),
        candidate.nSigTofPi1(),
        candidate.nSigTofKa1(),
        candidate.tpcTofNSigmaPi1(),
        candidate.tpcTofNSigmaKa1(),
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
        flagMcDecay,
        origin);
    } else {
      double cosThetaStar = candFlag == 0 ? HfHelper::cosThetaStarD0(candidate) : HfHelper::cosThetaStarD0bar(candidate);
      rowCandidateFull(
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
        candidate.nSigTpcPi0(),
        candidate.nSigTpcKa0(),
        candidate.nSigTofPi0(),
        candidate.nSigTofKa0(),
        candidate.tpcTofNSigmaPi0(),
        candidate.tpcTofNSigmaKa0(),
        candidate.nSigTpcPi1(),
        candidate.nSigTpcKa1(),
        candidate.nSigTofPi1(),
        candidate.nSigTofKa1(),
        candidate.tpcTofNSigmaPi1(),
        candidate.tpcTofNSigmaKa1(),
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
        flagMcDecay,
        origin);
    }
    if constexpr (ApplyMl) {
      if (candFlag == 0) {
        rowCandidateMl(
          candidate.mlProbD0()[0],
          candidate.mlProbD0()[1],
          candidate.mlProbD0()[2]);
      } else if (candFlag == 1) {
        rowCandidateMl(
          candidate.mlProbD0bar()[0],
          candidate.mlProbD0bar()[1],
          candidate.mlProbD0bar()[2]);
      }
    }
  }

  template <int ReconstructionType, bool ApplyMl, typename CandType>
  void processData(aod::Collisions const& collisions,
                   CandType const& candidates,
                   aod::Tracks const&, aod::BCs const&)
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
    if constexpr (ApplyMl) {
      rowCandidateMl.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      if (downSampleBkgFactor < 1.) {
        float const pseudoRndm = candidate.ptProng0() * 1000. - static_cast<int64_t>(candidate.ptProng0() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
      }
      double const yD = HfHelper::yD0(candidate);
      double const eD = HfHelper::eD0(candidate);
      double const ctD = HfHelper::ctD0(candidate);
      float massD0, massD0bar;
      float topolChi2PerNdf = -999.;
      if constexpr (ReconstructionType == aod::hf_cand::VertexerType::KfParticle) {
        massD0 = candidate.kfGeoMassD0();
        massD0bar = candidate.kfGeoMassD0bar();
        topolChi2PerNdf = candidate.kfTopolChi2OverNdf();
      } else {
        massD0 = HfHelper::invMassD0ToPiK(candidate);
        massD0bar = HfHelper::invMassD0barToKPi(candidate);
      }
      if (candidate.isSelD0()) {
        fillTable<ApplyMl>(candidate, 0, massD0, topolChi2PerNdf, ctD, yD, eD, 0, 0, 0);
      }
      if (candidate.isSelD0bar()) {
        fillTable<ApplyMl>(candidate, 1, massD0bar, topolChi2PerNdf, ctD, yD, eD, 0, 0, 0);
      }
    }
  }

  void processDataWithDCAFitterN(aod::Collisions const& collisions,
                                 soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0>> const& candidates,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    processData<aod::hf_cand::VertexerType::DCAFitter, false>(collisions, candidates, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processDataWithDCAFitterN, "Process data with DCAFitterN", true);

  void processDataWithDCAFitterNMl(aod::Collisions const& collisions,
                                   soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0, aod::HfMlD0>> const& candidates,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    processData<aod::hf_cand::VertexerType::DCAFitter, true>(collisions, candidates, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processDataWithDCAFitterNMl, "Process data with DCAFitterN and ML", false);

  void processDataWithKFParticle(aod::Collisions const& collisions,
                                 soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF, aod::HfSelD0>> const& candidates,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    processData<aod::hf_cand::VertexerType::KfParticle, false>(collisions, candidates, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processDataWithKFParticle, "Process data with KFParticle", false);

  void processDataWithKFParticleMl(aod::Collisions const& collisions,
                                   soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF, aod::HfSelD0, aod::HfMlD0>> const& candidates,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    processData<aod::hf_cand::VertexerType::KfParticle, true>(collisions, candidates, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processDataWithKFParticleMl, "Process data with KFParticle and ML", false);

  template <int ReconstructionType, bool OnlyBkg, bool OnlySig, bool ApplyMl, typename CandType>
  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const&,
                 CandType const& candidates,
                 MatchedGenCandidatesMc const& mcParticles,
                 aod::Tracks const&,
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
    if constexpr (ApplyMl) {
      rowCandidateMl.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      if constexpr (OnlyBkg) {
        if ((std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) || (fillCorrBkgs && (candidate.flagMcMatchRec() != 0))) {
          continue;
        }
        if (downSampleBkgFactor < 1.) {
          float const pseudoRndm = candidate.ptProng0() * 1000. - static_cast<int64_t>(candidate.ptProng0() * 1000);
          if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
            continue;
          }
        }
      }
      if constexpr (OnlySig) {
        if (fillCorrBkgs && candidate.flagMcMatchRec() == 0) {
          continue;
        }
        if (!fillCorrBkgs && std::abs(candidate.flagMcMatchRec()) != o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          continue;
        }
      }
      double const yD = HfHelper::yD0(candidate);
      double const eD = HfHelper::eD0(candidate);
      double const ctD = HfHelper::ctD0(candidate);
      float massD0, massD0bar;
      float topolChi2PerNdf = -999.;
      if constexpr (ReconstructionType == aod::hf_cand::VertexerType::KfParticle) {
        massD0 = candidate.kfGeoMassD0();
        massD0bar = candidate.kfGeoMassD0bar();
        topolChi2PerNdf = candidate.kfTopolChi2OverNdf();
      } else {
        massD0 = HfHelper::invMassD0ToPiK(candidate);
        massD0bar = HfHelper::invMassD0barToKPi(candidate);
      }
      if (candidate.isSelD0()) {
        fillTable<ApplyMl>(candidate, 0, massD0, topolChi2PerNdf, ctD, yD, eD, candidate.flagMcMatchRec(), candidate.flagMcDecayChanRec(), candidate.originMcRec());
      }
      if (candidate.isSelD0bar()) {
        fillTable<ApplyMl>(candidate, 1, massD0bar, topolChi2PerNdf, ctD, yD, eD, candidate.flagMcMatchRec(), candidate.flagMcDecayChanRec(), candidate.originMcRec());
      }
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(mcParticles.size());
    for (const auto& particle : mcParticles) {
      if ((std::abs(particle.flagMcMatchGen()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) || (fillCorrBkgs && particle.flagMcMatchGen() != 0)) {
        rowCandidateFullParticles(
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(particle.pVector(), o2::constants::physics::MassD0),
          particle.flagMcMatchGen(),
          particle.flagMcDecayChanGen(),
          particle.originMcGen());
      }
    }
  }

  void processMcWithDCAFitterOnlySig(aod::Collisions const& collisions,
                                     aod::McCollisions const& mcCollisions,
                                     SelectedCandidatesMc const&,
                                     MatchedGenCandidatesMc const& mcParticles,
                                     aod::Tracks const& tracks,
                                     aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, false, true, false>(collisions, mcCollisions, reconstructedCandSig, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithDCAFitterOnlySig, "Process MC with DCAFitterN only for signals", false);

  void processMcWithDCAFitterOnlySigMl(aod::Collisions const& collisions,
                                       aod::McCollisions const& mcCollisions,
                                       SelectedCandidatesMcMl const&,
                                       MatchedGenCandidatesMc const& mcParticles,
                                       aod::Tracks const& tracks,
                                       aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, false, true, true>(collisions, mcCollisions, reconstructedCandSigMl, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithDCAFitterOnlySigMl, "Process MC with DCAFitterN only for signals and ML", false);

  void processMcWithDCAFitterOnlyBkg(aod::Collisions const& collisions,
                                     aod::McCollisions const& mcCollisions,
                                     SelectedCandidatesMc const&,
                                     MatchedGenCandidatesMc const& mcParticles,
                                     aod::Tracks const& tracks,
                                     aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, true, false, false>(collisions, mcCollisions, reconstructedCandBkg, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithDCAFitterOnlyBkg, "Process MC with DCAFitterN only for background", false);

  void processMcWithDCAFitterOnlyBkgMl(aod::Collisions const& collisions,
                                       aod::McCollisions const& mcCollisions,
                                       SelectedCandidatesMcMl const&,
                                       MatchedGenCandidatesMc const& mcParticles,
                                       aod::Tracks const& tracks,
                                       aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, true, false, true>(collisions, mcCollisions, reconstructedCandBkgMl, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithDCAFitterOnlyBkgMl, "Process MC with DCAFitterN only for background with ML", false);

  void processMcWithDCAFitterAll(aod::Collisions const& collisions,
                                 aod::McCollisions const& mcCollisions,
                                 SelectedCandidatesMc const& candidates,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, false, false, false>(collisions, mcCollisions, candidates, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithDCAFitterAll, "Process MC with DCAFitterN", false);

  void processMcWithDCAFitterAllMl(aod::Collisions const& collisions,
                                   aod::McCollisions const& mcCollisions,
                                   SelectedCandidatesMcMl const& candidates,
                                   MatchedGenCandidatesMc const& mcParticles,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, false, false, true>(collisions, mcCollisions, candidates, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithDCAFitterAllMl, "Process MC with DCAFitterN with ML", false);

  void processMcWithKFParticleOnlySig(aod::Collisions const& collisions,
                                      aod::McCollisions const& mcCollisions,
                                      SelectedCandidatesMcKf const&,
                                      MatchedGenCandidatesMc const& mcParticles,
                                      aod::Tracks const& tracks,
                                      aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, false, true, false>(collisions, mcCollisions, reconstructedCandSigKF, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithKFParticleOnlySig, "Process MC with KFParticle only for signals", false);

  void processMcWithKFParticleOnlySigMl(aod::Collisions const& collisions,
                                        aod::McCollisions const& mcCollisions,
                                        SelectedCandidatesMcKfMl const&,
                                        MatchedGenCandidatesMc const& mcParticles,
                                        aod::Tracks const& tracks,
                                        aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, false, true, true>(collisions, mcCollisions, reconstructedCandSigKFMl, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithKFParticleOnlySigMl, "Process MC with KFParticle only for signals with ML", false);

  void processMcWithKFParticleOnlyBkg(aod::Collisions const& collisions,
                                      aod::McCollisions const& mcCollisions,
                                      SelectedCandidatesMcKf const&,
                                      MatchedGenCandidatesMc const& mcParticles,
                                      aod::Tracks const& tracks,
                                      aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, true, false, false>(collisions, mcCollisions, reconstructedCandBkgKF, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithKFParticleOnlyBkg, "Process MC with KFParticle only for background", false);

  void processMcWithKFParticleOnlyBkgMl(aod::Collisions const& collisions,
                                        aod::McCollisions const& mcCollisions,
                                        SelectedCandidatesMcKfMl const&,
                                        MatchedGenCandidatesMc const& mcParticles,
                                        aod::Tracks const& tracks,
                                        aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, true, false, true>(collisions, mcCollisions, reconstructedCandBkgKFMl, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithKFParticleOnlyBkgMl, "Process MC with KFParticle only for background with ML", false);

  void processMcWithKFParticleAll(aod::Collisions const& collisions,
                                  aod::McCollisions const& mcCollisions,
                                  SelectedCandidatesMcKf const& candidates,
                                  MatchedGenCandidatesMc const& mcParticles,
                                  aod::Tracks const& tracks,
                                  aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, true, false, false>(collisions, mcCollisions, candidates, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithKFParticleAll, "Process MC with KFParticle", false);

  void processMcWithKFParticleAllMl(aod::Collisions const& collisions,
                                    aod::McCollisions const& mcCollisions,
                                    SelectedCandidatesMcKfMl const& candidates,
                                    MatchedGenCandidatesMc const& mcParticles,
                                    aod::Tracks const& tracks,
                                    aod::BCs const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, false, false, true>(collisions, mcCollisions, candidates, mcParticles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorD0ToKPi, processMcWithKFParticleAllMl, "Process MC with KFParticle with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorD0ToKPi>(cfgc));
  return workflow;
}
