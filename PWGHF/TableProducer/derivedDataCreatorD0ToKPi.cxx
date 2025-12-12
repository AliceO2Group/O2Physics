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
/// \note Based on treeCreatorD0ToKPi.cxx
///
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGHF/Utils/utilsDerivedData.h"
#include "PWGHF/Utils/utilsPid.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <map>
#include <numeric>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::pid_tpc_tof_utils;
using namespace o2::analysis::hf_derived;

/// Writes the full information in an output TTree
struct HfDerivedDataCreatorD0ToKPi {
  HfProducesDerivedData<
    o2::aod::HfD0Bases,
    o2::aod::HfD0CollBases,
    o2::aod::HfD0CollIds,
    o2::aod::HfD0McCollBases,
    o2::aod::HfD0McCollIds,
    o2::aod::HfD0McRCollIds,
    o2::aod::HfD0PBases,
    o2::aod::HfD0PIds>
    rowsCommon;
  // Candidates
  Produces<o2::aod::HfD0Pars> rowCandidatePar;
  Produces<o2::aod::HfD0ParEs> rowCandidateParE;
  Produces<o2::aod::HfD0Sels> rowCandidateSel;
  Produces<o2::aod::HfD0Mls> rowCandidateMl;
  Produces<o2::aod::HfD0Ids> rowCandidateId;
  Produces<o2::aod::HfD0Mcs> rowCandidateMc;

  // Switches for filling tables
  HfConfigurableDerivedData confDerData;
  Configurable<bool> fillCandidatePar{"fillCandidatePar", true, "Fill candidate parameters"};
  Configurable<bool> fillCandidateParE{"fillCandidateParE", true, "Fill candidate extended parameters"};
  Configurable<bool> fillCandidateSel{"fillCandidateSel", true, "Fill candidate selection flags"};
  Configurable<bool> fillCandidateMl{"fillCandidateMl", true, "Fill candidate selection ML scores"};
  Configurable<bool> fillCandidateId{"fillCandidateId", true, "Fill original indices from the candidate table"};
  Configurable<bool> fillCandidateMc{"fillCandidateMc", true, "Fill candidate MC info"};
  // Parameters for production of training samples
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  SliceCache cache;
  static constexpr double Mass{o2::constants::physics::MassD0};

  using CollisionsWCentMult = soa::Join<aod::Collisions, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  using CollisionsWMcCentMult = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMultZeqs>;
  // using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;
  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0>>;
  using SelectedCandidatesKf = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF, aod::HfSelD0>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngMcRec, aod::HfSelD0>>;
  using SelectedCandidatesMcKf = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF, aod::HfCand2ProngMcRec, aod::HfSelD0>>;
  using SelectedCandidatesMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0, aod::HfMlD0>>;
  using SelectedCandidatesKfMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF, aod::HfSelD0, aod::HfMlD0>>;
  using SelectedCandidatesMcMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngMcRec, aod::HfSelD0, aod::HfMlD0>>;
  using SelectedCandidatesMcKfMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF, aod::HfCand2ProngMcRec, aod::HfSelD0, aod::HfMlD0>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>;
  using TypeMcCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_d0::isSelD0 >= 1 || aod::hf_sel_candidate_d0::isSelD0bar >= 1;
  Filter filterMcGenMatching = nabs(aod::hf_cand_2prong::flagMcMatchGen) == static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);

  Preslice<SelectedCandidates> candidatesPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesKf> candidatesKfPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMc> candidatesMcPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMcKf> candidatesMcKfPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMl> candidatesMlPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesKfMl> candidatesKfMlPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMcMl> candidatesMcMlPerCollision = aod::hf_cand::collisionId;
  Preslice<SelectedCandidatesMcKfMl> candidatesMcKfMlPerCollision = aod::hf_cand::collisionId;
  Preslice<MatchedGenCandidatesMc> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  // trivial partitions for all candidates to allow "->sliceByCached" inside processCandidates
  Partition<SelectedCandidates> candidatesAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesKf> candidatesKfAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesMc> candidatesMcAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesMcKf> candidatesMcKfAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesMl> candidatesMlAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesKfMl> candidatesKfMlAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesMcMl> candidatesMcMlAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  Partition<SelectedCandidatesMcKfMl> candidatesMcKfMlAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  // partitions for signal and background
  Partition<SelectedCandidatesMc> candidatesMcSig = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);
  Partition<SelectedCandidatesMc> candidatesMcBkg = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);
  Partition<SelectedCandidatesMcKf> candidatesMcKfSig = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);
  Partition<SelectedCandidatesMcKf> candidatesMcKfBkg = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);
  Partition<SelectedCandidatesMcMl> candidatesMcMlSig = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);
  Partition<SelectedCandidatesMcMl> candidatesMcMlBkg = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);
  Partition<SelectedCandidatesMcKfMl> candidatesMcKfMlSig = nabs(aod::hf_cand_2prong::flagMcMatchRec) == static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);
  Partition<SelectedCandidatesMcKfMl> candidatesMcKfMlBkg = nabs(aod::hf_cand_2prong::flagMcMatchRec) != static_cast<int8_t>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);

  void init(InitContext const&)
  {
    std::array<bool, 17> doprocess{doprocessDataWithDCAFitterN, doprocessDataWithKFParticle, doprocessMcWithDCAFitterSig, doprocessMcWithDCAFitterBkg, doprocessMcWithDCAFitterAll, doprocessMcWithKFParticleSig, doprocessMcWithKFParticleBkg, doprocessMcWithKFParticleAll,
                                   doprocessDataWithDCAFitterNMl, doprocessDataWithKFParticleMl, doprocessMcWithDCAFitterMlSig, doprocessMcWithDCAFitterMlBkg, doprocessMcWithDCAFitterMlAll, doprocessMcWithKFParticleMlSig, doprocessMcWithKFParticleMlBkg, doprocessMcWithKFParticleMlAll, doprocessMcGenOnly};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
    rowsCommon.init(confDerData);
  }

  template <typename T>
  void fillTablesCandidate(const T& candidate, int candFlag, double invMass, double cosThetaStar, double topoChi2,
                           double ct, double y, int8_t flagMc, int8_t origin, const std::vector<float>& mlScores)
  {
    rowsCommon.fillTablesCandidate(candidate, invMass, y);
    if (fillCandidatePar) {
      std::array<std::array<std::array<float, 3>, 2>, 2> sigmas{}; // PID nSigma [Expected][Hypothesis][TPC/TOF/TPC+TOF]
      if (candFlag == 0) {
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion], candidate, 0, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon], candidate, 0, Ka)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion], candidate, 1, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon], candidate, 1, Ka)
      } else if (candFlag == 1) {
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion], candidate, 1, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon], candidate, 1, Ka)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion], candidate, 0, Pi)
        GET_N_SIGMA_PRONG(sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon], candidate, 0, Ka)
      }
      rowCandidatePar(
        candidate.chi2PCA(),
        candidate.cpa(),
        candidate.cpaXY(),
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
        sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion][0],
        sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion][1],
        sigmas[HfProngSpecies::Pion][HfProngSpecies::Pion][2],
        sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon][0],
        sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon][1],
        sigmas[HfProngSpecies::Pion][HfProngSpecies::Kaon][2],
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion][0],
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion][1],
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Pion][2],
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon][0],
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon][1],
        sigmas[HfProngSpecies::Kaon][HfProngSpecies::Kaon][2],
        candidate.maxNormalisedDeltaIP(),
        candidate.impactParameterProduct());
    }
    if (fillCandidateParE) {
      rowCandidateParE(
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
    if (fillCandidateMl) {
      rowCandidateMl(
        mlScores);
    }
    if (fillCandidateId) {
      rowCandidateId(
        candidate.collisionId(),
        candidate.prong0Id(),
        candidate.prong1Id());
    }
    if (fillCandidateMc) {
      rowCandidateMc(
        flagMc,
        origin);
    }
  }

  template <aod::hf_cand::VertexerType ReconstructionType, bool IsMl, bool IsMc, bool OnlyBkg, bool OnlySig, typename CollType, typename CandType>
  void processCandidates(CollType const& collisions,
                         Partition<CandType>& candidates,
                         aod::Tracks const&,
                         aod::BCs const&)
  {
    // Fill collision properties
    if constexpr (IsMc) {
      if (confDerData.fillMcRCollId) {
        rowsCommon.matchedCollisions.clear();
      }
    }
    auto sizeTableColl = collisions.size();
    rowsCommon.reserveTablesColl(sizeTableColl);
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candidatesThisColl = candidates->sliceByCached(aod::hf_cand::collisionId, thisCollId, cache); // FIXME
      auto sizeTableCand = candidatesThisColl.size();
      LOGF(debug, "Rec. collision %d has %d candidates", thisCollId, sizeTableCand);
      // Skip collisions without HF candidates (and without HF particles in matched MC collisions if saving indices of reconstructed collisions matched to MC collisions)
      bool mcCollisionHasMcParticles{false};
      if constexpr (IsMc) {
        mcCollisionHasMcParticles = confDerData.fillMcRCollId && collision.has_mcCollision() && rowsCommon.hasMcParticles[collision.mcCollisionId()];
        LOGF(debug, "Rec. collision %d has MC collision %d with MC particles? %s", thisCollId, collision.mcCollisionId(), mcCollisionHasMcParticles ? "yes" : "no");
      }
      if (sizeTableCand == 0 && (!confDerData.fillMcRCollId || !mcCollisionHasMcParticles)) {
        LOGF(debug, "Skipping rec. collision %d", thisCollId);
        continue;
      }
      LOGF(debug, "Filling rec. collision %d at derived index %d", thisCollId, rowsCommon.rowCollBase.lastIndex() + 1);
      rowsCommon.fillTablesCollision<IsMc>(collision);

      // Fill candidate properties
      rowsCommon.reserveTablesCandidates(sizeTableCand);
      reserveTable(rowCandidatePar, fillCandidatePar, sizeTableCand);
      reserveTable(rowCandidateParE, fillCandidateParE, sizeTableCand);
      reserveTable(rowCandidateSel, fillCandidateSel, sizeTableCand);
      reserveTable(rowCandidateMl, fillCandidateMl, sizeTableCand);
      reserveTable(rowCandidateId, fillCandidateId, sizeTableCand);
      if constexpr (IsMc) {
        reserveTable(rowCandidateMc, fillCandidateMc, sizeTableCand);
      }
      int8_t flagMcRec = 0, origin = 0;
      for (const auto& candidate : candidatesThisColl) {
        if constexpr (IsMc) {
          flagMcRec = candidate.flagMcMatchRec();
          origin = candidate.originMcRec();
          if constexpr (OnlyBkg) {
            if (std::abs(flagMcRec) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
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
            if (std::abs(flagMcRec) != o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
              continue;
            }
          }
        } else {
          if (downSampleBkgFactor < 1.) {
            float const pseudoRndm = candidate.ptProng0() * 1000. - static_cast<int64_t>(candidate.ptProng0() * 1000);
            if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
              continue;
            }
          }
        }

        double const ct = HfHelper::ctD0(candidate);
        double const y = HfHelper::yD0(candidate);
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
        std::vector<float> mlScoresD0, mlScoresD0bar;
        if constexpr (IsMl) {
          std::copy(candidate.mlProbD0().begin(), candidate.mlProbD0().end(), std::back_inserter(mlScoresD0));
          std::copy(candidate.mlProbD0bar().begin(), candidate.mlProbD0bar().end(), std::back_inserter(mlScoresD0bar));
        }
        if (candidate.isSelD0()) {
          fillTablesCandidate(candidate, 0, massD0, HfHelper::cosThetaStarD0(candidate), topolChi2PerNdf, ct, y, flagMcRec, origin, mlScoresD0);
        }
        if (candidate.isSelD0bar()) {
          fillTablesCandidate(candidate, 1, massD0bar, HfHelper::cosThetaStarD0bar(candidate), topolChi2PerNdf, ct, y, flagMcRec, origin, mlScoresD0bar);
        }
      }
    }
  }
  void processDataWithDCAFitterN(CollisionsWCentMult const& collisions,
                                 SelectedCandidates const&,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, false, false, false, false>(collisions, candidatesAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processDataWithDCAFitterN, "Process data with DCAFitterN", true);

  void processDataWithKFParticle(CollisionsWCentMult const& collisions,
                                 SelectedCandidatesKf const&,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::KfParticle, false, false, false, false>(collisions, candidatesKfAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processDataWithKFParticle, "Process data with KFParticle", false);

  void processMcWithDCAFitterSig(CollisionsWMcCentMult const& collisions,
                                 SelectedCandidatesMc const&,
                                 TypeMcCollisions const& mcCollisions,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, false, true, false, true>(collisions, candidatesMcSig, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterSig, "Process MC with DCAFitterN only for signals", false);

  void processMcWithDCAFitterBkg(CollisionsWMcCentMult const& collisions,
                                 SelectedCandidatesMc const&,
                                 TypeMcCollisions const& mcCollisions,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, false, true, true, false>(collisions, candidatesMcBkg, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterBkg, "Process MC with DCAFitterN only for background", false);

  void processMcWithDCAFitterAll(CollisionsWMcCentMult const& collisions,
                                 SelectedCandidatesMc const&,
                                 TypeMcCollisions const& mcCollisions,
                                 MatchedGenCandidatesMc const& mcParticles,
                                 aod::Tracks const& tracks,
                                 aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, false, true, false, false>(collisions, candidatesMcAll, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterAll, "Process MC with DCAFitterN", false);

  void processMcWithKFParticleSig(CollisionsWMcCentMult const& collisions,
                                  SelectedCandidatesMcKf const&,
                                  TypeMcCollisions const& mcCollisions,
                                  MatchedGenCandidatesMc const& mcParticles,
                                  aod::Tracks const& tracks,
                                  aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::KfParticle, false, true, false, true>(collisions, candidatesMcKfSig, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleSig, "Process MC with KFParticle only for signals", false);

  void processMcWithKFParticleBkg(CollisionsWMcCentMult const& collisions,
                                  SelectedCandidatesMcKf const&,
                                  TypeMcCollisions const& mcCollisions,
                                  MatchedGenCandidatesMc const& mcParticles,
                                  aod::Tracks const& tracks,
                                  aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::KfParticle, false, true, true, false>(collisions, candidatesMcKfBkg, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleBkg, "Process MC with KFParticle only for background", false);

  void processMcWithKFParticleAll(CollisionsWMcCentMult const& collisions,
                                  SelectedCandidatesMcKf const&,
                                  TypeMcCollisions const& mcCollisions,
                                  MatchedGenCandidatesMc const& mcParticles,
                                  aod::Tracks const& tracks,
                                  aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::KfParticle, false, true, false, false>(collisions, candidatesMcKfAll, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleAll, "Process MC with KFParticle", false);

  // ML versions

  void processDataWithDCAFitterNMl(CollisionsWCentMult const& collisions,
                                   SelectedCandidatesMl const&,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, false, false, false>(collisions, candidatesMlAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processDataWithDCAFitterNMl, "Process data with DCAFitterN and ML", false);

  void processDataWithKFParticleMl(CollisionsWCentMult const& collisions,
                                   SelectedCandidatesKfMl const&,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, false, false, false>(collisions, candidatesKfMlAll, tracks, bcs);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processDataWithKFParticleMl, "Process data with KFParticle and ML", false);

  void processMcWithDCAFitterMlSig(CollisionsWMcCentMult const& collisions,
                                   SelectedCandidatesMcMl const&,
                                   TypeMcCollisions const& mcCollisions,
                                   MatchedGenCandidatesMc const& mcParticles,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, true, false, true>(collisions, candidatesMcMlSig, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterMlSig, "Process MC with DCAFitterN and ML only for signals", false);

  void processMcWithDCAFitterMlBkg(CollisionsWMcCentMult const& collisions,
                                   SelectedCandidatesMcMl const&,
                                   TypeMcCollisions const& mcCollisions,
                                   MatchedGenCandidatesMc const& mcParticles,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, true, true, false>(collisions, candidatesMcMlBkg, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterMlBkg, "Process MC with DCAFitterN and ML only for background", false);

  void processMcWithDCAFitterMlAll(CollisionsWMcCentMult const& collisions,
                                   SelectedCandidatesMcMl const&,
                                   TypeMcCollisions const& mcCollisions,
                                   MatchedGenCandidatesMc const& mcParticles,
                                   aod::Tracks const& tracks,
                                   aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::DCAFitter, true, true, false, false>(collisions, candidatesMcMlAll, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithDCAFitterMlAll, "Process MC with DCAFitterN and ML", false);

  void processMcWithKFParticleMlSig(CollisionsWMcCentMult const& collisions,
                                    SelectedCandidatesMcKfMl const&,
                                    TypeMcCollisions const& mcCollisions,
                                    MatchedGenCandidatesMc const& mcParticles,
                                    aod::Tracks const& tracks,
                                    aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, true, false, true>(collisions, candidatesMcKfMlSig, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleMlSig, "Process MC with KFParticle and ML only for signals", false);

  void processMcWithKFParticleMlBkg(CollisionsWMcCentMult const& collisions,
                                    SelectedCandidatesMcKfMl const&,
                                    TypeMcCollisions const& mcCollisions,
                                    MatchedGenCandidatesMc const& mcParticles,
                                    aod::Tracks const& tracks,
                                    aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, true, true, false>(collisions, candidatesMcKfMlBkg, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleMlBkg, "Process MC with KFParticle and ML only for background", false);

  void processMcWithKFParticleMlAll(CollisionsWMcCentMult const& collisions,
                                    SelectedCandidatesMcKfMl const&,
                                    TypeMcCollisions const& mcCollisions,
                                    MatchedGenCandidatesMc const& mcParticles,
                                    aod::Tracks const& tracks,
                                    aod::BCs const& bcs)
  {
    rowsCommon.preProcessMcCollisions(mcCollisions, mcParticlesPerMcCollision, mcParticles);
    processCandidates<aod::hf_cand::VertexerType::KfParticle, true, true, false, false>(collisions, candidatesMcKfMlAll, tracks, bcs);
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcWithKFParticleMlAll, "Process MC with KFParticle and ML", false);

  void processMcGenOnly(TypeMcCollisions const& mcCollisions,
                        MatchedGenCandidatesMc const& mcParticles)
  {
    rowsCommon.processMcParticles(mcCollisions, mcParticlesPerMcCollision, mcParticles, Mass);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorD0ToKPi, processMcGenOnly, "Process MC gen. only", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorD0ToKPi>(cfgc)};
}
