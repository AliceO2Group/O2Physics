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

/// \file correlatorDplusDplusReduced.cxx
/// \brief Writer of D+ → π+ K- π+ candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug, local optimization of analysis on small samples or ML training.
///        In this file are defined and filled the output tables
///
/// \author Valerio DI BELLA <valerio.di.bella@cern.ch>, IPHC Strasbourg
/// Based on the code of Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/ReducedDMesonPairsTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;

/// Writes the full information in an output TTree
struct HfCorrelatorDplusDplusReduced {
  Produces<o2::aod::HfCandDpFulls> rowCandidateFull;
  Produces<o2::aod::HfCandDpLites> rowCandidateLite;
  Produces<o2::aod::HfCandDpTinys> rowCandidateTiny;
  Produces<o2::aod::HfCandDpFullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandDpMls> rowCandidateMl;

  Produces<o2::aod::HfCandDpMcPs> rowCandidateMcParticles;
  Produces<o2::aod::HfCandDpMcEvs> rowCandidateMcCollisions;

  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 1, "Selection Flag for Dplus"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  Configurable<bool> fillCandidateTinyTable{"fillCandidateTinyTable", false, "Switch to fill tiny table with candidate properties"};
  // parameters for production of training samples
  Configurable<bool> fillCorrBkgs{"fillCorrBkgs", false, "Flag to fill derived tables with correlated background candidates"};
  Configurable<std::vector<int>> classMlIndexes{"classMlIndexes", {0, 2}, "Indexes of ML bkg and non-prompt scores."};
  Configurable<int> centEstimator{"centEstimator", 0, "Centrality estimation (None: 0, FT0C: 2, FT0M: 3)"};
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", true, "Enables processing of skimmed datasets"};
  Configurable<bool> skipSingleD{"skipSingleD", true, "Skip collisions with one or less D candidates"};

  HfHelper hfHelper;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKa, aod::HfSelDplusToPiKPi>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKa, aod::HfCand3ProngMcRec, aod::HfSelDplusToPiKPi>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;
  using SelectedCandidatesMcWithMl = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKa, aod::HfCand3ProngMcRec, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CollisionsCent = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::CentFT0Ms>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Filter filterMcGenMatching = (nabs(o2::aod::hf_cand_3prong::flagMcMatchGen) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi)) || (fillCorrBkgs && (nabs(o2::aod::hf_cand_3prong::flagMcMatchGen) != 0));

  Preslice<SelectedCandidates> tracksPerCollision = o2::aod::track::collisionId;
  Preslice<aod::McParticles> mcParticlesPerMcCollision = o2::aod::mcparticle::mcCollisionId;

  Partition<SelectedCandidatesMc> reconstructedCandSig = (nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi)) || (fillCorrBkgs && (nabs(o2::aod::hf_cand_3prong::flagMcMatchRec) != 0));
  Partition<SelectedCandidatesMc> reconstructedCandBkg = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi);
  Partition<SelectedCandidatesMcWithMl> reconstructedCandSigMl = (nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi)) || (fillCorrBkgs && (nabs(o2::aod::hf_cand_3prong::flagMcMatchRec) != 0));

  HistogramRegistry registry{"registry"};
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  void init(InitContext const&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    if (cfgSkimmedProcessing) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }
  }

  template <typename T>
  void fillEvent(const T& collision)
  {
    rowCandidateFullEvents(
      collision.numContrib(),
      collision.posX(),
      collision.posY(),
      collision.posZ());
  }

  template <typename Coll, bool doMc = false, bool doMl = false, typename T>
  void fillCandidateTable(const T& candidate, int localEvIdx = -1, int sign = 1)
  {
    int8_t flagMc = 0;
    int8_t originMc = 0;
    int8_t channelMc = 0;

    if constexpr (doMc) {
      flagMc = candidate.flagMcMatchRec();
      originMc = candidate.originMcRec();
      channelMc = candidate.flagMcDecayChanRec();
    }

    std::vector<float> outputMl = {-999., -999.};
    if constexpr (doMl) {
      for (unsigned int iclass = 0; iclass < classMlIndexes->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMlIndexes->at(iclass)];
      }
      rowCandidateMl(
        outputMl[0],
        outputMl[1]);
    }

    float cent{-1.};
    auto coll = candidate.template collision_as<Coll>();
    if (std::is_same_v<Coll, CollisionsCent> && centEstimator != CentralityEstimator::None) {
      cent = getCentralityColl(coll, centEstimator);
    }

    if (fillCandidateTinyTable) {
      rowCandidateTiny(
        candidate.isSelDplusToPiKPi(),
        hfHelper.invMassDplusToPiKPi(candidate),
        sign * candidate.pt(),
        candidate.eta(),
        candidate.phi(),
        localEvIdx,
        flagMc,
        originMc,
        channelMc);
    } else if (fillCandidateLiteTable) {
      rowCandidateLite(
        candidate.chi2PCA(),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXYNormalised(),
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.ptProng2(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impactParameter2(),
        candidate.impactParameterZ0(),
        candidate.impactParameterZ1(),
        candidate.impactParameterZ2(),
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
        candidate.nSigTpcPi2(),
        candidate.nSigTpcKa2(),
        candidate.nSigTofPi2(),
        candidate.nSigTofKa2(),
        candidate.tpcTofNSigmaPi2(),
        candidate.tpcTofNSigmaKa2(),
        candidate.isSelDplusToPiKPi(),
        hfHelper.invMassDplusToPiKPi(candidate),
        sign * candidate.pt(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        candidate.eta(),
        candidate.phi(),
        hfHelper.yDplus(candidate),
        cent,
        localEvIdx,
        flagMc,
        originMc,
        channelMc);
    } else {
      rowCandidateFull(
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
        candidate.impactParameterNormalised2(),
        candidate.ptProng2(),
        RecoDecay::p(candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()),
        candidate.pxProng0(),
        candidate.pyProng0(),
        candidate.pzProng0(),
        candidate.pxProng1(),
        candidate.pyProng1(),
        candidate.pzProng1(),
        candidate.pxProng2(),
        candidate.pyProng2(),
        candidate.pzProng2(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impactParameter2(),
        candidate.errorImpactParameter0(),
        candidate.errorImpactParameter1(),
        candidate.errorImpactParameter2(),
        candidate.impactParameterZ0(),
        candidate.impactParameterZ1(),
        candidate.impactParameterZ2(),
        candidate.errorImpactParameterZ0(),
        candidate.errorImpactParameterZ1(),
        candidate.errorImpactParameterZ2(),
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
        candidate.nSigTpcPi2(),
        candidate.nSigTpcKa2(),
        candidate.nSigTofPi2(),
        candidate.nSigTofKa2(),
        candidate.tpcTofNSigmaPi2(),
        candidate.tpcTofNSigmaKa2(),
        candidate.isSelDplusToPiKPi(),
        hfHelper.invMassDplusToPiKPi(candidate),
        sign * candidate.pt(),
        candidate.p(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        hfHelper.ctDplus(candidate),
        candidate.eta(),
        candidate.phi(),
        hfHelper.yDplus(candidate),
        hfHelper.eDplus(candidate),
        cent,
        localEvIdx,
        flagMc,
        originMc,
        channelMc);
    }
  }

  void processData(aod::Collisions const& collisions, SelectedCandidates const& candidates, aod::Tracks const&)
  {
    static int lastRunNumber = -1;
    // reserve memory
    rowCandidateFullEvents.reserve(collisions.size());
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidates.size());
    } else {
      rowCandidateFull.reserve(candidates.size());
    }

    for (const auto& collision : collisions) {
      if (cfgSkimmedProcessing) {
        auto bc = collision.bc_as<aod::BCsWithTimestamps>();
        int runNumber = bc.runNumber();
        if (lastRunNumber != runNumber) {
          lastRunNumber = runNumber;
          LOGF(info, "Initializing Zorro for run %d", runNumber);
          uint64_t currentTimestamp = bc.timestamp();
          zorro.initCCDB(ccdb.service, runNumber, currentTimestamp, "fHfDoubleCharm3P");
          zorro.populateHistRegistry(registry, runNumber);
        }
        zorro.isSelected(bc.globalBC());
      }

      const auto colId = collision.globalIndex();
      auto candidatesInThisCollision = candidates.sliceBy(tracksPerCollision, colId);
      if (skipSingleD)
        if (candidatesInThisCollision.size() < 2) // o2-linter: disable=magic-number (number of candidate must be larger than 1)
          continue;
      fillEvent(collision);
      for (const auto& candidate : candidatesInThisCollision) {
        auto prong_candidate = candidate.prong1_as<aod::Tracks>();
        auto candidate_sign = prong_candidate.sign();
        fillCandidateTable<aod::Collisions>(candidate, rowCandidateFullEvents.lastIndex(), candidate_sign);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusDplusReduced, processData, "Process data per collision", false);

  void processMcRec(aod::Collisions const& collisions,
                    SelectedCandidatesMc const& candidates)
  {
    // reserve memory
    rowCandidateFullEvents.reserve(collisions.size());
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidates.size());
    } else {
      rowCandidateFull.reserve(candidates.size());
    }

    for (const auto& collision : collisions) { // No skimming for MC data. No Zorro !
      const auto colId = collision.globalIndex();
      auto candidatesInThisCollision = candidates.sliceBy(tracksPerCollision, colId);
      if (skipSingleD)
        if (candidatesInThisCollision.size() < 2) // o2-linter: disable=magic-number (number of candidate must be larger than 1)
          continue;
      fillEvent(collision);
      for (const auto& candidate : candidatesInThisCollision) {
        auto prong_candidate = candidate.prong1_as<aod::Tracks>();
        auto candidate_sign = prong_candidate.sign();
        fillCandidateTable<aod::Collisions, true>(candidate, rowCandidateFullEvents.lastIndex(), candidate_sign);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusDplusReduced, processMcRec, "Process data per collision", false);

  void processMcGen(aod::McCollisions const& mccollisions, MatchedGenCandidatesMc const& mcparticles)
  {
    // reserve memory
    rowCandidateMcCollisions.reserve(mccollisions.size());
    rowCandidateMcParticles.reserve(mcparticles.size());

    for (const auto& mccollision : mccollisions) { // No skimming for MC data. No Zorro !
      const auto colId = mccollision.globalIndex();
      const auto particlesInThisCollision = mcparticles.sliceBy(mcParticlesPerMcCollision, colId);
      if (skipSingleD)
        if (particlesInThisCollision.size() < 2) // o2-linter: disable=magic-number (number of candidate must be larger than 1)
          continue;
      rowCandidateMcCollisions(
        mccollision.posX(),
        mccollision.posY(),
        mccollision.posZ());
      for (const auto& particle : particlesInThisCollision) {
        rowCandidateMcParticles(
          particle.pt(),
          particle.eta(),
          particle.phi(),
          particle.y(),
          rowCandidateMcCollisions.lastIndex(),
          particle.flagMcMatchGen(),
          particle.flagMcDecayChanGen(),
          particle.originMcGen());
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusDplusReduced, processMcGen, "Process MC data at the generator level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDplusDplusReduced>(cfgc)};
}
