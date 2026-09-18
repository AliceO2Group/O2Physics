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
#include "PWGHF/Core/HfMlResponseDplusToPiKPi.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/ReducedDMesonPairsTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
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

  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<bool> applySkimming{"applySkimming", false, "Flag to apply Skimming selections"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTDPlus"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_DPlusToKPiPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};

  Configurable<std::vector<double>> cutPtSkimming{"cutPtSkimming", {1, 5, 1000}, "pT bin limits for Skimming application"};
  Configurable<std::vector<double>> minM{"minM", {0.7, 0.7}, "Mass minimal for the cut for each pt bin"};
  Configurable<std::vector<double>> maxM{"maxM", {2.0, 2.1}, "Mass maximal for the cut for each pt bin"};
  Configurable<std::vector<double>> minCosTheta{"minCosTheta", {0.96, 0.98}, "CosTheta minimal for the cut for each pt bin"};
  Configurable<std::vector<double>> minDecayLength{"minDecayLength", {0.02, 0.03}, "DecayLength minimal for the cut for each pt bin"};
  Configurable<std::vector<double>> maxNsigmaTPC{"maxNsigmaTPC", {3, 3}, "NsigmaTPC maximal for the cut for each pt bin"};
  Configurable<std::vector<double>> maxNsigmaTOF{"maxNsigmaTOF", {3, 3}, "NsigmaTOF maximal for the cut for each pt bin"};

  Configurable<std::vector<double>> binsPtSkimming{"binsPtSkimming", {0}, "pT bin limits for Skimming application"};

  HfMlResponseDplusToPiKPi<float> hfMlResponse;

  std::vector<float> outputML;
  o2::ccdb::CcdbApi ccdbApi;

  HfHelper hfHelper;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKa, aod::HfSelDplusToPiKPi>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKa, aod::HfCand3ProngMcRec, aod::HfSelDplusToPiKPi>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;
  using SelectedCandidatesMcWithMl = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKa, aod::HfCand3ProngMcRec, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CollisionsCent = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::CentFT0Ms>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Filter filterMcGenMatching = (nabs(o2::aod::hf_cand_mc_flag::flagMcMatchGen) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi)) || (fillCorrBkgs && (nabs(o2::aod::hf_cand_mc_flag::flagMcMatchGen) != 0));

  Preslice<SelectedCandidates> tracksPerCollision = o2::aod::track::collisionId;
  Preslice<aod::McParticles> mcParticlesPerMcCollision = o2::aod::mcparticle::mcCollisionId;

  Partition<SelectedCandidatesMc> reconstructedCandSig = (nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi)) || (fillCorrBkgs && (nabs(o2::aod::hf_cand_mc_flag::flagMcMatchRec) != 0));
  Partition<SelectedCandidatesMc> reconstructedCandBkg = nabs(aod::hf_cand_mc_flag::flagMcMatchRec) != static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi);
  Partition<SelectedCandidatesMcWithMl> reconstructedCandSigMl = (nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi)) || (fillCorrBkgs && (nabs(o2::aod::hf_cand_mc_flag::flagMcMatchRec) != 0));

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

    if (applyMl) {
      hfMlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        hfMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        hfMlResponse.setModelPathsLocal(onnxFileNames);
      }
      hfMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      hfMlResponse.init();
    }
  }

  bool skimming(auto const& candidate,
                std::vector<double> const& PtcutSkimming,
                std::vector<double> const& Mmin,
                std::vector<double> const& Mmax,
                std::vector<double> const& CosThetamin,
                std::vector<double> const& DecayLengthmin,
                std::vector<double> const& NsigmaTPCmax,
                std::vector<double> const& NsigmaTOFmax)
  {
    if (candidate.pt() < PtcutSkimming[0] || candidate.pt() > PtcutSkimming[PtcutSkimming.size() - 1]) {
      return false;
    }
    for (size_t i = 1; i < PtcutSkimming.size(); i++) {
      if (candidate.pt() <= PtcutSkimming[i]) {
        if (hfHelper.invMassDplusToPiKPi(candidate) < Mmin[i - 1] ||
            hfHelper.invMassDplusToPiKPi(candidate) > Mmax[i - 1] ||
            candidate.cpa() < CosThetamin[i - 1] ||
            candidate.decayLength() < DecayLengthmin[i - 1] ||
            abs(candidate.nSigTofKa1()) > NsigmaTOFmax[i - 1] ||
            abs(candidate.nSigTpcKa1()) > NsigmaTPCmax[i - 1]) {
          return false;
        }
        return true;
      }
    }
    return false;
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

  template <typename Coll, bool DoMc = false, bool DoMl = false, typename T>
  void fillCandidateTable(const T& candidate, int localEvIdx = -1, int sign = 1)
  {
    int8_t flagMc = 0;
    int8_t originMc = 0;
    int8_t channelMc = 0;

    if constexpr (DoMc) {
      flagMc = candidate.flagMcMatchRec();
      originMc = candidate.originMcRec();
      channelMc = candidate.flagMcDecayChanRec();
    }

    std::vector<float> outML = {-999., -999.};
    if constexpr (DoMl) {
      for (unsigned int iclass = 0; iclass < classMlIndexes->size(); iclass++) {
        outML[iclass] = candidate.mlProbDplusToPiKPi()[classMlIndexes->at(iclass)];
      }
      rowCandidateMl(
        outML[0],
        outML[1]);
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

  void processData(aod::Collisions const& collisions,
                   SelectedCandidates const& candidates,
                   aod::Tracks const&,
                   aod::BCsWithTimestamps const&)
  {
    std::vector<double> skimmingCutPt = cutPtSkimming;
    std::vector<double> skimmingminM = minM;
    std::vector<double> skimmingmaxM = maxM;
    std::vector<double> skimmingminCosTheta = minCosTheta;
    std::vector<double> skimmingminDecayLength = minDecayLength;
    std::vector<double> skimmingmaxNsigmaTPC = maxNsigmaTPC;
    std::vector<double> skimmingmaxNsigmaTOF = maxNsigmaTOF;
    static int lastRunNumber = -1;
    // reserve memory
    rowCandidateFullEvents.reserve(collisions.size());
    if (fillCandidateTinyTable) {
      rowCandidateTiny.reserve(candidates.size());
    } else if (fillCandidateLiteTable) {
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
        auto prongCandidate = candidate.prong1_as<aod::Tracks>();
        auto candidateSign = -prongCandidate.sign();

        if (applySkimming &&
            !skimming(candidate,
                      skimmingCutPt,
                      skimmingminM,
                      skimmingmaxM,
                      skimmingminCosTheta,
                      skimmingminDecayLength,
                      skimmingmaxNsigmaTPC,
                      skimmingmaxNsigmaTOF)) {
          continue;
        }

        if (applyMl) {
          std::vector<float> inputFeatures = hfMlResponse.getInputFeatures(candidate);
          bool const isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, abs(candidate.pt()), outputML);
          if (!isSelectedMl) {
            continue;
          }
        }
        fillCandidateTable<aod::Collisions>(candidate, rowCandidateFullEvents.lastIndex(), candidateSign);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusDplusReduced, processData, "Process data per collision", false);

  void processMcRec(aod::Collisions const& collisions,
                    SelectedCandidatesMc const& candidates,
                    aod::Tracks const&)
  {
    std::vector<double> skimmingCutPt = cutPtSkimming;
    std::vector<double> skimmingminM = minM;
    std::vector<double> skimmingmaxM = maxM;
    std::vector<double> skimmingminCosTheta = minCosTheta;
    std::vector<double> skimmingminDecayLength = minDecayLength;
    std::vector<double> skimmingmaxNsigmaTPC = maxNsigmaTPC;
    std::vector<double> skimmingmaxNsigmaTOF = maxNsigmaTOF;
    // reserve memory
    rowCandidateFullEvents.reserve(collisions.size());
    if (fillCandidateTinyTable) {
      rowCandidateTiny.reserve(candidates.size());
    } else if (fillCandidateLiteTable) {
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
        auto prongCandidate = candidate.prong1_as<aod::Tracks>();
        auto candidateSign = -prongCandidate.sign();

        if (applySkimming &&
            !skimming(candidate,
                      skimmingCutPt,
                      skimmingminM,
                      skimmingmaxM,
                      skimmingminCosTheta,
                      skimmingminDecayLength,
                      skimmingmaxNsigmaTPC,
                      skimmingmaxNsigmaTOF)) {
          continue;
        }
        if (applyMl) {
          std::vector<float> inputFeatures = hfMlResponse.getInputFeatures(candidate);
          bool const isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, abs(candidate.pt()), outputML);
          if (!isSelectedMl) {
            continue;
          }
        }
        fillCandidateTable<aod::Collisions, true>(candidate, rowCandidateFullEvents.lastIndex(), candidateSign);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusDplusReduced, processMcRec, "Process data per collision", false);

  void processMcGen(aod::McCollisions const& mccollisions,
                    MatchedGenCandidatesMc const& mcparticles)
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
