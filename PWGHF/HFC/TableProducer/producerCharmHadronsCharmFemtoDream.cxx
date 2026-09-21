// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file producerCharmHadronsCharmFemtoDream.cxx
/// \brief Produce reduced D-meson candidate tables for charm-charm femtoscopy
/// \author Biao Zhang, Heidelberg University, biao.zhang@cern.ch

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfMlResponseD0ToKPi.h"
#include "PWGHF/Core/HfMlResponseDstarToD0Pi.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/ML/MlResponse.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfProducerCharmHadronsCharmFemtoDream {
  enum MlMode : uint8_t {
    NoMl = 0,
    FillMlFromSelector = 1,
    FillMlFromNewBDT = 2
  };
  // Each species needs its own model, feature order and pT-dependent cuts.
  struct MlConfig : ConfigurableGroup {
    std::string prefix;
    static inline const std::array<double, 3> DefaultCuts{1., 0., 0.};
    Configurable<int> applyMlMode{"applyMlMode", FillMlFromSelector, "0: no ML, 1: selector scores, 2: new BDT after selector"};
    Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{0., 36.}, "pT bin limits for new BDT"};
    Configurable<LabeledArray<double>> cutsMl{"cutsMl", {DefaultCuts.data(), 1, 3}, "New BDT cuts per pT bin: background, prompt, nonprompt"};
    Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{0, 1, 1}, "Reject scores above (0), below (1), or do not cut (2)"};
    Configurable<int> nClassesMl{"nClassesMl", 3, "Three output classes: background, prompt, nonprompt"};
    Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{}, "Ordered input feature names for new BDT"};
    Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{}, "Model files, one per pT bin"};
    Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{}, "CCDB model paths, one per pT bin"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "Timestamp used to retrieve models"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Load new BDT from CCDB instead of local files"};
  } mlD0{{}, "mlD0"}, mlDstar{{}, "mlDstar"};
  o2::analysis::HfMlResponseD0ToKPi<float> hfMlResponseD0;
  o2::analysis::HfMlResponseDstarToD0Pi<float> hfMlResponseDstar;
  o2::ccdb::CcdbApi ccdbApi;

  Produces<aod::FDCollisions> collisions;
  Produces<aod::FDColMasks> rowMasks;
  bool hasD0 = false;
  bool hasDstar = false;
  Produces<aod::FDHfCand2Prong> d0Rows;
  Produces<aod::FDHfCandDstar> dstarRows;
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Minimum D0 selector decision"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "CCDB URL"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "Run 3 magnetic field"};
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::hf_evsel::HfEventSelection hfEvSel;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  HistogramRegistry registry{"QA", {}, OutputObjHandlingPolicy::AnalysisObject};
  int runNumber = -1;
  using Collisions =
    soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using CollisionsWithFT0C = soa::Join<Collisions, aod::CentFT0Cs>;
  using D0s = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;
  using Dstars =
    soa::Join<aod::HfCandDstars, aod::HfD0FromDstar, aod::HfSelDstarToD0Pi>;

  using D0sMl = soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0, aod::HfMlD0>;
  using DstarsMl = soa::Join<aod::HfCandDstarsWPid, aod::HfD0FromDstar, aod::HfSelDstarToD0Pi, aod::HfMlDstarToD0Pi>;

  Filter filterSelectCandidateD0 = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0);
  Filter filterSelectCandidateDstar = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == true;

  template <typename Response>
  void initMl(MlConfig const& cfg, Response& response, bool withMl)
  {
    if (cfg.applyMlMode.value < NoMl || cfg.applyMlMode.value > FillMlFromNewBDT) {
      LOGP(fatal, "{}: invalid applyMlMode", cfg.prefix);
    }
    if (cfg.applyMlMode.value != FillMlFromNewBDT) {
      return;
    }
    if (!withMl) {
      LOGP(fatal, "{}: new BDT requires an Ml producer process", cfg.prefix);
    }
    auto const& edges = cfg.binsPtMl.value;
    if (edges.size() < 2 ||
        !std::all_of(edges.begin(), edges.end(), [](double x) { return std::isfinite(x); }) ||
        std::adjacent_find(edges.begin(), edges.end(), [](double a, double b) { return a >= b; }) != edges.end()) {
      LOGP(fatal, "{}: binsPtMl must be finite and strictly increasing", cfg.prefix);
    }
    const auto nBins = edges.size() - 1;
    if (cfg.nClassesMl.value != 3 || cfg.cutDirMl.value.size() != 3 ||
        cfg.cutsMl.value.rows() != nBins || cfg.cutsMl.value.cols() != 3 ||
        cfg.onnxFileNames.value.size() != nBins || cfg.namesInputFeatures.value.empty() ||
        (cfg.loadModelsFromCCDB.value && cfg.modelPathsCCDB.value.size() != nBins)) {
      LOGP(fatal, "{}: provide three classes, cuts/models for every pT bin, and input features", cfg.prefix);
    }
    for (auto direction : cfg.cutDirMl.value) {
      if (direction < o2::cuts_ml::CutGreater || direction > o2::cuts_ml::CutNot) {
        LOGP(fatal, "{}: invalid cutDirMl", cfg.prefix);
      }
    }
    for (unsigned int bin = 0; bin < nBins; ++bin) {
      for (unsigned int score = 0; score < 3; ++score) {
        if (!std::isfinite(cfg.cutsMl.value.get(bin, score))) {
          LOGP(fatal, "{}: cutsMl must be finite", cfg.prefix);
        }
      }
    }
    response.configure(edges, cfg.cutsMl.value, cfg.cutDirMl.value, cfg.nClassesMl.value);
    response.cacheInputFeaturesIndices(cfg.namesInputFeatures.value);
    if (cfg.loadModelsFromCCDB.value) {
      ccdbApi.init(ccdbUrl.value);
      response.setModelPathsCCDB(cfg.onnxFileNames.value, ccdbApi, cfg.modelPathsCCDB.value, cfg.timestampCCDB.value);
    } else {
      response.setModelPathsLocal(cfg.onnxFileNames.value);
    }
    response.init();
  }

  void init(InitContext const&)
  {
    if (static_cast<int>(doprocessD0D0) + doprocessD0Dstar + doprocessD0D0Ml +
          doprocessD0DstarMl + doprocessD0D0WithFT0C + doprocessD0DstarWithFT0C +
          doprocessD0D0MlWithFT0C + doprocessD0DstarMlWithFT0C !=
        1) {
      LOGP(fatal, "Enable exactly one charm-charm producer process");
    }
    if (selectionFlagD0 < 1) {
      LOGP(fatal, "selectionFlagD0 must be positive");
    }
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    const bool withMl = doprocessD0D0Ml || doprocessD0DstarMl || doprocessD0D0MlWithFT0C || doprocessD0DstarMlWithFT0C;
    initMl(mlD0, hfMlResponseD0, withMl);
    if (doprocessD0Dstar || doprocessD0DstarMl || doprocessD0DstarWithFT0C || doprocessD0DstarMlWithFT0C) {
      initMl(mlDstar, hfMlResponseDstar, withMl);
    }
    hfEvSel.init(registry, &zorroSummary);
    registry.add("events", ";stage (0=all,1=accepted);events", kTH1F,
                 {{2, -0.5, 1.5}});
    registry.add("candidates", ";species (0=D0,1=D0bar,2=Dstar+,3=Dstar-);rows",
                 kTH1F, {{4, -0.5, 3.5}});
  }

  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename Collision>
  bool acceptCollision(Collision const& col)
  {
    registry.fill(HIST("events"), 0);
    float cent = -1.f; // No centrality for pp MB.
    const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<
      true, CentEstimator,
      aod::BCsWithTimestamps>(col, cent, ccdb, registry);
    hfEvSel.fillHistograms(col, rejectionMask, cent);
    if (rejectionMask != 0) {
      return false;
    }
    hasD0 = false;
    hasDstar = false;
    auto bc = col.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag.value, nullptr, false);
    // Propagator field is in kG; FemtoDream uses Tesla.
    const float fieldTesla =
      0.1f * o2::base::Propagator::Instance()->getNominalBz();
    // MultV0M is the common reduced percentile column: FT0C or -1 (no estimator).
    collisions(col.posZ(), cent, col.multNTracksPV(), 2.f,
               fieldTesla);
    registry.fill(HIST("events"), 1);
    // Keep all accepted collisions, including zero/single-candidate events.
    return true;
  }

  template <typename Scores>
  std::array<float, 3> readScores(Scores const& scores)
  {
    if (scores.size() != 3) {
      LOGP(fatal,
           "Expected three selector ML scores: background, prompt, nonprompt");
      return {-1.f, -1.f, -1.f};
    }
    for (auto score : scores) {
      if (!std::isfinite(score)) {
        LOGP(fatal, "Non-finite selector ML score");
      }
    }
    return {scores[0], scores[1], scores[2]};
  }

  template <bool WithMl, typename Collision, typename Candidates>
  void fillD0(Collision const& col, Candidates const& candidates)
  {
    const auto timestamp = col.template bc_as<aod::BCsWithTimestamps>().timestamp();
    for (auto const& cand : candidates) {
      auto p0 = cand.template prong0_as<aod::Tracks>();
      auto p1 = cand.template prong1_as<aod::Tracks>();
      // The OR filter accepts the row if either hypothesis passes. Write only
      // the passing hypotheses, each with its own flavour and ML scores.
      for (int hypothesis = 0; hypothesis < 2; ++hypothesis) {
        if ((hypothesis == 0 ? cand.isSelD0() : cand.isSelD0bar()) <
            selectionFlagD0) {
          continue;
        }
        std::array<float, 3> scores{-1.f, -1.f, -1.f};
        if constexpr (WithMl) {
          if (mlD0.applyMlMode != NoMl) {
            scores = hypothesis == 0 ? readScores(cand.mlProbD0())
                                     : readScores(cand.mlProbD0bar());
          }
          if (mlD0.applyMlMode == FillMlFromNewBDT) {
            // Do not call the ML response with an out-of-range model index.
            if (!std::isfinite(cand.pt()) || cand.pt() < mlD0.binsPtMl.value.front() || cand.pt() >= mlD0.binsPtMl.value.back()) {
              continue;
            }
            const int pdgCode = hypothesis == 0 ? o2::constants::physics::kD0 : -o2::constants::physics::kD0;
            auto features = hfMlResponseD0.getInputFeatures<true>(cand, pdgCode);
            std::vector<float> output;
            if (!hfMlResponseD0.isSelectedMl(features, cand.pt(), output)) {
              continue; // Reject this hypothesis only.
            }
            scores = readScores(output);
          }
        }
        d0Rows(collisions.lastIndex(), timestamp, hypothesis == 0 ? 1 : -1,
               p0.globalIndex(), p1.globalIndex(), p0.pt(), p1.pt(), p0.eta(),
               p1.eta(), p0.phi(), p1.phi(), 1 << hypothesis, scores[0],
               scores[1], scores[2]);
        hasD0 = true;
        registry.fill(HIST("candidates"), hypothesis);
      }
    }
  }

  template <bool WithMl, typename Collision, typename Candidates>
  void fillDstar(Collision const& col,
                 Candidates const& candidates)
  {
    const auto timestamp = col.template bc_as<aod::BCsWithTimestamps>().timestamp();
    for (auto const& cand : candidates) {
      auto p0 = cand.template prong0_as<aod::Tracks>();
      auto p1 = cand.template prong1_as<aod::Tracks>();
      auto soft = cand.template prongPi_as<aod::Tracks>();
      std::array<float, 3> scores{-1.f, -1.f, -1.f};
      if constexpr (WithMl) {
        if (mlDstar.applyMlMode != NoMl) {
          scores = readScores(cand.mlProbDstarToD0Pi());
        }
        if (mlDstar.applyMlMode == FillMlFromNewBDT) {
          if (!std::isfinite(cand.pt()) || cand.pt() < mlDstar.binsPtMl.value.front() || cand.pt() >= mlDstar.binsPtMl.value.back()) {
            continue;
          }
          // Match the unswapped daughter convention used by D+track.
          auto features = hfMlResponseDstar.getInputFeatures(cand, false);
          std::vector<float> output;
          if (!hfMlResponseDstar.isSelectedMl(features, cand.pt(), output)) {
            continue;
          }
          scores = readScores(output);
        }
      }
      dstarRows(collisions.lastIndex(), timestamp, soft.sign(),
                p0.globalIndex(), p1.globalIndex(), soft.globalIndex(), p0.pt(),
                p1.pt(), soft.pt(), p0.eta(), p1.eta(), soft.eta(), p0.phi(),
                p1.phi(), soft.phi(), 1, scores[0], scores[1], scores[2]);
      hasDstar = true;
      registry.fill(HIST("candidates"), soft.sign() > 0 ? 2 : 3);
    }
  }

  void processD0D0(Collisions::iterator const& col,
                   aod::BCsWithTimestamps const&, aod::Tracks const&,
                   soa::Filtered<D0s> const& d0s)
  {
    if (acceptCollision<o2::hf_centrality::CentralityEstimator::None>(col)) {
      fillD0<false>(col, d0s);
      rowMasks(hasD0 ? 1 : 0, hasDstar ? 1 : 0, 0);
    }
  }
  PROCESS_SWITCH(HfProducerCharmHadronsCharmFemtoDream, processD0D0, "D0 only, data", true);

  void processD0Dstar(Collisions::iterator const& col,
                      aod::BCsWithTimestamps const&, aod::Tracks const&,
                      soa::Filtered<D0s> const& d0s, soa::Filtered<Dstars> const& dstars)
  {
    if (acceptCollision<o2::hf_centrality::CentralityEstimator::None>(col)) {
      fillD0<false>(col, d0s);
      fillDstar<false>(col, dstars);
      rowMasks(hasD0 ? 1 : 0, hasDstar ? 1 : 0, 0);
    }
  }
  PROCESS_SWITCH(HfProducerCharmHadronsCharmFemtoDream, processD0Dstar, "D0 and Dstar, data", false);

  void processD0D0Ml(Collisions::iterator const& col,
                     aod::BCsWithTimestamps const&, aod::Tracks const&,
                     soa::Filtered<D0sMl> const& d0s)
  {
    if (acceptCollision<o2::hf_centrality::CentralityEstimator::None>(col)) {
      fillD0<true>(col, d0s);
      rowMasks(hasD0 ? 1 : 0, hasDstar ? 1 : 0, 0);
    }
  }
  PROCESS_SWITCH(HfProducerCharmHadronsCharmFemtoDream, processD0D0Ml, "D0 with selector ML scores", false);

  void processD0DstarMl(Collisions::iterator const& col,
                        aod::BCsWithTimestamps const&, aod::Tracks const&,
                        soa::Filtered<D0sMl> const& d0s,
                        soa::Filtered<DstarsMl> const& dstars)
  {
    if (acceptCollision<o2::hf_centrality::CentralityEstimator::None>(col)) {
      fillD0<true>(col, d0s);
      fillDstar<true>(col, dstars);
      rowMasks(hasD0 ? 1 : 0, hasDstar ? 1 : 0, 0);
    }
  }
  PROCESS_SWITCH(HfProducerCharmHadronsCharmFemtoDream, processD0DstarMl, "D0 and Dstar with selector ML scores", false);

  void processD0D0WithFT0C(CollisionsWithFT0C::iterator const& col,
                           aod::BCsWithTimestamps const&, aod::Tracks const&,
                           soa::Filtered<D0s> const& d0s)
  {
    if (acceptCollision<o2::hf_centrality::CentralityEstimator::FT0C>(col)) {
      fillD0<false>(col, d0s);
      rowMasks(hasD0 ? 1 : 0, hasDstar ? 1 : 0, 0);
    }
  }
  PROCESS_SWITCH(HfProducerCharmHadronsCharmFemtoDream, processD0D0WithFT0C, "D0 only, data with FT0C centrality", false);

  void processD0DstarWithFT0C(CollisionsWithFT0C::iterator const& col,
                              aod::BCsWithTimestamps const&, aod::Tracks const&,
                              soa::Filtered<D0s> const& d0s, soa::Filtered<Dstars> const& dstars)
  {
    if (acceptCollision<o2::hf_centrality::CentralityEstimator::FT0C>(col)) {
      fillD0<false>(col, d0s);
      fillDstar<false>(col, dstars);
      rowMasks(hasD0 ? 1 : 0, hasDstar ? 1 : 0, 0);
    }
  }
  PROCESS_SWITCH(HfProducerCharmHadronsCharmFemtoDream, processD0DstarWithFT0C, "D0 and Dstar, data with FT0C centrality", false);

  void processD0D0MlWithFT0C(CollisionsWithFT0C::iterator const& col,
                             aod::BCsWithTimestamps const&, aod::Tracks const&,
                             soa::Filtered<D0sMl> const& d0s)
  {
    if (acceptCollision<o2::hf_centrality::CentralityEstimator::FT0C>(col)) {
      fillD0<true>(col, d0s);
      rowMasks(hasD0 ? 1 : 0, hasDstar ? 1 : 0, 0);
    }
  }
  PROCESS_SWITCH(HfProducerCharmHadronsCharmFemtoDream, processD0D0MlWithFT0C, "D0 with selector ML scores and FT0C centrality", false);

  void processD0DstarMlWithFT0C(CollisionsWithFT0C::iterator const& col,
                                aod::BCsWithTimestamps const&, aod::Tracks const&,
                                soa::Filtered<D0sMl> const& d0s,
                                soa::Filtered<DstarsMl> const& dstars)
  {
    if (acceptCollision<o2::hf_centrality::CentralityEstimator::FT0C>(col)) {
      fillD0<true>(col, d0s);
      fillDstar<true>(col, dstars);
      rowMasks(hasD0 ? 1 : 0, hasDstar ? 1 : 0, 0);
    }
  }
  PROCESS_SWITCH(HfProducerCharmHadronsCharmFemtoDream, processD0DstarMlWithFT0C, "D0 and Dstar with selector ML scores and FT0C centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfProducerCharmHadronsCharmFemtoDream>(cfgc)};
}
