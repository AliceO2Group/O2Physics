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

/// \file taskPtFlucCharmHadrons.cxx
/// \brief Analysis task for charm hadron pt fluctuation
///
/// \author Prottay Das, prottay.das@cern.ch
/// \author Wu Chuntai, UNIPD, CCNU, and INFN Padova
/// \author Andrea Rossi, INFN Padova

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TF1.h>
#include <TH2.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;

enum DecayChannel { DplusToPiKPi = 0,
                    D0ToPiK,
                    D0ToKPi };

enum CandD0Type { PureD0 = 0,
                  PureD0bar,
                  ReflectedD0,
                  ReflectedD0bar
};

namespace o2::aod
{
namespace CandInfo
{
DECLARE_SOA_COLUMN(Pt, pt, float);                   //! D candidate pT (GeV/c)
DECLARE_SOA_COLUMN(Mass, mass, float);               //! D invariant mass (GeV/c^2)
DECLARE_SOA_COLUMN(Eta, eta, float);                 //! D candidate eta
DECLARE_SOA_COLUMN(BdtScoreBkg, bdtScoreBkg, float); //! BDT background score = ML prob @ classMl[0]
DECLARE_SOA_COLUMN(BdtScoreNp, bdtScoreNp, float);   //! BDT non-prompt score = ML prob @ classMl[1]
} // namespace CandInfo

DECLARE_SOA_TABLE(HfCandInfo, "AOD", "HFCANDINFO",
                  CandInfo::Pt,
                  CandInfo::Mass,
                  CandInfo::Eta,
                  CandInfo::BdtScoreBkg,
                  CandInfo::BdtScoreNp);
} // namespace o2::aod

struct HfTaskPtFlucCharmHadrons {
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4)"};
  Configurable<int> selectionFlag{"selectionFlag", 1, "Selection Flag for hadron (e.g. 1 for skimming, 3 for topo. and kine., 7 for PID)"};
  Configurable<float> centralityMin{"centralityMin", 0.f, "Minimum centrality accepted (additional local cut)"};
  Configurable<float> centralityMax{"centralityMax", 100.f, "Maximum centrality accepted (additional local cut)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indices of BDT scores to be stored. Two indexes max."};
  Configurable<bool> saveCharmBulkCorrelations{"saveCharmBulkCorrelations", false, "Flag to save charm-bulk correlations in THnSparse"};

  Configurable<bool> usePtDepDcaXY{"usePtDepDcaXY", true, "Use pT-dependent DCAxy cut"};
  Configurable<std::string> dcaXYFunc{"dcaXYFunc", "(0.0105+0.035/(x^1.1))", "pT-dependent DCAxy cut function, x = pT (GeV/c); globalTrack default"};
  Configurable<float> maxDcaZ{"maxDcaZ", 2.0f, "Maximum |DCAz| (cm)"};

  // MC track-efficiency map (charged-hadron reconstruction efficiency, one TH2F per selection).
  Configurable<bool> applyTrackEffCorrection{"applyTrackEffCorrection", false, "Enable 1/epsilon(pT,eta) track weighting for <pT> correction"};
  Configurable<std::string> ccdbPathEffBase{"ccdbPathEffBase", "Users/c/chuntai/TrackEff/", "CCDB path to the TH2F track-efficiency map"};
  Configurable<std::string> trackEffName{"trackEffName", "LHC25a5", "Anchored MC appended to the CCDB efficiency path"};
  Configurable<int64_t> ccdbNoLaterThanEff{"ccdbNoLaterThanEff", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable creation timestamp for the efficiency object"};
  Configurable<float> effMin{"effMin", 0.0001f, "Minimum accepted efficiency; bins below are treated as no-correction"};

  // candidate info tree
  Configurable<bool> fillCandInfoTree{"fillCandInfoTree", false, "Fill the per-candidate info table (pT, mass, eta, BDT bkg, BDT non-prompt)"};

  // Event selection
  // Mean-pT (charged hadrons) definition: subevents with eta gap
  Configurable<float> etaTrkAMin{"etaTrkAMin", -0.8f, "A: negative track eta min"};
  Configurable<float> etaTrkAMax{"etaTrkAMax", -0.2f, "A: negative track eta max"};
  Configurable<float> etaTrkBMin{"etaTrkBMin", 0.2f, "B: positive track eta min"};
  Configurable<float> etaTrkBMax{"etaTrkBMax", 0.8f, "B: positive track eta max"};
  Configurable<float> etaTrkMax{"etaTrkMax", 0.8f, "Track eta max"};
  Configurable<float> etaCandMax{"etaCandMax", 0.8f, "Candidate eta max"};

  Configurable<float> ptTrkMin{"ptTrkMin", 0.2f, "Track pT min for <pT> (charged hadrons)"};
  Configurable<float> ptTrkMax{"ptTrkMax", 5.0f, "Track pT max for <pT> (charged hadrons)"};
  Configurable<float> ptCandMin{"ptCandMin", 0.f, "Min pT for D candidates"};
  Configurable<float> ptCandMax{"ptCandMax", 12.0f, "Max pT for D candidates"};
  Configurable<int> minNTrk{"minNTrk", 5, "Min charged tracks in each subevent to compute <pT>"};

  // Use ML in denominator AND numerator consistently
  Configurable<bool> applyMlCut{"applyMlCut", true, "Apply ML window for candidate selection (also in ND denominator)"};

  // ML window (assumes you store 2 ML outputs: classMl[0], classMl[1])
  Configurable<float> mlOneMin{"mlOneMin", 0.f, "ML1 min"};
  Configurable<float> mlOneMax{"mlOneMax", 1.f, "ML1 max"};
  Configurable<float> mlTwoMin{"mlTwoMin", 0.f, "ML2 min"};
  Configurable<float> mlTwoMax{"mlTwoMax", 1.f, "ML2 max"};

  SliceCache cache;
  HfEventSelection hfEvSel; // event selection and monitoring
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  Produces<o2::aod::HfCandInfo> rowCandInfo;

  // Cached track-efficiency map (loaded once per run).
  TH2F* hTrackEffPtEta = nullptr;
  TF1* fPtDepDcaXY = nullptr;
  int64_t mLastRunNumber = -1;
  std::string mCcdbPathEff;

  using CandDplusDataWMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandD0DataWMl = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

  // MC tables for the track-efficiency
  using McRecCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::McCollisionLabels>;
  using McRecTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
  Preslice<McRecTracks> perMcRecCollision = aod::track::collisionId;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;
  Filter filterSelectD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Filter filterSelectTracks = (nabs(aod::track::eta) < etaTrkMax) && (aod::track::pt > ptTrkMin) && (aod::track::pt < ptTrkMax);
  Partition<CandD0DataWMl> selectedD0ToPiKWMl = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
  Partition<CandD0DataWMl> selectedD0ToKPiWMl = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;

  // -------------------------
  // Axes / histograms
  // -------------------------
  ConfigurableAxis axisInvMass{"axisInvMass", {500, 1.65, 2.15}, "Inv mass axis"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 12.0}, "Candidate pT axis"};
  ConfigurableAxis axisCent{"axisCent", {VARIABLE_WIDTH, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 80.0}, "Centrality axis"};
  ConfigurableAxis axisSign{"axisSign", {5, -1.5, 3.5}, "Sign axis (integer centers -1, 0, 1, 2, 3)"};
  ConfigurableAxis axisMlOne{"axisMlOne", {1000, 0., 1.}, "ML score 1 axis"};
  ConfigurableAxis axisMlTwo{"axisMlTwo", {100, 0., 1.}, "ML score 2 axis"};
  ConfigurableAxis axisCandEta{"axisCandEta", {16, -0.8, 0.8}, "Candidate eta axis"};
  ConfigurableAxis axisMPtTrkA{"axisMPtTrkA", {200, 0., 2.0}, "Mean pT of tracks in subevent A"};
  ConfigurableAxis axisMPtTrkB{"axisMPtTrkB", {200, 0., 2.0}, "Mean pT of tracks in subevent B"};
  // ConfigurableAxis axisPtCandProduct{"axisPtCandProduct", {150, 0., 15}, "Candidate pT * mean pT of tracks (A and B) product"}; // to be moved to the TTree
  ConfigurableAxis axisPtTrkProduct{"axisPtTrkProduct", {200, 0., 2.}, "Mean pT of tracks in A * Mean pT of tracks in B"};
  ConfigurableAxis axisNTrkA{"axisNTrkA", {2000, 0., 2000.}, "N_{tracks} in subevent A"};
  ConfigurableAxis axisNTrkB{"axisNTrkB", {2000, 0., 2000.}, "N_{tracks} in subevent B"};
  ConfigurableAxis axisEffPt{"axisEffPt", {500, 0., 5.0}, "Track pT axis for efficiency (must cover [ptTrkMin, ptTrkMax])"};
  ConfigurableAxis axisEffEta{"axisEffEta", {160, -0.8, 0.8}, "Track eta axis for efficiency (must cover |eta| < etaTrkMax)"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    /// check process functions
    std::array<int, 2> processes = {doprocessDplusMl, doprocessD0Ml};
    const int nProcesses = std::accumulate(processes.begin(), processes.end(), 0);
    if (nProcesses > 1) {
      LOGP(fatal, "Only one D-meson process function should be enabled at a time, please check your configuration");
    } else if (nProcesses == 0 && doprocessMcRec == 0 && doprocessMcGen == 0) {
      LOGP(fatal, "No process function enabled");
    }

    const AxisSpec aInvMass{axisInvMass, "Inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec aPt{axisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec aCent{axisCent, "Centrality"};
    const AxisSpec aSign{axisSign, "Sign"};
    const AxisSpec aMlOne{axisMlOne, "Bkg score"};
    const AxisSpec aMlTwo{axisMlTwo, "FD score"};
    const AxisSpec aCandEta{axisCandEta, "Eta"};
    const AxisSpec aMPtTrkA{axisMPtTrkA, "Mean pT of tracks in subevent A (GeV/#it{c})"};
    const AxisSpec aMPtTrkB{axisMPtTrkB, "Mean pT of tracks in subevent B (GeV/#it{c})"};
    // const AxisSpec aPtCandProduct{axisPtCandProduct, "#it{p}_{T}^{cand} * <#it{p}_{T}>^{trks} (GeV^{2}/#it{c}^{2})"}; // to be moved to the TTree
    const AxisSpec aPtTrkProduct{axisPtTrkProduct, "#it{p}_{T}^{trks}(A) * #it{p}_{T}^{trks}(B) (GeV^{2}/#it{c}^{2})"};
    const AxisSpec aNTrkA{axisNTrkA, "N_{tracks} in subevent A"};
    const AxisSpec aNTrkB{axisNTrkB, "N_{tracks} in subevent B"};
    const AxisSpec aEffPt{axisEffPt, "p_{T} (GeV/c)"};
    const AxisSpec aEffEta{axisEffEta, "#eta"};

    // Event-level accumulators (charged hadrons only!)
    registry.add("hEvents", "events vs cent", HistType::kTH1F, {aCent}, true);
    registry.add("pMeanPtA", "<[pT]_A> vs cent (charged hadrons)", HistType::kTProfile, {aCent}, true);
    registry.add("pMeanPtB", "<[pT]_B> vs cent (charged hadrons)", HistType::kTProfile, {aCent}, true);
    registry.add("pMeanPtAmeanPtB", "<[pT]_A*[pT]_B> vs cent (charged hadrons)", HistType::kTProfile, {aCent}, true);

    // QA: multiplicities in A/B
    registry.add("pNA", "<N_A> vs cent (charged hadrons)", HistType::kTProfile, {aCent}, true);
    registry.add("pNB", "<N_B> vs cent (charged hadrons)", HistType::kTProfile, {aCent}, true);

    // charm-bulk correlations
    if (saveCharmBulkCorrelations) {
      registry.add("hCharmBulkCorrelations", "Charm-bulk correlations", HistType::kTHnSparseF, {aInvMass, aCent, aPt, aSign, aMlOne, aMlTwo, aCandEta, aMPtTrkA, aMPtTrkB, aNTrkA, aNTrkB}, true);
      registry.add("hMeanPtTrkAllColls", "Mean pT of charged hadrons for all collisions", HistType::kTHnSparseF, {aCent, aMPtTrkA, aMPtTrkB, aPtTrkProduct, aNTrkA, aNTrkB}, true);
    }

    // MC track-efficiency numerator/denominator + QA
    registry.add("hEffRec", "reconstructed matched primary charged tracks", HistType::kTH2F, {aEffPt, aEffEta}, true);
    registry.add("hEffGen", "generated primary charged tracks", HistType::kTH2F, {aEffPt, aEffEta}, true);
    registry.add("hEffRecFake", "reconstructed fake/secondary tracks (QA)", HistType::kTH2F, {aEffPt, aEffEta}, true);
    registry.add("hEffOutOfRange", "tracks out of efficiency-map range (QA)", HistType::kTH2F, {aEffPt, aEffEta}, true);
    registry.add("hEffInvalid", "tracks with invalid efficiency (QA)", HistType::kTH2F, {aEffPt, aEffEta}, true);

    hfEvSel.init(registry); // collision monitoring

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(ccdbNoLaterThanEff);

    if (usePtDepDcaXY) {
      fPtDepDcaXY = new TF1("ptDepDcaXY", dcaXYFunc.value.c_str(), 0.001, 1000.);
    }

    mCcdbPathEff = ccdbPathEffBase.value + trackEffName.value;
  }; // end init

  // ---------------------------------------------------------------------------
  // Helper functions
  // ---------------------------------------------------------------------------
  /// Check event selections for collision and fill event selection histograms, and revalculate centrality
  /// \param collision is the collision
  template <typename Coll>
  bool isSelectedHfCollision(Coll const& collision, float& cent)
  {
    o2::hf_evsel::HfCollisionRejectionMask collRejMask{};
    if (centEstimator == CentralityEstimator::FT0A) {
      collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0A, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
    } else if (centEstimator == CentralityEstimator::FT0C) {
      collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
    } else if (centEstimator == CentralityEstimator::FT0M) {
      collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
    } else if (centEstimator == CentralityEstimator::FV0A) {
      collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FV0A, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
    } else {
      LOG(fatal) << "Centrality estimator not recognized for collision selection";
    }
    hfEvSel.fillHistograms(collision, collRejMask, cent);
    return collRejMask == 0 && cent >= centralityMin.value && cent < centralityMax.value;
  }

  /// Get candidate mass
  template <DecayChannel channel, typename CandT, typename TCands, typename Trk>
  std::pair<float, float> getCandMassAndSign(const CandT& cand)
  {
    if constexpr (std::is_same_v<TCands, CandDplusDataWMl>) {
      return {HfHelper::invMassDplusToPiKPi(cand), cand.template prong0_as<Trk>().sign()};
    }
    if constexpr (std::is_same_v<TCands, CandD0DataWMl>) {
      if (channel == DecayChannel::D0ToPiK) {
        return {HfHelper::invMassD0ToPiK(cand), cand.isSelD0bar() ? CandD0Type::ReflectedD0bar : CandD0Type::PureD0};
      }
      if (channel == DecayChannel::D0ToKPi) {
        return {HfHelper::invMassD0barToKPi(cand), cand.isSelD0() ? CandD0Type::ReflectedD0 : CandD0Type::PureD0bar};
      }
    }
    LOG(fatal) << "Unsupported candidate type in getCandMassAndSign";
    return {0., 0.}; // unreachable
  }

  /// PWGHF-standard track selection (isGlobalTrackWoDCA) with pT-dependent DCAxy and flat DCAz.
  template <typename T>
  bool selectionTrack(const T& track) const
  {
    if (!track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (usePtDepDcaXY.value && std::abs(track.dcaXY()) >= fPtDepDcaXY->Eval(track.pt())) {
      return false;
    }
    if (std::abs(track.dcaZ()) >= maxDcaZ.value) {
      return false;
    }
    return true;
  }

  /// remove candidate daughters from the mean pT of tracks in A and B (if they are in the respective subevent)
  template <DecayChannel Channel, typename CandT>
  std::pair<float, int> removeDaughtersFromMeanPt(const CandT& cand, float rawMeanPt, int n, const std::vector<int>& trkIDs)
  {
    int removedCount = 0;
    float removedSumPt = 0.f;
    auto removeDaug = [&](int daugID, float daugPt) {
      if (std::binary_search(trkIDs.begin(), trkIDs.end(), daugID)) {
        removedSumPt += daugPt;
        ++removedCount;
      }
    };
    if constexpr (Channel == DecayChannel::DplusToPiKPi) {
      removeDaug(cand.prong0Id(), cand.ptProng0());
      removeDaug(cand.prong1Id(), cand.ptProng1());
      removeDaug(cand.prong2Id(), cand.ptProng2());
    } else if constexpr (Channel == DecayChannel::D0ToPiK || Channel == DecayChannel::D0ToKPi) {
      removeDaug(cand.prong0Id(), cand.ptProng0());
      removeDaug(cand.prong1Id(), cand.ptProng1());
    }
    if (removedCount > 0) {
      double totalSum = static_cast<double>(rawMeanPt) * n;
      return {static_cast<float>((totalSum - removedSumPt) / (n - removedCount)), n - removedCount};
    }
    return {rawMeanPt, n};
  }

  /// Track-efficiency weight w = 1/epsilon(pT, eta), with the mandatory guards (w=1 on invalid).
  float effWeight(float pt, float eta)
  {
    if (!applyTrackEffCorrection || hTrackEffPtEta == nullptr) {
      return 1.f;
    }
    const int bPt = hTrackEffPtEta->GetXaxis()->FindFixBin(pt);
    const int bEta = hTrackEffPtEta->GetYaxis()->FindFixBin(eta);
    if (bPt < 1 || bPt > hTrackEffPtEta->GetNbinsX() ||
        bEta < 1 || bEta > hTrackEffPtEta->GetNbinsY()) {
      registry.fill(HIST("hEffOutOfRange"), pt, eta);
      return 1.f;
    }
    const float eps = hTrackEffPtEta->GetBinContent(bPt, bEta);
    if (!std::isfinite(eps) || eps <= 0.f || eps < effMin) {
      registry.fill(HIST("hEffInvalid"), pt, eta);
      return 1.f;
    }
    return 1.f / eps;
  }

  // ---------------------------------------------------------------------------
  // Event-wise mean pT from charged hadron tracks
  // ---------------------------------------------------------------------------
  template <typename Trk>
  std::tuple<float, int, int> computeMeanPtCharged(Trk const& tracks, float etaMin, float etaMax)
  {
    double sumW = 0.0;
    double sumWpt = 0.0;
    int rawN = 0;

    for (const auto& trk : tracks) {
      if (!selectionTrack(trk)) {
        continue;
      }
      float eta = trk.eta();
      if (eta <= etaMin || eta >= etaMax) {
        continue;
      }
      float pt = trk.pt();
      float w = effWeight(pt, eta);
      sumW += w;
      sumWpt += w * pt;
      ++rawN;
    }

    const int weightedN = static_cast<int>(sumW);
    if (sumW < static_cast<double>(minNTrk.value)) {
      return {std::numeric_limits<float>::quiet_NaN(), weightedN, rawN};
    }
    return {static_cast<float>(sumWpt / sumW), weightedN, rawN};
  }

  template <DecayChannel Channel, typename Cand, typename TCands>
  std::array<float, 2> getMlScores(const Cand& cand) const
  {
    std::array<float, 2> scores{-999.f, -999.f};
    gsl::span<const float> probs;
    if constexpr (std::is_same_v<TCands, CandDplusDataWMl>) {
      probs = cand.mlProbDplusToPiKPi();
    } else if constexpr (Channel == DecayChannel::D0ToPiK) {
      probs = cand.mlProbD0();
    } else if constexpr (Channel == DecayChannel::D0ToKPi) {
      probs = cand.mlProbD0bar();
    } else {
      return scores;
    }
    if (probs.empty()) {
      return scores;
    }
    const auto nScores = std::min(classMl->size(), scores.size());
    for (std::size_t i = 0; i < nScores; ++i) {
      const auto idx = static_cast<std::size_t>(classMl->at(i));
      scores[i] = (idx < probs.size()) ? probs[idx] : -999.f;
    }
    return scores;
  }

  inline bool passMlCut(float ml1, float ml2) const
  {
    if (!applyMlCut.value) {
      return true;
    }
    return (ml1 >= mlOneMin.value && ml1 < mlOneMax.value &&
            ml2 >= mlTwoMin.value && ml2 < mlTwoMax.value);
  }

  /// True if the (absolute) PDG code is a charged hadron/lepton counted in the denominator.
  inline bool isChargedPdg(int pdgCode) const
  {
    const int apdg = std::abs(pdgCode);
    return apdg == PDG_t::kPiPlus || apdg == PDG_t::kKPlus || apdg == PDG_t::kProton ||
           apdg == PDG_t::kElectron || apdg == PDG_t::kMuonMinus;
  }

  /// Load the TH2F track-efficiency map once per run (cached; null -> fatal).
  template <typename Coll>
  void loadTrackEfficiency(Coll const& collision)
  {
    if (!applyTrackEffCorrection) {
      return;
    }
    const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (bc.runNumber() == mLastRunNumber) {
      return;
    }
    hTrackEffPtEta = ccdb->getForTimeStamp<TH2F>(mCcdbPathEff, bc.timestamp());
    if (hTrackEffPtEta == nullptr) {
      LOGF(fatal, "Track-efficiency map '%s' not found in CCDB for run %i", mCcdbPathEff.data(), bc.runNumber());
    }
    mLastRunNumber = bc.runNumber();
  }

  /// Compute the mean pT of the event using charged tracks in two eta-separated subevents, and correlate with D candidates
  /// \param collision is the collision with centrality and event selection
  /// \param candidates are the D candidates to correlate with the mean pT of the event
  /// \param tracks are the tracks used to compute the mean pT of the event
  template <DecayChannel Channel, typename T1, typename Trk>
  void runPtFlucAnalysis(Colls::iterator const& collision,
                         T1 const& candidates,
                         Trk const& tracks)
  {
    float cent{0.f};
    if (!isSelectedHfCollision(collision, cent)) {
      return;
    }

    loadTrackEfficiency(collision);

    // Compute <pT> in two eta-separated subevents using CHARGED tracks only
    const auto [RawMeanPtA, nA, nARaw] = computeMeanPtCharged(tracks, etaTrkAMin.value, etaTrkAMax.value);
    const auto [RawMeanPtB, nB, nBRaw] = computeMeanPtCharged(tracks, etaTrkBMin.value, etaTrkBMax.value);

    // in case double counting for D0
    if constexpr (Channel != DecayChannel::D0ToKPi) {
      registry.fill(HIST("pNA"), cent, static_cast<float>(nARaw));
      registry.fill(HIST("pNB"), cent, static_cast<float>(nBRaw));
      if (std::isfinite(RawMeanPtA) && std::isfinite(RawMeanPtB)) {
        registry.fill(HIST("hEvents"), cent, 1.f);
        registry.fill(HIST("pMeanPtA"), cent, RawMeanPtA);
        registry.fill(HIST("pMeanPtB"), cent, RawMeanPtB);
        registry.fill(HIST("pMeanPtAmeanPtB"), cent, RawMeanPtA * RawMeanPtB);
        if (saveCharmBulkCorrelations) {
          registry.fill(HIST("hMeanPtTrkAllColls"), cent, RawMeanPtA, RawMeanPtB, RawMeanPtA * RawMeanPtB, nA, nB);
        }
      }
    }

    if (!std::isfinite(RawMeanPtA) || !std::isfinite(RawMeanPtB)) {
      return;
    }

    // Collect per-subevent track IDs once, used for daughter removal in the charm-bulk correlations
    std::vector<int> trkIDA;
    std::vector<int> trkIDB;
    trkIDA.reserve(tracks.size() / 2);
    trkIDB.reserve(tracks.size() / 2);
    for (const auto& trk : tracks) {
      if (!selectionTrack(trk)) {
        continue;
      }
      const float eta = trk.eta();
      if (eta > etaTrkAMin.value && eta < etaTrkAMax.value) {
        trkIDA.push_back(trk.globalIndex());
      } else if (eta > etaTrkBMin.value && eta < etaTrkBMax.value) {
        trkIDB.push_back(trk.globalIndex());
      }
    }
    std::sort(trkIDA.begin(), trkIDA.end());
    std::sort(trkIDB.begin(), trkIDB.end());

    // charm-bulk correlations (optional, both eta sides).
    if (saveCharmBulkCorrelations) {
      for (const auto& cand : candidates) {
        const auto [ml1, ml2] = getMlScores<Channel, decltype(cand), T1>(cand);
        if (!passMlCut(ml1, ml2)) {
          continue;
        }
        const float eta = cand.eta();
        const float pt = cand.pt();
        // float candPtProduct{0.f}; // to be moved to the TTree
        float meanPtA = RawMeanPtA;
        float meanPtB = RawMeanPtB;
        int nATrk = nA;
        int nBTrk = nB;
        if (eta < 0. && std::abs(eta) < etaCandMax.value) {
          std::tie(meanPtB, nBTrk) = removeDaughtersFromMeanPt<Channel, decltype(cand)>(cand, RawMeanPtB, nB, trkIDB);
          // candPtProduct = pt * meanPtB; // to be moved to the TTree
        } else if (eta > 0. && eta < etaCandMax.value) {
          std::tie(meanPtA, nATrk) = removeDaughtersFromMeanPt<Channel, decltype(cand)>(cand, RawMeanPtA, nA, trkIDA);
          // candPtProduct = pt * meanPtA; // to be moved to the TTree
        }
        auto [invMass, sign] = getCandMassAndSign<Channel, decltype(cand), T1, Trk>(cand);
        registry.fill(HIST("hCharmBulkCorrelations"), invMass, cent, pt, sign, ml1, ml2, eta, meanPtA, meanPtB, nATrk, nBTrk);
      }
    }

    // Fill candidate info tree (optional, both eta sides)
    if (fillCandInfoTree) {
      for (const auto& cand : candidates) {
        const auto [ml1, ml2] = getMlScores<Channel, decltype(cand), T1>(cand);
        if (!passMlCut(ml1, ml2)) {
          continue;
        }
        const float pt = cand.pt();
        if (pt < ptCandMin.value || pt >= ptCandMax.value) {
          continue;
        }
        const auto [mass, sign] = getCandMassAndSign<Channel, decltype(cand), T1, Trk>(cand);
        (void)sign;
        rowCandInfo(pt, mass, cand.eta(), ml1, ml2); // pt * meanPtA and pt * meanPtB to be moved to this TTree
      }
    }
  } // end runPtFlucAnalysis

  // D0 with ML
  void processD0Ml(Colls::iterator const& collision,
                   CandD0DataWMl const& /*candidatesD0*/,
                   soa::Filtered<TracksWithExtra> const& tracks,
                   aod::BCsWithTimestamps const&)
  {
    auto candsD0ToPiKWMl = selectedD0ToPiKWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsD0ToKPiWMl = selectedD0ToKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runPtFlucAnalysis<DecayChannel::D0ToPiK>(collision, candsD0ToPiKWMl, tracks);
    runPtFlucAnalysis<DecayChannel::D0ToKPi>(collision, candsD0ToKPiWMl, tracks);
  }
  PROCESS_SWITCH(HfTaskPtFlucCharmHadrons, processD0Ml, "Process D0 candidates with ML", false);

  // Dplus with ML
  void processDplusMl(Colls::iterator const& collision,
                      CandDplusDataWMl const& candidatesDplus,
                      soa::Filtered<TracksWithExtra> const& tracks,
                      aod::BCsWithTimestamps const&)
  {
    runPtFlucAnalysis<DecayChannel::DplusToPiKPi>(collision, candidatesDplus, tracks);
  }
  PROCESS_SWITCH(HfTaskPtFlucCharmHadrons, processDplusMl, "Process Dplus candidates with ML", true);

  // MC reconstructed charged tracks (efficiency numerator)
  void processMcRec(aod::McCollision const& mcCollision,
                    aod::McParticles const& mcParticles,
                    soa::SmallGroups<McRecCollisions> const& collisions,
                    McRecTracks const& tracks,
                    aod::BCsWithTimestamps const&)
  {
    (void)mcParticles;
    bool foundSelected = false;
    int64_t collGlobalIndex = -1;
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision() || collision.mcCollisionId() != mcCollision.globalIndex()) {
        continue;
      }
      float cent{0.f};
      if (!isSelectedHfCollision(collision, cent)) {
        continue;
      }
      foundSelected = true;
      collGlobalIndex = collision.globalIndex();
      break;
    }
    if (!foundSelected) {
      return;
    }

    const auto tracksThisCollision = tracks.sliceBy(perMcRecCollision, collGlobalIndex);
    for (const auto& track : tracksThisCollision) {
      if (!track.has_mcParticle()) {
        registry.fill(HIST("hEffRecFake"), track.pt(), track.eta());
        continue;
      }
      auto part = track.mcParticle();
      if (!part.has_mcCollision() || part.mcCollisionId() != mcCollision.globalIndex()) {
        continue; // pileup contamination
      }
      if (!part.isPhysicalPrimary()) {
        registry.fill(HIST("hEffRecFake"), part.pt(), part.eta());
        continue;
      }
      if (!isChargedPdg(part.pdgCode())) {
        continue;
      }
      if (std::abs(track.eta()) >= etaTrkMax || track.pt() <= ptTrkMin || track.pt() >= ptTrkMax ||
          !selectionTrack(track)) {
        continue;
      }
      registry.fill(HIST("hEffRec"), part.pt(), part.eta());
    }
  }
  PROCESS_SWITCH(HfTaskPtFlucCharmHadrons, processMcRec, "MC reconstructed charged tracks (efficiency numerator)", false);

  // MC generated charged tracks (efficiency denominator)
  void processMcGen(aod::McCollision const& mcCollision,
                    aod::McParticles const& mcParticles,
                    soa::SmallGroups<McRecCollisions> const& collisions,
                    aod::BCsWithTimestamps const&)
  {
    bool foundSelected = false;
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision() || collision.mcCollisionId() != mcCollision.globalIndex()) {
        continue;
      }
      float cent{0.f};
      if (!isSelectedHfCollision(collision, cent)) {
        continue;
      }
      foundSelected = true;
      break;
    }
    if (!foundSelected) {
      return;
    }

    // Slice generated particles to this MC collision (denominator must count exactly the
    // same event population as the numerator).
    const auto mcParticlesPerColl = mcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());
    for (const auto& mcParticle : mcParticlesPerColl) {
      if (!mcParticle.has_mcCollision() || mcParticle.mcCollisionId() != mcCollision.globalIndex()) {
        continue; // consistency guard, mirrors processMcRec pileup rejection
      }
      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }
      if (!isChargedPdg(mcParticle.pdgCode())) {
        continue;
      }
      if (mcParticle.pt() <= ptTrkMin || mcParticle.pt() >= ptTrkMax || std::abs(mcParticle.eta()) >= etaTrkMax) {
        continue;
      }
      registry.fill(HIST("hEffGen"), mcParticle.pt(), mcParticle.eta());
    }
  }
  PROCESS_SWITCH(HfTaskPtFlucCharmHadrons, processMcGen, "MC generated charged tracks (efficiency denominator)", false);

}; // End struct HfTaskPtFlucCharmHadrons

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskPtFlucCharmHadrons>(cfgc)};
}
