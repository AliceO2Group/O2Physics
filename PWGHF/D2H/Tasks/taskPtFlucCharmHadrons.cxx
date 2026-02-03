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

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

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

#include <array>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <string>
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

struct HfTaskPtFlucCharmHadrons {
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4)"};
  Configurable<int> selectionFlag{"selectionFlag", 1, "Selection Flag for hadron (e.g. 1 for skimming, 3 for topo. and kine., 7 for PID)"};
  Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality accepted in SP computation"};
  Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality accepted in SP computation"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indices of BDT scores to be stored. Two indexes max."};

  Configurable<int> cfgITScluster{"cfgITScluster", 5, "Number of ITS cluster"};
  Configurable<int> cfgITSbarrel{"cfgITSbarrel", 1, "Number of ITS inner barrel"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.1f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.1f, "DCAz range for tracks"};

  // Event selection
  // Mean-pT (charged hadrons) definition: subevents with eta gap
  Configurable<float> etaAMin{"etaAMin", -0.8f, "A: negative eta min"};
  Configurable<float> etaAMax{"etaAMax", 0.0f, "A: negative eta max (0)"};
  Configurable<float> etaBMin{"etaBMin", 0.4f, "B: positive eta min (eta_min)"}; // or set from etaMinGap in code
  Configurable<float> etaBMax{"etaBMax", 0.8f, "B: positive eta max"};

  Configurable<float> ptTrkMin{"ptTrkMin", 0.2f, "Track pT min for <pT> (charged hadrons)"};
  Configurable<float> ptTrkMax{"ptTrkMax", 2.0f, "Track pT max for <pT> (charged hadrons)"};
  Configurable<float> ptTrkMinD{"ptTrkMinD", 0.2f, "Track pT min for D"};
  Configurable<float> ptTrkMaxD{"ptTrkMaxD", 2.0f, "Track pT max for D"};
  Configurable<int> minNTrk{"minNTrk", 5, "Min charged tracks in each subevent to compute <pT>"};

  // Use ML in denominator AND numerator consistently
  Configurable<bool> applyMlCut{"applyMlCut", true, "Apply ML window for candidate selection (also in ND denominator)"};

  // ML window (assumes you store 2 ML outputs: classMl[0], classMl[1])
  Configurable<float> mlOneMin{"mlOneMin", 0.f, "ML1 min"};
  Configurable<float> mlOneMax{"mlOneMax", 1.f, "ML1 max"};
  Configurable<float> mlTwoMin{"mlTwoMin", 0.f, "ML2 min"};
  Configurable<float> mlTwoMax{"mlTwoMax", 1.f, "ML2 max"};

  // Optional: require D candidate to be in subevent A acceptance (recommended for long-range)
  Configurable<bool> requireCandInA{"requireCandInA", true, "Require D candidate eta inside subevent A window"};

  SliceCache cache;
  HfEventSelection hfEvSel; // event selection and monitoring
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  using CandDplusDataWMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandD0DataWMl = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;
  Filter filterSelectD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;

  Partition<CandD0DataWMl> selectedD0ToPiKWMl = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
  Partition<CandD0DataWMl> selectedD0ToKPiWMl = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;

  // -------------------------
  // Axes / histograms
  // -------------------------
  ConfigurableAxis axisInvMass{"axisInvMass", {100, 1.78, 2.05}, "Inv mass axis"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0}, "Candidate pT axis"};
  ConfigurableAxis axisCent{"axisCent", {VARIABLE_WIDTH, 0.0, 10.0, 40.0, 80.0}, "Centrality axis"};
  ConfigurableAxis axisSign{"axisSign", {VARIABLE_WIDTH, -1.0, 1.0}, "Sign axis"};
  ConfigurableAxis thnConfigAxisMlOne{"thnConfigAxisMlOne", {1000, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisMlTwo{"thnConfigAxisMlTwo", {1000, 0., 1.}, ""};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {

    /// check process functions
    std::array<int, 2> processes = {doprocessDplusMl, doprocessD0Ml};
    const int nProcesses = std::accumulate(processes.begin(), processes.end(), 0);
    if (nProcesses > 1) {
      LOGP(fatal, "Only one process function should be enabled at a time, please check your configuration");
    } else if (nProcesses == 0) {
      LOGP(fatal, "No process function enabled");
    }

    const AxisSpec aInvMass{axisInvMass, "Inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec aPt{axisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec aCent{axisCent, "Centrality"};
    const AxisSpec aSign{axisSign, "Sign"};
    const AxisSpec thnAxisMlOne{thnConfigAxisMlOne, "Bkg score"};
    const AxisSpec thnAxisMlTwo{thnConfigAxisMlTwo, "FD score"};

    // Event-level accumulators (charged hadrons only!)
    registry.add("hEvents", "events vs cent", HistType::kTH1F, {aCent}, true);
    registry.add("hEventsD", "eventsD vs cent", HistType::kTH1F, {aCent}, true);
    // These are sufficient to compute sigma_<pT> offline in each centrality bin.
    registry.add("pMeanPtA", "<[pT]_A> vs cent (charged hadrons)", HistType::kTProfile, {aCent}, true);
    registry.add("pMeanPtB", "<[pT]_B> vs cent (charged hadrons)", HistType::kTProfile, {aCent}, true);
    registry.add("pMeanPtAmeanPtB", "<[pT]_A*[pT]_B> vs cent (charged hadrons)", HistType::kTProfile, {aCent}, true);

    // QA: multiplicities in A/B
    registry.add("pNA", "<N_A> vs cent (charged hadrons)", HistType::kTProfile, {aCent}, true);
    registry.add("pNB", "<N_B> vs cent (charged hadrons)", HistType::kTProfile, {aCent}, true);

    // Candidate mass distributions
    registry.add("hD_mass", "D candidate mass (weight=1)", HistType::kTHnSparseF, {aInvMass, aCent, aPt, aSign, thnAxisMlOne, thnAxisMlTwo}, true);
    registry.add("hD_f", "hD_f", HistType::kTHnSparseF, {aInvMass, aCent, aPt, aSign, thnAxisMlOne, thnAxisMlTwo}, true);
    registry.add("hD_f_wPtB", "hD_f_wPtB", HistType::kTHnSparseF, {aInvMass, aCent, aPt, aSign, thnAxisMlOne, thnAxisMlTwo}, true);

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }; // end init

  /// Get the centrality
  /// \param collision is the collision with the centrality information
  double getCentrality(Colls::iterator const& collision)
  {
    double cent = -999.;
    switch (centEstimator) {
      case CentralityEstimator::FV0A:
        cent = collision.centFV0A();
        break;
      case CentralityEstimator::FT0M:
        cent = collision.centFT0M();
        break;
      case CentralityEstimator::FT0A:
        cent = collision.centFT0A();
        break;
      case CentralityEstimator::FT0C:
        cent = collision.centFT0C();
        break;
      default:
        LOG(warning) << "Centrality estimator not valid. Possible values are V0A, T0M, T0A, T0C. Fallback to V0A";
        cent = collision.centFV0A();
        break;
    }
    return cent;
  }

  template <typename T>
  bool selectionTrack(const T& candidate) const
  {
    if (!(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster.value && candidate.tpcNClsFound() > cfgTPCcluster.value && candidate.itsNClsInnerBarrel() >= cfgITSbarrel.value && std::abs(candidate.dcaXY()) < cfgCutDCAxy.value && std::abs(candidate.dcaZ()) < cfgCutDCAz.value)) {
      return false;
    }
    return true;
  }

  // -------------------------
  // Event-wise mean pT from charged hadron tracks
  // -------------------------
  template <typename Trk>
  std::pair<float, int> computeMeanPtCharged(Trk const& tracks, float etaMin, float etaMax) const
  {
    double sumPt = 0.0;
    int n = 0;

    for (const auto& trk : tracks) {
      if (!selectionTrack(trk))
        continue;
      float eta = trk.eta();
      if (eta <= etaMin || eta >= etaMax) {
        continue;
      }
      float pt = trk.pt();
      if (pt < ptTrkMin.value || pt >= ptTrkMax.value) {
        continue;
      }
      sumPt += pt;
      ++n;
    }

    if (n < minNTrk.value) {
      return {std::numeric_limits<float>::quiet_NaN(), n};
    }
    return {static_cast<float>(sumPt / n), n};
  }

  // Helper: candidate eta-in-A cut
  template <typename Cand>
  bool passCandInA(const Cand& cand) const
  {
    if (!requireCandInA.value) {
      return true;
    }
    const float eta = cand.eta();
    return (eta > etaAMin.value && eta < etaAMax.value);
  }

  template <DecayChannel Channel, typename Cand, typename TCands>
  std::array<float, 2> getMlScores(const Cand& cand) const
  {
    float ml1 = -999.f;
    float ml2 = -999.f;

    if constexpr (std::is_same_v<TCands, CandDplusDataWMl>) {
      // D+
      const auto& probs = cand.mlProbDplusToPiKPi();
      ml1 = probs[classMl->at(0)];
      ml2 = probs[classMl->at(1)];
    } else {
      // D0 / D0bar depends on Channel
      if constexpr (Channel == DecayChannel::D0ToPiK) {
        const auto& probs = cand.mlProbD0();
        ml1 = probs[classMl->at(0)];
        ml2 = probs[classMl->at(1)];
      } else if constexpr (Channel == DecayChannel::D0ToKPi) {
        const auto& probs = cand.mlProbD0bar();
        ml1 = probs[classMl->at(0)];
        ml2 = probs[classMl->at(1)];
      }
    }

    return {ml1, ml2};
  }

  inline bool passMlCut(float ml1, float ml2) const
  {
    if (!applyMlCut.value) {
      return true;
    }
    return (ml1 >= mlOneMin.value && ml1 < mlOneMax.value &&
            ml2 >= mlTwoMin.value && ml2 < mlTwoMax.value);
  }

  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param candidates are the selected candidates
  template <DecayChannel Channel, typename T1, typename Trk>
  void runPtFlucAnalysis(Colls::iterator const& collision,
                         T1 const& candidates,
                         Trk const& tracks)
  {
    float cent = getCentrality(collision);
    if (cent < centralityMin || cent > centralityMax) {
      return;
    }

    // Compute <pT> in two eta-separated subevents using CHARGED tracks only
    auto [meanPtA, nA] = computeMeanPtCharged(tracks, etaAMin.value, etaAMax.value);
    auto [meanPtB, nB] = computeMeanPtCharged(tracks, etaBMin.value, etaBMax.value);

    // QA
    registry.fill(HIST("pNA"), cent, static_cast<float>(nA));
    registry.fill(HIST("pNB"), cent, static_cast<float>(nB));

    if (!std::isfinite(meanPtA) || !std::isfinite(meanPtB)) {
      return;
    }

    registry.fill(HIST("hEvents"), cent, 1.f);

    registry.fill(HIST("pMeanPtA"), cent, meanPtA);
    registry.fill(HIST("pMeanPtB"), cent, meanPtB);
    registry.fill(HIST("pMeanPtAmeanPtB"), cent, meanPtA * meanPtB);

    int nDcandTotA = 0;

    for (const auto& cand : candidates) {
      if (!passCandInA(cand)) {
        continue;
      }

      // compute ML exactly like numerator
      const auto [ml1, ml2] = getMlScores<Channel, decltype(cand), T1>(cand);

      // apply same ML cut
      if (!passMlCut(ml1, ml2)) {
        continue;
      }
      float pt = cand.pt();
      if (pt < ptTrkMinD.value || pt >= ptTrkMaxD.value) {
        continue;
      }

      ++nDcandTotA;
    }

    if (nDcandTotA <= 0) {
      return; // cannot build fraction
    }

    registry.fill(HIST("hEventsD"), cent, 1.f);
    const float invND = 1.f / static_cast<float>(nDcandTotA);

    for (const auto& candidate : candidates) {
      if (!passCandInA(candidate)) {
        continue;
      }

      // ML first (same definition)
      const auto [ml1, ml2] = getMlScores<Channel, decltype(candidate), T1>(candidate);
      if (!passMlCut(ml1, ml2)) {
        continue;
      }

      // compute mass
      float massCand = 0.;
      float signCand = 0.;
      if constexpr (std::is_same_v<T1, CandDplusDataWMl>) {
        massCand = HfHelper::invMassDplusToPiKPi(candidate);
        auto trackprong0 = candidate.template prong0_as<Trk>();
        signCand = trackprong0.sign();
      } else {
        if constexpr (Channel == DecayChannel::D0ToPiK) {
          massCand = HfHelper::invMassD0ToPiK(candidate);
          signCand = candidate.isSelD0bar() ? 3 : 1; // 3: reflected D0bar, 1: pure D0 excluding reflected D0bar
        } else if constexpr (Channel == DecayChannel::D0ToKPi) {
          massCand = HfHelper::invMassD0barToKPi(candidate);
          signCand = candidate.isSelD0() ? 3 : 2; // 3: reflected D0, 2: pure D0bar excluding reflected D0
        }
      }

      const double ptCand = candidate.pt();
      if (ptCand < ptTrkMinD.value || ptCand >= ptTrkMaxD.value) {
        continue;
      }

      registry.fill(HIST("hD_mass"), massCand, cent, ptCand, signCand, ml1, ml2, 1.f);
      registry.fill(HIST("hD_f"), massCand, cent, ptCand, signCand, ml1, ml2, invND);
      registry.fill(HIST("hD_f_wPtB"), massCand, cent, ptCand, signCand, ml1, ml2, meanPtB * invND);
    }
  }

  // D0 with ML
  void processD0Ml(Colls::iterator const& collision,
                   CandD0DataWMl const& /*candidatesD0*/,
                   TracksWithExtra const& tracks)
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
                      TracksWithExtra const& tracks)
  {
    runPtFlucAnalysis<DecayChannel::DplusToPiKPi>(collision, candidatesDplus, tracks);
  }
  PROCESS_SWITCH(HfTaskPtFlucCharmHadrons, processDplusMl, "Process Dplus candidates with ML", true);

}; // End struct HfTaskPtFlucCharmHadrons

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskPtFlucCharmHadrons>(cfgc)};
}
