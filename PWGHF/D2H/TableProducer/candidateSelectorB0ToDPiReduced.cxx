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

/// \file candidateSelectorB0ToDPiReduced.cxx
/// \brief B0 → D- π+ candidate selector
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

struct HfCandidateSelectorB0ToDPiReduced {
  Produces<aod::HfSelB0ToDPi> hfSelB0ToDPiCandidate; // table defined in CandidateSelectionTables.h

  Configurable<float> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // Enable PID
  Configurable<bool> usePid{"usePid", true, "Switch for PID selection at track level"};
  Configurable<bool> acceptPIDNotApplicable{"acceptPIDNotApplicable", true, "Switch to accept Status::NotApplicable [(NotApplicable for one detector) and (NotApplicable or Conditional for the other)] in PID selection"};
  // TPC PID
  Configurable<float> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<float> ptPidTpcMax{"ptPidTpcMax", 20., "Upper bound of track pT for TPC PID"};
  Configurable<float> nSigmaTpcMax{"nSigmaTpcMax", 5., "Nsigma cut on TPC only"};
  Configurable<float> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<float> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<float> ptPidTofMax{"ptPidTofMax", 20., "Upper bound of track pT for TOF PID"};
  Configurable<float> nSigmaTofMax{"nSigmaTofMax", 5., "Nsigma cut on TOF only"};
  Configurable<float> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_b0_to_d_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_b0_to_d_pi::cuts[0], hf_cuts_b0_to_d_pi::nBinsPt, hf_cuts_b0_to_d_pi::nCutVars, hf_cuts_b0_to_d_pi::labelsPt, hf_cuts_b0_to_d_pi::labelsCutVar}, "B0 candidate selection per pT bin"};
  // D-meson ML cuts
  Configurable<std::vector<double>> binsPtDmesMl{"binsPtDmesMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "D-meson pT bin limits for ML cuts"};
  Configurable<LabeledArray<double>> cutsDmesMl{"cutsDmesMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsDmesCutScore}, "D-meson ML cuts per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};

  // check if selectionFlagD (defined in dataCreatorDplusPiReduced.cxx) and usePid configurables are in sync
  bool selectionFlagDAndUsePidInSync = true;
  // variable that will store the value of selectionFlagD (defined in dataCreatorDplusPiReduced.cxx)
  int mySelectionFlagD = -1;

  HfHelper hfHelper;
  TrackSelectorPi selectorPion;

  HistogramRegistry registry{"registry"};

  using TracksPion = soa::Join<HfRedTracks, HfRedTracksPid>;

  void init(InitContext const& initContext)
  {
    std::array<bool, 2> doprocess{doprocessSelection, doprocessSelectionWithDmesMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function for data should be enabled at a time.");
    }

    if (usePid) {
      selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
      selectorPion.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
      selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
      selectorPion.setRangePtTof(ptPidTofMin, ptPidTofMax);
      selectorPion.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
      selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    }

    if (activateQA) {
      constexpr int kNBinsSelections = 1 + SelectionStep::NSelectionSteps;
      std::string labels[kNBinsSelections];
      labels[0] = "No selection";
      labels[1 + SelectionStep::RecoSkims] = "Skims selection";
      labels[1 + SelectionStep::RecoTopol] = "Skims & Topological selections";
      labels[1 + SelectionStep::RecoPID] = "Skims & Topological & PID selections";
      static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
      registry.add("hSelections", "Selections;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSelections, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
        registry.get<TH2>(HIST("hSelections"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }
    }
  }

  /// Main function to perform B0 candidate creation
  /// \param withDmesMl is the flag to use the table with ML scores for the D- daughter (only possible if present in the derived data)
  /// \param hfCandsB0 B0 candidates
  /// \param pionTracks pion tracks
  /// \param configs config inherited from the Dpi data creator
  template <bool withDmesMl, typename Cands>
  void runSelection(Cands const& hfCandsB0,
                    TracksPion const& pionTracks,
                    HfCandB0Configs const& configs)
  {
    // get DplusPi creator configurable
    for (const auto& config : configs) {
      mySelectionFlagD = config.mySelectionFlagD();

      if (usePid && !TESTBIT(mySelectionFlagD, SelectionStep::RecoPID)) {
        selectionFlagDAndUsePidInSync = false;
        LOG(warning) << "PID selections required on B0 daughters (usePid=true) but no PID selections on D candidates were required a priori (selectionFlagD<7). Set selectionFlagD=7 in hf-candidate-creator-b0";
      }
      if (!usePid && TESTBIT(mySelectionFlagD, SelectionStep::RecoPID)) {
        selectionFlagDAndUsePidInSync = false;
        LOG(warning) << "No PID selections required on B0 daughters (usePid=false) but PID selections on D candidates were required a priori (selectionFlagD=7). Set selectionFlagD<7 in hf-candidate-creator-b0";
      }
    }

    for (const auto& hfCandB0 : hfCandsB0) {
      int statusB0ToDPi = 0;
      auto ptCandB0 = hfCandB0.pt();

      // check if flagged as B0 → D π
      if (!TESTBIT(hfCandB0.hfflag(), hf_cand_b0::DecayType::B0ToDPi)) {
        hfSelB0ToDPiCandidate(statusB0ToDPi);
        if (activateQA) {
          registry.fill(HIST("hSelections"), 1, ptCandB0);
        }
        // LOGF(info, "B0 candidate selection failed at hfflag check");
        continue;
      }
      SETBIT(statusB0ToDPi, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusB0ToDPi = 1
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoSkims, ptCandB0);
      }

      // topological cuts
      if (!hfHelper.selectionB0ToDPiTopol(hfCandB0, cuts, binsPt)) {
        hfSelB0ToDPiCandidate(statusB0ToDPi);
        // LOGF(info, "B0 candidate selection failed at topology selection");
        continue;
      }

      if constexpr (withDmesMl) { // we include it in the topological selections
        if (!hfHelper.selectionDmesMlScoresForB(hfCandB0, cutsDmesMl, binsPtDmesMl)) {
          hfSelB0ToDPiCandidate(statusB0ToDPi);
          // LOGF(info, "B0 candidate selection failed at D-meson ML selection");
          continue;
        }
      }

      SETBIT(statusB0ToDPi, SelectionStep::RecoTopol); // RecoTopol = 1 --> statusB0ToDPi = 3
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoTopol, ptCandB0);
      }

      // checking if selectionFlagD and usePid are in sync
      if (!selectionFlagDAndUsePidInSync) {
        hfSelB0ToDPiCandidate(statusB0ToDPi);
        continue;
      }
      // track-level PID selection
      if (usePid) {
        auto trackPi = hfCandB0.template prong1_as<TracksPion>();
        int pidTrackPi = selectorPion.statusTpcAndTof(trackPi);
        if (!hfHelper.selectionB0ToDPiPid(pidTrackPi, acceptPIDNotApplicable.value)) {
          // LOGF(info, "B0 candidate selection failed at PID selection");
          hfSelB0ToDPiCandidate(statusB0ToDPi);
          continue;
        }
        SETBIT(statusB0ToDPi, SelectionStep::RecoPID); // RecoPID = 2 --> statusB0ToDPi = 7
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoPID, ptCandB0);
        }
      }
      hfSelB0ToDPiCandidate(statusB0ToDPi);
      // LOGF(info, "B0 candidate selection passed all selections");
    }
  }

  void processSelection(HfRedCandB0 const& hfCandsB0,
                        TracksPion const& pionTracks,
                        HfCandB0Configs const& configs)
  {
    runSelection<false>(hfCandsB0, pionTracks, configs);
  } // processSelection

  PROCESS_SWITCH(HfCandidateSelectorB0ToDPiReduced, processSelection, "Process selection without ML scores of D mesons", true);

  void processSelectionWithDmesMl(soa::Join<HfRedCandB0, HfRedB0DpMls> const& hfCandsB0,
                                  TracksPion const& pionTracks,
                                  HfCandB0Configs const& configs)
  {
    runSelection<true>(hfCandsB0, pionTracks, configs);
  } // processSelectionWithDmesMl

  PROCESS_SWITCH(HfCandidateSelectorB0ToDPiReduced, processSelectionWithDmesMl, "Process selection with ML scores of D mesons", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorB0ToDPiReduced>(cfgc)};
}
