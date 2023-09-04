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

/// \file candidateSelectorBsToDsPi.cxx
/// \brief Bs → Ds- π+ candidate selector
/// \note adapted from candidateSelectorB0ToDPi.cxx
///
/// \author Phil Stahlhut <phil.lennart.stahlhut@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand_bs; // from CandidateReconstructionTables.h
using namespace o2::analysis;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_bs_to_ds_pi; // from SelectorCuts.h

struct HfCandidateSelectorBsToDsPi {
  Produces<aod::HfSelBsToDsPi> hfSelBsToDsPiCandidate; // table defined in CandidateSelectionTables.h

  // Enable PID
  Configurable<bool> usePid{"usePid", true, "Switch for PID selection at track level"};
  Configurable<bool> acceptPIDNotApplicable{"acceptPIDNotApplicable", true, "Switch to accept Status::NotApplicable [(NotApplicable for one detector) and (NotApplicable or Conditional for the other)] in PID selection"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 20., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 5., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 20., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 5., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_bs_to_ds_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_bs_to_ds_pi::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "Bs candidate selection per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};

  // check if selectionFlagDs (defined in candidateCreatorBs.cxx) and usePid configurables are in sync
  bool selectionFlagDsAndUsePidInSync = true;

  TrackSelectorPi selectorPion;

  using TracksPidWithSel = soa::Join<aod::TracksWExtra, aod::TracksPidPi, aod::TrackSelection>;

  HistogramRegistry registry{"registry"};

  void init(InitContext& initContext)
  {
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

    int selectionFlagDs = -1;
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name.compare("hf-candidate-creator-bs") == 0) {
        for (const auto& option : device.options) {
          if (option.name.compare("selectionFlagDs") == 0) {
            selectionFlagDs = option.defaultValue.get<int>();
            LOGF(info, "selectionFlagDs = %d", selectionFlagDs);
          }
        }
      }
    }

    if (usePid && !TESTBIT(selectionFlagDs, SelectionStep::RecoPID)) {
      selectionFlagDsAndUsePidInSync = false;
      LOG(warning) << "PID selections required on Bs daughters (usePid=true) but no PID selections on Ds candidates were required a priori (selectionFlagDs<7). Set selectionFlagDs=7 in hf-candidate-creator-bs";
    }
    if (!usePid && TESTBIT(selectionFlagDs, SelectionStep::RecoPID)) {
      selectionFlagDsAndUsePidInSync = false;
      LOG(warning) << "No PID selections required on Bs daughters (usePid=false) but PID selections on Ds candidates were required a priori (selectionFlagDs=7). Set selectionFlagDs<7 in hf-candidate-creator-bs";
    }
  }

  void process(aod::HfCandBs const& hfCandsBs,
               TracksPidWithSel const&)
  {
    for (const auto& hfCandBs : hfCandsBs) {
      int statusBsToDsPi = 0;
      auto ptCandBs = hfCandBs.pt();

      // check if flagged as Bs → Ds π
      if (!TESTBIT(hfCandBs.hfflag(), hf_cand_bs::DecayType::BsToDsPi)) {
        hfSelBsToDsPiCandidate(statusBsToDsPi);
        if (activateQA) {
          registry.fill(HIST("hSelections"), 1, ptCandBs);
        }
        continue;
      }
      SETBIT(statusBsToDsPi, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusBsToDsPi = 1
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoSkims, ptCandBs);
      }

      // topological cuts
      if (!hf_sel_candidate_bs::selectionTopol(hfCandBs, cuts, binsPt)) {
        hfSelBsToDsPiCandidate(statusBsToDsPi);
        continue;
      }
      SETBIT(statusBsToDsPi, SelectionStep::RecoTopol); // RecoTopol = 1 --> statusBsToDsPi = 3
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoTopol, ptCandBs);
      }

      // checking if selectionFlagDs and usePid are in sync
      if (!selectionFlagDsAndUsePidInSync) {
        hfSelBsToDsPiCandidate(statusBsToDsPi);
        continue;
      }
      // track-level PID selection
      if (usePid) {
        auto trackPi = hfCandBs.prong1_as<TracksPidWithSel>();
        int pidTrackPi = selectorPion.statusTpcAndTof(trackPi);
        if (!hf_sel_candidate_bs::selectionPID(pidTrackPi, acceptPIDNotApplicable.value)) {
          hfSelBsToDsPiCandidate(statusBsToDsPi);
          continue;
        }
        SETBIT(statusBsToDsPi, SelectionStep::RecoPID); // RecoPID = 2 --> statusBsToDsPi = 7
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoPID, ptCandBs);
        }
      }

      hfSelBsToDsPiCandidate(statusBsToDsPi);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorBsToDsPi>(cfgc)};
}
