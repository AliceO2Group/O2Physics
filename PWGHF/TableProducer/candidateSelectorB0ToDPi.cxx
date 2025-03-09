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

/// \file candidateSelectorB0ToDPi.cxx
/// \brief B0 → D- π+ candidate selector
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include <vector>
#include <string>

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

struct HfCandidateSelectorB0ToDPi {
  Produces<aod::HfSelB0ToDPi> hfSelB0ToDPiCandidate; // table defined in CandidateSelectionTables.h

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
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
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_b0_to_d_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_b0_to_d_pi::Cuts[0], hf_cuts_b0_to_d_pi::NBinsPt, hf_cuts_b0_to_d_pi::NCutVars, hf_cuts_b0_to_d_pi::labelsPt, hf_cuts_b0_to_d_pi::labelsCutVar}, "B0 candidate selection per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  // check if selectionFlagD (defined in candidateCreatorB0.cxx) and usePid configurables are in sync

  bool selectionFlagDAndUsePidInSync = true;
  TrackSelectorPi selectorPion;
  HfHelper hfHelper;

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

    int selectionFlagD = -1;
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name.compare("hf-candidate-creator-b0") == 0) {
        for (const auto& option : device.options) {
          if (option.name.compare("selectionFlagD") == 0) {
            selectionFlagD = option.defaultValue.get<int>();
            LOGF(info, "selectionFlagD = %d", selectionFlagD);
          }
        }
      }
    }

    if (usePid && !TESTBIT(selectionFlagD, SelectionStep::RecoPID)) {
      selectionFlagDAndUsePidInSync = false;
      LOG(warning) << "PID selections required on B0 daughters (usePid=true) but no PID selections on D candidates were required a priori (selectionFlagD<7). Set selectionFlagD=7 in hf-candidate-creator-b0";
    }
    if (!usePid && TESTBIT(selectionFlagD, SelectionStep::RecoPID)) {
      selectionFlagDAndUsePidInSync = false;
      LOG(warning) << "No PID selections required on B0 daughters (usePid=false) but PID selections on D candidates were required a priori (selectionFlagD=7). Set selectionFlagD<7 in hf-candidate-creator-b0";
    }
  }

  void process(aod::HfCandB0 const& hfCandsB0,
               TracksPidWithSel const&)
  {
    for (const auto& hfCandB0 : hfCandsB0) {
      int statusB0ToDPi = 0;
      auto ptCandB0 = hfCandB0.pt();

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
        auto trackPi = hfCandB0.prong1_as<TracksPidWithSel>();
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorB0ToDPi>(cfgc)};
}
