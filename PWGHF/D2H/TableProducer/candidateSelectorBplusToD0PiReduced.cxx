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

/// \file candidateSelectorBplusToD0PiReduced.cxx
/// \brief B+ → D0bar π+ candidate selector
///
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari

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

struct HfCandidateSelectorBplusToD0PiReduced {
  Produces<aod::HfSelBplusToD0Pi> hfSelBplusToD0PiCandidate; // table defined in CandidateSelectionTables.h

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
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_bplus_to_d0_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_bplus_to_d0_pi::cuts[0], hf_cuts_bplus_to_d0_pi::nBinsPt, hf_cuts_bplus_to_d0_pi::nCutVars, hf_cuts_bplus_to_d0_pi::labelsPt, hf_cuts_bplus_to_d0_pi::labelsCutVar}, "B+ candidate selection per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};

  bool selectionFlagDAndUsePidInSync = true;
  // variable that will store the value of selectionFlagD (defined in dataCreatorD0PiReduced.cxx)
  int mySelectionFlagD0 = -1;
  int mySelectionFlagD0bar = -1;

  HfHelper hfHelper;
  TrackSelectorPi selectorPion;

  HistogramRegistry registry{"registry"};

  using TracksPion = soa::Join<HfRedTracks, HfRedTracksPid>;

  void init(InitContext const& initContext)
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
  }

  void process(HfRedCandBplus const& hfCandBs,
               TracksPion const&,
               HfCandBpConfigs const& configs)
  {
    // get DplusPi creator configurable
    for (const auto& config : configs) {
      mySelectionFlagD0 = config.mySelectionFlagD0();
      mySelectionFlagD0bar = config.mySelectionFlagD0bar();

      if ((usePid && !mySelectionFlagD0) || (usePid && !mySelectionFlagD0bar)) {
        selectionFlagDAndUsePidInSync = false;
        LOG(warning) << "PID selections required on B+ daughters (usePid=true) but no PID selections on D candidates were required a priori.";
      }
      if ((!usePid && mySelectionFlagD0) || (!usePid && mySelectionFlagD0bar)) {
        selectionFlagDAndUsePidInSync = false;
        LOG(warning) << "No PID selections required on Bp daughters (usePid=false) but PID selections on D candidates were required a priori.";
      }
    }

    for (const auto& hfCandBp : hfCandBs) {
      int statusBplus = 0;
      auto ptCandBplus = hfCandBp.pt();

      // check if flagged as B+ → D π
      if (!TESTBIT(hfCandBp.hfflag(), hf_cand_bplus::DecayType::BplusToD0Pi)) {
        hfSelBplusToD0PiCandidate(statusBplus);
        if (activateQA) {
          registry.fill(HIST("hSelections"), 1, ptCandBplus);
        }
        // LOGF(info, "B+ candidate selection failed at hfflag check");
        continue;
      }
      SETBIT(statusBplus, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusBplus = 1
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoSkims, ptCandBplus);
      }

      // topological cuts
      if (!hfHelper.selectionBplusToD0PiTopol(hfCandBp, cuts, binsPt)) {
        hfSelBplusToD0PiCandidate(statusBplus);
        // LOGF(info, "B+ candidate selection failed at topology selection");
        continue;
      }
      SETBIT(statusBplus, SelectionStep::RecoTopol); // RecoTopol = 1 --> statusBplus = 3
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoTopol, ptCandBplus);
      }

      // checking if selectionFlagD and usePid are in sync
      if (!selectionFlagDAndUsePidInSync) {
        hfSelBplusToD0PiCandidate(statusBplus);
        continue;
      }
      // track-level PID selection
      if (usePid) {
        auto trackPi = hfCandBp.prong1_as<TracksPion>();
        int pidTrackPi = selectorPion.statusTpcAndTof(trackPi);
        if (!hfHelper.selectionBplusToD0PiPid(pidTrackPi, acceptPIDNotApplicable.value)) {
          // LOGF(info, "B+ candidate selection failed at PID selection");
          hfSelBplusToD0PiCandidate(statusBplus);
          continue;
        }
        SETBIT(statusBplus, SelectionStep::RecoPID); // RecoPID = 2 --> statusBplus = 7
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoPID, ptCandBplus);
        }
      }
      hfSelBplusToD0PiCandidate(statusBplus);
      // LOGF(info, "B+ candidate selection passed all selections");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorBplusToD0PiReduced>(cfgc)};
}
