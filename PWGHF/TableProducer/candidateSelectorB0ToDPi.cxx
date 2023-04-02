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

#include "Common/Core/TrackSelectorPID.h"
#include "Framework/runDataProcessing.h"
// #include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand_b0; // from CandidateReconstructionTables.h
using namespace o2::analysis;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_b0_to_d_pi; // from SelectorCuts.h

// FIXME: store B0 creator configurable (until https://alice.its.cern.ch/jira/browse/O2-3582 solved)
namespace o2::aod
{
namespace hf_cand_b0_config
{
DECLARE_SOA_COLUMN(MySelectionFlagD, mySelectionFlagD, int);
} // namespace hf_cand_b0_config

DECLARE_SOA_TABLE(HfCandB0Config, "AOD", "HFCANDB0CONFIG", //!
                  hf_cand_b0_config::MySelectionFlagD);
} // namespace o2::aod

struct HfCandidateSelectorB0ToDPi {
  Produces<aod::HfSelB0ToDPi> hfSelB0ToDPiCandidate; // table defined in CandidateSelectionTables.h

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // Enable PID
  Configurable<bool> usePid{"usePid", true, "Switch for PID selection at track level"};
  Configurable<bool> acceptPIDNotApplicable{"acceptPIDNotApplicable", true, "Switch to accept Status::PIDNotApplicable [(NotApplicable for one detector) and (NotApplicable or Conditional for the other)] in PID selection"};
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
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_b0_to_d_pi::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "B0 candidate selection per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  // check if selectionFlagD (defined in candidateCreatorB0.cxx) and usePid configurables are in sync
  bool selectionFlagDAndUsePidInSync = true;
  // FIXME: store B0 creator configurable (until https://alice.its.cern.ch/jira/browse/O2-3582 solved)
  int mySelectionFlagD = -1;

  TrackSelectorPID selectorPion;

  using TracksPIDWithSel = soa::Join<aod::BigTracksPIDExtended, aod::TrackSelection>;

  HistogramRegistry registry{"registry"};

  void init(InitContext const& initContext)
  {
    if (usePid) {
      selectorPion.setPDG(kPiPlus);
      selectorPion.setRangePtTPC(ptPidTpcMin, ptPidTpcMax);
      selectorPion.setRangeNSigmaTPC(-nSigmaTpcMax, nSigmaTpcMax);
      selectorPion.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
      selectorPion.setRangePtTOF(ptPidTofMin, ptPidTofMax);
      selectorPion.setRangeNSigmaTOF(-nSigmaTofMax, nSigmaTofMax);
      selectorPion.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
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

    // FIXME: will be uncommented once https://alice.its.cern.ch/jira/browse/O2-3582 is solved
    /*int selectionFlagD = -1;
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      if (device.name.compare("hf-candidate-creator-b0") == 0) {
        for (auto const& option : device.options) {
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
    }*/
  }

  /// Apply topological cuts as defined in SelectorCuts.h
  /// \param hfCandB0 is the B0 candidate
  /// \param hfCandD is prong1 of B0 candidate
  /// \param trackPi is prong1 of B0 candidate
  /// \return true if candidate passes all selections
  template <typename T1, typename T2, typename T3>
  bool selectionTopol(const T1& hfCandB0, const T2& hfCandD, const T3& trackPi)
  {
    auto candpT = hfCandB0.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      // LOGF(info, "B0 topol selection failed at getpTBin");
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptCandMin || candpT >= ptCandMax) {
      return false;
    }

    // B0 mass cut
    if (std::abs(invMassB0ToDPi(hfCandB0) - RecoDecay::getMassPDG(pdg::Code::kB0)) > cuts->get(pTBin, "m")) {
      // Printf("B0 topol selection failed at mass diff check");
      return false;
    }

    // pion pt
    if (trackPi.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    // D- pt
    if (hfCandD.pt() < cuts->get(pTBin, "pT D")) {
      return false;
    }

    /*
    // D mass cut | already applied in candidateSelectorDplusToPiKPi.cxx
    if (std::abs(hf_cand_3prong::invMassDplusToPiKPi(hfCandD) - RecoDecay::getMassPDG(pdg::Code::kDMinus)) > cuts->get(pTBin, "DeltaMD")) {
      return false;
    }
    */

    // B0 Decay length
    if (hfCandB0.decayLength() < cuts->get(pTBin, "B0 decLen")) {
      return false;
    }

    // B0 Decay length XY
    if (hfCandB0.decayLengthXY() < cuts->get(pTBin, "B0 decLenXY")) {
      return false;
    }

    // B0 chi2PCA cut
    if (hfCandB0.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
      return false;
    }

    // B0 CPA cut
    if (hfCandB0.cpa() < cuts->get(pTBin, "CPA")) {
      return false;
    }

    // d0 of pi
    if (std::abs(hfCandB0.impactParameter1()) < cuts->get(pTBin, "d0 Pi")) {
      return false;
    }

    // d0 of D
    if (std::abs(hfCandB0.impactParameter0()) < cuts->get(pTBin, "d0 D")) {
      return false;
    }

    return true;
  }

  /// Apply PID selection
  /// \param pidTrackPi is the PID status of trackPi (prong1 of B0 candidate)
  /// \return true if prong1 of B0 candidate passes all selections
  template <typename T = int>
  bool selectionPID(const T& pidTrackPi)
  {
    if (!acceptPIDNotApplicable && pidTrackPi != TrackSelectorPID::Status::PIDAccepted) {
      return false;
    }
    if (acceptPIDNotApplicable && pidTrackPi == TrackSelectorPID::Status::PIDRejected) {
      return false;
    }

    return true;
  }

  void process(aod::HfCandB0 const& hfCandsB0, soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi> const&, TracksPIDWithSel const&,
               HfCandB0Config const& configs)
  {
    // FIXME: get B0 creator configurable (until https://alice.its.cern.ch/jira/browse/O2-3582 solved)
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
          registry.fill(HIST("hSelections"), 0, ptCandB0);
        }
        // LOGF(info, "B0 candidate selection failed at hfflag check");
        continue;
      }
      SETBIT(statusB0ToDPi, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusB0ToDPi = 1
      if (activateQA) {
        registry.fill(HIST("hSelections"), 1 + SelectionStep::RecoSkims, ptCandB0);
      }

      auto candD = hfCandB0.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>();
      auto trackPi = hfCandB0.prong1_as<TracksPIDWithSel>();

      // topological cuts
      if (!selectionTopol(hfCandB0, candD, trackPi)) {
        hfSelB0ToDPiCandidate(statusB0ToDPi);
        // LOGF(info, "B0 candidate selection failed at topology selection");
        continue;
      }
      SETBIT(statusB0ToDPi, SelectionStep::RecoTopol); // RecoTopol = 1 --> statusB0ToDPi = 3
      if (activateQA) {
        registry.fill(HIST("hSelections"), 1 + SelectionStep::RecoTopol, ptCandB0);
      }

      // checking if selectionFlagD and usePid are in sync
      if (!selectionFlagDAndUsePidInSync) {
        hfSelB0ToDPiCandidate(statusB0ToDPi);
        continue;
      }
      // track-level PID selection
      if (usePid) {
        int pidTrackPi = selectorPion.getStatusTrackPIDTpcAndTof(trackPi);
        if (!selectionPID(pidTrackPi)) {
          // LOGF(info, "B0 candidate selection failed at PID selection");
          hfSelB0ToDPiCandidate(statusB0ToDPi);
          continue;
        }
        SETBIT(statusB0ToDPi, SelectionStep::RecoPID); // RecoPID = 2 --> statusB0ToDPi = 7
        if (activateQA) {
          registry.fill(HIST("hSelections"), 1 + SelectionStep::RecoPID, ptCandB0);
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
