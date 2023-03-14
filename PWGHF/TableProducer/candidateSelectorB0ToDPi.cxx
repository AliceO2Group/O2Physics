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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand_b0; // from CandidateReconstructionTables.h
using namespace o2::analysis;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_b0_to_d_pi; // from SelectorCuts.h

struct HfCandidateSelectorB0ToDPi {
  Produces<aod::HfSelB0ToDPi> hfSelB0ToDPiCandidate; // table defined in CandidateSelectionTables.h

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // Enable PID
  Configurable<bool> usePid{"usePid", true, "Switch for PID selection at track level"};
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

  // Apply topological cuts as defined in SelectorCuts.h; return true if candidate passes all selections
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

    // D mass cut
    if (std::abs(hf_cand_3prong::invMassDplusToPiKPi(hfCandD) - RecoDecay::getMassPDG(pdg::Code::kDMinus)) > cuts->get(pTBin, "DeltaMD")) {
      return false;
    }

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

  // Apply PID selection; return true if candidate passes all selections
  template <typename T>
  bool selectionPID(const T& pidTrackPi)
  {
    if (pidTrackPi != TrackSelectorPID::Status::PIDAccepted) {
      return false;
    }

    return true;
  }

  using TracksPIDWithSel = soa::Join<aod::BigTracksPIDExtended, aod::TrackSelection>;

  void process(aod::HfCandB0 const& hfCandsB0, soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi> const&, TracksPIDWithSel const&)
  {
    int statusB0ToDPi = 0;

    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTPC(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPion.setRangePtTOF(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTOF(-nSigmaTofMax, nSigmaTofMax);
    selectorPion.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);

    for (const auto& hfCandB0 : hfCandsB0) {
      // check if flagged as B0 → D π
      if (!TESTBIT(hfCandB0.hfflag(), hf_cand_b0::DecayType::B0ToDPi)) {
        hfSelB0ToDPiCandidate(statusB0ToDPi);
        // LOGF(info, "B0 candidate selection failed at hfflag check");
        continue;
      }
      SETBIT(statusB0ToDPi, aod::SelectionStep::RecoSkims); // RecoSkims = 0 --> statusB0ToDPi = 1

      auto candD = hfCandB0.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>();
      auto trackPi = hfCandB0.prong1_as<TracksPIDWithSel>();

      // topological cuts
      if (!selectionTopol(hfCandB0, candD, trackPi)) {
        hfSelB0ToDPiCandidate(statusB0ToDPi);
        // LOGF(info, "B0 candidate selection failed at topology selection");
        continue;
      }
      SETBIT(statusB0ToDPi, aod::SelectionStep::RecoTopol); // RecoTopol = 1 --> statusB0ToDPi = 3

      // track-level PID selection
      if (usePid) {
        if (!TESTBIT(candD.isSelDplusToPiKPi(), aod::SelectionStep::RecoPID)) { // safety
          LOG(warning) << "PID not enabled in DplusToPiKPi selector. Set selectionFlagD=7 in hf-candidate-creator-b0";
          hfSelB0ToDPiCandidate(statusB0ToDPi);
          continue;
        }
        int pidTrackPi = selectorPion.getStatusTrackPIDTpcOrTof(trackPi);
        if (!selectionPID(pidTrackPi)) {
          // LOGF(info, "B0 candidate selection failed at PID selection");
          hfSelB0ToDPiCandidate(statusB0ToDPi);
          continue;
        }
        SETBIT(statusB0ToDPi, aod::SelectionStep::RecoPID); // RecoPID = 2 --> statusB0ToDPi = 7
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
