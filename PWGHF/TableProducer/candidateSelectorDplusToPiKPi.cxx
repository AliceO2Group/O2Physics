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

/// \file candidateSelectorDplusToPiKPi.cxx
/// \brief D± → π± K∓ π± selection task
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, Politecnico and INFN Torino
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::analysis::hf_cuts_dplus_to_pi_k_pi;

/// Struct for applying Dplus to piKpi selection cuts
struct HfCandidateSelectorDplusToPiKPi {
  Produces<aod::HfSelDplusToPiKPi> hfSelDplusToPiKPiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 1., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  // PID option
  Configurable<bool> acceptPIDNotApplicable{"acceptPIDNotApplicable", true, "Switch to accept Status::PIDNotApplicable [(NotApplicable for one detector) and (NotApplicable or Conditional for the other)] in PID selection"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 20., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 20., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_dplus_to_pi_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_dplus_to_pi_k_pi::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "Dplus candidate selection per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};

  TrackSelectorPID selectorPion;
  TrackSelectorPID selectorKaon;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    selectorPion.setPDG(kPiPlus);
    selectorPion.setRangePtTPC(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTPC(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangePtTOF(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTOF(-nSigmaTofMax, nSigmaTofMax);

    selectorKaon = selectorPion;
    selectorKaon.setPDG(kKPlus);

    if (activateQA) {
      constexpr int kNBinsSelections = 1 + aod::SelectionStep::NSelectionSteps;
      std::string labels[kNBinsSelections];
      labels[0] = "No selection";
      labels[1 + aod::SelectionStep::RecoSkims] = "Skims selection";
      labels[1 + aod::SelectionStep::RecoTopol] = "Skims & Topological selections";
      labels[1 + aod::SelectionStep::RecoPID] = "Skims & Topological & PID selections";
      static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
      registry.add("hSelections", "Selections;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSelections, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
        registry.get<TH2>(HIST("hSelections"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }
    }
  }

  /*
  /// Selection on goodness of daughter tracks
  /// \note should be applied at candidate selection
  /// \param track is daughter track
  /// \return true if track is good
  template <typename T>
  bool daughterSelection(const T& track)
  {
    if (track.tpcNClsFound() == 0) {
      return false; //is it clusters findable or found - need to check
    }
    return true;
  }
  */

  /// Candidate selections
  /// \param candidate is candidate
  /// \param trackPion1 is the first track with the pion hypothesis
  /// \param trackKaon is the track with the kaon hypothesis
  /// \param trackPion2 is the second track with the pion hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selection(const T1& candidate, const T2& trackPion1, const T2& trackKaon, const T2& trackPion2)
  {
    auto candpT = candidate.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }
    // check that the candidate pT is within the analysis range
    if (candpT < ptCandMin || candpT > ptCandMax) {
      return false;
    }
    // cut on daughter pT
    if (trackPion1.pt() < cuts->get(pTBin, "pT Pi") || trackKaon.pt() < cuts->get(pTBin, "pT K") || trackPion2.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    // invariant-mass cut
    if (std::abs(invMassDplusToPiKPi(candidate) - RecoDecay::getMassPDG(pdg::Code::kDPlus)) > cuts->get(pTBin, "deltaM")) {
      return false;
    }
    if (candidate.decayLength() < cuts->get(pTBin, "decay length")) {
      return false;
    }
    if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "normalized decay length XY")) {
      return false;
    }
    if (candidate.cpa() < cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }
    if (candidate.cpaXY() < cuts->get(pTBin, "cos pointing angle XY")) {
      return false;
    }
    if (std::abs(candidate.maxNormalisedDeltaIP()) > cuts->get(pTBin, "max normalized deltaIP")) {
      return false;
    }
    return true;
  }

  /// Apply PID selection
  /// \param pidTrackPos1Pion is the PID status of trackPos1 (prong0 of D candidate)
  /// \param pidTrackNegKaon is the PID status of trackNeg (prong1 of D candidate)
  /// \param pidTrackPos2Pion is the PID status of trackPos2 (prong2 of D candidate)
  /// \return true if prongs pass all selections
  template <typename T = int>
  bool selectionPID(const T& pidTrackPos1Pion, const T& pidTrackNegKaon, const T& pidTrackPos2Pion)
  {
    if (!acceptPIDNotApplicable &&
        (pidTrackPos1Pion != TrackSelectorPID::Status::PIDAccepted ||
         pidTrackNegKaon != TrackSelectorPID::Status::PIDAccepted ||
         pidTrackPos2Pion != TrackSelectorPID::Status::PIDAccepted)) {
      return false;
    }
    if (acceptPIDNotApplicable &&
        (pidTrackPos1Pion == TrackSelectorPID::Status::PIDRejected ||
         pidTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
         pidTrackPos2Pion == TrackSelectorPID::Status::PIDRejected)) {
      return false;
    }

    return true;
  }

  void process(aod::HfCand3Prong const& candidates, aod::BigTracksPID const&)
  {
    // looping over 3-prong candidates
    for (auto& candidate : candidates) {

      // final selection flag:
      auto statusDplusToPiKPi = 0;

      auto ptCand = candidate.pt();

      if (!TESTBIT(candidate.hfflag(), DecayType::DplusToPiKPi)) {
        hfSelDplusToPiKPiCandidate(statusDplusToPiKPi);
        if (activateQA) {
          registry.fill(HIST("hSelections"), 0, ptCand);
        }
        continue;
      }
      SETBIT(statusDplusToPiKPi, aod::SelectionStep::RecoSkims);
      if (activateQA) {
        registry.fill(HIST("hSelections"), 1 + aod::SelectionStep::RecoSkims, ptCand);
      }

      auto trackPos1 = candidate.prong0_as<aod::BigTracksPID>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<aod::BigTracksPID>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<aod::BigTracksPID>(); // positive daughter (negative for the antiparticles)

      /*
      // daughter track validity selection
      if (!daughterSelection(trackPos1) ||
          !daughterSelection(trackNeg) ||
          !daughterSelection(trackPos2)) {
        hfSelDplusToPiKPiCandidate(statusDplusToPiKPi);
        continue;
      }
      */

      // topological selection
      if (!selection(candidate, trackPos1, trackNeg, trackPos2)) {
        hfSelDplusToPiKPiCandidate(statusDplusToPiKPi);
        continue;
      }
      SETBIT(statusDplusToPiKPi, aod::SelectionStep::RecoTopol);
      if (activateQA) {
        registry.fill(HIST("hSelections"), 1 + aod::SelectionStep::RecoTopol, ptCand);
      }

      // track-level PID selection
      int pidTrackPos1Pion = selectorPion.getStatusTrackPIDTpcAndTof(trackPos1);
      int pidTrackNegKaon = selectorKaon.getStatusTrackPIDTpcAndTof(trackNeg);
      int pidTrackPos2Pion = selectorPion.getStatusTrackPIDTpcAndTof(trackPos2);

      if (!selectionPID(pidTrackPos1Pion, pidTrackNegKaon, pidTrackPos2Pion)) { // exclude D±
        hfSelDplusToPiKPiCandidate(statusDplusToPiKPi);
        continue;
      }
      SETBIT(statusDplusToPiKPi, aod::SelectionStep::RecoPID);
      if (activateQA) {
        registry.fill(HIST("hSelections"), 1 + aod::SelectionStep::RecoPID, ptCand);
      }

      hfSelDplusToPiKPiCandidate(statusDplusToPiKPi);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorDplusToPiKPi>(cfgc)};
}
