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

/// \file candidateSelectorJpsi.cxx
/// \brief J/ψ → e+ e−, μ+ μ− selection task
///
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"

#include "ALICE3/DataModel/MID.h"
#include "ALICE3/DataModel/RICH.h"
#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

namespace o2::aod
{
namespace hf_track_index_alice3_pid
{
DECLARE_SOA_INDEX_COLUMN(Track, track); //!
DECLARE_SOA_INDEX_COLUMN(RICH, rich);   //!
DECLARE_SOA_INDEX_COLUMN(MID, mid);     //!
} // namespace hf_track_index_alice3_pid

DECLARE_SOA_INDEX_TABLE_USER(HfTrackIndexALICE3PID, Tracks, "HFTRKIDXA3PID", //!
                             hf_track_index_alice3_pid::TrackId,
                             hf_track_index_alice3_pid::RICHId,
                             hf_track_index_alice3_pid::MIDId);
} // namespace o2::aod

struct HfCandidateSelectorJpsiAlice3PidIndexBuilder {
  Builds<o2::aod::HfTrackIndexALICE3PID> index;
  void init(InitContext&) {}
};

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec isAlice3{"isAlice3", VariantType::Bool, true, {"Switch between ALICE 2 and ALICE 3 detector setup"}};
  workflowOptions.push_back(isAlice3);
}

#include "Framework/runDataProcessing.h"

/// Struct for applying J/ψ → e+ e−, μ+ μ− selection cuts
struct HfCandidateSelectorJpsi {
  Produces<aod::HfSelJpsi> hfSelJpsiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  Configurable<bool> selectENotPi{"selectENotPi", true, "Apply combined TOF + RICH e/π selection"};
  // TPC
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 10., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  // TOF
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // RICH
  Configurable<double> ptPidRichMin{"ptPidRichMin", 0.15, "Lower bound of track pT for RICH PID"};
  Configurable<double> ptPidRichMax{"ptPidRichMax", 10., "Upper bound of track pT for RICH PID"};
  Configurable<double> nSigmaRichMax{"nSigmaRichMax", 3., "Nsigma cut on RICH only"};
  Configurable<double> nSigmaRichCombinedTofMax{"nSigmaRichCombinedTofMax", 5., "Nsigma cut on RICH combined with TOF"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_jpsi_to_e_e::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_jpsi_to_e_e::cuts[0], hf_cuts_jpsi_to_e_e::nBinsPt, hf_cuts_jpsi_to_e_e::nCutVars, hf_cuts_jpsi_to_e_e::labelsPt, hf_cuts_jpsi_to_e_e::labelsCutVar}, "Jpsi candidate selection per pT bin"};

  HfHelper hfHelper;
  TrackSelectorEl selectorElectron;
  TrackSelectorMu selectorMuon;

  using TracksSelAlice2 = soa::Join<aod::TracksWDca, aod::TracksPidEl>;
  using TracksSelAlice3 = soa::Join<aod::TracksWDca, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::HfTrackIndexALICE3PID>;

  void init(InitContext const&)
  {
    selectorElectron.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorElectron.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
    selectorElectron.setRangePtTof(ptPidTofMin, ptPidTofMax);
    selectorElectron.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
    selectorElectron.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    selectorElectron.setRangePtRich(ptPidRichMin, ptPidRichMax);
    selectorElectron.setRangeNSigmaRich(-nSigmaRichMax, nSigmaRichMax);
    selectorElectron.setRangeNSigmaRichCondTof(-nSigmaRichCombinedTofMax, nSigmaRichCombinedTofMax);
    selectorMuon = selectorElectron;
  }

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \param trackPos is the positive track
  /// \param trackNeg is the negative track
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selectionTopol(const T1& candidate, const T2& trackPos, const T2& trackNeg, int& selEE, int& selMuMu)
  {
    auto candpT = candidate.pt();
    auto pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptCandMin || candpT >= ptCandMax) {
      return false;
    }

    // cut on e+ e− invariant mass
    if (std::abs(hfHelper.invMassJpsiToEE(candidate) - o2::constants::physics::MassJPsi) > cuts->get(pTBin, "m")) {
      selEE = 0;
    }

    // cut on μ+ μ− invariant mass
    if (std::abs(hfHelper.invMassJpsiToMuMu(candidate) - o2::constants::physics::MassJPsi) > cuts->get(pTBin, "m")) {
      selMuMu = 0;
    }

    if (selEE == 0 && selMuMu == 0) {
      return false;
    }

    // cut on daughter pT (same cut used for both channels)
    if (trackNeg.pt() < cuts->get(pTBin, "pT El") || trackPos.pt() < cuts->get(pTBin, "pT El")) {
      return false;
    }

    // cut on daughter DCA - need to add secondary vertex constraint here
    if (std::abs(trackNeg.dcaXY()) > cuts->get(pTBin, "DCA_xy") || std::abs(trackPos.dcaXY()) > cuts->get(pTBin, "DCA_xy")) {
      return false;
    }

    // cut on daughter DCA - need to add secondary vertex constraint here
    if (std::abs(trackNeg.dcaZ()) > cuts->get(pTBin, "DCA_z") || std::abs(trackPos.dcaZ()) > cuts->get(pTBin, "DCA_z")) {
      return false;
    }

    // cut on chi2 point of closest approach
    if (std::abs(candidate.chi2PCA()) > cuts->get(pTBin, "chi2PCA")) {
      return false;
    }
    return true;
  }

  void processAlice2(aod::HfCand2Prong const& candidates,
                     TracksSelAlice2 const&)
  {
    // looping over 2-prong candidates
    for (const auto& candidate : candidates) {

      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::JpsiToEE) && !(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::JpsiToMuMu)) {
        hfSelJpsiCandidate(0, 0, 0, 0, 0, 0, 0, 0, 0);
        // hfSelJpsiCandidate(0, 0);
        continue;
      }

      auto trackPos = candidate.prong0_as<TracksSelAlice2>(); // positive daughter
      auto trackNeg = candidate.prong1_as<TracksSelAlice2>(); // negative daughter

      int selectedEETopol = 1;
      int selectedEETpc = 1;
      int selectedEETof = 1;
      int selectedEERich = 1;
      int selectedEETofRich = 1;
      int selectedMuMuTopol = 1;
      int selectedMuMuMid = 1;
      int selectedEE = 1;
      int selectedMuMu = 1;

      // track selection level need to add special cuts (additional cuts on decay length and d0 norm)

      if (!selectionTopol(candidate, trackPos, trackNeg, selectedEETopol, selectedMuMuTopol)) {
        selectedEETopol = 0;
        selectedMuMuTopol = 0;
        selectedEE = 0;
        selectedMuMu = 0;
        // if (!selectionTopol(candidate, trackPos, trackNeg, selectedEE, selectedMuMu)) {
        // hfSelJpsiCandidate(0, 0);
        // continue;
      }

      // track-level electron PID TOF selection
      if (selectorElectron.statusTof(trackPos) == TrackSelectorPID::Rejected ||
          selectorElectron.statusTof(trackNeg) == TrackSelectorPID::Rejected) {
        selectedEETof = 0;
        selectedEE = 0;
        // if (selectedMuMu == 0) {
        //   hfSelJpsiCandidate(0, 0);
        //   continue;
        // }
      }

      // track-level electron PID TPC selection
      if (selectorElectron.statusTpc(trackPos) == TrackSelectorPID::Rejected ||
          selectorElectron.statusTpc(trackNeg) == TrackSelectorPID::Rejected) {
        selectedEETpc = 0;
        selectedEE = 0;
      }

      hfSelJpsiCandidate(selectedEE,
                         selectedMuMu,
                         selectedEETopol,
                         selectedEETpc,
                         selectedEETof,
                         selectedEERich,
                         selectedEETofRich,
                         selectedMuMuTopol,
                         selectedMuMuMid);
      // hfSelJpsiCandidate(selectedEE, selectedMuMu);
    }
  }

  PROCESS_SWITCH(HfCandidateSelectorJpsi, processAlice2, "Use ALICE 2 detector setup", true);

  void processAlice3(aod::HfCand2Prong const& candidates,
                     TracksSelAlice3 const&,
                     aod::RICHs const&,
                     aod::MIDs const&)
  {
    // looping over 2-prong candidates
    for (const auto& candidate : candidates) {

      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::JpsiToEE) && !(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::JpsiToMuMu)) {
        hfSelJpsiCandidate(0, 0, 0, 0, 0, 0, 0, 0, 0);
        // hfSelJpsiCandidate(0, 0);
        continue;
      }

      auto trackPos = candidate.prong0_as<TracksSelAlice3>(); // positive daughter
      auto trackNeg = candidate.prong1_as<TracksSelAlice3>(); // negative daughter

      int selectedEETopol = 1;
      int selectedEETpc = 1;
      int selectedEETof = 1;
      int selectedEERich = 1;
      int selectedEETofRich = 1;
      int selectedMuMuTopol = 1;
      int selectedMuMuMid = 1;
      int selectedEE = 1;
      int selectedMuMu = 1;

      // track selection level need to add special cuts (additional cuts on decay length and d0 norm)

      if (!selectionTopol(candidate, trackPos, trackNeg, selectedEETopol, selectedMuMuTopol)) {
        selectedEETopol = 0;
        selectedMuMuTopol = 0;
        selectedEE = 0;
        selectedMuMu = 0;
        // if (!selectionTopol(candidate, trackPos, trackNeg, selectedEE, selectedMuMu)) {
        // hfSelJpsiCandidate(0, 0);
        // continue;
      }

      // if (selectENotPi) {
      //  combined TOF + RICH e selection with π rejection
      if (!selectorElectron.isElectronAndNotPion(trackPos) ||
          !selectorElectron.isElectronAndNotPion(trackNeg)) {
        selectedEETofRich = 0;
        selectedEE = 0;
      }
      //} else {
      // track-level electron PID TOF selection
      if (selectorElectron.statusTof(trackPos) == TrackSelectorPID::Rejected ||
          selectorElectron.statusTof(trackNeg) == TrackSelectorPID::Rejected) {
        selectedEETof = 0;
        selectedEE = 0;
        // if (selectedMuMu == 0) {
        //   hfSelJpsiCandidate(0, 0);
        //   continue;
        // }
      }

      // track-level electron PID RICH selection
      if (selectorElectron.statusRich(trackPos) == TrackSelectorPID::Rejected ||
          selectorElectron.statusRich(trackNeg) == TrackSelectorPID::Rejected) {
        selectedEERich = 0;
        selectedEE = 0;
      }
      //}

      // if (selectedEE == 0 && selectedMuMu == 0) {
      //   hfSelJpsiCandidate(0, 0);
      //   continue;
      // }

      // track-level muon PID MID selection
      if (selectorMuon.statusMid(trackPos) != TrackSelectorPID::Accepted ||
          selectorMuon.statusMid(trackNeg) != TrackSelectorPID::Accepted) {
        selectedMuMuMid = 0;
        selectedMuMu = 0;
      }

      hfSelJpsiCandidate(selectedEE,
                         selectedMuMu,
                         selectedEETopol,
                         selectedEETpc,
                         selectedEETof,
                         selectedEERich,
                         selectedEETofRich,
                         selectedMuMuTopol,
                         selectedMuMuMid);
      // hfSelJpsiCandidate(selectedEE, selectedMuMu);
    }
  }

  PROCESS_SWITCH(HfCandidateSelectorJpsi, processAlice3, "Use ALICE 3 detector setup", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  const bool isAlice3 = cfgc.options().get<bool>("isAlice3");
  if (isAlice3) {
    workflow.push_back(adaptAnalysisTask<HfCandidateSelectorJpsiAlice3PidIndexBuilder>(cfgc));
    workflow.push_back(adaptAnalysisTask<HfCandidateSelectorJpsi>(cfgc, SetDefaultProcesses{{{"processAlice2", false}, {"processAlice3", true}}}));
  } else {
    workflow.push_back(adaptAnalysisTask<HfCandidateSelectorJpsi>(cfgc));
  }
  return workflow;
}
