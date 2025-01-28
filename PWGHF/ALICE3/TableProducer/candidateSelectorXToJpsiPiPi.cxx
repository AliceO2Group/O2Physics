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

/// \file candidateSelectorXToJpsiPiPi.cxx
/// \brief X(3872) selection task.
/// \note Adapted from candidateSelectorJpsi.cxx
///
/// \author Rik Spijkers <r.spijkers@students.uu.nl>, Utrecht University
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

/// Struct for applying Jpsi selection cuts
struct HfCandidateSelectorXToJpsiPiPi {
  Produces<aod::HfSelXToJpsiPiPi> hfSelXToJpsiPiPiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 10., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  // Configurable<double> TPCNClsFindableMin{"TPCNClsFindableMin", 70., "Lower bound of TPC findable clusters for good PID"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 10., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_x_to_jpsi_pi_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_x_to_jpsi_pi_pi::Cuts[0], hf_cuts_x_to_jpsi_pi_pi::NBinsPt, hf_cuts_x_to_jpsi_pi_pi::NCutVars, hf_cuts_x_to_jpsi_pi_pi::labelsPt, hf_cuts_x_to_jpsi_pi_pi::labelsCutVar}, "Jpsi candidate selection per pT bin"};

  HfHelper hfHelper;

  using TracksSel = soa::Join<aod::Tracks, aod::TracksPidPi>;

  /// Selection on goodness of daughter tracks
  /// \note should be applied at candidate selection
  /// \param track is daughter track
  /// \return true if track is good
  template <typename T>
  bool daughterSelection(const T& /*track*/)
  {
    return true;
  }

  /// Conjugate independent toplogical cuts
  /// \param hfCandX is candidate
  /// \param trackNeutral is the track with the pi+ hypothesis
  /// \param trackPos is the track with the pi+ hypothesis
  /// \param trackNeg is the track with the pi- hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2, typename T3>
  bool selectionTopol(const T1& hfCandX, const T2& hfCandJpsi, const T3& trackPos, const T3& trackNeg)
  {
    auto candpT = hfCandX.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      // LOGF(debug, "X topol selection failed at getpTBin");
      return false;
    }

    if (candpT < ptCandMin || candpT >= ptCandMax) {
      // LOGF(debug, "X topol selection failed at cand pT check");
      return false; // check that the candidate pT is within the analysis range
    }

    if (std::abs(hfHelper.invMassXToJpsiPiPi(hfCandX) - o2::constants::physics::MassX3872) > cuts->get(pTBin, "m")) {
      // LOGF(debug, "X topol selection failed at mass diff check");
      return false; // check that mass difference is within bounds
    }

    if ((hfCandJpsi.pt() < cuts->get(pTBin, "pT Jpsi")) || (trackNeg.pt() < cuts->get(pTBin, "pT Pi")) || (trackPos.pt() < cuts->get(pTBin, "pT Pi"))) {
      // LOGF(debug, "X topol selection failed at daughter pT check");
      return false; // cut on daughter pT
    }

    if (hfCandX.cpa() < cuts->get(pTBin, "CPA")) {
      return false; // CPA check
    }

    if ((std::abs(hfCandX.impactParameter0()) > cuts->get(pTBin, "d0 Jpsi")) ||
        (std::abs(hfCandX.impactParameter1()) > cuts->get(pTBin, "d0 Pi")) ||
        (std::abs(hfCandX.impactParameter2()) > cuts->get(pTBin, "d0 Pi"))) {
      return false; // DCA check on daughters
    }

    // add more cuts: d0 product? PCA?

    return true;
  }

  /// Check if track is ok for TPC PID
  /// \param track is the track
  /// \note function to be expanded
  /// \return true if track is ok for TPC PID
  template <typename T>
  bool validTPCPID(const T& track)
  {
    if (std::abs(track.pt()) < ptPidTpcMin || std::abs(track.pt()) >= ptPidTpcMax) {
      return false;
    }
    // if (track.TPCNClsFindable() < TPCNClsFindableMin) return false;
    return true;
  }

  /// Check if track is ok for TOF PID
  /// \param track is the track
  /// \note function to be expanded
  /// \return true if track is ok for TOF PID
  template <typename T>
  bool validTofPid(const T& track)
  {
    if (std::abs(track.pt()) < ptPidTofMin || std::abs(track.pt()) >= ptPidTofMax) {
      return false;
    }
    return true;
  }

  /// Check if track is compatible with given TPC Nsigma cut for the pion hypothesis
  /// \param track is the track
  /// \param nSigmaCut is the nsigma threshold to test against
  /// \return true if track satisfies TPC pion hypothesis for given Nsigma cut
  template <typename T>
  bool selectionPIDTPC(const T& track, int nSigmaCut)
  {
    if (nSigmaCut > 999.) {
      return true;
    }
    return track.tpcNSigmaPi() < nSigmaCut;
  }

  /// Check if track is compatible with given TOF NSigma cut for the pion hypothesis
  /// \param track is the track
  /// \param nSigmaCut is the nSigma threshold to test against
  /// \return true if track satisfies TOF pion hypothesis for given NSigma cut
  template <typename T>
  bool selectionPIDTOF(const T& track, double nSigmaCut)
  {
    if (nSigmaCut > 999.) {
      return true;
    }
    return track.tofNSigmaPi() < nSigmaCut;
  }

  /// PID selection on daughter track
  /// \param track is the daughter track
  /// \return 1 if successful PID match, 0 if successful PID rejection, -1 if no PID info
  template <typename T>
  int selectionPID(const T& track)
  { // use both TPC and TOF here; in run5 only TOF makes sense. add some flag for run3/run5 data later?
    if (validTofPid(track)) {
      if (!selectionPIDTOF(track, nSigmaTofMax)) {
        return 0; // rejected by PID
      } else {
        return 1; // positive PID
      }
    } else {
      return -1; // no PID info
    }

    // no tpc in run5, so for now comment it out
    // if (validTPCPID(track)) {
    //   if (!selectionPIDTPC(track, nSigmaTpcMax)) {
    //     return 0; //rejected by PID
    //   } else {
    //     return 1; //positive PID
    //   }
    // } else {
    //   return -1; //no PID info
    // }
  }

  void process(aod::HfCandX const& hfCandXs,
               aod::HfCand2Prong const&,
               TracksSel const&)
  {
    for (const auto& hfCandX : hfCandXs) { // looping over X candidates
      // note the difference between Jpsi (index0) and pions (index1,2)
      auto candJpsi = hfCandX.prong0();
      auto trackPos = hfCandX.prong1_as<TracksSel>(); // positive daughter
      auto trackNeg = hfCandX.prong2_as<TracksSel>(); // negative daughter

      int selJpsiToEE = 1;
      int selJpsiToMuMu = 1;

      // check if flagged as X --> Jpsi Pi Pi
      if (!(hfCandX.hfflag() & 1 << hf_cand_x::DecayType::XToJpsiToEEPiPi)) {
        selJpsiToEE = 0;
      }

      if (!(hfCandX.hfflag() & 1 << hf_cand_x::DecayType::XToJpsiToMuMuPiPi)) {
        selJpsiToMuMu = 0;
      }

      if (selJpsiToEE == 0 && selJpsiToMuMu == 0) {
        hfSelXToJpsiPiPiCandidate(0, 0);
        continue;
      }

      // daughter track validity selection
      if (!daughterSelection(trackPos) || !daughterSelection(trackNeg)) {
        hfSelXToJpsiPiPiCandidate(0, 0);
        // LOGF(debug, "X candidate selection failed at daughter selection");
        continue;
      }

      // implement filter bit 4 cut - should be done before this task at the track selection level
      // need to add special cuts (additional cuts on decay length and d0 norm)

      if (!selectionTopol(hfCandX, candJpsi, trackPos, trackNeg)) {
        hfSelXToJpsiPiPiCandidate(0, 0);
        // LOGF(debug, "X candidate selection failed at selection topology");
        continue;
      }

      if (selectionPID(trackPos) == 0 || selectionPID(trackNeg) == 0) {
        hfSelXToJpsiPiPiCandidate(0, 0);
        // LOGF(debug, "X candidate selection failed at selection PID");
        continue;
      }

      hfSelXToJpsiPiPiCandidate(selJpsiToEE, selJpsiToMuMu);
      // LOGF(debug, "X candidate selection successful, candidate should be selected");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorXToJpsiPiPi>(cfgc)};
}
