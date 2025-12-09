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

/// \file candidateSelectorChicToJpsiGamma.cxx
/// \brief chi_c selection task.
/// \note Adapted from candidateSelectorXToJpsiPiPi.cxx
///
/// \author Alessandro De Falco <alessandro.de.falco@ca.infn.it>, Università/INFN Cagliari

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

/// Struct for applying Jpsi selection cuts
struct HfCandidateSelectorChicToJpsiGamma {
  Produces<aod::HfSelChicToJpsiGamma> hfSelChicToJpsiGammaCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // Configurable<double> TPCNClsFindableMin{"TPCNClsFindableMin", 70., "Lower bound of TPC findable clusters for good PID"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 10., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  // TPC PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 10., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_chic_to_jpsi_gamma::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_chic_to_jpsi_gamma::Cuts[0], hf_cuts_chic_to_jpsi_gamma::NBinsPt, hf_cuts_chic_to_jpsi_gamma::NCutVars, hf_cuts_chic_to_jpsi_gamma::labelsPt, hf_cuts_chic_to_jpsi_gamma::labelsCutVar}, "Jpsi candidate selection per pT bin"};

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
  /// \param hfCandChic is candidate
  /// \param trackNeutral is the track with the pi+ hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2, typename T3>
  bool selectionTopol(const T1& hfCandChic, const T2& hfCandJpsi, const T3& /*ecal*/)
  {
    auto candpT = hfCandChic.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    if (candpT < ptCandMin || candpT >= ptCandMax) {
      return false; // check that the candidate pT is within the analysis range
    }

    auto mchic = o2::constants::physics::MassChiC1; // chi_c1(1p)
    if (std::abs(HfHelper::invMassChicToJpsiGamma(hfCandChic) - mchic) > cuts->get(pTBin, "m")) {
      // LOGF(debug, "Chic topol selection failed at mass diff check");
      return false; // check that mass difference is within bounds
    }

    if ((hfCandJpsi.pt() < cuts->get(pTBin, "pT Jpsi"))) { // adf: Warning: no cut on photon
      return false;                                        // cut on daughter pT
    }

    // if ((hfCandJpsi.pt() < cuts->get(pTBin, "pT Jpsi")) || (trackNeg.pt() < cuts->get(pTBin, "pT Pi")) || (trackPos.pt() < cuts->get(pTBin, "pT Pi"))) {
    //   return false; //cut on daughter pT
    // }

    if (hfCandChic.cpa() < cuts->get(pTBin, "CPA")) {
      return false; // CPA check
    }

    if ((std::abs(hfCandChic.impactParameter0()) > cuts->get(pTBin, "d0 Jpsi"))) { // adf: Warning: no cut on photon
      return false;                                                                // DCA check on daughters
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

  //------------------------------------------------------------------------------------
  /// Check if track is compatible with given TPC Nsigma cut for the pion hypothesis
  /// \param track is the track
  /// \param nSigmaCut is the nsigma threshold to test against
  /// \return true if track satisfies TPC pion hypothesis for given Nsigma cut
  template <typename T>
  bool selectionPIDTPC(const T& /*track*/, int nSigmaCut)
  {
    if (nSigmaCut > 999.) {
      return true;
    }
    //    return track.tpcNSigmaPi() < nSigmaCut;
    return true;
  }

  // Check if track is compatible with given TOF NSigma cut for the pion hypothesis
  // \param track is the track
  // \param nSigmaCut is the nSigma threshold to test against
  // \return true if track satisfies TOF pion hypothesis for given NSigma cut
  template <typename T>
  bool selectionPIDTOF(const T& /*track*/, double nSigmaCut)
  {
    if (nSigmaCut > 999.) {
      return true;
    }
    //    return track.tofNSigmaPi() < nSigmaCut;
    return true;
  }

  // PID selection on daughter track
  // \param track is the daughter track
  // \return 1 if successful PID match, 0 if successful PID rejection, -1 if no PID info
  template <typename T>
  int selectionPID(const T& /*track*/)
  { // use both TPC and TOF here; in run5 only TOF makes sense. add some flag for run3/run5 data later?
    // if (validTofPid(track)) {
    //   if (!selectionPIDTOF(track, nSigmaTofMax)) {
    //     return 0; //rejected by PID
    //   } else {
    //     return 1; //positive PID
    //   }
    // } else {
    //   return -1; //no PID info

    return true;
  }

  //---------------------------------------------------------------

  void process(aod::HfCandChic const& hfCandChics,
               aod::HfCand2Prong const&,
               aod::ECALs const&)
  {
    for (const auto& hfCandChic : hfCandChics) { // looping over chi_c candidates
      // note the difference between Jpsi (index0) and pions (index1,2)
      auto candJpsi = hfCandChic.prong0();
      auto gamma = hfCandChic.prong1_as<aod::ECALs>();

      int selJpsiToEE = 1;
      int selJpsiToMuMu = 1;

      // check if flagged as chic --> Jpsi gamma
      if (!(hfCandChic.hfflag() & 1 << hf_cand_chic::DecayType::ChicToJpsiToEEGamma)) {
        selJpsiToEE = 0;
      }

      if (!(hfCandChic.hfflag() & 1 << hf_cand_chic::DecayType::ChicToJpsiToMuMuGamma)) {
        selJpsiToMuMu = 0;
      }

      if (selJpsiToEE == 0 && selJpsiToMuMu == 0) {
        hfSelChicToJpsiGammaCandidate(0, 0);
        continue;
      }

      // daughter track validity selection
      if (!daughterSelection(gamma)) { // no selection at present
        hfSelChicToJpsiGammaCandidate(0, 0);
        continue;
      }

      // implement filter bit 4 cut - should be done before this task at the track selection level
      // need to add special cuts (additional cuts on decay length and d0 norm)

      if (!selectionTopol(hfCandChic, candJpsi, gamma)) { // check selections
        hfSelChicToJpsiGammaCandidate(0, 0);
        continue;
      }

      if (selectionPID(gamma) == 0) { // no selection at present
        hfSelChicToJpsiGammaCandidate(0, 0);
        continue;
      }

      hfSelChicToJpsiGammaCandidate(selJpsiToEE, selJpsiToMuMu);
      // LOGF(debug, "Chi_c candidate selection successful, candidate should be selected");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorChicToJpsiGamma>(cfgc)};
}
