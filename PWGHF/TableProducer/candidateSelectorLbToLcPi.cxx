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

/// \file candidateSelectorLbToLcPi.cxx
/// \brief Λb0 → Λc+ π- candidate selector
///
/// \author Panos Christakoglou <panos.christakoglou@cern.ch>, Nikhef

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

struct HfCandidateSelectorLbToLcPi {
  Produces<aod::HfSelLbToLcPi> hfSelLbToLcPiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 10., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 5., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 10., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 5., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lb_to_lc_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lb_to_lc_pi::cuts[0], hf_cuts_lb_to_lc_pi::nBinsPt, hf_cuts_lb_to_lc_pi::nCutVars, hf_cuts_lb_to_lc_pi::labelsPt, hf_cuts_lb_to_lc_pi::labelsCutVar}, "Lb0 candidate selection per pT bin"};
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc+"};

  HfHelper hfHelper;

  // Apply topological cuts as defined in SelectorCuts.h; return true if candidate passes all cuts
  template <typename T1, typename T2, typename T3>
  bool selectionTopol(const T1& hfCandLb, const T2& hfCandLc, const T3& trackPi)
  {
    auto candpT = hfCandLb.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      // LOGF(debug, "Lb topol selection failed at getpTBin");
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptCandMin || candpT >= ptCandMax) {
      return false;
    }

    //Λb0 mass cut
    if (std::abs(hfHelper.invMassLbToLcPi(hfCandLb) - o2::analysis::pdg::MassLambdaB0) > cuts->get(pTBin, "m")) {
      // LOGF(debug, "Lb topol selection failed at mass diff check");
      return false;
    }

    // pion pt
    if (trackPi.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    // Lc+ pt
    if (hfCandLc.pt() < cuts->get(pTBin, "pT Lc+")) {
      return false;
    }

    // Lc mass
    // if (trackPi.sign() < 0) {
    // if (std::abs(hfHelper.invMassLcToPKPi(hfCandLc) - o2::analysis::pdg::MassLambdaCPlus) > cuts->get(pTBin, "DeltaMLc")) {
    // return false;
    // }
    // }

    // Lb Decay length
    if (hfCandLb.decayLength() < cuts->get(pTBin, "Lb decLen")) {
      return false;
    }

    // Lb Decay length XY
    if (hfCandLb.decayLengthXY() < cuts->get(pTBin, "Lb decLenXY")) {
      return false;
    }

    // Lb chi2PCA cut
    if (hfCandLb.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
      // LOGF(debug, "Lb selection failed at chi2PCA");
      return false;
    }

    // Lb CPA cut
    if (hfCandLb.cpa() < cuts->get(pTBin, "CPA")) {
      return false;
    }

    // d0 of pi
    if (std::abs(hfCandLb.impactParameter1()) < cuts->get(pTBin, "d0 Pi")) {
      return false;
    }

    // d0 of Lc+
    if (std::abs(hfCandLb.impactParameter0()) < cuts->get(pTBin, "d0 Lc+")) {
      return false;
    }

    return true;
  }

  void process(aod::HfCandLb const& hfCandLbs,
               soa::Join<aod::HfCand3Prong, aod::HfSelLc> const&,
               aod::Tracks const&)
  {
    for (const auto& hfCandLb : hfCandLbs) { // looping over Lb candidates

      int statusLb = 0;

      // check if flagged as Λb --> Λc+ π-
      if (!(hfCandLb.hfflag() & 1 << hf_cand_lb::DecayType::LbToLcPi)) {
        hfSelLbToLcPiCandidate(statusLb);
        // LOGF(debug, "Lb candidate selection failed at hfflag check");
        continue;
      }

      // Lc is always index0 and pi is index1 by default
      // auto candLc = hfCandLb.prong0();
      auto candLc = hfCandLb.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>();
      auto trackPi = hfCandLb.prong1();

      // topological cuts
      if (!selectionTopol(hfCandLb, candLc, trackPi)) {
        hfSelLbToLcPiCandidate(statusLb);
        // LOGF(debug, "Lb candidate selection failed at selection topology");
        continue;
      }

      hfSelLbToLcPiCandidate(1);
      // LOGF(debug, "Lb candidate selection successful, candidate should be selected");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfCandidateSelectorLbToLcPi>(cfgc));
  return workflow;
}
