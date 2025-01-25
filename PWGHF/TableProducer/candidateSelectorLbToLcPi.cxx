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

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

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

  Configurable<float> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // TPC PID
  Configurable<float> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<float> ptPidTpcMax{"ptPidTpcMax", 10., "Upper bound of track pT for TPC PID"};
  Configurable<float> nSigmaTpcMax{"nSigmaTpcMax", 5., "Nsigma cut on TPC only"};
  Configurable<float> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<float> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<float> ptPidTofMax{"ptPidTofMax", 10., "Upper bound of track pT for TOF PID"};
  Configurable<float> nSigmaTofMax{"nSigmaTofMax", 5., "Nsigma cut on TOF only"};
  Configurable<float> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lb_to_lc_pi::vecBinsPt}, "pT bin limits"};
  Configurable<float> impactParameterMaximum{"impactParameterMaximum", 0.2, "Maximum impact parameter for single tracks"};
  Configurable<float> maxDecayLengthError{"maxDecayLengthError", 0.015, "decay length error quality selection"};
  Configurable<float> maxDecayLengthXYError{"maxDecayLengthXYError", 0.01, "decay length xy error quality selection"};
  Configurable<float> maxVertexDistanceLbLc{"maxVertexDistanceLbLc", 0.05, "maximum distance between Lb and Lc vertex"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lb_to_lc_pi::cuts[0], hf_cuts_lb_to_lc_pi::nBinsPt, hf_cuts_lb_to_lc_pi::nCutVars, hf_cuts_lb_to_lc_pi::labelsPt, hf_cuts_lb_to_lc_pi::labelsCutVar}, "Lb0 candidate selection per pT bin"};
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc+"};

  HfHelper hfHelper;

  using TracksWExt = soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, aod::TrackSelection, o2::aod::TrackSelectionExtension, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::pidTPCFullPi>;

  bool passesImpactParameterResolution(float pT, float d0Resolution)
  {
    float expectedResolution(0.001 + 0.0052 * exp(-0.655 * pT));
    if (d0Resolution > expectedResolution * 1.5)
      return false;
    else
      return true;
  } // Compares to pT dependent cut on impact parameter resolution

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

    // Λb0 mass cut
    if (std::abs(hfHelper.invMassLbToLcPi(hfCandLb) - o2::constants::physics::MassLambdaB0) > cuts->get(pTBin, "m")) {
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

    float lcMass = 0.;
    if (hfCandLc.isSelLcToPKPi())
      lcMass = hfHelper.invMassLcToPKPi(hfCandLc);
    if (hfCandLc.isSelLcToPiKP())
      lcMass = hfHelper.invMassLcToPiKP(hfCandLc);
    if (std::abs(lcMass - o2::constants::physics::MassLambdaCPlus) > cuts->get(pTBin, "DeltaMLc"))
      return false;

    if (hfCandLb.errorDecayLengthXY() > maxDecayLengthXYError)
      return false;
    if (hfCandLb.errorDecayLength() > maxDecayLengthError)
      return false;

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

    // distance between Lb and Lc decay
    float diffXVert = hfCandLb.xSecondaryVertex() - hfCandLc.xSecondaryVertex();
    float diffYVert = hfCandLb.ySecondaryVertex() - hfCandLc.ySecondaryVertex();
    float diffZVert = hfCandLb.zSecondaryVertex() - hfCandLc.zSecondaryVertex();
    float vertexDistance = sqrt(diffXVert * diffXVert + diffYVert * diffYVert + diffZVert * diffZVert);
    if (vertexDistance > maxVertexDistanceLbLc) {
      return false;
    }

    return true;
  }

  void process(aod::HfCandLb const& hfCandLbs,
               soa::Join<aod::HfCand3Prong, aod::HfSelLc> const&,
               TracksWExt const&)
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
      auto trackPi = hfCandLb.prong1_as<TracksWExt>();

      // Check that the impact parameter resolution is not too far from the typical
      auto track0 = candLc.prong0_as<TracksWExt>();
      auto track1 = candLc.prong1_as<TracksWExt>();
      auto track2 = candLc.prong2_as<TracksWExt>();
      float reso0 = candLc.errorImpactParameter0();
      float reso1 = candLc.errorImpactParameter1();
      float reso2 = candLc.errorImpactParameter2();
      if (!passesImpactParameterResolution(track0.pt(), reso0) || !passesImpactParameterResolution(track1.pt(), reso1) || !passesImpactParameterResolution(track2.pt(), reso2) || !passesImpactParameterResolution(trackPi.pt(), hfCandLb.errorImpactParameter1())) {
        hfSelLbToLcPiCandidate(statusLb);
        continue;
      }

      // Maximum single-track impact parameter selection to suppress strange background
      if (std::abs(hfCandLb.impactParameter1()) > impactParameterMaximum || candLc.impactParameter0() > impactParameterMaximum || candLc.impactParameter1() > impactParameterMaximum || candLc.impactParameter2() > impactParameterMaximum) {
        hfSelLbToLcPiCandidate(statusLb);
        continue;
      }

      // topological cuts
      if (!selectionTopol(hfCandLb, candLc, trackPi)) {
        hfSelLbToLcPiCandidate(statusLb);
        // LOGF(debug, "Lb candidate selection failed at selection topology");
        continue;
      }

      // PID selection for pion
      if (trackPi.pt() > ptPidTpcMin && trackPi.pt() < ptPidTpcMax) {
        if (std::abs(trackPi.tpcNSigmaPi()) > nSigmaTpcMax) {
          hfSelLbToLcPiCandidate(statusLb);
          continue;
        }
      }
      if (trackPi.pt() > ptPidTofMin && trackPi.pt() < ptPidTofMax) {
        if (std::abs(trackPi.tofNSigmaPi()) > nSigmaTofMax) {
          hfSelLbToLcPiCandidate(statusLb);
          continue;
        }
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
