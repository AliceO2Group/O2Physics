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

/// \file candidateSelectorD0Alice3Forward.cxx
/// \brief D0(bar) → π± K∓ selection task
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include <vector>
#include <algorithm>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "ALICE3/DataModel/RICH.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

namespace indices
{
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_INDEX_COLUMN(FRICH, frich);
} // namespace indices
DECLARE_SOA_INDEX_TABLE_USER(FRICHTracksIndex, Tracks, "FRICHTRK", indices::TrackId, indices::FRICHId);
} // namespace o2::aod

struct HfCandidateSelectorD0Alice3ForwardRichIndexBuilder { // Builder of the RICH-track index linkage
  Builds<o2::aod::FRICHTracksIndex> indF;

  void init(InitContext&) {}
};

/// Struct for applying D0 selection cuts
struct HfCandidateSelectorD0Alice3Forward {
  Produces<aod::HfSelD0Alice3Forward> hfSelD0CandidateALICE3Forward;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_d0_to_pi_k::cuts[0], hf_cuts_d0_to_pi_k::nBinsPt, hf_cuts_d0_to_pi_k::nCutVars, hf_cuts_d0_to_pi_k::labelsPt, hf_cuts_d0_to_pi_k::labelsCutVar}, "D0 candidate selection per pT bin"};

  HfHelper hfHelper;

  using TracksSel = soa::Join<aod::TracksWDca, aod::FRICHTracksIndex>;

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T>
  bool selectionTopol(const T& candidate)
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
    // product of daughter impact parameters
    if (candidate.impactParameterProduct() > cuts->get(pTBin, "d0d0")) {
      return false;
    }
    // cosine of pointing angle
    if (candidate.cpa() < cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }
    // cosine of pointing angle XY
    if (candidate.cpaXY() < cuts->get(pTBin, "cos pointing angle xy")) {
      return false;
    }
    // normalised decay length in XY plane
    if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "normalized decay length XY")) {
      return false;
    }
    // candidate DCA
    // if (candidate.chi2PCA() > cuts[pTBin][1]) return false;

    // decay exponentail law, with tau = beta*gamma*ctau
    // decay length > ctau retains (1-1/e)
    if (std::abs(candidate.impactParameterNormalised0()) < 0.5 || std::abs(candidate.impactParameterNormalised1()) < 0.5) {
      return false;
    }
    double decayLengthCut = std::min((candidate.p() * 0.0066) + 0.01, cuts->get(pTBin, "minimum decay length"));
    if (candidate.decayLength() * candidate.decayLength() < decayLengthCut * decayLengthCut) {
      return false;
    }
    if (candidate.decayLength() > cuts->get(pTBin, "decay length")) {
      return false;
    }
    if (candidate.decayLengthXY() > cuts->get(pTBin, "decay length XY")) {
      return false;
    }
    if (candidate.decayLengthNormalised() * candidate.decayLengthNormalised() < 1.0) {
      // return false; // add back when getter fixed
    }
    return true;
  }

  /// Conjugate-dependent topological cuts
  /// \param candidate is candidate
  /// \param trackPion is the track with the pion hypothesis
  /// \param trackKaon is the track with the kaon hypothesis
  /// \note trackPion = positive and trackKaon = negative for D0 selection and inverse for D0bar
  /// \return true if candidate passes all cuts for the given Conjugate
  template <typename T1, typename T2>
  bool selectionTopolConjugate(const T1& candidate, const T2& trackPion, const T2& trackKaon)
  {
    auto candpT = candidate.pt();
    auto pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // invariant-mass cut
    if (trackPion.sign() > 0) {
      if (std::abs(hfHelper.invMassD0ToPiK(candidate) - o2::constants::physics::MassD0) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(hfHelper.invMassD0barToKPi(candidate) - o2::constants::physics::MassD0) > cuts->get(pTBin, "m")) {
        return false;
      }
    }

    // cut on daughter pT
    if (trackPion.pt() < cuts->get(pTBin, "pT Pi") || trackKaon.pt() < cuts->get(pTBin, "pT K")) {
      return false;
    }

    // cut on daughter DCA - need to add secondary vertex constraint here
    if (std::abs(trackPion.dcaXY()) > cuts->get(pTBin, "d0pi") || std::abs(trackKaon.dcaXY()) > cuts->get(pTBin, "d0K")) {
      return false;
    }

    // cut on cos(theta*)
    if (trackPion.sign() > 0) {
      if (std::abs(hfHelper.cosThetaStarD0(candidate)) > cuts->get(pTBin, "cos theta*")) {
        return false;
      }
    } else {
      if (std::abs(hfHelper.cosThetaStarD0bar(candidate)) > cuts->get(pTBin, "cos theta*")) {
        return false;
      }
    }

    return true;
  }

  void process(aod::HfCand2Prong const& candidates,
               TracksSel const&,
               aod::McParticles const&,
               aod::RICHs const&,
               aod::FRICHs const&)
  {

    for (const auto& candidate : candidates) {

      // selection flag
      int statusHFFlag = 0;
      int statusD0NoPid = 0;
      int statusD0RICHPID = 0;

      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        hfSelD0CandidateALICE3Forward(statusHFFlag, statusD0NoPid, statusD0RICHPID);
        continue;
      }
      statusHFFlag = 1;

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelD0CandidateALICE3Forward(statusHFFlag, statusD0NoPid, statusD0RICHPID);
        continue;
      }

      auto trackPos = candidate.prong0_as<TracksSel>();
      auto trackNeg = candidate.prong1_as<TracksSel>();

      // auto momentumPosTrack = trackPos.p();
      // auto momentumNegTrack = trackNeg.p();

      bool topolD0 = selectionTopolConjugate(candidate, trackPos, trackNeg);
      bool topolD0bar = selectionTopolConjugate(candidate, trackNeg, trackPos);

      if (!topolD0 && !topolD0bar) {
        hfSelD0CandidateALICE3Forward(statusHFFlag, statusD0NoPid, statusD0RICHPID);
        continue;
      }

      // float nsigmaTOFNegKaon = -5000.0;
      float nsigmaRICHNegKaon = -5000.0;
      // float nsigmaTOFPosPion = -5000.0;
      float nsigmaRICHPosPion = -5000.0;

      if (trackPos.has_frich()) {
        nsigmaRICHPosPion = trackPos.frich().frichNsigmaPi();
      }
      if (trackNeg.has_frich()) {
        nsigmaRICHNegKaon = trackNeg.frich().frichNsigmaKa();
      }

      // bool selectPionTOFplusRICH = false;
      // bool selectKaonTOFplusRICH = false;

      if (topolD0) {
        statusD0NoPid = 1;
        if ((std::abs(nsigmaRICHPosPion) < 3.0 && std::abs(nsigmaRICHNegKaon) < 3.0))
          statusD0RICHPID = 1;
      }
      hfSelD0CandidateALICE3Forward(statusHFFlag, statusD0NoPid, statusD0RICHPID);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfCandidateSelectorD0Alice3ForwardRichIndexBuilder>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfCandidateSelectorD0Alice3Forward>(cfgc));
  return workflow;
}
