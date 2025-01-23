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

/// \file candidateSelectorD0ParametrizedPid.cxx
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
DECLARE_SOA_INDEX_COLUMN(RICH, rich);
} // namespace indices
DECLARE_SOA_INDEX_TABLE_USER(RICHTracksIndex, Tracks, "RICHTRK", indices::TrackId, indices::RICHId);
} // namespace o2::aod

struct HfCandidateSelectorD0ParametrizedPidRichIndexBuilder { // Builder of the RICH-track index linkage
  Builds<o2::aod::RICHTracksIndex> indB;
  void init(InitContext&) {}
};

/// Struct for applying D0 selection cuts
struct HfCandidateSelectorD0ParametrizedPid {
  Produces<aod::HfSelD0ParametrizedPid> hfSelD0CandidateparametrizedPID;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  Configurable<double> etaPerfectPidMax{"etaPerfectPidMax", 1.75, "Eta cut for perfect PID"};
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

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::RICHTracksIndex, aod::McTrackLabels>;

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
      int statusD0NoPid = 0;
      int statusD0PerfectPid = 0;
      int statusD0 = 0;
      int statusD0barNoPid = 0;
      int statusD0barPerfectPid = 0;
      int statusD0bar = 0;

      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        hfSelD0CandidateparametrizedPID(statusD0NoPid, statusD0PerfectPid, statusD0, statusD0barNoPid, statusD0barPerfectPid, statusD0bar);
        continue;
      }

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelD0CandidateparametrizedPID(statusD0NoPid, statusD0PerfectPid, statusD0, statusD0barNoPid, statusD0barPerfectPid, statusD0bar);
        continue;
      }

      auto trackPos = candidate.prong0_as<TracksSel>();
      auto trackNeg = candidate.prong1_as<TracksSel>();

      // auto momentumPosTrack = trackPos.p();
      // auto momentumNegTrack = trackNeg.p();
      auto ptPosTrack = trackPos.pt();
      auto ptNegTrack = trackNeg.pt();
      auto etaPosTrack = std::abs(trackPos.eta());
      auto etaNegTrack = std::abs(trackNeg.eta());

      bool topolD0 = selectionTopolConjugate(candidate, trackPos, trackNeg);
      bool topolD0bar = selectionTopolConjugate(candidate, trackNeg, trackPos);

      if (!topolD0 && !topolD0bar) {
        hfSelD0CandidateparametrizedPID(statusD0NoPid, statusD0PerfectPid, statusD0, statusD0barNoPid, statusD0barPerfectPid, statusD0bar);
        continue;
      }

      int pdgPositive = 0;
      int pdgNegative = 0;
      if (trackPos.has_mcParticle()) {
        pdgPositive = trackPos.mcParticle().pdgCode();
      }
      if (trackNeg.has_mcParticle()) {
        pdgNegative = trackNeg.mcParticle().pdgCode();
      }

      bool selectPosPion = false;
      bool selectNegKaon = false;
      bool selectNegPion = false;
      bool selectPosKaon = false;

      if (etaPosTrack >= etaPerfectPidMax) {
        if (ptPosTrack < (19.58 / std::cosh(etaPosTrack))) {
          if (pdgPositive == 321)
            selectPosKaon = true;
          if (ptPosTrack < (11.65 / std::cosh(etaPosTrack)) && pdgPositive == 211)
            selectPosPion = true;
        } else {
          selectPosPion = true;
          selectPosKaon = true;
        }
      } else {
        if (trackPos.hasTOF()) {
          if (std::abs(trackPos.tofNSigmaPi()) < 3.0) {
            selectPosPion = true;
          }
          if (std::abs(trackPos.tofNSigmaKa()) < 3.0) {
            selectPosKaon = true;
          }
        }

        if (trackPos.has_rich() && !trackPos.hasTOF()) {
          if (std::abs(trackPos.rich().richNsigmaPi()) < 3.0) {
            selectPosPion = true;
          }
          if (std::abs(trackPos.rich().richNsigmaKa()) < 3.0) {
            selectPosKaon = true;
          }
        }

        if (trackPos.has_rich() && trackPos.hasTOF()) {
          if ((trackPos.rich().richNsigmaPi() * trackPos.rich().richNsigmaPi() + trackPos.tofNSigmaPi() * trackPos.tofNSigmaPi()) < 9.0) {
            selectPosPion = true;
          }
          if ((trackPos.rich().richNsigmaKa() * trackPos.rich().richNsigmaKa() + trackPos.tofNSigmaKa() * trackPos.tofNSigmaKa()) < 9.0) {
            selectPosKaon = true;
          }
        }
      }

      if (etaNegTrack >= etaPerfectPidMax) {
        if (ptNegTrack < (19.58 / std::cosh(etaNegTrack))) {
          if (pdgNegative == -321)
            selectNegKaon = true;
          if (ptNegTrack < (11.65 / std::cosh(etaNegTrack)) && pdgNegative == -211)
            selectNegPion = true;
        } else {
          selectNegPion = true;
          selectNegKaon = true;
        }
      } else {
        if (trackNeg.hasTOF()) {
          if (std::abs(trackNeg.tofNSigmaPi()) < 3.0) {
            selectNegPion = true;
          }
          if (std::abs(trackNeg.tofNSigmaKa()) < 3.0) {
            selectNegKaon = true;
          }
        }

        if (trackNeg.has_rich() && !trackNeg.hasTOF()) {
          if (std::abs(trackNeg.rich().richNsigmaPi()) < 3.0) {
            selectNegPion = true;
          }
          if (std::abs(trackNeg.rich().richNsigmaKa()) < 3.0) {
            selectNegKaon = true;
          }
        }

        if (trackNeg.has_rich() && trackNeg.hasTOF()) {
          if ((trackNeg.rich().richNsigmaPi() * trackNeg.rich().richNsigmaPi() + trackNeg.tofNSigmaPi() * trackNeg.tofNSigmaPi()) < 9.0) {
            selectNegPion = true;
          }
          if ((trackNeg.rich().richNsigmaKa() * trackNeg.rich().richNsigmaKa() + trackNeg.tofNSigmaKa() * trackNeg.tofNSigmaKa()) < 9.0) {
            selectNegKaon = true;
          }
        }
      }

      if (topolD0) {
        statusD0NoPid = 1;
        if (pdgPositive == 211 && pdgNegative == -321) {
          statusD0PerfectPid = 1;
        }
        if (selectPosPion && selectNegKaon) {
          statusD0 = 1;
        }
      }
      if (topolD0bar) {
        statusD0barNoPid = 1;
        if (pdgPositive == 321 && pdgNegative == -211) {
          statusD0barPerfectPid = 1;
        }
        if (selectNegPion && selectPosKaon) {
          statusD0bar = 1;
        }
      }
      hfSelD0CandidateparametrizedPID(statusD0NoPid, statusD0PerfectPid, statusD0, statusD0barNoPid, statusD0barPerfectPid, statusD0bar);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfCandidateSelectorD0ParametrizedPidRichIndexBuilder>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfCandidateSelectorD0ParametrizedPid>(cfgc));
  return workflow;
}
