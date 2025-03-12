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

/// \file candidateSelectorLcAlice3.cxx
/// \brief Λc± → p± K∓ π± selection task
///
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch>, University and INFN SALERNO
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include <vector>

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

struct HfCandidateSelectorLcAlice3RichIndexBuilder { // Builder of the RICH-track index linkage
  Builds<o2::aod::RICHTracksIndex> indB;

  void init(InitContext&) {}
};

/// Struct for applying Lc selection cuts
struct HfCandidateSelectorLcAlice3 {
  Produces<aod::HfSelLcAlice3> hfSelLcCandidateALICE3;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.1, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 1., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.5, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 2.5, "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<double> decayLengthXYNormalisedMin{"decayLengthXYNormalisedMin", 3., "Min. normalised decay length"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_p_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lc_to_p_k_pi::Cuts[0], hf_cuts_lc_to_p_k_pi::NBinsPt, hf_cuts_lc_to_p_k_pi::NCutVars, hf_cuts_lc_to_p_k_pi::labelsPt, hf_cuts_lc_to_p_k_pi::labelsCutVar}, "Lc candidate selection per pT bin"};

  HfHelper hfHelper;

  using TracksSel = soa::Join<aod::TracksWExtra, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::RICHTracksIndex, aod::McTrackLabels>;

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T>
  bool selectionTopol(const T& candidate)
  {
    auto candpT = candidate.pt();

    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptCandMin || candpT >= ptCandMax) {
      return false;
    }

    // cosine of pointing angle
    if (candidate.cpa() <= cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }

    // candidate chi2PCA
    if (candidate.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
      return false;
    }

    if (candidate.decayLength() <= cuts->get(pTBin, "decay length")) {
      return false;
    }

    if (candidate.decayLengthXYNormalised() < decayLengthXYNormalisedMin) {
      return false;
    }

    return true;
  }

  /// Conjugate-dependent topological cuts
  /// \param candidate is candidate
  /// \param trackProton is the track with the proton hypothesis
  /// \param trackPion is the track with the pion hypothesis
  /// \param trackKaon is the track with the kaon hypothesis
  /// \return true if candidate passes all cuts for the given Conjugate
  template <typename T1, typename T2>
  bool selectionTopolConjugate(const T1& candidate, const T2& trackProton, const T2& trackKaon, const T2& trackPion)
  {

    auto candpT = candidate.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // cut on daughter pT
    if (trackProton.pt() < cuts->get(pTBin, "pT p") || trackKaon.pt() < cuts->get(pTBin, "pT K") || trackPion.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    if (trackProton.globalIndex() == candidate.prong0Id()) {
      if (std::abs(hfHelper.invMassLcToPKPi(candidate) - o2::constants::physics::MassLambdaCPlus) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(hfHelper.invMassLcToPiKP(candidate) - o2::constants::physics::MassLambdaCPlus) > cuts->get(pTBin, "m")) {
        return false;
      }
    }

    return true;
  }

  void process(aod::HfCand3Prong const& candidates,
               TracksSel const&,
               aod::McParticles const&,
               aod::RICHs const&,
               aod::FRICHs const&)
  {
    for (const auto& candidate : candidates) {

      // selection flag

      int statusLcToPKPiNoPid = 0;
      int statusLcToPKPiPerfectPid = 0;
      int statusLcToPKPiTofPid = 0;
      int statusLcToPKPiTofPlusRichPid = 0;
      int statusLcToPiKPNoPid = 0;
      int statusLcToPiKPPerfectPid = 0;
      int statusLcToPiKPTofPid = 0;
      int statusLcToPiKPTofPlusRichPid = 0;

      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        hfSelLcCandidateALICE3(statusLcToPKPiNoPid, statusLcToPKPiPerfectPid, statusLcToPKPiTofPid, statusLcToPKPiTofPlusRichPid, statusLcToPiKPNoPid, statusLcToPiKPPerfectPid, statusLcToPiKPTofPid, statusLcToPiKPTofPlusRichPid);
        continue;
      }

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelLcCandidateALICE3(statusLcToPKPiNoPid, statusLcToPKPiPerfectPid, statusLcToPKPiTofPid, statusLcToPKPiTofPlusRichPid, statusLcToPiKPNoPid, statusLcToPiKPPerfectPid, statusLcToPiKPTofPid, statusLcToPiKPTofPlusRichPid);
        continue;
      }

      auto trackPos1 = candidate.prong0_as<TracksSel>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<TracksSel>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<TracksSel>(); // positive daughter (negative for the antiparticles)

      auto momentumPos1Track = trackPos1.p();
      auto momentumNegTrack = trackNeg.p();
      auto momentumPos2Track = trackPos2.p();

      bool topolLcToPKPi = selectionTopolConjugate(candidate, trackPos1, trackNeg, trackPos2);
      bool topolLcToPiKP = selectionTopolConjugate(candidate, trackPos2, trackNeg, trackPos1);

      if (!topolLcToPKPi && !topolLcToPiKP) {
        hfSelLcCandidateALICE3(statusLcToPKPiNoPid, statusLcToPKPiPerfectPid, statusLcToPKPiTofPid, statusLcToPKPiTofPlusRichPid, statusLcToPiKPNoPid, statusLcToPiKPPerfectPid, statusLcToPiKPTofPid, statusLcToPiKPTofPlusRichPid);
        continue;
      }

      const auto mcParticlePositive1 = trackPos1.mcParticle();
      const auto mcParticleNegative = trackNeg.mcParticle();
      const auto mcParticlePositive2 = trackPos2.mcParticle();

      int pdgPositive1 = mcParticlePositive1.pdgCode();
      int pdgNegative = mcParticleNegative.pdgCode();
      int pdgPositive2 = mcParticlePositive2.pdgCode();

      float nSigmaTOFPos1Proton = -5000.0;
      float nSigmaRICHPos1Proton = -5000.0;
      float nSigmaTOFPos2Proton = -5000.0;
      float nSigmaRICHPos2Proton = -5000.0;
      float nSigmaTOFNegKaon = -5000.0;
      float nSigmaRICHNegKaon = -5000.0;
      float nSigmaTOFPos1Pion = -5000.0;
      float nSigmaRICHPos1Pion = -5000.0;
      float nSigmaTOFPos2Pion = -5000.0;
      float nSigmaRICHPos2Pion = -5000.0;

      if (trackPos1.hasTOF()) {
        nSigmaTOFPos1Proton = trackPos1.tofNSigmaPr();
        nSigmaTOFPos1Pion = trackPos1.tofNSigmaPi();
      }
      if (trackPos2.hasTOF()) {
        nSigmaTOFPos2Pion = trackPos2.tofNSigmaPi();
        nSigmaTOFPos2Proton = trackPos2.tofNSigmaPr();
      }
      if (trackNeg.hasTOF()) {
        nSigmaTOFNegKaon = trackNeg.tofNSigmaKa();
      }

      if (trackPos1.has_rich()) {
        nSigmaRICHPos1Proton = trackPos1.rich().richNsigmaPr();
        nSigmaRICHPos1Pion = trackPos1.rich().richNsigmaPi();
      }
      if (trackPos2.has_rich()) {
        nSigmaRICHPos2Pion = trackPos2.rich().richNsigmaPi();
        nSigmaRICHPos2Proton = trackPos2.rich().richNsigmaPr();
      }
      if (trackNeg.has_rich()) {
        nSigmaRICHNegKaon = trackNeg.rich().richNsigmaKa();
      }

      bool selectProtonPos1TOFplusRICH = false;
      bool selectProtonPos2TOFplusRICH = false;
      bool selectPionPos1TOFplusRICH = false;
      bool selectPionPos2TOFplusRICH = false;
      bool selectKaonTOFplusRICH = false;

      if ((momentumPos1Track < 4.0 && std::abs(nSigmaTOFPos1Proton) < 3.0)) {
        selectProtonPos1TOFplusRICH = true;
      } else if ((momentumPos1Track > 4.0 && trackPos1.has_rich() && (nSigmaRICHPos1Proton * nSigmaRICHPos1Proton + nSigmaTOFPos1Proton * nSigmaTOFPos1Proton) < 9.0)) {
        selectProtonPos1TOFplusRICH = true;
      }
      if ((momentumPos2Track < 4.0 && std::abs(nSigmaTOFPos2Proton) < 3.0)) {
        selectProtonPos2TOFplusRICH = true;
      } else if ((momentumPos2Track > 4.0 && trackPos2.has_rich() && (nSigmaRICHPos2Proton * nSigmaRICHPos2Proton + nSigmaTOFPos2Proton * nSigmaTOFPos2Proton) < 9.0)) {
        selectProtonPos2TOFplusRICH = true;
      }

      if ((momentumPos1Track < 0.6 && std::abs(nSigmaTOFPos1Pion) < 3.0)) {
        selectPionPos1TOFplusRICH = true;
      } else if ((momentumPos1Track > 0.6 && trackPos1.has_rich() && (nSigmaRICHPos1Pion * nSigmaRICHPos1Pion + nSigmaTOFPos1Pion * nSigmaTOFPos1Pion) < 9.0)) {
        selectPionPos1TOFplusRICH = true;
      }
      if ((momentumPos2Track < 0.6 && std::abs(nSigmaTOFPos2Pion) < 3.0)) {
        selectPionPos2TOFplusRICH = true;
      } else if ((momentumPos2Track > 0.6 && trackPos2.has_rich() && (nSigmaRICHPos2Pion * nSigmaRICHPos2Pion + nSigmaTOFPos2Pion * nSigmaTOFPos2Pion) < 9.0)) {
        selectPionPos2TOFplusRICH = true;
      }

      if ((momentumNegTrack < 2.0 && std::abs(nSigmaTOFNegKaon) < 3.0)) {
        selectKaonTOFplusRICH = true;
      } else if ((momentumNegTrack > 2.0 && trackNeg.has_rich() && (nSigmaRICHNegKaon * nSigmaRICHNegKaon + nSigmaTOFNegKaon * nSigmaTOFNegKaon) < 9.0)) {
        selectKaonTOFplusRICH = true;
      }

      if (topolLcToPKPi) {
        statusLcToPKPiNoPid = 1;
        if (pdgPositive1 == kProton && pdgPositive2 == kPiPlus && pdgNegative == kKMinus) {
          statusLcToPKPiPerfectPid = 1;
        }
        if (std::abs(nSigmaTOFPos1Proton) < 3.0 && std::abs(nSigmaTOFPos2Pion) < 3.0 && std::abs(nSigmaTOFNegKaon) < 3.0) {
          statusLcToPKPiTofPid = 1;
        }
        if (selectProtonPos1TOFplusRICH && selectPionPos2TOFplusRICH && selectKaonTOFplusRICH) {
          statusLcToPKPiTofPlusRichPid = 1;
        }
      }

      if (topolLcToPiKP) {
        statusLcToPiKPNoPid = 1;
        if (pdgPositive2 == kProton && pdgPositive1 == kPiPlus && pdgNegative == kKMinus) {
          statusLcToPiKPPerfectPid = 1;
        }
        if (std::abs(nSigmaTOFPos2Proton) < 3.0 && std::abs(nSigmaTOFPos1Pion) < 3.0 && std::abs(nSigmaTOFNegKaon) < 3.0) {
          statusLcToPiKPTofPid = 1;
        }
        if (selectProtonPos2TOFplusRICH && selectPionPos1TOFplusRICH && selectKaonTOFplusRICH) {
          statusLcToPiKPTofPlusRichPid = 1;
        }
      }
      hfSelLcCandidateALICE3(statusLcToPKPiNoPid, statusLcToPKPiPerfectPid, statusLcToPKPiTofPid, statusLcToPKPiTofPlusRichPid, statusLcToPiKPNoPid, statusLcToPiKPPerfectPid, statusLcToPiKPTofPid, statusLcToPiKPTofPlusRichPid);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfCandidateSelectorLcAlice3RichIndexBuilder>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfCandidateSelectorLcAlice3>(cfgc));
  return workflow;
}
