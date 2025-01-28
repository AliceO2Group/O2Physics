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

/// \file candidateSelectorLcParametrizedPid.cxx
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

struct HfCandidateSelectorLcParametrizedPidRichIndexBuilder { // Builder of the RICH-track index linkage
  Builds<o2::aod::RICHTracksIndex> indB;

  void init(InitContext&) {}
};

/// Struct for applying Lc selection cuts
struct HfCandidateSelectorLcParametrizedPid {
  Produces<aod::HfSelLcParametrizedPid> hfSelLcCandidateparametrizedPID;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  Configurable<double> etaPerfectPidMax{"etaPerfectPidMax", 1.75, "Eta cut for perfect PID"};
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
  Configurable<double> decayLengthXYNormalisedMin{"decayLengthXYNormalisedMin", 3., "Normalised decay length"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_p_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lc_to_p_k_pi::cuts[0], hf_cuts_lc_to_p_k_pi::nBinsPt, hf_cuts_lc_to_p_k_pi::nCutVars, hf_cuts_lc_to_p_k_pi::labelsPt, hf_cuts_lc_to_p_k_pi::labelsCutVar}, "Lc candidate selection per pT bin"};

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
      int statusLcToPKPi = 0;
      int statusLcToPiKPNoPid = 0;
      int statusLcToPiKPPerfectPid = 0;
      int statusLcToPiKP = 0;

      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        hfSelLcCandidateparametrizedPID(statusLcToPKPiNoPid, statusLcToPKPiPerfectPid, statusLcToPKPi, statusLcToPiKPNoPid, statusLcToPiKPPerfectPid, statusLcToPiKP);
        continue;
      }

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelLcCandidateparametrizedPID(statusLcToPKPiNoPid, statusLcToPKPiPerfectPid, statusLcToPKPi, statusLcToPiKPNoPid, statusLcToPiKPPerfectPid, statusLcToPiKP);
        continue;
      }

      auto trackPos1 = candidate.prong0_as<TracksSel>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<TracksSel>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<TracksSel>(); // positive daughter (negative for the antiparticles)

      bool topolLcToPKPi = selectionTopolConjugate(candidate, trackPos1, trackNeg, trackPos2);
      bool topolLcToPiKP = selectionTopolConjugate(candidate, trackPos2, trackNeg, trackPos1);
      if (!topolLcToPKPi && !topolLcToPiKP) {
        hfSelLcCandidateparametrizedPID(statusLcToPKPiNoPid, statusLcToPKPiPerfectPid, statusLcToPKPi, statusLcToPiKPNoPid, statusLcToPiKPPerfectPid, statusLcToPiKP);
        continue;
      }

      auto ptPos1Track = trackPos1.pt();
      auto ptPos2Track = trackPos2.pt();
      auto ptNegTrack = trackNeg.pt();

      auto etaPos1Track = std::abs(trackPos1.eta());
      auto etaPos2Track = std::abs(trackPos2.eta());
      auto etaNegTrack = std::abs(trackNeg.eta());

      int pdgPositive1 = 0;
      int pdgPositive2 = 0;
      int pdgNegative = 0;
      if (trackPos1.has_mcParticle()) {
        pdgPositive1 = trackPos1.mcParticle().pdgCode();
      }
      if (trackPos2.has_mcParticle()) {
        pdgPositive2 = trackPos2.mcParticle().pdgCode();
      }
      if (trackNeg.has_mcParticle()) {
        pdgNegative = trackNeg.mcParticle().pdgCode();
      }

      bool selectPos1Proton = false;
      bool selectPos1Pion = false;
      bool selectPos2Pion = false;
      bool selectPos2Proton = false;
      bool selectNegKaon = false;

      if (etaPos1Track >= etaPerfectPidMax) {
        if (ptPos1Track < (19.58 / std::cosh(etaPos1Track))) {
          if (pdgPositive1 == 2212)
            selectPos1Proton = true;
          if (ptPos1Track < (11.65 / std::cosh(etaPos1Track)) && pdgPositive1 == 211)
            selectPos1Pion = true;
        } else {
          selectPos1Proton = true;
          selectPos1Pion = true;
        }
      } else {
        if (trackPos1.hasTOF()) {
          if (std::abs(trackPos1.tofNSigmaPi()) < 3.0) {
            selectPos1Pion = true;
          }
          if (std::abs(trackPos1.tofNSigmaPr()) < 3.0) {
            selectPos1Proton = true;
          }
        }

        if (trackPos1.has_rich() && !trackPos1.hasTOF()) {
          if (std::abs(trackPos1.rich().richNsigmaPi()) < 3.0) {
            selectPos1Pion = true;
          }
          if (std::abs(trackPos1.rich().richNsigmaPr()) < 3.0) {
            selectPos1Proton = true;
          }
        }

        if (trackPos1.has_rich() && trackPos1.hasTOF()) {
          if ((trackPos1.rich().richNsigmaPi() * trackPos1.rich().richNsigmaPi() + trackPos1.tofNSigmaPi() * trackPos1.tofNSigmaPi()) < 9.0) {
            selectPos1Pion = true;
          }
          if ((trackPos1.rich().richNsigmaPr() * trackPos1.rich().richNsigmaPr() + trackPos1.tofNSigmaPr() * trackPos1.tofNSigmaPr()) < 9.0) {
            selectPos1Proton = true;
          }
        }
      }

      if (etaPos2Track >= etaPerfectPidMax) {
        if (ptPos2Track < (19.58 / std::cosh(etaPos2Track))) {
          if (pdgPositive2 == 2212)
            selectPos2Proton = true;
          if (ptPos2Track < (11.65 / std::cosh(etaPos2Track)) && pdgPositive2 == 211)
            selectPos2Pion = true;
        } else {
          selectPos2Proton = true;
          selectPos2Pion = true;
        }
      } else {
        if (trackPos2.hasTOF()) {
          if (std::abs(trackPos2.tofNSigmaPi()) < 3.0) {
            selectPos2Pion = true;
          }
          if (std::abs(trackPos2.tofNSigmaPr()) < 3.0) {
            selectPos2Proton = true;
          }
        }

        if (trackPos2.has_rich() && !trackPos2.hasTOF()) {
          if (std::abs(trackPos2.rich().richNsigmaPi()) < 3.0) {
            selectPos2Pion = true;
          }
          if (std::abs(trackPos2.rich().richNsigmaPr()) < 3.0) {
            selectPos2Proton = true;
          }
        }

        if (trackPos2.has_rich() && trackPos2.hasTOF()) {
          if ((trackPos2.rich().richNsigmaPi() * trackPos2.rich().richNsigmaPi() + trackPos2.tofNSigmaPi() * trackPos2.tofNSigmaPi()) < 9.0) {
            selectPos2Pion = true;
          }
          if ((trackPos2.rich().richNsigmaPr() * trackPos2.rich().richNsigmaPr() + trackPos2.tofNSigmaPr() * trackPos2.tofNSigmaPr()) < 9.0) {
            selectPos2Proton = true;
          }
        }
      }

      if (etaNegTrack >= etaPerfectPidMax) {
        if (ptNegTrack < (19.58 / std::cosh(etaNegTrack))) {
          if (pdgNegative == -321)
            selectNegKaon = true;
        } else {
          selectNegKaon = true;
        }
      } else {
        if (trackNeg.hasTOF() && std::abs(trackNeg.tofNSigmaKa()) < 3.0) {
          selectNegKaon = true;
        }

        if (trackNeg.has_rich() && !trackNeg.hasTOF() && std::abs(trackNeg.rich().richNsigmaKa()) < 3.0) {
          selectNegKaon = true;
        }

        if (trackNeg.has_rich() && trackNeg.hasTOF() && (trackNeg.rich().richNsigmaKa() * trackNeg.rich().richNsigmaKa() + trackNeg.tofNSigmaKa() * trackNeg.tofNSigmaKa()) < 9.0) {
          selectNegKaon = true;
        }
      }

      if (topolLcToPKPi) {
        statusLcToPKPiNoPid = 1;
        if (pdgPositive1 == kProton && pdgPositive2 == kPiPlus && pdgNegative == kKMinus) {
          statusLcToPKPiPerfectPid = 1;
        }
        if (selectPos1Proton && selectPos2Pion && selectNegKaon) {
          statusLcToPKPi = 1;
        }
      }

      if (topolLcToPiKP) {
        statusLcToPiKPNoPid = 1;
        if (pdgPositive2 == kProton && pdgPositive1 == kPiPlus && pdgNegative == kKMinus) {
          statusLcToPiKPPerfectPid = 1;
        }
        if (selectPos2Proton && selectPos1Pion && selectNegKaon) {
          statusLcToPiKP = 1;
        }
      }
      hfSelLcCandidateparametrizedPID(statusLcToPKPiNoPid, statusLcToPKPiPerfectPid, statusLcToPKPi, statusLcToPiKPNoPid, statusLcToPiKPPerfectPid, statusLcToPiKP);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfCandidateSelectorLcParametrizedPidRichIndexBuilder>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfCandidateSelectorLcParametrizedPid>(cfgc));
  return workflow;
}
