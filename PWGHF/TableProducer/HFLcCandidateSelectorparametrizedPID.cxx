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

/// \file HFLcCandidateSelector.cxx
/// \brief Λc± → p± K∓ π± selection task
///
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch>, University and INFN SALERNO
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
//#include <iostream>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"
#include "ALICE3/DataModel/RICH.h"
#include "Common/Core/MC.h"
#include "Common/Core/PID/PIDResponse.h"
#include "ReconstructionDataFormats/PID.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::analysis::hf_cuts_lc_topkpi;

namespace o2::aod
{

namespace indices
{
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_INDEX_COLUMN(RICH, rich);
} // namespace indices
DECLARE_SOA_INDEX_TABLE_USER(RICHTracksIndex, Tracks, "RICHTRK", indices::TrackId, indices::RICHId);
} // namespace o2::aod

struct richIndexBuilder { // Builder of the RICH-track index linkage
  Builds<o2::aod::RICHTracksIndex> indB;
  void init(o2::framework::InitContext&) {}
};

/// Struct for applying Lc selection cuts
struct HFLcCandidateSelectorparametrizedPID {
  Produces<aod::HFSelLcCandidateparametrizedPID> hfSelLcCandidateparametrizedPID;
  Configurable<double> d_etaperfectPID{"d_etaperfectPID", 1.75, "Eta cut for perfect PID"};
  Configurable<double> d_normaliseddecaylengthxyCand{"d_normaliseddecaylengthxyCand", 3., "Normalised decay length"};
  Configurable<double> d_pTCandMin{"d_pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> d_pTCandMax{"d_pTCandMax", 36., "Upper bound of candidate pT"};
  Configurable<bool> d_FilterPID{"d_FilterPID", true, "Bool to use or not the PID at filtering level"};
  // TPC
  Configurable<double> d_pidTPCMinpT{"d_pidTPCMinpT", 0.1, "Lower bound of track pT for TPC PID"};
  Configurable<double> d_pidTPCMaxpT{"d_pidTPCMaxpT", 1., "Upper bound of track pT for TPC PID"};
  Configurable<double> d_nSigmaTPC{"d_nSigmaTPC", 3., "Nsigma cut on TPC only"};
  Configurable<double> d_nSigmaTPCCombined{"d_nSigmaTPCCombined", 5., "Nsigma cut on TPC combined with TOF"};
  //Configurable<double> d_TPCNClsFindablePIDCut{"d_TPCNClsFindablePIDCut", 70., "Lower bound of TPC findable clusters for good PID"};
  // TOF
  Configurable<double> d_pidTOFMinpT{"d_pidTOFMinpT", 0.5, "Lower bound of track pT for TOF PID"};
  Configurable<double> d_pidTOFMaxpT{"d_pidTOFMaxpT", 2.5, "Upper bound of track pT for TOF PID"};
  Configurable<double> d_nSigmaTOF{"d_nSigmaTOF", 3., "Nsigma cut on TOF only"};
  Configurable<double> d_nSigmaTOFCombined{"d_nSigmaTOFCombined", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_lc_topkpi::pTBins_v}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"Lc_to_p_K_pi_cuts", {hf_cuts_lc_topkpi::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "Lc candidate selection per pT bin"};

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

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T>
  bool selectionTopol(const T& candidate)
  {
    auto candpT = candidate.pt();

    int pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < d_pTCandMin || candpT >= d_pTCandMax) {
      return false;
    }

    // cosine of pointing angle
    if (candidate.cpa() <= cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }

    //candidate chi2PCA
    if (candidate.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
      return false;
    }

    if (candidate.decayLength() <= cuts->get(pTBin, "decay length")) {
      return false;
    }

    if (candidate.decayLengthXYNormalised() < d_normaliseddecaylengthxyCand) {
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
    int pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }

    // cut on daughter pT
    if (trackProton.pt() < cuts->get(pTBin, "pT p") || trackKaon.pt() < cuts->get(pTBin, "pT K") || trackPion.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    if (trackProton.globalIndex() == candidate.index0Id()) {
      if (std::abs(InvMassLcpKpi(candidate) - RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus)) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(InvMassLcpiKp(candidate) - RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus)) > cuts->get(pTBin, "m")) {
        return false;
      }
    }

    return true;
  }

  using Trks = soa::Join<aod::BigTracksPID, aod::Tracks, aod::RICHTracksIndex, aod::McTrackLabels, aod::TracksExtra>;
  void process(aod::HfCandProng3 const& candidates, Trks const& barreltracks, const aod::McParticles& mcParticles, const aod::RICHs&, const aod::FRICHs&)
  {
    for (auto& candidate : candidates) {

      // selection flag
      int statusLcPKPiNoPID = 0;
      int statusLcPKPiPerfectPID = 0;
      int statusLcPKPi = 0;
      int statusLcPiKPNoPID = 0;
      int statusLcPiKPPerfectPID = 0;
      int statusLcPiKP = 0;

      if (!(candidate.hfflag() & 1 << DecayType::LcToPKPi)) {
        hfSelLcCandidateparametrizedPID(statusLcPKPiNoPID, statusLcPKPiPerfectPID, statusLcPKPi, statusLcPiKPNoPID, statusLcPiKPPerfectPID, statusLcPiKP);
        continue;
      }

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelLcCandidateparametrizedPID(statusLcPKPiNoPID, statusLcPKPiPerfectPID, statusLcPKPi, statusLcPiKPNoPID, statusLcPiKPPerfectPID, statusLcPiKP);
        continue;
      }

      auto trackPos1 = candidate.index0_as<Trks>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.index1_as<Trks>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.index2_as<Trks>(); // positive daughter (negative for the antiparticles)

      bool topolLcpKpi = selectionTopolConjugate(candidate, trackPos1, trackNeg, trackPos2);
      bool topolLcpiKp = selectionTopolConjugate(candidate, trackPos2, trackNeg, trackPos1);
      if (!topolLcpKpi && !topolLcpiKp) {
        hfSelLcCandidateparametrizedPID(statusLcPKPiNoPID, statusLcPKPiPerfectPID, statusLcPKPi, statusLcPiKPNoPID, statusLcPiKPPerfectPID, statusLcPiKP);
        continue;
      }

      auto ptPos1Track = trackPos1.pt();
      auto ptPos2Track = trackPos2.pt();
      auto ptNegTrack = trackNeg.pt();

      auto etaPos1Track = std::abs(trackPos1.eta());
      auto etaPos2Track = std::abs(trackPos2.eta());
      auto etaNegTrack = std::abs(trackNeg.eta());

      const auto mcParticlePositive1 = trackPos1.mcParticle();
      const auto mcParticlePositive2 = trackPos2.mcParticle();
      const auto mcParticleNegative = trackNeg.mcParticle();
      int pdgPositive1 = mcParticlePositive1.pdgCode();
      int pdgPositive2 = mcParticlePositive2.pdgCode();
      int pdgNegative = mcParticleNegative.pdgCode();

      bool selectPos1Proton = false;
      bool selectPos1Pion = false;
      bool selectPos2Pion = false;
      bool selectPos2Proton = false;
      bool selectNegKaon = false;

      if (etaPos1Track >= d_etaperfectPID) {
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
          } else if (std::abs(trackPos1.tofNSigmaPr()) < 3.0) {
            selectPos1Proton = true;
          }
        }

        if (trackPos1.has_rich() && !trackPos1.hasTOF()) {
          if (std::abs(trackPos1.rich().richNsigmaPi()) < 3.0) {
            selectPos1Pion = true;
          } else if (std::abs(trackPos1.rich().richNsigmaPr()) < 3.0) {
            selectPos1Proton = true;
          }
        }

        if (trackPos1.has_rich() && trackPos1.hasTOF()) {
          if ((trackPos1.rich().richNsigmaPi() * trackPos1.rich().richNsigmaPi() + trackPos1.tofNSigmaPi() * trackPos1.tofNSigmaPi()) < 9.0) {
            selectPos1Pion = true;
          } else if ((trackPos1.rich().richNsigmaPr() * trackPos1.rich().richNsigmaPr() + trackPos1.tofNSigmaPr() * trackPos1.tofNSigmaPr()) < 9.0) {
            selectPos1Proton = true;
          }
        }
      }

      if (etaPos2Track >= d_etaperfectPID) {
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
          } else if (std::abs(trackPos2.tofNSigmaPr()) < 3.0) {
            selectPos2Proton = true;
          }
        }

        if (trackPos2.has_rich() && !trackPos2.hasTOF()) {
          if (std::abs(trackPos2.rich().richNsigmaPi()) < 3.0) {
            selectPos2Pion = true;
          } else if (std::abs(trackPos2.rich().richNsigmaPr()) < 3.0) {
            selectPos2Proton = true;
          }
        }

        if (trackPos2.has_rich() && trackPos2.hasTOF()) {
          if ((trackPos2.rich().richNsigmaPi() * trackPos2.rich().richNsigmaPi() + trackPos2.tofNSigmaPi() * trackPos2.tofNSigmaPi()) < 9.0) {
            selectPos2Pion = true;
          } else if ((trackPos2.rich().richNsigmaPr() * trackPos2.rich().richNsigmaPr() + trackPos2.tofNSigmaPr() * trackPos2.tofNSigmaPr()) < 9.0) {
            selectPos2Proton = true;
          }
        }
      }

      if (etaNegTrack >= d_etaperfectPID) {
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

        else if (trackNeg.has_rich() && !trackNeg.hasTOF() && std::abs(trackNeg.rich().richNsigmaKa()) < 3.0) {
          selectNegKaon = true;
        }

        else if (trackNeg.has_rich() && trackNeg.hasTOF() && (trackNeg.rich().richNsigmaKa() * trackNeg.rich().richNsigmaKa() + trackNeg.tofNSigmaKa() * trackNeg.tofNSigmaKa()) < 9.0) {
          selectNegKaon = true;
        }
      }

      if (topolLcpKpi) {
        statusLcPKPiNoPID = 1;
        if (pdgPositive1 == kProton && pdgPositive2 == kPiPlus && pdgNegative == kKMinus) {
          statusLcPKPiPerfectPID = 1;
        }
        if (selectPos1Proton && selectPos2Pion && selectNegKaon) {
          statusLcPKPi = 1;
        }
      }

      if (topolLcpiKp) {
        statusLcPiKPNoPID = 1;
        if (pdgPositive2 == kProton && pdgPositive1 == kPiPlus && pdgNegative == kKMinus) {
          statusLcPiKPPerfectPID = 1;
        }
        if (selectPos2Proton && selectPos1Pion && selectNegKaon) {
          statusLcPiKP = 1;
        }
      }
      hfSelLcCandidateparametrizedPID(statusLcPKPiNoPID, statusLcPKPiPerfectPID, statusLcPKPi, statusLcPiKPNoPID, statusLcPiKPPerfectPID, statusLcPiKP);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<richIndexBuilder>(cfgc));
  workflow.push_back(adaptAnalysisTask<HFLcCandidateSelectorparametrizedPID>(cfgc, TaskName{"hf-lc-candidate-selector-parametrizedPID"}));
  return workflow;
}
