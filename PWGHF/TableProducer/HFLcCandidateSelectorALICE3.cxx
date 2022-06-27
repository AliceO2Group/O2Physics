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
#include "Common/DataModel/PIDResponse.h"
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
struct HFLcCandidateSelectorALICE3 {
  Produces<aod::HFSelLcCandidateALICE3> hfSelLcCandidateALICE3;

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
      int statusLcPKPiTOFPID = 0;
      int statusLcPKPiTOFplusRICHPID = 0;
      int statusLcPiKPNoPID = 0;
      int statusLcPiKPPerfectPID = 0;
      int statusLcPiKPTOFPID = 0;
      int statusLcPiKPTOFplusRICHPID = 0;

      if (!(candidate.hfflag() & 1 << DecayType::LcToPKPi)) {
        hfSelLcCandidateALICE3(statusLcPKPiNoPID, statusLcPKPiPerfectPID, statusLcPKPiTOFPID, statusLcPKPiTOFplusRICHPID, statusLcPiKPNoPID, statusLcPiKPPerfectPID, statusLcPiKPTOFPID, statusLcPiKPTOFplusRICHPID);
        continue;
      }

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelLcCandidateALICE3(statusLcPKPiNoPID, statusLcPKPiPerfectPID, statusLcPKPiTOFPID, statusLcPKPiTOFplusRICHPID, statusLcPiKPNoPID, statusLcPiKPPerfectPID, statusLcPiKPTOFPID, statusLcPiKPTOFplusRICHPID);
        continue;
      }

      auto trackPos1 = candidate.index0_as<Trks>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.index1_as<Trks>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.index2_as<Trks>(); // positive daughter (negative for the antiparticles)

      auto momentumPos1Track = trackPos1.p();
      auto momentumNegTrack = trackNeg.p();
      auto momentumPos2Track = trackPos2.p();

      bool topolLcpKpi = selectionTopolConjugate(candidate, trackPos1, trackNeg, trackPos2);
      bool topolLcpiKp = selectionTopolConjugate(candidate, trackPos2, trackNeg, trackPos1);

      if (!topolLcpKpi && !topolLcpiKp) {
        hfSelLcCandidateALICE3(statusLcPKPiNoPID, statusLcPKPiPerfectPID, statusLcPKPiTOFPID, statusLcPKPiTOFplusRICHPID, statusLcPiKPNoPID, statusLcPiKPPerfectPID, statusLcPiKPTOFPID, statusLcPiKPTOFplusRICHPID);
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

      if (topolLcpKpi) {
        statusLcPKPiNoPID = 1;
        if (pdgPositive1 == kProton && pdgPositive2 == kPiPlus && pdgNegative == kKMinus) {
          statusLcPKPiPerfectPID = 1;
        }
        if (std::abs(nSigmaTOFPos1Proton) < 3.0 && std::abs(nSigmaTOFPos2Pion) < 3.0 && std::abs(nSigmaTOFNegKaon) < 3.0) {
          statusLcPKPiTOFPID = 1;
        }
        if (selectProtonPos1TOFplusRICH && selectPionPos2TOFplusRICH && selectKaonTOFplusRICH) {
          statusLcPKPiTOFplusRICHPID = 1;
        }
      }

      if (topolLcpiKp) {
        statusLcPiKPNoPID = 1;
        if (pdgPositive2 == kProton && pdgPositive1 == kPiPlus && pdgNegative == kKMinus) {
          statusLcPiKPPerfectPID = 1;
        }
        if (std::abs(nSigmaTOFPos2Proton) < 3.0 && std::abs(nSigmaTOFPos1Pion) < 3.0 && std::abs(nSigmaTOFNegKaon) < 3.0) {
          statusLcPiKPTOFPID = 1;
        }
        if (selectProtonPos2TOFplusRICH && selectPionPos1TOFplusRICH && selectKaonTOFplusRICH) {
          statusLcPiKPTOFplusRICHPID = 1;
        }
      }
      //std::cout << "status = " << statusLcPKPiNoPID << "\t" << statusLcPiKPNoPID << "\n";
      hfSelLcCandidateALICE3(statusLcPKPiNoPID, statusLcPKPiPerfectPID, statusLcPKPiTOFPID, statusLcPKPiTOFplusRICHPID, statusLcPiKPNoPID, statusLcPiKPPerfectPID, statusLcPiKPTOFPID, statusLcPiKPTOFplusRICHPID);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<richIndexBuilder>(cfgc));
  workflow.push_back(adaptAnalysisTask<HFLcCandidateSelectorALICE3>(cfgc, TaskName{"hf-lc-candidate-selector-ALICE3"}));
  return workflow;
}
