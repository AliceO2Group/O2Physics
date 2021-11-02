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

/// \file HFD0CandidateSelector.cxx
/// \brief D0(bar) → π± K∓ selection task
///
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
using namespace o2::aod::hf_cand_prong2;
using namespace o2::analysis::hf_cuts_d0_topik;

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

/// Struct for applying D0 selection cuts
struct HFD0CandidateSelectorparametrizedPID {
  Produces<aod::HFSelD0CandidateparametrizedPID> hfSelD0CandidateparametrizedPID;
  Configurable<double> d_etaperfectPID{"d_etaperfectPID", 1.75, "Eta cut for perfect PID"};
  Configurable<double> d_pTCandMin{"d_pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> d_pTCandMax{"d_pTCandMax", 50., "Upper bound of candidate pT"};
  // TPC
  Configurable<double> d_pidTPCMinpT{"d_pidTPCMinpT", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> d_pidTPCMaxpT{"d_pidTPCMaxpT", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> d_nSigmaTPC{"d_nSigmaTPC", 3., "Nsigma cut on TPC only"};
  Configurable<double> d_nSigmaTPCCombined{"d_nSigmaTPCCombined", 5., "Nsigma cut on TPC combined with TOF"};
  //Configurable<double> d_TPCNClsFindablePIDCut{"d_TPCNClsFindablePIDCut", 50., "Lower bound of TPC findable clusters for good PID"};
  // TOF
  Configurable<double> d_pidTOFMinpT{"d_pidTOFMinpT", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> d_pidTOFMaxpT{"d_pidTOFMaxpT", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> d_nSigmaTOF{"d_nSigmaTOF", 3., "Nsigma cut on TOF only"};
  Configurable<double> d_nSigmaTOFCombined{"d_nSigmaTOFCombined", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_d0_topik::pTBins_v}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"D0_to_pi_K_cuts", {hf_cuts_d0_topik::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "D0 candidate selection per pT bin"};

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
    auto pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < d_pTCandMin || candpT >= d_pTCandMax) {
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
    //if (candidate.chi2PCA() > cuts[pTBin][1]) return false;

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
      //return false; // add back when getter fixed
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
    auto pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }

    // invariant-mass cut
    if (trackPion.sign() > 0) {
      if (std::abs(InvMassD0(candidate) - RecoDecay::getMassPDG(pdg::Code::kD0)) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(InvMassD0bar(candidate) - RecoDecay::getMassPDG(pdg::Code::kD0)) > cuts->get(pTBin, "m")) {
        return false;
      }
    }

    // cut on daughter pT
    if (trackPion.pt() < cuts->get(pTBin, "pT Pi") || trackKaon.pt() < cuts->get(pTBin, "pT K")) {
      return false;
    }

    // cut on daughter DCA - need to add secondary vertex constraint here
    if (std::abs(trackPion.dcaPrim0()) > cuts->get(pTBin, "d0pi") || std::abs(trackKaon.dcaPrim0()) > cuts->get(pTBin, "d0K")) {
      return false;
    }

    // cut on cos(theta*)
    if (trackPion.sign() > 0) {
      if (std::abs(CosThetaStarD0(candidate)) > cuts->get(pTBin, "cos theta*")) {
        return false;
      }
    } else {
      if (std::abs(CosThetaStarD0bar(candidate)) > cuts->get(pTBin, "cos theta*")) {
        return false;
      }
    }

    return true;
  }

  using Trks = soa::Join<aod::BigTracksPID, aod::Tracks, aod::RICHTracksIndex, aod::McTrackLabels, aod::TracksExtra>;
  void process(aod::HfCandProng2 const& candidates, Trks const& barreltracks, const aod::McParticles& mcParticles, const aod::RICHs&, const aod::FRICHs&)
  {

    for (auto& candidate : candidates) {

      // selection flag
      int statusD0NoPID = 0;
      int statusD0PerfectPID = 0;
      int statusD0 = 0;
      int statusD0barNoPID = 0;
      int statusD0barPerfectPID = 0;
      int statusD0bar = 0;

      if (!(candidate.hfflag() & 1 << DecayType::D0ToPiK)) {
        hfSelD0CandidateparametrizedPID(statusD0NoPID, statusD0PerfectPID, statusD0, statusD0barNoPID, statusD0barPerfectPID, statusD0bar);
        continue;
      }

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelD0CandidateparametrizedPID(statusD0NoPID, statusD0PerfectPID, statusD0, statusD0barNoPID, statusD0barPerfectPID, statusD0bar);
        continue;
      }

      auto trackPos = candidate.index0_as<Trks>();
      auto trackNeg = candidate.index1_as<Trks>();

      auto momentumPosTrack = trackPos.p();
      auto momentumNegTrack = trackNeg.p();
      auto ptPosTrack = trackPos.pt();
      auto ptNegTrack = trackNeg.pt();
      auto etaPosTrack = std::abs(trackPos.eta());
      auto etaNegTrack = std::abs(trackNeg.eta());

      bool topolD0 = selectionTopolConjugate(candidate, trackPos, trackNeg);
      bool topolD0bar = selectionTopolConjugate(candidate, trackNeg, trackPos);

      if (!topolD0 && !topolD0bar) {
        hfSelD0CandidateparametrizedPID(statusD0NoPID, statusD0PerfectPID, statusD0, statusD0barNoPID, statusD0barPerfectPID, statusD0bar);
        continue;
      }

      const auto mcParticlePositive = trackPos.mcParticle();
      const auto mcParticleNegative = trackNeg.mcParticle();
      int pdgPositive = mcParticlePositive.pdgCode();
      int pdgNegative = mcParticleNegative.pdgCode();

      bool selectPosPion = false;
      bool selectNegKaon = false;
      bool selectNegPion = false;
      bool selectPosKaon = false;

      if (etaPosTrack >= d_etaperfectPID) {
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
          } else if (std::abs(trackPos.tofNSigmaKa()) < 3.0) {
            selectPosKaon = true;
          }
        }

        if (trackPos.has_rich() && !trackPos.hasTOF()) {
          if (std::abs(trackPos.rich().richNsigmaPi()) < 3.0) {
            selectPosPion = true;
          } else if (std::abs(trackPos.rich().richNsigmaKa()) < 3.0) {
            selectPosKaon = true;
          }
        }

        if (trackPos.has_rich() && trackPos.hasTOF()) {
          if ((trackPos.rich().richNsigmaPi() * trackPos.rich().richNsigmaPi() + trackPos.tofNSigmaPi() * trackPos.tofNSigmaPi()) < 9.0) {
            selectPosPion = true;
          } else if ((trackPos.rich().richNsigmaKa() * trackPos.rich().richNsigmaKa() + trackPos.tofNSigmaKa() * trackPos.tofNSigmaKa()) < 9.0) {
            selectPosKaon = true;
          }
        }
      }

      if (etaNegTrack >= d_etaperfectPID) {
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
          } else if (std::abs(trackNeg.tofNSigmaKa()) < 3.0) {
            selectNegKaon = true;
          }
        }

        if (trackNeg.has_rich() && !trackNeg.hasTOF()) {
          if (std::abs(trackNeg.rich().richNsigmaPi()) < 3.0) {
            selectNegPion = true;
          } else if (std::abs(trackNeg.rich().richNsigmaKa()) < 3.0) {
            selectNegKaon = true;
          }
        }

        if (trackNeg.has_rich() && trackNeg.hasTOF()) {
          if ((trackNeg.rich().richNsigmaPi() * trackNeg.rich().richNsigmaPi() + trackNeg.tofNSigmaPi() * trackNeg.tofNSigmaPi()) < 9.0) {
            selectNegPion = true;
          } else if ((trackNeg.rich().richNsigmaKa() * trackNeg.rich().richNsigmaKa() + trackNeg.tofNSigmaKa() * trackNeg.tofNSigmaKa()) < 9.0) {
            selectNegKaon = true;
          }
        }
      }

      if (topolD0) {
        statusD0NoPID = 1;
        if (pdgPositive == 211 && pdgNegative == -321) {
          statusD0PerfectPID = 1;
        }
        if (selectPosPion && selectNegKaon) {
          statusD0 = 1;
        }
      } else if (topolD0bar) {
        statusD0barNoPID = 1;
        if (pdgPositive == 321 && pdgNegative == -211) {
          statusD0barPerfectPID = 1;
        }
        if (selectNegPion && selectPosKaon) {
          statusD0bar = 1;
        }
      }
      hfSelD0CandidateparametrizedPID(statusD0NoPID, statusD0PerfectPID, statusD0, statusD0barNoPID, statusD0barPerfectPID, statusD0bar);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<richIndexBuilder>(cfgc));
  workflow.push_back(adaptAnalysisTask<HFD0CandidateSelectorparametrizedPID>(cfgc, TaskName{"hf-d0-candidate-selector-parametrizedPID"}));
  return workflow;
}
