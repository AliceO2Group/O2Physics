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

/// \file taskMini.cxx
/// \brief Mini version of the D0 analysis chain
///
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#include "Tutorials/PWGHF/DataModelMini.h"
//
#include "PWGHF/Core/HfHelper.h"
//
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TString.h>

#include <array>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Candidate creation =====================================================================

/// Candidate creator
/// Reconstruction of heavy-flavour 2-prong decay candidates
struct HfTaskMiniCandidateCreator2Prong {
  Produces<aod::HfTCand2ProngBase> rowCandidateBase;

  // vertexing parameters
  Configurable<float> magneticField{"magneticField", 5., "magnetic field [kG]"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<float> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<float> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<float> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<float> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};

  o2::vertexing::DCAFitterN<2> fitter; // 2-prong vertex fitter
  double massPiK{0.};
  double massKPi{0.};

  using TracksWithCov = soa::Join<Tracks, TracksCov>;

  OutputObj<TH1F> hMass{TH1F("hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", 500, 0., 5.)};

  void init(InitContext&)
  {
    // Configure the vertexer
    fitter.setBz(magneticField);
    fitter.setPropagateToPCA(propagateToPCA);
    fitter.setMaxR(maxR);
    fitter.setMaxDZIni(maxDZIni);
    fitter.setMinParamChange(minParamChange);
    fitter.setMinRelChi2Change(minRelChi2Change);
    fitter.setUseAbsDCA(useAbsDCA);
  }

  void process(aod::Collisions const&,
               aod::HfT2Prongs const& rowsTrackIndexProng2,
               TracksWithCov const&)
  {
    // loop over pairs of track indices
    for (const auto& rowTrackIndexProng2 : rowsTrackIndexProng2) {
      auto track0 = rowTrackIndexProng2.prong0_as<TracksWithCov>();
      auto track1 = rowTrackIndexProng2.prong1_as<TracksWithCov>();
      auto trackParVarPos1 = getTrackParCov(track0);
      auto trackParVarNeg1 = getTrackParCov(track1);
      auto collision = track0.collision();

      // reconstruct the 2-prong secondary vertex
      if (fitter.process(trackParVarPos1, trackParVarNeg1) == 0) {
        continue;
      }
      const auto& secondaryVertex = fitter.getPCACandidate();
      auto trackParVar0 = fitter.getTrack(0);
      auto trackParVar1 = fitter.getTrack(1);

      // get track momenta
      std::array<float, 3> pVec0;
      std::array<float, 3> pVec1;
      trackParVar0.getPxPyPzGlo(pVec0);
      trackParVar1.getPxPyPzGlo(pVec1);

      // fill candidate table rows
      rowCandidateBase(collision.globalIndex(),
                       collision.posX(), collision.posY(), collision.posZ(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       pVec0[0], pVec0[1], pVec0[2],
                       pVec1[0], pVec1[1], pVec1[2],
                       rowTrackIndexProng2.prong0Id(), rowTrackIndexProng2.prong1Id());

      // fill histograms
      // calculate invariant masses
      auto arrayMomenta = std::array{pVec0, pVec1};
      massPiK = RecoDecay::m(arrayMomenta, std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus});
      massKPi = RecoDecay::m(arrayMomenta, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
      hMass->Fill(massPiK);
      // hMass->Fill(massKPi);
    }
  }
};

/// Helper extension task
/// Extends the base table with expression columns (see the HfTCand2ProngExt table).
struct HfTaskMiniCandidateCreator2ProngExpressions {
  Spawns<aod::HfTCand2ProngExt> rowCandidateProng2;
  void init(InitContext const&) {}
};

// Candidate selection =====================================================================

/// D0 candidate selector
struct HfTaskMiniCandidateSelectorD0 {
  Produces<aod::HfTSelD0> hfSelD0Candidate;

  Configurable<float> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // TPC
  Configurable<float> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<float> ptPidTpcMax{"ptPidTpcMax", 5., "Upper bound of track pT for TPC PID"};
  Configurable<float> nSigmaTpc{"nSigmaTpc", 3., "Nsigma cut on TPC only"};
  // topological cuts
  Configurable<float> cpaMin{"cpaMin", 0.98, "Min. cosine of pointing angle"};
  Configurable<float> massWindow{"massWindow", 0.4, "Half-width of the invariant-mass window"};

  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;

  using TracksWithPid = soa::Join<Tracks,
                                  aod::pidTPCFullPi, aod::pidTPCFullKa,
                                  aod::pidTOFFullPi, aod::pidTOFFullKa>;

  void init(InitContext const&)
  {
    selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpc, nSigmaTpc);
    selectorKaon = selectorPion;
  }

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T>
  bool selectionTopol(const T& candidate)
  {
    // check that the candidate pT is within the analysis range
    if (candidate.pt() < ptCandMin || candidate.pt() >= ptCandMax) {
      return false;
    }
    // cosine of pointing angle
    if (candidate.cpa() < cpaMin) {
      return false;
    }
    return true;
  }

  /// Conjugate-dependent topological cuts
  /// \param candidate candidate
  /// \param trackPion the track with the pion hypothesis
  /// \param trackKaon the track with the kaon hypothesis
  /// \note trackPion = positive and trackKaon = negative for D0 selection and inverse for D0bar
  /// \return true if candidate passes all cuts for the given conjugate
  template <typename T1, typename T2>
  bool selectionTopolConjugate(const T1& candidate, const T2& trackPion, const T2& /*trackKaon*/)
  {
    // invariant-mass cut
    if (trackPion.sign() > 0) {
      if (std::abs(HfHelper::invMassD0ToPiK(candidate) - o2::constants::physics::MassD0) > massWindow) {
        return false;
      }
    } else {
      if (std::abs(HfHelper::invMassD0barToKPi(candidate) - o2::constants::physics::MassD0) > massWindow) {
        return false;
      }
    }
    return true;
  }

  void process(aod::HfTCand2Prong const& candidates,
               TracksWithPid const&)
  {
    // looping over 2-prong candidates
    for (const auto& candidate : candidates) {

      // final selection flag: 0 - rejected, 1 - accepted
      int statusD0 = 0;
      int statusD0bar = 0;

      auto trackPos = candidate.prong0_as<TracksWithPid>(); // positive daughter
      auto trackNeg = candidate.prong1_as<TracksWithPid>(); // negative daughter

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelD0Candidate(statusD0, statusD0bar);
        continue;
      }

      // conjugate-dependent topological selection for D0
      bool topolD0 = selectionTopolConjugate(candidate, trackPos, trackNeg);
      // conjugate-dependent topological selection for D0bar
      bool topolD0bar = selectionTopolConjugate(candidate, trackNeg, trackPos);

      if (!topolD0 && !topolD0bar) {
        hfSelD0Candidate(statusD0, statusD0bar);
        continue;
      }

      // track-level PID selection
      int pidTrackPosKaon = selectorKaon.statusTpcOrTof(trackPos);
      int pidTrackPosPion = selectorPion.statusTpcOrTof(trackPos);
      int pidTrackNegKaon = selectorKaon.statusTpcOrTof(trackNeg);
      int pidTrackNegPion = selectorPion.statusTpcOrTof(trackNeg);

      int pidD0 = -1;
      int pidD0bar = -1;

      if (pidTrackPosPion == TrackSelectorPID::Accepted &&
          pidTrackNegKaon == TrackSelectorPID::Accepted) {
        pidD0 = 1; // accept D0
      } else if (pidTrackPosPion == TrackSelectorPID::Rejected ||
                 pidTrackNegKaon == TrackSelectorPID::Rejected ||
                 pidTrackNegPion == TrackSelectorPID::Accepted ||
                 pidTrackPosKaon == TrackSelectorPID::Accepted) {
        pidD0 = 0; // exclude D0
      }

      if (pidTrackNegPion == TrackSelectorPID::Accepted &&
          pidTrackPosKaon == TrackSelectorPID::Accepted) {
        pidD0bar = 1; // accept D0bar
      } else if (pidTrackPosPion == TrackSelectorPID::Accepted ||
                 pidTrackNegKaon == TrackSelectorPID::Accepted ||
                 pidTrackNegPion == TrackSelectorPID::Rejected ||
                 pidTrackPosKaon == TrackSelectorPID::Rejected) {
        pidD0bar = 0; // exclude D0bar
      }

      if (pidD0 == 0 && pidD0bar == 0) {
        hfSelD0Candidate(statusD0, statusD0bar);
        continue;
      }

      if ((pidD0 == -1 || pidD0 == 1) && topolD0) {
        statusD0 = 1; // identified as D0
      }
      if ((pidD0bar == -1 || pidD0bar == 1) && topolD0bar) {
        statusD0bar = 1; // identified as D0bar
      }

      hfSelD0Candidate(statusD0, statusD0bar);
    }
  }
};

// Analysis task =====================================================================

/// D0 analysis task
struct HfTaskMiniD0 {
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection flag for D0 bar"};

  Partition<soa::Join<aod::HfTCand2Prong, aod::HfTSelD0>> selectedD0Candidates = aod::hf_selcandidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_selcandidate_d0::isSelD0bar >= selectionFlagD0bar;

  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    const TString strTitle = "D^{0} candidates";
    const TString strPt = "#it{p}_{T} (GeV/#it{c})";
    const TString strEntries = "entries";
    registry.add("hPtCand", strTitle + ";" + strPt + ";" + strEntries, {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hMass", strTitle + ";" + "inv. mass (#pi K) (GeV/#it{c}^{2})" + ";" + strEntries, {HistType::kTH1F, {{500, 0., 5.}}});
    registry.add("hCpaVsPtCand", strTitle + ";" + "cosine of pointing angle" + ";" + strPt + ";" + strEntries, {HistType::kTH2F, {{110, -1.1, 1.1}, {100, 0., 10.}}});
  }

  void process(soa::Join<aod::HfTCand2Prong, aod::HfTSelD0> const& /*candidates*/)
  {
    for (const auto& candidate : selectedD0Candidates) {
      if (candidate.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMass"), HfHelper::invMassD0ToPiK(candidate));
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMass"), HfHelper::invMassD0barToKPi(candidate));
      }
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hCpaVsPtCand"), candidate.cpa(), candidate.pt());
    }
  }
};

// Add all tasks in the workflow specification.
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskMiniCandidateCreator2Prong>(cfgc, TaskName{"hf-task-mini-candidate-creator-2prong"}),                        // o2-linter: disable=name/o2-task (wrong hyphenation)
    adaptAnalysisTask<HfTaskMiniCandidateCreator2ProngExpressions>(cfgc, TaskName{"hf-task-mini-candidate-creator-2prong-expressions"}), // o2-linter: disable=name/o2-task (wrong hyphenation)
    adaptAnalysisTask<HfTaskMiniCandidateSelectorD0>(cfgc),
    adaptAnalysisTask<HfTaskMiniD0>(cfgc)};
}
