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

/// \file k892pmanalysis.cxx
/// \brief Reconstruction of track-V0 decay resonance candidates
///
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>, Alessandro Sturniolo <a.sturniolo@cern.ch>

#include <TLorentzVector.h>

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "DataFormatsParameters/GRPObject.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct k892pmanalysis {
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ///// Configurables
  /// Histograms
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Binning of the centrality axis"};
  Configurable<float> cInvMassStart{"cInvMassStart", 0.6, "Invariant mass start"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 1.5, "Invariant mass end"};
  Configurable<int> cInvMassBins{"cInvMassBins", 900, "Invariant mass binning"};
  Configurable<float> cK0shortMassStart{"cK0shortMassStart", 0.4, "K0Short mass start"};
  Configurable<float> cK0shortMassEnd{"cK0shortMassEnd", 0.6, "K0Short mass end"};
  Configurable<int> cK0shortMassBins{"cK0shortMassBins", 50, "K0Short mass binning"};
  Configurable<int> cDCABins{"cDCABins", 150, "DCA binning"};
  /// Pre-selection cuts
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minimum pt cut"};
  /// DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  /// PID Selections
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"};                // TPC
  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};                // TOF
  Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"};  // Combined
  Configurable<bool> cUseOnlyTOFTrackPi{"cUseOnlyTOFTrackPi", false, "Use only TOF track for PID selection"}; // Use only TOF track for Pion PID selection
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection
  // Track selections
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};           // PV Contributor
  // V0 selections
  Configurable<double> cV0MinCosPA{"cV0MinCosPA", 0.97, "V0 minimum pointing angle cosine"};
  Configurable<double> cV0MaxDaughDCA{"cV0MaxDaughDCA", 1.0, "V0 daughter DCA Maximum"};

  void init(o2::framework::InitContext&)
  {
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec dcaxyAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec k0sMassAxis = {cK0shortMassBins, cK0shortMassStart, cK0shortMassEnd, "K^{0}_{S} Mass (GeV/#it{c}^2)"};

    // Mass QA (quick check)
    histos.add("k892pminvmass", "Invariant mass of charged K*(892)", kTH1F, {invMassAxis});
    histos.add("QAafter/k0shortmass", "Invariant mass of K0Short", kTH1F, {k0sMassAxis});
    // DCA QA
    histos.add("QAbefore/collMult", "Collision multiplicity", HistType::kTH1F, {binsCent});
    histos.add("QAbefore/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAbefore/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
    // histos.add("QAbefore/trkDCAz_k0s", "DCAz distribution of k0short track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAafter/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAafter/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
    // histos.add("QAafter/trkDCAz_k0s", "DCAz distribution of k0short track candidates", HistType::kTH1F, {dcazAxis});
    //  pT QA
    histos.add("QAbefore/trkpT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxisQA});
    histos.add("QAbefore/trkpT_k0s", "pT distribution of k0short track candidates", kTH1F, {ptAxisQA});
    histos.add("QAafter/trkpT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxisQA});
    histos.add("QAafter/trkpT_k0s", "pT distribution of k0short track candidates", kTH1F, {ptAxisQA});
    // Mass vs Pt vs Multiplicity 3-dimensional histogram
    histos.add("k892pmMassPtMult3d", "Charged K*(892) mass vs pT vs V0multiplicity distribution", kTH3F, {invMassAxis, ptAxis, centAxis});

    if (doprocessMCLight) {
      // MC QA
      histos.add("k892pmPtGen", "pT distribution of True MC charged K*(892)", kTH1F, {ptAxis});
      histos.add("k892pmPtRec", "pT distribution of Reconstructed MC charged K*(892)", kTH1F, {ptAxis});
    }
    // Print output histograms statistics
    LOG(info) << "Size of the histograms in spectraTOF";
    histos.print();
  }

  double massK0 = MassK0Short;
  double massPi = MassPionCharged;

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;

    return true;
  }

  template <typename V0Type>
  bool V0Cut(const V0Type v0)
  {
    // V0 track cuts
    if (std::abs(v0.eta()) > 0.8)
      return false;
    if (v0.v0CosPA() < cV0MinCosPA)
      return false;
    if (v0.daughDCA() > cV0MaxDaughDCA)
      return false;

    return true;
  }

  // Primary PID selection tools
  template <typename T>
  bool selectionPIDPrimaryPion(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPion) {
      tpcPIDPassed = true;
    }
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPi()) < cMaxTOFnSigmaPion) {
        tofPIDPassed = true;
      }
      if ((nsigmaCutCombinedPion > 0) && (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi() + candidate.tofNSigmaPi() * candidate.tofNSigmaPi() < nsigmaCutCombinedPion * nsigmaCutCombinedPion)) {
        tofPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType, typename V0sType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks, const V0sType& dV0s)
  {
    // auto multiplicity = collision.cent();
    auto multiplicity = collision.cent();
    histos.fill(HIST("QAbefore/collMult"), multiplicity);
    TLorentzVector lDecayDaughter, lDecayV0, lResonance;

    bool IsV0Processed = false;
    bool IsV0QAFilled = false;

    for (auto& trk : dTracks) {
      // Full index policy is needed to consider all possible combinations
      //// Initialize variables
      // trk: Pion, v0: K0s

      auto trkId = trk.index();
      auto trkptPi = trk.pt();

      // DCA QA (before cuts)
      histos.fill(HIST("QAbefore/trkDCAxy_pi"), trk.dcaXY());
      histos.fill(HIST("QAbefore/trkDCAz_pi"), trk.dcaZ());
      // pT QA (before cuts)
      histos.fill(HIST("QAbefore/trkpT_pi"), trkptPi);

      // apply the track cut
      if (!trackCut(trk) || !selectionPIDPrimaryPion(trk))
        continue;

      // DCA QA (QAafter cuts)
      histos.fill(HIST("QAafter/trkDCAxy_pi"), trk.dcaXY());
      histos.fill(HIST("QAafter/trkDCAz_pi"), trk.dcaZ());
      // pT QA (after cuts)
      histos.fill(HIST("QAafter/trkpT_pi"), trk.pt());

      for (auto& v0 : dV0s) {
        // Full index policy is needed to consider all possible combinations
        if (v0.indices()[0] == trkId || v0.indices()[1] == trkId)
          continue; // To avoid combining secondary and primary pions
        //// Initialize variables
        // trk: Pion, v0: K0s

        auto v0ptK0s = v0.pt();

        if (!IsV0QAFilled) {
          // pT QA (before cuts)
          histos.fill(HIST("QAbefore/trkpT_k0s"), v0ptK0s);
        }

        // apply the track cut
        if (!V0Cut(v0))
          continue;

        if (!IsV0QAFilled) {
          // pt QA (after cuts)
          histos.fill(HIST("QAafter/trkpT_k0s"), v0ptK0s);
          // K0s mass QA (after cuts)
          histos.fill(HIST("QAafter/k0shortmass"), v0.mK0Short());
        }

        lDecayDaughter.SetXYZM(trk.px(), trk.py(), trk.pz(), massPi);
        lDecayV0.SetXYZM(v0.px(), v0.py(), v0.pz(), massK0);
        lResonance = lDecayDaughter + lDecayV0;
        // Filling invariant mass histograms
        // K*(892)pm mass
        histos.fill(HIST("k892pminvmass"), lResonance.M());
        // K*(892)pm 3d mass, pt, multiplicity histogram
        histos.fill(HIST("k892pmMassPtMult3d"), lResonance.M(), lResonance.Pt(), multiplicity);
        if constexpr (IsMC) {
          if (abs(trk.pdgCode()) != 211 || abs(v0.pdgCode()) != 310) // Skip to next iteration if duaghters are not charged pion + K0s/AntiK0s
            continue;
          if (trk.motherPDG() != v0.motherPDG())
            continue;
          if (trk.motherPDG() != 323)
            continue;
          histos.fill(HIST("k892pmPtRec"), lResonance.Pt());
        }
        IsV0Processed = true;
      }
      if (IsV0Processed) {
        IsV0QAFilled = true;
      }
    }
  }

  void processDataLight(aod::ResoCollision& collision,
                        aod::ResoTracks const& resotracks,
                        aod::ResoV0s const& resov0s)
  {
    // LOG(info) << "new collision, zvtx: " << collision.posZ();
    fillHistograms<false, false>(collision, resotracks, resov0s);
  }
  PROCESS_SWITCH(k892pmanalysis, processDataLight, "Process Event for data", false);

  void processMCLight(aod::ResoCollision& collision,
                      soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks,
                      soa::Join<aod::ResoV0s, aod::ResoMCV0s> const& resov0s)
  {
    fillHistograms<true, false>(collision, resotracks, resov0s);
  }
  PROCESS_SWITCH(k892pmanalysis, processMCLight, "Process Event for MC", false);

  void processMCTrue(aod::ResoMCParents& resoParents)
  {
    for (auto& part : resoParents) {  // loop over all pre-filtered MC particles
      if (abs(part.pdgCode()) != 323) // K*892(pm)
        continue;
      if (abs(part.y()) > 0.5) // rapidity cut
        continue;
      bool pass1 = false;
      bool pass2 = false;
      if (part.daughterPDG1() == 211 && part.daughterPDG2() == 310) { // One decay to K0s and the other to pi+ (K*(892)+ mother) - Particle pass
        pass1 = true;
      }
      if (part.daughterPDG1() == -211 && part.daughterPDG2() == -310) { // One decay to AntiK0s and the other to pi- (K*(892)- mother) - Antiparticle pass
        pass2 = true;
      }
      if (!pass1 || !pass2) // Go on only if we have both decay products, else skip to next iteration
        continue;
      histos.fill(HIST("k892pmPtGen"), part.pt());
    }
  }
  PROCESS_SWITCH(k892pmanalysis, processMCTrue, "Process Event for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<k892pmanalysis>(cfgc, TaskName{"lf-k892pmanalysis"})};
}
