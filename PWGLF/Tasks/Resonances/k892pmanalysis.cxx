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

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <TLorentzVector.h>

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
  ConfigurableAxis binsEtaQA{"binsEtaQA", {VARIABLE_WIDTH, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, "Binning of pseudo-rapidity eta"};
  ConfigurableAxis binsV0CosPointAngleQA{"binsV0CosPointAngleQA", {VARIABLE_WIDTH, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0}, "Binning of V0 cosine of pointing angle"};

  Configurable<float> cInvMassStart{"cInvMassStart", 0.6, "Invariant mass start"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 1.5, "Invariant mass end"};
  Configurable<int> cInvMassBins{"cInvMassBins", 900, "Invariant mass binning"};

  Configurable<float> cK0shortMassStart{"cK0shortMassStart", 0.4, "K0Short mass start"};
  Configurable<float> cK0shortMassEnd{"cK0shortMassEnd", 0.6, "K0Short mass end"};
  Configurable<int> cK0shortMassBins{"cK0shortMassBins", 50, "K0Short mass binning"};

  Configurable<float> cLambdaAntiLambdaMassStart{"cLambdaAntiLambdaMassStart", 1.0, "V0 mass (in the (Anti)Lambda0 hypothesis) start"};
  Configurable<float> cLambdaAntiLambdaMassEnd{"cLambdaAntiLambdaMassEnd", 2.0, "V0 mass (in the (Anti)Lambda0 hypothesis) end"};
  Configurable<int> cLambdaAntiLambdaMassBins{"cLambdaAntiLambdaMassBins", 250, "V0 mass (in the (Anti)Lambda0 hypothesis) binning"};

  Configurable<int> cDCABinsQA{"cDCABinsQA", 150, "DCA binning"};

  Configurable<int> cTpcNsigmaPionBinsQA{"cTpcNsigmaPionBinsQA", 140, "tpcNSigmaPi binning"};

  // Pre-selection cuts
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minimum pt cut"};
  Configurable<double> cMaxEtacut{"cMaxEtacut", 0.8, "Track maximum eta cut"};
  Configurable<double> cMaxV0Etacut{"cMaxV0Etacut", 0.8, "V0 maximum eta cut"};
  // DCA Selections
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
  // Competing V0 rejection
  Configurable<double> cV0MassWindow{"cV0MassWindow", 0.0043, "Mass window for competing Lambda0 rejection"};
  // Resonance selection
  // Configurable<double> cMaxResRapidity{"cMaxResRapidity", 0.5, "Maximum pseudo-rapidity value of reconstructed K*(892)pm resonance"};
  // Event mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Mixing bins - multiplicity"};

  void init(o2::framework::InitContext&)
  {
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec etaAxisQA = {binsEtaQA, "#eta"};
    AxisSpec dcaxyAxisQA = {cDCABinsQA, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxisQA = {cDCABinsQA, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec daughdcaAxisQa = {cDCABinsQA, 0.0, 3.0, "V0 daughters DCA (cm)"};
    AxisSpec CosPointAngleAxisQA = {binsV0CosPointAngleQA, "CosPA"};
    AxisSpec tpcNSigmaPiAxisQA = {cTpcNsigmaPionBinsQA, -7.0, 7.0, "N#sigma_{TPC}"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec k0sMassAxisQA = {cK0shortMassBins, cK0shortMassStart, cK0shortMassEnd, "K^{0}_{S} Mass (GeV/#it{c}^2)"};
    AxisSpec lambdaAntilambdaMassAxisQA = {cLambdaAntiLambdaMassBins, cLambdaAntiLambdaMassStart, cLambdaAntiLambdaMassEnd, "(Anti-)#Lambda^{0} Mass (GeV/#it{c}^2)"};
    AxisSpec k892pmCountAxis = {2, 0., 2., "K*^{+}(892) = 1, K*^{-}(892) = 2"};
    AxisSpec goodTrackCountAxis = {3, 0., 3., "Passed track = 1, Passed V0 = 2, Passed track and V0 = 3"};

    // Mass QA (quick check)
    // QA before
    histos.add("QAbefore/k0shortmassPt", "Invariant mass of K0Short vs K0Short Pt", kTH2F, {ptAxisQA, k0sMassAxisQA});
    histos.add("QAbefore/lambda0mass", "Invariant mass of V0 in Lambda0 hypothesis", kTH1F, {lambdaAntilambdaMassAxisQA});
    histos.add("QAbefore/antilambda0mass", "Invariant mass of V0 in AntiLambda0 hypothesis", kTH1F, {lambdaAntilambdaMassAxisQA});
    // QA after
    histos.add("QAafter/k0shortmassPt", "Invariant mass of K0Short vs K0Short Pt", kTH2F, {ptAxisQA, k0sMassAxisQA});
    histos.add("QAafter/lambda0mass", "Invariant mass of V0 in Lambda0 hypothesis", kTH1F, {lambdaAntilambdaMassAxisQA});
    histos.add("QAafter/antilambda0mass", "Invariant mass of V0 in AntiLambda0 hypothesis", kTH1F, {lambdaAntilambdaMassAxisQA});
    histos.add("k892pminvmass", "Invariant mass of charged K*(892)", kTH1F, {invMassAxis});
    // Multiplicity and accepted events QA
    histos.add("QAbefore/collMult", "Collision multiplicity", HistType::kTH1F, {centAxis});
    // DCA QA
    // QA before
    histos.add("QAbefore/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxisQA});
    histos.add("QAbefore/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxisQA});
    histos.add("QAbefore/trkV0DaughDCA", "DCA distribution of V0 daughters", kTH1F, {daughdcaAxisQa});
    // QA after
    histos.add("QAafter/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxisQA});
    histos.add("QAafter/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxisQA});
    histos.add("QAafter/trkV0DaughDCA", "DCA distribution of V0 daughters", kTH1F, {daughdcaAxisQa});
    // Pseudo-rapidity QA
    // QA before
    histos.add("QAbefore/trkEta_pi", "Primary pion track pseudo-rapidity eta", kTH1F, {etaAxisQA});
    histos.add("QAbefore/trkEta_k0s", "K0short track pseudo-rapidity eta", kTH1F, {etaAxisQA});
    histos.add("QAbefore/k892pmRapidity", "Reconstructed K*(892)^{#pm} rapidity", kTH1F, {etaAxisQA});
    // QA after
    histos.add("QAafter/trkEta_pi", "Primary pion track pseudo-rapidity eta", kTH1F, {etaAxisQA});
    histos.add("QAafter/trkEta_k0s", "K0short track pseudo-rapidity eta", kTH1F, {etaAxisQA});
    histos.add("QAafter/k892pmRapidity", "Reconstructed K*(892)^{#pm} rapidity", kTH1F, {etaAxisQA});
    // Cosine pointing angle QA
    // QA before
    histos.add("QAbefore/V0CosPA", "V0 cosine of pointing angle", kTH1F, {CosPointAngleAxisQA});
    // QA after
    histos.add("QAafter/V0CosPA", "V0 cosine of pointing angle", kTH1F, {CosPointAngleAxisQA});
    // pT QA
    // QA before
    histos.add("QAbefore/trkpT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxisQA});
    histos.add("QAbefore/trkpT_k0s", "pT distribution of k0short track candidates", kTH1F, {ptAxisQA});
    // QA after
    histos.add("QAafter/trkpT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxisQA});
    histos.add("QAafter/trkpT_k0s", "pT distribution of k0short track candidates", kTH1F, {ptAxisQA});
    // Primary pion TPC PID
    // QA before
    histos.add("QAbefore/tpcNsigmaPionQA", "NsigmaTPC distribution of primary pion candidates", kTH2F, {ptAxisQA, tpcNSigmaPiAxisQA});
    // QA after
    histos.add("QAafter/tpcNsigmaPionQA", "NsigmaTPC distribution of primary pion candidates", kTH2F, {ptAxisQA, tpcNSigmaPiAxisQA});
    // Good tracks and V0 counts QA
    histos.add("QAafter/hGoodTracksV0s", "Number of good track and V0 passed", kTH1F, {goodTrackCountAxis});
    // Mass vs Pt vs Multiplicity 3-dimensional histogram
    histos.add("k892pmMassPtMult3d", "Charged K*(892) mass vs pT vs V0 multiplicity distribution", kTH3F, {invMassAxis, ptAxis, centAxis});

    if (doprocessMELight) {
      // Event-mixing background
      histos.add("k892pmMEbackground", "Mixed-event reconstructed K*(892) invariant mass distribution", kTH1F, {invMassAxis});
      histos.add("k892pmMEMassPtMult3d", "Mixed-event reconstructed K*(892) mass vs pT vs V0 multiplicity distribution", kTH3F, {invMassAxis, ptAxis, centAxis});
    }

    if (doprocessMCLight) {
      // MC QA
      histos.add("k892pmPtGen", "pT distribution of True MC charged K*(892)", kTH1F, {ptAxis});
      histos.add("k892pmPtRec", "pT distribution of Reconstructed MC charged K*(892)", kTH1F, {ptAxis});
    }

    if (doprocessMCTrue) {
      // DEBUG HISTOGRAMS
      histos.add("hK892pmCounter", "Generated MC resonances", kTH1F, {k892pmCountAxis});
      // histos.add("hDaughterCounter", "Generated MC resonance daughters", kTH1F, {daughterCountAxis});
    }
    // Print output histograms statistics
    LOG(info) << "Size of the histograms in spectraTOF";
    histos.print();
  }

  double massK0 = MassK0Short;
  double massPi = MassPionCharged;
  double massLambda0 = MassLambda;
  double massAntiLambda0 = MassLambda0Bar;

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.eta()) > cMaxEtacut)
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
    if (std::abs(v0.eta()) > cMaxV0Etacut)
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

  /*// Resonance candidate selection
  template <typename ResoCandidate>
  bool selectionResoK892pm(const ResoCandidate& resoCandidate)
  {
    // Rapidity cut
    if (resoCandidate.Rapidity() > cMaxResRapidity) {
      continue;
    }
    return true;
  }*/

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
      auto trketaPi = trk.eta();
      // LOG(INFO) << "Primary Pion Eta before cuts is: " << trketaPi;

      // Pseudo-rapidity QA (before cuts)
      histos.fill(HIST("QAbefore/trkEta_pi"), trketaPi);

      if (!IsMix) {
        // DCA QA (before cuts)
        histos.fill(HIST("QAbefore/trkDCAxy_pi"), trk.dcaXY());
        histos.fill(HIST("QAbefore/trkDCAz_pi"), trk.dcaZ());
        // Pseudo-rapidity QA (before cuts)
        histos.fill(HIST("QAbefore/trkEta_pi"), trketaPi);
        // pT QA (before cuts)
        histos.fill(HIST("QAbefore/trkpT_pi"), trkptPi);
        // TPC PID (before cuts)
        histos.fill(HIST("QAbefore/tpcNsigmaPionQA"), trkptPi, trk.tpcNSigmaPi());
      }

      // apply the track cut
      if (!trackCut(trk) || !selectionPIDPrimaryPion(trk))
        continue;

      histos.fill(HIST("QAafter/hGoodTracksV0s"), 0.5);

      // LOG(INFO) << "Primary Pion Eta after cuts is: " << trketaPi;

      if (!IsMix) {
        // DCA QA (QAafter cuts)
        histos.fill(HIST("QAafter/trkDCAxy_pi"), trk.dcaXY());
        histos.fill(HIST("QAafter/trkDCAz_pi"), trk.dcaZ());
        // Pseudo-rapidity QA (after cuts)
        histos.fill(HIST("QAafter/trkEta_pi"), trketaPi);
        // pT QA (after cuts)
        histos.fill(HIST("QAafter/trkpT_pi"), trk.pt());
        // TPC PID (after cuts)
        histos.fill(HIST("QAafter/tpcNsigmaPionQA"), trkptPi, trk.tpcNSigmaPi());
      }

      for (auto& v0 : dV0s) {
        // Full index policy is needed to consider all possible combinations
        if (v0.indices()[0] == trkId || v0.indices()[1] == trkId)
          continue; // To avoid combining secondary and primary pions
        //// Initialize variables
        // trk: Pion, v0: K0s

        auto v0ptK0s = v0.pt();

        if (!IsMix && !IsV0QAFilled) {
          // V0 daughter DCA (before cuts)
          histos.fill(HIST("QAbefore/trkV0DaughDCA"), v0.daughDCA());
          // Pseudo-rapidity QA (before cuts)
          histos.fill(HIST("QAbefore/trkEta_k0s"), v0.eta());
          // CosPA (before cuts)
          histos.fill(HIST("QAbefore/V0CosPA"), v0.v0CosPA());
          // pT QA (before cuts)
          histos.fill(HIST("QAbefore/trkpT_k0s"), v0ptK0s);
          // K0s mass QA (before cuts)
          histos.fill(HIST("QAbefore/k0shortmassPt"), v0ptK0s, v0.mK0Short());
        }

        // apply the track cut
        if (!V0Cut(v0))
          continue;

        histos.fill(HIST("QAafter/hGoodTracksV0s"), 1.5);

        // QA (anti)Lamda0 mass before competing V0 rejection cut
        if (!IsMix && !IsV0QAFilled) {
          histos.fill(HIST("QAbefore/lambda0mass"), v0.mLambda());
          histos.fill(HIST("QAbefore/antilambda0mass"), v0.mAntiLambda());
        }

        // apply the competing V0 rejection cut (excluding Lambda0 candidates, massLambdaPDG = 1115.683 MeV/c2)
        if (std::abs(v0.mLambda() - massLambda0) < cV0MassWindow)
          continue;
        if (std::abs(v0.mAntiLambda() - massAntiLambda0) < cV0MassWindow)
          continue;

        if (!IsMix && !IsV0QAFilled) {
          // V0 daughter DCA (after cuts)
          histos.fill(HIST("QAafter/trkV0DaughDCA"), v0.daughDCA());
          // Pseudo-rapidity QA (before cuts)
          histos.fill(HIST("QAafter/trkEta_k0s"), v0.eta());
          // CosPA (after cuts)
          histos.fill(HIST("QAafter/V0CosPA"), v0.v0CosPA());
          // pt QA (after cuts)
          histos.fill(HIST("QAafter/trkpT_k0s"), v0ptK0s);
          // K0s mass QA (after cuts)
          histos.fill(HIST("QAafter/k0shortmassPt"), v0ptK0s, v0.mK0Short());
          histos.fill(HIST("QAafter/lambda0mass"), v0.mLambda());
          histos.fill(HIST("QAafter/antilambda0mass"), v0.mAntiLambda());
        }

        lDecayDaughter.SetXYZM(trk.px(), trk.py(), trk.pz(), massPi);
        lDecayV0.SetXYZM(v0.px(), v0.py(), v0.pz(), massK0);
        lResonance = lDecayDaughter + lDecayV0;

        if (!IsMix && !IsV0QAFilled) {
          // Pseudo-rapidity QA (before eta cut)
          histos.fill(HIST("QAbefore/k892pmRapidity"), lResonance.Rapidity());
        }

        // Checking whether the mid-rapidity condition is met
        if (std::abs(lResonance.Rapidity()) > 0.5) {
          continue;
        }

        if (!IsMix && !IsV0QAFilled) {
          // Pseudo-rapidity QA (before eta cut)
          histos.fill(HIST("QAafter/k892pmRapidity"), lResonance.Rapidity());
        }

        // Counting how many resonances passed
        histos.fill(HIST("QAafter/hGoodTracksV0s"), 2.5);
        // Filling invariant mass histograms
        if constexpr (!IsMix) {
          // Reconstructed K*(892)pm mass
          histos.fill(HIST("k892pminvmass"), lResonance.M());
          // Reconstructed K*(892)pm 3d mass, pt, multiplicity histogram
          histos.fill(HIST("k892pmMassPtMult3d"), lResonance.M(), lResonance.Pt(), multiplicity);
          if constexpr (IsMC) {
            // LOG(info) << "track PDG:\t" << trk.pdgCode() << "\tV0 PDG:\t" << v0.pdgCode();
            if (std::abs(trk.pdgCode()) != 211 || std::abs(v0.pdgCode()) != 310) // Skip to next iteration if daughters are not charged pion + K0s/AntiK0s
              continue;
            if (trk.motherPDG() != v0.motherPDG())
              continue;
            // LOG(info) << "track PDG:\t" << trk.pdgCode() << "\tV0 PDG:\t" << v0.pdgCode();
            if (trk.motherPDG() != 323)
              continue;
            histos.fill(HIST("k892pmPtRec"), lResonance.Pt());
          }
        } else {
          // Mixed-event reeconstructed K*(892)pm mass
          histos.fill(HIST("k892pmMEbackground"), lResonance.M());
          // Mixed-event reconstructed K*(892)pm 3d mass, pt, multiplicity histogram
          histos.fill(HIST("k892pmMEMassPtMult3d"), lResonance.M(), lResonance.Pt(), multiplicity);
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
    for (auto& part : resoParents) {       // loop over all pre-filtered MC particles
      if (std::abs(part.pdgCode()) != 323) // K*892(pm)
        continue;
      if (std::abs(part.y()) > 0.5) // rapidity cut
        continue;
      bool pass1 = false;
      bool pass2 = false;
      /*// Sanity check: looking for K*0 resonances for sanity check
      if (std::abs(part.pdgCode()) == 323) {
        LOG(info) << "Found charged K*: " << part.pdgCode() << ". Daughters' PDG are " << part.daughterPDG1() << " and " << part.daughterPDG2();
      }
      if (std::abs(part.pdgCode()) == 313) {
        LOG(info) << "Found non-charged K*: " << part.pdgCode() << ". Daughters' PDG are " << part.daughterPDG1() << " and " << part.daughterPDG2();
      }*/

      if (part.daughterPDG1() == 211 && part.daughterPDG2() == 310) { // One decay to K0s and the other to pi+ (K*(892)+ mother) - Particle pass
        pass1 = true;
        histos.fill(HIST("hK892pmCounter"), 0.5);
      }
      if (part.daughterPDG1() == -211 && part.daughterPDG2() == -310) { // One decay to AntiK0s and the other to pi- (K*(892)- mother) - Antiparticle pass
        pass2 = true;
        histos.fill(HIST("hK892pmCounter"), 1.5);
      }
      /*if (std::abs(part.daughterPDG1()) == 211)
        histos.fill(HIST("hDaughterCounter"), 0.5);
      if (std::abs(part.daughterPDG2()) == 310)
        histos.fill(HIST("hDaughterCounter"), 1.5);
      if (std::abs(part.daughterPDG1()) == 211 && std::abs(part.daughterPDG2()) == 310)
        histos.fill(HIST("hDaughterCounter"), 2.5);*/
      // if (!pass1 || !pass2) // Go on only if we have both decay products, else skip to next iteration
      if (!pass1 && !pass2) // Go on only if we have both decay products, else skip to next iteration
        continue;
      histos.fill(HIST("k892pmPtGen"), part.pt());
    }
  }
  PROCESS_SWITCH(k892pmanalysis, processMCTrue, "Process Event for MC", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMELight(o2::aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks, aod::ResoV0s const& resov0s)
  {
    auto tracksV0sTuple = std::make_tuple(resotracks, resov0s);
    // auto V0sTuple = std::make_tuple(resov0s);
    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    Pair<aod::ResoCollisions, aod::ResoTracks, aod::ResoV0s, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksV0sTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, resotracks1, collision2, resov0s2] : pairs) {
      fillHistograms<false, true>(collision1, resotracks1, resov0s2);
    }
  };
  PROCESS_SWITCH(k892pmanalysis, processMELight, "Process Event Mixing light without partition", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<k892pmanalysis>(cfgc, TaskName{"lf-k892pmanalysis"})};
}
