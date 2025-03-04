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

/// \file chargedkstaranalysis.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
///
/// \author Protay

#include "TF1.h"
// #include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
// #include <TPDGCode.h>
#include <string>

#include <array>
#include <cmath>
#include <cstdlib>

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"

// For charged kstarpp analysis
#include "PWGLF/DataModel/LFResonanceTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

struct chargedkstaranalysis {

  // Connect to ccdb
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  // For charged Kstarpp analysis use Resonance Initalizer and THnSparse
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Binning of the centrality axis"};
  ConfigurableAxis binsPt{"binsPt",
                          {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0},
                          "Binning of the pT axis"};
  ConfigurableAxis etabins{"etabins",
                           {VARIABLE_WIDTH, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0},
                           "Eta Binning"};
  Configurable<int> cDCABinsQA{"cDCABinsQA", 150, "DCA binning"};
  ConfigurableAxis binsPtQA{"binsPtQA",
                            {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0},
                            "Binning of the pT axis"};

  HistogramRegistry histos1{
    "histos1",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Pre-selection cuts
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minimum pt cut"};
  Configurable<bool> confevtcollintimerangestandard{"confevtcollintimerangestandard", true, "Evt sel: apply NoCollInTimeRangeStandard"};
  /// PID Selections
  Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999,
                                             "Combined nSigma cut for Pion"}; // Combined

  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5,
                                       "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0,
                                       "Track DCAz cut to PV Maximum"};
  // Track selections
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true,
                                     "Primary track selection"}; // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true,
                                         "Global track selection without DCA"}; // kQualityTracks (kTrackType |
                                                                                // kTPCNCls | kTPCCrossedRows |
                                                                                // kTPCCrossedRowsOverNCls |
                                                                                // kTPCChi2NDF | kTPCRefit |
                                                                                // kITSNCls | kITSChi2NDF |
                                                                                // kITSRefit | kITSHits) |
                                                                                // kInAcceptanceTracks (kPtRange |
                                                                                // kEtaRange)
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true,
                                      "PV contributor track selection"}; // PV Contributor
  // V0 selections
  Configurable<double> cV0MinCosPA{"cV0MinCosPA", 0.97,
                                   "V0 minimum pointing angle cosine"};
  Configurable<double> cV0MaxDaughDCA{"cV0MaxDaughDCA", 1.0,
                                      "V0 daughter DCA Maximum"};
  // Competing V0 rejection
  Configurable<double> cV0MassWindow{"cV0MassWindow", 0.0043, "Mass window for competing Lambda0 rejection"};
  Configurable<float> cInvMassStart{"cInvMassStart", 0.6, "Invariant mass start"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 1.5, "Invariant mass end"};
  Configurable<int> cInvMassBins{"cInvMassBins", 900, "Invariant mass binning"};

  // Event mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis cfgvtxbins{"cfgvtxbins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgmultbins{"cfgmultbins", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Mixing bins - multiplicity"};
  Configurable<int> cTpcNsigmaPionBinsQA{"cTpcNsigmaPionBinsQA", 140, "tpcNSigmaPi binning"};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurable parameters for V0 selection
  Configurable<float> confdaugheta{"confdaugheta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> nSigmaCutTPC{"nSigmaCutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nSigmaCutTOF{"nSigmaCutTOF", 3.0,
                                   "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the Combined Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  void init(InitContext const&)
  {
    AxisSpec dcaxyAxisQA = {cDCABinsQA, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxisQA = {cDCABinsQA, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec tpcNSigmaPiAxisQA = {cTpcNsigmaPionBinsQA, -7.0, 7.0,
                                  "N#sigma_{TPC}"};
    AxisSpec pidQAAxis = {130, -6.5, 6.5};
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd,
                            "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec etaAxis = {etabins, "#eta"};
    AxisSpec goodTrackCountAxis = {
      3, 0., 3., "Passed track = 1, Passed V0 = 2, Passed track and V0 = 3"};
    // register histograms
    histos1.add("hVertexZ", "hVertexZ", HistType::kTH1F, {{nBins, -15., 15.}});
    // Multiplicity and accepted events QA
    histos1.add("QAbefore/collMult", "Collision multiplicity", HistType::kTH1F,
                {centAxis});
    // QA before
    histos1.add("QAbefore/pi_Eta", "Primary pion track eta", kTH1F, {etaAxis});
    histos1.add("QAbefore/k0s_Eta", "K0short track eta", kTH1F, {etaAxis});
    histos1.add("QAbefore/chargedkstarpmRapidity",
                "Reconstructed K*^{#pm} rapidity", kTH1F, {etaAxis});
    histos1.add("QAbefore/trkpionTOFPID", "TOF PID of bachelor pion candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos1.add("QAbefore/trkpionTPCTOFPID", "TPC-TOF PID map of bachelor pion candidates", HistType::kTH2D, {pidQAAxis, pidQAAxis});

    histos1.add("QAbefore/DCAxy_pi",
                "DCAxy distribution of pion track candidates", HistType::kTH1F,
                {dcaxyAxisQA});
    histos1.add("QAbefore/DCAz_pi",
                "DCAz distribution of pion track candidates", HistType::kTH1F,
                {dcazAxisQA});
    histos1.add("QAbefore/pT_pi", "pT distribution of pion track candidates",
                kTH1F, {ptAxisQA});
    histos1.add("QAbefore/tpcNsigmaPionQA",
                "NsigmaTPC distribution of primary pion candidates", kTH2F,
                {ptAxisQA, tpcNSigmaPiAxisQA});

    // QA after
    histos1.add("QAAfter/DCAxy_pi",
                "DCAxy distribution of pion track candidates", HistType::kTH1F,
                {dcaxyAxisQA});
    histos1.add("QAAfter/DCAz_pi", "DCAz distribution of pion track candidates",
                HistType::kTH1F, {dcazAxisQA});
    histos1.add("QAAfter/pT_pi", "pT distribution of pion track candidates",
                kTH1F, {ptAxisQA});
    histos1.add("QAAfter/tpcNsigmaPionQA",
                "NsigmaTPC distribution of primary pion candidates", kTH2F,
                {ptAxisQA, tpcNSigmaPiAxisQA});
    histos1.add("QAAfter/pi_Eta", "Primary pion track eta", kTH1F, {etaAxis});

    // Good tracks and V0 counts QA
    histos1.add("QAafter/hGoodTracksV0s", "Number of good track and V0 passed",
                kTH1F, {goodTrackCountAxis});
    histos1.add("chargedkstarinvmassUlikeSign",
                "Invariant mass of charged K*(892)", kTH1F, {invMassAxis});
    histos1.add("chargedkstarinvmassMixedEvent",
                "Invariant mass of charged K*(892)", kTH1F, {invMassAxis});

    // Mass vs Pt vs Multiplicity 3-dimensional histogram
    //    histos1.add("chargekstarMassPtMult", "Charged K*(892) mass vs pT vs V0
    //    multiplicity distribution", kTH3F, {invMassAxis, ptAxis, centAxis});

    histos1.add("chargekstarMassPtMultPtUnlikeSign",
                "Invariant mass of CKS meson Unlike Sign", kTHnSparseF,
                {invMassAxis, ptAxis, centAxis}, true);
    histos1.add("chargekstarMassPtMultPtMixedEvent",
                "Invariant mass of CKS meson MixedEvent Sign", kTHnSparseF,
                {invMassAxis, ptAxis, centAxis}, true);
  }
  double massPi = o2::constants::physics::MassPionCharged;
  double massK0s = o2::constants::physics::MassK0Short;
  double massKa = o2::constants::physics::MassKPlus;
  ROOT::Math::PtEtaPhiMVector cksvector;

  double massK0 = o2::constants::physics::MassK0Short;
  double massPicharged = o2::constants::physics::MassPionCharged;
  double massLambda0 = o2::constants::physics::MassLambda;
  double massAntiLambda0 = o2::constants::physics::MassLambda0Bar;
  // Fill histograms (main function)
  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType,
            typename V0sType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks,
                      const V0sType& dV0s)
  {
    // auto multiplicity = collision.cent();
    auto multiplicity = collision.cent();
    histos1.fill(HIST("QAbefore/collMult"), multiplicity);
    TLorentzVector lDecayDaughter, lDecayV0, lResonance;

    for (const auto& track : dTracks) { // loop over all dTracks1 to find the bachelor pion
      auto trackId = track.index();
      auto trackptPi = track.pt();
      auto tracketaPi = track.eta();
      auto istrkhasTOF = track.hasTOF();
      auto trkNSigmaPiTPC = track.tpcNSigmaPi();
      auto trkNSigmaPiTOF = (istrkhasTOF) ? track.tofNSigmaPi() : -999.;

      histos1.fill(HIST("QAbefore/pi_Eta"), tracketaPi);

      if (!IsMix) {
        // TPC PID (before cuts)
        histos1.fill(HIST("QAbefore/tpcNsigmaPionQA"), trackptPi, trkNSigmaPiTPC);
        if (istrkhasTOF) {
          histos1.fill(HIST("QAbefore/trkpionTOFPID"), trackptPi, trkNSigmaPiTOF);
          histos1.fill(HIST("QAbefore/trkpionTPCTOFPID"), trkNSigmaPiTPC, trkNSigmaPiTOF);
        }
        // DCA QA (before cuts)
        histos1.fill(HIST("QAbefore/DCAxy_pi"), track.dcaXY());
        histos1.fill(HIST("QAbefore/DCAz_pi"), track.dcaZ());
        // Pseudo-rapidity QA (before cuts)
        histos1.fill(HIST("QAbefore/pi_Eta"), tracketaPi);
        // pT QA (before cuts)
        histos1.fill(HIST("QAbefore/pT_pi"), trackptPi);
      }

      // apply the track cut
      if (!trackCutpp(track) || !selectionPIDpp(track))
        continue;

      histos1.fill(HIST("QAafter/hGoodTracksV0s"), 0.5);

      if (!IsMix) {
        // DCA QA (before cuts)
        histos1.fill(HIST("QAAfter/DCAxy_pi"), track.dcaXY());
        histos1.fill(HIST("QAAfter/DCAz_pi"), track.dcaZ());
        // Pseudo-rapidity QA (before cuts)
        histos1.fill(HIST("QAAfter/pi_Eta"), tracketaPi);
        // pT QA (before cuts)
        histos1.fill(HIST("QAAfter/pT_pi"), trackptPi);
        // TPC PID (before cuts)
        histos1.fill(HIST("QAAfter/tpcNsigmaPionQA"), trackptPi,
                     track.tpcNSigmaPi());
      }

      for (const auto& v0 : dV0s) {

        // Full index policy is needed to consider all possible combinations
        if (v0.indices()[0] == trackId || v0.indices()[1] == trackId)
          continue; // To avoid combining secondary and primary pions
                    //// Initialize variables
        // trk: Pion, v0: K0s
        // apply the track cut
        if (!v0cut(v0))
          continue;
        histos1.fill(HIST("QAafter/hGoodTracksV0s"), 1.5);

        lDecayDaughter.SetXYZM(track.px(), track.py(), track.pz(), massPi);
        lDecayV0.SetXYZM(v0.px(), v0.py(), v0.pz(), massK0);
        lResonance = lDecayDaughter + lDecayV0;
        // Counting how many resonances passed
        histos1.fill(HIST("QAafter/hGoodTracksV0s"), 2.5);

        // Checking whether the mid-rapidity condition is met
        if (std::abs(lResonance.Rapidity()) > 0.5)
          continue;
        if constexpr (!IsMix) {
          histos1.fill(HIST("chargedkstarinvmassUlikeSign"), lResonance.M());
          // Reconstructed K*(892)pm 3d mass, pt, multiplicity histogram
          histos1.fill(HIST("chargekstarMassPtMultPtUnlikeSign"),
                       lResonance.M(), lResonance.Pt(), multiplicity);

        } else {
          histos1.fill(HIST("chargedkstarinvmassMixedEvent"), lResonance.M());
          // Reconstructed K*(892)pm 3d mass, pt, multiplicity histogram
          histos1.fill(HIST("chargekstarMassPtMultPtMixedEvent"),
                       lResonance.M(), lResonance.Pt(), multiplicity);
        }
      }
    }
  }

  template <typename T>
  bool selectionPIDpp(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaPi()) < nSigmaCutTPC) {
      tpcPIDPassed = true;
    }
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPi()) < nSigmaCutTOF) {
        tofPIDPassed = true;
      }
      if ((nsigmaCutCombinedPion > 0) &&
          (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi() +
             candidate.tofNSigmaPi() * candidate.tofNSigmaPi() <
           nsigmaCutCombinedPion * nsigmaCutCombinedPion)) {
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

  template <typename TrackType>
  bool trackCutpp(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.eta()) > confdaugheta)
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
  bool v0cut(const V0Type v0)
  {
    // V0 track cuts
    if (std::abs(v0.eta()) > confdaugheta)
      return false;
    if (v0.v0CosPA() < cV0MinCosPA)
      return false;
    if (v0.daughDCA() > cV0MaxDaughDCA)
      return false;

    // apply the competing V0 rejection cut (excluding Lambda0 candidates,
    // massLambdaPDG = 1115.683 MeV/c2)

    if (std::abs(v0.mLambda() - massLambda0) < cV0MassWindow)
      return false;
    if (std::abs(v0.mAntiLambda() - massAntiLambda0) < cV0MassWindow)
      return false;

    return true;
  }

  /*
  SameKindPair<EventCandidates, TrackCandidates,
       BinningTypeVertexContributor>
    pair{binningOnPositions, cfgNoMixedEvents, -1, &cache};
  */

  void processSEnew(aod::ResoCollision const& collision,
                    aod::ResoTracks const& resotracks,
                    aod::ResoV0s const& resov0s)
  {
    // Fill the event counter
    histos1.fill(HIST("hVertexZ"), collision.posZ());
    fillHistograms<false, false>(collision, resotracks,
                                 resov0s); // Fill histograms, no MC, no mixing
  }
  PROCESS_SWITCH(chargedkstaranalysis, processSEnew, "Process Same event new",
                 true);

  using BinningTypeVtxZT0M =
    ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMEnew(aod::ResoCollisions const& collisions,
                    aod::ResoTracks const& resotracks,
                    aod::ResoV0s const& resov0s)
  {
    auto tracksV0sTuple = std::make_tuple(resotracks, resov0s);
    auto v0stuple = std::make_tuple(resov0s);
    BinningTypeVtxZT0M colBinning{{cfgvtxbins, cfgmultbins}, true};
    Pair<aod::ResoCollisions, aod::ResoTracks, aod::ResoV0s, BinningTypeVtxZT0M>
      pairs{colBinning, nEvtMixing, -1, collisions,
            tracksV0sTuple, &cache}; // -1 is the number of the bin to skip
    for (const auto& [c1, restrk1, c2, resov0s2] : pairs) {
      fillHistograms<false, true>(c1, restrk1, resov0s2);
    }
  }
  PROCESS_SWITCH(chargedkstaranalysis, processMEnew, "Process Mixed events new",
                 true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<chargedkstaranalysis>(cfgc)};
}
