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
///
/// \brief This task serves to do hadron-(strange hadron) correlation studies.
///  The yield will be calculated using the two-particle correlation method.
///  Trigger particle : Hadrons
///  Associated Particles : V0s or Cascades
///  this task requires the hStrangeCorrelationFilter to have been run before.
///
/// \author Kai Cui
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFHStrangeCorrelationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/Multiplicity.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksComplete = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;

struct correlateSameEvents {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  /// Function to aid in calculating delta-phi
  /// \param phi1 first phi value
  /// \param phi2 second phi value
  Double_t ComputeDeltaPhi(Double_t phi1, Double_t phi2)
  {
    // To be completely sure, use inner products
    Double_t x1, y1, x2, y2;
    x1 = TMath::Cos(phi1);
    y1 = TMath::Sin(phi1);
    x2 = TMath::Cos(phi2);
    y2 = TMath::Sin(phi2);
    Double_t lInnerProd = x1 * x2 + y1 * y2;
    Double_t lVectorProd = x1 * y2 - x2 * y1;

    Double_t lReturnVal = 0;
    if (lVectorProd > 1e-8) {
      lReturnVal = TMath::ACos(lInnerProd);
    }
    if (lVectorProd < -1e-8) {
      lReturnVal = -TMath::ACos(lInnerProd);
    }

    if (lReturnVal < -TMath::Pi() / 2.) {
      lReturnVal += 2. * TMath::Pi();
    }

    return lReturnVal;
  }

  void init(InitContext const&)
  {
    // define usual axes to be used
    const AxisSpec axisPhi{40, -0.5 * M_PI, 1.5 * M_PI, "#phi"};
    const AxisSpec axisEta{400, -1, 1, "#eta"};
    const AxisSpec axisDeltaEta{400, -2, 2, "#Delta#eta"};
    const AxisSpec axisPt{100, 0, 10, "#it{p}_{T} (GeV/c)"};
    const AxisSpec axisK0ShortMass{400, 0.400f, 0.600f, "Inv. Mass (GeV/c^{2})"};
    const AxisSpec axisLambdaMass{400, 1.01f, 1.21f, "Inv. Mass (GeV/c^{2})"};

    // correlation histograms in phi alone (warning: not mixed-event-corrected)
    histos.add("correlationPhiHadronK0short", "correlationPhiHadronK0short", kTH1F, {axisPhi});
    histos.add("correlationPhiHadronLambda", "correlationPhiHadronLambda", kTH1F, {axisPhi});
    histos.add("correlationPhiHadronAntiLambda", "correlationPhiHadronAntiLambda", kTH1F, {axisPhi});
    histos.add("correlationPhiHadronXiMinus", "correlationPhiHadronXiMinus", kTH1F, {axisPhi});
    histos.add("correlationPhiHadronXiPlus", "correlationPhiHadronXiPlus", kTH1F, {axisPhi});
    histos.add("correlationPhiHadronOmegaMinus", "correlationPhiHadronOmegaMinus", kTH1F, {axisPhi});
    histos.add("correlationPhiHadronOmegaPlus", "correlationPhiHadronOmegaPlus", kTH1F, {axisPhi});

    // full correlation functions
    histos.add("correlationFullHadronK0Short", "correlationFullHadronK0Short", kTH2F, {axisPhi, axisDeltaEta});
    histos.add("correlationFullHadronLambda", "correlationFullHadronLambda", kTH2F, {axisPhi, axisDeltaEta});
    histos.add("correlationFullHadronAntiLambda", "correlationFullHadronAntiLambda", kTH2F, {axisPhi, axisDeltaEta});
    histos.add("correlationFullHadronXiMinus", "correlationFullHadronXiMinus", kTH2F, {axisPhi, axisDeltaEta});
    histos.add("correlationFullHadronXiPlus", "correlationFullHadronXiPlus", kTH2F, {axisPhi, axisDeltaEta});
    histos.add("correlationFullHadronOmegaMinus", "correlationFullHadronOmegaMinus", kTH2F, {axisPhi, axisDeltaEta});
    histos.add("correlationFullHadronOmegaPlus", "correlationFullHadronOmegaPlus", kTH2F, {axisPhi, axisDeltaEta});

    // Some QA plots
    histos.add("hMassK0Short", "hMassK0Short", kTH1F, {axisK0ShortMass});
    histos.add("hMassLambda", "hMassLambda", kTH1F, {axisLambdaMass});
    histos.add("hMassAntiLambda", "hMassAntiLambda", kTH1F, {axisLambdaMass});
    histos.add("hTrackEta", "hTrackEta", kTH1F, {axisEta});
    histos.add("hV0Eta", "hV0Eta", kTH1F, {axisEta});
    histos.add("hPtK0Short", "hPtK0Short", kTH1F, {axisPt});
    histos.add("hPtLambda", "hPtLambda", kTH1F, {axisPt});
    histos.add("hPtAntiLambda", "hPtAntiLambda", kTH1F, {axisPt});
  }

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator const& collision,
               aod::TriggerTracks const& triggerTracks, aod::AssocV0s const& associatedV0s, aod::AssocCascades const& associatedCascades, aod::V0Datas const&, TracksComplete const&, soa::Filtered<aod::V0Datas> const& V0s, aod::CascDataExt const& Cascades, aod::V0sLinked const&)
  {
    // ________________________________________________
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    if (TMath::Abs(collision.posZ()) > 10.0) {
      return;
    }
    // ________________________________________________
    // Do basic QA
    for (auto const& v0 : associatedV0s) {
      auto v0Data = v0.v0Data();
      histos.fill(HIST("hV0Eta"), v0Data.eta());

      if (v0.compatibleK0Short()) {
        // K0Short compatible
        histos.fill(HIST("hMassK0Short"), v0Data.mK0Short());
        histos.fill(HIST("hPtK0Short"), v0Data.pt());
      }

      if (v0.compatibleLambda()) {
        // Lambda compatible
        histos.fill(HIST("hMassLambda"), v0Data.mLambda());
        histos.fill(HIST("hPtLambda"), v0Data.pt());
      }
      if (v0.compatibleAntiLambda()) {
        // AntiLambda compatible
        histos.fill(HIST("hMassAntiLambda"), v0Data.mAntiLambda());
        histos.fill(HIST("hPtAntiLambda"), v0Data.pt());
      }
    }
    for (auto const& triggerTrack : triggerTracks) {
      auto track = triggerTrack.track_as<aod::Tracks>();
      histos.fill(HIST("hTrackEta"), track.eta());
    }

    // ________________________________________________
    // Do hadron - V0 correlations
    for (auto const& [triggerTrackRef, assocTrackRef] : combinations(o2::soa::CombinationsFullIndexPolicy(triggerTracks, associatedV0s))) {
      // De-reference
      auto triggerTrack = triggerTrackRef.track_as<aod::Tracks>();
      auto assocTrack = assocTrackRef.v0Data();
      // Correlate
      if (assocTrackRef.compatibleK0Short()) {
        histos.fill(HIST("correlationPhiHadronK0Short"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()));
        histos.fill(HIST("correlationFullHadronK0Short"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
      }
      if (assocTrackRef.compatibleLambda()) {
        histos.fill(HIST("correlationPhiHadronLambda"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()));
        histos.fill(HIST("correlationFullHadronLambda"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
      }
      if (assocTrackRef.compatibleAntiLambda()) {
        histos.fill(HIST("correlationPhiHadronAntiLambda"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()));
        histos.fill(HIST("correlationFullHadronAntiLambda"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
      }
    }

    // ________________________________________________
    // Do hadron - cascade correlations
    for (auto const& [triggerTrackRef, assocTrackRef] : combinations(o2::soa::CombinationsFullIndexPolicy(triggerTracks, associatedCascades))) {
      // De-reference
      auto triggerTrack = triggerTrackRef.track_as<aod::Tracks>();
      auto assocTrack = assocTrackRef.cascData();
      // Correlate
      if (assocTrackRef.compatibleXiMinus()) {
        histos.fill(HIST("correlationPhiHadronXiMinus"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()));
        histos.fill(HIST("correlationFullHadronXiMinus"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
      }
      if (assocTrackRef.compatibleXiPlus()) {
        histos.fill(HIST("correlationPhiHadronXiPlus"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()));
        histos.fill(HIST("correlationFullHadronXiPlus"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
      }
      if (assocTrackRef.compatibleOmegaMinus()) {
        histos.fill(HIST("correlationPhiHadronOmegaMinus"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()));
        histos.fill(HIST("correlationFullHadronOmegaMinus"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
      }
      if (assocTrackRef.compatibleOmegaPlus()) {
        histos.fill(HIST("correlationPhiHadronOmegaPlus"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()));
        histos.fill(HIST("correlationFullHadronOmegaPlus"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
      }
    }
  }
};

struct correlateMixedEvents {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  // using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0A>;
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
  BinningType colBinning{{ConfVtxBins, ConfMultBins}, true};
  Pair<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>, aod::TriggerTracks, aod::AssocV0s, BinningType> pairV0s{colBinning, 5, -1}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored
  Pair<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>, aod::TriggerTracks, aod::AssocCascades, BinningType> pairCascades{colBinning, 5, -1}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  /// Function to aid in calculating delta-phi
  /// \param phi1 first phi value
  /// \param phi2 second phi value
  Double_t ComputeDeltaPhi(Double_t phi1, Double_t phi2)
  {
    // To be completely sure, use inner products
    Double_t x1, y1, x2, y2;
    x1 = TMath::Cos(phi1);
    y1 = TMath::Sin(phi1);
    x2 = TMath::Cos(phi2);
    y2 = TMath::Sin(phi2);
    Double_t lInnerProd = x1 * x2 + y1 * y2;
    Double_t lVectorProd = x1 * y2 - x2 * y1;

    Double_t lReturnVal = 0;
    if (lVectorProd > 1e-8) {
      lReturnVal = TMath::ACos(lInnerProd);
    }
    if (lVectorProd < -1e-8) {
      lReturnVal = -TMath::ACos(lInnerProd);
    }

    if (lReturnVal < -TMath::Pi() / 2.) {
      lReturnVal += 2. * TMath::Pi();
    }

    return lReturnVal;
  }

  void init(InitContext const&)
  {
    // define usual axes to be used
    const AxisSpec axisPhi{40, -0.5 * M_PI, 1.5 * M_PI, "#phi"};
    const AxisSpec axisEta{400, -1, 1, "#eta"};
    const AxisSpec axisDeltaEta{400, -2, 2, "#Delta#eta"};

    // correlation histograms in phi alone (warning: not mixed-event-corrected)
    histos.add("deltaPhiDistribution", "deltaPhiDistribution", kTH1F, {axisPhi});
    histos.add("deltaEtaDistribution", "deltaEtaDistribution", kTH1F, {axisEta});

    // full correlation functions
    histos.add("mixedFullHadronK0Short", "mixedFullHadronK0Short", kTH2F, {axisPhi, axisEta});
    histos.add("mixedFullHadronLambda", "mixedFullHadronLambda", kTH2F, {axisPhi, axisEta});
    histos.add("mixedFullHadronAntiLambda", "mixedFullHadronAntiLambda", kTH2F, {axisPhi, axisEta});
    histos.add("mixedFullHadronXiMinus", "mixedFullHadronXiMinus", kTH2F, {axisPhi, axisEta});
    histos.add("mixedFullHadronXiPlus", "mixedFullHadronXiPlus", kTH2F, {axisPhi, axisEta});
    histos.add("mixedFullHadronOmegaMinus", "mixedFullHadronOmegaMinus", kTH2F, {axisPhi, axisEta});
    histos.add("mixedFullHadronOmegaPlus", "mixedFullHadronOmegaPlus", kTH2F, {axisPhi, axisEta});
  }

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults> const& collisions, aod::TriggerTracks const& tracks, aod::AssocV0s const& assocv0s, aod::AssocCascades const& associatedCascades, aod::V0Datas const&, aod::CascDatas const&)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d V0s %d", collisions.size(), tracks.size(), assocv0s.size());
    // ________________________________________________
    // Mixed event loop for associated = V0s
    for (auto& [c1, tracks1, c2, tracks2] : pairV0s) {
      LOGF(info, "Mixed event collisions: (%d, %d)", c1.globalIndex(), c2.globalIndex());
      for (auto& [t1, t2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", t1.collisionId(), t2.collisionId(), c1.index(), c2.index(), t1.collision().index(), t2.collision().index());
        // De-reference
        auto triggerTrack = t1.track_as<aod::Tracks>();
        auto assocTrack = t2.v0Data();
        // do correlations
        histos.fill(HIST("deltaPhiDistribution"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()));
        histos.fill(HIST("deltaEtaDistribution"), triggerTrack.eta() - assocTrack.eta());
        // Correlate
        if (t2.compatibleK0Short()) {
          histos.fill(HIST("mixedHadronK0Short"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
        }
        if (t2.compatibleLambda()) {
          histos.fill(HIST("mixedHadronLambda"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
        }
        if (t2.compatibleAntiLambda()) {
          histos.fill(HIST("mixedHadronAntiLambda"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
        }
      }
    }
    // ________________________________________________
    // Mixed event loop for associated = Cascades
    for (auto& [c1, tracks1, c2, tracks2] : pairCascades) {
      LOGF(info, "Mixed event collisions: (%d, %d)", c1.globalIndex(), c2.globalIndex());
      for (auto& [t1, t2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", t1.collisionId(), t2.collisionId(), c1.index(), c2.index(), t1.collision().index(), t2.collision().index());
        // De-reference
        auto triggerTrack = t1.track_as<aod::Tracks>();
        auto assocTrack = t2.cascData();
        // Correlate
        if (t2.compatibleXiMinus()) {
          histos.fill(HIST("mixedHadronXiMinus"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
        }
        if (t2.compatibleXiPlus()) {
          histos.fill(HIST("mixedHadronXiPlus"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
        }
        if (t2.compatibleOmegaMinus()) {
          histos.fill(HIST("mixedHadronOmegaMinus"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
        }
        if (t2.compatibleOmegaPlus()) {
          histos.fill(HIST("mixedHadronOmegaPlus"), ComputeDeltaPhi(triggerTrack.phi(), assocTrack.phi()), triggerTrack.eta() - assocTrack.eta());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<correlateSameEvents>(cfgc),
    adaptAnalysisTask<correlateMixedEvents>(cfgc)};
}
