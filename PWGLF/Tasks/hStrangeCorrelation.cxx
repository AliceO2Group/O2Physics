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
/// \author Kai Cui (kaicui@mails.ccnu.edu.cn)
/// \author Lucia Anna Tarasovicova (lucia.anna.husova@cern.ch)
/// \author David Dobrigkeit Chinellato (david.dobrigkeit.chinellato@cern.ch)
/// \author Zhongbao Yin (Zhong-Bao.Yin@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFHStrangeCorrelationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/Multiplicity.h"
#include "Framework/StaticFor.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksComplete = soa::Join<aod::Tracks, aod::TracksExtra>;

struct correlateStrangeness {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsPhi{"nBinsPhi", 72, "Number of phi bins"};
  Configurable<int> nBinsEta{"nBinsEta", 80, "Number of eta bins"};
  Configurable<int> nBinsDeltaEta{"nBinsDeltaEta", 160, "Number of delta-eta bins"};
  Configurable<int> nBinsDeltaPhi{"nBinsDeltaPhi", 144, "Number of delta-phi bins"};
  Configurable<int> nBinsMass{"nBinsMass", 200, "Number of mass bins"};
  Configurable<int> zVertexCut{"zVertexCut", 10, "Cut on PV position"};

  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;

  // collision slicing for mixed events
  Preslice<aod::TriggerTracks> collisionSliceTracks = aod::triggerTracks::collisionId;
  Preslice<aod::AssocV0s> collisionSliceV0s = aod::assocV0s::collisionId;
  Preslice<aod::AssocCascades> collisionSliceCascades = aod::assocCascades::collisionId;

  static constexpr std::string_view v0names[] = {"K0Short", "Lambda", "AntiLambda"};
  static constexpr std::string_view cascadenames[] = {"XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};

  /// Function to aid in calculating delta-phi
  /// \param phi1 first phi value
  /// \param phi2 second phi value
  Double_t ComputeDeltaPhi(Double_t phi1, Double_t phi2)
  {
    Double_t deltaPhi = phi1 - phi2;
    if (deltaPhi < -TMath::Pi() / 2.) {
      deltaPhi += 2. * TMath::Pi();
    }
    if (deltaPhi > 3 * TMath::Pi() / 2.) {
      deltaPhi -= 2. * TMath::Pi();
    }
    return deltaPhi;
  }

  void fillCorrelationsV0(aod::TriggerTracks const& triggers, aod::AssocV0s const& assocs, bool mixing)
  {
    for (auto& triggerTrack : triggers) {
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!mixing)
        histos.fill(HIST("sameEvent/TriggerParticlesV0"), trigg.pt());
      for (auto& assocCandidate : assocs) {
        auto assoc = assocCandidate.v0Data();

        //---] removing autocorrelations [---
        auto postrack = assoc.posTrack_as<TracksComplete>();
        auto negtrack = assoc.negTrack_as<TracksComplete>();
        if (trigg.globalIndex() == postrack.globalIndex())
          continue;
        if (trigg.globalIndex() == negtrack.globalIndex())
          continue;
        // TODO: add histogram checking how many pairs are rejected (should be small!)

        float deltaphi = ComputeDeltaPhi(trigg.phi(), assoc.phi());
        float deltaeta = trigg.eta() - assoc.eta();
        float ptassoc = assoc.pt();
        static_for<0, 2>([&](auto i) {
          constexpr int index = i.value;
          if (assocCandidate.compatible(index) && !mixing && assocCandidate.inMassRegionCheck(index, 1))
            histos.fill(HIST("sameEvent/LeftBg") + HIST(v0names[index]), deltaphi, deltaeta, ptassoc);
          if (assocCandidate.compatible(index) && !mixing && assocCandidate.inMassRegionCheck(index, 2))
            histos.fill(HIST("sameEvent/Signal") + HIST(v0names[index]), deltaphi, deltaeta, ptassoc);
          if (assocCandidate.compatible(index) && !mixing && assocCandidate.inMassRegionCheck(index, 3))
            histos.fill(HIST("sameEvent/RightBg") + HIST(v0names[index]), deltaphi, deltaeta, ptassoc);
          if (assocCandidate.compatible(index) && mixing && assocCandidate.inMassRegionCheck(index, 1))
            histos.fill(HIST("mixedEvent/LeftBg") + HIST(v0names[index]), deltaphi, deltaeta, ptassoc);
          if (assocCandidate.compatible(index) && mixing && assocCandidate.inMassRegionCheck(index, 2))
            histos.fill(HIST("mixedEvent/Signal") + HIST(v0names[index]), deltaphi, deltaeta, ptassoc);
          if (assocCandidate.compatible(index) && mixing && assocCandidate.inMassRegionCheck(index, 3))
            histos.fill(HIST("mixedEvent/RightBg") + HIST(v0names[index]), deltaphi, deltaeta, ptassoc);
        });
      }
    }
  }

  void fillCorrelationsCascade(aod::TriggerTracks const& triggers, aod::AssocCascades const& assocs, bool mixing)
  {
    for (auto& triggerTrack : triggers) {
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!mixing)
        histos.fill(HIST("sameEvent/TriggerParticlesCascade"), trigg.pt());
      for (auto& assocCandidate : assocs) {
        auto assoc = assocCandidate.cascData();

        //---] removing autocorrelations [---
        auto v0index = assoc.v0_as<o2::aod::V0sLinked>();
        if (!(v0index.has_v0Data()))
          continue;                      // this should not happen - included for safety
        auto assocV0 = v0index.v0Data(); // de-reference index to correct v0data in case it exists
        auto postrack = assocV0.posTrack_as<TracksComplete>();
        auto negtrack = assocV0.negTrack_as<TracksComplete>();
        auto bachtrack = assoc.bachelor_as<TracksComplete>();
        if (trigg.globalIndex() == postrack.globalIndex())
          continue;
        if (trigg.globalIndex() == negtrack.globalIndex())
          continue;
        if (trigg.globalIndex() == bachtrack.globalIndex())
          continue;
        // TODO: add histogram checking how many pairs are rejected (should be small!)

        float deltaphi = ComputeDeltaPhi(trigg.phi(), assoc.phi());
        float deltaeta = trigg.eta() - assoc.eta();
        float ptassoc = assoc.pt();
        static_for<0, 3>([&](auto i) {
          constexpr int index = i.value;
          if (assocCandidate.compatible(index) && !mixing && assocCandidate.inMassRegionCheck(index, 1))
            histos.fill(HIST("sameEvent/LeftBg") + HIST(cascadenames[index]), deltaphi, deltaeta, ptassoc);
          if (assocCandidate.compatible(index) && !mixing && assocCandidate.inMassRegionCheck(index, 2))
            histos.fill(HIST("sameEvent/Signal") + HIST(cascadenames[index]), deltaphi, deltaeta, ptassoc);
          if (assocCandidate.compatible(index) && !mixing && assocCandidate.inMassRegionCheck(index, 3))
            histos.fill(HIST("sameEvent/RightBg") + HIST(cascadenames[index]), deltaphi, deltaeta, ptassoc);
          if (assocCandidate.compatible(index) && mixing && assocCandidate.inMassRegionCheck(index, 1))
            histos.fill(HIST("mixedEvent/LeftBg") + HIST(cascadenames[index]), deltaphi, deltaeta, ptassoc);
          if (assocCandidate.compatible(index) && mixing && assocCandidate.inMassRegionCheck(index, 2))
            histos.fill(HIST("mixedEvent/Signal") + HIST(cascadenames[index]), deltaphi, deltaeta, ptassoc);
          if (assocCandidate.compatible(index) && mixing && assocCandidate.inMassRegionCheck(index, 3))
            histos.fill(HIST("mixedEvent/RightBg") + HIST(cascadenames[index]), deltaphi, deltaeta, ptassoc);
        });
      }
    }
  }

  void init(InitContext const&)
  {
    // define usual axes to be used
    const AxisSpec axisPhi{nBinsPhi, -0.5 * M_PI, 1.5 * M_PI, "#phi"};
    const AxisSpec axisEta{nBinsEta, -0.8, +0.8, "#eta"};
    const AxisSpec axisDeltaEta{nBinsDeltaEta, -1.6, 1.6, "#Delta#eta"};
    const AxisSpec axisDeltaPhi{nBinsDeltaPhi, -0.5 * M_PI, 1.5 * M_PI, "#Delta#phi"};
    const AxisSpec axisPtFine{100, 0, 10, "#it{p}_{T} (GeV/c)"};
    const AxisSpec axisPt{10, 0, 5, "#it{p}_{T} (GeV/c)"};
    const AxisSpec axisK0ShortMass{nBinsMass, 0.400f, 0.600f, "Inv. Mass (GeV/c^{2})"};
    const AxisSpec axisLambdaMass{nBinsMass, 1.01f, 1.21f, "Inv. Mass (GeV/c^{2})"};
    const AxisSpec axisXiMass{nBinsMass, 1.22f, 1.42f, "Inv. Mass (GeV/c^{2})"};
    const AxisSpec axisOmegaMass{nBinsMass, 1.57f, 1.77f, "Inv. Mass (GeV/c^{2})"};

    // same-event correlation functions
    histos.add("sameEvent/SignalK0Short", "SignalK0Short", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/SignalLambda", "SignalLambda", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/SignalAntiLambda", "SignalAntiLambda", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/SignalXiMinus", "SignalXiMinus", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/SignalXiPlus", "SignalXiPlus", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/SignalOmegaMinus", "SignalOmegaMinus", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/SignalOmegaPlus", "SignalOmegaPlus", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});

    histos.add("sameEvent/LeftBgK0Short", "LeftBgK0Short", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/LeftBgLambda", "LeftBgLambda", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/LeftBgAntiLambda", "LeftBgAntiLambda", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/LeftBgXiMinus", "LeftBgXiMinus", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/LeftBgXiPlus", "LeftBgXiPlus", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/LeftBgOmegaMinus", "LeftBgOmegaMinus", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/LeftBgOmegaPlus", "LeftBgOmegaPlus", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});

    histos.add("sameEvent/RightBgK0Short", "RightBgK0Short", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/RightBgLambda", "RightBgLambda", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/RightBgAntiLambda", "RightBgAntiLambda", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/RightBgXiMinus", "RightBgXiMinus", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/RightBgXiPlus", "RightBgXiPlus", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/RightBgOmegaMinus", "RightBgOmegaMinus", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});
    histos.add("sameEvent/RightBgOmegaPlus", "RightBgOmegaPlus", kTH3F, {axisDeltaPhi, axisDeltaEta, axisPt});

    // mixed-event correlation functions
    histos.addClone("sameEvent/", "mixedEvent/");

    // Some QA plots
    histos.add("h2dMassK0Short", "h2dMassK0Short", kTH2F, {axisPtFine, axisK0ShortMass});
    histos.add("h2dMassLambda", "h2dMassLambda", kTH2F, {axisPtFine, axisLambdaMass});
    histos.add("h2dMassAntiLambda", "h2dMassAntiLambda", kTH2F, {axisPtFine, axisLambdaMass});
    histos.add("h2dMassXiMinus", "h2dMassXiMinus", kTH2F, {axisPtFine, axisXiMass});
    histos.add("h2dMassXiPlus", "h2dMassXiPlus", kTH2F, {axisPtFine, axisXiMass});
    histos.add("h2dMassOmegaMinus", "h2dMassOmegaMinus", kTH2F, {axisPtFine, axisOmegaMass});
    histos.add("h2dMassOmegaPlus", "h2dMassOmegaPlus", kTH2F, {axisPtFine, axisOmegaMass});
    histos.add("hTrackEtaVsPt", "hTrackEtaVsPt", kTH2F, {axisPt, axisEta});
    histos.add("hV0EtaVsPt", "hV0EtaVsPt", kTH2F, {axisPt, axisEta});
    histos.add("hCascEtaVsPt", "hCascEtaVsPt", kTH2F, {axisPt, axisEta});

    histos.add("sameEvent/TriggerParticlesV0", "TriggersV0", kTH1F, {{400, 0, 20}});
    histos.add("sameEvent/TriggerParticlesCascade", "TriggersCascade", kTH1F, {{400, 0, 20}});
  }

  void processSameEvent(soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator const& collision,
                        aod::AssocV0s const& associatedV0s, aod::AssocCascades const& associatedCascades, aod::TriggerTracks const& triggerTracks,
                        aod::V0Datas const&, aod::V0sLinked const&, aod::CascDatas const&, TracksComplete const&)
  {
    // ________________________________________________
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    if (TMath::Abs(collision.posZ()) > zVertexCut) {
      return;
    }
    // ________________________________________________
    // Do basic QA
    for (auto const& v0 : associatedV0s) {
      auto v0Data = v0.v0Data();
      histos.fill(HIST("hV0EtaVsPt"), v0Data.pt(), v0Data.eta());
      static_for<0, 2>([&](auto i) {
        constexpr int index = i.value;
        if (v0.compatible(index))
          histos.fill(HIST("h2dMass") + HIST(v0names[index]), v0Data.pt(), v0Data.m(index));
      });
    }
    for (auto const& casc : associatedCascades) {
      auto cascData = casc.cascData();
      histos.fill(HIST("hCascEtaVsPt"), cascData.pt(), cascData.eta());
      static_for<0, 3>([&](auto i) {
        constexpr int index = i.value;
        if (casc.compatible(index))
          histos.fill(HIST("h2dMass") + HIST(cascadenames[index]), cascData.pt(), cascData.m(index));
      });
    }
    for (auto const& triggerTrack : triggerTracks) {
      auto track = triggerTrack.track_as<TracksComplete>();
      histos.fill(HIST("hTrackEtaVsPt"), track.pt(), track.eta());
    }

    // ________________________________________________
    // Do hadron - V0 correlations
    fillCorrelationsV0(triggerTracks, associatedV0s, false);

    // ________________________________________________
    // Do hadron - cascade correlations
    fillCorrelationsCascade(triggerTracks, associatedCascades, false);
  }
  PROCESS_SWITCH(correlateStrangeness, processSameEvent, "Process same events", true);

  void processMixedEvent(soa::Join<aod::Collisions, aod::EvSels, aod::Mults> const& collisions,
                         aod::AssocV0s const& associatedV0s, aod::AssocCascades const& associatedCascades, aod::TriggerTracks const& triggerTracks,
                         aod::V0Datas const&, aod::V0sLinked const&, aod::CascDatas const&, TracksComplete const&)
  {

    BinningType colBinning{{ConfVtxBins, ConfMultBins}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, collisions, collisions)) {
      // ________________________________________________
      // Perform basic event selection on both collisions
      if (!collision1.sel8() || !collision2.sel8())
        continue;
      if (TMath::Abs(collision1.posZ()) > zVertexCut && TMath::Abs(collision2.posZ()) > zVertexCut)
        continue;

      // ________________________________________________
      // Do slicing
      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());
      auto slicedAssocV0s = associatedV0s.sliceBy(collisionSliceV0s, collision2.globalIndex());
      auto slicedAssocCascades = associatedCascades.sliceBy(collisionSliceCascades, collision2.globalIndex());

      // ________________________________________________
      // Do hadron - V0 correlations
      fillCorrelationsV0(slicedTriggerTracks, slicedAssocV0s, true);

      // ________________________________________________
      // Do hadron - cascade correlations
      fillCorrelationsCascade(slicedTriggerTracks, slicedAssocCascades, true);
    }
  }
  PROCESS_SWITCH(correlateStrangeness, processMixedEvent, "Process mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<correlateStrangeness>(cfgc)};
}