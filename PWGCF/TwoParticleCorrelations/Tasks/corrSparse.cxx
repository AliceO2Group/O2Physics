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

/// \file corrSparse.cxx
/// \brief Provides a sparse with usefull two particle correlation info
/// \author Thor Jensen (djt288@alumni.ku.dk)

#include <vector>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "CommonConstants/MathConstants.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct CorrSparse {
  O2_DEFINE_CONFIGURABLE(cfgZVtxCut, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMin, float, 0.2f, "minimum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMax, float, 10.0f, "maximum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgEtaCut, float, 0.8f, "Eta cut")
  O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -2, 2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity / centrality axis for histograms"};
  ConfigurableAxis vtxMix{"vtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis multMix{"multMix", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity / centrality axis for mixed event histograms"};

  // make the filters and cuts.

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgZVtxCut && (aod::evsel::sel8) == true;

  Filter trackFilter = (nabs(aod::track::eta) < cfgEtaCut) && (aod::track::pt > cfgPtCutMin) && (aod::track::pt < cfgPtCutMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  // define the filtered collisions and tracks

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel, o2::aod::BarrelMults>>; // aod::CentFT0Cs
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>>;

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    LOGF(info, "Starting init");
    // Make histograms to check the distributions after cuts
    registry.add("deltaEta_deltaPhi_same", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}}); // check to see the delta eta and delta phi distribution
    registry.add("deltaEta_deltaPhi_mixed", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}});
    registry.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
    registry.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
    registry.add("pT", "pT", {HistType::kTH1D, {axisPtTrigger}});
    registry.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
    registry.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});

    registry.add("Sparse_mixed", "", {HistType::kTHnSparseF, {{axisMultiplicity, axisVertex, axisPtTrigger, axisPtAssoc, axisDeltaPhi, axisDeltaEta}}}); // Make the output sparse
    registry.add("Sparse_same", "", {HistType::kTHnSparseF, {{axisMultiplicity, axisVertex, axisPtTrigger, axisPtAssoc, axisDeltaPhi, axisDeltaEta}}});
    registry.add("Trig_Hist", "", {HistType::kTH3F, {{axisMultiplicity, axisVertex, axisPtTrigger}}});

    registry.add("eventcount", "bin", {HistType::kTH1F, {{3, 0, 3, "bin"}}}); // histogram to see how many events are in the same and mixed event
  }
  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };

  // fill multiple histograms
  template <typename TCollision, typename TTracks>
  void fillYield(TCollision collision, TTracks tracks) // function to fill the yield and etaphi histograms.
  {
    registry.fill(HIST("Nch"), tracks.size());
    registry.fill(HIST("zVtx"), collision.posZ());

    for (auto const& track1 : tracks) {
      registry.fill(HIST("Phi"), track1.phi());
      registry.fill(HIST("Eta"), track1.eta());
      registry.fill(HIST("pT"), track1.pt());
    }
  }

  template <typename TTracks>
  void fillCorrelations(TTracks tracks1, TTracks tracks2, float posZ, int system, float Nch) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (system == SameEvent) {
        registry.fill(HIST("Trig_Hist"), Nch, posZ, track1.pt());
      }

      for (auto const& track2 : tracks2) {

        if (track1.pt() <= track2.pt())
          continue; // skip if the trigger pt is less than the associate p

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        // fill the right sparse and histograms
        if (system == SameEvent) {
          registry.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta);
          registry.fill(HIST("Sparse_same"), Nch, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta);
        } else if (system == MixedEvent) {
          registry.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta);
          registry.fill(HIST("Sparse_mixed"), Nch, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta);
        }
      }
    }
  }

  void processSame(AodCollisions::iterator const& collision, AodTracks const& tracks)
  {

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin
    fillYield(collision, tracks);
    fillCorrelations(tracks, tracks, collision.posZ(), SameEvent, tracks.size()); // fill the SE histogram and Sparse
  }
  PROCESS_SWITCH(CorrSparse, processSame, "Process same event", true);

  // event mixing

  SliceCache cache;
  using MixedBinning = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;

  // the process for filling the mixed events
  void processMixed(AodCollisions const& collisions, AodTracks const& tracks)
  {
    MixedBinning binningOnVtxAndMult{{vtxMix, multMix}, true}; // true is for 'ignore overflows' (true by default)
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<AodCollisions, AodTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto const& [collision1, tracks1, collision2, tracks2] : pairs) {
      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      fillCorrelations(tracks1, tracks2, collision1.posZ(), MixedEvent, tracks1.size());
    }
  }
  PROCESS_SWITCH(CorrSparse, processMixed, "Process mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CorrSparse>(cfgc),
  };
}
