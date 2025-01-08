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
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct CorrSparse {
  Configurable<float> cfgZVtxCut = {"cfgZVtxCut", 10.0, "Vertex z cut. Default 10 cm"};
  Configurable<float> cfgPtCutMin = {"cfgPtCutMin", 0.2, "Minimum accepted track pT. Default 0.2 GeV"};
  Configurable<float> cfgPtCutMax = {"cfgPtCutMax", 5.0, "Maximum accepted track pT. Default 5.0 GeV"};
  Configurable<float> cfgEtaCut = {"cfgEtaCut", 0.8, "Eta cut. Default 0.8"};
  ConfigurableAxis axisVertex{"axisVertex", {7, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -2, 2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity / centrality axis for histograms"};
  HistogramRegistry registry{"registry"};
  int logcolls = 0;
  int logcollpairs = 0;

  void init(InitContext&)
  {
    LOGF(info, "Starting init");
    registry.add("Yield", "pT vs eta vs Nch", {HistType::kTH3F, {{40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}, {100, 0, 100, "Nch"}}}); // check to see total number of tracks
    registry.add("etaphi_Trigger", "eta vs phi vs Nch", {HistType::kTH3F, {{100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}, {100, 0, 100, "Nch"}}});

    // Make histograms to check the distributions after cuts
    registry.add("deltaEta_deltaPhi_same", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}}); // check to see the delta eta and delta phi distribution
    registry.add("deltaEta_deltaPhi_mixed", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}});
    registry.add("Phi", "Phi", {HistType::kTH1D, {{72, 0, TwoPI, "Phi"}}});
    registry.add("Eta", "Eta", {HistType::kTH1D, {{40, -2, 2, "Eta"}}});
    registry.add("pT", "pT", {HistType::kTH1D, {axisPtTrigger}});
    registry.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});

    registry.add("Sparse_mixed", "", {HistType::kTHnSparseF, {{axisVertex, axisPtTrigger, axisPtAssoc, axisMultiplicity, axisDeltaPhi, axisDeltaEta}}}); // Make the output sparse
    registry.add("Sparse_same", "", {HistType::kTHnSparseF, {{axisVertex, axisPtTrigger, axisPtAssoc, axisMultiplicity, axisDeltaPhi, axisDeltaEta}}});

    const int maxMixBin = axisMultiplicity->size() * axisVertex->size();
    registry.add("eventcount", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}}); // histogram to see how many events are in the same and mixed event
  }

  // fill multiple histograms
  template <typename TCollision, typename TTracks>
  void fillYield(TCollision /*collision*/, float /*centrality*/, TTracks tracks) // function to fill the yield and etaphi histograms.
  {
    registry.fill(HIST("Nch"), tracks.size());
    for (auto const& track1 : tracks) {
      registry.fill(HIST("Yield"), track1.pt(), track1.eta(), track1.size());
      registry.fill(HIST("etaphi_Trigger"), track1.eta(), track1.phi(), track1.size());
      registry.fill(HIST("Phi"), track1.phi());
      registry.fill(HIST("Eta"), track1.eta());
      registry.fill(HIST("pT"), track1.pt());
    }
  }

  template <typename TCollision>
  bool fillCollision(TCollision collision, float /*centrality*/)
  {

    if (!collision.sel8()) {
      return false;
    }

    return true;
  }

  template <typename TTracks>
  void fillCorrelations(TTracks tracks1, TTracks tracks2, float posZ, int system, float Nch) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      for (auto const& track2 : tracks2) {
        if (track1 == track2) {
          continue;
        }

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        // fill the right sparse and histograms
        if (system == 1) {
          registry.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta);
          registry.fill(HIST("Sparse_same"), posZ, track1.pt(), track2.pt(), Nch, deltaPhi, deltaEta);
        } else if (system == 2) {
          registry.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta);
          registry.fill(HIST("Sparse_mixed"), posZ, track1.pt(), track2.pt(), Nch, deltaPhi, deltaEta);
        }
      }
    }
  }

  // make the filters and cuts.

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgZVtxCut;

  Filter trackFilter = (nabs(aod::track::eta) < cfgEtaCut) && (aod::track::pt > cfgPtCutMin) && (aod::track::pt < cfgPtCutMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  // define the filtered collisions and tracks

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;

  void processSame(AodCollisions::iterator const& collision, AodTracks const& tracks)
  {
    const auto centrality = collision.centFT0C();

    registry.fill(HIST("eventcount"), -2); // because its same event i put it in the -2 bin
    fillYield(collision, centrality, tracks);
    fillCorrelations(tracks, tracks, collision.posZ(), 1, tracks.size()); // fill the SE histogram and Sparse
  }
  PROCESS_SWITCH(CorrSparse, processSame, "Process same event", true);

  // i do the event mixing (i have not changed this from the tutorial i got).
  std::vector<double> vtxBinsEdges{VARIABLE_WIDTH, -7.0f, -5.0f, -3.0f, -1.0f, 1.0f, 3.0f, 5.0f, 7.0f};
  std::vector<double> multBinsEdges{VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0, 100.1f};
  SliceCache cache;

  ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>
    bindingOnVtxAndMult{{vtxBinsEdges, multBinsEdges}, true}; // true is for 'ignore overflows' (true by default)
  SameKindPair<AodCollisions,
               AodTracks,
               ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>>
    pair{bindingOnVtxAndMult, 5, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  // the process for filling the mixed events
  void processMixed(AodCollisions const& /*collisions*/, AodTracks const& /*tracks*/)
  {
    for (auto const& [collision1, tracks1, collision2, tracks2] : pair) {

      if (fillCollision(collision1, collision1.centFT0C()) == false) {
        continue;
      }

      registry.fill(HIST("eventcount"), 1); // fill the mixed event in the 1 bin

      fillCorrelations(tracks1, tracks2, collision1.posZ(), 2, tracks1.size());
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
