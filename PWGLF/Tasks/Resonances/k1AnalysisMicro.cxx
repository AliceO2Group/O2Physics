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
/// \file k1AnalysisMicro.cxx
/// \brief Reconstruction of track-track decay resonance candidates
/// \author Su-Jeong Ji <su-jeong.ji@cern.ch>, Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///

#include "PWGLF/Core/K1AnalysisMicroCore.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <tuple>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::analysis::k1micro;

struct K1AnalysisMicro {
  // Module-initializer v001 tables; full tracks keep their unversioned schema as a fallback.
  using ResoCollisions = aod::ResoCollisions_001;
  using ResoMCCols = soa::Join<ResoCollisions, aod::ResoMCCollisions_001>;
  using ResoTracks = aod::ResoTracks; // no v001 exists; K1 does not need ResoTrackTracks (trackId unused)
  using ResoMicroTracks = aod::ResoMicroTracks_001;
  using ResoMCTracks = soa::Join<ResoTracks, aod::ResoMCTracks>;
  using ResoMCMicroTracks = soa::Join<ResoMicroTracks, aod::ResoMCMicroTracks_001>;
  using ResoMCParents = aod::ResoMCParents_001;

  SliceCache cache;
  // Registered only to enable the slice cache that SameKindPair (event mixing) needs, as in Xi1820Analysis
  Preslice<ResoTracks> perResoCollisionTrack = aod::resodaughter::resoCollisionId;
  Preslice<ResoMicroTracks> perResoCollisionMicroTrack = aod::resodaughter::resoCollisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Selection shared with the K1 training-table task (plain JSON keys, no group prefix)
  EventCuts eventCuts;
  TrackCuts trackCuts;
  PionPidCuts pionPID;
  KaonPidCuts kaonPID;
  SecondaryCuts secondaryCuts;
  CandidateCuts candidateCuts;
  HistogramOptions histogramOptions;

  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};

  K1AnalysisMicroCore core;

  void init(InitContext&)
  {
    const int sameEventModes = static_cast<int>(doprocessResoTracks) + static_cast<int>(doprocessResoMicroTracks) +
                               static_cast<int>(doprocessMC) + static_cast<int>(doprocessMCMicro);
    const int mixedEventModes = static_cast<int>(doprocessME) + static_cast<int>(doprocessMEMicro);
    if (sameEventModes > 1 || mixedEventModes > 1) {
      LOG(fatal) << "Enable at most one same-event mode and one mixing mode";
    }

    ProcessModes modes;
    modes.microTracks = doprocessResoMicroTracks || doprocessMCMicro || doprocessMEMicro;
    modes.mcReco = doprocessMC || doprocessMCMicro;
    modes.mcRecoMicro = doprocessMCMicro;
    modes.mcGen = doprocessMCTrue;
    core.init(histos, eventCuts, trackCuts, pionPID, kaonPID, secondaryCuts, candidateCuts, histogramOptions, modes);

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in K1 Analysis Task";
    histos.print();
  }

  void processResoTracks(ResoCollisions::iterator const& collision,
                         ResoTracks const& resotracks)
  {
    if (!core.passesEventCuts(collision)) {
      return;
    }
    core.fillHistograms<false, false, false>(histos, collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processResoTracks, "Process ResoTracks", false);

  void processResoMicroTracks(ResoCollisions::iterator const& collision,
                              ResoMicroTracks const& resomicrotracks)
  {
    if (!core.passesEventCuts(collision)) {
      return;
    }
    core.fillHistograms<false, false, true>(histos, collision, resomicrotracks, resomicrotracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processResoMicroTracks, "Process ResoMicroTracks", true);

  void processMC(ResoMCCols::iterator const& collision,
                 ResoMCTracks const& resotracks)
  {
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    histos.fill(HIST("MCReco/collisions"), 0.5);
    core.fillHistograms<true, false, false>(histos, collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processMC, "Process Event for MC", false);

  void processMCMicro(ResoMCCols::iterator const& collision, ResoMCMicroTracks const& tracks)
  {
    // The modular producer already selected these reconstructed collisions.
    // Apply precisely the same reconstruction loop as the data baseline.
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    histos.fill(HIST("MCReco/collisions"), 0.5);
    histos.fill(HIST("MCReco/microTracks"), 0.5, tracks.size());
    core.fillHistograms<true, false, true>(histos, collision, tracks, tracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processMCMicro, "Process reconstructed MC with micro v001 tables", false);

  void processMCTrue(ResoMCCols::iterator const& collision, ResoMCParents const& resoParents)
  {
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    core.fillGenerated(histos, resoParents);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processMCTrue, "Process generated K1 in selected events with v001 parents", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processME(ResoCollisions const& collisions, ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<ResoCollisions, ResoTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (!core.passesEventCuts(collision1) || !core.passesEventCuts(collision2)) {
        continue;
      }
      core.fillHistograms<false, true, false>(histos, collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(K1AnalysisMicro, processME, "Process EventMixing light without partition", false);

  // Processing Event Mixing -- Micro
  void processMEMicro(ResoCollisions const& collisions, ResoMicroTracks const& resomicrotracks)
  {
    auto tracksTuple = std::make_tuple(resomicrotracks);
    BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<ResoCollisions, ResoMicroTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (!core.passesEventCuts(collision1) || !core.passesEventCuts(collision2)) {
        continue;
      }
      core.fillHistograms<false, true, true>(histos, collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(K1AnalysisMicro, processMEMicro, "Process EventMixing light without partition", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<K1AnalysisMicro>(cfgc)};
}
