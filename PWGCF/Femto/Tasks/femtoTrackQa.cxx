// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoTrackQa.cxx
/// \brief QA task for tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/OutputObjHeader.h"
#include "Framework/runDataProcessing.h"

using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

struct FemtoTrackQa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<FCols, FColMasks, FColPos, FColSphericities, FColMults>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<FTracks, FTrackMasks, FTrackDcas, FTrackExtras, FTrackPids>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, FTrackLabels>;

  SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;
  colhistmanager::CollisionHistManager colHistManager;

  // setup tracks
  trackbuilder::ConfTrackSelection1 confTrackSelection;
  trackhistmanager::ConfTrackBinning1 confTrackBinning;
  trackhistmanager::ConfTrackQaBinning1 confTrackQaBinning;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrackQa> trackHistManager;

  Partition<FemtoTracks> trackPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  Preslice<FemtoTracks> perColReco = femtobase::stored::fColId;

  particlecleaner::ConfTrackCleaner1 confTrackCleaner;
  particlecleaner::ParticleCleaner trackCleaner;

  Partition<FemtoTracksWithLabel> trackWithLabelPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  Preslice<FemtoTracksWithLabel> perColRecoWithLabel = femtobase::stored::fColId;

  HistogramRegistry hRegistry{"FemtoTrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  template <modes::Mode mode,
            typename MakeColSpec,
            typename MakeTrackSpec>
  void initMode(MakeColSpec&& makeColSpec,
                MakeTrackSpec&& makeTrackSpec)
  {
    auto colHistSpec = makeColSpec(confCollisionBinning, confCollisionQaBinning);
    colHistManager.init<mode>(&hRegistry, colHistSpec, confCollisionQaBinning);
    auto trackHistSpec = makeTrackSpec(confTrackBinning, confTrackQaBinning);

    trackHistManager.init<mode>(
      &hRegistry,
      trackHistSpec,
      confTrackSelection.chargeAbs.value,
      confTrackSelection.chargeSign.value,
      confTrackSelection.pdgCodeAbs.value,
      confTrackQaBinning);
  }

  void init(InitContext&)
  {
    if ((doprocessData + doprocessMc) > 1) {
      LOG(fatal) << "More than 1 process function is activated. Breaking...";
    }
    bool processData = doprocessData;
    trackCleaner.init(confTrackCleaner);
    if (processData) {
      auto colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, colHistSpec, confCollisionQaBinning);
      auto trackHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confTrackBinning, confTrackQaBinning);
      trackHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, trackHistSpec, confTrackSelection, confTrackQaBinning);
    } else {
      auto colHistSpec = colhistmanager::makeColMcQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, colHistSpec, confCollisionQaBinning);
      auto trackHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confTrackBinning, confTrackQaBinning);
      trackHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, trackHistSpec, confTrackSelection, confTrackQaBinning);
    }
    hRegistry.print();
  };

  void processData(FilteredFemtoCollision const& col, FemtoTracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col);
    auto trackSlice = trackPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& track : trackSlice) {
      trackHistManager.fill<modes::Mode::kAnalysis_Qa>(track, tracks);
    }
  };
  PROCESS_SWITCH(FemtoTrackQa, processData, "Track QA in Data", true);

  void processMc(FilteredFemtoCollisionWithLabel const& col, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(col, mcCols);
    auto trackSlice = trackWithLabelPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& track : trackSlice) {
      if (!trackCleaner.isClean(track, mcParticles, mcMothers, mcPartonicMothers)) {
        continue;
      }
      trackHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(track, tracks, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  PROCESS_SWITCH(FemtoTrackQa, processMc, "Track QA in Monte Carlo", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoTrackQa>(cfgc),
  };
  return workflow;
}
