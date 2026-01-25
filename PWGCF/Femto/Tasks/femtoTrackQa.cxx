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

#include <map>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoTrackQa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks, o2::aod::FColPos, o2::aod::FColSphericities, o2::aod::FColMults>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks, o2::aod::FTrackDcas, o2::aod::FTrackExtras, o2::aod::FTrackPids>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;
  colhistmanager::CollisionHistManager colHistManager;

  // setup tracks
  trackbuilder::ConfTrackSelection1 confTrackSelection;
  trackhistmanager::ConfTrackBinning1 confTrackBinning;
  trackhistmanager::ConfTrackQaBinning1 confTrackQaBinning;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrackQa> trackHistManager;

  o2::framework::Partition<FemtoTracks> trackPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  o2::framework::Preslice<FemtoTracks> perColReco = o2::aod::femtobase::stored::fColId;

  particlecleaner::ConfTrackCleaner1 confTrackCleaner;
  particlecleaner::ParticleCleaner trackCleaner;

  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  o2::framework::Preslice<FemtoTracksWithLabel> perColRecoWithLabel = o2::aod::femtobase::stored::fColId;

  o2::framework::HistogramRegistry hRegistry{"FemtoTrackQA", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    if ((doprocessData + doprocessMc) > 1) {
      LOG(fatal) << "More than 1 process function is activated. Breaking...";
    }
    bool processData = doprocessData;
    trackCleaner.init(confTrackCleaner);

    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec;

    if (processData) {
      colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, colHistSpec, confCollisionQaBinning);
      trackHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confTrackBinning, confTrackQaBinning);
      trackHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, trackHistSpec, confTrackSelection, confTrackQaBinning);
    } else {
      colHistSpec = colhistmanager::makeColMcQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, colHistSpec, confCollisionQaBinning);
      trackHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confTrackBinning, confTrackQaBinning);
      trackHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, trackHistSpec, confTrackSelection, confTrackQaBinning);
    }
    hRegistry.print();
  };

  void processData(FilteredFemtoCollision const& col, FemtoTracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col);
    auto trackSlice = trackPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& track : trackSlice) {
      trackHistManager.fill<modes::Mode::kAnalysis_Qa>(track, tracks);
    }
  };
  PROCESS_SWITCH(FemtoTrackQa, processData, "Track QA in Data", true);

  void processMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(col, mcCols);
    auto trackSlice = trackWithLabelPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& track : trackSlice) {
      if (!trackCleaner.isClean(track, mcParticles, mcMothers, mcPartonicMothers)) {
        continue;
      }
      trackHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(track, tracks, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  PROCESS_SWITCH(FemtoTrackQa, processMc, "Track QA in Monte Carlo", false);
};

o2::framework::WorkflowSpec
  defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoTrackQa>(cfgc),
  };
  return workflow;
}
