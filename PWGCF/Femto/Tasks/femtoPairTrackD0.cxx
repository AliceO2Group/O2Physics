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

/// \file femtoPairTrackD0.cxx
/// \brief Tasks that computes correlation between tracks and D0 mesons
/// \author igor.ptak.stud@pw.edu.pl, WUT, igor.ptak.stud@pw.edu.pl

#include "PWGCF/Femto/Core/charmHadronBuilder.h"
#include "PWGCF/Femto/Core/charmHadronHistManager.h"
#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <map>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoPairTrackD0 {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoD0s = o2::soa::Join<o2::aod::FD0s, o2::aod::FD0Masks>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup tracks
  trackbuilder::ConfTrackSelection1 confTrackSelection;
  trackhistmanager::ConfTrackBinning1 confTrackBinning;
  particlecleaner::ConfTrackCleaner1 confTrackCleaner;

  o2::framework::Partition<FemtoTracks> trackPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  o2::framework::Preslice<FemtoTracks> perColTracks = o2::aod::femtobase::stored::fColId;

  // setup for D0 daughters
  trackhistmanager::ConfD0PosDauBinning confPosDauBinning;
  trackhistmanager::ConfD0NegDauBinning confNegDauBinning;

  // setup D0s
  charmhadronbuilder::ConfD0Selection d0Selection;
  charmhadronhistmanager::ConfD0Binning1 confD0Binning;
  particlecleaner::ConfD0Cleaner1 confD0Cleaner;

  o2::framework::Partition<FemtoD0s> d0Partition = MAKE_D0_PARTITION(d0Selection);
  o2::framework::Preslice<FemtoD0s> perColD0s = o2::aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;

  pairbuilder::PairTrackD0Builder<
    trackhistmanager::PrefixTrack1,
    charmhadronhistmanager::PrefixD01,
    trackhistmanager::PrefixD01PosDaughter,
    trackhistmanager::PrefixD01NegDaughter,
    pairhistmanager::PrefixTrackD0Se,
    pairhistmanager::PrefixTrackD0Me,
    closepairrejection::PrefixTrackD0DaughterSe,
    closepairrejection::PrefixTrackD0DaughterMe,
    modes::CharmHadron::kD0>
    pairTrackD0Builder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoTrackD0", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  closepairrejection::ConfCprTrackD0Daughter confCpr;

  void init(o2::framework::InitContext&)
  {
    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
    std::map<charmhadronhistmanager::CharmHadronHist, std::vector<o2::framework::AxisSpec>> d0HistSpec = charmhadronhistmanager::makeD0HistSpecMap(confD0Binning);
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairTrackD0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);

    pairTrackD0Builder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, d0Selection, confD0Cleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, d0HistSpec, posDauSpec, negDauSpec, pairTrackD0HistSpec, cprHistSpec);

    hRegistry.print();
  };

  void processSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoD0s const& d0s)
  {
    pairTrackD0Builder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition, d0s, d0Partition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackD0, processSameEvent, "Enable processing same event processing for tracks and D0s", true);

  void processMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoD0s const& /*d0s*/)
  {
    pairTrackD0Builder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, trackPartition, d0Partition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackD0, processMixedEvent, "Enable processing mixed event processing for tracks and D0s", true);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackD0>(context),
  };
  return workflow;
}
