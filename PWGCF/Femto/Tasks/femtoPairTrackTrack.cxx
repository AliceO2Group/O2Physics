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

/// \file femtoPairTrackTrack.cxx
/// \brief Tasks that computes correlation between two tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/collisionTagger.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/pairCleaner.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>

#include <map>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoPairTrackTrack {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  // opt-in: carries FColTags (see TableProducer/femtoProducerCollisionTagger.cxx), only needed by
  // processMixedEventWithTagging. collisionFilter below (declared on the base FemtoCollisions
  // columns) auto-applies to this joined/filtered type too, exactly as it already does for
  // FilteredFemtoCollisionsWithLabel above; the extra tag-bit condition is then applied via a
  // Partition declared locally inside processMixedEventWithTagging itself (see there), the same
  // "Partition declared and bound inside process()" idiom as Tutorials/src/eventMixing.cxx's
  // MixedEventsPartitionedTracks -- so nothing here needs a second top-level Filter.
  using FemtoCollisionsWithTags = o2::soa::Join<FemtoCollisions, o2::aod::FColTags>;
  using FilteredFemtoCollisionsWithTags = o2::soa::Filtered<FemtoCollisionsWithTags>;

  // collision table with event shape, same filter as standard collisions
  using FemtoCollisionsWithEventShape = o2::soa::Join<FemtoCollisions, o2::aod::FColShapes>;
  using FilteredFemtoCollisionsWithEventShape = o2::soa::Filtered<FemtoCollisionsWithEventShape>;
  using FilteredFemtoCollisionWithEventShape = FilteredFemtoCollisionsWithEventShape::iterator;

  // opt-in: carries both FColShapes and FColTags, for processMixedEventWithTaggingAndEventShape --
  // same collisionFilter auto-apply reasoning as FemtoCollisionsWithTags above.
  using FemtoCollisionsWithEventShapeAndTags = o2::soa::Join<FemtoCollisionsWithEventShape, o2::aod::FColTags>;
  using FilteredFemtoCollisionsWithEventShapeAndTags = o2::soa::Filtered<FemtoCollisionsWithEventShapeAndTags>;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;

  // for analysis which require particles at high pt, add tof mass so sidebands can be used
  using FemtoTracksWithMass = o2::soa::Join<FemtoTracks, o2::aod::FTrackMass>;

  // opt-in: adds the packed purity tables (nSigma, see TrackSelection1/2.puritySlotMask) and FTrackDcas
  // (fraction estimation) for the FitVars histograms. Needs the producer's
  // TrackTables.produceItsPurity/produceTpcPurity/produceTofPurity/produceTrackDcas enabled.
  using FemtoTracksWithFitVars = o2::soa::Join<FemtoTracks, o2::aod::FItsPurity, o2::aod::FTpcPurity, o2::aod::FTofPurity, o2::aod::FTrackDcas>;

  // the Mc process functions always carry the purity/fraction tables (not opt-in like the data
  // path's WithFitVars variant) -- requires the producer's TrackTables.produceItsPurity/
  // produceTpcPurity/produceTofPurity/produceTrackDcas enabled whenever MC processing is used
  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels, o2::aod::FItsPurity, o2::aod::FTpcPurity, o2::aod::FTofPurity, o2::aod::FTrackDcas>;

  using FemtoMcParticlesWithLabel = o2::soa::Join<o2::aod::FMcParticles, o2::aod::FMcMotherLabels>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // opt-in, only needed by processMixedEventWithTagging -- see processMixedEventWithTagging for
  // how it's used.
  collisiontagger::ConfTagging<collisiontagger::role::TrackTrackPair> confTagging;

  // setup tracks
  trackbuilder::ConfTrackSelection1 confTrackSelections1;
  particlecleaner::ConfTrackCleaner1 confTrackCleaner1;
  trackhistmanager::ConfTrackBinning1 confTrackBinning1;

  trackbuilder::ConfTrackSelection2 confTrackSelections2;
  particlecleaner::ConfTrackCleaner2 confTrackCleaner2;
  trackhistmanager::ConfTrackBinning2 confTrackBinning2;

  o2::framework::Partition<FemtoTracks> trackPartition1 = MAKE_TRACK_PARTITION(confTrackSelections1);
  o2::framework::Partition<FemtoTracks> trackPartition2 = MAKE_TRACK_PARTITION(confTrackSelections2);
  o2::framework::Preslice<FemtoTracks> perColtracks = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoTracksWithMass> trackWithMassPartition1 = MAKE_TRACK_PARTITION_WITH_MASS(confTrackSelections1);
  o2::framework::Partition<FemtoTracksWithMass> trackWithMassPartition2 = MAKE_TRACK_PARTITION_WITH_MASS(confTrackSelections2);
  o2::framework::Preslice<FemtoTracksWithMass> perColtracksWithMass = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoTracksWithFitVars> trackWithFitVarsPartition1 = MAKE_TRACK_PARTITION(confTrackSelections1);
  o2::framework::Partition<FemtoTracksWithFitVars> trackWithFitVarsPartition2 = MAKE_TRACK_PARTITION(confTrackSelections2);
  o2::framework::Preslice<FemtoTracksWithFitVars> perColtracksWithFitVars = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition1 = MAKE_TRACK_PARTITION(confTrackSelections1);
  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition2 = MAKE_TRACK_PARTITION(confTrackSelections2);
  o2::framework::Preslice<FemtoTracksWithLabel> perColtracksWithLabel = o2::aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;
  paircleaner::ConfPairCleanerBinning confPairCleaner;

  closepairrejection::ConfCprTrackTrack confCpr;

  pairbuilder::PairTrackTrackBuilder<
    trackhistmanager::PrefixTrack1Se,
    trackhistmanager::PrefixTrack2Se,
    trackhistmanager::PrefixTrack1Me,
    trackhistmanager::PrefixTrack2Me,
    pairhistmanager::PrefixTrackTrackSe,
    pairhistmanager::PrefixTrackTrackMe,
    closepairrejection::PrefixTrackTrackSe,
    closepairrejection::PrefixTrackTrackMe>
    pairTrackTrackBuilder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  std::vector<double> defaultEventPlaneAngleBins{10, 0, o2::constants::math::TwoPI};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent, o2::aod::femtocollisions::EventPlaneAngle> mixBinsVtxCentEventPlaneAngle{{defaultVtxBins, defaultCentBins, defaultEventPlaneAngleBins}};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoTrackTrack", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext& context)
  {
    if ((static_cast<int>(doprocessSameEvent) + static_cast<int>(doprocessSameEventWithMass) + static_cast<int>(doprocessSameEventWithEventShape) + static_cast<int>(doprocessSameEventWithFitVars) + static_cast<int>(doprocessSameEventMc)) > 1 ||
        (static_cast<int>(doprocessMixedEvent) + static_cast<int>(doprocessMixedEventWithMass) + static_cast<int>(doprocessMixedEventWithEventShape) + static_cast<int>(doprocessMixedEventWithFitVars) + static_cast<int>(doprocessMixedEventMc) + static_cast<int>(doprocessMixedEventWithTagging) + static_cast<int>(doprocessMixedEventWithTaggingAndMass) + static_cast<int>(doprocessMixedEventWithTaggingAndEventShape) + static_cast<int>(doprocessMixedEventWithTaggingAndFitVars)) > 1) {
      LOG(fatal) << "More than 1 same or mixed event process function is activated. Breaking...";
    }

    if (doprocessMixedEventWithTagging || doprocessMixedEventWithTaggingAndMass || doprocessMixedEventWithTaggingAndEventShape || doprocessMixedEventWithTaggingAndFitVars) {
      // discoverTaggedRoles() is the same canonical device-list walk femtoProducerCollisionTagger.cxx uses
      // to assign bits -- replaying it here (rather than asking the tagger directly, which has no
      // way to communicate a runtime-computed value back out through RunningWorkflowInfo) is how
      // this task finds out which bit(s) ended up being its own.
      auto const& workflows = context.services().get<o2::framework::RunningWorkflowInfo const>();
      auto const& ownDevice = context.services().get<o2::framework::DeviceSpec const>();
      auto roles = collisiontagger::discoverTaggedRoles(workflows);
      int bit = collisiontagger::findOwnBit(roles, ownDevice.name, confTagging.Role);
      if (bit < 0) {
        LOG(fatal) << "A processMixedEventWithTagging* function is enabled, but no femto-collision-tagger "
                   << "device (with processTrackTrack enabled) tagging this device (" << ownDevice.name
                   << ") was found in this workflow.";
      }
      confTagging.internalTagMask.value = static_cast<datatypes::CollisionTagType>(1) << bit;
      LOG(info) << "Collision tagger bit for " << ownDevice.name << ": " << bit;
    }
    bool processData = doprocessSameEvent || doprocessMixedEvent || doprocessSameEventWithMass || doprocessMixedEventWithMass || doprocessSameEventWithEventShape || doprocessMixedEventWithEventShape || doprocessMixedEventWithTagging || doprocessMixedEventWithTaggingAndMass || doprocessMixedEventWithTaggingAndEventShape;
    bool processFitVars = doprocessSameEventWithFitVars || doprocessMixedEventWithFitVars || doprocessMixedEventWithTaggingAndFitVars;
    bool processMc = doprocessSameEventMc || doprocessMixedEventMc;
    if (static_cast<int>(processData) + static_cast<int>(processFitVars) + static_cast<int>(processMc) > 1) {
      LOG(fatal) << "More than one of data, purity or mc processing is activated. Breaking...";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins.value, confMixing.multBins.value}};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}};
    mixBinsVtxCentEventPlaneAngle = {{confMixing.vtxBins.value, confMixing.centBins.value, confMixing.eventPlaneAngle.value}};

    // setup histogram specs
    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec1;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec2;
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairHistSpec;
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);
    std::map<paircleaner::PairCleanerHist, std::vector<o2::framework::AxisSpec>> pairCleanerHistSpec = paircleaner::makePairCleanerHistSpecMap(confPairCleaner);

    if (processFitVars) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      trackHistSpec1 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning1);
      trackHistSpec2 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning2);
      pairHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
      pairTrackTrackBuilder.init<modes::Mode::kSe_Reco_FitVars, modes::Mode::kMe_Reco_FitVars>(&hRegistry, confCollisionBinning, confTrackSelections1, confTrackSelections2, confTrackCleaner1, confTrackCleaner2, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec1, trackHistSpec2, pairHistSpec, cprHistSpec, pairCleanerHistSpec);
    } else if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      trackHistSpec1 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning1);
      trackHistSpec2 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning2);
      pairHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
      pairTrackTrackBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confTrackSelections1, confTrackSelections2, confTrackCleaner1, confTrackCleaner2, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec1, trackHistSpec2, pairHistSpec, cprHistSpec, pairCleanerHistSpec);
    } else {
      colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      trackHistSpec1 = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning1);
      trackHistSpec2 = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning2);
      pairHistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning, confMixing);
      pairTrackTrackBuilder.init<modes::Mode::kSe_Reco_Mc_FitVars, modes::Mode::kMe_Reco_Mc_FitVars>(&hRegistry, confCollisionBinning, confTrackSelections1, confTrackSelections2, confTrackCleaner1, confTrackCleaner2, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec1, trackHistSpec2, pairHistSpec, cprHistSpec, pairCleanerHistSpec);
    }
    hRegistry.print();
  };

  void processSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks)
  {
    pairTrackTrackBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition1, trackPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processSameEvent, "Enable processing same event processing", true);

  void processSameEventWithMass(FilteredFemtoCollision const& col, FemtoTracksWithMass const& tracks)
  {
    pairTrackTrackBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackWithMassPartition1, trackWithMassPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processSameEventWithMass, "Enable processing same event processing (with track masses)", false);

  void processSameEventWithFitVars(FilteredFemtoCollision const& col, FemtoTracksWithFitVars const& tracks)
  {
    pairTrackTrackBuilder.processSameEvent<modes::Mode::kSe_Reco_FitVars>(col, tracks, trackWithFitVarsPartition1, trackWithFitVarsPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processSameEventWithFitVars, "Enable processing same event processing (with FitVars tables)", false);

  void processSameEventWithEventShape(FilteredFemtoCollisionWithEventShape const& col, FemtoTracks const& tracks)
  {
    pairTrackTrackBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition1, trackPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processSameEventWithEventShape, "Enable processing same event processing with event shape information", false);

  void processSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackTrackBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc_FitVars>(col, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processSameEventMc, "Enable processing same event processing", false);

  // Opt-in counterpart of processMixedEvent() below: `cols` already gets the standard
  // collisionFilter applied (Filtered<FemtoCollisionsWithTags> picks it up automatically, exactly
  // as FilteredFemtoCollisionsWithLabel above does). On top of that, a Partition is declared and
  // bound right here, locally, to pick out only the collisions the tagger found a qualifying pair
  // in -- the same "declare + bindTable() inside process()" idiom as
  // Tutorials/src/eventMixing.cxx's MixedEventsPartitionedTracks, so this needs no change to
  // pairProcessHelpers.h/pairBuilder.h: processMixedEvent is called exactly as below, just with a
  // further-filtered collision view. Requires a femto-collision-tagger device (processTrackTrack
  // enabled) in the same workflow -- see init().
  void processMixedEventWithTagging(FilteredFemtoCollisionsWithTags const& cols, FemtoTracks const& tracks)
  {
    o2::framework::Partition<FilteredFemtoCollisionsWithTags> taggedCollisions = ncheckbit(o2::aod::femtocollisions::collisionTag, confTagging.internalTagMask);
    taggedCollisions.bindTable(cols);
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kMe_Reco>(taggedCollisions, tracks, trackPartition1, trackPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent, mixBinsVtxCentEventPlaneAngle);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEventWithTagging, "Enable processing mixed event processing, pre-filtered using the collision tagger", false);

  // Same tagging filter as processMixedEventWithTagging above, carrying FemtoTracksWithMass
  // instead of plain FemtoTracks -- the tag bit means the same thing regardless of which track
  // table variant is used downstream, since it's computed from the same TrackSelection1/2 cut
  // values either way.
  void processMixedEventWithTaggingAndMass(FilteredFemtoCollisionsWithTags const& cols, FemtoTracksWithMass const& tracks)
  {
    o2::framework::Partition<FilteredFemtoCollisionsWithTags> taggedCollisions = ncheckbit(o2::aod::femtocollisions::collisionTag, confTagging.internalTagMask);
    taggedCollisions.bindTable(cols);
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kMe_Reco>(taggedCollisions, tracks, trackWithMassPartition1, trackWithMassPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent, mixBinsVtxCentEventPlaneAngle);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEventWithTaggingAndMass, "Enable processing mixed event processing, pre-filtered using the collision tagger (with track masses)", false);

  // Same tagging filter as processMixedEventWithTagging above, on a collision table that also
  // carries FColShapes, needed for the kVtxCentEventPlaneAngle mixing policy.
  void processMixedEventWithTaggingAndEventShape(FilteredFemtoCollisionsWithEventShapeAndTags const& cols, FemtoTracks const& tracks)
  {
    o2::framework::Partition<FilteredFemtoCollisionsWithEventShapeAndTags> taggedCollisions = ncheckbit(o2::aod::femtocollisions::collisionTag, confTagging.internalTagMask);
    taggedCollisions.bindTable(cols);
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kMe_Reco>(taggedCollisions, tracks, trackPartition1, trackPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent, mixBinsVtxCentEventPlaneAngle);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEventWithTaggingAndEventShape, "Enable processing mixed event processing, pre-filtered using the collision tagger, with event shape information", false);

  // Same tagging filter as processMixedEventWithTagging above, carrying FemtoTracksWithFitVars
  // instead of plain FemtoTracks -- same rationale as processMixedEventWithTaggingAndMass above.
  void processMixedEventWithTaggingAndFitVars(FilteredFemtoCollisionsWithTags const& cols, FemtoTracksWithFitVars const& tracks)
  {
    o2::framework::Partition<FilteredFemtoCollisionsWithTags> taggedCollisions = ncheckbit(o2::aod::femtocollisions::collisionTag, confTagging.internalTagMask);
    taggedCollisions.bindTable(cols);
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kMe_Reco_FitVars>(taggedCollisions, tracks, trackWithFitVarsPartition1, trackWithFitVarsPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent, mixBinsVtxCentEventPlaneAngle);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEventWithTaggingAndFitVars, "Enable processing mixed event processing, pre-filtered using the collision tagger (with FitVars tables)", false);

  void processMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks)
  {
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, trackPartition1, trackPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent, mixBinsVtxCentEventPlaneAngle);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEvent, "Enable processing mixed event processing", true);

  void processMixedEventWithMass(FilteredFemtoCollisions const& cols, FemtoTracksWithMass const& tracks)
  {
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, trackWithMassPartition1, trackWithMassPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent, mixBinsVtxCentEventPlaneAngle);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEventWithMass, "Enable processing mixed event processing (with track masses)", false);

  void processMixedEventWithFitVars(FilteredFemtoCollisions const& cols, FemtoTracksWithFitVars const& tracks)
  {
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kMe_Reco_FitVars>(cols, tracks, trackWithFitVarsPartition1, trackWithFitVarsPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent, mixBinsVtxCentEventPlaneAngle);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEventWithFitVars, "Enable processing mixed event processing (with FitVars tables)", false);

  void processMixedEventWithEventShape(FilteredFemtoCollisionsWithEventShape const& cols, FemtoTracks const& tracks)
  {
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, trackPartition1, trackPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent, mixBinsVtxCentEventPlaneAngle);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEventWithEventShape, "Enable processing mixed event processing with event shape information", false);

  void processMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc_FitVars>(cols, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent, mixBinsVtxCentEventPlaneAngle);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEventMc, "Enable processing mixed event processing", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackTrack>(context),
  };
  return workflow;
}
