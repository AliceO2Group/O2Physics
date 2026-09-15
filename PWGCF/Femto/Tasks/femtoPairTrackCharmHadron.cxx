// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoPairTrackCharmHadron.cxx
/// \brief Tasks that computes correlation between tracks and charm hadrons
/// \author Igor Ptak, WUT, igor.tomasz.ptak@cern.ch

#include "PWGCF/Femto/Core/charmHadronBuilder.h"
#include "PWGCF/Femto/Core/charmHadronHistManager.h"
#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/pairCleaner.h"
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

struct FemtoPairTrackCharmHadron {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoD0s = o2::soa::Join<o2::aod::FD0s, o2::aod::FD0Masks>;
  using FemtoLcs = o2::soa::Join<o2::aod::FLcs, o2::aod::FLcMasks>;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;
  using FemtoD0sWithLabel = o2::soa::Join<FemtoD0s, o2::aod::FD0Labels>;
  using FemtoLcsWithLabel = o2::soa::Join<FemtoLcs, o2::aod::FLcLabels>;
  using FemtoMcParticlesWithLabel = o2::soa::Join<o2::aod::FMcParticles, o2::aod::FMcMotherLabels>;

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

  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  o2::framework::Preslice<FemtoTracksWithLabel> perColTracksWithLabel = o2::aod::femtobase::stored::fColId;

  // setup for D0 daughters
  trackhistmanager::ConfD01PosDauBinning confPosDauBinning;
  trackhistmanager::ConfD01NegDauBinning confNegDauBinning;

  // setup D0s
  charmhadronbuilder::ConfD0Selection1 d0Selection;
  charmhadronhistmanager::ConfD0Binning1 confD0Binning;
  particlecleaner::ConfD0Cleaner1 confD0Cleaner;

  o2::framework::Partition<FemtoD0s> d0Partition = MAKE_D0_PARTITION(d0Selection);
  o2::framework::Preslice<FemtoD0s> perColD0s = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoD0sWithLabel> d0WithLabelPartition = MAKE_D0_PARTITION(d0Selection);
  o2::framework::Preslice<FemtoD0sWithLabel> perColD0sWithLabel = o2::aod::femtobase::stored::fColId;

  // setup for Lc daughters
  trackhistmanager::ConfLc1ProtonDauBinning confProtonDauBinning;
  trackhistmanager::ConfLc1KaonDauBinning confKaonDauBinning;
  trackhistmanager::ConfLc1PionDauBinning confPionDauBinning;

  // setup Lcs
  charmhadronbuilder::ConfLcSelection1 lcSelection;
  charmhadronhistmanager::ConfLcBinning1 confLcBinning;
  particlecleaner::ConfLcCleaner1 confLcCleaner;

  o2::framework::Partition<FemtoLcs> lcPartition = MAKE_CHARM3PRONG_PARTITION(lcSelection);
  o2::framework::Preslice<FemtoLcs> perColLcs = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoLcsWithLabel> lcWithLabelPartition = MAKE_CHARM3PRONG_PARTITION(lcSelection);
  o2::framework::Preslice<FemtoLcsWithLabel> perColLcsWithLabel = o2::aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;
  paircleaner::ConfPairCleanerBinning confPairCleaner;

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

  pairbuilder::PairTrackLcBuilder<
    trackhistmanager::PrefixTrack1,
    charmhadronhistmanager::PrefixLc1,
    trackhistmanager::PrefixLc1ProtonDaughter,
    trackhistmanager::PrefixLc1KaonDaughter,
    trackhistmanager::PrefixLc1PionDaughter,
    pairhistmanager::PrefixTrackLcSe,
    pairhistmanager::PrefixTrackLcMe,
    closepairrejection::PrefixTrackLcProtonSe,
    closepairrejection::PrefixTrackLcKaonSe,
    closepairrejection::PrefixTrackLcPionSe,
    closepairrejection::PrefixTrackLcProtonMe,
    closepairrejection::PrefixTrackLcKaonMe,
    closepairrejection::PrefixTrackLcPionMe,
    modes::CharmHadron::kLc>
    pairTrackLcBuilder;

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
  closepairrejection::ConfCprTrackLcProton confCprLcProton;
  closepairrejection::ConfCprTrackLcKaon confCprLcKaon;
  closepairrejection::ConfCprTrackLcPion confCprLcPion;

  void init(o2::framework::InitContext&)
  {
    bool processD0 = doprocessD0SameEvent || doprocessD0MixedEvent;
    bool processD0Mc = doprocessD0SameEventMc || doprocessD0MixedEventMc;
    bool processLc = doprocessLcSameEvent || doprocessLcMixedEvent;
    bool processLcMc = doprocessLcSameEventMc || doprocessLcMixedEventMc;

    if ((static_cast<int>(processD0) + static_cast<int>(processD0Mc) + static_cast<int>(processLc) + static_cast<int>(processLcMc)) > 1) {
      LOG(fatal) << "Only one charm hadron species and data/mc mode can be processed at a time. Breaking...";
    }

    bool processData = processD0 || processLc;

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins.value, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> protonDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> kaonDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> pionDauSpec;
    std::map<charmhadronhistmanager::CharmHadronHist, std::vector<o2::framework::AxisSpec>> d0HistSpec;
    std::map<charmhadronhistmanager::CharmHadronHist, std::vector<o2::framework::AxisSpec>> lcHistSpec;
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairTrackD0HistSpec;
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairTrackLcHistSpec;
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecLc = closepairrejection::makeCprHistSpecMap(confCprLcProton);
    std::map<paircleaner::PairCleanerHist, std::vector<o2::framework::AxisSpec>> pairCleanerHistSpec = paircleaner::makePairCleanerHistSpecMap(confPairCleaner);

    if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      trackHistSpec = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
      if (processD0) {
        posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
        negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
        d0HistSpec = charmhadronhistmanager::makeCharmHadronHistSpecMap(confD0Binning);
        pairTrackD0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
        pairTrackD0Builder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, d0Selection, confD0Cleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, d0HistSpec, posDauSpec, negDauSpec, pairTrackD0HistSpec, cprHistSpec, pairCleanerHistSpec);
      }
      if (processLc) {
        protonDauSpec = trackhistmanager::makeTrackHistSpecMap(confProtonDauBinning);
        kaonDauSpec = trackhistmanager::makeTrackHistSpecMap(confKaonDauBinning);
        pionDauSpec = trackhistmanager::makeTrackHistSpecMap(confPionDauBinning);
        lcHistSpec = charmhadronhistmanager::makeCharmHadronHistSpecMap(confLcBinning);
        pairTrackLcHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
        pairTrackLcBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, lcSelection, confLcCleaner, confCprLcProton, confCprLcKaon, confCprLcPion, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, lcHistSpec, protonDauSpec, kaonDauSpec, pionDauSpec, pairTrackLcHistSpec, cprHistSpecLc, pairCleanerHistSpec);
      }
    } else {
      colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      trackHistSpec = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning);
      if (processD0Mc) {
        posDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confPosDauBinning);
        negDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confNegDauBinning);
        d0HistSpec = charmhadronhistmanager::makeCharmHadronMcHistSpecMap(confD0Binning);
        pairTrackD0HistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning, confMixing);
        pairTrackD0Builder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, d0Selection, confD0Cleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, d0HistSpec, posDauSpec, negDauSpec, pairTrackD0HistSpec, cprHistSpec, pairCleanerHistSpec);
      }
      if (processLcMc) {
        protonDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confProtonDauBinning);
        kaonDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confKaonDauBinning);
        pionDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confPionDauBinning);
        lcHistSpec = charmhadronhistmanager::makeCharmHadronMcHistSpecMap(confLcBinning);
        pairTrackLcHistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning, confMixing);
        pairTrackLcBuilder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, lcSelection, confLcCleaner, confCprLcProton, confCprLcKaon, confCprLcPion, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, lcHistSpec, protonDauSpec, kaonDauSpec, pionDauSpec, pairTrackLcHistSpec, cprHistSpecLc, pairCleanerHistSpec);
      }
    }

    hRegistry.print();
  };

  void processD0SameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoD0s const& d0s)
  {
    pairTrackD0Builder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition, d0s, d0Partition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCharmHadron, processD0SameEvent, "Enable processing same event processing for tracks and D0s", true);

  void processD0MixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoD0s const& /*d0s*/)
  {
    pairTrackD0Builder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, trackPartition, d0Partition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCharmHadron, processD0MixedEvent, "Enable processing mixed event processing for tracks and D0s", true);

  void processD0SameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoD0sWithLabel const& d0s, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackD0Builder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, trackWithLabelPartition, d0s, d0WithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCharmHadron, processD0SameEventMc, "Enable processing same event processing for tracks and D0s with MC information", false);

  void processD0MixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoD0sWithLabel const& /*d0s*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackD0Builder.processMixedEvent<modes::Mode::kMe_Reco_Mc>(cols, mcCols, tracks, trackWithLabelPartition, d0WithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCharmHadron, processD0MixedEventMc, "Enable processing mixed event processing for tracks and D0s with MC information", false);

  void processLcSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLcs const& lcs)
  {
    pairTrackLcBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition, lcs, lcPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCharmHadron, processLcSameEvent, "Enable processing same event processing for tracks and Lcs", false);

  void processLcMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoLcs const& /*lcs*/)
  {
    pairTrackLcBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, trackPartition, lcPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCharmHadron, processLcMixedEvent, "Enable processing mixed event processing for tracks and Lcs", false);

  void processLcSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLcsWithLabel const& lcs, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackLcBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, trackWithLabelPartition, lcs, lcWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCharmHadron, processLcSameEventMc, "Enable processing same event processing for tracks and Lcs with MC information", false);

  void processLcMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLcsWithLabel const& /*lcs*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackLcBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc>(cols, mcCols, tracks, trackWithLabelPartition, lcWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCharmHadron, processLcMixedEventMc, "Enable processing mixed event processing for tracks and Lcs with MC information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackCharmHadron>(context),
  };
  return workflow;
}
