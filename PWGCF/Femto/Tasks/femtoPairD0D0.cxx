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

/// \file femtoPairD0D0.cxx
/// \brief Tasks that computes correlation between two D0 mesons
/// \author Igor Ptak, WUT, igor.tomasz.ptak@cern.ch

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

// Each pair slot (1 and 2) has its own selection, cleaner, binning and partition.
// The sign of the selection picks the decay hypothesis:
//   d0-d0       : D0Selection1.sign = +1, D0Selection2.sign = +1, Mixing.sameSpecies = true
//   d0bar-d0bar : D0Selection1.sign = -1, D0Selection2.sign = -1, Mixing.sameSpecies = true
//   d0-d0bar    : D0Selection1.sign = +1, D0Selection2.sign = -1, Mixing.sameSpecies = false
struct FemtoPairD0D0 {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::aod::FTracks;
  using FemtoD0s = o2::soa::Join<o2::aod::FD0s, o2::aod::FD0Masks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;
  using FemtoD0sWithLabel = o2::soa::Join<FemtoD0s, o2::aod::FD0Labels>;

  using FemtoMcParticlesWithLabel = o2::soa::Join<o2::aod::FMcParticles, o2::aod::FMcMotherLabels>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup for daughters
  // separate binning per d0 slot, so that both hypotheses can be binned independently
  trackhistmanager::ConfD01PosDauBinning confD01PosDauBinning;
  trackhistmanager::ConfD01NegDauBinning confD01NegDauBinning;
  trackhistmanager::ConfD02PosDauBinning confD02PosDauBinning;
  trackhistmanager::ConfD02NegDauBinning confD02NegDauBinning;

  // setup d0s
  charmhadronbuilder::ConfD0Selection1 confD0Selection1;
  charmhadronbuilder::ConfD0Selection2 confD0Selection2;
  particlecleaner::ConfD0Cleaner1 confD0Cleaner1;
  particlecleaner::ConfD0Cleaner2 confD0Cleaner2;
  charmhadronhistmanager::ConfD0Binning1 confD0Binning1;
  charmhadronhistmanager::ConfD0Binning2 confD0Binning2;

  o2::framework::Partition<FemtoD0s> d0Partition1 = MAKE_D0_PARTITION(confD0Selection1);
  o2::framework::Partition<FemtoD0s> d0Partition2 = MAKE_D0_PARTITION(confD0Selection2);
  o2::framework::Preslice<FemtoD0s> perColD0s = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoD0sWithLabel> d0WithLabelPartition1 = MAKE_D0_PARTITION(confD0Selection1);
  o2::framework::Partition<FemtoD0sWithLabel> d0WithLabelPartition2 = MAKE_D0_PARTITION(confD0Selection2);
  o2::framework::Preslice<FemtoD0sWithLabel> perColD0sWithLabel = o2::aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;

  pairbuilder::PairD0D0Builder<
    charmhadronhistmanager::PrefixD01,
    trackhistmanager::PrefixD01PosDaughter,
    trackhistmanager::PrefixD01NegDaughter,
    charmhadronhistmanager::PrefixD02,
    trackhistmanager::PrefixD02PosDaughter,
    trackhistmanager::PrefixD02NegDaughter,
    pairhistmanager::PrefixD0D0Se,
    pairhistmanager::PrefixD0D0Me,
    closepairrejection::PrefixD0D0PosSe,
    closepairrejection::PrefixD0D0NegSe,
    closepairrejection::PrefixD0D0PosMe,
    closepairrejection::PrefixD0D0NegMe,
    modes::CharmHadron::kD0,
    modes::CharmHadron::kD0>
    pairD0D0Builder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoD0D0", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  closepairrejection::ConfCprD0DaugherD0DaughterPos confCprPos;
  closepairrejection::ConfCprD0DaugherD0DaughterNeg confCprNeg;

  void init(o2::framework::InitContext&)
  {
    bool processData = doprocessSameEvent || doprocessMixedEvent;
    bool processMc = doprocessSameEventMc || doprocessMixedEventMc;

    if (processData && processMc) {
      LOG(fatal) << "Both data and mc processing is enabled. Breaking...";
    }
    if (!processData && !processMc) {
      LOG(fatal) << "Neither data nor mc processing is enabled. Breaking...";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins.value, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histograms
    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDauSpec1;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDauSpec1;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDauSpec2;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDauSpec2;
    std::map<charmhadronhistmanager::CharmHadronHist, std::vector<o2::framework::AxisSpec>> d0HistSpec1;
    std::map<charmhadronhistmanager::CharmHadronHist, std::vector<o2::framework::AxisSpec>> d0HistSpec2;
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairD0D0HistSpec;
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecPos = closepairrejection::makeCprHistSpecMap(confCprPos);
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecNeg = closepairrejection::makeCprHistSpecMap(confCprNeg);

    if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      posDauSpec1 = trackhistmanager::makeTrackHistSpecMap(confD01PosDauBinning);
      negDauSpec1 = trackhistmanager::makeTrackHistSpecMap(confD01NegDauBinning);
      posDauSpec2 = trackhistmanager::makeTrackHistSpecMap(confD02PosDauBinning);
      negDauSpec2 = trackhistmanager::makeTrackHistSpecMap(confD02NegDauBinning);
      d0HistSpec1 = charmhadronhistmanager::makeD0HistSpecMap(confD0Binning1);
      d0HistSpec2 = charmhadronhistmanager::makeD0HistSpecMap(confD0Binning2);
      pairD0D0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
      pairD0D0Builder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confD0Selection1, confD0Selection2, confD0Cleaner1, confD0Cleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, d0HistSpec1, d0HistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairD0D0HistSpec, cprHistSpecPos, cprHistSpecNeg);
    } else {
      colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      posDauSpec1 = trackhistmanager::makeTrackMcHistSpecMap(confD01PosDauBinning);
      negDauSpec1 = trackhistmanager::makeTrackMcHistSpecMap(confD01NegDauBinning);
      posDauSpec2 = trackhistmanager::makeTrackMcHistSpecMap(confD02PosDauBinning);
      negDauSpec2 = trackhistmanager::makeTrackMcHistSpecMap(confD02NegDauBinning);
      d0HistSpec1 = charmhadronhistmanager::makeD0McHistSpecMap(confD0Binning1);
      d0HistSpec2 = charmhadronhistmanager::makeD0McHistSpecMap(confD0Binning2);
      pairD0D0HistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning, confMixing);
      pairD0D0Builder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistry, confCollisionBinning, confD0Selection1, confD0Selection2, confD0Cleaner1, confD0Cleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, d0HistSpec1, d0HistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairD0D0HistSpec, cprHistSpecPos, cprHistSpecNeg);
    }
  }

  void processSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoD0s const& d0s)
  {
    pairD0D0Builder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, d0s, d0Partition1, d0Partition2, cache);
  }
  PROCESS_SWITCH(FemtoPairD0D0, processSameEvent, "Enable processing same event processing for d0-d0", true);

  void processSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoD0sWithLabel const& d0s, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairD0D0Builder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, d0s, d0WithLabelPartition1, d0WithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairD0D0, processSameEventMc, "Enable processing same event processing for d0-d0 with mc information", false);

  void processMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoD0s const& /*d0s*/)
  {
    pairD0D0Builder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, d0Partition1, d0Partition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairD0D0, processMixedEvent, "Enable processing mixed event processing for d0-d0", true);

  void processMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoD0sWithLabel const& /*d0s*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairD0D0Builder.processMixedEvent<modes::Mode::kMe_Reco_Mc>(cols, mcCols, tracks, d0WithLabelPartition1, d0WithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairD0D0, processMixedEventMc, "Enable processing mixed event processing for d0-d0 with mc information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairD0D0>(context),
  };
  return workflow;
}
