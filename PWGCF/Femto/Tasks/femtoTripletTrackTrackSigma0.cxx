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

/// \file femtoTripletTrackTrackV0.cxx
/// \brief Tasks that computes correlation between two tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/closeTripletRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairCleaner.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/tripletBuilder.h"
#include "PWGCF/Femto/Core/tripletHistManager.h"
#include "PWGCF/Femto/Core/sigma0Builder.h"
#include "PWGCF/Femto/Core/sigma0HistManager.h"
#include "PWGCF/Femto/Core/v0HistManager.h"
#include "PWGCF/Femto/Core/photonHistManager.h"
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
#include <string>

using namespace o2::analysis::femto;

struct FemtoTripletTrackTrackSigma0 {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoSigma0s = o2::soa::Join<o2::aod::FSigma0s, o2::aod::FSigma0Masks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;
  using FemtoSigma0sWithLabel = o2::soa::Join<FemtoSigma0s, o2::aod::FSigma0Labels>;

  using FemtoMcParticlesWithLabel = o2::soa::Join<o2::aod::FMcParticles, o2::aod::FMcMotherLabels>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup tracks
  trackbuilder::ConfTrackSelection1 confTrackSelections1;
  trackhistmanager::ConfTrackBinning1 confTrackBinning1;
  trackbuilder::ConfTrackSelection2 confTrackSelections2;
  trackhistmanager::ConfTrackBinning2 confTrackBinning2;
  inline static constexpr std::string_view  PrefixTrack1 = "Track1/";
  inline static constexpr std::string_view  PrefixTrack2 = "Track2/";

  // particle cleaners (mc only)
  particlecleaner::ConfTrackCleaner1 confTrackCleaner1;
  particlecleaner::ConfTrackCleaner2 confTrackCleaner2;
  particlecleaner::ConfSigma0Cleaner confSigma0Cleaner;

  o2::framework::Partition<FemtoTracks> trackPartition1 = MAKE_TRACK_PARTITION(confTrackSelections1);
  o2::framework::Partition<FemtoTracks> trackPartition2 = MAKE_TRACK_PARTITION(confTrackSelections2);
  o2::framework::Preslice<FemtoTracks> perColtracks = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition1 = MAKE_TRACK_PARTITION(confTrackSelections1);
  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition2 = MAKE_TRACK_PARTITION(confTrackSelections2);
  o2::framework::Preslice<FemtoTracksWithLabel> perColtracksWithLabel = o2::aod::femtobase::stored::fColId;

  // setup for photon and lambda
  v0histmanager::ConfLambdaBinning1 confLambdaBinning;
  inline static constexpr std::string_view  PrefixPhoton = "Sigma0/DauPhoton/";
  photonhistmanager::ConfPhotonBinning confPhotonBinning;
  inline static constexpr std::string_view  PrefixLambda = "Sigma0/DauLambda/";

  // setup sigma0s
  sigma0builder::ConfSigma0SelectionDefault confSigma0Selection;
  sigma0histmanager::ConfSigma0BinningDefault confSigma0Binning;

  o2::framework::Partition<FemtoSigma0s> sigma0Partition = MAKE_SIGMA0_PARTITION(confSigma0Selection);
  o2::framework::Preslice<FemtoSigma0s> perColSigma0s = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoSigma0sWithLabel> sigma0WithLabelPartition = MAKE_SIGMA0_PARTITION(confSigma0Selection);
  o2::framework::Preslice<FemtoSigma0sWithLabel> perColSigma0sWithLabel = o2::aod::femtobase::stored::fColId;

  // setup triplets
  triplethistmanager::ConfTripletBinning confTripletBinning;
  triplethistmanager::ConfTripletCuts confTripletCuts;
  paircleaner::ConfPairCleanerBinning confTripletCleanerBinning;

  closetripletrejection::ConfCtrTrackTrackV0 confCtr;
  closepairrejection::ConfCprTrackV0Daughter confCpr;

  tripletbuilder::TripletTrackTrackSigma0Builder<
    PrefixTrack1,
    PrefixTrack2,
    sigma0histmanager::PrefixSigma0,
    PrefixLambda,
    PrefixPhoton,
    triplethistmanager::PrefixTrackTrackSigma0Se,
    triplethistmanager::PrefixTrackTrackSigma0Me,
    closetripletrejection::PrefixTrack1Track2Se,
    closetripletrejection::PrefixTrack1Sigma0Se,
    closetripletrejection::PrefixTrack2Sigma0Se,
    closetripletrejection::PrefixTrack1Track2Me,
    closetripletrejection::PrefixTrack1Sigma0Me,
    closetripletrejection::PrefixTrack2Sigma0Me>
    tripletTrackTrackSigma0Builder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}};
  triplethistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoTrackTrackSigma0", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    if ((static_cast<int>(doprocessSameEvent) + static_cast<int>(doprocessSameEventMc)) > 1 || (static_cast<int>(doprocessMixedEvent) + static_cast<int>(doprocessMixedEventMc)) > 1) {
      LOG(fatal) << "More than 1 same or mixed event process function is activated. Breaking...";
    }
    bool processData = doprocessSameEvent || doprocessMixedEvent;
    bool processMc = doprocessSameEventMc || doprocessMixedEventMc;
    if (processData && processMc) {
      LOG(fatal) << "Both data and mc processing is activated. Breaking...";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins.value, confMixing.multBins.value}};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}};

    // setup histogram specs
    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec1;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec2;

    std::map<v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> lambdaSpec;
    std::map<photonhistmanager::PhotonHist, std::vector<o2::framework::AxisSpec>> photonSpec;
    std::map<sigma0histmanager::Sigma0Hist, std::vector<o2::framework::AxisSpec>> sigma0HistSpec;
    std::map<triplethistmanager::TripletHist, std::vector<o2::framework::AxisSpec>> tripletHistSpec;
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> ctrHistSpec = closepairrejection::makeCprHistSpecMap(confCtr);
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);
    std::map<paircleaner::PairCleanerHist, std::vector<o2::framework::AxisSpec>> tripletCleanerHistSpec = paircleaner::makePairCleanerHistSpecMap(confTripletCleanerBinning);

    if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      trackHistSpec1 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning1);
      trackHistSpec2 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning2);
      sigma0HistSpec = sigma0histmanager::makeSigma0HistSpecMap(confSigma0Binning);
      lambdaSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
      photonSpec = photonhistmanager::makePhotonHistSpecMap(confPhotonBinning);
      tripletHistSpec = triplethistmanager::makeTripletHistSpecMap(confTripletBinning, confMixing);
      tripletTrackTrackSigma0Builder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confTrackSelections1, confTrackSelections2, confSigma0Selection, confTrackCleaner1, confTrackCleaner2, confSigma0Cleaner, confCtr, confMixing, confTripletBinning, confTripletCuts, colHistSpec, trackHistSpec1, trackHistSpec2, sigma0HistSpec, tripletHistSpec, ctrHistSpec, tripletCleanerHistSpec);
    } else {
    }
    hRegistry.print();
  };

  void processSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoSigma0s const& /*sigmas*/)
  {
    tripletTrackTrackSigma0Builder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition1, trackPartition2, sigma0Partition, cache);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackSigma0, processSameEvent, "Enable processing same event processing", true);

  void processSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigma0sWithLabel const& /*sigmas*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    tripletTrackTrackSigma0Builder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, sigma0WithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackSigma0, processSameEventMc, "Enable processing same event processing", false);

  void processMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoSigma0s const& /*sigmas*/)
  {
    tripletTrackTrackSigma0Builder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, trackPartition1, trackPartition2, sigma0Partition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackSigma0, processMixedEvent, "Enable processing mixed event processing", true);

  void processMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigma0sWithLabel const& /*sigmas*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    tripletTrackTrackSigma0Builder.processMixedEvent<modes::Mode::kMe_Reco_Mc>(cols, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, sigma0WithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackSigma0, processMixedEventMc, "Enable processing mixed event processing", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoTripletTrackTrackSigma0>(context),
  };
  return workflow;
}
