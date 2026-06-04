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

/// \file femtoTripletTrackTrackCascade.cxx
/// \brief Tasks that computes correlation between two tracks and a cascade
/// \author Roberta Ferioli roberta.ferioli@cern.ch, Raffaele Del Grande raffaele.del.grande@cern.ch, CTU Prague

#include "PWGCF/Femto/Core/cascadeBuilder.h"
#include "PWGCF/Femto/Core/cascadeHistManager.h"
#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/closeTripletRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/tripletBuilder.h"
#include "PWGCF/Femto/Core/tripletHistManager.h"
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

struct FemtoTripletTrackTrackCascade {
  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoXis = o2::soa::Join<o2::aod::FXis, o2::aod::FXiMasks>;
  using FemtoOmegas = o2::soa::Join<o2::aod::FOmegas, o2::aod::FOmegaMasks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;
  using FemtoXisWithLabel = o2::soa::Join<FemtoXis, o2::aod::FLambdaLabels>;
  using FemtoOmegasWithLabel = o2::soa::Join<FemtoOmegas, o2::aod::FK0shortLabels>;

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

  o2::framework::Partition<FemtoTracks> trackPartition1 = MAKE_TRACK_PARTITION(confTrackSelections1);
  o2::framework::Partition<FemtoTracks> trackPartition2 = MAKE_TRACK_PARTITION(confTrackSelections2);
  o2::framework::Preslice<FemtoTracks> perColtracks = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition1 = MAKE_TRACK_PARTITION(confTrackSelections1);
  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition2 = MAKE_TRACK_PARTITION(confTrackSelections2);
  o2::framework::Preslice<FemtoTracksWithLabel> perColtracksWithLabel = o2::aod::femtobase::stored::fColId;

  // setup for daughters
  trackhistmanager::ConfCascadePosDauBinning confPosDauBinning;
  trackhistmanager::ConfCascadeNegDauBinning confNegDauBinning;
  trackhistmanager::ConfCascadeBachelorBinning confBachelorBinning;

  // setup xis
  cascadebuilder::ConfXiSelection confXiSelection;
  cascadehistmanager::ConfXiBinning confXiBinning;
  particlecleaner::ConfXiCleaner1 confXiCleaner;
  particlecleaner::ConfTrackCleaner1 confTrackCleaner1;
  particlecleaner::ConfTrackCleaner2 confTrackCleaner2;

  o2::framework::Partition<FemtoXis> xiPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  o2::framework::Preslice<FemtoXis> perColXis = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoXisWithLabel> xiWithLabelPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  o2::framework::Preslice<FemtoXisWithLabel> perColXisWithLabel = o2::aod::femtobase::stored::fColId;

  // setup omegas
  cascadebuilder::ConfOmegaSelection confOmegaSelection;
  cascadehistmanager::ConfOmegaBinning confOmegaBinning;
  particlecleaner::ConfOmegaCleaner1 confOmegaCleaner;

  o2::framework::Partition<FemtoOmegas> omegaPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  o2::framework::Preslice<FemtoOmegas> perColOmegas = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoOmegasWithLabel> omegaWithLabelPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  o2::framework::Preslice<FemtoOmegasWithLabel> perColOmegasWithLabel = o2::aod::femtobase::stored::fColId;

  // setup triplets
  triplethistmanager::ConfTripletBinning confTripletBinning;
  triplethistmanager::ConfTripletCuts confTripletCuts;

  closetripletrejection::ConfCtrTrackTrackTrack confCtr;

  tripletbuilder::TripletTrackTrackCascadeBuilder<
    modes::Cascade::kXi,
    trackhistmanager::PrefixTrack1,
    trackhistmanager::PrefixTrack2,
    cascadehistmanager::PrefixXi,
    trackhistmanager::PrefixCascadeBachelor,
    trackhistmanager::PrefixCascadePosDaughter,
    trackhistmanager::PrefixCascadeNegDaughter,
    triplethistmanager::PrefixTrackTrackCascadeSe,
    triplethistmanager::PrefixTrackTrackCascadeMe,
    closetripletrejection::PrefixTrack1Track2Se,
    closetripletrejection::PrefixTrack1CascadeBachelorSe,
    closetripletrejection::PrefixTrack2CascadeBachelorSe,
    closetripletrejection::PrefixTrack1V0DaughterSe,
    closetripletrejection::PrefixTrack2V0DaughterSe,
    closetripletrejection::PrefixTrack1Track2Me,
    closetripletrejection::PrefixTrack1CascadeBachelorMe,
    closetripletrejection::PrefixTrack2CascadeBachelorMe,
    closetripletrejection::PrefixTrack1V0DaughterMe,
    closetripletrejection::PrefixTrack2V0DaughterMe>
    tripletTrackTrackXiBuilder;

  tripletbuilder::TripletTrackTrackCascadeBuilder<
    modes::Cascade::kOmega,
    trackhistmanager::PrefixTrack1,
    trackhistmanager::PrefixTrack2,
    cascadehistmanager::PrefixOmega,
    trackhistmanager::PrefixCascadeBachelor,
    trackhistmanager::PrefixCascadePosDaughter,
    trackhistmanager::PrefixCascadeNegDaughter,
    triplethistmanager::PrefixTrackTrackCascadeSe,
    triplethistmanager::PrefixTrackTrackCascadeMe,
    closetripletrejection::PrefixTrack1Track2Se,
    closetripletrejection::PrefixTrack1CascadeBachelorSe,
    closetripletrejection::PrefixTrack2CascadeBachelorSe,
    closetripletrejection::PrefixTrack1V0DaughterSe,
    closetripletrejection::PrefixTrack2V0DaughterSe,
    closetripletrejection::PrefixTrack1Track2Me,
    closetripletrejection::PrefixTrack1CascadeBachelorMe,
    closetripletrejection::PrefixTrack2CascadeBachelorMe,
    closetripletrejection::PrefixTrack1V0DaughterMe,
    closetripletrejection::PrefixTrack2V0DaughterMe>
    tripletTrackTrackOmegaBuilder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  triplethistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoTrackTrackCascade", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  closepairrejection::ConfCprTrackCascadeBachelor confCprBachelor;
  closepairrejection::ConfCprTrackV0Daughter confCprV0Daughter;

  void init(o2::framework::InitContext&)
  {
    if ((doprocessXiSameEvent + doprocessXiSameEventMc) > 1 || (doprocessXiMixedEvent + doprocessXiMixedEventMc) > 1 || (doprocessOmegaSameEvent + doprocessOmegaSameEventMc) > 1 || (doprocessOmegaMixedEvent + doprocessOmegaMixedEventMc) > 1) {
      LOG(fatal) << "More than 1 same or mixed event process function is activated. Breaking...";
    }
    bool processData = doprocessXiSameEvent || doprocessXiMixedEvent || doprocessOmegaSameEvent || doprocessOmegaMixedEvent;
    bool processMc = doprocessXiSameEventMc || doprocessXiMixedEventMc || doprocessOmegaSameEventMc || doprocessOmegaMixedEventMc;
    if (processData && processMc) {
      LOG(fatal) << "Both data and mc processing is activated. Breaking...";
    }

    bool processXi = doprocessXiSameEvent || doprocessXiSameEventMc || doprocessXiMixedEvent || doprocessXiMixedEventMc;
    bool processOmega = doprocessOmegaSameEvent || doprocessOmegaSameEventMc || doprocessOmegaMixedEvent || doprocessOmegaMixedEventMc;

    if (processXi && processOmega) {
      LOG(fatal) << "Both xi-track and omega-track processing is enabled. Breaking...";
    }
    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histogram specs
    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec1;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec2;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> bachelorHistSpec;

    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDauSpec;
    std::map<cascadehistmanager::CascadeHist, std::vector<o2::framework::AxisSpec>> xiHistSpec;
    std::map<cascadehistmanager::CascadeHist, std::vector<o2::framework::AxisSpec>> omegaHistSpec;
    std::map<triplethistmanager::TripletHist, std::vector<o2::framework::AxisSpec>> tripletTrackTrackCascadeHistSpec;

    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecBachelor = closepairrejection::makeCprHistSpecMap(confCprBachelor);
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecV0Daughter = closepairrejection::makeCprHistSpecMap(confCprV0Daughter);
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> ctrHistSpec = closepairrejection::makeCprHistSpecMap(confCtr);

    if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      trackHistSpec1 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning1);
      trackHistSpec2 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning2);
      bachelorHistSpec = trackhistmanager::makeTrackHistSpecMap(confBachelorBinning);
      posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
      negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
      tripletTrackTrackCascadeHistSpec = triplethistmanager::makeTripletHistSpecMap(confTripletBinning);

      if (processXi) {
        xiHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confXiBinning);
        tripletTrackTrackXiBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confCollisionBinning, confTrackSelections1, confTrackSelections2, confTrackCleaner1, confCtr, confXiSelection, confXiCleaner, confCprBachelor, confCprV0Daughter, confMixing, confTripletBinning, confTripletCuts, colHistSpec, trackHistSpec1, trackHistSpec2, xiHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, tripletTrackTrackCascadeHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter, ctrHistSpec);
      } else {
        omegaHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confOmegaBinning);
        tripletTrackTrackOmegaBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confCollisionBinning, confTrackSelections1, confTrackSelections2, confTrackCleaner1, confCtr, confOmegaSelection, confOmegaCleaner, confCprBachelor, confCprV0Daughter, confMixing, confTripletBinning, confTripletCuts, colHistSpec, trackHistSpec1, trackHistSpec2, omegaHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, tripletTrackTrackCascadeHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter, ctrHistSpec);
      }
    } else {
      colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      trackHistSpec1 = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning1);
      trackHistSpec2 = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning2);
      bachelorHistSpec = trackhistmanager::makeTrackMcHistSpecMap(confBachelorBinning);
      posDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confPosDauBinning);
      negDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confNegDauBinning);
      tripletTrackTrackCascadeHistSpec = triplethistmanager::makeTripletMcHistSpecMap(confTripletBinning);
      if (processXi) {
        xiHistSpec = cascadehistmanager::makeCascadeMcHistSpecMap(confXiBinning);
        tripletTrackTrackXiBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, confCollisionBinning, confTrackSelections1, confTrackSelections2, confTrackCleaner1, confCtr, confXiSelection, confXiCleaner, confCprBachelor, confCprV0Daughter, confMixing, confTripletBinning, confTripletCuts, colHistSpec, trackHistSpec1, trackHistSpec2, xiHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, tripletTrackTrackCascadeHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter, ctrHistSpec);
      } else {
        omegaHistSpec = cascadehistmanager::makeCascadeMcHistSpecMap(confOmegaBinning);
        tripletTrackTrackOmegaBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, confCollisionBinning, confTrackSelections1, confTrackSelections2, confTrackCleaner1, confCtr, confOmegaSelection, confOmegaCleaner, confCprBachelor, confCprV0Daughter, confMixing, confTripletBinning, confTripletCuts, colHistSpec, trackHistSpec1, trackHistSpec2, omegaHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, tripletTrackTrackCascadeHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter, ctrHistSpec);
      }
    }
    hRegistry.print();
  };

  void processXiSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoXis const& /*xis*/)
  {
    tripletTrackTrackXiBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition1, trackPartition2, xiPartition, cache);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackCascade, processXiSameEvent, "Enable processing same event processing for tracks and xis", true);

  void processXiSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoXisWithLabel const& /*xis*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    tripletTrackTrackXiBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, xiWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackCascade, processXiSameEventMc, "Enable processing same event processing for tracks and xis with mc information", false);

  void processXiMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoXis const& /*xis*/)
  {
    tripletTrackTrackXiBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition1, trackPartition2, xiPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackCascade, processXiMixedEvent, "Enable processing mixed event processing for tracks and xis", true);

  void processXiMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoXisWithLabel const& /*xis*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& /*mcMothers*/, o2::aod::FMcPartMoths const& /*mcPartonicMothers*/)
  {
    tripletTrackTrackXiBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, xiWithLabelPartition, mcParticles, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackCascade, processXiMixedEventMc, "Enable processing mixed event processing for tracks and xis with mc information", false);

  void processOmegaSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoOmegas const& /*omega*/)
  {
    tripletTrackTrackOmegaBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition1, trackPartition2, omegaPartition, cache);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackCascade, processOmegaSameEvent, "Enable processing same event processing for tracks and omegas", false);

  void processOmegaSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoOmegasWithLabel const& /*omegas*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    tripletTrackTrackOmegaBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, omegaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackCascade, processOmegaSameEventMc, "Enable processing same event processing for tracks and omegas with mc information", false);

  void processOmegaMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoOmegas const& /*omegas*/)
  {
    tripletTrackTrackOmegaBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition1, trackPartition2, omegaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackCascade, processOmegaMixedEvent, "Enable processing mixed event processing for tracks and omegas", false);

  void processOmegaMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoXisWithLabel const& /*xis*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& /*mcMothers*/, o2::aod::FMcPartMoths const& /*mcPartonicMothers*/)
  {
    tripletTrackTrackOmegaBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, omegaWithLabelPartition, mcParticles, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackCascade, processOmegaMixedEventMc, "Enable processing mixed event processing for tracks and omegas with mc information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoTripletTrackTrackCascade>(cfgc),
  };
  return workflow;
}
