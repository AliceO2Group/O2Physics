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

/// \file femtoDndetaTriplet.cxx
/// \brief charged-particle pseudorapidity density in events triggered on femto triplets
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
///
/// An event is triggered if at least one triplet passes the particle selections, the triplet cleaner,
/// the close triplet rejection and the triplet cuts (e.g. TripletCuts.q3Max).
/// There is no generator-level triplet builder yet, so the event/signal loss is not conditioned on a generator-level trigger.
/// Needs a femto producer run in pass-through mode.

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/closeTripletRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/dndetaBuilder.h"
#include "PWGCF/Femto/Core/dndetaHistManager.h"
#include "PWGCF/Femto/Core/mcBuilder.h"
#include "PWGCF/Femto/Core/mcParticleHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairCleaner.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/tripletBuilder.h"
#include "PWGCF/Femto/Core/tripletHistManager.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGCF/Femto/Core/v0HistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <map>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoDndetaTriplet {

  // setup tables
  // FCols provides the default centrality (FT0M in pp, FT0C in PbPb), FColMults the INEL>0 and correlation estimators
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks, o2::aod::FColMults>;
  using FemtoCollision = FemtoCollisions::iterator;
  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;

  // same join as femtoTrackQa, the track histogram managers (QA) need mass, DCAs, extras and PID
  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks, o2::aod::FTrackMass, o2::aod::FTrackDcas, o2::aod::FTrackExtras, o2::aod::FTrackPids>;
  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;

  using FemtoLambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks>;
  using FemtoLambdasWithLabel = o2::soa::Join<FemtoLambdas, o2::aod::FLambdaLabels>;
  using FemtoK0shorts = o2::soa::Join<o2::aod::FK0shorts, o2::aod::FK0shortMasks>;

  using FemtoMcParticles = o2::soa::Join<o2::aod::FMcParticles, o2::aod::FMcMotherLabels>;

  o2::framework::SliceCache cache;
  o2::framework::Service<o2::framework::O2DatabasePDG> pdg = {};

  // setup dndeta
  dndetabuilder::ConfDndeta confDndeta;
  dndetabuilder::ConfDndetaAcceptance confDndetaAcceptance;
  dndetahistmanager::ConfDndetaBinning confDndetaBinning;
  dndetabuilder::DndetaBuilder dndetaBuilder;

  // pass-through mode does not guarantee fMcColId is written in sorted order (reco-driven label
  // rows can be created out of mc-collision order), so this must be an unsorted grouping
  o2::framework::PresliceUnsorted<FemtoMcParticles> perMcColParticles = o2::aod::femtomcparticle::fMcColId;

  // setup collisions (applied in the dndeta builder, no filter, so the cutflow sees every collision)
  collisionbuilder::ConfCollisionSelection collisionSelection;
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup tracks for dndeta (bitmask with track quality bits only, see dndetaBuilder.h)
  trackbuilder::ConfTrackSelectionDndetaGlobal confDndetaTracksGlobal;
  trackhistmanager::ConfDndetaTrackGlobalBinning confDndetaTrackGlobalBinning;
  trackhistmanager::ConfDndetaTrackGlobalQaBinning confDndetaTrackGlobalQaBinning;

  trackbuilder::ConfTrackSelectionDndetaItsOnly confDndetaTracksItsOnly;
  trackhistmanager::ConfDndetaTrackItsOnlyBinning confDndetaTrackItsOnlyBinning;
  trackhistmanager::ConfDndetaTrackItsOnlyQaBinning confDndetaTrackItsOnlyQaBinning;

  o2::framework::Partition<FemtoTracks> dndetaTrackPartitionGlobal = MAKE_TRACK_PARTITION(confDndetaTracksGlobal);
  o2::framework::Partition<FemtoTracks> dndetaTrackPartitionItsOnly = MAKE_TRACK_PARTITION(confDndetaTracksItsOnly);

  o2::framework::Partition<FemtoTracksWithLabel> dndetaTrackWithLabelPartitionGlobal = MAKE_TRACK_PARTITION(confDndetaTracksGlobal);
  o2::framework::Partition<FemtoTracksWithLabel> dndetaTrackWithLabelPartitionItsOnly = MAKE_TRACK_PARTITION(confDndetaTracksItsOnly);

  // the *Mc process functions take the ungrouped tables (no ::iterator first argument), so DPL never
  // auto-registers the fColId slicing cache for FemtoTracksWithLabel; this unused Preslice is the only
  // thing that requests/enables it, otherwise sliceByCached() throws "Disabled cache ... is requested"
  // once the first accepted reco collision actually reaches it
  o2::framework::Preslice<FemtoTracksWithLabel> perColTracksWithLabel = o2::aod::femtobase::stored::fColId;

  // setup tracks (trigger)
  trackbuilder::ConfTrackSelection1 confTrackSelection1;
  trackhistmanager::ConfTrackBinning1 confTrackBinning1;
  particlecleaner::ConfTrackCleaner1 confTrackCleaner1;

  trackbuilder::ConfTrackSelection2 confTrackSelection2;
  trackhistmanager::ConfTrackBinning2 confTrackBinning2;
  particlecleaner::ConfTrackCleaner2 confTrackCleaner2;

  trackbuilder::ConfTrackSelection3 confTrackSelection3;
  trackhistmanager::ConfTrackBinning3 confTrackBinning3;
  particlecleaner::ConfTrackCleaner3 confTrackCleaner3;

  o2::framework::Partition<FemtoTracks> trackPartition1 = MAKE_TRACK_PARTITION(confTrackSelection1);
  o2::framework::Partition<FemtoTracks> trackPartition2 = MAKE_TRACK_PARTITION(confTrackSelection2);
  o2::framework::Partition<FemtoTracks> trackPartition3 = MAKE_TRACK_PARTITION(confTrackSelection3);

  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition1 = MAKE_TRACK_PARTITION(confTrackSelection1);
  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition2 = MAKE_TRACK_PARTITION(confTrackSelection2);
  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition3 = MAKE_TRACK_PARTITION(confTrackSelection3);

  // setup for daughters
  trackhistmanager::ConfV01PosDauBinning confPosDauBinning;
  trackhistmanager::ConfV01NegDauBinning confNegDauBinning;

  // setup lambdas (trigger)
  v0builder::ConfLambdaSelection1 confLambdaSelection;
  v0histmanager::ConfLambdaBinning1 confLambdaBinning;
  particlecleaner::ConfLambdaCleaner1 confLambdaCleaner;

  o2::framework::Partition<FemtoLambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  o2::framework::Partition<FemtoLambdasWithLabel> lambdaWithLabelPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  // same reason as perColTracksWithLabel above, needed by processTrackTrackV0Mc
  o2::framework::Preslice<FemtoLambdasWithLabel> perColLambdasWithLabel = o2::aod::femtobase::stored::fColId;

  // setup strangeness yields (dedicated selections, independent of the trigger)
  v0builder::ConfLambdaSelectionStrangeness confStrangeLambdaSelection;
  v0builder::ConfK0shortSelectionStrangeness confStrangeK0shortSelection;

  o2::framework::Partition<FemtoLambdas> strangeLambdaPartition = MAKE_LAMBDA_PARTITION(confStrangeLambdaSelection);
  o2::framework::Partition<FemtoK0shorts> strangeK0shortPartition = MAKE_K0SHORT_PARTITION(confStrangeK0shortSelection);

  // setup mc particles (generator-level trigger)
  // FMcParticles also contains rows of reconstructed secondaries (created via labels), so set
  // McParticleSelection1/2/3.requireOrigin = true (origin = 2, physical primary) for an unbiased generator-level trigger
  mcbuilder::ConfMcParticleSelection1 confMcParticleSelection1;
  mcparticlehistmanager::ConfMcParticleBinning1 confMcParticleBinning1;
  particlecleaner::ConfMcParticleCleaner1 confMcParticleCleaner1;

  mcbuilder::ConfMcParticleSelection2 confMcParticleSelection2;
  mcparticlehistmanager::ConfMcParticleBinning2 confMcParticleBinning2;
  particlecleaner::ConfMcParticleCleaner2 confMcParticleCleaner2;

  mcbuilder::ConfMcParticleSelection3 confMcParticleSelection3;
  mcparticlehistmanager::ConfMcParticleBinning3 confMcParticleBinning3;
  particlecleaner::ConfMcParticleCleaner3 confMcParticleCleaner3;

  o2::framework::Partition<FemtoMcParticles> mcParticlePartition1 = MAKE_MC_PARTICLE_PARTITION(confMcParticleSelection1);
  o2::framework::Partition<FemtoMcParticles> mcParticlePartition2 = MAKE_MC_PARTICLE_PARTITION(confMcParticleSelection2);
  o2::framework::Partition<FemtoMcParticles> mcParticlePartition3 = MAKE_MC_PARTICLE_PARTITION(confMcParticleSelection3);

  // setup triplets
  triplethistmanager::ConfTripletBinning confTripletBinning;
  triplethistmanager::ConfTripletCuts confTripletCuts;
  paircleaner::ConfPairCleanerBinning confTripletCleanerBinning;
  triplethistmanager::ConfMixing confMixing; // only same event is processed, used for particle12(3)AreSameSpecies/seed

  closetripletrejection::ConfCtrTrackTrackTrack confCtrTrackTrackTrack;
  closetripletrejection::ConfCtrTrackTrackV0 confCtrTrackTrackV0;
  closetripletrejection::ConfCtrMcParticleMcParticleMcParticle confCtrMcParticleMcParticleMcParticle; // generator-level trigger, no cut by default

  tripletbuilder::TripletTrackTrackTrackBuilder<
    trackhistmanager::PrefixTrack1,
    trackhistmanager::PrefixTrack2,
    trackhistmanager::PrefixTrack3,
    triplethistmanager::PrefixTrackTrackTrackSe,
    triplethistmanager::PrefixTrackTrackTrackMe,
    closetripletrejection::PrefixTrack1Track2Se,
    closetripletrejection::PrefixTrack2Track3Se,
    closetripletrejection::PrefixTrack1Track3Se,
    closetripletrejection::PrefixTrack1Track2Me,
    closetripletrejection::PrefixTrack2Track3Me,
    closetripletrejection::PrefixTrack1Track3Me>
    tripletTrackTrackTrackBuilder;

  tripletbuilder::TripletTrackTrackV0Builder<
    modes::V0::kLambda,
    trackhistmanager::PrefixTrack1,
    trackhistmanager::PrefixTrack2,
    v0histmanager::PrefixLambda1,
    trackhistmanager::PrefixV01PosDaughter,
    trackhistmanager::PrefixV01NegDaughter,
    triplethistmanager::PrefixTrackTrackLambdaSe,
    triplethistmanager::PrefixTrackTrackLambdaMe,
    closetripletrejection::PrefixTrack1Track2Se,
    closetripletrejection::PrefixTrack1V0Se,
    closetripletrejection::PrefixTrack2V0Se,
    closetripletrejection::PrefixTrack1Track2Me,
    closetripletrejection::PrefixTrack1V0Me,
    closetripletrejection::PrefixTrack2V0Me>
    tripletTrackTrackLambdaBuilder;

  tripletbuilder::TripletMcParticleMcParticleMcParticleBuilder<
    mcparticlehistmanager::PrefixMcParticle1,
    mcparticlehistmanager::PrefixMcParticle2,
    mcparticlehistmanager::PrefixMcParticle3,
    triplethistmanager::PrefixMcParticleMcParticleMcParticleSe,
    triplethistmanager::PrefixMcParticleMcParticleMcParticleMe,
    closetripletrejection::PrefixMcParticle1McParticle2Se,
    closetripletrejection::PrefixMcParticle2McParticle3Se,
    closetripletrejection::PrefixMcParticle1McParticle3Se,
    closetripletrejection::PrefixMcParticle1McParticle2Me,
    closetripletrejection::PrefixMcParticle2McParticle3Me,
    closetripletrejection::PrefixMcParticle1McParticle3Me>
    tripletMcParticleMcParticleMcParticleBuilder;

  // dndeta histograms and trigger QA are kept in separate registries, so the fixed histogram names of the builders cannot clash
  o2::framework::HistogramRegistry hRegistry{"FemtoDndetaTriplet", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};
  o2::framework::HistogramRegistry hRegistryTrigger{"FemtoDndetaTripletTrigger", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};
  o2::framework::HistogramRegistry hRegistryGenTrigger{"FemtoDndetaTripletGenTrigger", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    const bool processDataWithStrangeness = doprocessTrackTrackTrackWithStrangeness || doprocessTrackTrackV0WithStrangeness;
    const bool processMc = doprocessTrackTrackTrackMc || doprocessTrackTrackV0Mc;
    const int nProcess = static_cast<int>(doprocessTrackTrackTrack) + static_cast<int>(doprocessTrackTrackV0) +
                         static_cast<int>(doprocessTrackTrackTrackWithStrangeness) + static_cast<int>(doprocessTrackTrackV0WithStrangeness) +
                         static_cast<int>(doprocessTrackTrackTrackMc) + static_cast<int>(doprocessTrackTrackV0Mc);
    if (nProcess != 1) {
      LOG(fatal) << "Exactly one process function has to be activated (found " << nProcess << "). Breaking...";
    }

    const modes::Trigger trigger = (doprocessTrackTrackTrack || doprocessTrackTrackTrackWithStrangeness || doprocessTrackTrackTrackMc) ? modes::Trigger::kTrackTrackTrack : modes::Trigger::kTrackTrackV0;
    LOG(info) << "Trigger: " << modes::triggerToString(trigger);

    dndetaBuilder.init(&hRegistry, confDndeta, confDndetaAcceptance, collisionSelection,
                       confDndetaTracksGlobal, confDndetaTrackGlobalBinning, confDndetaTrackGlobalQaBinning,
                       confDndetaTracksItsOnly, confDndetaTrackItsOnlyBinning, confDndetaTrackItsOnlyQaBinning,
                       confDndetaBinning, processMc, processDataWithStrangeness);
    if (dndetaBuilder.useGenTrigger() && !processMc) {
      LOG(warn) << "Dndeta.useGenTrigger is only used for mc processing. Ignoring it.";
    }

    // setup histogram specs
    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec1;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec2;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec3;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDauSpec;
    std::map<v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> lambdaHistSpec;
    std::map<triplethistmanager::TripletHist, std::vector<o2::framework::AxisSpec>> tripletHistSpec;
    std::map<paircleaner::PairCleanerHist, std::vector<o2::framework::AxisSpec>> tripletCleanerHistSpec = paircleaner::makePairCleanerHistSpecMap(confTripletCleanerBinning);

    if (!processMc) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      trackHistSpec1 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning1);
      trackHistSpec2 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning2);
      trackHistSpec3 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning3);
      posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
      negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
      lambdaHistSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
      tripletHistSpec = triplethistmanager::makeTripletHistSpecMap(confTripletBinning, confMixing);
    } else {
      colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      trackHistSpec1 = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning1);
      trackHistSpec2 = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning2);
      trackHistSpec3 = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning3);
      posDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confPosDauBinning);
      negDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confNegDauBinning);
      lambdaHistSpec = v0histmanager::makeV0McHistSpecMap(confLambdaBinning);
      tripletHistSpec = triplethistmanager::makeTripletMcHistSpecMap(confTripletBinning, confMixing);
    }

    // only the builder of the active trigger is initialized
    // only same event processing is used; the mixed event histograms are booked but stay empty
    if (trigger == modes::Trigger::kTrackTrackTrack) {
      std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> ctrHistSpec = closepairrejection::makeCprHistSpecMap(confCtrTrackTrackTrack);
      if (!processMc) {
        tripletTrackTrackTrackBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistryTrigger, confCollisionBinning, confTrackSelection1, confTrackSelection2, confTrackSelection3, confTrackCleaner1, confTrackCleaner2, confTrackCleaner3, confCtrTrackTrackTrack, confMixing, confTripletBinning, confTripletCuts, colHistSpec, trackHistSpec1, trackHistSpec2, trackHistSpec3, tripletHistSpec, ctrHistSpec, tripletCleanerHistSpec);
      } else {
        tripletTrackTrackTrackBuilder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistryTrigger, confCollisionBinning, confTrackSelection1, confTrackSelection2, confTrackSelection3, confTrackCleaner1, confTrackCleaner2, confTrackCleaner3, confCtrTrackTrackTrack, confMixing, confTripletBinning, confTripletCuts, colHistSpec, trackHistSpec1, trackHistSpec2, trackHistSpec3, tripletHistSpec, ctrHistSpec, tripletCleanerHistSpec);
      }
    } else {
      std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> ctrHistSpec = closepairrejection::makeCprHistSpecMap(confCtrTrackTrackV0);
      if (!processMc) {
        tripletTrackTrackLambdaBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistryTrigger, confCollisionBinning, confTrackSelection1, confTrackSelection2, confLambdaSelection, confTrackCleaner1, confTrackCleaner2, confLambdaCleaner, confCtrTrackTrackV0, confMixing, confTripletBinning, confTripletCuts, colHistSpec, trackHistSpec1, trackHistSpec2, lambdaHistSpec, posDauSpec, negDauSpec, tripletHistSpec, ctrHistSpec, tripletCleanerHistSpec);
      } else {
        tripletTrackTrackLambdaBuilder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistryTrigger, confCollisionBinning, confTrackSelection1, confTrackSelection2, confLambdaSelection, confTrackCleaner1, confTrackCleaner2, confLambdaCleaner, confCtrTrackTrackV0, confMixing, confTripletBinning, confTripletCuts, colHistSpec, trackHistSpec1, trackHistSpec2, lambdaHistSpec, posDauSpec, negDauSpec, tripletHistSpec, ctrHistSpec, tripletCleanerHistSpec);
      }
    }

    // generator-level trigger; the species come from McParticleSelection1/2/3, so it covers
    // both triplet triggers (a lambda is just another mc particle)
    if (processMc && dndetaBuilder.useGenTrigger()) {
      std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> genColHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      std::map<mcparticlehistmanager::McParticleHist, std::vector<o2::framework::AxisSpec>> mcParticleHistSpec1 = mcparticlehistmanager::makeMcParticleHistSpecMap(confMcParticleBinning1);
      std::map<mcparticlehistmanager::McParticleHist, std::vector<o2::framework::AxisSpec>> mcParticleHistSpec2 = mcparticlehistmanager::makeMcParticleHistSpecMap(confMcParticleBinning2);
      std::map<mcparticlehistmanager::McParticleHist, std::vector<o2::framework::AxisSpec>> mcParticleHistSpec3 = mcparticlehistmanager::makeMcParticleHistSpecMap(confMcParticleBinning3);
      std::map<triplethistmanager::TripletHist, std::vector<o2::framework::AxisSpec>> genTripletHistSpec = triplethistmanager::makeTripletMcTruthHistSpecMap(confTripletBinning, confMixing);
      std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> genCtrHistSpec = closepairrejection::makeCprHistSpecMap(confCtrMcParticleMcParticleMcParticle);
      tripletMcParticleMcParticleMcParticleBuilder.init<modes::Mode::kSe_Mc, modes::Mode::kMe_Mc>(&hRegistryGenTrigger, confCollisionBinning, confMcParticleSelection1, confMcParticleSelection2, confMcParticleSelection3, confMcParticleBinning1, confMcParticleBinning2, confMcParticleBinning3, confMcParticleCleaner1, confMcParticleCleaner2, confMcParticleCleaner3, confCtrMcParticleMcParticleMcParticle, confMixing, confTripletBinning, confTripletCuts, genColHistSpec, mcParticleHistSpec1, mcParticleHistSpec2, mcParticleHistSpec3, genTripletHistSpec, genCtrHistSpec, tripletCleanerHistSpec);
    }

    hRegistry.print();
    hRegistryTrigger.print();
    hRegistryGenTrigger.print();
  };

  /// generator-level triplet trigger, or always true if not configured
  template <typename T1, typename T2, typename T3, typename T4>
  bool evaluateGenTrigger(T1 const& mcCol, T2 const& mcParticles, T3 const& mcMothers, T4 const& mcPartonicMothers)
  {
    if (!dndetaBuilder.useGenTrigger()) {
      return true;
    }
    return tripletMcParticleMcParticleMcParticleBuilder.processSameEvent<modes::Mode::kSe_Mc>(mcCol, mcParticles, mcMothers, mcPartonicMothers, mcParticlePartition1, mcParticlePartition2, mcParticlePartition3, cache);
  }

  template <typename T>
  void fillStrangeness(T const& col)
  {
    auto lambdas = strangeLambdaPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto k0shorts = strangeK0shortPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    dndetaBuilder.processStrangeness(col, lambdas, k0shorts);
  }

  // data
  void processTrackTrackTrack(FemtoCollision const& col, FemtoTracks const& tracks)
  {
    dndetaBuilder.processData(col, dndetaTrackPartitionGlobal, dndetaTrackPartitionItsOnly, cache, [&]() {
      return tripletTrackTrackTrackBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition1, trackPartition2, trackPartition3, cache);
    });
  }
  PROCESS_SWITCH(FemtoDndetaTriplet, processTrackTrackTrack, "Data, track-track-track trigger", true);

  void processTrackTrackV0(FemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/)
  {
    dndetaBuilder.processData(col, dndetaTrackPartitionGlobal, dndetaTrackPartitionItsOnly, cache, [&]() {
      return tripletTrackTrackLambdaBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition1, trackPartition2, lambdaPartition, cache);
    });
  }
  PROCESS_SWITCH(FemtoDndetaTriplet, processTrackTrackV0, "Data, track-track-lambda trigger", false);

  // data with strangeness yields
  void processTrackTrackTrackWithStrangeness(FemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoK0shorts const& /*k0shorts*/)
  {
    if (dndetaBuilder.processData(col, dndetaTrackPartitionGlobal, dndetaTrackPartitionItsOnly, cache, [&]() {
          return tripletTrackTrackTrackBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition1, trackPartition2, trackPartition3, cache);
        })) {
      fillStrangeness(col);
    }
  }
  PROCESS_SWITCH(FemtoDndetaTriplet, processTrackTrackTrackWithStrangeness, "Data with strangeness yields, track-track-track trigger", false);

  void processTrackTrackV0WithStrangeness(FemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoK0shorts const& /*k0shorts*/)
  {
    if (dndetaBuilder.processData(col, dndetaTrackPartitionGlobal, dndetaTrackPartitionItsOnly, cache, [&]() {
          return tripletTrackTrackLambdaBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition1, trackPartition2, lambdaPartition, cache);
        })) {
      fillStrangeness(col);
    }
  }
  PROCESS_SWITCH(FemtoDndetaTriplet, processTrackTrackV0WithStrangeness, "Data with strangeness yields, track-track-lambda trigger", false);

  // mc
  void processTrackTrackTrackMc(FemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    dndetaBuilder.processMc(
      cols, mcCols, dndetaTrackWithLabelPartitionGlobal, dndetaTrackWithLabelPartitionItsOnly, cache, mcParticles, mcMothers, perMcColParticles, pdg,
      [&](auto const& col) {
        return tripletTrackTrackTrackBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, trackWithLabelPartition3, mcParticles, mcMothers, mcPartonicMothers, cache);
      },
      [&](auto const& mcCol) { return evaluateGenTrigger(mcCol, mcParticles, mcMothers, mcPartonicMothers); });
  }
  PROCESS_SWITCH(FemtoDndetaTriplet, processTrackTrackTrackMc, "MC, track-track-track trigger", false);

  void processTrackTrackV0Mc(FemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, FemtoMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    dndetaBuilder.processMc(
      cols, mcCols, dndetaTrackWithLabelPartitionGlobal, dndetaTrackWithLabelPartitionItsOnly, cache, mcParticles, mcMothers, perMcColParticles, pdg,
      [&](auto const& col) {
        return tripletTrackTrackLambdaBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, lambdaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
      },
      [&](auto const& mcCol) { return evaluateGenTrigger(mcCol, mcParticles, mcMothers, mcPartonicMothers); });
  }
  PROCESS_SWITCH(FemtoDndetaTriplet, processTrackTrackV0Mc, "MC, track-track-lambda trigger", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoDndetaTriplet>(context),
  };
  return workflow;
}
