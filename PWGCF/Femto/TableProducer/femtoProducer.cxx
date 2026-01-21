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

/// \file femtoProducer.cxx
/// \brief Tasks that produces all femto tables
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/Femto/Core/cascadeBuilder.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/kinkBuilder.h"
#include "PWGCF/Femto/Core/mcBuilder.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGLF/DataModel/LFKinkDecayTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/OutputObjHeader.h"
#include "Framework/runDataProcessing.h"

#include "fairlogger/Logger.h"

#include <chrono>
#include <cstdint>

using namespace o2::analysis::femto;

namespace o2::analysis::femto
{
namespace consumeddata
{
using Run3PpCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::Mults, o2::aod::CentFT0As, o2::aod::CentFT0Cs, o2::aod::CentFT0Ms>;
using Run3PpMcRecoCollisions = o2::soa::Join<Run3PpCollisions, o2::aod::McCollisionLabels>;
using Run3PpMcGenCollisions = o2::soa::Join<o2::aod::McCollisions, o2::aod::MultsExtraMC, o2::aod::McCentFT0Ms>;

using Run3FullPidTracks =
  soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA,
            o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr, o2::aod::pidTPCFullDe, o2::aod::pidTPCFullTr, o2::aod::pidTPCFullHe,
            o2::aod::pidTOFFullEl, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr, o2::aod::pidTOFFullDe, o2::aod::pidTOFFullTr, o2::aod::pidTOFFullHe,
            o2::aod::pidTOFbeta, o2::aod::pidTOFmass>;
using Run3McRecoTracks = soa::Join<Run3FullPidTracks, o2::aod::McTrackLabels>;

using Run3Vzeros = o2::aod::V0Datas;
using Run3RecoVzeros = o2::soa::Join<o2::aod::V0Datas, o2::aod::McV0Labels>;

using Run3Cascades = o2::aod::CascDatas;
using Run3RecoCascades = o2::soa::Join<Run3Cascades, o2::aod::McCascLabels>;

using Run3Kinks = o2::aod::KinkCands;

using Run3McGenParticles = o2::aod::McParticles;

} // namespace consumeddata
} // namespace o2::analysis::femto

struct FemtoProducer {

  // ccdb
  collisionbuilder::ConfCcdb confCcdb;

  // collision builder
  collisionbuilder::CollisionBuilderProducts collisionBuilderProducts;
  collisionbuilder::ConfCollisionTables confCollisionTables;
  collisionbuilder::ConfCollisionFilters confCollisionFilters;
  collisionbuilder::ConfCollisionBits confCollisionBits;
  collisionbuilder::ConfCollisionRctFlags confCollisionRctFlags;
  collisionbuilder::CollisionBuilder<collisionbuilder::ColSelHistName> collisionBuilder;

  // track builder
  trackbuilder::TrackBuilderProducts trackBuilderProducts;
  trackbuilder::ConfTrackTables confTrackTables;
  trackbuilder::TrackBuilder<trackbuilder::TrackSelHistName> trackBuilder;
  trackbuilder::ConfTrackBits confTrackBits;
  trackbuilder::ConfTrackFilters confTrackFilters;

  // v0 builders
  v0builder::V0BuilderProducts v0builderProducts;
  v0builder::ConfV0Tables confV0Tables;
  v0builder::ConfV0Filters confV0Filters;
  v0builder::ConfK0shortBits confK0shortBits;
  v0builder::V0Builder<modes::V0::kK0short, v0builder::K0shortSelHistName> k0shortBuilder;
  v0builder::ConfLambdaBits confLambdaBits;
  v0builder::V0Builder<modes::V0::kLambda, v0builder::LambdaSelHistName> lambdaBuilder;
  v0builder::V0Builder<modes::V0::kAntiLambda, v0builder::AntilambdaSelHistName> antilambdaBuilder;

  // cascade builder
  cascadebuilder::CascadeBuilderProducts cascadeBuilderProducts;
  cascadebuilder::ConfCascadeTables confCascadeTables;
  cascadebuilder::ConfCascadeFilters confCascadeFilters;
  cascadebuilder::ConfXiBits confXiBits;
  cascadebuilder::CascadeBuilder<modes::Cascade::kXi, cascadebuilder::XiSelHistName> xiBuilder;
  cascadebuilder::ConfOmegaBits confOmegaBits;
  cascadebuilder::CascadeBuilder<modes::Cascade::kOmega, cascadebuilder::OmegaSelHistName> omegaBuilder;

  // kink builder
  kinkbuilder::KinkBuilderProducts kinkBuilderProducts;
  kinkbuilder::ConfKinkTables confKinkTables;
  kinkbuilder::ConfKinkFilters confKinkFilters;
  kinkbuilder::ConfSigmaBits confSigmaBits;
  kinkbuilder::KinkBuilder<modes::Kink::kSigma, kinkbuilder::SigmaSelHistName> sigmaBuilder;
  kinkbuilder::ConfSigmaPlusBits confSigmaPlusBits;
  kinkbuilder::KinkBuilder<modes::Kink::kSigmaPlus, kinkbuilder::SigmaPlusSelHistName> sigmaPlusBuilder;

  // mc builder
  mcbuilder::ConfMc confMc;
  mcbuilder::ConfMcTables confMcTables;
  mcbuilder::McBuilderProducts mcProducts;
  mcbuilder::McBuilder mcBuilder;

  // histogramming
  // add histograms in next iteration
  o2::framework::HistogramRegistry hRegistry{"FemtoProducer", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // data members
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(o2::framework::InitContext& context)
  {
    if ((xiBuilder.fillAnyTable() || omegaBuilder.fillAnyTable()) && (!doprocessTracksV0sCascadesRun3pp && !doprocessTracksV0sCascadesKinksRun3pp)) {
      LOG(fatal) << "At least one cascade table is enabled, but wrong process function is enabled. Breaking...";
    }
    if ((lambdaBuilder.fillAnyTable() || antilambdaBuilder.fillAnyTable() || k0shortBuilder.fillAnyTable()) && (!doprocessTracksV0sCascadesRun3pp && !doprocessTracksV0sRun3pp && !doprocessTracksV0sCascadesKinksRun3pp && !doprocessTracksV0sRun3ppMc)) {
      LOG(fatal) << "At least one v0 table is enabled, but wrong process function is enabled. Breaking...";
    }
    if ((sigmaBuilder.fillAnyTable() || sigmaPlusBuilder.fillAnyTable()) && (!doprocessTracksKinksRun3pp && !doprocessTracksV0sCascadesKinksRun3pp && !doprocessTracksKinksRun3ppMc)) {
      LOG(fatal) << "At least one kink table is enabled, but wrong process function is enabled. Breaking...";
    }
    if (mcBuilder.fillAnyTable() && (!doprocessTracksV0sRun3ppMc && !doprocessTracksKinksRun3ppMc)) {
      LOG(fatal) << "At least one mc table is enabled, but wrong process function is enabled. Breaking...";
    }

    // init ccdb
    ccdb->setURL(confCcdb.ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    // collision selection
    collisionBuilder.init(&hRegistry, confCollisionFilters, confCollisionBits, confCollisionRctFlags, confCcdb, confCollisionTables, context);

    // configure track builder
    trackBuilder.init(&hRegistry, confTrackBits, confTrackFilters, confTrackTables, context);

    // configure v0 builder
    k0shortBuilder.init(&hRegistry, confK0shortBits, confV0Filters, confV0Tables, context);
    lambdaBuilder.init(&hRegistry, confLambdaBits, confV0Filters, confV0Tables, context);
    antilambdaBuilder.init(&hRegistry, confLambdaBits, confV0Filters, confV0Tables, context);

    // configure kink builder
    sigmaBuilder.init(&hRegistry, confSigmaBits, confKinkFilters, confKinkTables, context);
    sigmaPlusBuilder.init(&hRegistry, confSigmaPlusBits, confKinkFilters, confKinkTables, context);

    // cascade selections
    xiBuilder.init(&hRegistry, confXiBits, confCascadeFilters, confCascadeTables, context);
    omegaBuilder.init(&hRegistry, confOmegaBits, confCascadeFilters, confCascadeTables, context);

    // configure mcBuilder
    mcBuilder.init(confMc, confMcTables, context);

    hRegistry.print();
  }

  // processing collisions
  template <modes::System system, typename T1, typename T2, typename T3>
  bool processCollisions(T1 const& col, T2 const& /* bcs*/, T3 const& tracks)
  {
    collisionBuilder.reset();
    auto bc = col.template bc_as<T2>();
    collisionBuilder.initCollision<system>(bc, col, tracks, ccdb, hRegistry);
    if (!collisionBuilder.checkCollision(col)) {
      return false;
    }
    // do not fill collsions into the table here
    // collisions are filled if at least one partilce is found in the collisions
    return true;
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  bool processMcCollisions(T1 const& col, T2 const& mcCols, T3 const& /* bcs*/, T4 const& tracks, T5 const& mcParticles)
  {
    collisionBuilder.reset();
    mcBuilder.reset(mcCols, mcParticles); // we call this function always first so we reset the mcBuilder here
    auto bc = col.template bc_as<T3>();
    collisionBuilder.initCollision<system>(bc, col, tracks, ccdb, hRegistry);
    if (!collisionBuilder.checkCollision(col, mcCols)) {
      return false;
    }
    // do not fill collsions into the table here
    // collisions are filled if at least one partilce is found in the collisions
    return true;
  }

  // processing tracks
  template <modes::System system, typename T1, typename T2>
  void processTracks(T1 const& col, T2 const& tracksWithItsPid)
  {
    trackBuilder.reset(tracksWithItsPid);
    trackBuilder.fillTracks<system>(col, collisionBuilder, collisionBuilderProducts, tracksWithItsPid, trackBuilderProducts);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processMcTracks(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& tracksWithItsPid, T5 const& mcParticles)
  {
    trackBuilder.reset(tracksWithItsPid);
    trackBuilder.fillMcTracks<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, tracks, tracksWithItsPid, trackBuilderProducts, mcParticles, mcBuilder, mcProducts);
  }

  // processing v0s
  template <modes::System system, typename T1, typename T2, typename T3, typename T4>
  void processV0s(T1 const& col, T2 const& tracks, T3 const& tracksWithItsPid, T4 const& v0s)
  {
    lambdaBuilder.fillV0s<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, v0builderProducts, v0s, tracks, tracksWithItsPid, trackBuilder);
    antilambdaBuilder.fillV0s<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, v0builderProducts, v0s, tracks, tracksWithItsPid, trackBuilder);
    k0shortBuilder.fillV0s<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, v0builderProducts, v0s, tracks, tracksWithItsPid, trackBuilder);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processMcV0s(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& tracksWithItsPid, T5 const& v0s, T6 const& mcParticles)
  {
    lambdaBuilder.fillMcV0s<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, v0builderProducts, v0s, tracks, tracksWithItsPid, trackBuilder, mcParticles, mcBuilder, mcProducts);
    antilambdaBuilder.fillMcV0s<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, v0builderProducts, v0s, tracks, tracksWithItsPid, trackBuilder, mcParticles, mcBuilder, mcProducts);
    k0shortBuilder.fillMcV0s<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, v0builderProducts, v0s, tracks, tracksWithItsPid, trackBuilder, mcParticles, mcBuilder, mcProducts);
  }

  // processing kinks
  template <modes::System system, typename T1, typename T2, typename T3, typename T4>
  void processKinks(T1 const& col, T2 const& tracks, T3 const& tracksWithItsPid, T4 const& kinks)
  {
    sigmaBuilder.fillKinks<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, tracksWithItsPid, trackBuilder);
    sigmaPlusBuilder.fillKinks<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, tracksWithItsPid, trackBuilder);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processMcKinks(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& tracksWithItsPid, T5 const& kinks, T6 const& mcParticles)
  {
    sigmaBuilder.fillMcKinks<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, tracksWithItsPid, trackBuilder, mcParticles, mcBuilder, mcProducts);
    sigmaPlusBuilder.fillMcKinks<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, tracksWithItsPid, trackBuilder, mcParticles, mcBuilder, mcProducts);
  }

  // processing cascades
  template <modes::System system, typename T1, typename T2, typename T3, typename T4>
  void processCascades(T1 const& col, T2 const& tracks, T3 const& tracksWithItsPid, T4 const& cascades)
  {
    xiBuilder.fillCascades<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, cascadeBuilderProducts,
                                   cascades, tracks, tracksWithItsPid, trackBuilder);
    omegaBuilder.fillCascades<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, cascadeBuilderProducts,
                                      cascades, tracks, tracksWithItsPid, trackBuilder);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processMcCascades(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& tracksWithItsPid, T5 const& cascades, T6 const& mcParticles)
  {
    xiBuilder.fillMcCascades<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, cascadeBuilderProducts, cascades, tracks, tracksWithItsPid, trackBuilder, mcParticles, mcBuilder, mcProducts);
    omegaBuilder.fillMcCascades<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, cascadeBuilderProducts, cascades, tracks, tracksWithItsPid, trackBuilder, mcParticles, mcBuilder, mcProducts);
  }

  // proccess functions
  void processTracksRun3pp(consumeddata::Run3PpCollisions::iterator const& col,
                           o2::aod::BCsWithTimestamps const& bcs,
                           consumeddata::Run3FullPidTracks const& tracks)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksRun3pp, "Process tracks", true);

  // process tracks and v0s
  void processTracksV0sRun3pp(consumeddata::Run3PpCollisions::iterator const& col,
                              o2::aod::BCsWithTimestamps const& bcs,
                              consumeddata::Run3FullPidTracks const& tracks,
                              consumeddata::Run3Vzeros const& v0s)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processV0s<modes::System::kPP_Run3>(col, tracks, tracksWithItsPid, v0s);
  };
  PROCESS_SWITCH(FemtoProducer, processTracksV0sRun3pp, "Process tracks and v0s", false);

  // process tracks and kinks
  void processTracksKinksRun3pp(consumeddata::Run3PpCollisions::iterator const& col,
                                o2::aod::BCsWithTimestamps const& bcs,
                                consumeddata::Run3FullPidTracks const& tracks,
                                consumeddata::Run3Kinks const& kinks)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processKinks<modes::System::kPP_Run3>(col, tracks, tracksWithItsPid, kinks);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksKinksRun3pp, "Process tracks and kinks", false);

  // process tracks, v0s and cascades
  void processTracksV0sCascadesRun3pp(consumeddata::Run3PpCollisions::iterator const& col,
                                      o2::aod::BCsWithTimestamps const& bcs,
                                      consumeddata::Run3FullPidTracks const& tracks,
                                      consumeddata::Run3Vzeros const& v0s,
                                      consumeddata::Run3Cascades const& cascades)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processV0s<modes::System::kPP_Run3>(col, tracks, tracksWithItsPid, v0s);
    processCascades<modes::System::kPP_Run3>(col, tracks, tracksWithItsPid, cascades);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesRun3pp, "Provide Tracks, V0s and Cascades for Run3", false);

  // process tracks, v0s, cascades and kinks
  void processTracksV0sCascadesKinksRun3pp(consumeddata::Run3PpCollisions::iterator const& col,
                                           o2::aod::BCsWithTimestamps const& bcs,
                                           consumeddata::Run3FullPidTracks const& tracks,
                                           consumeddata::Run3Vzeros const& v0s,
                                           consumeddata::Run3Cascades const& cascades,
                                           consumeddata::Run3Kinks const& kinks)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processV0s<modes::System::kPP_Run3>(col, tracks, tracksWithItsPid, v0s);
    processKinks<modes::System::kPP_Run3>(col, tracks, tracksWithItsPid, kinks);
    processCascades<modes::System::kPP_Run3>(col, tracks, tracksWithItsPid, cascades);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesKinksRun3pp, "Provide Tracks, V0s and Cascades for Run3", false);

  // process monte carlo tracks
  void processTracksRun3ppMc(consumeddata::Run3PpMcRecoCollisions::iterator const& col,
                             consumeddata::Run3PpMcGenCollisions const& mcCols,
                             o2::aod::BCsWithTimestamps const& bcs,
                             consumeddata::Run3McRecoTracks const& tracks,
                             consumeddata::Run3McGenParticles const& mcParticles)
  {
    if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracks, mcParticles)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracks, tracksWithItsPid, mcParticles);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksRun3ppMc, "Provide reconstructed and generated Tracks", false);

  // process monte carlo tracks and v0s
  void processTracksV0sRun3ppMc(consumeddata::Run3PpMcRecoCollisions::iterator const& col,
                                consumeddata::Run3PpMcGenCollisions const& mcCols,
                                o2::aod::BCsWithTimestamps const& bcs,
                                consumeddata::Run3McRecoTracks const& tracks,
                                consumeddata::Run3RecoVzeros const& v0s,
                                consumeddata::Run3McGenParticles const& mcParticles)
  {
    if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracks, mcParticles)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracks, tracksWithItsPid, mcParticles);
    processMcV0s<modes::System::kPP_Run3_MC>(col, mcCols, tracks, tracksWithItsPid, v0s, mcParticles);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sRun3ppMc, "Provide reconstructed and generated tracks and v0s", false);

  // process monte carlo tracks and kinks
  void processTracksKinksRun3ppMc(consumeddata::Run3PpMcRecoCollisions::iterator const& col,
                                  consumeddata::Run3PpMcGenCollisions const& mcCols,
                                  o2::aod::BCsWithTimestamps const& bcs,
                                  consumeddata::Run3McRecoTracks const& tracks,
                                  consumeddata::Run3Kinks const& kinks,
                                  consumeddata::Run3McGenParticles const& mcParticles)
  {
    if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracks, mcParticles)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracks, tracksWithItsPid, mcParticles);
    processMcKinks<modes::System::kPP_Run3_MC>(col, mcCols, tracks, tracksWithItsPid, kinks, mcParticles);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksKinksRun3ppMc, "Provide reconstructed and generated tracks and kinks", false);

  // process monte carlo tracks and v0s and kinks (adding cascades later here)
  void processTracksV0sKinksRun3ppMc(consumeddata::Run3PpMcRecoCollisions::iterator const& col,
                                     consumeddata::Run3PpMcGenCollisions const& mcCols,
                                     o2::aod::BCsWithTimestamps const& bcs,
                                     consumeddata::Run3McRecoTracks const& tracks,
                                     consumeddata::Run3RecoVzeros const& v0s,
                                     consumeddata::Run3Kinks const& kinks,
                                     consumeddata::Run3McGenParticles const& mcParticles)
  {
    if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracks, mcParticles)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracks, tracksWithItsPid, mcParticles);
    processMcV0s<modes::System::kPP_Run3_MC>(col, mcCols, tracks, tracksWithItsPid, v0s, mcParticles);
    processMcKinks<modes::System::kPP_Run3_MC>(col, mcCols, tracks, tracksWithItsPid, kinks, mcParticles);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sKinksRun3ppMc, "Provide reconstructed and generated tracks and v0s and kinks", false);

  // process monte carlo tracks and v0s
  void processTracksV0sCascadesRun3ppMc(consumeddata::Run3PpMcRecoCollisions::iterator const& col,
                                        consumeddata::Run3PpMcGenCollisions const& mcCols,
                                        o2::aod::BCsWithTimestamps const& bcs,
                                        consumeddata::Run3McRecoTracks const& tracks,
                                        consumeddata::Run3RecoVzeros const& v0s,
                                        consumeddata::Run3RecoCascades const& cascades,
                                        consumeddata::Run3McGenParticles const& mcParticles)
  {
    if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracks, mcParticles)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracks, tracksWithItsPid, mcParticles);
    processMcV0s<modes::System::kPP_Run3_MC>(col, mcCols, tracks, tracksWithItsPid, v0s, mcParticles);
    processMcCascades<modes::System::kPP_Run3_MC>(col, mcCols, tracks, tracksWithItsPid, cascades, mcParticles);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesRun3ppMc, "Provide reconstructed and generated tracks and v0s", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{adaptAnalysisTask<FemtoProducer>(cfgc)};
  return workflow;
}
