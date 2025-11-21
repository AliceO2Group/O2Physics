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
/// \brief Tasks that produces the all femto tables
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/Femto/Core/cascadeBuilder.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/kinkBuilder.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/twoTrackResonanceBuilder.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGLF/DataModel/LFKinkDecayTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

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
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/OutputObjHeader.h"
#include "Framework/runDataProcessing.h"

#include "fairlogger/Logger.h"

#include <chrono>
#include <cstdint>
#include <string>
#include <unordered_map>

using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

namespace o2::analysis::femto
{
namespace consumeddata
{
using Run3PpCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms>;

using Run3FullPidTracks =
  soa::Join<Tracks, TracksExtra, TracksDCA,
            pidTPCFullEl, pidTPCFullPi, pidTPCFullKa, pidTPCFullPr, pidTPCFullDe, pidTPCFullTr, pidTPCFullHe,
            pidTOFFullEl, pidTOFFullPi, pidTOFFullKa, pidTOFFullPr, pidTOFFullDe, pidTOFFullTr, pidTOFFullHe,
            pidTOFbeta, pidTOFmass>;

using Run3PpVzeros = V0Datas;

using Run3PpCascades = CascDatas;

using Run3PpKinks = KinkCands;

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

  // resonance daughter filters and partitions
  twotrackresonancebuilder::ConfTwoTrackResonanceDaughterFilters confResonanceDaughterFilters;
  // caching and preslicing
  SliceCache cache;
  Preslice<Tracks> perColTracks = track::collisionId;
  Partition<consumeddata::Run3FullPidTracks> partitionPositiveDaughters =
    (track::signed1Pt > 0.f) &&
    (track::pt > confResonanceDaughterFilters.ptMin && track::pt < confResonanceDaughterFilters.ptMax) &&
    (track::eta > confResonanceDaughterFilters.etaMin && track::eta < confResonanceDaughterFilters.etaMax) &&
    (track::phi > confResonanceDaughterFilters.phiMin && track::phi < confResonanceDaughterFilters.phiMax);
  Partition<consumeddata::Run3FullPidTracks> partitionNegativeDaughters =
    (track::signed1Pt < 0.f) &&
    (track::pt > confResonanceDaughterFilters.ptMin && track::pt < confResonanceDaughterFilters.ptMax) &&
    (track::eta > confResonanceDaughterFilters.etaMin && track::eta < confResonanceDaughterFilters.etaMax) &&
    (track::phi > confResonanceDaughterFilters.phiMin && track::phi < confResonanceDaughterFilters.phiMax);

  // resonance builders
  twotrackresonancebuilder::TwoTrackResonanceBuilderProducts twoTrackResonanceBuilderProducts;
  twotrackresonancebuilder::ConfTwoTrackResonanceTables confTwoTrackResonanceTables;
  twotrackresonancebuilder::ConfRhoFilters confRhoFilters;
  twotrackresonancebuilder::ConfRho0Bits confRho0Bits;
  twotrackresonancebuilder::TwoTrackResonanceBuilder<modes::TwoTrackResonance::kRho0, twotrackresonancebuilder::RhoSelHistName> rho0Builder;
  twotrackresonancebuilder::ConfPhiFilters confPhiFilters;
  twotrackresonancebuilder::ConfPhiBits confPhiBits;
  twotrackresonancebuilder::TwoTrackResonanceBuilder<modes::TwoTrackResonance::kPhi, twotrackresonancebuilder::PhiSelHistName> phiBuilder;
  twotrackresonancebuilder::ConfKstarFilters confKstarFilters;
  twotrackresonancebuilder::ConfKstar0Bits confKstar0Bits;
  twotrackresonancebuilder::TwoTrackResonanceBuilder<modes::TwoTrackResonance::kKstar0, twotrackresonancebuilder::Kstar0SelHistName> kstar0Builder;
  twotrackresonancebuilder::TwoTrackResonanceBuilder<modes::TwoTrackResonance::kKstar0Bar, twotrackresonancebuilder::Kstar0barSelHistName> kstar0barBuilder;

  // histogramming
  // add histograms in next iteration
  HistogramRegistry hRegistry{"FemtoProducer", {}, OutputObjHandlingPolicy::AnalysisObject};

  // data members
  Service<o2::ccdb::BasicCCDBManager> ccdb;            /// Accessing the CCDB
  std::unordered_map<int64_t, int64_t> indexMapTracks; // for mapping tracks to lambdas, cascades and resonances

  void init(InitContext& context)
  {
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

    // configure resonance selections
    rho0Builder.init(&hRegistry, confRho0Bits, confRhoFilters, confResonanceDaughterFilters, confTwoTrackResonanceTables, context);
    phiBuilder.init(&hRegistry, confPhiBits, confPhiFilters, confResonanceDaughterFilters, confTwoTrackResonanceTables, context);
    kstar0Builder.init(&hRegistry, confKstar0Bits, confKstarFilters, confResonanceDaughterFilters, confTwoTrackResonanceTables, context);
    kstar0barBuilder.init(&hRegistry, confKstar0Bits, confKstarFilters, confResonanceDaughterFilters, confTwoTrackResonanceTables, context);

    if ((xiBuilder.fillAnyTable() || omegaBuilder.fillAnyTable()) && (!doprocessTracksV0sCascadesRun3pp && !doprocessTracksV0sCascadesKinksRun3pp)) {
      LOG(fatal) << "At least one cascade table is enabled, but wrong process function is enabled. Breaking...";
    }
    if ((lambdaBuilder.fillAnyTable() || antilambdaBuilder.fillAnyTable() || k0shortBuilder.fillAnyTable()) && (!doprocessTracksV0sCascadesRun3pp && !doprocessTracksV0sRun3pp && !doprocessTracksV0sCascadesKinksRun3pp)) {
      LOG(fatal) << "At least one v0 table is enabled, but wrong process function is enabled. Breaking...";
    }
    if ((sigmaBuilder.fillAnyTable() || sigmaPlusBuilder.fillAnyTable()) && (!doprocessTracksKinksRun3pp && !doprocessTracksV0sCascadesKinksRun3pp)) {
      LOG(fatal) << "At least one kink table is enabled, but wrong process function is enabled. Breaking...";
    }
  }

  // Core implementations
  template <modes::System system, typename T1, typename T2, typename T3>
  bool processCollisions(T1 const& col, T2 const& /* bcs*/, T3 const& tracks)
  {
    collisionBuilder.reset();
    auto bc = col.template bc_as<T2>();
    collisionBuilder.initCollision<system>(bc, col, tracks, ccdb, hRegistry);
    if (!collisionBuilder.checkCollision(col)) {
      return false;
    }
    return true;
  }

  template <modes::System system, typename T1, typename T2>
  void processTracks(T1 const& col, T2 const& tracksWithItsPid)
  {
    trackBuilder.fillTracks<system>(col, collisionBuilder, collisionBuilderProducts, tracksWithItsPid, trackBuilderProducts, indexMapTracks);
  }

  template <modes::System system, typename T1, typename T2>
  void processResonances(T1 const& col, T2 const& /*tracks*/)
  {
    auto groupPositiveTracks = partitionPositiveDaughters->sliceByCached(track::collisionId, col.globalIndex(), cache);
    auto groupNegativeTracks = partitionNegativeDaughters->sliceByCached(track::collisionId, col.globalIndex(), cache);
    rho0Builder.fillResonances<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, twoTrackResonanceBuilderProducts, groupPositiveTracks, groupNegativeTracks, trackBuilder, indexMapTracks);
    phiBuilder.fillResonances<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, twoTrackResonanceBuilderProducts, groupPositiveTracks, groupNegativeTracks, trackBuilder, indexMapTracks);
    kstar0Builder.fillResonances<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, twoTrackResonanceBuilderProducts, groupPositiveTracks, groupNegativeTracks, trackBuilder, indexMapTracks);
    kstar0barBuilder.fillResonances<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, twoTrackResonanceBuilderProducts, groupPositiveTracks, groupNegativeTracks, trackBuilder, indexMapTracks);
  }

  // add v0s
  template <modes::System system, typename T1, typename T2, typename T3>
  void processV0s(T1 const& col, T2 const& tracks, T3 const& v0s)
  {
    lambdaBuilder.fillV0s<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, v0builderProducts, v0s, tracks, trackBuilder, indexMapTracks);
    antilambdaBuilder.fillV0s<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, v0builderProducts, v0s, tracks, trackBuilder, indexMapTracks);
    k0shortBuilder.fillV0s<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, v0builderProducts, v0s, tracks, trackBuilder, indexMapTracks);
  }

  // add kinks
  template <modes::System system, typename T1, typename T2, typename T3>
  void processKinks(T1 const& col, T2 const& tracks, T3 const& kinks)
  {
    sigmaBuilder.fillKinks<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, trackBuilder, indexMapTracks);
    sigmaPlusBuilder.fillKinks<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, trackBuilder, indexMapTracks);
  }

  // add cascades
  template <modes::System system, typename T1, typename T2, typename T3>
  void processCascades(T1 const& col, T2 const& tracks, T3 const& cascades)
  {
    xiBuilder.fillCascades<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, cascadeBuilderProducts,
                                   cascades, tracks, trackBuilder, indexMapTracks);
    omegaBuilder.fillCascades<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, cascadeBuilderProducts,
                                      cascades, tracks, trackBuilder, indexMapTracks);
  }

  // proccess functions
  void processTracksRun3pp(consumeddata::Run3PpCollisions::iterator const& col,
                           BCsWithTimestamps const& bcs,
                           consumeddata::Run3FullPidTracks const& tracks)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    indexMapTracks.clear();
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processResonances<modes::System::kPP_Run3>(col, tracks);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksRun3pp, "Process tracks", true);

  // process tracks and v0s
  void processTracksV0sRun3pp(consumeddata::Run3PpCollisions::iterator const& col,
                              BCsWithTimestamps const& bcs,
                              consumeddata::Run3FullPidTracks const& tracks,
                              consumeddata::Run3PpVzeros const& v0s)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    indexMapTracks.clear();
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processResonances<modes::System::kPP_Run3>(col, tracks);
    processV0s<modes::System::kPP_Run3>(col, tracks, v0s);
  };
  PROCESS_SWITCH(FemtoProducer, processTracksV0sRun3pp, "Process tracks and v0s", false);

  // process tracks and kinks
  void processTracksKinksRun3pp(consumeddata::Run3PpCollisions::iterator const& col,
                                BCsWithTimestamps const& bcs,
                                consumeddata::Run3FullPidTracks const& tracks,
                                consumeddata::Run3PpKinks const& kinks)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    indexMapTracks.clear();
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processResonances<modes::System::kPP_Run3>(col, tracks);
    processKinks<modes::System::kPP_Run3>(col, tracks, kinks);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksKinksRun3pp, "Process tracks and kinks", false);

  // process tracks, v0s and cascades
  void processTracksV0sCascadesRun3pp(consumeddata::Run3PpCollisions::iterator const& col,
                                      BCsWithTimestamps const& bcs,
                                      consumeddata::Run3FullPidTracks const& tracks,
                                      consumeddata::Run3PpVzeros const& v0s,
                                      consumeddata::Run3PpCascades const& cascades)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    indexMapTracks.clear();
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processResonances<modes::System::kPP_Run3>(col, tracks);
    processV0s<modes::System::kPP_Run3>(col, tracks, v0s);
    processCascades<modes::System::kPP_Run3>(col, tracks, cascades);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesRun3pp, "Provide Tracks, V0s and Cascades for Run3", false);

  // process tracks, v0s, cascades and kinks
  void processTracksV0sCascadesKinksRun3pp(consumeddata::Run3PpCollisions::iterator const& col,
                                           BCsWithTimestamps const& bcs,
                                           consumeddata::Run3FullPidTracks const& tracks,
                                           consumeddata::Run3PpVzeros const& v0s,
                                           consumeddata::Run3PpCascades const& cascades,
                                           consumeddata::Run3PpKinks const& kinks)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    indexMapTracks.clear();
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processResonances<modes::System::kPP_Run3>(col, tracks);
    processV0s<modes::System::kPP_Run3>(col, tracks, v0s);
    processKinks<modes::System::kPP_Run3>(col, tracks, kinks);
    processCascades<modes::System::kPP_Run3>(col, tracks, cascades);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesKinksRun3pp, "Provide Tracks, V0s and Cascades for Run3", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoProducer>(cfgc)};
  return workflow;
}
