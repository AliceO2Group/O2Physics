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
/// \brief Tasks that produces most femto tables
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/Femto/Core/cascadeBuilder.h"
#include "PWGCF/Femto/Core/charmHadronBuilder.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/kinkBuilder.h"
#include "PWGCF/Femto/Core/mcBuilder.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGLF/DataModel/LFKinkDecayTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

using namespace o2::analysis::femto;

namespace o2::analysis::femto::rawinputs
{
using Run3PpCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::Mults, o2::aod::CentFT0Ms, o2::aod::CentFT0Cs, o2::aod::CentFT0As>;
using Run3PpMcRecoCollisions = o2::soa::Join<Run3PpCollisions, o2::aod::McCollisionLabels>;
using Run3PpMcGenCollisions = o2::soa::Join<o2::aod::McCollisions, o2::aod::MultsExtraMC, o2::aod::McCentFT0Ms>;

using Run3PbPbCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::Mults, o2::aod::CentFT0Ms, o2::aod::CentFT0Cs, o2::aod::CentFT0As>;
using Run3PbPbCollisionsWithEventShape = o2::soa::Join<Run3PbPbCollisions, o2::aod::QvectorFT0CVecs, o2::aod::QvectorFT0AVecs>;

using Run3PbPbMcRecoCollisions = o2::soa::Join<Run3PbPbCollisions, o2::aod::McCollisionLabels>;
using Run3PbPbMcGenCollisions = o2::soa::Join<o2::aod::McCollisions, o2::aod::MultsExtraMC, o2::aod::McCentFT0Cs>;

using Run3FullPidTracks =
  soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA,
            o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr, o2::aod::pidTPCFullDe, o2::aod::pidTPCFullTr, o2::aod::pidTPCFullHe,
            o2::aod::pidTOFFullEl, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr, o2::aod::pidTOFFullDe, o2::aod::pidTOFFullTr, o2::aod::pidTOFFullHe,
            o2::aod::pidTOFbeta, o2::aod::pidTOFmass>;
using Run3McRecoTracks = soa::Join<Run3FullPidTracks, o2::aod::McTrackLabels>;

using Run3Vzeros = o2::soa::Join<o2::aod::V0Datas, o2::aod::V0TOFPIDs, o2::aod::V0TOFNSigmas>;
using Run3RecoVzeros = o2::soa::Join<Run3Vzeros, o2::aod::McV0Labels>;

using Run3D0s = soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>;
using Run3RecoD0s = soa::Join<Run3D0s, aod::HfCand2ProngMcRec>;

using Run3Cascades = o2::soa::Join<o2::aod::CascDatas, o2::aod::CascTOFPIDs, o2::aod::CascTOFNSigmas>;
using Run3RecoCascades = o2::soa::Join<Run3Cascades, o2::aod::McCascLabels>;

using Run3Kinks = o2::aod::KinkCands;

using Run3McGenParticles = o2::aod::McParticles;

} // namespace o2::analysis::femto::rawinputs

struct FemtoProducer {

  o2::framework::Preslice<rawinputs::Run3McGenParticles> perMcCollision = o2::aod::mcparticle::mcCollisionId;
  o2::framework::Preslice<rawinputs::Run3McRecoTracks> perColRecoTracks = o2::aod::track::collisionId;
  o2::framework::Preslice<rawinputs::Run3RecoVzeros> perColRecoV0s = o2::aod::v0data::collisionId;
  o2::framework::Preslice<rawinputs::Run3RecoCascades> perColRecoCascades = o2::aod::cascdata::collisionId;
  o2::framework::Preslice<rawinputs::Run3Kinks> perColRecoKinks = o2::aod::kinkcand::collisionId;
  o2::framework::Preslice<rawinputs::Run3RecoD0s> perColRecoD0s = o2::aod::hf_cand::collisionId;

  // ccdb config
  collisionbuilder::ConfCcdb confCcdb;

  // collision builder
  collisionbuilder::CollisionBuilderProducts collisionBuilderProducts;
  collisionbuilder::ConfCollisionTables confCollisionTables;
  collisionbuilder::ConfCollisionFilters confCollisionFilters;
  collisionbuilder::ConfCollisionBits confCollisionBits;
  collisionbuilder::ConfCollisionRctFlags confCollisionRctFlags;
  collisionbuilder::CollisionBuilder<collisionbuilder::ColSelHistName, collisionbuilder::CollisionFilterHistName> collisionBuilder;

  // track builder
  trackbuilder::TrackBuilderProducts trackBuilderProducts;
  trackbuilder::ConfTrackTables confTrackTables;
  trackbuilder::TrackBuilder<trackbuilder::TrackSelHistName, trackbuilder::TrackFilterHistName> trackBuilder;
  trackbuilder::ConfTrackBits confTrackBits;
  trackbuilder::ConfTrackFilters confTrackFilters;

  // v0 builders
  v0builder::V0BuilderProducts v0builderProducts;
  v0builder::ConfV0Tables confV0Tables;
  v0builder::ConfV0Filters confV0Filters;
  v0builder::ConfK0shortBits confK0shortBits;
  v0builder::V0Builder<modes::V0::kK0short, v0builder::K0shortSelHistName, v0builder::K0shortFilterHistName> k0shortBuilder;
  v0builder::ConfLambdaBits confLambdaBits;
  v0builder::V0Builder<modes::V0::kLambda, v0builder::LambdaSelHistName, v0builder::LambdaFilterHistName> lambdaBuilder;
  v0builder::V0Builder<modes::V0::kAntiLambda, v0builder::AntilambdaSelHistName, v0builder::AntiLambdaFilterHistName> antilambdaBuilder;

  // charm hadron builder
  charmhadronbuilder::CharmHadronBuilderProducts charmHadronBuilderProducts;
  charmhadronbuilder::ConfD0Filters confD0Filters;
  charmhadronbuilder::ConfD0Bits confD0Bits;
  charmhadronbuilder::ConfD0Tables confD0Tables;
  charmhadronbuilder::CharmHadronBuilder<modes::CharmHadron::kD0, charmhadronbuilder::D0SelHistName, charmhadronbuilder::D0FilterHistName> d0Builder;
  charmhadronbuilder::CharmHadronBuilder<modes::CharmHadron::kD0Bar, charmhadronbuilder::D0barSelHistName, charmhadronbuilder::D0barFilterHistName> d0barBuilder;

  // cascade builder
  cascadebuilder::CascadeBuilderProducts cascadeBuilderProducts;
  cascadebuilder::ConfCascadeTables confCascadeTables;
  cascadebuilder::ConfCascadeFilters confCascadeFilters;
  cascadebuilder::ConfXiBits confXiBits;
  cascadebuilder::CascadeBuilder<modes::Cascade::kXi, cascadebuilder::XiSelHistName, cascadebuilder::XiFilterHistName> xiBuilder;
  cascadebuilder::ConfOmegaBits confOmegaBits;
  cascadebuilder::CascadeBuilder<modes::Cascade::kOmega, cascadebuilder::OmegaSelHistName, cascadebuilder::OmegaFilterHistName> omegaBuilder;

  // kink builder
  kinkbuilder::KinkBuilderProducts kinkBuilderProducts;
  kinkbuilder::ConfKinkTables confKinkTables;
  kinkbuilder::ConfKinkFilters confKinkFilters;
  kinkbuilder::ConfSigmaBits confSigmaBits;
  kinkbuilder::KinkBuilder<modes::Kink::kSigma, kinkbuilder::SigmaSelHistName, kinkbuilder::SigmaFilterHistName> sigmaBuilder;
  kinkbuilder::ConfSigmaPlusBits confSigmaPlusBits;
  kinkbuilder::KinkBuilder<modes::Kink::kSigmaPlus, kinkbuilder::SigmaPlusSelHistName, kinkbuilder::SigmaPlusFilterHistName> sigmaPlusBuilder;

  // mc builder
  mcbuilder::ConfMc confMc;
  mcbuilder::ConfMcTables confMcTables;
  mcbuilder::McBuilderProducts mcProducts;
  mcbuilder::McBuilder mcBuilder;

  // histogramming
  o2::framework::HistogramRegistry hRegistry{"FemtoProducer", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // data members
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb = {}; /// Accessing the CCDB
  o2::framework::Service<o2::framework::O2DatabasePDG> pdgDb = {};

  void init(o2::framework::InitContext& context)
  {
    // ---- ccdb ---------------------------------------------------------------
    LOG(info) << "Setting up connection to CCDB with URL: " << confCcdb.ccdbUrl.value;
    ccdb->setURL(confCcdb.ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    // ---- builders -----------------------------------------------------------
    // Order matters here: CollisionBuilder must init before TrackBuilder (which checks
    // collisionBuilder.producingCollisions()/producingLiteCollisions()), and TrackBuilder
    // must init before the V0/Cascade/Kink builders (which check
    // trackBuilder.producingTracks()/producingLiteTracks()).
    collisionBuilder.init(&hRegistry, confCollisionFilters, confCollisionBits, confCollisionRctFlags, confCcdb, confCollisionTables, context);
    trackBuilder.init(&hRegistry, confTrackBits, confTrackFilters, confTrackTables, context, collisionBuilder);

    k0shortBuilder.init(&hRegistry, confK0shortBits, confV0Filters, confV0Tables, context, trackBuilder);
    lambdaBuilder.init(&hRegistry, confLambdaBits, confV0Filters, confV0Tables, context, trackBuilder);
    antilambdaBuilder.init(&hRegistry, confLambdaBits, confV0Filters, confV0Tables, context, trackBuilder);

    d0Builder.init(&hRegistry, confD0Bits, confD0Filters, confD0Tables, context);
    d0barBuilder.init(&hRegistry, confD0Bits, confD0Filters, confD0Tables, context);

    sigmaBuilder.init(&hRegistry, confSigmaBits, confKinkFilters, confKinkTables, context, trackBuilder);
    sigmaPlusBuilder.init(&hRegistry, confSigmaPlusBits, confKinkFilters, confKinkTables, context, trackBuilder);

    xiBuilder.init(&hRegistry, confXiBits, confCascadeFilters, confCascadeTables, context, trackBuilder);
    omegaBuilder.init(&hRegistry, confOmegaBits, confCascadeFilters, confCascadeTables, context, trackBuilder);

    mcBuilder.init(confMc, confMcTables, context);

    // ---- guard: enabled tables vs. enabled process function ------------------
    if ((xiBuilder.fillAnyTable() || omegaBuilder.fillAnyTable()) &&
        (!doprocessTracksV0sCascadesRun3pp && !doprocessTracksV0sCascadesRun3PbPb &&
         !doprocessTracksV0sCascadesKinksRun3pp && !doprocessTracksV0sCascadesRun3ppMc &&
         !doprocessTracksV0sCascadesRun3PbPbMc)) {
      LOG(fatal) << "At least one cascade table is enabled, but wrong process function is enabled. Breaking...";
    }
    if ((lambdaBuilder.fillAnyTable() || antilambdaBuilder.fillAnyTable() || k0shortBuilder.fillAnyTable()) &&
        (!doprocessTracksV0sCascadesRun3pp && !doprocessTracksV0sCascadesRun3PbPb &&
         !doprocessTracksV0sRun3pp && !doprocessTracksV0sCascadesKinksRun3pp &&
         !doprocessTracksV0sRun3ppMc && !doprocessTracksV0sRun3PbPb && !doprocessTracksV0sRun3PbPbMc &&
         !doprocessTracksV0sCascadesRun3ppMc && !doprocessTracksV0sCascadesRun3PbPbMc &&
         !doprocessTracksV0sKinksRun3ppMc)) {
      LOG(fatal) << "At least one v0 table is enabled, but wrong process function is enabled. Breaking...";
    }
    if ((sigmaBuilder.fillAnyTable() || sigmaPlusBuilder.fillAnyTable()) &&
        (!doprocessTracksKinksRun3pp && !doprocessTracksV0sCascadesKinksRun3pp &&
         !doprocessTracksKinksRun3ppMc && !doprocessTracksV0sKinksRun3ppMc)) {
      LOG(fatal) << "At least one kink table is enabled, but wrong process function is enabled. Breaking...";
    }
    if ((d0Builder.fillAnyTable() || d0barBuilder.fillAnyTable()) &&
        (!doprocessTracksD0sRun3pp && !doprocessTracksD0sRun3PbPb &&
         !doprocessTracksD0sRun3ppMc && !doprocessTracksD0sRun3PbPbMc)) {
      LOG(fatal) << "At least one d0 table is enabled, but wrong process function is enabled. Breaking...";
    }
    if (mcBuilder.fillAnyTable() &&
        (!doprocessTracksRun3ppMc && !doprocessTracksRun3PbPbMc &&
         !doprocessTracksV0sRun3ppMc && !doprocessTracksV0sRun3PbPbMc &&
         !doprocessTracksV0sCascadesRun3ppMc && !doprocessTracksV0sCascadesRun3PbPbMc &&
         !doprocessTracksKinksRun3ppMc && !doprocessTracksV0sKinksRun3ppMc &&
         !doprocessTracksD0sRun3ppMc && !doprocessTracksD0sRun3PbPbMc &&
         !doprocessMcOnly)) {
      LOG(fatal) << "At least one mc table is enabled, but wrong process function is enabled. Breaking...";
    }

    // ---- guard: pass-through consistency across builders ---------------------
    {
      std::vector<std::pair<std::string, bool>> passThroughFlags;
      auto add = [&passThroughFlags](const std::string& name, bool enabled, bool passThrough) {
        if (enabled) {
          passThroughFlags.emplace_back(name, passThrough);
        }
      };
      // substitute the real accessors if they live on the builders rather than the conf groups
      add("collision", collisionBuilder.fillAnyTable(), collisionBuilder.isPassThrough());
      add("track", trackBuilder.fillAnyTable(), trackBuilder.isPassThrough());
      add("k0short", k0shortBuilder.fillAnyTable(), k0shortBuilder.isPassThrough());
      add("lambda", lambdaBuilder.fillAnyTable(), lambdaBuilder.isPassThrough());
      add("antilambda", antilambdaBuilder.fillAnyTable(), antilambdaBuilder.isPassThrough());
      add("d0", d0Builder.fillAnyTable(), d0Builder.isPassThrough());
      add("d0bar", d0barBuilder.fillAnyTable(), d0barBuilder.isPassThrough());
      add("sigma", sigmaBuilder.fillAnyTable(), sigmaBuilder.isPassThrough());
      add("sigmaplus", sigmaPlusBuilder.fillAnyTable(), sigmaPlusBuilder.isPassThrough());
      add("xi", xiBuilder.fillAnyTable(), xiBuilder.isPassThrough());
      add("omega", omegaBuilder.fillAnyTable(), omegaBuilder.isPassThrough());
      add("mc", mcBuilder.fillAnyTable(), mcBuilder.isPassThrough());

      const bool anyPassThrough = std::any_of(passThroughFlags.begin(), passThroughFlags.end(),
                                              [](const auto& f) { return f.second; });
      const bool allPassThrough = std::all_of(passThroughFlags.begin(), passThroughFlags.end(),
                                              [](const auto& f) { return f.second; });

      if (anyPassThrough && !allPassThrough) {
        LOG(warn) << "Builders disagree about pass-through mode:";
        for (const auto& [name, passThrough] : passThroughFlags) {
          LOG(warn) << "  " << name << " builder: pass-through = " << (passThrough ? "true" : "false");
        }
        LOG(warn) << "This is legal, but check it is intended:";
        LOG(warn) << "  - collisions in pass-through, particles not: every collision is written, "
                  << "but the particle tables only cover collisions with a candidate.";
        LOG(warn) << "  - particles in pass-through, collisions not: the extra particle rows hang off "
                  << "a candidate-biased event sample, which will bias event-normalised observables.";
      }
    }

    // ---- guard: exactly one process function ---------------------------------
    const int nProcesses =
      static_cast<int>(doprocessTracksRun3pp) +
      static_cast<int>(doprocessTracksRun3PbPb) +
      static_cast<int>(doprocessTracksRun3PbPbWithEventShape) +
      static_cast<int>(doprocessTracksV0sRun3pp) +
      static_cast<int>(doprocessTracksV0sRun3PbPb) +
      static_cast<int>(doprocessTracksV0sCascadesRun3pp) +
      static_cast<int>(doprocessTracksV0sCascadesRun3PbPb) +
      static_cast<int>(doprocessTracksKinksRun3pp) +
      static_cast<int>(doprocessTracksV0sCascadesKinksRun3pp) +
      static_cast<int>(doprocessTracksD0sRun3pp) +
      static_cast<int>(doprocessTracksD0sRun3PbPb) +
      static_cast<int>(doprocessTracksRun3ppMc) +
      static_cast<int>(doprocessTracksRun3PbPbMc) +
      static_cast<int>(doprocessTracksV0sRun3ppMc) +
      static_cast<int>(doprocessTracksV0sRun3PbPbMc) +
      static_cast<int>(doprocessTracksV0sCascadesRun3ppMc) +
      static_cast<int>(doprocessTracksV0sCascadesRun3PbPbMc) +
      static_cast<int>(doprocessTracksKinksRun3ppMc) +
      static_cast<int>(doprocessTracksV0sKinksRun3ppMc) +
      static_cast<int>(doprocessTracksD0sRun3ppMc) +
      static_cast<int>(doprocessTracksD0sRun3PbPbMc) +
      static_cast<int>(doprocessMcOnly);

    if (nProcesses != 1) {
      LOG(fatal) << "Exactly one process function must be enabled, found " << nProcesses << ". Breaking...";
    }

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
    // Normally the collision row is written lazily, once a candidate is found. In
    // pass-through we need every accepted collision, including those with no candidate,
    // otherwise event-normalised observables lose exactly those events.
    if (collisionBuilder.isPassThrough()) {
      collisionBuilder.fillCollision<system>(collisionBuilderProducts, col);
    }
    return true;
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  bool processMcCollisions(T1 const& col, T2 const& mcCols, T3 const& /* bcs*/, T4 const& tracks, T5 const& mcParticles)
  {
    collisionBuilder.reset();
    // in pass-through the maps are built once per timeframe
    // resetting here would wipe them and duplicate every mc collision and particle
    if (!mcBuilder.isPassThrough()) {
      mcBuilder.reset(mcCols, mcParticles);
    }
    auto bc = col.template bc_as<T3>();
    collisionBuilder.initCollision<system>(bc, col, tracks, ccdb, hRegistry);
    if (!collisionBuilder.checkCollision(col, mcCols)) {
      return false;
    }
    // NOTE: must be fillMcCollision, not fillCollision - fillMcCollision returns early on
    //       mCollisionAlreadyFilled, so filling the collision first would permanently skip
    //       the collision label
    if (collisionBuilder.isPassThrough()) {
      collisionBuilder.fillMcCollision<system>(collisionBuilderProducts, col, mcCols, mcProducts, mcBuilder);
    }
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
  template <modes::System system, typename T1, typename T2, typename T3>
  void processV0s(T1 const& col, T2 const& tracks, T3 const& v0s)
  {
    lambdaBuilder.fillV0s<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, v0builderProducts, v0s, tracks, trackBuilder);
    antilambdaBuilder.fillV0s<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, v0builderProducts, v0s, tracks, trackBuilder);
    k0shortBuilder.fillV0s<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, v0builderProducts, v0s, tracks, trackBuilder);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processMcV0s(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& v0s, T5 const& mcParticles)
  {
    lambdaBuilder.fillMcV0s<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, v0builderProducts, v0s, tracks, trackBuilder, mcParticles, mcBuilder, mcProducts);
    antilambdaBuilder.fillMcV0s<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, v0builderProducts, v0s, tracks, trackBuilder, mcParticles, mcBuilder, mcProducts);
    k0shortBuilder.fillMcV0s<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, v0builderProducts, v0s, tracks, trackBuilder, mcParticles, mcBuilder, mcProducts);
  }

  // processing kinks
  template <modes::System system, typename T1, typename T2, typename T3>
  void processKinks(T1 const& col, T2 const& tracks, T3 const& kinks)
  {
    sigmaBuilder.fillKinks<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, trackBuilder);
    sigmaPlusBuilder.fillKinks<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, trackBuilder);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processMcKinks(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& kinks, T5 const& mcParticles)
  {
    sigmaBuilder.fillMcKinks<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, trackBuilder, mcParticles, mcBuilder, mcProducts);
    sigmaPlusBuilder.fillMcKinks<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, trackBuilder, mcParticles, mcBuilder, mcProducts);
  }

  // processing cascades
  template <modes::System system, typename T1, typename T2, typename T3>
  void processCascades(T1 const& col, T2 const& tracks, T3 const& cascades)
  {
    xiBuilder.fillCascades<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, cascadeBuilderProducts,
                                   cascades, tracks, trackBuilder);
    omegaBuilder.fillCascades<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, cascadeBuilderProducts,
                                      cascades, tracks, trackBuilder);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processMcCascades(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& cascades, T5 const& mcParticles)
  {
    xiBuilder.fillMcCascades<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, cascadeBuilderProducts, cascades, tracks, trackBuilder, mcParticles, mcBuilder, mcProducts);
    omegaBuilder.fillMcCascades<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, cascadeBuilderProducts, cascades, tracks, trackBuilder, mcParticles, mcBuilder, mcProducts);
  }

  // processing D0s
  template <modes::System system, typename T1, typename T2, typename T3>
  void processD0s(T1 const& col, T2 const& tracks, T3 const& candidates)
  {
    d0Builder.fillD0s<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, charmHadronBuilderProducts, candidates, tracks, trackBuilder);
    d0barBuilder.fillD0s<system>(col, collisionBuilder, collisionBuilderProducts, trackBuilderProducts, charmHadronBuilderProducts, candidates, tracks, trackBuilder);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processMcD0s(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& candidates, T5 const& mcParticles)
  {
    d0Builder.fillMcD0s<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, charmHadronBuilderProducts, candidates, tracks, trackBuilder, mcParticles, mcBuilder, mcProducts);
    d0barBuilder.fillMcD0s<system>(col, collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, charmHadronBuilderProducts, candidates, tracks, trackBuilder, mcParticles, mcBuilder, mcProducts);
  }

  // ==========================================================================
  //  process functions - data
  // ==========================================================================

  void processTracksRun3pp(rawinputs::Run3PpCollisions::iterator const& col,
                           o2::aod::BCsWithTimestamps const& bcs,
                           rawinputs::Run3FullPidTracks const& tracks)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksRun3pp, "Provide tracks", true);

  void processTracksRun3PbPb(rawinputs::Run3PbPbCollisions::iterator const& col,
                             o2::aod::BCsWithTimestamps const& bcs,
                             rawinputs::Run3FullPidTracks const& tracks)
  {
    if (!processCollisions<modes::System::kPbPb_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPbPb_Run3>(col, tracksWithItsPid);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksRun3PbPb, "Provide tracks in PbPb collisions", false);

  void processTracksRun3PbPbWithEventShape(rawinputs::Run3PbPbCollisionsWithEventShape::iterator const& col,
                                           o2::aod::BCsWithTimestamps const& bcs,
                                           rawinputs::Run3FullPidTracks const& tracks)
  {
    if (!processCollisions<modes::System::kPbPb_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPbPb_Run3>(col, tracksWithItsPid);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksRun3PbPbWithEventShape, "Provide tracks in PbPb collisions with event shape information", false);

  void processTracksV0sRun3pp(rawinputs::Run3PpCollisions::iterator const& col,
                              o2::aod::BCsWithTimestamps const& bcs,
                              rawinputs::Run3FullPidTracks const& tracks,
                              rawinputs::Run3Vzeros const& v0s)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processV0s<modes::System::kPP_Run3>(col, tracks, v0s);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sRun3pp, "Provide tracks and v0s", false);

  void processTracksV0sRun3PbPb(rawinputs::Run3PbPbCollisions::iterator const& col,
                                o2::aod::BCsWithTimestamps const& bcs,
                                rawinputs::Run3FullPidTracks const& tracks,
                                rawinputs::Run3Vzeros const& v0s)
  {
    if (!processCollisions<modes::System::kPbPb_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPbPb_Run3>(col, tracksWithItsPid);
    processV0s<modes::System::kPbPb_Run3>(col, tracks, v0s);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sRun3PbPb, "Provide tracks and v0s in PbPb collisions", false);

  void processTracksV0sCascadesRun3pp(rawinputs::Run3PpCollisions::iterator const& col,
                                      o2::aod::BCsWithTimestamps const& bcs,
                                      rawinputs::Run3FullPidTracks const& tracks,
                                      rawinputs::Run3Vzeros const& v0s,
                                      rawinputs::Run3Cascades const& cascades)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processV0s<modes::System::kPP_Run3>(col, tracks, v0s);
    processCascades<modes::System::kPP_Run3>(col, tracks, cascades);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesRun3pp, "Provide tracks, v0s and cascades", false);

  void processTracksV0sCascadesRun3PbPb(rawinputs::Run3PbPbCollisions::iterator const& col,
                                        o2::aod::BCsWithTimestamps const& bcs,
                                        rawinputs::Run3FullPidTracks const& tracks,
                                        rawinputs::Run3Vzeros const& v0s,
                                        rawinputs::Run3Cascades const& cascades)
  {
    if (!processCollisions<modes::System::kPbPb_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPbPb_Run3>(col, tracksWithItsPid);
    processV0s<modes::System::kPbPb_Run3>(col, tracks, v0s);
    processCascades<modes::System::kPbPb_Run3>(col, tracks, cascades);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesRun3PbPb, "Provide tracks, v0s and cascades in PbPb collisions", false);

  void processTracksKinksRun3pp(rawinputs::Run3PpCollisions::iterator const& col,
                                o2::aod::BCsWithTimestamps const& bcs,
                                rawinputs::Run3FullPidTracks const& tracks,
                                rawinputs::Run3Kinks const& kinks)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processKinks<modes::System::kPP_Run3>(col, tracks, kinks);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksKinksRun3pp, "Provide tracks and kinks", false);

  void processTracksV0sCascadesKinksRun3pp(rawinputs::Run3PpCollisions::iterator const& col,
                                           o2::aod::BCsWithTimestamps const& bcs,
                                           rawinputs::Run3FullPidTracks const& tracks,
                                           rawinputs::Run3Vzeros const& v0s,
                                           rawinputs::Run3Cascades const& cascades,
                                           rawinputs::Run3Kinks const& kinks)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processV0s<modes::System::kPP_Run3>(col, tracks, v0s);
    processCascades<modes::System::kPP_Run3>(col, tracks, cascades);
    processKinks<modes::System::kPP_Run3>(col, tracks, kinks);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesKinksRun3pp, "Provide tracks, v0s, cascades and kinks", false);

  void processTracksD0sRun3pp(rawinputs::Run3PpCollisions::iterator const& col,
                              o2::aod::BCsWithTimestamps const& bcs,
                              rawinputs::Run3FullPidTracks const& tracks,
                              rawinputs::Run3D0s const& candidates)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processD0s<modes::System::kPP_Run3>(col, tracks, candidates);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksD0sRun3pp, "Provide tracks and D0s", false);

  void processTracksD0sRun3PbPb(rawinputs::Run3PbPbCollisions::iterator const& col,
                                o2::aod::BCsWithTimestamps const& bcs,
                                rawinputs::Run3FullPidTracks const& tracks,
                                rawinputs::Run3D0s const& candidates)
  {
    if (!processCollisions<modes::System::kPbPb_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPbPb_Run3>(col, tracksWithItsPid);
    processD0s<modes::System::kPbPb_Run3>(col, tracks, candidates);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksD0sRun3PbPb, "Provide tracks and D0s in PbPb collisions", false);

  // ==========================================================================
  //  process functions - monte carlo (ungrouped)
  //
  //  No argument is an ::iterator, so DPL hands over the full timeframe tables and
  //  does no auto-grouping. Tables that are iterated over (tracks, v0s, cascades,
  //  kinks) are sliced per collision; the full track table is passed alongside so
  //  the v0/cascade/kink builders can resolve their daughter indices.
  //  These functions own the per-timeframe reset of the mc builder maps, which is
  //  why processMcCollisions skips its own reset in pass-through mode.
  // ==========================================================================

  void processTracksRun3ppMc(rawinputs::Run3PpMcGenCollisions const& mcCols,
                             rawinputs::Run3PpMcRecoCollisions const& cols,
                             o2::aod::BCsWithTimestamps const& bcs,
                             rawinputs::Run3McRecoTracks const& tracks,
                             rawinputs::Run3McGenParticles const& mcParticles)
  {
    if (mcBuilder.isPassThrough()) {
      mcBuilder.reset(mcCols, mcParticles);
      mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, perMcCollision, mcProducts, pdgDb);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol, mcParticles)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksRun3ppMc, "Provide reconstructed and generated tracks", false);

  void processTracksRun3PbPbMc(rawinputs::Run3PbPbMcGenCollisions const& mcCols,
                               rawinputs::Run3PbPbMcRecoCollisions const& cols,
                               o2::aod::BCsWithTimestamps const& bcs,
                               rawinputs::Run3McRecoTracks const& tracks,
                               rawinputs::Run3McGenParticles const& mcParticles)
  {
    if (mcBuilder.isPassThrough()) {
      mcBuilder.reset(mcCols, mcParticles);
      mcBuilder.fillMcPassThrough<modes::System::kPbPb_Run3_MC>(mcCols, mcParticles, perMcCollision, mcProducts, pdgDb);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPbPb_Run3_MC>(col, mcCols, bcs, tracksThisCol, mcParticles)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPbPb_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksRun3PbPbMc, "Provide reconstructed and generated tracks in PbPb collisions", false);

  void processTracksV0sRun3ppMc(rawinputs::Run3PpMcGenCollisions const& mcCols,
                                rawinputs::Run3PpMcRecoCollisions const& cols,
                                o2::aod::BCsWithTimestamps const& bcs,
                                rawinputs::Run3McRecoTracks const& tracks,
                                rawinputs::Run3RecoVzeros const& v0s,
                                rawinputs::Run3McGenParticles const& mcParticles)
  {
    if (mcBuilder.isPassThrough()) {
      mcBuilder.reset(mcCols, mcParticles);
      mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, perMcCollision, mcProducts, pdgDb);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol, mcParticles)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto v0sThisCol = v0s.sliceBy(perColRecoV0s, col.globalIndex());
      processMcV0s<modes::System::kPP_Run3_MC>(col, mcCols, tracks, v0sThisCol, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sRun3ppMc, "Provide reconstructed and generated tracks and v0s", false);

  void processTracksV0sRun3PbPbMc(rawinputs::Run3PbPbMcGenCollisions const& mcCols,
                                  rawinputs::Run3PbPbMcRecoCollisions const& cols,
                                  o2::aod::BCsWithTimestamps const& bcs,
                                  rawinputs::Run3McRecoTracks const& tracks,
                                  rawinputs::Run3RecoVzeros const& v0s,
                                  rawinputs::Run3McGenParticles const& mcParticles)
  {
    if (mcBuilder.isPassThrough()) {
      mcBuilder.reset(mcCols, mcParticles);
      mcBuilder.fillMcPassThrough<modes::System::kPbPb_Run3_MC>(mcCols, mcParticles, perMcCollision, mcProducts, pdgDb);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPbPb_Run3_MC>(col, mcCols, bcs, tracksThisCol, mcParticles)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPbPb_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto v0sThisCol = v0s.sliceBy(perColRecoV0s, col.globalIndex());
      processMcV0s<modes::System::kPbPb_Run3_MC>(col, mcCols, tracks, v0sThisCol, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sRun3PbPbMc, "Provide reconstructed and generated tracks and v0s in PbPb collisions", false);

  void processTracksV0sCascadesRun3ppMc(rawinputs::Run3PpMcGenCollisions const& mcCols,
                                        rawinputs::Run3PpMcRecoCollisions const& cols,
                                        o2::aod::BCsWithTimestamps const& bcs,
                                        rawinputs::Run3McRecoTracks const& tracks,
                                        rawinputs::Run3RecoVzeros const& v0s,
                                        rawinputs::Run3RecoCascades const& cascades,
                                        rawinputs::Run3McGenParticles const& mcParticles)
  {
    if (mcBuilder.isPassThrough()) {
      mcBuilder.reset(mcCols, mcParticles);
      mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, perMcCollision, mcProducts, pdgDb);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol, mcParticles)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto v0sThisCol = v0s.sliceBy(perColRecoV0s, col.globalIndex());
      processMcV0s<modes::System::kPP_Run3_MC>(col, mcCols, tracks, v0sThisCol, mcParticles);

      auto cascadesThisCol = cascades.sliceBy(perColRecoCascades, col.globalIndex());
      processMcCascades<modes::System::kPP_Run3_MC>(col, mcCols, tracks, cascadesThisCol, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesRun3ppMc, "Provide reconstructed and generated tracks, v0s and cascades", false);

  void processTracksV0sCascadesRun3PbPbMc(rawinputs::Run3PbPbMcGenCollisions const& mcCols,
                                          rawinputs::Run3PbPbMcRecoCollisions const& cols,
                                          o2::aod::BCsWithTimestamps const& bcs,
                                          rawinputs::Run3McRecoTracks const& tracks,
                                          rawinputs::Run3RecoVzeros const& v0s,
                                          rawinputs::Run3RecoCascades const& cascades,
                                          rawinputs::Run3McGenParticles const& mcParticles)
  {
    if (mcBuilder.isPassThrough()) {
      mcBuilder.reset(mcCols, mcParticles);
      mcBuilder.fillMcPassThrough<modes::System::kPbPb_Run3_MC>(mcCols, mcParticles, perMcCollision, mcProducts, pdgDb);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPbPb_Run3_MC>(col, mcCols, bcs, tracksThisCol, mcParticles)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPbPb_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto v0sThisCol = v0s.sliceBy(perColRecoV0s, col.globalIndex());
      processMcV0s<modes::System::kPbPb_Run3_MC>(col, mcCols, tracks, v0sThisCol, mcParticles);

      auto cascadesThisCol = cascades.sliceBy(perColRecoCascades, col.globalIndex());
      processMcCascades<modes::System::kPbPb_Run3_MC>(col, mcCols, tracks, cascadesThisCol, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesRun3PbPbMc, "Provide reconstructed and generated tracks, v0s and cascades in PbPb collisions", false);

  void processTracksKinksRun3ppMc(rawinputs::Run3PpMcGenCollisions const& mcCols,
                                  rawinputs::Run3PpMcRecoCollisions const& cols,
                                  o2::aod::BCsWithTimestamps const& bcs,
                                  rawinputs::Run3McRecoTracks const& tracks,
                                  rawinputs::Run3Kinks const& kinks,
                                  rawinputs::Run3McGenParticles const& mcParticles)
  {
    if (mcBuilder.isPassThrough()) {
      mcBuilder.reset(mcCols, mcParticles);
      mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, perMcCollision, mcProducts, pdgDb);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol, mcParticles)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto kinksThisCol = kinks.sliceBy(perColRecoKinks, col.globalIndex());
      processMcKinks<modes::System::kPP_Run3_MC>(col, mcCols, tracks, kinksThisCol, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksKinksRun3ppMc, "Provide reconstructed and generated tracks and kinks", false);

  void processTracksV0sKinksRun3ppMc(rawinputs::Run3PpMcGenCollisions const& mcCols,
                                     rawinputs::Run3PpMcRecoCollisions const& cols,
                                     o2::aod::BCsWithTimestamps const& bcs,
                                     rawinputs::Run3McRecoTracks const& tracks,
                                     rawinputs::Run3RecoVzeros const& v0s,
                                     rawinputs::Run3Kinks const& kinks,
                                     rawinputs::Run3McGenParticles const& mcParticles)
  {
    if (mcBuilder.isPassThrough()) {
      mcBuilder.reset(mcCols, mcParticles);
      mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, perMcCollision, mcProducts, pdgDb);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol, mcParticles)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto v0sThisCol = v0s.sliceBy(perColRecoV0s, col.globalIndex());
      processMcV0s<modes::System::kPP_Run3_MC>(col, mcCols, tracks, v0sThisCol, mcParticles);

      auto kinksThisCol = kinks.sliceBy(perColRecoKinks, col.globalIndex());
      processMcKinks<modes::System::kPP_Run3_MC>(col, mcCols, tracks, kinksThisCol, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sKinksRun3ppMc, "Provide reconstructed and generated tracks, v0s and kinks", false);

  void processTracksD0sRun3ppMc(rawinputs::Run3PpMcGenCollisions const& mcCols,
                                rawinputs::Run3PpMcRecoCollisions const& cols,
                                o2::aod::BCsWithTimestamps const& bcs,
                                rawinputs::Run3McRecoTracks const& tracks,
                                rawinputs::Run3RecoD0s const& d0s,
                                rawinputs::Run3McGenParticles const& mcParticles)
  {
    if (mcBuilder.isPassThrough()) {
      mcBuilder.reset(mcCols, mcParticles);
      mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, perMcCollision, mcProducts, pdgDb);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol, mcParticles)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto d0sThisCol = d0s.sliceBy(perColRecoD0s, col.globalIndex());
      processMcD0s<modes::System::kPP_Run3_MC>(col, mcCols, tracks, d0sThisCol, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksD0sRun3ppMc, "Provide reconstructed and generated tracks and D0s", false);

  void processTracksD0sRun3PbPbMc(rawinputs::Run3PbPbMcGenCollisions const& mcCols,
                                  rawinputs::Run3PbPbMcRecoCollisions const& cols,
                                  o2::aod::BCsWithTimestamps const& bcs,
                                  rawinputs::Run3McRecoTracks const& tracks,
                                  rawinputs::Run3RecoD0s const& d0s,
                                  rawinputs::Run3McGenParticles const& mcParticles)
  {
    if (mcBuilder.isPassThrough()) {
      mcBuilder.reset(mcCols, mcParticles);
      mcBuilder.fillMcPassThrough<modes::System::kPbPb_Run3_MC>(mcCols, mcParticles, perMcCollision, mcProducts, pdgDb);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPbPb_Run3_MC>(col, mcCols, bcs, tracksThisCol, mcParticles)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPbPb_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto d0sThisCol = d0s.sliceBy(perColRecoD0s, col.globalIndex());
      processMcD0s<modes::System::kPbPb_Run3_MC>(col, mcCols, tracks, d0sThisCol, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksD0sRun3PbPbMc, "Provide reconstructed and generated tracks and D0s in PbPb collisions", false);

  // ==========================================================================
  //  process functions - generator level only (for MCGEN datasets)
  // ==========================================================================

  // do not preslice the mcParticles table passed to fillMcParticle, otherwise the finding of the partonic mother can fail
  void processMcOnly(o2::aod::McCollisions const& mcCols, o2::aod::McParticles const& mcParticles)
  {
    mcBuilder.reset(mcCols, mcParticles);
    for (const auto& mcCol : mcCols) {
      auto particlesThisCollision = mcParticles.sliceBy(perMcCollision, mcCol.globalIndex());
      mcBuilder.fillMcCollision(mcCol, particlesThisCollision, mcProducts, pdgDb);
      for (const auto& mcParticle : particlesThisCollision) {
        mcBuilder.fillMcParticle<modes::System::kMC>(mcParticle, mcParticles, mcCol, mcProducts);
      }
    }
  }
  PROCESS_SWITCH(FemtoProducer, processMcOnly, "Provide generated particles", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{adaptAnalysisTask<FemtoProducer>(context)};
  return workflow;
}
