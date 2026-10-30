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
#include "PWGCF/Femto/Core/photonBuilder.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
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
#include "Common/DataModel/Qvectors.h"
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

using Run3Vzeros = o2::soa::Join<o2::aod::V0Datas, o2::aod::V0TOFNSigmas>;
using Run3RecoVzeros = o2::soa::Join<Run3Vzeros, o2::aod::McV0Labels>;

using Run3D0s = soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>;
using Run3RecoD0s = soa::Join<Run3D0s, aod::HfCand2ProngMcRec>;

using Run3Lcs = soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi>;
using Run3RecoLcs = soa::Join<Run3Lcs, aod::HfCand3ProngMcRec>;

using Run3Cascades = o2::soa::Join<o2::aod::CascDatas, o2::aod::CascTOFPIDs, o2::aod::CascTOFNSigmas>;
using Run3RecoCascades = o2::soa::Join<Run3Cascades, o2::aod::McCascLabels>;

using Run3Kinks = o2::aod::KinkCands;

using Run3McGenParticles = o2::aod::McParticles;

} // namespace o2::analysis::femto::rawinputs

struct FemtoProducer {

  // all Preslice declarations, grouped into a single data member (see o2::framework::PresliceGroup's
  // own doc comment) to stay well under the 100-direct-data-member cap on task structs
  struct : public o2::framework::PresliceGroup {
    o2::framework::Preslice<rawinputs::Run3McGenParticles> perMcCollision = o2::aod::mcparticle::mcCollisionId;
    o2::framework::Preslice<rawinputs::Run3McRecoTracks> perColRecoTracks = o2::aod::track::collisionId;
    o2::framework::Preslice<rawinputs::Run3RecoVzeros> perColRecoV0s = o2::aod::v0data::collisionId;
    o2::framework::Preslice<o2::aod::V0PhotonsKF> perColRecoPCMPhotons = o2::aod::v0photonkf::collisionId;
    o2::framework::Preslice<rawinputs::Run3RecoCascades> perColRecoCascades = o2::aod::cascdata::collisionId;
    o2::framework::Preslice<rawinputs::Run3Kinks> perColRecoKinks = o2::aod::kinkcand::collisionId;
    o2::framework::Preslice<rawinputs::Run3RecoD0s> perColRecoD0s = o2::aod::hf_cand::collisionId;
    o2::framework::Preslice<rawinputs::Run3RecoLcs> perColRecoLcs = o2::aod::hf_cand::collisionId;
  } preslices;

  // plain builder helper objects, grouped into a single data member (not Configurable/Partition/
  // Preslice/Produces, so the framework never needs to discover them individually -- free to nest)
  struct {
    collisionbuilder::CollisionBuilder<collisionbuilder::ColSelHistName, collisionbuilder::CollisionFilterHistName> collisionBuilder;
    trackbuilder::TrackBuilder<trackbuilder::TrackSelHistName, trackbuilder::TrackFilterHistName> trackBuilder;
    v0builder::V0Builder<modes::V0::kK0short, v0builder::K0shortSelHistName, v0builder::K0shortFilterHistName> k0shortBuilder;
    v0builder::V0Builder<modes::V0::kLambda, v0builder::LambdaSelHistName, v0builder::LambdaFilterHistName> lambdaBuilder;
    v0builder::V0Builder<modes::V0::kAntiLambda, v0builder::AntilambdaSelHistName, v0builder::AntiLambdaFilterHistName> antilambdaBuilder;
    photonbuilder::PhotonBuilder photonBuilder;
    charmhadronbuilder::CharmHadronBuilder<modes::CharmHadron::kD0, charmhadronbuilder::D0SelHistName, charmhadronbuilder::D0FilterHistName> d0Builder;
    charmhadronbuilder::CharmHadronBuilder<modes::CharmHadron::kD0Bar, charmhadronbuilder::D0barSelHistName, charmhadronbuilder::D0barFilterHistName> d0barBuilder;
    charmhadronbuilder::CharmHadronBuilder<modes::CharmHadron::kLc, charmhadronbuilder::LcSelHistName, charmhadronbuilder::LcFilterHistName> lcBuilder;
    charmhadronbuilder::CharmHadronBuilder<modes::CharmHadron::kLcBar, charmhadronbuilder::LcBarSelHistName, charmhadronbuilder::LcBarFilterHistName> lcBarBuilder;
    cascadebuilder::CascadeBuilder<modes::Cascade::kXi, cascadebuilder::XiSelHistName, cascadebuilder::XiFilterHistName> xiBuilder;
    cascadebuilder::CascadeBuilder<modes::Cascade::kOmega, cascadebuilder::OmegaSelHistName, cascadebuilder::OmegaFilterHistName> omegaBuilder;
    kinkbuilder::KinkBuilder<modes::Kink::kSigma, kinkbuilder::SigmaSelHistName, kinkbuilder::SigmaFilterHistName> sigmaBuilder;
    kinkbuilder::KinkBuilder<modes::Kink::kSigmaPlus, kinkbuilder::SigmaPlusSelHistName, kinkbuilder::SigmaPlusFilterHistName> sigmaPlusBuilder;
    mcbuilder::McBuilder mcBuilder;
  } builders;

  // ccdb config
  collisionbuilder::ConfCcdb confCcdb;

  // collision builder
  collisionbuilder::CollisionBuilderProducts collisionBuilderProducts;
  collisionbuilder::ConfCollisionTables confCollisionTables;
  collisionbuilder::ConfCollisionFilters confCollisionFilters;
  collisionbuilder::ConfCollisionBits confCollisionBits;
  collisionbuilder::ConfCollisionRctFlags confCollisionRctFlags;

  // track builder
  trackbuilder::TrackBuilderProducts trackBuilderProducts;
  trackbuilder::ConfTrackTables confTrackTables;
  trackbuilder::ConfTrackBits confTrackBits;
  trackbuilder::ConfTrackFilters confTrackFilters;

  // v0 builders
  v0builder::V0BuilderProducts v0builderProducts;
  v0builder::ConfV0Tables confV0Tables;
  v0builder::ConfV0Filters confV0Filters;
  v0builder::ConfK0shortBits confK0shortBits;
  v0builder::ConfLambdaBits confLambdaBits;

  // photon (PCM) builder
  photonbuilder::PhotonBuilderProducts photonBuilderProducts;
  photonbuilder::ConfPhotonTables confPhotonTables;
  photonbuilder::ConfPhotonFilters confPhotonFilters;
  photonbuilder::ConfPhotonBits confPhotonBits;

  // charm hadron builder
  charmhadronbuilder::CharmHadronBuilderProducts charmHadronBuilderProducts;
  charmhadronbuilder::ConfCharmHadronFilters confCharmHadronFilters;
  charmhadronbuilder::ConfCharmHadronTables confCharmHadronTables;
  charmhadronbuilder::ConfD0Bits confD0Bits;
  charmhadronbuilder::ConfLcBits confLcBits;

  // cascade builder
  cascadebuilder::CascadeBuilderProducts cascadeBuilderProducts;
  cascadebuilder::ConfCascadeTables confCascadeTables;
  cascadebuilder::ConfCascadeFilters confCascadeFilters;
  cascadebuilder::ConfXiBits confXiBits;
  cascadebuilder::ConfOmegaBits confOmegaBits;

  // kink builder
  kinkbuilder::KinkBuilderProducts kinkBuilderProducts;
  kinkbuilder::ConfKinkTables confKinkTables;
  kinkbuilder::ConfKinkFilters confKinkFilters;
  kinkbuilder::ConfSigmaBits confSigmaBits;
  kinkbuilder::ConfSigmaPlusBits confSigmaPlusBits;

  // mc builder
  mcbuilder::ConfMc confMc;
  mcbuilder::ConfMcTables confMcTables;
  mcbuilder::McBuilderProducts mcProducts;

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
    // builders.collisionBuilder.producingCollisions()/producingLiteCollisions()), and TrackBuilder
    // must init before the V0/Cascade/Kink builders (which check
    // builders.trackBuilder.producingTracks()/producingLiteTracks()).
    builders.collisionBuilder.init(&hRegistry, confCollisionFilters, confCollisionBits, confCollisionRctFlags, confCcdb, confCollisionTables, context);
    builders.trackBuilder.init(&hRegistry, confTrackBits, confTrackFilters, confTrackTables, context, builders.collisionBuilder);

    builders.k0shortBuilder.init(&hRegistry, confK0shortBits, confV0Filters, confV0Tables, context, builders.trackBuilder);
    builders.lambdaBuilder.init(&hRegistry, confLambdaBits, confV0Filters, confV0Tables, context, builders.trackBuilder);
    builders.antilambdaBuilder.init(&hRegistry, confLambdaBits, confV0Filters, confV0Tables, context, builders.trackBuilder);

    builders.photonBuilder.init(&hRegistry, confPhotonBits, confPhotonFilters, confPhotonTables, context);

    builders.d0Builder.init(&hRegistry, confD0Bits, confCharmHadronFilters, confCharmHadronTables, context);
    builders.d0barBuilder.init(&hRegistry, confD0Bits, confCharmHadronFilters, confCharmHadronTables, context);

    builders.lcBuilder.init(&hRegistry, confLcBits, confCharmHadronFilters, confCharmHadronTables, context);
    builders.lcBarBuilder.init(&hRegistry, confLcBits, confCharmHadronFilters, confCharmHadronTables, context);

    builders.sigmaBuilder.init(&hRegistry, confSigmaBits, confKinkFilters, confKinkTables, context, builders.trackBuilder);
    builders.sigmaPlusBuilder.init(&hRegistry, confSigmaPlusBits, confKinkFilters, confKinkTables, context, builders.trackBuilder);

    builders.xiBuilder.init(&hRegistry, confXiBits, confCascadeFilters, confCascadeTables, context, builders.trackBuilder);
    builders.omegaBuilder.init(&hRegistry, confOmegaBits, confCascadeFilters, confCascadeTables, context, builders.trackBuilder);

    builders.mcBuilder.init(confMc, confMcTables, context);

    // ---- guard: enabled tables vs. enabled process function ------------------
    if ((builders.xiBuilder.fillAnyTable() || builders.omegaBuilder.fillAnyTable()) &&
        (!doprocessTracksV0sCascadesRun3pp && !doprocessTracksV0sCascadesRun3PbPb &&
         !doprocessTracksV0sCascadesRun3PbPbWithEventShape &&
         !doprocessTracksV0sCascadesPhotonsRun3PbPbWithEventShape &&
         !doprocessTracksV0sCascadesKinksRun3pp && !doprocessTracksV0sCascadesKinksPhotonsRun3pp &&
         !doprocessTracksV0sCascadesRun3ppMc &&
         !doprocessTracksV0sCascadesRun3PbPbMc)) {
      LOG(fatal) << "At least one cascade table is enabled, but wrong process function is enabled. Breaking...";
    }
    if ((builders.lambdaBuilder.fillAnyTable() || builders.antilambdaBuilder.fillAnyTable() || builders.k0shortBuilder.fillAnyTable()) &&
        (!doprocessTracksV0sCascadesRun3pp && !doprocessTracksV0sCascadesRun3PbPb &&
         !doprocessTracksV0sRun3pp && !doprocessTracksV0sCascadesKinksRun3pp && !doprocessTracksV0sCascadesKinksPhotonsRun3pp &&
         !doprocessTracksV0sCascadesRun3PbPbWithEventShape && !doprocessTracksV0sCascadesPhotonsRun3PbPbWithEventShape &&
         !doprocessTracksV0sRun3ppMc && !doprocessTracksV0sRun3PbPb && !doprocessTracksV0sRun3PbPbMc &&
         !doprocessTracksV0sCascadesRun3ppMc && !doprocessTracksV0sCascadesRun3PbPbMc &&
         !doprocessTracksV0sKinksRun3ppMc)) {
      LOG(fatal) << "At least one v0 table is enabled, but wrong process function is enabled. Breaking...";
    }
    if (builders.photonBuilder.fillAnyTable() &&
        (!doprocessTracksPhotonsRun3pp && !doprocessTracksPhotonsRun3PbPb &&
         !doprocessTracksV0sCascadesKinksPhotonsRun3pp && !doprocessTracksV0sCascadesPhotonsRun3PbPbWithEventShape)) {
      LOG(fatal) << "At least one photon table is enabled, but wrong process function is enabled. Breaking...";
    }
    if ((builders.sigmaBuilder.fillAnyTable() || builders.sigmaPlusBuilder.fillAnyTable()) &&
        (!doprocessTracksKinksRun3pp && !doprocessTracksV0sCascadesKinksRun3pp && !doprocessTracksV0sCascadesKinksPhotonsRun3pp &&
         !doprocessTracksKinksRun3ppMc && !doprocessTracksV0sKinksRun3ppMc)) {
      LOG(fatal) << "At least one kink table is enabled, but wrong process function is enabled. Breaking...";
    }
    if ((builders.d0Builder.fillAnyTable() || builders.d0barBuilder.fillAnyTable()) &&
        (!doprocessTracksD0sRun3pp && !doprocessTracksD0sRun3PbPb &&
         !doprocessTracksD0sRun3ppMc && !doprocessTracksD0sRun3PbPbMc)) {
      LOG(fatal) << "At least one d0 table is enabled, but wrong process function is enabled. Breaking...";
    }
    if ((builders.lcBuilder.fillAnyTable() || builders.lcBarBuilder.fillAnyTable()) &&
        (!doprocessTracksLcsRun3pp && !doprocessTracksLcsRun3PbPb &&
         !doprocessTracksLcsRun3ppMc && !doprocessTracksLcsRun3PbPbMc)) {
      LOG(fatal) << "At least one lc table is enabled, but wrong process function is enabled. Breaking...";
    }
    if (builders.mcBuilder.fillAnyTable() &&
        (!doprocessTracksRun3ppMc && !doprocessTracksRun3PbPbMc &&
         !doprocessTracksV0sRun3ppMc && !doprocessTracksV0sRun3PbPbMc &&
         !doprocessTracksV0sCascadesRun3ppMc && !doprocessTracksV0sCascadesRun3PbPbMc &&
         !doprocessTracksKinksRun3ppMc && !doprocessTracksV0sKinksRun3ppMc &&
         !doprocessTracksD0sRun3ppMc && !doprocessTracksD0sRun3PbPbMc &&
         !doprocessTracksLcsRun3ppMc && !doprocessTracksLcsRun3PbPbMc &&
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
      add("collision", builders.collisionBuilder.fillAnyTable(), builders.collisionBuilder.isPassThrough());
      add("track", builders.trackBuilder.fillAnyTable(), builders.trackBuilder.isPassThrough());
      add("k0short", builders.k0shortBuilder.fillAnyTable(), builders.k0shortBuilder.isPassThrough());
      add("lambda", builders.lambdaBuilder.fillAnyTable(), builders.lambdaBuilder.isPassThrough());
      add("antilambda", builders.antilambdaBuilder.fillAnyTable(), builders.antilambdaBuilder.isPassThrough());
      add("photon", builders.photonBuilder.fillAnyTable(), builders.photonBuilder.isPassThrough());
      add("d0", builders.d0Builder.fillAnyTable(), builders.d0Builder.isPassThrough());
      add("d0bar", builders.d0barBuilder.fillAnyTable(), builders.d0barBuilder.isPassThrough());
      add("lc", builders.lcBuilder.fillAnyTable(), builders.lcBuilder.isPassThrough());
      add("lcbar", builders.lcBarBuilder.fillAnyTable(), builders.lcBarBuilder.isPassThrough());
      add("sigma", builders.sigmaBuilder.fillAnyTable(), builders.sigmaBuilder.isPassThrough());
      add("sigmaplus", builders.sigmaPlusBuilder.fillAnyTable(), builders.sigmaPlusBuilder.isPassThrough());
      add("xi", builders.xiBuilder.fillAnyTable(), builders.xiBuilder.isPassThrough());
      add("omega", builders.omegaBuilder.fillAnyTable(), builders.omegaBuilder.isPassThrough());
      add("mc", builders.mcBuilder.fillAnyTable(), builders.mcBuilder.isPassThrough());

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

      // in mc pass-through every generated collision is written, but reco collisions of other
      // sub-generators are still rejected in checkCollision -> biased event-loss denominator
      if (builders.mcBuilder.isPassThrough() && builders.collisionBuilder.isPassThrough() && builders.collisionBuilder.subGeneratorId() >= 0) {
        LOG(warn) << "MC and collision pass-through are enabled, but CollisionFilter.subGeneratorId = "
                  << builders.collisionBuilder.subGeneratorId() << " still rejects reconstructed collisions of other sub-generators, "
                  << "while all generated collisions are written. Set it to -1 for unbiased event/signal-loss corrections.";
      }

      // daughter tracks can be written a second time under the collision of their mother candidate;
      // without FTrackExtras (fillType) these rows cannot be told apart from selected tracks downstream
      const bool anyDaughterBuilder = builders.k0shortBuilder.fillAnyTable() || builders.lambdaBuilder.fillAnyTable() || builders.antilambdaBuilder.fillAnyTable() ||
                                      builders.xiBuilder.fillAnyTable() || builders.omegaBuilder.fillAnyTable() ||
                                      builders.sigmaBuilder.fillAnyTable() || builders.sigmaPlusBuilder.fillAnyTable() ||
                                      builders.d0Builder.fillAnyTable() || builders.d0barBuilder.fillAnyTable() || builders.lcBuilder.fillAnyTable() || builders.lcBarBuilder.fillAnyTable();
      if (builders.trackBuilder.isPassThrough() && anyDaughterBuilder && !builders.trackBuilder.producingTrackExtras()) {
        LOG(warn) << "Track pass-through is enabled together with candidates that register daughter tracks, but FTrackExtras is not produced. "
                  << "Daughter-only track rows cannot be identified downstream (enable TrackTables.produceTrackExtras).";
      }
    }

    // ---- guard: exactly one process function ---------------------------------
    const int nProcesses =
      static_cast<int>(doprocessTracksRun3pp) +
      static_cast<int>(doprocessTracksRun3PbPb) +
      static_cast<int>(doprocessTracksRun3PbPbWithEventShape) +
      static_cast<int>(doprocessTracksV0sRun3pp) +
      static_cast<int>(doprocessTracksPhotonsRun3pp) +
      static_cast<int>(doprocessTracksPhotonsRun3PbPb) +
      static_cast<int>(doprocessTracksV0sRun3PbPb) +
      static_cast<int>(doprocessTracksV0sCascadesRun3pp) +
      static_cast<int>(doprocessTracksV0sCascadesRun3PbPb) +
      static_cast<int>(doprocessTracksV0sCascadesRun3PbPbWithEventShape) +
      static_cast<int>(doprocessTracksV0sCascadesPhotonsRun3PbPbWithEventShape) +
      static_cast<int>(doprocessTracksKinksRun3pp) +
      static_cast<int>(doprocessTracksV0sCascadesKinksRun3pp) +
      static_cast<int>(doprocessTracksV0sCascadesKinksPhotonsRun3pp) +
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
      static_cast<int>(doprocessTracksLcsRun3pp) +
      static_cast<int>(doprocessTracksLcsRun3PbPb) +
      static_cast<int>(doprocessTracksLcsRun3ppMc) +
      static_cast<int>(doprocessTracksLcsRun3PbPbMc) +
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
    builders.collisionBuilder.reset();
    auto bc = col.template bc_as<T2>();
    builders.collisionBuilder.initCollision<system>(bc, col, tracks, ccdb, hRegistry);
    if (!builders.collisionBuilder.checkCollision(col)) {
      return false;
    }
    // Normally the collision row is written lazily, once a candidate is found. In
    // pass-through we need every accepted collision, including those with no candidate,
    // otherwise event-normalised observables lose exactly those events.
    if (builders.collisionBuilder.isPassThrough()) {
      builders.collisionBuilder.fillCollision<system>(collisionBuilderProducts, col);
    }
    return true;
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4>
  bool processMcCollisions(T1 const& col, T2 const& mcCols, T3 const& /* bcs*/, T4 const& tracks)
  {
    builders.collisionBuilder.reset();
    // builders.mcBuilder.reset() runs once per timeframe in the calling processTracksXXXMc function
    // (both pass-through and lazy-fill mode); resetting here on every reco collision would
    // wipe the dedup maps and duplicate every mc collision/particle referenced by more than
    // one reco collision (e.g. split-vertex reconstruction under pileup)
    auto bc = col.template bc_as<T3>();
    builders.collisionBuilder.initCollision<system>(bc, col, tracks, ccdb, hRegistry);
    if (!builders.collisionBuilder.checkCollision(col, mcCols)) {
      return false;
    }
    // NOTE: must be fillMcCollision, not fillCollision - fillMcCollision returns early on
    //       mCollisionAlreadyFilled, so filling the collision first would permanently skip
    //       the collision label
    if (builders.collisionBuilder.isPassThrough()) {
      builders.collisionBuilder.fillMcCollision<system>(collisionBuilderProducts, col, mcCols, mcProducts, builders.mcBuilder);
    }
    return true;
  }

  // processing tracks
  template <modes::System system, typename T1, typename T2>
  void processTracks(T1 const& col, T2 const& tracksWithItsPid)
  {
    builders.trackBuilder.reset(tracksWithItsPid);
    builders.trackBuilder.fillTracks<system>(col, builders.collisionBuilder, collisionBuilderProducts, tracksWithItsPid, trackBuilderProducts);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processMcTracks(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& tracksWithItsPid, T5 const& mcParticles)
  {
    builders.trackBuilder.reset(tracksWithItsPid);
    builders.trackBuilder.fillMcTracks<system>(col, builders.collisionBuilder, collisionBuilderProducts, mcCols, tracks, tracksWithItsPid, trackBuilderProducts, mcParticles, builders.mcBuilder, mcProducts);
  }

  // processing v0s
  template <modes::System system, typename T1, typename T2, typename T3>
  void processV0s(T1 const& col, T2 const& tracks, T3 const& v0s)
  {
    builders.lambdaBuilder.fillV0s<system>(col, builders.collisionBuilder, collisionBuilderProducts, trackBuilderProducts, v0builderProducts, v0s, tracks, builders.trackBuilder);
    builders.antilambdaBuilder.fillV0s<system>(col, builders.collisionBuilder, collisionBuilderProducts, trackBuilderProducts, v0builderProducts, v0s, tracks, builders.trackBuilder);
    builders.k0shortBuilder.fillV0s<system>(col, builders.collisionBuilder, collisionBuilderProducts, trackBuilderProducts, v0builderProducts, v0s, tracks, builders.trackBuilder);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processMcV0s(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& v0s, T5 const& mcParticles)
  {
    builders.lambdaBuilder.fillMcV0s<system>(col, builders.collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, v0builderProducts, v0s, tracks, builders.trackBuilder, mcParticles, builders.mcBuilder, mcProducts);
    builders.antilambdaBuilder.fillMcV0s<system>(col, builders.collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, v0builderProducts, v0s, tracks, builders.trackBuilder, mcParticles, builders.mcBuilder, mcProducts);
    builders.k0shortBuilder.fillMcV0s<system>(col, builders.collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, v0builderProducts, v0s, tracks, builders.trackBuilder, mcParticles, builders.mcBuilder, mcProducts);
  }

  // processing photons (PCM)
  template <modes::System system, typename T1, typename T2, typename T3>
  void processPhotons(T1 const& col, T2 const& photons, T3 const& v0legs)
  {
    builders.photonBuilder.fillPhotons<system>(col, builders.collisionBuilder, collisionBuilderProducts, photonBuilderProducts, photons, v0legs);
  }

  // processing kinks
  template <modes::System system, typename T1, typename T2, typename T3>
  void processKinks(T1 const& col, T2 const& tracks, T3 const& kinks)
  {
    builders.sigmaBuilder.fillKinks<system>(col, builders.collisionBuilder, collisionBuilderProducts, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, builders.trackBuilder);
    builders.sigmaPlusBuilder.fillKinks<system>(col, builders.collisionBuilder, collisionBuilderProducts, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, builders.trackBuilder);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processMcKinks(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& kinks, T5 const& mcParticles)
  {
    builders.sigmaBuilder.fillMcKinks<system>(col, builders.collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, builders.trackBuilder, mcParticles, builders.mcBuilder, mcProducts);
    builders.sigmaPlusBuilder.fillMcKinks<system>(col, builders.collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, kinkBuilderProducts, kinks, tracks, builders.trackBuilder, mcParticles, builders.mcBuilder, mcProducts);
  }

  // processing cascades
  template <modes::System system, typename T1, typename T2, typename T3>
  void processCascades(T1 const& col, T2 const& tracks, T3 const& cascades)
  {
    builders.xiBuilder.fillCascades<system>(col, builders.collisionBuilder, collisionBuilderProducts, trackBuilderProducts, cascadeBuilderProducts,
                                            cascades, tracks, builders.trackBuilder);
    builders.omegaBuilder.fillCascades<system>(col, builders.collisionBuilder, collisionBuilderProducts, trackBuilderProducts, cascadeBuilderProducts,
                                               cascades, tracks, builders.trackBuilder);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processMcCascades(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& cascades, T5 const& mcParticles)
  {
    builders.xiBuilder.fillMcCascades<system>(col, builders.collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, cascadeBuilderProducts, cascades, tracks, builders.trackBuilder, mcParticles, builders.mcBuilder, mcProducts);
    builders.omegaBuilder.fillMcCascades<system>(col, builders.collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, cascadeBuilderProducts, cascades, tracks, builders.trackBuilder, mcParticles, builders.mcBuilder, mcProducts);
  }

  // processing D0s
  template <modes::System system, typename T1, typename T2, typename T3>
  void processD0s(T1 const& col, T2 const& tracks, T3 const& candidates)
  {
    builders.d0Builder.fillCharmHadrons<system>(col, builders.collisionBuilder, collisionBuilderProducts, trackBuilderProducts, charmHadronBuilderProducts, candidates, tracks, builders.trackBuilder);
    builders.d0barBuilder.fillCharmHadrons<system>(col, builders.collisionBuilder, collisionBuilderProducts, trackBuilderProducts, charmHadronBuilderProducts, candidates, tracks, builders.trackBuilder);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processMcD0s(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& candidates, T5 const& mcParticles)
  {
    builders.d0Builder.fillMcCharmHadrons<system>(col, builders.collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, charmHadronBuilderProducts, candidates, tracks, builders.trackBuilder, mcParticles, builders.mcBuilder, mcProducts);
    builders.d0barBuilder.fillMcCharmHadrons<system>(col, builders.collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, charmHadronBuilderProducts, candidates, tracks, builders.trackBuilder, mcParticles, builders.mcBuilder, mcProducts);
  }

  // processing Lcs
  template <modes::System system, typename T1, typename T2, typename T3>
  void processLcs(T1 const& col, T2 const& tracks, T3 const& candidates)
  {
    builders.lcBuilder.fillCharmHadrons<system>(col, builders.collisionBuilder, collisionBuilderProducts, trackBuilderProducts, charmHadronBuilderProducts, candidates, tracks, builders.trackBuilder);
    builders.lcBarBuilder.fillCharmHadrons<system>(col, builders.collisionBuilder, collisionBuilderProducts, trackBuilderProducts, charmHadronBuilderProducts, candidates, tracks, builders.trackBuilder);
  }
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processMcLcs(T1 const& col, T2 const& mcCols, T3 const& tracks, T4 const& candidates, T5 const& mcParticles)
  {
    builders.lcBuilder.fillMcCharmHadrons<system>(col, builders.collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, charmHadronBuilderProducts, candidates, tracks, builders.trackBuilder, mcParticles, builders.mcBuilder, mcProducts);
    builders.lcBarBuilder.fillMcCharmHadrons<system>(col, builders.collisionBuilder, collisionBuilderProducts, mcCols, trackBuilderProducts, charmHadronBuilderProducts, candidates, tracks, builders.trackBuilder, mcParticles, builders.mcBuilder, mcProducts);
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

  void processTracksPhotonsRun3pp(rawinputs::Run3PpCollisions::iterator const& col,
                                  o2::aod::BCsWithTimestamps const& bcs,
                                  rawinputs::Run3FullPidTracks const& tracks,
                                  o2::aod::V0PhotonsKF const& photons,
                                  o2::aod::V0Legs const& v0legs)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    // V0PhotonsKF/V0Legs carry a plain CollisionId column (not a framework index column), so they
    // are not auto-grouped by the DPL binder and must be sliced manually, same as EMPhotonFilter.cxx
    auto photonsThisCol = photons.sliceBy(preslices.perColRecoPCMPhotons, col.globalIndex());
    processPhotons<modes::System::kPP_Run3>(col, photonsThisCol, v0legs);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksPhotonsRun3pp, "Provide tracks and PCM photons", false);

  void processTracksPhotonsRun3PbPb(rawinputs::Run3PbPbCollisions::iterator const& col,
                                    o2::aod::BCsWithTimestamps const& bcs,
                                    rawinputs::Run3FullPidTracks const& tracks,
                                    o2::aod::V0PhotonsKF const& photons,
                                    o2::aod::V0Legs const& v0legs)
  {
    if (!processCollisions<modes::System::kPbPb_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPbPb_Run3>(col, tracksWithItsPid);
    auto photonsThisCol = photons.sliceBy(preslices.perColRecoPCMPhotons, col.globalIndex());
    processPhotons<modes::System::kPbPb_Run3>(col, photonsThisCol, v0legs);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksPhotonsRun3PbPb, "Provide tracks and PCM photons in PbPb collisions", false);

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

  void processTracksV0sCascadesRun3PbPbWithEventShape(rawinputs::Run3PbPbCollisionsWithEventShape::iterator const& col,
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
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesRun3PbPbWithEventShape, "Provide tracks, v0s and cascades in PbPb collisions with event shape information", false);

  void processTracksV0sCascadesPhotonsRun3PbPbWithEventShape(rawinputs::Run3PbPbCollisionsWithEventShape::iterator const& col,
                                                             o2::aod::BCsWithTimestamps const& bcs,
                                                             rawinputs::Run3FullPidTracks const& tracks,
                                                             rawinputs::Run3Vzeros const& v0s,
                                                             rawinputs::Run3Cascades const& cascades,
                                                             o2::aod::V0PhotonsKF const& photons,
                                                             o2::aod::V0Legs const& v0legs)
  {
    if (!processCollisions<modes::System::kPbPb_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPbPb_Run3>(col, tracksWithItsPid);
    processV0s<modes::System::kPbPb_Run3>(col, tracks, v0s);
    processCascades<modes::System::kPbPb_Run3>(col, tracks, cascades);
    // V0PhotonsKF/V0Legs carry a plain CollisionId column (not a framework index column), so they
    // are not auto-grouped by the DPL binder and must be sliced manually, same as EMPhotonFilter.cxx
    auto photonsThisCol = photons.sliceBy(preslices.perColRecoPCMPhotons, col.globalIndex());
    processPhotons<modes::System::kPbPb_Run3>(col, photonsThisCol, v0legs);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesPhotonsRun3PbPbWithEventShape, "Provide tracks, v0s, cascades and PCM photons in PbPb collisions with event shape information", false);

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

  void processTracksV0sCascadesKinksPhotonsRun3pp(rawinputs::Run3PpCollisions::iterator const& col,
                                                  o2::aod::BCsWithTimestamps const& bcs,
                                                  rawinputs::Run3FullPidTracks const& tracks,
                                                  rawinputs::Run3Vzeros const& v0s,
                                                  rawinputs::Run3Cascades const& cascades,
                                                  rawinputs::Run3Kinks const& kinks,
                                                  o2::aod::V0PhotonsKF const& photons,
                                                  o2::aod::V0Legs const& v0legs)
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
    // V0PhotonsKF/V0Legs carry a plain CollisionId column (not a framework index column), so they
    // are not auto-grouped by the DPL binder and must be sliced manually, same as EMPhotonFilter.cxx
    auto photonsThisCol = photons.sliceBy(preslices.perColRecoPCMPhotons, col.globalIndex());
    processPhotons<modes::System::kPP_Run3>(col, photonsThisCol, v0legs);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksV0sCascadesKinksPhotonsRun3pp, "Provide tracks, v0s, cascades, kinks and PCM photons", false);

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
  //  These functions own the per-timeframe reset of the mc builder maps (once, before the
  //  reco-collision loop, in both pass-through and lazy-fill mode); processMcCollisions never
  //  resets them itself, since it runs once per reco collision and the dedup maps have to
  //  survive across all reco collisions in the timeframe.
  // ==========================================================================

  void processTracksRun3ppMc(rawinputs::Run3PpMcGenCollisions const& mcCols,
                             rawinputs::Run3PpMcRecoCollisions const& cols,
                             o2::aod::BCsWithTimestamps const& bcs,
                             rawinputs::Run3McRecoTracks const& tracks,
                             rawinputs::Run3McGenParticles const& mcParticles)
  {
    builders.mcBuilder.reset(mcCols, mcParticles);
    if (builders.mcBuilder.isPassThrough()) {
      builders.mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, preslices.perMcCollision, mcProducts);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(preslices.perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol)) {
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
    builders.mcBuilder.reset(mcCols, mcParticles);
    if (builders.mcBuilder.isPassThrough()) {
      builders.mcBuilder.fillMcPassThrough<modes::System::kPbPb_Run3_MC>(mcCols, mcParticles, preslices.perMcCollision, mcProducts);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(preslices.perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPbPb_Run3_MC>(col, mcCols, bcs, tracksThisCol)) {
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
    builders.mcBuilder.reset(mcCols, mcParticles);
    if (builders.mcBuilder.isPassThrough()) {
      builders.mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, preslices.perMcCollision, mcProducts);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(preslices.perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto v0sThisCol = v0s.sliceBy(preslices.perColRecoV0s, col.globalIndex());
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
    builders.mcBuilder.reset(mcCols, mcParticles);
    if (builders.mcBuilder.isPassThrough()) {
      builders.mcBuilder.fillMcPassThrough<modes::System::kPbPb_Run3_MC>(mcCols, mcParticles, preslices.perMcCollision, mcProducts);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(preslices.perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPbPb_Run3_MC>(col, mcCols, bcs, tracksThisCol)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPbPb_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto v0sThisCol = v0s.sliceBy(preslices.perColRecoV0s, col.globalIndex());
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
    builders.mcBuilder.reset(mcCols, mcParticles);
    if (builders.mcBuilder.isPassThrough()) {
      builders.mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, preslices.perMcCollision, mcProducts);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(preslices.perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto v0sThisCol = v0s.sliceBy(preslices.perColRecoV0s, col.globalIndex());
      processMcV0s<modes::System::kPP_Run3_MC>(col, mcCols, tracks, v0sThisCol, mcParticles);

      auto cascadesThisCol = cascades.sliceBy(preslices.perColRecoCascades, col.globalIndex());
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
    builders.mcBuilder.reset(mcCols, mcParticles);
    if (builders.mcBuilder.isPassThrough()) {
      builders.mcBuilder.fillMcPassThrough<modes::System::kPbPb_Run3_MC>(mcCols, mcParticles, preslices.perMcCollision, mcProducts);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(preslices.perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPbPb_Run3_MC>(col, mcCols, bcs, tracksThisCol)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPbPb_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto v0sThisCol = v0s.sliceBy(preslices.perColRecoV0s, col.globalIndex());
      processMcV0s<modes::System::kPbPb_Run3_MC>(col, mcCols, tracks, v0sThisCol, mcParticles);

      auto cascadesThisCol = cascades.sliceBy(preslices.perColRecoCascades, col.globalIndex());
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
    builders.mcBuilder.reset(mcCols, mcParticles);
    if (builders.mcBuilder.isPassThrough()) {
      builders.mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, preslices.perMcCollision, mcProducts);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(preslices.perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto kinksThisCol = kinks.sliceBy(preslices.perColRecoKinks, col.globalIndex());
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
    builders.mcBuilder.reset(mcCols, mcParticles);
    if (builders.mcBuilder.isPassThrough()) {
      builders.mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, preslices.perMcCollision, mcProducts);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(preslices.perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto v0sThisCol = v0s.sliceBy(preslices.perColRecoV0s, col.globalIndex());
      processMcV0s<modes::System::kPP_Run3_MC>(col, mcCols, tracks, v0sThisCol, mcParticles);

      auto kinksThisCol = kinks.sliceBy(preslices.perColRecoKinks, col.globalIndex());
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
    builders.mcBuilder.reset(mcCols, mcParticles);
    if (builders.mcBuilder.isPassThrough()) {
      builders.mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, preslices.perMcCollision, mcProducts);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(preslices.perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto d0sThisCol = d0s.sliceBy(preslices.perColRecoD0s, col.globalIndex());
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
    builders.mcBuilder.reset(mcCols, mcParticles);
    if (builders.mcBuilder.isPassThrough()) {
      builders.mcBuilder.fillMcPassThrough<modes::System::kPbPb_Run3_MC>(mcCols, mcParticles, preslices.perMcCollision, mcProducts);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(preslices.perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPbPb_Run3_MC>(col, mcCols, bcs, tracksThisCol)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPbPb_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto d0sThisCol = d0s.sliceBy(preslices.perColRecoD0s, col.globalIndex());
      processMcD0s<modes::System::kPbPb_Run3_MC>(col, mcCols, tracks, d0sThisCol, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksD0sRun3PbPbMc, "Provide reconstructed and generated tracks and D0s in PbPb collisions", false);

  void processTracksLcsRun3pp(rawinputs::Run3PpCollisions::iterator const& col,
                              o2::aod::BCsWithTimestamps const& bcs,
                              rawinputs::Run3FullPidTracks const& tracks,
                              rawinputs::Run3Lcs const& candidates)
  {
    if (!processCollisions<modes::System::kPP_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3>(col, tracksWithItsPid);
    processLcs<modes::System::kPP_Run3>(col, tracks, candidates);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksLcsRun3pp, "Provide tracks and Lcs", false);

  void processTracksLcsRun3PbPb(rawinputs::Run3PbPbCollisions::iterator const& col,
                                o2::aod::BCsWithTimestamps const& bcs,
                                rawinputs::Run3FullPidTracks const& tracks,
                                rawinputs::Run3Lcs const& candidates)
  {
    if (!processCollisions<modes::System::kPbPb_Run3>(col, bcs, tracks)) {
      return;
    }
    auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                            o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPbPb_Run3>(col, tracksWithItsPid);
    processLcs<modes::System::kPbPb_Run3>(col, tracks, candidates);
  }
  PROCESS_SWITCH(FemtoProducer, processTracksLcsRun3PbPb, "Provide tracks and Lcs in PbPb collisions", false);

  void processTracksLcsRun3ppMc(rawinputs::Run3PpMcGenCollisions const& mcCols,
                                rawinputs::Run3PpMcRecoCollisions const& cols,
                                o2::aod::BCsWithTimestamps const& bcs,
                                rawinputs::Run3McRecoTracks const& tracks,
                                rawinputs::Run3RecoLcs const& lcs,
                                rawinputs::Run3McGenParticles const& mcParticles)
  {
    builders.mcBuilder.reset(mcCols, mcParticles);
    if (builders.mcBuilder.isPassThrough()) {
      builders.mcBuilder.fillMcPassThrough<modes::System::kPP_Run3_MC>(mcCols, mcParticles, preslices.perMcCollision, mcProducts);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(preslices.perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPP_Run3_MC>(col, mcCols, bcs, tracksThisCol)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPP_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto lcsThisCol = lcs.sliceBy(preslices.perColRecoLcs, col.globalIndex());
      processMcLcs<modes::System::kPP_Run3_MC>(col, mcCols, tracks, lcsThisCol, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksLcsRun3ppMc, "Provide reconstructed and generated tracks and Lcs", false);

  void processTracksLcsRun3PbPbMc(rawinputs::Run3PbPbMcGenCollisions const& mcCols,
                                  rawinputs::Run3PbPbMcRecoCollisions const& cols,
                                  o2::aod::BCsWithTimestamps const& bcs,
                                  rawinputs::Run3McRecoTracks const& tracks,
                                  rawinputs::Run3RecoLcs const& lcs,
                                  rawinputs::Run3McGenParticles const& mcParticles)
  {
    builders.mcBuilder.reset(mcCols, mcParticles);
    if (builders.mcBuilder.isPassThrough()) {
      builders.mcBuilder.fillMcPassThrough<modes::System::kPbPb_Run3_MC>(mcCols, mcParticles, preslices.perMcCollision, mcProducts);
    }

    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(preslices.perColRecoTracks, col.globalIndex());
      if (!processMcCollisions<modes::System::kPbPb_Run3_MC>(col, mcCols, bcs, tracksThisCol)) {
        continue;
      }
      auto tracksWithItsPid = o2::soa::Attach<rawinputs::Run3McRecoTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      processMcTracks<modes::System::kPbPb_Run3_MC>(col, mcCols, tracksThisCol, tracksWithItsPid, mcParticles);

      auto lcsThisCol = lcs.sliceBy(preslices.perColRecoLcs, col.globalIndex());
      processMcLcs<modes::System::kPbPb_Run3_MC>(col, mcCols, tracks, lcsThisCol, mcParticles);
    }
  }
  PROCESS_SWITCH(FemtoProducer, processTracksLcsRun3PbPbMc, "Provide reconstructed and generated tracks and Lcs in PbPb collisions", false);

  // ==========================================================================
  //  process functions - generator level only (for MCGEN datasets)
  // ==========================================================================

  // do not preslice the mcParticles table passed to fillMcParticle, otherwise the finding of the partonic mother can fail
  void processMcOnly(o2::aod::McCollisions const& mcCols, o2::aod::McParticles const& mcParticles)
  {
    builders.mcBuilder.reset(mcCols, mcParticles);
    for (const auto& mcCol : mcCols) {
      auto particlesThisCollision = mcParticles.sliceBy(preslices.perMcCollision, mcCol.globalIndex());
      builders.mcBuilder.fillMcCollision(mcCol, particlesThisCollision, mcProducts, pdgDb);
      for (const auto& mcParticle : particlesThisCollision) {
        builders.mcBuilder.fillMcParticle<modes::System::kMC>(mcParticle, mcParticles, mcCol, mcProducts);
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
