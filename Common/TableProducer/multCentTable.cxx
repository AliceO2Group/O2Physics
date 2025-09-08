// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file multCentTable.cxx
/// \brief unified, self-configuring mult/cent provider
/// \author ALICE

//===============================================================
//
// Unified, self-configuring multiplicity+centrality task
// still work in progress: use at your own discretion
//
//===============================================================

#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/Core/MetadataHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Tools/Multiplicity/MultModule.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

o2::common::core::MetadataHelper metadataInfo; // Metadata helper

struct MultCentTable {
  o2::common::multiplicity::standardConfigurables opts;
  o2::common::multiplicity::products products;
  o2::common::multiplicity::MultModule module;

  // CCDB boilerplate declarations
  o2::framework::Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  // hold multiplicity values for layover to centrality calculation
  std::vector<o2::common::multiplicity::multEntry> mults;

  // slicers
  Preslice<soa::Join<aod::TracksIU, aod::TracksExtra>> slicerTracksIU = o2::aod::track::collisionId;
  Preslice<soa::Join<aod::TracksIU, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension>> slicerTracksIUwithSelections = o2::aod::track::collisionId;
  Preslice<soa::Join<aod::Tracks, aod::TracksExtra>> slicerTrackRun2 = o2::aod::track::collisionId;

  void init(o2::framework::InitContext& initContext)
  {
    // CCDB boilerplate init
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setURL(ccdburl.value);
    ccdb->setFatalWhenNull(false); // please never crash on your own, all exceptions captured (as they always should)

    // task-specific
    module.init(metadataInfo, opts, initContext);
  }

  void processRun2(soa::Join<aod::Collisions, aod::Run2MatchedSparse> const& collisions,
                   soa::Join<aod::Tracks, aod::TracksExtra> const& tracks,
                   soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps> const& bcs,
                   aod::Zdcs const&,
                   aod::FV0As const&,
                   aod::FV0Cs const&,
                   aod::FT0s const&)
  {
    mults.clear();
    for (auto const& collision : collisions) {
      o2::common::multiplicity::multEntry mult;
      const auto& bc = bcs.rawIteratorAt(collision.getId<aod::indices::BCId>());
      const uint64_t collIdx = collision.globalIndex();
      auto tracksThisCollision = tracks.sliceBy(slicerTrackRun2, collIdx);
      mult = module.collisionProcessRun2(collision, tracksThisCollision, bc, products);
      mults.push_back(mult);
    }
  }

  void processRun3(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                   soa::Join<aod::TracksIU, aod::TracksExtra> const& tracks,
                   soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse> const&,
                   aod::Zdcs const&,
                   aod::FV0As const&,
                   aod::FT0s const&,
                   aod::FDDs const&)
  {
    mults.clear();
    for (auto const& collision : collisions) {
      o2::common::multiplicity::multEntry mult;
      const auto& bc = collision.bc_as<soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>>();
      const uint64_t collIdx = collision.globalIndex();
      auto tracksThisCollision = tracks.sliceBy(slicerTracksIU, collIdx);
      mult = module.collisionProcessRun3(ccdb, metadataInfo, collision, tracksThisCollision, bc, products);
      mults.push_back(mult);
    }
  }

  void processRun3WithGlobalCounters(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                     soa::Join<aod::TracksIU, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension> const& tracks,
                                     soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse> const&,
                                     aod::Zdcs const&,
                                     aod::FV0As const&,
                                     aod::FT0s const&,
                                     aod::FDDs const&)
  {
    mults.clear();
    for (auto const& collision : collisions) {
      o2::common::multiplicity::multEntry mult;
      const auto& bc = collision.bc_as<soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>>();
      const uint64_t collIdx = collision.globalIndex();
      auto tracksThisCollision = tracks.sliceBy(slicerTracksIUwithSelections, collIdx);
      mult = module.collisionProcessRun3(ccdb, metadataInfo, collision, tracksThisCollision, bc, products);
      mults.push_back(mult);
    }
  }
  void processMFT(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                  o2::aod::MFTTracks const& mfttracks,
                  soa::SmallGroups<aod::BestCollisionsFwd> const& retracks)
  {
    if (opts.mEnabledTables[o2::common::multiplicity::kMFTMults]) {
      // populates MFT information in the mults buffer (in addition to filling table)
      module.collisionProcessMFT(collision, mfttracks, retracks, mults, products);
    }
  }
  void processMonteCarlo(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    if (opts.mEnabledTables[o2::common::multiplicity::kMultMCExtras]) {
      module.collisionProcessMonteCarlo(mcCollision, mcParticles, pdg, products);
    }
  }
  void processMonteCarlo2Mults(soa::Join<aod::McCollisionLabels, aod::Collisions>::iterator const& collision)
  {
    if (opts.mEnabledTables[o2::common::multiplicity::kMult2MCExtras]) {
      // establish simple interlink for posterior analysis (derived data)
      products.tableExtraMult2MCExtras(collision.mcCollisionId());
    }
  }
  void processCentralityRun2(aod::Collisions const& collisions, soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps> const& bcs)
  {
    // it is important that this function is at the end of the other process functions.
    // it requires `mults` to be properly set, which will only happen after the other process
    // functions have been called.

    // internally, the function below will do nothing if no centrality is requested.
    // it is thus safer to always keep the actual process function for centrality
    // generation to true, since the requisites for being in this context are
    // always fulfilled
    if (collisions.size() != static_cast<int64_t>(mults.size())) {
      LOGF(fatal, "Size of collisions doesn't match size of multiplicity buffer!");
    }
    module.generateCentralitiesRun2(ccdb, metadataInfo, bcs, mults, products);
  }
  void processCentralityRun3(aod::Collisions const& collisions, soa::Join<aod::BCs, aod::BcSels, aod::Timestamps> const& bcs, aod::FT0s const&)
  {
    // it is important that this function is at the end of the other process functions.
    // it requires `mults` to be properly set, which will only happen after the other process
    // functions have been called.

    // internally, the function below will do nothing if no centrality is requested.
    // it is thus safer to always keep the actual process function for centrality
    // generation to true, since the requisites for being in this context are
    // always fulfilled
    if (collisions.size() != static_cast<int64_t>(mults.size())) {
      LOGF(fatal, "Size of collisions doesn't match size of multiplicity buffer!");
    }
    module.generateCentralitiesRun3(ccdb, metadataInfo, bcs, mults, products);
  }

  PROCESS_SWITCH(MultCentTable, processRun2, "Process Run 2", false);
  PROCESS_SWITCH(MultCentTable, processRun3, "Process Run 3", true);
  PROCESS_SWITCH(MultCentTable, processRun3WithGlobalCounters, "Process Run 3 + global tracking counters", false);
  PROCESS_SWITCH(MultCentTable, processMFT, "Process MFT info", false);
  PROCESS_SWITCH(MultCentTable, processMonteCarlo, "Process Monte Carlo information", false);
  PROCESS_SWITCH(MultCentTable, processMonteCarlo2Mults, "Process Monte Carlo information", false);
  PROCESS_SWITCH(MultCentTable, processCentralityRun2, "Generate Run 2 centralities", false);
  PROCESS_SWITCH(MultCentTable, processCentralityRun3, "Generate Run 3 centralities", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  WorkflowSpec workflow{adaptAnalysisTask<MultCentTable>(cfgc)};
  return workflow;
}
