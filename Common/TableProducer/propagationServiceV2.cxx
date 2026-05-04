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

/// \file propagationServiceV2.cxx
/// \brief V2: GRPMagField and MeanVertexObject sourced from aod::GloCCDBObjects declarative CCDB table.
/// \author ALICE

//===============================================================
//
// Merged track propagation + strangeness building task
//
// Provides a common task to deal with track propagation and
// strangeness building in a single DPL device that is particularly
// adequate for pipelining.
//
//===============================================================

#include "PWGLF/Utils/strangenessBuilderModule.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/GloCCDBObjects.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/Tools/StandardCCDBLoader.h"
#include "Common/Tools/TrackPropagationModule.h"
#include "Common/Tools/TrackTuner.h"

#include <CCDB/BasicCCDBManager.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TObject.h>

#include <string>

using namespace o2;
using namespace o2::framework;

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using FullTracksExtWithPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTPCFullHe>;
using FullTracksExtIUWithPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTPCFullHe>;
using FullTracksExtLabeled = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::McTrackLabels>;
using FullTracksExtLabeledIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::McTrackLabels>;
using FullTracksExtLabeledWithPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTPCFullHe, aod::McTrackLabels>;
using FullTracksExtLabeledIUWithPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTPCFullHe, aod::McTrackLabels>;
using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;

// For dE/dx association in pre-selection
using TracksExtraWithPID = soa::Join<aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTPCFullHe>;

struct propagationServiceV2 {
  // Service<BasicCCDBManager> kept for MatLUT (rectifyPtrFromFile) and
  // strangenessBuilderModule (V-drift via ccdb->instance()).
  // GRPMagField and MeanVertex are sourced from CCDB columns instead.
  o2::framework::Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // propagation stuff — ccdbLoader used only for lut + mMeanVtx (set from column) + runNumber
  o2::common::StandardCCDBLoaderConfigurables standardCCDBLoaderConfigurables;
  o2::common::StandardCCDBLoader ccdbLoader;

  // Declarative CCDB path overrides (replace grpmagPath / mVtxPath in StandardCCDBLoaderConfigurables)
  // Option names: "ccdb:fGRPMagField" and "ccdb:fMeanVertex" respectively.
  ConfigurableCCDBPath<o2::aod::ccdbGlo::GRPMagField> grpmagPath;
  ConfigurableCCDBPath<o2::aod::ccdbGlo::MeanVertex> mVtxPath;

  // boilerplate: strangeness builder stuff
  o2::pwglf::strangenessbuilder::products products;
  o2::pwglf::strangenessbuilder::coreConfigurables baseOpts;
  o2::pwglf::strangenessbuilder::v0Configurables v0BuilderOpts;
  o2::pwglf::strangenessbuilder::cascadeConfigurables cascadeBuilderOpts;
  o2::pwglf::strangenessbuilder::preSelectOpts preSelectOpts;
  o2::pwglf::strangenessbuilder::BuilderModule strangenessBuilderModule;

  // the track tuner object -> needs to be here as it inherits from ConfigurableGroup (+ has its own copy of ccdbApi)
  TrackTuner trackTunerObj;

  // track propagation
  o2::common::TrackPropagationProducts trackPropagationProducts;
  o2::common::TrackPropagationConfigurables trackPropagationConfigurables;
  o2::common::TrackPropagationModule trackPropagation;

  using BCsWithCCDB = soa::Join<aod::BCsWithTimestamps, aod::GloCCDBObjects>;

  // registry
  HistogramRegistry histos{"histos"};

  void init(o2::framework::InitContext& initContext)
  {
    // Only needed for MatLUT fetch and strangenessBuilderModule V-drift
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setURL(ccdburl.value);

    // task-specific
    trackPropagation.init(trackPropagationConfigurables, trackTunerObj, histos, initContext);
    strangenessBuilderModule.init(baseOpts, v0BuilderOpts, cascadeBuilderOpts, preSelectOpts, histos, initContext);
  }

  // Load MatLUT once (needs rectifyPtrFromFile, kept manual), set B-field per run from
  // GRPMagField CCDB column, and refresh mMeanVtx pointer every call (pointer into current
  // BC table, valid only for the duration of this process() invocation).
  template <typename TBC>
  void initCCDB(TBC const& bc0)
  {
    if (!ccdbLoader.lut) {
      LOG(info) << "Loading material look-up table for run: " << bc0.runNumber();
      ccdbLoader.lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(
        ccdb->template getForRun<o2::base::MatLayerCylSet>(standardCCDBLoaderConfigurables.lutPath.value, bc0.runNumber()));
      o2::base::Propagator::Instance()->setMatLUT(ccdbLoader.lut);
    }
    // Always refresh: pointer into current BC table, invalidated after process() returns
    ccdbLoader.mMeanVtx = &bc0.meanVertex();
    if (ccdbLoader.runNumber != bc0.runNumber()) {
      const auto& grpmag = bc0.grpMagField(); // from declarative CCDB column
      LOG(info) << "Setting B-field to current " << grpmag.getL3Current() << " A for run " << bc0.runNumber() << " from GRPMagField CCDB column";
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      ccdbLoader.runNumber = bc0.runNumber();
    }
  }

  void processRealData(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtIU const& tracks, BCsWithCCDB const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());
    trackPropagation.fillTrackTables<false>(trackPropagationConfigurables, trackTunerObj, ccdbLoader, collisions, tracks, trackPropagationProducts, histos);
    strangenessBuilderModule.dataProcess(ccdb, histos, collisions, static_cast<TObject*>(nullptr), v0s, cascades, trackedCascades, tracks, bcs, static_cast<TObject*>(nullptr), products);
  }

  void processMonteCarlo(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions, aod::McCollisions const& mccollisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtLabeledIU const& tracks, BCsWithCCDB const& bcs, aod::McParticles const& mcParticles)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());
    trackPropagation.fillTrackTables<true>(trackPropagationConfigurables, trackTunerObj, ccdbLoader, collisions, tracks, trackPropagationProducts, histos);
    strangenessBuilderModule.dataProcess(ccdb, histos, collisions, mccollisions, v0s, cascades, trackedCascades, tracks, bcs, mcParticles, products);
  }

  void processRealDataWithPID(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtIUWithPID const& tracks, BCsWithCCDB const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());
    trackPropagation.fillTrackTables<false>(trackPropagationConfigurables, trackTunerObj, ccdbLoader, collisions, tracks, trackPropagationProducts, histos);
    strangenessBuilderModule.dataProcess(ccdb, histos, collisions, static_cast<TObject*>(nullptr), v0s, cascades, trackedCascades, tracks, bcs, static_cast<TObject*>(nullptr), products);
  }

  void processMonteCarloWithPID(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions, aod::McCollisions const& mccollisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtLabeledIUWithPID const& tracks, BCsWithCCDB const& bcs, aod::McParticles const& mcParticles)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());
    trackPropagation.fillTrackTables<true>(trackPropagationConfigurables, trackTunerObj, ccdbLoader, collisions, tracks, trackPropagationProducts, histos);
    strangenessBuilderModule.dataProcess(ccdb, histos, collisions, mccollisions, v0s, cascades, trackedCascades, tracks, bcs, mcParticles, products);
  }

  PROCESS_SWITCH(propagationServiceV2, processRealData, "process real data", true);
  PROCESS_SWITCH(propagationServiceV2, processMonteCarlo, "process monte carlo", false);
  PROCESS_SWITCH(propagationServiceV2, processRealDataWithPID, "process real data", false);
  PROCESS_SWITCH(propagationServiceV2, processMonteCarloWithPID, "process monte carlo", false);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<propagationServiceV2>(cfgc)};
  return workflow;
}
