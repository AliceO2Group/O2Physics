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

/// \file trackPropagationTester.cxx
/// \brief testing ground for track propagation
/// \author ALICE

//===============================================================
//
// Experimental version of the track propagation task
// this utilizes an analysis task module that can be employed elsewhere
// and allows for the re-utilization of a material LUT
//
// candidate approach for core service approach
//
//===============================================================

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"
#include "PWGLF/Utils/strangenessBuilderModule.h"
#include "Common/Tools/StandardCCDBLoader.h"

// The Run 3 AO2D stores the tracks at the point of innermost update. For a track with ITS this is the innermost (or second innermost)
// ITS layer. For a track without ITS, this is the TPC inner wall or for loopers in the TPC even a radius beyond that.
// In order to use the track parameters, the tracks have to be propagated to the collision vertex which is done by this task.
// The task consumes the TracksIU and TracksCovIU tables and produces Tracks and TracksCov to which then the user analysis can subscribe.
//
// This task is not needed for Run 2 converted data.
// There are two versions of the task (see process flags), one producing also the covariance matrix and the other only the tracks table.

using namespace o2;
using namespace o2::framework;
// using namespace o2::framework::expressions;

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

struct StrangenessBuilderTester {
  // CCDB boilerplate declarations
  o2::framework::Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // propagation stuff
  o2::common::StandardCCDBLoaderConfigurables standardCCDBLoaderConfigurables;
  o2::common::StandardCCDBLoader ccdbLoader;

  // boilerplate: strangeness builder stuff
  o2::pwglf::strangenessbuilder::products products;
  o2::pwglf::strangenessbuilder::coreConfigurables baseOpts;
  o2::pwglf::strangenessbuilder::v0Configurables v0BuilderOpts;
  o2::pwglf::strangenessbuilder::cascadeConfigurables cascadeBuilderOpts;
  o2::pwglf::strangenessbuilder::preSelectOpts preSelectOpts;
  o2::pwglf::strangenessbuilder::BuilderModule strangenessBuilderModule;

  // registry
  HistogramRegistry histos{"histos"};

  void init(o2::framework::InitContext& initContext)
  {
    // CCDB boilerplate init
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setURL(ccdburl.value);

    // task-specific
    strangenessBuilderModule.init(baseOpts, v0BuilderOpts, cascadeBuilderOpts, preSelectOpts, histos, initContext);
  }

  void processRealData(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    ccdbLoader.initCCDBfromBCs(standardCCDBLoaderConfigurables, ccdb, bcs);
    strangenessBuilderModule.dataProcess(ccdb, histos, collisions, static_cast<TObject*>(nullptr), v0s, cascades, trackedCascades, tracks, bcs, static_cast<TObject*>(nullptr), products);
  }

  void processMonteCarlo(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions, aod::McCollisions const& mccollisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtLabeledIU const& tracks, aod::BCsWithTimestamps const& bcs, aod::McParticles const& mcParticles)
  {
    ccdbLoader.initCCDBfromBCs(standardCCDBLoaderConfigurables, ccdb, bcs);
    strangenessBuilderModule.dataProcess(ccdb, histos, collisions, mccollisions, v0s, cascades, trackedCascades, tracks, bcs, mcParticles, products);
  }

  void processRealDataWithPID(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtIUWithPID const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    ccdbLoader.initCCDBfromBCs(standardCCDBLoaderConfigurables, ccdb, bcs);
    strangenessBuilderModule.dataProcess(ccdb, histos, collisions, static_cast<TObject*>(nullptr), v0s, cascades, trackedCascades, tracks, bcs, static_cast<TObject*>(nullptr), products);
  }

  void processMonteCarloWithPID(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions, aod::McCollisions const& mccollisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtLabeledIUWithPID const& tracks, aod::BCsWithTimestamps const& bcs, aod::McParticles const& mcParticles)
  {
    ccdbLoader.initCCDBfromBCs(standardCCDBLoaderConfigurables, ccdb, bcs);
    strangenessBuilderModule.dataProcess(ccdb, histos, collisions, mccollisions, v0s, cascades, trackedCascades, tracks, bcs, mcParticles, products);
  }

  PROCESS_SWITCH(StrangenessBuilderTester, processRealData, "process real data", true);
  PROCESS_SWITCH(StrangenessBuilderTester, processMonteCarlo, "process monte carlo", false);
  PROCESS_SWITCH(StrangenessBuilderTester, processRealDataWithPID, "process real data", false);
  PROCESS_SWITCH(StrangenessBuilderTester, processMonteCarloWithPID, "process monte carlo", false);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<StrangenessBuilderTester>(cfgc)};
  return workflow;
}
