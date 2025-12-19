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

/// \file propagationService.cxx
/// \brief
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

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Tools/StandardCCDBLoader.h"
#include "Common/Tools/TrackPropagationModule.h"
#include "Common/Tools/TrackTuner.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/GeomConstants.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include <string>

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

struct propagationServiceRun2 {
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
    ccdb->setFatalWhenNull(false);

    // task-specific
    strangenessBuilderModule.init(baseOpts, v0BuilderOpts, cascadeBuilderOpts, preSelectOpts, histos, initContext);
  }

  void processRealData(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, FullTracksExt const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    ccdbLoader.initCCDBfromBCs(standardCCDBLoaderConfigurables, ccdb, bcs);
    strangenessBuilderModule.dataProcess(ccdb, histos, collisions, static_cast<TObject*>(nullptr), v0s, cascades, static_cast<TObject*>(nullptr), tracks, bcs, static_cast<TObject*>(nullptr), products);
  }

  void processMonteCarlo(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions, aod::McCollisions const& mccollisions, aod::V0s const& v0s, aod::Cascades const& cascades, FullTracksExtLabeled const& tracks, aod::BCsWithTimestamps const& bcs, aod::McParticles const& mcParticles)
  {
    ccdbLoader.initCCDBfromBCs(standardCCDBLoaderConfigurables, ccdb, bcs);
    strangenessBuilderModule.dataProcess(ccdb, histos, collisions, mccollisions, v0s, cascades, static_cast<TObject*>(nullptr), tracks, bcs, mcParticles, products);
  }

  void processRealDataWithPID(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, FullTracksExtWithPID const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    ccdbLoader.initCCDBfromBCs(standardCCDBLoaderConfigurables, ccdb, bcs);
    strangenessBuilderModule.dataProcess(ccdb, histos, collisions, static_cast<TObject*>(nullptr), v0s, cascades, static_cast<TObject*>(nullptr), tracks, bcs, static_cast<TObject*>(nullptr), products);
  }

  void processMonteCarloWithPID(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions, aod::McCollisions const& mccollisions, aod::V0s const& v0s, aod::Cascades const& cascades, FullTracksExtLabeledWithPID const& tracks, aod::BCsWithTimestamps const& bcs, aod::McParticles const& mcParticles)
  {
    ccdbLoader.initCCDBfromBCs(standardCCDBLoaderConfigurables, ccdb, bcs);
    strangenessBuilderModule.dataProcess(ccdb, histos, collisions, mccollisions, v0s, cascades, static_cast<TObject*>(nullptr), tracks, bcs, mcParticles, products);
  }

  PROCESS_SWITCH(propagationServiceRun2, processRealData, "process real data", true);
  PROCESS_SWITCH(propagationServiceRun2, processMonteCarlo, "process monte carlo", false);
  PROCESS_SWITCH(propagationServiceRun2, processRealDataWithPID, "process real data", false);
  PROCESS_SWITCH(propagationServiceRun2, processMonteCarloWithPID, "process monte carlo", false);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<propagationServiceRun2>(cfgc)};
  return workflow;
}
