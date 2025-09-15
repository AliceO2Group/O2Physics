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
// Modularized version of TPC PID task
//
//===============================================================

#include "Common/Core/MetadataHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Tools/PID/pidTPCModule.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TObject.h>

#include <chrono>
#include <string>

using namespace o2;
using namespace o2::framework;

o2::common::core::MetadataHelper metadataInfo; // Metadata helper

struct pidTpcServiceRun2 {

  // CCDB boilerplate declarations
  o2::framework::Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  o2::aod::pid::pidTPCProducts products;
  o2::aod::pid::pidTPCConfigurables pidTPCopts;
  o2::aod::pid::pidTPCModule pidTPC;

  void init(o2::framework::InitContext& initContext)
  {
    // CCDB boilerplate init
    ccdb->setURL(ccdburl.value);
    ccdb->setFatalWhenNull(false); // manual fallback in case ccdb entry empty
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(ccdburl.value);

    // task-specific
    pidTPC.init(ccdb, ccdbApi, initContext, pidTPCopts, metadataInfo);
  }

  void processTracks(soa::Join<aod::Collisions, aod::EvSels> const& collisions, soa::Join<aod::Tracks, aod::TracksExtra> const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    pidTPC.process(ccdb, ccdbApi, bcs, collisions, tracks, static_cast<TObject*>(nullptr), products);
  }

  void processTracksMC(soa::Join<aod::Collisions, aod::EvSels> const& collisions, soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels> const& tracks, aod::BCsWithTimestamps const& bcs, aod::McParticles const&)
  {
    pidTPC.process(ccdb, ccdbApi, bcs, collisions, tracks, static_cast<TObject*>(nullptr), products);
  }

  PROCESS_SWITCH(pidTpcServiceRun2, processTracks, "Process Tracks", true);
  PROCESS_SWITCH(pidTpcServiceRun2, processTracksMC, "Process Tracks in MC (enables tune-on-data)", false);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata for later too
  metadataInfo.initMetadata(cfgc);

  WorkflowSpec workflow{adaptAnalysisTask<pidTpcServiceRun2>(cfgc)};
  return workflow;
}
