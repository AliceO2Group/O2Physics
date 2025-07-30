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

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Tools/StandardCCDBLoader.h"
#include "Common/Tools/TrackPropagationModule.h"
#include "Common/Tools/TrackTuner.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/GeomConstants.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
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

struct TrackPropagationTester {
  o2::common::StandardCCDBLoaderConfigurables standardCCDBLoaderConfigurables;
  o2::common::TrackPropagationProducts trackPropagationProducts;
  o2::common::TrackPropagationConfigurables trackPropagationConfigurables;

  // the track tuner object -> needs to be here as it inherits from ConfigurableGroup (+ has its own copy of ccdbApi)
  TrackTuner trackTunerObj;

  // CCDB boilerplate declarations
  o2::framework::Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  o2::common::StandardCCDBLoader ccdbLoader;
  o2::common::TrackPropagationModule trackPropagation;

  HistogramRegistry registry{"registry"};

  void init(o2::framework::InitContext& initContext)
  {
    // CCDB boilerplate init
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setURL(ccdburl.value);

    // task-specific
    trackPropagation.init(trackPropagationConfigurables, trackTunerObj, registry, initContext);
  }

  void processReal(aod::Collisions const& collisions, soa::Join<aod::StoredTracksIU, aod::TracksCovIU, aod::TracksExtra> const& tracks, aod::Collisions const&, aod::BCs const& bcs)
  {
    // task-specific
    ccdbLoader.initCCDBfromBCs(standardCCDBLoaderConfigurables, ccdb, bcs);
    trackPropagation.fillTrackTables<false>(trackPropagationConfigurables, trackTunerObj, ccdbLoader, collisions, tracks, trackPropagationProducts, registry);
  }
  PROCESS_SWITCH(TrackPropagationTester, processReal, "Process Real Data", true);

  // -----------------------
  void processMc(aod::Collisions const& collisions, soa::Join<aod::StoredTracksIU, aod::McTrackLabels, aod::TracksCovIU, aod::TracksExtra> const& tracks, aod::McParticles const&, aod::Collisions const&, aod::BCs const& bcs)
  {
    ccdbLoader.initCCDBfromBCs(standardCCDBLoaderConfigurables, ccdb, bcs);
    trackPropagation.fillTrackTables<false>(trackPropagationConfigurables, trackTunerObj, ccdbLoader, collisions, tracks, trackPropagationProducts, registry);
  }
  PROCESS_SWITCH(TrackPropagationTester, processMc, "Process Monte Carlo", false);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<TrackPropagationTester>(cfgc)};
  return workflow;
}
