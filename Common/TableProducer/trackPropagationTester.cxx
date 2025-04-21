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

//===============================================================
//
// Experimental version of the track propagation task
// this utilizes an analysis task module that can be employed elsewhere
// and allows for the re-utilization of a material LUT
//
// candidate approach for core service approach
//
//===============================================================

#include "TableHelper.h"
#include "Common/Tools/TrackTuner.h"
#include "Common/Tools/TrackPropagationModule.h"
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

struct TrackPropagationTester {
  // produces group to be passed to track propagation module
  struct : ProducesGroup {
    Produces<aod::StoredTracks> tracksParPropagated;
    Produces<aod::TracksExtension> tracksParExtensionPropagated;
    Produces<aod::StoredTracksCov> tracksParCovPropagated;
    Produces<aod::TracksCovExtension> tracksParCovExtensionPropagated;
    Produces<aod::TracksDCA> tracksDCA;
    Produces<aod::TracksDCACov> tracksDCACov;
    Produces<aod::TrackTunerTable> tunertable;
  } trackPropagationProducts;

  // Configurables
  struct : ConfigurableGroup {
    std::string prefix = "ccdb";
    Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  } ccdbConfigurables;

  struct : ConfigurableGroup {
    std::string prefix = "trackPropagation";
    Configurable<float> minPropagationRadius{"minPropagationDistance", o2::constants::geom::XTPCInnerRef + 0.1, "Only tracks which are at a smaller radius will be propagated, defaults to TPC inner wall"};
    // for TrackTuner only (MC smearing)
    Configurable<bool> useTrackTuner{"useTrackTuner", false, "Apply track tuner corrections to MC"};
    Configurable<bool> useTrkPid{"useTrkPid", false, "use pid in tracking"};
    Configurable<bool> fillTrackTunerTable{"fillTrackTunerTable", false, "flag to fill track tuner table"};
    Configurable<int> trackTunerConfigSource{"trackTunerConfigSource", aod::track_tuner::InputString, "1: input string; 2: TrackTuner Configurables"};
    Configurable<std::string> trackTunerParams{"trackTunerParams", "debugInfo=0|updateTrackDCAs=1|updateTrackCovMat=1|updateCurvature=0|updateCurvatureIU=0|updatePulls=0|isInputFileFromCCDB=1|pathInputFile=Users/m/mfaggin/test/inputsTrackTuner/PbPb2022|nameInputFile=trackTuner_DataLHC22sPass5_McLHC22l1b2_run529397.root|pathFileQoverPt=Users/h/hsharma/qOverPtGraphs|nameFileQoverPt=D0sigma_Data_removal_itstps_MC_LHC22b1b.root|usePvRefitCorrections=0|qOverPtMC=-1.|qOverPtData=-1.", "TrackTuner parameter initialization (format: <name>=<value>|<name>=<value>)"};
    ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  } trackPropagationConfigurables;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  StandardCCDBLoader ccdbLoader;
  TrackPropagationModule trackPropagationMod;

  // registry
  HistogramRegistry registry{"registry"};

  void init(o2::framework::InitContext& initContext)
  {
    // configure ccdb
    ccdb->setURL(ccdbConfigurables.ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    ccdbLoader.readConfiguration(ccdbConfigurables);
    trackPropagationMod.init(initContext);
    trackPropagationMod.initHistograms(registry, trackPropagationConfigurables);
  }

  void processReal(soa::Join<aod::StoredTracksIU, aod::TracksCovIU, aod::TracksExtra> const& tracks, aod::Collisions const&, aod::BCs const& bcs)
  {
    ccdbLoader.initCCDBfromBCs(ccdb, bcs);
    trackPropagationMod.getFromCCDBLoader(ccdbLoader);
    trackPropagationMod.fillTrackTables<false>(tracks, trackPropagationProducts, registry);
  }
  PROCESS_SWITCH(TrackPropagationTester, processReal, "Process Real Data", true);

  // -----------------------
  void processMc(soa::Join<aod::StoredTracksIU, aod::McTrackLabels, aod::TracksCovIU, aod::TracksExtra> const& tracks, aod::McParticles const&, aod::Collisions const&, aod::BCs const& bcs)
  {
    ccdbLoader.initCCDBfromBCs(ccdb, bcs);
    trackPropagationMod.getFromCCDBLoader(ccdbLoader);
    trackPropagationMod.fillTrackTables<true>(tracks, trackPropagationProducts, registry);
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
 