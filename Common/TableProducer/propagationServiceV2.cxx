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
/// \brief V2: GRPMagField, MeanVertexObject and the material LUT sourced from declarative CCDB tables.
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
#include "Common/DataModel/TpcCCDBObjects.h"
#include "Common/DataModel/TrackTunerCCDBObjects.h"
#include "Common/Tools/TrackPropagationModule.h"
#include "Common/Tools/TrackTuner.h"

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
  // No CCDB client of any kind: GRPMagField, MeanVertex, the material LUT and the TPC
  // drift correction all arrive as declarative CCDB columns.
  // NB: TrackTuner still holds its own CcdbApi for the DCA calibration files when
  // useTrackTuner is enabled; its path is derived from the run number, which a CCDB
  // column cannot express today.

  // Keep every CCDB path user-overridable, as ccdb.lutPath / ccdb.grpmagPath /
  // ccdb.mVtxPath were before the migration. The option name and default are derived
  // from the column itself ("ccdb:fMatLUT" and friends), so they cannot drift from the
  // declaration the way the old per-task Configurables did. Declaring them is enough —
  // the accessors stay bc.matLUT() / bc.grpMagField() / bc.meanVertex().
  o2::framework::ConfigurableCCDBPath<aod::ccdbGlo::MatLUT> matLUTPath;
  o2::framework::ConfigurableCCDBPath<aod::ccdbGlo::GRPMagField> grpMagFieldPath;
  o2::framework::ConfigurableCCDBPath<aod::ccdbGlo::MeanVertex> meanVertexPath;

  // Everything this task needs is read straight off the CCDB columns at the point of
  // use; no StandardCCDBLoader, and no CCDB query of its own. The single piece of
  // retained state is the run number, needed only to avoid re-installing the magnetic
  // field: that one is not a lookup but global state in the Propagator /
  // TGeoGlobalMagField singletons, and installing it rebuilds or rescales the field map.
  int mRunNumber = -1;

  // boilerplate: strangeness builder stuff
  o2::pwglf::strangenessbuilder::products products;
  o2::pwglf::strangenessbuilder::coreConfigurables baseOpts;
  o2::pwglf::strangenessbuilder::v0Configurables v0BuilderOpts;
  o2::pwglf::strangenessbuilder::cascadeConfigurables cascadeBuilderOpts;
  o2::pwglf::strangenessbuilder::preSelectOpts preSelectOpts;
  o2::pwglf::strangenessbuilder::eventSelectOpts eventSelectOpts;
  o2::pwglf::strangenessbuilder::BuilderModule strangenessBuilderModule;

  // the track tuner object -> needs to be here as it inherits from ConfigurableGroup (+ has its own copy of ccdbApi)
  TrackTuner trackTunerObj;

  // track propagation
  o2::common::TrackPropagationProducts trackPropagationProducts;
  o2::common::TrackPropagationConfigurables trackPropagationConfigurables;
  o2::common::TrackPropagationModule trackPropagation;

  using BCsWithCCDB = soa::Join<aod::BCsWithTimestamps, aod::GloCCDBObjects, aod::GeomCCDBObjects, aod::TpcCalibCCDBObjects, aod::TrackTunerCCDBObjects>;

  // registry
  HistogramRegistry histos{"histos"};

  void init(o2::framework::InitContext& initContext)
  {
    // task-specific
    trackPropagation.init(trackPropagationConfigurables, trackTunerObj, histos, initContext, /*calibFromCCDBColumns=*/true);
    strangenessBuilderModule.init(baseOpts, v0BuilderOpts, cascadeBuilderOpts, preSelectOpts, eventSelectOpts, histos, initContext);
  }

  /// Install into the Propagator the two things which are global state rather than
  /// values: the magnetic field and the material LUT.
  template <typename TBC>
  void initPropagator(TBC const& bc0)
  {
    if (mRunNumber != bc0.runNumber()) {
      LOG(info) << "Setting B-field to current " << bc0.grpMagField().getL3Current() << " A for run " << bc0.runNumber() << " from GRPMagField CCDB column";
      o2::base::Propagator::initFieldFromGRP(&bc0.grpMagField());
      mRunNumber = bc0.runNumber();
    }
    // A pointer store, so it costs nothing to redo every timeframe — and doing so
    // means a relocated column buffer is picked up for free instead of dangling.
    // The column's finaliser has already run MatLayerCylSet::rectifyPtrFromFile.
    o2::base::Propagator::Instance()->setMatLUT(&bc0.matLUT());
  }

  void processRealData(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtIU const& tracks, BCsWithCCDB const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    auto bc0 = bcs.begin();
    initPropagator(bc0);
    trackPropagation.fillTrackTables<false>(trackPropagationConfigurables, trackTunerObj, bc0.runNumber(), &bc0.meanVertex(), &bc0.trackTunerDca(), &bc0.trackTunerQOverPt(), collisions, tracks, trackPropagationProducts, histos);
    strangenessBuilderModule.dataProcess(histos, collisions, static_cast<TObject*>(nullptr), v0s, cascades, trackedCascades, tracks, bcs, static_cast<TObject*>(nullptr), products);
  }

  void processMonteCarlo(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions, aod::McCollisions const& mccollisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtLabeledIU const& tracks, BCsWithCCDB const& bcs, aod::McParticles const& mcParticles)
  {
    if (bcs.size() == 0) {
      return;
    }
    auto bc0 = bcs.begin();
    initPropagator(bc0);
    trackPropagation.fillTrackTables<true>(trackPropagationConfigurables, trackTunerObj, bc0.runNumber(), &bc0.meanVertex(), &bc0.trackTunerDca(), &bc0.trackTunerQOverPt(), collisions, tracks, trackPropagationProducts, histos);
    strangenessBuilderModule.dataProcess(histos, collisions, mccollisions, v0s, cascades, trackedCascades, tracks, bcs, mcParticles, products);
  }

  void processRealDataWithPID(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtIUWithPID const& tracks, BCsWithCCDB const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    auto bc0 = bcs.begin();
    initPropagator(bc0);
    trackPropagation.fillTrackTables<false>(trackPropagationConfigurables, trackTunerObj, bc0.runNumber(), &bc0.meanVertex(), &bc0.trackTunerDca(), &bc0.trackTunerQOverPt(), collisions, tracks, trackPropagationProducts, histos);
    strangenessBuilderModule.dataProcess(histos, collisions, static_cast<TObject*>(nullptr), v0s, cascades, trackedCascades, tracks, bcs, static_cast<TObject*>(nullptr), products);
  }

  void processMonteCarloWithPID(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions, aod::McCollisions const& mccollisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtLabeledIUWithPID const& tracks, BCsWithCCDB const& bcs, aod::McParticles const& mcParticles)
  {
    if (bcs.size() == 0) {
      return;
    }
    auto bc0 = bcs.begin();
    initPropagator(bc0);
    trackPropagation.fillTrackTables<true>(trackPropagationConfigurables, trackTunerObj, bc0.runNumber(), &bc0.meanVertex(), &bc0.trackTunerDca(), &bc0.trackTunerQOverPt(), collisions, tracks, trackPropagationProducts, histos);
    strangenessBuilderModule.dataProcess(histos, collisions, mccollisions, v0s, cascades, trackedCascades, tracks, bcs, mcParticles, products);
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
