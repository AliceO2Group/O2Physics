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

// O2:
#include <CCDB/BasicCCDBManager.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/DataTypes.h"
#include "Common/DataModel/TrackSelectionTables.h" // needed for aod::TracksDCA table
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

using namespace o2;
using namespace o2::framework;

// *) Run 3:
using EventSelection = soa::Join<aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs>;
using CollisionRec = soa::Join<aod::Collisions, EventSelection>::iterator;
using CollisionRecSim = soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelection>::iterator;
using CollisionSim = aod::McCollision; // TBI 20240120 add support for centrality also for this case
using TracksRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
using TrackRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>::iterator;
using TracksRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
using TrackRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>::iterator;
using TracksSim = aod::McParticles;
using TrackSim = aod::McParticles::iterator;

// *) Converted Run 2:
using EventSelection_Run2 = soa::Join<aod::EvSels, aod::Mults, aod::CentRun2V0Ms, aod::CentRun2SPDTrks>;
using CollisionRec_Run2 = soa::Join<aod::Collisions, EventSelection_Run2>::iterator;
using CollisionRecSim_Run2 = soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelection_Run2>::iterator;
using TracksRecSim_Run2 = soa::Join<aod::Tracks, aod::McTrackLabels>;
using TrackRecSim_Run2 = soa::Join<aod::Tracks, aod::McTrackLabels>::iterator;

// *) Converted Run 1:
//    TBI 20240205 Since centrality calibration is not available in converted Run 1 data, I cannot treat it for the time being in the same way as converted Run 2.
//                 Once calibration is available, just use Run 2 above both for Run 2 and Run 1
using EventSelection_Run1 = soa::Join<aod::EvSels, aod::Mults>; // TBI 20240205 no calibration for centrality in converted LHC10h and LHC11h, at the time of writing
using CollisionRec_Run1 = soa::Join<aod::Collisions, EventSelection_Run1>::iterator;

// *) ROOT:
#include <TList.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGrid.h>
#include <Riostream.h>
#include <TRandom3.h>
#include <TComplex.h>
#include <TStopwatch.h>
#include <TF1.h>
#include <TF3.h>
using namespace std;

// *) Enums:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-Enums.h"

// *) Global constants:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-GlobalConstants.h"

// *) Main task:
struct MultiparticleCorrelationsAB // this name is used in lower-case format to name the TDirectoryFile in AnalysisResults.root
{

  // *) CCDB:
  Service<ccdb::BasicCCDBManager> ccdb;

// *) Configurables (cuts):
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-Configurables.h"

// *) Data members:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-DataMembers.h"

// *) Member functions:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-MemberFunctions.h"

  // -------------------------------------------

  // *) Initialize and book all objects:
  void init(InitContext const&)
  {
    // *) Trick to avoid name clashes, part 1;
    // *) Default configuration, booking, binning and cuts;
    // *) Insanity checks;
    // *) Book random generator;
    // *) Book base list;
    // *) Book all remaining objects;
    // ...
    // *) Trick to avoid name clashes, part 2;

    // *) Trick to avoid name clashes, part 1:
    Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    // *) Default configuration, booking, binning and cuts:
    DefaultConfiguration();
    DefaultBooking(); // here I decide only which histograms are booked, not details like binning, etc. That's done later in Book* member functions.
    DefaultBinning();
    DefaultCuts(); // Remark: has to be called after DefaultBinning(), since some default cuts are defined through default binning, to ease bookeeping

    // *) Set what to process - only rec, both rec and sim, only sim:
    // WhatToProcess(); // yes, this can be called here, after calling all Default* member functions above, because this has an effect only on Book* members functions, and the ones called afterward

    // *) Insanity checks:
    InsanityChecks();

    // *) Book random generator:
    delete gRandom;
    gRandom = new TRandom3(tc.fRandomSeed); // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID

    // *) Book base list:
    BookBaseList();

    // *) Book all remaining objects;
    BookAndNestAllLists();
    BookResultsHistograms(); // yes, this one has to be booked first, because it defines the commong binning for other groups of histograms
    BookEventHistograms();
    BookEventCutsHistograms();
    BookParticleHistograms();
    BookParticleCutsHistograms();
    BookQvectorHistograms();
    BookCorrelationsHistograms();
    BookWeightsHistograms();
    BookNestedLoopsHistograms();
    BookInternalValidationHistograms();
    BookTest0Histograms();
    BookTheRest(); // here I book everything that was not sorted (yet) in the specific functions above

    // *) Trick to avoid name clashes, part 2:
    TH1::AddDirectory(oldHistAddStatus);

  } // void init(InitContext const&)

  // -------------------------------------------

  // Since I am subscribing to different tables in each case, there are 3 separate implementations of process(...)
  // A) Process only reconstructed data;
  // B) Process both reconstructed and corresponding MC truth simulated data;
  // C) Process only simulated data.

  // For Run 2 converted data, I have the following implementations of process(...)
  // D) Process only converted reconstructed Run 2 data;
  // E) Process both converted reconstructed and corresponding MC truth simulated Run 2 data;
  // F) Process only converted simulated Run 2 data.

  // For Run 1 converted data, I have the following implementations of process(...)
  // G) Process only converted reconstructed Run 1 data;
  // H) Process both converted reconstructed and corresponding MC truth simulated Run 1 data;
  // I) Process only converted simulated Run 1 data.

  // For testing purposes I have processTest(...)
  // J) Process data with minimum subscription to the tables.

  // -------------------------------------------

  // A) Process only reconstructed data:
  void processRec(CollisionRec const& collision, aod::BCs const&, TracksRec const& tracks)
  {
    // Remark: Do not use here LOGF(fatal, ...) or LOGF(info, ...), because their stdout/stderr is suppressed. Use them in regular member functions instead.

    // *) Steer all analysis steps:
    Steer<eRec>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRec, "process only reconstructed data", true); // yes, keep always one process switch "true", so that I have default running version

  // -------------------------------------------

  // B) Process both reconstructed and corresponding MC truth simulated data:
  void processRecSim(CollisionRecSim const& collision, aod::BCs const&, TracksRecSim const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    Steer<eRecAndSim>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRecSim, "process both reconstructed and corresponding MC truth simulated data", false);

  // -------------------------------------------

  // C) Process only simulated data:
  void processSim(CollisionSim const& collision, aod::BCs const&, TracksSim const& tracks)
  {
    Steer<eSim>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processSim, "process only simulated data", false);

  // -------------------------------------------

  // D) Process only converted reconstructed Run 2 data:
  void processRec_Run2(CollisionRec_Run2 const& collision, aod::BCs const&, aod::Tracks const& tracks)
  {
    Steer<eRec_Run2>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRec_Run2, "process only converted reconstructed Run 2 data", false);

  // -------------------------------------------

  // E) Process both converted reconstructed and corresponding MC truth simulated Run 2 data:
  void processRecSim_Run2(CollisionRecSim_Run2 const& collision, aod::BCs const&, TracksRecSim_Run2 const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    Steer<eRecAndSim_Run2>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRecSim_Run2, "process both converted reconstructed and simulated Run 2 data", false);

  // -------------------------------------------

  // F) Process only converted simulated Run 2 data:
  void processSim_Run2(aod::Collision const&) // TBI 20240424 not ready yet, just a dummy to version to get later "doprocess..." variable.
  {
    // Steer<eSim_Run2>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processSim_Run2, "process converted only simulated Run 2 data", false);

  // -------------------------------------------

  // G) Process only converted reconstructed Run 1 data:
  void processRec_Run1(CollisionRec_Run1 const& collision, aod::BCs const&, aod::Tracks const& tracks)
  {
    Steer<eRec_Run1>(collision, tracks);
  }

  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRec_Run1, "process only converted reconstructed Run 1 data", false);

  // -------------------------------------------

  // H) Process both converted reconstructed and corresponding MC truth simulated Run 1 data;
  void processRecSim_Run1(aod::Collision const&) // TBI 20240424 not ready yet, just a dummy to version to get later "doprocess..." variable.
  {
    // Steer<eRecSim_Run1>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRecSim_Run1, "process both converted reconstructed and simulated Run 1 data", false);

  // -------------------------------------------

  // I) Process only converted simulated Run 1 data.
  void processSim_Run1(aod::Collision const&) // TBI 20240424 not ready yet, just a dummy to version to get later "doprocess..." variable.
  {
    // Steer<eSim_Run1>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processSim_Run1, "process only converted simulated Run 1 data", false);

  // -------------------------------------------

  // J) Process data with minimum subscription to the tables, for testing purposes:
  void processTest(aod::Collision const& collision, aod::BCs const&, aod::Tracks const& tracks)
  {
    Steer<eTest>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processTest, "test processing", false);

}; // struct MultiparticleCorrelationsAB

// -------------------------------------------

// *) ...
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiparticleCorrelationsAB>(cfgc),
  };
} // WorkflowSpec...
