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
#include "Common/DataModel/TrackSelectionTables.h" // needed for aod::TracksDCA table
using namespace o2;
using namespace o2::framework;

using CollisionRec = aod::Collision;
using CollisionSim = aod::McCollision;

using TracksRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
using TrackRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>::iterator;

using TracksRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
using TrackRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>::iterator;

using TracksSim = aod::McParticles;
using TrackSim = aod::McParticles::iterator;

// ROOT:
#include "TList.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGrid.h"
#include "Riostream.h"
#include "TRandom3.h"
#include <TComplex.h>
using namespace std;

// *) Enums:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-Enums.h"

// *) Global constants:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-GlobalConstants.h"

// *) These are indended flags for PROCESS_SWITCH, have to be global, at least for the time being...
//    TBI 20231017 check this further, it doesn't work yet this way. It seems I have to pass to PROCESS_SWITCH(  ) only literals 'true' or 'false'
//    TBI 20231020 I could as well re-define them as data members...
bool gProcessRec = false;
bool gProcessRecSim = false;
bool gProcessSim = false;

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
    // *) Book base list;
    // *) Default configuration, booking, binning and cuts;
    // *) Configure the task with setters and getters;
    // *) Book all remaining objects;
    // ...
    // *) Trick to avoid name clashes, part 2;

    // *) Trick to avoid name clashes, part 1:
    Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    // *) Default configuration, booking, binning and cuts:
    DefaultConfiguration();
    DefaultBooking();
    DefaultBinning();
    DefaultCuts(); // Remark: has to be called after DefaultBinning(), since some default cuts are defined through default binning, to ease bookeeping

    // *) Set what to process - only rec, both rec and sim, only sim:
    WhatToProcess(); // yes, this can be called here, after calling all Default* member functions above, because this has an effect only on Book* members functions

    // *) Book random generator:
    delete gRandom;
    gRandom = new TRandom3(fRandomSeed); // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID

    // *) Book base list:
    BookBaseList();

    // *) Book all remaining objects;
    BookAndNestAllLists();
    BookEventHistograms();
    BookParticleHistograms();
    BookQvectorHistograms();
    BookCorrelationsHistograms();
    BookWeightsHistograms();
    BookNestedLoopsHistograms();
    BookTest0Histograms();
    BookResultsHistograms();

    // *) Trick to avoid name clashes, part 2:
    TH1::AddDirectory(oldHistAddStatus);

  } // void init(InitContext const&)

  // -------------------------------------------

  // Since I am subscribing to different tables in each case, there are 3 separate implementations of process(...)
  // A) Process only reconstructed data;
  // B) Process both reconstructed and simulated data;
  // C) Process only simulated data.

  // -------------------------------------------

  // A) Process only reconstructed data:
  void processRec(CollisionRec const& collision, aod::BCs const&, TracksRec const& tracks)
  {
    // ...

    // *) TBI 20231020 Temporary here (use configurable 'cfWhatToProcess' to set this flag corerectly):
    if (!gProcessRec) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
    }

    // *) If I reached max number of events, ignore the remaining collisions:
    if (!fProcessRemainingEvents) {
      return; // TBI 20231008 Temporarily implemented this way. But what I really need here is a graceful exit
              // from subsequent processing (which will also dump the output file, etc.). When that's possible,
              // move this to a member function Steer*(...)
    }

    // *) Steer all analysis steps:
    SteerRec(collision, tracks);
    // TBI 20231021 If I want to do some postprocessing after Steer(...), re-define Steer from void to bool, so that event cuts have effect, etc,

  } // void processRec(...)
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRec, "process only reconstructed information", false);

  // -------------------------------------------

  // B) Process both reconstructed and simulated data:
  void processRecSim(CollisionRec const& collision, aod::BCs const&, TracksRecSim const& tracks, aod::McParticles const&)
  {
    // ...

    // *) TBI 20231020 Temporary here (use configurable 'cfWhatToProcess' to set this flag corerectly):
    if (!gProcessRecSim) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
    }

    // *) If I reached max number of events, ignore the remaining collisions:
    if (!fProcessRemainingEvents) {
      return; // TBI 20231008 Temporarily implemented this way. But what I really need here is a graceful exit
              // from subsequent processing (which will also dump the output file, etc.). When that's possible,
              // move this to a member function Steer(...)
    }

    // *) Steer all analysis steps:
    SteerRecSim(collision, tracks);
    // TBI 20231021 If I want to do some postprocessing after Steer(...), re-define Steer from void to bool, so that event cuts have effect, etc,

  } // void processRecSim(...)
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRecSim, "process both reconstructed and simulated information", true);

  // -------------------------------------------

  // C) Process only simulated data:
  void processSim(CollisionSim const& collision, aod::BCs const&, TracksSim const& tracks)
  {
    // ...

    // *) If I reached max number of events, ignore the remaining collisions:
    if (!fProcessRemainingEvents) {
      return; // TBI 20231008 Temporarily implemented this way. But what I really need here is a graceful exit from subsequent processing (which will also dump the output file, etc.)
    }

    // *) TBI 20231020 Temporary here (use configurable 'cfWhatToProcess' to set this flag corerectly):
    if (!gProcessSim) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
    }

    // *) Steer all analysis steps:
    SteerSim(collision, tracks);
    // TBI 20231021 If I want to do some postprocessing after Steer(...), re-define Steer from void to bool, so that event cuts have effect, etc,

  } // void processSim(...)
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processSim, "process only simulated information", false);

}; // struct MultiparticleCorrelationsAB

// -------------------------------------------

// *) ...
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiparticleCorrelationsAB>(cfgc),
  };
} // WorkflowSpec...
