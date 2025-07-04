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

/// \file multiparticle-correlations-ab.cxx
/// \brief ... TBI 20250425
/// \author Ante.Bilandzic@cern.ch

// O2:
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h" // needed for aod::TracksDCA table

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

using namespace o2;
using namespace o2::framework;

// *) Run 3:
using BCs_Run3 = soa::Join<aod::BCs, aod::Timestamps>; // TBI 20241126 under testing
// Remark 1: I have already timestamp in workflow, due to track-propagation.
// Remark 2: For consistency with notation below, drop _Run3 and instead use _Run2 and _Run1 TBI 20250401 not sure any longer what I wanted to say here...

// using EventSelection = soa::Join<aod::EvSels, aod::Mults, aod::MultsGlobal, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs>; // TBI 20241209 validating "MultsGlobal"
//  for using collision.multNTracksGlobal() TBI 20250128 do i still need this?
using EventSelection = soa::Join<aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs>;
// TBI 20250128 I can't join here directly aod::CentNGlobals, see email from DDC from 20250127 if this one requires a special treatment
//              See in https://github.com/AliceO2Group/O2Physics/blob/master/Common/DataModel/Centrality.h how centrality tables are named exactly
using CollisionRec = soa::Join<aod::Collisions, EventSelection>::iterator; // use in json "isMC": "true" for "event-selection-task"
using CollisionRecSim = soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelection>::iterator;
// using CollisionRecSim = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::MultsExtraMC, EventSelection>::iterator; // TBI 20241210 validating "MultsExtraMC" for multMCNParticlesEta08
using CollisionSim = aod::McCollision;
using TracksRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
// using TrackRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>::iterator;
using TracksRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>; // + use in json "isMC" : "true"
using TrackRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>::iterator;
using TracksSim = aod::McParticles;
using TrackSim = aod::McParticles::iterator;

// *) Run 2:
using EventSelection_Run2 = soa::Join<aod::EvSels, aod::Mults, aod::CentRun2V0Ms, aod::CentRun2SPDTrks>; // TBI 20240517 do not subscribe to CentRun2CL0s and CentRun2CL1s => see enum
using CollisionRec_Run2 = soa::Join<aod::Collisions, EventSelection_Run2>::iterator;                     // use in json "isRun2MC" : "true" for "event-selection-task"
using CollisionRecSim_Run2 = soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelection_Run2>::iterator;
// Remark: For tracks, I can use everything same as in Run 3

// *) Run 1:
//    TBI 20240205 Since centrality calibration is not available in converted Run 1 data, I cannot treat it for the time being in the same way as converted Run 2.
//                 Once calibration is available, just use Run 2 above both for Run 2 and Run 1
using EventSelection_Run1 = soa::Join<aod::EvSels, aod::Mults>; // TBI 20240205 no calibration for centrality in converted LHC10h and LHC11h, at the time of writing
using CollisionRec_Run1 = soa::Join<aod::Collisions, EventSelection_Run1>::iterator;
using CollisionRecSim_Run1 = soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelection_Run1>::iterator;
// Remark: For tracks, I can use everything same as in Run 3

// *) QA:
//    Remark: This is Run 3 "Rec" + subscription to additional few tables (otherwise unnecessary in my analysis, e.g. some specific detector tables), used only for QA purposes.
//            Therefore, I start all definitions from what I have defined for Run 3 "Rec", and on top of it join these additional tables for QA.
using BCs_QA = soa::Join<BCs_Run3, aod::BcSels, aod::Run3MatchedToBCSparse>;
//             *) BcSels => bc.has_foundFT0(), etc.
//             *) Run3MatchedToBCSparse => bc.has_zdc(), etc. TBI 20250401 at the moment, I do not use this one
using Collision_QA = CollisionRec; // if I would need additional tables for QA, just join 'em here with CollisionRec
using TracksRec_QA = TracksRec;    // if I would need additional tables for QA, just join 'em here with TracksRec

// *) ROOT:
#include <TComplex.h>
#include <TDatabasePDG.h>
#include <TExMap.h>
#include <TF1.h>
#include <TF3.h>
#include <TFile.h>
#include <TFormula.h>
#include <TGrid.h>
#include <TH1D.h>
#include <THnSparse.h>
#include <TList.h>
#include <TObjString.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TSystem.h>

#include <Riostream.h>

#include <complex>
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
  ctpRateFetcher mRateFetcher; // see email from MP on 20240508 and example usage in O2Physics/PWGLF/TableProducer/Common/zdcSP.cxx

  // *) O2DatabsePDG service shared service between different tasks (do not use TDatabasePDG directly, because it's not shared)
  //    See also Tutorials/src/usingPDGCervice.cxx
  //    TBI 20250625 enable the line below and switch to O2DatabsePDG when memory consumption with O2DatabsePDG is resolved, and then replace in all functions
  //      tc.fDatabasePDG->GetParticle(track.pdgCode()) with pdg->GetParticle(track.pdgCode()) + same for mcParticle + remove TDatabasePDG.h
  // Service<o2::framework::O2DatabasePDG> pdg;

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
    // *) Insanity checks before booking;
    // *) Book random generator;
    // *) Book base list;
    // *) Book all remaining objects;
    // *) Insanity checks after booking;
    // *) Trick to avoid name clashes, part 2;

    // *) Trick to avoid name clashes, part 1:
    bool oldHistAddStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    // *) Default configuration, booking, binning and cuts:
    insanityChecksOnDefinitionsOfConfigurables(); // values passed via configurables are insanitized here. Nothing is initialized yet via configurables in this method
    defaultConfiguration();                       // here default values from configurables are taken into account
    defaultBooking();                             // here I decide only which histograms are booked, not details like binning, etc.
    defaultBinning();                             // here default values for bins are either hardwired, or values for bins provided via configurables are taken into account
    defaultCuts();                                // here default values for cuts are either hardwired, or defined through default binning to ease bookeeping,
                                                  // or values for cuts provided via configurables are taken into account
                                                  // Remark: defaultCuts() has to be called after defaultBinning()

    // *) Specific cuts:
    if (tc.fUseSpecificCuts) {
      specificCuts(tc.fWhichSpecificCuts); // after default cuts are applied, on top of them apply analysis-specific cuts. Has to be called after defaultBinning() and defaultCuts()
    }

    // *) Insanity checks before booking:
    insanityChecksBeforeBooking(); // check only hardwired values and the ones obtained from configurables

    // *) Book random generator:
    delete gRandom;
    gRandom = new TRandom3(tc.fRandomSeed); // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID

    // *) Book base list:
    bookBaseList();

    // *) Book all remaining objects;
    bookAndNestAllLists();
    bookResultsHistograms(); // yes, this one has to be booked first, because it defines the common binning for other groups of histograms, w/ or w/o clonning
    bookQAHistograms();
    bookEventHistograms();
    bookEventCutsHistograms();
    bookParticleHistograms();
    bookParticleCutsHistograms(); // memStatus: 50913
    bookQvectorHistograms();      // memStatus: 50913 (without differential q-vectors and eta separations)
    bookCorrelationsHistograms();
    bookWeightsHistograms();
    bookCentralityWeightsHistograms();
    bookNestedLoopsHistograms();
    bookNUAHistograms();
    bookInternalValidationHistograms();
    bookTest0Histograms();
    bookEtaSeparationsHistograms();
    bookTheRest(); // I book everything that was not sorted (yet) in the specific functions above
    // memStatus: 50913 (without differential q-vectors and eta separations)

    // *) I can purge a few objects used for common consistent booking across different groups of histograms:
    purgeAfterBooking();

    // *) Insanity checks after booking:
    insanityChecksAfterBooking(); // pointers of all local histograms, etc., are available, so I can do insanity checks directly on all booked objects

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
  void processRec(CollisionRec const& collision, BCs_Run3 const& bcs, TracksRec const& tracks)
  {
    // Remark: Do not use here LOGF(fatal, ...) or LOGF(info, ...), because their stdout/stderr is suppressed. Use them in regular member functions instead.

    // *) Steer all analysis steps:
    Steer<eRec>(collision, bcs, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRec, "process only reconstructed data", true); // yes, keep always one process switch "true", so that I have default running version

  // -------------------------------------------

  // B) Process both reconstructed and corresponding MC truth simulated data:
  void processRecSim(CollisionRecSim const& collision, aod::BCs const& bcs, TracksRecSim const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    Steer<eRecAndSim>(collision, bcs, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRecSim, "process both reconstructed and corresponding MC truth simulated data", false);

  // -------------------------------------------

  // C) Process only simulated data:
  void processSim(CollisionSim const& /*collision*/, aod::BCs const& /*bcs*/, TracksSim const& /*tracks*/)
  {
    //    Steer<eSim>(collision, bcs, tracks); // TBI 20240517 not ready yet, but I do not really need this one urgently, since RecSim is working, and I need that one for efficiencies...
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processSim, "process only simulated data", false);

  // -------------------------------------------

  // D) Process only converted reconstructed Run 2 data:
  void processRec_Run2(CollisionRec_Run2 const& collision, aod::BCs const& bcs, TracksRec const& tracks)
  {
    Steer<eRec_Run2>(collision, bcs, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRec_Run2, "process only converted reconstructed Run 2 data", false);

  // -------------------------------------------

  // E) Process both converted reconstructed and corresponding MC truth simulated Run 2 data:
  void processRecSim_Run2(CollisionRecSim_Run2 const& collision, aod::BCs const& bcs, TracksRecSim const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    Steer<eRecAndSim_Run2>(collision, bcs, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRecSim_Run2, "process both converted reconstructed and simulated Run 2 data", false);

  // -------------------------------------------

  // F) Process only converted simulated Run 2 data:
  void processSim_Run2(CollisionSim const& /*collision*/) // TBI 20240517 extend this subscription eventually
  {
    // Steer<eSim_Run2>(collision, tracks); // TBI 20240517 not ready yet, but I do not really need this one urgently, since RecSim_Run2 is working, and I need that one for efficiencies...
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processSim_Run2, "process only converted simulated Run 2 data", false);

  // -------------------------------------------

  // G) Process only converted reconstructed Run 1 data:
  void processRec_Run1(CollisionRec_Run1 const& collision, aod::BCs const& bcs, TracksRec const& tracks)
  {
    Steer<eRec_Run1>(collision, bcs, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRec_Run1, "process only converted reconstructed Run 1 data", false);

  // -------------------------------------------

  // H) Process both converted reconstructed and corresponding MC truth simulated Run 1 data;
  void processRecSim_Run1(CollisionRecSim_Run1 const& /*collision*/, aod::BCs const& /*bcs*/, TracksRecSim const& /*tracks*/, aod::McParticles const&, aod::McCollisions const&)
  {
    // Steer<eRecAndSim_Run1>(collision, bcs, tracks); // TBI 20240517 not ready yet, but for benchmarking in any case I need only "Rec"
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processRecSim_Run1, "process both converted reconstructed and simulated Run 1 data", false);

  // -------------------------------------------

  // I) Process only converted simulated Run 1 data.
  void processSim_Run1(aod::Collision const&) // TBI 20240424 not ready yet, just a dummy to version to get later "doprocess..." variable.
  {
    // Steer<eSim_Run1>(collision, tracks); // TBI 20240517 not ready yet, but for benchmarking in any case I need only "Rec"
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processSim_Run1, "process only converted simulated Run 1 data", false);

  // -------------------------------------------

  // J) Process data with minimum subscription to the tables, for testing purposes:
  //    Remark: To keep this branch as simple as possible, I do not subscribe to centrality table. Therefore, when running with "processTest": "true" in JSON,
  //            I have to remove "| o2-analysis-centrality-table $JsonFile \" from workflow (yes, remove, not comment out!)
  void processTest(aod::Collision const& collision, aod::BCs const& bcs, aod::Tracks const& tracks)
  {
    Steer<eTest>(collision, bcs, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processTest, "test processing", false);

  // -------------------------------------------

  // K) Process data with more than necessary subscriptions to the tables, only for QA purposes:
  //    Remark 1: This is basically the main "processRec" switch, merely enhanced with subscription to few more tables (e.g. detector specific), only for QA purposes.
  //    Remark 2: Ideally, i use the same workflow for "processRec" and "processQA", but most likely at some point I will have to establish separate workflow for "processQA"
  void processQA(Collision_QA const& collision, BCs_QA const& bcs, TracksRec_QA const& tracks, aod::FT0s const&)
  {
    // Summary for additional tables subscribed to directly here:
    // *) FT0s => bc.foundFT0().sumAmpC(), etc.
    Steer<eQA>(collision, bcs, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsAB, processQA, "QA processing", false);

  // -------------------------------------------

  // L) Process extra Monte Carlo info the from table HepMCHeavyIons:
  //    Remark 1: Under testing, merge eventually this process switch with processRecSim, processRecSim_Run2, and processRecSim_Run1 above;
  //              This switch does everything the same as processRecSim (so it works only for Run 3 at the moment), except extra info is processed from table HepMCHeavyIons with dedicated function call
  //              ProcessHepMCHeavyIons(hepMChi). I use this dedicated function, in order not to modify call to Steer(...) by adding the fourth argument.
  //              TBI 20250429 see if I can circumvent this with templates (i can NOT join HepMCHeavyIons and McParticles), in order to keep call to Steer(...) as simple as it is now
  //    Remark 2: In MC LHC24g3 and LHC24e2c, for HepMCHeavyIons only impact parameter is filled;
  //              As soon as HepMCHeavyIons is correctly filled in MC productions, merge this switch with processRecSim, processRecSim_Run2, and processRecSim_Run1 above.
  //              Most notably, I will need hep.centrality() (centrality at generated level), and hep.sigmaInelNN()
  void processHepMChi(CollisionRecSim const& collision, aod::BCs const& bcs, TracksRecSim const& tracks, aod::HepMCHeavyIons const& hepMChi, aod::McParticles const&, aod::McCollisions const&)
  {
    // Comment the weather here...

    // *) Check if this collision has the corresponding MC collision:
    if (!collision.has_mcCollision()) {
      LOGF(warning, "\033[1;31m%s at line %d : No MC collision for this collision, skip... \033[0m", __FUNCTION__, __LINE__);
      return;
    }

    // *) For this collision, get the corresponding mcCollision, and then profit from the fact that HepMCHeavyIons have index to mcCollision by default (no need to join with McCollisionLabels):
    auto hep = hepMChi.iteratorAt(collision.mcCollision().globalIndex());

    // *) Quick insanity check that HepMCHeavyIons and McCollisions refer to the same MC collision:
    //    Since both of them provide getter impactParameter(), i simply check if it gives the same value in both cases:
    if (std::abs(hep.impactParameter() - collision.mcCollision().impactParameter()) > tc.fFloatingPointPrecision) {
      LOGF(fatal, "\033[1;31m%s at line %d : impactParameter accessed from HepMCHeavyIons and McCollisions is not the same, they do not correspond to the same MC event \033[0m", __FUNCTION__, __LINE__);
    }

    // *) Okay, extract all extra info from HepMCHeavyIons:
    ProcessHepMCHeavyIons(hep);

    // *) Call the Steer(...)
    //    TBI 20250429 For the time being, only Run 3 call for Steer(...) is supported. When generalizing to Run 2 and Run 1 process switches, perhaps the better strategy is
    //                 just to inject ProcessHepMCHeavyIons(hep); , and keep call to Steer(...) as it is now?
    Steer<eRecAndSim>(collision, bcs, tracks); // TBI 20250429 remember that I have hardwired here eRecAndSim, so this now works only for Run 3

  } // void processHepMChi( ... )

  PROCESS_SWITCH(MultiparticleCorrelationsAB, processHepMChi, "HepMCHeavyIons processing", false);

  // -------------------------------------------

  // ... ctd. here with further process switches ...

  // -------------------------------------------

}; // struct MultiparticleCorrelationsAB

// -------------------------------------------

// *) ...
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiparticleCorrelationsAB>(cfgc),
  };
} // WorkflowSpec...
