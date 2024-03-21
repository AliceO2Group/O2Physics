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

#include <cmath>
#include <algorithm>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TableHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "PWGCF/DataModel/DptDptFiltered.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <TList.h>
#include <TDirectory.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile3D.h>

#include "PWGCF/TableProducer/dptdptfilter.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis;

#define DPTDPTFILTERLOGCOLLISIONS debug
#define DPTDPTFILTERLOGTRACKS debug

namespace o2::analysis::dptdptfilter
{
using DptDptFullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using DptDptFullTracksAmbiguous = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::TrackCompColls>;
using DptDptTracksPID = soa::Join<aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFbeta>;
using DptDptFullTracksPID = soa::Join<DptDptFullTracks, DptDptTracksPID>;
using DptDptFullTracksPIDAmbiguous = soa::Join<DptDptFullTracksAmbiguous, DptDptTracksPID>;
using DptDptFullTracksDetLevel = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using DptDptFullTracksDetLevelAmbiguous = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::TrackCompColls>;
using DptDptFullTracksPIDDetLevel = soa::Join<DptDptFullTracksDetLevel, DptDptTracksPID>;
using DptDptFullTracksPIDDetLevelAmbiguous = soa::Join<DptDptFullTracksDetLevelAmbiguous, DptDptTracksPID>;

bool fullDerivedData = false; /* produce full derived data for its external storage */

/// \enum MatchRecoGenSpecies
/// \brief The species considered by the matching test
enum MatchRecoGenSpecies {
  kDptDptCharged = 0, ///< charged particle/track
  kDptDptElectron,    ///< electron
  kDptDptMuon,        ///< muon
  kDptDptPion,        ///< pion
  kDptDptKaon,        ///< kaon
  kDptDptProton,      ///< proton
  kDptDptNoOfSpecies, ///< the number of considered species
  kWrongSpecies = -1
};

const char* speciesName[kDptDptNoOfSpecies] = {"h", "e", "mu", "pi", "ka", "p"};

const char* speciesTitle[kDptDptNoOfSpecies] = {"", "e", "#mu", "#pi", "K", "p"};

//============================================================================================
// The DptDptFilter histogram objects
// TODO: consider registering in the histogram registry
//============================================================================================
TH1F* fhCentMultB = nullptr;
TH1F* fhCentMultA = nullptr;
TH1F* fhVertexZB = nullptr;
TH1F* fhVertexZA = nullptr;
TH1F* fhMultB = nullptr;
TH1F* fhMultA = nullptr;
TH1F* fhPB = nullptr;
std::vector<TH1F*> fhPA;
TH1F* fhPtB = nullptr;
std::vector<TH1F*> fhPtA;
TH1F* fhPtPosB = nullptr;
std::vector<TH1F*> fhPtPosA;
TH1F* fhPtNegB = nullptr;
std::vector<TH1F*> fhPtNegA;
std::vector<TH2F*> fhNPosNegA;
std::vector<TH1F*> fhDeltaNA;

TH1F* fhEtaB = nullptr;
TH1F* fhEtaA = nullptr;

TH1F* fhPhiB = nullptr;
TH1F* fhPhiA = nullptr;

TH1F* fhDCAxyB = nullptr;
TH1F* fhDCAxyA = nullptr;
TH1F* fhFineDCAxyA = nullptr;
TH1F* fhDCAzB = nullptr;
TH1F* fhDCAzA = nullptr;
TH1F* fhFineDCAzA = nullptr;

TH2D* fhAmbiguousTrackType = nullptr;
TH2F* fhAmbiguousTrackPt = nullptr;
TH2F* fhAmbiguityDegree = nullptr;
TH2F* fhCompatibleCollisionsZVtxRms = nullptr;

TH1F* fhTrueCentMultB = nullptr;
TH1F* fhTrueCentMultA = nullptr;
TH1F* fhTrueVertexZB = nullptr;
TH1F* fhTrueVertexZA = nullptr;
TH1F* fhTrueVertexZAA = nullptr;
TH1F* fhTruePB = nullptr;
std::vector<TH1F*> fhTruePA;
TH1F* fhTruePtB = nullptr;
std::vector<TH1F*> fhTruePtA;
TH1F* fhTruePtPosB = nullptr;
std::vector<TH1F*> fhTruePtPosA;
TH1F* fhTruePtNegB = nullptr;
std::vector<TH1F*> fhTruePtNegA;
std::vector<TH2F*> fhTrueNPosNegA;
std::vector<TH1F*> fhTrueDeltaNA;

TH1F* fhTrueEtaB = nullptr;
TH1F* fhTrueEtaA = nullptr;

TH1F* fhTruePhiB = nullptr;
TH1F* fhTruePhiA = nullptr;

TH1F* fhTrueDCAxyB = nullptr;
TH1F* fhTrueDCAxyA = nullptr;
TH1F* fhTrueDCAzB = nullptr;
TH1F* fhTrueDCAxyBid = nullptr;
TH1F* fhTrueDCAzA = nullptr;

//============================================================================================
// The DptDptFilter multiplicity counters
//============================================================================================
std::vector<int> trkMultPos;  // multiplicity of positive tracks
std::vector<int> trkMultNeg;  // multiplicity of negative tracks
std::vector<int> partMultPos; // multiplicity of positive particles
std::vector<int> partMultNeg; // multiplicity of negative particles
} // namespace o2::analysis::dptdptfilter

using namespace dptdptfilter;

struct DptDptFilter {
  struct : ConfigurableGroup {
    Configurable<std::string> cfgCCDBUrl{"input_ccdburl", "http://ccdb-test.cern.ch:8080", "The CCDB url for the input file"};
    Configurable<std::string> cfgCCDBPathName{"input_ccdbpath", "", "The CCDB path for the input file. Default \"\", i.e. don't load from CCDB"};
    Configurable<std::string> cfgCCDBDate{"input_ccdbdate", "20220307", "The CCDB date for the input file"};
    Configurable<std::string> cfgCCDBPeriod{"input_ccdbperiod", "LHC22o", "The CCDB dataset period for the input file"};
  } cfginputfile;
  Configurable<bool> cfgFullDerivedData{"fullderiveddata", false, "Produce the full derived data for external storage. Default false"};
  Configurable<std::string> cfgCentMultEstimator{"centmultestimator", "V0M", "Centrality/multiplicity estimator detector: V0M,CL0,CL1,FV0A,FT0M,FT0A,FT0C,NTPV,NOCM: none. Default V0M"};
  Configurable<std::string> cfgSystem{"syst", "PbPb", "System: pp, PbPb, Pbp, pPb, XeXe, ppRun3, PbPbRun3. Default PbPb"};
  Configurable<std::string> cfgDataType{"datatype", "data", "Data type: data, datanoevsel, MC, FastMC, OnTheFlyMC. Default data"};
  Configurable<std::string> cfgTriggSel{"triggsel", "MB", "Trigger selection: MB, None. Default MB"};
  Configurable<std::string> cfgCentSpec{"centralities", "00-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80", "Centrality/multiplicity ranges in min-max separated by commas"};
  Configurable<o2::analysis::DptDptBinningCuts> cfgBinning{"binning",
                                                           {28, -7.0, 7.0, 18, 0.2, 2.0, 16, -0.8, 0.8, 72, 0.5},
                                                           "triplets - nbins, min, max - for z_vtx, pT, eta and phi, binning plus bin fraction of phi origin shift"};
  Configurable<bool> cfgTraceCollId0{"tracecollid0", false, "Trace particles in collisions id 0. Default false"};

  OutputObj<TList> fOutput{"DptDptFilterCollisionsInfo", OutputObjHandlingPolicy::AnalysisObject};

  Produces<aod::DptDptCFAcceptedCollisions> acceptedcollisions;
  Produces<aod::DptDptCFCollisionsInfo> collisionsinfo;
  Produces<aod::DptDptCFAcceptedTrueCollisions> acceptedtrueevents;
  Produces<aod::DptDptCFGenCollisionsInfo> gencollisionsinfo;

  void init(InitContext const&)
  {
    using namespace dptdptfilter;

    LOGF(info, "DptDptFilterTask::init()");

    fullDerivedData = cfgFullDerivedData;

    /* update with the configurable values */
    /* the binning */
    ptbins = cfgBinning->mPTbins;
    ptlow = cfgBinning->mPTmin;
    ptup = cfgBinning->mPTmax;
    etabins = cfgBinning->mEtabins;
    etalow = cfgBinning->mEtamin;
    etaup = cfgBinning->mEtamax;
    zvtxbins = cfgBinning->mZVtxbins;
    zvtxlow = cfgBinning->mZVtxmin;
    zvtxup = cfgBinning->mZVtxmax;
    /* the track types and combinations */
    /* the centrality/multiplicity estimation */
    fCentMultEstimator = getCentMultEstimator(cfgCentMultEstimator);
    /* the trigger selection */
    fTriggerSelection = getTriggerSelection(cfgTriggSel);
    traceCollId0 = cfgTraceCollId0;

    /* if the system type is not known at this time, we have to put the initialization somewhere else */
    fSystem = getSystemType(cfgSystem);
    fDataType = getDataType(cfgDataType);

    /* create the output list which will own the task histograms */
    TList* fOutputList = new TList();
    fOutputList->SetOwner(true);
    fOutput.setObject(fOutputList);

    if ((fDataType == kData) || (fDataType == kDataNoEvtSel) || (fDataType == kMC)) {
      /* create the reconstructed data histograms */
      /* TODO: proper axes and axes titles according to the system; still incomplete */
      std::string multestimator = getCentMultEstimatorName(fCentMultEstimator);
      if (fSystem > kPbp) {
        fhCentMultB = new TH1F("CentralityB", "Centrality before cut; centrality (%)", 100, 0, 100);
        fhCentMultA = new TH1F("CentralityA", "Centrality; centrality (%)", 100, 0, 100);
        fhMultB = new TH1F("MultB", TString::Format("%s Multiplicity before cut;%s Multiplicity;Collisions", multestimator.c_str(), multestimator.c_str()), 4001, -0.5, 4000.5);
        fhMultA = new TH1F("MultA", TString::Format("%s Multiplicity;%s Multiplicity;Collisions", multestimator.c_str(), multestimator.c_str()), 4001, -0.5, 4000.5);
      } else {
        /* for pp, pPb and Pbp systems use multiplicity instead */
        fhCentMultB = new TH1F("MultiplicityB", "Multiplicity before cut; multiplicity (%)", 100, 0, 100);
        fhCentMultA = new TH1F("MultiplicityA", "Multiplicity; multiplicity (%)", 100, 0, 100);
        fhMultB = new TH1F("MultB", TString::Format("%s Multiplicity before cut;%s Multiplicity;Collisions", multestimator.c_str(), multestimator.c_str()), 601, -0.5, 600.5);
        fhMultA = new TH1F("MultA", TString::Format("%s Multiplicity;%s Multiplicity;Collisions", multestimator.c_str(), multestimator.c_str()), 601, -0.5, 600.5);
      }

      fhVertexZB = new TH1F("VertexZB", "Vertex Z; z_{vtx}", 60, -15, 15);
      fhVertexZA = new TH1F("VertexZA", "Vertex Z; z_{vtx}", zvtxbins, zvtxlow, zvtxup);

      /* add the hstograms to the output list */
      fOutputList->Add(fhCentMultB);
      fOutputList->Add(fhCentMultA);
      fOutputList->Add(fhMultB);
      fOutputList->Add(fhMultA);
      fOutputList->Add(fhVertexZB);
      fOutputList->Add(fhVertexZA);
    }

    if ((fDataType != kData) && (fDataType != kDataNoEvtSel)) {
      /* create the true data histograms */
      /* TODO: proper axes and axes titles according to the system; still incomplete */
      if (fSystem > kPbp) {
        fhTrueCentMultB = new TH1F("TrueCentralityB", "Centrality before (truth); centrality (%)", 100, 0, 100);
        fhTrueCentMultA = new TH1F("TrueCentralityA", "Centrality (truth); centrality (%)", 100, 0, 100);
      } else {
        /* for pp, pPb and Pbp systems use multiplicity instead */
        fhTrueCentMultB = new TH1F("TrueMultiplicityB", "Multiplicity before (truth); multiplicity (%)", 100, 0, 100);
        fhTrueCentMultA = new TH1F("TrueMultiplicityA", "Multiplicity (truth); multiplicity (%)", 100, 0, 100);
      }

      fhTrueVertexZB = new TH1F("TrueVertexZB", "Vertex Z before (truth); z_{vtx}", 60, -15, 15);
      fhTrueVertexZA = new TH1F("TrueVertexZA", "Vertex Z (truth); z_{vtx}", zvtxbins, zvtxlow, zvtxup);
      fhTrueVertexZAA = new TH1F("TrueVertexZAA", "Vertex Z (truth rec associated); z_{vtx}", zvtxbins, zvtxlow, zvtxup);

      /* add the hstograms to the output list */
      fOutputList->Add(fhTrueCentMultB);
      fOutputList->Add(fhTrueCentMultA);
      fOutputList->Add(fhTrueVertexZB);
      fOutputList->Add(fhTrueVertexZA);
      fOutputList->Add(fhTrueVertexZAA);
    }
  }

  template <typename CollisionObject, typename TracksObject>
  void processReconstructed(CollisionObject const& collision, TracksObject const& ftracks, float centormult);

  void processWithCent(aod::CollisionEvSelCent const& collision, DptDptFullTracks const& ftracks);
  PROCESS_SWITCH(DptDptFilter, processWithCent, "Process reco with centrality", false);

  void processWithRun2Cent(aod::CollisionEvSelRun2Cent const& collision, DptDptFullTracks const& ftracks);
  PROCESS_SWITCH(DptDptFilter, processWithRun2Cent, "Process reco with Run !/2 centrality", false);

  void processWithoutCent(aod::CollisionEvSel const& collision, DptDptFullTracks const& ftracks);
  PROCESS_SWITCH(DptDptFilter, processWithoutCent, "Process reco without centrality", false);

  void processWithCentPID(aod::CollisionEvSelCent const& collision, DptDptFullTracksPID const& ftracks);
  PROCESS_SWITCH(DptDptFilter, processWithCentPID, "Process PID reco with centrality", false);

  void processWithRun2CentPID(aod::CollisionEvSelRun2Cent const& collision, DptDptFullTracksPID const& ftracks);
  PROCESS_SWITCH(DptDptFilter, processWithRun2CentPID, "Process PID reco with centrality", false);

  void processWithoutCentPID(aod::CollisionEvSel const& collision, DptDptFullTracksPID const& ftracks);
  PROCESS_SWITCH(DptDptFilter, processWithoutCentPID, "Process PID reco without centrality", false);

  void processWithCentDetectorLevel(aod::CollisionEvSelCent const& collision, DptDptFullTracksDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(DptDptFilter, processWithCentDetectorLevel, "Process MC detector level with centrality", false);

  void processWithRun2CentDetectorLevel(aod::CollisionEvSelRun2Cent const& collision, DptDptFullTracksDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(DptDptFilter, processWithRun2CentDetectorLevel, "Process MC detector level with centrality", false);

  void processWithoutCentDetectorLevel(aod::CollisionEvSel const& collision, DptDptFullTracksDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(DptDptFilter, processWithoutCentDetectorLevel, "Process MC detector level without centrality", false);

  void processWithCentPIDDetectorLevel(aod::CollisionEvSelCent const& collision, DptDptFullTracksPIDDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(DptDptFilter, processWithCentPIDDetectorLevel, "Process PID MC detector level with centrality", false);

  void processWithRun2CentPIDDetectorLevel(aod::CollisionEvSelRun2Cent const& collision, DptDptFullTracksPIDDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(DptDptFilter, processWithRun2CentPIDDetectorLevel, "Process PID MC detector level with centrality", false);

  void processWithoutCentPIDDetectorLevel(aod::CollisionEvSel const& collision, DptDptFullTracksPIDDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(DptDptFilter, processWithoutCentPIDDetectorLevel, "Process PID MC detector level without centrality", false);

  template <typename CollisionObject, typename ParticlesList>
  void processGenerated(CollisionObject const& mccollision, ParticlesList const& mcparticles, float centormult);

  template <typename CollisionsGroup, typename AllCollisions>
  void processGeneratorLevel(aod::McCollision const& mccollision,
                             CollisionsGroup const& collisions,
                             aod::McParticles const& mcparticles,
                             AllCollisions const& allcollisions,
                             float defaultcent);

  void processWithCentGeneratorLevel(aod::McCollision const& mccollision,
                                     soa::SmallGroups<soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>> const& collisions,
                                     aod::McParticles const& mcparticles,
                                     aod::CollisionsEvSelCent const& allcollisions);
  PROCESS_SWITCH(DptDptFilter, processWithCentGeneratorLevel, "Process generated with centrality", false);

  void processWithRun2CentGeneratorLevel(aod::McCollision const& mccollision,
                                         soa::SmallGroups<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>> const& collisions,
                                         aod::McParticles const& mcparticles,
                                         aod::CollisionsEvSelRun2Cent const& allcollisions);
  PROCESS_SWITCH(DptDptFilter, processWithRun2CentGeneratorLevel, "Process generated with centrality", false);

  void processWithoutCentGeneratorLevel(aod::McCollision const& mccollision,
                                        soa::SmallGroups<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>> const& collisions,
                                        aod::McParticles const& mcparticles,
                                        aod::CollisionsEvSel const& allcollisions);
  PROCESS_SWITCH(DptDptFilter, processWithoutCentGeneratorLevel, "Process generated without centrality", false);

  void processVertexGenerated(aod::McCollisions const&);
  PROCESS_SWITCH(DptDptFilter, processVertexGenerated, "Process vertex generator level", false);
};

template <typename CollisionObject, typename TracksObject>
void DptDptFilter::processReconstructed(CollisionObject const& collision, TracksObject const& ftracks, float tentativecentmult)
{
  using namespace dptdptfilter;

  LOGF(DPTDPTFILTERLOGCOLLISIONS, "DptDptFilterTask::processReconstructed(). New collision with %d tracks", ftracks.size());

  float mult = extractMultiplicity(collision, fCentMultEstimator);

  fhCentMultB->Fill(tentativecentmult);
  fhMultB->Fill(mult);
  fhVertexZB->Fill(collision.posZ());
  uint8_t acceptedevent = uint8_t(false);
  float centormult = tentativecentmult;
  if (IsEvtSelected(collision, centormult)) {
    acceptedevent = true;
    fhCentMultA->Fill(centormult);
    fhMultA->Fill(mult);
    fhVertexZA->Fill(collision.posZ());
    if (fullDerivedData) {
      acceptedcollisions(collision.bcId(), collision.posZ(), acceptedevent, centormult);
    } else {
      collisionsinfo(acceptedevent, centormult);
    }
  } else {
    if (!fullDerivedData) {
      /* the tracks are done at a different level */
      collisionsinfo(uint8_t(false), 105.0);
    }
  }
}

void DptDptFilter::processWithCent(aod::CollisionEvSelCent const& collision, DptDptFullTracks const& ftracks)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void DptDptFilter::processWithRun2Cent(aod::CollisionEvSelRun2Cent const& collision, DptDptFullTracks const& ftracks)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void DptDptFilter::processWithoutCent(aod::CollisionEvSel const& collision, DptDptFullTracks const& ftracks)
{
  processReconstructed(collision, ftracks, 50.0);
}

void DptDptFilter::processWithCentPID(aod::CollisionEvSelCent const& collision, DptDptFullTracksPID const& ftracks)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void DptDptFilter::processWithRun2CentPID(aod::CollisionEvSelRun2Cent const& collision, DptDptFullTracksPID const& ftracks)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void DptDptFilter::processWithoutCentPID(aod::CollisionEvSel const& collision, DptDptFullTracksPID const& ftracks)
{
  processReconstructed(collision, ftracks, 50.0);
}

void DptDptFilter::processWithCentDetectorLevel(aod::CollisionEvSelCent const& collision, DptDptFullTracksDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void DptDptFilter::processWithRun2CentDetectorLevel(aod::CollisionEvSelRun2Cent const& collision, DptDptFullTracksDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void DptDptFilter::processWithoutCentDetectorLevel(aod::CollisionEvSel const& collision, DptDptFullTracksDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, 50.0);
}

void DptDptFilter::processWithCentPIDDetectorLevel(aod::CollisionEvSelCent const& collision, DptDptFullTracksPIDDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void DptDptFilter::processWithRun2CentPIDDetectorLevel(aod::CollisionEvSelRun2Cent const& collision, DptDptFullTracksPIDDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void DptDptFilter::processWithoutCentPIDDetectorLevel(aod::CollisionEvSel const& collision, DptDptFullTracksPIDDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, 50.0);
}

template <typename CollisionObject, typename ParticlesList>
void DptDptFilter::processGenerated(CollisionObject const& mccollision, ParticlesList const& mcparticles, float centormult)
{
  using namespace dptdptfilter;

  uint8_t acceptedevent = uint8_t(false);
  if (IsEvtSelected(mccollision, centormult)) {
    acceptedevent = uint8_t(true);
  }
  if (fullDerivedData) {
    acceptedtrueevents(mccollision.bcId(), mccollision.posZ(), acceptedevent, centormult);
  } else {
    gencollisionsinfo(acceptedevent, centormult);
  }
}

template <typename CollisionsGroup, typename AllCollisions>
void DptDptFilter::processGeneratorLevel(aod::McCollision const& mccollision,
                                         CollisionsGroup const& collisions,
                                         aod::McParticles const& mcparticles,
                                         AllCollisions const& allcollisions,
                                         float defaultcent)
{
  using namespace dptdptfilter;

  LOGF(DPTDPTFILTERLOGCOLLISIONS, "DptDptFilterTask::processGeneratorLevel(). New generated collision with %d reconstructed collisions and %d particles", collisions.size(), mcparticles.size());

  if (collisions.size() > 1) {
    LOGF(DPTDPTFILTERLOGCOLLISIONS, "DptDptFilterTask::processGeneratorLevel(). Generated collision with more than one reconstructed collisions. Processing only the first accepted for centrality/multiplicity classes extraction");
  }

  bool processed = false;
  for (auto& tmpcollision : collisions) {
    if (tmpcollision.has_mcCollision()) {
      if (tmpcollision.mcCollisionId() == mccollision.globalIndex()) {
        typename AllCollisions::iterator const& collision = allcollisions.iteratorAt(tmpcollision.globalIndex());
        if (IsEvtSelected(collision, defaultcent)) {
          fhTrueVertexZAA->Fill((mccollision.posZ()));
          processGenerated(mccollision, mcparticles, defaultcent);
          processed = true;
          break; /* TODO: only processing the first reconstructed accepted collision */
        }
      }
    }
  }
  if (!processed && !fullDerivedData) {
    gencollisionsinfo(uint8_t(false), 105.0);
  }
}

void DptDptFilter::processWithCentGeneratorLevel(aod::McCollision const& mccollision,
                                                 soa::SmallGroups<soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>> const& collisions,
                                                 aod::McParticles const& mcparticles,
                                                 aod::CollisionsEvSelCent const& allcollisions)
{
  processGeneratorLevel(mccollision, collisions, mcparticles, allcollisions, 50.0);
}

void DptDptFilter::processWithRun2CentGeneratorLevel(aod::McCollision const& mccollision,
                                                     soa::SmallGroups<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>> const& collisions,
                                                     aod::McParticles const& mcparticles,
                                                     aod::CollisionsEvSelRun2Cent const& allcollisions)
{
  processGeneratorLevel(mccollision, collisions, mcparticles, allcollisions, 50.0);
}

void DptDptFilter::processWithoutCentGeneratorLevel(aod::McCollision const& mccollision,
                                                    soa::SmallGroups<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>> const& collisions,
                                                    aod::McParticles const& mcparticles,
                                                    aod::CollisionsEvSel const& allcollisions)
{
  processGeneratorLevel(mccollision, collisions, mcparticles, allcollisions, 50.0);
}

void DptDptFilter::processVertexGenerated(aod::McCollisions const& mccollisions)
{
  for (aod::McCollision const& mccollision : mccollisions) {
    fhTrueVertexZB->Fill(mccollision.posZ());
    /* we assign a default value */
    float centmult = 50.0f;
    if (IsEvtSelected(mccollision, centmult)) {
      fhTrueVertexZA->Fill((mccollision.posZ()));
    }
  }
}

/// RMS calculation. Taken from PWGHF/Tasks/taskMcValidation.cxx
/// \param vec  vector of values to compute RMS
template <typename T>
T computeRMS(std::vector<T>& vec)
{
  T sum = std::accumulate(vec.begin(), vec.end(), 0.0);
  T mean = sum / vec.size();

  std::vector<T> diff(vec.size());
  std::transform(vec.begin(), vec.end(), diff.begin(), [mean](T x) { return x - mean; });
  T sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  T stdev = std::sqrt(sq_sum / vec.size());

  return stdev;
}

struct DptDptFilterTracks {

  Produces<aod::ScannedTracks> scannedtracks;
  Produces<aod::DptDptCFTracksInfo> tracksinfo;
  Produces<aod::ScannedTrueTracks> scannedgentracks;
  Produces<aod::DptDptCFGenTracksInfo> gentracksinfo;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  bool storedccdbinfo = false;

  std::string cfgCCDBUrl{"http://ccdb-test.cern.ch:8080"};
  std::string cfgCCDBPathName{""};
  std::string cfgCCDBDate{"20220307"};
  std::string cfgCCDBPeriod{"LHC22o"};

  Configurable<bool> cfgFullDerivedData{"fullderiveddata", false, "Produce the full derived data for external storage. Default false"};
  Configurable<int> cfgTrackType{"trktype", 4, "Type of selected tracks: 0 = no selection, 1 = Run2 global tracks FB96, 3 = Run3 tracks, 4 = Run3 tracks MM sel, 5 = Run2 TPC only tracks, 7 = Run 3 TPC only tracks. Default 4"};
  Configurable<o2::analysis::CheckRangeCfg> cfgTraceDCAOutliers{"trackdcaoutliers", {false, 0.0, 0.0}, "Track the generator level DCAxy outliers: false/true, low dcaxy, up dcaxy. Default {false,0.0,0.0}"};
  Configurable<float> cfgTraceOutOfSpeciesParticles{"trackoutparticles", false, "Track the particles which are not e,mu,pi,K,p: false/true. Default false"};
  Configurable<int> cfgRecoIdMethod{"recoidmethod", 0, "Method for identifying reconstructed tracks: 0 No PID, 1 PID, 2 mcparticle. Default 0"};
  Configurable<o2::analysis::TrackSelectionCfg> cfgTrackSelection{"tracksel", {false, false, 0, 70, 0.8, 2.4, 3.2}, "Track selection: {useit: true/false, ongen: true/false, tpccls, tpcxrws, tpcxrfc, dcaxy, dcaz}. Default {false,0.70.0.8,2.4,3.2}"};
  Configurable<o2::analysis::TrackSelectionPIDCfg> cfgPionPIDSelection{"pipidsel",
                                                                       {},
                                                                       "PID criteria for pions"};
  Configurable<o2::analysis::TrackSelectionPIDCfg> cfgKaonPIDSelection{"kapidsel",
                                                                       {},
                                                                       "PID criteria for kaons"};
  Configurable<o2::analysis::TrackSelectionPIDCfg> cfgProtonPIDSelection{"prpidsel",
                                                                         {},
                                                                         "PID criteria for protons"};
  Configurable<o2::analysis::TrackSelectionPIDCfg> cfgElectronPIDSelection{"elpidsel",
                                                                           {},
                                                                           "PID criteria for electrons"};
  Configurable<o2::analysis::TrackSelectionPIDCfg> cfgMuonPIDSelection{"mupidsel",
                                                                       {},
                                                                       "PID criteria for muons"};

  OutputObj<TList> fOutput{"DptDptFilterTracksInfo", OutputObjHandlingPolicy::AnalysisObject};
  PIDSpeciesSelection pidselector;
  bool checkAmbiguousTracks = false;

  void init(InitContext& initContext)
  {
    LOGF(info, "DptDptFilterTracks::init()");

    fullDerivedData = cfgFullDerivedData;

    /* update with the configurable values */
    /* self configure the binning */
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mZVtxbins", zvtxbins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mZVtxmin", zvtxlow, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mZVtxmax", zvtxup, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPTbins", ptbins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPTmin", ptlow, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPTmax", ptup, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mEtabins", etabins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mEtamin", etalow, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mEtamax", etaup, false);

    /* self configure the CCDB access to the input file */
    getTaskOptionValue(initContext, "dpt-dpt-filter", "input_ccdburl", cfgCCDBUrl, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "input_ccdbpath", cfgCCDBPathName, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "input_ccdbdate", cfgCCDBDate, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "input_ccdbperiod", cfgCCDBPeriod, false);

    /* the track types and combinations */
    tracktype = cfgTrackType.value;
    initializeTrackSelection();
    traceDCAOutliers = cfgTraceDCAOutliers;
    traceOutOfSpeciesParticles = cfgTraceOutOfSpeciesParticles;
    recoIdMethod = cfgRecoIdMethod;
    if (cfgTrackSelection->mUseIt) {
      useOwnTrackSelection = true;
      if (cfgTrackSelection->mOnGen) {
        useOwnParticleSelection = true;
        particleMaxDCAxy = cfgTrackSelection->mDCAxy;
        particleMaxDCAZ = cfgTrackSelection->mDCAz;
      }
      ownTrackSelection.ResetITSRequirements();
      ownTrackSelection.SetMinNClustersTPC(cfgTrackSelection->mTPCclusters);
      ownTrackSelection.SetMinNCrossedRowsTPC(cfgTrackSelection->mTPCxRows);
      ownTrackSelection.SetMaxDcaXYPtDep([&](float) { return cfgTrackSelection->mDCAxy; });
      ownTrackSelection.SetMaxDcaXY(cfgTrackSelection->mDCAxy);
      ownTrackSelection.SetMaxDcaZ(cfgTrackSelection->mDCAz);
      o2::aod::track::TrackTypeEnum ttype;
      switch (tracktype) {
        case 1:
          ttype = o2::aod::track::Run2Track;
          break;
        case 3:
          ttype = o2::aod::track::Track;
          break;
        default:
          ttype = o2::aod::track::Track;
          break;
      }
      ownTrackSelection.SetTrackType(ttype);
    } else {
      useOwnTrackSelection = false;
    }

    /* self configure system type and data type */
    /* if the system type is not known at this time, we have to put the initialization somewhere else */
    std::string tmpstr;
    getTaskOptionValue(initContext, "dpt-dpt-filter", "syst", tmpstr, false);
    fSystem = getSystemType(tmpstr);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "datatype", tmpstr, false);
    fDataType = getDataType(tmpstr);
    fPDG = TDatabasePDG::Instance();

    /* required ambiguous tracks checks? */
    if (dofilterDetectorLevelWithoutPIDAmbiguous || dofilterDetectorLevelWithPIDAmbiguous || dofilterRecoWithoutPIDAmbiguous || dofilterRecoWithPIDAmbiguous) {
      checkAmbiguousTracks = true;
    }

    /* configure the PID selection */
    auto insertInPIDselector = [&](auto cfg, uint sp) {
      if (cfg.value.mUseIt) {
        if (cfg.value.mExclude) {
          pidselector.AddExclude(sp, &(cfg.value));
          LOGF(info, "Incorporated species: %s to PID selection for exclusion", pidselector.spnames[sp].data());
        } else {
          pidselector.Add(sp, &(cfg.value));
          LOGF(info, "Incorporated species: %s to PID selection", pidselector.spnames[sp].data());
        }
      }
    };
    insertInPIDselector(cfgElectronPIDSelection, 0);
    insertInPIDselector(cfgMuonPIDSelection, 1);
    insertInPIDselector(cfgPionPIDSelection, 2);
    insertInPIDselector(cfgKaonPIDSelection, 3);
    insertInPIDselector(cfgProtonPIDSelection, 4);

    /* create the output list which will own the task histograms */
    TList* fOutputList = new TList();
    fOutputList->SetOwner(true);
    fOutput.setObject(fOutputList);

    /* incorporate configuration parameters to the output */
    fOutputList->Add(new TParameter<Int_t>("TrackType", cfgTrackType, 'f'));
    fOutputList->Add(new TParameter<Int_t>("TrackOneCharge", 1, 'f'));
    fOutputList->Add(new TParameter<Int_t>("TrackTwoCharge", -1, 'f'));

    if ((fDataType == kData) || (fDataType == kDataNoEvtSel) || (fDataType == kMC)) {
      /* create the reconstructed data histograms */
      fhPB = new TH1F("fHistPB", "p distribution for reconstructed before;p (GeV/c);dN/dp (c/GeV)", 100, 0.0, 15.0);
      fhPtB = new TH1F("fHistPtB", "p_{T} distribution for reconstructed before;p_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhPtPosB = new TH1F("fHistPtPosB", "P_{T} distribution for reconstructed (#plus) before;P_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhPtNegB = new TH1F("fHistPtNegB", "P_{T} distribution for reconstructed (#minus) before;P_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhEtaB = new TH1F("fHistEtaB", "#eta distribution for reconstructed before;#eta;counts", 40, -2.0, 2.0);
      fhEtaA = new TH1F("fHistEtaA", "#eta distribution for reconstructed;#eta;counts", etabins, etalow, etaup);
      fhPhiB = new TH1F("fHistPhiB", "#phi distribution for reconstructed before;#phi;counts", 360, 0.0, constants::math::TwoPI);
      fhPhiA = new TH1F("fHistPhiA", "#phi distribution for reconstructed;#phi;counts", 360, 0.0, constants::math::TwoPI);
      fhDCAxyB = new TH1F("DCAxyB", "DCA_{xy} distribution for reconstructed before;DCA_{xy} (cm);counts", 1000, -4.0, 4.0);
      fhDCAxyA = new TH1F("DCAxyA", "DCA_{xy} distribution for reconstructed;DCA_{xy} (cm);counts", 1000, -4., 4.0);
      fhFineDCAxyA = new TH1F("FineDCAxyA", "DCA_{xy} distribution for reconstructed;DCA_{xy} (cm);counts", 4000, -1.0, 1.0);
      fhDCAzB = new TH1F("DCAzB", "DCA_{z} distribution for reconstructed before;DCA_{z} (cm);counts", 1000, -4.0, 4.0);
      fhDCAzA = new TH1F("DCAzA", "DCA_{z} distribution for reconstructed;DCA_{z} (cm);counts", 1000, -4.0, 4.0);
      fhFineDCAzA = new TH1F("FineDCAzA", "DCA_{z} distribution for reconstructed;DCA_{z} (cm);counts", 4000, -1.0, 1.0);

      if (checkAmbiguousTracks) {
        /* let's allocate the ambigous tracks tracking histograms*/
        fhAmbiguousTrackType = new TH2D("fHistAmbiguousTracksType", "Ambiguous tracks type vs. multiplicity class;Ambiguous track type;Multiplicity (%);counts", 4, -0.5, 3.5, 101, -0.5, 100.5);
        fhAmbiguousTrackPt = new TH2F("fHistAmbiguousTracksPt", "Ambiguous tracks #it{p}_{T} vs. multiplicity class;#it{p}_{T} (GeV/#it{c});Multiplicity (%);counts", 100, 0.0, 15.0, 101, -0.5, 100.5);
        fhAmbiguityDegree = new TH2F("fHistAmbiguityDegree", "Ambiguity degree vs. multiplicity class;Ambiguity degree;Multiplicity (%);counts", 31, -0.5, 30.5, 101, -0.5, 100.5);
        fhCompatibleCollisionsZVtxRms = new TH2F("fHistCompatibleCollisionsZVtxRms", "Compatible collisions #it{z}_{vtx} RMS;#sigma_{#it{z}_{vtx}};Multiplicity (%);counts", 100, -10.0, 10.0, 101, -0.5, 100.5);
      }

      uint nspecies = pidselector.getNSpecies();
      auto reserveHistos = [&](uint n) {
        fhPA.reserve(n);
        fhPtA.reserve(n);
        fhPtPosA.reserve(n);
        fhPtNegA.reserve(n);
        fhNPosNegA.reserve(n);
        fhDeltaNA.reserve(n);
        trkMultPos.reserve(n);
        trkMultNeg.reserve(n);
      };
      auto createHistos = [&](uint ix, auto name, auto title) {
        fhPA[ix] = new TH1F(TString::Format("fHistPA_%s", name).Data(),
                            TString::Format("#it{p} distribution for reconstructed %s;#it{p} (GeV/#it{c});d#it{N}/d#it{p} (#it{c}/GeV)", title).Data(),
                            ptbins, ptlow, ptup);
        fhPtA[ix] = new TH1F(TString::Format("fHistPtA_%s", name).Data(),
                             TString::Format("#it{p}_{T} distribution for reconstructed %s;#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{P}_{T} (#it{c}/GeV)", title).Data(),
                             ptbins, ptlow, ptup);
        fhPtPosA[ix] = new TH1F(TString::Format("fHistPtPosA_%s", name).Data(),
                                TString::Format("#it{p}_{T} distribution for reconstructed  %s^{#plus};#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T} (#it{c}/GeV)", title).Data(),
                                ptbins, ptlow, ptup);
        fhPtNegA[ix] = new TH1F(TString::Format("fHistPtNegA_%s", name).Data(),
                                TString::Format("#it{p}_{T} distribution for reconstructed  %s^{#minus};#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T} (#it{c}/GeV)", title).Data(),
                                ptbins, ptlow, ptup);
        fhNPosNegA[ix] = new TH2F(TString::Format("fhNPosNegA_%s", name).Data(),
                                  TString::Format("#it{N}(%s^{#plus}) #it{N}(%s^{#minus}) distribution for reconstructed;#it{N}(%s^{#plus});#it{N}(%s^{#minus})", title, title, title, title).Data(),
                                  40, -0.5, 39.5, 40, -0.5, 39.5);
        fhDeltaNA[ix] = new TH1F(TString::Format("fhDeltaNA_%s", name).Data(),
                                 TString::Format("#it{N}(%s^{#plus}) #minus #it{N}(%s^{#minus}) distribution for reconstructed;#it{N}(%s^{#plus}) #minus #it{N}(%s^{#minus})", title, title, title, title).Data(),
                                 79, -39.5, 39.5);
        trkMultPos[ix] = 0;
        trkMultNeg[ix] = 0;
      };
      if (nspecies == 0) {
        /* no species so charged tracks analysis*/
        LOGF(info, "Unidentified analysis");
        reserveHistos(1);
        createHistos(0, pidselector.getHadFName(), pidselector.getHadTitle());
      } else {
        reserveHistos(nspecies);
        LOGF(info, "Identified analysis with %d species", nspecies);
        for (uint sp = 0; sp < nspecies; ++sp) {
          LOGF(info, "Adding species %s", pidselector.getSpeciesFName(sp));
          createHistos(sp, pidselector.getSpeciesFName(sp), pidselector.getSpeciesTitle(sp));
        }
      }

      /* add the hstograms to the output list */
      fOutputList->Add(fhPB);
      fOutputList->Add(fhPtB);
      fOutputList->Add(fhPtPosB);
      fOutputList->Add(fhPtNegB);
      fOutputList->Add(fhEtaB);
      fOutputList->Add(fhEtaA);
      fOutputList->Add(fhPhiB);
      fOutputList->Add(fhPhiA);
      fOutputList->Add(fhDCAxyB);
      fOutputList->Add(fhDCAxyA);
      fOutputList->Add(fhFineDCAxyA);
      fOutputList->Add(fhDCAzB);
      fOutputList->Add(fhDCAzA);
      fOutputList->Add(fhFineDCAzA);
      if (checkAmbiguousTracks) {
        fOutputList->Add(fhAmbiguousTrackType);
        fOutputList->Add(fhAmbiguousTrackPt);
        fOutputList->Add(fhAmbiguityDegree);
        fOutputList->Add(fhCompatibleCollisionsZVtxRms);
      }
      uint nhsets = (nspecies > 0) ? nspecies : 1;
      for (uint sp = 0; sp < nhsets; ++sp) {
        fOutputList->Add(fhPA[sp]);
        fOutputList->Add(fhPtA[sp]);
        fOutputList->Add(fhPtPosA[sp]);
        fOutputList->Add(fhPtNegA[sp]);
        fOutputList->Add(fhNPosNegA[sp]);
        fOutputList->Add(fhDeltaNA[sp]);
      }
    }

    if ((fDataType != kData) && (fDataType != kDataNoEvtSel)) {
      /* create the true data histograms */
      fhTruePB = new TH1F("fTrueHistPB", "#it{p} distribution before (truth);#it{p} (GeV/#it{c});d#it{N}/d#it{p} (#it{c}/GeV)", 100, 0.0, 15.0);
      fhTruePtB = new TH1F("fTrueHistPtB", "#it{p}_{T} distribution before (truth);#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T} (#it{c}/GeV)", 100, 0.0, 15.0);
      fhTruePtPosB = new TH1F("fTrueHistPtPosB", "#it{p}_{T} distribution (#plus) before (truth);#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T} (#it{c}/GeV)", 100, 0.0, 15.0);
      fhTruePtNegB = new TH1F("fTrueHistPtNegB", "#it{p}_{T} distribution (#minus) before (truth);#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T} (#it{c}/GeV)", 100, 0.0, 15.0);
      fhTrueEtaB = new TH1F("fTrueHistEtaB", "#eta distribution before (truth);#eta;counts", 40, -2.0, 2.0);
      fhTrueEtaA = new TH1F("fTrueHistEtaA", "#eta distribution (truth);#eta;counts", etabins, etalow, etaup);
      fhTruePhiB = new TH1F("fTrueHistPhiB", "#phi distribution before (truth);#phi;counts", 360, 0.0, constants::math::TwoPI);
      fhTruePhiA = new TH1F("fTrueHistPhiA", "#phi distribution (truth);#phi;counts", 360, 0.0, constants::math::TwoPI);
      fhTrueDCAxyB = new TH1F("TrueDCAxyB", "DCA_{#it{xy}} distribution for generated before;DCA_{#it{xy}} (cm);counts", 1000, -4.0, 4.0);
      if (traceDCAOutliers.mDoIt) {
        fhTrueDCAxyBid = new TH1F("PDGCodeDCAxyB",
                                  TString::Format("PDG code within %.2f<|DCA_{#it{xy}}|<%.2f; PDG code", traceDCAOutliers.mLowValue, traceDCAOutliers.mUpValue).Data(),
                                  100, 0.5, 100.5);
      }
      fhTrueDCAxyA = new TH1F("TrueDCAxyA", "DCA_{#it{xy}} distribution for generated;DCA_{#it{xy}};counts (cm)", 1000, -4., 4.0);
      fhTrueDCAzB = new TH1F("TrueDCAzB", "DCA_{#it{z}} distribution for generated before;DCA_{#it{z}} (cm);counts", 1000, -4.0, 4.0);
      fhTrueDCAzA = new TH1F("TrueDCAzA", "DCA_{#it{z}} distribution for generated;DCA_{#it{z}} (cm);counts", 1000, -4.0, 4.0);

      auto reserveTruthHistos = [&](uint n) {
        fhTruePA.reserve(n);
        fhTruePtA.reserve(n);
        fhTruePtPosA.reserve(n);
        fhTruePtNegA.reserve(n);
        fhTrueNPosNegA.reserve(n);
        fhTrueDeltaNA.reserve(n);
        partMultPos.reserve(n);
        partMultNeg.reserve(n);
      };
      auto createTruthHistos = [&](uint ix, auto name, auto title) {
        fhTruePA[ix] = new TH1F(TString::Format("fTrueHistPA_%s", name).Data(),
                                TString::Format("#it{p} distribution %s (truth);#it{p} (GeV/#it{c});d#it{N}/d#it{p} (#it{c}/GeV)", title).Data(),
                                ptbins, ptlow, ptup);
        fhTruePtA[ix] = new TH1F(TString::Format("fTrueHistPtA_%s", name).Data(),
                                 TString::Format("#it{p}_{T} distribution %s (truth);#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T} (#it{c}/GeV)", title).Data(),
                                 ptbins, ptlow, ptup);
        fhTruePtPosA[ix] = new TH1F(TString::Format("fTrueHistPtPosA_%s", name).Data(),
                                    TString::Format("#it{p}_{T} distribution %s^{#plus} (truth);#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T} (#it{c}/GeV)", title).Data(),
                                    ptbins, ptlow, ptup);
        fhTruePtNegA[ix] = new TH1F(TString::Format("fTrueHistPtNegA_%s", name).Data(),
                                    TString::Format("#it{p}_{T} distribution %s^{#minus} (truth);#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T} (#it{c}/GeV)", title).Data(),
                                    ptbins, ptlow, ptup);
        fhTrueNPosNegA[ix] = new TH2F(TString::Format("fhTrueNPosNegA_%s", name).Data(),
                                      TString::Format("#it{N}(%s^{#plus}) #it{N}(%s^{#minus}) distribution (truth);#it{N}(%s^{#plus});#it{N}(%s^{#minus})", title, title, title, title).Data(),
                                      40, -0.5, 39.5, 40, -0.5, 39.5);
        fhTrueDeltaNA[ix] = new TH1F(TString::Format("fhTrueDeltaNA_%s", name).Data(),
                                     TString::Format("#it{N}(%s^{#plus}) #minus #it{N}(%s^{#minus}) distribution (truth);#it{N}(%s^{#plus}) #minus #it{N}(%s^{#minus})", title, title, title, title).Data(),
                                     79, -39.5, 39.5);
        partMultPos[ix] = 0;
        partMultNeg[ix] = 0;
      };
      uint nspecies = pidselector.getNSpecies();
      if (nspecies == 0) {
        reserveTruthHistos(1);
        createTruthHistos(0, pidselector.getHadFName(), pidselector.getHadTitle());
      } else {
        reserveTruthHistos(nspecies);
        for (uint sp = 0; sp < nspecies; ++sp) {
          createTruthHistos(sp, pidselector.getSpeciesFName(sp), pidselector.getSpeciesTitle(sp));
        }
      }

      /* add the hstograms to the output list */
      fOutputList->Add(fhTruePB);
      fOutputList->Add(fhTruePtB);
      fOutputList->Add(fhTruePtPosB);
      fOutputList->Add(fhTruePtNegB);
      fOutputList->Add(fhTrueEtaB);
      fOutputList->Add(fhTrueEtaA);
      fOutputList->Add(fhTruePhiB);
      fOutputList->Add(fhTruePhiA);
      fOutputList->Add(fhTrueDCAxyB);
      if (traceDCAOutliers.mDoIt) {
        fOutputList->Add(fhTrueDCAxyBid);
      }
      fOutputList->Add(fhTrueDCAxyA);
      fOutputList->Add(fhTrueDCAzB);
      fOutputList->Add(fhTrueDCAzA);
      uint nhsets = (nspecies > 0) ? nspecies : 1;
      for (uint sp = 0; sp < nhsets; ++sp) {
        fOutputList->Add(fhTruePA[sp]);
        fOutputList->Add(fhTruePtA[sp]);
        fOutputList->Add(fhTruePtPosA[sp]);
        fOutputList->Add(fhTruePtNegA[sp]);
        fOutputList->Add(fhTrueNPosNegA[sp]);
        fOutputList->Add(fhTrueDeltaNA[sp]);
      }
    }
    /* initialize access to the CCDB */
    ccdb->setURL(cfgCCDBUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void getCCDBInformation()
  {
    using namespace analysis::dptdptfilter;

    /* let's get a potential PID adjustment */
    if (cfgCCDBPathName.length() > 0 && !storedccdbinfo) {
      LOGF(info, "Getting information for PID adjustment from %s, at %s, for %s", cfgCCDBPathName.c_str(), cfgCCDBDate.c_str(), cfgCCDBPeriod.c_str());
      TList* pidinfo = getCCDBInput(ccdb, cfgCCDBPathName.c_str(), cfgCCDBDate.c_str(), cfgCCDBPeriod.c_str());
      if (pidinfo != nullptr) {
        pidselector.storePIDAdjustments(pidinfo);
      }
      storedccdbinfo = true;
    }
  }

  template <typename TrackObject>
  int8_t trackIdentification(TrackObject const& track);
  template <typename TrackObject>
  int8_t selectTrack(TrackObject const& track);
  template <typename CollisionObjects, typename TrackObject>
  int8_t selectTrackAmbiguousCheck(CollisionObjects const& collisions, TrackObject const& track);
  template <typename ParticleObject>
  int8_t identifyParticle(ParticleObject const& particle);
  template <typename ParticleObject, typename MCCollisionObject>
  int8_t selectParticle(ParticleObject const& particle, MCCollisionObject const& mccollision);
  template <typename TrackObject>
  void fillTrackHistosBeforeSelection(TrackObject const& track);
  template <typename TrackObject>
  void fillTrackHistosAfterSelection(TrackObject const& track, int8_t sp);
  template <typename ParticleObject, typename MCCollisionObject>
  void fillParticleHistosBeforeSelection(ParticleObject const& particle,
                                         MCCollisionObject const& collision,
                                         float charge);
  template <typename ParticleObject, typename MCCollisionObject>
  void fillParticleHistosAfterSelection(ParticleObject const& particle,
                                        MCCollisionObject const& collision,
                                        float charge,
                                        int8_t sp);

  /* TODO: as it is now when the derived data is stored (fullDerivedData = true) */
  /* the collision index stored with the track is wrong. This has to be fixed    */
  template <typename passedtracks>
  void filterTracks(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const& collisions,
                    passedtracks const& tracks)
  {
    /* do check for special adjustments */
    getCCDBInformation();

    int naccepted = 0;
    int ncollaccepted = 0;
    if (!fullDerivedData) {
      tracksinfo.reserve(tracks.size());
    }
    for (auto collision : collisions) {
      if (collision.collisionaccepted()) {
        ncollaccepted++;
      }
    }
    for (auto track : tracks) {
      int8_t pid = -1;
      if (track.has_collision() && (track.template collision_as<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>()).collisionaccepted()) {
        pid = selectTrackAmbiguousCheck(collisions, track);
        if (!(pid < 0)) {
          naccepted++;
          if (fullDerivedData) {
            LOGF(fatal, "Stored derived data not prepared for saving the proper new collision id");
            scannedtracks((track.template collision_as<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>()).globalIndex(), pid, track.pt(), track.eta(), track.phi());
          } else {
            tracksinfo(pid);
          }
        } else {
          if (!fullDerivedData) {
            /* track not accepted */
            tracksinfo(pid);
          }
        }
      } else {
        if (!fullDerivedData) {
          /* track collision not accepted */
          tracksinfo(pid);
        }
      }
    }
    LOGF(DPTDPTFILTERLOGCOLLISIONS,
         "Processed %d accepted collisions out of a total of %d with  %d accepted tracks out of a "
         "total of %d",
         ncollaccepted,
         collisions.size(),
         naccepted,
         tracks.size());
  }

  /* TODO: for the time being the full derived data is still not supported  */
  /* for doing that we need to get the index of the associated mc collision */
  void filterParticles(soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo> const& gencollisions, aod::McParticles const& particles)
  {
    using namespace dptdptfilter;

    int acceptedparticles = 0;
    int acceptedcollisions = 0;
    if (!fullDerivedData) {
      gentracksinfo.reserve(particles.size());
    }

    for (auto gencoll : gencollisions) {
      if (gencoll.collisionaccepted()) {
        acceptedcollisions++;
      }
    }

    for (auto& particle : particles) {
      float charge = getCharge(particle);

      int8_t pid = -1;
      if (charge != 0) {
        if (particle.has_mcCollision() && (particle.template mcCollision_as<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>()).collisionaccepted()) {
          auto mccollision = particle.template mcCollision_as<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>();
          pid = selectParticle(particle, mccollision);
          if (!(pid < 0)) {
            acceptedparticles++;
          }
        }
      } else {
        if ((particle.mcCollisionId() == 0) && traceCollId0) {
          LOGF(DPTDPTFILTERLOGTRACKS, "Particle %d with fractional charge or equal to zero", particle.globalIndex());
        }
      }
      if (!fullDerivedData) {
        gentracksinfo(pid);
      }
    }
    LOGF(DPTDPTFILTERLOGCOLLISIONS,
         "Processed %d accepted generated collisions out of a total of %d with  %d accepted particles out of a "
         "total of %d",
         acceptedcollisions,
         gencollisions.size(),
         acceptedparticles,
         particles.size());
  }

  void filterRecoWithPID(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksPID const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterRecoWithPID, "Not stored derived data track filtering", false)

  void filterRecoWithPIDAmbiguous(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksPIDAmbiguous const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterRecoWithPIDAmbiguous, "Not stored derived data track filtering with ambiguous tracks check", false)

  void filterDetectorLevelWithPID(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksPIDDetLevel const& tracks, aod::McParticles const&)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterDetectorLevelWithPID, "Not stored derived data detector level track filtering", false)

  void filterDetectorLevelWithPIDAmbiguous(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksPIDDetLevelAmbiguous const& tracks, aod::McParticles const&)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterDetectorLevelWithPIDAmbiguous, "Not stored derived data detector level track filtering with ambiguous tracks check", false)

  void filterRecoWithoutPID(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const& collisions, DptDptFullTracks const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterRecoWithoutPID, "Track filtering without PID information", true)

  void filterRecoWithoutPIDAmbiguous(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const& collisions, DptDptFullTracksAmbiguous const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterRecoWithoutPIDAmbiguous, "Track filtering without PID information with ambiguous tracks check", false)

  void filterDetectorLevelWithoutPID(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const& collisions, DptDptFullTracksDetLevel const& tracks, aod::McParticles const&)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterDetectorLevelWithoutPID, "Detector level track filtering without PID information", false)

  void filterDetectorLevelWithoutPIDAmbiguous(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const& collisions, DptDptFullTracksDetLevelAmbiguous const& tracks, aod::McParticles const&)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterDetectorLevelWithoutPIDAmbiguous, "Detector level track filtering without PID information with ambiguous tracks check", false)

  void filterGenerated(soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo> const& gencollisions, aod::McParticles const& particles)
  {
    filterParticles(gencollisions, particles);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterGenerated, "Generated particles filtering", true)
};

template <typename TrackObject>
int8_t DptDptFilterTracks::trackIdentification(TrackObject const& track)
{
  using namespace dptdptfilter;

  int8_t sp = -127;
  if (recoIdMethod == 0) {
    sp = 0;
  } else if (recoIdMethod == 1) {
    if constexpr (framework::has_type_v<aod::pidtpc_tiny::TPCNSigmaStorePi, typename TrackObject::all_columns>) {
      sp = pidselector.whichSpecies(track);
    } else {
      LOGF(fatal, "Track identification required but PID information not present");
    }
  } else if (recoIdMethod == 2) {
    if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
      sp = identifyParticle(track.template mcParticle_as<aod::McParticles>());
    } else {
      LOGF(fatal, "Track identification required from MC particle but MC information not present");
    }
  }
  return sp;
}

template <typename TrackObject>
int8_t DptDptFilterTracks::selectTrack(TrackObject const& track)
{
  using namespace dptdptfilter;

  /* before track selection */
  fillTrackHistosBeforeSelection(track);

  /* track selection */
  int8_t sp = -127;
  if (AcceptTrack(track)) {
    /* the track has been accepted */
    /* let's identify it */
    sp = trackIdentification(track);
    if (!(sp < 0)) {
      /* fill the species histograms */
      fillTrackHistosAfterSelection(track, sp);
      /* update species multiplicities */
      if (track.sign() > 0) {
        trkMultPos[sp]++;
        /* positive tracks even pid */
        sp = sp * 2;
      } else if (track.sign() < 0) {
        trkMultNeg[sp]++;
        /* negative tracks odd pid */
        sp = sp * 2 + 1;
      }
    }
  }
  return sp;
}

template <typename CollisionObjects, typename TrackObject>
int8_t DptDptFilterTracks::selectTrackAmbiguousCheck(CollisionObjects const& collisions, TrackObject const& track)
{
  bool ambiguoustrack = false;
  int tracktype = 0; /* no ambiguous */
  std::vector<double> zvertexes{};
  /* ambiguous tracks checks if required */
  if constexpr (has_type_v<aod::track_association::CollisionIds, typename TrackObject::all_columns>) {
    if (track.compatibleCollIds().size() > 0) {
      if (track.compatibleCollIds().size() == 1) {
        if (track.collisionId() != track.compatibleCollIds()[0]) {
          /* ambiguous track! */
          ambiguoustrack = true;
          /* in principle we should not be here because the track is associated to two collisions at least */
          tracktype = 2;
          zvertexes.push_back(collisions.iteratorAt(track.collisionId()).posZ());
          zvertexes.push_back(collisions.iteratorAt(track.compatibleCollIds()[0]).posZ());
        } else {
          /* we consider the track as no ambiguous */
          tracktype = 1;
        }
      } else {
        /* ambiguous track! */
        ambiguoustrack = true;
        tracktype = 3;
        /* the track is associated to more than one collision */
        for (const auto& collIdx : track.compatibleCollIds()) {
          zvertexes.push_back(collisions.iteratorAt(collIdx).posZ());
        }
      }
    }
  }

  float multiplicityclass = (track.template collision_as<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>()).centmult();
  if (ambiguoustrack) {
    /* keep track of ambiguous tracks */
    fhAmbiguousTrackType->Fill(tracktype, multiplicityclass);
    fhAmbiguousTrackPt->Fill(track.pt(), multiplicityclass);
    fhAmbiguityDegree->Fill(zvertexes.size(), multiplicityclass);
    if (tracktype == 2) {
      fhCompatibleCollisionsZVtxRms->Fill(-computeRMS(zvertexes), multiplicityclass);
    } else {
      fhCompatibleCollisionsZVtxRms->Fill(computeRMS(zvertexes), multiplicityclass);
    }
    return -1;
  } else {
    if (checkAmbiguousTracks) {
      /* feedback of no ambiguous tracks only if checks required */
      fhAmbiguousTrackType->Fill(tracktype, multiplicityclass);
    }
    return selectTrack(track);
  }
}

template <typename TrackObject>
void DptDptFilterTracks::fillTrackHistosBeforeSelection(TrackObject const& track)
{
  fhPB->Fill(track.p());
  fhPtB->Fill(track.pt());
  fhEtaB->Fill(track.eta());
  fhPhiB->Fill(track.phi());
  if (track.sign() > 0) {
    fhPtPosB->Fill(track.pt());
  } else {
    fhPtNegB->Fill(track.pt());
  }
  fhDCAxyB->Fill(track.dcaXY());
  fhDCAzB->Fill(track.dcaZ());
}

template <typename TrackObject>
void DptDptFilterTracks::fillTrackHistosAfterSelection(TrackObject const& track, int8_t sp)
{
  fhEtaA->Fill(track.eta());
  fhPhiA->Fill(track.phi());
  fhDCAxyA->Fill(track.dcaXY());
  fhDCAzA->Fill(track.dcaZ());
  if (track.dcaXY() < 1.0) {
    fhFineDCAxyA->Fill(track.dcaXY());
  }
  if (track.dcaZ() < 1.0) {
    fhFineDCAzA->Fill(track.dcaZ());
  }
  fhPA[sp]->Fill(track.p());
  fhPtA[sp]->Fill(track.pt());
  if (track.sign() > 0) {
    fhPtPosA[sp]->Fill(track.pt());
  } else {
    fhPtNegA[sp]->Fill(track.pt());
  }
}

template <typename ParticleObject>
inline int8_t DptDptFilterTracks::identifyParticle(ParticleObject const& particle)
{
  using namespace dptdptfilter;
  return pidselector.whichTruthSpecies(particle);
}

template <typename ParticleObject, typename MCCollisionObject>
inline int8_t DptDptFilterTracks::selectParticle(ParticleObject const& particle, MCCollisionObject const& mccollision)
{
  float charge = getCharge(particle);
  int8_t sp = -127;
  if (charge != 0) {
    /* before particle selection */
    fillParticleHistosBeforeSelection(particle, mccollision, charge);

    /* track selection */
    if (AcceptParticle(particle, mccollision)) {
      /* the particle has been accepted */
      /* let's identify the particle */
      sp = identifyParticle(particle);
      if (!(sp < 0)) {
        /* fill the species  histograms */
        fillParticleHistosAfterSelection(particle, mccollision, charge, sp);
        /* update species multiplicities */
        if (charge > 0) {
          partMultPos[sp]++;
          /* positive species even pid */
          sp = sp * 2;
        } else {
          partMultNeg[sp]++;
          /* negative species odd pid */
          sp = sp * 2 + 1;
        }
      }
    }
  }
  return sp;
}

template <typename ParticleObject, typename MCCollisionObject>
void DptDptFilterTracks::fillParticleHistosBeforeSelection(ParticleObject const& particle, MCCollisionObject const& collision, float charge)
{
  fhTruePB->Fill(particle.p());
  fhTruePtB->Fill(particle.pt());
  fhTrueEtaB->Fill(particle.eta());
  fhTruePhiB->Fill(particle.phi());
  if (charge > 0) {
    fhTruePtPosB->Fill(particle.pt());
  } else if (charge < 0) {
    fhTruePtNegB->Fill(particle.pt());
  }

  float dcaxy = TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                            (particle.vy() - collision.posY()) * (particle.vy() - collision.posY()));
  if (traceDCAOutliers.mDoIt && (traceDCAOutliers.mLowValue < dcaxy) && (dcaxy < traceDCAOutliers.mUpValue)) {
    fhTrueDCAxyBid->Fill(TString::Format("%d", particle.pdgCode()).Data(), 1.0);
  }

  fhTrueDCAxyB->Fill(TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                                 (particle.vy() - collision.posY()) * (particle.vy() - collision.posY())));
  fhTrueDCAzB->Fill((particle.vz() - collision.posZ()));
}

template <typename ParticleObject, typename MCCollisionObject>
void DptDptFilterTracks::fillParticleHistosAfterSelection(ParticleObject const& particle, MCCollisionObject const& collision, float charge, int8_t sp)
{
  fhTrueEtaA->Fill(particle.eta());
  fhTruePhiA->Fill(particle.phi());
  float dcaxy = TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                            (particle.vy() - collision.posY()) * (particle.vy() - collision.posY()));
  if (traceDCAOutliers.mDoIt && (traceDCAOutliers.mLowValue < dcaxy) && (dcaxy < traceDCAOutliers.mUpValue)) {
    LOGF(info, "DCAxy outlier: Particle with index %d and pdg code %d assigned to MC collision %d, pT: %f, phi: %f, eta: %f",
         particle.globalIndex(), particle.pdgCode(), particle.mcCollisionId(), particle.pt(), particle.phi(), particle.eta());
    LOGF(info, "               With status %d and flags %0X", particle.statusCode(), particle.flags());
  }

  fhTrueDCAxyA->Fill(TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                                 (particle.vy() - collision.posY()) * (particle.vy() - collision.posY())));
  fhTrueDCAzA->Fill((particle.vz() - collision.posZ()));
  fhTruePA[sp]->Fill(particle.p());
  fhTruePtA[sp]->Fill(particle.pt());
  if (charge > 0) {
    fhTruePtPosA[sp]->Fill(particle.pt());
  } else {
    fhTruePtNegA[sp]->Fill(particle.pt());
  }
}

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<DptDptFilter>(cfgc,
                                                        SetDefaultProcesses{
                                                          {{"processWithoutCent", true},
                                                           {"processWithoutCentMC", true}}}),
                        adaptAnalysisTask<DptDptFilterTracks>(cfgc)};
  return workflow;
}
