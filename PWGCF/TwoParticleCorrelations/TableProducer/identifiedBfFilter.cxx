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
#include "PWGCF/TwoParticleCorrelations/TableProducer/identifiedBfFilter.h"

#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/IdentifiedBfFiltered.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Framework/runDataProcessing.h"
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

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis;

#define IDENTIFIEDBFFILTERLOGCOLLISIONS debug
#define IDENTIFIEDBFFILTERLOGTRACKS debug

namespace o2::analysis::identifiedbffilter
{
using IdBfFullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
using IdBfFullTracksAmbiguous = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackCompColls>;
using IdBfTracksPID = soa::Join<aod::pidTPCEl, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFEl, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
using IdBfTracksFullPID = soa::Join<aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
using IdBfFullTracksPID = soa::Join<IdBfFullTracks, IdBfTracksPID>;
using IdBfFullTracksFullPID = soa::Join<IdBfFullTracks, IdBfTracksFullPID>;
using IdBfFullTracksPIDAmbiguous = soa::Join<IdBfFullTracksAmbiguous, IdBfTracksPID>;
using IdBfFullTracksFullPIDAmbiguous = soa::Join<IdBfFullTracksAmbiguous, IdBfTracksFullPID>;
using IdBfFullTracksDetLevel = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TracksDCA>;
using IdBfFullTracksDetLevelAmbiguous = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TracksDCA, aod::TrackCompColls>;
using IdBfFullTracksPIDDetLevel = soa::Join<IdBfFullTracksDetLevel, IdBfTracksPID>;
using IdBfFullTracksFullPIDDetLevel = soa::Join<IdBfFullTracksDetLevel, IdBfTracksFullPID>;
using IdBfFullTracksPIDDetLevelAmbiguous = soa::Join<IdBfFullTracksDetLevelAmbiguous, IdBfTracksPID>;
using IdBfFullTracksFullPIDDetLevelAmbiguous = soa::Join<IdBfFullTracksDetLevelAmbiguous, IdBfTracksFullPID>;

bool fullDerivedData = false; /* produce full derived data for its external storage */

TList* ccdblst = nullptr;
bool loadfromccdb = false;

//============================================================================================
// The IdentifiedBfFilter histogram objects
// TODO: consider registering in the histogram registry
//============================================================================================
TH1F* fhCentMultB = nullptr;
TH1F* fhCentMultA = nullptr;
TH1F* fhVertexZB = nullptr;
TH1F* fhVertexZA = nullptr;
TH1F* fhMultB = nullptr;
TH1F* fhMultA = nullptr;
TH2F* fhYZB = nullptr;
TH2F* fhXYB = nullptr;
TH2F* fhYZA = nullptr;
TH2F* fhXYA = nullptr;
TH1F* fhPB = nullptr;
TH1F* fhPA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhPtB = nullptr;
TH1F* fhPtA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhPtPosB = nullptr;
TH1F* fhPtPosA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhPtNegB = nullptr;
TH1F* fhPtNegA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhNPosNegA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhDeltaNA[kIdBfNoOfSpecies + 1] = {nullptr};

TH1I* fhNClustersB = nullptr;
TH2F* fhPhiYB = nullptr;
TH2F* fhPtYB = nullptr;
TH1F* fhChi2B = nullptr;
TH1I* fhITSNclB = nullptr;

TH1I* fhNClustersA = nullptr;
TH2F* fhPhiYA = nullptr;
TH2F* fhPtYA = nullptr;
TH1F* fhChi2A = nullptr;
TH1I* fhITSNclA = nullptr;

TH2F* fhNSigmaTPC[kIdBfNoOfSpecies] = {nullptr};
TH2F* fhNSigmaTOF[kIdBfNoOfSpecies] = {nullptr};
TH2F* fhNSigmaCombo[kIdBfNoOfSpecies] = {nullptr};
TH2F* fhNSigmaTPC_IdTrks[kIdBfNoOfSpecies] = {nullptr};

TH1F* fhNSigmaCorrection[kIdBfNoOfSpecies] = {nullptr};

TH1F* fhEtaB = nullptr;
TH1F* fhEtaA = nullptr;

TH1F* fhPhiB = nullptr;
TH1F* fhPhiA = nullptr;

TH2F* fhdEdxB = nullptr;
TH2F* fhdEdxIPTPCB = nullptr;
TH2F* fhdEdxA[kIdBfNoOfSpecies + 2] = {nullptr};
TH2F* fhdEdxIPTPCA[kIdBfNoOfSpecies + 2] = {nullptr};

TH1F* fhMassB = nullptr;
TH1F* fhMassA[kIdBfNoOfSpecies + 1] = {nullptr};

TH2S* fhDoublePID = nullptr;

TH1F* fhDCAxyB = nullptr;
TH1F* fhDCAxyA = nullptr;
TH1F* fhFineDCAxyA = nullptr;
TH1F* fhDCAzB = nullptr;
TH2F* fhDCAxyzB = nullptr;
TH1F* fhDCAzA = nullptr;
TH1F* fhFineDCAzA = nullptr;
TH2F* fhDCAxyzA = nullptr;

TH1F* fhWrongTrackID = nullptr;

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
TH1F* fhTrueCharge = nullptr;
TH1F* fhTruePA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhTruePtB = nullptr;
TH1F* fhTruePtA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhTruePtPosB = nullptr;
TH1F* fhTruePtPosA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhTruePtNegB = nullptr;
TH1F* fhTruePtNegA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhTrueNPosNegA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhTrueDeltaNA[kIdBfNoOfSpecies + 1] = {nullptr};

TH2F* fhTruePhiYB = nullptr;
TH2F* fhTruePtYB = nullptr;

TH2F* fhTruePhiYA = nullptr;
TH2F* fhTruePtYA = nullptr;

TH1F* fhTrueEtaB = nullptr;
TH1F* fhTrueEtaA = nullptr;

TH1F* fhTruePhiB = nullptr;
TH1F* fhTruePhiA = nullptr;

TH1F* fhTrueDCAxyB = nullptr;
TH1F* fhTrueDCAxyA = nullptr;
TH1F* fhTrueDCAzB = nullptr;
TH1F* fhTrueDCAxyBid = nullptr;
TH1F* fhTrueDCAzA = nullptr;
TH2F* fhTrueDCAxyzB = nullptr;
TH2F* fhTrueDCAxyzA = nullptr;

//============================================================================================
// The IdentifiedBfFilter multiplicity counters
//============================================================================================
int trkMultPos[kIdBfNoOfSpecies + 1];  // multiplicity of positive tracks
int trkMultNeg[kIdBfNoOfSpecies + 1];  // multiplicity of negative tracks
int partMultPos[kIdBfNoOfSpecies + 1]; // multiplicity of positive particles
int partMultNeg[kIdBfNoOfSpecies + 1]; // multiplicity of negative particles
} // namespace o2::analysis::identifiedbffilter

using namespace identifiedbffilter;

struct IdentifiedBfFilter {
  Configurable<bool> cfgFullDerivedData{"fullderiveddata", false, "Produce the full derived data for external storage. Default false"};
  Configurable<std::string> cfgCentMultEstimator{"centmultestimator", "V0M", "Centrality/multiplicity estimator detector: V0M,CL0,CL1,FV0A,FT0M,FT0A,FT0C,NTPV,NOCM: none. Default V0M"};
  Configurable<std::string> cfgSystem{"syst", "PbPb", "System: pp, PbPb, Pbp, pPb, XeXe, ppRun3, PbPbRun3. Default PbPb"};
  Configurable<std::string> cfgDataType{"datatype", "data", "Data type: data, datanoevsel, MC, FastMC, OnTheFlyMC. Default data"};
  Configurable<std::string> cfgTriggSel{"triggsel", "MB", "Trigger selection: MB, None. Default MB"};
  Configurable<o2::analysis::DptDptBinningCuts> cfgBinning{"binning",
                                                           {28, -7.0, 7.0, 18, 0.2, 2.0, 16, -0.8, 0.8, 72, 0.5},
                                                           "triplets - nbins, min, max - for z_vtx, pT, eta and phi, binning plus bin fraction of phi origin shift"};
  Configurable<bool> cfgTraceCollId0{"tracecollid0", false, "Trace particles in collisions id 0. Default false"};

  OutputObj<TList> fOutput{"IdentifiedBfFilterCollisionsInfo", OutputObjHandlingPolicy::AnalysisObject};

  Produces<aod::IdentifiedBfCFAcceptedCollisions> acceptedcollisions;
  Produces<aod::IdentifiedBfCFCollisionsInfo> collisionsinfo;
  Produces<aod::IdentifiedBfCFAcceptedTrueCollisions> acceptedtrueevents;
  Produces<aod::IdentifiedBfCFGenCollisionsInfo> gencollisionsinfo;

  void init(InitContext const&)
  {
    using namespace identifiedbffilter;

    LOGF(info, "IdentifiedBfFilterTask::init()");

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

  void processWithCent(aod::CollisionEvSelCent const& collision, IdBfFullTracks const& ftracks);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithCent, "Process reco with centrality", false);

  void processWithRun2Cent(aod::CollisionEvSelRun2Cent const& collision, IdBfFullTracks const& ftracks);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithRun2Cent, "Process reco with Run !/2 centrality", false);

  void processWithoutCent(aod::CollisionEvSel const& collision, IdBfFullTracks const& ftracks);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithoutCent, "Process reco without centrality", false);

  void processWithCentPID(aod::CollisionEvSelCent const& collision, IdBfFullTracksPID const& ftracks);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithCentPID, "Process PID reco with centrality", false);

  void processWithRun2CentPID(aod::CollisionEvSelRun2Cent const& collision, IdBfFullTracksPID const& ftracks);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithRun2CentPID, "Process PID reco with centrality", false);

  void processWithoutCentPID(aod::CollisionEvSel const& collision, IdBfFullTracksPID const& ftracks);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithoutCentPID, "Process PID reco without centrality", false);

  void processWithCentFullPID(aod::CollisionEvSelCent const& collision, IdBfFullTracksFullPID const& ftracks);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithCentFullPID, "Process Full PID reco with centrality", false);

  void processWithRun2CentFullPID(aod::CollisionEvSelRun2Cent const& collision, IdBfFullTracksFullPID const& ftracks);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithRun2CentFullPID, "Process Full PID reco with centrality", false);

  void processWithoutCentFullPID(aod::CollisionEvSel const& collision, IdBfFullTracksFullPID const& ftracks);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithoutCentFullPID, "Process Full PID reco without centrality", false);

  void processWithCentDetectorLevel(aod::CollisionEvSelCent const& collision, IdBfFullTracksDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithCentDetectorLevel, "Process MC detector level with centrality", false);

  void processWithRun2CentDetectorLevel(aod::CollisionEvSelRun2Cent const& collision, IdBfFullTracksDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithRun2CentDetectorLevel, "Process MC detector level with centrality", false);

  void processWithoutCentDetectorLevel(aod::CollisionEvSel const& collision, IdBfFullTracksDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithoutCentDetectorLevel, "Process MC detector level without centrality", false);

  void processWithCentPIDDetectorLevel(aod::CollisionEvSelCent const& collision, IdBfFullTracksPIDDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithCentPIDDetectorLevel, "Process PID MC detector level with centrality", false);

  void processWithRun2CentPIDDetectorLevel(aod::CollisionEvSelRun2Cent const& collision, IdBfFullTracksPIDDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithRun2CentPIDDetectorLevel, "Process PID MC detector level with centrality", false);

  void processWithoutCentPIDDetectorLevel(aod::CollisionEvSel const& collision, IdBfFullTracksPIDDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithoutCentPIDDetectorLevel, "Process PID MC detector level without centrality", false);

  void processWithCentFullPIDDetectorLevel(aod::CollisionEvSelCent const& collision, IdBfFullTracksFullPIDDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithCentFullPIDDetectorLevel, "Process Full PID MC detector level with centrality", false);

  void processWithRun2CentFullPIDDetectorLevel(aod::CollisionEvSelRun2Cent const& collision, IdBfFullTracksFullPIDDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithRun2CentFullPIDDetectorLevel, "Process Full PID MC detector level with centrality", false);

  void processWithoutCentFullPIDDetectorLevel(aod::CollisionEvSel const& collision, IdBfFullTracksFullPIDDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithoutCentFullPIDDetectorLevel, "Process Full PID MC detector level without centrality", false);

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
  PROCESS_SWITCH(IdentifiedBfFilter, processWithCentGeneratorLevel, "Process generated with centrality", false);

  void processWithRun2CentGeneratorLevel(aod::McCollision const& mccollision,
                                         soa::SmallGroups<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>> const& collisions,
                                         aod::McParticles const& mcparticles,
                                         aod::CollisionsEvSelRun2Cent const& allcollisions);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithRun2CentGeneratorLevel, "Process generated with centrality", false);

  void processWithoutCentGeneratorLevel(aod::McCollision const& mccollision,
                                        soa::SmallGroups<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>> const& collisions,
                                        aod::McParticles const& mcparticles,
                                        aod::CollisionsEvSel const& allcollisions);
  PROCESS_SWITCH(IdentifiedBfFilter, processWithoutCentGeneratorLevel, "Process generated without centrality", false);

  void processVertexGenerated(aod::McCollisions const&);
  PROCESS_SWITCH(IdentifiedBfFilter, processVertexGenerated, "Process vertex generator level", false);
};

template <typename CollisionObject, typename TracksObject>
void IdentifiedBfFilter::processReconstructed(CollisionObject const& collision, TracksObject const& ftracks, float tentativecentmult)
{
  using namespace identifiedbffilter;

  LOGF(IDENTIFIEDBFFILTERLOGCOLLISIONS, "IdentifiedBfFilterTask::processReconstructed(). New collision with %d tracks", ftracks.size());

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

void IdentifiedBfFilter::processWithCent(aod::CollisionEvSelCent const& collision, IdBfFullTracks const& ftracks)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void IdentifiedBfFilter::processWithRun2Cent(aod::CollisionEvSelRun2Cent const& collision, IdBfFullTracks const& ftracks)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void IdentifiedBfFilter::processWithoutCent(aod::CollisionEvSel const& collision, IdBfFullTracks const& ftracks)
{
  processReconstructed(collision, ftracks, 50.0);
}

void IdentifiedBfFilter::processWithCentPID(aod::CollisionEvSelCent const& collision, IdBfFullTracksPID const& ftracks)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void IdentifiedBfFilter::processWithRun2CentPID(aod::CollisionEvSelRun2Cent const& collision, IdBfFullTracksPID const& ftracks)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void IdentifiedBfFilter::processWithoutCentPID(aod::CollisionEvSel const& collision, IdBfFullTracksPID const& ftracks)
{
  processReconstructed(collision, ftracks, 50.0);
}

void IdentifiedBfFilter::processWithCentFullPID(aod::CollisionEvSelCent const& collision, IdBfFullTracksFullPID const& ftracks)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void IdentifiedBfFilter::processWithRun2CentFullPID(aod::CollisionEvSelRun2Cent const& collision, IdBfFullTracksFullPID const& ftracks)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void IdentifiedBfFilter::processWithoutCentFullPID(aod::CollisionEvSel const& collision, IdBfFullTracksFullPID const& ftracks)
{
  processReconstructed(collision, ftracks, 50.0);
}

void IdentifiedBfFilter::processWithCentDetectorLevel(aod::CollisionEvSelCent const& collision, IdBfFullTracksDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void IdentifiedBfFilter::processWithRun2CentDetectorLevel(aod::CollisionEvSelRun2Cent const& collision, IdBfFullTracksDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void IdentifiedBfFilter::processWithoutCentDetectorLevel(aod::CollisionEvSel const& collision, IdBfFullTracksDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, 50.0);
}

void IdentifiedBfFilter::processWithCentPIDDetectorLevel(aod::CollisionEvSelCent const& collision, IdBfFullTracksPIDDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void IdentifiedBfFilter::processWithRun2CentPIDDetectorLevel(aod::CollisionEvSelRun2Cent const& collision, IdBfFullTracksPIDDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void IdentifiedBfFilter::processWithoutCentPIDDetectorLevel(aod::CollisionEvSel const& collision, IdBfFullTracksPIDDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, 50.0);
}

void IdentifiedBfFilter::processWithCentFullPIDDetectorLevel(aod::CollisionEvSelCent const& collision, IdBfFullTracksFullPIDDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void IdentifiedBfFilter::processWithRun2CentFullPIDDetectorLevel(aod::CollisionEvSelRun2Cent const& collision, IdBfFullTracksFullPIDDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, getCentMultPercentile(collision));
}

void IdentifiedBfFilter::processWithoutCentFullPIDDetectorLevel(aod::CollisionEvSel const& collision, IdBfFullTracksFullPIDDetLevel const& ftracks, aod::McParticles const&)
{
  processReconstructed(collision, ftracks, 50.0);
}

template <typename CollisionObject, typename ParticlesList>
void IdentifiedBfFilter::processGenerated(CollisionObject const& mccollision, ParticlesList const&, float centormult)
{
  using namespace identifiedbffilter;

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
void IdentifiedBfFilter::processGeneratorLevel(aod::McCollision const& mccollision,
                                               CollisionsGroup const& collisions,
                                               aod::McParticles const& mcparticles,
                                               AllCollisions const& allcollisions,
                                               float defaultcent)
{
  using namespace identifiedbffilter;

  LOGF(IDENTIFIEDBFFILTERLOGCOLLISIONS, "IdentifiedBfFilterTask::processGeneratorLevel(). New generated collision with %d reconstructed collisions and %d particles", collisions.size(), mcparticles.size());

  if (collisions.size() > 1) {
    LOGF(IDENTIFIEDBFFILTERLOGCOLLISIONS, "IdentifiedBfFilterTask::processGeneratorLevel(). Generated collision with more than one reconstructed collisions. Processing only the first accepted for centrality/multiplicity classes extraction");
  }

  bool processed = false;
  for (auto& tmpcollision : collisions) {
    if (tmpcollision.has_mcCollision()) {
      if (tmpcollision.mcCollisionId() == mccollision.globalIndex()) {
        typename AllCollisions::iterator const& collision = allcollisions.iteratorAt(tmpcollision.globalIndex());
        if (IsEvtSelected(collision, defaultcent)) {
          fhTrueVertexZAA->Fill(mccollision.posZ());
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

void IdentifiedBfFilter::processWithCentGeneratorLevel(aod::McCollision const& mccollision,
                                                       soa::SmallGroups<soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>> const& collisions,
                                                       aod::McParticles const& mcparticles,
                                                       aod::CollisionsEvSelCent const& allcollisions)
{
  processGeneratorLevel(mccollision, collisions, mcparticles, allcollisions, 50.0);
}

void IdentifiedBfFilter::processWithRun2CentGeneratorLevel(aod::McCollision const& mccollision,
                                                           soa::SmallGroups<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>> const& collisions,
                                                           aod::McParticles const& mcparticles,
                                                           aod::CollisionsEvSelRun2Cent const& allcollisions)
{
  processGeneratorLevel(mccollision, collisions, mcparticles, allcollisions, 50.0);
}

void IdentifiedBfFilter::processWithoutCentGeneratorLevel(aod::McCollision const& mccollision,
                                                          soa::SmallGroups<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>> const& collisions,
                                                          aod::McParticles const& mcparticles,
                                                          aod::CollisionsEvSel const& allcollisions)
{
  processGeneratorLevel(mccollision, collisions, mcparticles, allcollisions, 50.0);
}

void IdentifiedBfFilter::processVertexGenerated(aod::McCollisions const& mccollisions)
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

struct IdentifiedBfFilterTracks {

  struct : ConfigurableGroup {
    Configurable<std::string> cfgCCDBUrl{"input_ccdburl", "http://ccdb-test.cern.ch:8080", "The CCDB url for the input file"};
    Configurable<std::string> cfgCCDBPathName{"input_ccdbpath", "", "The CCDB path for the input file. Default \"\", i.e. don't load from CCDB"};
    Configurable<std::string> cfgCCDBDate{"input_ccdbdate", "20220307", "The CCDB date for the input file"};
  } cfgcentersinputfile;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  TList* getCCDBInput(const char* ccdbpath, const char* ccdbdate)
  {

    std::tm cfgtm = {};
    std::stringstream ss(ccdbdate);
    ss >> std::get_time(&cfgtm, "%Y%m%d");
    cfgtm.tm_hour = 12;
    int64_t timestamp = std::mktime(&cfgtm) * 1000;

    TList* lst = ccdb->getForTimeStamp<TList>(ccdbpath, timestamp);
    if (lst != nullptr) {
      LOGF(info, "Correctly loaded CCDB input object");
    } else {
      LOGF(error, "CCDB input object could not be loaded");
    }
    return lst;
  }

  Produces<aod::ScannedTracks> scannedtracks;
  Produces<aod::IdentifiedBfCFTracksInfo> tracksinfo;
  Produces<aod::ScannedTrueTracks> scannedgentracks;
  Produces<aod::IdentifiedBfCFGenTracksInfo> gentracksinfo;

  Configurable<bool> cfgFullDerivedData{"fullderiveddata", false, "Produce the full derived data for external storage. Default false"};
  Configurable<int> cfgTrackType{"trktype", 1, "Type of selected tracks: 0 = no selection, 1 = Run2 global tracks FB96, 3 = Run3 tracks, 5 = Run2 TPC only tracks, 7 = Run 3 TPC only tracks. Default 1"};
  Configurable<std::string> cfgSystem{"syst", "PbPb", "System: pp, PbPb, Pbp, pPb, XeXe, ppRun3, PbPbRun3. Default PbPb"};
  Configurable<std::string> cfgDataType{"datatype", "data", "Data type: data, datanoevsel, MC, FastMC, OnTheFlyMC. Default data"};
  Configurable<o2::analysis::DptDptBinningCuts> cfgBinning{"binning",
                                                           {28, -7.0, 7.0, 18, 0.2, 2.0, 16, -0.8, 0.8, 72, 0.5},
                                                           "triplets - nbins, min, max - for z_vtx, pT, eta and phi, binning plus bin fraction of phi origin shift"};
  Configurable<o2::analysis::CheckRangeCfg> cfgTraceDCAOutliers{"trackdcaoutliers", {false, 0.0, 0.0}, "Track the generator level DCAxy outliers: false/true, low dcaxy, up dcaxy. Default {false,0.0,0.0}"};
  Configurable<float> cfgTraceOutOfSpeciesParticles{"trackoutparticles", false, "Track the particles which are not e,mu,pi,K,p: false/true. Default false"};
  Configurable<int> cfgRecoIdMethod{"recoidmethod", 0, "Method for identifying reconstructed tracks: 0 No PID, 1 PID, 2 mcparticle. Default 0"};
  Configurable<o2::analysis::TrackSelectionCfg> cfgTrackSelection{"tracksel", {false, false, 0, 70, 0.8, 2.4, 3.2}, "Track selection: {useit: true/false, ongen: true/false, tpccls, tpcxrws, tpcxrfc, dcaxy, dcaz}. Default {false,0.70.0.8,2.4,3.2}"};
  Configurable<bool> reqTOF{"requireTOF", false, "Require TOF data for PID. Default false"};
  Configurable<bool> onlyTOF{"onlyTOF", false, "Only use TOF data for PID. Default false"};

  Configurable<int> pidEl{"pidEl", -1, "Identify Electron Tracks"};
  Configurable<int> pidPi{"pidPi", -1, "Identify Pion Tracks"};
  Configurable<int> pidKa{"pidKa", -1, "Identify Kaon Tracks"};
  Configurable<int> pidPr{"pidPr", -1, "Identify Proton Tracks"};

  Configurable<float> minPIDSigma{"minpidsigma", -3.0, "Minimum required sigma for PID Acceptance"};
  Configurable<float> maxPIDSigma{"maxpidsigma", 3.0, "Maximum required sigma for PID Acceptance"};

  Configurable<float> minRejectSigma{"minrejectsigma", -1.0, "Minimum required sigma for PID double match rejection"};
  Configurable<float> maxRejectSigma{"maxrejectsigma", 1.0, "Maximum required sigma for PID double match rejection"};

  Configurable<float> tofCut{"TOFCutoff", 0.8, "Momentum under which we don't use TOF PID data"};

  OutputObj<TList> fOutput{"IdentifiedBfFilterTracksInfo", OutputObjHandlingPolicy::AnalysisObject};
  bool checkAmbiguousTracks = false;

  void init(InitContext&)
  {
    LOGF(info, "IdentifiedBfFilterTracks::init()");

    // ccdb info
    ccdb->setURL(cfgcentersinputfile.cfgCCDBUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    LOGF(info, "Initizalized CCDB");

    loadfromccdb = cfgcentersinputfile.cfgCCDBPathName->length() > 0;

    if (ccdblst == nullptr) {
      if (loadfromccdb) {
        LOGF(info, "Loading CCDB Objects");

        ccdblst = getCCDBInput(cfgcentersinputfile.cfgCCDBPathName->c_str(), cfgcentersinputfile.cfgCCDBDate->c_str());
        for (int i = 0; i < kIdBfNoOfSpecies; i++) {
          fhNSigmaCorrection[i] = reinterpret_cast<TH1F*>(ccdblst->FindObject(Form("centerBin_%s", speciesName[i])));
        }
      }
    }
    LOGF(info, "Loaded CCDB Objects");

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

    /* if the system type is not known at this time, we have to put the initialization somewhere else */
    fSystem = getSystemType(cfgSystem);
    fDataType = getDataType(cfgDataType);
    fPDG = TDatabasePDG::Instance();

    /* required ambiguous tracks checks? */
    if (dofilterDetectorLevelWithoutPIDAmbiguous || dofilterDetectorLevelWithPIDAmbiguous || dofilterRecoWithoutPIDAmbiguous || dofilterRecoWithPIDAmbiguous) {
      checkAmbiguousTracks = true;
    }

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
      fhXYB = new TH2F("fHistXYB", "x and y distribution for reconstructed before", 1000, -10.0, 10.0, 1000, -10.0, 10.0);
      fhYZB = new TH2F("fHistYZB", "y and z distribution for reconstructed before", 1000, -10.0, 10.0, 1000, -10.0, 10.0);
      fhXYA = new TH2F("fHistXYA", "x and y distribution for reconstructed after", 1000, -10.0, 10.0, 1000, -10.0, 10.0);
      fhYZA = new TH2F("fHistYZA", "y and z distribution for reconstructed after", 1000, -10.0, 10.0, 1000, -10.0, 10.0);
      fhPB = new TH1F("fHistPB", "p distribution for reconstructed before;p (GeV/c);dN/dp (c/GeV)", 100, 0.0, 15.0);
      fhPtB = new TH1F("fHistPtB", "p_{T} distribution for reconstructed before;p_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhPtPosB = new TH1F("fHistPtPosB", "P_{T} distribution for reconstructed (#plus) before;P_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhPtNegB = new TH1F("fHistPtNegB", "P_{T} distribution for reconstructed (#minus) before;P_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhEtaB = new TH1F("fHistEtaB", "#eta distribution for reconstructed before;#eta;counts", 40, -2.0, 2.0);
      fhEtaA = new TH1F("fHistEtaA", "#eta distribution for reconstructed;#eta;counts", etabins, etalow, etaup);
      fhPhiB = new TH1F("fHistPhiB", "#phi distribution for reconstructed before;#phi;counts", 360, 0.0, constants::math::TwoPI);
      fhdEdxB = new TH2F("fHistdEdxB", "dE/dx vs P before; P (GeV/c); dE/dx (a.u.)", ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
      fhdEdxIPTPCB = new TH2F("fHistdEdxIPTPCB", "dE/dx vs P_{IP} before; P (GeV/c); dE/dx (a.u.)", ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
      fhPhiA = new TH1F("fHistPhiA", "#phi distribution for reconstructed;#phi;counts", 360, 0.0, constants::math::TwoPI);
      fhDCAxyB = new TH1F("DCAxyB", "DCA_{xy} distribution for reconstructed before;DCA_{xy} (cm);counts", 1000, -4.0, 4.0);
      fhDCAxyA = new TH1F("DCAxyA", "DCA_{xy} distribution for reconstructed;DCA_{xy} (cm);counts", 1000, -4., 4.0);
      fhDCAxyzB = new TH2F("DCAxyzB", "DCA_{xy} vs DCA_{z} distribution for reconstructed before;DCA_{xy} (cm); DCA_{z} (cm);counts", 1000, -4.0, 4.0, 1000, -4.0, 4.0);
      fhFineDCAxyA = new TH1F("FineDCAxyA", "DCA_{xy} distribution for reconstructed;DCA_{xy} (cm);counts", 4000, -1.0, 1.0);
      fhDCAzB = new TH1F("DCAzB", "DCA_{z} distribution for reconstructed before;DCA_{z} (cm);counts", 1000, -4.0, 4.0);
      fhDCAzA = new TH1F("DCAzA", "DCA_{z} distribution for reconstructed;DCA_{z} (cm);counts", 1000, -4.0, 4.0);
      fhFineDCAzA = new TH1F("FineDCAzA", "DCA_{z} distribution for reconstructed;DCA_{z} (cm);counts", 4000, -1.0, 1.0);
      fhDCAxyzA = new TH2F("DCAxyzA", "DCA_{xy} vs DCA_{z} distribution for reconstructed;DCA_{xy} (cm); DCA_{z} (cm);counts", 1000, -4.0, 4.0, 1000, -4.0, 4.0);

      fhNClustersB = new TH1I("fHistNClB", "TPC NClusters distribution for reconstructed before;counts", 201, 0, 200);
      fhPhiYB = new TH2F("fHistPhiYB", "#phi vs #eta distribution for reconstructed before;#phi;#eta;counts", 360, 0.0, constants::math::TwoPI, 40, -2.0, 2.0);
      fhPtYB = new TH2F("fHistPtYB", "p_{T} vs #eta distribution for reconstructed before;p_{T} (GeV/c);#eta;counts", 100, 0.0, 15.0, 40, -2.0, 2.0);
      fhChi2B = new TH1F("fHistChi2B", "#chi^{2}/Ncl TPC distribution for reconstructed before;", 100, 0.0, 20.0);
      fhITSNclB = new TH1I("fHistITSNClB", "ITS NClusters distribution for reconstructed before;counts", 21, 0, 20);

      fhNClustersA = new TH1I("fHistNClA", "TPC NClusters distribution for reconstructed after;counts", 201, 0, 200);
      fhPhiYA = new TH2F("fHistPhiYA", "#phi vs #eta distribution for reconstructed after;#phi;#eta;counts", 360, 0.0, constants::math::TwoPI, 40, -2.0, 2.0);
      fhPtYA = new TH2F("fHistPtYA", "p_{T} vs #eta distribution for reconstructed after;p_{T} (GeV/c);#eta;counts", 100, 0.0, 15.0, 40, -2.0, 2.0);
      fhChi2A = new TH1F("fHistChi2A", "#chi^{2}/Ncl TPC distribution for reconstructed after;", 100, 0.0, 20.0);
      fhITSNclA = new TH1I("fHistITSNClA", "ITS NClusters distribution for reconstructed after;counts", 21, 0, 20);

      fhDoublePID = new TH2S("DoublePID", "PIDs for double match;Original Species;Secondary Species", kIdBfNoOfSpecies, 0, kIdBfNoOfSpecies, kIdBfNoOfSpecies, 0, kIdBfNoOfSpecies);

      fhWrongTrackID = new TH1F("WrongTrackId", "Wrong Tracks From Double Track Id distribution in p", ptbins, ptlow, ptup);
      if (checkAmbiguousTracks) {
        /* let's allocate the ambigous tracks tracking histograms*/
        fhAmbiguousTrackType = new TH2D("fHistAmbiguousTracksType", "Ambiguous tracks type vs. multiplicity class;Ambiguous track type;Multiplicity (%);counts", 4, -0.5, 3.5, 101, -0.5, 100.5);
        fhAmbiguousTrackPt = new TH2F("fHistAmbiguousTracksPt", "Ambiguous tracks #it{p}_{T} vs. multiplicity class;#it{p}_{T} (GeV/#it{c});Multiplicity (%);counts", 100, 0.0, 15.0, 101, -0.5, 100.5);
        fhAmbiguityDegree = new TH2F("fHistAmbiguityDegree", "Ambiguity degree vs. multiplicity class;Ambiguity degree;Multiplicity (%);counts", 31, -0.5, 30.5, 101, -0.5, 100.5);
        fhCompatibleCollisionsZVtxRms = new TH2F("fHistCompatibleCollisionsZVtxRms", "Compatible collisions #it{z}_{vtx} RMS;#sigma_{#it{z}_{vtx}};Multiplicity (%);counts", 100, -10.0, 10.0, 101, -0.5, 100.5);
      }

      for (int sp = 0; sp < kIdBfNoOfSpecies; ++sp) {
        fhNSigmaTPC[sp] = new TH2F(TString::Format("fhNSigmaTPC_%s", speciesName[sp]).Data(),
                                   TString::Format("N Sigma from TPC vs P for %s;N #sigma;p (GeV/c)", speciesTitle[sp]).Data(),
                                   48, -6, 6,
                                   ptbins, ptlow, ptup);
        fhNSigmaTOF[sp] = new TH2F(TString::Format("fhNSigmaTOF_%s", speciesName[sp]).Data(),
                                   TString::Format("N Sigma from TOF vs P for %s;N #sigma;p (GeV/c)", speciesTitle[sp]).Data(),
                                   48, -6, 6,
                                   ptbins, ptlow, ptup);
        fhNSigmaCombo[sp] = new TH2F(TString::Format("fhNSigmaCombo_%s", speciesName[sp]).Data(),
                                     TString::Format("N Sigma from Combo vs P for %s;N #sigma;p (GeV/c)", speciesTitle[sp]).Data(),
                                     48, -6, 6,
                                     ptbins, ptlow, ptup);
        fhNSigmaTPC_IdTrks[sp] = new TH2F(TString::Format("fhNSigmaTPC_IdTrks_%s", speciesName[sp]).Data(),
                                          TString::Format("N Sigma from TPC vs P for Identified %s;N #sigma;p (GeV/c)", speciesTitle[sp]).Data(),
                                          48, -6, 6,
                                          ptbins, ptlow, ptup);
      }
      LOGF(info, "Making histos");

      for (int sp = 0; sp < kIdBfNoOfSpecies + 1; ++sp) {
        fhPA[sp] = new TH1F(TString::Format("fHistPA_%s", speciesName[sp]).Data(),
                            TString::Format("p distribution for reconstructed %s;p (GeV/c);dN/dp (c/GeV)", speciesTitle[sp]).Data(),
                            ptbins, ptlow, ptup);
        fhPtA[sp] = new TH1F(TString::Format("fHistPtA_%s", speciesName[sp]),
                             TString::Format("p_{T} distribution for reconstructed %s;p_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                             ptbins, ptlow, ptup);
        fhPtPosA[sp] = new TH1F(TString::Format("fHistPtPosA_%s", speciesName[sp]),
                                TString::Format("P_{T} distribution for reconstructed  %s^{#plus};P_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                ptbins, ptlow, ptup);
        fhPtNegA[sp] = new TH1F(TString::Format("fHistPtNegA_%s", speciesName[sp]),
                                TString::Format("P_{T} distribution for reconstructed  %s^{#minus};P_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                ptbins, ptlow, ptup);
        fhNPosNegA[sp] = new TH2F(TString::Format("fhNPosNegA_%s", speciesName[sp]).Data(),
                                  TString::Format("N(%s^{#plus}) N(%s^{#minus}) distribution for reconstructed;N(%s^{#plus});N(%s^{#minus})", speciesTitle[sp], speciesTitle[sp], speciesTitle[sp], speciesTitle[sp]).Data(),
                                  40, -0.5, 39.5, 40, -0.5, 39.5);
        fhDeltaNA[sp] = new TH1F(TString::Format("fhDeltaNA_%s", speciesName[sp]).Data(),
                                 TString::Format("N(%s^{#plus}) #minus N(%s^{#minus}) distribution for reconstructed;N(%s^{#plus}) #minus N(%s^{#minus})", speciesTitle[sp], speciesTitle[sp], speciesTitle[sp], speciesTitle[sp]).Data(),
                                 79, -39.5, 39.5);
        fhdEdxA[sp] = new TH2F(TString::Format("fhdEdxA_%s", speciesName[sp]).Data(),
                               TString::Format("dE/dx vs P reconstructed %s; P (GeV/c); dE/dx (a.u.)", speciesTitle[sp]).Data(),
                               ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
        fhdEdxIPTPCA[sp] = new TH2F(TString::Format("fhdEdxIPTPCA_%s", speciesName[sp]).Data(),
                                    TString::Format("dE/dx vs P_{IP} reconstructed %s; P (GeV/c); dE/dx (a.u.)", speciesTitle[sp]).Data(),
                                    ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
      }
      fhdEdxA[kIdBfNoOfSpecies + 1] = new TH2F(TString::Format("fhdEdxA_WrongSpecies").Data(),
                                               TString::Format("dE/dx vs P reconstructed Wrong Species; P (GeV/c); dE/dx (a.u.)").Data(),
                                               ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
      fhdEdxIPTPCA[kIdBfNoOfSpecies + 1] = new TH2F(TString::Format("fhdEdxIPTPCA_WrongSpecies").Data(),
                                                    TString::Format("dE/dx vs P_{IP} reconstructed Wrong Species; P (GeV/c); dE/dx (a.u.)").Data(),
                                                    ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
      /* add the hstograms to the output list */
      fOutputList->Add(fhXYB);
      fOutputList->Add(fhYZB);
      fOutputList->Add(fhXYA);
      fOutputList->Add(fhYZA);
      fOutputList->Add(fhNClustersB);
      fOutputList->Add(fhNClustersA);
      fOutputList->Add(fhPhiYB);
      fOutputList->Add(fhPhiYA);
      fOutputList->Add(fhPtYB);
      fOutputList->Add(fhPtYA);
      fOutputList->Add(fhChi2B);
      fOutputList->Add(fhChi2A);
      fOutputList->Add(fhITSNclB);
      fOutputList->Add(fhITSNclA);
      fOutputList->Add(fhPB);
      fOutputList->Add(fhPtB);
      fOutputList->Add(fhPtPosB);
      fOutputList->Add(fhPtNegB);
      fOutputList->Add(fhEtaB);
      fOutputList->Add(fhEtaA);
      fOutputList->Add(fhPhiB);
      fOutputList->Add(fhPhiA);
      fOutputList->Add(fhdEdxB);
      fOutputList->Add(fhdEdxIPTPCB);
      fOutputList->Add(fhDCAxyB);
      fOutputList->Add(fhDCAxyA);
      fOutputList->Add(fhDCAxyzB);
      fOutputList->Add(fhDCAxyzA);
      fOutputList->Add(fhWrongTrackID);
      fOutputList->Add(fhDoublePID);
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

      for (int sp = 0; sp < kIdBfNoOfSpecies; ++sp) {
        fOutputList->Add(fhNSigmaTPC[sp]);
        fOutputList->Add(fhNSigmaTOF[sp]);
        fOutputList->Add(fhNSigmaCombo[sp]);
        fOutputList->Add(fhNSigmaTPC_IdTrks[sp]);
      }

      LOGF(info, "Adding Histos to list");
      for (int sp = 0; sp < kIdBfNoOfSpecies + 1; ++sp) {
        fOutputList->Add(fhPA[sp]);
        fOutputList->Add(fhPtA[sp]);
        fOutputList->Add(fhPtPosA[sp]);
        fOutputList->Add(fhPtNegA[sp]);
        fOutputList->Add(fhNPosNegA[sp]);
        fOutputList->Add(fhDeltaNA[sp]);
        fOutputList->Add(fhdEdxA[sp]);
        fOutputList->Add(fhdEdxIPTPCA[sp]);
      }
      LOGF(info, "Adding Additional Histos to list");
      fOutputList->Add(fhdEdxA[kIdBfNoOfSpecies + 1]);
      fOutputList->Add(fhdEdxIPTPCA[kIdBfNoOfSpecies + 1]);
      LOGF(info, "Added additional histos to list ");
    }

    if ((fDataType != kData) && (fDataType != kDataNoEvtSel)) {
      /* create the true data histograms */

      fhTruePB = new TH1F("fTrueHistPB", "p distribution before (truth);p (GeV/c);dN/dp (c/GeV)", 100, 0.0, 15.0);
      fhTrueCharge = new TH1F("fTrueHistCharge", "Charge distribution before (truth);charge;count", 3, -1.0, 1.0);

      fhTruePhiYB = new TH2F("fTrueHistPhiYB", "#phi vs #eta distribution before (truth);#phi;#eta;counts", 360, 0.0, constants::math::TwoPI, 40, -2.0, 2.0);
      fhTruePtYB = new TH2F("fTrueHistPtYB", "p_{T} vs #eta distribution before (truth);p_{T} (GeV/c);#eta;counts", 100, 0.0, 15.0, 40, -2.0, 2.0);

      fhTruePhiYA = new TH2F("fTrueHistPhiYA", "#phi vs #eta distribution after (truth);#phi;#eta;counts", 360, 0.0, constants::math::TwoPI, 40, -2.0, 2.0);
      fhTruePtYA = new TH2F("fTrueHistPtYA", "p_{T} vs #eta distribution after (truth);p_{T} (GeV/c);#eta;counts", 100, 0.0, 15.0, 40, -2.0, 2.0);

      fhTruePtB = new TH1F("fTrueHistPtB", "p_{T} distribution before (truth);p_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhTruePtPosB = new TH1F("fTrueHistPtPosB", "P_{T} distribution (#plus) before (truth);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhTruePtNegB = new TH1F("fTrueHistPtNegB", "P_{T} distribution (#minus) before (truth);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhTrueEtaB = new TH1F("fTrueHistEtaB", "#eta distribution before (truth);#eta;counts", 40, -2.0, 2.0);
      fhTrueEtaA = new TH1F("fTrueHistEtaA", "#eta distribution (truth);#eta;counts", etabins, etalow, etaup);
      fhTruePhiB = new TH1F("fTrueHistPhiB", "#phi distribution before (truth);#phi;counts", 360, 0.0, constants::math::TwoPI);
      fhTruePhiA = new TH1F("fTrueHistPhiA", "#phi distribution (truth);#phi;counts", 360, 0.0, constants::math::TwoPI);
      fhTrueDCAxyB = new TH1F("TrueDCAxyB", "DCA_{xy} distribution for generated before;DCA_{xy} (cm);counts", 1000, -4.0, 4.0);
      if (traceDCAOutliers.mDoIt) {
        fhTrueDCAxyBid = new TH1F("PDGCodeDCAxyB",
                                  TString::Format("PDG code within %.2f<|DCA_{#it{xy}}|<%.2f; PDG code", traceDCAOutliers.mLowValue, traceDCAOutliers.mUpValue).Data(),
                                  100, 0.5, 100.5);
      }
      fhTrueDCAxyA = new TH1F("TrueDCAxyA", "DCA_{xy} distribution for generated;DCA_{xy};counts (cm)", 1000, -4., 4.0);
      fhTrueDCAzB = new TH1F("TrueDCAzB", "DCA_{z} distribution for generated before;DCA_{z} (cm);counts", 1000, -4.0, 4.0);
      fhTrueDCAzA = new TH1F("TrueDCAzA", "DCA_{z} distribution for generated;DCA_{z} (cm);counts", 1000, -4.0, 4.0);
      fhTrueDCAxyzB = new TH2F("TrueDCAxyzB", "DCA_{xy} vs DCA_{z} distribution for generated before;DCA_{xy} (cm);DCA_{z} (cm);counts", 1000, -4.0, 4.0, 1000, -4.0, 4.0);
      fhTrueDCAxyzA = new TH2F("TrueDCAxyzA", "DCA_{xy} vs DCA_{z} distribution for generated after;DCA_{xy} (cm);DCA_{z} (cm);counts", 1000, -4.0, 4.0, 1000, -4.0, 4.0);

      for (int sp = 0; sp < kIdBfNoOfSpecies + 1; ++sp) {
        fhTruePA[sp] = new TH1F(TString::Format("fTrueHistPA_%s", speciesName[sp]).Data(),
                                TString::Format("p distribution %s (truth);p (GeV/c);dN/dp (c/GeV)", speciesTitle[sp]).Data(),
                                ptbins, ptlow, ptup);
        fhTruePtA[sp] = new TH1F(TString::Format("fTrueHistPtA_%s", speciesName[sp]),
                                 TString::Format("p_{T} distribution %s (truth);p_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                 ptbins, ptlow, ptup);
        fhTruePtPosA[sp] = new TH1F(TString::Format("fTrueHistPtPosA_%s", speciesName[sp]),
                                    TString::Format("P_{T} distribution %s^{#plus} (truth);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                    ptbins, ptlow, ptup);
        fhTruePtNegA[sp] = new TH1F(TString::Format("fTrueHistPtNegA_%s", speciesName[sp]),
                                    TString::Format("P_{T} distribution %s^{#minus} (truth);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                    ptbins, ptlow, ptup);
        fhTrueNPosNegA[sp] = new TH2F(TString::Format("fhTrueNPosNegA_%s", speciesName[sp]).Data(),
                                      TString::Format("N(%s^{#plus}) N(%s^{#minus}) distribution (truth);N(%s^{#plus});N(%s^{#minus})", speciesTitle[sp], speciesTitle[sp], speciesTitle[sp], speciesTitle[sp]).Data(),
                                      40, -0.5, 39.5, 40, -0.5, 39.5);
        fhTrueDeltaNA[sp] = new TH1F(TString::Format("fhTrueDeltaNA_%s", speciesName[sp]).Data(),
                                     TString::Format("N(%s^{#plus}) #minus N(%s^{#minus}) distribution (truth);N(%s^{#plus}) #minus N(%s^{#minus})", speciesTitle[sp], speciesTitle[sp], speciesTitle[sp], speciesTitle[sp]).Data(),
                                     79, -39.5, 39.5);
      }

      /* add the hstograms to the output list */
      fOutputList->Add(fhTruePhiYB);
      fOutputList->Add(fhTruePtYB);
      fOutputList->Add(fhTruePhiYA);
      fOutputList->Add(fhTruePtYA);
      fOutputList->Add(fhTruePB);
      fOutputList->Add(fhTruePtB);
      fOutputList->Add(fhTruePtPosB);
      fOutputList->Add(fhTruePtNegB);
      fOutputList->Add(fhTrueCharge);
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
      fOutputList->Add(fhTrueDCAxyzB);
      fOutputList->Add(fhTrueDCAxyzA);

      for (int sp = 0; sp < kIdBfNoOfSpecies + 1; ++sp) {
        fOutputList->Add(fhTruePA[sp]);
        fOutputList->Add(fhTruePtA[sp]);
        fOutputList->Add(fhTruePtPosA[sp]);
        fOutputList->Add(fhTruePtNegA[sp]);
        fOutputList->Add(fhTrueNPosNegA[sp]);
        fOutputList->Add(fhTrueDeltaNA[sp]);
      }
    }
    /* initialize access to the CCDB */
  }

  template <typename TrackObject>
  inline MatchRecoGenSpecies IdentifyTrack(TrackObject const& track);
  template <typename TrackObject>
  int8_t AcceptTrack(TrackObject const& track);
  template <typename ParticleObject, typename MCCollisionObject>
  int8_t AcceptParticle(ParticleObject& particle, MCCollisionObject const& mccollision);
  template <typename CollisionObjects, typename TrackObject>
  int8_t selectTrackAmbiguousCheck(CollisionObjects const& collisions, TrackObject const& track);
  template <typename ParticleObject>
  inline MatchRecoGenSpecies IdentifyParticle(ParticleObject const& particle);
  template <typename TrackObject>
  void fillTrackHistosBeforeSelection(TrackObject const& track);
  template <typename TrackObject>
  void fillTrackHistosAfterSelection(TrackObject const& track, MatchRecoGenSpecies sp);
  template <typename ParticleObject, typename MCCollisionObject>
  void fillParticleHistosBeforeSelection(ParticleObject const& particle,
                                         MCCollisionObject const& collision,
                                         float charge);
  template <typename ParticleObject, typename MCCollisionObject>
  void fillParticleHistosAfterSelection(ParticleObject const& particle,
                                        MCCollisionObject const& collision,
                                        float charge,
                                        MatchRecoGenSpecies sp);

  /* TODO: as it is now when the derived data is stored (fullDerivedData = true) */
  /* the collision index stored with the track is wrong. This has to be fixed    */
  template <typename passedtracks>
  void filterTracks(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo> const& collisions,
                    passedtracks const& tracks)
  {
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
      if (track.has_collision() && (track.template collision_as<soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>>()).collisionaccepted()) {
        pid = selectTrackAmbiguousCheck(collisions, track);
        if (!(pid < 0)) {
          naccepted++;
          /* update charged multiplicities */
          if (pid % 2 == 0) {
            trkMultPos[kIdBfCharged]++;
          }
          if (pid % 2 == 1) {
            trkMultNeg[kIdBfCharged]++;
          }
          if (fullDerivedData) {
            LOGF(fatal, "Stored derived data not prepared for saving the proper new collision id");
            scannedtracks((track.template collision_as<soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>>()).globalIndex(), pid, track.pt(), track.eta(), track.phi());
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
    LOGF(IDENTIFIEDBFFILTERLOGCOLLISIONS,
         "Processed %d accepted collisions out of a total of %d with  %d accepted tracks out of a "
         "total of %d",
         ncollaccepted,
         collisions.size(),
         naccepted,
         tracks.size());
  }

  /* TODO: for the time being the full derived data is still not supported  */
  /* for doing that we need to get the index of the associated mc collision */
  void filterParticles(soa::Join<aod::McCollisions, aod::IdentifiedBfCFGenCollisionsInfo> const& gencollisions, aod::McParticles const& particles)
  {
    using namespace identifiedbffilter;
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
      int8_t pid = -1;
      if (particle.isPhysicalPrimary()) {
        TParticlePDG* pdgpart = fPDG->GetParticle(particle.pdgCode());
        float charge = 0;
        if (pdgpart != nullptr) {
          charge = getCharge(pdgpart->Charge());
          // print charge
        }
        fhTrueCharge->Fill(charge);
        if (charge != 0) {
          if (particle.has_mcCollision() && (particle.template mcCollision_as<soa::Join<aod::McCollisions, aod::IdentifiedBfCFGenCollisionsInfo>>()).collisionaccepted()) {
            auto mccollision = particle.template mcCollision_as<soa::Join<aod::McCollisions, aod::IdentifiedBfCFGenCollisionsInfo>>();
            /* before particle selection */
            fillParticleHistosBeforeSelection(particle, mccollision, charge);

            /* track selection */
            pid = AcceptParticle(particle, mccollision);
            if (!(pid < 0)) { // if PID isn't negative
              acceptedparticles++;
            }
          }
        } else {
          if ((particle.mcCollisionId() == 0) && traceCollId0) {
            LOGF(IDENTIFIEDBFFILTERLOGTRACKS, "Particle %d with fractional charge or equal to zero", particle.globalIndex());
          }
        }

      } else {
        if ((particle.mcCollisionId() == 0) && traceCollId0) {
          LOGF(IDENTIFIEDBFFILTERLOGTRACKS, "Particle %d not Physical Primary", particle.globalIndex());
        }
      }
      if (!fullDerivedData) {
        gentracksinfo(pid);
      }
    }
    LOGF(info,
         "Processed %d accepted generated collisions out of a total of %d with  %d accepted particles out of a "
         "total of %d",
         acceptedcollisions,
         gencollisions.size(),
         acceptedparticles,
         particles.size());
  }

  void filterRecoWithPID(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksPID const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterRecoWithPID, "Not stored derived data track filtering", false)

  void filterRecoWithPIDAmbiguous(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksPIDAmbiguous const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterRecoWithPIDAmbiguous, "Not stored derived data track filtering with ambiguous tracks check", false)

  void filterDetectorLevelWithPID(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksPIDDetLevel const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterDetectorLevelWithPID, "Not stored derived data detector level track filtering", false)

  void filterDetectorLevelWithPIDAmbiguous(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksPIDDetLevelAmbiguous const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterDetectorLevelWithPIDAmbiguous, "Not stored derived data detector level track filtering with ambiguous tracks check", false)

  void filterRecoWithFullPID(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksFullPID const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterRecoWithFullPID, "Not stored derived data track filtering", false)

  void filterRecoWithFullPIDAmbiguous(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksFullPIDAmbiguous const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterRecoWithFullPIDAmbiguous, "Not stored derived data track filtering with ambiguous tracks check", false)

  void filterDetectorLevelWithFullPID(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksFullPIDDetLevel const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterDetectorLevelWithFullPID, "Not stored derived data detector level track filtering", false)

  void filterDetectorLevelWithFullPIDAmbiguous(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksFullPIDDetLevelAmbiguous const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterDetectorLevelWithFullPIDAmbiguous, "Not stored derived data detector level track filtering with ambiguous tracks check", false)

  void filterRecoWithoutPID(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo> const& collisions, IdBfFullTracks const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterRecoWithoutPID, "Track filtering without PID information", true)

  void filterRecoWithoutPIDAmbiguous(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo> const& collisions, IdBfFullTracksAmbiguous const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterRecoWithoutPIDAmbiguous, "Track filtering without PID information with ambiguous tracks check", false)

  void filterDetectorLevelWithoutPID(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo> const& collisions, IdBfFullTracksDetLevel const& tracks, aod::McParticles const&)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterDetectorLevelWithoutPID, "Detector level track filtering without PID information", false)

  void filterDetectorLevelWithoutPIDAmbiguous(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo> const& collisions, IdBfFullTracksDetLevelAmbiguous const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterDetectorLevelWithoutPIDAmbiguous, "Detector level track filtering without PID information with ambiguous tracks check", false)

  void filterGenerated(soa::Join<aod::McCollisions, aod::IdentifiedBfCFGenCollisionsInfo> const& gencollisions, aod::McParticles const& particles)
  {
    filterParticles(gencollisions, particles);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterGenerated, "Generated particles filtering", true)
};

template <typename ParticleObject>
inline MatchRecoGenSpecies IdentifiedBfFilterTracks::IdentifyParticle(ParticleObject const& particle)
{
  using namespace identifiedbffilter;

  int pdgcode = fabs(particle.pdgCode());

  switch (pdgcode) {
    case pdgcodeEl:
      return kIdBfElectron;
      break;

    case pdgcodePi:
      return kIdBfPion;
      break;
    case pdgcodeKa:
      return kIdBfKaon;
      break;
    case pdgcodePr:
      return kIdBfProton;
      break;

    default:
      if (traceOutOfSpeciesParticles) {
        LOGF(info, "Wrong particle passed selection cuts. PDG code: %d", pdgcode);
      }
      return kWrongSpecies;
      break;
  }
}

template <typename TrackObject>
void fillNSigmaHistos(TrackObject const& track)
{

  float actualTPCNSigmaEl = track.tpcNSigmaEl();
  float actualTPCNSigmaPi = track.tpcNSigmaPi();
  float actualTPCNSigmaKa = track.tpcNSigmaKa();
  float actualTPCNSigmaPr = track.tpcNSigmaPr();

  if (loadfromccdb) {
    actualTPCNSigmaEl = actualTPCNSigmaEl - fhNSigmaCorrection[kIdBfElectron]->GetBinContent(fhNSigmaCorrection[kIdBfElectron]->FindBin(track.tpcInnerParam()));
    actualTPCNSigmaPi = actualTPCNSigmaPi - fhNSigmaCorrection[kIdBfPion]->GetBinContent(fhNSigmaCorrection[kIdBfPion]->FindBin(track.tpcInnerParam()));
    actualTPCNSigmaKa = actualTPCNSigmaKa - fhNSigmaCorrection[kIdBfKaon]->GetBinContent(fhNSigmaCorrection[kIdBfKaon]->FindBin(track.tpcInnerParam()));
    actualTPCNSigmaPr = actualTPCNSigmaPr - fhNSigmaCorrection[kIdBfProton]->GetBinContent(fhNSigmaCorrection[kIdBfProton]->FindBin(track.tpcInnerParam()));
  }

  fhNSigmaTPC[kIdBfElectron]->Fill(actualTPCNSigmaEl, track.tpcInnerParam());
  fhNSigmaTPC[kIdBfPion]->Fill(actualTPCNSigmaPi, track.tpcInnerParam());
  fhNSigmaTPC[kIdBfKaon]->Fill(actualTPCNSigmaKa, track.tpcInnerParam());
  fhNSigmaTPC[kIdBfProton]->Fill(actualTPCNSigmaPr, track.tpcInnerParam());

  fhNSigmaTOF[kIdBfElectron]->Fill(track.tofNSigmaEl(), track.tpcInnerParam());
  fhNSigmaTOF[kIdBfPion]->Fill(track.tofNSigmaPi(), track.tpcInnerParam());
  fhNSigmaTOF[kIdBfKaon]->Fill(track.tofNSigmaKa(), track.tpcInnerParam());
  fhNSigmaTOF[kIdBfProton]->Fill(track.tofNSigmaPr(), track.tpcInnerParam());

  fhNSigmaCombo[kIdBfElectron]->Fill(sqrtf(track.tofNSigmaEl() * track.tofNSigmaEl() + actualTPCNSigmaEl * actualTPCNSigmaEl), track.tpcInnerParam());
  fhNSigmaCombo[kIdBfPion]->Fill(sqrtf(track.tofNSigmaPi() * track.tofNSigmaPi() + actualTPCNSigmaPi * actualTPCNSigmaPi), track.tpcInnerParam());
  fhNSigmaCombo[kIdBfKaon]->Fill(sqrtf(track.tofNSigmaKa() * track.tofNSigmaKa() + actualTPCNSigmaKa * actualTPCNSigmaKa), track.tpcInnerParam());
  fhNSigmaCombo[kIdBfProton]->Fill(sqrtf(track.tofNSigmaPr() * track.tofNSigmaPr() + actualTPCNSigmaPr * actualTPCNSigmaPr), track.tpcInnerParam());
}

template <typename TrackObject>
inline MatchRecoGenSpecies IdentifiedBfFilterTracks::IdentifyTrack(TrackObject const& track)
{
  using namespace o2::analysis::identifiedbffilter;

  fillNSigmaHistos(track);

  float actualTPCNSigmaEl = track.tpcNSigmaEl();
  float actualTPCNSigmaPi = track.tpcNSigmaPi();
  float actualTPCNSigmaKa = track.tpcNSigmaKa();
  float actualTPCNSigmaPr = track.tpcNSigmaPr();

  float nsigmas[kIdBfNoOfSpecies];

  if (loadfromccdb) {
    actualTPCNSigmaEl = actualTPCNSigmaEl - fhNSigmaCorrection[kIdBfElectron]->GetBinContent(fhNSigmaCorrection[kIdBfElectron]->FindBin(track.tpcInnerParam()));
    actualTPCNSigmaPi = actualTPCNSigmaPi - fhNSigmaCorrection[kIdBfPion]->GetBinContent(fhNSigmaCorrection[kIdBfPion]->FindBin(track.tpcInnerParam()));
    actualTPCNSigmaKa = actualTPCNSigmaKa - fhNSigmaCorrection[kIdBfKaon]->GetBinContent(fhNSigmaCorrection[kIdBfKaon]->FindBin(track.tpcInnerParam()));
    actualTPCNSigmaPr = actualTPCNSigmaPr - fhNSigmaCorrection[kIdBfProton]->GetBinContent(fhNSigmaCorrection[kIdBfProton]->FindBin(track.tpcInnerParam()));
  }

  if (track.tpcInnerParam() < tofCut && !reqTOF && !onlyTOF) {

    nsigmas[kIdBfElectron] = actualTPCNSigmaEl;
    nsigmas[kIdBfPion] = actualTPCNSigmaPi;
    nsigmas[kIdBfKaon] = actualTPCNSigmaKa;
    nsigmas[kIdBfProton] = actualTPCNSigmaPr;

  } else {
    /* introduce require TOF flag */
    if (track.hasTOF() && !onlyTOF) {
      nsigmas[kIdBfElectron] = sqrtf(actualTPCNSigmaEl * actualTPCNSigmaEl + track.tofNSigmaEl() * track.tofNSigmaEl());
      nsigmas[kIdBfPion] = sqrtf(actualTPCNSigmaPi * actualTPCNSigmaPi + track.tofNSigmaPi() * track.tofNSigmaPi());
      nsigmas[kIdBfKaon] = sqrtf(actualTPCNSigmaKa * actualTPCNSigmaKa + track.tofNSigmaKa() * track.tofNSigmaKa());
      nsigmas[kIdBfProton] = sqrtf(actualTPCNSigmaPr * actualTPCNSigmaPr + track.tofNSigmaPr() * track.tofNSigmaPr());

    } else if (!reqTOF || !onlyTOF) {
      nsigmas[kIdBfElectron] = actualTPCNSigmaEl;
      nsigmas[kIdBfPion] = actualTPCNSigmaPi;
      nsigmas[kIdBfKaon] = actualTPCNSigmaKa;
      nsigmas[kIdBfProton] = actualTPCNSigmaPr;

    } else if (track.hasTOF() && onlyTOF) {
      nsigmas[kIdBfElectron] = track.tofNSigmaEl();
      nsigmas[kIdBfPion] = track.tofNSigmaPi();
      nsigmas[kIdBfKaon] = track.tofNSigmaKa();
      nsigmas[kIdBfProton] = track.tofNSigmaPr();
    } else {
      return kWrongSpecies;
    }
  }

  if (!pidEl) {
    nsigmas[kIdBfElectron] = 999.0f;
  }
  if (!pidPi) {
    nsigmas[kIdBfPion] = 999.0f;
  }
  if (!pidKa) {
    nsigmas[kIdBfKaon] = 999.0f;
  }
  if (!pidPr) {
    nsigmas[kIdBfProton] = 999.0f;
  }

  float min_nsigma = 999.0f;
  MatchRecoGenSpecies sp_min_nsigma = kWrongSpecies;
  for (int sp = 0; sp < kIdBfNoOfSpecies; ++sp) {
    if (fabs(nsigmas[sp]) < fabs(min_nsigma)) { // Check if species nsigma is less than current nsigma
      min_nsigma = nsigmas[sp];                // If yes, set species nsigma to current nsigma
      sp_min_nsigma = MatchRecoGenSpecies(sp); // set current species sp number to current sp
    }
  }
  bool doublematch = false;
  MatchRecoGenSpecies spDouble = kWrongSpecies;
  if (min_nsigma < maxPIDSigma && min_nsigma > minPIDSigma) {         // Check that current nsigma is in accpetance range
    for (int sp = 0; (sp < kIdBfNoOfSpecies) && !doublematch; ++sp) { // iterate over all species while there's no double match and we're in the list
      if (sp != sp_min_nsigma) {                                      // for species not current minimum nsigma species
        if (nsigmas[sp] < maxRejectSigma && nsigmas[sp] > minRejectSigma) { // If secondary species is in rejection range
          doublematch = true;                                         // Set double match true
          spDouble = MatchRecoGenSpecies(sp);
        }
      }
    }
    if (doublematch) { // if double match true
      fhWrongTrackID->Fill(track.p());
      fhdEdxA[kIdBfNoOfSpecies]->Fill(track.p(), track.tpcSignal());
      fhdEdxIPTPCA[kIdBfNoOfSpecies]->Fill(track.tpcInnerParam(), track.tpcSignal());
      fhDoublePID->Fill(sp_min_nsigma, spDouble);
      return kWrongSpecies; // Return wrong species value
    } else {
      if (sp_min_nsigma == 0) {
        fhNSigmaTPC_IdTrks[sp_min_nsigma]->Fill(actualTPCNSigmaEl, track.tpcInnerParam());
      }
      if (sp_min_nsigma == 1) {
        fhNSigmaTPC_IdTrks[sp_min_nsigma]->Fill(actualTPCNSigmaPi, track.tpcInnerParam());
      }
      if (sp_min_nsigma == 2) {
        fhNSigmaTPC_IdTrks[sp_min_nsigma]->Fill(actualTPCNSigmaKa, track.tpcInnerParam());
      }
      if (sp_min_nsigma == 3) {
        fhNSigmaTPC_IdTrks[sp_min_nsigma]->Fill(actualTPCNSigmaPr, track.tpcInnerParam());
      }

      return sp_min_nsigma;
    }
  } else {
    return kWrongSpecies;
  }
}

/// \brief Accepts or not the passed track
/// \param track the track of interest
/// \return the internal track id, -1 if not accepted
/// TODO: the PID implementation
/// For the time being we keep the convention
/// - positive track pid even
/// - negative track pid odd
template <typename TrackObject>
inline int8_t IdentifiedBfFilterTracks::AcceptTrack(TrackObject const& track)
{
  fillTrackHistosBeforeSelection(track); // <Fill "before selection" histo

  /* TODO: incorporate a mask in the scanned tracks table for the rejecting track reason */
  if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
    if (track.mcParticleId() < 0) {
      return -1;
    }
  }

  if (matchTrackType(track)) {
    if (ptlow < track.pt() && track.pt() < ptup && etalow < track.eta() && track.eta() < etaup) {
      fillTrackHistosAfterSelection(track, kIdBfCharged);
      MatchRecoGenSpecies sp = kWrongSpecies;
      if (recoIdMethod == 0) {
        sp = kIdBfCharged;
      } else if (recoIdMethod == 1) {

        if constexpr (framework::has_type_v<aod::pidtpc_tiny::TPCNSigmaStorePi, typename TrackObject::all_columns> || framework::has_type_v<aod::pidtpc::TPCNSigmaPi, typename TrackObject::all_columns>) {
          sp = IdentifyTrack(track);
        } else {
          LOGF(fatal, "Track identification required but PID information not present");
        }
      } else if (recoIdMethod == 2) {
        if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
          sp = IdentifyParticle(track.template mcParticle_as<aod::McParticles>());
        } else {
          LOGF(fatal, "Track identification required from MC particle but MC information not present");
        }
      }
      if (sp == kWrongSpecies) {
        return -1;
      }
      if (!(sp < 0)) {
        fillTrackHistosAfterSelection(track, sp); //<Fill accepted track histo with PID
        if (track.sign() > 0) {                   // if positive
          trkMultPos[sp]++; //<< Update Particle Multiplicity
          return speciesChargeValue1[sp];
        }
        if (track.sign() < 0) { // if negative
          trkMultNeg[sp]++; //<< Update Particle Multiplicity
          return speciesChargeValue1[sp] + 1;
        }
      }
    }
  }
  return -1;
}

/// \brief Accepts or not the passed generated particle
/// \param track the particle of interest
/// \return `true` if the particle is accepted, `false` otherwise
template <typename ParticleObject, typename MCCollisionObject>
inline int8_t IdentifiedBfFilterTracks::AcceptParticle(ParticleObject& particle, MCCollisionObject const& mccollision)
{
  /* overall momentum cut */
  if (!(overallminp < particle.p())) {
    return kWrongSpecies;
  }
  TParticlePDG* pdgpart = fPDG->GetParticle(particle.pdgCode());
  float charge = 0;
  if (pdgpart != nullptr) {
    charge = getCharge(pdgpart->Charge());
  }

  if (particle.isPhysicalPrimary() && fabs(charge) > 0.0) {
    if ((particle.mcCollisionId() == 0) && traceCollId0) {
      LOGF(info, "Particle %d passed isPhysicalPrimary", particle.globalIndex());
    }

    if (ptlow < particle.pt() && particle.pt() < ptup && etalow < particle.eta() && particle.eta() < etaup) {
      MatchRecoGenSpecies sp = IdentifyParticle(particle);
      if (sp != kWrongSpecies) {
        if (sp != kIdBfCharged) {
          /* fill the charged particle histograms */
          fillParticleHistosAfterSelection(particle, mccollision, charge, kIdBfCharged);
          /* update charged multiplicities */
          if (charge == 1) {
            partMultPos[kIdBfCharged]++;
          } else if (charge == -1) {
            partMultNeg[kIdBfCharged]++;
          }
        }
        /* fill the species  histograms */
        fillParticleHistosAfterSelection(particle, mccollision, charge, sp);
        /* update species multiplicities */
        if (charge == 1) {
          partMultPos[sp]++;
        } else if (charge == -1) {
          partMultNeg[sp]++;
        }
      }
      if (charge == 1) {
        return speciesChargeValue1[sp];

      } else if (charge == -1) {
        return speciesChargeValue1[sp] + 1;
      }
    }
  } else {
    if ((particle.mcCollisionId() == 0) && traceCollId0) {
      LOGF(info, "Particle %d NOT passed isPhysicalPrimary", particle.globalIndex());
    }
  }
  return kWrongSpecies;
}

template <typename CollisionObjects, typename TrackObject>
int8_t IdentifiedBfFilterTracks::selectTrackAmbiguousCheck(CollisionObjects const& collisions, TrackObject const& track)
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

  float multiplicityclass = (track.template collision_as<soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>>()).centmult();
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
    return AcceptTrack(track);
  }
}

template <typename TrackObject>
void IdentifiedBfFilterTracks::fillTrackHistosBeforeSelection(TrackObject const& track)
{
  fhXYB->Fill(track.x(), track.y());
  fhYZB->Fill(track.y(), track.z());
  fhNClustersB->Fill(track.tpcNClsFound());
  fhPhiYB->Fill(track.phi(), track.eta());
  fhPtYB->Fill(track.pt(), track.eta());
  fhChi2B->Fill(track.tpcChi2NCl());
  fhITSNclB->Fill(track.itsNCls());

  fhPB->Fill(track.p());
  fhPtB->Fill(track.pt());
  fhEtaB->Fill(track.eta());
  fhPhiB->Fill(track.phi());
  fhdEdxB->Fill(track.p(), track.tpcSignal());
  fhdEdxIPTPCB->Fill(track.tpcInnerParam(), track.tpcSignal());
  if (track.sign() > 0) {
    fhPtPosB->Fill(track.pt());
  } else {
    fhPtNegB->Fill(track.pt());
  }
  fhDCAxyB->Fill(track.dcaXY());
  fhDCAzB->Fill(track.dcaZ());
  fhDCAxyzB->Fill(track.dcaXY(), track.dcaZ());
}

template <typename TrackObject>
void IdentifiedBfFilterTracks::fillTrackHistosAfterSelection(TrackObject const& track, MatchRecoGenSpecies sp)
{
  /* the charged species should have been called first so avoid double counting */
  if (sp == kIdBfCharged) {
    fhEtaA->Fill(track.eta());
    fhPhiA->Fill(track.phi());
    fhXYA->Fill(track.x(), track.y());
    fhYZA->Fill(track.y(), track.z());
    fhNClustersA->Fill(track.tpcNClsFound());
    fhPhiYA->Fill(track.phi(), track.eta());
    fhPtYA->Fill(track.pt(), track.eta());
    fhChi2A->Fill(track.tpcChi2NCl());
    fhITSNclA->Fill(track.itsNCls());
    fhDCAxyA->Fill(track.dcaXY());
    fhDCAzA->Fill(track.dcaZ());
    fhDCAxyzA->Fill(track.dcaXY(), track.dcaZ());

    if (track.dcaXY() < 1.0) {
      fhFineDCAxyA->Fill(track.dcaXY());
    }
    if (track.dcaZ() < 1.0) {
      fhFineDCAzA->Fill(track.dcaZ());
    }
  }
  fhPA[sp]->Fill(track.p());
  fhPtA[sp]->Fill(track.pt());
  fhdEdxA[sp]->Fill(track.p(), track.tpcSignal());
  fhdEdxIPTPCA[sp]->Fill(track.tpcInnerParam(), track.tpcSignal());
  if (track.sign() > 0) {
    fhPtPosA[sp]->Fill(track.pt());
  } else {
    fhPtNegA[sp]->Fill(track.pt());
  }
}

template <typename ParticleObject, typename MCCollisionObject>
void IdentifiedBfFilterTracks::fillParticleHistosBeforeSelection(ParticleObject const& particle, MCCollisionObject const& collision, float charge)
{
  fhTruePhiYB->Fill(particle.phi(), particle.eta());
  fhTruePtYB->Fill(particle.pt(), particle.eta());
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
  fhTrueDCAxyzB->Fill(TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                                  (particle.vy() - collision.posY()) * (particle.vy() - collision.posY())),
                      (particle.vz() - collision.posZ()));
  fhTrueDCAzB->Fill((particle.vz() - collision.posZ()));
}

template <typename ParticleObject, typename MCCollisionObject>
void IdentifiedBfFilterTracks::fillParticleHistosAfterSelection(ParticleObject const& particle, MCCollisionObject const& collision, float charge, MatchRecoGenSpecies sp)
{
  /* the charged species should have been called first so avoid double counting */
  if (sp == kIdBfCharged) {
    fhTruePhiYA->Fill(particle.phi(), particle.eta());
    fhTruePtYA->Fill(particle.pt(), particle.eta());
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
    fhTrueDCAxyzA->Fill(TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                                    (particle.vy() - collision.posY()) * (particle.vy() - collision.posY())),
                        (particle.vz() - collision.posZ()));
  }
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
  WorkflowSpec workflow{adaptAnalysisTask<IdentifiedBfFilter>(cfgc,
                                                              SetDefaultProcesses{
                                                                {{"processWithoutCent", true},
                                                                 {"processWithoutCentGeneratorLevel", true}}}),
                        adaptAnalysisTask<IdentifiedBfFilterTracks>(cfgc)};
  return workflow;
}
