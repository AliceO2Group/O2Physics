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

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "PWGCF/DataModel/DptDptFiltered.h"
#include "Common/DataModel/TrackSelectionTables.h"
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
using DptDptFullTracksPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
using DptDptFullTracksPIDDetLevel = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TracksDCA, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
using DptDptFullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
using DptDptFullTracksDetLevel = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TracksDCA>;

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
// The DptDptFilter output objects
//============================================================================================
TH1F* fhCentMultB = nullptr;
TH1F* fhCentMultA = nullptr;
TH1F* fhVertexZB = nullptr;
TH1F* fhVertexZA = nullptr;
TH1F* fhMultB = nullptr;
TH1F* fhMultA = nullptr;
TH1F* fhPB = nullptr;
TH1F* fhPA[kDptDptNoOfSpecies] = {nullptr};
TH1F* fhPtB = nullptr;
TH1F* fhPtA[kDptDptNoOfSpecies] = {nullptr};
TH1F* fhPtPosB = nullptr;
TH1F* fhPtPosA[kDptDptNoOfSpecies] = {nullptr};
TH1F* fhPtNegB = nullptr;
TH1F* fhPtNegA[kDptDptNoOfSpecies] = {nullptr};
TH2F* fhNPosNegA[kDptDptNoOfSpecies] = {nullptr};
TH1F* fhDeltaNA[kDptDptNoOfSpecies] = {nullptr};

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

TH1F* fhTrueCentMultB = nullptr;
TH1F* fhTrueCentMultA = nullptr;
TH1F* fhTrueVertexZB = nullptr;
TH1F* fhTrueVertexZA = nullptr;
TH1F* fhTrueVertexZAA = nullptr;
TH1F* fhTruePB = nullptr;
TH1F* fhTruePA[kDptDptNoOfSpecies] = {nullptr};
TH1F* fhTruePtB = nullptr;
TH1F* fhTruePtA[kDptDptNoOfSpecies] = {nullptr};
TH1F* fhTruePtPosB = nullptr;
TH1F* fhTruePtPosA[kDptDptNoOfSpecies] = {nullptr};
TH1F* fhTruePtNegB = nullptr;
TH1F* fhTruePtNegA[kDptDptNoOfSpecies] = {nullptr};
TH2F* fhTrueNPosNegA[kDptDptNoOfSpecies] = {nullptr};
TH1F* fhTrueDeltaNA[kDptDptNoOfSpecies] = {nullptr};

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
int trkMultPos[kDptDptNoOfSpecies];  // multiplicity of positive tracks
int trkMultNeg[kDptDptNoOfSpecies];  // multiplicity of negative tracks
int partMultPos[kDptDptNoOfSpecies]; // multiplicity of positive particles
int partMultNeg[kDptDptNoOfSpecies]; // multiplicity of negative particles
} // namespace o2::analysis::dptdptfilter

using namespace dptdptfilter;

struct DptDptFilter {
  Configurable<bool> cfgFullDerivedData{"fullderiveddata", false, "Produce the full derived data for external storage. Default false"};
  Configurable<std::string> cfgCentMultEstimator{"centmultestimator", "V0M", "Centrality/multiplicity estimator detector: V0M,CL0,CL1,FV0A,FT0M,FT0A,FT0C,NTPV,NOCM: none. Default V0M"};
  Configurable<std::string> cfgSystem{"syst", "PbPb", "System: pp, PbPb, Pbp, pPb, XeXe, ppRun3, PbPbRun3. Default PbPb"};
  Configurable<std::string> cfgDataType{"datatype", "data", "Data type: data, datanoevsel, MC, FastMC, OnTheFlyMC. Default data"};
  Configurable<std::string> cfgTriggSel{"triggsel", "MB", "Trigger selection: MB, None. Default MB"};
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

struct DptDptFilterTracks {
  Produces<aod::ScannedTracks> scannedtracks;
  Produces<aod::DptDptCFTracksInfo> tracksinfo;
  Produces<aod::ScannedTrueTracks> scannedgentracks;
  Produces<aod::DptDptCFGenTracksInfo> gentracksinfo;

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

  OutputObj<TList> fOutput{"DptDptFilterTracksInfo", OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    LOGF(info, "DptDptFilterTracks::init()");

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

      for (int sp = 0; sp < kDptDptNoOfSpecies; ++sp) {
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

      for (int sp = 0; sp < kDptDptNoOfSpecies; ++sp) {
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
      fhTruePB = new TH1F("fTrueHistPB", "p distribution before (truth);p (GeV/c);dN/dp (c/GeV)", 100, 0.0, 15.0);
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

      for (int sp = 0; sp < kDptDptNoOfSpecies; ++sp) {
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

      for (int sp = 0; sp < kDptDptNoOfSpecies; ++sp) {
        fOutputList->Add(fhTruePA[sp]);
        fOutputList->Add(fhTruePtA[sp]);
        fOutputList->Add(fhTruePtPosA[sp]);
        fOutputList->Add(fhTruePtNegA[sp]);
        fOutputList->Add(fhTrueNPosNegA[sp]);
        fOutputList->Add(fhTrueDeltaNA[sp]);
      }
    }
  }

  template <typename TrackObject>
  inline MatchRecoGenSpecies IdentifyTrack(TrackObject const& track);
  template <typename TrackObject>
  MatchRecoGenSpecies trackIdentification(TrackObject const& track);
  template <typename TrackObject>
  int8_t selectTrack(TrackObject const& track);
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
  void filterTracks(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const& collisions,
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
      if (track.has_collision() && (track.template collision_as<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>()).collisionaccepted()) {
        pid = selectTrack(track);
        if (!(pid < 0)) {
          naccepted++;
          if (fullDerivedData) {
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
      float charge = 0.0;
      TParticlePDG* pdgparticle = fPDG->GetParticle(particle.pdgCode());
      if (pdgparticle != nullptr) {
        charge = (pdgparticle->Charge() / 3 >= 1) ? 1.0 : ((pdgparticle->Charge() / 3 <= -1) ? -1.0 : 0.0);
      }

      int8_t pid = -1;

      if (charge != 0) {
        if (particle.has_mcCollision() && (particle.template mcCollision_as<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>()).collisionaccepted()) {
          auto mccollision = particle.template mcCollision_as<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>();
          /* before particle selection */
          fillParticleHistosBeforeSelection(particle, mccollision, charge);

          /* track selection */
          /* TODO: at some point the pid has to be substituted by the identified species */
          pid = AcceptParticle(particle, mccollision);
          if (!(pid < 0)) {
            /* the particle has been accepted */
            /* let's identify the particle */
            /* TODO: probably this needs to go to AcceptParticle */
            MatchRecoGenSpecies sp = IdentifyParticle(particle);
            if (sp != kWrongSpecies) {
              if (sp != kDptDptCharged) {
                /* fill the charged particle histograms */
                fillParticleHistosAfterSelection(particle, mccollision, charge, kDptDptCharged);
                /* update charged multiplicities */
                if (pid % 2 == 0) {
                  partMultPos[kDptDptCharged]++;
                }
                if (pid % 2 == 1) {
                  partMultNeg[kDptDptCharged]++;
                }
              }
              /* fill the species  histograms */
              fillParticleHistosAfterSelection(particle, mccollision, charge, sp);
              /* update species multiplicities */
              if (pid % 2 == 0) {
                partMultPos[sp]++;
              }
              if (pid % 2 == 1) {
                partMultNeg[sp]++;
              }
              acceptedparticles++;
            } else {
              pid = -1;
            }
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

  void filterRecoWithoutPID(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const& collisions, DptDptFullTracks const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterRecoWithoutPID, "Track filtering without PID information", true)

  void filterGenerated(soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo> const& gencollisions, aod::McParticles const& particles)
  {
    filterParticles(gencollisions, particles);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterGenerated, "Generated particles filering", true)
};

template <typename ParticleObject>
inline MatchRecoGenSpecies DptDptFilterTracks::IdentifyParticle(ParticleObject const& particle)
{
  using namespace dptdptfilter;

  constexpr int pdgcodeEl = 11;
  constexpr int pdgcodeMu = 13;
  constexpr int pdgcodePi = 211;
  constexpr int pdgcodeKa = 321;
  constexpr int pdgcodePr = 2212;

  int pdgcode = abs(particle.pdgCode());

  switch (pdgcode) {
    case pdgcodeEl:
      return kDptDptElectron;
      break;
    case pdgcodeMu:
      return kDptDptMuon;
      break;
    case pdgcodePi:
      return kDptDptPion;
      break;
    case pdgcodeKa:
      return kDptDptKaon;
      break;
    case pdgcodePr:
      return kDptDptProton;
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
inline MatchRecoGenSpecies DptDptFilterTracks::IdentifyTrack(TrackObject const& track)
{
  using namespace o2::analysis::dptdptfilter;

  float nsigmas[kDptDptNoOfSpecies];
  if (track.p() < 0.8) {
    nsigmas[kDptDptCharged] = 999.0f;
    nsigmas[kDptDptElectron] = track.tpcNSigmaEl();
    nsigmas[kDptDptMuon] = track.tpcNSigmaMu();
    nsigmas[kDptDptPion] = track.tpcNSigmaPi();
    nsigmas[kDptDptKaon] = track.tpcNSigmaKa();
    nsigmas[kDptDptProton] = track.tpcNSigmaPr();
  } else {
    /* introduce require TOF flag */
    if (track.hasTOF()) {
      nsigmas[kDptDptCharged] = 999.0f;
      nsigmas[kDptDptElectron] = sqrtf(track.tpcNSigmaEl() * track.tpcNSigmaEl() + track.tofNSigmaEl() * track.tofNSigmaEl());
      nsigmas[kDptDptMuon] = sqrtf(track.tpcNSigmaMu() * track.tpcNSigmaMu() + track.tofNSigmaMu() * track.tofNSigmaMu());
      nsigmas[kDptDptPion] = sqrtf(track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi());
      nsigmas[kDptDptKaon] = sqrtf(track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa());
      nsigmas[kDptDptProton] = sqrtf(track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr());
    } else {
      nsigmas[kDptDptCharged] = 999.0f;
      nsigmas[kDptDptElectron] = track.tpcNSigmaEl();
      nsigmas[kDptDptMuon] = track.tpcNSigmaMu();
      nsigmas[kDptDptPion] = track.tpcNSigmaPi();
      nsigmas[kDptDptKaon] = track.tpcNSigmaKa();
      nsigmas[kDptDptProton] = track.tpcNSigmaPr();
    }
  }
  float min_nsigma = 999.0f;
  MatchRecoGenSpecies sp_min_nsigma = kWrongSpecies;
  for (int sp = 0; sp < kDptDptNoOfSpecies; ++sp) {
    if (nsigmas[sp] < min_nsigma) {
      min_nsigma = nsigmas[sp];
      sp_min_nsigma = MatchRecoGenSpecies(sp);
    }
  }
  bool doublematch = false;
  if (min_nsigma < 3.0) {
    for (int sp = 0; (sp < kDptDptNoOfSpecies) && !doublematch; ++sp) {
      if (sp != sp_min_nsigma) {
        if (nsigmas[sp] < 3.0) {
          doublematch = true;
        }
      }
    }
    if (doublematch) {
      return kWrongSpecies;
    } else {
      return sp_min_nsigma;
    }
  } else {
    return kWrongSpecies;
  }
}

template <typename TrackObject>
MatchRecoGenSpecies DptDptFilterTracks::trackIdentification(TrackObject const& track)
{
  using namespace dptdptfilter;

  MatchRecoGenSpecies sp = kWrongSpecies;
  if (recoIdMethod == 0) {
    sp = kDptDptCharged;
  } else if (recoIdMethod == 1) {
    if constexpr (framework::has_type_v<aod::pidtpc_tiny::TPCNSigmaStorePi, typename TrackObject::all_columns>) {
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
  return sp;
}

template <typename TrackObject>
int8_t DptDptFilterTracks::selectTrack(TrackObject const& track)
{
  using namespace dptdptfilter;

  /* before track selection */
  fillTrackHistosBeforeSelection(track);

  /* track selection */
  int8_t pid = AcceptTrack(track);
  if (!(pid < 0)) {
    /* the track has been accepted */
    /* let's identify it */
    /* TODO: probably this needs to go inside AcceptTrack */
    MatchRecoGenSpecies sp = trackIdentification(track);
    if (sp != kWrongSpecies) {
      if (sp != kDptDptCharged) {
        /* fill the charged histograms */
        fillTrackHistosAfterSelection(track, kDptDptCharged);
        /* update charged multiplicities */
        if (pid % 2 == 0) {
          trkMultPos[kDptDptCharged]++;
        }
        if (pid % 2 == 1) {
          trkMultNeg[kDptDptCharged]++;
        }
      }
      /* fill the species histograms */
      fillTrackHistosAfterSelection(track, sp);
      /* update species multiplicities */
      if (pid % 2 == 0) {
        trkMultPos[sp]++;
      }
      if (pid % 2 == 1) {
        trkMultNeg[sp]++;
      }
    }
  }
  return pid;
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
void DptDptFilterTracks::fillTrackHistosAfterSelection(TrackObject const& track, MatchRecoGenSpecies sp)
{
  /* the charged species should have been called first so avoid double counting */
  if (sp == kDptDptCharged) {
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
  }
  fhPA[sp]->Fill(track.p());
  fhPtA[sp]->Fill(track.pt());
  if (track.sign() > 0) {
    fhPtPosA[sp]->Fill(track.pt());
  } else {
    fhPtNegA[sp]->Fill(track.pt());
  }
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
void DptDptFilterTracks::fillParticleHistosAfterSelection(ParticleObject const& particle, MCCollisionObject const& collision, float charge, MatchRecoGenSpecies sp)
{
  /* the charged species should have been called first so avoid double counting */
  if (sp == kDptDptCharged) {
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
  WorkflowSpec workflow{adaptAnalysisTask<DptDptFilter>(cfgc,
                                                        SetDefaultProcesses{
                                                          {{"processWithoutCent", true},
                                                           {"processWithoutCentMC", true}}}),
                        adaptAnalysisTask<DptDptFilterTracks>(cfgc)};
  return workflow;
}
