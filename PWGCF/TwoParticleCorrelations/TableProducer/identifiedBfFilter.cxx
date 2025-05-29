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

/// \file identifiedBfFilter.cxx
/// \brief Filters collisions and tracks according to selection criteria
/// \author bghanley1995@gmail.com

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
#include "Framework/O2DatabasePDGPlugin.h"
#include <TROOT.h>
#include <TParameter.h>
#include <TPDGCode.h>
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
using IdBfTracksTOF = soa::Join<aod::TOFSignal, aod::pidTOFFlags, aod::pidTOFbeta, aod::pidTOFmass>;
using IdBfTracksFullPID = soa::Join<aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
using IdBfFullTracksPID = soa::Join<IdBfFullTracks, IdBfTracksPID, IdBfTracksTOF>;
using IdBfFullTracksFullPID = soa::Join<IdBfFullTracks, IdBfTracksFullPID>;
using IdBfFullTracksPIDAmbiguous = soa::Join<IdBfFullTracksAmbiguous, IdBfTracksPID, IdBfTracksTOF>;
using IdBfFullTracksFullPIDAmbiguous = soa::Join<IdBfFullTracksAmbiguous, IdBfTracksFullPID>;
using IdBfFullTracksDetLevel = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TracksDCA>;
using IdBfFullTracksDetLevelAmbiguous = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TracksDCA, aod::TrackCompColls>;
using IdBfFullTracksPIDDetLevel = soa::Join<IdBfFullTracksDetLevel, IdBfTracksPID, IdBfTracksTOF>;
using IdBfFullTracksFullPIDDetLevel = soa::Join<IdBfFullTracksDetLevel, IdBfTracksFullPID>;
using IdBfFullTracksPIDDetLevelAmbiguous = soa::Join<IdBfFullTracksDetLevelAmbiguous, IdBfTracksPID, IdBfTracksTOF>;
using IdBfFullTracksFullPIDDetLevelAmbiguous = soa::Join<IdBfFullTracksDetLevelAmbiguous, IdBfTracksFullPID>;

bool fullDerivedData = false; /* produce full derived data for its external storage */

TList* ccdblst = nullptr;
bool loadfromccdb = false;

std::vector<int> recoIdMethods = {0, 1, 2}; // Reconstructed PID Methods, 0 is no PID, 1 is calculated PID, 2 is MC PID
std::vector<int> trackTypes = {0, 1, 2, 3};
const int twoDenom = 2; // Used to test if a value is even or odd

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

TH2F* fhPtEtaPosA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPtEtaNegA[kIdBfNoOfSpecies + 1] = {nullptr};

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
TH2F* fhNSigmaTPCIdTrks[kIdBfNoOfSpecies] = {nullptr};

TH1F* fhNSigmaCorrection[kIdBfNoOfSpecies] = {nullptr};

TH1F* fhEtaB = nullptr;
TH1F* fhEtaA = nullptr;

TH1F* fhPhiB = nullptr;
TH1F* fhPhiA = nullptr;

TH1F* fhTrackLengthB = nullptr;
TH1F* fhTrackLengthTOFB = nullptr;
TH2F* fhTrackTimeB = nullptr;
TH2F* fhTrackBetaInvB = nullptr;
TH2F* fhTrackTimeIPB = nullptr;
TH2F* fhTrackBetaInvIPB = nullptr;
TH2F* fhdEdxB = nullptr;
TH2F* fhdEdxIPTPCB = nullptr;
TH1F* fhTrackLengthA[kIdBfNoOfSpecies + 2] = {nullptr};
TH1F* fhTrackLengthTOFA[kIdBfNoOfSpecies + 2] = {nullptr};
TH2F* fhdEdxA[kIdBfNoOfSpecies + 2] = {nullptr};
TH2F* fhdEdxIPTPCA[kIdBfNoOfSpecies + 2] = {nullptr};
TH2F* fhTrackTimeA[kIdBfNoOfSpecies + 2] = {nullptr};
TH2F* fhTrackBetaInvA[kIdBfNoOfSpecies + 2] = {nullptr};
TH2F* fhTrackTimeIPA[kIdBfNoOfSpecies + 2] = {nullptr};
TH2F* fhTrackBetaInvIPA[kIdBfNoOfSpecies + 2] = {nullptr};

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

TH2S* fhTruePIDMismatch = nullptr;
TH1S* fhTruePIDCorrect = nullptr;

std::vector<std::vector<TH2F*>> fhTrueNSigmaTPC = {o2::analysis::identifiedbffilter::kIdBfNoOfSpecies, {o2::analysis::identifiedbffilter::kIdBfNoOfSpecies, nullptr}};
std::vector<std::vector<TH2F*>> fhTrueNSigmaTOF = {o2::analysis::identifiedbffilter::kIdBfNoOfSpecies, {o2::analysis::identifiedbffilter::kIdBfNoOfSpecies, nullptr}};

std::vector<std::vector<TH2F*>> fhPrimaryNSigmaTPC = {o2::analysis::identifiedbffilter::kIdBfNoOfSpecies, {o2::analysis::identifiedbffilter::kIdBfNoOfSpecies, nullptr}};
std::vector<std::vector<TH2F*>> fhPrimaryNSigmaTOF = {o2::analysis::identifiedbffilter::kIdBfNoOfSpecies, {o2::analysis::identifiedbffilter::kIdBfNoOfSpecies, nullptr}};

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

TH2F* fhTruedEdxA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhTruedEdxIPTPCA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhTrueTrackLengthA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhTrueTrackLengthTOFA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhTrueTrackTimeA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhTrueTrackBetaInvA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhTrueTrackTimeIPA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhTrueTrackBetaInvIPA[kIdBfNoOfSpecies + 1] = {nullptr};

TH1F* fhPrimaryPA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhPrimaryPtA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPrimarydEdxA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPrimarydEdxIPTPCA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhPrimaryTrackLengthA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhPrimaryTrackLengthTOFA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPrimaryTrackTimeA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPrimaryTrackBetaInvA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPrimaryTrackTimeIPA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPrimaryTrackBetaInvIPA[kIdBfNoOfSpecies + 1] = {nullptr};

TH1F* fhPurePA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhPurePtA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPuredEdxA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPuredEdxIPTPCA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhPureTrackLengthA[kIdBfNoOfSpecies + 1] = {nullptr};
TH1F* fhPureTrackLengthTOFA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPureTrackTimeA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPureTrackBetaInvA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPureTrackTimeIPA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhPureTrackBetaInvIPA[kIdBfNoOfSpecies + 1] = {nullptr};

TH2F* fhTruePtEtaPosA[kIdBfNoOfSpecies + 1] = {nullptr};
TH2F* fhTruePtEtaNegA[kIdBfNoOfSpecies + 1] = {nullptr};

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

TH1F* fhPrimaryPB = nullptr;
TH1F* fhPrimaryPtB = nullptr;
TH2F* fhPrimarydEdxB = nullptr;
TH2F* fhPrimarydEdxIPTPCB = nullptr;
TH1F* fhPrimaryTrackLengthB = nullptr;
TH1F* fhPrimaryTrackLengthTOFB = nullptr;
TH2F* fhPrimaryTrackTimeB = nullptr;
TH2F* fhPrimaryTrackBetaInvB = nullptr;
TH2F* fhPrimaryTrackTimeIPB = nullptr;
TH2F* fhPrimaryTrackBetaInvIPB = nullptr;

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
  Configurable<bool> cfgFullDerivedData{"cfgFullDerivedData", false, "Produce the full derived data for external storage. Default false"};
  Configurable<std::string> cfgCentMultEstimator{"cfgCentMultEstimator", "V0M", "Centrality/multiplicity estimator detector: V0M,CL0,CL1,FV0A,FT0M,FT0A,FT0C,NTPV,NOCM: none. Default V0M"};
  Configurable<std::string> cfgSystem{"cfgSystem", "PbPb", "System: pp, PbPb, Pbp, pPb, XeXe, ppRun3, PbPbRun3. Default PbPb"};
  Configurable<std::string> cfgDataType{"cfgDataType", "data", "Data type: data, datanoevsel, MC, FastMC, OnTheFlyMC. Default data"};
  Configurable<std::string> cfgTriggSel{"cfgTriggSel", "MB", "Trigger selection: MB, None. Default MB"};
  Configurable<o2::analysis::DptDptBinningCuts> cfgBinning{"cfgBinning",
                                                           {28, -7.0, 7.0, 18, 0.2, 2.0, 16, -0.8, 0.8, 72, 0.5},
                                                           "triplets - nbins, min, max - for z_vtx, pT, eta and phi, binning plus bin fraction of phi origin shift"};
  Configurable<bool> cfgTraceCollId0{"cfgTraceCollId0", false, "Trace particles in collisions id 0. Default false"};

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
  if (isEvtSelected(collision, centormult)) {
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
  if (isEvtSelected(mccollision, centormult)) {
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
  for (const auto& tmpcollision : collisions) {
    if (tmpcollision.has_mcCollision()) {
      if (tmpcollision.mcCollisionId() == mccollision.globalIndex()) {
        typename AllCollisions::iterator const& collision = allcollisions.iteratorAt(tmpcollision.globalIndex());
        if (isEvtSelected(collision, defaultcent)) {
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
    if (isEvtSelected(mccollision, centmult)) {
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
  T sqSum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  T stdev = std::sqrt(sqSum / vec.size());

  return stdev;
}

struct IdentifiedBfFilterTracks {

  struct : ConfigurableGroup {
    Configurable<std::string> inputCCDBUrl{"inputCCDBUrl", "http://ccdb-test.cern.ch:8080", "The CCDB url for the input file"};
    Configurable<std::string> inputCCDBPathName{"inputCCDBPathName", "", "The CCDB path for the input file. Default \"\", i.e. don't load from CCDB"};
    Configurable<std::string> inputCCDBDate{"inputCCDBDate", "20220307", "The CCDB date for the input file"};
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

  Service<o2::framework::O2DatabasePDG> fPDG;
  Produces<aod::ScannedTracks> scannedtracks;
  Produces<aod::IdentifiedBfCFTracksInfo> tracksinfo;
  Produces<aod::ScannedTrueTracks> scannedgentracks;
  Produces<aod::IdentifiedBfCFGenTracksInfo> gentracksinfo;

  Configurable<bool> cfgFullDerivedData{"cfgFullDerivedData", false, "Produce the full derived data for external storage. Default false"};
  Configurable<int> cfgTrackType{"cfgTrackType", 1, "Type of selected tracks: 0 = no selection, 1 = Run2 global tracks FB96, 3 = Run3 tracks, 5 = Run2 TPC only tracks, 7 = Run 3 TPC only tracks. Default 1"};
  Configurable<std::string> cfgSystem{"cfgSystem", "PbPb", "System: pp, PbPb, Pbp, pPb, XeXe, ppRun3, PbPbRun3. Default PbPb"};
  Configurable<std::string> cfgDataType{"cfgDataType", "data", "Data type: data, datanoevsel, MC, FastMC, OnTheFlyMC. Default data"};
  Configurable<o2::analysis::DptDptBinningCuts> cfgBinning{"cfgBinning",
                                                           {28, -7.0, 7.0, 18, 0.2, 2.0, 16, -0.8, 0.8, 72, 0.5},
                                                           "triplets - nbins, min, max - for z_vtx, pT, eta and phi, binning plus bin fraction of phi origin shift"};
  Configurable<o2::analysis::CheckRangeCfg> cfgTraceDCAOutliers{"cfgTraceDCAOutliers", {false, 0.0, 0.0}, "Track the generator level DCAxy outliers: false/true, low dcaxy, up dcaxy. Default {false,0.0,0.0}"};
  Configurable<float> cfgTraceOutOfSpeciesParticles{"cfgTraceOutOfSpeciesParticles", false, "Track the particles which are not e,mu,pi,K,p: false/true. Default false"};
  Configurable<int> cfgRecoIdMethod{"cfgRecoIdMethod", 0, "Method for identifying reconstructed tracks: 0 No PID, 1 PID, 2 mcparticle. Default 0"};
  Configurable<o2::analysis::TrackSelectionCfg> cfgTrackSelection{"cfgTrackSelection", {false, false, 0, 70, 0.8, 2.4, 3.2}, "Track selection: {useit: true/false, ongen: true/false, tpccls, tpcxrws, tpcxrfc, dcaxy, dcaz}. Default {false,0.70.0.8,2.4,3.2}"};
  Configurable<bool> reqTOF{"reqTOF", false, "Require TOF data for PID. Default false"};
  Configurable<bool> onlyTOF{"onlyTOF", false, "Only use TOF data for PID. Default false"};

  // Configurable<int> pidEl{"pidEl", -1, "Identify Electron Tracks"};
  // Configurable<int> pidPi{"pidPi", -1, "Identify Pion Tracks"};
  // Configurable<int> pidKa{"pidKa", -1, "Identify Kaon Tracks"};
  // Configurable<int> pidPr{"pidPr", -1, "Identify Proton Tracks"};

  // Configurable<float> minPIDSigma{"minPIDSigma", -3.0, "Minimum required sigma for PID Acceptance"};
  // Configurable<float> maxPIDSigma{"maxPIDSigma", 3.0, "Maximum required sigma for PID Acceptance"};

  // Configurable<std::vector<float>> minPIDSigmas{"minPIDSigmas", {-3.,-3.,-3.,-3.},"Minimum required sigma for PID Acceptance, {e, pi, K, p}"};
  // Configurable<std::vector<std::vector<float>>> acceptPIDSigmas{"acceptPIDSigmas", {{-3.,3.},{-3.,3.},{-3.,3.},{-3.,3.}},"Sigma range for PID Acceptance, {e, pi, K, p}"};
  // Configurable<std::vector<float>> minRejectSigmas{"minRejectSigmas", {-1.,-1.,-1.,-1.},"Minimum required sigma for PID double match rejection, {e, pi, K, p}"};

  // Configurable<float> minRejectSigma{"minRejectSigma", -1.0, "Minimum required sigma for PID double match rejection"};
  // Configurable<float> maxRejectSigma{"maxRejectSigma", 1.0, "Maximum required sigma for PID double match rejection"};

  Configurable<std::vector<int>> cfgDoPID{"cfgDoPID", {-1, -1, -1, -1}, "Do PID for particle, {e, pi, K, p}"};
  Configurable<std::vector<float>> cfgTOFCut{"cfgTOFCut", {0.8, 0.8, 0.8, 0.8}, "Momentum under which we don't use TOF PID data, {e, pi, K, p}"};
  Configurable<std::vector<float>> cfgTPCCut{"cfgTPCCut", {1.2, 1.2, 1.2, 1.2}, "Momentum over which we don't use TPC PID data, {e, pi, K, p}"};

  Configurable<bool> makeNSigmaPlots{"makeNSigmaPlots", false, "Produce the N Sigma Plots for external storage. Default false"};

  struct : ConfigurableGroup {
    Configurable<std::vector<float>> rejectPIDSigmasEl{"rejectPIDSigmasEl", {0., 1., 1., 1.}, "Required sigma for PID double match rejection of electrons for {e, pi, K, p}"};
    Configurable<std::vector<float>> rejectPIDSigmasPi{"rejectPIDSigmasPi", {1., 0., 1., 1.}, "Required sigma for PID double match rejection of pions for {e, pi, K, p}"};
    Configurable<std::vector<float>> rejectPIDSigmasKa{"rejectPIDSigmasKa", {1., 1., 0., 1.}, "Required sigma for PID double match rejection of kaons for {e, pi, K, p}"};
    Configurable<std::vector<float>> rejectPIDSigmasPr{"rejectPIDSigmasPr", {1., 1., 1., 0.}, "Required sigma for PID double match rejection of protons for {e, pi, K, p}"};
  } rejectPIDSigmas;

  struct : ConfigurableGroup {
    Configurable<std::vector<float>> acceptPIDSigmasEl{"acceptPIDSigmasEl", {-3., 3.}, "Sigma range for PID Acceptance for electrons"};
    Configurable<std::vector<float>> acceptPIDSigmasPi{"acceptPIDSigmasPi", {-3., 3.}, "Sigma range for PID Acceptance for pions"};
    Configurable<std::vector<float>> acceptPIDSigmasKa{"acceptPIDSigmasKa", {-3., 3.}, "Sigma range for PID Acceptance for kaons"};
    Configurable<std::vector<float>> acceptPIDSigmasPr{"acceptPIDSigmasPr", {-3., 3.}, "Sigma range for PID Acceptance for protons"};
  } acceptPIDSigmas;

  OutputObj<TList> fOutput{"IdentifiedBfFilterTracksInfo", OutputObjHandlingPolicy::AnalysisObject};
  bool checkAmbiguousTracks = false;

  void init(InitContext&)
  {
    LOGF(info, "IdentifiedBfFilterTracks::init()");

    // ccdb info
    ccdb->setURL(cfgcentersinputfile.inputCCDBUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    LOGF(info, "Initizalized CCDB");

    loadfromccdb = cfgcentersinputfile.inputCCDBPathName->length() > 0;

    if (ccdblst == nullptr) {
      if (loadfromccdb) {
        LOGF(info, "Loading CCDB Objects");

        ccdblst = getCCDBInput(cfgcentersinputfile.inputCCDBPathName->c_str(), cfgcentersinputfile.inputCCDBDate->c_str());
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

    LOGF(info, "Initializing ranges");
    acceptRange.push_back(acceptPIDSigmas.acceptPIDSigmasEl);
    acceptRange.push_back(acceptPIDSigmas.acceptPIDSigmasPi);
    acceptRange.push_back(acceptPIDSigmas.acceptPIDSigmasKa);
    acceptRange.push_back(acceptPIDSigmas.acceptPIDSigmasPr);

    rejectRange.push_back(rejectPIDSigmas.rejectPIDSigmasEl);
    rejectRange.push_back(rejectPIDSigmas.rejectPIDSigmasPi);
    rejectRange.push_back(rejectPIDSigmas.rejectPIDSigmasKa);
    rejectRange.push_back(rejectPIDSigmas.rejectPIDSigmasPr);

    doPID = cfgDoPID;
    tofCut = cfgTOFCut;
    tpcCut = cfgTPCCut;
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

    /* required ambiguous tracks checks? */
    if (dofilterDetectorLevelWithoutPIDAmbiguous || dofilterDetectorLevelWithPIDAmbiguous || dofilterRecoWithoutPIDAmbiguous || dofilterRecoWithPIDAmbiguous) {
      checkAmbiguousTracks = true;
    }

    /* create the output list which will own the task histograms */
    TList* fOutputList = new TList();
    fOutputList->SetOwner(true);
    fOutput.setObject(fOutputList);

    /* incorporate configuration parameters to the output */
    fOutputList->Add(new TParameter<int>("TrackType", cfgTrackType, 'f'));
    fOutputList->Add(new TParameter<int>("TrackOneCharge", 1, 'f'));
    fOutputList->Add(new TParameter<int>("TrackTwoCharge", -1, 'f'));

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
      fhTrackLengthB = new TH1F(TString::Format("fhTrackLengthB").Data(),
                                TString::Format("Track Length; L (cm)").Data(),
                                1000, 0., 1000.0);
      fhTrackLengthTOFB = new TH1F(TString::Format("fhTrackLengthTOFB").Data(),
                                   TString::Format("Track Length with TOF; L (cm)").Data(),
                                   1000, 0.0, 1000.0);
      fhTrackTimeB = new TH2F(TString::Format("fhTrackTimeB").Data(),
                              TString::Format("Track Time vs P; P (GeV/c); Track Time(ns)").Data(),
                              ptbins, ptlow, ptup, 1000, 0.0, 10.0);
      fhTrackBetaInvB = new TH2F(TString::Format("fhTrackBetaInvB").Data(),
                                 TString::Format("1/#Beta vs P; P (GeV/c); 1/#Beta(ns/m)").Data(),
                                 ptbins, ptlow, ptup, 1000, 0.0, 10.0);
      fhTrackTimeIPB = new TH2F(TString::Format("fhTrackTimeIPB").Data(),
                                TString::Format("Track Time vs P_{IP}; P (GeV/c); Track Time(ns)").Data(),
                                ptbins, ptlow, ptup, 1000, 0.0, 10.0);
      fhTrackBetaInvIPB = new TH2F(TString::Format("fhTrackBetaInvIPB").Data(),
                                   TString::Format("1/#Beta vs P_{IP}; P (GeV/c); 1/#Beta(ns/m)").Data(),
                                   ptbins, ptlow, ptup, 1000, 0.0, 10.0);
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
        fhNSigmaTPCIdTrks[sp] = new TH2F(TString::Format("fhNSigmaTPC_IdTrks_%s", speciesName[sp]).Data(),
                                         TString::Format("N Sigma from TPC vs P for Identified %s;N #sigma;p (GeV/c)", speciesTitle[sp]).Data(),
                                         48, -6, 6,
                                         ptbins, ptlow, ptup);
      }

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
        fhPtEtaPosA[sp] = new TH2F(TString::Format("fHistPtEtaPosA_%s", speciesName[sp]),
                                   TString::Format("P_{T} vs #eta distribution for reconstructed  %s^{#plus};P_{T} (GeV/c);#eta;dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                   ptbins, ptlow, ptup,
                                   etabins, etalow, etaup);
        fhPtEtaNegA[sp] = new TH2F(TString::Format("fHistPtEtaNegA_%s", speciesName[sp]),
                                   TString::Format("P_{T} vs #eta distribution for reconstructed  %s^{#minus};P_{T} (GeV/c);#eta;dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                   ptbins, ptlow, ptup,
                                   etabins, etalow, etaup);
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
        fhTrackLengthA[sp] = new TH1F(TString::Format("fhTrackLengthA_%s", speciesName[sp]).Data(),
                                      TString::Format("Track Length of reconstructed %s; L (cm)", speciesTitle[sp]).Data(),
                                      1000, 0.0, 1000.0);
        fhTrackLengthTOFA[sp] = new TH1F(TString::Format("fhTrackLengthTOFA_%s", speciesName[sp]).Data(),
                                         TString::Format("Track Length of reconstructed %s with TOF; L (cm)", speciesTitle[sp]).Data(),
                                         1000, 0.0, 1000.0);
        fhTrackTimeA[sp] = new TH2F(TString::Format("fhTrackTimeA_%s", speciesName[sp]).Data(),
                                    TString::Format("Track Time vs P reconstructed %s; P (GeV/c); Track Time(ns)", speciesTitle[sp]).Data(), ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        fhTrackBetaInvA[sp] = new TH2F(TString::Format("fhTrackBetaInvA_%s", speciesName[sp]).Data(),
                                       TString::Format("1/#Beta vs P reconstructed %s; P (GeV/c); 1/#Beta(ns/m)", speciesTitle[sp]).Data(),
                                       ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        fhTrackTimeIPA[sp] = new TH2F(TString::Format("fhTrackTimeIPA_%s", speciesName[sp]).Data(),
                                      TString::Format("Track Time vs P_{IP} reconstructed %s; P (GeV/c); Track Time(ns)", speciesTitle[sp]).Data(),
                                      ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        fhTrackBetaInvIPA[sp] = new TH2F(TString::Format("fhTrackBetaInvIPA_%s", speciesName[sp]).Data(),
                                         TString::Format("1/#Beta vs P_{IP} reconstructed %s; P (GeV/c); 1/#Beta(ns/m)", speciesTitle[sp]).Data(),
                                         ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        LOGF(info, "Made Histos");
      }
      fhdEdxA[kIdBfNoOfSpecies + 1] = new TH2F(TString::Format("fhdEdxA_WrongSpecies").Data(),
                                               TString::Format("dE/dx vs P reconstructed Wrong Species; P (GeV/c); dE/dx (a.u.)").Data(),
                                               ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
      fhdEdxIPTPCA[kIdBfNoOfSpecies + 1] = new TH2F(TString::Format("fhdEdxIPTPCA_WrongSpecies").Data(),
                                                    TString::Format("dE/dx vs P_{IP} reconstructed Wrong Species; P (GeV/c); dE/dx (a.u.)").Data(),
                                                    ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
      fhTrackLengthA[kIdBfNoOfSpecies + 1] = new TH1F(TString::Format("fhTrackLengthA_WrongSpecies").Data(),
                                                      TString::Format("Track Length of reconstructed Wrong Species; L (cm)").Data(),
                                                      1000, 0.0, 1000.0);
      fhTrackLengthTOFA[kIdBfNoOfSpecies + 1] = new TH1F(TString::Format("fhTrackLengthTOFA_WrongSpecies").Data(),
                                                         TString::Format("Track Length of reconstructed Wrong Species with TOF; L (cm)").Data(),
                                                         1000, 0.0, 1000.0);
      fhTrackTimeA[kIdBfNoOfSpecies + 1] = new TH2F(TString::Format("fhTrackTimeA_WrongSpecies").Data(),
                                                    TString::Format("Track Time vs P reconstructed Wrong Species; P (GeV/c); Track Time(ns)").Data(),
                                                    ptbins, ptlow, ptup, 1000, 0.0, 10.0);
      fhTrackBetaInvA[kIdBfNoOfSpecies + 1] = new TH2F(TString::Format("fhTrackBetaInvA_WrongSpecies").Data(),
                                                       TString::Format("1/#Beta vs P reconstructed Wrong Species; P (GeV/c); 1/#Beta(ns/m)").Data(),
                                                       ptbins, ptlow, ptup, 1000, 0.0, 10.0);
      fhTrackTimeIPA[kIdBfNoOfSpecies + 1] = new TH2F(TString::Format("fhTrackTimeIPA_WrongSpecies").Data(),
                                                      TString::Format("Track Time vs P_{IP} reconstructed Wrong Species; P (GeV/c); Track Time(ns)").Data(),
                                                      ptbins, ptlow, ptup, 1000, 0.0, 10.0);
      fhTrackBetaInvIPA[kIdBfNoOfSpecies + 1] = new TH2F(TString::Format("fhTrackBetaInvIPA_WrongSpecies").Data(),
                                                         TString::Format("1/#Beta vs P_{IP} reconstructed Wrong Species; P (GeV/c); 1/#Beta(ns/m)").Data(),
                                                         ptbins, ptlow, ptup, 1000, 0.0, 10.0);

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
      fOutputList->Add(fhTrackLengthB);
      fOutputList->Add(fhTrackLengthTOFB);
      fOutputList->Add(fhTrackTimeB);
      fOutputList->Add(fhTrackBetaInvB);
      fOutputList->Add(fhTrackTimeIPB);
      fOutputList->Add(fhTrackBetaInvIPB);
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
        fOutputList->Add(fhNSigmaTPCIdTrks[sp]);
      }

      for (int sp = 0; sp < kIdBfNoOfSpecies + 1; ++sp) {
        fOutputList->Add(fhPA[sp]);
        fOutputList->Add(fhPtA[sp]);
        fOutputList->Add(fhPtPosA[sp]);
        fOutputList->Add(fhPtNegA[sp]);
        fOutputList->Add(fhPtEtaPosA[sp]);
        fOutputList->Add(fhPtEtaNegA[sp]);
        fOutputList->Add(fhNPosNegA[sp]);
        fOutputList->Add(fhDeltaNA[sp]);
        fOutputList->Add(fhdEdxA[sp]);
        fOutputList->Add(fhdEdxIPTPCA[sp]);
        fOutputList->Add(fhTrackLengthA[sp]);
        fOutputList->Add(fhTrackLengthTOFA[sp]);
        fOutputList->Add(fhTrackTimeA[sp]);
        fOutputList->Add(fhTrackBetaInvA[sp]);
        fOutputList->Add(fhTrackTimeIPA[sp]);
        fOutputList->Add(fhTrackBetaInvIPA[sp]);
        LOGF(info, "Added histos");
      }
      fOutputList->Add(fhdEdxA[kIdBfNoOfSpecies + 1]);
      fOutputList->Add(fhdEdxIPTPCA[kIdBfNoOfSpecies + 1]);
      fOutputList->Add(fhTrackLengthA[kIdBfNoOfSpecies + 1]);
      fOutputList->Add(fhTrackLengthTOFA[kIdBfNoOfSpecies + 1]);
      fOutputList->Add(fhTrackTimeA[kIdBfNoOfSpecies + 1]);
      fOutputList->Add(fhTrackBetaInvA[kIdBfNoOfSpecies + 1]);
      fOutputList->Add(fhTrackTimeIPA[kIdBfNoOfSpecies + 1]);
      fOutputList->Add(fhTrackBetaInvIPA[kIdBfNoOfSpecies + 1]);
    }

    if ((fDataType != kData) && (fDataType != kDataNoEvtSel)) {
      /* create the true data histograms */

      fhTruePIDMismatch = new TH2S("fHistTruePIDMismatch", "Mismatched Generated and Reconstructed PID;Generated Species;Reconstructed Species", kIdBfNoOfSpecies, 0, kIdBfNoOfSpecies, kIdBfNoOfSpecies, 0, kIdBfNoOfSpecies);
      fhTruePIDCorrect = new TH1S("fHistTruePIDCorrect", "Correct PID between Generated and Reconstructed PID;Species", kIdBfNoOfSpecies, 0, kIdBfNoOfSpecies);

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
      fhPrimaryPB = new TH1F(TString::Format("fhPrimaryPB").Data(),
                             TString::Format("p distribution Primary Before Selection;p (GeV/c);dN/dp (c/GeV)").Data(),
                             ptbins, ptlow, ptup);
      fhPrimaryPtB = new TH1F(TString::Format("fhPrimaryPtB"),
                              TString::Format("p_{T} distribution Primary Before Selection ;p_{T} (GeV/c);dN/dP_{T} (c/GeV)").Data(),
                              ptbins, ptlow, ptup);
      fhPrimarydEdxB = new TH2F(TString::Format("fhPrimarydEdxB").Data(),
                                TString::Format("dE/dx vs P Primary Before Selection; P (GeV/c); dE/dx (a.u.)").Data(),
                                ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
      fhPrimarydEdxIPTPCB = new TH2F(TString::Format("fhPrimarydEdxIPTPCB").Data(),
                                     TString::Format("dE/dx vs P_{IP} Primary Before Selection; P (GeV/c); dE/dx (a.u.)").Data(),
                                     ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
      fhPrimaryTrackLengthB = new TH1F(TString::Format("fhPrimaryTrackLengthB").Data(),
                                       TString::Format("Track Length of Primary Before Selection; L (cm)").Data(),
                                       1000, 0.0, 1000.0);
      fhPrimaryTrackLengthTOFB = new TH1F(TString::Format("fhPrimaryTrackLengthTOFB").Data(),
                                          TString::Format("Track Length of Primary Before Selection with TOF; L (cm)").Data(),
                                          1000, 0.0, 1000.0);
      fhPrimaryTrackTimeB = new TH2F(TString::Format("fhPrimaryTrackTimeB").Data(),
                                     TString::Format("Track Time vs P Primary Before Selection; P (GeV/c); Track Time(ns)").Data(),
                                     ptbins, ptlow, ptup, 1000, 0.0, 10.0);
      fhPrimaryTrackBetaInvB = new TH2F(TString::Format("fhPrimaryTrackBetaInvB").Data(),
                                        TString::Format("1/#Beta vs P Primary Before Selection; P (GeV/c); 1/#Beta(ns/m)").Data(),
                                        ptbins, ptlow, ptup, 1000, 0.0, 10.0);
      fhPrimaryTrackTimeIPB = new TH2F(TString::Format("fhPrimaryTrackTimeIPB").Data(),
                                       TString::Format("Track Time vs P_{IP} Primary Before Selection; P (GeV/c); Track Time(ns)").Data(),
                                       ptbins, ptlow, ptup, 1000, 0.0, 10.0);
      fhPrimaryTrackBetaInvIPB = new TH2F(TString::Format("fhPrimaryTrackBetaInvIPB").Data(),
                                          TString::Format("1/#Beta vs P_{IP} Primary Before Selection; P (GeV/c); 1/#Beta(ns/m)").Data(),
                                          ptbins, ptlow, ptup, 1000, 0.0, 10.0);
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
        fhTruePtEtaPosA[sp] = new TH2F(TString::Format("fTrueHistPtEtaPosA_%s", speciesName[sp]),
                                       TString::Format("P_{T} vs #eta distribution %s^{#plus} (truth);P_{T} (GeV/c);#eta;dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                       ptbins, ptlow, ptup,
                                       etabins, etalow, etaup);
        fhTruePtEtaNegA[sp] = new TH2F(TString::Format("fTrueHistPtEtaNegA_%s", speciesName[sp]),
                                       TString::Format("P_{T} vs #eta distribution %s^{#minus} (truth);P_{T} (GeV/c);#eta;dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                       ptbins, ptlow, ptup,
                                       etabins, etalow, etaup);
        fhTrueNPosNegA[sp] = new TH2F(TString::Format("fhTrueNPosNegA_%s", speciesName[sp]).Data(),
                                      TString::Format("N(%s^{#plus}) N(%s^{#minus}) distribution (truth);N(%s^{#plus});N(%s^{#minus})", speciesTitle[sp], speciesTitle[sp], speciesTitle[sp], speciesTitle[sp]).Data(),
                                      40, -0.5, 39.5, 40, -0.5, 39.5);
        fhTrueDeltaNA[sp] = new TH1F(TString::Format("fhTrueDeltaNA_%s", speciesName[sp]).Data(),
                                     TString::Format("N(%s^{#plus}) #minus N(%s^{#minus}) distribution (truth);N(%s^{#plus}) #minus N(%s^{#minus})", speciesTitle[sp], speciesTitle[sp], speciesTitle[sp], speciesTitle[sp]).Data(),
                                     79, -39.5, 39.5);
        fhTruedEdxA[sp] = new TH2F(TString::Format("fhTruedEdxA_%s", speciesName[sp]).Data(),
                                   TString::Format("dE/dx vs P generated %s; P (GeV/c); dE/dx (a.u.)", speciesTitle[sp]).Data(),
                                   ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
        fhTruedEdxIPTPCA[sp] = new TH2F(TString::Format("fhTruedEdxIPTPCA_%s", speciesName[sp]).Data(),
                                        TString::Format("dE/dx vs P_{IP} generated %s; P (GeV/c); dE/dx (a.u.)", speciesTitle[sp]).Data(),
                                        ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
        fhTrueTrackLengthA[sp] = new TH1F(TString::Format("fhTrueTrackLengthA_%s", speciesName[sp]).Data(),
                                          TString::Format("Track Length of generated %s; L (cm)", speciesTitle[sp]).Data(),
                                          1000, 0.0, 1000.0);
        fhTrueTrackLengthTOFA[sp] = new TH1F(TString::Format("fhTrueTrackLengthTOFA_%s", speciesName[sp]).Data(),
                                             TString::Format("Track Length of generated %s with TOF; L (cm)", speciesTitle[sp]).Data(),
                                             1000, 0.0, 1000.0);
        fhTrueTrackTimeA[sp] = new TH2F(TString::Format("fhTrueTrackTimeA_%s", speciesName[sp]).Data(),
                                        TString::Format("Track Time vs P generated %s; P (GeV/c); Track Time(ns)", speciesTitle[sp]).Data(),
                                        ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        fhTrueTrackBetaInvA[sp] = new TH2F(TString::Format("fhTrueTrackBetaInvA_%s", speciesName[sp]).Data(),
                                           TString::Format("1/#Beta vs P generated %s; P (GeV/c); 1/#Beta(ns/m)", speciesTitle[sp]).Data(),
                                           ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        fhTrueTrackTimeIPA[sp] = new TH2F(TString::Format("fhTrueTrackTimeIPA_%s", speciesName[sp]).Data(),
                                          TString::Format("Track Time vs P_{IP} generated %s; P (GeV/c); Track Time(ns)", speciesTitle[sp]).Data(),
                                          ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        fhTrueTrackBetaInvIPA[sp] = new TH2F(TString::Format("fhTrueTrackBetaInvIPA_%s", speciesName[sp]).Data(),
                                             TString::Format("1/#Beta vs P_{IP} generated %s; P (GeV/c); 1/#Beta(ns/m)", speciesTitle[sp]).Data(),
                                             ptbins, ptlow, ptup, 1000, 0.0, 10.0);

        fhPrimaryPA[sp] = new TH1F(TString::Format("fhPrimaryPA_%s", speciesName[sp]).Data(),
                                   TString::Format("p distribution Primary %s;p (GeV/c);dN/dp (c/GeV)", speciesTitle[sp]).Data(),
                                   ptbins, ptlow, ptup);
        fhPrimaryPtA[sp] = new TH1F(TString::Format("fhPrimaryPtA_%s", speciesName[sp]),
                                    TString::Format("p_{T} distribution Primary %s ;p_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                    ptbins, ptlow, ptup);
        fhPrimarydEdxA[sp] = new TH2F(TString::Format("fhPrimarydEdxA_%s", speciesName[sp]).Data(),
                                      TString::Format("dE/dx vs P Primary %s; P (GeV/c); dE/dx (a.u.)", speciesTitle[sp]).Data(),
                                      ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
        fhPrimarydEdxIPTPCA[sp] = new TH2F(TString::Format("fhPrimarydEdxIPTPCA_%s", speciesName[sp]).Data(),
                                           TString::Format("dE/dx vs P_{IP} Primary %s; P (GeV/c); dE/dx (a.u.)", speciesTitle[sp]).Data(),
                                           ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
        fhPrimaryTrackLengthA[sp] = new TH1F(TString::Format("fhPrimaryTrackLengthA_%s", speciesName[sp]).Data(),
                                             TString::Format("Track Length of Primary %s; L (cm)", speciesTitle[sp]).Data(),
                                             1000, 0.0, 1000.0);
        fhPrimaryTrackLengthTOFA[sp] = new TH1F(TString::Format("fhPrimaryTrackLengthTOFA_%s", speciesName[sp]).Data(),
                                                TString::Format("Track Length of Primary %s with TOF; L (cm)", speciesTitle[sp]).Data(),
                                                1000, 0.0, 1000.0);
        fhPrimaryTrackTimeA[sp] = new TH2F(TString::Format("fhPrimaryTrackTimeA_%s", speciesName[sp]).Data(),
                                           TString::Format("Track Time vs P Primary %s; P (GeV/c); Track Time(ns)", speciesTitle[sp]).Data(),
                                           ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        fhPrimaryTrackBetaInvA[sp] = new TH2F(TString::Format("fhPrimaryTrackBetaInvA_%s", speciesName[sp]).Data(),
                                              TString::Format("1/#Beta vs P Primary %s; P (GeV/c); 1/#Beta(ns/m)", speciesTitle[sp]).Data(),
                                              ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        fhPrimaryTrackTimeIPA[sp] = new TH2F(TString::Format("fhPrimaryTrackTimeIPA_%s", speciesName[sp]).Data(),
                                             TString::Format("Track Time vs P_{IP} Primary %s; P (GeV/c); Track Time(ns)", speciesTitle[sp]).Data(),
                                             ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        fhPrimaryTrackBetaInvIPA[sp] = new TH2F(TString::Format("fhPrimaryTrackBetaInvIPA_%s", speciesName[sp]).Data(),
                                                TString::Format("1/#Beta vs P_{IP} Primary %s; P (GeV/c); 1/#Beta(ns/m)", speciesTitle[sp]).Data(),
                                                ptbins, ptlow, ptup, 1000, 0.0, 10.0);

        fhPurePA[sp] = new TH1F(TString::Format("fhPurePA_%s", speciesName[sp]).Data(),
                                TString::Format("p distribution Pure %s;p (GeV/c);dN/dp (c/GeV)", speciesTitle[sp]).Data(),
                                ptbins, ptlow, ptup);
        fhPurePtA[sp] = new TH1F(TString::Format("fhPurePtA_%s", speciesName[sp]),
                                 TString::Format("p_{T} distribution Pure %s ;p_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                 ptbins, ptlow, ptup);
        fhPuredEdxA[sp] = new TH2F(TString::Format("fhPuredEdxA_%s", speciesName[sp]).Data(),
                                   TString::Format("dE/dx vs P Pure %s; P (GeV/c); dE/dx (a.u.)", speciesTitle[sp]).Data(),
                                   ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
        fhPuredEdxIPTPCA[sp] = new TH2F(TString::Format("fhPuredEdxIPTPCA_%s", speciesName[sp]).Data(),
                                        TString::Format("dE/dx vs P_{IP} Pure %s; P (GeV/c); dE/dx (a.u.)", speciesTitle[sp]).Data(),
                                        ptbins, ptlow, ptup, 1000, 0.0, 1000.0);
        fhPureTrackLengthA[sp] = new TH1F(TString::Format("fhPureTrackLengthA_%s", speciesName[sp]).Data(),
                                          TString::Format("Track Length of Pure %s; L (cm)", speciesTitle[sp]).Data(),
                                          1000, 0.0, 1000.0);
        fhPureTrackLengthTOFA[sp] = new TH1F(TString::Format("fhPureTrackLengthTOFA_%s", speciesName[sp]).Data(),
                                             TString::Format("Track Length of Pure %s with TOF; L (cm)", speciesTitle[sp]).Data(),
                                             1000, 0.0, 1000.0);
        fhPureTrackTimeA[sp] = new TH2F(TString::Format("fhPureTrackTimeA_%s", speciesName[sp]).Data(),
                                        TString::Format("Track Time vs P Pure %s; P (GeV/c); Track Time(ns)", speciesTitle[sp]).Data(),
                                        ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        fhPureTrackBetaInvA[sp] = new TH2F(TString::Format("fhPureTrackBetaInvA_%s", speciesName[sp]).Data(),
                                           TString::Format("1/#Beta vs P Pure %s; P (GeV/c); 1/#Beta(ns/m)", speciesTitle[sp]).Data(),
                                           ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        fhPureTrackTimeIPA[sp] = new TH2F(TString::Format("fhPureTrackTimeIPA_%s", speciesName[sp]).Data(),
                                          TString::Format("Track Time vs P_{IP} Pure %s; P (GeV/c); Track Time(ns)", speciesTitle[sp]).Data(),
                                          ptbins, ptlow, ptup, 1000, 0.0, 10.0);
        fhPureTrackBetaInvIPA[sp] = new TH2F(TString::Format("fhPureTrackBetaInvIPA_%s", speciesName[sp]).Data(),
                                             TString::Format("1/#Beta vs P_{IP} Pure %s; P (GeV/c); 1/#Beta(ns/m)", speciesTitle[sp]).Data(),
                                             ptbins, ptlow, ptup, 1000, 0.0, 10.0);
      }
      if (makeNSigmaPlots) {
        for (int sp1 = 0; sp1 < kIdBfNoOfSpecies; ++sp1) {
          for (int sp2 = 0; sp2 < kIdBfNoOfSpecies; ++sp2) {
            fhTrueNSigmaTPC[sp1][sp2] = new TH2F(TString::Format("fhTrueNSigmaTPC%s_%s", speciesName[sp1], speciesName[sp2]).Data(),
                                                 TString::Format("N #sigma %s from TPC vs P for generated %s;N #sigma;p (GeV/c)", speciesTitle[sp1], speciesTitle[sp2]).Data(),
                                                 48, -6, 6,
                                                 ptbins, ptlow, ptup);

            fhTrueNSigmaTOF[sp1][sp2] = new TH2F(TString::Format("fhTrueNSigmaTOF%s_%s", speciesName[sp1], speciesName[sp2]).Data(),
                                                 TString::Format("N #sigma %s from TOF vs P for generated %s;N #sigma;p (GeV/c)", speciesTitle[sp1], speciesTitle[sp2]).Data(),
                                                 48, -6, 6,
                                                 ptbins, ptlow, ptup);
            fhPrimaryNSigmaTPC[sp1][sp2] = new TH2F(TString::Format("fhPrimaryNSigmaTPC%s_%s", speciesName[sp1], speciesName[sp2]).Data(),
                                                    TString::Format("N #sigma %s from TPC vs P for primary %s;N #sigma;p (GeV/c)", speciesTitle[sp1], speciesTitle[sp2]).Data(),
                                                    48, -6, 6,
                                                    ptbins, ptlow, ptup);

            fhPrimaryNSigmaTOF[sp1][sp2] = new TH2F(TString::Format("fhPrimaryNSigmaTOF%s_%s", speciesName[sp1], speciesName[sp2]).Data(),
                                                    TString::Format("N #sigma %s from TOF vs P for primary %s;N #sigma;p (GeV/c)", speciesTitle[sp1], speciesTitle[sp2]).Data(),
                                                    48, -6, 6,
                                                    ptbins, ptlow, ptup);
          }
        }
      }

      /* add the hstograms to the output list */
      fOutputList->Add(fhTruePIDMismatch);
      fOutputList->Add(fhTruePIDCorrect);
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
      fOutputList->Add(fhPrimaryPB);
      fOutputList->Add(fhPrimaryPtB);
      fOutputList->Add(fhPrimarydEdxB);
      fOutputList->Add(fhPrimarydEdxIPTPCB);
      fOutputList->Add(fhPrimaryTrackLengthB);
      fOutputList->Add(fhPrimaryTrackLengthTOFB);
      fOutputList->Add(fhPrimaryTrackTimeB);
      fOutputList->Add(fhPrimaryTrackBetaInvB);
      fOutputList->Add(fhPrimaryTrackTimeIPB);
      fOutputList->Add(fhPrimaryTrackBetaInvIPB);
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
        fOutputList->Add(fhTruePtEtaPosA[sp]);
        fOutputList->Add(fhTruePtEtaNegA[sp]);
        fOutputList->Add(fhTrueNPosNegA[sp]);
        fOutputList->Add(fhTrueDeltaNA[sp]);

        fOutputList->Add(fhTruedEdxA[sp]);
        fOutputList->Add(fhTruedEdxIPTPCA[sp]);
        fOutputList->Add(fhTrueTrackLengthA[sp]);
        fOutputList->Add(fhTrueTrackLengthTOFA[sp]);
        fOutputList->Add(fhTrueTrackTimeA[sp]);
        fOutputList->Add(fhTrueTrackBetaInvA[sp]);
        fOutputList->Add(fhTrueTrackTimeIPA[sp]);
        fOutputList->Add(fhTrueTrackBetaInvIPA[sp]);

        fOutputList->Add(fhPrimaryPA[sp]);
        fOutputList->Add(fhPrimaryPtA[sp]);
        fOutputList->Add(fhPrimarydEdxA[sp]);
        fOutputList->Add(fhPrimarydEdxIPTPCA[sp]);
        fOutputList->Add(fhPrimaryTrackLengthA[sp]);
        fOutputList->Add(fhPrimaryTrackLengthTOFA[sp]);
        fOutputList->Add(fhPrimaryTrackTimeA[sp]);
        fOutputList->Add(fhPrimaryTrackBetaInvA[sp]);
        fOutputList->Add(fhPrimaryTrackTimeIPA[sp]);
        fOutputList->Add(fhPrimaryTrackBetaInvIPA[sp]);

        fOutputList->Add(fhPurePA[sp]);
        fOutputList->Add(fhPurePtA[sp]);
        fOutputList->Add(fhPuredEdxA[sp]);
        fOutputList->Add(fhPuredEdxIPTPCA[sp]);
        fOutputList->Add(fhPureTrackLengthA[sp]);
        fOutputList->Add(fhPureTrackLengthTOFA[sp]);
        fOutputList->Add(fhPureTrackTimeA[sp]);
        fOutputList->Add(fhPureTrackBetaInvA[sp]);
        fOutputList->Add(fhPureTrackTimeIPA[sp]);
        fOutputList->Add(fhPureTrackBetaInvIPA[sp]);
      }
      if (makeNSigmaPlots) {
        for (int sp1 = 0; sp1 < kIdBfNoOfSpecies; ++sp1) {
          for (int sp2 = 0; sp2 < kIdBfNoOfSpecies; ++sp2) {
            fOutputList->Add(fhTrueNSigmaTPC[sp1][sp2]);
            fOutputList->Add(fhTrueNSigmaTOF[sp1][sp2]);
            fOutputList->Add(fhPrimaryNSigmaTPC[sp1][sp2]);
            fOutputList->Add(fhPrimaryNSigmaTOF[sp1][sp2]);
          }
        }
      }
    }
    /* initialize access to the CCDB */
  }

  template <typename TrackObject>
  inline MatchRecoGenSpecies identifyTrack(TrackObject const& track);
  template <typename TrackObject>
  int8_t acceptTrack(TrackObject const& track);
  template <typename ParticleObject, typename MCCollisionObject>
  int8_t acceptParticle(ParticleObject& particle, MCCollisionObject const& mccollision);
  template <typename CollisionObjects, typename TrackObject>
  int8_t selectTrackAmbiguousCheck(CollisionObjects const& collisions, TrackObject const& track);
  template <typename ParticleObject>
  inline void identifyPIDMismatch(ParticleObject const& particle, MatchRecoGenSpecies const& trkId);
  template <typename ParticleObject>
  inline void identifyRealNSigma(ParticleObject const& particle, std::vector<float> tpcNSigma, std::vector<float> tofNSigma, float tpcInnerParam);
  template <typename ParticleObject>
  inline MatchRecoGenSpecies identifyParticle(ParticleObject const& particle);
  template <typename TrackObject>
  void fillTrackHistosBeforeSelection(TrackObject const& track);
  template <typename TrackObject>
  void fillTrackHistosAfterSelection(TrackObject const& track, MatchRecoGenSpecies sp);
  template <typename ParticleObject>
  bool isPrimary(ParticleObject const& particle);
  template <typename TrackObject>
  void fillRealPIDTrackHistosAfter(TrackObject const& track, MatchRecoGenSpecies sp);
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
    // LOGF(info, "Top of filterTracks");
    int naccepted = 0;
    int ncollaccepted = 0;
    if (!fullDerivedData) {
      tracksinfo.reserve(tracks.size());
    }
    for (const auto& collision : collisions) {
      if (collision.collisionaccepted()) {
        ncollaccepted++;
      }
    }
    for (const auto& track : tracks) {
      int8_t pid = -1;
      if (track.has_collision() && (track.template collision_as<soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>>()).collisionaccepted()) {
        pid = selectTrackAmbiguousCheck(collisions, track);
        if (!(pid < 0)) {
          naccepted++;
          /* update charged multiplicities */
          if (pid % twoDenom == trackTypes[0]) {
            trkMultPos[kIdBfCharged]++;
          }
          if (pid % twoDenom == trackTypes[1]) {
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

    for (const auto& gencoll : gencollisions) {
      if (gencoll.collisionaccepted()) {
        acceptedcollisions++;
      }
    }

    for (const auto& particle : particles) {
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
            pid = acceptParticle(particle, mccollision);
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

  void filterDetectorLevelWithPID(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksPIDDetLevel const& tracks, aod::McParticles const&)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterDetectorLevelWithPID, "Not stored derived data detector level track filtering", false)

  void filterDetectorLevelWithPIDAmbiguous(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksPIDDetLevelAmbiguous const& tracks, aod::McParticles const&)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterDetectorLevelWithPIDAmbiguous, "Not stored derived data detector level track filtering with ambiguous tracks check", false)

  void filterRecoWithFullPID(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksFullPID const& tracks, aod::McParticles const&)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterRecoWithFullPID, "Not stored derived data track filtering", false)

  void filterRecoWithFullPIDAmbiguous(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksFullPIDAmbiguous const& tracks)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterRecoWithFullPIDAmbiguous, "Not stored derived data track filtering with ambiguous tracks check", false)

  void filterDetectorLevelWithFullPID(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksFullPIDDetLevel const& tracks, aod::McParticles const&)
  {
    filterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(IdentifiedBfFilterTracks, filterDetectorLevelWithFullPID, "Not stored derived data detector level track filtering", false)

  void filterDetectorLevelWithFullPIDAmbiguous(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo>& collisions, IdBfFullTracksFullPIDDetLevelAmbiguous const& tracks, aod::McParticles const&)
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

  void filterDetectorLevelWithoutPIDAmbiguous(soa::Join<aod::Collisions, aod::IdentifiedBfCFCollisionsInfo> const& collisions, IdBfFullTracksDetLevelAmbiguous const& tracks, aod::McParticles const&)
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
inline MatchRecoGenSpecies IdentifiedBfFilterTracks::identifyParticle(ParticleObject const& particle)
{
  using namespace identifiedbffilter;
  // LOGF(info, "Top of identifyParticle");

  int pdgcode = std::fabs(particle.pdgCode());

  switch (pdgcode) {
    case kPositron:
      return kIdBfElectron;
      break;

    case kPiPlus:
      return kIdBfPion;
      break;
    case kKPlus:
      return kIdBfKaon;
      break;
    case kProton:
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

template <typename ParticleObject>
inline void IdentifiedBfFilterTracks::identifyPIDMismatch(ParticleObject const& particle, MatchRecoGenSpecies const& trkId)
{
  MatchRecoGenSpecies realPID = identifyParticle(particle);
  if (!(realPID < 0)) {
    if (realPID == trkId) {
      fhTruePIDCorrect->Fill(realPID);
    } else {
      fhTruePIDMismatch->Fill(realPID, trkId);
    }
  }
}

template <typename ParticleObject>
inline void IdentifiedBfFilterTracks::identifyRealNSigma(ParticleObject const& particle, std::vector<float> tpcNSigma, std::vector<float> tofNSigma, float tpcInnerParam)
{

  MatchRecoGenSpecies realPID = identifyParticle(particle);
  if (!(realPID < 0)) {
    fhTrueNSigmaTPC[kIdBfElectron][realPID]->Fill(tpcNSigma[kIdBfElectron], tpcInnerParam);
    fhTrueNSigmaTOF[kIdBfElectron][realPID]->Fill(tofNSigma[kIdBfElectron], tpcInnerParam);
    fhTrueNSigmaTPC[kIdBfPion][realPID]->Fill(tpcNSigma[kIdBfPion], tpcInnerParam);
    fhTrueNSigmaTOF[kIdBfPion][realPID]->Fill(tofNSigma[kIdBfPion], tpcInnerParam);
    fhTrueNSigmaTPC[kIdBfKaon][realPID]->Fill(tpcNSigma[kIdBfKaon], tpcInnerParam);
    fhTrueNSigmaTOF[kIdBfKaon][realPID]->Fill(tofNSigma[kIdBfKaon], tpcInnerParam);
    fhTrueNSigmaTPC[kIdBfProton][realPID]->Fill(tpcNSigma[kIdBfProton], tpcInnerParam);
    fhTrueNSigmaTOF[kIdBfProton][realPID]->Fill(tofNSigma[kIdBfProton], tpcInnerParam);

    if (particle.isPhysicalPrimary()) {
      fhPrimaryNSigmaTPC[kIdBfElectron][realPID]->Fill(tpcNSigma[kIdBfElectron], tpcInnerParam);
      fhPrimaryNSigmaTOF[kIdBfElectron][realPID]->Fill(tofNSigma[kIdBfElectron], tpcInnerParam);
      fhPrimaryNSigmaTPC[kIdBfPion][realPID]->Fill(tpcNSigma[kIdBfPion], tpcInnerParam);
      fhPrimaryNSigmaTOF[kIdBfPion][realPID]->Fill(tofNSigma[kIdBfPion], tpcInnerParam);
      fhPrimaryNSigmaTPC[kIdBfKaon][realPID]->Fill(tpcNSigma[kIdBfKaon], tpcInnerParam);
      fhPrimaryNSigmaTOF[kIdBfKaon][realPID]->Fill(tofNSigma[kIdBfKaon], tpcInnerParam);
      fhPrimaryNSigmaTPC[kIdBfProton][realPID]->Fill(tpcNSigma[kIdBfProton], tpcInnerParam);
      fhPrimaryNSigmaTOF[kIdBfProton][realPID]->Fill(tofNSigma[kIdBfProton], tpcInnerParam);
    }
  }
}

template <typename TrackObject>
void fillNSigmaHistos(TrackObject const& track)
{

  float actualTPCNSigma[kIdBfNoOfSpecies];

  actualTPCNSigma[kIdBfElectron] = track.tpcNSigmaEl();
  actualTPCNSigma[kIdBfPion] = track.tpcNSigmaPi();
  actualTPCNSigma[kIdBfKaon] = track.tpcNSigmaKa();
  actualTPCNSigma[kIdBfProton] = track.tpcNSigmaPr();

  float actualTOFNSigma[kIdBfNoOfSpecies];

  actualTOFNSigma[kIdBfElectron] = track.tofNSigmaEl();
  actualTOFNSigma[kIdBfPion] = track.tofNSigmaPi();
  actualTOFNSigma[kIdBfKaon] = track.tofNSigmaKa();
  actualTOFNSigma[kIdBfProton] = track.tofNSigmaPr();

  if (loadfromccdb) {
    actualTPCNSigma[kIdBfElectron] = actualTPCNSigma[kIdBfElectron] - fhNSigmaCorrection[kIdBfElectron]->GetBinContent(fhNSigmaCorrection[kIdBfElectron]->FindBin(track.tpcInnerParam()));
    actualTPCNSigma[kIdBfPion] = actualTPCNSigma[kIdBfPion] - fhNSigmaCorrection[kIdBfPion]->GetBinContent(fhNSigmaCorrection[kIdBfPion]->FindBin(track.tpcInnerParam()));
    actualTPCNSigma[kIdBfKaon] = actualTPCNSigma[kIdBfKaon] - fhNSigmaCorrection[kIdBfKaon]->GetBinContent(fhNSigmaCorrection[kIdBfKaon]->FindBin(track.tpcInnerParam()));
    actualTPCNSigma[kIdBfProton] = actualTPCNSigma[kIdBfProton] - fhNSigmaCorrection[kIdBfProton]->GetBinContent(fhNSigmaCorrection[kIdBfProton]->FindBin(track.tpcInnerParam()));
  }

  fhNSigmaTPC[kIdBfElectron]->Fill(actualTPCNSigma[kIdBfElectron], track.tpcInnerParam());
  fhNSigmaTPC[kIdBfPion]->Fill(actualTPCNSigma[kIdBfPion], track.tpcInnerParam());
  fhNSigmaTPC[kIdBfKaon]->Fill(actualTPCNSigma[kIdBfKaon], track.tpcInnerParam());
  fhNSigmaTPC[kIdBfProton]->Fill(actualTPCNSigma[kIdBfProton], track.tpcInnerParam());

  fhNSigmaTOF[kIdBfElectron]->Fill(actualTOFNSigma[kIdBfElectron], track.tpcInnerParam());
  fhNSigmaTOF[kIdBfPion]->Fill(actualTOFNSigma[kIdBfPion], track.tpcInnerParam());
  fhNSigmaTOF[kIdBfKaon]->Fill(actualTOFNSigma[kIdBfKaon], track.tpcInnerParam());
  fhNSigmaTOF[kIdBfProton]->Fill(actualTOFNSigma[kIdBfProton], track.tpcInnerParam());

  fhNSigmaCombo[kIdBfElectron]->Fill(sqrtf(actualTOFNSigma[kIdBfElectron] * actualTOFNSigma[kIdBfElectron] + actualTPCNSigma[kIdBfElectron] * actualTPCNSigma[kIdBfElectron]), track.tpcInnerParam());
  fhNSigmaCombo[kIdBfPion]->Fill(sqrtf(actualTOFNSigma[kIdBfPion] * actualTOFNSigma[kIdBfPion] + actualTPCNSigma[kIdBfPion] * actualTPCNSigma[kIdBfPion]), track.tpcInnerParam());
  fhNSigmaCombo[kIdBfKaon]->Fill(sqrtf(actualTOFNSigma[kIdBfKaon] * actualTOFNSigma[kIdBfKaon] + actualTPCNSigma[kIdBfKaon] * actualTPCNSigma[kIdBfKaon]), track.tpcInnerParam());
  fhNSigmaCombo[kIdBfProton]->Fill(sqrtf(actualTOFNSigma[kIdBfProton] * actualTOFNSigma[kIdBfProton] + actualTPCNSigma[kIdBfProton] * actualTPCNSigma[kIdBfProton]), track.tpcInnerParam());
}

/// \brief Identifies the passed track with TPC and TOF data
/// \param track the track of interest
/// \return the internal track id, -1 if not accepted

template <typename TrackObject>
inline MatchRecoGenSpecies IdentifiedBfFilterTracks::identifyTrack(TrackObject const& track)
{
  using namespace o2::analysis::identifiedbffilter;

  fillNSigmaHistos(track);

  std::vector<float> actualTPCNSigma(kIdBfNoOfSpecies, 0.);

  actualTPCNSigma[kIdBfElectron] = track.tpcNSigmaEl();
  actualTPCNSigma[kIdBfPion] = track.tpcNSigmaPi();
  actualTPCNSigma[kIdBfKaon] = track.tpcNSigmaKa();
  actualTPCNSigma[kIdBfProton] = track.tpcNSigmaPr();

  std::vector<float> actualTOFNSigma(kIdBfNoOfSpecies, 0.);

  actualTOFNSigma[kIdBfElectron] = track.tofNSigmaEl();
  actualTOFNSigma[kIdBfPion] = track.tofNSigmaPi();
  actualTOFNSigma[kIdBfKaon] = track.tofNSigmaKa();
  actualTOFNSigma[kIdBfProton] = track.tofNSigmaPr();

  float nsigmas[kIdBfNoOfSpecies];

  if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
    if (makeNSigmaPlots) {
      identifyRealNSigma(track.template mcParticle_as<aod::McParticles>(), actualTPCNSigma, actualTOFNSigma, track.tpcInnerParam());
    }
  }

  if (loadfromccdb) {
    for (int iSp = 0; iSp < kIdBfNoOfSpecies; iSp++) {
      actualTPCNSigma[iSp] = actualTPCNSigma[iSp] - fhNSigmaCorrection[iSp]->GetBinContent(fhNSigmaCorrection[iSp]->FindBin(track.tpcInnerParam()));
    }
  }

  for (int iSp = 0; iSp < kIdBfNoOfSpecies; iSp++) {

    if (track.tpcInnerParam() < tofCut[iSp] && track.tpcInnerParam() < tpcCut[iSp] && !onlyTOF) {
      nsigmas[iSp] = actualTPCNSigma[iSp];
    } else if (track.tpcInnerParam() > tofCut[iSp] && track.tpcInnerParam() < tpcCut[iSp] && !onlyTOF && track.hasTOF()) {
      nsigmas[iSp] = sqrtf(actualTPCNSigma[iSp] * actualTPCNSigma[iSp] + actualTOFNSigma[iSp] * actualTOFNSigma[iSp]);
    } else if (track.hasTOF() && ((track.tpcInnerParam() > tofCut[iSp] && track.tpcInnerParam() > tpcCut[iSp]) || onlyTOF)) {
      nsigmas[iSp] = actualTOFNSigma[iSp];
    } else {
      return kWrongSpecies;
    }
  }

  float minNSigma = 999.0f;
  MatchRecoGenSpecies spMinNSigma = kWrongSpecies;
  for (int sp = 0; sp < kIdBfNoOfSpecies; ++sp) {
    if (doPID[sp]) {                                       // Check if we're IDing PID for this species
      if (std::fabs(nsigmas[sp]) < std::fabs(minNSigma)) { // Check if species nsigma is less than current nsigma
        minNSigma = nsigmas[sp];                           // If yes, set species nsigma to current nsigma
        spMinNSigma = MatchRecoGenSpecies(sp);             // set current species sp number to current sp
      }
    }
  }
  bool doublematch = false;
  MatchRecoGenSpecies spDouble = kWrongSpecies;
  // LOGF(info,"Looking at accept range");
  if (minNSigma < acceptRange[spMinNSigma][1] && minNSigma > acceptRange[spMinNSigma][0]) { // Check that current nsigma is in accpetance range
    // LOGF(info,"In accept Range");
    for (int sp = 0; (sp < kIdBfNoOfSpecies) && !doublematch; ++sp) { // iterate over all species while there's no double match and we're in the list
      if (sp != spMinNSigma) {                                        // for species not current minimum nsigma species
        // LOGF(info, "looking at Reject Range");
        if (std::fabs(nsigmas[sp]) < rejectRange[spMinNSigma][sp]) {  // If secondary species is in rejection range
          doublematch = true;                                         // Set double match true
          spDouble = MatchRecoGenSpecies(sp);
        }
      }
    }
    if (doublematch) { // if double match true
      fhWrongTrackID->Fill(track.p());
      fhdEdxA[kIdBfNoOfSpecies]->Fill(track.p(), track.tpcSignal());
      fhdEdxIPTPCA[kIdBfNoOfSpecies]->Fill(track.tpcInnerParam(), track.tpcSignal());
      fhTrackLengthA[kIdBfNoOfSpecies]->Fill(track.length());
      fhTrackTimeA[kIdBfNoOfSpecies]->Fill(track.p(), (track.trackTime()));
      fhTrackTimeIPA[kIdBfNoOfSpecies]->Fill(track.tpcInnerParam(), (track.trackTime()));

      if constexpr (framework::has_type_v<aod::pidtofbeta::Beta, typename TrackObject::all_columns>) {
        fhTrackBetaInvA[kIdBfNoOfSpecies]->Fill(track.p(), 1 / track.beta());
        fhTrackBetaInvIPA[kIdBfNoOfSpecies]->Fill(track.tpcInnerParam(), 1 / track.beta());
      }
      fhDoublePID->Fill(spMinNSigma, spDouble);
      return kWrongSpecies; // Return wrong species value
    } else {
      fhNSigmaTPCIdTrks[spMinNSigma]->Fill(actualTPCNSigma[spMinNSigma], track.tpcInnerParam());

      if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
        identifyPIDMismatch(track.template mcParticle_as<aod::McParticles>(), spMinNSigma);
      }
      return spMinNSigma;
    }
  } else {
    return kWrongSpecies;
  }
}

/// \brief Accepts or not the passed track
/// \param track the track of interest
/// \return the internal track id, -1 if not accepted

template <typename TrackObject>
inline int8_t IdentifiedBfFilterTracks::acceptTrack(TrackObject const& track)
{
  // LOGF(info,"Top of acceptTrack");
  fillTrackHistosBeforeSelection(track); // <Fill "before selection" histo

  /* TODO: incorporate a mask in the scanned tracks table for the rejecting track reason */
  if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
    if (track.mcParticleId() < 0) {
      // LOGF(info,"No matching MC particle");
      return -1;
    }
  }

  if (matchTrackType(track)) {
    // LOGF(info, "Track type match");
    if (ptlow < track.pt() && track.pt() < ptup && etalow < track.eta() && track.eta() < etaup) {
      // LOGF(info, "Track Accepted");
      fillTrackHistosAfterSelection(track, kIdBfCharged);
      MatchRecoGenSpecies sp = kWrongSpecies;
      if (recoIdMethod == recoIdMethods[0]) {
        sp = kIdBfCharged;
      } else if (recoIdMethod == recoIdMethods[1]) {
        if constexpr (framework::has_type_v<aod::pidtpc_tiny::TPCNSigmaStorePi, typename TrackObject::all_columns> || framework::has_type_v<aod::pidtpc::TPCNSigmaPi, typename TrackObject::all_columns>) {
          sp = identifyTrack(track);
        } else {
          LOGF(fatal, "Track identification required but PID information not present");
        }
      } else if (recoIdMethod == recoIdMethods[2]) {
        if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
          sp = identifyParticle(track.template mcParticle_as<aod::McParticles>());
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
inline int8_t IdentifiedBfFilterTracks::acceptParticle(ParticleObject& particle, MCCollisionObject const& mccollision)
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
  if ((particle.flags() & 0x8) != 0x8) {
    if (particle.isPhysicalPrimary() && std::fabs(charge) > 0.0) {
      if ((particle.mcCollisionId() == 0) && traceCollId0) {
        LOGF(info, "Particle %d passed isPhysicalPrimary", particle.globalIndex());
      }

      if (ptlow < particle.pt() && particle.pt() < ptup && etalow < particle.eta() && particle.eta() < etaup) {
        MatchRecoGenSpecies sp = kWrongSpecies;
        if (recoIdMethod == recoIdMethods[0]) {
          sp = kIdBfCharged;
        }
        if (recoIdMethod == recoIdMethods[1] || recoIdMethod == recoIdMethods[2]) {
          sp = identifyParticle(particle);
        }

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
  } else {
    if ((particle.mcCollisionId() == 0) && traceCollId0) {
      LOGF(info, "Particle %d Out of Bunch Pileup", particle.globalIndex());
    }
  }
  return kWrongSpecies;
}

template <typename CollisionObjects, typename TrackObject>
int8_t IdentifiedBfFilterTracks::selectTrackAmbiguousCheck(CollisionObjects const& collisions, TrackObject const& track)
{
  // LOGF(info,"Top of AmbiguousCheck");
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
    if (tracktype == trackTypes[2]) {
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
    return acceptTrack(track);
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

  if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {

    if (isPrimary(track.template mcParticle_as<aod::McParticles>())) {
      fhPrimaryPB->Fill(track.p());
      fhPrimaryPtB->Fill(track.pt());
      fhPrimarydEdxB->Fill(track.p(), track.tpcSignal());
      fhPrimarydEdxIPTPCB->Fill(track.tpcInnerParam(), track.tpcSignal());
      fhPrimaryTrackLengthB->Fill(track.length());
      if (track.hasTOF() && track.p() > tofCut[0]) {
        fhPrimaryTrackLengthTOFB->Fill(track.length());
        fhPrimaryTrackTimeB->Fill(track.p(), (track.trackTime()));
        fhPrimaryTrackTimeIPB->Fill(track.tpcInnerParam(), (track.trackTime()));

        if constexpr (framework::has_type_v<aod::pidtofbeta::Beta, typename TrackObject::all_columns>) {
          fhPrimaryTrackBetaInvB->Fill(track.p(), 1 / track.beta());
          fhPrimaryTrackBetaInvIPB->Fill(track.tpcInnerParam(), 1 / track.beta());
        }
      }
    }
  }
}
template <typename ParticleObject>
bool IdentifiedBfFilterTracks::isPrimary(ParticleObject const& particle)
{
  return particle.isPhysicalPrimary();
}
template <typename TrackObject>
void IdentifiedBfFilterTracks::fillRealPIDTrackHistosAfter(TrackObject const& track, MatchRecoGenSpecies sp)
{

  if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
    MatchRecoGenSpecies realPID = identifyParticle(track.template mcParticle_as<aod::McParticles>());
    if (!(realPID < 0)) {
      fhTruedEdxA[realPID]->Fill(track.p(), track.tpcSignal());
      fhTruedEdxIPTPCA[realPID]->Fill(track.tpcInnerParam(), track.tpcSignal());
      fhTrueTrackLengthA[realPID]->Fill(track.length());
      if (track.hasTOF() && track.p() > tofCut[realPID]) {
        fhTrueTrackLengthTOFA[realPID]->Fill(track.length());
        fhTrueTrackTimeA[realPID]->Fill(track.p(), (track.trackTime()));
        fhTrueTrackTimeIPA[realPID]->Fill(track.tpcInnerParam(), (track.trackTime()));
        if constexpr (framework::has_type_v<aod::pidtofbeta::Beta, typename TrackObject::all_columns>) {
          fhTrueTrackBetaInvA[realPID]->Fill(track.p(), 1 / track.beta());
          fhTrueTrackBetaInvIPA[realPID]->Fill(track.tpcInnerParam(), 1 / track.beta());
        }
      }
    }

    if (isPrimary(track.template mcParticle_as<aod::McParticles>())) {
      fhPrimaryPA[sp]->Fill(track.p());
      fhPrimaryPtA[sp]->Fill(track.pt());
      fhPrimarydEdxA[sp]->Fill(track.p(), track.tpcSignal());
      fhPrimarydEdxIPTPCA[sp]->Fill(track.tpcInnerParam(), track.tpcSignal());
      fhPrimaryTrackLengthA[sp]->Fill(track.length());
      if (track.hasTOF() && track.p() > tofCut[sp]) {
        fhPrimaryTrackLengthTOFA[sp]->Fill(track.length());
        fhPrimaryTrackTimeA[sp]->Fill(track.p(), (track.trackTime()));
        fhPrimaryTrackTimeIPA[sp]->Fill(track.tpcInnerParam(), (track.trackTime()));

        if constexpr (framework::has_type_v<aod::pidtofbeta::Beta, typename TrackObject::all_columns>) {
          fhPrimaryTrackBetaInvA[sp]->Fill(track.p(), 1 / track.beta());
          fhPrimaryTrackBetaInvIPA[sp]->Fill(track.tpcInnerParam(), 1 / track.beta());
        }
      }

      if (sp == realPID) {
        fhPurePA[realPID]->Fill(track.p());
        fhPurePtA[realPID]->Fill(track.pt());
        fhPuredEdxA[realPID]->Fill(track.p(), track.tpcSignal());
        fhPuredEdxIPTPCA[realPID]->Fill(track.tpcInnerParam(), track.tpcSignal());
        fhPureTrackLengthA[realPID]->Fill(track.length());
        if (track.hasTOF() && track.p() > tofCut[realPID]) {
          fhPureTrackLengthTOFA[realPID]->Fill(track.length());
          fhPureTrackTimeA[realPID]->Fill(track.p(), (track.trackTime()));
          fhPureTrackTimeIPA[realPID]->Fill(track.tpcInnerParam(), (track.trackTime()));

          if constexpr (framework::has_type_v<aod::pidtofbeta::Beta, typename TrackObject::all_columns>) {
            fhPureTrackBetaInvA[realPID]->Fill(track.p(), 1 / track.beta());
            fhPureTrackBetaInvIPA[realPID]->Fill(track.tpcInnerParam(), 1 / track.beta());
          }
        }
      }
    }
  }
}
template <typename TrackObject>
void IdentifiedBfFilterTracks::fillTrackHistosAfterSelection(TrackObject const& track, MatchRecoGenSpecies sp)
{
  /* the charged species should have been called first so avoid double counting */
  // LOGF(info,"Top of AfterSelection");
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
  fhTrackLengthA[sp]->Fill(track.length());
  if (track.hasTOF() && track.p() > tofCut[sp]) {
    fhTrackLengthTOFA[sp]->Fill(track.length());
    fhTrackTimeA[sp]->Fill(track.p(), (track.trackTime()));
    fhTrackTimeIPA[sp]->Fill(track.tpcInnerParam(), (track.trackTime()));

    if constexpr (framework::has_type_v<aod::pidtofbeta::Beta, typename TrackObject::all_columns>) {
      fhTrackBetaInvA[sp]->Fill(track.p(), 1 / track.beta());
      fhTrackBetaInvIPA[sp]->Fill(track.tpcInnerParam(), 1 / track.beta());
    }
  }
  if (track.sign() > 0) {
    fhPtPosA[sp]->Fill(track.pt());
    fhPtEtaPosA[sp]->Fill(track.pt(), track.eta());
  } else {
    fhPtNegA[sp]->Fill(track.pt());
    fhPtEtaNegA[sp]->Fill(track.pt(), track.eta());
  }
  if ((fDataType != kData) && (fDataType != kDataNoEvtSel)) {
    fillRealPIDTrackHistosAfter(track, sp);
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

  float dcaxy = std::sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                          (particle.vy() - collision.posY()) * (particle.vy() - collision.posY()));
  if (traceDCAOutliers.mDoIt && (traceDCAOutliers.mLowValue < dcaxy) && (dcaxy < traceDCAOutliers.mUpValue)) {
    fhTrueDCAxyBid->Fill(TString::Format("%d", particle.pdgCode()).Data(), 1.0);
  }

  fhTrueDCAxyB->Fill(std::sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                               (particle.vy() - collision.posY()) * (particle.vy() - collision.posY())));
  fhTrueDCAxyzB->Fill(std::sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
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
    float dcaxy = std::sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                            (particle.vy() - collision.posY()) * (particle.vy() - collision.posY()));
    if (traceDCAOutliers.mDoIt && (traceDCAOutliers.mLowValue < dcaxy) && (dcaxy < traceDCAOutliers.mUpValue)) {
      LOGF(info, "DCAxy outlier: Particle with index %d and pdg code %d assigned to MC collision %d, pT: %f, phi: %f, eta: %f",
           particle.globalIndex(), particle.pdgCode(), particle.mcCollisionId(), particle.pt(), particle.phi(), particle.eta());
      LOGF(info, "               With status %d and flags %0X", particle.statusCode(), particle.flags());
    }

    fhTrueDCAxyA->Fill(std::sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                                 (particle.vy() - collision.posY()) * (particle.vy() - collision.posY())));
    fhTrueDCAzA->Fill((particle.vz() - collision.posZ()));
    fhTrueDCAxyzA->Fill(std::sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                                  (particle.vy() - collision.posY()) * (particle.vy() - collision.posY())),
                        (particle.vz() - collision.posZ()));
  }
  fhTruePA[sp]->Fill(particle.p());
  fhTruePtA[sp]->Fill(particle.pt());
  if (charge > 0) {
    fhTruePtPosA[sp]->Fill(particle.pt());
    fhTruePtEtaPosA[sp]->Fill(particle.pt(), particle.eta());
  } else {
    fhTruePtNegA[sp]->Fill(particle.pt());
    fhTruePtEtaNegA[sp]->Fill(particle.pt(), particle.eta());
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
