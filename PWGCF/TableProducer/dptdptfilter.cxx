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

/// \file dptdptfilter.cxx
/// \brief Filters collisions and tracks according to selection criteria
/// \author victor.gonzalez.sebastian@gmail.com

#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "CommonConstants/PhysicsConstants.h"
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
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include <TROOT.h>
#include <TPDGCode.h>
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
using DptDptTracksFullPID = soa::Join<aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using DptDptFullTracksPID = soa::Join<DptDptFullTracks, DptDptTracksPID>;
using DptDptFullTracksPIDAmbiguous = soa::Join<DptDptFullTracksAmbiguous, DptDptTracksPID>;
using DptDptFullTracksFullPID = soa::Join<DptDptFullTracks, DptDptTracksFullPID>;
using DptDptFullTracksFullPIDAmbiguous = soa::Join<DptDptFullTracksAmbiguous, DptDptTracksFullPID>;
using DptDptFullTracksDetLevel = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using DptDptFullTracksDetLevelAmbiguous = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::TrackCompColls>;
using DptDptFullTracksPIDDetLevel = soa::Join<DptDptFullTracksDetLevel, DptDptTracksPID>;
using DptDptFullTracksPIDDetLevelAmbiguous = soa::Join<DptDptFullTracksDetLevelAmbiguous, DptDptTracksPID>;
using DptDptFullTracksFullPIDDetLevel = soa::Join<DptDptFullTracksDetLevel, DptDptTracksFullPID>;
using DptDptFullTracksFullPIDDetLevelAmbiguous = soa::Join<DptDptFullTracksDetLevelAmbiguous, DptDptTracksFullPID>;

bool fullDerivedData = false; /* produce full derived data for its external storage */
TpcExcludeTrack tpcExcluder;  ///< the TPC excluder object instance

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

const char* eventSelectionSteps[knCollisionSelectionFlags] = {
  "IN",
  "MB",
  "INT7",
  "SEL7",
  "SEL8",
  "NOSAMEBUNCHPUP",
  "ISGOODZVTXFT0VSPV",
  "ISVERTEXITSTPC",
  "ISVERTEXTOFMATCHED",
  "ISVERTEXTRDMATCHED",
  "NOCOLLINTIMERANGE",
  "NOCOLLINROF",
  "OCCUPANCY",
  "ISGOODITSLAYER3",
  "ISGOODITSLAYER0123",
  "ISGOODITSLAYERALL",
  "CENTRALITY",
  "ZVERTEX",
  "SELECTED"};

//============================================================================================
// The DptDptFilter histogram objects
// TODO: consider registering in the histogram registry
//============================================================================================
TH1D* fhEventSelection = nullptr;
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

//////////////////////////////////////////////////////////////////////////////
// Multiplicity in principle for on the fly generated events
//////////////////////////////////////////////////////////////////////////////

struct Multiplicity {
  enum MultEst {
    kV0M,
    kCL1,
    kCL1GAP
  };

  float getMultiplicityClass() { return multiplicityClass; }
  float getMultiplicity() { return multiplicity; }

  MultEst classestimator = kV0M;

  float multiplicityClass = -1.0;
  float multiplicity = 0.0;
  bool inelgth0 = false;
  int v0am = 0;
  int v0cm = 0;
  int cl1m = 0;
  int cl1EtaGapM = 0;
  int dNchdEta = 0;
  int nPart = 0;
  TH1F* fhNPartTot = nullptr;           ///< total number of particles analyzed
  TH1F* fhMultiplicity;                 ///< the multiplicity distribution
  TH2F* fhV0Multiplicity;               ///< the V0M multiplicity histogram
  TH2F* fhCL1Multiplicity;              ///< the CL1 multiplicity histogram
  TH2F* fhCL1EtaGapMultiplicity;        ///< the CL1 with an eta gap multiplicity histogram
  const TH1* fhV0MMultPercentile;       ///< the V0M Centrality / Multiplicity percentile estimation histogram
  const TH1* fhCL1MultPercentile;       ///< the CL1 Centrality / Multiplicity percentile estimation histogram
  const TH1* fhCL1EtaGapMultPercentile; ///< the CL1 with an eta gap Centrality / Multiplicity percentile estimation histogram

  void init(TList* hlist)
  {
    fhNPartTot = new TH1F("CollisionNpart", "Collision analyzed particles;number of particles;counts", 8000, -0.5, 8000 - 0.5);
    fhMultiplicity = new TH1F("CollisionMultiplicity", "Event multiplicity;multiplicity (%);counts", 101, -0.5, 101 - 0.5);
    fhV0Multiplicity = new TH2F("V0Multiplicity", "V0M;V0M;d#it{N}/d#eta;counts", 3000, -9.5, 3000 - 9.5, 2500, -9.5, 2500 - 9.5);
    fhCL1Multiplicity = new TH2F("CL1Multiplicity", "CL1M;CL1M;d#it{N}/d#eta;counts", 3000, -9.5, 3000 - 9.5, 2500, -9.5, 2500 - 9.5);
    fhCL1EtaGapMultiplicity = new TH2F("CL1EtaGapMultiplicity", "CL1M (excl |#eta|<0.8);CL1M;d#it{N}/d#eta;counts", 3000, -9.5, 3000 - 9.5, 2500, -9.5, 2500 - 9.5);

    hlist->Add(fhNPartTot);
    hlist->Add(fhMultiplicity);
    hlist->Add(fhV0Multiplicity);
    hlist->Add(fhCL1Multiplicity);
    hlist->Add(fhCL1EtaGapMultiplicity);
  }

  void setMultiplicityPercentiles(TList* list)
  {
    LOGF(info, "setMultiplicityPercentiles()", "From list %s", list->GetName());
    fhV0MMultPercentile = reinterpret_cast<TH1*>(list->FindObject("V0MCentMult"));
    fhCL1MultPercentile = reinterpret_cast<TH1*>(list->FindObject("CL1MCentMult"));
    fhCL1EtaGapMultPercentile = reinterpret_cast<TH1*>(list->FindObject("CL1EtaGapMCentMult"));

    if (fhV0MMultPercentile == nullptr || fhCL1MultPercentile == nullptr || fhCL1EtaGapMultPercentile == nullptr) {
      LOGF(fatal, "setMultiplicityPercentiles()", "Percentiles histograms not correctly loaded. ABORTING!!!");
      return;
    }
  }

  template <typename Particle>
  bool addParticleToMultiplicity(const Particle& p)
  {
    /* on the fly MC production */
    /* get event multiplicity according to the passed eta range  */
    /* event multiplicity as number of primary charged particles */
    /* based on AliAnalysisTaskPhiCorrelations implementation    */
    int pdgcode = std::abs(p.pdgCode());
    auto addTo = [](const Particle& p, int& est, float etamin, float etamax) {
      if (p.eta() < etamax && etamin < p.eta()) {
        est = est + 1;
      }
    };

    /* pdg checks */
    switch (pdgcode) {
      case kPiPlus:
      case kKPlus:
      case kProton:
        /* not clear if we should use IsPhysicalPrimary here */
        /* TODO: adapt to FT0M Run 3 and other estimators */
        if (0.001 < p.pt() && p.pt() < 50.0) {
          if (p.eta() < 1.0 && -1.0 < p.eta()) {
            inelgth0 = true;
          }
          addTo(p, v0am, 2.8, 5.1);
          addTo(p, v0cm, -3.7, -1.7);
          addTo(p, cl1m, -1.4, 1.4);
          addTo(p, cl1EtaGapM, -1.4, -0.8);
          addTo(p, cl1EtaGapM, 0.8, 1.4);
          addTo(p, dNchdEta, -0.5, 0.5);
          nPart++;
        }
        break;
      default:
        break;
    }
    return true;
  }

  template <typename CollisionParticles>
  void extractMultiplicity(const CollisionParticles& particles)
  {
    multiplicityClass = 105;
    multiplicity = 0;
    inelgth0 = false;
    nPart = 0;
    v0am = 0;
    v0cm = 0;
    cl1m = 0;
    cl1EtaGapM = 0;
    dNchdEta = 0;

    for (auto const& particle : particles) {
      addParticleToMultiplicity(particle);
    }

    if (inelgth0) {
      if (fhNPartTot != nullptr) {
        fhNPartTot->Fill(nPart);
      }
      if (fhV0Multiplicity != nullptr) {
        fhV0Multiplicity->Fill(v0am + v0cm, dNchdEta);
      }
      if (fhCL1Multiplicity != nullptr) {
        fhCL1Multiplicity->Fill(cl1m, dNchdEta);
      }
      if (fhCL1EtaGapMultiplicity != nullptr) {
        fhCL1EtaGapMultiplicity->Fill(cl1EtaGapM, dNchdEta);
      }
      switch (classestimator) {
        case kV0M:
          if (fhV0MMultPercentile != nullptr) {
            multiplicityClass = fhV0MMultPercentile->GetBinContent(fhV0MMultPercentile->FindFixBin(v0am + v0cm));
            multiplicity = v0am + v0cm;
          }
          break;
        case kCL1:
          if (fhCL1MultPercentile != nullptr) {
            multiplicityClass = fhCL1MultPercentile->GetBinContent(fhCL1MultPercentile->FindFixBin(cl1m));
            multiplicity = cl1m;
          }
          break;
        case kCL1GAP:
          if (fhCL1EtaGapMultPercentile != nullptr) {
            multiplicityClass = fhCL1EtaGapMultPercentile->GetBinContent(fhCL1EtaGapMultPercentile->FindFixBin(cl1EtaGapM));
            multiplicity = cl1EtaGapM;
          }
          break;
        default:
          break;
      }
      fhMultiplicity->Fill(multiplicityClass);
    }
  }
};

//////////////////////////////////////////////////////////////////////////////
// The filter class
//////////////////////////////////////////////////////////////////////////////

struct DptDptFilter {
  struct : ConfigurableGroup {
    Configurable<std::string> cfgCCDBUrl{"input_ccdburl", "http://ccdb-test.cern.ch:8080", "The CCDB url for the input file"};
    Configurable<std::string> cfgCCDBPathName{"input_ccdbpath", "", "The CCDB path for the input file. Default \"\", i.e. don't load from CCDB"};
    Configurable<std::string> cfgCCDBDate{"input_ccdbdate", "20220307", "The CCDB date for the input file"};
    Configurable<std::string> cfgCCDBPeriod{"input_ccdbperiod", "LHC22o", "The CCDB dataset period for the input file"};
  } cfginputfile;
  Configurable<bool> cfgFullDerivedData{"fullderiveddata", false, "Produce the full derived data for external storage. Default false"};
  Configurable<std::string> cfgCentMultEstimator{"centmultestimator", "V0M", "Centrality/multiplicity estimator detector: V0M,CL0,CL1,FV0A,FT0M,FT0A,FT0C,NTPV,NOCM: none. Default V0M"};

  struct : ConfigurableGroup {
    std::string prefix = "cfgEventSelection";
    Configurable<std::string> itsDeadMaps{"itsDeadMaps", "", "Level of inactive chips: nocheck(empty), goodIts3, goodIts0123, goodItsAll. Default empty"};
    struct : ConfigurableGroup {
      std::string prefix = "cfgOccupancySelection";
      Configurable<std::string> cfgOccupancyEstimation{"cfgOccupancyEstimation", "None", "Occupancy estimation: None, Tracks, FT0C. Default None"};
      Configurable<float> cfgMinOccupancy{"cfgMinOccupancy", 0.0f, "Minimum allowed occupancy. Depends on the occupancy estimation"};
      Configurable<float> cfgMaxOccupancy{"cfgMaxOccupancy", 1e6f, "Maximum allowed occupancy. Depends on the occupancy estimation"};
    } cfgOccupancySelection;
  } cfgEventSelection;
  Configurable<std::string> cfgSystem{"syst", "PbPb", "System: pp, PbPb, Pbp, pPb, XeXe, ppRun3, PbPbRun3. Default PbPb"};
  Configurable<std::string> cfgDataType{"datatype", "data", "Data type: data, datanoevsel, MC, FastMC, OnTheFlyMC. Default data"};
  Configurable<std::string> cfgTriggSel{"triggsel", "MB", "Trigger selection: MB,VTXTOFMATCHED,VTXTRDMATCHED,VTXTRDTOFMATCHED,None. Default MB"};
  Configurable<std::string> cfgCentSpec{"centralities", "00-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80", "Centrality/multiplicity ranges in min-max separated by commas"};
  Configurable<float> cfgOverallMinP{"overallminp", 0.0f, "The overall minimum momentum for the analysis. Default: 0.0"};
  struct : ConfigurableGroup {
    std::string prefix = "cfgTpcExclusion";
    Configurable<int> method{"method", 0, "The method for excluding tracks within the TPC. 0: no exclusion; 1: static; 2: dynamic. Default: 0"};
    Configurable<std::string> positiveLowCut{"positiveLowCut", "0.0787/x - 0.0236", "The lower cut function for positive tracks"};
    Configurable<std::string> positiveUpCut{"positiveUpCut", "0.0892/x + 0.0251", "The upper cut function for positive tracks"};
    Configurable<std::string> negativeLowCut{"negativeLowCut", "pi/9.0 - (0.0892/x + 0.0251)", "The lower cut function for negative tracks"};
    Configurable<std::string> negativeUpCut{"negativeUpCut", "pi/9 - (0.0787/x - 0.0236)", "The upper cut function for negative tracks"};
  } cfgTpcExclusion;
  Configurable<o2::analysis::DptDptBinningCuts> cfgBinning{"binning",
                                                           {28, -7.0, 7.0, 18, 0.2, 2.0, 16, -0.8, 0.8, 72, 0.5},
                                                           "triplets - nbins, min, max - for z_vtx, pT, eta and phi, binning plus bin fraction of phi origin shift"};
  Configurable<bool> cfgTraceCollId0{"tracecollid0", false, "Trace particles in collisions id 0. Default false"};

  OutputObj<TList> fOutput{"DptDptFilterCollisionsInfo", OutputObjHandlingPolicy::AnalysisObject};

  Produces<aod::DptDptCFAcceptedCollisions> acceptedcollisions;
  Produces<aod::DptDptCFCollisionsInfo> collisionsinfo;
  Produces<aod::DptDptCFAcceptedTrueCollisions> acceptedtrueevents;
  Produces<aod::DptDptCFGenCollisionsInfo> gencollisionsinfo;

  Multiplicity multiplicity;

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
    if (doprocessWithoutCent || doprocessWithoutCentDetectorLevel || doprocessWithoutCentGeneratorLevel) {
      fCentMultEstimator = kNOCM;
    } else if (doprocessOnTheFlyGeneratorLevel) {
      fCentMultEstimator = kFV0A;
    } else {
      fCentMultEstimator = getCentMultEstimator(cfgCentMultEstimator);
    }
    /* the occupancy selection */
    fOccupancyEstimation = getOccupancyEstimator(cfgEventSelection.cfgOccupancySelection.cfgOccupancyEstimation);
    fMinOccupancy = cfgEventSelection.cfgOccupancySelection.cfgMinOccupancy;
    fMaxOccupancy = cfgEventSelection.cfgOccupancySelection.cfgMaxOccupancy;
    /* the ITS dead map check */
    fItsDeadMapCheck = getItsDeadMapCheck(cfgEventSelection.itsDeadMaps);

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
      fhEventSelection = new TH1D("EventSelection", ";;counts", knCollisionSelectionFlags, -0.5f, static_cast<float>(knCollisionSelectionFlags) - 0.5f);
      for (int ix = 0; ix < knCollisionSelectionFlags; ++ix) {
        fhEventSelection->GetXaxis()->SetBinLabel(ix + 1, eventSelectionSteps[ix]);
      }
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
      fOutputList->Add(fhEventSelection);
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
      if (!doprocessOnTheFlyGeneratorLevel) {
        fhTrueVertexZAA = new TH1F("TrueVertexZAA", "Vertex Z (truth rec associated); z_{vtx}", zvtxbins, zvtxlow, zvtxup);
      }

      /* add the hstograms to the output list */
      fOutputList->Add(fhTrueCentMultB);
      fOutputList->Add(fhTrueCentMultA);
      fOutputList->Add(fhTrueVertexZB);
      fOutputList->Add(fhTrueVertexZA);
      if (doprocessOnTheFlyGeneratorLevel) {
        multiplicity.init(fOutputList);
      } else {
        fOutputList->Add(fhTrueVertexZAA);
      }
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

  void processWithCentDetectorLevel(aod::CollisionEvSelCent const& collision, DptDptFullTracksDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(DptDptFilter, processWithCentDetectorLevel, "Process MC detector level with centrality", false);

  void processWithRun2CentDetectorLevel(aod::CollisionEvSelRun2Cent const& collision, DptDptFullTracksDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(DptDptFilter, processWithRun2CentDetectorLevel, "Process MC detector level with centrality", false);

  void processWithoutCentDetectorLevel(aod::CollisionEvSel const& collision, DptDptFullTracksDetLevel const& ftracks, aod::McParticles const&);
  PROCESS_SWITCH(DptDptFilter, processWithoutCentDetectorLevel, "Process MC detector level without centrality", false);

  template <typename CollisionObject, typename ParticlesList>
  bool processGenerated(CollisionObject const& mccollision, ParticlesList const& mcparticles, float centormult);

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

  void processOnTheFlyGeneratorLevel(aod::McCollision const& mccollision,
                                     aod::McParticles const& mcparticles);
  PROCESS_SWITCH(DptDptFilter, processOnTheFlyGeneratorLevel, "Process on the fly generated events", false);

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
  if (isEventSelected(collision, centormult)) {
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
  /* report the event selection */
  for (int iflag = 0; iflag < knCollisionSelectionFlags; ++iflag) {
    if (collisionFlags.test(iflag)) {
      fhEventSelection->Fill(iflag);
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

template <typename CollisionObject, typename ParticlesList>
bool DptDptFilter::processGenerated(CollisionObject const& mccollision, ParticlesList const&, float centormult)
{
  using namespace dptdptfilter;

  uint8_t acceptedevent = uint8_t(false);
  if (isEventSelected(mccollision, centormult)) {
    acceptedevent = uint8_t(true);
  }
  if (fullDerivedData) {
    acceptedtrueevents(mccollision.bcId(), mccollision.posZ(), acceptedevent, centormult);
  } else {
    gencollisionsinfo(acceptedevent, centormult);
  }
  return static_cast<bool>(acceptedevent);
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
  for (auto const& tmpcollision : collisions) {
    if (tmpcollision.has_mcCollision()) {
      if (tmpcollision.mcCollisionId() == mccollision.globalIndex()) {
        typename AllCollisions::iterator const& collision = allcollisions.iteratorAt(tmpcollision.globalIndex());
        if (isEventSelected(collision, defaultcent)) {
          if (processGenerated(mccollision, mcparticles, defaultcent)) {
            fhTrueVertexZAA->Fill((mccollision.posZ()));
          }
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

void DptDptFilter::processOnTheFlyGeneratorLevel(aod::McCollision const& mccollision,
                                                 aod::McParticles const& mcparticles)
{
  uint8_t acceptedEvent = uint8_t(false);
  fhTrueVertexZB->Fill(mccollision.posZ());
  /* we assign a default value for the time being */
  float centormult = 50.0f;
  if (isEventSelected(mccollision, centormult)) {
    acceptedEvent = true;
    multiplicity.extractMultiplicity(mcparticles);
    fhTrueVertexZA->Fill((mccollision.posZ()));
    centormult = multiplicity.getMultiplicityClass();
  }
  if (fullDerivedData) {
    acceptedtrueevents(mccollision.bcId(), mccollision.posZ(), acceptedEvent, centormult);
  } else {
    gencollisionsinfo(acceptedEvent, centormult);
  }
}

void DptDptFilter::processVertexGenerated(aod::McCollisions const& mccollisions)
{
  for (aod::McCollision const& mccollision : mccollisions) {
    fhTrueVertexZB->Fill(mccollision.posZ());
    /* we assign a default value */
    float centmult = 50.0f;
    if (isEventSelected(mccollision, centmult)) {
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
  T stdDev = std::sqrt(sqSum / vec.size());

  return stdDev;
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

  Configurable<bool> cfgOutDebugInfo{"outdebuginfo", false, "Out detailed debug information per track into a text file. Default false"};
  Configurable<bool> cfgFullDerivedData{"fullderiveddata", false, "Produce the full derived data for external storage. Default false"};
  Configurable<int> cfgTrackType{"trktype", 4, "Type of selected tracks: 0 = no selection;1 = Run2 global tracks FB96;3 = Run3 tracks;4 = Run3 tracks MM sel;5 = Run2 TPC only tracks;7 = Run 3 TPC only tracks;30-33 = any/two on 3 ITS,any/all in 7 ITS;40-43 same as 30-33 w tighter DCAxy;50-53 w tighter pT DCAz. Default 4"};
  Configurable<bool> cfgOnlyInOneSide{"onlyinoneside", false, "select tracks that don't cross the TPC central membrane. Default false"};
  Configurable<o2::analysis::CheckRangeCfg> cfgTraceDCAOutliers{"trackdcaoutliers", {false, 0.0, 0.0}, "Track the generator level DCAxy outliers: false/true, low dcaxy, up dcaxy. Default {false,0.0,0.0}"};
  Configurable<float> cfgTraceOutOfSpeciesParticles{"trackoutparticles", false, "Track the particles which are not e,mu,pi,K,p: false/true. Default false"};
  Configurable<int> cfgRecoIdMethod{"recoidmethod", 0, "Method for identifying reconstructed tracks: 0 No PID, 1 PID, 2 mcparticle, 3 mcparticle only primaries, 4 mcparticle only sec, 5 mcparicle only sec from decays, 6 mcparticle only sec from material. Default 0"};
  Configurable<o2::analysis::TrackSelectionTuneCfg> cfgTuneTrackSelection{"tunetracksel", {}, "Track selection: {useit: true/false, tpccls-useit, tpcxrws-useit, tpcxrfc-useit, tpcshcls-useit, dcaxy-useit, dcaz-useit}. Default {false,0.70,false,0.8,false,0.4,false,2.4,false,3.2,false}"};
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
  Service<o2::framework::O2DatabasePDG> fPDG;
  PIDSpeciesSelection pidselector;
  bool checkAmbiguousTracks = false;

  std::vector<bool> particleReconstructed;

  void init(InitContext& initContext)
  {
    LOGF(info, "DptDptFilterTracks::init()");

    fullDerivedData = cfgFullDerivedData;

    /* update with the configurable values */
    /* self configure the binning */
    getTaskOptionValue(initContext, "dpt-dpt-filter", "overallminp", overallminp, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mZVtxbins", zvtxbins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mZVtxmin", zvtxlow, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mZVtxmax", zvtxup, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPTbins", ptbins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPTmin", ptlow, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPTmax", ptup, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mEtabins", etabins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mEtamin", etalow, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mEtamax", etaup, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPhibins", phibins, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPhibinshift", phibinshift, false);

    TpcExclusionMethod tpcExclude = kNOEXCLUSION; ///< exclude tracks within the TPC according to this method
    std::string pLowCut;
    std::string pUpCut;
    std::string nLowCut;
    std::string nUpCut;
    {
      int tmpTpcExclude = 0;
      getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgTpcExclusion.method", tmpTpcExclude, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgTpcExclusion.positiveLowCut", pLowCut, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgTpcExclusion.positiveUpCut", pUpCut, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgTpcExclusion.negativeLowCut", nLowCut, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "cfgTpcExclusion.negativeUpCut", nUpCut, false);
      tpcExclude = static_cast<TpcExclusionMethod>(tmpTpcExclude);
    }
    /* self configure the CCDB access to the input file */
    getTaskOptionValue(initContext, "dpt-dpt-filter", "input_ccdburl", cfgCCDBUrl, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "input_ccdbpath", cfgCCDBPathName, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "input_ccdbdate", cfgCCDBDate, false);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "input_ccdbperiod", cfgCCDBPeriod, false);

    /* the track types and combinations */
    tracktype = cfgTrackType.value;
    initializeTrackSelection(cfgTuneTrackSelection.value);
    traceDCAOutliers = cfgTraceDCAOutliers;
    traceOutOfSpeciesParticles = cfgTraceOutOfSpeciesParticles;
    recoIdMethod = cfgRecoIdMethod;
    onlyInOneSide = cfgOnlyInOneSide.value;

    /* the TPC excluder object instance */
    tpcExcluder = TpcExcludeTrack(tpcExclude);
    tpcExcluder.setCuts(pLowCut, pUpCut, nLowCut, nUpCut);

    /* self configure system type and data type */
    /* if the system type is not known at this time, we have to put the initialization somewhere else */
    std::string tmpstr;
    getTaskOptionValue(initContext, "dpt-dpt-filter", "syst", tmpstr, false);
    fSystem = getSystemType(tmpstr);
    getTaskOptionValue(initContext, "dpt-dpt-filter", "datatype", tmpstr, false);
    fDataType = getDataType(tmpstr);

    /* required ambiguous tracks checks? */
    if (dofilterDetectorLevelWithoutPIDAmbiguous || dofilterDetectorLevelWithPIDAmbiguous || dofilterDetectorLevelWithFullPIDAmbiguous ||
        dofilterRecoWithoutPIDAmbiguous || dofilterRecoWithPIDAmbiguous || dofilterRecoWithFullPIDAmbiguous) {
      checkAmbiguousTracks = true;
    }

    /* configure the PID selection */
    auto insertInPIDselector = [&](auto cfg, uint sp) {
      if (cfg.value.mUseIt) {
        if (cfg.value.mExclude) {
          pidselector.addExcludedSpecies(sp, &(cfg.value));
          LOGF(info, "Incorporated species: %s to PID selection for exclusion", pidselector.spnames[sp].data());
        } else {
          pidselector.addSpecies(sp, &(cfg.value));
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
    fOutputList->Add(new TParameter<int>("TrackType", cfgTrackType, 'f'));
    fOutputList->Add(new TParameter<int>("TrackOneCharge", 1, 'f'));
    fOutputList->Add(new TParameter<int>("TrackTwoCharge", -1, 'f'));

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
    /* the debug info output file if required */
    if (cfgOutDebugInfo) {
      debugstream.open("tracings.csv");
      debugstream << "p,piw,pt,hastof,dEdx,beta,tpcnEl,tpcnMu,tpcnPi,tpcnKa,tpcnPr,tofnEl,tofnMu,tofnPi,tofnKa,tofnPr,tpcnElSft,tpcnMuSft,tpcnPiSft,tpcnKaSft,tpcnPrSft,tofnElSft,tofnMuSft,tofnPiSft,tofnKaSft,tofnPrSft,idcode,pid,pid2,truepid,phprim,process\n";
      debugstream.close();
    }
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

  template <StrongDebugging outdebug, typename TrackObject>
  int8_t trackIdentification(TrackObject const& track);
  template <StrongDebugging outdebug, typename CollisionsObject, typename TrackObject>
  int8_t selectTrack(TrackObject const& track);
  template <StrongDebugging outdebug, typename CollisionObjects, typename TrackObject>
  int8_t selectTrackAmbiguousCheck(CollisionObjects const& collisions, TrackObject const& track);
  template <typename ParticleObject>
  int8_t identifyParticle(ParticleObject const& particle);
  template <typename ParticleObject>
  int8_t identifyPrimaryParticle(ParticleObject const& particle);
  template <typename ParticleObject>
  int8_t identifySecondaryParticle(ParticleObject const& particle);
  template <typename ParticleObject>
  int8_t identifySecFromDecayParticle(ParticleObject const& particle);
  template <typename ParticleObject>
  int8_t identifySecFromMaterialParticle(ParticleObject const& particle);
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
  template <StrongDebugging outdebug, typename passedtracks>
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
    for (auto const& collision : collisions) {
      if (collision.collisionaccepted()) {
        ncollaccepted++;
      }
    }
    for (auto const& track : tracks) {
      int8_t pid = -1;
      if (track.has_collision() && (track.template collision_as<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>()).collisionaccepted()) {
        pid = selectTrackAmbiguousCheck<outdebug>(collisions, track);
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

  /* filter the tracks but not creating the filtered tracks table */
  /* the aim is to fill the structure of the generated particles  */
  /* that were reconstructed                                      */
  template <typename passedtracks>
  void filterTracksSpecial(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const&, passedtracks const& tracks)
  {
    /* do check for special adjustments */
    getCCDBInformation();

    for (auto const& track : tracks) {
      int8_t pid = -1;
      if (track.has_collision() && (track.template collision_as<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>()).collisionaccepted()) {
        pid = selectTrack<kNODEBUG, soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>(track);
        if (!(pid < 0)) {
          particleReconstructed[track.mcParticleId()] = true;
        }
      }
    }
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

    for (auto const& gencoll : gencollisions) {
      if (gencoll.collisionaccepted()) {
        acceptedcollisions++;
      }
    }

    for (auto const& particle : particles) {
      int8_t pid = -1;
      auto pdgpart = fPDG->GetParticle(particle.pdgCode());
      float charge = pdgpart != nullptr ? getCharge(pdgpart->Charge()) : 0;

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

  /* we produce the derived particle table incoporating only the particles that were accepted but not were reconstructed */
  void filterParticlesSpecial(soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo> const& gencollisions, aod::McParticles const& particles)
  {
    using namespace dptdptfilter;

    int acceptedparticles = 0;
    int acceptedcollisions = 0;
    if (!fullDerivedData) {
      gentracksinfo.reserve(particles.size());
    }

    for (auto const& gencoll : gencollisions) {
      if (gencoll.collisionaccepted()) {
        acceptedcollisions++;
      }
    }

    for (auto const& particle : particles) {
      int8_t pid = -1;
      auto pdgpart = fPDG->GetParticle(particle.pdgCode());
      float charge = pdgpart != nullptr ? getCharge(pdgpart->Charge()) : 0;

      if (charge != 0) {
        if (particle.has_mcCollision() && (particle.template mcCollision_as<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>()).collisionaccepted()) {
          auto mccollision = particle.template mcCollision_as<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>();
          pid = selectParticle(particle, mccollision);
          if (!(pid < 0)) {
            if (particleReconstructed[particle.globalIndex()]) {
              /* the particle was reconstructed and accepted, reject it */
              pid = -1;
            } else {
              acceptedparticles++;
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

  template <typename passedtracks>
  void doFilterTracks(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const& collisions, passedtracks const& tracks)
  {
    if (cfgOutDebugInfo) {
      debugstream.open("tracings.csv", std::ios::app);
      filterTracks<kDEBUG>(collisions, tracks);
      debugstream.close();
    } else {
      filterTracks<kNODEBUG>(collisions, tracks);
    }
  }

  void filterRecoWithPID(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksPID const& tracks)
  {
    doFilterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterRecoWithPID, "Not stored derived data track filtering", false)

  void filterRecoWithFullPID(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksFullPID const& tracks)
  {
    doFilterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterRecoWithFullPID, "Not stored derived data track filtering", false)

  void filterRecoWithPIDAmbiguous(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksPIDAmbiguous const& tracks)
  {
    doFilterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterRecoWithPIDAmbiguous, "Not stored derived data track filtering with ambiguous tracks check", false)

  void filterRecoWithFullPIDAmbiguous(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksFullPIDAmbiguous const& tracks)
  {
    doFilterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterRecoWithFullPIDAmbiguous, "Not stored derived data track filtering with ambiguous tracks check", false)

  void filterDetectorLevelWithPID(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksPIDDetLevel const& tracks, aod::McParticles const&)
  {
    doFilterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterDetectorLevelWithPID, "Not stored derived data detector level track filtering", false)

  void filterDetectorLevelWithFullPID(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksFullPIDDetLevel const& tracks, aod::McParticles const&)
  {
    doFilterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterDetectorLevelWithFullPID, "Not stored derived data detector level track filtering", false)

  void filterDetectorLevelWithPIDAmbiguous(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksPIDDetLevelAmbiguous const& tracks, aod::McParticles const&)
  {
    doFilterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterDetectorLevelWithPIDAmbiguous, "Not stored derived data detector level track filtering with ambiguous tracks check", false)

  void filterDetectorLevelWithFullPIDAmbiguous(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksFullPIDDetLevelAmbiguous const& tracks, aod::McParticles const&)
  {
    doFilterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterDetectorLevelWithFullPIDAmbiguous, "Not stored derived data detector level track filtering with ambiguous tracks check", false)

  void filterRecoWithoutPID(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const& collisions, DptDptFullTracks const& tracks)
  {
    doFilterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterRecoWithoutPID, "Track filtering without PID information", true)

  void filterRecoWithoutPIDAmbiguous(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const& collisions, DptDptFullTracksAmbiguous const& tracks)
  {
    doFilterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterRecoWithoutPIDAmbiguous, "Track filtering without PID information with ambiguous tracks check", false)

  void filterDetectorLevelWithoutPID(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const& collisions, DptDptFullTracksDetLevel const& tracks, aod::McParticles const&)
  {
    doFilterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterDetectorLevelWithoutPID, "Detector level track filtering without PID information", false)

  void filterDetectorLevelWithoutPIDAmbiguous(soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo> const& collisions, DptDptFullTracksDetLevelAmbiguous const& tracks, aod::McParticles const&)
  {
    doFilterTracks(collisions, tracks);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterDetectorLevelWithoutPIDAmbiguous, "Detector level track filtering without PID information with ambiguous tracks check", false)

  void filterGenerated(soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo> const& gencollisions, aod::McParticles const& particles)
  {
    filterParticles(gencollisions, particles);
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterGenerated, "Generated particles filtering", true)

  void filterGeneratedNotReconstructedWithPID(soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo> const& gencollisions, aod::McParticles const& particles,
                                              soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>& collisions, DptDptFullTracksPIDDetLevel const& tracks)
  {
    particleReconstructed.resize(particles.size());
    filterTracksSpecial(collisions, tracks);
    filterParticlesSpecial(gencollisions, particles);
    particleReconstructed.clear();
  }
  PROCESS_SWITCH(DptDptFilterTracks, filterGeneratedNotReconstructedWithPID, "Generated particles filtering", false)
};

template <StrongDebugging outdebug, typename TrackObject>
int8_t DptDptFilterTracks::trackIdentification(TrackObject const& track)
{
  using namespace dptdptfilter;

  int8_t sp = -127;
  if (recoIdMethod == 0) {
    sp = 0;
  } else if (recoIdMethod == 1) {
    if constexpr (framework::has_type_v<aod::pidtpc_tiny::TPCNSigmaStorePi, typename TrackObject::all_columns> || framework::has_type_v<aod::pidtpc::TPCNSigmaPi, typename TrackObject::all_columns>) {
      sp = pidselector.whichSpecies<outdebug>(track);
    } else {
      LOGF(fatal, "Track identification required but PID information not present");
    }
  } else {
    if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
      switch (recoIdMethod) {
        case 2:
          sp = identifyParticle(track.template mcParticle_as<aod::McParticles>());
          break;
        case 3:
          sp = identifyPrimaryParticle(track.template mcParticle_as<aod::McParticles>());
          break;
        case 4:
          sp = identifySecondaryParticle(track.template mcParticle_as<aod::McParticles>());
          break;
        case 5:
          sp = identifySecFromDecayParticle(track.template mcParticle_as<aod::McParticles>());
          break;
        case 6:
          sp = identifySecFromMaterialParticle(track.template mcParticle_as<aod::McParticles>());
          break;
        default:
          LOGF(fatal, "Track identification method %d not recognized. Fix it!", recoIdMethod);
      }
    } else {
      LOGF(fatal, "Track identification required from MC particle but MC information not present");
    }
  }
  return sp;
}

template <StrongDebugging outdebug, typename CollisionsObject, typename TrackObject>
int8_t DptDptFilterTracks::selectTrack(TrackObject const& track)
{
  using namespace dptdptfilter;

  /* before track selection */
  fillTrackHistosBeforeSelection(track);

  /* track selection */
  int8_t sp = -127;
  if (acceptTrack<CollisionsObject>(track)) {
    /* the track has been accepted */
    /* let's identify it */
    sp = trackIdentification<outdebug>(track);
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

template <StrongDebugging outdebug, typename CollisionObjects, typename TrackObject>
int8_t DptDptFilterTracks::selectTrackAmbiguousCheck(CollisionObjects const& collisions, TrackObject const& track)
{
  bool ambiguoustrack = false;
  int ambtracktype = 0; /* no ambiguous */
  std::vector<double> zvertexes{};
  /* ambiguous tracks checks if required */
  if constexpr (has_type_v<aod::track_association::CollisionIds, typename TrackObject::all_columns>) {
    if (track.compatibleCollIds().size() > 0) {
      if (track.compatibleCollIds().size() == 1) {
        if (track.collisionId() != track.compatibleCollIds()[0]) {
          /* ambiguous track! */
          ambiguoustrack = true;
          /* in principle we should not be here because the track is associated to two collisions at least */
          ambtracktype = 2;
          zvertexes.push_back(collisions.iteratorAt(track.collisionId()).posZ());
          zvertexes.push_back(collisions.iteratorAt(track.compatibleCollIds()[0]).posZ());
        } else {
          /* we consider the track as no ambiguous */
          ambtracktype = 1;
        }
      } else {
        /* ambiguous track! */
        ambiguoustrack = true;
        ambtracktype = 3;
        /* the track is associated to more than one collision */
        for (const auto& collIdx : track.compatibleCollIds()) {
          zvertexes.push_back(collisions.iteratorAt(collIdx).posZ());
        }
      }
    }
  }

  float multiplicityClass = (track.template collision_as<CollisionObjects>()).centmult();
  if (ambiguoustrack) {
    /* keep track of ambiguous tracks */
    fhAmbiguousTrackType->Fill(ambtracktype, multiplicityClass);
    fhAmbiguousTrackPt->Fill(track.pt(), multiplicityClass);
    fhAmbiguityDegree->Fill(zvertexes.size(), multiplicityClass);
    if (ambtracktype == 2) {
      fhCompatibleCollisionsZVtxRms->Fill(-computeRMS(zvertexes), multiplicityClass);
    } else {
      fhCompatibleCollisionsZVtxRms->Fill(computeRMS(zvertexes), multiplicityClass);
    }
    return -1;
  } else {
    if (checkAmbiguousTracks) {
      /* feedback of no ambiguous tracks only if checks required */
      fhAmbiguousTrackType->Fill(ambtracktype, multiplicityClass);
    }
    return selectTrack<outdebug, CollisionObjects>(track);
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

template <typename ParticleObject>
inline int8_t DptDptFilterTracks::identifyPrimaryParticle(ParticleObject const& particle)
{
  using namespace dptdptfilter;
  return pidselector.whichTruthPrimarySpecies(particle);
}

template <typename ParticleObject>
inline int8_t DptDptFilterTracks::identifySecondaryParticle(ParticleObject const& particle)
{
  using namespace dptdptfilter;
  return pidselector.whichTruthSecondarySpecies(particle);
}

template <typename ParticleObject>
inline int8_t DptDptFilterTracks::identifySecFromDecayParticle(ParticleObject const& particle)
{
  using namespace dptdptfilter;
  return pidselector.whichTruthSecFromDecaySpecies(particle);
}

template <typename ParticleObject>
inline int8_t DptDptFilterTracks::identifySecFromMaterialParticle(ParticleObject const& particle)
{
  using namespace dptdptfilter;
  return pidselector.whichTruthSecFromMaterialSpecies(particle);
}

template <typename ParticleObject, typename MCCollisionObject>
inline int8_t DptDptFilterTracks::selectParticle(ParticleObject const& particle, MCCollisionObject const& mccollision)
{
  int8_t sp = -127;
  auto pdgpart = fPDG->GetParticle(particle.pdgCode());
  float charge = pdgpart != nullptr ? getCharge(pdgpart->Charge()) : 0;
  if (charge != 0) {
    /* before particle selection */
    fillParticleHistosBeforeSelection(particle, mccollision, charge);

    /* track selection */
    if (acceptParticle(particle, mccollision)) {
      /* the particle has been accepted */
      /* the particle is only accepted if it is a primary particle */
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

  float dcaxy = std::sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                          (particle.vy() - collision.posY()) * (particle.vy() - collision.posY()));
  if (traceDCAOutliers.mDoIt && (traceDCAOutliers.mLowValue < dcaxy) && (dcaxy < traceDCAOutliers.mUpValue)) {
    fhTrueDCAxyBid->Fill(TString::Format("%d", particle.pdgCode()).Data(), 1.0);
  }

  fhTrueDCAxyB->Fill(std::sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                               (particle.vy() - collision.posY()) * (particle.vy() - collision.posY())));
  fhTrueDCAzB->Fill((particle.vz() - collision.posZ()));
}

template <typename ParticleObject, typename MCCollisionObject>
void DptDptFilterTracks::fillParticleHistosAfterSelection(ParticleObject const& particle, MCCollisionObject const& collision, float charge, int8_t sp)
{
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
