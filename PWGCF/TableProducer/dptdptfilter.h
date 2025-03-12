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

/// \file dptdptfilter.h
/// \brief Filters collisions and tracks according to selection criteria
/// \author victor.gonzalez.sebastian@gmail.com

#ifndef PWGCF_TABLEPRODUCER_DPTDPTFILTER_H_
#define PWGCF_TABLEPRODUCER_DPTDPTFILTER_H_

#include <CCDB/BasicCCDBManager.h>
#include <TF1.h>
#include <TList.h>
#include <vector>
#include <bitset>
#include <string>
#include <iomanip>
#include <fstream>
#include <locale>
#include <sstream>
#include <ranges>
#include <functional>
#include <map>

#include "ReconstructionDataFormats/PID.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "PWGCF/Core/AnalysisConfigurableCuts.h"

namespace o2
{
namespace aod
{
using CollisionsEvSelCent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>;
using CollisionEvSelCent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>::iterator;
using CollisionsEvSelRun2Cent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2CL0s, aod::CentRun2CL1s>;
using CollisionEvSelRun2Cent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2CL0s, aod::CentRun2CL1s>::iterator;
using CollisionsEvSel = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
using CollisionEvSel = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>::iterator;
using TrackData = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>::iterator;
} // namespace aod
namespace analysis
{
namespace dptdptfilter
{
/// \enum SystemType
/// \brief The type of the system under analysis
enum SystemType {
  kNoSystem = 0, ///< no system defined
  kpp,           ///< **p-p** system
  kpPb,          ///< **p-Pb** system
  kPbp,          ///< **Pb-p** system
  kPbPb,         ///< **Pb-Pb** system
  kXeXe,         ///< **Xe-Xe** system
  kppRun3,       ///< **p-p Run 3** system
  kPbPbRun3,     ///< **Pb-Pb Run 3** system
  knSystems      ///< number of handled systems
};

/// \enum DataType
/// \brief Which kind of data is the task addressing
enum DataType {
  kData = 0,     ///< actual data, not generated
  kMC,           ///< Generator level and detector level
  kFastMC,       ///< Gererator level but stored dataset
  kOnTheFly,     ///< On the fly generator level data
  kDataNoEvtSel, ///< actual data but not event selection available yet
  knGenData      ///< number of different generator data types
};

/// \enum CentMultEstimatorType
/// \brief The detector used to estimate centrality/multiplicity
enum CentMultEstimatorType {
  kNOCM = 0,           ///< do not use centrality/multiplicity estimator
  kV0M,                ///< V0M centrality/multiplicity estimator Run 1/2
  kCL0,                ///< CL0 centrality/multiplicity estimator Run 1/2
  kCL1,                ///< CL1 centrality/multiplicity estimator Run 1/2
  kFV0A,               ///< FV0A centrality/multiplicity estimator Run 3
  kFT0M,               ///< FT0M centrality/multiplicity estimator Run 3
  kFT0A,               ///< FT0A centrality/multiplicity estimator Run 3
  kFT0C,               ///< FT0C centrality/multiplicity estimator Run 3
  kNTPV,               ///< NTPV centrality/multiplicity estimator Run 3
  knCentMultEstimators ///< number of centrality/mutiplicity estimator
};

/// \enum TriggerSelectionType
/// \brief The type of trigger to apply for event selection
enum TriggerSelectionType {
  kNONE = 0,              ///< do not use trigger selection
  kMB,                    ///< Minimum bias trigger
  kMBEXTRA,               ///< Additional Run3 event quality
  kVTXTOFMATCHED,         ///< at least one vertex contributor is matched to TOF
  kVTXTRDMATCHED,         ///< at least one vertex contributor is matched to TRD
  kVTXTRDTOFMATCHED,      ///< at least one vertex contributor is matched to TRD and TOF
  kEXTRAVTXTOFMATCHED,    ///< Additional Run3 event quality and at least one vertex contributor is matched to TOF
  kEXTRAVTXTRDMATCHED,    ///< Additional Run3 event quality and at least one vertex contributor is matched to TRD
  kEXTRAVTXTRDTOFMATCHED, ///< Additional Run3 event quality and at least one vertex contributor is matched to TRD and TOF
  knEventSelection        ///< number of triggers for event selection
};

/// \enum OccupancyEstimationType
/// \brief The type of occupancy estimation
enum OccupancyEstimationType {
  kNOOCC = 0,           ///< do not use occupancy estimation
  kTRACKSOCC,           ///< occupancy estimated using tracks
  kFT0COCC,             ///< occupancy estimated using the FT0C
  knOccupancyEstimators ///< the number of occupancy estimators
};

/// \enum CollisionSelectionFlags
/// \brief The different criteria for selecting/rejecting collisions
enum CollisionSelectionFlags {
  kIN = 0,                  ///< new unhandled, yet, event
  kMBBIT,                   ///< minimum  bias
  kINT7BIT,                 ///< INT7 Run 1/2
  kSEL7BIT,                 ///< Sel7 Run 1/2
  kSEL8BIT,                 ///< Sel8
  kNOSAMEBUNCHPUPBIT,       ///< no same bunch pile up
  kISGOODZVTXFT0VSPVBIT,    ///< good zvtx FT0 vs PV
  kISVERTEXITSTPCBIT,       ///< is vertex TPC and ITS
  kISVERTEXTOFMATCHEDBIT,   ///< vertex contributor with TOF matched
  kISVERTEXTRDMATCHEDBIT,   ///< vertex contributor with TRD matche
  kNOCOLLINTIMERANGEBIT,    ///< no collision in time range
  kNOCOLLINROFBIT,          ///< no collision in readout
  kOCCUPANCYBIT,            ///< occupancy within limits
  kISGOODITSLAYER3BIT,      ///< right level of inactive chips for ITS layer 3
  kISGOODITSLAYER0123BIT,   ///< right level of inactive chips for ITS layers 0,1,2, and 3
  kISGOODITSLAYERALLBIT,    ///< right level of inactive chips for all seven ITS layers
  kCENTRALITYBIT,           ///< centrality cut passed
  kZVERTEXBIT,              ///< zvtx cut passed
  kSELECTED,                ///< the event has passed all selections
  knCollisionSelectionFlags ///< number of flags
};

/// \enum StrongDebugging
/// \brief Enable a per track information debugging. Only for local analyses
enum StrongDebugging {
  kNODEBUG = 0, ///< do not debug
  kDEBUG        ///< output debugging information on a per track basis to a text file
};

/// \enum TpcExclusionMethod
/// \brief Methods for excluding tracks witin the TPC
enum TpcExclusionMethod {
  kNOEXCLUSION = 0, ///< do not exclude tracks within the TPC
  kSTATIC,          ///< exclude tracks statically on the bins of the TPC sector borders; only valid if 72 bins and origin shifted by 0.5
  kDYNAMIC          ///< pT dependent exclusion matching the sector borders a la Alex Dobrin
};

/// \enum ItsDeadMapsCheckType
/// \brief Check for the right level of ITS dead chips
enum ItsDeadMapsCheckType {
  kNOCHECK = 0,        ///< no check
  kGOODITSLAYER3,      ///< check good the 3 ITS layer
  kGOODITSLAYER0123,   ///< check good the 0,1,2,and 3 ITS layers
  kGOODITSLAYERALL,    ///< check good all ITS layers
  kNOGOODITSLAYER3,    ///< check no good the  3 ITS layer
  kNOGOODITSLAYER0123, ///< check no good the 0,1,2,and 3 ITS layers
  kNOGOODITSLAYERALL   ///< check no good all ITS layers
};

//============================================================================================
// The debug output stream
//============================================================================================
std::ofstream debugstream;

//============================================================================================
// The overall minimum momentum
//============================================================================================
float overallminp = 0.0f;

//============================================================================================
// The collision selection flags and configuration objects
//============================================================================================
std::bitset<32> collisionFlags;

//============================================================================================
// The DptDptFilter configuration objects
//============================================================================================
int ptbins = 18;
float ptlow = 0.2, ptup = 2.0;
int etabins = 16;
float etalow = -0.8, etaup = 0.8;
int zvtxbins = 40;
float zvtxlow = -10.0, zvtxup = 10.0;
int phibins = 72;
float philow = 0.0f;
float phiup = constants::math::TwoPI;
float phibinshift = 0.0f;

struct TpcExcludeTrack;             ///< forward declaration of the excluder object
bool onlyInOneSide = false;         ///< select only tracks that don't cross the TPC central membrane
extern TpcExcludeTrack tpcExcluder; ///< the TPC excluder object instance

/* selection criteria from PWGMM */
// default quality criteria for tracks with ITS contribution
static constexpr o2::aod::track::TrackSelectionFlags::flagtype TrackSelectionITS =
  o2::aod::track::TrackSelectionFlags::kITSNCls | o2::aod::track::TrackSelectionFlags::kITSChi2NDF |
  o2::aod::track::TrackSelectionFlags::kITSHits;

// default quality criteria for tracks with TPC contribution
static constexpr o2::aod::track::TrackSelectionFlags::flagtype TrackSelectionTPC =
  o2::aod::track::TrackSelectionFlags::kTPCNCls |
  o2::aod::track::TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  o2::aod::track::TrackSelectionFlags::kTPCChi2NDF;

// default standard DCA cuts
static constexpr o2::aod::track::TrackSelectionFlags::flagtype TrackSelectionDCA =
  o2::aod::track::TrackSelectionFlags::kDCAz | o2::aod::track::TrackSelectionFlags::kDCAxy;

int tracktype = 1;
std::vector<TrackSelection*> trackFilters = {};              // the vector of track selectors
std::vector<std::function<float(float)>> maxDcaZPtDeps = {}; // max dca in z axis as function of pT for each track selector
TrackSelection* extraTrackFilter = nullptr;
bool dca2Dcut = false;
float sharedTpcClusters = 1.0; ///< max fraction of shared TPC clusters
float maxDCAz = 1e6f;
float maxDCAxy = 1e6f;

inline TList* getCCDBInput(auto& ccdb, const char* ccdbpath, const char* ccdbdate, const char* period = "")
{
  std::tm cfgtm = {};
  std::stringstream ss(ccdbdate);
  ss >> std::get_time(&cfgtm, "%Y%m%d");
  cfgtm.tm_hour = 12;
  int64_t timestamp = std::mktime(&cfgtm) * 1000;

  TList* lst = nullptr;
  if (std::strlen(period) != 0) {
    std::map<std::string, std::string> metadata{{"Period", period}};
    lst = ccdb->template getSpecific<TList>(ccdbpath, timestamp, metadata);
  } else {
    lst = ccdb->template getForTimeStamp<TList>(ccdbpath, timestamp);
  }
  if (lst != nullptr) {
    LOGF(info, "Correctly loaded CCDB input object");
  } else {
    LOGF(error, "CCDB input object could not be loaded");
  }
  return lst;
}

inline void initializeTrackSelection(TrackSelectionTuneCfg& tune)
{
  auto addTrackFilter = [](auto filter, std::function<float(float)> dcaZCutFunc = {}) {
    trackFilters.push_back(filter);
    maxDcaZPtDeps.push_back(dcaZCutFunc);
  };
  auto highQualityTpcTrack = []() {
    TrackSelection* tpcTrack = new TrackSelection(getGlobalTrackSelection());
    tpcTrack->ResetITSRequirements();
    tpcTrack->SetRequireITSRefit(false);
    tpcTrack->SetMinNClustersTPC(120);
    sharedTpcClusters = 0.2;
    return tpcTrack;
  };
  switch (tracktype) {
    case 1: { /* Run2 global track */
      TrackSelection* globalRun2 = new TrackSelection(getGlobalTrackSelection());
      globalRun2->SetTrackType(o2::aod::track::Run2Track); // Run 2 track asked by default
      globalRun2->SetMaxChi2PerClusterTPC(2.5f);
      TrackSelection* globalSDDRun2 = new TrackSelection(getGlobalTrackSelectionSDD());
      globalSDDRun2->SetTrackType(o2::aod::track::Run2Track); // Run 2 track asked by default
      globalSDDRun2->SetMaxChi2PerClusterTPC(2.5f);
      addTrackFilter(globalRun2);
      addTrackFilter(globalSDDRun2);
    } break;
    case 3: { /* Run3 track */
      TrackSelection* globalRun3 = new TrackSelection(getGlobalTrackSelection());
      globalRun3->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
      globalRun3->ResetITSRequirements();
      globalRun3->SetRequireHitsInITSLayers(1, {0, 1, 2});
      TrackSelection* globalSDDRun3 = new TrackSelection(getGlobalTrackSelection());
      globalSDDRun3->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
      globalSDDRun3->ResetITSRequirements();
      globalSDDRun3->SetRequireNoHitsInITSLayers({0, 1, 2});
      globalSDDRun3->SetRequireHitsInITSLayers(1, {3});
      addTrackFilter(globalRun3);
      addTrackFilter(globalSDDRun3);
    } break;
    case 5: { /* Run2 TPC only track */
      TrackSelection* tpcOnly = new TrackSelection;
      tpcOnly->SetTrackType(o2::aod::track::Run2Track); // Run 2 track asked by default
      tpcOnly->SetMinNClustersTPC(50);
      tpcOnly->SetMaxChi2PerClusterTPC(4);
      tpcOnly->SetMaxDcaZ(3.2f);
      maxDCAz = 3.2;
      tpcOnly->SetMaxDcaXY(2.4f);
      maxDCAxy = 2.4;
      dca2Dcut = true;
      addTrackFilter(tpcOnly);
    } break;
    case 7: { /* Run3 TPC only track */
      TrackSelection* tpcOnly = new TrackSelection;
      tpcOnly->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
      tpcOnly->SetMinNClustersTPC(50);
      tpcOnly->SetMaxChi2PerClusterTPC(4);
      tpcOnly->SetMaxDcaZ(3.2f);
      maxDCAz = 3.2;
      tpcOnly->SetMaxDcaXY(2.4f);
      maxDCAxy = 2.4;
      dca2Dcut = true;
      addTrackFilter(tpcOnly);
    } break;
    case 10: { /* Run3 track primary vertex contributor */
      TrackSelection* globalRun3 = new TrackSelection(getGlobalTrackSelection());
      globalRun3->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
      globalRun3->ResetITSRequirements();
      globalRun3->SetRequireHitsInITSLayers(1, {0, 1, 2});
      TrackSelection* globalSDDRun3 = new TrackSelection(getGlobalTrackSelection());
      globalSDDRun3->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
      globalSDDRun3->ResetITSRequirements();
      globalSDDRun3->SetRequireNoHitsInITSLayers({0, 1, 2});
      globalSDDRun3->SetRequireHitsInITSLayers(1, {3});
      addTrackFilter(globalRun3);
      addTrackFilter(globalSDDRun3);
    } break;
    case 30: { /* Run 3 default global track: kAny on 3 IB layers of ITS */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default)));
    } break;
    case 31: { /* Run 3 global track: kTwo on 3 IB layers of ITS */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibTwo, TrackSelection::GlobalTrackRun3DCAxyCut::Default)));
    } break;
    case 32: { /* Run 3 global track: kAny on all 7 layers of ITS */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSallAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default)));
    } break;
    case 33: { /* Run 3 global track: kAll on all 7 layers of ITS */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::Default)));
    } break;
    case 40: { /* Run 3 global track: kAny on 3 IB layers of ITS, tighter DCAxy */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)));
    } break;
    case 41: { /* Run 3 global track: kTwo on 3 IB layers of ITS, tighter DCAxy */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibTwo, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)));
    } break;
    case 42: { /* Run 3 global track: kAny on all 7 layers of ITS, tighter DCAxy */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSallAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)));
    } break;
    case 43: { /* Run 3 global track: kAll on all 7 layers of ITS, tighter DCAxy */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)));
    } break;
    case 50: { /* Run 3 global track: kAny on 3 IB layers of ITS, tighter DCAxy, tighter pT dep DCAz */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; });
    } break;
    case 51: { /* Run 3 global track: kTwo on 3 IB layers of ITS, tighter DCAxy, tighter pT dep DCAz */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibTwo, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; });
    } break;
    case 52: { /* Run 3 global track: kAny on all 7 layers of ITS, tighter DCAxy, tighter pT dep DCAz */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSallAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; });
    } break;
    case 53: { /* Run 3 global track: kAll on all 7 layers of ITS, tighter DCAxy, tighter pT dep DCAz */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; });
    } break;
    case 60: { /* Run 3 global track: kAny on 3 IB layers of ITS, tighter DCAxy, tighter pT dep DCAz, plus TPC+TOF only tracks */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; });
      extraTrackFilter = highQualityTpcTrack();
    } break;
    case 61: { /* Run 3 global track: kTwo on 3 IB layers of ITS, tighter DCAxy, tighter pT dep DCAz, plus TPC+TOF only tracks */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibTwo, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; });
      extraTrackFilter = highQualityTpcTrack();
    } break;
    case 62: { /* Run 3 global track: kAny on all 7 layers of ITS, tighter DCAxy, tighter pT dep DCAz, plus TPC+TOF only tracks */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSallAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; });
      extraTrackFilter = highQualityTpcTrack();
    } break;
    case 63: { /* Run 3 global track: kAll on all 7 layers of ITS, tighter DCAxy, tighter pT dep DCAz, plus TPC+TOF only tracks */
      addTrackFilter(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; });
      extraTrackFilter = highQualityTpcTrack();
    } break;
    default:
      break;
  }
  if (tune.mUseIt) {
    for (auto const& filter : trackFilters) {
      if (tune.mUseTPCclusters) {
        filter->SetMinNClustersTPC(tune.mTPCclusters);
      }
      if (tune.mUseTPCxRows) {
        filter->SetMinNCrossedRowsTPC(tune.mTPCxRows);
      }
      if (tune.mUseTPCXRoFClusters) {
        filter->SetMinNCrossedRowsOverFindableClustersTPC(tune.mTPCXRoFClusters);
      }
      if (tune.mUseDCAxy) {
        /* DCAxy is tricky due to how the pT dependence is implemented */
        filter->SetMaxDcaXYPtDep([&tune](float) { return tune.mDCAxy; });
        filter->SetMaxDcaXY(tune.mDCAxy);
      }
      if (tune.mUseDCAz) {
        filter->SetMaxDcaZ(tune.mDCAz);
      }
    }
    if (tune.mUseDCAz) {
      for (auto dcaZCut : maxDcaZPtDeps) { // o2-linter: disable=const-ref-in-for-loop
        dcaZCut = [&tune](float) { return tune.mDCAz; };
      }
    }
    if (tune.mUseFractionTpcSharedClusters) {
      sharedTpcClusters = tune.mFractionTpcSharedClusters;
    }
  }
}

SystemType fSystem = kNoSystem;
DataType fDataType = kData;
CentMultEstimatorType fCentMultEstimator = kV0M;
TriggerSelectionType fTriggerSelection = kMB;
OccupancyEstimationType fOccupancyEstimation = kNOOCC; /* the occupancy estimator to use */
ItsDeadMapsCheckType fItsDeadMapCheck = kNOCHECK;      /* the check of the ITS dead maps to use */

float fMinOccupancy = 0.0f; /* the minimum allowed occupancy */
float fMaxOccupancy = 1e6f; /* the maximum allowed occupancy */

/* adaptations for the pp nightly checks */
analysis::CheckRangeCfg traceDCAOutliers;
bool traceOutOfSpeciesParticles = false;
int recoIdMethod = 0;
float particleMaxDCAxy = 999.9f;
float particleMaxDCAZ = 999.9f;
bool traceCollId0 = false;

inline TriggerSelectionType getTriggerSelection(std::string const& triggstr)
{
  if (triggstr.empty() || triggstr == "MB") {
    return kMB;
  } else if (triggstr == "MBEXTRA") {
    return kMBEXTRA;
  } else if (triggstr == "VTXTOFMATCHED") {
    return kVTXTOFMATCHED;
  } else if (triggstr == "VTXTRDMATCHED") {
    return kVTXTRDMATCHED;
  } else if (triggstr == "VTXTRDTOFMATCHED") {
    return kVTXTRDTOFMATCHED;
  } else if (triggstr == "EXTRAVTXTOFMATCHED") {
    return kEXTRAVTXTOFMATCHED;
  } else if (triggstr == "EXTRAVTXTRDMATCHED") {
    return kEXTRAVTXTRDMATCHED;
  } else if (triggstr == "EXTRAVTXTRDTOFMATCHED") {
    return kEXTRAVTXTRDTOFMATCHED;
  } else if (triggstr == "None") {
    return kNONE;
  } else {
    LOGF(fatal, "Wrong trigger selection: %s", triggstr.c_str());
    return kNONE;
  }
}

inline SystemType getSystemType(std::string const& sysstr)
{
  /* we have to figure out how extract the system type */
  if (sysstr.empty() || (sysstr == "PbPb")) {
    return kPbPb;
  } else if (sysstr == "pp") {
    return kpp;
  } else if (sysstr == "pPb") {
    return kpPb;
  } else if (sysstr == "Pbp") {
    return kPbp;
  } else if (sysstr == "pPb") {
    return kpPb;
  } else if (sysstr == "XeXe") {
    return kXeXe;
  } else if (sysstr == "ppRun3") {
    return kppRun3;
  } else if (sysstr == "PbPbRun3") {
    return kPbPbRun3;
  } else {
    LOGF(fatal, "DptDptCorrelations::getSystemType(). Wrong system type: %s", sysstr.c_str());
  }
  return kPbPb;
}

/// \brief Type of data according to the configuration string
/// \param datastr The data type configuration string
/// \return Internal code for the passed kind of data string
inline DataType getDataType(std::string const& datastr)
{
  /* we have to figure out how extract the type of data*/
  if (datastr.empty() || (datastr == "data")) {
    return kData;
  } else if (datastr == "datanoevsel") {
    return kDataNoEvtSel;
  } else if (datastr == "MC") {
    return kMC;
  } else if (datastr == "FastMC") {
    return kFastMC;
  } else if (datastr == "OnTheFlyMC") {
    return kOnTheFly;
  } else {
    LOGF(fatal, "DptDptCorrelations::getDataType(). Wrong type of dat: %d", datastr.c_str());
  }
  return kData;
}

inline CentMultEstimatorType getCentMultEstimator(std::string const& datastr)
{
  if (datastr == "V0M") {
    return kV0M;
  } else if (datastr == "CL0") {
    return kCL0;
  } else if (datastr == "CL1") {
    return kCL1;
  } else if (datastr == "FV0A") {
    return kFV0A;
  } else if (datastr == "FT0M") {
    return kFT0M;
  } else if (datastr == "FT0A") {
    return kFT0A;
  } else if (datastr == "FT0C") {
    return kFT0C;
  } else if (datastr == "NTPV") {
    return kNTPV;
  } else if (datastr == "NOCM") {
    return kNOCM;
  } else {
    LOGF(fatal, "Centrality/Multiplicity estimator %s not supported yet", datastr.c_str());
  }
  return kNOCM;
}

inline std::string getCentMultEstimatorName(CentMultEstimatorType est)
{
  switch (est) {
    case kV0M:
      return "V0M";
      break;
    case kCL0:
      return "CL0";
      break;
    case kCL1:
      return "CL1";
      break;
    case kFV0A:
      return "FV0A";
      break;
    case kFT0M:
      return "FT0M";
      break;
    case kFT0A:
      return "FT0A";
      break;
    case kFT0C:
      return "FT0C";
      break;
    case kNTPV:
      return "NTPV";
      break;
    case kNOCM:
      return "NOCM";
      break;
    default:
      LOGF(fatal, "Centrality/Multiplicity estimator %d not supported yet", (int)est);
      return "WRONG";
      break;
  }
}

inline OccupancyEstimationType getOccupancyEstimator(const std::string& estimator)
{
  if (estimator == "None") {
    return kNOOCC;
  } else if (estimator == "Tracks") {
    return kTRACKSOCC;
  } else if (estimator == "FT0C") {
    return kFT0COCC;
  } else {
    LOGF(fatal, "Occupancy estimator %s not supported yet", estimator.c_str());
    return kNOOCC;
  }
}

inline ItsDeadMapsCheckType getItsDeadMapCheck(const std::string& check)
{
  if (check.length() == 0 || check == "None") {
    return kNOCHECK;
  } else if (check == "goodIts3") {
    return kGOODITSLAYER3;
  } else if (check == "goodIts0123") {
    return kGOODITSLAYER0123;
  } else if (check == "goodItsAll") {
    return kGOODITSLAYERALL;
  } else if (check == "noGoodIts3") {
    return kNOGOODITSLAYER3;
  } else if (check == "noGoodIts0123") {
    return kNOGOODITSLAYER0123;
  } else if (check == "noGoodItsAll") {
    return kNOGOODITSLAYERALL;
  } else {
    LOGF(fatal, "ITS dead map check %s not implemented", check.c_str());
    return kNOCHECK;
  }
}

//////////////////////////////////////////////////////////////////////////////////
/// Trigger selection
//////////////////////////////////////////////////////////////////////////////////

/// \brief Trigger selection for reconstructed and detector level collision tables
template <typename CollisionObject>
inline bool triggerSelectionReco(CollisionObject const& collision)
{
  bool trigsel = false;
  switch (fSystem) {
    case kpp:
    case kpPb:
    case kPbp:
    case kPbPb:
    case kXeXe:
      switch (fTriggerSelection) {
        case kMB:
          switch (fDataType) {
            case kData:
              if (collision.alias_bit(kINT7)) {
                collisionFlags.set(kINT7BIT);
                if (collision.sel7()) {
                  trigsel = true;
                  collisionFlags.set(kSEL7BIT);
                }
              }
              break;
            case kMC:
              if (collision.sel7()) {
                trigsel = true;
                collisionFlags.set(kSEL7BIT);
              }
              break;
            default:
              collisionFlags.set(kMBBIT);
              trigsel = true;
              break;
          }
          break;
        case kNONE:
          trigsel = true;
          break;
        default:
          break;
      }
      break;
    case kppRun3:
    case kPbPbRun3: {
      auto run3Accepted = [](auto const& coll) {
        return coll.sel8() &&
               coll.selection_bit(aod::evsel::kNoCollInTimeRangeStandard) &&
               coll.selection_bit(aod::evsel::kNoCollInRofStandard);
      };
      auto run3ExtraAccepted = [](auto const& coll) {
        return coll.sel8() &&
               coll.selection_bit(aod::evsel::kNoSameBunchPileup) &&
               coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) &&
               coll.selection_bit(aod::evsel::kIsVertexITSTPC);
      };
      auto setCollisionFlags = [](auto& flags, auto const& coll) {
        flags.set(kSEL8BIT, coll.sel8() != 0);
        flags.set(kNOSAMEBUNCHPUPBIT, coll.selection_bit(aod::evsel::kNoSameBunchPileup));
        flags.set(kISGOODZVTXFT0VSPVBIT, coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV));
        flags.set(kISVERTEXITSTPCBIT, coll.selection_bit(aod::evsel::kIsVertexITSTPC));
        flags.set(kISVERTEXTOFMATCHEDBIT, coll.selection_bit(aod::evsel::kIsVertexTOFmatched));
        flags.set(kISVERTEXTRDMATCHEDBIT, coll.selection_bit(aod::evsel::kIsVertexTRDmatched));
        flags.set(kNOCOLLINTIMERANGEBIT, coll.selection_bit(aod::evsel::kNoCollInTimeRangeStandard));
        flags.set(kNOCOLLINROFBIT, coll.selection_bit(aod::evsel::kNoCollInRofStandard));
        flags.set(kISGOODITSLAYER3BIT, coll.selection_bit(aod::evsel::kIsGoodITSLayer3));
        flags.set(kISGOODITSLAYER0123BIT, coll.selection_bit(aod::evsel::kIsGoodITSLayer0123));
        flags.set(kISGOODITSLAYERALLBIT, coll.selection_bit(aod::evsel::kIsGoodITSLayersAll));
      };
      switch (fTriggerSelection) {
        case kMB:
          if (run3Accepted(collision)) {
            trigsel = true;
          }
          break;
        case kMBEXTRA:
          if (run3ExtraAccepted(collision)) {
            trigsel = true;
          }
          break;
        case kVTXTOFMATCHED:
          if (run3Accepted(collision) && collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
            trigsel = true;
          }
          break;
        case kVTXTRDMATCHED:
          if (run3Accepted(collision) && collision.selection_bit(aod::evsel::kIsVertexTRDmatched)) {
            trigsel = true;
          }
          break;
        case kVTXTRDTOFMATCHED:
          if (run3Accepted(collision) && collision.selection_bit(aod::evsel::kIsVertexTRDmatched) && collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
            trigsel = true;
          }
          break;
        case kEXTRAVTXTOFMATCHED:
          if (run3ExtraAccepted(collision) && collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
            trigsel = true;
          }
          break;
        case kEXTRAVTXTRDMATCHED:
          if (run3ExtraAccepted(collision) && collision.selection_bit(aod::evsel::kIsVertexTRDmatched)) {
            trigsel = true;
          }
          break;
        case kEXTRAVTXTRDTOFMATCHED:
          if (run3ExtraAccepted(collision) && collision.selection_bit(aod::evsel::kIsVertexTRDmatched) && collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
            trigsel = true;
          }
          break;
        case kNONE:
          trigsel = true;
          break;
        default:
          break;
      }
      setCollisionFlags(collisionFlags, collision);
    } break;
    default:
      break;
  }
  return trigsel;
}

/// \brief Trigger selection by default: unknow subscribed collision table
template <typename CollisionObject>
inline bool triggerSelection(CollisionObject const&)
{
  LOGF(fatal, "Trigger selection not implemented for this kind of collisions");
  return false;
}

/// \brief Trigger selection for reconstructed collision tables without centrality/multiplicity
template <>
inline bool triggerSelection<aod::CollisionEvSel>(aod::CollisionEvSel const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for reconstructed collision tables with Run 2 centrality/multiplicity
template <>
inline bool triggerSelection<aod::CollisionEvSelRun2Cent>(aod::CollisionEvSelRun2Cent const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for reconstructed collision tables with centrality/multiplicity
template <>
inline bool triggerSelection<aod::CollisionEvSelCent>(aod::CollisionEvSelCent const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for detector level collision tables without centrality/multiplicity
template <>
inline bool triggerSelection<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for detector level collision tables with Run 2 centrality/multiplicity
template <>
inline bool triggerSelection<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for detector level collision tables with centrality/multiplicity
template <>
inline bool triggerSelection<soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for generator level collison table
template <>
inline bool triggerSelection<aod::McCollision>(aod::McCollision const&)
{
  return true;
}

//////////////////////////////////////////////////////////////////////////////////
/// Multiplicity extraction
//////////////////////////////////////////////////////////////////////////////////

/// \brief Extract the collision multiplicity from the event selection information
template <typename CollisionObject>
inline float extractMultiplicity(CollisionObject const& collision, CentMultEstimatorType est)
{
  switch (est) {
    case kV0M:
      return collision.multFV0M();
      break;
    case kCL0:
      return collision.multTracklets();
      break;
    case kCL1:
      return collision.multTracklets();
      break;
    case kFV0A:
      return collision.multFV0A();
      break;
    case kFT0M:
      return collision.multFT0M();
      break;
    case kFT0A:
      return collision.multFT0A();
      break;
    case kFT0C:
      return collision.multFT0M();
      break;
    case kNTPV:
      return collision.multNTracksPV();
      break;
    case kNOCM:
      return collision.multFT0M();
      break;
    default:
      LOGF(fatal, "Centrality/Multiplicity estimator %d not supported yet", (int)est);
      return collision.multFT0M();
      break;
  }
}

//////////////////////////////////////////////////////////////////////////////////
/// Centrality selection
//////////////////////////////////////////////////////////////////////////////////

/// \brief Centrality/multiplicity percentile
template <typename CollisionObject>
  requires(o2::aod::HasRun2Centrality<CollisionObject>)
float getCentMultPercentile(CollisionObject collision)
{
  switch (fCentMultEstimator) {
    case kV0M:
      return collision.centRun2V0M();
    case kCL0:
      return collision.centRun2CL0();
    case kCL1:
      return collision.centRun2CL1();
    default:
      return 105.0;
  }
}

template <typename CollisionObject>
  requires(o2::aod::HasCentrality<CollisionObject>)
float getCentMultPercentile(CollisionObject collision)
{
  switch (fCentMultEstimator) {
    case kFV0A:
      return collision.centFV0A();
    case kFT0M:
      return collision.centFT0M();
    case kFT0A:
      return collision.centFT0A();
    case kFT0C:
      return collision.centFT0C();
    case kNTPV:
      return collision.centNTPV();
    default:
      return 105.0;
  }
}

/// \brief Centrality selection when there is centrality/multiplicity information
template <typename CollisionObject>
inline bool centralitySelectionMult(CollisionObject collision, float& centmult)
{
  float mult = getCentMultPercentile(collision);
  if (mult < 100 && 0 < mult) {
    centmult = mult;
    collisionFlags.set(kCENTRALITYBIT);
    return true;
  }
  return false;
}

/// \brief Centrality selection when there is not centrality/multiplicity information
template <typename CollisionObject>
inline bool centralitySelectionNoMult(CollisionObject const&, float& centmult)
{
  bool centmultsel = false;
  switch (fCentMultEstimator) {
    case kNOCM:
      centmult = 50.0;
      centmultsel = true;
      collisionFlags.set(kCENTRALITYBIT);
      break;
    default:
      break;
  }
  return centmultsel;
}

/// \brief Centrality selection by default: unknown subscribed collision table
template <typename CollisionObject>
inline bool centralitySelection(CollisionObject const&, float&)
{
  LOGF(fatal, "Centrality selection not implemented for this kind of collisions");
  return false;
}

/// \brief Centrality selection for reconstructed and detector level collision tables with centrality/multiplicity information
template <>
inline bool centralitySelection<aod::CollisionEvSelCent>(aod::CollisionEvSelCent const& collision, float& centmult)
{
  return centralitySelectionMult(collision, centmult);
}

/// \brief Centrality selection for reconstructed and detector level collision tables with Run 2 centrality/multiplicity information
template <>
inline bool centralitySelection<aod::CollisionEvSelRun2Cent>(aod::CollisionEvSelRun2Cent const& collision, float& centmult)
{
  return centralitySelectionMult(collision, centmult);
}

/// \brief Centrality selection for reconstructed and detector level collision tables without centrality/multiplicity information
template <>
inline bool centralitySelection<aod::CollisionEvSel>(aod::CollisionEvSel const& collision, float& centmult)
{
  return centralitySelectionNoMult(collision, centmult);
}

/// \brief Centrality selection for detector level collision tables without centrality/multiplicity
template <>
inline bool centralitySelection<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator const& collision, float& centmult)
{
  return centralitySelectionNoMult(collision, centmult);
}

/// \brief Centrality selection for detector level collision tables with centrality/multiplicity
template <>
inline bool centralitySelection<soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator const& collision, float& centmult)
{
  return centralitySelectionMult(collision, centmult);
}

/// \brief Centrality selection for detector level collision tables with Run 2 centrality/multiplicity
template <>
inline bool centralitySelection<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator const& collision, float& centmult)
{
  return centralitySelectionMult(collision, centmult);
}

/// \brief Centrality selection for generator level collision table
template <>
inline bool centralitySelection<aod::McCollision>(aod::McCollision const&, float& centmult)
{
  if (centmult < 100 && 0 < centmult) {
    return true;
  } else {
    return false;
  }
}

//////////////////////////////////////////////////////////////////////////////////
/// Occupancy selection
//////////////////////////////////////////////////////////////////////////////////

/// \brief select on the collision occupancy
/// \return true if collison passes the occupancy cut false otherwise
template <typename CollisionObject>
inline bool selectOnOccupancy(CollisionObject collision)
{
  switch (fOccupancyEstimation) {
    case kNOOCC:
      collisionFlags.set(kOCCUPANCYBIT);
      return true;
    case kTRACKSOCC:
      if ((fMinOccupancy <= collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < fMaxOccupancy)) {
        collisionFlags.set(kOCCUPANCYBIT);
        return true;
      } else {
        return false;
      }
    case kFT0COCC:
      if ((fMinOccupancy <= collision.ft0cOccupancyInTimeRange()) && (collision.ft0cOccupancyInTimeRange() < fMaxOccupancy)) {
        collisionFlags.set(kOCCUPANCYBIT);
        return true;
      } else {
        return false;
      }
    default:
      return false;
  }
}

/// \brief Occupancy selection by default: unknown subscribed collision table
template <typename CollisionObject>
inline bool occupancySelection(CollisionObject const&)
{
  LOGF(fatal, "Occupancy selection not implemented for this kind of collisions");
  return false;
}

/// \brief Occupancy selection for reconstructed and detector level collision tables with centrality/multiplicity information
template <>
inline bool occupancySelection<aod::CollisionEvSelCent>(aod::CollisionEvSelCent const& collision)
{
  return selectOnOccupancy(collision);
}

/// \brief Occupancy selection for reconstructed and detector level collision tables with Run 2 centrality/multiplicity information
template <>
inline bool occupancySelection<aod::CollisionEvSelRun2Cent>(aod::CollisionEvSelRun2Cent const&)
{
  return true;
}

/// \brief Occupancy selection for reconstructed and detector level collision tables without centrality/multiplicity information
template <>
inline bool occupancySelection<aod::CollisionEvSel>(aod::CollisionEvSel const& collision)
{
  return selectOnOccupancy(collision);
}

/// \brief Occupancy selection for detector level collision tables without centrality/multiplicity
template <>
inline bool occupancySelection<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator const& collision)
{
  return selectOnOccupancy(collision);
}

/// \brief Occupancy selection for detector level collision tables with centrality/multiplicity
template <>
inline bool occupancySelection<soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator const& collision)
{
  return selectOnOccupancy(collision);
}

/// \brief Occupancy selection for detector level collision tables with Run 2 centrality/multiplicity
template <>
inline bool occupancySelection<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator const&)
{
  return true;
}

/// \brief Occupancy selection for generator level collision table
template <>
inline bool occupancySelection<aod::McCollision>(aod::McCollision const&)
{
  return true;
}

//////////////////////////////////////////////////////////////////////////////////
/// ITS dead maps selection
//////////////////////////////////////////////////////////////////////////////////

/// \brief select on the ITS dead maps
/// \return true if the collision passes the ITS dead maps check false otherwise
template <typename CollisionObject>
inline bool selectOnItsDeadMaps(CollisionObject coll)
{
  auto checkFlag = [](auto flag, bool invert = false) {
    return invert ? !flag : flag;
  };
  switch (fItsDeadMapCheck) {
    case kNOCHECK:
      return true;
    case kGOODITSLAYER3:
      return checkFlag(coll.selection_bit(aod::evsel::kIsGoodITSLayer3));
    case kGOODITSLAYER0123:
      return checkFlag(coll.selection_bit(aod::evsel::kIsGoodITSLayer0123));
    case kGOODITSLAYERALL:
      return checkFlag(coll.selection_bit(aod::evsel::kIsGoodITSLayersAll));
    case kNOGOODITSLAYER3:
      return checkFlag(coll.selection_bit(aod::evsel::kIsGoodITSLayer3), true);
    case kNOGOODITSLAYER0123:
      return checkFlag(coll.selection_bit(aod::evsel::kIsGoodITSLayer0123), true);
    case kNOGOODITSLAYERALL:
      return checkFlag(coll.selection_bit(aod::evsel::kIsGoodITSLayersAll), true);
    default:
      return false;
  }
}

/// \brief Occupancy selection by default: unknown subscribed collision table
template <typename CollisionObject>
inline bool itsDeadMapsSelection(CollisionObject const&)
{
  LOGF(fatal, "Occupancy selection not implemented for this kind of collisions");
  return false;
}

/// \brief Occupancy selection for reconstructed and detector level collision tables with centrality/multiplicity information
template <>
inline bool itsDeadMapsSelection<aod::CollisionEvSelCent>(aod::CollisionEvSelCent const& collision)
{
  return selectOnItsDeadMaps(collision);
}

/// \brief Occupancy selection for reconstructed and detector level collision tables with Run 2 centrality/multiplicity information
template <>
inline bool itsDeadMapsSelection<aod::CollisionEvSelRun2Cent>(aod::CollisionEvSelRun2Cent const&)
{
  return true;
}

/// \brief Occupancy selection for reconstructed and detector level collision tables without centrality/multiplicity information
template <>
inline bool itsDeadMapsSelection<aod::CollisionEvSel>(aod::CollisionEvSel const& collision)
{
  return selectOnItsDeadMaps(collision);
}

/// \brief Occupancy selection for detector level collision tables without centrality/multiplicity
template <>
inline bool itsDeadMapsSelection<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator const& collision)
{
  return selectOnItsDeadMaps(collision);
}

/// \brief Occupancy selection for detector level collision tables with centrality/multiplicity
template <>
inline bool itsDeadMapsSelection<soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator const& collision)
{
  return selectOnItsDeadMaps(collision);
}

/// \brief Occupancy selection for detector level collision tables with Run 2 centrality/multiplicity
template <>
inline bool itsDeadMapsSelection<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator const&)
{
  return true;
}

/// \brief Occupancy selection for generator level collision table
template <>
inline bool itsDeadMapsSelection<aod::McCollision>(aod::McCollision const&)
{
  return true;
}

//////////////////////////////////////////////////////////////////////////////////
/// Event selection
//////////////////////////////////////////////////////////////////////////////////

template <typename CollisionObject>
inline bool isEventSelected(CollisionObject const& collision, float& centormult)
{
  collisionFlags.reset();
  collisionFlags.set(kIN);

  bool trigsel = triggerSelection(collision);

  bool occupancysel = occupancySelection(collision);

  bool itsdeadmapssel = itsDeadMapsSelection(collision);

  bool zvtxsel = false;
  /* TODO: vertex quality checks */
  if (zvtxlow < collision.posZ() && collision.posZ() < zvtxup) {
    if (onlyInOneSide) {
      if (collision.posZ() != 0.0) {
        /* if only one side, we accept collisions which have zvtx different than zero */
        zvtxsel = true;
        collisionFlags.set(kZVERTEXBIT);
      }
    } else {
      zvtxsel = true;
      collisionFlags.set(kZVERTEXBIT);
    }
  }

  bool centmultsel = centralitySelection(collision, centormult);

  bool accepted = trigsel && occupancysel && itsdeadmapssel && zvtxsel && centmultsel;

  if (accepted) {
    collisionFlags.set(kSELECTED);
  }

  return accepted;
}

//////////////////////////////////////////////////////////////////////////////////
/// Track selection
//////////////////////////////////////////////////////////////////////////////////

struct TpcExcludeTrack {
  TpcExcludeTrack()
  {
    method = kNOEXCLUSION;
  }
  explicit TpcExcludeTrack(TpcExclusionMethod m) // o2-linter: disable=name/function-variable
  {
    switch (m) {
      case kNOEXCLUSION:
        method = m;
        break;
      case kSTATIC:
        if (phibinshift == 0.5f && phibins == 72) {
          method = m;
        } else {
          LOGF(fatal, "Static TPC exclusion method with bin shift: %.2f and number of bins %d. Please fix it", phibinshift, phibins);
        }
        break;
      case kDYNAMIC:
        method = m;
        break;
      default:
        LOGF(fatal, "Wrong TPC tracks exclusion method %d. Please, fix it", static_cast<int>(m));
    }
    philow = 0.0f;
    phiup = constants::math::TwoPI;
    phibinwidth = (phiup - philow) / static_cast<float>(phibins);
    phiup = phiup - phibinwidth * phibinshift;
    philow = philow - phibinwidth * phibinshift;
  }

  template <typename TrackObject>
  int getPhiBinIx(TrackObject const& track)
  {
    float phi = RecoDecay::constrainAngle(track.phi(), philow);
    return static_cast<int>((phi - philow) / phibinwidth);
  }

  template <typename TrackObject>
  bool exclude(TrackObject const& track)
  {
    constexpr int kNoOfTpcSectors = 18;
    constexpr float kTpcPhiSectorWidth = (constants::math::TwoPI) / kNoOfTpcSectors;

    switch (method) {
      case kNOEXCLUSION: {
        return false;
      } break;
      case kSTATIC: {
        int phiBinIx = getPhiBinIx(track);
        /* bins multiple of four have got sector border */
        if ((phiBinIx % 4) != 0) {
          return false;
        } else {
          return true;
        }
      } break;
      case kDYNAMIC: {
        float phiInTpcSector = std::fmod(track.phi(), kTpcPhiSectorWidth);
        if (track.sign() > 0) {
          return (phiInTpcSector < positiveUpCut->Eval(track.pt())) && (positiveLowCut->Eval(track.pt()) < phiInTpcSector);
        } else {
          return (phiInTpcSector < negativeUpCut->Eval(track.pt())) && (negativeLowCut->Eval(track.pt()) < phiInTpcSector);
        }
      } break;
      default:
        return false;
    }
  }

  void setCuts(std::string pLowCut, std::string pUpCut, std::string nLowCut, std::string nUpCut)
  {
    LOGF(info, "Setting the TPC exclusion cuts: pLow=%s, pUp=%s, nLow=%s, nUp=%s", pLowCut, pUpCut, nLowCut, nUpCut);
    positiveLowCut = new TF1("posLowCut", pLowCut.c_str(), ptlow, ptup);
    positiveUpCut = new TF1("posUpCut", pUpCut.c_str(), ptlow, ptup);
    negativeLowCut = new TF1("negLowCut", nLowCut.c_str(), ptlow, ptup);
    negativeUpCut = new TF1("negUpCut", nUpCut.c_str(), ptlow, ptup);
  }

  TpcExclusionMethod method = kNOEXCLUSION;
  float phibinwidth = 0.0;
  TF1* positiveLowCut = nullptr;
  TF1* positiveUpCut = nullptr;
  TF1* negativeLowCut = nullptr;
  TF1* negativeUpCut = nullptr;
};

template <typename TrackObject>
inline bool matchTrackType(TrackObject const& track)
{
  using namespace o2::aod::track;

  if (tracktype == 4) {
    // under tests MM track selection
    // see: https://indico.cern.ch/event/1383788/contributions/5816953/attachments/2805905/4896281/TrackSel_GlobalTracks_vs_MMTrackSel.pdf
    // it should be equivalent to this
    //       (track.passedDCAxy && track.passedDCAz && track.passedGoldenChi2) &&
    //       (track.passedITSNCls && track.passedITSChi2NDF && track.passedITSHits) &&
    //       (!track.hasTPC || (track.passedTPCNCls && track.passedTPCChi2NDF && track.passedTPCCrossedRowsOverNCls));
    return track.hasITS() && ((track.trackCutFlag() & TrackSelectionITS) == TrackSelectionITS) &&
           (!track.hasTPC() || ((track.trackCutFlag() & TrackSelectionTPC) == TrackSelectionTPC)) &&
           ((track.trackCutFlag() & TrackSelectionDCA) == TrackSelectionDCA);
  } else {
    for (size_t i = 0; i < trackFilters.size(); ++i) {
      auto const& filter = trackFilters[i];
      auto const& maxDcaZPtDep = maxDcaZPtDeps[i];
      if (filter->IsSelected(track)) {
        /* additional more stringent track cuts if needed */
        auto checkDca2Dcut = [&](auto const& track) {
          if (dca2Dcut) {
            if (track.dcaXY() * track.dcaXY() / maxDCAxy / maxDCAxy + track.dcaZ() * track.dcaZ() / maxDCAz / maxDCAz > 1) {
              return false;
            } else {
              return true;
            }
          } else {
            return true;
          }
        };
        auto checkDcaZcut = [&](auto const& track) {
          return ((maxDcaZPtDep) ? std::fabs(track.dcaZ()) <= maxDcaZPtDep(track.pt()) : true);
        };

        /* tight pT dependent DCAz cut */
        if (!checkDcaZcut(track)) {
          return false;
        }
        /* 2D DCA xy-o-z cut */
        if (!checkDca2Dcut(track)) {
          continue;
        }
        /* shared fraction of TPC clusters */
        if (!(track.tpcFractionSharedCls() < sharedTpcClusters)) {
          continue;
        }
        /* primary vertex contributor */
        if (tracktype == 10) {
          if (!track.isPVContributor()) {
            continue;
          }
        }
        return true;
      }
    }
    /* check extra less stringent selections */
    if (extraTrackFilter != nullptr) {
      if (extraTrackFilter->IsSelected(track)) {
        /* shared fraction of TPC clusters */
        if (!(track.tpcFractionSharedCls() < sharedTpcClusters)) {
          return false;
        }
        /* we require TOF */
        if (!(track.hasTOF())) {
          return false;
        }
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }
}

/// \brief Checks if the passed track is within the acceptance conditions of the analysis
/// \param track the track of interest
/// \return true if the track is in the acceptance, otherwise false
template <typename CollisionsObject, typename TrackObject>
inline bool inTheAcceptance(TrackObject const& track)
{
  /* the side on which the collision happened */
  float side = track.template collision_as<CollisionsObject>().posZ();

  /* overall minimum momentum cut for the analysis */
  if (!(overallminp < track.p())) {
    return false;
  }

  if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
    if (track.mcParticleId() < 0) {
      return false;
    }
  }

  /* check the side of the collision and decide if one side */
  if (onlyInOneSide) {
    if (side < 0) {
      if (track.eta() >= 0.0) {
        return false;
      }
    } else if (side > 0) {
      if (track.eta() <= 0.0) {
        return false;
      }
    } else {
      /* if only one side we should not have collisions with zvtx = 0 */
      LOGF(fatal, "Selecting one side tracks with a zvtx zero collision");
    }
  }

  if (ptlow < track.pt() && track.pt() < ptup && etalow < track.eta() && track.eta() < etaup) {
    return !tpcExcluder.exclude(track);
  }
  return false;
}

/// \brief Accepts or not the passed track
/// \param track the track of interest
/// \return true if the track is accepted, otherwise false
template <typename CollisionsObject, typename TrackObject>
inline bool acceptTrack(TrackObject const& track)
{
  if (inTheAcceptance<CollisionsObject>(track)) {
    if (matchTrackType(track)) {
      return true;
    }
  }
  return false;
}

template <typename ParticleObject, typename MCCollisionObject>
void exploreMothers(ParticleObject& particle, MCCollisionObject& collision)
{
  for (const auto& m : particle.template mothers_as<aod::McParticles>()) {
    LOGF(info, "   mother index: %d", m.globalIndex());
    LOGF(info, "   Tracking back mother");
    LOGF(info, "   assigned collision Id: %d, looping on collision Id: %d", m.mcCollisionId(), collision.globalIndex());
    LOGF(info, "   index: %d, pdg code: %d", m.globalIndex(), m.pdgCode());
    LOGF(info, "   Passed  isPhysicalPrimary(): %s", m.isPhysicalPrimary() ? "YES" : "NO");

    exploreMothers(m, collision);
  }
}

inline float getCharge(float pdgCharge)
{
  float charge = (pdgCharge / 3 >= 1) ? 1.0 : ((pdgCharge / 3 <= -1) ? -1.0 : 0);
  return charge;
}

/// \brief Accepts or not the passed generated particle
/// \param track the particle of interest
/// \return `true` if the particle is accepted, `false` otherwise
template <typename ParticleObject, typename MCCollisionObject>
inline bool acceptParticle(ParticleObject& particle, MCCollisionObject const&)
{
  /* overall momentum cut */
  if (!(overallminp < particle.p())) {
    return false;
  }

  if (particle.isPhysicalPrimary()) {
    if ((particle.mcCollisionId() == 0) && traceCollId0) {
      LOGF(info, "Particle %d passed isPhysicalPrimary", particle.globalIndex());
    }

    if (ptlow < particle.pt() && particle.pt() < ptup && etalow < particle.eta() && particle.eta() < etaup) {
      return true;
    }
  } else {
    if ((particle.mcCollisionId() == 0) && traceCollId0) {
      LOGF(info, "Particle %d NOT passed isPhysicalPrimary", particle.globalIndex());
    }
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////////
/// PID
//////////////////////////////////////////////////////////////////////////////////

struct PIDSpeciesSelection {
  const std::vector<int> pdgcodes = {11, 13, 211, 321, 2212};
  const std::vector<std::string_view> spnames = {"e", "mu", "pi", "ka", "p"};
  const std::vector<std::string_view> sptitles = {"e", "#mu", "#pi", "K", "p"};
  const std::vector<std::string_view> spfnames = {"E", "Mu", "Pi", "Ka", "Pr"};
  const std::vector<std::string_view> spadjnames = {"Electron", "Muon", "Pion", "Kaon", "Proton"};
  const std::vector<double> spmasses = {o2::constants::physics::MassElectron, o2::constants::physics::MassMuon, o2::constants::physics::MassPionCharged, o2::constants::physics::MassKaonCharged, o2::constants::physics::MassProton};
  const std::vector<std::string_view> chadjnames = {"P", "M"};
  const char* hadname = "h";
  const char* hadtitle = "h";
  const char* hadfname = "Ha";
  uint getNSpecies() { return config.size(); }
  const char* getSpeciesName(uint8_t ix) { return spnames[species[ix]].data(); }
  const char* getSpeciesTitle(uint8_t ix) { return sptitles[species[ix]].data(); }
  const char* getSpeciesFName(uint8_t ix) { return spfnames[species[ix]].data(); }
  double getSpeciesMass(uint8_t ix) { return spmasses[species[ix]]; }
  const char* getHadName() { return hadname; }
  const char* getHadTitle() { return hadtitle; }
  const char* getHadFName() { return hadfname; }
  bool isSpeciesBeingSelected(uint8_t sp) { return std::find(species.begin(), species.end(), sp) != species.end(); }
  bool isGlobalSpecies(uint8_t isp, o2::track::PID::ID glsp) { return species[isp] == glsp; }
  void storePIDAdjustments(TList* lst)
  {
    auto storedetectorwithcharge = [&](auto& detectorstore, auto detectorname, auto charge) {
      for (uint isp = 0; isp < spadjnames.size(); ++isp) {
        TString fullhname = TString::Format("%s%s%s_Difference", detectorname, spadjnames[isp].data(), charge);
        detectorstore[isp] = static_cast<TH1*>(lst->FindObject(fullhname.Data()));
      }
    };
    auto reportadjdetectorwithcharge = [&](auto& detectorstore, auto detectorname, auto charge) {
      for (uint isp = 0; isp < spadjnames.size(); ++isp) {
        if (detectorstore[isp] != nullptr) {
          LOGF(info, "Stored nsigmas adjust for detector %s and species %s%s in histogram %s", detectorname, spadjnames[isp].data(), charge, detectorstore[isp]->GetName());
        }
      }
    };
    storedetectorwithcharge(tpcnsigmasshiftpos, "TPC", "P");
    storedetectorwithcharge(tofnsigmasshiftpos, "TOF", "P");
    storedetectorwithcharge(tpcnsigmasshiftneg, "TPC", "M");
    storedetectorwithcharge(tofnsigmasshiftneg, "TOF", "M");

    reportadjdetectorwithcharge(tpcnsigmasshiftpos, "TPC", "P");
    reportadjdetectorwithcharge(tofnsigmasshiftpos, "TOF", "P");
    reportadjdetectorwithcharge(tpcnsigmasshiftneg, "TPC", "M");
    reportadjdetectorwithcharge(tofnsigmasshiftneg, "TOF", "M");
  }
  void addSpecies(uint8_t sp, o2::analysis::TrackSelectionPIDCfg* incfg)
  {
    o2::analysis::TrackSelectionPIDCfg* cfg = new o2::analysis::TrackSelectionPIDCfg(*incfg);
    config.push_back(cfg);
    species.push_back(sp);

    auto last = config[config.size() - 1];
    uint8_t lastsp = species[config.size() - 1];
    LOGF(info, "Inserted species %d with", lastsp);
    LOGF(info, "  minTPC nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMinNSigmasTPC[0], last->mMinNSigmasTPC[1], last->mMinNSigmasTPC[2], last->mMinNSigmasTPC[3], last->mMinNSigmasTPC[4]);
    LOGF(info, "  maxTPC nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMaxNSigmasTPC[0], last->mMaxNSigmasTPC[1], last->mMaxNSigmasTPC[2], last->mMaxNSigmasTPC[3], last->mMaxNSigmasTPC[4]);
    LOGF(info, "  minTOF nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMinNSigmasTOF[0], last->mMinNSigmasTOF[1], last->mMinNSigmasTOF[2], last->mMinNSigmasTOF[3], last->mMinNSigmasTOF[4]);
    LOGF(info, "  maxTOF nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMaxNSigmasTOF[0], last->mMaxNSigmasTOF[1], last->mMaxNSigmasTOF[2], last->mMaxNSigmasTOF[3], last->mMaxNSigmasTOF[4]);
    LOGF(info, "  %.1f < pT < %.1f", last->mPtMin, last->mPtMax);
  }
  void addExcludedSpecies(uint8_t sp, const o2::analysis::TrackSelectionPIDCfg* incfg)
  {
    o2::analysis::TrackSelectionPIDCfg* cfg = new o2::analysis::TrackSelectionPIDCfg(*incfg);
    configexclude.push_back(cfg);
    speciesexclude.push_back(sp);
    auto last = configexclude[configexclude.size() - 1];
    uint8_t lastsp = speciesexclude[configexclude.size() - 1];

    LOGF(info, "Inserted species %d for exclusion with", lastsp);
    LOGF(info, "  minTPC nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMinNSigmasTPC[0], last->mMinNSigmasTPC[1], last->mMinNSigmasTPC[2], last->mMinNSigmasTPC[3], last->mMinNSigmasTPC[4]);
    LOGF(info, "  maxTPC nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMaxNSigmasTPC[0], last->mMaxNSigmasTPC[1], last->mMaxNSigmasTPC[2], last->mMaxNSigmasTPC[3], last->mMaxNSigmasTPC[4]);
    LOGF(info, "  minTOF nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMinNSigmasTOF[0], last->mMinNSigmasTOF[1], last->mMinNSigmasTOF[2], last->mMinNSigmasTOF[3], last->mMinNSigmasTOF[4]);
    LOGF(info, "  maxTOF nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMaxNSigmasTOF[0], last->mMaxNSigmasTOF[1], last->mMaxNSigmasTOF[2], last->mMaxNSigmasTOF[3], last->mMaxNSigmasTOF[4]);
    LOGF(info, "  %.1f < pT < %.1f", last->mPtMin, last->mPtMax);
  }
  template <StrongDebugging outdebug, typename TrackObject>
  int8_t whichSpecies(TrackObject const& track)
  {
    TString debuginfo;
    std::vector<float> tpcnsigmas = {track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::vector<float> tofnsigmas = {track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};

    auto outmomentumdebug = [&]() {
      if constexpr (outdebug != 0) {
        debuginfo += TString::Format("%.5f,%.5f,%.5f,%d,%.2f,%.4f,", track.p(), track.tpcInnerParam(), track.pt(), track.hasTOF() ? 1 : 0, track.tpcSignal(), track.beta());
      }
    };
    auto outnsigmasdebug = [&]() {
      if constexpr (outdebug != 0) {
        for (auto const& tpcn : tpcnsigmas) {
          debuginfo += TString::Format("%.4f,", tpcn);
        }
        for (auto const& tofn : tofnsigmas) {
          debuginfo += TString::Format("%.4f,", tofn);
        }
      }
    };

    /* out debug if needed */
    outmomentumdebug();
    /* out debug if needed */
    outnsigmasdebug();

    auto closeTo = [](auto& values, auto& mindet, auto& maxdet, uint8_t sp) {
      if (mindet[sp] <= values[sp] && values[sp] < maxdet[sp]) {
        return true;
      } else {
        return false;
      }
    };
    auto awayFrom = [](auto& values, auto& mindet, auto& maxdet, uint8_t sp) {
      for (int ix = 0; ix < 5; ix++) {
        if (ix != sp) {
          if (mindet[ix] <= values[ix] && values[ix] < maxdet[ix]) {
            return false;
          }
        } else {
          continue;
        }
      }
      return true;
    };
    auto closeToTPC = [&](auto& config, uint8_t sp) {
      return closeTo(tpcnsigmas, config->mMinNSigmasTPC, config->mMaxNSigmasTPC, sp);
    };
    auto awayFromTPC = [&](auto& config, uint8_t sp) {
      return awayFrom(tpcnsigmas, config->mMinNSigmasTPC, config->mMaxNSigmasTPC, sp);
    };
    auto closeToTOF = [&](auto& config, uint8_t sp) {
      return closeTo(tofnsigmas, config->mMinNSigmasTOF, config->mMaxNSigmasTOF, sp);
    };
    auto closeToTPCTOF = [&](auto& config, uint8_t sp) {
      if (config->m2Dcut) {
        float a = (config->mMaxNSigmasTPC[sp] - config->mMinNSigmasTPC[sp]) / 2.0;
        float b = (config->mMaxNSigmasTOF[sp] - config->mMinNSigmasTOF[sp]) / 2.0;
        float oa = (config->mMaxNSigmasTPC[sp] + config->mMinNSigmasTPC[sp]) / 2.0;
        float ob = (config->mMaxNSigmasTOF[sp] + config->mMinNSigmasTOF[sp]) / 2.0;
        float vtpc = tpcnsigmas[sp] - oa;
        float vtof = tofnsigmas[sp] - ob;
        return ((vtpc * vtpc / a / a + vtof * vtof / b / b) < 1);
      } else {
        return closeToTPC(config, sp) && closeToTOF(config, sp);
      }
    };
    auto awayFromTPCTOF = [&](auto& config, uint8_t sp) {
      for (uint8_t ix = 0; ix < 5; ++ix) {
        if (ix != sp) {
          if (closeToTPCTOF(config, ix)) {
            return false;
          }
        } else {
          continue;
        }
      }
      return true;
    };
    auto aboveThreshold = [&](auto& config) {
      return ((config->mPThreshold > 0.0) && (config->mPThreshold < track.p()));
    };
    auto isA = [&](auto& config, uint8_t sp) {
      if (track.hasTOF()) {
        return closeToTPCTOF(config, sp) && awayFromTPCTOF(config, sp);
      } else {
        if (aboveThreshold(config)) {
          if (config->mRequireTOF) {
            return false;
          }
        }
      }
      /* we are here                                                */
      /* - below the threshold                                      */
      /* - above the threshold without TOF and without requiring it */
      /* so we check only the TPC information                       */
      return closeToTPC(config, sp) && awayFromTPC(config, sp);
    };

    auto adjustnsigmas = [&]() {
      if (track.sign() > 0) {
        for (uint isp = 0; isp < spnames.size(); ++isp) {
          if (tpcnsigmasshiftpos[isp] != nullptr) {
            TH1* h = tpcnsigmasshiftpos[isp];
            tpcnsigmas[isp] -= h->GetBinContent(h->GetXaxis()->FindFixBin(track.p()));
          }
          if (tofnsigmasshiftpos[isp] != nullptr) {
            TH1* h = tofnsigmasshiftpos[isp];
            tofnsigmas[isp] -= h->GetBinContent(h->GetXaxis()->FindFixBin(track.p()));
          }
        }
      } else {
        for (uint isp = 0; isp < spnames.size(); ++isp) {
          if (tpcnsigmasshiftneg[isp] != nullptr) {
            TH1* h = tpcnsigmasshiftneg[isp];
            tpcnsigmas[isp] -= h->GetBinContent(h->GetXaxis()->FindFixBin(track.p()));
          }
          if (tofnsigmasshiftneg[isp] != nullptr) {
            TH1* h = tofnsigmasshiftneg[isp];
            tofnsigmas[isp] -= h->GetBinContent(h->GetXaxis()->FindFixBin(track.p()));
          }
        }
      }
    };

    auto outpiddebug = [&](int code, int pid, int pid2) {
      if constexpr (outdebug != 0) {
        int truepid = -1;
        int isphysicalprimary = -1;
        int process = -1;
        if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
          if (!(track.mcParticleId() < 0)) {
            auto particle = track.template mcParticle_as<aod::McParticles>();
            truepid = particle.pdgCode();
            isphysicalprimary = particle.isPhysicalPrimary() ? 1 : 0;
            process = particle.getProcess();
          }
        }
        debuginfo += TString::Format("%d,%d,%d,%d,%d,%d\n", code, pid, pid2, truepid, isphysicalprimary, process);
        debugstream << debuginfo;
      }
    };

    /* adjust the nsigmas values if appropriate */
    adjustnsigmas();
    /* out debug info if needed */
    outnsigmasdebug();

    /* let's first check the exclusion from the analysis */
    for (uint8_t ix = 0; ix < configexclude.size(); ++ix) {
      if (isA(configexclude[ix], speciesexclude[ix])) {
        /* out debug info if needed */
        outpiddebug(1, speciesexclude[ix], -1);
        return -(ix + 1);
      }
    }
    /* we don't exclude it so check which species if any required */
    if (config.size() > 0) {
      int8_t id = -127;
      for (uint8_t ix = 0; ix < config.size(); ++ix) {
        if (isA(config[ix], species[ix])) {
          if (id < 0) {
            id = ix;
          } else {
            /* out debug info if needed */
            outpiddebug(2, species[id], species[ix]);
            /* already identified once */
            return -127;
          }
        }
      }
      /* check a species transverse momentum cut */
      if (!(id < 0)) {
        /* identified */
        if ((track.pt() < config[id]->mPtMin) || config[id]->mPtMax < track.pt()) {
          /* rejected */
          outpiddebug(4, species[id], -1);
          return -127;
        }
      }
      /* out debug info if needed */
      if (id < 0) {
        /* not identified */
        outpiddebug(3, -1, -1);
      } else {
        /* identified */
        outpiddebug(0, species[id], -1);
      }
      return id;
    } else {
      outpiddebug(0, 0, -1);
      /* charged hadron */
      return 0;
    }
  }
  template <typename ParticleObject>
  int8_t whichTruthSpecies(ParticleObject part)
  {
    int pdgcode = std::abs(part.pdgCode());
    /* let's first check the exclusion from the analysis */
    for (uint8_t ix = 0; ix < configexclude.size(); ++ix) {
      if (pdgcode == pdgcodes[speciesexclude[ix]]) {
        return -(ix + 1);
      }
    }
    /* we don't exclude it so check which species if any required */
    if (config.size() > 0) {
      for (uint8_t ix = 0; ix < config.size(); ++ix) {
        if (pdgcode == pdgcodes[species[ix]]) {
          /* check a species transverse momentum cut */
          if ((part.pt() < config[ix]->mPtMin) || config[ix]->mPtMax < part.pt()) {
            /* rejected */
            return -127;
          }
          return ix;
        }
      }
      return -127;
    } else {
      return 0;
    }
  }

  template <typename ParticleObject>
  int8_t whichTruthPrimarySpecies(ParticleObject part)
  {
    if (part.isPhysicalPrimary()) {
      return whichTruthSpecies(part);
    } else {
      return -127;
    }
  }

  template <typename ParticleObject>
  int8_t whichTruthSecondarySpecies(ParticleObject part)
  {
    if (!(part.isPhysicalPrimary())) {
      return whichTruthSpecies(part);
    } else {
      return -127;
    }
  }

  template <typename ParticleObject>
  int8_t whichTruthSecFromDecaySpecies(ParticleObject part)
  {
    if (!(part.isPhysicalPrimary()) && (part.getProcess() == 4)) {
      return whichTruthSpecies(part);
    } else {
      return -127;
    }
  }

  template <typename ParticleObject>
  int8_t whichTruthSecFromMaterialSpecies(ParticleObject part)
  {
    if (!(part.isPhysicalPrimary()) && (part.getProcess() != 4)) {
      return whichTruthSpecies(part);
    } else {
      return -127;
    }
  }

  std::vector<const o2::analysis::TrackSelectionPIDCfg*> config;        ///< the PID selection configuration of the species to include in the analysis
  std::vector<uint8_t> species;                                         ///< the species index of the species to include in the analysis
  std::vector<const o2::analysis::TrackSelectionPIDCfg*> configexclude; ///< the PID selection configuration of the species to exclude from the analysis
  std::vector<uint8_t> speciesexclude;                                  ///< the species index of teh species to exclude from the analysis
  std::vector<TH1*> tpcnsigmasshiftpos{spnames.size(), nullptr};
  std::vector<TH1*> tpcnsigmasshiftneg{spnames.size(), nullptr};
  std::vector<TH1*> tofnsigmasshiftpos{spnames.size(), nullptr};
  std::vector<TH1*> tofnsigmasshiftneg{spnames.size(), nullptr};
};

} // namespace dptdptfilter
} // namespace analysis
} // namespace o2

#endif // PWGCF_TABLEPRODUCER_DPTDPTFILTER_H_
