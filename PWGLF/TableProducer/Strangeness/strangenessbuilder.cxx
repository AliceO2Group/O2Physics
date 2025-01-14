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
//
// Strangeness builder task 
// ========================
// 
// This task produces all tables that may be necessary for 
// strangeness analyses. A single device is provided to 
// ensure better computing resource (memory) management. 
//
//  process functions: 
//  -- processPreselectTPCPID ..: pre-selects TPC dE/dx-compatible candidates.
//  -- processRealData .........: use this OR processSimulation but NOT both
//  -- processSimulation .......: use this OR processRealData but NOT both
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "TableHelper.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/Utils/strangenessBuilderHelper.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

using namespace o2;
using namespace o2::framework;

static constexpr int nParameters = 1;
static const std::vector<std::string> tableNames{"V0Indices",           //.0 (standard analysis: V0Data)
                                                 "V0CoresBase",         //.1 (standard analyses: V0Data)
                                                 "V0Covs",              //.2
                                                 "CascIndices",         //.3 (standard analyses: CascData)
                                                 "KFCascIndices",       // 4 (standard analyses: KFCascData)
                                                 "TraCascIndices",      //.5 (standard analyses: TraCascData)
                                                 "StoredCascCores",     //.6 (standard analyses: CascData)
                                                 "StoredKFCascCores",   // 7 (standard analyses: KFCascData)
                                                 "StoredTraCascCores",  //.8 (standard analyses: TraCascData)
                                                 "CascCovs",            //.9
                                                 "KFCascCovs",          // 10
                                                 "TraCascCovs",         //.11
                                                 "V0TrackXs",           //.12
                                                 "CascTrackXs",         //.13
                                                 "CascBBs",             //.14
                                                 "V0DauCovs",           //.15 (requested: tracking studies)
                                                 "V0DauCovIUs",         //.16 (requested: tracking studies)
                                                 "V0TraPosAtDCAs",      //.17 (requested: tracking studies)
                                                 "V0TraPosAtIUs",       //.18 (requested: tracking studies)
                                                 "V0Ivanovs",           //.19 (requested: tracking studies)
                                                 "McV0Labels",          // 20 (MC/standard analysis)
                                                 "V0MCCores",           // 21 (MC)
                                                 "V0CoreMCLabels",      // 22 (MC)
                                                 "V0MCCollRefs",        // 23 (MC)
                                                 "McCascLabels",        // 24 (MC/standard analysis)
                                                 "McKFCascLabels",      // 25 (MC)
                                                 "McTraCascLabels",     // 26 (MC)
                                                 "McCascBBTags",        // 27 (MC)
                                                 "CascMCCores",         // 28 (MC)
                                                 "CascCoreMCLabels",    // 29 (MC)
                                                 "CascMCCollRefs",      // 30 (MC)
                                                 "StraCollision",       // 31 (derived)
                                                 "StraCollLabels",      // 32 (derived)
                                                 "StraMCCollisions",    // 33 (MC/derived)
                                                 "StraMCCollMults",     // 34 (MC/derived)
                                                 "StraCents",           // 35 (derived)
                                                 "StraEvSels",          // 36 (derived)
                                                 "StraStamps",          // 37 (derived)
                                                 "V0CollRefs",          // 38 (derived)
                                                 "CascCollRefs",        // 39 (derived)
                                                 "KFCascCollRefs",      // 40 (derived)
                                                 "TraCascCollRefs",     // 41 (derived)
                                                 "DauTrackExtras",      // 42 (derived)
                                                 "DauTrackMCIds",       // 43 (MC/derived)
                                                 "DauTrackTPCPIDs",     // 44 (derived)
                                                 "DauTrackTOFPIDs",     // 45 (derived)
                                                 "V0Extras",            // 46 (derived)
                                                 "CascExtras",          // 47 (derived)
                                                 "StraTrackExtras",     // 48 (derived)
                                                 "CascToTraRefs",       // 49 (derived)
                                                 "CascToKFRefs",        // 50 (derived)
                                                 "TraToCascRefs",       // 51 (derived)
                                                 "KFToCascRefs",        // 52 (derived)
                                                 "V0MCMothers",         // 53 (MC/derived)
                                                 "CascMCMothers",       // 54 (MC/derived)
                                                 "MotherMCParts",       // 55 (MC/derived)
                                                 "StraFT0AQVs",         // 56 (derived)
                                                 "StraFT0CQVs",         // 57 (derived)
                                                 "StraFT0MQVs",         // 58 (derived)
                                                 "StraFV0AQVs",         // 59 (derived)
                                                 "StraTPCQVs",          // 60 (derived)
                                                 "StraFT0CQVsEv",       // 61 (derived)
                                                 "StraZDCSP",           // 62 (derived)
                                                 "GeK0Short",           // 63 (MC/derived)
                                                 "GeLambda",            // 64 (MC/derived)
                                                 "GeAntiLambda",        // 65 (MC/derived)
                                                 "GeXiMinus",           // 66 (MC/derived)
                                                 "GeXiPlus",            // 67 (MC/derived)
                                                 "GeOmegaMinus",        // 68 (MC/derived)
                                                 "GeOmegaPlus",         // 69 (MC/derived)
                                                 "V0FoundTags",         // 70 (MC/derived)
                                                 "CascFoundTags",       // 71 (MC/derived)
                                                 "StraOrigins"          // 72 (derived)
                                                 };

static constexpr int nTablesConst = 73;

static const std::vector<std::string> parameterNames{"enable"};
static const int defaultParameters[nTablesConst][nParameters]{
  {-1},  {1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, //0-9
  {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, //10-19
  {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, //20-29
  {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, //30-39
  {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, //40-49
  {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, //50-59
  {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, //60-69
  {-1}, {-1}, {-1}                                            //70-72
  };

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;

// For dE/dx association in pre-selection
using TracksExtraWithPID = soa::Join<aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullHe>;

struct StrangenessBuilder {
  // helper object
  o2::pwglf::strangenessBuilderHelper straHelper; 

  // table index : match order above
  enum tableIndex { kV0Indices = 0,
                    kV0CoresBase,
                    kV0Covs,
                    kCascIndices,
                    kKFCascIndices,
                    kTraCascIndices,
                    kStoredCascCores,
                    kStoredKFCascCores,
                    kStoredTraCascCores,
                    kCascCovs,
                    kKFCascCovs,
                    kTraCascCovs,
                    kV0TrackXs,
                    kCascTrackXs,
                    kCascBBs,
                    kV0DauCovs,
                    kV0DauCovIUs,
                    kV0TraPosAtDCAs,
                    kV0TraPosAtIUs,
                    kV0Ivanovs,
                    kMcV0Labels,
                    kV0MCCores,
                    kV0CoreMCLabels,
                    kV0MCCollRefs,
                    kMcCascLabels,
                    kMcKFCascLabels,
                    kMcTraCascLabels,
                    kMcCascBBTags,
                    kCascMCCores,
                    kCascCoreMCLabels,
                    kCascMCCollRefs,
                    kStraCollision,
                    kStraCollLabels,
                    kStraMCCollisions,
                    kStraMCCollMults,
                    kStraCents,
                    kStraEvSels,
                    kStraStamps,
                    kV0CollRefs,
                    kCascCollRefs,
                    kKFCascCollRefs,
                    kTraCascCollRefs,
                    kDauTrackExtras,
                    kDauTrackMCIds,
                    kDauTrackTPCPIDs,
                    kDauTrackTOFPIDs,
                    kV0Extras,
                    kCascExtras,
                    kStraTrackExtras,
                    kCascToTraRefs,
                    kCascToKFRefs,
                    kTraToCascRefs,
                    kKFToCascRefs,
                    kV0MCMothers,
                    kCascMCMothers,
                    kMotherMCParts,
                    kStraFT0AQVs,
                    kStraFT0CQVs,
                    kStraFT0MQVs,
                    kStraFV0AQVs,
                    kStraTPCQVs,
                    kStraFT0CQVsEv,
                    kStraZDCSP,
                    kGeK0Short,
                    kGeLambda,
                    kGeAntiLambda,
                    kGeXiMinus,
                    kGeXiPlus,
                    kGeOmegaMinus,
                    kGeOmegaPlus,
                    kV0FoundTags,
                    kCascFoundTags,
                    kStraOrigins,
                    nTables};

  //__________________________________________________
  // V0 tables
  Produces<aod::V0Indices> v0indices;         // standard part of V0Datas
  Produces<aod::V0CoresBase> v0cores;         // standard part of V0Datas
  Produces<aod::V0Covs> v0covs;               // for decay chain reco

  //__________________________________________________
  // cascade tables
  Produces<aod::CascIndices> cascidx;                 // standard part of CascDatas
  Produces<aod::KFCascIndices> kfcascidx;             // standard part of KFCascDatas
  Produces<aod::TraCascIndices> trackedcascidx;       // standard part of TraCascDatas
  Produces<aod::StoredCascCores> cascdata;            // standard part of CascDatas
  Produces<aod::StoredKFCascCores> kfcascdata;        // standard part of KFCascDatas
  Produces<aod::StoredTraCascCores> trackedcascdata;  // standard part of TraCascDatas
  Produces<aod::CascCovs> casccovs;                   // for decay chain reco
  Produces<aod::KFCascCovs> kfcasccovs;               // for decay chain reco
  Produces<aod::TraCascCovs> tracasccovs;             // for decay chain reco

  //__________________________________________________
  // interlink tables
  Produces<aod::V0DataLink> v0dataLink;           // de-refs V0s -> V0Data
  Produces<aod::CascDataLink> cascdataLink;       // de-refs Cascades -> CascData
  Produces<aod::KFCascDataLink> kfcascdataLink;   // de-refs Cascades -> KFCascData
  Produces<aod::TraCascDataLink> tracascdataLink; // de-refs Cascades -> TraCascData

  //__________________________________________________
  // secondary auxiliary tables
  Produces<aod::V0TrackXs> v0trackXs;       // for decay chain reco
  Produces<aod::CascTrackXs> cascTrackXs;   // for decay chain reco

  // further auxiliary / optional if desired
  Produces<aod::CascBBs> cascbb;      
  Produces<aod::V0DauCovs> v0daucovs;            // covariances of daughter tracks
  Produces<aod::V0DauCovIUs> v0daucovIUs;        // covariances of daughter tracks
  Produces<aod::V0TraPosAtDCAs> v0dauPositions;  // auxiliary debug information
  Produces<aod::V0TraPosAtIUs> v0dauPositionsIU; // auxiliary debug information
  Produces<aod::V0Ivanovs> v0ivanovs;            // information for Marian's tests

  //__________________________________________________
  // MC information: V0
  Produces<aod::McV0Labels> v0labels;           // MC labels for V0s
  Produces<aod::V0MCCores> v0mccores;           // mc info storage
  Produces<aod::V0CoreMCLabels> v0CoreMCLabels; // interlink V0Cores -> V0MCCores
  Produces<aod::V0MCCollRefs> v0mccollref;      // references collisions from V0MCCores
  
  // MC information: Cascades
  Produces<aod::McCascLabels> casclabels;            // MC labels for cascades
  Produces<aod::McKFCascLabels> kfcasclabels;        // MC labels for KF cascades
  Produces<aod::McTraCascLabels> tracasclabels;      // MC labels for tracked cascades
  Produces<aod::McCascBBTags> bbtags;                // bb tags (inv structure tagging in mc)
  Produces<aod::CascMCCores> cascmccores;            // mc info storage
  Produces<aod::CascCoreMCLabels> cascCoreMClabels;  // interlink CascCores -> CascMCCores
  Produces<aod::CascMCCollRefs> cascmccollrefs;      // references MC collisions from MC cascades

  //__________________________________________________
  // fundamental building blocks of derived data
  Produces<aod::StraCollision> strangeColl;        // characterises collisions
  Produces<aod::StraCollLabels> strangeCollLabels; // characterises collisions
  Produces<aod::StraMCCollisions> strangeMCColl;   // characterises collisions / MC
  Produces<aod::StraMCCollMults> strangeMCMults;   // characterises collisions / MC mults
  Produces<aod::StraCents> strangeCents;           // characterises collisions / centrality
  Produces<aod::StraEvSels> strangeEvSels;         // characterises collisions / centrality / sel8 selection
  Produces<aod::StraStamps> strangeStamps;         // provides timestamps, run numbers
  Produces<aod::V0CollRefs> v0collref;             // references collisions from V0s
  Produces<aod::CascCollRefs> casccollref;         // references collisions from cascades
  Produces<aod::KFCascCollRefs> kfcasccollref;     // references collisions from KF cascades
  Produces<aod::TraCascCollRefs> tracasccollref;   // references collisions from tracked cascades

  //__________________________________________________
  // track extra references
  Produces<aod::DauTrackExtras> dauTrackExtras;   // daughter track detector properties
  Produces<aod::DauTrackMCIds> dauTrackMCIds;     // daughter track MC Particle ID
  Produces<aod::DauTrackTPCPIDs> dauTrackTPCPIDs; // daughter track TPC PID
  Produces<aod::DauTrackTOFPIDs> dauTrackTOFPIDs; // daughter track TOF PID
  Produces<aod::V0Extras> v0Extras;               // references DauTracks from V0s
  Produces<aod::CascExtras> cascExtras;           // references DauTracks from cascades
  Produces<aod::StraTrackExtras> straTrackExtras; // references DauTracks from tracked cascades (for the actual tracked cascade, not its daughters)

  //__________________________________________________
  // cascade interlinks
  Produces<aod::CascToTraRefs> cascToTraRefs; // cascades -> tracked
  Produces<aod::CascToKFRefs> cascToKFRefs;   // cascades -> KF
  Produces<aod::TraToCascRefs> traToCascRefs; // tracked -> cascades
  Produces<aod::KFToCascRefs> kfToCascRefs;   // KF -> cascades

  //__________________________________________________
  // mother information
  Produces<aod::V0MCMothers> v0mothers;       // V0 mother references
  Produces<aod::CascMCMothers> cascmothers;   // casc mother references
  Produces<aod::MotherMCParts> motherMCParts; // mc particles for mothers

  //__________________________________________________
  // Q-vectors
  Produces<aod::StraFT0AQVs> StraFT0AQVs;     // FT0A Q-vector
  Produces<aod::StraFT0CQVs> StraFT0CQVs;     // FT0C Q-vector
  Produces<aod::StraFT0MQVs> StraFT0MQVs;     // FT0M Q-vector
  Produces<aod::StraFV0AQVs> StraFV0AQVs;     // FV0A Q-vector
  Produces<aod::StraTPCQVs> StraTPCQVs;       // TPC Q-vector
  Produces<aod::StraFT0CQVsEv> StraFT0CQVsEv; // events used to compute FT0C Q-vector (LF)
  Produces<aod::StraZDCSP> StraZDCSP;         // ZDC Sums and Products

  //__________________________________________________
  // Generated binned data
  // this is a hack while the system does not do better
  Produces<aod::GeK0Short> geK0Short;
  Produces<aod::GeLambda> geLambda;
  Produces<aod::GeAntiLambda> geAntiLambda;
  Produces<aod::GeXiMinus> geXiMinus;
  Produces<aod::GeXiPlus> geXiPlus;
  Produces<aod::GeOmegaMinus> geOmegaMinus;
  Produces<aod::GeOmegaPlus> geOmegaPlus;

  //__________________________________________________
  // Found tags for findable exercise
  Produces<aod::V0FoundTags> v0FoundTags;
  Produces<aod::CascFoundTags> cascFoundTags;

  //__________________________________________________
  // Debug
  Produces<aod::StraOrigins> straOrigin;

  Configurable<LabeledArray<int>> enabledTables{"enabledTables",
                                                {defaultParameters[0], nTables, nParameters, tableNames, parameterNames},
                                                "Produce this table: -1 for autodetect; otherwise, 0/1 is false/true"};
  std::vector<int> mEnabledTables; // Vector of enabled tables

  // CCDB options
  struct : ConfigurableGroup {
    std::string prefix = "ccdb";
    Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  } ccdbConfigurations;

  // V0 building options
  struct : ConfigurableGroup {
    std::string prefix = "v0BuilderOpts";
    Configurable<bool> generatePhotonCandidates{"generatePhotonCandidates", false, "generate gamma conversion candidates (V0s using TPC-only tracks)"};
  } v0BuilderOpts;

  // cascade building options
  struct : ConfigurableGroup {
    std::string prefix = "cascadeBuilderOpts";
    Configurable<bool> useCascadeMomentumAtPrimVtx{"useCascadeMomentumAtPrimVtx", false, "use cascade momentum at PV"};

    // KF building specific
    Configurable<bool> kfTuneForOmega{"kfTuneForOmega", false, "if enabled, take main cascade properties from Omega fit instead of Xi fit (= default)"};
    Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "KF Construct Method"};
    Configurable<bool> kfUseV0MassConstraint{"kfUseV0MassConstraint", true, "KF: use Lambda mass constraint"};
    Configurable<bool> kfUseCascadeMassConstraint{"kfUseCascadeMassConstraint", false, "KF: use Cascade mass constraint - WARNING: not adequate for inv mass analysis of Xi"};
    Configurable<bool> kfDoDCAFitterPreMinimV0{"kfDoDCAFitterPreMinimV0", true, "KF: do DCAFitter pre-optimization before KF fit to include material corrections for V0"};
    Configurable<bool> kfDoDCAFitterPreMinimCasc{"kfDoDCAFitterPreMinimCasc", true, "KF: do DCAFitter pre-optimization before KF fit to include material corrections for Xi"};
  } cascadeBuilderOpts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int mRunNumber;
  o2::base::MatLayerCylSet* lut = nullptr;

  // for tagging V0s used in cascades
  std::vector<o2::pwglf::v0candidate> v0sFromCascades; // Vector of v0 candidates used in cascades
  std::vector<int> v0Map;                              // index to relate V0s -> v0sFromCascades

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    mRunNumber = 0;
    
    mEnabledTables.resize(nTables, 0);

    for (int i = 0; i < nTables; i++) {
      int f = enabledTables->get(tableNames[i].c_str(), "enable");
      if (f == 1) {
        mEnabledTables[i] = 1;
        LOGF(info, "Enabled table: %s", tableNames[i].c_str());
      }
    }

    ccdb->setURL(ccdbConfigurations.ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  bool initCCDB(aod::BCsWithTimestamps const& bcs, aod::Collisions const& collisions)
  {
    auto bc = collisions.size() ? collisions.begin().bc_as<aod::BCsWithTimestamps>() : bcs.begin();
    if (!bcs.size()) {
      LOGF(warn, "No BC found, skipping this DF.");
      return false; // signal to skip this DF
    }

    if (mRunNumber == bc.runNumber()) {
      return true;
    }

    auto timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (doprocessRealDataRun2) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbConfigurations.grpPath, timestamp);
      if (!grpo) {
        LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpPath << " of object GRPObject for timestamp " << timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbConfigurations.grpmagPath, timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpmagPath << " of object GRPMagField for timestamp " << timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
    }
    // Fetch magnetic field from ccdb for current collision
    auto magneticField = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << magneticField << " kG";

    // Set magnetic field value once known
    straHelper.fitter.setBz(magneticField);

    // acquire LUT for this timestamp
    LOG(info) << "Loading material look-up table for timestamp: " << timestamp;
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(ccdbConfigurations.lutPath, timestamp));
    o2::base::Propagator::Instance()->setMatLUT(lut);
    straHelper.lut = lut;

    LOG(info) << "Fully configured for run: " << bc.runNumber();
    // mmark this run as configured  
    mRunNumber = bc.runNumber();

    return true;
  }

  //__________________________________________________
  template <typename TV0s, typename TCascades>
  void markV0sUsedInCascades(TV0s const& v0s, TCascades const& cascades)
  {
    v0sFromCascades.clear();
    v0Map.resize(v0s.size(), -2); // marks not used
    for (auto& cascade : cascades) {
      v0Map[cascade.v0Id()] = -1; // marks used (but isn't the index of a properly built V0, which would be >= 0)
    }
  }

  template <class TTracks, typename TCollisions, typename TV0s>
  void buildV0s(TCollisions const& collisions, TV0s const& v0s)
  {
    int nV0s = 0;
    // Loops over all V0s in the time frame
    for (auto& v0 : v0s) {
      if(!mEnabledTables[kV0CoresBase] && v0Map[v0.globalIndex()] == -2){
        // this v0 hasn't been used by cascades and we're not generating V0s, so skip it
        v0dataLink(-1, -1);
        continue; 
      }

      // Get tracks and generate candidate
      auto const& collision = v0.collision();
      auto const& posTrack = v0.template posTrack_as<TTracks>();
      auto const& negTrack = v0.template negTrack_as<TTracks>();
      if(!straHelper.buildV0Candidate(collision, posTrack, negTrack, v0.isCollinearV0())){ 
        v0dataLink(-1, -1);
        continue;
      }
      nV0s++;
      if(v0Map[v0.globalIndex()]==-1){
        v0Map[v0.globalIndex()] = v0sFromCascades.size(); // provide actual valid index in buffer
        v0sFromCascades.push_back(straHelper.v0);
      }
      // fill requested cursors only if type is not 0
      if (v0.v0Type() == 1 || (v0.v0Type() == 1 && v0BuilderOpts.generatePhotonCandidates)) {
        if(mEnabledTables[kV0Indices]){
          // for referencing (especially - but not only - when using derived data)
          v0indices(v0.posTrackId(), v0.negTrackId(),
                    v0.collisionId(), v0.globalIndex());
        }
        if(mEnabledTables[kV0TrackXs]){
          // further decay chains may need this
          v0trackXs(straHelper.v0.positiveTrackX, straHelper.v0.negativeTrackX);
        }
        if(mEnabledTables[kV0CoresBase]){
          // standard analysis
          v0cores(straHelper.v0.position[0], straHelper.v0.position[1], straHelper.v0.position[2],
                  straHelper.v0.positiveMomentum[0], straHelper.v0.positiveMomentum[1], straHelper.v0.positiveMomentum[2],
                  straHelper.v0.negativeMomentum[0], straHelper.v0.negativeMomentum[1], straHelper.v0.negativeMomentum[2],
                  straHelper.v0.daughterDCA,
                  straHelper.v0.positiveDCAxy,
                  straHelper.v0.negativeDCAxy,
                  TMath::Cos(straHelper.v0.pointingAngle),
                  straHelper.v0.dcaXY,
                  v0.v0Type());
          v0dataLink(v0cores.lastIndex(), -1);
        }
        if(mEnabledTables[kV0TraPosAtDCAs]){
          // for tracking studies
          v0dauPositions(straHelper.v0.positivePosition[0], straHelper.v0.positivePosition[1], straHelper.v0.positivePosition[2],
                         straHelper.v0.negativePosition[0], straHelper.v0.negativePosition[1], straHelper.v0.negativePosition[2]);

        }
        if(mEnabledTables[kV0TraPosAtIUs]){
          // for tracking studies
          std::array<float, 3> positivePositionIU;
          std::array<float, 3> negativePositionIU;
          o2::track::TrackPar positiveTrackParam = getTrackPar(posTrack);
          o2::track::TrackPar negativeTrackParam = getTrackPar(negTrack);
          positiveTrackParam.getXYZGlo(positivePositionIU);
          negativeTrackParam.getXYZGlo(negativePositionIU);
          v0dauPositionsIU(positivePositionIU[0], positivePositionIU[1], positivePositionIU[2],
                           negativePositionIU[0], negativePositionIU[1], negativePositionIU[2]);
        }
      }
    }
    LOGF(info, "V0s in DF: %i, V0s built: %i, V0s built and buffered for cascades: %i.", v0s.size(), nV0s, v0sFromCascades.size());
  }

  template <class TTracks, typename TCollisions, typename TCascades>
  void buildCascades(TCollisions const& collisions, TCascades const& cascades)
  {
    if(!mEnabledTables[kStoredCascCores]){ 
      return; // don't do if no request for cascades in place
    }
    int nCascades = 0;
    // Loops over all V0s in the time frame
    for (auto& cascade : cascades) {
      // Get tracks and generate candidate
      auto const& collision = cascade.collision();
      auto const& v0 = cascade.v0();
      auto const& posTrack = v0.template posTrack_as<TTracks>();
      auto const& negTrack = v0.template negTrack_as<TTracks>();
      auto const& bachTrack = cascade.template bachelor_as<TTracks>();
      if(!straHelper.buildCascadeCandidate(collision, 
                                           v0sFromCascades[v0Map[v0.globalIndex()]],
                                           posTrack, 
                                           negTrack, 
                                           bachTrack, 
                                           mEnabledTables[kCascBBs],
                                           cascadeBuilderOpts.useCascadeMomentumAtPrimVtx,  
                                           mEnabledTables[kCascCovs])){
        cascdataLink(-1);
        continue; // didn't work out, skip
      }
      nCascades++;

      // generate analysis tables as required
      if(mEnabledTables[kCascIndices]){
        cascidx(cascade.globalIndex(),
                straHelper.cascade.positiveTrack, straHelper.cascade.negativeTrack,
                straHelper.cascade.bachelorTrack, straHelper.cascade.collisionId);
      }
      if (mEnabledTables[kStoredCascCores]){
        cascdata(straHelper.cascade.charge, straHelper.cascade.massXi, straHelper.cascade.massOmega,
                 straHelper.cascade.cascadePosition[0], straHelper.cascade.cascadePosition[1], straHelper.cascade.cascadePosition[2],
                 straHelper.cascade.v0Position[0], straHelper.cascade.v0Position[1], straHelper.cascade.v0Position[2],
                 straHelper.cascade.positiveMomentum[0], straHelper.cascade.positiveMomentum[1], straHelper.cascade.positiveMomentum[2],
                 straHelper.cascade.negativeMomentum[0], straHelper.cascade.negativeMomentum[1], straHelper.cascade.negativeMomentum[2],
                 straHelper.cascade.bachelorMomentum[0], straHelper.cascade.bachelorMomentum[1], straHelper.cascade.bachelorMomentum[2],
                 straHelper.cascade.cascadeMomentum[0], straHelper.cascade.cascadeMomentum[1], straHelper.cascade.cascadeMomentum[2],
                 straHelper.cascade.v0DaughterDCA, straHelper.cascade.cascadeDaughterDCA,
                 straHelper.cascade.positiveDCAxy, straHelper.cascade.negativeDCAxy,
                 straHelper.cascade.bachelorDCAxy, straHelper.cascade.cascadeDCAxy, straHelper.cascade.cascadeDCAz);
        // interlink always produced if cascades generated
        cascdataLink(cascdata.lastIndex());
      }

      if (mEnabledTables[kCascTrackXs]) {
        cascTrackXs(straHelper.cascade.positiveTrackX, straHelper.cascade.negativeTrackX, straHelper.cascade.bachelorTrackX);
      }
      if( mEnabledTables[kCascBBs]){
        cascbb(straHelper.cascade.bachBaryonCosPA, straHelper.cascade.bachBaryonDCAxyToPV);
      }
      if (mEnabledTables[kCascCovs]) {
        casccovs(straHelper.cascade.covariance);
      }
    }
      
    LOGF(info, "Cascades in DF: %i, cascades built: %i", cascades.size(), nCascades);
  }


  template <class TTracks, typename TCollisions, typename TCascades>
  void buildKFCascades(TCollisions const& collisions, TCascades const& cascades)
  {
    if(!mEnabledTables[kStoredKFCascCores]){ 
      return; // don't do if no request for cascades in place
    }
    int nCascades = 0;
    // Loops over all V0s in the time frame
    for (auto& cascade : cascades) {
      // Get tracks and generate candidate
      auto const& collision = cascade.collision();
      auto const& v0 = cascade.v0();
      auto const& posTrack = v0.template posTrack_as<TTracks>();
      auto const& negTrack = v0.template negTrack_as<TTracks>();
      auto const& bachTrack = cascade.template bachelor_as<TTracks>();
      if(!straHelper.buildCascadeCandidateWithKF(collision, 
                                                 posTrack, 
                                                 negTrack, 
                                                 bachTrack, 
                                                 mEnabledTables[kCascBBs], 
                                                 cascadeBuilderOpts.kfConstructMethod, 
                                                 cascadeBuilderOpts.kfTuneForOmega, 
                                                 cascadeBuilderOpts.kfUseV0MassConstraint, 
                                                 cascadeBuilderOpts.kfUseCascadeMassConstraint, 
                                                 cascadeBuilderOpts.kfDoDCAFitterPreMinimV0, 
                                                 cascadeBuilderOpts.kfDoDCAFitterPreMinimCasc)){
        kfcascdataLink(-1);
        continue; // didn't work out, skip
      }
      nCascades++;

      // generate analysis tables as required
      if(mEnabledTables[kKFCascIndices]){
        kfcascidx(cascade.globalIndex(),
                  straHelper.cascade.positiveTrack, straHelper.cascade.negativeTrack,
                  straHelper.cascade.bachelorTrack, straHelper.cascade.collisionId);
      }
      if (mEnabledTables[kStoredKFCascCores]){
        kfcascdata(straHelper.cascade.charge, straHelper.cascade.massXi, straHelper.cascade.massOmega,
                   straHelper.cascade.cascadePosition[0], straHelper.cascade.cascadePosition[1], straHelper.cascade.cascadePosition[2],
                   straHelper.cascade.v0Position[0], straHelper.cascade.v0Position[1], straHelper.cascade.v0Position[2],
                   straHelper.cascade.positivePosition[0], straHelper.cascade.positivePosition[1], straHelper.cascade.positivePosition[2],
                   straHelper.cascade.negativePosition[0], straHelper.cascade.negativePosition[1], straHelper.cascade.negativePosition[2],
                   straHelper.cascade.positiveMomentum[0], straHelper.cascade.positiveMomentum[1], straHelper.cascade.positiveMomentum[2],
                   straHelper.cascade.negativeMomentum[0], straHelper.cascade.negativeMomentum[1], straHelper.cascade.negativeMomentum[2],
                   straHelper.cascade.bachelorMomentum[0], straHelper.cascade.bachelorMomentum[1], straHelper.cascade.bachelorMomentum[2],
                   straHelper.cascade.v0Momentum[0], straHelper.cascade.v0Momentum[1], straHelper.cascade.v0Momentum[2],
                   straHelper.cascade.cascadeMomentum[0], straHelper.cascade.cascadeMomentum[1], straHelper.cascade.cascadeMomentum[2],
                   straHelper.cascade.v0DaughterDCA, straHelper.cascade.cascadeDaughterDCA,
                   straHelper.cascade.positiveDCAxy, straHelper.cascade.negativeDCAxy,
                   straHelper.cascade.bachelorDCAxy, straHelper.cascade.cascadeDCAxy, straHelper.cascade.cascadeDCAz,
                   straHelper.cascade.kfMLambda, straHelper.cascade.kfV0Chi2, straHelper.cascade.kfCascadeChi2);
        // interlink always produced if cascades generated
        kfcascdataLink(kfcascdata.lastIndex());
      }
      if (mEnabledTables[kKFCascCovs]) {
        kfcasccovs(straHelper.cascade.covariance, straHelper.cascade.kfTrackCovarianceV0, straHelper.cascade.kfTrackCovariancePos, straHelper.cascade.kfTrackCovarianceNeg);
      }
    }
      
    LOGF(info, "KF Cascades in DF: %i, KF cascades built: %i", cascades.size(), nCascades);
  }

  template <class TTracks, typename TCollisions, typename TStrangeTracks>
  void buildTrackedCascades(TCollisions const& collisions, TStrangeTracks const& cascadeTracks)
  {
    if(!mEnabledTables[kStoredTraCascCores]){ 
      return; // don't do if no request for cascades in place
    }
    int nCascades = 0;
    // Loops over all V0s in the time frame
    for (auto& cascadeTrack : cascadeTracks) {
      // Get tracks and generate candidate
      if (!cascadeTrack.has_track())
        continue; // safety (should be fine but depends on future stratrack dev)

      auto const& strangeTrack = cascadeTrack.template track_as<TTracks>();
      auto const& collision = strangeTrack.collision();
      auto const& cascade = strangeTrack.cascade();
      auto const& v0 = cascade.v0();
      auto const& posTrack = v0.template posTrack_as<TTracks>();
      auto const& negTrack = v0.template negTrack_as<TTracks>();
      auto const& bachTrack = cascade.template bachelor_as<TTracks>();
      if(!straHelper.buildCascadeCandidate(collision, 
                                           posTrack, 
                                           negTrack, 
                                           bachTrack, 
                                           mEnabledTables[kCascBBs],
                                           cascadeBuilderOpts.useCascadeMomentumAtPrimVtx,
                                           mEnabledTables[kCascCovs])){
        tracascdataLink(-1);
        continue; // didn't work out, skip
      }

      // recalculate DCAxy, DCAz with strange track
      auto strangeTrackParCov = getTrackParCov(strangeTrack);
      gpu::gpustd::array<float, 2> dcaInfo;
      strangeTrackParCov.setPID(o2::track::PID::XiMinus); // FIXME: not OK for omegas
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, strangeTrackParCov, 2.f, straHelper.fitter.getMatCorrType(), &dcaInfo);
      straHelper.cascade.cascadeDCAxy = dcaInfo[0]; 
      straHelper.cascade.cascadeDCAz = dcaInfo[1]; 

      // get momentum from strange track (should not be very different)
      strangeTrackParCov.getPxPyPzGlo(straHelper.cascade.cascadeMomentum);

      // accounting
      nCascades++;

      // generate analysis tables as required
      if(mEnabledTables[kTraCascIndices]){
        cascidx(cascade.globalIndex(),
                straHelper.cascade.positiveTrack, straHelper.cascade.negativeTrack,
                straHelper.cascade.bachelorTrack, straHelper.cascade.collisionId);
      }
      if (mEnabledTables[kStoredTraCascCores]){
        cascdata(straHelper.cascade.charge, cascadeTrack.xiMass(), cascadeTrack.omegaMass(),
                 cascadeTrack.decayX(), cascadeTrack.decayY(), cascadeTrack.decayZ(),
                 straHelper.cascade.v0Position[0], straHelper.cascade.v0Position[1], straHelper.cascade.v0Position[2],
                 straHelper.cascade.positiveMomentum[0], straHelper.cascade.positiveMomentum[1], straHelper.cascade.positiveMomentum[2],
                 straHelper.cascade.negativeMomentum[0], straHelper.cascade.negativeMomentum[1], straHelper.cascade.negativeMomentum[2],
                 straHelper.cascade.bachelorMomentum[0], straHelper.cascade.bachelorMomentum[1], straHelper.cascade.bachelorMomentum[2],
                 straHelper.cascade.cascadeMomentum[0], straHelper.cascade.cascadeMomentum[1], straHelper.cascade.cascadeMomentum[2],
                 straHelper.cascade.v0DaughterDCA, straHelper.cascade.cascadeDaughterDCA,
                 straHelper.cascade.positiveDCAxy, straHelper.cascade.negativeDCAxy,
                 straHelper.cascade.bachelorDCAxy, straHelper.cascade.cascadeDCAxy, straHelper.cascade.cascadeDCAz, 
                 cascadeTrack.matchingChi2(), cascadeTrack.topologyChi2(), cascadeTrack.itsClsSize());
        // interlink always produced if base core table generated
        tracascdataLink(cascdata.lastIndex());
      }
      if (mEnabledTables[kCascCovs]) {
        std::array<float, 21> traCovMat = {0.};
        strangeTrackParCov.getCovXYZPxPyPzGlo(traCovMat);
        float traCovMatArray[21];
        for (int ii = 0; ii < 21; ii++) {
          traCovMatArray[ii] = traCovMat[ii];
        }
        tracasccovs(traCovMatArray);
      }
    }
    LOGF(info, "Tracked cascades in DF: %i, tracked cascades built: %i", cascadeTracks.size(), nCascades);
  }

  //__________________________________________________
  template <typename TCollisions, typename TV0s, typename TCascades, typename TTrackedCascades, typename TTracks, typename TBCs>
  void dataProcess(TCollisions const& collisions, TV0s const& v0s, TCascades const& cascades, TTrackedCascades const& trackedCascades, TTracks const&, TBCs const& bcs)
  {
    if(!initCCDB(bcs, collisions)) return;

    // mark V0s that will be buffered for the cascade building
    markV0sUsedInCascades(v0s, cascades);

    // build V0s 
    buildV0s<TTracks>(collisions, v0s);
    
    // build cascades
    buildCascades<TTracks>(collisions, cascades);
    buildKFCascades<TTracks>(collisions, cascades);

    // build tracked cascades only if subscription is Run 3 like (doesn't exist in Run 2)
    if constexpr (requires { TTrackedCascades::iterator; }) {
      buildTrackedCascades<TTracks>(collisions, trackedCascades);
    }
  }

  void processRealData(aod::Collisions const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    dataProcess(collisions, v0s, cascades, trackedCascades, tracks, bcs);
  }

  void processRealDataRun2(aod::Collisions const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExt const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    dataProcess(collisions, v0s, cascades, (TObject*) nullptr, tracks, bcs);
  }

  void processSimulationFindable(aod::Collisions const& collisions, aod::FindableV0s const& v0s, aod::Cascades const& cascades, FullTracksExtIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    dataProcess(collisions, v0s, cascades, (TObject*) nullptr, tracks, bcs);
  }

  PROCESS_SWITCH(StrangenessBuilder, processRealData, "process real data", true);
  PROCESS_SWITCH(StrangenessBuilder, processRealDataRun2, "process real data (Run 2)", false);
  PROCESS_SWITCH(StrangenessBuilder, processSimulationFindable, "process simulation findable (requires lambdakzeromcfinder)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<StrangenessBuilder>(cfgc)};
}
