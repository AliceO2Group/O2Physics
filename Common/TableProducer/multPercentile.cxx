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
/// \file multCentrality.cxx
/// \brief Produces multiplicity and centrality percentiles tables
///
/// \author ALICE
///


#include <vector>
#include <algorithm>
#include <map>
#include <string>

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "TableHelper.h"
#include "MetadataHelper.h"
#include "TList.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

MetadataHelper metadataInfo; // Metadata helper

namespace multiplicity{
static constexpr int kFV0Mults = 0;
static constexpr int kFT0Mults = 1;
static constexpr int kFDDMults = 2;
static constexpr int kZDCMults = 3;
static constexpr int kTrackletMults = 4;
static constexpr int kTPCMults = 5;
static constexpr int kPVMults = 6;
static constexpr int kMultsExtra = 7;
static constexpr int kMultSelections = 8;
static constexpr int kFV0MultZeqs = 9;
static constexpr int kFT0MultZeqs = 10;
static constexpr int kFDDMultZeqs = 11;
static constexpr int kPVMultZeqs = 12;
static constexpr int kMultMCExtras = 13;
static constexpr int kNTables = 14;

// Checking that the Zeq tables are after the normal ones
static_assert(kFV0Mults < kFV0MultZeqs);
static_assert(kFT0Mults < kFT0MultZeqs);
static_assert(kFDDMults < kFDDMultZeqs);
static_assert(kPVMults < kPVMultZeqs);


static const std::vector<std::string> tableNames{"FV0Mults",       // 0
    "FT0Mults",       // 1
    "FDDMults",       // 2
    "ZDCMults",       // 3
    "TrackletMults",  // 4
    "TPCMults",       // 5
    "PVMults",        // 6
    "MultsExtra",     // 7
    "MultSelections", // 8
    "FV0MultZeqs",    // 9
    "FT0MultZeqs",    // 10
    "FDDMultZeqs",    // 11
    "PVMultZeqs",     // 12
    "MultMCExtras"};  // 13
    static const int defaultParameters[kNTables][1]{{-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}};
}

namespace centrality{
    static constexpr int kRun2V0Ms = 0;
static constexpr int kRun2V0As = 1;
static constexpr int kRun2SPDTrks = 2;
static constexpr int kRun2SPDClss = 3;
static constexpr int kRun2CL0s = 4;
static constexpr int kRun2CL1s = 5;
static constexpr int kFV0As = 6;
static constexpr int kFT0Ms = 7;
static constexpr int kFT0As = 8;
static constexpr int kFT0Cs = 9;
static constexpr int kFT0CVariant1s = 10;
static constexpr int kFDDMs = 11;
static constexpr int kNTPVs = 12;
static constexpr int kNGlobals = 13;
static constexpr int kMFTs = 14;
static constexpr int kNTables = 15;
static const std::vector<std::string> tableNames{"CentRun2V0Ms", // 0
    "CentRun2V0As", // 1
    "CentRun2SPDTrks", // 2
    "CentRun2SPDClss", // 3
    "CentRun2CL0s", // 4
    "CentRun2CL1s", // 5
    "CentFV0As", // 6 
    "CentFT0Ms", // 7
    "CentFT0As", // 8
    "CentFT0Cs", // 9
    "CentFT0CVariant1s", // 10
    "CentFDDMs", // 11
    "CentNTPVs", // 12
    "CentNGlobals", // 13
    "CentMFTs"};

static const int defaultParameters[kNTables][1]{{-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}};

}

static const std::vector<std::string> parameterNames{"Enable"};




struct multiplicityPercentile {
    SliceCache cache;
    // Services
    Service<o2::ccdb::BasicCCDBManager> ccdb;
    Service<o2::framework::O2DatabasePDG> pdg;



    // Multiplicity tables
    Produces<aod::FV0Mults> tableFV0;             // 0
    Produces<aod::FV0AOuterMults> tableFV0AOuter; // 0-bis (produced with FV0)
    Produces<aod::FT0Mults> tableFT0;             // 1
    Produces<aod::FDDMults> tableFDD;             // 2
    Produces<aod::ZDCMults> tableZDC;             // 3
    Produces<aod::TrackletMults> tableTracklet;   // 4
    Produces<aod::TPCMults> tableTpc;             // 5
    Produces<aod::PVMults> tablePv;               // 6
    Produces<aod::MultsExtra> tableExtra;         // 7
    Produces<aod::MultSelections> multSelections; // 8
    Produces<aod::FV0MultZeqs> tableFV0Zeqs;      // 9
    Produces<aod::FT0MultZeqs> tableFT0Zeqs;      // 10
    Produces<aod::FDDMultZeqs> tableFDDZeqs;      // 11
    Produces<aod::PVMultZeqs> tablePVZeqs;        // 12
    Produces<aod::MultMCExtras> tableExtraMc;     // 13
    Produces<aod::Mult2MCExtras> tableExtraMult2MCExtras;
    Produces<aod::MFTMults> mftMults;             // Not accounted for, produced using custom process function to avoid dependencies
    Produces<aod::MultsGlobal> multsGlobal;       // Not accounted for, produced based on process function processGlobalTrackingCounters
    Configurable<LabeledArray<int>> enabledMultiplicityTables{"enabledMultiplicityTables",
        {defaultParameters[0], multiplicity::kNTables, 1, multiplicity::tableNames, parameterNames},
        "Produce multiplicity tables depending on needs. Values different than -1 override the automatic setup: the corresponding table can be set off (0) or on (1)"};


    // Centrality tables
    Produces<aod::CentRun2V0Ms> centRun2V0M;
    Produces<aod::CentRun2V0As> centRun2V0A;
    Produces<aod::CentRun2SPDTrks> centRun2SPDTracklets;
    Produces<aod::CentRun2SPDClss> centRun2SPDClusters;
    Produces<aod::CentRun2CL0s> centRun2CL0;
    Produces<aod::CentRun2CL1s> centRun2CL1;
    Produces<aod::CentFV0As> centFV0A;
    Produces<aod::CentFT0Ms> centFT0M;
    Produces<aod::CentFT0As> centFT0A;
    Produces<aod::CentFT0Cs> centFT0C;
    Produces<aod::CentFT0CVariant1s> centFT0CVariant1;
    Produces<aod::CentFDDMs> centFDDM;
    Produces<aod::CentNTPVs> centNTPV;
    Produces<aod::CentNGlobals> centNGlobals;
    Produces<aod::CentMFTs> centMFTs;
    Configurable<LabeledArray<int>> enabledCentralityTables{"enabledCentralityTables",
                                                  {defaultParameters[0], centrality::kNTables, 1, centrality::tableNames, parameterNames},
                                                  "Produce centrality tables depending on needs. Values different than -1 override the automatic setup: the corresponding table can be set off (0) or on (1)"};
  

// Configuration
Configurable<bool> produceHistograms{"produceHistograms", false, {"Option to produce debug histograms"}};
Configurable<bool> autoSetupFromMetadata{"autoSetupFromMetadata", true, {"Autosetup the Run 2 and Run 3 processing from the metadata"}};


struct : ConfigurableGroup {
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "The CCDB endpoint url address"};
    Configurable<std::string> ccdbPath{"ccdbPath", "Centrality/Estimators", "The CCDB path for centrality/multiplicity information"};
    Configurable<std::string> genName{"genName", "", "Genearator name: HIJING, PYTHIA8, ... Default: \"\""};
    Configurable<bool> doNotCrashOnNull{"doNotCrashOnNull", false, {"Option to not crash on null and instead fill required tables with dummy info"}};
    Configurable<std::string> reconstructionPass{"reconstructionPass", "", {"Apass to use when fetching the calibration tables. Empty (default) does not check for any pass. Use `metadata` to fetch it from the AO2D metadata. Otherwise it will override the metadata."}};
  } ccdbConfig;


  // Configurable multiplicity
  struct : ConfigurableGroup {
  Configurable<int> doVertexZeq{"doVertexZeq", 1, "if 1: do vertex Z eq mult table"};
  Configurable<float> fractionOfEvents{"fractionOfEvents", 2.0, "Fractions of events to keep in case the QA is used"};                                         
} multiplicityConfig;





  // Configurable multiplicity
  struct : ConfigurableGroup {

  } centralityConfig;

  int mRunNumber;
  bool lCalibLoaded;
  TList* lCalibObjects;
  TProfile* hVtxZFV0A;
  TProfile* hVtxZFT0A;
  TProfile* hVtxZFT0C;
  TProfile* hVtxZFDDA;

  TProfile* hVtxZFDDC;
  TProfile* hVtxZNTracks;
  std::vector<int> mEnabledTables; // Vector of enabled tables

  // Debug output
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};
  OutputObj<TList> listCalib{"calib-list", OutputObjHandlingPolicy::QAObject};

  unsigned int randomSeed = 0;
  void init(InitContext& context)
  {
  }

  using Run2Tracks = soa::Join<aod::Tracks, aod::TracksExtra>;
  Partition<Run2Tracks> run2tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));
  Partition<Run2Tracks> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  Partition<Run2Tracks> pvContribTracks = (nabs(aod::track::eta) < 0.8f) && ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Partition<Run2Tracks> pvContribTracksEta1 = (nabs(aod::track::eta) < 1.0f) && ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::TracksIU> perColIU = aod::track::collisionId;
  Preslice<aod::MFTTracks> perCollisionMFT = o2::aod::fwdtrack::collisionId;

  using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;


}