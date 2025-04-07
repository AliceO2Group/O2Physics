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
    ccdb->setURL(ccdbConfig.ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false); // don't fatal, please - exception is caught explicitly (as it should)
    


      if (autoSetupFromMetadata && metadataInfo.isFullyDefined()) {
        LOG(info) << "Autosetting the processing from the metadata";
          if (metadataInfo.isRun3()) {
            doprocessRun2.value = false;
          } else {
            doprocessRun3.value = false;
          }
      }
  
      randomSeed = static_cast<unsigned int>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      if (doprocessRun2 == false && doprocessRun3 == false) {
        LOGF(fatal, "Neither processRun2 nor processRun3 enabled. Please choose one.");
      }
      if (doprocessRun2 == true && doprocessRun3 == true) {
        LOGF(fatal, "Cannot enable processRun2 and processRun3 at the same time. Please choose one.");
      }
  
      // Enable or disable the multiplicity tables
      bool enabledMultiplicities[multiplicity::kNTables] = {false};
      for (int i = 0; i < multiplicity::kNTables; i++) {
        int f = enabledTables->get(tableNames[i].c_str(), "Enable");
        enableFlagIfTableRequired(context, tableNames[i], f);
        if (f == 1) {
          enabledMultiplicities[i] = true;
          mEnabledTables.push_back(i);
          if (fractionOfEvents <= 1.f && (tableNames[i] != "MultsExtra")) {
            LOG(fatal) << "Cannot have a fraction of events <= 1 and multiplicity table consumed.";
          }
        }
      }
      // Handle the custom cases.
      if (enabledMultiplicities[kMultMCExtras]) {
        if (enabledTables->get(tableNames[kMultMCExtras].c_str(), "Enable") == -1) {
          doprocessMC.value = true;
          LOG(info) << "Enabling MC processing due to " << tableNames[kMultMCExtras] << " table being enabled.";
        }
      }
  
      // Check that the tables are enabled consistenly
      if (enabledMultiplicities[kFV0MultZeqs] && !enabledMultiplicities[kFV0Mults]) { // FV0
        mEnabledTables.push_back(kFV0Mults);
        LOG(info) << "Cannot have the " << tableNames[kFV0MultZeqs] << " table enabled and not the one on " << tableNames[kFV0Mults] << ". Enabling it.";
      }
      if (enabledMultiplicities[kFT0MultZeqs] && !enabledMultiplicities[kFT0Mults]) { // FT0
        mEnabledTables.push_back(kFT0Mults);
        LOG(info) << "Cannot have the " << tableNames[kFT0MultZeqs] << " table enabled and not the one on " << tableNames[kFT0Mults] << ". Enabling it.";
      }
      if (enabledMultiplicities[kFDDMultZeqs] && !enabledMultiplicities[kFDDMults]) { // FDD
        mEnabledTables.push_back(kFDDMults);
        LOG(info) << "Cannot have the " << tableNames[kFDDMultZeqs] << " table enabled and not the one on " << tableNames[kFDDMults] << ". Enabling it.";
      }
      if (enabledMultiplicities[kPVMultZeqs] && !enabledMultiplicities[kPVMults]) { // PV
        mEnabledTables.push_back(kPVMults);
        LOG(info) << "Cannot have the " << tableNames[kPVMultZeqs] << " table enabled and not the one on " << tableNames[kPVMults] << ". Enabling it.";
      }
      std::sort(mEnabledTables.begin(), mEnabledTables.end());
  
      mRunNumber = 0;
      lCalibLoaded = false;
      lCalibObjects = nullptr;
      hVtxZFV0A = nullptr;
      hVtxZFT0A = nullptr;
      hVtxZFT0C = nullptr;
      hVtxZFDDA = nullptr;
      hVtxZFDDC = nullptr;
      hVtxZNTracks = nullptr;
      listCalib.setObject(new TList);

      if (!produceHistograms.value) {
        return;
      }
      histos.add("FT0A", "FT0A vs FT0A eq.", HistType::kTH2D, {{1000, 0, 1000, "FT0A multiplicity"}, {1000, 0, 1000, "FT0A multiplicity eq."}});
      histos.add("FT0C", "FT0C vs FT0C eq.", HistType::kTH2D, {{1000, 0, 1000, "FT0C multiplicity"}, {1000, 0, 1000, "FT0C multiplicity eq."}});
      histos.add("FT0CMultvsPV", "FT0C vs mult.", HistType::kTH2D, {{1000, 0, 1000, "FT0C mult."}, {100, 0, 100, "PV mult."}});
      histos.add("FT0AMultvsPV", "FT0A vs mult.", HistType::kTH2D, {{1000, 0, 1000, "FT0A mult."}, {100, 0, 100, "PV mult."}});
    }
  
    /// Dummy process function for BCs, needed in case both Run2 and Run3 process functions are disabled
    void process(aod::BCs const&) {}

  using Run2Tracks = soa::Join<aod::Tracks, aod::TracksExtra>;
  Partition<Run2Tracks> run2tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));
  Partition<Run2Tracks> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  Partition<Run2Tracks> pvContribTracks = (nabs(aod::track::eta) < 0.8f) && ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Partition<Run2Tracks> pvContribTracksEta1 = (nabs(aod::track::eta) < 1.0f) && ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::TracksIU> perColIU = aod::track::collisionId;
  Preslice<aod::MFTTracks> perCollisionMFT = o2::aod::fwdtrack::collisionId;

  using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;



  using Run3TracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;
  Partition<Run3TracksIU> tracksIUWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  Partition<Run3TracksIU> pvAllContribTracksIU = ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Partition<Run3TracksIU> pvContribTracksIU = (nabs(aod::track::eta) < 0.8f) && ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Partition<Run3TracksIU> pvContribTracksIUEta1 = (nabs(aod::track::eta) < 1.0f) && ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Partition<Run3TracksIU> pvContribTracksIUEtaHalf = (nabs(aod::track::eta) < 0.5f) && ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));

  void processRun3(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                   Run3TracksIU const&,
                   BCsWithRun3Matchings const&,
                   aod::Zdcs const&,
                   aod::FV0As const&,
                   aod::FT0s const&,
                   aod::FDDs const&)
  {
    // reserve memory for multiplicity tables
    for (const auto& i : mEnabledTables) {
      switch (i) {
        case kFV0Mults: // FV0
          tableFV0.reserve(collisions.size());
          tableFV0AOuter.reserve(collisions.size());
          break;
        case kFT0Mults: // FT0
          tableFT0.reserve(collisions.size());
          break;
        case kFDDMults: // FDD
          tableFDD.reserve(collisions.size());
          break;
        case kZDCMults: // ZDC
          tableZDC.reserve(collisions.size());
          break;
        case kTrackletMults: // Tracklets (Run 2 only, nothing to do) (to be removed!)
          tableTracklet.reserve(collisions.size());
          break;
        case kTPCMults: // TPC
          tableTpc.reserve(collisions.size());
          break;
        case kPVMults: // PV multiplicity
          tablePv.reserve(collisions.size());
          break;
        case kMultsExtra: // Extra information
          tableExtra.reserve(collisions.size());
          break;
        case kMultSelections: // Extra information
          multSelections.reserve(collisions.size());
          break;
        case kFV0MultZeqs: // Equalized multiplicity for FV0
          tableFV0Zeqs.reserve(collisions.size());
          break;
        case kFT0MultZeqs: // Equalized multiplicity for FT0
          tableFT0Zeqs.reserve(collisions.size());
          break;
        case kFDDMultZeqs: // Equalized multiplicity for FDD
          tableFDDZeqs.reserve(collisions.size());
          break;
        case kPVMultZeqs: // Equalized multiplicity for PV
          tablePVZeqs.reserve(collisions.size());
          break;
        case kMultMCExtras: // MC extra information (nothing to do in the data)
          break;
        default:
          LOG(fatal) << "Unknown table requested: " << i;
          break;
      }
    }

    // reserve memory for centrality tables

    // Initializing multiplicity values
    float multFV0A = 0.f;
    float multFV0AOuter = 0.f;
    float multFV0C = 0.f;
    float multFT0A = 0.f;
    float multFT0C = 0.f;
    float multFDDA = 0.f;
    float multFDDC = 0.f;
    float multZNA = -1.f;
    float multZNC = -1.f;
    float multZEM1 = -1.f;
    float multZEM2 = -1.f;
    float multZPA = -1.f;
    float multZPC = -1.f;

    float multZeqFV0A = 0.f;
    float multZeqFT0A = 0.f;
    float multZeqFT0C = 0.f;
    float multZeqFDDA = 0.f;
    float multZeqFDDC = 0.f;
    float multZeqNContribs = 0.f;

    for (auto const& collision : collisions) {
      if ((fractionOfEvents < 1.f) && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > fractionOfEvents) { // Skip events that are not sampled (only for the QA)
        return;
      }
      int multNContribs = 0;
      int multNContribsEta1 = 0;
      int multNContribsEtaHalf = 0;

      /* check the previous run number */
      const auto& bc = collision.bc_as<BCsWithRun3Matchings>();
      if (doVertexZeq > 0) {
        if (bc.runNumber() != mRunNumber) {
          mRunNumber = bc.runNumber(); // mark this run as at least tried
          if (ccdbConfig.reconstructionPass.value == "") {
            lCalibObjects = ccdb->getForRun<TList>(ccdbConfig.ccdbPath, mRunNumber);
          } else if (ccdbConfig.reconstructionPass.value == "metadata") {
            std::map<std::string, std::string> metadata;
            metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
            LOGF(info, "Loading CCDB for reconstruction pass (from metadata): %s", metadataInfo.get("RecoPassName"));
            lCalibObjects = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, mRunNumber, metadata);
          } else {
            std::map<std::string, std::string> metadata;
            metadata["RecoPassName"] = ccdbConfig.reconstructionPass.value;
            LOGF(info, "Loading CCDB for reconstruction pass (from provided argument): %s", ccdbConfig.reconstructionPass.value);
            lCalibObjects = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, mRunNumber, metadata);
          }

          if (lCalibObjects) {
            if (produceHistograms) {
              listCalib->Add(lCalibObjects->Clone(Form("%i", bc.runNumber())));
            }

            hVtxZFV0A = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFV0A"));
            hVtxZFT0A = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFT0A"));
            hVtxZFT0C = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFT0C"));
            hVtxZFDDA = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFDDA"));
            hVtxZFDDC = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFDDC"));
            hVtxZNTracks = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZNTracksPV"));
            lCalibLoaded = true;
            // Capture error
            if (!hVtxZFV0A || !hVtxZFT0A || !hVtxZFT0C || !hVtxZFDDA || !hVtxZFDDC || !hVtxZNTracks) {
              LOGF(error, "Problem loading CCDB objects! Please check");
              lCalibLoaded = false;
            }
          } else {
            LOGF(error, "Problem loading CCDB object! Please check");
            lCalibLoaded = false;
          }
        }
      }

      for (const auto& i : mEnabledTables) {
        switch (i) {
          case kFV0Mults: // FV0
          {
            multFV0A = 0.f;
            multFV0AOuter = 0.f;
            multFV0C = 0.f;
            // using FV0 row index from event selection task
            if (collision.has_foundFV0()) {
              const auto& fv0 = collision.foundFV0();
              for (size_t ii = 0; ii < fv0.amplitude().size(); ii++) {
                auto amplitude = fv0.amplitude()[ii];
                auto channel = fv0.channel()[ii];
                multFV0A += amplitude;
                if (channel > 7) {
                  multFV0AOuter += amplitude;
                }
              }
            } else {
              multFV0A = -999.f;
              multFV0C = -999.f;
            }
            tableFV0(multFV0A, multFV0C);
            tableFV0AOuter(multFV0AOuter);
            LOGF(debug, "multFV0A=%5.0f multFV0C=%5.0f", multFV0A, multFV0C);
          } break;
          case kFT0Mults: // FT0
          {
            multFT0A = 0.f;
            multFT0C = 0.f;
            // using FT0 row index from event selection task
            if (collision.has_foundFT0()) {
              const auto& ft0 = collision.foundFT0();
              for (const auto& amplitude : ft0.amplitudeA()) {
                multFT0A += amplitude;
              }
              for (const auto& amplitude : ft0.amplitudeC()) {
                multFT0C += amplitude;
              }
            } else {
              multFT0A = -999.f;
              multFT0C = -999.f;
            }
            tableFT0(multFT0A, multFT0C);
            LOGF(debug, "multFT0A=%5.0f multFT0C=%5.0f", multFV0A, multFV0C);
          } break;
          case kFDDMults: // FDD
          {
            multFDDA = 0.f;
            multFDDC = 0.f;
            // using FDD row index from event selection task
            if (collision.has_foundFDD()) {
              const auto& fdd = collision.foundFDD();
              for (const auto& amplitude : fdd.chargeA()) {
                multFDDA += amplitude;
              }
              for (const auto& amplitude : fdd.chargeC()) {
                multFDDC += amplitude;
              }
            } else {
              multFDDA = -999.f;
              multFDDC = -999.f;
            }
            tableFDD(multFDDA, multFDDC);
            LOGF(debug, "multFDDA=%5.0f multFDDC=%5.0f", multFDDA, multFDDC);
          } break;
          case kZDCMults: // ZDC
          {
            multZNA = -1.f;
            multZNC = -1.f;
            multZEM1 = -1.f;
            multZEM2 = -1.f;
            multZPA = -1.f;
            multZPC = -1.f;
            if (bc.has_zdc()) {
              multZNA = bc.zdc().amplitudeZNA();
              multZNC = bc.zdc().amplitudeZNC();
              multZEM1 = bc.zdc().amplitudeZEM1();
              multZEM2 = bc.zdc().amplitudeZEM2();
              multZPA = bc.zdc().amplitudeZPA();
              multZPC = bc.zdc().amplitudeZPC();
            } else {
              multZNA = -999.f;
              multZNC = -999.f;
              multZEM1 = -999.f;
              multZEM2 = -999.f;
              multZPA = -999.f;
              multZPC = -999.f;
            }
            tableZDC(multZNA, multZNC, multZEM1, multZEM2, multZPA, multZPC);
            LOGF(debug, "multZNA=%6.0f multZNC=%6.0f", multZNA, multZNC);
          } break;
          case kTrackletMults: // Tracklets (only Run2) nothing to do (to be removed!)
          {
            tableTracklet(0);
          } break;
          case kTPCMults: // TPC
          {
            const auto& tracksGrouped = tracksIUWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
            const int multTPC = tracksGrouped.size();
            tableTpc(multTPC);
            LOGF(debug, "multTPC=%i", multTPC);
          } break;
          case kPVMults: // PV multiplicity
          {
            // use only one single grouping operation, then do loop
            const auto& tracksThisCollision = pvContribTracksIUEta1.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
            multNContribsEta1 = tracksThisCollision.size();
            for (const auto& track : tracksThisCollision) {
              if (std::abs(track.eta()) < 0.8) {
                multNContribs++;
              }
              if (std::abs(track.eta()) < 0.5) {
                multNContribsEtaHalf++;
              }
            }

            tablePv(multNContribs, multNContribsEta1, multNContribsEtaHalf);
            LOGF(debug, "multNContribs=%i, multNContribsEta1=%i, multNContribsEtaHalf=%i", multNContribs, multNContribsEta1, multNContribsEtaHalf);
          } break;
          case kMultsExtra: // Extra
          {
            int nHasITS = 0, nHasTPC = 0, nHasTOF = 0, nHasTRD = 0;
            int nITSonly = 0, nTPConly = 0, nITSTPC = 0;
            const auto& pvAllContribsGrouped = pvAllContribTracksIU->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
            const auto& tpcTracksGrouped = tracksIUWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

            for (const auto& track : pvAllContribsGrouped) {
              if (track.hasITS()) {
                nHasITS++;
                if (track.hasTPC())
                  nITSTPC++;
                if (!track.hasTPC() && !track.hasTOF() && !track.hasTRD())
                  nITSonly++;
              }
              if (track.hasTPC()) {
                nHasTPC++;
                if (!track.hasITS() && !track.hasTOF() && !track.hasTRD())
                  nTPConly++;
              }
              if (track.hasTOF())
                nHasTOF++;
              if (track.hasTRD())
                nHasTRD++;
            }

            int nAllTracksTPCOnly = 0;
            int nAllTracksITSTPC = 0;
            for (const auto& track : tpcTracksGrouped) {
              if (track.hasITS()) {
                nAllTracksITSTPC++;
              } else {
                nAllTracksTPCOnly++;
              }
            }

            tableExtra(collision.numContrib(), collision.chi2(), collision.collisionTimeRes(),
                       mRunNumber, collision.posZ(), collision.sel8(),
                       nHasITS, nHasTPC, nHasTOF, nHasTRD, nITSonly, nTPConly, nITSTPC,
                       nAllTracksTPCOnly, nAllTracksITSTPC,
                       collision.trackOccupancyInTimeRange(),
                       collision.ft0cOccupancyInTimeRange(),
                       collision.flags());
          } break;
          case kMultSelections: // Multiplicity selections
          {
            multSelections(collision.selection_raw());
          } break;
          case kFV0MultZeqs: // Z equalized FV0
          {
            if (std::fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
              multZeqFV0A = hVtxZFV0A->Interpolate(0.0) * multFV0A / hVtxZFV0A->Interpolate(collision.posZ());
            }
            tableFV0Zeqs(multZeqFV0A);
          } break;
          case kFT0MultZeqs: // Z equalized FT0
          {
            if (std::fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
              multZeqFT0A = hVtxZFT0A->Interpolate(0.0) * multFT0A / hVtxZFT0A->Interpolate(collision.posZ());
              multZeqFT0C = hVtxZFT0C->Interpolate(0.0) * multFT0C / hVtxZFT0C->Interpolate(collision.posZ());
            }
            if (produceHistograms.value) {
              histos.fill(HIST("FT0A"), multFT0A, multZeqFT0A);
              histos.fill(HIST("FT0C"), multFT0C, multZeqFT0C);
              histos.fill(HIST("FT0AMultvsPV"), multZeqFT0A, multNContribs);
              histos.fill(HIST("FT0CMultvsPV"), multZeqFT0C, multNContribs);
            }
            tableFT0Zeqs(multZeqFT0A, multZeqFT0C);
          } break;
          case kFDDMultZeqs: // Z equalized FDD
          {
            if (std::fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
              multZeqFDDA = hVtxZFDDA->Interpolate(0.0) * multFDDA / hVtxZFDDA->Interpolate(collision.posZ());
              multZeqFDDC = hVtxZFDDC->Interpolate(0.0) * multFDDC / hVtxZFDDC->Interpolate(collision.posZ());
            }
            tableFDDZeqs(multZeqFDDA, multZeqFDDC);
          } break;
          case kPVMultZeqs: // Z equalized PV
          {
            if (std::fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
              multZeqNContribs = hVtxZNTracks->Interpolate(0.0) * multNContribs / hVtxZNTracks->Interpolate(collision.posZ());
            }
            tablePVZeqs(multZeqNContribs);
          } break;
          case kMultMCExtras: // MC only (nothing to do)
          {
          } break;
          default: // Default
          {
            LOG(fatal) << "Unknown table requested: " << i;
          } break;
        }
      }
    }
  }


}