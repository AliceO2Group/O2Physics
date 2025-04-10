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

namespace multiplicity
{
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
} // namespace multiplicity

namespace centrality
{
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
static const std::vector<std::string> tableNames{"CentRun2V0Ms",      // 0
                                                 "CentRun2V0As",      // 1
                                                 "CentRun2SPDTrks",   // 2
                                                 "CentRun2SPDClss",   // 3
                                                 "CentRun2CL0s",      // 4
                                                 "CentRun2CL1s",      // 5
                                                 "CentFV0As",         // 6
                                                 "CentFT0Ms",         // 7
                                                 "CentFT0As",         // 8
                                                 "CentFT0Cs",         // 9
                                                 "CentFT0CVariant1s", // 10
                                                 "CentFDDMs",         // 11
                                                 "CentNTPVs",         // 12
                                                 "CentNGlobals",      // 13
                                                 "CentMFTs"};

static const int defaultParameters[kNTables][1]{{-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}};

} // namespace centrality

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
  Produces<aod::MFTMults> mftMults;       // Not accounted for, produced using custom process function to avoid dependencies
  Produces<aod::MultsGlobal> multsGlobal; // Not accounted for, produced based on process function processGlobalTrackingCounters
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

  struct TagRun2V0MCalibration {
    bool lCalibLoaded;
    TList* lCalibObjectsMultiplicity;
    TProfile* hVtxZFV0A;
    TProfile* hVtxZFT0A;
    TProfile* hVtxZFT0C;
    TProfile* hVtxZFDDA;
    TProfile* hVtxZFDDC;
    TProfile* hVtxZNTracks;
  } centralityCalib;

  int mRunNumber;
  std::vector<int> mEnabledMultiplicityTables; // Vector of enabled multiplicity tables
  std::vector<int> mEnabledCentralityTables;   // Vector of enabled centrality tables

  struct TagRun2V0MCalibration {
    bool mCalibrationStored = false;
    TFormula* mMCScale = nullptr;
    float mMCScalePars[6] = {0.0};
    TH1* mhVtxAmpCorrV0A = nullptr;
    TH1* mhVtxAmpCorrV0C = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2V0MInfo;
  struct TagRun2V0ACalibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorrV0A = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2V0AInfo;
  struct TagRun2SPDTrackletsCalibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorr = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2SPDTksInfo;
  struct TagRun2SPDClustersCalibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorrCL0 = nullptr;
    TH1* mhVtxAmpCorrCL1 = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2SPDClsInfo;
  struct TagRun2CL0Calibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorr = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2CL0Info;
  struct TagRun2CL1Calibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorr = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2CL1Info;

  struct CentralityCalibration {
    std::string name = "";
    bool mCalibrationStored = false;
    TH1* mhMultSelCalib = nullptr;
    float mMCScalePars[6] = {0.0};
    TFormula* mMCScale = nullptr;
    explicit CentralityCalibration(std::string name)
      : name(name),
        mCalibrationStored(false),
        mhMultSelCalib(nullptr),
        mMCScalePars{0.0},
        mMCScale(nullptr)
    {
    }
    bool isSane(bool fatalize = false)
    {
      if (!mhMultSelCalib) {
        return true;
      }
      for (int i = 1; i < mhMultSelCalib->GetNbinsX() + 1; i++) {
        if (mhMultSelCalib->GetXaxis()->GetBinLowEdge(i) > mhMultSelCalib->GetXaxis()->GetBinUpEdge(i)) {
          if (fatalize) {
            LOG(fatal) << "Centrality calibration table " << name << " has bins with low edge > up edge";
          }
          LOG(warning) << "Centrality calibration table " << name << " has bins with low edge > up edge";
          return false;
        }
      }
      return true;
    }
  };
  CentralityCalibration fv0aInfo = CentralityCalibration("FV0");
  CentralityCalibration ft0mInfo = CentralityCalibration("FT0");
  CentralityCalibration ft0aInfo = CentralityCalibration("FT0A");
  CentralityCalibration ft0cInfo = CentralityCalibration("FT0C");
  CentralityCalibration ft0cVariant1Info = CentralityCalibration("FT0Cvar1");
  CentralityCalibration fddmInfo = CentralityCalibration("FDD");
  CentralityCalibration ntpvInfo = CentralityCalibration("NTracksPV");
  CentralityCalibration nGlobalInfo = CentralityCalibration("NGlobal");
  CentralityCalibration mftInfo = CentralityCalibration("MFT");

  // Debug output
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};
  OutputObj<TList> listCalibMultiplicity{"calib-list", OutputObjHandlingPolicy::QAObject};

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
        mEnabledMultiplicityTables.push_back(i);
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
      mEnabledMultiplicityTables.push_back(kFV0Mults);
      LOG(info) << "Cannot have the " << tableNames[kFV0MultZeqs] << " table enabled and not the one on " << tableNames[kFV0Mults] << ". Enabling it.";
    }
    if (enabledMultiplicities[kFT0MultZeqs] && !enabledMultiplicities[kFT0Mults]) { // FT0
      mEnabledMultiplicityTables.push_back(kFT0Mults);
      LOG(info) << "Cannot have the " << tableNames[kFT0MultZeqs] << " table enabled and not the one on " << tableNames[kFT0Mults] << ". Enabling it.";
    }
    if (enabledMultiplicities[kFDDMultZeqs] && !enabledMultiplicities[kFDDMults]) { // FDD
      mEnabledMultiplicityTables.push_back(kFDDMults);
      LOG(info) << "Cannot have the " << tableNames[kFDDMultZeqs] << " table enabled and not the one on " << tableNames[kFDDMults] << ". Enabling it.";
    }
    if (enabledMultiplicities[kPVMultZeqs] && !enabledMultiplicities[kPVMults]) { // PV
      mEnabledMultiplicityTables.push_back(kPVMults);
      LOG(info) << "Cannot have the " << tableNames[kPVMultZeqs] << " table enabled and not the one on " << tableNames[kPVMults] << ". Enabling it.";
    }
    std::sort(mEnabledMultiplicityTables.begin(), mEnabledMultiplicityTables.end());

    /* Checking the tables which are requested in the workflow and enabling them */
    for (int i = 0; i < centrality::kNTables; i++) {
      int f = enabledCentralityTables->get(centrality::tableNames[i].c_str(), "Enable");
      enableFlagIfTableRequired(context, centrality::tableNames[i], f);
      if (f == 1) {
        if (centrality::tableNames[i].find("Run2") != std::string::npos) {
          if (doprocessRun3) {
            LOG(fatal) << "Cannot enable Run2 table `" << centrality::tableNames[i] << "` while running in Run3 mode. Please check and disable them.";
          }
        } else {
          if (doprocessRun2) {
            LOG(fatal) << "Cannot enable Run3 table `" << centrality::tableNames[i] << "` while running in Run2 mode. Please check and disable them.";
          }
        }
        isTableEnabled[i] = true;
        mEnabledCentralityTables.push_back(i);
      }
    }

    if (mEnabledCentralityTables.size() == 0) {
      LOGF(fatal, "No table enabled. Please enable at least one table.");
    }
    std::sort(mEnabledCentralityTables.begin(), mEnabledCentralityTables.end());

    mRunNumber = 0;
    lCalibLoaded = false;
    lCalibObjectsMultiplicity = nullptr;
    hVtxZFV0A = nullptr;
    hVtxZFT0A = nullptr;
    hVtxZFT0C = nullptr;
    hVtxZFDDA = nullptr;
    hVtxZFDDC = nullptr;
    hVtxZNTracks = nullptr;
    listCalibMultiplicity.setObject(new TList);

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
    for (const auto tableId : mEnabledMultiplicityTables) {
      switch (tableId) {
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
          LOG(fatal) << "Unknown table requested: " << tableId;
          break;
      }
    }

    // reserve memory for centrality tables
    for (auto const tableId : mEnabledCentralityTables) {
      switch (tableId) {
        case centrality::kFV0As:
          centFV0A.reserve(collisions.size());
          break;
        case centrality::kFT0Ms:
          centFT0M.reserve(collisions.size());
          break;
        case centrality::kFT0As:
          centFT0A.reserve(collisions.size());
          break;
        case centrality::kFT0Cs:
          centFT0C.reserve(collisions.size());
          break;
        case centrality::kFT0CVariant1s:
          centFT0CVariant1.reserve(collisions.size());
          break;
        case centrality::kFDDMs:
          centFDDM.reserve(collisions.size());
          break;
        case centrality::kNTPVs:
          centNTPV.reserve(collisions.size());
          break;
        case centrality::kNGlobals:
          centNGlobals.reserve(collisions.size());
          break;
        case centrality::kMFTs:
          centMFTs.reserve(collisions.size());
          break;
        default:
          LOGF(fatal, "Table %d not supported in Run3", tableId);
          break;
      }
    }

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
            lCalibObjectsMultiplicity = ccdb->getForRun<TList>(ccdbConfig.ccdbPath, mRunNumber);
          } else if (ccdbConfig.reconstructionPass.value == "metadata") {
            std::map<std::string, std::string> metadata;
            metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
            LOGF(info, "Loading CCDB for reconstruction pass (from metadata): %s", metadataInfo.get("RecoPassName"));
            lCalibObjectsMultiplicity = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, mRunNumber, metadata);
          } else {
            std::map<std::string, std::string> metadata;
            metadata["RecoPassName"] = ccdbConfig.reconstructionPass.value;
            LOGF(info, "Loading CCDB for reconstruction pass (from provided argument): %s", ccdbConfig.reconstructionPass.value);
            lCalibObjectsMultiplicity = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, mRunNumber, metadata);
          }

          if (lCalibObjectsMultiplicity) {
            if (produceHistograms) {
              listCalibMultiplicity->Add(lCalibObjectsMultiplicity->Clone(Form("%i", bc.runNumber())));
            }

            hVtxZFV0A = static_cast<TProfile*>(lCalibObjectsMultiplicity->FindObject("hVtxZFV0A"));
            hVtxZFT0A = static_cast<TProfile*>(lCalibObjectsMultiplicity->FindObject("hVtxZFT0A"));
            hVtxZFT0C = static_cast<TProfile*>(lCalibObjectsMultiplicity->FindObject("hVtxZFT0C"));
            hVtxZFDDA = static_cast<TProfile*>(lCalibObjectsMultiplicity->FindObject("hVtxZFDDA"));
            hVtxZFDDC = static_cast<TProfile*>(lCalibObjectsMultiplicity->FindObject("hVtxZFDDC"));
            hVtxZNTracks = static_cast<TProfile*>(lCalibObjectsMultiplicity->FindObject("hVtxZNTracksPV"));
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

      // First we compute the multiplicity
      for (const auto tableId : mEnabledMultiplicityTables) {
        switch (tableId) {
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
            LOG(fatal) << "Unknown table requested: " << tableId;
          } break;
        }
      }

      if (mEnabledCentralityTables.size() == 0) { // If no centrality table is required skip the rest
        continue;
      }

      if (bc.runNumber() != mRunNumber) {
        mRunNumber = bc.runNumber(); // mark that this run has been attempted already regardless of outcome
        LOGF(info, "timestamp=%llu, run number=%d", bc.timestamp(), bc.runNumber());
        TList* calibrationList = nullptr;
        // Check if the ccdb path is a root file
        if (ccdbConfig.ccdbPath.value.find(".root") != std::string::npos) { // File
          LOG(info) << "Fetching centrality calibration from TFile" << ccdbConfig.ccdbPath.value.c_str() << " and pass '" << ccdbConfig.reconstructionPass.value << "'";
          TFile f(ccdbConfig.ccdbPath.value.c_str(), "READ");
          f.GetObject(ccdbConfig.reconstructionPass.value.c_str(), calibrationList);
          if (!calibrationList) {
            f.ls();
            LOG(fatal) << "No calibration list " << ccdbConfig.reconstructionPass.value << " found in the file " << ccdbConfig.ccdbPath.value;
          }
        } else { // CCDB
          LOG(info) << "Fetching centrality calibration from ccdb" << ccdbConfig.ccdbPath.value << " and pass '" << ccdbConfig.reconstructionPass.value << "'";
          std::map<std::string, std::string> metadata;
          if (ccdbConfig.reconstructionPass.value == "") {
            LOG(info) << "No pass required";
          } else {
            if (ccdbConfig.reconstructionPass.value == "metadata") {
              LOGF(info, "Loading CCDB for reconstruction pass (from metadata): %s", metadataInfo.get("RecoPassName"));
              ccdbConfig.reconstructionPass.value = metadataInfo.get("RecoPassName");
            }
            metadata["RecoPassName"] = ccdbConfig.reconstructionPass.value;
          }
          calibrationList = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, bc.runNumber(), metadata);
        }

        if (calibrationList != nullptr) {
          if (produceHistograms) {
            listCalib->Add(calibrationList->Clone(Form("%i", bc.runNumber())));
          }
          LOGF(info, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          auto getccdb = [calibrationList, bc](struct CentralityCalibration& estimator, const Configurable<std::string> generatorName) { // TODO: to consider the name inside the estimator structure
            estimator.mhMultSelCalib = reinterpret_cast<TH1*>(calibrationList->FindObject(TString::Format("hCalibZeq%s", estimator.name.c_str()).Data()));
            estimator.mMCScale = reinterpret_cast<TFormula*>(calibrationList->FindObject(TString::Format("%s-%s", generatorName->c_str(), estimator.name.c_str()).Data()));
            if (estimator.mhMultSelCalib != nullptr) {
              if (generatorName->length() != 0) {
                LOGF(info, "Retrieving MC calibration for %d, generator name: %s", bc.runNumber(), generatorName->c_str());
                if (estimator.mMCScale != nullptr) {
                  for (int ixpar = 0; ixpar < 6; ++ixpar) {
                    estimator.mMCScalePars[ixpar] = estimator.mMCScale->GetParameter(ixpar);
                    LOGF(info, "Parameter index %i value %.5f", ixpar, estimator.mMCScalePars[ixpar]);
                  }
                } else {
                  LOGF(warning, "MC Scale information from %s for run %d not available", estimator.name.c_str(), bc.runNumber());
                }
              }
              estimator.mCalibrationStored = true;
              estimator.isSane();
            } else {
              LOGF(info, "Calibration information from %s for run %d not available, will fill this estimator with invalid values and continue (no crash).", estimator.name.c_str(), bc.runNumber());
            }
          };

          for (auto const& table : mEnabledCentralityTables) {
            switch (table) {
              case centrality::kFV0As:
                getccdb(fv0aInfo, ccdbConfig.genName);
                break;
              case centrality::kFT0Ms:
                getccdb(ft0mInfo, ccdbConfig.genName);
                break;
              case centrality::kFT0As:
                getccdb(ft0aInfo, ccdbConfig.genName);
                break;
              case centrality::kFT0Cs:
                getccdb(ft0cInfo, ccdbConfig.genName);
                break;
              case centrality::kFT0CVariant1s:
                getccdb(ft0cVariant1Info, ccdbConfig.genName);
                break;
              case centrality::kFDDMs:
                getccdb(fddmInfo, ccdbConfig.genName);
                break;
              case centrality::kNTPVs:
                getccdb(ntpvInfo, ccdbConfig.genName);
                break;
              case centrality::kNGlobals:
                getccdb(nGlobalInfo, ccdbConfig.genName);
                break;
              case centrality::kMFTs:
                getccdb(mftInfo, ccdbConfig.genName);
                break;
              default:
                LOGF(fatal, "Table %d not supported in Run3", table);
                break;
            }
          }
        } else {
          if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
            LOGF(fatal, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
          } else { // only if asked: continue filling with non-valid values (105)
            LOGF(info, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu, will fill tables with dummy values", bc.runNumber(), bc.timestamp());
          }
        }
      }

      /**
       * @brief Populates a table with data based on the given calibration information and multiplicity.
       *
       * @param table The table to populate.
       * @param estimator The calibration information.
       * @param multiplicity The multiplicity value.
       */

      auto populateCentralityTable = [&](auto& table,
                                         struct CentralityCalibration& estimator,
                                         float multiplicity) {
        const bool assignOutOfRange = embedINELgtZEROselection && !collision.isInelGt0();
        auto scaleMC = [](float x, float pars[6]) {
          return std::pow(((pars[0] + pars[1] * std::pow(x, pars[2])) - pars[3]) / pars[4], 1.0f / pars[5]);
        };

        float percentile = 105.0f;
        float scaledMultiplicity = multiplicity;
        if (estimator.mCalibrationStored) {
          if (estimator.mMCScale != nullptr) {
            scaledMultiplicity = scaleMC(multiplicity, estimator.mMCScalePars);
            LOGF(debug, "Unscaled %s multiplicity: %f, scaled %s multiplicity: %f", estimator.name.c_str(), multiplicity, estimator.name.c_str(), scaledMultiplicity);
          }
          percentile = estimator.mhMultSelCalib->GetBinContent(estimator.mhMultSelCalib->FindFixBin(scaledMultiplicity));
          if (assignOutOfRange)
            percentile = 100.5f;
        }
        LOGF(debug, "%s centrality/multiplicity percentile = %.0f for a zvtx eq %s value %.0f", estimator.name.c_str(), percentile, estimator.name.c_str(), scaledMultiplicity);
        table(percentile);
        return percentile;
      };

      for (auto const& table : mEnabledCentralityTables) {
        switch (table) {
          case centrality::kFV0As:
            populateCentralityTable(centFV0A, fv0aInfo, multZeqFV0A);
            break;
          case centrality::kFT0Ms:
            const float perC = populateCentralityTable(centFT0M, ft0mInfo, multZeqFT0A + multZeqFT0C);
            break;
          case centrality::kFT0As:
            const float perC = populateCentralityTable(centFT0A, ft0aInfo, multZeqFT0A);
            break;
          case centrality::kFT0Cs:
            const float perC = populateCentralityTable(centFT0C, ft0cInfo, multZeqFT0C);
            break;
          case centrality::kFT0CVariant1s:
            populateCentralityTable(centFT0CVariant1, ft0cVariant1Info, multZeqFT0C);
            break;
          case centrality::kFDDMs:
            populateCentralityTable(centFDDM, fddmInfo, multZeqFDDA + multZeqFDDC);
            break;
          case centrality::kNTPVs:
            populateCentralityTable(centNTPV, ntpvInfo, multZeqNTracksPV);
            break;
          case centrality::kNGlobals:
            populateCentralityTable(centNGlobals, nGlobalInfo, multNTracksGlobal);
            break;
          case centrality::kMFTs:
            populateCentralityTable(centMFTs, mftInfo, mftNtracks);
            break;
          default:
            LOGF(fatal, "Table %d not supported in Run3", table);
            break;
        }
      }
    }
  }
  PROCESS_SWITCH(CentralityTable, processRun3, "Provide Run3 calibrated centrality/multiplicity percentiles tables", true);
}

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<multPercentile>(cfgc)};
}
