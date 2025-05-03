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
/// \file multiplicityTable.cxx
/// \brief Produces multiplicity tables
///
/// \author ALICE
///

#include <string>
#include <vector>
#include <algorithm>
#include <map>

#include <TH1F.h>
#include <TFormula.h>
#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "TableHelper.h"
#include "MetadataHelper.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"
#include "TableHelper.h"
#include "TList.h"
#include "TFormula.h"

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
static const int defaultParameters[kNTables][1]{{-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}};
} // namespace multiplicity

static const std::vector<std::string> parameterNames{"Enable"};

struct MultiplicityTable {
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
                                                            {multiplicity::defaultParameters[0], multiplicity::kNTables, 1, multiplicity::tableNames, parameterNames},
                                                            "Produce multiplicity tables depending on needs. Values different than -1 override the automatic setup: the corresponding table can be set off (0) or on (1)"};

  // Configuration
  Configurable<bool> produceHistograms{"produceHistograms", false, {"Option to produce debug histograms"}};
  Configurable<bool> autoSetupFromMetadata{"autoSetupFromMetadata", true, {"Autosetup the Run 2 and Run 3 processing from the metadata"}};

  struct : ConfigurableGroup {
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "The CCDB endpoint url address"};
    Configurable<bool> doNotCrashOnNull{"doNotCrashOnNull", false, {"Option to not crash on null and instead fill required tables with dummy info"}};
    Configurable<std::string> reconstructionPass{"reconstructionPass", "", {"Apass to use when fetching the calibration tables. Empty (default) does not check for any pass. Use `metadata` to fetch it from the AO2D metadata. Otherwise it will override the metadata."}};
  } ccdbConfig;

  // Configurable multiplicity
  struct : ConfigurableGroup {
    Configurable<std::string> ccdbPath{"ccdbPath", "Centrality/Calibration", "The CCDB path for centrality/multiplicity information"};
    Configurable<bool> doVertexZeq{"doVertexZeq", 1, "if 1: do vertex Z eq mult table"};
    Configurable<float> fractionOfEvents{"fractionOfEvents", 2.0, "Fractions of events to keep in case the QA is used"};
    Configurable<int> minNclsITSGlobalTrack{"minNclsITSGlobalTrack", 5, "min. number of ITS clusters for global tracks"};
    Configurable<int> minNclsITSibGlobalTrack{"minNclsITSibGlobalTrack", 1, "min. number of ITSib clusters for global tracks"};
  } multiplicityConfig;
  Configurable<float> minPtGlobalTrack{"minPtGlobalTrack", 0.15, "min. pT for global tracks"};
  Configurable<float> maxPtGlobalTrack{"maxPtGlobalTrack", 1e+10, "max. pT for global tracks"};

  struct MultiplicityCalibration {
    bool isCalibLoaded = false;
    TList* mCalibObjectsMultiplicity = nullptr;
    TProfile* hVtxZFV0A = nullptr;
    TProfile* hVtxZFT0A = nullptr;
    TProfile* hVtxZFT0C = nullptr;
    TProfile* hVtxZFDDA = nullptr;
    TProfile* hVtxZFDDC = nullptr;
    TProfile* hVtxZNTracks = nullptr;
    int mRunNumber = 0;
  } multCalib;

  template <typename T>
  void initMultCalib(const T& bc)
  {
    if (multiplicityConfig.doVertexZeq == false) {
      return;
    }
    if (bc.runNumber() == multCalib.mRunNumber) {
      return;
    }
    multCalib.mRunNumber = bc.runNumber(); // mark this run as at least tried
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
    multCalib.mCalibObjectsMultiplicity = ccdb->getSpecificForRun<TList>(multiplicityConfig.ccdbPath, bc.runNumber(), metadata);
    if (multCalib.mCalibObjectsMultiplicity) {
      if (produceHistograms) {
        listCalib->Add(multCalib.mCalibObjectsMultiplicity->Clone(Form("multCalib_run%i", bc.runNumber())));
      }

      auto findProfile = [&](const char* name) {
        TProfile* p = static_cast<TProfile*>(multCalib.mCalibObjectsMultiplicity->FindObject(name));
        if (!p) {
          multCalib.mCalibObjectsMultiplicity->ls();
          LOGF(error, "CCDB object %s not found in list!", name);
          multCalib.isCalibLoaded = false;
          return static_cast<TProfile*>(nullptr);
        }
        return p;
      };
      multCalib.isCalibLoaded = true;
      multCalib.hVtxZFV0A = findProfile("hVtxZFV0A");
      multCalib.hVtxZFT0A = findProfile("hVtxZFT0A");
      multCalib.hVtxZFT0C = findProfile("hVtxZFT0C");
      multCalib.hVtxZFDDA = findProfile("hVtxZFDDA");
      multCalib.hVtxZFDDC = findProfile("hVtxZFDDC");
      multCalib.hVtxZNTracks = findProfile("hVtxZNTracksPV");
    } else {
      LOGF(error, "Problem loading CCDB object! Please check");
      multCalib.isCalibLoaded = false;
    }
  }

  // List of tables enabled
  std::vector<int> mEnabledMultiplicityTables; // Vector of enabled multiplicity tables
  bool isMultiplicityTableEnabled(int table) const
  {
    for (int i : mEnabledMultiplicityTables) {
      if (i == table) {
        return true;
      }
    }
    return false;
  }

  // Debug output
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TList> listCalib{"calib-list", OutputObjHandlingPolicy::QAObject};

  unsigned int randomSeed = 0;
  void init(InitContext& context)
  {
    ccdb->setURL(ccdbConfig.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false); // don't fatal, please - exception is caught explicitly (as it should)

    // If both Run 2 and Run 3 data process flags are enabled then we check the metadata
    if (autoSetupFromMetadata && metadataInfo.isFullyDefined()) {
      LOG(info) << "Autosetting the processing from the metadata";
      if (doprocessRun2 == true && doprocessRun3 == true) {
        if (metadataInfo.isRun3()) {
          doprocessRun2.value = false;
        } else {
          doprocessRun3.value = false;
        }
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
    for (int i = 0; i < multiplicity::kNTables; i++) {
      int f = enabledMultiplicityTables->get(multiplicity::tableNames[i].c_str(), "Enable");
      enableFlagIfTableRequired(context, multiplicity::tableNames[i], f);
      if (f == 1) {
        mEnabledMultiplicityTables.push_back(i);
        if (multiplicityConfig.fractionOfEvents <= 1.f && (multiplicity::tableNames[i] != "MultsExtra")) {
          LOG(fatal) << "Cannot have a fraction of events <= 1 and multiplicity table consumed.";
        }
      }
    }
    // Handle the custom cases.
    if (isMultiplicityTableEnabled(multiplicity::kMultMCExtras)) {
      if (enabledMultiplicityTables->get(multiplicity::tableNames[multiplicity::kMultMCExtras].c_str(), "Enable") == -1) {
        doprocessMC.value = true;
        LOG(info) << "Enabling MC processing due to " << multiplicity::tableNames[multiplicity::kMultMCExtras] << " table being enabled.";
      }
    }

    // Check that the tables are enabled consistenly
    if (isMultiplicityTableEnabled(multiplicity::kFV0MultZeqs) && !isMultiplicityTableEnabled(multiplicity::kFV0Mults)) { // FV0
      mEnabledMultiplicityTables.push_back(multiplicity::kFV0Mults);
      LOG(info) << "Cannot have the " << multiplicity::tableNames[multiplicity::kFV0MultZeqs] << " table enabled and not the one on " << multiplicity::tableNames[multiplicity::kFV0Mults] << ". Enabling it.";
    }
    if (isMultiplicityTableEnabled(multiplicity::kFT0MultZeqs) && !isMultiplicityTableEnabled(multiplicity::kFT0Mults)) { // FT0
      mEnabledMultiplicityTables.push_back(multiplicity::kFT0Mults);
      LOG(info) << "Cannot have the " << multiplicity::tableNames[multiplicity::kFT0MultZeqs] << " table enabled and not the one on " << multiplicity::tableNames[multiplicity::kFT0Mults] << ". Enabling it.";
    }
    if (isMultiplicityTableEnabled(multiplicity::kFDDMultZeqs) && !isMultiplicityTableEnabled(multiplicity::kFDDMults)) { // FDD
      mEnabledMultiplicityTables.push_back(multiplicity::kFDDMults);
      LOG(info) << "Cannot have the " << multiplicity::tableNames[multiplicity::kFDDMultZeqs] << " table enabled and not the one on " << multiplicity::tableNames[multiplicity::kFDDMults] << ". Enabling it.";
    }
    if (isMultiplicityTableEnabled(multiplicity::kPVMultZeqs) && !isMultiplicityTableEnabled(multiplicity::kPVMults)) { // PV
      mEnabledMultiplicityTables.push_back(multiplicity::kPVMults);
      LOG(info) << "Cannot have the " << multiplicity::tableNames[multiplicity::kPVMultZeqs] << " table enabled and not the one on " << multiplicity::tableNames[multiplicity::kPVMults] << ". Enabling it.";
    }
    std::sort(mEnabledMultiplicityTables.begin(), mEnabledMultiplicityTables.end());

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

  using BCsWithTimestampsAndRun2Infos = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;
  using Run2Tracks = soa::Join<aod::Tracks, aod::TracksExtra>;
  Partition<Run2Tracks> run2tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));
  Partition<Run2Tracks> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  Partition<Run2Tracks> pvContribTracks = (nabs(aod::track::eta) < 0.8f) && ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Partition<Run2Tracks> pvContribTracksEta1 = (nabs(aod::track::eta) < 1.0f) && ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::TracksIU> perColIU = aod::track::collisionId;
  Preslice<aod::MFTTracks> perCollisionMFT = o2::aod::fwdtrack::collisionId;

  void processRun2(aod::CollisionMatchedRun2Sparse const& collision,
                   Run2Tracks const&,
                   aod::Zdcs const&,
                   aod::FV0As const&,
                   aod::FV0Cs const&,
                   aod::FT0s const&,
                   BCsWithTimestampsAndRun2Infos const& bcs)
  {
    float multFV0A = 0.f;
    float multFV0C = 0.f;
    float multFT0A = 0.f;
    float multFT0C = 0.f;
    float multFDDA = 0.f;
    float multFDDC = 0.f;
    float multZNA = 0.f;
    float multZNC = 0.f;

    auto trackletsGrouped = run2tracklets->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksGrouped = tracksWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    int multTracklets = trackletsGrouped.size();
    int multTPC = tracksGrouped.size();
    int multNContribs = 0;
    int multNContribsEta1 = 0;
    int multNContribsEtaHalf = 0;

    if (collision.has_fv0a()) {
      for (const auto& amplitude : collision.fv0a().amplitude()) {
        multFV0A += amplitude;
      }
    }
    if (collision.has_fv0c()) {
      for (const auto& amplitude : collision.fv0c().amplitude()) {
        multFV0C += amplitude;
      }
    }
    if (collision.has_ft0()) {
      auto ft0 = collision.ft0();
      for (const auto& amplitude : ft0.amplitudeA()) {
        multFT0A += amplitude;
      }
      for (const auto& amplitude : ft0.amplitudeC()) {
        multFT0C += amplitude;
      }
    }
    if (collision.has_zdc()) {
      auto zdc = collision.zdc();
      multZNA = zdc.energyCommonZNA();
      multZNC = zdc.energyCommonZNC();
    }

    // Try to do something Similar to https://github.com/alisw/AliPhysics/blob/22862a945004f719f8e9664c0264db46e7186a48/OADB/AliPPVsMultUtils.cxx#L541C26-L541C37
    for (const auto& tracklet : trackletsGrouped) {
      if (std::abs(tracklet.eta()) < 1.0) {
        multNContribsEta1++;
      }
      if (std::abs(tracklet.eta()) < 0.8) {
        multNContribs++;
      }
      if (std::abs(tracklet.eta()) < 0.5) {
        multNContribsEtaHalf++;
      }
    }

    LOGF(debug, "multFV0A=%5.0f multFV0C=%5.0f multFT0A=%5.0f multFT0C=%5.0f multFDDA=%5.0f multFDDC=%5.0f multZNA=%6.0f multZNC=%6.0f multTracklets=%i multTPC=%i multNContribsEta1=%i multNContribs=%i multNContribsEtaHalf=%i", multFV0A, multFV0C, multFT0A, multFT0C, multFDDA, multFDDC, multZNA, multZNC, multTracklets, multTPC, multNContribs, multNContribsEta1, multNContribsEtaHalf);
    tableFV0(multFV0A, multFV0C);
    tableFT0(multFT0A, multFT0C);
    tableFDD(multFDDA, multFDDC);
    tableZDC(multZNA, multZNC, 0.0f, 0.0f, 0.0f, 0.0f);
    tableTracklet(multTracklets);
    tableTpc(multTPC);
    tablePv(multNContribs, multNContribsEta1, multNContribsEtaHalf);
  }
  PROCESS_SWITCH(MultiplicityTable, processRun2, "Produce Run 2 multiplicity tables. Autoset if both processRun2 and processRun3 are enabled", true);

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
        case multiplicity::kFV0Mults: // FV0
          tableFV0.reserve(collisions.size());
          tableFV0AOuter.reserve(collisions.size());
          break;
        case multiplicity::kFT0Mults: // FT0
          tableFT0.reserve(collisions.size());
          break;
        case multiplicity::kFDDMults: // FDD
          tableFDD.reserve(collisions.size());
          break;
        case multiplicity::kZDCMults: // ZDC
          tableZDC.reserve(collisions.size());
          break;
        case multiplicity::kTrackletMults: // Tracklets (Run 2 only, nothing to do) (to be removed!)
          tableTracklet.reserve(collisions.size());
          break;
        case multiplicity::kTPCMults: // TPC
          tableTpc.reserve(collisions.size());
          break;
        case multiplicity::kPVMults: // PV multiplicity
          tablePv.reserve(collisions.size());
          break;
        case multiplicity::kMultsExtra: // Extra information
          tableExtra.reserve(collisions.size());
          break;
        case multiplicity::kMultSelections: // Extra information
          multSelections.reserve(collisions.size());
          break;
        case multiplicity::kFV0MultZeqs: // Equalized multiplicity for FV0
          tableFV0Zeqs.reserve(collisions.size());
          break;
        case multiplicity::kFT0MultZeqs: // Equalized multiplicity for FT0
          tableFT0Zeqs.reserve(collisions.size());
          break;
        case multiplicity::kFDDMultZeqs: // Equalized multiplicity for FDD
          tableFDDZeqs.reserve(collisions.size());
          break;
        case multiplicity::kPVMultZeqs: // Equalized multiplicity for PV
          tablePVZeqs.reserve(collisions.size());
          break;
        case multiplicity::kMultMCExtras: // MC extra information (nothing to do in the data)
          break;
        default:
          LOG(fatal) << "Unknown table requested: " << tableId;
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
      if ((multiplicityConfig.fractionOfEvents < 1.f) && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > multiplicityConfig.fractionOfEvents) { // Skip events that are not sampled (only for the QA)
        return;
      }
      int multNContribs = 0;
      int multNContribsEta1 = 0;
      int multNContribsEtaHalf = 0;

      /* check the previous run number */
      const auto& bc = collision.bc_as<BCsWithRun3Matchings>();
      initMultCalib(bc);

      // First we compute the multiplicity
      for (const auto tableId : mEnabledMultiplicityTables) {
        switch (tableId) {
          case multiplicity::kFV0Mults: // FV0
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
          case multiplicity::kFT0Mults: // FT0
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
          case multiplicity::kFDDMults: // FDD
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
          case multiplicity::kZDCMults: // ZDC
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
          case multiplicity::kTrackletMults: // Tracklets (only Run2) nothing to do (to be removed!)
          {
            tableTracklet(0);
          } break;
          case multiplicity::kTPCMults: // TPC
          {
            const auto& tracksGrouped = tracksIUWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
            const int multTPC = tracksGrouped.size();
            tableTpc(multTPC);
            LOGF(debug, "multTPC=%i", multTPC);
          } break;
          case multiplicity::kPVMults: // PV multiplicity
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
          case multiplicity::kMultsExtra: // Extra
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
                       bc.runNumber(), collision.posZ(), collision.sel8(),
                       nHasITS, nHasTPC, nHasTOF, nHasTRD, nITSonly, nTPConly, nITSTPC,
                       nAllTracksTPCOnly, nAllTracksITSTPC,
                       collision.trackOccupancyInTimeRange(),
                       collision.ft0cOccupancyInTimeRange(),
                       collision.flags());
          } break;
          case multiplicity::kMultSelections: // Multiplicity selections
          {
            multSelections(collision.selection_raw());
          } break;
          case multiplicity::kFV0MultZeqs: // Z equalized FV0
          {
            if (std::fabs(collision.posZ()) < 15.0f && multCalib.isCalibLoaded) {
              multZeqFV0A = multCalib.hVtxZFV0A->Interpolate(0.0) * multFV0A / multCalib.hVtxZFV0A->Interpolate(collision.posZ());
            }
            tableFV0Zeqs(multZeqFV0A);
          } break;
          case multiplicity::kFT0MultZeqs: // Z equalized FT0
          {
            if (std::fabs(collision.posZ()) < 15.0f && multCalib.isCalibLoaded) {
              multZeqFT0A = multCalib.hVtxZFT0A->Interpolate(0.0) * multFT0A / multCalib.hVtxZFT0A->Interpolate(collision.posZ());
              multZeqFT0C = multCalib.hVtxZFT0C->Interpolate(0.0) * multFT0C / multCalib.hVtxZFT0C->Interpolate(collision.posZ());
            }
            if (produceHistograms.value) {
              histos.fill(HIST("FT0A"), multFT0A, multZeqFT0A);
              histos.fill(HIST("FT0C"), multFT0C, multZeqFT0C);
              histos.fill(HIST("FT0AMultvsPV"), multZeqFT0A, multNContribs);
              histos.fill(HIST("FT0CMultvsPV"), multZeqFT0C, multNContribs);
            }
            tableFT0Zeqs(multZeqFT0A, multZeqFT0C);
          } break;
          case multiplicity::kFDDMultZeqs: // Z equalized FDD
          {
            if (std::fabs(collision.posZ()) < 15.0f && multCalib.isCalibLoaded) {
              multZeqFDDA = multCalib.hVtxZFDDA->Interpolate(0.0) * multFDDA / multCalib.hVtxZFDDA->Interpolate(collision.posZ());
              multZeqFDDC = multCalib.hVtxZFDDC->Interpolate(0.0) * multFDDC / multCalib.hVtxZFDDC->Interpolate(collision.posZ());
            }
            tableFDDZeqs(multZeqFDDA, multZeqFDDC);
          } break;
          case multiplicity::kPVMultZeqs: // Z equalized PV
          {
            if (std::fabs(collision.posZ()) < 15.0f && multCalib.isCalibLoaded) {
              multZeqNContribs = multCalib.hVtxZNTracks->Interpolate(0.0) * multNContribs / multCalib.hVtxZNTracks->Interpolate(collision.posZ());
            }
            tablePVZeqs(multZeqNContribs);
          } break;
          case multiplicity::kMultMCExtras: // MC only (nothing to do)
          {
          } break;
          default: // Default
          {
            LOG(fatal) << "Unknown table requested: " << tableId;
          } break;
        }
      }
    }
  }
  PROCESS_SWITCH(MultiplicityTable, processRun3, "Produce Run 3 multiplicity tables. Autoset if both processRun2 and processRun3 are enabled", true);

  // one loop better than multiple sliceby calls
  // FIT FT0C: -3.3 < η < -2.1
  // FOT FT0A:  3.5 < η <  4.9
  Filter mcParticleFilter = (aod::mcparticle::eta < 7.0f) && (aod::mcparticle::eta > -7.0f);
  using McParticlesFiltered = soa::Filtered<aod::McParticles>;
  void processMC(aod::McCollision const& mcCollision, McParticlesFiltered const& mcParticles)
  {
    int multFT0A = 0;
    int multFV0A = 0;
    int multFT0C = 0;
    int multFDDA = 0;
    int multFDDC = 0;
    int multBarrelEta05 = 0;
    int multBarrelEta08 = 0;
    int multBarrelEta10 = 0;
    for (auto const& mcPart : mcParticles) {
      if (!mcPart.isPhysicalPrimary()) {
        continue;
      }

      auto charge = 0.;
      auto* p = pdg->GetParticle(mcPart.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 1e-3) {
        continue; // reject neutral particles in counters
      }

      if (std::abs(mcPart.eta()) < 1.0) {
        multBarrelEta10++;
        if (std::abs(mcPart.eta()) < 0.8) {
          multBarrelEta08++;
          if (std::abs(mcPart.eta()) < 0.5) {
            multBarrelEta05++;
          }
        }
      }
      if (-3.3 < mcPart.eta() && mcPart.eta() < -2.1)
        multFT0C++;
      if (3.5 < mcPart.eta() && mcPart.eta() < 4.9)
        multFT0A++;
      if (2.2 < mcPart.eta() && mcPart.eta() < 5.0)
        multFV0A++;
      if (-6.9 < mcPart.eta() && mcPart.eta() < -4.9)
        multFDDC++;
      if (4.7 < mcPart.eta() && mcPart.eta() < 6.3)
        multFDDA++;
    }
    tableExtraMc(multFT0A, multFT0C, multFV0A, multFDDA, multFDDC, multBarrelEta05, multBarrelEta08, multBarrelEta10, mcCollision.posZ());
  }
  PROCESS_SWITCH(MultiplicityTable, processMC, "Produce MC multiplicity tables", false);

  void processMC2Mults(soa::Join<aod::McCollisionLabels, aod::Collisions>::iterator const& collision)
  {
    tableExtraMult2MCExtras(collision.mcCollisionId()); // interlink
  }
  PROCESS_SWITCH(MultiplicityTable, processMC2Mults, "Produce MC -> Mult map", false);

  using Run3Tracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;
  Partition<Run3Tracks> pvContribGlobalTracksEta1 = (minPtGlobalTrack < aod::track::pt && aod::track::pt < maxPtGlobalTrack) && (nabs(aod::track::eta) < 1.0f) && ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor)) && requireQualityTracksInFilter();
  void processGlobalTrackingCounters(aod::Collision const& collision, soa::Join<Run3TracksIU, aod::TrackSelection, aod::TrackSelectionExtension> const& tracksIU, Run3Tracks const&)
  {
    // counter from Igor
    int nGlobalTracks = 0;
    int multNbrContribsEta05GlobalTrackWoDCA = 0;
    int multNbrContribsEta08GlobalTrackWoDCA = 0;
    int multNbrContribsEta10GlobalTrackWoDCA = 0;

    auto pvContribGlobalTracksEta1PerCollision = pvContribGlobalTracksEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (const auto& track : pvContribGlobalTracksEta1PerCollision) {
      if (track.itsNCls() < multiplicityConfig.minNclsITSGlobalTrack || track.itsNClsInnerBarrel() < multiplicityConfig.minNclsITSibGlobalTrack) {
        continue;
      }
      multNbrContribsEta10GlobalTrackWoDCA++;

      if (std::abs(track.eta()) < 0.8) {
        multNbrContribsEta08GlobalTrackWoDCA++;
      }
      if (std::abs(track.eta()) < 0.5) {
        multNbrContribsEta05GlobalTrackWoDCA++;
      }
    }

    for (const auto& track : tracksIU) {
      if (std::fabs(track.eta()) < 0.8 && track.tpcNClsFound() >= 80 && track.tpcNClsCrossedRows() >= 100) {
        if (track.isGlobalTrack()) {
          nGlobalTracks++;
        }
      }
    }

    LOGF(debug, "nGlobalTracks = %d, multNbrContribsEta08GlobalTrackWoDCA = %d, multNbrContribsEta10GlobalTrackWoDCA = %d, multNbrContribsEta05GlobalTrackWoDCA = %d", nGlobalTracks, multNbrContribsEta08GlobalTrackWoDCA, multNbrContribsEta10GlobalTrackWoDCA, multNbrContribsEta05GlobalTrackWoDCA);

    multsGlobal(nGlobalTracks, multNbrContribsEta08GlobalTrackWoDCA, multNbrContribsEta10GlobalTrackWoDCA, multNbrContribsEta05GlobalTrackWoDCA);
  }
  PROCESS_SWITCH(MultiplicityTable, processGlobalTrackingCounters, "Produce Run 3 global counters", false);

  void processRun3MFT(soa::Join<aod::Collisions, aod::EvSels>::iterator const&,
                      o2::aod::MFTTracks const& mftTracks,
                      soa::SmallGroups<aod::BestCollisionsFwd> const& retracks)
  {
    int nAllTracks = 0;
    int nTracks = 0;

    for (const auto& track : mftTracks) {
      if (track.nClusters() >= 5) { // hardcoded for now
        nAllTracks++;
      }
    }

    if (retracks.size() > 0) {
      for (const auto& retrack : retracks) {
        auto track = retrack.mfttrack();
        if (track.nClusters() < 5) {
          continue; // min cluster requirement
        }
        if ((track.eta() > -2.0f) && (track.eta() < -3.9f)) {
          continue; // too far to be of true interest
        }
        if (std::abs(retrack.bestDCAXY()) > 2.0f) {
          continue; // does not point to PV properly
        }
        nTracks++;
      }
    }
    mftMults(nAllTracks, nTracks);
  }
  PROCESS_SWITCH(MultiplicityTable, processRun3MFT, "Produce MFT mult tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<MultiplicityTable>(cfgc)};
}
