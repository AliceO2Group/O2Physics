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
static constexpr int nTables = 14;

// Checking that the Zeq tables are after the normal ones
static_assert(kFV0Mults < kFV0MultZeqs);
static_assert(kFT0Mults < kFT0MultZeqs);
static_assert(kFDDMults < kFDDMultZeqs);
static_assert(kPVMults < kPVMultZeqs);

static constexpr int nParameters = 1;
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
static const std::vector<std::string> parameterNames{"Enable"};
static const int defaultParameters[nTables][nParameters]{{-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}};

struct MultiplicityTable {
  SliceCache cache;
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

  // For vertex-Z corrections in calibration
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  using Run2Tracks = soa::Join<aod::Tracks, aod::TracksExtra>;
  Partition<Run2Tracks> run2tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));
  Partition<Run2Tracks> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  Partition<Run2Tracks> pvContribTracks = (nabs(aod::track::eta) < 0.8f) && ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);
  Partition<Run2Tracks> pvContribTracksEta1 = (nabs(aod::track::eta) < 1.0f) && ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);
  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::TracksIU> perColIU = aod::track::collisionId;
  Preslice<aod::MFTTracks> perCollisionMFT = o2::aod::fwdtrack::collisionId;

  using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;

  // Configurable
  Configurable<int> doVertexZeq{"doVertexZeq", 1, "if 1: do vertex Z eq mult table"};
  Configurable<float> fractionOfEvents{"fractionOfEvents", 2.0, "Fractions of events to keep in case the QA is used"};
  Configurable<LabeledArray<int>> enabledTables{"enabledTables",
                                                {defaultParameters[0], nTables, nParameters, tableNames, parameterNames},
                                                "Produce tables depending on needs. Values different than -1 override the automatic setup: the corresponding table can be set off (0) or on (1)"};

  struct : ConfigurableGroup {
    Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "The CCDB endpoint url address"};
    Configurable<std::string> ccdbPath{"ccdbpath", "Centrality/Calibration", "The CCDB path for centrality/multiplicity information"};
    Configurable<std::string> reconstructionPass{"reconstructionPass", "", {"Apass to use when fetching the calibration tables. Empty (default) does not check for any pass. Use `metadata` to fetch it from the AO2D metadata. Otherwise it will override the metadata."}};
  } ccdbConfig;

  Configurable<bool> produceHistograms{"produceHistograms", false, {"Option to produce debug histograms"}};

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
    if (metadataInfo.isFullyDefined() && !doprocessRun2 && !doprocessRun3) { // Check if the metadata is initialized (only if not forced from the workflow configuration)
      if (metadataInfo.isRun3()) {
        doprocessRun3.value = true;
      } else {
        doprocessRun2.value = false;
      }
    }

    randomSeed = static_cast<unsigned int>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    if (doprocessRun2 == false && doprocessRun3 == false) {
      LOGF(fatal, "Neither processRun2 nor processRun3 enabled. Please choose one.");
    }
    if (doprocessRun2 == true && doprocessRun3 == true) {
      LOGF(fatal, "Cannot enable processRun2 and processRun3 at the same time. Please choose one.");
    }

    bool tEnabled[nTables] = {false};
    for (int i = 0; i < nTables; i++) {
      int f = enabledTables->get(tableNames[i].c_str(), "Enable");
      enableFlagIfTableRequired(context, tableNames[i], f);
      if (f == 1) {
        tEnabled[i] = true;
        mEnabledTables.push_back(i);
        if (fractionOfEvents <= 1.f && (tableNames[i] != "MultsExtra")) {
          LOG(fatal) << "Cannot have a fraction of events <= 1 and multiplicity table consumed.";
        }
      }
    }
    // Handle the custom cases.
    if (tEnabled[kMultMCExtras]) {
      if (enabledTables->get(tableNames[kMultMCExtras].c_str(), "Enable") == -1) {
        doprocessMC.value = true;
        LOG(info) << "Enabling MC processing due to " << tableNames[kMultMCExtras] << " table being enabled.";
      }
    }

    // Check that the tables are enabled consistenly
    if (tEnabled[kFV0MultZeqs] && !tEnabled[kFV0Mults]) { // FV0
      mEnabledTables.push_back(kFV0Mults);
      LOG(info) << "Cannot have the " << tableNames[kFV0MultZeqs] << " table enabled and not the one on " << tableNames[kFV0Mults] << ". Enabling it.";
    }
    if (tEnabled[kFT0MultZeqs] && !tEnabled[kFT0Mults]) { // FT0
      mEnabledTables.push_back(kFT0Mults);
      LOG(info) << "Cannot have the " << tableNames[kFT0MultZeqs] << " table enabled and not the one on " << tableNames[kFT0Mults] << ". Enabling it.";
    }
    if (tEnabled[kFDDMultZeqs] && !tEnabled[kFDDMults]) { // FDD
      mEnabledTables.push_back(kFDDMults);
      LOG(info) << "Cannot have the " << tableNames[kFDDMultZeqs] << " table enabled and not the one on " << tableNames[kFDDMults] << ". Enabling it.";
    }
    if (tEnabled[kPVMultZeqs] && !tEnabled[kPVMults]) { // PV
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

    ccdb->setURL(ccdbConfig.ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false); // don't fatal, please - exception is caught explicitly (as it should)

    if (!produceHistograms.value) {
      return;
    }
    histos.add("FT0A", "FT0A vs FT0A eq.", HistType::kTH2D, {{1000, 0, 1000, "FT0A multiplicity"}, {1000, 0, 1000, "FT0A multiplicity eq."}});
    histos.add("FT0C", "FT0C vs FT0C eq.", HistType::kTH2D, {{1000, 0, 1000, "FT0C multiplicity"}, {1000, 0, 1000, "FT0C multiplicity eq."}});
    histos.add("FT0CMultvsPV", "FT0C vs mult.", HistType::kTH2D, {{1000, 0, 1000, "FT0C mult."}, {100, 0, 100, "PV mult."}});
    histos.add("FT0AMultvsPV", "FT0A vs mult.", HistType::kTH2D, {{1000, 0, 1000, "FT0A mult."}, {100, 0, 100, "PV mult."}});

    listCalib.setObject(new TList);
  }

  /// Dummy process function for BCs, needed in case both Run2 and Run3 process functions are disabled
  void process(aod::BCs const&) {}

  void processRun2(aod::Run2MatchedSparse::iterator const& collision,
                   Run2Tracks const&,
                   aod::BCs const&,
                   aod::Zdcs const&,
                   aod::FV0As const&,
                   aod::FV0Cs const&,
                   aod::FT0s const&)
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
      for (auto amplitude : collision.fv0a().amplitude()) {
        multFV0A += amplitude;
      }
    }
    if (collision.has_fv0c()) {
      for (auto amplitude : collision.fv0c().amplitude()) {
        multFV0C += amplitude;
      }
    }
    if (collision.has_ft0()) {
      auto ft0 = collision.ft0();
      for (auto amplitude : ft0.amplitudeA()) {
        multFT0A += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        multFT0C += amplitude;
      }
    }
    if (collision.has_zdc()) {
      auto zdc = collision.zdc();
      multZNA = zdc.energyCommonZNA();
      multZNC = zdc.energyCommonZNC();
    }

    LOGF(debug, "multFV0A=%5.0f multFV0C=%5.0f multFT0A=%5.0f multFT0C=%5.0f multFDDA=%5.0f multFDDC=%5.0f multZNA=%6.0f multZNC=%6.0f multTracklets=%i multTPC=%i", multFV0A, multFV0C, multFT0A, multFT0C, multFDDA, multFDDC, multZNA, multZNC, multTracklets, multTPC);
    tableFV0(multFV0A, multFV0C);
    tableFT0(multFT0A, multFT0C);
    tableFDD(multFDDA, multFDDC);
    tableZDC(multZNA, multZNC, 0.0f, 0.0f, 0.0f, 0.0f);
    tableTracklet(multTracklets);
    tableTpc(multTPC);
    tablePv(multNContribs, multNContribsEta1, multNContribsEtaHalf);
  }

  using Run3TracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;
  Partition<Run3TracksIU> tracksIUWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  Partition<Run3TracksIU> pvAllContribTracksIU = ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);
  Partition<Run3TracksIU> pvContribTracksIU = (nabs(aod::track::eta) < 0.8f) && ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);
  Partition<Run3TracksIU> pvContribTracksIUEta1 = (nabs(aod::track::eta) < 1.0f) && ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);
  Partition<Run3TracksIU> pvContribTracksIUEtaHalf = (nabs(aod::track::eta) < 0.5f) && ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);

  void processRun3(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                   Run3TracksIU const&,
                   BCsWithRun3Matchings const&,
                   aod::Zdcs const&,
                   aod::FV0As const&,
                   aod::FT0s const&,
                   aod::FDDs const&)
  {
    // reserve memory
    for (auto i : mEnabledTables) {
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

      for (auto i : mEnabledTables) {
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
              for (auto amplitude : ft0.amplitudeA()) {
                multFT0A += amplitude;
              }
              for (auto amplitude : ft0.amplitudeC()) {
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
              for (auto amplitude : fdd.chargeA()) {
                multFDDA += amplitude;
              }
              for (auto amplitude : fdd.chargeC()) {
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
            for (auto track : tracksThisCollision) {
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

            for (auto track : pvAllContribsGrouped) {
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
            for (auto track : tpcTracksGrouped) {
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
            if (fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
              multZeqFV0A = hVtxZFV0A->Interpolate(0.0) * multFV0A / hVtxZFV0A->Interpolate(collision.posZ());
            }
            tableFV0Zeqs(multZeqFV0A);
          } break;
          case kFT0MultZeqs: // Z equalized FT0
          {
            if (fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
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
            if (fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
              multZeqFDDA = hVtxZFDDA->Interpolate(0.0) * multFDDA / hVtxZFDDA->Interpolate(collision.posZ());
              multZeqFDDC = hVtxZFDDC->Interpolate(0.0) * multFDDC / hVtxZFDDC->Interpolate(collision.posZ());
            }
            tableFDDZeqs(multZeqFDDA, multZeqFDDC);
          } break;
          case kPVMultZeqs: // Z equalized PV
          {
            if (fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
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

  // one loop better than multiple sliceby calls
  // FIT FT0C: -3.3 < η < -2.1
  // FOT FT0A:  3.5 < η <  4.9
  Filter mcParticleFilter = (aod::mcparticle::eta < 7.0f) && (aod::mcparticle::eta > -7.0f);
  using mcParticlesFiltered = soa::Filtered<aod::McParticles>;

  void processMC(aod::McCollision const& mcCollision, mcParticlesFiltered const& mcParticles)
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

  void processMC2Mults(soa::Join<aod::McCollisionLabels, aod::Collisions>::iterator const& collision)
  {
    tableExtraMult2MCExtras(collision.mcCollisionId()); // interlink
  }

  Configurable<float> min_pt_globaltrack{"min_pt_globaltrack", 0.15, "min. pT for global tracks"};
  Configurable<float> max_pt_globaltrack{"max_pt_globaltrack", 1e+10, "max. pT for global tracks"};
  Configurable<int> min_ncluster_its_globaltrack{"min_ncluster_its_globaltrack", 5, "min. number of ITS clusters for global tracks"};
  Configurable<int> min_ncluster_itsib_globaltrack{"min_ncluster_itsib_globaltrack", 1, "min. number of ITSib clusters for global tracks"};

  using Run3Tracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;
  Partition<Run3Tracks> pvContribGlobalTracksEta1 = (min_pt_globaltrack < aod::track::pt && aod::track::pt < max_pt_globaltrack) && (nabs(aod::track::eta) < 1.0f) && ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor) && requireQualityTracksInFilter();

  void processGlobalTrackingCounters(aod::Collision const& collision, soa::Join<Run3TracksIU, aod::TrackSelection, aod::TrackSelectionExtension> const& tracksIU, Run3Tracks const&)
  {
    // counter from Igor
    int nGlobalTracks = 0;
    int multNContribsEta05_kGlobalTrackWoDCA = 0;
    int multNContribsEta08_kGlobalTrackWoDCA = 0;
    int multNContribsEta10_kGlobalTrackWoDCA = 0;

    auto pvContribGlobalTracksEta1_per_collision = pvContribGlobalTracksEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (auto& track : pvContribGlobalTracksEta1_per_collision) {
      if (track.itsNCls() < min_ncluster_its_globaltrack || track.itsNClsInnerBarrel() < min_ncluster_itsib_globaltrack) {
        continue;
      }
      multNContribsEta10_kGlobalTrackWoDCA++;

      if (std::abs(track.eta()) < 0.8) {
        multNContribsEta08_kGlobalTrackWoDCA++;
      }
      if (std::abs(track.eta()) < 0.5) {
        multNContribsEta05_kGlobalTrackWoDCA++;
      }
    }

    for (auto& track : tracksIU) {
      if (fabs(track.eta()) < 0.8 && track.tpcNClsFound() >= 80 && track.tpcNClsCrossedRows() >= 100) {
        if (track.isGlobalTrack()) {
          nGlobalTracks++;
        }
      }
    }

    LOGF(debug, "nGlobalTracks = %d, multNContribsEta08_kGlobalTrackWoDCA = %d, multNContribsEta10_kGlobalTrackWoDCA = %d, multNContribsEta05_kGlobalTrackWoDCA = %d", nGlobalTracks, multNContribsEta08_kGlobalTrackWoDCA, multNContribsEta10_kGlobalTrackWoDCA, multNContribsEta05_kGlobalTrackWoDCA);

    multsGlobal(nGlobalTracks, multNContribsEta08_kGlobalTrackWoDCA, multNContribsEta10_kGlobalTrackWoDCA, multNContribsEta05_kGlobalTrackWoDCA);
  }

  void processRun3MFT(soa::Join<aod::Collisions, aod::EvSels>::iterator const&,
                      o2::aod::MFTTracks const& mftTracks,
                      soa::SmallGroups<aod::BestCollisionsFwd> const& retracks)
  {
    int nAllTracks = 0;
    int nTracks = 0;

    for (auto& track : mftTracks) {
      if (track.nClusters() >= 5) { // hardcoded for now
        nAllTracks++;
      }
    }

    if (retracks.size() > 0) {
      for (auto& retrack : retracks) {
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

  // Process switches
  PROCESS_SWITCH(MultiplicityTable, processRun2, "Produce Run 2 multiplicity tables", false);
  PROCESS_SWITCH(MultiplicityTable, processRun3, "Produce Run 3 multiplicity tables", true);
  PROCESS_SWITCH(MultiplicityTable, processGlobalTrackingCounters, "Produce Run 3 global counters", false);
  PROCESS_SWITCH(MultiplicityTable, processMC, "Produce MC multiplicity tables", false);
  PROCESS_SWITCH(MultiplicityTable, processMC2Mults, "Produce MC -> Mult map", false);
  PROCESS_SWITCH(MultiplicityTable, processRun3MFT, "Produce MFT mult tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<MultiplicityTable>(cfgc)};
}
