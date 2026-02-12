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

/// \file MultModule.h
/// \brief combined multiplicity + centrality module with autodetect features
/// \author ALICE

#ifndef COMMON_TOOLS_MULTIPLICITY_MULTMODULE_H_
#define COMMON_TOOLS_MULTIPLICITY_MULTMODULE_H_

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/DataSpecUtils.h>
#include <Framework/DataTypes.h>
#include <Framework/DeviceSpec.h>
#include <Framework/Logger.h>
#include <Framework/RunningWorkflowInfo.h>

#include <TFile.h>
#include <TFormula.h>
#include <TH1.h>
#include <TList.h>
#include <TProfile.h>
#include <TString.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

//__________________________________________
// MultModule

namespace o2
{
namespace common
{
namespace multiplicity
{

// statics necessary for the configurables in this namespace
static constexpr int nParameters = 1;
static const std::vector<std::string> tableNames{
  // multiplicity subcomponent
  "FV0Mults",
  "FV0AOuterMults",
  "FT0Mults",
  "FDDMults",
  "ZDCMults",
  "TrackletMults",
  "TPCMults",
  "PVMults",
  "MultsExtra",
  "MultSelections",
  "FV0MultZeqs",
  "FT0MultZeqs",
  "FDDMultZeqs",
  "PVMultZeqs",
  "GlobalMultZeqs",
  "MFTMultZeqs",
  "MultMCExtras",
  "Mult2MCExtras",
  "MFTMults",
  "MultsGlobal",

  // centrality subcomponent
  "CentRun2V0Ms",
  "CentRun2V0As",
  "CentRun2SPDTrks",
  "CentRun2SPDClss",
  "CentRun2CL0s",
  "CentRun2CL1s",
  "CentFV0As",
  "CentFT0Ms",
  "CentFT0As",
  "CentFT0Cs",
  "CentFT0CVariant1s",
  "CentFT0CVariant2s",
  "CentFDDMs",
  "CentNTPVs",
  "CentNGlobals",
  "CentMFTs",
  "BCCentFT0Ms",
  "BCCentFT0As",
  "BCCentFT0Cs"};

static constexpr int nTablesConst = 39;

static const std::vector<std::string> parameterNames{"enable"};
static const int defaultParameters[nTablesConst][nParameters]{
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1}};

// table index : match order above
enum tableIndex { kFV0Mults,       // standard
                  kFV0AOuterMults, // standard
                  kFT0Mults,       // standard
                  kFDDMults,       // standard
                  kZDCMults,       // standard
                  kTrackletMults,  // Run 2
                  kTPCMults,       // standard
                  kPVMults,        // standard
                  kMultsExtra,     // standard
                  kMultSelections, // event selection
                  kFV0MultZeqs,    // zeq calib, standard
                  kFT0MultZeqs,    // zeq calib, standard
                  kFDDMultZeqs,    // zeq calib, standard
                  kPVMultZeqs,     // zeq calib, standard
                  kGlobalMultZeqs, // zeq calib, extra
                  kMFTMultZeqs,    // zeq calib, extra
                  kMultMCExtras,   // MC exclusive
                  kMult2MCExtras,  // MC exclusive
                  kMFTMults,       // requires MFT task
                  kMultsGlobal,    // requires track selection task

                  // centrality subcomponent
                  kCentRun2V0Ms,      // Run 2
                  kCentRun2V0As,      // Run 2
                  kCentRun2SPDTrks,   // Run 2
                  kCentRun2SPDClss,   // Run 2
                  kCentRun2CL0s,      // Run 2
                  kCentRun2CL1s,      // Run 2
                  kCentFV0As,         // standard Run 3
                  kCentFT0Ms,         // standard Run 3
                  kCentFT0As,         // standard Run 3
                  kCentFT0Cs,         // standard Run 3
                  kCentFT0CVariant1s, // standard Run 3
                  kCentFT0CVariant2s, // standard Run 3
                  kCentFDDMs,         // standard Run 3
                  kCentNTPVs,         // standard Run 3
                  kCentNGlobals,      // requires track selection task
                  kCentMFTs,          // requires MFT task
                  kBCCentFT0Ms,       // bc centrality
                  kBCCentFT0As,       // bc centrality
                  kBCCentFT0Cs,       // bc centrality
                  kNTables };

struct products : o2::framework::ProducesGroup {
  //__________________________________________________
  // multiplicity tables
  o2::framework::Produces<aod::FV0Mults> tableFV0;
  o2::framework::Produces<aod::FV0AOuterMults> tableFV0AOuter;
  o2::framework::Produces<aod::FT0Mults> tableFT0;
  o2::framework::Produces<aod::FDDMults> tableFDD;
  o2::framework::Produces<aod::ZDCMults> tableZDC;
  o2::framework::Produces<aod::TrackletMults> tableTracklet;
  o2::framework::Produces<aod::TPCMults> tableTpc;
  o2::framework::Produces<aod::PVMults> tablePv;
  o2::framework::Produces<aod::MultsExtra> tableExtra;
  o2::framework::Produces<aod::MultSelections> multSelections;
  o2::framework::Produces<aod::FV0MultZeqs> tableFV0Zeqs;
  o2::framework::Produces<aod::FT0MultZeqs> tableFT0Zeqs;
  o2::framework::Produces<aod::FDDMultZeqs> tableFDDZeqs;
  o2::framework::Produces<aod::PVMultZeqs> tablePVZeqs;
  o2::framework::Produces<aod::GlobalMultZeqs> tableNGlobalZeqs;
  o2::framework::Produces<aod::MFTMultZeqs> tableNMFTZeqs;
  o2::framework::Produces<aod::MultMCExtras> tableExtraMc;
  o2::framework::Produces<aod::Mult2MCExtras> tableExtraMult2MCExtras;
  o2::framework::Produces<aod::MFTMults> mftMults;
  o2::framework::Produces<aod::MultsGlobal> multsGlobal;

  //__________________________________________________
  // centrality tables (per collision / default)
  o2::framework::Produces<aod::CentRun2V0Ms> centRun2V0M;
  o2::framework::Produces<aod::CentRun2V0As> centRun2V0A;
  o2::framework::Produces<aod::CentRun2SPDTrks> centRun2SPDTracklets;
  o2::framework::Produces<aod::CentRun2SPDClss> centRun2SPDClusters;
  o2::framework::Produces<aod::CentRun2CL0s> centRun2CL0;
  o2::framework::Produces<aod::CentRun2CL1s> centRun2CL1;
  o2::framework::Produces<aod::CentFV0As> centFV0A;
  o2::framework::Produces<aod::CentFT0Ms> centFT0M;
  o2::framework::Produces<aod::CentFT0As> centFT0A;
  o2::framework::Produces<aod::CentFT0Cs> centFT0C;
  o2::framework::Produces<aod::CentFT0CVariant1s> centFT0CVariant1;
  o2::framework::Produces<aod::CentFT0CVariant2s> centFT0CVariant2;
  o2::framework::Produces<aod::CentFDDMs> centFDDM;
  o2::framework::Produces<aod::CentNTPVs> centNTPV;
  o2::framework::Produces<aod::CentNGlobals> centNGlobals;
  o2::framework::Produces<aod::CentMFTs> centMFTs;
  o2::framework::Produces<aod::BCCentFT0As> bcCentFT0A;
  o2::framework::Produces<aod::BCCentFT0Cs> bcCentFT0C;
  o2::framework::Produces<aod::BCCentFT0Ms> bcCentFT0M;

  //__________________________________________________
  // centrality tables per BC
  // FIXME - future development
};

// for providing temporary buffer
// FIXME ideally cursors could be readable
// to avoid duplicate memory allocation but ok
struct multEntry {
  float multFV0A = 0.0f;
  float multFV0C = 0.0f;
  float multFV0AOuter = 0.0f;
  float multFT0A = 0.0f;
  float multFT0C = 0.0f;
  float multFDDA = 0.0f;
  float multFDDC = 0.0f;
  float multZNA = 0.0f;
  float multZNC = 0.0f;
  float multZEM1 = 0.0f;
  float multZEM2 = 0.0f;
  float multZPA = 0.0f;
  float multZPC = 0.0f;
  int multTracklets = 0;

  int multNContribs = 0;        // PVMult 0.8
  int multNContribsEta1 = 0;    // PVMult 1.0
  int multNContribsEtaHalf = 0; // PVMult 0.5
  int multTPC = 0;              // all TPC (PV contrib unchecked)
  int multHasTPC = 0;           // extras
  int multHasITS = 0;           // extras
  int multHasTOF = 0;           // extras
  int multHasTRD = 0;           // extras
  int multITSOnly = 0;          // extras
  int multTPCOnly = 0;          // extras
  int multITSTPC = 0;           // extras
  int multAllTracksTPCOnly = 0; // extras
  int multAllTracksITSTPC = 0;  // extras

  float multFV0AZeq = -999.0f;
  float multFV0CZeq = -999.0f;
  float multFT0AZeq = -999.0f;
  float multFT0CZeq = -999.0f;
  float multFDDAZeq = -999.0f;
  float multFDDCZeq = -999.0f;
  float multNContribsZeq = 0;
  float multMFTTracksZeq = 0;
  float multGlobalTracksZeq = 0;

  int multGlobalTracks = 0;                     // multsGlobal
  int multNbrContribsEta05GlobalTrackWoDCA = 0; // multsGlobal
  int multNbrContribsEta08GlobalTrackWoDCA = 0; // multsGlobal
  int multNbrContribsEta10GlobalTrackWoDCA = 0; // multsGlobal

  int multMFTAllTracks = 0; // mft
  int multMFTTracks = 0;    // mft

  // For Run2 only
  float posZ = -999.0f;
  uint16_t spdClustersL0 = 0;
  uint16_t spdClustersL1 = 0;
};

// strangenessBuilder: 1st-order configurables
struct standardConfigurables : o2::framework::ConfigurableGroup {
  // self-configuration configurables
  o2::framework::Configurable<o2::framework::LabeledArray<int>> enabledTables{"enabledTables",
                                                                              {defaultParameters[0], nTablesConst, nParameters, tableNames, parameterNames},
                                                                              "Produce this table: -1 for autodetect; otherwise, 0/1 is false/true"};
  std::vector<int> mEnabledTables; // Vector of enabled tables

  // Autoconfigure process functions
  o2::framework::Configurable<bool> autoConfigureProcess{"autoConfigureProcess", false, "if true, will configure process function switches based on metadata"};

  // do vertex-Z equalized or not
  o2::framework::Configurable<int> doVertexZeq{"doVertexZeq", 1, "if 1: do vertex Z eq mult table"};

  // global track counter configurables
  o2::framework::Configurable<float> minPtGlobalTrack{"minPtGlobalTrack", 0.15, "min. pT for global tracks"};
  o2::framework::Configurable<float> maxPtGlobalTrack{"maxPtGlobalTrack", 1e+10, "max. pT for global tracks"};
  o2::framework::Configurable<int> minNclsITSGlobalTrack{"minNclsITSGlobalTrack", 5, "min. number of ITS clusters for global tracks"};
  o2::framework::Configurable<int> minNclsITSibGlobalTrack{"minNclsITSibGlobalTrack", 1, "min. number of ITSib clusters for global tracks"};

  // MFT track counter configurables
  o2::framework::Configurable<int> minNclsMFTTrack{"minNclsMFTTrack", 5, "min. number of MFT clusters for MFT tracks"};
  o2::framework::Configurable<float> maxDCAxyToPVMFTTrack{"maxDCAxyToPVMFTTrack", 2.0f, "max DCAxy to PV for MFT tracks (cm)"};
  o2::framework::Configurable<float> minEtaMFTTrack{"minEtaMFTTrack", -1e+09f, "min. pseudorapidity for MFT tracks (nominal: -3.6)"};
  o2::framework::Configurable<float> maxEtaMFTTrack{"maxEtaMFTTrack", 1e+09f, "max. pseudorapidity for MFT tracks (nominal: -2.45)"};

  // ccdb information
  o2::framework::Configurable<std::string> ccdbPathVtxZ{"ccdbPathVtxZ", "Centrality/Calibration", "The CCDB path for vertex-Z calibration"};
  o2::framework::Configurable<std::string> ccdbPathCentrality{"ccdbPathCentrality", "Centrality/Estimators", "The CCDB path for centrality information"};
  o2::framework::Configurable<std::string> reconstructionPass{"reconstructionPass", "", {"Apass to use when fetching the calibration tables. Empty (default) does not check for any pass. Use `metadata` to fetch it from the AO2D metadata. Otherwise it will override the metadata."}};

  // centrality operation
  o2::framework::Configurable<std::string> generatorName{"generatorName", "", {"Specify if and only if this is MC. Typical: PYTHIA"}};
  o2::framework::Configurable<bool> embedINELgtZEROselection{"embedINELgtZEROselection", false, {"Option to do percentile 100.5 if not INELgtZERO"}};
};

class MultModule
{
 public:
  MultModule()
  {
    // constructor
    mRunNumber = 0;
    mRunNumberCentrality = 0;
    lCalibLoaded = false;
    lCalibObjects = nullptr;
    hVtxZFV0A = nullptr;
    hVtxZFT0A = nullptr;
    hVtxZFT0C = nullptr;
    hVtxZFDDA = nullptr;
    hVtxZFDDC = nullptr;
    hVtxZNTracks = nullptr;
    hVtxZNMFTTracks = nullptr;
    hVtxZNGlobalTracks = nullptr;
  }

  // internal: calib related, vtx-z profiles
  int mRunNumber;
  int mRunNumberCentrality;
  bool lCalibLoaded;
  TList* lCalibObjects;
  TProfile* hVtxZFV0A;
  TProfile* hVtxZFT0A;
  TProfile* hVtxZFT0C;
  TProfile* hVtxZFDDA;
  TProfile* hVtxZFDDC;
  TProfile* hVtxZNTracks;
  TProfile* hVtxZNMFTTracks;    // non-legacy, added August/2025
  TProfile* hVtxZNGlobalTracks; // non-legacy, added August/2025

  // declaration of structs here
  // (N.B.: will be invisible to the outside, create your own copies)
  o2::common::multiplicity::standardConfigurables internalOpts;

  //_________________________________________________
  // centrality-related objects
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
  struct CalibrationInfo {
    std::string name = "";
    bool mCalibrationStored = false;
    TH1* mhMultSelCalib = nullptr;
    float mMCScalePars[6] = {0.0};
    TFormula* mMCScale = nullptr;
    explicit CalibrationInfo(std::string name)
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

  CalibrationInfo fv0aInfo = CalibrationInfo("FV0");
  CalibrationInfo ft0mInfo = CalibrationInfo("FT0");
  CalibrationInfo ft0aInfo = CalibrationInfo("FT0A");
  CalibrationInfo ft0cInfo = CalibrationInfo("FT0C");
  CalibrationInfo ft0cVariant1Info = CalibrationInfo("FT0Cvar1");
  CalibrationInfo ft0cVariant2Info = CalibrationInfo("FT0Cvar2");
  CalibrationInfo fddmInfo = CalibrationInfo("FDD");
  CalibrationInfo ntpvInfo = CalibrationInfo("NTracksPV");
  CalibrationInfo nGlobalInfo = CalibrationInfo("NGlobal");
  CalibrationInfo mftInfo = CalibrationInfo("MFT");

  template <typename TMetadatainfo, typename TConfigurables, typename TInitContext>
  void init(TMetadatainfo const& metadataInfo, TConfigurables& opts, TInitContext& context)
  {
    // read in configurations from the task where it's used
    internalOpts = opts;
    internalOpts.mEnabledTables.resize(nTablesConst, 0);

    LOGF(info, "Configuring tables to generate");
    LOGF(info, "Metadata information: isMC? %i", metadataInfo.isMC());
    const auto& workflows = context.services().template get<o2::framework::RunningWorkflowInfo const>();

    TString listOfRequestors[nTablesConst];
    for (int i = 0; i < nTablesConst; i++) {
      int f = internalOpts.enabledTables->get(tableNames[i].c_str(), "enable");
      if (f == 1) {
        internalOpts.mEnabledTables[i] = 1;
        listOfRequestors[i] = "manual enabling";
      }
      if (f == -1) {
        // autodetect this table in other devices
        for (o2::framework::DeviceSpec const& device : workflows.devices) {
          // Step 1: check if this device subscribed to the V0data table
          for (auto const& input : device.inputs) {
            if (o2::framework::DataSpecUtils::partialMatch(input.matcher, o2::header::DataOrigin("AOD"))) {
              auto&& [origin, description, version] = o2::framework::DataSpecUtils::asConcreteDataMatcher(input.matcher);
              std::string tableNameWithVersion = tableNames[i];
              if (version > 0) {
                tableNameWithVersion += Form("_%03d", version);
              }
              if (input.matcher.binding == tableNameWithVersion) {
                LOGF(info, "Device %s has subscribed to %s (version %i)", device.name, tableNames[i], version);
                listOfRequestors[i].Append(Form("%s ", device.name.c_str()));
                internalOpts.mEnabledTables[i] = 1;
              }
            }
          }
        }
      }
    }

    // dependency checker
    if (internalOpts.mEnabledTables[kCentFV0As] && !internalOpts.mEnabledTables[kFV0MultZeqs]) {
      internalOpts.mEnabledTables[kFV0MultZeqs] = 1;
      listOfRequestors[kFV0MultZeqs].Append(Form("%s ", "dependency check"));
    }
    if ((internalOpts.mEnabledTables[kCentFT0As] || internalOpts.mEnabledTables[kCentFT0Cs] || internalOpts.mEnabledTables[kCentFT0Ms] || internalOpts.mEnabledTables[kCentFT0CVariant1s]) && !internalOpts.mEnabledTables[kFT0MultZeqs]) {
      internalOpts.mEnabledTables[kFT0MultZeqs] = 1;
      listOfRequestors[kFT0MultZeqs].Append(Form("%s ", "dependency check"));
    }
    if (internalOpts.mEnabledTables[kCentFDDMs] && !internalOpts.mEnabledTables[kFDDMultZeqs]) {
      internalOpts.mEnabledTables[kFDDMultZeqs] = 1;
      listOfRequestors[kFDDMultZeqs].Append(Form("%s ", "dependency check"));
    }
    if (internalOpts.mEnabledTables[kCentMFTs] && !internalOpts.mEnabledTables[kMFTMults]) {
      internalOpts.mEnabledTables[kMFTMults] = 1;
      listOfRequestors[kMFTMults].Append(Form("%s ", "dependency check"));
    }
    if (internalOpts.mEnabledTables[kCentMFTs] && !internalOpts.mEnabledTables[kMFTMultZeqs]) {
      internalOpts.mEnabledTables[kMFTMultZeqs] = 1;
      listOfRequestors[kMFTMultZeqs].Append(Form("%s ", "dependency check"));
    }
    if (internalOpts.mEnabledTables[kCentNGlobals] && !internalOpts.mEnabledTables[kMultsGlobal]) {
      internalOpts.mEnabledTables[kMultsGlobal] = 1;
      listOfRequestors[kMultsGlobal].Append(Form("%s ", "dependency check"));
    }
    if (internalOpts.mEnabledTables[kCentNGlobals] && !internalOpts.mEnabledTables[kGlobalMultZeqs]) {
      internalOpts.mEnabledTables[kGlobalMultZeqs] = 1;
      listOfRequestors[kGlobalMultZeqs].Append(Form("%s ", "dependency check"));
    }
    if (internalOpts.embedINELgtZEROselection.value > 0 && !internalOpts.mEnabledTables[kPVMults]) {
      internalOpts.mEnabledTables[kPVMults] = 1;
      listOfRequestors[kPVMults].Append(Form("%s ", "dependency check"));
    }

    // list enabled tables
    for (int i = 0; i < nTablesConst; i++) {
      // printout to be improved in the future
      if (internalOpts.mEnabledTables[i]) {
        LOGF(info, " -~> Table enabled: %s, requested by %s", tableNames[i], listOfRequestors[i].Data());
      }
    }

    mRunNumber = 0;
    mRunNumberCentrality = 0;
    lCalibLoaded = false;
    hVtxZFV0A = nullptr;
    hVtxZFT0A = nullptr;
    hVtxZFT0C = nullptr;
    hVtxZFDDA = nullptr;
    hVtxZFDDC = nullptr;
    hVtxZNTracks = nullptr;
    hVtxZNMFTTracks = nullptr;
    hVtxZNGlobalTracks = nullptr;

    opts = internalOpts;
  }

  //__________________________________________________
  template <typename TCollision, typename TTracks, typename TBC, typename TOutputGroup>
  o2::common::multiplicity::multEntry collisionProcessRun2(TCollision const& collision, TTracks const& tracks, TBC const& bc, TOutputGroup& cursors)
  {
    // initialize properties
    o2::common::multiplicity::multEntry mults;

    mults.posZ = collision.posZ();
    mults.spdClustersL0 = bc.spdClustersL0();
    mults.spdClustersL1 = bc.spdClustersL1();
    //_______________________________________________________________________
    // forward detector signals, raw
    if (collision.has_fv0a()) {
      for (const auto& amplitude : collision.fv0a().amplitude()) {
        mults.multFV0A += amplitude;
      }
    }
    if (collision.has_fv0c()) {
      for (const auto& amplitude : collision.fv0c().amplitude()) {
        mults.multFV0C += amplitude;
      }
    }
    if (collision.has_ft0()) {
      auto ft0 = collision.ft0();
      for (const auto& amplitude : ft0.amplitudeA()) {
        mults.multFT0A += amplitude;
      }
      for (const auto& amplitude : ft0.amplitudeC()) {
        mults.multFT0C += amplitude;
      }
    }
    if (collision.has_zdc()) {
      auto zdc = collision.zdc();
      mults.multZNA = zdc.energyCommonZNA();
      mults.multZNC = zdc.energyCommonZNC();
    }

    //_______________________________________________________________________
    // determine if barrel track loop is required, do it (once!) if so but save CPU if not
    if (internalOpts.mEnabledTables[kPVMults] || internalOpts.mEnabledTables[kTPCMults] || internalOpts.mEnabledTables[kTrackletMults]) {
      // Try to do something Similar to https://github.com/alisw/AliPhysics/blob/22862a945004f719f8e9664c0264db46e7186a48/OADB/AliPPVsMultUtils.cxx#L541C26-L541C37
      for (const auto& track : tracks) {
        // check whether the track is a tracklet
        if (track.trackType() == o2::aod::track::Run2Tracklet) {
          if (internalOpts.mEnabledTables[kTrackletMults]) {
            mults.multTracklets++;
          }
          if (internalOpts.mEnabledTables[kPVMults]) {
            if (std::abs(track.eta()) < 1.0) {
              mults.multNContribsEta1++; // pvmults
              if (std::abs(track.eta()) < 0.8) {
                mults.multNContribs++; // pvmults
                if (std::abs(track.eta()) < 0.5) {
                  mults.multNContribsEtaHalf++; // pvmults
                }
              }
            }
          }
        }
        // check whether the track is a global ITS-TPC track
        if (track.tpcNClsFindable() > 0) {
          if (internalOpts.mEnabledTables[kTPCMults]) {
            mults.multTPC++;
          }
        }
      }
    }

    // fill standard cursors if required
    if (internalOpts.mEnabledTables[kFV0Mults]) {
      cursors.tableFV0(mults.multFV0A, mults.multFV0C);
    }
    if (internalOpts.mEnabledTables[kFT0Mults]) {
      cursors.tableFT0(mults.multFT0A, mults.multFT0C);
    }
    if (internalOpts.mEnabledTables[kFDDMults]) {
      cursors.tableFDD(mults.multFDDA, mults.multFDDC);
    }
    if (internalOpts.mEnabledTables[kZDCMults]) {
      cursors.tableZDC(mults.multZNA, mults.multZNC, 0.0f, 0.0f, 0.0f, 0.0f);
    }
    if (internalOpts.mEnabledTables[kTrackletMults]) { // Tracklets only Run2
      cursors.tableTracklet(mults.multTracklets);
    }
    if (internalOpts.mEnabledTables[kTPCMults]) {
      cursors.tableTpc(mults.multTPC);
    }
    if (internalOpts.mEnabledTables[kPVMults]) {
      cursors.tablePv(mults.multNContribs, mults.multNContribsEta1, mults.multNContribsEtaHalf);
    }

    return mults;
  }

  //__________________________________________________
  template <typename TCCDB, typename TMetadataInfo, typename TCollision, typename TTracks, typename TBC, typename TOutputGroup>
  o2::common::multiplicity::multEntry collisionProcessRun3(TCCDB const& ccdb, TMetadataInfo const& metadataInfo, TCollision const& collision, TTracks const& tracks, TBC const& bc, TOutputGroup& cursors)
  {
    // initialize properties
    o2::common::multiplicity::multEntry mults;

    //_______________________________________________________________________
    // preparatory steps
    if (internalOpts.doVertexZeq > 0) {
      if (bc.runNumber() != mRunNumber) {
        mRunNumber = bc.runNumber(); // mark this run as at least tried
        if (internalOpts.reconstructionPass.value == "") {
          lCalibObjects = ccdb->template getForRun<TList>(internalOpts.ccdbPathVtxZ, mRunNumber);
        } else if (internalOpts.reconstructionPass.value == "metadata") {
          std::map<std::string, std::string> metadata;
          metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
          LOGF(info, "Loading CCDB for reconstruction pass (from metadata): %s", metadataInfo.get("RecoPassName"));
          lCalibObjects = ccdb->template getSpecificForRun<TList>(internalOpts.ccdbPathVtxZ, mRunNumber, metadata);
        } else {
          std::map<std::string, std::string> metadata;
          metadata["RecoPassName"] = internalOpts.reconstructionPass.value;
          LOGF(info, "Loading CCDB for reconstruction pass (from provided argument): %s", internalOpts.reconstructionPass.value);
          lCalibObjects = ccdb->template getSpecificForRun<TList>(internalOpts.ccdbPathVtxZ, mRunNumber, metadata);
        }

        if (lCalibObjects) {
          hVtxZFV0A = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFV0A"));
          hVtxZFT0A = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFT0A"));
          hVtxZFT0C = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFT0C"));
          hVtxZFDDA = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFDDA"));
          hVtxZFDDC = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFDDC"));
          hVtxZNTracks = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZNTracksPV"));
          hVtxZNMFTTracks = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZMFT"));
          hVtxZNGlobalTracks = static_cast<TProfile*>(lCalibObjects->FindObject("hVtxZNGlobals"));
          lCalibLoaded = true;
          // Capture error
          if (!hVtxZFV0A || !hVtxZFT0A || !hVtxZFT0C || !hVtxZFDDA || !hVtxZFDDC || !hVtxZNTracks) {
            LOGF(error, "Problem loading CCDB objects! Please check");
            lCalibLoaded = false;
          }
          if (!hVtxZNMFTTracks) {
            LOGF(info, "MFT track counter: vertex-Z calibration not loaded, will run without.");
          }
          if (!hVtxZNGlobalTracks) {
            LOGF(info, "Global track counter: vertex-Z calibration not loaded, will run without.");
          }
        } else {
          LOGF(error, "Problem loading CCDB object! Please check");
          lCalibLoaded = false;
        }
      }
    }

    //_______________________________________________________________________
    // forward detector signals, raw
    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();
      for (size_t ii = 0; ii < fv0.amplitude().size(); ii++) {
        auto amplitude = fv0.amplitude()[ii];
        auto channel = fv0.channel()[ii];
        mults.multFV0A += amplitude;
        if (channel > 7) {
          mults.multFV0AOuter += amplitude;
        }
      }
    } else {
      mults.multFV0A = -999.f;
      mults.multFV0AOuter = -999.f;
    }
    if (collision.has_foundFT0()) {
      const auto& ft0 = collision.foundFT0();
      for (const auto& amplitude : ft0.amplitudeA()) {
        mults.multFT0A += amplitude;
      }
      for (const auto& amplitude : ft0.amplitudeC()) {
        mults.multFT0C += amplitude;
      }
    } else {
      mults.multFT0A = -999.f;
      mults.multFT0C = -999.f;
    }
    if (collision.has_foundFDD()) {
      const auto& fdd = collision.foundFDD();
      for (const auto& amplitude : fdd.chargeA()) {
        mults.multFDDA += amplitude;
      }
      for (const auto& amplitude : fdd.chargeC()) {
        mults.multFDDC += amplitude;
      }
    } else {
      mults.multFDDA = -999.f;
      mults.multFDDC = -999.f;
    }
    if (bc.has_zdc()) {
      mults.multZNA = bc.zdc().amplitudeZNA();
      mults.multZNC = bc.zdc().amplitudeZNC();
      mults.multZEM1 = bc.zdc().amplitudeZEM1();
      mults.multZEM2 = bc.zdc().amplitudeZEM2();
      mults.multZPA = bc.zdc().amplitudeZPA();
      mults.multZPC = bc.zdc().amplitudeZPC();
    } else {
      mults.multZNA = -999.f;
      mults.multZNC = -999.f;
      mults.multZEM1 = -999.f;
      mults.multZEM2 = -999.f;
      mults.multZPA = -999.f;
      mults.multZPC = -999.f;
    }

    // fill standard cursors if required
    if (internalOpts.mEnabledTables[kTrackletMults]) { // Tracklets (only Run2) nothing to do (to be removed!)
      cursors.tableTracklet(0);
    }
    if (internalOpts.mEnabledTables[kFV0Mults]) {
      cursors.tableFV0(mults.multFV0A, mults.multFV0C);
    }
    if (internalOpts.mEnabledTables[kFV0AOuterMults]) {
      cursors.tableFV0AOuter(mults.multFV0AOuter);
    }
    if (internalOpts.mEnabledTables[kFT0Mults]) {
      cursors.tableFT0(mults.multFT0A, mults.multFT0C);
    }
    if (internalOpts.mEnabledTables[kFDDMults]) {
      cursors.tableFDD(mults.multFDDA, mults.multFDDC);
    }
    if (internalOpts.mEnabledTables[kZDCMults]) {
      cursors.tableZDC(mults.multZNA, mults.multZNC, mults.multZEM1, mults.multZEM2, mults.multZPA, mults.multZPC);
    }

    //_______________________________________________________________________
    // fill selections (for posterior derived analysis if requested)
    if (internalOpts.mEnabledTables[kMultSelections]) {
      cursors.multSelections(collision.selection_raw());
    }

    //_______________________________________________________________________
    // vertex-Z equalized signals
    if (internalOpts.mEnabledTables[kFV0MultZeqs]) {
      if (mults.multFV0A > -1.0f && std::fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
        mults.multFV0AZeq = hVtxZFV0A->Interpolate(0.0) * mults.multFV0A / hVtxZFV0A->Interpolate(collision.posZ());
      } else {
        mults.multFV0AZeq = 0.0f;
      }
      cursors.tableFV0Zeqs(mults.multFV0AZeq);
    }
    if (internalOpts.mEnabledTables[kFT0MultZeqs]) {
      if (mults.multFT0A > -1.0f && std::fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
        mults.multFT0AZeq = hVtxZFT0A->Interpolate(0.0) * mults.multFT0A / hVtxZFT0A->Interpolate(collision.posZ());
      } else {
        mults.multFT0AZeq = 0.0f;
      }
      if (mults.multFT0C > -1.0f && std::fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
        mults.multFT0CZeq = hVtxZFT0C->Interpolate(0.0) * mults.multFT0C / hVtxZFT0C->Interpolate(collision.posZ());
      } else {
        mults.multFT0CZeq = 0.0f;
      }
      cursors.tableFT0Zeqs(mults.multFT0AZeq, mults.multFT0CZeq);
    }
    if (internalOpts.mEnabledTables[kFDDMultZeqs]) {
      if (mults.multFDDA > -1.0f && std::fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
        mults.multFDDAZeq = hVtxZFDDA->Interpolate(0.0) * mults.multFDDA / hVtxZFDDA->Interpolate(collision.posZ());
      } else {
        mults.multFDDAZeq = 0.0f;
      }
      if (mults.multFDDC > -1.0f && std::fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
        mults.multFDDCZeq = hVtxZFDDC->Interpolate(0.0) * mults.multFDDC / hVtxZFDDC->Interpolate(collision.posZ());
      } else {
        mults.multFDDCZeq = 0.0f;
      }
      cursors.tableFDDZeqs(mults.multFDDAZeq, mults.multFDDCZeq);
    }

    //_______________________________________________________________________
    // determine if barrel track loop is required, do it (once!) if so but save CPU if not
    if (internalOpts.mEnabledTables[kTPCMults] || internalOpts.mEnabledTables[kPVMults] || internalOpts.mEnabledTables[kMultsExtra] || internalOpts.mEnabledTables[kPVMultZeqs] || internalOpts.mEnabledTables[kMultsGlobal] || internalOpts.mEnabledTables[kGlobalMultZeqs]) {
      // single loop to calculate all
      for (const auto& track : tracks) {
        if (track.hasTPC()) {
          mults.multTPC++;
          if (track.hasITS()) {
            mults.multAllTracksITSTPC++; // multsextra
          } else {
            mults.multAllTracksTPCOnly++; // multsextra
          }
        }
        // PV contributor checked explicitly
        if (track.isPVContributor()) {
          if (std::abs(track.eta()) < 1.0) {
            mults.multNContribsEta1++; // pvmults
            if (std::abs(track.eta()) < 0.8) {
              mults.multNContribs++; // pvmults
              if (std::abs(track.eta()) < 0.5) {
                mults.multNContribsEtaHalf++; // pvmults
              }
            }
          }
          if (track.hasITS()) {
            mults.multHasITS++; // multsextra
            if (track.hasTPC())
              mults.multITSTPC++; // multsextra
            if (!track.hasTPC() && !track.hasTOF() && !track.hasTRD()) {
              mults.multITSOnly++; // multsextra
            }
          }
          if (track.hasTPC()) {
            mults.multHasTPC++; // multsextra
            if (!track.hasITS() && !track.hasTOF() && !track.hasTRD()) {
              mults.multTPCOnly++; // multsextra
            }
          }
          if (track.hasTOF()) {
            mults.multHasTOF++; // multsextra
          }
          if (track.hasTRD()) {
            mults.multHasTRD++; // multsextra
          }
        }

        // global counters: do them only in case information is provided in tracks table
        if constexpr (requires { track.isQualityTrack(); }) {
          if (track.pt() < internalOpts.maxPtGlobalTrack.value && track.pt() > internalOpts.minPtGlobalTrack.value && std::fabs(track.eta()) < 1.0f && track.isPVContributor() && track.isQualityTrack()) {
            if (track.itsNCls() < internalOpts.minNclsITSGlobalTrack || track.itsNClsInnerBarrel() < internalOpts.minNclsITSibGlobalTrack) {
              continue;
            }
            mults.multNbrContribsEta10GlobalTrackWoDCA++;

            if (std::abs(track.eta()) < 0.8) {
              mults.multNbrContribsEta08GlobalTrackWoDCA++;
            }
            if (std::abs(track.eta()) < 0.5) {
              mults.multNbrContribsEta05GlobalTrackWoDCA++;
            }
          }
          if (std::fabs(track.eta()) < 0.8 && track.tpcNClsFound() >= 80 && track.tpcNClsCrossedRows() >= 100) {
            if (track.isGlobalTrack()) {
              mults.multGlobalTracks++;
            }
          }
        } // end constexpr requires track selection stuff
      }

      cursors.multsGlobal(mults.multGlobalTracks, mults.multNbrContribsEta08GlobalTrackWoDCA, mults.multNbrContribsEta10GlobalTrackWoDCA, mults.multNbrContribsEta05GlobalTrackWoDCA);

      if (!hVtxZNGlobalTracks || std::fabs(collision.posZ()) > 15.0f) {
        mults.multGlobalTracksZeq = mults.multGlobalTracks; // if no equalization available, don't do it
      } else {
        mults.multGlobalTracksZeq = hVtxZNGlobalTracks->Interpolate(0.0) * mults.multGlobalTracks / hVtxZNGlobalTracks->Interpolate(collision.posZ());
      }

      // provide vertex-Z equalized Nglobals (or non-equalized if missing or beyond range)
      if (internalOpts.mEnabledTables[kGlobalMultZeqs]) {
        cursors.tableNGlobalZeqs(mults.multGlobalTracksZeq);
      }
    }

    // fill track counters at this stage if requested
    if (internalOpts.mEnabledTables[kTPCMults]) {
      cursors.tableTpc(mults.multTPC);
    }
    if (internalOpts.mEnabledTables[kPVMults]) {
      cursors.tablePv(mults.multNContribs, mults.multNContribsEta1, mults.multNContribsEtaHalf);
    }
    if (internalOpts.mEnabledTables[kMultsExtra]) {
      cursors.tableExtra(collision.numContrib(), collision.chi2(), collision.collisionTimeRes(),
                         bc.runNumber(), collision.posZ(), collision.sel8(),
                         mults.multHasITS, mults.multHasTPC, mults.multHasTOF, mults.multHasTRD,
                         mults.multITSOnly, mults.multTPCOnly, mults.multITSTPC,
                         mults.multAllTracksTPCOnly, mults.multAllTracksITSTPC,
                         collision.trackOccupancyInTimeRange(),
                         collision.ft0cOccupancyInTimeRange(),
                         collision.flags());
    }
    if (internalOpts.mEnabledTables[kPVMultZeqs]) {
      if (std::fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
        mults.multNContribsZeq = hVtxZNTracks->Interpolate(0.0) * mults.multNContribs / hVtxZNTracks->Interpolate(collision.posZ());
      } else {
        mults.multNContribsZeq = 0.0f;
      }
      cursors.tablePVZeqs(mults.multNContribsZeq);
    }

    // return multiplicity object such that it is handled properly when computing centrality
    return mults;
  }

  //__________________________________________________
  template <typename TMCCollision, typename TMCParticles, typename TPDGService, typename TOutputGroup>
  void collisionProcessMonteCarlo(TMCCollision const& mccollision, TMCParticles const& mcparticles, TPDGService const& pdg, TOutputGroup& cursors)
  {
    int multFT0A = 0;
    int multFV0A = 0;
    int multFT0C = 0;
    int multFDDA = 0;
    int multFDDC = 0;
    int multBarrelEta05 = 0;
    int multBarrelEta08 = 0;
    int multBarrelEta10 = 0;
    for (auto const& mcPart : mcparticles) {
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
    cursors.tableExtraMc(multFT0A, multFT0C, multFV0A, multFDDA, multFDDC, multBarrelEta05, multBarrelEta08, multBarrelEta10, mccollision.posZ());
  }

  //__________________________________________________
  template <typename TCollision, typename TMFTTracks, typename TBestCollisionsFwd, typename TMultBuffer, typename TOutputGroup>
  void collisionProcessMFT(TCollision const& collision, TMFTTracks const& mfttracks, TBestCollisionsFwd const& retracks, TMultBuffer& mults, TOutputGroup& cursors)
  {
    int nAllTracks = 0;
    int nTracks = 0;

    for (const auto& track : mfttracks) {
      if (track.nClusters() >= minNclsMFTTrack) {
        nAllTracks++;
      }
    }

    if (retracks.size() > 0) {
      for (const auto& retrack : retracks) {
        auto track = retrack.mfttrack();
        if (track.nClusters() < minNclsMFTTrack) {
          continue; // min cluster requirement
        }
        if (track.eta() > maxEtaMFTTrack || track.eta() < minEtaMFTTrack) {
          continue; // too far to be of true interest
        }
        if (std::abs(retrack.bestDCAXY()) > maxDCAxyToPVMFTTrack) {
          continue; // does not point to PV properly
        }
        nTracks++;
      }
    }
    cursors.mftMults(nAllTracks, nTracks);
    mults[collision.globalIndex()].multMFTAllTracks = nAllTracks;
    mults[collision.globalIndex()].multMFTTracks = nTracks;

    // vertex-Z equalized MFT
    if (!hVtxZNMFTTracks || std::fabs(collision.posZ()) > 15.0f) {
      mults[collision.globalIndex()].multMFTTracksZeq = mults[collision.globalIndex()].multMFTTracks; // if no equalization available, don't do it
    } else {
      mults[collision.globalIndex()].multMFTTracksZeq = hVtxZNMFTTracks->Interpolate(0.0) * mults[collision.globalIndex()].multMFTTracks / hVtxZNMFTTracks->Interpolate(collision.posZ());
    }

    // provide vertex-Z equalized Nglobals (or non-equalized if missing or beyond range)
    if (internalOpts.mEnabledTables[kMFTMultZeqs]) {
      cursors.tableNMFTZeqs(mults[collision.globalIndex()].multMFTTracksZeq);
    }
  }

  //__________________________________________________
  template <typename TCCDB, typename TMetadata, typename TBC>
  void ConfigureCentralityRun2(TCCDB& ccdb, TMetadata const& metadataInfo, TBC const& bc)
  {
    if (bc.runNumber() != mRunNumberCentrality) {
      mRunNumberCentrality = bc.runNumber(); // mark that this run has been attempted already regardless of outcome
      LOGF(info, "centrality loading procedure for timestamp=%llu, run number=%d", bc.timestamp(), bc.runNumber());
      TList* callst = nullptr;
      // Check if the ccdb path is a root file
      if (internalOpts.ccdbPathCentrality.value.find(".root") != std::string::npos) {
        TFile f(internalOpts.ccdbPathCentrality.value.c_str(), "READ");
        f.GetObject(internalOpts.reconstructionPass.value.c_str(), callst);
        if (!callst) {
          f.ls();
          LOG(fatal) << "No calibration list " << internalOpts.reconstructionPass.value << " found.";
        }
      } else {
        if (internalOpts.reconstructionPass.value == "") {
          callst = ccdb->template getForRun<TList>(internalOpts.ccdbPathCentrality, bc.runNumber());
        } else if (internalOpts.reconstructionPass.value == "metadata") {
          std::map<std::string, std::string> metadata;
          metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
          LOGF(info, "Loading CCDB for reconstruction pass (from metadata): %s", metadataInfo.get("RecoPassName"));
          callst = ccdb->template getSpecificForRun<TList>(internalOpts.ccdbPathCentrality, bc.runNumber(), metadata);
        } else {
          std::map<std::string, std::string> metadata;
          metadata["RecoPassName"] = internalOpts.reconstructionPass.value;
          LOGF(info, "Loading CCDB for reconstruction pass (from provided argument): %s", internalOpts.reconstructionPass.value);
          callst = ccdb->template getSpecificForRun<TList>(internalOpts.ccdbPathCentrality, bc.runNumber(), metadata);
        }
      }

      Run2V0MInfo.mCalibrationStored = false;
      Run2V0AInfo.mCalibrationStored = false;
      Run2SPDTksInfo.mCalibrationStored = false;
      Run2SPDClsInfo.mCalibrationStored = false;
      Run2CL0Info.mCalibrationStored = false;
      Run2CL1Info.mCalibrationStored = false;
      if (callst != nullptr) {
        auto getccdb = [callst](const char* ccdbhname) {
          TH1* h = reinterpret_cast<TH1*>(callst->FindObject(ccdbhname));
          return h;
        };
        auto getformulaccdb = [callst](const char* ccdbhname) {
          TFormula* f = reinterpret_cast<TFormula*>(callst->FindObject(ccdbhname));
          return f;
        };

        if (internalOpts.mEnabledTables[kCentRun2V0Ms]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2V0MInfo.mhVtxAmpCorrV0A = getccdb("hVtx_fAmplitude_V0A_Normalized");
          Run2V0MInfo.mhVtxAmpCorrV0C = getccdb("hVtx_fAmplitude_V0C_Normalized");
          Run2V0MInfo.mhMultSelCalib = getccdb("hMultSelCalib_V0M");
          Run2V0MInfo.mMCScale = getformulaccdb(TString::Format("%s-V0M", internalOpts.generatorName->c_str()).Data());
          if ((Run2V0MInfo.mhVtxAmpCorrV0A != nullptr) && (Run2V0MInfo.mhVtxAmpCorrV0C != nullptr) && (Run2V0MInfo.mhMultSelCalib != nullptr)) {
            if (internalOpts.generatorName->length() != 0) {
              if (Run2V0MInfo.mMCScale != nullptr) {
                for (int ixpar = 0; ixpar < 6; ++ixpar) {
                  Run2V0MInfo.mMCScalePars[ixpar] = Run2V0MInfo.mMCScale->GetParameter(ixpar);
                }
              } else {
                // continue filling with non-valid values (105)
                LOGF(info, "MC Scale information from V0M for run %d not available", bc.runNumber());
              }
            }
            Run2V0MInfo.mCalibrationStored = true;
          } else {
            // continue filling with non-valid values (105)
            LOGF(info, "Calibration information from V0M for run %d corrupted, will fill V0M tables with dummy values", bc.runNumber());
          }
        }
        if (internalOpts.mEnabledTables[kCentRun2V0As]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2V0AInfo.mhVtxAmpCorrV0A = getccdb("hVtx_fAmplitude_V0A_Normalized");
          Run2V0AInfo.mhMultSelCalib = getccdb("hMultSelCalib_V0A");
          if ((Run2V0AInfo.mhVtxAmpCorrV0A != nullptr) && (Run2V0AInfo.mhMultSelCalib != nullptr)) {
            Run2V0AInfo.mCalibrationStored = true;
          } else {
            // continue filling with non-valid values (105)
            LOGF(info, "Calibration information from V0A for run %d corrupted, will fill V0A tables with dummy values", bc.runNumber());
          }
        }
        if (internalOpts.mEnabledTables[kCentRun2SPDTrks]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2SPDTksInfo.mhVtxAmpCorr = getccdb("hVtx_fnTracklets_Normalized");
          Run2SPDTksInfo.mhMultSelCalib = getccdb("hMultSelCalib_SPDTracklets");
          if ((Run2SPDTksInfo.mhVtxAmpCorr != nullptr) && (Run2SPDTksInfo.mhMultSelCalib != nullptr)) {
            Run2SPDTksInfo.mCalibrationStored = true;
          } else {
            // continue filling with non-valid values (105)
            LOGF(info, "Calibration information from SPD tracklets for run %d corrupted, will fill SPD tracklets tables with dummy values", bc.runNumber());
          }
        }
        if (internalOpts.mEnabledTables[kCentRun2SPDClss]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2SPDClsInfo.mhVtxAmpCorrCL0 = getccdb("hVtx_fnSPDClusters0_Normalized");
          Run2SPDClsInfo.mhVtxAmpCorrCL1 = getccdb("hVtx_fnSPDClusters1_Normalized");
          Run2SPDClsInfo.mhMultSelCalib = getccdb("hMultSelCalib_SPDClusters");
          if ((Run2SPDClsInfo.mhVtxAmpCorrCL0 != nullptr) && (Run2SPDClsInfo.mhVtxAmpCorrCL1 != nullptr) && (Run2SPDClsInfo.mhMultSelCalib != nullptr)) {
            Run2SPDClsInfo.mCalibrationStored = true;
          } else {
            // continue filling with non-valid values (105)
            LOGF(info, "Calibration information from SPD clusters for run %d corrupted, will fill SPD clusters tables with dummy values", bc.runNumber());
          }
        }
        if (internalOpts.mEnabledTables[kCentRun2CL0s]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2CL0Info.mhVtxAmpCorr = getccdb("hVtx_fnSPDClusters0_Normalized");
          Run2CL0Info.mhMultSelCalib = getccdb("hMultSelCalib_CL0");
          if ((Run2CL0Info.mhVtxAmpCorr != nullptr) && (Run2CL0Info.mhMultSelCalib != nullptr)) {
            Run2CL0Info.mCalibrationStored = true;
          } else {
            // continue filling with non-valid values (105)
            LOGF(info, "Calibration information from CL0 multiplicity for run %d corrupted, will fill CL0 multiplicity tables with dummy values", bc.runNumber());
          }
        }
        if (internalOpts.mEnabledTables[kCentRun2CL1s]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2CL1Info.mhVtxAmpCorr = getccdb("hVtx_fnSPDClusters1_Normalized");
          Run2CL1Info.mhMultSelCalib = getccdb("hMultSelCalib_CL1");
          if ((Run2CL1Info.mhVtxAmpCorr != nullptr) && (Run2CL1Info.mhMultSelCalib != nullptr)) {
            Run2CL1Info.mCalibrationStored = true;
          } else {
            // continue filling with non-valid values (105)
            LOGF(info, "Calibration information from CL1 multiplicity for run %d corrupted, will fill CL1 multiplicity tables with dummy values", bc.runNumber());
          }
        }
      } else {
        LOGF(info, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu, will fill tables with dummy values", bc.runNumber(), bc.timestamp());
      }
    }
  }

  //__________________________________________________
  template <typename TCCDB, typename TMetadata, typename TBC>
  void ConfigureCentralityRun3(TCCDB& ccdb, TMetadata const& metadataInfo, TBC const& bc)
  {
    if (bc.runNumber() != mRunNumberCentrality) {
      mRunNumberCentrality = bc.runNumber(); // mark that this run has been attempted already regardless of outcome
      LOGF(info, "centrality loading procedure for timestamp=%llu, run number=%d", bc.timestamp(), bc.runNumber());

      // capture the need for PYTHIA calibration in Pb-Pb runs
      if (metadataInfo.isMC() && mRunNumber >= 544013 && mRunNumber <= 545367) {
        LOGF(info, "This is MC for Pb-Pb. Setting generatorName automatically to PYTHIA");
        internalOpts.generatorName.value = "PYTHIA";
      }

      // capture the need for PYTHIA calibration in light ion runs automatically
      if (metadataInfo.isMC() && mRunNumber >= 564250 && mRunNumber <= 564472) {
        LOGF(info, "This is MC for light ion runs. Setting generatorName automatically to PYTHIA");
        internalOpts.generatorName.value = "PYTHIA";
      }

      LOGF(info, "centrality loading procedure for timestamp=%llu, run number=%d", bc.timestamp(), bc.runNumber());
      TList* callst = nullptr;
      // Check if the ccdb path is a root file
      if (internalOpts.ccdbPathCentrality.value.find(".root") != std::string::npos) {
        TFile f(internalOpts.ccdbPathCentrality.value.c_str(), "READ");
        f.GetObject(internalOpts.reconstructionPass.value.c_str(), callst);
        if (!callst) {
          f.ls();
          LOG(fatal) << "No calibration list " << internalOpts.reconstructionPass.value << " found.";
        }
      } else {
        if (internalOpts.reconstructionPass.value == "") {
          callst = ccdb->template getForRun<TList>(internalOpts.ccdbPathCentrality, bc.runNumber());
        } else if (internalOpts.reconstructionPass.value == "metadata") {
          std::map<std::string, std::string> metadata;
          metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
          LOGF(info, "Loading CCDB for reconstruction pass (from metadata): %s", metadataInfo.get("RecoPassName"));
          callst = ccdb->template getSpecificForRun<TList>(internalOpts.ccdbPathCentrality, bc.runNumber(), metadata);
        } else {
          std::map<std::string, std::string> metadata;
          metadata["RecoPassName"] = internalOpts.reconstructionPass.value;
          LOGF(info, "Loading CCDB for reconstruction pass (from provided argument): %s", internalOpts.reconstructionPass.value);
          callst = ccdb->template getSpecificForRun<TList>(internalOpts.ccdbPathCentrality, bc.runNumber(), metadata);
        }
      }

      fv0aInfo.mCalibrationStored = false;
      ft0mInfo.mCalibrationStored = false;
      ft0aInfo.mCalibrationStored = false;
      ft0cInfo.mCalibrationStored = false;
      ft0cVariant1Info.mCalibrationStored = false;
      ft0cVariant2Info.mCalibrationStored = false;
      fddmInfo.mCalibrationStored = false;
      ntpvInfo.mCalibrationStored = false;
      nGlobalInfo.mCalibrationStored = false;
      mftInfo.mCalibrationStored = false;
      if (callst != nullptr) {
        LOGF(info, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
        auto getccdb = [callst, bc](struct CalibrationInfo& estimator, const o2::framework::Configurable<std::string> generatorName) { // TODO: to consider the name inside the estimator structure
          estimator.mhMultSelCalib = reinterpret_cast<TH1*>(callst->FindObject(TString::Format("hCalibZeq%s", estimator.name.c_str()).Data()));
          estimator.mMCScale = reinterpret_cast<TFormula*>(callst->FindObject(TString::Format("%s-%s", generatorName->c_str(), estimator.name.c_str()).Data()));
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

        // invoke loading only for requested centralities
        if (internalOpts.mEnabledTables[kCentFV0As])
          getccdb(fv0aInfo, internalOpts.generatorName);
        if (internalOpts.mEnabledTables[kCentFT0Ms])
          getccdb(ft0mInfo, internalOpts.generatorName);
        if (internalOpts.mEnabledTables[kCentFT0As])
          getccdb(ft0aInfo, internalOpts.generatorName);
        if (internalOpts.mEnabledTables[kCentFT0Cs])
          getccdb(ft0cInfo, internalOpts.generatorName);
        if (internalOpts.mEnabledTables[kCentFT0CVariant1s])
          getccdb(ft0cVariant1Info, internalOpts.generatorName);
        if (internalOpts.mEnabledTables[kCentFT0CVariant2s])
          getccdb(ft0cVariant2Info, internalOpts.generatorName);
        if (internalOpts.mEnabledTables[kCentFDDMs])
          getccdb(fddmInfo, internalOpts.generatorName);
        if (internalOpts.mEnabledTables[kCentNTPVs])
          getccdb(ntpvInfo, internalOpts.generatorName);
        if (internalOpts.mEnabledTables[kCentNGlobals])
          getccdb(nGlobalInfo, internalOpts.generatorName);
        if (internalOpts.mEnabledTables[kCentMFTs])
          getccdb(mftInfo, internalOpts.generatorName);
      } else {
        LOGF(info, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu, will fill tables with dummy values", bc.runNumber(), bc.timestamp());
      }
    }
  }

  //__________________________________________________
  template <typename TCCDB, typename TMetadata, typename TBCs, typename TMultBuffer, typename TOutputGroup>
  void generateCentralitiesRun3(TCCDB& ccdb, TMetadata const& metadataInfo, TBCs const& bcs, TMultBuffer const& mults, TOutputGroup& cursors)
  {
    // takes multiplicity buffer and generates the desirable centrality values (if any)

    // first step: did someone actually ask for it? Otherwise, go home
    if (
      internalOpts.mEnabledTables[kCentFV0As] || internalOpts.mEnabledTables[kCentFT0Ms] ||
      internalOpts.mEnabledTables[kCentFT0As] || internalOpts.mEnabledTables[kCentFT0Cs] ||
      internalOpts.mEnabledTables[kCentFT0CVariant1s] ||
      internalOpts.mEnabledTables[kCentFT0CVariant2s] ||
      internalOpts.mEnabledTables[kCentFDDMs] ||
      internalOpts.mEnabledTables[kCentNTPVs] || internalOpts.mEnabledTables[kCentNGlobals] ||
      internalOpts.mEnabledTables[kCentMFTs] || internalOpts.mEnabledTables[kBCCentFT0Ms] ||
      internalOpts.mEnabledTables[kBCCentFT0As] || internalOpts.mEnabledTables[kBCCentFT0Cs]) {
      // check and update centrality calibration objects for Run 3
      const auto& firstbc = bcs.begin();
      ConfigureCentralityRun3(ccdb, metadataInfo, firstbc);

      /************************************************************
       * @brief Populates a table with data based on the given calibration information and multiplicity.
       *
       * @param table The table to populate.
       * @param estimator The calibration information.
       * @param multiplicity The multiplicity value.
       *************************************************************/

      auto populateTable = [&](auto& table, struct CalibrationInfo& estimator, float multiplicity, bool isInelGt0) {
        const bool assignOutOfRange = internalOpts.embedINELgtZEROselection && !isInelGt0;
        auto scaleMC = [](float x, const float pars[6]) {
          float core = ((pars[0] + pars[1] * std::pow(x, pars[2])) - pars[3]) / pars[4];
          if (core < 0.0f) {
            return 0.0f; // this should be marked as low multiplicity and not mapped, core^pars[5] would be NaN
          }
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

      // populate centralities per event
      for (size_t iEv = 0; iEv < mults.size(); iEv++) {
        bool isInelGt0 = (mults[iEv].multNContribsEta1 > 0);
        if (internalOpts.mEnabledTables[kCentFV0As])
          populateTable(cursors.centFV0A, fv0aInfo, mults[iEv].multFV0AZeq, isInelGt0);
        if (internalOpts.mEnabledTables[kCentFT0Ms])
          populateTable(cursors.centFT0M, ft0mInfo, mults[iEv].multFT0AZeq + mults[iEv].multFT0CZeq, isInelGt0);
        if (internalOpts.mEnabledTables[kCentFT0As])
          populateTable(cursors.centFT0A, ft0aInfo, mults[iEv].multFT0AZeq, isInelGt0);
        if (internalOpts.mEnabledTables[kCentFT0Cs])
          populateTable(cursors.centFT0C, ft0cInfo, mults[iEv].multFT0CZeq, isInelGt0);
        if (internalOpts.mEnabledTables[kCentFT0CVariant1s])
          populateTable(cursors.centFT0CVariant1, ft0cVariant1Info, mults[iEv].multFT0CZeq, isInelGt0);
        if (internalOpts.mEnabledTables[kCentFT0CVariant2s])
          populateTable(cursors.centFT0CVariant2, ft0cVariant2Info, mults[iEv].multFT0CZeq, isInelGt0);
        if (internalOpts.mEnabledTables[kCentFDDMs])
          populateTable(cursors.centFDDM, fddmInfo, mults[iEv].multFDDAZeq + mults[iEv].multFDDCZeq, isInelGt0);
        if (internalOpts.mEnabledTables[kCentNTPVs])
          populateTable(cursors.centNTPV, ntpvInfo, mults[iEv].multNContribs, isInelGt0);
        if (internalOpts.mEnabledTables[kCentNGlobals])
          populateTable(cursors.centNGlobals, nGlobalInfo, mults[iEv].multGlobalTracksZeq, isInelGt0);
        if (internalOpts.mEnabledTables[kCentMFTs])
          populateTable(cursors.centMFTs, mftInfo, mults[iEv].multMFTTracksZeq, isInelGt0);
      }

      // populate centralities per BC
      for (size_t ibc = 0; ibc < static_cast<size_t>(bcs.size()); ibc++) {
        float bcMultFT0A = 0;
        float bcMultFT0C = 0;

        const auto& bc = bcs.rawIteratorAt(ibc);
        if (bc.has_foundFT0()) {
          const auto& ft0 = bc.foundFT0();
          for (const auto& amplitude : ft0.amplitudeA()) {
            bcMultFT0A += amplitude;
          }
          for (const auto& amplitude : ft0.amplitudeC()) {
            bcMultFT0C += amplitude;
          }
        } else {
          bcMultFT0A = -999.f;
          bcMultFT0C = -999.f;
        }

        if (internalOpts.mEnabledTables[kBCCentFT0Ms])
          populateTable(cursors.bcCentFT0M, ft0mInfo, bcMultFT0A + bcMultFT0C, true);
        if (internalOpts.mEnabledTables[kBCCentFT0As])
          populateTable(cursors.bcCentFT0A, ft0aInfo, bcMultFT0A, true);
        if (internalOpts.mEnabledTables[kBCCentFT0Cs])
          populateTable(cursors.bcCentFT0C, ft0cInfo, bcMultFT0C, true);
      }
    }
  }
  //__________________________________________________
  template <typename TCCDB, typename TMetadata, typename TBCs, typename TMultBuffer, typename TOutputGroup>
  void generateCentralitiesRun2(TCCDB& ccdb, TMetadata const& metadataInfo, TBCs const& bcs, TMultBuffer const& mults, TOutputGroup& cursors)
  {
    // takes multiplicity buffer and generates the desirable centrality values (if any)
    // For Run 2
    if (
      internalOpts.mEnabledTables[kCentRun2V0Ms] || internalOpts.mEnabledTables[kCentRun2V0As] ||
      internalOpts.mEnabledTables[kCentRun2SPDTrks] || internalOpts.mEnabledTables[kCentRun2SPDClss] ||
      internalOpts.mEnabledTables[kCentRun2CL0s] || internalOpts.mEnabledTables[kCentRun2CL1s]) {
      // check and update centrality calibration objects for Run 3
      const auto& firstbc = bcs.begin();
      ConfigureCentralityRun2(ccdb, metadataInfo, firstbc);

      auto scaleMC = [](float x, const float pars[6]) {
        float core = ((pars[0] + pars[1] * std::pow(x, pars[2])) - pars[3]) / pars[4];
        if (core < 0.0f) {
          return 0.0f; // this should be marked as low multiplicity and not mapped, core^pars[5] would be NaN
        }
        return std::pow(((pars[0] + pars[1] * std::pow(x, pars[2])) - pars[3]) / pars[4], 1.0f / pars[5]);
      };

      // populate centralities per event
      for (size_t iEv = 0; iEv < mults.size(); iEv++) {
        if (internalOpts.mEnabledTables[kCentRun2V0Ms]) {
          float cV0M = 105.0f;
          if (Run2V0MInfo.mCalibrationStored) {
            float v0m;
            if (Run2V0MInfo.mMCScale != nullptr) {
              v0m = scaleMC(mults[iEv].multFV0A + mults[iEv].multFV0C, Run2V0MInfo.mMCScalePars);
              LOGF(debug, "Unscaled v0m: %f, scaled v0m: %f", mults[iEv].multFV0A + mults[iEv].multFV0C, v0m);
            } else {
              v0m = mults[iEv].multFV0A * Run2V0MInfo.mhVtxAmpCorrV0A->GetBinContent(Run2V0MInfo.mhVtxAmpCorrV0A->FindFixBin(mults[iEv].posZ)) +
                    mults[iEv].multFV0C * Run2V0MInfo.mhVtxAmpCorrV0C->GetBinContent(Run2V0MInfo.mhVtxAmpCorrV0C->FindFixBin(mults[iEv].posZ));
            }
            cV0M = Run2V0MInfo.mhMultSelCalib->GetBinContent(Run2V0MInfo.mhMultSelCalib->FindFixBin(v0m));
          }
          LOGF(debug, "centRun2V0M=%.0f", cV0M);
          // fill centrality columns
          cursors.centRun2V0M(cV0M);
        }
        if (internalOpts.mEnabledTables[kCentRun2V0As]) {
          float cV0A = 105.0f;
          if (Run2V0AInfo.mCalibrationStored) {
            float v0a = mults[iEv].multFV0A * Run2V0AInfo.mhVtxAmpCorrV0A->GetBinContent(Run2V0AInfo.mhVtxAmpCorrV0A->FindFixBin(mults[iEv].posZ));
            cV0A = Run2V0AInfo.mhMultSelCalib->GetBinContent(Run2V0AInfo.mhMultSelCalib->FindFixBin(v0a));
          }
          LOGF(debug, "centRun2V0A=%.0f", cV0A);
          // fill centrality columns
          cursors.centRun2V0A(cV0A);
        }
        if (internalOpts.mEnabledTables[kCentRun2SPDTrks]) {
          float cSPD = 105.0f;
          if (Run2SPDTksInfo.mCalibrationStored) {
            float spdm = mults[iEv].multTracklets * Run2SPDTksInfo.mhVtxAmpCorr->GetBinContent(Run2SPDTksInfo.mhVtxAmpCorr->FindFixBin(mults[iEv].posZ));
            cSPD = Run2SPDTksInfo.mhMultSelCalib->GetBinContent(Run2SPDTksInfo.mhMultSelCalib->FindFixBin(spdm));
          }
          LOGF(debug, "centSPDTracklets=%.0f", cSPD);
          cursors.centRun2SPDTracklets(cSPD);
        }
        if (internalOpts.mEnabledTables[kCentRun2SPDClss]) {
          float cSPD = 105.0f;
          if (Run2SPDClsInfo.mCalibrationStored) {
            float spdm = mults[iEv].spdClustersL0 * Run2SPDClsInfo.mhVtxAmpCorrCL0->GetBinContent(Run2SPDClsInfo.mhVtxAmpCorrCL0->FindFixBin(mults[iEv].posZ)) +
                         mults[iEv].spdClustersL1 * Run2SPDClsInfo.mhVtxAmpCorrCL1->GetBinContent(Run2SPDClsInfo.mhVtxAmpCorrCL1->FindFixBin(mults[iEv].posZ));
            cSPD = Run2SPDClsInfo.mhMultSelCalib->GetBinContent(Run2SPDClsInfo.mhMultSelCalib->FindFixBin(spdm));
          }
          LOGF(debug, "centSPDClusters=%.0f", cSPD);
          cursors.centRun2SPDClusters(cSPD);
        }
        if (internalOpts.mEnabledTables[kCentRun2CL0s]) {
          float cCL0 = 105.0f;
          if (Run2CL0Info.mCalibrationStored) {
            float cl0m = mults[iEv].spdClustersL0 * Run2CL0Info.mhVtxAmpCorr->GetBinContent(Run2CL0Info.mhVtxAmpCorr->FindFixBin(mults[iEv].posZ));
            cCL0 = Run2CL0Info.mhMultSelCalib->GetBinContent(Run2CL0Info.mhMultSelCalib->FindFixBin(cl0m));
          }
          LOGF(debug, "centCL0=%.0f", cCL0);
          cursors.centRun2CL0(cCL0);
        }
        if (internalOpts.mEnabledTables[kCentRun2CL1s]) {
          float cCL1 = 105.0f;
          if (Run2CL1Info.mCalibrationStored) {
            float cl1m = mults[iEv].spdClustersL1 * Run2CL1Info.mhVtxAmpCorr->GetBinContent(Run2CL1Info.mhVtxAmpCorr->FindFixBin(mults[iEv].posZ));
            cCL1 = Run2CL1Info.mhMultSelCalib->GetBinContent(Run2CL1Info.mhMultSelCalib->FindFixBin(cl1m));
          }
          LOGF(debug, "centCL1=%.0f", cCL1);
          cursors.centRun2CL1(cCL1);
        }
      }
    }
  }
}; // end BuilderModule

} // namespace multiplicity
} // namespace common
} // namespace o2

#endif // COMMON_TOOLS_MULTIPLICITY_MULTMODULE_H_
