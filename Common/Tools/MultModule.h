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

#ifndef COMMON_TOOLS_MULTMODULE_H_
#define COMMON_TOOLS_MULTMODULE_H_

#include <memory>
#include <cstdlib>
#include <cmath>
#include <array>
#include <string>
#include "Framework/AnalysisDataModel.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramSpec.h"
#include "TableHelper.h"
#include "Common/Core/TPCVDriftManager.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

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
  "MultMCExtras",
  "kMult2MCExtras",
  "kMFTMults",
  "kMultsGlobal",

  //centrality subcomponent
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
  "CentFDDMs",
  "CentNTPVs",
  "CentNGlobals",
  "CentMFTs"};

static constexpr int nTablesConst = 32;

static const std::vector<std::string> parameterNames{"enable"};
static const int defaultParameters[nTablesConst][nParameters]{
  {-1}, {-1}, {-1}, {-1}, {-1},
  {-1}, {-1}, {-1}, {-1}, {-1},
  {-1}, {-1}, {-1}, {-1}, {-1},
  {-1}, {-1}, {-1}, {-1}, {-1},
  {-1}, {-1}, {-1}, {-1}, {-1},
  {-1}, {-1}, {-1}, {-1}, {-1},
  {-1}, {-1}}; 

// table index : match order above
enum tableIndex {   kFV0Mults,        // standard
                    kFV0AOuterMults,  // standard
                    kFT0Mults,        // standard 
                    kFDDMults,        // standard
                    kZDCMults,        // standard 
                    kTrackletMults,   // Run 2 
                    kTPCMults,        // standard
                    kPVMults,         // standard
                    kMultsExtra,      // standard
                    kMultSelections,  // event selection
                    kFV0MultZeqs,     // zeq calib, standard
                    kFT0MultZeqs,     // zeq calib, standard
                    kFDDMultZeqs,     // zeq calib, standard
                    kPVMultZeqs,      // zeq calib, standard
                    kMultMCExtras,    // MC exclusive
                    kMult2MCExtras,   // MC exclusive
                    kMFTMults,        // requires MFT task
                    kMultsGlobal,     // requires track selection task

                    //centrality subcomponent
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
                    kCentFDDMs,         // standard Run 3
                    kCentNTPVs,         // standard Run 3
                    kCentNGlobals,      // requires track selection task 
                    kCentMFTs,          // requires MFT task
                    kNTables};

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
  o2::framework::Produces<aod::CentFDDMs> centFDDM;
  o2::framework::Produces<aod::CentNTPVs> centNTPV;
  o2::framework::Produces<aod::CentNGlobals> centNGlobals;
  o2::framework::Produces<aod::CentMFTs> centMFTs;

  //__________________________________________________
  // centrality tables per BC
  // FIXME - future development
};

// for providing temporary buffer 
// FIXME ideally cursors could be readable 
// to avoid duplicate memory allocation but ok
struct multEntry {
  float multFV0A = -999.0f;
  float multFV0C = -999.0f;
  float multFV0AOuter = -999.0f;
  float multFT0A = -999.0f;
  float multFT0C = -999.0f;
  float multFDDA = -999.0f;
  float multFDDC = -999.0f;
  float multZNA = -999.0f;
  float multZNC = -999.0f;
  float multZEM1 = -999.0f;
  float multZEM2 = -999.0f;
  float multZPA = -999.0f;
  float multZPC = -999.0f;
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

  int multGlobalTracks = 0; // multsGlobal                 
  int multNbrContribsEta05GlobalTrackWoDCA = 0; // multsGlobal
  int multNbrContribsEta08GlobalTrackWoDCA = 0; // multsGlobal
  int multNbrContribsEta10GlobalTrackWoDCA = 0; // multsGlobal

  int multMFTAllTracks = 0; // mft
  int multMFTTracks = 0; // mft
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

  // ccdb information
  o2::framework::Configurable<std::string> ccdbPathVtxZ{"ccdbPathVtxZ", "Centrality/Calibration", "The CCDB path for centrality/multiplicity information"};
  o2::framework::Configurable<std::string> reconstructionPass{"reconstructionPass", "", {"Apass to use when fetching the calibration tables. Empty (default) does not check for any pass. Use `metadata` to fetch it from the AO2D metadata. Otherwise it will override the metadata."}};

};

class MultModule
{
 public:
  MultModule()
  {
    // constructor
    mRunNumber = 0;
    lCalibLoaded = false;
    lCalibObjects = nullptr;
    hVtxZFV0A = nullptr;
    hVtxZFT0A = nullptr;
    hVtxZFT0C = nullptr;
    hVtxZFDDA = nullptr;
    hVtxZFDDC = nullptr;
    hVtxZNTracks = nullptr;
  }

  // internal: calib related, vtx-z profiles
  int mRunNumber;
  bool lCalibLoaded;
  TList* lCalibObjects;
  TProfile* hVtxZFV0A;
  TProfile* hVtxZFT0A;
  TProfile* hVtxZFT0C;
  TProfile* hVtxZFDDA;
  TProfile* hVtxZFDDC;
  TProfile* hVtxZNTracks;

  // declaration of structs here
  // (N.B.: will be invisible to the outside, create your own copies)
  o2::common::multiplicity::standardConfigurables internalOpts;

  template <typename TConfigurables, typename TInitContext>
  void init(TConfigurables const& opts, TInitContext& context)
  {
    // read in configurations from the task where it's used
    internalOpts = opts;

    mRunNumber = 0;

    internalOpts.mEnabledTables.resize(nTablesConst, 0);

    LOGF(info, "Configuring tables to generate");
    auto& workflows = context.services().template get<o2::framework::RunningWorkflowInfo const>();

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

    // list enabled tables
    for (int i = 0; i < nTablesConst; i++) {
      // printout to be improved in the future
      if (internalOpts.mEnabledTables[i]) {
        LOGF(info, " -~> Table enabled: %s, requested by %s", tableNames[i], listOfRequestors[i].Data());
      }
    }

    mRunNumber = 0;
    lCalibLoaded = false;
    hVtxZFV0A = nullptr;
    hVtxZFT0A = nullptr;
    hVtxZFT0C = nullptr;
    hVtxZFDDA = nullptr;
    hVtxZFDDC = nullptr;
    hVtxZNTracks = nullptr;
  }

  //__________________________________________________
  template <typename TCollision, typename TTracks, typename TBCs, typename TZdc, typename TFV0A, typename TFV0C, typename TFT0>
  void collisionProcessRun2(TCollision const& collision, TTracks const& tracks, TBCs const& bcs, TZdc const& zdc, TFV0A const& fv0a, TFV0C const& fv0c, TFT0 const& ft0 )
  {
    // initialize properties
    o2::common::multiplicity::multEntry mults; 
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
    }
    if (collision.has_foundFT0()) {
      const auto& ft0 = collision.foundFT0();
      for (const auto& amplitude : ft0.amplitudeA()) {
        mults.multFT0A += amplitude;
      }
      for (const auto& amplitude : ft0.amplitudeC()) {
        mults.multFT0C += amplitude;
      }
    }
    if (collision.has_foundFDD()) {
      const auto& fdd = collision.foundFDD();
      for (const auto& amplitude : fdd.chargeA()) {
        mults.multFDDA += amplitude;
      }
      for (const auto& amplitude : fdd.chargeC()) {
        mults.multFDDC += amplitude;
      }
    }
    if (bc.has_zdc()) {
      mults.multZNA = bc.zdc().amplitudeZNA();
      mults.multZNC = bc.zdc().amplitudeZNC();
      mults.multZEM1 = bc.zdc().amplitudeZEM1();
      mults.multZEM2 = bc.zdc().amplitudeZEM2();
      mults.multZPA = bc.zdc().amplitudeZPA();
      mults.multZPC = bc.zdc().amplitudeZPC();
    } 

    // fill standard cursors if required
    if(internalOpts.mEnabledTables[kFV0Mults]){
      cursors.tableFV0(mults.multFV0A, mults.multFV0C);
    }
    if(internalOpts.mEnabledTables[kFV0AOuterMults]){
      cursors.tableFV0AOuter(mults.multFV0AOuter);
    }
    if(internalOpts.mEnabledTables[kFT0Mults]){
      cursors.tableFT0(mults.multFT0A, mults.multFT0C);
    }
    if(internalOpts.mEnabledTables[kFDDMults]){
      cursors.tableFDD(mults.multFDDA, mults.multFDDC);
    }
    if(internalOpts.mEnabledTables[kZDCMults]){
      cursors.tableZDC(mults.multZNA, mults.multZNC, mults.multZEM1, mults.multZEM2, mults.multZPA, mults.multZPC);
    }

    //_______________________________________________________________________
    // forward detector signals, vertex-Z equalized
    if(internalOpts.mEnabledTables[kFV0MultZeqs]){
      if(std::fabs(collision.posZ() && lCalibLoaded)){
        mults.multFV0AZeq = hVtxZFV0A->Interpolate(0.0) * mults.multFV0A / hVtxZFV0A->Interpolate(collision.posZ());
      }else{
        mults.multFV0AZeq = 0.0f;
      }
      cursors.tableFV0Zeqs(mults.multFV0AZeq);
    }
    if(internalOpts.mEnabledTables[kFT0MultZeqs]){
      if(std::fabs(collision.posZ() && lCalibLoaded)){
        mults.multFT0AZeq = hVtxZFT0A->Interpolate(0.0) * mults.multFT0A / hVtxZFT0A->Interpolate(collision.posZ());
        mults.multFT0CZeq = hVtxZFT0C->Interpolate(0.0) * mults.multFT0C / hVtxZFT0C->Interpolate(collision.posZ());
      }else{
        mults.multFT0AZeq = 0.0f;
        mults.multFT0CZeq = 0.0f;
      }
      cursors.tableFT0Zeqs(mults.multFT0AZeq, mults.multFT0CZeq);
    }
    if(internalOpts.mEnabledTables[kFDDMultZeqs]){
      if(std::fabs(collision.posZ() && lCalibLoaded)){
        mults.multFDDAZeq = hVtxZFDDA->Interpolate(0.0) * mults.multFDDA / hVtxZFDDA->Interpolate(collision.posZ());
        mults.multFDDCZeq = hVtxZFDDC->Interpolate(0.0) * mults.multFDDC / hVtxZFDDC->Interpolate(collision.posZ());
      }else{
        mults.multFDDAZeq = 0.0f;
        mults.multFDDCZeq = 0.0f;
      }
      cursors.tableFDDZeqs(mults.multFDDAZeq, mults.multFDDCZeq);
    }

    //_______________________________________________________________________
    // determine if barrel track loop is required, do it (once!) if so but save CPU if not
    if(internalOpts.mEnabledTables[kTPCMults] || internalOpts.mEnabledTables[kPVMults] || internalOpts.mEnabledTables[kMultsExtra] || internalOpts.mEnabledTables[kPVMultZeqs] || internalOpts.mEnabledTables[kMultsGlobal]){
      // single loop to calculate all
      for (const auto& track : tracks) {
        if(track.hasTPC()){ 
          mults.multTPC++;
          if(track.hasITS()){
            mults.multAllTracksITSTPC++; // multsextra
          }else{
            mults.multAllTracksTPCOnly++; // multsextra
          }
        }
        // PV contributor checked explicitly
        if(track.isPVContributor()){
          if(std::abs(track.eta())<1.0){
            mults.multNContribsEta1++; // pvmults
            if(std::abs(track.eta())<0.8){
              mults.multNContribs++; // pvmults
              if(std::abs(track.eta())<0.5){
                mults.multNContribsEtaHalf++; // pvmults
              }
            }
          }
          if(track.hasITS()){ 
            mults.multHasITS++; // multsextra
            if (track.hasTPC())
              mults.multITSTPC++; // multsextra
            if (!track.hasTPC() && !track.hasTOF() && !track.hasTRD()){
              mults.multITSOnly++;  // multsextra
            }
          }
          if(track.hasTPC()){ 
            mults.multHasTPC++; // multsextra
            if (!track.hasITS() && !track.hasTOF() && !track.hasTRD()){
              mults.multTPCOnly++;  // multsextra
            }
          }
          if(track.hasTOF()){ 
            mults.multHasTOF++; // multsextra
          }
          if(track.hasTRD()){ 
            mults.multHasTRD++; // multsextra
          }
        }

        // global counters: do them only in case information is provided in tracks table
        if constexpr (requires { tracks.isQualityTrack(); }) {
          if(track.pt()<internalOpts.maxPtGlobalTrack.value && track.pt()>internalOpts.minPtGlobalTrack.value && std::fabs(track.eta())<1.0f && track.isPVContributor() && tracks.isQualityTrack()){ 
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
    }

    // fill track counters at this stage if requested
    if(internalOpts.mEnabledTables[kTPCMults]){
      cursors.tableTpc(mults.multTPC);
    } 
    if(internalOpts.mEnabledTables[kPVMults]){ 
      cursors.tablePv(mults.multNContribs, mults.multNContribsEta1, mults.multNContribsEtaHalf);
    }
    if(internalOpts.mEnabledTables[kMultsExtra]){
      cursors.tableExtra(collision.numContrib(), collision.chi2(), collision.collisionTimeRes(),
                         bc.runNumber(), collision.posZ(), collision.sel8(),
                         mults.multHasITS, mults.multHasTPC, mults.multHasTOF, mults.multHasTRD, 
                         mults.multITSOnly, mults.multTPCOnly, mults.multITSTPC,
                         mults.multAllTracksTPCOnly, mults.multAllTracksITSTPC, 
                         collision.trackOccupancyInTimeRange(),
                         collision.ft0cOccupancyInTimeRange(),
                         collision.flags());
    }
    if(internalOpts.mEnabledTables[kPVMultZeqs]){
      if(std::fabs(collision.posZ()) && lCalibLoaded){
        mults.multNContribsZeq = hVtxZNTracks->Interpolate(0.0) * mults.multNContribs / hVtxZNTracks->Interpolate(collision.posZ());
      }else{
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
    cursors.mftMults(nAllTracks, nTracks);
    mults[collision.globalIndex()].multMFTAllTracks = nAllTracks; 
    mults[collision.globalIndex()].multMFTTracks = nTracks; 
  }
}; // end BuilderModule

} // namespace multiplicity
} // namespace tools
} // namespace o2

#endif // COMMON_TOOLS_MULTMODULE_H_