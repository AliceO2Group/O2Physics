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
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//   Configurable workflow for running several DQ or other PWG analyses

#include <iostream>
#include <vector>
#include <algorithm>
#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TString.h>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "Common/Core/TableHelper.h"

using std::cout;
using std::endl;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Some definitions
namespace o2::aod
{
namespace dqanalysisflags
{
DECLARE_SOA_COLUMN(MixingHash, mixingHash, int);                                     //! Hash used in event mixing
DECLARE_SOA_BITMAP_COLUMN(IsEventSelected, isEventSelected, 32);                     //! Event decision
DECLARE_SOA_BITMAP_COLUMN(IsBarrelSelected, isBarrelSelected, 32);                   //! Barrel track decisions
DECLARE_SOA_COLUMN(BarrelAmbiguityInBunch, barrelAmbiguityInBunch, int8_t);          //! Barrel track in-bunch ambiguity
DECLARE_SOA_COLUMN(BarrelAmbiguityOutOfBunch, barrelAmbiguityOutOfBunch, int8_t);    //! Barrel track out of bunch ambiguity
DECLARE_SOA_BITMAP_COLUMN(IsMuonSelected, isMuonSelected, 32);                       //! Muon track decisions (joinable to ReducedMuonsAssoc)
DECLARE_SOA_BITMAP_COLUMN(IsBarrelSelectedPrefilter, isBarrelSelectedPrefilter, 32); //! Barrel prefilter decisions (joinable to ReducedTracksAssoc)
// Bcandidate columns for ML analysis of B->Jpsi+K
DECLARE_SOA_COLUMN(massBcandidate, MBcandidate, float);
DECLARE_SOA_COLUMN(pTBcandidate, PtBcandidate, float);
DECLARE_SOA_COLUMN(LxyBcandidate, lxyBcandidate, float);
DECLARE_SOA_COLUMN(LxyzBcandidate, lxyzBcandidate, float);
DECLARE_SOA_COLUMN(LzBcandidate, lzBcandidate, float);
DECLARE_SOA_COLUMN(TauxyBcandidate, tauxyBcandidate, float);
DECLARE_SOA_COLUMN(TauzBcandidate, tauzBcandidate, float);
DECLARE_SOA_COLUMN(CosPBcandidate, cosPBcandidate, float);
DECLARE_SOA_COLUMN(Chi2Bcandidate, chi2Bcandidate, float);
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);                                                            //!  joinable to ReducedEvents
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASH", dqanalysisflags::MixingHash);                                                             //!  joinable to ReducedEvents
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTS", dqanalysisflags::IsBarrelSelected);                                                    //!  joinable to ReducedTracksAssoc
DECLARE_SOA_TABLE(BarrelAmbiguities, "AOD", "DQBARRELAMB", dqanalysisflags::BarrelAmbiguityInBunch, dqanalysisflags::BarrelAmbiguityOutOfBunch); //!  joinable to ReducedBarrelTracks
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);                                                       //!  joinable to ReducedMuonsAssoc
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTER", dqanalysisflags::IsBarrelSelectedPrefilter);                                                  //!  joinable to ReducedTracksAssoc
DECLARE_SOA_TABLE(BmesonCandidates, "AOD", "DQBMESONS", dqanalysisflags::massBcandidate, dqanalysisflags::pTBcandidate, dqanalysisflags::LxyBcandidate, dqanalysisflags::LxyzBcandidate, dqanalysisflags::LzBcandidate, dqanalysisflags::TauxyBcandidate, dqanalysisflags::TauzBcandidate, dqanalysisflags::CosPBcandidate, dqanalysisflags::Chi2Bcandidate);
} // namespace o2::aod

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsHashSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts>;
using MyEventsVtxCovSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvector>;
using MyEventsQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector>;
using MyEventsHashSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvector>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksWithAmbiguities = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelAmbiguities>;
using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksWithCovWithAmbiguities = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelAmbiguities>;
using MyDielectronCandidates = soa::Join<aod::Dielectrons, aod::DielectronsExtra>;
using MyDimuonCandidates = soa::Join<aod::Dimuons, aod::DimuonsExtra>;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyMuonTracksSelectedWithColl = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsInfo, aod::MuonTrackCuts>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;
// constexpr static uint32_t gkEventFillMapWithQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector;
// constexpr static uint32_t gkEventFillMapWithCovQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
// constexpr static uint32_t gkTrackFillMapWithColl = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID | VarManager::ObjTypes::ReducedTrackCollInfo;

// constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCov;
// constexpr static uint32_t gkMuonFillMapWithColl = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCollInfo;

constexpr static int pairTypeEE = VarManager::kDecayToEE;
constexpr static int pairTypeMuMu = VarManager::kDecayToMuMu;
// constexpr static int pairTypeEMu = VarManager::kElectronMuon;

// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses, const char* histGroups); // defines histograms for all tasks

template <typename TMap>
void PrintBitMap(TMap map, int nbits)
{
  for (int i = 0; i < nbits; i++) {
    cout << ((map & (TMap(1) << i)) > 0 ? "1" : "0");
  }
}

// Analysis task that produces event decisions and the Hash table used in event mixing
struct AnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};
  // TODO: Provide the mixing variables and binning directly via configurables (e.g. vectors of float)
  Configurable<string> fConfigMixingVariables{"cfgMixingVars", "", "Mixing configs separated by a comma, default no mixing"};
  Configurable<string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};

  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;
  AnalysisCompositeCut* fEventCut;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

  std::map<int64_t, bool> fSelMap;                     // key: reduced event global index, value: event selection decision
  std::map<uint64_t, std::vector<int64_t>> fBCCollMap; // key: global BC, value: vector of reduced event global indices
  std::map<string, string> fMetadataRCT, fHeader;
  int fCurrentRun;

  void init(o2::framework::InitContext&)
  {
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    VarManager::SetDefaultVarNames();
    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;SameBunchCorrelations", fConfigAddEventHistogram.value.data()); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                                                             // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    TString mixVarsString = fConfigMixingVariables.value;
    std::unique_ptr<TObjArray> objArray(mixVarsString.Tokenize(","));
    if (objArray->GetEntries() > 0) {
      fMixHandler = new MixingHandler("mixingHandler", "mixing handler");
      fMixHandler->Init();
      for (int iVar = 0; iVar < objArray->GetEntries(); ++iVar) {
        dqmixing::SetUpMixing(fMixHandler, objArray->At(iVar)->GetName());
      }
    }

    fCurrentRun = -1;
    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
    fCCDBApi.init(fConfigCcdbUrl.value);
  }

  template <uint32_t TEventFillMap, typename TEvents>
  void runEventSelection(TEvents const& events)
  {
    if (events.size() > 0 && events.begin().runNumber() != fCurrentRun) {
      fHeader = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", events.begin().runNumber()), fMetadataRCT, -1);
      uint64_t sor = std::atol(fHeader["SOR"].c_str());
      uint64_t eor = std::atol(fHeader["EOR"].c_str());
      cout << "=========================== SOR / EOR is " << sor << " / " << eor << endl;
    }

    fSelMap.clear();
    fBCCollMap.clear();

    for (auto& event : events) {
      // Reset the fValues array and fill event observables
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEvent<TEventFillMap>(event);

      bool decision = false;
      // if QA is requested fill histograms before event selections
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues); // automatically fill all the histograms in the class Event
      }
      if (fEventCut->IsSelected(VarManager::fgValues)) {
        if (fConfigQA) {
          fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
        }
        decision = true;
      }
      fSelMap[event.globalIndex()] = decision;
      if (fBCCollMap.find(event.globalBC()) == fBCCollMap.end()) {
        std::vector<int64_t> evIndices = {event.globalIndex()};
        fBCCollMap[event.globalBC()] = evIndices;
      } else {
        auto& evIndices = fBCCollMap[event.globalBC()];
        evIndices.push_back(event.globalIndex());
      }

      if (fMixHandler != nullptr) {
        int hh = fMixHandler->FindEventCategory(VarManager::fgValues);
        hash(hh);
      }
    }
  }

  template <uint32_t TEventFillMap, typename TEvents>
  void publishSelections(TEvents const& events)
  {

    std::map<int64_t, bool> collisionSplittingMap; // key: event global index, value: whether pileup event is a possible splitting

    // Reset the fValues array and fill event observables
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    // loop over the BC map, find BCs with more than one collision and compute 2-event correlation quantities
    for (auto& [bc, evIndices] : fBCCollMap) {
      if (evIndices.size() < 2) {
        continue;
      }
      for (auto ev1Idx = evIndices.begin(); ev1Idx != evIndices.end(); ++ev1Idx) {
        if (!fSelMap[*ev1Idx]) {
          continue;
        }
        auto ev1 = events.rawIteratorAt(*ev1Idx);
        for (auto ev2Idx = std::next(ev1Idx); ev2Idx != evIndices.end(); ++ev2Idx) {
          if (!fSelMap[*ev2Idx]) {
            continue;
          }
          auto ev2 = events.rawIteratorAt(*ev2Idx);
          VarManager::FillTwoEvents(ev1, ev2);
          if (TMath::Abs(VarManager::fgValues[VarManager::kTwoEvDeltaZ]) < 1.0) { // this is a possible collision split
            collisionSplittingMap[*ev1Idx] = true;
            collisionSplittingMap[*ev2Idx] = true;
          }
          fHistMan->FillHistClass("SameBunchCorrelations", VarManager::fgValues);
        }
      }
    }

    // publish the table
    uint32_t evSel = 0;
    for (auto& event : events) {
      evSel = 0;
      if (fSelMap[event.globalIndex()]) { // event passed the user cuts
        evSel |= (uint32_t(1) << 0);
      }
      std::vector<int64_t> sameBunchEvents = fBCCollMap[event.globalBC()];
      if (sameBunchEvents.size() > 1) { // event with in-bunch pileup
        evSel |= (uint32_t(1) << 1);
      }
      if (collisionSplittingMap.find(event.globalIndex()) != collisionSplittingMap.end()) { // event with possible fake in-bunch pileup (collision splitting)
        evSel |= (uint32_t(1) << 2);
      }
      eventSel(evSel);
    }
  }

  void processSkimmed(MyEvents const& events)
  {
    runEventSelection<gkEventFillMap>(events);
    publishSelections<gkEventFillMap>(events);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDummy, "Dummy function", false);
};

// Produces a table with barrel track decisions (joinable to the ReducedTracksAssociations)
// Here one should add all the track cuts needed through the workflow (e.g. cuts for same-even pairing, electron prefiltering, track for dilepton-track correlations)
struct AnalysisTrackSelection {
  Produces<aod::BarrelTrackCuts> trackSel;
  Produces<aod::BarrelAmbiguities> trackAmbiguities;
  OutputObj<THashList> fOutputList{"output"};

  Configurable<string> fConfigCuts{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<string> fConfigCcdbPathTPC{"ccdb-path-tpc", "Users/z/zhxiong/TPCPID/PostCalib", "base path to the ccdb object"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<bool> fConfigComputeTPCpostCalib{"cfgTPCpostCalib", false, "If true, compute TPC post-calibrated n-sigmas"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> fConfigDummyRunlist{"cfgDummyRunlist", false, "If true, use dummy runlist"};
  Configurable<int> fConfigInitRunNumber{"cfgInitRunNumber", 543215, "Initial run number used in run by run checks"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fTrackCuts;

  int fCurrentRun; // current run (needed to detect run changes for loading CCDB parameters)

  std::map<int64_t, std::vector<int64_t>> fNAssocsInBunch;    // key: track global index, value: vector of global index for events associated in-bunch (events that have in-bunch pileup or splitting)
  std::map<int64_t, std::vector<int64_t>> fNAssocsOutOfBunch; // key: track global index, value: vector of global index for events associated out-of-bunch (events that have no in-bunch pileup)

  void init(o2::framework::InitContext&)
  {
    fCurrentRun = 0;

    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      // set one histogram directory for each defined track cut
      TString histDirNames = "TrackBarrel_BeforeCuts;";
      for (auto& cut : fTrackCuts) {
        histDirNames += Form("TrackBarrel_%s;", cut.GetName());
      }
      histDirNames += "TrackBarrel_AmbiguityInBunch;TrackBarrel_AmbiguityOutOfBunch;";

      DefineHistograms(fHistMan, histDirNames.Data(), fConfigAddTrackHistogram.value.data()); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                        // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
    if (fConfigDummyRunlist) {
      VarManager::SetDummyRunlist(fConfigInitRunNumber);
    }
    if (fConfigComputeTPCpostCalib) {
      fCCDB->setURL(fConfigCcdbUrl.value);
      fCCDB->setCaching(true);
      fCCDB->setLocalObjectValidityChecking();
      fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
    }
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTracks>
  void runTrackSelection(ReducedTracksAssoc const& assocs, TEvents const& events, TTracks const& tracks)
  {
    fNAssocsInBunch.clear();
    fNAssocsOutOfBunch.clear();

    if (fConfigComputeTPCpostCalib && events.size() > 0 && fCurrentRun != events.begin().runNumber()) {
      auto calibList = fCCDB->getForTimeStamp<TList>(fConfigCcdbPathTPC.value, events.begin().timestamp());
      VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
      VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
      VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
      VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
      VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
      VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));

      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, events.begin().timestamp());
      if (grpmag != nullptr) {
        VarManager::SetMagneticField(grpmag->getNominalL3Field());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", events.begin().timestamp());
      }

      fCurrentRun = events.begin().runNumber();
    }

    trackSel.reserve(assocs.size());
    trackAmbiguities.reserve(tracks.size());
    uint32_t filterMap = 0;
    int iCut = 0;

    for (auto& assoc : assocs) {
      auto event = assoc.template reducedevent_as<TEvents>();
      if (!event.isEventSelected_bit(0)) {
        trackSel(0);
        continue;
      }
      VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
      // fill event information which might be needed in histograms/cuts that combine track and event properties
      VarManager::FillEvent<TEventFillMap>(event);

      auto track = assoc.template reducedtrack_as<TTracks>();
      filterMap = 0;
      VarManager::FillTrack<TTrackFillMap>(track);
      // compute quantities which depend on the associated collision, such as DCA
      VarManager::FillTrackCollision<TTrackFillMap>(track, event);
      if (fConfigQA) {
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }
      iCut = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, iCut++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << iCut);
          if (fConfigQA) {
            fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
          }
        }
      } // end loop over cuts

      trackSel(filterMap);

      // count the number of associations per track
      if (filterMap > 0) {
        if (event.isEventSelected_bit(1)) {
          if (fNAssocsInBunch.find(track.globalIndex()) == fNAssocsInBunch.end()) {
            std::vector<int64_t> evVector = {event.globalIndex()};
            fNAssocsInBunch[track.globalIndex()] = evVector;
          } else {
            auto& evVector = fNAssocsInBunch[track.globalIndex()];
            evVector.push_back(event.globalIndex());
          }
        } else {
          if (fNAssocsOutOfBunch.find(track.globalIndex()) == fNAssocsOutOfBunch.end()) {
            std::vector<int64_t> evVector = {event.globalIndex()};
            fNAssocsOutOfBunch[track.globalIndex()] = evVector;
          } else {
            auto& evVector = fNAssocsOutOfBunch[track.globalIndex()];
            evVector.push_back(event.globalIndex());
          }
        }
      }
    } // end loop over associations

    // QA the collision-track associations
    for (auto& [trackIdx, evIndices] : fNAssocsInBunch) {
      if (evIndices.size() == 1) {
        continue;
      }
      auto track = tracks.rawIteratorAt(trackIdx);
      VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
      VarManager::FillTrack<TTrackFillMap>(track);
      VarManager::fgValues[VarManager::kBarrelNAssocsInBunch] = static_cast<float>(evIndices.size());
      fHistMan->FillHistClass("TrackBarrel_AmbiguityInBunch", VarManager::fgValues);
    } // end loop over in-bunch ambiguous tracks

    for (auto& [trackIdx, evIndices] : fNAssocsOutOfBunch) {
      if (evIndices.size() == 1) {
        continue;
      }
      auto track = tracks.rawIteratorAt(trackIdx);
      VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
      VarManager::FillTrack<TTrackFillMap>(track);
      VarManager::fgValues[VarManager::kBarrelNAssocsOutOfBunch] = static_cast<float>(evIndices.size());
      fHistMan->FillHistClass("TrackBarrel_AmbiguityOutOfBunch", VarManager::fgValues);
    } // end loop over out-of-bunch ambiguous tracks

    // publish the ambiguity table
    for (auto& track : tracks) {
      int8_t nInBunch = 0;
      if (fNAssocsInBunch.find(track.globalIndex()) != fNAssocsInBunch.end()) {
        nInBunch = fNAssocsInBunch[track.globalIndex()].size();
      }
      int8_t nOutOfBunch = 0;
      if (fNAssocsOutOfBunch.find(track.globalIndex()) != fNAssocsOutOfBunch.end()) {
        nOutOfBunch = fNAssocsOutOfBunch[track.globalIndex()].size();
      }
      trackAmbiguities(nInBunch, nOutOfBunch);
    }

  } // end runTrackSelection()

  void processSkimmed(ReducedTracksAssoc const& assocs, MyEventsSelected const& events, MyBarrelTracks const& tracks)
  {
    runTrackSelection<gkEventFillMap, gkTrackFillMap>(assocs, events, tracks);
  }
  void processSkimmedWithCov(ReducedTracksAssoc const& assocs, MyEventsVtxCovSelected const& events, MyBarrelTracksWithCov const& tracks)
  {
    runTrackSelection<gkEventFillMapWithCov, gkTrackFillMapWithCov>(assocs, events, tracks);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmed, "Run barrel track selection on DQ skimmed track associations", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmedWithCov, "Run barrel track selection on DQ skimmed tracks w/ cov matrix associations", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummy, "Dummy function", false);
};

// Produces a table with muon decisions (joinable to the ReducedMuonsAssociations)
// Here one should add all the track cuts needed through the workflow (e.g. cuts for same-event pairing, track for dilepton-track correlations)
struct AnalysisMuonSelection {
  Produces<aod::MuonTrackCuts> muonSel;
  OutputObj<THashList> fOutputList{"output"};

  Configurable<string> fConfigCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddMuonHistogram{"cfgAddMuonHistogram", "", "Comma separated list of histograms"};
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> fConfigGeoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fMuonCuts;

  int fCurrentRun; // current run kept to detect run changes and trigger loading params from CCDB

  void init(o2::framework::InitContext&)
  {
    fCurrentRun = 0;
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fMuonCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      // set one histogram directory for each defined track cut
      TString histDirNames = "TrackMuon_BeforeCuts;";
      for (auto& cut : fMuonCuts) {
        histDirNames += Form("TrackMuon_%s;", cut.GetName());
      }

      DefineHistograms(fHistMan, histDirNames.Data(), fConfigAddMuonHistogram.value.data()); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                       // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      fCCDB->get<TGeoManager>(fConfigGeoPath);
    }
  }

  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvents, typename TMuons>
  void runMuonSelection(ReducedMuonsAssoc const& assocs, TEvents const& events, TMuons const& /*muons*/)
  {
    if (events.size() > 0 && fCurrentRun != events.begin().runNumber()) {
      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, events.begin().timestamp());
      if (grpmag != nullptr) {
        o2::base::Propagator::initFieldFromGRP(grpmag);
        VarManager::SetMagneticField(grpmag->getNominalL3Field());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", events.begin().timestamp());
      }
      fCurrentRun = events.begin().runNumber();
    }

    muonSel.reserve(assocs.size());
    uint32_t filterMap = 0;
    int iCut = 0;

    for (auto& assoc : assocs) {
      auto event = assoc.template reducedevent_as<TEvents>();
      if (!event.isEventSelected_bit(0)) {
        muonSel(0);
        continue;
      }
      VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
      // fill event information which might be needed in histograms/cuts that combine track and event properties
      VarManager::FillEvent<TEventFillMap>(event);

      auto track = assoc.template reducedmuon_as<TMuons>();
      filterMap = 0;
      VarManager::FillTrack<TMuonFillMap>(track);
      // compute quantities which depend on the associated collision
      VarManager::FillPropagateMuon<TMuonFillMap>(track, event);
      if (fConfigQA) {
        fHistMan->FillHistClass("TrackMuon_BeforeCuts", VarManager::fgValues);
      }
      iCut = 0;
      for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, iCut++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << iCut);
          if (fConfigQA) {
            fHistMan->FillHistClass(Form("TrackMuon_%s", (*cut).GetName()), VarManager::fgValues);
          }
        }
      } // end loop over cuts
      muonSel(filterMap);
    } // end loop over assocs
  }

  void processSkimmed(ReducedMuonsAssoc const& assocs, MyEventsVtxCovSelected const& events, MyMuonTracksWithCov const& muons)
  {
    runMuonSelection<gkEventFillMapWithCov, gkMuonFillMapWithCov>(assocs, events, muons);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisMuonSelection, processSkimmed, "Run muon selection on DQ skimmed muons", false);
  PROCESS_SWITCH(AnalysisMuonSelection, processDummy, "Dummy function", false);
};

// Run the prefilter selection (e.g. electron prefiltering for photon conversions)
// This takes uses a sample of tracks selected with loose cuts (fConfigPrefilterTrackCut) and combines them
//  with the sample of tracks to be used in downstream analysis (fConfigTrackCuts). If a pair is found to pass
//  the pair prefilter cut (cfgPrefilterPairCut), the analysis track is tagged to be removed from analysis.
struct AnalysisPrefilterSelection {
  Produces<aod::Prefilter> prefilter; // joinable with ReducedTracksAssoc

  // Configurables
  Configurable<std::string> fConfigPrefilterTrackCut{"cfgPrefilterTrackCut", "", "Prefilter track cut"};
  Configurable<std::string> fConfigPrefilterPairCut{"cfgPrefilterPairCut", "", "Prefilter pair cut"};
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Track cuts for which to run the prefilter"};

  std::map<uint32_t, uint32_t> fPrefilterMap;
  AnalysisCompositeCut* fPairCut;
  uint32_t fPrefilterMask;
  int fPrefilterCutBit;

  Preslice<aod::ReducedTracksAssoc> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;

  void init(o2::framework::InitContext& initContext)
  {
    // get the list of track cuts to be prefiltered
    TString trackCutsStr = fConfigTrackCuts.value;
    TObjArray* objArrayTrackCuts = nullptr;
    if (!trackCutsStr.IsNull()) {
      objArrayTrackCuts = trackCutsStr.Tokenize(",");
    }
    if (objArrayTrackCuts->GetEntries() == 0) {
      LOG(fatal) << " No track cuts to prefilter!";
    }

    // get the list of cuts that were computed in the barrel track-selection task and create a bit mask
    //  to mark just the ones we want to apply a prefilter on
    fPrefilterMask = 0;
    fPrefilterCutBit = -1;
    string trackCuts;
    getTaskOptionValue<string>(initContext, "analysis-track-selection", "cfgTrackCuts", trackCuts, false);
    TString allTrackCutsStr = trackCuts;
    TString prefilterTrackCutStr = fConfigPrefilterTrackCut.value;
    if (!trackCutsStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(allTrackCutsStr.Tokenize(","));
      if (objArray->FindObject(prefilterTrackCutStr.Data()) == nullptr) {
        LOG(fatal) << " Prefilter track cut not among the cuts calculated by the track-selection task! ";
      }
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        TString tempStr = objArray->At(icut)->GetName();
        if (objArrayTrackCuts->FindObject(tempStr.Data()) != nullptr) {
          fPrefilterMask |= (uint32_t(1) << icut);
        }
        if (tempStr.CompareTo(fConfigPrefilterTrackCut.value) == 0) {
          fPrefilterCutBit = icut;
        }
      }
    }
    if (fPrefilterMask == 0) {
      LOG(fatal) << "No specified track cuts for prefiltering";
    }
    if (fPrefilterCutBit < 0) {
      LOG(fatal) << "No or incorrectly specified loose track prefilter cut";
    }

    // setup the prefilter pair cut
    fPairCut = new AnalysisCompositeCut(true);
    TString pairCutStr = fConfigPrefilterPairCut.value;
    if (!pairCutStr.IsNull()) {
      fPairCut = dqcuts::GetCompositeCut(pairCutStr.Data());
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    VarManager::SetDefaultVarNames();

    VarManager::SetupTwoProngDCAFitter(5.0f, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, true); // TODO: get these parameters from Configurables
    VarManager::SetupTwoProngFwdDCAFitter(5.0f, true, 200.0f, 1.0e-3f, 0.9f, true);
  }

  template <uint32_t TTrackFillMap, typename TTracks>
  void runPrefilter(soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& assocs, TTracks const& /*tracks*/)
  {

    for (auto& [assoc1, assoc2] : o2::soa::combinations(assocs, assocs)) {
      auto track1 = assoc1.template reducedtrack_as<TTracks>();
      auto track2 = assoc2.template reducedtrack_as<TTracks>();

      // NOTE: here we restrict to just pairs of opposite sign (conversions), but in principle this can be made
      // a configurable and check also same-sign pairs (track splitting)
      if (track1.sign() * track2.sign() > 0) {
        continue;
      }

      // here we check the cuts fulfilled by both tracks, for both the tight and loose selections
      uint32_t track1Candidate = (assoc1.isBarrelSelected_raw() & fPrefilterMask);
      uint32_t track2Candidate = (assoc2.isBarrelSelected_raw() & fPrefilterMask);
      bool track1Loose = assoc1.isBarrelSelected_bit(fPrefilterCutBit);
      bool track2Loose = assoc2.isBarrelSelected_bit(fPrefilterCutBit);

      if (!((track1Candidate > 0 && track2Loose) || (track2Candidate > 0 && track1Loose))) {
        continue;
      }

      // compute pair quantities
      VarManager::FillPair<VarManager::kDecayToEE, TTrackFillMap>(track1, track2);
      // if the pair fullfils the criteria, add an entry into the prefilter map for the two tracks
      if (fPairCut->IsSelected(VarManager::fgValues)) {
        if (fPrefilterMap.find(track1.globalIndex()) == fPrefilterMap.end() && track1Candidate > 0) {
          fPrefilterMap[track1.globalIndex()] = track1Candidate;
        }
        if (fPrefilterMap.find(track2.globalIndex()) == fPrefilterMap.end() && track2Candidate > 0) {
          fPrefilterMap[track2.globalIndex()] = track2Candidate;
        }
      }
    } // end loop over combinations
  }

  void processBarrelSkimmed(MyEvents const& events, soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& assocs, MyBarrelTracks const& tracks)
  {

    fPrefilterMap.clear();

    for (auto& event : events) {
      auto groupedAssocs = assocs.sliceBy(trackAssocsPerCollision, event.globalIndex());
      if (groupedAssocs.size() > 1) {
        runPrefilter<gkTrackFillMap>(groupedAssocs, tracks);
      }
    }
    uint32_t mymap = -1;
    for (auto& assoc : assocs) {
      auto track = assoc.template reducedtrack_as<MyBarrelTracks>();
      mymap = -1;
      if (fPrefilterMap.find(track.globalIndex()) != fPrefilterMap.end()) {
        // NOTE: publish the bitwise negated bits (~), so there will be zeroes for cuts that failed the prefiltering and 1 everywhere else
        mymap = ~fPrefilterMap[track.globalIndex()];
      }
      prefilter(mymap);
    }
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisPrefilterSelection, processBarrelSkimmed, "Run Prefilter selection on reduced tracks", false);
  PROCESS_SWITCH(AnalysisPrefilterSelection, processDummy, "Do nothing", false);
};

// Run the same-event pairing
// This task assumes that both legs of the resonance fulfill the same cuts (symmetric decay channel)
//  Runs combinatorics for barrel-barrel, muon-muon and barrel-muon combinations
// TODO: implement properly the barrel-muon combinations
// The task implements also process functions for running event mixing
struct AnalysisSameEventPairing {

  Produces<aod::Dielectrons> dielectronList;
  Produces<aod::Dimuons> dimuonList;
  Produces<aod::DielectronsExtra> dielectronsExtraList;
  Produces<aod::DimuonsExtra> dimuonsExtraList;
  Produces<aod::DimuonsAll> dimuonAllList;
  Produces<aod::DileptonFlow> dileptonFlowList;
  Produces<aod::DileptonsInfo> dileptonInfoList;

  o2::base::MatLayerCylSet* fLUT = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  OutputObj<THashList> fOutputList{"output"};

  Configurable<string> fConfigTrackCuts{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  Configurable<string> fConfigMuonCuts{"cfgMuonCuts", "", "Comma separated list of muon cuts"};
  Configurable<string> fConfigPairCuts{"cfgPairCuts", "", "Comma separated list of pair cuts"};

  Configurable<int> fConfigMixingDepth{"cfgMixingDepth", 100, "Number of Events stored for event mixing"};
  Configurable<std::string> fConfigAddEventMixingHistogram{"cfgAddEventMixingHistogram", "", "Comma separated list of histograms"};

  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<string> fConfigCcdbPath{"ccdb-path", "Users/lm", "base path to the ccdb object"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> fConfigGRPMagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> fConfigUseRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};

  Configurable<std::string> fConfigAddSEPHistogram{"cfgAddSEPHistogram", "", "Comma separated list of histograms"};
  Configurable<bool> fConfigFlatTables{"cfgFlatTables", false, "Produce a single flat tables with all relevant information of the pairs and single tracks"};
  Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
  Configurable<bool> fConfigUseAbsDCA{"cfgUseAbsDCA", false, "Use absolute DCA minimization instead of chi^2 minimization in secondary vertexing"};
  Configurable<bool> fConfigPropToPCA{"cfgPropToPCA", false, "Propagate tracks to secondary vertex"};
  Configurable<bool> fConfigCorrFullGeo{"cfgCorrFullGeo", false, "Use full geometry to correct for MCS effects in track propagation"};
  Configurable<bool> fConfigNoCorr{"cfgNoCorrFwdProp", false, "Do not correct for MCS effects in track propagation"};
  Configurable<std::string> fConfigLutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> fConfigGeoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> fConfigCollisionSystem{"syst", "pp", "Collision system, pp or PbPb"};
  Configurable<float> fConfigCenterMassEnergy{"energy", 13600, "Center of mass energy in GeV"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected == uint32_t(1);

  HistogramManager* fHistMan;

  // keep histogram class names in maps, so we don't have to buld their names in the pair loops
  std::map<int, std::vector<TString>> fTrackHistNames;
  std::map<int, std::vector<TString>> fMuonHistNames;
  std::map<int, std::vector<TString>> fTrackMuonHistNames;
  std::vector<AnalysisCompositeCut> fPairCuts;

  uint32_t fTrackFilterMask; // mask for the track cuts required in this task to be applied on the barrel cuts produced upstream
  uint32_t fMuonFilterMask;  // mask for the muon cuts required in this task to be applied on the muon cuts produced upstream
  int fNCutsBarrel;
  int fNCutsMuon;
  int fNPairCuts;

  bool fEnableBarrelMixingHistos;
  bool fEnableBarrelHistos;
  bool fEnableMuonHistos;
  bool fEnableMuonMixingHistos;

  NoBinningPolicy<aod::dqanalysisflags::MixingHash> hashBin;

  Preslice<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter>> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;
  Preslice<soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts>> muonAssocsPerCollision = aod::reducedtrack_association::reducedeventId;

  void init(o2::framework::InitContext& context)
  {
    fEnableBarrelHistos = context.mOptions.get<bool>("processAllSkimmed") || context.mOptions.get<bool>("processBarrelOnlySkimmed");
    fEnableBarrelMixingHistos = context.mOptions.get<bool>("processMixingAllSkimmed");
    fEnableMuonHistos = context.mOptions.get<bool>("processAllSkimmed") || context.mOptions.get<bool>("processMuonOnlySkimmed");
    fEnableMuonMixingHistos = context.mOptions.get<bool>("processMixingAllSkimmed");
    bool isDummy = context.mOptions.get<bool>("processDummy");
    if (isDummy) {
      if (fEnableBarrelHistos || fEnableBarrelMixingHistos || fEnableMuonHistos || fEnableMuonMixingHistos) {
        LOG(warning) << "The dummy process function is enabled while you have enabled normal process function. Check your configuration file!" << endl;
      } else {
        LOG(info) << "Dummy function enabled. Skipping the rest of init()" << endl;
        return;
      }
    }

    // Keep track of all the histogram class names to avoid composing strings in the pairing loop
    TString histNames = "";
    std::vector<TString> names;

    TString cutNamesStr = fConfigPairCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // get the list of cuts for tracks/muons, check that they were played by the barrel/muon selection tasks
    //   and make a mask for active cuts (barrel and muon selection tasks may run more cuts, needed for other analyses)
    TString trackCutsStr = fConfigTrackCuts.value;
    TObjArray* objArrayTrackCuts = nullptr;
    if (!trackCutsStr.IsNull()) {
      objArrayTrackCuts = trackCutsStr.Tokenize(",");
    }
    TString muonCutsStr = fConfigMuonCuts.value;
    TObjArray* objArrayMuonCuts = nullptr;
    if (!muonCutsStr.IsNull()) {
      objArrayMuonCuts = muonCutsStr.Tokenize(",");
    }

    // get the barrel track selection cuts
    string tempCuts;
    getTaskOptionValue<string>(context, "analysis-track-selection", "cfgTrackCuts", tempCuts, false);
    TString tempCutsStr = tempCuts;
    if (!trackCutsStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(tempCutsStr.Tokenize(","));
      fNCutsBarrel = objArray->GetEntries();
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        TString tempStr = objArray->At(icut)->GetName();
        if (objArrayTrackCuts->FindObject(tempStr.Data()) != nullptr) {
          fTrackFilterMask |= (uint32_t(1) << icut);

          if (fEnableBarrelHistos) {
            names = {
              Form("PairsBarrelSEPM_%s", objArray->At(icut)->GetName()),
              Form("PairsBarrelSEPP_%s", objArray->At(icut)->GetName()),
              Form("PairsBarrelSEMM_%s", objArray->At(icut)->GetName())};
            histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
            if (fEnableBarrelMixingHistos) {
              names.push_back(Form("PairsBarrelMEPM_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelMEPP_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelMEMM_%s", objArray->At(icut)->GetName()));
              histNames += Form("%s;%s;%s;", names[3].Data(), names[4].Data(), names[5].Data());
            }
            names.push_back(Form("PairsBarrelSEPM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
            names.push_back(Form("PairsBarrelSEPP_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
            names.push_back(Form("PairsBarrelSEMM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
            names.push_back(Form("PairsBarrelSEPM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
            names.push_back(Form("PairsBarrelSEPP_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
            names.push_back(Form("PairsBarrelSEMM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
            histNames += Form("%s;%s;%s;", names[(fEnableBarrelMixingHistos ? 6 : 3)].Data(), names[(fEnableBarrelMixingHistos ? 7 : 4)].Data(), names[(fEnableBarrelMixingHistos ? 8 : 5)].Data());
            histNames += Form("%s;%s;%s;", names[(fEnableBarrelMixingHistos ? 9 : 6)].Data(), names[(fEnableBarrelMixingHistos ? 10 : 7)].Data(), names[(fEnableBarrelMixingHistos ? 11 : 8)].Data());
            fTrackHistNames[icut] = names;

            TString cutNamesStr = fConfigPairCuts.value;
            if (!cutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
              fNPairCuts = objArrayPair->GetEntries();
              for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
                names = {
                  Form("PairsBarrelSEPM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsBarrelSEPP_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsBarrelSEMM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName())};
                histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
                fTrackHistNames[fNCutsBarrel + icut * fNPairCuts + iPairCut] = names;
              } // end loop (pair cuts)
            }   // end if (pair cuts)
          }     // end if enableBarrelHistos
        }
      }
    }
    // get the muon track selection cuts
    getTaskOptionValue<string>(context, "analysis-muon-selection", "cfgMuonCuts", tempCuts, false);
    tempCutsStr = tempCuts;
    if (!muonCutsStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(tempCutsStr.Tokenize(","));
      fNCutsMuon = objArray->GetEntries();
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        TString tempStr = objArray->At(icut)->GetName();
        if (objArrayMuonCuts->FindObject(tempStr.Data()) != nullptr) {
          fMuonFilterMask |= (uint32_t(1) << icut);

          if (fEnableMuonHistos) {
            // no pair cuts
            names = {
              Form("PairsMuonSEPM_%s", objArray->At(icut)->GetName()),
              Form("PairsMuonSEPP_%s", objArray->At(icut)->GetName()),
              Form("PairsMuonSEMM_%s", objArray->At(icut)->GetName())};
            histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
            if (fEnableMuonMixingHistos) {
              names.push_back(Form("PairsMuonMEPM_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonMEPP_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonMEMM_%s", objArray->At(icut)->GetName()));
              histNames += Form("%s;%s;%s;", names[3].Data(), names[4].Data(), names[5].Data());
            }
            fMuonHistNames[icut] = names;

            TString cutNamesStr = fConfigPairCuts.value;
            if (!cutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
              fNPairCuts = objArrayPair->GetEntries();
              for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
                names = {
                  Form("PairsMuonSEPM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsMuonSEPP_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsMuonSEMM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName())};
                histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
                fMuonHistNames[fNCutsMuon + icut * fNCutsMuon + iPairCut] = names;
              } // end loop (pair cuts)
            }   // end if (pair cuts)
          }
        }
      }
    }

    fCurrentRun = 0;

    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();

    if (fConfigNoCorr) {
      VarManager::SetupFwdDCAFitterNoCorr();
    } else if (fConfigCorrFullGeo || (fConfigUseKFVertexing && fConfigPropToPCA)) {
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        fCCDB->get<TGeoManager>(fConfigGeoPath);
      }
    } else {
      fLUT = o2::base::MatLayerCylSet::rectifyPtrFromFile(fCCDB->get<o2::base::MatLayerCylSet>(fConfigLutPath));
      VarManager::SetupMatLUTFwdDCAFitter(fLUT);
    }

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    /*if (context.mOptions.get<bool>("processElectronMuonSkimmed") || context.mOptions.get<bool>("processAllSkimmed")) {
      TString cutNamesBarrel = fConfigTrackCuts.value;
      TString cutNamesMuon = fConfigMuonCuts.value;
      if (!cutNamesBarrel.IsNull() && !cutNamesMuon.IsNull()) {
        std::unique_ptr<TObjArray> objArrayBarrel(cutNamesBarrel.Tokenize(","));
        std::unique_ptr<TObjArray> objArrayMuon(cutNamesMuon.Tokenize(","));
        if (objArrayBarrel->GetEntries() == objArrayMuon->GetEntries()) {   // one must specify equal number of barrel and muon cuts
          for (int icut = 0; icut < objArrayBarrel->GetEntries(); ++icut) { // loop over track cuts
            // no pair cuts
            names = {
              Form("PairsEleMuSEPM_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName()),
              Form("PairsEleMuSEPP_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName()),
              Form("PairsEleMuSEMM_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName())};
            histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
            fTrackMuonHistNames.push_back(names);

            TString cutNamesStr = fConfigPairCuts.value;
            if (!cutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
              for (int iPairCut = 0; iPairCut < objArrayPair->GetEntries(); ++iPairCut) { // loop over pair cuts
                std::vector<TString> names = {
                  Form("PairsEleMuSEPM_%s_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsEleMuSEPP_%s_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsEleMuSEMM_%s_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName())};
                histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
                fTrackMuonHistNames.push_back(names);
              } // end loop (pair cuts)
            }   // end if (pair cuts)
          }     // end loop (track cuts)
        }       // end if (equal number of cuts)
      }         // end if (track cuts)
    }*/

    VarManager::SetCollisionSystem((TString)fConfigCollisionSystem, fConfigCenterMassEnergy); // set collision system and center of mass energy

    DefineHistograms(fHistMan, histNames.Data(), fConfigAddSEPHistogram.value.data()); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());                                   // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  void initParamsFromCCDB(uint64_t timestamp, bool withTwoProngFitter = true)
  {

    if (fConfigUseRemoteField.value) {
      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fConfigGRPMagPath, timestamp);
      float magField = 0.0;
      if (grpmag != nullptr) {
        magField = grpmag->getNominalL3Field();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", timestamp);
      }
      if (withTwoProngFitter) {
        if (fConfigUseKFVertexing.value) {
          VarManager::SetupTwoProngKFParticle(magField);
        } else {
          VarManager::SetupTwoProngDCAFitter(magField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value); // TODO: get these parameters from Configurables
          VarManager::SetupTwoProngFwdDCAFitter(magField, true, 200.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value);
        }
      } else {
        VarManager::SetupTwoProngDCAFitter(magField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value); // needed because take in varmanager Bz from fgFitterTwoProngBarrel for PhiV calculations
      }
    } else {
      if (withTwoProngFitter) {
        if (fConfigUseKFVertexing.value) {
          VarManager::SetupTwoProngKFParticle(fConfigMagField.value);
        } else {
          VarManager::SetupTwoProngDCAFitter(fConfigMagField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value); // TODO: get these parameters from Configurables
          VarManager::SetupTwoProngFwdDCAFitter(fConfigMagField.value, true, 200.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value);
        }
      } else {
        VarManager::SetupTwoProngDCAFitter(fConfigMagField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value); // needed because take in varmanager Bz from fgFitterTwoProngBarrel for PhiV calculations
      }
    }
  }

  // Template function to run same event pairing (barrel-barrel, muon-muon, barrel-muon)
  template <bool TTwoProngFitter, int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTrackAssocs, typename TTracks>
  void runSameEventPairing(TEvents const& events, Preslice<TTrackAssocs>& preslice, TTrackAssocs const& assocs, TTracks const& /*tracks*/)
  {
    if (fCurrentRun != events.begin().runNumber()) {
      initParamsFromCCDB(events.begin().timestamp(), TTwoProngFitter);
      fCurrentRun = events.begin().runNumber();
    }

    TString cutNames = fConfigTrackCuts.value;
    std::map<int, std::vector<TString>> histNames = fTrackHistNames;
    int ncuts = fNCutsBarrel;
    if constexpr (TPairType == pairTypeMuMu) {
      cutNames = fConfigMuonCuts.value;
      histNames = fMuonHistNames;
      ncuts = fNCutsMuon;
    }
    int histIdxOffset = 0;
    if constexpr (TPairType == pairTypeEE) {
      if (fEnableBarrelMixingHistos) {
        histIdxOffset = 3;
      }
    }
    /*if constexpr (TPairType == pairTypeEMu) {
      cutNames = fConfigMuonCuts.value;
      histNames = fTrackMuonHistNames;
    }*/

    uint32_t twoTrackFilter = 0;
    int sign1 = 0;
    int sign2 = 0;
    dielectronList.reserve(1);
    dimuonList.reserve(1);
    dielectronsExtraList.reserve(1);
    dimuonsExtraList.reserve(1);
    dileptonInfoList.reserve(1);
    dileptonFlowList.reserve(1);
    if (fConfigFlatTables.value) {
      dimuonAllList.reserve(1);
    }
    constexpr bool eventHasQvector = ((TEventFillMap & VarManager::ObjTypes::ReducedEventQvector) > 0);
    constexpr bool eventHasQvectorCentr = ((TEventFillMap & VarManager::ObjTypes::CollisionQvect) > 0);
    constexpr bool trackHasCov = ((TTrackFillMap & VarManager::ObjTypes::TrackCov) > 0 || (TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelCov) > 0);

    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);

      auto groupedAssocs = assocs.sliceBy(preslice, event.globalIndex());
      if (groupedAssocs.size() == 0) {
        continue;
      }

      for (auto& [a1, a2] : o2::soa::combinations(groupedAssocs, groupedAssocs)) {

        if constexpr (TPairType == VarManager::kDecayToEE || TPairType == VarManager::kDecayToPiPi) {
          twoTrackFilter = a1.isBarrelSelected_raw() & a2.isBarrelSelected_raw() & a1.isBarrelSelectedPrefilter_raw() & a2.isBarrelSelectedPrefilter_raw() & fTrackFilterMask;

          if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }

          auto t1 = a1.template reducedtrack_as<TTracks>();
          auto t2 = a2.template reducedtrack_as<TTracks>();
          sign1 = t1.sign();
          sign2 = t2.sign();
          // store the ambiguity number of the two dilepton legs in the last 4 digits of the two-track filter
          if (t1.barrelAmbiguityInBunch() > 1) {
            twoTrackFilter |= (uint32_t(1) << 28);
          }
          if (t2.barrelAmbiguityInBunch() > 1) {
            twoTrackFilter |= (uint32_t(1) << 29);
          }
          if (t1.barrelAmbiguityOutOfBunch() > 1) {
            twoTrackFilter |= (uint32_t(1) << 30);
          }
          if (t2.barrelAmbiguityOutOfBunch() > 1) {
            twoTrackFilter |= (uint32_t(1) << 31);
          }

          VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
          if constexpr (TTwoProngFitter) {
            VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, fConfigPropToPCA);
          }
          if constexpr (eventHasQvector) {
            VarManager::FillPairVn<TPairType>(t1, t2);
          }

          dielectronList(event.globalIndex(), VarManager::fgValues[VarManager::kMass],
                         VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi],
                         t1.sign() + t2.sign(), twoTrackFilter, 0);

          if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrackCollInfo) > 0) {
            dileptonInfoList(t1.collisionId(), event.posX(), event.posY(), event.posZ());
          }
          if constexpr (trackHasCov && TTwoProngFitter) {
            dielectronsExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauzProjected], VarManager::fgValues[VarManager::kVertexingLzProjected], VarManager::fgValues[VarManager::kVertexingLxyProjected]);
          }
        }

        if constexpr (TPairType == VarManager::kDecayToMuMu) {
          twoTrackFilter = a1.isMuonSelected_raw() & a2.isMuonSelected_raw() & fMuonFilterMask;
          if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }

          auto t1 = a1.template reducedmuon_as<TTracks>();
          auto t2 = a2.template reducedmuon_as<TTracks>();
          sign1 = t1.sign();
          sign2 = t2.sign();

          VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
          if constexpr (TTwoProngFitter) {
            VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, fConfigPropToPCA);
          }
          if constexpr (eventHasQvector) {
            VarManager::FillPairVn<TPairType>(t1, t2);
          }

          dimuonList(event.globalIndex(), VarManager::fgValues[VarManager::kMass],
                     VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi],
                     t1.sign() + t2.sign(), twoTrackFilter, 0);
          if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedMuonCollInfo) > 0) {
            dileptonInfoList(t1.collisionId(), event.posX(), event.posY(), event.posZ());
          }

          if constexpr (TTwoProngFitter) {
            dimuonsExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingLz], VarManager::fgValues[VarManager::kVertexingLxy]);
            if (fConfigFlatTables.value) {
              dimuonAllList(event.posX(), event.posY(), event.posZ(), event.numContrib(),
                            -999., -999., -999.,
                            VarManager::fgValues[VarManager::kMass],
                            false,
                            VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), VarManager::fgValues[VarManager::kVertexingChi2PCA],
                            VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingTauzErr],
                            VarManager::fgValues[VarManager::kVertexingTauxy], VarManager::fgValues[VarManager::kVertexingTauxyErr],
                            VarManager::fgValues[VarManager::kCosPointingAngle],
                            VarManager::fgValues[VarManager::kPt1], VarManager::fgValues[VarManager::kEta1], VarManager::fgValues[VarManager::kPhi1], t1.sign(),
                            VarManager::fgValues[VarManager::kPt2], VarManager::fgValues[VarManager::kEta2], VarManager::fgValues[VarManager::kPhi2], t2.sign(),
                            t1.fwdDcaX(), t1.fwdDcaY(), t2.fwdDcaX(), t2.fwdDcaY(),
                            0., 0.,
                            t1.chi2MatchMCHMID(), t2.chi2MatchMCHMID(),
                            t1.chi2MatchMCHMFT(), t2.chi2MatchMCHMFT(),
                            t1.chi2(), t2.chi2(),
                            -999., -999., -999., -999.,
                            -999., -999., -999., -999.,
                            -999., -999., -999., -999.,
                            -999., -999., -999., -999.,
                            t1.isAmbiguous(), t2.isAmbiguous(),
                            VarManager::fgValues[VarManager::kU2Q2], VarManager::fgValues[VarManager::kU3Q3],
                            VarManager::fgValues[VarManager::kR2EP_AB], VarManager::fgValues[VarManager::kR2SP_AB], VarManager::fgValues[VarManager::kCentFT0C],
                            VarManager::fgValues[VarManager::kCos2DeltaPhi], VarManager::fgValues[VarManager::kCos3DeltaPhi],
                            VarManager::fgValues[VarManager::kCORR4POI], VarManager::fgValues[VarManager::kCORR2POI], VarManager::fgValues[VarManager::kM01POI], VarManager::fgValues[VarManager::kM0111POI], VarManager::fgValues[VarManager::kMultDimuons],
                            VarManager::fgValues[VarManager::kVertexingPz], VarManager::fgValues[VarManager::kVertexingSV]);
            }
            if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedMuonCollInfo) > 0) {
              if constexpr (eventHasQvector == true || eventHasQvectorCentr == true) {
                dileptonFlowList(t1.collisionId(), VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kCentFT0C],
                                 VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(),
                                 VarManager::fgValues[VarManager::kU2Q2], VarManager::fgValues[VarManager::kR2SP_AB], VarManager::fgValues[VarManager::kR2SP_AC], VarManager::fgValues[VarManager::kR2SP_BC],
                                 VarManager::fgValues[VarManager::kU3Q3], VarManager::fgValues[VarManager::kR3SP],
                                 VarManager::fgValues[VarManager::kCos2DeltaPhi], VarManager::fgValues[VarManager::kR2EP_AB], VarManager::fgValues[VarManager::kR2EP_AC], VarManager::fgValues[VarManager::kR2EP_BC],
                                 VarManager::fgValues[VarManager::kCos3DeltaPhi], VarManager::fgValues[VarManager::kR3EP],
                                 VarManager::fgValues[VarManager::kCORR4POI], VarManager::fgValues[VarManager::kCORR2POI], VarManager::fgValues[VarManager::kM01POI], VarManager::fgValues[VarManager::kM0111POI],
                                 VarManager::fgValues[VarManager::kCORR2REF], VarManager::fgValues[VarManager::kCORR4REF], VarManager::fgValues[VarManager::kM11REF], VarManager::fgValues[VarManager::kM11REF],
                                 VarManager::fgValues[VarManager::kMultDimuons], VarManager::fgValues[VarManager::kMultA]);
              }
            }
          }
        }
        // TODO: the model for the electron-muon combination has to be thought through
        /*if constexpr (TPairType == VarManager::kElectronMuon) {
          twoTrackFilter = a1.isBarrelSelected_raw() & a1.isBarrelSelectedPrefilter_raw() & a2.isMuonSelected_raw() & fTwoTrackFilterMask;
        }*/

        // Fill histograms
        bool isAmbiInBunch = false;
        bool isAmbiOutOfBunch = false;
        for (int icut = 0; icut < ncuts; icut++) {
          if (twoTrackFilter & (uint32_t(1) << icut)) {
            isAmbiInBunch = (twoTrackFilter & (uint32_t(1) << 28)) || (twoTrackFilter & (uint32_t(1) << 29));
            isAmbiOutOfBunch = (twoTrackFilter & (uint32_t(1) << 30)) || (twoTrackFilter & (uint32_t(1) << 31));
            if (sign1 * sign2 < 0) {
              fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
              if (isAmbiInBunch) {
                fHistMan->FillHistClass(histNames[icut][3 + histIdxOffset].Data(), VarManager::fgValues);
              }
              if (isAmbiOutOfBunch) {
                fHistMan->FillHistClass(histNames[icut][3 + histIdxOffset + 3].Data(), VarManager::fgValues);
              }
            } else {
              if (sign1 > 0) {
                fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
                if (isAmbiInBunch) {
                  fHistMan->FillHistClass(histNames[icut][4 + histIdxOffset].Data(), VarManager::fgValues);
                }
                if (isAmbiOutOfBunch) {
                  fHistMan->FillHistClass(histNames[icut][4 + histIdxOffset + 3].Data(), VarManager::fgValues);
                }
              } else {
                fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
                if (isAmbiInBunch) {
                  fHistMan->FillHistClass(histNames[icut][5 + histIdxOffset].Data(), VarManager::fgValues);
                }
                if (isAmbiOutOfBunch) {
                  fHistMan->FillHistClass(histNames[icut][5 + histIdxOffset + 3].Data(), VarManager::fgValues);
                }
              }
            }
            for (unsigned int iPairCut = 0; iPairCut < fPairCuts.size(); iPairCut++) {
              AnalysisCompositeCut cut = fPairCuts.at(iPairCut);
              if (!(cut.IsSelected(VarManager::fgValues))) // apply pair cuts
                continue;
              if (sign1 * sign2 < 0) {
                fHistMan->FillHistClass(histNames[ncuts + icut * ncuts + iPairCut][0].Data(), VarManager::fgValues);
              } else {
                if (sign1 > 0) {
                  fHistMan->FillHistClass(histNames[ncuts + icut * ncuts + iPairCut][1].Data(), VarManager::fgValues);
                } else {
                  fHistMan->FillHistClass(histNames[ncuts + icut * ncuts + iPairCut][2].Data(), VarManager::fgValues);
                }
              }
            } // end loop (pair cuts)
          }
        } // end loop (cuts)
      }   // end loop over pairs of track associations
    }     // end loop over events
  }

  template <int TPairType, uint32_t TEventFillMap, typename TAssoc1, typename TAssoc2, typename TTracks1, typename TTracks2>
  void runMixedPairing(TAssoc1 const& assocs1, TAssoc2 const& assocs2, TTracks1 const& /*tracks1*/, TTracks2 const& /*tracks2*/)
  {
    std::map<int, std::vector<TString>> histNames = fTrackHistNames;
    int pairSign = 0;
    int ncuts = 0;
    uint32_t twoTrackFilter = 0;
    for (auto& a1 : assocs1) {
      for (auto& a2 : assocs2) {
        if constexpr (TPairType == VarManager::kDecayToEE) {
          twoTrackFilter = a1.isBarrelSelected_raw() & a2.isBarrelSelected_raw() & a1.isBarrelSelectedPrefilter_raw() & a2.isBarrelSelectedPrefilter_raw() & fTrackFilterMask;
          if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }
          auto t1 = a1.template reducedtrack_as<TTracks1>();
          auto t2 = a1.template reducedtrack_as<TTracks2>();
          VarManager::FillPairME<TPairType>(t1, t2);
          if constexpr ((TEventFillMap & VarManager::ObjTypes::ReducedEventQvector) > 0) {
            VarManager::FillPairVn<TPairType>(t1, t2);
          }
          pairSign = t1.sign() + t2.sign();
          ncuts = fNCutsBarrel;
        }
        if constexpr (TPairType == VarManager::kDecayToMuMu) {
          twoTrackFilter = a1.isMuonSelected_raw() & a2.isMuonSelected_raw() & fMuonFilterMask;
          if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }
          auto t1 = a1.template reducedmuon_as<TTracks1>();
          auto t2 = a1.template reducedmuon_as<TTracks2>();
          VarManager::FillPairME<TPairType>(t1, t2);
          if constexpr ((TEventFillMap & VarManager::ObjTypes::ReducedEventQvector) > 0) {
            VarManager::FillPairVn<TPairType>(t1, t2);
          }
          ncuts = fNCutsMuon;
          histNames = fMuonHistNames;
        }
        /*if constexpr (TPairType == VarManager::kElectronMuon) {
          twoTrackFilter = a1.isBarrelSelected_raw() & a1.isBarrelSelectedPrefilter_raw() & a2.isMuonSelected_raw() & fTrackFilterMask;
        }*/

        for (int icut = 0; icut < ncuts; icut++) {
          if (!(twoTrackFilter & (uint32_t(1) << icut))) {
            continue; // cut not passed
          }
          if (pairSign == 0) {
            fHistMan->FillHistClass(histNames[icut][3].Data(), VarManager::fgValues);
          } else {
            if (pairSign > 0) {
              fHistMan->FillHistClass(histNames[icut][4].Data(), VarManager::fgValues);
            } else {
              fHistMan->FillHistClass(histNames[icut][5].Data(), VarManager::fgValues);
            }
          }
        } // end for (cuts)
      }   // end for (track2)
    }     // end for (track1)
  }

  // barrel-barrel and muon-muon event mixing
  template <int TPairType, uint32_t TEventFillMap, typename TEvents, typename TAssocs, typename TTracks>
  void runSameSideMixing(TEvents& events, TAssocs const& assocs, TTracks const& tracks, Preslice<TAssocs>& preSlice)
  {
    events.bindExternalIndices(&assocs);
    int mixingDepth = fConfigMixingDepth.value;
    for (auto& [event1, event2] : selfCombinations(hashBin, mixingDepth, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event1, VarManager::fgValues);

      auto assocs1 = assocs.sliceBy(preSlice, event1.globalIndex());
      assocs1.bindExternalIndices(&events);

      auto assocs2 = assocs.sliceBy(preSlice, event2.globalIndex());
      assocs2.bindExternalIndices(&events);

      runMixedPairing<TPairType, TEventFillMap>(assocs1, assocs2, tracks, tracks);
    } // end event loop
  }

  void processAllSkimmed(MyEventsVtxCovSelected const& events,
                         soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs, MyBarrelTracksWithCovWithAmbiguities const& barrelTracks,
                         soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCov const& muons)
  {
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithCov, gkMuonFillMapWithCov>(events, muonAssocsPerCollision, muonAssocs, muons);
    // runSameEventPairing<true, VarManager::kElectronMuon, gkEventFillMap, gkTrackFillMap>(event, tracks, muons);
  }

  void processBarrelOnlySkimmed(MyEventsVtxCovSelected const& events,
                                soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                                MyBarrelTracksWithCovWithAmbiguities const& barrelTracks)
  {
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processMuonOnlySkimmed(MyEventsVtxCovSelected const& events,
                              soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCov const& muons)
  {
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithCov, gkMuonFillMapWithCov>(events, muonAssocsPerCollision, muonAssocs, muons);
  }

  void processMixingAllSkimmed(soa::Filtered<MyEventsHashSelected>& events,
                               soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& trackAssocs, MyBarrelTracksWithCov const& tracks,
                               soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCov const& muons)
  {
    runSameSideMixing<pairTypeEE, gkEventFillMap>(events, trackAssocs, tracks, trackAssocsPerCollision);
    runSameSideMixing<pairTypeMuMu, gkEventFillMap>(events, muonAssocs, muons, muonAssocsPerCollision);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processAllSkimmed, "Run all types of pairing, with skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlySkimmed, "Run barrel only pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMuonOnlySkimmed, "Run muon only pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMixingAllSkimmed, "Run all types of mixed pairing, with skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", false);
};

// Combines dileptons with barrel or muon tracks for either resonance or correlation analyses
// Dileptons produced with all the selection cuts specified in the same-event pairing task are combined with the
//   tracks passing the fConfigTrackCut cut. The dileptons cuts from the same-event pairing task are auto-detected
struct AnalysisDileptonTrack {
  Produces<aod::BmesonCandidates> BmesonsTable;
  OutputObj<THashList> fOutputList{"output"};

  Configurable<string> fConfigTrackCut{"cfgTrackCut", "kaonPID", "Cut for the track to be correlated with the dileptons"};
  Configurable<float> fConfigDileptonLowMass{"cfgDileptonLowMass", 2.8, "Low mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonHighMass{"cfgDileptonHighMass", 3.2, "High mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonpTCut{"cfgDileptonpTCut", 0.0, "pT cut for dileptons used in the triplet vertexing"};
  Configurable<float> fConfigDileptonLxyCut{"cfgDileptonLxyCut", 0.0, "Lxy cut for dileptons used in the triplet vertexing"};
  Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};

  Configurable<std::string> fConfigHistogramSubgroups{"cfgDileptonTrackHistogramsSubgroups", "invmass,vertexing", "Comma separated list of dilepton-track histogram subgroups"};
  Configurable<int> fConfigMixingDepth{"cfgMixingDepth", 5, "Event mixing pool depth"};

  Configurable<bool> fConfigUseRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<std::string> fConfigGRPmagPath{"cfgGrpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};

  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.
  int fNCuts;
  int fTrackCutBit;
  std::map<int, TString> fHistNamesDileptonTrack;
  std::map<int, TString> fHistNamesDileptons;
  std::map<int, TString> fHistNamesME;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  // TODO: The filter expressions seem to always use the default value of configurables, not the values from the actual configuration file
  Filter eventFilter = aod::dqanalysisflags::isEventSelected == uint32_t(1);
  Filter dileptonFilter = aod::reducedpair::pt > fConfigDileptonpTCut&& aod::reducedpair::mass > fConfigDileptonLowMass&& aod::reducedpair::mass<fConfigDileptonHighMass && aod::reducedpair::sign == 0 && aod::reducedpair::lxy> fConfigDileptonLxyCut;
  Filter filterBarrel = aod::dqanalysisflags::isBarrelSelected > uint32_t(0);
  Filter filterMuon = aod::dqanalysisflags::isMuonSelected > uint32_t(0);

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair; // fill map

  // use two values array to avoid mixing up the quantities
  float* fValuesDilepton;
  float* fValuesHadron;
  HistogramManager* fHistMan;

  NoBinningPolicy<aod::dqanalysisflags::MixingHash> fHashBin;

  void init(o2::framework::InitContext& context)
  {
    bool isBarrel = context.mOptions.get<bool>("processBarrelSkimmed");
    bool isBarrelME = context.mOptions.get<bool>("processBarrelMixedEvent");
    bool isMuon = context.mOptions.get<bool>("processMuonSkimmed");
    bool isMuonME = context.mOptions.get<bool>("processMuonMixedEvent");
    bool isAnyProcessEnabled = isBarrel || isBarrelME || isMuon || isMuonME;
    bool isDummy = context.mOptions.get<bool>("processDummy");
    if (isDummy) {
      if (isAnyProcessEnabled) {
        LOG(warning) << "Dummy function is enabled even if there are normal process functions running! Fix your config!" << endl;
      } else {
        LOG(info) << "Dummy function is enabled. Skipping the rest of the init function" << endl;
        return;
      }
    }

    fCurrentRun = 0;
    fValuesDilepton = new float[VarManager::kNVars];
    fValuesHadron = new float[VarManager::kNVars];
    fTrackCutBit = -1;
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // For each track/muon selection used to produce dileptons, create a separate histogram directory using the
    // name of the track/muon cut.
    // Also, create a map which will hold the name of the histogram directories so they can be accessed directly in the pairing loop
    if (isBarrel || isMuon) {
      // get the list of single track and muon cuts computed in the dedicated tasks upstream
      string tempCutsSingle;
      if (isBarrel) {
        getTaskOptionValue<string>(context, "analysis-track-selection", "cfgTrackCuts", tempCutsSingle, false);
      } else {
        getTaskOptionValue<string>(context, "analysis-muon-selection", "cfgMuonCuts", tempCutsSingle, false);
      }
      TString tempCutsSingleStr = tempCutsSingle;
      TObjArray* objArraySingleCuts = nullptr;
      if (!tempCutsSingleStr.IsNull()) {
        objArraySingleCuts = tempCutsSingleStr.Tokenize(",");
      }
      if (objArraySingleCuts->FindObject(fConfigTrackCut.value.data()) == nullptr) {
        LOG(fatal) << " Track cut chosen for the correlation task was not computed in the single-track task! Check it out!";
      }
      // get the cuts employed for same-event pairing
      string tempCutsPair;
      if (isBarrel) {
        getTaskOptionValue<string>(context, "analysis-same-event-pairing", "cfgTrackCuts", tempCutsPair, false);
      } else {
        getTaskOptionValue<string>(context, "analysis-same-event-pairing", "cfgMuonCuts", tempCutsPair, false);
      }
      TString tempCutsPairStr = tempCutsPair;
      if (!tempCutsSingleStr.IsNull() && !tempCutsPairStr.IsNull()) {
        std::unique_ptr<TObjArray> objArray(tempCutsPairStr.Tokenize(","));
        fNCuts = objArray->GetEntries();
        for (int icut = 0; icut < objArraySingleCuts->GetEntries(); ++icut) {
          TString tempStr = objArraySingleCuts->At(icut)->GetName();
          if (objArray->FindObject(tempStr.Data()) != nullptr) {
            fHistNamesDileptonTrack[icut] = Form("DileptonTrack_%s_%s", tempStr.Data(), fConfigTrackCut.value.data());
            fHistNamesDileptons[icut] = Form("DileptonsSelected_%s", tempStr.Data());
            DefineHistograms(fHistMan, fHistNamesDileptonTrack[icut], fConfigHistogramSubgroups.value.data()); // define dilepton-track histograms
            DefineHistograms(fHistMan, fHistNamesDileptons[icut], "barrel,vertexing");                         // define dilepton histograms
            if (isBarrelME || isMuonME) {
              fHistNamesME[icut] = Form("DileptonTrackME_%s", tempStr.Data());
              DefineHistograms(fHistMan, fHistNamesME[icut], "mixedevent"); // define ME histograms
            }
          }
          if (tempStr.CompareTo(fConfigTrackCut.value.data()) == 0) {
            fTrackCutBit = icut; // the bit correspoding to the track to be combined with dileptons
          }
        }
      }
    }
    if (fHistNamesDileptons.size() == 0) {
      LOG(fatal) << " No valid dilepton cuts ";
    }

    if (context.mOptions.get<bool>("processBarrelMixedEvent")) {
      DefineHistograms(fHistMan, "DileptonTrackME", "mixedevent"); // define all histograms
    }

    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  // init parameters from CCDB
  void initParamsFromCCDB(uint64_t timestamp)
  {
    if (fConfigUseRemoteField.value) {
      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fConfigGRPmagPath.value, timestamp);
      float magField = 0.0;
      if (grpmag != nullptr) {
        magField = grpmag->getNominalL3Field();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", timestamp);
      }
      if (fConfigUseKFVertexing.value) {
        VarManager::SetupThreeProngKFParticle(magField);
      } else {
        VarManager::SetupThreeProngDCAFitter(); // TODO: get these parameters from Configurables
      }
    } else {
      if (fConfigUseKFVertexing.value) {
        VarManager::SetupThreeProngKFParticle(fConfigMagField.value);
      } else {
        VarManager::SetupThreeProngDCAFitter(); // TODO: get these parameters from Configurables
      }
    }
  }

  // Template function to run pair - hadron combinations
  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks, typename TTrackAssocs, typename TDileptons>
  void runDileptonHadron(TEvent const& event, TTrackAssocs const& assocs, TTracks const& tracks, TDileptons const& dileptons)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesHadron);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    VarManager::FillEvent<TEventFillMap>(event, fValuesHadron);
    VarManager::FillEvent<TEventFillMap>(event, fValuesDilepton);

    for (auto dilepton : dileptons) {
      // get full track info of tracks based on the index
      auto lepton1 = tracks.rawIteratorAt(dilepton.index0Id());
      auto lepton2 = tracks.rawIteratorAt(dilepton.index1Id());
      // Check that the dilepton has zero charge
      if (dilepton.sign() != 0) {
        continue;
      }

      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesDilepton);
      for (int icut = 0; icut < fNCuts; icut++) {
        if (dilepton.filterMap_bit(icut)) {
          fHistMan->FillHistClass(fHistNamesDileptons[icut].Data(), fValuesDilepton);
        }
      }

      // loop over hadrons
      for (auto& assoc : assocs) {
        if constexpr (TCandidateType == VarManager::kBtoJpsiEEK) {
          if (!assoc.isBarrelSelected_bit(fTrackCutBit)) {
            continue;
          }
          auto track = assoc.template reducedtrack_as<TTracks>();
          if (track.globalIndex() == dilepton.index0Id() || track.globalIndex() == dilepton.index1Id()) {
            continue;
          }
          VarManager::FillDileptonHadron(dilepton, track, fValuesHadron);
          VarManager::FillDileptonTrackVertexing<TCandidateType, TEventFillMap, TTrackFillMap>(event, lepton1, lepton2, track, fValuesHadron);
        }
        if constexpr (TCandidateType == VarManager::kBcToThreeMuons) {
          if (!assoc.isMuonSelected_bit(fTrackCutBit)) {
            continue;
          }
          auto track = assoc.template reducedmuon_as<TTracks>();
          if (track.globalIndex() == dilepton.index0Id() || track.globalIndex() == dilepton.index1Id()) {
            continue;
          }

          VarManager::FillDileptonHadron(dilepton, track, fValuesHadron);
          VarManager::FillDileptonTrackVertexing<TCandidateType, TEventFillMap, TTrackFillMap>(event, lepton1, lepton2, track, fValuesHadron);
        }

        for (int icut = 0; icut < fNCuts; icut++) {
          if (dilepton.filterMap_bit(icut)) {
            fHistMan->FillHistClass(fHistNamesDileptonTrack[icut].Data(), fValuesHadron);
          }
        }
        // table to be written out for ML analysis
        BmesonsTable(fValuesHadron[VarManager::kPairMass], fValuesHadron[VarManager::kPairPt], fValuesHadron[VarManager::kVertexingLxy], fValuesHadron[VarManager::kVertexingLxyz], fValuesHadron[VarManager::kVertexingLz], fValuesHadron[VarManager::kVertexingTauxy], fValuesHadron[VarManager::kVertexingTauz], fValuesHadron[VarManager::kCosPointingAngle], fValuesHadron[VarManager::kVertexingChi2PCA]);
      }
    }
  }

  Preslice<aod::ReducedTracksAssoc> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;
  Preslice<MyDielectronCandidates> dielectronsPerCollision = aod::reducedpair::reducedeventId;

  void processBarrelSkimmed(soa::Filtered<MyEventsVtxCovSelected> const& events,
                            soa::Filtered<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> const& assocs,
                            MyBarrelTracksWithCov const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {
    // set up KF or DCAfitter
    if (events.size() == 0) {
      return;
    }
    if (fCurrentRun != events.begin().runNumber()) { // start: runNumber
      initParamsFromCCDB(events.begin().timestamp());
      fCurrentRun = events.begin().runNumber();
    } // end: runNumber
    for (auto& event : events) {
      auto groupedBarrelAssocs = assocs.sliceBy(trackAssocsPerCollision, event.globalIndex());
      auto groupedDielectrons = dileptons.sliceBy(dielectronsPerCollision, event.globalIndex());
      runDileptonHadron<VarManager::kBtoJpsiEEK, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, groupedBarrelAssocs, tracks, groupedDielectrons);
    }
  }

  Preslice<aod::ReducedMuonsAssoc> muonAssocsPerCollision = aod::reducedtrack_association::reducedeventId;
  Preslice<MyDimuonCandidates> dimuonsPerCollision = aod::reducedpair::reducedeventId;

  void processMuonSkimmed(soa::Filtered<MyEventsVtxCovSelected> const& events,
                          soa::Filtered<soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts>> const& assocs,
                          MyMuonTracksWithCov const& tracks, soa::Filtered<MyDimuonCandidates> const& dileptons)
  {
    // set up KF or DCAfitter
    if (events.size() == 0) {
      return;
    }
    if (fCurrentRun != events.begin().runNumber()) { // start: runNumber
      initParamsFromCCDB(events.begin().timestamp());
      fCurrentRun = events.begin().runNumber();
    } // end: runNumber
    for (auto& event : events) {
      auto groupedMuonAssocs = assocs.sliceBy(muonAssocsPerCollision, event.globalIndex());
      auto groupedDimuons = dileptons.sliceBy(dimuonsPerCollision, event.globalIndex());
      runDileptonHadron<VarManager::kBcToThreeMuons, gkEventFillMapWithCov, gkMuonFillMapWithCov>(event, groupedMuonAssocs, tracks, groupedDimuons);
    }
  }

  void processBarrelMixedEvent(soa::Filtered<MyEventsHashSelected>& events,
                               soa::Filtered<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> const& assocs,
                               MyBarrelTracksWithCov const&, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {
    if (events.size() == 0) {
      return;
    }
    events.bindExternalIndices(&dileptons);
    events.bindExternalIndices(&assocs);

    for (auto& [event1, event2] : selfCombinations(fHashBin, fConfigMixingDepth.value, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<gkEventFillMap>(event1, VarManager::fgValues);

      auto evDileptons = dileptons.sliceBy(dielectronsPerCollision, event1.globalIndex());
      evDileptons.bindExternalIndices(&events);

      auto evAssocs = assocs.sliceBy(trackAssocsPerCollision, event2.globalIndex());
      evAssocs.bindExternalIndices(&events);

      for (auto& assoc : evAssocs) {
        if (!assoc.isBarrelSelected_bit(fTrackCutBit)) {
          continue;
        }
        auto track = assoc.template reducedtrack_as<MyBarrelTracksWithCov>();

        for (auto dilepton : evDileptons) {
          VarManager::FillDileptonHadron(dilepton, track, VarManager::fgValues);
          for (int icut = 0; icut < fNCuts; icut++) {
            if (dilepton.filterMap_bit(icut)) {
              fHistMan->FillHistClass(fHistNamesME[icut].Data(), VarManager::fgValues);
            }
          }
        } // end for (dileptons)
      }   // end for (assocs)
    }     // end event loop
  }

  void processMuonMixedEvent(soa::Filtered<MyEventsHashSelected>& events,
                             soa::Filtered<soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts>> const& assocs,
                             MyMuonTracksWithCov const&, soa::Filtered<MyDimuonCandidates> const& dileptons)
  {
    if (events.size() == 0) {
      return;
    }
    events.bindExternalIndices(&dileptons);
    events.bindExternalIndices(&assocs);

    for (auto& [event1, event2] : selfCombinations(fHashBin, fConfigMixingDepth.value, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<gkEventFillMap>(event1, VarManager::fgValues);

      auto evDileptons = dileptons.sliceBy(dimuonsPerCollision, event1.globalIndex());
      evDileptons.bindExternalIndices(&events);

      auto evAssocs = assocs.sliceBy(muonAssocsPerCollision, event2.globalIndex());
      evAssocs.bindExternalIndices(&events);

      for (auto& assoc : evAssocs) {

        if (!assoc.isMuonSelected_bit(fTrackCutBit)) {
          continue;
        }
        auto track = assoc.template reducedmuon_as<MyMuonTracksWithCov>();

        for (auto dilepton : evDileptons) {
          VarManager::FillDileptonHadron(dilepton, track, VarManager::fgValues);
          for (int icut = 0; icut < fNCuts; icut++) {
            if (dilepton.filterMap_bit(icut)) {
              fHistMan->FillHistClass(fHistNamesME[icut].Data(), VarManager::fgValues);
            }
          }
        } // end for (dileptons)
      }   // end for (assocs)
    }     // end event loop
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisDileptonTrack, processBarrelSkimmed, "Run barrel dilepton-track pairing, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processMuonSkimmed, "Run muon dilepton-track pairing, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processBarrelMixedEvent, "Run barrel dilepton-hadron mixed event pairing", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processMuonMixedEvent, "Run muon dilepton-hadron mixed event pairing", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
    adaptAnalysisTask<AnalysisMuonSelection>(cfgc),
    adaptAnalysisTask<AnalysisPrefilterSelection>(cfgc),
    adaptAnalysisTask<AnalysisSameEventPairing>(cfgc),
    adaptAnalysisTask<AnalysisDileptonTrack>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses, const char* histGroups)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //  The histogram classes are provided in the histClasses string, separated by semicolon ";"
  //  The histogram classes and their components histograms are defined below depending on the name of the histogram class
  //
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    TString histName = histGroups;
    // NOTE: The level of detail for histogramming can be controlled via configurables
    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", histName);
    }

    if (classStr.Contains("SameBunchCorrelations")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "two-collisions", histName);
    }

    if (classStr.Contains("Track") && !classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
        if (classStr.Contains("PIDCalibElectron")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_electron");
        }
        if (classStr.Contains("PIDCalibPion")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_pion");
        }
        if (classStr.Contains("PIDCalibProton")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_proton");
        }
        if (classStr.Contains("Ambiguity")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", Form("%s,ambiguity", histName.Data()));
        }
      }
      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
      }
    }

    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }

    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "barrel,vertexing");
    }

    if (classStr.Contains("DileptonTrack") && !classStr.Contains("ME")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-track", histName);
    }

    if (classStr.Contains("DileptonTrackME")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-track", "mixedevent");
    }

    if (classStr.Contains("HadronsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
    }

    if (classStr.Contains("DileptonHadronInvMass")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-mass");
    }

    if (classStr.Contains("DileptonHadronCorrelation")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-correlation");
    }
  } // end loop over histogram classes
}
