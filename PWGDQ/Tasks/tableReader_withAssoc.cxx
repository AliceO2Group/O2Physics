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

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <numeric>
#include <vector>
#include <algorithm>
#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TString.h>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/Configurable.h"
#include "Framework/OutputObjHeader.h"
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
#include "ITSMFTBase/DPLAlpideParam.h"
#include "Common/CCDB/EventSelectionParams.h"

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
DECLARE_SOA_COLUMN(MuonAmbiguityInBunch, muonAmbiguityInBunch, int8_t);              //! Muon track in-bunch ambiguity
DECLARE_SOA_COLUMN(MuonAmbiguityOutOfBunch, muonAmbiguityOutOfBunch, int8_t);        //! Muon track out of bunch ambiguity
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
DECLARE_SOA_TABLE(MuonAmbiguities, "AOD", "DQMUONAMB", dqanalysisflags::MuonAmbiguityInBunch, dqanalysisflags::MuonAmbiguityOutOfBunch);         //!  joinable to ReducedMuonTracks
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTER", dqanalysisflags::IsBarrelSelectedPrefilter);                                                  //!  joinable to ReducedTracksAssoc
DECLARE_SOA_TABLE(BmesonCandidates, "AOD", "DQBMESONS", dqanalysisflags::massBcandidate, dqanalysisflags::pTBcandidate, dqanalysisflags::LxyBcandidate, dqanalysisflags::LxyzBcandidate, dqanalysisflags::LzBcandidate, dqanalysisflags::TauxyBcandidate, dqanalysisflags::TauzBcandidate, dqanalysisflags::CosPBcandidate, dqanalysisflags::Chi2Bcandidate);
} // namespace o2::aod

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsMultExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;
using MyEventsZdc = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedZdcs>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsMultExtraSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll, aod::EventCuts>;
using MyEventsHashSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts>;
using MyEventsVtxCovSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvector>;
using MyEventsVtxCovZdcSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::ReducedZdcs, aod::EventCuts>;
using MyEventsQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector>;
using MyEventsHashSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvector>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksWithAmbiguities = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelAmbiguities>;
using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksWithCovWithAmbiguities = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelAmbiguities>;
using MyBarrelTracksWithCovWithAmbiguitiesWithColl = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelAmbiguities, aod::ReducedTracksBarrelInfo>;
using MyDielectronCandidates = soa::Join<aod::Dielectrons, aod::DielectronsExtra>;
using MyDitrackCandidates = soa::Join<aod::Ditracks, aod::DitracksExtra>;
using MyDimuonCandidates = soa::Join<aod::Dimuons, aod::DimuonsExtra>;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyMuonTracksWithCovWithAmbiguities = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonAmbiguities>;
using MyMuonTracksSelectedWithColl = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsInfo, aod::MuonTrackCuts>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkEventFillMapWithZdc = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedZdc;
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;
constexpr static uint32_t gkEventFillMapWithCovZdc = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ReducedZdc;
constexpr static uint32_t gkEventFillMapWithMultExtra = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventMultExtra;
// constexpr static uint32_t gkEventFillMapWithQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector;
// constexpr static uint32_t gkEventFillMapWithCovQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithCovWithColl = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID | VarManager::ObjTypes::ReducedTrackCollInfo;

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
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
  Configurable<int> fConfigITSROFrameStartBorderMargin{"ITSROFrameStartBorderMargin", -1, "Number of bcs at the start of ITS RO Frame border. Take from CCDB if -1"};
  Configurable<int> fConfigITSROFrameEndBorderMargin{"ITSROFrameEndBorderMargin", -1, "Number of bcs at the end of ITS RO Frame border. Take from CCDB if -1"};
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;
  AnalysisCompositeCut* fEventCut;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

  std::map<int64_t, bool> fSelMap;                     // key: reduced event global index, value: event selection decision
  std::map<uint64_t, std::vector<int64_t>> fBCCollMap; // key: global BC, value: vector of reduced event global indices
  int fCurrentRun;

  void init(o2::framework::InitContext&)
  {
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
    DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;SameBunchCorrelations", fConfigAddEventHistogram.value.data()); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());                                                                             // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

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
      std::map<string, string> metadataRCT, header;
      header = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", events.begin().runNumber()), metadataRCT, -1);
      uint64_t sor = std::atol(header["SOR"].c_str());
      uint64_t eor = std::atol(header["EOR"].c_str());
      VarManager::SetSORandEOR(sor, eor);
      LOG(info) << "============================= SOR / EOR :: " << sor << " / " << eor;

      auto alppar = fCCDB->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", events.begin().timestamp());
      EventSelectionParams* par = fCCDB->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", events.begin().timestamp());
      int itsROFrameStartBorderMargin = fConfigITSROFrameStartBorderMargin < 0 ? par->fITSROFrameStartBorderMargin : fConfigITSROFrameStartBorderMargin;
      int itsROFrameEndBorderMargin = fConfigITSROFrameEndBorderMargin < 0 ? par->fITSROFrameEndBorderMargin : fConfigITSROFrameEndBorderMargin;
      VarManager::SetITSROFBorderselection(alppar->roFrameBiasInBC, alppar->roFrameLengthInBC, itsROFrameStartBorderMargin, itsROFrameEndBorderMargin);
      LOGP(debug, "==============++++++++++++========== roBias / roLength / start / end :: {} / {} / {} / {}", alppar->roFrameBiasInBC, alppar->roFrameLengthInBC, itsROFrameStartBorderMargin, itsROFrameEndBorderMargin);

      fCurrentRun = events.begin().runNumber();
    }

    fSelMap.clear();
    fBCCollMap.clear();

    for (auto& event : events) {
      // Reset the fValues array and fill event observables
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEvent<TEventFillMap>(event);

      bool decision = false;
      fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues); // automatically fill all the histograms in the class Event
      if (fEventCut->IsSelected(VarManager::fgValues)) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
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
  void processSkimmedWithZdc(MyEventsZdc const& events)
  {
    runEventSelection<gkEventFillMapWithZdc>(events);
    publishSelections<gkEventFillMapWithZdc>(events);
  }
  void processSkimmedWithMultExtra(MyEventsMultExtra const& events)
  {
    runEventSelection<gkEventFillMapWithMultExtra>(events);
    publishSelections<gkEventFillMapWithMultExtra>(events);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedWithZdc, "Run event selection on DQ skimmed events, with ZDC", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedWithMultExtra, "Run event selection on DQ skimmed events, with mult extra", false);
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
  // Track related options
  Configurable<bool> fPropTrack{"cfgPropTrack", true, "Propgate tracks to associated collision to recalculate DCA and momentum vector"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fTrackCuts;

  int fCurrentRun; // current run kept to detect run changes and trigger loading params from CCDB

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
    fCCDBApi.init(fConfigCcdbUrl.value);
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTracks>
  void runTrackSelection(ReducedTracksAssoc const& assocs, TEvents const& events, TTracks const& tracks)
  {
    fNAssocsInBunch.clear();
    fNAssocsOutOfBunch.clear();

    if (events.size() > 0 && fCurrentRun != events.begin().runNumber()) {
      if (fConfigComputeTPCpostCalib) {
        auto calibList = fCCDB->getForTimeStamp<TList>(fConfigCcdbPathTPC.value, events.begin().timestamp());
        VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
      }

      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, events.begin().timestamp());
      if (grpmag != nullptr) {
        VarManager::SetMagneticField(grpmag->getNominalL3Field());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", events.begin().timestamp());
      }

      std::map<string, string> metadataRCT, header;
      header = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", events.begin().runNumber()), metadataRCT, -1);
      uint64_t sor = std::atol(header["SOR"].c_str());
      uint64_t eor = std::atol(header["EOR"].c_str());
      VarManager::SetSORandEOR(sor, eor);

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
      if (fPropTrack) {
        VarManager::FillTrackCollision<TTrackFillMap>(track, event);
      }
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
  void processSkimmedWithMultExtra(ReducedTracksAssoc const& assocs, MyEventsMultExtraSelected const& events, MyBarrelTracks const& tracks)
  {
    runTrackSelection<gkEventFillMapWithMultExtra, gkTrackFillMap>(assocs, events, tracks);
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
  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmedWithMultExtra, "Run barrel track selection on DQ skimmed track associations, with extra multiplicity tables", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmedWithCov, "Run barrel track selection on DQ skimmed tracks w/ cov matrix associations", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummy, "Dummy function", false);
};

// Produces a table with muon decisions (joinable to the ReducedMuonsAssociations)
// Here one should add all the track cuts needed through the workflow (e.g. cuts for same-event pairing, track for dilepton-track correlations)
struct AnalysisMuonSelection {
  Produces<aod::MuonTrackCuts> muonSel;
  Produces<aod::MuonAmbiguities> muonAmbiguities;
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

  std::map<int64_t, std::vector<int64_t>> fNAssocsInBunch;    // key: muon global index, value: vector of global index for events associated in-bunch (events that have in-bunch pileup or splitting)
  std::map<int64_t, std::vector<int64_t>> fNAssocsOutOfBunch; // key: muon global index, value: vector of global index for events associated out-of-bunch (events that have no in-bunch pileup)

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
      histDirNames += "TrackMuon_AmbiguityInBunch;TrackMuon_AmbiguityOutOfBunch;";

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
  void runMuonSelection(ReducedMuonsAssoc const& assocs, TEvents const& events, TMuons const& muons)
  {
    fNAssocsInBunch.clear();
    fNAssocsOutOfBunch.clear();

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
    muonAmbiguities.reserve(muons.size());
    uint32_t filterMap = 0;
    int iCut = 0;

    for (auto& assoc : assocs) {
      auto event = assoc.template reducedevent_as<TEvents>();
      if (!event.isEventSelected_bit(0)) {
        muonSel(0);
        continue;
      }
      VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
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
    } // end loop over assocs

    // QA the collision-track associations
    if (fConfigQA) {
      for (auto& [trackIdx, evIndices] : fNAssocsInBunch) {
        if (evIndices.size() == 1) {
          continue;
        }
        auto track = muons.rawIteratorAt(trackIdx);
        VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
        VarManager::FillTrack<TMuonFillMap>(track);
        VarManager::fgValues[VarManager::kMuonNAssocsInBunch] = static_cast<float>(evIndices.size());
        fHistMan->FillHistClass("TrackMuon_AmbiguityInBunch", VarManager::fgValues);
      } // end loop over in-bunch ambiguous tracks

      for (auto& [trackIdx, evIndices] : fNAssocsOutOfBunch) {
        if (evIndices.size() == 1) {
          continue;
        }
        auto track = muons.rawIteratorAt(trackIdx);
        VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
        VarManager::FillTrack<TMuonFillMap>(track);
        VarManager::fgValues[VarManager::kMuonNAssocsOutOfBunch] = static_cast<float>(evIndices.size());
        fHistMan->FillHistClass("TrackMuon_AmbiguityOutOfBunch", VarManager::fgValues);
      } // end loop over out-of-bunch ambiguous tracks
    }

    // publish the ambiguity table
    for (auto& track : muons) {
      int8_t nInBunch = 0;
      if (fNAssocsInBunch.find(track.globalIndex()) != fNAssocsInBunch.end()) {
        nInBunch = fNAssocsInBunch[track.globalIndex()].size();
      }
      int8_t nOutOfBunch = 0;
      if (fNAssocsOutOfBunch.find(track.globalIndex()) != fNAssocsOutOfBunch.end()) {
        nOutOfBunch = fNAssocsOutOfBunch[track.globalIndex()].size();
      }
      muonAmbiguities(nInBunch, nOutOfBunch);
    }
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
    if (initContext.mOptions.get<bool>("processDummy")) {
      LOG(info) << "Dummy function enabled. Skipping the rest of init()" << endl;
      return;
    }
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
  Produces<aod::DielectronsInfo> dielectronInfoList;
  Produces<aod::DimuonsExtra> dimuonsExtraList;
  Produces<aod::DielectronsAll> dielectronAllList;
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
  // Configurable<std::string> fConfigAddEventMixingHistogram{"cfgAddEventMixingHistogram", "", "Comma separated list of histograms"};

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
  // Track related options
  Configurable<bool> fPropTrack{"cfgPropTrack", true, "Propgate tracks to associated collision to recalculate DCA and momentum vector"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

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
    fEnableBarrelHistos = context.mOptions.get<bool>("processAllSkimmed") || context.mOptions.get<bool>("processBarrelOnlySkimmed") || context.mOptions.get<bool>("processBarrelOnlyWithCollSkimmed") || context.mOptions.get<bool>("processBarrelOnlySkimmedNoCov");
    fEnableBarrelMixingHistos = context.mOptions.get<bool>("processMixingAllSkimmed") || context.mOptions.get<bool>("processMixingBarrelSkimmed");
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
            names.push_back(Form("PairsMuonSEPM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
            names.push_back(Form("PairsMuonSEPP_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
            names.push_back(Form("PairsMuonSEMM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
            names.push_back(Form("PairsMuonSEPM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
            names.push_back(Form("PairsMuonSEPP_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
            names.push_back(Form("PairsMuonSEMM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
            histNames += Form("%s;%s;%s;", names[(fEnableMuonMixingHistos ? 6 : 3)].Data(), names[(fEnableMuonMixingHistos ? 7 : 4)].Data(), names[(fEnableMuonMixingHistos ? 8 : 5)].Data());
            histNames += Form("%s;%s;%s;", names[(fEnableMuonMixingHistos ? 9 : 6)].Data(), names[(fEnableMuonMixingHistos ? 10 : 7)].Data(), names[(fEnableMuonMixingHistos ? 11 : 8)].Data());
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
    fCCDBApi.init(fConfigCcdbUrl.value);

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

  void initParamsFromCCDB(uint64_t timestamp, int runNumber, bool withTwoProngFitter = true)
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

    std::map<string, string> metadataRCT, header;
    header = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", runNumber), metadataRCT, -1);
    uint64_t sor = std::atol(header["SOR"].c_str());
    uint64_t eor = std::atol(header["EOR"].c_str());
    VarManager::SetSORandEOR(sor, eor);
  }

  // Template function to run same event pairing (barrel-barrel, muon-muon, barrel-muon)
  template <bool TTwoProngFitter, int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTrackAssocs, typename TTracks>
  void runSameEventPairing(TEvents const& events, Preslice<TTrackAssocs>& preslice, TTrackAssocs const& assocs, TTracks const& /*tracks*/)
  {
    if (events.size() > 0) { // Additional protection to avoid crashing of events.begin().runNumber()
      if (fCurrentRun != events.begin().runNumber()) {
        initParamsFromCCDB(events.begin().timestamp(), events.begin().runNumber(), TTwoProngFitter);
        fCurrentRun = events.begin().runNumber();
      }
    }

    TString cutNames = fConfigTrackCuts.value;
    std::map<int, std::vector<TString>> histNames = fTrackHistNames;
    int ncuts = fNCutsBarrel;
    int histIdxOffset = 0;
    if constexpr (TPairType == pairTypeMuMu) {
      cutNames = fConfigMuonCuts.value;
      histNames = fMuonHistNames;
      ncuts = fNCutsMuon;
      if (fEnableMuonMixingHistos) {
        histIdxOffset = 3;
      }
    }
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
    uint32_t dileptonMcDecision = 0; // placeholder, copy of the dqEfficiency.cxx one
    int sign1 = 0;
    int sign2 = 0;
    dielectronList.reserve(1);
    dimuonList.reserve(1);
    dielectronsExtraList.reserve(1);
    dielectronInfoList.reserve(1);
    dimuonsExtraList.reserve(1);
    dileptonInfoList.reserve(1);
    dileptonFlowList.reserve(1);
    if (fConfigFlatTables.value) {
      dielectronAllList.reserve(1);
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

      bool isFirst = true;
      for (auto& [a1, a2] : o2::soa::combinations(groupedAssocs, groupedAssocs)) {
        if constexpr (TPairType == VarManager::kDecayToEE) {
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
          // compute quantities which depend on the associated collision, such as DCA
          if (fPropTrack) {
            VarManager::FillPairCollision<TPairType, TTrackFillMap>(event, t1, t2);
          }
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
            dielectronInfoList(t1.collisionId(), t1.trackId(), t2.trackId());
            dileptonInfoList(t1.collisionId(), event.posX(), event.posY(), event.posZ());
          }
          if constexpr (trackHasCov && TTwoProngFitter) {
            dielectronsExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauzProjected], VarManager::fgValues[VarManager::kVertexingLzProjected], VarManager::fgValues[VarManager::kVertexingLxyProjected]);
            if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelPID) > 0) {
              if (fConfigFlatTables.value) {
                dielectronAllList(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), twoTrackFilter, dileptonMcDecision,
                                  t1.pt(), t1.eta(), t1.phi(), t1.tpcNClsCrossedRows(), t1.tpcNClsFound(), t1.tpcChi2NCl(), t1.dcaXY(), t1.dcaZ(), t1.tpcSignal(), t1.tpcNSigmaEl(), t1.tpcNSigmaPi(), t1.tpcNSigmaPr(), t1.beta(), t1.tofNSigmaEl(), t1.tofNSigmaPi(), t1.tofNSigmaPr(),
                                  t2.pt(), t2.eta(), t2.phi(), t2.tpcNClsCrossedRows(), t2.tpcNClsFound(), t2.tpcChi2NCl(), t2.dcaXY(), t2.dcaZ(), t2.tpcSignal(), t2.tpcNSigmaEl(), t2.tpcNSigmaPi(), t2.tpcNSigmaPr(), t2.beta(), t2.tofNSigmaEl(), t2.tofNSigmaPi(), t2.tofNSigmaPr(),
                                  VarManager::fgValues[VarManager::kKFTrack0DCAxyz], VarManager::fgValues[VarManager::kKFTrack1DCAxyz], VarManager::fgValues[VarManager::kKFDCAxyzBetweenProngs], VarManager::fgValues[VarManager::kKFTrack0DCAxy], VarManager::fgValues[VarManager::kKFTrack1DCAxy], VarManager::fgValues[VarManager::kKFDCAxyBetweenProngs],
                                  VarManager::fgValues[VarManager::kKFTrack0DeviationFromPV], VarManager::fgValues[VarManager::kKFTrack1DeviationFromPV], VarManager::fgValues[VarManager::kKFTrack0DeviationxyFromPV], VarManager::fgValues[VarManager::kKFTrack1DeviationxyFromPV],
                                  VarManager::fgValues[VarManager::kKFMass], VarManager::fgValues[VarManager::kKFChi2OverNDFGeo], VarManager::fgValues[VarManager::kVertexingLxyz], VarManager::fgValues[VarManager::kVertexingLxyzOverErr], VarManager::fgValues[VarManager::kVertexingLxy], VarManager::fgValues[VarManager::kVertexingLxyOverErr], VarManager::fgValues[VarManager::kVertexingTauxy], VarManager::fgValues[VarManager::kVertexingTauxyErr], VarManager::fgValues[VarManager::kKFCosPA], VarManager::fgValues[VarManager::kKFJpsiDCAxyz], VarManager::fgValues[VarManager::kKFJpsiDCAxy],
                                  VarManager::fgValues[VarManager::kKFPairDeviationFromPV], VarManager::fgValues[VarManager::kKFPairDeviationxyFromPV],
                                  VarManager::fgValues[VarManager::kKFMassGeoTop], VarManager::fgValues[VarManager::kKFChi2OverNDFGeoTop]);
              }
            }
          }
        }

        if constexpr (TPairType == VarManager::kDecayToMuMu) {
          twoTrackFilter = a1.isMuonSelected_raw() & a2.isMuonSelected_raw() & fMuonFilterMask;
          if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }

          auto t1 = a1.template reducedmuon_as<TTracks>();
          auto t2 = a2.template reducedmuon_as<TTracks>();
          if (t1.matchMCHTrackId() == t2.matchMCHTrackId())
            continue;
          if (t1.matchMFTTrackId() == t2.matchMFTTrackId())
            continue;
          sign1 = t1.sign();
          sign2 = t2.sign();
          // store the ambiguity number of the two dilepton legs in the last 4 digits of the two-track filter
          if (t1.muonAmbiguityInBunch() > 1) {
            twoTrackFilter |= (uint32_t(1) << 28);
          }
          if (t2.muonAmbiguityInBunch() > 1) {
            twoTrackFilter |= (uint32_t(1) << 29);
          }
          if (t1.muonAmbiguityOutOfBunch() > 1) {
            twoTrackFilter |= (uint32_t(1) << 30);
          }
          if (t2.muonAmbiguityOutOfBunch() > 1) {
            twoTrackFilter |= (uint32_t(1) << 31);
          }

          VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
          // compute quantities which depend on the associated collision, such as DCA
          if (fPropTrack) {
            VarManager::FillPairCollision<TPairType, TTrackFillMap>(event, t1, t2);
          }
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
                            (twoTrackFilter & (uint32_t(1) << 28)) || (twoTrackFilter & (uint32_t(1) << 30)), (twoTrackFilter & (uint32_t(1) << 29)) || (twoTrackFilter & (uint32_t(1) << 31)),
                            VarManager::fgValues[VarManager::kU2Q2], VarManager::fgValues[VarManager::kU3Q3],
                            VarManager::fgValues[VarManager::kR2EP_AB], VarManager::fgValues[VarManager::kR2SP_AB], VarManager::fgValues[VarManager::kCentFT0C],
                            VarManager::fgValues[VarManager::kCos2DeltaPhi], VarManager::fgValues[VarManager::kCos3DeltaPhi],
                            VarManager::fgValues[VarManager::kCORR2POI], VarManager::fgValues[VarManager::kCORR4POI], VarManager::fgValues[VarManager::kM01POI], VarManager::fgValues[VarManager::kM0111POI], VarManager::fgValues[VarManager::kMultDimuons],
                            VarManager::fgValues[VarManager::kVertexingPz], VarManager::fgValues[VarManager::kVertexingSV]);
            }
            if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedMuonCollInfo) > 0) {
              if constexpr (eventHasQvector == true || eventHasQvectorCentr == true) {
                dileptonFlowList(t1.collisionId(), VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kCentFT0C],
                                 VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), isFirst,
                                 VarManager::fgValues[VarManager::kU2Q2], VarManager::fgValues[VarManager::kR2SP_AB], VarManager::fgValues[VarManager::kR2SP_AC], VarManager::fgValues[VarManager::kR2SP_BC],
                                 VarManager::fgValues[VarManager::kU3Q3], VarManager::fgValues[VarManager::kR3SP],
                                 VarManager::fgValues[VarManager::kCos2DeltaPhi], VarManager::fgValues[VarManager::kR2EP_AB], VarManager::fgValues[VarManager::kR2EP_AC], VarManager::fgValues[VarManager::kR2EP_BC],
                                 VarManager::fgValues[VarManager::kCos3DeltaPhi], VarManager::fgValues[VarManager::kR3EP],
                                 VarManager::fgValues[VarManager::kCORR2POI], VarManager::fgValues[VarManager::kCORR4POI], VarManager::fgValues[VarManager::kM01POI], VarManager::fgValues[VarManager::kM0111POI],
                                 VarManager::fgValues[VarManager::kCORR2REF], VarManager::fgValues[VarManager::kCORR4REF], VarManager::fgValues[VarManager::kM11REF], VarManager::fgValues[VarManager::kM1111REF],
                                 VarManager::fgValues[VarManager::kMultDimuons], VarManager::fgValues[VarManager::kMultA]);
              }
            }
          }
          if (t1.sign() != t2.sign()) {
            isFirst = false;
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
          auto t2 = a2.template reducedtrack_as<TTracks2>();
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
          auto t2 = a2.template reducedmuon_as<TTracks2>();
          if (t1.matchMCHTrackId() == t2.matchMCHTrackId())
            continue;
          if (t1.matchMFTTrackId() == t2.matchMFTTrackId())
            continue;
          VarManager::FillPairME<TPairType>(t1, t2);
          if constexpr ((TEventFillMap & VarManager::ObjTypes::ReducedEventQvector) > 0) {
            VarManager::FillPairVn<TPairType>(t1, t2);
          }
          pairSign = t1.sign() + t2.sign();
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
                         soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCovWithAmbiguities const& muons)
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

  void processBarrelOnlySkimmedNoCov(MyEventsSelected const& events,
                                     soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                                     MyBarrelTracksWithAmbiguities const& barrelTracks)
  {
    runSameEventPairing<false, VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMap>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processBarrelOnlyWithCollSkimmed(MyEventsVtxCovSelected const& events,
                                        soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                                        MyBarrelTracksWithCovWithAmbiguitiesWithColl const& barrelTracks)
  {
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCovWithColl>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processMuonOnlySkimmed(MyEventsVtxCovSelected const& events,
                              soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCovWithAmbiguities const& muons)
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

  void processMixingBarrelSkimmed(soa::Filtered<MyEventsHashSelected>& events,
                                  soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& trackAssocs, MyBarrelTracksWithCov const& tracks)
  {
    runSameSideMixing<pairTypeEE, gkEventFillMap>(events, trackAssocs, tracks, trackAssocsPerCollision);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processAllSkimmed, "Run all types of pairing, with skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlySkimmed, "Run barrel only pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlyWithCollSkimmed, "Run barrel only pairing, with skimmed tracks and with collision information", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlySkimmedNoCov, "Run barrel only pairing (no covariances), with skimmed tracks and with collision information", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMuonOnlySkimmed, "Run muon only pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMixingAllSkimmed, "Run all types of mixed pairing, with skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMixingBarrelSkimmed, "Run barrel type mixing pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", false);
};

// Run pairing for resonance with legs fulfilling separate cuts (asymmetric decay channel)
struct AnalysisAsymmetricPairing {

  Produces<aod::Ditracks> ditrackList;
  Produces<aod::DitracksExtra> ditrackExtraList;

  o2::base::MatLayerCylSet* fLUT = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  // Output objects
  OutputObj<THashList> fOutputList{"output"};

  // Configurables
  Configurable<string> fConfigLegCuts{"cfgLegCuts", "", "<leg-A-1>:<leg-B-1>[:<leg-C-1>],[<leg-A-2>:<leg-B-2>[:<leg-C-1>],...]"};
  Configurable<uint32_t> fConfigLegAFilterMask{"cfgLegAFilterMask", 0, "Filter mask corresponding to cuts in event-selection"};
  Configurable<uint32_t> fConfigLegBFilterMask{"cfgLegBFilterMask", 0, "Filter mask corresponding to cuts in event-selection"};
  Configurable<uint32_t> fConfigLegCFilterMask{"cfgLegCFilterMask", 0, "Filter mask corresponding to cuts in event-selection"};
  Configurable<string> fConfigCommonTrackCuts{"cfgCommonTrackCuts", "", "Comma separated list of cuts to be applied to all legs"};
  Configurable<string> fConfigPairCuts{"cfgPairCuts", "", "Comma separated list of pair cuts"};
  Configurable<bool> fConfigSkipAmbiguousIdCombinations{"cfgSkipAmbiguousIdCombinations", true, "Choose whether to skip pairs/triples which pass a stricter combination of cuts, e.g. KKPi triplets for D+ -> KPiPi"};

  Configurable<std::string> fConfigHistogramSubgroups{"cfgAsymmetricPairingHistogramsSubgroups", "barrel,vertexing", "Comma separated list of asymmetric-pairing histogram subgroups"};
  Configurable<bool> fConfigSameSignHistograms{"cfgSameSignHistograms", false, "Include same sign pair histograms for 2-prong decays"};
  Configurable<bool> fConfigAmbiguousHistograms{"cfgAmbiguousHistograms", false, "Include separate histograms for pairs/triplets with ambiguous tracks"};

  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> fConfigGRPMagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> fConfigUseRemoteField{"cfgUseRemoteField", false, "Choose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};

  Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
  Configurable<bool> fConfigUseAbsDCA{"cfgUseAbsDCA", false, "Use absolute DCA minimization instead of chi^2 minimization in secondary vertexing"};
  Configurable<bool> fConfigPropToPCA{"cfgPropToPCA", false, "Propagate tracks to secondary vertex"};
  Configurable<std::string> fConfigLutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  HistogramManager* fHistMan;

  std::map<int, std::vector<TString>> fTrackHistNames;
  std::vector<AnalysisCompositeCut> fPairCuts;

  // Filter masks to find legs in BarrelTrackCuts table
  uint32_t fLegAFilterMask;
  uint32_t fLegBFilterMask;
  uint32_t fLegCFilterMask;
  // Map tracking which pair of leg cuts the track cuts participate in
  std::map<int, uint32_t> fTrackCutFilterMasks;
  // Filter map for common track cuts
  uint32_t fCommonTrackCutMask;
  // Map tracking which common track cut the track cuts correspond to
  std::map<int, uint32_t> fCommonTrackCutFilterMasks;

  int fNLegCuts;
  int fNPairCuts;
  int fNCommonTrackCuts;

  Preslice<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;

  // Partitions for triplets and asymmetric pairs
  Partition<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> legACandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & fConfigLegAFilterMask) > uint32_t(0);
  Partition<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> legBCandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & fConfigLegBFilterMask) > uint32_t(0);
  Partition<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> legCCandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & fConfigLegCFilterMask) > uint32_t(0);

  void init(o2::framework::InitContext& context)
  {
    bool isDummy = context.mOptions.get<bool>("processDummy");
    if (isDummy) {
      LOG(info) << "Dummy function enabled. Skipping the rest of init()" << endl;
      return;
    }

    TString histNames = "";
    std::vector<TString> names;

    // Get the leg cut filter maps
    fLegAFilterMask = fConfigLegAFilterMask.value;
    fLegBFilterMask = fConfigLegBFilterMask.value;
    fLegCFilterMask = fConfigLegCFilterMask.value;
    // Get the pair cuts
    TString cutNamesStr = fConfigPairCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // Get the barrel track selection cuts
    string tempCuts;
    getTaskOptionValue<string>(context, "analysis-track-selection", "cfgTrackCuts", tempCuts, false);
    TString tempCutsStr = tempCuts;
    std::unique_ptr<TObjArray> objArray(tempCutsStr.Tokenize(","));
    // Get the common leg cuts
    int commonCutIdx;
    TString commonNamesStr = fConfigCommonTrackCuts.value;
    if (!commonNamesStr.IsNull()) { // if common track cuts
      std::unique_ptr<TObjArray> objArrayCommon(commonNamesStr.Tokenize(","));
      fNCommonTrackCuts = objArrayCommon->GetEntries();
      for (int icut = 0; icut < fNCommonTrackCuts; ++icut) {
        commonCutIdx = objArray->IndexOf(objArrayCommon->At(icut));
        if (commonCutIdx >= 0) {
          fCommonTrackCutMask |= uint32_t(1) << objArray->IndexOf(objArrayCommon->At(icut));
          fCommonTrackCutFilterMasks[icut] = uint32_t(1) << objArray->IndexOf(objArrayCommon->At(icut));
        } else {
          LOGF(fatal, "Common track cut %s was not calculated upstream. Check the config!", objArrayCommon->At(icut)->GetName());
        }
      }
    }
    // Check that the leg cut masks make sense
    if (static_cast<int>(std::floor(TMath::Log2(fLegAFilterMask))) + 1 > objArray->GetEntries()) {
      LOGF(fatal, "fConfigLegAFilterMask has highest bit at position %d, but track-selection only has %d cuts!", static_cast<int>(std::floor(TMath::Log2(fLegAFilterMask))) + 1, objArray->GetEntries());
    }
    if (static_cast<int>(std::floor(TMath::Log2(fLegBFilterMask))) + 1 > objArray->GetEntries()) {
      LOGF(fatal, "fConfigLegBFilterMask has highest bit at position %d, but track-selection only has %d cuts!", static_cast<int>(std::floor(TMath::Log2(fLegBFilterMask))) + 1, objArray->GetEntries());
    }
    if (static_cast<int>(std::floor(TMath::Log2(fLegCFilterMask))) + 1 > objArray->GetEntries()) {
      LOGF(fatal, "fConfigLegCFilterMask has highest bit at position %d, but track-selection only has %d cuts!", static_cast<int>(std::floor(TMath::Log2(fLegCFilterMask))) + 1, objArray->GetEntries());
    }

    // Get the cuts defining the legs
    uint32_t fConstructedLegAFilterMask = 0;
    uint32_t fConstructedLegBFilterMask = 0;
    uint32_t fConstructedLegCFilterMask = 0;
    TString legCutsStr = fConfigLegCuts.value;
    std::unique_ptr<TObjArray> objArrayLegs(legCutsStr.Tokenize(","));
    if (objArrayLegs->GetEntries() == 0) {
      LOG(fatal) << "No cuts defining legs. Check the config!";
    }
    fNLegCuts = objArrayLegs->GetEntries();
    std::vector<bool> isThreeProng;
    int legAIdx;
    int legBIdx;
    int legCIdx;
    // Loop over leg defining cuts
    for (int icut = 0; icut < fNLegCuts; ++icut) {
      TString legsStr = objArrayLegs->At(icut)->GetName();
      std::unique_ptr<TObjArray> legs(legsStr.Tokenize(":"));
      if (legs->GetEntries() == 3) {
        isThreeProng.push_back(true);
      } else if (legs->GetEntries() == 2) {
        isThreeProng.push_back(false);
      } else {
        LOGF(fatal, "Leg cuts %s has the wrong format and could not be parsed!", legsStr.Data());
        continue;
      }
      // Find leg cuts in the track selection cuts
      legAIdx = objArray->IndexOf(legs->At(0));
      if (legAIdx >= 0) {
        fConstructedLegAFilterMask |= (uint32_t(1) << legAIdx);
        fTrackCutFilterMasks[icut] |= uint32_t(1) << legAIdx;
      } else {
        LOGF(fatal, "Leg A cut %s was not calculated upstream. Check the config!", legs->At(0)->GetName());
        continue;
      }
      legBIdx = objArray->IndexOf(legs->At(1));
      if (legBIdx >= 0) {
        fConstructedLegBFilterMask |= (uint32_t(1) << legBIdx);
        fTrackCutFilterMasks[icut] |= uint32_t(1) << legBIdx;
      } else {
        LOGF(fatal, "Leg B cut %s was not calculated upstream. Check the config!", legs->At(1)->GetName());
        continue;
      }
      if (isThreeProng[icut]) {
        legCIdx = objArray->IndexOf(legs->At(2));
        if (legCIdx >= 0) {
          fConstructedLegCFilterMask |= (uint32_t(1) << legCIdx);
          fTrackCutFilterMasks[icut] |= uint32_t(1) << legCIdx;
        } else {
          LOGF(fatal, "Leg C cut %s was not calculated upstream. Check the config!", legs->At(2)->GetName());
          continue;
        }
      }
      if (isThreeProng[icut]) {
        names = {
          Form("TripletsBarrelSE_%s", legsStr.Data())};
        histNames += Form("%s;", names[0].Data());
        if (fConfigAmbiguousHistograms.value) {
          names.push_back(Form("TripletsBarrelSE_ambiguous_%s", legsStr.Data()));
          histNames += Form("%s;", names[1].Data());
        }
        fTrackHistNames[icut] = names;

        std::unique_ptr<TObjArray> objArrayCommon(commonNamesStr.Tokenize(","));
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
          names = {};
          names.push_back(Form("TripletsBarrelSE_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName()));
          histNames += Form("%s;", names[0].Data());
          fTrackHistNames[fNLegCuts + icut * fNCommonTrackCuts + iCommonCut] = names;
        }

        TString cutNamesStr = fConfigPairCuts.value;
        if (!cutNamesStr.IsNull()) { // if pair cuts
          std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
          fNPairCuts = objArrayPair->GetEntries();
          for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
            names = {};
            names.push_back(Form("TripletsBarrelSE_%s_%s", legsStr.Data(), objArrayPair->At(iPairCut)->GetName()));
            histNames += Form("%s;", names[0].Data());
            fTrackHistNames[fNLegCuts * (fNCommonTrackCuts + 1) + icut * fNPairCuts + iPairCut] = names;
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              names = {};
              names.push_back(Form("TripletsBarrelSE_%s_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName()));
              histNames += Form("%s;", names[0].Data());
              fTrackHistNames[(fNLegCuts * (fNCommonTrackCuts + 1) + fNLegCuts * fNPairCuts) + icut * (fNPairCuts * fNCommonTrackCuts + 1) + iCommonCut * (1 + fNPairCuts) + iPairCut] = names;
            } // end loop (common cuts)
          } // end loop (pair cuts)
        }   // end if (pair cuts)
      } else {
        names = {};
        std::vector<TString> pairHistPrefixes = {"PairsBarrelSEPM"};
        if (fConfigSameSignHistograms.value) {
          pairHistPrefixes.push_back("PairsBarrelSEPP");
          pairHistPrefixes.push_back("PairsBarrelSEMM");
        }
        int fNPairHistPrefixes = pairHistPrefixes.size();

        for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
          names.push_back(Form("%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data()));
          histNames += Form("%s;", names[iPrefix].Data());
        }
        if (fConfigAmbiguousHistograms.value) {
          for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
            names.push_back(Form("%s_ambiguous_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data()));
            histNames += Form("%s;", names[fNPairHistPrefixes + iPrefix].Data());
          }
        }
        fTrackHistNames[icut] = names;

        std::unique_ptr<TObjArray> objArrayCommon(commonNamesStr.Tokenize(","));
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
          names = {};
          for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
            names.push_back(Form("%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName()));
            histNames += Form("%s;", names[iPrefix].Data());
          }
          fTrackHistNames[fNLegCuts + icut * fNCommonTrackCuts + iCommonCut] = names;
        }

        if (!cutNamesStr.IsNull()) { // if pair cuts
          std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
          fNPairCuts = objArrayPair->GetEntries();
          for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
            names = {};
            for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
              names.push_back(Form("%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayPair->At(iPairCut)->GetName()));
              histNames += Form("%s;", names[iPrefix].Data());
            }
            fTrackHistNames[fNLegCuts * (fNCommonTrackCuts + 1) + icut * fNPairCuts + iPairCut] = names;
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              names = {};
              for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                names.push_back(Form("%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName()));
                histNames += Form("%s;", names[iPrefix].Data());
              }
              fTrackHistNames[(fNLegCuts * (fNCommonTrackCuts + 1) + fNLegCuts * fNPairCuts) + icut * (fNPairCuts * fNCommonTrackCuts + 1) + iCommonCut * (1 + fNPairCuts) + iPairCut] = names;
            } // end loop (common cuts)
          } // end loop (pair cuts)
        }   // end if (pair cuts)
      }
    }
    // Make sure the leg cuts are covered by the configured filter masks
    if (fLegAFilterMask != fConstructedLegAFilterMask) {
      LOGF(fatal, "cfgLegAFilterMask (%d) is not equal to the mask constructed by the cuts specified in cfgLegCuts (%d)!", fLegAFilterMask, fConstructedLegAFilterMask);
    }
    if (fLegBFilterMask != fConstructedLegBFilterMask) {
      LOGF(fatal, "cfgLegBFilterMask (%d) is not equal to the mask constructed by the cuts specified in cfgLegCuts (%d)!", fLegBFilterMask, fConstructedLegBFilterMask);
    }
    if (fLegCFilterMask != fConstructedLegCFilterMask) {
      LOGF(fatal, "cfgLegCFilterMask (%d) is not equal to the mask constructed by the cuts specified in cfgLegCuts (%d)!", fLegCFilterMask, fConstructedLegCFilterMask);
    }
    // Make sure only pairs or only triplets of leg cuts were given
    int tripletCheckSum = std::count(isThreeProng.begin(), isThreeProng.end(), true);
    if (tripletCheckSum != 0 && tripletCheckSum != fNLegCuts) {
      LOGF(fatal, "A mix of pairs and triplets was given as leg cuts. Check your config!");
    }

    fCurrentRun = 0;

    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();

    fLUT = o2::base::MatLayerCylSet::rectifyPtrFromFile(fCCDB->get<o2::base::MatLayerCylSet>(fConfigLutPath));
    VarManager::SetupMatLUTFwdDCAFitter(fLUT);

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    DefineHistograms(fHistMan, histNames.Data(), fConfigHistogramSubgroups.value.data()); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());                                      // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  void initParamsFromCCDB(uint64_t timestamp, bool isTriplets)
  {
    if (fConfigUseRemoteField.value) {
      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fConfigGRPMagPath, timestamp);
      float magField = 0.0;
      if (grpmag != nullptr) {
        magField = grpmag->getNominalL3Field();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", timestamp);
      }
      if (isTriplets) {
        if (fConfigUseKFVertexing.value) {
          VarManager::SetupThreeProngKFParticle(magField);
        } else {
          VarManager::SetupThreeProngDCAFitter(magField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value);
        }
      } else {
        if (fConfigUseKFVertexing.value) {
          VarManager::SetupTwoProngKFParticle(magField);
        } else {
          VarManager::SetupTwoProngDCAFitter(magField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value); // TODO: get these parameters from Configurables
        }
      }
    } else {
      if (isTriplets) {
        if (fConfigUseKFVertexing.value) {
          VarManager::SetupThreeProngKFParticle(fConfigMagField.value);
        } else {
          VarManager::SetupThreeProngDCAFitter(fConfigMagField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value);
        }
      } else {
        if (fConfigUseKFVertexing.value) {
          VarManager::SetupTwoProngKFParticle(fConfigMagField.value);
        } else {
          VarManager::SetupTwoProngDCAFitter(fConfigMagField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value); // TODO: get these parameters from Configurables
        }
      }
    }
  }

  // Template function to run same event pairing with asymmetric pairs (e.g. kaon-pion)
  template <bool TTwoProngFitter, int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTrackAssocs, typename TTracks>
  void runAsymmetricPairing(TEvents const& events, Preslice<TTrackAssocs>& preslice, TTrackAssocs const& /*assocs*/, TTracks const& /*tracks*/)
  {
    if (events.size() > 0) { // Additional protection to avoid crashing of events.begin().runNumber()
      if (fCurrentRun != events.begin().runNumber()) {
        initParamsFromCCDB(events.begin().timestamp(), false);
        fCurrentRun = events.begin().runNumber();
      }
    }

    std::map<int, std::vector<TString>> histNames = fTrackHistNames;

    int sign1 = 0;
    int sign2 = 0;
    ditrackList.reserve(1);
    ditrackExtraList.reserve(1);

    constexpr bool trackHasCov = ((TTrackFillMap & VarManager::ObjTypes::TrackCov) > 0 || (TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelCov) > 0);

    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event, VarManager::fgValues);

      auto groupedLegAAssocs = legACandidateAssocs.sliceBy(preslice, event.globalIndex());
      if (groupedLegAAssocs.size() == 0) {
        continue;
      }
      auto groupedLegBAssocs = legBCandidateAssocs.sliceBy(preslice, event.globalIndex());
      if (groupedLegBAssocs.size() == 0) {
        continue;
      }

      // TODO: Think about double counting
      std::set<std::pair<int, int>> globIdxPairs;
      for (auto& [a1, a2] : combinations(soa::CombinationsFullIndexPolicy(groupedLegAAssocs, groupedLegBAssocs))) {

        uint32_t twoTrackFilter = 0;
        uint32_t twoTrackCommonFilter = 0;
        uint32_t pairFilter = 0;
        bool isPairIdWrong = false;
        for (int icut = 0; icut < fNLegCuts; ++icut) {
          // Find leg pair definitions both candidates participate in
          if ((((a1.isBarrelSelected_raw() & fLegAFilterMask) | (a2.isBarrelSelected_raw() & fLegBFilterMask)) & fTrackCutFilterMasks[icut]) == fTrackCutFilterMasks[icut]) {
            twoTrackFilter |= (uint32_t(1) << icut);
            // If the supposed pion passes a kaon cut, this is a K+K-. Skip it.
            if (TPairType == VarManager::kDecayToKPi && fConfigSkipAmbiguousIdCombinations.value) {
              if (a2.isBarrelSelected_raw() & fLegAFilterMask) {
                isPairIdWrong = true;
              }
            }
          }
        }

        if (!twoTrackFilter || isPairIdWrong) {
          continue;
        }

        // Find common track cuts both candidates pass
        twoTrackCommonFilter |= a1.isBarrelSelected_raw() & a2.isBarrelSelected_raw() & fCommonTrackCutMask;

        auto t1 = a1.template reducedtrack_as<TTracks>();
        auto t2 = a2.template reducedtrack_as<TTracks>();

        // Avoid self-pairs
        if (t1.globalIndex() == t2.globalIndex()) {
          continue;
        }

        sign1 = t1.sign();
        sign2 = t2.sign();
        // store the ambiguity number of the two dilepton legs in the last 4 digits of the two-track filter
        if (t1.barrelAmbiguityInBunch() > 1 || t1.barrelAmbiguityOutOfBunch() > 1) {
          twoTrackFilter |= (uint32_t(1) << 30);
        }
        if (t2.barrelAmbiguityInBunch() > 1 || t2.barrelAmbiguityOutOfBunch() > 1) {
          twoTrackFilter |= (uint32_t(1) << 31);
        }

        VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
        if constexpr (TTwoProngFitter) {
          VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, fConfigPropToPCA);
        }

        // Fill histograms
        bool isAmbi = false;
        for (int icut = 0; icut < fNLegCuts; icut++) {
          if (twoTrackFilter & (uint32_t(1) << icut)) {
            isAmbi = (twoTrackFilter & (uint32_t(1) << 30)) || (twoTrackFilter & (uint32_t(1) << 31));
            if (sign1 * sign2 < 0) {
              fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
              if (isAmbi && fConfigAmbiguousHistograms.value) {
                fHistMan->FillHistClass(histNames[icut][3].Data(), VarManager::fgValues);
              }
            } else if (fConfigSameSignHistograms.value) {
              if (sign1 > 0) {
                fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
                if (isAmbi && fConfigAmbiguousHistograms.value) {
                  fHistMan->FillHistClass(histNames[icut][4].Data(), VarManager::fgValues);
                }
              } else {
                fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
                if (isAmbi && fConfigAmbiguousHistograms.value) {
                  fHistMan->FillHistClass(histNames[icut][5].Data(), VarManager::fgValues);
                }
              }
            }
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
              if (twoTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
                if (sign1 * sign2 < 0) {
                  fHistMan->FillHistClass(histNames[fNLegCuts + icut * fNCommonTrackCuts + iCommonCut][0].Data(), VarManager::fgValues);
                } else if (fConfigSameSignHistograms.value) {
                  if (sign1 > 0) {
                    fHistMan->FillHistClass(histNames[fNLegCuts + icut * fNCommonTrackCuts + iCommonCut][1].Data(), VarManager::fgValues);
                  } else {
                    fHistMan->FillHistClass(histNames[fNLegCuts + icut * fNCommonTrackCuts + iCommonCut][2].Data(), VarManager::fgValues);
                  }
                }
              }
            } // end loop (common cuts)
            for (unsigned int iPairCut = 0; iPairCut < fPairCuts.size(); iPairCut++) {
              AnalysisCompositeCut cut = fPairCuts.at(iPairCut);
              if (!(cut.IsSelected(VarManager::fgValues))) // apply pair cuts
                continue;
              pairFilter |= (uint32_t(1) << iPairCut);
              // Histograms with pair cuts
              if (sign1 * sign2 < 0) {
                fHistMan->FillHistClass(histNames[fNLegCuts * (fNCommonTrackCuts + 1) + icut * fNPairCuts + iPairCut][0].Data(), VarManager::fgValues);
              } else if (fConfigSameSignHistograms.value) {
                if (sign1 > 0) {
                  fHistMan->FillHistClass(histNames[fNLegCuts * (fNCommonTrackCuts + 1) + icut * fNPairCuts + iPairCut][1].Data(), VarManager::fgValues);
                } else {
                  fHistMan->FillHistClass(histNames[fNLegCuts * (fNCommonTrackCuts + 1) + icut * fNPairCuts + iPairCut][2].Data(), VarManager::fgValues);
                }
              }
              // Histograms with pair cuts and common track cuts
              for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                if (twoTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
                  if (sign1 * sign2 < 0) {
                    fHistMan->FillHistClass(histNames[(fNLegCuts * (fNCommonTrackCuts + 1) + fNLegCuts * fNPairCuts) + icut * (fNPairCuts * fNCommonTrackCuts + 1) + iCommonCut * (1 + fNPairCuts) + iPairCut][0].Data(), VarManager::fgValues);
                  } else if (fConfigSameSignHistograms.value) {
                    if (sign1 > 0) {
                      fHistMan->FillHistClass(histNames[(fNLegCuts * (fNCommonTrackCuts + 1) + fNLegCuts * fNPairCuts) + icut * (fNPairCuts * fNCommonTrackCuts + 1) + iCommonCut * (1 + fNPairCuts) + iPairCut][1].Data(), VarManager::fgValues);
                    } else {
                      fHistMan->FillHistClass(histNames[(fNLegCuts * (fNCommonTrackCuts + 1) + fNLegCuts * fNPairCuts) + icut * (fNPairCuts * fNCommonTrackCuts + 1) + iCommonCut * (1 + fNPairCuts) + iPairCut][2].Data(), VarManager::fgValues);
                    }
                  }
                }
              }
            } // end loop (pair cuts)
          }
        } // end loop (cuts)
        ditrackList(event.globalIndex(), VarManager::fgValues[VarManager::kMass],
                    VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi],
                    t1.sign() + t2.sign(), twoTrackFilter, pairFilter, twoTrackCommonFilter);
        if constexpr (trackHasCov && TTwoProngFitter) {
          ditrackExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauzProjected], VarManager::fgValues[VarManager::kVertexingLzProjected], VarManager::fgValues[VarManager::kVertexingLxyProjected]);
        }
      }   // end inner assoc loop (leg A)
    }     // end event loop
  }

  // Template function to run same event triplets (e.g. D+->K-pi+pi+)
  template <bool TThreeProngFitter, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTrackAssocs, typename TTracks>
  void runThreeProng(TEvents const& events, Preslice<TTrackAssocs>& preslice, TTrackAssocs const& /*assocs*/, TTracks const& tracks, VarManager::PairCandidateType tripletType)
  {
    if (events.size() > 0) { // Additional protection to avoid crashing of events.begin().runNumber()
      if (fCurrentRun != events.begin().runNumber()) {
        initParamsFromCCDB(events.begin().timestamp(), true);
        fCurrentRun = events.begin().runNumber();
      }
    }

    std::map<int, std::vector<TString>> histNames = fTrackHistNames;

    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event, VarManager::fgValues);

      auto groupedLegAAssocs = legACandidateAssocs.sliceBy(preslice, event.globalIndex());
      if (groupedLegAAssocs.size() == 0) {
        continue;
      }
      auto groupedLegBAssocs = legBCandidateAssocs.sliceBy(preslice, event.globalIndex());
      if (groupedLegBAssocs.size() == 0) {
        continue;
      }
      auto groupedLegCAssocs = legCCandidateAssocs.sliceBy(preslice, event.globalIndex());
      if (groupedLegCAssocs.size() == 0) {
        continue;
      }

      std::set<std::tuple<int64_t, int64_t, int64_t>> globIdxTriplets;
      // Based on triplet type, make suitable combinations of the partitions
      if (tripletType == VarManager::kTripleCandidateToPKPi) {
        for (auto& [a1, a2, a3] : combinations(soa::CombinationsFullIndexPolicy(groupedLegAAssocs, groupedLegBAssocs, groupedLegCAssocs))) {
          readTriplet<TThreeProngFitter, TEventFillMap, TTrackFillMap>(a1, a2, a3, tracks, event, tripletType, histNames);
        }
      } else if (tripletType == VarManager::kTripleCandidateToKPiPi) {
        for (auto& a1 : groupedLegAAssocs) {
          for (auto& [a2, a3] : combinations(groupedLegBAssocs, groupedLegCAssocs)) {
            readTriplet<TThreeProngFitter, TEventFillMap, TTrackFillMap>(a1, a2, a3, tracks, event, tripletType, histNames);
          }
        }
      } else {
        LOG(fatal) << "Given tripletType not recognized. Don't know how to make combinations!" << endl;
      }
    } // end event loop
  }

  // Helper function to process triplet
  template <bool TThreeProngFitter, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TTrackAssoc, typename TTracks, typename TEvent>
  void readTriplet(TTrackAssoc const& a1, TTrackAssoc const& a2, TTrackAssoc const& a3, TTracks const& /*tracks*/, TEvent const& event, VarManager::PairCandidateType tripletType, std::map<int, std::vector<TString>> histNames)
  {
    uint32_t threeTrackFilter = 0;
    uint32_t threeTrackCommonFilter = 0;
    for (int icut = 0; icut < fNLegCuts; ++icut) {
      // Find out which leg cut combination the triplet passes
      if ((((a1.isBarrelSelected_raw() & fLegAFilterMask) | (a2.isBarrelSelected_raw() & fLegBFilterMask) | (a3.isBarrelSelected_raw() & fLegCFilterMask)) & fTrackCutFilterMasks[icut]) == fTrackCutFilterMasks[icut]) {
        threeTrackFilter |= (uint32_t(1) << icut);
        if (tripletType == VarManager::kTripleCandidateToPKPi && fConfigSkipAmbiguousIdCombinations.value) {
          // Check if the supposed pion passes as a proton or kaon, if so, skip this triplet. It is pKp or pKK.
          if ((a3.isBarrelSelected_raw() & fLegAFilterMask) || (a3.isBarrelSelected_raw() & fLegBFilterMask)) {
            return;
          }
          // Check if the supposed kaon passes as a proton, if so, skip this triplet. It is ppPi.
          if (a2.isBarrelSelected_raw() & fLegAFilterMask) {
            return;
          }
        }
        if (tripletType == VarManager::kTripleCandidateToKPiPi && fConfigSkipAmbiguousIdCombinations.value) {
          // Check if one of the supposed pions pass as a kaon, if so, skip this triplet. It is KKPi.
          if ((a2.isBarrelSelected_raw() & fLegAFilterMask) || (a3.isBarrelSelected_raw() & fLegAFilterMask)) {
            return;
          }
        }
      }
    }
    if (!threeTrackFilter) {
      return;
    }

    // Find common track cuts all candidates pass
    threeTrackCommonFilter |= a1.isBarrelSelected_raw() & a2.isBarrelSelected_raw() & a3.isBarrelSelected_raw() & fCommonTrackCutMask;

    auto t1 = a1.template reducedtrack_as<TTracks>();
    auto t2 = a2.template reducedtrack_as<TTracks>();
    auto t3 = a3.template reducedtrack_as<TTracks>();

    // Avoid self-pairs
    if (t1 == t2 || t1 == t3 || t2 == t3) {
      return;
    }

    // store the ambiguity of the three legs in the last 3 digits of the two-track filter
    if (t1.barrelAmbiguityInBunch() > 1 || t1.barrelAmbiguityOutOfBunch() > 1) {
      threeTrackFilter |= (uint32_t(1) << 29);
    }
    if (t2.barrelAmbiguityInBunch() > 1 || t2.barrelAmbiguityOutOfBunch() > 1) {
      threeTrackFilter |= (uint32_t(1) << 30);
    }
    if (t3.barrelAmbiguityInBunch() > 1 || t3.barrelAmbiguityOutOfBunch() > 1) {
      threeTrackFilter |= (uint32_t(1) << 31);
    }

    VarManager::FillTriple(t1, t2, t3, VarManager::fgValues, tripletType);
    if constexpr (TThreeProngFitter) {
      VarManager::FillTripletVertexing<TEventFillMap, TTrackFillMap>(event, t1, t2, t3, tripletType);
    }

    // Fill histograms
    for (int icut = 0; icut < fNLegCuts; icut++) {
      if (threeTrackFilter & (uint32_t(1) << icut)) {
        fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
        if (fConfigAmbiguousHistograms.value && ((threeTrackFilter & (uint32_t(1) << 29)) || (threeTrackFilter & (uint32_t(1) << 30)) || (threeTrackFilter & (uint32_t(1) << 31)))) {
          fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
        }
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
          if (threeTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
            fHistMan->FillHistClass(histNames[fNLegCuts + icut * fNCommonTrackCuts + iCommonCut][0].Data(), VarManager::fgValues);
          }
        } // end loop (common cuts)
        for (unsigned int iPairCut = 0; iPairCut < fPairCuts.size(); iPairCut++) {
          AnalysisCompositeCut cut = fPairCuts.at(iPairCut);
          if (!(cut.IsSelected(VarManager::fgValues))) { // apply pair cuts
            continue;
          }
          // Histograms with pair cuts
          fHistMan->FillHistClass(histNames[fNLegCuts * (fNCommonTrackCuts + 1) + icut * fNPairCuts + iPairCut][0].Data(), VarManager::fgValues);
          // Histograms with pair cuts and common track cuts
          for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
            if (threeTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
              fHistMan->FillHistClass(histNames[(fNLegCuts * (fNCommonTrackCuts + 1) + fNLegCuts * fNPairCuts) + icut * (fNPairCuts * fNCommonTrackCuts + 1) + iCommonCut * (1 + fNPairCuts) + iPairCut][0].Data(), VarManager::fgValues);
            }
          }
        } // end loop (pair cuts)
      }
    } // end loop (cuts)
  }

  void processKaonPionSkimmed(MyEventsVtxCovZdcSelected const& events,
                              soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& barrelAssocs,
                              MyBarrelTracksWithCovWithAmbiguities const& barrelTracks)
  {
    runAsymmetricPairing<true, VarManager::kDecayToKPi, gkEventFillMapWithCovZdc, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processKaonPionPionSkimmed(MyEventsVtxCovZdcSelected const& events,
                                  soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& barrelAssocs,
                                  MyBarrelTracksWithCovWithAmbiguities const& barrelTracks)
  {
    runThreeProng<true, gkEventFillMapWithCovZdc, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, VarManager::kTripleCandidateToKPiPi);
  }

  void processProtonKaonPionSkimmed(MyEventsVtxCovZdcSelected const& events,
                                    soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& barrelAssocs,
                                    MyBarrelTracksWithCovWithAmbiguities const& barrelTracks)
  {
    runThreeProng<true, gkEventFillMapWithCovZdc, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, VarManager::kTripleCandidateToPKPi);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisAsymmetricPairing, processKaonPionSkimmed, "Run kaon pion pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processKaonPionPionSkimmed, "Run kaon pion pion triplets, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processProtonKaonPionSkimmed, "Run proton kaon pion triplets, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", true);
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
  int fNPairCuts;
  int fNCommonTrackCuts;
  std::map<int, int> fCommonTrackCutMap;
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
    bool isBarrelAsymmetric = context.mOptions.get<bool>("processDstarToD0Pi");
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
    TString histNames = "";

    // For each track/muon selection used to produce dileptons, create a separate histogram directory using the
    // name of the track/muon cut.
    // Also, create a map which will hold the name of the histogram directories so they can be accessed directly in the pairing loop
    if (isBarrel || isMuon || isBarrelAsymmetric) {
      // get the list of single track and muon cuts computed in the dedicated tasks upstream
      string tempCutsSingle;
      if (isBarrel || isBarrelAsymmetric) {
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
      string tempCutsSinglePair;
      string pairCuts;
      string pairCommonCuts;
      string tempCutsTrack;
      if (isBarrel) {
        getTaskOptionValue<string>(context, "analysis-same-event-pairing", "cfgTrackCuts", tempCutsSinglePair, false);
        getTaskOptionValue<string>(context, "analysis-same-event-pairing", "cfgPairCuts", pairCuts, false);
      } else if (isMuon) {
        getTaskOptionValue<string>(context, "analysis-same-event-pairing", "cfgMuonCuts", tempCutsSinglePair, false);
        getTaskOptionValue<string>(context, "analysis-same-event-pairing", "cfgPairCuts", pairCuts, false);
      } else if (isBarrelAsymmetric) {
        getTaskOptionValue<string>(context, "analysis-track-selection", "cfgTrackCuts", tempCutsTrack, false);
        getTaskOptionValue<string>(context, "analysis-asymmetric-pairing", "cfgLegCuts", tempCutsSinglePair, false);
        getTaskOptionValue<string>(context, "analysis-asymmetric-pairing", "cfgPairCuts", pairCuts, false);
        getTaskOptionValue<string>(context, "analysis-asymmetric-pairing", "cfgCommonTrackCuts", pairCommonCuts, false);
      }

      // If asymmetric pair is used, it may have common track cuts
      TString pairCommonCutsStr = pairCommonCuts;
      if (!pairCommonCutsStr.IsNull()) { // if common track cuts
        TString tempCutsTrackStr = tempCutsTrack;
        std::unique_ptr<TObjArray> objArrayTempTrack(tempCutsTrackStr.Tokenize(","));
        int fNTempTrackCuts = objArrayTempTrack->GetEntries();
        std::unique_ptr<TObjArray> objArrayCommon(pairCommonCutsStr.Tokenize(","));
        fNCommonTrackCuts = objArrayCommon->GetEntries();
        for (int icut = 0; icut < fNCommonTrackCuts; ++icut) {
          for (int iicut = 0; iicut < fNTempTrackCuts; ++iicut) {
            if (std::strcmp(objArrayCommon->At(icut)->GetName(), objArrayTempTrack->At(iicut)->GetName()) == 0) {
              fCommonTrackCutMap[icut] = iicut;
            }
          }
        }
      }

      TString tempCutsSinglePairStr = tempCutsSinglePair;
      bool cutFound;
      if (!tempCutsSingleStr.IsNull() && !tempCutsSinglePairStr.IsNull()) {
        std::unique_ptr<TObjArray> objArray(tempCutsSinglePairStr.Tokenize(","));
        fNCuts = objArray->GetEntries();
        for (int icut = 0; icut < fNCuts; ++icut) {
          TString tempStr = objArray->At(icut)->GetName();
          if (!isBarrelAsymmetric) {
            cutFound = objArraySingleCuts->FindObject(tempStr.Data()) != nullptr;
          } else {
            std::unique_ptr<TObjArray> legObjArray(tempStr.Tokenize(":"));
            cutFound = true;
            for (int iicut = 0; iicut < legObjArray->GetEntries(); ++iicut) {
              TString tempLegStr = legObjArray->At(iicut)->GetName();
              if (objArraySingleCuts->FindObject(tempLegStr.Data()) == nullptr) {
                cutFound = false;
              }
            }
          }
          if (cutFound) {
            fHistNamesDileptonTrack[icut] = Form("DileptonTrack_%s_%s", tempStr.Data(), fConfigTrackCut.value.data());
            fHistNamesDileptons[icut] = Form("DileptonsSelected_%s", tempStr.Data());
            TString pairCutsStr = pairCuts;
            DefineHistograms(fHistMan, fHistNamesDileptonTrack[icut], fConfigHistogramSubgroups.value.data()); // define dilepton-track histograms
            DefineHistograms(fHistMan, fHistNamesDileptons[icut], "barrel,vertexing");                         // define dilepton histograms
            if (!pairCommonCutsStr.IsNull()) {
              std::unique_ptr<TObjArray> objArrayCommon(pairCommonCutsStr.Tokenize(","));
              for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                // fTrackHistNames[fNLegCuts + icut * fNCommonTrackCuts + iCommonCut] = names;
                fHistNamesDileptonTrack[fNCuts + icut * fNCommonTrackCuts + iCommonCut] = Form("DileptonTrack_%s_%s_%s", tempStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), fConfigTrackCut.value.data());
                fHistNamesDileptons[fNCuts + icut * fNCommonTrackCuts + iCommonCut] = Form("DileptonsSelected_%s_%s", tempStr.Data(), objArrayCommon->At(iCommonCut)->GetName());
                DefineHistograms(fHistMan, fHistNamesDileptonTrack[fNCuts + icut * fNCommonTrackCuts + iCommonCut], fConfigHistogramSubgroups.value.data()); // define dilepton-track histograms
                DefineHistograms(fHistMan, fHistNamesDileptons[fNCuts + icut * fNCommonTrackCuts + iCommonCut], "barrel,vertexing");                         // define dilepton histograms
              }
            }
            if (!pairCutsStr.IsNull()) {
              std::unique_ptr<TObjArray> objArrayPairCuts(pairCutsStr.Tokenize(","));
              fNPairCuts = objArrayPairCuts->GetEntries();
              for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) {
                fHistNamesDileptonTrack[fNCuts * (fNCommonTrackCuts + 1) + icut * fNPairCuts + iPairCut] = Form("DileptonTrack_%s_%s_%s", tempStr.Data(), objArrayPairCuts->At(iPairCut)->GetName(), fConfigTrackCut.value.data());
                fHistNamesDileptons[fNCuts * (fNCommonTrackCuts + 1) + icut * fNPairCuts + iPairCut] = Form("DileptonsSelected_%s_%s", tempStr.Data(), objArrayPairCuts->At(iPairCut)->GetName());
                DefineHistograms(fHistMan, fHistNamesDileptonTrack[fNCuts * (fNCommonTrackCuts + 1) + icut * fNPairCuts + iPairCut], fConfigHistogramSubgroups.value.data()); // define dilepton-track histograms
                DefineHistograms(fHistMan, fHistNamesDileptons[fNCuts * (fNCommonTrackCuts + 1) + icut * fNPairCuts + iPairCut], "barrel,vertexing");                         // define dilepton histograms
                if (!pairCommonCutsStr.IsNull()) {
                  std::unique_ptr<TObjArray> objArrayCommon(pairCommonCutsStr.Tokenize(","));
                  for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                    fHistNamesDileptonTrack[(fNCuts * (fNCommonTrackCuts + 1) + fNCuts * fNPairCuts) + icut * (fNPairCuts * fNCommonTrackCuts + 1) + iCommonCut * (1 + fNPairCuts) + iPairCut] = Form("DileptonTrack_%s_%s_%s_%s", tempStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPairCuts->At(iPairCut)->GetName(), fConfigTrackCut.value.data());
                    fHistNamesDileptons[(fNCuts * (fNCommonTrackCuts + 1) + fNCuts * fNPairCuts) + icut * (fNPairCuts * fNCommonTrackCuts + 1) + iCommonCut * (1 + fNPairCuts) + iPairCut] = Form("DileptonsSelected_%s_%s_%s", tempStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPairCuts->At(iPairCut)->GetName());
                    DefineHistograms(fHistMan, fHistNamesDileptonTrack[(fNCuts * (fNCommonTrackCuts + 1) + fNCuts * fNPairCuts) + icut * (fNPairCuts * fNCommonTrackCuts + 1) + iCommonCut * (1 + fNPairCuts) + iPairCut], fConfigHistogramSubgroups.value.data()); // define dilepton-track histograms
                    DefineHistograms(fHistMan, fHistNamesDileptons[(fNCuts * (fNCommonTrackCuts + 1) + fNCuts * fNPairCuts) + icut * (fNPairCuts * fNCommonTrackCuts + 1) + iCommonCut * (1 + fNPairCuts) + iPairCut], "barrel,vertexing");                         // define dilepton histograms
                  }
                }
              }
            }
            if (isBarrelME || isMuonME) {
              fHistNamesME[icut] = Form("DileptonTrackME_%s", tempStr.Data());
              DefineHistograms(fHistMan, fHistNamesME[icut], "mixedevent"); // define ME histograms
            }
          }
        }
        for (int icut = 0; icut < objArraySingleCuts->GetEntries(); ++icut) {
          TString tempStr = objArraySingleCuts->At(icut)->GetName();
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
        VarManager::SetupThreeProngDCAFitter(magField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, false); // TODO: get these parameters from Configurables
      }
    } else {
      if (fConfigUseKFVertexing.value) {
        VarManager::SetupThreeProngKFParticle(fConfigMagField.value);
      } else {
        VarManager::SetupThreeProngDCAFitter(fConfigMagField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, false); // TODO: get these parameters from Configurables
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
          if constexpr (TCandidateType == VarManager::kDstarToD0KPiPi) { // Dielectrons and Dimuons don't have the PairFilterMap column
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
              if (dilepton.commonFilterMap_bit(fCommonTrackCutMap[iCommonCut])) {
                fHistMan->FillHistClass(fHistNamesDileptons[fNCuts + icut * fNCommonTrackCuts + iCommonCut].Data(), fValuesDilepton);
              }
            }
            for (int iPairCut = 0; iPairCut < fNPairCuts; iPairCut++) {
              if (dilepton.pairFilterMap_bit(iPairCut)) {
                fHistMan->FillHistClass(fHistNamesDileptons[fNCuts * (fNCommonTrackCuts + 1) + icut * fNPairCuts + iPairCut].Data(), fValuesDilepton);
                for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
                  if (dilepton.commonFilterMap_bit(fCommonTrackCutMap[iCommonCut])) {
                    fHistMan->FillHistClass(fHistNamesDileptons[(fNCuts * (fNCommonTrackCuts + 1) + fNCuts * fNPairCuts) + icut * (fNPairCuts * fNCommonTrackCuts + 1) + iCommonCut * (1 + fNPairCuts) + iPairCut].Data(), fValuesDilepton);
                  }
                }
              }
            }
          }
        }
      }

      // loop over hadrons
      for (auto& assoc : assocs) {
        if constexpr (TCandidateType == VarManager::kBtoJpsiEEK || TCandidateType == VarManager::kDstarToD0KPiPi) {
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
            if constexpr (TCandidateType == VarManager::kDstarToD0KPiPi) { // Dielectrons and Dimuons don't have the PairFilterMap column
              for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
                if (dilepton.commonFilterMap_bit(fCommonTrackCutMap[iCommonCut])) {
                  fHistMan->FillHistClass(fHistNamesDileptonTrack[fNCuts + icut * fNCommonTrackCuts + iCommonCut].Data(), fValuesHadron);
                }
              }
              for (int iPairCut = 0; iPairCut < fNPairCuts; iPairCut++) {
                if (dilepton.pairFilterMap_bit(iPairCut)) {
                  fHistMan->FillHistClass(fHistNamesDileptonTrack[fNCuts * (fNCommonTrackCuts + 1) + icut * fNPairCuts + iPairCut].Data(), fValuesHadron);
                  for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
                    if (dilepton.commonFilterMap_bit(fCommonTrackCutMap[iCommonCut])) {
                      fHistMan->FillHistClass(fHistNamesDileptonTrack[(fNCuts * (fNCommonTrackCuts + 1) + fNCuts * fNPairCuts) + icut * (fNPairCuts * fNCommonTrackCuts + 1) + iCommonCut * (1 + fNPairCuts) + iPairCut].Data(), fValuesHadron);
                    }
                  }
                }
              }
            }
          }
        }
        // table to be written out for ML analysis
        BmesonsTable(fValuesHadron[VarManager::kPairMass], fValuesHadron[VarManager::kPairPt], fValuesHadron[VarManager::kVertexingLxy], fValuesHadron[VarManager::kVertexingLxyz], fValuesHadron[VarManager::kVertexingLz], fValuesHadron[VarManager::kVertexingTauxy], fValuesHadron[VarManager::kVertexingTauz], fValuesHadron[VarManager::kCosPointingAngle], fValuesHadron[VarManager::kVertexingChi2PCA]);
      }
    }
  }

  Preslice<aod::ReducedTracksAssoc> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;
  Preslice<MyDielectronCandidates> dielectronsPerCollision = aod::reducedpair::reducedeventId;
  Preslice<MyDitrackCandidates> ditracksPerCollision = aod::reducedpair::reducedeventId;

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

  void processDstarToD0Pi(soa::Filtered<MyEventsVtxCovSelected> const& events,
                          soa::Filtered<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> const& assocs,
                          MyBarrelTracksWithCov const& tracks, soa::Filtered<MyDitrackCandidates> const& ditracks)
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
      auto groupedDitracks = ditracks.sliceBy(ditracksPerCollision, event.globalIndex());
      runDileptonHadron<VarManager::kDstarToD0KPiPi, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, groupedBarrelAssocs, tracks, groupedDitracks);
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
  PROCESS_SWITCH(AnalysisDileptonTrack, processDstarToD0Pi, "Run barrel pairing of D0 daughters with pion candidate, using skimmed data", false);
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
    adaptAnalysisTask<AnalysisAsymmetricPairing>(cfgc),
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

    if (classStr.Contains("Triplets")) {
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
