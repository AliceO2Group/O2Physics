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
//
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <memory>
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
// TODO: the barrel amd muon selection columns are bit maps so unsigned types should be used, however, for now this is not supported in Filter expressions
// TODO: For now in the tasks we just statically convert from unsigned int to int, which should be fine as long as we do
//      not use a large number of bits (>=30)
// Bcandidate columns for ML analysis of B->Jpsi+K
DECLARE_SOA_COLUMN(MixingHash, mixingHash, int);
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelected, isBarrelSelected, int);
DECLARE_SOA_COLUMN(IsMuonSelected, isMuonSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelectedPrefilter, isBarrelSelectedPrefilter, int);
DECLARE_SOA_COLUMN(IsPrefilterVetoed, isPrefilterVetoed, int);
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

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASH", dqanalysisflags::MixingHash);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTS", dqanalysisflags::IsBarrelSelected, dqanalysisflags::IsBarrelSelectedPrefilter);
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTER", dqanalysisflags::IsPrefilterVetoed);
DECLARE_SOA_TABLE(BmesonCandidates, "AOD", "DQBMESONS", dqanalysisflags::massBcandidate, dqanalysisflags::pTBcandidate, dqanalysisflags::LxyBcandidate, dqanalysisflags::LxyzBcandidate, dqanalysisflags::LzBcandidate, dqanalysisflags::TauxyBcandidate, dqanalysisflags::TauzBcandidate, dqanalysisflags::CosPBcandidate, dqanalysisflags::Chi2Bcandidate);
} // namespace o2::aod

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsHashSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts>;
using MyEventsVtxCovSelectedMultExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;
using MyEventsVtxCovSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvector>;
using MyEventsVtxCovQvectorExtraWithRefFlow = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::ReducedEventsQvector, aod::ReducedEventsQvectorExtra, aod::ReducedEventsRefFlow>;
using MyEventsVtxCovSelectedQvectorExtraWithRefFlow = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvector, aod::ReducedEventsQvectorExtra, aod::ReducedEventsRefFlow>;
using MyEventsVtxCovSelectedQvectorCentr = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvectorCentr>;
using MyEventsQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector>;
using MyEventsQvectorMultExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;
using MyEventsQvectorCentr = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvectorCentr>;
using MyEventsQvectorCentrMultExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvectorCentr, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;
using MyEventsQvectorExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector, aod::ReducedEventsQvectorExtra>;
using MyEventsHashSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvector>;
using MyEventsHashSelectedQvectorExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, ReducedEventVtxCov, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvector, aod::ReducedEventsQvectorExtra, aod::ReducedEventsRefFlow>;
using MyEventsHashSelectedQvectorCentr = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvectorCentr>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;
using MyBarrelTracksSelectedWithPrefilter = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::Prefilter>;
using MyBarrelTracksSelectedWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;
using MyBarrelTracksSelectedWithColl = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::ReducedTracksBarrelInfo>;
using MyDielectronCandidates = soa::Join<aod::Dielectrons, aod::DielectronsExtra>;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonTracksSelected = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::MuonTrackCuts>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyMuonTracksSelectedWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonTrackCuts>;
using MyMuonTracksSelectedWithColl = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsInfo, aod::MuonTrackCuts>;
using MyMftTracks = soa::Join<aod::ReducedMFTs, aod::ReducedMFTsExtra>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;
constexpr static uint32_t gkEventFillMapWithQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkEventFillMapWithMultExtra = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventMultExtra;
constexpr static uint32_t gkEventFillMapWithQvectorMultExtra = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector | VarManager::ObjTypes::ReducedEventMultExtra;
constexpr static uint32_t gkEventFillMapWithQvectorCentr = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::CollisionQvect;
constexpr static uint32_t gkEventFillMapWithQvectorCentrMultExtra = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::CollisionQvect | VarManager::ObjTypes::ReducedEventMultExtra;
constexpr static uint32_t gkEventFillMapWithCovQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkEventFillMapWithCovQvectorExtraWithRefFlow = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector | VarManager::ObjTypes::ReducedEventQvectorExtra | VarManager::ObjTypes::ReducedEventRefFlow;
constexpr static uint32_t gkEventFillMapWithCovQvectorCentr = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::CollisionQvect;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithColl = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID | VarManager::ObjTypes::ReducedTrackCollInfo;

constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCov;
constexpr static uint32_t gkMuonFillMapWithColl = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCollInfo;

constexpr static int pairTypeEE = VarManager::kDecayToEE;
constexpr static int pairTypeMuMu = VarManager::kDecayToMuMu;
constexpr static int pairTypeEMu = VarManager::kElectronMuon;

// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar); // defines histograms for all tasks

struct AnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};
  // TODO: Provide the mixing variables and binning directly via configurables (e.g. vectors of float)
  Configurable<string> fConfigMixingVariables{"cfgMixingVars", "", "Mixing configs separated by a comma, default no mixing"};
  Configurable<string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<bool> fConfigRunZorro{"cfgRunZorro", false, "Enable event selection with zorro [WARNING: under debug, do not enable!]"};
  Configurable<int> fConfigITSROFrameStartBorderMargin{"ITSROFrameStartBorderMargin", -1, "Number of bcs at the start of ITS RO Frame border. Take from CCDB if -1"};
  Configurable<int> fConfigITSROFrameEndBorderMargin{"ITSROFrameEndBorderMargin", -1, "Number of bcs at the end of ITS RO Frame border. Take from CCDB if -1"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;
  AnalysisCompositeCut* fEventCut;
  int fLastRun;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    VarManager::SetDefaultVarNames();
    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;", fConfigAddEventHistogram); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                           // provide the list of required variables so that VarManager knows what to fill
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

    // CCDB configuration
    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    // Not later than now objects
    // fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);

    fLastRun = -1;
  }

  template <uint32_t TEventFillMap, typename TEvent>
  void runEventSelection(TEvent const& event)
  {
    if (event.runNumber() != fLastRun) {
      auto alppar = fCCDB->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", event.timestamp());
      EventSelectionParams* par = fCCDB->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", event.timestamp());
      int itsROFrameStartBorderMargin = fConfigITSROFrameStartBorderMargin < 0 ? par->fITSROFrameStartBorderMargin : fConfigITSROFrameStartBorderMargin;
      int itsROFrameEndBorderMargin = fConfigITSROFrameEndBorderMargin < 0 ? par->fITSROFrameEndBorderMargin : fConfigITSROFrameEndBorderMargin;
      VarManager::SetITSROFBorderselection(alppar->roFrameBiasInBC, alppar->roFrameLengthInBC, itsROFrameStartBorderMargin, itsROFrameEndBorderMargin);
      fLastRun = event.runNumber();
    }

    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);

    VarManager::FillEvent<TEventFillMap>(event);
    // TODO: make this condition at compile time
    if (fConfigQA) {
      fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues); // automatically fill all the histograms in the class Event
    }

    if (!fConfigRunZorro) {
      if (fEventCut->IsSelected(VarManager::fgValues)) {
        if (fConfigQA) {
          fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
        }
        eventSel(1);
      } else {
        eventSel(0);
      }
    } else {
      if (fEventCut->IsSelected(VarManager::fgValues) && event.tag_bit(56)) { // This is the bit used for the software trigger event selections [TO BE DONE: find a more clear way to use it]
        if (fConfigQA) {
          fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
        }
        eventSel(1);
      } else {
        eventSel(0);
      }
    }

    if (fMixHandler != nullptr) {
      int hh = fMixHandler->FindEventCategory(VarManager::fgValues);
      hash(hh);
    }
  }

  void processSkimmed(MyEvents::iterator const& event)
  {
    runEventSelection<gkEventFillMap>(event);
  }
  void processSkimmedQVector(MyEventsQvector::iterator const& event)
  {
    runEventSelection<gkEventFillMapWithQvector>(event);
  }
  void processSkimmedQVectorCentr(MyEventsQvectorCentr::iterator const& event)
  {
    runEventSelection<gkEventFillMapWithQvectorCentr>(event);
  }
  void processSkimmedQVectorMultExtra(MyEventsQvectorMultExtra::iterator const& event)
  {
    runEventSelection<gkEventFillMapWithQvectorMultExtra>(event);
  }
  void processSkimmedQVectorCentrMultExtra(MyEventsQvectorCentrMultExtra::iterator const& event)
  {
    runEventSelection<gkEventFillMapWithQvectorCentrMultExtra>(event);
  }
  void processSkimmedQVectorExtraRef(MyEventsVtxCovQvectorExtraWithRefFlow::iterator const& event)
  {
    runEventSelection<gkEventFillMapWithCovQvectorExtraWithRefFlow>(event);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedQVector, "Run event selection on DQ skimmed events with Q vector from GFW", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedQVectorCentr, "Run event selection on DQ skimmed events with Q vector from CFW", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedQVectorMultExtra, "Run event selection on DQ skimmed events with Q vector from GFW and MultPV", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedQVectorCentrMultExtra, "Run event selection on DQ skimmed events with Q vector from CFW and MultPV", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedQVectorExtraRef, "Run event selection on DQ skimmed events with Q vector and subscribing to reference flow table", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDummy, "Dummy function", false);
  // TODO: Add process functions subscribing to Framework Collision
};

struct AnalysisTrackSelection {
  Produces<aod::BarrelTrackCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  // The list of cuts should contain all the track cuts needed later in analysis, including
  //  for candidate electron selection (+ eventual prefilter cuts) and other needs like quarkonium - hadron correlations
  // The user must ensure using them properly in the tasks downstream
  // NOTE: For now, the candidate electron cuts must be provided first, then followed by any other needed selections
  Configurable<string> fConfigCuts{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
  Configurable<int> fConfigPrefilterCutId{"cfgPrefilterCutId", 32, "Id of the Prefilter track cut (starting at 0)"}; // In order to create another column prefilter (should be temporary before improving cut selection in configurables, then displaced to AnalysisPrefilterSelection)
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<string> fConfigCcdbPathTPC{"ccdb-path-tpc", "Users/z/zhxiong/TPCPID/PostCalib", "base path to the ccdb object"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<bool> fConfigComputeTPCpostCalib{"cfgTPCpostCalib", false, "If true, compute TPC post-calibrated n-sigmas"};
  Configurable<std::string> fConfigRunPeriods{"cfgRunPeriods", "LHC22f", "run periods for used data"};
  Configurable<bool> fConfigDummyRunlist{"cfgDummyRunlist", false, "If true, use dummy runlist"};
  Configurable<int> fConfigInitRunNumber{"cfgInitRunNumber", 543215, "Initial run number used in run by run checks"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fTrackCuts;

  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

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

      DefineHistograms(fHistMan, histDirNames.Data(), fConfigAddTrackHistogram); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                           // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
    if (fConfigDummyRunlist) {
      VarManager::SetDummyRunlist(fConfigInitRunNumber);
    }
    if (fConfigComputeTPCpostCalib) {
      // CCDB configuration
      fCCDB->setURL(fConfigCcdbUrl.value);
      fCCDB->setCaching(true);
      fCCDB->setLocalObjectValidityChecking();
      // Not later than now objects
      fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
    }
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runTrackSelection(TEvent const& event, TTracks const& tracks)
  {
    VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
    // fill event information which might be needed in histograms/cuts that combine track and event properties
    VarManager::FillEvent<TEventFillMap>(event);

    // check whether the run changed, and if so, update calibrations in the VarManager
    // TODO: Here, for the run number and timestamp we assume the function runs with the
    //      DQ skimmed model. However, we need a compile time check so to make this compatible
    //      also with the full data model.
    if (fConfigComputeTPCpostCalib && fCurrentRun != event.runNumber()) {
      auto calibList = fCCDB->getForTimeStamp<TList>(fConfigCcdbPathTPC.value, event.timestamp());
      VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
      VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
      VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
      VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
      VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
      VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
      fCurrentRun = event.runNumber();
    }

    trackSel.reserve(tracks.size());
    uint32_t filterMap = 0;
    bool prefilterSelected = false;
    int iCut = 0;

    for (auto& track : tracks) {
      filterMap = 0;
      prefilterSelected = false;
      VarManager::FillTrack<TTrackFillMap>(track);
      if (fConfigQA) { // TODO: make this compile time
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }
      iCut = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, iCut++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          if (iCut != fConfigPrefilterCutId) {
            filterMap |= (static_cast<uint32_t>(1) << iCut);
          }
          if (iCut == fConfigPrefilterCutId) {
            prefilterSelected = true;
          }
          if (fConfigQA) { // TODO: make this compile time
            fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
          }
        }
      }

      trackSel(static_cast<int>(filterMap), static_cast<int>(prefilterSelected));
    } // end loop over tracks
  }

  void processSkimmed(MyEvents::iterator const& event, MyBarrelTracks const& tracks)
  {
    runTrackSelection<gkEventFillMap, gkTrackFillMap>(event, tracks);
  }
  void processSkimmedWithCov(MyEventsVtxCov::iterator const& event, MyBarrelTracksWithCov const& tracks)
  {
    runTrackSelection<gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, tracks);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmed, "Run barrel track selection on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmedWithCov, "Run barrel track selection on DQ skimmed tracks w/ cov matrix", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummy, "Dummy function", false);
};

struct AnalysisMuonSelection {
  Produces<aod::MuonTrackCuts> muonSel;
  OutputObj<THashList> fOutputList{"output"};
  Configurable<string> fConfigCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddMuonHistogram{"cfgAddMuonHistogram", "", "Comma separated list of histograms"};

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fMuonCuts;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

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

      DefineHistograms(fHistMan, histDirNames.Data(), fConfigAddMuonHistogram); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                          // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvent, typename TMuons>
  void runMuonSelection(TEvent const& event, TMuons const& muons)
  {
    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
    VarManager::FillEvent<TEventFillMap>(event);

    muonSel.reserve(muons.size());
    uint32_t filterMap = 0;
    int iCut = 0;

    for (auto& muon : muons) {
      filterMap = 0;
      VarManager::FillTrack<TMuonFillMap>(muon);
      if (fConfigQA) { // TODO: make this compile time
        fHistMan->FillHistClass("TrackMuon_BeforeCuts", VarManager::fgValues);
      }

      iCut = 0;
      for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, iCut++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (static_cast<uint32_t>(1) << iCut);
          if (fConfigQA) { // TODO: make this compile time
            fHistMan->FillHistClass(Form("TrackMuon_%s", (*cut).GetName()), VarManager::fgValues);
          }
        }
      }
      muonSel(static_cast<int>(filterMap));
    } // end loop over tracks
  }

  void processSkimmed(MyEvents::iterator const& event, MyMuonTracks const& muons)
  {
    runMuonSelection<gkEventFillMap, gkMuonFillMap>(event, muons);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisMuonSelection, processSkimmed, "Run muon selection on DQ skimmed muons", false);
  PROCESS_SWITCH(AnalysisMuonSelection, processDummy, "Dummy function", false);
};

struct AnalysisPrefilterSelection {
  SliceCache cache;
  Produces<aod::Prefilter> prefilter;
  Preslice<MyBarrelTracks> perCollision = aod::reducedtrack::reducedeventId;

  // Configurables
  Configurable<std::string> fConfigPrefilterPairCut{"cfgPrefilterPairCut", "", "Prefilter pair cut"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  o2::parameters::GRPMagField* grpmag = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  Filter barrelTracksSelectedPrefilter = aod::dqanalysisflags::isBarrelSelectedPrefilter > 0;

  Partition<soa::Filtered<MyBarrelTracksSelected>> barrelTracksSelected = aod::dqanalysisflags::isBarrelSelected > 0;

  std::map<int, bool> fPrefiltermap;
  AnalysisCompositeCut* fPairCut;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    fCurrentRun = 0;

    ccdb->setURL(ccdburl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    fPairCut = new AnalysisCompositeCut(true);
    TString pairCutStr = fConfigPrefilterPairCut.value;
    if (!pairCutStr.IsNull()) {
      fPairCut = dqcuts::GetCompositeCut(pairCutStr.Data());
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    VarManager::SetDefaultVarNames();

    // Do we need this ?
    VarManager::SetupTwoProngDCAFitter(5.0f, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, true); // TODO: get these parameters from Configurables
    VarManager::SetupTwoProngFwdDCAFitter(5.0f, true, 200.0f, 1.0e-3f, 0.9f, true);
  }

  template <int TPairType, uint32_t TTrackFillMap, typename TTracks1, typename TTracks2>
  void runPrefilterPairing(TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
      if (track1.sign() * track2.sign() > 0) {
        continue;
      }

      // pairing
      VarManager::FillPair<TPairType, TTrackFillMap>(track1, track2);

      if (fPairCut->IsSelected(VarManager::fgValues)) {
        fPrefiltermap[track1.globalIndex()] = true;
        fPrefiltermap[track2.globalIndex()] = true;
      }
    }
  }

  void processBarrelSkimmed(MyEventsSelected const& events, soa::Filtered<MyBarrelTracksSelected> const& filteredTracks, MyBarrelTracks const& tracks)
  {
    const int pairType = VarManager::kDecayToEE;
    fPrefiltermap.clear();

    if (events.size() > 0 && fCurrentRun != events.begin().runNumber()) {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, events.begin().timestamp());
      if (grpmag != nullptr) {
        VarManager::SetMagneticField(grpmag->getNominalL3Field());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", events.begin().timestamp());
      }
      fCurrentRun = events.begin().runNumber();
    }

    for (auto& event : events) {
      if (event.isEventSelected()) {
        auto groupedPrefilterCandidates = filteredTracks.sliceBy(perCollision, event.globalIndex());
        auto groupedBarrelCandidates = barrelTracksSelected->sliceByCached(aod::reducedtrack::reducedeventId, event.globalIndex(), cache);
        runPrefilterPairing<pairType, gkTrackFillMap>(groupedPrefilterCandidates, groupedBarrelCandidates);
      }
    } // end loop events

    // Fill Prefilter bits for all tracks to have something joinable to MyBarrelTracksSelected
    for (auto& track : tracks) {
      prefilter(static_cast<int>(fPrefiltermap[track.globalIndex()]));
    }
  }

  void processDummy(MyEvents&)
  {
  }

  PROCESS_SWITCH(AnalysisPrefilterSelection, processBarrelSkimmed, "Run Prefilter selection on reduced tracks", false);
  PROCESS_SWITCH(AnalysisPrefilterSelection, processDummy, "Do nothing", false);
};

struct AnalysisEventMixing {
  OutputObj<THashList> fOutputList{"output"};
  // Here one should provide the list of electron and muon candidate cuts in the same order as specified in the above
  // single particle selection tasks to preserve the correspondence between the track cut name and its
  //  bit position in the cuts bitmap
  // TODO: Create a configurable to specify exactly on which of the bits one should run the event mixing
  Configurable<string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<string> fConfigMuonCuts{"cfgMuonCuts", "", "Comma separated list of muon cuts"};
  Configurable<int> fConfigMixingDepth{"cfgMixingDepth", 100, "Number of Events stored for event mixing"};
  Configurable<std::string> fConfigAddEventMixingHistogram{"cfgAddEventMixingHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> fConfigAmbiguousHist{"cfgAmbiHist", false, "Enable Ambiguous histograms for time association studies"};
  Configurable<string> ccdbPathFlow{"ccdb-path-flow", "Users/c/chizh/FlowResolution", "path to the ccdb object for flow resolution factors"};
  Configurable<bool> fConfigFlowReso{"cfgFlowReso", false, "Enable loading of flow resolution factors from CCDB"};
  Configurable<bool> fConfigSingleMuCumulants{"cfgSingleMuCumulants", false, "Enable loading of flow resolution factors from CCDB"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  o2::parameters::GRPMagField* grpmag = nullptr;
  TH1D* ResoFlowSP = nullptr; // Resolution factors for flow analysis, this will be loaded from CCDB
  TH1D* ResoFlowEP = nullptr; // Resolution factors for flow analysis, this will be loaded from CCDB
  TH2D* SingleMuv22m = nullptr; // Single muon v22, loaded from CCDB
  TH2D* SingleMuv24m = nullptr; // Single muon v24, loaded from CCDB
  TH2D* SingleMuv22p = nullptr; // Single antimuon v22, loaded from CCDB
  TH2D* SingleMuv24p = nullptr; // Single antimuon v24, loaded from CCDB
  int fCurrentRun;            // needed to detect if the run changed and trigger update of calibrations etc.

  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected == 1;
  Filter filterTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;
  Filter filterMuonTrackSelected = aod::dqanalysisflags::isMuonSelected > 0;

  HistogramManager* fHistMan;
  // NOTE: The bit mask is required to run pairing just based on the desired electron/muon candidate cuts
  uint32_t fTwoTrackFilterMask = 0;
  uint32_t fTwoMuonFilterMask = 0;
  std::vector<std::vector<TString>> fTrackHistNames;
  std::vector<std::vector<TString>> fMuonHistNames;
  std::vector<std::vector<TString>> fTrackMuonHistNames;

  NoBinningPolicy<aod::dqanalysisflags::MixingHash> hashBin;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    fCurrentRun = 0;

    ccdb->setURL(ccdburl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Keep track of all the histogram class names to avoid composing strings in the event mixing pairing
    TString histNames = "";
    if (context.mOptions.get<bool>("processBarrelSkimmed") || context.mOptions.get<bool>("processBarrelVnSkimmed")) {
      TString cutNames = fConfigTrackCuts.value;
      if (!cutNames.IsNull()) {
        std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
        for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
          std::vector<TString> names = {
            Form("PairsBarrelMEPM_%s", objArray->At(icut)->GetName()),
            Form("PairsBarrelMEPP_%s", objArray->At(icut)->GetName()),
            Form("PairsBarrelMEMM_%s", objArray->At(icut)->GetName())};
          histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          fTrackHistNames.push_back(names);
          fTwoTrackFilterMask |= (static_cast<uint32_t>(1) << icut);
        }
      }
    }
    if (context.mOptions.get<bool>("processMuonSkimmed") || context.mOptions.get<bool>("processMuonVnSkimmed") || context.mOptions.get<bool>("processMuonVnCentrSkimmed") || context.mOptions.get<bool>("processMuonVnExtraSkimmed")) {
      TString cutNames = fConfigMuonCuts.value;
      if (!cutNames.IsNull()) {
        std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
        for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
          std::vector<TString> names = {
            Form("PairsMuonMEPM_%s", objArray->At(icut)->GetName()),
            Form("PairsMuonMEPP_%s", objArray->At(icut)->GetName()),
            Form("PairsMuonMEMM_%s", objArray->At(icut)->GetName())};
          if (fConfigAmbiguousHist) {
            histNames += Form("%s;%s;%s;%s_unambiguous;%s_unambiguous;%s_unambiguous;", names[0].Data(), names[1].Data(), names[2].Data(), names[0].Data(), names[1].Data(), names[2].Data());
          } else {
            histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          }
          fMuonHistNames.push_back(names);
          fTwoMuonFilterMask |= (static_cast<uint32_t>(1) << icut);
        }
      }
    }
    if (context.mOptions.get<bool>("processBarrelMuonSkimmed")) {
      TString cutNamesBarrel = fConfigTrackCuts.value;
      TString cutNamesMuon = fConfigMuonCuts.value;
      if (!cutNamesBarrel.IsNull() && !cutNamesMuon.IsNull()) {
        std::unique_ptr<TObjArray> objArrayBarrel(cutNamesBarrel.Tokenize(","));
        std::unique_ptr<TObjArray> objArrayMuon(cutNamesMuon.Tokenize(","));
        if (objArrayBarrel->GetEntries() == objArrayMuon->GetEntries()) { // one must specify equal number of barrel and muon cuts
          for (int icut = 0; icut < objArrayBarrel->GetEntries(); ++icut) {
            std::vector<TString> names = {
              Form("PairsEleMuMEPM_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName()),
              Form("PairsEleMuMEPP_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName()),
              Form("PairsEleMuMEMM_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName())};
            histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
            fTrackMuonHistNames.push_back(names);
            fTwoTrackFilterMask |= (static_cast<uint32_t>(1) << icut);
            fTwoMuonFilterMask |= (static_cast<uint32_t>(1) << icut);
          }
        }
      }
    }

    DefineHistograms(fHistMan, histNames.Data(), fConfigAddEventMixingHistogram); // define all histograms
    // Additional histograms via JSON
    dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str());
    VarManager::SetUseVars(fHistMan->GetUsedVars());                              // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  template <uint32_t TEventFillMap, int TPairType, typename TTracks1, typename TTracks2>
  void runMixedPairing(TTracks1 const& tracks1, TTracks2 const& tracks2)
  {

    unsigned int ncuts = fTrackHistNames.size();
    std::vector<std::vector<TString>> histNames = fTrackHistNames;
    if constexpr (TPairType == pairTypeMuMu) {
      ncuts = fMuonHistNames.size();
      histNames = fMuonHistNames;
    }
    if constexpr (TPairType == pairTypeEMu) {
      ncuts = fTrackMuonHistNames.size();
      histNames = fTrackMuonHistNames;
    }

    uint32_t twoTrackFilter = 0;
    uint32_t mult_dimuons = 0;
    for (auto& track1 : tracks1) {
      for (auto& track2 : tracks2) {
        if constexpr (TPairType == VarManager::kDecayToEE) {
          twoTrackFilter = static_cast<uint32_t>(track1.isBarrelSelected()) & static_cast<uint32_t>(track2.isBarrelSelected()) & fTwoTrackFilterMask;
        }
        if constexpr (TPairType == VarManager::kDecayToMuMu) {
          twoTrackFilter = static_cast<uint32_t>(track1.isMuonSelected()) & static_cast<uint32_t>(track2.isMuonSelected()) & fTwoMuonFilterMask;
          if (fConfigSingleMuCumulants) {
            VarManager::FillTwoMixEventsCumulants(SingleMuv22m, SingleMuv24m, SingleMuv22p, SingleMuv24p, track1, track2);
          }
        }
        if constexpr (TPairType == VarManager::kElectronMuon) {
          twoTrackFilter = static_cast<uint32_t>(track1.isBarrelSelected()) & static_cast<uint32_t>(track2.isMuonSelected()) & fTwoTrackFilterMask;
        }

        if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
          continue;
        }
        VarManager::FillPairME<TEventFillMap, TPairType>(track1, track2);

        for (unsigned int icut = 0; icut < ncuts; icut++) {
          if (twoTrackFilter & (static_cast<uint32_t>(1) << icut)) {
            if (track1.sign() * track2.sign() < 0) {
              mult_dimuons++;
              fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
              if (fConfigAmbiguousHist && !(track1.isAmbiguous() || track2.isAmbiguous())) {
                fHistMan->FillHistClass(Form("%s_unambiguous", histNames[icut][0].Data()), VarManager::fgValues);
              }
            } else {
              if (track1.sign() > 0) {
                fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
                if (fConfigAmbiguousHist && !(track1.isAmbiguous() || track2.isAmbiguous())) {
                  fHistMan->FillHistClass(Form("%s_unambiguous", histNames[icut][1].Data()), VarManager::fgValues);
                }
              } else {
                fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
                if (fConfigAmbiguousHist && !(track1.isAmbiguous() || track2.isAmbiguous())) {
                  fHistMan->FillHistClass(Form("%s_unambiguous", histNames[icut][2].Data()), VarManager::fgValues);
                }
              }
            }
          } // end if (filter bits)
        } // end for (cuts)
      } // end for (track2)
    } // end for (track1)
    VarManager::fgValues[VarManager::kMultDimuonsME] = mult_dimuons;
  }

  // barrel-barrel and muon-muon event mixing
  template <int TPairType, uint32_t TEventFillMap, typename TEvents, typename TTracks>
  void runSameSide(TEvents& events, TTracks const& tracks, Preslice<TTracks>& preSlice)
  {
    if (events.size() > 0 && fCurrentRun != events.begin().runNumber()) {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, events.begin().timestamp());
      if (grpmag != nullptr) {
        VarManager::SetMagneticField(grpmag->getNominalL3Field());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", events.begin().timestamp());
      }
      if (fConfigFlowReso) {
        TString PathFlow = ccdbPathFlow.value;
        TString ccdbPathFlowSP = Form("%s/ScalarProduct", PathFlow.Data());
        TString ccdbPathFlowEP = Form("%s/EventPlane", PathFlow.Data());
        ResoFlowSP = ccdb->getForTimeStamp<TH1D>(ccdbPathFlowSP.Data(), events.begin().timestamp());
        ResoFlowEP = ccdb->getForTimeStamp<TH1D>(ccdbPathFlowEP.Data(), events.begin().timestamp());
        if (ResoFlowSP == nullptr || ResoFlowEP == nullptr) {
          LOGF(fatal, "Resolution factor is not available in CCDB at timestamp=%llu", events.begin().timestamp());
        }
      }
      if (fConfigSingleMuCumulants) {
        TString PathFlow = ccdbPathFlow.value;
        TString ccdbPathMuv22m = Form("%s/SingleMuv22m", PathFlow.Data());
        TString ccdbPathMuv24m = Form("%s/SingleMuv24m", PathFlow.Data());
        TString ccdbPathMuv22p = Form("%s/SingleMuv22p", PathFlow.Data());
        TString ccdbPathMuv24p = Form("%s/SingleMuv24p", PathFlow.Data());
        SingleMuv22m = ccdb->getForTimeStamp<TH2D>(ccdbPathMuv22m.Data(), events.begin().timestamp());
        SingleMuv24m = ccdb->getForTimeStamp<TH2D>(ccdbPathMuv24m.Data(), events.begin().timestamp());
        SingleMuv22p = ccdb->getForTimeStamp<TH2D>(ccdbPathMuv22p.Data(), events.begin().timestamp());
        SingleMuv24p = ccdb->getForTimeStamp<TH2D>(ccdbPathMuv24p.Data(), events.begin().timestamp());
        if (SingleMuv22m == nullptr || SingleMuv24m == nullptr || SingleMuv22p == nullptr || SingleMuv24p == nullptr) {
          LOGF(fatal, "Single muon cumulants are not available in CCDB at timestamp=%llu", events.begin().timestamp());
        }
      }
      fCurrentRun = events.begin().runNumber();
    }

    events.bindExternalIndices(&tracks);
    int mixingDepth = fConfigMixingDepth.value;
    for (auto& [event1, event2] : selfCombinations(hashBin, mixingDepth, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event1, VarManager::fgValues);

      auto tracks1 = tracks.sliceBy(preSlice, event1.globalIndex());
      tracks1.bindExternalIndices(&events);

      auto tracks2 = tracks.sliceBy(preSlice, event2.globalIndex());
      tracks2.bindExternalIndices(&events);

      VarManager::FillTwoMixEvents<TEventFillMap>(event1, event2, tracks1, tracks2);
      if (fConfigFlowReso) {
        VarManager::FillTwoMixEventsFlowResoFactor(ResoFlowSP, ResoFlowEP);
      }
      runMixedPairing<TEventFillMap, TPairType>(tracks1, tracks2);
    } // end event loop
  }

  // barrel-muon event mixing
  template <uint32_t TEventFillMap, typename TEvents, typename TTracks, typename TMuons>
  void runBarrelMuon(TEvents& events, TTracks const& tracks, TMuons const& muons)
  {
    if (events.size() > 0 && fCurrentRun != events.begin().runNumber()) {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, events.begin().timestamp());
      if (grpmag != nullptr) {
        VarManager::SetMagneticField(grpmag->getNominalL3Field());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", events.begin().timestamp());
      }
      fCurrentRun = events.begin().runNumber();
    }

    events.bindExternalIndices(&muons);

    for (auto& [event1, event2] : selfCombinations(hashBin, 100, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event1, VarManager::fgValues);

      auto tracks1 = tracks.sliceBy(perEventsSelectedT, event1.globalIndex());
      tracks1.bindExternalIndices(&events);

      auto muons2 = muons.sliceBy(perEventsSelectedM, event2.globalIndex());
      muons2.bindExternalIndices(&events);

      runMixedPairing<TEventFillMap, pairTypeEMu>(tracks1, muons2);
    } // end event loop
  }

  Preslice<soa::Filtered<MyBarrelTracksSelected>> perEventsSelectedT = aod::reducedtrack::reducedeventId;
  Preslice<soa::Filtered<MyMuonTracksSelected>> perEventsSelectedM = aod::reducedmuon::reducedeventId;

  void processBarrelSkimmed(soa::Filtered<MyEventsHashSelected>& events, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    runSameSide<pairTypeEE, gkEventFillMap>(events, tracks, perEventsSelectedT);
  }
  void processMuonSkimmed(soa::Filtered<MyEventsHashSelected>& events, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    runSameSide<pairTypeMuMu, gkEventFillMap>(events, muons, perEventsSelectedM);
  }
  void processBarrelMuonSkimmed(soa::Filtered<MyEventsHashSelected>& events, soa::Filtered<MyBarrelTracksSelected> const& tracks, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    runBarrelMuon<gkEventFillMap>(events, tracks, muons);
  }
  void processBarrelVnSkimmed(soa::Filtered<MyEventsHashSelectedQvector>& events, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    runSameSide<pairTypeEE, gkEventFillMapWithQvector>(events, tracks, perEventsSelectedT);
  }
  void processMuonVnSkimmed(soa::Filtered<MyEventsHashSelectedQvector>& events, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    runSameSide<pairTypeMuMu, gkEventFillMapWithQvector>(events, muons, perEventsSelectedM);
  }
  void processMuonVnCentrSkimmed(soa::Filtered<MyEventsHashSelectedQvectorCentr>& events, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    runSameSide<pairTypeMuMu, gkEventFillMapWithQvectorCentr>(events, muons, perEventsSelectedM);
  }
  void processMuonVnExtraSkimmed(soa::Filtered<MyEventsHashSelectedQvectorExtra>& events, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    runSameSide<pairTypeMuMu, gkEventFillMapWithCovQvectorExtraWithRefFlow>(events, muons, perEventsSelectedM);
  }
  // TODO: This is a dummy process function for the case when the user does not want to run any of the process functions (no event mixing)
  //    If there is no process function enabled, the workflow hangs
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventMixing, processBarrelSkimmed, "Run barrel-barrel mixing on skimmed tracks", false);
  PROCESS_SWITCH(AnalysisEventMixing, processMuonSkimmed, "Run muon-muon mixing on skimmed muons", false);
  PROCESS_SWITCH(AnalysisEventMixing, processBarrelMuonSkimmed, "Run barrel-muon mixing on skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisEventMixing, processBarrelVnSkimmed, "Run barrel-barrel vn mixing on skimmed tracks", false);
  PROCESS_SWITCH(AnalysisEventMixing, processMuonVnSkimmed, "Run muon-muon vn mixing on skimmed tracks", false);
  PROCESS_SWITCH(AnalysisEventMixing, processMuonVnCentrSkimmed, "Run muon-muon vn mixing on skimmed tracks from central framework", false);
  PROCESS_SWITCH(AnalysisEventMixing, processMuonVnExtraSkimmed, "Run muon-muon vn mixing on skimmed tracks from GFW", false);
  PROCESS_SWITCH(AnalysisEventMixing, processDummy, "Dummy function", false);
};

struct AnalysisSameEventPairing {

  Produces<aod::Dielectrons> dielectronList;
  Produces<aod::Dimuons> dimuonList;
  Produces<aod::DielectronsExtra> dielectronExtraList;
  Produces<aod::DielectronsInfo> dielectronInfoList;
  Produces<aod::DimuonsExtra> dimuonExtraList;
  Produces<aod::DielectronsAll> dielectronAllList;
  Produces<aod::DimuonsAll> dimuonAllList;
  Produces<aod::DileptonFlow> dileptonFlowList;
  Produces<aod::DileptonsInfo> dileptonInfoList;
  Produces<aod::DileptonsMiniTree> dileptonMiniTree;
  float mMagField = 0.0;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  TH1D* ResoFlowSP = nullptr; // Resolution factors for flow analysis, this will be loaded from CCDB
  TH1D* ResoFlowEP = nullptr; // Resolution factors for flow analysis, this will be loaded from CCDB
  int fCurrentRun;            // needed to detect if the run changed and trigger update of calibrations etc.

  OutputObj<THashList> fOutputList{"output"};
  Configurable<string> fConfigTrackCuts{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  Configurable<string> fConfigMuonCuts{"cfgMuonCuts", "", "Comma separated list of muon cuts"};
  Configurable<string> fConfigPairCuts{"cfgPairCuts", "", "Comma separated list of pair cuts"};
  Configurable<string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<string> ccdbPath{"ccdb-path", "Users/lm", "base path to the ccdb object"};
  Configurable<string> ccdbPathFlow{"ccdb-path-flow", "Users/c/chizh/FlowResolution", "path to the ccdb object for flow resolution factors"};
  Configurable<bool> fConfigFlowReso{"cfgFlowReso", false, "Enable loading of flow resolution factors from CCDB"};
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> fConfigAddSEPHistogram{"cfgAddSEPHistogram", "", "Comma separated list of histograms"};
  Configurable<bool> fConfigFlatTables{"cfgFlatTables", false, "Produce a single flat tables with all relevant information of the pairs and single tracks"};
  Configurable<bool> fConfigAmbiguousHist{"cfgAmbiHist", false, "Enable Ambiguous histograms for time association studies"};
  Configurable<bool> fConfigMultDimuons{"cfgMultDimuons", false, "Multiplicity for Unlike Dimuons"};
  Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
  Configurable<bool> fUseRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> fUseAbsDCA{"cfgUseAbsDCA", false, "Use absolute DCA minimization instead of chi^2 minimization in secondary vertexing"};
  Configurable<bool> fPropToPCA{"cfgPropToPCA", false, "Propagate tracks to secondary vertex"};
  Configurable<bool> fCorrFullGeo{"cfgCorrFullGeo", false, "Use full geometry to correct for MCS effects in track propagation"};
  Configurable<bool> fNoCorr{"cfgNoCorrFwdProp", false, "Do not correct for MCS effects in track propagation"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> fCollisionSystem{"syst", "pp", "Collision system, pp or PbPb"};
  Configurable<float> fCenterMassEnergy{"energy", 13600, "Center of mass energy in GeV"};
  Configurable<bool> fConfigCumulants{"cfgCumulants", false, "If true, fill Cumulants with Weights different than 0"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};

  // Configurables to create output tree (flat tables or minitree)
  struct : ConfigurableGroup {
    Configurable<bool> fConfigMiniTree{"useMiniTree.cfgMiniTree", false, "Produce a single flat table with minimal information for analysis"};
    Configurable<float> fConfigMiniTreeMinMass{"useMiniTree.cfgMiniTreeMinMass", 2, "Min. mass cut for minitree"};
    Configurable<float> fConfigMiniTreeMaxMass{"useMiniTree.cfgMiniTreeMaxMass", 5, "Max. mass cut for minitree"};
  } useMiniTree;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected == 1;
  // NOTE: the barrel filter map contains decisions for both electrons and hadrons used in the correlation task
  Filter filterBarrelTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;
  Filter filterMuonTrackSelected = aod::dqanalysisflags::isMuonSelected > 0;
  Filter prefilter = aod::dqanalysisflags::isPrefilterVetoed == 0;

  HistogramManager* fHistMan;

  // NOTE: The track filter produced by the barrel track selection contain a number of electron cut decisions and one last cut for hadrons used in the
  //           dilepton - hadron task downstream. So the bit mask is required to select pairs just based on the electron cuts
  // TODO: provide as Configurable the list and names of the cuts which should be used in pairing
  uint32_t fTwoTrackFilterMask = 0;
  uint32_t fTwoMuonFilterMask = 0;
  std::vector<std::vector<TString>> fTrackHistNames;
  std::vector<std::vector<TString>> fMuonHistNames;
  std::vector<std::vector<TString>> fTrackMuonHistNames;
  std::vector<AnalysisCompositeCut> fPairCuts;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    fCurrentRun = 0;

    ccdb->setURL(ccdburl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    if (fNoCorr) {
      VarManager::SetupFwdDCAFitterNoCorr();
    } else if (fCorrFullGeo || (fConfigUseKFVertexing && fPropToPCA)) {
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
    } else {
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
      VarManager::SetupMatLUTFwdDCAFitter(lut);
    }

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Keep track of all the histogram class names to avoid composing strings in the event mixing pairing
    TString histNames = "";
    std::vector<TString> names;

    TString cutNamesStr = fConfigPairCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    if (context.mOptions.get<bool>("processDecayToEESkimmed") || context.mOptions.get<bool>("processDecayToEESkimmedNoTwoProngFitter") || context.mOptions.get<bool>("processDecayToEESkimmedWithCov") || context.mOptions.get<bool>("processDecayToEESkimmedWithCovNoTwoProngFitter") || context.mOptions.get<bool>("processDecayToEEVertexingSkimmed") || context.mOptions.get<bool>("processVnDecayToEESkimmed") || context.mOptions.get<bool>("processDecayToEEPrefilterSkimmed") || context.mOptions.get<bool>("processDecayToEEPrefilterSkimmedNoTwoProngFitter") || context.mOptions.get<bool>("processDecayToEESkimmedWithColl") || context.mOptions.get<bool>("processDecayToEESkimmedWithCollNoTwoProngFitter") || context.mOptions.get<bool>("processDecayToPiPiSkimmed") || context.mOptions.get<bool>("processAllSkimmed")) {
      TString cutNames = fConfigTrackCuts.value;
      if (!cutNames.IsNull()) { // if track cuts
        std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
        for (int icut = 0; icut < objArray->GetEntries(); ++icut) { // loop over track cuts
          fTwoTrackFilterMask |= (static_cast<uint32_t>(1) << icut);
          // no pair cuts
          names = {
            Form("PairsBarrelSEPM_%s", objArray->At(icut)->GetName()),
            Form("PairsBarrelSEPP_%s", objArray->At(icut)->GetName()),
            Form("PairsBarrelSEMM_%s", objArray->At(icut)->GetName())};
          histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          fTrackHistNames.push_back(names);

          TString cutNamesStr = fConfigPairCuts.value;
          if (!cutNamesStr.IsNull()) { // if pair cuts
            std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
            for (int iPairCut = 0; iPairCut < objArrayPair->GetEntries(); ++iPairCut) { // loop over pair cuts
              names = {
                Form("PairsBarrelSEPM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                Form("PairsBarrelSEPP_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                Form("PairsBarrelSEMM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName())};
              if (fConfigAmbiguousHist) {
                histNames += Form("%s;%s;%s;%s_unambiguous;%s_unambiguous;%s_unambiguous;", names[0].Data(), names[1].Data(), names[2].Data(), names[0].Data(), names[1].Data(), names[2].Data());
              } else {
                histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
              }
              fTrackHistNames.push_back(names);
            } // end loop (pair cuts)
          } // end if (pair cuts)
        } // end loop (track cuts)
      } // end if (track cuts)
    }

    if (context.mOptions.get<bool>("processDecayToMuMuSkimmed") || context.mOptions.get<bool>("processDecayToMuMuSkimmedWithMult") || context.mOptions.get<bool>("processDecayToMuMuVertexingSkimmed") || context.mOptions.get<bool>("processDecayToMuMuSkimmedWithColl") || context.mOptions.get<bool>("processVnDecayToMuMuSkimmed") || context.mOptions.get<bool>("processVnDecayToMuMuSkimmedWithWeights") || context.mOptions.get<bool>("processVnDecayToMuMuSkimmedWithWeightsAndColl") || context.mOptions.get<bool>("processVnCentrDecayToMuMuSkimmed") || context.mOptions.get<bool>("processAllSkimmed")) {
      TString cutNames = fConfigMuonCuts.value;
      if (!cutNames.IsNull()) {
        std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
        for (int icut = 0; icut < objArray->GetEntries(); ++icut) { // loop over track cuts
          fTwoMuonFilterMask |= (static_cast<uint32_t>(1) << icut);
          // no pair cuts
          names = {
            Form("PairsMuonSEPM_%s", objArray->At(icut)->GetName()),
            Form("PairsMuonSEPP_%s", objArray->At(icut)->GetName()),
            Form("PairsMuonSEMM_%s", objArray->At(icut)->GetName())};
          if (fConfigAmbiguousHist) {
            histNames += Form("%s;%s;%s;%s_unambiguous;%s_unambiguous;%s_unambiguous;", names[0].Data(), names[1].Data(), names[2].Data(), names[0].Data(), names[1].Data(), names[2].Data());
          } else {
            histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          }
          fMuonHistNames.push_back(names);

          TString cutNamesStr = fConfigPairCuts.value;
          if (!cutNamesStr.IsNull()) { // if pair cuts
            std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
            for (int iPairCut = 0; iPairCut < objArrayPair->GetEntries(); ++iPairCut) { // loop over pair cuts
              names = {
                Form("PairsMuonSEPM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                Form("PairsMuonSEPP_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                Form("PairsMuonSEMM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName())};
              histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
              fMuonHistNames.push_back(names);
            } // end loop (pair cuts)
          } // end if (pair cuts)
        } // end loop (track cuts)
      } // end if (track cuts)
    }
    if (context.mOptions.get<bool>("processElectronMuonSkimmed") || context.mOptions.get<bool>("processAllSkimmed")) {
      TString cutNamesBarrel = fConfigTrackCuts.value;
      TString cutNamesMuon = fConfigMuonCuts.value;
      if (!cutNamesBarrel.IsNull() && !cutNamesMuon.IsNull()) {
        std::unique_ptr<TObjArray> objArrayBarrel(cutNamesBarrel.Tokenize(","));
        std::unique_ptr<TObjArray> objArrayMuon(cutNamesMuon.Tokenize(","));
        if (objArrayBarrel->GetEntries() == objArrayMuon->GetEntries()) {   // one must specify equal number of barrel and muon cuts
          for (int icut = 0; icut < objArrayBarrel->GetEntries(); ++icut) { // loop over track cuts
            fTwoTrackFilterMask |= (static_cast<uint32_t>(1) << icut);
            fTwoMuonFilterMask |= (static_cast<uint32_t>(1) << icut);
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
            } // end if (pair cuts)
          } // end loop (track cuts)
        } // end if (equal number of cuts)
      } // end if (track cuts)
    }

    // Usage example of ccdb
    // ccdb->setURL(url.value);
    // ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // ccdb->setCreatedNotAfter(nolaterthan.value);

    VarManager::SetCollisionSystem((TString)fCollisionSystem, fCenterMassEnergy); // set collision system and center of mass energy

    DefineHistograms(fHistMan, histNames.Data(), fConfigAddSEPHistogram); // define all histograms
    dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
    VarManager::SetUseVars(fHistMan->GetUsedVars());                      // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  // Template function to run same event pairing (barrel-barrel, muon-muon, barrel-muon)
  template <bool TTwoProngFitter, int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks1, typename TTracks2>
  void runSameEventPairing(TEvent const& event, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    if (fConfigCumulants && VarManager::fgValues[VarManager::kM11REF] == 0) {

      return;
    }
    if (fCurrentRun != event.runNumber()) {
      if (fUseRemoteField) {
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, event.timestamp());
        if (grpmag != nullptr) {
          mMagField = grpmag->getNominalL3Field();
        } else {
          LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", event.timestamp());
        }
        if constexpr (TTwoProngFitter == true) {
          if (fConfigUseKFVertexing) {
            VarManager::SetupTwoProngKFParticle(mMagField);
          } else {
            VarManager::SetupTwoProngDCAFitter(mMagField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fUseAbsDCA.value);
            VarManager::SetupTwoProngFwdDCAFitter(mMagField, true, 200.0f, 1.0e-3f, 0.9f, fUseAbsDCA.value);
          }
        }
        VarManager::SetMagneticField(mMagField);
      } else {
        if constexpr (TTwoProngFitter == true) {
          if (fConfigUseKFVertexing) {
            VarManager::SetupTwoProngKFParticle(fConfigMagField.value);
          } else {
            VarManager::SetupTwoProngDCAFitter(fConfigMagField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fUseAbsDCA.value);
            VarManager::SetupTwoProngFwdDCAFitter(fConfigMagField.value, true, 200.0f, 1.0e-3f, 0.9f, fUseAbsDCA.value);
          }
        }
        VarManager::SetMagneticField(fConfigMagField.value);
      }
      if (fConfigFlowReso) {
        TString PathFlow = ccdbPathFlow.value;
        TString ccdbPathFlowSP = Form("%s/ScalarProduct", PathFlow.Data());
        TString ccdbPathFlowEP = Form("%s/EventPlane", PathFlow.Data());
        ResoFlowSP = ccdb->getForTimeStamp<TH1D>(ccdbPathFlowSP.Data(), event.timestamp());
        ResoFlowEP = ccdb->getForTimeStamp<TH1D>(ccdbPathFlowEP.Data(), event.timestamp());
        if (ResoFlowSP == nullptr || ResoFlowEP == nullptr) {
          LOGF(fatal, "Resolution factor is not available in CCDB at timestamp=%llu", event.timestamp());
        }
      }
      fCurrentRun = event.runNumber();
    }

    TString cutNames = fConfigTrackCuts.value;
    std::vector<std::vector<TString>> histNames = fTrackHistNames;
    if constexpr (TPairType == pairTypeMuMu) {
      cutNames = fConfigMuonCuts.value;
      histNames = fMuonHistNames;
    }
    if constexpr (TPairType == pairTypeEMu) {
      cutNames = fConfigMuonCuts.value;
      histNames = fTrackMuonHistNames;
    }
    std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
    int ncuts = objArray->GetEntries();

    uint32_t twoTrackFilter = 0;
    uint32_t dileptonFilterMap = 0;
    uint32_t dileptonMcDecision = 0; // placeholder, copy of the dqEfficiency.cxx one
    dielectronList.reserve(1);
    dimuonList.reserve(1);
    dielectronExtraList.reserve(1);
    dielectronInfoList.reserve(1);
    dimuonExtraList.reserve(1);
    dileptonInfoList.reserve(1);
    dileptonFlowList.reserve(1);
    if (fConfigFlatTables.value) {
      dielectronAllList.reserve(1);
      dimuonAllList.reserve(1);
    }
    if (useMiniTree.fConfigMiniTree) {
      dileptonMiniTree.reserve(1);
    }

    if (fConfigMultDimuons.value) {

      uint32_t mult_dimuons = 0;
      uint32_t mult_antimuons = 0;
      uint32_t mult_muons = 0;

      for (auto& t : tracks1) {
        if constexpr (TPairType == VarManager::kDecayToMuMu) {
          if (static_cast<uint32_t>(t.isMuonSelected()) & fTwoMuonFilterMask) {
            if (t.sign() < 0) {
              mult_muons++;
            } else {
              mult_antimuons++;
            }
          }
        }
      }

      for (auto& [t1, t2] : combinations(tracks1, tracks2)) {
        if constexpr (TPairType == VarManager::kDecayToMuMu) {
          twoTrackFilter = static_cast<uint32_t>(t1.isMuonSelected()) & static_cast<uint32_t>(t2.isMuonSelected()) & fTwoMuonFilterMask;
        }

        if (twoTrackFilter && (t1.sign() != t2.sign())) {
          mult_dimuons++;
        }
      }

      VarManager::fgValues[VarManager::kMultDimuons] = mult_dimuons;
      VarManager::fgValues[VarManager::kMultMuons] = mult_muons;
      VarManager::fgValues[VarManager::kMultAntiMuons] = mult_antimuons;
    }

    if (fConfigFlowReso) {
      VarManager::FillEventFlowResoFactor(ResoFlowSP, ResoFlowEP);
    }

    bool isFirst = true;
    for (auto& [t1, t2] : combinations(tracks1, tracks2)) {
      if constexpr (TPairType == VarManager::kDecayToEE || TPairType == VarManager::kDecayToPiPi) {
        twoTrackFilter = static_cast<uint32_t>(t1.isBarrelSelected()) & static_cast<uint32_t>(t2.isBarrelSelected()) & fTwoTrackFilterMask;
      }
      if constexpr (TPairType == VarManager::kDecayToMuMu) {
        twoTrackFilter = static_cast<uint32_t>(t1.isMuonSelected()) & static_cast<uint32_t>(t2.isMuonSelected()) & fTwoMuonFilterMask;
      }
      if constexpr (TPairType == VarManager::kElectronMuon) {
        twoTrackFilter = static_cast<uint32_t>(t1.isBarrelSelected()) & static_cast<uint32_t>(t2.isMuonSelected()) & fTwoTrackFilterMask;
      }
      if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }
      constexpr bool eventHasQvector = ((TEventFillMap & VarManager::ObjTypes::ReducedEventQvector) > 0);
      constexpr bool eventHasQvectorCentr = ((TEventFillMap & VarManager::ObjTypes::CollisionQvect) > 0);

      // TODO: FillPair functions need to provide a template argument to discriminate between cases when cov matrix is available or not
      VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
      if constexpr ((TPairType == pairTypeEE) || (TPairType == pairTypeMuMu)) { // call this just for ee or mumu pairs
        if constexpr (TTwoProngFitter == true) {
          VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, fPropToPCA);
        }
        if constexpr (eventHasQvector) {
          VarManager::FillPairVn<TEventFillMap, TPairType>(t1, t2);
        }
        if constexpr (eventHasQvectorCentr) {
          VarManager::FillPairVn<TEventFillMap, TPairType>(t1, t2);
        }
      }

      // TODO: provide the type of pair to the dilepton table (e.g. ee, mumu, emu...)
      dileptonFilterMap = twoTrackFilter;

      if constexpr (TPairType == pairTypeEE) {
        dielectronList(event, VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), dileptonFilterMap, dileptonMcDecision);
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrackCollInfo) > 0) {
          dielectronInfoList(t1.collisionId(), t1.trackId(), t2.trackId());
        }
      }
      if constexpr ((TPairType == pairTypeEE) && (TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelPID) > 0) {
        if (fConfigFlatTables.value) {
          dielectronAllList(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), dileptonFilterMap, dileptonMcDecision,
                            t1.pt(), t1.eta(), t1.phi(), t1.itsClusterMap(), t1.itsChi2NCl(), t1.tpcNClsCrossedRows(), t1.tpcNClsFound(), t1.tpcChi2NCl(), t1.dcaXY(), t1.dcaZ(), t1.tpcSignal(), t1.tpcNSigmaEl(), t1.tpcNSigmaPi(), t1.tpcNSigmaPr(), t1.beta(), t1.tofNSigmaEl(), t1.tofNSigmaPi(), t1.tofNSigmaPr(),
                            t2.pt(), t2.eta(), t2.phi(), t2.itsClusterMap(), t2.itsChi2NCl(), t2.tpcNClsCrossedRows(), t2.tpcNClsFound(), t2.tpcChi2NCl(), t2.dcaXY(), t2.dcaZ(), t2.tpcSignal(), t2.tpcNSigmaEl(), t2.tpcNSigmaPi(), t2.tpcNSigmaPr(), t2.beta(), t2.tofNSigmaEl(), t2.tofNSigmaPi(), t2.tofNSigmaPr(),
                            VarManager::fgValues[VarManager::kKFTrack0DCAxyz], VarManager::fgValues[VarManager::kKFTrack1DCAxyz], VarManager::fgValues[VarManager::kKFDCAxyzBetweenProngs], VarManager::fgValues[VarManager::kKFTrack0DCAxy], VarManager::fgValues[VarManager::kKFTrack1DCAxy], VarManager::fgValues[VarManager::kKFDCAxyBetweenProngs],
                            VarManager::fgValues[VarManager::kKFTrack0DeviationFromPV], VarManager::fgValues[VarManager::kKFTrack1DeviationFromPV], VarManager::fgValues[VarManager::kKFTrack0DeviationxyFromPV], VarManager::fgValues[VarManager::kKFTrack1DeviationxyFromPV],
                            VarManager::fgValues[VarManager::kKFMass], VarManager::fgValues[VarManager::kKFChi2OverNDFGeo], VarManager::fgValues[VarManager::kVertexingLxyz], VarManager::fgValues[VarManager::kVertexingLxyzOverErr], VarManager::fgValues[VarManager::kVertexingLxy], VarManager::fgValues[VarManager::kVertexingLxyOverErr], VarManager::fgValues[VarManager::kVertexingTauxy], VarManager::fgValues[VarManager::kVertexingTauxyErr], VarManager::fgValues[VarManager::kKFCosPA], VarManager::fgValues[VarManager::kKFJpsiDCAxyz], VarManager::fgValues[VarManager::kKFJpsiDCAxy],
                            VarManager::fgValues[VarManager::kKFPairDeviationFromPV], VarManager::fgValues[VarManager::kKFPairDeviationxyFromPV],
                            VarManager::fgValues[VarManager::kKFMassGeoTop], VarManager::fgValues[VarManager::kKFChi2OverNDFGeoTop]);
        }
      }
      if constexpr (TPairType == pairTypeMuMu) {
        dimuonList(event, VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), dileptonFilterMap, dileptonMcDecision);
      }
      if constexpr ((TPairType == pairTypeMuMu && ((TTrackFillMap & VarManager::ObjTypes::ReducedMuonCollInfo) > 0)) || (TPairType == pairTypeEE && ((TTrackFillMap & VarManager::ObjTypes::ReducedTrackCollInfo) > 0))) {
        dileptonInfoList(t1.collisionId(), event.posX(), event.posY(), event.posZ());
      }

      constexpr bool trackHasCov = ((TTrackFillMap & VarManager::ObjTypes::TrackCov) > 0 || (TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelCov) > 0);
      if constexpr ((TPairType == pairTypeEE) && trackHasCov && (TTwoProngFitter == true)) {
        dielectronExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingLz], VarManager::fgValues[VarManager::kVertexingLxy]);
      }
      if constexpr ((TPairType == pairTypeMuMu) && (TTwoProngFitter == true)) {
        // LOGP(info, "mu1 collId = {}, mu2 collId = {}", t1.collisionId(), t2.collisionId());
        dimuonExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingLz], VarManager::fgValues[VarManager::kVertexingLxy]);
        if (fConfigFlatTables.value) {
          dimuonAllList(event.posX(), event.posY(), event.posZ(), event.numContrib(),
                        event.selection_raw(), 0,
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
                        VarManager::fgValues[VarManager::kCORR2POI], VarManager::fgValues[VarManager::kCORR4POI], VarManager::fgValues[VarManager::kM01POI], VarManager::fgValues[VarManager::kM0111POI], VarManager::fgValues[VarManager::kMultDimuons],
                        VarManager::fgValues[VarManager::kVertexingPz], VarManager::fgValues[VarManager::kVertexingSV]);
        }
        auto collId = 0;
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedMuonCollInfo) > 0) {
          collId = t1.collisionId();
        }
        if constexpr (eventHasQvector == true || eventHasQvectorCentr == true) {
          dileptonFlowList(collId, VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kCentFT0C],
                           VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), isFirst,
                           VarManager::fgValues[VarManager::kU2Q2], VarManager::fgValues[VarManager::kR2SP_AB], VarManager::fgValues[VarManager::kR2SP_AC], VarManager::fgValues[VarManager::kR2SP_BC],
                           VarManager::fgValues[VarManager::kU3Q3], VarManager::fgValues[VarManager::kR3SP],
                           VarManager::fgValues[VarManager::kCos2DeltaPhi], VarManager::fgValues[VarManager::kR2EP_AB], VarManager::fgValues[VarManager::kR2EP_AC], VarManager::fgValues[VarManager::kR2EP_BC],
                           VarManager::fgValues[VarManager::kCos3DeltaPhi], VarManager::fgValues[VarManager::kR3EP],
                           VarManager::fgValues[VarManager::kCORR2POI], VarManager::fgValues[VarManager::kCORR4POI], VarManager::fgValues[VarManager::kM01POI], VarManager::fgValues[VarManager::kM0111POI],
                           VarManager::fgValues[VarManager::kCORR2REF], VarManager::fgValues[VarManager::kCORR4REF], VarManager::fgValues[VarManager::kM11REF], VarManager::fgValues[VarManager::kM1111REF],
                           VarManager::fgValues[VarManager::kMultDimuons], VarManager::fgValues[VarManager::kMultA]);
        }
        if (t1.sign() != t2.sign()) {
          isFirst = false;
        }
      }

      int iCut = 0;
      for (int icut = 0; icut < ncuts; icut++) {
        if (twoTrackFilter & (static_cast<uint32_t>(1) << icut)) {
          if (t1.sign() * t2.sign() < 0) {
            fHistMan->FillHistClass(histNames[iCut][0].Data(), VarManager::fgValues);
            if (fConfigAmbiguousHist && !(t1.isAmbiguous() || t2.isAmbiguous())) {
              fHistMan->FillHistClass(Form("%s_unambiguous", histNames[iCut][0].Data()), VarManager::fgValues);
            }
            if (useMiniTree.fConfigMiniTree) {
              float dileptonMass = VarManager::fgValues[VarManager::kMass];
              if (dileptonMass > useMiniTree.fConfigMiniTreeMinMass && dileptonMass < useMiniTree.fConfigMiniTreeMaxMass) {
                dileptonMiniTree(VarManager::fgValues[VarManager::kMass],
                                 VarManager::fgValues[VarManager::kPt],
                                 VarManager::fgValues[VarManager::kRap],
                                 VarManager::fgValues[VarManager::kCentFT0C],
                                 VarManager::fgValues[VarManager::kCos2DeltaPhi],
                                 VarManager::fgValues[VarManager::kPt1],
                                 VarManager::fgValues[VarManager::kEta1],
                                 VarManager::fgValues[VarManager::kPhi1],
                                 VarManager::fgValues[VarManager::kPt2],
                                 VarManager::fgValues[VarManager::kEta2],
                                 VarManager::fgValues[VarManager::kPhi2]);
              }
            }
          } else {
            if (t1.sign() > 0) {
              fHistMan->FillHistClass(histNames[iCut][1].Data(), VarManager::fgValues);
              if (fConfigAmbiguousHist && !(t1.isAmbiguous() || t2.isAmbiguous())) {
                fHistMan->FillHistClass(Form("%s_unambiguous", histNames[iCut][1].Data()), VarManager::fgValues);
              }
            } else {
              fHistMan->FillHistClass(histNames[iCut][2].Data(), VarManager::fgValues);
              if (fConfigAmbiguousHist && !(t1.isAmbiguous() || t2.isAmbiguous())) {
                fHistMan->FillHistClass(Form("%s_unambiguous", histNames[iCut][2].Data()), VarManager::fgValues);
              }
            }
          }
          iCut++;
          for (unsigned int iPairCut = 0; iPairCut < fPairCuts.size(); iPairCut++, iCut++) {
            AnalysisCompositeCut cut = fPairCuts.at(iPairCut);
            if (!(cut.IsSelected(VarManager::fgValues))) // apply pair cuts
              continue;
            if (t1.sign() * t2.sign() < 0) {
              fHistMan->FillHistClass(histNames[iCut][0].Data(), VarManager::fgValues);
              if (fConfigAmbiguousHist && !(t1.isAmbiguous() || t2.isAmbiguous())) {
                fHistMan->FillHistClass(Form("%s_unambiguous", histNames[iCut][0].Data()), VarManager::fgValues);
              }
            } else {
              if (t1.sign() > 0) {
                fHistMan->FillHistClass(histNames[iCut][1].Data(), VarManager::fgValues);
                if (fConfigAmbiguousHist && !(t1.isAmbiguous() || t2.isAmbiguous())) {
                  fHistMan->FillHistClass(Form("%s_unambiguous", histNames[iCut][1].Data()), VarManager::fgValues);
                }
              } else {
                fHistMan->FillHistClass(histNames[iCut][2].Data(), VarManager::fgValues);
                if (fConfigAmbiguousHist && !(t1.isAmbiguous() || t2.isAmbiguous())) {
                  fHistMan->FillHistClass(Form("%s_unambiguous", histNames[iCut][2].Data()), VarManager::fgValues);
                }
              }
            }
          } // end loop (pair cuts)
        } else { // end if (filter bits)
          iCut = iCut + 1 + fPairCuts.size();
        }
      } // end loop (cuts)
    } // end loop over pairs
  }

  void processDecayToEESkimmed(soa::Filtered<MyEventsSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMap>(event, tracks, tracks);
  }
  void processDecayToEESkimmedNoTwoProngFitter(soa::Filtered<MyEventsSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<false, VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMap>(event, tracks, tracks);
  }
  void processDecayToEESkimmedWithCov(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelectedWithCov> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCov>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, tracks, tracks);
  }
  void processDecayToEESkimmedWithCovNoTwoProngFitter(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelectedWithCov> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCov>(event, VarManager::fgValues);
    runSameEventPairing<false, VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, tracks, tracks);
  }
  void processDecayToEEVertexingSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelectedWithCov> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, tracks, tracks);
  }
  void processDecayToEEPrefilterSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelectedWithPrefilter> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMap>(event, tracks, tracks);
  }
  void processDecayToEEPrefilterSkimmedNoTwoProngFitter(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelectedWithPrefilter> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<false, VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMap>(event, tracks, tracks);
  }
  void processDecayToEESkimmedWithColl(soa::Filtered<MyEventsSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelectedWithColl> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMapWithColl>(event, tracks, tracks);
  }
  void processDecayToEESkimmedWithCollNoTwoProngFitter(soa::Filtered<MyEventsSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelectedWithColl> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<false, VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMapWithColl>(event, tracks, tracks);
  }
  void processDecayToMuMuSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMap, gkMuonFillMap>(event, muons, muons);
  }
  void processDecayToMuMuSkimmedWithMult(soa::Filtered<MyEventsVtxCovSelectedMultExtra>::iterator const& event, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithMultExtra>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithMultExtra, gkMuonFillMap>(event, muons, muons);
  }
  void processDecayToMuMuVertexingSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyMuonTracksSelectedWithCov> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithCov, gkMuonFillMapWithCov>(event, muons, muons);
  }
  void processDecayToMuMuSkimmedWithColl(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyMuonTracksSelectedWithColl> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMap, gkMuonFillMapWithColl>(event, muons, muons);
  }
  void processVnDecayToEESkimmed(soa::Filtered<MyEventsVtxCovSelectedQvector>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCovQvector>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithCovQvector, gkTrackFillMap>(event, tracks, tracks);
  }
  void processVnDecayToMuMuSkimmed(soa::Filtered<MyEventsVtxCovSelectedQvector>::iterator const& event, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCovQvector>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithCovQvector, gkMuonFillMap>(event, muons, muons);
  }
  void processVnDecayToMuMuSkimmedWithWeights(soa::Filtered<MyEventsVtxCovSelectedQvectorExtraWithRefFlow>::iterator const& event, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCovQvectorExtraWithRefFlow>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithCovQvectorExtraWithRefFlow, gkMuonFillMap>(event, muons, muons);
  }

  void processVnDecayToMuMuSkimmedWithWeightsAndColl(soa::Filtered<MyEventsVtxCovSelectedQvectorExtraWithRefFlow>::iterator const& event, soa::Filtered<MyMuonTracksSelectedWithColl> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCovQvectorExtraWithRefFlow>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithCovQvectorExtraWithRefFlow, gkMuonFillMapWithColl>(event, muons, muons);
  }

  void processVnCentrDecayToMuMuSkimmed(soa::Filtered<MyEventsVtxCovSelectedQvectorCentr>::iterator const& event, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCovQvectorCentr>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithCovQvectorCentr, gkMuonFillMap>(event, muons, muons);
  }
  void processElectronMuonSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kElectronMuon, gkEventFillMap, gkTrackFillMap>(event, tracks, muons);
  }
  void processDecayToPiPiSkimmed(soa::Filtered<MyEventsSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToPiPi, gkEventFillMap, gkTrackFillMap>(event, tracks, tracks);
  }
  void processAllSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMap>(event, tracks, tracks);
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMap, gkMuonFillMap>(event, muons, muons);
    runSameEventPairing<true, VarManager::kElectronMuon, gkEventFillMap, gkTrackFillMap>(event, tracks, muons);
  }
  // TODO: dummy function for the case when no process function is enabled
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmed, "Run electron-electron pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmedNoTwoProngFitter, "Run electron-electron pairing, with skimmed tracks but no two prong fitter", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmedWithCov, "Run electron-electron pairing, with skimmed covariant tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmedWithCovNoTwoProngFitter, "Run electron-electron pairing, with skimmed covariant tracks but no two prong fitter", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEEVertexingSkimmed, "Run electron-electron pairing and vertexing, with skimmed electrons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEEPrefilterSkimmed, "Run electron-electron pairing, with skimmed tracks and prefilter from AnalysisPrefilterSelection", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEEPrefilterSkimmedNoTwoProngFitter, "Run electron-electron pairing, with skimmed tracks and prefilter from AnalysisPrefilterSelection but no two prong fitter", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmedWithColl, "Run electron-electron pairing, with skimmed tracks and with collision information", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmedWithCollNoTwoProngFitter, "Run electron-electron pairing, with skimmed tracks and with collision information but no two prong fitter", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToMuMuSkimmed, "Run muon-muon pairing, with skimmed muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToMuMuSkimmedWithMult, "Run muon-muon pairing, with skimmed muons and multiplicity", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToMuMuVertexingSkimmed, "Run muon-muon pairing and vertexing, with skimmed muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToMuMuSkimmedWithColl, "Run muon-muon pairing keeping the info of AO2D collision, with skimmed muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processVnDecayToEESkimmed, "Run electron-electron pairing, with skimmed tracks for vn", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processVnDecayToMuMuSkimmed, "Run muon-muon pairing, with skimmed tracks for vn", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processVnDecayToMuMuSkimmedWithWeights, "Run muon-muon pairing, with skimmed tracks for vn using weights for Q-vectors", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processVnDecayToMuMuSkimmedWithWeightsAndColl, "Run muon-muon pairing, with skimmed tracks for vn using weights for Q-vectors with collisionId provided", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processVnCentrDecayToMuMuSkimmed, "Run muon-muon pairing, with skimmed tracks for vn from central framework", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processElectronMuonSkimmed, "Run electron-muon pairing, with skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToPiPiSkimmed, "Run pion-pion pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processAllSkimmed, "Run all types of pairing, with skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", false);
};

struct AnalysisDileptonHadron {
  //
  // This task combines dilepton candidates with a track and could be used for example
  //  in analyses with the dilepton as one of the decay products of a higher mass resonance (e.g. B -> Jpsi + K)
  //    or in dilepton + hadron correlations, etc.
  //
  //  The barrel and muon track filtering tasks can produce multiple parallel decisions, which are used to produce
  //   dileptons which inherit the logical intersection of the track level decisions (see the AnalysisSameEventPairing task).
  //  This can be used also in the dilepton-hadron correlation analysis. However, in this model of the task, we use all the dileptons produced in the
  //    lepton pairing task to combine them with the hadrons selected by the barrel track selection.
  //  To be modified/adapted if new requirements appear
  float mMagField = 0.0;
  o2::parameters::GRPMagField* grpmag = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  OutputObj<THashList> fOutputList{"output"};
  // TODO: For now this is only used to determine the position in the filter bit map for the hadron cut
  Configurable<string> fConfigTrackCuts{"cfgLeptonCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  // comment: add list of subgroups (must define subgroups under )
  Configurable<std::string> fConfigAddDileptonHadHistogram{"cfgAddDileptonHadHistogram", "", "Comma separated list of histograms"};
  Configurable<int> fConfigMixingDepth{"cfgMixingDepth", 5, "Event mixing pool depth"};
  Configurable<float> fConfigDileptonLowMass{"cfgDileptonLowMass", 0.0, "Low mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonHighMass{"cfgDileptonHighMass", 5.0, "High mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonpTCut{"cfgDileptonpTCut", 0.0, "pT cut for dileptons used in the triplet vertexing"};
  Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
  Configurable<bool> fUseRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Produces<aod::BmesonCandidates> BmesonsTable;

  Filter eventFilter = aod::dqanalysisflags::isEventSelected == 1;
  Filter dileptonFilter = aod::reducedpair::pt > fConfigDileptonpTCut&& aod::reducedpair::mass > fConfigDileptonLowMass&& aod::reducedpair::mass < fConfigDileptonHighMass&& aod::reducedpair::sign == 0;
  Filter filterBarrelTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair; // fill map

  // use two values array to avoid mixing up the quantities
  float* fValuesDilepton;
  float* fValuesHadron;
  HistogramManager* fHistMan;

  // NOTE: the barrel track filter is shared between the filters for dilepton electron candidates (first n-bits)
  //       and the associated hadrons (n+1 bit) --> see the barrel track selection task
  //      The current condition should be replaced when bitwise operators will become available in Filter expressions
  int fNHadronCutBit;

  NoBinningPolicy<aod::dqanalysisflags::MixingHash> hashBin;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    fCurrentRun = 0;
    fValuesDilepton = new float[VarManager::kNVars];
    fValuesHadron = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // TODO: Create separate histogram directories for each selection used in the creation of the dileptons
    // TODO: Implement possibly multiple selections for the associated track ?
    if (context.mOptions.get<bool>("processSkimmed")) {
      DefineHistograms(fHistMan, "DileptonsSelected;DileptonHadronInvMass;DileptonHadronCorrelationSE", fConfigAddDileptonHadHistogram); // define all histograms
    }
    if (context.mOptions.get<bool>("processMixedEvent")) {
      DefineHistograms(fHistMan, "DileptonHadronInvMassME;DileptonHadronCorrelationME", fConfigAddDileptonHadHistogram); // define all histograms
    }

    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    TString configCutNamesStr = fConfigTrackCuts.value;
    if (!configCutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(configCutNamesStr.Tokenize(","));
      fNHadronCutBit = objArray->GetEntries() - 1;
    } else {
      fNHadronCutBit = 0;
    }
  }

  // Template function to run pair - hadron combinations
  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runDileptonHadron(TEvent const& event, TTracks const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {

    // set up KF or DCAfitter
    if (fCurrentRun != event.runNumber()) { // start: runNumber
      if (fUseRemoteField.value) {
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, event.timestamp());
        if (grpmag != nullptr) {
          mMagField = grpmag->getNominalL3Field();
        } else {
          LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", event.timestamp());
        }
        if (fConfigUseKFVertexing) {
          VarManager::SetupThreeProngKFParticle(mMagField);
        } else {
          VarManager::SetupThreeProngDCAFitter(mMagField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, false); // TODO: get these parameters from Configurables
        }
      } else {
        if (fConfigUseKFVertexing) {
          VarManager::SetupThreeProngKFParticle(fConfigMagField.value);
        } else {
          VarManager::SetupThreeProngDCAFitter(fConfigMagField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, false); // TODO: get these parameters from Configurables
        }
      }
      fCurrentRun = event.runNumber();
    } // end: runNumber

    VarManager::ResetValues(0, VarManager::kNVars, fValuesHadron);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    VarManager::FillEvent<TEventFillMap>(event, fValuesHadron);
    VarManager::FillEvent<TEventFillMap>(event, fValuesDilepton);

    // Set the global index offset to find the proper lepton
    // TO DO: remove it once the issue with lepton index is solved
    int indexOffset = -999;
    std::vector<int> trackGlobalIndexes;

    if (dileptons.size() > 0) {
      for (auto track : tracks) {
        trackGlobalIndexes.push_back(track.globalIndex());
        // std::cout << track.index() << " " << track.globalIndex() << std::endl;
      }
    }
    // loop once over dileptons for QA purposes
    for (auto dilepton : dileptons) {
      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesDilepton);
      fHistMan->FillHistClass("DileptonsSelected", fValuesDilepton);

      // get the index of the electron legs
      int indexLepton1 = dilepton.index0Id();
      int indexLepton2 = dilepton.index1Id();

      if (indexOffset == -999) {
        indexOffset = trackGlobalIndexes.at(0);
      }
      trackGlobalIndexes.clear();

      // get full track info of tracks based on the index
      auto lepton1 = tracks.iteratorAt(indexLepton1 - indexOffset);
      auto lepton2 = tracks.iteratorAt(indexLepton2 - indexOffset);

      // Check that the dilepton has zero charge
      if (dilepton.sign() != 0) {
        continue;
      }

      // Check that the leptons are opposite sign
      if (lepton1.sign() * lepton2.sign() > 0) {
        continue;
      }

      // loop over hadrons
      for (auto& hadron : tracks) {
        if (!(static_cast<uint32_t>(hadron.isBarrelSelected()) & (static_cast<uint32_t>(1) << fNHadronCutBit))) {
          continue;
        }

        // if the hadron is either of the electron legs, continue
        int index = hadron.globalIndex();
        if (index == indexLepton1 || index == indexLepton2) {
          continue;
        }

        VarManager::FillDileptonHadron(dilepton, hadron, fValuesHadron);
        // VarManager::FillDileptonTrackVertexing<TCandidateType, TEventFillMap, TTrackFillMap>(event, lepton1, lepton2, hadron, fValuesHadron);
        fHistMan->FillHistClass("DileptonHadronInvMass", fValuesHadron);
        fHistMan->FillHistClass("DileptonHadronCorrelationSE", fValuesHadron);
        // table to be written out for ML analysis
        BmesonsTable(fValuesHadron[VarManager::kPairMass], fValuesHadron[VarManager::kPairPt], fValuesHadron[VarManager::kVertexingLxy], fValuesHadron[VarManager::kVertexingLxyz], fValuesHadron[VarManager::kVertexingLz], fValuesHadron[VarManager::kVertexingTauxy], fValuesHadron[VarManager::kVertexingTauz], fValuesHadron[VarManager::kCosPointingAngle], fValuesHadron[VarManager::kVertexingChi2PCA]);
      }
    }
  }

  void processSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyBarrelTracksSelectedWithCov const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {
    runDileptonHadron<VarManager::kBtoJpsiEEK, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, tracks, dileptons);
  }

  Preslice<soa::Filtered<MyDielectronCandidates>> perEventPairs = aod::reducedpair::reducedeventId;
  Preslice<soa::Filtered<MyBarrelTracksSelected>> perEventTracks = aod::reducedtrack::reducedeventId;

  void processMixedEvent(soa::Filtered<MyEventsHashSelected>& events, soa::Filtered<MyDielectronCandidates> const& dileptons, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    events.bindExternalIndices(&dileptons);
    events.bindExternalIndices(&tracks);

    for (auto& [event1, event2] : selfCombinations(hashBin, fConfigMixingDepth.value, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<gkEventFillMap>(event1, VarManager::fgValues);

      auto evDileptons = dileptons.sliceBy(perEventPairs, event1.globalIndex());
      evDileptons.bindExternalIndices(&events);

      auto evTracks = tracks.sliceBy(perEventTracks, event2.globalIndex());
      evTracks.bindExternalIndices(&events);

      for (auto dilepton : evDileptons) {
        for (auto& track : evTracks) {

          if (!(static_cast<uint32_t>(track.isBarrelSelected()) & (static_cast<uint32_t>(1) << fNHadronCutBit))) {
            continue;
          }

          VarManager::FillDileptonHadron(dilepton, track, VarManager::fgValues);
          fHistMan->FillHistClass("DileptonHadronInvMassME", VarManager::fgValues);
          fHistMan->FillHistClass("DileptonHadronCorrelationME", VarManager::fgValues);
        } // end for (track)
      } // end for (dilepton)

    } // end event loop
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisDileptonHadron, processSkimmed, "Run dilepton-hadron pairing, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonHadron, processMixedEvent, "Run dilepton-hadron mixed event pairing", false);
  PROCESS_SWITCH(AnalysisDileptonHadron, processDummy, "Dummy function", false);
};

struct AnalysisDileptonTrackTrack {
  OutputObj<THashList> fOutputList{"output"};

  Configurable<std::string> fConfigTrackCut1{"cfgTrackCut1", "pionPIDCut1", "track1 cut"}; // used for select the tracks from SelectedTracks
  Configurable<std::string> fConfigTrackCut2{"cfgTrackCut2", "pionPIDCut2", "track2 cut"}; // used for select the tracks from SelectedTracks
  Configurable<std::string> fConfigDileptonCut{"cfgDiLeptonCut", "pairJpsi2", "Dilepton cut"};
  Configurable<std::string> fConfigQuadrupletCuts{"cfgQuadrupletCuts", "pairX3872Cut1", "Comma separated list of Dilepton-Track-Track cut"};
  Configurable<std::string> fConfigAddDileptonHistogram{"cfgAddDileptonHistogram", "barrel", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddQuadrupletHistogram{"cfgAddQuadrupletHistogram", "xtojpsipipi", "Comma separated list of histograms"};

  Produces<aod::DileptonTrackTrackCandidates> DileptonTrackTrackTable;

  Filter eventFilter = aod::dqanalysisflags::isEventSelected == 1;
  Filter dileptonFilter = aod::reducedpair::sign == 0;
  Filter filterBarrelTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair; // fill map

  // use some values array to avoid mixing up the quantities
  float* fValuesDitrack;
  float* fValuesQuadruplet;
  HistogramManager* fHistMan;

  // cut name setting
  TString fTrackCutName1;
  TString fTrackCutName2;
  bool fIsSameTrackCut = false;
  AnalysisCompositeCut fDileptonCut;
  std::vector<TString> fQuadrupletCutNames;
  std::vector<AnalysisCompositeCut> fQuadrupletCuts;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    fValuesQuadruplet = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // define cuts
    fTrackCutName1 = fConfigTrackCut1.value;
    fTrackCutName2 = fConfigTrackCut2.value;
    if (fTrackCutName1 == fTrackCutName2) {
      fIsSameTrackCut = true;
    }
    TString configDileptonCutNamesStr = fConfigDileptonCut.value;
    fDileptonCut = *dqcuts::GetCompositeCut(configDileptonCutNamesStr.Data());
    TString configQuadruletCutNamesStr = fConfigQuadrupletCuts.value;
    std::unique_ptr<TObjArray> objArray(configQuadruletCutNamesStr.Tokenize(","));
    for (Int_t icut = 0; icut < objArray->GetEntries(); ++icut) {
      TString cutName = objArray->At(icut)->GetName();
      fQuadrupletCutNames.push_back(cutName);
      fQuadrupletCuts.push_back(*dqcuts::GetCompositeCut(cutName.Data()));
    }

    if (!context.mOptions.get<bool>("processDummy")) {
      DefineHistograms(fHistMan, Form("Dileptons_%s", configDileptonCutNamesStr.Data()), fConfigAddDileptonHistogram);
      if (!configQuadruletCutNamesStr.IsNull()) {
        for (std::size_t icut = 0; icut < fQuadrupletCutNames.size(); ++icut) {
          if (fIsSameTrackCut) {
            DefineHistograms(fHistMan, Form("QuadrupletSEPM_%s", fQuadrupletCutNames[icut].Data()), fConfigAddQuadrupletHistogram);
          } else {
            DefineHistograms(fHistMan, Form("QuadrupletSEPM_%s", fQuadrupletCutNames[icut].Data()), fConfigAddQuadrupletHistogram);
            DefineHistograms(fHistMan, Form("QuadrupletSEMP_%s", fQuadrupletCutNames[icut].Data()), fConfigAddQuadrupletHistogram);
          }
          DefineHistograms(fHistMan, Form("QuadrupletSEPP_%s", fQuadrupletCutNames[icut].Data()), fConfigAddQuadrupletHistogram);
          DefineHistograms(fHistMan, Form("QuadrupletSEMM_%s", fQuadrupletCutNames[icut].Data()), fConfigAddQuadrupletHistogram);
        }
      }
    }

    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }
  // Template function to run pair - track - track combinations
  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runDileptonTrackTrack(TEvent const& event, TTracks const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesQuadruplet);
    VarManager::FillEvent<TEventFillMap>(event, fValuesQuadruplet);

    // LOGF(info, "Number of dileptons: %d", dileptons.size());

    // loop over dileptons
    for (auto dilepton : dileptons) {
      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesQuadruplet);

      // apply the dilepton cut
      if (!fDileptonCut.IsSelected(fValuesQuadruplet))
        continue;

      fHistMan->FillHistClass(Form("Dileptons_%s", fDileptonCut.GetName()), fValuesQuadruplet);

      // get the index of the electron legs
      int indexLepton1 = dilepton.index0Id();
      int indexLepton2 = dilepton.index1Id();

      // loop over hadrons pairs
      for (auto& [t1, t2] : combinations(tracks, tracks)) {
        // avoid self-combinations
        if (t1.globalIndex() == indexLepton1 || t1.globalIndex() == indexLepton2 || t2.globalIndex() == indexLepton1 || t2.globalIndex() == indexLepton2) {
          continue;
        }

        if (fIsSameTrackCut) {
          if (!(static_cast<uint32_t>(t1.isBarrelSelected()) & (static_cast<uint32_t>(1) << 1)) || !(static_cast<uint32_t>(t2.isBarrelSelected()) & (static_cast<uint32_t>(1) << 1))) {
            continue;
          }
        } else {
          if (!(static_cast<uint32_t>(t1.isBarrelSelected()) & (static_cast<uint32_t>(1) << 1)) || !(static_cast<uint32_t>(t2.isBarrelSelected()) & (static_cast<uint32_t>(1) << 2))) {
            continue;
          }
        }

        // fill variables
        VarManager::FillDileptonTrackTrack<TCandidateType>(dilepton, t1, t2, fValuesQuadruplet);

        int iCut = 0;
        uint32_t CutDecision = 0;
        for (auto cutname = fQuadrupletCutNames.begin(); cutname != fQuadrupletCutNames.end(); cutname++, iCut++) {
          // apply dilepton-track-track cut
          if (fQuadrupletCuts[iCut].IsSelected(fValuesQuadruplet)) {
            CutDecision |= (1 << iCut);
            if (fIsSameTrackCut) {
              if (t1.sign() * t2.sign() < 0) {
                fHistMan->FillHistClass(Form("QuadrupletSEPM_%s", fQuadrupletCutNames[iCut].Data()), fValuesQuadruplet);
              }
            } else {
              if ((t1.sign() < 0) && (t2.sign() > 0)) {
                fHistMan->FillHistClass(Form("QuadrupletSEMP_%s", fQuadrupletCutNames[iCut].Data()), fValuesQuadruplet);
              } else if ((t1.sign() > 0) && (t2.sign() < 0)) {
                fHistMan->FillHistClass(Form("QuadrupletSEPM_%s", fQuadrupletCutNames[iCut].Data()), fValuesQuadruplet);
              }
            }
            if ((t1.sign() > 0) && (t2.sign() > 0)) {
              fHistMan->FillHistClass(Form("QuadrupletSEPP_%s", fQuadrupletCutNames[iCut].Data()), fValuesQuadruplet);
            } else if ((t1.sign() < 0) && (t2.sign() < 0)) {
              fHistMan->FillHistClass(Form("QuadrupletSEMM_%s", fQuadrupletCutNames[iCut].Data()), fValuesQuadruplet);
            }
          }
        } // loop over dilepton-track-track cuts

        // table to be written out for ML analysis
        if (!CutDecision)
          continue;
        DileptonTrackTrackTable(fValuesQuadruplet[VarManager::kQuadMass], fValuesQuadruplet[VarManager::kQuadPt], fValuesQuadruplet[VarManager::kQuadEta], fValuesQuadruplet[VarManager::kQuadPhi], fValuesQuadruplet[VarManager::kRap],
                                fValuesQuadruplet[VarManager::kQ], fValuesQuadruplet[VarManager::kDeltaR1], fValuesQuadruplet[VarManager::kDeltaR2], fValuesQuadruplet[VarManager::kDeltaR],
                                dilepton.mass(), dilepton.pt(), dilepton.eta(), dilepton.phi(), dilepton.sign(),
                                fValuesQuadruplet[VarManager::kDitrackMass], fValuesQuadruplet[VarManager::kDitrackPt], t1.pt(), t2.pt(), t1.eta(), t2.eta(), t1.phi(), t2.phi(), t1.sign(), t2.sign());
      } // end loop over track-track pairs
    } // end loop over dileptons
  }

  void processJpsiPiPi(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyBarrelTracksSelectedWithCov const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {
    runDileptonTrackTrack<VarManager::kXtoJpsiPiPi, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, tracks, dileptons);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisDileptonTrackTrack, processJpsiPiPi, "Run dilepton-dihadron pairing to study X(3872), using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonTrackTrack, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
    adaptAnalysisTask<AnalysisMuonSelection>(cfgc),
    adaptAnalysisTask<AnalysisPrefilterSelection>(cfgc),
    adaptAnalysisTask<AnalysisEventMixing>(cfgc),
    adaptAnalysisTask<AnalysisSameEventPairing>(cfgc),
    adaptAnalysisTask<AnalysisDileptonHadron>(cfgc),
    adaptAnalysisTask<AnalysisDileptonTrackTrack>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar)
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

    TString histName = configVar.value;
    // NOTE: The level of detail for histogramming can be controlled via configurables
    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", histName);
    }

    if (classStr.Contains("Track") && !classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
      }
      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
      }
    }

    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }

    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "barrel");
    }

    if (classStr.Contains("HadronsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
    }

    if (classStr.Contains("DileptonHadronInvMass")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-mass");
    }

    if (classStr.Contains("DileptonHadronCorrelation")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-array-correlation");
    }

    if (classStr.Contains("DitrackSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }

    if (classStr.Contains("Quadruplet")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-dihadron", histName);
    }

  } // end loop over histogram classes
}
