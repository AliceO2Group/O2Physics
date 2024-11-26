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

#include <map>
#include <string>
#include <memory>
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
DECLARE_SOA_COLUMN(IsBarrelSelectedPrefilter, isBarrelSelectedPrefilter, int);
DECLARE_SOA_COLUMN(IsPrefilterVetoed, isPrefilterVetoed, int);
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASH", dqanalysisflags::MixingHash);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTS", dqanalysisflags::IsBarrelSelected, dqanalysisflags::IsBarrelSelectedPrefilter);
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTER", dqanalysisflags::IsPrefilterVetoed);
} // namespace o2::aod

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsHashSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts>;
using MyEventsVtxCovSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvector>;
// using MyEventsVtxCovSelectedQvectorCentr = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvectorCentr>;
using MyEventsQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector>;
using MyEventsQvectorMultExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;
using MyEventsQvectorCentr = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvectorCentr>;
using MyEventsQvectorCentrMultExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvectorCentr, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;
using MyEventsHashSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvector>;
// using MyEventsHashSelectedQvectorCentr = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvectorCentr>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;
using MyBarrelTracksSelectedWithPrefilter = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::Prefilter>;
using MyBarrelTracksSelectedWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;
using MyBarrelTracksSelectedWithColl = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::ReducedTracksBarrelInfo>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;
constexpr static uint32_t gkEventFillMapWithQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkEventFillMapWithQvectorMultExtra = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector | VarManager::ObjTypes::ReducedEventMultExtra;
constexpr static uint32_t gkEventFillMapWithQvectorCentr = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::CollisionQvect;
constexpr static uint32_t gkEventFillMapWithQvectorCentrMultExtra = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::CollisionQvect | VarManager::ObjTypes::ReducedEventMultExtra;
constexpr static uint32_t gkEventFillMapWithCovQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector;
// constexpr static uint32_t gkEventFillMapWithCovQvectorCentr = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::CollisionQvect;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithColl = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID | VarManager::ObjTypes::ReducedTrackCollInfo;

constexpr static int pairTypeEE = VarManager::kDecayToEE;

// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar); // defines histograms for all tasks

struct AnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};
  // TODO: Provide the mixing variables and binning directly via configurables (e.g. vectors of float)
  Configurable<string> fConfigMixingVariables{"cfgMixingVars", "", "Mixing configs separated by a comma, default no mixing"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<int> fConfigITSROFrameStartBorderMargin{"ITSROFrameStartBorderMargin", -1, "Number of bcs at the start of ITS RO Frame border. Take from CCDB if -1"};
  Configurable<int> fConfigITSROFrameEndBorderMargin{"ITSROFrameEndBorderMargin", -1, "Number of bcs at the end of ITS RO Frame border. Take from CCDB if -1"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  AnalysisCompositeCut* fEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<float> cfgCentFT0CMin{"cfgCentralityMin", -1000000000.f, "min. centrality"};
    Configurable<float> cfgCentFT0CMax{"cfgCentralityMax", 1000000000.f, "max. centrality"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -1000000000, "min. occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
  } eventcuts;

  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;
  int fLastRun;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  void init(o2::framework::InitContext&)
  {
    fEventCut = new AnalysisCompositeCut(true);
    fEventCut->AddCut(GetEventCut());
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
    if (fEventCut->IsSelected(VarManager::fgValues)) {
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
      }
      eventSel(1);
    } else {
      eventSel(0);
    }

    if (fMixHandler != nullptr) {
      int hh = fMixHandler->FindEventCategory(VarManager::fgValues);
      hash(hh);
    }
  }

  AnalysisCut* GetEventCut()
  {
    AnalysisCut* cut = new AnalysisCut("eventcut", "eventcut");
    cut->AddCut(VarManager::kVtxZ, -eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    if (eventcuts.cfgRequireSel8)
      cut->AddCut(VarManager::kIsSel8, 0.5, 1.5);
    if (eventcuts.cfgRequireNoTFB)
      cut->AddCut(VarManager::kIsNoTFBorder, 0.5, 1.5);
    if (eventcuts.cfgRequireNoITSROFB)
      cut->AddCut(VarManager::kIsNoITSROFBorder, 0.5, 1.5);
    if (eventcuts.cfgRequireNoSameBunchPileup)
      cut->AddCut(VarManager::kIsNoSameBunch, 0.5, 1.5);
    if (eventcuts.cfgRequireGoodZvtxFT0vsPV)
      cut->AddCut(VarManager::kIsGoodZvtxFT0vsPV, 0.5, 1.5);
    cut->AddCut(VarManager::kCentFT0C, eventcuts.cfgCentFT0CMin, eventcuts.cfgCentFT0CMax);
    cut->AddCut(VarManager::kTrackOccupancyInTimeRange, eventcuts.cfgTrackOccupancyMin, eventcuts.cfgTrackOccupancyMax);
    return cut;
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
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedQVector, "Run event selection on DQ skimmed events with Q vector from GFW", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedQVectorCentr, "Run event selection on DQ skimmed events with Q vector from CFW", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedQVectorMultExtra, "Run event selection on DQ skimmed events with Q vector from GFW and MultPV", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedQVectorCentrMultExtra, "Run event selection on DQ skimmed events with Q vector from CFW and MultPV", false);
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
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<bool> fConfigQAIfSelEvt{"cfgQAIfSelEvt", true, "If true, fill QA only for selected events"};
  Configurable<string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
  Configurable<int> fConfigPrefilterCutId{"cfgPrefilterCutId", 32, "Id of the Prefilter track cut (starting at 0)"}; // In order to create another column prefilter (should be temporary before improving cut selection in configurables, then displaced to AnalysisPrefilterSelection)
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<string> fConfigCcdbPathTPC{"ccdb-path-tpc", "Users/z/zhxiong/TPCPID/PostCalib", "base path to the ccdb object"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<bool> fConfigComputeTPCpostCalib{"cfgTPCpostCalib", false, "If true, compute TPC post-calibrated n-sigmas"};
  Configurable<std::string> fConfigRunPeriods{"cfgRunPeriods", "LHC22f", "run periods for used data"};
  Configurable<bool> fConfigDummyRunlist{"cfgDummyRunlist", false, "If true, use dummy runlist"};
  Configurable<int> fConfigInitRunNumber{"cfgInitRunNumber", 543215, "Initial run number used in run by run checks"};
  Configurable<std::size_t> fConfigNbTrackCut{"cfgNbTrackCut", 1, "Number of cuts including prefilter cut, need to be below 30"};

  std::vector<AnalysisCompositeCut> fTrackCuts;
  struct : ConfigurableGroup {
    std::string prefix = "trackingcut_group";
    Configurable<std::vector<float>> cfgEtaMax{"cfgEtaMax", {0.8f}, "Comma separated list of eta ranges"};
    Configurable<std::vector<float>> cfgPtMin{"cfgPtMin", {0.2f}, "Comma separated list of pt min"};
    Configurable<std::vector<float>> cfgPtMax{"cfgPtMax", {20.f}, "Comma separated list of pt max"};
    Configurable<std::vector<float>> cfgDCAxyMax{"cfgDCAxyMax", {3.f}, "Comma separated list of dcaxy max"};
    Configurable<std::vector<float>> cfgDCAzMax{"cfgDCAzMax", {1.f}, "Comma separated list of dcaz max"};
    Configurable<std::vector<int>> cfgIsSPDfirst{"cfgIsSPDfirst", {1}, "Comma separated list of if one requires one hit in the first ITS layer"};
    Configurable<std::vector<int>> cfgIsSPDany{"cfgIsSPDany", {-1}, "Comma separated list of if one requires one hit in the first two ITS layers"};
    Configurable<std::vector<int>> cfgIsITSibAny{"cfgIsITSibAny", {-1}, "Comma separated list of if one requires one hit in the first three ITS layers"};
    Configurable<std::vector<float>> cfgITSchi2Max{"cfgITSchi2Max", {5.f}, "Comma separated list of its chi2 max"};
    Configurable<std::vector<float>> cfgITSnclsMin{"cfgITSnclsMin", {4.5f}, "Comma separated list of min number of ITS clusters"};
    Configurable<std::vector<float>> cfgITSnclsMax{"cfgITSnclsMax", {7.5f}, "Comma separated list of max number of ITS clusters"};
    Configurable<std::vector<float>> cfgTPCchi2Max{"cfgTPCchi2Max", {4.f}, "Comma separated list of tpc chi2 max"};
    Configurable<std::vector<float>> cfgTPCnclsMin{"cfgTPCnclsMin", {90.f}, "Comma separated list of min number of TPC clusters"};
    Configurable<std::vector<float>> cfgTPCnclsMax{"cfgTPCnclsMax", {170.f}, "Comma separated list of max number of TPC clusters"};
    Configurable<std::vector<float>> cfgTPCnclsCRMin{"cfgTPCnclsCRMin", {80.f}, "Comma separated list of min number of TPC crossed rows"};
    Configurable<std::vector<float>> cfgTPCnclsCRMax{"cfgTPCnclsCRMax", {170.f}, "Comma separated list of max number of TPC crossed rows"};
    Configurable<std::vector<int>> cfgIsDalitzLeg{"cfgIsDalitzLeg", {-1}, "Comma separated list of if one requires hit for prefilter done during skimming, should be between 1 and 8"};
  } trackcuts;

  struct : ConfigurableGroup {
    std::string prefix = "pidcut_group";
    Configurable<std::vector<int>> cfgPIDmode{"cfgPIDmode", {1}, "List of PID mode: 1 TPChadrrejection, 2 TOFreq, 3 OR between both"};
    Configurable<std::vector<int>> cfgRejBadTOF{"cfgRejBadTOF", {1}, "List of reject bad TOF: 1 yes, 0 no"};
    Configurable<std::vector<float>> cfgTPCNSigmaElMin{"cfgTPCNSigmaElMin", {-1.f}, "Comma separated list of min TPC nsigma e for inclusion"};
    Configurable<std::vector<float>> cfgTPCNSigmaElMax{"cfgTPCNSigmaElMax", {3.f}, "Comma separated list of max TPC nsigma e for inclusion"};
    Configurable<std::vector<float>> cfgTPCNSigmaPiMin{"cfgTPCNSigmaPiMin", {-3.f}, "Comma separated list of min TPC nsigma pion for exclusion"};
    Configurable<std::vector<float>> cfgTPCNSigmaPiMax{"cfgTPCNSigmaPiMax", {4.f}, "Comma separated list of max TPC nsigma pion for exclusion"};
    Configurable<std::vector<float>> cfgTPCNSigmaPrMin{"cfgTPCNSigmaPrMin", {-3.f}, "Comma separated list of min TPC nsigma proton for exclusion"};
    Configurable<std::vector<float>> cfgTPCNSigmaPrMax{"cfgTPCNSigmaPrMax", {4.f}, "Comma separated list of max TPC nsigma proton for exclusion"};
    Configurable<std::vector<float>> cfgTPCNSigmaKaMin{"cfgTPCNSigmaKaMin", {-3.f}, "Comma separated list of min TPC nsigma kaon for exclusion"};
    Configurable<std::vector<float>> cfgTPCNSigmaKaMax{"cfgTPCNSigmaKaMax", {4.f}, "Comma separated list of max TPC nsigma kaon for exclusion"};
    Configurable<std::vector<float>> cfgTPCNSigmaMuMin{"cfgTPCNSigmaMuMin", {0.f}, "Comma separated list of min TPC nsigma muon for exclusion"};
    Configurable<std::vector<float>> cfgTPCNSigmaMuMax{"cfgTPCNSigmaMuMax", {0.f}, "Comma separated list of max TPC nsigma muon for exclusion"};
    Configurable<std::vector<float>> cfgTOFNSigmaElMin{"cfgTOFNSigmaElMin", {-3.f}, "Comma separated list of min TOF nsigma e for inclusion"};
    Configurable<std::vector<float>> cfgTOFNSigmaElMax{"cfgTOFNSigmaElMax", {2.f}, "Comma separated list of max TOF nsigma e for inclusion"};
  } pidcuts;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  HistogramManager* fHistMan;

  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  void init(o2::framework::InitContext&)
  {
    fCurrentRun = 0;

    if (fConfigNbTrackCut > 0 && CheckSize()) {
      for (std::size_t icut = 0; icut < fConfigNbTrackCut; ++icut) {
        AnalysisCompositeCut* cut = new AnalysisCompositeCut(Form("trackcut%zu", icut), Form("trackcut%zu", icut));
        cut->AddCut(GetTrackCut(icut));
        cut->AddCut(GetPIDCut(icut));
        fTrackCuts.push_back(*cut);
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

  bool CheckSize()
  {
    auto veceta = (std::vector<float>)trackcuts.cfgEtaMax;
    auto vecptmin = (std::vector<float>)trackcuts.cfgPtMin;
    auto vecptmax = (std::vector<float>)trackcuts.cfgPtMax;
    auto vecDCAxymax = (std::vector<float>)trackcuts.cfgDCAxyMax;
    auto vecDCAzmax = (std::vector<float>)trackcuts.cfgDCAzMax;
    auto vecIsSPDfirst = (std::vector<int>)trackcuts.cfgIsSPDfirst;
    auto vecIsSPDany = (std::vector<int>)trackcuts.cfgIsSPDany;
    auto vecIsITSibAny = (std::vector<int>)trackcuts.cfgIsITSibAny;
    auto vecITSchi2max = (std::vector<float>)trackcuts.cfgITSchi2Max;
    auto vecITSnclsmin = (std::vector<float>)trackcuts.cfgITSnclsMin;
    auto vecITSnclsmax = (std::vector<float>)trackcuts.cfgITSnclsMax;
    auto vecTPCchi2max = (std::vector<float>)trackcuts.cfgTPCchi2Max;
    auto vecTPCnclsmin = (std::vector<float>)trackcuts.cfgTPCnclsMin;
    auto vecTPCnclsmax = (std::vector<float>)trackcuts.cfgTPCnclsMax;
    auto vecTPCnclsCRmin = (std::vector<float>)trackcuts.cfgTPCnclsCRMin;
    auto vecTPCnclsCRmax = (std::vector<float>)trackcuts.cfgTPCnclsCRMax;
    auto vecIsDalitzLeg = (std::vector<int>)trackcuts.cfgIsDalitzLeg;
    auto vecPIDmode = (std::vector<int>)pidcuts.cfgPIDmode;
    auto vecrejbadtof = (std::vector<int>)pidcuts.cfgRejBadTOF;
    auto vecTPCnsigmaelmin = (std::vector<float>)pidcuts.cfgTPCNSigmaElMin;
    auto vecTPCnsigmaelmax = (std::vector<float>)pidcuts.cfgTPCNSigmaElMax;
    auto vecTPCnsigmapimin = (std::vector<float>)pidcuts.cfgTPCNSigmaPiMin;
    auto vecTPCnsigmapimax = (std::vector<float>)pidcuts.cfgTPCNSigmaPiMax;
    auto vecTPCnsigmaprmin = (std::vector<float>)pidcuts.cfgTPCNSigmaPrMin;
    auto vecTPCnsigmaprmax = (std::vector<float>)pidcuts.cfgTPCNSigmaPrMax;
    auto vecTPCnsigmakamin = (std::vector<float>)pidcuts.cfgTPCNSigmaKaMin;
    auto vecTPCnsigmakamax = (std::vector<float>)pidcuts.cfgTPCNSigmaKaMax;
    auto vecTPCnsigmamumin = (std::vector<float>)pidcuts.cfgTPCNSigmaMuMin;
    auto vecTPCnsigmamumax = (std::vector<float>)pidcuts.cfgTPCNSigmaMuMax;
    auto vecTOFnsigmaelmin = (std::vector<float>)pidcuts.cfgTOFNSigmaElMin;
    auto vecTOFnsigmaelmax = (std::vector<float>)pidcuts.cfgTOFNSigmaElMax;

    if (veceta.size() != fConfigNbTrackCut)
      return false;
    if (vecptmin.size() != fConfigNbTrackCut)
      return false;
    if (vecptmax.size() != fConfigNbTrackCut)
      return false;
    if (vecDCAxymax.size() != fConfigNbTrackCut)
      return false;
    if (vecDCAzmax.size() != fConfigNbTrackCut)
      return false;
    if (vecIsSPDfirst.size() != fConfigNbTrackCut)
      return false;
    if (vecIsSPDany.size() != fConfigNbTrackCut)
      return false;
    if (vecIsITSibAny.size() != fConfigNbTrackCut)
      return false;
    if (vecITSchi2max.size() != fConfigNbTrackCut)
      return false;
    if (vecITSnclsmin.size() != fConfigNbTrackCut)
      return false;
    if (vecITSnclsmax.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCchi2max.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnclsmin.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnclsmax.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnclsCRmin.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnclsCRmax.size() != fConfigNbTrackCut)
      return false;
    if (vecIsDalitzLeg.size() != fConfigNbTrackCut)
      return false;
    if (vecPIDmode.size() != fConfigNbTrackCut)
      return false;
    if (vecrejbadtof.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnsigmaelmin.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnsigmaelmax.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnsigmapimin.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnsigmapimax.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnsigmaprmin.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnsigmaprmax.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnsigmakamin.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnsigmakamax.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnsigmamumin.size() != fConfigNbTrackCut)
      return false;
    if (vecTPCnsigmamumax.size() != fConfigNbTrackCut)
      return false;
    if (vecTOFnsigmaelmin.size() != fConfigNbTrackCut)
      return false;
    if (vecTOFnsigmaelmax.size() != fConfigNbTrackCut)
      return false;
    return true;
  }

  AnalysisCompositeCut* GetPIDCut(unsigned int i)
  {
    auto vecPIDmode = (std::vector<int>)pidcuts.cfgPIDmode;
    auto vecrejbadtof = (std::vector<int>)pidcuts.cfgRejBadTOF;
    auto vecTPCnsigmaelmin = (std::vector<float>)pidcuts.cfgTPCNSigmaElMin;
    auto vecTPCnsigmaelmax = (std::vector<float>)pidcuts.cfgTPCNSigmaElMax;
    auto vecTPCnsigmapimin = (std::vector<float>)pidcuts.cfgTPCNSigmaPiMin;
    auto vecTPCnsigmapimax = (std::vector<float>)pidcuts.cfgTPCNSigmaPiMax;
    auto vecTPCnsigmaprmin = (std::vector<float>)pidcuts.cfgTPCNSigmaPrMin;
    auto vecTPCnsigmaprmax = (std::vector<float>)pidcuts.cfgTPCNSigmaPrMax;
    auto vecTPCnsigmakamin = (std::vector<float>)pidcuts.cfgTPCNSigmaKaMin;
    auto vecTPCnsigmakamax = (std::vector<float>)pidcuts.cfgTPCNSigmaKaMax;
    auto vecTPCnsigmamumin = (std::vector<float>)pidcuts.cfgTPCNSigmaMuMin;
    auto vecTPCnsigmamumax = (std::vector<float>)pidcuts.cfgTPCNSigmaMuMax;
    auto vecTOFnsigmaelmin = (std::vector<float>)pidcuts.cfgTOFNSigmaElMin;
    auto vecTOFnsigmaelmax = (std::vector<float>)pidcuts.cfgTOFNSigmaElMax;

    AnalysisCompositeCut* cut_tpc_nSigma = new AnalysisCompositeCut(Form("pid_TPCnSigma_%d", i), Form("pid_TPCnSigma_%d", i), kTRUE);
    AnalysisCut* cuttpc = new AnalysisCut(Form("pidcuttpc_%d", i), Form("pidcuttpc_%d", i));
    if (!fConfigComputeTPCpostCalib) {
      cuttpc->AddCut(VarManager::kTPCnSigmaEl, vecTPCnsigmaelmin.at(i), vecTPCnsigmaelmax.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
      cuttpc->AddCut(VarManager::kTPCnSigmaPi, vecTPCnsigmapimin.at(i), vecTPCnsigmapimax.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
      cuttpc->AddCut(VarManager::kTPCnSigmaKa, vecTPCnsigmakamin.at(i), vecTPCnsigmakamax.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
      cuttpc->AddCut(VarManager::kTPCnSigmaPr, vecTPCnsigmaprmin.at(i), vecTPCnsigmaprmax.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
      cuttpc->AddCut(VarManager::kTPCnSigmaMu, vecTPCnsigmamumin.at(i), vecTPCnsigmamumax.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
    } else {
      cuttpc->AddCut(VarManager::kTPCnSigmaEl_Corr, vecTPCnsigmaelmin.at(i), vecTPCnsigmaelmax.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
      cuttpc->AddCut(VarManager::kTPCnSigmaPi_Corr, vecTPCnsigmapimin.at(i), vecTPCnsigmapimax.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
      cuttpc->AddCut(VarManager::kTPCnSigmaKa, vecTPCnsigmakamin.at(i), vecTPCnsigmakamax.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
      cuttpc->AddCut(VarManager::kTPCnSigmaPr_Corr, vecTPCnsigmaprmin.at(i), vecTPCnsigmaprmax.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
      cuttpc->AddCut(VarManager::kTPCnSigmaMu, vecTPCnsigmamumin.at(i), vecTPCnsigmamumax.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
    }

    if (vecrejbadtof.at(i) > 0) {
      cuttpc->AddCut(VarManager::kTOFbeta, 0., 0.985, true, VarManager::kPin, 0.0, 1e+10, false);
      cuttpc->AddCut(VarManager::kTOFbeta, 1.015, 999999999., true, VarManager::kPin, 0.0, 1e+10, false);
    }
    cut_tpc_nSigma->AddCut(cuttpc);

    if (vecPIDmode.at(i) == 1)
      return cut_tpc_nSigma;

    AnalysisCompositeCut* cut_tof_nSigma = new AnalysisCompositeCut(Form("pid_TOFnSigma_%d", i), Form("pid_TOFnSigma_%d", i), kTRUE);
    AnalysisCut* cuttof = new AnalysisCut(Form("pidcuttof_%d", i), Form("pidcuttof_%d", i));
    if (!fConfigComputeTPCpostCalib) {
      cuttof->AddCut(VarManager::kTPCnSigmaEl, vecTPCnsigmaelmin.at(i), vecTPCnsigmaelmax.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
      cuttof->AddCut(VarManager::kTPCnSigmaPi, vecTPCnsigmapimin.at(i), vecTPCnsigmapimax.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
    } else {
      cuttof->AddCut(VarManager::kTPCnSigmaEl_Corr, vecTPCnsigmaelmin.at(i), vecTPCnsigmaelmax.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
      cuttof->AddCut(VarManager::kTPCnSigmaPi_Corr, vecTPCnsigmapimin.at(i), vecTPCnsigmapimax.at(i), true, VarManager::kPin, 0.0, 1e+10, false);
    }
    cuttof->AddCut(VarManager::kTOFnSigmaEl, vecTOFnsigmaelmin.at(i), vecTOFnsigmaelmax.at(i), false, VarManager::kPin, 0.0, 1e+10, false);
    cut_tof_nSigma->AddCut(cuttof);

    if (vecPIDmode.at(i) == 2)
      return cut_tof_nSigma;

    AnalysisCompositeCut* cut_pid_OR = new AnalysisCompositeCut(Form("e_NSigma_%d", i), Form("e_NSigma_%d", i), kFALSE);
    cut_pid_OR->AddCut(cut_tpc_nSigma);
    cut_pid_OR->AddCut(cut_tof_nSigma);
    return cut_pid_OR;
  }

  AnalysisCut* GetTrackCut(unsigned int i)
  {
    auto veceta = (std::vector<float>)trackcuts.cfgEtaMax;
    auto vecptmin = (std::vector<float>)trackcuts.cfgPtMin;
    auto vecptmax = (std::vector<float>)trackcuts.cfgPtMax;

    auto vecDCAxymax = (std::vector<float>)trackcuts.cfgDCAxyMax;
    auto vecDCAzmax = (std::vector<float>)trackcuts.cfgDCAzMax;

    auto vecIsSPDfirst = (std::vector<int>)trackcuts.cfgIsSPDfirst;
    auto vecIsSPDany = (std::vector<int>)trackcuts.cfgIsSPDany;
    auto vecIsITSibAny = (std::vector<int>)trackcuts.cfgIsITSibAny;

    auto vecITSchi2max = (std::vector<float>)trackcuts.cfgITSchi2Max;
    auto vecITSnclsmin = (std::vector<float>)trackcuts.cfgITSnclsMin;
    auto vecITSnclsmax = (std::vector<float>)trackcuts.cfgITSnclsMax;

    auto vecTPCchi2max = (std::vector<float>)trackcuts.cfgTPCchi2Max;
    auto vecTPCnclsmin = (std::vector<float>)trackcuts.cfgTPCnclsMin;
    auto vecTPCnclsmax = (std::vector<float>)trackcuts.cfgTPCnclsMax;
    auto vecTPCnclsCRmin = (std::vector<float>)trackcuts.cfgTPCnclsCRMin;
    auto vecTPCnclsCRmax = (std::vector<float>)trackcuts.cfgTPCnclsCRMax;

    auto vecIsDalitzLeg = (std::vector<int>)trackcuts.cfgIsDalitzLeg;

    AnalysisCut* cut = new AnalysisCut(Form("tracksel%d", i), Form("tracksel%d", i));
    cut->AddCut(VarManager::kPt, vecptmin.at(i), vecptmax.at(i));
    cut->AddCut(VarManager::kEta, -veceta.at(i), veceta.at(i));

    cut->AddCut(VarManager::kTrackDCAxy, -vecDCAxymax.at(i), vecDCAxymax.at(i));
    cut->AddCut(VarManager::kTrackDCAz, -vecDCAzmax.at(i), vecDCAzmax.at(i));

    if (vecIsSPDfirst.at(i) > 0)
      cut->AddCut(VarManager::kIsSPDfirst, 0.5, 1.5);
    if (vecIsSPDany.at(i) > 0)
      cut->AddCut(VarManager::kIsSPDany, 0.5, 1.5);
    if (vecIsITSibAny.at(i) > 0)
      cut->AddCut(VarManager::kIsITSibAny, 0.5, 1.5);

    cut->AddCut(VarManager::kITSchi2, 0.0, vecITSchi2max.at(i));
    cut->AddCut(VarManager::kITSncls, vecITSnclsmin.at(i), vecITSnclsmax.at(i));

    cut->AddCut(VarManager::kTPCchi2, 0.0, vecTPCchi2max.at(i));
    cut->AddCut(VarManager::kTPCnclsCR, vecTPCnclsCRmin.at(i), vecTPCnclsCRmax.at(i));
    cut->AddCut(VarManager::kTPCncls, vecTPCnclsmin.at(i), vecTPCnclsmax.at(i));

    if (vecIsDalitzLeg.at(i) > 0 && vecIsDalitzLeg.at(i) <= 8)
      cut->AddCut(VarManager::kIsDalitzLeg + vecIsDalitzLeg.at(i) - 1, 0.5, 1.5);
    return cut;
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
        if (fConfigQAIfSelEvt) {
          if (event.isEventSelected()) {
            fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
          }
        } else {
          fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
        }
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
            if (fConfigQAIfSelEvt) {
              if (event.isEventSelected()) {
                fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
              }
            } else {
              fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
            }
          }
        }
      }

      trackSel(static_cast<int>(filterMap), static_cast<int>(prefilterSelected));
    } // end loop over tracks
  }

  void processSkimmed(MyEventsSelected::iterator const& event, MyBarrelTracks const& tracks)
  {
    runTrackSelection<gkEventFillMap, gkTrackFillMap>(event, tracks);
  }
  void processSkimmedWithCov(MyEventsVtxCovSelected::iterator const& event, MyBarrelTracksWithCov const& tracks)
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

struct AnalysisPrefilterSelection {
  SliceCache cache;
  Produces<aod::Prefilter> prefilter;
  Preslice<MyBarrelTracks> perCollision = aod::reducedtrack::reducedeventId;

  // Configurables
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  o2::parameters::GRPMagField* grpmag = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  Filter barrelTracksSelectedPrefilter = aod::dqanalysisflags::isBarrelSelectedPrefilter > 0;

  Partition<soa::Filtered<MyBarrelTracksSelected>> barrelTracksSelected = aod::dqanalysisflags::isBarrelSelected > 0;

  std::map<int, bool> fPrefiltermap;
  AnalysisCompositeCut* fPairCut;
  struct : ConfigurableGroup {
    std::string prefix = "pairprecut_group";
    Configurable<float> cfgMassMax{"cfgMassMax", 0.05f, "max. mass"};
    Configurable<float> cfgOpAngMax{"cfgOpAngMax", 0.05f, "max. opening angle"};
    Configurable<float> cfgPhiVMin{"cfgPhiVMin", 2.f, "min. phiv"};

  } pairprecuts;

  AnalysisCut* GetPairPreCut()
  {
    AnalysisCut* cut = new AnalysisCut("PairPre", "PairPre");
    cut->AddCut(VarManager::kMass, 0.0, pairprecuts.cfgMassMax);
    cut->AddCut(VarManager::kOpeningAngle, 0.0, pairprecuts.cfgOpAngMax);
    cut->AddCut(VarManager::kPairPhiv, pairprecuts.cfgPhiVMin, 3.2);
    return cut;
  }

  void init(o2::framework::InitContext&)
  {
    fCurrentRun = 0;

    ccdb->setURL(ccdburl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    fPairCut = new AnalysisCompositeCut(true);
    fPairCut->AddCut(GetPairPreCut());

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    VarManager::SetDefaultVarNames();
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
  Configurable<std::size_t> fConfigNbTrackCut{"cfgNbTrackCut", 1, "Number of cuts without prefilter cut, need to be consistent with the track selection"};
  Configurable<int> fConfigMixingDepth{"cfgMixingDepth", 100, "Number of Events stored for event mixing"};
  Configurable<std::string> fConfigAddEventMixingHistogram{"cfgAddEventMixingHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  o2::parameters::GRPMagField* grpmag = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected == 1;
  Filter filterTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;

  HistogramManager* fHistMan;
  // NOTE: The bit mask is required to run pairing just based on the desired electron/muon candidate cuts
  uint32_t fTwoTrackFilterMask = 0;
  std::vector<std::vector<TString>> fTrackHistNames;

  NoBinningPolicy<aod::dqanalysisflags::MixingHash> hashBin;

  void init(o2::framework::InitContext& /*context*/)
  {
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
    if (fConfigNbTrackCut > 0 && fConfigNbTrackCut < 31) {
      for (std::size_t icut = 0; icut < fConfigNbTrackCut; ++icut) {
        std::vector<TString> names = {
          Form("PairsBarrelMEPM_trackcut%zu", icut),
          Form("PairsBarrelMEPP_trackcut%zu", icut),
          Form("PairsBarrelMEMM_trackcut%zu", icut)};
        histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
        fTrackHistNames.push_back(names);
        fTwoTrackFilterMask |= (static_cast<uint32_t>(1) << icut);
      }
    }

    DefineHistograms(fHistMan, histNames.Data(), fConfigAddEventMixingHistogram); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());                              // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  template <uint32_t TEventFillMap, int TPairType, typename TTracks1, typename TTracks2>
  void runMixedPairing(TTracks1 const& tracks1, TTracks2 const& tracks2)
  {

    unsigned int ncuts = fTrackHistNames.size();
    std::vector<std::vector<TString>> histNames = fTrackHistNames;

    uint32_t twoTrackFilter = 0;
    for (auto& track1 : tracks1) {
      for (auto& track2 : tracks2) {
        twoTrackFilter = static_cast<uint32_t>(track1.isBarrelSelected()) & static_cast<uint32_t>(track2.isBarrelSelected()) & fTwoTrackFilterMask;

        if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
          continue;
        }
        VarManager::FillPairME<TEventFillMap, TPairType>(track1, track2);

        constexpr bool eventHasQvector = (VarManager::ObjTypes::ReducedEventQvector > 0);
        if constexpr (eventHasQvector) {
          VarManager::FillPairVn<TEventFillMap, TPairType>(track1, track2);
        }
        constexpr bool eventHasQvectorCentr = (VarManager::ObjTypes::CollisionQvect > 0);
        if constexpr (eventHasQvectorCentr) {
          VarManager::FillPairVn<TEventFillMap, TPairType>(track1, track2);
        }

        for (unsigned int icut = 0; icut < ncuts; icut++) {
          if (twoTrackFilter & (static_cast<uint32_t>(1) << icut)) {
            if (track1.sign() * track2.sign() < 0) {
              fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
            } else {
              if (track1.sign() > 0) {
                fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
              } else {
                fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
              }
            }
          } // end if (filter bits)
        } // end for (cuts)
      } // end for (track2)
    } // end for (track1)
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
      runMixedPairing<TEventFillMap, TPairType>(tracks1, tracks2);
    } // end event loop
  }

  Preslice<soa::Filtered<MyBarrelTracksSelected>> perEventsSelectedT = aod::reducedtrack::reducedeventId;

  void processBarrelSkimmed(soa::Filtered<MyEventsHashSelected>& events, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    runSameSide<pairTypeEE, gkEventFillMap>(events, tracks, perEventsSelectedT);
  }
  void processBarrelVnSkimmed(soa::Filtered<MyEventsHashSelectedQvector>& events, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    runSameSide<pairTypeEE, gkEventFillMapWithQvector>(events, tracks, perEventsSelectedT);
  }
  // TODO: This is a dummy process function for the case when the user does not want to run any of the process functions (no event mixing)
  //    If there is no process function enabled, the workflow hangs
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventMixing, processBarrelSkimmed, "Run barrel-barrel mixing on skimmed tracks", false);
  PROCESS_SWITCH(AnalysisEventMixing, processBarrelVnSkimmed, "Run barrel-barrel vn mixing on skimmed tracks", false);
  PROCESS_SWITCH(AnalysisEventMixing, processDummy, "Dummy function", false);
};

struct AnalysisSameEventPairing {

  float mMagField = 0.0;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  OutputObj<THashList> fOutputList{"output"};
  Configurable<std::size_t> fConfigNbTrackCut{"cfgNbTrackCut", 1, "Number of track cuts without prefilter cut, need to be consistent with the track selection"};
  Configurable<std::size_t> fConfigNbPairCut{"cfgNbPairCut", 1, "Number of pair cuts, need to be below 4 right now"};
  Configurable<string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<string> ccdbPath{"ccdb-path", "Users/lm", "base path to the ccdb object"};
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> fConfigAddSEPHistogram{"cfgAddSEPHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> fCollisionSystem{"syst", "pp", "Collision system, pp or PbPb"};
  Configurable<float> fCenterMassEnergy{"energy", 13600, "Center of mass energy in GeV"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected == 1;
  // NOTE: the barrel filter map contains decisions for both electrons and hadrons used in the correlation task
  Filter filterBarrelTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;
  Filter prefilter = aod::dqanalysisflags::isPrefilterVetoed == 0;

  HistogramManager* fHistMan;

  // NOTE: The track filter produced by the barrel track selection contain a number of electron cut decisions and one last cut for hadrons used in the
  //           dilepton - hadron task downstream. So the bit mask is required to select pairs just based on the electron cuts
  // TODO: provide as Configurable the list and names of the cuts which should be used in pairing
  uint32_t fTwoTrackFilterMask = 0;
  std::vector<std::vector<TString>> fTrackHistNames;
  std::vector<AnalysisCompositeCut> fPairCuts;
  struct : ConfigurableGroup {
    std::string prefix = "paircut_group";
    Configurable<std::vector<int>> cfgRej{"cfgRej", {1}, "max. mass"};
    Configurable<std::vector<float>> cfgMassMax{"cfgMassMax", {0.1f}, "max. mass"};
    Configurable<std::vector<float>> cfgMassMin{"cfgMassMin", {0.f}, "min. mass"};
    Configurable<std::vector<float>> cfgOpAngMax{"cfgOpAngMax", {-1.f}, "max. opening angle"};
    Configurable<std::vector<float>> cfgPhiVMin{"cfgPhiVMin", {2.f}, "min. phiv"};

  } paircuts;

  bool CheckSize()
  {
    auto vecRej = (std::vector<int>)paircuts.cfgRej;
    auto vecmassmin = (std::vector<float>)paircuts.cfgMassMin;
    auto vecmassmax = (std::vector<float>)paircuts.cfgMassMax;
    auto vecopamax = (std::vector<float>)paircuts.cfgOpAngMax;
    auto vecphivmin = (std::vector<float>)paircuts.cfgPhiVMin;

    if (vecRej.size() != fConfigNbPairCut)
      return false;
    if (vecmassmin.size() != fConfigNbPairCut)
      return false;
    if (vecmassmax.size() != fConfigNbPairCut)
      return false;
    if (vecopamax.size() != fConfigNbPairCut)
      return false;
    if (vecphivmin.size() != fConfigNbPairCut)
      return false;
    return true;
  }

  AnalysisCompositeCut* GetPairCut(unsigned int i)
  {
    auto vecRej = (std::vector<int>)paircuts.cfgRej;
    auto vecmassmin = (std::vector<float>)paircuts.cfgMassMin;
    auto vecmassmax = (std::vector<float>)paircuts.cfgMassMax;
    auto vecopamax = (std::vector<float>)paircuts.cfgOpAngMax;
    auto vecphivmin = (std::vector<float>)paircuts.cfgPhiVMin;

    if (vecRej.at(i)) {
      AnalysisCompositeCut* cut_pair = new AnalysisCompositeCut(Form("cut_pair_%d", i), Form("cut_pair_%d", i), kFALSE);
      if (vecphivmin.at(i) > 0) {
        AnalysisCut* cutphiv = new AnalysisCut(Form("phiv_%d", i), Form("phiv_%d", i));
        cutphiv->AddCut(VarManager::kPairPhiv, vecphivmin.at(i), 3.2, true);
        cut_pair->AddCut(cutphiv);
      }
      if (vecmassmin.at(i) < vecmassmax.at(i)) {
        AnalysisCut* cutmass = new AnalysisCut(Form("mass_%d", i), Form("mass_%d", i));
        cutmass->AddCut(VarManager::kMass, vecmassmin.at(i), vecmassmax.at(i), true);
        cut_pair->AddCut(cutmass);
      }
      if (vecopamax.at(i) > 0) {
        AnalysisCut* cutopa = new AnalysisCut(Form("opa_%d", i), Form("opa_%d", i));
        cutopa->AddCut(VarManager::kOpeningAngle, 0.0, vecopamax.at(i), true);
        cut_pair->AddCut(cutopa);
      }
      return cut_pair;
    } else {
      AnalysisCompositeCut* cut_pair = new AnalysisCompositeCut(Form("cut_pair_%d", i), Form("cut_pair_%d", i), kTRUE);
      if (vecphivmin.at(i) > 0) {
        AnalysisCut* cutphiv = new AnalysisCut(Form("phiv_%d", i), Form("phiv_%d", i));
        cutphiv->AddCut(VarManager::kPairPhiv, vecphivmin.at(i), 3.2);
        cut_pair->AddCut(cutphiv);
      }
      if (vecmassmin.at(i) < vecmassmax.at(i)) {
        AnalysisCut* cutmass = new AnalysisCut(Form("mass_%d", i), Form("mass_%d", i));
        cutmass->AddCut(VarManager::kMass, vecmassmin.at(i), vecmassmax.at(i));
        cut_pair->AddCut(cutmass);
      }
      if (vecopamax.at(i) > 0) {
        AnalysisCut* cutopa = new AnalysisCut(Form("opa_%d", i), Form("opa_%d", i));
        cutopa->AddCut(VarManager::kOpeningAngle, 0.0, vecopamax.at(i));
        cut_pair->AddCut(cutopa);
      }
      return cut_pair;
    }
  }

  void init(o2::framework::InitContext& /*context*/)
  {
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
    std::vector<TString> names;

    if (fConfigNbPairCut > 0 && CheckSize()) {
      for (std::size_t icut = 0; icut < fConfigNbPairCut; ++icut) {
        fPairCuts.push_back(*GetPairCut(icut));
      }
    }

    if (fConfigNbTrackCut > 0 && fConfigNbTrackCut < 31) {   // if track cuts
      for (std::size_t icut = 0; icut < fConfigNbTrackCut; ++icut) { // loop over track cuts
        fTwoTrackFilterMask |= (static_cast<uint32_t>(1) << icut);
        // no pair cuts
        names = {
          Form("PairsBarrelSEPM_trackcut%zu", icut),
          Form("PairsBarrelSEPP_trackcut%zu", icut),
          Form("PairsBarrelSEMM_trackcut%zu", icut)};
        histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
        fTrackHistNames.push_back(names);

        for (std::size_t iPairCut = 0; iPairCut < fPairCuts.size(); ++iPairCut) { // loop over pair cuts
          names = {
            Form("PairsBarrelSEPM_trackcut%zu_paircut%zu", icut, iPairCut),
            Form("PairsBarrelSEPP_trackcut%zu_paircut%zu", icut, iPairCut),
            Form("PairsBarrelSEMM_trackcut%zu_paircut%zu", icut, iPairCut)};

          histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());

          fTrackHistNames.push_back(names);
        } // end loop (pair cuts)
      } // end loop (track cuts)
    } // end if (track cuts)

    VarManager::SetCollisionSystem((TString)fCollisionSystem, fCenterMassEnergy); // set collision system and center of mass energy

    DefineHistograms(fHistMan, histNames.Data(), fConfigAddSEPHistogram); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());                      // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  // Template function to run same event pairing (barrel-barrel, muon-muon, barrel-muon)
  template <int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks1, typename TTracks2>
  void runSameEventPairing(TEvent const& event, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    if (fCurrentRun != event.runNumber()) {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, event.timestamp());
      if (grpmag != nullptr) {
        mMagField = grpmag->getNominalL3Field();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", event.timestamp());
      }
      VarManager::SetMagneticField(mMagField);
      fCurrentRun = event.runNumber();
    }

    std::vector<std::vector<TString>> histNames = fTrackHistNames;
    uint32_t twoTrackFilter = 0;

    for (auto& [t1, t2] : combinations(tracks1, tracks2)) {
      twoTrackFilter = static_cast<uint32_t>(t1.isBarrelSelected()) & static_cast<uint32_t>(t2.isBarrelSelected()) & fTwoTrackFilterMask;

      if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }
      constexpr bool eventHasQvector = ((TEventFillMap & VarManager::ObjTypes::ReducedEventQvector) > 0);
      constexpr bool eventHasQvectorCentr = ((TEventFillMap & VarManager::ObjTypes::CollisionQvect) > 0);

      // TODO: FillPair functions need to provide a template argument to discriminate between cases when cov matrix is available or not
      VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
      if constexpr (eventHasQvector) {
        VarManager::FillPairVn<TEventFillMap, TPairType>(t1, t2);
      }
      if constexpr (eventHasQvectorCentr) {
        VarManager::FillPairVn<TEventFillMap, TPairType>(t1, t2);
      }

      int iCut = 0;
      for (std::size_t icut = 0; icut < fConfigNbTrackCut; icut++) {
        if (twoTrackFilter & (static_cast<uint32_t>(1) << icut)) {
          if (t1.sign() * t2.sign() < 0) {
            fHistMan->FillHistClass(histNames[iCut][0].Data(), VarManager::fgValues);
          } else {
            if (t1.sign() > 0) {
              fHistMan->FillHistClass(histNames[iCut][1].Data(), VarManager::fgValues);
            } else {
              fHistMan->FillHistClass(histNames[iCut][2].Data(), VarManager::fgValues);
            }
          }
          iCut++;
          for (unsigned int iPairCut = 0; iPairCut < fPairCuts.size(); iPairCut++, iCut++) {
            AnalysisCompositeCut cut = fPairCuts.at(iPairCut);
            if (!(cut.IsSelected(VarManager::fgValues))) // apply pair cuts
              continue;
            if (t1.sign() * t2.sign() < 0) {
              fHistMan->FillHistClass(histNames[iCut][0].Data(), VarManager::fgValues);
            } else {
              if (t1.sign() > 0) {
                fHistMan->FillHistClass(histNames[iCut][1].Data(), VarManager::fgValues);
              } else {
                fHistMan->FillHistClass(histNames[iCut][2].Data(), VarManager::fgValues);
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
    runSameEventPairing<VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMap>(event, tracks, tracks);
  }
  void processDecayToEESkimmedWithCov(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelectedWithCov> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCov>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, tracks, tracks);
  }
  void processDecayToEEPrefilterSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelectedWithPrefilter> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMap>(event, tracks, tracks);
  }
  void processDecayToEESkimmedWithColl(soa::Filtered<MyEventsSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelectedWithColl> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMapWithColl>(event, tracks, tracks);
  }
  void processVnDecayToEESkimmed(soa::Filtered<MyEventsVtxCovSelectedQvector>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCovQvector>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kDecayToEE, gkEventFillMapWithCovQvector, gkTrackFillMap>(event, tracks, tracks);
  }
  // TODO: dummy function for the case when no process function is enabled
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmed, "Run electron-electron pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmedWithCov, "Run electron-electron pairing, with skimmed covariant tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEEPrefilterSkimmed, "Run electron-electron pairing, with skimmed tracks and prefilter from AnalysisPrefilterSelection", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmedWithColl, "Run electron-electron pairing, with skimmed tracks and with collision information", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processVnDecayToEESkimmed, "Run electron-electron pairing, with skimmed tracks for vn", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
    adaptAnalysisTask<AnalysisPrefilterSelection>(cfgc),
    adaptAnalysisTask<AnalysisEventMixing>(cfgc),
    adaptAnalysisTask<AnalysisSameEventPairing>(cfgc)};
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
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
    }

    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }
  } // end loop over histogram classes
}
