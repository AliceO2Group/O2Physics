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
#include "CCDB/BasicCCDBManager.h"
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
#include <TH1F.h>
#include <THashList.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <algorithm>

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
DECLARE_SOA_COLUMN(MixingHash, mixingHash, int);
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelected, isBarrelSelected, int);
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASH", dqanalysisflags::MixingHash);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTS", dqanalysisflags::IsBarrelSelected);
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
using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;
using MyBarrelTracksSelectedWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkEventFillMapWithQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkEventFillMapWithCovQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector;

constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
// constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;

constexpr static int pairTypeEE = VarManager::kDecayToEE;

// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses); // defines histograms for all tasks

struct AnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};
  // TODO: Provide the mixing variables and binning directly via configurables (e.g. vectors of float)
  Configurable<string> fConfigMixingVariables{"cfgMixingVars", "", "Mixing configs separated by a coma, default no mixing"};
  Configurable<string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};

  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;
  AnalysisCompositeCut* fEventCut;

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
      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;"); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                 // provide the list of required variables so that VarManager knows what to fill
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
  }

  template <uint32_t TEventFillMap, typename TEvent>
  void runEventSelection(TEvent const& event)
  {
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

  void processSkimmed(MyEvents::iterator const& event)
  {
    runEventSelection<gkEventFillMap>(event);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
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
  Configurable<string> fConfigCuts{"cfgTrackCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fTrackCuts;

  void init(o2::framework::InitContext&)
  {
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

      DefineHistograms(fHistMan, histDirNames.Data()); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runTrackSelection(TEvent const& event, TTracks const& tracks)
  {
    VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
    // fill event information which might be needed in histograms/cuts that combine track and event properties
    VarManager::FillEvent<TEventFillMap>(event);

    trackSel.reserve(tracks.size());
    uint32_t filterMap = 0;
    int iCut = 0;

    for (auto& track : tracks) {
      filterMap = 0;
      VarManager::FillTrack<TTrackFillMap>(track);
      if (fConfigQA) { // TODO: make this compile time
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }

      iCut = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, iCut++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << iCut);
          if (fConfigQA) { // TODO: make this compile time
            fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
          }
        }
      }
      trackSel(static_cast<int>(filterMap));
    } // end loop over tracks
  }

  void processSkimmed(MyEvents::iterator const& event, MyBarrelTracks const& tracks)
  {
    runTrackSelection<gkEventFillMap, gkTrackFillMap>(event, tracks);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmed, "Run barrel track selection on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummy, "Dummy function", false);
};

struct AnalysisEventMixing {
  OutputObj<THashList> fOutputList{"output"};
  // Here one should provide the list of electron and muon candidate cuts in the same order as specified in the above
  // single particle selection tasks to preserve the correspondence between the track cut name and its
  //  bit position in the cuts bitmap
  // TODO: Create a configurable to specify exactly on which of the bits one should run the event mixing
  Configurable<string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};

  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected == 1;
  Filter filterTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;

  HistogramManager* fHistMan;
  // NOTE: The bit mask is required to run pairing just based on the desired electron candidate cuts
  uint32_t fTwoTrackFilterMask = 0;
  std::vector<std::vector<TString>> fTrackHistNames;

  NoBinningPolicy<aod::dqanalysisflags::MixingHash> hashBin;

  void init(o2::framework::InitContext& context)
  {
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
          fTwoTrackFilterMask |= (uint32_t(1) << icut);
        }
      }
    }

    DefineHistograms(fHistMan, histNames.Data());    // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  template <int TPairType, typename TTracks1, typename TTracks2>
  void runMixedPairing(TTracks1 const& tracks1, TTracks2 const& tracks2)
  {

    unsigned int ncuts = fTrackHistNames.size();
    std::vector<std::vector<TString>> histNames = fTrackHistNames;

    uint32_t twoTrackFilter = 0;
    for (auto& track1 : tracks1) {
      for (auto& track2 : tracks2) {
        if constexpr (TPairType == VarManager::kDecayToEE) {
          twoTrackFilter = uint32_t(track1.isBarrelSelected()) & uint32_t(track2.isBarrelSelected()) & fTwoTrackFilterMask;
        }

        if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
          continue;
        }
        VarManager::FillPairME<TPairType>(track1, track2);

        constexpr bool eventHasQvector = (VarManager::ObjTypes::ReducedEventQvector > 0);
        if constexpr (eventHasQvector) {
          VarManager::FillPairVn<TPairType>(track1, track2);
        }

        for (unsigned int icut = 0; icut < ncuts; icut++) {
          if (twoTrackFilter & (uint32_t(1) << icut)) {
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
        }   // end for (cuts)
      }     // end for (track2)
    }       // end for (track1)
  }

  // barrel-barrel and muon-muon event mixing
  template <int TPairType, uint32_t TEventFillMap, typename TEvents, typename TTracks>
  void runSameSide(TEvents& events, TTracks const& tracks)
  {
    events.bindExternalIndices(&tracks);
    auto tracksTuple = std::make_tuple(tracks);
    GroupSlicer slicerTracks(events, tracksTuple);
    for (auto& [event1, event2] : selfCombinations(hashBin, 100, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event1, VarManager::fgValues);
      auto it1 = slicerTracks.begin();
      for (auto& slice : slicerTracks) {
        if (slice.groupingElement().index() == event1.index()) {
          it1 = slice;
          break;
        }
      }
      auto tracks1 = std::get<TTracks>(it1.associatedTables());
      tracks1.bindExternalIndices(&events);
      auto it2 = slicerTracks.begin();
      for (auto& slice : slicerTracks) {
        if (slice.groupingElement().index() == event2.index()) {
          it2 = slice;
          break;
        }
      }
      auto tracks2 = std::get<TTracks>(it2.associatedTables());
      tracks2.bindExternalIndices(&events);

      runMixedPairing<TPairType>(tracks1, tracks2);
    } // end event loop
  }

  void processBarrelSkimmed(soa::Filtered<MyEventsHashSelected>& events, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    runSameSide<pairTypeEE, gkEventFillMap>(events, tracks);
  }
  void processBarrelVnSkimmed(soa::Filtered<MyEventsHashSelectedQvector>& events, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    runSameSide<pairTypeEE, gkEventFillMapWithQvector>(events, tracks);
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

  Produces<aod::Dileptons> dileptonList;
  Produces<aod::DileptonsExtra> dileptonExtraList;
  Produces<aod::DileptonFlow> dileptonFlowList;

  OutputObj<THashList> fOutputList{"output"};
  Configurable<string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};
  Configurable<string> ccdbPath{"ccdb-path", "Users/lm", "base path to the ccdb object"};
  Configurable<long> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected == 1;
  // NOTE: the barrel filter map contains decisions for both electrons and hadrons used in the correlation task
  Filter filterBarrelTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;

  HistogramManager* fHistMan;

  // NOTE: The track filter produced by the barrel track selection contain a number of electron cut decisions and one last cut for hadrons used in the
  //           dilepton - hadron task downstream. So the bit mask is required to select pairs just based on the electron cuts
  // TODO: provide as Configurable the list and names of the cuts which should be used in pairing
  uint32_t fTwoTrackFilterMask = 0;
  std::vector<std::vector<TString>> fTrackHistNames;

  void init(o2::framework::InitContext& context)
  {
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Keep track of all the histogram class names to avoid composing strings in the event mixing pairing
    TString histNames = "";
    if (context.mOptions.get<bool>("processDecayToEESkimmed") || context.mOptions.get<bool>("processVnDecayToEESkimmed") || context.mOptions.get<bool>("processAllSkimmed")) {
      TString cutNames = fConfigTrackCuts.value;
      if (!cutNames.IsNull()) {
        std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
        for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
          std::vector<TString> names = {
            Form("PairsBarrelSEPM_%s", objArray->At(icut)->GetName()),
            Form("PairsBarrelSEPP_%s", objArray->At(icut)->GetName()),
            Form("PairsBarrelSEMM_%s", objArray->At(icut)->GetName())};
          histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          fTrackHistNames.push_back(names);
          fTwoTrackFilterMask |= (uint32_t(1) << icut);
        }
      }
    }

    DefineHistograms(fHistMan, histNames.Data());    // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    VarManager::SetupTwoProngDCAFitter(5.0f, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, true); // TODO: get these parameters from Configurables
    VarManager::SetupTwoProngFwdDCAFitter(5.0f, true, 200.0f, 1.0e-3f, 0.9f, true);
  }

  // Template function to run same event pairing (barrel-barrel, muon-muon, barrel-muon)
  template <int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks1, typename TTracks2>
  void runSameEventPairing(TEvent const& event, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {

    unsigned int ncuts = fTrackHistNames.size();
    std::vector<std::vector<TString>> histNames = fTrackHistNames;

    uint32_t twoTrackFilter = 0;
    uint32_t dileptonFilterMap = 0;
    uint32_t dileptonMcDecision = 0; // placeholder, copy of the dqEfficiency.cxx one
    dileptonList.reserve(1);
    dileptonExtraList.reserve(1);
    for (auto& [t1, t2] : combinations(tracks1, tracks2)) {
      if constexpr (TPairType == VarManager::kDecayToEE) {
        twoTrackFilter = uint32_t(t1.isBarrelSelected()) & uint32_t(t2.isBarrelSelected()) & fTwoTrackFilterMask;
      }
      if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }
      constexpr bool eventHasQvector = ((TEventFillMap & VarManager::ObjTypes::ReducedEventQvector) > 0);

      // TODO: FillPair functions need to provide a template argument to discriminate between cases when cov matrix is available or not
      VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);

      // TODO: provide the type of pair to the dilepton table (e.g. ee, mumu, emu...)
      dileptonFilterMap = twoTrackFilter;
      dileptonList(event, VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), dileptonFilterMap, dileptonMcDecision);

      if constexpr (eventHasQvector) {
        dileptonFlowList(VarManager::fgValues[VarManager::kU2Q2], VarManager::fgValues[VarManager::kU3Q3], VarManager::fgValues[VarManager::kCos2DeltaPhi], VarManager::fgValues[VarManager::kCos3DeltaPhi]);
      }

      for (unsigned int icut = 0; icut < ncuts; icut++) {
        if (twoTrackFilter & (uint32_t(1) << icut)) {
          if (t1.sign() * t2.sign() < 0) {
            fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
          } else {
            if (t1.sign() > 0) {
              fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
            } else {
              fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
            }
          }
        } // end if (filter bits)
      }   // end for (cuts)
    }     // end loop over pairs
  }

  void processDecayToEESkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMap>(event, tracks, tracks);
  }
  void processVnDecayToEESkimmed(soa::Filtered<MyEventsVtxCovSelectedQvector>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCovQvector>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kDecayToEE, gkEventFillMapWithCovQvector, gkTrackFillMap>(event, tracks, tracks);
  }
  void processAllSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMap>(event, tracks, tracks);
  }
  // TODO: dummy function for the case when no process function is enabled
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmed, "Run electron-electron pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processVnDecayToEESkimmed, "Run electron-electron pairing, with skimmed tracks for vn", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processAllSkimmed, "Run all types of pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
    adaptAnalysisTask<AnalysisEventMixing>(cfgc),
    adaptAnalysisTask<AnalysisSameEventPairing>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses)
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

    // NOTE: The level of detail for histogramming can be controlled via configurables
    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "trigger,cent");
    }

    if (classStr.Contains("Track")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "its,tpcpid,dca,tofpid");
      }
    }

    if (classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_lmee", "");
      }
    }

    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_barrel");
    }

  } // end loop over histogram classes
}
