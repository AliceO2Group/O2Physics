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
//          alexander.tiekoetter@cern.ch
//   Configurable workflow for running several DQ or other PWG analyses

#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/DataModel/ReducedTablesAlice3.h"

#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/collisionAlice3.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/Core/TableHelper.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TObjString.h>
#include <TString.h>

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

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
DECLARE_SOA_COLUMN(MixingHash, mixingHash, int);                                     //! Hash used in event mixing //need to understand
DECLARE_SOA_BITMAP_COLUMN(IsEventSelected, isEventSelected, 32);                     //! Event decision
DECLARE_SOA_BITMAP_COLUMN(IsBarrelSelected, isBarrelSelected, 32);                   //! Barrel track decisions
DECLARE_SOA_COLUMN(BarrelAmbiguityInBunch, barrelAmbiguityInBunch, int8_t);          //! Barrel track in-bunch ambiguity
DECLARE_SOA_COLUMN(BarrelAmbiguityOutOfBunch, barrelAmbiguityOutOfBunch, int8_t);    //! Barrel track out of bunch ambiguity
DECLARE_SOA_BITMAP_COLUMN(IsBarrelSelectedPrefilter, isBarrelSelectedPrefilter, 32); //! Barrel prefilter decisions (joinable to ReducedA3TracksAssoc)

DECLARE_SOA_COLUMN(Massee, massee, float);
DECLARE_SOA_COLUMN(Etaee, etaee, float);
DECLARE_SOA_COLUMN(Rapee, rapee, float);
DECLARE_SOA_COLUMN(Phiee, phiee, float);
DECLARE_SOA_COLUMN(Ptee, ptee, float);
DECLARE_SOA_COLUMN(Lxyee, lxyee, float);
DECLARE_SOA_COLUMN(LxyeePoleMass, lxyeepolemass, float);
DECLARE_SOA_COLUMN(Lzee, lzee, float);
DECLARE_SOA_COLUMN(MultiplicityFT0A, multiplicityFT0AJPsi2ee, float);
DECLARE_SOA_COLUMN(MultiplicityFT0C, multiplicityFT0CJPsi2ee, float);
DECLARE_SOA_COLUMN(PercentileFT0M, percentileFT0MJPsi2ee, float);
DECLARE_SOA_COLUMN(MultiplicityNContrib, multiplicityNContribJPsi2ee, float);
DECLARE_SOA_COLUMN(AmbiguousInBunchPairs, AmbiguousJpsiPairsInBunch, bool);
DECLARE_SOA_COLUMN(AmbiguousOutOfBunchPairs, AmbiguousJpsiPairsOutOfBunch, bool);
DECLARE_SOA_COLUMN(Corrassoc, corrassoc, bool);
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);
// Candidate columns efficiency calculation for prompt-non-prompt JPsi separation
DECLARE_SOA_COLUMN(OniaPt, oniaPt, float);
DECLARE_SOA_COLUMN(OniaY, oniaY, float);
DECLARE_SOA_COLUMN(OniaEta, oniaEta, float);
DECLARE_SOA_COLUMN(OniaPhi, oniaPhi, float);
DECLARE_SOA_COLUMN(OniaVz, oniaVz, float);
DECLARE_SOA_COLUMN(OniaVtxZ, oniaVtxZ, float);
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);                                                            //!  joinable to ReducedEvents
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASHA", dqanalysisflags::MixingHash);                                                            //!  joinable to ReducedEvents
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTS", dqanalysisflags::IsBarrelSelected);                                                    //!  joinable to ReducedA3TracksAssoc
DECLARE_SOA_TABLE(BarrelAmbiguities, "AOD", "DQBARRELAMB", dqanalysisflags::BarrelAmbiguityInBunch, dqanalysisflags::BarrelAmbiguityOutOfBunch); //!  joinable to ReducedBarrelTracks
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTER", dqanalysisflags::IsBarrelSelectedPrefilter);                                                  //!  joinable to ReducedA3TracksAssoc

DECLARE_SOA_TABLE(OniaMCTruth, "AOD", "MCTRUTHONIA", dqanalysisflags::OniaPt, dqanalysisflags::OniaEta, dqanalysisflags::OniaY, dqanalysisflags::OniaPhi, dqanalysisflags::OniaVz, dqanalysisflags::OniaVtxZ, dqanalysisflags::MultiplicityFT0A, dqanalysisflags::MultiplicityFT0C, dqanalysisflags::PercentileFT0M, dqanalysisflags::MultiplicityNContrib);

} // namespace o2::aod

// TODO: USE PROPER TABLES

using MyEvents = soa::Join<aod::ReducedA3Events, aod::ReducedA3MCEventLabels>;
using MyEventsSelected = soa::Join<aod::ReducedA3Events, aod::EventCuts, aod::ReducedA3MCEventLabels>;
using MyEventsVtxCov = soa::Join<aod::ReducedA3Events, aod::ReducedA3EventsVtxCov, aod::ReducedA3MCEventLabels>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedA3Events, aod::ReducedA3EventsVtxCov, aod::EventCuts, aod::ReducedA3MCEventLabels>;

using MyBarrelAssocs = soa::Join<aod::ReducedA3TracksAssoc, aod::BarrelTrackCuts>;
using MyBarrelAssocsPrefilter = soa::Join<aod::ReducedA3TracksAssoc, aod::BarrelTrackCuts, aod::Prefilter>;

using MyBarrelTracks = soa::Join<aod::ReducedA3Tracks, aod::ReducedA3TracksBarrel,
                                 aod::ReducedA3PIDTOF, aod::ReducedA3PIDRich, aod::ReducedA3PIDRichSignals,
                                 aod::ReducedA3TracksBarrelLabels>;

using MyBarrelTracksWithCov = soa::Join<aod::ReducedA3Tracks, aod::ReducedA3TracksBarrel, aod::ReducedA3TracksBarrelCov,
                                        aod::ReducedA3PIDTOF, aod::ReducedA3PIDRich, aod::ReducedA3PIDRichSignals,
                                        aod::ReducedA3TracksBarrelLabels>;

using MyBarrelTracksWithCovWithAmbiguities = soa::Join<aod::ReducedA3Tracks, aod::ReducedA3TracksBarrel, aod::ReducedA3TracksBarrelCov,
                                                       aod::ReducedA3PIDTOF, aod::ReducedA3PIDRich, aod::ReducedA3PIDRichSignals,
                                                       aod::BarrelAmbiguities, aod::ReducedA3TracksBarrelLabels>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent;
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventVtxCov;

constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;

// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses, const char* histGroups); // defines histograms for all tasks

constexpr int TWO_PRONG = 2;
constexpr int THREE_PRONG = 3;

// Analysis task that produces event decisions and the Hash table used in event mixing
struct AnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};
  Configurable<std::string> fConfigMixingVariables{"cfgMixingVars", "", "Mixing configs separated by a comma, default no mixing"};
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigEventCutsJSON{"cfgEventCutsJSON", "", "Additional event cuts specified in JSON format"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddEventMCHistogram{"cfgAddEventMCHistogram", "generator", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Add event histograms defined via JSON formatting (see HistogramsLibrary)"};

  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;

  AnalysisCompositeCut* fEventCut;

  std::map<int64_t, bool> fSelMap; // key: reduced event global index, value: event selection decision

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }
    VarManager::SetDefaultVarNames();

    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    if (eventCutStr != "") {
      AnalysisCut* cut = dqcuts::GetAnalysisCut(eventCutStr.Data());
      if (cut != nullptr) {
        fEventCut->AddCut(cut);
      }
    }
    // Additional cuts via JSON
    TString eventCutJSONStr = fConfigEventCutsJSON.value;
    if (eventCutJSONStr != "") {
      std::vector<AnalysisCut*> jsonCuts = dqcuts::GetCutsFromJSON(eventCutJSONStr.Data());
      for (const auto& cutIt : jsonCuts) {
        fEventCut->AddCut(cutIt);
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(true);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;", fConfigAddEventHistogram.value.data());
      DefineHistograms(fHistMan, "EventsMC", fConfigAddEventMCHistogram.value.data());
      dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // aditional histograms via JSON
      VarManager::SetUseVars(fHistMan->GetUsedVars());
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

  void runEventSelection(MyEventsVtxCov const& events, ReducedA3MCEvents const& mcEvents)
  {
    fSelMap.clear();

    for (const auto& event : events) {
      // Reset the fValues array and fill event observables
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEventAlice3<gkEventFillMap>(event);
      if (event.has_reducedA3MCEvent()) {
        VarManager::FillEventAlice3<VarManager::ObjTypes::ReducedEventMC>(event.reducedA3MCEvent());
      }

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
      if (fMixHandler != nullptr) {
        int hh = fMixHandler->FindEventCategory(VarManager::fgValues);
        hash(hh);
      }
    }

    for (const auto& event : mcEvents) {
      // Reset the fValues array and fill event observables
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEventAlice3<VarManager::ObjTypes::ReducedEventMC>(event);
      if (fConfigQA) {
        fHistMan->FillHistClass("EventsMC", VarManager::fgValues);
      }
    }
  }

  void publishSelections(MyEventsVtxCov const& events)
  {
    // publish the table
    uint32_t evSel = static_cast<uint32_t>(0);
    for (const auto& event : events) {
      evSel = 0;
      if (fSelMap[event.globalIndex()]) { // event passed the user cuts
        evSel |= (static_cast<uint32_t>(1) << 0);
      }
      eventSel(evSel);
    }
  }

  void processSkimmed(MyEventsVtxCov const& events, aod::ReducedA3MCEvents const& mcEvents)
  {
    runEventSelection(events, mcEvents);
    publishSelections(events);
  }

  void processDummy(MyEventsVtxCov const&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDummy, "Dummy function", true);
};

// Produces a table with barrel track decisions (joinable to the ReducedA3TracksAssociations)
// Here one should add all the track cuts needed through the workflow (e.g. cuts for same-even pairing, electron prefiltering, track for dilepton-track correlations)
struct AnalysisTrackSelection {
  Produces<aod::BarrelTrackCuts> trackSel;
  Produces<aod::BarrelAmbiguities> trackAmbiguities;
  OutputObj<THashList> fOutputList{"output"};

  Configurable<std::string> fConfigCuts{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigCutsJSON{"cfgBarrelTrackCutsJSON", "", "Additional list of barrel track cuts in JSON format"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};
  Configurable<bool> fConfigPublishAmbiguity{"cfgPublishAmbiguity", true, "If true, publish ambiguity table and fill QA histograms"};

  Configurable<std::string> fConfigMCSignals{"cfgTrackMCSignals", "", "Comma separated list of MC signals"};
  Configurable<std::string> fConfigMCSignalsJSON{"cfgTrackMCsignalsJSON", "", "Additional list of MC signals via JSON"};

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut*> fTrackCuts;
  std::vector<MCSignal*> fMCSignals; // list of signals to be checked
  std::vector<TString> fHistNamesReco;
  std::vector<TString> fHistNamesMCMatched;

  std::map<int64_t, std::vector<int64_t>> fNAssocsInBunch;    // key: track global index, value: vector of global index for events associated in-bunch (events that have in-bunch pileup or splitting)
  std::map<int64_t, std::vector<int64_t>> fNAssocsOutOfBunch; // key: track global index, value: vector of global index for events associated out-of-bunch (events that have no in-bunch pileup)

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }
    VarManager::SetDefaultVarNames();

    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // add extra cuts from JSON
    TString addTrackCutsStr = fConfigCutsJSON.value;
    if (addTrackCutsStr != "") {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
      for (const auto& t : addTrackCuts) {
        fTrackCuts.push_back(reinterpret_cast<AnalysisCompositeCut*>(t));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    TString configSigNamesStr = fConfigMCSignals.value;
    std::unique_ptr<TObjArray> sigNamesArray(configSigNamesStr.Tokenize(","));
    // Setting the MC signals
    for (int isig = 0; isig < sigNamesArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(sigNamesArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() != 1) { // NOTE: only 1 prong signals
          continue;
        }
        fMCSignals.push_back(sig);
      }
    }
    // Add the MCSignals from the JSON config
    TString addMCSignalsStr = fConfigMCSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (const auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() != 1) { // NOTE: only 1 prong signals
          continue;
        }
        fMCSignals.push_back(mcIt);
      }
    }

    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(true);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      // Configure histogram classes for each track cut;
      // Add histogram classes for each track cut and for each requested MC signal (reconstructed tracks with MC truth)
      TString histClasses = "AssocsBarrel_BeforeCuts;";
      for (const auto& cut : fTrackCuts) {
        TString nameStr = Form("AssocsBarrel_%s", cut->GetName());
        fHistNamesReco.push_back(nameStr);
        histClasses += Form("%s;", nameStr.Data());
        for (const auto& sig : fMCSignals) {
          TString nameStr2 = Form("AssocsCorrectBarrel_%s_%s", cut->GetName(), sig->GetName());
          fHistNamesMCMatched.push_back(nameStr2);
          histClasses += Form("%s;", nameStr2.Data());
          nameStr2 = Form("AssocsIncorrectBarrel_%s_%s", cut->GetName(), sig->GetName());
          fHistNamesMCMatched.push_back(nameStr2);
          histClasses += Form("%s;", nameStr2.Data());
        }
      }

      DefineHistograms(fHistMan, histClasses.Data(), fConfigAddTrackHistogram.value.data());
      if (fConfigPublishAmbiguity) {
        DefineHistograms(fHistMan, "TrackBarrel_AmbiguityInBunch;TrackBarrel_AmbiguityOutOfBunch;", "ambiguity");
      }
      dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                       // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  void runTrackSelection(ReducedA3TracksAssoc const& assocs, MyEventsVtxCovSelected const& /*events*/, MyBarrelTracksWithCov const& tracks, ReducedA3MCEvents const& /*eventsMC*/, ReducedA3MCTracks const& tracksMC)
  {
    fNAssocsInBunch.clear();
    fNAssocsOutOfBunch.clear();

    trackSel.reserve(assocs.size());
    trackAmbiguities.reserve(tracks.size());

    // Loop over associations
    for (const auto& assoc : assocs) {
      auto event = assoc.template reducedA3event_as<MyEventsVtxCovSelected>();
      if (!event.isEventSelected_bit(0)) {
        trackSel(0);
        continue;
      }

      VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
      // fill event information which might be needed in histograms/cuts that combine track and event properties
      VarManager::FillEventAlice3<gkEventFillMapWithCov>(event);
      if (event.has_reducedA3MCEvent()) {
        VarManager::FillEventAlice3<VarManager::ObjTypes::ReducedEventMC>(event.reducedA3MCEvent());
      }

      auto track = assoc.template reducedA3track_as<MyBarrelTracksWithCov>();
      VarManager::FillTrackAlice3<gkTrackFillMapWithCov>(track);
      // compute quantities which depend on the associated collision, such as DCA
      VarManager::FillTrackCollision<gkTrackFillMapWithCov>(track, event);

      bool isCorrectAssoc = false;
      if (track.has_reducedA3MCTrack()) {
        auto trackMC = track.reducedA3MCTrack();
        auto eventMCfromTrack = trackMC.reducedA3MCEvent();
        if (event.has_reducedA3MCEvent()) {
          isCorrectAssoc = (eventMCfromTrack.globalIndex() == event.reducedA3MCEvent().globalIndex());
        }
        VarManager::FillTrackMC(tracksMC, trackMC);
        VarManager::FillResolutions(trackMC, track);
      }

      if (fConfigQA) {
        fHistMan->FillHistClass("AssocsBarrel_BeforeCuts", VarManager::fgValues);
      }

      int iCut = 0;
      uint32_t filterMap = static_cast<uint32_t>(0);
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, iCut++) {
        if ((*cut)->IsSelected(VarManager::fgValues)) {
          filterMap |= (static_cast<uint32_t>(1) << iCut);
          if (fConfigQA) {
            fHistMan->FillHistClass(fHistNamesReco[iCut], VarManager::fgValues);
          }
        }
      } // end loop over cuts
      trackSel(filterMap);

      // compute MC matching decisions and fill histograms for matched associations
      int isig = 0;
      if (filterMap > 0 && track.has_reducedA3MCTrack()) {
        // loop over all MC signals
        for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {
          // check if this MC signal is matched
          if ((*sig)->CheckSignal(true, track.reducedA3MCTrack())) {
            // mcDecision |= (static_cast<uint32_t>(1) << isig);
            //  loop over cuts and fill histograms for the cuts that are fulfilled
            for (unsigned int icut = 0; icut < fTrackCuts.size(); icut++) {
              if (filterMap & (static_cast<uint32_t>(1) << icut)) {
                if (isCorrectAssoc) {
                  fHistMan->FillHistClass(fHistNamesMCMatched[icut * 2 * fMCSignals.size() + 2 * isig].Data(), VarManager::fgValues);
                } else {
                  fHistMan->FillHistClass(fHistNamesMCMatched[icut * 2 * fMCSignals.size() + 2 * isig + 1].Data(), VarManager::fgValues);
                }
              }
            } // end loop over cuts
          }
        }
      } // end if (filterMap > 0)

      // count the number of associations per track
      if (fConfigPublishAmbiguity && filterMap > 0) {
        if (event.isEventSelected_bit(1)) {
          // for this track, count the number of associated collisions with in-bunch pileup and out of bunch associations
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
    // TODO: some tracks can be associated to both collisions that have in bunch pileup and collisions from different bunches
    //       So one could QA these tracks separately
    if (fConfigPublishAmbiguity) {
      if (fConfigQA) {
        for (const auto& [trackIdx, evIndices] : fNAssocsInBunch) {
          if (evIndices.size() == 1) {
            continue;
          }
          auto track = tracks.rawIteratorAt(trackIdx);
          VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
          VarManager::FillTrackAlice3<gkTrackFillMapWithCov>(track);
          VarManager::fgValues[VarManager::kBarrelNAssocsInBunch] = static_cast<float>(evIndices.size());
          fHistMan->FillHistClass("TrackBarrel_AmbiguityInBunch", VarManager::fgValues);
        } // end loop over in-bunch ambiguous tracks

        for (const auto& [trackIdx, evIndices] : fNAssocsOutOfBunch) {
          if (evIndices.size() == 1) {
            continue;
          }
          auto track = tracks.rawIteratorAt(trackIdx);
          VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
          VarManager::FillTrackAlice3<gkTrackFillMapWithCov>(track);
          VarManager::fgValues[VarManager::kBarrelNAssocsOutOfBunch] = static_cast<float>(evIndices.size());
          fHistMan->FillHistClass("TrackBarrel_AmbiguityOutOfBunch", VarManager::fgValues);
        } // end loop over out-of-bunch ambiguous tracks
      }

      // publish the ambiguity table
      for (const auto& track : tracks) {
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
    }
  } // end runTrackSelection()

  void processSkimmedWithCov(ReducedA3TracksAssoc const& assocs, MyEventsVtxCovSelected const& events, MyBarrelTracksWithCov const& tracks, ReducedA3MCEvents const& eventsMC, ReducedA3MCTracks const& tracksMC)
  {
    runTrackSelection(assocs, events, tracks, eventsMC, tracksMC);
  }

  void processDummy(MyEvents const&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmedWithCov, "Run barrel track selection on DQ skimmed tracks w/ cov matrix associations", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummy, "Dummy function", true);
};

struct AnalysisPrefilterSelection {
  Produces<aod::Prefilter> prefilter; // joinable with ReducedA3TracksAssoc

  // Configurables
  Configurable<std::string> fConfigPrefilterTrackCut{"cfgPrefilterTrackCut", "", "Prefilter track cut"};
  Configurable<std::string> fConfigPrefilterPairCut{"cfgPrefilterPairCut", "", "Prefilter pair cut"};
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Track cuts for which to run the prefilter"};
  // Track related options
  Configurable<bool> fPropTrack{"cfgPropTrack", false, "Propgate tracks to associated collision to recalculate DCA and momentum vector"};

  std::map<uint32_t, uint32_t> fPrefilterMap;
  AnalysisCompositeCut* fPairCut;
  uint32_t fPrefilterMask;
  int fPrefilterCutBit;

  PresliceUnsorted<aod::ReducedA3TracksAssoc> trackAssocsPerCollision = aod::reducedA3track_association::reducedA3eventId;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    bool runPrefilter = true;
    // get the list of track cuts to be prefiltered
    TString trackCutsStr = fConfigTrackCuts.value;
    TObjArray* objArrayTrackCuts = nullptr;
    if (!trackCutsStr.IsNull()) {
      objArrayTrackCuts = trackCutsStr.Tokenize(",");
      if (objArrayTrackCuts == nullptr) {
        runPrefilter = false;
      }
    } else {
      LOG(warn) << " No track cuts to prefilter! Prefilter will not be run";
      runPrefilter = false;
    }
    // get the cut to be used as loose selection
    TString prefilterTrackCutStr = fConfigPrefilterTrackCut.value;
    if (prefilterTrackCutStr.IsNull()) {
      LOG(warn) << " No prefilter loose selection specified! Prefilter will not be run";
      runPrefilter = false;
    }

    fPrefilterMask = 0;
    fPrefilterCutBit = -1;
    if (runPrefilter) {
      // get the list of cuts that were computed in the barrel track-selection task and create a bit mask
      //  to mark just the ones we want to apply a prefilter on
      string trackCuts;
      getTaskOptionValue<string>(context, "analysis-track-selection", "cfgTrackCuts", trackCuts, false);
      TString allTrackCutsStr = trackCuts;
      // check also the cuts added via JSON and add them to the string of cuts
      getTaskOptionValue<string>(context, "analysis-track-selection", "cfgBarrelTrackCutsJSON", trackCuts, false);
      TString addTrackCutsStr = trackCuts;
      if (addTrackCutsStr != "") {
        std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
        for (const auto& t : addTrackCuts) {
          allTrackCutsStr += Form(",%s", t->GetName());
        }
      }

      std::unique_ptr<TObjArray> objArray(allTrackCutsStr.Tokenize(","));
      if (objArray == nullptr) {
        LOG(fatal) << " Not getting any track cuts from the barrel-track-selection ";
      }
      if (objArray->FindObject(prefilterTrackCutStr.Data()) == nullptr) {
        LOG(fatal) << " Prefilter track cut not among the cuts calculated by the track-selection task! ";
      }
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        TString tempStr = objArray->At(icut)->GetName();
        if (objArrayTrackCuts->FindObject(tempStr.Data()) != nullptr) {
          fPrefilterMask |= (static_cast<uint32_t>(1) << icut);
        }
        if (tempStr.CompareTo(fConfigPrefilterTrackCut.value) == 0) {
          fPrefilterCutBit = icut;
        }
      }
      // setup the prefilter pair cut
      fPairCut = new AnalysisCompositeCut(true);
      TString pairCutStr = fConfigPrefilterPairCut.value;
      if (!pairCutStr.IsNull()) {
        fPairCut = dqcuts::GetCompositeCut(pairCutStr.Data());
      }
    }
    if (fPrefilterMask == static_cast<uint32_t>(0) || fPrefilterCutBit < 0) {
      LOG(warn) << "No specified loose cut or track cuts for prefiltering. This task will do nothing.";
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    VarManager::SetDefaultVarNames();

    VarManager::SetupTwoProngDCAFitter(5.0f, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, true); // TODO: get these parameters from Configurables
    VarManager::SetupTwoProngFwdDCAFitter(5.0f, true, 200.0f, 1.0e-3f, 0.9f, true);
  }

  void runPrefilter(MyEvents::iterator const& event, MyBarrelAssocs const& assocs, MyBarrelTracks const& /*tracks*/)
  {
    if (fPrefilterCutBit < 0 || fPrefilterMask == 0) {
      return;
    }

    for (const auto& [assoc1, assoc2] : o2::soa::combinations(assocs, assocs)) {
      auto track1 = assoc1.template reducedA3track_as<MyBarrelTracks>();
      auto track2 = assoc2.template reducedA3track_as<MyBarrelTracks>();

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
      VarManager::FillPairAlice3<VarManager::kDecayToEE, gkTrackFillMap>(track1, track2);
      if (fPropTrack) {
        VarManager::FillPairCollision<VarManager::kDecayToEE, gkTrackFillMap>(event, track1, track2);
      }
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

  void processBarrelSkimmed(MyEvents const& events, MyBarrelAssocs const& assocs, MyBarrelTracks const& tracks)
  {

    fPrefilterMap.clear();

    for (const auto& event : events) {
      auto groupedAssocs = assocs.sliceBy(trackAssocsPerCollision, event.globalIndex());
      if (groupedAssocs.size() > 1) {
        runPrefilter(event, groupedAssocs, tracks);
      }
    }

    uint32_t mymap = -1;
    // If cuts were not configured, then produce a map with all 1's and publish it for all associations
    if (fPrefilterCutBit < 0 || fPrefilterMask == 0) {
      for (int i = 0; i < assocs.size(); ++i) {
        prefilter(mymap);
      }
    } else {
      for (const auto& assoc : assocs) {
        // TODO: just use the index from the assoc (no need to cast the whole track)
        auto track = assoc.template reducedA3track_as<MyBarrelTracks>();
        mymap = -1;
        if (fPrefilterMap.find(track.globalIndex()) != fPrefilterMap.end()) {
          // NOTE: publish the bitwise negated bits (~), so there will be zeroes for cuts that failed the prefiltering and 1 everywhere else
          mymap = ~fPrefilterMap[track.globalIndex()];
          prefilter(mymap);
        } else {
          prefilter(mymap); // track did not pass the prefilter selections, so publish just 1's
        }
      }
    }
  }

  void processDummy(MyEvents const&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisPrefilterSelection, processBarrelSkimmed, "Run Prefilter selection on reduced tracks", false);
  PROCESS_SWITCH(AnalysisPrefilterSelection, processDummy, "Do nothing", true);
};

// Run the same-event pairing
// This task assumes that both legs of the resonance fulfill the same cuts (symmetric decay channel)
//  Runs combinatorics for barrel-barrel combinations
// The task implements also process functions for running event mixing
struct AnalysisSameEventPairing {

  Produces<aod::Dielectrons> dielectronList;
  Produces<aod::DielectronsExtra> dielectronsExtraList;
  Produces<aod::DielectronsAll> dielectronAllList;
  Produces<aod::OniaMCTruth> MCTruthTableEffi;

  o2::base::MatLayerCylSet* fLUT = nullptr;
  OutputObj<THashList> fOutputList{"output"};

  struct : ConfigurableGroup {
    Configurable<std::string> track{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
    Configurable<std::string> pair{"cfgPairCuts", "", "Comma separated list of pair cuts, !!! Use only if you know what you are doing, otherwise leave empty"};
    Configurable<std::string> MCgenAcc{"cfgMCGenAccCut", "", "cut for MC generated particles acceptance"};
    // TODO: Add pair cuts via JSON
  } fConfigCuts;

  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddSEPHistogram{"cfgAddSEPHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};

  struct : ConfigurableGroup {
    Configurable<bool> flatTables{"cfgFlatTables", false, "Produce a single flat tables with all relevant information of the pairs and single tracks"};
    Configurable<bool> propToPCA{"cfgPropToPCA", false, "Propagate tracks to secondary vertex"};
    Configurable<std::string> collisionSystem{"syst", "pp", "Collision system, pp or PbPb"};
    Configurable<float> centerMassEnergy{"energy", 13600, "Center of mass energy in GeV"};
  } fConfigOptions;

  struct : ConfigurableGroup {
    Configurable<std::string> genSignals{"cfgBarrelMCGenSignals", "", "Comma separated list of MC signals (generated)"};
    Configurable<std::string> genSignalsJSON{"cfgMCGenSignalsJSON", "", "Additional list of MC signals (generated) via JSON"};
    Configurable<std::string> recSignals{"cfgBarrelMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};
    Configurable<std::string> recSignalsJSON{"cfgMCRecSignalsJSON", "", "Comma separated list of MC signals (reconstructed) via JSON"};
    Configurable<bool> skimSignalOnly{"cfgSkimSignalOnly", false, "Configurable to select only matched candidates"};
  } fConfigMC;

  // Track related options
  Configurable<bool> fPropTrack{"cfgPropTrack", true, "Propgate tracks to associated collision to recalculate DCA and momentum vector"};

  // Filter filterEventSelected = aod::dqanalysisflags::isEventSelected & uint32_t(1);
  Filter eventFilter = aod::dqanalysisflags::isEventSelected > static_cast<uint32_t>(0);

  HistogramManager* fHistMan;

  // keep histogram class names in maps, so we don't have to buld their names in the pair loops
  std::map<int, std::vector<TString>> fTrackHistNames;
  std::map<int, std::vector<TString>> fBarrelHistNamesMCmatched;
  std::vector<MCSignal*> fRecMCSignals;
  std::vector<MCSignal*> fGenMCSignals;

  std::vector<AnalysisCompositeCut> fPairCuts;
  AnalysisCompositeCut fMCGenAccCut;
  bool fUseMCGenAccCut = false;

  uint32_t fTrackFilterMask; // mask for the track cuts required in this task to be applied on the barrel cuts produced upstream
  int fNCutsBarrel;
  int fNPairCuts;
  bool fHasTwoProngGenMCsignals = false;

  bool fEnableBarrelHistos;

  PresliceUnsorted<MyBarrelAssocsPrefilter> trackAssocsPerCollision = aod::reducedA3track_association::reducedA3eventId;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }
    bool isMCGen = context.mOptions.get<bool>("processMCGen") || context.mOptions.get<bool>("processMCGenWithGrouping") || context.mOptions.get<bool>("processBarrelOnlySkimmed");
    VarManager::SetDefaultVarNames();

    fEnableBarrelHistos = context.mOptions.get<bool>("processBarrelOnlySkimmed");

    // Keep track of all the histogram class names to avoid composing strings in the pairing loop
    TString histNames = "";
    TString cutNamesStr = fConfigCuts.pair.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // get the list of cuts for tracks, check that they were played by the barrel selection tasks
    //   and make a mask for active cuts (barrel selection tasks may run more cuts, needed for other analyses)
    TString trackCutsStr = fConfigCuts.track.value;
    TObjArray* objArrayTrackCuts = nullptr;
    if (!trackCutsStr.IsNull()) {
      objArrayTrackCuts = trackCutsStr.Tokenize(",");
    }

    // Setting the MC rec signal names
    TString sigNamesStr = fConfigMC.recSignals.value;
    std::unique_ptr<TObjArray> objRecSigArray(sigNamesStr.Tokenize(","));
    for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() != TWO_PRONG) { // NOTE: 2-prong signals required
          continue;
        }
        fRecMCSignals.push_back(sig);
      }
    }

    // Add the MCSignals from the JSON config
    TString addMCSignalsStr = fConfigMC.recSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (const auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() != TWO_PRONG) { // NOTE: only 2 prong signals
          continue;
        }
        fRecMCSignals.push_back(mcIt);
      }
    }

    // get the barrel track selection cuts
    string tempCuts;
    getTaskOptionValue<string>(context, "analysis-track-selection", "cfgTrackCuts", tempCuts, false);
    TString tempCutsStr = tempCuts;
    // check also the cuts added via JSON and add them to the string of cuts
    getTaskOptionValue<string>(context, "analysis-track-selection", "cfgBarrelTrackCutsJSON", tempCuts, false);
    TString addTrackCutsStr = tempCuts;
    if (addTrackCutsStr != "") {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
      for (const auto& t : addTrackCuts) {
        tempCutsStr += Form(",%s", t->GetName());
      }
    }

    // get the mc generated acceptance cut
    TString mcGenAccCutStr = fConfigCuts.MCgenAcc.value;
    if (mcGenAccCutStr != "") {
      AnalysisCut* cut = dqcuts::GetAnalysisCut(mcGenAccCutStr.Data());
      if (cut != nullptr) {
        fMCGenAccCut.AddCut(cut);
      }
      fUseMCGenAccCut = true;
    }

    // check that the barrel track cuts array required in this task is not empty
    if (!trackCutsStr.IsNull()) {
      // tokenize and loop over the barrel cuts produced by the barrel track selection task
      std::unique_ptr<TObjArray> objArray(tempCutsStr.Tokenize(","));
      fNCutsBarrel = objArray->GetEntries();
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        TString tempStr = objArray->At(icut)->GetName();
        // if the current barrel selection cut is required in this task, then switch on the corresponding bit in the mask
        // and assign histogram directories
        if (objArrayTrackCuts->FindObject(tempStr.Data()) != nullptr) {
          fTrackFilterMask |= (static_cast<uint32_t>(1) << icut);

          if (fEnableBarrelHistos) {
            // assign the pair hist directories for the current cut
            std::vector<TString> names = {
              Form("PairsBarrelSEPM_%s", objArray->At(icut)->GetName()),
              Form("PairsBarrelSEPP_%s", objArray->At(icut)->GetName()),
              Form("PairsBarrelSEMM_%s", objArray->At(icut)->GetName())};
            if (fConfigQA) {
              // assign separate hist directories for ambiguous tracks
              names.push_back(Form("PairsBarrelSEPM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelSEPP_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelSEMM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelSEPM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelSEPP_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelSEMM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
            }
            for (const auto& n : names) {
              histNames += Form("%s;", n.Data());
            }
            fTrackHistNames[icut] = names;

            // if there are pair cuts specified, assign hist directories for each barrel cut - pair cut combination
            // NOTE: This could possibly lead to large histogram outputs. It is strongly advised to use pair cuts only
            //   if you know what you are doing.
            TString cutNamesStr = fConfigCuts.pair.value;
            if (!cutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
              fNPairCuts = objArrayPair->GetEntries();
              for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
                names = {
                  Form("PairsBarrelSEPM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsBarrelSEPP_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsBarrelSEMM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName())};
                histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
                // NOTE: In the numbering scheme for the map key, we use the number of barrel cuts in the barrel-track selection task
                fTrackHistNames[fNCutsBarrel + icut * fNPairCuts + iPairCut] = names;
              } // end loop (pair cuts)
            } // end if (pair cuts)

            // assign hist directories for the MC matched pairs for each (track cut,MCsignal) combination
            if (!sigNamesStr.IsNull()) {
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
                auto sig = fRecMCSignals.at(isig);
                names = {
                  Form("PairsBarrelSEPM_%s_%s", objArray->At(icut)->GetName(), sig->GetName()),
                  Form("PairsBarrelSEPP_%s_%s", objArray->At(icut)->GetName(), sig->GetName()),
                  Form("PairsBarrelSEMM_%s_%s", objArray->At(icut)->GetName(), sig->GetName())};
                if (fConfigQA) {
                  names.push_back(Form("PairsBarrelSEPMCorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPMIncorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPM_ambiguousInBunch_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPM_ambiguousInBunchCorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPM_ambiguousInBunchIncorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPM_ambiguousOutOfBunch_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPM_ambiguousOutOfBunchCorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPM_ambiguousOutOfBunchIncorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                }
                for (const auto& n : names) {
                  histNames += Form("%s;", n.Data());
                }
                fBarrelHistNamesMCmatched.try_emplace(icut * fRecMCSignals.size() + isig, names);
              } // end loop over MC signals
            }
          } // end if enableBarrelHistos
        }
      }
    }

    // Add histogram classes for each specified MCsignal at the generator level
    // TODO: create a std::vector of hist classes to be used at Fill time, to avoid using Form in the process function
    TString sigGenNamesStr = fConfigMC.genSignals.value;
    std::unique_ptr<TObjArray> objGenSigArray(sigGenNamesStr.Tokenize(","));
    for (int isig = 0; isig < objGenSigArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objGenSigArray->At(isig)->GetName());
      if (sig) {
        fGenMCSignals.push_back(sig);
      }
    }

    // Add the MCSignals from the JSON config
    TString addMCSignalsGenStr = fConfigMC.genSignalsJSON.value;
    if (addMCSignalsGenStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsGenStr.Data());
      for (const auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() > TWO_PRONG) { // NOTE: only 2 prong signals
          continue;
        }
        fGenMCSignals.push_back(mcIt);
      }
    }

    if (isMCGen) {
      for (const auto& sig : fGenMCSignals) {
        if (sig->GetNProngs() == 1) {
          histNames += Form("MCTruthGen_%s;", sig->GetName()); // TODO: Add these names to a std::vector to avoid using Form in the process function
          histNames += Form("MCTruthGenSel_%s;", sig->GetName());
        } else if (sig->GetNProngs() == TWO_PRONG) {
          histNames += Form("MCTruthGenPair_%s;", sig->GetName());
          histNames += Form("MCTruthGenPairSel_%s;", sig->GetName());
          fHasTwoProngGenMCsignals = true;
        }
      }
    }

    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(true);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    VarManager::SetCollisionSystem((TString)fConfigOptions.collisionSystem, fConfigOptions.centerMassEnergy); // set collision system and center of mass energy

    DefineHistograms(fHistMan, histNames.Data(), fConfigAddSEPHistogram.value.data());     // define all histograms
    dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
    VarManager::SetUseVars(fHistMan->GetUsedVars());                                       // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  // Template function to run same event pairing (barrel-barrel)
  void runSameEventPairing(MyEventsVtxCovSelected const& events, PresliceUnsorted<MyBarrelAssocsPrefilter>& preslice, MyBarrelAssocsPrefilter const& assocs, MyBarrelTracksWithCovWithAmbiguities const& /*tracks*/, ReducedA3MCEvents const& /*mcEvents*/, ReducedA3MCTracks const& /*mcTracks*/)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }

    TString cutNames = fConfigCuts.track.value;
    std::map<int, std::vector<TString>> histNames = fTrackHistNames;
    std::map<int, std::vector<TString>> histNamesMC = fBarrelHistNamesMCmatched;
    int ncuts = fNCutsBarrel;

    uint32_t twoTrackFilter = static_cast<uint32_t>(0);
    int sign1 = 0;
    int sign2 = 0;
    uint32_t mcDecision = static_cast<uint32_t>(0);
    bool isCorrectAssoc_leg1 = false;
    bool isCorrectAssoc_leg2 = false;
    dielectronList.reserve(1);
    dielectronsExtraList.reserve(1);

    if (fConfigOptions.flatTables.value) {
      dielectronAllList.reserve(1);
    }

    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEventAlice3<gkEventFillMap>(event, VarManager::fgValues);
      VarManager::FillEventAlice3<VarManager::ObjTypes::ReducedEventMC>(event.reducedA3MCEvent(), VarManager::fgValues);

      auto groupedAssocs = assocs.sliceBy(preslice, event.globalIndex());
      if (groupedAssocs.size() == 0) {
        continue;
      }

      for (const auto& [a1, a2] : o2::soa::combinations(groupedAssocs, groupedAssocs)) {

        twoTrackFilter = a1.isBarrelSelected_raw() & a2.isBarrelSelected_raw() & a1.isBarrelSelectedPrefilter_raw() & a2.isBarrelSelectedPrefilter_raw() & fTrackFilterMask;

        if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
          continue;
        }

        auto t1 = a1.template reducedA3track_as<MyBarrelTracksWithCovWithAmbiguities>();
        auto t2 = a2.template reducedA3track_as<MyBarrelTracksWithCovWithAmbiguities>();
        sign1 = t1.sign();
        sign2 = t2.sign();
        // store the ambiguity number of the two dilepton legs in the last 4 digits of the two-track filter
        if (t1.barrelAmbiguityInBunch() > 1) {
          twoTrackFilter |= (static_cast<uint32_t>(1) << 28);
        }
        if (t2.barrelAmbiguityInBunch() > 1) {
          twoTrackFilter |= (static_cast<uint32_t>(1) << 29);
        }
        if (t1.barrelAmbiguityOutOfBunch() > 1) {
          twoTrackFilter |= (static_cast<uint32_t>(1) << 30);
        }
        if (t2.barrelAmbiguityOutOfBunch() > 1) {
          twoTrackFilter |= (static_cast<uint32_t>(1) << 31);
        }

        // run MC matching for this pair
        int isig = 0;
        mcDecision = 0;
        for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
          if (t1.has_reducedA3MCTrack() && t2.has_reducedA3MCTrack()) {
            if ((*sig)->CheckSignal(true, t1.reducedA3MCTrack(), t2.reducedA3MCTrack())) {
              mcDecision |= (static_cast<uint32_t>(1) << isig);
            }
          }
        } // end loop over MC signals
        if (t1.has_reducedA3MCTrack() && t2.has_reducedA3MCTrack()) {
          isCorrectAssoc_leg1 = (t1.reducedA3MCTrack().reducedA3MCEvent() == event.reducedA3MCEvent());
          isCorrectAssoc_leg2 = (t2.reducedA3MCTrack().reducedA3MCEvent() == event.reducedA3MCEvent());
        }

        VarManager::FillPairAlice3<VarManager::kDecayToEE, gkTrackFillMap>(t1, t2);
        if (fPropTrack) {
          VarManager::FillPairCollision<VarManager::kDecayToEE, gkTrackFillMap>(event, t1, t2);
        }
        /* TODO: Reimplement Pair vertexing when secondary vertexing is available
        if constexpr (TTwoProngFitter) {
          // VarManager::FillPairVertexing<VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMap>(event, t1, t2, fConfigOptions.propToPCA);
        }*/
        if (!fConfigMC.skimSignalOnly || (fConfigMC.skimSignalOnly && mcDecision > 0)) {
          dielectronList(event.globalIndex(), VarManager::fgValues[VarManager::kMass],
                         VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi],
                         t1.sign() + t2.sign(), twoTrackFilter, mcDecision);
        }

        // Fill histograms
        bool isAmbiInBunch = false;
        bool isAmbiOutOfBunch = false;

        for (int icut = 0; icut < ncuts; icut++) {
          if (twoTrackFilter & (static_cast<uint32_t>(1) << icut)) {
            isAmbiInBunch = (twoTrackFilter & (static_cast<uint32_t>(1) << 28)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 29));
            isAmbiOutOfBunch = (twoTrackFilter & (static_cast<uint32_t>(1) << 30)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 31));
            if (sign1 * sign2 < 0) {                                                    // +- pairs
              fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues); // reconstructed, unmatched
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {        // loop over MC signals
                if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                  fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][0].Data(), VarManager::fgValues); // matched signal
                  if (fConfigQA) {
                    if (isCorrectAssoc_leg1 && isCorrectAssoc_leg2) { // correct track-collision association
                      fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][3].Data(), VarManager::fgValues);
                    } else { // incorrect track-collision association
                      fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][4].Data(), VarManager::fgValues);
                    }
                    if (isAmbiInBunch) { // ambiguous in bunch
                      fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][5].Data(), VarManager::fgValues);
                      if (isCorrectAssoc_leg1 && isCorrectAssoc_leg2) {
                        fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][6].Data(), VarManager::fgValues);
                      } else {
                        fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][7].Data(), VarManager::fgValues);
                      }
                    }
                    if (isAmbiOutOfBunch) { // ambiguous out of bunch
                      fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][8].Data(), VarManager::fgValues);
                      if (isCorrectAssoc_leg1 && isCorrectAssoc_leg2) {
                        fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][9].Data(), VarManager::fgValues);
                      } else {
                        fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][10].Data(), VarManager::fgValues);
                      }
                    }
                  }
                }
                if (fConfigQA) {
                  if (isAmbiInBunch) {
                    fHistMan->FillHistClass(histNames[icut][3].Data(), VarManager::fgValues);
                  }
                  if (isAmbiOutOfBunch) {
                    fHistMan->FillHistClass(histNames[icut][3 + 3].Data(), VarManager::fgValues);
                  }
                }
              }
            } else {
              if (sign1 > 0) { // ++ pairs
                fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
                for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                  if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                    fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][1].Data(), VarManager::fgValues);
                  }
                }
                if (fConfigQA) {
                  if (isAmbiInBunch) {
                    fHistMan->FillHistClass(histNames[icut][4].Data(), VarManager::fgValues);
                  }
                  if (isAmbiOutOfBunch) {
                    fHistMan->FillHistClass(histNames[icut][4 + 3].Data(), VarManager::fgValues);
                  }
                }
              } else { // -- pairs
                fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
                for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                  if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                    fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][2].Data(), VarManager::fgValues);
                  }
                }
                if (fConfigQA) {
                  if (isAmbiInBunch) {
                    fHistMan->FillHistClass(histNames[icut][5].Data(), VarManager::fgValues);
                  }
                  if (isAmbiOutOfBunch) {
                    fHistMan->FillHistClass(histNames[icut][5 + 3].Data(), VarManager::fgValues);
                  }
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
      } // end loop over pairs of track associations
    } // end loop over events
  }

  PresliceUnsorted<ReducedA3MCTracks> perReducedMcEvent = aod::reducedA3trackMC::reducedA3MCEventId;

  void runMCGenWithGrouping(MyEventsVtxCovSelected const& events, ReducedA3MCEvents const& /*mcEvents*/, ReducedA3MCTracks const& mcTracks)
  {
    [[maybe_unused]] uint32_t mcDecision = 0;
    int isig = 0;

    for (const auto& mctrack : mcTracks) {
      VarManager::FillTrackMC(mcTracks, mctrack);
      // if we have a mc generated acceptance cut, apply it here
      if (fUseMCGenAccCut) {
        if (!fMCGenAccCut.IsSelected(VarManager::fgValues)) {
          continue;
        }
      }
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (const auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), VarManager::fgValues);
        }
      }
    }
    // Fill Generated histograms taking into account selected collisions
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedA3MCEvent()) {
        continue;
      }

      for (const auto& track : mcTracks) {
        if (track.reducedA3MCEventId() != event.reducedA3MCEventId()) {
          continue;
        }
        VarManager::FillTrackMC(mcTracks, track);
        // if we have a mc generated acceptance cut, apply it here
        if (fUseMCGenAccCut) {
          if (!fMCGenAccCut.IsSelected(VarManager::fgValues)) {
            continue;
          }
        }
        auto track_raw = mcTracks.rawIteratorAt(track.globalIndex());
        mcDecision = 0;
        isig = 0;
        for (const auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, track_raw)) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), VarManager::fgValues);
            MCTruthTableEffi(VarManager::fgValues[VarManager::kMCPt], VarManager::fgValues[VarManager::kMCEta], VarManager::fgValues[VarManager::kMCY], VarManager::fgValues[VarManager::kMCPhi], VarManager::fgValues[VarManager::kMCVz], VarManager::fgValues[VarManager::kMCVtxZ], VarManager::fgValues[VarManager::kMultFT0A], VarManager::fgValues[VarManager::kMultFT0C], VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kVtxNcontribReal]);
          }
          isig++;
        }
      }
    } // end loop over reconstructed events

    if (fHasTwoProngGenMCsignals) {
      for (const auto& [t1, t2] : combinations(mcTracks, mcTracks)) {
        auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
        auto t2_raw = mcTracks.rawIteratorAt(t2.globalIndex());
        if (t1_raw.reducedA3MCEventId() == t2_raw.reducedA3MCEventId()) {
          for (const auto& sig : fGenMCSignals) {
            if (sig->GetNProngs() != TWO_PRONG) { // NOTE: 2-prong signals required here
              continue;
            }
            if (sig->CheckSignal(true, t1_raw, t2_raw)) {
              VarManager::FillPairMC<VarManager::kDecayToEE>(t1, t2);
              if (fUseMCGenAccCut) {
                if (!fMCGenAccCut.IsSelected(VarManager::fgValues)) {
                  continue;
                }
              }
              fHistMan->FillHistClass(Form("MCTruthGenPair_%s", sig->GetName()), VarManager::fgValues);
            }
          }
        }
      }
    }
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedA3MCEvent()) {
        continue;
      }
      // CURRENTLY ONLY FOR 1-GENERATION 2-PRONG SIGNALS
      if (fHasTwoProngGenMCsignals) {
        auto groupedMCTracks = mcTracks.sliceBy(perReducedMcEvent, event.reducedA3MCEventId());
        groupedMCTracks.bindInternalIndicesTo(&mcTracks);
        for (const auto& [t1, t2] : combinations(groupedMCTracks, groupedMCTracks)) {
          auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
          auto t2_raw = mcTracks.rawIteratorAt(t2.globalIndex());
          if (t1_raw.reducedA3MCEventId() == t2_raw.reducedA3MCEventId()) {
            mcDecision = 0;
            isig = 0;
            for (const auto& sig : fGenMCSignals) {
              if (sig->GetNProngs() != TWO_PRONG) { // NOTE: 2-prong signals required here
                continue;
              }
              if (sig->CheckSignal(true, t1_raw, t2_raw)) {
                mcDecision |= (static_cast<uint32_t>(1) << isig);
                VarManager::FillPairMC<VarManager::kDecayToEE>(t1, t2);
                if (fUseMCGenAccCut) {
                  if (!fMCGenAccCut.IsSelected(VarManager::fgValues)) {
                    continue;
                  }
                }
                fHistMan->FillHistClass(Form("MCTruthGenPairSel_%s", sig->GetName()), VarManager::fgValues);
              }
              isig++;
            }
          }
        }
      } // end loop over reconstructed events
    }
  }

  void processBarrelOnlySkimmed(MyEventsVtxCovSelected const& events,
                                MyBarrelAssocsPrefilter const& barrelAssocs,
                                MyBarrelTracksWithCovWithAmbiguities const& barrelTracks, ReducedA3MCEvents const& mcEvents, ReducedA3MCTracks const& mcTracks)
  {
    runSameEventPairing(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks);
    runMCGenWithGrouping(events, mcEvents, mcTracks);
  }

  PresliceUnsorted<ReducedA3MCTracks> perReducedMcGenEvent = aod::reducedA3trackMC::reducedA3MCEventId;

  void processMCGen(soa::Filtered<MyEventsVtxCovSelected> const& events, ReducedA3MCEvents const& /*mcEvents*/, ReducedA3MCTracks const& mcTracks)
  {
    // Fill Generated histograms taking into account all generated tracks
    [[maybe_unused]] uint32_t mcDecision = 0;
    int isig = 0;

    for (const auto& mctrack : mcTracks) {
      VarManager::FillTrackMC(mcTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (const auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), VarManager::fgValues);
        }
      }
    }

    // Fill Generated histograms taking into account selected collisions
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedA3MCEvent()) {
        continue;
      }
      VarManager::FillEventAlice3<gkEventFillMap>(event, VarManager::fgValues);
      VarManager::FillEventAlice3<VarManager::ObjTypes::ReducedEventMC>(event.reducedA3MCEvent(), VarManager::fgValues);
      // auto groupedMCTracks = mcTracks.sliceBy(perReducedMcGenEvent, event.reducedA3MCEventId());
      // groupedMCTracks.bindInternalIndicesTo(&mcTracks);
      // for (const auto& track : groupedMCTracks) {
      for (const auto& track : mcTracks) {
        if (track.reducedA3MCEventId() != event.reducedA3MCEventId()) {
          continue;
        }
        VarManager::FillTrackMC(mcTracks, track);
        auto track_raw = mcTracks.rawIteratorAt(track.globalIndex());
        // auto track_raw = groupedMCTracks.rawIteratorAt(track.globalIndex());
        mcDecision = 0;
        isig = 0;
        for (const auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, track_raw)) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), VarManager::fgValues);
            MCTruthTableEffi(VarManager::fgValues[VarManager::kMCPt], VarManager::fgValues[VarManager::kMCEta], VarManager::fgValues[VarManager::kMCY], VarManager::fgValues[VarManager::kMCPhi], VarManager::fgValues[VarManager::kMCVz], VarManager::fgValues[VarManager::kMCVtxZ], VarManager::fgValues[VarManager::kMultFT0A], VarManager::fgValues[VarManager::kMultFT0C], VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kVtxNcontribReal]);
          }
          isig++;
        }
      }
    } // end loop over reconstructed events
    if (fHasTwoProngGenMCsignals) {
      for (const auto& [t1, t2] : combinations(mcTracks, mcTracks)) {
        auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
        auto t2_raw = mcTracks.rawIteratorAt(t2.globalIndex());
        if (t1_raw.reducedA3MCEventId() == t2_raw.reducedA3MCEventId()) {
          for (const auto& sig : fGenMCSignals) {
            if (sig->GetNProngs() != TWO_PRONG) { // NOTE: 2-prong signals required here
              continue;
            }
            if (sig->CheckSignal(true, t1_raw, t2_raw)) {
              fHistMan->FillHistClass(Form("MCTruthGenPair_%s", sig->GetName()), VarManager::fgValues);
            }
          }
        }
      }
    }
    // Fill Generated PAIR histograms taking into account selected collisions
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedA3MCEvent()) {
        continue;
      }

      if (fHasTwoProngGenMCsignals) {
        for (const auto& [t1, t2] : combinations(mcTracks, mcTracks)) {
          if (t1.reducedA3MCEventId() != event.reducedA3MCEventId()) {
            continue;
          }
          if (t2.reducedA3MCEventId() != event.reducedA3MCEventId()) {
            continue;
          }
          auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
          auto t2_raw = mcTracks.rawIteratorAt(t2.globalIndex());
          if (t1_raw.reducedA3MCEventId() == t2_raw.reducedA3MCEventId()) {
            mcDecision = 0;
            isig = 0;
            for (const auto& sig : fGenMCSignals) {
              if (sig->GetNProngs() != TWO_PRONG) { // NOTE: 2-prong signals required here
                continue;
              }
              if (sig->CheckSignal(true, t1_raw, t2_raw)) {
                mcDecision |= (static_cast<uint32_t>(1) << isig);
                fHistMan->FillHistClass(Form("MCTruthGenPairSel_%s", sig->GetName()), VarManager::fgValues);
              }
              isig++;
            }
          }
        }
      }
    } // end loop over reconstructed events
  }

  void processMCGenWithGrouping(soa::Filtered<MyEventsVtxCovSelected> const& events, ReducedA3MCEvents const& /*mcEvents*/, ReducedA3MCTracks const& mcTracks)
  {
    [[maybe_unused]] uint32_t mcDecision = 0;
    int isig = 0;

    for (const auto& mctrack : mcTracks) {
      VarManager::FillTrackMC(mcTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (const auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), VarManager::fgValues);
        }
      }
    }
    // Fill Generated histograms taking into account selected collisions
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedA3MCEvent()) {
        continue;
      }

      for (const auto& track : mcTracks) {
        if (track.reducedA3MCEventId() != event.reducedA3MCEventId()) {
          continue;
        }
        VarManager::FillTrackMC(mcTracks, track);
        auto track_raw = mcTracks.rawIteratorAt(track.globalIndex());
        mcDecision = 0;
        isig = 0;
        for (const auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, track_raw)) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), VarManager::fgValues);
            MCTruthTableEffi(VarManager::fgValues[VarManager::kMCPt], VarManager::fgValues[VarManager::kMCEta], VarManager::fgValues[VarManager::kMCY], VarManager::fgValues[VarManager::kMCPhi], VarManager::fgValues[VarManager::kMCVz], VarManager::fgValues[VarManager::kMCVtxZ], VarManager::fgValues[VarManager::kMultFT0A], VarManager::fgValues[VarManager::kMultFT0C], VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kVtxNcontribReal]);
          }
          isig++;
        }
      }
    } // end loop over reconstructed events
    if (fHasTwoProngGenMCsignals) {
      for (const auto& [t1, t2] : combinations(mcTracks, mcTracks)) {
        auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
        auto t2_raw = mcTracks.rawIteratorAt(t2.globalIndex());
        if (t1_raw.reducedA3MCEventId() == t2_raw.reducedA3MCEventId()) {
          for (const auto& sig : fGenMCSignals) {
            if (sig->GetNProngs() != TWO_PRONG) { // NOTE: 2-prong signals required here
              continue;
            }
            if (sig->CheckSignal(true, t1_raw, t2_raw)) {
              fHistMan->FillHistClass(Form("MCTruthGenPair_%s", sig->GetName()), VarManager::fgValues);
            }
          }
        }
      }
    }
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedA3MCEvent()) {
        continue;
      }
      // CURRENTLY ONLY FOR 1-GENERATION 2-PRONG SIGNALS
      if (fHasTwoProngGenMCsignals) {
        auto groupedMCTracks = mcTracks.sliceBy(perReducedMcEvent, event.reducedA3MCEventId());
        groupedMCTracks.bindInternalIndicesTo(&mcTracks);
        for (const auto& [t1, t2] : combinations(groupedMCTracks, groupedMCTracks)) {
          auto t1_raw = groupedMCTracks.rawIteratorAt(t1.globalIndex());
          auto t2_raw = groupedMCTracks.rawIteratorAt(t2.globalIndex());
          if (t1_raw.reducedA3MCEventId() == t2_raw.reducedA3MCEventId()) {
            mcDecision = 0;
            isig = 0;
            for (const auto& sig : fGenMCSignals) {
              if (sig->GetNProngs() != TWO_PRONG) { // NOTE: 2-prong signals required here
                continue;
              }
              if (sig->CheckSignal(true, t1_raw, t2_raw)) {
                mcDecision |= (static_cast<uint32_t>(1) << isig);
                fHistMan->FillHistClass(Form("MCTruthGenPairSel_%s", sig->GetName()), VarManager::fgValues);
              }
              isig++;
            }
          }
        }
      } // end loop over reconstructed events
    }
  }

  void processDummy(MyEvents const&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlySkimmed, "Run barrel only pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMCGen, "Loop over MC particle stack and fill generator level histograms", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMCGenWithGrouping, "Loop over MC particle stack (grouped MCTracks) and fill generator level histograms", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", true);
};

struct AnalysisAsymmetricPairing {

  Produces<aod::Ditracks> ditrackList;
  Produces<aod::DitracksExtra> ditrackExtraList;

  // Output objects
  OutputObj<THashList> fOutputList{"output"};

  // Configurables
  Configurable<std::string> fConfigLegCuts{"cfgLegCuts", "", "<leg-A-1>:<leg-B-1>[:<leg-C-1>],[<leg-A-2>:<leg-B-2>[:<leg-C-1>],...]"};
  Configurable<uint32_t> fConfigLegAFilterMask{"cfgLegAFilterMask", 0, "Filter mask corresponding to cuts in event-selection"};
  Configurable<uint32_t> fConfigLegBFilterMask{"cfgLegBFilterMask", 0, "Filter mask corresponding to cuts in event-selection"};
  Configurable<uint32_t> fConfigLegCFilterMask{"cfgLegCFilterMask", 0, "Filter mask corresponding to cuts in event-selection"};
  Configurable<std::string> fConfigCommonTrackCuts{"cfgCommonTrackCuts", "", "Comma separated list of cuts to be applied to all legs"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "", "Comma separated list of pair cuts"};
  Configurable<std::string> fConfigPairCutsJSON{"cfgPairCutsJSON", "", "Additional list of pair cuts in JSON format"};
  Configurable<bool> fConfigSkipAmbiguousIdCombinations{"cfgSkipAmbiguousIdCombinations", true, "Choose whether to skip pairs/triples which pass a stricter combination of cuts, e.g. KKPi triplets for D+ -> KPiPi"};

  Configurable<std::string> fConfigHistogramSubgroups{"cfgAsymmetricPairingHistogramsSubgroups", "barrel,vertexing", "Comma separated list of asymmetric-pairing histogram subgroups"};
  Configurable<bool> fConfigSameSignHistograms{"cfgSameSignHistograms", false, "Include same sign pair histograms for 2-prong decays"};
  Configurable<bool> fConfigReflectedHistograms{"cfgReflectedHistograms", false, "Include separate histograms for pairs which are reflections of previously counted pairs"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};

  Configurable<std::string> fConfigMCGenSignals{"cfgBarrelMCGenSignals", "", "Comma separated list of MC signals (generated)"};
  Configurable<std::string> fConfigMCRecSignals{"cfgBarrelMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};
  Configurable<std::string> fConfigMCRecSignalsJSON{"cfgMCRecSignalsJSON", "", "Additional list of MC signals (reconstructed) via JSON"};
  Configurable<std::string> fConfigMCGenSignalsJSON{"cfgMCGenSignalsJSON", "", "Comma separated list of MC signals (generated) via JSON"};

  HistogramManager* fHistMan;

  std::vector<AnalysisCompositeCut*> fPairCuts;
  int fNPairHistPrefixes;

  std::vector<MCSignal*> fRecMCSignals;
  std::vector<MCSignal*> fGenMCSignals;

  // Filter masks to find legs in BarrelTrackCuts table
  uint32_t fLegAFilterMask;
  uint32_t fLegBFilterMask;
  uint32_t fLegCFilterMask;
  // Maps tracking which combination of leg cuts the track cuts participate in
  std::map<int, uint32_t> fConstructedLegAFilterMasksMap;
  std::map<int, uint32_t> fConstructedLegBFilterMasksMap;
  std::map<int, uint32_t> fConstructedLegCFilterMasksMap;
  // Filter map for common track cuts
  uint32_t fCommonTrackCutMask;
  // Map tracking which common track cut the track cuts correspond to
  std::map<int, uint32_t> fCommonTrackCutFilterMasks;

  int fNLegCuts;
  int fNPairCuts = 0;
  int fNCommonTrackCuts;
  // vectors for cut names and signal names, for easy access when calling FillHistogramList()
  std::vector<TString> fLegCutNames;
  std::vector<TString> fPairCutNames;
  std::vector<TString> fCommonCutNames;
  std::vector<TString> fRecMCSignalNames;

  Filter eventFilter = aod::dqanalysisflags::isEventSelected > static_cast<uint32_t>(0);

  PresliceUnsorted<MyBarrelAssocs> trackAssocsPerCollision = aod::reducedA3track_association::reducedA3eventId;
  // PresliceUnsorted<aod::ReducedA3TracksAssoc> trackAssocsPerCollision = aod::reducedA3track_association::reducedA3eventId;

  // Partitions for triplets and asymmetric pairs
  Partition<MyBarrelAssocs> legACandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & fConfigLegAFilterMask) > static_cast<uint32_t>(0);
  Partition<MyBarrelAssocs> legBCandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & fConfigLegBFilterMask) > static_cast<uint32_t>(0);
  Partition<MyBarrelAssocs> legCCandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & fConfigLegCFilterMask) > static_cast<uint32_t>(0);

  // Map to track how many times a pair of tracks has been encountered
  std::map<std::pair<int32_t, int32_t>, int8_t> fPairCount;

  void init(o2::framework::InitContext& context)
  {
    bool isMCGen = context.mOptions.get<bool>("processMCGen") || context.mOptions.get<bool>("processMCGenWithEventSelection");
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(true);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Get the leg cut filter masks
    fLegAFilterMask = fConfigLegAFilterMask.value;
    fLegBFilterMask = fConfigLegBFilterMask.value;
    fLegCFilterMask = fConfigLegCFilterMask.value;

    // Get the pair cuts
    TString cutNamesStr = fConfigPairCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // Extra pair cuts via JSON
    TString addPairCutsStr = fConfigPairCutsJSON.value;
    if (addPairCutsStr != "") {
      std::vector<AnalysisCut*> addPairCuts = dqcuts::GetCutsFromJSON(addPairCutsStr.Data());
      for (const auto& t : addPairCuts) {
        fPairCuts.push_back(reinterpret_cast<AnalysisCompositeCut*>(t));
        cutNamesStr += Form(",%s", t->GetName());
      }
    }
    std::unique_ptr<TObjArray> objArrayPairCuts(cutNamesStr.Tokenize(","));
    fNPairCuts = objArrayPairCuts->GetEntries();
    for (int j = 0; j < fNPairCuts; j++) {
      fPairCutNames.push_back(objArrayPairCuts->At(j)->GetName());
    }

    // Setting the MC rec signal names
    TString sigNamesStr = fConfigMCRecSignals.value;
    std::unique_ptr<TObjArray> objRecSigArray(sigNamesStr.Tokenize(","));
    for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
      if (sig) {
        fRecMCSignals.push_back(sig);
      }
    }
    // Add the reco MCSignals from the JSON config
    TString addMCSignalsStr = fConfigMCRecSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (const auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() != TWO_PRONG && mcIt->GetNProngs() != THREE_PRONG) {
          LOG(fatal) << "Signal at reconstructed level requested (" << mcIt->GetName() << ") " << "does not have 2 or 3 prongs! Fix it";
        }
        fRecMCSignals.push_back(mcIt);
        sigNamesStr += Form(",%s", mcIt->GetName());
      }
    }
    // Put all the reco MCSignal names in the vector for histogram naming
    std::unique_ptr<TObjArray> objArrayRecMCSignals(sigNamesStr.Tokenize(","));
    for (int i = 0; i < objArrayRecMCSignals->GetEntries(); i++) {
      fRecMCSignalNames.push_back(objArrayRecMCSignals->At(i)->GetName());
    }

    // Get the barrel track selection cuts
    string tempCuts;
    getTaskOptionValue<string>(context, "analysis-track-selection", "cfgTrackCuts", tempCuts, false);
    TString tempCutsStr = tempCuts;
    // check also the cuts added via JSON and add them to the string of cuts
    getTaskOptionValue<string>(context, "analysis-track-selection", "cfgBarrelTrackCutsJSON", tempCuts, false);
    TString addTrackCutsStr = tempCuts;
    if (addTrackCutsStr != "") {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
      for (const auto& t : addTrackCuts) {
        tempCutsStr += Form(",%s", t->GetName());
      }
    }
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
          fCommonTrackCutMask |= static_cast<uint32_t>(1) << objArray->IndexOf(objArrayCommon->At(icut));
          fCommonTrackCutFilterMasks[icut] = static_cast<uint32_t>(1) << objArray->IndexOf(objArrayCommon->At(icut));
          fCommonCutNames.push_back(objArrayCommon->At(icut)->GetName());
        } else {
          LOGF(fatal, "Common track cut %s was not calculated upstream. Check the config!", objArrayCommon->At(icut)->GetName());
        }
      }
    }
    // Check that the leg cut masks make sense
    if (static_cast<int>(std::floor(std::log2(fLegAFilterMask))) + 1 > objArray->GetEntries()) {
      LOGF(fatal, "fConfigLegAFilterMask has highest bit at position %d, but track-selection only has %d cuts!", static_cast<int>(std::floor(std::log2(fLegAFilterMask))) + 1, objArray->GetEntries());
    }
    if (static_cast<int>(std::floor(std::log2(fLegBFilterMask))) + 1 > objArray->GetEntries()) {
      LOGF(fatal, "fConfigLegBFilterMask has highest bit at position %d, but track-selection only has %d cuts!", static_cast<int>(std::floor(std::log2(fLegBFilterMask))) + 1, objArray->GetEntries());
    }
    if (static_cast<int>(std::floor(std::log2(fLegCFilterMask))) + 1 > objArray->GetEntries()) {
      LOGF(fatal, "fConfigLegCFilterMask has highest bit at position %d, but track-selection only has %d cuts!", static_cast<int>(std::floor(std::log2(fLegCFilterMask))) + 1, objArray->GetEntries());
    }

    // Get the cuts defining the legs
    uint32_t fConstructedLegAFilterMask = 0;
    uint32_t fConstructedLegBFilterMask = 0;
    uint32_t fConstructedLegCFilterMask = 0;
    TString legCutsStr = fConfigLegCuts.value;
    std::unique_ptr<TObjArray> objArrayLegs(legCutsStr.Tokenize(","));
    if (objArrayLegs->GetEntries() == 0 && !isMCGen) {
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
      if (legs->GetEntries() == THREE_PRONG) {
        isThreeProng.push_back(true);
      } else if (legs->GetEntries() == TWO_PRONG) {
        isThreeProng.push_back(false);
      } else {
        LOGF(fatal, "Leg cuts %s has the wrong format and could not be parsed!", legsStr.Data());
        continue;
      }
      // Find leg cuts in the track selection cuts
      legAIdx = objArray->IndexOf(legs->At(0));
      if (legAIdx >= 0) {
        fConstructedLegAFilterMask |= static_cast<uint32_t>(1) << legAIdx;
        fConstructedLegAFilterMasksMap[icut] |= static_cast<uint32_t>(1) << legAIdx;
      } else {
        LOGF(fatal, "Leg A cut %s was not calculated upstream. Check the config!", legs->At(0)->GetName());
        continue;
      }
      legBIdx = objArray->IndexOf(legs->At(1));
      if (legBIdx >= 0) {
        fConstructedLegBFilterMask |= static_cast<uint32_t>(1) << legBIdx;
        fConstructedLegBFilterMasksMap[icut] |= static_cast<uint32_t>(1) << legBIdx;
      } else {
        LOGF(fatal, "Leg B cut %s was not calculated upstream. Check the config!", legs->At(1)->GetName());
        continue;
      }
      if (isThreeProng[icut]) {
        legCIdx = objArray->IndexOf(legs->At(2));
        if (legCIdx >= 0) {
          fConstructedLegCFilterMask |= static_cast<uint32_t>(1) << legCIdx;
          fConstructedLegCFilterMasksMap[icut] |= static_cast<uint32_t>(1) << legCIdx;
        } else {
          LOGF(fatal, "Leg C cut %s was not calculated upstream. Check the config!", legs->At(2)->GetName());
          continue;
        }
      }
      // Leg cut config is fine, store the leg cut name in a vector
      fLegCutNames.push_back(legsStr);

      // Define histogram and histogram directory names
      if (isThreeProng[icut]) {
        DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s", legsStr.Data()), fConfigHistogramSubgroups.value.data());
        if (fConfigQA) {
          DefineHistograms(fHistMan, Form("TripletsBarrelSE_ambiguous_%s", legsStr.Data()), fConfigHistogramSubgroups.value.data());
        }

        std::unique_ptr<TObjArray> objArrayCommon(commonNamesStr.Tokenize(","));
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
          DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName()), fConfigHistogramSubgroups.value.data());
        }

        TString cutNamesStr = fConfigPairCuts.value;
        if (!cutNamesStr.IsNull()) { // if pair cuts
          std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
          fNPairCuts = objArrayPair->GetEntries();
          for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
            DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s", legsStr.Data(), objArrayPair->At(iPairCut)->GetName()), fConfigHistogramSubgroups.value.data());
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName()), fConfigHistogramSubgroups.value.data());
            } // end loop (common cuts)
          } // end loop (pair cuts)
        } // end if (pair cuts)

        // TODO: assign hist directories for the MC matched triplets for each (leg cut combo,MCsignal) combination
        if (!sigNamesStr.IsNull()) {
          for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
            auto sig = fRecMCSignals.at(isig);
            DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s", legsStr.Data(), sig->GetName()), fConfigHistogramSubgroups.value.data());

            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), sig->GetName()), fConfigHistogramSubgroups.value.data());
            }

            if (!cutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
              for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
                DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s_%s", legsStr.Data(), objArrayPair->At(iPairCut)->GetName(), sig->GetName()), fConfigHistogramSubgroups.value.data());
                for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                  DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName(), sig->GetName()), fConfigHistogramSubgroups.value.data());
                } // end loop (common cuts)
              } // end loop (pair cuts)
            } // end if (pair cuts)
          } // end loop over MC signals
        } // end if (MC signals)
      } else {
        std::vector<TString> pairHistPrefixes = {"PairsBarrelSEPM"};
        if (fConfigSameSignHistograms.value) {
          pairHistPrefixes.push_back("PairsBarrelSEPP");
          pairHistPrefixes.push_back("PairsBarrelSEMM");
        }
        fNPairHistPrefixes = pairHistPrefixes.size();

        for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
          DefineHistograms(fHistMan, Form("%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data()), fConfigHistogramSubgroups.value.data());
        }
        if (fConfigQA) {
          for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
            DefineHistograms(fHistMan, Form("%s_ambiguous_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data()), fConfigHistogramSubgroups.value.data());
          }
        }
        if (fConfigReflectedHistograms.value) {
          for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
            DefineHistograms(fHistMan, Form("%s_reflected_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data()), fConfigHistogramSubgroups.value.data());
          }
        }

        std::unique_ptr<TObjArray> objArrayCommon(commonNamesStr.Tokenize(","));
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
          for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
            DefineHistograms(fHistMan, Form("%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName()), fConfigHistogramSubgroups.value.data());
          }
        }

        if (!cutNamesStr.IsNull()) { // if pair cuts
          std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
          fNPairCuts = objArrayPair->GetEntries();
          for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
            for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
              DefineHistograms(fHistMan, Form("%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayPair->At(iPairCut)->GetName()), fConfigHistogramSubgroups.value.data());
            }
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                DefineHistograms(fHistMan, Form("%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName()), fConfigHistogramSubgroups.value.data());
              }
            } // end loop (common cuts)
          } // end loop (pair cuts)
        } // end if (pair cuts)

        // assign hist directories for the MC matched triplets for each (leg cut combo,MCsignal) combination
        if (!sigNamesStr.IsNull()) {
          for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
            auto sig = fRecMCSignals.at(isig);
            for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
              DefineHistograms(fHistMan, Form("%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), sig->GetName()), fConfigHistogramSubgroups.value.data());
            }
            if (fConfigReflectedHistograms.value) {
              for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                DefineHistograms(fHistMan, Form("%s_reflected_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), sig->GetName()), fConfigHistogramSubgroups.value.data());
              }
            }

            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                DefineHistograms(fHistMan, Form("%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), sig->GetName()), fConfigHistogramSubgroups.value.data());
              }
            }

            if (!cutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
              for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
                for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                  DefineHistograms(fHistMan, Form("%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayPair->At(iPairCut)->GetName(), sig->GetName()), fConfigHistogramSubgroups.value.data());
                }
                for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                  for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                    DefineHistograms(fHistMan, Form("%s_%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName(), sig->GetName()), fConfigHistogramSubgroups.value.data());
                  }
                } // end loop (common cuts)
              } // end loop (pair cuts)
            } // end if (pair cuts)
          } // end loop over MC signals
        } // end if (MC signals)
      }
    }

    // Add histogram classes for each specified MCsignal at the generator level
    // TODO: create a std::vector of hist classes to be used at Fill time, to avoid using Form in the process function
    TString sigGenNamesStr = fConfigMCGenSignals.value;
    std::unique_ptr<TObjArray> objGenSigArray(sigGenNamesStr.Tokenize(","));
    for (int isig = 0; isig < objGenSigArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objGenSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() == 1) { // NOTE: 1-prong signals required
          fGenMCSignals.push_back(sig);
          DefineHistograms(fHistMan, Form("MCTruthGen_%s;", sig->GetName()), fConfigHistogramSubgroups.value.data());    // TODO: Add these names to a std::vector to avoid using Form in the process function
          DefineHistograms(fHistMan, Form("MCTruthGenSel_%s;", sig->GetName()), fConfigHistogramSubgroups.value.data()); // TODO: Add these names to a std::vector to avoid using Form in the process function
        }
      }
    }

    // Add the gen MCSignals from the JSON config
    addMCSignalsStr = fConfigMCGenSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (const auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() == 1) {
          fGenMCSignals.push_back(mcIt);
          DefineHistograms(fHistMan, Form("MCTruthGen_%s;", mcIt->GetName()), fConfigHistogramSubgroups.value.data());    // TODO: Add these names to a std::vector to avoid using Form in the process function
          DefineHistograms(fHistMan, Form("MCTruthGenSel_%s;", mcIt->GetName()), fConfigHistogramSubgroups.value.data()); // TODO: Add these names to a std::vector to avoid using Form in the process function
        }
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

    dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
    VarManager::SetUseVars(fHistMan->GetUsedVars());                                       // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  // Function to run same event pairing with asymmetric pairs (e.g. kaon-pion)
  void runAsymmetricPairing(MyEventsVtxCovSelected const& events, PresliceUnsorted<MyBarrelAssocs>& preslice, MyBarrelAssocs const& /*assocs*/, MyBarrelTracksWithCovWithAmbiguities const& /*tracks*/, ReducedA3MCEvents const& /*mcEvents*/, ReducedA3MCTracks const& /*mcTracks*/)
  {
    fPairCount.clear();

    int sign1 = 0;
    int sign2 = 0;
    uint32_t mcDecision = 0;
    ditrackList.reserve(1);
    ditrackExtraList.reserve(1);

    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEventAlice3<gkEventFillMapWithCov>(event, VarManager::fgValues);

      auto groupedLegAAssocs = legACandidateAssocs.sliceBy(preslice, event.globalIndex());
      if (groupedLegAAssocs.size() == 0) {
        continue;
      }
      auto groupedLegBAssocs = legBCandidateAssocs.sliceBy(preslice, event.globalIndex());
      if (groupedLegBAssocs.size() == 0) {
        continue;
      }

      for (const auto& [a1, a2] : combinations(soa::CombinationsFullIndexPolicy(groupedLegAAssocs, groupedLegBAssocs))) {

        uint32_t twoTrackFilter = 0;
        uint32_t twoTrackCommonFilter = 0;
        uint32_t pairFilter = 0;
        bool isPairIdWrong = false;
        for (int icut = 0; icut < fNLegCuts; ++icut) {
          // Find leg pair definitions both candidates participate in
          if ((a1.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut]) && (a2.isBarrelSelected_raw() & fConstructedLegBFilterMasksMap[icut])) {
            twoTrackFilter |= static_cast<uint32_t>(1) << icut;
            // If the supposed pion passes a kaon cut, this is a K+K-. Skip it.
            if (fConfigSkipAmbiguousIdCombinations.value) {
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

        auto t1 = a1.template reducedA3track_as<MyBarrelTracksWithCovWithAmbiguities>();
        auto t2 = a2.template reducedA3track_as<MyBarrelTracksWithCovWithAmbiguities>();

        // Avoid self-pairs
        if (t1.globalIndex() == t2.globalIndex()) {
          continue;
        }

        bool isReflected = false;
        std::pair<int32_t, int32_t> trackIds(t1.globalIndex(), t2.globalIndex());
        if (fPairCount.find(trackIds) != fPairCount.end()) {
          // Double counting is possible due to track-collision ambiguity. Skip pairs which were counted before
          fPairCount[trackIds] += 1;
          continue;
        }
        if (fPairCount.find(std::pair(trackIds.second, trackIds.first)) != fPairCount.end()) {
          isReflected = true;
        }
        fPairCount[trackIds] += 1;

        sign1 = t1.sign();
        sign2 = t2.sign();
        // store the ambiguity number of the two dilepton legs in the last 4 digits of the two-track filter
        if (t1.barrelAmbiguityInBunch() > 1 || t1.barrelAmbiguityOutOfBunch() > 1) {
          twoTrackFilter |= static_cast<uint32_t>(1) << 30;
        }
        if (t2.barrelAmbiguityInBunch() > 1 || t2.barrelAmbiguityOutOfBunch() > 1) {
          twoTrackFilter |= static_cast<uint32_t>(1) << 31;
        }

        // run MC matching for this pair
        int isig = 0;
        mcDecision = 0;
        for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
          if (t1.has_reducedA3MCTrack() && t2.has_reducedA3MCTrack()) {
            VarManager::FillPairMC<VarManager::kDecayToKPi>(t1.reducedA3MCTrack(), t2.reducedA3MCTrack());
            if ((*sig)->CheckSignal(true, t1.reducedA3MCTrack(), t2.reducedA3MCTrack())) {
              mcDecision |= static_cast<uint32_t>(1) << isig;
            }
          }
        } // end loop over MC signals

        VarManager::FillPairAlice3<VarManager::kDecayToKPi, gkTrackFillMapWithCov>(t1, t2);
        /*TODO: Reimplement when secondary vertexing is available
        if constexpr (TTwoProngFitter) {
          VarManager::FillPairVertexing<VarManager::kDecayToKPi, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, t1, t2, fConfigPropToPCA);
        }*/

        // Fill histograms
        bool isAmbi = false;
        for (int icut = 0; icut < fNLegCuts; icut++) {
          if (twoTrackFilter & (static_cast<uint32_t>(1) << icut)) {
            isAmbi = (twoTrackFilter & (static_cast<uint32_t>(1) << 30)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 31));
            if (sign1 * sign2 < 0) {                                                                                // +- pairs
              fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s", fLegCutNames[icut].Data()), VarManager::fgValues); // reconstructed, unmatched
              if (isAmbi && fConfigQA) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPM_ambiguous_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
              }
              if (isReflected && fConfigReflectedHistograms.value) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPM_reflected_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
              }
            } else if (fConfigSameSignHistograms.value) {
              if (sign1 > 0) { // ++ pairs
                fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                if (isAmbi && fConfigQA) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPP_ambiguous_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                }
                if (isReflected && fConfigReflectedHistograms.value) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPP_reflected_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                }
              } else { // -- pairs
                fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                if (isAmbi && fConfigQA) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEMM_ambiguous_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                }
                if (isReflected && fConfigReflectedHistograms) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEMM_reflected_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                }
              }
            }
            for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
              if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                if (sign1 * sign2 < 0) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                  if (isReflected && fConfigReflectedHistograms.value) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPM_reflected_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                  }
                } else if (fConfigSameSignHistograms.value) {
                  if (sign1 > 0) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    if (isReflected && fConfigReflectedHistograms.value) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEPP_reflected_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    }
                  } else {
                    fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    if (isReflected && fConfigReflectedHistograms.value) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEMM_reflected_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    }
                  }
                }
              }
            }
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
              if (twoTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
                if (sign1 * sign2 < 0) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), VarManager::fgValues);
                } else if (fConfigSameSignHistograms.value) {
                  if (sign1 > 0) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), VarManager::fgValues);
                  } else {
                    fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), VarManager::fgValues);
                  }
                }
                for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                  if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                    if (sign1 * sign2 < 0) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    } else if (fConfigSameSignHistograms.value) {
                      if (sign1 > 0) {
                        fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                      } else {
                        fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                      }
                    }
                  }
                }
              }
            } // end loop (common cuts)
            int iPairCut = 0;
            for (auto cut = fPairCuts.begin(); cut != fPairCuts.end(); cut++, iPairCut++) {
              if (!((*cut)->IsSelected(VarManager::fgValues))) // apply pair cuts
                continue;
              pairFilter |= (static_cast<uint32_t>(1) << iPairCut);
              // Histograms with pair cuts
              if (sign1 * sign2 < 0) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
              } else if (fConfigSameSignHistograms.value) {
                if (sign1 > 0) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
                } else {
                  fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
                }
              }
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                  if (sign1 * sign2 < 0) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                  } else if (fConfigSameSignHistograms.value) {
                    if (sign1 > 0) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    } else {
                      fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    }
                  }
                }
              }
              // Histograms with pair cuts and common track cuts
              for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                if (twoTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
                  if (sign1 * sign2 < 0) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
                  } else if (fConfigSameSignHistograms.value) {
                    if (sign1 > 0) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
                    } else {
                      fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
                    }
                  }
                  for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                    if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                      if (sign1 * sign2 < 0) {
                        fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                      } else if (fConfigSameSignHistograms.value) {
                        if (sign1 > 0) {
                          fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                        } else {
                          fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                        }
                      }
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
      } // end inner assoc loop (leg A)
    } // end event loop
  }

  // Function to run same event triplets (e.g. D+->K-pi+pi+)
  void runThreeProng(MyEventsVtxCovSelected const& events, PresliceUnsorted<MyBarrelAssocs>& preslice, MyBarrelAssocs const& /*assocs*/, MyBarrelTracksWithCovWithAmbiguities const& tracks, ReducedA3MCEvents const& /*mcEvents*/, ReducedA3MCTracks const& /*mcTracks*/, VarManager::PairCandidateType tripletType)
  {
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEventAlice3<gkEventFillMapWithCov>(event, VarManager::fgValues);

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

      // Based on triplet type, make suitable combinations of the partitions
      if (tripletType == VarManager::kTripleCandidateToPKPi) {
        for (const auto& [a1, a2, a3] : combinations(soa::CombinationsFullIndexPolicy(groupedLegAAssocs, groupedLegBAssocs, groupedLegCAssocs))) {
          readTriplet(a1, a2, a3, tracks, event, tripletType);
        }
      } else if (tripletType == VarManager::kTripleCandidateToKPiPi) {
        for (const auto& a1 : groupedLegAAssocs) {
          for (const auto& [a2, a3] : combinations(groupedLegBAssocs, groupedLegCAssocs)) {
            readTriplet(a1, a2, a3, tracks, event, tripletType);
          }
        }
      } else {
        LOG(fatal) << "Given tripletType not recognized. Don't know how to make combinations!\n";
      }
    } // end event loop
  }

  // Helper function to process triplet
  void readTriplet(MyBarrelAssocs::iterator const& a1, MyBarrelAssocs::iterator const& a2, MyBarrelAssocs::iterator const& a3, MyBarrelTracksWithCovWithAmbiguities const& /*tracks*/, MyEventsVtxCovSelected::iterator const& /*event*/, VarManager::PairCandidateType tripletType)
  {
    uint32_t mcDecision = 0;

    uint32_t threeTrackFilter = 0;
    uint32_t threeTrackCommonFilter = 0;
    for (int icut = 0; icut < fNLegCuts; ++icut) {
      // Find out which leg cut combinations the triplet passes
      if ((a1.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut]) && (a2.isBarrelSelected_raw() & fConstructedLegBFilterMasksMap[icut]) && (a3.isBarrelSelected_raw() & fConstructedLegCFilterMasksMap[icut])) {
        threeTrackFilter |= (static_cast<uint32_t>(1) << icut);
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

    auto t1 = a1.template reducedA3track_as<MyBarrelTracksWithCovWithAmbiguities>();
    auto t2 = a2.template reducedA3track_as<MyBarrelTracksWithCovWithAmbiguities>();
    auto t3 = a3.template reducedA3track_as<MyBarrelTracksWithCovWithAmbiguities>();

    // Avoid self-pairs
    if (t1 == t2 || t1 == t3 || t2 == t3) {
      return;
    }
    // Check charge
    if (tripletType == VarManager::kTripleCandidateToKPiPi) {
      if (!((t1.sign() == -1 && t2.sign() == 1 && t3.sign() == 1) || (t1.sign() == 1 && t2.sign() == -1 && t3.sign() == -1))) {
        return;
      }
    }
    if (tripletType == VarManager::kTripleCandidateToPKPi) {
      if (!((t1.sign() == 1 && t2.sign() == -1 && t3.sign() == 1) || (t1.sign() == -1 && t2.sign() == 1 && t3.sign() == -1))) {
        return;
      }
    }

    // store the ambiguity of the three legs in the last 3 digits of the two-track filter
    if (t1.barrelAmbiguityInBunch() > 1 || t1.barrelAmbiguityOutOfBunch() > 1) {
      threeTrackFilter |= (static_cast<uint32_t>(1) << 29);
    }
    if (t2.barrelAmbiguityInBunch() > 1 || t2.barrelAmbiguityOutOfBunch() > 1) {
      threeTrackFilter |= (static_cast<uint32_t>(1) << 30);
    }
    if (t3.barrelAmbiguityInBunch() > 1 || t3.barrelAmbiguityOutOfBunch() > 1) {
      threeTrackFilter |= (static_cast<uint32_t>(1) << 31);
    }

    // run MC matching for this triplet
    int isig = 0;
    mcDecision = 0;
    for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
      if (t1.has_reducedA3MCTrack() && t2.has_reducedA3MCTrack() && t3.has_reducedA3MCTrack()) {
        if ((*sig)->CheckSignal(true, t1.reducedA3MCTrack(), t2.reducedA3MCTrack(), t3.reducedA3MCTrack())) {
          mcDecision |= (static_cast<uint32_t>(1) << isig);
        }
      }
    } // end loop over MC signals

    VarManager::FillTriple(t1, t2, t3, VarManager::fgValues, tripletType);
    /* TODO: Reimplement when secondary vertexing is available
    if constexpr (TThreeProngFitter) {
      VarManager::FillTripletVertexing<gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, t1, t2, t3, tripletType);
    }*/

    // Fill histograms
    bool isAmbi = false;
    for (int icut = 0; icut < fNLegCuts; icut++) {
      isAmbi = (threeTrackFilter & (static_cast<uint32_t>(1) << 29)) || (threeTrackFilter & (static_cast<uint32_t>(1) << 30)) || (threeTrackFilter & (static_cast<uint32_t>(1) << 31));
      if (threeTrackFilter & (static_cast<uint32_t>(1) << icut)) {
        fHistMan->FillHistClass(Form("TripletsBarrelSE_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
        for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
          if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
            fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues); // matched signal
          }
        } // end loop (MC signals)
        if (fConfigQA && isAmbi) {
          fHistMan->FillHistClass(Form("TripletsBarrelSE_ambiguous_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
        }
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
          if (threeTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
            fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), VarManager::fgValues);
            for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
              if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues); // matched signal
              }
            } // end loop (MC signals)
          }
        } // end loop (common cuts)
        int iPairCut = 0;
        for (auto cut = fPairCuts.begin(); cut != fPairCuts.end(); cut++, iPairCut++) {
          if (!((*cut)->IsSelected(VarManager::fgValues))) // apply pair cuts
            continue;
          // Histograms with pair cuts
          fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
          for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
            if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
              fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues); // matched signal
            }
          } // end loop (MC signals)
          // Histograms with pair cuts and common track cuts
          for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
            if (threeTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
              fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                  fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues); // matched signal
                }
              } // end loop (MC signals)
            }
          }
        } // end loop (pair cuts)
      }
    } // end loop (cuts)
  }

  void processKaonPionSkimmed(MyEventsVtxCovSelected const& events,
                              MyBarrelAssocs const& barrelAssocs,
                              MyBarrelTracksWithCovWithAmbiguities const& barrelTracks,
                              ReducedA3MCEvents const& mcEvents, ReducedA3MCTracks const& mcTracks)
  {
    runAsymmetricPairing(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks);
  }

  void processKaonPionPionSkimmed(MyEventsVtxCovSelected const& events,
                                  MyBarrelAssocs const& barrelAssocs,
                                  MyBarrelTracksWithCovWithAmbiguities const& barrelTracks,
                                  ReducedA3MCEvents const& mcEvents, ReducedA3MCTracks const& mcTracks)
  {
    runThreeProng(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks, VarManager::kTripleCandidateToKPiPi);
  }

  void processMCGen(ReducedA3MCTracks const& mcTracks)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // group all the MC tracks which belong to the MC event corresponding to the current reconstructed event
    // auto groupedMCTracks = tracksMC.sliceBy(aod::reducedA3trackMC::reducedA3MCEventId, event.reducedMCevent().globalIndex());
    for (const auto& mctrack : mcTracks) {

      VarManager::FillTrackMC(mcTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (const auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), VarManager::fgValues);
        }
      }
    }
  }

  PresliceUnsorted<ReducedA3MCTracks> perReducedMcEvent = aod::reducedA3trackMC::reducedA3MCEventId;

  void processMCGenWithEventSelection(soa::Filtered<MyEventsVtxCovSelected> const& events,
                                      ReducedA3MCEvents const& /*mcEvents*/, ReducedA3MCTracks const& mcTracks)
  {
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedA3MCEvent()) {
        continue;
      }

      auto groupedMCTracks = mcTracks.sliceBy(perReducedMcEvent, event.reducedA3MCEventId());
      groupedMCTracks.bindInternalIndicesTo(&mcTracks);
      for (const auto& track : groupedMCTracks) {

        VarManager::FillTrackMC(mcTracks, track);

        auto track_raw = groupedMCTracks.rawIteratorAt(track.globalIndex());
        for (const auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, track_raw)) {
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), VarManager::fgValues);
          }
        }
      }
    } // end loop over reconstructed events
  }

  void processDummy(MyEvents const&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisAsymmetricPairing, processKaonPionSkimmed, "Run kaon pion pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processKaonPionPionSkimmed, "Run kaon pion pion triplets, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processMCGen, "Loop over MC particle stack and fill generator level histograms", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processMCGenWithEventSelection, "Loop over MC particle stack and fill generator level histograms", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
    adaptAnalysisTask<AnalysisPrefilterSelection>(cfgc),
    adaptAnalysisTask<AnalysisSameEventPairing>(cfgc),
    adaptAnalysisTask<AnalysisAsymmetricPairing>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses, const char* histGroups)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //  The histogram classes are provided in the histClasses string, separated by semicolon ";"
  //  The histogram classes and their components histograms are defined below depending on the name of the histogram class
  //
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (int iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    TString histName = histGroups;
    // NOTE: The level of detail for histogramming can be controlled via configurables
    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", histName);
    }

    if (classStr.Contains("SameBunchCorrelations") || classStr.Contains("OutOfBunchCorrelations")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "two-collisions", histName);
    }

    // TODO: CHANGE TO PROPER PID

    if ((classStr.Contains("Track") || classStr.Contains("Assoc")) && !classStr.Contains("Pairs")) {
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
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "ambiguity");
        }
      }
    }

    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }

    if (classStr.Contains("Triplets")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }

    if (classStr.Contains("MCTruthGenPair")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_pair", histName);
    }

    if (classStr.Contains("MCTruthGenSelBR")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_triple");
    } else if (classStr.Contains("MCTruthGen")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_track");
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

  } // end loop over histogram classes
}
