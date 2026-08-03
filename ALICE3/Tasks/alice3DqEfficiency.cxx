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
/// \author Ionut Cristian Arsene <iarsene@cern.ch>, Oslo
/// \author Alexander Tiekoetter <alexander.tiekoetter@cern.ch>, Muenster
/// \brief Skimming Task for DQ Table Maker
/// \file alice3DqEfficiency.cxx

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

#include "ALICE3/DataModel/ReducedTablesAlice3.h"
#include "Common/Core/TableHelper.h"

#include <DetectorsBase/MatLayerCylSet.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <THashList.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
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
using namespace o2::common::core;

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
DECLARE_SOA_COLUMN(LxyeePoleMass, lxyeePoleMass, float);
DECLARE_SOA_COLUMN(Lzee, lzee, float);
DECLARE_SOA_COLUMN(MultiplicityFT0A, multiplicityFT0A, float);
DECLARE_SOA_COLUMN(MultiplicityFT0C, multiplicityFT0C, float);
DECLARE_SOA_COLUMN(PercentileFT0M, percentileFT0M, float);
DECLARE_SOA_COLUMN(MultiplicityNContrib, multiplicityNContrib, float);
DECLARE_SOA_COLUMN(AmbiguousInBunchPairs, ambiguousInBunchPairs, bool);
DECLARE_SOA_COLUMN(AmbiguousOutOfBunchPairs, ambiguousOutOfBunchPairs, bool);
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

using MyEvents = soa::Join<aod::ReA3Events, aod::ReducedA3MCEventLabels>;
using MyEventsSelected = soa::Join<aod::ReA3Events, aod::EventCuts, aod::ReducedA3MCEventLabels>;
using MyEventsVtxCov = soa::Join<aod::ReA3Events, aod::ReducedA3EventsVtxCov, aod::ReducedA3MCEventLabels>;
using MyEventsVtxCovSelected = soa::Join<aod::ReA3Events, aod::ReducedA3EventsVtxCov, aod::EventCuts, aod::ReducedA3MCEventLabels>;

using MyBarrelAssocs = soa::Join<aod::ReducedA3TracksAssoc, aod::BarrelTrackCuts>;
using MyBarrelAssocsPrefilter = soa::Join<aod::ReducedA3TracksAssoc, aod::BarrelTrackCuts, aod::Prefilter>;

using MyBarrelTracks = soa::Join<aod::ReA3Tracks, aod::ReducedA3TracksBarrel,
                                 aod::ReducedA3PIDTOF, aod::ReducedA3PIDRich, aod::ReducedA3PIDRichSignals,
                                 aod::ReducedA3TracksBarrelLabels>;

using MyBarrelTracksWithCov = soa::Join<aod::ReA3Tracks, aod::ReducedA3TracksBarrel, aod::ReducedA3TracksBarrelCov,
                                        aod::ReducedA3PIDTOF, aod::ReducedA3PIDRich, aod::ReducedA3PIDRichSignals,
                                        aod::ReducedA3TracksBarrelLabels>;

using MyBarrelTracksWithCovWithAmbiguities = soa::Join<aod::ReA3Tracks, aod::ReducedA3TracksBarrel, aod::ReducedA3TracksBarrelCov,
                                                       aod::ReducedA3PIDTOF, aod::ReducedA3PIDRich, aod::ReducedA3PIDRichSignals,
                                                       aod::BarrelAmbiguities, aod::ReducedA3TracksBarrelLabels>;

constexpr static uint32_t GkEventFillMap = VarManager::ObjTypes::ReducedEvent;
constexpr static uint32_t GkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventVtxCov;

constexpr static uint32_t GkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t GkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;

namespace dqefficiency_helpers
{
inline float* varValues() { return static_cast<float*>(VarManager::fgValues); }
inline TString* varNames() { return static_cast<TString*>(VarManager::fgVariableNames); }
inline TString* varUnits() { return static_cast<TString*>(VarManager::fgVariableUnits); }
} // namespace dqefficiency_helpers

// Global function used to define needed histogram classes
void defineHistograms(HistogramManager* histMan, const TString& histClasses, const char* histGroups); // defines histograms for all tasks

constexpr int TwoProng = 2;
constexpr int ThreeProng = 3;

// Analysis task that produces event decisions and the Hash table used in event mixing
struct Alice3DqEfficiencyAnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};
  Configurable<std::string> cfgMixingVars{"cfgMixingVars", "", "Mixing configs separated by a comma, default no mixing"};
  Configurable<std::string> cfgEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> cfgEventCutsJSON{"cfgEventCutsJSON", "", "Additional event cuts specified in JSON format"};
  Configurable<bool> cfgQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> cfgAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> cfgAddEventMCHistogram{"cfgAddEventMCHistogram", "generator", "Comma separated list of histograms"};
  Configurable<std::string> cfgAddJSONHistograms{"cfgAddJSONHistograms", "", "Add event histograms defined via JSON formatting (see HistogramsLibrary)"};

  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;

  AnalysisCompositeCut* fEventCut = nullptr;

  std::map<int64_t, bool> fSelMap; // key: reduced event global index, value: event selection decision

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }
    VarManager::SetDefaultVarNames();

    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = cfgEventCuts.value;
    if (eventCutStr != "") {
      AnalysisCut* cut = dqcuts::GetAnalysisCut(eventCutStr.Data());
      if (cut != nullptr) {
        fEventCut->AddCut(cut);
      }
    }
    // Additional cuts via JSON
    TString eventCutJSONStr = cfgEventCutsJSON.value;
    if (eventCutJSONStr != "") {
      std::vector<AnalysisCut*> jsonCuts = dqcuts::GetCutsFromJSON(eventCutJSONStr.Data());
      for (const auto& cutIt : jsonCuts) {
        fEventCut->AddCut(cutIt);
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (cfgQA) {
      fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(true);
      fHistMan->SetDefaultVarNames(dqefficiency_helpers::varNames(), dqefficiency_helpers::varUnits());
      defineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;", cfgAddEventHistogram.value.data());
      defineHistograms(fHistMan, "EventsMC", cfgAddEventMCHistogram.value.data());
      dqhistograms::AddHistogramsFromJSON(fHistMan, cfgAddJSONHistograms.value.c_str()); // aditional histograms via JSON
      VarManager::SetUseVars(fHistMan->GetUsedVars());
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    TString mixVarsString = cfgMixingVars.value;
    std::unique_ptr<TObjArray> objArray(mixVarsString.Tokenize(","));
    if (objArray->GetEntries() > 0) {
      fMixHandler = new MixingHandler("mixingHandler", "mixing handler");
      fMixHandler->Init();
      for (int iVar = 0; iVar < objArray->GetEntries(); ++iVar) {
        dqmixing::SetUpMixing(fMixHandler, objArray->At(iVar)->GetName());
      }
    }
  }

  void runEventSelection(MyEventsVtxCov const& events, ReA3MCEvents const& mcEvents)
  {
    fSelMap.clear();

    for (const auto& event : events) {
      // Reset the fValues array and fill event observables
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEventAlice3<GkEventFillMap>(event);
      if (event.has_reA3MCEvent()) {
        VarManager::FillEventAlice3<VarManager::ObjTypes::ReducedEventMC>(event.reA3MCEvent());
      }

      bool decision = false;
      // if QA is requested fill histograms before event selections
      if (cfgQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", dqefficiency_helpers::varValues()); // automatically fill all the histograms in the class Event
      }
      if (fEventCut->IsSelected(dqefficiency_helpers::varValues())) {
        if (cfgQA) {
          fHistMan->FillHistClass("Event_AfterCuts", dqefficiency_helpers::varValues());
        }
        decision = true;
      }
      fSelMap[event.globalIndex()] = decision;
      if (fMixHandler != nullptr) {
        int hh = fMixHandler->FindEventCategory(dqefficiency_helpers::varValues());
        hash(hh);
      }
    }

    for (const auto& event : mcEvents) {
      // Reset the fValues array and fill event observables
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEventAlice3<VarManager::ObjTypes::ReducedEventMC>(event);
      if (cfgQA) {
        fHistMan->FillHistClass("EventsMC", dqefficiency_helpers::varValues());
      }
    }
  }

  void publishSelections(MyEventsVtxCov const& events)
  {
    // publish the table
    auto evSel = static_cast<uint32_t>(0);
    for (const auto& event : events) {
      evSel = 0;
      if (fSelMap[event.globalIndex()]) { // event passed the user cuts
        evSel |= (static_cast<uint32_t>(1) << 0);
      }
      eventSel(evSel);
    }
  }

  void processSkimmed(MyEventsVtxCov const& events, aod::ReA3MCEvents const& mcEvents)
  {
    runEventSelection(events, mcEvents);
    publishSelections(events);
  }

  void processDummy(MyEventsVtxCov const&)
  {
    // do nothing
  }

  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisEventSelection, processDummy, "Dummy function", true);
};

// Produces a table with barrel track decisions (joinable to the ReducedA3TracksAssociations)
// Here one should add all the track cuts needed through the workflow (e.g. cuts for same-even pairing, electron prefiltering, track for dilepton-track correlations)
struct Alice3DqEfficiencyAnalysisTrackSelection {
  Produces<aod::BarrelTrackCuts> trackSel;
  Produces<aod::BarrelAmbiguities> trackAmbiguities;
  OutputObj<THashList> fOutputList{"output"};

  Configurable<std::string> cfgTrackCuts{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  Configurable<std::string> cfgBarrelTrackCutsJSON{"cfgBarrelTrackCutsJSON", "", "Additional list of barrel track cuts in JSON format"};
  Configurable<bool> cfgQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> cfgAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> cfgAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};
  Configurable<bool> cfgPublishAmbiguity{"cfgPublishAmbiguity", true, "If true, publish ambiguity table and fill QA histograms"};

  Configurable<std::string> cfgTrackMCSignals{"cfgTrackMCSignals", "", "Comma separated list of MC signals"};
  Configurable<std::string> cfgTrackMCsignalsJSON{"cfgTrackMCsignalsJSON", "", "Additional list of MC signals via JSON"};

  HistogramManager* fHistMan = nullptr;
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

    TString cutNamesStr = cfgTrackCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // add extra cuts from JSON
    TString addTrackCutsStr = cfgBarrelTrackCutsJSON.value;
    if (addTrackCutsStr != "") {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
      for (const auto& t : addTrackCuts) {
        fTrackCuts.push_back(static_cast<AnalysisCompositeCut*>(t));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    TString configSigNamesStr = cfgTrackMCSignals.value;
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
    TString addMCSignalsStr = cfgTrackMCsignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (const auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() != 1) { // NOTE: only 1 prong signals
          continue;
        }
        fMCSignals.push_back(mcIt);
      }
    }

    if (cfgQA) {
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(true);
      fHistMan->SetDefaultVarNames(dqefficiency_helpers::varNames(), dqefficiency_helpers::varUnits());

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

      defineHistograms(fHistMan, histClasses.Data(), cfgAddTrackHistogram.value.data());
      if (cfgPublishAmbiguity) {
        defineHistograms(fHistMan, "TrackBarrel_AmbiguityInBunch;TrackBarrel_AmbiguityOutOfBunch;", "ambiguity");
      }
      dqhistograms::AddHistogramsFromJSON(fHistMan, cfgAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                   // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  void runTrackSelection(ReducedA3TracksAssoc const& assocs, MyEventsVtxCovSelected const& /*events*/, MyBarrelTracksWithCov const& tracks, ReA3MCEvents const& /*eventsMC*/, ReA3MCTracks const& tracksMC)
  {
    fNAssocsInBunch.clear();
    fNAssocsOutOfBunch.clear();

    trackSel.reserve(assocs.size());
    trackAmbiguities.reserve(tracks.size());

    // Loop over associations
    for (const auto& assoc : assocs) {
      auto event = assoc.template reA3event_as<MyEventsVtxCovSelected>();
      if (!event.isEventSelected_bit(0)) {
        trackSel(0);
        continue;
      }

      VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
      // fill event information which might be needed in histograms/cuts that combine track and event properties
      VarManager::FillEventAlice3<GkEventFillMapWithCov>(event);
      if (event.has_reA3MCEvent()) {
        VarManager::FillEventAlice3<VarManager::ObjTypes::ReducedEventMC>(event.reA3MCEvent());
      }

      auto track = assoc.template reA3track_as<MyBarrelTracksWithCov>();
      VarManager::FillTrackAlice3<GkTrackFillMapWithCov>(track);
      // compute quantities which depend on the associated collision, such as DCA
      VarManager::FillTrackCollision<GkTrackFillMapWithCov>(track, event);

      bool isCorrectAssoc = false;
      if (track.has_reA3MCTrack()) {
        auto trackMC = track.reA3MCTrack();
        auto eventMCfromTrack = trackMC.reA3MCEvent();
        if (event.has_reA3MCEvent()) {
          isCorrectAssoc = (eventMCfromTrack.globalIndex() == event.reA3MCEvent().globalIndex());
        }
        VarManager::FillTrackMC(tracksMC, trackMC);
        VarManager::FillResolutions(trackMC, track);
      }

      if (cfgQA) {
        fHistMan->FillHistClass("AssocsBarrel_BeforeCuts", dqefficiency_helpers::varValues());
      }

      int iCut = 0;
      auto filterMap = static_cast<uint32_t>(0);
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, iCut++) {
        if ((*cut)->IsSelected(dqefficiency_helpers::varValues())) {
          filterMap |= (static_cast<uint32_t>(1) << iCut);
          if (cfgQA) {
            fHistMan->FillHistClass(fHistNamesReco[iCut], dqefficiency_helpers::varValues());
          }
        }
      } // end loop over cuts
      trackSel(filterMap);

      // compute MC matching decisions and fill histograms for matched associations
      int isig = 0;
      if (filterMap > 0 && track.has_reA3MCTrack()) {
        // loop over all MC signals
        for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {
          // check if this MC signal is matched
          if ((*sig)->CheckSignal(true, track.reA3MCTrack())) {
            // mcDecision |= (static_cast<uint32_t>(1) << isig);
            //  loop over cuts and fill histograms for the cuts that are fulfilled
            for (unsigned int icut = 0; icut < fTrackCuts.size(); icut++) {
              if ((filterMap & (static_cast<uint32_t>(1) << icut)) != 0u) {
                if (isCorrectAssoc) {
                  fHistMan->FillHistClass(fHistNamesMCMatched[icut * 2 * fMCSignals.size() + 2 * isig].Data(), dqefficiency_helpers::varValues());
                } else {
                  fHistMan->FillHistClass(fHistNamesMCMatched[icut * 2 * fMCSignals.size() + 2 * isig + 1].Data(), dqefficiency_helpers::varValues());
                }
              }
            } // end loop over cuts
          }
        }
      } // end if (filterMap > 0)

      // count the number of associations per track
      if (cfgPublishAmbiguity && filterMap > 0) {
        if (event.isEventSelected_bit(1)) {
          // for this track, count the number of associated collisions with in-bunch pileup and out of bunch associations
          if (!fNAssocsInBunch.contains(track.globalIndex())) {
            std::vector<int64_t> evVector = {event.globalIndex()};
            fNAssocsInBunch[track.globalIndex()] = evVector;
          } else {
            auto& evVector = fNAssocsInBunch[track.globalIndex()];
            evVector.push_back(event.globalIndex());
          }
        } else {
          if (!fNAssocsOutOfBunch.contains(track.globalIndex())) {
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
    if (cfgPublishAmbiguity) {
      if (cfgQA) {
        for (const auto& [trackIdx, evIndices] : fNAssocsInBunch) {
          if (evIndices.size() == 1) {
            continue;
          }
          auto track = tracks.rawIteratorAt(trackIdx);
          VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
          VarManager::FillTrackAlice3<GkTrackFillMapWithCov>(track);
          VarManager::fgValues[VarManager::kBarrelNAssocsInBunch] = static_cast<float>(evIndices.size());
          fHistMan->FillHistClass("TrackBarrel_AmbiguityInBunch", dqefficiency_helpers::varValues());
        } // end loop over in-bunch ambiguous tracks

        for (const auto& [trackIdx, evIndices] : fNAssocsOutOfBunch) {
          if (evIndices.size() == 1) {
            continue;
          }
          auto track = tracks.rawIteratorAt(trackIdx);
          VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
          VarManager::FillTrackAlice3<GkTrackFillMapWithCov>(track);
          VarManager::fgValues[VarManager::kBarrelNAssocsOutOfBunch] = static_cast<float>(evIndices.size());
          fHistMan->FillHistClass("TrackBarrel_AmbiguityOutOfBunch", dqefficiency_helpers::varValues());
        } // end loop over out-of-bunch ambiguous tracks
      }

      // publish the ambiguity table
      for (const auto& track : tracks) {
        int8_t nInBunch = 0;
        if (!fNAssocsInBunch.contains(track.globalIndex())) {
          nInBunch = fNAssocsInBunch[track.globalIndex()].size();
        }
        int8_t nOutOfBunch = 0;
        if (!fNAssocsOutOfBunch.contains(track.globalIndex())) {
          nOutOfBunch = fNAssocsOutOfBunch[track.globalIndex()].size();
        }
        trackAmbiguities(nInBunch, nOutOfBunch);
      }
    }
  } // end runTrackSelection()

  void processSkimmedWithCov(ReducedA3TracksAssoc const& assocs, MyEventsVtxCovSelected const& events, MyBarrelTracksWithCov const& tracks, ReA3MCEvents const& eventsMC, ReA3MCTracks const& tracksMC)
  {
    runTrackSelection(assocs, events, tracks, eventsMC, tracksMC);
  }

  void processDummy(MyEvents const&)
  {
    // do nothing
  }

  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisTrackSelection, processSkimmedWithCov, "Run barrel track selection on DQ skimmed tracks w/ cov matrix associations", false);
  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisTrackSelection, processDummy, "Dummy function", true);
};

struct Alice3DqEfficiencyAnalysisPrefilterSelection {
  Produces<aod::Prefilter> prefilter; // joinable with ReducedA3TracksAssoc

  // Configurables
  Configurable<std::string> cfgPrefilterTrackCut{"cfgPrefilterTrackCut", "", "Prefilter track cut"};
  Configurable<std::string> cfgPrefilterPairCut{"cfgPrefilterPairCut", "", "Prefilter pair cut"};
  Configurable<std::string> cfgTrackCuts{"cfgTrackCuts", "", "Track cuts for which to run the prefilter"};
  // Track related options
  Configurable<bool> cfgPropTrack{"cfgPropTrack", false, "Propgate tracks to associated collision to recalculate DCA and momentum vector"};

  std::map<uint32_t, uint32_t> fPrefilterMap;
  AnalysisCompositeCut* fPairCut = nullptr;
  uint32_t fPrefilterMask = 0;
  int fPrefilterCutBit = -1;

  PresliceUnsorted<aod::ReducedA3TracksAssoc> trackAssocsPerCollision = aod::reducedA3track_association::reA3eventId;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    bool runPrefilter = true;
    // get the list of track cuts to be prefiltered
    TString trackCutsStr = cfgTrackCuts.value;
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
    TString prefilterTrackCutStr = cfgPrefilterTrackCut.value;
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
        if (tempStr.CompareTo(cfgPrefilterTrackCut.value) == 0) {
          fPrefilterCutBit = icut;
        }
      }
      // setup the prefilter pair cut
      fPairCut = new AnalysisCompositeCut(true);
      TString pairCutStr = cfgPrefilterPairCut.value;
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
      auto track1 = assoc1.template reA3track_as<MyBarrelTracks>();
      auto track2 = assoc2.template reA3track_as<MyBarrelTracks>();

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

      if ((track1Candidate == 0 || !track2Loose) && (track2Candidate == 0 || !track1Loose)) {
        continue;
      }

      // compute pair quantities
      VarManager::FillPairAlice3<VarManager::kDecayToEE, GkTrackFillMap>(track1, track2);
      if (cfgPropTrack) {
        VarManager::FillPairCollision<VarManager::kDecayToEE, GkTrackFillMap>(event, track1, track2);
      }
      // if the pair fullfils the criteria, add an entry into the prefilter map for the two tracks
      if (fPairCut->IsSelected(dqefficiency_helpers::varValues())) {
        if (!fPrefilterMap.contains(track1.globalIndex()) && track1Candidate > 0) {
          fPrefilterMap[track1.globalIndex()] = track1Candidate;
        }
        if (!fPrefilterMap.contains(track2.globalIndex()) && track2Candidate > 0) {
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
        auto track = assoc.template reA3track_as<MyBarrelTracks>();
        mymap = -1;
        if (!fPrefilterMap.contains(track.globalIndex())) {
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

  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisPrefilterSelection, processBarrelSkimmed, "Run Prefilter selection on reduced tracks", false);
  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisPrefilterSelection, processDummy, "Do nothing", true);
};

// Run the same-event pairing
// This task assumes that both legs of the resonance fulfill the same cuts (symmetric decay channel)
//  Runs combinatorics for barrel-barrel combinations
// The task implements also process functions for running event mixing
struct Alice3DqEfficiencyAnalysisSameEventPairing {

  Produces<aod::Dielectrons> dielectronList;
  Produces<aod::DielectronsExtra> dielectronsExtraList;
  Produces<aod::DielectronsAll> dielectronAllList;
  Produces<aod::OniaMCTruth> mcTruthTableEffi;

  o2::base::MatLayerCylSet* fLUT = nullptr;
  OutputObj<THashList> fOutputList{"output"};

  struct : ConfigurableGroup {
    Configurable<std::string> cfgTrackCuts{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
    Configurable<std::string> cfgPairCuts{"cfgPairCuts", "", "Comma separated list of pair cuts, !!! Use only if you know what you are doing, otherwise leave empty"};
    Configurable<std::string> cfgMCGenAccCut{"cfgMCGenAccCut", "", "cut for MC generated particles acceptance"};
    // TODO: Add pair cuts via JSON
  } cfgTrackCuts;

  Configurable<bool> cfgQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> cfgAddSEPHistogram{"cfgAddSEPHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> cfgAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};

  struct : ConfigurableGroup {
    Configurable<bool> cfgFlatTables{"cfgFlatTables", false, "Produce a single flat tables with all relevant information of the pairs and single tracks"};
    Configurable<bool> cfgPropToPCA{"cfgPropToPCA", false, "Propagate tracks to secondary vertex"};
    Configurable<std::string> syst{"syst", "pp", "Collision system, pp or PbPb"};
    Configurable<float> energy{"energy", 13600, "Center of mass energy in GeV"};
  } fConfigOptions;

  struct : ConfigurableGroup {
    Configurable<std::string> cfgBarrelMCGenSignals{"cfgBarrelMCGenSignals", "", "Comma separated list of MC signals (generated)"};
    Configurable<std::string> cfgMCGenSignalsJSON{"cfgMCGenSignalsJSON", "", "Additional list of MC signals (generated) via JSON"};
    Configurable<std::string> cfgBarrelMCRecSignals{"cfgBarrelMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};
    Configurable<std::string> cfgMCRecSignalsJSON{"cfgMCRecSignalsJSON", "", "Comma separated list of MC signals (reconstructed) via JSON"};
    Configurable<bool> cfgSkimSignalOnly{"cfgSkimSignalOnly", false, "Configurable to select only matched candidates"};
  } fConfigMC;

  // Track related options
  Configurable<bool> cfgPropTrack{"cfgPropTrack", true, "Propgate tracks to associated collision to recalculate DCA and momentum vector"};

  // Filter filterEventSelected = aod::dqanalysisflags::isEventSelected & uint32_t(1);
  Filter eventFilter = aod::dqanalysisflags::isEventSelected > static_cast<uint32_t>(0);

  HistogramManager* fHistMan = nullptr;

  // keep histogram class names in maps, so we don't have to buld their names in the pair loops
  std::map<int, std::vector<TString>> fTrackHistNames;
  std::map<int, std::vector<TString>> fBarrelHistNamesMCmatched;
  std::vector<MCSignal*> fRecMCSignals;
  std::vector<MCSignal*> fGenMCSignals;

  std::vector<AnalysisCompositeCut> fPairCuts;
  AnalysisCompositeCut fMCGenAccCut;
  bool fUseMCGenAccCut = false;

  uint32_t fTrackFilterMask = 0; // mask for the track cuts required in this task to be applied on the barrel cuts produced upstream
  int fNCutsBarrel = 0;
  int fNPairCuts = 0;
  bool fHasTwoProngGenMCsignals = false;

  bool fEnableBarrelHistos = false;

  PresliceUnsorted<MyBarrelAssocsPrefilter> trackAssocsPerCollision = aod::reducedA3track_association::reA3eventId;

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
    TString pairCutNamesStr = cfgTrackCuts.cfgPairCuts.value;
    if (!pairCutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(pairCutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // get the list of cuts for tracks, check that they were played by the barrel selection tasks
    //   and make a mask for active cuts (barrel selection tasks may run more cuts, needed for other analyses)
    TString trackCutsStr = cfgTrackCuts.cfgTrackCuts.value;
    TObjArray* objArrayTrackCuts = nullptr;
    if (!trackCutsStr.IsNull()) {
      objArrayTrackCuts = trackCutsStr.Tokenize(",");
    }

    // Setting the MC rec signal names
    TString sigNamesStr = fConfigMC.cfgBarrelMCRecSignals.value;
    std::unique_ptr<TObjArray> objRecSigArray(sigNamesStr.Tokenize(","));
    for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() != TwoProng) { // NOTE: 2-prong signals required
          continue;
        }
        fRecMCSignals.push_back(sig);
      }
    }

    // Add the MCSignals from the JSON config
    TString addMCSignalsStr = fConfigMC.cfgMCRecSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (const auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() != TwoProng) { // NOTE: only 2 prong signals
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
    TString mcGenAccCutStr = cfgTrackCuts.cfgMCGenAccCut.value;
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
            if (cfgQA) {
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
            TString pairCutHistNamesStr = cfgTrackCuts.cfgPairCuts.value;
            if (!pairCutHistNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(pairCutHistNamesStr.Tokenize(","));
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
                if (cfgQA) {
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
    TString sigGenNamesStr = fConfigMC.cfgBarrelMCGenSignals.value;
    std::unique_ptr<TObjArray> objGenSigArray(sigGenNamesStr.Tokenize(","));
    for (int isig = 0; isig < objGenSigArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objGenSigArray->At(isig)->GetName());
      if (sig) {
        fGenMCSignals.push_back(sig);
      }
    }

    // Add the MCSignals from the JSON config
    TString addMCSignalsGenStr = fConfigMC.cfgMCGenSignalsJSON.value;
    if (addMCSignalsGenStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsGenStr.Data());
      for (const auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() > TwoProng) { // NOTE: only 2 prong signals
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
        } else if (sig->GetNProngs() == TwoProng) {
          histNames += Form("MCTruthGenPair_%s;", sig->GetName());
          histNames += Form("MCTruthGenPairSel_%s;", sig->GetName());
          fHasTwoProngGenMCsignals = true;
        }
      }
    }

    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(true);
    fHistMan->SetDefaultVarNames(dqefficiency_helpers::varNames(), dqefficiency_helpers::varUnits());

    VarManager::SetCollisionSystem((TString)fConfigOptions.syst, fConfigOptions.energy); // set collision system and center of mass energy

    defineHistograms(fHistMan, histNames.Data(), cfgAddSEPHistogram.value.data());     // define all histograms
    dqhistograms::AddHistogramsFromJSON(fHistMan, cfgAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
    VarManager::SetUseVars(fHistMan->GetUsedVars());                                   // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  // Template function to run same event pairing (barrel-barrel)
  void runSameEventPairing(MyEventsVtxCovSelected const& events, PresliceUnsorted<MyBarrelAssocsPrefilter>& preslice, MyBarrelAssocsPrefilter const& assocs, MyBarrelTracksWithCovWithAmbiguities const& /*tracks*/, ReA3MCEvents const& /*mcEvents*/, ReA3MCTracks const& /*mcTracks*/)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }

    TString cutNames = cfgTrackCuts.cfgTrackCuts.value;
    std::map<int, std::vector<TString>> histNames = fTrackHistNames;
    std::map<int, std::vector<TString>> histNamesMC = fBarrelHistNamesMCmatched;
    int ncuts = fNCutsBarrel;

    auto twoTrackFilter = static_cast<uint32_t>(0);
    int sign1 = 0;
    int sign2 = 0;
    auto mcDecision = static_cast<uint32_t>(0);
    bool isCorrectAssocLeg1 = false;
    bool isCorrectAssocLeg2 = false;

    int64_t reserveSize = 0;
    for (auto const& event : events) {
      if (event.isEventSelected_bit(0)) {
        auto groupedAssocs = assocs.sliceBy(preslice, event.globalIndex());
        size_t nGood = 0;
        for (auto const& t : groupedAssocs) {
          if (t.isBarrelSelected_raw() > 0u) {
            nGood++;
          }
        }
        reserveSize += nGood * (nGood - 1) / 2;
      }
    }

    dielectronList.reserve(reserveSize);
    dielectronsExtraList.reserve(reserveSize);

    if (fConfigOptions.cfgFlatTables.value) {
      dielectronAllList.reserve(reserveSize);
    }

    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEventAlice3<GkEventFillMap>(event, dqefficiency_helpers::varValues());
      VarManager::FillEventAlice3<VarManager::ObjTypes::ReducedEventMC>(event.reA3MCEvent(), dqefficiency_helpers::varValues());

      auto groupedAssocs = assocs.sliceBy(preslice, event.globalIndex());
      if (groupedAssocs.size() == 0) {
        continue;
      }

      for (const auto& [a1, a2] : o2::soa::combinations(groupedAssocs, groupedAssocs)) {

        twoTrackFilter = a1.isBarrelSelected_raw() & a2.isBarrelSelected_raw() & a1.isBarrelSelectedPrefilter_raw() & a2.isBarrelSelectedPrefilter_raw() & fTrackFilterMask;

        if (twoTrackFilter == 0u) { // the tracks must have at least one filter bit in common to continue
          continue;
        }

        auto t1 = a1.template reA3track_as<MyBarrelTracksWithCovWithAmbiguities>();
        auto t2 = a2.template reA3track_as<MyBarrelTracksWithCovWithAmbiguities>();
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
        int iSigMc = 0;
        mcDecision = 0;
        for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, iSigMc++) {
          if (t1.has_reA3MCTrack() && t2.has_reA3MCTrack()) {
            if ((*sig)->CheckSignal(true, t1.reA3MCTrack(), t2.reA3MCTrack())) {
              mcDecision |= (static_cast<uint32_t>(1) << iSigMc);
            }
          }
        } // end loop over MC signals
        if (t1.has_reA3MCTrack() && t2.has_reA3MCTrack()) {
          isCorrectAssocLeg1 = (t1.reA3MCTrack().reA3MCEvent() == event.reA3MCEvent());
          isCorrectAssocLeg2 = (t2.reA3MCTrack().reA3MCEvent() == event.reA3MCEvent());
        }

        VarManager::FillPairAlice3<VarManager::kDecayToEE, GkTrackFillMap>(t1, t2);
        if (cfgPropTrack) {
          VarManager::FillPairCollision<VarManager::kDecayToEE, GkTrackFillMap>(event, t1, t2);
        }

        VarManager::FillPairVertexingAlice3<VarManager::kDecayToEE, GkEventFillMap, GkTrackFillMap>(event, t1, t2, fConfigOptions.cfgPropToPCA);
        if (!fConfigMC.cfgSkimSignalOnly || mcDecision > 0) {
          dielectronList(event.globalIndex(), VarManager::fgValues[VarManager::kMass],
                         VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi],
                         t1.sign() + t2.sign(), twoTrackFilter, mcDecision);
        }

        // Fill histograms
        bool isAmbiInBunch = false;
        bool isAmbiOutOfBunch = false;

        for (int icut = 0; icut < ncuts; icut++) {
          if ((twoTrackFilter & (static_cast<uint32_t>(1) << icut)) != 0u) {
            isAmbiInBunch = ((twoTrackFilter & (static_cast<uint32_t>(1) << 28)) != 0u) || ((twoTrackFilter & (static_cast<uint32_t>(1) << 29)) != 0u);
            isAmbiOutOfBunch = ((twoTrackFilter & (static_cast<uint32_t>(1) << 30)) != 0u) || ((twoTrackFilter & (static_cast<uint32_t>(1) << 31)) != 0u);
            if (sign1 * sign2 < 0) {                                                                 // +- pairs
              fHistMan->FillHistClass(histNames[icut][0].Data(), dqefficiency_helpers::varValues()); // reconstructed, unmatched
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {                     // loop over MC signals
                if ((mcDecision & (static_cast<uint32_t>(1) << isig)) != 0u) {
                  fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][0].Data(), dqefficiency_helpers::varValues()); // matched signal
                  if (cfgQA) {
                    if (isCorrectAssocLeg1 && isCorrectAssocLeg2) { // correct track-collision association
                      fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][3].Data(), dqefficiency_helpers::varValues());
                    } else { // incorrect track-collision association
                      fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][4].Data(), dqefficiency_helpers::varValues());
                    }
                    if (isAmbiInBunch) { // ambiguous in bunch
                      fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][5].Data(), dqefficiency_helpers::varValues());
                      if (isCorrectAssocLeg1 && isCorrectAssocLeg2) {
                        fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][6].Data(), dqefficiency_helpers::varValues());
                      } else {
                        fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][7].Data(), dqefficiency_helpers::varValues());
                      }
                    }
                    if (isAmbiOutOfBunch) { // ambiguous out of bunch
                      fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][8].Data(), dqefficiency_helpers::varValues());
                      if (isCorrectAssocLeg1 && isCorrectAssocLeg2) {
                        fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][9].Data(), dqefficiency_helpers::varValues());
                      } else {
                        fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][10].Data(), dqefficiency_helpers::varValues());
                      }
                    }
                  }
                }
                if (cfgQA) {
                  if (isAmbiInBunch) {
                    fHistMan->FillHistClass(histNames[icut][3].Data(), dqefficiency_helpers::varValues());
                  }
                  if (isAmbiOutOfBunch) {
                    fHistMan->FillHistClass(histNames[icut][3 + 3].Data(), dqefficiency_helpers::varValues());
                  }
                }
              }
            } else {
              if (sign1 > 0) { // ++ pairs
                fHistMan->FillHistClass(histNames[icut][1].Data(), dqefficiency_helpers::varValues());
                for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                  if ((mcDecision & (static_cast<uint32_t>(1) << isig)) != 0u) {
                    fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][1].Data(), dqefficiency_helpers::varValues());
                  }
                }
                if (cfgQA) {
                  if (isAmbiInBunch) {
                    fHistMan->FillHistClass(histNames[icut][4].Data(), dqefficiency_helpers::varValues());
                  }
                  if (isAmbiOutOfBunch) {
                    fHistMan->FillHistClass(histNames[icut][4 + 3].Data(), dqefficiency_helpers::varValues());
                  }
                }
              } else { // -- pairs
                fHistMan->FillHistClass(histNames[icut][2].Data(), dqefficiency_helpers::varValues());
                for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                  if ((mcDecision & (static_cast<uint32_t>(1) << isig)) != 0u) {
                    fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][2].Data(), dqefficiency_helpers::varValues());
                  }
                }
                if (cfgQA) {
                  if (isAmbiInBunch) {
                    fHistMan->FillHistClass(histNames[icut][5].Data(), dqefficiency_helpers::varValues());
                  }
                  if (isAmbiOutOfBunch) {
                    fHistMan->FillHistClass(histNames[icut][5 + 3].Data(), dqefficiency_helpers::varValues());
                  }
                }
              }
            }
            for (unsigned int iPairCut = 0; iPairCut < fPairCuts.size(); iPairCut++) {
              AnalysisCompositeCut cut = fPairCuts.at(iPairCut);
              if (!(cut.IsSelected(dqefficiency_helpers::varValues()))) { // apply pair cuts
                continue;
              }
              if (sign1 * sign2 < 0) {
                fHistMan->FillHistClass(histNames[ncuts + icut * ncuts + iPairCut][0].Data(), dqefficiency_helpers::varValues());
              } else {
                if (sign1 > 0) {
                  fHistMan->FillHistClass(histNames[ncuts + icut * ncuts + iPairCut][1].Data(), dqefficiency_helpers::varValues());
                } else {
                  fHistMan->FillHistClass(histNames[ncuts + icut * ncuts + iPairCut][2].Data(), dqefficiency_helpers::varValues());
                }
              }
            } // end loop (pair cuts)
          }
        } // end loop (cuts)
      } // end loop over pairs of track associations
    } // end loop over events
  }

  PresliceUnsorted<ReA3MCTracks> perReducedMcEvent = aod::reducedA3trackMC::reA3MCEventId;

  void runMCGenWithGrouping(MyEventsVtxCovSelected const& events, ReA3MCEvents const& /*mcEvents*/, ReA3MCTracks const& mcTracks)
  {
    for (const auto& mctrack : mcTracks) {
      VarManager::FillTrackMC(mcTracks, mctrack);
      // if we have a mc generated acceptance cut, apply it here
      if (fUseMCGenAccCut) {
        if (!fMCGenAccCut.IsSelected(dqefficiency_helpers::varValues())) {
          continue;
        }
      }
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (const auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), dqefficiency_helpers::varValues());
        }
      }
    }
    // Fill Generated histograms taking into account selected collisions
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reA3MCEvent()) {
        continue;
      }

      for (const auto& track : mcTracks) {
        if (track.reA3MCEventId() != event.reA3MCEventId()) {
          continue;
        }
        VarManager::FillTrackMC(mcTracks, track);
        // if we have a mc generated acceptance cut, apply it here
        if (fUseMCGenAccCut) {
          if (!fMCGenAccCut.IsSelected(dqefficiency_helpers::varValues())) {
            continue;
          }
        }
        auto trackRaw = mcTracks.rawIteratorAt(track.globalIndex());
        for (const auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, trackRaw)) {
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), dqefficiency_helpers::varValues());
            mcTruthTableEffi(VarManager::fgValues[VarManager::kMCPt], VarManager::fgValues[VarManager::kMCEta], VarManager::fgValues[VarManager::kMCY], VarManager::fgValues[VarManager::kMCPhi], VarManager::fgValues[VarManager::kMCVz], VarManager::fgValues[VarManager::kMCVtxZ], VarManager::fgValues[VarManager::kMultFT0A], VarManager::fgValues[VarManager::kMultFT0C], VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kVtxNcontribReal]);
          }
        }
      }
    } // end loop over reconstructed events

    if (fHasTwoProngGenMCsignals) {
      for (const auto& [t1, t2] : combinations(mcTracks, mcTracks)) {
        auto t1Raw = mcTracks.rawIteratorAt(t1.globalIndex());
        auto t2Raw = mcTracks.rawIteratorAt(t2.globalIndex());
        if (t1Raw.reA3MCEventId() == t2Raw.reA3MCEventId()) {
          for (const auto& sig : fGenMCSignals) {
            if (sig->GetNProngs() != TwoProng) { // NOTE: 2-prong signals required here
              continue;
            }
            if (sig->CheckSignal(true, t1Raw, t2Raw)) {
              VarManager::FillPairMC<VarManager::kDecayToEE>(t1, t2);
              if (fUseMCGenAccCut) {
                if (!fMCGenAccCut.IsSelected(dqefficiency_helpers::varValues())) {
                  continue;
                }
              }
              fHistMan->FillHistClass(Form("MCTruthGenPair_%s", sig->GetName()), dqefficiency_helpers::varValues());
            }
          }
        }
      }
    }
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reA3MCEvent()) {
        continue;
      }
      // CURRENTLY ONLY FOR 1-GENERATION 2-PRONG SIGNALS
      if (fHasTwoProngGenMCsignals) {
        auto groupedMCTracks = mcTracks.sliceBy(perReducedMcEvent, event.reA3MCEventId());
        groupedMCTracks.bindInternalIndicesTo(&mcTracks);
        for (const auto& [t1, t2] : combinations(groupedMCTracks, groupedMCTracks)) {
          auto t1Raw = mcTracks.rawIteratorAt(t1.globalIndex());
          auto t2Raw = mcTracks.rawIteratorAt(t2.globalIndex());
          if (t1Raw.reA3MCEventId() == t2Raw.reA3MCEventId()) {
            for (const auto& sig : fGenMCSignals) {
              if (sig->GetNProngs() != TwoProng) { // NOTE: 2-prong signals required here
                continue;
              }
              if (sig->CheckSignal(true, t1Raw, t2Raw)) {
                VarManager::FillPairMC<VarManager::kDecayToEE>(t1, t2);
                if (fUseMCGenAccCut) {
                  if (!fMCGenAccCut.IsSelected(dqefficiency_helpers::varValues())) {
                    continue;
                  }
                }
                fHistMan->FillHistClass(Form("MCTruthGenPairSel_%s", sig->GetName()), dqefficiency_helpers::varValues());
              }
            }
          }
        }
      } // end loop over reconstructed events
    }
  }

  void processBarrelOnlySkimmed(MyEventsVtxCovSelected const& events,
                                MyBarrelAssocsPrefilter const& barrelAssocs,
                                MyBarrelTracksWithCovWithAmbiguities const& barrelTracks, ReA3MCEvents const& mcEvents, ReA3MCTracks const& mcTracks)
  {
    runSameEventPairing(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks);
    runMCGenWithGrouping(events, mcEvents, mcTracks);
  }

  PresliceUnsorted<ReA3MCTracks> perReducedMcGenEvent = aod::reducedA3trackMC::reA3MCEventId;

  void processMCGen(soa::Filtered<MyEventsVtxCovSelected> const& events, ReA3MCEvents const& /*mcEvents*/, ReA3MCTracks const& mcTracks)
  {
    for (const auto& mctrack : mcTracks) {
      VarManager::FillTrackMC(mcTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (const auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), dqefficiency_helpers::varValues());
        }
      }
    }

    // Fill Generated histograms taking into account selected collisions
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reA3MCEvent()) {
        continue;
      }
      VarManager::FillEventAlice3<GkEventFillMap>(event, dqefficiency_helpers::varValues());
      VarManager::FillEventAlice3<VarManager::ObjTypes::ReducedEventMC>(event.reA3MCEvent(), dqefficiency_helpers::varValues());

      for (const auto& track : mcTracks) {
        if (track.reA3MCEventId() != event.reA3MCEventId()) {
          continue;
        }
        VarManager::FillTrackMC(mcTracks, track);
        auto trackRaw = mcTracks.rawIteratorAt(track.globalIndex());
        for (const auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, trackRaw)) {
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), dqefficiency_helpers::varValues());
            mcTruthTableEffi(VarManager::fgValues[VarManager::kMCPt], VarManager::fgValues[VarManager::kMCEta], VarManager::fgValues[VarManager::kMCY], VarManager::fgValues[VarManager::kMCPhi], VarManager::fgValues[VarManager::kMCVz], VarManager::fgValues[VarManager::kMCVtxZ], VarManager::fgValues[VarManager::kMultFT0A], VarManager::fgValues[VarManager::kMultFT0C], VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kVtxNcontribReal]);
          }
        }
      }
    } // end loop over reconstructed events
    if (fHasTwoProngGenMCsignals) {
      for (const auto& [t1, t2] : combinations(mcTracks, mcTracks)) {
        auto t1Raw = mcTracks.rawIteratorAt(t1.globalIndex());
        auto t2Raw = mcTracks.rawIteratorAt(t2.globalIndex());
        if (t1Raw.reA3MCEventId() == t2Raw.reA3MCEventId()) {
          for (const auto& sig : fGenMCSignals) {
            if (sig->GetNProngs() != TwoProng) { // NOTE: 2-prong signals required here
              continue;
            }
            if (sig->CheckSignal(true, t1Raw, t2Raw)) {
              fHistMan->FillHistClass(Form("MCTruthGenPair_%s", sig->GetName()), dqefficiency_helpers::varValues());
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
      if (!event.has_reA3MCEvent()) {
        continue;
      }

      if (fHasTwoProngGenMCsignals) {
        for (const auto& [t1, t2] : combinations(mcTracks, mcTracks)) {
          if (t1.reA3MCEventId() != event.reA3MCEventId()) {
            continue;
          }
          if (t2.reA3MCEventId() != event.reA3MCEventId()) {
            continue;
          }
          auto t1Raw = mcTracks.rawIteratorAt(t1.globalIndex());
          auto t2Raw = mcTracks.rawIteratorAt(t2.globalIndex());
          if (t1Raw.reA3MCEventId() == t2Raw.reA3MCEventId()) {
            for (const auto& sig : fGenMCSignals) {
              if (sig->GetNProngs() != TwoProng) { // NOTE: 2-prong signals required here
                continue;
              }
              if (sig->CheckSignal(true, t1Raw, t2Raw)) {
                fHistMan->FillHistClass(Form("MCTruthGenPairSel_%s", sig->GetName()), dqefficiency_helpers::varValues());
              }
            }
          }
        }
      }
    } // end loop over reconstructed events
  }

  void processMCGenWithGrouping(soa::Filtered<MyEventsVtxCovSelected> const& events, ReA3MCEvents const& /*mcEvents*/, ReA3MCTracks const& mcTracks)
  {
    for (const auto& mctrack : mcTracks) {
      VarManager::FillTrackMC(mcTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (const auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), dqefficiency_helpers::varValues());
        }
      }
    }
    // Fill Generated histograms taking into account selected collisions
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reA3MCEvent()) {
        continue;
      }

      for (const auto& track : mcTracks) {
        if (track.reA3MCEventId() != event.reA3MCEventId()) {
          continue;
        }
        VarManager::FillTrackMC(mcTracks, track);
        auto trackRaw = mcTracks.rawIteratorAt(track.globalIndex());
        for (const auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, trackRaw)) {
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), dqefficiency_helpers::varValues());
            mcTruthTableEffi(VarManager::fgValues[VarManager::kMCPt], VarManager::fgValues[VarManager::kMCEta], VarManager::fgValues[VarManager::kMCY], VarManager::fgValues[VarManager::kMCPhi], VarManager::fgValues[VarManager::kMCVz], VarManager::fgValues[VarManager::kMCVtxZ], VarManager::fgValues[VarManager::kMultFT0A], VarManager::fgValues[VarManager::kMultFT0C], VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kVtxNcontribReal]);
          }
        }
      }
    } // end loop over reconstructed events
    if (fHasTwoProngGenMCsignals) {
      for (const auto& [t1, t2] : combinations(mcTracks, mcTracks)) {
        auto t1Raw = mcTracks.rawIteratorAt(t1.globalIndex());
        auto t2Raw = mcTracks.rawIteratorAt(t2.globalIndex());
        if (t1Raw.reA3MCEventId() == t2Raw.reA3MCEventId()) {
          for (const auto& sig : fGenMCSignals) {
            if (sig->GetNProngs() != TwoProng) { // NOTE: 2-prong signals required here
              continue;
            }
            if (sig->CheckSignal(true, t1Raw, t2Raw)) {
              fHistMan->FillHistClass(Form("MCTruthGenPair_%s", sig->GetName()), dqefficiency_helpers::varValues());
            }
          }
        }
      }
    }
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reA3MCEvent()) {
        continue;
      }
      // CURRENTLY ONLY FOR 1-GENERATION 2-PRONG SIGNALS
      if (fHasTwoProngGenMCsignals) {
        auto groupedMCTracks = mcTracks.sliceBy(perReducedMcEvent, event.reA3MCEventId());
        groupedMCTracks.bindInternalIndicesTo(&mcTracks);
        for (const auto& [t1, t2] : combinations(groupedMCTracks, groupedMCTracks)) {
          auto t1Raw = groupedMCTracks.rawIteratorAt(t1.globalIndex());
          auto t2Raw = groupedMCTracks.rawIteratorAt(t2.globalIndex());
          if (t1Raw.reA3MCEventId() == t2Raw.reA3MCEventId()) {
            for (const auto& sig : fGenMCSignals) {
              if (sig->GetNProngs() != TwoProng) { // NOTE: 2-prong signals required here
                continue;
              }
              if (sig->CheckSignal(true, t1Raw, t2Raw)) {
                fHistMan->FillHistClass(Form("MCTruthGenPairSel_%s", sig->GetName()), dqefficiency_helpers::varValues());
              }
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

  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisSameEventPairing, processBarrelOnlySkimmed, "Run barrel only pairing, with skimmed tracks", false);
  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisSameEventPairing, processMCGen, "Loop over MC particle stack and fill generator level histograms", false);
  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisSameEventPairing, processMCGenWithGrouping, "Loop over MC particle stack (grouped MCTracks) and fill generator level histograms", false);
  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisSameEventPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", true);
};

struct Alice3DqEfficiencyAnalysisAsymmetricPairing {

  Produces<aod::Ditracks> ditrackList;
  Produces<aod::DitracksExtra> ditrackExtraList;

  // Output objects
  OutputObj<THashList> fOutputList{"output"};

  // Configurables
  Configurable<std::string> cfgLegCuts{"cfgLegCuts", "", "<leg-A-1>:<leg-B-1>[:<leg-C-1>],[<leg-A-2>:<leg-B-2>[:<leg-C-1>],...]"};
  Configurable<uint32_t> cfgLegAFilterMask{"cfgLegAFilterMask", 0, "Filter mask corresponding to cuts in event-selection"};
  Configurable<uint32_t> cfgLegBFilterMask{"cfgLegBFilterMask", 0, "Filter mask corresponding to cuts in event-selection"};
  Configurable<uint32_t> cfgLegCFilterMask{"cfgLegCFilterMask", 0, "Filter mask corresponding to cuts in event-selection"};
  Configurable<std::string> cfgCommonTrackCuts{"cfgCommonTrackCuts", "", "Comma separated list of cuts to be applied to all legs"};
  Configurable<std::string> cfgPairCuts{"cfgPairCuts", "", "Comma separated list of pair cuts"};
  Configurable<std::string> cfgPairCutsJSON{"cfgPairCutsJSON", "", "Additional list of pair cuts in JSON format"};
  Configurable<bool> cfgSkipAmbiguousIdCombinations{"cfgSkipAmbiguousIdCombinations", true, "Choose whether to skip pairs/triples which pass a stricter combination of cuts, e.g. KKPi triplets for D+ -> KPiPi"};

  Configurable<std::string> cfgAsymmetricPairingHistogramsSubgroups{"cfgAsymmetricPairingHistogramsSubgroups", "barrel,vertexing", "Comma separated list of asymmetric-pairing histogram subgroups"};
  Configurable<bool> cfgSameSignHistograms{"cfgSameSignHistograms", false, "Include same sign pair histograms for 2-prong decays"};
  Configurable<bool> cfgReflectedHistograms{"cfgReflectedHistograms", false, "Include separate histograms for pairs which are reflections of previously counted pairs"};
  Configurable<bool> cfgQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> cfgAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};

  Configurable<std::string> cfgBarrelMCGenSignals{"cfgBarrelMCGenSignals", "", "Comma separated list of MC signals (generated)"};
  Configurable<std::string> cfgBarrelMCRecSignals{"cfgBarrelMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};
  Configurable<std::string> cfgMCRecSignalsJSON{"cfgMCRecSignalsJSON", "", "Additional list of MC signals (reconstructed) via JSON"};
  Configurable<std::string> cfgMCGenSignalsJSON{"cfgMCGenSignalsJSON", "", "Comma separated list of MC signals (generated) via JSON"};

  HistogramManager* fHistMan = nullptr;

  std::vector<AnalysisCompositeCut*> fPairCuts;
  int fNPairHistPrefixes = 0;

  std::vector<MCSignal*> fRecMCSignals;
  std::vector<MCSignal*> fGenMCSignals;

  // Filter masks to find legs in BarrelTrackCuts table
  uint32_t fLegAFilterMask = 0;
  uint32_t fLegBFilterMask = 0;
  uint32_t fLegCFilterMask = 0;
  // Maps tracking which combination of leg cuts the track cuts participate in
  std::map<int, uint32_t> fConstructedLegAFilterMasksMap;
  std::map<int, uint32_t> fConstructedLegBFilterMasksMap;
  std::map<int, uint32_t> fConstructedLegCFilterMasksMap;
  // Filter map for common track cuts
  uint32_t fCommonTrackCutMask = 0;
  // Map tracking which common track cut the track cuts correspond to
  std::map<int, uint32_t> fCommonTrackCutFilterMasks;

  int fNLegCuts = 0;
  int fNPairCuts = 0;
  int fNCommonTrackCuts = 0;
  // vectors for cut names and signal names, for easy access when calling FillHistogramList()
  std::vector<TString> fLegCutNames;
  std::vector<TString> fPairCutNames;
  std::vector<TString> fCommonCutNames;
  std::vector<TString> fRecMCSignalNames;

  Filter eventFilter = aod::dqanalysisflags::isEventSelected > static_cast<uint32_t>(0);

  PresliceUnsorted<MyBarrelAssocs> trackAssocsPerCollision = aod::reducedA3track_association::reA3eventId;
  // PresliceUnsorted<aod::ReducedA3TracksAssoc> trackAssocsPerCollision = aod::reducedA3track_association::reA3eventId;

  // Partitions for triplets and asymmetric pairs
  Partition<MyBarrelAssocs> legACandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & cfgLegAFilterMask) > static_cast<uint32_t>(0);
  Partition<MyBarrelAssocs> legBCandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & cfgLegBFilterMask) > static_cast<uint32_t>(0);
  Partition<MyBarrelAssocs> legCCandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & cfgLegCFilterMask) > static_cast<uint32_t>(0);

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
    fHistMan->SetDefaultVarNames(dqefficiency_helpers::varNames(), dqefficiency_helpers::varUnits());

    // Get the leg cut filter masks
    fLegAFilterMask = cfgLegAFilterMask.value;
    fLegBFilterMask = cfgLegBFilterMask.value;
    fLegCFilterMask = cfgLegCFilterMask.value;

    // Get the pair cuts
    TString pairCutNamesStr = cfgPairCuts.value;
    if (!pairCutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(pairCutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // Extra pair cuts via JSON
    TString addPairCutsStr = cfgPairCutsJSON.value;
    if (addPairCutsStr != "") {
      std::vector<AnalysisCut*> addPairCuts = dqcuts::GetCutsFromJSON(addPairCutsStr.Data());
      for (const auto& t : addPairCuts) {
        fPairCuts.push_back(static_cast<AnalysisCompositeCut*>(t));
        pairCutNamesStr += Form(",%s", t->GetName());
      }
    }
    std::unique_ptr<TObjArray> objArrayPairCuts(pairCutNamesStr.Tokenize(","));
    fNPairCuts = objArrayPairCuts->GetEntries();
    for (int j = 0; j < fNPairCuts; j++) {
      fPairCutNames.push_back(objArrayPairCuts->At(j)->GetName());
    }

    // Setting the MC rec signal names
    TString sigNamesStr = cfgBarrelMCRecSignals.value;
    std::unique_ptr<TObjArray> objRecSigArray(sigNamesStr.Tokenize(","));
    for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
      if (sig) {
        fRecMCSignals.push_back(sig);
      }
    }
    // Add the reco MCSignals from the JSON config
    TString addMCSignalsStr = cfgMCRecSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (const auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() != TwoProng && mcIt->GetNProngs() != ThreeProng) {
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
    int commonCutIdx = -1;
    TString commonNamesStr = cfgCommonTrackCuts.value;
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
      LOGF(fatal, "cfgLegAFilterMask has highest bit at position %d, but track-selection only has %d cuts!", static_cast<int>(std::floor(std::log2(fLegAFilterMask))) + 1, objArray->GetEntries());
    }
    if (static_cast<int>(std::floor(std::log2(fLegBFilterMask))) + 1 > objArray->GetEntries()) {
      LOGF(fatal, "cfgLegBFilterMask has highest bit at position %d, but track-selection only has %d cuts!", static_cast<int>(std::floor(std::log2(fLegBFilterMask))) + 1, objArray->GetEntries());
    }
    if (static_cast<int>(std::floor(std::log2(fLegCFilterMask))) + 1 > objArray->GetEntries()) {
      LOGF(fatal, "cfgLegCFilterMask has highest bit at position %d, but track-selection only has %d cuts!", static_cast<int>(std::floor(std::log2(fLegCFilterMask))) + 1, objArray->GetEntries());
    }

    // Get the cuts defining the legs
    uint32_t fConstructedLegAFilterMask = 0;
    uint32_t fConstructedLegBFilterMask = 0;
    uint32_t fConstructedLegCFilterMask = 0;
    TString legCutsStr = cfgLegCuts.value;
    std::unique_ptr<TObjArray> objArrayLegs(legCutsStr.Tokenize(","));
    if (objArrayLegs->GetEntries() == 0 && !isMCGen) {
      LOG(fatal) << "No cuts defining legs. Check the config!";
    }
    fNLegCuts = objArrayLegs->GetEntries();
    std::vector<bool> isThreeProng;
    int legAIdx = -1;
    int legBIdx = -1;
    int legCIdx = -1;
    // Loop over leg defining cuts
    for (int icut = 0; icut < fNLegCuts; ++icut) {
      TString legsStr = objArrayLegs->At(icut)->GetName();
      std::unique_ptr<TObjArray> legs(legsStr.Tokenize(":"));
      if (legs->GetEntries() == ThreeProng) {
        isThreeProng.push_back(true);
      } else if (legs->GetEntries() == TwoProng) {
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
        defineHistograms(fHistMan, Form("TripletsBarrelSE_%s", legsStr.Data()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
        if (cfgQA) {
          defineHistograms(fHistMan, Form("TripletsBarrelSE_ambiguous_%s", legsStr.Data()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
        }

        std::unique_ptr<TObjArray> objArrayCommon(commonNamesStr.Tokenize(","));
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
          defineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
        }

        TString pairCutHistNamesStr = cfgPairCuts.value;
        if (!pairCutHistNamesStr.IsNull()) { // if pair cuts
          std::unique_ptr<TObjArray> objArrayPair(pairCutHistNamesStr.Tokenize(","));
          fNPairCuts = objArrayPair->GetEntries();
          for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
            defineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s", legsStr.Data(), objArrayPair->At(iPairCut)->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              defineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
            } // end loop (common cuts)
          } // end loop (pair cuts)
        } // end if (pair cuts)

        // TODO: assign hist directories for the MC matched triplets for each (leg cut combo,MCsignal) combination
        if (!sigNamesStr.IsNull()) {
          for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
            auto sig = fRecMCSignals.at(isig);
            defineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s", legsStr.Data(), sig->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());

            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              defineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), sig->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
            }

            if (!pairCutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(pairCutNamesStr.Tokenize(","));
              for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
                defineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s_%s", legsStr.Data(), objArrayPair->At(iPairCut)->GetName(), sig->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
                for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                  defineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName(), sig->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
                } // end loop (common cuts)
              } // end loop (pair cuts)
            } // end if (pair cuts)
          } // end loop over MC signals
        } // end if (MC signals)
      } else {
        std::vector<TString> pairHistPrefixes = {"PairsBarrelSEPM"};
        if (cfgSameSignHistograms.value) {
          pairHistPrefixes.push_back("PairsBarrelSEPP");
          pairHistPrefixes.push_back("PairsBarrelSEMM");
        }
        fNPairHistPrefixes = pairHistPrefixes.size();

        for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
          defineHistograms(fHistMan, Form("%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
        }
        if (cfgQA) {
          for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
            defineHistograms(fHistMan, Form("%s_ambiguous_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
          }
        }
        if (cfgReflectedHistograms.value) {
          for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
            defineHistograms(fHistMan, Form("%s_reflected_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
          }
        }

        std::unique_ptr<TObjArray> objArrayCommon(commonNamesStr.Tokenize(","));
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
          for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
            defineHistograms(fHistMan, Form("%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
          }
        }

        if (!pairCutNamesStr.IsNull()) { // if pair cuts
          std::unique_ptr<TObjArray> objArrayPair(pairCutNamesStr.Tokenize(","));
          fNPairCuts = objArrayPair->GetEntries();
          for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
            for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
              defineHistograms(fHistMan, Form("%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayPair->At(iPairCut)->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
            }
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                defineHistograms(fHistMan, Form("%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
              }
            } // end loop (common cuts)
          } // end loop (pair cuts)
        } // end if (pair cuts)

        // assign hist directories for the MC matched triplets for each (leg cut combo,MCsignal) combination
        if (!sigNamesStr.IsNull()) {
          for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
            auto sig = fRecMCSignals.at(isig);
            for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
              defineHistograms(fHistMan, Form("%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), sig->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
            }
            if (cfgReflectedHistograms.value) {
              for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                defineHistograms(fHistMan, Form("%s_reflected_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), sig->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
              }
            }

            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                defineHistograms(fHistMan, Form("%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), sig->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
              }
            }

            if (!pairCutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(pairCutNamesStr.Tokenize(","));
              for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
                for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                  defineHistograms(fHistMan, Form("%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayPair->At(iPairCut)->GetName(), sig->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
                }
                for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                  for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                    defineHistograms(fHistMan, Form("%s_%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName(), sig->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());
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
    TString sigGenNamesStr = cfgBarrelMCGenSignals.value;
    std::unique_ptr<TObjArray> objGenSigArray(sigGenNamesStr.Tokenize(","));
    for (int isig = 0; isig < objGenSigArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objGenSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() == 1) { // NOTE: 1-prong signals required
          fGenMCSignals.push_back(sig);
          defineHistograms(fHistMan, Form("MCTruthGen_%s;", sig->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());    // TODO: Add these names to a std::vector to avoid using Form in the process function
          defineHistograms(fHistMan, Form("MCTruthGenSel_%s;", sig->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data()); // TODO: Add these names to a std::vector to avoid using Form in the process function
        }
      }
    }

    // Add the gen MCSignals from the JSON config
    addMCSignalsStr = cfgMCGenSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (const auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() == 1) {
          fGenMCSignals.push_back(mcIt);
          defineHistograms(fHistMan, Form("MCTruthGen_%s;", mcIt->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data());    // TODO: Add these names to a std::vector to avoid using Form in the process function
          defineHistograms(fHistMan, Form("MCTruthGenSel_%s;", mcIt->GetName()), cfgAsymmetricPairingHistogramsSubgroups.value.data()); // TODO: Add these names to a std::vector to avoid using Form in the process function
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

    dqhistograms::AddHistogramsFromJSON(fHistMan, cfgAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
    VarManager::SetUseVars(fHistMan->GetUsedVars());                                   // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  // Function to run same event pairing with asymmetric pairs (e.g. kaon-pion)
  void runAsymmetricPairing(MyEventsVtxCovSelected const& events, PresliceUnsorted<MyBarrelAssocs>& preslice, MyBarrelAssocs const& /*assocs*/, MyBarrelTracksWithCovWithAmbiguities const& /*tracks*/, ReA3MCEvents const& /*mcEvents*/, ReA3MCTracks const& /*mcTracks*/)
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
      VarManager::FillEventAlice3<GkEventFillMapWithCov>(event, dqefficiency_helpers::varValues());

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
          if (((a1.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut]) != 0u) && ((a2.isBarrelSelected_raw() & fConstructedLegBFilterMasksMap[icut]) != 0u)) {
            twoTrackFilter |= static_cast<uint32_t>(1) << icut;
            // If the supposed pion passes a kaon cut, this is a K+K-. Skip it.
            if (cfgSkipAmbiguousIdCombinations.value) {
              if ((a2.isBarrelSelected_raw() & fLegAFilterMask) != 0u) {
                isPairIdWrong = true;
              }
            }
          }
        }

        if (twoTrackFilter == 0u || isPairIdWrong) {
          continue;
        }

        // Find common track cuts both candidates pass
        twoTrackCommonFilter |= a1.isBarrelSelected_raw() & a2.isBarrelSelected_raw() & fCommonTrackCutMask;

        auto t1 = a1.template reA3track_as<MyBarrelTracksWithCovWithAmbiguities>();
        auto t2 = a2.template reA3track_as<MyBarrelTracksWithCovWithAmbiguities>();

        // Avoid self-pairs
        if (t1.globalIndex() == t2.globalIndex()) {
          continue;
        }

        bool isReflected = false;
        std::pair<int32_t, int32_t> trackIds(t1.globalIndex(), t2.globalIndex());
        if (fPairCount.contains(trackIds)) {
          // Double counting is possible due to track-collision ambiguity. Skip pairs which were counted before
          fPairCount[trackIds] += 1;
          continue;
        }
        if (fPairCount.contains(std::pair(trackIds.second, trackIds.first))) {
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
        int iSigMc = 0;
        mcDecision = 0;
        for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, iSigMc++) {
          if (t1.has_reA3MCTrack() && t2.has_reA3MCTrack()) {
            VarManager::FillPairMC<VarManager::kDecayToKPi>(t1.reA3MCTrack(), t2.reA3MCTrack());
            if ((*sig)->CheckSignal(true, t1.reA3MCTrack(), t2.reA3MCTrack())) {
              mcDecision |= static_cast<uint32_t>(1) << iSigMc;
            }
          }
        } // end loop over MC signals

        VarManager::FillPairAlice3<VarManager::kDecayToKPi, GkTrackFillMapWithCov>(t1, t2);
        VarManager::FillPairVertexingAlice3<VarManager::kDecayToKPi, GkEventFillMapWithCov, GkTrackFillMapWithCov>(event, t1, t2, true);

        // Fill histograms
        bool isAmbi = false;
        for (int icut = 0; icut < fNLegCuts; icut++) {
          if ((twoTrackFilter & (static_cast<uint32_t>(1) << icut)) != 0u) {
            isAmbi = ((twoTrackFilter & (static_cast<uint32_t>(1) << 30)) != 0u) || ((twoTrackFilter & (static_cast<uint32_t>(1) << 31)) != 0u);
            if (sign1 * sign2 < 0) {                                                                                             // +- pairs
              fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s", fLegCutNames[icut].Data()), dqefficiency_helpers::varValues()); // reconstructed, unmatched
              if (isAmbi && cfgQA) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPM_ambiguous_%s", fLegCutNames[icut].Data()), dqefficiency_helpers::varValues());
              }
              if (isReflected && cfgReflectedHistograms.value) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPM_reflected_%s", fLegCutNames[icut].Data()), dqefficiency_helpers::varValues());
              }
            } else if (cfgSameSignHistograms.value) {
              if (sign1 > 0) { // ++ pairs
                fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s", fLegCutNames[icut].Data()), dqefficiency_helpers::varValues());
                if (isAmbi && cfgQA) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPP_ambiguous_%s", fLegCutNames[icut].Data()), dqefficiency_helpers::varValues());
                }
                if (isReflected && cfgReflectedHistograms.value) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPP_reflected_%s", fLegCutNames[icut].Data()), dqefficiency_helpers::varValues());
                }
              } else { // -- pairs
                fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s", fLegCutNames[icut].Data()), dqefficiency_helpers::varValues());
                if (isAmbi && cfgQA) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEMM_ambiguous_%s", fLegCutNames[icut].Data()), dqefficiency_helpers::varValues());
                }
                if (isReflected && cfgReflectedHistograms) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEMM_reflected_%s", fLegCutNames[icut].Data()), dqefficiency_helpers::varValues());
                }
              }
            }
            for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
              if ((mcDecision & (static_cast<uint32_t>(1) << isig)) != 0u) {
                if (sign1 * sign2 < 0) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                  if (isReflected && cfgReflectedHistograms.value) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPM_reflected_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                  }
                } else if (cfgSameSignHistograms.value) {
                  if (sign1 > 0) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                    if (isReflected && cfgReflectedHistograms.value) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEPP_reflected_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                    }
                  } else {
                    fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                    if (isReflected && cfgReflectedHistograms.value) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEMM_reflected_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                    }
                  }
                }
              }
            }
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
              if ((twoTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) != 0u) {
                if (sign1 * sign2 < 0) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), dqefficiency_helpers::varValues());
                } else if (cfgSameSignHistograms.value) {
                  if (sign1 > 0) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), dqefficiency_helpers::varValues());
                  } else {
                    fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), dqefficiency_helpers::varValues());
                  }
                }
                for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                  if ((mcDecision & (static_cast<uint32_t>(1) << isig)) != 0u) {
                    if (sign1 * sign2 < 0) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                    } else if (cfgSameSignHistograms.value) {
                      if (sign1 > 0) {
                        fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                      } else {
                        fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                      }
                    }
                  }
                }
              }
            } // end loop (common cuts)
            int iPairCut = 0;
            for (auto cut = fPairCuts.begin(); cut != fPairCuts.end(); cut++, iPairCut++) {
              if (!((*cut)->IsSelected(dqefficiency_helpers::varValues()))) { // apply pair cuts
                continue;
              }
              pairFilter |= (static_cast<uint32_t>(1) << iPairCut);
              // Histograms with pair cuts
              if (sign1 * sign2 < 0) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), dqefficiency_helpers::varValues());
              } else if (cfgSameSignHistograms.value) {
                if (sign1 > 0) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), dqefficiency_helpers::varValues());
                } else {
                  fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), dqefficiency_helpers::varValues());
                }
              }
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                if ((mcDecision & (static_cast<uint32_t>(1) << isig)) != 0u) {
                  if (sign1 * sign2 < 0) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                  } else if (cfgSameSignHistograms.value) {
                    if (sign1 > 0) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                    } else {
                      fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                    }
                  }
                }
              }
              // Histograms with pair cuts and common track cuts
              for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                if ((twoTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) != 0u) {
                  if (sign1 * sign2 < 0) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), dqefficiency_helpers::varValues());
                  } else if (cfgSameSignHistograms.value) {
                    if (sign1 > 0) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), dqefficiency_helpers::varValues());
                    } else {
                      fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), dqefficiency_helpers::varValues());
                    }
                  }
                  for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                    if ((mcDecision & (static_cast<uint32_t>(1) << isig)) != 0u) {
                      if (sign1 * sign2 < 0) {
                        fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                      } else if (cfgSameSignHistograms.value) {
                        if (sign1 > 0) {
                          fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
                        } else {
                          fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues());
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
  template <bool TThreeProngFitter, typename TEvents, typename TTrackAssocs, typename TTracks>
  void runThreeProng(TEvents const& events, PresliceUnsorted<TTrackAssocs>& preslice, TTrackAssocs const& /*assocs*/, TTracks const& tracks, ReA3MCEvents const& /*mcEvents*/, ReA3MCTracks const& /*mcTracks*/, VarManager::PairCandidateType tripletType)
  {
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEventAlice3<GkEventFillMapWithCov>(event, dqefficiency_helpers::varValues());

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
          readTriplet<TThreeProngFitter>(a1, a2, a3, tracks, event, tripletType);
        }
      } else if (tripletType == VarManager::kTripleCandidateToKPiPi) {
        for (const auto& a1 : groupedLegAAssocs) {
          for (const auto& [a2, a3] : combinations(groupedLegBAssocs, groupedLegCAssocs)) {
            readTriplet<TThreeProngFitter>(a1, a2, a3, tracks, event, tripletType);
          }
        }
      } else {
        LOG(fatal) << "Given tripletType not recognized. Don't know how to make combinations!\n";
      }
    } // end event loop
  }

  template <bool TThreeProngFitter, typename TTrackAssoc, typename TTracks, typename TEvent>
  void readTriplet(TTrackAssoc const& a1, TTrackAssoc const& a2, TTrackAssoc const& a3, TTracks const& /*tracks*/, TEvent const& event, VarManager::PairCandidateType tripletType)
  {
    uint32_t mcDecision = 0;

    uint32_t threeTrackFilter = 0;
    uint32_t threeTrackCommonFilter = 0;
    for (int icut = 0; icut < fNLegCuts; ++icut) {
      // Find out which leg cut combinations the triplet passes
      if ((a1.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut]) && (a2.isBarrelSelected_raw() & fConstructedLegBFilterMasksMap[icut]) && (a3.isBarrelSelected_raw() & fConstructedLegCFilterMasksMap[icut])) {
        threeTrackFilter |= (static_cast<uint32_t>(1) << icut);
        if (tripletType == VarManager::kTripleCandidateToPKPi && cfgSkipAmbiguousIdCombinations.value) {
          // Check if the supposed pion passes as a proton or kaon, if so, skip this triplet. It is pKp or pKK.
          if ((a3.isBarrelSelected_raw() & fLegAFilterMask) || (a3.isBarrelSelected_raw() & fLegBFilterMask)) {
            return;
          }
          // Check if the supposed kaon passes as a proton, if so, skip this triplet. It is ppPi.
          if (a2.isBarrelSelected_raw() & fLegAFilterMask) {
            return;
          }
        }
        if (tripletType == VarManager::kTripleCandidateToKPiPi && cfgSkipAmbiguousIdCombinations.value) {
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

    auto t1 = a1.template reA3track_as<TTracks>();
    auto t2 = a2.template reA3track_as<TTracks>();
    auto t3 = a3.template reA3track_as<TTracks>();

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
    int iSigMc = 0;
    mcDecision = 0;
    for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, iSigMc++) {
      if (t1.has_reA3MCTrack() && t2.has_reA3MCTrack() && t3.has_reA3MCTrack()) {
        if ((*sig)->CheckSignal(true, t1.reA3MCTrack(), t2.reA3MCTrack(), t3.reA3MCTrack())) {
          mcDecision |= (static_cast<uint32_t>(1) << iSigMc);
        }
      }
    } // end loop over MC signals

    VarManager::FillTriple(t1, t2, t3, dqefficiency_helpers::varValues(), tripletType);
    if constexpr (TThreeProngFitter) {
      VarManager::FillTripletVertexingALICE3<GkEventFillMapWithCov, GkTrackFillMapWithCov>(event, t1, t2, t3, tripletType);
    }

    // Fill histograms
    bool isAmbi = false;
    for (int icut = 0; icut < fNLegCuts; icut++) {
      isAmbi = (threeTrackFilter & (static_cast<uint32_t>(1) << 29)) || (threeTrackFilter & (static_cast<uint32_t>(1) << 30)) || (threeTrackFilter & (static_cast<uint32_t>(1) << 31));
      if (threeTrackFilter & (static_cast<uint32_t>(1) << icut)) {
        fHistMan->FillHistClass(Form("TripletsBarrelSE_%s", fLegCutNames[icut].Data()), dqefficiency_helpers::varValues());
        for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
          if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
            fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues()); // matched signal
          }
        } // end loop (MC signals)
        if (cfgQA && isAmbi) {
          fHistMan->FillHistClass(Form("TripletsBarrelSE_ambiguous_%s", fLegCutNames[icut].Data()), dqefficiency_helpers::varValues());
        }
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
          if (threeTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
            fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), dqefficiency_helpers::varValues());
            for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
              if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues()); // matched signal
              }
            } // end loop (MC signals)
          }
        } // end loop (common cuts)
        int iPairCut = 0;
        for (auto cut = fPairCuts.begin(); cut != fPairCuts.end(); cut++, iPairCut++) {
          if (!((*cut)->IsSelected(dqefficiency_helpers::varValues()))) { // apply pair cuts
            continue;
          }
          // Histograms with pair cuts
          fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), dqefficiency_helpers::varValues());
          for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
            if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
              fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues()); // matched signal
            }
          } // end loop (MC signals)
          // Histograms with pair cuts and common track cuts
          for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
            if (threeTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
              fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), dqefficiency_helpers::varValues());
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                  fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), dqefficiency_helpers::varValues()); // matched signal
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
                              ReA3MCEvents const& mcEvents, ReA3MCTracks const& mcTracks)
  {
    runAsymmetricPairing(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks);
  }

  void processKaonPionPionSkimmed(MyEventsVtxCovSelected const& events,
                                  MyBarrelAssocs const& barrelAssocs,
                                  MyBarrelTracksWithCovWithAmbiguities const& barrelTracks,
                                  ReA3MCEvents const& mcEvents, ReA3MCTracks const& mcTracks)
  {
    runThreeProng<true>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks, VarManager::kTripleCandidateToKPiPi);
  }

  void processProtonKaonPionSkimmed(MyEventsVtxCovSelected const& events,
                                    MyBarrelAssocs const& barrelAssocs,
                                    MyBarrelTracksWithCovWithAmbiguities const& barrelTracks,
                                    ReA3MCEvents const& mcEvents, ReA3MCTracks const& mcTracks)
  {
    runThreeProng<true>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks, VarManager::kTripleCandidateToPKPi);
  }

  void processMCGen(ReA3MCTracks const& mcTracks)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // group all the MC tracks which belong to the MC event corresponding to the current reconstructed event
    // auto groupedMCTracks = tracksMC.sliceBy(aod::reducedA3trackMC::reA3MCEventId, event.reducedMCevent().globalIndex());
    for (const auto& mctrack : mcTracks) {

      VarManager::FillTrackMC(mcTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (const auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), dqefficiency_helpers::varValues());
        }
      }
    }
  }

  PresliceUnsorted<ReA3MCTracks> perReducedMcEvent = aod::reducedA3trackMC::reA3MCEventId;

  void processMCGenWithEventSelection(soa::Filtered<MyEventsVtxCovSelected> const& events,
                                      ReA3MCEvents const& /*mcEvents*/, ReA3MCTracks const& mcTracks)
  {
    for (const auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reA3MCEvent()) {
        continue;
      }

      auto groupedMCTracks = mcTracks.sliceBy(perReducedMcEvent, event.reA3MCEventId());
      groupedMCTracks.bindInternalIndicesTo(&mcTracks);
      for (const auto& track : groupedMCTracks) {

        VarManager::FillTrackMC(mcTracks, track);

        auto trackRaw = groupedMCTracks.rawIteratorAt(track.globalIndex());
        for (const auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, trackRaw)) {
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), dqefficiency_helpers::varValues());
          }
        }
      }
    } // end loop over reconstructed events
  }

  void processDummy(MyEvents const&)
  {
    // do nothing
  }

  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisAsymmetricPairing, processKaonPionSkimmed, "Run kaon pion pairing, with skimmed tracks", false);
  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisAsymmetricPairing, processKaonPionPionSkimmed, "Run kaon pion pion triplets, with skimmed tracks", false);
  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisAsymmetricPairing, processProtonKaonPionSkimmed, "Run proton kaon pion triplets, with skimmed tracks", false);
  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisAsymmetricPairing, processMCGen, "Loop over MC particle stack and fill generator level histograms", false);
  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisAsymmetricPairing, processMCGenWithEventSelection, "Loop over MC particle stack and fill generator level histograms", false);
  PROCESS_SWITCH(Alice3DqEfficiencyAnalysisAsymmetricPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Alice3DqEfficiencyAnalysisEventSelection>(cfgc),
    adaptAnalysisTask<Alice3DqEfficiencyAnalysisTrackSelection>(cfgc),
    adaptAnalysisTask<Alice3DqEfficiencyAnalysisPrefilterSelection>(cfgc),
    adaptAnalysisTask<Alice3DqEfficiencyAnalysisSameEventPairing>(cfgc),
    adaptAnalysisTask<Alice3DqEfficiencyAnalysisAsymmetricPairing>(cfgc)};
}

void defineHistograms(HistogramManager* histMan, const TString& histClasses, const char* histGroups)
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
