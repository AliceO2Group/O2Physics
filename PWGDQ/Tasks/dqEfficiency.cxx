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
// Analysis task for processing O2::DQ MC skimmed AODs
//
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include <TMath.h>
#include <TH1F.h>
#include <THashList.h>
#include <TString.h>
#include <iostream>
#include <vector>

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
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelected, isBarrelSelected, int);
DECLARE_SOA_COLUMN(IsMuonSelected, isMuonSelected, int);
} // namespace reducedevent

DECLARE_SOA_TABLE(EventCuts, "AOD", "EVENTCUTS", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "BARRELTRACKCUTS", dqanalysisflags::IsBarrelSelected);
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);
} // namespace o2::aod

//using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsMC>;
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedMCEventLabels>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::ReducedMCEventLabels>;
// TODO: make secondary vertexing optional
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::ReducedMCEventLabels>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedMCEventLabels>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::ReducedTracksBarrelLabels>;
using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::ReducedTracksBarrelLabels>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::ReducedTracksBarrelLabels>;
using MyBarrelTracksSelectedWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::ReducedTracksBarrelLabels>;

using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsLabels>;
using MyMuonTracksSelected = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::MuonTrackCuts, aod::ReducedMuonsLabels>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::ReducedMuonsLabels>;
using MyMuonTracksSelectedWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonTrackCuts, aod::ReducedMuonsLabels>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkMCEventFillMap = VarManager::ObjTypes::ReducedEventMC;
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
//constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCov;
constexpr static uint32_t gkParticleMCFillMap = VarManager::ObjTypes::ParticleMC;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

struct AnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  OutputObj<THashList> fOutputList{"output"};
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};

  HistogramManager* fHistMan;
  AnalysisCompositeCut* fEventCut;

  void init(o2::framework::InitContext&)
  {
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    VarManager::SetDefaultVarNames();
    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;"); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                 // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, typename TEvent, typename TEventsMC>
  void runSelection(TEvent const& event, TEventsMC const& mcEvents)
  {
    // Reset the values array
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);

    VarManager::FillEvent<TEventFillMap>(event);
    if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
      VarManager::FillEvent<TEventMCFillMap>(event.reducedMCevent());
    }
    if constexpr ((TEventMCFillMap & VarManager::ObjTypes::CollisionMC) > 0) {
      VarManager::FillEvent<TEventMCFillMap>(event.mcCollision());
    }
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
  }

  void processSkimmed(MyEvents::iterator const& event, aod::ReducedMCEvents const& mcEvents)
  {
    runSelection<gkEventFillMap, gkMCEventFillMap>(event, mcEvents);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDummy, "Dummy process function", false);
};

struct AnalysisTrackSelection {
  Produces<aod::BarrelTrackCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  Configurable<std::string> fConfigCuts{"cfgTrackCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMCSignals{"cfgTrackMCSignals", "", "Comma separated list of MC signals"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fTrackCuts;
  std::vector<MCSignal> fMCSignals; // list of signals to be checked
  std::vector<TString> fHistNamesReco;
  std::vector<std::vector<TString>> fHistNamesMCMatched;

  void init(o2::framework::InitContext&)
  {
    // Setting the cut names
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    TString configSigNamesStr = fConfigMCSignals.value;
    std::unique_ptr<TObjArray> sigNamesArray(configSigNamesStr.Tokenize(","));

    // Setting the MC signal names
    for (int isig = 0; isig < sigNamesArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(sigNamesArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() != 1) { // NOTE: only 1 prong signals
          continue;
        }
        fMCSignals.push_back(*sig);
      }
    }

    // Configure histogram classes for each track cut;
    // Add histogram classes for each track cut and for each requested MC signal (reconstructed tracks with MC truth)
    TString histClasses = "TrackBarrel_BeforeCuts;";
    for (auto& cut : fTrackCuts) {
      TString nameStr = Form("TrackBarrel_%s", cut.GetName());
      fHistNamesReco.push_back(nameStr);
      histClasses += Form("%s;", nameStr.Data());
      std::vector<TString> mcnames;
      for (auto& sig : fMCSignals) {
        TString nameStr2 = Form("TrackBarrel_%s_%s", cut.GetName(), sig.GetName());
        printf("Adding my histogram class %s\n", nameStr2.Data());
        mcnames.push_back(nameStr2);
        histClasses += Form("%s;", nameStr2.Data());
      }
      fHistNamesMCMatched.push_back(mcnames);
    }

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, histClasses.Data());  // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackFillMap, uint32_t TTrackMCFillMap, typename TEvent, typename TTracks, typename TEventsMC, typename TTracksMC>
  void runSelection(TEvent const& event, TTracks const& tracks, TEventsMC const& eventsMC, TTracksMC const& tracksMC)
  {
    VarManager::ResetValues(0, VarManager::kNMCParticleVariables);
    // fill event information which might be needed in histograms that combine track and event properties
    VarManager::FillEvent<TEventFillMap>(event);
    if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
      VarManager::FillEvent<TEventMCFillMap>(event.reducedMCevent());
    }
    if constexpr ((TEventMCFillMap & VarManager::ObjTypes::CollisionMC) > 0) {
      VarManager::FillEvent<TEventMCFillMap>(event.mcCollision());
    }

    uint32_t filterMap = 0;
    trackSel.reserve(tracks.size());
    for (auto& track : tracks) {
      filterMap = 0;

      VarManager::FillTrack<TTrackFillMap>(track); // compute track quantities
      // compute MC matched quantities
      if constexpr (TTrackMCFillMap & VarManager::ObjTypes::ReducedTrack) {
        VarManager::FillTrack<gkParticleMCFillMap>(track.reducedMCTrack());
      }
      if constexpr (TTrackMCFillMap & VarManager::ObjTypes::Track) {
        VarManager::FillTrack<gkParticleMCFillMap>(track.mcParticle());
      }

      if (fConfigQA) {
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }

      // compute track selection and publish the bit map
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << i);
          if (fConfigQA) {
            fHistMan->FillHistClass(fHistNamesReco[i].Data(), VarManager::fgValues);
          }
        }
      }
      trackSel(static_cast<int>(filterMap));
      if (!filterMap) {
        continue;
      }

      if (!fConfigQA) {
        continue;
      }

      // compute MC matching decisions
      uint32_t mcDecision = 0;
      int isig = 0;
      for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
          if ((*sig).CheckSignal(false, tracksMC, track.reducedMCTrack())) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
          if ((*sig).CheckSignal(false, tracksMC, track.template mcParticle_as<aod::McParticles_001>())) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }
      }

      // fill histograms
      for (unsigned int i = 0; i < fMCSignals.size(); i++) {
        if (!(mcDecision & (uint32_t(1) << i))) {
          continue;
        }
        for (unsigned int j = 0; j < fTrackCuts.size(); j++) {
          if (filterMap & (uint8_t(1) << j)) {
            fHistMan->FillHistClass(fHistNamesMCMatched[j][i].Data(), VarManager::fgValues);
          }
        } // end loop over cuts
      }   // end loop over MC signals
    }     // end loop over tracks
  }

  void processSkimmed(MyEventsSelected::iterator const& event, MyBarrelTracks const& tracks, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runSelection<gkEventFillMap, gkMCEventFillMap, gkTrackFillMap, gkParticleMCFillMap>(event, tracks, eventsMC, tracksMC);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmed, "Run barrel track selection on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummy, "Dummy process function", false);
};

struct AnalysisMuonSelection {
  Produces<aod::MuonTrackCuts> muonSel;
  OutputObj<THashList> fOutputList{"output"};
  Configurable<std::string> fConfigCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};
  Configurable<std::string> fConfigMCSignals{"cfgMuonMCSignals", "", "Comma separated list of MC signals"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fTrackCuts;
  std::vector<MCSignal> fMCSignals; // list of signals to be checked
  std::vector<TString> fHistNamesReco;
  std::vector<std::vector<TString>> fHistNamesMCMatched;

  void init(o2::framework::InitContext&)
  {
    // Setting the cut names
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    TString configSigNamesStr = fConfigMCSignals.value;
    std::unique_ptr<TObjArray> sigNamesArray(configSigNamesStr.Tokenize(","));

    // Setting the MC signal names
    for (int isig = 0; isig < sigNamesArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(sigNamesArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() != 1) { // NOTE: only 1 prong signals
          continue;
        }
        fMCSignals.push_back(*sig);
      }
    }

    // Configure histogram classes for each track cut;
    // Add histogram classes for each track cut and for each requested MC signal (reconstructed tracks with MC truth)
    TString histClasses = "Muon_BeforeCuts;";
    for (auto& cut : fTrackCuts) {
      TString nameStr = Form("Muon_%s;", cut.GetName());
      fHistNamesReco.push_back(nameStr);
      histClasses += Form("%s;", nameStr.Data());
      std::vector<TString> mcnames;
      for (auto& sig : fMCSignals) {
        TString nameStr2 = Form("TrackBarrel_%s_%s", cut.GetName(), sig.GetName());
        printf("Adding my histogram class %s\n", nameStr2.Data());
        mcnames.push_back(nameStr2);
        histClasses += Form("%s;", nameStr2.Data());
      }
      fHistNamesMCMatched.push_back(mcnames);
    }

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, histClasses.Data());  // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TMuonFillMap, uint32_t TMuonMCFillMap, typename TEvent, typename TMuons, typename TEventsMC, typename TMuonsMC>
  void runSelection(TEvent const& event, TMuons const& muons, TEventsMC const& eventsMC, TMuonsMC const& muonsMC)
  {
    //cout << "Event ######################################" << endl;
    VarManager::ResetValues(0, VarManager::kNMCParticleVariables);
    // fill event information which might be needed in histograms that combine track and event properties
    VarManager::FillEvent<TEventFillMap>(event);
    if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
      VarManager::FillEvent<TEventMCFillMap>(event.reducedMCevent());
    }
    if constexpr ((TEventMCFillMap & VarManager::ObjTypes::CollisionMC) > 0) {
      VarManager::FillEvent<TEventMCFillMap>(event.mcCollision());
    }

    uint32_t filterMap = 0;
    muonSel.reserve(muons.size());
    for (auto& muon : muons) {
      filterMap = 0;
      VarManager::FillTrack<TMuonFillMap>(muon); // compute muon quantities

      // compute MC matched quantities using either the DQ skimmed or the Framework data models
      if constexpr ((TMuonFillMap & VarManager::ObjTypes::ReducedMuon) > 0) {
        VarManager::FillTrack<gkParticleMCFillMap>(muon.reducedMCTrack());
      }
      if constexpr ((TMuonFillMap & VarManager::ObjTypes::Muon) > 0) {
        VarManager::FillTrack<gkParticleMCFillMap>(muon.template mcParticle_as<aod::McParticles_001>());
      }

      if (fConfigQA) {
        fHistMan->FillHistClass("Muon_BeforeCuts", VarManager::fgValues);
      }

      // compute the cut selections and publish the filter bit map
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << i);
          if (fConfigQA) {
            fHistMan->FillHistClass(fHistNamesReco[i].Data(), VarManager::fgValues);
          }
        }
      }
      muonSel(static_cast<int>(filterMap));

      // if no filter fulfilled, continue
      if (!filterMap) {
        continue;
      }

      // everything below is related to filling QA histograms
      if (!fConfigQA) {
        continue;
      }

      // compute MC matching decisions
      uint32_t mcDecision = 0;
      int isig = 0;
      for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {
        if constexpr ((TMuonFillMap & VarManager::ObjTypes::ReducedMuon) > 0) {
          if ((*sig).CheckSignal(false, muonsMC, muon.reducedMCTrack())) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }
        if constexpr ((TMuonFillMap & VarManager::ObjTypes::Muon) > 0) {
          if ((*sig).CheckSignal(false, muonsMC, muon.template mcParticle_as<aod::McParticles_001>())) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }
      }

      // fill histograms
      for (unsigned int i = 0; i < fMCSignals.size(); i++) {
        if (!(mcDecision & (uint32_t(1) << i))) {
          continue;
        }
        for (unsigned int j = 0; j < fTrackCuts.size(); j++) {
          if (filterMap & (uint8_t(1) << j)) {
            fHistMan->FillHistClass(fHistNamesMCMatched[j][i].Data(), VarManager::fgValues);
          }
        } // end loop over cuts
      }   // end loop over MC signals
    }     // end loop over muons
  }

  void processSkimmed(MyEventsSelected::iterator const& event, MyMuonTracks const& muons, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runSelection<gkEventFillMap, gkMCEventFillMap, gkMuonFillMap, gkParticleMCFillMap>(event, muons, eventsMC, tracksMC);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisMuonSelection, processSkimmed, "Run muon selection on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisMuonSelection, processDummy, "Dummy process function", false);
};

struct AnalysisSameEventPairing {
  Produces<aod::Dileptons> dileptonList;
  Produces<aod::DileptonsExtra> dileptonExtraList;
  Produces<aod::DimuonsAll> dimuonAllList;
  OutputObj<THashList> fOutputList{"output"};
  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected == 1;
  Filter filterBarrelTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;
  Filter filterMuonSelected = aod::dqanalysisflags::isMuonSelected > 0;
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMCRecSignals{"cfgBarrelMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};
  Configurable<std::string> fConfigMCGenSignals{"cfgBarrelMCGenSignals", "", "Comma separated list of MC signals (generated)"};
  Configurable<bool> fConfigFlatTables{"cfgFlatTables", false, "Produce a single flat tables with all relevant information of the pairs and single tracks"};
  // TODO: here we specify signals, however signal decisions are precomputed and stored in mcReducedFlags
  // TODO: The tasks based on skimmed MC could/should rely ideally just on these flags
  // TODO:   special AnalysisCuts to be prepared in this direction
  // TODO: cuts on the MC truth information to be added if needed

  HistogramManager* fHistMan;
  std::vector<std::vector<TString>> fBarrelHistNames;
  std::vector<std::vector<TString>> fBarrelHistNamesMCmatched;
  std::vector<std::vector<TString>> fMuonHistNames;
  std::vector<std::vector<TString>> fMuonHistNamesMCmatched;
  std::vector<std::vector<TString>> fBarrelMuonHistNames;
  std::vector<std::vector<TString>> fBarrelMuonHistNamesMCmatched;
  std::vector<MCSignal> fRecMCSignals;
  std::vector<MCSignal> fGenMCSignals;

  void init(o2::framework::InitContext& context)
  {
    bool enableBarrelHistos = context.mOptions.get<bool>("processJpsiToEESkimmed");
    bool enableMuonHistos = context.mOptions.get<bool>("processJpsiToMuMuSkimmed") || context.mOptions.get<bool>("processJpsiToMuMuVertexingSkimmed");
    //bool enableBarrelMuonHistos = context.mOptions.get<bool>("processElectronMuonSkimmed");

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // For barrel-barrel, muon-muon or barrel-muon create:
    // a) 3 histogram classes with SEPM,SEPP and SEMM pairing
    // b) 1 histogram class for each specified MCsignal  (in total we have n X m histogram classes for each track and MCsignal combination)
    //    For the MC matching, for now we create histogram classes for just the PM pairs
    TString sigNamesStr = fConfigMCRecSignals.value;
    std::unique_ptr<TObjArray> objRecSigArray(sigNamesStr.Tokenize(","));
    TString histNames = "";

    // Setting the MC rec signal names
    for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() != 2) { // NOTE: 2-prong signals required
          continue;
        }
        fRecMCSignals.push_back(*sig);
      }
    }

    if (enableBarrelHistos) {
      TString cutNames = fConfigTrackCuts.value;
      if (!cutNames.IsNull()) {
        std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
        for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
          std::vector<TString> names = {
            Form("PairsBarrelSEPM_%s", objArray->At(icut)->GetName()),
            Form("PairsBarrelSEPP_%s", objArray->At(icut)->GetName()),
            Form("PairsBarrelSEMM_%s", objArray->At(icut)->GetName())};
          histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          fBarrelHistNames.push_back(names);
          std::vector<TString> mcSigClasses;
          if (!sigNamesStr.IsNull()) {
            for (auto& sig : fRecMCSignals) {
              TString histName = Form("PairsBarrelSEPM_%s_%s", objArray->At(icut)->GetName(), sig.GetName());
              histNames += Form("%s;", histName.Data());
              mcSigClasses.push_back(histName);
            } // end loop over MC signals
          }
          fBarrelHistNamesMCmatched.push_back(mcSigClasses);
        } // end loop over cuts
      }   // end if(cutNames.IsNull())
    }     // end if processBarrel

    if (enableMuonHistos) {
      TString cutNames = fConfigMuonCuts.value;
      if (!cutNames.IsNull()) {
        std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
        for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
          std::vector<TString> names = {
            Form("PairsMuonSEPM_%s", objArray->At(icut)->GetName()),
            Form("PairsMuonSEPP_%s", objArray->At(icut)->GetName()),
            Form("PairsMuonSEMM_%s", objArray->At(icut)->GetName())
          };
          histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          fMuonHistNames.push_back(names);
          std::vector<TString> mcSigClasses;
          if (!sigNamesStr.IsNull()) {
            for (auto& sig : fRecMCSignals) {
              TString histName = Form("PairsMuonSEPM_%s_%s", objArray->At(icut)->GetName(), sig.GetName());
              histNames += Form("%s;", histName.Data());
              mcSigClasses.push_back(histName);
            } // end loop over MC signals
          }
          fMuonHistNamesMCmatched.push_back(mcSigClasses);
        }  // end loop over cuts
      }    // end if(cutNames.IsNull())
    }      // end if processMuon

    // NOTE: For the electron-muon pairing, the policy is that the user specifies n track and n muon cuts via configurables
    //     So for each barrel cut there is a corresponding muon cut
    /*if (enableBarrelMuonHistos) {
      TString cutNamesBarrel = fConfigTrackCuts.value;
      TString cutNamesMuon = fConfigMuonCuts.value;
      if (!cutNamesBarrel.IsNull()) {
        std::unique_ptr<TObjArray> objArrayBarrel(cutNamesBarrel.Tokenize(","));
        std::unique_ptr<TObjArray> objArrayMuon(cutNamesMuon.Tokenize(","));
        for (int icut = 0; icut < objArrayBarrel->GetEntries(); ++icut) {
          if (icut >= objArrayMuon->GetEntries()) {
            // there are fewer muon cuts specified wrt barrel cuts
            break;
          }
          std::vector<TString> names = {
            Form("PairsEleMuSEPM_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName()),
            Form("PairsEleMuSEPP_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName()),
            Form("PairsEleMuSEMM_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName())
          };
          histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          fBarrelMuonHistNames.push_back(names);
          std::vector<TString> mcSigClasses;
          if (!sigNamesStr.IsNull()) {
            for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
              MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
              if (sig) {
                if (sig->GetNProngs() != 2) { // NOTE: 2-prong signals required
                  continue;
                }
                fRecMCSignals.push_back(*sig);
                TString histName = Form("PairsEleMuSEPM_%s_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName(), sig->GetName());
                histNames += Form("%s;", histName.Data());
                mcSigClasses.push_back(histName);
              }
            } // end loop over MC signals
          }
          fBarrelMuonHistNamesMCmatched.push_back(mcSigClasses);
        }  // end loop over cuts
      }  // end if(cutNames.IsNull())
    }  // end if processBarrelMuon
    */

    // Add histogram classes for each specified MCsignal at the generator level
    // TODO: create a std::vector of hist classes to be used at Fill time, to avoid using Form in the process function
    TString sigGenNamesStr = fConfigMCGenSignals.value;
    std::unique_ptr<TObjArray> objGenSigArray(sigGenNamesStr.Tokenize(","));
    for (int isig = 0; isig < objGenSigArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objGenSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() == 1) { // NOTE: 1-prong signals required
          fGenMCSignals.push_back(*sig);
          histNames += Form("MCTruthGen_%s;", sig->GetName()); // TODO: Add these names to a std::vector to avoid using Form in the process function
        } else if (sig->GetNProngs() == 2) { // NOTE: 2-prong signals required
          fGenMCSignals.push_back(*sig);
          histNames += Form("MCTruthGenPair_%s;", sig->GetName());
        }
      }
    }

    DefineHistograms(fHistMan, histNames.Data());    // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    VarManager::SetupTwoProngDCAFitter(5.0f, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, true); // TODO: get these parameters from Configurables
    VarManager::SetupTwoProngFwdDCAFitter(5.0f, true, 200.0f, 1.0e-3f, 0.9f, true);
  }

  template <int TPairType, uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks1, typename TTracks2, typename TEventsMC, typename TTracksMC>
  void runPairing(TEvent const& event, TTracks1 const& tracks1, TTracks2 const& tracks2, TEventsMC const& eventsMC, TTracksMC const& tracksMC)
  {
    // establish the right histogram classes to be filled depending on TPairType (ee,mumu,emu)
    unsigned int ncuts = fBarrelHistNames.size();
    std::vector<std::vector<TString>> histNames = fBarrelHistNames;
    std::vector<std::vector<TString>> histNamesMCmatched = fBarrelHistNamesMCmatched;
    if constexpr (TPairType == VarManager::kJpsiToMuMu) {
      ncuts = fMuonHistNames.size();
      histNames = fMuonHistNames;
      histNamesMCmatched = fMuonHistNamesMCmatched;
    }
    if constexpr (TPairType == VarManager::kElectronMuon) {
      ncuts = fBarrelMuonHistNames.size();
      histNames = fBarrelMuonHistNames;
      histNamesMCmatched = fBarrelMuonHistNamesMCmatched;
    }

    // Loop over two track combinations
    uint8_t twoTrackFilter = 0;
    uint32_t dileptonFilterMap = 0;
    uint32_t dileptonMcDecision = 0;
    dileptonList.reserve(1);
    dileptonExtraList.reserve(1);
    if (fConfigFlatTables.value) {
      dimuonAllList.reserve(1);
    }

    for (auto& [t1, t2] : combinations(tracks1, tracks2)) {
      if constexpr (TPairType == VarManager::kJpsiToEE) {
        twoTrackFilter = uint32_t(t1.isBarrelSelected()) & uint32_t(t2.isBarrelSelected());
      }
      if constexpr (TPairType == VarManager::kJpsiToMuMu) {
        twoTrackFilter = uint32_t(t1.isMuonSelected()) & uint32_t(t2.isMuonSelected());
      }
      if constexpr (TPairType == VarManager::kElectronMuon) {
        twoTrackFilter = uint32_t(t1.isBarrelSelected()) & uint32_t(t2.isMuonSelected());
      }
      if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }
      VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
      // secondary vertexing is not implemented for e-mu pairs so we need to hide this function from the e-mu analysis for now
      if constexpr ((TPairType == VarManager::kJpsiToEE) || (TPairType == VarManager::kJpsiToMuMu)) {
        VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, VarManager::fgValues);
      }

      // run MC matching for this pair
      uint32_t mcDecision = 0;
      int isig = 0;
      for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
        if constexpr (TTrackFillMap & VarManager::ObjTypes::ReducedTrack || TTrackFillMap & VarManager::ObjTypes::ReducedMuon) { // for skimmed DQ model
          if ((*sig).CheckSignal(false, tracksMC, t1.reducedMCTrack(), t2.reducedMCTrack())) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }
        if constexpr (TTrackFillMap & VarManager::ObjTypes::Track || TTrackFillMap & VarManager::ObjTypes::Muon) { // for Framework data model
          if ((*sig).CheckSignal(false, tracksMC, t1.template mcParticle_as<aod::McParticles_001>(), t2.template mcParticle_as<aod::McParticles_001>())) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }
      } // end loop over MC signals

      dileptonFilterMap = twoTrackFilter;
      dileptonMcDecision = mcDecision;
      dileptonList(event, VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), dileptonFilterMap, dileptonMcDecision);
      constexpr bool muonHasCov = ((TTrackFillMap & VarManager::ObjTypes::MuonCov) > 0 || (TTrackFillMap & VarManager::ObjTypes::ReducedMuonCov) > 0);
      if constexpr ((TPairType == VarManager::kJpsiToMuMu) && muonHasCov) {
        dileptonExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingLz], VarManager::fgValues[VarManager::kVertexingLxy]);
        if (fConfigFlatTables.value) {
          dimuonAllList(event.posX(), event.posY(), event.posZ(), event.reducedMCevent().mcPosX(), event.reducedMCevent().mcPosY(), event.reducedMCevent().mcPosZ(), VarManager::fgValues[VarManager::kMass], dileptonMcDecision, VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingTauzErr], VarManager::fgValues[VarManager::kVertexingTauxy], VarManager::fgValues[VarManager::kVertexingTauxyErr], t1.pt(), t1.eta(), t1.phi(), t1.sign(), t2.pt(), t2.eta(), t2.phi(), t2.sign(), t1.mcMask(), t2.mcMask(), t1.chi2MatchMCHMID(), t2.chi2MatchMCHMID(), t1.chi2MatchMCHMFT(), t2.chi2MatchMCHMFT(), t1.reducedMCTrack().pt(), t1.reducedMCTrack().eta(), t1.reducedMCTrack().phi(), t1.reducedMCTrack().e(), t2.reducedMCTrack().pt(), t2.reducedMCTrack().eta(), t2.reducedMCTrack().phi(), t2.reducedMCTrack().e(), t1.reducedMCTrack().vx(), t1.reducedMCTrack().vy(), t1.reducedMCTrack().vz(), t1.reducedMCTrack().vt(), t2.reducedMCTrack().vx(), t2.reducedMCTrack().vy(), t2.reducedMCTrack().vz(), t2.reducedMCTrack().vt());
        }
      }

      // Loop over all fulfilled cuts and fill pair histograms
      for (unsigned int icut = 0; icut < ncuts; icut++) {
        if (twoTrackFilter & (uint8_t(1) << icut)) {
          if (t1.sign() * t2.sign() < 0) {
            fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
            for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
              if (mcDecision & (uint32_t(1) << isig)) {
                fHistMan->FillHistClass(histNamesMCmatched[icut][isig].Data(), VarManager::fgValues);
              }
            }
          } else {
            if (t1.sign() > 0) {
              fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
            } else {
              fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
            }
          }
        }
      }
    } // end loop over barrel track pairs
  }   // end runPairing

  template <typename TTracksMC>
  void runMCGen(TTracksMC const& groupedMCTracks)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // group all the MC tracks which belong to the MC event corresponding to the current reconstructed event
    //auto groupedMCTracks = tracksMC.sliceBy(aod::reducedtrackMC::reducedMCeventId, event.reducedMCevent().globalIndex());
    for (auto& mctrack : groupedMCTracks) {
      VarManager::FillTrack<gkParticleMCFillMap>(mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (auto& sig : fGenMCSignals) {
        if (sig.GetNProngs() != 1) { // NOTE: 1-prong signals required
          continue;
        }
        if (sig.CheckSignal(false, groupedMCTracks, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig.GetName()), VarManager::fgValues);
        }
      }
    }

    //    // loop over mc stack and fill histograms for pure MC truth signals
    for (auto& sig : fGenMCSignals) {
      if (sig.GetNProngs() != 2) { // NOTE: 2-prong signals required
        continue;
      }
      for (auto& [t1, t2] : combinations(groupedMCTracks, groupedMCTracks)) {
        if (sig.CheckSignal(false, groupedMCTracks, t1, t2)) {
          VarManager::FillPairMC(t1, t2);
          fHistMan->FillHistClass(Form("MCTruthGenPair_%s", sig.GetName()), VarManager::fgValues);
        }
      }
    } //end of true pairing loop
  }   // end runMCGen

  Preslice<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;

  void processJpsiToEESkimmed(soa::Filtered<MyEventsSelected>::iterator const& event,
                              soa::Filtered<MyBarrelTracksSelected> const& tracks,
                              ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event);
    VarManager::FillEvent<gkMCEventFillMap>(event.reducedMCevent());

    runPairing<VarManager::kJpsiToEE, gkEventFillMap, gkMCEventFillMap, gkTrackFillMap>(event, tracks, tracks, eventsMC, tracksMC);
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    runMCGen(groupedMCTracks);
  }

  void processJpsiToMuMuSkimmed(soa::Filtered<MyEventsSelected>::iterator const& event,
                                soa::Filtered<MyMuonTracksSelected> const& muons,
                                ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event);
    VarManager::FillEvent<gkMCEventFillMap>(event.reducedMCevent());

    runPairing<VarManager::kJpsiToMuMu, gkEventFillMap, gkMCEventFillMap, gkMuonFillMap>(event, muons, muons, eventsMC, tracksMC);
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    runMCGen(groupedMCTracks);
  }

  void processJpsiToMuMuVertexingSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event,
                                         soa::Filtered<MyMuonTracksSelectedWithCov> const& muons,
                                         ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event);
    VarManager::FillEvent<gkMCEventFillMap>(event.reducedMCevent());

    runPairing<VarManager::kJpsiToMuMu, gkEventFillMapWithCov, gkMCEventFillMap, gkMuonFillMapWithCov>(event, muons, muons, eventsMC, tracksMC);
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    runMCGen(groupedMCTracks);
  }

  /*void processElectronMuonSkimmed(soa::Filtered<MyEventsSelected>::iterator const& event,
                                  soa::Filtered<MyBarrelTracksSelected> const& tracks, soa::Filtered<MyMuonTracksSelected> const& muons,
                                  ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC) {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event);
    VarManager::FillEvent<gkMCEventFillMap>(event.reducedMCevent());

    runPairing<VarManager::kElectronMuon, gkEventFillMap, gkMCEventFillMap, gkTrackFillMap>(event, tracks, muons, eventsMC, tracksMC);
    auto groupedMCTracks = tracksMC.sliceBy(aod::reducedtrackMC::reducedMCeventId, event.reducedMCevent().globalIndex());
    runMCGen(groupedMCTracks);
  }*/
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processJpsiToEESkimmed, "Run barrel barrel pairing on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processJpsiToMuMuSkimmed, "Run muon muon pairing on DQ skimmed muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processJpsiToMuMuVertexingSkimmed, "Run muon muon pairing on DQ skimmed muons including vertexing", false);
  //PROCESS_SWITCH(AnalysisSameEventPairing, processElectronMuonSkimmed, "Run barrel muon pairing on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy process function", false);
};

struct AnalysisDileptonTrack {
  Produces<aod::DileptonTrackCandidates> dileptontrackcandidatesList;
  OutputObj<THashList> fOutputList{"output"};
  // TODO: For now this is only used to determine the position in the filter bit map for the hadron cut
  Configurable<string> fConfigTrackCuts{"cfgLeptonCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<bool> fConfigFillCandidateTable{"cfgFillCandidateTable", false, "Produce a single flat tables with all relevant information dilepton-track candidates"};
  Filter eventFilter = aod::dqanalysisflags::isEventSelected == 1;
  //Filter dileptonFilter = aod::reducedpair::mass > 2.92f && aod::reducedpair::mass < 3.16f && aod::reducedpair::sign == 0;
  //Filter dileptonFilter = aod::reducedpair::mass > 2.6f && aod::reducedpair::mass < 3.5f && aod::reducedpair::sign == 0;

  Configurable<std::string> fConfigMCRecSignals{"cfgBarrelMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};
  Configurable<std::string> fConfigMCGenSignals{"cfgBarrelMCGenSignals", "", "Comma separated list of MC signals (generated)"};

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair; // fill map

  // use two values array to avoid mixing up the quantities
  float* fValuesDilepton;
  float* fValuesTrack;
  HistogramManager* fHistMan;

  std::vector<std::vector<TString>> fMuonHistNames;
  std::vector<std::vector<TString>> fMuonHistNamesMCmatched;
  std::vector<TString> fRecMCSignalsNames;

  std::vector<MCSignal> fRecMCSignals;
  std::vector<MCSignal> fGenMCSignals;

  // NOTE: the barrel track filter is shared between the filters for dilepton electron candidates (first n-bits)
  //       and the associated hadrons (n+1 bit) --> see the barrel track selection task
  //      The current condition should be replaced when bitwise operators will become available in Filter expresions
  int fNHadronCutBit;

  void init(o2::framework::InitContext& context)
  {
    TString sigNamesStr = fConfigMCRecSignals.value;
    std::unique_ptr<TObjArray> objRecSigArray(sigNamesStr.Tokenize(","));
    TString histNames;

    fValuesDilepton = new float[VarManager::kNVars];
    fValuesTrack = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // TODO: Create separate histogram directories for each selection used in the creation of the dileptons
    // TODO: Implement possibly multiple selections for the associated track ?
    if (context.mOptions.get<bool>("processDimuonMuonSkimmed")) {
      //DefineHistograms(fHistMan, "DileptonsSelected;DileptonTrackInvMass;DileptonsSelected_matchedMC;DileptonTrackInvMass_matchedMC;"); // define all histograms
      //VarManager::SetUseVars(fHistMan->GetUsedVars());
      //fOutputList.setObject(fHistMan->GetMainHistogramList());

      histNames += "DileptonsSelected;DileptonTrackInvMass;";

      if (!sigNamesStr.IsNull()) {
        for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
          MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
          if (sig) {
            if (sig->GetNProngs() == 1) {
              fRecMCSignals.push_back(*sig);
              TString histName = Form("LeptonsSelected_matchedMC_%s", sig->GetName());
              histNames += Form("%s;", histName.Data());
              fRecMCSignalsNames.push_back(sig->GetName());
            }
            if (sig->GetNProngs() == 2) {
              fRecMCSignals.push_back(*sig);
              TString histName = Form("DileptonsSelected_matchedMC_%s", sig->GetName());
              histNames += Form("%s;", histName.Data());
              fRecMCSignalsNames.push_back(sig->GetName());
            }
            if (sig->GetNProngs() == 3) {
              fRecMCSignals.push_back(*sig);
              TString histName = Form("DileptonTrackInvMass_matchedMC_%s", sig->GetName());
              histNames += Form("%s;", histName.Data());
              fRecMCSignalsNames.push_back(sig->GetName());
            }
          }
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
            fGenMCSignals.push_back(*sig);
            histNames += Form("MCTruthGen_%s;", sig->GetName()); // TODO: Add these names to a std::vector to avoid using Form in the process function
          }
        }
      }

      DefineHistograms(fHistMan, histNames.Data()); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    TString configCutNamesStr = fConfigTrackCuts.value;
    if (!configCutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(configCutNamesStr.Tokenize(","));
      fNHadronCutBit = objArray->GetEntries();
    } else {
      fNHadronCutBit = 0;
    }
  }

  // Template function to run pair - track combinations
  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks, typename TEventsMC, typename TTracksMC>
  void runDileptonTrack(TEvent const& event, TTracks const& tracks, soa::Join<aod::Dileptons, aod::DileptonsExtra> const& dileptons, TEventsMC const& eventsMC, TTracksMC const& tracksMC)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesTrack);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    VarManager::FillEvent<TEventFillMap>(event, fValuesTrack);
    VarManager::FillEvent<TEventFillMap>(event, fValuesDilepton);

    // Set the global index offset to find the proper lepton
    // TO DO: remove it once the issue with lepton index is solved
    int indexOffset = -999;
    for (auto dilepton : dileptons) {

      int indexLepton1 = dilepton.index0Id();
      int indexLepton2 = dilepton.index1Id();

      if (indexOffset == -999) {
        indexOffset = indexLepton1;
      }

      auto lepton1 = tracks.iteratorAt(indexLepton1 - indexOffset);
      auto lepton2 = tracks.iteratorAt(indexLepton2 - indexOffset);

      // Check that the dilepton has zero charge
      if (dilepton.sign() != 0) {
        continue;
      }

      // Check that the muons are opposite sign
      if (lepton1.sign() * lepton2.sign() > 0) {
        continue;
      }

      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesDilepton);
      fHistMan->FillHistClass("DileptonsSelected", fValuesDilepton);

      auto lepton1MC = lepton1.reducedMCTrack();
      auto lepton2MC = lepton2.reducedMCTrack();

      // run MC matching for this pair
      uint32_t mcDecision = 0;
      int isig = 0;
      for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
        if constexpr (TTrackFillMap & VarManager::ObjTypes::ReducedTrack || TTrackFillMap & VarManager::ObjTypes::ReducedMuon) { // for skimmed DQ model
          if ((*sig).CheckSignal(false, tracksMC, lepton1MC, lepton2MC)) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }
      } // end loop over MC signals

      for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
        if (mcDecision & (uint32_t(1) << isig)) {
          fHistMan->FillHistClass(Form("DileptonsSelected_matchedMC_%s", fRecMCSignalsNames[isig].Data()), fValuesDilepton);
        }
      }

      if (fConfigFillCandidateTable.value) {
        dileptontrackcandidatesList.reserve(1);
      }
      for (auto& track : tracks) {
        auto trackMC = track.reducedMCTrack();
        int index = track.globalIndex();

        // Remove combinations in which the track index is the same as the dilepton legs indices
        if (index == indexLepton1 || index == indexLepton2) {
          continue;
        }

        VarManager::FillDileptonTrackVertexing<TCandidateType, TEventFillMap, TTrackFillMap>(event, lepton1, lepton2, track, fValuesTrack);
        fHistMan->FillHistClass("DileptonTrackInvMass", fValuesTrack);

        mcDecision = 0;
        isig = 0;
        for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
          if constexpr (TTrackFillMap & VarManager::ObjTypes::ReducedTrack || TTrackFillMap & VarManager::ObjTypes::ReducedMuon || TTrackFillMap & VarManager::ObjTypes::ReducedMuon) { // for skimmed DQ model
            if ((*sig).CheckSignal(false, tracksMC, lepton1MC, lepton2MC, trackMC)) {
              mcDecision |= (uint32_t(1) << isig);
            }
          }
        }

        if (fConfigFillCandidateTable.value) {
          dileptontrackcandidatesList(mcDecision, fValuesTrack[VarManager::kPairMass], fValuesTrack[VarManager::kPairPt], fValuesTrack[VarManager::kPairEta], fValuesTrack[VarManager::kVertexingTauz], fValuesTrack[VarManager::kVertexingTauxy], fValuesTrack[VarManager::kVertexingLz], fValuesTrack[VarManager::kVertexingLxy]);
        }

        for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
          if (mcDecision & (uint32_t(1) << isig)) {
            fHistMan->FillHistClass(Form("DileptonTrackInvMass_matchedMC_%s", fRecMCSignalsNames[isig].Data()), fValuesTrack);
          }
        }
      }
    }
  }

  template <typename TTracksMC>
  void runMCGen(TTracksMC const& groupedMCTracks)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // group all the MC tracks which belong to the MC event corresponding to the current reconstructed event
    // auto groupedMCTracks = tracksMC.sliceBy(aod::reducedtrackMC::reducedMCeventId, event.reducedMCevent().globalIndex());
    for (auto& mctrack : groupedMCTracks) {
      VarManager::FillTrack<gkParticleMCFillMap>(mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (auto& sig : fGenMCSignals) {
        if (sig.GetNProngs() != 1) { // NOTE: 1-prong signals required
          continue;
        }
        if (sig.CheckSignal(false, groupedMCTracks, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig.GetName()), VarManager::fgValues);
        }
      }
    }
  }

  Preslice<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;

  void processDimuonMuonSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyMuonTracksSelectedWithCov const& tracks, soa::Join<aod::Dileptons, aod::DileptonsExtra> const& dileptons, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runDileptonTrack<VarManager::kBcToThreeMuons, gkEventFillMapWithCov, gkMCEventFillMap, gkMuonFillMapWithCov>(event, tracks, dileptons, eventsMC, tracksMC);
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    runMCGen(groupedMCTracks);
  }
  void processDielectronKaonSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyBarrelTracksSelectedWithCov const& tracks, soa::Join<aod::Dileptons, aod::DileptonsExtra> const& dileptons, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runDileptonTrack<VarManager::kBtoJpsiEEK, gkEventFillMapWithCov, gkMCEventFillMap, gkTrackFillMap>(event, tracks, dileptons, eventsMC, tracksMC);
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    runMCGen(groupedMCTracks);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisDileptonTrack, processDimuonMuonSkimmed, "Run dimuon-muon pairing, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processDielectronKaonSkimmed, "Run dielectron-kaon pairing, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
    adaptAnalysisTask<AnalysisMuonSelection>(cfgc),
    adaptAnalysisTask<AnalysisSameEventPairing>(cfgc),
    adaptAnalysisTask<AnalysisDileptonTrack>(cfgc)};
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
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "trigger,cent,mc");
    }

    if (classStr.Contains("Track")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "its,tpcpid,dca,tofpid,mc");
      }
    }

    if (classStr.Contains("Muon")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "muon");
    }

    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_barrel", "vertexing-barrel");
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_dimuon", "vertexing-forward");
    }

    if (classStr.Contains("MCTruthGenPair")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_pair");
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Pt", "MC generator p_{T} distribution", false, 200, 0.0, 20.0, VarManager::kMCPt);
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Eta", "MC generator #eta distribution", false, 500, -5.0, 5.0, VarManager::kMCEta);
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Phi", "MC generator #varphi distribution", false, 500, -6.3, 6.3, VarManager::kMCPhi);
    }
    if (classStr.Contains("MCTruthGen")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth");
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Pt", "MC generator p_{T} distribution", false, 200, 0.0, 20.0, VarManager::kMCPt);
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Eta", "MC generator #eta distribution", false, 500, -5.0, 5.0, VarManager::kMCEta);
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Phi", "MC generator #varphi distribution", false, 500, -6.3, 6.3, VarManager::kMCPhi);
    }
    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_barrel");
    }
    if (classStr.Contains("LeptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "muon");
    }
    if (classStr.Contains("DileptonTrackInvMass")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-track-mass");
    }

  } // end loop over histogram classes
}
