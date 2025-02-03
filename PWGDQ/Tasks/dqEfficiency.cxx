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
#include <iostream>
#include <vector>
#include <TMath.h>
#include <TH1F.h>
#include <THashList.h>
#include <TString.h>
#include <cstdio>
#include <string>
#include <memory>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"

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
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "EVENTCUTS", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "BARRELTRACKCUTS", dqanalysisflags::IsBarrelSelected);
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);
} // namespace o2::aod

// using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsMC>;
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
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
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
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;"); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                 // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, typename TEvent, typename TEventsMC>
  void runSelection(TEvent const& event, TEventsMC const& /*mcEvents*/)
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

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

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
  void runSelection(TEvent const& event, TTracks const& tracks, TEventsMC const& /*eventsMC*/, TTracksMC const& /*tracksMC*/)
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
          filterMap |= (static_cast<uint32_t>(1) << i);
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
          if ((*sig).CheckSignal(false, track.reducedMCTrack())) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
          }
        }
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
          if ((*sig).CheckSignal(false, track.template mcParticle_as<aod::McParticles_001>())) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
          }
        }
      }

      // fill histograms
      for (unsigned int i = 0; i < fMCSignals.size(); i++) {
        if (!(mcDecision & (static_cast<uint32_t>(1) << i))) {
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
  void processSkimmedWithCov(MyEventsSelected::iterator const& event, MyBarrelTracksWithCov const& tracks, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runSelection<gkEventFillMap, gkMCEventFillMap, gkTrackFillMapWithCov, gkParticleMCFillMap>(event, tracks, eventsMC, tracksMC);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmedWithCov, "Run barrel track selection on DQ skimmed tracks with covariance", false);
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

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

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
    TString histClasses = "TrackMuon_BeforeCuts;";
    for (auto& cut : fTrackCuts) {
      TString nameStr = Form("TrackMuon_%s", cut.GetName());
      fHistNamesReco.push_back(nameStr);
      histClasses += Form("%s;", nameStr.Data());
      std::vector<TString> mcnames;
      for (auto& sig : fMCSignals) {
        TString nameStr2 = Form("TrackMuon_%s_%s", cut.GetName(), sig.GetName());
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
  void runSelection(TEvent const& event, TMuons const& muons, TEventsMC const& /*eventsMC*/, TMuonsMC const& /*muonsMC*/)
  {
    // cout << "Event ######################################" << endl;
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
        fHistMan->FillHistClass("TrackMuon_BeforeCuts", VarManager::fgValues);
      }

      // compute the cut selections and publish the filter bit map
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (static_cast<uint32_t>(1) << i);
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
          if ((*sig).CheckSignal(false, muon.reducedMCTrack())) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
          }
        }
        if constexpr ((TMuonFillMap & VarManager::ObjTypes::Muon) > 0) {
          if ((*sig).CheckSignal(false, muon.template mcParticle_as<aod::McParticles_001>())) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
          }
        }
      }

      // fill histograms
      for (unsigned int i = 0; i < fMCSignals.size(); i++) {
        if (!(mcDecision & (static_cast<uint32_t>(1) << i))) {
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
  Produces<aod::Dielectrons> dielectronList;
  Produces<aod::Dimuons> dimuonList;
  Produces<aod::DielectronsExtra> dielectronExtraList;
  Produces<aod::DielectronsInfo> dielectronInfoList;
  Produces<aod::DimuonsExtra> dimuonExtraList;
  Produces<aod::DielectronsAll> dielectronAllList;
  Produces<aod::DimuonsAll> dimuonAllList;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  float mMagField = 0.0;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  OutputObj<THashList> fOutputList{"output"};
  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected == 1;
  Filter filterBarrelTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;
  Filter filterMuonSelected = aod::dqanalysisflags::isMuonSelected > 0;
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMCRecSignals{"cfgMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};
  Configurable<std::string> fConfigMCGenSignals{"cfgMCGenSignals", "", "Comma separated list of MC signals (generated)"};
  Configurable<bool> fConfigSkimSignalOnly{"fConfigSkimSignalOnly", false, "Configurable to select only matched candidates"};
  Configurable<bool> fConfigFlatTables{"cfgFlatTables", false, "Produce a single flat tables with all relevant information of the pairs and single tracks"};
  Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
  Configurable<bool> fUseRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> fConfigAmbiguousHist{"cfgAmbiHist", false, "Enable Ambiguous histograms for time association studies"};
  Configurable<bool> fUseAbsDCA{"cfgUseAbsDCA", false, "Use absolute DCA minimization instead of chi^2 minimization in secondary vertexing"};
  Configurable<bool> fPropToPCA{"cfgPropToPCA", false, "Propagate tracks to secondary vertex"};
  Configurable<bool> fCorrFullGeo{"cfgCorrFullGeo", false, "Use full geometry to correct for MCS effects in track propagation"};
  Configurable<bool> fNoCorr{"cfgNoCorrFwdProp", false, "Do not correct for MCS effects in track propagation"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

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
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    fCurrentRun = 0;

    ccdb->setURL(ccdburl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    if (fNoCorr) {
      VarManager::SetupFwdDCAFitterNoCorr();
    } else if (fCorrFullGeo) {
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
    } else {
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
      VarManager::SetupMatLUTFwdDCAFitter(lut);
    }

    bool enableBarrelHistos = context.mOptions.get<bool>("processDecayToEESkimmed") || context.mOptions.get<bool>("processDecayToEEVertexingSkimmed");
    bool enableMuonHistos = context.mOptions.get<bool>("processDecayToMuMuSkimmed") || context.mOptions.get<bool>("processDecayToMuMuVertexingSkimmed");
    // bool enableBarrelMuonHistos = context.mOptions.get<bool>("processElectronMuonSkimmed");

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
            Form("PairsMuonSEMM_%s", objArray->At(icut)->GetName())};
          if (fConfigAmbiguousHist) {
            histNames += Form("%s;%s;%s;%s_unambiguous;%s_unambiguous;%s_unambiguous;", names[0].Data(), names[1].Data(), names[2].Data(), names[0].Data(), names[1].Data(), names[2].Data());
          } else {
            histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          }
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
        } // end loop over cuts
      }   // end if(cutNames.IsNull())
    }     // end if processMuon

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
        } else if (sig->GetNProngs() == 2) {                   // NOTE: 2-prong signals required
          fGenMCSignals.push_back(*sig);
          histNames += Form("MCTruthGenPair_%s;", sig->GetName());
        }
      }
    }

    DefineHistograms(fHistMan, histNames.Data());    // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  template <int TPairType, uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks1, typename TTracks2, typename TEventsMC, typename TTracksMC>
  void runPairing(TEvent const& event, TTracks1 const& tracks1, TTracks2 const& tracks2, TEventsMC const& /*eventsMC*/, TTracksMC const& /*tracksMC*/)
  {
    if (fCurrentRun != event.runNumber()) {
      if (fUseRemoteField) {
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, event.timestamp());
        if (grpmag != nullptr) {
          mMagField = grpmag->getNominalL3Field();
        } else {
          LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", event.timestamp());
        }
        if (fConfigUseKFVertexing) {
          VarManager::SetupTwoProngKFParticle(mMagField);
        } else {
          VarManager::SetupTwoProngDCAFitter(mMagField, fPropToPCA.value, 200.0f, 4.0f, 1.0e-3f, 0.9f, fUseAbsDCA.value); // TODO: get these parameters from Configurables
          VarManager::SetupTwoProngFwdDCAFitter(mMagField, fPropToPCA.value, 200.0f, 1.0e-3f, 0.9f, fUseAbsDCA.value);
        }
      } else {
        if (fConfigUseKFVertexing) {
          VarManager::SetupTwoProngKFParticle(fConfigMagField);
        } else {
          VarManager::SetupTwoProngDCAFitter(fConfigMagField.value, fPropToPCA.value, 200.0f, 4.0f, 1.0e-3f, 0.9f, fUseAbsDCA.value); // TODO: get these parameters from Configurables
          VarManager::SetupTwoProngFwdDCAFitter(fConfigMagField.value, fPropToPCA.value, 200.0f, 1.0e-3f, 0.9f, fUseAbsDCA.value);
        }
      }
      fCurrentRun = event.runNumber();
    }

    // establish the right histogram classes to be filled depending on TPairType (ee,mumu,emu)
    unsigned int ncuts = fBarrelHistNames.size();
    std::vector<std::vector<TString>> histNames = fBarrelHistNames;
    std::vector<std::vector<TString>> histNamesMCmatched = fBarrelHistNamesMCmatched;
    if constexpr (TPairType == VarManager::kDecayToMuMu) {
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
    dielectronList.reserve(1);
    dielectronInfoList.reserve(1);
    dimuonList.reserve(1);
    dielectronExtraList.reserve(1);
    dimuonExtraList.reserve(1);
    if (fConfigFlatTables.value) {
      dielectronAllList.reserve(1);
      dimuonAllList.reserve(1);
    }

    for (auto& [t1, t2] : combinations(tracks1, tracks2)) {
      if constexpr (TPairType == VarManager::kDecayToEE) {
        twoTrackFilter = static_cast<uint32_t>(t1.isBarrelSelected()) & static_cast<uint32_t>(t2.isBarrelSelected());
      }
      if constexpr (TPairType == VarManager::kDecayToMuMu) {
        twoTrackFilter = static_cast<uint32_t>(t1.isMuonSelected()) & static_cast<uint32_t>(t2.isMuonSelected());
      }
      if constexpr (TPairType == VarManager::kElectronMuon) {
        twoTrackFilter = static_cast<uint32_t>(t1.isBarrelSelected()) & static_cast<uint32_t>(t2.isMuonSelected());
      }
      if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }
      VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
      // secondary vertexing is not implemented for e-mu pairs so we need to hide this function from the e-mu analysis for now
      if constexpr ((TPairType == VarManager::kDecayToEE) || (TPairType == VarManager::kDecayToMuMu)) {
        VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, fPropToPCA, VarManager::fgValues);
      }

      // run MC matching for this pair
      uint32_t mcDecision = 0;
      int isig = 0;
      for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
        if constexpr (TTrackFillMap & VarManager::ObjTypes::ReducedTrack || TTrackFillMap & VarManager::ObjTypes::ReducedMuon) { // for skimmed DQ model
          if ((*sig).CheckSignal(false, t1.reducedMCTrack(), t2.reducedMCTrack())) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
          }
        }
        if constexpr (TTrackFillMap & VarManager::ObjTypes::Track || TTrackFillMap & VarManager::ObjTypes::Muon) { // for Framework data model
          if ((*sig).CheckSignal(false, t1.template mcParticle_as<aod::McParticles_001>(), t2.template mcParticle_as<aod::McParticles_001>())) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
          }
        }
      } // end loop over MC signals

      dileptonFilterMap = twoTrackFilter;
      dileptonMcDecision = mcDecision;
      if (!fConfigSkimSignalOnly || (fConfigSkimSignalOnly && mcDecision > 0)) {
        if constexpr (TPairType == VarManager::kDecayToEE) {
          dielectronList(event, VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), dileptonFilterMap, dileptonMcDecision);
          if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrackCollInfo) > 0) {
            dielectronInfoList(t1.collisionId(), t1.trackId(), t2.trackId());
          }
          dielectronExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingLz], VarManager::fgValues[VarManager::kVertexingLxy]);
        }
        if constexpr ((TPairType == VarManager::kDecayToEE) && (TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelPID) > 0) {
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
      }
      if constexpr (TPairType == VarManager::kDecayToMuMu) {
        dimuonList(event, VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), dileptonFilterMap, dileptonMcDecision);
        dimuonExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingLz], VarManager::fgValues[VarManager::kVertexingLxy]);
      }

      constexpr bool muonHasCov = ((TTrackFillMap & VarManager::ObjTypes::MuonCov) > 0 || (TTrackFillMap & VarManager::ObjTypes::ReducedMuonCov) > 0);
      if constexpr ((TPairType == VarManager::kDecayToMuMu) && muonHasCov) {
        if (fConfigFlatTables.value) {
          dimuonAllList(event.posX(), event.posY(), event.posZ(), event.numContrib(),
                        event.selection_raw(), 0,
                        event.reducedMCevent().mcPosX(), event.reducedMCevent().mcPosY(), event.reducedMCevent().mcPosZ(),
                        VarManager::fgValues[VarManager::kMass],
                        dileptonMcDecision,
                        VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), VarManager::fgValues[VarManager::kVertexingChi2PCA],
                        VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingTauzErr],
                        VarManager::fgValues[VarManager::kVertexingTauxy], VarManager::fgValues[VarManager::kVertexingTauxyErr],
                        VarManager::fgValues[VarManager::kCosPointingAngle],
                        VarManager::fgValues[VarManager::kPt1], VarManager::fgValues[VarManager::kEta1], VarManager::fgValues[VarManager::kPhi1], t1.sign(),
                        VarManager::fgValues[VarManager::kPt2], VarManager::fgValues[VarManager::kEta2], VarManager::fgValues[VarManager::kPhi2], t2.sign(),
                        t1.fwdDcaX(), t1.fwdDcaY(), t2.fwdDcaX(), t2.fwdDcaY(),
                        t1.mcMask(), t2.mcMask(),
                        t1.chi2MatchMCHMID(), t2.chi2MatchMCHMID(),
                        t1.chi2MatchMCHMFT(), t2.chi2MatchMCHMFT(),
                        t1.chi2(), t2.chi2(),
                        t1.reducedMCTrack().pt(), t1.reducedMCTrack().eta(), t1.reducedMCTrack().phi(), t1.reducedMCTrack().e(),
                        t2.reducedMCTrack().pt(), t2.reducedMCTrack().eta(), t2.reducedMCTrack().phi(), t2.reducedMCTrack().e(),
                        t1.reducedMCTrack().vx(), t1.reducedMCTrack().vy(), t1.reducedMCTrack().vz(), t1.reducedMCTrack().vt(),
                        t2.reducedMCTrack().vx(), t2.reducedMCTrack().vy(), t2.reducedMCTrack().vz(), t2.reducedMCTrack().vt(),
                        t1.isAmbiguous(), t2.isAmbiguous(), -999., -999., -999., -999., -999., -999., -999., -999., -999.,
                        -999., -999., -999., VarManager::fgValues[VarManager::kVertexingPz],
                        VarManager::fgValues[VarManager::kVertexingSV]);
        }
      }

      // Loop over all fulfilled cuts and fill pair histograms
      for (unsigned int icut = 0; icut < ncuts; icut++) {
        if (twoTrackFilter & (uint8_t(1) << icut)) {
          if (t1.sign() * t2.sign() < 0) {
            fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
            if (fConfigAmbiguousHist && !(t1.isAmbiguous() || t2.isAmbiguous())) {
              fHistMan->FillHistClass(Form("%s_unambiguous", histNames[icut][0].Data()), VarManager::fgValues);
            }
            for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
              if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                fHistMan->FillHistClass(histNamesMCmatched[icut][isig].Data(), VarManager::fgValues);
                if (fConfigAmbiguousHist && !(t1.isAmbiguous() || t2.isAmbiguous())) {
                  fHistMan->FillHistClass(Form("%s_unambiguous", histNamesMCmatched[icut][isig].Data()), VarManager::fgValues);
                }
              }
            }
          } else {
            if (t1.sign() > 0) {
              fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
              if (fConfigAmbiguousHist && !(t1.isAmbiguous() || t2.isAmbiguous())) {
                fHistMan->FillHistClass(Form("%s_unambiguous", histNames[icut][1].Data()), VarManager::fgValues);
              }
            } else {
              fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
              if (fConfigAmbiguousHist && !(t1.isAmbiguous() || t2.isAmbiguous())) {
                fHistMan->FillHistClass(Form("%s_unambiguous", histNames[icut][2].Data()), VarManager::fgValues);
              }
            }
          }
        }
      }
    } // end loop over barrel track pairs
  }   // end runPairing

  template <typename TTracksMC>
  void runMCGen(TTracksMC& groupedMCTracks)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // group all the MC tracks which belong to the MC event corresponding to the current reconstructed event
    // auto groupedMCTracks = tracksMC.sliceBy(aod::reducedtrackMC::reducedMCeventId, event.reducedMCevent().globalIndex());
    for (auto& mctrack : groupedMCTracks) {
      VarManager::FillTrackMC(groupedMCTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (auto& sig : fGenMCSignals) {
        if (sig.GetNProngs() != 1) { // NOTE: 1-prong signals required
          continue;
        }
        bool checked = false;
        if constexpr (soa::is_soa_filtered_v<TTracksMC>) {
          auto mctrack_raw = groupedMCTracks.rawIteratorAt(mctrack.globalIndex());
          checked = sig.CheckSignal(false, mctrack_raw);
        } else {
          checked = sig.CheckSignal(false, mctrack);
        }
        if (checked) {
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
        bool checked = false;
        if constexpr (soa::is_soa_filtered_v<TTracksMC>) {
          auto t1_raw = groupedMCTracks.rawIteratorAt(t1.globalIndex());
          auto t2_raw = groupedMCTracks.rawIteratorAt(t2.globalIndex());
          checked = sig.CheckSignal(false, t1_raw, t2_raw);
        } else {
          checked = sig.CheckSignal(false, t1, t2);
        }
        if (checked) {
          VarManager::FillPairMC(t1, t2);
          fHistMan->FillHistClass(Form("MCTruthGenPair_%s", sig.GetName()), VarManager::fgValues);
        }
      }
    } // end of true pairing loop
  }   // end runMCGen

  // Preslice<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;
  PresliceUnsorted<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;

  void processDecayToEESkimmed(soa::Filtered<MyEventsSelected>::iterator const& event,
                               soa::Filtered<MyBarrelTracksSelected> const& tracks,
                               ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event);
    VarManager::FillEvent<gkMCEventFillMap>(event.reducedMCevent());

    runPairing<VarManager::kDecayToEE, gkEventFillMap, gkMCEventFillMap, gkTrackFillMap>(event, tracks, tracks, eventsMC, tracksMC);
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    runMCGen(groupedMCTracks);
  }

  void processDecayToEEVertexingSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event,
                                        soa::Filtered<MyBarrelTracksSelectedWithCov> const& tracks,
                                        ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCov>(event);
    VarManager::FillEvent<gkMCEventFillMap>(event.reducedMCevent());

    runPairing<VarManager::kDecayToEE, gkEventFillMapWithCov, gkMCEventFillMap, gkTrackFillMapWithCov>(event, tracks, tracks, eventsMC, tracksMC);
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    runMCGen(groupedMCTracks);
  }

  void processDecayToMuMuSkimmed(soa::Filtered<MyEventsSelected>::iterator const& event,
                                 soa::Filtered<MyMuonTracksSelected> const& muons,
                                 ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event);
    VarManager::FillEvent<gkMCEventFillMap>(event.reducedMCevent());

    runPairing<VarManager::kDecayToMuMu, gkEventFillMap, gkMCEventFillMap, gkMuonFillMap>(event, muons, muons, eventsMC, tracksMC);
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    runMCGen(groupedMCTracks);
  }

  void processDecayToMuMuVertexingSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event,
                                          soa::Filtered<MyMuonTracksSelectedWithCov> const& muons,
                                          ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event);
    VarManager::FillEvent<gkMCEventFillMap>(event.reducedMCevent());

    runPairing<VarManager::kDecayToMuMu, gkEventFillMapWithCov, gkMCEventFillMap, gkMuonFillMapWithCov>(event, muons, muons, eventsMC, tracksMC);
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

  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmed, "Run barrel barrel pairing on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEEVertexingSkimmed, "Run barrel barrel pairing on DQ skimmed tracks including vertexing", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToMuMuSkimmed, "Run muon muon pairing on DQ skimmed muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToMuMuVertexingSkimmed, "Run muon muon pairing on DQ skimmed muons including vertexing", false);
  // PROCESS_SWITCH(AnalysisSameEventPairing, processElectronMuonSkimmed, "Run barrel muon pairing on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy process function", false);
};

struct AnalysisDileptonTrack {
  Produces<aod::DileptonTrackCandidates> dileptontrackcandidatesList;
  OutputObj<THashList> fOutputList{"output"};
  // TODO: For now this is only used to determine the position in the filter bit map for the hadron cut
  Configurable<string> fConfigTrackCuts{"cfgLeptonCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<bool> fConfigFillCandidateTable{"cfgFillCandidateTable", false, "Produce a single flat tables with all relevant information dilepton-track candidates"};
  Filter eventFilter = aod::dqanalysisflags::isEventSelected == 1;
  // Filter dileptonFilter = aod::reducedpair::mass > 2.92f && aod::reducedpair::mass < 3.16f && aod::reducedpair::sign == 0;
  // Filter dileptonFilter = aod::reducedpair::mass > 2.6f && aod::reducedpair::mass < 3.5f && aod::reducedpair::sign == 0;

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
  //      The current condition should be replaced when bitwise operators will become available in Filter expressions
  int fNHadronCutBit;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

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
    if (context.mOptions.get<bool>("processDimuonMuonSkimmed") || context.mOptions.get<bool>("processDielectronKaonSkimmed")) {
      // DefineHistograms(fHistMan, "DileptonsSelected;DileptonTrackInvMass;DileptonsSelected_matchedMC;DileptonTrackInvMass_matchedMC;"); // define all histograms
      // VarManager::SetUseVars(fHistMan->GetUsedVars());
      // fOutputList.setObject(fHistMan->GetMainHistogramList());

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
  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks, typename TDileptons, typename TEventsMC, typename TTracksMC>
  void runDileptonTrack(TEvent const& event, TTracks const& tracks, TDileptons const& dileptons, TEventsMC const& /*eventsMC*/, TTracksMC const& /*tracksMC*/)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesTrack);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    VarManager::FillEvent<TEventFillMap>(event, fValuesTrack);
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
    for (auto dilepton : dileptons) {

      int indexLepton1 = dilepton.index0Id();
      int indexLepton2 = dilepton.index1Id();

      if (indexOffset == -999) {
        indexOffset = trackGlobalIndexes.at(0);
      }
      trackGlobalIndexes.clear();

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
          if ((*sig).CheckSignal(false, lepton1MC, lepton2MC)) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
          }
        }
      } // end loop over MC signals

      for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
        if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
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

        VarManager::FillDileptonHadron(dilepton, track, fValuesTrack);
        VarManager::FillDileptonTrackVertexing<TCandidateType, TEventFillMap, TTrackFillMap>(event, lepton1, lepton2, track, fValuesTrack);
        fHistMan->FillHistClass("DileptonTrackInvMass", fValuesTrack);

        mcDecision = 0;
        isig = 0;
        for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
          if constexpr (TTrackFillMap & VarManager::ObjTypes::ReducedTrack || TTrackFillMap & VarManager::ObjTypes::ReducedMuon || TTrackFillMap & VarManager::ObjTypes::ReducedMuon) { // for skimmed DQ model
            if ((*sig).CheckSignal(false, lepton1MC, lepton2MC, trackMC)) {
              mcDecision |= (static_cast<uint32_t>(1) << isig);
            }
          }
        }

        if (fConfigFillCandidateTable.value) {
          dileptontrackcandidatesList(mcDecision, fValuesTrack[VarManager::kPairMass], fValuesTrack[VarManager::kPairPt], fValuesTrack[VarManager::kPairEta], fValuesTrack[VarManager::kVertexingTauz], fValuesTrack[VarManager::kVertexingTauxy], fValuesTrack[VarManager::kVertexingLz], fValuesTrack[VarManager::kVertexingLxy]);
        }

        for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
          if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
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
      VarManager::FillTrackMC(groupedMCTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (auto& sig : fGenMCSignals) {
        if (sig.GetNProngs() != 1) { // NOTE: 1-prong signals required
          continue;
        }
        bool checked = false;
        if constexpr (soa::is_soa_filtered_v<TTracksMC>) {
          auto mctrack_raw = groupedMCTracks.rawIteratorAt(mctrack.globalIndex());
          checked = sig.CheckSignal(false, mctrack_raw);
        } else {
          checked = sig.CheckSignal(false, mctrack);
        }
        if (checked) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig.GetName()), VarManager::fgValues);
        }
      }
    }
  }

  // Preslice<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;
  PresliceUnsorted<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;

  void processDimuonMuonSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyMuonTracksSelectedWithCov const& tracks, soa::Join<aod::Dimuons, aod::DimuonsExtra> const& dileptons, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runDileptonTrack<VarManager::kBcToThreeMuons, gkEventFillMapWithCov, gkMCEventFillMap, gkMuonFillMapWithCov>(event, tracks, dileptons, eventsMC, tracksMC);
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    runMCGen(groupedMCTracks);
  }
  void processDielectronKaonSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyBarrelTracksSelectedWithCov const& tracks, soa::Join<aod::Dielectrons, aod::DielectronsExtra> const& dileptons, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runDileptonTrack<VarManager::kBtoJpsiEEK, gkEventFillMapWithCov, gkMCEventFillMap, gkTrackFillMapWithCov>(event, tracks, dileptons, eventsMC, tracksMC);
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

struct AnalysisDileptonTrackTrack {
  OutputObj<THashList> fOutputList{"output"};

  Configurable<std::string> fConfigTrackCut1{"cfgTrackCut1", "pionPIDCut1", "track1 cut"}; // used for select the tracks from SelectedTracks
  Configurable<std::string> fConfigTrackCut2{"cfgTrackCut2", "pionPIDCut2", "track2 cut"}; // used for select the tracks from SelectedTracks
  Configurable<std::string> fConfigDileptonCut{"cfgDileptonCut", "pairJpsi2", "Dilepton cut"};
  Configurable<std::string> fConfigQuadrupletCuts{"cfgQuadrupletCuts", "pairX3872Cut1", "Comma separated list of Dilepton-Track-Track cut"};
  Configurable<std::string> fConfigMCRecSignals{"cfgBarrelMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};
  Configurable<std::string> fConfigMCGenSignals{"cfgBarrelMCGenSignals", "", "Comma separated list of MC signals (generated)"};
  Configurable<std::string> fConfigDileptonMCRecSignal{"cfgDileptonMCRecSignal", "", "Comma separated list of MC signals (reconstructed)"};

  Produces<aod::DileptonTrackTrackCandidates> DileptonTrackTrackTable;
  HistogramManager* fHistMan;

  std::vector<TString> fRecMCSignalsNames;
  std::vector<MCSignal> fRecMCSignals;
  std::vector<MCSignal> fGenMCSignals;
  std::vector<MCSignal> fDileptonMCRecSignals;

  Filter eventFilter = aod::dqanalysisflags::isEventSelected == 1;
  Filter dileptonFilter = aod::reducedpair::sign == 0;
  Filter filterBarrelTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;

  float* fValuesQuadruplet;

  std::vector<TString> fQuadrupletCutNames;
  AnalysisCompositeCut fDileptonCut;
  std::vector<AnalysisCompositeCut> fQuadrupletCuts;
  TString fTrackCutName1;
  TString fTrackCutName2;
  bool fIsSameTrackCut = false;

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair; // fill map

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

    TString histNames;

    // check if the same track cuts are used for both tracks
    if (fConfigTrackCut1.value == fConfigTrackCut2.value) {
      fIsSameTrackCut = true;
    }
    fTrackCutName1 = fConfigTrackCut1.value;
    fTrackCutName2 = fConfigTrackCut2.value;

    // dilepton cut
    TString configDileptonCutNamesStr = fConfigDileptonCut.value;
    fDileptonCut = *dqcuts::GetCompositeCut(configDileptonCutNamesStr.Data());

    // dilepton-track-track cuts
    TString configQuadruletCutNamesStr = fConfigQuadrupletCuts.value;
    std::unique_ptr<TObjArray> objArray(configQuadruletCutNamesStr.Tokenize(","));
    for (Int_t icut = 0; icut < objArray->GetEntries(); ++icut) {
      TString cutName = objArray->At(icut)->GetName();
      fQuadrupletCutNames.push_back(cutName);
      fQuadrupletCuts.push_back(*dqcuts::GetCompositeCut(cutName.Data()));
    }

    // reconstructible MC signals
    TString sigNamesStr = fConfigMCRecSignals.value;
    std::unique_ptr<TObjArray> objRecSigArray(sigNamesStr.Tokenize(","));
    for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() == 3) {
          fRecMCSignals.push_back(*sig);
        }
      }
    }

    // fill the histogram names
    if (!context.mOptions.get<bool>("processDummy")) {
      // Title_DileptonTrackTrackCutName
      if (!configQuadruletCutNamesStr.IsNull()) {
        for (std::size_t icut = 0; icut < fQuadrupletCutNames.size(); ++icut) {
          if (fIsSameTrackCut) {
            histNames += Form("QuadrupletSEPM_%s;", fQuadrupletCutNames[icut].Data());
          } else {
            histNames += Form("QuadrupletSEPM_%s;", fQuadrupletCutNames[icut].Data());
            histNames += Form("QuadrupletSEMP_%s;", fQuadrupletCutNames[icut].Data());
          }
          histNames += Form("QuadrupletSEPP_%s;", fQuadrupletCutNames[icut].Data());
          histNames += Form("QuadrupletSEMM_%s;", fQuadrupletCutNames[icut].Data());
          // Reco MC signals
          for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
            MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
            if (sig) {
              if (sig->GetNProngs() == 4) {
                fRecMCSignals.push_back(*sig);
                fRecMCSignalsNames.push_back(sig->GetName());
                histNames += Form("MCTruthRecQaud_%s_%s;", fQuadrupletCutNames[icut].Data(), sig->GetName());
              }
            }
          }
        } // loop over dilepton-track-track cuts
      }
    }

    // Genarate MC signals
    TString sigGenNamesStr = fConfigMCGenSignals.value;
    std::unique_ptr<TObjArray> objGenSigArray(sigGenNamesStr.Tokenize(","));
    for (int isig = 0; isig < objGenSigArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objGenSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() == 1) { // NOTE: 1-prong signals required
          fGenMCSignals.push_back(*sig);
          histNames += Form("MCTruthGenQaud_%s;", sig->GetName()); // TODO: Add these names to a std::vector to avoid using Form in the process function
        }
      }
    }

    DefineHistograms(fHistMan, histNames.Data()); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    // dilepton MC signal
    TString configDileptonMCRecSignalStr = fConfigDileptonMCRecSignal.value;
    std::unique_ptr<TObjArray> objDileptonMCRecSignalArray(configDileptonMCRecSignalStr.Tokenize(","));
    for (int isig = 0; isig < objDileptonMCRecSignalArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objDileptonMCRecSignalArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() == 2) {
          fDileptonMCRecSignals.push_back(*sig);
        }
      }
    }
  }

  // Template function to run pair - track - track combinations
  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks, typename TDileptons, typename TEventsMC, typename TTracksMC>
  void runDileptonTrackTrack(TEvent const& event, TTracks const& tracks, TDileptons const& dileptons, TEventsMC const& /*eventsMC*/, TTracksMC const& /*tracksMC*/)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesQuadruplet);
    VarManager::FillEvent<TEventFillMap>(event, fValuesQuadruplet);

    // LOGF(info, "Number of dileptons: %d", dileptons.size());
    int indexOffset = -999;
    std::vector<int> trackGlobalIndexes;

    if (dileptons.size() > 0) {
      for (auto track : tracks) {
        trackGlobalIndexes.push_back(track.globalIndex());
        // std::cout << track.index() << " " << track.globalIndex() << std::endl;
      }
    }

    // loop over dileptons
    for (auto dilepton : dileptons) {
      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesQuadruplet);

      // apply the dilepton cut
      if (!fDileptonCut.IsSelected(fValuesQuadruplet))
        continue;

      // get the index of the electron legs
      int indexLepton1 = dilepton.index0Id();
      int indexLepton2 = dilepton.index1Id();

      if (indexOffset == -999) {
        indexOffset = trackGlobalIndexes.at(0);
      }
      trackGlobalIndexes.clear();

      auto lepton1 = tracks.iteratorAt(indexLepton1 - indexOffset);
      auto lepton2 = tracks.iteratorAt(indexLepton2 - indexOffset);

      auto lepton1MC = lepton1.reducedMCTrack();
      auto lepton2MC = lepton2.reducedMCTrack();

      // loop over track - track combinations
      for (auto& [t1, t2] : combinations(tracks, tracks)) {
        // avoid self-combinations
        if (t1.globalIndex() == indexLepton1 || t1.globalIndex() == indexLepton2 || t2.globalIndex() == indexLepton1 || t2.globalIndex() == indexLepton2) {
          // LOGF(info, "self-combination: %d=%d %d=%d", t1.globalIndex(), indexLepton1, indexLepton2, t2.globalIndex());
          continue;
        }

        // dilepton combinate with two same particles
        if ((fIsSameTrackCut && (t1.isBarrelSelected() & (static_cast<uint32_t>(1) << 1)) && (t2.isBarrelSelected() & (static_cast<uint32_t>(1) << 1))) ||
            (!fIsSameTrackCut && (t1.isBarrelSelected() & (static_cast<uint32_t>(1) << 1)) && (t2.isBarrelSelected() & (static_cast<uint32_t>(1) << 2)))) {
        } else {
          continue;
        }

        // Fill the Histograms
        VarManager::FillDileptonTrackTrack<TCandidateType>(dilepton, t1, t2, fValuesQuadruplet);
        uint32_t CutDecision = 0;
        uint32_t mcDecision = 0;
        int iCut = 0;
        for (auto cut = fQuadrupletCuts.begin(); cut != fQuadrupletCuts.end(); cut++, iCut++) {
          if (cut->IsSelected(fValuesQuadruplet)) {
            CutDecision |= (static_cast<uint32_t>(1) << iCut);
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

            // Reco MC signals
            if (fRecMCSignals.size() > 0) {
              int isig = 0;
              for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
                if ((*sig).CheckSignal(true, lepton1MC, lepton2MC, t1.reducedMCTrack(), t2.reducedMCTrack())) {
                  mcDecision |= (static_cast<uint32_t>(1) << isig);
                }
              }
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
                if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                  fHistMan->FillHistClass(Form("MCTruthRecQaud_%s_%s", fQuadrupletCutNames[iCut].Data(), fRecMCSignalsNames[isig].Data()), fValuesQuadruplet);
                }
              }
            }
          }
        } // end loop over cuts

        // Fill the DileptonTrackTrackCandidates table
        if (!CutDecision)
          continue;
        if (!mcDecision)
          continue;
        DileptonTrackTrackTable(fValuesQuadruplet[VarManager::kQuadMass], fValuesQuadruplet[VarManager::kQuadPt], fValuesQuadruplet[VarManager::kQuadEta], fValuesQuadruplet[VarManager::kQuadPhi], fValuesQuadruplet[VarManager::kRap],
                                fValuesQuadruplet[VarManager::kQ], fValuesQuadruplet[VarManager::kDeltaR1], fValuesQuadruplet[VarManager::kDeltaR2], fValuesQuadruplet[VarManager::kDeltaR],
                                dilepton.mass(), dilepton.pt(), dilepton.eta(), dilepton.phi(), dilepton.sign(),
                                fValuesQuadruplet[VarManager::kDitrackMass], fValuesQuadruplet[VarManager::kDitrackPt], t1.pt(), t2.pt(), t1.eta(), t2.eta(), t1.phi(), t2.phi(), t1.sign(), t2.sign());
      } // end loop over track - track combinations
    } // end loop over dileptons
  };

  template <int TCandidateType>
  void runMCGen(ReducedMCTracks const& mcTracks)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    for (auto& track : mcTracks) {
      VarManager::FillTrackMC(mcTracks, track);
      for (auto& sig : fGenMCSignals) {
        if (sig.CheckSignal(true, track)) {
          int daughterIdFirst = track.daughtersIds()[0];
          int daughterIdEnd = track.daughtersIds()[1];
          int Ndaughters = daughterIdEnd - daughterIdFirst + 1;
          if (Ndaughters == 3) {
            auto dilepton = mcTracks.rawIteratorAt(daughterIdFirst);
            auto track1 = mcTracks.rawIteratorAt(daughterIdFirst + 1);
            auto track2 = mcTracks.rawIteratorAt(daughterIdFirst + 2);
            VarManager::FillQaudMC<TCandidateType>(dilepton, track1, track2);
          }
          fHistMan->FillHistClass(Form("MCTruthGenQaud_%s", sig.GetName()), VarManager::fgValues);
        }
      }
    }
  };

  PresliceUnsorted<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;

  void processX3872(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyBarrelTracksSelectedWithCov const& tracks, soa::Join<aod::Dielectrons, aod::DielectronsExtra> const& dileptons, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runDileptonTrackTrack<VarManager::kXtoJpsiPiPi, gkEventFillMapWithCov, gkMCEventFillMap, gkTrackFillMapWithCov>(event, tracks, dileptons, eventsMC, tracksMC);
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    runMCGen<VarManager::kXtoJpsiPiPi>(groupedMCTracks);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisDileptonTrackTrack, processX3872, "Run X3872 reconstruction", false);
  PROCESS_SWITCH(AnalysisDileptonTrackTrack, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
    adaptAnalysisTask<AnalysisMuonSelection>(cfgc),
    adaptAnalysisTask<AnalysisSameEventPairing>(cfgc),
    adaptAnalysisTask<AnalysisDileptonTrack>(cfgc),
    adaptAnalysisTask<AnalysisDileptonTrackTrack>(cfgc)};
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

    if (classStr.Contains("Track") && !classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "its,tpcpid,dca,tofpid,mc,mcMother");
      }
      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "muon");
      }
    }

    if (classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "barrel,vertexing-barrel,kalman-filter");
      }
      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "dimuon,dimuon-multi-diff,vertexing-forward");
      }
    }

    if (classStr.Contains("MCTruthGenPair")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_pair");
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Pt_Rapidity", "MC generator p_{T}, y distribution", false, 120, 0.0, 30.0, VarManager::kMCPt, 150, 2.5, 4.0, VarManager::kMCY);
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Eta", "MC generator #eta distribution", false, 200, 2.5, 4.0, VarManager::kMCEta);
      // histMan->AddHistogram(objArray->At(iclass)->GetName(), "Rapidity", "MC generator y distribution", false, 150, 2.5, 4.0, VarManager::kMCY);
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Phi", "MC generator #varphi distribution", false, 50, 0.0, 2. * TMath::Pi(), VarManager::kMCPhi);
    }
    if (classStr.Contains("MCTruthGenQaud")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_quad");
    }
    if (classStr.Contains("MCTruthGen")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth");
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Pt_Rapidity", "MC generator p_{T}, y distribution", false, 120, 0.0, 30.0, VarManager::kMCPt, 150, 2.5, 4.0, VarManager::kMCY);
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Eta", "MC generator #eta distribution", false, 200, 2.5, 4.0, VarManager::kMCEta);
      // histMan->AddHistogram(objArray->At(iclass)->GetName(), "Rapidity", "MC generator y distribution", false, 150, 2.5, 4.0, VarManager::kMCY);
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Phi", "MC generator #varphi distribution", false, 50, 0.0, 2. * TMath::Pi(), VarManager::kMCPhi);
    }
    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "barrel,dimuon");
    }
    if (classStr.Contains("LeptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "muon");
    }
    if (classStr.Contains("DileptonTrackInvMass")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-mass");
    }
    if (classStr.Contains("Quadruplet") || classStr.Contains("MCTruthRecQaud")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-dihadron", "xtojpsipipi");
    }

  } // end loop over histogram classes
}
