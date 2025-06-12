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
// Contact: florian.eisenhut@cern.ch
//
// Analysis task to generate Monte Carlo templates of different heavy-flavour dielectron sources
//
#include <iostream>
#include <vector>
#include <TMath.h>
#include <TH1F.h>
#include <THashList.h>
#include <TString.h>
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
#include "Common/DataModel/TrackSelectionTables.h"

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

namespace emanalysisflags
{
DECLARE_SOA_COLUMN(IsMCEventSelected, isMCEventSelected, int);
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelected, isBarrelSelected, int);
} // namespace emanalysisflags

DECLARE_SOA_TABLE(EventMCCuts, "AOD", "EVENTMCCUTS", emanalysisflags::IsMCEventSelected);
DECLARE_SOA_TABLE(EventCuts, "AOD", "EVENTCUTS", emanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "BARRELTRACKCUTS", emanalysisflags::IsBarrelSelected);
} // namespace o2::aod

// No skimming: works for events and single tracks
using MyEventsAOD = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using MyEventsSelectedAOD = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::EventCuts>;
// using MyMCEventsSelectedAOD = soa::Join<aod::McCollisions, aod::EventMCCuts>;
using MyBarrelTracksAOD = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                    aod::pidTPCFullKa, aod::pidTPCFullPr,
                                    aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                    aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                    aod::McTrackLabels>;
using MyBarrelTracksSelectedAOD = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection,
                                            aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                            aod::pidTPCFullKa, aod::pidTPCFullPr,
                                            aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                            aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                            aod::BarrelTrackCuts, aod::McTrackLabels>;
// using MyMCTrackAOD = soa::Join<aod::McParticles, aod::SmearedTracks>;

constexpr static uint32_t gkEventFillMapAOD = VarManager::ObjTypes::Collision;
constexpr static uint32_t gkMCEventFillMapAOD = VarManager::ObjTypes::CollisionMC;
constexpr static uint32_t gkTrackFillMapAOD = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkParticleMCFillMapAOD = VarManager::ObjTypes::ParticleMC;

// Skimmed data
// using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsMC>;
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedMCEventLabels>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::ReducedMCEventLabels>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedMCEventLabels>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::ReducedTracksBarrelLabels>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::ReducedTracksBarrelLabels>;
using MyBarrelTracksSelectedWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::ReducedTracksBarrelLabels>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkMCEventFillMap = VarManager::ObjTypes::ReducedEventMC;
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkParticleMCFillMap = VarManager::ObjTypes::ParticleMC;

void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar); // defines histograms for all tasks

struct AnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  OutputObj<THashList> fOutputList{"output"};
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
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
      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;", fConfigAddEventHistogram); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                           // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, typename TEvent, typename TEventsMC>
  void runEventSelection(TEvent const& event, TEventsMC const&)
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
    runEventSelection<gkEventFillMap, gkMCEventFillMap>(event, mcEvents);
  }

  void processAOD(MyEventsAOD::iterator const& event, aod::McCollisions const& mcEvents)
  {
    runEventSelection<gkEventFillMapAOD, gkMCEventFillMapAOD>(event, mcEvents);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }
  void processDummyAOD(MyEventsAOD&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventSelection, processAOD, "Run event selection without skimming", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDummy, "Dummy process function", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDummyAOD, "Dummy process function", false);
};

struct AnalysisTrackSelection {
  Produces<aod::BarrelTrackCuts> trackSel;
  Filter filterEventSelected = aod::emanalysisflags::isEventSelected == 1;
  Filter filterMCEventSelected = aod::emanalysisflags::isMCEventSelected == 1;
  OutputObj<THashList> fOutputList{"output"};
  Configurable<std::string> fConfigCuts{"cfgTrackCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMCSignals{"cfgTrackMCSignals", "", "Comma separated list of MC signals"};
  Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<bool> fConfigMCTruthGen{"cfgMCTruthGen", false, "If true, fill MCTruthGen histograms"};

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
        fTrackCuts.emplace_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
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
        fMCSignals.emplace_back(*sig);
      }
    }

    // Configure histogram classes for each track cut;
    // Add histogram classes for each track cut and for each requested MC signal (reconstructed tracks with MC truth)
    TString histClasses = "TrackBarrel_BeforeCuts;";
    for (auto& cut : fTrackCuts) {
      TString nameStr = Form("TrackBarrel_%s", cut.GetName());
      fHistNamesReco.emplace_back(nameStr);
      histClasses += Form("%s;", nameStr.Data());
      std::vector<TString> mcnames;
      for (auto& sig : fMCSignals) {
        TString nameStr2 = Form("TrackBarrel_%s_%s", cut.GetName(), sig.GetName());
        printf("Adding my histogram class %s\n", nameStr2.Data());
        mcnames.emplace_back(nameStr2);
        histClasses += Form("%s;", nameStr2.Data());
      }
      fHistNamesMCMatched.emplace_back(mcnames);
    }

    if (fConfigQA) {
      if (fConfigMCTruthGen) {
        // Add histogram classes for each MC signal at generated level
        // std::vector<TString> mcnamesgen;
        for (int isig = 0; isig < sigNamesArray->GetEntries(); ++isig) {
          MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(sigNamesArray->At(isig)->GetName());
          if (sig) {
            if (sig->GetNProngs() == 1) { // NOTE: only 1 prong signals
              TString nameStr2 = Form("MCTruthGenTrack_%s", sig->GetName());
              // mcnamesgen.push_back(nameStr2);
              histClasses += Form("%s;", nameStr2.Data()); // TODO: Add these names to a std::vector to avoid using Form in the process function
            }
          }
        }
      }

      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("SingleElectronQA", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, histClasses.Data(), fConfigAddTrackHistogram); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                          // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TTrackFillMap, uint32_t TTrackMCFillMap, typename TTracks, typename TTracksMC>
  void runTrackSelection(TTracks const& tracks, TTracksMC const&, bool passEvFilter, bool writeTable)
  {
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

      if (fConfigQA && passEvFilter) {
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }

      // compute track selection and publish the bit map
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << i);
          if (fConfigQA && passEvFilter) {
            fHistMan->FillHistClass(fHistNamesReco[i].Data(), VarManager::fgValues);
          }
        }
      }
      if (writeTable)
        trackSel(static_cast<int>(filterMap));
      if (!filterMap || !passEvFilter) {
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
          if ((*sig).CheckSignal(true, track.reducedMCTrack())) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
          if ((*sig).CheckSignal(true, track.template mcParticle_as<aod::McParticles_001>())) {
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
      } // end loop over MC signals
    } // end loop over tracks
  } // end runTrackSelection

  template <typename TTracksMC>
  void runMCGenTrack(TTracksMC const& groupedMCTracks)
  {
    for (auto& mctrack : groupedMCTracks) {
      if (abs(mctrack.pdgCode()) != 11)
        continue;
      VarManager::FillTrackMC(groupedMCTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (auto& sig : fMCSignals) {
        if (sig.GetNProngs() != 1) { // NOTE: 1-prong signals required
          continue;
        }
        bool checked = false;
        if constexpr (soa::is_soa_filtered_v<TTracksMC>) {
          auto mctrack_raw = groupedMCTracks.rawIteratorAt(mctrack.globalIndex());
          checked = sig.CheckSignal(true, mctrack_raw);
        } else {
          checked = sig.CheckSignal(true, mctrack);
        }
        if (checked && fConfigQA) {
          fHistMan->FillHistClass(Form("MCTruthGenTrack_%s", sig.GetName()), VarManager::fgValues);
        }
      }
    }
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TTrackMCFillMap, typename TEvent, typename TTracks, typename TTracksMC>
  void runDataFill(TEvent const& event, TTracks const& tracks, TTracksMC const& tracksMC, bool writeTable)
  {
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::ResetValues(0, VarManager::kNMCParticleVariables);
    VarManager::FillEvent<TEventFillMap>(event);

    runTrackSelection<TTrackFillMap, TTrackMCFillMap>(tracks, tracksMC, true, writeTable);
    if (fConfigMCTruthGen)
      runMCGenTrack(tracksMC);
  }

  void processSkimmed(soa::Filtered<MyEventsSelected>::iterator const& event, MyBarrelTracks const& tracks, ReducedMCTracks const& tracksMC)
  {
    runDataFill<gkEventFillMap, gkTrackFillMap, gkParticleMCFillMap>(event, tracks, tracksMC, true);
  }
  void processAOD(MyBarrelTracksAOD const& tracks, aod::McParticles const& tracksMC)
  {
    runTrackSelection<gkTrackFillMapAOD, gkParticleMCFillMapAOD>(tracks, tracksMC, false, true);
  }
  void processAODFillHist(soa::Filtered<MyEventsSelectedAOD>::iterator const& event, MyBarrelTracksAOD const& tracks, aod::McParticles const& tracksMC)
  {
    runDataFill<gkEventFillMapAOD, gkTrackFillMapAOD, gkParticleMCFillMapAOD>(event, tracks, tracksMC, false);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }
  void processDummyAOD(MyEventsAOD&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmed, "Run barrel track selection on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processAOD, "Run barrel track selection without skimming", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processAODFillHist, "Run barrel track selection without skimming to fill track histograms", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummy, "Dummy process function", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummyAOD, "Dummy process function", false);
};

struct AnalysisSameEventPairing {
  OutputObj<THashList> fOutputList{"output"};
  Filter filterEventSelected = aod::emanalysisflags::isEventSelected == 1;
  Filter filterBarrelTrackSelected = aod::emanalysisflags::isBarrelSelected > 0;
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMCRecSignals{"cfgBarrelMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};
  Configurable<std::string> fConfigMCGenSignals{"cfgBarrelMCGenSignals", "", "Comma separated list of MC signals (generated)"};
  Configurable<std::string> fConfigAddSEPHistogram{"cfgAddSEPHistogram", "", "Comma separated list of histograms"};
  Configurable<bool> fConfigRunMCGenPair{"cfgRunMCGenPair", false, "Do pairing of true MC particles"};
  Configurable<bool> fPropToPCA{"cfgPropToPCA", false, "Propagate tracks to secondary vertex"};
  Configurable<bool> fConfigDoSecVtxProp{"cfgDoSecVtxProp", false, "Propagate tracks to secondary vertex"};
  // TODO: here we specify signals, however signal decisions are precomputed and stored in mcReducedFlags
  // TODO: The tasks based on skimmed MC could/should rely ideally just on these flags
  // TODO:   special AnalysisCuts to be prepared in this direction
  // TODO: cuts on the MC truth information to be added if needed

  HistogramManager* fHistMan;
  std::vector<std::vector<TString>> fBarrelHistNames;
  std::vector<std::vector<std::vector<TString>>> fBarrelHistNamesMCmatched;
  std::vector<MCSignal> fRecMCSignals;
  std::vector<MCSignal> fGenMCSignals;

  void init(o2::framework::InitContext& context)
  {
    bool enableBarrelHistos = context.mOptions.get<bool>("processDecayToEESkimmed") || context.mOptions.get<bool>("processDecayToEESkimmedWithCov") || context.mOptions.get<bool>("processDecayToEEAOD");

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // For barrel-barrel create:
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
        fRecMCSignals.emplace_back(*sig);
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
          fBarrelHistNames.emplace_back(names);
          std::vector<std::vector<TString>> mcSigClasses;
          if (!sigNamesStr.IsNull()) {
            for (auto& sig : fRecMCSignals) {
              std::vector<TString> names = {
                Form("PairsBarrelSEPM_%s_%s", objArray->At(icut)->GetName(), sig.GetName()),
                Form("PairsBarrelSEPP_%s_%s", objArray->At(icut)->GetName(), sig.GetName()),
                Form("PairsBarrelSEMM_%s_%s", objArray->At(icut)->GetName(), sig.GetName())};
              histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
              mcSigClasses.emplace_back(names);
            } // end loop over MC signals
          }
          fBarrelHistNamesMCmatched.emplace_back(mcSigClasses);
        } // end loop over cuts
      }   // end if(cutNames.IsNull())
    }     // end if processBarrel

    // Add histogram classes for each specified MCsignal at the generator level
    // TODO: create a std::vector of hist classes to be used at Fill time, to avoid using Form in the process function
    TString sigGenNamesStr = fConfigMCGenSignals.value;
    std::unique_ptr<TObjArray> objGenSigArray(sigGenNamesStr.Tokenize(","));
    for (int isig = 0; isig < objGenSigArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objGenSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() == 2) { // NOTE: 2-prong signals required
          fGenMCSignals.emplace_back(*sig);
          histNames += Form("MCTruthGenPair_%s;", sig->GetName());
        }
      }
    }

    DefineHistograms(fHistMan, histNames.Data(), fConfigAddSEPHistogram); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());                      // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    VarManager::SetMagneticField(5.0f);
    VarManager::SetupTwoProngDCAFitter(5.0f, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, true); // TODO: get these parameters from Configurables
    VarManager::SetupTwoProngFwdDCAFitter(5.0f, true, 200.0f, 1.0e-3f, 0.9f, true);
  }

  template <int TPairType, uint32_t TTrackFillMap, typename TTracks1, typename TTracks2>
  void runPairing(TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    // establish the right histogram classes to be filled depending on TPairType (ee,mumu,emu)
    unsigned int ncuts = fBarrelHistNames.size();
    std::vector<std::vector<TString>> histNames = fBarrelHistNames;
    std::vector<std::vector<std::vector<TString>>> histNamesMCmatched = fBarrelHistNamesMCmatched;

    // Loop over two track combinations
    uint8_t twoTrackFilter = 0;
    // uint32_t dileptonFilterMap = 0;
    // uint32_t dileptonMcDecision = 0;

    for (auto& [t1, t2] : combinations(tracks1, tracks2)) {
      if constexpr (TPairType == VarManager::kDecayToEE) {
        twoTrackFilter = uint32_t(t1.isBarrelSelected()) & uint32_t(t2.isBarrelSelected());
      }
      if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }
      VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);

      // run MC matching for this pair
      uint32_t mcDecision = 0;
      int isig = 0;
      for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) { // for skimmed DQ model
          if ((*sig).CheckSignal(true, t1.reducedMCTrack(), t2.reducedMCTrack())) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) { // for Framework data model
          if ((*sig).CheckSignal(true, t1.template mcParticle_as<aod::McParticles_001>(), t2.template mcParticle_as<aod::McParticles_001>())) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }
      } // end loop over MC signals

      // dileptonFilterMap = twoTrackFilter;
      // dileptonMcDecision = mcDecision;

      // Loop over all fulfilled cuts and fill pair histograms
      for (unsigned int icut = 0; icut < ncuts; icut++) {
        if (twoTrackFilter & (uint8_t(1) << icut)) {
          if (t1.sign() * t2.sign() < 0) {
            fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
            for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
              if (mcDecision & (uint32_t(1) << isig)) {
                fHistMan->FillHistClass(histNamesMCmatched[icut][isig][0].Data(), VarManager::fgValues);
              }
            }
          } else {
            if (t1.sign() > 0) {
              fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
                if (mcDecision & (uint32_t(1) << isig)) {
                  fHistMan->FillHistClass(histNamesMCmatched[icut][isig][1].Data(), VarManager::fgValues);
                }
              }
            } else {
              fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
                if (mcDecision & (uint32_t(1) << isig)) {
                  fHistMan->FillHistClass(histNamesMCmatched[icut][isig][2].Data(), VarManager::fgValues);
                }
              }
            }
          }
        }
      }
    } // end loop over barrel track pairs
  }   // end runPairing

  template <typename TTracksMC>
  void runMCGenPair(TTracksMC const& groupedMCTracks)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    for (auto& sig : fGenMCSignals) {
      if (sig.GetNProngs() != 2) { // NOTE: 2-prong signals required
        continue;
      }
      for (auto& [t1, t2] : combinations(groupedMCTracks, groupedMCTracks)) {
        // only select electrons
        if ((abs(t1.pdgCode()) != 11) || (abs(t2.pdgCode()) != 11))
          continue;
        // only select physical primary particles
        if (!t1.isPhysicalPrimary() || !t2.isPhysicalPrimary())
          continue;

        bool checked = false;
        if constexpr (soa::is_soa_filtered_v<TTracksMC>) {
          auto t1_raw = groupedMCTracks.rawIteratorAt(t1.globalIndex());
          auto t2_raw = groupedMCTracks.rawIteratorAt(t2.globalIndex());
          checked = sig.CheckSignal(true, t1_raw, t2_raw);
        } else {
          checked = sig.CheckSignal(true, t1, t2);
        }
        if (checked) {
          VarManager::FillPairMC<VarManager::kDecayToEE>(t1, t2);
          fHistMan->FillHistClass(Form("MCTruthGenPair_%s", sig.GetName()), VarManager::fgValues);
        }
      }
    } // end of true pairing loop
  }   // end runMCGen

  // skimmed
  PresliceUnsorted<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;
  // non skimmed
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  void processDecayToEESkimmed(soa::Filtered<MyEventsSelected>::iterator const& event,
                               soa::Filtered<MyBarrelTracksSelected> const& tracks,
                               ReducedMCEvents const&, ReducedMCTracks const& tracksMC)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event);
    VarManager::FillEvent<gkMCEventFillMap>(event.reducedMCevent());

    runPairing<VarManager::kDecayToEE, gkTrackFillMap>(tracks, tracks);
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    if (fConfigRunMCGenPair)
      runMCGenPair(groupedMCTracks);
  }

  void processDecayToEESkimmedWithCov(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event,
                                      soa::Filtered<MyBarrelTracksSelectedWithCov> const& tracks,
                                      ReducedMCEvents const&, ReducedMCTracks const& tracksMC)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCov>(event);
    VarManager::FillEvent<gkMCEventFillMap>(event.reducedMCevent());

    runPairing<VarManager::kDecayToEE, gkTrackFillMapWithCov>(tracks, tracks);
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    if (fConfigRunMCGenPair)
      runMCGenPair(groupedMCTracks);
  }

  void processDecayToEEAOD(soa::Filtered<MyEventsSelectedAOD>::iterator const& event,
                           soa::Filtered<MyBarrelTracksSelectedAOD> const& tracks,
                           aod::McCollisions const&, aod::McParticles const& tracksMC)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapAOD>(event);
    VarManager::FillEvent<gkMCEventFillMapAOD>(event.mcCollision());

    runPairing<VarManager::kDecayToEE, gkTrackFillMapAOD>(tracks, tracks);
    auto groupedMCTracks = tracksMC.sliceBy(perMcCollision, event.mcCollision().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    if (fConfigRunMCGenPair)
      runMCGenPair(groupedMCTracks);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }
  void processDummyAOD(MyEventsAOD&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmed, "Run barrel barrel pairing on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEESkimmedWithCov, "Run barrel barrel pairing on DQ skimmed covariant tracks including vertexing", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToEEAOD, "Run barrel barrel pairing on non skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy process function", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummyAOD, "Dummy process function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
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
      }
    }

    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }

    if (classStr.Contains("MCTruthGenPair")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_pair");
    }
    if (classStr.Contains("MCTruthGenTrack")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_track");
    }
  } // end loop over histogram classes
}
