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
//
// Analysis task for calculating single electron and dielectron efficiency
//
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Field/MagneticField.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "TGeoGlobalMagField.h"
#include <TH1F.h>
#include <THashList.h>
#include <TLorentzVector.h>
#include <TMath.h>
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
using MyEventsNoSkimmed = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using MyEventsSelectedNoSkimmed = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::EventCuts>;
using MyMCEventsSelectedNoSkimmed = soa::Join<aod::McCollisions, aod::EventMCCuts>;
using MyBarrelTracksNoSkimmed = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection,
                                          aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                          aod::pidTPCFullKa, aod::pidTPCFullPr,
                                          aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                          aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                          aod::McTrackLabels>;
using MyBarrelTracksSelectedNoSkimmed = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection,
                                                  aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                                  aod::pidTPCFullKa, aod::pidTPCFullPr,
                                                  aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                                  aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                                  aod::BarrelTrackCuts, aod::McTrackLabels>;
using MyMCTrackNoSkimmed = soa::Join<aod::McParticles, aod::SmearedElectrons>;

constexpr static uint32_t gkEventFillMapNoSkimmed = VarManager::ObjTypes::Collision;
constexpr static uint32_t gkMCEventFillMapNoSkimmed = VarManager::ObjTypes::CollisionMC;
constexpr static uint32_t gkTrackFillMapNoSkimmed = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID | VarManager::ObjTypes::TrackPIDExtra;
constexpr static uint32_t gkParticleMCFillMapNoSkimmed = VarManager::ObjTypes::ParticleMC;

// Skimmed data: works up to dielectron efficiency
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedMCEventLabels>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::ReducedMCEventLabels>;
using MyMCEventsSelected = soa::Join<aod::ReducedMCEvents, aod::EventMCCuts>;
using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::ReducedTracksBarrelLabels>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::ReducedTracksBarrelLabels>;
using MyMCReducedTracks = soa::Join<ReducedMCTracks, aod::SmearedElectrons>;
//
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkMCEventFillMap = VarManager::ObjTypes::ReducedEventMC;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkParticleMCFillMap = VarManager::ObjTypes::ParticleMC;

void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar); // defines histograms for all tasks
void SetBinsLinear(std::vector<double>& fBins, const double min, const double max, const unsigned int steps);

struct AnalysisEventSelection {

  Produces<aod::EventCuts> eventSel;
  OutputObj<THashList> fOutputList{"QA"};
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
  Configurable<bool> fConfigOnlyInjectedEvents{"cfgOnlyInjectedEvents", false, "Use only on Non-skimmed data! If true, select only injected events"};
  Configurable<std::vector<int>> fSubGenIDs{"cfgSubGenIDs", {0, 1, 2, 3}, "Use only on Non-skimmed data! Provide a comma separated list of subGenIDs to select, e.g. 0,1,2,3"};

  HistogramManager* fHistMan;
  AnalysisCompositeCut* fEventCut;
  HistogramRegistry registry{"HistoAnalysisEvent", {}, OutputObjHandlingPolicy::AnalysisObject};

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

      AxisSpec axisSubGen = {4, -0.5, 3.5, "MC SubGenerator ID"};
      registry.add<TH1>("Generator/SubGenerator_BeforeCuts", "", HistType::kTH1D, {axisSubGen}, true);
      registry.add<TH1>("Generator/SubGenerator_SelectedInjected", "", HistType::kTH1D, {axisSubGen}, true);
      registry.add<TH1>("Generator/SubGenerator_AfterCuts", "", HistType::kTH1D, {axisSubGen}, true);
    }
  }

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, typename TEvent, typename TEventsMC>
  void runSelection(TEvent const& event, TEventsMC const& /*mcEvents*/)
  {
    // Reset the values array
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    bool pass = true;

    int32_t subGeneratorID = -999;
    VarManager::FillEvent<TEventFillMap>(event);
    if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
      VarManager::FillEvent<TEventMCFillMap>(event.reducedMCevent());
      // TODO: Get access to subgenerator ID in skimmed data
      // generatorID = event.reducedMCevent().generatorsID();
    }
    if constexpr ((TEventMCFillMap & VarManager::ObjTypes::CollisionMC) > 0) {
      if (!event.has_mcCollision()) {
        pass = false;
      } else {
        VarManager::FillEvent<TEventMCFillMap>(event.mcCollision());
        subGeneratorID = event.mcCollision().getSubGeneratorId();
      }
    }

    registry.fill(HIST("Generator/SubGenerator_BeforeCuts"), subGeneratorID);

    // check if SubGeneratorID is part of list:
    // if SubGenerator is not part of it, reject event, return
    // fill event histos only if event is from SubGenerator
    if (fConfigOnlyInjectedEvents && !(std::find(fSubGenIDs->begin(), fSubGenIDs->end(), subGeneratorID) != fSubGenIDs->end())) {
      eventSel(0);
      return;
    }

    if (fConfigQA) {
      fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues); // automatically fill all the histograms in the class Event
      registry.fill(HIST("Generator/SubGenerator_SelectedInjected"), subGeneratorID);
    }
    if (fEventCut->IsSelected(VarManager::fgValues) && pass) {
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
        registry.fill(HIST("Generator/SubGenerator_AfterCuts"), subGeneratorID);
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

  void processNoSkimmed(MyEventsNoSkimmed::iterator const& event, aod::McCollisions const& mcEvents)
  {
    runSelection<gkEventFillMapNoSkimmed, gkMCEventFillMapNoSkimmed>(event, mcEvents);
  }

  void processDummyNoSkimmed(MyEventsNoSkimmed&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processNoSkimmed, "Run event selection without skimming", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDummy, "Dummy process function", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDummyNoSkimmed, "Dummy process function", false);
};

struct AnalysisEventQa {

  Filter filterEventSelected = aod::emanalysisflags::isEventSelected == 1;
  HistogramRegistry registry{"HistoAnalysisEventQa", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {

    VarManager::SetUseVars(AnalysisCut::fgUsedVars);
    VarManager::SetDefaultVarNames();

    registry.add("MCEvRecEv", "", HistType::kTH1D, {{20, 0., 20.}}, true);
    registry.add("mctrack", "", HistType::kTProfile, {{20, 0., 20.}}, true);
    registry.add("MCEvent", "", HistType::kTH1D, {{1, 0., 1.}}, true);
    registry.add("RecEvent", "", HistType::kTH1D, {{1, 0., 1.}}, true);
  }

  PresliceUnsorted<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackMCFillMap, typename TEvents, typename TEventsMC, typename TTracksMC>
  void runSelection(TEvents const& events, TEventsMC const& /*eventsMC*/, TTracksMC const& tracksMC)
  {

    uint8_t eventFilter = 0;
    std::map<uint64_t, int> fMCEventNbmctrack;
    std::map<uint64_t, int> fMCEventNbReco;
    std::map<uint64_t, int> fMCEventLabels;

    int fMCCounters = 0;
    // int fEvCounters = 0;

    // First loop

    for (auto& event : events) {
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEvent<TEventFillMap>(event);
      eventFilter = uint32_t(event.isEventSelected());
      if (!eventFilter)
        continue;
      // fEvCounters++;
      registry.fill(HIST("RecEvent"), 0.5);

      Int_t midrap = 0;
      Int_t globalindexmc = -1;

      // skimmed data
      if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
        auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
        midrap = dNdetach(groupedMCTracks);

        auto mcEvent = event.reducedMCevent();
        globalindexmc = mcEvent.globalIndex();
      }
      // Not skimmed data
      if constexpr ((TEventMCFillMap & VarManager::ObjTypes::CollisionMC) > 0) {
        if (!event.has_mcCollision()) {
          continue;
        }
        auto groupedMCTracks = tracksMC.sliceBy(perMcCollision, event.mcCollision().globalIndex());
        midrap = dNdetach(groupedMCTracks);

        auto mcEvent = event.mcCollision();
        globalindexmc = mcEvent.globalIndex();
      }

      if (!(fMCEventLabels.find(globalindexmc) != fMCEventLabels.end())) {
        fMCEventLabels[globalindexmc] = fMCCounters;
        fMCEventNbReco[globalindexmc] = 1;
        fMCEventNbmctrack[globalindexmc] = midrap;
        registry.fill(HIST("MCEvent"), 0.5);
        fMCCounters++;
      } else {
        fMCEventNbReco[globalindexmc] = fMCEventNbReco.find(globalindexmc)->second + 1;
      }

    } // end loop over events

    for (const auto& [mcEv, NbRecEv] : fMCEventNbReco) {
      registry.fill(HIST("MCEvRecEv"), fMCEventNbReco.find(mcEv)->second);
    }

    for (const auto& [mcEv, NbRecEv] : fMCEventNbmctrack) {
      registry.fill(HIST("mctrack"), fMCEventNbReco.find(mcEv)->second, fMCEventNbmctrack.find(mcEv)->second);
    }
  }

  template <typename TTracksMC>
  Int_t dNdetach(TTracksMC const& groupedMCTracks)
  {

    Int_t midrap = 0;
    for (auto mctrack : groupedMCTracks) {
      if (TMath::Abs(mctrack.eta()) < 0.5 && mctrack.isPhysicalPrimary() && (TMath::Abs(mctrack.pdgCode()) == 211 || mctrack.pdgCode() == 111)) {
        midrap++;
      }
    }
    return midrap;
  }

  void processSkimmed(soa::Filtered<MyEventsSelected> const& events, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runSelection<gkEventFillMap, gkMCEventFillMap, gkParticleMCFillMap>(events, eventsMC, tracksMC);
  }

  void processNoSkimmed(soa::Filtered<MyEventsSelectedNoSkimmed> const& events, aod::McCollisions const& eventsMC, aod::McParticles const& tracksMC)
  {
    runSelection<gkEventFillMapNoSkimmed, gkMCEventFillMapNoSkimmed, gkParticleMCFillMapNoSkimmed>(events, eventsMC, tracksMC);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  void processDummyNoSkimmed(MyEventsNoSkimmed&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventQa, processNoSkimmed, "Run event QA without skimming", false);
  PROCESS_SWITCH(AnalysisEventQa, processSkimmed, "Run event QA on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventQa, processDummy, "Dummy process function", false);
  PROCESS_SWITCH(AnalysisEventQa, processDummyNoSkimmed, "Dummy process function", false);
};

struct AnalysisMCEvent {

  Produces<aod::EventMCCuts> eventMCSel;

  void init(o2::framework::InitContext&)
  {

    VarManager::SetUseVars(AnalysisCut::fgUsedVars);
    VarManager::SetDefaultVarNames();
  }

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, typename TEvents, typename TEventsMC>
  void runSelection(TEvents const& events, TEventsMC const& eventsMC)
  {

    uint8_t eventFilter = 0;
    Int_t globalindex = -1;

    for (auto& eventMC : eventsMC) {
      bool pass = false;
      globalindex = eventMC.globalIndex();

      for (auto& event : events) {
        Int_t globalindexmc = -1;
        eventFilter = uint32_t(event.isEventSelected());
        if (!eventFilter)
          continue;

        // skimmed data
        if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
          auto mcEvent = event.reducedMCevent();
          globalindexmc = mcEvent.globalIndex();
        }
        // Not skimmed data
        if constexpr ((TEventMCFillMap & VarManager::ObjTypes::CollisionMC) > 0) {
          if (!event.has_mcCollision()) {
            continue;
          }
          auto mcEvent = event.mcCollision();
          globalindexmc = mcEvent.globalIndex();
        }
        if (globalindexmc == globalindex) {
          pass = true;
        }
      }
      if (pass) {
        eventMCSel(1);
      } else {
        eventMCSel(0);
      }
    }
  }

  void processSkimmed(MyEventsSelected const& events, ReducedMCEvents const& eventsMC)
  {
    runSelection<gkEventFillMap, gkMCEventFillMap>(events, eventsMC);
  }

  void processNoSkimmed(MyEventsSelectedNoSkimmed const& events, aod::McCollisions const& eventsMC)
  {
    runSelection<gkEventFillMapNoSkimmed, gkMCEventFillMapNoSkimmed>(events, eventsMC);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  void processDummyNoSkimmed(MyEventsNoSkimmed&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisMCEvent, processNoSkimmed, "Run MC event selection without skimming", false);
  PROCESS_SWITCH(AnalysisMCEvent, processSkimmed, "Run MC event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisMCEvent, processDummy, "Dummy process function", false);
  PROCESS_SWITCH(AnalysisMCEvent, processDummyNoSkimmed, "Dummy process function", false);
};

struct AnalysisTrackSelection {

  Produces<aod::BarrelTrackCuts> trackSel;
  Filter filterEventSelected = aod::emanalysisflags::isEventSelected == 1;
  Filter filterMCEventSelected = aod::emanalysisflags::isMCEventSelected == 1;

  // configurables
  Configurable<std::string> fConfigCuts{"cfgTrackCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMCSignals{"cfgTrackMCSignals", "", "Comma separated list of MC signals"};

  // 3D histos for efficiency
  Configurable<bool> fConfigUsePtVec{"cfgUsePtVecEff", true, "If true, non-linear pt bins but vector pt bins"};
  ConfigurableAxis ptBinsVec{"ptBinsVec", {0., 0., 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.20, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285, 0.29, 0.295, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80, 2.00, 2.40, 2.80, 3.20, 3.70, 4.50, 6.00, 8.00, 10., 12., 14., 16., 18., 20.}, "Pt binning vector"};
  ConfigurableAxis ptBins{"ptBins", {10, 0.f, 10.f}, "Pt binning"};
  ConfigurableAxis etaBins{"etaBins", {16, -0.8f, 0.8f}, "Eta binning"};
  ConfigurableAxis phiBins{"phiBins", {63, -0.f, 6.3f}, "Phi binning"};
  Configurable<bool> fConfigRecWithMC{"cfgEffRecWithMCVars", false, "If true, fill also 3D histograms at reconstructed level with mc variables"};
  Configurable<bool> fConfigMCCollz{"cfgMCCollz", false, "If true, look only at reconstructed track associated to mc track from a MC collision within 10cm"};

  // Resolution histos
  Configurable<bool> fConfigResolutionOn{"cfgResolution", false, "If true, fill resolution histograms"};
  Configurable<bool> fConfigUsePtVecRes{"cfgUsePtVecRes", true, "If true, non-linear pt bins predefined in res histos"};
  ConfigurableAxis ptResBinsVec{"ptResBinsVec", {0., 0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.20, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285, 0.29, 0.295, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80, 2.00, 2.40, 2.80, 3.20, 3.70, 4.50, 6.00, 8.00, 10., 12., 14., 16., 18., 20.}, "Pt binning vector for resolution"};
  ConfigurableAxis ptResBins{"ptResBins", {20, 0.f, 20.f}, "Pt binning for resolution"};
  ConfigurableAxis deltaptResBins{"deltaptResBins", {500, -0.5f, 0.5f}, "DeltaPt binning for resolution"};
  ConfigurableAxis deltaetaResBins{"deltaetaResBins", {500, -0.5f, 0.5f}, "DeltaEta binning for resolution"};
  ConfigurableAxis deltaphiResBins{"deltaphiResBins", {500, -0.5f, 0.5f}, "DeltaPhi binning for resolution"};

  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};

  // output lists
  OutputObj<THashList> fOutputQA{"SingleElectronQA"};
  HistogramRegistry registry{"HistoAnalysisTrackSelection", {}, OutputObjHandlingPolicy::AnalysisObject};
  THashList* fQASingleElectronList; // QA in case on with histo from manager outputs

  // Cuts and signals
  // AnalysisCompositeCut* fEventCut; // Taken from event selection part
  std::vector<AnalysisCompositeCut> fTrackCuts; // list of track cuts
  AnalysisCompositeCut* fTrackCutsRes;          // track cut for resolution map
  std::vector<MCSignal> fMCSignals;             // list of signals to be checked
  MCSignal* fMCSignalRes;                       // signal for res

  // 3D histos
  std::vector<std::shared_ptr<TH3>> fHistGenPosPart;
  std::vector<std::shared_ptr<TH3>> fHistGenNegPart;
  std::vector<std::shared_ptr<TH3>> fHistGenSmearedPosPart;
  std::vector<std::shared_ptr<TH3>> fHistGenSmearedNegPart;
  std::vector<std::shared_ptr<TH3>> fHistRecPosPart;
  std::vector<std::shared_ptr<TH3>> fHistRecNegPart;
  std::vector<std::shared_ptr<TH3>> fHistRecPosPartMC;
  std::vector<std::shared_ptr<TH3>> fHistRecNegPartMC;
  //
  std::vector<std::shared_ptr<TH3>> fHistRecPosSingleRecPartMC;
  std::vector<std::shared_ptr<TH3>> fHistRecNegSingleRecPartMC;
  std::vector<std::shared_ptr<TH3>> fHistRecPosClassCollDoubleCountPartMC;
  std::vector<std::shared_ptr<TH3>> fHistRecNegClassCollDoubleCountPartMC;
  std::vector<std::shared_ptr<TH3>> fHistRecPosClassAmbigCollDoubleCountPartMC;
  std::vector<std::shared_ptr<TH3>> fHistRecNegClassAmbigCollDoubleCountPartMC;

  // Res histos
  std::vector<std::shared_ptr<TH2>> fHistRes;

  // QA
  HistogramManager* fHistManQA;                            // histo manager
  std::vector<TString> fHistNamesRecoQA;                   // list of histo names for all reconstructed tracks in histo manager
  std::vector<std::vector<TString>> fHistNamesMCMatchedQA; // list of histo names for reconstructed signals in histo manager
  std::vector<TString> fHistNamesMCQA;                     // list of histo names for generated signals in histo manager

  void init(o2::framework::InitContext&)
  {

    // Create list output for QA
    fQASingleElectronList = new THashList;
    fQASingleElectronList->SetOwner(kTRUE);
    fQASingleElectronList->SetName("SEQA");

    // Binning 3D histos for single electron efficiency
    AxisSpec axisEta{etaBins, "#it{#eta}_{e}"};
    AxisSpec axisPhi{phiBins, "#it{#varphi}_{e} (rad)"};
    AxisSpec axisPt{ptBins, "#it{p}_{T,e} (GeV/#it{c})"};
    AxisSpec axisMCColl = {3, -0.5, 2.5, "MCcoll info"};
    AxisSpec axisDoubleCount = {2, -0.5, 1.5, "Double count info"};
    AxisSpec axisAmbig = {2, -0.5, 1.5, "Ambiguous info"};

    // List of track cuts
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    VarManager::SetDefaultVarNames();

    // List of MC signals
    TString configSigNamesStr = fConfigMCSignals.value;
    std::unique_ptr<TObjArray> sigNamesArray(configSigNamesStr.Tokenize(","));
    for (int isig = 0; isig < sigNamesArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(sigNamesArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() != 1) { // NOTE: only 1 prong signals
          continue;
        }
        // List of signal to be checked
        fMCSignals.push_back(*sig);
      }
    }

    // Efficiency histograms
    // Generated histograms
    // Generated true
    for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
      if (!fConfigUsePtVec) {
        fHistGenPosPart.push_back(registry.add<TH3>(Form("SingleElectron/Generated/Ngen_Pos_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
        fHistGenNegPart.push_back(registry.add<TH3>(Form("SingleElectron/Generated/Ngen_Neg_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
      } else {
        fHistGenPosPart.push_back(registry.add<TH3>(Form("SingleElectron/Generated/Ngen_Pos_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
        fHistGenNegPart.push_back(registry.add<TH3>(Form("SingleElectron/Generated/Ngen_Neg_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
      }
    }

    // Generated smeared
    for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
      if (!fConfigUsePtVec) {
        fHistGenSmearedPosPart.push_back(registry.add<TH3>(Form("SingleElectron/GeneratedSmeared/Ngen_Pos_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
        fHistGenSmearedNegPart.push_back(registry.add<TH3>(Form("SingleElectron/GeneratedSmeared/Ngen_Neg_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
      } else {
        fHistGenSmearedPosPart.push_back(registry.add<TH3>(Form("SingleElectron/GeneratedSmeared/Ngen_Pos_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
        fHistGenSmearedNegPart.push_back(registry.add<TH3>(Form("SingleElectron/GeneratedSmeared/Ngen_Neg_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
      }
    }

    // Reconstruction level: for every cutsetting one list and every MCsignal 2 histograms with pos and neg charge
    for (unsigned int list_i = 0; list_i < fTrackCuts.size(); ++list_i) {
      for (unsigned int i = 0; i < fMCSignals.size(); ++i) {

        if (!fConfigUsePtVec) {
          fHistRecPosPart.push_back(registry.add<TH3>(Form("SingleElectron/%s/Nrec_Pos_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
          fHistRecNegPart.push_back(registry.add<TH3>(Form("SingleElectron/%s/Nrec_Neg_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
        } else {
          fHistRecPosPart.push_back(registry.add<TH3>(Form("SingleElectron/%s/Nrec_Pos_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
          fHistRecNegPart.push_back(registry.add<TH3>(Form("SingleElectron/%s/Nrec_Neg_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
        }
      }
    }

    // Reconstruction level: for every cutsetting one list and every MCsignal 2 histograms with pos and neg charge, filled with mc variables
    if (fConfigRecWithMC) {
      for (unsigned int list_i = 0; list_i < fTrackCuts.size(); ++list_i) {
        for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
          if (!fConfigUsePtVec) {
            fHistRecPosPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Pos_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
            fHistRecNegPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Neg_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
            fHistRecPosSingleRecPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Pos_SingleRec_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
            fHistRecNegSingleRecPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Neg_SingleRec_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
            fHistRecPosClassCollDoubleCountPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Pos_ClassCollDoubleCount_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisMCColl, axisDoubleCount}, true));
            fHistRecNegClassCollDoubleCountPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Neg_ClassCollDoubleCount_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisMCColl, axisDoubleCount}, true));
          } else {
            fHistRecPosPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Pos_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
            fHistRecNegPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Neg_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
            fHistRecPosSingleRecPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Pos_SingleRec_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
            fHistRecNegSingleRecPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Neg_SingleRec_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
            fHistRecPosClassCollDoubleCountPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Pos_ClassCollDoubleCount_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisMCColl, axisDoubleCount}, true));
            fHistRecNegClassCollDoubleCountPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Neg_ClassCollDoubleCount_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisMCColl, axisDoubleCount}, true));
          }
          fHistRecPosClassAmbigCollDoubleCountPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Pos_ClassAmigCollDoubleCount_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisAmbig, axisMCColl, axisDoubleCount}, true));
          fHistRecNegClassAmbigCollDoubleCountPartMC.push_back(registry.add<TH3>(Form("SingleElectron/%s_MCVars/Nrec_Neg_ClassAmigCollDoubleCount_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisAmbig, axisMCColl, axisDoubleCount}, true));
        }
      }
      // Histo without track cut
      for (unsigned int i = 0; i < fMCSignals.size(); ++i) {

        if (!fConfigUsePtVec) {
          fHistRecPosPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Pos_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
          fHistRecNegPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Neg_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
          fHistRecPosSingleRecPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Pos_SingleRec_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
          fHistRecNegSingleRecPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Neg_SingleRec_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisEta, axisPhi}, true));
          fHistRecPosClassCollDoubleCountPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Pos_ClassCollDoubleCount_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisMCColl, axisDoubleCount}, true));
          fHistRecNegClassCollDoubleCountPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Neg_ClassCollDoubleCount_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisPt, axisMCColl, axisDoubleCount}, true));

        } else {
          fHistRecPosPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Pos_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
          fHistRecNegPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Neg_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
          fHistRecPosSingleRecPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Pos_SingleRec_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
          fHistRecNegSingleRecPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Neg_SingleRec_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisEta, axisPhi}, true));
          fHistRecPosClassCollDoubleCountPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Pos_ClassCollDoubleCount_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisMCColl, axisDoubleCount}, true));
          fHistRecNegClassCollDoubleCountPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Neg_ClassCollDoubleCount_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {{ptBinsVec, "#it{p}_{T,e} (GeV/#it{c})"}, axisMCColl, axisDoubleCount}, true));
        }
        fHistRecPosClassAmbigCollDoubleCountPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Pos_ClassAmigCollDoubleCount_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisAmbig, axisMCColl, axisDoubleCount}, true));
        fHistRecNegClassAmbigCollDoubleCountPartMC.push_back(registry.add<TH3>(Form("SingleElectron/NoCut_MCVars/Nrec_Neg_ClassAmigCollDoubleCount_%s", fMCSignals.at(i).GetName()), "", HistType::kTH3D, {axisAmbig, axisMCColl, axisDoubleCount}, true));
      }
    }

    // Resolution histogramms
    if (fConfigResolutionOn) {
      // Binning for resolution
      AxisSpec axisPtRes{ptResBins, "#it{p}^{gen}_{T,e} (GeV/#it{c})"};
      AxisSpec axisDeltaptRes{deltaptResBins, "(p^{gen}_{T} - p^{rec}_{T}) / p^{gen}_{T} (GeV/c)"};
      AxisSpec axisDeltaetaRes{deltaetaResBins, "#eta^{gen} - #eta^{rec}"};
      AxisSpec axisDeltaphiRes{deltaphiResBins, "#varphi^{gen} - #varphi^{rec} (rad)"};

      // Create the histos
      if (!fConfigUsePtVecRes) {
        fHistRes.push_back(registry.add<TH2>("Resolution/PtGen_DeltaPtOverPtGen", "", HistType::kTH2D, {axisPtRes, axisDeltaptRes}, true));
        fHistRes.push_back(registry.add<TH2>("Resolution/PtGen_DeltaEta", "", HistType::kTH2D, {axisPtRes, axisDeltaetaRes}, true));
        fHistRes.push_back(registry.add<TH2>("Resolution/PtGen_DeltaPhi_Ele", "", HistType::kTH2D, {axisPtRes, axisDeltaphiRes}, true));
        fHistRes.push_back(registry.add<TH2>("Resolution/PtGen_DeltaPhi_Pos", "", HistType::kTH2D, {axisPtRes, axisDeltaphiRes}, true));
      } else {
        fHistRes.push_back(registry.add<TH2>("Resolution/PtGen_DeltaPtOverPtGen", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,e} (GeV/#it{c})"}, axisDeltaptRes}, true));
        fHistRes.push_back(registry.add<TH2>("Resolution/PtGen_DeltaEta", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,e} (GeV/#it{c})"}, axisDeltaetaRes}, true));
        fHistRes.push_back(registry.add<TH2>("Resolution/PtGen_DeltaPhi_Ele", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,e} (GeV/#it{c})"}, axisDeltaphiRes}, true));
        fHistRes.push_back(registry.add<TH2>("Resolution/PtGen_DeltaPhi_Pos", "", HistType::kTH2D, {{ptResBinsVec, "#it{p}^{gen}_{T,e} (GeV/#it{c})"}, axisDeltaphiRes}, true));
      }
    }

    // Configure QA histogram classes
    if (fConfigQA) {
      TString histClassesQA = "TrackBarrel_BeforeCuts;";
      for (auto& cut : fTrackCuts) {

        // All reconstructed leptons
        TString nameStr = Form("TrackBarrel_%s", cut.GetName());
        fHistNamesRecoQA.push_back(nameStr);
        histClassesQA += Form("%s;", nameStr.Data());

        // All reconstructed leptons matched to a 1 prong signal or MC 1 prong signal directly
        std::vector<TString> mcnamesreco;
        for (unsigned int isig = 0; isig < fMCSignals.size(); ++isig) {
          TString nameStr2 = Form("TrackBarrel_%s_%s", cut.GetName(), fMCSignals.at(isig).GetName());
          mcnamesreco.push_back(nameStr2);
          histClassesQA += Form("%s;", nameStr2.Data());
        }
        fHistNamesMCMatchedQA.push_back(mcnamesreco);
      }

      // Add histogram classes for each MC signal at generated level
      // std::vector<TString> mcnamesgen;
      for (unsigned int isig = 0; isig < fMCSignals.size(); ++isig) {
        TString nameStr2 = Form("MCTruthGen_%s", fMCSignals.at(isig).GetName());
        fHistNamesMCQA.push_back(nameStr2);
        // mcnamesgen.push_back(nameStr2);
        histClassesQA += Form("%s;", nameStr2.Data());
      }
      // fHistNamesMCQA.push_back(mcnamesgen);

      fHistManQA = new HistogramManager("SingleElectronQA", "aa", VarManager::kNVars);
      fHistManQA->SetUseDefaultVariableNames(kTRUE);
      fHistManQA->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistManQA, histClassesQA.Data(), fConfigAddTrackHistogram); // define all histograms
      VarManager::SetUseVars(fHistManQA->GetUsedVars());                            // provide the list of required variables so that VarManager knows what to fill
      fQASingleElectronList = fHistManQA->GetMainHistogramList();
    }
    fOutputQA.setObject(fQASingleElectronList);
  }

  PresliceUnsorted<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;
  Preslice<MyBarrelTracks> perReducedEventTracks = aod::reducedtrack::reducedeventId;

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<MyBarrelTracksNoSkimmed> perCollisionTracks = aod::track::collisionId;

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackFillMap, uint32_t TTrackMCFillMap, typename TEvents, typename TTracks, typename TEventsMC, typename TTracksMC>
  void runSelection(TEvents const& events, TTracks const& tracks, TEventsMC const& /*eventsMC*/, TTracksMC const& tracksMC, bool write)
  {

    uint8_t eventFilter = 0;
    bool pass = true;
    std::map<uint64_t, int> fMCEventLabels;
    int fCounters = 0; //! [0] - particle counter, [1] - event counter

    for (auto& event : events) {
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::ResetValues(0, VarManager::kNMCParticleVariables);
      // fill event information which might be needed in histograms that combine track and event properties
      VarManager::FillEvent<TEventFillMap>(event);
      // if(!fEventCut->IsSelected(VarManager::fgValues)) continue;
      eventFilter = uint32_t(event.isEventSelected());
      if (!eventFilter) {
        pass = false;
      }
      if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
        VarManager::FillEvent<TEventMCFillMap>(event.reducedMCevent());
      }
      if constexpr ((TEventMCFillMap & VarManager::ObjTypes::CollisionMC) > 0) {
        VarManager::FillEvent<TEventMCFillMap>(event.mcCollision());
      }

      // Look if we did not already saw the collision and fill the denominator of the single electron efficiency
      Int_t globalindexmc = -1;
      if (pass) {
        if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
          auto mcEvent = event.reducedMCevent();
          globalindexmc = mcEvent.globalIndex();
        }
        if constexpr ((TEventMCFillMap & VarManager::ObjTypes::CollisionMC) > 0) {
          auto mcEvent = event.mcCollision();
          globalindexmc = mcEvent.globalIndex();
        }
        if (!(fMCEventLabels.find(globalindexmc) != fMCEventLabels.end())) {
          fMCEventLabels[globalindexmc] = fCounters;
          fCounters++;
          // skimmed data
          if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
            auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
            groupedMCTracks.bindInternalIndicesTo(&tracksMC);
            runMCGenTrack<false>(groupedMCTracks);
          }
          // Not skimmed data
          if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
            auto groupedMCTracks = tracksMC.sliceBy(perMcCollision, event.mcCollision().globalIndex());
            groupedMCTracks.bindInternalIndicesTo(&tracksMC);
            runMCGenTrack<false>(groupedMCTracks);
          }
        }
      }

      // Loop over reconstructed tracks belonging to the event and fill the numerator of the efficiency as well as the resolution map
      if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
        auto groupedTracks = tracks.sliceBy(perReducedEventTracks, event.globalIndex());
        runRecTrack<gkTrackFillMap>(groupedTracks, tracksMC, pass, write);
      }
      if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
        auto groupedTracks = tracks.sliceBy(perCollisionTracks, event.globalIndex());
        runRecTrack<gkTrackFillMapNoSkimmed>(groupedTracks, tracksMC, pass, write);
      }
    } // end loop over events
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TTrackMCFillMap, typename TEvent, typename TTracks, typename TTracksMC>
  void runDataFill(TEvent const& event, TTracks const& tracks, TTracksMC const& tracksMC, bool write)
  {

    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::ResetValues(0, VarManager::kNMCParticleVariables);
    VarManager::FillEvent<TEventFillMap>(event);

    runRecTrack<TTrackFillMap>(tracks, tracksMC, true, write);
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TTrackMCFillMap, typename TEvents, typename TEventsMC, typename TTracks, typename TTracksMC, typename TAmbigTracks>
  void runDataFillMore(TEvents const& events, TEventsMC eventsMC, TTracks const& tracks, TTracksMC const& tracksMC, TAmbigTracks const& ambiTracksMid)
  {

    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::ResetValues(0, VarManager::kNMCParticleVariables);

    runRecTrackMore<TEventFillMap, TTrackFillMap>(events, eventsMC, tracks, tracksMC, ambiTracksMid);
  }

  template <uint32_t TEventMCFillMap, uint32_t TTrackMCFillMap, typename TEventsMC, typename TTracksMC>
  void runMCFill(TEventsMC const& eventMC, TTracksMC const& tracksMC)
  {
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::FillEvent<TEventMCFillMap>(eventMC);
    runMCGenTrack<true>(tracksMC);
  }

  template <uint32_t TEventMCFillMap, uint32_t TTrackMCFillMap, typename TEventsMC, typename TTracksMC>
  void runMCFillMore(TEventsMC const& eventMC, TTracksMC const& tracksMC)
  {
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::FillEvent<TEventMCFillMap>(eventMC);
    runMCGenTrackMore<true>(tracksMC, eventMC);
  }

  template <uint32_t TTrackFillMap, uint32_t TTrackMCFillMap, typename TTracks, typename TTracksMC>
  void runDataSelection(TTracks const& tracks, TTracksMC const& tracksMC)
  {

    runRecTrack<TTrackFillMap>(tracks, tracksMC, false, true);
  }

  template <bool smeared, typename TTracksMC>
  void runMCGenTrack(TTracksMC const& groupedMCTracks)
  {
    for (auto& mctrack : groupedMCTracks) {
      VarManager::ResetValues(0, VarManager::kNMCParticleVariables);
      VarManager::FillTrackMC(groupedMCTracks, mctrack);
      int isig = 0;
      for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {
        bool checked = false;
        if constexpr (soa::is_soa_filtered_v<TTracksMC>) {
          auto mctrack_raw = groupedMCTracks.rawIteratorAt(mctrack.globalIndex());
          checked = (*sig).CheckSignal(true, mctrack_raw);
        } else {
          checked = (*sig).CheckSignal(true, mctrack);
        }
        if (checked) {
          if (mctrack.pdgCode() > 0) {
            fHistGenNegPart[isig]->Fill(mctrack.pt(), mctrack.eta(), mctrack.phi());
            if constexpr (smeared)
              fHistGenSmearedNegPart[isig]->Fill(mctrack.ptSmeared(), mctrack.etaSmeared(), mctrack.phiSmeared());
          } else {
            fHistGenPosPart[isig]->Fill(mctrack.pt(), mctrack.eta(), mctrack.phi());
            if constexpr (smeared)
              fHistGenSmearedPosPart[isig]->Fill(mctrack.ptSmeared(), mctrack.etaSmeared(), mctrack.phiSmeared());
          }
          if (fConfigQA) {
            fHistManQA->FillHistClass(fHistNamesMCQA[isig].Data(), VarManager::fgValues);
          }
        }
      }
    }
  }

  template <bool smeared, typename TTracksMC, typename TEventsMC>
  void runMCGenTrackMore(TTracksMC const& groupedMCTracks, TEventsMC const& /*eventMC*/)
  {

    for (auto& mctrack : groupedMCTracks) {

      bool mccollisionwithin10 = false;
      auto mccollision = mctrack.mcCollision();
      Double_t zmc = mccollision.posZ();
      if (TMath::Abs(zmc) < 10.)
        mccollisionwithin10 = true;

      if (!mccollisionwithin10 && fConfigMCCollz)
        continue;

      VarManager::ResetValues(0, VarManager::kNMCParticleVariables);
      VarManager::FillTrackMC(groupedMCTracks, mctrack);
      int isig = 0;
      for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {
        bool checked = false;
        if constexpr (soa::is_soa_filtered_v<TTracksMC>) {
          auto mctrack_raw = groupedMCTracks.rawIteratorAt(mctrack.globalIndex());
          checked = (*sig).CheckSignal(true, mctrack_raw);
        } else {
          checked = (*sig).CheckSignal(true, mctrack);
        }
        if (checked) {
          if (mctrack.pdgCode() > 0) {
            fHistGenNegPart[isig]->Fill(mctrack.pt(), mctrack.eta(), mctrack.phi());
            if constexpr (smeared)
              fHistGenSmearedNegPart[isig]->Fill(mctrack.ptSmeared(), mctrack.etaSmeared(), mctrack.phiSmeared());
          } else {
            fHistGenPosPart[isig]->Fill(mctrack.pt(), mctrack.eta(), mctrack.phi());
            if constexpr (smeared)
              fHistGenSmearedPosPart[isig]->Fill(mctrack.ptSmeared(), mctrack.etaSmeared(), mctrack.phiSmeared());
          }
          if (fConfigQA)
            // fHistManQA->FillHistClass(Form("MCTruthGen_%s", (*sig).GetName()), VarManager::fgValues);
            fHistManQA->FillHistClass(fHistNamesMCQA[isig].Data(), VarManager::fgValues);
        }
      }
    }
  }

  template <uint32_t TTrackFillMap, typename TTracks, typename TTracksMC>
  void runRecTrack(TTracks const& groupedTracks, TTracksMC const& tracksMC, bool pass, bool write)
  {

    uint32_t filterMap = 0;
    trackSel.reserve(groupedTracks.size());

    for (auto& track : groupedTracks) {
      filterMap = 0;

      VarManager::ResetValues(0, VarManager::kNMCParticleVariables);
      VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
      VarManager::FillTrack<TTrackFillMap>(track); // compute track quantities

      // compute MC matched quantities
      if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
        VarManager::FillTrackMC(tracksMC, track.reducedMCTrack());
      }
      if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
        // If no MC particle is found, skip the track
        if (track.has_mcParticle()) {
          auto mctrack = track.template mcParticle_as<aod::McParticles>();
          VarManager::FillTrackMC(tracksMC, mctrack);
        }
      }

      // no track cut
      if (fConfigQA && pass) {
        fHistManQA->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }

      // compute track selection and publish the bit map
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << i);
          if (fConfigQA && pass) {
            fHistManQA->FillHistClass(fHistNamesRecoQA[i].Data(), VarManager::fgValues);
          }
        }
      }
      if (write)
        trackSel(static_cast<int>(filterMap));
      if (!filterMap || !pass) {
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
          if (track.has_mcParticle()) {
            auto mctrack = track.template mcParticle_as<aod::McParticles>();
            if ((*sig).CheckSignal(true, mctrack)) {
              mcDecision |= (uint32_t(1) << isig);
            }
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
            if (track.sign() < 0) {
              fHistRecNegPart[j * fMCSignals.size() + i]->Fill(track.pt(), track.eta(), track.phi());
            } else {
              fHistRecPosPart[j * fMCSignals.size() + i]->Fill(track.pt(), track.eta(), track.phi());
            }

            if (fConfigRecWithMC) {

              Double_t mcpt = -10000.;
              Double_t mceta = -10000.;
              Double_t mcphi = -1000.;

              if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
                auto mctrack = track.reducedMCTrack();
                mcpt = mctrack.pt();
                mceta = mctrack.eta();
                mcphi = mctrack.phi();
              }
              if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
                if (track.has_mcParticle()) {
                  auto mctrack = track.template mcParticle_as<aod::McParticles>();
                  mcpt = mctrack.pt();
                  mceta = mctrack.eta();
                  mcphi = mctrack.phi();
                }
              }

              if (track.sign() < 0) {
                fHistRecNegPartMC[j * fMCSignals.size() + i]->Fill(mcpt, mceta, mcphi);
              } else {
                fHistRecPosPartMC[j * fMCSignals.size() + i]->Fill(mcpt, mceta, mcphi);
              }
            }

            if (fConfigResolutionOn && (i == 0) && (j == 0)) {

              Double_t mcpt = -10000.;
              Double_t mceta = -10000.;
              Double_t mcphi = -1000.;
              Int_t mcpdg = -10000.;

              if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
                auto mctrack = track.reducedMCTrack();
                mcpt = mctrack.pt();
                mceta = mctrack.eta();
                mcphi = mctrack.phi();
                mcpdg = mctrack.pdgCode();
              }
              if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
                if (track.has_mcParticle()) {
                  auto mctrack = track.template mcParticle_as<aod::McParticles>();
                  mcpt = mctrack.pt();
                  mceta = mctrack.eta();
                  mcphi = mctrack.phi();
                  mcpdg = mctrack.pdgCode();
                }
              }
              Double_t deltaptoverpt = -1000.;
              if (mcpt > 0.)
                deltaptoverpt = (mcpt - track.pt()) / mcpt;
              Double_t deltaeta = mceta - track.eta();
              Double_t deltaphi = mcphi - track.phi();
              fHistRes[0]->Fill(mcpt, deltaptoverpt);
              fHistRes[1]->Fill(mcpt, deltaeta);
              if (mcpdg < 0) {
                fHistRes[2]->Fill(mcpt, deltaphi);
              } else {
                fHistRes[3]->Fill(mcpt, deltaphi);
              }
            }
            if (fConfigQA)
              fHistManQA->FillHistClass(fHistNamesMCMatchedQA[j][i].Data(), VarManager::fgValues);
          }
        } // end loop over cuts
      } // end loop over MC signals
    } // end loop over reconstructed track belonging to the events
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TEventsMC, typename TTracks, typename TTracksMC, typename TAmbigTracks>
  void runRecTrackMore(TEvents const& events, TEventsMC const& /*eventsMC*/, TTracks const& groupedTracks, TTracksMC const& tracksMC, TAmbigTracks const& ambiTracksMid)
  {

    std::map<uint64_t, int> fRecTrackLabels[fTrackCuts.size() + 1];

    uint32_t filterMap = 0;
    trackSel.reserve(groupedTracks.size());

    for (auto& track : groupedTracks) {

      // How many time the associated MC track was seen for this cut
      Int_t fRecCounters[fTrackCuts.size() + 1];
      for (unsigned int k = 0; k < fTrackCuts.size() + 1; k++) {
        fRecCounters[k] = 0;
      }

      filterMap = 0;
      Int_t ambiguousinfo = 0;
      Int_t collisioninfo = -1;
      Int_t mcCollisionIddmctrack = -999;
      Int_t mcCollisionIdrectrack = -999;
      bool mccollisionwithin10 = false;

      VarManager::FillTrack<TTrackFillMap>(track); // compute track quantities

      // Do ambiguous tracks
      for (auto& ambiTrackMid : ambiTracksMid) {
        if (ambiTrackMid.trackId() == track.globalIndex()) {
          ambiguousinfo = 1;
          break;
        }
      }

      // Do matching
      // if (!track.has_collision()) printf("CollisionId %d\n",track.collisionId());
      if (track.has_collision()) {
        Int_t reccollisionid = track.collisionId();
        if (ambiguousinfo == 1)
          printf("Has reccollision but is ambiguous\n");
        // printf("Look for the reconstructed collision %d\n",reccollisionid);
        bool pass = 0;
        for (auto& event : events) {
          if (event.isEventSelected() == 1) {
            VarManager::FillEvent<TEventFillMap>(event);
            // printf("Global index of collision %d\n",event.globalIndex());
            if (reccollisionid == event.globalIndex()) {
              pass = 1;
              // printf("Found a collision with the same id %d and %d\n",reccollisionid,event.globalIndex());
              if (ambiguousinfo == 1)
                printf("Has reccollision and found it in the list but is ambiguous\n");
              if (event.has_mcCollision()) {
                mcCollisionIdrectrack = event.mcCollisionId();
                if (ambiguousinfo == 1)
                  printf("Has reccollision with mccollision but is ambiguous\n");
              } else {
                if (ambiguousinfo == 1)
                  printf("Has reccollision but without mccollision and is ambiguous\n");
              }
              break;
            }
          }
        }
        if (!pass) // rec collision of track is not selected by isSelected
          continue;
        // else       rec collision of track is selected by isSelected
      } else {
        // printf("Not attached to a reconstructed collision\n");
      }

      if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
        // If no MC particle is found, skip the track
        if (track.has_mcParticle()) {
          // if (ambiguousinfo == 1) printf("Has mcparticle but is ambiguous\n");
          // printf("Found a mc track\n");
          auto mctrack = track.template mcParticle_as<aod::McParticles>();
          mcCollisionIddmctrack = mctrack.mcCollisionId();
          auto mccollision = mctrack.mcCollision();
          Double_t zmc = mccollision.posZ();
          if (TMath::Abs(zmc) < 10.)
            mccollisionwithin10 = true;
          VarManager::FillTrackMC(tracksMC, mctrack);
        }
      }

      // no track cut
      if (fConfigQA) {
        fHistManQA->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }

      // if (ambiguousinfo == 1) printf("Values are %d and %d\n",mcCollisionIddmctrack,mcCollisionIdrectrack);
      //  compute collision and mccollision info
      if (mcCollisionIddmctrack != -999 && mcCollisionIdrectrack != -999) {
        if (mcCollisionIddmctrack == mcCollisionIdrectrack)
          collisioninfo = 0;
        else
          collisioninfo = 1;
      } else {
        collisioninfo = 2;
      }
      // printf("collision info %d\n",collisioninfo);

      // compute track selection and publish the bit map
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << i);
          if (fConfigQA) {
            fHistManQA->FillHistClass(fHistNamesRecoQA[i].Data(), VarManager::fgValues);
          }
        }
      }

      // compute MC matching decisions
      uint32_t mcDecision = 0;
      int isig = 0;
      Int_t mctrackindex = -999;
      Int_t doublereconstructedtrack[fTrackCuts.size() + 1];
      for (unsigned int k = 0; k < fTrackCuts.size() + 1; k++) {
        doublereconstructedtrack[k] = 0;
      }
      for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
          if (track.has_mcParticle()) {
            auto mctrack = track.template mcParticle_as<aod::McParticles>();
            if ((*sig).CheckSignal(true, mctrack)) {
              mcDecision |= (uint32_t(1) << isig);
              mctrackindex = mctrack.globalIndex();
            }
          }
        }
      }

      // // Double reconstructed track only for the signal (they should not be redundant or crossing!!)
      // for (unsigned int i = 0; i < fMCSignals.size(); i++) {
      //   if (!(mcDecision & (uint32_t(1) << i))) {
      //     continue;
      //   }

      // no track cuts
      if (!(fRecTrackLabels[fTrackCuts.size()].find(mctrackindex) != fRecTrackLabels[fTrackCuts.size()].end())) {
        fRecTrackLabels[fTrackCuts.size()][mctrackindex] = fRecCounters[fTrackCuts.size()];
        fRecCounters[fTrackCuts.size()]++;
      } else {
        // printf("For cut %d, found a mc collision track already reconstructed %d for selected collision with the same mc collision %d\n",j,mctrackindex,mcCollisionId);
        doublereconstructedtrack[fTrackCuts.size()] = 1;
        fRecTrackLabels[fTrackCuts.size()][mctrackindex] = fRecTrackLabels[fTrackCuts.size()].find(mctrackindex)->second + 1;
      }
      // track cuts
      for (unsigned int j = 0; j < fTrackCuts.size(); j++) {
        if (filterMap & (uint8_t(1) << j)) {

          if (!(fRecTrackLabels[j].find(mctrackindex) != fRecTrackLabels[j].end())) {
            fRecTrackLabels[j][mctrackindex] = fRecCounters[j];
            fRecCounters[j]++;
          } else {
            // printf("For cut %d, found a mc collision track already reconstructed %d for selected collision with the same mc collision %d\n",j,mctrackindex,mcCollisionId);
            doublereconstructedtrack[j] = 1;
            fRecTrackLabels[j][mctrackindex] = fRecTrackLabels[j].find(mctrackindex)->second + 1;
          }
        }
      }
      // }

      // fill histograms
      for (unsigned int i = 0; i < fMCSignals.size(); i++) {
        if (!(mcDecision & (uint32_t(1) << i))) {
          continue;
        }
        Double_t mmcpt = -10000.;
        Double_t mmceta = -10000.;
        Double_t mmcphi = -1000.;
        // No track cut
        if (fConfigRecWithMC) {
          if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
            auto mctrack = track.reducedMCTrack();
            mmcpt = mctrack.pt();
            mmceta = mctrack.eta();
            mmcphi = mctrack.phi();
          }
          if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
            if (track.has_mcParticle()) {
              auto mctrack = track.template mcParticle_as<aod::McParticles>();
              mmcpt = mctrack.pt();
              mmceta = mctrack.eta();
              mmcphi = mctrack.phi();
            }
          }

          if ((mccollisionwithin10 && fConfigMCCollz) || (!fConfigMCCollz)) {
            if (track.sign() < 0) {
              fHistRecNegPartMC[fTrackCuts.size() * fMCSignals.size() + i]->Fill(mmcpt, mmceta, mmcphi);
              if (doublereconstructedtrack[fTrackCuts.size()] == 0)
                fHistRecNegSingleRecPartMC[fTrackCuts.size() * fMCSignals.size() + i]->Fill(mmcpt, mmceta, mmcphi);
              if (TMath::Abs(mmceta) < 0.8) {
                fHistRecNegClassCollDoubleCountPartMC[fTrackCuts.size() * fMCSignals.size() + i]->Fill(mmcpt, collisioninfo, doublereconstructedtrack[fTrackCuts.size()]);
                fHistRecNegClassAmbigCollDoubleCountPartMC[fTrackCuts.size() * fMCSignals.size() + i]->Fill(ambiguousinfo, collisioninfo, doublereconstructedtrack[fTrackCuts.size()]);
              }
            } else {
              fHistRecPosPartMC[fTrackCuts.size() * fMCSignals.size() + i]->Fill(mmcpt, mmceta, mmcphi);
              if (doublereconstructedtrack[fTrackCuts.size()] == 0)
                fHistRecPosSingleRecPartMC[fTrackCuts.size() * fMCSignals.size() + i]->Fill(mmcpt, mmceta, mmcphi);
              if (TMath::Abs(mmceta) < 0.8) {
                fHistRecPosClassCollDoubleCountPartMC[fTrackCuts.size() * fMCSignals.size() + i]->Fill(mmcpt, collisioninfo, doublereconstructedtrack[fTrackCuts.size()]);
                fHistRecPosClassAmbigCollDoubleCountPartMC[fTrackCuts.size() * fMCSignals.size() + i]->Fill(ambiguousinfo, collisioninfo, doublereconstructedtrack[fTrackCuts.size()]);
              }
            }
          }
        }
        // track cuts
        for (unsigned int j = 0; j < fTrackCuts.size(); j++) {
          if (filterMap & (uint8_t(1) << j)) {
            if (track.sign() < 0) {
              fHistRecNegPart[j * fMCSignals.size() + i]->Fill(track.pt(), track.eta(), track.phi());
            } else {
              fHistRecPosPart[j * fMCSignals.size() + i]->Fill(track.pt(), track.eta(), track.phi());
            }

            if (fConfigRecWithMC) {
              if ((mccollisionwithin10 && fConfigMCCollz) || (!fConfigMCCollz)) {
                if (track.sign() < 0) {
                  fHistRecNegPartMC[j * fMCSignals.size() + i]->Fill(mmcpt, mmceta, mmcphi);
                  if (doublereconstructedtrack[j] == 0)
                    fHistRecNegSingleRecPartMC[j * fMCSignals.size() + i]->Fill(mmcpt, mmceta, mmcphi);
                  if (TMath::Abs(mmceta) < 0.8) {
                    fHistRecNegClassCollDoubleCountPartMC[j * fMCSignals.size() + i]->Fill(mmcpt, collisioninfo, doublereconstructedtrack[j]);
                    fHistRecNegClassAmbigCollDoubleCountPartMC[j * fMCSignals.size() + i]->Fill(ambiguousinfo, collisioninfo, doublereconstructedtrack[fTrackCuts.size()]);
                  }
                } else {
                  fHistRecPosPartMC[j * fMCSignals.size() + i]->Fill(mmcpt, mmceta, mmcphi);
                  if (doublereconstructedtrack[j] == 0)
                    fHistRecPosSingleRecPartMC[j * fMCSignals.size() + i]->Fill(mmcpt, mmceta, mmcphi);
                  if (TMath::Abs(mmceta) < 0.8) {
                    fHistRecPosClassCollDoubleCountPartMC[j * fMCSignals.size() + i]->Fill(mmcpt, collisioninfo, doublereconstructedtrack[j]);
                    fHistRecPosClassAmbigCollDoubleCountPartMC[j * fMCSignals.size() + i]->Fill(ambiguousinfo, collisioninfo, doublereconstructedtrack[fTrackCuts.size()]);
                  }
                }
              }
            }

            if (fConfigQA)
              fHistManQA->FillHistClass(fHistNamesMCMatchedQA[j][i].Data(), VarManager::fgValues);
          }
        } // end loop over cuts
      } // end loop over MC signals
    } // end loop over reconstructed track belonging to the events
  }

  void processSkimmed(soa::Filtered<MyEventsSelected> const& events, MyBarrelTracks const& tracks, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runSelection<gkEventFillMap, gkMCEventFillMap, gkTrackFillMap, gkParticleMCFillMap>(events, tracks, eventsMC, tracksMC, true);
  }

  void processDataSkimmed(soa::Filtered<MyEventsSelected>::iterator const& event, MyBarrelTracks const& tracks, ReducedMCTracks const& tracksMC)
  {
    runDataFill<gkEventFillMap, gkTrackFillMap, gkParticleMCFillMap>(event, tracks, tracksMC, true);
  }

  void processMCSkimmed(soa::Filtered<MyMCEventsSelected> const& eventsMC, MyMCReducedTracks const& tracksMC)
  {
    for (auto& eventMC : eventsMC) {
      auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, eventMC.globalIndex());
      groupedMCTracks.bindInternalIndicesTo(&tracksMC);
      runMCFill<gkMCEventFillMap, gkParticleMCFillMap>(eventMC, groupedMCTracks);
    }
  }

  void processDataSelectionNoSkimmed(MyBarrelTracksNoSkimmed const& tracks, aod::McParticles const& tracksMC)
  {
    runDataSelection<gkTrackFillMapNoSkimmed, gkParticleMCFillMapNoSkimmed>(tracks, tracksMC);
  }

  void processMCNoSkimmed(soa::Filtered<MyMCEventsSelectedNoSkimmed>::iterator const& eventMC, MyMCTrackNoSkimmed const& tracksMC)
  // void processMCNoSkimmed(aod::McCollisions::iterator const& eventMC, MyMCTrackNoSkimmed const& tracksMC)
  {
    runMCFill<gkMCEventFillMapNoSkimmed, gkParticleMCFillMapNoSkimmed>(eventMC, tracksMC);
  }

  // void processMCNoSkimmed(MyMCTrackNoSkimmed const& tracksMC)
  // {
  //   runMCGenTrack<true>(tracksMC);
  // }

  void processMCNoSkimmedMore(soa::Filtered<MyMCEventsSelectedNoSkimmed>::iterator const& eventMC, MyMCTrackNoSkimmed const& tracksMC)
  // void processMCNoSkimmedMore(aod::McCollisions::iterator const& eventMC, MyMCTrackNoSkimmed const& tracksMC)
  {
    runMCFillMore<gkMCEventFillMapNoSkimmed, gkParticleMCFillMapNoSkimmed>(eventMC, tracksMC);
  }

  // void processMCNoSkimmedMore(MyMCTrackNoSkimmed const& tracksMC, aod::McCollisions const& eventsMC)
  // {
  //   runMCGenTrackMore<true>(tracksMC, eventsMC);
  // }

  void processDataNoSkimmed(soa::Filtered<MyEventsSelectedNoSkimmed>::iterator const& event, MyBarrelTracksNoSkimmed const& tracks, aod::McParticles const& tracksMC)
  {
    runDataFill<gkEventFillMapNoSkimmed, gkTrackFillMapNoSkimmed, gkParticleMCFillMapNoSkimmed>(event, tracks, tracksMC, false);
  }

  // void processDataNoSkimmed(MyBarrelTracksNoSkimmed const& tracks, aod::McParticles const& tracksMC)
  // {
  //   runRecTrack<gkTrackFillMapNoSkimmed>(tracks, tracksMC, true, false);
  //   // runDataFill<gkEventFillMapNoSkimmed, gkTrackFillMapNoSkimmed, gkParticleMCFillMapNoSkimmed>(event, tracks, tracksMC, false);
  // }

  void processDataNoSkimmedMore(soa::Filtered<MyEventsSelectedNoSkimmed> const& events, aod::McCollisions const& eventsMC, MyBarrelTracksNoSkimmed const& tracks, aod::McParticles const& tracksMC, aod::AmbiguousTracks const& ambiTracksMid)
  {
    runDataFillMore<gkEventFillMapNoSkimmed, gkTrackFillMapNoSkimmed, gkParticleMCFillMapNoSkimmed>(events, eventsMC, tracks, tracksMC, ambiTracksMid);
  }

  // void processDataNoSkimmedMore(MyBarrelTracksNoSkimmed const& tracks, aod::McParticles const& tracksMC, aod::AmbiguousTracks const& ambiTracksMid, MyEventsNoSkimmed const& events, aod::McCollisions const& eventsMC)
  // {
  //   runRecTrackMore<gkTrackFillMapNoSkimmed>(events, eventsMC, tracks, tracksMC, ambiTracksMid);
  // }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  void processDummyNoSkimmed(MyEventsNoSkimmed&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisTrackSelection, processDataSelectionNoSkimmed, "Run barrel track selection without skimming", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processMCNoSkimmed, "Run the barrel MC track filling without skimming", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processMCNoSkimmedMore, "Run the barrel MC track filling without skimming", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDataNoSkimmed, "Run the barrel data track filling without skimming", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDataNoSkimmedMore, "Run the barrel data track filling without skimming", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmed, "Run the data and mc barrel data track selection and filling on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDataSkimmed, "Run the barrel data track selection and filling on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processMCSkimmed, "Run the barrel mc track filling on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummy, "Dummy process function", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummyNoSkimmed, "Dummy process function", false);
};

struct AnalysisSameEventPairing {

  // Filter based on previous components in the task
  Filter filterEventSelected = aod::emanalysisflags::isEventSelected == 1;
  Filter filterMCEventSelected = aod::emanalysisflags::isMCEventSelected == 1;
  Filter filterBarrelTrackSelected = aod::emanalysisflags::isBarrelSelected > 0;
  // Filter filterMCTrackSelected = nabs(o2::aod::mcparticle::pdgCode) == 11;

  // Partition
  // Partition<ReducedMCTracks> mcSkimmedElectrons = nabs(o2::aod::mcparticle::pdgCode) == 11;
  // Partition<aod::McParticles> mcNotSkimmedElectrons = nabs(o2::aod::mcparticle::pdgCode) == 11;

  float mMagField = 0.0;
  o2::parameters::GRPMagField* grpmag = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMCSignals{"cfgBarrelMCSignals", "", "Comma separated list of MC signals"};
  Configurable<double> fConfigMinPt{"cfgMinPtFidCut", 0., "Fiducial min Pt for MC signal"};
  Configurable<double> fConfigMaxPt{"cfgMaxPtFidCut", 10., "Fiducial max Pt for MC signal"};
  Configurable<double> fConfigMinEta{"cfgMinEtaFidCut", -0.8, "Fiducial min eta for MC signal"};
  Configurable<double> fConfigMaxEta{"cfgMaxEtaFidCut", 0.8, "Fiducial max eta for MC signal"};
  Configurable<bool> fConfigRecWithMC{"cfgEffRecWithMCVars", false, "If true, fill also 3D histograms at reconstructed level with mc variables"};
  Configurable<bool> fConfigFillLS{"cfgEffFillLS", false, "If true, fill LS histograms"};
  Configurable<bool> fConfigIsPrimary{"cfgIsPrimary", false, "If true, fill only with primary particles"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddSEPHistogram{"cfgAddSEPHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> fUseRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // 2D histos: mee and ptee
  ConfigurableAxis pteeBins{"pteeBins", {100, 0.f, 10.f}, "Ptee binning"};
  ConfigurableAxis meeBins{"meeBins", {600, 0.f, 3.5f}, "Mee binning"};

  // output lists
  OutputObj<THashList> fOutputQA{"DielectronQA"};
  HistogramRegistry registry{"HistoSameEventSelection", {}, OutputObjHandlingPolicy::AnalysisObject};
  THashList* fQADielectronList; // QA in case on with histo from manager outputs

  // Cuts and signals
  std::vector<AnalysisCompositeCut> fTrackCuts; // list of track cuts to init properly the histos
  std::vector<MCSignal> fMCSignals;             // list of signals with one prong to be checked: ULS 2D histos

  // 2D histo vectors
  std::vector<std::shared_ptr<TH2>> fHistGenPair;
  std::vector<std::shared_ptr<TH2>> fHistGenSmearedPair;
  std::vector<std::shared_ptr<TH2>> fHistRecPair;
  std::vector<std::shared_ptr<TH2>> fHistRecPairMC;

  // QA
  HistogramManager* fHistManQA; // histo manager

  void init(o2::framework::InitContext&)
  {
    fCurrentRun = 0;

    ccdb->setURL(ccdburl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    // Create list output for QA
    fQADielectronList = new THashList;
    fQADielectronList->SetOwner(kTRUE);
    fQADielectronList->SetName("DEQA");

    // Binning 2D histos
    AxisSpec axisPtee{pteeBins, "#it{p}_{T,ee} (GeV/#it{c})"};
    AxisSpec axisMee{meeBins, "#it{m}_{ee} (GeV/#it{c}^{2})"};

    // List of track cuts
    TString cutNamesStr = fConfigTrackCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    VarManager::SetDefaultVarNames();

    // List of MC signals
    TString configSigNamesStr = fConfigMCSignals.value;
    std::unique_ptr<TObjArray> sigNamesArray(configSigNamesStr.Tokenize(","));
    for (int isig = 0; isig < sigNamesArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(sigNamesArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() == 2) { // only 2 prong signals
          fMCSignals.push_back(*sig);
        }
        // List of signal to be checked
      }
    }

    // Configure 2D histograms
    // Generated particles
    for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
      fHistGenPair.push_back(registry.add<TH2>(Form("Dielectron/Generated/Ngen_Pair_ULS_%s", fMCSignals.at(i).GetName()), "", HistType::kTH2D, {axisMee, axisPtee}, true));
      if (fConfigFillLS) {
        fHistGenPair.push_back(registry.add<TH2>(Form("Dielectron/Generated/Ngen_Pair_LS_%s", fMCSignals.at(i).GetName()), "", HistType::kTH2D, {axisMee, axisPtee}, true));
      }
    }

    // Generated+smeared particles
    for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
      fHistGenSmearedPair.push_back(registry.add<TH2>(Form("Dielectron/GeneratedSmeared/NgenSmeared_Pair_ULS_%s", fMCSignals.at(i).GetName()), "", HistType::kTH2D, {axisMee, axisPtee}, true));
      if (fConfigFillLS) {
        fHistGenSmearedPair.push_back(registry.add<TH2>(Form("Dielectron/GeneratedSmeared/NgenSmeared_Pair_LS_%s", fMCSignals.at(i).GetName()), "", HistType::kTH2D, {axisMee, axisPtee}, true));
      }
    }

    // Rec with reconstructed variables
    for (unsigned int list_i = 0; list_i < fTrackCuts.size(); ++list_i) {
      for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
        fHistRecPair.push_back(registry.add<TH2>(Form("Dielectron/%s/Nrec_Pair_ULS_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH2D, {axisMee, axisPtee}, true));
        if (fConfigFillLS) {
          fHistRecPair.push_back(registry.add<TH2>(Form("Dielectron/%s/Nrec_Pair_LS_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH2D, {axisMee, axisPtee}, true));
        }
      }
    }

    // Rec with MC variables
    if (fConfigRecWithMC) {
      for (unsigned int list_i = 0; list_i < fTrackCuts.size(); ++list_i) {
        for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
          fHistRecPairMC.push_back(registry.add<TH2>(Form("Dielectron/%s_MCVars/Nrec_Pair_ULS_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH2D, {axisMee, axisPtee}, true));
          if (fConfigFillLS) {
            fHistRecPairMC.push_back(registry.add<TH2>(Form("Dielectron/%s_MCVars/Nrec_Pair_LS_%s", fTrackCuts.at(list_i).GetName(), fMCSignals.at(i).GetName()), "", HistType::kTH2D, {axisMee, axisPtee}, true));
          }
        }
      }
    }

    // Configure QA histogram classes
    if (fConfigQA) {
      std::vector<TString> names;
      TString histClassesQA = "";
      for (auto& cut : fTrackCuts) {

        // All passing dileptons
        names = {
          Form("PairsBarrelSEPM_%s", cut.GetName()),
          Form("PairsBarrelSEPP_%s", cut.GetName()),
          Form("PairsBarrelSEMM_%s", cut.GetName())};
        histClassesQA += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());

        // All reconstructed dileptons matched to a MC signal
        std::vector<TString> mcnamesreco;
        for (unsigned int isig = 0; isig < fMCSignals.size(); ++isig) {
          names = {
            Form("PairsBarrelSEPM_%s_%s", cut.GetName(), fMCSignals.at(isig).GetName()),
            Form("PairsBarrelSEPP_%s_%s", cut.GetName(), fMCSignals.at(isig).GetName()),
            Form("PairsBarrelSEMM_%s_%s", cut.GetName(), fMCSignals.at(isig).GetName())};
          histClassesQA += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
        }
      }

      fHistManQA = new HistogramManager("DielectronQA", "aa", VarManager::kNVars);
      fHistManQA->SetUseDefaultVariableNames(kTRUE);
      fHistManQA->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistManQA, histClassesQA.Data(), fConfigAddSEPHistogram); // define all histograms
      VarManager::SetUseVars(fHistManQA->GetUsedVars());                          // provide the list of required variables so that VarManager knows what to fill
      fQADielectronList = fHistManQA->GetMainHistogramList();
    }
    fOutputQA.setObject(fQADielectronList);
  }

  PresliceUnsorted<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;
  Preslice<MyBarrelTracks> perReducedEventTracks = aod::reducedtrack::reducedeventId;

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<MyBarrelTracksNoSkimmed> perCollisionTracks = aod::track::collisionId;

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTracks, typename TEventsMC, typename TTracksMC>
  void runPairing(TEvents const& events, TTracks const& tracks, TEventsMC const& /*eventsMC*/, TTracksMC const& tracksMC)
  {

    std::map<uint64_t, int> fMCEventLabels;
    int fCounters = 0; //! [0] - particle counter, [1] - event counter

    for (auto& event : events) {

      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::ResetValues(0, VarManager::kNMCParticleVariables);
      // fill event information which might be needed in histograms that combine track and event properties
      VarManager::FillEvent<TEventFillMap>(event);
      if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
        VarManager::FillEvent<TEventMCFillMap>(event.reducedMCevent());
      }
      if constexpr ((TEventMCFillMap & VarManager::ObjTypes::CollisionMC) > 0) {
        VarManager::FillEvent<TEventMCFillMap>(event.mcCollision());
      }

      // Look if we did not already saw the collision and fill the denominator of the single electron efficiency
      Int_t globalindexmc = -1;
      if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
        auto mcEvent = event.reducedMCevent();
        globalindexmc = mcEvent.globalIndex();
      }
      if constexpr ((TEventMCFillMap & VarManager::ObjTypes::CollisionMC) > 0) {
        auto mcEvent = event.mcCollision();
        globalindexmc = mcEvent.globalIndex();
      }
      if (!(fMCEventLabels.find(globalindexmc) != fMCEventLabels.end())) {
        fMCEventLabels[globalindexmc] = fCounters;
        fCounters++;
        // skimmed data
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
          auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
          groupedMCTracks.bindInternalIndicesTo(&tracksMC);
          runMCGenPair<false>(groupedMCTracks);
        }
        // Not skimmed data
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
          auto groupedMCTracks = tracksMC.sliceBy(perMcCollision, event.mcCollision().globalIndex());
          groupedMCTracks.bindInternalIndicesTo(&tracksMC);
          runMCGenPair<false>(groupedMCTracks);
        }
      }

      if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
        auto groupedTracks = tracks.sliceBy(perReducedEventTracks, event.globalIndex());
        groupedTracks.bindInternalIndicesTo(&tracks);
        runRecPair<TTrackFillMap>(groupedTracks, tracksMC);
      }
      if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
        auto groupedTracks = tracks.sliceBy(perCollisionTracks, event.globalIndex());
        runRecPair<TTrackFillMap>(groupedTracks, tracksMC);
      }
    } // end loop over reconstructed event
  } // end loop pairing function

  template <uint32_t TEventMCFillMap, uint32_t TTrackFillMap, typename TEventMC, typename TTracksMC>
  void runMCPairing(TEventMC const& /*eventMC*/, TTracksMC const& tracksMC)
  {

    runMCGenPair<true>(tracksMC);

  } // end loop pairing function

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks, typename TTracksMC>
  void runDataPairing(TEvent const& event, TTracks const& tracks, TTracksMC const& tracksMC)
  {

    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    // fill event information which might be needed in histograms that combine track and event properties
    VarManager::FillEvent<TEventFillMap>(event);
    runRecPair<TTrackFillMap>(tracks, tracksMC);

  } // end loop pairing function

  template <bool smeared, typename TTracksMC>
  void runMCGenPair(TTracksMC const& groupedMCTracks)
  {
    //
    Double_t masse = 0.00051099895; // 0.5 MeV/c2 -> 0.0005 GeV/c2

    for (auto& [t1, t2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(groupedMCTracks, groupedMCTracks))) {
      // for (auto& [t1, t2] : combinations(groupedMCTracks, groupedMCTracks)) {

      if ((abs(t1.pdgCode()) != 11) || (abs(t2.pdgCode()) != 11))
        continue;
      // if (!t1.producedByGenerator() || !t2.producedByGenerator())
      //   continue;
      if (fConfigIsPrimary) {
        if (!t1.isPhysicalPrimary() || !t2.isPhysicalPrimary())
          continue;
      }
      if (!fConfigFillLS && (t1.pdgCode() * t2.pdgCode() > 0))
        continue; // ULS only

      // True MC values
      TLorentzVector Lvec1;
      TLorentzVector Lvec2;
      Lvec1.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), masse);
      Lvec2.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), masse);
      TLorentzVector LvecM = Lvec1 + Lvec2;
      double mass = LvecM.M();
      double pairpt = LvecM.Pt();
      // Smeared MC values
      double masssmeared = -1.;
      double pairptsmeared = -1.;
      if constexpr (smeared) {
        TLorentzVector Lvec1smeared;
        TLorentzVector Lvec2smeared;
        Lvec1smeared.SetPtEtaPhiM(t1.ptSmeared(), t1.etaSmeared(), t1.phiSmeared(), masse);
        Lvec2smeared.SetPtEtaPhiM(t2.ptSmeared(), t2.etaSmeared(), t2.phiSmeared(), masse);
        TLorentzVector LvecMsmeared = Lvec1smeared + Lvec2smeared;
        masssmeared = LvecMsmeared.M();
        pairptsmeared = LvecMsmeared.Pt();
      }

      // Fiducial cut MC value
      Bool_t genfidcut = kTRUE;
      if ((t1.eta() > fConfigMaxEta) || (t2.eta() > fConfigMaxEta) || (t1.eta() < fConfigMinEta) || (t2.eta() < fConfigMinEta) || (t1.pt() > fConfigMaxPt) || (t2.pt() > fConfigMaxPt) || (t1.pt() < fConfigMinPt) || (t2.pt() < fConfigMinPt))
        genfidcut = kFALSE;

      // Fiducial cut Smeared values
      Bool_t genfidcutsmeared = kTRUE;
      if constexpr (smeared) {
        if ((t1.etaSmeared() > fConfigMaxEta) || (t2.etaSmeared() > fConfigMaxEta) || (t1.etaSmeared() < fConfigMinEta) || (t2.etaSmeared() < fConfigMinEta) || (t1.ptSmeared() > fConfigMaxPt) || (t2.ptSmeared() > fConfigMaxPt) || (t1.ptSmeared() < fConfigMinPt) || (t2.ptSmeared() < fConfigMinPt))
          genfidcutsmeared = kFALSE;
      }

      int isig = 0;
      for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {
        bool checked = false;
        if constexpr (soa::is_soa_filtered_v<TTracksMC>) {
          auto t1_raw = groupedMCTracks.rawIteratorAt(t1.globalIndex());
          auto t2_raw = groupedMCTracks.rawIteratorAt(t2.globalIndex());
          checked = (*sig).CheckSignal(true, t1_raw, t2_raw);
        } else {
          checked = (*sig).CheckSignal(true, t1, t2);
        }
        if (checked) {

          // not smeared after fiducial cuts
          if (genfidcut) {
            if (!fConfigFillLS) {
              fHistGenPair[isig]->Fill(mass, pairpt);
            } else {
              if (t1.pdgCode() * t2.pdgCode() < 0) {
                fHistGenPair[isig * 2]->Fill(mass, pairpt);
              } else {
                fHistGenPair[isig * 2 + 1]->Fill(mass, pairpt);
              }
            }
          }
          // Smeared
          if constexpr (smeared) {
            if (genfidcutsmeared) {
              if (!fConfigFillLS) {
                fHistGenSmearedPair[isig]->Fill(masssmeared, pairptsmeared);
              } else {
                if (t1.pdgCode() * t2.pdgCode() < 0) {
                  fHistGenSmearedPair[isig * 2]->Fill(masssmeared, pairptsmeared);
                } else {
                  fHistGenSmearedPair[isig * 2 + 1]->Fill(masssmeared, pairptsmeared);
                }
              }
            }
          }
        }
      }
    } // end of true pairing loop
  } // end runMCGen

  template <uint32_t TTrackFillMap, typename TTracks, typename TTracksMC>
  void runRecPair(TTracks const& tracks, TTracksMC const& /*tracksMC*/)
  {
    //
    Double_t masse = 0.00051099895; // 0.5 MeV/c2 -> 0.0005 GeV/c2 , to be taken differently
    Bool_t uls = kTRUE;

    // Loop over two track combinations
    uint8_t twoTrackFilter = 0;

    for (auto& [t1, t2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {
      // for (auto& [t1, t2] : combinations(tracks, tracks)) {

      twoTrackFilter = uint32_t(t1.isBarrelSelected()) & uint32_t(t2.isBarrelSelected());

      if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }

      if (t1.sign() * t2.sign() > 0)
        uls = kFALSE;
      else
        uls = kTRUE;

      if (!fConfigFillLS && !uls)
        continue;

      // Fiducial cut in case they are not done in the reconstructed cuts
      Bool_t recfidcut = kTRUE;
      if ((t1.eta() > fConfigMaxEta) || (t2.eta() > fConfigMaxEta) || (t1.eta() < fConfigMinEta) || (t2.eta() < fConfigMinEta) || (t1.pt() > fConfigMaxPt) || (t2.pt() > fConfigMaxPt) || (t1.pt() < fConfigMinPt) || (t2.pt() < fConfigMinPt))
        recfidcut = kFALSE;

      //
      VarManager::ResetValues(0, VarManager::kNPairVariables);
      VarManager::FillPair<VarManager::kDecayToEE, TTrackFillMap>(t1, t2);

      // Fill the QA for all passing tracks
      if (fConfigQA) {
        for (unsigned int j = 0; j < fTrackCuts.size(); j++) {
          if (twoTrackFilter & (uint8_t(1) << j)) {
            if (!fConfigFillLS) {
              fHistManQA->FillHistClass(Form("PairsBarrelSEPM_%s", fTrackCuts.at(j).GetName()), VarManager::fgValues);
            } else {
              if (uls) {
                fHistManQA->FillHistClass(Form("PairsBarrelSEPM_%s", fTrackCuts.at(j).GetName()), VarManager::fgValues);
              } else {
                if (t1.sign() > 0) {
                  fHistManQA->FillHistClass(Form("PairsBarrelSEPP_%s", fTrackCuts.at(j).GetName()), VarManager::fgValues);
                } else {
                  fHistManQA->FillHistClass(Form("PairsBarrelSEMM_%s", fTrackCuts.at(j).GetName()), VarManager::fgValues);
                }
              }
            }
          }
        }
      }

      // run MC matching for this pair
      uint32_t mcDecision = 0;
      int isig = 0;
      for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {

        if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) { // for skimmed DQ model
          if (fConfigIsPrimary) {
            if (!t1.reducedMCTrack().isPhysicalPrimary() || !t2.reducedMCTrack().isPhysicalPrimary())
              continue;
          }
          if ((*sig).CheckSignal(true, t1.reducedMCTrack(), t2.reducedMCTrack())) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }

        if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
          if (t1.has_mcParticle() && t2.has_mcParticle()) {
            auto mct1 = t1.template mcParticle_as<aod::McParticles>();
            auto mct2 = t2.template mcParticle_as<aod::McParticles>();
            if (fConfigIsPrimary) {
              if (!mct1.isPhysicalPrimary() || !mct2.isPhysicalPrimary())
                continue;
            }
            if ((*sig).CheckSignal(true, mct1, mct2)) {
              mcDecision |= (uint32_t(1) << isig);
            }
          }
        }

      } // end of loop MC signals

      for (unsigned int i = 0; i < fMCSignals.size(); i++) {
        if (!(mcDecision & (uint32_t(1) << i))) {
          continue;
        }
        if (recfidcut) {
          for (unsigned int j = 0; j < fTrackCuts.size(); j++) {
            if (twoTrackFilter & (uint8_t(1) << j)) {
              if (!fConfigFillLS) {
                fHistRecPair[j * fMCSignals.size() + i]->Fill(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt]);
                if (fConfigQA)
                  fHistManQA->FillHistClass(Form("PairsBarrelSEPM_%s_%s", fTrackCuts.at(j).GetName(), fMCSignals.at(i).GetName()), VarManager::fgValues);
              } else {
                if (uls) {
                  fHistRecPair[j * (2 * fMCSignals.size()) + 2 * i]->Fill(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt]);
                  if (fConfigQA)
                    fHistManQA->FillHistClass(Form("PairsBarrelSEPM_%s_%s", fTrackCuts.at(j).GetName(), fMCSignals.at(i).GetName()), VarManager::fgValues);
                } else {
                  fHistRecPair[j * (2 * fMCSignals.size()) + 2 * i + 1]->Fill(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt]);
                  if (fConfigQA) {
                    if (t1.sign() > 0) {
                      fHistManQA->FillHistClass(Form("PairsBarrelSEPP_%s_%s", fTrackCuts.at(j).GetName(), fMCSignals.at(i).GetName()), VarManager::fgValues);
                    } else {
                      fHistManQA->FillHistClass(Form("PairsBarrelSEMM_%s_%s", fTrackCuts.at(j).GetName(), fMCSignals.at(i).GetName()), VarManager::fgValues);
                    }
                  }
                }
              }
            }
          }
        }

        if (fConfigRecWithMC) {
          Double_t mass = -1000.;
          Double_t pairpt = 1000.;
          // Fiducial cut in case they are not done in the reconstructed cuts
          Bool_t genfidcut = kTRUE;
          if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
            auto mctrack1 = t1.reducedMCTrack();
            auto mctrack2 = t2.reducedMCTrack();
            TLorentzVector Lvec2, Lvec1;
            Lvec1.SetPtEtaPhiM(mctrack1.pt(), mctrack1.eta(), mctrack1.phi(), masse);
            Lvec2.SetPtEtaPhiM(mctrack2.pt(), mctrack2.eta(), mctrack2.phi(), masse);
            TLorentzVector LvecM = Lvec1 + Lvec2;
            mass = LvecM.M();
            pairpt = LvecM.Pt();
            if ((mctrack1.eta() > fConfigMaxEta) || (mctrack2.eta() > fConfigMaxEta) || (mctrack1.eta() < fConfigMinEta) || (mctrack2.eta() < fConfigMinEta) || (mctrack1.pt() > fConfigMaxPt) || (mctrack2.pt() > fConfigMaxPt) || (mctrack1.pt() < fConfigMinPt) || (mctrack2.pt() < fConfigMinPt))
              genfidcut = kFALSE;
          }
          if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
            if (t1.has_mcParticle() && t2.has_mcParticle()) {
              auto mctrack1 = t1.template mcParticle_as<aod::McParticles>();
              auto mctrack2 = t2.template mcParticle_as<aod::McParticles>();
              TLorentzVector Lvec2, Lvec1;
              Lvec1.SetPtEtaPhiM(mctrack1.pt(), mctrack1.eta(), mctrack1.phi(), masse);
              Lvec2.SetPtEtaPhiM(mctrack2.pt(), mctrack2.eta(), mctrack2.phi(), masse);
              TLorentzVector LvecM = Lvec1 + Lvec2;
              mass = LvecM.M();
              pairpt = LvecM.Pt();
              if ((mctrack1.eta() > fConfigMaxEta) || (mctrack2.eta() > fConfigMaxEta) || (mctrack1.eta() < fConfigMinEta) || (mctrack2.eta() < fConfigMinEta) || (mctrack1.pt() > fConfigMaxPt) || (mctrack2.pt() > fConfigMaxPt) || (mctrack1.pt() < fConfigMinPt) || (mctrack2.pt() < fConfigMinPt))
                genfidcut = kFALSE;
            } else {
              genfidcut = kFALSE;
            }
          }
          if (genfidcut) {
            for (unsigned int j = 0; j < fTrackCuts.size(); j++) {
              if (twoTrackFilter & (uint8_t(1) << j)) {
                if (!fConfigFillLS) {
                  fHistRecPairMC[j * fMCSignals.size() + i]->Fill(mass, pairpt);
                } else {
                  if (uls) {
                    fHistRecPairMC[j * (2 * fMCSignals.size()) + 2 * i]->Fill(mass, pairpt);
                  } else {
                    fHistRecPairMC[j * (2 * fMCSignals.size()) + 2 * i + 1]->Fill(mass, pairpt);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  void processToEESkimmed(soa::Filtered<MyEventsSelected> const& events,
                          soa::Filtered<MyBarrelTracksSelected> const& tracks,
                          ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    if (events.size() > 0 && fCurrentRun != events.begin().runNumber()) {
      if (fUseRemoteField.value) {
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, events.begin().timestamp());
        if (grpmag != nullptr) {
          mMagField = grpmag->getNominalL3Field();
        } else {
          LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", events.begin().timestamp());
        }
        VarManager::SetMagneticField(mMagField);
      } else {
        VarManager::SetMagneticField(fConfigMagField.value);
      }
      fCurrentRun = events.begin().runNumber();
    }

    runPairing<gkEventFillMap, gkMCEventFillMap, gkTrackFillMap>(events, tracks, eventsMC, tracksMC);
  }

  void processDataToEESkimmed(soa::Filtered<MyEventsSelected>::iterator const& event,
                              soa::Filtered<MyBarrelTracksSelected> const& tracks,
                              ReducedMCTracks const& tracksMC)
  {

    if (fCurrentRun != event.runNumber()) {
      if (fUseRemoteField.value) {
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, event.timestamp());
        if (grpmag != nullptr) {
          mMagField = grpmag->getNominalL3Field();
        } else {
          LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", event.timestamp());
        }
        VarManager::SetMagneticField(mMagField);
      } else {
        VarManager::SetMagneticField(fConfigMagField.value);
      }
      fCurrentRun = event.runNumber();
    }

    runDataPairing<gkEventFillMap, gkTrackFillMap>(event, tracks, tracksMC);
  }

  void processMCToEESkimmed(soa::Filtered<MyMCEventsSelected> const& eventsMC,
                            MyMCReducedTracks const& tracksMC)
  {
    for (auto& eventMC : eventsMC) {
      auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, eventMC.globalIndex());
      groupedMCTracks.bindInternalIndicesTo(&tracksMC);
      runMCPairing<gkMCEventFillMap, gkTrackFillMap>(eventMC, groupedMCTracks);
    }
  }

  void processDataToEENoSkimmed(soa::Filtered<MyEventsSelectedNoSkimmed>::iterator const& event,
                                aod::BCsWithTimestamps const&,
                                soa::Filtered<MyBarrelTracksSelectedNoSkimmed> const& tracks,
                                aod::McParticles const& tracksMC)
  {
    auto bc = event.template bc_as<aod::BCsWithTimestamps>();
    if (fCurrentRun != bc.runNumber()) {
      if (fUseRemoteField.value) {
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
        if (grpmag != nullptr) {
          mMagField = grpmag->getNominalL3Field();
        } else {
          LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", bc.timestamp());
        }
        VarManager::SetMagneticField(mMagField);
      } else {
        VarManager::SetMagneticField(fConfigMagField.value);
      }
      fCurrentRun = bc.runNumber();
    }

    runDataPairing<gkEventFillMapNoSkimmed, gkTrackFillMapNoSkimmed>(event, tracks, tracksMC);
  }

  void processMCToEENoSkimmed(soa::Filtered<MyMCEventsSelectedNoSkimmed>::iterator const& eventMC, MyMCTrackNoSkimmed const& tracksMC)
  {
    runMCPairing<gkMCEventFillMapNoSkimmed, gkTrackFillMapNoSkimmed>(eventMC, tracksMC);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  void processDummyNoSkimmed(MyEventsNoSkimmed&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processToEESkimmed, "Run the data and mc barrel barrel pairing on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDataToEESkimmed, "Run the data barrel barrel pairing on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMCToEESkimmed, "Run the MC barrel barrel pairing on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDataToEENoSkimmed, "Run the data for the barrel pairing on not skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMCToEENoSkimmed, "Run the MC for the barrel pairing on not skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy process function", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummyNoSkimmed, "Dummy process function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisEventQa>(cfgc),
    adaptAnalysisTask<AnalysisMCEvent>(cfgc),
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
      }
    }

    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }
    if (classStr.Contains("MCTruthGenPair")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_pair");
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Pt", "MC generator p_{T} distribution", false, 200, 0.0, 20.0, VarManager::kMCPt);
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Eta", "MC generator #eta distribution", false, 500, -5.0, 5.0, VarManager::kMCEta);
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Phi", "MC generator #varphi distribution", false, 500, -6.3, 6.3, VarManager::kMCPhi);
    }
    if (classStr.Contains("MCTruthGen")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_track");
    }
    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_barrel");
    }
  } // end loop over histogram classes
}
