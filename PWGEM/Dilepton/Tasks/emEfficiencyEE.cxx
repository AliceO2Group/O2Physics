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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include <TMath.h>
#include <TH1F.h>
#include <THashList.h>
#include <TLorentzVector.h>
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
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelected, isBarrelSelected, int);
} // namespace emanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "EVENTCUTS", emanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "BARRELTRACKCUTS", emanalysisflags::IsBarrelSelected);
} // namespace o2::aod

// No skimming: works for events and single tracks
using MyEventsNoSkimmed = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using MyEventsSelectedNoSkimmed = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::EventCuts>;
using MyBarrelTracksNoSkimmed = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                          aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                          aod::pidTPCFullKa, aod::pidTPCFullPr,
                                          aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                          aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                          aod::McTrackLabels>;
using MyBarrelTracksSelectedNoSkimmed = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                  aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                                  aod::pidTPCFullKa, aod::pidTPCFullPr,
                                                  aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                                  aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                                  aod::BarrelTrackCuts, aod::McTrackLabels>;

constexpr static uint32_t gkEventFillMapNoSkimmed = VarManager::ObjTypes::Collision;
constexpr static uint32_t gkMCEventFillMapNoSkimmed = VarManager::ObjTypes::CollisionMC;
constexpr static uint32_t gkTrackFillMapNoSkimmed = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkParticleMCFillMapNoSkimmed = VarManager::ObjTypes::ParticleMC;

// Skimmed data: works up to dielectron efficiency
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedMCEventLabels>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::ReducedMCEventLabels>;
using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::ReducedTracksBarrelLabels>;
// using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::ReducedTracksBarrelLabels>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::ReducedTracksBarrelLabels>;
// using MyBarrelTracksSelectedWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::ReducedTracksBarrelLabels>;

//
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkMCEventFillMap = VarManager::ObjTypes::ReducedEventMC;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
// constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkParticleMCFillMap = VarManager::ObjTypes::ParticleMC;

void DefineHistograms(HistogramManager* histMan, TString histClasses);
void SetBinsLinear(std::vector<double>& fBins, const double min, const double max, const unsigned int steps);

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
    bool pass = true;

    VarManager::FillEvent<TEventFillMap>(event);
    if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
      VarManager::FillEvent<TEventMCFillMap>(event.reducedMCevent());
    }
    if constexpr ((TEventMCFillMap & VarManager::ObjTypes::CollisionMC) > 0) {
      if (!event.has_mcCollision()) {
        pass = false;
      } else {
        VarManager::FillEvent<TEventMCFillMap>(event.mcCollision());
      }
    }
    if (fConfigQA) {
      fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues); // automatically fill all the histograms in the class Event
    }
    if (fEventCut->IsSelected(VarManager::fgValues) && pass) {
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
  OutputObj<THashList> fOutputList{"output"};
  THashList* fMainList;
  TH1D* fNbRecCollisionPerMCCollision;
  TProfile* fNbmctrack;
  TH1D* fNbMcEvent;
  TH1D* fNbRecEvent;

  void init(o2::framework::InitContext&)
  {

    VarManager::SetUseVars(AnalysisCut::fgUsedVars);
    VarManager::SetDefaultVarNames();

    fMainList = new THashList;
    fMainList->SetOwner(kTRUE);
    fMainList->SetName("OAnalysisEventQA");
    fNbRecCollisionPerMCCollision = new TH1D("MCEvRecEv", "", 20, 0., 20.);
    fNbmctrack = new TProfile("mctrack", "", 20, 0., 20.);
    fNbMcEvent = new TH1D("MCEvent", "", 1, 0., 1.);
    fNbRecEvent = new TH1D("RecEvent", "", 1, 0., 1.);
    fMainList->Add(fNbRecCollisionPerMCCollision);
    fMainList->Add(fNbmctrack);
    fMainList->Add(fNbMcEvent);
    fMainList->Add(fNbRecEvent);
    fOutputList.setObject(fMainList);
  }

  Preslice<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;
  Preslice<aod::McParticles_001> perMcCollision = aod::mcparticle::mcCollisionId;

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackMCFillMap, typename TEvents, typename TEventsMC, typename TTracksMC>
  void runSelection(TEvents const& events, TEventsMC const& eventsMC, TTracksMC const& tracksMC)
  {

    uint8_t eventFilter = 0;
    std::map<uint64_t, int> fMCEventNbmctrack;
    std::map<uint64_t, int> fMCEventNbReco;
    std::map<uint64_t, int> fMCEventLabels;

    int fMCCounters = 0;
    int fEvCounters = 0;

    // First loop

    for (auto& event : events) {
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEvent<TEventFillMap>(event);
      eventFilter = uint32_t(event.isEventSelected());
      if (!eventFilter)
        continue;
      fEvCounters++;
      fNbRecEvent->Fill(0.5);

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
        fNbMcEvent->Fill(0.5);
        fMCCounters++;
      } else {
        fMCEventNbReco[globalindexmc] = fMCEventNbReco.find(globalindexmc)->second + 1;
      }

    } // end loop over events

    for (const auto& [mcEv, NbRecEv] : fMCEventNbReco) {
      fNbRecCollisionPerMCCollision->Fill(fMCEventNbReco.find(mcEv)->second);
    }

    for (const auto& [mcEv, NbRecEv] : fMCEventNbmctrack) {
      fNbmctrack->Fill(fMCEventNbReco.find(mcEv)->second, fMCEventNbmctrack.find(mcEv)->second);
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

  void processNoSkimmed(soa::Filtered<MyEventsSelectedNoSkimmed> const& events, aod::McCollisions const& eventsMC, aod::McParticles_001 const& tracksMC)
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

  PROCESS_SWITCH(AnalysisEventQa, processNoSkimmed, "Run event selection without skimming", false);
  PROCESS_SWITCH(AnalysisEventQa, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventQa, processDummy, "Dummy process function", false);
  PROCESS_SWITCH(AnalysisEventQa, processDummyNoSkimmed, "Dummy process function", false);
};

struct AnalysisTrackSelection {

  Produces<aod::BarrelTrackCuts> trackSel;
  Filter filterEventSelected = aod::emanalysisflags::isEventSelected == 1;

  // configurables
  // Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandardNoINT7", "Event selection"};
  Configurable<std::string> fConfigCuts{"cfgTrackCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMCSignals{"cfgTrackMCSignals", "", "Comma separated list of MC signals"};

  // 3D histos for efficiency
  Configurable<bool> fConfigUsePtVec{"cfgUsePtVecEff", true, "If true, non-linear pt bins predefined"};
  Configurable<double> fConfigMinPt{"cfgMinPtEff", 0., "min Pt in 3D histos"};
  Configurable<double> fConfigMaxPt{"cfgMaxPtEff", 10., "max Pt in 3D histos"};
  Configurable<int> fConfigStepPt{"cfgStepPtEff", 1, "Nb of steps in pt in 3D histos"};
  Configurable<double> fConfigMinEta{"cfgMinEtaEff", -0.8, "min Eta in 3D histos"};
  Configurable<double> fConfigMaxEta{"cfgMaxEtaEff", 0.8, "max Eta in 3D histos"};
  Configurable<int> fConfigStepEta{"cfgStepEtaEff", 16, "Nb of steps in Eta in 3D histos"};
  Configurable<double> fConfigMinPhi{"cfgMinPhiEff", 0., "min Phi in 3D histos"};
  Configurable<double> fConfigMaxPhi{"cfgMaxPhiEff", 6.3, "max Phi in 3D histos"};
  Configurable<int> fConfigStepPhi{"cfgStepPhiEff", 63, "Nb of steps in Phi in 3D histos"};

  Configurable<bool> fConfigRecWithMC{"cfgEffRecWithMCVars", false, "If true, fill also 3D histograms at reconstructed level with mc variables"};

  // Resolution histos
  Configurable<bool> fConfigResolutionOn{"cfgResolution", false, "If true, fill resolution histograms"};
  Configurable<bool> fConfigUsePtVecRes{"cfgUsePtVecRes", true, "If true, non-linear pt bins predefined in res histos"};
  Configurable<double> fConfigMinPtRes{"cfgMinPtRes", 0., "min Pt in res histos"};
  Configurable<double> fConfigMaxPtRes{"cfgMaxPtRes", 20., "max Pt in res histos"};
  Configurable<int> fConfigStepPtRes{"cfgStepPtRes", 1, "Nb of steps in pt in res histos"};
  Configurable<int> fConfigStepDeltaPt{"cfgStepDeltaPtRes", 1, "Nb of steps in delta pt in res histos"};
  Configurable<double> fConfigMinDeltaEta{"cfgMinDeltaEtaRes", -0.5, "min delta Eta in res histos"};
  Configurable<double> fConfigMaxDeltaEta{"cfgMaxDeltaEtaRes", 0.5, "max delta Eta in res histos"};
  Configurable<int> fConfigStepDeltaEta{"cfgStepDeltaEtaRes", 500, "Nb of steps in detla Eta in res histos"};
  Configurable<double> fConfigMinDeltaPhi{"cfgMinDeltaPhiRes", -0.5, "min delta Phi in res histos"};
  Configurable<double> fConfigMaxDeltaPhi{"cfgMaxDeltaPhiRes", 0.5, "max delta Phi in res histos"};
  Configurable<int> fConfigStepDeltaPhi{"cfgStepDeltaPhiRes", 500, "Nb of steps in delta Phi in res histos"};

  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};

  // output lists
  OutputObj<THashList> fOutputList{"output"};
  THashList* fMainList;             // Main list
  THashList* fSingleElectronList;   // 3D histos for MC and reconstructed signals
  THashList* fResolutionList;       // Resolution histograms
  THashList* fQASingleElectronList; // QA in case on with histo manager outputs

  // Cuts and signals
  // AnalysisCompositeCut* fEventCut; // Taken from event selection part
  std::vector<AnalysisCompositeCut> fTrackCuts; // list of track cuts
  AnalysisCompositeCut* fTrackCutsRes;          // track cut for resolution map
  std::vector<MCSignal> fMCSignals;             // list of signals to be checked
  MCSignal* fMCSignalRes;                       // signal for res

  // 3D histos
  std::vector<TH3D*> fHistGenPosPart;
  std::vector<TH3D*> fHistGenNegPart;
  std::vector<TH3D*> fHistGenSmearedPosPart;
  std::vector<TH3D*> fHistGenSmearedNegPart;
  std::vector<TH3D*> fHistRecPosPart;
  std::vector<TH3D*> fHistRecNegPart;
  std::vector<TH3D*> fHistRecPosPartMC;
  std::vector<TH3D*> fHistRecNegPartMC;
  // Binning
  std::vector<double> fPtBins;
  std::vector<double> fEtaBins;
  std::vector<double> fPhiBins;

  // Res histos
  std::vector<TH2D*> fHistRes;
  // Binning
  std::vector<double> fPtResBins;
  std::vector<double> fDeltaEtaBins;
  std::vector<double> fDeltaPhiBins;

  // QA
  HistogramManager* fHistManQA;                            // histo manager
  std::vector<TString> fHistNamesRecoQA;                   // list of histo names for all reconstructed tracks in histo manager
  std::vector<std::vector<TString>> fHistNamesMCMatchedQA; // list of histo names for reconstructed signals in histo manager
  std::vector<std::vector<TString>> fHistNamesMCQA;        // list of histo names for generated signals in histo manager

  void init(o2::framework::InitContext&)
  {

    // Create list output
    fMainList = new THashList;
    fMainList->SetOwner(kTRUE);
    fMainList->SetName("trackselection");

    // Create list output for 3D eta,phi,pt
    fSingleElectronList = new THashList;
    fSingleElectronList->SetOwner(kTRUE);
    fSingleElectronList->SetName("SingleElectron");

    // Binning 3D histos
    if (fConfigUsePtVec) {
      const Int_t Npt = 68;
      Double_t pte[Npt] = {0.00, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.20, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285, 0.29, 0.295, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80, 2.00, 2.40, 2.80, 3.20, 3.70, 4.50, 6.00, 8.00, 10.};
      std::vector<double> v_pte(pte, std::end(pte));
      fPtBins = v_pte;
    } else {
      SetBinsLinear(fPtBins, fConfigMinPt, fConfigMaxPt, fConfigStepPt);
    }
    SetBinsLinear(fEtaBins, fConfigMinEta, fConfigMaxEta, fConfigStepEta);
    SetBinsLinear(fPhiBins, fConfigMinPhi, fConfigMaxPhi, fConfigStepPhi);
    const int fNptBins = fPtBins.size() - 1;
    const int fNetaBins = fEtaBins.size() - 1;
    const int fNphiBins = fPhiBins.size() - 1;

    // Event cut: taken from event selection part
    // fEventCut = new AnalysisCompositeCut(true);
    // TString eventCutStr = fConfigEventCuts.value;
    // fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));

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

    // Configure 3D histograms
    // Create List with generated particles
    TList* Generated = new TList();
    Generated->SetOwner();
    Generated->SetName("Generated");
    for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
      TH3D* th3_tmp_pos = new TH3D(Form("Ngen_Pos_%s", fMCSignals.at(i).GetName()), ";p_{T};#eta;#varphi", fNptBins, fPtBins.data(), fNetaBins, fEtaBins.data(), fNphiBins, fPhiBins.data());
      th3_tmp_pos->Sumw2();
      fHistGenPosPart.push_back(th3_tmp_pos);
      Generated->Add(th3_tmp_pos);
      TH3D* th3_tmp_neg = new TH3D(Form("Ngen_Neg_%s", fMCSignals.at(i).GetName()), ";p_{T};#eta;#varphi", fNptBins, fPtBins.data(), fNetaBins, fEtaBins.data(), fNphiBins, fPhiBins.data());
      th3_tmp_neg->Sumw2();
      fHistGenNegPart.push_back(th3_tmp_neg);
      Generated->Add(th3_tmp_neg);
    }
    // Create List with generated+smeared particles
    TList* GeneratedSmeared = new TList();
    GeneratedSmeared->SetName("GeneratedSmeared");
    GeneratedSmeared->SetOwner();
    for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
      TH3D* th3_tmp_pos = new TH3D(Form("Ngen_Pos_%s", fMCSignals.at(i).GetName()), ";p_{T};#eta;#varphi", fNptBins, fPtBins.data(), fNetaBins, fEtaBins.data(), fNphiBins, fPhiBins.data());
      th3_tmp_pos->Sumw2();
      fHistGenSmearedPosPart.push_back(th3_tmp_pos);
      GeneratedSmeared->Add(th3_tmp_pos);
      TH3D* th3_tmp_neg = new TH3D(Form("Ngen_Neg_%s", fMCSignals.at(i).GetName()), ";p_{T};#eta;#varphi", fNptBins, fPtBins.data(), fNetaBins, fEtaBins.data(), fNphiBins, fPhiBins.data());
      th3_tmp_neg->Sumw2();
      fHistGenSmearedNegPart.push_back(th3_tmp_neg);
      GeneratedSmeared->Add(th3_tmp_neg);
    }

    fSingleElectronList->Add(Generated);
    fSingleElectronList->Add(GeneratedSmeared);

    // For every cutsetting one list and every MCsignal 2 histograms with pos and neg charge
    for (unsigned int list_i = 0; list_i < fTrackCuts.size(); ++list_i) {
      TList* list = new TList();
      list->SetName(fTrackCuts.at(list_i).GetName());
      list->SetOwner();

      for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
        TH3D* th3_tmp_pos = new TH3D(Form("Nrec_Pos_%s", fMCSignals.at(i).GetName()), ";p_{T};#eta;#varphi", fNptBins, fPtBins.data(), fNetaBins, fEtaBins.data(), fNphiBins, fPhiBins.data());
        th3_tmp_pos->Sumw2();
        th3_tmp_pos->SetDirectory(0x0);
        fHistRecPosPart.push_back(th3_tmp_pos);
        list->Add(th3_tmp_pos);
        TH3D* th3_tmp_neg = new TH3D(Form("Nrec_Neg_%s", fMCSignals.at(i).GetName()), ";p_{T};#eta;#varphi", fNptBins, fPtBins.data(), fNetaBins, fEtaBins.data(), fNphiBins, fPhiBins.data());
        th3_tmp_neg->Sumw2();
        th3_tmp_neg->SetDirectory(0x0);
        fHistRecNegPart.push_back(th3_tmp_neg);
        list->Add(th3_tmp_neg);
      }
      fSingleElectronList->Add(list);
    }
    // For every cutsetting one list and every MCsignal 2 histograms with pos and neg charge, filled with mc variables
    if (fConfigRecWithMC) {
      for (unsigned int list_i = 0; list_i < fTrackCuts.size(); ++list_i) {
        TList* list = new TList();
        list->SetName(Form("%s_MCVars", fTrackCuts.at(list_i).GetName()));
        list->SetOwner();

        for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
          TH3D* th3_tmp_pos = new TH3D(Form("Nrec_Pos_%s", fMCSignals.at(i).GetName()), ";p_{T};#eta;#varphi", fNptBins, fPtBins.data(), fNetaBins, fEtaBins.data(), fNphiBins, fPhiBins.data());
          th3_tmp_pos->Sumw2();
          th3_tmp_pos->SetDirectory(0x0);
          fHistRecPosPartMC.push_back(th3_tmp_pos);
          list->Add(th3_tmp_pos);
          TH3D* th3_tmp_neg = new TH3D(Form("Nrec_Neg_%s", fMCSignals.at(i).GetName()), ";p_{T};#eta;#varphi", fNptBins, fPtBins.data(), fNetaBins, fEtaBins.data(), fNphiBins, fPhiBins.data());
          th3_tmp_neg->Sumw2();
          th3_tmp_neg->SetDirectory(0x0);
          fHistRecNegPartMC.push_back(th3_tmp_neg);
          list->Add(th3_tmp_neg);
        }
        fSingleElectronList->Add(list);
      }
    }
    fMainList->Add(fSingleElectronList);

    // Resolution histogramms
    if (fConfigResolutionOn) {

      // Binning 3D histos
      if (fConfigUsePtVecRes) {
        const Int_t Npt = 73;
        Double_t pte[Npt] = {0.00, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.20, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285, 0.29, 0.295, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80, 2.00, 2.40, 2.80, 3.20, 3.70, 4.50, 6.00, 8.00, 10., 12.0, 14., 16., 18., 20.};
        std::vector<double> v_pte(pte, std::end(pte));
        fPtResBins = v_pte;
      } else {
        SetBinsLinear(fPtResBins, fConfigMinPtRes, fConfigMaxPtRes, fConfigStepPtRes);
      }
      SetBinsLinear(fDeltaEtaBins, fConfigMinDeltaEta, fConfigMaxDeltaEta, fConfigStepDeltaEta);
      SetBinsLinear(fDeltaPhiBins, fConfigMinDeltaPhi, fConfigMaxDeltaPhi, fConfigStepDeltaPhi);
      const int fNptresBins = fPtResBins.size() - 1;
      const int fNDeltaetaBins = fDeltaEtaBins.size() - 1;
      const int fNDeltaphiBins = fDeltaPhiBins.size() - 1;

      fResolutionList = new THashList;
      fResolutionList->SetOwner(kTRUE);
      fResolutionList->SetName("Resolution");

      TH2D* thPtGen_DeltaPtOverPtGen = new TH2D("PtGen_DeltaPtOverPtGen", "", fNptresBins, fPtResBins.data(), fConfigStepDeltaPt, -1., +1.);
      TH2D* thPtGen_DeltaEta = new TH2D("PtGen_DeltaEta", "", fNptresBins, fPtResBins.data(), fNDeltaetaBins, fDeltaEtaBins.data());
      TH2D* thPtGen_DeltaPhi_Ele = new TH2D("PtGen_DeltaPhi_Ele", "", fNptresBins, fPtResBins.data(), fNDeltaphiBins, fDeltaPhiBins.data());
      TH2D* thPtGen_DeltaPhi_Pos = new TH2D("PtGen_DeltaPhi_Pos", "", fNptresBins, fPtResBins.data(), fNDeltaphiBins, fDeltaPhiBins.data());

      thPtGen_DeltaPtOverPtGen->Sumw2();
      thPtGen_DeltaEta->Sumw2();
      thPtGen_DeltaPhi_Ele->Sumw2();
      thPtGen_DeltaPhi_Pos->Sumw2();

      thPtGen_DeltaPtOverPtGen->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      thPtGen_DeltaPtOverPtGen->GetYaxis()->SetTitle("(p^{gen}_{T} - p^{rec}_{T}) / p^{gen}_{T} (GeV/c)");
      thPtGen_DeltaEta->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      thPtGen_DeltaEta->GetYaxis()->SetTitle("#eta^{gen} - #eta^{rec}");
      thPtGen_DeltaPhi_Ele->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      thPtGen_DeltaPhi_Ele->GetYaxis()->SetTitle("#varphi^{gen} - #varphi^{rec} (rad)");
      thPtGen_DeltaPhi_Pos->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      thPtGen_DeltaPhi_Pos->GetYaxis()->SetTitle("#varphi^{gen} - #varphi^{rec} (rad)");

      fHistRes.push_back(thPtGen_DeltaPtOverPtGen);
      fHistRes.push_back(thPtGen_DeltaEta);
      fHistRes.push_back(thPtGen_DeltaPhi_Ele);
      fHistRes.push_back(thPtGen_DeltaPhi_Pos);

      fResolutionList->Add(thPtGen_DeltaPtOverPtGen);
      fResolutionList->Add(thPtGen_DeltaEta);
      fResolutionList->Add(thPtGen_DeltaPhi_Ele);
      fResolutionList->Add(thPtGen_DeltaPhi_Pos);

      fMainList->Add(fResolutionList);
    }

    // Configure QA histogram classes
    if (fConfigQA) {
      // Create list output for QA
      fQASingleElectronList = new THashList;
      fQASingleElectronList->SetOwner(kTRUE);
      fQASingleElectronList->SetName("SingleElectronQA");
      TString histClassesQA = "";
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
      std::vector<TString> mcnamesgen;
      for (unsigned int isig = 0; isig < fMCSignals.size(); ++isig) {
        TString nameStr2 = Form("MCTruthGen_%s", fMCSignals.at(isig).GetName());
        mcnamesgen.push_back(nameStr2);
        histClassesQA += Form("%s;", nameStr2.Data());
      }
      fHistNamesMCQA.push_back(mcnamesgen);

      fHistManQA = new HistogramManager("SingleElectronQA", "aa", VarManager::kNVars);
      fHistManQA->SetUseDefaultVariableNames(kTRUE);
      fHistManQA->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistManQA, histClassesQA.Data()); // define all histograms
      VarManager::SetUseVars(fHistManQA->GetUsedVars());  // provide the list of required variables so that VarManager knows what to fill
      fQASingleElectronList = fHistManQA->GetMainHistogramList();
      fMainList->Add(fQASingleElectronList);
    }
    fOutputList.setObject(fMainList);
  }

  Preslice<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;
  Preslice<MyBarrelTracks> perReducedEventTracks = aod::reducedtrack::reducedeventId;

  Preslice<aod::McParticles_001> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<MyBarrelTracksNoSkimmed> perCollisionTracks = aod::track::collisionId;

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackFillMap, uint32_t TTrackMCFillMap, typename TEvents, typename TTracks, typename TEventsMC, typename TTracksMC>
  void runSelection(TEvents const& events, TTracks const& tracks, TEventsMC const& eventsMC, TTracksMC const& tracksMC)
  {

    uint8_t eventFilter = 0;
    std::map<uint64_t, int> fMCEventLabels;
    int fCounters = 0; //! [0] - particle counter, [1] - event counter

    for (auto& event : events) {
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::ResetValues(0, VarManager::kNMCParticleVariables);
      // fill event information which might be needed in histograms that combine track and event properties
      VarManager::FillEvent<TEventFillMap>(event);
      // if(!fEventCut->IsSelected(VarManager::fgValues)) continue;
      eventFilter = uint32_t(event.isEventSelected());
      if (!eventFilter)
        continue;
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
          runMCGenTrack(groupedMCTracks);
        }
        // Not skimmed data
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
          auto groupedMCTracks = tracksMC.sliceBy(perMcCollision, event.mcCollision().globalIndex());
          groupedMCTracks.bindInternalIndicesTo(&tracksMC);
          runMCGenTrack(groupedMCTracks);
        }
      }

      // Loop over reconstructed tracks belonging to the event and fill the numerator of the efficiency as well as the resolution map
      if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
        auto groupedTracks = tracks.sliceBy(perReducedEventTracks, event.globalIndex());
        runRecTrack<gkTrackFillMap>(groupedTracks, tracksMC);
      }
      if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
        auto groupedTracks = tracks.sliceBy(perCollisionTracks, event.globalIndex());
        runRecTrack<gkTrackFillMapNoSkimmed>(groupedTracks, tracksMC);
      }
    } // end loop over events
  }

  template <typename TTracksMC>
  void runMCGenTrack(TTracksMC const& groupedMCTracks)
  {

    // loop over mc stack and fill histograms for pure MC truth signals
    // group all the MC tracks which belong to the MC event corresponding to the current reconstructed event
    // auto groupedMCTracks = tracksMC.sliceBy(aod::reducedtrackMC::reducedMCeventId, event.reducedMCevent().globalIndex());
    for (auto& mctrack : groupedMCTracks) {
      VarManager::FillTrack<gkParticleMCFillMap>(mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      int isig = 0;
      for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {
        if ((*sig).CheckSignal(true, groupedMCTracks, mctrack)) {
          if (mctrack.pdgCode() > 0) {
            dynamic_cast<TH3D*>(fHistGenNegPart.at(isig))->Fill(mctrack.pt(), mctrack.eta(), mctrack.phi());
          } else {
            dynamic_cast<TH3D*>(fHistGenPosPart.at(isig))->Fill(mctrack.pt(), mctrack.eta(), mctrack.phi());
          }
          if (fConfigQA)
            fHistManQA->FillHistClass(Form("MCTruthGen_%s", (*sig).GetName()), VarManager::fgValues);
        }
      }
    }
  }

  template <uint32_t TTrackFillMap, typename TTracks, typename TTracksMC>
  void runRecTrack(TTracks const& groupedTracks, TTracksMC const& tracksMC)
  {

    uint32_t filterMap = 0;
    trackSel.reserve(groupedTracks.size());
    for (auto& track : groupedTracks) {
      filterMap = 0;

      VarManager::FillTrack<TTrackFillMap>(track); // compute track quantities

      // compute MC matched quantities
      if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
        VarManager::FillTrack<gkParticleMCFillMap>(track.reducedMCTrack());
      }
      if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
        // If no MC particle is found, skip the track
        if (track.has_mcParticle()) {
          auto mctrack = track.template mcParticle_as<aod::McParticles_001>();
          VarManager::FillTrack<gkParticleMCFillMapNoSkimmed>(mctrack);
        }
      }

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
      trackSel(static_cast<int>(filterMap));
      if (!filterMap) {
        continue;
      }

      // compute MC matching decisions
      uint32_t mcDecision = 0;
      int isig = 0;
      for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {

        if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
          if ((*sig).CheckSignal(true, tracksMC, track.reducedMCTrack())) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::Track) > 0) {
          if (track.has_mcParticle()) {
            auto mctrack = track.template mcParticle_as<aod::McParticles_001>();
            if ((*sig).CheckSignal(true, tracksMC, mctrack)) {
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
              dynamic_cast<TH3D*>(fHistRecNegPart.at(j * fMCSignals.size() + i))->Fill(track.pt(), track.eta(), track.phi());
            } else {
              dynamic_cast<TH3D*>(fHistRecPosPart.at(j * fMCSignals.size() + i))->Fill(track.pt(), track.eta(), track.phi());
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
                  auto mctrack = track.template mcParticle_as<aod::McParticles_001>();
                  mcpt = mctrack.pt();
                  mceta = mctrack.eta();
                  mcphi = mctrack.phi();
                }
              }

              if (track.sign() < 0) {
                dynamic_cast<TH3D*>(fHistRecNegPartMC.at(j * fMCSignals.size() + i))->Fill(mcpt, mceta, mcphi);
              } else {
                dynamic_cast<TH3D*>(fHistRecPosPartMC.at(j * fMCSignals.size() + i))->Fill(mcpt, mceta, mcphi);
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
                  auto mctrack = track.template mcParticle_as<aod::McParticles_001>();
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
              dynamic_cast<TH2D*>(fHistRes.at(0))->Fill(mcpt, deltaptoverpt);
              dynamic_cast<TH2D*>(fHistRes.at(1))->Fill(mcpt, deltaeta);
              if (mcpdg < 0) {
                dynamic_cast<TH2D*>(fHistRes.at(2))->Fill(mcpt, deltaphi);
              } else {
                dynamic_cast<TH2D*>(fHistRes.at(3))->Fill(mcpt, deltaphi);
              }
            }
            if (fConfigQA)
              fHistManQA->FillHistClass(fHistNamesMCMatchedQA[j][i].Data(), VarManager::fgValues);
          }
        } // end loop over cuts
      }   // end loop over MC signals
    }     // end loop over reconstructed track belonging to the events
  }

  void processSkimmed(soa::Filtered<MyEventsSelected> const& events, MyBarrelTracks const& tracks, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runSelection<gkEventFillMap, gkMCEventFillMap, gkTrackFillMap, gkParticleMCFillMap>(events, tracks, eventsMC, tracksMC);
  }

  void processNoSkimmed(soa::Filtered<MyEventsSelectedNoSkimmed> const& events, MyBarrelTracksNoSkimmed const& tracks, aod::McCollisions const& eventsMC, aod::McParticles_001 const& tracksMC)
  {
    runSelection<gkEventFillMapNoSkimmed, gkMCEventFillMapNoSkimmed, gkTrackFillMapNoSkimmed, gkParticleMCFillMapNoSkimmed>(events, tracks, eventsMC, tracksMC);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  void processDummyNoSkimmed(MyEventsNoSkimmed&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisTrackSelection, processNoSkimmed, "Run barrel track selection without skimming", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmed, "Run barrel track selection on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummy, "Dummy process function", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummyNoSkimmed, "Dummy process function", false);
};

struct AnalysisSameEventPairing {

  // Produces<aod::Dileptons> dileptonList;

  // Filter based on previous components in the task
  Filter filterEventSelected = aod::emanalysisflags::isEventSelected == 1;
  Filter filterBarrelTrackSelected = aod::emanalysisflags::isBarrelSelected > 0;

  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMCSignals{"cfgBarrelMCSignals", "", "Comma separated list of MC signals"};
  Configurable<double> fConfigMinPt{"cfgMinPtFidCut", 0., "Fiducial min Pt for MC signal"};
  Configurable<double> fConfigMaxPt{"cfgMaxPtFidCut", 10., "Fiducial max Pt for MC signal"};
  Configurable<double> fConfigMinEta{"cfgMinEtaFidCut", -0.8, "Fiducial min eta for MC signal"};
  Configurable<double> fConfigMaxEta{"cfgMaxEtaFidCut", 0.8, "Fiducial max eta for MC signal"};
  Configurable<bool> fConfigFlatTables{"cfgFlatTables", false, "Produce a single flat tables with all relevant information of the pairs and single tracks"};
  Configurable<bool> fConfigRecWithMC{"cfgEffRecWithMCVars", false, "If true, fill also 3D histograms at reconstructed level with mc variables"};
  Configurable<bool> fConfigFillLS{"cfgEffFillLS", false, "If true, fill LS histograms"};

  // TODO: here we specify signals, however signal decisions are precomputed and stored in mcReducedFlags
  // TODO: The tasks based on skimmed MC could/should rely ideally just on these flags
  // TODO:   special AnalysisCuts to be prepared in this direction
  // TODO: cuts on the MC truth information to be added if needed

  // 2D histos: mee and ptee
  Configurable<double> fConfigMinPtee{"cfgMinPtee", 0., "min Ptee in 2D histos"};
  Configurable<double> fConfigMaxPtee{"cfgMaxPtee", 10., "max Ptee in 2D histos"};
  Configurable<int> fConfigStepPtee{"cfgStepPtee", 100, "Nb of steps in ptee in 2D histos"};
  Configurable<double> fConfigMinMee{"cfgMinMee", 0., "min Mee in 2D histos"};
  Configurable<double> fConfigMaxMee{"cfgMaxMee", 3.5, "max Mee in 2D histos"};
  Configurable<int> fConfigStepMee{"cfgStepMee", 600, "Nb of steps in Mee in 2D histos"};

  // output lists
  OutputObj<THashList> fOutputList{"output"};
  THashList* fMainList;   // Main list
  THashList* fPairList;   // 2D histos for MC and reconstructed signals
  THashList* fQAPairList; // QA in case on with histo manager outputs

  // Cuts and signals
  // AnalysisCompositeCut* fEventCut; // Taken from event selection part
  std::vector<AnalysisCompositeCut> fTrackCuts; // list of track cuts
  std::vector<MCSignal> fMCSignals;             // list of signals with one prong to be checked: ULS 2D histos

  // 2D histo vectors
  std::vector<TH2D*> fHistGenPair;
  std::vector<TH2D*> fHistGenSmearedPair;
  std::vector<TH2D*> fHistRecPair;
  std::vector<TH2D*> fHistRecPairMC;
  // Binning
  std::vector<double> fPteeBins;
  std::vector<double> fMeeBins;

  // QA: to be defined

  void init(o2::framework::InitContext& context)
  {

    // Create list output
    fMainList = new THashList;
    fMainList->SetOwner(kTRUE);
    fMainList->SetName("pairselection");

    // Create list output for 3D eta,phi,pt
    fPairList = new THashList;
    fPairList->SetOwner(kTRUE);
    fPairList->SetName("Dielectron");

    // Binning 2D histos
    SetBinsLinear(fPteeBins, fConfigMinPtee, fConfigMaxPtee, fConfigStepPtee);
    SetBinsLinear(fMeeBins, fConfigMinMee, fConfigMaxMee, fConfigStepMee);
    const int fNpteeBins = fPteeBins.size() - 1;
    const int fNmeeBins = fMeeBins.size() - 1;

    // Event cut: taken from event selection part
    // fEventCut = new AnalysisCompositeCut(true);
    // TString eventCutStr = fConfigEventCuts.value;
    // fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));

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
    // Create List with generated particles
    TList* Generated = new TList();
    Generated->SetOwner();
    Generated->SetName("Generated");
    for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
      TH2D* th2_tmp = new TH2D(Form("Ngen_Pair_ULS_%s", fMCSignals.at(i).GetName()), ";m_{ee};p_{T,ee}", fNmeeBins, fMeeBins.data(), fNpteeBins, fPteeBins.data());
      th2_tmp->Sumw2();
      fHistGenPair.push_back(th2_tmp);
      Generated->Add(th2_tmp);
      if (fConfigFillLS) {
        TH2D* th3_tmp = new TH2D(Form("Ngen_Pair_LS_%s", fMCSignals.at(i).GetName()), ";m_{ee};p_{T,ee}", fNmeeBins, fMeeBins.data(), fNpteeBins, fPteeBins.data());
        th3_tmp->Sumw2();
        fHistGenPair.push_back(th3_tmp);
        Generated->Add(th3_tmp);
      }
    }

    // Create List with generated+smeared particles
    TList* GeneratedSmeared = new TList();
    GeneratedSmeared->SetName("GeneratedSmeared");
    GeneratedSmeared->SetOwner();
    for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
      TH2D* th2_tmp = new TH2D(Form("NgenSmeared_Pair_ULS_%s", fMCSignals.at(i).GetName()), ";m_{ee};p_{T,ee}", fNmeeBins, fMeeBins.data(), fNpteeBins, fPteeBins.data());
      th2_tmp->Sumw2();
      fHistGenSmearedPair.push_back(th2_tmp);
      GeneratedSmeared->Add(th2_tmp);
      if (fConfigFillLS) {
        TH2D* th3_tmp = new TH2D(Form("NgenSmeared_Pair_LS_%s", fMCSignals.at(i).GetName()), ";m_{ee};p_{T,ee}", fNmeeBins, fMeeBins.data(), fNpteeBins, fPteeBins.data());
        th3_tmp->Sumw2();
        fHistGenSmearedPair.push_back(th3_tmp);
        GeneratedSmeared->Add(th3_tmp);
      }
    }

    fPairList->Add(Generated);
    fPairList->Add(GeneratedSmeared);

    // Rec with reconstructed variables
    for (unsigned int list_i = 0; list_i < fTrackCuts.size(); ++list_i) {
      TList* list = new TList();
      list->SetName(fTrackCuts.at(list_i).GetName());
      list->SetOwner();

      for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
        TH2D* th2_tmp = new TH2D(Form("Nrec_Pair_ULS_%s", fMCSignals.at(i).GetName()), ";m_{ee};p_{T,ee}", fNmeeBins, fMeeBins.data(), fNpteeBins, fPteeBins.data());
        th2_tmp->Sumw2();
        th2_tmp->SetDirectory(0x0);
        fHistRecPair.push_back(th2_tmp);
        list->Add(th2_tmp);
        if (fConfigFillLS) {
          TH2D* th3_tmp = new TH2D(Form("Nrec_Pair_LS_%s", fMCSignals.at(i).GetName()), ";m_{ee};p_{T,ee}", fNmeeBins, fMeeBins.data(), fNpteeBins, fPteeBins.data());
          th3_tmp->Sumw2();
          th3_tmp->SetDirectory(0x0);
          fHistRecPair.push_back(th3_tmp);
          list->Add(th3_tmp);
        }
      }
      fPairList->Add(list);
    }

    // Rec with MC variables
    if (fConfigRecWithMC) {
      for (unsigned int list_i = 0; list_i < fTrackCuts.size(); ++list_i) {
        TList* list = new TList();
        list->SetName(Form("%s_MCVars", fTrackCuts.at(list_i).GetName()));
        list->SetOwner();

        for (unsigned int i = 0; i < fMCSignals.size(); ++i) {
          TH2D* th2_tmp = new TH2D(Form("Nrec_Pair_ULS_%s", fMCSignals.at(i).GetName()), ";m_{ee};p_{T,ee}", fNmeeBins, fMeeBins.data(), fNpteeBins, fPteeBins.data());
          th2_tmp->Sumw2();
          th2_tmp->SetDirectory(0x0);
          fHistRecPairMC.push_back(th2_tmp);
          list->Add(th2_tmp);
          if (fConfigFillLS) {
            TH2D* th3_tmp = new TH2D(Form("Nrec_Pair_LS_%s", fMCSignals.at(i).GetName()), ";m_{ee};p_{T,ee}", fNmeeBins, fMeeBins.data(), fNpteeBins, fPteeBins.data());
            th3_tmp->Sumw2();
            th3_tmp->SetDirectory(0x0);
            fHistRecPairMC.push_back(th3_tmp);
            list->Add(th3_tmp);
          }
        }
        fPairList->Add(list);
      }
    }
    fMainList->Add(fPairList);
    fOutputList.setObject(fMainList);
  }

  Preslice<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;
  Preslice<MyBarrelTracks> perReducedEventTracks = aod::reducedtrack::reducedeventId;

  template <uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTracks, typename TEventsMC, typename TTracksMC>
  void runPairing(TEvents const& events, TTracks const& tracks, TEventsMC const& eventsMC, TTracksMC const& tracksMC)
  {

    std::map<uint64_t, int> fMCEventLabels;
    int fCounters = 0; //! [0] - particle counter, [1] - event counter

    for (auto& event : events) {

      if (!event.isEventSelected()) {
        return;
      }

      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::ResetValues(0, VarManager::kNMCParticleVariables);
      // fill event information which might be needed in histograms that combine track and event properties
      VarManager::FillEvent<TEventFillMap>(event);
      if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
        VarManager::FillEvent<TEventMCFillMap>(event.reducedMCevent());
      }

      // Look if we did not already saw the collision and fill the denominator of the single electron efficiency
      Int_t globalindexmc = -1;
      if constexpr ((TEventMCFillMap & VarManager::ObjTypes::ReducedEventMC) > 0) {
        auto mcEvent = event.reducedMCevent();
        globalindexmc = mcEvent.globalIndex();
      }
      if (!(fMCEventLabels.find(globalindexmc) != fMCEventLabels.end())) {
        fMCEventLabels[globalindexmc] = fCounters;
        fCounters++;
        // skimmed data
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
          auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
          groupedMCTracks.bindInternalIndicesTo(&tracksMC);
          runMCGenPair(groupedMCTracks);
        }
      }

      // auto groupedTrackspos = posTracks->sliceByCached(aod::reducedtrack::reducedeventId, event.globalIndex());
      // groupedTrackspos.bindInternalIndicesTo(&tracks);
      // auto groupedTracksneg = negTracks->sliceByCached(aod::reducedtrack::reducedeventId, event.globalIndex());
      // groupedTracksneg.bindInternalIndicesTo(&tracks);

      if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) {
        auto groupedTracks = tracks.sliceBy(perReducedEventTracks, event.globalIndex());
        groupedTracks.bindInternalIndicesTo(&tracks);
        runRecPair<TTrackFillMap>(groupedTracks, tracksMC);
      }
    } // end loop over reconstructed event
  }   // end loop pairing function

  template <typename TTracksMC>
  void runMCGenPair(TTracksMC const& groupedMCTracks)
  {
    //
    Double_t masse = 0.00051099895; // 0.5 MeV/c2 -> 0.0005 GeV/c2

    for (auto& [t1, t2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(groupedMCTracks, groupedMCTracks))) {
      // for (auto& [t1, t2] : combinations(groupedMCTracks, groupedMCTracks)) {

      if ((abs(t1.pdgCode()) != 11) || (abs(t2.pdgCode()) != 11))
        continue;
      if (!fConfigFillLS && (t1.pdgCode() * t2.pdgCode() > 0))
        continue; // ULS only

      TLorentzVector Lvec1;
      TLorentzVector Lvec2;
      Lvec1.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), masse);
      Lvec2.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), masse);
      TLorentzVector LvecM = Lvec1 + Lvec2;
      double mass = LvecM.M();
      double pairpt = LvecM.Pt();
      // double opangle = Lvec1.Angle(Lvec2.Vect());

      // Fiducial cut
      Bool_t genfidcut = kTRUE;
      if ((t1.eta() > fConfigMaxEta) || (t2.eta() > fConfigMaxEta) || (t1.eta() < fConfigMinEta) || (t2.eta() < fConfigMinEta) || (t1.pt() > fConfigMaxPt) || (t2.pt() > fConfigMaxPt) || (t1.pt() < fConfigMinPt) || (t2.pt() < fConfigMinPt))
        genfidcut = kFALSE;

      int isig = 0;
      for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {
        if ((*sig).CheckSignal(true, groupedMCTracks, t1, t2)) {

          // not smeared after fiducial cuts
          if (genfidcut) {
            if (!fConfigFillLS)
              dynamic_cast<TH2D*>(fHistGenPair.at(isig))->Fill(mass, pairpt);
            else {
              if (t1.pdgCode() * t2.pdgCode() < 0) {
                dynamic_cast<TH2D*>(fHistGenPair.at(isig * 2))->Fill(mass, pairpt);
              } else {
                dynamic_cast<TH2D*>(fHistGenPair.at(isig * 2 + 1))->Fill(mass, pairpt);
              }
            }
          }
          // need to implement smeared
        }
      }
    } // end of true pairing loop
  }   // end runMCGen
  template <uint32_t TTrackFillMap, typename TTracks, typename TTracksMC>
  void runRecPair(TTracks const& tracks, TTracksMC const& tracksMC)
  {
    //
    Double_t masse = 0.00051099895; // 0.5 MeV/c2 -> 0.0005 GeV/c2 , to be taken differently
    Bool_t uls = kTRUE;

    // Loop over two track combinations
    uint8_t twoTrackFilter = 0;
    // uint32_t dileptonFilterMap = 0;
    // uint32_t dileptonMcDecision = 0;
    // dileptonList.reserve(1);

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
      VarManager::FillPair<VarManager::kJpsiToEE, TTrackFillMap>(t1, t2);

      // run MC matching for this pair
      uint32_t mcDecision = 0;
      int isig = 0;
      for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {

        if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrack) > 0) { // for skimmed DQ model
          if ((*sig).CheckSignal(true, tracksMC, t1.reducedMCTrack(), t2.reducedMCTrack())) {
            mcDecision |= (uint32_t(1) << isig);
          }
        }

      } // end of loop MC signals

      // dileptonFilterMap = twoTrackFilter;
      // dileptonMcDecision = mcDecision;
      // dileptonList(event, VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), dileptonFilterMap, dileptonMcDecision);

      for (unsigned int i = 0; i < fMCSignals.size(); i++) {
        if (!(mcDecision & (uint32_t(1) << i))) {
          continue;
        }
        if (recfidcut) {
          for (unsigned int j = 0; j < fTrackCuts.size(); j++) {
            if (twoTrackFilter & (uint8_t(1) << j)) {
              if (!fConfigFillLS) {
                dynamic_cast<TH2D*>(fHistRecPair.at(j * fMCSignals.size() + i))->Fill(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt]);
              } else {
                if (uls) {
                  dynamic_cast<TH2D*>(fHistRecPair.at(j * (2 * fMCSignals.size()) + 2 * i))->Fill(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt]);
                } else {
                  dynamic_cast<TH2D*>(fHistRecPair.at(j * (2 * fMCSignals.size()) + 2 * i + 1))->Fill(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt]);
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
          if (genfidcut) {
            for (unsigned int j = 0; j < fTrackCuts.size(); j++) {
              if (twoTrackFilter & (uint8_t(1) << j)) {
                if (!fConfigFillLS) {
                  dynamic_cast<TH2D*>(fHistRecPairMC.at(j * fMCSignals.size() + i))->Fill(mass, pairpt);
                } else {
                  if (uls) {
                    dynamic_cast<TH2D*>(fHistRecPairMC.at(j * (2 * fMCSignals.size()) + 2 * i))->Fill(mass, pairpt);
                  } else {
                    dynamic_cast<TH2D*>(fHistRecPairMC.at(j * (2 * fMCSignals.size()) + 2 * i + 1))->Fill(mass, pairpt);
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
    runPairing<gkEventFillMap, gkMCEventFillMap, gkTrackFillMap>(events, tracks, eventsMC, tracksMC);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  void processDummyNoSkimmed(MyEventsNoSkimmed&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processToEESkimmed, "Run barrel barrel pairing on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy process function", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummyNoSkimmed, "Dummy process function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisEventQa>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
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
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "trigger,cent,mc");
    }

    if (classStr.Contains("Track")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "its,tpcpid,dca,tofpid,mc");
      }
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
    }
    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_barrel");
    }
  }
} // end loop over histogram classes

void SetBinsLinear(std::vector<double>& fBins, const double min, const double max, const unsigned int steps)
{
  fBins.clear();
  const double stepSize = (max - min) / steps;
  for (unsigned int i = 0; i < steps + 1; ++i) {
    fBins.push_back(i * stepSize + min);
  }
}
