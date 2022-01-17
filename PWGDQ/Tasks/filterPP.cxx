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
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "EventFiltering/filterTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include <TH1F.h>
#include <TH2I.h>
#include <THashList.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <memory>
#include <cstring>

using std::cout;
using std::endl;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Some definitions
namespace
{
enum DQTriggers {
  kSingleE = 0,
  kSingleMuLow,
  kSingleMuHigh,
  kDiElectron,
  kDiMuon,
  kNTriggersDQ
};
} // namespace
namespace o2::aod
{
namespace dqppfilter
{
DECLARE_SOA_COLUMN(IsDQEventSelected, isDQEventSelected, int);
DECLARE_SOA_COLUMN(IsDQBarrelSelected, isDQBarrelSelected, uint8_t);
DECLARE_SOA_COLUMN(IsDQMuonSelected, isDQMuonSelected, uint8_t);
} // namespace dqppfilter

DECLARE_SOA_TABLE(DQEventCuts, "AOD", "DQEVENTCUTS", dqppfilter::IsDQEventSelected);
DECLARE_SOA_TABLE(DQBarrelTrackCuts, "AOD", "DQBARRELCUTS", dqppfilter::IsDQBarrelSelected);
DECLARE_SOA_TABLE(DQMuonsCuts, "AOD", "DQMUONCUTS", dqppfilter::IsDQMuonSelected);
} // namespace o2::aod

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsSelected = soa::Join<aod::Collisions, aod::EvSels, aod::DQEventCuts>;
using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyBarrelTracksSelected = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::TrackSelection,
                                         aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                         aod::pidTPCFullKa, aod::pidTPCFullPr,
                                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                         aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                         aod::DQBarrelTrackCuts>;
using MyMuons = aod::FwdTracks;
using MyMuonsSelected = soa::Join<aod::FwdTracks, aod::DQMuonsCuts>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::Muon;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

struct DQEventSelectionTask {
  Produces<aod::DQEventCuts> eventSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan = nullptr;
  AnalysisCompositeCut* fEventCut;

  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Comma separated list of event cuts; multiple cuts are applied with a logical AND"};

  void init(o2::framework::InitContext&)
  {
    DefineCuts();

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;"); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());                 // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  void DefineCuts()
  {
    // default cut is "eventStandard" (kINT7 and vtxZ selection)
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));
    if (!eventCutStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(eventCutStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fEventCut->AddCut(dqcuts::GetAnalysisCut(objArray->At(icut)->GetName()));
      }
    }

    // NOTE: Additional cuts to those specified via the Configurable may still be added

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  template <uint32_t TEventFillMap, typename TEvent>
  void runDQEventSelection(TEvent const& collision, aod::BCs const& bcs)
  {
    // Reset the Values array
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);

    VarManager::FillEvent<TEventFillMap>(collision);
    fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues); // automatically fill all the histograms in the class Event
    if (fEventCut->IsSelected(VarManager::fgValues)) {
      fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
      eventSel(1);
    } else {
      eventSel(0);
    }
  }
  void processDQEventSelection(MyEvents::iterator const& collision, aod::BCs const& bcs)
  {
    runDQEventSelection<gkEventFillMap>(collision, bcs);
  }
  void processDQDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQEventSelectionTask, processDQEventSelection, "Run event selection", false);
  PROCESS_SWITCH(DQEventSelectionTask, processDQDummy, "Dummy function", false);
};

struct DQBarrelTrackSelectionTask {
  Produces<aod::DQBarrelTrackCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fTrackCuts;

  Configurable<std::string> fConfigCuts{"cfgBarrelTrackCuts", "jpsiPID2", "Comma separated list of ADDITIONAL barrel track cuts"};

  void init(o2::framework::InitContext&)
  {
    DefineCuts();

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    TString cutNames = "TrackBarrel_BeforeCuts;";
    for (auto& cut : fTrackCuts) {
      cutNames += Form("TrackBarrel_%s;", cut.GetName());
    }

    DefineHistograms(fHistMan, cutNames.Data());     // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  void DefineCuts()
  {
    // available cuts: jpsiKineAndQuality, jpsiPID1, jpsiPID2
    TString cutNamesStr = "jpsiPID1," + fConfigCuts.value; // "jpsiPID1" is the fixed standard cut provided to the Central Event Filtering task
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // NOTE: Additional cuts to those specified via the Configurable may still be added

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runDQTrackSelection(TEvent const& collision, aod::BCs const& bcs, TTracks const& tracksBarrel)
  {
    uint8_t filterMap = uint8_t(0);
    trackSel.reserve(tracksBarrel.size());

    VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
    // fill event information which might be needed in histograms that combine track and event properties
    VarManager::FillEvent<TEventFillMap>(collision);

    for (auto& track : tracksBarrel) {
      filterMap = uint8_t(0);
      // TODO: if a condition on the event selection is applied, the histogram output is missing
      VarManager::FillTrack<TTrackFillMap>(track);
      fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); ++cut, ++i) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint8_t(1) << i);
          fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
        }
      }
      trackSel(filterMap);
    } // end loop over tracks
  }

  void processDQTrackSelection(MyEvents::iterator const& collision, aod::BCs const& bcs, MyBarrelTracks const& tracks)
  {
    runDQTrackSelection<gkEventFillMap, gkTrackFillMap>(collision, bcs, tracks);
  }
  void processDQDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQBarrelTrackSelectionTask, processDQTrackSelection, "Run barrel track selection", false);
  PROCESS_SWITCH(DQBarrelTrackSelectionTask, processDQDummy, "Dummy function", false);
};

struct DQMuonsSelection {
  Produces<aod::DQMuonsCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fTrackCuts;

  Configurable<std::string> fConfigCuts{"cfgMuonsCuts", "muonQualityCuts", "Comma separated list of ADDITIONAL muon track cuts"};

  void init(o2::framework::InitContext&)
  {
    DefineCuts();

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    TString cutNames = "TrackMuons_BeforeCuts;";
    for (unsigned int i = 0; i < fTrackCuts.size(); i++) {
      cutNames += Form("TrackMuons_%s;", fTrackCuts[i].GetName());
    }

    DefineHistograms(fHistMan, cutNames.Data());     // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  void DefineCuts()
  {
    // available cuts: muonQualityCuts
    TString cutNamesStr = "muonLowPt,muonHighPt," + fConfigCuts.value; // "muonLowPt" & "muonHighPt" are the fixed standard cut provided to the Central Event Filtering task
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // NOTE: Additional cuts to those specified via the Configurable may still be added

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvent, typename TMuons>
  void runDQMuonSelection(TEvent const& collision, aod::BCs const& bcs, TMuons const& muons)
  {
    uint8_t filterMap = uint8_t(0);
    trackSel.reserve(muons.size());

    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
    // fill event information which might be needed in histograms that combine track and event properties
    VarManager::FillEvent<TEventFillMap>(collision);

    for (auto& muon : muons) {
      filterMap = uint8_t(0);
      // TODO: if a condition on the event selection is applied, the histogram output is missing
      VarManager::FillTrack<TMuonFillMap>(muon);
      fHistMan->FillHistClass("TrackMuons_BeforeCuts", VarManager::fgValues);
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); ++cut, ++i) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint8_t(1) << i);
          fHistMan->FillHistClass(Form("TrackMuons_%s", (*cut).GetName()), VarManager::fgValues);
        }
      }
      trackSel(filterMap);
    }
  }

  void processDQMuonSelection(MyEvents::iterator const& collision, aod::BCs const& bcs, MyMuons const& muons)
  {
    runDQMuonSelection<gkEventFillMap, gkMuonFillMap>(collision, bcs, muons);
  }
  void processDQDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQMuonsSelection, processDQMuonSelection, "Run muon selection", false);
  PROCESS_SWITCH(DQMuonsSelection, processDQDummy, "Dummy function", false);
};

struct DQFilterPPTask {
  Produces<aod::DQEventFilter> eventFilter;
  OutputObj<THashList> fOutputList{"output"};
  OutputObj<TH1I> fStatsSingleBarrel{"statsSingleBarrel"};
  OutputObj<TH1I> fStatsSingleMuon{"statsSingleMuon"};
  OutputObj<TH2I> fStatsBarrel{"statsBarrel"};
  OutputObj<TH2I> fStatsMuon{"statsMuon"};
  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fPairCuts;

  Configurable<std::string> fConfigTrackCuts{"cfgBarrelTrackCuts", "jpsiPID2", "Comma separated list of ADDITIONAL barrel track cuts"};
  Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of ADDITIONAL muon track cuts"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "pairMassLow,pairJpsi", "Comma separated list of ADDITIONAL pair cuts"};

  int fNTrackCuts;
  int fNMuonCuts;
  int fNPairCuts;
  TObjArray* fTrkCutsNameArray;
  TObjArray* fMuonCutsNameArray;

  void DefineCuts()
  {
    // available pair cuts in CutsLibrary: pairNoCut,pairMassLow,pairJpsi,pairPsi2S,pairUpsilon,pairJpsiPtLow1, pairJpsiPtLow2
    TString pairCutNamesStr = "pairNoCut," + fConfigPairCuts.value; // "pairNoCut" is the fixed standard cut provided to the Central Event Filtering task
    std::unique_ptr<TObjArray> objArray(pairCutNamesStr.Tokenize(","));
    fNPairCuts = objArray->GetEntries();
    if (fNPairCuts) {
      for (int icut = 0; icut < fNPairCuts; ++icut) {
        fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // NOTE: Additional cuts to those specified via the Configurable may still be added

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  void init(o2::framework::InitContext&)
  {
    DefineCuts();

    // initialize the variable manager
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // configure histograms
    TString trackCutNamesStr = "jpsiPID1," + fConfigTrackCuts.value;
    fTrkCutsNameArray = trackCutNamesStr.Tokenize(",");
    fNTrackCuts = fTrkCutsNameArray->GetEntries();
    TString muonCutNamesStr = "muonLowPt,muonHighPt," + fConfigMuonCuts.value;
    fMuonCutsNameArray = muonCutNamesStr.Tokenize(",");
    fNMuonCuts = fMuonCutsNameArray->GetEntries();
    TString histNames = "";
    for (int i = 0; i < fNTrackCuts; i++) {
      histNames += Form("TrackBarrel_%s;PairsBarrelPM_%s;", fTrkCutsNameArray->At(i)->GetName(), fTrkCutsNameArray->At(i)->GetName());
    }
    for (int i = 0; i < fNMuonCuts; i++) {
      histNames += Form("TrackMuons_%s;PairsMuonsPM_%s;PairsMuonsPP_%s;PairsMuonsMM_%s;", fMuonCutsNameArray->At(i)->GetName(), fMuonCutsNameArray->At(i)->GetName(), fMuonCutsNameArray->At(i)->GetName(), fMuonCutsNameArray->At(i)->GetName());
    }
    DefineHistograms(fHistMan, histNames.Data());    // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    // configure the stats histogram
    fStatsSingleBarrel.setObject(new TH1I("statsSingleBarrel", "statsSingleBarrel", fNTrackCuts, -0.5, -0.5 + fNTrackCuts));
    for (int i = 0; i < fNTrackCuts; ++i) {
      fStatsSingleBarrel->GetXaxis()->SetBinLabel(i + 1, fTrkCutsNameArray->At(i)->GetName());
    }

    fStatsSingleMuon.setObject(new TH1I("statsSingleMuon", "statsSingleMuon", fNMuonCuts, -0.5, -0.5 + fNMuonCuts));
    for (int i = 0; i < fNMuonCuts; ++i) {
      fStatsSingleMuon->GetXaxis()->SetBinLabel(i + 1, fMuonCutsNameArray->At(i)->GetName());
    }

    fStatsBarrel.setObject(new TH2I("statsBarrel", "statsBarrel", fNPairCuts, -0.5, -0.5 + fNPairCuts, fNTrackCuts, -0.5, -0.5 + fNTrackCuts));
    for (int i = 0; i < fNPairCuts; ++i) {
      fStatsBarrel->GetXaxis()->SetBinLabel(i + 1, fPairCuts[i].GetName());
    }
    for (int i = 0; i < fNTrackCuts; ++i) {
      fStatsBarrel->GetYaxis()->SetBinLabel(i + 1, fTrkCutsNameArray->At(i)->GetName());
    }

    fStatsMuon.setObject(new TH2I("statsMuon", "statsMuon", fNPairCuts, -0.5, -0.5 + fNPairCuts, fNMuonCuts, -0.5, -0.5 + fNMuonCuts));
    for (int i = 0; i < fNPairCuts; ++i) {
      fStatsMuon->GetXaxis()->SetBinLabel(i + 1, fPairCuts[i].GetName());
    }
    for (int i = 0; i < fNMuonCuts; ++i) {
      fStatsMuon->GetYaxis()->SetBinLabel(i + 1, fMuonCutsNameArray->At(i)->GetName());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, typename TEvent, typename TTracks, typename TMuons>
  void runDQFilterPP(TEvent const& collision, aod::BCs const& bcs, TTracks const& tracksBarrel, TMuons const& muons)
  {
    uint64_t filter = 0; // first
    constexpr int pairTypeEE = VarManager::kJpsiToEE;
    constexpr int pairTypeMuMu = VarManager::kJpsiToMuMu;

    if (!collision.isDQEventSelected()) {
      return;
    }
    // Reset the Values array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<TEventFillMap>(collision);
    uint8_t* singleBarrelCount = new uint8_t[fNTrackCuts];                 // single electron trigger
    uint8_t* singleMuonsCount = new uint8_t[fNMuonCuts];                   // single muon trigger
    uint8_t* pairsBarrelCount = new uint8_t[fNPairCuts * fNTrackCuts];     // dielectron trigger
    uint8_t* pairsUnlikeMuonsCount = new uint8_t[fNPairCuts * fNMuonCuts]; // dimuon trigger Unlike
    uint8_t* pairsLikeMuonsCount = new uint8_t[fNPairCuts * fNMuonCuts];   // dimuon trigger like

    for (int i = 0; i < fNPairCuts * fNTrackCuts; i++) {
      pairsBarrelCount[i] = 0;
    }
    for (int i = 0; i < fNPairCuts * fNMuonCuts; i++) {
      pairsUnlikeMuonsCount[i] = 0;
      pairsLikeMuonsCount[i] = 0;
    }
    for (int i = 0; i < fNMuonCuts; i++) {
      singleMuonsCount[i] = 0;
    }
    for (int i = 0; i < fNTrackCuts; i++) {
      singleBarrelCount[i] = 0;
    }

    uint8_t cutSingleFilter = 0;
    for (auto track : tracksBarrel) {
      cutSingleFilter = track.isDQBarrelSelected();
      if (!cutSingleFilter) {
        continue;
      }
      for (int i = 0; i < fNTrackCuts; ++i) {
        if (!(cutSingleFilter & (uint8_t(1) << i))) {
          continue;
        }
        VarManager::FillTrack<TTrackFillMap>(track);
        fHistMan->FillHistClass(Form("TrackBarrel_%s", fTrkCutsNameArray->At(i)->GetName()), VarManager::fgValues);
        singleBarrelCount[i] += 1;
      }
    }

    uint8_t cutFilter = 0;
    for (auto& [t1, t2] : combinations(tracksBarrel, tracksBarrel)) {
      //        for (auto tneg : negTracks) { // +- pairs
      cutFilter = t1.isDQBarrelSelected() & t2.isDQBarrelSelected();
      if (!cutFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }
      VarManager::FillPair<pairTypeEE, TTrackFillMap>(t1, t2); // compute pair quantities
      if (t1.sign() * t2.sign() < 0) {
        for (int i = 0; i < fNTrackCuts; ++i) {
          if (!(cutFilter & (uint8_t(1) << i))) {
            continue;
          }
          fHistMan->FillHistClass(Form("PairsBarrelPM_%s", fTrkCutsNameArray->At(i)->GetName()), VarManager::fgValues);
          int j = 0;
          for (auto cut = fPairCuts.begin(); cut != fPairCuts.end(); cut++, j++) {
            if ((*cut).IsSelected(VarManager::fgValues)) {
              pairsBarrelCount[j + i * fNPairCuts] += 1;
            }
          }
        }
      }
    }

    uint8_t cutSingleMuonFilter = 0;
    for (auto muon1 : muons) {
      cutSingleMuonFilter = muon1.isDQMuonSelected();
      if (!cutSingleMuonFilter) {
        continue;
      }
      for (int i = 0; i < fNMuonCuts; ++i) {
        if (!(cutSingleMuonFilter & (uint8_t(1) << i))) {
          continue;
        }
        VarManager::FillTrack<TMuonFillMap>(muon1);
        fHistMan->FillHistClass(Form("TrackMuons_%s", fMuonCutsNameArray->At(i)->GetName()), VarManager::fgValues);
        singleMuonsCount[i] += 1;
      }
    }

    uint8_t cutMuonFilter = 0;
    for (auto& [muon1, muon2] : combinations(muons, muons)) {
      cutMuonFilter = muon1.isDQMuonSelected() & muon2.isDQMuonSelected();
      if (!cutMuonFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }
      VarManager::FillPair<pairTypeMuMu, TMuonFillMap>(muon1, muon2); // compute pair quantities
      for (int i = 0; i < fNMuonCuts; ++i) {
        if (!(cutMuonFilter & (uint8_t(1) << i))) {
          continue;
        }
        if (muon1.sign() * muon2.sign() < 0) {
          fHistMan->FillHistClass(Form("PairsMuonsPM_%s", fMuonCutsNameArray->At(i)->GetName()), VarManager::fgValues);
        } else {
          if (muon1.sign() > 0) {
            fHistMan->FillHistClass(Form("PairsMuonsPP_%s", fMuonCutsNameArray->At(i)->GetName()), VarManager::fgValues);
          } else {
            fHistMan->FillHistClass(Form("PairsMuonsMM_%s", fMuonCutsNameArray->At(i)->GetName()), VarManager::fgValues);
          }
        }

        int j = 0;
        for (auto cut = fPairCuts.begin(); cut != fPairCuts.end(); cut++, j++) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            if (muon1.sign() * muon2.sign() < 0) {
              pairsUnlikeMuonsCount[j + i * fNPairCuts] += 1;
            } else {
              pairsLikeMuonsCount[j + i * fNPairCuts] += 1;
            }
          }
        }
      }
    }

    // Fill DQ bit map
    for (int i = 0; i < fNTrackCuts; i++) {
      if (singleBarrelCount[i] > 0) {
        filter |= (uint64_t(1) << i);
        fStatsSingleBarrel->Fill(i);
      }
      for (int j = 0; j < fNPairCuts; j++) {
        if (pairsBarrelCount[j + i * fNPairCuts] > 0) {
          filter |= (uint64_t(1) << (8 + j + i * fNPairCuts));
          fStatsBarrel->Fill(j, i);
        }
      }
    }
    for (int i = 0; i < fNMuonCuts; i++) {
      if (singleMuonsCount[i] > 0) {
        filter |= (uint64_t(1) << (32 + i));
        fStatsSingleMuon->Fill(i);
      }
      for (int j = 0; j < fNPairCuts; j++) {
        if (pairsUnlikeMuonsCount[j + i * fNPairCuts] > 0) {
          filter |= (uint64_t(1) << (40 + j + i * fNPairCuts));
          fStatsMuon->Fill(j, i);
        }
        if (pairsLikeMuonsCount[j + i * fNPairCuts] > 0) {
          filter |= (uint64_t(1) << (52 + j + i * fNPairCuts));
          //            fStatsMuon->Fill(j, i);
        }
      }
    }

    delete[] pairsBarrelCount;
    delete[] pairsUnlikeMuonsCount;
    delete[] pairsLikeMuonsCount;
    delete[] singleMuonsCount;
    delete[] singleBarrelCount;

    eventFilter(filter);
  }
  void processDQFilterPP(MyEventsSelected::iterator const& collision, aod::BCs const& bcs, MyBarrelTracksSelected const& tracks, MyMuonsSelected const& muons)
  {
    runDQFilterPP<gkEventFillMap, gkTrackFillMap, gkMuonFillMap>(collision, bcs, tracks, muons);
  }
  // TODO: dummy function for the case when no process function is enabled
  void processDQDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQFilterPPTask, processDQFilterPP, "Run filter task", false);
  PROCESS_SWITCH(DQFilterPPTask, processDQDummy, "Dummy function", false);
};

struct DQCentralFilterPPTask {
  Produces<aod::DqFilters> dqtable;
  OutputObj<THashList> fOutputList{"output"};
  OutputObj<TH1I> fStatsTriggers{"statsTriggers"};
  HistogramManager* fHistMan;

  int fNTrackCuts;
  int fNMuonCuts;
  int fNPairCuts;
  TObjArray* fTrkCutsNameArray;
  TObjArray* fMuonCutsNameArray;

  void init(o2::framework::InitContext&)
  {

    // initialize the variable manager
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // configure histograms
    TString trackCutNamesStr = "jpsiPID1";
    fTrkCutsNameArray = trackCutNamesStr.Tokenize(",");
    fNTrackCuts = fTrkCutsNameArray->GetEntries();
    TString muonCutNamesStr = "muonLowPt,muonHighPt";
    fMuonCutsNameArray = muonCutNamesStr.Tokenize(",");
    fNMuonCuts = fMuonCutsNameArray->GetEntries();

    // configure the stats histogram
    fStatsTriggers.setObject(new TH1I("statsTriggers", "statsTriggers", kNTriggersDQ, -0.5, -0.5 + kNTriggersDQ));
    fStatsTriggers->GetXaxis()->SetBinLabel(1, "Single Electron");
    fStatsTriggers->GetXaxis()->SetBinLabel(2, "Single Low pT Muon");
    fStatsTriggers->GetXaxis()->SetBinLabel(3, "Single High pT Muon");
    fStatsTriggers->GetXaxis()->SetBinLabel(4, "Di-Electron");
    fStatsTriggers->GetXaxis()->SetBinLabel(5, "Di-Muon");
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, typename TEvent, typename TTracks, typename TMuons>
  void runCentralDQFilterPP(TEvent const& collision, aod::BCs const& bcs, TTracks const& tracksBarrel, TMuons const& muons)
  {
    bool keepEvent[kNTriggersDQ]{false};

    if (!collision.isDQEventSelected()) {
      return;
    }
    // Reset the Values array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<TEventFillMap>(collision);

    uint8_t cutSingleFilter = 0;
    for (auto track : tracksBarrel) {
      cutSingleFilter = track.isDQBarrelSelected();
      if (!cutSingleFilter) {
        continue;
      }
      for (int i = 0; i < 1; ++i) { // Single Barrel count : bit 0 corresponds to "jpsiPID1" cut
        if (!(cutSingleFilter & (uint8_t(1) << i))) {
          continue;
        }
        keepEvent[0] = true; // Single Electron
      }
    }

    uint8_t cutFilter = 0;
    for (auto& [t1, t2] : combinations(tracksBarrel, tracksBarrel)) {
      //        for (auto tneg : negTracks) { // +- pairs
      cutFilter = t1.isDQBarrelSelected() & t2.isDQBarrelSelected();
      if (!cutFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }
      if (t1.sign() * t2.sign() < 0) {
        for (int i = 0; i < 1; ++i) {
          if (!(cutFilter & (uint8_t(1) << i))) {
            continue;
          }
          keepEvent[3] = true; // Di-Electron
        }
      }
    }

    uint8_t cutSingleMuonFilter = 0;
    for (auto muon1 : muons) {
      cutSingleMuonFilter = muon1.isDQMuonSelected();
      if (!cutSingleMuonFilter) {
        continue;
      }
      for (int i = 0; i < 2; ++i) {
        if (!(cutSingleMuonFilter & (uint8_t(1) << i))) {
          continue;
        }
        keepEvent[i + 1] = true; // Low pT muon (keepEvent[1]) ang High pT muons (keepEvent[2])
      }
    }

    uint8_t cutMuonFilter = 0;
    for (auto& [muon1, muon2] : combinations(muons, muons)) {
      cutMuonFilter = muon1.isDQMuonSelected() & muon2.isDQMuonSelected();
      if (!cutMuonFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }
      for (int i = 0; i < 1; ++i) {
        if (!(cutMuonFilter & (uint8_t(1) << i))) {
          continue;
        }
        if (muon1.sign() * muon2.sign() < 0) {
          keepEvent[4] = true; // Di-Muon (low pT muons, no pair cut)
        }
      }
    }
    for (int i = 0; i < kNTriggersDQ; ++i) {
      if (keepEvent[i] > 0) { // the tracks must have at least one filter bit in common to continue
        fStatsTriggers->Fill(i);
      }
    }
    // Filling the table
    dqtable(keepEvent[kSingleE], keepEvent[kSingleMuLow], keepEvent[kSingleMuHigh], keepEvent[kDiElectron], keepEvent[kDiMuon]);
  }
  void processCentralDQFilterPP(MyEventsSelected::iterator const& collision, aod::BCs const& bcs, MyBarrelTracksSelected const& tracks, MyMuonsSelected const& muons)
  {
    runCentralDQFilterPP<gkEventFillMap, gkTrackFillMap, gkMuonFillMap>(collision, bcs, tracks, muons);
  }
  // TODO: dummy function for the case when no process function is enabled
  void processDQDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQCentralFilterPPTask, processCentralDQFilterPP, "Run central filter task", false);
  PROCESS_SWITCH(DQCentralFilterPPTask, processDQDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DQEventSelectionTask>(cfgc),
    adaptAnalysisTask<DQBarrelTrackSelectionTask>(cfgc),
    adaptAnalysisTask<DQMuonsSelection>(cfgc),
    adaptAnalysisTask<DQFilterPPTask>(cfgc),
    adaptAnalysisTask<DQCentralFilterPPTask>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //  The histogram classes are provided in the histClasses string, separated by semicolon ";"
  //  NOTE: Additional histograms to those predefined in the library can be added here !!
  const int kNRuns = 25;
  int runs[kNRuns] = {
    285009, 285011, 285012, 285013, 285014, 285015, 285064, 285065, 285066, 285106,
    285108, 285125, 285127, 285165, 285200, 285202, 285203, 285222, 285224, 285327,
    285328, 285347, 285364, 285365, 285396};
  VarManager::SetRunNumbers(kNRuns, runs);

  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "trigger,vtxPbPb");
    }

    if (classStr.Contains("Track")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "its,tpcpid,dca");
    }

    if (classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_barrel", "vertexing-barrel");
      }
      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_dimuon");
      }
    }
  }
}
