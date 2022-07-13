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
#include "Common/DataModel/PIDResponse.h"
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
  kDiElectron,
  kSingleMuLow,
  kSingleMuHigh,
  kDiMuon,
  kNTriggersDQ
};
} // namespace
namespace o2::aod
{
namespace dqppfilter
{
DECLARE_SOA_COLUMN(IsDQEventSelected, isDQEventSelected, int);
DECLARE_SOA_COLUMN(IsDQBarrelSelected, isDQBarrelSelected, uint32_t);
DECLARE_SOA_COLUMN(IsDQMuonSelected, isDQMuonSelected, uint32_t);
} // namespace dqppfilter

DECLARE_SOA_TABLE(DQEventCuts, "AOD", "DQEVENTCUTS", dqppfilter::IsDQEventSelected);
DECLARE_SOA_TABLE(DQBarrelTrackCuts, "AOD", "DQBARRELCUTS", dqppfilter::IsDQBarrelSelected);
DECLARE_SOA_TABLE(DQMuonsCuts, "AOD", "DQMUONCUTS", dqppfilter::IsDQMuonSelected);
} // namespace o2::aod

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsSelected = soa::Join<aod::Collisions, aod::EvSels, aod::DQEventCuts>;
// TODO: subscribe to the bare needed minimum, in particular for the CEFP task
using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyBarrelTracksTiny = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                     aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
                                     aod::pidTPCKa, aod::pidTPCPr,
                                     aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi,
                                     aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFbeta>;

using MyBarrelTracksSelected = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                         aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                         aod::pidTPCFullKa, aod::pidTPCFullPr,
                                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                         aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                         aod::DQBarrelTrackCuts>;
using MyBarrelTracksSelectedTiny = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                             aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
                                             aod::pidTPCKa, aod::pidTPCPr,
                                             aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi,
                                             aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFbeta,
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
  Configurable<bool> fConfigQA{"cfgWithQA", false, "If true, fill QA histograms"};
  // TODO: configure the histogram classes to be filled by QA

  void init(o2::framework::InitContext&)
  {
    // Construct the composite cut out of the user provided event selection cut(s)
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut = new AnalysisCompositeCut(true);
    if (!eventCutStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(eventCutStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fEventCut->AddCut(dqcuts::GetAnalysisCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;"); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                 // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  void processEventSelection(MyEvents::iterator const& collision, aod::BCs const& bcs)
  {
    // Reset the Values array
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);

    VarManager::FillEvent<gkEventFillMap>(collision);
    if (fConfigQA) {
      fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
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

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQEventSelectionTask, processEventSelection, "Run event selection", false);
  PROCESS_SWITCH(DQEventSelectionTask, processDummy, "Dummy function", false);
};

struct DQBarrelTrackSelection {
  Produces<aod::DQBarrelTrackCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;

  Configurable<std::string> fConfigCuts{"cfgBarrelTrackCuts", "jpsiPID1", "Comma separated list of ADDITIONAL barrel track cuts"};
  Configurable<bool> fConfigQA{"cfgWithQA", false, "If true, fill QA histograms"};
  // TODO: configure the histogram classes to be filled by QA

  std::vector<AnalysisCompositeCut> fTrackCuts;
  std::vector<TString> fCutHistNames;

  void init(o2::framework::InitContext&)
  {
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        AnalysisCompositeCut* cut = dqcuts::GetCompositeCut(objArray->At(icut)->GetName());
        if (cut) {
          fTrackCuts.push_back(*cut);
        } else {
          LOGF(fatal, "Invalid barrel track cut provided: %s", objArray->At(icut)->GetName());
        }
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      TString cutNames = "TrackBarrel_BeforeCuts;";
      for (auto& cut : fTrackCuts) {
        cutNames += Form("TrackBarrel_%s;", cut.GetName());
        fCutHistNames.push_back(Form("TrackBarrel_%s", cut.GetName()));
      }

      DefineHistograms(fHistMan, cutNames.Data());     // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runTrackSelection(TEvent const& collisions, aod::BCs const& bcs, TTracks const& tracksBarrel)
  {
    uint32_t filterMap = uint32_t(0);
    trackSel.reserve(tracksBarrel.size());
    int CollisionId = -1;

    VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);

    for (auto& track : tracksBarrel) {
      filterMap = uint32_t(0);
      if (!track.has_collision()) {
        trackSel(uint32_t(0));
      } else {
        // fill event information which might be needed in histograms or cuts that combine track and event properties
        if (track.collisionId() != CollisionId) { // check if the track belongs to a different event than the previous one
          CollisionId = track.collisionId();
          auto collision = track.template collision_as<TEvent>();
          VarManager::FillEvent<TEventFillMap>(collision);
        }
        VarManager::FillTrack<TTrackFillMap>(track);
        if (fConfigQA) {
          fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
        }
        int i = 0;
        for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); ++cut, ++i) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            filterMap |= (uint32_t(1) << i);
            if (fConfigQA) {
              fHistMan->FillHistClass(fCutHistNames[i].Data(), VarManager::fgValues);
            }
          }
        }
        trackSel(filterMap);
      }
    } // end loop over tracks
  }

  void processSelection(MyEvents const& collisions, aod::BCs const& bcs, MyBarrelTracks const& tracks)
  {
    runTrackSelection<gkEventFillMap, gkTrackFillMap>(collisions, bcs, tracks);
  }
  void processSelectionTiny(MyEvents const& collisions, aod::BCs const& bcs, MyBarrelTracksTiny const& tracks)
  {
    runTrackSelection<gkEventFillMap, gkTrackFillMap>(collisions, bcs, tracks);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQBarrelTrackSelection, processSelection, "Run barrel track selection", false);
  PROCESS_SWITCH(DQBarrelTrackSelection, processSelectionTiny, "Run barrel track selection", false);
  PROCESS_SWITCH(DQBarrelTrackSelection, processDummy, "Dummy function", false);
};

struct DQMuonsSelection {
  Produces<aod::DQMuonsCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;

  Configurable<std::string> fConfigCuts{"cfgMuonsCuts", "muonQualityCuts", "Comma separated list of ADDITIONAL muon track cuts"};
  Configurable<bool> fConfigQA{"cfgWithQA", false, "If true, fill QA histograms"};
  // TODO: configure the histogram classes to be filled by QA

  std::vector<AnalysisCompositeCut> fTrackCuts;
  std::vector<TString> fCutHistNames;

  void init(o2::framework::InitContext&)
  {
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars);

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      TString cutNames = "Muon_BeforeCuts;";
      for (unsigned int i = 0; i < fTrackCuts.size(); i++) {
        cutNames += Form("Muon_%s;", fTrackCuts[i].GetName());
        fCutHistNames.push_back(Form("Muon_%s", fTrackCuts[i].GetName()));
      }

      DefineHistograms(fHistMan, cutNames.Data());     // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvent, typename TMuons>
  void runMuonSelection(TEvent const& collisions, aod::BCs const& bcs, TMuons const& muons)
  {
    uint32_t filterMap = uint32_t(0);
    trackSel.reserve(muons.size());
    int CollisionId = -1;

    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);

    for (auto& muon : muons) {
      filterMap = uint32_t(0);
      if (!muon.has_collision()) {
        trackSel(uint32_t(0));
      } else {
        // fill event information which might be needed in histograms or cuts that combine track and event properties
        if (muon.collisionId() != CollisionId) { // check if the track belongs to a different event than the previous one
          CollisionId = muon.collisionId();
          auto collision = muon.template collision_as<TEvent>();
          VarManager::FillEvent<TEventFillMap>(collision);
        }
        VarManager::FillTrack<TMuonFillMap>(muon);
        if (fConfigQA) {
          fHistMan->FillHistClass("Muon_BeforeCuts", VarManager::fgValues);
        }
        int i = 0;
        for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); ++cut, ++i) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            filterMap |= (uint32_t(1) << i);
            if (fConfigQA) {
              fHistMan->FillHistClass(fCutHistNames[i].Data(), VarManager::fgValues);
            }
          }
        }
        trackSel(filterMap);
      }
    } // end loop over tracks
  }

  void processSelection(MyEvents const& collisions, aod::BCs const& bcs, MyMuons const& muons)
  {
    runMuonSelection<gkEventFillMap, gkMuonFillMap>(collisions, bcs, muons);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQMuonsSelection, processSelection, "Run muon selection", false);
  PROCESS_SWITCH(DQMuonsSelection, processDummy, "Dummy function", false);
};

struct DQFilterPPTask {
  Produces<aod::DQEventFilter> eventFilter;
  Produces<aod::DqFilters> dqtable;
  OutputObj<THashList> fOutputList{"output"};
  OutputObj<TH1I> fStats{"Statistics"};
  HistogramManager* fHistMan;

  Configurable<std::string> fConfigBarrelSelections{"cfgBarrelSels", "jpsiPID1:pairMassLow:1", "<track-cut>:[<pair-cut>]:<n>,[<track-cut>:[<pair-cut>]:<n>],..."};
  Configurable<std::string> fConfigMuonSelections{"cfgMuonSels", "muonQualityCuts:pairNoCut:1", "<muon-cut>:[<pair-cut>]:<n>"};
  Configurable<bool> fConfigQA{"cfgWithQA", false, "If true, fill QA histograms"};

  Filter filterBarrelTrackSelected = aod::dqppfilter::isDQBarrelSelected > uint32_t(0);
  Filter filterMuonTrackSelected = aod::dqppfilter::isDQMuonSelected > uint32_t(0);

  int fNBarrelCuts;                                    // number of barrel selections
  int fNMuonCuts;                                      // number of muon selections
  std::vector<bool> fBarrelRunPairing;                 // bit map on whether the selections require pairing (barrel)
  std::vector<bool> fMuonRunPairing;                   // bit map on whether the selections require pairing (muon)
  std::vector<int> fBarrelNreqObjs;                    // minimal number of tracks/pairs required (barrel)
  std::vector<int> fMuonNreqObjs;                      // minimal number of tracks/pairs required (muon)
  std::map<int, AnalysisCompositeCut> fBarrelPairCuts; // map of barrel pair cuts
  std::map<int, AnalysisCompositeCut> fMuonPairCuts;   // map of muon pair cuts
  std::map<int, TString> fBarrelPairHistNames;         // map with names of the barrel pairing histogram directories
  std::map<int, TString> fMuonPairHistNames;           // map with names of the muon pairing histogram directories

  void DefineCuts()
  {
    TString barrelSelsStr = fConfigBarrelSelections.value;
    std::unique_ptr<TObjArray> objArray(barrelSelsStr.Tokenize(","));
    fNBarrelCuts = objArray->GetEntries();
    if (fNBarrelCuts) {
      for (int icut = 0; icut < fNBarrelCuts; ++icut) {
        TString selStr = objArray->At(icut)->GetName();
        std::unique_ptr<TObjArray> sel(selStr.Tokenize(":"));
        if (sel->GetEntries() < 2 || sel->GetEntries() > 3) {
          continue;
        }
        // if "sel" contains 3 entries, it means the user provided both the track and pair cuts
        if (sel->GetEntries() == 3) {
          fBarrelPairCuts[icut] = (*dqcuts::GetCompositeCut(sel->At(1)->GetName()));
          fBarrelRunPairing.push_back(true);
          fBarrelNreqObjs.push_back(std::atoi(sel->At(2)->GetName()));
          fBarrelPairHistNames[icut] = Form("PairsBarrelSEPM_%s_%s", sel->At(0)->GetName(), sel->At(1)->GetName());
        } else {
          fBarrelNreqObjs.push_back(std::atoi(sel->At(1)->GetName()));
          fBarrelRunPairing.push_back(false);
        }
      }
    }
    TString muonSelsStr = fConfigMuonSelections.value;
    std::unique_ptr<TObjArray> objArray2(muonSelsStr.Tokenize(","));
    fNMuonCuts = objArray2->GetEntries();
    if (fNMuonCuts) {
      for (int icut = 0; icut < fNMuonCuts; ++icut) {
        TString selStr = objArray2->At(icut)->GetName();
        std::unique_ptr<TObjArray> sel(selStr.Tokenize(":"));
        if (sel->GetEntries() < 2 || sel->GetEntries() > 3) {
          continue;
        }
        // if "sel" contains 3 entries, it means the user provided both the track and pair cuts
        if (sel->GetEntries() == 3) {
          fMuonPairCuts[icut] = (*dqcuts::GetCompositeCut(sel->At(1)->GetName()));
          fMuonRunPairing.push_back(true);
          fMuonNreqObjs.push_back(std::atoi(sel->At(2)->GetName()));
          fMuonPairHistNames[icut] = Form("PairsMuonSEPM_%s_%s", sel->At(0)->GetName(), sel->At(1)->GetName());
        } else {
          fMuonNreqObjs.push_back(std::atoi(sel->At(1)->GetName()));
          fMuonRunPairing.push_back(false);
        }
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars);

    // setup the Stats histogram
    fStats.setObject(new TH1I("Statistics", "Stats for DQ triggers", fNBarrelCuts + fNMuonCuts + 2, -2.5, -0.5 + fNBarrelCuts + fNMuonCuts));
    fStats->GetXaxis()->SetBinLabel(1, "Events inspected");
    fStats->GetXaxis()->SetBinLabel(2, "Events selected");
    if (fNBarrelCuts) {
      for (int ib = 3; ib < 3 + fNBarrelCuts; ib++) {
        fStats->GetXaxis()->SetBinLabel(ib, objArray->At(ib - 3)->GetName());
      }
    }
    if (fNMuonCuts) {
      for (int ib = 3 + fNBarrelCuts; ib < 3 + fNBarrelCuts + fNMuonCuts; ib++) {
        fStats->GetXaxis()->SetBinLabel(ib, objArray2->At(ib - 3 - fNBarrelCuts)->GetName());
      }
    }
  }

  void init(o2::framework::InitContext&)
  {
    DefineCuts();

    if (fConfigQA) {
      // initialize the variable manager
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      TString histNames = "";
      for (const auto& [key, value] : fBarrelPairHistNames) {
        histNames += value;
        histNames += ";";
      }
      for (const auto& [key, value] : fMuonPairHistNames) {
        histNames += value;
        histNames += ";";
      }
      DefineHistograms(fHistMan, histNames.Data());
      VarManager::SetUseVars(fHistMan->GetUsedVars());
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, typename TEvent, typename TTracks, typename TMuons>
  void runFilterPP(TEvent const& collision, aod::BCs const& bcs, TTracks const& tracksBarrel, TMuons const& muons)
  {
    fStats->Fill(-2.0);
    // if the event is not selected produce tables and return
    if (!collision.isDQEventSelected()) {
      eventFilter(0);
      dqtable(false, false, false, false, false);
      return;
    }
    fStats->Fill(-1.0);

    if (tracksBarrel.size() == 0 && muons.size() == 0) {
      eventFilter(0);
      dqtable(false, false, false, false, false);
      return;
    }

    // Reset the values array and compute event quantities
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<TEventFillMap>(collision);

    std::vector<int> objCountersBarrel(fNBarrelCuts, 0); // init all counters to zero
    // count the number of barrel tracks fulfilling each cut
    for (auto track : tracksBarrel) {
      for (int i = 0; i < fNBarrelCuts; ++i) {
        if (track.isDQBarrelSelected() & (uint32_t(1) << i)) {
          objCountersBarrel[i] += 1;
        }
      }
    }

    // check which selections require pairing
    uint32_t pairingMask = 0; // in order to know which of the selections actually require pairing
    for (int i = 0; i < fNBarrelCuts; i++) {
      if (fBarrelRunPairing[i]) {
        if (objCountersBarrel[i] > 1) { // pairing has to be enabled and at least two tracks are needed
          pairingMask |= (uint32_t(1) << i);
        }
        objCountersBarrel[i] = 0; // reset counters for selections where pairing is needed (count pairs instead)
      }
    }

    // run pairing if there is at least one selection that requires it
    uint32_t pairFilter = 0;
    if (pairingMask > 0) {
      for (auto& [t1, t2] : combinations(tracksBarrel, tracksBarrel)) {
        // keep just opposite-sign pairs
        if (t1.sign() * t2.sign() > 0) {
          continue;
        }
        // check the pairing mask and that the tracks share a cut bit
        pairFilter = pairingMask & t1.isDQBarrelSelected() & t2.isDQBarrelSelected();
        if (pairFilter == 0) {
          continue;
        }
        // construct the pair and apply pair cuts
        VarManager::FillPair<VarManager::kJpsiToEE, TTrackFillMap>(t1, t2); // compute pair quantities
        for (int icut = 0; icut < fNBarrelCuts; icut++) {
          if (!(pairFilter & (uint32_t(1) << icut))) {
            continue;
          }
          if (!fBarrelPairCuts[icut].IsSelected(VarManager::fgValues)) {
            continue;
          }
          objCountersBarrel[icut] += 1; // count the pair
          if (fConfigQA) {              // fill histograms if QA is enabled
            fHistMan->FillHistClass(fBarrelPairHistNames[icut].Data(), VarManager::fgValues);
          }
        }
      }
    }

    std::vector<int> objCountersMuon(fNMuonCuts, 0); // init all counters to zero
    // count the number of muon tracks fulfilling each selection
    for (auto muon : muons) {
      for (int i = 0; i < fNMuonCuts; ++i) {
        if (muon.isDQMuonSelected() & (uint32_t(1) << i)) {
          objCountersMuon[i] += 1;
        }
      }
    }

    // check which muon selections require pairing
    pairingMask = 0; // reset the mask for the muons
    for (int i = 0; i < fNMuonCuts; i++) {
      if (fMuonRunPairing[i]) { // pairing has to be enabled and at least two tracks are needed
        if (objCountersMuon[i] > 1) {
          pairingMask |= (uint32_t(1) << i);
        }
        objCountersMuon[i] = 0; // reset counters for selections where pairing is needed (count pairs instead)
      }
    }

    // run pairing if there is at least one selection that requires it
    pairFilter = 0;
    if (pairingMask > 0) {
      for (auto& [t1, t2] : combinations(muons, muons)) {
        // keep just opposite-sign pairs
        if (t1.sign() * t2.sign() > 0) {
          continue;
        }
        // check the pairing mask and that the tracks share a cut bit
        pairFilter = pairingMask & t1.isDQMuonSelected() & t2.isDQMuonSelected();
        if (pairFilter == 0) {
          continue;
        }
        // construct the pair and apply cuts
        VarManager::FillPair<VarManager::kJpsiToMuMu, TTrackFillMap>(t1, t2); // compute pair quantities
        for (int icut = 0; icut < fNMuonCuts; icut++) {
          if (!(pairFilter & (uint32_t(1) << icut))) {
            continue;
          }
          if (!fMuonPairCuts[icut].IsSelected(VarManager::fgValues)) {
            continue;
          }
          objCountersMuon[icut] += 1;
          if (fConfigQA) {
            fHistMan->FillHistClass(fMuonPairHistNames[icut].Data(), VarManager::fgValues);
          }
        }
      }
    }

    // compute the decisions and publish
    // NOTE: For the CEFP decisions, decisions are placed in a vector of bool in an ordered way:
    //       start with all configured barrel selections and then continue with those from muons
    //       The configured order has to be in sync with that implemented in the cefp task and can be done
    //       by preparing a dedicated json configuration file
    std::vector<bool> decisions(kNTriggersDQ, false); // event decisions to be transmited to CEFP
    uint64_t filter = 0;
    for (int i = 0; i < fNBarrelCuts; i++) {
      if (objCountersBarrel[i] >= fBarrelNreqObjs[i]) {
        filter |= (uint64_t(1) << i);
        fStats->Fill(float(i));
        if (i < kNTriggersDQ) {
          decisions[i] = true;
        }
      }
    }
    for (int i = 0; i < fNMuonCuts; i++) {
      if (objCountersMuon[i] >= fMuonNreqObjs[i]) {
        filter |= (uint64_t(1) << (i + fNBarrelCuts));
        fStats->Fill(float(i + fNBarrelCuts));
        if (i + fNBarrelCuts < kNTriggersDQ) {
          decisions[i + fNBarrelCuts] = true;
        }
      }
    }
    eventFilter(filter);
    dqtable(decisions[0], decisions[1], decisions[2], decisions[3], decisions[4]);
  }

  void processFilterPP(MyEventsSelected::iterator const& collision, aod::BCs const& bcs,
                       soa::Filtered<MyBarrelTracksSelected> const& tracks, soa::Filtered<MyMuonsSelected> const& muons)
  {
    runFilterPP<gkEventFillMap, gkTrackFillMap, gkMuonFillMap>(collision, bcs, tracks, muons);
  }

  void processFilterPPTiny(MyEventsSelected::iterator const& collision, aod::BCs const& bcs,
                           soa::Filtered<MyBarrelTracksSelectedTiny> const& tracks, soa::Filtered<MyMuonsSelected> const& muons)
  {
    runFilterPP<gkEventFillMap, gkTrackFillMap, gkMuonFillMap>(collision, bcs, tracks, muons);
  }

  // TODO: dummy function for the case when no process function is enabled
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQFilterPPTask, processFilterPP, "Run filter task", false);
  PROCESS_SWITCH(DQFilterPPTask, processFilterPPTiny, "Run filter task", false);
  PROCESS_SWITCH(DQFilterPPTask, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DQEventSelectionTask>(cfgc),
    adaptAnalysisTask<DQBarrelTrackSelection>(cfgc),
    adaptAnalysisTask<DQMuonsSelection>(cfgc),
    adaptAnalysisTask<DQFilterPPTask>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //
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

    if (classStr.Contains("Muon")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "muon");
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
