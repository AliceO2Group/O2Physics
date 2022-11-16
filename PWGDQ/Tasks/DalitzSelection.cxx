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
//                Task to select electrons from dalitz decay
//        It can produce track and pair histograms for selected tracks
//      It creates a bitmap with selections to be stored during skimming
//
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/DataTypes.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;
using std::array;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyReducedEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;

using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA,
                                 aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullMu,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullMu,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using BarrelReducedTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::Collision;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkReducedEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkReducedTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;

struct dalitzPairing {
  Produces<o2::aod::DalitzBits> dalitzbits;
  Preslice<MyBarrelTracks> perCollision = aod::track::collisionId;

  // Configurables
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandardNoINT7", "Event selection"};
  Configurable<std::string> fConfigDalitzTrackCuts{"cfgDalitzTrackCuts", "", "Dalitz track selection cuts, separated by a comma"};
  Configurable<std::string> fConfigDalitzPairCuts{"cfgDalitzPairCuts", "", "Dalitz pair selection cuts"};
  Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
  Configurable<bool> fQA{"cfgQA", true, "QA histograms"};
  Configurable<float> fConfigBarrelTrackPINLow{"cfgBarrelLowPIN", 0.1f, "Low pt cut for Dalitz tracks in the barrel"};
  Configurable<float> fConfigEtaCut{"cfgEtaCut", 0.9f, "Eta cut for Dalitz tracks in the barrel"};
  Configurable<float> fConfigTPCNSigLow{"cfgTPCNSigElLow", -3.f, "Low TPCNSigEl cut for Dalitz tracks in the barrel"};
  Configurable<float> fConfigTPCNSigHigh{"cfgTPCNSigElHigh", 3.f, "High TPCNsigEl cut for Dalitz tracks in the barrel"};

  Filter filterBarrelTrack = o2::aod::track::tpcInnerParam >= fConfigBarrelTrackPINLow && nabs(o2::aod::track::eta) <= fConfigEtaCut && o2::aod::pidtpc::tpcNSigmaEl <= fConfigTPCNSigHigh && o2::aod::pidtpc::tpcNSigmaEl >= fConfigTPCNSigLow;

  OutputObj<THashList> fOutputList{"output"}; //! the histogram manager output list
  OutputObj<TList> fStatsList{"Statistics"};  //! skimming statistics

  std::map<int, uint8_t> trackmap;
  std::map<int, uint8_t> dalitzmap;

  AnalysisCompositeCut* fEventCut;
  std::vector<AnalysisCompositeCut> fTrackCuts;
  std::vector<AnalysisCompositeCut> fPairCuts;
  int nCuts = 0;

  HistogramManager* fHistMan;

  void init(o2::framework::InitContext&)
  {
    // Event cuts
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));

    // Barrel track cuts
    TString cutNamesStr = fConfigDalitzTrackCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // Pair cuts
    TString cutNamesPairStr = fConfigDalitzPairCuts.value;
    if (!cutNamesPairStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesPairStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    if (fTrackCuts.size() != fPairCuts.size()) {
      std::cout << "WARNING: YOU SHOULD PROVIDE THE SAME NUMBER OF TRACK AND PAIR CUTS" << std::endl;
    }
    nCuts = std::min(fTrackCuts.size(), fPairCuts.size());

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Create the histogram class names to be added to the histogram manager
    TString histClasses = "";

    if (fQA) {
      for (int icut = 0; icut < nCuts; icut++) {
        AnalysisCompositeCut trackCut = fTrackCuts.at(icut);
        AnalysisCompositeCut pairCut = fPairCuts.at(icut);
        histClasses += Form("TrackBarrel_%s_%s;", trackCut.GetName(), pairCut.GetName());
        histClasses += Form("Pair_%s_%s;", trackCut.GetName(), pairCut.GetName());
      }
    }

    // Define histograms

    std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
    for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
      TString classStr = objArray->At(iclass)->GetName();
      fHistMan->AddHistClass(classStr.Data());

      if (classStr.Contains("Event")) {
        dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", "");
      }
      TString histTrackName = fConfigAddTrackHistogram.value;
      if (classStr.Contains("Track")) {
        dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histTrackName);
      }

      if (classStr.Contains("Pair")) {
        dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "pair_barrel", "dalitz");
      }
    }

    fStatsList.setObject(new TList());
    fStatsList->SetOwner(kTRUE);

    // Dalitz selection statistics: one bin for each (track,pair) selection
    TH1I* histTracks = new TH1I("TrackStats", "Dalitz selection statistics", nCuts, -0.5, nCuts - 0.5);
    for (int icut = 0; icut < nCuts; icut++) {
      AnalysisCompositeCut trackCut = fTrackCuts.at(icut);
      AnalysisCompositeCut pairCut = fPairCuts.at(icut);
      histTracks->GetXaxis()->SetBinLabel(icut + 1, Form("%s_%s", trackCut.GetName(), pairCut.GetName()));
    }
    if (fQA) {
      fStatsList->Add(histTracks);
    }

    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  template <uint32_t TTrackFillMap, typename TTracks>
  void runTrackSelection(TTracks const& tracksBarrel)
  {
    for (auto& track : tracksBarrel) {
      uint8_t filterMap = uint8_t(0);
      VarManager::FillTrack<TTrackFillMap>(track);
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); ++cut, ++i) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint8_t(1) << i);
        }
      }
      if (filterMap) {
        trackmap[track.globalIndex()] = filterMap;
      }
    } // end loop over tracks
  }

  template <uint32_t TTrackFillMap, typename TTracks>
  void runDalitzPairing(TTracks const& tracks1, TTracks const& tracks2)
  {
    const int TPairType = VarManager::kJpsiToEE; // For dielectron
    for (auto& [track1, track2] : o2::soa::combinations(CombinationsStrictlyUpperIndexPolicy(tracks1, tracks2))) {
      if (track1.sign() * track2.sign() > 0) {
        continue;
      }

      uint8_t twoTracksFilterMap = trackmap[track1.globalIndex()] & trackmap[track2.globalIndex()];
      if (!twoTracksFilterMap)
        continue;

      // pairing
      VarManager::FillPair<TPairType, TTrackFillMap>(track1, track2);
      uint8_t track1Untagged = uint8_t(0);
      uint8_t track2Untagged = uint8_t(0);

      // Fill pair selection map and fill pair histogram
      for (int icut = 0; icut < nCuts; icut++) {
        if (!(twoTracksFilterMap & (uint8_t(1) << icut)))
          continue;
        AnalysisCompositeCut pairCut = fPairCuts.at(icut);
        if (pairCut.IsSelected(VarManager::fgValues)) {
          if (fQA) {
            AnalysisCompositeCut trackCut = fTrackCuts.at(icut);
            fHistMan->FillHistClass(Form("Pair_%s_%s", trackCut.GetName(), pairCut.GetName()), VarManager::fgValues);
          }

          // Check if tracks were already tagged
          bool b1 = dalitzmap[track1.globalIndex()] & (uint8_t(1) << icut);
          if (!b1) {
            track1Untagged |= (uint8_t(1) << icut);
            if (fQA) {
              ((TH1I*)fStatsList->At(0))->Fill(icut);
            }
          }
          bool b2 = dalitzmap[track2.globalIndex()] & (uint8_t(1) << icut);
          if (!b2) {
            track2Untagged |= (uint8_t(1) << icut);
            if (fQA) {
              ((TH1I*)fStatsList->At(0))->Fill(icut);
            }
          }
        }
      }

      // Tag tracks which are not already tagged
      dalitzmap[track1.globalIndex()] |= track1Untagged;
      dalitzmap[track2.globalIndex()] |= track2Untagged;

      // Fill track histograms if not already tagged
      if (fQA) {
        VarManager::FillTrack<TTrackFillMap>(track1);
        for (int icut = 0; icut < nCuts; icut++) {
          if (track1Untagged & (uint8_t(1) << icut)) {
            AnalysisCompositeCut trackCut = fTrackCuts.at(icut);
            AnalysisCompositeCut pairCut = fPairCuts.at(icut);
            fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", trackCut.GetName(), pairCut.GetName()), VarManager::fgValues);
          }
        }
        VarManager::FillTrack<TTrackFillMap>(track2);
        for (int icut = 0; icut < nCuts; icut++) {
          if (track2Untagged & (uint8_t(1) << icut)) {
            AnalysisCompositeCut trackCut = fTrackCuts.at(icut);
            AnalysisCompositeCut pairCut = fPairCuts.at(icut);
            fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", trackCut.GetName(), pairCut.GetName()), VarManager::fgValues);
          }
        }
      }
    } // end of tracksP,N loop
  }

  void processFullTracks(MyEvents const& collisions, soa::Filtered<MyBarrelTracks> const& filteredTracks, MyBarrelTracks const& tracks)
  {
    dalitzmap.clear();

    for (auto& collision : collisions) {
      trackmap.clear();
      VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
      VarManager::FillEvent<gkEventFillMap>(collision);
      bool isEventSelected = fEventCut->IsSelected(VarManager::fgValues);

      if (isEventSelected) {
        auto groupedFilteredTracks = filteredTracks.sliceBy(perCollision, collision.globalIndex());
        runTrackSelection<gkTrackFillMap>(groupedFilteredTracks);
        runDalitzPairing<gkTrackFillMap>(groupedFilteredTracks, groupedFilteredTracks);
      }
    }

    for (auto& track : tracks) {// Fill dalitz bits
      dalitzbits(static_cast<int>(dalitzmap[track.globalIndex()]));
    }
  }

  void processDummy(MyEvents&)
  {
  }

  PROCESS_SWITCH(dalitzPairing, processFullTracks, "Run Dalitz selection on AO2D tables", false);
  PROCESS_SWITCH(dalitzPairing, processDummy, "Do nothing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dalitzPairing>(cfgc)};
}
