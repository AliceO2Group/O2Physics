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
DECLARE_SOA_COLUMN(IsBarrelSelected, isBarrelSelected, uint8_t);
DECLARE_SOA_COLUMN(IsMuonSelected, isMuonSelected, uint8_t);
} // namespace dqppfilter

DECLARE_SOA_TABLE(DQEventCuts, "AOD", "DQEVENTCUTS", dqppfilter::IsDQEventSelected);
DECLARE_SOA_TABLE(DQBarrelTrackCuts, "AOD", "DQBARRELCUTS", dqppfilter::IsBarrelSelected);
DECLARE_SOA_TABLE(DQMuonsCuts, "AOD", "DQMUONCUTS", dqppfilter::IsMuonSelected);
} // namespace o2::aod

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsSelected = soa::Join<aod::Collisions, aod::EvSels, aod::DQEventCuts>;
using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksExtended, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyBarrelTracksSelected = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksExtended, aod::TrackSelection,
                                         aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                         aod::pidTPCFullKa, aod::pidTPCFullPr,
                                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                         aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                         aod::DQBarrelTrackCuts>;
using MyMuons = aod::FwdTracks;
using MyMuonsSelected = soa::Join<aod::FwdTracks, aod::DQMuonsCuts>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::Muon;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

struct DQEventSelectionTask {
  Produces<aod::DQEventCuts> eventSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;
  AnalysisCompositeCut fEventCut{true};

  float* fValues;

  Configurable<std::string> fConfigCuts{"cfgEventCuts", "eventStandard", "Comma separated list of event cuts; multiple cuts are applied with a logical AND"};

  void init(o2::framework::InitContext&)
  {
    fValues = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;"); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());                 // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    DefineCuts();
  }

  void DefineCuts()
  {
    // default cut is "eventStandard" (kINT7 and vtxZ selection)
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fEventCut.AddCut(dqcuts::GetAnalysisCut(objArray->At(icut)->GetName()));
      }
    }

    // NOTE: Additional cuts to those specified via the Configurable may still be added

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  void process(MyEvents::iterator const& event, aod::BCs const& bcs)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables, fValues);

    VarManager::FillEvent<gkEventFillMap>(event, fValues);
    fHistMan->FillHistClass("Event_BeforeCuts", fValues); // automatically fill all the histograms in the class Event
    if (fEventCut.IsSelected(fValues)) {
      fHistMan->FillHistClass("Event_AfterCuts", fValues);
      eventSel(1);
    } else {
      eventSel(0);
    }
  }
};

struct DQBarrelTrackSelectionTask {
  Produces<aod::DQBarrelTrackCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fTrackCuts;

  float* fValues; // array to be used by the VarManager

  Configurable<std::string> fConfigCuts{"cfgBarrelTrackCuts", "jpsiPID2", "Comma separated list of ADDITIONAL barrel track cuts"};

  void init(o2::framework::InitContext&)
  {
    DefineCuts();

    fValues = new float[VarManager::kNVars];
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

  void process(MyEvents::iterator const& event, MyBarrelTracks const& tracks, aod::BCs const& bcs)
  {
    uint8_t filterMap = uint8_t(0);
    trackSel.reserve(tracks.size());

    VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables, fValues);
    // fill event information which might be needed in histograms that combine track and event properties
    VarManager::FillEvent<gkEventFillMap>(event, fValues);

    for (auto& track : tracks) {
      filterMap = uint8_t(0);
      // TODO: if a condition on the event selection is applied, the histogram output is missing
      VarManager::FillTrack<gkTrackFillMap>(track, fValues);
      fHistMan->FillHistClass("TrackBarrel_BeforeCuts", fValues);
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); ++cut, ++i) {
        if ((*cut).IsSelected(fValues)) {
          filterMap |= (uint8_t(1) << i);
          fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), fValues);
        }
      }
      trackSel(filterMap);
    }
  }
};

struct DQMuonsSelection {
  Produces<aod::DQMuonsCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fTrackCuts;

  float* fValues;

  Configurable<std::string> fConfigCuts{"cfgMuonsCuts", "muonQualityCuts", "Comma separated list of ADDITIONAL muon track cuts"};

  void init(o2::framework::InitContext&)
  {
    DefineCuts();

    fValues = new float[VarManager::kNVars];
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

  void process(MyEvents::iterator const& event, MyMuons const& muons, aod::BCs const& bcs)
  {
    uint8_t filterMap = uint8_t(0);
    trackSel.reserve(muons.size());

    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables, fValues);
    // fill event information which might be needed in histograms that combine track and event properties
    VarManager::FillEvent<gkEventFillMap>(event, fValues);

    for (auto& muon : muons) {
      filterMap = uint8_t(0);
      // TODO: if a condition on the event selection is applied, the histogram output is missing
      VarManager::FillTrack<gkMuonFillMap>(muon, fValues);
      fHistMan->FillHistClass("TrackMuons_BeforeCuts", fValues);
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); ++cut, ++i) {
        if ((*cut).IsSelected(fValues)) {
          filterMap |= (uint8_t(1) << i);
          fHistMan->FillHistClass(Form("TrackMuons_%s", (*cut).GetName()), fValues);
        }
      }
      trackSel(filterMap);
    }
  }
};

struct DQFilterPPTask {
  Produces<aod::DQEventFilter> eventFilter;
  Produces<aod::DqFilters> dqtable;
  OutputObj<THashList> fOutputList{"output"};
  OutputObj<TH1I> fStatsSingleBarrel{"statsSingleBarrel"};
  OutputObj<TH1I> fStatsSingleMuon{"statsSingleMuon"};
  OutputObj<TH2I> fStatsBarrel{"statsBarrel"};
  OutputObj<TH2I> fStatsMuon{"statsMuon"};
  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fPairCuts;

  float* fValues;

  Partition<MyBarrelTracksSelected> posTracks = aod::track::signed1Pt > 0.0f && aod::dqppfilter::isBarrelSelected > uint8_t(0);
  Partition<MyBarrelTracksSelected> negTracks = aod::track::signed1Pt < 0.0f && aod::dqppfilter::isBarrelSelected > uint8_t(0);

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

    fValues = new float[VarManager::kNVars];

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

  void
    process(MyEventsSelected::iterator const& event, MyBarrelTracksSelected const& tracks, MyMuonsSelected const& muons, aod::BCs const& bcs)
  {
    uint64_t filter = 0; // first 32 bits for dielectron, last 32 for dimuons

    //[0] = SingleE, [1] = SingleMulow, [2] = SingleMuHigh, [3] = DiElectron, [4] = DiMuon
    bool keepEvent[kNTriggersDQ]{false};

    constexpr int pairTypeEE = VarManager::kJpsiToEE;
    constexpr int pairTypeMuMu = VarManager::kJpsiToMuMu;
    if (event.isDQEventSelected() == 1) {
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars, fValues);
      VarManager::FillEvent<gkEventFillMap>(event, fValues);
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
      for (auto track : tracks) {
        cutSingleFilter = track.isBarrelSelected();
        if (!cutSingleFilter) {
          continue;
        }
        for (int i = 0; i < fNTrackCuts; ++i) {
          if (!(cutSingleFilter & (uint8_t(1) << i))) {
            continue;
          }
          VarManager::FillTrack<gkTrackFillMap>(track, fValues);
          fHistMan->FillHistClass(Form("TrackBarrel_%s", fTrkCutsNameArray->At(i)->GetName()), fValues);
          singleBarrelCount[i] += 1;
        }
      }

      uint8_t cutFilter = 0;
      for (auto tpos : posTracks) {
        for (auto tneg : negTracks) { // +- pairs
          cutFilter = tpos.isBarrelSelected() & tneg.isBarrelSelected();
          if (!cutFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }
          VarManager::FillPair<pairTypeEE, gkTrackFillMap>(tpos, tneg, fValues); // compute pair quantities
          for (int i = 0; i < fNTrackCuts; ++i) {
            if (!(cutFilter & (uint8_t(1) << i))) {
              continue;
            }
            fHistMan->FillHistClass(Form("PairsBarrelPM_%s", fTrkCutsNameArray->At(i)->GetName()), fValues);
            int j = 0;
            for (auto cut = fPairCuts.begin(); cut != fPairCuts.end(); cut++, j++) {
              if ((*cut).IsSelected(fValues)) {
                pairsBarrelCount[j + i * fNPairCuts] += 1;
              }
            }
          }
        }
      }

      uint8_t cutSingleMuonFilter = 0;
      for (auto muon1 : muons) {
        cutSingleMuonFilter = muon1.isMuonSelected();
        if (!cutSingleMuonFilter) {
          continue;
        }
        for (int i = 0; i < fNMuonCuts; ++i) {
          if (!(cutSingleMuonFilter & (uint8_t(1) << i))) {
            continue;
          }
          VarManager::FillTrack<gkMuonFillMap>(muon1, fValues);
          fHistMan->FillHistClass(Form("TrackMuons_%s", fMuonCutsNameArray->At(i)->GetName()), fValues);
          singleMuonsCount[i] += 1;
        }
      }

      uint8_t cutMuonFilter = 0;
      for (auto& [muon1, muon2] : combinations(muons, muons)) {
        cutMuonFilter = muon1.isMuonSelected() & muon2.isMuonSelected();
        if (!cutMuonFilter) { // the tracks must have at least one filter bit in common to continue
          continue;
        }
        VarManager::FillPair<pairTypeMuMu, gkMuonFillMap>(muon1, muon2, fValues); // compute pair quantities
        for (int i = 0; i < fNMuonCuts; ++i) {
          if (!(cutMuonFilter & (uint8_t(1) << i))) {
            continue;
          }
          if (muon1.sign() * muon2.sign() < 0) {
            fHistMan->FillHistClass(Form("PairsMuonsPM_%s", fMuonCutsNameArray->At(i)->GetName()), fValues);
          } else {
            if (muon1.sign() > 0) {
              fHistMan->FillHistClass(Form("PairsMuonsPP_%s", fMuonCutsNameArray->At(i)->GetName()), fValues);
            } else {
              fHistMan->FillHistClass(Form("PairsMuonsMM_%s", fMuonCutsNameArray->At(i)->GetName()), fValues);
            }
          }

          int j = 0;
          for (auto cut = fPairCuts.begin(); cut != fPairCuts.end(); cut++, j++) {
            if ((*cut).IsSelected(fValues)) {
              if (muon1.sign() * muon2.sign() < 0) {
                pairsUnlikeMuonsCount[j + i * fNPairCuts] += 1;
              } else {
                pairsLikeMuonsCount[j + i * fNPairCuts] += 1;
              }
            }
          }
        }
      }

      // Fill event selections tag for central event filtering
      if (singleBarrelCount[0] > 0) { // Single Barrel count : bit 0 corresponds to "jpsiPID1" cut
        keepEvent[0] = true;
      }
      if (singleMuonsCount[0] > 0) { // Single Muon count : bit 0 corresponds to "muonLowPt" cut
        keepEvent[1] = true;
      }
      if (singleMuonsCount[1] > 0) { // Single Muon count : bit 1 corresponds to "muonHighPt" cut
        keepEvent[2] = true;
      }
      if (pairsBarrelCount[0] > 0) { // Pair Barrel count : bit 0 corresponds "jpsiPID1" & "pairNoCut"  [jPaircut + iTrackCut * fNPairCuts]
        keepEvent[3] = true;
      }
      if (pairsUnlikeMuonsCount[0] > 0) { // Pair Muon count : bit 0 corresponds "muonLowPt" & "pairNoCut"  [jPaircut + iMuonCut * fNPairCuts]
        keepEvent[4] = true;
      }
      // Filling the table
      dqtable(keepEvent[kSingleE], keepEvent[kSingleMuLow], keepEvent[kSingleMuHigh], keepEvent[kDiElectron], keepEvent[kDiMuon]);

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
    }
    eventFilter(filter);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DQEventSelectionTask>(cfgc),
    adaptAnalysisTask<DQBarrelTrackSelectionTask>(cfgc),
    adaptAnalysisTask<DQMuonsSelection>(cfgc),
    adaptAnalysisTask<DQFilterPPTask>(cfgc)};
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
