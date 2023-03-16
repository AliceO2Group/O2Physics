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

/// \file taskSingleMuon.cxx
/// \brief Task derived from the DQ framework and used to extract the observables on single muons needed for the HF-muon analysis.
/// \author Maolin Zhang <maolin.zhang@cern.ch>, CCNU

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;

namespace o2::aod
{
namespace dq_analysis_flags
{
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
}
DECLARE_SOA_TABLE(EventCuts, "AOD", "EVENTCUTS", dq_analysis_flags::IsEventSelected);
} // namespace o2::aod

using MyCollisions = o2::soa::Join<aod::Collisions, aod::EvSels>;
using MyMcCollisions = o2::soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
using MyEventsSelected = o2::soa::Join<aod::Collisions, aod::EvSels, o2::aod::EventCuts>;
using MyMcEventsSelected = o2::soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, o2::aod::EventCuts>;
using MyMuons = o2::soa::Join<aod::FwdTracks, aod::FwdTracksDCA>;
using MyMcMuons = o2::soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;

constexpr static uint32_t gMuonFillMap(VarManager::ObjTypes::Muon);
constexpr static uint32_t gEventFillMap(VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision);
constexpr static uint32_t gTrackMCFillMap(VarManager::ObjTypes::ParticleMC);

/// Defines the histograms for classes required in analysis.
void defineHistograms(HistogramManager*, TString);

struct HfTaskSingleMuonEventSelection {
  Produces<aod::EventCuts> eventSel;

  Configurable<bool> applySoftwareTrigger{"applySoftwareTrigger", false, "whether to apply the software trigger"};
  Configurable<int> softwareTrigger{"softwareTrigger", VarManager::kIsMuonSingleLowPt7, "software trigger flag"};
  Configurable<bool> applyCutZVtx{"applyCutZVtx", false, "whether to apply the VtxZ cut"};
  Configurable<float> zVtxMin{"zVtxMin", -10., "min. z of primary vertex [cm]"};
  Configurable<float> zVtxMax{"zVtxMax", 10., "max. z of primary vertex [cm]"};

  float* values;
  AnalysisCompositeCut* eventCut;
  HistogramManager* histMan;

  OutputObj<THashList> outputList{"output"};

  void init(o2::framework::InitContext&)
  {
    VarManager::SetDefaultVarNames();
    values = new float[VarManager::kNVars];
    eventCut = new AnalysisCompositeCut("event selection", "Event Selection", true);

    AnalysisCut cut;
    if (applyCutZVtx)
      cut.AddCut(VarManager::kVtxZ, zVtxMin, zVtxMax, false);
    if (applySoftwareTrigger)
      cut.AddCut(softwareTrigger, 0.5, 1.5, false);
    eventCut->AddCut(&cut);

    histMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    histMan->SetUseDefaultVariableNames(kTRUE);
    histMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    defineHistograms(histMan, "EventBeforeCuts;EventAfterCuts;");
    VarManager::SetUseVars(histMan->GetUsedVars());
    outputList.setObject(histMan->GetMainHistogramList());
  }

  template <uint32_t TEventFillMap, typename TEvent>
  void runEventSel(TEvent const& event, aod::BCs const& bcs)
  {
    // select events
    if (event.sel8()) {
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEvent<TEventFillMap>(event, values);
      histMan->FillHistClass("EventBeforeCuts", values);

      if (eventCut->IsSelected(values)) {
        histMan->FillHistClass("EventAfterCuts", values);
        eventSel(1);
      } else {
        eventSel(0);
      }
    } else {
      eventSel(0);
    }
  }

  void processEvent(MyCollisions::iterator const& event, aod::BCs const& bcs)
  {
    runEventSel<gEventFillMap>(event, bcs);
  }
  PROCESS_SWITCH(HfTaskSingleMuonEventSelection, processEvent, "run event selection with real data", true);

  void processEventMc(MyMcCollisions::iterator const& event, aod::BCs const& bcs)
  {
    runEventSel<gEventFillMap>(event, bcs);
  }
  PROCESS_SWITCH(HfTaskSingleMuonEventSelection, processEventMc, "run event selection with MC data", false);

  void processNoEventSelection(MyCollisions::iterator const& event, aod::BCs const& bcs)
  {
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::FillEvent<gEventFillMap>(event, values);
    eventSel(1);
  }
  PROCESS_SWITCH(HfTaskSingleMuonEventSelection, processNoEventSelection, "skip the event selection", false);
};

struct HfTaskSingleMuonSelection {
  Configurable<std::string> muonCuts{"muonCuts", "muonQualityCuts", "muon selection"};
  Configurable<bool> fillLiteHist{"fillLiteHist", true, "fill lite version muon histograms"};
  Configurable<bool> fillMcHist{"fillMcHist", false, "fill MC-related muon histograms"};

  float* values;
  AnalysisCompositeCut* trackCut;

  o2::framework::HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  void init(o2::framework::InitContext&)
  {
    AxisSpec axispT{200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{500, -5., 5., "#eta"};
    AxisSpec axisDCA{500, 0., 5., "DCA (cm)"};
    AxisSpec axisSign{5, -2.5, 2.5, "Charge"};
    AxisSpec axisP{500, 0., 500., "p (GeV/#it{c})"};
    AxisSpec axisVtxZ{400, -20., 20., "Vertex Z (cm)"};

    AxisSpec axisTrkType{5, -0.5, 4.5, "TrackType"};
    AxisSpec axisMCmask{200, -0.5, 199.5, "mcMask"};
    AxisSpec axisChi2MatchMCHMFT{100, 0., 100., "MCH-MFT matching #chi^{2}"};
    // kinematics for MC
    AxisSpec axispTGen{200, 0., 200., "#it{p}_{T} Truth (GeV/#it{c})"};
    AxisSpec axisEtaGen{500, -5., 5., "#eta_{Truth}"};
    AxisSpec axisPGen{500, 0., 500., "p_{Truth} (GeV/#it{c})"};
    AxisSpec axispTDif{200, -2., 2., "#it{p}_{T} diff (GeV/#it{c})"};
    AxisSpec axisEtaDif{200, -2., 2., "#eta diff"};
    AxisSpec axisPDif{200, -2., 2., "p diff (GeV/#it{c})"};

    HistogramConfigSpec hTHnMu{HistType::kTHnSparseD, {axispT, axisEta, axisDCA, axisSign, axisP, axisVtxZ, axisTrkType, axisChi2MatchMCHMFT, axisMCmask}, 9};
    HistogramConfigSpec hTHnPt{HistType::kTHnSparseD, {axispT, axispTGen, axispTDif, axisTrkType}, 4};
    HistogramConfigSpec hTHnEta{HistType::kTHnSparseD, {axisEta, axisEtaGen, axisEtaDif, axisTrkType}, 4};
    HistogramConfigSpec hTHnP{HistType::kTHnSparseD, {axisP, axisPGen, axisPDif, axisTrkType}, 4};
    HistogramConfigSpec hTHnMuLite{HistType::kTHnSparseD, {axispT, axisEta, axisDCA, axisTrkType, axisChi2MatchMCHMFT, axisMCmask}, 6};

    if (fillLiteHist) {
      registry.add("hMuLiteBeforeCuts", "", hTHnMuLite);
      registry.add("hMuLiteAfterCuts", "", hTHnMuLite);
    } else {
      registry.add("hMuBeforeCuts", "", hTHnMu);
      registry.add("hMuAfterCuts", "", hTHnMu);
    }
    if (fillMcHist) {
      registry.add("hPtBeforeCuts", "", hTHnPt);
      registry.add("hPtAfterCuts", "", hTHnPt);
      registry.add("hEtaBeforeCuts", "", hTHnEta);
      registry.add("hEtaAfterCuts", "", hTHnEta);
      registry.add("hPBeforeCuts", "", hTHnP);
      registry.add("hPAfterCuts", "", hTHnP);
    }

    VarManager::SetDefaultVarNames();
    values = new float[VarManager::kNVars];

    trackCut = new AnalysisCompositeCut(true);
    TString selectStr = muonCuts.value;
    trackCut->AddCut(dqcuts::GetAnalysisCut(selectStr.Data()));
    VarManager::SetUseVars(AnalysisCut::fgUsedVars);
  }

  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvent, typename TMuons>
  void runMuonSel(TEvent const& event, aod::BCs const& bcs, TMuons const& tracks)
  {
    // select muons in data
    if (event.isEventSelected() == 0) {
      return;
    }

    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables, values);
    VarManager::FillEvent<gEventFillMap>(event, values);

    // loop over muon tracks
    for (auto const& track : tracks) {
      VarManager::FillTrack<TMuonFillMap>(track, values);

      // compute DCAXY
      // Same for MC
      const auto dcaXY(std::sqrt(values[VarManager::kMuonDCAx] * values[VarManager::kMuonDCAx] + values[VarManager::kMuonDCAy] * values[VarManager::kMuonDCAy]));

      if (fillLiteHist) {
        // Before Muon Cuts
        registry.fill(HIST("hMuLiteBeforeCuts"),
                      values[VarManager::kPt],
                      values[VarManager::kEta], dcaXY,
                      values[VarManager::kMuonTrackType], values[VarManager::kMuonChi2MatchMCHMFT], 0);
        // After Muon Cuts
        if (trackCut->IsSelected(values)) {
          registry.fill(HIST("hMuLiteAfterCuts"),
                        values[VarManager::kPt],
                        values[VarManager::kEta], dcaXY,
                        values[VarManager::kMuonTrackType], values[VarManager::kMuonChi2MatchMCHMFT], 0);
        }
      } else {
        // Before Muon Cuts
        registry.fill(HIST("hMuBeforeCuts"),
                      values[VarManager::kPt],
                      values[VarManager::kEta], dcaXY,
                      values[VarManager::kCharge], track.p(),
                      values[VarManager::kVtxZ],
                      values[VarManager::kMuonTrackType], values[VarManager::kMuonChi2MatchMCHMFT], 0);
        // After Muon Cuts
        if (trackCut->IsSelected(values)) {
          registry.fill(HIST("hMuAfterCuts"),
                        values[VarManager::kPt],
                        values[VarManager::kEta], dcaXY,
                        values[VarManager::kCharge], track.p(),
                        values[VarManager::kVtxZ],
                        values[VarManager::kMuonTrackType], values[VarManager::kMuonChi2MatchMCHMFT], 0);
        }
      }
    } // end loop over muon tracks
  }

  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, uint32_t TTrackMCFillMap, typename TEvent, typename TMuons, typename TMC>
  void runMuonSelMc(TEvent const& event, aod::BCs const& bcs, TMuons const& tracks, TMC const& mc)
  {
    // select muons in MC
    if (event.isEventSelected() == 0) {
      return;
    }

    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables, values);
    VarManager::FillEvent<gEventFillMap>(event, values);

    // loop over muon tracks
    for (auto const& track : tracks) {
      VarManager::FillTrack<TMuonFillMap>(track, values);

      if (!track.has_mcParticle()) {
        continue;
      }
      auto mcParticle = track.mcParticle();

      VarManager::FillTrack<TTrackMCFillMap>(mcParticle, values);

      // compute DCAXY
      const auto dcaXY(std::sqrt(values[VarManager::kMuonDCAx] * values[VarManager::kMuonDCAx] + values[VarManager::kMuonDCAy] * values[VarManager::kMuonDCAy]));

      if (fillLiteHist) {
        // Before Muon Cuts
        registry.fill(HIST("hMuLiteBeforeCuts"),
                      values[VarManager::kPt],
                      values[VarManager::kEta], dcaXY,
                      values[VarManager::kMuonTrackType], values[VarManager::kMuonChi2MatchMCHMFT],
                      track.mcMask());
        if (fillMcHist) {
          registry.fill(HIST("hPtBeforeCuts"),
                        values[VarManager::kPt], values[VarManager::kMCPt], values[VarManager::kMCPt] - values[VarManager::kPt], values[VarManager::kMuonTrackType]);
          registry.fill(HIST("hEtaBeforeCuts"),
                        values[VarManager::kEta], values[VarManager::kMCEta], values[VarManager::kMCEta] - values[VarManager::kEta], values[VarManager::kMuonTrackType]);
          registry.fill(HIST("hPBeforeCuts"),
                        values[VarManager::kP], values[VarManager::kMCPt] * std::cosh(values[VarManager::kMCEta]), values[VarManager::kMCPt] * std::cosh(values[VarManager::kMCEta]) - values[VarManager::kP], values[VarManager::kMuonTrackType]);
        }
        // After Muon Cuts
        if (trackCut->IsSelected(values)) {
          registry.fill(HIST("hMuLiteAfterCuts"),
                        values[VarManager::kPt],
                        values[VarManager::kEta], dcaXY,
                        values[VarManager::kMuonTrackType], values[VarManager::kMuonChi2MatchMCHMFT],
                        track.mcMask());
          if (fillMcHist) {
            registry.fill(HIST("hPtAfterCuts"),
                          values[VarManager::kPt], values[VarManager::kMCPt], values[VarManager::kMCPt] - values[VarManager::kPt], values[VarManager::kMuonTrackType]);
            registry.fill(HIST("hEtaAfterCuts"),
                          values[VarManager::kEta], values[VarManager::kMCEta], values[VarManager::kMCEta] - values[VarManager::kEta], values[VarManager::kMuonTrackType]);
            registry.fill(HIST("hPAfterCuts"),
                          values[VarManager::kP], values[VarManager::kMCPt] * std::cosh(values[VarManager::kMCEta]), values[VarManager::kMCPt] * std::cosh(values[VarManager::kMCEta]) - values[VarManager::kP], values[VarManager::kMuonTrackType]);
          }
        }
      } else {
        // Before Muon Cuts
        registry.fill(HIST("hMuBeforeCuts"),
                      values[VarManager::kPt],
                      values[VarManager::kEta], dcaXY,
                      values[VarManager::kCharge], track.p(),
                      values[VarManager::kVtxZ],
                      values[VarManager::kMuonTrackType], values[VarManager::kMuonChi2MatchMCHMFT],
                      track.mcMask());
        if (fillMcHist) {
          registry.fill(HIST("hPtBeforeCuts"),
                        values[VarManager::kPt], values[VarManager::kMCPt], values[VarManager::kMCPt] - values[VarManager::kPt], values[VarManager::kMuonTrackType]);
          registry.fill(HIST("hEtaBeforeCuts"),
                        values[VarManager::kEta], values[VarManager::kMCEta], values[VarManager::kMCEta] - values[VarManager::kEta], values[VarManager::kMuonTrackType]);
          registry.fill(HIST("hPBeforeCuts"),
                        values[VarManager::kP], values[VarManager::kMCPt] * std::cosh(values[VarManager::kMCEta]), values[VarManager::kMCPt] * std::cosh(values[VarManager::kMCEta]) - values[VarManager::kP], values[VarManager::kMuonTrackType]);
        }
        // After Muon Cuts
        if (trackCut->IsSelected(values)) {
          registry.fill(HIST("hMuAfterCuts"),
                        values[VarManager::kPt],
                        values[VarManager::kEta], dcaXY,
                        values[VarManager::kCharge], track.p(),
                        values[VarManager::kVtxZ],
                        values[VarManager::kMuonTrackType], values[VarManager::kMuonChi2MatchMCHMFT],
                        track.mcMask());
          if (fillMcHist) {
            registry.fill(HIST("hPtAfterCuts"),
                          values[VarManager::kPt], values[VarManager::kMCPt], values[VarManager::kMCPt] - values[VarManager::kPt], values[VarManager::kMuonTrackType]);
            registry.fill(HIST("hEtaAfterCuts"),
                          values[VarManager::kEta], values[VarManager::kMCEta], values[VarManager::kMCEta] - values[VarManager::kEta], values[VarManager::kMuonTrackType]);
            registry.fill(HIST("hPAfterCuts"),
                          values[VarManager::kP], values[VarManager::kMCPt] * std::cosh(values[VarManager::kMCEta]), values[VarManager::kMCPt] * std::cosh(values[VarManager::kMCEta]) - values[VarManager::kP], values[VarManager::kMuonTrackType]);
          }
        }
      }
    } // end loop over muon tracks
  }

  void processDummy(MyEventsSelected&)
  {
    // do nothing
  }

  void processMuon(MyEventsSelected::iterator const& event, aod::BCs const& bcs,
                   MyMuons const& tracks)
  {
    runMuonSel<gEventFillMap, gMuonFillMap>(event, bcs, tracks);
  }

  void processMuonMc(MyMcEventsSelected::iterator const& event, aod::BCs const& bcs,
                     MyMcMuons const& tracks, aod::McParticles const& mc)
  {
    runMuonSelMc<gEventFillMap, gMuonFillMap, gTrackMCFillMap>(event, bcs, tracks, mc);
  }

  PROCESS_SWITCH(HfTaskSingleMuonSelection, processDummy, "do nothing", false);
  PROCESS_SWITCH(HfTaskSingleMuonSelection, processMuon, "run muon selection with real data", false);
  PROCESS_SWITCH(HfTaskSingleMuonSelection, processMuonMc, "run muon selection with MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskSingleMuonEventSelection>(cfgc),
    adaptAnalysisTask<HfTaskSingleMuonSelection>(cfgc),
  };
}

void defineHistograms(HistogramManager* histMan, TString histClasses)
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

    if (classStr.Contains("Event"))
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "trigger,all");
  }
}
