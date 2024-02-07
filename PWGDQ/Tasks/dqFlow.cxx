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
///       Task to compute Q vectors and other quantities related from the generic framework
///       Generic framework O2 version is a port of the AliPhysics version
///       To be used in the DQ analyses aiming for flow measurements
///       Run the standalone task with:
///       o2-analysis-timestamp --aod-file AO2D.root -b | o2-analysis-event-selection -b | o2-analysis-multiplicity-table -b | o2-analysis-centrality-table -b | o2-analysis-fdd-converter -b | o2-analysis-trackselection -b | o2-analysis-trackextension -b | o2-analysis-pid-tpc-full -b | o2-analysis-pid-tof-full -b | o2-analysis-pid-tof-base -b | o2-analysis-pid-tof-beta -b | o2-analysis-dq-flow -b
///       tested (June 2, 2022) on AO2D.root files from train production 242

#include <TH1F.h>
#include <THashList.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <TRandom3.h>
#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWCumulant.h"
#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"

using std::complex;
using std::cout;
using std::endl;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Declarations of various short names
using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;
using MyEventsWithCentRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;

using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyBarrelTracksWithCov = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection,
                                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                        aod::pidTPCFullKa, aod::pidTPCFullPr,
                                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                        aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;
using MyMuons = aod::FwdTracks;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCentRun2;
constexpr static uint32_t gkEventFillMapRun3 = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

struct DQEventQvector {
  Produces<ReducedEventsQvector> eventQvector;

  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<bool> fConfigQA{"cfgQA", true, "If true, fill QA histograms"};
  Configurable<bool> fConfigFlow{"cfgFlow", true, "If true, fill Flow histograms"};
  Configurable<float> fConfigCutPtMin{"cfgCutPtMin", 0.2f, "Minimal pT for tracks"};
  Configurable<float> fConfigCutPtMax{"cfgCutPtMax", 12.0f, "Maximal pT for tracks"};
  Configurable<float> fConfigCutEtaMin{"cfgCutEtaMin", -0.8f, "Eta min range for tracks"};
  Configurable<float> fConfigCutEtaMax{"cfgCutEtaMax", 0.8f, "Eta max range for tracks"};
  Configurable<float> fConfigEtaLimitMin{"cfgEtaLimitMin", -0.4f, "Eta gap min separation, only if using subEvents"};
  Configurable<float> fConfigEtaLimitMax{"cfgEtaLimitMax", 0.4f, "Eta gap max separation, only if using subEvents"};
  Configurable<uint8_t> fConfigNPow{"cfgNPow", 0, "Power of weights for Q vector"};

  // Access to the efficiencies and acceptances from CCDB
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> fConfigEfficiency{"ccdb-path-efficiency", "Users/r/rcaron/efficiency", "CCDB path to efficiency object"};
  Configurable<std::string> fConfigAcceptance{"ccdb-path-acceptance", "", "CCDB path to acceptance or GFWWeights object"};
  Configurable<std::string> fConfigURL{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  // Configurables for FlowContainer (e.g charged particles pt-differential v2{2}, v2{3}, ...)
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1}, "multiplicity / centrality axis for histograms"};

  // Define the filter for barrel tracks and forward tracks
  Filter trackFilter = (nabs(aod::track::eta) <= fConfigCutEtaMax) && (aod::track::pt > fConfigCutPtMin) && (aod::track::pt < fConfigCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));
  Filter fwdFilter = (aod::fwdtrack::eta < -2.45f) && (aod::fwdtrack::eta > -3.6f);

  // Histograms used for optionnal efficiency and non-uniform acceptance corrections
  struct Config {
    TH1D* mEfficiency = nullptr;
    GFWWeights* mAcceptance = nullptr;
  } cfg;

  // Define output
  HistogramManager* fHistMan = nullptr;
  AnalysisCompositeCut* fEventCut;
  OutputObj<THashList> fOutputList{"outputQA"};
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};

  // Define global variables for generic framework
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TRandom3* fRndm = new TRandom3(0);

  // Initialize CCDB, efficiencies and acceptances from CCDB, histograms, GFW, FlowContainer
  void init(o2::framework::InitContext&)
  {
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut = new AnalysisCompositeCut(true);
    if (!eventCutStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(eventCutStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fEventCut->AddCut(dqcuts::GetAnalysisCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    ccdb->setURL(fConfigURL.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(fConfigNoLaterThan.value);

    VarManager::SetDefaultVarNames();

    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts"); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    // Global effiencies
    if (fConfigEfficiency.value.empty() == false) {
      cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(fConfigEfficiency.value, fConfigNoLaterThan.value);
    }

    // Global acceptance or GFWeights to correct for NUA in the track loop
    if (fConfigAcceptance.value.empty() == false) {
      cfg.mAcceptance = ccdb->getForTimeStamp<GFWWeights>(fConfigAcceptance.value, fConfigNoLaterThan.value);
    }

    // Reference flow
    TObjArray* oba = new TObjArray();
    oba->Add(new TNamed("ChGap22", "ChGap22"));   // for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap24", "ChGap24"));   // for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChFull22", "ChFull22")); // no-gap case
    oba->Add(new TNamed("ChFull24", "ChFull24")); // no-gap case
    oba->Add(new TNamed("ChGap32", "ChGap32"));   // gap-case
    fFC->SetName("FlowContainer");
    fFC->Initialize(oba, axisMultiplicity, 10);
    delete oba;

    // Define regions of positive and negative eta in order to create gaps
    fGFW->AddRegion("refN", fConfigCutEtaMin, fConfigEtaLimitMin, 1, 1);
    fGFW->AddRegion("refP", fConfigEtaLimitMax, fConfigCutEtaMax, 1, 1);
    fGFW->AddRegion("full", fConfigCutEtaMin, fConfigCutEtaMax, 1, 2);
    // Defined the different charged particle correlations
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2} refN {-2}", "ChGap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2 2} refN {-2 -2}", "ChGap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {3} refN {-3}", "ChGap32", kFALSE));

    fGFW->CreateRegions();
  }

  // Fill the FlowContainer
  void FillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm, bool fillflag)
  {
    // Calculate the correlations from the GFW
    double dnx, dny, valx;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    dny = fGFW->Calculate(corrconf, 0, kTRUE).imag();
    if (dnx == 0) {
      return;
    }

    if (!corrconf.pTDif) {
      valx = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(valx) < 1) {
        fFC->FillProfile(corrconf.Head.c_str(), cent, valx, 1, rndm);
        if (dny == 0) {
          return;
        }
      }
      return;
    }
    uint8_t nAxisPtBins = 31;
    for (int i = 1; i <= nAxisPtBins; i++) {
      dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
      if (dnx == 0) {
        return;
      }
      valx = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(valx) < 1) {
        // Fill the charged particle correlation vs pT profiles
        fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, valx, 1., rndm);
      }
      return;
    }
  }

  // Templated function instantianed for all of the process functions
  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runFillQvector(TEvent const& collision, aod::BCsWithTimestamps const&, TTracks const& tracks1)
  {
    // Fill the event properties within the VarManager
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<TEventFillMap>(collision);

    // TODO: bc that could be used later to get timestamp for acceptance/GFWWeights
    // auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    // if (fConfigAcceptance.value.empty() == false) { cfg.mAcceptance = ccdb->getForTimeStamp<GFWWeights>(fConfigAcceptance.value, bc.timestamp());}

    fGFW->Clear();

    constexpr bool eventHasCentRun2 = ((TEventFillMap & VarManager::ObjTypes::CollisionCentRun2) > 0);
    constexpr bool eventHasCentRun3 = ((TEventFillMap & VarManager::ObjTypes::CollisionCent) > 0);

    // Acceptance and efficiency weights
    float weff = 1.0, wacc = 1.0;

    // Fill the GFW object in the track loop
    for (auto& track : tracks1) {
      if (cfg.mEfficiency) {
        weff = cfg.mEfficiency->GetBinContent(cfg.mEfficiency->FindBin(track.pt()));
      } else {
        weff = 1.0;
      }
      if (weff == 0) {
        continue;
      }
      weff = 1. / weff;
      if (cfg.mAcceptance) {
        wacc = cfg.mAcceptance->GetNUA(track.phi(), track.eta(), collision.posZ());
      } else {
        wacc = 1.0;
      }
      // Fill the GFW for each track to compute Q vector and correction using weights
      fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, 3); // using default values for ptin = 0 and mask = 3
    }

    float l_Random = fRndm->Rndm(); // used only to compute correlators
    bool fillFlag = kFALSE;         // could be used later
    for (uint64_t l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      if constexpr (eventHasCentRun2) {
        FillFC(corrconfigs.at(l_ind), VarManager::fgValues[VarManager::kCentVZERO], l_Random, fillFlag);
      }
      if constexpr (eventHasCentRun3) {
        FillFC(corrconfigs.at(l_ind), VarManager::fgValues[VarManager::kCentFT0C], l_Random, fillFlag);
      }
    }

    // Define quantities needed for the different eta regions
    uint8_t nentriesN = 0.0;
    uint8_t nentriesP = 0.0;
    uint8_t nentriesFull = 0.0;
    complex<double> Q2vecN;
    complex<double> Q2vecP;
    complex<double> Q2vecFull;
    complex<double> Q3vecN;
    complex<double> Q3vecP;
    complex<double> Q3vecFull;

    if (fGFW && (tracks1.size() > 0)) {
      // Obtain the GFWCumulant where Q is calculated (index=region, with different eta gaps)
      GFWCumulant gfwCumN = fGFW->GetCumulant(0);
      GFWCumulant gfwCumP = fGFW->GetCumulant(1);
      GFWCumulant gfwCumFull = fGFW->GetCumulant(2);

      // and the multiplicity of the event in each region
      nentriesN = gfwCumN.GetN();
      nentriesP = gfwCumP.GetN();
      nentriesFull = gfwCumFull.GetN();

      // Get the Q vector for selected harmonic, power (for minPt=0)
      Q2vecN = gfwCumN.Vec(2, fConfigNPow);
      Q2vecP = gfwCumP.Vec(2, fConfigNPow);
      Q2vecFull = gfwCumFull.Vec(2, fConfigNPow);
      Q3vecN = gfwCumN.Vec(3, fConfigNPow);
      Q3vecP = gfwCumP.Vec(3, fConfigNPow);
      Q3vecFull = gfwCumFull.Vec(3, fConfigNPow);
    }

    // Fill the VarManager::fgValues with the Q vector quantities
    VarManager::FillQVectorFromGFW(collision, Q2vecFull, Q2vecN, Q2vecP, Q3vecFull, Q3vecN, Q3vecP, nentriesFull, nentriesN, nentriesP);

    if (fConfigQA) {
      if ((tracks1.size() > 0) && (nentriesFull * nentriesN * nentriesP != 0.0)) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
        if (fEventCut->IsSelected(VarManager::fgValues)) {
          fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
        }
      }
    }

    // Fill the tree for the reduced event table with Q vector quantities
    if (fEventCut->IsSelected(VarManager::fgValues)) {
      eventQvector(VarManager::fgValues[VarManager::kQ2X0A], VarManager::fgValues[VarManager::kQ2Y0A], VarManager::fgValues[VarManager::kQ2X0B], VarManager::fgValues[VarManager::kQ2Y0B], VarManager::fgValues[VarManager::kQ2X0C], VarManager::fgValues[VarManager::kQ2Y0C], VarManager::fgValues[VarManager::kMultA], VarManager::fgValues[VarManager::kMultC], VarManager::fgValues[VarManager::kMultC], VarManager::fgValues[VarManager::kQ3X0A], VarManager::fgValues[VarManager::kQ3Y0A], VarManager::fgValues[VarManager::kQ3X0B], VarManager::fgValues[VarManager::kQ3Y0B], VarManager::fgValues[VarManager::kQ3X0C], VarManager::fgValues[VarManager::kQ3Y0C]);
    }
  }

  // Process to fill Q vector using barrel tracks in a reduced event table for barrel/muon tracks flow related analyses Run 2
  void processBarrelQvectorRun2(MyEventsWithCent::iterator const& collisions, aod::BCsWithTimestamps const& bcs, soa::Filtered<MyBarrelTracks> const& tracks)
  {
    runFillQvector<gkEventFillMap, gkTrackFillMap>(collisions, bcs, tracks);
  }

  // Process to fill Q vector using barrel tracks in a reduced event table for barrel/muon tracks flow related analyses Run 3
  void processBarrelQvector(MyEventsWithCentRun3::iterator const& collisions, aod::BCsWithTimestamps const& bcs, soa::Filtered<MyBarrelTracks> const& tracks)
  {
    runFillQvector<gkEventFillMapRun3, gkTrackFillMap>(collisions, bcs, tracks);
  }

  // Process to fill Q vector using forward tracks in a reduced event table for barrel/muon tracks flow related analyses Run 3
  void processForwardQvector(MyEventsWithCentRun3::iterator const& collisions, aod::BCsWithTimestamps const& bcs, soa::Filtered<aod::MFTTracks> const& tracks)
  {
    runFillQvector<gkEventFillMapRun3, 0u>(collisions, bcs, tracks);
  }

  // TODO: dummy function for the case when no process function is enabled
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQEventQvector, processBarrelQvectorRun2, "Run q-vector task on barrel tracks for Run2", false);
  PROCESS_SWITCH(DQEventQvector, processBarrelQvector, "Run q-vector task on barrel tracks for Run3", false);
  PROCESS_SWITCH(DQEventQvector, processForwardQvector, "Run q-vector task on forward tracks for Run3", false);
  PROCESS_SWITCH(DQEventQvector, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DQEventQvector>(cfgc)};
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
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "qvector,trigger,cent,res");
    }
  }
}
