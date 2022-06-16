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
///       Task to compute Q vectors and other quanitites related from the generic framework
///       Generic framework O2 version is a port of the AliPhysics version
///       To be used in the DQ analyses aiming for flow measurements
///       Run the standalone task with:
///       o2-analysis-timestamp --aod-file AO2D.root -b | o2-analysis-event-selection -b | o2-analysis-multiplicity-table -b | o2-analysis-centrality-table -b | o2-analysis-fdd-converter -b | o2-analysis-trackselection -b | o2-analysis-trackextension -b | o2-analysis-pid-tpc-full -b | o2-analysis-pid-tof-full -b | o2-analysis-pid-tof-base -b | o2-analysis-pid-tof-beta -b | o2-analysis-dq-flow -b
///       tested (June 2, 2022) on AO2D.root files from train production 242

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
#include "GFW.h"
#include "GFWCumulant.h"
#include "FlowContainer.h"
#include "GFWWeights.h"
#include <TRandom3.h>
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include <TH1F.h>
#include <THashList.h>
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

// Declarations of various short names
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

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;

using MyMuons = aod::FwdTracks;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

struct AnalysisQvector {
  Produces<ReducedEventsQvector> eventQvector;

  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "", "Comma separated list of muon cuts"};
  Configurable<bool> fConfigQA{"cfgQA", true, "If true, fill QA histograms"};

  // Configurable<float> fConfigVtxCut{"cfgVtxCut", 12.0, "Z vertex cut"};
  Configurable<float> fConfigCutPtMin{"cfgCutPtMin", 0.2f, "Minimal pT for tracks"};
  Configurable<float> fConfigCutPtMax{"cfgCutPtMax", 12.0f, "Maximal pT for tracks"};
  Configurable<float> fConfigCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> fConfigEtaLimit{"cfgEtaLimit", 0.4f, "Eta gap separation, only if using subEvents"};
  Configurable<int> fConfigNHarm{"cfgNHarm", 2, "Harmonic number of Q vector"};
  Configurable<int> fConfigNPow{"cfgNPow", 0, "Power of weights for Q vector"};

  // Access to the efficiencies and acceptances from CCDB
  Configurable<std::string> fConfigEfficiency{"cfgEfficiency", "", "CCDB path to efficiency object"};
  Configurable<std::string> fConfigAcceptance{"cfgAcceptance", "", "CCDB path to acceptance object"};
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> fConfigURL{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};
  Configurable<std::string> fConfigCCDBPath{"ccdb-path", "Users/lm", "base path to the ccdb object"};
  Configurable<long> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  // Configurables for FlowContainer (e.g charged particles pt-differential v22, v23, ...)
  //  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  //  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  //  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};
  //  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1}, "multiplicity / centrality axis for histograms"};
  //  AxisSpec axisCentBins{{0, 5., 10., 20., 30., 40., 50., 60., 70., 80.}, "centrality percentile"};

  // Filter collisionFilter = nabs(aod::collision::posZ) < fConfigVtxCut;
  Filter trackFilter = (nabs(aod::track::eta) < fConfigCutEta) && (aod::track::pt > fConfigCutPtMin) && (aod::track::pt < fConfigCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  struct Config {
    TH1D* mEfficiency = nullptr;
    GFWWeights* mAcceptance = nullptr;
  } cfg;

  // Define output
  HistogramManager* fHistMan = nullptr;
  AnalysisCompositeCut* fEventCut;
  OutputObj<THashList> fOutputList{"outputQA"};
  // OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};  // Need to add a dictionary for FlowContainer output

  // define global variables for generic framework
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
    auto histCCDB = ccdb->get<TH1F>(fConfigCCDBPath.value);
    if (!histCCDB) {
      LOGF(fatal, "CCDB histogram not found");
    }

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
      // cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(fConfigEfficiency.value, fConfigNoLaterThan.value);
      if (cfg.mEfficiency)
        LOGF(info, "Loaded efficiency histogram %s (%p)", fConfigEfficiency.value.c_str(), (void*)cfg.mEfficiency);
      else
        LOGF(info, "Could not load efficiency histogram from %s (%p)", fConfigEfficiency.value.c_str(), (void*)cfg.mEfficiency);
    }

    // Reference flow
    //    TObjArray* oba = new TObjArray();
    //    oba->Add(new TNamed("ChGap22", "ChGap22"));   // for gap (|eta|>0.4) case
    //    oba->Add(new TNamed("ChGap24", "ChGap24"));   // for gap (|eta|>0.4) case
    //    oba->Add(new TNamed("ChFull22", "ChFull22")); // no-gap case
    //    oba->Add(new TNamed("ChFull24", "ChFull24")); // no-gap case
    //    oba->Add(new TNamed("ChGap32", "ChGap32"));   // gap-case
    // fFC->SetName("FlowContainer");
    // fFC->Initialize(oba, axisMultiplicity, 10);
    // delete oba;

    int pows[] = {3, 0, 2, 2, 3, 3, 3};
    int powsFull[] = {5, 0, 4, 4, 3, 3, 3};
    // Define regions of positive and negative eta in order to create gaps
    fGFW->AddRegion("refN", 7, pows, -fConfigCutEta, -fConfigEtaLimit, 1, 1);
    fGFW->AddRegion("refP", 7, pows, fConfigEtaLimit, fConfigCutEta, 1, 1);
    fGFW->AddRegion("full", 7, powsFull, -fConfigCutEta, fConfigCutEta, 1, 2);

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2} refN {-2}", "ChGap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2 2} refN {-2 -2}", "ChGap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {3} refN {-3}", "ChGap32", kFALSE));
  }

  // TODO: make available the flowcontainer output (add a dictionary somewhere...)
  // Fill the FlowContainer
  void FillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm, bool fillflag, bool dqflag)
  {
    // Calculate the correlations from the GFW
    //    double dnx, dny, valx, valy;
    //    dnx = fGFW->Calculate(corrconf, 0, kTRUE).Re();
    //    dny = fGFW->Calculate(corrconf, 0, kTRUE).Im();
    //    if (dnx == 0) {
    //      return;
    //    }
    //
    //    if (!corrconf.pTDif) {
    //      valx = fGFW->Calculate(corrconf, 0, kFALSE).Re() / dnx;
    //      if (TMath::Abs(valx) < 1) {
    //        fFC->FillProfile(corrconf.Head.Data(), cent, valx, 1, rndm);
    //        if (dny == 0) {
    //          return;
    //        }
    //        valy = fGFW->Calculate(corrconf, 0, kFALSE).Re() / dny;
    //      }
    //      return;
    //    }
    //    bool DisableOverlap = kFALSE;
    //    int nAxisPtBins = 31;
    //    for (int i = 1; i <= nAxisPtBins; i++) {
    //      dnx = fGFW->Calculate(corrconf, 0, kTRUE, DisableOverlap).Re();
    //      if (dnx == 0) {
    //        return;
    //      }
    //      valx = fGFW->Calculate(corrconf, 0, kFALSE, DisableOverlap).Re() / dnx;
    //      if (TMath::Abs(valx) < 1) {
    //         fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.Data(), i), cent, valx, 1., rndm);
    //      }
    //      return;
    //    }
  }

  // Templated function instantianed for all of the process functions
  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runFillQvector(TEvent const& collision, aod::BCs const& bcs, TTracks const& tracks1)
  {
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<TEventFillMap>(collision);

    // TODO: properly access to config files from ccdb using bc.timestamp()
    // auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (fConfigAcceptance.value.empty() == false) {
      // cfg.mAcceptance = ccdb->getForTimeStamp<GFWWeights>(fConfigAcceptance.value, bc.timestamp());
      if (cfg.mAcceptance) {
        LOGF(info, "Loaded acceptance histogram from %s (%p)", fConfigAcceptance.value.c_str(), (void*)cfg.mAcceptance);
      } else {
        LOGF(warning, "Could not load acceptance histogram from %s (%p)", fConfigAcceptance.value.c_str(), (void*)cfg.mAcceptance);
      }
    }

    fGFW->Clear();

    // acceptance and efficiency weights
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
        wacc = 1;
      }
      // VarManager::FillTrack<TTrackFillMap>(track);

      // Fill the GFW for each track to compute Q vector
      fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, 3); // using default values for ptin=0 and mask=3
    }

    //    float l_Random = fRndm->Rndm(); // used only to compute correlators
    //    bool fillFlag = kFALSE;         // could be used later
    //    bool DQEventFlag = kFALSE;      // could be used later
    //    for (unsigned long int l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
    //      FillFC(corrconfigs.at(l_ind), collision.centRun2V0M(), l_Random, fillFlag, DQEventFlag);
    //    };

    // Obtain the GFWCumulant where Q is calculated (index=region, with different eta gaps)
    GFWCumulant gfwCumN = fGFW->GetCumulant(0);
    GFWCumulant gfwCumP = fGFW->GetCumulant(1);
    GFWCumulant gfwCum = fGFW->GetCumulant(2);

    // Get the multiplicity of the event in this region
    int nentriesN = gfwCumN.GetN();
    int nentriesP = gfwCumP.GetN();
    int nentries = gfwCum.GetN();

    // Get the Q vector for selected harmonic, power (for minPt=0)
    TComplex QvecN = gfwCumN.Vec(fConfigNHarm, fConfigNPow);
    TComplex QvecP = gfwCumP.Vec(fConfigNHarm, fConfigNPow);
    TComplex Qvec = gfwCum.Vec(fConfigNHarm, fConfigNPow);

    // TODO: move this part to later analysis development
    // compute the resolution from Q vectors estimated with different eta gaps
    //    Double_t resGap = 0.0;
    //    if (nentriesN > 0 && nentriesP > 0) {
    //      resGap = (QvecP.Re() * QvecN.Re() + QvecP.Im() * QvecN.Im()) / (nentriesN * nentriesP);
    //    }

    // Fill the VarManager::fgValues with the Q vector quantities
    VarManager::FillQVectorFromGFW(collision, Qvec, QvecN, QvecP, nentries, nentriesN, nentriesP);
    if (fConfigQA) {
      if (nentriesN * nentriesP * nentries != 0) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
        if (fEventCut->IsSelected(VarManager::fgValues)) {
          fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
        }
      }
    }

    // Fill the tree for the reduced event table with Q vector quantities
    if (fEventCut->IsSelected(VarManager::fgValues)) {
      eventQvector(VarManager::fgValues[VarManager::kQ2X0A], VarManager::fgValues[VarManager::kQ2Y0A], VarManager::fgValues[VarManager::kQ2X0B], VarManager::fgValues[VarManager::kQ2Y0B], VarManager::fgValues[VarManager::kQ2X0C], VarManager::fgValues[VarManager::kQ2Y0C], VarManager::fgValues[VarManager::kMultA], VarManager::fgValues[VarManager::kMultC], VarManager::fgValues[VarManager::kMultC]);
    }
  }

  // Process to fill Q vector in a reduced event table for barrel/muon tracks flow related analyses
  void process(MyEventsWithCent::iterator const& collisions, aod::BCs const& bcs, soa::Filtered<MyBarrelTracks> const& tracks)
  {
    runFillQvector<gkEventFillMap, gkTrackFillMap>(collisions, bcs, tracks);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisQvector>(cfgc)};
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
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "qvector,trigger,cent");
    }
  }
}
