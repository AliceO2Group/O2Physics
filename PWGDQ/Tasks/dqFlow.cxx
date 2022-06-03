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
using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;

using MyBarrelTracksWithCov = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksExtended, aod::TrackSelection,
                                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                        aod::pidTPCFullKa, aod::pidTPCFullPr,
                                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                        aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;

using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsWithFilter = soa::Join<aod::Collisions, aod::EvSels, aod::DQEventFilter>;
using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;

using MyMuons = aod::FwdTracks;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

struct AnalysisQvector {
  Produces<ReducedEventsQvector> eventQvector;

  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "", "Comma separated list of muon cuts"};
  Configurable<bool> fConfigQA{"cfgQA", true, "If true, fill QA histograms"};

  Configurable<float> fVtxCut{"VtxCut", 12.0, "Z vertex cut"};
  Configurable<float> cfgCutPtMin{"cfgCutPtMin", 0.2f, "Minimal pT for tracks"};
  Configurable<float> cfgCutPtMax{"cfgCutPtMax", 12.0f, "Maximal pT for tracks"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgEtaLimit{"cfgEtaLimit", 0.4f, "Eta gap separation, only if subEvents=true"};
  // Configurable<bool> bUseWeights{"UseWeights", true, "If true, fill Q vectors with weights for phi and p_T"};
  // Configurable<bool> bSubEvents{"SubEvents", true, "If true, fill use sub-events methods with different detector gaps"};
  Configurable<int> nHarm{"nHarm", 2, "Harmonic number of Q vector"};
  Configurable<int> nPow{"nPow", 0, "Power of weights for Q vector"};

  // Access to the efficiencies and acceptances from CCDB
  Configurable<std::string> fcfgEfficiency{"cfgEfficiency", "", "CCDB path to efficiency object"};
  Configurable<std::string> fcfgAcceptance{"cfgAcceptance", "", "CCDB path to acceptance object"};
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdb-path", "Users/lm", "base path to the ccdb object"};
  Configurable<long> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  // Configurables for FlowContainer (e.g charged particles pt-differential v22, v23, ...)
  //  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  //  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  //  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};
  //  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1}, "multiplicity / centrality axis for histograms"};
  //  AxisSpec axisCentBins{{0, 5., 10., 20., 30., 40., 50., 60., 70., 80.}, "centrality percentile"};

  Filter collisionFilter = nabs(aod::collision::posZ) < fVtxCut;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  struct Config {
    TH1D* mEfficiency = nullptr;
    GFWWeights* mAcceptance = nullptr;
  } cfg;

  // Define output
  HistogramManager* fHistMan = nullptr;
  AnalysisCompositeCut* fEventCut;
  std::vector<TString> fTrackHistNames;
  OutputObj<THashList> fOutputList{"output"};
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

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(nolaterthan.value);
    auto histCCDB = ccdb->get<TH1F>(ccdbPath.value);
    if (!histCCDB) {
      LOGF(fatal, "CCDB histogram not found");
    }

    TString histNames = "";
    VarManager::SetDefaultVarNames();

    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      std::vector<TString> names = {
        "Event_BeforeCuts",
        "Event_AfterCuts"};
      histNames += Form("%s;%s", names[0].Data(), names[1].Data());
      for (int i = 0; i < 2; i++) {
        fTrackHistNames.push_back(names[i]);
      }

      DefineHistograms(fHistMan, histNames.Data());    // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    // Global effiencies
    if (fcfgEfficiency.value.empty() == false) {
      // cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(fcfgEfficiency.value, nolaterthan.value);
      if (cfg.mEfficiency)
        LOGF(info, "Loaded efficiency histogram %s (%p)", fcfgEfficiency.value.c_str(), (void*)cfg.mEfficiency);
      else
        LOGF(info, "Could not load efficiency histogram from %s (%p)", fcfgEfficiency.value.c_str(), (void*)cfg.mEfficiency);
    }

    // Reference flow
    TObjArray* oba = new TObjArray();
    oba->Add(new TNamed("ChGap22", "ChGap22"));   // for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap24", "ChGap24"));   // for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChFull22", "ChFull22")); // no-gap case
    oba->Add(new TNamed("ChFull24", "ChFull24")); // no-gap case
    oba->Add(new TNamed("ChGap32", "ChGap32"));   // gap-case

    // fFC->SetName("FlowContainer");
    // fFC->Initialize(oba, axisMultiplicity, 10);
    delete oba;

    int pows[] = {3, 0, 2, 2, 3, 3, 3};
    int powsFull[] = {5, 0, 4, 4, 3, 3, 3};
    // Define regions of positive and negative eta in order to create gaps
    fGFW->AddRegion("refN", 7, pows, -cfgCutEta, -cfgEtaLimit, 1, 1);
    fGFW->AddRegion("refP", 7, pows, cfgEtaLimit, cfgCutEta, 1, 1);
    fGFW->AddRegion("full", 7, powsFull, -cfgCutEta, cfgCutEta, 1, 2);

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2} refN {-2}", "ChGap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2 2} refN {-2 -2}", "ChGap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {3} refN {-3}", "ChGap32", kFALSE));
  }

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
    // Reset the fValues and Qn vector array
    VarManager::ResetValues(0, VarManager::kNVars);
    // Fill the event information with the VarManager
    VarManager::FillEvent<TEventFillMap>(collision);

    uint64_t tag = 0;
    // store the selection decisions
    for (int i = 0; i < kNsel; i++) {
      if (collision.selection()[i] > 0) {
        tag |= (uint64_t(1) << i);
      }
    }
    if (collision.sel7()) {
      tag |= (uint64_t(1) << kNsel); //! SEL7 stored at position kNsel in the tag bit map
    }
    int mult = tracks1.size();

    // TODO: properly access to config files from ccdb using bc.timestamp()
    // auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (fcfgAcceptance.value.empty() == false) {
      // cfg.mAcceptance = ccdb->getForTimeStamp<GFWWeights>(fcfgAcceptance.value, bc.timestamp());
      if (cfg.mAcceptance) {
        LOGF(info, "Loaded acceptance histogram from %s (%p)", fcfgAcceptance.value.c_str(), (void*)cfg.mAcceptance);
      } else {
        LOGF(warning, "Could not load acceptance histogram from %s (%p)", fcfgAcceptance.value.c_str(), (void*)cfg.mAcceptance);
      }
    }
    if (mult < 1) {
      return;
    }
    if (!collision.sel7()) {
      return;
    }
    float vtxz = collision.posZ();

    fGFW->Clear();
    const auto centrality = collision.centRun2V0M();
    if (centrality > 100) {
      return;
    }

    float l_Random = fRndm->Rndm(); // used only to compute correlators
    float weff = 1.0, wacc = 1.0;   // acceptance and efficiency weights

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
        wacc = cfg.mAcceptance->GetNUA(track.phi(), track.eta(), vtxz);
      } else {
        wacc = 1;
      }
      // Fill the track information in the VarManager
      VarManager::FillTrack<TTrackFillMap>(track);

      int ptin = 0; // default value
      int mask = 3; // default value
      // Fill the GFW for each track in order to compute Q vector
      fGFW->Fill(track.eta(), ptin, track.phi(), wacc * weff, mask);
    }

    bool fillFlag = kFALSE;    // could be used later
    bool DQEventFlag = kFALSE; // could be used later
    for (unsigned long int l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      FillFC(corrconfigs.at(l_ind), centrality, l_Random, fillFlag, DQEventFlag);
    };

    // Obtain the GFWCumulant where Q is calculated (index=region, with different eta gaps)
    GFWCumulant fGFWCumulantN = fGFW->GetCumulant(0);
    GFWCumulant fGFWCumulantP = fGFW->GetCumulant(1);
    GFWCumulant fGFWCumulant = fGFW->GetCumulant(2);

    // Get the multiplicity of the event in this region
    int nentriesN = fGFWCumulant.GetN();
    int nentriesP = fGFWCumulant.GetN();
    int nentries = fGFWCumulant.GetN();

    // Get the Q vector for selected harmonic, power (for minPt=0)
    TComplex QvecN = fGFWCumulantN.Vec(nHarm, nPow);
    TComplex QvecP = fGFWCumulantP.Vec(nHarm, nPow);
    TComplex Qvec = fGFWCumulant.Vec(nHarm, nPow);

    Double_t resGap = 0.0;

    if (nentriesN > 0 && nentriesP > 0) {
      // TODO: provide other calculation of R
      // compute the resolution from Q vectors estimated with different eta gaps
      resGap = (QvecP.Re() * QvecN.Re() + QvecP.Im() * QvecN.Im()) / (nentriesN * nentriesP);

      // Fill the VarManager::fgValues with the Q vector quantities
      VarManager::FillQVectorFromGFW<TEventFillMap>(collision, Qvec, resGap, nentries);
      if (fConfigQA) {
        fHistMan->FillHistClass(fTrackHistNames[0].Data(), VarManager::fgValues);
        if (fEventCut->IsSelected(VarManager::fgValues)) {
          fHistMan->FillHistClass(fTrackHistNames[1].Data(), VarManager::fgValues);
        }
      }
    }

    // Fill the tree for the reduced event table with Q vector quantities
    eventQvector(tag, collision.bc().runNumber(), collision.posZ(), VarManager::fgValues[VarManager::kCentVZERO], VarManager::fgValues[VarManager::kQ2X0], VarManager::fgValues[VarManager::kQ2Y0], VarManager::fgValues[VarManager::kPsi2], VarManager::fgValues[VarManager::kRes]);
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
