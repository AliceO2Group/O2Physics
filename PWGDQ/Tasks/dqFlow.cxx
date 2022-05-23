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
using MyCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>;

using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsWithFilter = soa::Join<aod::Collisions, aod::EvSels, aod::DQEventFilter>;
using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;

using MyMuons = aod::FwdTracks;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;

struct AnalysisQvector {

  // TODO: Provide an access to the Q vector for basics dilepton flow related analyses (to be adapted)

  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "", "Comma separated list of muon cuts"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};

  Configurable<float> fVtxCut{"VtxCut", 14.0, "Z vertex cut"};
  Configurable<bool> bUseWeights{"UseWeights", true, "If true, fill Q vectors with weights for phi and p_T"};
  Configurable<bool> bSubEvents{"SubEvents", true, "If true, fill use sub-events methods with different detector gaps"};
  Configurable<float> fEtaLimit{"EtaLimit", 0.0, "Eta gap separation (e.g ITS=0.0, MFT=-3.05,...), only if subEvents=true"};
  Configurable<int> nHarm{"nHarm", 2, "Harmonic number of Q vector"};
  Configurable<int> nPow{"nPow", 0, "Power of weights for Q vector"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1}, "multiplicity / centrality axis for histograms"};

  // TODO: Obtain efficiencies and acceptances from CCDB
  Configurable<std::string> fcfgEfficiency{"cfgEfficiency", "", "CCDB path to efficiency object"};
  Configurable<std::string> fcfgAcceptance{"cfgAcceptance", "", "CCDB path to acceptance object"};
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdb-path", "Users/lm", "base path to the ccdb object"};
  Configurable<long> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  Filter collisionFilter = nabs(aod::collision::posZ) < fVtxCut;

  HistogramManager* fHistMan = nullptr;
  AnalysisCompositeCut* fEventCut;

  struct Config {
    TH1D* mEfficiency = nullptr;
    GFWWeights* mAcceptance = nullptr;
  } cfg;

  // Define output
  OutputObj<THashList> fOutputList{"output"};
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TRandom3* fRndm = new TRandom3(0);

  // Initialize CCDB, efficiencies and acceptances from CCDB, histograms, GFW, FlowContainer
  // in order to fill Q vector with or without corrections
  void init(o2::framework::InitContext&)
  {
    // TODO: initialize

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    // Global effiencies
    if (fcfgEfficiency.value.empty() == false) {
      cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(fcfgEfficiency.value, nolaterthan.value);
      if (cfg.mEfficiency)
        LOGF(info, "Loaded efficiency histogram %s (%p)", fcfgEfficiency.value.c_str(), (void*)cfg.mEfficiency);
      else
        LOGF(info, "Could not load efficiency histogram from %s (%p)", fcfgEfficiency.value.c_str(), (void*)cfg.mEfficiency);
    }

    TObjArray* oba = new TObjArray();
    // Reference flow
    oba->Add(new TNamed("ChGap22", "ChGap22"));   // for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap24", "ChGap24"));   // for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChFull22", "ChFull22")); // no-gap case
    oba->Add(new TNamed("ChFull24", "ChFull24")); // no-gap case
    oba->Add(new TNamed("ChGap32", "ChGap32"));   // gap-case
    fFC->SetName("FlowContainer");
    fFC->Initialize(oba, axisMultiplicity, 10);
    delete oba;

    int pows[] = {3, 0, 2, 2, 3, 3, 3};
    int powsFull[] = {5, 0, 4, 4, 3, 3, 3};
    fGFW->AddRegion("refN", 7, pows, -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP", 7, pows, 0.4, 0.8, 1, 1);
    fGFW->AddRegion("full", 7, powsFull, -0.8, 0.8, 1, 2);

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2} refN {-2}", "ChGap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2 2} refN {-2 -2}", "ChGap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {3} refN {-3}", "ChGap32", kFALSE));
  }

  void FillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    // Calculate the correlations from the GFW
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).Re();
    if (dnx == 0) {
      return;
    }
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).Re() / dnx;
      if (TMath::Abs(val) < 1) {
        fFC->FillProfile(corrconf.Head.Data(), cent, val, 1, rndm);
      }
      return;
    }
  }

  // Template function to run fill Q vector (alltracks-barrel, alltracks-muon)
  //  template <typename TCollision, typename BC, typename TTracks1, typename TTracks2 >
  //  void runFillQvector(TCollision const& collision, BC const& bcs, TTracks1 const& tracks1)
  void process(MyCollisions::iterator const& collision, aod::BCsWithTimestamps const& bcs, aod::Tracks const& tracks1, MyBarrelTracks const& barreltracks)
  {
    // Reset the fValues and Qn vector array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::ResetQvector();

    // TODO: implement main functions to fill Q vector with corrections for selected tracks

    // auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (fcfgAcceptance.value.empty() == false) {
      // cfg.mAcceptance = ccdb->getForTimeStamp<GFWWeights>(fcfgAcceptance.value, bc.timestamp());
      if (cfg.mAcceptance) {
        LOGF(info, "Loaded acceptance histogram from %s (%p)", fcfgAcceptance.value.c_str(), (void*)cfg.mAcceptance);
      } else {
        LOGF(warning, "Could not load acceptance histogram from %s (%p)", fcfgAcceptance.value.c_str(), (void*)cfg.mAcceptance);
      }
    }
    if (tracks1.size() < 1) {
      return;
    }
    if (!collision.sel7()) {
      return;
    }
    // LOGF(info, "Tracks for collision: %d | Vertex: %.1f | INT7: %d | V0M: %.1f", tracks.size(), collision.posZ(), collision.sel7(), collision.centV0M());
    float vtxz = collision.posZ();

    fGFW->Clear();
    const auto centrality = collision.centRun2V0M();
    if (centrality > 100) {
      return;
    }

    float l_Random = fRndm->Rndm();
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
        wacc = cfg.mAcceptance->GetNUA(track.phi(), track.eta(), vtxz);
      } else {
        wacc = 1;
      }
      fGFW->Fill(track.eta(), 1, track.phi(), wacc * weff, 3);
    }
    for (unsigned long int l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      FillFC(corrconfigs.at(l_ind), centrality, l_Random);
    };
  }

  //  // Process Q vector for Barrel tracks related analyses
  //  void processQvectorBarrelSkimmed(MyCollisions::iterator const& collision, aod::BCsWithTimestamps const& bcs, aod::Tracks const& tracks, MyBarrelTracks const& barreltracks)
  //  {
  //
  //    // TODO: fill Q vector with corrections on selected tracks for barrel analysis
  //    runFillQvector(collision, bcs, tracks, barreltracks);
  //  }
  //
  //  // Process Q vector for Muon tracks related analyses
  //  void processQvectorMuonsSkimmed(MyCollisions::iterator const& collision, aod::BCsWithTimestamps const& bcs, aod::Tracks const& tracks, MyMuons const& muons)
  //  {
  //
  //    // TODO: fill Q vector with corrections on selected tracks for muon analysis
  //    runFillQvector(collision, bcs, tracks, muons);
  //  }

  // Dummy function for the case when no process function is enabled
  //  void processDummy()
  //  {
  //    // do nothing
  //  }

  // PROCESS_SWITCH(AnalysisQvector, processQvectorBarrelSkimmed, "Fill Q vectors for selected events and tracks, for ee flow analyses", false);
  // PROCESS_SWITCH(AnalysisQvector, processQvectorMuonsSkimmed, "Fill Q vectors for selected events and tracks, for mumu flow analyses", false);
  // PROCESS_SWITCH(AnalysisQvector, processDummy, "Dummy function, enabled only if none of the others are enabled", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisQvector>(cfgc)};
}
