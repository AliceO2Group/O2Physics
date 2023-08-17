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

#include <CCDB/BasicCCDBManager.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"

#include "GFWPowerArray.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "FlowContainer.h"
#include "GFWWeights.h"
#include <TProfile.h>
#include <TRandom3.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct GenericFramework {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 3.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "Number of chi2 per TPC cluster")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 0.1f, "Maximum DCA xy")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")

  ConfigurableAxis axisVertex{"axisVertex", {22, -11, 11}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1}, "multiplicity / centrality axis for histograms"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls) && (aod::track::dcaXY < cfgCutDCAxy);
  using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  struct Config {
    TH1D* mEfficiency = nullptr;
    GFWWeights* mAcceptance = nullptr;
    bool correctionsLoaded = false;
  } cfg;

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TRandom3* fRndm = new TRandom3(0);
  TAxis* fPtAxis;

  void init(InitContext const&)
  {

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* PtBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, PtBins);

    if (doprocessWeights) {
      fWeights->SetPtBins(nPtBins, PtBins);
      fWeights->Init(true, false);
    }

    if (doprocessData) {
      registry.add("hPhi", "", {HistType::kTH1D, {axisPhi}});
      registry.add("hEta", "", {HistType::kTH1D, {axisEta}});
      registry.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
      registry.add("hPhiEtaVtxZ_corrected", "", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      registry.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});

      TObjArray* oba = new TObjArray();
      // Reference flow
      oba->Add(new TNamed("ChGap22", "ChGap22"));     // for gap (|eta|>0.4) case
      oba->Add(new TNamed("ChGap32", "ChGap32"));     //
      oba->Add(new TNamed("ChGap42", "ChGap42"));     //
      oba->Add(new TNamed("ChFull22", "ChFull22"));   // no-gap case
      oba->Add(new TNamed("ChFull24", "ChFull24"));   // no-gap case
      oba->Add(new TNamed("ChFull26", "ChFull26"));   // no-gap case
      oba->Add(new TNamed("ChFull28", "ChFull28"));   // no-gap case
      oba->Add(new TNamed("ChFull210", "ChFull210")); // no-gap case
      fFC->SetName("FlowContainer");
      fFC->Initialize(oba, axisMultiplicity, cfgNbootstrap);
      delete oba;

      fGFW->AddRegion("refN", -0.8, -0.4, 1, 1);
      fGFW->AddRegion("refP", 0.4, 0.8, 1, 1);
      fGFW->AddRegion("full", -0.8, 0.8, 1, 2);
      CreateCorrConfigs();
      fGFW->CreateRegions();
    }
  }

  void loadCorrections(uint64_t timestamp)
  {
    if (cfg.correctionsLoaded)
      return;
    if (cfgAcceptance.value.empty() == false) {
      cfg.mAcceptance = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance, timestamp);
      if (cfg.mAcceptance)
        LOGF(info, "Loaded acceptance weights from %s (%p)", cfgAcceptance.value.c_str(), (void*)cfg.mAcceptance);
      else
        LOGF(warning, "Could not load acceptance weights from %s (%p)", cfgAcceptance.value.c_str(), (void*)cfg.mAcceptance);
    }
    if (cfgEfficiency.value.empty() == false) {
      cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgEfficiency, timestamp);
      if (cfg.mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)cfg.mEfficiency);
    }
    cfg.correctionsLoaded = true;
  }

  void CreateCorrConfigs()
  {
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2} refN {-2}", "ChGap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {3} refN {-3}", "ChGap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {4} refN {-4}", "ChGap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 2 -2 -2 -2}", "ChFull26", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 2 2 -2 -2 -2 -2}", "ChFull28", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 2 2 2 -2 -2 -2 -2 -2}", "ChFull210", kFALSE));
  }

  void FillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        fFC->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      return;
    }
    return;
  }

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>::iterator const& collision, aod::BCsWithTimestamps const&, myTracks const& tracks)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    loadCorrections(bc.timestamp());

    if (tracks.size() < 1)
      return;

    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);

    fGFW->Clear();
    const auto centrality = collision.centFT0C();
    registry.fill(HIST("hCent"), centrality);
    if (centrality > 100)
      return;
    float l_Random = fRndm->Rndm();
    float weff = 1, wacc = 1;

    for (auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < 70)
        continue;
      if (track.dcaXY() > 7 * (0.0026 + 0.0050 / pow(track.pt(), 1.01)))
        continue;
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hEta"), track.eta());

      if (cfg.mEfficiency)
        weff = cfg.mEfficiency->GetBinContent(cfg.mEfficiency->FindBin(track.pt()));
      else
        weff = 1.0;
      if (weff == 0)
        continue;
      weff = 1. / weff;
      if (cfg.mAcceptance)
        wacc = cfg.mAcceptance->GetNUA(track.phi(), track.eta(), vtxz);
      else
        wacc = 1;
      registry.fill(HIST("hPhiEtaVtxZ_corrected"), track.phi(), track.eta(), vtxz, wacc);
      fGFW->Fill(track.eta(), 1, track.phi(), wacc * weff, 3);
    }
    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      FillFC(corrconfigs.at(l_ind), centrality, l_Random);
    }
  }
  PROCESS_SWITCH(GenericFramework, processData, "Process analysis for data", true);

  void processWeights(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>::iterator const& collision, myTracks const& tracks)
  {
    const auto centrality = collision.centFT0C();
    if (centrality > 100)
      return;
    for (auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < 70)
        continue;
      if (track.dcaXY() > 7 * (0.0026 + 0.0050 / pow(track.pt(), 1.01)))
        continue;
      fWeights->Fill(track.phi(), track.eta(), collision.posZ(), track.pt(), centrality, 0);
    }
  }
  PROCESS_SWITCH(GenericFramework, processWeights, "Process weights for acceptance corrections", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<GenericFramework>(cfgc),
  };
}
