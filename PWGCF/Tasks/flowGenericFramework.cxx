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
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMin, float, 0.2f, "Minimal pT for reference tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMax, float, 3.0f, "Maximal pT for reference tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "Number of chi2 per TPC cluster")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 0.1f, "Maximum DCA xy")
  O2_DEFINE_CONFIGURABLE(cfgEtaSep, float, 0.4f, "Eta gap for flow calculations")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")

  ConfigurableAxis axisVertex{"axisVertex", {22, -11, 11}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1}, "multiplicity / centrality axis for histograms"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls) && (aod::track::dcaXY < cfgCutDCAxy);
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
      registry.add("hPt", "", {HistType::kTH1D, {axisPt}});
      registry.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
      registry.add("hPhiEtaVtxZ_corrected", "", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      registry.add("hCent", "", {HistType::kTH1D, {axisMultiplicity}});

      TObjArray* oba = new TObjArray();
      // Reference flow
      oba->Add(new TNamed("ChGapP22", "ChGapP22")); // for positive gap case
      for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("ChGapP22_pt_%i", i + 1), "ChGapP22_pTDiff"));
      oba->Add(new TNamed("ChGapP32", "ChGapP32"));
      for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("ChGapP32_pt_%i", i + 1), "ChGapP32_pTDiff"));
      oba->Add(new TNamed("ChGapP42", "ChGapP42"));
      for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("ChGapP42_pt_%i", i + 1), "ChGapP42_pTDiff"));

      oba->Add(new TNamed("ChGapN22", "ChGapN22")); // for negative gap case
      for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("ChGapN22_pt_%i", i + 1), "ChGapN22_pTDiff"));
      oba->Add(new TNamed("ChGapN32", "ChGapN32"));
      for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("ChGapN32_pt_%i", i + 1), "ChGapN32_pTDiff"));
      oba->Add(new TNamed("ChGapN42", "ChGapN42"));
      for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("ChGapN42_pt_%i", i + 1), "ChGapN42_pTDiff"));

      oba->Add(new TNamed("ChFull22", "ChFull22")); // no-gap case
      for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("ChFull22_pt_%i", i + 1), "ChFull22_pTDiff"));
      oba->Add(new TNamed("ChFull24", "ChFull24")); // no-gap case
      for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("ChFull24_pt_%i", i + 1), "ChFull24_pTDiff"));
      oba->Add(new TNamed("ChFull26", "ChFull26"));   // no-gap case
      oba->Add(new TNamed("ChFull28", "ChFull28"));   // no-gap case
      oba->Add(new TNamed("ChFull210", "ChFull210")); // no-gap case
      fFC->SetName("FlowContainer");
      fFC->SetXAxis(fPtAxis);
      fFC->Initialize(oba, axisMultiplicity, cfgNbootstrap);
      delete oba;

      fGFW->AddRegion("refN", -cfgCutEta, cfgEtaSep, 1, 1);
      fGFW->AddRegion("refP", cfgEtaSep, cfgCutEta, 1, 1);
      fGFW->AddRegion("refFull", -cfgCutEta, cfgCutEta, 1, 1);

      fGFW->AddRegion("poiN", -cfgCutEta, -cfgEtaSep, nPtBins + 1, 2);
      fGFW->AddRegion("poiP", cfgEtaSep, cfgCutEta, nPtBins + 1, 2);
      fGFW->AddRegion("poiFull", -cfgCutEta, cfgCutEta, nPtBins + 1, 2);

      fGFW->AddRegion("olN", -cfgCutEta, -cfgEtaSep, nPtBins + 1, 4);
      fGFW->AddRegion("olP", cfgEtaSep, cfgCutEta, nPtBins + 1, 4);
      fGFW->AddRegion("olFull", -cfgCutEta, cfgCutEta, nPtBins + 1, 4);

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
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2} refN {-2}", "ChGapP22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiP refP | olP {2} refN {-2}", "ChGapP22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN {2} refP {-2}", "ChGapN22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN refN | olN {2} refP {-2}", "ChGapN22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {3} refN {-3}", "ChGapP32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiP refP | olP {3} refN {-3}", "ChGapP32", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN {3} refP {-3}", "ChGapN32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN refN | olN {3} refP {-3}", "ChGapN32", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {4} refN {-4}", "ChGapP42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiP refP | olP {4} refN {-4}", "ChGapP42", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN {4} refP {-4}", "ChGapN42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN refN | olN {4} refP {-4}", "ChGapN42", kTRUE));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refFull {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiFull refFull | olFull {2 -2}", "ChFull22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refFull {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiFull refFull | olFull {2 2 -2 -2}", "ChFull24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refFull {2 2 2 -2 -2 -2}", "ChFull26", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refFull {2 2 2 2 -2 -2 -2 -2}", "ChFull28", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refFull {2 2 2 2 2 -2 -2 -2 -2 -2}", "ChFull210", kFALSE));
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
    for (Int_t i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
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
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hEta"), track.eta());

      double pt = track.pt();
      if (cfg.mEfficiency)
        weff = cfg.mEfficiency->GetBinContent(cfg.mEfficiency->FindBin(pt));
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
      registry.fill(HIST("hPt"), pt);
      bool WithinPtPOI = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);       // within POI pT range
      bool WithinPtRef = (cfgCutPtRefMin < pt) && (pt < cfgCutPtRefMax); // within RF pT range
      if (WithinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), wacc * weff, 1);
      if (WithinPtPOI)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), wacc * weff, 2);
      if (WithinPtPOI && WithinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), wacc * weff, 4);
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
