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
#include <cmath>
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
#include "TList.h"
#include <TProfile.h>
#include <TRandom3.h>
#include <TF1.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowGFWPbPb {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW(); // GFW class used from main src
  std::vector<GFW::CorrConfig> corrconfigs;
  TRandom3* fRndm = new TRandom3(0);
  TAxis* fPtAxis;
  std::vector<std::vector<std::shared_ptr<TProfile>>> BootstrapArray; // TProfile is a shared pointer

  enum ExtraProfile {

    // here are TProfiles for vn-pt correlations that are not implemented in GFW
    kc22,
    kc24,
    kc26,
    kc28,
    kc22etagap,

    // Count the total number of enum
    kCount_ExtraProfile
  };

  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;   // collisions filter
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>>; // tracks filter

  void init(InitContext const&) // Initialization
  {
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    // Add some output objects to the histogram registry
    registry.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{1, 0, 1}}});
    registry.add("hPhi", "", {HistType::kTH1D, {axisPhi}});
    registry.add("hEta", "", {HistType::kTH1D, {axisEta}});
    registry.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
    registry.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    registry.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});
    registry.add("c22", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c24", ";Centrality  (%) ; C_{2}{4}", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c26", ";Centrality  (%) ; C_{2}{6}", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c28", ";Centrality  (%) ; C_{2}{8}", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c22etagap", ";Centrality  (%) ; C_{2}{2} (|#eta| < 0.8) ", {HistType::kTProfile, {axisMultiplicity}});

    // initial array
    BootstrapArray.resize(cfgNbootstrap);
    for (int i = 0; i < cfgNbootstrap; i++) {
      BootstrapArray[i].resize(kCount_ExtraProfile);
    }
    for (int i = 0; i < cfgNbootstrap; i++) {
      BootstrapArray[i][kc22] = registry.add<TProfile>(Form("BootstrapContainer_%d/c22", i), ";Centrality  (%) ; C_{2}{2}", {HistType::kTProfile, {axisMultiplicity}});
      BootstrapArray[i][kc24] = registry.add<TProfile>(Form("BootstrapContainer_%d/c24", i), ";Centrality  (%) ; C_{2}{4}", {HistType::kTProfile, {axisMultiplicity}});
      BootstrapArray[i][kc26] = registry.add<TProfile>(Form("BootstrapContainer_%d/c26", i), ";Centrality  (%) ; C_{2}{6}", {HistType::kTProfile, {axisMultiplicity}});
      BootstrapArray[i][kc28] = registry.add<TProfile>(Form("BootstrapContainer_%d/c28", i), ";Centrality  (%) ; C_{2}{8}", {HistType::kTProfile, {axisMultiplicity}});
      BootstrapArray[i][kc22etagap] = registry.add<TProfile>(Form("BootstrapContainer_%d/c22etagap", i), ";Centrality  (%) ; C_{2}{2} (|#eta| < 0.8)", {HistType::kTProfile, {axisMultiplicity}});
    }

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* PtBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, PtBins);

    // add in FlowContainer to Get boostrap sample automatically
    TObjArray* oba = new TObjArray();
    fFC->SetXAxis(fPtAxis);
    fFC->SetName("FlowContainer");
    fFC->Initialize(oba, axisMultiplicity, cfgNbootstrap);
    delete oba;

    fGFW->AddRegion("full", -0.8, 0.8, 1, 1); // eta region -0.8 to 0.8
    fGFW->AddRegion("refN10", -0.8, -0.5, 1, 1);
    fGFW->AddRegion("refP10", 0.5, 0.8, 1, 1);

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 2 -2 -2 -2}", "ChFull26", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 2 2  -2 -2 -2 -2}", "ChFull28", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {2} refP10 {-2}", "Ch10Gap22", kFALSE));

    fGFW->CreateRegions(); // finalize the initialization

  } // end of Initialization

  template <char... chars>
  void FillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        registry.fill(tarName, cent, val, dnx);
      return;
    }
    return;
  }

  void FillProfile(const GFW::CorrConfig& corrconf, std::shared_ptr<TProfile> tarName, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1) {
        tarName->Fill(cent, val, dnx);
      }
      return;
    }
    return;
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

  void process(aodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, aodTracks const& tracks)
  {
    registry.fill(HIST("hEventCount"), 0.5);
    if (!collision.sel8())
      return;

    int Ntot = tracks.size();
    if (Ntot < 1)
      return;

    // if (!collision.sel7())  return;

    float vtxz = collision.posZ();
    float l_Random = fRndm->Rndm();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), Ntot);
    registry.fill(HIST("hCent"), collision.centFT0C());
    fGFW->Clear();

    const auto cent = collision.centFT0C();
    float weff = 1, wacc = 1;

    for (auto& track : tracks) {
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hEta"), track.eta());

      fGFW->Fill(track.eta(), 1, track.phi(), wacc * weff, 1);
    }

    // Filling c22 with ROOT TProfile
    FillProfile(corrconfigs.at(0), HIST("c22"), cent);
    FillProfile(corrconfigs.at(1), HIST("c24"), cent);
    FillProfile(corrconfigs.at(2), HIST("c26"), cent);
    FillProfile(corrconfigs.at(3), HIST("c28"), cent);
    FillProfile(corrconfigs.at(4), HIST("c22etagap"), cent);

    // Filling Bootstrap Samples
    int SampleIndex = static_cast<int>(cfgNbootstrap * l_Random);
    FillProfile(corrconfigs.at(0), BootstrapArray[SampleIndex][kc22], cent);
    FillProfile(corrconfigs.at(1), BootstrapArray[SampleIndex][kc24], cent);
    FillProfile(corrconfigs.at(2), BootstrapArray[SampleIndex][kc26], cent);
    FillProfile(corrconfigs.at(3), BootstrapArray[SampleIndex][kc28], cent);
    FillProfile(corrconfigs.at(4), BootstrapArray[SampleIndex][kc22etagap], cent);

    // Filling Flow Container
    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      FillFC(corrconfigs.at(l_ind), cent, l_Random);
    }

  } // End of process
};  // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGFWPbPb>(cfgc)};
}
