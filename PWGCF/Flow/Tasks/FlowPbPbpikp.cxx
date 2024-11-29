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
#include <vector>
#include <iostream>
#include <utility>
#include <array>
#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/StepTHn.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Multiplicity.h"
#include "CommonConstants/PhysicsConstants.h"

#include "PWGCF/GenericFramework/Core/GFWPowerArray.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWCumulant.h"
#include "PWGCF/GenericFramework/Core/FlowContainer.h"

#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"

#include <TProfile.h>
#include <TRandom3.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct GfwPidflow {
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

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
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 5.00, 6.00, 8.00, 10.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisparticles{"axisparticles", {3, 0, 3}, "axis for different hadrons"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TAxis* fPtAxis;
  TRandom3* fRndm = new TRandom3(0);

  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::pidBayes, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;

  void init(InitContext const&)
  {
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    histos.add("hPhi", "", {HistType::kTH1D, {axisPhi}});
    histos.add("hEta", "", {HistType::kTH1D, {axisEta}});
    histos.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
    histos.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    histos.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});
    histos.add("hPt", "", {HistType::kTH1D, {axisPt}});
    histos.add("c22_gap08", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c22_gap08_pi", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c22_gap08_ka", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c22_gap08_pr", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c24_full", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("TofTpcNsigma", "", {HistType::kTHnSparseD, {{axisparticles, axisNsigmaTPC, axisNsigmaTOF}}});

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* PtBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, PtBins);

    TObjArray* oba = new TObjArray();
    oba->Add(new TNamed("Ch08Gap22", "Ch08Gap22"));
    for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Ch08Gap22_pt_%i", i + 1), "Ch08Gap22_pTDiff"));
    oba->Add(new TNamed("Pi08Gap22", "Pi08Gap22"));
    for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Pi08Gap22_pt_%i", i + 1), "Pi08Gap22_pTDiff"));
    oba->Add(new TNamed("Ka08Gap22", "Ka08Gap22"));
    for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Ka08Gap22_pt_%i", i + 1), "Ka08Gap22_pTDiff"));
    oba->Add(new TNamed("Pr08Gap22", "Pr08Gap22"));
    for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Pr08Gap22_pt_%i", i + 1), "Pr08Gap22_pTDiff"));
    oba->Add(new TNamed("ChFull24", "ChFull24"));
    for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("ChFull24_pt_%i", i + 1), "ChFull24_pTDiff"));

    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fPtAxis);
    fFC->Initialize(oba, axisMultiplicity, cfgNbootstrap);
    delete oba;

    fGFW->AddRegion("refN08", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP08", 0.4, 0.8, 1, 1);
    fGFW->AddRegion("full", -0.8, 0.8, 1, 512);

    // charged parts
    fGFW->AddRegion("poiN", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 128);
    fGFW->AddRegion("olN", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 256);

    // pion
    fGFW->AddRegion("poiNpi", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("olNpi", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 16);

    // kaon
    fGFW->AddRegion("poiNk", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 4);
    fGFW->AddRegion("olNk", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 32);

    // proton
    fGFW->AddRegion("poiNpr", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 8);
    fGFW->AddRegion("olNpr", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 64);

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Ch08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Pi08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Ka08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Pr08Gap22", kFALSE));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN refN08 | olN {2} refP08 {-2}", "Ch08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNpi refN08 | olNpi {2} refP08 {-2}", "Pi08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNk refN08 | olNk {2} refP08 {-2}", "Ka08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNpr refN08 | olNpr {2} refP08 {-2}", "Pr08Gap22", kTRUE));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    fGFW->CreateRegions();
  }

  template <typename TTrack>
  std::pair<int, int> GetBayesID(TTrack track)
  {
    std::array<int, 3> bayesprobs = {static_cast<int>(track.bayesPi()), static_cast<int>(track.bayesKa()), static_cast<int>(track.bayesPr())};
    int bayesid = -1;
    int prob = 0;

    for (int i = 0; i < 3; ++i) {
      if (bayesprobs[i] > prob && bayesprobs[i] > 80) {
        bayesid = i;
        prob = bayesprobs[i];
      }
    }
    return std::make_pair(bayesid, prob);
  }

  template <typename TTrack>
  int GetBayesPIDIndex(TTrack track)
  {
    int maxProb[3] = {80, 80, 80};
    int pidID = -1;
    std::pair<int, int> idprob = GetBayesID(track);
    if (idprob.first == 0 || idprob.first == 1 || idprob.first == 2) { // 0 = pion, 1 = kaon, 2 = proton
      pidID = idprob.first;
      float nsigmaTPC[3] = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
      if (idprob.second > maxProb[pidID]) {
        if (abs(nsigmaTPC[pidID]) > 3)
          return 0;
        return pidID + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
      } else {
        return 0;
      }
    }
    return 0;
  }

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
        histos.fill(tarName, cent, val, dnx);
      return;
    }
    for (Int_t i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        histos.fill(tarName, fPtAxis->GetBinCenter(i), val, dnx);
    }
    return;
  }

  void FillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0) {
      return;
    }
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1) {
        fFC->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      }
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
    int Ntot = tracks.size();
    if (Ntot < 1)
      return;
    if (!collision.sel8())
      return;
    float l_Random = fRndm->Rndm();

    float vtxz = collision.posZ();
    histos.fill(HIST("hVtxZ"), vtxz);
    histos.fill(HIST("hMult"), Ntot);
    histos.fill(HIST("hCent"), collision.centFT0C());
    fGFW->Clear();
    const auto cent = collision.centFT0C();
    float weff = 1, wacc = 1;
    int pidIndex;
    for (auto& track : tracks) {
      double pt = track.pt();
      histos.fill(HIST("hPhi"), track.phi());
      histos.fill(HIST("hEta"), track.eta());
      histos.fill(HIST("hPt"), pt);

      histos.fill(HIST("TofTpcNsigma"), 0, track.tpcNSigmaPi(), track.tofNSigmaPi());
      histos.fill(HIST("TofTpcNsigma"), 1, track.tpcNSigmaKa(), track.tofNSigmaKa());
      histos.fill(HIST("TofTpcNsigma"), 2, track.tpcNSigmaPr(), track.tofNSigmaPr());

      bool WithinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
      bool WithinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);       // within RF pT range

      pidIndex = GetBayesPIDIndex(track);
      if (WithinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), wacc * weff, 1);
      if (WithinPtPOI)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), wacc * weff, 128);
      if (WithinPtPOI && WithinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), wacc * weff, 256);
      fGFW->Fill(track.eta(), 1, track.phi(), wacc * weff, 512);

      if (pidIndex) {
        if (WithinPtPOI)
          fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), wacc * weff, 1 << (pidIndex));
        if (WithinPtPOI && WithinPtRef)
          fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), wacc * weff, 1 << (pidIndex + 3));
      }
    }

    // Filling c22 with ROOT TProfile
    FillProfile(corrconfigs.at(0), HIST("c22_gap08"), cent);
    FillProfile(corrconfigs.at(1), HIST("c22_gap08_pi"), cent);
    FillProfile(corrconfigs.at(2), HIST("c22_gap08_ka"), cent);
    FillProfile(corrconfigs.at(3), HIST("c22_gap08_pr"), cent);
    FillProfile(corrconfigs.at(4), HIST("c24_full"), cent);

    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      FillFC(corrconfigs.at(l_ind), cent, l_Random);
    }

  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<GfwPidflow>(cfgc)};
}
