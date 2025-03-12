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

/// \file   flowMc.cxx
/// \author Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  Feb/5/2025
/// \brief  QC of synthetic flow exercise

#include <CCDB/BasicCCDBManager.h>
#include <vector>
#include <string>
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGMM/Mult/DataModel/Index.h" // for Particles2Tracks table
#include "GFWPowerArray.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "GFWWeights.h"
#include "FlowContainer.h"
#include <TProfile.h>
#include <TRandom3.h>
#include <TPDGCode.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowMc {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> minB{"minB", 0.0f, "min impact parameter"};
  Configurable<float> maxB{"maxB", 20.0f, "max impact parameter"};
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, false, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgFlowAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgFlowEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgCentVsIPTruth, std::string, "", "CCDB path to centrality vs IP truth")
  O2_DEFINE_CONFIGURABLE(cfgFlowCumulantEnabled, bool, false, "switch of calculating flow")
  O2_DEFINE_CONFIGURABLE(cfgFlowCumulantNbootstrap, int, 30, "Number of subsamples")

  ConfigurableAxis axisB{"axisB", {100, 0.0f, 20.0f}, ""};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "X axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {100, 0.0f, constants::math::TwoPI}, ""};
  ConfigurableAxis axisNch{"axisNch", {300, 0.0f, 3000.0f}, "Nch in |eta|<0.8"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f}, "pt axis"};

  // Cent vs IP
  TH1D* mCentVsIPTruth = nullptr;
  bool centVsIPTruthLoaded = false;

  // Corrections
  TH1D* mEfficiency = nullptr;
  GFWWeights* mAcceptance = nullptr;
  bool correctionsLoaded = false;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  OutputObj<FlowContainer> fFCTrue{FlowContainer("FlowContainerTrue")};
  OutputObj<FlowContainer> fFCReco{FlowContainer("FlowContainerReco")};
  GFW* fGFWTrue = new GFW();
  GFW* fGFWReco = new GFW();
  TAxis* fPtAxis;
  std::vector<GFW::CorrConfig> corrconfigsTruth;
  std::vector<GFW::CorrConfig> corrconfigsReco;
  TRandom3* fRndm = new TRandom3(0);

  void init(InitContext&)
  {
    // pT histograms
    histos.add<TH1>("hImpactParameter", "hImpactParameter", HistType::kTH1D, {axisB});
    histos.add<TH2>("hNchVsImpactParameter", "hNchVsImpactParameter", HistType::kTH2D, {axisB, axisNch});
    histos.add<TH1>("hEventPlaneAngle", "hEventPlaneAngle", HistType::kTH1D, {axisPhi});
    histos.add<TH2>("hPtVsPhiGenerated", "hPtVsPhiGenerated", HistType::kTH2D, {axisPhi, axisPt});
    histos.add<TH2>("hPtVsPhiGlobal", "hPtVsPhiGlobal", HistType::kTH2D, {axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGenerated", "hBVsPtVsPhiGenerated", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGlobal", "hBVsPtVsPhiGlobal", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiAny", "hBVsPtVsPhiAny", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiTPCTrack", "hBVsPtVsPhiTPCTrack", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiITSTrack", "hBVsPtVsPhiITSTrack", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiITSABTrack", "hBVsPtVsPhiITSABTrack", HistType::kTH3D, {axisB, axisPhi, axisPt});

    histos.add<TH3>("hBVsPtVsPhiGeneratedK0Short", "hBVsPtVsPhiGeneratedK0Short", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGlobalK0Short", "hBVsPtVsPhiGlobalK0Short", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGeneratedLambda", "hBVsPtVsPhiGeneratedLambda", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGlobalLambda", "hBVsPtVsPhiGlobalLambda", HistType::kTH3D, {axisB, axisPhi, axisPt});

    histos.add<TH3>("hBVsPtVsPhiGeneratedXi", "hBVsPtVsPhiGeneratedXi", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGlobalXi", "hBVsPtVsPhiGlobalXi", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGeneratedOmega", "hBVsPtVsPhiGeneratedOmega", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGlobalOmega", "hBVsPtVsPhiGlobalOmega", HistType::kTH3D, {axisB, axisPhi, axisPt});

    histos.add<TH1>("hPhi", "#phi distribution", HistType::kTH1D, {axisPhi});
    histos.add<TH1>("hPhiWeighted", "corrected #phi distribution", HistType::kTH1D, {axisPhi});
    histos.add<TH2>("hEPVsPhiMC", "hEPVsPhiMC;Event Plane Angle; #varphi", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add<TH2>("hEPVsPhi", "hEPVsPhi;Event Plane Angle; #varphi", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add<TH2>("hPtNchGenerated", "Reco production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH2>("hPtNchGlobal", "Global production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH1>("hPtMCGen", "Monte Carlo Truth; pT (GeV/c);", {HistType::kTH1D, {axisPt}});
    histos.add<TH1>("hPtMCGlobal", "Monte Carlo Global; pT (GeV/c);", {HistType::kTH1D, {axisPt}});

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    if (cfgOutputNUAWeights) {
      fWeights->setPtBins(nPtBins, ptBins);
      fWeights->init(true, false);
    }

    if (cfgFlowCumulantEnabled) {
      TObjArray* oba = new TObjArray();
      oba->Add(new TNamed("ChFull22", "ChFull22"));
      for (auto i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("ChFull22_pt_%i", i + 1), "ChFull22_pTDiff"));
      oba->Add(new TNamed("Ch10Gap22", "Ch10Gap22"));
      for (auto i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("Ch10Gap22_pt_%i", i + 1), "Ch10Gap22_pTDiff"));
      fFCTrue->SetName("FlowContainerTrue");
      fFCTrue->SetXAxis(fPtAxis);
      fFCTrue->Initialize(oba, axisCentrality, cfgFlowCumulantNbootstrap);
      fFCReco->SetName("FlowContainerReco");
      fFCReco->SetXAxis(fPtAxis);
      fFCReco->Initialize(oba, axisCentrality, cfgFlowCumulantNbootstrap);
      delete oba;

      fGFWTrue->AddRegion("full", -0.8, 0.8, 1, 1);
      fGFWTrue->AddRegion("refN10", -0.8, -0.5, 1, 1);
      fGFWTrue->AddRegion("refP10", 0.5, 0.8, 1, 1);
      fGFWTrue->AddRegion("poiN10", -0.8, -0.5, 1 + fPtAxis->GetNbins(), 2);
      fGFWTrue->AddRegion("poifull", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 2);
      fGFWTrue->AddRegion("olN10", -0.8, -0.5, 1 + fPtAxis->GetNbins(), 4);
      fGFWTrue->AddRegion("olfull", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 4);
      corrconfigsTruth.push_back(fGFWTrue->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
      corrconfigsTruth.push_back(fGFWTrue->GetCorrelatorConfig("poifull full | olfull {2 -2}", "ChFull22", kTRUE));
      corrconfigsTruth.push_back(fGFWTrue->GetCorrelatorConfig("refN10 {2} refP10 {-2}", "Ch10Gap22", kFALSE));
      corrconfigsTruth.push_back(fGFWTrue->GetCorrelatorConfig("poiN10 refN10 | olN10 {2} refP10 {-2}", "Ch10Gap22", kTRUE));
      fGFWTrue->CreateRegions();

      fGFWReco->AddRegion("full", -0.8, 0.8, 1, 1);
      fGFWReco->AddRegion("refN10", -0.8, -0.5, 1, 1);
      fGFWReco->AddRegion("refP10", 0.5, 0.8, 1, 1);
      fGFWReco->AddRegion("poiN10", -0.8, -0.5, 1 + fPtAxis->GetNbins(), 2);
      fGFWReco->AddRegion("poifull", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 2);
      fGFWReco->AddRegion("olN10", -0.8, -0.5, 1 + fPtAxis->GetNbins(), 4);
      fGFWReco->AddRegion("olfull", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 4);
      corrconfigsReco.push_back(fGFWReco->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
      corrconfigsReco.push_back(fGFWReco->GetCorrelatorConfig("poifull full | olfull {2 -2}", "ChFull22", kTRUE));
      corrconfigsReco.push_back(fGFWReco->GetCorrelatorConfig("refN10 {2} refP10 {-2}", "Ch10Gap22", kFALSE));
      corrconfigsReco.push_back(fGFWReco->GetCorrelatorConfig("poiN10 refN10 | olN10 {2} refP10 {-2}", "Ch10Gap22", kTRUE));
      fGFWReco->CreateRegions();
    }
  }

  void loadCorrections(uint64_t timestamp)
  {
    if (correctionsLoaded)
      return;
    if (cfgFlowAcceptance.value.empty() == false) {
      mAcceptance = ccdb->getForTimeStamp<GFWWeights>(cfgFlowAcceptance, timestamp);
      if (mAcceptance)
        LOGF(info, "Loaded acceptance weights from %s (%p)", cfgFlowAcceptance.value.c_str(), (void*)mAcceptance);
      else
        LOGF(warning, "Could not load acceptance weights from %s (%p)", cfgFlowAcceptance.value.c_str(), (void*)mAcceptance);
    }
    if (cfgFlowEfficiency.value.empty() == false) {
      mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgFlowEfficiency, timestamp);
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgFlowEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgFlowEfficiency.value.c_str(), (void*)mEfficiency);
    }
    correctionsLoaded = true;
  }

  bool setCurrentParticleWeights(float& weight_nue, float& weight_nua, float phi, float eta, float pt, float vtxz)
  {
    float eff = 1.;
    if (mEfficiency)
      eff = mEfficiency->GetBinContent(mEfficiency->FindBin(pt));
    else
      eff = 1.0;
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    if (mAcceptance)
      weight_nua = mAcceptance->getNUA(phi, eta, vtxz);
    else
      weight_nua = 1;
    return true;
  }

  void fillFC(GFW* fGFW, bool isMCTruth, const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (!corrconf.pTDif) {
      if (dnx == 0)
        return;
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        if (isMCTruth)
          fFCTrue->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
        else
          fFCReco->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      }
      return;
    }
    for (auto i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        if (isMCTruth)
          fFCTrue->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
        else
          fFCReco->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
      }
    }
    return;
  }

  void loadCentVsIPTruth(uint64_t timestamp)
  {
    if (centVsIPTruthLoaded)
      return;
    if (cfgCentVsIPTruth.value.empty() == false) {
      mCentVsIPTruth = ccdb->getForTimeStamp<TH1D>(cfgCentVsIPTruth, timestamp);
      if (mCentVsIPTruth)
        LOGF(info, "Loaded CentVsIPTruth weights from %s (%p)", cfgCentVsIPTruth.value.c_str(), (void*)mCentVsIPTruth);
      else
        LOGF(fatal, "Failed to load CentVsIPTruth weights from %s", cfgCentVsIPTruth.value.c_str());

      centVsIPTruthLoaded = true;
    } else {
      LOGF(fatal, "when calculate flow, Cent Vs IP distribution must be provided");
    }
  }

  using RecoTracks = soa::Join<aod::TracksIU, aod::TracksExtra>;

  void process(aod::McCollision const& mcCollision, aod::BCsWithTimestamps const&, soa::Join<aod::McParticles, aod::ParticlesToTracks> const& mcParticles, RecoTracks const&)
  {

    float imp = mcCollision.impactParameter();
    float evPhi = mcCollision.eventPlaneAngle();
    float vtxz = mcCollision.posZ();
    evPhi = RecoDecay::constrainAngle(evPhi);

    int64_t nCh = 0;
    int64_t nChGlobal = 0;
    float centrality = 0;
    float lRandom = fRndm->Rndm();
    float weff = 1.;
    float wacc = 1.;
    auto bc = mcCollision.bc_as<aod::BCsWithTimestamps>();
    loadCorrections(bc.timestamp());

    if (imp > minB && imp < maxB) {
      // event within range
      histos.fill(HIST("hImpactParameter"), imp);
      histos.fill(HIST("hEventPlaneAngle"), evPhi);
      if (cfgFlowCumulantEnabled) {
        loadCentVsIPTruth(bc.timestamp());
        centrality = mCentVsIPTruth->GetBinContent(mCentVsIPTruth->GetXaxis()->FindBin(imp));
        fGFWTrue->Clear();
        fGFWReco->Clear();
      }

      for (auto const& mcParticle : mcParticles) {
        int pdgCode = std::abs(mcParticle.pdgCode());
        if (pdgCode != PDG_t::kElectron && pdgCode != PDG_t::kMuonMinus && pdgCode != PDG_t::kPiPlus && pdgCode != kKPlus && pdgCode != PDG_t::kProton)
          continue;
        if (!mcParticle.isPhysicalPrimary())
          continue;
        if (std::fabs(mcParticle.eta()) > 0.8) // main acceptance
          continue;
        if (mcParticle.has_tracks()) {
          auto const& tracks = mcParticle.tracks_as<RecoTracks>();
          for (auto const& track : tracks) {
            if (track.hasTPC() && track.hasITS())
              nChGlobal++;
          }
        }
      }

      for (auto const& mcParticle : mcParticles) {
        // focus on bulk: e, mu, pi, k, p
        int pdgCode = std::abs(mcParticle.pdgCode());
        if (pdgCode != PDG_t::kElectron && pdgCode != PDG_t::kMuonMinus && pdgCode != PDG_t::kPiPlus && pdgCode != kKPlus && pdgCode != PDG_t::kProton)
          continue;

        if (!mcParticle.isPhysicalPrimary())
          continue;
        if (std::fabs(mcParticle.eta()) > 0.8) // main acceptance
          continue;

        float deltaPhi = mcParticle.phi() - mcCollision.eventPlaneAngle();
        deltaPhi = RecoDecay::constrainAngle(deltaPhi);
        histos.fill(HIST("hPtVsPhiGenerated"), deltaPhi, mcParticle.pt());
        histos.fill(HIST("hBVsPtVsPhiGenerated"), imp, deltaPhi, mcParticle.pt());
        histos.fill(HIST("hPtNchGenerated"), mcParticle.pt(), nChGlobal);
        histos.fill(HIST("hPtMCGen"), mcParticle.pt());

        nCh++;

        bool validGlobal = false;
        bool validTrack = false;
        bool validTPCTrack = false;
        bool validITSTrack = false;
        bool validITSABTrack = false;
        if (mcParticle.has_tracks()) {
          auto const& tracks = mcParticle.tracks_as<RecoTracks>();
          for (auto const& track : tracks) {
            if (track.hasTPC() && track.hasITS()) {
              validGlobal = true;
            }
            if (track.hasTPC() || track.hasITS()) {
              validTrack = true;
            }
            if (track.hasTPC()) {
              validTPCTrack = true;
            }
            if (track.hasITS() && track.itsChi2NCl() > -1e-6) {
              validITSTrack = true;
            }
            if (track.hasITS() && track.itsChi2NCl() < -1e-6) {
              validITSABTrack = true;
            }
          }
        }

        bool withinPtRef = (cfgCutPtRefMin < mcParticle.pt()) && (mcParticle.pt() < cfgCutPtRefMax); // within RF pT range
        bool withinPtPOI = (cfgCutPtPOIMin < mcParticle.pt()) && (mcParticle.pt() < cfgCutPtPOIMax); // within POI pT range
        if (cfgOutputNUAWeights && withinPtRef)
          fWeights->fill(mcParticle.phi(), mcParticle.eta(), vtxz, mcParticle.pt(), 0, 0);
        if (!setCurrentParticleWeights(weff, wacc, mcParticle.phi(), mcParticle.eta(), mcParticle.pt(), vtxz))
          continue;

        if (cfgFlowCumulantEnabled) {
          if (withinPtRef)
            fGFWTrue->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 1);
          if (withinPtPOI)
            fGFWTrue->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 2);
          if (withinPtPOI && withinPtRef)
            fGFWTrue->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 4);

          if (validGlobal) {
            if (withinPtRef)
              fGFWReco->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 1);
            if (withinPtPOI)
              fGFWReco->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 2);
            if (withinPtPOI && withinPtRef)
              fGFWReco->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 4);
          }
        }

        if (withinPtRef) {
          histos.fill(HIST("hEPVsPhiMC"), evPhi, mcParticle.phi());
        }

        if (validGlobal && withinPtRef) {
          histos.fill(HIST("hPhi"), mcParticle.phi());
          histos.fill(HIST("hPhiWeighted"), mcParticle.phi(), wacc);
          histos.fill(HIST("hEPVsPhi"), evPhi, mcParticle.phi());
        }

        // if valid global, fill
        if (validGlobal) {
          histos.fill(HIST("hPtVsPhiGlobal"), deltaPhi, mcParticle.pt(), wacc * weff);
          histos.fill(HIST("hBVsPtVsPhiGlobal"), imp, deltaPhi, mcParticle.pt(), wacc * weff);
          histos.fill(HIST("hPtNchGlobal"), mcParticle.pt(), nChGlobal);
          histos.fill(HIST("hPtMCGlobal"), mcParticle.pt());
        }
        // if any track present, fill
        if (validTrack)
          histos.fill(HIST("hBVsPtVsPhiAny"), imp, deltaPhi, mcParticle.pt(), wacc * weff);
        if (validTPCTrack)
          histos.fill(HIST("hBVsPtVsPhiTPCTrack"), imp, deltaPhi, mcParticle.pt(), wacc * weff);
        if (validITSTrack)
          histos.fill(HIST("hBVsPtVsPhiITSTrack"), imp, deltaPhi, mcParticle.pt(), wacc * weff);
        if (validITSABTrack)
          histos.fill(HIST("hBVsPtVsPhiITSABTrack"), imp, deltaPhi, mcParticle.pt(), wacc * weff);
      }

      if (cfgFlowCumulantEnabled) {
        for (uint l_ind = 0; l_ind < corrconfigsTruth.size(); l_ind++) {
          fillFC(fGFWTrue, true, corrconfigsTruth.at(l_ind), centrality, lRandom);
        }
        for (uint l_ind = 0; l_ind < corrconfigsReco.size(); l_ind++) {
          fillFC(fGFWReco, false, corrconfigsReco.at(l_ind), centrality, lRandom);
        }
      }
    }
    histos.fill(HIST("hNchVsImpactParameter"), imp, nCh);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowMc>(cfgc)};
}
