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
#include <string>
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/DataModel/TrackSelectionTables.h"
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
  O2_DEFINE_CONFIGURABLE(cfgFlowAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgFlowEfficiency, std::string, "", "CCDB path to efficiency object")

  ConfigurableAxis axisB{"axisB", {100, 0.0f, 20.0f}, ""};
  ConfigurableAxis axisPhi{"axisPhi", {100, 0.0f, constants::math::TwoPI}, ""};
  ConfigurableAxis axisNch{"axisNch", {300, 0.0f, 3000.0f}, "Nch in |eta|<0.8"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f}, "pt axis"};

  // Corrections
  TH1D* mEfficiency = nullptr;
  GFWWeights* mAcceptance = nullptr;
  bool correctionsLoaded = false;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};

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

    if (cfgOutputNUAWeights) {
      o2::framework::AxisSpec axis = axisPt;
      int nPtBins = axis.binEdges.size() - 1;
      double* ptBins = &(axis.binEdges)[0];

      fWeights->setPtBins(nPtBins, ptBins);
      fWeights->init(true, false);
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

  using RecoTracks = soa::Join<aod::TracksIU, aod::TracksExtra>;

  void process(aod::McCollision const& mcCollision, aod::BCsWithTimestamps const&, soa::Join<aod::McParticles, aod::ParticlesToTracks> const& mcParticles, RecoTracks const&)
  {

    float imp = mcCollision.impactParameter();
    float evPhi = mcCollision.eventPlaneAngle();
    float vtxz = mcCollision.posZ();
    if (evPhi < 0)
      evPhi += constants::math::TwoPI;

    int64_t nCh = 0;
    float weff = 1.;
    float wacc = 1.;
    auto bc = mcCollision.bc_as<aod::BCsWithTimestamps>();
    loadCorrections(bc.timestamp());

    if (imp > minB && imp < maxB) {
      // event within range
      histos.fill(HIST("hImpactParameter"), imp);
      histos.fill(HIST("hEventPlaneAngle"), evPhi);

      for (auto const& mcParticle : mcParticles) {
        // focus on bulk: e, mu, pi, k, p
        int pdgCode = std::abs(mcParticle.pdgCode());
        if (pdgCode != 11 && pdgCode != 13 && pdgCode != 211 && pdgCode != 321 && pdgCode != 2212)
          continue;

        if (!mcParticle.isPhysicalPrimary())
          continue;
        if (std::fabs(mcParticle.eta()) > 0.8) // main acceptance
          continue;

        float deltaPhi = mcParticle.phi() - mcCollision.eventPlaneAngle();
        if (deltaPhi < 0)
          deltaPhi += constants::math::TwoPI;
        if (deltaPhi > constants::math::TwoPI)
          deltaPhi -= constants::math::TwoPI;
        histos.fill(HIST("hPtVsPhiGenerated"), deltaPhi, mcParticle.pt());
        histos.fill(HIST("hBVsPtVsPhiGenerated"), imp, deltaPhi, mcParticle.pt());

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
        if (cfgOutputNUAWeights && withinPtRef)
          fWeights->fill(mcParticle.phi(), mcParticle.eta(), vtxz, mcParticle.pt(), 0, 0);
        if (!setCurrentParticleWeights(weff, wacc, mcParticle.phi(), mcParticle.eta(), mcParticle.pt(), vtxz))
          continue;
        if (withinPtRef) {
          histos.fill(HIST("hPhi"), mcParticle.phi());
          histos.fill(HIST("hPhiWeighted"), mcParticle.phi(), wacc);
        }

        // if valid global, fill
        if (validGlobal) {
          histos.fill(HIST("hPtVsPhiGlobal"), deltaPhi, mcParticle.pt(), wacc * weff);
          histos.fill(HIST("hBVsPtVsPhiGlobal"), imp, deltaPhi, mcParticle.pt(), wacc * weff);
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
    }
    histos.fill(HIST("hNchVsImpactParameter"), imp, nCh);
  }

  using LabeledCascades = soa::Join<aod::CascDataExt, aod::McCascLabels>;

  void processCascade(aod::McParticle const& mcParticle, soa::SmallGroups<LabeledCascades> const& cascades, RecoTracks const&, aod::McCollisions const&)
  {
    auto mcCollision = mcParticle.mcCollision();
    float imp = mcCollision.impactParameter();

    int pdgCode = std::abs(mcParticle.pdgCode());
    if (pdgCode != 3312 && pdgCode != 3334)
      return;

    if (!mcParticle.isPhysicalPrimary())
      return;
    if (std::fabs(mcParticle.eta()) > 0.8)
      return;

    float deltaPhi = mcParticle.phi() - mcCollision.eventPlaneAngle();
    if (deltaPhi < 0)
      deltaPhi += constants::math::TwoPI;
    if (deltaPhi > constants::math::TwoPI)
      deltaPhi -= constants::math::TwoPI;
    if (pdgCode == 3312)
      histos.fill(HIST("hBVsPtVsPhiGeneratedXi"), imp, deltaPhi, mcParticle.pt());
    if (pdgCode == 3334)
      histos.fill(HIST("hBVsPtVsPhiGeneratedOmega"), imp, deltaPhi, mcParticle.pt());

    if (cascades.size() > 0) {
      if (pdgCode == 3312)
        histos.fill(HIST("hBVsPtVsPhiGlobalXi"), imp, deltaPhi, mcParticle.pt());
      if (pdgCode == 3334)
        histos.fill(HIST("hBVsPtVsPhiGlobalOmega"), imp, deltaPhi, mcParticle.pt());
    }
  }
  PROCESS_SWITCH(FlowMc, processCascade, "Process cascades", true);

  using LabeledV0s = soa::Join<aod::V0Datas, aod::McV0Labels>;

  void processV0s(aod::McParticle const& mcParticle, soa::SmallGroups<LabeledV0s> const& v0s, RecoTracks const&, aod::McCollisions const&)
  {
    auto mcCollision = mcParticle.mcCollision();
    float imp = mcCollision.impactParameter();

    int pdgCode = std::abs(mcParticle.pdgCode());
    if (pdgCode != 310 && pdgCode != 3122)
      return;

    if (!mcParticle.isPhysicalPrimary())
      return;
    if (std::fabs(mcParticle.eta()) > 0.8)
      return;

    float deltaPhi = mcParticle.phi() - mcCollision.eventPlaneAngle();
    if (deltaPhi < 0)
      deltaPhi += constants::math::TwoPI;
    if (deltaPhi > constants::math::TwoPI)
      deltaPhi -= constants::math::TwoPI;
    if (pdgCode == 310)
      histos.fill(HIST("hBVsPtVsPhiGeneratedK0Short"), imp, deltaPhi, mcParticle.pt());
    if (pdgCode == 3122)
      histos.fill(HIST("hBVsPtVsPhiGeneratedLambda"), imp, deltaPhi, mcParticle.pt());

    if (v0s.size() > 0) {
      if (pdgCode == 310)
        histos.fill(HIST("hBVsPtVsPhiGlobalK0Short"), imp, deltaPhi, mcParticle.pt());
      if (pdgCode == 3122)
        histos.fill(HIST("hBVsPtVsPhiGlobalLambda"), imp, deltaPhi, mcParticle.pt());
    }
  }
  PROCESS_SWITCH(FlowMc, processV0s, "Process V0s", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowMc>(cfgc)};
}
