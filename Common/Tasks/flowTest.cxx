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

// flow test for QC of synthetic flow exercise
// cross-PWG effort in tracking studies
// includes basic tracking, V0s and Cascades

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGMM/Mult/DataModel/Index.h" // for Particles2Tracks table

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct flowTest {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> minB{"minB", 0.0f, "min impact parameter"};
  Configurable<float> maxB{"maxB", 20.0f, "max impact parameter"};

  ConfigurableAxis axisB{"axisB", {100, 0.0f, 20.0f}, ""};
  ConfigurableAxis axisPhi{"axisPhi", {100, 0.0f, 2.0f * TMath::Pi()}, ""};
  ConfigurableAxis axisNch{"axisNch", {300, 0.0f, 3000.0f}, "Nch in |eta|<0.8"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f}, "pt axis"};

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
  }

  using recoTracks = soa::Join<aod::TracksIU, aod::TracksExtra>;

  void process(aod::McCollision const& mcCollision, soa::Join<aod::McParticles, aod::ParticlesToTracks> const& mcParticles, recoTracks const&)
  {

    float imp = mcCollision.impactParameter();
    float evPhi = mcCollision.eventPlaneAngle();
    if (evPhi < 0)
      evPhi += 2. * TMath::Pi();

    long nCh = 0;

    if (imp > minB && imp < maxB) {
      // event within range
      histos.fill(HIST("hImpactParameter"), imp);
      histos.fill(HIST("hEventPlaneAngle"), evPhi);

      for (auto const& mcParticle : mcParticles) {
        // focus on bulk: e, mu, pi, k, p
        int pdgCode = TMath::Abs(mcParticle.pdgCode());
        if (pdgCode != 11 && pdgCode != 13 && pdgCode != 211 && pdgCode != 321 && pdgCode != 2212)
          continue;

        if (!mcParticle.isPhysicalPrimary())
          continue;
        if (TMath::Abs(mcParticle.eta()) > 0.8) // main acceptance
          continue;

        float deltaPhi = mcParticle.phi() - mcCollision.eventPlaneAngle();
        if (deltaPhi < 0)
          deltaPhi += 2. * TMath::Pi();
        if (deltaPhi > 2. * TMath::Pi())
          deltaPhi -= 2. * TMath::Pi();
        histos.fill(HIST("hPtVsPhiGenerated"), deltaPhi, mcParticle.pt());
        histos.fill(HIST("hBVsPtVsPhiGenerated"), imp, deltaPhi, mcParticle.pt());

        nCh++;

        bool validGlobal = false;
        bool validTrack = false;
        bool validTPCTrack = false;
        bool validITSTrack = false;
        bool validITSABTrack = false;
        if (mcParticle.has_tracks()) {
          auto const& tracks = mcParticle.tracks_as<recoTracks>();
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

        // if valid global, fill
        if (validGlobal) {
          histos.fill(HIST("hPtVsPhiGlobal"), deltaPhi, mcParticle.pt());
          histos.fill(HIST("hBVsPtVsPhiGlobal"), imp, deltaPhi, mcParticle.pt());
        }
        // if any track present, fill
        if (validTrack)
          histos.fill(HIST("hBVsPtVsPhiAny"), imp, deltaPhi, mcParticle.pt());
        if (validTPCTrack)
          histos.fill(HIST("hBVsPtVsPhiTPCTrack"), imp, deltaPhi, mcParticle.pt());
        if (validITSTrack)
          histos.fill(HIST("hBVsPtVsPhiITSTrack"), imp, deltaPhi, mcParticle.pt());
        if (validITSABTrack)
          histos.fill(HIST("hBVsPtVsPhiITSABTrack"), imp, deltaPhi, mcParticle.pt());
      }
    }
    histos.fill(HIST("hNchVsImpactParameter"), imp, nCh);
  }

  using LabeledCascades = soa::Join<aod::CascDataExt, aod::McCascLabels>;

  void processCascade(aod::McParticle const& mcParticle, soa::SmallGroups<LabeledCascades> const& cascades, recoTracks const&, aod::McCollisions const&)
  {
    auto mcCollision = mcParticle.mcCollision();
    float imp = mcCollision.impactParameter();

    int pdgCode = TMath::Abs(mcParticle.pdgCode());
    if (pdgCode != 3312 && pdgCode != 3334)
      return;

    if (!mcParticle.isPhysicalPrimary())
      return;
    if (TMath::Abs(mcParticle.eta()) > 0.8)
      return;

    float deltaPhi = mcParticle.phi() - mcCollision.eventPlaneAngle();
    if (deltaPhi < 0)
      deltaPhi += 2. * TMath::Pi();
    if (deltaPhi > 2. * TMath::Pi())
      deltaPhi -= 2. * TMath::Pi();
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
  PROCESS_SWITCH(flowTest, processCascade, "Process cascades", true);

  using LabeledV0s = soa::Join<aod::V0Datas, aod::McV0Labels>;

  void processV0s(aod::McParticle const& mcParticle, soa::SmallGroups<LabeledV0s> const& v0s, recoTracks const&, aod::McCollisions const&)
  {
    auto mcCollision = mcParticle.mcCollision();
    float imp = mcCollision.impactParameter();

    int pdgCode = TMath::Abs(mcParticle.pdgCode());
    if (pdgCode != 310 && pdgCode != 3122)
      return;

    if (!mcParticle.isPhysicalPrimary())
      return;
    if (TMath::Abs(mcParticle.eta()) > 0.8)
      return;

    float deltaPhi = mcParticle.phi() - mcCollision.eventPlaneAngle();
    if (deltaPhi < 0)
      deltaPhi += 2. * TMath::Pi();
    if (deltaPhi > 2. * TMath::Pi())
      deltaPhi -= 2. * TMath::Pi();
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
  PROCESS_SWITCH(flowTest, processV0s, "Process V0s", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowTest>(cfgc)};
}
