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

///
/// \file   qaEventTrackLite.cxx
/// \author Mario Krüger <mario.kruger@cern.ch>
/// \author Mattia Faggin <mattia.faggin@cern.ch>
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>
/// \brief  Light version of the task to produce QA objects for the track and the event properties.
///         This task runs on prefiltered data
///

#include "qaEventTrack.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

// Lite version of the QA task to run on skimmed dataset
struct qaEventTrackLite {
  // Binning
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0}, ""};
  ConfigurableAxis binsImpPar{"binsImpPar", {200, -0.15, 0.15}, "Impact parameter binning (cm)"};
  ConfigurableAxis binsEta{"binsEta", {800, -2., 2.}, "Eta binning"};
  ConfigurableAxis binsPhi{"binsPhi", {180, 0., 2 * M_PI}, "Phi binning"};

  HistogramRegistry histos;

  Configurable<bool> bItsStandalone{"bItsStandalone", false, "Select only ITS standalone DPG tracks"};
  Configurable<bool> bTpcOnly{"bTpcOnly", false, "Select only TPC only DPG tracks"};
  Configurable<bool> bItsTpcMatched{"bItsTpcMatched", false, "Select ITS-TPC matched DPG tracks"};
  // Kinematic selections
  Configurable<float> ptMin{"ptMin", 0., "Minimum track pt"};
  Configurable<float> etaMin{"etaMin", -10., "Minimum eta for DPG tracks"};
  Configurable<float> etaMax{"etaMax", 10., "Maximum eta for DPG tracks"};
  // ITS selections
  Configurable<float> chi2ItsMin{"chi2ItsMin", -1001.f, "Max ITS chi2"};
  Configurable<float> chi2ItsMax{"chi2ItsMax", 1000.f, "Max ITS chi2"};
  // TPC selections
  Configurable<int> nClusterTpcMin{"nClusterTpcMin", -1001, "Minimum number of TPC clusters"};
  Configurable<int> nCrossedRowsTpcMin{"nCrossedRowsTpcMin", -1001, "Minimum number of TPC crossed rows"};
  Configurable<float> nCrossedRowsTpcOverFindableClustersTpcMin{"nCrossedRowsTpcOverFindableClustersTpcMin", -1, "Minimum ratio between TPC crossed rows and findable clusters"};
  Configurable<float> chi2TpcMin{"chi2TpcMin", -1001.f, "Max TPC chi2"};
  Configurable<float> chi2TpcMax{"chi2TpcMax", 1000.f, "Max TPC chi2"};
  // TOF selections
  Configurable<float> chi2TofMin{"chi2TofMin", -1001.f, "Max TOF chi2"};
  Configurable<float> lengthMin{"lengthMin", -1001.f, "Min length"};

  void init(InitContext const&)
  {
    const AxisSpec axisPt{binsPt, "#it{p}_{T} [GeV/c]"};
    const AxisSpec axisEta{binsEta, "#it{#eta}"};
    const AxisSpec axisPhi{binsPhi, "#it{#phi} [rad]"};

    // kine histograms
    histos.add("Tracks/VertexPositionZ", "", kTH1D, {{100, -20.f, 20.f, "Vertex Z [cm]"}});
    histos.add("Tracks/Kine/pt", "#it{p}_{T}", kTH1D, {axisPt});
    histos.add("Tracks/Kine/eta", "#eta", kTH1D, {axisEta});
    histos.add("Tracks/Kine/phi", "#phi", kTH1D, {axisPhi});
    histos.add("Tracks/length", "track length in cm;#it{Length} [cm];", kTH1D, {{400, 0, 1000}});
    const AxisSpec axisImpParRPhi{binsImpPar, "#it{d}_{r#it{#varphi}} (#cm)"};
    const AxisSpec axisImpParZAxis{binsImpPar, "#it{d}_{z} (#cm)"};
    histos.add("Tracks/dcaXY", "distance of closest approach in #it{xy} plane", kTH1D, {axisImpParRPhi});
    histos.add("Tracks/dcaZ", "distance of closest approach in #it{z}", kTH1D, {axisImpParZAxis});
    histos.add("Tracks/dcaXYvsPt", "d_#it{xy} vs. #it{p}_{T}", kTH2D, {axisPt, axisImpParRPhi});
    histos.add("Tracks/dcaZvsPt", "d_#it{z} vs. #it{p}_{T}", kTH2D, {axisPt, axisImpParRPhi});

    // its histograms
    histos.add("Tracks/ITS/itsChi2NCl", "chi2 per ITS cluster;chi2 / cluster ITS", kTH1D, {{100, 0, 40}});
    histos.add("Tracks/ITS/itsNCl", "ITS number of clusters;# clusters ITS", kTH1D, {{8, -0.5, 7.5}});
    histos.add("Tracks/ITS/itsNClvsItsHitmap", "ITS number of clusters vs. ITS hitmap;# clusters ITS; ITS hitmap", kTH2D, {{8, -0.5, 7.5, "# clusters ITS"}, {128, 0, 128, "ITS hitmap"}});
    // tpc histograms
    histos.add("Tracks/TPC/tpcChi2NCl", "chi2 per cluster in TPC;chi2 / cluster TPC", kTH1D, {{100, 0, 10}});
    histos.add("Tracks/TPC/tpcNClsFound", "number of found TPC clusters;# clusters TPC", kTH1D, {{165, -0.5, 164.5}});
    histos.add("Tracks/TPC/tpcCrossedRows", "number of crossed TPC rows;# crossed rows TPC", kTH1D, {{165, -0.5, 164.5}});
    histos.add("Tracks/TPC/tpcCrossedRowsOverFindableCls", "crossed TPC rows over findable clusters;crossed rows / findable clusters TPC", kTH1D, {{60, 0.7, 1.3}});
    histos.add("Tracks/TPC/tpcNClsFoundvsPt", "", kTH2D, {axisPt, {165, -0.5, 164.5, "# clusters TPC"}});
    histos.add("Tracks/TPC/tpcCrossedRowsvsPt", "", kTH2D, {axisPt, {165, -0.5, 164.5, "# crossed rows TPC"}});
    histos.add("Tracks/TPC/tpcCrossedRowsOverFindableClsvsPt", "", kTH2D, {axisPt, {60, 0.7, 1.3, "crossed rows / findable clusters TPC"}});
    // trd histograms
    histos.add("Tracks/TRD/trdChi2", "chi2 in TRD", kTH1D, {{100, 0, 10, "chi2 / cluster TRD"}});
    // tof histograms
    histos.add("Tracks/TOF/tofChi2", "chi2 in TOF", kTH1D, {{100, 0, 10, "chi2 / cluster TOF"}});
    // matching histogram
    histos.add("Tracks/matchedDet", "matched detectors", kTH1D, {{4, 0.5, 4.5, ""}});
    histos.get<TH1>(HIST("Tracks/matchedDet"))->GetXaxis()->SetBinLabel(1, "hasTPC");
    histos.get<TH1>(HIST("Tracks/matchedDet"))->GetXaxis()->SetBinLabel(2, "hasITS");
    histos.get<TH1>(HIST("Tracks/matchedDet"))->GetXaxis()->SetBinLabel(3, "hasTRD");
    histos.get<TH1>(HIST("Tracks/matchedDet"))->GetXaxis()->SetBinLabel(4, "hasTOF");
    // kinematics
    histos.add("Tracks/relativeResoPt", "relative #it{p}_{T} resolution;#sigma(#it{p}_{T})/#it{p}_{T};#it{p}_{T};", kTH2D, {{axisPt, {500, 0., 1, "#sigma{#it{p}}/#it{p}_{T}"}}});
    histos.add("Tracks/relativeResoPtMean", "mean relative #it{p}_{T} resolution;#LT(#it{p}_{T})/#it{p}_{T}#GT;#it{p}_{T};", kTProfile, {axisPt});

    // MC histograms
    if (doprocessMCLite) {
      histos.add("Particles/PDGs", "Particle PDGs;PDG Code", kTH1D, {{100, 0.f, 100.f}});
      histos.add("Particles/Kine/pt", "Particle #it{p}_{T}", kTH1D, {axisPt});
      histos.add("Particles/Kine/eta", "Particle #eta", kTH1D, {axisEta});
      histos.add("Particles/Kine/phi", "Particle #phi", kTH1D, {axisPhi});
    }
  }

  ///////////////
  /// Filters ///
  ///////////////
  // Kinematics
  Filter ptCut = o2::aod::dpgtrack::pt > ptMin;
  Filter etaCut = etaMin < o2::aod::dpgtrack::eta && o2::aod::dpgtrack::eta < etaMax;
  // Detector matching
  Filter itsStandaloneTracks = (bItsStandalone.node() == false) || (o2::aod::dpgtrack::hasITS == true && o2::aod::dpgtrack::hasTPC == false);
  Filter tpcOnlyTracks = (bTpcOnly.node() == false) || (o2::aod::dpgtrack::hasITS == false && o2::aod::dpgtrack::hasTPC == true);
  Filter itsTpcMatchedTracks = (bItsTpcMatched.node() == false) || (o2::aod::dpgtrack::hasITS == true && o2::aod::dpgtrack::hasTPC == true);
  // ITS
  Filter itsChi2 = (bTpcOnly.node() == true) || (chi2ItsMin < o2::aod::track::itsChi2NCl && o2::aod::track::itsChi2NCl < chi2ItsMax);
  // TPC
  Filter tpcChi2s = (bItsStandalone.node() == true) || (chi2TpcMin < o2::aod::track::tpcChi2NCl && o2::aod::track::tpcChi2NCl < chi2TpcMax);
  Filter tpcNclusters = (bItsStandalone.node() == true) || (o2::aod::dpgtrack::tpcNClsFound > (int16_t)nClusterTpcMin);
  Filter tpcNcrossedRows = (bItsStandalone.node() == true) || (o2::aod::dpgtrack::tpcNClsCrossedRows > (int16_t)nCrossedRowsTpcMin);
  Filter tpcNcrossedRowsOverFindableClusters = (bItsStandalone.node() == true) || (o2::aod::dpgtrack::tpcCrossedRowsOverFindableCls > nCrossedRowsTpcOverFindableClustersTpcMin);
  // TOF
  Filter tofChi = o2::aod::track::tofChi2 > chi2TofMin;
  Filter length = o2::aod::track::length > lengthMin;

  // Process data
  void processDataLite(o2::soa::Filtered<aod::DPGTracks> const& tracks, aod::DPGCollisions const&)
  {
    for (const auto& track : tracks) {
      histos.fill(HIST("Tracks/VertexPositionZ"), track.dpgCollision().posZ());
      histos.fill(HIST("Tracks/Kine/pt"), track.pt());
      histos.fill(HIST("Tracks/Kine/eta"), track.eta());
      histos.fill(HIST("Tracks/Kine/phi"), track.phi());
      histos.fill(HIST("Tracks/dcaXY"), track.dcaXY());
      histos.fill(HIST("Tracks/dcaZ"), track.dcaZ());
      histos.fill(HIST("Tracks/dcaXYvsPt"), track.pt(), track.dcaXY());
      histos.fill(HIST("Tracks/dcaZvsPt"), track.pt(), track.dcaZ());
      histos.fill(HIST("Tracks/length"), track.length());
      histos.fill(HIST("Tracks/ITS/itsChi2NCl"), track.itsChi2NCl());
      histos.fill(HIST("Tracks/ITS/itsNCl"), track.itsNCls());
      histos.fill(HIST("Tracks/ITS/itsNClvsItsHitmap"), track.itsNCls(), track.itsClusterMap());
      histos.fill(HIST("Tracks/TPC/tpcChi2NCl"), track.tpcChi2NCl());
      histos.fill(HIST("Tracks/TPC/tpcNClsFound"), track.tpcNClsFound());
      histos.fill(HIST("Tracks/TPC/tpcCrossedRows"), track.tpcNClsCrossedRows());
      histos.fill(HIST("Tracks/TPC/tpcCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("Tracks/TPC/tpcNClsFoundvsPt"), track.pt(), track.tpcNClsFound());
      histos.fill(HIST("Tracks/TPC/tpcCrossedRowsvsPt"), track.pt(), track.tpcNClsCrossedRows());
      histos.fill(HIST("Tracks/TPC/tpcCrossedRowsOverFindableClsvsPt"), track.pt(), track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("Tracks/TRD/trdChi2"), track.trdChi2());
      histos.fill(HIST("Tracks/TOF/tofChi2"), track.tofChi2());
      if (track.hasTPC()) {
        histos.fill(HIST("Tracks/matchedDet"), 1);
      }
      if (track.hasITS()) {
        histos.fill(HIST("Tracks/matchedDet"), 2);
      }
      if (track.hasTRD()) {
        histos.fill(HIST("Tracks/matchedDet"), 3);
      }
      if (track.hasTOF()) {
        histos.fill(HIST("Tracks/matchedDet"), 4);
      }
      histos.fill(HIST("Tracks/relativeResoPt"), track.pt(), track.ptReso());
      histos.fill(HIST("Tracks/relativeResoPtMean"), track.pt(), track.ptReso());
    }
  }
  PROCESS_SWITCH(qaEventTrackLite, processDataLite, "process data lite", true);

  // Process MC
  void processMCLite(o2::soa::Filtered<soa::Join<aod::DPGTracks, aod::DPGRecoParticles>> const& tracks, aod::DPGCollisions const&, aod::DPGNonRecoParticles const& particles)
  {
    for (const auto& track : tracks) {
      if (track.productionMode() == 0) {
        histos.get<TH1>(HIST("Particles/PDGs"))->Fill(Form("%i", track.pdgCode()), 1);
      }

      histos.fill(HIST("Particles/Kine/pt"), track.ptMC());
      histos.fill(HIST("Particles/Kine/eta"), track.etaMC());
      histos.fill(HIST("Particles/Kine/phi"), track.phiMC());

      histos.fill(HIST("Tracks/Kine/pt"), track.pt());
      histos.fill(HIST("Tracks/Kine/eta"), track.eta());
      histos.fill(HIST("Tracks/Kine/phi"), track.phi());
      histos.fill(HIST("Tracks/dcaXY"), track.dcaXY());
      histos.fill(HIST("Tracks/dcaZ"), track.dcaZ());
      histos.fill(HIST("Tracks/dcaXYvsPt"), track.pt(), track.dcaXY());
      histos.fill(HIST("Tracks/dcaZvsPt"), track.pt(), track.dcaZ());
      histos.fill(HIST("Tracks/length"), track.length());
      histos.fill(HIST("Tracks/ITS/itsChi2NCl"), track.itsChi2NCl());
      histos.fill(HIST("Tracks/TPC/tpcChi2NCl"), track.tpcChi2NCl());
      histos.fill(HIST("Tracks/TPC/tpcNClsFound"), track.tpcNClsFound());
      histos.fill(HIST("Tracks/TPC/tpcCrossedRows"), track.tpcNClsCrossedRows());
      histos.fill(HIST("Tracks/TPC/tpcCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("Tracks/TPC/tpcNClsFoundvsPt"), track.pt(), track.tpcNClsFound());
      histos.fill(HIST("Tracks/TPC/tpcCrossedRowsvsPt"), track.pt(), track.tpcNClsCrossedRows());
      histos.fill(HIST("Tracks/TPC/tpcCrossedRowsOverFindableClsvsPt"), track.pt(), track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("Tracks/TRD/trdChi2"), track.trdChi2());
      histos.fill(HIST("Tracks/TOF/tofChi2"), track.tofChi2());
      if (track.hasTPC()) {
        histos.fill(HIST("Tracks/matchedDet"), 1);
      }
      if (track.hasITS()) {
        histos.fill(HIST("Tracks/matchedDet"), 2);
      }
      if (track.hasTRD()) {
        histos.fill(HIST("Tracks/matchedDet"), 3);
      }
      if (track.hasTOF()) {
        histos.fill(HIST("Tracks/matchedDet"), 4);
      }
    }

    for (const auto& particle : particles) {
      histos.fill(HIST("Particles/Kine/pt"), particle.ptMC());
      histos.fill(HIST("Particles/Kine/eta"), particle.etaMC());
      histos.fill(HIST("Particles/Kine/phi"), particle.phiMC());
    }
  }
  PROCESS_SWITCH(qaEventTrackLite, processMCLite, "process MC lite", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qaEventTrackLite>(cfgc)};
}
