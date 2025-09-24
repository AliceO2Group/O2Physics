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

/// \file systematicsMapping.cxx
/// \brief Task to perform a systematics study for K0s and charged Kaons
/// \author Nicolò Jacazio, Universita del Piemonte Orientale (IT)
/// \since September 22, 2025

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <algorithm>

using namespace o2;
using namespace o2::framework;

struct SystematicsMapping {
  // Returns a unique index for the combination of cuts
  ConfigurableAxis ptBins{"ptBins", {100, 0.f, 10.f}, "Binning for pT (GeV/c)"};
  ConfigurableAxis etaBins{"etaBins", {40, -1.0f, 1.0f}, "Binning for #eta"};
  ConfigurableAxis phiBins{"phiBins", {36, 0.f, o2::constants::math::TwoPI}, "Binning for #phi (rad)"};
  // Define the Signal axis
  ConfigurableAxis invariantMassBins{"invariantMassBins", {100, -0.1f, 0.1f}, "Binning for the invariant mass (GeV/c^2)"};
  ConfigurableAxis nsigmaBins{"nsigmaBins", {100, -10.f, 10.f}, "Binning for nSigma"};
  // Selection bins
  ConfigurableAxis tpcCrossedRowsBins{"tpcCrossedRowsBins", {5, 70, 100, 120, 135, 150}, "Min TPC clusters for tracks"};
  ConfigurableAxis itsClustersBins{"itsClustersBins", {5, 0, 6}, "Min ITS clusters for tracks"};
  ConfigurableAxis dcaBins{"dcaBins", {100, 0.f, 5.f}, "Binning for DCA (cm)"};
  ConfigurableAxis chi2Bins{"chi2Bins", {100, 0.f, 100.f}, "Binning for chi2"};
  // Selection configurables
  Configurable<float> selectionPosZ{"selectionPosZ", 10.f, "Max |z| of the primary vertex"};
  // V0 selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 10, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.0, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.0, "DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 0.0, "Radius"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};

  HistogramRegistry registry{"registry"};

  template <typename T>
  bool isCollisionSelected(T const& collision)
  {
    return collision.sel8() && std::abs(collision.posZ()) <= selectionPosZ;
  }

  void init(InitContext const&)
  {
    const AxisSpec ptAxis{ptBins, "#it{p}_{T} (GeV/c)"};
    const AxisSpec etaAxis{etaBins, "#eta"};
    const AxisSpec phiAxis{phiBins, "#phi (rad)"};
    const AxisSpec invariantMassAxis{invariantMassBins, "Invariant Mass (GeV/c^{2})"};
    const AxisSpec nsigmaAxisTPC{nsigmaBins, "nSigma TPC"};
    const AxisSpec nsigmaAxisTOF{nsigmaBins, "nSigma TOF"};
    const AxisSpec tpcCrossedRowsAxis{tpcCrossedRowsBins, "TPC crossed rows"};
    const AxisSpec itsClustersAxis{itsClustersBins, "ITS clusters"};
    const AxisSpec dcaXYAxis{dcaBins, "DCAxy (cm)"};
    const AxisSpec dcaZAxis{dcaBins, "DCAz (cm)"};
    const AxisSpec chi2TPCAxis{chi2Bins, "TPC Chi2"};
    const AxisSpec chi2ITSAxis{chi2Bins, "ITS Chi2"};

    if (doprocessData) {

      // First we define the histograms on which we are cutting (tpc clusters, its clusters, ..)
      registry.add("K/hTPCCrossedRows", "", HistType::kTH1F, {{100, 0, 200}});
      registry.add("K/hITSClusters", "", HistType::kTH1F, {{10, 0, 10}});
      registry.add("K/hDCAxy", "", HistType::kTH1F, {dcaBins});
      registry.add("K/hDCAz", "", HistType::kTH1F, {dcaBins});
      registry.add("K/hChi2OverNCLsTPC", "", HistType::kTH1F, {chi2Bins});
      registry.add("K/hChi2OverNCLsITS", "", HistType::kTH1F, {chi2Bins});
      registry.addClone("K/", "K0s/");

      // Add the signal histograms
      registry.add("K/SignalPositive", "", HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis, nsigmaAxisTPC, nsigmaAxisTOF, tpcCrossedRowsAxis, itsClustersAxis, dcaXYAxis, dcaZAxis, chi2TPCAxis, chi2ITSAxis});
      registry.add("K/SignalNegative", "", HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis, nsigmaAxisTPC, nsigmaAxisTOF, tpcCrossedRowsAxis, itsClustersAxis, dcaXYAxis, dcaZAxis, chi2TPCAxis, chi2ITSAxis});
      registry.add("K0s/Signal", "", HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis, invariantMassBins, nsigmaAxisTPC, nsigmaAxisTOF, tpcCrossedRowsAxis, itsClustersAxis, dcaXYAxis, dcaZAxis, chi2TPCAxis, chi2ITSAxis});
    }

    if (doprocessMc) {
      registry.add("K/GeneratedPositive", "", HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis});
      registry.add("K/GeneratedNegative", "", HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis});
      registry.add("K0s/Generated", "", HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis});
    }
  }

  using TrackType = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::TracksDCA>;
  using CollisionType = soa::Join<aod::Collisions, aod::EvSels>;

  void processData(CollisionType const& collisions,
                   TrackType const& tracks,
                   aod::V0Datas const& v0s)
  {
    for (const auto& collision : collisions) {
      if (isCollisionSelected(collision))
        continue; // MB selection

      // Kaon loop
      for (const auto& track : tracks) {
        if (track.collisionId() != collision.globalIndex())
          continue;
        registry.fill(HIST("K/hTPCCrossedRows"), track.tpcNClsCrossedRows());
        registry.fill(HIST("K/hITSClusters"), track.itsNCls());
        registry.fill(HIST("K/hDCAxy"), track.dcaXY());
        registry.fill(HIST("K/hDCAz"), track.dcaZ());
        registry.fill(HIST("K/hChi2OverNCLsTPC"), track.tpcChi2NCl());
        registry.fill(HIST("K/hChi2OverNCLsITS"), track.itsChi2NCl());
        if (track.sign() > 0)
          registry.fill(HIST("K/SignalPositive"), track.pt(), track.eta(), track.phi(), track.tpcNSigmaKa(), track.tofNSigmaKa(), track.tpcNClsCrossedRows(), track.itsNCls(), track.dcaXY(), track.dcaZ(), track.tpcChi2NCl(), track.itsChi2NCl());
        else
          registry.fill(HIST("K/SignalNegative"), track.pt(), track.eta(), track.phi(), track.tpcNSigmaKa(), track.tofNSigmaKa(), track.tpcNClsCrossedRows(), track.itsNCls(), track.dcaXY(), track.dcaZ(), track.tpcChi2NCl(), track.itsChi2NCl());
      }

      // K0s loop
      for (const auto& v0 : v0s) {
        if (v0.collisionId() != collision.globalIndex())
          continue;

        const auto& posTrack = v0.posTrack_as<TrackType>();
        const auto& negTrack = v0.negTrack_as<TrackType>();
        if (v0.v0radius() < v0radius ||
            v0.v0cosPA() < v0cospa ||
            std::abs(posTrack.eta()) > etadau ||
            std::abs(negTrack.eta()) > etadau)
          continue;
        registry.fill(HIST("K0s/hTPCCrossedRows"), std::min(posTrack.tpcNClsCrossedRows(), negTrack.tpcNClsCrossedRows()));
        registry.fill(HIST("K0s/hITSClusters"), std::min(posTrack.itsNCls(), negTrack.itsNCls()));
        registry.fill(HIST("K0s/hDCAxy"), std::min(posTrack.dcaXY(), negTrack.dcaXY()));
        registry.fill(HIST("K0s/hDCAz"), std::min(posTrack.dcaZ(), negTrack.dcaZ()));
        registry.fill(HIST("K0s/hChi2OverNCLsTPC"), std::min(posTrack.tpcChi2NCl(), negTrack.tpcChi2NCl()));
        registry.fill(HIST("K0s/hChi2OverNCLsITS"), std::min(posTrack.itsChi2NCl(), negTrack.itsChi2NCl()));
        registry.fill(HIST("K0s/Signal"), v0.pt(), v0.eta(), v0.phi(), v0.mK0Short() - constants::physics::MassK0Short,
                      std::max(posTrack.tpcNSigmaPi(), negTrack.tpcNSigmaPi()),
                      std::max(posTrack.tofNSigmaPi(), negTrack.tofNSigmaPi()),
                      std::min(posTrack.tpcNClsCrossedRows(), negTrack.tpcNClsCrossedRows()),
                      std::min(posTrack.itsNCls(), negTrack.itsNCls()),
                      std::min(posTrack.dcaXY(), negTrack.dcaXY()),
                      std::min(posTrack.dcaZ(), negTrack.dcaZ()),
                      std::min(posTrack.tpcChi2NCl(), negTrack.tpcChi2NCl()),
                      std::min(posTrack.itsChi2NCl(), negTrack.itsChi2NCl()));
      }
    }
  }
  PROCESS_SWITCH(SystematicsMapping, processData, "Systematics study for K0s and charged Kaons", true);

  void processMc(soa::Join<CollisionType, aod::McCollisionLabels> const& collisions,
                 aod::McParticles const& particles,
                 aod::McCollisions const&)
  {
    for (const auto& collision : collisions) {
      if (!isCollisionSelected(collision))
        continue; // MB selection
      if (!collision.has_mcCollision())
        continue;
      const auto& mcCollision = collision.mcCollision();

      for (const auto& particle : particles) {
        if (particle.mcCollisionId() != mcCollision.globalIndex())
          continue;
        if (!particle.isPhysicalPrimary())
          continue;
        switch (particle.pdgCode()) {
          case 321: // K+
            registry.fill(HIST("K/GeneratedPositive"), particle.pt(), particle.eta(), particle.phi());
            break;
          case -321: // K-
            registry.fill(HIST("K/GeneratedNegative"), particle.pt(), particle.eta(), particle.phi());
            break;
          case 310: // K0s
            registry.fill(HIST("K0s/Generated"), particle.pt(), particle.eta(), particle.phi());
            break;
          default:
            break;
        }
      }
    }
  }
  PROCESS_SWITCH(SystematicsMapping, processMc, "Systematics study for K0s and charged Kaons on MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<SystematicsMapping>(cfgc)}; }
