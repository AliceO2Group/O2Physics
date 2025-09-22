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
  ConfigurableAxis tpcClusterBins{"tpcClusterBins", {5, 70, 100, 120, 135, 150}, "Min TPC clusters for tracks"};
  ConfigurableAxis itsClustersBins{"itsClustersBins", {5, 0, 6}, "Min ITS clusters for tracks"};
  // Selection configurables
  Configurable<float> selectionPosZ{"selectionPosZ", 10.f, "Max |z| of the primary vertex"};

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

    if (doprocessData) {

      // First we define the histograms on which we are cutting (tpc clusters, its clusters, ..)
      registry.add("K/hTPCClusters", "", HistType::kTH1F, {{100, 0, 200}});
      registry.add("K/hITSClusters", "", HistType::kTH1F, {{10, 0, 10}});
      registry.addClone("K/", "K0s/");

      // Add the signal histograms
      registry.add("K/SignalPositive", "", HistType::kTHnSparseF, {ptBins, etaBins, phiBins, nsigmaBins, tpcClusterBins, itsClustersBins});
      registry.add("K/SignalNegative", "", HistType::kTHnSparseF, {ptBins, etaBins, phiBins, nsigmaBins, tpcClusterBins, itsClustersBins});
      registry.add("K0s/Signal", "", HistType::kTHnSparseF, {ptBins, etaBins, phiBins, invariantMassBins, tpcClusterBins, itsClustersBins});
    }

    if (doprocessMc) {
      registry.add("K/GeneratedPositive", "", HistType::kTHnSparseF, {ptBins, etaBins, phiBins});
      registry.add("K/GeneratedNegative", "", HistType::kTHnSparseF, {ptBins, etaBins, phiBins});
      registry.add("K0s/Generated", "", HistType::kTHnSparseF, {ptBins, etaBins, phiBins});
    }
  }

  using TrackType = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullKa, aod::pidTOFFullPi>;
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
        registry.fill(HIST("hTPCClusters"), track.tpcNClsFound());
        registry.fill(HIST("hITSClusters"), track.itsNCls());
        if (track.sign() > 0)
          registry.fill(HIST("K/SignalPositive"), track.pt(), track.eta(), track.phi(), track.tpcNSigmaKa(), track.tpcNClsFound(), track.itsNCls());
        else
          registry.fill(HIST("K/SignalNegative"), track.pt(), track.eta(), track.phi(), track.tpcNSigmaKa(), track.tpcNClsFound(), track.itsNCls());
      }

      // K0s loop
      for (const auto& v0 : v0s) {
        if (v0.collisionId() != collision.globalIndex())
          continue;
        const auto& posTrack = v0.posTrack_as<TrackType>();
        const auto& negTrack = v0.negTrack_as<TrackType>();
        registry.fill(HIST("K0s/Signal"), v0.pt(), v0.eta(), v0.phi(), v0.mK0Short() - constants::physics::MassK0Short, std::min(posTrack.tpcNClsFound(), negTrack.tpcNClsFound()), std::min(posTrack.itsNCls(), negTrack.itsNCls()));
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
