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
  ConfigurableAxis dcaBins{"dcaBins", {VARIABLE_WIDTH, -3.0, -2.95, -2.9, -2.85, -2.8, -2.75, -2.7, -2.65, -2.6, -2.55, -2.5, -2.45, -2.4, -2.35, -2.3, -2.25, -2.2, -2.15, -2.1, -2.05, -2.0, -1.975, -1.95, -1.925, -1.9, -1.875, -1.85, -1.825, -1.8, -1.775, -1.75, -1.725, -1.7, -1.675, -1.65, -1.625, -1.6, -1.575, -1.55, -1.525, -1.5, -1.475, -1.45, -1.425, -1.4, -1.375, -1.35, -1.325, -1.3, -1.275, -1.25, -1.225, -1.2, -1.175, -1.15, -1.125, -1.1, -1.075, -1.05, -1.025, -1.0, -0.99, -0.98, -0.97, -0.96, -0.95, -0.94, -0.93, -0.92, -0.91, -0.9, -0.89, -0.88, -0.87, -0.86, -0.85, -0.84, -0.83, -0.82, -0.81, -0.8, -0.79, -0.78, -0.77, -0.76, -0.75, -0.74, -0.73, -0.72, -0.71, -0.7, -0.69, -0.68, -0.67, -0.66, -0.65, -0.64, -0.63, -0.62, -0.61, -0.6, -0.59, -0.58, -0.57, -0.56, -0.55, -0.54, -0.53, -0.52, -0.51, -0.5, -0.49, -0.48, -0.47, -0.46, -0.45, -0.44, -0.43, -0.42, -0.41, -0.4, -0.396, -0.392, -0.388, -0.384, -0.38, -0.376, -0.372, -0.368, -0.364, -0.36, -0.356, -0.352, -0.348, -0.344, -0.34, -0.336, -0.332, -0.328, -0.324, -0.32, -0.316, -0.312, -0.308, -0.304, -0.3, -0.296, -0.292, -0.288, -0.284, -0.28, -0.276, -0.272, -0.268, -0.264, -0.26, -0.256, -0.252, -0.248, -0.244, -0.24, -0.236, -0.232, -0.228, -0.224, -0.22, -0.216, -0.212, -0.208, -0.204, -0.2, -0.198, -0.196, -0.194, -0.192, -0.19, -0.188, -0.186, -0.184, -0.182, -0.18, -0.178, -0.176, -0.174, -0.172, -0.17, -0.168, -0.166, -0.164, -0.162, -0.16, -0.158, -0.156, -0.154, -0.152, -0.15, -0.148, -0.146, -0.144, -0.142, -0.14, -0.138, -0.136, -0.134, -0.132, -0.13, -0.128, -0.126, -0.124, -0.122, -0.12, -0.118, -0.116, -0.114, -0.112, -0.11, -0.108, -0.106, -0.104, -0.102, -0.1, -0.099, -0.098, -0.097, -0.096, -0.095, -0.094, -0.093, -0.092, -0.091, -0.09, -0.089, -0.088, -0.087, -0.086, -0.085, -0.084, -0.083, -0.082, -0.081, -0.08, -0.079, -0.078, -0.077, -0.076, -0.075, -0.074, -0.073, -0.072, -0.071, -0.07, -0.069, -0.068, -0.067, -0.066, -0.065, -0.064, -0.063, -0.062, -0.061, -0.06, -0.059, -0.058, -0.057, -0.056, -0.055, -0.054, -0.053, -0.052, -0.051, -0.05, -0.049, -0.048, -0.047, -0.046, -0.045, -0.044, -0.043, -0.042, -0.041, -0.04, -0.039, -0.038, -0.037, -0.036, -0.035, -0.034, -0.033, -0.032, -0.031, -0.03, -0.029, -0.028, -0.027, -0.026, -0.025, -0.024, -0.023, -0.022, -0.021, -0.02, -0.019, -0.018, -0.017, -0.016, -0.015, -0.014, -0.013, -0.012, -0.011, -0.01, -0.009, -0.008, -0.007, -0.006, -0.005, -0.004, -0.003, -0.002, -0.001, -0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.04, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.05, 0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 0.058, 0.059, 0.06, 0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.07, 0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099, 0.1, 0.102, 0.104, 0.106, 0.108, 0.11, 0.112, 0.114, 0.116, 0.118, 0.12, 0.122, 0.124, 0.126, 0.128, 0.13, 0.132, 0.134, 0.136, 0.138, 0.14, 0.142, 0.144, 0.146, 0.148, 0.15, 0.152, 0.154, 0.156, 0.158, 0.16, 0.162, 0.164, 0.166, 0.168, 0.17, 0.172, 0.174, 0.176, 0.178, 0.18, 0.182, 0.184, 0.186, 0.188, 0.19, 0.192, 0.194, 0.196, 0.198, 0.2, 0.204, 0.208, 0.212, 0.216, 0.22, 0.224, 0.228, 0.232, 0.236, 0.24, 0.244, 0.248, 0.252, 0.256, 0.26, 0.264, 0.268, 0.272, 0.276, 0.28, 0.284, 0.288, 0.292, 0.296, 0.3, 0.304, 0.308, 0.312, 0.316, 0.32, 0.324, 0.328, 0.332, 0.336, 0.34, 0.344, 0.348, 0.352, 0.356, 0.36, 0.364, 0.368, 0.372, 0.376, 0.38, 0.384, 0.388, 0.392, 0.396, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0}, "Binning for DCA (cm)"};
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
      registry.add("K/ReconstructedPositive", "", HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis, nsigmaAxisTPC, nsigmaAxisTOF, tpcCrossedRowsAxis, itsClustersAxis, dcaXYAxis, dcaZAxis, chi2TPCAxis, chi2ITSAxis});
      registry.add("K/ReconstructedNegative", "", HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis, nsigmaAxisTPC, nsigmaAxisTOF, tpcCrossedRowsAxis, itsClustersAxis, dcaXYAxis, dcaZAxis, chi2TPCAxis, chi2ITSAxis});
      registry.add("K0s/Reconstructed", "", HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis, invariantMassBins, nsigmaAxisTPC, nsigmaAxisTOF, tpcCrossedRowsAxis, itsClustersAxis, dcaXYAxis, dcaZAxis, chi2TPCAxis, chi2ITSAxis});
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
                 soa::Join<TrackType, aod::McTrackLabels> const& tracks,
                 soa::Join<aod::V0Datas, aod::McV0Labels> const& v0s,
                 aod::McParticles const& particles,
                 aod::McCollisions const&)
  {
    for (const auto& collision : collisions) {
      if (!isCollisionSelected(collision))
        continue; // MB selection
      if (!collision.has_mcCollision())
        continue;
      const auto& mcCollision = collision.mcCollision();

      for (const auto& track : tracks) {
        if (track.collisionId() != collision.globalIndex())
          continue;
        if (!track.has_mcParticle())
          continue;
        const auto& mcParticle = track.mcParticle();
        if (mcParticle.mcCollisionId() != mcCollision.globalIndex())
          continue;
        if (!mcParticle.isPhysicalPrimary())
          continue;
        switch (mcParticle.pdgCode()) {
          case 321: // K+
            registry.fill(HIST("K/ReconstructedPositive"), track.pt(), track.eta(), track.phi(), track.tpcNSigmaKa(), track.tofNSigmaKa(), track.tpcNClsCrossedRows(), track.itsNCls(), track.dcaXY(), track.dcaZ(), track.tpcChi2NCl(), track.itsChi2NCl());
            break;
          case -321: // K-
            registry.fill(HIST("K/ReconstructedNegative"), track.pt(), track.eta(), track.phi(), track.tpcNSigmaKa(), track.tofNSigmaKa(), track.tpcNClsCrossedRows(), track.itsNCls(), track.dcaXY(), track.dcaZ(), track.tpcChi2NCl(), track.itsChi2NCl());
            break;
          default:
            break;
        }
      }

      for (const auto& v0 : v0s) {
        if (v0.collisionId() != collision.globalIndex())
          continue;
        if (!v0.has_mcParticle())
          continue;
        const auto& mcParticle = v0.mcParticle();
        if (mcParticle.mcCollisionId() != mcCollision.globalIndex())
          continue;
        if (!mcParticle.isPhysicalPrimary())
          continue;
        if (std::abs(mcParticle.pdgCode()) != 310)
          continue;
        const auto& posTrack = v0.posTrack_as<TrackType>();
        const auto& negTrack = v0.negTrack_as<TrackType>();
        if (v0.v0radius() < v0radius ||
            v0.v0cosPA() < v0cospa ||
            std::abs(posTrack.eta()) > etadau ||
            std::abs(negTrack.eta()) > etadau)
          continue;

        registry.fill(HIST("K0s/Reconstructed"), v0.pt(), v0.eta(), v0.phi(), v0.mK0Short() - constants::physics::MassK0Short,
                      std::max(posTrack.tpcNSigmaPi(), negTrack.tpcNSigmaPi()),
                      std::max(posTrack.tofNSigmaPi(), negTrack.tofNSigmaPi()),
                      std::min(posTrack.tpcNClsCrossedRows(), negTrack.tpcNClsCrossedRows()),
                      std::min(posTrack.itsNCls(), negTrack.itsNCls()),
                      std::min(posTrack.dcaXY(), negTrack.dcaXY()),
                      std::min(posTrack.dcaZ(), negTrack.dcaZ()),
                      std::min(posTrack.tpcChi2NCl(), negTrack.tpcChi2NCl()),
                      std::min(posTrack.itsChi2NCl(), negTrack.itsChi2NCl()));
      }

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
