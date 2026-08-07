/// \file flattenicityTask.cxx
/// \brief Flattenicity analysis task for UE studies
/// \author Eisha Rani
/// \since August 2026

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

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct FlattenicityTask {

  // --- Flattenicity constants ---
  static constexpr int NCH_A = 96;
  static constexpr int NCH_C = 112;
  static constexpr int NCELL = NCH_A + NCH_C;
  static constexpr int NPHISECTORS = 8;
  static constexpr int NETA_A = NCH_A / NPHISECTORS;
  static constexpr int NETA_C = NCH_C / NPHISECTORS;

  // FT0 acceptance
  static constexpr float FT0A_ETA_MIN = 3.5f;
  static constexpr float FT0A_ETA_MAX = 4.9f;
  static constexpr float FT0C_ETA_MIN = -3.3f;
  static constexpr float FT0C_ETA_MAX = -2.1f;

  // --- Event selection constants ---
  static constexpr float VERTEX_CUT = 10.0f;
  static constexpr float FLAT_MIN = 0.0f;
  static constexpr float INEL_ETA_CUT = 1.0f;
  static constexpr float MIDRAP_ETA_CUT = 0.8f;

  // --- Configurables ---
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum pT"};

  Configurable<bool> isRun3{"isRun3", true, "is Run3 dataset"};
  Configurable<bool> timeEvsel{"timeEvsel", true, "TPC Time frame boundary cut"};
  Configurable<bool> piluprejection{"piluprejection", true, "Pileup rejection"};
  Configurable<bool> goodzvertex{"goodzvertex", true, "Good Z vertex"};

  // --- Track selection ---
  TrackSelection mySelectionPrim;

  // --- Histograms ---
  HistogramRegistry registry;

  // --- Init ---
  void init(InitContext const&) override
  {
    // Initialize track selection
    mySelectionPrim = myTrackSelectionPrim();

    // Define histograms
    AxisSpec flatBins = {40, 0.0, 1.0, "#rho"};
    AxisSpec nchBins = {100, -0.5, 99.5, "N_{ch}"};

    registry.add("hFlattenicityTruth", "Truth flattenicity; 1-#rho; Events",
                 HistType::kTH1D, {flatBins});
    registry.add("hFlattenicityReco", "Reco flattenicity; 1-#rho; Events",
                 HistType::kTH1D, {flatBins});
    registry.add("hFlattenicityCorrelation", "Truth vs Reco; 1-#rho_{truth}; 1-#rho_{reco}",
                 HistType::kTH2D, {flatBins, flatBins});
    registry.add("hNch", "Reco Nch distribution; N_{ch}; Events",
                 HistType::kTH1D, {nchBins});
    registry.add("hNchTruth", "Truth Nch distribution; N_{ch}; Events",
                 HistType::kTH1D, {nchBins});
  }

  // --- Track selection function ---
  TrackSelection myTrackSelectionPrim()
  {
    TrackSelection selectedTracks;
    selectedTracks.SetPtRange(0.1f, 1e10f);
    selectedTracks.SetEtaRange(-0.8f, 0.8f);
    selectedTracks.SetRequireITSRefit(true);
    selectedTracks.SetRequireTPCRefit(true);
    selectedTracks.SetMinNCrossedRowsTPC(70);
    selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.4f);
    selectedTracks.SetMaxChi2PerClusterTPC(4.0f);
    selectedTracks.SetRequireHitsInITSLayers(1, {0, 1});
    selectedTracks.SetMaxChi2PerClusterITS(36.0f);
    selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / std::pow(pt, 1.1f); });
    selectedTracks.SetMaxDcaZ(2.0f);
    return selectedTracks;
  }

  // --- Helper: Get cell ID for a particle in FT0 acceptance ---
  int getCellId(float eta, float phi)
  {
    // Check if in FT0-A acceptance
    if (eta > FT0A_ETA_MIN && eta < FT0A_ETA_MAX) {
      int phiBin = static_cast<int>(std::floor(phi / (2.0f * o2::constants::math::PI / static_cast<float>(NPHISECTORS))));
      phiBin = std::clamp(phiBin, 0, NPHISECTORS - 1);
      int etaBin = static_cast<int>(std::floor((eta - FT0A_ETA_MIN) / ((FT0A_ETA_MAX - FT0A_ETA_MIN) / static_cast<float>(NETA_A))));
      etaBin = std::clamp(etaBin, 0, NETA_A - 1);
      return etaBin * NPHISECTORS + phiBin;
    }

    // Check if in FT0-C acceptance
    if (eta > FT0C_ETA_MIN && eta < FT0C_ETA_MAX) {
      int phiBin = static_cast<int>(std::floor(phi / (2.0f * o2::constants::math::PI / static_cast<float>(NPHISECTORS))));
      phiBin = std::clamp(phiBin, 0, NPHISECTORS - 1);
      int etaBin = static_cast<int>(std::floor((eta - FT0C_ETA_MIN) / ((FT0C_ETA_MAX - FT0C_ETA_MIN) / static_cast<float>(NETA_C))));
      etaBin = std::clamp(etaBin, 0, NETA_C - 1);
      return NCH_A + etaBin * NPHISECTORS + phiBin;
    }

    return -1; // Not in FT0 acceptance
  }

  // --- Flattenicity calculation ---
  float calculateFlattenicity(const std::vector<float>& counts)
  {
    if (counts.size() != static_cast<size_t>(NCELL)) {
      return -1.0f;
    }

    float total = 0.0f;
    for (const auto& c : counts) {
      total += c;
    }
    if (total <= 0.0f) {
      return -1.0f;
    }

    float mean = total / static_cast<float>(NCELL);
    if (mean <= 0.0f) {
      return -1.0f;
    }

    float sumSq = 0.0f;
    for (const auto& c : counts) {
      sumSq += (c - mean) * (c - mean);
    }

    float rho = std::sqrt(sumSq / (static_cast<float>(NCELL) * static_cast<float>(NCELL))) / mean;
    return 1.0f - rho;
  }

  // --- Process Data ---
  void processData(aod::Collision const& collision,
                   soa::Filtered<aod::Tracks> const& tracks,
                   aod::FT0s const& ft0s)
  {
    // Event selection (Paola/Jesus)
    if (!collision.sel8()) {
      return;
    }
    if (std::abs(collision.posZ()) >= VERTEX_CUT) {
      return;
    }

    // Track loop for Nch
    int nch = 0;
    for (const auto& track : tracks) {
      if (!mySelectionPrim.IsSelected(track)) {
        continue;
      }
      nch++;
    }
    registry.fill(HIST("hNch"), nch);

    // FT0 flattenicity
    auto ft0 = collision.ft0();
    if (ft0.hasAmplitudeA() && ft0.hasAmplitudeC()) {
      auto ampA = ft0.amplitudeA();
      auto ampC = ft0.amplitudeC();

      std::vector<float> counts(NCELL, 0.0f);
      for (int i = 0; i < static_cast<int>(ampA.size()) && i < NCH_A; ++i) {
        counts[i] = ampA[i];
      }
      for (int i = 0; i < static_cast<int>(ampC.size()) && i < NCH_C; ++i) {
        counts[NCH_A + i] = ampC[i];
      }

      float flat = calculateFlattenicity(counts);
      if (flat >= FLAT_MIN) {
        registry.fill(HIST("hFlattenicityReco"), flat);
      }
    }
  }
  PROCESS_SWITCH(FlattenicityTask, processData, "Process data", true);

  // --- Process MC ---
  void processMC(aod::McCollision const& mcCollision,
                 aod::McParticles const& particles,
                 soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions,
                 aod::FT0s const& ft0s,
                 aod::BCs const& /*bcs*/)
  {
    // ---- Truth-level processing ----
    bool inel = false;
    int nchTruth = 0;
    std::vector<float> truthCounts(NCELL, 0.0f);

    for (const auto& particle : particles) {
      // Check if physical primary
      if (!particle.isPhysicalPrimary()) {
        continue;
      }

      // Check if charged
      if (std::abs(particle.pdgCode()) == 0) {
        continue;
      }

      // Check pT > 0
      if (particle.pt() <= 0.0f) {
        continue;
      }

      // INEL>0 check: primary charged with |eta| < INEL_ETA_CUT
      if (std::abs(particle.eta()) < INEL_ETA_CUT) {
        inel = true;
      }

      // Nch at midrapidity: |eta| < MIDRAP_ETA_CUT, pT > cfgTrkLowPtCut
      if (std::abs(particle.eta()) < MIDRAP_ETA_CUT && particle.pt() > cfgTrkLowPtCut) {
        nchTruth++;
      }

      // Flattenicity: particles in FT0 acceptance
      int cellId = getCellId(particle.eta(), particle.phi());
      if (cellId >= 0 && cellId < NCELL) {
        truthCounts[cellId] += 1.0f;
      }
    }

    // Apply truth-level event selection (Paola/Jesus)
    if (!inel) {
      return;
    }
    if (std::abs(mcCollision.posZ()) >= VERTEX_CUT) {
      return;
    }

    // Fill truth multiplicity
    registry.fill(HIST("hNchTruth"), nchTruth);
    registry.fill(HIST("hNch"), nchTruth);

    // Calculate truth flattenicity
    float truthFlat = calculateFlattenicity(truthCounts);
    if (truthFlat >= FLAT_MIN) {
      registry.fill(HIST("hFlattenicityTruth"), truthFlat);
    }

    // ---- Reconstructed-level processing for matched collisions ----
    for (const auto& collision : collisions) {
      // Apply reconstruction-level event selection
      if (!collision.sel8()) {
        continue;
      }
      if (std::abs(collision.posZ()) >= VERTEX_CUT) {
        continue;
      }

      // Get FT0 flattenicity for this collision
      auto ft0 = collision.ft0();
      if (!ft0.hasAmplitudeA() || !ft0.hasAmplitudeC()) {
        continue;
      }

      auto ampA = ft0.amplitudeA();
      auto ampC = ft0.amplitudeC();

      std::vector<float> recoCounts(NCELL, 0.0f);
      for (int i = 0; i < static_cast<int>(ampA.size()) && i < NCH_A; ++i) {
        recoCounts[i] = ampA[i];
      }
      for (int i = 0; i < static_cast<int>(ampC.size()) && i < NCH_C; ++i) {
        recoCounts[NCH_A + i] = ampC[i];
      }

      float recoFlat = calculateFlattenicity(recoCounts);
      if (recoFlat >= FLAT_MIN && truthFlat >= FLAT_MIN) {
        registry.fill(HIST("hFlattenicityReco"), recoFlat);
        registry.fill(HIST("hFlattenicityCorrelation"), truthFlat, recoFlat);
      }
    }
  }
  PROCESS_SWITCH(FlattenicityTask, processMC, "Process MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<FlattenicityTask>(cfgc));
  return workflow;
}
