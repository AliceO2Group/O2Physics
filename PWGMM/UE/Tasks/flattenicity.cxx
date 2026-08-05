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

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct flattenicityTask {

  // --- Flattenicity constants ---
  static constexpr int N_CH_A = 96;
  static constexpr int N_CH_C = 112;
  static constexpr int N_CELL = N_CH_A + N_CH_C;
  static constexpr int N_PHI_SECTORS = 8;
  static constexpr int N_ETA_A = N_CH_A / N_PHI_SECTORS;
  static constexpr int N_ETA_C = N_CH_C / N_PHI_SECTORS;

  // FT0 acceptance
  static constexpr float FT0A_ETA_MIN = 3.5;
  static constexpr float FT0A_ETA_MAX = 4.9;
  static constexpr float FT0C_ETA_MIN = -3.3;
  static constexpr float FT0C_ETA_MAX = -2.1;

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
  void init(InitContext const&)
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
    registry.add("hNch", "Nch distribution; N_{ch}; Events",
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
    selectedTracks.SetMaxChi2PerClusterTPC(4.f);
    selectedTracks.SetRequireHitsInITSLayers(1, {0, 1});
    selectedTracks.SetMaxChi2PerClusterITS(36.f);
    selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / pow(pt, 1.1f); });
    selectedTracks.SetMaxDcaZ(2.f);
    return selectedTracks;
  }

  // --- Flattenicity calculation ---
  float calculateFlattenicity(const std::vector<float>& counts)
  {
    if (counts.size() != N_CELL)
      return -1;

    float total = 0;
    for (auto c : counts)
      total += c;
    if (total <= 0)
      return -1;

    float mean = total / N_CELL;
    if (mean <= 0)
      return -1;

    float sumSq = 0;
    for (auto c : counts) {
      sumSq += (c - mean) * (c - mean);
    }

    float rho = std::sqrt(sumSq / (N_CELL * N_CELL)) / mean;
    return 1.0f - rho;
  }

  // --- Process Data ---
  void processData(aod::Collision const& collision,
                   soa::Filtered<aod::Tracks> const& tracks,
                   aod::FT0s const& ft0s)
  {
    // Event selection (Paola/Jesus)
    if (!collision.sel8())
      return;
    if (std::abs(collision.posZ()) >= 10.0f)
      return;

    // Track loop for Nch
    int nch = 0;
    for (auto& track : tracks) {
      if (!mySelectionPrim.IsSelected(track))
        continue;
      nch++;
    }
    registry.fill(HIST("hNch"), nch);

    // FT0 flattenicity
    auto ft0 = collision.ft0();
    if (ft0.hasAmplitudeA() && ft0.hasAmplitudeC()) {
      auto ampA = ft0.amplitudeA();
      auto ampC = ft0.amplitudeC();

      std::vector<float> counts(N_CELL, 0.0f);
      for (int i = 0; i < ampA.size() && i < N_CH_A; ++i) {
        counts[i] = ampA[i];
      }
      for (int i = 0; i < ampC.size() && i < N_CH_C; ++i) {
        counts[N_CH_A + i] = ampC[i];
      }

      float flat = calculateFlattenicity(counts);
      if (flat >= 0) {
        registry.fill(HIST("hFlattenicityReco"), flat);
      }
    }
  }
  PROCESS_SWITCH(flattenicityTask, processData, "Process data", true);

  // --- Process MC ---
  void processMC(aod::McCollision const& mcCollision,
                 aod::McParticles const& particles,
                 soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions,
                 aod::FT0s const& ft0s)
  {
    // MC processing - to be implemented
  }
  PROCESS_SWITCH(flattenicityTask, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<flattenicityTask>(cfgc));
  return workflow;
}
