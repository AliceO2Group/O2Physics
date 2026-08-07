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

/// \file flattenicityTask.cxx
/// \brief Flattenicity analysis task for UE studies
/// \author Eisha Rani
/// \since August 2026

#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

struct FlattenicityTask {

  // ============================================
  // FT0 Constants
  // ============================================
  static constexpr int N_PHI_SECTORS = 8;
  static constexpr int N_ETA_A = 12;
  static constexpr int N_ETA_C = 14;
  static constexpr int N_CH_A = 96;
  static constexpr int N_CH_C = 112;
  static constexpr int N_CELL = N_CH_A + N_CH_C; // 208

  static constexpr float FT0A_ETA_MIN = 3.5;
  static constexpr float FT0A_ETA_MAX = 4.9;
  static constexpr float FT0C_ETA_MIN = -3.3;
  static constexpr float FT0C_ETA_MAX = -2.1;

  static constexpr int PHYSICAL_PRIMARY_BIT = 0x4;

  // ============================================
  // Histogram Definitions
  // ============================================
  HistogramRegistry histos{
    "histos",
    {
      // Event-level
      {"hEvents", "Event selection;;Counts", {HistType::kTH1F, {{5, 0, 5}}}},

      // dNch/deta
      {"hNch_INEL", "Nch distribution (INEL>0);N_{ch};Entries", {HistType::kTH1F, {{100, -0.5, 99.5}}}},
      {"hNch_FT0", "Nch distribution (INEL>0 & FT0);N_{ch};Entries", {HistType::kTH1F, {{100, -0.5, 99.5}}}},

      // Flattenicity
      {"hFlattenicity", "Flattenicity distribution;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicity_vs_Nch", "Flattenicity vs Nch;N_{ch};1-#rho", {HistType::kTH2F, {{50, -0.5, 99.5}, {50, 0.0, 1.0}}}},

      // FT0 cell occupancy
      {"hCellOccupancy", "FT0 cell occupancy;Cell ID;Entries", {HistType::kTH1F, {{N_CELL, 0, N_CELL}}}},
      {"hCellOccupancyFT0A", "FT0-A cell occupancy;Cell ID;Entries", {HistType::kTH1F, {{N_CH_A, 0, N_CH_A}}}},
      {"hCellOccupancyFT0C", "FT0-C cell occupancy;Cell ID;Entries", {HistType::kTH1F, {{N_CH_C, 0, N_CH_C}}}},

      // Multiplicity classes
      {"hFlattenicity_0_10", "Flattenicity 0-10%;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicity_10_20", "Flattenicity 10-20%;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicity_20_30", "Flattenicity 20-30%;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicity_30_40", "Flattenicity 30-40%;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicity_40_50", "Flattenicity 40-50%;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicity_50_60", "Flattenicity 50-60%;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicity_60_70", "Flattenicity 60-70%;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicity_70_80", "Flattenicity 70-80%;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicity_80_90", "Flattenicity 80-90%;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicity_90_100", "Flattenicity 90-100%;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
    }};

  // ============================================
  // Configurables
  // ============================================
  Configurable<float> cfgPtMin{"cfgPtMin", 0.1, "Minimum pT for tracks"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8, "Maximum |eta| for tracks"};
  Configurable<float> cfgVzMax{"cfgVzMax", 10.0, "Maximum |vz| for collisions"};
  Configurable<int> cfgNCrossedRowsTPC{"cfgNCrossedRowsTPC", 70, "Minimum TPC crossed rows"};
  Configurable<float> cfgChi2PerClusterTPC{"cfgChi2PerClusterTPC", 4.0, "Maximum TPC chi2 per cluster"};
  Configurable<float> cfgChi2PerClusterITS{"cfgChi2PerClusterITS", 36.0, "Maximum ITS chi2 per cluster"};
  Configurable<float> cfgDCAZ{"cfgDCAZ", 0.1, "Maximum DCA z"};
  Configurable<bool> cfgRequireGoldenChi2{"cfgRequireGoldenChi2", true, "Require golden chi2"};

  // ============================================
  // Particle charge function
  // ============================================
  int getCharge(int pdgCode)
  {
    switch (std::abs(pdgCode)) {
      case 211:  // pion
      case 321:  // kaon
      case 2212: // proton
        return 1;
      case 11: // electron
      case 13: // muon
        return -1;
      default:
        return 0;
    }
  }

  // ============================================
  // Flattenicity calculation
  // ============================================
  float computeFlattenicity(const std::array<float, N_CELL>& counts)
  {
    float total = 0.0;
    for (int i = 0; i < N_CELL; i++) {
      total += counts[i];
    }

    if (total <= 0)
      return -1.0;

    float mean = total / N_CELL;
    if (mean <= 0)
      return -1.0;

    float sumSq = 0.0;
    for (int i = 0; i < N_CELL; i++) {
      sumSq += (counts[i] - mean) * (counts[i] - mean);
    }

    float rho = std::sqrt(sumSq) / (N_CELL * mean);
    return rho;
  }

  // ============================================
  // Assign particle to FT0 cell
  // ============================================
  int assignToFT0Cell(float eta, float phi, bool& isFT0A)
  {
    // Check if in FT0 acceptance
    bool inFT0A = (eta > FT0A_ETA_MIN && eta < FT0A_ETA_MAX);
    bool inFT0C = (eta > FT0C_ETA_MIN && eta < FT0C_ETA_MAX);

    if (!inFT0A && !inFT0C)
      return -1;

    isFT0A = inFT0A;

    // Phi bin
    int phiBin = static_cast<int>(std::floor(phi / (2 * TMath::Pi() / N_PHI_SECTORS)));
    phiBin = std::max(0, std::min(phiBin, N_PHI_SECTORS - 1));

    int cellId = -1;

    if (inFT0A) {
      // FT0-A: cells 0-95
      float etaWidth = (FT0A_ETA_MAX - FT0A_ETA_MIN) / N_ETA_A;
      int etaBin = static_cast<int>(std::floor((eta - FT0A_ETA_MIN) / etaWidth));
      etaBin = std::max(0, std::min(etaBin, N_ETA_A - 1));
      cellId = etaBin * N_PHI_SECTORS + phiBin;
    } else if (inFT0C) {
      // FT0-C: cells 96-207
      float etaWidth = (FT0C_ETA_MAX - FT0C_ETA_MIN) / N_ETA_C;
      int etaBin = static_cast<int>(std::floor((eta - FT0C_ETA_MIN) / etaWidth));
      etaBin = std::max(0, std::min(etaBin, N_ETA_C - 1));
      cellId = N_CH_A + etaBin * N_PHI_SECTORS + phiBin;
    }

    return cellId;
  }

  // ============================================
  // Track selection (Paola/Jesus criteria)
  // ============================================
  template <typename T>
  bool isSelectedTrack(const T& track)
  {
    // pT selection
    if (track.pt() < cfgPtMin)
      return false;

    // Eta selection
    if (std::abs(track.eta()) > cfgEtaMax)
      return false;

    // TPC crossed rows
    if (track.tpcNClsCrossedRows() < cfgNCrossedRowsTPC)
      return false;

    // TPC chi2 per cluster
    if (track.tpcChi2NCl() > cfgChi2PerClusterTPC)
      return false;

    // ITS chi2 per cluster
    if (track.itsChi2NCl() > cfgChi2PerClusterITS)
      return false;

    // DCA z
    if (std::abs(track.dcaZ()) > cfgDCAZ)
      return false;

    // Golden chi2 (global track)
    if (cfgRequireGoldenChi2 && !track.isGlobalTrack())
      return false;

    return true;
  }

  // ============================================
  // Process MC collisions
  // ============================================
  void processMC(
    aod::McCollision const& /* mcCollision */,
    aod::McParticles const& mcParticles)
  {
    // Initialize counters
    std::array<float, N_CELL> truthCounts;
    truthCounts.fill(0.0);

    int nch_INEL = 0;
    int nch_FT0 = 0;
    bool hasFT0A = false;
    bool hasFT0C = false;

    // Loop over MC particles
    for (const auto& particle : mcParticles) {
      // Check if primary
      if (!(particle.flags() & PHYSICAL_PRIMARY_BIT))
        continue;

      // Check if charged
      int charge = getCharge(particle.pdgCode());
      if (charge == 0)
        continue;

      // pT > 0.1
      if (particle.pt() < cfgPtMin)
        continue;

      // INEL>0: |eta| < 1
      if (std::abs(particle.eta()) < 1.0) {
        nch_INEL++;
      }

      // dNch/deta: |eta| < 0.8
      if (std::abs(particle.eta()) < cfgEtaMax) {
        nch_FT0++;
      }

      // FT0 acceptance
      bool inFT0A = (particle.eta() > FT0A_ETA_MIN && particle.eta() < FT0A_ETA_MAX);
      bool inFT0C = (particle.eta() > FT0C_ETA_MIN && particle.eta() < FT0C_ETA_MAX);

      if (inFT0A)
        hasFT0A = true;
      if (inFT0C)
        hasFT0C = true;

      // Assign to FT0 cell
      bool isFT0A = false;
      int cellId = assignToFT0Cell(particle.eta(), particle.phi(), isFT0A);

      if (cellId >= 0 && cellId < N_CELL) {
        truthCounts[cellId] += 1.0;
        histos.fill(HIST("hCellOccupancy"), cellId);
        if (isFT0A) {
          histos.fill(HIST("hCellOccupancyFT0A"), cellId);
        } else {
          histos.fill(HIST("hCellOccupancyFT0C"), cellId - N_CH_A);
        }
      }
    }

    // Event selection
    bool isINEL = (nch_INEL > 0);
    bool isFT0 = (hasFT0A && hasFT0C);

    histos.fill(HIST("hEvents"), 0); // All events

    if (isINEL) {
      histos.fill(HIST("hEvents"), 1); // INEL>0
      histos.fill(HIST("hNch_INEL"), nch_FT0);
    }

    if (isINEL && isFT0) {
      histos.fill(HIST("hEvents"), 2); // INEL>0 & FT0
      histos.fill(HIST("hNch_FT0"), nch_FT0);

      // Compute flattenicity
      float rho = computeFlattenicity(truthCounts);
      if (rho > 0) {
        float flattenicity = 1.0 - rho;
        histos.fill(HIST("hFlattenicity"), flattenicity);
        histos.fill(HIST("hFlattenicity_vs_Nch"), nch_FT0, flattenicity);

        // Multiplicity classes (based on Nch)
        if (nch_FT0 < 5) {
          histos.fill(HIST("hFlattenicity_0_10"), flattenicity);
        } else if (nch_FT0 < 8) {
          histos.fill(HIST("hFlattenicity_10_20"), flattenicity);
        } else if (nch_FT0 < 11) {
          histos.fill(HIST("hFlattenicity_20_30"), flattenicity);
        } else if (nch_FT0 < 14) {
          histos.fill(HIST("hFlattenicity_30_40"), flattenicity);
        } else if (nch_FT0 < 17) {
          histos.fill(HIST("hFlattenicity_40_50"), flattenicity);
        } else if (nch_FT0 < 20) {
          histos.fill(HIST("hFlattenicity_50_60"), flattenicity);
        } else if (nch_FT0 < 24) {
          histos.fill(HIST("hFlattenicity_60_70"), flattenicity);
        } else if (nch_FT0 < 28) {
          histos.fill(HIST("hFlattenicity_70_80"), flattenicity);
        } else if (nch_FT0 < 33) {
          histos.fill(HIST("hFlattenicity_80_90"), flattenicity);
        } else {
          histos.fill(HIST("hFlattenicity_90_100"), flattenicity);
        }
      }
    }
  }

  PROCESS_SWITCH(FlattenicityTask, processMC, "Process MC events", true);

  // ============================================
  // Process data collisions
  // ============================================
  void processData(
    aod::Collision const& collision,
    aod::FT0s const& ft0s,
    FullTracks const& tracks)
  {
    // Event selection: |vz| < 10 cm
    if (std::abs(collision.posZ()) > cfgVzMax)
      return;

    // Find FT0 matching this collision's BC
    auto ft0 = ft0s.begin();
    bool foundFT0 = false;
    for (const auto& f : ft0s) {
      if (f.bcId() == collision.bcId()) {
        ft0 = f;
        foundFT0 = true;
        break;
      }
    }
    if (!foundFT0)
      return;

    // Track selection and counting
    // int nTracks = 0;
    std::array<float, N_CELL> recoCounts;
    recoCounts.fill(0.0);

    for (const auto& track : tracks) {
      if (!isSelectedTrack(track))
        continue;
      // nTracks++;

      // Assign to FT0 cell using track extrapolation
      bool isFT0A = false;
      int cellId = assignToFT0Cell(track.eta(), track.phi(), isFT0A);
      if (cellId >= 0 && cellId < N_CELL) {
        recoCounts[cellId] += 1.0;
      }
    }

    histos.fill(HIST("hEvents"), 3); // Data events

    // Compute flattenicity from FT0 amplitudes
    // The FT0 channels are stored as arrays of (channel, amplitude) pairs
    // The channelA() and channelC() return vectors of pairs
    std::array<float, N_CELL> ft0Counts;
    ft0Counts.fill(0.0);

    // FT0-A channels (0-95)
    for (int i = 0; i < N_CH_A; i++) {
      ft0Counts[i] = ft0.channelA()[i];
    }

    // FT0-C channels (96-207)
    for (int i = 0; i < N_CH_C; i++) {
      ft0Counts[N_CH_A + i] = ft0.channelC()[i];
    }

    float rho = computeFlattenicity(ft0Counts);
    if (rho > 0) {
      histos.fill(HIST("hFlattenicity"), 1.0 - rho);
    }
  }

  PROCESS_SWITCH(FlattenicityTask, processData, "Process data events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<FlattenicityTask>(cfgc));
  return workflow;
}
