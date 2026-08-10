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

#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
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
  static constexpr int NPhiSectors = 8;
  static constexpr int NEtaA = 12;
  static constexpr int NEtaC = 14;
  static constexpr int NchA = 96;
  static constexpr int NchC = 112;
  static constexpr int NCell = NchA + NchC; // 208

  static constexpr float FT0AEtaMin = 3.5;
  static constexpr float FT0AEtaMax = 4.9;
  static constexpr float FT0CEtaMin = -3.3;
  static constexpr float FT0CEtaMax = -2.1;

  static constexpr int NPhysicalPrimaryBit = 0x4;

  // ============================================
  // Multiplicity class boundaries
  // ============================================
  static constexpr int NchBound10 = 5;
  static constexpr int NchBound20 = 8;
  static constexpr int NchBound30 = 11;
  static constexpr int NchBound40 = 14;
  static constexpr int NchBound50 = 17;
  static constexpr int NchBound60 = 20;
  static constexpr int NchBound70 = 24;
  static constexpr int NchBound80 = 28;
  static constexpr int NchBound90 = 33;

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
      {"hCellOccupancy", "FT0 cell occupancy;Cell ID;Entries", {HistType::kTH1F, {{NCell, 0, NCell}}}},
      {"hCellOccupancyFT0A", "FT0-A cell occupancy;Cell ID;Entries", {HistType::kTH1F, {{NchA, 0, NchA}}}},
      {"hCellOccupancyFT0C", "FT0-C cell occupancy;Cell ID;Entries", {HistType::kTH1F, {{NchC, 0, NchC}}}},

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
      case PDG_t::kPiPlus: // 211
      case PDG_t::kKPlus:  // 321
      case PDG_t::kProton: // 2212
        return 1;
      case PDG_t::kElectron: // 11
      case PDG_t::kMuonPlus: // 13
        return -1;
      default:
        return 0;
    }
  }

  // ============================================
  // Flattenicity calculation
  // ============================================
  float computeFlattenicity(const std::array<float, NCell>& counts)
  {
    float total = 0.0;
    for (int i = 0; i < NCell; i++) {
      total += counts[i];
    }

    if (total <= 0) {
      return -1.0;
    }

    float mean = total / NCell;
    if (mean <= 0) {
      return -1.0;
    }

    float sumSq = 0.0;
    for (int i = 0; i < NCell; i++) {
      sumSq += (counts[i] - mean) * (counts[i] - mean);
    }

    float rho = std::sqrt(sumSq) / (NCell * mean);
    return rho;
  }

  // ============================================
  // Assign particle to FT0 cell
  // ============================================
  int assignToFT0Cell(float eta, float phi, bool& isFT0A)
  {
    // Check if in FT0 acceptance
    bool inFT0A = (eta > FT0AEtaMin && eta < FT0AEtaMax);
    bool inFT0C = (eta > FT0CEtaMin && eta < FT0CEtaMax);

    if (!inFT0A && !inFT0C) {
      return -1;
    }

    isFT0A = inFT0A;

    // Phi bin
    int phiBin = static_cast<int>(std::floor(phi / (o2::constants::math::TwoPI / NPhiSectors)));
    phiBin = std::max(0, std::min(phiBin, NPhiSectors - 1));

    int cellId = -1;

    if (inFT0A) {
      // FT0-A: cells 0-95
      float etaWidth = (FT0AEtaMax - FT0AEtaMin) / NEtaA;
      int etaBin = static_cast<int>(std::floor((eta - FT0AEtaMin) / etaWidth));
      etaBin = std::max(0, std::min(etaBin, NEtaA - 1));
      cellId = etaBin * NPhiSectors + phiBin;
    } else if (inFT0C) {
      // FT0-C: cells 96-207
      float etaWidth = (FT0CEtaMax - FT0CEtaMin) / NEtaC;
      int etaBin = static_cast<int>(std::floor((eta - FT0CEtaMin) / etaWidth));
      etaBin = std::max(0, std::min(etaBin, NEtaC - 1));
      cellId = NchA + etaBin * NPhiSectors + phiBin;
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
    if (track.pt() < cfgPtMin) {
      return false;
    }

    // Eta selection
    if (std::abs(track.eta()) > cfgEtaMax) {
      return false;
    }

    // TPC crossed rows
    if (track.tpcNClsCrossedRows() < cfgNCrossedRowsTPC) {
      return false;
    }

    // TPC chi2 per cluster
    if (track.tpcChi2NCl() > cfgChi2PerClusterTPC) {
      return false;
    }

    // ITS chi2 per cluster
    if (track.itsChi2NCl() > cfgChi2PerClusterITS) {
      return false;
    }

    // DCA z
    if (std::abs(track.dcaZ()) > cfgDCAZ) {
      return false;
    }

    // Golden chi2 (global track)
    if (cfgRequireGoldenChi2 && !track.isGlobalTrack()) {
      return false;
    }

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
    std::array<float, NCell> truthCounts{};

    int nchINEL = 0;
    int nchFT0 = 0;
    bool hasFT0A = false;
    bool hasFT0C = false;

    // Loop over MC particles
    for (const auto& particle : mcParticles) {
      // Check if primary
      if ((particle.flags() & NPhysicalPrimaryBit) == 0) {
        continue;
      }

      // Check if charged
      int charge = getCharge(particle.pdgCode());
      if (charge == 0) {
        continue;
      }

      // pT > 0.1
      if (particle.pt() < cfgPtMin) {
        continue;
      }

      // INEL>0: |eta| < 1
      if (std::abs(particle.eta()) < 1.0) {
        nchINEL++;
      }

      // dNch/deta: |eta| < 0.8
      if (std::abs(particle.eta()) < cfgEtaMax) {
        nchFT0++;
      }

      // FT0 acceptance
      bool inFT0A = (particle.eta() > FT0AEtaMin && particle.eta() < FT0AEtaMax);
      bool inFT0C = (particle.eta() > FT0CEtaMin && particle.eta() < FT0CEtaMax);

      if (inFT0A) {
        hasFT0A = true;
      }
      if (inFT0C) {
        hasFT0C = true;
      }

      // Assign to FT0 cell
      bool isFT0A = false;
      int cellId = assignToFT0Cell(particle.eta(), particle.phi(), isFT0A);

      if (cellId >= 0 && cellId < NCell) {
        truthCounts[cellId] += 1.0;
        histos.fill(HIST("hCellOccupancy"), cellId);
        if (isFT0A) {
          histos.fill(HIST("hCellOccupancyFT0A"), cellId);
        } else {
          histos.fill(HIST("hCellOccupancyFT0C"), cellId - NchA);
        }
      }
    }

    // Event selection
    bool isINEL = (nchINEL > 0);
    bool isFT0 = (hasFT0A && hasFT0C);

    histos.fill(HIST("hEvents"), 0); // All events

    if (isINEL) {
      histos.fill(HIST("hEvents"), 1); // INEL>0
      histos.fill(HIST("hNch_INEL"), nchFT0);
    }

    if (isINEL && isFT0) {
      histos.fill(HIST("hEvents"), 2); // INEL>0 & FT0
      histos.fill(HIST("hNch_FT0"), nchFT0);

      // Compute flattenicity
      float rho = computeFlattenicity(truthCounts);
      if (rho > 0) {
        float flattenicity = 1.0 - rho;
        histos.fill(HIST("hFlattenicity"), flattenicity);
        histos.fill(HIST("hFlattenicity_vs_Nch"), nchFT0, flattenicity);

        // Multiplicity classes (based on Nch)
        if (nchFT0 < NchBound10) {
          histos.fill(HIST("hFlattenicity_0_10"), flattenicity);
        } else if (nchFT0 < NchBound20) {
          histos.fill(HIST("hFlattenicity_10_20"), flattenicity);
        } else if (nchFT0 < NchBound30) {
          histos.fill(HIST("hFlattenicity_20_30"), flattenicity);
        } else if (nchFT0 < NchBound40) {
          histos.fill(HIST("hFlattenicity_30_40"), flattenicity);
        } else if (nchFT0 < NchBound50) {
          histos.fill(HIST("hFlattenicity_40_50"), flattenicity);
        } else if (nchFT0 < NchBound60) {
          histos.fill(HIST("hFlattenicity_50_60"), flattenicity);
        } else if (nchFT0 < NchBound70) {
          histos.fill(HIST("hFlattenicity_60_70"), flattenicity);
        } else if (nchFT0 < NchBound80) {
          histos.fill(HIST("hFlattenicity_70_80"), flattenicity);
        } else if (nchFT0 < NchBound90) {
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
    if (std::abs(collision.posZ()) > cfgVzMax) {
      return;
    }

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
    if (!foundFT0) {
      return;
    }

    // Track selection and counting
    std::array<float, NCell> recoCounts{};

    for (const auto& track : tracks) {
      if (!isSelectedTrack(track)) {
        continue;
      }

      // Assign to FT0 cell using track extrapolation
      bool isFT0A = false;
      int cellId = assignToFT0Cell(track.eta(), track.phi(), isFT0A);
      if (cellId >= 0 && cellId < NCell) {
        recoCounts[cellId] += 1.0;
      }
    }

    histos.fill(HIST("hEvents"), 3); // Data events

    // Compute flattenicity from FT0 amplitudes
    std::array<float, NCell> ft0Counts{};

    // FT0-A channels (0-95)
    for (int i = 0; i < NchA; i++) {
      ft0Counts[i] = ft0.channelA()[i];
    }

    // FT0-C channels (96-207)
    for (int i = 0; i < NchC; i++) {
      ft0Counts[NchA + i] = ft0.channelC()[i];
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
