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
/// \brief Flattenicity analysis task for UE studies, with two flattenicity
///        definitions: particle-based (charged particles/tracks mapped into
///        FT0 cells) and FT0-detector-amplitude-based (real FT0 channel
///        signals), for Hyperloop
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

      // Definition 1: flattenicity from charged particles/tracks mapped into FT0 cells
      {"hFlattenicityParticles", "Flattenicity from charged particles in FT0 acceptance;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityParticles_vs_Nch", "Flattenicity (particles) vs Nch;N_{ch};1-#rho", {HistType::kTH2F, {{50, -0.5, 99.5}, {50, 0.0, 1.0}}}},

      // Definition 2: flattenicity from FT0 detector amplitudes only
      {"hFlattenicityFT0", "Flattenicity from FT0 detector amplitudes;1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},

      // FT0 cell occupancy (from particle-level mapping)
      {"hCellOccupancy", "FT0 cell occupancy;Cell ID;Entries", {HistType::kTH1F, {{NCell, 0, NCell}}}},
      {"hCellOccupancyFT0A", "FT0-A cell occupancy;Cell ID;Entries", {HistType::kTH1F, {{NchA, 0, NchA}}}},
      {"hCellOccupancyFT0C", "FT0-C cell occupancy;Cell ID;Entries", {HistType::kTH1F, {{NchC, 0, NchC}}}},

      // Percentile classes (I-VIII: 0-1,1-5,5-10,10-20,20-30,30-40,40-50,50-100%)
      // -- particle-based definition
      {"hFlattenicityParticles_0_1", "Flattenicity (particles) class I (0-1%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityParticles_1_5", "Flattenicity (particles) class II (1-5%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityParticles_5_10", "Flattenicity (particles) class III (5-10%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityParticles_10_20", "Flattenicity (particles) class IV (10-20%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityParticles_20_30", "Flattenicity (particles) class V (20-30%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityParticles_30_40", "Flattenicity (particles) class VI (30-40%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityParticles_40_50", "Flattenicity (particles) class VII (40-50%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityParticles_50_100", "Flattenicity (particles) class VIII (50-100%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},

      // -- FT0-amplitude-based definition
      {"hFlattenicityFT0_0_1", "Flattenicity (FT0) class I (0-1%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityFT0_1_5", "Flattenicity (FT0) class II (1-5%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityFT0_5_10", "Flattenicity (FT0) class III (5-10%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityFT0_10_20", "Flattenicity (FT0) class IV (10-20%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityFT0_20_30", "Flattenicity (FT0) class V (20-30%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityFT0_30_40", "Flattenicity (FT0) class VI (30-40%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityFT0_40_50", "Flattenicity (FT0) class VII (40-50%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
      {"hFlattenicityFT0_50_100", "Flattenicity (FT0) class VIII (50-100%);1-#rho;Entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
    }};

  // ============================================
  // Configurables
  // ============================================
  Configurable<float> cfgPtMin{"cfgPtMin", 0.0, "Minimum pT for tracks/particles (GeV/c)"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8, "Maximum |eta| for tracks"};
  Configurable<float> cfgVzMax{"cfgVzMax", 10.0, "Maximum |vz| for collisions"};
  Configurable<int> cfgNCrossedRowsTPC{"cfgNCrossedRowsTPC", 70, "Minimum TPC crossed rows"};
  Configurable<float> cfgChi2PerClusterTPC{"cfgChi2PerClusterTPC", 4.0, "Maximum TPC chi2 per cluster"};
  Configurable<float> cfgChi2PerClusterITS{"cfgChi2PerClusterITS", 36.0, "Maximum ITS chi2 per cluster"};
  Configurable<float> cfgDCAZ{"cfgDCAZ", 0.1, "Maximum DCA z"};
  Configurable<bool> cfgRequireGoldenChi2{"cfgRequireGoldenChi2", true, "Require golden chi2"};

  // Percentile-class boundaries (in 1-rho), one full set of 7 edges per
  // definition (classes I-VIII: 0-1,1-5,5-10,10-20,20-30,30-40,40-50,50-100%).
  // Kept as individually named Configurables for Hyperloop -- these are NOT
  // computed on the fly; they must be set from an actual percentile
  // calibration of the corresponding 1-rho distribution (particle-based and
  // FT0-amplitude-based distributions are different quantities and will not
  // share the same boundary values). The defaults below are placeholders.

  // -- particle-based definition boundaries
  Configurable<float> cfgParticlesCut0To1{"cfgParticlesCut0To1", 0.90, "Particles 1-rho lower edge for 0-1%"};
  Configurable<float> cfgParticlesCut1To5{"cfgParticlesCut1To5", 0.80, "Particles 1-rho lower edge for 1-5%"};
  Configurable<float> cfgParticlesCut5To10{"cfgParticlesCut5To10", 0.75, "Particles 1-rho lower edge for 5-10%"};
  Configurable<float> cfgParticlesCut10To20{"cfgParticlesCut10To20", 0.65, "Particles 1-rho lower edge for 10-20%"};
  Configurable<float> cfgParticlesCut20To30{"cfgParticlesCut20To30", 0.55, "Particles 1-rho lower edge for 20-30%"};
  Configurable<float> cfgParticlesCut30To40{"cfgParticlesCut30To40", 0.45, "Particles 1-rho lower edge for 30-40%"};
  Configurable<float> cfgParticlesCut40To50{"cfgParticlesCut40To50", 0.35, "Particles 1-rho lower edge for 40-50%"};

  // -- FT0-amplitude-based definition boundaries
  Configurable<float> cfgFT0Cut0To1{"cfgFT0Cut0To1", 0.904, "FT0 1-rho lower edge for 0-1%"};
  Configurable<float> cfgFT0Cut1To5{"cfgFT0Cut1To5", 0.888, "FT0 1-rho lower edge for 1-5%"};
  Configurable<float> cfgFT0Cut5To10{"cfgFT0Cut5To10", 0.840, "FT0 1-rho lower edge for 5-10%"};
  Configurable<float> cfgFT0Cut10To20{"cfgFT0Cut10To20", 0.780, "FT0 1-rho lower edge for 10-20%"};
  Configurable<float> cfgFT0Cut20To30{"cfgFT0Cut20To30", 0.720, "FT0 1-rho lower edge for 20-30%"};
  Configurable<float> cfgFT0Cut30To40{"cfgFT0Cut30To40", 0.660, "FT0 1-rho lower edge for 30-40%"};
  Configurable<float> cfgFT0Cut40To50{"cfgFT0Cut40To50", 0.600, "FT0 1-rho lower edge for 40-50%"};

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
  // Generic ALICE flattenicity formula, shared by both definitions
  // ============================================
  float computeRhoFromCells(const std::array<float, NCell>& counts)
  {
    float total = 0.0f;
    for (int i = 0; i < NCell; ++i) {
      total += counts[i];
    }
    if (total <= 0.0f) {
      return -1.0f;
    }

    float mean = total / NCell;
    if (mean <= 0.0f) {
      return -1.0f;
    }

    float sumSq = 0.0f;
    for (int i = 0; i < NCell; ++i) {
      sumSq += (counts[i] - mean) * (counts[i] - mean);
    }

    return std::sqrt(sumSq) / (NCell * mean);
  }

  // Definition 1: charged-particle/track flattenicity in FT0 acceptance.
  float computeFlattenicityParticles(const std::array<float, NCell>& particleCounts)
  {
    float rho = computeRhoFromCells(particleCounts);
    return (rho >= 0.0f) ? (1.0f - rho) : -1.0f;
  }

  // Definition 2: FT0-detector-only flattenicity from amplitudes.
  float computeFlattenicityFT0(const std::array<float, NCell>& ft0Counts)
  {
    float rho = computeRhoFromCells(ft0Counts);
    return (rho >= 0.0f) ? (1.0f - rho) : -1.0f;
  }

  // ============================================
  // Assign particle/track to FT0 cell
  // ============================================
  int assignToFT0Cell(float eta, float phi, bool& isFT0A)
  {
    bool inFT0A = (eta > FT0AEtaMin && eta < FT0AEtaMax);
    bool inFT0C = (eta > FT0CEtaMin && eta < FT0CEtaMax);
    if (!inFT0A && !inFT0C) {
      return -1;
    }

    isFT0A = inFT0A;

    int phiBin = static_cast<int>(std::floor(phi / (o2::constants::math::TwoPI / NPhiSectors)));
    phiBin = std::max(0, std::min(phiBin, NPhiSectors - 1));

    if (inFT0A) {
      // FT0-A: cells 0-95
      float etaWidth = (FT0AEtaMax - FT0AEtaMin) / NEtaA;
      int etaBin = static_cast<int>(std::floor((eta - FT0AEtaMin) / etaWidth));
      etaBin = std::max(0, std::min(etaBin, NEtaA - 1));
      return etaBin * NPhiSectors + phiBin;
    }

    // FT0-C: cells 96-207
    float etaWidth = (FT0CEtaMax - FT0CEtaMin) / NEtaC;
    int etaBin = static_cast<int>(std::floor((eta - FT0CEtaMin) / etaWidth));
    etaBin = std::max(0, std::min(etaBin, NEtaC - 1));
    return NchA + etaBin * NPhiSectors + phiBin;
  }

  // ============================================
  // Percentile-class filling helpers
  // ============================================
  void fillParticlesPercentileHistograms(float flat)
  {
    if (flat >= cfgParticlesCut0To1) {
      histos.fill(HIST("hFlattenicityParticles_0_1"), flat);
    } else if (flat >= cfgParticlesCut1To5) {
      histos.fill(HIST("hFlattenicityParticles_1_5"), flat);
    } else if (flat >= cfgParticlesCut5To10) {
      histos.fill(HIST("hFlattenicityParticles_5_10"), flat);
    } else if (flat >= cfgParticlesCut10To20) {
      histos.fill(HIST("hFlattenicityParticles_10_20"), flat);
    } else if (flat >= cfgParticlesCut20To30) {
      histos.fill(HIST("hFlattenicityParticles_20_30"), flat);
    } else if (flat >= cfgParticlesCut30To40) {
      histos.fill(HIST("hFlattenicityParticles_30_40"), flat);
    } else if (flat >= cfgParticlesCut40To50) {
      histos.fill(HIST("hFlattenicityParticles_40_50"), flat);
    } else {
      histos.fill(HIST("hFlattenicityParticles_50_100"), flat);
    }
  }

  void fillFT0PercentileHistograms(float flat)
  {
    if (flat >= cfgFT0Cut0To1) {
      histos.fill(HIST("hFlattenicityFT0_0_1"), flat);
    } else if (flat >= cfgFT0Cut1To5) {
      histos.fill(HIST("hFlattenicityFT0_1_5"), flat);
    } else if (flat >= cfgFT0Cut5To10) {
      histos.fill(HIST("hFlattenicityFT0_5_10"), flat);
    } else if (flat >= cfgFT0Cut10To20) {
      histos.fill(HIST("hFlattenicityFT0_10_20"), flat);
    } else if (flat >= cfgFT0Cut20To30) {
      histos.fill(HIST("hFlattenicityFT0_20_30"), flat);
    } else if (flat >= cfgFT0Cut30To40) {
      histos.fill(HIST("hFlattenicityFT0_30_40"), flat);
    } else if (flat >= cfgFT0Cut40To50) {
      histos.fill(HIST("hFlattenicityFT0_40_50"), flat);
    } else {
      histos.fill(HIST("hFlattenicityFT0_50_100"), flat);
    }
  }

  // ============================================
  // Track selection (Paola/Jesus criteria)
  // ============================================
  template <typename T>
  bool isSelectedTrack(const T& track)
  {
    // pT > 0
    if (track.pt() <= 0.0f || track.pt() < cfgPtMin) {
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
  // (only the particle-based definition is available here: there is no
  // FT0 detector object at MC-truth level in this process signature)
  // ============================================
  void processMC(
    aod::McCollision const& /* mcCollision */,
    aod::McParticles const& mcParticles)
  {
    std::array<float, NCell> particleCounts{};

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

      // pT > 0
      if (particle.pt() <= 0.0f || particle.pt() < cfgPtMin) {
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
        particleCounts[cellId] += 1.0f;
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

      float flatParticles = computeFlattenicityParticles(particleCounts);
      if (flatParticles >= 0.0f) {
        histos.fill(HIST("hFlattenicityParticles"), flatParticles);
        histos.fill(HIST("hFlattenicityParticles_vs_Nch"), nchFT0, flatParticles);
        fillParticlesPercentileHistograms(flatParticles);
      }
    }
  }

  PROCESS_SWITCH(FlattenicityTask, processMC, "Process MC events", true);

  // ============================================
  // Process data collisions
  // (both definitions are available here: particle-based from selected
  // reconstructed tracks, and FT0-amplitude-based from the real FT0 signal)
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

    // Definition 1: selected charged tracks mapped into FT0 cells
    std::array<float, NCell> recoParticleCounts{};
    for (const auto& track : tracks) {
      if (!isSelectedTrack(track)) {
        continue;
      }
      bool isFT0A = false;
      int cellId = assignToFT0Cell(track.eta(), track.phi(), isFT0A);
      if (cellId >= 0 && cellId < NCell) {
        recoParticleCounts[cellId] += 1.0f;
      }
    }

    // Definition 2: FT0 detector amplitudes only
    std::array<float, NCell> ft0Counts{};
    for (int i = 0; i < NchA; ++i) {
      ft0Counts[i] = ft0.channelA()[i];
    }
    for (int i = 0; i < NchC; ++i) {
      ft0Counts[NchA + i] = ft0.channelC()[i];
    }

    histos.fill(HIST("hEvents"), 3); // Data events

    float flatParticles = computeFlattenicityParticles(recoParticleCounts);
    if (flatParticles >= 0.0f) {
      histos.fill(HIST("hFlattenicityParticles"), flatParticles);
      fillParticlesPercentileHistograms(flatParticles);
    }

    float flatFT0 = computeFlattenicityFT0(ft0Counts);
    if (flatFT0 >= 0.0f) {
      histos.fill(HIST("hFlattenicityFT0"), flatFT0);
      fillFT0PercentileHistograms(flatFT0);
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
