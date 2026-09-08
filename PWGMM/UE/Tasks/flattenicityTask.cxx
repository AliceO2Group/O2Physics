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

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
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
using CollisionsWithCent = soa::Join<aod::Collisions, aod::CentFT0Ms>;
using CollisionsWithCentAndMcLabel = soa::Join<aod::Collisions, aod::CentFT0Ms, aod::McCollisionLabels>;

struct FlattenicityTask {
  // ============================================
  // FT0 Constants (Based on multFilter.cxx from Antonio)
  // ============================================
  static constexpr int NPhiSectors = 8;
  static constexpr int NSectorsA = 24;                       // FT0-A sectors
  static constexpr int NSectorsC = 28;                       // FT0-C sectors
  static constexpr int ChannelsPerSector = 4;                // Channels per sector
  static constexpr int NchA = NSectorsA * ChannelsPerSector; // 96
  static constexpr int NchC = NSectorsC * ChannelsPerSector; // 112
  static constexpr int NCell = NchA + NchC;                  // 208

  // Eta-bin counts for the particle/track cell grid (assignToFT0Cell).
  // Chosen so that NEtaA*NPhiSectors == NchA and NEtaC*NPhiSectors == NchC
  // exactly, so FT0-A and FT0-C cell IDs cannot overlap. This is separate
  // from NSectorsA/NSectorsC above, which are only used for the FT0
  // detector-amplitude flattenicity (hFlattenicityFT0/FT0A/FT0C).
  static constexpr int NEtaA = NchA / NPhiSectors; // 12
  static constexpr int NEtaC = NchC / NPhiSectors; // 14

  static constexpr float FT0AEtaMin = 3.5;
  static constexpr float FT0AEtaMax = 4.9;
  static constexpr float FT0CEtaMin = -3.3;
  static constexpr float FT0CEtaMax = -2.1;

  static constexpr int NPhysicalPrimaryBit = 0x4;

  // Centrality percentile class boundaries
  static constexpr float CentBound1 = 1.0f;
  static constexpr float CentBound5 = 5.0f;
  static constexpr float CentBound10 = 10.0f;
  static constexpr float CentBound20 = 20.0f;
  static constexpr float CentBound30 = 30.0f;
  static constexpr float CentBound40 = 40.0f;
  static constexpr float CentBound50 = 50.0f;
  static constexpr float CentBound60 = 60.0f;
  static constexpr float CentBound70 = 70.0f;
  static constexpr float CentBound80 = 80.0f;
  static constexpr float CentBound90 = 90.0f;
  static constexpr float CentBound95 = 95.0f;

  // ============================================
  // Helper functions to map FT0 channels to sectors (from multFilter.cxx)
  // ============================================
  int getT0ASector(int i_ch)
  {
    for (int i_sec = 0; i_sec < NSectorsA; ++i_sec) {
      if (i_ch >= ChannelsPerSector * i_sec && i_ch <= (ChannelsPerSector - 1) + ChannelsPerSector * i_sec) {
        return i_sec;
      }
    }
    return -1;
  }

  int getT0CSector(int i_ch)
  {
    for (int i_sec = 0; i_sec < NSectorsC; ++i_sec) {
      if (i_ch >= ChannelsPerSector * i_sec && i_ch <= (ChannelsPerSector - 1) + ChannelsPerSector * i_sec) {
        return i_sec;
      }
    }
    return -1;
  }

  // ============================================
  // Histogram Definitions - 100 BINS (bin width = 0.01)
  // ============================================
  HistogramRegistry histos{
    "histos",
    {
      // Event-level
      {"hEvents", "Event selection;;Counts", {HistType::kTH1F, {{5, 0, 5}}}},

      // dNch/deta
      {"hNch_INEL", "Nch distribution (INEL>0);N_{ch};Entries", {HistType::kTH1F, {{100, -0.5, 99.5}}}},
      {"hNch_FT0", "Nch distribution (INEL>0 & FT0);N_{ch};Entries", {HistType::kTH1F, {{100, -0.5, 99.5}}}},

      // ============================================
      // Flattenicity Histograms - 100 BINS!
      // Bin width = 0.01, allowing 0-1% and 1-5% classes
      // ============================================
      {"hFlattenicityParticles", "Flattenicity from charged particles;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityParticles_vs_Nch", "Flattenicity (particles) vs Nch;N_{ch};1-#rho", {HistType::kTH2F, {{50, -0.5, 99.5}, {100, 0.0, 1.0}}}},
      {"hFlattenicityFT0", "Flattenicity from FT0 detector amplitudes (avg of FT0-A and FT0-C);1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0A", "Flattenicity from FT0-A amplitudes (24 sectors);1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0C", "Flattenicity from FT0-C amplitudes (28 sectors);1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},

      // FT0 cell occupancy
      {"hCellOccupancy", "FT0 cell occupancy;Cell ID;Entries", {HistType::kTH1F, {{NCell, 0, NCell}}}},
      {"hCellOccupancyFT0A", "FT0-A cell occupancy;Cell ID;Entries", {HistType::kTH1F, {{NchA, 0, NchA}}}},
      {"hCellOccupancyFT0C", "FT0-C cell occupancy;Cell ID;Entries", {HistType::kTH1F, {{NchC, 0, NchC}}}},

      // ============================================
      // MULTIPLICITY CLASSES USING PERCENTILES
      // All classes have 100 bins (bin width = 0.01)
      // ============================================
      // 0-1%
      {"hFlattenicityParticles_0_1", "Flattenicity (particles) class 0-1%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_0_1", "Flattenicity (FT0) class 0-1%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      // 1-5%
      {"hFlattenicityParticles_1_5", "Flattenicity (particles) class 1-5%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_1_5", "Flattenicity (FT0) class 1-5%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      // 5-10%
      {"hFlattenicityParticles_5_10", "Flattenicity (particles) class 5-10%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_5_10", "Flattenicity (FT0) class 5-10%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      // 10-20%
      {"hFlattenicityParticles_10_20", "Flattenicity (particles) class 10-20%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_10_20", "Flattenicity (FT0) class 10-20%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      // 20-30%
      {"hFlattenicityParticles_20_30", "Flattenicity (particles) class 20-30%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_20_30", "Flattenicity (FT0) class 20-30%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      // 30-40%
      {"hFlattenicityParticles_30_40", "Flattenicity (particles) class 30-40%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_30_40", "Flattenicity (FT0) class 30-40%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      // 40-50%
      {"hFlattenicityParticles_40_50", "Flattenicity (particles) class 40-50%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_40_50", "Flattenicity (FT0) class 40-50%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      // 50-60%
      {"hFlattenicityParticles_50_60", "Flattenicity (particles) class 50-60%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_50_60", "Flattenicity (FT0) class 50-60%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      // 60-70%
      {"hFlattenicityParticles_60_70", "Flattenicity (particles) class 60-70%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_60_70", "Flattenicity (FT0) class 60-70%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      // 70-80%
      {"hFlattenicityParticles_70_80", "Flattenicity (particles) class 70-80%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_70_80", "Flattenicity (FT0) class 70-80%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      // 80-90%
      {"hFlattenicityParticles_80_90", "Flattenicity (particles) class 80-90%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_80_90", "Flattenicity (FT0) class 80-90%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      // 90-95%
      {"hFlattenicityParticles_90_95", "Flattenicity (particles) class 90-95%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_90_95", "Flattenicity (FT0) class 90-95%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      // 95-100%
      {"hFlattenicityParticles_95_100", "Flattenicity (particles) class 95-100%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"hFlattenicityFT0_95_100", "Flattenicity (FT0) class 95-100%;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},

      // ============================================
      // GEN vs RECO CORRELATION (3 plots)
      // ============================================
      {"hFlatParticles_Gen_vs_Rec", "Flattenicity (particles): gen vs reco;1-#rho (gen);1-#rho (reco)", {HistType::kTH2F, {{100, 0.0, 1.0}, {100, 0.0, 1.0}}}},
      {"hFlatFT0_Gen_vs_Rec", "Flattenicity: gen (particles) vs reco (FT0 amplitudes);1-#rho (gen);1-#rho (reco FT0)", {HistType::kTH2F, {{100, 0.0, 1.0}, {100, 0.0, 1.0}}}},
      {"hFlattenicityGen_MatchedToReco", "Gen-level flattenicity (only events with a reco match);1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
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
  // Templated on array size so it can be reused for the 208-cell
  // particle-level grid as well as the 24-cell (FT0-A) and
  // 28-cell (FT0-C) detector-amplitude arrays.
  // ============================================
  template <std::size_t N>
  float computeFlattenicity(const std::array<float, N>& counts)
  {
    float total = 0.0;
    for (std::size_t i = 0; i < N; i++) {
      total += counts[i];
    }
    if (total <= 0) {
      return -1.0;
    }
    float mean = total / N;
    if (mean <= 0) {
      return -1.0;
    }
    float sumSq = 0.0;
    for (std::size_t i = 0; i < N; i++) {
      sumSq += (counts[i] - mean) * (counts[i] - mean);
    }
    float rho = std::sqrt(sumSq) / (N * mean);
    return rho;
  }

  // ============================================
  // Assign particle to FT0 cell (Updated to match multFilter.cxx)
  // ============================================
  int assignToFT0Cell(float eta, float phi, bool& isFT0A)
  {
    bool inFT0A = (eta > FT0AEtaMin && eta < FT0AEtaMax);
    bool inFT0C = (eta > FT0CEtaMin && eta < FT0CEtaMax);
    if (!inFT0A && !inFT0C) {
      return -1;
    }
    isFT0A = inFT0A;

    // Map phi to sector (0-7)
    int phiBin = static_cast<int>(std::floor(phi / (o2::constants::math::TwoPI / NPhiSectors)));
    phiBin = std::max(0, std::min(phiBin, NPhiSectors - 1));

    int cellId = -1;
    if (inFT0A) {
      // FT0-A: NEtaA eta bins x NPhiSectors phi bins = NchA cells (0..NchA-1)
      float etaWidth = (FT0AEtaMax - FT0AEtaMin) / NEtaA;
      int sector = static_cast<int>(std::floor((eta - FT0AEtaMin) / etaWidth));
      sector = std::max(0, std::min(sector, NEtaA - 1));
      cellId = sector * NPhiSectors + phiBin;
    } else if (inFT0C) {
      // FT0-C: NEtaC eta bins x NPhiSectors phi bins = NchC cells,
      // offset by NchA so FT0-A and FT0-C cell IDs never overlap.
      float etaWidth = (FT0CEtaMax - FT0CEtaMin) / NEtaC;
      int sector = static_cast<int>(std::floor((eta - FT0CEtaMin) / etaWidth));
      sector = std::max(0, std::min(sector, NEtaC - 1));
      cellId = NchA + sector * NPhiSectors + phiBin;
    }
    return cellId;
  }

  // ============================================
  // Track selection
  // ============================================
  template <typename T>
  bool isSelectedTrack(const T& track)
  {
    if (track.pt() < cfgPtMin)
      return false;
    if (std::abs(track.eta()) > cfgEtaMax)
      return false;
    if (track.tpcNClsCrossedRows() < cfgNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > cfgChi2PerClusterTPC)
      return false;
    if (track.itsChi2NCl() > cfgChi2PerClusterITS)
      return false;
    if (std::abs(track.dcaZ()) > cfgDCAZ)
      return false;
    if (cfgRequireGoldenChi2 && !track.isGlobalTrack())
      return false;
    return true;
  }

  // ============================================
  // Helper function to fill multiplicity class histograms
  // ============================================
  void fillMultiplicityClass(float flattenicity, float centrality, bool isParticle)
  {
    if (isParticle) {
      if (centrality < CentBound1) {
        histos.fill(HIST("hFlattenicityParticles_0_1"), flattenicity);
      } else if (centrality < CentBound5) {
        histos.fill(HIST("hFlattenicityParticles_1_5"), flattenicity);
      } else if (centrality < CentBound10) {
        histos.fill(HIST("hFlattenicityParticles_5_10"), flattenicity);
      } else if (centrality < CentBound20) {
        histos.fill(HIST("hFlattenicityParticles_10_20"), flattenicity);
      } else if (centrality < CentBound30) {
        histos.fill(HIST("hFlattenicityParticles_20_30"), flattenicity);
      } else if (centrality < CentBound40) {
        histos.fill(HIST("hFlattenicityParticles_30_40"), flattenicity);
      } else if (centrality < CentBound50) {
        histos.fill(HIST("hFlattenicityParticles_40_50"), flattenicity);
      } else if (centrality < CentBound60) {
        histos.fill(HIST("hFlattenicityParticles_50_60"), flattenicity);
      } else if (centrality < CentBound70) {
        histos.fill(HIST("hFlattenicityParticles_60_70"), flattenicity);
      } else if (centrality < CentBound80) {
        histos.fill(HIST("hFlattenicityParticles_70_80"), flattenicity);
      } else if (centrality < CentBound90) {
        histos.fill(HIST("hFlattenicityParticles_80_90"), flattenicity);
      } else if (centrality < CentBound95) {
        histos.fill(HIST("hFlattenicityParticles_90_95"), flattenicity);
      } else {
        histos.fill(HIST("hFlattenicityParticles_95_100"), flattenicity);
      }
    } else {
      if (centrality < CentBound1) {
        histos.fill(HIST("hFlattenicityFT0_0_1"), flattenicity);
      } else if (centrality < CentBound5) {
        histos.fill(HIST("hFlattenicityFT0_1_5"), flattenicity);
      } else if (centrality < CentBound10) {
        histos.fill(HIST("hFlattenicityFT0_5_10"), flattenicity);
      } else if (centrality < CentBound20) {
        histos.fill(HIST("hFlattenicityFT0_10_20"), flattenicity);
      } else if (centrality < CentBound30) {
        histos.fill(HIST("hFlattenicityFT0_20_30"), flattenicity);
      } else if (centrality < CentBound40) {
        histos.fill(HIST("hFlattenicityFT0_30_40"), flattenicity);
      } else if (centrality < CentBound50) {
        histos.fill(HIST("hFlattenicityFT0_40_50"), flattenicity);
      } else if (centrality < CentBound60) {
        histos.fill(HIST("hFlattenicityFT0_50_60"), flattenicity);
      } else if (centrality < CentBound70) {
        histos.fill(HIST("hFlattenicityFT0_60_70"), flattenicity);
      } else if (centrality < CentBound80) {
        histos.fill(HIST("hFlattenicityFT0_70_80"), flattenicity);
      } else if (centrality < CentBound90) {
        histos.fill(HIST("hFlattenicityFT0_80_90"), flattenicity);
      } else if (centrality < CentBound95) {
        histos.fill(HIST("hFlattenicityFT0_90_95"), flattenicity);
      } else {
        histos.fill(HIST("hFlattenicityFT0_95_100"), flattenicity);
      }
    }
  }

  // ============================================
  // Process MC collisions (generator level only)
  // ============================================
  void processMC(aod::McCollisions const& /* mcCollisions */, aod::McParticles const& mcParticles)
  {
    if (mcParticles.size() == 0) {
      LOG(warning) << "No MC particles found in this data frame";
      return;
    }
    std::array<float, NCell> truthCounts{};
    int nchINEL = 0;
    int nchFT0 = 0;
    bool hasFT0A = false;
    bool hasFT0C = false;

    for (const auto& particle : mcParticles) {
      if ((particle.flags() & NPhysicalPrimaryBit) == 0)
        continue;
      if (getCharge(particle.pdgCode()) == 0)
        continue;
      if (particle.pt() < cfgPtMin)
        continue;

      if (std::abs(particle.eta()) < 1.0)
        nchINEL++;
      if (std::abs(particle.eta()) < cfgEtaMax)
        nchFT0++;

      bool inFT0A = (particle.eta() > FT0AEtaMin && particle.eta() < FT0AEtaMax);
      bool inFT0C = (particle.eta() > FT0CEtaMin && particle.eta() < FT0CEtaMax);
      if (inFT0A)
        hasFT0A = true;
      if (inFT0C)
        hasFT0C = true;

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

    bool isINEL = (nchINEL > 0);
    bool isFT0 = (hasFT0A && hasFT0C);
    histos.fill(HIST("hEvents"), 0);
    if (isINEL) {
      histos.fill(HIST("hEvents"), 1);
      histos.fill(HIST("hNch_INEL"), nchFT0);
    }
    if (isINEL && isFT0) {
      histos.fill(HIST("hEvents"), 2);
      histos.fill(HIST("hNch_FT0"), nchFT0);
      float rho = computeFlattenicity(truthCounts);
      if (rho > 0) {
        float flattenicity = 1.0 - rho;
        histos.fill(HIST("hFlattenicityParticles"), flattenicity);
        histos.fill(HIST("hFlattenicityParticles_vs_Nch"), nchFT0, flattenicity);
      }
    }
  }
  PROCESS_SWITCH(FlattenicityTask, processMC, "Process MC events", true);

  // ============================================
  // Process data collisions (reconstruction level only)
  // ============================================
  void processData(CollisionsWithCent::iterator const& collision, aod::FT0s const& ft0s, FullTracks const& tracks)
  {
    if (std::abs(collision.posZ()) > cfgVzMax)
      return;
    float centrality = collision.centFT0M();

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

    std::array<float, NCell> recoCounts{};
    for (const auto& track : tracks) {
      if (!isSelectedTrack(track))
        continue;
      bool isFT0A = false;
      int cellId = assignToFT0Cell(track.eta(), track.phi(), isFT0A);
      if (cellId >= 0 && cellId < NCell) {
        recoCounts[cellId] += 1.0;
      }
    }
    histos.fill(HIST("hEvents"), 3);

    float rhoParticles = computeFlattenicity(recoCounts);
    if (rhoParticles > 0) {
      float flattenicity = 1.0 - rhoParticles;
      histos.fill(HIST("hFlattenicityParticles"), flattenicity);
      fillMultiplicityClass(flattenicity, centrality, true);
    }

    // ============================================
    // FT0 detector-amplitude flattenicity.
    // Matches multFilter.cxx (Antonio): FT0-A (24 sectors) and
    // FT0-C (28 sectors) are each turned into their own flattenicity
    // value, then averaged. No further phi subdivision is applied
    // here, unlike the particle-level cell grid above.
    // ============================================
    std::array<float, NSectorsA> ft0CountsA{};
    std::array<float, NSectorsC> ft0CountsC{};
    if (ft0.amplitudeA().size() > 0) {
      for (std::size_t i = 0; i < ft0.amplitudeA().size(); i++) {
        uint8_t channel = ft0.channelA()[i];
        int sector = getT0ASector(channel);
        if (sector >= 0 && sector < NSectorsA) {
          ft0CountsA[sector] += ft0.amplitudeA()[i];
        }
      }
    }
    if (ft0.amplitudeC().size() > 0) {
      for (std::size_t i = 0; i < ft0.amplitudeC().size(); i++) {
        uint8_t channel = ft0.channelC()[i];
        int sector = getT0CSector(channel);
        if (sector >= 0 && sector < NSectorsC) {
          ft0CountsC[sector] += ft0.amplitudeC()[i];
        }
      }
    }
    float rhoFT0A = computeFlattenicity(ft0CountsA);
    float rhoFT0C = computeFlattenicity(ft0CountsC);
    if (rhoFT0A > 0) {
      histos.fill(HIST("hFlattenicityFT0A"), 1.0 - rhoFT0A);
    }
    if (rhoFT0C > 0) {
      histos.fill(HIST("hFlattenicityFT0C"), 1.0 - rhoFT0C);
    }
    if (rhoFT0A > 0 && rhoFT0C > 0) {
      float rhoFT0 = (rhoFT0A + rhoFT0C) / 2.0;
      float flattenicity = 1.0 - rhoFT0;
      histos.fill(HIST("hFlattenicityFT0"), flattenicity);
      fillMultiplicityClass(flattenicity, centrality, false);
    }
  }
  PROCESS_SWITCH(FlattenicityTask, processData, "Process data events", true);

  // ============================================
  // Process GEN vs RECO correlation (3 additional plots)
  // ============================================
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;
  SliceCache cache;

  void processGenRecCorrelation(
    CollisionsWithCentAndMcLabel::iterator const& collision,
    aod::McCollisions const&,
    aod::McParticles const& mcParticles,
    aod::FT0s const& ft0s,
    FullTracks const& tracks)
  {
    if (std::abs(collision.posZ()) > cfgVzMax)
      return;
    if (!collision.has_mcCollision())
      return;

    const auto& mcCollision = collision.mcCollision();
    const auto& particlesInThisCollision = mcParticles.sliceBy(perMCCol, mcCollision.globalIndex());

    std::array<float, NCell> truthCounts{};
    for (const auto& particle : particlesInThisCollision) {
      if ((particle.flags() & NPhysicalPrimaryBit) == 0)
        continue;
      if (getCharge(particle.pdgCode()) == 0)
        continue;
      if (particle.pt() < cfgPtMin)
        continue;
      bool isFT0A = false;
      int cellId = assignToFT0Cell(particle.eta(), particle.phi(), isFT0A);
      if (cellId >= 0 && cellId < NCell) {
        truthCounts[cellId] += 1.0;
      }
    }
    float rhoGen = computeFlattenicity(truthCounts);
    if (rhoGen <= 0)
      return;
    float flatGen = 1.0 - rhoGen;

    auto ft0 = ft0s.begin();
    bool foundFT0 = false;
    for (const auto& f : ft0s) {
      if (f.bcId() == collision.bcId()) {
        ft0 = f;
        foundFT0 = true;
        break;
      }
    }

    std::array<float, NCell> recoCounts{};
    for (const auto& track : tracks) {
      if (!isSelectedTrack(track))
        continue;
      bool isFT0A = false;
      int cellId = assignToFT0Cell(track.eta(), track.phi(), isFT0A);
      if (cellId >= 0 && cellId < NCell) {
        recoCounts[cellId] += 1.0;
      }
    }
    float rhoRecParticles = computeFlattenicity(recoCounts);
    histos.fill(HIST("hFlattenicityGen_MatchedToReco"), flatGen);

    if (rhoRecParticles > 0) {
      float flatRecParticles = 1.0 - rhoRecParticles;
      histos.fill(HIST("hFlatParticles_Gen_vs_Rec"), flatGen, flatRecParticles);
    }

    if (foundFT0) {
      // Same FT0-A/FT0-C split-then-average fix as processData above.
      std::array<float, NSectorsA> ft0CountsA{};
      std::array<float, NSectorsC> ft0CountsC{};
      if (ft0.amplitudeA().size() > 0) {
        for (std::size_t i = 0; i < ft0.amplitudeA().size(); i++) {
          uint8_t channel = ft0.channelA()[i];
          int sector = getT0ASector(channel);
          if (sector >= 0 && sector < NSectorsA) {
            ft0CountsA[sector] += ft0.amplitudeA()[i];
          }
        }
      }
      if (ft0.amplitudeC().size() > 0) {
        for (std::size_t i = 0; i < ft0.amplitudeC().size(); i++) {
          uint8_t channel = ft0.channelC()[i];
          int sector = getT0CSector(channel);
          if (sector >= 0 && sector < NSectorsC) {
            ft0CountsC[sector] += ft0.amplitudeC()[i];
          }
        }
      }
      float rhoRecFT0A = computeFlattenicity(ft0CountsA);
      float rhoRecFT0C = computeFlattenicity(ft0CountsC);
      if (rhoRecFT0A > 0 && rhoRecFT0C > 0) {
        float rhoRecFT0 = (rhoRecFT0A + rhoRecFT0C) / 2.0;
        float flatRecFT0 = 1.0 - rhoRecFT0;
        histos.fill(HIST("hFlatFT0_Gen_vs_Rec"), flatGen, flatRecFT0);
      }
    }
  }
  PROCESS_SWITCH(FlattenicityTask, processGenRecCorrelation, "Process gen-vs-reco correlation (needs MC file)", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<FlattenicityTask>(cfgc));
  return workflow;
}
