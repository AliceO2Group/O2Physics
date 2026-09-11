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
/// \brief Run-1 calibration task for FT0 flattenicity (fine-binned histograms for percentile extraction)
/// \author Eisha Rani
/// \since September 2026

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using FullTracks =
  soa::Join<aod::Tracks,
            aod::TracksExtra,
            aod::TracksDCA,
            aod::TrackSelection>;

using CollisionsWithSelection =
  soa::Join<aod::Collisions, aod::EvSels>;

struct FlattenicityTask {

  // NOLINTNEXTLINE(misc-include-cleaner) - Service<> is provided by Framework/AnalysisTask.h; there is no separate standalone header for it in this O2 version.
  Service<o2::framework::O2DatabasePDG> pdg{};

  // ========================================================================
  // FT0 geometry
  // ========================================================================

  static constexpr int NPhiSectors = 8;
  static constexpr int NSectorsA = 24;
  static constexpr int NSectorsC = 28;
  static constexpr int NSectorsCombined = NSectorsA + NSectorsC; // 52
  static constexpr int ChannelsPerSector = 4;
  static constexpr int NchA = NSectorsA * ChannelsPerSector; // 96
  static constexpr int NchC = NSectorsC * ChannelsPerSector; // 112
  static constexpr int NCell = NchA + NchC;                  // 208
  static constexpr int NEtaA = NchA / NPhiSectors;           // 12
  static constexpr int NEtaC = NchC / NPhiSectors;           // 14
  static constexpr float FT0AEtaMin = 3.5;
  static constexpr float FT0AEtaMax = 4.9;
  static constexpr float FT0CEtaMin = -3.3;
  static constexpr float FT0CEtaMax = -2.1;
  static constexpr int NPhysicalPrimaryBit = 0x4;
  static constexpr float ChargeEpsilon = 1e-6f;

  // ========================================================================
  // Run-1 fine binning
  // ========================================================================
  static constexpr int NFlatBins = 10000;

  // ========================================================================
  // Histograms
  // ========================================================================
  HistogramRegistry histos{
    "histos",
    {
      {"QA/hEvents", "Event selection;selection;Entries", {HistType::kTH1F, {{5, 0, 5}}}},
      {"QA/hNchINEL", "INEL>0 multiplicity (|eta|<0.8);N_{ch};Entries", {HistType::kTH1F, {{100, -0.5, 99.5}}}},

      {"Truth/hFlattenicityTruthAverage_fine", "MC truth flattenicity - average of FT0A-like/FT0C-like grids;1-#rho;Entries", {HistType::kTH1F, {{NFlatBins, 0.0, 1.0}}}},
      {"Truth/hFlattenicityTruthCombined_fine", "MC truth flattenicity - combined 208-cell grid;1-#rho;Entries", {HistType::kTH1F, {{NFlatBins, 0.0, 1.0}}}},
      {"Reco/hFlattenicityFT0Average_fine", "FT0 flattenicity - average of FT0A and FT0C;1-#rho;Entries", {HistType::kTH1F, {{NFlatBins, 0.0, 1.0}}}},
      {"Reco/hFlattenicityFT0Combined_fine", "FT0 flattenicity - combined FT0A+FT0C;1-#rho;Entries", {HistType::kTH1F, {{NFlatBins, 0.0, 1.0}}}},
      {"Reco/hFT0Mraw_fine", "Raw FT0M amplitude sum;FT0A+FT0C amplitude;Entries", {HistType::kTH1F, {{NFlatBins, 0.0, 20000.0}}}},

      {"Truth/hFlattenicityTruthAverage", "MC truth flattenicity - average definition;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"Truth/hFlattenicityTruthCombined", "MC truth flattenicity - combined definition;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"Reco/hFlattenicityFT0Average", "FT0 flattenicity - average definition;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"Reco/hFlattenicityFT0Combined", "FT0 flattenicity - combined definition;1-#rho;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}}},
      {"Reco/hFT0Mraw", "Raw FT0M amplitude sum;FT0A+FT0C amplitude;Entries", {HistType::kTH1F, {{200, 0.0, 20000.0}}}},

      {"Truth/hCellOccupancy", "Truth FT0 cell occupancy;cell ID;Entries", {HistType::kTH1F, {{NCell, 0, NCell}}}},
      {"Truth/hCellOccupancyFT0A", "Truth FT0A cell occupancy;cell ID;Entries", {HistType::kTH1F, {{NchA, 0, NchA}}}},
      {"Truth/hCellOccupancyFT0C", "Truth FT0C cell occupancy;cell ID;Entries", {HistType::kTH1F, {{NchC, 0, NchC}}}},
    }};

  // ========================================================================
  // Configurables
  // ========================================================================

  Configurable<float> cfgPtMin{"cfgPtMin", 0.0, "Minimum pT for charged particles/tracks (GeV/c) - Antonio's instruction: require pT > 0, not a physical pT-threshold cut"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8, "Maximum |eta| for INEL>0 and track selection"};
  Configurable<float> cfgVzMax{"cfgVzMax", 10.0, "Maximum absolute collision vertex z (cm)"};
  Configurable<bool> cfgApplySel8{"cfgApplySel8", true, "Apply sel8 event selection (set to true for all samples)"};
  Configurable<int> cfgNCrossedRowsTPC{"cfgNCrossedRowsTPC", 70, "Minimum TPC crossed rows"};
  Configurable<float> cfgChi2PerClusterTPC{"cfgChi2PerClusterTPC", 4.0, "Maximum TPC chi2 per cluster"};
  Configurable<float> cfgChi2PerClusterITS{"cfgChi2PerClusterITS", 36.0, "Maximum ITS chi2 per cluster"};
  Configurable<float> cfgDCAZ{"cfgDCAZ", 0.1, "Maximum absolute DCAz (cm)"};
  Configurable<bool> cfgRequireGoldenChi2{"cfgRequireGoldenChi2", true, "Require global/golden track selection"};

  // ========================================================================
  // PDG charge helper
  // ========================================================================
  bool isChargedParticle(int pdgCode)
  {
    const auto* pdgParticle = pdg->GetParticle(pdgCode);
    if (!pdgParticle) {
      return false;
    }
    if (std::abs(pdgParticle->Charge()) < ChargeEpsilon) {
      return false;
    }
    return true;
  }

  // ========================================================================
  // Rho / flattenicity calculation
  // ========================================================================
  template <std::size_t N>
  float computeRho(const std::array<float, N>& counts)
  {
    float total = 0.0f;
    for (std::size_t i = 0; i < N; ++i) {
      total += counts[i];
    }
    if (total <= 0.0f) {
      return -1.0f;
    }
    const float mean = total / static_cast<float>(N);
    if (mean <= 0.0f) {
      return -1.0f;
    }
    float sumSq = 0.0f;
    for (std::size_t i = 0; i < N; ++i) {
      const float diff = counts[i] - mean;
      sumSq += diff * diff;
    }
    return std::sqrt(sumSq) / (static_cast<float>(N) * mean);
  }

  template <std::size_t N>
  float computeFlattenicity(const std::array<float, N>& counts)
  {
    const float rho = computeRho(counts);
    if (rho < 0.0f) {
      return -1.0f;
    }
    return 1.0f - rho;
  }

  // ========================================================================
  // FT0 channel -> sector
  // ========================================================================
  int getFT0ASector(int channel)
  {
    if (channel < 0 || channel >= NchA) {
      return -1;
    }
    return channel / ChannelsPerSector;
  }
  int getFT0CSector(int channel)
  {
    if (channel < 0 || channel >= NchC) {
      return -1;
    }
    return channel / ChannelsPerSector;
  }

  // ========================================================================
  // Generator particle -> 208-cell grid
  // ========================================================================
  int assignToFT0Cell(float eta, float phi, bool& isFT0A)
  {
    const bool inFT0A = (eta > FT0AEtaMin && eta < FT0AEtaMax);
    const bool inFT0C = (eta > FT0CEtaMin && eta < FT0CEtaMax);
    if (!inFT0A && !inFT0C) {
      return -1;
    }
    isFT0A = inFT0A;

    // Wrap phi into [0, 2pi) using the O2-standard helper - avoids the
    // manual +=/-= TwoPI pattern the linter flags (two-pi-add-subtract).
    const float wrappedPhi = RecoDecay::constrainAngle(phi, 0.0f);

    const float phiWidth = o2::constants::math::TwoPI / static_cast<float>(NPhiSectors);
    int phiBin = static_cast<int>(std::floor(wrappedPhi / phiWidth));
    phiBin = std::max(0, std::min(phiBin, NPhiSectors - 1));

    if (inFT0A) {
      const float etaWidth = (FT0AEtaMax - FT0AEtaMin) / static_cast<float>(NEtaA);
      int etaBin = static_cast<int>(std::floor((eta - FT0AEtaMin) / etaWidth));
      etaBin = std::max(0, std::min(etaBin, NEtaA - 1));
      return etaBin * NPhiSectors + phiBin;
    }
    const float etaWidth = (FT0CEtaMax - FT0CEtaMin) / static_cast<float>(NEtaC);
    int etaBin = static_cast<int>(std::floor((eta - FT0CEtaMin) / etaWidth));
    etaBin = std::max(0, std::min(etaBin, NEtaC - 1));
    return NchA + etaBin * NPhiSectors + phiBin;
  }

  // ========================================================================
  // Compute both reco definitions
  // ========================================================================
  template <typename FT0>
  bool computeFT0Flattenicities(const FT0& ft0, float& flattenicityAverage, float& flattenicityCombined, float& sumAmplitude)
  {
    std::array<float, NSectorsA> countsA{};
    std::array<float, NSectorsC> countsC{};
    std::array<float, NSectorsCombined> countsCombined{};
    sumAmplitude = 0.0f;
    for (std::size_t i = 0; i < ft0.amplitudeA().size(); ++i) {
      const int channel = ft0.channelA()[i];
      const int sector = getFT0ASector(channel);
      const float amplitude = ft0.amplitudeA()[i];
      if (sector >= 0 && sector < NSectorsA) {
        countsA[sector] += amplitude;
        countsCombined[sector] += amplitude;
      }
      sumAmplitude += amplitude;
    }
    for (std::size_t i = 0; i < ft0.amplitudeC().size(); ++i) {
      const int channel = ft0.channelC()[i];
      const int sector = getFT0CSector(channel);
      const float amplitude = ft0.amplitudeC()[i];
      if (sector >= 0 && sector < NSectorsC) {
        countsC[sector] += amplitude;
        countsCombined[NSectorsA + sector] += amplitude;
      }
      sumAmplitude += amplitude;
    }
    const float rhoA = computeRho(countsA);
    const float rhoC = computeRho(countsC);
    if (rhoA < 0.0f || rhoC < 0.0f) {
      return false;
    }
    flattenicityAverage = 1.0f - 0.5f * (rhoA + rhoC);
    const float rhoCombined = computeRho(countsCombined);
    if (rhoCombined < 0.0f) {
      return false;
    }
    flattenicityCombined = 1.0f - rhoCombined;
    return true;
  }

  // ========================================================================
  // processMC
  // ========================================================================
  void processMC(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    if (std::abs(mcCollision.posZ()) > cfgVzMax) {
      return;
    }
    // Combined 208-cell truth grid (cells 0..NchA-1 = FT0A-like, NchA..NCell-1
    // = FT0C-like) AND the two separate sub-arrays, filled in parallel, so we
    // can compute BOTH truth definitions - average and combined - the same
    // way the reco side already does for the real FT0-A/FT0-C sectors.
    std::array<float, NCell> truthCounts{};
    std::array<float, NchA> truthCountsA{};
    std::array<float, NchC> truthCountsC{};
    int nChINEL = 0;
    bool hasFT0A = false, hasFT0C = false;

    for (const auto& particle : mcParticles) {
      if ((particle.flags() & NPhysicalPrimaryBit) == 0) {
        continue;
      }
      if (!isChargedParticle(particle.pdgCode())) {
        continue;
      }
      if (particle.pt() <= cfgPtMin) {
        continue; // Antonio: require pT > 0
      }
      if (std::abs(particle.eta()) < cfgEtaMax) {
        ++nChINEL;
      }

      const bool inFT0A = (particle.eta() > FT0AEtaMin && particle.eta() < FT0AEtaMax);
      const bool inFT0C = (particle.eta() > FT0CEtaMin && particle.eta() < FT0CEtaMax);
      if (inFT0A) {
        hasFT0A = true;
      }
      if (inFT0C) {
        hasFT0C = true;
      }

      bool isFT0A = false;
      const int cellId = assignToFT0Cell(particle.eta(), particle.phi(), isFT0A);
      if (cellId < 0 || cellId >= NCell) {
        continue;
      }
      truthCounts[cellId] += 1.0f;
      if (isFT0A) {
        truthCountsA[cellId] += 1.0f;
      } else {
        truthCountsC[cellId - NchA] += 1.0f;
      }

      histos.fill(HIST("Truth/hCellOccupancy"), cellId);
      if (isFT0A) {
        histos.fill(HIST("Truth/hCellOccupancyFT0A"), cellId);
      } else {
        histos.fill(HIST("Truth/hCellOccupancyFT0C"), cellId - NchA);
      }
    }

    if (nChINEL == 0) {
      return;
    }
    histos.fill(HIST("QA/hEvents"), 1);
    histos.fill(HIST("QA/hNchINEL"), nChINEL);

    // Require activity on BOTH sides of the detector (Antonio's instruction),
    // for MC truth exactly as already required on the reco side.
    if (!hasFT0A || !hasFT0C) {
      return;
    }
    histos.fill(HIST("QA/hEvents"), 2);

    // Definition 1: average of two separately-normalized rho's (mirrors the
    // reco "average" definition, using the truth-grid A/C sub-arrays).
    const float rhoATruth = computeRho(truthCountsA);
    const float rhoCTruth = computeRho(truthCountsC);
    if (rhoATruth >= 0.0f && rhoCTruth >= 0.0f) {
      const float flattenicityTruthAverage = 1.0f - 0.5f * (rhoATruth + rhoCTruth);
      if (flattenicityTruthAverage >= 0.0f && flattenicityTruthAverage <= 1.0f) {
        histos.fill(HIST("Truth/hFlattenicityTruthAverage_fine"), flattenicityTruthAverage);
        histos.fill(HIST("Truth/hFlattenicityTruthAverage"), flattenicityTruthAverage);
      }
    }

    // Definition 2: one combined rho over the full 208-cell truth grid.
    const float flattenicityTruthCombined = computeFlattenicity(truthCounts);
    if (flattenicityTruthCombined >= 0.0f) {
      histos.fill(HIST("Truth/hFlattenicityTruthCombined_fine"), flattenicityTruthCombined);
      histos.fill(HIST("Truth/hFlattenicityTruthCombined"), flattenicityTruthCombined);
    }
  }

  PROCESS_SWITCH(FlattenicityTask, processMC, "Process MC truth for Run-1 calibration", true);

  // ========================================================================
  // processData
  // ========================================================================
  void processData(CollisionsWithSelection::iterator const& collision,
                   aod::FT0s const& ft0s,
                   FullTracks const& tracks)
  {
    if (cfgApplySel8 && !collision.sel8()) {
      return;
    }
    if (std::abs(collision.posZ()) > cfgVzMax) {
      return;
    }

    int nChINEL = 0;
    for (const auto& track : tracks) {
      if (track.collisionId() != collision.globalIndex()) {
        continue;
      }
      if (track.pt() <= cfgPtMin) {
        continue;
      }
      if (std::abs(track.eta()) < cfgEtaMax) {
        ++nChINEL;
      }
    }
    if (nChINEL == 0) {
      return;
    }

    histos.fill(HIST("QA/hEvents"), 3);
    histos.fill(HIST("QA/hNchINEL"), nChINEL);

    bool foundFT0 = false;
    aod::FT0s::iterator ft0;
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

    float flattenicityAverage = -1.0f, flattenicityCombined = -1.0f, sumFT0M = 0.0f;
    if (!computeFT0Flattenicities(ft0, flattenicityAverage, flattenicityCombined, sumFT0M)) {
      return;
    }

    histos.fill(HIST("Reco/hFT0Mraw"), sumFT0M);
    histos.fill(HIST("Reco/hFT0Mraw_fine"), sumFT0M);

    if (flattenicityAverage >= 0.0f && flattenicityAverage <= 1.0f) {
      histos.fill(HIST("Reco/hFlattenicityFT0Average_fine"), flattenicityAverage);
      histos.fill(HIST("Reco/hFlattenicityFT0Average"), flattenicityAverage);
    }
    if (flattenicityCombined >= 0.0f && flattenicityCombined <= 1.0f) {
      histos.fill(HIST("Reco/hFlattenicityFT0Combined_fine"), flattenicityCombined);
      histos.fill(HIST("Reco/hFlattenicityFT0Combined"), flattenicityCombined);
    }
  }

  PROCESS_SWITCH(FlattenicityTask, processData, "Process data / MC-reco for Run-1 calibration", true);
};

// ============================================================================
// Workflow
// ============================================================================
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FlattenicityTask>(cfgc)};
}
