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
///
/// \file jFlucEfficiencyTask.cxx
/// \brief Task to calculate the efficiency of the cf-derived tracks/particles
/// \author DongJo Kim, Jasper Parkkila, Bong-Hwi Lim (djkim@cern.ch, jparkkil@cern.ch, bong-hwi.lim@cern.ch)
/// \since March 2024

#include <vector>
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JFlucEfficiencyTask {
  // Add the pT binning array as a static member
  static constexpr std::array<double, 94> PttJacek = {
    0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
    0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
    1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
    2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
    4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0,
    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
    26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0,
    70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0,
    170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0,
    270.0, 280.0, 290.0, 300.0};

  // Update the axisPt configuration with proper vector initialization
  ConfigurableAxis axisPt{"axisPt", std::vector<double>(PttJacek.begin(), PttJacek.end()), "pT axis"};

  // Configurable for track selection
  Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Minimum transverse momentum"};
  Configurable<float> cfgPtMax{"cfgPtMax", 300.0f, "Maximum transverse momentum"};
  Configurable<float> cfgEtaMin{"cfgEtaMin", -1.0f, "Minimum pseudorapidity"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", 1.0f, "Maximum pseudorapidity"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Vertex cut"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0.0f, "Min centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 100.0f, "Max centrality"};
  Configurable<uint8_t> cfgTrackBitMask{"cfgTrackBitMask", 0, "BitMask for track selection systematics"};

  // Configurable axes
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "multiplicity / centrality axis"};

  // Filter declarations
  Filter cfCollisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex) &&
                             (aod::cfcollision::multiplicity > cfgCentMin) &&
                             (aod::cfcollision::multiplicity < cfgCentMax);
  Filter cfMCCollisionFilter = (nabs(aod::mccollision::posZ) < cfgCutVertex) &&
                               (aod::cfmccollision::multiplicity > cfgCentMin) &&
                               (aod::cfmccollision::multiplicity < cfgCentMax);
  Filter cfMCParticleFilter = (aod::cfmcparticle::pt >= cfgPtMin) &&
                              (aod::cfmcparticle::pt <= cfgPtMax) &&
                              (aod::cfmcparticle::eta >= cfgEtaMin) &&
                              (aod::cfmcparticle::eta <= cfgEtaMax);
  Filter cfTrackFilter = (aod::cftrack::pt >= cfgPtMin) &&
                         (aod::cftrack::pt <= cfgPtMax) &&
                         (aod::cftrack::eta >= cfgEtaMin) &&
                         (aod::cftrack::eta <= cfgEtaMax) &&
                         ((aod::track::trackType & (uint8_t)cfgTrackBitMask) == (uint8_t)cfgTrackBitMask);

  Filter trackFilter = (nabs(aod::track::eta) < cfgEtaMax) &&
                       (aod::track::pt > cfgPtMin) &&
                       ((requireGlobalTrackInFilter()) ||
                        (aod::track::isGlobalTrackSDD == (uint8_t) true));

  Configurable<bool> cfgEfficiencyFromData{"cfgEfficiencyFromData", false, "Calculate efficiency using data events as reference"};
  Configurable<int> cfgVerbosity{"cfgVerbosity", 0, "Verbosity level"};

  // Histogram Registry
  HistogramRegistry registry{
    "registry",
    {{"hEventCounterMC", "Event counter MC;Counter;Counts", {HistType::kTH1F, {{3, -0.5, 2.5}}}},
     {"hEventCounterReco", "Event counter Reco;Counter;Counts", {HistType::kTH1F, {{3, -0.5, 2.5}}}},
     {"hZVertexMC", "MC Z vertex distribution;Z vertex (cm);Centrality (%)", {HistType::kTH2F, {{200, -20, 20}, {axisMultiplicity}}}},
     {"hZVertexReco", "Reconstructed Z vertex distribution;Z vertex (cm);Centrality (%)", {HistType::kTH2F, {{200, -20, 20}, {axisMultiplicity}}}},
     {"hZVertexCorrelation", "Z vertex correlation;MC Z vertex (cm);Reco Z vertex (cm)", {HistType::kTH2F, {{200, -20, 20}, {200, -20, 20}}}}}};

  // Configurable for debugging
  Configurable<bool> debugMode{"debugMode", false, "Debug mode"};

  void init(InitContext const&)
  {
    if (debugMode) {
      LOGF(info, "Initializing JFlucEfficiencyTask");
    }

    if (doprocessMC) {
      registry.add("hPtGen", "Generated p_{T} (all);p_{T} (GeV/c);Centrality (%);Counts",
                   o2::framework::HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)});
      registry.add("hEtaGen", "Generated #eta (all);#eta;Centrality (%);Counts",
                   o2::framework::HistType::kTH2F, {AxisSpec(100, -1, 1), AxisSpec(axisMultiplicity)});
      registry.add("hPtGenPos", "Generated p_{T} (positive);p_{T} (GeV/c);Centrality (%);Counts",
                   o2::framework::HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)});

      registry.add("hPtGenNeg", "Generated p_{T} (negative);p_{T} (GeV/c);Centrality (%);Counts",
                   o2::framework::HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)});
    }

    if (doprocessData) {
      registry.add("hPtRec", "Reconstructed p_{T} (all);p_{T} (GeV/c);Centrality (%);Counts",
                   o2::framework::HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)});

      registry.add("hEtaRec", "Reconstructed #eta (all);#eta;Centrality (%);Counts",
                   o2::framework::HistType::kTH2F, {AxisSpec(100, -1, 1), AxisSpec(axisMultiplicity)});

      registry.add("hPtRecPos", "Reconstructed p_{T} (positive);p_{T} (GeV/c);Centrality (%);Counts",
                   o2::framework::HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)});

      registry.add("hPtRecNeg", "Reconstructed p_{T} (negative);p_{T} (GeV/c);Centrality (%);Counts",
                   o2::framework::HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)});
    }

    if (cfgEfficiencyFromData) {
      registry.add("hPtGenData", "Generated p_{T} from data events (all);p_{T} (GeV/c);Centrality (%);Counts",
                   {HistType::kTH2F, {axisPt, axisMultiplicity}});
      registry.add("hEtaGenData", "Generated #eta from data events (all);#eta;Centrality (%);Counts",
                   {HistType::kTH2F, {AxisSpec(100, -1, 1), axisMultiplicity}});
      registry.add("hPtGenDataPos", "Generated p_{T} from data events (positive);p_{T} (GeV/c);Centrality (%);Counts",
                   {HistType::kTH2F, {axisPt, axisMultiplicity}});
      registry.add("hPtGenDataNeg", "Generated p_{T} from data events (negative);p_{T} (GeV/c);Centrality (%);Counts",
                   {HistType::kTH2F, {axisPt, axisMultiplicity}});
      registry.add("hPtRecData", "Reconstructed p_{T} (all);p_{T} (GeV/c);Centrality (%);Counts",
                   o2::framework::HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)});
      registry.add("hEtaRecData", "Reconstructed #eta (all);#eta;Centrality (%);Counts",
                   o2::framework::HistType::kTH2F, {AxisSpec(100, -1, 1), AxisSpec(axisMultiplicity)});
      registry.add("hPtRecDataPos", "Reconstructed p_{T} (positive);p_{T} (GeV/c);Centrality (%);Counts",
                   o2::framework::HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)});
      registry.add("hPtRecDataNeg", "Reconstructed p_{T} (negative);p_{T} (GeV/c);Centrality (%);Counts",
                   o2::framework::HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)});
    }

    // Initialize histogram labels
    auto h1 = registry.get<TH1>(HIST("hEventCounterMC"));
    auto h2 = registry.get<TH1>(HIST("hEventCounterReco"));

    if (h1 && h2) {
      h1->GetXaxis()->SetBinLabel(1, "All MC Events");
      h1->GetXaxis()->SetBinLabel(2, "Selected MC Events");
      h1->GetXaxis()->SetBinLabel(3, "Analyzed MC Events");

      h2->GetXaxis()->SetBinLabel(1, "All Reco Events");
      h2->GetXaxis()->SetBinLabel(2, "Selected Reco Events");
      h2->GetXaxis()->SetBinLabel(3, "Analyzed Reco Events");
    } else {
      LOGF(error, "Failed to get histograms from registry");
    }
  }

  void processMC(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, soa::Filtered<aod::CFMcParticles> const& mcParticles)
  {
    float centrality = mcCollision.multiplicity();

    for (const auto& particle : mcParticles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }

      registry.fill(HIST("hPtGen"), particle.pt(), centrality);
      registry.fill(HIST("hEtaGen"), particle.eta(), centrality);

      if (particle.sign() > 0) { // Positive particles
        registry.fill(HIST("hPtGenPos"), particle.pt(), centrality);
      } else if (particle.sign() < 0) { // Negative particles
        registry.fill(HIST("hPtGenNeg"), particle.pt(), centrality);
      }
    }
  }

  void processData(soa::Filtered<aod::CFCollisions>::iterator const& cfCollision, soa::Filtered<aod::CFTracks> const& cfTracks)
  {
    float centrality = cfCollision.multiplicity();

    if (centrality < cfgCentMin || centrality > cfgCentMax) {
      return;
    }

    for (const auto& track : cfTracks) {
      registry.fill(HIST("hPtRec"), track.pt(), centrality);
      registry.fill(HIST("hEtaRec"), track.eta(), centrality);

      if (track.sign() > 0) { // Positive tracks
        registry.fill(HIST("hPtRecPos"), track.pt(), centrality);
      } else if (track.sign() < 0) { // Negative tracks
        registry.fill(HIST("hPtRecNeg"), track.pt(), centrality);
      }
    }
  }

  template <typename TCollision, typename TTracks>
  void fillQA(const TCollision& /*collision*/, float multiplicity, const TTracks& tracks)
  {
    registry.fill(HIST("multiplicity"), multiplicity);
    for (const auto& track : tracks) {
      registry.fill(HIST("yields"), multiplicity, track.pt(), track.eta());
      registry.fill(HIST("etaphi"), track.eta(), track.phi());
    }
  }

  // NOTE SmallGroups includes soa::Filtered always
  Preslice<aod::CFTracksWithLabel> perCollision = aod::cftrack::cfCollisionId;
  void processEfficiency(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision,
                         aod::CFMcParticles const& mcParticles,
                         soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions,
                         aod::CFTracksWithLabel const& tracks)
  {
    try {
      // Count MC events and fill MC z-vertex with centrality
      registry.fill(HIST("hEventCounterMC"), 0);
      registry.fill(HIST("hZVertexMC"), mcCollision.posZ(), mcCollision.multiplicity());

      if (debugMode) {
        LOGF(info, "Processing MC collision %d at z = %.3f", mcCollision.globalIndex(), mcCollision.posZ());
      }

      // Fill MC particle histograms
      for (const auto& mcParticle : mcParticles) {
        if (!mcParticle.isPhysicalPrimary())
          continue;

        // Fill generated particle histograms
        registry.fill(HIST("hPtGen"), mcParticle.pt(), mcCollision.multiplicity());
        registry.fill(HIST("hEtaGen"), mcParticle.eta(), mcCollision.multiplicity());

        if (mcParticle.sign() > 0) {
          registry.fill(HIST("hPtGenPos"), mcParticle.pt(), mcCollision.multiplicity());
        } else if (mcParticle.sign() < 0) {
          registry.fill(HIST("hPtGenNeg"), mcParticle.pt(), mcCollision.multiplicity());
        }

        if (cfgEfficiencyFromData) {
          registry.fill(HIST("hPtGenData"), mcParticle.pt(), mcCollision.multiplicity());
          registry.fill(HIST("hEtaGenData"), mcParticle.eta(), mcCollision.multiplicity());
          if (mcParticle.sign() > 0) {
            registry.fill(HIST("hPtGenDataPos"), mcParticle.pt(), mcCollision.multiplicity());
          } else if (mcParticle.sign() < 0) {
            registry.fill(HIST("hPtGenDataNeg"), mcParticle.pt(), mcCollision.multiplicity());
          }
        }
      }

      registry.fill(HIST("hEventCounterMC"), 1);

      // Check reconstructed collisions
      if (collisions.size() == 0) {
        if (debugMode) {
          LOGF(info, "No reconstructed collisions found for MC collision %d", mcCollision.globalIndex());
        }
        return;
      }

      // Process reconstructed events
      for (const auto& collision : collisions) {
        registry.fill(HIST("hEventCounterReco"), 0);
        registry.fill(HIST("hZVertexReco"), collision.posZ(), collision.multiplicity());
        registry.fill(HIST("hZVertexCorrelation"), mcCollision.posZ(), collision.posZ());

        if (debugMode) {
          LOGF(info, "Processing reconstructed collision %d at z = %.3f",
               collision.globalIndex(), collision.posZ());
        }

        // Fill track histograms
        auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
        for (const auto& track : groupedTracks) {
          if (!track.has_cfMCParticle()) {
            if (debugMode) {
              LOGF(debug, "Track without MC particle found");
            }
            continue;
          }

          registry.fill(HIST("hPtRecData"), track.pt(), collision.multiplicity());
          registry.fill(HIST("hEtaRecData"), track.eta(), collision.multiplicity());

          if (track.sign() > 0) {
            registry.fill(HIST("hPtRecDataPos"), track.pt(), collision.multiplicity());
          } else if (track.sign() < 0) {
            registry.fill(HIST("hPtRecDataNeg"), track.pt(), collision.multiplicity());
          }
        }

        // Count selected and analyzed events
        registry.fill(HIST("hEventCounterReco"), 1);
        registry.fill(HIST("hEventCounterReco"), 2);
      }

      registry.fill(HIST("hEventCounterMC"), 2);

    } catch (const std::exception& e) {
      LOGF(error, "Exception caught in processEfficiency: %s", e.what());
    } catch (...) {
      LOGF(error, "Unknown exception caught in processEfficiency");
    }
  }

  PROCESS_SWITCH(JFlucEfficiencyTask, processMC, "Process MC only", false);
  PROCESS_SWITCH(JFlucEfficiencyTask, processData, "Process data only", false);
  PROCESS_SWITCH(JFlucEfficiencyTask, processEfficiency, "Process efficiency task", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JFlucEfficiencyTask>(cfgc)};
}
