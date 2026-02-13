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

// **Hadronic Correction in the EMCAL framework: to avoid the double counting of the charged particles' contribution in jets**
/// \author Archita Rani Dash <archita.rani.dash@cern.ch>

#include "PWGJE/DataModel/EMCALClusterDefinition.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include "TVector2.h"
#include <TF1.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct EmcalClusterHadronicCorrectionTask {
  Produces<aod::JClustersCorrectedEnergies> clusterEnergyCorrectedTable;

  HistogramRegistry registry;
  // Configurables for Histogram Binning
  PresliceUnsorted<aod::JEMCTracks> perTrackMatchedTrack = aod::jemctrack::trackId;

  // define configurables here
  Configurable<float> minTrackPt{"minTrackPt", 0.15, "Minimum pT for tracks"};
  Configurable<double> hadCorr1{"hadCorr1", 1., "hadronic correction fraction for complete cluster energy subtraction for one matched track"};                // 100% - default
  Configurable<double> hadCorr2{"hadCorr2", 0.7, "hadronic correction fraction for systematic studies for one matched track"};                                // 70%
  Configurable<double> hadCorralltrks1{"hadCorralltrks1", 1., "hadronic correction fraction for complete cluster energy subtraction for all matched tracks"}; // 100% - all tracks
  Configurable<double> hadCorralltrks2{"hadCorralltrks2", 0.7, "hadronic correction fraction for systematic studies for all matched tracks"};                 // 70%
  Configurable<float> minDEta{"minDEta", 0.01, "Minimum dEta between track and cluster"};
  Configurable<float> minDPhi{"minDPhi", 0.01, "Minimum dPhi between track and cluster"};
  Configurable<double> constantSubtractionValue{"constantSubtractionValue", 0.236, "Value to be used for constant MIP subtraction (only applicable if using constant subtraction in M02 scheme)"};

  // pT-dependent track-matching configurables
  Configurable<float> eta0{"eta0", 0.04, "Param 0 in eta for pt-dependent matching"};
  Configurable<float> eta1{"eta1", 0.010, "Param 1 in eta for pt-dependent matching"};
  Configurable<float> eta2{"eta2", 2.5, "Param 2 in eta for pt-dependent matching"};
  Configurable<float> phi0{"phi0", 0.09, "Param 0 in phi for pt-dependent matching"};
  Configurable<float> phi1{"phi1", 0.015, "Param 1 in phi for pt-dependent matching"};
  Configurable<float> phi2{"phi2", 2.0, "Param 2 in phi for pt-dependent matching"};

  Configurable<bool> doHadCorrSyst{"doHadCorrSyst", false, "Do hadronic correction for systematic studies"};
  Configurable<bool> doMomDepMatching{"doMomDepMatching", true, "Do momentum dependent track matching"}; // to be always set to true in Run 3
  Configurable<bool> useM02SubtractionScheme1{"useM02SubtractionScheme1", false, "Flag to enable hadronic correction scheme using cluster M02 value for clusterE1 and clusterEAll1"};
  Configurable<bool> useM02SubtractionScheme2{"useM02SubtractionScheme2", false, "Flag to enable hadronic correction scheme using cluster M02 value for clusterE2 and clusterEAll2"};
  Configurable<bool> useFraction1{"useFraction1", false, "Fractional momentum subtraction for clusterE1 and clusterEAll1"};
  Configurable<bool> useFraction2{"useFraction2", false, "Fractional momentum subtraction for clusterE2 and clusterEAll2"};

  void init(o2::framework::InitContext&)
  {
    // Event histograms
    registry.add("h_allcollisions", "Total events; event status;entries", {HistType::kTH1F, {{1, 0.5, 1.5}}});

    // Matched-Cluster histograms
    registry.add("h_matchedclusters", "Total matched clusters; cluster status;entries", {HistType::kTH1F, {{1, 0.5, 1.5}}});
    registry.add("h_ClsE", "; Cls E w/o correction (GeV); entries", {HistType::kTH1F, {{350, 0., 350.}}});
    registry.add("h_Ecluster1", "; Ecluster1 (GeV); entries", {HistType::kTH1F, {{350, 0., 350.}}});
    registry.add("h_Ecluster2", "; Ecluster2 (GeV); entries", {HistType::kTH1F, {{350, 0., 350.}}});
    registry.add("h_EclusterAll1", "; EclusterAll1 (GeV); entries", {HistType::kTH1F, {{350, 0., 350.}}});
    registry.add("h_EclusterAll2", "; EclusterAll2 (GeV); entries", {HistType::kTH1F, {{350, 0., 350.}}});
    registry.add("h_ClsTime", "Cluster time distribution of uncorrected cluster E; #it{t}_{cls} (ns); entries", {HistType::kTH1F, {{1500, -600., 900.}}});
    registry.add("h_ClsM02", "Cluster M02 distribution of uncorrected cluster E; #it{M}_{02}; entries", {HistType::kTH1F, {{400, 0., 5.}}});
    registry.add("h2_ClsEvsNmatches", "Original cluster energy vs Nmatches; Cls E w/o correction (GeV); Nmatches", {HistType::kTH2F, {{350, 0., 350.}, {100, -0.5, 21.}}});
    registry.add("h2_ClsEvsEcluster1", "; Cls E w/o correction (GeV); Ecluster1 (GeV)", {HistType::kTH2F, {{350, 0., 350.}, {350, 0., 350.}}});
    registry.add("h2_ClsEvsEcluster2", "; Cls E w/o correction (GeV); Ecluster2 (GeV)", {HistType::kTH2F, {{350, 0., 350.}, {350, 0., 350.}}});
    registry.add("h2_ClsEvsEclusterAll1", "; Cls E w/o correction (GeV); EclusterAll1 (GeV)", {HistType::kTH2F, {{350, 0., 350.}, {350, 0., 350.}}});
    registry.add("h2_ClsEvsEclusterAll2", "; Cls E w/o correction (GeV); EclusterAll2 (GeV)", {HistType::kTH2F, {{350, 0., 350.}, {350, 0., 350.}}});

    // Matched-Track histograms
    registry.add("h_matchedtracks", "Total matched tracks; track status;entries", {HistType::kTH1F, {{1, 0.5, 1.5}}});
  }

  // The matching of clusters and tracks is already centralised in the EMCAL framework.
  // One only needs to apply a filter on matched clusters
  // Here looping over all collisions matched to EMCAL clusters
  void processMatchedCollisions(aod::JetCollision const&, soa::Join<aod::JClusters, aod::JClusterTracks> const& clusters, aod::JEMCTracks const& emcTracks, aod::JetTracks const&)
  {
    registry.fill(HIST("h_allcollisions"), 1);

    // skip events with no clusters
    if (clusters.size() == 0) {
      return;
    }

    // Looping over all clusters matched to the collision
    for (const auto& cluster : clusters) {

      registry.fill(HIST("h_matchedclusters"), 1);

      double clusterE1;
      double clusterE2;
      double clusterEAll1;
      double clusterEAll2;
      clusterE1 = clusterE2 = clusterEAll1 = clusterEAll2 = cluster.energy();

      registry.fill(HIST("h_ClsE"), cluster.energy());
      registry.fill(HIST("h_ClsM02"), cluster.m02());
      registry.fill(HIST("h_ClsTime"), cluster.time());

      int nMatches = 0;         // counter for closest matched track
      double closestTrkP = 0.0; // closest track momentum
      double totalTrkP = 0.0;   // total track momentum

      // pT-dependent track-matching instead of PID based track-matching to be adapted from Run 2 - suggested by Markus Fasel

      TF1 funcPtDepEta("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      funcPtDepEta.SetParameters(eta0, eta1, eta2);
      TF1 funcPtDepPhi("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      funcPtDepPhi.SetParameters(phi0, phi1, phi2);

      // No matched tracks (trackless case)
      if (cluster.matchedTracks().size() == 0) {
        // Use original cluster energy values, no subtraction needed.
        registry.fill(HIST("h2_ClsEvsNmatches"), cluster.energy(), 0);
        registry.fill(HIST("h_Ecluster1"), clusterE1);
        registry.fill(HIST("h_Ecluster2"), clusterE2);
        registry.fill(HIST("h_EclusterAll1"), clusterEAll1);
        registry.fill(HIST("h_EclusterAll2"), clusterEAll2);
        registry.fill(HIST("h2_ClsEvsEcluster1"), cluster.energy(), clusterE1);
        registry.fill(HIST("h2_ClsEvsEcluster2"), cluster.energy(), clusterE2);
        registry.fill(HIST("h2_ClsEvsEclusterAll1"), cluster.energy(), clusterEAll1);
        registry.fill(HIST("h2_ClsEvsEclusterAll2"), cluster.energy(), clusterEAll2);
        clusterEnergyCorrectedTable(clusterE1, clusterE2, clusterEAll1, clusterEAll2);
        continue;
      }

      // Looping over all matched tracks for the cluster
      // Total number of matched tracks = 20 (hard-coded)
      for (const auto& matchedTrack : cluster.matchedTracks_as<aod::JetTracks>()) {
        if (matchedTrack.pt() < minTrackPt) {
          continue;
        }
        double mom = std::abs(matchedTrack.p());
        registry.fill(HIST("h_matchedtracks"), 1);

        // CASE 1: skip tracks with a very low pT
        constexpr double kMinMom = 1e-6;
        if (mom < kMinMom) {
          continue;
        } // end CASE 1

        // CASE 2:
        //  a) If one matched track -> 100% energy subtraction
        //  b) If more than one matched track -> 100% energy subtraction
        //  c) If you want to do systematic studies -> perform the above two checks a) and b), and then subtract 70% energy instead of 100%

        // Perform dEta/dPhi matching
        auto emcTrack = (emcTracks.sliceBy(perTrackMatchedTrack, matchedTrack.globalIndex())).iteratorAt(0);
        double dEta = emcTrack.etaDiff();
        double dPhi = emcTrack.phiDiff();

        // Apply the eta and phi matching thresholds
        // dEta and dPhi cut : ensures that the matched track is within the desired eta/phi window

        // Do pT-dependent track matching
        if (doMomDepMatching) {
          auto trackEtaMax = funcPtDepEta.Eval(mom);
          auto trackPhiHigh = +funcPtDepPhi.Eval(mom);
          auto trackPhiLow = -funcPtDepPhi.Eval(mom);

          if ((dPhi < trackPhiHigh && dPhi > trackPhiLow) && std::fabs(dEta) < trackEtaMax) {
            if (nMatches == 0) {
              closestTrkP = mom;
            }
            totalTrkP += mom;
            nMatches++;
          }
        } else {
          // Do fixed dEta/dPhi matching (non-pT dependent)
          if (std::fabs(dEta) >= minDEta || std::fabs(dPhi) >= minDPhi) {
            continue; // Skip this track if outside the fixed cut region
          }

          // If track passes fixed dEta/dPhi cuts, process it
          if (nMatches == 0) {
            closestTrkP = mom; // Closest track match
          }
          totalTrkP += mom; // Accumulate momentum
          nMatches++;       // Count this match
        }

      } // End of track loop
      registry.fill(HIST("h2_ClsEvsNmatches"), cluster.energy(), nMatches);

      if (nMatches >= 1) {
        if (useM02SubtractionScheme1) {
          // Do M02-based correction if enabled
          clusterE1 = subtractM02ClusterEnergy(cluster.m02(), clusterE1, nMatches, closestTrkP, hadCorr1, useFraction1);
          clusterEAll1 = subtractM02ClusterEnergy(cluster.m02(), clusterEAll1, nMatches, totalTrkP, hadCorralltrks1, useFraction1);
        } else {
          // Default energy subtraction (100% and 70%)
          clusterE1 = subtractClusterEnergy(clusterE1, closestTrkP, hadCorr1, nMatches, useFraction1);
          clusterEAll1 = subtractClusterEnergy(clusterEAll1, totalTrkP, hadCorralltrks1, nMatches, useFraction1);
        }

        if (useM02SubtractionScheme2) {
          // Do M02-based correction if enabled
          clusterE2 = subtractM02ClusterEnergy(cluster.m02(), clusterE2, nMatches, closestTrkP, hadCorr2, useFraction2);
          clusterEAll2 = subtractM02ClusterEnergy(cluster.m02(), clusterEAll2, nMatches, totalTrkP, hadCorralltrks2, useFraction2);
        } else {
          // Default energy subtraction (100% and 70%)
          clusterE2 = subtractClusterEnergy(clusterE2, closestTrkP, hadCorr2, nMatches, useFraction2);
          clusterEAll2 = subtractClusterEnergy(clusterEAll2, totalTrkP, hadCorralltrks2, nMatches, useFraction2);
        }
      }
      registry.fill(HIST("h_Ecluster1"), clusterE1);
      registry.fill(HIST("h_Ecluster2"), clusterE2);
      registry.fill(HIST("h_EclusterAll1"), clusterEAll1);
      registry.fill(HIST("h_EclusterAll2"), clusterEAll2);
      registry.fill(HIST("h2_ClsEvsEcluster1"), cluster.energy(), clusterE1);
      registry.fill(HIST("h2_ClsEvsEcluster2"), cluster.energy(), clusterE2);
      registry.fill(HIST("h2_ClsEvsEclusterAll1"), cluster.energy(), clusterEAll1);
      registry.fill(HIST("h2_ClsEvsEclusterAll2"), cluster.energy(), clusterEAll2);

      // Fill the table with all four corrected energies
      clusterEnergyCorrectedTable(clusterE1, clusterE2, clusterEAll1, clusterEAll2);
    } // End of cluster loop
  } // process function ends
  PROCESS_SWITCH(EmcalClusterHadronicCorrectionTask, processMatchedCollisions, "hadronic correction", true);

  // Helper function to prevent negative energy subtraction
  double subtractClusterEnergy(double clusterE, double mom, double correctionFactor, int nMatches, bool useFraction)
  {
    double corrE = clusterE;

    // if (UseConstantSubtractionValue) {
    if (!useFraction) {
      corrE = clusterE - constantSubtractionValue * nMatches; // Use constant value for MIP-subtraction (regardless of the cluster-shape)
    } else {
      corrE = clusterE - correctionFactor * mom; // Fractional momentum subtraction
    }
    return (corrE < 0) ? 0 : corrE;
  }

  // Helper function for M02-based energy subtraction (based on cluster-shape)
  double subtractM02ClusterEnergy(double m02, double clusterE, int nMatches, double totalTrkP, double correctionFactor, bool useFraction)
  {
    double corrE = clusterE;

    // For M02 in the single photon region, the signal is primarily: Single photons, single electrons, single MIPs
    if (m02 > 0.1 && m02 < 0.4) { // circular clusters(electron/photon)
      corrE = 0;                  // Single electron, single MIP
    } else if (m02 > 0.4) {
      // Large M02 region (M02 > 0.4), more complex overlaps and hadronic showers.
      // The signal is primarily: Single hadronic shower, photon-photon overlap, photon-MIP overlap, MIP-MIP overlap,
      // MIP-hadronic shower overlap, hadronic shower - hadronic shower overlap)

      if (!useFraction) {
        corrE = clusterE - constantSubtractionValue * nMatches; // Use constant value for MIP-subtraction (regardless of the cluster-shape)
      } else {
        corrE = clusterE - correctionFactor * totalTrkP; // Fractional momentum subtraction
      }
    }
    return (corrE < 0) ? 0 : corrE; // Prevent negative energy
  }

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EmcalClusterHadronicCorrectionTask>(cfgc, TaskName{"emcal-cluster-hadronic-correction-task"})}; }
