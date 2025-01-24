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

#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <TF1.h>
#include <vector>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "DetectorsBase/GeometryManager.h"

#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/EMCALMatchedTracks.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"
#include "PWGJE/DataModel/emcalClusterHadronicCorrectionTask.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"
#include "EMCALBase/Geometry.h"
#include "EMCALBase/ClusterFactory.h"
#include "EMCALBase/NonlinearityHandler.h"
#include "EMCALReconstruction/Clusterizer.h"
#include "PWGJE/Core/JetUtilities.h"
#include "TVector2.h"

#include "CommonDataFormat/InteractionRecord.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using collisionEvSelIt = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator;
using myTracks = o2::soa::Filtered<o2::soa::Join<o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::FullTracks, o2::aod::TrackSelection>>;

struct EmcalClusterHadronicCorrectionTask {
  Produces<o2::aod::EmcalHCs> hadroniccorrectedclusters;

  HistogramRegistry registry;
  // Configurables for Histogram Binning
  Preslice<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalclustercell::emcalclusterId; // looking at clusterID column in the EMCALMatchedTracks for every cluster

  // define configurables here
  Configurable<int> mClusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  Configurable<float> minTime{"minTime", -25., "Minimum cluster time for time cut"};
  Configurable<float> maxTime{"maxTime", 20., "Maximum cluster time for time cut"};
  Configurable<float> minM02{"minM02", 0.1, "Minimum M02 for M02 cut"};
  Configurable<float> maxM02{"maxM02", 0.9, "Maximum M02 for M02 cut"};
  Configurable<float> minTrackPt{"minTrackPt", 0.15, "Minimum pT for tracks"};
  Configurable<double> fHadCorr1{"HadCorr1", 1., "hadronic correction fraction for complete cluster energy subtraction for one matched track"};                // 100% - default
  Configurable<double> fHadCorr2{"HadCorr2", 0.7, "hadronic correction fraction for systematic studies for one matched track"};                                // 70%
  Configurable<double> fHadCorralltrks1{"HadCorralltrks1", 1., "hadronic correction fraction for complete cluster energy subtraction for all matched tracks"}; // 100% - all tracks
  Configurable<double> fHadCorralltrks2{"HadCorralltrks2", 0.7, "hadronic correction fraction for systematic studies for all matched tracks"};                 // 70%
  Configurable<float> minDEta{"minDEta", 0.01, "Minimum dEta between track and cluster"};
  Configurable<float> minDPhi{"minDPhi", 0.01, "Minimum dPhi between track and cluster"};
  Configurable<double> fConstantSubtractionValue{"ConstantSubtractionValue", 0.236, "Value to be used for constant MIP subtraction (only applicable if using constant subtraction in M02 scheme)"};

  // pT-dependent track-matching configurables
  Configurable<float> Eta0{"eta0", 0.04, "Param 0 in eta for pt-dependent matching"};
  Configurable<float> Eta1{"eta1", 0.010, "Param 1 in eta for pt-dependent matching"};
  Configurable<float> Eta2{"eta2", 2.5, "Param 2 in eta for pt-dependent matching"};
  Configurable<float> Phi0{"phi0", 0.09, "Param 0 in phi for pt-dependent matching"};
  Configurable<float> Phi1{"phi1", 0.015, "Param 1 in phi for pt-dependent matching"};
  Configurable<float> Phi2{"phi2", 2.0, "Param 2 in phi for pt-dependent matching"};
  Configurable<double> PhiMatch{"phiMatch", 0.050, "phi match value in pp"};
  Configurable<double> EtaMatch{"etaMatch", 0.025, "eta match value in pp"};

  Configurable<bool> doHadCorrSyst{"doHadCorrSyst", false, "Do hadronic correction for systematic studies"};
  Configurable<bool> doMomDepMatching{"doMomDepMatching", true, "Do momentum dependent track matching"}; // to be always set to true in Run 3
  Configurable<bool> UseM02SubtractionScheme1{"UseM02SubtractionScheme1", false, "Flag to enable hadronic correction scheme using cluster M02 value for Ecluster1 and EclusterAll1"};
  Configurable<bool> UseM02SubtractionScheme2{"UseM02SubtractionScheme2", false, "Flag to enable hadronic correction scheme using cluster M02 value for Ecluster2 and EclusterAll2"};
  Configurable<bool> UseFraction1{"UseFraction1", false, "Fractional momentum subtraction for Ecluster1 and EclusterAll1"};
  Configurable<bool> UseFraction2{"UseFraction2", false, "Fractional momentum subtraction for Ecluster2 and EclusterAll2"};

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
    registry.add("h_ClsTime", "Cluster time distribution of uncorrected cluster E; #it{t}_{cls} (ns); entries", {HistType::kTH1F, {{500, -250., 250.}}});
    registry.add("h_ClsM02", "Cluster M02 distribution of uncorrected cluster E; #it{M}_{02}; entries", {HistType::kTH1F, {{400, 0., 5.}}});
    registry.add("h2_ClsEvsNmatches", "Original cluster energy vs Nmatches; Cls E w/o correction (GeV); Nmatches", {HistType::kTH2F, {{350, 0., 350.}, {100, -0.5, 21.}}});
    registry.add("h2_ClsEvsEcluster1", "; Cls E w/o correction (GeV); Ecluster1 (GeV)", {HistType::kTH2F, {{350, 0., 350.}, {350, 0., 350.}}});
    registry.add("h2_ClsEvsEcluster2", "; Cls E w/o correction (GeV); Ecluster2 (GeV)", {HistType::kTH2F, {{350, 0., 350.}, {350, 0., 350.}}});
    registry.add("h2_ClsEvsEclusterAll1", "; Cls E w/o correction (GeV); EclusterAll1 (GeV)", {HistType::kTH2F, {{350, 0., 350.}, {350, 0., 350.}}});
    registry.add("h2_ClsEvsEclusterAll2", "; Cls E w/o correction (GeV); EclusterAll2 (GeV)", {HistType::kTH2F, {{350, 0., 350.}, {350, 0., 350.}}});

    // Matched-Track histograms
    registry.add("h_matchedtracks", "Total matched tracks; track status;entries", {HistType::kTH1F, {{1, 0.5, 1.5}}});
  }

  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == mClusterDefinition);
  Filter trackSelection = (o2::aod::track::pt >= minTrackPt);
  // The matching of clusters and tracks is already centralised in the EMCAL framework.
  // One only needs to apply a filter on matched clusters
  // Here looping over all collisions matched to EMCAL clusters
  void processMatchedCollisions(collisionEvSelIt const&, o2::aod::EMCALClusters const& clusters, o2::aod::EMCALMatchedTracks const& matchedtracks, myTracks const&)
  {
    registry.fill(HIST("h_allcollisions"), 1);

    // skip events with no clusters
    if (clusters.size() == 0) {
      return;
    }

    // Looping over all clusters matched to the collision
    for (const auto& cluster : clusters) {
      registry.fill(HIST("h_matchedclusters"), 1);

      double Ecluster1;
      double Ecluster2;
      double EclusterAll1;
      double EclusterAll2;
      Ecluster1 = Ecluster2 = EclusterAll1 = EclusterAll2 = cluster.energy();

      registry.fill(HIST("h_ClsE"), cluster.energy());
      registry.fill(HIST("h_ClsM02"), cluster.m02());
      registry.fill(HIST("h_ClsTime"), cluster.time());
      // selecting ALL MATCHED TRACKS after slicing all entries in perClusterMatchedTracks by the cluster globalIndex
      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex());

      int Nmatches = 0;         // counter for closest matched track
      double closestTrkP = 0.0; // closest track momentum
      double totalTrkP = 0.0;   // total track momentum

      // pT-dependent track-matching instead of PID based track-matching to be adapted from Run 2 - suggested by Markus Fasel

      TF1 funcPtDepEta("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      funcPtDepEta.SetParameters(Eta0, Eta1, Eta2);
      TF1 funcPtDepPhi("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      funcPtDepEta.SetParameters(Phi0, Phi1, Phi2);

      // No matched tracks (trackless case)
      if (tracksofcluster.size() == 0) {
        // Use original cluster energy values, no subtraction needed.
        registry.fill(HIST("h2_ClsEvsNmatches"), cluster.energy(), 0);
        registry.fill(HIST("h_Ecluster1"), Ecluster1);
        registry.fill(HIST("h_Ecluster2"), Ecluster2);
        registry.fill(HIST("h_EclusterAll1"), EclusterAll1);
        registry.fill(HIST("h_EclusterAll2"), EclusterAll2);
        registry.fill(HIST("h2_ClsEvsEcluster1"), cluster.energy(), Ecluster1);
        registry.fill(HIST("h2_ClsEvsEcluster2"), cluster.energy(), Ecluster2);
        registry.fill(HIST("h2_ClsEvsEclusterAll1"), cluster.energy(), EclusterAll1);
        registry.fill(HIST("h2_ClsEvsEclusterAll2"), cluster.energy(), EclusterAll2);
        hadroniccorrectedclusters(Ecluster1, Ecluster2, EclusterAll1, EclusterAll2);
        continue;
      }

      // Looping over all matched tracks for the cluster
      // Total number of matched tracks = 20 (hard-coded)
      for (const auto& match : tracksofcluster) {

        double mom = abs(match.track_as<myTracks>().p());
        registry.fill(HIST("h_matchedtracks"), 1);

        // CASE 1: skip tracks with a very low pT
        if (mom < 1e-6) {
          continue;
        } // end CASE 1

        // CASE 2:
        //  a) If one matched track -> 100% energy subtraction
        //  b) If more than one matched track -> 100% energy subtraction
        //  c) If you want to do systematic studies -> perform the above two checks a) and b), and then subtract 70% energy instead of 100%

        // Perform dEta/dPhi matching
        double dEta = match.track_as<myTracks>().trackEtaEmcal() - cluster.eta();
        double dPhi = TVector2::Phi_mpi_pi(match.track_as<myTracks>().trackPhiEmcal() - cluster.phi());

        // Apply the eta and phi matching thresholds
        // dEta and dPhi cut : ensures that the matched track is within the desired eta/phi window

        // Do pT-dependent track matching
        if (doMomDepMatching) {
          auto trackEtaMax = funcPtDepEta.Eval(mom);
          auto trackPhiHigh = +funcPtDepPhi.Eval(mom);
          auto trackPhiLow = -funcPtDepPhi.Eval(mom);

          if ((dPhi < trackPhiHigh && dPhi > trackPhiLow) && fabs(dEta) < trackEtaMax) {
            if (Nmatches == 0) {
              closestTrkP = mom;
            }
            totalTrkP += mom;
            Nmatches++;
          }
        } else {
          // Do fixed dEta/dPhi matching (non-pT dependent)
          if (fabs(dEta) >= minDEta || fabs(dPhi) >= minDPhi) {
            continue; // Skip this track if outside the fixed cut region
          }

          // If track passes fixed dEta/dPhi cuts, process it
          if (Nmatches == 0) {
            closestTrkP = mom; // Closest track match
          }
          totalTrkP += mom; // Accumulate momentum
          Nmatches++;       // Count this match
        }

      } // End of track loop
      registry.fill(HIST("h2_ClsEvsNmatches"), cluster.energy(), Nmatches);

      if (Nmatches >= 1) {
        if (UseM02SubtractionScheme1) {
          // Do M02-based correction if enabled
          Ecluster1 = subtractM02ClusterEnergy(cluster.m02(), Ecluster1, Nmatches, closestTrkP, fHadCorr1, UseFraction1);
          EclusterAll1 = subtractM02ClusterEnergy(cluster.m02(), EclusterAll1, Nmatches, totalTrkP, fHadCorralltrks1, UseFraction1);
        } else {
          // Default energy subtraction (100% and 70%)
          Ecluster1 = subtractClusterEnergy(Ecluster1, closestTrkP, fHadCorr1, Nmatches, UseFraction1);
          EclusterAll1 = subtractClusterEnergy(EclusterAll1, totalTrkP, fHadCorralltrks1, Nmatches, UseFraction1);
        }

        if (UseM02SubtractionScheme2) {
          // Do M02-based correction if enabled
          Ecluster2 = subtractM02ClusterEnergy(cluster.m02(), Ecluster2, Nmatches, closestTrkP, fHadCorr2, UseFraction2);
          EclusterAll2 = subtractM02ClusterEnergy(cluster.m02(), EclusterAll2, Nmatches, totalTrkP, fHadCorralltrks2, UseFraction2);
        } else {
          // Default energy subtraction (100% and 70%)
          Ecluster2 = subtractClusterEnergy(Ecluster2, closestTrkP, fHadCorr2, Nmatches, UseFraction2);
          EclusterAll2 = subtractClusterEnergy(EclusterAll2, totalTrkP, fHadCorralltrks2, Nmatches, UseFraction2);
        }
      }
      registry.fill(HIST("h_Ecluster1"), Ecluster1);
      registry.fill(HIST("h_Ecluster2"), Ecluster2);
      registry.fill(HIST("h_EclusterAll1"), EclusterAll1);
      registry.fill(HIST("h_EclusterAll2"), EclusterAll2);
      registry.fill(HIST("h2_ClsEvsEcluster1"), cluster.energy(), Ecluster1);
      registry.fill(HIST("h2_ClsEvsEcluster2"), cluster.energy(), Ecluster2);
      registry.fill(HIST("h2_ClsEvsEclusterAll1"), cluster.energy(), EclusterAll1);
      registry.fill(HIST("h2_ClsEvsEclusterAll2"), cluster.energy(), EclusterAll2);

      // Fill the table with all four corrected energies
      hadroniccorrectedclusters(Ecluster1, Ecluster2, EclusterAll1, EclusterAll2);

    } // End of cluster loop
  } // process function ends

  // Helper function to prevent negative energy subtraction
  double subtractClusterEnergy(double Ecluster, double mom, double correctionFactor, int Nmatches, bool UseFraction)
  {
    double Ecorr = Ecluster;

    // if (UseConstantSubtractionValue) {
    if (!UseFraction) {
      Ecorr = Ecluster - fConstantSubtractionValue * Nmatches; // Use constant value for MIP-subtraction (regardless of the cluster-shape)
    } else {
      Ecorr = Ecluster - correctionFactor * mom; // Fractional momentum subtraction
    }
    return (Ecorr < 0) ? 0 : Ecorr;
  }

  // Helper function for M02-based energy subtraction (based on cluster-shape)
  double subtractM02ClusterEnergy(double m02, double Ecluster, int Nmatches, double totalTrkP, double correctionFactor, bool UseFraction)
  {
    double Ecorr = Ecluster;

    // For M02 in the single photon region, the signal is primarily: Single photons, single electrons, single MIPs
    if (m02 > 0.1 && m02 < 0.4) { // circular clusters(electron/photon)
      Ecorr = 0;           // Single electron, single MIP
    } else if (m02 > 0.4) {
      // Large M02 region (M02 > 0.4), more complex overlaps and hadronic showers.
      // The signal is primarily: Single hadronic shower, photon-photon overlap, photon-MIP overlap, MIP-MIP overlap,
      // MIP-hadronic shower overlap, hadronic shower - hadronic shower overlap)

      if (!UseFraction) {
        Ecorr = Ecluster - fConstantSubtractionValue * Nmatches; // Use constant value for MIP-subtraction (regardless of the cluster-shape)
      } else {
        Ecorr = Ecluster - correctionFactor * totalTrkP; // Fractional momentum subtraction
      }
    }
    return (Ecorr < 0) ? 0 : Ecorr; // Prevent negative energy
  }

  PROCESS_SWITCH(EmcalClusterHadronicCorrectionTask, processMatchedCollisions, "Process matched clusters from collision", true);
}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EmcalClusterHadronicCorrectionTask>(cfgc, TaskName{"emcal-cluster-hadronic-correction-task"})}; }
