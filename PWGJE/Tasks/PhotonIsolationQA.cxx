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

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"

#include "CommonDataFormat/InteractionRecord.h"

// \struct PhotonIsolationQA
/// \brief Task to select emcal clusters originating from promt photons
/// \author Berend van Beuzekom <b.d.vanbeuzekom@students.uu.nl>
/// \since 30-05-2024
///
/// This task is designed to select EMCal clusters originating from prompt photons.
/// Cluster selection is performed using Pt_Iso (density of particles around the cluster minus UE density) and cluster shape.
/// The task can also be run over Monte Carlo data to obtain a correction factor for the purity measurement using the ABCD method.

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using myGlobTracks = o2::soa::Join<o2::aod::FullTracks, o2::aod::TrackSelection>;
using collisionEvSelIt = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
using mcCells = o2::soa::Join<aod::Calos, aod::McCaloLabels_001>;
using filteredMCCells = o2::soa::Filtered<mcCells>;
using MCClusters = o2::soa::Join<o2::aod::EMCALMCClusters, o2::aod::EMCALClusters>;

struct PhotonIsolationQA {
  HistogramRegistry Data_Info{"Data_Info", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_Info{"MC_Info", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  using o2HistType = HistType;
  using o2Axis = AxisSpec;

  Configurable<float> maxPosZ{"maxPosZ", 10.0f, "maximum z position of collision in cm"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<int> mClusterDefinition{"clusterDefinition", 0, "cluster definition to be selected, e.g. 10=kV3Default, 0 = kV1Default"};
  Configurable<float> minTime{"minTime", -30., "Minimum cluster time for time cut"};
  Configurable<float> maxTime{"maxTime", +35., "Maximum cluster time for time cut"};
  Configurable<float> minClusterEnergy{"minClusterEnergy", 0.7f, "Minimal cluster energy"};
  Configurable<int> minNCells{"minNCelss", 2, "Minimal amount of cells per cluster"};
  Configurable<int> maxNLM{"maxNLM", 2, "Maximal amount of local maxima per cluster"};
  Configurable<bool> ExoticContribution{"ExoticContribution", false, "Exotic cluster in the data"};
  Configurable<float> minDPhi{"minDPhi", 0.01, "Minimum dPhi between track and cluster"};

  Filter PosZFilter = nabs(aod::collision::posZ) < maxPosZ;
  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == mClusterDefinition) && (o2::aod::emcalcluster::time >= minTime) && (o2::aod::emcalcluster::time <= maxTime) && (o2::aod::emcalcluster::energy > minClusterEnergy) && (o2::aod::emcalcluster::nCells >= minNCells) && (o2::aod::emcalcluster::nlm <= maxNLM) && (o2::aod::emcalcluster::isExotic == ExoticContribution);
  Filter emccellfilter = aod::calo::caloType == 1; // mc emcal cell

  using selectedCollisions = soa::Filtered<collisionEvSelIt>;
  using selectedClusters = soa::Filtered<aod::EMCALClusters>;
  using selectedMCClusters = soa::Filtered<MCClusters>;

  Preslice<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalclustercell::emcalclusterId;

  int eventSelection = -1;
  int trackSelection = -1;
  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    const o2Axis PosZ_Axis{100, -15, 15, "Z Positions (cm)"};
    const o2Axis Num_Cluster_Axis{10, 0, 10, "N Clusters"};
    const o2Axis Shower_Shape_Long_Axis{100, 0, 4, "#sigma^{2}_{long}"};
    const o2Axis Shower_Shape_Short_Axis{100, 0, 4, "#sigma^{2}_{short}"};
    const o2Axis Phi_Axis{100, 1, 6, "#phi"};
    const o2Axis Eta_Axis{100, -0.9, 0.9, "#eta"};
    const o2Axis Energy_Axis{100, 0, 100, "E (GeV/c)"};
    const o2Axis NCells_Axis{50, 0, 50, "N Cells"};
    const o2Axis Pt_Axis{100, 0, 100, "P_{T} (GeV/c)"};
    const o2Axis P_Axis{100, 0, 100, "P (GeV/c)"};
    const o2Axis RatioP_Axis{50, 0, 10, "E/P"};
    const o2Axis Num_Track_Axis{20, 0, 20, "N Tracks"};
    const o2Axis PtIso_Axis{100, -10, 15, "P_{T, Iso} (GeV/c)"};
    const o2Axis Rho_Axis{100, 0, 100, "#rho (#frac{GeV/c}{A})"};
    const o2Axis ABCD_Axis{5, 0, 5, "ABCD"};
    const o2Axis NLM_Axis{3, 0, 3, "NLM"};
    const o2Axis Fraction_Axis{50, 0, 2};

    Data_Info.add("hPosZ", "Z Position of collision", o2HistType::kTH1F, {PosZ_Axis});
    Data_Info.add("hNumClusters", "Number of cluster per collision", o2HistType::kTH1F, {Num_Cluster_Axis});
    Data_Info.add("hClusterLocation", "Location of shower in eta phi plane", o2HistType::kTH2F, {{Eta_Axis}, {Phi_Axis}});
    Data_Info.add("hEnergy", "Energy of the cluster", o2HistType::kTH1F, {Energy_Axis});
    Data_Info.add("hEnergy_ShowerShapeLong", "Energy vs Shower shape long axis", o2HistType::kTH2F, {{Energy_Axis}, {Shower_Shape_Long_Axis}});
    Data_Info.add("hEnergy_ShowerShapeShort", "Energy vs Shower shape short axis", o2HistType::kTH2F, {{Energy_Axis}, {Shower_Shape_Short_Axis}});
    Data_Info.add("hEnergy_m02_m20", "Energy cluster Vs m02 vs m20", o2HistType::kTHnSparseL, {{Energy_Axis}, {Shower_Shape_Long_Axis}, {Shower_Shape_Short_Axis}});
    Data_Info.add("hEnergy_NCells", "Energy vs Number of cells in cluster", o2HistType::kTH2F, {{Energy_Axis}, {NCells_Axis}});

    Data_Info.add("hEvsNumTracks", "Energy of cluster vs matched tracks", o2HistType::kTH2F, {{Energy_Axis}, {Num_Track_Axis}});
    Data_Info.add("hTrackPt", "Pt of matched track", o2HistType::kTH1F, {Pt_Axis});
    Data_Info.add("hTrackEta", "Eta of matched track", o2HistType::kTH1F, {Eta_Axis});
    Data_Info.add("hTrackPhi", "Phi of matched track", o2HistType::kTH1F, {Phi_Axis});
    Data_Info.add("hRatioClusterETrackP", "Ratio between energy of cluster and P of tracks", o2HistType::kTH1F, {RatioP_Axis});
    Data_Info.add("hRatioClusterETrackPrackPt", "Ratio between E and P vs Pt of matched tracks", o2HistType::kTH2F, {{RatioP_Axis}, {Pt_Axis}});

    Data_Info.add("hEvsPtIso", "Pt_Iso", o2HistType::kTH2F, {{Energy_Axis}, {PtIso_Axis}});
    Data_Info.add("hRho_Perpen_Cone", "Density of perpendicular cone", o2HistType::kTH1F, {Rho_Axis});
    Data_Info.add("hShowerShape", "Shower shape", o2HistType::kTH2F, {{Shower_Shape_Long_Axis}, {Shower_Shape_Short_Axis}});
    Data_Info.add("hSigmaLongvsPtIso", "Long shower shape vs Pt_Iso", o2HistType::kTH2F, {{Shower_Shape_Long_Axis}, {PtIso_Axis}});
    Data_Info.add("hABCDControlRegion", "Yield Control Regions", o2HistType::kTH2F, {{ABCD_Axis}, {Energy_Axis}});
    Data_Info.add("hE_M02_M20_NLM_NCells_PtIso", "Energy vs M02 vs M20 vs NLM vs NCells vs PtIso", o2HistType::kTHnSparseL, {{Energy_Axis}, {Shower_Shape_Long_Axis}, {Shower_Shape_Short_Axis}, {NLM_Axis}, {NCells_Axis}, {PtIso_Axis}});

    MC_Info.add("hPosZ", "Z Position of collision", o2HistType::kTH1F, {PosZ_Axis});
    MC_Info.add("hNumClusters", "Number of cluster per collision", o2HistType::kTH1F, {Num_Cluster_Axis});
    MC_Info.add("hClusterLocation", "Location of shower in eta phi plane", o2HistType::kTH2F, {{Eta_Axis}, {Phi_Axis}});
    MC_Info.add("hEnergy", "Energy of the cluster", o2HistType::kTH1F, {Energy_Axis});
    MC_Info.add("hEnergy_ShowerShapeLong", "Energy vs Shower shape long axis", o2HistType::kTH2F, {{Energy_Axis}, {Shower_Shape_Long_Axis}});
    MC_Info.add("hEnergy_ShowerShapeShort", "Energy vs Shower shape short axis", o2HistType::kTH2F, {{Energy_Axis}, {Shower_Shape_Short_Axis}});
    MC_Info.add("hEnergy_m02_m20", "Energy cluster Vs m02 vs m20", o2HistType::kTHnSparseL, {{Energy_Axis}, {Shower_Shape_Long_Axis}, {Shower_Shape_Short_Axis}});
    MC_Info.add("hEnergy_NCells", "Energy vs Number of cells in cluster", o2HistType::kTH2F, {{Energy_Axis}, {NCells_Axis}});

    MC_Info.add("hEvsNumTracks", "Energy of cluster vs matched tracks", o2HistType::kTH2F, {{Energy_Axis}, {Num_Track_Axis}});
    MC_Info.add("hTrackPt", "Pt of matched track", o2HistType::kTH1F, {Pt_Axis});
    MC_Info.add("hTrackEta", "Eta of matched track", o2HistType::kTH1F, {Eta_Axis});
    MC_Info.add("hTrackPhi", "Phi of matched track", o2HistType::kTH1F, {Phi_Axis});
    MC_Info.add("hRatioClusterETrackP", "Ratio between energy of cluster and P of tracks", o2HistType::kTH1F, {RatioP_Axis});
    MC_Info.add("hRatioClusterETrackPrackPt", "Ratio between E and P vs Pt of matched tracks", o2HistType::kTH2F, {{RatioP_Axis}, {Pt_Axis}});
    MC_Info.add("hE_M02_M20_NLM_NCells_PtIso", "Energy vs M02 vs M20 vs NLM vs NCells vs PtIso", o2HistType::kTHnSparseL, {{Energy_Axis}, {Shower_Shape_Long_Axis}, {Shower_Shape_Short_Axis}, {NLM_Axis}, {NCells_Axis}, {PtIso_Axis}});

    MC_Info.add("hEvsPtIso", "Pt_Iso", o2HistType::kTH2F, {{Energy_Axis}, {PtIso_Axis}});
    MC_Info.add("hRho_Perpen_Cone", "Density of perpendicular cone", o2HistType::kTH1F, {Rho_Axis});
    MC_Info.add("hShowerShape", "Shower shape", o2HistType::kTH2F, {{Shower_Shape_Long_Axis}, {Shower_Shape_Short_Axis}});
    MC_Info.add("hSigmaLongvsPtIso", "Long shower shape vs Pt_Iso", o2HistType::kTH2F, {{Shower_Shape_Long_Axis}, {PtIso_Axis}});
    MC_Info.add("hABCDControlRegion", "Yield Control Regions", o2HistType::kTH2F, {{ABCD_Axis}, {Energy_Axis}});
    MC_Info.add("hMcParticlesToCluster", "Energy of particles linked to mc cluster", o2HistType::kTH1F, {{100, 0, 100}});
    MC_Info.add("hMcParticleStatusCode", "generator status code of mc particle linked to mc cluster", o2HistType::kTH1F, {{200, 0, 200}});
    MC_Info.add("hMcParticleProcessCode", "physical process code of mc particle linked to mc cluster", o2HistType::kTH1F, {{200, 0, 200}});

    std::vector<std::string> bin_names = {"A", "B", "C", "D", "True Bckgr A"};
    for (size_t i = 0; i < bin_names.size(); i++) {
      MC_Info.get<TH2>(HIST("hABCDControlRegion"))->GetXaxis()->SetBinLabel(i + 1, bin_names[i].c_str());
      Data_Info.get<TH2>(HIST("hABCDControlRegion"))->GetXaxis()->SetBinLabel(i + 1, bin_names[i].c_str());
    }
  }

  // boolian returns true if a track is matched to the cluster with E_Cluster/P_track < 1.75. Otherwise returns false
  bool track_matching(const auto& cluster, o2::aod::EMCALMatchedTracks const& matched_tracks)
  {
    for (const auto& match : matched_tracks) {
      double abs_pt = abs(match.track_as<myGlobTracks>().pt());
      if ((cluster.energy() / abs_pt) < 1.75) {
        return true;
      }
    }
    return false;
  }

  // sums the pt of all tracks around the cluster within a radius of R = 0.4
  double sum_Pt_tracks_in_cone(const auto& cluster, myGlobTracks const& tracks, double Cone_Radius = 0.4)
  {
    double sum_Pt = 0.;
    for (const auto& track : tracks) {
      double dphi = cluster.phi() - track.phi();
      if (abs(dphi) > M_PI) {
        dphi = 2. * M_PI - abs(dphi);
      }
      double distance = sqrt(pow((cluster.eta() - track.eta()), 2) + pow(dphi, 2));
      if (distance < Cone_Radius) {
        sum_Pt += track.pt();
      }
    }
    return sum_Pt;
  }

  // Calculates the Pt density of tracks of the UE with the perpendicular cone method
  double Rho_Perpendicular_Cone(const auto& cluster, myGlobTracks const& tracks, double Perpendicular_Cone_Radius = 0.4)
  {
    double sum_Pt = 0.;
    double Perpendicular_Phi1 = cluster.phi() + (1. / 2.) * M_PI;
    double Perpendicular_Phi2 = cluster.phi() - (1. / 2.) * M_PI;
    if (Perpendicular_Phi1 > 2. * M_PI) {
      Perpendicular_Phi1 = Perpendicular_Phi1 - 2. * M_PI;
    }
    if (Perpendicular_Phi2 < 0.) {
      Perpendicular_Phi2 = 2. * M_PI + Perpendicular_Phi2;
    }
    for (const auto& track : tracks) {
      double dphi1 = Perpendicular_Phi1 - track.phi();
      double dphi2 = Perpendicular_Phi2 - track.phi();
      if (abs(dphi1) > M_PI) {
        dphi1 = 2. * M_PI - abs(dphi1);
      }
      if (abs(dphi2) > M_PI) {
        dphi2 = 2. * M_PI - abs(dphi2);
      }
      double distance1 = sqrt(pow((cluster.eta() - track.eta()), 2) + pow(dphi1, 2));
      double distance2 = sqrt(pow((cluster.eta() - track.eta()), 2) + pow(dphi2, 2));
      if (distance1 < Perpendicular_Cone_Radius) {
        sum_Pt += track.pt();
      }
      if (distance2 < Perpendicular_Cone_Radius) {
        sum_Pt += track.pt();
      }
    }
    double Rho = sum_Pt / (2. * M_PI * pow(Perpendicular_Cone_Radius, 2));
    return Rho;
  }

  // Calculates the Pt_Isolation. (Check if the photon candidate is isolated)
  double Pt_Iso(double sum_Pt_tracks_in_cone, double rho, double Cone_Radius = 0.4)
  {
    double Pt_Iso = sum_Pt_tracks_in_cone - rho * M_PI * pow(Cone_Radius, 2);
    return Pt_Iso;
  }

  void fillclusterhistos(const auto cluster, HistogramRegistry registry)
  {
    registry.fill(HIST("hClusterLocation"), cluster.eta(), cluster.phi());
    registry.fill(HIST("hEnergy"), cluster.energy());
    registry.fill(HIST("hEnergy_ShowerShapeLong"), cluster.energy(), cluster.m02());
    registry.fill(HIST("hEnergy_ShowerShapeShort"), cluster.energy(), cluster.m20());
    registry.fill(HIST("hEnergy_NCells"), cluster.energy(), cluster.nCells());
    registry.fill(HIST("hEnergy_m02_m20"), cluster.energy(), cluster.m02(), cluster.m20());
  }

  void fillABCDHisto(HistogramRegistry registry, const auto& cluster, double Pt_iso)
  {
    if ((Pt_iso < 1.5) && (cluster.m02() < 0.3) && (cluster.m02() > 0.1)) {
      registry.fill(HIST("hABCDControlRegion"), 0.5, cluster.energy());
    }
    if ((Pt_iso > 4.0) && (cluster.m02() < 0.3) && (cluster.m02() > 0.1)) {
      registry.fill(HIST("hABCDControlRegion"), 1.5, cluster.energy());
    }
    if ((Pt_iso < 1.5) && (cluster.m02() < 2.0) && (cluster.m02() > 0.4)) {
      registry.fill(HIST("hABCDControlRegion"), 2.5, cluster.energy());
    }
    if ((Pt_iso > 4.0) && (cluster.m02() < 2.0) && (cluster.m02() > 0.4)) {
      registry.fill(HIST("hABCDControlRegion"), 3.5, cluster.energy());
    }
  }

  void fillMatchedTrackHistos(HistogramRegistry registry, const auto& cluster, o2::aod::EMCALMatchedTracks const& matched_tracks)
  {
    for (const auto& match : matched_tracks) {
      registry.fill(HIST("hTrackPt"), match.track_as<myGlobTracks>().pt());
      registry.fill(HIST("hTrackEta"), match.track_as<myGlobTracks>().eta());
      registry.fill(HIST("hTrackPhi"), match.track_as<myGlobTracks>().phi());
      registry.fill(HIST("hRatioClusterETrackP"), cluster.energy() / abs(match.track_as<myGlobTracks>().p()));
      registry.fill(HIST("hRatioClusterETrackPrackPt"), cluster.energy() / abs(match.track_as<myGlobTracks>().p()), match.track_as<myGlobTracks>().pt());
    }
  }

  // process monte carlo data
  void processMC(selectedCollisions::iterator const& theCollision, selectedMCClusters const& mcclusters, aod::StoredMcParticles_001 const&, myGlobTracks const& tracks, o2::aod::EMCALClusterCells const& /*emccluscells*/, o2::aod::EMCALMatchedTracks const& matchedtracks)
  {
    MC_Info.fill(HIST("hPosZ"), theCollision.posZ());

    if (mcclusters.size() > 0) {
      MC_Info.fill(HIST("hNumClusters"), mcclusters.size());
    }

    for (auto& mccluster : mcclusters) {
      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, mccluster.globalIndex());
      fillclusterhistos(mccluster, MC_Info);
      fillMatchedTrackHistos(MC_Info, mccluster, tracksofcluster);
      MC_Info.fill(HIST("hEvsNumTracks"), mccluster.energy(), tracksofcluster.size());

      if (!track_matching(mccluster, tracksofcluster)) { // no track with significant momentum is matched to cluster
        double Pt_Cone = sum_Pt_tracks_in_cone(mccluster, tracks);
        double Rho_Perpen_Cone = Rho_Perpendicular_Cone(mccluster, tracks);
        double Pt_iso = Pt_Iso(Pt_Cone, Rho_Perpen_Cone);

        MC_Info.fill(HIST("hEvsPtIso"), mccluster.energy(), Pt_iso);
        MC_Info.fill(HIST("hRho_Perpen_Cone"), Rho_Perpen_Cone);
        MC_Info.fill(HIST("hShowerShape"), mccluster.m02(), mccluster.m20());
        MC_Info.fill(HIST("hSigmaLongvsPtIso"), mccluster.m02(), Pt_iso);
        MC_Info.fill(HIST("hE_M02_M20_NLM_NCells_PtIso"), mccluster.energy(), mccluster.m02(), mccluster.m20(), mccluster.nlm(), mccluster.nCells(), Pt_iso);
        fillABCDHisto(MC_Info, mccluster, Pt_iso);
        // acces mc true info
        auto ClusterParticles = mccluster.mcParticle_as<aod::StoredMcParticles_001>();
        bool background = true;
        for (auto& clusterparticle : ClusterParticles) {
          if (clusterparticle.pdgCode() == 22) {
            MC_Info.fill(HIST("hMcParticlesToCluster"), clusterparticle.e());
            MC_Info.fill(HIST("hMcParticleProcessCode"), clusterparticle.getProcess());
            MC_Info.fill(HIST("hMcParticleStatusCode"), clusterparticle.getGenStatusCode());
            if ((clusterparticle.getProcess() >= 21) && (clusterparticle.getProcess() <= 29))
              background = false; // Particles from the hardest subprocess
          }
        }

        if (background) {
          if ((Pt_iso < 1.5) && (mccluster.m02() < 0.3) && (mccluster.m02() > 0.1)) {
            MC_Info.fill(HIST("hABCDControlRegion"), 4.5, mccluster.energy());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(PhotonIsolationQA, processMC, "proces MC data", true);

  void processData(selectedCollisions::iterator const& theCollision, selectedClusters const& clusters, o2::aod::EMCALClusterCells const& /*emccluscells*/, o2::aod::EMCALMatchedTracks const& matchedtracks, myGlobTracks const& alltracks)
  {
    Data_Info.fill(HIST("hPosZ"), theCollision.posZ());

    if (clusters.size() > 0) {
      Data_Info.fill(HIST("hNumClusters"), clusters.size());
    }

    for (const auto& cluster : clusters) {
      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex());
      fillclusterhistos(cluster, Data_Info);
      fillMatchedTrackHistos(Data_Info, cluster, tracksofcluster);
      Data_Info.fill(HIST("hEvsNumTracks"), cluster.energy(), tracksofcluster.size());

      if (!track_matching(cluster, tracksofcluster)) { // no track with significant momentum is matched to cluster
        double Pt_Cone = sum_Pt_tracks_in_cone(cluster, alltracks);
        double Rho_Perpen_Cone = Rho_Perpendicular_Cone(cluster, alltracks);
        double Pt_iso = Pt_Iso(Pt_Cone, Rho_Perpen_Cone);
        Data_Info.fill(HIST("hEvsPtIso"), cluster.energy(), Pt_iso);
        Data_Info.fill(HIST("hRho_Perpen_Cone"), Rho_Perpen_Cone);
        Data_Info.fill(HIST("hShowerShape"), cluster.m02(), cluster.m20());
        Data_Info.fill(HIST("hSigmaLongvsPtIso"), cluster.m02(), Pt_iso);
        Data_Info.fill(HIST("hE_M02_M20_NLM_NCells_PtIso"), cluster.energy(), cluster.m02(), cluster.m20(), cluster.nlm(), cluster.nCells(), Pt_iso);
        fillABCDHisto(Data_Info, cluster, Pt_iso);
      }
    }
  }

  PROCESS_SWITCH(PhotonIsolationQA, processData, "proces data", true);
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PhotonIsolationQA>(cfgc, TaskName{"photon-isolation-qa"})};
}
