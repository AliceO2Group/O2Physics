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
#include <set>
#include <utility>

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
using collisionEvSelItMC = o2::aod::McCollisions;
using MCClusters = o2::soa::Join<o2::aod::EMCALMCClusters, o2::aod::EMCALClusters>;

struct PhotonIsolationQA {
  HistogramRegistry Data_Info{"Data_Info", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_Info{"MC_Info", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  using o2HistType = HistType;
  using o2Axis = AxisSpec;

  o2::emcal::Geometry* mGeometry = nullptr;

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
  Configurable<float> Track_matching_Radius{"Track_matching_Radius", 0.05, "Radius for which a high energetic track is matched to a cluster"};
  Configurable<bool> isMC{"isMC", true, "should be set to true if the data set is monte carlo"};

  Filter PosZFilter = nabs(aod::collision::posZ) < maxPosZ;
  Filter PosZFilterMC = nabs(aod::mccollision::posZ) < maxPosZ;
  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == mClusterDefinition) && (o2::aod::emcalcluster::time >= minTime) && (o2::aod::emcalcluster::time <= maxTime) && (o2::aod::emcalcluster::energy > minClusterEnergy) && (o2::aod::emcalcluster::nCells >= minNCells) && (o2::aod::emcalcluster::nlm <= maxNLM) && (o2::aod::emcalcluster::isExotic == ExoticContribution);
  Filter emccellfilter = aod::calo::caloType == 1;

  using selectedCollisions = soa::Filtered<collisionEvSelIt>;
  using selectedMcCollisions = soa::Filtered<collisionEvSelItMC>;
  using selectedClusters = soa::Filtered<aod::EMCALClusters>;
  using selectedMCClusters = soa::Filtered<MCClusters>;

  // Preslices
  Preslice<aod::Collisions> collisionsPerBC = aod::collision::bcId;
  Preslice<aod::McCollisions> McCollisionsPerBC = aod::mccollision::bcId;
  Preslice<aod::Tracks> TracksPercollision = aod::track::collisionId;
  Preslice<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<selectedMCClusters> ClustersPerCol = aod::emcalcluster::collisionId;
  Preslice<aod::EMCALClusterCells> CellsPerCluster = aod::emcalclustercell::emcalclusterId;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(300000);
    auto nCells = mGeometry->GetNCells();

    const o2Axis PosZ_Axis{100, -15, 15, "Z Positions (cm)"};
    const o2Axis Num_Cluster_Axis{10, 0, 10, "N Clusters"};
    const o2Axis Shower_Shape_Long_Axis{100, 0, 4, "#sigma^{2}_{long}"};
    const o2Axis Shower_Shape_Short_Axis{100, 0, 4, "#sigma^{2}_{short}"};
    const o2Axis Phi_Axis{100, 1, 6, "#phi"};
    const o2Axis Eta_Axis{100, -0.9, 0.9, "#eta"};
    const o2Axis Energy_Axis{100, 0, 100, "E (GeV/c)"};
    const o2Axis NCells_Axis{50, 0, 50, "N Cells"};
    const o2Axis Num_Track_Axis{20, -0.5, 19.5, "N Tracks"};
    const o2Axis PtIso_Axis{100, -10, 15, "P_{T, Iso} (GeV/c)"};
    const o2Axis Rho_Axis{200, 0, 200, "#rho (#frac{GeV/c}{A})"};
    const o2Axis ABCD_Axis{5, 0, 5, "ABCD"};
    const o2Axis NLM_Axis{50, -0.5, 49.5, "NLM"};
    const o2Axis PDG_Axis{20000, -10000.5, 9999.5, "PDG Code"};
    const o2Axis Status_Code_Axis{400, -200.5, 199.5, "Status Code"};
    const o2Axis BC_Axis{100, -0.5, 99.5, "Col per BC"};
    const o2Axis CellNumber_Axis{nCells, -0.5, nCells - 0.5, "Cell Number"};
    const o2Axis SM_Flag_Axis{2, -0.5, 1.5};

    Data_Info.add("hPosZ", "Z Position of collision", o2HistType::kTH1F, {PosZ_Axis});
    Data_Info.add("hNumClusters", "Number of cluster per collision", o2HistType::kTH1F, {Num_Cluster_Axis});
    Data_Info.add("hClusterLocation", "Location of shower in eta phi plane", o2HistType::kTH2F, {{Eta_Axis}, {Phi_Axis}});
    Data_Info.add("hEnergy_ShowerShapeLong", "Energy vs Shower shape long axis", o2HistType::kTH2F, {{Energy_Axis}, {Shower_Shape_Long_Axis}});
    Data_Info.add("hEnergy_ShowerShapeShort", "Energy vs Shower shape short axis", o2HistType::kTH2F, {{Energy_Axis}, {Shower_Shape_Short_Axis}});
    Data_Info.add("hEnergy_m02_m20", "Energy cluster Vs m02 vs m20", o2HistType::kTH3F, {{Energy_Axis}, {Shower_Shape_Long_Axis}, {Shower_Shape_Short_Axis}});
    Data_Info.add("hEnergy_NCells", "Energy vs Number of cells in cluster", o2HistType::kTH2F, {{Energy_Axis}, {NCells_Axis}});
    Data_Info.add("hEvsNumTracks", "Energy of cluster vs matched tracks", o2HistType::kTH2F, {{Energy_Axis}, {Num_Track_Axis}});

    Data_Info.add("hEvsPtIso", "Pt_Iso", o2HistType::kTH2F, {{Energy_Axis}, {PtIso_Axis}});
    Data_Info.add("hRho_Perpen_Cone", "Energy vs Density of perpendicular cone", o2HistType::kTH2F, {{Energy_Axis}, {Rho_Axis}});
    Data_Info.add("hShowerShape", "Shower shape", o2HistType::kTH2F, {{Shower_Shape_Long_Axis}, {Shower_Shape_Short_Axis}});
    Data_Info.add("hSigmaLongvsPtIso", "Long shower shape vs Pt_Iso", o2HistType::kTH3F, {{Shower_Shape_Long_Axis}, {PtIso_Axis}, {Energy_Axis}});
    Data_Info.add("hABCDControlRegion", "Yield Control Regions", o2HistType::kTH2F, {{ABCD_Axis}, {Energy_Axis}});
    Data_Info.add("hCollperBC", "collisions per BC", o2HistType::kTH1F, {BC_Axis});
    Data_Info.add("hEnergy_NLM_Flag", "Energy vs NLM", o2HistType::kTH3F, {{Energy_Axis}, {NLM_Axis}, {SM_Flag_Axis}});
    Data_Info.add("hNCells_NLM_Flag", "Energy vs NLM", o2HistType::kTH3F, {{NCells_Axis}, {NLM_Axis}, {SM_Flag_Axis}});

    MC_Info.add("hPosZ", "Z Position of collision", o2HistType::kTH1F, {PosZ_Axis});
    MC_Info.get<TH1>(HIST("hPosZ"))->Sumw2();
    MC_Info.add("hNumClusters", "Number of cluster per collision", o2HistType::kTH1F, {Num_Cluster_Axis});
    MC_Info.get<TH1>(HIST("hNumClusters"))->Sumw2();
    MC_Info.add("hClusterLocation", "Location of shower in eta phi plane", o2HistType::kTH2F, {{Eta_Axis}, {Phi_Axis}});
    MC_Info.add("hEnergy_ShowerShapeLong", "Energy vs Shower shape long axis", o2HistType::kTH2F, {{Energy_Axis}, {Shower_Shape_Long_Axis}});
    MC_Info.get<TH2>(HIST("hEnergy_ShowerShapeLong"))->Sumw2();
    MC_Info.add("hEnergy_ShowerShapeShort", "Energy vs Shower shape short axis", o2HistType::kTH2F, {{Energy_Axis}, {Shower_Shape_Short_Axis}});
    MC_Info.get<TH2>(HIST("hEnergy_ShowerShapeShort"))->Sumw2();
    MC_Info.add("hEnergy_m02_m20", "Energy cluster Vs m02 vs m20", o2HistType::kTH3F, {{Energy_Axis}, {Shower_Shape_Long_Axis}, {Shower_Shape_Short_Axis}});
    MC_Info.get<TH3>(HIST("hEnergy_m02_m20"))->Sumw2();
    MC_Info.add("hEnergy_NCells", "Energy vs Number of cells in cluster", o2HistType::kTH2F, {{Energy_Axis}, {NCells_Axis}});
    MC_Info.get<TH2>(HIST("hEnergy_NCells"))->Sumw2();

    MC_Info.add("hEvsNumTracks", "Energy of cluster vs matched tracks", o2HistType::kTH2F, {{Energy_Axis}, {Num_Track_Axis}});
    MC_Info.get<TH2>(HIST("hEvsNumTracks"))->Sumw2();
    MC_Info.add("hEvsPtIso", "Pt_Iso", o2HistType::kTH2F, {{Energy_Axis}, {PtIso_Axis}});
    MC_Info.get<TH2>(HIST("hEvsPtIso"))->Sumw2();
    MC_Info.add("hRho_Perpen_Cone", "Energy vs Density of perpendicular cone", o2HistType::kTH2F, {{Energy_Axis}, {Rho_Axis}});
    MC_Info.get<TH2>(HIST("hRho_Perpen_Cone"))->Sumw2();
    MC_Info.add("hShowerShape", "Shower shape", o2HistType::kTH2F, {{Shower_Shape_Long_Axis}, {Shower_Shape_Short_Axis}});
    MC_Info.get<TH2>(HIST("hShowerShape"))->Sumw2();
    MC_Info.add("hSigmaLongvsPtIso", "Long shower shape vs Pt_Iso", o2HistType::kTH3F, {{Shower_Shape_Long_Axis}, {PtIso_Axis}, {Energy_Axis}});
    MC_Info.get<TH3>(HIST("hSigmaLongvsPtIso"))->Sumw2();
    MC_Info.add("hABCDControlRegion", "Yield Control Regions", o2HistType::kTH2F, {{ABCD_Axis}, {Energy_Axis}});
    MC_Info.get<TH2>(HIST("hABCDControlRegion"))->Sumw2();
    MC_Info.add("hClusterEnergy_MCParticleEnergy", "Energy cluster vs energy particle of cluster", o2HistType::kTH2F, {{Energy_Axis}, {Energy_Axis}});
    MC_Info.get<TH2>(HIST("hClusterEnergy_MCParticleEnergy"))->Sumw2();
    MC_Info.add("hMotherPDG", "PDG code of candidate photons mother", o2HistType::kTH1F, {{2000, -1000.5, 999.5}});
    MC_Info.add("hMotherStatusCode", "Statuscode of candidate photons mother", o2HistType::kTH1F, {{400, -200.5, 199.5}});
    MC_Info.add("hMotherStatusCodeVsPDG", "Statuscode of candidate photons mother", o2HistType::kTH2F, {{Status_Code_Axis}, {PDG_Axis}});
    MC_Info.add("hCollperBC", "collisions per BC", o2HistType::kTH1F, {BC_Axis});
    MC_Info.add("hEnergy_NLM_Flag", "Energy vs NLM", o2HistType::kTH3F, {{Energy_Axis}, {NLM_Axis}, {SM_Flag_Axis}});
    MC_Info.get<TH3>(HIST("hEnergy_NLM_Flag"))->Sumw2();
    MC_Info.add("hNCells_NLM_Flag", "Energy vs NLM", o2HistType::kTH3F, {{NCells_Axis}, {NLM_Axis}, {SM_Flag_Axis}});
    MC_Info.get<TH3>(HIST("hNCells_NLM_Flag"))->Sumw2();
    MC_Info.add("hPromtPhoton", "Energy vs m02 vs NCells, PtIso", o2HistType::kTHnSparseF, {{Energy_Axis}, {Shower_Shape_Long_Axis}, {NCells_Axis}, {PtIso_Axis}});

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
      double dphi = cluster.phi() - match.track_as<myGlobTracks>().phi();
      if (abs(dphi) > M_PI) {
        dphi = 2. * M_PI - abs(dphi);
      }
      double distance = sqrt(pow((cluster.eta() - match.track_as<myGlobTracks>().eta()), 2) + pow(dphi, 2));
      if (distance < Track_matching_Radius) {
        double abs_pt = abs(match.track_as<myGlobTracks>().pt());
        if ((cluster.energy() / abs_pt) < 1.75) {
          return true;
        }
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

  void fillclusterhistos(const auto cluster, HistogramRegistry registry, double weight = 1.0)
  {
    registry.fill(HIST("hClusterLocation"), cluster.eta(), cluster.phi());
    if (isMC == true) {
      registry.fill(HIST("hEnergy_ShowerShapeLong"), cluster.energy(), cluster.m02(), weight);
      registry.fill(HIST("hEnergy_ShowerShapeShort"), cluster.energy(), cluster.m20(), weight);
      registry.fill(HIST("hEnergy_NCells"), cluster.energy(), cluster.nCells(), weight);
      registry.fill(HIST("hEnergy_m02_m20"), cluster.energy(), cluster.m02(), cluster.m20(), weight);
      registry.fill(HIST("hShowerShape"), cluster.m02(), cluster.m20(), weight);
    } else {
      registry.fill(HIST("hEnergy_ShowerShapeLong"), cluster.energy(), cluster.m02());
      registry.fill(HIST("hEnergy_ShowerShapeShort"), cluster.energy(), cluster.m20());
      registry.fill(HIST("hEnergy_NCells"), cluster.energy(), cluster.nCells());
      registry.fill(HIST("hEnergy_m02_m20"), cluster.energy(), cluster.m02(), cluster.m20());
      registry.fill(HIST("hShowerShape"), cluster.m02(), cluster.m20());
    }
  }

  void fillABCDHisto(HistogramRegistry registry, const auto& cluster, double Pt_iso, double weight = 1.0)
  {
    if (isMC == true) {
      if ((Pt_iso < 1.5) && (cluster.m02() < 0.3) && (cluster.m02() > 0.1)) {
        registry.fill(HIST("hABCDControlRegion"), 0.5, cluster.energy(), weight);
      }
      if ((Pt_iso > 4.0) && (cluster.m02() < 0.3) && (cluster.m02() > 0.1)) {
        registry.fill(HIST("hABCDControlRegion"), 1.5, cluster.energy(), weight);
      }
      if ((Pt_iso < 1.5) && (cluster.m02() < 2.0) && (cluster.m02() > 0.4)) {
        registry.fill(HIST("hABCDControlRegion"), 2.5, cluster.energy(), weight);
      }
      if ((Pt_iso > 4.0) && (cluster.m02() < 2.0) && (cluster.m02() > 0.4)) {
        registry.fill(HIST("hABCDControlRegion"), 3.5, cluster.energy(), weight);
      }
    } else {
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
  }

  // iterates over all mothers to check if photon originated from hard scattering (statuscode = abs(23))
  template <typename T>
  int getOriginalMotherIndex(const typename T::iterator& particle)
  {
    if (abs(particle.getGenStatusCode()) == 23) {
      return particle.getGenStatusCode();
    }
    auto mother = particle;

    while (mother.has_mothers()) {
      mother = mother.template mothers_first_as<T>();

      MC_Info.fill(HIST("hMotherPDG"), mother.pdgCode());
      MC_Info.fill(HIST("hMotherStatusCode"), mother.getGenStatusCode());
      MC_Info.fill(HIST("hMotherStatusCodeVsPDG"), mother.getGenStatusCode(), mother.pdgCode());

      int motherStatusCode = mother.getGenStatusCode();
      int motherPDGCode = mother.pdgCode();

      if (abs(motherStatusCode) == 23 && motherPDGCode == 22) {
        return motherStatusCode;
      }
    }
    return -1.0;
  }

  // Calculates the number of local maxima within a cluster
  std::pair<int, int> CalculateNLM(const auto& ClusterCells)
  {
    std::vector<std::vector<float>> Cell_Info(ClusterCells.size(), std::vector<float>(3));
    std::vector<int> supermodules;

    int idx = 0;
    for (auto& Cell : ClusterCells) {
      auto [supermodule, module, phiInModule, etaInModule] = mGeometry->GetCellIndex(Cell.calo().cellNumber());
      auto [row, col] = mGeometry->GetCellPhiEtaIndexInSModule(supermodule, module, phiInModule, etaInModule);
      supermodules.push_back(supermodule);
      Cell_Info[idx++] = {static_cast<float>(row), static_cast<float>(col), Cell.calo().amplitude()};
    }

    std::vector<std::vector<float>> local_max_vector(ClusterCells.size(), std::vector<float>(3));
    for (size_t i = 0; i < Cell_Info.size(); ++i) {
      float row = Cell_Info[i][0];
      float col = Cell_Info[i][1];
      float amp = Cell_Info[i][2];

      bool updated;
      do {
        updated = false;
        for (const auto& cell : Cell_Info) {
          float dr = cell[0] - row;
          float dc = cell[1] - col;
          if (std::abs(dr) <= 1 && std::abs(dc) <= 1 && cell[2] > amp) {
            row += dr;
            col += dc;
            amp = cell[2];
            updated = true;
          }
        }
      } while (updated);

      local_max_vector[i] = {row, col, amp};
    }

    std::set<int> uniqueRows;
    for (const auto& localMax : local_max_vector) {
      uniqueRows.insert(localMax[0]);
    }

    int NLM = uniqueRows.size();

    // flag = 0 if cluster falls in 1 supermodule. flag = 1 if cluster falls in multiple supermodules and will have automatically more local maxima
    int flag = (std::unordered_set<int>(supermodules.begin(), supermodules.end()).size() > 1) ? 1 : 0;
    return std::make_pair(NLM, flag);
  }

  // process monte carlo data
  void processMC(aod::BCs const& bcs, selectedMcCollisions const& Collisions, selectedMCClusters const& mcclusters, aod::McParticles const&, myGlobTracks const& tracks, o2::aod::EMCALMatchedTracks const& matchedtracks, aod::Calos const&, aod::EMCALClusterCells const& ClusterCells)
  {
    for (auto bc : bcs) {
      auto collisionsInBC = Collisions.sliceBy(McCollisionsPerBC, bc.globalIndex());
      MC_Info.fill(HIST("hCollperBC"), collisionsInBC.size());
      if (collisionsInBC.size() == 1) {
        for (const auto& Collision : collisionsInBC) {
          MC_Info.fill(HIST("hPosZ"), Collision.posZ(), Collision.weight());
          auto ClustersInCol = mcclusters.sliceBy(ClustersPerCol, Collision.globalIndex());
          auto tracksInCol = tracks.sliceBy(TracksPercollision, Collision.globalIndex());

          if (ClustersInCol.size() > 0) {
            MC_Info.fill(HIST("hNumClusters"), ClustersInCol.size(), Collision.weight());
          }

          for (auto& mccluster : ClustersInCol) {
            auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, mccluster.globalIndex());
            fillclusterhistos(mccluster, MC_Info, Collision.weight());
            MC_Info.fill(HIST("hEvsNumTracks"), mccluster.energy(), tracksofcluster.size(), Collision.weight());

            auto CellsInCluster = ClusterCells.sliceBy(CellsPerCluster, mccluster.globalIndex());
            auto [NLM, flag] = CalculateNLM(CellsInCluster);
            MC_Info.fill(HIST("hEnergy_NLM_Flag"), mccluster.energy(), NLM, flag, Collision.weight());
            MC_Info.fill(HIST("hNCells_NLM_Flag"), mccluster.nCells(), NLM, flag, Collision.weight());

            if (!track_matching(mccluster, tracksofcluster)) { // no track with significant momentum is matched to cluster
              if (NLM <= maxNLM) {
                double Pt_Cone = sum_Pt_tracks_in_cone(mccluster, tracksInCol);
                double Rho_Perpen_Cone = Rho_Perpendicular_Cone(mccluster, tracksInCol);
                double Pt_iso = Pt_Iso(Pt_Cone, Rho_Perpen_Cone);

                MC_Info.fill(HIST("hEvsPtIso"), mccluster.energy(), Pt_iso, Collision.weight());
                MC_Info.fill(HIST("hRho_Perpen_Cone"), mccluster.energy(), Rho_Perpen_Cone, Collision.weight());
                MC_Info.fill(HIST("hSigmaLongvsPtIso"), mccluster.m02(), Pt_iso, mccluster.energy(), Collision.weight());
                fillABCDHisto(MC_Info, mccluster, Pt_iso, Collision.weight());

                // acces mc true info
                auto ClusterParticles = mccluster.mcParticle_as<aod::McParticles>();
                bool background = true;
                for (auto& clusterparticle : ClusterParticles) {
                  if (clusterparticle.pdgCode() == 22) {
                    MC_Info.fill(HIST("hClusterEnergy_MCParticleEnergy"), mccluster.energy(), clusterparticle.e(), Collision.weight());
                    int first_mother_status_code = getOriginalMotherIndex<aod::McParticles>(clusterparticle);
                    if (abs(first_mother_status_code) == 23) {
                      background = false;
                      MC_Info.fill(HIST("hPromtPhoton"), mccluster.energy(), mccluster.m02(), mccluster.nCells(), Pt_iso);
                    }
                  }
                }
                if (background) {
                  if ((Pt_iso < 1.5) && (mccluster.m02() < 0.3) && (mccluster.m02() > 0.1)) {
                    MC_Info.fill(HIST("hABCDControlRegion"), 4.5, mccluster.energy(), Collision.weight());
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(PhotonIsolationQA, processMC, "proces MC data", true);

  void processData(aod::BCs const& bcs, selectedCollisions const& Collisions, selectedClusters const& clusters, o2::aod::EMCALMatchedTracks const& matchedtracks, myGlobTracks const& tracks, aod::Calos const&, aod::EMCALClusterCells const& ClusterCells)
  {
    for (auto bc : bcs) {
      auto collisionsInBC = Collisions.sliceBy(collisionsPerBC, bc.globalIndex());
      Data_Info.fill(HIST("hCollperBC"), collisionsInBC.size());
      if (collisionsInBC.size() == 1) {
        for (const auto& Collision : collisionsInBC) {
          Data_Info.fill(HIST("hPosZ"), Collision.posZ());
          auto ClustersInCol = clusters.sliceBy(ClustersPerCol, Collision.globalIndex());
          auto tracksInCol = tracks.sliceBy(TracksPercollision, Collision.globalIndex());

          if (ClustersInCol.size() > 0) {
            Data_Info.fill(HIST("hNumClusters"), ClustersInCol.size());
          }

          for (auto& cluster : ClustersInCol) {
            auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex());
            fillclusterhistos(cluster, Data_Info);
            Data_Info.fill(HIST("hEvsNumTracks"), cluster.energy(), tracksofcluster.size());

            auto CellsInCluster = ClusterCells.sliceBy(CellsPerCluster, cluster.globalIndex());
            auto [NLM, flag] = CalculateNLM(CellsInCluster);
            Data_Info.fill(HIST("hEnergy_NLM_Flag"), cluster.energy(), NLM, flag);
            Data_Info.fill(HIST("hNCells_NLM_Flag"), cluster.nCells(), NLM, flag);

            if (!track_matching(cluster, tracksofcluster)) { // no track with significant momentum is matched to cluster
              if (NLM < maxNLM) {
                double Pt_Cone = sum_Pt_tracks_in_cone(cluster, tracksInCol);
                double Rho_Perpen_Cone = Rho_Perpendicular_Cone(cluster, tracksInCol);
                double Pt_iso = Pt_Iso(Pt_Cone, Rho_Perpen_Cone);

                Data_Info.fill(HIST("hEvsPtIso"), cluster.energy(), Pt_iso);
                Data_Info.fill(HIST("hRho_Perpen_Cone"), cluster.energy(), Rho_Perpen_Cone);
                Data_Info.fill(HIST("hSigmaLongvsPtIso"), cluster.m02(), Pt_iso, cluster.energy());
                fillABCDHisto(Data_Info, cluster, Pt_iso);
              }
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(PhotonIsolationQA, processData, "proces data", true);
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PhotonIsolationQA>(cfgc, TaskName{"photon-isolation-qa"})};
}
