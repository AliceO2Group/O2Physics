// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
///
/// \file emcalBcWisePi0.cxx
///
/// \brief Task that extracts pi0s from BC wise derived data of EMCal clusters
///
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) Goethe University Frankfurt
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "TString.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/LorentzRotation.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"

#include "CommonConstants/MathConstants.h"
#include "EMCALBase/Geometry.h"
#include "PWGEM/PhotonMeson/DataModel/bcWiseTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct EmcalBcWisePi0 {
  HistogramRegistry mHistManager{"EmcalGammaGammaBcWiseHistograms"};

  Configurable<int> cfgClusterDefinition{"cfgClusterDefinition", 10, "Clusterizer to be selected, e.g. 13 for kV3MostSplitLowSeed"};
  Configurable<float> cfgMinClusterEnergy{"cfgMinClusterEnergy", 0.7, "Minimum energy of selected clusters (GeV)"};
  Configurable<float> cfgMinM02{"cfgMinM02", 0.1, "Minimum M02 of selected clusters"};
  Configurable<float> cfgMaxM02{"cfgMaxM02", 0.7, "Maximum M02 of selected clusters"};
  Configurable<float> cfgMinTime{"cfgMinTime", -15, "Minimum time of selected clusters (ns)"};
  Configurable<float> cfgMaxTime{"cfgMaxTime", 15, "Maximum time of selected clusters (ns)"};
  Configurable<float> cfgRapidityCut{"cfgRapidityCut", 0.8f, "Maximum absolute rapidity of counted particles"};
  Configurable<float> cfgMinOpenAngle{"cfgMinOpenAngle", 0.0202, "Minimum opening angle between photons"};
  Configurable<int> cfgDistanceToEdge{"cfgDistanceToEdge", 1, "Distance to edge in cells required for rotated cluster to be accepted"};
  Configurable<int> cfgBGEventDownsampling{"cfgBGEventDownsampling", 1, "Only calculate background for every n-th event (performance reasons in PbPb)"};

  static constexpr float DefaultCentralityWindow[2] = {-1., 101.};
  Configurable<LabeledArray<float>> cfgCentralityWindow{"cfgCentralityWindow", {DefaultCentralityWindow, 2, {"Min_Centrality", "Max_Centrality"}}, "Select centrality window (also requires unique collision)"};

  Filter clusterDefinitionFilter = aod::bcwisecluster::storedDefinition == static_cast<uint8_t>(cfgClusterDefinition);
  Filter energyFilter = aod::bcwisecluster::storedE > static_cast<uint16_t>(cfgMinClusterEnergy* aod::emdownscaling::downscalingFactors[aod::emdownscaling::kEnergy]);
  Filter m02Filter = (aod::bcwisecluster::storedNCells == static_cast<uint8_t>(1) || (aod::bcwisecluster::storedM02 > static_cast<uint16_t>(cfgMinM02 * aod::emdownscaling::downscalingFactors[aod::emdownscaling::kM02]) && aod::bcwisecluster::storedM02 < static_cast<uint16_t>(cfgMaxM02 * aod::emdownscaling::downscalingFactors[aod::emdownscaling::kM02])));
  Filter timeFilter = (aod::bcwisecluster::storedTime > static_cast<int16_t>(cfgMinTime * aod::emdownscaling::downscalingFactors[aod::emdownscaling::kTime]) && aod::bcwisecluster::storedTime < static_cast<int16_t>(cfgMaxTime * aod::emdownscaling::downscalingFactors[aod::emdownscaling::kTime]));

  emcal::Geometry* emcalGeom;

  void init(InitContext const&)
  {
    emcalGeom = emcal::Geometry::GetInstanceFromRunNumber(300000);
    const int nEventBins = 8;
    mHistManager.add("nBCs", "Number of BCs;;#bf{#it{N}_{BC}}", HistType::kTH1F, {{nEventBins, -0.5, 7.5}});
    const TString binLabels[nEventBins] = {"All", "FT0", "TVX", "kTVXinEMC", "Cell", "Cluster", "NoBorder", "Collision"};
    for (int iBin = 0; iBin < nEventBins; iBin++)
      mHistManager.get<TH1>(HIST("nBCs"))->GetXaxis()->SetBinLabel(iBin + 1, binLabels[iBin]);

    mHistManager.add("Centrality", "FT0M centrality;FT0M centrality (%);#bf{#it{N}_{BC} (only BCs containing exactly 1 collision)}", HistType::kTH1F, {{100, 0., 100.}});

    mHistManager.add("clusterE", "Energy of cluster;#bf{#it{E} (GeV)};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, 0, 20}});
    mHistManager.add("clusterM02", "Shape of cluster;#bf{#it{M}_{02}};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, 0, 2}});
    mHistManager.add("clusterTime", "Time of cluster;#bf{#it{t} (ns)};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, -100, 100}});
    mHistManager.add("clusterNCells", "Number of cells per cluster;#bf{#it{N}_{cells}};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{51, 0., 50.5}});
    mHistManager.add("clusterExotic", "Is cluster exotic?;#bf{Exotic?};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{2, -0.5, 1.5}});
    mHistManager.add("clusterEtaPhi", "Eta/Phi distribution of clusters;#eta;#phi", HistType::kTH2F, {{400, -0.8, 0.8}, {400, 0, constants::math::TwoPI}});

    mHistManager.add("invMassVsPt", "Invariant mass and pT of meson candidates", HistType::kTH2F, {{400, 0., 0.8}, {200, 0., 20.}});
    mHistManager.add("invMassVsPtBackground", "Invariant mass and pT of background meson candidates", HistType::kTH2F, {{400, 0., 0.8}, {200, 0., 20.}});
  }

  /// \brief returns if cluster is too close to edge of EMCal
  bool isTooCloseToEdge(const int cellID, const int DistanceToBorder = 1)
  {
    if (DistanceToBorder <= 0)
      return false;
    if (cellID < 0)
      return true;

    // check distance to border in case the cell is okay
    auto [iSupMod, iMod, iPhi, iEta] = emcalGeom->GetCellIndex(cellID);
    auto [irow, icol] = emcalGeom->GetCellPhiEtaIndexInSModule(iSupMod, iMod, iPhi, iEta);
    int iRowLast = (emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::EMCAL_THIRD || emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::DCAL_EXT) ? 8 : 24;

    return (irow < DistanceToBorder || (iRowLast - irow) <= DistanceToBorder);
  }

  void fillEventHists(const auto& bc, const auto& collisions, const auto& clusters)
  {
    mHistManager.fill(HIST("nBCs"), 0);
    if (bc.hasFT0())
      mHistManager.fill(HIST("nBCs"), 1);
    if (bc.hasTVX())
      mHistManager.fill(HIST("nBCs"), 2);
    if (bc.haskTVXinEMC())
      mHistManager.fill(HIST("nBCs"), 3);
    if (bc.hasEMCCell())
      mHistManager.fill(HIST("nBCs"), 4);
    if (clusters.size() > 0)
      mHistManager.fill(HIST("nBCs"), 5);
    if (bc.hasNoTFROFBorder())
      mHistManager.fill(HIST("nBCs"), 6);
    if (collisions.size() > 0)
      mHistManager.fill(HIST("nBCs"), 7);

    if (collisions.size() == 1) {
      for (const auto& collision : collisions)
        mHistManager.fill(HIST("Centrality"), collision.centrality());
    }
  }

  void fillClusterHists(const auto& clusters)
  {
    for (const auto& cluster : clusters) {
      mHistManager.fill(HIST("clusterE"), cluster.e());
      mHistManager.fill(HIST("clusterM02"), cluster.m02());
      mHistManager.fill(HIST("clusterTime"), cluster.time());
      mHistManager.fill(HIST("clusterNCells"), cluster.nCells());
      mHistManager.fill(HIST("clusterEtaPhi"), cluster.eta(), cluster.phi());
      mHistManager.fill(HIST("clusterExotic"), cluster.isExotic());
    }
  }

  void reconstructMesons(const auto& clusters, int bcId)
  {
    for (const auto& [g1, g2] : soa::combinations(soa::CombinationsStrictlyUpperIndexPolicy(clusters, clusters))) {
      ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      if (std::fabs(v12.Rapidity()) > cfgRapidityCut)
        continue;

      float openingAngle12 = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));
      if (openingAngle12 < cfgMinOpenAngle)
        continue;

      mHistManager.fill(HIST("invMassVsPt"), v12.M(), v12.Pt());

      if (clusters.size() < 3)
        continue;

      if (bcId % cfgBGEventDownsampling != 0)
        continue;

      // "else: Calculate background"

      ROOT::Math::AxisAngle rotationAxis(v12.Vect(), constants::math::PIHalf);
      ROOT::Math::Rotation3D rotationMatrix(rotationAxis);
      for (ROOT::Math::PtEtaPhiMVector vi : {v1, v2}) {

        vi = rotationMatrix * vi;

        try {
          int iCellID = emcalGeom->GetAbsCellIdFromEtaPhi(vi.Eta(), vi.Phi());
          if (isTooCloseToEdge(iCellID, cfgDistanceToEdge))
            continue;
        } catch (o2::emcal::InvalidPositionException& e) {
          continue;
        }

        for (const auto& g3 : clusters) {
          if (g3.globalIndex() == g1.globalIndex() || g3.globalIndex() == g2.globalIndex())
            continue;

          ROOT::Math::PtEtaPhiMVector v3(g3.pt(), g3.eta(), g3.phi(), 0.);

          float openingAnglei3 = std::acos(vi.Vect().Dot(v3.Vect()) / (vi.P() * v3.P()));
          if (openingAnglei3 < cfgMinOpenAngle)
            continue;

          ROOT::Math::PtEtaPhiMVector vBG = v3 + vi;

          mHistManager.fill(HIST("invMassVsPtBackground"), vBG.M(), vBG.Pt());
        }
      }
    }
  }

  bool isCentralitySelected(const auto& collisions)
  {
    if (cfgCentralityWindow->get("Min_Centrality") > -0.5 || cfgCentralityWindow->get("Max_Centrality") < 100.5) { // Centrality window is set
      if (collisions.size() != 1)
        return false;
      for (const auto& collision : collisions) {
        if (collision.centrality() < cfgCentralityWindow->get("Min_Centrality") || collision.centrality() > cfgCentralityWindow->get("Max_Centrality"))
          return false;
      }
    }
    return true;
  }

  void process(aod::BCWiseBCs::iterator const& bc, aod::BCWiseCollisions const& collisions, soa::Filtered<aod::BCWiseClusters> const& clusters)
  {
    if (!isCentralitySelected(collisions))
      return;

    fillEventHists(bc, collisions, clusters);

    fillClusterHists(clusters);

    reconstructMesons(clusters, bc.globalIndex());
  }
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EmcalBcWisePi0>(cfgc)}; }
