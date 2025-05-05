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

using SelectedClusters = soa::Filtered<aod::BCWiseClusters>;
using SelectedMCClusters = soa::Filtered<soa::Join<aod::BCWiseClusters, aod::BCWiseMCClusters>>;

struct EmcalBcWisePi0 {
  HistogramRegistry mHistManager{"EmcalGammaGammaBcWiseHistograms"};

  Configurable<bool> cfgRequirekTVXinEMC{"cfgRequirekTVXinEMC", true, "Reconstruct pi0s only in kTVXinEMC triggered BCs"};
  Configurable<int> cfgSelectOnlyUniqueAmbiguous{"cfgSelectOnlyUniqueAmbiguous", 0, "0: all clusters, 1: only unique clusters, 2: only ambiguous clusters"};

  Configurable<int> cfgClusterDefinition{"cfgClusterDefinition", 13, "Clusterizer to be selected, e.g. 13 for kV3MostSplitLowSeed"};
  Configurable<int16_t> cfgMinClusterEnergy{"cfgMinClusterEnergy", 700, "Minimum energy of selected clusters (MeV)"};
  Configurable<int16_t> cfgMinM02{"cfgMinM02", 1000, "Minimum M02 of selected clusters (x1000)"};
  Configurable<int16_t> cfgMaxM02{"cfgMaxM02", 7000, "Maximum M02 of selected clusters (x1000)"};
  Configurable<int16_t> cfgMinTime{"cfgMinTime", -1500, "Minimum time of selected clusters (10 ps)"};
  Configurable<int16_t> cfgMaxTime{"cfgMaxTime", 1500, "Maximum time of selected clusters (10 ps)"};
  Configurable<float> cfgRapidityCut{"cfgRapidityCut", 0.8f, "Maximum absolute rapidity of counted particles"};
  Configurable<float> cfgMinOpenAngle{"cfgMinOpenAngle", 0.0202, "Minimum opening angle between photons"};
  Configurable<int> cfgDistanceToEdge{"cfgDistanceToEdge", 1, "Distance to edge in cells required for rotated cluster to be accepted"};
  Configurable<int> cfgBGEventDownsampling{"cfgBGEventDownsampling", 1, "Only calculate background for every n-th event (performance reasons in PbPb)"};

  static constexpr float DefaultCentralityWindow[2] = {-1., 101.};
  Configurable<LabeledArray<float>> cfgCentralityWindow{"cfgCentralityWindow", {DefaultCentralityWindow, 2, {"Min_Centrality", "Max_Centrality"}}, "Select centrality window (also requires unique collision)"};

  Configurable<bool> cfgIsMC{"cfgIsMC", false, "Flag to indicate if the task is running on MC data and should fill MC histograms"};

  Filter clusterDefinitionFilter = aod::bcwisecluster::storedDefinition == cfgClusterDefinition;
  Filter energyFilter = aod::bcwisecluster::storedE > cfgMinClusterEnergy;
  Filter m02Filter = (aod::bcwisecluster::storedNCells == 1 || (aod::bcwisecluster::storedM02 > cfgMinM02 && aod::bcwisecluster::storedM02 < cfgMaxM02));
  Filter timeFilter = (aod::bcwisecluster::storedTime > cfgMinTime && aod::bcwisecluster::storedTime < cfgMaxTime);

  emcal::Geometry* emcalGeom;

  void init(InitContext const&)
  {
    emcalGeom = emcal::Geometry::GetInstanceFromRunNumber(300000);
    const int nEventBins = 8;
    mHistManager.add("Event/nBCs", "Number of BCs;;#bf{#it{N}_{BC}}", HistType::kTH1F, {{nEventBins, -0.5, 7.5}});
    const TString binLabels[nEventBins] = {"All", "FT0", "TVX", "kTVXinEMC", "Cell", "Cluster", "NoBorder", "Collision"};
    for (int iBin = 0; iBin < nEventBins; iBin++)
      mHistManager.get<TH1>(HIST("Event/nBCs"))->GetXaxis()->SetBinLabel(iBin + 1, binLabels[iBin]);

    mHistManager.add("Event/nCollPerBC", "Number of collisions per BC;#bf{#it{N}_{coll}};#bf{#it{N}_{BC}}", HistType::kTH1F, {{5, -0.5, 4.5}});
    mHistManager.add("Event/Z1VsZ2", "Z vertex positions for BCs with two collisions;#bf{#it{z}_{vtx}^{1} (cm)};#bf{#it{z}_{vtx}^{2} (cm)}", HistType::kTH2F, {{150, -15, 15}, {150, -15, 15}});
    mHistManager.add("Event/dZ", "Distance between vertices for BCs with two collisions;#bf{#Delta #it{z}_{vtx} (cm)};#bf{#it{N}_{BC}}", HistType::kTH1F, {{600, -30, 30}});

    mHistManager.add("Event/Centrality", "FT0M centrality;FT0M centrality (%);#bf{#it{N}_{BC} (only BCs containing exactly 1 collision)}", HistType::kTH1F, {{100, 0., 100.}});
    mHistManager.add("Event/CentralityVsAmplitude", "FT0M AmplitudeVsCentrality;FT0M Centrality;FT0M Amplitude", HistType::kTH2F, {{105, 0, 105}, {600, 0, 300000}});

    mHistManager.add("Cluster/E", "Energy of cluster;#bf{#it{E} (GeV)};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, 0, 20}});
    mHistManager.add("Cluster/M02", "Shape of cluster;#bf{#it{M}_{02}};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, 0, 2}});
    mHistManager.add("Cluster/Time", "Time of cluster;#bf{#it{t} (ns)};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, -100, 100}});
    mHistManager.add("Cluster/NCells", "Number of cells per cluster;#bf{#it{N}_{cells}};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{51, 0., 50.5}});
    mHistManager.add("Cluster/Exotic", "Is cluster exotic?;#bf{Exotic?};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{2, -0.5, 1.5}});
    mHistManager.add("Cluster/EtaPhi", "Eta/Phi distribution of clusters;#eta;#phi", HistType::kTH2F, {{400, -0.8, 0.8}, {400, 0, constants::math::TwoPI}});

    mHistManager.add("GG/invMassVsPt", "Invariant mass and pT of meson candidates;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})}", HistType::kTH2F, {{200, 0., 0.4}, {200, 0., 20.}});
    mHistManager.add("GG/invMassVsPtBackground", "Invariant mass and pT of background meson candidates;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})}", HistType::kTH2F, {{200, 0., 0.4}, {200, 0., 20.}});

    if (cfgIsMC) {
      mHistManager.add("True/clusterERecVsETrue", "True vs reconstructed energy of cluster inducing particle;#bf{#it{E}_{rec} (GeV)};#bf{#it{E}_{true}^{cls inducing part} (GeV)}", HistType::kTH2F, {{200, 0, 20}, {200, 0, 20}});
      mHistManager.add("True/PtRecVsPtTrue", "True vs reconstructed pT of true pi0s;#bf{#it{p}_{T}^{rec} (GeV/#it{c})};#bf{#it{p}_{T}^{true} (GeV/#it{c})}", HistType::kTH2F, {{200, 0., 20.}, {200, 0., 20.}});
      mHistManager.add("True/invMassVsPt_Primary", "Invariant mass and pT of meson candidates", HistType::kTH2F, {{200, 0., 0.4}, {200, 0., 20.}});
      mHistManager.add("True/invMassVsPt_Secondary", "Invariant mass and pT of meson candidates", HistType::kTH2F, {{200, 0., 0.4}, {200, 0., 20.}});
      mHistManager.add("True/invMassVsPt_HadronicShower", "Invariant mass and pT of meson candidates", HistType::kTH2F, {{200, 0., 0.4}, {200, 0., 20.}});

      mHistManager.add("Generated/pi0_AllBCs", "pT spectrum of generated pi0s in all BCs;#bf{#it{p}_{T} (GeV/#it{c})};#bf{#it{N}_{#pi^{0}}^{gen}}", HistType::kTH1F, {{200, 0, 20}});
      mHistManager.add("Generated/pi0_TVX", "pT spectrum of generated pi0s in TVX triggered BCs;#bf{#it{p}_{T} (GeV/#it{c})};#bf{#it{N}_{#pi^{0}}^{gen}}", HistType::kTH1F, {{200, 0, 20}});
      mHistManager.add("Generated/pi0_kTVXinEMC", "pT spectrum of generated pi0s in kTVXinEMC triggered BCs;#bf{#it{p}_{T} (GeV/#it{c})};#bf{#it{N}_{#pi^{0}}^{gen}}", HistType::kTH1F, {{200, 0, 20}});
    }
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
    mHistManager.fill(HIST("Event/nBCs"), 0);
    if (bc.hasFT0())
      mHistManager.fill(HIST("Event/nBCs"), 1);
    if (bc.hasTVX())
      mHistManager.fill(HIST("Event/nBCs"), 2);
    if (bc.haskTVXinEMC())
      mHistManager.fill(HIST("Event/nBCs"), 3);
    if (bc.hasEMCCell())
      mHistManager.fill(HIST("Event/nBCs"), 4);
    if (clusters.size() > 0)
      mHistManager.fill(HIST("Event/nBCs"), 5);
    if (bc.hasNoTFROFBorder())
      mHistManager.fill(HIST("Event/nBCs"), 6);
    if (collisions.size() > 0)
      mHistManager.fill(HIST("Event/nBCs"), 7);

    if (collisions.size() == 1) {
      mHistManager.fill(HIST("Event/Centrality"), collisions.iteratorAt(0).centrality());
      mHistManager.fill(HIST("Event/CentralityVsAmplitude"), collisions.iteratorAt(0).centrality(), bc.ft0Amplitude());
    } else {
      mHistManager.fill(HIST("Event/CentralityVsAmplitude"), 103, bc.ft0Amplitude());
    }

    mHistManager.fill(HIST("Event/nCollPerBC"), collisions.size());
    if (collisions.size() == 2) {
      mHistManager.fill(HIST("Event/Z1VsZ2"), collisions.iteratorAt(0).zVtx(), collisions.iteratorAt(1).zVtx());
      mHistManager.fill(HIST("Event/dZ"), collisions.iteratorAt(0).zVtx() - collisions.iteratorAt(1).zVtx());
    }
  }

  void fillClusterHists(const auto& clusters)
  {
    for (const auto& cluster : clusters) {
      mHistManager.fill(HIST("Cluster/E"), cluster.e());
      mHistManager.fill(HIST("Cluster/M02"), cluster.m02());
      mHistManager.fill(HIST("Cluster/Time"), cluster.time());
      mHistManager.fill(HIST("Cluster/NCells"), cluster.nCells());
      mHistManager.fill(HIST("Cluster/EtaPhi"), cluster.eta(), cluster.phi());
      mHistManager.fill(HIST("Cluster/Exotic"), cluster.isExotic());
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

      mHistManager.fill(HIST("GG/invMassVsPt"), v12.M(), v12.Pt());

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

          mHistManager.fill(HIST("GG/invMassVsPtBackground"), vBG.M(), vBG.Pt());
        }
      }
    }
  }
  void reconstructTrueMesons(const auto& clusters, const auto& mcPi0s)
  {
    for (const auto& [g1, g2] : soa::combinations(soa::CombinationsStrictlyUpperIndexPolicy(clusters, clusters))) {
      if (g1.pi0ID() == -1 || g1.pi0ID() != g2.pi0ID())
        continue;

      ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      if (std::fabs(v12.Rapidity()) > cfgRapidityCut)
        continue;

      float openingAngle12 = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));
      if (openingAngle12 < cfgMinOpenAngle)
        continue;

      const auto& mcPi0 = mcPi0s.iteratorAt(g1.pi0ID() - mcPi0s.offset());

      mHistManager.fill(HIST("True/PtRecVsPtTrue"), v12.Pt(), mcPi0.pt());

      if (mcPi0.isPrimary())
        mHistManager.fill(HIST("True/invMassVsPt_Primary"), v12.M(), v12.Pt());
      else if (mcPi0.isFromWD())
        mHistManager.fill(HIST("True/invMassVsPt_Secondary"), v12.M(), v12.Pt());
      else
        mHistManager.fill(HIST("True/invMassVsPt_HadronicShower"), v12.M(), v12.Pt());
    }
  }

  bool isBCSelected(const auto& bc, const auto& collisions)
  {
    if (cfgRequirekTVXinEMC && !bc.haskTVXinEMC())
      return false;

    if (cfgCentralityWindow->get("Min_Centrality") > -0.5 || cfgCentralityWindow->get("Max_Centrality") < 100.5) { // Centrality window is set
      if (collisions.size() != 1)
        return false;
      if (collisions.iteratorAt(0).centrality() < cfgCentralityWindow->get("Min_Centrality") || collisions.iteratorAt(0).centrality() > cfgCentralityWindow->get("Max_Centrality"))
        return false;
    }

    if (cfgSelectOnlyUniqueAmbiguous == 1 && collisions.size() != 1)
      return false;
    if (cfgSelectOnlyUniqueAmbiguous == 2 && collisions.size() == 1)
      return false;

    return true;
  }

  void fillGeneratedPi0Hists(const auto& mcPi0s, const auto& bc)
  {
    for (const auto& mcPi0 : mcPi0s) {
      mHistManager.fill(HIST("Generated/pi0_AllBCs"), mcPi0.pt());
      if (bc.hasTVX())
        mHistManager.fill(HIST("Generated/pi0_TVX"), mcPi0.pt());
      if (bc.haskTVXinEMC())
        mHistManager.fill(HIST("Generated/pi0_kTVXinEMC"), mcPi0.pt());
    }
  }

  void process(aod::BCWiseBCs::iterator const& bc, aod::BCWiseCollisions const& collisions, SelectedClusters const& clusters)
  {
    if (!isBCSelected(bc, collisions))
      return;

    fillEventHists(bc, collisions, clusters);

    fillClusterHists(clusters);

    reconstructMesons(clusters, bc.globalIndex());
  }

  void processMCInfo(aod::BCWiseBCs::iterator const& bc, aod::BCWiseCollisions const& collisions, SelectedMCClusters const& clusters, aod::BCWiseMCPi0s const& mcPi0s)
  {
    if (!cfgIsMC) {
      LOG(fatal) << "MC processing is not enabled, but the task is running on MC data. Please set cfgIsMC to true.";
      return;
    }

    fillGeneratedPi0Hists(mcPi0s, bc); // Fill before BC selection to also store pi0s in BCs that were not triggered

    if (!isBCSelected(bc, collisions))
      return;

    for (const auto& cluster : clusters)
      mHistManager.fill(HIST("True/clusterERecVsETrue"), cluster.e(), cluster.trueE());

    reconstructTrueMesons(clusters, mcPi0s);
  }
  PROCESS_SWITCH(EmcalBcWisePi0, processMCInfo, "Run true and gen", false);
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EmcalBcWisePi0>(cfgc)}; }
