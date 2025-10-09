// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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
/// \file emcalBcWiseGammaGamma.cxx
///
/// \brief Task that extracts pi0s and eta mesons from BC wise derived data of EMCal clusters
///
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) Goethe University Frankfurt
///

#include "PWGEM/PhotonMeson/DataModel/bcWiseTables.h"

#include "CommonConstants/MathConstants.h"
#include "EMCALBase/Geometry.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Math/AxisAngle.h"
#include "Math/LorentzRotation.h"
#include "Math/Rotation3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TString.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using SelectedClusters = soa::Filtered<aod::BCWiseClusters>;
using SelectedMCClusters = soa::Filtered<soa::Join<aod::BCWiseClusters, aod::BCWiseMCClusters>>;

struct EmcalBcWiseGammaGamma {
  HistogramRegistry mHistManager{"EmcalGammaGammaBcWiseHistograms"};

  Configurable<bool> cfgRequirekTVXinEMC{"cfgRequirekTVXinEMC", true, "Reconstruct mesons only in kTVXinEMC triggered BCs"};
  Configurable<bool> cfgRequireEMCCell{"cfgRequireEMCCell", true, "Reconstruct mesons only in BCs containing at least one EMCal cell (workaround for kTVXinEMC trigger)"};
  Configurable<int> cfgSelectOnlyUniqueAmbiguous{"cfgSelectOnlyUniqueAmbiguous", 0, "0: all clusters, 1: only unique clusters, 2: only ambiguous clusters"};
  Configurable<int> cfgCentralityEstimator{"cfgCentralityEstimator", 0, "0: FT0C, 1: FT0M"};

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
  Configurable<int> cfgMinTimeSinceSOF{"cfgMinTimeSinceSOF", -1, "Only analyze events with a time since start of fill larger than this value (in minutes)"};
  Configurable<int> cfgMaxTimeSinceSOF{"cfgMaxTimeSinceSOF", 100000, "Only analyze events with a time since start of fill smaller than this value (in minutes)"};

  ConfigurableAxis cfgCentralityBinning{"cfgCentralityBinning", {VARIABLE_WIDTH, 0.f, 5.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f, 90.f, 100.f, 101.f, 102.f}, "FT0M centrality (%)"};

  Configurable<bool> cfgIsMC{"cfgIsMC", false, "Flag to indicate if the task is running on MC data and should fill MC histograms"};

  Filter clusterDefinitionFilter = aod::bcwisecluster::storedDefinition == cfgClusterDefinition;
  Filter energyFilter = aod::bcwisecluster::storedE > cfgMinClusterEnergy;
  Filter m02Filter = (aod::bcwisecluster::storedNCells == 1 || (aod::bcwisecluster::storedM02 > cfgMinM02 && aod::bcwisecluster::storedM02 < cfgMaxM02));
  Filter timeFilter = (aod::bcwisecluster::storedTime > cfgMinTime && aod::bcwisecluster::storedTime < cfgMaxTime);

  emcal::Geometry* emcalGeom;

  void init(InitContext const&)
  {
    emcalGeom = emcal::Geometry::GetInstanceFromRunNumber(300000);
    const int nEventBins = 6;
    mHistManager.add("Event/nBCs", "Number of BCs;;#bf{FT0M centrality (%)};#bf{#it{N}_{BC}}", HistType::kTH2F, {{nEventBins, -0.5, 5.5}, cfgCentralityBinning});
    mHistManager.add("Event/nCollisions", "Number of Collisions (BCs x P(mu));;#bf{FT0M centrality (%)};#bf{#it{N}_{coll}}", HistType::kTH2F, {{nEventBins, -0.5, 5.5}, cfgCentralityBinning});
    const TString binLabels[nEventBins] = {"All", "FT0", "TVX", "kTVXinEMC", "Cell", "Cluster"};
    for (int iBin = 0; iBin < nEventBins; iBin++) {
      mHistManager.get<TH2>(HIST("Event/nBCs"))->GetXaxis()->SetBinLabel(iBin + 1, binLabels[iBin]);
      mHistManager.get<TH2>(HIST("Event/nCollisions"))->GetXaxis()->SetBinLabel(iBin + 1, binLabels[iBin]);
    }

    mHistManager.add("Event/nCollPerBC", "Number of collisions per BC;#bf{#it{N}_{coll}};#bf{FT0M centrality (%)};#bf{#it{N}_{BC}}", HistType::kTH2F, {{5, -0.5, 4.5}, cfgCentralityBinning});
    mHistManager.add("Event/Z1VsZ2", "Z vertex positions for BCs with two collisions;#bf{#it{z}_{vtx}^{1} (cm)};#bf{#it{z}_{vtx}^{2} (cm)}", HistType::kTH2F, {{150, -15, 15}, {150, -15, 15}});
    mHistManager.add("Event/dZ", "Distance between vertices for BCs with two collisions;#bf{#Delta #it{z}_{vtx} (cm)};#bf{#it{N}_{BC}}", HistType::kTH1F, {{600, -30, 30}});
    mHistManager.add("Event/Mu", "Probablity of a collision in the BC;#bf{#mu};#bf{#it{N}_{BC}}", HistType::kTH1F, {{2000, 0., 0.4}});
    mHistManager.add("Event/TimeSinceSOF", "Time of BC since the start of fill;#bf{t-t_{SOF} (min)};#bf{#it{N}_{BC}}", HistType::kTH1F, {{2400, 0., 1200}});

    mHistManager.add("Event/Centrality", "FT0M centrality;FT0M centrality (%);#bf{#it{N}_{BC}}", HistType::kTH1F, {cfgCentralityBinning});
    mHistManager.add("Event/CentralityVsAmplitude", "FT0M AmplitudeVsCentrality;FT0M Centrality;FT0M Amplitude", HistType::kTH2F, {cfgCentralityBinning, {600, 0, 300000}});

    mHistManager.add("Cluster/E", "Energy of cluster;#bf{#it{E} (GeV)};#bf{FT0M centrality (%)};#bf{#it{N}_{clusters}}", HistType::kTH2F, {{200, 0, 20}, cfgCentralityBinning});
    mHistManager.add("Cluster/M02", "Shape of cluster;#bf{#it{M}_{02}};#bf{FT0M centrality (%)};#bf{#it{N}_{clusters}}", HistType::kTH2F, {{200, 0, 2}, cfgCentralityBinning});
    mHistManager.add("Cluster/Time", "Time of cluster;#bf{#it{t} (ns)};#bf{FT0M centrality (%)};#bf{#it{N}_{clusters}}", HistType::kTH2F, {{200, -100, 100}, cfgCentralityBinning});
    mHistManager.add("Cluster/NCells", "Number of cells per cluster;#bf{#it{N}_{cells}};#bf{FT0M centrality (%)};#bf{#it{N}_{clusters}}", HistType::kTH2F, {{51, 0., 50.5}, cfgCentralityBinning});
    mHistManager.add("Cluster/Exotic", "Is cluster exotic?;#bf{Exotic?};#bf{FT0M centrality (%)};#bf{#it{N}_{clusters}}", HistType::kTH2F, {{2, -0.5, 1.5}, cfgCentralityBinning});
    mHistManager.add("Cluster/EtaPhi", "Eta/Phi distribution of clusters;#eta;#phi;#bf{FT0M centrality (%)};#bf{#it{N}_{clusters}}", HistType::kTH3F, {{400, -0.8, 0.8}, {400, 0, constants::math::TwoPI}, cfgCentralityBinning});

    mHistManager.add("GG/invMassVsPt", "Invariant mass and pT of meson candidates;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})};#bf{FT0M centrality (%)}", HistType::kTH3F, {{400, 0., 0.8}, {300, 0, 30}, cfgCentralityBinning});
    mHistManager.add("GG/invMassVsPtBackground", "Invariant mass and pT of background meson candidates;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})};#bf{FT0M centrality (%)}", HistType::kTH3F, {{400, 0., 0.8}, {300, 0, 30}, cfgCentralityBinning});

    if (cfgIsMC) {
      mHistManager.add("True/clusterERecVsETrue", "True vs reconstructed energy of cluster inducing particle;#bf{#it{E}_{rec} (GeV)};#bf{#it{E}_{true}^{cls inducing part} (GeV)};#bf{FT0M centrality (%)}", HistType::kTH3F, {{200, 0, 20}, {200, 0, 20}, cfgCentralityBinning});
      mHistManager.add("True/pi0_PtRecVsPtTrue", "True vs reconstructed pT of true pi0s;#bf{#it{p}_{T}^{rec} (GeV/#it{c})};#bf{#it{p}_{T}^{true} (GeV/#it{c})};#bf{FT0M centrality (%)}", HistType::kTH3F, {{300, 0, 30}, {300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("True/pi0_invMassVsPt_Primary", "Reconstructed validated primary pi0;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})};#bf{FT0M centrality (%)}", HistType::kTH3F, {{400, 0., 0.8}, {300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("True/pi0_invMassVsPt_Secondary", "Reconstructed validated pi0 from secondary decay;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})};#bf{FT0M centrality (%)}", HistType::kTH3F, {{400, 0., 0.8}, {300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("True/pi0_invMassVsPt_HadronicShower", "Reconstructed validated pi0 from hadronic shower;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})};#bf{FT0M centrality (%)}", HistType::kTH3F, {{400, 0., 0.8}, {300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("True/eta_PtRecVsPtTrue", "True vs reconstructed pT of true eta meson;#bf{#it{p}_{T}^{rec} (GeV/#it{c})};#bf{#it{p}_{T}^{true} (GeV/#it{c})};#bf{FT0M centrality (%)}", HistType::kTH3F, {{300, 0, 30}, {300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("True/eta_invMassVsPt_Primary", "Reconstructed validated primary eta meson;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})};#bf{FT0M centrality (%)}", HistType::kTH3F, {{400, 0., 0.8}, {300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("True/eta_invMassVsPt_Secondary", "Reconstructed validated eta meson from secondary decay;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})};#bf{FT0M centrality (%)}", HistType::kTH3F, {{400, 0., 0.8}, {300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("True/eta_invMassVsPt_HadronicShower", "Reconstructed validated eta meson from hadronic shower;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})};#bf{FT0M centrality (%)}", HistType::kTH3F, {{400, 0., 0.8}, {300, 0, 30}, cfgCentralityBinning});

      mHistManager.add("Generated/pi0_AllBCs", "pT spectrum of generated pi0s in all BCs;#bf{#it{p}_{T} (GeV/#it{c})};#bf{FT0M centrality (%)};#bf{#it{N}_{#pi^{0}}^{gen}}", HistType::kTH2F, {{300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("Generated/pi0_FT0", "pT spectrum of generated pi0s in BCs with found FT0;#bf{#it{p}_{T} (GeV/#it{c})};#bf{FT0M centrality (%)};#bf{#it{N}_{#pi^{0}}^{gen}}", HistType::kTH2F, {{300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("Generated/pi0_TVX", "pT spectrum of generated pi0s in TVX triggered BCs;#bf{#it{p}_{T} (GeV/#it{c})};#bf{FT0M centrality (%)};#bf{#it{N}_{#pi^{0}}^{gen}}", HistType::kTH2F, {{300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("Generated/pi0_kTVXinEMC", "pT spectrum of generated pi0s in kTVXinEMC triggered BCs;#bf{#it{p}_{T} (GeV/#it{c})};#bf{FT0M centrality (%)};#bf{#it{N}_{#pi^{0}}^{gen}}", HistType::kTH2F, {{300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("Accepted/pi0_kTVXinEMC", "pT spectrum of accepted pi0s in kTVXinEMC triggered BCs;#bf{#it{p}_{T} (GeV/#it{c})};#bf{FT0M centrality (%)};#bf{#it{N}_{#pi^{0}}^{acc}}", HistType::kTH2F, {{300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("Generated/eta_AllBCs", "pT spectrum of generated eta mesons in all BCs;#bf{#it{p}_{T} (GeV/#it{c})};#bf{FT0M centrality (%)};#bf{#it{N}_{#eta}^{gen}}", HistType::kTH2F, {{300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("Generated/eta_FT0", "pT spectrum of generated eta mesons in BCs with found FT0;#bf{#it{p}_{T} (GeV/#it{c})};#bf{FT0M centrality (%)};#bf{#it{N}_{#eta}^{gen}}", HistType::kTH2F, {{300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("Generated/eta_TVX", "pT spectrum of generated eta mesons in TVX triggered BCs;#bf{#it{p}_{T} (GeV/#it{c})};#bf{FT0M centrality (%)};#bf{#it{N}_{#eta}^{gen}}", HistType::kTH2F, {{300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("Generated/eta_kTVXinEMC", "pT spectrum of generated eta mesons in kTVXinEMC triggered BCs;#bf{#it{p}_{T} (GeV/#it{c})};#bf{FT0M centrality (%)};#bf{#it{N}_{#eta}^{gen}}", HistType::kTH2F, {{300, 0, 30}, cfgCentralityBinning});
      mHistManager.add("Accepted/eta_kTVXinEMC", "pT spectrum of accepted eta mesons in kTVXinEMC triggered BCs;#bf{#it{p}_{T} (GeV/#it{c})};#bf{FT0M centrality (%)};#bf{#it{N}_{#eta}^{acc}}", HistType::kTH2F, {{300, 0, 30}, cfgCentralityBinning});
    }
  }

  float getCentrality(const auto& bc)
  {
    if (cfgCentralityEstimator == 0)
      return bc.ft0cCentrality();
    else if (cfgCentralityEstimator == 1)
      return bc.ft0mCentrality();
    else
      throw std::runtime_error("Unknown centrality estimator selected");
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
    mHistManager.fill(HIST("Event/nBCs"), 0, getCentrality(bc));
    float mu = bc.mu();
    mHistManager.fill(HIST("Event/Mu"), mu);
    mHistManager.fill(HIST("Event/TimeSinceSOF"), bc.timeSinceSOF() / 60.);
    double p = mu > 0.001 ? mu / (1 - std::exp(-mu)) : 1.; // No pile-up for small mu (protection against division by zero)
    mHistManager.fill(HIST("Event/nCollisions"), 0, getCentrality(bc), p);
    if (bc.hasFT0()) {
      mHistManager.fill(HIST("Event/nBCs"), 1, getCentrality(bc));
      mHistManager.fill(HIST("Event/nCollisions"), 1, getCentrality(bc), p);
    }
    if (bc.hasTVX()) {
      mHistManager.fill(HIST("Event/nBCs"), 2, getCentrality(bc));
      mHistManager.fill(HIST("Event/nCollisions"), 2, getCentrality(bc), p);
    }
    if (bc.haskTVXinEMC()) {
      mHistManager.fill(HIST("Event/nBCs"), 3, getCentrality(bc));
      mHistManager.fill(HIST("Event/nCollisions"), 3, getCentrality(bc), p);
    }
    if (bc.hasEMCCell()) {
      mHistManager.fill(HIST("Event/nBCs"), 4, getCentrality(bc));
      mHistManager.fill(HIST("Event/nCollisions"), 4, getCentrality(bc), p);
    }
    if (clusters.size() > 0) {
      mHistManager.fill(HIST("Event/nBCs"), 5, getCentrality(bc));
      mHistManager.fill(HIST("Event/nCollisions"), 5, getCentrality(bc), p);
    }

    mHistManager.fill(HIST("Event/Centrality"), getCentrality(bc));
    mHistManager.fill(HIST("Event/CentralityVsAmplitude"), getCentrality(bc), bc.ft0Amplitude());

    mHistManager.fill(HIST("Event/nCollPerBC"), collisions.size(), getCentrality(bc));
    if (collisions.size() == 2) {
      mHistManager.fill(HIST("Event/Z1VsZ2"), collisions.iteratorAt(0).zVtx(), collisions.iteratorAt(1).zVtx());
      mHistManager.fill(HIST("Event/dZ"), collisions.iteratorAt(0).zVtx() - collisions.iteratorAt(1).zVtx());
    }
  }

  void fillClusterHists(const auto& clusters, float centrality)
  {
    for (const auto& cluster : clusters) {
      mHistManager.fill(HIST("Cluster/E"), cluster.e(), centrality);
      mHistManager.fill(HIST("Cluster/M02"), cluster.m02(), centrality);
      mHistManager.fill(HIST("Cluster/Time"), cluster.time(), centrality);
      mHistManager.fill(HIST("Cluster/NCells"), cluster.nCells(), centrality);
      mHistManager.fill(HIST("Cluster/EtaPhi"), cluster.eta(), cluster.phi(), centrality);
      mHistManager.fill(HIST("Cluster/Exotic"), cluster.isExotic(), centrality);
    }
  }

  void reconstructMesons(const auto& clusters, const auto& bc)
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

      mHistManager.fill(HIST("GG/invMassVsPt"), v12.M(), v12.Pt(), getCentrality(bc));

      if (clusters.size() < 3)
        continue;

      if (bc.globalIndex() % cfgBGEventDownsampling != 0)
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

          mHistManager.fill(HIST("GG/invMassVsPtBackground"), vBG.M(), vBG.Pt(), getCentrality(bc));
        }
      }
    }
  }
  void reconstructTrueMesons(const auto& clusters, const auto& mcPi0s, const auto& mcEtas, const auto& bc)
  {
    for (const auto& [g1, g2] : soa::combinations(soa::CombinationsStrictlyUpperIndexPolicy(clusters, clusters))) {
      if (g1.mesonID() != g2.mesonID() || g1.mesonID() == -1)
        continue;

      ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      if (std::fabs(v12.Rapidity()) > cfgRapidityCut)
        continue;

      float openingAngle12 = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));
      if (openingAngle12 < cfgMinOpenAngle)
        continue;

      if (!g1.isEta()) {
        const auto& mcPi0 = mcPi0s.iteratorAt(g1.mesonID() - mcPi0s.offset());

        mHistManager.fill(HIST("True/pi0_PtRecVsPtTrue"), v12.Pt(), mcPi0.pt(), getCentrality(bc));

        if (mcPi0.isPrimary())
          mHistManager.fill(HIST("True/pi0_invMassVsPt_Primary"), v12.M(), v12.Pt(), getCentrality(bc));
        else if (mcPi0.isFromWD())
          mHistManager.fill(HIST("True/pi0_invMassVsPt_Secondary"), v12.M(), v12.Pt(), getCentrality(bc));
        else
          mHistManager.fill(HIST("True/pi0_invMassVsPt_HadronicShower"), v12.M(), v12.Pt(), getCentrality(bc));
      } else {
        const auto& mcEta = mcEtas.iteratorAt(g1.mesonID() - mcEtas.offset());

        mHistManager.fill(HIST("True/eta_PtRecVsPtTrue"), v12.Pt(), mcEta.pt(), getCentrality(bc));

        if (mcEta.isPrimary())
          mHistManager.fill(HIST("True/eta_invMassVsPt_Primary"), v12.M(), v12.Pt(), getCentrality(bc));
        else if (mcEta.isFromWD())
          mHistManager.fill(HIST("True/eta_invMassVsPt_Secondary"), v12.M(), v12.Pt(), getCentrality(bc));
        else
          mHistManager.fill(HIST("True/eta_invMassVsPt_HadronicShower"), v12.M(), v12.Pt(), getCentrality(bc));
      }
    }
  }

  bool isBCSelected(const auto& bc, const auto& collisions)
  {
    if (cfgRequirekTVXinEMC && !bc.haskTVXinEMC())
      return false;
    if (cfgRequireEMCCell && !bc.hasEMCCell())
      return false;
    if (cfgSelectOnlyUniqueAmbiguous == 1 && collisions.size() != 1)
      return false;
    if (cfgSelectOnlyUniqueAmbiguous == 2 && collisions.size() == 1)
      return false;
    if (cfgMinTimeSinceSOF > bc.timeSinceSOF() / 60 || cfgMaxTimeSinceSOF < bc.timeSinceSOF() / 60)
      return false;
    return true;
  }

  void fillGeneratedMesonHists(const auto& mcPi0s, const auto& mcEtas, const auto& bc)
  {
    for (const auto& mcPi0 : mcPi0s) {
      if (mcPi0.isPrimary()) {
        mHistManager.fill(HIST("Generated/pi0_AllBCs"), mcPi0.pt(), getCentrality(bc));
        if (bc.hasFT0())
          mHistManager.fill(HIST("Generated/pi0_FT0"), mcPi0.pt(), getCentrality(bc));
        if (bc.hasTVX())
          mHistManager.fill(HIST("Generated/pi0_TVX"), mcPi0.pt(), getCentrality(bc));
        if (bc.haskTVXinEMC())
          mHistManager.fill(HIST("Generated/pi0_kTVXinEMC"), mcPi0.pt(), getCentrality(bc));
        if (mcPi0.isAccepted() && bc.haskTVXinEMC())
          mHistManager.fill(HIST("Accepted/pi0_kTVXinEMC"), mcPi0.pt(), getCentrality(bc));
      }
    }
    for (const auto& mcEta : mcEtas) {
      if (mcEta.isPrimary()) {
        mHistManager.fill(HIST("Generated/eta_AllBCs"), mcEta.pt(), getCentrality(bc));
        if (bc.hasFT0())
          mHistManager.fill(HIST("Generated/eta_FT0"), mcEta.pt(), getCentrality(bc));
        if (bc.hasTVX())
          mHistManager.fill(HIST("Generated/eta_TVX"), mcEta.pt(), getCentrality(bc));
        if (bc.haskTVXinEMC())
          mHistManager.fill(HIST("Generated/eta_kTVXinEMC"), mcEta.pt(), getCentrality(bc));
        if (mcEta.isAccepted() && bc.haskTVXinEMC())
          mHistManager.fill(HIST("Accepted/eta_kTVXinEMC"), mcEta.pt(), getCentrality(bc));
      }
    }
  }

  void process(aod::BCWiseBCs::iterator const& bc, aod::BCWiseCollisions const& collisions, SelectedClusters const& clusters)
  {
    if (!isBCSelected(bc, collisions))
      return;

    fillEventHists(bc, collisions, clusters);

    fillClusterHists(clusters, getCentrality(bc));

    reconstructMesons(clusters, bc);
  }

  void processMCInfo(aod::BCWiseBCs::iterator const& bc, aod::BCWiseCollisions const& collisions, SelectedMCClusters const& clusters, aod::BCWiseMCPi0s const& mcPi0s, aod::BCWiseMCEtas const& mcEtas)
  {
    if (!cfgIsMC)
      LOG(fatal) << "MC processing is not enabled, but the task is running on MC data. Please set cfgIsMC to true.";

    fillGeneratedMesonHists(mcPi0s, mcEtas, bc); // Fill before BC selection to also store pi0s and eta mesons in BCs that were not triggered

    if (!isBCSelected(bc, collisions))
      return;

    for (const auto& cluster : clusters)
      mHistManager.fill(HIST("True/clusterERecVsETrue"), cluster.e(), cluster.trueE(), getCentrality(bc));

    reconstructTrueMesons(clusters, mcPi0s, mcEtas, bc);
  }
  PROCESS_SWITCH(EmcalBcWiseGammaGamma, processMCInfo, "Run true and gen", false);
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EmcalBcWiseGammaGamma>(cfgc)}; }
