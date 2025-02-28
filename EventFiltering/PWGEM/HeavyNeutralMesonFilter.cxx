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
/// \file HeavyNeutralMesonFilter.cxx
///
/// \brief This code loops over collisions to filter events contaning heavy mesons (omega or eta') using EMCal clusters
///
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt
///

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <vector>
#include "TVector3.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/EventSelection.h"
#include "EMCALBase/Geometry.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyBCs = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps>;
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;

using SelectedTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TrackSelection>>;
using SelectedClusters = soa::Filtered<aod::EMCALClusters>; // Clusters from collisions with only one collision in the BC

struct Photon { // Struct to store photons (unique and ambiguous clusters that passed the cuts)
  Photon(float eta, float phi, float energy) : eta(eta), phi(phi), e(energy), theta(2 * std::atan2(std::exp(-eta), 1)), px(e * std::sin(theta) * std::cos(phi)), py(e * std::sin(theta) * std::sin(phi)), pz(e * std::cos(theta)), pt(std::sqrt(px * px + py * py))
  {
    photon.SetPxPyPzE(px, py, pz, e);
  }

  ROOT::Math::PxPyPzEVector photon;
  float eta, phi, e;
  float theta;
  float px, py, pz, pt;
};

namespace HeavyNeutralMeson
{
enum MesonType {
  kPi0 = 1 << 0,     // 0001
  kEta = 1 << 1,     // 0010
  kOmega = 1 << 2,   // 0100
  kEtaPrime = 1 << 3 // 1000
};
}

struct Meson {
  Meson(Photon p1, Photon p2) : p1(p1), p2(p2), mesonType(0)
  {
    pMeson = p1.photon + p2.photon;
  }
  Meson(Meson gg, float eTracks, float pxTracks, float pyTracks, float pzTracks) : p1(gg.p1), p2(gg.p2), mesonType(0)
  {
    pMeson.SetPxPyPzE(gg.pMeson.Px() + pxTracks, gg.pMeson.Py() + pyTracks, gg.pMeson.Pz() + pzTracks, gg.pMeson.e() + eTracks);
  }
  Photon p1, p2;
  ROOT::Math::PxPyPzEVector pMeson;

  uint32_t mesonType;

  void addMesonType(uint32_t type) { mesonType |= type; }
  bool hasMesonType(uint32_t type) const { return mesonType & type; }

  float m() const { return pMeson.M(); }
  float pT() const { return pMeson.Pt(); }
};

struct HeavyNeutralMesonFilter {

  Produces<aod::HeavyNeutralMesonFilters> tags;

  HistogramRegistry mHistManager{"HeavyNeutralMesonFilterHistograms", {}, OutputObjHandlingPolicy::QAObject};

  Configurable<int> cfgClusterDefinition{"cfgClusterDefinition", 13, "Clusterizer to be selected, e.g. 13 for kV3MostSplitLowSeed"};
  Configurable<float> cfgMinClusterEnergy{"cfgMinClusterEnergy", 0.5, "Minimum energy of selected clusters (GeV)"};
  Configurable<float> cfgMinM02{"cfgMinM02", 0., "Minimum M02 of selected clusters"};
  Configurable<float> cfgMaxM02{"cfgMaxM02", 1, "Maximum M02 of selected clusters"};
  Configurable<float> cfgMinTime{"cfgMinTime", -25, "Minimum time of selected clusters (ns)"};
  Configurable<float> cfgMaxTime{"cfgMaxTime", 25, "Maximum time of selected clusters (ns)"};

  Configurable<float> cfgTrackMinPt{"cfgTrackMinPt", 0.1, "Minimum momentum of tracks (GeV/c)"};

  static constexpr float defaultMassWindows[2][4] = {{0., 0.4, 0.6, 1.}, {0.4, 0.8, 0.8, 1.2}};
  Configurable<LabeledArray<float>> massWindowOmega{"massWindowOmega", {defaultMassWindows[0], 4, {"pi0_min", "pi0_max", "omega_min", "omega_max"}}, "Mass window for selected omegas and their decay pi0"};
  Configurable<LabeledArray<float>> massWindowEtaPrime{"massWindowEtaPrime", {defaultMassWindows[1], 4, {"eta_min", "eta_max", "etaprime_min", "etaprime_max"}}, "Mass window for selected eta' and their decay eta"};

  static constexpr float defaultMinPts[4] = {0., 5., 0., 5.}; // LowPtOmega, HighPtOmega, LowPtEtaPrime, HighPtEtaPrime
  Configurable<LabeledArray<float>> minHNMPt{"minHNMPt", {defaultMinPts, 4, {"omega_lowPt", "omega_highPt", "etaPrime_lowPt", "etaPrime_highPt"}}, "Minimum pT values for the heavy neutral mesons of the low and high pT trigger"};

  Filter clusterDefinitionFilter = aod::emcalcluster::definition == cfgClusterDefinition;
  Filter energyFilter = aod::emcalcluster::energy > cfgMinClusterEnergy;
  Filter m02Filter = (aod::emcalcluster::nCells == 1 || (aod::emcalcluster::m02 > cfgMinM02 && aod::emcalcluster::m02 < cfgMaxM02));
  Filter timeFilter = (aod::emcalcluster::time > cfgMinTime && aod::emcalcluster::time < cfgMaxTime);

  Filter trackPtFilter = aod::track::pt > cfgTrackMinPt;

  std::vector<Photon> mPhotons;
  std::vector<Meson> mGGs;

  bool colContainsLowPtOmega = false;
  bool colContainsHighPtOmega = false;
  bool colContainsLowPtEtaPrime = false;
  bool colContainsHighPtEtaPrime = false;

  emcal::Geometry* emcalGeom;

  void init(InitContext const&)
  {
    emcalGeom = emcal::Geometry::GetInstanceFromRunNumber(300000);
    auto hCollisionCounter = mHistManager.add<TH1>("Event/hCollisionCounter", "Number of collisions;;#bf{#it{N}_{Coll}}", HistType::kTH1F, {{8, -0.5, 7.5}});
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "kTVXinEMC");
    hCollisionCounter->GetXaxis()->SetBinLabel(3, "L0Triggered");
    hCollisionCounter->GetXaxis()->SetBinLabel(4, "2+ tracks&clusters");
    hCollisionCounter->GetXaxis()->SetBinLabel(5, "Low #it{p}_{T} #omega");
    hCollisionCounter->GetXaxis()->SetBinLabel(6, "High #it{p}_{T} #omega");
    hCollisionCounter->GetXaxis()->SetBinLabel(7, "Low #it{p}_{T} #eta'");
    hCollisionCounter->GetXaxis()->SetBinLabel(8, "High #it{p}_{T} #eta'");

    mHistManager.add("Event/nClusters", "Number of BCs (with (k)TVX);#bf{TVX triggered};#bf{#it{N}_{BCs}}", HistType::kTH1F, {{51, -0.5, 50.5}});
    mHistManager.add("Event/nGGs", "Number of BCs (with (k)TVX);#bf{TVX triggered};#bf{#it{N}_{BCs}}", HistType::kTH1F, {{51, -0.5, 50.5}});
    mHistManager.add("Event/nTracks", "Number of BCs (with (k)TVX);#bf{TVX triggered};#bf{#it{N}_{BCs}}", HistType::kTH1F, {{51, -0.5, 50.5}});
    mHistManager.add("Event/nHeavyNeutralMesons", "Number of BCs (with (k)TVX);#bf{TVX triggered};#bf{#it{N}_{BCs}}", HistType::kTH1F, {{51, -0.5, 50.5}});

    mHistManager.add("nCollisionsVsClusters", "Number of collisions vs number of clusters;N_{Collisions};N_{Clusters}", HistType::kTH2F, {{4, -0.5, 3.5}, {21, -0.5, 20.5}});

    mHistManager.add("Track/trackP", "Momentum of tracks;#bf{#it{p} (GeV/#it{c})};#bf{#it{N}_{tracks}}", HistType::kTH1F, {{200, 0, 20}});

    mHistManager.add("Cluster/clusterE", "Energy of cluster;#bf{#it{E} (GeV)};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, 0, 20}});
    mHistManager.add("Cluster/clusterM02", "Shape of cluster;#bf{#it{M}_{02}};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, 0, 2}});
    mHistManager.add("Cluster/clusterTime", "Time of cluster;#bf{#it{t} (ns)};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, -100, 100}});

    mHistManager.add("GG/invMassVsPt", "Invariant mass and pT of meson candidates", HistType::kTH2F, {{400, massWindowOmega->get("pi0_min"), massWindowEtaPrime->get("eta_max")}, {250, 0., 25.}});

    mHistManager.add("HeavyNeutralMeson/invMassVsPt", "Invariant mass and pT of meson candidates", HistType::kTH2F, {{600, massWindowOmega->get("omega_min"), massWindowEtaPrime->get("etaprime_max")}, {250, 0., 25.}});
  }

  void process(MyCollisions::iterator const& collision, MyBCs const&, SelectedClusters const& clusters, SelectedTracks const& tracks)
  {
    mHistManager.fill(HIST("Event/hCollisionCounter"), 0.);

    colContainsLowPtOmega = false;
    colContainsHighPtOmega = false;
    colContainsLowPtEtaPrime = false;
    colContainsHighPtEtaPrime = false;

    if (isBCSelected(collision.foundBC_as<MyBCs>(), clusters.size(), tracks.size() > 2)) {
      processClusters(clusters);

      processGGs();
      processHeavyNeutralMesons(tracks);

      mPhotons.clear();
      mGGs.clear();
    }

    if (colContainsLowPtOmega)
      mHistManager.fill(HIST("Event/hCollisionCounter"), 4.);
    if (colContainsHighPtOmega)
      mHistManager.fill(HIST("Event/hCollisionCounter"), 5.);
    if (colContainsLowPtEtaPrime)
      mHistManager.fill(HIST("Event/hCollisionCounter"), 6.);
    if (colContainsHighPtEtaPrime)
      mHistManager.fill(HIST("Event/hCollisionCounter"), 7.);
    tags(colContainsLowPtOmega, colContainsHighPtOmega, colContainsLowPtEtaPrime, colContainsHighPtEtaPrime);
  }

  bool isBCSelected(MyBCs::iterator const& bc, int nClusters, bool enoughTracks)
  {
    if (bc.alias_bit(kTVXinEMC))
      mHistManager.fill(HIST("Event/hCollisionCounter"), 1.);
    else if (bc.alias_bit(kEMC7) || bc.alias_bit(kEG1) || bc.alias_bit(kEG2) || bc.alias_bit(kDG1) || bc.alias_bit(kDG2) || bc.alias_bit(kEJ1) || bc.alias_bit(kEJ2) || bc.alias_bit(kDJ1) || bc.alias_bit(kDJ2))
      mHistManager.fill(HIST("Event/hCollisionCounter"), 2.);
    else
      return false;
    if (nClusters < 2 || !enoughTracks)
      return false;
    mHistManager.fill(HIST("Event/hCollisionCounter"), 3.);

    return true;
  }

  void processClusters(SelectedClusters const& clusters)
  {
    for (const auto& cluster : clusters) {
      mHistManager.fill(HIST("Cluster/clusterE"), cluster.energy());
      mHistManager.fill(HIST("Cluster/clusterM02"), cluster.m02());
      mHistManager.fill(HIST("Cluster/clusterTime"), cluster.time());
      mPhotons.push_back(Photon(cluster.eta(), cluster.phi(), cluster.energy()));
    }
    mHistManager.fill(HIST("Event/nClusters"), clusters.size());
  }

  /// \brief Process light neutral meson candidates (pi0 or eta), calculate invariant mass and pT and fill histograms
  void processGGs()
  {
    // loop over all photon combinations and build meson candidates
    for (unsigned int ig1 = 0; ig1 < mPhotons.size(); ++ig1) {
      for (unsigned int ig2 = ig1 + 1; ig2 < mPhotons.size(); ++ig2) {
        Meson lightMeson(mPhotons[ig1], mPhotons[ig2]); // build lightMeson from photons
        mHistManager.fill(HIST("GG/invMassVsPt"), lightMeson.m(), lightMeson.pT());
        if (lightMeson.m() > massWindowOmega->get("pi0_min") && lightMeson.m() < massWindowOmega->get("pi0_max"))
          lightMeson.addMesonType(HeavyNeutralMeson::kPi0);
        else if (lightMeson.m() > massWindowEtaPrime->get("eta_min") && lightMeson.m() < massWindowEtaPrime->get("eta_max"))
          lightMeson.addMesonType(HeavyNeutralMeson::kEta);
        else
          continue;
        mGGs.push_back(lightMeson);
      }
    }
    mHistManager.fill(HIST("Event/nGGs"), mGGs.size());
  }

  void processHeavyNeutralMesons(SelectedTracks const& tracks)
  {
    for (const auto& posTrack : tracks) {
      if (!posTrack.isGlobalTrack() || posTrack.sign() < 0)
        continue;
      for (const auto& negTrack : tracks) {
        if (!negTrack.isGlobalTrack() || negTrack.sign() > 0)
          continue;
        for (unsigned long iGG = 0; iGG < mGGs.size(); iGG++) {
          Meson heavyNeutralMeson(mGGs.at(iGG), posTrack.energy(constants::physics::MassPiPlus) + negTrack.energy(constants::physics::MassPiPlus), posTrack.px() + negTrack.px(), posTrack.py() + negTrack.py(), posTrack.pz() + negTrack.pz());
          mHistManager.fill(HIST("HeavyNeutralMeson/invMassVsPt"), heavyNeutralMeson.m(), heavyNeutralMeson.pT());
          if (mGGs.at(iGG).hasMesonType(HeavyNeutralMeson::kPi0) && heavyNeutralMeson.m() > massWindowOmega->get("omega_min") && heavyNeutralMeson.m() < massWindowOmega->get("omega_max")) {
            if(minHNMPt->get("omega_lowPt") < heavyNeutralMeson.pT())
              colContainsLowPtOmega = true;
            if(minHNMPt->get("omega_highPt") < heavyNeutralMeson.pT())
              colContainsHighPtOmega = true;
          }
          if (mGGs.at(iGG).hasMesonType(HeavyNeutralMeson::kEta) && heavyNeutralMeson.m() > massWindowEtaPrime->get("etaprime_min") && heavyNeutralMeson.m() < massWindowEtaPrime->get("etaprime_max")) {
            if(minHNMPt->get("etaPrime_lowPt") < heavyNeutralMeson.pT())
            colContainsLowPtEtaPrime = true;
          if(minHNMPt->get("etaPrime_highPt") < heavyNeutralMeson.pT())
            colContainsHighPtEtaPrime = true;
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HeavyNeutralMesonFilter>(cfgc)};
}
