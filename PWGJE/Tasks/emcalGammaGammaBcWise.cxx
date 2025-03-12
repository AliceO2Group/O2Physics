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
/// \file emcalGammaGammaBcWise.cxx
///
/// \brief This code loops over BCs to pair photons using only EMCal + TVX (no central barrel, no collisions)
///
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt
///

#include <CommonConstants/MathConstants.h>
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/EventSelection.h"
#include "EMCALBase/Geometry.h"
#include "PWGJE/DataModel/EMCALClusters.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyBCs = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps, aod::Run3MatchedToBCSparse>;
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;

using SelectedUniqueClusters = soa::Filtered<aod::EMCALClusters>;             // Clusters from collisions with only one collision in the BC
using SelectedAmbiguousClusters = soa::Filtered<aod::EMCALAmbiguousClusters>; // Clusters from BCs with multiple collisions (no vertex assignment possible)

struct Photon { // Struct to store photons (unique and ambiguous clusters that passed the cuts)
  Photon(float eta, float phi, float energy) : eta(eta), phi(phi), e(energy), theta(2 * std::atan2(std::exp(-eta), 1)), px(e * std::sin(theta) * std::cos(phi)), py(e * std::sin(theta) * std::sin(phi)), pz(e * std::cos(theta)), pt(std::sqrt(px * px + py * py))
  {
    photon.SetPxPyPzE(px, py, pz, e);
  }

  TLorentzVector photon;
  float eta, phi, e;
  float theta;
  float px, py, pz, pt;
};

/// \brief returns if cluster is too close to edge of EMCal (using rotation background method only for EMCal!)
bool IsTooCloseToEdge(const int cellID, const int DistanceToBorder = 1, emcal::Geometry* emcalGeom = nullptr)
{
  if (DistanceToBorder <= 0) {
    return false;
  }
  if (cellID < 0) {
    return true;
  }

  int iBadCell = -1;

  // check distance to border in case the cell is okay
  auto [iSupMod, iMod, iPhi, iEta] = emcalGeom->GetCellIndex(cellID);
  auto [irow, icol] = emcalGeom->GetCellPhiEtaIndexInSModule(iSupMod, iMod, iPhi, iEta);

  // Check rows/phi
  int iRowLast = 24;
  if (emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::EMCAL_HALF) {
    iRowLast /= 2; // 2/3 sm case
  } else if (emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::EMCAL_THIRD) {
    iRowLast /= 3; // 1/3 sm case
  } else if (emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::DCAL_EXT) {
    iRowLast /= 3; // 1/3 sm case
  }

  if (irow < DistanceToBorder || (iRowLast - irow) <= DistanceToBorder) {
    iBadCell = 1;
  }

  if (iBadCell > 0) {
    return true;
  }
  return false;
}

bool isPhotonAccepted(Photon const& p, emcal::Geometry* emcalGeom = nullptr)
{
  int cellID = 0;
  try {
    cellID = emcalGeom->GetAbsCellIdFromEtaPhi(p.eta, p.phi);
    if (IsTooCloseToEdge(cellID, 1, emcalGeom))
      cellID = -1;
  } catch (o2::emcal::InvalidPositionException& e) {
    cellID = -1;
  }

  if (cellID == -1)
    return false;

  return true;
}

struct Meson {
  Meson(Photon p1, Photon p2) : p1(p1), p2(p2)
  {
    pMeson = p1.photon + p2.photon;
  }
  Photon p1, p2;
  TLorentzVector pMeson;

  float m() const { return pMeson.M(); }
  float pT() const { return pMeson.Pt(); }
  float openingAngle() const { return p1.photon.Angle(p2.photon.Vect()); }
};

struct EmcalGammaGammaBcWise {
  HistogramRegistry mHistManager{"EmcalGammaGammaBcWiseHistograms"};

  Configurable<int> cfgClusterDefinition{"cfgClusterDefinition", 13, "Clusterizer to be selected, e.g. 13 for kV3MostSplitLowSeed"};
  Configurable<float> cfgMinClusterEnergy{"cfgMinClusterEnergy", 0.7, "Minimum energy of selected clusters (GeV)"};
  Configurable<float> cfgMinM02{"cfgMinM02", 0.1, "Minimum M02 of selected clusters"};
  Configurable<float> cfgMaxM02{"cfgMaxM02", 0.7, "Maximum M02 of selected clusters"};
  Configurable<float> cfgMinTime{"cfgMinTime", -15, "Minimum time of selected clusters (ns)"};
  Configurable<float> cfgMaxTime{"cfgMaxTime", 15, "Maximum time of selected clusters (ns)"};
  Configurable<float> cfgMinOpenAngle{"cfgMinOpenAngle", 0.0202, "Minimum opening angle between photons"};

  Filter clusterDefinitionFilter = aod::emcalcluster::definition == cfgClusterDefinition;
  Filter energyFilter = aod::emcalcluster::energy > cfgMinClusterEnergy;
  Filter m02Filter = (aod::emcalcluster::nCells == 1 || (aod::emcalcluster::m02 > cfgMinM02 && aod::emcalcluster::m02 < cfgMaxM02));
  Filter timeFilter = (aod::emcalcluster::time > cfgMinTime && aod::emcalcluster::time < cfgMaxTime);

  std::vector<Photon> mPhotons;

  emcal::Geometry* emcalGeom;

  void init(InitContext const&)
  {
    emcalGeom = emcal::Geometry::GetInstanceFromRunNumber(300000);
    mHistManager.add("nBCs", "Number of BCs (with (k)TVX);#bf{TVX triggered};#bf{#it{N}_{BCs}}", HistType::kTH1F, {{3, -0.5, 2.5}});

    mHistManager.add("nCollisionsVsClusters", "Number of collisions vs number of clusters;N_{Collisions};N_{Clusters}", HistType::kTH2F, {{4, -0.5, 3.5}, {21, -0.5, 20.5}});

    mHistManager.add("clusterE", "Energy of cluster;#bf{#it{E} (GeV)};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, 0, 20}});
    mHistManager.add("clusterM02", "Shape of cluster;#bf{#it{M}_{02}};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, 0, 2}});
    mHistManager.add("clusterTime", "Time of cluster;#bf{#it{t} (ns)};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, -100, 100}});

    mHistManager.add("invMassVsPt", "Invariant mass and pT of meson candidates", HistType::kTH2F, {{400, 0., 0.8}, {200, 0., 20.}});
    mHistManager.add("invMassVsPtBackground", "Invariant mass and pT of background meson candidates", HistType::kTH2F, {{400, 0., 0.8}, {200, 0., 20.}});
  }

  PresliceUnsorted<MyCollisions> perFoundBC = aod::evsel::foundBCId;
  Preslice<aod::EMCALClusters> perCol = aod::emcalcluster::collisionId;
  Preslice<aod::EMCALAmbiguousClusters> perBC = aod::emcalcluster::bcId;

  void process(MyBCs const& bcs, MyCollisions const& collisions, SelectedUniqueClusters const& uClusters, SelectedAmbiguousClusters const& aClusters)
  {
    for (const auto& bc : bcs) {
      mHistManager.fill(HIST("nBCs"), 0.);
      if (!bc.selection_bit(aod::evsel::kIsTriggerTVX))
        continue;
      mHistManager.fill(HIST("nBCs"), 1.);
      if (!bc.alias_bit(kTVXinEMC))
        continue;
      mHistManager.fill(HIST("nBCs"), 2.);

      auto collisionsInFoundBC = collisions.sliceBy(perFoundBC, bc.globalIndex());

      if (collisionsInFoundBC.size() == 1) { // Unique
        auto clustersInCollision = uClusters.sliceBy(perCol, collisionsInFoundBC.begin().globalIndex());
        processClusters(clustersInCollision);
      } else { // Ambiguous
        auto clustersInBC = aClusters.sliceBy(perBC, bc.globalIndex());
        processClusters(clustersInBC);
      }

      mHistManager.fill(HIST("nCollisionsVsClusters"), collisionsInFoundBC.size(), mPhotons.size());

      processMesons();
    }
  }

  /// \brief Process EMCAL clusters (either ambigous or unique)
  template <typename Clusters>
  void processClusters(Clusters const& clusters)
  {
    mPhotons.clear();

    // loop over all clusters from accepted collision
    for (const auto& cluster : clusters) {

      mHistManager.fill(HIST("clusterE"), cluster.energy());
      mHistManager.fill(HIST("clusterM02"), cluster.m02());
      mHistManager.fill(HIST("clusterTime"), cluster.time());

      mPhotons.push_back(Photon(cluster.eta(), cluster.phi(), cluster.energy()));
    }
  }

  /// \brief Process meson candidates, calculate invariant mass and pT and fill histograms
  void processMesons()
  {
    if (mPhotons.size() < 2) // if less then 2 clusters are found, skip event
      return;

    // loop over all photon combinations and build meson candidates
    for (unsigned int ig1 = 0; ig1 < mPhotons.size(); ++ig1) {
      for (unsigned int ig2 = ig1 + 1; ig2 < mPhotons.size(); ++ig2) {
        // build meson from photons
        if (mPhotons[ig1].photon.Angle(mPhotons[ig2].photon.Vect()) < cfgMinOpenAngle)
          continue;
        Meson meson(mPhotons[ig1], mPhotons[ig2]);
        mHistManager.fill(HIST("invMassVsPt"), meson.m(), meson.pT());

        calculateBackground(meson, ig1, ig2); // calculate background candidates (rotation background)
      }
    }
  }

  /// \brief Calculate background (using rotation background method)
  void calculateBackground(const Meson& meson, const unsigned int ig1, const unsigned int ig2)
  {
    if (mPhotons.size() < 3) // if less than 3 clusters are present, skip event
      return;

    TVector3 lvRotationPion = (meson.pMeson).Vect(); // calculate rotation axis
    for (unsigned int ig3 = 0; ig3 < mPhotons.size(); ++ig3) {
      if (ig3 == ig1 || ig3 == ig2) // Skip if photon is one of the meson constituents
        continue;
      for (const unsigned int ig : {ig1, ig2}) {
        TLorentzVector lvRotationPhoton(mPhotons[ig].px, mPhotons[ig].py, mPhotons[ig].pz, mPhotons[ig].e);
        lvRotationPhoton.Rotate(constants::math::PIHalf, lvRotationPion);
        Photon rotPhoton(lvRotationPhoton.Eta(), lvRotationPhoton.Phi(), lvRotationPhoton.E());
        if (!isPhotonAccepted(rotPhoton, emcalGeom))
          continue;
        if (rotPhoton.photon.Angle(mPhotons[ig3].photon.Vect()) < cfgMinOpenAngle)
          continue;
        Meson mesonRotated(rotPhoton, mPhotons[ig3]);
        mHistManager.fill(HIST("invMassVsPtBackground"), mesonRotated.m(), mesonRotated.pT());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<EmcalGammaGammaBcWise>(cfgc)};
}
