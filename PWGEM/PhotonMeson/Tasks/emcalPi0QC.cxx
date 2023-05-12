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

#include "EMCALBase/Geometry.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"

#include "CommonDataFormat/InteractionRecord.h"

#include "TLorentzVector.h"
#include "TVector3.h"

// \struct Pi0QCTask
/// \brief Simple monitoring task for EMCal clusters
/// \author Joshua Koenig <joshua.konig@cern.ch>, Goethe University Frankfurt
/// \since 25.05.2022
///
/// This task is meant to be used for QC for the emcal using properties of the pi0
/// - ...
/// Simple event selection using the flag doEventSel is provided, which selects INT7 events if set to 1
/// For pilot beam data, instead of relying on the event selection, one can veto specific BC IDS using the flag
/// fDoVetoBCID.

using namespace o2::framework;
using namespace o2::framework::expressions;
using collisionEvSelIt = o2::aod::Collision;
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
using selectedCluster = o2::soa::Filtered<o2::aod::EMCALCluster>;
using selectedAmbiguousClusters = o2::soa::Filtered<o2::aod::EMCALAmbiguousClusters>;
using selectedAmbiguousCluster = o2::soa::Filtered<o2::aod::EMCALAmbiguousCluster>;

struct Photon {
  Photon(float eta_tmp, float phi_tmp, float energy_tmp, int clusid = 0)
  {
    eta = eta_tmp;
    phi = phi_tmp;
    energy = energy_tmp;
    theta = 2 * std::atan2(std::exp(-eta), 1);
    px = energy * std::sin(theta) * std::cos(phi);
    py = energy * std::sin(theta) * std::sin(phi);
    pz = energy * std::cos(theta);
    pt = std::sqrt(px * px + py * py);
    photon.SetPxPyPzE(px, py, pz, energy);
    id = clusid;
  }

  TLorentzVector photon;
  float pt;
  float px;
  float py;
  float pz;
  float eta;
  float phi;
  float energy;
  float theta;
  int id;
};

struct Meson {
  Meson(Photon p1, Photon p2) : pgamma1(p1),
                                pgamma2(p2)
  {
    pMeson = p1.photon + p2.photon;
  }
  Photon pgamma1;
  Photon pgamma2;
  TLorentzVector pMeson;

  float getMass() const { return pMeson.M(); }
  float getPt() const { return pMeson.Pt(); }
  float getOpeningAngle() const { return pgamma1.photon.Angle(pgamma2.photon.Vect()); }
};

struct EventMixVec {

  void AddEvent(std::vector<Photon> vecGamma)
  {
    if (vecEvtMix.size() < nEVtMixSize) {
      vecEvtMix.push_back(vecGamma);
    } else {
      vecEvtMix.erase(vecEvtMix.begin() + nEVtMixSize - 1);
      vecEvtMix.push_back(vecGamma);
    }
  }
  Photon* getPhoton(unsigned int iEvt, unsigned int iGamma)
  {
    if (vecEvtMix.size() >= iEvt)
      return nullptr;
    if (vecEvtMix[iEvt].size() >= iGamma)
      return nullptr;
    return &vecEvtMix[iEvt][iGamma];
  }

  std::vector<std::vector<Photon>> vecEvtMix;
  unsigned int nEVtMixSize = 20;
};

struct Pi0QCTask {
  HistogramRegistry mHistManager{"NeutralMesonHistograms"};
  o2::emcal::Geometry* mGeometry = nullptr;

  // configurable parameters
  // TODO adapt mDoEventSel switch to also allow selection of other triggers (e.g. EMC7)
  Configurable<bool> mDoEventSel{"doEventSel", 0, "demand kINT7"};
  Configurable<std::string> mVetoBCID{"vetoBCID", "", "BC ID(s) to be excluded, this should be used as an alternative to the event selection"};
  Configurable<std::string> mSelectBCID{"selectBCID", "all", "BC ID(s) to be included, this should be used as an alternative to the event selection"};
  Configurable<double> mVertexCut{"vertexCut", -1, "apply z-vertex cut with value in cm"};
  Configurable<int> mTimeMin{"TimeMinCut", -600, "apply min timing cut (in ns)"};
  Configurable<int> mTimeMax{"TimeMaxCut", 900, "apply min timing cut (in ns)"};
  Configurable<float> mClusterMinM02Cut{"MinM02Cut", 0.1, "apply min M02 cut"};
  Configurable<float> mClusterMaxM02Cut{"MaxM02Cut", 0.7, "apply max M02 cut"};
  Configurable<float> mMinEnergyCut{"MinEnergyCut", 0.7, "apply min cluster energy cut"};
  Configurable<int> mMinNCellsCut{"MinNCellsCut", 1, "apply min cluster number of cell cut"};
  Configurable<float> mMinOpenAngleCut{"OpeningAngleCut", 0.0202, "apply min opening angle cut"};
  Configurable<std::string> mClusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  std::vector<int> mVetoBCIDs;
  std::vector<int> mSelectBCIDs;

  // define cluster filter. It selects only those clusters which are of the type
  // specified in the string mClusterDefinition,e.g. kV3Default, which is V3 clusterizer with default
  // clusterization parameters
  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);

  // define container for photons
  std::vector<Photon> mPhotons;
  // define container for photons for each collision
  std::map<int, std::vector<Photon>> mapPhotons;

  // event mixing class
  EventMixVec evtMix;

  /// \brief Create output histograms and initialize geometry
  void init(InitContext const&)
  {
    // create histograms
    using o2HistType = HistType;
    using o2Axis = AxisSpec;

    // load geometry just in case we need it
    mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(300000);

    // create common axes
    LOG(info) << "Creating histograms";
    const o2Axis bcAxis{3501, -0.5, 3500.5};
    const o2Axis energyAxis{makeClusterBinning(), "#it{p}_{T} (GeV)"};

    // event properties
    mHistManager.add("eventsAll", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsEMCTrigg", "Number of EMC triggered events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsSelected", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventBCAll", "Bunch crossing ID of event (all events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventBCSelected", "Bunch crossing ID of event (selected events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventVertexZAll", "z-vertex of event (all events)", o2HistType::kTH1F, {{200, -20, 20}});
    mHistManager.add("eventVertexZSelected", "z-vertex of event (selected events)", o2HistType::kTH1F, {{200, -20, 20}});

    // cluster properties
    mHistManager.add("clusterE", "Energy of cluster", o2HistType::kTH1F, {energyAxis});
    mHistManager.add("clusterE_SimpleBinning", "Energy of cluster", o2HistType::kTH1F, {{400, 0, 100}});
    mHistManager.add("clusterTime", "Time of cluster", o2HistType::kTH1F, {{500, -250, 250}});
    mHistManager.add("clusterEtaPhi", "Eta and phi of cluster", o2HistType::kTH2F, {{100, -1, 1}, {100, 0, 2 * TMath::Pi()}});
    mHistManager.add("clusterM02", "M02 of cluster", o2HistType::kTH1F, {{400, 0, 5}});
    mHistManager.add("clusterM20", "M20 of cluster", o2HistType::kTH1F, {{400, 0, 2.5}});
    mHistManager.add("clusterNLM", "Number of local maxima of cluster", o2HistType::kTH1I, {{10, 0, 10}});
    mHistManager.add("clusterNCells", "Number of cells in cluster", o2HistType::kTH1I, {{50, 0, 50}});
    mHistManager.add("clusterDistanceToBadChannel", "Distance to bad channel", o2HistType::kTH1F, {{100, 0, 100}});

    // meson related histograms
    mHistManager.add("invMassVsPt", "invariant mass and pT of meson candidates", o2HistType::kTH2F, {{400, 0, 0.8}, {energyAxis}});
    mHistManager.add("invMassVsPtBackground", "invariant mass and pT of background meson candidates", o2HistType::kTH2F, {{400, 0, 0.8}, {energyAxis}});
    mHistManager.add("invMassVsPtMixedBackground", "invariant mass and pT of mixed background meson candidates", o2HistType::kTH2F, {{400, 0, 0.8}, {energyAxis}});

    if (mVetoBCID->length()) {
      std::stringstream parser(mVetoBCID.value);
      std::string token;
      int bcid;
      while (std::getline(parser, token, ',')) {
        bcid = std::stoi(token);
        LOG(info) << "Veto BCID " << bcid;
        mVetoBCIDs.push_back(bcid);
      }
    }
    if (mSelectBCID.value != "all") {
      std::stringstream parser(mSelectBCID.value);
      std::string token;
      int bcid;
      while (std::getline(parser, token, ',')) {
        bcid = std::stoi(token);
        LOG(info) << "Select BCID " << bcid;
        mSelectBCIDs.push_back(bcid);
      }
    }
  }
  /// \brief Process EMCAL clusters that are matched to a collisions

  // void processCollisions(collisionEvSelIt const& collision, selectedClusters const& clusters)
  void processCollisions(o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator const& collision, selectedClusters const& clusters)
  {
    // for(const auto & collision : theCollisions){
    mHistManager.fill(HIST("eventsAll"), 1);
    LOG(debug) << "processCollisions";
    // do event selection if mDoEventSel is specified
    // currently the event selection is hard coded to kINT7
    // but other selections are possible that are defined in TriggerAliases.h
    if (mDoEventSel && (!collision.alias()[kINT7])) {
      LOG(debug) << "Event not selected becaus it is not kINT7, skipping";
      return;
    }
    mHistManager.fill(HIST("eventVertexZAll"), collision.posZ());
    if (mVertexCut > 0 && std::abs(collision.posZ()) > mVertexCut) {
      LOG(debug) << "Event not selected because of z-vertex cut z= " << collision.posZ() << " > " << mVertexCut << " cm, skipping";
      return;
    }
    mHistManager.fill(HIST("eventsSelected"), 1);
    mHistManager.fill(HIST("eventVertexZSelected"), collision.posZ());

    ProcessClusters(clusters);
    ProcessMesons(clusters);
  }
  PROCESS_SWITCH(Pi0QCTask, processCollisions, "Process clusters from collision", false);

  /// \brief Process EMCAL clusters that are not matched to a collision
  /// This is not needed for most users
  void processAmbiguous(o2::aod::BCs::iterator const& bc, selectedAmbiguousClusters const& clusters)
  {
    LOG(debug) << "processAmbiguous";
    // TODO: remove this loop and put it in separate process function that only takes care of ambiguous clusters
    o2::InteractionRecord eventIR;
    eventIR.setFromLong(bc.globalBC());
    mHistManager.fill(HIST("eventBCAll"), eventIR.bc);
    if (std::find(mVetoBCIDs.begin(), mVetoBCIDs.end(), eventIR.bc) != mVetoBCIDs.end()) {
      LOG(info) << "Event rejected because of veto BCID " << eventIR.bc;
      return;
    }
    if (mSelectBCIDs.size() && (std::find(mSelectBCIDs.begin(), mSelectBCIDs.end(), eventIR.bc) == mSelectBCIDs.end())) {
      return;
    }
    mHistManager.fill(HIST("eventBCSelected"), eventIR.bc);

    ProcessAmbigousClusters(clusters);
    ProcessMesons(clusters);
  }
  PROCESS_SWITCH(Pi0QCTask, processAmbiguous, "Process Ambiguous clusters", false);

  /// \brief Process EMCAL clusters that are matched to a collisions
  template <typename Clusters>
  void ProcessClusters(Clusters const& clusters)
  {
    LOG(debug) << "ProcessClusters";
    // clear photon vector
    mPhotons.clear();
    mapPhotons.clear();

    int globalCollID = -1000;

    // loop over all clusters from accepted collision
    // auto eventClusters = clusters.select(o2::aod::emcalcluster::bcId == theCollision.bc().globalBC());
    for (const auto& cluster : clusters) {

      // o2::InteractionRecord eventIR;
      auto collID = cluster.collisionId();
      if (globalCollID == -1000)
        globalCollID = collID;

      if (globalCollID != collID) {
        LOG(info) << "Something went wrong with the collision ID";
      }

      // fill histograms of cluster properties
      // in this implementation the cluster properties are directly
      // loaded from the flat table, in the future one should
      // consider using the AnalysisCluster object to work with
      // after loading.
      LOG(debug) << "Cluster energy: " << cluster.energy();
      LOG(debug) << "Cluster time: " << cluster.time();
      LOG(debug) << "Cluster M02: " << cluster.m02();
      mHistManager.fill(HIST("clusterE"), cluster.energy());
      mHistManager.fill(HIST("clusterTime"), cluster.time());
      mHistManager.fill(HIST("clusterE_SimpleBinning"), cluster.energy());
      mHistManager.fill(HIST("clusterEtaPhi"), cluster.eta(), cluster.phi());
      mHistManager.fill(HIST("clusterM02"), cluster.m02());
      mHistManager.fill(HIST("clusterM20"), cluster.m20());
      mHistManager.fill(HIST("clusterNLM"), cluster.nlm());
      mHistManager.fill(HIST("clusterNCells"), cluster.nCells());
      mHistManager.fill(HIST("clusterDistanceToBadChannel"), cluster.distanceToBadChannel());

      // apply basic cluster cuts
      if (cluster.energy() < mMinEnergyCut) {
        LOG(debug) << "Cluster rejected because of energy cut";
        continue;
      }
      if (cluster.nCells() <= mMinNCellsCut) {
        LOG(debug) << "Cluster rejected because of nCells cut";
        continue;
      }
      if (cluster.m02() < mClusterMinM02Cut || cluster.m02() > mClusterMaxM02Cut) {
        LOG(debug) << "Cluster rejected because of m02 cut";
        continue;
      }
      if (cluster.time() < mTimeMin || cluster.time() > mTimeMax) {
        LOG(debug) << "Cluster rejected because of time cut";
        continue;
      }

      // put clusters in photon vector
      mPhotons.push_back(Photon(cluster.eta(), cluster.phi(), cluster.energy(), cluster.id()));
    }
  }

  /// \brief Process EMCAL clusters that are matched to a collisions
  template <typename Clusters>
  void ProcessAmbigousClusters(Clusters const& clusters)
  {
    LOG(debug) << "ProcessClusters";
    // clear photon vector
    mPhotons.clear();

    // loop over all clusters from accepted collision
    for (const auto& cluster : clusters) {

      // fill histograms of cluster properties
      // in this implementation the cluster properties are directly
      // loaded from the flat table, in the future one should
      // consider using the AnalysisCluster object to work with
      // after loading.
      LOG(debug) << "Cluster energy: " << cluster.energy();
      LOG(debug) << "Cluster time: " << cluster.time();
      LOG(debug) << "Cluster M02: " << cluster.m02();
      mHistManager.fill(HIST("clusterE"), cluster.energy());
      mHistManager.fill(HIST("clusterTime"), cluster.time());
      mHistManager.fill(HIST("clusterE_SimpleBinning"), cluster.energy());
      mHistManager.fill(HIST("clusterEtaPhi"), cluster.eta(), cluster.phi());
      mHistManager.fill(HIST("clusterM02"), cluster.m02());
      mHistManager.fill(HIST("clusterM20"), cluster.m20());
      mHistManager.fill(HIST("clusterNLM"), cluster.nlm());
      mHistManager.fill(HIST("clusterNCells"), cluster.nCells());
      mHistManager.fill(HIST("clusterDistanceToBadChannel"), cluster.distanceToBadChannel());

      // apply basic cluster cuts
      if (cluster.energy() < mMinEnergyCut) {
        LOG(debug) << "Cluster rejected because of energy cut";
        continue;
      }
      if (cluster.nCells() <= mMinNCellsCut) {
        LOG(debug) << "Cluster rejected because of nCells cut";
        continue;
      }
      if (cluster.m02() < mClusterMinM02Cut || cluster.m02() > mClusterMaxM02Cut) {
        LOG(debug) << "Cluster rejected because of m02 cut";
        continue;
      }
      if (cluster.time() < mTimeMin || cluster.time() > mTimeMax) {
        LOG(debug) << "Cluster rejected because of time cut";
        continue;
      }

      // put clusters in photon vector
      mPhotons.push_back(Photon(cluster.eta(), cluster.phi(), cluster.energy(), cluster.id()));
    }
  }

  /// \brief Process meson candidates, calculate invariant mass and pT and fill histograms
  template <typename Clusters>
  void ProcessMesons(Clusters const& clusters)
  {
    LOG(debug) << "ProcessMesons " << mPhotons.size();

    mHistManager.fill(HIST("eventsEMCTrigg"), 1);

    // if less then 2 clusters are found, skip event
    if (mPhotons.size() < 2)
      return;

    // loop over all photon combinations and build meson candidates
    for (unsigned int ig1 = 0; ig1 < mPhotons.size(); ++ig1) {
      for (unsigned int ig2 = ig1 + 1; ig2 < mPhotons.size(); ++ig2) {

        // build meson from photons
        Meson meson(mPhotons[ig1], mPhotons[ig2]);
        if (meson.getOpeningAngle() > mMinOpenAngleCut) {
          mHistManager.fill(HIST("invMassVsPt"), meson.getMass(), meson.getPt());
        }

        // calculate background candidates (rotation background)
        CalculateBackground(meson, ig1, ig2);
      }
      CalculateMixedBack(mPhotons[ig1]);
    }

    evtMix.AddEvent(mPhotons);
  }

  /// \brief Calculate background (using rotation background method)
  void CalculateBackground(const Meson& meson, unsigned int ig1, unsigned int ig2)
  {
    // if less than 3 clusters are present, skip event
    if (mPhotons.size() < 3) {
      return;
    }
    const double rotationAngle = M_PI / 2.0; //0.78539816339; // rotaion angle 90°

    TLorentzVector lvRotationPhoton1; // photon candidates which get rotated
    TLorentzVector lvRotationPhoton2; // photon candidates which get rotated
    TVector3 lvRotationPion;          // rotation axis
    for (unsigned int ig3 = 0; ig3 < mPhotons.size(); ++ig3) {
      // continue if photons are identical
      if (ig3 == ig1 || ig3 == ig2) {
        continue;
      }
      // calculate rotation axis
      lvRotationPion = (meson.pMeson).Vect();

      // initialize photons for rotation
      lvRotationPhoton1.SetPxPyPzE(mPhotons[ig1].px, mPhotons[ig1].py, mPhotons[ig1].pz, mPhotons[ig1].energy);
      lvRotationPhoton2.SetPxPyPzE(mPhotons[ig2].px, mPhotons[ig2].py, mPhotons[ig2].pz, mPhotons[ig2].energy);

      // rotate photons around rotation axis
      lvRotationPhoton1.Rotate(rotationAngle, lvRotationPion);
      lvRotationPhoton2.Rotate(rotationAngle, lvRotationPion);

      // initialize Photon objects for rotated photons
      Photon rotPhoton1(lvRotationPhoton1.Eta(), lvRotationPhoton1.Phi(), lvRotationPhoton1.E(), mPhotons[ig1].id);
      Photon rotPhoton2(lvRotationPhoton2.Eta(), lvRotationPhoton2.Phi(), lvRotationPhoton2.E(), mPhotons[ig2].id);

      // build meson from rotated photons
      Meson mesonRotated1(rotPhoton1, mPhotons[ig3]);
      Meson mesonRotated2(rotPhoton2, mPhotons[ig3]);

      // Fill histograms
      if (mesonRotated1.getOpeningAngle() > mMinOpenAngleCut) {
        mHistManager.fill(HIST("invMassVsPtBackground"), mesonRotated1.getMass(), mesonRotated1.getPt());
      }
      if (mesonRotated2.getOpeningAngle() > mMinOpenAngleCut) {
        mHistManager.fill(HIST("invMassVsPtBackground"), mesonRotated2.getMass(), mesonRotated2.getPt());
      }
    }
  }

  void CalculateMixedBack(Photon gamma)
  {
    for (unsigned int i = 0; i < evtMix.vecEvtMix.size(); ++i) {
      for (unsigned int ig1 = 0; ig1 < evtMix.vecEvtMix[i].size(); ++ig1) {
        Meson meson(gamma, evtMix.vecEvtMix[i][ig1]);
        if (meson.getOpeningAngle() > mMinOpenAngleCut) {
          mHistManager.fill(HIST("invMassVsPtMixedBackground"), meson.getMass(), meson.getPt());
        }
      }
    }
  }

  /// \brief Create binning for cluster energy/pT axis (variable bin size)
  /// direct port from binning often used in AliPhysics for debugging
  /// \return vector with bin limits
  std::vector<double> makeClusterBinning() const
  {

    std::vector<double> result;
    int nBinsPt = 179;
    double maxPt = 60;
    for (Int_t i = 0; i < nBinsPt + 1; i++) {
      if (i < 100) {
        result.emplace_back(0.10 * i);
      } else if (i < 140) {
        result.emplace_back(10. + 0.25 * (i - 100));
      } else if (i < 180) {
        result.emplace_back(20. + 1.00 * (i - 140));
      } else {
        result.emplace_back(maxPt);
      }
    }
    return result;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<Pi0QCTask>(cfgc, TaskName{"EMCPi0QCTask"}, SetDefaultProcesses{{{"processCollisions", true}, {"processAmbiguous", false}}}),
    adaptAnalysisTask<Pi0QCTask>(cfgc, TaskName{"EMCPi0QCTaskAmbiguous"}, SetDefaultProcesses{{{"processCollisions", false}, {"processAmbiguous", true}}})};
  return workflow;
}
