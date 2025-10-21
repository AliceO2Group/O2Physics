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

/// \file emcalPi0Qc.cxx
/// \brief Simple monitoring task for EMCal clusters
/// \author Joshua Koenig <joshua.konig@cern.ch>, Goethe University Frankfurt
/// \struct EmcalPi0Qc
/// \since 25.05.2022
///
/// This task is meant to be used for QC for the emcal using properties of the pi0
/// - ...
/// Simple event selection using the flag doEventSel is provided, which selects INT7 events if set to 1
/// For pilot beam data, instead of relying on the event selection, one can veto specific BC IDS using the flag
/// fDoVetoBCID.

#include "PWGJE/DataModel/EMCALClusterDefinition.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/MathConstants.h>
#include <CommonDataFormat/InteractionRecord.h>
#include <EMCALBase/Geometry.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TVector3.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2::framework;
using namespace o2::framework::expressions;
using SelectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
using MyCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::EMCALMatchedCollisions>;
using MyBCs = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;
using SelectedAmbiguousClusters = o2::soa::Filtered<o2::aod::EMCALAmbiguousClusters>;

struct Photon {
  Photon(float eta_tmp, float phi_tmp, float energy_tmp, int clusid = 0)
  {
    eta = eta_tmp;
    phi = phi_tmp;
    onDCal = (phi < 6 && phi > 4);
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
  bool onDCal; // Checks whether photon is in phi region of the DCal, otherwise: EMCal
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

  void addEvent(std::vector<Photon> vecGamma)
  {
    if (vecEvtMix.size() < nEVtMixSize) {
      vecEvtMix.push_back(vecGamma);
    } else {
      vecEvtMix.erase(vecEvtMix.begin() + nEVtMixSize - 1);
      vecEvtMix.push_back(vecGamma);
    }
  }

  std::vector<std::vector<Photon>> vecEvtMix;
  unsigned int nEVtMixSize = 20;
};

struct EmcalPi0Qc {
  HistogramRegistry mHistManager{"NeutralMesonHistograms"};
  o2::emcal::Geometry* mGeometry = nullptr;

  Filter emccellfilter = o2::aod::calo::caloType == 1;

  // configurable parameters
  // TODO adapt mDoEventSel switch to also allow selection of other triggers (e.g. EMC7)
  Configurable<bool> mDoEventSel{"mDoEventSel", 0, "demand kINT7"};
  Configurable<bool> mRequireCaloReadout{"mRequireCaloReadout", 0, "require kTVXinEMC"};
  Configurable<bool> mRequireEMCalCells{"mRequireEMCalCells", 0, "require at least one EMC cell in each collision"};
  Configurable<std::string> mVetoBCID{"mVetoBCID", "", "BC ID(s) to be excluded, this should be used as an alternative to the event selection"};
  Configurable<std::string> mSelectBCID{"mSelectBCID", "all", "BC ID(s) to be included, this should be used as an alternative to the event selection"};
  Configurable<double> mVertexCut{"mVertexCut", -1, "apply z-vertex cut with value in cm"};
  Configurable<int> mTimeMin{"mTimeMin", -600, "apply min timing cut (in ns)"};
  Configurable<int> mTimeMax{"mTimeMax", 900, "apply min timing cut (in ns)"};
  Configurable<float> mClusterMinM02Cut{"mClusterMinM02Cut", 0.1, "apply min M02 cut"};
  Configurable<float> mClusterMaxM02Cut{"mClusterMaxM02Cut", 0.7, "apply max M02 cut"};
  Configurable<float> mMinEnergyCut{"mMinEnergyCut", 0.7, "apply min cluster energy cut"};
  Configurable<int> mMinNCellsCut{"mMinNCellsCut", 1, "apply min cluster number of cell cut"};
  Configurable<float> mMinOpenAngleCut{"mMinOpenAngleCut", 0.0202, "apply min opening angle cut"};
  Configurable<std::string> mClusterDefinition{"mClusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<bool> mSplitEMCalDCal{"mSplitEMCalDCal", 0, "Create and fill inv mass histograms for photons on EMCal and DCal individually"};
  std::vector<int> mVetoBCIDs;
  std::vector<int> mSelectBCIDs;

  ConfigurableAxis pTBinning{"pTBinning", {500, 0.0f, 50.0f}, "Binning used along pT axis for inv mass histograms"};
  ConfigurableAxis invmassBinning{"invmassBinning", {400, 0.0f, 0.8f}, "Binning used for inv mass axis in inv mass - pT histograms"};

  // define cluster filter. It selects only those clusters which are of the type
  // specified in the string mClusterDefinition,e.g. kV3Default, which is V3 clusterizer with default
  // clusterization parameters
  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);

  // define container for photons
  std::vector<Photon> mPhotons;

  // event mixing class
  EventMixVec evtMix;

  /// \brief Create output histograms and initialize geometry
  void init(InitContext const&)
  {
    // load geometry just in case we need it
    mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(300000);

    // create common axes
    LOG(info) << "Creating histograms";
    const AxisSpec bcAxis{3501, -0.5, 3500.5};
    const AxisSpec energyAxis{makeClusterBinning(), "#it{E} (GeV)"};

    mHistManager.add("events", "events;;#it{count}", HistType::kTH1F, {{6, 0.5, 6.5}});
    auto heventType = mHistManager.get<TH1>(HIST("events"));
    heventType->GetXaxis()->SetBinLabel(1, "All events");
    heventType->GetXaxis()->SetBinLabel(2, "sel8 + readout");
    heventType->GetXaxis()->SetBinLabel(3, "1+ Contributor");
    heventType->GetXaxis()->SetBinLabel(4, "z<10cm");
    heventType->GetXaxis()->SetBinLabel(5, "unique col");
    heventType->GetXaxis()->SetBinLabel(6, "EMCAL cell>0");
    mHistManager.add("eventBCAll", "Bunch crossing ID of event (all events)", HistType::kTH1F, {bcAxis});
    mHistManager.add("eventBCSelected", "Bunch crossing ID of event (selected events)", HistType::kTH1F, {bcAxis});
    mHistManager.add("eventVertexZAll", "z-vertex of event (all events)", HistType::kTH1F, {{200, -20, 20}});
    mHistManager.add("eventVertexZSelected", "z-vertex of event (selected events)", HistType::kTH1F, {{200, -20, 20}});

    // cluster properties
    for (const bool& iBeforeCuts : {false, true}) {
      const char* clusterDirectory = iBeforeCuts ? "ClustersBeforeCuts" : "ClustersAfterCuts";
      mHistManager.add(Form("%s/clusterE", clusterDirectory), "Energy of cluster", HistType::kTH1F, {energyAxis});
      mHistManager.add(Form("%s/clusterE_SimpleBinning", clusterDirectory), "Energy of cluster", HistType::kTH1F, {{400, 0, 100, "#it{E} (GeV)"}});
      mHistManager.add(Form("%s/clusterTime", clusterDirectory), "Time of cluster", HistType::kTH1F, {{500, -250, 250, "#it{t}_{cls} (ns)"}});
      mHistManager.add(Form("%s/clusterEtaPhi", clusterDirectory), "Eta and phi of cluster", HistType::kTH2F, {{100, -1, 1, "#eta"}, {100, 0, o2::constants::math::TwoPI, "#phi"}});
      mHistManager.add(Form("%s/clusterM02", clusterDirectory), "M02 of cluster", HistType::kTH1F, {{400, 0, 5, "#it{M}_{02}"}});
      mHistManager.add(Form("%s/clusterM20", clusterDirectory), "M20 of cluster", HistType::kTH1F, {{400, 0, 2.5, "#it{M}_{20}"}});
      mHistManager.add(Form("%s/clusterNLM", clusterDirectory), "Number of local maxima of cluster", HistType::kTH1I, {{10, 0, 10, "#it{N}_{local maxima}"}});
      mHistManager.add(Form("%s/clusterNCells", clusterDirectory), "Number of cells in cluster", HistType::kTH1I, {{50, 0, 50, "#it{N}_{cells}"}});
      mHistManager.add(Form("%s/clusterDistanceToBadChannel", clusterDirectory), "Distance to bad channel", HistType::kTH1F, {{100, 0, 100, "#it{d}"}});
    }

    // meson related histograms
    mHistManager.add("invMassVsPt", "invariant mass and pT of meson candidates", HistType::kTH2F, {invmassBinning, pTBinning});
    mHistManager.add("invMassVsPtBackground", "invariant mass and pT of background meson candidates", HistType::kTH2F, {invmassBinning, pTBinning});
    mHistManager.add("invMassVsPtMixedBackground", "invariant mass and pT of mixed background meson candidates", HistType::kTH2F, {invmassBinning, pTBinning});

    if (mSplitEMCalDCal) {
      mHistManager.add("invMassVsPt_EMCal", "invariant mass and pT of meson candidates with both clusters on EMCal", HistType::kTH2F, {invmassBinning, pTBinning});
      mHistManager.add("invMassVsPtBackground_EMCal", "invariant mass and pT of background meson candidates with both clusters on EMCal", HistType::kTH2F, {invmassBinning, pTBinning});
      mHistManager.add("invMassVsPtMixedBackground_EMCal", "invariant mass and pT of mixed background meson candidates with both clusters on EMCal", HistType::kTH2F, {invmassBinning, pTBinning});
      mHistManager.add("invMassVsPt_DCal", "invariant mass and pT of meson candidates with both clusters on DCal", HistType::kTH2F, {invmassBinning, pTBinning});
      mHistManager.add("invMassVsPtBackground_DCal", "invariant mass and pT of background meson candidates with both clusters on DCal", HistType::kTH2F, {invmassBinning, pTBinning});
      mHistManager.add("invMassVsPtMixedBackground_DCal", "invariant mass and pT of mixed background meson candidates with both clusters on DCal", HistType::kTH2F, {invmassBinning, pTBinning});
    }

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

  PresliceUnsorted<SelectedClusters> perCollision = o2::aod::emcalcluster::collisionId;

  /// \brief Process EMCAL clusters that are matched to a collisions
  void processCollision(MyBCs const&, MyCollisions const& collisions, SelectedClusters const& clusters, o2::soa::Filtered<o2::aod::Calos> const& cells)
  {
    std::unordered_map<uint64_t, int> cellGlobalBCs;
    // Build map of number of cells for corrected BCs using global BCs
    // used later in the determination whether a BC has EMC cell content (for speed reason)
    for (const auto& cell : cells) {
      auto globalbcid = cell.bc_as<MyBCs>().globalBC();
      auto found = cellGlobalBCs.find(globalbcid);
      if (found != cellGlobalBCs.end()) {
        found->second++;
      } else {
        cellGlobalBCs.insert(std::pair<uint64_t, int>(globalbcid, 1));
      }
    }

    for (const auto& collision : collisions) {
      mHistManager.fill(HIST("events"), 1); // Fill "All events" bin of event histogram

      if (mDoEventSel && (!collision.sel8() || (mRequireCaloReadout && !collision.alias_bit(kTVXinEMC)))) { // Check sel8 and whether EMC was read out
        continue;
      }
      mHistManager.fill(HIST("events"), 2); // Fill sel8 + readout

      if (mDoEventSel && collision.numContrib() < 0.5) { // Skip collisions without contributors
        continue;
      }
      mHistManager.fill(HIST("events"), 3); // Fill >1 vtx contr. bin of event histogram

      mHistManager.fill(HIST("eventVertexZAll"), collision.posZ());
      if (mVertexCut > 0 && std::abs(collision.posZ()) > mVertexCut) {
        continue;
      }
      mHistManager.fill(HIST("events"), 4); // Fill z-Vertex selected bin of event histogram
      mHistManager.fill(HIST("eventVertexZSelected"), collision.posZ());

      if (mDoEventSel && collision.ambiguous()) { // Skip ambiguous collisions (those that are in BCs including multiple collisions)
        continue;
      }
      mHistManager.fill(HIST("events"), 5); // Fill "One collision in BC" bin of event histogram

      if (mDoEventSel) {
        auto found = cellGlobalBCs.find(collision.foundBC_as<MyBCs>().globalBC());
        if (mRequireEMCalCells && (found == cellGlobalBCs.end() || found->second == 0)) { // Skip collisions without any readout EMCal cells
          continue;
        }
      }
      mHistManager.fill(HIST("events"), 6); // Fill at least one non0 cell in EMCal of event histogram (Selected)

      auto clustersPerColl = clusters.sliceBy(perCollision, collision.collisionId());
      processClusters(clustersPerColl);
      processMesons();
    }
  }
  PROCESS_SWITCH(EmcalPi0Qc, processCollision, "Process clusters from collision", false);

  /// \brief Process EMCAL clusters that are not matched to a collision
  /// This is not needed for most users
  void processAmbiguous(o2::aod::BCs::iterator const& bc, SelectedAmbiguousClusters const& clusters)
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

    processAmbiguousClusters(clusters);
    processMesons();
  }
  PROCESS_SWITCH(EmcalPi0Qc, processAmbiguous, "Process Ambiguous clusters", false);

  /// \brief Process EMCAL clusters that are matched to a collisions
  template <typename Clusters>
  void processClusters(Clusters const& clusters)
  {
    LOG(debug) << "processClusters";
    // clear photon vector
    mPhotons.clear();

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

      fillClusterQAHistos<decltype(cluster), 0>(cluster);

      if (clusterRejectedByCut(cluster))
        continue;

      fillClusterQAHistos<decltype(cluster), 1>(cluster);

      // put clusters in photon vector
      mPhotons.push_back(Photon(cluster.eta(), cluster.phi(), cluster.energy(), cluster.id()));
    }
  }

  /// \brief Process EMCAL clusters that are not matched to a collisions
  template <typename Clusters>
  void processAmbiguousClusters(Clusters const& clusters)
  {
    LOG(debug) << "processClusters";
    // clear photon vector
    mPhotons.clear();

    // loop over all clusters from accepted collision
    for (const auto& cluster : clusters) {

      fillClusterQAHistos<decltype(cluster), 0>(cluster);

      if (clusterRejectedByCut(cluster))
        continue;

      fillClusterQAHistos<decltype(cluster), 1>(cluster);

      // put clusters in photon vector
      mPhotons.push_back(Photon(cluster.eta(), cluster.phi(), cluster.energy(), cluster.id()));
    }
  }

  /// \brief Fills the standard QA histograms for a given cluster
  template <typename Cluster, int BeforeCuts>
  void fillClusterQAHistos(Cluster const& cluster)
  {
    // In this implementation the cluster properties are directly loaded from the flat table,
    // in the future one should consider using the AnalysisCluster object to work with after loading.
    static constexpr std::string_view kClusterQAHistEnergy[2] = {"ClustersBeforeCuts/clusterE", "ClustersAfterCuts/clusterE"};
    static constexpr std::string_view kClusterQAHistEnergySimpleBinning[2] = {"ClustersBeforeCuts/clusterE_SimpleBinning", "ClustersAfterCuts/clusterE_SimpleBinning"};
    static constexpr std::string_view kClusterQAHistTime[2] = {"ClustersBeforeCuts/clusterTime", "ClustersAfterCuts/clusterTime"};
    static constexpr std::string_view kClusterQAHistEtaPhi[2] = {"ClustersBeforeCuts/clusterEtaPhi", "ClustersAfterCuts/clusterEtaPhi"};
    static constexpr std::string_view kClusterQAHistM02[2] = {"ClustersBeforeCuts/clusterM02", "ClustersAfterCuts/clusterM02"};
    static constexpr std::string_view kClusterQAHistM20[2] = {"ClustersBeforeCuts/clusterM20", "ClustersAfterCuts/clusterM20"};
    static constexpr std::string_view kClusterQAHistNLM[2] = {"ClustersBeforeCuts/clusterNLM", "ClustersAfterCuts/clusterNLM"};
    static constexpr std::string_view kClusterQAHistNCells[2] = {"ClustersBeforeCuts/clusterNCells", "ClustersAfterCuts/clusterNCells"};
    static constexpr std::string_view kClusterQAHistDistanceToBadChannel[2] = {"ClustersBeforeCuts/clusterDistanceToBadChannel", "ClustersAfterCuts/clusterDistanceToBadChannel"};
    mHistManager.fill(HIST(kClusterQAHistEnergy[BeforeCuts]), cluster.energy());
    mHistManager.fill(HIST(kClusterQAHistEnergySimpleBinning[BeforeCuts]), cluster.energy());
    mHistManager.fill(HIST(kClusterQAHistTime[BeforeCuts]), cluster.time());
    mHistManager.fill(HIST(kClusterQAHistEtaPhi[BeforeCuts]), cluster.eta(), cluster.phi());
    mHistManager.fill(HIST(kClusterQAHistM02[BeforeCuts]), cluster.m02());
    mHistManager.fill(HIST(kClusterQAHistM20[BeforeCuts]), cluster.m20());
    mHistManager.fill(HIST(kClusterQAHistNLM[BeforeCuts]), cluster.nlm());
    mHistManager.fill(HIST(kClusterQAHistNCells[BeforeCuts]), cluster.nCells());
    mHistManager.fill(HIST(kClusterQAHistDistanceToBadChannel[BeforeCuts]), cluster.distanceToBadChannel());
  }

  /// \brief Return a boolean that states, whether a cluster should be rejected by the applied cluster cuts
  template <typename Cluster>
  bool clusterRejectedByCut(Cluster const& cluster)
  {
    // apply basic cluster cuts
    if (cluster.energy() < mMinEnergyCut) {
      LOG(debug) << "Cluster rejected because of energy cut";
      return true;
    }
    if (cluster.nCells() < mMinNCellsCut) {
      LOG(debug) << "Cluster rejected because of nCells cut";
      return true;
    }
    // Only apply M02 cut when cluster contains more than one cell
    if (cluster.nCells() > 1) {
      if (cluster.m02() < mClusterMinM02Cut || cluster.m02() > mClusterMaxM02Cut) {
        LOG(debug) << "Cluster rejected because of m02 cut";
        return true;
      }
    }
    if (cluster.time() < mTimeMin || cluster.time() > mTimeMax) {
      LOG(debug) << "Cluster rejected because of time cut";
      return true;
    }
    return false;
  }

  /// \brief Process meson candidates, calculate invariant mass and pT and fill histograms
  void processMesons()
  {
    LOG(debug) << "processMesons " << mPhotons.size();

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

          if (mSplitEMCalDCal) {
            if (!mPhotons[ig1].onDCal && !mPhotons[ig2].onDCal) {
              mHistManager.fill(HIST("invMassVsPt_EMCal"), meson.getMass(), meson.getPt());
            } else if (mPhotons[ig1].onDCal && mPhotons[ig2].onDCal) {
              mHistManager.fill(HIST("invMassVsPt_DCal"), meson.getMass(), meson.getPt());
            }
          }
        }

        // calculate background candidates (rotation background)
        calculateBackground(meson, ig1, ig2);
      }
      calculateMixedBack(mPhotons[ig1]);
    }

    evtMix.addEvent(mPhotons);
  }

  /// \brief Calculate background (using rotation background method)
  void calculateBackground(const Meson& meson, unsigned int ig1, unsigned int ig2)
  {
    // if less than 3 clusters are present, skip event
    if (mPhotons.size() < 3) {
      return;
    }
    const double rotationAngle = o2::constants::math::PIHalf; // 0.78539816339; // rotaion angle 90°

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
        if (mSplitEMCalDCal) {
          if (!mPhotons[ig1].onDCal && !mPhotons[ig2].onDCal && !mPhotons[ig3].onDCal) {
            mHistManager.fill(HIST("invMassVsPtBackground_EMCal"), mesonRotated1.getMass(), mesonRotated1.getPt());
          } else if (mPhotons[ig1].onDCal && mPhotons[ig2].onDCal && mPhotons[ig3].onDCal) {
            mHistManager.fill(HIST("invMassVsPtBackground_DCal"), mesonRotated1.getMass(), mesonRotated1.getPt());
          }
        }
      }
      if (mesonRotated2.getOpeningAngle() > mMinOpenAngleCut) {
        mHistManager.fill(HIST("invMassVsPtBackground"), mesonRotated2.getMass(), mesonRotated2.getPt());
        if (mSplitEMCalDCal) {
          if (!mPhotons[ig1].onDCal && !mPhotons[ig2].onDCal && !mPhotons[ig3].onDCal) {
            mHistManager.fill(HIST("invMassVsPtBackground_EMCal"), mesonRotated2.getMass(), mesonRotated2.getPt());
          } else if (mPhotons[ig1].onDCal && mPhotons[ig2].onDCal && mPhotons[ig3].onDCal) {
            mHistManager.fill(HIST("invMassVsPtBackground_DCal"), mesonRotated2.getMass(), mesonRotated2.getPt());
          }
        }
      }
    }
  }

  void calculateMixedBack(Photon gamma)
  {
    for (unsigned int i = 0; i < evtMix.vecEvtMix.size(); ++i) {
      for (unsigned int ig1 = 0; ig1 < evtMix.vecEvtMix[i].size(); ++ig1) {
        Meson meson(gamma, evtMix.vecEvtMix[i][ig1]);
        if (meson.getOpeningAngle() > mMinOpenAngleCut) {
          mHistManager.fill(HIST("invMassVsPtMixedBackground"), meson.getMass(), meson.getPt());
          if (mSplitEMCalDCal) {
            if (!gamma.onDCal && !evtMix.vecEvtMix[i][ig1].onDCal) {
              mHistManager.fill(HIST("invMassVsPtMixedBackground_EMCal"), meson.getMass(), meson.getPt());
            } else if (gamma.onDCal && evtMix.vecEvtMix[i][ig1].onDCal) {
              mHistManager.fill(HIST("invMassVsPtMixedBackground_DCal"), meson.getMass(), meson.getPt());
            }
          }
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
    for (int i = 0; i < nBinsPt + 1; i++) {
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
    adaptAnalysisTask<EmcalPi0Qc>(cfgc, TaskName{"EmcalPi0QcAssociate"}, SetDefaultProcesses{{{"processCollision", true}, {"processAmbiguous", false}}}),  // o2-linter: disable=name/o2-task (adapted multiple times)
    adaptAnalysisTask<EmcalPi0Qc>(cfgc, TaskName{"EmcalPi0QcAmbiguous"}, SetDefaultProcesses{{{"processCollision", false}, {"processAmbiguous", true}}})}; // o2-linter: disable=name/o2-task (adapted multiple times)
  return workflow;
}
