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
#include <numeric>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"

#include "CommonDataFormat/InteractionRecord.h"

// \struct ClusterMonitor
/// \brief Simple monitoring task for EMCal clusters
/// \author Florian Jonas <florian.jonas@cern.ch>, Oak Ridge National Laoratory
/// \since 30.03.2022
///
/// This task is meant to be used for monitoring EMCal clusters, allowing to track simple cluster
/// properties, such as:
/// - cluster energy
/// - cluster position
/// - cluster time
/// - cluster shape
/// - ...
/// Simple event selection using the flag doEventSel is provided, which selects INT7 events if set to 1
/// For pilot beam data, instead of relying on the event selection, one can veto specific BC IDS using the flag
/// fDoVetoBCID.
using namespace o2::framework;
using namespace o2::framework::expressions;
using collisionEvSelIt = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator;
using bcEvSelIt = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>::iterator;
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
using selectedAmbiguousClusters = o2::soa::Filtered<o2::aod::EMCALAmbiguousClusters>;
struct ClusterMonitor {
  HistogramRegistry mHistManager{"ClusterMonitorHistograms"};
  o2::emcal::Geometry* mGeometry = nullptr;

  Preslice<o2::aod::EMCALClusterCells> perCluster = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<o2::aod::EMCALAmbiguousClusterCells> perClusterAmb = o2::aod::emcalclustercell::emcalambiguousclusterId;
  // configurable parameters
  // TODO adapt mDoEventSel switch to also allow selection of other triggers (e.g. EMC7)
  Configurable<bool> mDoEventSel{"doEventSel", 0, "demand kINT7"};
  Configurable<std::string> mVetoBCID{"vetoBCID", "", "BC ID(s) to be excluded, this should be used as an alternative to the event selection"};
  Configurable<std::string> mSelectBCID{"selectBCID", "all", "BC ID(s) to be included, this should be used as an alternative to the event selection"};
  Configurable<double> mVertexCut{"vertexCut", -1, "apply z-vertex cut with value in cm"};
  Configurable<int> mClusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  ConfigurableAxis mClusterTimeBinning{"clustertime-binning", {1500, -600, 900}, ""};
  ConfigurableAxis mNumberClusterBinning{"numclusters-binning", {201, -0.5, 200.5}, ""};

  std::vector<int> mVetoBCIDs;
  std::vector<int> mSelectBCIDs;
  std::vector<float> mCellTime;

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
    const o2Axis energyAxis{makeEnergyBinningAliPhysics(), "E_{clus} (GeV)"};
    const o2Axis amplitudeAxisLarge{1000, 0., 100., "amplitudeLarge", "Amplitude (GeV)"};
    const o2Axis supermoduleAxis{20, -0.5, 19.5, "Supermodule ID"};
    o2Axis timeAxis{mClusterTimeBinning, "t_{cl} (ns)"};
    o2Axis numberClustersAxis{mNumberClusterBinning, "Number of clusters / event"};
    const AxisSpec thAxisCellTimeDiff{3000, -1500, 1500, "#Delta#it{t}_{cell} (ns)"};
    const AxisSpec thAxisCellTimeMean{1500, -600, 900, "#LT#it{t}_{cell}#GT (ns)"};

    // event properties
    mHistManager.add("eventsAll", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsSelected", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventBCAll", "Bunch crossing ID of event (all events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventBCSelected", "Bunch crossing ID of event (selected events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventVertexZAll", "z-vertex of event (all events)", o2HistType::kTH1F, {{200, -20, 20}});
    mHistManager.add("eventVertexZSelected", "z-vertex of event (selected events)", o2HistType::kTH1F, {{200, -20, 20}});
    mHistManager.add("numberOfClustersEvents", "number of clusters per event (selected events)", o2HistType::kTH1F, {numberClustersAxis});
    mHistManager.add("numberOfClustersBC", "number of clusters per bunch crossing (ambiguous BCs)", o2HistType::kTH1F, {numberClustersAxis});
    mHistManager.add("numberOfClustersEventsRejected", "number of clusters per event (rejected events)", o2HistType::kTH1F, {numberClustersAxis});
    mHistManager.add("numberOfClustersSMEvents", "number of clusters per supermodule per event (selected events)", o2HistType::kTH2F, {numberClustersAxis, {20, -0.5, 19.5, "SupermoduleID"}});
    mHistManager.add("numberOfClustersSMBC", "number of clusters per supermodule per bunch crossing (ambiguous BCs)", o2HistType::kTH2F, {numberClustersAxis, {20, -0.5, 19.5, "SupermoduleID"}});

    // cluster properties (matched clusters)
    mHistManager.add("clusterE", "Energy of cluster", o2HistType::kTH1F, {energyAxis});
    mHistManager.add("clusterEMatched", "Energy of cluster (with match)", o2HistType::kTH1F, {energyAxis});
    mHistManager.add("clusterESupermodule", "Energy of the cluster vs. supermoduleID", o2HistType::kTH2F, {energyAxis, supermoduleAxis});
    mHistManager.add("clusterE_SimpleBinning", "Energy of cluster", o2HistType::kTH1F, {{2000, 0, 200}});
    mHistManager.add("clusterEtaPhi", "Eta and phi of cluster", o2HistType::kTH2F, {{100, -1, 1}, {100, 0, 2 * TMath::Pi()}});
    mHistManager.add("clusterM02", "M02 of cluster", o2HistType::kTH1F, {{400, 0, 5}});
    mHistManager.add("clusterM20", "M20 of cluster", o2HistType::kTH1F, {{400, 0, 2.5}});
    mHistManager.add("clusterNLM", "Number of local maxima of cluster", o2HistType::kTH1I, {{10, 0, 10}});
    mHistManager.add("clusterNCells", "Number of cells in cluster", o2HistType::kTH1I, {{50, 0, 50}});
    mHistManager.add("clusterDistanceToBadChannel", "Distance to bad channel", o2HistType::kTH1F, {{100, 0, 100}});
    mHistManager.add("clusterTimeVsE", "Cluster time vs energy", o2HistType::kTH2F, {timeAxis, energyAxis});
    mHistManager.add("clusterAmpFractionLeadingCell", "Fraction of energy in leading cell", o2HistType::kTH1F, {{100, 0, 1}});
    mHistManager.add("clusterCellTimeDiff", "Cell time difference in clusters", o2HistType::kTH1D, {thAxisCellTimeDiff});
    mHistManager.add("clusterCellTimeMean", "Mean cell time per cluster", o2HistType::kTH1D, {thAxisCellTimeMean});

    // add histograms per supermodule
    for (int ism = 0; ism < 20; ++ism) {
      mHistManager.add(Form("clusterTimeVsESM/clusterTimeVsESM%d", ism), Form("Cluster time vs energy in Supermodule %d", ism), o2HistType::kTH2F, {timeAxis, amplitudeAxisLarge});
      mHistManager.add(Form("clusterNcellVsESM/clusterNCellVsESM%d", ism), Form("Cluster number of cells vs energy in Supermodule %d", ism), o2HistType::kTH2F, {{50, 0, 50}, amplitudeAxisLarge});
      mHistManager.add(Form("clusterM02VsESM/clusterM02VsESM%d", ism), Form("Cluster M02 vs energy in Supermodule %d", ism), o2HistType::kTH2F, {{400, 0, 5}, amplitudeAxisLarge});
      mHistManager.add(Form("clusterM20VsESM/clusterM20VsESM%d", ism), Form("Cluster M20 vs energy in Supermodule %d", ism), o2HistType::kTH2F, {{400, 0, 2.5}, amplitudeAxisLarge});
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

  // define cluster filter. It selects only those clusters which are of the type
  // sadly passing of the string at runtime is not possible for technical region so cluster definition is
  // an integer instead
  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == mClusterDefinition);

  /// \brief Process EMCAL clusters that are matched to a collisions
  void processCollisions(collisionEvSelIt const& theCollision, selectedClusters const& clusters, o2::aod::EMCALClusterCells const& emccluscells, o2::aod::Calos const&)
  {
    mHistManager.fill(HIST("eventsAll"), 1);

    // do event selection if mDoEventSel is specified
    // currently the event selection is hard coded to kINT7
    // but other selections are possible that are defined in TriggerAliases.h
    bool isSelected = true;
    if (mDoEventSel) {
      if (theCollision.bc().runNumber() < 300000) {
        if (!theCollision.alias_bit(kINT7)) {
          isSelected = false;
        }
      } else {
        if (!theCollision.alias_bit(kTVXinEMC)) {
          isSelected = false;
        }
      }
    }
    if (!isSelected) {
      LOG(debug) << "Event not selected because it is not kINT7 or does not have EMCAL in readout, skipping";
      mHistManager.fill(HIST("numberOfClustersEventsRejected"), clusters.size());
      return;
    }
    mHistManager.fill(HIST("eventVertexZAll"), theCollision.posZ());
    if (mVertexCut > 0 && TMath::Abs(theCollision.posZ()) > mVertexCut) {
      LOG(debug) << "Event not selected because of z-vertex cut z= " << theCollision.posZ() << " > " << mVertexCut << " cm, skipping";
      return;
    }
    mHistManager.fill(HIST("eventsSelected"), 1);
    mHistManager.fill(HIST("eventVertexZSelected"), theCollision.posZ());
    mHistManager.fill(HIST("numberOfClustersEvents"), clusters.size());

    LOG(debug) << "bunch crossing ID" << theCollision.bcId();
    std::array<int, 20> numberOfClustersSM;
    std::fill(numberOfClustersSM.begin(), numberOfClustersSM.end(), 0);
    // loop over all clusters from accepted collision
    // auto eventClusters = clusters.select(o2::aod::emcalcluster::bcId == theCollision.bc().globalBC());
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
      mHistManager.fill(HIST("clusterE_SimpleBinning"), cluster.energy());
      mHistManager.fill(HIST("clusterEtaPhi"), cluster.eta(), cluster.phi());
      mHistManager.fill(HIST("clusterM02"), cluster.m02());
      mHistManager.fill(HIST("clusterM20"), cluster.m20());
      mHistManager.fill(HIST("clusterTimeVsE"), cluster.time(), cluster.energy());
      mHistManager.fill(HIST("clusterNLM"), cluster.nlm());
      mHistManager.fill(HIST("clusterNCells"), cluster.nCells());
      mHistManager.fill(HIST("clusterDistanceToBadChannel"), cluster.distanceToBadChannel());
      // loop over cells in cluster
      LOG(debug) << "Cluster energy: " << cluster.energy();
      LOG(debug) << "Cluster index: " << cluster.index();
      LOG(debug) << "ncells in cluster: " << cluster.nCells();
      LOG(debug) << "real cluster id" << cluster.id();

      // monitor observables per supermodule
      try {
        auto supermoduleID = mGeometry->SuperModuleNumberFromEtaPhi(cluster.eta(), cluster.phi());
        mHistManager.fill(HIST("clusterESupermodule"), cluster.energy(), supermoduleID);
        fillSupermoduleHistograms(supermoduleID, cluster.energy(), cluster.time(), cluster.nCells(), cluster.m02(), cluster.m20());
        numberOfClustersSM[supermoduleID]++;
      } catch (o2::emcal::InvalidPositionException& e) {
        // Imprecision of the position at the sector boundaries, mostly due to
        // vertex imprecision. Skip these clusters for the now.
      }

      // example of loop over all cells of current cluster
      // cell.calo() allows to access the cell properties as defined in AnalysisDataModel
      // In this exammple, we loop over all cells and find the cell of maximum energy and plot the fraction
      // it carries of the whole cluster
      auto cellsofcluster = emccluscells.sliceBy(perCluster, cluster.globalIndex());
      double maxamp = 0;
      double ampfraction = 0;
      mCellTime.clear();
      mCellTime.reserve(cellsofcluster.size());
      for (const auto& cell : cellsofcluster) {
        // example how to get any information of the cell associated with cluster
        LOG(debug) << "Cell ID:" << cell.calo().amplitude() << " Time " << cell.calo().time();
        if (cell.calo().amplitude() > maxamp) {
          maxamp = cell.calo().amplitude();
        }
        mCellTime.push_back(cell.calo().time());
      } // end of loop over cells
      mHistManager.fill(HIST("clusterCellTimeMean"), std::accumulate(mCellTime.begin(), mCellTime.end(), 0.0f) / mCellTime.size());
      for (std::size_t iCell1 = 0; iCell1 < mCellTime.size() - 1; iCell1++) {
        for (std::size_t iCell2 = iCell1 + 1; iCell2 < mCellTime.size(); iCell2++) {
          mHistManager.fill(HIST("clusterCellTimeDiff"), mCellTime[iCell1] - mCellTime[iCell2]);
        }
      }
      ampfraction = maxamp / cluster.energy();
      mHistManager.fill(HIST("clusterAmpFractionLeadingCell"), ampfraction);
    }
    for (int supermoduleID = 0; supermoduleID < 20; supermoduleID++) {
      mHistManager.fill(HIST("numberOfClustersSMEvents"), numberOfClustersSM[supermoduleID], supermoduleID);
    }
  }
  PROCESS_SWITCH(ClusterMonitor, processCollisions, "Process clusters from collision", false);

  /// \brief Process EMCAL clusters that are not matched to a collision
  /// This is not needed for most users

  void processAmbiguous(bcEvSelIt const& bc, selectedAmbiguousClusters const& clusters)
  {
    // loop over bc , if requested (mVetoBCID >= 0), reject everything from a certain BC
    // this can be used as alternative to event selection (e.g. for pilot beam data)
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
    bool isSelected = true;
    if (mDoEventSel) {
      if (bc.runNumber() < 300000) {
        if (!bc.alias_bit(kINT7)) {
          isSelected = false;
        }
      } else {
        if (!bc.alias_bit(kTVXinEMC)) {
          isSelected = false;
        }
      }
    }
    if (!isSelected) {
      return;
    }
    mHistManager.fill(HIST("eventBCSelected"), eventIR.bc);
    mHistManager.fill(HIST("numberOfClustersBC"), clusters.size());

    std::array<int, 20> numberOfClustersSM;
    std::fill(numberOfClustersSM.begin(), numberOfClustersSM.end(), 0);
    // loop over ambiguous clusters
    for (const auto& cluster : clusters) {
      mHistManager.fill(HIST("clusterE"), cluster.energy());
      mHistManager.fill(HIST("clusterE_SimpleBinning"), cluster.energy());
      mHistManager.fill(HIST("clusterEtaPhi"), cluster.eta(), cluster.phi());
      mHistManager.fill(HIST("clusterM02"), cluster.m02());
      mHistManager.fill(HIST("clusterM20"), cluster.m20());
      mHistManager.fill(HIST("clusterTimeVsE"), cluster.time(), cluster.energy());
      mHistManager.fill(HIST("clusterNLM"), cluster.nlm());
      mHistManager.fill(HIST("clusterNCells"), cluster.nCells());
      mHistManager.fill(HIST("clusterDistanceToBadChannel"), cluster.distanceToBadChannel());

      try {
        auto supermoduleID = mGeometry->SuperModuleNumberFromEtaPhi(cluster.eta(), cluster.phi());
        mHistManager.fill(HIST("clusterESupermodule"), cluster.energy(), supermoduleID);
        fillSupermoduleHistograms(supermoduleID, cluster.energy(), cluster.time(), cluster.nCells(), cluster.m02(), cluster.m20());
        numberOfClustersSM[supermoduleID]++;
      } catch (o2::emcal::InvalidPositionException& e) {
        // Imprecision of the position at the sector boundaries, mostly due to
        // vertex imprecision. Skip these clusters for the now.
      }
    }
    for (int supermoduleID = 0; supermoduleID < 20; supermoduleID++) {
      mHistManager.fill(HIST("numberOfClustersSMBC"), numberOfClustersSM[supermoduleID], supermoduleID);
    }
  }
  PROCESS_SWITCH(ClusterMonitor, processAmbiguous, "Process Ambiguous clusters", false);

  void fillSupermoduleHistograms(int supermoduleID, double energy, double time, int ncell, double m02, double m20)
  {
    // workaround to have the histogram names per supermodule
    // handled as constexpr
    switch (supermoduleID) {
      case 0:
        supermoduleHistHelper<0>(energy, time, ncell, m02, m20);
        break;
      case 1:
        supermoduleHistHelper<1>(energy, time, ncell, m02, m20);
        break;
      case 2:
        supermoduleHistHelper<2>(energy, time, ncell, m02, m20);
        break;
      case 3:
        supermoduleHistHelper<3>(energy, time, ncell, m02, m20);
        break;
      case 4:
        supermoduleHistHelper<4>(energy, time, ncell, m02, m20);
        break;
      case 5:
        supermoduleHistHelper<5>(energy, time, ncell, m02, m20);
        break;
      case 6:
        supermoduleHistHelper<6>(energy, time, ncell, m02, m20);
        break;
      case 7:
        supermoduleHistHelper<7>(energy, time, ncell, m02, m20);
        break;
      case 8:
        supermoduleHistHelper<8>(energy, time, ncell, m02, m20);
        break;
      case 9:
        supermoduleHistHelper<9>(energy, time, ncell, m02, m20);
        break;
      case 10:
        supermoduleHistHelper<10>(energy, time, ncell, m02, m20);
        break;
      case 11:
        supermoduleHistHelper<11>(energy, time, ncell, m02, m20);
        break;
      case 12:
        supermoduleHistHelper<12>(energy, time, ncell, m02, m20);
        break;
      case 13:
        supermoduleHistHelper<13>(energy, time, ncell, m02, m20);
        break;
      case 14:
        supermoduleHistHelper<14>(energy, time, ncell, m02, m20);
        break;
      case 15:
        supermoduleHistHelper<15>(energy, time, ncell, m02, m20);
        break;
      case 16:
        supermoduleHistHelper<16>(energy, time, ncell, m02, m20);
        break;
      case 17:
        supermoduleHistHelper<17>(energy, time, ncell, m02, m20);
        break;
      case 18:
        supermoduleHistHelper<18>(energy, time, ncell, m02, m20);
        break;
      case 19:
        supermoduleHistHelper<19>(energy, time, ncell, m02, m20);
        break;
      default:
        break;
    }
  }

  template <int supermoduleID>
  void supermoduleHistHelper(double energy, double time, int ncells, double m02, double m20)
  {
    static constexpr std::string_view clusterEnergyTimeHist[20] = {"clusterTimeVsESM/clusterTimeVsESM0", "clusterTimeVsESM/clusterTimeVsESM1", "clusterTimeVsESM/clusterTimeVsESM2", "clusterTimeVsESM/clusterTimeVsESM3", "clusterTimeVsESM/clusterTimeVsESM4", "clusterTimeVsESM/clusterTimeVsESM5", "clusterTimeVsESM/clusterTimeVsESM6", "clusterTimeVsESM/clusterTimeVsESM7", "clusterTimeVsESM/clusterTimeVsESM8", "clusterTimeVsESM/clusterTimeVsESM9", "clusterTimeVsESM/clusterTimeVsESM10", "clusterTimeVsESM/clusterTimeVsESM11", "clusterTimeVsESM/clusterTimeVsESM12", "clusterTimeVsESM/clusterTimeVsESM13", "clusterTimeVsESM/clusterTimeVsESM14", "clusterTimeVsESM/clusterTimeVsESM15", "clusterTimeVsESM/clusterTimeVsESM16", "clusterTimeVsESM/clusterTimeVsESM17", "clusterTimeVsESM/clusterTimeVsESM18", "clusterTimeVsESM/clusterTimeVsESM19"};
    static constexpr std::string_view clusterEnergyNcellHistSM[20] = {"clusterNcellVsESM/clusterNCellVsESM0", "clusterNcellVsESM/clusterNCellVsESM1", "clusterNcellVsESM/clusterNCellVsESM2", "clusterNcellVsESM/clusterNCellVsESM3", "clusterNcellVsESM/clusterNCellVsESM4", "clusterNcellVsESM/clusterNCellVsESM5", "clusterNcellVsESM/clusterNCellVsESM6", "clusterNcellVsESM/clusterNCellVsESM7", "clusterNcellVsESM/clusterNCellVsESM8", "clusterNcellVsESM/clusterNCellVsESM9", "clusterNcellVsESM/clusterNCellVsESM10", "clusterNcellVsESM/clusterNCellVsESM11", "clusterNcellVsESM/clusterNCellVsESM12", "clusterNcellVsESM/clusterNCellVsESM13", "clusterNcellVsESM/clusterNCellVsESM14", "clusterNcellVsESM/clusterNCellVsESM15", "clusterNcellVsESM/clusterNCellVsESM16", "clusterNcellVsESM/clusterNCellVsESM17", "clusterNcellVsESM/clusterNCellVsESM18", "clusterNcellVsESM/clusterNCellVsESM19"};
    static constexpr std::string_view clusterEnergyM02HistSM[20] = {"clusterM02VsESM/clusterM02VsESM0", "clusterM02VsESM/clusterM02VsESM1", "clusterM02VsESM/clusterM02VsESM2", "clusterM02VsESM/clusterM02VsESM3", "clusterM02VsESM/clusterM02VsESM4", "clusterM02VsESM/clusterM02VsESM5", "clusterM02VsESM/clusterM02VsESM6", "clusterM02VsESM/clusterM02VsESM7", "clusterM02VsESM/clusterM02VsESM8", "clusterM02VsESM/clusterM02VsESM9", "clusterM02VsESM/clusterM02VsESM10", "clusterM02VsESM/clusterM02VsESM11", "clusterM02VsESM/clusterM02VsESM12", "clusterM02VsESM/clusterM02VsESM13", "clusterM02VsESM/clusterM02VsESM14", "clusterM02VsESM/clusterM02VsESM15", "clusterM02VsESM/clusterM02VsESM16", "clusterM02VsESM/clusterM02VsESM17", "clusterM02VsESM/clusterM02VsESM18", "clusterM02VsESM/clusterM02VsESM19"};
    static constexpr std::string_view clusterEnergyM20HistSM[20] = {"clusterM20VsESM/clusterM20VsESM0", "clusterM20VsESM/clusterM20VsESM1", "clusterM20VsESM/clusterM20VsESM2", "clusterM20VsESM/clusterM20VsESM3", "clusterM20VsESM/clusterM20VsESM4", "clusterM20VsESM/clusterM20VsESM5", "clusterM20VsESM/clusterM20VsESM6", "clusterM20VsESM/clusterM20VsESM7", "clusterM20VsESM/clusterM20VsESM8", "clusterM20VsESM/clusterM20VsESM9", "clusterM20VsESM/clusterM20VsESM10", "clusterM20VsESM/clusterM20VsESM11", "clusterM20VsESM/clusterM20VsESM12", "clusterM20VsESM/clusterM20VsESM13", "clusterM20VsESM/clusterM20VsESM14", "clusterM20VsESM/clusterM20VsESM15", "clusterM20VsESM/clusterM20VsESM16", "clusterM20VsESM/clusterM20VsESM17", "clusterM20VsESM/clusterM20VsESM18", "clusterM20VsESM/clusterM20VsESM19"};
    mHistManager.fill(HIST(clusterEnergyTimeHist[supermoduleID]), time, energy);
    mHistManager.fill(HIST(clusterEnergyNcellHistSM[supermoduleID]), ncells, energy);
    mHistManager.fill(HIST(clusterEnergyM02HistSM[supermoduleID]), m02, energy);
    mHistManager.fill(HIST(clusterEnergyM20HistSM[supermoduleID]), m20, energy);
  }

  /// \brief Create binning for cluster energy axis (variable bin size)
  /// \return vector with bin limits
  std::vector<double>
    makeEnergyBinning() const
  {
    auto fillBinLimits = [](std::vector<double>& binlimits, double max, double binwidth) {
      auto current = *binlimits.rbegin();
      while (current < max) {
        current += binwidth;
        binlimits.emplace_back(current);
      }
    };
    std::vector<double> result = {0.};
    fillBinLimits(result, 2., 0.1);
    fillBinLimits(result, 5., 0.2);
    fillBinLimits(result, 10., 0.5);
    fillBinLimits(result, 20., 1.);
    fillBinLimits(result, 50., 2.);
    fillBinLimits(result, 100., 5.);
    fillBinLimits(result, 200., 10.);
    return result;
  }

  /// \brief Create binning for cluster energy axis (variable bin size)
  /// direct port from binning often used in AliPhysics for debugging
  /// \return vector with bin limits
  std::vector<double> makeEnergyBinningAliPhysics() const
  {

    std::vector<double> result;
    Int_t nBinsClusterE = 235;
    for (Int_t i = 0; i < nBinsClusterE + 1; i++) {
      if (i < 1)
        result.emplace_back(0.3 * i);
      else if (i < 55)
        result.emplace_back(0.3 + 0.05 * (i - 1));
      else if (i < 105)
        result.emplace_back(3. + 0.1 * (i - 55));
      else if (i < 140)
        result.emplace_back(8. + 0.2 * (i - 105));
      else if (i < 170)
        result.emplace_back(15. + 0.5 * (i - 140));
      else if (i < 190)
        result.emplace_back(30. + 1.0 * (i - 170));
      else if (i < 215)
        result.emplace_back(50. + 2.0 * (i - 190));
      else if (i < 235)
        result.emplace_back(100. + 5.0 * (i - 215));
      else if (i < 245)
        result.emplace_back(200. + 10.0 * (i - 235));
    }
    return result;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<ClusterMonitor>(cfgc, TaskName{"EMCClusterMonitorTask"}, SetDefaultProcesses{{{"processCollisions", true}, {"processAmbiguous", false}}}),
    adaptAnalysisTask<ClusterMonitor>(cfgc, TaskName{"EMCClusterMonitorTaskAmbiguous"}, SetDefaultProcesses{{{"processCollisions", false}, {"processAmbiguous", true}}})};
  return workflow;
}
