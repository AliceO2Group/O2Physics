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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"

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
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
using selectedAmbiguousClusters = o2::soa::Filtered<o2::aod::EMCALAmbiguousClusters>;
struct ClusterMonitor {
  HistogramRegistry mHistManager{"ClusterMonitorHistograms"};
  o2::emcal::Geometry* mGeometry = nullptr;

  // configurable parameters
  // TODO adapt mDoEventSel switch to also allow selection of other triggers (e.g. EMC7)
  Configurable<bool> mDoEventSel{"doEventSel", 0, "demand kINT7"};
  Configurable<std::string> mVetoBCID{"vetoBCID", "", "BC ID(s) to be excluded, this should be used as an alternative to the event selection"};
  Configurable<std::string> mSelectBCID{"selectBCID", "all", "BC ID(s) to be included, this should be used as an alternative to the event selection"};
  Configurable<double> mVertexCut{"vertexCut", -1, "apply z-vertex cut with value in cm"};
  Configurable<std::string> mClusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  std::vector<int> mVetoBCIDs;
  std::vector<int> mSelectBCIDs;

  // define cluster filter. It selects only those clusters which are of the type
  // specified in the string mClusterDefinition,e.g. kV3Default, which is V3 clusterizer with default
  // clusterization parameters
  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);

  /// \brief Create output histograms and initialize geometry
  void init(InitContext const&)
  {
    // create histograms
    using o2HistType = HistType;
    using o2Axis = AxisSpec;

    // load geometry just in case we need it
    mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(300000);

    // create histograms for cluster QA
    Double_t timeMin = -600;
    Double_t timeMax = 900;

    // create common axes
    LOG(info) << "Creating histograms";
    const o2Axis bcAxis{3501, -0.5, 3500.5};
    const o2Axis energyAxis{makeEnergyBinningAliPhysics(), "E_{clus} (GeV)"};
    const o2Axis timeAxis{800, timeMin, timeMax};

    // event properties
    mHistManager.add("eventsAll", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsSelected", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventBCAll", "Bunch crossing ID of event (all events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventBCSelected", "Bunch crossing ID of event (selected events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventVertexZAll", "z-vertex of event (all events)", o2HistType::kTH1F, {{200, -20, 20}});
    mHistManager.add("eventVertexZSelected", "z-vertex of event (selected events)", o2HistType::kTH1F, {{200, -20, 20}});

    // cluster properties (matched clusters)
    mHistManager.add("clusterE", "Energy of cluster", o2HistType::kTH1F, {energyAxis});
    mHistManager.add("clusterE_SimpleBinning", "Energy of cluster", o2HistType::kTH1F, {{400, 0, 100}});
    mHistManager.add("clusterEtaPhi", "Eta and phi of cluster", o2HistType::kTH2F, {{100, -1, 1}, {100, 0, 2 * TMath::Pi()}});
    mHistManager.add("clusterM02", "M02 of cluster", o2HistType::kTH1F, {{400, 0, 5}});
    mHistManager.add("clusterM20", "M20 of cluster", o2HistType::kTH1F, {{400, 0, 2.5}});
    mHistManager.add("clusterNLM", "Number of local maxima of cluster", o2HistType::kTH1I, {{10, 0, 10}});
    mHistManager.add("clusterNCells", "Number of cells in cluster", o2HistType::kTH1I, {{50, 0, 50}});
    mHistManager.add("clusterDistanceToBadChannel", "Distance to bad channel", o2HistType::kTH1F, {{100, 0, 100}});
    mHistManager.add("clusterTimeVsE", "Cluster time vs energy", o2HistType::kTH2F, {timeAxis, energyAxis});

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
  void processCollisions(collisionEvSelIt const& theCollision, selectedClusters const& clusters, o2::aod::BCs const& bcs)
  {
    mHistManager.fill(HIST("eventsAll"), 1);

    // do event selection if mDoEventSel is specified
    // currently the event selection is hard coded to kINT7
    // but other selections are possible that are defined in TriggerAliases.h
    if (mDoEventSel && (!theCollision.alias()[kINT7])) {
      LOG(debug) << "Event not selected becaus it is not kINT7, skipping";
      return;
    }
    mHistManager.fill(HIST("eventVertexZAll"), theCollision.posZ());
    if (mVertexCut > 0 && TMath::Abs(theCollision.posZ()) > mVertexCut) {
      LOG(debug) << "Event not selected because of z-vertex cut z= " << theCollision.posZ() << " > " << mVertexCut << " cm, skipping";
      return;
    }
    mHistManager.fill(HIST("eventVertexZAll"), theCollision.posZ());
    if (mVertexCut > 0 && TMath::Abs(theCollision.posZ()) > mVertexCut) {
      LOG(debug) << "Event not selected because of z-vertex cut z= " << theCollision.posZ() << " > " << mVertexCut << " cm, skipping";
      return;
    }
    mHistManager.fill(HIST("eventsSelected"), 1);
    mHistManager.fill(HIST("eventVertexZSelected"), theCollision.posZ());

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
    }
  }
  PROCESS_SWITCH(ClusterMonitor, processCollisions, "Process clusters from collision", false);

  /// \brief Process EMCAL clusters that are not matched to a collision
  /// This is not needed for most users

  void processAmbiguous(o2::aod::BC const bc, selectedAmbiguousClusters const& clusters)
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
    mHistManager.fill(HIST("eventBCSelected"), eventIR.bc);
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
    }
  }
  PROCESS_SWITCH(ClusterMonitor, processAmbiguous, "Process Ambiguous clusters", false);

  /// \brief Create binning for cluster energy axis (variable bin size)
  /// \return vector with bin limits
  std::vector<double> makeEnergyBinning() const
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