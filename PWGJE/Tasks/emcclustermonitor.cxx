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

using collisionEvSelIt = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator;

struct ClusterMonitor {
  o2::framework::HistogramRegistry mHistManager{"ClusterMonitorHistograms"};
  o2::emcal::Geometry* mGeometry = nullptr;

  // configurable parameters
  // TODO adapt mDoEventSel switch to also allow selection of other triggers (e.g. EMC7)
  o2::framework::Configurable<bool> mDoEventSel{"doEventSel", 0, "demand kINT7"};
  o2::framework::Configurable<std::string> mVetoBCID{"vetoBCID", "", "BC ID(s) to be excluded, this should be used as an alternative to the event selection"};
  o2::framework::Configurable<std::string> mSelectBCID{"selectBCID", "all", "BC ID(s) to be included, this should be used as an alternative to the event selection"};
  std::vector<int> mVetoBCIDs;
  std::vector<int> mSelectBCIDs;

  // TODO add multiple configurables that allow to cut on cluster properties

  /// \brief Create output histograms and initialize geometry
  void init(o2::framework::InitContext const&)
  {
    // create histograms
    using o2HistType = o2::framework::HistType;
    using o2Axis = o2::framework::AxisSpec;

    // load geometry just in case we need it
    mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(300000);

    // create histograms for cluster QA
    Double_t timeMin = -600;
    Double_t timeMax = 900;

    // create common axes
    LOG(info) << "Creating histograms";
    const o2Axis bcAxis{3501, -0.5, 3500.5};
    const o2Axis energyAxis{makeEnergyBinning(), "E_{clus} (GeV)"};
    const o2Axis timeAxis{800, timeMin, timeMax};

    // event properties
    mHistManager.add("eventsAll", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsSelected", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventBCAll", "Bunch crossing ID of event (all events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventBCSelected", "Bunch crossing ID of event (selected events)", o2HistType::kTH1F, {bcAxis});

    // cluster properties
    mHistManager.add("clusterE", "Energy of cluster", o2HistType::kTH1F, {{600, 0, 100}});
    mHistManager.add("clusterEtaPhi", "Eta and phi of cluster", o2HistType::kTH2F, {{100, -1, 1}, {100, 0, 2 * TMath::Pi()}});
    mHistManager.add("clusterM02", "M02 of cluster", o2HistType::kTH1F, {{100, 0, 3}});
    mHistManager.add("clusterM20", "M20 of cluster", o2HistType::kTH1F, {{100, 0, 3}});
    mHistManager.add("clusterNLM", "Number of local maxima of cluster", o2HistType::kTH1I, {{20, 0, 20}});
    mHistManager.add("clusterNCells", "Number of cells in cluster", o2HistType::kTH1I, {{100, 0, 100}});
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
  /// \brief Process EMCAL clusters
  void process(collisionEvSelIt const& theCollision, o2::aod::EMCALClusters const& clusters, o2::aod::BCs const& bcs)
  {

    mHistManager.fill(HIST("eventsAll"), 1);

    // do event selection if mDoEventSel is specified
    // currently the event selection is hard coded to kINT7
    // but other selections are possible that are defined in TriggerAliases.h
    if (mDoEventSel && (!theCollision.alias()[kINT7])) {
      LOG(info) << "Event not selected becaus it is not kINT7, skipping";
      return;
    }
    mHistManager.fill(HIST("eventsSelected"), 1);

    // loop over bc , if requested (mVetoBCID >= 0), reject everything from a certain BC
    // this can be used as alternative to event selection (e.g. for pilot beam data)
    for (const auto& bc : bcs) {
      o2::InteractionRecord eventIR;
      eventIR.setFromLong(bc.globalBC());
      mHistManager.fill(HIST("eventBCAll"), eventIR.bc);
      if (std::find(mVetoBCIDs.begin(), mVetoBCIDs.end(), eventIR.bc) != mVetoBCIDs.end()) {
        LOG(info) << "Event rejected because of veto BCID " << eventIR.bc;
        continue;
      }
      if (mSelectBCIDs.size() && (std::find(mSelectBCIDs.begin(), mSelectBCIDs.end(), eventIR.bc) == mSelectBCIDs.end())) {
        continue;
      }
      mHistManager.fill(HIST("eventBCSelected"), eventIR.bc);
    }

    // loop over all clusters from accepted collision
    for (const auto& cluster : clusters) {
      // fill histograms of cluster properties
      // in this implementation the cluster properties are directly
      // loaded from the flat table, in the future one should
      // consider using the AnalysisCluster object to work with
      // after loading.
      o2::InteractionRecord eventIR;
      eventIR.setFromLong(cluster.bc().globalBC());
      if (std::find(mVetoBCIDs.begin(), mVetoBCIDs.end(), eventIR.bc) != mVetoBCIDs.end()) {
        continue;
      }
      if (mSelectBCIDs.size() && (std::find(mSelectBCIDs.begin(), mSelectBCIDs.end(), eventIR.bc) == mSelectBCIDs.end())) {
        continue;
      }
      LOG(debug) << "Cluster energy: " << cluster.energy();
      LOG(debug) << "Cluster time: " << cluster.time();
      LOG(debug) << "Cluster M02: " << cluster.m02();
      mHistManager.fill(HIST("clusterE"), cluster.energy());
      mHistManager.fill(HIST("clusterEtaPhi"), cluster.eta(), cluster.phi());
      mHistManager.fill(HIST("clusterM02"), cluster.m02());
      mHistManager.fill(HIST("clusterM20"), cluster.m20());
      mHistManager.fill(HIST("clusterTimeVsE"), cluster.time(), cluster.energy());
      mHistManager.fill(HIST("clusterNLM"), cluster.nlm());
      mHistManager.fill(HIST("clusterNCells"), cluster.nCells());
      mHistManager.fill(HIST("clusterDistanceToBadChannel"), cluster.distanceToBadChannel());
    }
  }

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
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<ClusterMonitor>(cfgc)};
}