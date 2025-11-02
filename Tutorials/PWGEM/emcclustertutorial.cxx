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

#include "PWGJE/DataModel/EMCALClusters.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// \struct ClusterTutorial
/// \brief Skeleton task for EMCal cluster tutorial
/// \author Florian Jonas <florian.jonas@cern.ch>, Oak Ridge National Laoratory; Joshua König <joshua.konig@cern.ch>, Goethe-University Frankfurt; Marvin Hemmer <marvin.hemmer@cern.ch>, Goethe-University Frankfurt
/// \since 09.04.2023
using namespace o2::framework;
using namespace o2::framework::expressions;
using collisionEvSelIt = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator;

// we need to filter the emcal clusters, since we only want clusters from one specific clusterizer algorithm
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
struct ClusterTutorial {
  HistogramRegistry mHistManager{"ClusterMonitorHistograms"};

  Preslice<o2::aod::EMCALClusterCells> perCluster = o2::aod::emcalclustercell::emcalclusterId;
  // configurable parameters
  Configurable<bool> mDoEventSel{"doEventSel", 0, "demand kINT7"};
  Configurable<double> mVertexCut{"vertexCut", -1, "apply z-vertex cut with abs value in cm"};
  Configurable<int> mClusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};

  // define cluster filter. It selects only those clusters which are of the type.
  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == mClusterDefinition);

  /// \brief Create output histograms and initialize geometry
  void init(InitContext const&)
  {
    mHistManager.add("clusterE", "Energy of cluster;'it{E} (GeV)", HistType::kTH1F, {{200, 0., 40.}});
  }

  /// \brief Process EMCAL clusters that are matched to a collisions
  void processCalo(collisionEvSelIt const& theCollision, selectedClusters const& clusters, o2::aod::EMCALClusterCells const& /*emccluscells*/)
  {

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
      return;
    }
    if (mVertexCut > 0 && TMath::Abs(theCollision.posZ()) > mVertexCut) {
      return;
    }
    // loop over all clusters from accepted collision
    for (const auto& cluster : clusters) {
      mHistManager.fill(HIST("clusterE"), cluster.energy());
    }
  }
  PROCESS_SWITCH(ClusterTutorial, processCalo, "process only reconstructed info", true);

  void processDummy(collisionEvSelIt const&)
  {
    // do nothing
  }
  PROCESS_SWITCH(ClusterTutorial, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ClusterTutorial>(cfgc)};
}
