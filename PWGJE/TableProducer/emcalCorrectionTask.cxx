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

// EMCAL Correction Task
//
// Author: Raymond Ehlers & Florian Jonas

#include <algorithm>
#include <cmath>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "DetectorsBase/GeometryManager.h"

#include "PWGJE/DataModel/EMCALClusters.h"

#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"
#include "EMCALBase/Geometry.h"
#include "EMCALBase/ClusterFactory.h"
#include "EMCALReconstruction/Clusterizer.h"
#include "TVector2.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct EmcalCorrectionTask {
  Produces<o2::aod::EMCALClusters> clusters;
  Produces<o2::aod::EMCALAmbiguousClusters> clustersAmbiguous;

  // Options for the clusterization
  // 1 corresponds to EMCAL cells based on the Run2 definition.
  Configurable<int> selectedCellType{"selectedCellType", 1, "EMCAL Cell type"};
  Configurable<std::string> clusterDefinitions{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default. Multiple definitions can be specified separated by comma"};
  // CDB service (for geometry)
  Service<o2::ccdb::BasicCCDBManager> mCcdbManager;

  // Clusterizer and related
  // Apparently streaming these objects really doesn't work, and causes problems for setting up the workflow.
  // So we use unique_ptr and define them below.
  std::vector<std::unique_ptr<o2::emcal::Clusterizer<o2::emcal::Cell>>> mClusterizers;
  std::vector<std::unique_ptr<o2::emcal::ClusterFactory<o2::emcal::Cell>>> mClusterFactories;
  // Cells and clusters
  std::vector<o2::emcal::Cell> mEmcalCells;
  std::vector<o2::emcal::AnalysisCluster> mAnalysisClusters;

  std::vector<o2::aod::EMCALClusterDefinition> mClusterDefinitions;
  // QA
  // NOTE: This is not comprehensive.
  OutputObj<TH1F> hCellE{"hCellE"};
  OutputObj<TH1I> hCellTowerID{"hCellTowerID"};
  OutputObj<TH2F> hCellEtaPhi{"hCellEtaPhi"};
  OutputObj<TH2I> hCellRowCol{"hCellRowCol"};
  OutputObj<TH1F> hClusterE{"hClusterE"};
  OutputObj<TH2F> hClusterEtaPhi{"hClusterEtaPhi"};

  void init(InitContext const&)
  {
    LOG(debug) << "Start init!";
    // NOTE: The geometry manager isn't necessary just to load the EMCAL geometry.
    //       However, it _is_ necessary for loading the misalignment matrices as of September 2020
    //       Eventually, those matrices will be moved to the CCDB, but it's not yet ready.

    // The geomatrices from the CCDB are needed in order to calculate the cluster kinematics
    const char* ccdburl = "http://alice-ccdb.cern.ch"; // Possibly the parameter might become configurable, in order to be able to load new EMCAL calibration objects
    mCcdbManager->setURL(ccdburl);
    mCcdbManager->setCaching(true);
    mCcdbManager->setLocalObjectValidityChecking();
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      mCcdbManager->get<TGeoManager>("GLO/Config/Geometry");
    }
    LOG(debug) << "After load geometry!";
    o2::emcal::Geometry* geometry = o2::emcal::Geometry::GetInstanceFromRunNumber(223409);
    if (!geometry) {
      LOG(error) << "Failure accessing geometry";
    }

    // read all the cluster definitions specified in the options
    if (clusterDefinitions->length()) {
      std::stringstream parser(clusterDefinitions.value);
      std::string token;
      o2::aod::EMCALClusterDefinition clusDef;
      while (std::getline(parser, token, ',')) {
        clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(token);
        mClusterDefinitions.push_back(clusDef);
      }
    }
    for (auto& clusterDefinition : mClusterDefinitions) {
      mClusterizers.emplace_back(std::make_unique<o2::emcal::Clusterizer<o2::emcal::Cell>>(1, clusterDefinition.timeMin, clusterDefinition.timeMax, clusterDefinition.gradientCut, true, clusterDefinition.seedEnergy, clusterDefinition.minCellEnergy));
      mClusterFactories.emplace_back(std::make_unique<o2::emcal::ClusterFactory<o2::emcal::Cell>>());
      LOG(info) << "Cluster definition initialized: " << clusterDefinition.toString();
      LOG(info) << "timeMin: " << clusterDefinition.timeMin;
      LOG(info) << "timeMax: " << clusterDefinition.timeMax;
      LOG(info) << "gradientCut: " << clusterDefinition.gradientCut;
      LOG(info) << "seedEnergy: " << clusterDefinition.seedEnergy;
      LOG(info) << "minCellEnergy: " << clusterDefinition.minCellEnergy;
    }
    for (auto& clusterizer : mClusterizers) {
      clusterizer->setGeometry(geometry);
    }

    if (mClusterizers.size() == 0) {
      LOG(error) << "No cluster definitions specified!";
    }

    LOG(debug) << "Completed init!";

    // Setup QA hists.
    hCellE.setObject(new TH1F("hCellE", "hCellE", 200, 0.0, 100));
    hCellTowerID.setObject(new TH1I("hCellTowerID", "hCellTowerID", 20000, 0, 20000));
    hCellEtaPhi.setObject(new TH2F("hCellEtaPhi", "hCellEtaPhi", 160, -0.8, 0.8, 72, 0, 2 * 3.14159));
    // NOTE: Reversed column and row because it's more natural for presentatin.
    hCellRowCol.setObject(new TH2I("hCellRowCol", "hCellRowCol;Column;Row", 97, 0, 97, 600, 0, 600));
    hClusterE.setObject(new TH1F("hClusterE", "hClusterE", 200, 0.0, 100));
    hClusterEtaPhi.setObject(new TH2F("hClusterEtaPhi", "hClusterEtaPhi", 160, -0.8, 0.8, 72, 0, 2 * 3.14159));
  }

  // void process(aod::Collision const& collision, soa::Filtered<aod::Tracks> const& fullTracks, aod::Calos const& cells)
  // void process(aod::Collision const& collision, aod::Tracks const& tracks, aod::Calos const& cells)
  // void process(aod::BCs const& bcs, aod::Collision const& collision, aod::Calos const& cells)

  //  Appears to need the BC to be accessed to be available in the collision table...
  void process(aod::BC const& bc, aod::Collisions const& collisions, aod::Calos const& cells)
  {
    LOG(debug) << "Starting process.";
    // Convert aod::Calo to o2::emcal::Cell which can be used with the clusterizer.
    // In particular, we need to filter only EMCAL cells.
    mEmcalCells.clear();
    for (auto& cell : cells) {
      if (cell.caloType() != selectedCellType) {
        // LOG(debug) << "Rejected";
        continue;
      }
      // LOG(debug) << "Cell E: " << cell.getEnergy();
      // LOG(debug) << "Cell E: " << cell;
      mEmcalCells.emplace_back(o2::emcal::Cell(
        cell.cellNumber(),
        cell.amplitude(),
        cell.time(),
        o2::emcal::intToChannelType(cell.cellType())));
    }

    // Cell QA
    // For convenience, use the clusterizer stored geometry to get the eta-phi
    for (auto& cell : mEmcalCells) {
      hCellE->Fill(cell.getEnergy());
      hCellTowerID->Fill(cell.getTower());
      auto res = mClusterizers.at(0)->getGeometry()->EtaPhiFromIndex(cell.getTower());
      hCellEtaPhi->Fill(std::get<0>(res), TVector2::Phi_0_2pi(std::get<1>(res)));
      res = mClusterizers.at(0)->getGeometry()->GlobalRowColFromIndex(cell.getTower());
      // NOTE: Reversed column and row because it's more natural for presentatin.
      hCellRowCol->Fill(std::get<1>(res), std::get<0>(res));
    }

    // TODO: Helpful for now, but should be removed.
    LOG(debug) << "Converted EMCAL cells";
    for (auto& cell : mEmcalCells) {
      LOG(debug) << cell.getTower() << ": E: " << cell.getEnergy() << ", time: " << cell.getTimeStamp() << ", type: " << cell.getType();
    }

    LOG(debug) << "Converted cells. Contains: " << mEmcalCells.size() << ". Originally " << cells.size() << ". About to run clusterizer.";
    // this is a test
    // Run the clusterizers
    LOG(debug) << "Running clusterizers";
    Int_t i = 0;
    for (auto& clusterizer : mClusterizers) {
      clusterizer->findClusters(mEmcalCells);

      auto emcalClusters = clusterizer->getFoundClusters();
      auto emcalClustersInputIndices = clusterizer->getFoundClustersInputIndices();
      LOG(debug) << "Retrieved results. About to setup cluster factory.";

      // Convert to analysis clusters.
      // First, the cluster factory requires cluster and cell information in order to build the clusters.
      mAnalysisClusters.clear();
      mClusterFactories.at(i)->reset();
      mClusterFactories.at(i)->setClustersContainer(*emcalClusters);
      mClusterFactories.at(i)->setCellsContainer(mEmcalCells);
      mClusterFactories.at(i)->setCellsIndicesContainer(*emcalClustersInputIndices);

      LOG(debug) << "Cluster factory set up.";
      // Convert to analysis clusters.
      for (int icl = 0; icl < mClusterFactories.at(i)->getNumberOfClusters(); icl++) {
        auto analysisCluster = mClusterFactories.at(i)->buildCluster(icl);
        mAnalysisClusters.emplace_back(analysisCluster);
      }
      LOG(debug) << "Converted to analysis clusters.";

      float vx = 0, vy = 0, vz = 0;
      bool hasCollision = false;
      if (collisions.size() > 1) {
        LOG(error) << "More than one collision in the bc. This is not supported.";
      } else {
        // dummy loop to get the first collision
        for (const auto& col : collisions) {
          vx = col.posX();
          vy = col.posY();
          vz = col.posZ();
          hasCollision = true;

          // we found a collision, put the clusters into the none ambiguous table
          clusters.reserve(mAnalysisClusters.size());
          for (const auto& cluster : mAnalysisClusters) {

            // Determine the cluster eta, phi, correcting for the vertex position.
            auto pos = cluster.getGlobalPosition();
            pos = pos - math_utils::Point3D<float>{vx, vy, vz};
            // Normalize the vector and rescale by energy.
            pos *= (cluster.E() / std::sqrt(pos.Mag2()));

            // save to table
            LOG(debug) << "Writing cluster definition " << static_cast<int>(mClusterDefinitions.at(i)) << " to table.";
            clusters(col, cluster.getID(), cluster.E(), cluster.getCoreEnergy(), pos.Eta(), TVector2::Phi_0_2pi(pos.Phi()),
                     cluster.getM02(), cluster.getM20(), cluster.getNCells(), cluster.getClusterTime(),
                     cluster.getIsExotic(), cluster.getDistanceToBadChannel(), cluster.getNExMax(), static_cast<int>(mClusterDefinitions.at(i)));

            // fill histograms
            hClusterE->Fill(cluster.E());
            hClusterEtaPhi->Fill(pos.Eta(), TVector2::Phi_0_2pi(pos.Phi()));
          } // end of cluster loop
        }   // end of collision loop
      }
      if (!hasCollision) {
        LOG(warning) << "No vertex found for event. Assuming (0,0,0).";
      }

      // Store the clusters in the table where a mathcing collision could
      // be identified.
      if (!hasCollision) { // ambiguous
        clustersAmbiguous.reserve(mAnalysisClusters.size());
        for (const auto& cluster : mAnalysisClusters) {
          auto pos = cluster.getGlobalPosition();
          pos = pos - math_utils::Point3D<float>{vx, vy, vz};
          // Normalize the vector and rescale by energy.
          pos *= (cluster.E() / std::sqrt(pos.Mag2()));

          // We have our necessary properties. Now we store outputs

          // LOG(debug) << "Cluster E: " << cluster.E();
          clustersAmbiguous(bc, cluster.getID(), cluster.E(), cluster.getCoreEnergy(), pos.Eta(), TVector2::Phi_0_2pi(pos.Phi()),
                            cluster.getM02(), cluster.getM20(), cluster.getNCells(), cluster.getClusterTime(),
                            cluster.getIsExotic(), cluster.getDistanceToBadChannel(), cluster.getNExMax(), static_cast<int>(mClusterDefinitions.at(i)));
        }
      }
      LOG(debug) << "Cluster loop done for clusterizer " << i;
      i++;
    } // end of clusterizer loop
    LOG(debug) << "Done with process.";
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<EmcalCorrectionTask>(cfgc, TaskName{"emcal-correction-task"})};
}
