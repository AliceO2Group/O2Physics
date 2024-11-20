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
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
/// \author Florian Jonas <florian.jonas@cern.ch>

#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <string>
#include <tuple>
#include <vector>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "DetectorsBase/GeometryManager.h"

#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/CellLabel.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"
#include "EMCALBase/Geometry.h"
#include "EMCALBase/ClusterFactory.h"
#include "EMCALBase/NonlinearityHandler.h"
#include "EMCALReconstruction/Clusterizer.h"
#include "PWGJE/Core/JetUtilities.h"
#include "TVector2.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using myGlobTracks = o2::soa::Join<o2::aod::FullTracks, o2::aod::TrackSelection>;
using bcEvSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;
using collEventSels = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
using filteredCells = o2::soa::Filtered<aod::Calos>;
using mcCells = o2::soa::Join<aod::Calos, aod::McCaloLabels_001>;
using filteredMCCells = o2::soa::Filtered<mcCells>;

struct EmcalCorrectionTask {
  Produces<o2::aod::EMCALClusters> clusters;
  Produces<o2::aod::EMCALMCClusters> mcclusters;
  Produces<o2::aod::EMCALAmbiguousClusters> clustersAmbiguous;
  Produces<o2::aod::EMCALClusterCells> clustercells; // cells belonging to given cluster
  Produces<o2::aod::EMCALAmbiguousClusterCells> clustercellsambiguous;
  Produces<o2::aod::EMCALMatchedTracks> matchedTracks;
  Produces<o2::aod::EMCALMatchedCollisions> emcalcollisionmatch;

  // Preslices
  Preslice<myGlobTracks> perCollision = o2::aod::track::collisionId;
  PresliceUnsorted<collEventSels> collisionsPerFoundBC = aod::evsel::foundBCId;
  Preslice<aod::Collisions> collisionsPerBC = aod::collision::bcId;
  Preslice<filteredCells> cellsPerFoundBC = aod::calo::bcId;
  Preslice<filteredMCCells> mcCellsPerFoundBC = aod::calo::bcId;

  // Options for the clusterization
  // 1 corresponds to EMCAL cells based on the Run2 definition.
  Configurable<int> selectedCellType{"selectedCellType", 1, "EMCAL Cell type"};
  Configurable<std::string> clusterDefinitions{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default. Multiple definitions can be specified separated by comma"};
  Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.4f, "Max matching distance track-cluster"};
  Configurable<std::string> nonlinearityFunction{"nonlinearityFunction", "DATA_TestbeamFinal", "Nonlinearity correction at cluster level"};
  Configurable<bool> disableNonLin{"disableNonLin", false, "Disable NonLin correction if set to true"};
  Configurable<bool> hasShaperCorrection{"hasShaperCorrection", true, "Apply correction for shaper saturation"};
  Configurable<int> applyCellAbsScale{"applyCellAbsScale", 0, "Enable absolute cell energy scale to correct for energy loss in material in front of EMCal"};
  Configurable<std::vector<float>> vCellAbsScaleFactor{"cellAbsScaleFactor", {1.f}, "values for absolute cell energy calibration. Different values correspond to different regions or SM types of EMCal"};
  Configurable<float> logWeight{"logWeight", 4.5, "logarithmic weight for the cluster center of gravity calculation"};
  Configurable<float> exoticCellFraction{"exoticCellFraction", 0.97, "Good cell if fraction < 1-ecross/ecell"};
  Configurable<float> exoticCellDiffTime{"exoticCellDiffTime", 1.e6, "If time of candidate to exotic and close cell is larger than exoticCellDiffTime (in ns), it must be noisy, set amp to 0"};
  Configurable<float> exoticCellMinAmplitude{"exoticCellMinAmplitude", 4, "Check for exotic only if amplitud is larger than this value"};
  Configurable<float> exoticCellInCrossMinAmplitude{"exoticCellInCrossMinAmplitude", 0.1, "Minimum energy of cells in cross, if lower not considered in cross"};
  Configurable<bool> useWeightExotic{"useWeightExotic", false, "States if weights should be used for exotic cell cut"};
  Configurable<bool> isMC{"isMC", false, "States if run over MC"};
  Configurable<int> applyCellTimeShift{"applyCellTimeShift", 0, "apply shift to the cell time for data and MC; For data: 0 = off; non-zero = log function extracted from data - For MC: 0 = off; 1 = const shift; 2 = eta-dependent shift"};

  // Require EMCAL cells (CALO type 1)
  Filter emccellfilter = aod::calo::caloType == selectedCellType;

  // CDB service (for geometry)
  Service<o2::ccdb::BasicCCDBManager> mCcdbManager;

  // Clusterizer and related
  // Apparently streaming these objects really doesn't work, and causes problems for setting up the workflow.
  // So we use unique_ptr and define them below.
  std::vector<std::unique_ptr<o2::emcal::Clusterizer<o2::emcal::Cell>>> mClusterizers;
  o2::emcal::ClusterFactory<o2::emcal::Cell> mClusterFactories;
  o2::emcal::NonlinearityHandler mNonlinearityHandler;
  // Cells and clusters
  std::vector<o2::emcal::AnalysisCluster> mAnalysisClusters;
  std::vector<o2::emcal::ClusterLabel> mClusterLabels;

  std::vector<o2::aod::EMCALClusterDefinition> mClusterDefinitions;
  // QA
  o2::framework::HistogramRegistry mHistManager{"EMCALCorrectionTaskQAHistograms"};

  // EMCal geometry
  o2::emcal::Geometry* geometry;

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
    geometry = o2::emcal::Geometry::GetInstanceFromRunNumber(223409);
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
    mClusterFactories.setGeometry(geometry);
    mClusterFactories.SetECALogWeight(logWeight);
    mClusterFactories.setExoticCellFraction(exoticCellFraction);
    mClusterFactories.setExoticCellDiffTime(exoticCellDiffTime);
    mClusterFactories.setExoticCellMinAmplitude(exoticCellMinAmplitude);
    mClusterFactories.setExoticCellInCrossMinAmplitude(exoticCellInCrossMinAmplitude);
    mClusterFactories.setUseWeightExotic(useWeightExotic);
    for (auto& clusterDefinition : mClusterDefinitions) {
      mClusterizers.emplace_back(std::make_unique<o2::emcal::Clusterizer<o2::emcal::Cell>>(1E9, clusterDefinition.timeMin, clusterDefinition.timeMax, clusterDefinition.gradientCut, clusterDefinition.doGradientCut, clusterDefinition.seedEnergy, clusterDefinition.minCellEnergy));
      LOG(info) << "Cluster definition initialized: " << clusterDefinition.toString();
      LOG(info) << "timeMin: " << clusterDefinition.timeMin;
      LOG(info) << "timeMax: " << clusterDefinition.timeMax;
      LOG(info) << "gradientCut: " << clusterDefinition.gradientCut;
      LOG(info) << "seedEnergy: " << clusterDefinition.seedEnergy;
      LOG(info) << "minCellEnergy: " << clusterDefinition.minCellEnergy;
      LOG(info) << "storageID" << clusterDefinition.storageID;
    }
    for (auto& clusterizer : mClusterizers) {
      clusterizer->setGeometry(geometry);
    }

    if (mClusterizers.size() == 0) {
      LOG(error) << "No cluster definitions specified!";
    }

    mNonlinearityHandler = o2::emcal::NonlinearityFactory::getInstance().getNonlinearity(static_cast<std::string>(nonlinearityFunction));
    LOG(info) << "Using nonlinearity parameterisation: " << nonlinearityFunction.value;
    LOG(info) << "Apply shaper saturation correction:  " << (hasShaperCorrection.value ? "yes" : "no");

    LOG(debug) << "Completed init!";

    // Setup QA hists.
    // NOTE: This is not comprehensive.
    using o2HistType = o2::framework::HistType;
    using o2Axis = o2::framework::AxisSpec;
    o2Axis energyAxis{200, 0., 100., "E (GeV)"},
      timeAxis{300, -100, 200., "t (ns)"},
      etaAxis{160, -0.8, 0.8, "#eta"},
      phiAxis{72, 0, 2 * 3.14159, "phi"};
    mHistManager.add("hCellE", "hCellE", o2HistType::kTH1F, {energyAxis});
    mHistManager.add("hCellTowerID", "hCellTowerID", o2HistType::kTH1D, {{20000, 0, 20000}});
    mHistManager.add("hCellEtaPhi", "hCellEtaPhi", o2HistType::kTH2F, {etaAxis, phiAxis});
    // NOTE: Reversed column and row because it's more natural for presentation.
    mHistManager.add("hCellRowCol", "hCellRowCol;Column;Row", o2HistType::kTH2D, {{97, 0, 97}, {600, 0, 600}});
    mHistManager.add("hClusterE", "hClusterE", o2HistType::kTH1F, {energyAxis});
    mHistManager.add("hClusterEtaPhi", "hClusterEtaPhi", o2HistType::kTH2F, {etaAxis, phiAxis});
    mHistManager.add("hClusterTime", "hClusterTime", o2HistType::kTH1F, {timeAxis});
    mHistManager.add("hGlobalTrackEtaPhi", "hGlobalTrackEtaPhi", o2HistType::kTH2F, {etaAxis, phiAxis});
    mHistManager.add("hGlobalTrackMult", "hGlobalTrackMult", o2HistType::kTH1D, {{200, -0.5, 199.5, "N_{trk}"}});
    mHistManager.add("hCollisionType", "hCollisionType;;#it{count}", o2HistType::kTH1D, {{3, -0.5, 2.5}});
    auto hCollisionType = mHistManager.get<TH1>(HIST("hCollisionType"));
    hCollisionType->GetXaxis()->SetBinLabel(1, "no collision");
    hCollisionType->GetXaxis()->SetBinLabel(2, "normal collision");
    hCollisionType->GetXaxis()->SetBinLabel(3, "mult. collisions");
    mHistManager.add("hClusterType", "hClusterType;;#it{count}", o2HistType::kTH1D, {{3, -0.5, 2.5}});
    auto hClusterType = mHistManager.get<TH1>(HIST("hClusterType"));
    hClusterType->GetXaxis()->SetBinLabel(1, "no collision");
    hClusterType->GetXaxis()->SetBinLabel(2, "normal collision");
    hClusterType->GetXaxis()->SetBinLabel(3, "mult. collisions");
    mHistManager.add("hCollPerBC", "hCollPerBC;#it{N}_{coll.};#it{count}", o2HistType::kTH1D, {{100, -0.5, 99.5}});
    mHistManager.add("hBC", "hBC;;#it{count}", o2HistType::kTH1D, {{8, -0.5, 7.5}});
    mHistManager.add("hCollisionTimeReso", "hCollisionTimeReso;#Delta t_{coll};#it{count}", o2HistType::kTH1D, {{2000, 0, 2000}});
    auto hBC = mHistManager.get<TH1>(HIST("hBC"));
    hBC->GetXaxis()->SetBinLabel(1, "with EMCal cells");
    hBC->GetXaxis()->SetBinLabel(2, "with EMCal cells but no collision");
    hBC->GetXaxis()->SetBinLabel(3, "with EMCal cells and collision");
    hBC->GetXaxis()->SetBinLabel(4, "with EMCal cells and mult. collisions");
    hBC->GetXaxis()->SetBinLabel(5, "no EMCal cells and no collision");
    hBC->GetXaxis()->SetBinLabel(6, "no EMCal cells and with collision");
    hBC->GetXaxis()->SetBinLabel(7, "no EMCal cells and mult. collisions");
    hBC->GetXaxis()->SetBinLabel(8, "all BC");
    if (isMC) {
      mHistManager.add("hContributors", "hContributors;contributor per cell hit;#it{counts}", o2HistType::kTH1I, {{20, 0, 20}});
      mHistManager.add("hMCParticleEnergy", "hMCParticleEnergy;#it{E} (GeV/#it{c});#it{counts}", o2HistType::kTH1F, {energyAxis});
    }
  }

  // void process(aod::Collision const& collision, soa::Filtered<aod::Tracks> const& fullTracks, aod::Calos const& cells)
  // void process(aod::Collision const& collision, aod::Tracks const& tracks, aod::Calos const& cells)
  // void process(aod::BCs const& bcs, aod::Collision const& collision, aod::Calos const& cells)

  //  Appears to need the BC to be accessed to be available in the collision table...
  void processFull(bcEvSels const& bcs, collEventSels const& collisions, myGlobTracks const& tracks, filteredCells const& cells)
  {
    LOG(debug) << "Starting process full.";

    int nBCsProcessed = 0;
    int nCellsProcessed = 0;
    std::unordered_map<uint64_t, int> numberCollsInBC; // Number of collisions mapped to the global BC index of all BCs
    std::unordered_map<uint64_t, int> numberCellsInBC; // Number of cells mapped to the global BC index of all BCs to check whether EMCal was readout
    for (auto bc : bcs) {
      LOG(debug) << "Next BC";
      // Convert aod::Calo to o2::emcal::Cell which can be used with the clusterizer.
      // In particular, we need to filter only EMCAL cells.

      // Get the collisions matched to the BC using foundBCId of the collision
      auto collisionsInFoundBC = collisions.sliceBy(collisionsPerFoundBC, bc.globalIndex());
      auto cellsInBC = cells.sliceBy(cellsPerFoundBC, bc.globalIndex());

      numberCollsInBC.insert(std::pair<uint64_t, int>(bc.globalIndex(), collisionsInFoundBC.size()));
      numberCellsInBC.insert(std::pair<uint64_t, int>(bc.globalIndex(), cellsInBC.size()));

      if (!cellsInBC.size()) {
        LOG(debug) << "No cells found for BC";
        countBC(collisionsInFoundBC.size(), false);
        continue;
      }
      // Counters for BCs with matched collisions
      countBC(collisionsInFoundBC.size(), true);
      std::vector<o2::emcal::Cell> cellsBC;
      std::vector<int64_t> cellIndicesBC;
      for (auto& cell : cellsInBC) {
        auto amplitude = cell.amplitude();
        if (static_cast<bool>(hasShaperCorrection)) {
          amplitude = o2::emcal::NonlinearityHandler::evaluateShaperCorrectionCellEnergy(amplitude);
        }
        if (applyCellAbsScale) {
          amplitude *= GetAbsCellScale(cell.cellNumber());
        }
        cellsBC.emplace_back(cell.cellNumber(),
                             amplitude,
                             cell.time() + getCellTimeShift(cell.cellNumber(), amplitude),
                             o2::emcal::intToChannelType(cell.cellType()));
        cellIndicesBC.emplace_back(cell.globalIndex());
      }
      LOG(detail) << "Number of cells for BC (CF): " << cellsBC.size();
      nCellsProcessed += cellsBC.size();

      fillQAHistogram(cellsBC);

      // TODO: Helpful for now, but should be removed.
      LOG(debug) << "Converted EMCAL cells";
      for (auto& cell : cellsBC) {
        LOG(debug) << cell.getTower() << ": E: " << cell.getEnergy() << ", time: " << cell.getTimeStamp() << ", type: " << cell.getType();
      }

      LOG(debug) << "Converted cells. Contains: " << cellsBC.size() << ". Originally " << cellsInBC.size() << ". About to run clusterizer.";
      //  this is a test
      //  Run the clusterizers
      LOG(debug) << "Running clusterizers";
      for (size_t iClusterizer = 0; iClusterizer < mClusterizers.size(); iClusterizer++) {
        cellsToCluster(iClusterizer, cellsBC);

        if (collisionsInFoundBC.size() == 1) {
          // dummy loop to get the first collision
          for (const auto& col : collisionsInFoundBC) {
            if (col.foundBCId() == bc.globalIndex()) {
              mHistManager.fill(HIST("hCollisionTimeReso"), col.collisionTimeRes());
              mHistManager.fill(HIST("hCollPerBC"), 1);
              mHistManager.fill(HIST("hCollisionType"), 1);
              math_utils::Point3D<float> vertex_pos = {col.posX(), col.posY(), col.posZ()};

              std::vector<std::vector<int>> clusterToTrackIndexMap;
              std::vector<std::vector<int>> trackToClusterIndexMap;
              std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>> IndexMapPair{clusterToTrackIndexMap, trackToClusterIndexMap};
              std::vector<int64_t> trackGlobalIndex;
              doTrackMatching<collEventSels::filtered_iterator>(col, tracks, IndexMapPair, vertex_pos, trackGlobalIndex);

              // Store the clusters in the table where a matching collision could
              // be identified.
              FillClusterTable<collEventSels::filtered_iterator>(col, vertex_pos, iClusterizer, cellIndicesBC, IndexMapPair, trackGlobalIndex);
            }
          }
        } else { // ambiguous
          // LOG(warning) << "No vertex found for event. Assuming (0,0,0).";
          bool hasCollision = false;
          mHistManager.fill(HIST("hCollPerBC"), collisionsInFoundBC.size());
          if (collisionsInFoundBC.size() == 0) {
            mHistManager.fill(HIST("hCollisionType"), 0);
          } else {
            hasCollision = true;
            mHistManager.fill(HIST("hCollisionType"), 2);
          }
          FillAmbigousClusterTable<bcEvSels::iterator>(bc, iClusterizer, cellIndicesBC, hasCollision);
        }

        LOG(debug) << "Cluster loop done for clusterizer " << iClusterizer;
      } // end of clusterizer loop
      LOG(debug) << "Done with process BC.";
      nBCsProcessed++;
    } // end of bc loop

    // Loop through all collisions and fill emcalcollisionmatch with a boolean stating, whether the collision was ambiguous (not the only collision in its BC)
    for (const auto& collision : collisions) {
      auto globalbcid = collision.foundBC_as<bcEvSels>().globalIndex();
      auto foundColls = numberCollsInBC.find(globalbcid);
      auto foundCells = numberCellsInBC.find(globalbcid);
      if (foundColls != numberCollsInBC.end() && foundCells != numberCellsInBC.end()) {
        emcalcollisionmatch(collision.globalIndex(), foundColls->second != 1, foundCells->second > 0);
      } else {
        LOG(warning) << "BC not found in map of number of collisions.";
      }
    } // end of collision loop

    LOG(detail) << "Processed " << nBCsProcessed << " BCs with " << nCellsProcessed << " cells";
  }
  PROCESS_SWITCH(EmcalCorrectionTask, processFull, "run full analysis", true);

  void processMCFull(bcEvSels const& bcs, collEventSels const& collisions, myGlobTracks const& tracks, filteredMCCells const& cells, aod::StoredMcParticles_001 const&)
  {
    LOG(debug) << "Starting process full.";

    int nBCsProcessed = 0;
    int nCellsProcessed = 0;
    std::unordered_map<uint64_t, int> numberCollsInBC; // Number of collisions mapped to the global BC index of all BCs
    std::unordered_map<uint64_t, int> numberCellsInBC; // Number of cells mapped to the global BC index of all BCs to check whether EMCal was readout
    for (auto bc : bcs) {
      LOG(debug) << "Next BC";
      // Convert aod::Calo to o2::emcal::Cell which can be used with the clusterizer.
      // In particular, we need to filter only EMCAL cells.

      // Get the collisions matched to the BC using foundBCId of the collision
      auto collisionsInFoundBC = collisions.sliceBy(collisionsPerFoundBC, bc.globalIndex());
      auto cellsInBC = cells.sliceBy(mcCellsPerFoundBC, bc.globalIndex());

      numberCollsInBC.insert(std::pair<uint64_t, int>(bc.globalIndex(), collisionsInFoundBC.size()));
      numberCellsInBC.insert(std::pair<uint64_t, int>(bc.globalIndex(), cellsInBC.size()));

      if (!cellsInBC.size()) {
        LOG(debug) << "No cells found for BC";
        countBC(collisionsInFoundBC.size(), false);
        continue;
      }
      // Counters for BCs with matched collisions
      countBC(collisionsInFoundBC.size(), true);
      std::vector<o2::emcal::Cell> cellsBC;
      std::vector<int64_t> cellIndicesBC;
      std::vector<o2::emcal::CellLabel> cellLabels;
      for (auto& cell : cellsInBC) {
        mHistManager.fill(HIST("hContributors"), cell.mcParticle_as<aod::StoredMcParticles_001>().size());
        auto cellParticles = cell.mcParticle_as<aod::StoredMcParticles_001>();
        for (auto& cellparticle : cellParticles) {
          mHistManager.fill(HIST("hMCParticleEnergy"), cellparticle.e());
        }
        auto amplitude = cell.amplitude();
        if (static_cast<bool>(hasShaperCorrection)) {
          amplitude = o2::emcal::NonlinearityHandler::evaluateShaperCorrectionCellEnergy(amplitude);
        }
        cellsBC.emplace_back(cell.cellNumber(),
                             amplitude,
                             cell.time() + getCellTimeShift(cell.cellNumber(), amplitude),
                             o2::emcal::intToChannelType(cell.cellType()));
        cellIndicesBC.emplace_back(cell.globalIndex());
        cellLabels.emplace_back(cell.mcParticleIds(), cell.amplitudeA());
      }
      LOG(detail) << "Number of cells for BC (CF): " << cellsBC.size();
      nCellsProcessed += cellsBC.size();

      fillQAHistogram(cellsBC);

      // TODO: Helpful for now, but should be removed.
      LOG(debug) << "Converted EMCAL cells";
      for (auto& cell : cellsBC) {
        LOG(debug) << cell.getTower() << ": E: " << cell.getEnergy() << ", time: " << cell.getTimeStamp() << ", type: " << cell.getType();
      }

      LOG(debug) << "Converted cells. Contains: " << cellsBC.size() << ". Originally " << cellsInBC.size() << ". About to run clusterizer.";
      //  this is a test
      //  Run the clusterizers
      LOG(debug) << "Running clusterizers";
      for (size_t iClusterizer = 0; iClusterizer < mClusterizers.size(); iClusterizer++) {
        cellsToCluster(iClusterizer, cellsBC, cellLabels);

        if (collisionsInFoundBC.size() == 1) {
          // dummy loop to get the first collision
          for (const auto& col : collisionsInFoundBC) {
            if (col.foundBCId() == bc.globalIndex()) {
              mHistManager.fill(HIST("hCollPerBC"), 1);
              mHistManager.fill(HIST("hCollisionType"), 1);
              math_utils::Point3D<float> vertex_pos = {col.posX(), col.posY(), col.posZ()};

              std::vector<std::vector<int>> clusterToTrackIndexMap;
              std::vector<std::vector<int>> trackToClusterIndexMap;
              std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>> IndexMapPair{clusterToTrackIndexMap, trackToClusterIndexMap};
              std::vector<int64_t> trackGlobalIndex;
              doTrackMatching<collEventSels::filtered_iterator>(col, tracks, IndexMapPair, vertex_pos, trackGlobalIndex);

              // Store the clusters in the table where a matching collision could
              // be identified.
              FillClusterTable<collEventSels::filtered_iterator>(col, vertex_pos, iClusterizer, cellIndicesBC, IndexMapPair, trackGlobalIndex);
            }
          }
        } else { // ambiguous
          // LOG(warning) << "No vertex found for event. Assuming (0,0,0).";
          bool hasCollision = false;
          mHistManager.fill(HIST("hCollPerBC"), collisionsInFoundBC.size());
          if (collisionsInFoundBC.size() == 0) {
            mHistManager.fill(HIST("hCollisionType"), 0);
          } else {
            hasCollision = true;
            mHistManager.fill(HIST("hCollisionType"), 2);
          }
          FillAmbigousClusterTable<bcEvSels::iterator>(bc, iClusterizer, cellIndicesBC, hasCollision);
        }
        LOG(debug) << "Cluster loop done for clusterizer " << iClusterizer;
      } // end of clusterizer loop
      LOG(debug) << "Done with process BC.";
      nBCsProcessed++;
    } // end of bc loop

    // Loop through all collisions and fill emcalcollisionmatch with a boolean stating, whether the collision was ambiguous (not the only collision in its BC)
    for (const auto& collision : collisions) {
      auto globalbcid = collision.foundBC_as<bcEvSels>().globalIndex();
      auto foundColls = numberCollsInBC.find(globalbcid);
      auto foundCells = numberCellsInBC.find(globalbcid);
      if (foundColls != numberCollsInBC.end() && foundCells != numberCellsInBC.end()) {
        emcalcollisionmatch(collision.globalIndex(), foundColls->second != 1, foundCells->second > 0);
      } else {
        LOG(warning) << "BC not found in map of number of collisions.";
      }
    } // end of collision loop

    LOG(detail) << "Processed " << nBCsProcessed << " BCs with " << nCellsProcessed << " cells";
  }
  PROCESS_SWITCH(EmcalCorrectionTask, processMCFull, "run full analysis with MC info", false);
  void processStandalone(aod::BCs const& bcs, aod::Collisions const& collisions, filteredCells const& cells)
  {
    LOG(debug) << "Starting process standalone.";
    int nBCsProcessed = 0;
    int nCellsProcessed = 0;
    for (auto bc : bcs) {
      LOG(debug) << "Next BC";
      // Convert aod::Calo to o2::emcal::Cell which can be used with the clusterizer.
      // In particular, we need to filter only EMCAL cells.

      // Get the collisions matched to the BC using global bc index of the collision
      // since we do not have event selection available here!
      auto collisionsInBC = collisions.sliceBy(collisionsPerBC, bc.globalIndex());
      auto cellsInBC = cells.sliceBy(cellsPerFoundBC, bc.globalIndex());

      if (!cellsInBC.size()) {
        LOG(debug) << "No cells found for BC";
        countBC(collisionsInBC.size(), false);
        continue;
      }
      // Counters for BCs with matched collisions
      countBC(collisionsInBC.size(), true);
      std::vector<o2::emcal::Cell> cellsBC;
      std::vector<int64_t> cellIndicesBC;
      for (auto& cell : cellsInBC) {
        cellsBC.emplace_back(cell.cellNumber(),
                             cell.amplitude(),
                             cell.time() + getCellTimeShift(cell.cellNumber(), cell.amplitude()),
                             o2::emcal::intToChannelType(cell.cellType()));
        cellIndicesBC.emplace_back(cell.globalIndex());
      }
      LOG(detail) << "Number of cells for BC (CF): " << cellsBC.size();
      nCellsProcessed += cellsBC.size();

      fillQAHistogram(cellsBC);

      // TODO: Helpful for now, but should be removed.
      LOG(debug) << "Converted EMCAL cells";
      for (auto& cell : cellsBC) {
        LOG(debug) << cell.getTower() << ": E: " << cell.getEnergy() << ", time: " << cell.getTimeStamp() << ", type: " << cell.getType();
      }

      LOG(debug) << "Converted cells. Contains: " << cellsBC.size() << ". Originally " << cellsInBC.size() << ". About to run clusterizer.";

      //  this is a test
      //  Run the clusterizers
      LOG(debug) << "Running clusterizers";
      for (size_t iClusterizer = 0; iClusterizer < mClusterizers.size(); iClusterizer++) {
        cellsToCluster(iClusterizer, cellsBC);

        if (collisionsInBC.size() == 1) {
          // dummy loop to get the first collision
          for (const auto& col : collisionsInBC) {
            mHistManager.fill(HIST("hCollPerBC"), 1);
            mHistManager.fill(HIST("hCollisionType"), 1);
            math_utils::Point3D<float> vertex_pos = {col.posX(), col.posY(), col.posZ()};

            // Store the clusters in the table where a matching collision could
            // be identified.
            FillClusterTable<aod::Collision>(col, vertex_pos, iClusterizer, cellIndicesBC);
          }
        } else { // ambiguous
          // LOG(warning) << "No vertex found for event. Assuming (0,0,0).";
          bool hasCollision = false;
          mHistManager.fill(HIST("hCollPerBC"), collisionsInBC.size());
          if (collisionsInBC.size() == 0) {
            mHistManager.fill(HIST("hCollisionType"), 0);
          } else {
            hasCollision = true;
            mHistManager.fill(HIST("hCollisionType"), 2);
          }
          FillAmbigousClusterTable<aod::BC>(bc, iClusterizer, cellIndicesBC, hasCollision);
        }

        LOG(debug) << "Cluster loop done for clusterizer " << iClusterizer;
      } // end of clusterizer loop
      LOG(detail) << "Processed " << nBCsProcessed << " BCs with " << nCellsProcessed << " cells";
      nBCsProcessed++;
    } // end of bc loop
    LOG(debug) << "Done with process BC.";
  }
  PROCESS_SWITCH(EmcalCorrectionTask, processStandalone, "run stand alone analysis", false);

  void cellsToCluster(size_t iClusterizer, const gsl::span<o2::emcal::Cell> cellsBC, std::optional<const gsl::span<o2::emcal::CellLabel>> cellLabels = std::nullopt)
  {
    mClusterizers.at(iClusterizer)->findClusters(cellsBC);

    auto emcalClusters = mClusterizers.at(iClusterizer)->getFoundClusters();
    auto emcalClustersInputIndices = mClusterizers.at(iClusterizer)->getFoundClustersInputIndices();
    LOG(debug) << "Retrieved results. About to setup cluster factory.";

    // Convert to analysis clusters.
    // First, the cluster factory requires cluster and cell information in order
    // to build the clusters.
    mAnalysisClusters.clear();
    mClusterLabels.clear();
    mClusterFactories.reset();
    if (cellLabels) {
      mClusterFactories.setContainer(*emcalClusters, cellsBC, *emcalClustersInputIndices, cellLabels);
    } else {
      mClusterFactories.setContainer(*emcalClusters, cellsBC, *emcalClustersInputIndices);
    }

    LOG(debug) << "Cluster factory set up.";
    // Convert to analysis clusters.
    for (int icl = 0; icl < mClusterFactories.getNumberOfClusters(); icl++) {
      o2::emcal::ClusterLabel clusterLabel;
      auto analysisCluster = mClusterFactories.buildCluster(icl, &clusterLabel);
      mAnalysisClusters.emplace_back(analysisCluster);
      mClusterLabels.push_back(clusterLabel);
      LOG(debug) << "Cluster " << icl << ": E: " << analysisCluster.E()
                 << ", NCells " << analysisCluster.getNCells();
    }
    LOG(debug) << "Converted to analysis clusters.";
  }

  template <typename Collision>
  void FillClusterTable(Collision const& col, math_utils::Point3D<float> const& vertex_pos, size_t iClusterizer, const gsl::span<int64_t> cellIndicesBC, std::optional<std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>>> const& IndexMapPair = std::nullopt, std::optional<std::vector<int64_t>> const& trackGlobalIndex = std::nullopt)
  {
    // we found a collision, put the clusters into the none ambiguous table
    clusters.reserve(mAnalysisClusters.size());
    if (mClusterLabels.size() > 0) {
      mcclusters.reserve(mClusterLabels.size());
    }
    int cellindex = -1;
    unsigned int iCluster = 0;
    for (const auto& cluster : mAnalysisClusters) {
      // Determine the cluster eta, phi, correcting for the vertex position.
      auto pos = cluster.getGlobalPosition();
      pos = pos - vertex_pos;
      // Normalize the vector and rescale by energy.
      pos *= (cluster.E() / std::sqrt(pos.Mag2()));

      // Correct for nonlinear behaviour
      float nonlinCorrEnergy = cluster.E();
      if (!disableNonLin) {
        try {
          nonlinCorrEnergy = mNonlinearityHandler.getCorrectedClusterEnergy(cluster);
        } catch (o2::emcal::NonlinearityHandler::UninitException& e) {
          LOG(error) << e.what();
        }
      }

      // save to table
      LOG(debug) << "Writing cluster definition "
                 << static_cast<int>(mClusterDefinitions.at(iClusterizer))
                 << " to table.";
      mHistManager.fill(HIST("hClusterType"), 1);
      clusters(col, cluster.getID(), nonlinCorrEnergy, cluster.getCoreEnergy(), cluster.E(),
               pos.Eta(), TVector2::Phi_0_2pi(pos.Phi()), cluster.getM02(),
               cluster.getM20(), cluster.getNCells(),
               cluster.getClusterTime(), cluster.getIsExotic(),
               cluster.getDistanceToBadChannel(), cluster.getNExMax(),
               static_cast<int>(mClusterDefinitions.at(iClusterizer)));
      if (mClusterLabels.size() > 0) {
        mcclusters(mClusterLabels[iCluster].getLabels(), mClusterLabels[iCluster].getEnergyFractions());
      }
      clustercells.reserve(cluster.getNCells());
      // loop over cells in cluster and save to table
      for (int ncell = 0; ncell < cluster.getNCells(); ncell++) {
        cellindex = cluster.getCellIndex(ncell);
        LOG(debug) << "trying to find cell index " << cellindex << " in map";
        clustercells(clusters.lastIndex(), cellIndicesBC[cellindex]);
      } // end of cells of cluser loop
      // fill histograms
      mHistManager.fill(HIST("hClusterE"), cluster.E());
      mHistManager.fill(HIST("hClusterTime"), cluster.getClusterTime());
      mHistManager.fill(HIST("hClusterEtaPhi"), pos.Eta(), TVector2::Phi_0_2pi(pos.Phi()));
      if (IndexMapPair && trackGlobalIndex) {
        for (unsigned int iTrack = 0; iTrack < std::get<0>(*IndexMapPair)[iCluster].size(); iTrack++) {
          if (std::get<0>(*IndexMapPair)[iCluster][iTrack] >= 0) {
            LOG(debug) << "Found track " << (*trackGlobalIndex)[std::get<0>(*IndexMapPair)[iCluster][iTrack]] << " in cluster " << cluster.getID();
            matchedTracks(clusters.lastIndex(), (*trackGlobalIndex)[std::get<0>(*IndexMapPair)[iCluster][iTrack]]);
          }
        }
      }
      iCluster++;
    } // end of cluster loop
  }

  template <typename BC>
  void FillAmbigousClusterTable(BC const& bc, size_t iClusterizer, const gsl::span<int64_t> cellIndicesBC, bool hasCollision)
  {
    int cellindex = -1;
    clustersAmbiguous.reserve(mAnalysisClusters.size());
    for (const auto& cluster : mAnalysisClusters) {
      auto pos = cluster.getGlobalPosition();
      pos = pos - math_utils::Point3D<float>{0., 0., 0.};
      // Normalize the vector and rescale by energy.
      pos *= (cluster.E() / std::sqrt(pos.Mag2()));

      // Correct for nonlinear behaviour
      float nonlinCorrEnergy = cluster.E();
      try {
        nonlinCorrEnergy = mNonlinearityHandler.getCorrectedClusterEnergy(cluster);
      } catch (o2::emcal::NonlinearityHandler::UninitException& e) {
        LOG(error) << e.what();
      }

      // We have our necessary properties. Now we store outputs

      // LOG(debug) << "Cluster E: " << cluster.E();
      if (!hasCollision) {
        mHistManager.fill(HIST("hClusterType"), 0);
      } else {
        mHistManager.fill(HIST("hClusterType"), 2);
      }
      clustersAmbiguous(
        bc, cluster.getID(), nonlinCorrEnergy, cluster.getCoreEnergy(), cluster.E(),
        pos.Eta(), TVector2::Phi_0_2pi(pos.Phi()), cluster.getM02(),
        cluster.getM20(), cluster.getNCells(), cluster.getClusterTime(),
        cluster.getIsExotic(), cluster.getDistanceToBadChannel(),
        cluster.getNExMax(), static_cast<int>(mClusterDefinitions.at(iClusterizer)));
      clustercellsambiguous.reserve(cluster.getNCells());
      for (int ncell = 0; ncell < cluster.getNCells(); ncell++) {
        cellindex = cluster.getCellIndex(ncell);
        clustercellsambiguous(clustersAmbiguous.lastIndex(),
                              cellIndicesBC[cellindex]);
      } // end of cells of cluster loop
    } // end of cluster loop
  }

  template <typename Collision>
  void doTrackMatching(Collision const& col, myGlobTracks const& tracks, std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>>& IndexMapPair, math_utils::Point3D<float>& vertex_pos, std::vector<int64_t>& trackGlobalIndex)
  {
    auto groupedTracks = tracks.sliceBy(perCollision, col.globalIndex());
    int NTracksInCol = groupedTracks.size();
    std::vector<double> trackPhi;
    std::vector<double> trackEta;
    // reserve memory to reduce on the fly memory allocation
    trackPhi.reserve(NTracksInCol);
    trackEta.reserve(NTracksInCol);
    trackGlobalIndex.reserve(NTracksInCol);
    FillTrackInfo<decltype(groupedTracks)>(groupedTracks, trackPhi, trackEta, trackGlobalIndex);

    int NClusterInCol = mAnalysisClusters.size();
    std::vector<double> clusterPhi;
    std::vector<double> clusterEta;
    clusterPhi.reserve(NClusterInCol);
    clusterEta.reserve(NClusterInCol);

    // TODO one loop that could in principle be combined with the other
    // loop to improve performance
    for (const auto& cluster : mAnalysisClusters) {
      // Determine the cluster eta, phi, correcting for the vertex
      // position.
      auto pos = cluster.getGlobalPosition();
      pos = pos - vertex_pos;
      // Normalize the vector and rescale by energy.
      pos *= (cluster.E() / std::sqrt(pos.Mag2()));
      clusterPhi.emplace_back(TVector2::Phi_0_2pi(pos.Phi()));
      clusterEta.emplace_back(pos.Eta());
    }
    IndexMapPair =
      jetutilities::MatchClustersAndTracks(clusterPhi, clusterEta,
                                           trackPhi, trackEta,
                                           maxMatchingDistance, 20);
  }

  template <typename Tracks>
  void FillTrackInfo(Tracks const& tracks, std::vector<double>& trackPhi, std::vector<double>& trackEta, std::vector<int64_t>& trackGlobalIndex)
  {
    int NTrack = 0;
    for (auto& track : tracks) {
      // TODO only consider tracks in current emcal/dcal acceptanc
      if (!track.isGlobalTrack()) { // only global tracks
        continue;
      }
      NTrack++;
      trackPhi.emplace_back(TVector2::Phi_0_2pi(track.trackPhiEmcal()));
      trackEta.emplace_back(track.trackEtaEmcal());
      mHistManager.fill(HIST("hGlobalTrackEtaPhi"), track.trackEtaEmcal(),
                        TVector2::Phi_0_2pi(track.trackPhiEmcal()));
      trackGlobalIndex.emplace_back(track.globalIndex());
    }
    mHistManager.fill(HIST("hGlobalTrackMult"), NTrack);
  }

  void countBC(int numberOfCollisions, bool hasEMCcells)
  {
    int emcDataOffset = hasEMCcells ? 0 : 3;
    int collisionOffset = 2;
    switch (numberOfCollisions) {
      case 0:
        collisionOffset = 0;
        break;
      case 1:
        collisionOffset = 1;
        break;
      default:
        collisionOffset = 2;
        break;
    }
    mHistManager.fill(HIST("hBC"), 7); // All collisions
    if (hasEMCcells) {
      mHistManager.fill(HIST("hBC"), 0);
    }
    mHistManager.fill(HIST("hBC"), 1 + emcDataOffset + collisionOffset);
  }

  void fillQAHistogram(const gsl::span<o2::emcal::Cell> cellsBC)
  {
    // Cell QA
    // For convenience, use the clusterizer stored geometry to get the eta-phi
    for (auto& cell : cellsBC) {
      mHistManager.fill(HIST("hCellE"), cell.getEnergy());
      mHistManager.fill(HIST("hCellTowerID"), cell.getTower());
      auto res = mClusterizers.at(0)->getGeometry()->EtaPhiFromIndex(cell.getTower());
      mHistManager.fill(HIST("hCellEtaPhi"), std::get<0>(res), TVector2::Phi_0_2pi(std::get<1>(res)));
      res = mClusterizers.at(0)->getGeometry()->GlobalRowColFromIndex(cell.getTower());
      // NOTE: Reversed column and row because it's more natural for presentation.
      mHistManager.fill(HIST("hCellRowCol"), std::get<1>(res), std::get<0>(res));
    }
  }

  float GetAbsCellScale(const int cellID)
  {
    // Apply cell scale based on SM types (Full, Half (not used), EMC 1/3, DCal, DCal 1/3)
    // Same as in Run2 data
    if (applyCellAbsScale == 1) {
      int iSM = mClusterizers.at(0)->getGeometry()->GetSuperModuleNumber(cellID);
      return vCellAbsScaleFactor.value[mClusterizers.at(0)->getGeometry()->GetSMType(iSM)];

      // Apply cell scale based on columns to accoutn for material of TRD structures
    } else if (applyCellAbsScale == 2) {
      auto res = mClusterizers.at(0)->getGeometry()->GlobalRowColFromIndex(cellID);
      return vCellAbsScaleFactor.value[std::get<1>(res)];
    } else {
      return 1.f;
    }
  }

  // Apply shift of the cell time in data and MC
  // In MC this has to be done to shift the cell time, which is not calibrated to 0 due to the flight time of the particles to the EMCal surface (~15ns)
  // In data this is done to correct for the time walk effect
  float getCellTimeShift(const int16_t cellID, const float cellEnergy)
  {
    if (isMC) {
      if (applyCellTimeShift == 1) { // constant shift
        LOG(debug) << "shift the cell time by 15ns";
        return -15.f;                       // roughly calculated by assuming particles travel with v=c (photons) and EMCal is 4.4m away from vertex
      } else if (applyCellTimeShift == 2) { // eta dependent shift ( as larger eta values are further away from collision point)
        // Use distance between vertex and EMCal (at eta = 0) and distance on EMCal surface (cell size times column) to calculate distance to cell
        // 0.2 is cell size in m (0.06) divided by the speed of light in m/ns (0.3)
        // 47.5 is the "middle" of the EMCal (2*48 cells in one column)
        float timeCol = 0.2f * (geometry->GlobalCol(cellID) - 47.5f); // calculate time to get to specific column
        float time = -sqrt(215.f + timeCol * timeCol);                // 215 is 14.67ns^2 (time it takes to get the cell at eta = 0)
        LOG(debug) << "shift the cell time by " << time << "  applyCellTimeShift " << applyCellTimeShift;
        return time;
      } else {
        return 0.f;
      }
    } else { // data
      if (applyCellTimeShift != 0) {
        if (cellEnergy < 0.3) // Cells with tless than 300 MeV cannot be the leading cell in the cluster, so their time does not require precise calibration
          return 0.f;
        else if (cellEnergy < 4.)                                         // Low energy regime
          return (0.57284 + 0.82194 * TMath::Log(1.30651 * cellEnergy));  // Parameters extracted from LHC22o (pp), but also usable for other periods
        else                                                              // High energy regime
          return (-0.05858 + 1.50593 * TMath::Log(0.97591 * cellEnergy)); // Parameters extracted from LHC22o (pp), but also usable for other periods
      } else {                                                            // Dont apply cell time shift if applyCellTimeShift == 0
        return 0.f;
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<EmcalCorrectionTask>(cfgc, TaskName{"emcal-correction-task"})};
}
