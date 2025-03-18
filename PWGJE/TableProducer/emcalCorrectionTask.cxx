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
/// EMCAL Correction Task
///
/// \file emcalCorrectionTask.cxx
///
/// \brief Task that provides EMCal clusters and applies necessary corrections
///
/// \author Raymond Ehlers (raymond.ehlers@cern.ch) ORNL, Florian Jonas (florian.jonas@cern.ch)
///

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <string>
#include <tuple>
#include <vector>
#include <random>

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
using MyGlobTracks = o2::soa::Join<o2::aod::FullTracks, o2::aod::TrackSelection>;
using BcEvSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;
using CollEventSels = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
using FilteredCells = o2::soa::Filtered<aod::Calos>;
using McCells = o2::soa::Join<aod::Calos, aod::McCaloLabels_001>;
using FilteredMcCells = o2::soa::Filtered<McCells>;

struct EmcalCorrectionTask {
  Produces<o2::aod::EMCALClusters> clusters;
  Produces<o2::aod::EMCALMCClusters> mcclusters;
  Produces<o2::aod::EMCALAmbiguousClusters> clustersAmbiguous;
  Produces<o2::aod::EMCALClusterCells> clustercells; // cells belonging to given cluster
  Produces<o2::aod::EMCALAmbiguousClusterCells> clustercellsambiguous;
  Produces<o2::aod::EMCALMatchedTracks> matchedTracks;
  Produces<o2::aod::EMCALMatchedCollisions> emcalcollisionmatch;

  // Preslices
  Preslice<MyGlobTracks> perCollision = o2::aod::track::collisionId;
  PresliceUnsorted<CollEventSels> collisionsPerFoundBC = aod::evsel::foundBCId;
  Preslice<aod::Collisions> collisionsPerBC = aod::collision::bcId;
  Preslice<FilteredCells> cellsPerFoundBC = aod::calo::bcId;
  Preslice<FilteredMcCells> mcCellsPerFoundBC = aod::calo::bcId;

  // Options for the clusterization
  // 1 corresponds to EMCAL cells based on the Run2 definition.
  Configurable<int> selectedCellType{"selectedCellType", 1, "EMCAL Cell type"};
  Configurable<std::string> clusterDefinitions{"clusterDefinitions", "kV3Default", "cluster definition to be selected, e.g. V3Default. Multiple definitions can be specified separated by comma"};
  Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.4f, "Max matching distance track-cluster"};
  Configurable<std::string> nonlinearityFunction{"nonlinearityFunction", "DATA_TestbeamFinal_NoScale", "Nonlinearity correction at cluster level. Default for data should be DATA_TestbeamFinal_NoScale. Default for MC should be MC_TestbeamFinal."};
  Configurable<bool> disableNonLin{"disableNonLin", false, "Disable NonLin correction if set to true"};
  Configurable<bool> hasShaperCorrection{"hasShaperCorrection", true, "Apply correction for shaper saturation"};
  Configurable<int> applyCellAbsScale{"applyCellAbsScale", 0, "Enable absolute cell energy scale to correct for energy loss in material in front of EMCal"};
  Configurable<std::vector<float>> cellAbsScaleFactors{"cellAbsScaleFactors", {1.f}, "values for absolute cell energy calibration. Different values correspond to different regions or SM types of EMCal"};
  Configurable<float> logWeight{"logWeight", 4.5, "logarithmic weight for the cluster center of gravity calculation"};
  Configurable<float> exoticCellFraction{"exoticCellFraction", 0.97, "Good cell if fraction < 1-ecross/ecell"};
  Configurable<float> exoticCellDiffTime{"exoticCellDiffTime", 1.e6, "If time of candidate to exotic and close cell is larger than exoticCellDiffTime (in ns), it must be noisy, set amp to 0"};
  Configurable<float> exoticCellMinAmplitude{"exoticCellMinAmplitude", 4, "Check for exotic only if amplitud is larger than this value"};
  Configurable<float> exoticCellInCrossMinAmplitude{"exoticCellInCrossMinAmplitude", 0.1, "Minimum energy of cells in cross, if lower not considered in cross"};
  Configurable<bool> useWeightExotic{"useWeightExotic", false, "States if weights should be used for exotic cell cut"};
  Configurable<bool> isMC{"isMC", false, "States if run over MC"};
  Configurable<bool> applyCellTimeCorrection{"applyCellTimeCorrection", true, "apply a correction to the cell time for data and MC: Shift both average cell times to 0 and smear MC time distribution to fit data better. For MC requires isMC to be true"};
  Configurable<float> trackMinPt{"trackMinPt", 0.3, "Minimum pT for tracks to perform track matching, to reduce computing time. Tracks below a certain pT will be loopers anyway."};

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

  // Random number generator to draw cell time smearing for MC
  std::random_device rd{};
  std::mt19937_64 rdgen{rd()};
  std::normal_distribution<> normalgaus{0, 1}; // mean = 0, stddev = 1 (apply amplitude of smearing after drawing random for performance reasons)

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
    for (const auto& clusterDefinition : mClusterDefinitions) {
      mClusterizers.emplace_back(std::make_unique<o2::emcal::Clusterizer<o2::emcal::Cell>>(clusterDefinition.timeDiff, clusterDefinition.timeMin, clusterDefinition.timeMax, clusterDefinition.gradientCut, clusterDefinition.doGradientCut, clusterDefinition.seedEnergy, clusterDefinition.minCellEnergy));
      LOG(info) << "Cluster definition initialized: " << clusterDefinition.toString();
      LOG(info) << "timeMin: " << clusterDefinition.timeMin;
      LOG(info) << "timeMax: " << clusterDefinition.timeMax;
      LOG(info) << "timeDiff: " << clusterDefinition.timeDiff;
      LOG(info) << "gradientCut: " << clusterDefinition.gradientCut;
      LOG(info) << "seedEnergy: " << clusterDefinition.seedEnergy;
      LOG(info) << "minCellEnergy: " << clusterDefinition.minCellEnergy;
      LOG(info) << "storageID: " << clusterDefinition.storageID;
    }
    for (const auto& clusterizer : mClusterizers) {
      clusterizer->setGeometry(geometry);
    }

    if (mClusterizers.size() == 0) {
      LOG(error) << "No cluster definitions specified!";
    }

    mNonlinearityHandler = o2::emcal::NonlinearityFactory::getInstance().getNonlinearity(static_cast<std::string>(nonlinearityFunction));
    LOG(info) << "Using nonlinearity parameterisation: " << nonlinearityFunction.value;
    LOG(info) << "Apply shaper saturation correction:  " << (hasShaperCorrection.value ? "yes" : "no");

    LOG(debug) << "Completed init!";

    // Define the cell energy binning
    std::vector<double> cellEnergyBins;
    for (int i = 0; i < 51; i++)
      cellEnergyBins.emplace_back(0.1 * (i - 0) + 0.0); // from 0 to 5 GeV/c, every 0.1 GeV
    for (int i = 51; i < 76; i++)
      cellEnergyBins.emplace_back(0.2 * (i - 51) + 5.2); // from 5.2 to 10.0 GeV, every 0.2 GeV
    for (int i = 76; i < 166; i++)
      cellEnergyBins.emplace_back(1. * (i - 76) + 11.); // from 11.0 to 100. GeV, every 1 GeV

    // Setup QA hists.
    // NOTE: This is not comprehensive.
    using O2HistType = o2::framework::HistType;
    o2::framework::AxisSpec energyAxis{200, 0., 100., "E (GeV)"},
      timeAxis{300, -100, 200., "t (ns)"},
      etaAxis{160, -0.8, 0.8, "#eta"},
      phiAxis{72, 0, 2 * 3.14159, "phi"},
      nlmAxis{50, -0.5, 49.5, "NLM"};
    mHistManager.add("hCellE", "hCellE", O2HistType::kTH1F, {energyAxis});
    mHistManager.add("hCellTowerID", "hCellTowerID", O2HistType::kTH1D, {{20000, 0, 20000}});
    mHistManager.add("hCellEtaPhi", "hCellEtaPhi", O2HistType::kTH2F, {etaAxis, phiAxis});
    mHistManager.add("hHGCellTimeEnergy", "hCellTime", O2HistType::kTH2F, {{300, -30, 30}, cellEnergyBins}); // Cell time vs energy for high gain cells (low energies)
    mHistManager.add("hLGCellTimeEnergy", "hCellTime", O2HistType::kTH2F, {{300, -30, 30}, cellEnergyBins}); // Cell time vs energy for low gain cells (high energies)
    // NOTE: Reversed column and row because it's more natural for presentation.
    mHistManager.add("hCellRowCol", "hCellRowCol;Column;Row", O2HistType::kTH2D, {{96, -0.5, 95.5}, {208, -0.5, 207.5}});
    mHistManager.add("hClusterE", "hClusterE", O2HistType::kTH1F, {energyAxis});
    mHistManager.add("hClusterNLM", "hClusterNLM", O2HistType::kTH1F, {nlmAxis});
    mHistManager.add("hClusterEtaPhi", "hClusterEtaPhi", O2HistType::kTH2F, {etaAxis, phiAxis});
    mHistManager.add("hClusterTime", "hClusterTime", O2HistType::kTH1F, {timeAxis});
    mHistManager.add("hGlobalTrackEtaPhi", "hGlobalTrackEtaPhi", O2HistType::kTH2F, {etaAxis, phiAxis});
    mHistManager.add("hGlobalTrackMult", "hGlobalTrackMult", O2HistType::kTH1D, {{200, -0.5, 199.5, "N_{trk}"}});
    mHistManager.add("hCollisionType", "hCollisionType;;#it{count}", O2HistType::kTH1D, {{3, -0.5, 2.5}});
    auto hCollisionType = mHistManager.get<TH1>(HIST("hCollisionType"));
    hCollisionType->GetXaxis()->SetBinLabel(1, "no collision");
    hCollisionType->GetXaxis()->SetBinLabel(2, "normal collision");
    hCollisionType->GetXaxis()->SetBinLabel(3, "mult. collisions");
    mHistManager.add("hBCMatchErrors", "hBCMatchErrors;;#it{N}_{BC}", O2HistType::kTH1D, {{3, -0.5, 2.5}});
    auto hBCMatchErrors = mHistManager.get<TH1>(HIST("hBCMatchErrors"));
    hBCMatchErrors->GetXaxis()->SetBinLabel(1, "Normal");
    hBCMatchErrors->GetXaxis()->SetBinLabel(2, "Wrong collisionID order");
    hBCMatchErrors->GetXaxis()->SetBinLabel(3, "foundBCId != globalIndex");
    mHistManager.add("hClusterType", "hClusterType;;#it{count}", O2HistType::kTH1D, {{3, -0.5, 2.5}});
    auto hClusterType = mHistManager.get<TH1>(HIST("hClusterType"));
    hClusterType->GetXaxis()->SetBinLabel(1, "no collision");
    hClusterType->GetXaxis()->SetBinLabel(2, "normal collision");
    hClusterType->GetXaxis()->SetBinLabel(3, "mult. collisions");
    mHistManager.add("hCollPerBC", "hCollPerBC;#it{N}_{coll.};#it{count}", O2HistType::kTH1D, {{100, -0.5, 99.5}});
    mHistManager.add("hBC", "hBC;;#it{count}", O2HistType::kTH1D, {{8, -0.5, 7.5}});
    mHistManager.add("hCollisionTimeReso", "hCollisionTimeReso;#Delta t_{coll};#it{count}", O2HistType::kTH1D, {{2000, 0, 2000}});
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
      mHistManager.add("hContributors", "hContributors;contributor per cell hit;#it{counts}", O2HistType::kTH1I, {{20, 0, 20}});
      mHistManager.add("hMCParticleEnergy", "hMCParticleEnergy;#it{E} (GeV/#it{c});#it{counts}", O2HistType::kTH1F, {energyAxis});
    }
  }

  // void process(aod::Collision const& collision, soa::Filtered<aod::Tracks> const& fullTracks, aod::Calos const& cells)
  // void process(aod::Collision const& collision, aod::Tracks const& tracks, aod::Calos const& cells)
  // void process(aod::BCs const& bcs, aod::Collision const& collision, aod::Calos const& cells)

  //  Appears to need the BC to be accessed to be available in the collision table...
  void processFull(BcEvSels const& bcs, CollEventSels const& collisions, MyGlobTracks const& tracks, FilteredCells const& cells)
  {
    LOG(debug) << "Starting process full.";

    int previousCollisionId = 0; // Collision ID of the last unique BC. Needed to skip unordered collisions to ensure ordered collisionIds in the cluster table
    int nBCsProcessed = 0;
    int nCellsProcessed = 0;
    std::unordered_map<uint64_t, int> numberCollsInBC; // Number of collisions mapped to the global BC index of all BCs
    std::unordered_map<uint64_t, int> numberCellsInBC; // Number of cells mapped to the global BC index of all BCs to check whether EMCal was readout
    for (const auto& bc : bcs) {
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
      for (const auto& cell : cellsInBC) {
        auto amplitude = cell.amplitude();
        if (static_cast<bool>(hasShaperCorrection) && emcal::intToChannelType(cell.cellType()) == emcal::ChannelType_t::LOW_GAIN) { // Apply shaper correction to LG cells
          amplitude = o2::emcal::NonlinearityHandler::evaluateShaperCorrectionCellEnergy(amplitude);
        }
        if (applyCellAbsScale) {
          amplitude *= getAbsCellScale(cell.cellNumber());
        }
        cellsBC.emplace_back(cell.cellNumber(),
                             amplitude,
                             cell.time() + getCellTimeShift(cell.cellNumber(), amplitude, o2::emcal::intToChannelType(cell.cellType())),
                             o2::emcal::intToChannelType(cell.cellType()));
        cellIndicesBC.emplace_back(cell.globalIndex());
      }
      LOG(detail) << "Number of cells for BC (CF): " << cellsBC.size();
      nCellsProcessed += cellsBC.size();

      fillQAHistogram(cellsBC);

      // TODO: Helpful for now, but should be removed.
      LOG(debug) << "Converted EMCAL cells";
      for (const auto& cell : cellsBC) {
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
            if (previousCollisionId > col.globalIndex()) {
              mHistManager.fill(HIST("hBCMatchErrors"), 1);
              continue;
            }
            previousCollisionId = col.globalIndex();
            if (col.foundBCId() == bc.globalIndex()) {
              mHistManager.fill(HIST("hBCMatchErrors"), 0); // CollisionID ordered and foundBC matches -> Fill as healthy
              mHistManager.fill(HIST("hCollisionTimeReso"), col.collisionTimeRes());
              mHistManager.fill(HIST("hCollPerBC"), 1);
              mHistManager.fill(HIST("hCollisionType"), 1);
              math_utils::Point3D<float> vertexPos = {col.posX(), col.posY(), col.posZ()};

              std::vector<std::vector<int>> clusterToTrackIndexMap;
              std::vector<std::vector<int>> trackToClusterIndexMap;
              std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>> indexMapPair{clusterToTrackIndexMap, trackToClusterIndexMap};
              std::vector<int64_t> trackGlobalIndex;
              doTrackMatching<CollEventSels::filtered_iterator>(col, tracks, indexMapPair, vertexPos, trackGlobalIndex);

              // Store the clusters in the table where a matching collision could
              // be identified.
              fillClusterTable<CollEventSels::filtered_iterator>(col, vertexPos, iClusterizer, cellIndicesBC, indexMapPair, trackGlobalIndex);
            } else {
              mHistManager.fill(HIST("hBCMatchErrors"), 2);
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
          fillAmbigousClusterTable<BcEvSels::iterator>(bc, iClusterizer, cellIndicesBC, hasCollision);
        }

        LOG(debug) << "Cluster loop done for clusterizer " << iClusterizer;
      } // end of clusterizer loop
      LOG(debug) << "Done with process BC.";
      nBCsProcessed++;
    } // end of bc loop

    // Loop through all collisions and fill emcalcollisionmatch with a boolean stating, whether the collision was ambiguous (not the only collision in its BC)
    for (const auto& collision : collisions) {
      auto globalbcid = collision.foundBC_as<BcEvSels>().globalIndex();
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

  void processMCFull(BcEvSels const& bcs, CollEventSels const& collisions, MyGlobTracks const& tracks, FilteredMcCells const& cells, aod::StoredMcParticles_001 const&)
  {
    LOG(debug) << "Starting process full.";

    int previousCollisionId = 0; // Collision ID of the last unique BC. Needed to skip unordered collisions to ensure ordered collisionIds in the cluster table
    int nBCsProcessed = 0;
    int nCellsProcessed = 0;
    std::unordered_map<uint64_t, int> numberCollsInBC; // Number of collisions mapped to the global BC index of all BCs
    std::unordered_map<uint64_t, int> numberCellsInBC; // Number of cells mapped to the global BC index of all BCs to check whether EMCal was readout
    for (const auto& bc : bcs) {
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
      for (const auto& cell : cellsInBC) {
        mHistManager.fill(HIST("hContributors"), cell.mcParticle_as<aod::StoredMcParticles_001>().size());
        auto cellParticles = cell.mcParticle_as<aod::StoredMcParticles_001>();
        for (const auto& cellparticle : cellParticles) {
          mHistManager.fill(HIST("hMCParticleEnergy"), cellparticle.e());
        }
        auto amplitude = cell.amplitude();
        if (static_cast<bool>(hasShaperCorrection) && emcal::intToChannelType(cell.cellType()) == emcal::ChannelType_t::LOW_GAIN) { // Apply shaper correction to LG cells
          amplitude = o2::emcal::NonlinearityHandler::evaluateShaperCorrectionCellEnergy(amplitude);
        }
        cellsBC.emplace_back(cell.cellNumber(),
                             amplitude,
                             cell.time() + getCellTimeShift(cell.cellNumber(), amplitude, o2::emcal::intToChannelType(cell.cellType())),
                             o2::emcal::intToChannelType(cell.cellType()));
        cellIndicesBC.emplace_back(cell.globalIndex());
        cellLabels.emplace_back(cell.mcParticleIds(), cell.amplitudeA());
      }
      LOG(detail) << "Number of cells for BC (CF): " << cellsBC.size();
      nCellsProcessed += cellsBC.size();

      fillQAHistogram(cellsBC);

      // TODO: Helpful for now, but should be removed.
      LOG(debug) << "Converted EMCAL cells";
      for (const auto& cell : cellsBC) {
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
            if (previousCollisionId > col.globalIndex()) {
              mHistManager.fill(HIST("hBCMatchErrors"), 1);
              continue;
            }
            previousCollisionId = col.globalIndex();
            if (col.foundBCId() == bc.globalIndex()) {
              mHistManager.fill(HIST("hBCMatchErrors"), 0); // CollisionID ordered and foundBC matches -> Fill as healthy
              mHistManager.fill(HIST("hCollPerBC"), 1);
              mHistManager.fill(HIST("hCollisionType"), 1);
              math_utils::Point3D<float> vertexPos = {col.posX(), col.posY(), col.posZ()};

              std::vector<std::vector<int>> clusterToTrackIndexMap;
              std::vector<std::vector<int>> trackToClusterIndexMap;
              std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>> indexMapPair{clusterToTrackIndexMap, trackToClusterIndexMap};
              std::vector<int64_t> trackGlobalIndex;
              doTrackMatching<CollEventSels::filtered_iterator>(col, tracks, indexMapPair, vertexPos, trackGlobalIndex);

              // Store the clusters in the table where a matching collision could
              // be identified.
              fillClusterTable<CollEventSels::filtered_iterator>(col, vertexPos, iClusterizer, cellIndicesBC, indexMapPair, trackGlobalIndex);
            } else {
              mHistManager.fill(HIST("hBCMatchErrors"), 2);
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
          fillAmbigousClusterTable<BcEvSels::iterator>(bc, iClusterizer, cellIndicesBC, hasCollision);
        }
        LOG(debug) << "Cluster loop done for clusterizer " << iClusterizer;
      } // end of clusterizer loop
      LOG(debug) << "Done with process BC.";
      nBCsProcessed++;
    } // end of bc loop

    // Loop through all collisions and fill emcalcollisionmatch with a boolean stating, whether the collision was ambiguous (not the only collision in its BC)
    for (const auto& collision : collisions) {
      auto globalbcid = collision.foundBC_as<BcEvSels>().globalIndex();
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
  void processStandalone(aod::BCs const& bcs, aod::Collisions const& collisions, FilteredCells const& cells)
  {
    LOG(debug) << "Starting process standalone.";
    int previousCollisionId = 0; // Collision ID of the last unique BC. Needed to skip unordered collisions to ensure ordered collisionIds in the cluster table
    int nBCsProcessed = 0;
    int nCellsProcessed = 0;
    for (const auto& bc : bcs) {
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
      for (const auto& cell : cellsInBC) {
        cellsBC.emplace_back(cell.cellNumber(),
                             cell.amplitude(),
                             cell.time() + getCellTimeShift(cell.cellNumber(), cell.amplitude(), o2::emcal::intToChannelType(cell.cellType())),
                             o2::emcal::intToChannelType(cell.cellType()));
        cellIndicesBC.emplace_back(cell.globalIndex());
      }
      LOG(detail) << "Number of cells for BC (CF): " << cellsBC.size();
      nCellsProcessed += cellsBC.size();

      fillQAHistogram(cellsBC);

      // TODO: Helpful for now, but should be removed.
      LOG(debug) << "Converted EMCAL cells";
      for (const auto& cell : cellsBC) {
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
            if (previousCollisionId > col.globalIndex()) {
              mHistManager.fill(HIST("hBCMatchErrors"), 1);
              continue;
            }
            previousCollisionId = col.globalIndex();
            mHistManager.fill(HIST("hBCMatchErrors"), 0); // CollisionID ordered and foundBC matches -> Fill as healthy
            mHistManager.fill(HIST("hCollPerBC"), 1);
            mHistManager.fill(HIST("hCollisionType"), 1);
            math_utils::Point3D<float> vertexPos = {col.posX(), col.posY(), col.posZ()};

            // Store the clusters in the table where a matching collision could
            // be identified.
            fillClusterTable<aod::Collision>(col, vertexPos, iClusterizer, cellIndicesBC);
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
          fillAmbigousClusterTable<aod::BC>(bc, iClusterizer, cellIndicesBC, hasCollision);
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
    // in preparation for future O2 changes
    // mClusterFactories.setClusterizerSettings(mClusterDefinitions.at(iClusterizer).minCellEnergy, mClusterDefinitions.at(iClusterizer).timeMin, mClusterDefinitions.at(iClusterizer).timeMax, mClusterDefinitions.at(iClusterizer).recalcShowerShape5x5);
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
  void fillClusterTable(Collision const& col, math_utils::Point3D<float> const& vertexPos, size_t iClusterizer, const gsl::span<int64_t> cellIndicesBC, std::optional<std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>>> const& indexMapPair = std::nullopt, std::optional<std::vector<int64_t>> const& trackGlobalIndex = std::nullopt)
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
      pos = pos - vertexPos;
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
      mHistManager.fill(HIST("hClusterNLM"), cluster.getNExMax());
      mHistManager.fill(HIST("hClusterTime"), cluster.getClusterTime());
      mHistManager.fill(HIST("hClusterEtaPhi"), pos.Eta(), TVector2::Phi_0_2pi(pos.Phi()));
      if (indexMapPair && trackGlobalIndex) {
        for (unsigned int iTrack = 0; iTrack < std::get<0>(*indexMapPair)[iCluster].size(); iTrack++) {
          if (std::get<0>(*indexMapPair)[iCluster][iTrack] >= 0) {
            LOG(debug) << "Found track " << (*trackGlobalIndex)[std::get<0>(*indexMapPair)[iCluster][iTrack]] << " in cluster " << cluster.getID();
            matchedTracks(clusters.lastIndex(), (*trackGlobalIndex)[std::get<0>(*indexMapPair)[iCluster][iTrack]]);
          }
        }
      }
      iCluster++;
    } // end of cluster loop
  }

  template <typename BC>
  void fillAmbigousClusterTable(BC const& bc, size_t iClusterizer, const gsl::span<int64_t> cellIndicesBC, bool hasCollision)
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
  void doTrackMatching(Collision const& col, MyGlobTracks const& tracks, std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>>& indexMapPair, math_utils::Point3D<float>& vertexPos, std::vector<int64_t>& trackGlobalIndex)
  {
    auto groupedTracks = tracks.sliceBy(perCollision, col.globalIndex());
    int nTracksInCol = groupedTracks.size();
    std::vector<double> trackPhi;
    std::vector<double> trackEta;
    // reserve memory to reduce on the fly memory allocation
    trackPhi.reserve(nTracksInCol);
    trackEta.reserve(nTracksInCol);
    trackGlobalIndex.reserve(nTracksInCol);
    fillTrackInfo<decltype(groupedTracks)>(groupedTracks, trackPhi, trackEta, trackGlobalIndex);

    int nClusterInCol = mAnalysisClusters.size();
    std::vector<double> clusterPhi;
    std::vector<double> clusterEta;
    clusterPhi.reserve(nClusterInCol);
    clusterEta.reserve(nClusterInCol);

    // TODO one loop that could in principle be combined with the other
    // loop to improve performance
    for (const auto& cluster : mAnalysisClusters) {
      // Determine the cluster eta, phi, correcting for the vertex
      // position.
      auto pos = cluster.getGlobalPosition();
      pos = pos - vertexPos;
      // Normalize the vector and rescale by energy.
      pos *= (cluster.E() / std::sqrt(pos.Mag2()));
      clusterPhi.emplace_back(TVector2::Phi_0_2pi(pos.Phi()));
      clusterEta.emplace_back(pos.Eta());
    }
    indexMapPair =
      jetutilities::MatchClustersAndTracks(clusterPhi, clusterEta,
                                           trackPhi, trackEta,
                                           maxMatchingDistance, 20);
  }

  template <typename Tracks>
  void fillTrackInfo(Tracks const& tracks, std::vector<double>& trackPhi, std::vector<double>& trackEta, std::vector<int64_t>& trackGlobalIndex)
  {
    int nTrack = 0;
    for (const auto& track : tracks) {
      // TODO only consider tracks in current emcal/dcal acceptanc
      if (!track.isGlobalTrack()) { // only global tracks
        continue;
      }
      // Tracks that do not point to the EMCal/DCal/PHOS get default values of -999
      // This way we can cut out tracks that do not point to the EMCal+DCal
      if (track.trackEtaEmcal() < -900 || track.trackPhiEmcal() < -900) {
        continue;
      }
      if (trackMinPt > 0 && track.pt() < trackMinPt) {
        continue;
      }
      nTrack++;
      trackPhi.emplace_back(TVector2::Phi_0_2pi(track.trackPhiEmcal()));
      trackEta.emplace_back(track.trackEtaEmcal());
      mHistManager.fill(HIST("hGlobalTrackEtaPhi"), track.trackEtaEmcal(),
                        TVector2::Phi_0_2pi(track.trackPhiEmcal()));
      trackGlobalIndex.emplace_back(track.globalIndex());
    }
    mHistManager.fill(HIST("hGlobalTrackMult"), nTrack);
  }

  void countBC(int numberOfCollisions, bool hasEMCCells)
  {
    int emcDataOffset = hasEMCCells ? 0 : 3;
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
    if (hasEMCCells) {
      mHistManager.fill(HIST("hBC"), 0);
    }
    mHistManager.fill(HIST("hBC"), 1 + emcDataOffset + collisionOffset);
  }

  void fillQAHistogram(const gsl::span<o2::emcal::Cell> cellsBC)
  {
    // Cell QA
    // For convenience, use the clusterizer stored geometry to get the eta-phi
    for (const auto& cell : cellsBC) {
      mHistManager.fill(HIST("hCellE"), cell.getEnergy());
      if (cell.getLowGain())
        mHistManager.fill(HIST("hLGCellTimeEnergy"), cell.getTimeStamp(), cell.getEnergy());
      else if (cell.getHighGain())
        mHistManager.fill(HIST("hHGCellTimeEnergy"), cell.getTimeStamp(), cell.getEnergy());
      mHistManager.fill(HIST("hCellTowerID"), cell.getTower());
      auto res = mClusterizers.at(0)->getGeometry()->EtaPhiFromIndex(cell.getTower());
      mHistManager.fill(HIST("hCellEtaPhi"), std::get<0>(res), TVector2::Phi_0_2pi(std::get<1>(res)));
      res = mClusterizers.at(0)->getGeometry()->GlobalRowColFromIndex(cell.getTower());
      // NOTE: Reversed column and row because it's more natural for presentation.
      mHistManager.fill(HIST("hCellRowCol"), std::get<1>(res), std::get<0>(res));
    }
  }

  float getAbsCellScale(const int cellID)
  {
    // Apply cell scale based on SM types (Full, Half (not used), EMC 1/3, DCal, DCal 1/3)
    // Same as in Run2 data
    if (applyCellAbsScale == 1) {
      int iSM = mClusterizers.at(0)->getGeometry()->GetSuperModuleNumber(cellID);
      return cellAbsScaleFactors.value[mClusterizers.at(0)->getGeometry()->GetSMType(iSM)];

      // Apply cell scale based on columns to accoutn for material of TRD structures
    } else if (applyCellAbsScale == 2) {
      auto res = mClusterizers.at(0)->getGeometry()->GlobalRowColFromIndex(cellID);
      return cellAbsScaleFactors.value[std::get<1>(res)];
    } else {
      return 1.f;
    }
  }

  // Apply shift of the cell time in data and MC
  // In MC this has to be done to shift the cell time, which is not calibrated to 0 due to the flight time of the particles to the EMCal surface (~15ns)
  // In data this is done to correct for the time walk effect
  float getCellTimeShift(const int16_t cellID, const float cellEnergy, const emcal::ChannelType_t cellType)
  {
    if (!applyCellTimeCorrection) {
      return 0.f;
    }
    float timeshift = 0.f;
    float timesmear = 0.f;
    if (isMC) { // ---> MC
      // Shift the time to 0, as the TOF was simulated -> eta dependent shift (as larger eta values are further away from collision point)
      // Use distance between vertex and EMCal (at eta = 0) and distance on EMCal surface (cell size times column) to calculate distance to cell
      // 0.2 is cell size in m (0.06) divided by the speed of light in m/ns (0.3) - 47.5 is the "middle" of the EMCal (2*48 cells in one column)
      float timeCol = 0.2f * (geometry->GlobalCol(cellID) - 47.5f); // calculate time to get to specific column
      timeshift = -std::sqrt(215.f + timeCol * timeCol);            // 215 is 14.67ns^2 (time it takes to get the cell at eta = 0)

      // Also smear the time to account for the broader time resolution in data than in MC
      if (cellEnergy < 0.3)                                                       // Cells with tless than 300 MeV cannot be the leading cell in the cluster, so their time does not require precise calibration
        timesmear = 0.;                                                           // They will therefore not be smeared and only get their shift
      else if (cellType == emcal::ChannelType_t::HIGH_GAIN)                       // High gain cells -> Low energies
        timesmear = normalgaus(rdgen) * (1.6 + 9.5 * std::exp(-3. * cellEnergy)); // Parameters extracted from LHC24f3b & LHC22o (pp), but also usable for other periods
      else if (cellType == emcal::ChannelType_t::LOW_GAIN)                        // Low gain cells -> High energies
        timesmear = normalgaus(rdgen) * (5.0);                                    // Parameters extracted from LHC24g4 & LHC24aj (pp), but also usable for other periods

    } else {                                                    // ---> Data
      if (cellEnergy < 0.3) {                                   // Cells with tless than 300 MeV cannot be the leading cell in the cluster, so their time does not require precise calibration
        timeshift = 0.;                                         // In data they will not be shifted (they are close to 0 anyways)
      } else if (cellType == emcal::ChannelType_t::HIGH_GAIN) { // High gain cells -> Low energies
        if (cellEnergy < 4.)                                    // Low energy regime
          timeshift = 0.8 * std::log(2.7 * cellEnergy);         // Parameters extracted from LHC22o (pp), but also usable for other periods
        else                                                    // Medium energy regime
          timeshift = 1.5 * std::log(0.9 * cellEnergy);         // Parameters extracted from LHC22o (pp), but also usable for other periods
      } else if (cellType == emcal::ChannelType_t::LOW_GAIN) {  // Low gain cells -> High energies
        if (cellEnergy < 30.)                                   // High energy regime
          timeshift = 1.9 * std::log(0.09 * cellEnergy);        // Parameters extracted from LHC24aj (pp), but also usable for other periods
        else                                                    // Very high energy regime
          timeshift = 1.9;                                      // Parameters extracted from LHC24aj (pp), but also usable for other periods
      }
      LOG(debug) << "Shift the cell time by " << timeshift << " + " << timesmear << " ns";
    }
    return timeshift + timesmear;
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<EmcalCorrectionTask>(cfgc)};
}
