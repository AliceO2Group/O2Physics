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
/// \file emcalCorrectionTask.cxx
/// \brief Task that provides EMCal clusters and applies necessary corrections
/// \author Raymond Ehlers (raymond.ehlers@cern.ch) ORNL, Florian Jonas (florian.jonas@cern.ch), Marvin Hemmer (marvin.hemmer@cern.ch)

#include "PWGJE/Core/emcalCrossTalkEmulation.h"
#include "PWGJE/Core/utilsTrackMatchingEMC.h"
#include "PWGJE/DataModel/EMCALClusterDefinition.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"
//
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h" // for EM V0 legs

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsEMCAL/AnalysisCluster.h>
#include <DataFormatsEMCAL/Cell.h>
#include <DataFormatsEMCAL/CellLabel.h>
#include <DataFormatsEMCAL/ClusterLabel.h>
#include <DataFormatsEMCAL/Constants.h>
#include <DetectorsBase/GeometryManager.h>
#include <EMCALBase/ClusterFactory.h>
#include <EMCALBase/Geometry.h>
#include <EMCALBase/NonlinearityHandler.h>
#include <EMCALCalib/GainCalibrationFactors.h>
#include <EMCALCalibration/EMCALTempCalibExtractor.h>
#include <EMCALReconstruction/Clusterizer.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/WorkflowSpec.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <fairlogger/Logger.h>

#include <GPUROOTCartesianFwd.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <gsl/span>
#include <memory>
#include <random>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::emccrosstalk;
using namespace tmemcutilities;
using MyGlobTracks = o2::soa::Join<o2::aod::FullTracks, o2::aod::TrackSelection>;
using BcEvSels = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
using CollEventSels = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
using FilteredCells = o2::soa::Filtered<aod::Calos>;
using McCells = o2::soa::Join<aod::Calos, aod::McCaloLabels_001>;
using FilteredMcCells = o2::soa::Filtered<McCells>;
using EMV0Legs = aod::V0Legs;

enum CellScaleMode {
  ModeNone = 0,
  ModeSMWise = 1,
  ModeColumnWise = 2,
  NumberModes = 3
};

struct EmcalCorrectionTask {
  Produces<o2::aod::EMCALClusters> clusters;
  Produces<o2::aod::EMCALMCClusters> mcclusters;
  Produces<o2::aod::EMCALAmbiguousClusters> clustersAmbiguous;
  Produces<o2::aod::EMCALAmbiguousMCClusters> mcclustersAmbiguous;
  Produces<o2::aod::EMCALClusterCells> clustercells; // cells belonging to given cluster
  Produces<o2::aod::EMCALAmbiguousClusterCells> clustercellsambiguous;
  Produces<o2::aod::EMCALMatchedTracks> matchedTracks;
  Produces<o2::aod::EMCMatchSecs> matchedSecondaries;
  Produces<o2::aod::EMCALMatchedCollisions> emcalcollisionmatch;

  // Preslices
  Preslice<MyGlobTracks> perCollision = o2::aod::track::collisionId;
  PresliceUnsorted<CollEventSels> collisionsPerFoundBC = aod::evsel::foundBCId;
  Preslice<aod::Collisions> collisionsPerBC = aod::collision::bcId;
  Preslice<FilteredCells> cellsPerFoundBC = aod::calo::bcId;
  Preslice<FilteredMcCells> mcCellsPerFoundBC = aod::calo::bcId;
  PresliceUnsorted<EMV0Legs> perCollisionEMV0Legs = aod::v0leg::collisionId;

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
  Configurable<bool> fillQA{"fillQA", false, "Switch to turn on QA histograms."};
  Configurable<bool> useCCDBAlignment{"useCCDBAlignment", false, "EXPERTS ONLY! Switch to use the alignment object stored in CCDB instead of using the default alignment from the global geometry object."};
  Configurable<bool> applyTempCalib{"applyTempCalib", false, "Switch to turn on Temperature calibration."};
  Configurable<std::string> pathTempCalibCCDB{"pathTempCalibCCDB", "Users/j/jokonig/EMCalTempCalibParams", "Path in the ccdb where slope and intercept for each cell are stored"}; // change to official path as soon as it is available
  Configurable<bool> useTempCalibMean{"useTempCalibMean", false, "Switch to turn on Temperature mean calculation instead of median."};
  Configurable<float> mcCellEnergyShift{"mcCellEnergyShift", 1., "Relative shift of the MC cell energy. 1.1 for 10% shift to higher mass, etc. Only applied to MC."};
  Configurable<float> mcCellEnergyResolutionBroadening{"mcCellEnergyResolutionBroadening", 0., "Relative widening of the MC cell energy resolution. 0 for no widening, 0.1 for 10% widening, etc. Only applied to MC."};
  Configurable<bool> applyGainCalibShift{"applyGainCalibShift", false, "Apply shift for cell gain calibration to use values before cell format change (Sept. 2023)"};

  // cross talk emulation configs
  EmcCrossTalkConf emcCrossTalkConf;

  // cross talk emulation class for handling the cross talk
  emccrosstalk::EMCCrossTalk emcCrossTalk;

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

  // Cluster Eta and Phi used for track matching later
  std::vector<float> mClusterPhi;
  std::vector<float> mClusterEta;

  std::vector<o2::aod::EMCALClusterDefinition> mClusterDefinitions;
  // QA
  o2::framework::HistogramRegistry mHistManager{"EMCALCorrectionTaskQAHistograms"};

  // Random number generator to draw cell time smearing for MC
  std::random_device rd{};
  std::mt19937_64 rdgen{rd()};
  std::normal_distribution<> normalgaus{0, 1}; // mean = 0, stddev = 1 (apply amplitude of smearing after drawing random for performance reasons)

  // EMCal geometry
  o2::emcal::Geometry* geometry;

  // EMCal cell temperature calibrator
  std::unique_ptr<o2::emcal::EMCALTempCalibExtractor> mTempCalibExtractor;
  bool mIsTempCalibInitialized = false;

  // Gain calibration
  std::array<float, 17664> mArrGainCalibDiff;

  std::vector<std::pair<int, int>> mExtraTimeShiftRunRanges;

  // Current run number
  int runNumber{0};

  static constexpr float TrackNotOnEMCal = -900.f;
  static constexpr int kMaxMatchesPerCluster = 20; // Maximum number of tracks to match per cluster

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
    if (useCCDBAlignment.value) {
      geometry->SetMisalMatrixFromCcdb();
    }

    if (applyTempCalib) {
      mTempCalibExtractor = std::make_unique<o2::emcal::EMCALTempCalibExtractor>();
    }

    // gain calibration shift initialization
    if (applyGainCalibShift) {
      initializeGainCalibShift();
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

    // 500 clusters per event is a good upper limit
    mClusterPhi.reserve(500 * mClusterizers.size());
    mClusterEta.reserve(500 * mClusterizers.size());

    mNonlinearityHandler = o2::emcal::NonlinearityFactory::getInstance().getNonlinearity(static_cast<std::string>(nonlinearityFunction));
    LOG(info) << "Using nonlinearity parameterisation: " << nonlinearityFunction.value;
    LOG(info) << "Apply shaper saturation correction:  " << (hasShaperCorrection.value ? "yes" : "no");

    LOG(debug) << "Completed init!";

    // Define the cell energy binning
    std::vector<double> cellEnergyBins;
    for (int i = 0; i < 51; i++) {                      // o2-linter: disable=magic-number (just numbers for binning)
      cellEnergyBins.emplace_back(0.1 * (i - 0) + 0.0); // from 0 to 5 GeV/c, every 0.1 GeV
    }
    for (int i = 51; i < 76; i++) {                      // o2-linter: disable=magic-number (just numbers for binning)
      cellEnergyBins.emplace_back(0.2 * (i - 51) + 5.2); // from 5.2 to 10.0 GeV, every 0.2 GeV
    }
    for (int i = 76; i < 166; i++) {                    // o2-linter: disable=magic-number (just numbers for binning)
      cellEnergyBins.emplace_back(1. * (i - 76) + 11.); // from 11.0 to 100. GeV, every 1 GeV
    }

    // Setup QA hists.
    // NOTE: This is not comprehensive.
    using O2HistType = o2::framework::HistType;
    o2::framework::AxisSpec energyAxis{200, 0., 100., "#it{E} (GeV)"},
      timeAxis{300, -100, 200., "#it{t} (ns)"},
      etaAxis{160, -0.8, 0.8, "#it{#eta}"},
      phiAxis{72, 0, 2 * 3.14159, "#it{#varphi} (rad)"},
      nlmAxis{50, -0.5, 49.5, "NLM"},
      fCrossAxis{100, 0., 1., "F_{+}"},
      sigmaLongAxis{100, 0., 1.0, "#sigma^{2}_{long}"},
      sigmaShortAxis{100, 0., 1.0, "#sigma^{2}_{short}"},
      nCellAxis{60, -0.5, 59.5, "#it{n}_{cells}"},
      energyDenseAxis = {7000, 0.f, 70.f, "#it{E}_{cell} (GeV)"};
    o2::framework::AxisSpec axisDeltaEta{400, -0.2, 0.2, "#Delta#eta"};
    o2::framework::AxisSpec axisDeltaPhi{400, -0.2, 0.2, "#Delta#varphi (rad)"};
    o2::framework::AxisSpec axisNCluster{1000, 0, 1000, "#it{N}_{clus.}"};
    mHistManager.add("hCellE", "hCellE", O2HistType::kTH1D, {energyAxis});
    mHistManager.add("hCellTowerID", "hCellTowerID", O2HistType::kTH1D, {{20000, 0, 20000}});
    mHistManager.add("hCellEtaPhi", "hCellEtaPhi", O2HistType::kTH2F, {etaAxis, phiAxis});
    mHistManager.add("hHGCellTimeEnergy", "hCellTime", O2HistType::kTH2F, {{300, -30, 30}, cellEnergyBins}); // Cell time vs energy for high gain cells (low energies)
    mHistManager.add("hLGCellTimeEnergy", "hCellTime", O2HistType::kTH2F, {{300, -30, 30}, cellEnergyBins}); // Cell time vs energy for low gain cells (high energies)
    mHistManager.add("hTempCalibCorrection", "hTempCalibCorrection", O2HistType::kTH1F, {{5000, 0.5, 1.5}});
    // NOTE: Reversed column and row because it's more natural for presentation.
    mHistManager.add("hCellRowCol", "hCellRowCol;Column;Row", O2HistType::kTH2D, {{96, -0.5, 95.5}, {208, -0.5, 207.5}});
    mHistManager.add("hClusterE", "hClusterE", O2HistType::kTH1D, {energyAxis});
    mHistManager.add("hNCluster", "hNCluster", O2HistType::kTH1D, {axisNCluster});
    mHistManager.add("hClusterNLM", "hClusterNLM", O2HistType::kTH1D, {nlmAxis});
    mHistManager.add("hClusterEtaPhi", "hClusterEtaPhi", O2HistType::kTH2F, {etaAxis, phiAxis});
    mHistManager.add("hClusterTime", "hClusterTime", O2HistType::kTH1D, {timeAxis});
    mHistManager.add("hCollisionType", "hCollisionType;;#it{count}", O2HistType::kTH1D, {{3, -0.5, 2.5}});
    mHistManager.add("hMatchedPrimaryTracks", "hMatchedPrimaryTracks", O2HistType::kTH2F, {axisDeltaEta, axisDeltaPhi});
    mHistManager.add("hMatchedSecondaries", "hMatchedSecondaries", O2HistType::kTH2F, {axisDeltaEta, axisDeltaPhi});
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
    if (isMC.value) {
      mHistManager.add("hContributors", "hContributors;contributor per cell hit;#it{counts}", O2HistType::kTH1D, {{20, 0, 20}});
      mHistManager.add("hMCParticleEnergy", "hMCParticleEnergy;#it{E} (GeV/#it{c});#it{counts}", O2HistType::kTH1D, {energyAxis});
    }
    if (fillQA.value) {
      mHistManager.add("hClusterNCellE", "hClusterNCellE", O2HistType::kTH2D, {energyAxis, nCellAxis});
      mHistManager.add("hClusterFCrossE", "hClusterFCrossE", O2HistType::kTH2D, {energyAxis, fCrossAxis});
      mHistManager.add("hClusterFCrossSigmaLongE", "hClusterFCrossSigmaLongE", O2HistType::kTH3F, {energyAxis, fCrossAxis, sigmaLongAxis});
      mHistManager.add("hClusterFCrossSigmaShortE", "hClusterFCrossSigmaShortE", O2HistType::kTH3F, {energyAxis, fCrossAxis, sigmaShortAxis});
    }

    if (isMC.value && emcCrossTalkConf.enableCrossTalk.value) {
      emcCrossTalk.initObjects(emcCrossTalkConf);
      if (emcCrossTalkConf.createHistograms.value) {
        mHistManager.add<TH1>("hCellEnergyDistBefore", "Cell energy before cross-talk emulation", {o2::framework::HistType::kTH1D, {energyDenseAxis}});
        mHistManager.add<TH1>("hCellEnergyDistAfter", "Cell energy after cross-talk emulation", {o2::framework::HistType::kTH1D, {energyDenseAxis}});
      }
    }

    // For some runs, LG cells require an extra time shift of 2 * 8.8ns due to problems in the time calibration
    // Affected run ranges (inclusive) are initialised here (min,max)
    mExtraTimeShiftRunRanges.emplace_back(535365, 535645); // LHC23g-LHC23h
    mExtraTimeShiftRunRanges.emplace_back(535725, 536126); // LHC23h-LHC23l
    mExtraTimeShiftRunRanges.emplace_back(536199, 536202); // LHC23l-LHC23m
    mExtraTimeShiftRunRanges.emplace_back(536239, 536346); // LHC23m-LHC23n
    mExtraTimeShiftRunRanges.emplace_back(536565, 536590); // Commisioning-LHC23r
    mExtraTimeShiftRunRanges.emplace_back(542280, 543854); // LHC23zv-LHC23zy
    mExtraTimeShiftRunRanges.emplace_back(559544, 559856); // PbPb 2024
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

      // get run number
      runNumber = bc.runNumber();

      if (applyTempCalib && !mIsTempCalibInitialized) { // needs to be called once
        mTempCalibExtractor->InitializeFromCCDB(pathTempCalibCCDB, static_cast<uint64_t>(runNumber));
        mIsTempCalibInitialized = true;
      }

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
        if (applyGainCalibShift) {
          amplitude *= mArrGainCalibDiff[cell.cellNumber()];
        }
        if (applyTempCalib) {
          float tempCalibFactor = mTempCalibExtractor->getGainCalibFactor(static_cast<uint16_t>(cell.cellNumber()));
          amplitude /= tempCalibFactor;
          mHistManager.fill(HIST("hTempCalibCorrection"), tempCalibFactor);
        }
        cellsBC.emplace_back(cell.cellNumber(),
                             amplitude,
                             cell.time() + getCellTimeShift(cell.cellNumber(), amplitude, o2::emcal::intToChannelType(cell.cellType()), runNumber),
                             o2::emcal::intToChannelType(cell.cellType()));
        cellIndicesBC.emplace_back(cell.globalIndex());
      }
      LOG(detail) << "Number of cells for BC (CF): " << cellsBC.size();
      nCellsProcessed += cellsBC.size();

      fillQAHistogram(cellsBC);

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

              MatchResult indexMapPair;
              std::vector<int64_t> trackGlobalIndex;
              doTrackMatching<CollEventSels::filtered_iterator>(col, tracks, indexMapPair, trackGlobalIndex);

              // Store the clusters in the table where a matching collision could
              // be identified.
              fillClusterTable<CollEventSels::filtered_iterator>(col, vertexPos, iClusterizer, cellIndicesBC, &indexMapPair, &trackGlobalIndex);
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

        mClusterPhi.clear();
        mClusterEta.clear();
        LOG(debug) << "Cluster loop done for clusterizer " << iClusterizer;
      } // end of clusterizer loop
      LOG(debug) << "Done with process BC.";
      nBCsProcessed++;
    } // end of bc loop

    // Loop through all collisions and fill emcalcollisionmatch with a boolean stating, whether the collision was ambiguous (not the only collision in its BC)
    for (const auto& collision : collisions) {
      if (applySoftwareTriggerSelection) {
        if (!zorro.isSelected(collision.foundBC_as<BcEvSels>().globalBC())) {
          continue;
        }
      }
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

  void processWithSecondaries(BcEvSels const& bcs, CollEventSels const& collisions, MyGlobTracks const& tracks, FilteredCells const& cells, EMV0Legs const& v0legs)
  {
    LOG(debug) << "Starting process full.";

    int previousCollisionId = 0; // Collision ID of the last unique BC. Needed to skip unordered collisions to ensure ordered collisionIds in the cluster table
    int nBCsProcessed = 0;
    int nCellsProcessed = 0;
    std::unordered_map<uint64_t, int> numberCollsInBC; // Number of collisions mapped to the global BC index of all BCs
    std::unordered_map<uint64_t, int> numberCellsInBC; // Number of cells mapped to the global BC index of all BCs to check whether EMCal was readout
    for (const auto& bc : bcs) {
      LOG(debug) << "Next BC";

      // get run number
      runNumber = bc.runNumber();

      if (applyTempCalib && !mIsTempCalibInitialized) { // needs to be called once
        mTempCalibExtractor->InitializeFromCCDB(pathTempCalibCCDB, static_cast<uint64_t>(runNumber));
        mIsTempCalibInitialized = true;
      }

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
        if (applyGainCalibShift) {
          amplitude *= mArrGainCalibDiff[cell.cellNumber()];
        }
        if (applyTempCalib) {
          float tempCalibFactor = mTempCalibExtractor->getGainCalibFactor(static_cast<uint16_t>(cell.cellNumber()));
          amplitude /= tempCalibFactor;
          mHistManager.fill(HIST("hTempCalibCorrection"), tempCalibFactor);
        }
        cellsBC.emplace_back(cell.cellNumber(),
                             amplitude,
                             cell.time() + getCellTimeShift(cell.cellNumber(), amplitude, o2::emcal::intToChannelType(cell.cellType()), runNumber),
                             o2::emcal::intToChannelType(cell.cellType()));
        cellIndicesBC.emplace_back(cell.globalIndex());
      }
      LOG(detail) << "Number of cells for BC (CF): " << cellsBC.size();
      nCellsProcessed += cellsBC.size();

      fillQAHistogram(cellsBC);

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

              MatchResult indexMapPair;
              std::vector<int64_t> trackGlobalIndex;
              doTrackMatching<CollEventSels::filtered_iterator>(col, tracks, indexMapPair, trackGlobalIndex);

              MatchResult indexMapPairSecondary;
              std::vector<int64_t> secondaryGlobalIndex;
              doSecondaryTrackMatching<CollEventSels::filtered_iterator>(col, v0legs, indexMapPairSecondary, secondaryGlobalIndex, tracks);

              // Store the clusters in the table where a matching collision could
              // be identified.
              fillClusterTable<CollEventSels::filtered_iterator>(col, vertexPos, iClusterizer, cellIndicesBC, &indexMapPair, &trackGlobalIndex, &indexMapPairSecondary, &secondaryGlobalIndex);
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

        mClusterPhi.clear();
        mClusterEta.clear();
        LOG(debug) << "Cluster loop done for clusterizer " << iClusterizer;
      } // end of clusterizer loop
      LOG(debug) << "Done with process BC.";
      nBCsProcessed++;
    } // end of bc loop

    // Loop through all collisions and fill emcalcollisionmatch with a boolean stating, whether the collision was ambiguous (not the only collision in its BC)
    for (const auto& collision : collisions) {
      if (applySoftwareTriggerSelection) {
        if (!zorro.isSelected(collision.foundBC_as<BcEvSels>().globalBC())) {
          continue;
        }
      }
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
  PROCESS_SWITCH(EmcalCorrectionTask, processWithSecondaries, "run full analysis with secondary track matching", false);

  void processMCFull(BcEvSels const& bcs, CollEventSels const& collisions, MyGlobTracks const& tracks, FilteredMcCells const& cells, aod::StoredMcParticles_001 const&)
  {
    LOG(debug) << "Starting processMCFull.";

    int previousCollisionId = 0; // Collision ID of the last unique BC. Needed to skip unordered collisions to ensure ordered collisionIds in the cluster table
    int nBCsProcessed = 0;
    int nCellsProcessed = 0;
    std::unordered_map<uint64_t, int> numberCollsInBC; // Number of collisions mapped to the global BC index of all BCs
    std::unordered_map<uint64_t, int> numberCellsInBC; // Number of cells mapped to the global BC index of all BCs to check whether EMCal was readout
    for (const auto& bc : bcs) {
      LOG(debug) << "Next BC";
      // Convert aod::Calo to o2::emcal::Cell which can be used with the clusterizer.
      // In particular, we need to filter only EMCAL cells.

      // get run number
      runNumber = bc.runNumber();

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
          const auto& ids = cell.mcParticleIds();
          const auto& amps = cell.amplitudeA();

          if (ids.empty() || amps.empty()) {
            LOGF(warning, "Skipping cell with empty MC info: absId=%d", cell.cellNumber());
            continue;
          }
          mHistManager.fill(HIST("hMCParticleEnergy"), cellparticle.e());
        }
        auto amplitude = cell.amplitude();
        if (static_cast<bool>(hasShaperCorrection) && emcal::intToChannelType(cell.cellType()) == emcal::ChannelType_t::LOW_GAIN) { // Apply shaper correction to LG cells
          amplitude = o2::emcal::NonlinearityHandler::evaluateShaperCorrectionCellEnergy(amplitude);
        }
        if (mcCellEnergyShift != 1.) {
          amplitude *= mcCellEnergyShift; // Fine tune the MC cell energy
        }
        if (mcCellEnergyResolutionBroadening != 0.) {
          amplitude *= (1. + normalgaus(rdgen) * mcCellEnergyResolutionBroadening); // Fine tune the MC cell energy resolution
        }
        cellsBC.emplace_back(cell.cellNumber(),
                             amplitude,
                             cell.time() + getCellTimeShift(cell.cellNumber(), amplitude, o2::emcal::intToChannelType(cell.cellType()), runNumber),
                             o2::emcal::intToChannelType(cell.cellType()));
        cellIndicesBC.emplace_back(cell.globalIndex());
        cellLabels.emplace_back(std::vector<int>{cell.mcParticleIds().begin(), cell.mcParticleIds().end()}, std::vector<float>{cell.amplitudeA().begin(), cell.amplitudeA().end()});
      }
      if (isMC.value && emcCrossTalkConf.enableCrossTalk.value) {
        if (emcCrossTalkConf.createHistograms.value) {
          for (const auto& cell : cellsBC) {
            mHistManager.fill(HIST("hCellEnergyDistBefore"), cell.getAmplitude());
          }
        }
        emcCrossTalk.setCells(cellsBC, cellLabels);
        bool isOkCrossTalk = emcCrossTalk.run();
        if (!isOkCrossTalk) {
          LOG(info) << "Cross talk emulation failed!";
        } else {
          // When we get new cells we also need to add additional entries into cellIndicesBC.
          // Adding -1 and later when filling the clusterID<->cellID table skip all cases where this is -1
          if (cellIndicesBC.size() < cellsBC.size()) {
            cellIndicesBC.reserve(cellsBC.size());
            size_t nMissing = cellsBC.size() - cellIndicesBC.size();
            cellIndicesBC.insert(cellIndicesBC.end(), nMissing, -1);
          }
          if (emcCrossTalkConf.createHistograms.value) {
            for (const auto& cell : cellsBC) {
              mHistManager.fill(HIST("hCellEnergyDistAfter"), cell.getAmplitude());
            }
          }
        } // cross talk emulation was okay
      } // if (isMC.value && emcCrossTalkConf.enableCrossTalk.value)
      // shaper correction has to come AFTER cross talk
      for (auto& cell : cellsBC) { // o2-linter: disable=const-ref-in-for-loop (we are changing a value here)
        if (cell.getLowGain()) {
          cell.setAmplitude(o2::emcal::NonlinearityHandler::evaluateShaperCorrectionCellEnergy(cell.getAmplitude()));
        }
      }
      LOG(detail) << "Number of cells for BC (CF): " << cellsBC.size();
      nCellsProcessed += cellsBC.size();

      fillQAHistogram(cellsBC);

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

              MatchResult indexMapPair;
              std::vector<int64_t> trackGlobalIndex;
              doTrackMatching<CollEventSels::filtered_iterator>(col, tracks, indexMapPair, trackGlobalIndex);

              // Store the clusters in the table where a matching collision could
              // be identified.
              fillClusterTable<CollEventSels::filtered_iterator>(col, vertexPos, iClusterizer, cellIndicesBC, &indexMapPair, &trackGlobalIndex);
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
        mClusterPhi.clear();
        mClusterEta.clear();
        LOG(debug) << "Cluster loop done for clusterizer " << iClusterizer;
      } // end of clusterizer loop
      LOG(debug) << "Done with process BC.";
      nBCsProcessed++;
    } // end of bc loop

    // Loop through all collisions and fill emcalcollisionmatch with a boolean stating, whether the collision was ambiguous (not the only collision in its BC)
    for (const auto& collision : collisions) {
      if (applySoftwareTriggerSelection) {
        if (!zorro.isSelected(collision.foundBC_as<BcEvSels>().globalBC())) {
          continue;
        }
      }
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

  void processMCWithSecondaries(BcEvSels const& bcs, CollEventSels const& collisions, MyGlobTracks const& tracks, FilteredMcCells const& cells, aod::StoredMcParticles_001 const&, EMV0Legs const& v0legs)
  {
    LOG(debug) << "Starting processMCWithSecondaries.";

    int previousCollisionId = 0; // Collision ID of the last unique BC. Needed to skip unordered collisions to ensure ordered collisionIds in the cluster table
    int nBCsProcessed = 0;
    int nCellsProcessed = 0;
    std::unordered_map<uint64_t, int> numberCollsInBC; // Number of collisions mapped to the global BC index of all BCs
    std::unordered_map<uint64_t, int> numberCellsInBC; // Number of cells mapped to the global BC index of all BCs to check whether EMCal was readout
    for (const auto& bc : bcs) {
      LOG(debug) << "Next BC";
      // Convert aod::Calo to o2::emcal::Cell which can be used with the clusterizer.
      // In particular, we need to filter only EMCAL cells.

      // get run number
      runNumber = bc.runNumber();

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
          const auto& ids = cell.mcParticleIds();
          const auto& amps = cell.amplitudeA();

          if (ids.empty() || amps.empty()) {
            LOGF(warning, "Skipping cell with empty MC info: absId=%d", cell.cellNumber());
            continue;
          }
          mHistManager.fill(HIST("hMCParticleEnergy"), cellparticle.e());
        }
        auto amplitude = cell.amplitude();
        if (static_cast<bool>(hasShaperCorrection) && emcal::intToChannelType(cell.cellType()) == emcal::ChannelType_t::LOW_GAIN) { // Apply shaper correction to LG cells
          amplitude = o2::emcal::NonlinearityHandler::evaluateShaperCorrectionCellEnergy(amplitude);
        }
        if (mcCellEnergyShift != 1.) {
          amplitude *= mcCellEnergyShift; // Fine tune the MC cell energy
        }
        if (mcCellEnergyResolutionBroadening != 0.) {
          amplitude *= (1. + normalgaus(rdgen) * mcCellEnergyResolutionBroadening); // Fine tune the MC cell energy resolution
        }
        cellsBC.emplace_back(cell.cellNumber(),
                             amplitude,
                             cell.time() + getCellTimeShift(cell.cellNumber(), amplitude, o2::emcal::intToChannelType(cell.cellType()), runNumber),
                             o2::emcal::intToChannelType(cell.cellType()));
        cellIndicesBC.emplace_back(cell.globalIndex());
        cellLabels.emplace_back(std::vector<int>{cell.mcParticleIds().begin(), cell.mcParticleIds().end()}, std::vector<float>{cell.amplitudeA().begin(), cell.amplitudeA().end()});
      }
      if (isMC.value && emcCrossTalkConf.enableCrossTalk.value) {
        if (emcCrossTalkConf.createHistograms.value) {
          for (const auto& cell : cellsBC) {
            mHistManager.fill(HIST("hCellEnergyDistBefore"), cell.getAmplitude());
          }
        }
        emcCrossTalk.setCells(cellsBC, cellLabels);
        bool isOkCrossTalk = emcCrossTalk.run();
        if (!isOkCrossTalk) {
          LOG(info) << "Cross talk emulation failed!";
        } else {
          // When we get new cells we also need to add additional entries into cellIndicesBC.
          // Adding -1 and later when filling the clusterID<->cellID table skip all cases where this is -1
          if (cellIndicesBC.size() < cellsBC.size()) {
            cellIndicesBC.reserve(cellsBC.size());
            size_t nMissing = cellsBC.size() - cellIndicesBC.size();
            cellIndicesBC.insert(cellIndicesBC.end(), nMissing, -1);
          }
          if (emcCrossTalkConf.createHistograms.value) {
            for (const auto& cell : cellsBC) {
              mHistManager.fill(HIST("hCellEnergyDistAfter"), cell.getAmplitude());
            }
          }
        } // cross talk emulation was okay
      } // if (isMC.value && emcCrossTalkConf.enableCrossTalk.value)
      // shaper correction has to come AFTER cross talk
      for (auto& cell : cellsBC) { // o2-linter: disable=const-ref-in-for-loop (we are changing a value here)
        if (cell.getLowGain()) {
          cell.setAmplitude(o2::emcal::NonlinearityHandler::evaluateShaperCorrectionCellEnergy(cell.getAmplitude()));
        }
      }
      LOG(detail) << "Number of cells for BC (CF): " << cellsBC.size();
      nCellsProcessed += cellsBC.size();

      fillQAHistogram(cellsBC);

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

              MatchResult indexMapPair;
              std::vector<int64_t> trackGlobalIndex;
              doTrackMatching<CollEventSels::filtered_iterator>(col, tracks, indexMapPair, trackGlobalIndex);

              MatchResult indexMapPairSecondary;
              std::vector<int64_t> secondaryGlobalIndex;
              doSecondaryTrackMatching<CollEventSels::filtered_iterator>(col, v0legs, indexMapPairSecondary, secondaryGlobalIndex, tracks);

              // Store the clusters in the table where a matching collision could
              // be identified.
              fillClusterTable<CollEventSels::filtered_iterator>(col, vertexPos, iClusterizer, cellIndicesBC, &indexMapPair, &trackGlobalIndex, &indexMapPairSecondary, &secondaryGlobalIndex);
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
        mClusterPhi.clear();
        mClusterEta.clear();
        LOG(debug) << "Cluster loop done for clusterizer " << iClusterizer;
      } // end of clusterizer loop
      LOG(debug) << "Done with process BC.";
      nBCsProcessed++;
    } // end of bc loop

    // Loop through all collisions and fill emcalcollisionmatch with a boolean stating, whether the collision was ambiguous (not the only collision in its BC)
    for (const auto& collision : collisions) {
      if (applySoftwareTriggerSelection) {
        if (!zorro.isSelected(collision.foundBC_as<BcEvSels>().globalBC())) {
          continue;
        }
      }
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
  PROCESS_SWITCH(EmcalCorrectionTask, processMCWithSecondaries, "run full analysis with MC info", false);

  void processStandalone(BcEvSels const& bcs, aod::Collisions const& collisions, FilteredCells const& cells)
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

      // get run number
      runNumber = bc.runNumber();

      if (applyTempCalib && !mIsTempCalibInitialized) { // needs to be called once
        mTempCalibExtractor->InitializeFromCCDB(pathTempCalibCCDB, static_cast<uint64_t>(runNumber));
        mIsTempCalibInitialized = true;
      }

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
        auto amplitude = cell.amplitude();
        if (static_cast<bool>(hasShaperCorrection) && emcal::intToChannelType(cell.cellType()) == emcal::ChannelType_t::LOW_GAIN) { // Apply shaper correction to LG cells
          amplitude = o2::emcal::NonlinearityHandler::evaluateShaperCorrectionCellEnergy(amplitude);
        }
        if (applyGainCalibShift) {
          amplitude *= mArrGainCalibDiff[cell.cellNumber()];
        }
        if (applyTempCalib) {
          float tempCalibFactor = mTempCalibExtractor->getGainCalibFactor(static_cast<uint16_t>(cell.cellNumber()));
          amplitude /= tempCalibFactor;
          mHistManager.fill(HIST("hTempCalibCorrection"), tempCalibFactor);
        }
        cellsBC.emplace_back(cell.cellNumber(),
                             amplitude,
                             cell.time() + getCellTimeShift(cell.cellNumber(), amplitude, o2::emcal::intToChannelType(cell.cellType()), runNumber),
                             o2::emcal::intToChannelType(cell.cellType()));
        cellIndicesBC.emplace_back(cell.globalIndex());
      }
      LOG(detail) << "Number of cells for BC (CF): " << cellsBC.size();
      nCellsProcessed += cellsBC.size();

      fillQAHistogram(cellsBC);

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
          fillAmbigousClusterTable<BcEvSels::iterator>(bc, iClusterizer, cellIndicesBC, hasCollision);
        }

        mClusterPhi.clear();
        mClusterEta.clear();
        LOG(debug) << "Cluster loop done for clusterizer " << iClusterizer;
      } // end of clusterizer loop
      LOG(detail) << "Processed " << nBCsProcessed << " BCs with " << nCellsProcessed << " cells";
      nBCsProcessed++;
    } // end of bc loop
    LOG(debug) << "Done with process BC.";
  }
  PROCESS_SWITCH(EmcalCorrectionTask, processStandalone, "run stand alone analysis", false);

  void cellsToCluster(size_t iClusterizer, const gsl::span<o2::emcal::Cell> cellsBC, gsl::span<const o2::emcal::CellLabel> cellLabels = {})
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
    if (cellLabels.empty()) {
      mClusterFactories.setContainer(*emcalClusters, cellsBC, *emcalClustersInputIndices);
    } else {
      mClusterFactories.setContainer(*emcalClusters, cellsBC, *emcalClustersInputIndices, cellLabels);
    }

    LOG(debug) << "Cluster factory set up.";
    // Convert to analysis clusters.
    for (int icl = 0; icl < mClusterFactories.getNumberOfClusters(); icl++) {
      o2::emcal::ClusterLabel clusterLabel;
      auto analysisCluster = mClusterFactories.buildCluster(icl, &clusterLabel);
      mAnalysisClusters.emplace_back(analysisCluster);
      mClusterLabels.push_back(clusterLabel);
      auto pos = analysisCluster.getGlobalPosition();
      mClusterPhi.emplace_back(RecoDecay::constrainAngle(pos.Phi()));
      mClusterEta.emplace_back(pos.Eta());
      LOG(debug) << "Cluster " << icl << ": E: " << analysisCluster.E() << ", NCells " << analysisCluster.getNCells();
    }
    mHistManager.fill(HIST("hNCluster"), mAnalysisClusters.size());
    LOG(debug) << "Converted to analysis clusters.";
  }

  template <typename Collision>
  void fillClusterTable(Collision const& col, math_utils::Point3D<float> const& vertexPos, size_t iClusterizer, const gsl::span<int64_t> cellIndicesBC, MatchResult* indexMapPair = nullptr, const std::vector<int64_t>* trackGlobalIndex = nullptr, MatchResult* indexMapPairSecondaries = nullptr, const std::vector<int64_t>* secondariesGlobalIndex = nullptr)
  {
    // average number of cells per cluster, only used the reseve a reasonable amount for the clustercells table
    const size_t nAvgNcells = 3;
    // we found a collision, put the clusters into the none ambiguous table
    clusters.reserve(mAnalysisClusters.size());
    if (!mClusterLabels.empty()) {
      mcclusters.reserve(mClusterLabels.size());
    }
    clustercells.reserve(mAnalysisClusters.size() * nAvgNcells);

    // get the clusterType once
    const auto clusterType = static_cast<int>(mClusterDefinitions[iClusterizer]);

    int cellindex = -1;
    unsigned int iCluster = 0;
    float energy = 0.f;
    for (const auto& cluster : mAnalysisClusters) {
      energy = cluster.E();
      // Determine the cluster eta, phi, correcting for the vertex position.
      auto pos = cluster.getGlobalPosition();
      pos = pos - vertexPos;
      // Normalize the vector and rescale by energy.
      pos *= (energy / std::sqrt(pos.Mag2()));

      // Correct for nonlinear behaviour
      float nonlinCorrEnergy = energy;
      if (!disableNonLin) {
        try {
          nonlinCorrEnergy = mNonlinearityHandler.getCorrectedClusterEnergy(cluster);
        } catch (o2::emcal::NonlinearityHandler::UninitException& e) {
          LOG(error) << e.what();
        }
      }

      // save to table
      LOG(debug) << "Writing cluster definition "
                 << clusterType
                 << " to table.";
      mHistManager.fill(HIST("hClusterType"), 1);
      clusters(col, cluster.getID(), nonlinCorrEnergy, cluster.getCoreEnergy(), energy,
               pos.Eta(), RecoDecay::constrainAngle(pos.Phi()), cluster.getM02(),
               cluster.getM20(), cluster.getNCells(),
               cluster.getClusterTime(), cluster.getIsExotic(),
               cluster.getDistanceToBadChannel(), cluster.getNExMax(),
               clusterType);
      if (!mClusterLabels.empty()) {
        mcclusters(mClusterLabels[iCluster].getLabels(), mClusterLabels[iCluster].getEnergyFractions());
      }
      // loop over cells in cluster and save to table
      for (int ncell = 0; ncell < cluster.getNCells(); ncell++) {
        cellindex = cluster.getCellIndex(ncell);
        LOG(debug) << "trying to find cell index " << cellindex << " in map";
        if (cellIndicesBC[cellindex] >= 0) {
          clustercells(clusters.lastIndex(), cellIndicesBC[cellindex]);
        }
      } // end of cells of cluser loop
      // fill histograms
      mHistManager.fill(HIST("hClusterE"), energy);
      mHistManager.fill(HIST("hClusterNLM"), cluster.getNExMax());
      mHistManager.fill(HIST("hClusterTime"), cluster.getClusterTime());
      mHistManager.fill(HIST("hClusterEtaPhi"), pos.Eta(), RecoDecay::constrainAngle(pos.Phi()));
      if (fillQA.value) {
        mHistManager.fill(HIST("hClusterNCellE"), cluster.E(), cluster.getNCells());
        mHistManager.fill(HIST("hClusterFCrossE"), cluster.E(), cluster.getFCross());
        mHistManager.fill(HIST("hClusterFCrossSigmaLongE"), cluster.E(), cluster.getFCross(), cluster.getM02());
        mHistManager.fill(HIST("hClusterFCrossSigmaShortE"), cluster.E(), cluster.getFCross(), cluster.getM20());
      }
      if (indexMapPair && trackGlobalIndex) {
        if (iCluster < indexMapPair->matchIndexTrack.size() && indexMapPair->matchIndexTrack.size() > 0) {
          for (unsigned int iTrack = 0; iTrack < indexMapPair->matchIndexTrack[iCluster].size(); iTrack++) {
            if (indexMapPair->matchIndexTrack[iCluster][iTrack] >= 0) {
              LOG(debug) << "Found track " << (*trackGlobalIndex)[indexMapPair->matchIndexTrack[iCluster][iTrack]] << " in cluster " << cluster.getID();
              matchedTracks(clusters.lastIndex(), (*trackGlobalIndex)[indexMapPair->matchIndexTrack[iCluster][iTrack]], indexMapPair->matchDeltaPhi[iCluster][iTrack], indexMapPair->matchDeltaEta[iCluster][iTrack]);
              mHistManager.fill(HIST("hMatchedPrimaryTracks"), indexMapPair->matchDeltaEta[iCluster][iTrack], indexMapPair->matchDeltaPhi[iCluster][iTrack]);
            }
          }
        }
      }
      if (indexMapPairSecondaries && secondariesGlobalIndex) {
        if (iCluster < indexMapPairSecondaries->matchIndexTrack.size() && indexMapPairSecondaries->matchIndexTrack.size() > 0) {
          for (unsigned int iTrack = 0; iTrack < indexMapPairSecondaries->matchIndexTrack[iCluster].size(); iTrack++) {
            if (indexMapPairSecondaries->matchIndexTrack[iCluster][iTrack] >= 0) {
              LOG(debug) << "Found secondary track " << (*secondariesGlobalIndex)[indexMapPairSecondaries->matchIndexTrack[iCluster][iTrack]] << " in cluster " << cluster.getID();
              matchedSecondaries(clusters.lastIndex(), (*secondariesGlobalIndex)[indexMapPairSecondaries->matchIndexTrack[iCluster][iTrack]], indexMapPairSecondaries->matchDeltaPhi[iCluster][iTrack], indexMapPairSecondaries->matchDeltaEta[iCluster][iTrack]);
              mHistManager.fill(HIST("hMatchedSecondaries"), indexMapPairSecondaries->matchDeltaEta[iCluster][iTrack], indexMapPairSecondaries->matchDeltaPhi[iCluster][iTrack]);
            }
          }
        }
      }
      iCluster++;
    } // end of cluster loop
  }

  template <typename BC>
  void fillAmbigousClusterTable(BC const& bc, size_t iClusterizer, const gsl::span<int64_t> cellIndicesBC, bool hasCollision)
  {
    // average number of cells per cluster, only used the reseve a reasonable amount for the clustercells table
    const size_t nAvgNcells = 3;
    int cellindex = -1;
    clustersAmbiguous.reserve(mAnalysisClusters.size());
    if (mClusterLabels.size() > 0) {
      mcclustersAmbiguous.reserve(mClusterLabels.size());
    }
    clustercellsambiguous.reserve(mAnalysisClusters.size() * nAvgNcells);
    unsigned int iCluster = 0;
    float energy = 0.f;
    for (const auto& cluster : mAnalysisClusters) {
      energy = cluster.E();
      auto pos = cluster.getGlobalPosition();
      pos = pos - math_utils::Point3D<float>{0., 0., 0.};
      // Normalize the vector and rescale by energy.
      pos *= (energy / std::sqrt(pos.Mag2()));

      // Correct for nonlinear behaviour
      float nonlinCorrEnergy = energy;
      try {
        nonlinCorrEnergy = mNonlinearityHandler.getCorrectedClusterEnergy(cluster);
      } catch (o2::emcal::NonlinearityHandler::UninitException& e) {
        LOG(error) << e.what();
      }

      // We have our necessary properties. Now we store outputs

      // LOG(debug) << "Cluster E: " << energy;
      if (!hasCollision) {
        mHistManager.fill(HIST("hClusterType"), 0);
      } else {
        mHistManager.fill(HIST("hClusterType"), 2);
      }
      clustersAmbiguous(
        bc, cluster.getID(), nonlinCorrEnergy, cluster.getCoreEnergy(), energy,
        pos.Eta(), RecoDecay::constrainAngle(pos.Phi()), cluster.getM02(),
        cluster.getM20(), cluster.getNCells(), cluster.getClusterTime(),
        cluster.getIsExotic(), cluster.getDistanceToBadChannel(),
        cluster.getNExMax(), static_cast<int>(mClusterDefinitions.at(iClusterizer)));
      if (mClusterLabels.size() > 0) {
        mcclustersAmbiguous(mClusterLabels[iCluster].getLabels(), mClusterLabels[iCluster].getEnergyFractions());
      }
      for (int ncell = 0; ncell < cluster.getNCells(); ncell++) {
        cellindex = cluster.getCellIndex(ncell);
        clustercellsambiguous(clustersAmbiguous.lastIndex(),
                              cellIndicesBC[cellindex]);
      } // end of cells of cluster loop
      iCluster++;
    } // end of cluster loop
  }

  template <typename Collision>
  void doTrackMatching(Collision const& col, MyGlobTracks const& tracks, MatchResult& indexMapPair, std::vector<int64_t>& trackGlobalIndex)
  {
    auto groupedTracks = tracks.sliceBy(perCollision, col.globalIndex());
    int nTracksInCol = groupedTracks.size();
    std::vector<float> trackPhi;
    std::vector<float> trackEta;
    // reserve memory to reduce on the fly memory allocation
    trackPhi.reserve(nTracksInCol);
    trackEta.reserve(nTracksInCol);
    trackGlobalIndex.reserve(nTracksInCol);
    fillTrackInfo<decltype(groupedTracks)>(groupedTracks, trackPhi, trackEta, trackGlobalIndex);

    indexMapPair = matchTracksToCluster(mClusterPhi, mClusterEta, trackPhi, trackEta, maxMatchingDistance, kMaxMatchesPerCluster);
  }

  template <typename Collision>
  void doSecondaryTrackMatching(Collision const& col, EMV0Legs const& v0legs, MatchResult& indexMapPair, std::vector<int64_t>& trackGlobalIndex, MyGlobTracks const& tracks)
  {
    auto groupedV0Legs = v0legs.sliceBy(perCollisionEMV0Legs, col.globalIndex());
    int nLegsInCol = groupedV0Legs.size();
    std::vector<float> trackPhi;
    std::vector<float> trackEta;
    // reserve memory to reduce on the fly memory allocation
    trackPhi.reserve(nLegsInCol);
    trackEta.reserve(nLegsInCol);
    trackGlobalIndex.reserve(nLegsInCol);

    float trackEtaEmcal = 0.f;
    float trackPhiEmcal = 0.f;
    for (const auto& leg : groupedV0Legs) {
      if (leg.trackId() < 0 || leg.trackId() > tracks.size()) {
        continue;
      }
      auto track = tracks.iteratorAt(leg.trackId());
      trackEtaEmcal = track.trackEtaEmcal();
      trackPhiEmcal = track.trackPhiEmcal();
      // Tracks that do not point to the EMCal/DCal/PHOS get default values of -999
      // This way we can cut out tracks that do not point to the EMCal+DCal
      if (trackEtaEmcal < TrackNotOnEMCal || trackPhiEmcal < TrackNotOnEMCal) {
        continue;
      }
      if (trackMinPt > 0 && track.pt() < trackMinPt) {
        continue;
      }
      trackPhi.emplace_back(RecoDecay::constrainAngle(trackPhiEmcal));
      trackEta.emplace_back(trackEtaEmcal);
      trackGlobalIndex.emplace_back(track.globalIndex());
    }
    indexMapPair = matchTracksToCluster(mClusterPhi, mClusterEta, trackPhi, trackEta, maxMatchingDistance, kMaxMatchesPerCluster);
  }

  template <typename Tracks>
  void fillTrackInfo(Tracks const& tracks, std::vector<float>& trackPhi, std::vector<float>& trackEta, std::vector<int64_t>& trackGlobalIndex)
  {
    for (const auto& track : tracks) {
      if (!track.isGlobalTrack()) { // only global tracks
        continue;
      }
      // Tracks that do not point to the EMCal/DCal/PHOS get default values of -999
      // This way we can cut out tracks that do not point to the EMCal+DCal
      if (track.trackEtaEmcal() < TrackNotOnEMCal || track.trackPhiEmcal() < TrackNotOnEMCal) {
        continue;
      }
      if (trackMinPt > 0 && track.pt() < trackMinPt) {
        continue;
      }
      trackPhi.emplace_back(RecoDecay::constrainAngle(track.trackPhiEmcal()));
      trackEta.emplace_back(track.trackEtaEmcal());
      trackGlobalIndex.emplace_back(track.globalIndex());
    }
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
      mHistManager.fill(HIST("hCellEtaPhi"), std::get<0>(res), RecoDecay::constrainAngle(std::get<1>(res)));
      res = mClusterizers.at(0)->getGeometry()->GlobalRowColFromIndex(cell.getTower());
      // NOTE: Reversed column and row because it's more natural for presentation.
      mHistManager.fill(HIST("hCellRowCol"), std::get<1>(res), std::get<0>(res));
    }
  }

  float getAbsCellScale(const int cellID)
  {
    // Apply cell scale based on SM types (Full, Half (not used), EMC 1/3, DCal, DCal 1/3)
    // Same as in Run2 data
    if (applyCellAbsScale == CellScaleMode::ModeSMWise) {
      int iSM = mClusterizers.at(0)->getGeometry()->GetSuperModuleNumber(cellID);
      return cellAbsScaleFactors.value[mClusterizers.at(0)->getGeometry()->GetSMType(iSM)];

      // Apply cell scale based on columns to accoutn for material of TRD structures
    } else if (applyCellAbsScale == CellScaleMode::ModeColumnWise) {
      auto res = mClusterizers.at(0)->getGeometry()->GlobalRowColFromIndex(cellID);
      return cellAbsScaleFactors.value[std::get<1>(res)];
    } else {
      return 1.f;
    }
  }

  // Apply shift of the cell time in data and MC
  // In MC this has to be done to shift the cell time, which is not calibrated to 0 due to the flight time of the particles to the EMCal surface (~15ns)
  // In data this is done to correct for the time walk effect
  float getCellTimeShift(const int16_t cellID, const float cellEnergy, const emcal::ChannelType_t cellType, const int runNumber)
  {
    if (!applyCellTimeCorrection) {
      return 0.f;
    }
    float timeshift = 0.f;
    float timesmear = 0.f;
    const float minLeaderEnergy = 0.3f;
    const float lowEnergyRegime = 4.f;
    const float highEnergyRegime = 30.f;
    if (isMC) { // ---> MC
      // Shift the time to 0, as the TOF was simulated -> eta dependent shift (as larger eta values are further away from collision point)
      // Use distance between vertex and EMCal (at eta = 0) and distance on EMCal surface (cell size times column) to calculate distance to cell
      // 0.2 is cell size in m (0.06) divided by the speed of light in m/ns (0.3) - 47.5 is the "middle" of the EMCal (2*48 cells in one column)
      float timeCol = 0.2f * (geometry->GlobalCol(cellID) - 47.5f); // calculate time to get to specific column
      timeshift = -std::sqrt(215.f + timeCol * timeCol);            // 215 is 14.67ns^2 (time it takes to get the cell at eta = 0)

      // Also smear the time to account for the broader time resolution in data than in MC
      if (cellEnergy < minLeaderEnergy)                                           // Cells with tless than 300 MeV cannot be the leading cell in the cluster, so their time does not require precise calibration
        timesmear = 0.;                                                           // They will therefore not be smeared and only get their shift
      else if (cellType == emcal::ChannelType_t::HIGH_GAIN)                       // High gain cells -> Low energies
        timesmear = normalgaus(rdgen) * (1.6 + 9.5 * std::exp(-3. * cellEnergy)); // Parameters extracted from LHC24f3b & LHC22o (pp), but also usable for other periods
      else if (cellType == emcal::ChannelType_t::LOW_GAIN)                        // Low gain cells -> High energies
        timesmear = normalgaus(rdgen) * (5.0);                                    // Parameters extracted from LHC24g4 & LHC24aj (pp), but also usable for other periods

    } else {                                                    // ---> Data
      if (cellEnergy < minLeaderEnergy) {                       // Cells with tless than 300 MeV cannot be the leading cell in the cluster, so their time does not require precise calibration
        timeshift = 0.;                                         // In data they will not be shifted (they are close to 0 anyways)
      } else if (cellType == emcal::ChannelType_t::HIGH_GAIN) { // High gain cells -> Low energies
        if (cellEnergy < lowEnergyRegime)                       // Low energy regime
          timeshift = 0.8 * std::log(2.7 * cellEnergy);         // Parameters extracted from LHC22o (pp), but also usable for other periods
        else                                                    // Medium energy regime
          timeshift = 1.5 * std::log(0.9 * cellEnergy);         // Parameters extracted from LHC22o (pp), but also usable for other periods
      } else if (cellType == emcal::ChannelType_t::LOW_GAIN) {  // Low gain cells -> High energies
        if (cellEnergy < highEnergyRegime)                      // High energy regime
          timeshift = 1.9 * std::log(0.09 * cellEnergy);        // Parameters extracted from LHC24aj (pp), but also usable for other periods
        else                                                    // Very high energy regime
          timeshift = 1.9;                                      // Parameters extracted from LHC24aj (pp), but also usable for other periods
      }
      // Temporary extra shift for bug in time calibraiton of apass4 Pb-Pb 2024, requires pos shift of 2*8.8 ns for low gain cells
      if (cellType == emcal::ChannelType_t::LOW_GAIN) {
        for (const auto& range : mExtraTimeShiftRunRanges) {
          if (runNumber >= range.first && runNumber <= range.second) {
            timeshift += 2 * 8.8;
          }
        }
      }
      LOG(debug) << "Shift the cell time by " << timeshift << " + " << timesmear << " ns";
    }
    return timeshift + timesmear;
  };

  void initializeGainCalibShift()
  {
    auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
    uint64_t tsOld = 1634853602000; // timestamp corresponding to LHC22o old gain calib object
    o2::emcal::GainCalibrationFactors* paramsOld = ccdbMgr.getForTimeStamp<o2::emcal::GainCalibrationFactors>("EMC/Calib/GainCalibFactors", tsOld);
    uint64_t tsNew = 1734853602000; // timestamp corresponding to new gain calib object (new cell compression)
    o2::emcal::GainCalibrationFactors* paramsNew = ccdbMgr.getForTimeStamp<o2::emcal::GainCalibrationFactors>("EMC/Calib/GainCalibFactors", tsNew);
    for (uint16_t i = 0; i < mArrGainCalibDiff.size(); ++i) {
      mArrGainCalibDiff[i] = paramsNew->getGainCalibFactors(i) == 0 ? 1. : paramsOld->getGainCalibFactors(i) / paramsNew->getGainCalibFactors(i);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<EmcalCorrectionTask>(cfgc)};
}
