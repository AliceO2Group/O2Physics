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

/// \file taskPi0FlowEMC.cxx
/// \brief Analysis task for neutral pion flow with EMCal
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/PhotonMeson/Core/EMBitFlags.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <EMCALBase/Geometry.h>
#include <EMCALBase/GeometryBase.h>
#include <EMCALCalib/BadChannelMap.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/AxisAngle.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TF1.h>
#include <TH1.h>
#include <TString.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;

enum QvecEstimator {
  FT0M = 0,
  FT0A = 1,
  FT0C,
  TPCPos,
  TPCNeg,
  TPCTot,
  FV0A
};

enum CentralityEstimator {
  None = 0,
  CFT0A = 1,
  CFT0C,
  CFT0M,
  NCentralityEstimators
};

enum Harmonics {
  kNone = 0,
  kDirect = 1,
  kElliptic = 2,
  kTriangluar = 3,
  kQuadrangular = 4,
  kPentagonal = 5,
  kHexagonal = 6,
  kHeptagonal = 7,
  kOctagonal = 8
};

enum class MapLevel {
  kGood = 1,
  kNoBad = 2,
  kInEMC = 3,
  kAll = 4
};

struct TaskPi0FlowEMC {
  static constexpr float MinEnergy = 0.7f;

  // configurable for flow
  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<int> qvecDetector{"qvecDetector", 0, "Detector for Q vector estimation (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5, FV0A: 6)"};
  Configurable<int> qvecSubADetector{"qvecSubADetector", 3, "Sub A Detector for Q vector estimation for resolution (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5, FV0A: 6)"};
  Configurable<int> qvecSubBDetector{"qvecSubBDetector", 4, "Sub B Detector for Q vector estimation for resolution (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5, FV0A: 6)"};
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int> cfgEMCalMapLevelBackground{"cfgEMCalMapLevelBackground", 4, "Different levels of correction for the background, the smaller number includes the level of the higher number (4: none, 3: only inside EMCal, 2: remove edges, 1: exclude bad channels)"};
  Configurable<int> cfgEMCalMapLevelSameEvent{"cfgEMCalMapLevelSameEvent", 4, "Different levels of correction for the same event, the smaller number includes the level of the higher number (4: none, 3: only inside EMCal, 2: remove edges, 1: exclude bad channels)"};
  Configurable<int> cfgDistanceToEdge{"cfgDistanceToEdge", 1, "Distance to edge in cells required for rotated cluster to be accepted"};
  Configurable<bool> cfgDoM02{"cfgDoM02", false, "Flag to enable flow vs M02 for single photons"};
  Configurable<bool> cfgDoReverseScaling{"cfgDoReverseScaling", false, "Flag to reverse the scaling that is possibly applied during NonLin"};
  Configurable<bool> cfgDoPlaneQA{"cfgDoPlaneQA", false, "Flag to enable QA plots comparing in and out of plane"};
  Configurable<float> cfgMaxQVector{"cfgMaxQVector", 20.f, "Maximum allowed absolute QVector value."};
  Configurable<float> cfgMaxAsymmetry{"cfgMaxAsymmetry", 0.1f, "Maximum allowed asymmetry for photon pairs used in calibration."};

  // configurable axis
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {400, 0.0, 0.8}, "invariant mass axis for the neutral meson"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0., 20.}, "pT axis for the neutral meson"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {20, 0., 100.}, "centrality axis for the current event"};
  ConfigurableAxis thnConfigAxisCosNPhi{"thnConfigAxisCosNPhi", {100, -1., 1.}, "cos(n*phi) axis for the current event"};
  ConfigurableAxis thnConfigAxisCosDeltaPhi{"thnConfigAxisCosDeltaPhi", {8, -1., 1.}, "cos(delta phi) axis for the current event"};
  ConfigurableAxis thnConfigAxisM02{"thnConfigAxisM02", {200, 0., 5.}, "M02 axis for the EMCal cluster"};
  ConfigurableAxis thnConfigAxisEnergyCalib{"thnConfigAxisEnergyCalib", {200, 0., 20.}, "energy axis for the emcal clusters for the calibration process"};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcuts";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<bool> cfgRequireEMCReadoutInMB{"cfgRequireEMCReadoutInMB", true, "require the EMC to be read out in an MB collision (kTVXinEMC)"};
    Configurable<bool> cfgRequireEMCHardwareTriggered{"cfgRequireEMCHardwareTriggered", false, "require the EMC to be hardware triggered (kEMC7 or kDMC7)"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -1, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<float> cfgMinCent{"cfgMinCent", 0, "min. centrality (%)"};
    Configurable<float> cfgMaxCent{"cfgMaxCent", 90, "max. centrality (%)"};
    Configurable<bool> onlyKeepWeightedEvents{"onlyKeepWeightedEvents", false, "flag to keep only weighted events (for JJ MCs) and remove all MB events (with weight = 1)"};
  } eventcuts;

  EMCPhotonCut fEMCCut;
  struct : ConfigurableGroup {
    std::string prefix = "emccuts";
    Configurable<std::string> clusterDefinition{"clusterDefinition", "kV3MostSplitSmallestTimeDiff", "Clusterizer to be selected, e.g. V3Default"};
    Configurable<float> cfgEMCminTime{"cfgEMCminTime", -25., "Minimum cluster time for EMCal time cut"};
    Configurable<float> cfgEMCmaxTime{"cfgEMCmaxTime", +30., "Maximum cluster time for EMCal time cut"};
    Configurable<float> cfgEMCminM02{"cfgEMCminM02", 0.1, "Minimum M02 for EMCal M02 cut"};
    Configurable<float> cfgEMCmaxM02{"cfgEMCmaxM02", 0.7, "Maximum M02 for EMCal M02 cut"};
    Configurable<float> cfgEMCminE{"cfgEMCminE", 0.7, "Minimum cluster energy for EMCal energy cut"};
    Configurable<int> cfgEMCminNCell{"cfgEMCminNCell", 1, "Minimum number of cells per cluster for EMCal NCell cut"};
    Configurable<std::vector<float>> cfgEMCTMEta{"cfgEMCTMEta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> cfgEMCTMPhi{"cfgEMCTMPhi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> emcSecTMEta{"emcSecTMEta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> emcSecTMPhi{"emcSecTMPhi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<float> cfgEMCEoverp{"cfgEMCEoverp", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
    Configurable<bool> cfgEMCUseExoticCut{"cfgEMCUseExoticCut", true, "FLag to use the EMCal exotic cluster cut"};
    Configurable<bool> cfgEMCUseTM{"cfgEMCUseTM", false, "flag to use EMCal track matching cut or not"};
    Configurable<bool> emcUseSecondaryTM{"emcUseSecondaryTM", false, "flag to use EMCal secondary track matching cut or not"};
    Configurable<bool> cfgEnableQA{"cfgEnableQA", false, "flag to turn QA plots on/off"};
  } emccuts;

  V0PhotonCut fV0PhotonCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "PCMcuts";
    o2::framework::Configurable<bool> requireV0WithITSTPC{"requireV0WithITSTPC", false, "flag to enforce V0s have ITS and TPC"};
    o2::framework::Configurable<bool> requireV0WithITSOnly{"requireV0WithITSOnly", false, "flag to select V0s with ITSonly tracks"};
    o2::framework::Configurable<bool> requireV0WithTPCOnly{"requireV0WithTPCOnly", false, "flag to select V0s with TPConly tracks"};
    o2::framework::Configurable<float> minPtV0{"minPtV0", 0.1, "min pT for v0 photons at PV"};
    o2::framework::Configurable<float> maxPtV0{"maxPtV0", 1e+10, "max pT for v0 photons at PV"};
    o2::framework::Configurable<float> minEtaV0{"minEtaV0", -0.8, "min eta for v0 photons at PV"};
    o2::framework::Configurable<float> maxEtaV0{"maxEtaV0", 0.8, "max eta for v0 photons at PV"};
    o2::framework::Configurable<float> minRV0{"minRV0", 4.0, "min v0 radius"};
    o2::framework::Configurable<float> maxRV0{"maxRV0", 90.0, "max v0 radius"};
    o2::framework::Configurable<float> maxAlphaAP{"maxAlphaAP", 0.95, "max alpha for AP cut"};
    o2::framework::Configurable<float> maxQtAP{"maxQtAP", 0.01, "max qT for AP cut"};
    o2::framework::Configurable<float> minCosPA{"minCosPA", 0.999, "min V0 CosPA"};
    o2::framework::Configurable<float> maxPCA{"maxPCA", 1.5, "max distance btween 2 legs"};
    o2::framework::Configurable<float> maxChi2KF{"maxChi2KF", 1.e+10f, "max chi2/ndf with KF"};
    o2::framework::Configurable<bool> rejectV0onITSib{"rejectV0onITSib", true, "flag to reject V0s on ITSib"};
    o2::framework::Configurable<bool> applyPrefilter{"applyPrefilter", false, "flag to apply prefilter to V0"};

    o2::framework::Configurable<int> minNClusterTPC{"minNClusterTPC", 0, "min NCluster TPC"};
    o2::framework::Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 40, "min ncrossed rows in TPC"};
    o2::framework::Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "min fraction of crossed rows over findable clusters in TPC"};
    o2::framework::Configurable<float> maxFracSharedClustersTPC{"maxFracSharedClustersTPC", 999.f, "max fraction of shared clusters in TPC"};
    o2::framework::Configurable<int> minNClusterITS{"minNClusterITS", 0, "min NCluster ITS"};
    o2::framework::Configurable<float> minMeanClusterSizeITSob{"minMeanClusterSizeITSob", 0.f, "min <cluster size> ITSob"};
    o2::framework::Configurable<float> maxMeanClusterSizeITSob{"maxMeanClusterSizeITSob", 16.f, "max <cluster size> ITSob"};
    o2::framework::Configurable<float> macChi2TPC{"macChi2TPC", 4.f, "max chi2/NclsTPC"};
    o2::framework::Configurable<float> macChi2ITS{"macChi2ITS", 36.f, "max chi2/NclsITS"};
    o2::framework::Configurable<float> minTPCNSigmaEl{"minTPCNSigmaEl", -3.0f, "min. TPC n sigma for electron"};
    o2::framework::Configurable<float> maxTPCNSigmaEl{"maxTPCNSigmaEl", +3.0f, "max. TPC n sigma for electron"};
    o2::framework::Configurable<bool> disableITSOnly{"disableITSOnly", false, "flag to disable ITSonly tracks"};
    o2::framework::Configurable<bool> disableTPCOnly{"disableTPCOnly", false, "flag to disable TPConly tracks"};
    o2::framework::Configurable<bool> doQA{"doQA", false, "flag to set QA flag."};
  } pcmcuts;

  struct : ConfigurableGroup {
    std::string prefix = "mesonConfig";
    Configurable<float> minOpenAngle{"minOpenAngle", 0.0202, "apply min opening angle. Default value one EMCal cell"};
    Configurable<bool> enableTanThetadPhi{"enableTanThetadPhi", false, "flag to turn cut opening angle in delta theta delta phi on/off"};
    Configurable<float> minTanThetadPhi{"minTanThetadPhi", 4., "apply min opening angle in delta theta delta phi to cut on late conversion"};
    Configurable<float> maxEnergyAsymmetry{"maxEnergyAsymmetry", 1., "apply max energy asymmetry for meson candidate"};
    Configurable<bool> cfgEnableQA{"cfgEnableQA", false, "flag to turn QA plots on/off"};
    ConfigurableAxis thConfigAxisTanThetaPhi{"thConfigAxisTanThetaPhi", {180, -90.f, 90.f}, ""};
  } mesonConfig;

  struct : ConfigurableGroup {
    std::string prefix = "mixingConfig";
    ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis cfgCentBins{"cfgCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "Mixing bins - centrality"};
    ConfigurableAxis cfgEPBins{"cfgEPBins", {8, o2::constants::math::PIHalf, o2::constants::math::PIHalf}, "Mixing bins - event plane angle"};
    ConfigurableAxis cfgOccupancyBins{"cfgOccupancyBins", {VARIABLE_WIDTH, 0, 100, 500, 1000, 2000}, "Mixing bins - occupancy"};
    Configurable<int> cfgMixingDepth{"cfgMixingDepth", 2, "Mixing depth"};
  } mixingConfig;

  struct : ConfigurableGroup {
    std::string prefix = "rotationConfig";
    Configurable<bool> cfgDoRotation{"cfgDoRotation", false, "Flag to enable rotation background method."};
    Configurable<int> cfgDownsampling{"cfgDownsampling", 1, "Calculate rotation background only for every <value> collision."};
    Configurable<float> cfgRotAngle{"cfgRotAngle", std::move(const_cast<float&>(o2::constants::math::PIHalf)), "Angle used for the rotation method."};
    Configurable<bool> cfgUseWeights{"cfgUseWeights", false, "Flag to enable weights for rotation background method."};
  } rotationConfig;

  struct : ConfigurableGroup {
    std::string prefix = "correctionConfig";
    Configurable<std::string> cfgSpresoPath{"cfgSpresoPath", "Users/m/mhemmer/EM/Flow/Resolution", "Path to SP resolution file"};
    Configurable<int> cfgApplySPresolution{"cfgApplySPresolution", 0, "Apply resolution correction"};
    Configurable<bool> doEMCalCalib{"doEMCalCalib", 0, "Produce output for EMCal calibration"};
    Configurable<bool> cfgEnableNonLin{"cfgEnableNonLin", false, "flag to turn extra non linear energy calibration on/off"};
  } correctionConfig;

  SliceCache cache;
  EventPlaneHelper epHelper;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  int runNow = 0;
  int runBefore = -1;

  // Filter clusterFilter = aod::skimmedcluster::time >= emccuts.cfgEMCminTime && aod::skimmedcluster::time <= emccuts.cfgEMCmaxTime && aod::skimmedcluster::m02 >= emccuts.cfgEMCminM02 && aod::skimmedcluster::m02 <= emccuts.cfgEMCmaxM02 && aod::skimmedcluster::e >= emccuts.cfgEMCminE;
  Filter collisionFilter = (nabs(aod::collision::posZ) <= eventcuts.cfgZvtxMax) && (aod::evsel::ft0cOccupancyInTimeRange <= eventcuts.cfgFT0COccupancyMax) && (aod::evsel::ft0cOccupancyInTimeRange >= eventcuts.cfgFT0COccupancyMin);
  // using FilteredEMCalPhotons = soa::Filtered<soa::Join<aod::EMCEMEventIds, aod::MinClusters>>;
  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::MinClusters>;
  using PCMPhotons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
  using FilteredCollsWithQvecs = soa::Filtered<soa::Join<aod::EMEvents, aod::EMEventsAlias, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>>;
  using CollsWithQvecs = soa::Join<aod::EMEvents, aod::EMEventsAlias, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
  using Colls = soa::Join<aod::EMEvents, aod::EMEventsAlias, aod::EMEventsMult, aod::EMEventsCent>;

  static constexpr std::size_t NQVecEntries = 6;

  PresliceOptional<EMCalPhotons> perCollisionEMC = o2::aod::emccluster::emeventId;
  PresliceOptional<PCMPhotons> perCollisionPCM = aod::v0photonkf::emeventId;
  PresliceOptional<MinMTracks> perEMCClusterMT = o2::aod::mintm::minClusterId;
  PresliceOptional<MinMSTracks> perEMCClusterMS = o2::aod::mintm::minClusterId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  o2::emcal::Geometry* emcalGeom;
  o2::emcal::BadChannelMap* mBadChannels;
  TH1D* h1SPResolution = nullptr;
  // Constants for eta and phi ranges for the look up table
  static constexpr double EtaMin = -0.75, etaMax = 0.75;
  static constexpr int NBinsEta = 150; // 150 bins for eta

  static constexpr double PhiMin = 1.35, phiMax = 5.75;
  static constexpr int NBinsPhi = 440; // (440 bins = 0.01 step size covering most regions)

  std::array<int8_t, NBinsEta * NBinsPhi> lookupTable1D;
  float epsilon = 1.e-8;

  // static constexpr
  static constexpr int8_t NMinPhotonRotBkg = 3;

  // Usage when cfgEnableNonLin is enabled
  std::unique_ptr<TF1> fEMCalCorrectionFactor; // ("fEMCalCorrectionFactor","(1 + [0]/x + [1]/x^2) / (1 + [2]/x)", 0.3, 100.);
  float energyCorrectionFactor = 1.f;
  uint8_t nColl = 1;

  // To access the 1D array
  inline int getIndex(int iEta, int iPhi)
  {
    return iEta * NBinsPhi + iPhi;
  }

  // Function to access the lookup table
  inline int8_t checkEtaPhi1D(double eta, double phi)
  {
    if (eta < EtaMin || eta > etaMax || phi < PhiMin || phi > phiMax) {
      return 3; // Out of bounds
    }

    // Compute indices directly
    int iEta = static_cast<int>((eta - EtaMin) / ((etaMax - EtaMin) / NBinsEta));
    int iPhi = static_cast<int>((phi - PhiMin) / ((phiMax - PhiMin) / NBinsPhi));

    return lookupTable1D[getIndex(iEta, iPhi)];
  }

  void defineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetRequireEMCReadoutInMB(eventcuts.cfgRequireEMCReadoutInMB);
    fEMEventCut.SetRequireEMCHardwareTriggered(eventcuts.cfgRequireEMCHardwareTriggered);
  }

  void defineEMCCut()
  {
    fEMCCut = EMCPhotonCut("fEMCCut", "fEMCCut");

    fEMCCut.SetTrackMatchingEtaParams(emccuts.cfgEMCTMEta->at(0), emccuts.cfgEMCTMEta->at(1), emccuts.cfgEMCTMEta->at(2));
    fEMCCut.SetTrackMatchingPhiParams(emccuts.cfgEMCTMPhi->at(0), emccuts.cfgEMCTMPhi->at(1), emccuts.cfgEMCTMPhi->at(2));

    fEMCCut.SetSecTrackMatchingEtaParams(emccuts.emcSecTMEta->at(0), emccuts.emcSecTMEta->at(1), emccuts.emcSecTMEta->at(2));
    fEMCCut.SetSecTrackMatchingPhiParams(emccuts.emcSecTMPhi->at(0), emccuts.emcSecTMPhi->at(1), emccuts.emcSecTMPhi->at(2));
    fEMCCut.SetMinEoverP(emccuts.cfgEMCEoverp);

    fEMCCut.SetMinE(emccuts.cfgEMCminE);
    fEMCCut.SetMinNCell(emccuts.cfgEMCminNCell);
    fEMCCut.SetM02Range(emccuts.cfgEMCminM02, emccuts.cfgEMCmaxM02);
    fEMCCut.SetTimeRange(emccuts.cfgEMCminTime, emccuts.cfgEMCmaxTime);
    fEMCCut.SetUseExoticCut(emccuts.cfgEMCUseExoticCut);
    fEMCCut.SetClusterizer(emccuts.clusterDefinition);
    fEMCCut.SetUseTM(emccuts.cfgEMCUseTM.value);                // disables or enables TM
    fEMCCut.SetUseSecondaryTM(emccuts.emcUseSecondaryTM.value); // disables or enables secondary TM
    fEMCCut.SetDoQA(emccuts.cfgEnableQA.value);
  }

  void definePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    // for v0
    fV0PhotonCut.SetV0PtRange(pcmcuts.minPtV0, pcmcuts.maxPtV0);
    fV0PhotonCut.SetV0EtaRange(pcmcuts.minEtaV0, pcmcuts.maxEtaV0);
    fV0PhotonCut.SetMinCosPA(pcmcuts.minCosPA);
    fV0PhotonCut.SetMaxPCA(pcmcuts.maxPCA);
    fV0PhotonCut.SetMaxChi2KF(pcmcuts.maxChi2KF);
    fV0PhotonCut.SetRxyRange(pcmcuts.minRV0, pcmcuts.maxRV0);
    fV0PhotonCut.SetAPRange(pcmcuts.maxAlphaAP, pcmcuts.maxQtAP);
    fV0PhotonCut.RejectITSib(pcmcuts.rejectV0onITSib);

    // for track
    fV0PhotonCut.SetMinNClustersTPC(pcmcuts.minNClusterTPC);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.minNCrossedRowsTPC);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(pcmcuts.minNCrossedRowsOverFindableClustersTPC);
    fV0PhotonCut.SetMaxFracSharedClustersTPC(pcmcuts.maxFracSharedClustersTPC);
    fV0PhotonCut.SetChi2PerClusterTPC(0.f, pcmcuts.macChi2TPC);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.minTPCNSigmaEl, pcmcuts.maxTPCNSigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(0.f, pcmcuts.macChi2ITS);
    fV0PhotonCut.SetNClustersITS(pcmcuts.minNClusterITS, 7);
    fV0PhotonCut.SetMeanClusterSizeITSob(pcmcuts.minMeanClusterSizeITSob, pcmcuts.maxMeanClusterSizeITSob);
    fV0PhotonCut.SetDisableITSonly(pcmcuts.disableITSOnly);
    fV0PhotonCut.SetDisableTPConly(pcmcuts.disableTPCOnly);
    fV0PhotonCut.SetRequireITSTPC(pcmcuts.requireV0WithITSTPC);
    fV0PhotonCut.SetRequireITSonly(pcmcuts.requireV0WithITSOnly);
    fV0PhotonCut.SetRequireTPConly(pcmcuts.requireV0WithTPCOnly);

    fV0PhotonCut.setDoQA(pcmcuts.doQA.value);
  }

  void init(InitContext&)
  {
    if (harmonic != kElliptic && harmonic != kTriangluar) {
      LOG(info) << "Harmonic was set to " << harmonic << " but can only be 2 or 3!";
    }

    defineEMEventCut();
    defineEMCCut();
    fEMCCut.addQAHistograms(&registry);
    definePCMCut();
    fV0PhotonCut.addQAHistograms(&registry);
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&registry);

    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality (%)"};
    const AxisSpec thnAxisCosDeltaPhi{thnConfigAxisCosDeltaPhi, Form("cos(%d(#varphi - #Psi_{sub}))", harmonic.value)};
    const AxisSpec thnAxisM02{thnConfigAxisM02, "M_{02}"};
    const AxisSpec thAxisTanThetaPhi{mesonConfig.thConfigAxisTanThetaPhi, "atan(#Delta#theta/#Delta#varphi)"};
    const AxisSpec thAxisClusterEnergy{thnConfigAxisPt, "#it{E} (GeV)"};
    const AxisSpec thAxisEnergyCalib{thnConfigAxisEnergyCalib, "#it{E}_{clus} (GeV)"};
    const AxisSpec thAxisAlpha{100, -1., +1, "#alpha"};
    const AxisSpec thAxisEnergy{1000, 0., 100., "#it{E}_{clus} (GeV)"};
    const AxisSpec thAxisEta{320, -0.8, 0.8, "#eta"};
    const AxisSpec thAxisPhi{500, 0, 2 * 3.14159, "phi"};

    const AxisSpec thnAxisMixingVtx{mixingConfig.cfgVtxBins, "#it{z} (cm)"};
    const AxisSpec thnAxisMixingCent{mixingConfig.cfgCentBins, "Centrality (%)"};
    const AxisSpec thnAxisMixingEP{mixingConfig.cfgEPBins, Form("cos(%d#varphi)", harmonic.value)};

    registry.add("hSparsePi0Flow", "<v_n> vs m_{inv} vs p_T vs cent for same event", HistType::kTProfile3D, {thnAxisInvMass, thnAxisPt, thnAxisCent});
    registry.add("hSparsePi0", "m_{inv} vs p_T vs cent for same event", HistType::kTH3D, {thnAxisInvMass, thnAxisPt, thnAxisCent});

    registry.add("hSparseBkgMixFlow", "<v_n> vs m_{inv} vs p_T vs cent for mixed event", HistType::kTProfile3D, {thnAxisInvMass, thnAxisPt, thnAxisCent});
    registry.add("hSparseBkgMix", "m_{inv} vs p_T vs cent for mixed event", HistType::kTH3D, {thnAxisInvMass, thnAxisPt, thnAxisCent});

    if (rotationConfig.cfgDoRotation.value) {
      registry.add("hSparseBkgRotFlow", "<v_n> vs m_{inv} vs p_T vs cent for rotation background", HistType::kTProfile3D, {thnAxisInvMass, thnAxisPt, thnAxisCent});
      registry.add("hSparseBkgRot", "m_{inv} vs p_T vs cent for rotation background", HistType::kTH3D, {thnAxisInvMass, thnAxisPt, thnAxisCent});
    }

    registry.add("h3DMixingCount", "THn Event Mixing QA", HistType::kTH3D, {thnAxisMixingVtx, thnAxisMixingCent, thnAxisMixingEP});
    if (cfgDoPlaneQA.value) {
      registry.add("hSparsePi0FlowPlane", "THn for SP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCent, thnAxisCosDeltaPhi});
    }
    auto hMesonCuts = registry.add<TH1>("hMesonCuts", "hMesonCuts;;Counts", kTH1D, {{6, 0.5, 6.5}}, false);
    hMesonCuts->GetXaxis()->SetBinLabel(1, "in");
    hMesonCuts->GetXaxis()->SetBinLabel(2, "opening angle");
    hMesonCuts->GetXaxis()->SetBinLabel(3, "#it{M}_{#gamma#gamma}");
    hMesonCuts->GetXaxis()->SetBinLabel(4, "#it{p}_{T}");
    hMesonCuts->GetXaxis()->SetBinLabel(5, "conversion cut");
    hMesonCuts->GetXaxis()->SetBinLabel(6, "out");

    auto hMesonCutsMixed = registry.add<TH1>("hMesonCutsMixed", "hMesonCutsMixed;;Counts", kTH1D, {{6, 0.5, 6.5}}, false);
    hMesonCutsMixed->GetXaxis()->SetBinLabel(1, "in");
    hMesonCutsMixed->GetXaxis()->SetBinLabel(2, "opening angle");
    hMesonCutsMixed->GetXaxis()->SetBinLabel(3, "#it{M}_{#gamma#gamma}");
    hMesonCutsMixed->GetXaxis()->SetBinLabel(4, "#it{p}_{T}");
    hMesonCutsMixed->GetXaxis()->SetBinLabel(5, "conversion cut");
    hMesonCutsMixed->GetXaxis()->SetBinLabel(6, "out");

    if (emccuts.cfgEnableQA.value) {
      registry.add("clusterQA/hEClusterBefore", "Histo for cluster energy before cuts", HistType::kTH1D, {thAxisClusterEnergy});
      registry.add("clusterQA/hEClusterAfter", "Histo for cluster energy after cuts", HistType::kTH1D, {thAxisClusterEnergy});
      registry.add("clusterQA/hClusterEtaPhiBefore", "hClusterEtaPhiBefore", HistType::kTH2D, {thAxisPhi, thAxisEta});
      registry.add("clusterQA/hClusterEtaPhiAfter", "hClusterEtaPhiAfter", HistType::kTH2D, {thAxisPhi, thAxisEta});
      if (rotationConfig.cfgDoRotation.value) {
        registry.add("clusterQA/hClusterBackEtaPhiBefore", "hClusterBackEtaPhiBefore", HistType::kTH2D, {thAxisPhi, thAxisEta});
        registry.add("clusterQA/hClusterBackEtaPhiAfter", "hClusterBackEtaPhiAfter", HistType::kTH2D, {thAxisPhi, thAxisEta});
      }
    }

    if (mesonConfig.cfgEnableQA.value) {
      registry.add("mesonQA/hInvMassPt", "Histo for inv pair mass vs pt", HistType::kTH2D, {thnAxisInvMass, thnAxisPt});
      registry.add("mesonQA/hTanThetaPhi", "Histo for identification of conversion cluster", HistType::kTH2D, {thnAxisInvMass, thAxisTanThetaPhi});
      registry.add("mesonQA/hAlphaPt", "Histo of meson asymmetry vs pT", HistType::kTH2D, {thAxisAlpha, thnAxisPt});
      registry.add("mesonQA/hInvMassPtMixed", "Histo for inv pair mass vs pt for mixed event", HistType::kTH2D, {thnAxisInvMass, thnAxisPt});
      registry.add("mesonQA/hTanThetaPhiMixed", "Histo for identification of conversion cluster for mixed event", HistType::kTH2D, {thnAxisInvMass, thAxisTanThetaPhi});
      registry.add("mesonQA/hAlphaPtMixed", "Histo of meson asymmetry vs pT for mixed event", HistType::kTH2D, {thAxisAlpha, thnAxisPt});
    }

    if (correctionConfig.doEMCalCalib.value) {
      registry.add("hSparseCalibSE", "THn for Calib same event", HistType::kTHnSparseF, {thnAxisInvMass, thAxisEnergyCalib, thnAxisCent});
      registry.add("hSparseCalibBack", "THn for Calib background", HistType::kTHnSparseF, {thnAxisInvMass, thAxisEnergyCalib, thnAxisCent});
    }

    if (cfgDoM02.value) {
      registry.add("p3DM02Flow", "<v_n> vs M_{02} vs p_T vs cent", HistType::kTProfile3D, {thnAxisM02, thnAxisPt, thnAxisCent});
      registry.add("h3DSparsePi0", "M_{02} vs p_T vs cent", HistType::kTH3D, {thnAxisM02, thnAxisPt, thnAxisCent});
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    LOG(info) << "thnConfigAxisInvMass.value[1] = " << thnConfigAxisInvMass.value[1] << " thnConfigAxisInvMass.value.back() = " << thnConfigAxisInvMass.value.back();
    LOG(info) << "thnConfigAxisPt.value[1] = " << thnConfigAxisPt.value[1] << " thnConfigAxisPt.value.back() = " << thnConfigAxisPt.value.back();

    fEMCalCorrectionFactor = std::make_unique<TF1>("fEMCalCorrectionFactor", "(1 + [0]/x + [1]/x^2) / (1 + [2]/x)", 0.3, 100.);
    fEMCalCorrectionFactor->SetParameters(-5.33426e-01, 1.40144e-02, -5.24434e-01);
  }; // end init

  /// Change radians to degree
  /// \param angle in radians
  /// \return angle in degree
  constexpr float getAngleDegree(float angle)
  {
    return angle * o2::constants::math::Rad2Deg;
  }

  /// Compute the delta psi in the range [0, pi/harmonic]
  /// \param psi1 is the first angle
  /// \param psi2 is the second angle
  float getDeltaPsiInRange(float psi1, float psi2)
  {
    float deltaPsi = psi1 - psi2;
    return RecoDecay::constrainAngle(deltaPsi, 0.f, harmonic);
  }

  /// Fill THnSparse
  /// \param mass is the invariant mass of the candidate
  /// \param pt is the transverse momentum of the candidate
  /// \param cent is the centrality of the collision
  /// \param sp is the scalar product
  template <const int histType>
  void fillThn(const float mass, const float pt, const float cent, const float sp)
  {
    static constexpr std::string_view FlowHistTypes[3] = {"hSparsePi0Flow", "hSparseBkgRotFlow", "hSparseBkgMixFlow"};
    static constexpr std::string_view HistTypes[3] = {"hSparsePi0", "hSparseBkgRot", "hSparseBkgMix"};
    registry.fill(HIST(FlowHistTypes[histType]), mass, pt, cent, sp);
    registry.fill(HIST(HistTypes[histType]), mass, pt, cent);
  }

  /// Get the centrality
  /// \param collision is the collision with the centrality information
  template <o2::soa::is_iterator TCollision>
  float getCentrality(TCollision const& collision)
  {
    float cent = -999.;
    switch (centEstimator) {
      case CentralityEstimator::CFT0M:
        cent = collision.centFT0M();
        break;
      case CentralityEstimator::CFT0A:
        cent = collision.centFT0A();
        break;
      case CentralityEstimator::CFT0C:
        cent = collision.centFT0C();
        break;
      default:
        LOG(warning) << "Centrality estimator not valid. Possible values are T0M, T0A, T0C. Fallback to T0C";
        cent = collision.centFT0C();
        break;
    }
    return cent;
  }

  /// Get all used Q vector
  /// \param collision is the collision with the Q vector information
  template <o2::soa::is_iterator TCollision>
  std::array<float, NQVecEntries> getAllQvec(TCollision const& collision)
  {
    // Retrieve the Q vectors using the helper function for each detector
    auto [xQVecMain, yQVecMain] = getQvec(collision, qvecDetector);
    auto [xQVecSubA, yQVecSubA] = getQvec(collision, qvecSubADetector);
    auto [xQVecSubB, yQVecSubB] = getQvec(collision, qvecSubBDetector);

    return {xQVecMain, yQVecMain, xQVecSubA, yQVecSubA, xQVecSubB, yQVecSubB};
  }

  /// Get the Q vector
  /// \param collision is the collision with the Q vector information
  template <o2::soa::is_iterator TCollision>
  std::pair<float, float> getQvec(TCollision const& collision, int detector)
  {
    float xQVec = -999.f;
    float yQVec = -999.f;

    switch (detector) {
      case QvecEstimator::FT0M:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xft0m();
          yQVec = collision.q2yft0m();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xft0m();
          yQVec = collision.q3yft0m();
        }
        break;
      case QvecEstimator::FT0A:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xft0a();
          yQVec = collision.q2yft0a();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xft0a();
          yQVec = collision.q3yft0a();
        }
        break;
      case QvecEstimator::FT0C:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xft0c();
          yQVec = collision.q2yft0c();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xft0c();
          yQVec = collision.q3yft0c();
        }
        break;
      case QvecEstimator::TPCPos:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xbpos();
          yQVec = collision.q2ybpos();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xbpos();
          yQVec = collision.q3ybpos();
        }
        break;
      case QvecEstimator::TPCNeg:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xbneg();
          yQVec = collision.q2ybneg();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xbneg();
          yQVec = collision.q3ybneg();
        }
        break;
      case QvecEstimator::TPCTot:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xbtot();
          yQVec = collision.q2ybtot();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xbtot();
          yQVec = collision.q3ybtot();
        }
        break;
      case QvecEstimator::FV0A:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xfv0a();
          yQVec = collision.q2yfv0a();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xfv0a();
          yQVec = collision.q3yfv0a();
        }
        break;
      default:
        LOG(warning) << "Q vector estimator not valid. Falling back to FT0M";
        if (harmonic == kElliptic) {
          xQVec = collision.q2xft0m();
          yQVec = collision.q2yft0m();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xft0m();
          yQVec = collision.q3yft0m();
        }
        break;
    }
    return {xQVec, yQVec};
  }

  /// Check if the QVector values are within reasonable range
  /// \param collision is the collision with the Q vector information
  bool isQvecGood(std::array<float, NQVecEntries> const& QVecs)
  {
    bool isgood = true;
    for (const auto& QVec : QVecs) {
      if (std::fabs(QVec) > cfgMaxQVector) {
        isgood = false;
        break;
      }
    }
    return isgood;
  }

  bool isTooCloseToEdge(const int cellID, const int DistanceToBorder = 1)
  {
    if (DistanceToBorder <= 0) {
      return false;
    }
    if (cellID < 0) {
      return true;
    }

    int iBadCell = -1;

    // check distance to border in case the cell is okay
    auto [iSupMod, iMod, iPhi, iEta] = emcalGeom->GetCellIndex(cellID);
    auto [irow, icol] = emcalGeom->GetCellPhiEtaIndexInSModule(iSupMod, iMod, iPhi, iEta);

    // Check rows/phi
    int iRowLast = 24;
    if (emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::EMCAL_HALF) {
      iRowLast /= 2; // 2/3 sm case
    } else if (emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::EMCAL_THIRD) {
      iRowLast /= 3; // 1/3 sm case
    } else if (emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::DCAL_EXT) {
      iRowLast /= 3; // 1/3 sm case
    }

    if (irow < DistanceToBorder || (iRowLast - irow) <= DistanceToBorder) {
      iBadCell = 1;
    }

    if (iBadCell > 0) {
      return true;
    }
    return false;
  }

  bool isCellMasked(int cellID)
  {
    bool masked = false;
    if (mBadChannels) {
      auto maskStatus = mBadChannels->getChannelStatus(cellID);
      masked = (maskStatus != o2::emcal::BadChannelMap::MaskType_t::GOOD_CELL);
    }
    return masked;
  }

  template <o2::soa::is_iterator TCollision>
  void initCCDB(TCollision const& collision)
  {
    // Load EMCal geometry
    emcalGeom = o2::emcal::Geometry::GetInstanceFromRunNumber(collision.runNumber());
    // Load Bad Channel map
    mBadChannels = ccdb->getForTimeStamp<o2::emcal::BadChannelMap>("EMC/Calib/BadChannelMap", collision.timestamp());
    lookupTable1D.fill(-1);
    double binWidthEta = (etaMax - EtaMin) / NBinsEta;
    double binWidthPhi = (phiMax - PhiMin) / NBinsPhi;

    if (cfgEMCalMapLevelBackground.value >= static_cast<int>(MapLevel::kAll) && cfgEMCalMapLevelSameEvent >= static_cast<int>(MapLevel::kAll)) {
      // in this case we do not want to check the clusters, so just say thery are all good.
      lookupTable1D.fill(0); // good
    } else {
      for (int iEta = 0; iEta < NBinsEta; ++iEta) {
        double etaCenter = EtaMin + (iEta + 0.5) * binWidthEta;
        for (int iPhi = 0; iPhi < NBinsPhi; ++iPhi) {
          double phiCenter = PhiMin + (iPhi + 0.5) * binWidthPhi;
          try {
            // Get the cell ID
            int cellID = emcalGeom->GetAbsCellIdFromEtaPhi(etaCenter, phiCenter);

            // Check conditions for the cell
            if (isTooCloseToEdge(cellID, 1)) {
              lookupTable1D[getIndex(iEta, iPhi)] = 2; // Edge
            } else if (isCellMasked(cellID)) {
              lookupTable1D[getIndex(iEta, iPhi)] = 1; // Bad
            } else {
              lookupTable1D[getIndex(iEta, iPhi)] = 0; // Good
            }
          } catch (o2::emcal::InvalidPositionException& e) {
            lookupTable1D[getIndex(iEta, iPhi)] = 3; // Outside geometry
          }
        }
      }
    }
    if (correctionConfig.cfgApplySPresolution.value) {
      h1SPResolution = ccdb->getForTimeStamp<TH1D>(correctionConfig.cfgSpresoPath.value, collision.timestamp());
    }
  }

  /// \brief Calculate background using rotation background method
  template <typename TPhotons, o2::soa::is_iterator TCollision>
  void rotationBackground(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photon1, ROOT::Math::PtEtaPhiMVector photon2, TPhotons const& photonsColl, unsigned int ig1, unsigned int ig2, TCollision const& collision)
  {
    // if less than 3 clusters are present skip event since we need at least 3 clusters
    if (photonsColl.size() < NMinPhotonRotBkg) {
      return;
    }
    float weight = 1.f;
    if (rotationConfig.cfgUseWeights) {
      weight = 1.f / (2.f * static_cast<float>(photonsColl.size() - 2));
    }

    auto [xQVec, yQVec] = getQvec(collision, qvecDetector);
    float cent = getCentrality(collision);
    int iCellIDPhoton1 = 0;
    int iCellIDPhoton2 = 0;

    ROOT::Math::AxisAngle rotationAxis(meson.Vect(), rotationConfig.cfgRotAngle);
    ROOT::Math::Rotation3D rotationMatrix(rotationAxis);
    photon1 = rotationMatrix * photon1;
    photon2 = rotationMatrix * photon2;

    if (emccuts.cfgEnableQA.value) {
      registry.fill(HIST("clusterQA/hClusterBackEtaPhiBefore"), RecoDecay::constrainAngle(photon1.Phi()), photon1.Eta()); // before check but after rotation
      registry.fill(HIST("clusterQA/hClusterBackEtaPhiBefore"), RecoDecay::constrainAngle(photon2.Phi()), photon2.Eta()); // before check but after rotation
    }

    if (checkEtaPhi1D(photon1.Eta(), RecoDecay::constrainAngle(photon1.Phi())) >= cfgEMCalMapLevelBackground.value) {
      iCellIDPhoton1 = -1;
    } else if (emccuts.cfgEnableQA.value) {
      registry.fill(HIST("clusterQA/hClusterBackEtaPhiAfter"), RecoDecay::constrainAngle(photon1.Phi()), photon1.Eta()); // after check
    }
    if (checkEtaPhi1D(photon2.Eta(), RecoDecay::constrainAngle(photon2.Phi())) >= cfgEMCalMapLevelBackground.value) {
      iCellIDPhoton2 = -1;
    } else if (emccuts.cfgEnableQA.value) {
      registry.fill(HIST("clusterQA/hClusterBackEtaPhiAfter"), RecoDecay::constrainAngle(photon2.Phi()), photon2.Eta()); // after check
    }
    if (iCellIDPhoton1 == -1 && iCellIDPhoton2 == -1) {
      return;
    }
    for (const auto& photon : photonsColl) {
      if (photon.globalIndex() == ig1 || photon.globalIndex() == ig2) {
        // only combine rotated photons with other photons
        continue;
      }
      if (!(fEMCCut.IsSelected(photon))) {
        continue;
      }
      if (checkEtaPhi1D(photon.eta(), RecoDecay::constrainAngle(photon.phi())) >= cfgEMCalMapLevelBackground.value) {
        continue;
      }
      energyCorrectionFactor = 1.f;
      if (correctionConfig.cfgEnableNonLin.value) {
        energyCorrectionFactor = fEMCalCorrectionFactor->Eval(photon.e() > MinEnergy ? photon.e() : MinEnergy);
      }
      ROOT::Math::PtEtaPhiMVector photon3(energyCorrectionFactor * photon.pt(), photon.eta(), photon.phi(), 0.);
      if (iCellIDPhoton1 >= 0) {
        ROOT::Math::PtEtaPhiMVector mother1 = photon1 + photon3;
        float openingAngle1 = std::acos(photon1.Vect().Dot(photon3.Vect()) / (photon1.P() * photon3.P()));
        float cosNPhi1 = std::cos(harmonic * mother1.Phi());
        float sinNPhi1 = std::sin(harmonic * mother1.Phi());
        float scalprodCand1 = cosNPhi1 * xQVec + sinNPhi1 * yQVec;

        if (correctionConfig.cfgApplySPresolution.value) {
          scalprodCand1 = scalprodCand1 / h1SPResolution->GetBinContent(h1SPResolution->FindBin(cent + epsilon));
        }

        if (openingAngle1 > mesonConfig.minOpenAngle && thnConfigAxisInvMass.value[1] <= mother1.M() && thnConfigAxisInvMass.value.back() >= mother1.M() && thnConfigAxisPt.value[1] <= mother1.Pt() && thnConfigAxisPt.value.back() >= mother1.Pt()) {
          if (mesonConfig.enableTanThetadPhi) {
            float dTheta = photon1.Theta() - photon3.Theta();
            float dPhi = photon1.Phi() - photon3.Phi();
            if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi <= std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
              continue;
            }
          }
          registry.fill(HIST("hSparseBkgRotFlow"), mother1.M(), mother1.Pt(), cent, scalprodCand1, weight);
          registry.fill(HIST("hSparseBkgRot"), mother1.M(), mother1.Pt(), cent, weight);
        }
        if (iCellIDPhoton2 >= 0) {
          ROOT::Math::PtEtaPhiMVector mother2 = photon2 + photon3;
          float openingAngle2 = std::acos(photon2.Vect().Dot(photon3.Vect()) / (photon2.P() * photon3.P()));
          float cosNPhi2 = std::cos(harmonic * mother2.Phi());
          float sinNPhi2 = std::sin(harmonic * mother2.Phi());
          float scalprodCand2 = cosNPhi2 * xQVec + sinNPhi2 * yQVec;

          if (correctionConfig.cfgApplySPresolution.value) {
            scalprodCand2 = scalprodCand2 / h1SPResolution->GetBinContent(h1SPResolution->FindBin(cent + epsilon));
          }

          if (openingAngle2 > mesonConfig.minOpenAngle && thnConfigAxisInvMass.value[1] <= mother2.M() && thnConfigAxisInvMass.value.back() >= mother2.M() && thnConfigAxisPt.value[1] <= mother2.Pt() && thnConfigAxisPt.value.back() >= mother2.Pt()) {
            if (mesonConfig.enableTanThetadPhi) {
              float dTheta = photon2.Theta() - photon3.Theta();
              float dPhi = photon2.Phi() - photon3.Phi();
              if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi <= std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
                continue;
              }
            }
            registry.fill(HIST("hSparseBkgRotFlow"), mother2.M(), mother2.Pt(), cent, scalprodCand2, weight);
            registry.fill(HIST("hSparseBkgRot"), mother2.M(), mother2.Pt(), cent, weight);
          }
        } // end of loop over third photon
      }
    }
    return;
  }

  /// \brief Calculate background using rotation background method
  /// \param meson meson candidate whos momentum is used as rotation axis
  /// \param photonPCM first photon candidate (PCM)
  /// \param photonEMC second photon candidate (EMC)
  /// \param photonsCollPCM subtable of current collisions PCM photons
  /// \param photonsCollEMC subtable of current collisions EMCal photons
  /// \param giPCM global index of photonPCM
  /// \param giEMC global index of photonEMC
  /// \param collision current collision iterator
  template <typename TPCMPhotons, typename TEMCPhotons, o2::soa::is_iterator TCollision>
  void rotationBackgroundPCMEMC(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photonPCM, ROOT::Math::PtEtaPhiMVector photonEMC, TPCMPhotons const& photonsCollPCM, TEMCPhotons const& photonsCollEMC, unsigned int giPCM, unsigned int giEMC, TCollision const& collision)
  {
    // if less than 3 photon candidates are present skip event since we need at least 3 candidates
    if ((photonsCollEMC.size() + photonsCollPCM.size()) < NMinPhotonRotBkg) {
      return;
    }
    float weight = 1.f;
    if (rotationConfig.cfgUseWeights) {
      weight = 1.f / static_cast<float>(photonsCollPCM.size() + photonsCollEMC.size() - 2);
    }

    auto [xQVec, yQVec] = getQvec(collision, qvecDetector);
    float cent = getCentrality(collision);
    int iCellIDphotonPCM = 0;
    int iCellIDphotonEMC = 0;

    ROOT::Math::AxisAngle rotationAxis(meson.Vect(), rotationConfig.cfgRotAngle);
    ROOT::Math::Rotation3D rotationMatrix(rotationAxis);
    photonPCM = rotationMatrix * photonPCM;
    photonEMC = rotationMatrix * photonEMC;

    if (emccuts.cfgEnableQA.value) {
      registry.fill(HIST("clusterQA/hClusterBackEtaPhiBefore"), RecoDecay::constrainAngle(photonEMC.Phi()), photonEMC.Eta()); // before check but after rotation
    }

    if (photonPCM.Eta() < pcmcuts.minEtaV0 || photonPCM.Eta() > pcmcuts.maxEtaV0) {
      iCellIDphotonPCM = -1;
    }
    if (checkEtaPhi1D(photonEMC.Eta(), RecoDecay::constrainAngle(photonEMC.Phi())) >= cfgEMCalMapLevelBackground.value) {
      iCellIDphotonEMC = -1;
    } else if (emccuts.cfgEnableQA.value) {
      registry.fill(HIST("clusterQA/hClusterBackEtaPhiAfter"), RecoDecay::constrainAngle(photonEMC.Phi()), photonEMC.Eta()); // after check
    }
    if (iCellIDphotonPCM == -1 && iCellIDphotonEMC == -1) {
      return;
    }
    if (iCellIDphotonPCM > -1) {
      for (const auto& photon : photonsCollEMC) {
        if (photon.globalIndex() == giEMC) {
          // only combine rotated photons with other photons
          continue;
        }
        if (!(fEMCCut.IsSelected(photon))) {
          continue;
        }
        if (checkEtaPhi1D(photon.eta(), RecoDecay::constrainAngle(photon.phi())) >= cfgEMCalMapLevelBackground.value) {
          continue;
        }
        energyCorrectionFactor = 1.f;
        if (correctionConfig.cfgEnableNonLin.value) {
          energyCorrectionFactor = fEMCalCorrectionFactor->Eval(photon.e() > MinEnergy ? photon.e() : MinEnergy);
        }
        ROOT::Math::PtEtaPhiMVector photon3(energyCorrectionFactor * photon.pt(), photon.eta(), photon.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector mother1 = photonPCM + photon3;
        float openingAngle = std::acos(photonPCM.Vect().Dot(photon3.Vect()) / (photonPCM.P() * photon3.P()));
        float cosNPhi = std::cos(harmonic * mother1.Phi());
        float sinNPhi = std::sin(harmonic * mother1.Phi());
        float scalprodCand = cosNPhi * xQVec + sinNPhi * yQVec;

        if (correctionConfig.cfgApplySPresolution.value) {
          scalprodCand = scalprodCand / h1SPResolution->GetBinContent(h1SPResolution->FindBin(cent + epsilon));
        }
        if (openingAngle > mesonConfig.minOpenAngle && thnConfigAxisInvMass.value[1] <= mother1.M() && thnConfigAxisInvMass.value.back() >= mother1.M() && thnConfigAxisPt.value[1] <= mother1.Pt() && thnConfigAxisPt.value.back() >= mother1.Pt()) {
          if (mesonConfig.enableTanThetadPhi) {
            float dTheta = photonPCM.Theta() - photon3.Theta();
            float dPhi = photonPCM.Phi() - photon3.Phi();
            if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi <= std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
              continue;
            }
          }
          registry.fill(HIST("hSparseBkgRotFlow"), mother1.M(), mother1.Pt(), cent, scalprodCand, weight);
          registry.fill(HIST("hSparseBkgRot"), mother1.M(), mother1.Pt(), cent, weight);
        }
      } // end of loop over EMC photons
    } // if(iCellIDphotonPCM > -1)
    if (iCellIDphotonEMC > -1) {
      for (const auto& photon : photonsCollPCM) {
        if (photon.globalIndex() == giPCM) {
          // only combine rotated photons with other photons
          continue;
        }
        if (!(fV0PhotonCut.IsSelected<decltype(photon), aod::V0Legs>(photon))) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector photon3(photon.pt(), photon.eta(), photon.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector mother1 = photonEMC + photon3;
        float openingAngle = std::acos(photonEMC.Vect().Dot(photon3.Vect()) / (photonEMC.P() * photon3.P()));
        float cosNPhi = std::cos(harmonic * mother1.Phi());
        float sinNPhi = std::sin(harmonic * mother1.Phi());
        float scalprodCand = cosNPhi * xQVec + sinNPhi * yQVec;

        if (correctionConfig.cfgApplySPresolution.value) {
          scalprodCand = scalprodCand / h1SPResolution->GetBinContent(h1SPResolution->FindBin(cent + epsilon));
        }
        if (openingAngle > mesonConfig.minOpenAngle && thnConfigAxisInvMass.value[1] <= mother1.M() && thnConfigAxisInvMass.value.back() >= mother1.M() && thnConfigAxisPt.value[1] <= mother1.Pt() && thnConfigAxisPt.value.back() >= mother1.Pt()) {
          if (mesonConfig.enableTanThetadPhi) {
            float dTheta = photonEMC.Theta() - photon3.Theta();
            float dPhi = photonEMC.Phi() - photon3.Phi();
            if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi <= std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
              continue;
            }
          }
          registry.fill(HIST("hSparseBkgRotFlow"), mother1.M(), mother1.Pt(), cent, scalprodCand, weight);
          registry.fill(HIST("hSparseBkgRot"), mother1.M(), mother1.Pt(), cent, weight);
        }
      } // end of loop over PCM photons
    } // if(iCellIDphotonEMC > -1)
    return;
  }

  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param meson are the selected candidates
  template <const int histType, o2::soa::is_iterator TCollision>
  void runFlowAnalysis(TCollision const& collision, ROOT::Math::PtEtaPhiMVector const& meson)
  {
    auto [xQVec, yQVec] = getQvec(collision, qvecDetector);
    float cent = getCentrality(collision);

    float massCand = meson.M();
    float ptCand = meson.Pt();
    float phiCand = meson.Phi();

    float cosNPhi = std::cos(harmonic * phiCand);
    float sinNPhi = std::sin(harmonic * phiCand);
    float scalprodCand = cosNPhi * xQVec + sinNPhi * yQVec;

    if (correctionConfig.cfgApplySPresolution.value) {
      scalprodCand = scalprodCand / h1SPResolution->GetBinContent(h1SPResolution->FindBin(cent + epsilon));
    }

    if (cfgDoPlaneQA.value && histType == 0) {
      float epAngle = epHelper.GetEventPlane(xQVec, yQVec, harmonic);
      float cosDeltaPhi = std::cos(harmonic * getDeltaPsiInRange(phiCand, epAngle));
      registry.fill(HIST("hSparsePi0FlowPlane"), massCand, ptCand, cent, cosDeltaPhi);
    }

    fillThn<histType>(massCand, ptCand, cent, scalprodCand);
    return;
  }

  /// \brief check if standard event cuts + FT0 occupancy + centrality + QVec good is
  /// \param collision collision that will be checked
  /// \return true if collision survives all checks, otherwise false
  template <typename TCollision>
  bool isFullEventSelected(TCollision const& collision, bool fillHisto = false)
  {
    if (fillHisto) {
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&registry, collision);
    }
    if (!(fEMEventCut.IsSelected(collision))) {
      // general event selection
      return false;
    }
    if (!(eventcuts.cfgFT0COccupancyMin <= collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax)) {
      // occupancy selection
      return false;
    }
    float cent = getCentrality(collision);
    if (cent < eventcuts.cfgMinCent || cent > eventcuts.cfgMaxCent) {
      // event selection
      return false;
    }
    if (!isQvecGood(getAllQvec(collision))) {
      // selection based on QVector
      return false;
    }
    if (fillHisto) {
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&registry, collision);
      registry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // accepted
      registry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // accepted
    }
    return true;
  }

  template <typename TCollision, typename TPhotons1, typename TPhotons2>
  void runPairingLoop(TCollision const& collision, TPhotons1 const& photons1, TPhotons2 const& photons2, EMBitFlags const& flags1, EMBitFlags const& flags2)
  {
    for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1, photons2))) {
      if (!(flags1.test(g1.globalIndex())) || !(flags2.test(g2.globalIndex()))) {
        continue;
      }

      // Cut edge clusters away, similar to rotation method to ensure same acceptance is used
      if (cfgDistanceToEdge.value) {
        if (checkEtaPhi1D(g1.eta(), RecoDecay::constrainAngle(g1.phi())) >= cfgEMCalMapLevelSameEvent.value) {
          continue;
        }
        if (checkEtaPhi1D(g2.eta(), RecoDecay::constrainAngle(g2.phi())) >= cfgEMCalMapLevelSameEvent.value) {
          continue;
        }
      }
      if (correctionConfig.cfgEnableNonLin.value) {
        energyCorrectionFactor = fEMCalCorrectionFactor->Eval(g1.e() > MinEnergy ? g1.e() : MinEnergy);
      }
      if (correctionConfig.cfgEnableNonLin.value) {
        energyCorrectionFactor = fEMCalCorrectionFactor->Eval(g2.e() > MinEnergy ? g2.e() : MinEnergy);
      }

      ROOT::Math::PtEtaPhiMVector v1(energyCorrectionFactor * g1.pt(), g1.eta(), g1.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector v2(energyCorrectionFactor * g2.pt(), g2.eta(), g2.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;

      float dTheta = v1.Theta() - v2.Theta();
      float dPhi = v1.Phi() - v2.Phi();
      float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));

      registry.fill(HIST("hMesonCuts"), 1);
      if (openingAngle <= mesonConfig.minOpenAngle) {
        registry.fill(HIST("hMesonCuts"), 2);
        continue;
      }
      if (rotationConfig.cfgDoRotation.value && nColl % rotationConfig.cfgDownsampling == 0) {
        rotationBackground<EMCalPhotons>(vMeson, v1, v2, photons1, g1.globalIndex(), g2.globalIndex(), collision);
      }
      if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
        registry.fill(HIST("hMesonCuts"), 3);
        continue;
      }
      if (thnConfigAxisPt.value[1] > vMeson.Pt() || thnConfigAxisPt.value.back() < vMeson.Pt()) {
        registry.fill(HIST("hMesonCuts"), 4);
        continue;
      }
      if (mesonConfig.cfgEnableQA.value) {
        registry.fill(HIST("mesonQA/hInvMassPt"), vMeson.M(), vMeson.Pt());
        registry.fill(HIST("mesonQA/hTanThetaPhi"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
        registry.fill(HIST("mesonQA/hAlphaPt"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
      }
      if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
        registry.fill(HIST("hMesonCuts"), 5);
        continue;
      }
      registry.fill(HIST("hMesonCuts"), 6);
      runFlowAnalysis<0>(collision, vMeson);
    }
  }

  // Pi0 from EMCal
  void processEMCal(CollsWithQvecs const& collisions, EMCalPhotons const& clusters, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds)
  {
    if (clusters.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    EMBitFlags flags(clusters.size());
    fEMCCut.AreSelectedRunning(flags, clusters, matchedPrims, matchedSeconds, &registry);

    energyCorrectionFactor = 1.f;
    if (cfgDoReverseScaling.value) {
      energyCorrectionFactor = 1.0505f;
    }
    for (const auto& collision : collisions) {

      if (!isFullEventSelected(collision, true)) {
        continue;
      }
      runNow = collision.runNumber();
      if (runNow != runBefore) {
        initCCDB(collision);
        runBefore = runNow;
      }

      auto photonsPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());

      if (emccuts.cfgEnableQA.value) {
        for (const auto& photon : photonsPerCollision) {
          registry.fill(HIST("clusterQA/hEClusterBefore"), photon.e());                      // before cuts
          registry.fill(HIST("clusterQA/hClusterEtaPhiBefore"), photon.phi(), photon.eta()); // before cuts
          if (!(fEMCCut.IsSelected(photon))) {
            continue;
          }
          if (cfgDistanceToEdge.value && (checkEtaPhi1D(photon.eta(), RecoDecay::constrainAngle(photon.phi())) >= cfgEMCalMapLevelSameEvent.value)) {
            continue;
          }
          registry.fill(HIST("clusterQA/hEClusterAfter"), photon.e());                      // accepted after cuts
          registry.fill(HIST("clusterQA/hClusterEtaPhiAfter"), photon.phi(), photon.eta()); // after cuts
        }
      }
      runPairingLoop(collision, photonsPerCollision, photonsPerCollision, flags, flags);
      if (rotationConfig.cfgDoRotation.value) {
        if (nColl % rotationConfig.cfgDownsampling == 0) {
          nColl = 1; // reset counter
        } else {
          nColl++;
        }
      }
    } // end of loop over collisions
  }
  PROCESS_SWITCH(TaskPi0FlowEMC, processEMCal, "Process EMCal Pi0 candidates", true);

  // Pi0 from EMCal
  void processEMCalMixed(FilteredCollsWithQvecs const& collisions, EMCalPhotons const& clusters, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds)
  {
    if (clusters.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    EMBitFlags flags(clusters.size());
    fEMCCut.AreSelectedRunning(flags, clusters, matchedPrims, matchedSeconds);
    energyCorrectionFactor = 1.f;
    if (cfgDoReverseScaling.value) {
      energyCorrectionFactor = 1.0505f;
    }

    using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, emevent::EP2FT0C<emevent::Q2xFT0C, emevent::Q2yFT0C>>;
    BinningType binningMixedEvent{{mixingConfig.cfgVtxBins, mixingConfig.cfgCentBins, mixingConfig.cfgEPBins}, true};

    auto clustersTuple = std::make_tuple(clusters);
    SameKindPair<FilteredCollsWithQvecs, EMCalPhotons, BinningType> pair{binningMixedEvent, mixingConfig.cfgMixingDepth, -1, collisions, clustersTuple, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

    for (const auto& [c1, clusters1, c2, clusters2] : pair) {
      if (!(fEMEventCut.IsSelected(c1)) || !(fEMEventCut.IsSelected(c2))) {
        // general event selection
        continue;
      }
      if (!(eventcuts.cfgFT0COccupancyMin <= c1.ft0cOccupancyInTimeRange() && c1.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax) || !(eventcuts.cfgFT0COccupancyMin <= c2.ft0cOccupancyInTimeRange() && c2.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax)) {
        // occupancy selection
        continue;
      }
      if (getCentrality(c1) < eventcuts.cfgMinCent || getCentrality(c1) > eventcuts.cfgMaxCent || getCentrality(c2) < eventcuts.cfgMinCent || getCentrality(c2) > eventcuts.cfgMaxCent) {
        // event selection
        continue;
      }
      if (!isQvecGood(getAllQvec(c1)) || !isQvecGood(getAllQvec(c2))) {
        // selection based on QVector
        continue;
      }
      runNow = c1.runNumber();
      if (runNow != runBefore) {
        initCCDB(c1);
        runBefore = runNow;
      }
      registry.fill(HIST("h3DMixingCount"), c1.posZ(), getCentrality(c1), c1.ep2ft0m());
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(clusters1, clusters2))) {
        if (!(flags.test(g1.globalIndex())) || !(flags.test(g2.globalIndex()))) {
          continue;
        }
        // Cut edge clusters away, similar to rotation method to ensure same acceptance is used
        if (cfgDistanceToEdge.value) {
          if (checkEtaPhi1D(g1.eta(), RecoDecay::constrainAngle(g1.phi())) >= cfgEMCalMapLevelBackground.value) {
            continue;
          }
          if (checkEtaPhi1D(g2.eta(), RecoDecay::constrainAngle(g2.phi())) >= cfgEMCalMapLevelBackground.value) {
            continue;
          }
        }
        if (correctionConfig.cfgEnableNonLin.value) {
          energyCorrectionFactor = fEMCalCorrectionFactor->Eval(g1.e() > MinEnergy ? g1.e() : MinEnergy);
        }
        ROOT::Math::PtEtaPhiMVector v1(energyCorrectionFactor * g1.pt(), g1.eta(), g1.phi(), 0.);
        if (correctionConfig.cfgEnableNonLin.value) {
          energyCorrectionFactor = fEMCalCorrectionFactor->Eval(g2.e() > MinEnergy ? g2.e() : MinEnergy);
        }
        ROOT::Math::PtEtaPhiMVector v2(energyCorrectionFactor * g2.pt(), g2.eta(), g2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;

        float dTheta = v1.Theta() - v2.Theta();
        float dPhi = v1.Phi() - v2.Phi();
        float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));

        registry.fill(HIST("hMesonCutsMixed"), 1);
        if (openingAngle <= mesonConfig.minOpenAngle) {
          registry.fill(HIST("hMesonCutsMixed"), 2);
          continue;
        }
        if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
          registry.fill(HIST("hMesonCutsMixed"), 3);
          continue;
        }
        if (thnConfigAxisPt.value[1] > vMeson.Pt() || thnConfigAxisPt.value.back() < vMeson.Pt()) {
          registry.fill(HIST("hMesonCutsMixed"), 4);
          continue;
        }
        if (mesonConfig.cfgEnableQA.value) {
          registry.fill(HIST("mesonQA/hInvMassPtMixed"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("mesonQA/hTanThetaPhiMixed"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
          registry.fill(HIST("mesonQA/hAlphaPtMixed"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
          registry.fill(HIST("hMesonCutsMixed"), 5);
          continue;
        }
        registry.fill(HIST("hMesonCutsMixed"), 6);
        runFlowAnalysis<2>(c1, vMeson);
      }
    }
  }
  PROCESS_SWITCH(TaskPi0FlowEMC, processEMCalMixed, "Process EMCal Pi0 mixed event candidates", false);

  // PCM-EMCal same event
  void processEMCalPCMC(CollsWithQvecs const& collisions, EMCalPhotons const& clusters, PCMPhotons const& photons, aod::V0Legs const&, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds)
  {
    if (clusters.size() <= 0 && photons.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    EMBitFlags emcFlags(clusters.size());
    if (clusters.size() > 0) {
      fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds, &registry);
    }

    EMBitFlags v0flags(photons.size());
    if (photons.size() > 0) {
      fV0PhotonCut.AreSelectedRunning<decltype(photons), aod::V0Legs>(v0flags, photons, &registry);
    }

    energyCorrectionFactor = 1.f;
    if (cfgDoReverseScaling.value) {
      energyCorrectionFactor = 1.0505f;
    }
    for (const auto& collision : collisions) {

      if (!isFullEventSelected(collision, true)) {
        continue;
      }
      runNow = collision.runNumber();
      if (runNow != runBefore) {
        initCCDB(collision);
        runBefore = runNow;
      }

      auto photonsEMCPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());
      auto photonsPCMPerCollision = photons.sliceBy(perCollisionPCM, collision.globalIndex());

      if (emccuts.cfgEnableQA) {
        for (const auto& photon : photonsEMCPerCollision) {
          registry.fill(HIST("clusterQA/hEClusterBefore"), photon.e()); // before cuts
          if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(photon))) {
            continue;
          }
          registry.fill(HIST("clusterQA/hEClusterAfter"), photon.e()); // accepted after cuts
        }
      }
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photonsEMCPerCollision, photonsPCMPerCollision))) {
        if (!(emcFlags.test(g1.globalIndex())) || !(v0flags.test(g2.globalIndex()))) {
          continue;
        }

        // Cut edge clusters away, similar to rotation method to ensure same acceptance is used
        if (cfgDistanceToEdge.value) {
          if (checkEtaPhi1D(g1.eta(), RecoDecay::constrainAngle(g1.phi())) >= cfgEMCalMapLevelSameEvent.value) {
            continue;
          }
        }

        // EMCal photon v1
        energyCorrectionFactor = fEMCalCorrectionFactor->Eval(g1.e() > MinEnergy ? g1.e() : MinEnergy);
        ROOT::Math::PtEtaPhiMVector v1(energyCorrectionFactor * g1.pt(), g1.eta(), g1.phi(), 0.);
        // PCM photon v2s
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;

        float dTheta = v1.Theta() - v2.Theta();
        float dPhi = v1.Phi() - v2.Phi();
        float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));
        registry.fill(HIST("hMesonCuts"), 1);
        if (openingAngle <= mesonConfig.minOpenAngle) {
          registry.fill(HIST("hMesonCuts"), 2);
          continue;
        }
        if (rotationConfig.cfgDoRotation.value) {
          rotationBackgroundPCMEMC<decltype(photonsPCMPerCollision), decltype(photonsEMCPerCollision)>(vMeson, v2, v1, photonsPCMPerCollision, photonsEMCPerCollision, g2.globalIndex(), g1.globalIndex(), collision);
        }
        if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
          registry.fill(HIST("hMesonCuts"), 3);
          continue;
        }
        if (thnConfigAxisPt.value[1] > vMeson.Pt() || thnConfigAxisPt.value.back() < vMeson.Pt()) {
          registry.fill(HIST("hMesonCuts"), 4);
          continue;
        }
        if (mesonConfig.cfgEnableQA.value) {
          registry.fill(HIST("mesonQA/hInvMassPt"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("mesonQA/hTanThetaPhi"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
          registry.fill(HIST("mesonQA/hAlphaPt"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
          registry.fill(HIST("hMesonCuts"), 5);
          continue;
        }
        runFlowAnalysis<0>(collision, vMeson);
      }
    }
  }
  PROCESS_SWITCH(TaskPi0FlowEMC, processEMCalPCMC, "Process neutral meson flow using PCM-EMC same event", false);

  // PCM-EMCal mixed event
  void processEMCalPCMMixed(CollsWithQvecs const& collisions, EMCalPhotons const& clusters, PCMPhotons const& pcmPhotons, aod::V0Legs const&, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds)
  {
    if (clusters.size() <= 0 && pcmPhotons.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }

    using BinningTypeMixed = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, emevent::EP2FT0C<emevent::Q2xFT0C, emevent::Q2yFT0C>>;
    BinningTypeMixed binningOnPositions{{mixingConfig.cfgVtxBins, mixingConfig.cfgCentBins, mixingConfig.cfgEPBins}, true};

    auto associatedTables = std::make_tuple(clusters, pcmPhotons);

    Pair<CollsWithQvecs, EMCalPhotons, PCMPhotons, BinningTypeMixed> pairPCMEMC{binningOnPositions, mixingConfig.cfgMixingDepth, -1, collisions, associatedTables, &cache}; // indicates that mixingConfig.cfgMixingDepth events should be mixed and under/overflow (-1) to be ignored

    energyCorrectionFactor = 1.f;
    EMBitFlags emcFlags(clusters.size());
    if (clusters.size() > 0) {
      fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds, &registry);
    }

    EMBitFlags v0flags(pcmPhotons.size());
    if (pcmPhotons.size() > 0) {
      fV0PhotonCut.AreSelectedRunning<decltype(pcmPhotons), aod::V0Legs>(v0flags, pcmPhotons, &registry);
    }

    for (const auto& [c1, photonEMC, c2, photonPCM] : pairPCMEMC) {
      if (!(fEMEventCut.IsSelected(c1)) || !(fEMEventCut.IsSelected(c2))) {
        // general event selection
        continue;
      }
      if (!(eventcuts.cfgFT0COccupancyMin <= c1.ft0cOccupancyInTimeRange() && c1.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax) || !(eventcuts.cfgFT0COccupancyMin <= c2.ft0cOccupancyInTimeRange() && c2.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax)) {
        // occupancy selection
        continue;
      }
      if (getCentrality(c1) < eventcuts.cfgMinCent || getCentrality(c1) > eventcuts.cfgMaxCent || getCentrality(c2) < eventcuts.cfgMinCent || getCentrality(c2) > eventcuts.cfgMaxCent) {
        // event selection
        continue;
      }
      runNow = c1.runNumber();
      if (runNow != runBefore) {
        initCCDB(c1);
        runBefore = runNow;
      }
      registry.fill(HIST("h3DMixingCount"), c1.posZ(), getCentrality(c1), c1.ep2ft0m());
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photonEMC, photonPCM))) {
        if (!(emcFlags.test(g1.globalIndex())) || !(v0flags.test(g2.globalIndex()))) {
          continue;
        }
        // Cut edge clusters away, similar to rotation method to ensure same acceptance is used
        if (cfgDistanceToEdge.value) {
          if (checkEtaPhi1D(g1.eta(), RecoDecay::constrainAngle(g1.phi())) >= cfgEMCalMapLevelBackground.value) {
            continue;
          }
        }
        energyCorrectionFactor = fEMCalCorrectionFactor->Eval(g1.e() > MinEnergy ? g1.e() : MinEnergy);
        ROOT::Math::PtEtaPhiMVector v1(energyCorrectionFactor * g1.pt(), g1.eta(), g1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;

        float dTheta = v1.Theta() - v2.Theta();
        float dPhi = v1.Phi() - v2.Phi();
        float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));

        registry.fill(HIST("hMesonCutsMixed"), 1);
        if (openingAngle <= mesonConfig.minOpenAngle) {
          registry.fill(HIST("hMesonCutsMixed"), 2);
          continue;
        }
        if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
          registry.fill(HIST("hMesonCutsMixed"), 3);
          continue;
        }
        if (thnConfigAxisPt.value[1] > vMeson.Pt() || thnConfigAxisPt.value.back() < vMeson.Pt()) {
          registry.fill(HIST("hMesonCutsMixed"), 4);
          continue;
        }
        if (mesonConfig.cfgEnableQA.value) {
          registry.fill(HIST("mesonQA/hInvMassPtMixed"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("mesonQA/hTanThetaPhiMixed"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
          registry.fill(HIST("mesonQA/hAlphaPtMixed"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
          registry.fill(HIST("hMesonCutsMixed"), 5);
          continue;
        }
        registry.fill(HIST("hMesonCutsMixed"), 6);
        runFlowAnalysis<2>(c1, vMeson);
      }
    }
  }
  PROCESS_SWITCH(TaskPi0FlowEMC, processEMCalPCMMixed, "Process neutral meson flow using PCM-EMC mixed event", false);

  // Pi0 from EMCal
  void processM02(CollsWithQvecs const& collisions, EMCalPhotons const& clusters, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds)
  {
    if (clusters.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    EMBitFlags emcFlags(clusters.size());
    fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds);

    for (const auto& collision : collisions) {
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&registry, collision);
      if (!(fEMEventCut.IsSelected(collision))) {
        // general event selection
        continue;
      }
      if (!(eventcuts.cfgFT0COccupancyMin <= collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax)) {
        // occupancy selection
        continue;
      }
      float cent = getCentrality(collision);
      if (cent < eventcuts.cfgMinCent || cent > eventcuts.cfgMaxCent) {
        // event selection
        continue;
      }
      if (!isQvecGood(getAllQvec(collision))) {
        // selection based on QVector
        continue;
      }
      runNow = collision.runNumber();
      if (runNow != runBefore) {
        initCCDB(collision);
        runBefore = runNow;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&registry, collision);
      registry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // accepted
      registry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // accepted

      auto photonsPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());

      for (const auto& photon : photonsPerCollision) {
        if (emccuts.cfgEnableQA.value) {
          registry.fill(HIST("clusterQA/hEClusterBefore"), photon.e());                      // before cuts
          registry.fill(HIST("clusterQA/hClusterEtaPhiBefore"), photon.phi(), photon.eta()); // before cuts
        }
        auto matchedPrimsPerCluster = matchedPrims.sliceBy(perEMCClusterMT, photon.globalIndex());
        auto matchedSecondsPerCluster = matchedSeconds.sliceBy(perEMCClusterMS, photon.globalIndex());
        if (!(emcFlags.test(photon.globalIndex()))) {
          continue;
        }
        if (cfgDistanceToEdge.value && (checkEtaPhi1D(photon.eta(), RecoDecay::constrainAngle(photon.phi())) >= cfgEMCalMapLevelSameEvent.value)) {
          continue;
        }
        if (emccuts.cfgEnableQA.value) {
          registry.fill(HIST("clusterQA/hEClusterAfter"), photon.e());                      // accepted after cuts
          registry.fill(HIST("clusterQA/hClusterEtaPhiAfter"), photon.phi(), photon.eta()); // after cuts
        }

        auto [xQVec, yQVec] = getQvec(collision, qvecDetector);
        float cent = getCentrality(collision);

        float phiCand = photon.phi();

        float cosNPhi = std::cos(harmonic * phiCand);
        float sinNPhi = std::sin(harmonic * phiCand);
        float scalprodCand = cosNPhi * xQVec + sinNPhi * yQVec;

        if (correctionConfig.cfgApplySPresolution.value) {
          scalprodCand = scalprodCand / h1SPResolution->GetBinContent(h1SPResolution->FindBin(cent + epsilon));
        }
        if (cfgDoM02.value) {
          registry.fill(HIST("p3DM02Flow"), photon.m02(), photon.pt(), cent, scalprodCand);
          registry.fill(HIST("h3DSparsePi0"), photon.m02(), photon.pt(), cent);
        }
        return;
      } // end of loop over single cluster
    } // end of loop over collisions
  } // processM02
  PROCESS_SWITCH(TaskPi0FlowEMC, processM02, "Process single EMCal clusters as function of M02", false);

  // Pi0 from EMCal
  void processPCM(CollsWithQvecs const& collisions, PCMPhotons const& photons, aod::V0Legs const&)
  {
    if (photons.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    EMBitFlags v0flags(photons.size());
    fV0PhotonCut.AreSelectedRunning<decltype(photons), aod::V0Legs>(v0flags, photons, &registry);
    for (const auto& collision : collisions) {

      if (!isFullEventSelected(collision, true)) {
        continue;
      }
      runNow = collision.runNumber();
      if (runNow != runBefore) {
        initCCDB(collision);
        runBefore = runNow;
      }

      auto photonsPerCollision = photons.sliceBy(perCollisionPCM, collision.globalIndex());
      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photonsPerCollision, photonsPerCollision))) {
        if (!(v0flags.test(g1.globalIndex())) || !(v0flags.test(g2.globalIndex()))) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;

        float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));

        registry.fill(HIST("hMesonCuts"), 1);
        if (openingAngle <= mesonConfig.minOpenAngle) {
          registry.fill(HIST("hMesonCuts"), 2);
          continue;
        }
        if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
          registry.fill(HIST("hMesonCuts"), 3);
          continue;
        }
        if (thnConfigAxisPt.value[1] > vMeson.Pt() || thnConfigAxisPt.value.back() < vMeson.Pt()) {
          registry.fill(HIST("hMesonCuts"), 4);
          continue;
        }
        if (mesonConfig.cfgEnableQA.value) {
          registry.fill(HIST("mesonQA/hInvMassPt"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("mesonQA/hAlphaPt"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        registry.fill(HIST("hMesonCuts"), 6);
        runFlowAnalysis<0>(collision, vMeson);
      }
    } // end of loop over collisions
  }
  PROCESS_SWITCH(TaskPi0FlowEMC, processPCM, "Process PCM Pi0 candidates", false);

  // PCM-EMCal mixed event
  void processPCMMixed(FilteredCollsWithQvecs const& collisions, PCMPhotons const& pcmPhotons, aod::V0Legs const&)
  {

    if (pcmPhotons.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }

    using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, emevent::EP2FT0C<emevent::Q2xFT0C, emevent::Q2yFT0C>>;
    BinningType binningMixedEvent{{mixingConfig.cfgVtxBins, mixingConfig.cfgCentBins, mixingConfig.cfgEPBins}, true};

    auto pcmPhotonTuple = std::make_tuple(pcmPhotons);
    SameKindPair<FilteredCollsWithQvecs, PCMPhotons, BinningType> pair{binningMixedEvent, mixingConfig.cfgMixingDepth, -1, collisions, pcmPhotonTuple, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

    energyCorrectionFactor = 1.f;

    EMBitFlags v0flags(pcmPhotons.size());
    fV0PhotonCut.AreSelectedRunning<decltype(pcmPhotons), aod::V0Legs>(v0flags, pcmPhotons);

    for (const auto& [c1, photon1, c2, photon2] : pair) {
      if (!(fEMEventCut.IsSelected(c1)) || !(fEMEventCut.IsSelected(c2))) {
        // general event selection
        continue;
      }
      if (!(eventcuts.cfgFT0COccupancyMin <= c1.ft0cOccupancyInTimeRange() && c1.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax) || !(eventcuts.cfgFT0COccupancyMin <= c2.ft0cOccupancyInTimeRange() && c2.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax)) {
        // occupancy selection
        continue;
      }
      if (getCentrality(c1) < eventcuts.cfgMinCent || getCentrality(c1) > eventcuts.cfgMaxCent || getCentrality(c2) < eventcuts.cfgMinCent || getCentrality(c2) > eventcuts.cfgMaxCent) {
        // event selection
        continue;
      }
      runNow = c1.runNumber();
      if (runNow != runBefore) {
        initCCDB(c1);
        runBefore = runNow;
      }
      registry.fill(HIST("h3DMixingCount"), c1.posZ(), getCentrality(c1), c1.ep2ft0m());
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photon1, photon2))) {
        if (!(v0flags.test(g1.globalIndex())) || !(v0flags.test(g2.globalIndex()))) {
          continue;
        }
        energyCorrectionFactor = fEMCalCorrectionFactor->Eval(g1.e() > MinEnergy ? g1.e() : MinEnergy);
        ROOT::Math::PtEtaPhiMVector v1(energyCorrectionFactor * g1.pt(), g1.eta(), g1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;

        float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));

        registry.fill(HIST("hMesonCutsMixed"), 1);
        if (openingAngle <= mesonConfig.minOpenAngle) {
          registry.fill(HIST("hMesonCutsMixed"), 2);
          continue;
        }
        if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
          registry.fill(HIST("hMesonCutsMixed"), 3);
          continue;
        }
        if (thnConfigAxisPt.value[1] > vMeson.Pt() || thnConfigAxisPt.value.back() < vMeson.Pt()) {
          registry.fill(HIST("hMesonCutsMixed"), 4);
          continue;
        }
        if (mesonConfig.cfgEnableQA.value) {
          registry.fill(HIST("mesonQA/hInvMassPtMixed"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("mesonQA/hAlphaPtMixed"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        registry.fill(HIST("hMesonCutsMixed"), 6);
        runFlowAnalysis<2>(c1, vMeson);
      }
    }
  }
  PROCESS_SWITCH(TaskPi0FlowEMC, processPCMMixed, "Process neutral meson flow using PCM-EMC mixed event", false);

}; // End struct TaskPi0FlowEMC

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TaskPi0FlowEMC>(cfgc)};
}
