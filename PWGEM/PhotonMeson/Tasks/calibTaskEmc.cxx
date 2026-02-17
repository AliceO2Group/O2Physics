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

/// \file calibTaskEmc.cxx
/// \brief Task to produce calibration values for EMCal
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/PhotonMeson/Core/EMBitFlags.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/emcalHistoDefinitions.h"

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

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
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

struct CalibTaskEmc {
  static constexpr float MinEnergy = 0.7f;

  // configurable for flow
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int> cfgEMCalMapLevelBackground{"cfgEMCalMapLevelBackground", 4, "Different levels of correction for the background, the smaller number includes the level of the higher number (4: none, 3: only inside EMCal, 2: exclude bad channels, 1: remove edges)"};
  Configurable<int> cfgEMCalMapLevelSameEvent{"cfgEMCalMapLevelSameEvent", 4, "Different levels of correction for the same event, the smaller number includes the level of the higher number (4: none, 3: only inside EMCal, 2: exclude bad channels, 1: remove edges)"};
  Configurable<int> cfgDistanceToEdge{"cfgDistanceToEdge", 1, "Distance to edge in cells required for rotated cluster to be accepted"};
  Configurable<float> cfgMaxAsymmetry{"cfgMaxAsymmetry", 0.1f, "Maximum allowed asymmetry for photon pairs used in calibration."};

  // configurable axis
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {400, 0.0, 0.8}, ""};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {20, 0., 100.}, ""};
  ConfigurableAxis thnConfigAxisEnergyCalib{"thnConfigAxisEnergyCalib", {200, 0., 20.}, "Cluster energy axis for calibration"};
  ConfigurableAxis thnConfigAxisPtCalib{"thnConfigAxisPtCalib", {200, 0., 20.}, "Pt axis for calibration, mostly for PCM"};

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
    Configurable<bool> cfgEnableQA{"cfgEnableQA", false, "flag to turn QA plots on/off"};
    ConfigurableAxis thConfigAxisTanThetaPhi{"thConfigAxisTanThetaPhi", {180, -90.f, 90.f}, ""};
  } mesonConfig;

  struct : ConfigurableGroup {
    std::string prefix = "mixingConfig";
    ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis cfgCentBins{"cfgCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "Mixing bins - centrality"};
    Configurable<int> cfgMixingDepth{"cfgMixingDepth", 2, "Mixing depth"};
  } mixingConfig;

  struct : ConfigurableGroup {
    std::string prefix = "correctionConfig";
    Configurable<bool> cfgEnableNonLin{"cfgEnableNonLin", false, "flag to turn extra non linear energy calibration on/off"};
  } correctionConfig;

  SliceCache cache;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  int runNow = 0;
  int runBefore = -1;

  Filter collisionFilter = (nabs(aod::collision::posZ) <= eventcuts.cfgZvtxMax) && (aod::evsel::ft0cOccupancyInTimeRange <= eventcuts.cfgFT0COccupancyMax) && (aod::evsel::ft0cOccupancyInTimeRange >= eventcuts.cfgFT0COccupancyMin);
  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::MinClusters>;
  using PCMPhotons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
  using FilteredCollsWithQvecs = soa::Filtered<soa::Join<aod::EMEvents, aod::EMEventsAlias, aod::EMEventsMult, aod::EMEventsCent>>;
  using CollsWithQvecs = soa::Join<aod::EMEvents, aod::EMEventsAlias, aod::EMEventsMult, aod::EMEventsCent>;
  using Colls = soa::Join<aod::EMEvents, aod::EMEventsAlias, aod::EMEventsMult, aod::EMEventsCent>;

  static constexpr std::size_t NQVecEntries = 6;

  PresliceOptional<EMCalPhotons> perCollisionEMC = o2::aod::emccluster::emeventId;
  PresliceOptional<PCMPhotons> perCollisionPCM = aod::v0photonkf::emeventId;
  PresliceOptional<MinMTracks> perEMCClusterMT = o2::aod::mintm::minClusterId;
  PresliceOptional<MinMSTracks> perEMCClusterMS = o2::aod::mintm::minClusterId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  o2::emcal::Geometry* emcalGeom;
  o2::emcal::BadChannelMap* mBadChannels;
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

    defineEMEventCut();
    defineEMCCut();
    fEMCCut.addQAHistograms(&registry);
    definePCMCut();
    fV0PhotonCut.addQAHistograms(&registry);
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&registry);

    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality (%)"};
    const AxisSpec thAxisTanThetaPhi{mesonConfig.thConfigAxisTanThetaPhi, "atan(#Delta#theta/#Delta#varphi)"};
    const AxisSpec thAxisEnergyCalib{thnConfigAxisEnergyCalib, "#it{E}_{clus} (GeV)"};
    const AxisSpec thnAxisPtCalib{thnConfigAxisPtCalib, "#it{p}_{T} (GeV)"};
    const AxisSpec thAxisAlpha{100, -1., +1, "#alpha"};
    const AxisSpec thAxisEnergy{1000, 0., 100., "#it{E}_{clus} (GeV)"};
    const AxisSpec thAxisEta{320, -0.8, 0.8, "#eta"};
    const AxisSpec thAxisPhi{500, 0, 2 * 3.14159, "phi"};

    const AxisSpec thnAxisMixingVtx{mixingConfig.cfgVtxBins, "#it{z} (cm)"};
    const AxisSpec thnAxisMixingCent{mixingConfig.cfgCentBins, "Centrality (%)"};

    if (doprocessEMCal || doprocessEMCalPCMC) {
      registry.add("hSparsePi0", "m_{inv} vs p_T vs cent for same event", HistType::kTH3D, {thnAxisInvMass, thAxisEnergyCalib, thnAxisCent});
    } else if (doprocessPCM) {
      registry.add("hSparsePi0", "m_{inv} vs p_T vs cent for same event", HistType::kTH3D, {thnAxisInvMass, thnAxisPtCalib, thnAxisCent});
    }

    if (doprocessEMCalMixed || doprocessEMCalPCMMixed) {
      registry.add("hSparseBkgMix", "m_{inv} vs p_T vs cent for mixed event", HistType::kTH3D, {thnAxisInvMass, thAxisEnergyCalib, thnAxisCent});
    } else if (doprocessPCMMixed) {
      registry.add("hSparseBkgMix", "m_{inv} vs p_T vs cent for mixed event", HistType::kTH3D, {thnAxisInvMass, thnAxisPtCalib, thnAxisCent});
    }

    if (doprocessEMCalMixed || doprocessEMCalPCMMixed || doprocessPCMMixed) {
      registry.add("hMixingCount", "THn Event Mixing QA", HistType::kTH2D, {thnAxisMixingVtx, thnAxisMixingCent});
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

    if (mesonConfig.cfgEnableQA.value) {
      registry.add("mesonQA/hInvMassPt", "Histo for inv pair mass vs pt", HistType::kTH2D, {thnAxisInvMass, thnAxisPtCalib});
      registry.add("mesonQA/hTanThetaPhi", "Histo for identification of conversion cluster", HistType::kTH2D, {thnAxisInvMass, thAxisTanThetaPhi});
      registry.add("mesonQA/hAlphaPt", "Histo of meson asymmetry vs pT", HistType::kTH2D, {thAxisAlpha, thnAxisPtCalib});
      registry.add("mesonQA/hInvMassPtMixed", "Histo for inv pair mass vs pt for mixed event", HistType::kTH2D, {thnAxisInvMass, thnAxisPtCalib});
      registry.add("mesonQA/hTanThetaPhiMixed", "Histo for identification of conversion cluster for mixed event", HistType::kTH2D, {thnAxisInvMass, thAxisTanThetaPhi});
      registry.add("mesonQA/hAlphaPtMixed", "Histo of meson asymmetry vs pT for mixed event", HistType::kTH2D, {thAxisAlpha, thnAxisPtCalib});
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

  }; // end init

  /// Change radians to degree
  /// \param angle in radians
  /// \return angle in degree
  constexpr float getAngleDegree(float angle)
  {
    return angle * o2::constants::math::Rad2Deg;
  }

  /// Fill THnSparse
  /// \param mass is the invariant mass of the candidate
  /// \param pt is the transverse momentum of the candidate
  /// \param cent is the centrality of the collision
  /// \param sp is the scalar product
  template <const int histType>
  void fillThn(const float mass, const float pt, const float cent)
  {
    static constexpr std::string_view HistTypes[3] = {"hSparsePi0", "hSparseBkgRot", "hSparseBkgMix"};
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
  }

  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param meson are the selected candidates
  template <const int histType, o2::soa::is_iterator TCollision>
  void runFlowAnalysis(TCollision const& collision, ROOT::Math::PtEtaPhiMVector const& meson, float photonCalibValue = 0.f)
  {
    float cent = getCentrality(collision);
    float massCand = meson.M();

    fillThn<histType>(massCand, photonCalibValue, cent);
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
      float asymmetry = (g1.e() - g2.e()) / (g1.e() + g2.e());
      if (std::fabs(asymmetry) > cfgMaxAsymmetry) { // only use symmetric decays
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

      ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
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
      if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
        registry.fill(HIST("hMesonCuts"), 3);
        continue;
      }
      if (mesonConfig.cfgEnableQA.value) {
        registry.fill(HIST("mesonQA/hInvMassPt"), vMeson.M(), vMeson.Pt());
        registry.fill(HIST("mesonQA/hTanThetaPhi"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
        registry.fill(HIST("mesonQA/hAlphaPt"), asymmetry, vMeson.Pt());
      }
      if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
        registry.fill(HIST("hMesonCuts"), 5);
        continue;
      }
      registry.fill(HIST("hMesonCuts"), 6);
      runFlowAnalysis<0>(collision, vMeson, g1.e());
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
      runPairingLoop(collision, photonsPerCollision, photonsPerCollision, flags, flags);
    } // end of loop over collisions
  }
  PROCESS_SWITCH(CalibTaskEmc, processEMCal, "Process EMCal Pi0 candidates", true);

  // Pi0 from EMCal
  void processEMCalMixed(FilteredCollsWithQvecs const& collisions, EMCalPhotons const& clusters, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds)
  {
    if (clusters.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    EMBitFlags flags(clusters.size());
    fEMCCut.AreSelectedRunning(flags, clusters, matchedPrims, matchedSeconds);

    using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
    BinningType binningMixedEvent{{mixingConfig.cfgVtxBins, mixingConfig.cfgCentBins}, true};

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
      runNow = c1.runNumber();
      if (runNow != runBefore) {
        initCCDB(c1);
        runBefore = runNow;
      }
      registry.fill(HIST("hMixingCount"), c1.posZ(), getCentrality(c1));
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(clusters1, clusters2))) {
        if (!(flags.test(g1.globalIndex())) || !(flags.test(g2.globalIndex()))) {
          continue;
        }
        float asymmetry = (g1.e() - g2.e()) / (g1.e() + g2.e());
        if (std::fabs(asymmetry) > cfgMaxAsymmetry) { // only use symmetric decays
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
        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
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
        if (mesonConfig.cfgEnableQA.value) {
          registry.fill(HIST("mesonQA/hInvMassPtMixed"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("mesonQA/hTanThetaPhiMixed"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
          registry.fill(HIST("mesonQA/hAlphaPtMixed"), asymmetry, vMeson.Pt());
        }
        if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
          registry.fill(HIST("hMesonCutsMixed"), 5);
          continue;
        }
        registry.fill(HIST("hMesonCutsMixed"), 6);
        runFlowAnalysis<2>(c1, vMeson, g1.e());
      }
    }
  }
  PROCESS_SWITCH(CalibTaskEmc, processEMCalMixed, "Process EMCal Pi0 mixed event candidates", false);

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
        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
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
        if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
          registry.fill(HIST("hMesonCuts"), 3);
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
        runFlowAnalysis<0>(collision, vMeson, g1.pt());
      }
    }
  }
  PROCESS_SWITCH(CalibTaskEmc, processEMCalPCMC, "Process neutral meson flow using PCM-EMC same event", false);

  // PCM-EMCal mixed event
  void processEMCalPCMMixed(CollsWithQvecs const& collisions, EMCalPhotons const& clusters, PCMPhotons const& pcmPhotons, aod::V0Legs const&, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds)
  {
    if (clusters.size() <= 0 && pcmPhotons.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    using BinningTypeMixed = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
    BinningTypeMixed binningOnPositions{{mixingConfig.cfgVtxBins, mixingConfig.cfgCentBins}, true};

    auto associatedTables = std::make_tuple(clusters, pcmPhotons);

    Pair<CollsWithQvecs, EMCalPhotons, PCMPhotons, BinningTypeMixed> pairPCMEMC{binningOnPositions, mixingConfig.cfgMixingDepth, -1, collisions, associatedTables, &cache}; // indicates that mixingConfig.cfgMixingDepth events should be mixed and under/overflow (-1) to be ignored

    EMBitFlags emcFlags(clusters.size());
    if (clusters.size() > 0) {
      fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds);
    }

    EMBitFlags v0flags(pcmPhotons.size());
    if (pcmPhotons.size() > 0) {
      fV0PhotonCut.AreSelectedRunning<decltype(pcmPhotons), aod::V0Legs>(v0flags, pcmPhotons);
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
      registry.fill(HIST("hMixingCount"), c1.posZ(), getCentrality(c1));
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
        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
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
        runFlowAnalysis<2>(c1, vMeson, g1.pt());
      }
    }
  }
  PROCESS_SWITCH(CalibTaskEmc, processEMCalPCMMixed, "Process neutral meson flow using PCM-EMC mixed event", false);

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

        float asymmetry = (g1.pt() - g2.pt()) / (g1.pt() + g2.pt());
        if (std::fabs(asymmetry) > cfgMaxAsymmetry) { // only use symmetric decays
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
        if (mesonConfig.cfgEnableQA.value) {
          registry.fill(HIST("mesonQA/hInvMassPt"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("mesonQA/hAlphaPt"), asymmetry, vMeson.Pt());
        }
        registry.fill(HIST("hMesonCuts"), 6);
        runFlowAnalysis<0>(collision, vMeson, g1.pt());
      }
    } // end of loop over collisions
  }
  PROCESS_SWITCH(CalibTaskEmc, processPCM, "Process PCM Pi0 candidates", false);

  // PCM-EMCal mixed event
  void processPCMMixed(FilteredCollsWithQvecs const& collisions, PCMPhotons const& pcmPhotons, aod::V0Legs const&)
  {
    if (pcmPhotons.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
    BinningType binningMixedEvent{{mixingConfig.cfgVtxBins, mixingConfig.cfgCentBins}, true};

    auto pcmPhotonTuple = std::make_tuple(pcmPhotons);
    SameKindPair<FilteredCollsWithQvecs, PCMPhotons, BinningType> pair{binningMixedEvent, mixingConfig.cfgMixingDepth, -1, collisions, pcmPhotonTuple, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

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
      registry.fill(HIST("hMixingCount"), c1.posZ(), getCentrality(c1));
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photon1, photon2))) {
        if (!(v0flags.test(g1.globalIndex())) || !(v0flags.test(g2.globalIndex()))) {
          continue;
        }

        float asymmetry = (g1.pt() - g2.pt()) / (g1.pt() + g2.pt());
        if (std::fabs(asymmetry) > cfgMaxAsymmetry) { // only use symmetric decays
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
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
        if (mesonConfig.cfgEnableQA.value) {
          registry.fill(HIST("mesonQA/hInvMassPtMixed"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("mesonQA/hAlphaPtMixed"), asymmetry, vMeson.Pt());
        }
        registry.fill(HIST("hMesonCutsMixed"), 6);
        runFlowAnalysis<2>(c1, vMeson, g1.pt());
      }
    }
  }
  PROCESS_SWITCH(CalibTaskEmc, processPCMMixed, "Process neutral meson flow using PCM-EMC mixed event", false);

}; // End struct CalibTaskEmc

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<CalibTaskEmc>(cfgc)};
}
