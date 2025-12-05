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

#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/emcalHistoDefinitions.h"
//
#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"

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
#include <memory>
#include <numbers>
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

struct CalibTaskEmc {
  static constexpr float MinEnergy = 0.7f;

  // configurable for flow
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> cfgDoRotation{"cfgDoRotation", false, "Flag to enable rotation background method"};
  Configurable<int> cfgDownsampling{"cfgDownsampling", 1, "Calculate rotation background only for every <value> collision"};
  Configurable<int> cfgEMCalMapLevelBackground{"cfgEMCalMapLevelBackground", 4, "Different levels of correction for the background, the smaller number includes the level of the higher number (4: none, 3: only inside EMCal, 2: exclude bad channels, 1: remove edges)"};
  Configurable<int> cfgEMCalMapLevelSameEvent{"cfgEMCalMapLevelSameEvent", 4, "Different levels of correction for the same event, the smaller number includes the level of the higher number (4: none, 3: only inside EMCal, 2: exclude bad channels, 1: remove edges)"};
  Configurable<float> cfgRotAngle{"cfgRotAngle", std::move(const_cast<float&>(o2::constants::math::PIHalf)), "Angle used for the rotation method"};
  Configurable<int> cfgDistanceToEdge{"cfgDistanceToEdge", 1, "Distance to edge in cells required for rotated cluster to be accepted"};
  Configurable<float> cfgMaxAsymmetry{"cfgMaxAsymmetry", 0.1f, "Maximum allowed asymmetry for photon pairs used in calibration when using EMC-EMC."};

  // configurable axis
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {400, 0.0, 0.8}, ""};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0., 20.}, ""};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {20, 0., 100.}, ""};
  ConfigurableAxis thnConfigAxisEnergyCalib{"thnConfigAxisEnergyCalib", {200, 0., 20.}, ""};

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
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -1, "min. track occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -1, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<float> cfgMinCent{"cfgMinCent", 0, "min. centrality (%)"};
    Configurable<float> cfgMaxCent{"cfgMaxCent", 90, "max. centrality (%)"};
    Configurable<bool> onlyKeepWeightedEvents{"onlyKeepWeightedEvents", false, "flag to keep only weighted events (for JJ MCs) and remove all MB events (with weight = 1)"};
    Configurable<bool> cfgEnableQA{"cfgEnableQA", false, "flag to turn QA plots on/off"};
  } eventcuts;

  EMCPhotonCut fEMCCut;
  struct : ConfigurableGroup {
    std::string prefix = "emccuts";
    Configurable<std::string> clusterDefinition{"clusterDefinition", "kV3Default", "Clusterizer to be selected, e.g. V3Default"};
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
    Configurable<bool> cfgEMCUseTM{"cfgEMCUseTM", true, "flag to use EMCal track matching cut or not"};
    Configurable<bool> emcUseSecondaryTM{"emcUseSecondaryTM", false, "flag to use EMCal secondary track matching cut or not"};
    Configurable<bool> cfgEnableQA{"cfgEnableQA", false, "flag to turn QA plots on/off"};
  } emccuts;

  V0PhotonCut fPCMPhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcuts";
    Configurable<bool> cfgRequireV0WithITSTPC{"cfgRequireV0WithITSTPC", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfgRequireV0WithITSonly{"cfgRequireV0WithITSonly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfgRequireV0WithTPConly{"cfgRequireV0WithTPConly", false, "flag to select V0s with TPConly tracks"};
    Configurable<float> cfgMinPtV0{"cfgMinPtV0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfgMaxPtV0{"cfgMaxPtV0", 1e+10, "max pT for v0 photons at PV"};
    Configurable<float> cfgMinEtaV0{"cfgMinEtaV0", -0.8, "min eta for v0 photons at PV"};
    Configurable<float> cfgMaxEtaV0{"cfgMaxEtaV0", 0.8, "max eta for v0 photons at PV"};
    Configurable<float> cfgMinV0Radius{"cfgMinV0Radius", 4.0, "min v0 radius"};
    Configurable<float> cfgMaxV0Radius{"cfgMaxV0Radius", 90.0, "max v0 radius"};
    Configurable<float> cfgMaxAlphaAP{"cfgMaxAlphaAP", 0.95, "max alpha for AP cut"};
    Configurable<float> cfgMaxQtAP{"cfgMaxQtAP", 0.01, "max qT for AP cut"};
    Configurable<float> cfgMinCosPA{"cfgMinCosPA", 0.999, "min V0 CosPA"};
    Configurable<float> cfgMaxPCA{"cfgMaxPCA", 1.5, "max distance btween 2 legs"};
    Configurable<float> cfgMaxChi2KF{"cfgMaxChi2KF", 1e+10, "max chi2/ndf with KF"};
    Configurable<bool> cfgRejectV0OnITSib{"cfgRejectV0OnITSib", true, "flag to reject V0s on ITSib"};

    Configurable<int> cfgMinNClusterTPC{"cfgMinNClusterTPC", 0, "min ncluster tpc"};
    Configurable<int> cfgMinNCrossedRows{"cfgMinNCrossedRows", 40, "min ncrossed rows"};
    Configurable<float> cfgMaxFracSharedClusterTPC{"cfgMaxFracSharedClusterTPC", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfgMaxChi2TPC{"cfgMaxChi2TPC", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfgMaxChi2ITS{"cfgMaxChi2ITS", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfgMinTPCNSigmaEl{"cfgMinTPCNSigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfgMaxTPCNSigmaEl{"cfgMaxTPCNSigmaEl", +3.0, "max. TPC n sigma for electron"};
    Configurable<bool> cfgDisableITSOnly{"cfgDisableITSOnly", false, "flag to disable ITSonly tracks"};
    Configurable<bool> cfgDisableTPCOnly{"cfgDisableTPCOnly", false, "flag to disable TPConly tracks"};
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
    ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -9.0f, -8.f, -7.0f, -6.f, -5.0f, -4.f, -3.f, -2.f, -1.f, 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis cfgCentBins{"cfgCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "Mixing bins - centrality"};
    Configurable<int> cfgMixingDepth{"cfgMixingDepth", 2, "Mixing depth"};
  } mixingConfig;

  SliceCache cache;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  int runNow = 0;
  int runBefore = -1;

  // Filter clusterFilter = aod::skimmedcluster::time >= emccuts.cfgEMCminTime && aod::skimmedcluster::time <= emccuts.cfgEMCmaxTime && aod::skimmedcluster::m02 >= emccuts.cfgEMCminM02 && aod::skimmedcluster::m02 <= emccuts.cfgEMCmaxM02 && skimmedcluster::e >= emccuts.cfgEMCminE;
  // Filter collisionFilter = (nabs(aod::collision::posZ) <= eventcuts.cfgZvtxMax) && (aod::evsel::ft0cOccupancyInTimeRange <= eventcuts.cfgFT0COccupancyMax) && (aod::evsel::ft0cOccupancyInTimeRange >= eventcuts.cfgFT0COccupancyMin);
  // using FilteredEMCalPhotons = soa::Filtered<soa::Join<aod::EMCEMEventIds, aod::SkimEMCClusters>>;
  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::SkimEMCClusters>;
  using PCMPhotons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
  using Colls = soa::Join<aod::EMEvents, aod::EMEventsAlias, aod::EMEventsMult, aod::EMEventsCent>;

  // for event mixing
  using MyEMH = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, o2::aod::pwgem::dilepton::utils::EMTrack>;
  MyEMH* emh1 = nullptr;
  MyEMH* emh2 = nullptr;

  Preslice<EMCalPhotons> perCollisionEMC = aod::emccluster::emeventId;
  Preslice<PCMPhotons> perCollisionPCM = aod::v0photonkf::emeventId;

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningType binningOnPositions{{mixingConfig.cfgVtxBins, mixingConfig.cfgCentBins}, true};
  Pair<Colls, EMCalPhotons, PCMPhotons, BinningType> pairPCMEMC{binningOnPositions, mixingConfig.cfgMixingDepth, -1, &cache}; // indicates that mixingConfig.cfgMixingDepth events should be mixed and under/overflow (-1) to be ignored

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  o2::emcal::Geometry* emcalGeom;
  o2::emcal::BadChannelMap* mBadChannels;
  // Constants for eta and phi ranges
  double etaMin = -0.75, etaMax = 0.75;
  int nBinsEta = 150; // 150 bins for eta

  double phiMin = 1.35, phiMax = 5.75;
  int nBinsPhi = 440; // (440 bins = 0.01 step size covering most regions)

  std::vector<int8_t> lookupTable1D;
  float epsilon = 1.e-8;

  // static constexpr
  static constexpr int64_t NMinPhotonRotBkg = 3;
  static constexpr int64_t NMinPhotonRotBkgMixed = 2;

  // Usage when cfgEnableNonLin is enabled
  std::unique_ptr<TF1> fEMCalCorrectionFactor; // ("fEMCalCorrectionFactor","(1 + [0]/x + [1]/x^2) / (1 + [2]/x)", 0.3, 100.);

  // To access the 1D array
  inline int getIndex(int iEta, int iPhi)
  {
    return iEta * nBinsPhi + iPhi;
  }

  // Function to access the lookup table
  inline int8_t checkEtaPhi1D(double eta, double phi)
  {
    if (eta < etaMin || eta > etaMax || phi < phiMin || phi > phiMax) {
      return 3; // Out of bounds
    }

    // Compute indices directly
    int iEta = static_cast<int>((eta - etaMin) / ((etaMax - etaMin) / nBinsEta));
    int iPhi = static_cast<int>((phi - phiMin) / ((phiMax - phiMin) / nBinsPhi));

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
  }

  void DefinePCMCut()
  {
    fPCMPhotonCut = V0PhotonCut("fPCMPhotonCut", "fPCMPhotonCut");

    // for v0
    fPCMPhotonCut.SetV0PtRange(pcmcuts.cfgMinPtV0, pcmcuts.cfgMaxPtV0);
    fPCMPhotonCut.SetV0EtaRange(pcmcuts.cfgMinEtaV0, pcmcuts.cfgMaxEtaV0);
    fPCMPhotonCut.SetMinCosPA(pcmcuts.cfgMinCosPA);
    fPCMPhotonCut.SetMaxPCA(pcmcuts.cfgMaxPCA);
    fPCMPhotonCut.SetMaxChi2KF(pcmcuts.cfgMaxChi2KF);
    fPCMPhotonCut.SetRxyRange(pcmcuts.cfgMinV0Radius, pcmcuts.cfgMaxV0Radius);
    fPCMPhotonCut.SetAPRange(pcmcuts.cfgMaxAlphaAP, pcmcuts.cfgMaxQtAP);
    fPCMPhotonCut.RejectITSib(pcmcuts.cfgRejectV0OnITSib);

    // for track
    fPCMPhotonCut.SetMinNClustersTPC(pcmcuts.cfgMinNClusterTPC);
    fPCMPhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfgMinNCrossedRows);
    fPCMPhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fPCMPhotonCut.SetMaxFracSharedClustersTPC(pcmcuts.cfgMaxFracSharedClusterTPC);
    fPCMPhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfgMaxChi2TPC);
    fPCMPhotonCut.SetTPCNsigmaElRange(pcmcuts.cfgMinTPCNSigmaEl, pcmcuts.cfgMaxTPCNSigmaEl);
    fPCMPhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfgMaxChi2ITS);
    fPCMPhotonCut.SetNClustersITS(0, 7);
    fPCMPhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fPCMPhotonCut.SetDisableITSonly(pcmcuts.cfgDisableITSOnly);
    fPCMPhotonCut.SetDisableTPConly(pcmcuts.cfgDisableTPCOnly);
    fPCMPhotonCut.SetRequireITSTPC(pcmcuts.cfgRequireV0WithITSTPC);
    fPCMPhotonCut.SetRequireITSonly(pcmcuts.cfgRequireV0WithITSonly);
    fPCMPhotonCut.SetRequireTPConly(pcmcuts.cfgRequireV0WithTPConly);
  }

  void init(InitContext&)
  {
    defineEMEventCut();
    defineEMCCut();
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&registry);

    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality (%)"};
    const AxisSpec thAxisTanThetaPhi{mesonConfig.thConfigAxisTanThetaPhi, "atan(#Delta#theta/#Delta#varphi)"};
    const AxisSpec thAxisClusterEnergy{thnConfigAxisPt, "#it{E} (GeV)"};
    const AxisSpec thAxisEnergyCalib{thnConfigAxisEnergyCalib, "#it{E}_{clus} (GeV)"};
    const AxisSpec thAxisAlpha{100, -1., +1, "#alpha"};
    const AxisSpec thAxisEnergy{1000, 0., 100., "#it{E}_{clus} (GeV)"};

    registry.add("hSparseCalibSE", "THn for Calib same event", HistType::kTHnSparseF, {thnAxisInvMass, thAxisEnergyCalib, thnAxisCent});
    registry.add("hSparseCalibBack", "THn for Calib background", HistType::kTHnSparseF, {thnAxisInvMass, thAxisEnergyCalib, thnAxisCent});

    auto hClusterCuts = registry.add<TH1>("hClusterCuts", "hClusterCuts;;Counts", kTH1D, {{6, 0.5, 6.5}}, false);
    hClusterCuts->GetXaxis()->SetBinLabel(1, "in");
    hClusterCuts->GetXaxis()->SetBinLabel(2, "opening angle");
    hClusterCuts->GetXaxis()->SetBinLabel(3, "#it{M}_{#gamma#gamma}");
    hClusterCuts->GetXaxis()->SetBinLabel(4, "#it{p}_{T}");
    hClusterCuts->GetXaxis()->SetBinLabel(5, "conversion cut");
    hClusterCuts->GetXaxis()->SetBinLabel(6, "out");

    auto hClusterCutsMixed = registry.add<TH1>("hClusterCutsMixed", "hClusterCutsMixed;;Counts", kTH1D, {{6, 0.5, 6.5}}, false);
    hClusterCutsMixed->GetXaxis()->SetBinLabel(1, "in");
    hClusterCutsMixed->GetXaxis()->SetBinLabel(2, "opening angle");
    hClusterCutsMixed->GetXaxis()->SetBinLabel(3, "#it{M}_{#gamma#gamma}");
    hClusterCutsMixed->GetXaxis()->SetBinLabel(4, "#it{p}_{T}");
    hClusterCutsMixed->GetXaxis()->SetBinLabel(5, "conversion cut");
    hClusterCutsMixed->GetXaxis()->SetBinLabel(6, "out");

    if (eventcuts.cfgEnableQA) {
      auto hCollisionEMCCheck = registry.add<TH1>("hCollisionEMCCheck", "collision counter;;Counts", kTH1D, {{7, 0.5, 7.5}}, false);
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(1, "all");
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(2, "EMC MB Readout");
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(3, "has clusters");
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(4, "EMC MB Readout & has clusters");
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(5, "EMC MB Readout but no clusters");
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(6, "No EMC MB Readout but has clusters");
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(7, "No EMC MB Readout and no clusters");
    }

    if (emccuts.cfgEnableQA) {
      registry.add("hEClusterBefore", "Histo for cluster energy before cuts", HistType::kTH1D, {thAxisClusterEnergy});
      registry.add("hEClusterAfter", "Histo for cluster energy after cuts", HistType::kTH1D, {thAxisClusterEnergy});
    }

    if (mesonConfig.cfgEnableQA) {
      registry.add("hInvMassPt", "Histo for inv pair mass vs pt", HistType::kTH2D, {thnAxisInvMass, thnAxisPt});
      registry.add("hTanThetaPhi", "Histo for identification of conversion cluster", HistType::kTH2D, {thnAxisInvMass, thAxisTanThetaPhi});
      registry.add("hAlphaPt", "Histo of meson asymmetry vs pT", HistType::kTH2D, {thAxisAlpha, thnAxisPt});
      registry.add("hInvMassPtMixed", "Histo for inv pair mass vs pt for mixed event", HistType::kTH2D, {thnAxisInvMass, thnAxisPt});
      registry.add("hTanThetaPhiMixed", "Histo for identification of conversion cluster for mixed event", HistType::kTH2D, {thnAxisInvMass, thAxisTanThetaPhi});
      registry.add("hAlphaPtMixed", "Histo of meson asymmetry vs pT for mixed event", HistType::kTH2D, {thAxisAlpha, thnAxisPt});
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
  float getAngleDegree(float angle)
  {
    return angle * 180.f * std::numbers::inv_pi_v<float>;
  }

  /// Get the centrality
  /// \param collision is the collision with the centrality information
  template <typename TCollision>
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

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    // Load EMCal geometry
    emcalGeom = o2::emcal::Geometry::GetInstanceFromRunNumber(collision.runNumber());
    // Load Bad Channel map
    mBadChannels = ccdb->getForTimeStamp<o2::emcal::BadChannelMap>("EMC/Calib/BadChannelMap", collision.timestamp());
    lookupTable1D = std::vector<int8_t>(nBinsEta * nBinsPhi, -1);
    double binWidthEta = (etaMax - etaMin) / nBinsEta;
    double binWidthPhi = (phiMax - phiMin) / nBinsPhi;

    for (int iEta = 0; iEta < nBinsEta; ++iEta) {
      double etaCenter = etaMin + (iEta + 0.5) * binWidthEta;
      for (int iPhi = 0; iPhi < nBinsPhi; ++iPhi) {
        double phiCenter = phiMin + (iPhi + 0.5) * binWidthPhi;
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

  /// \brief Calculate background using rotation background method for calib
  template <typename TPhotons, typename TCollision>
  void rotationBackgroundCalib(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photon1, ROOT::Math::PtEtaPhiMVector photon2, TPhotons const& photons_coll, unsigned int ig1, unsigned int ig2, TCollision const& collision)
  {
    // if less than 3 clusters are present skip event since we need at least 3 clusters
    if (photons_coll.size() < NMinPhotonRotBkg) {
      return;
    }
    float cent = getCentrality(collision);
    int iCellIDPhoton1 = 0;
    int iCellIDPhoton2 = 0;

    ROOT::Math::AxisAngle rotationAxis(meson.Vect(), cfgRotAngle.value);
    ROOT::Math::Rotation3D rotationMatrix(rotationAxis);
    photon1 = rotationMatrix * photon1;
    photon2 = rotationMatrix * photon2;

    if (checkEtaPhi1D(photon1.Eta(), RecoDecay::constrainAngle(photon1.Phi())) >= cfgEMCalMapLevelBackground.value) {
      iCellIDPhoton1 = -1;
    }
    if (checkEtaPhi1D(photon2.Eta(), RecoDecay::constrainAngle(photon2.Phi())) >= cfgEMCalMapLevelBackground.value) {
      iCellIDPhoton2 = -1;
    }

    if (iCellIDPhoton1 == -1 && iCellIDPhoton2 == -1) {
      return;
    }
    for (const auto& photon : photons_coll) {
      if (photon.globalIndex() == ig1 || photon.globalIndex() == ig2) {
        // only combine rotated photons with other photons
        continue;
      }
      if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(photon))) {
        continue;
      }
      if (checkEtaPhi1D(photon.eta(), RecoDecay::constrainAngle(photon.phi())) >= cfgEMCalMapLevelBackground.value) {
        continue;
      }
      float energyCorrectionFactor = 1.f;
      energyCorrectionFactor = fEMCalCorrectionFactor->Eval(photon.e() > MinEnergy ? photon.e() : MinEnergy);
      ROOT::Math::PtEtaPhiMVector photon3(energyCorrectionFactor * photon.pt(), photon.eta(), photon.phi(), 0.);
      if (iCellIDPhoton1 >= 0) {
        if (std::fabs((photon1.E() - photon3.E()) / (photon1.E() + photon3.E()) < cfgMaxAsymmetry)) { // only use symmetric decays
          ROOT::Math::PtEtaPhiMVector mother1 = photon1 + photon3;
          float openingAngle1 = std::acos(photon1.Vect().Dot(photon3.Vect()) / (photon1.P() * photon3.P()));

          if (openingAngle1 > mesonConfig.minOpenAngle && thnConfigAxisInvMass.value[1] <= mother1.M() && thnConfigAxisInvMass.value.back() >= mother1.M() && thnConfigAxisPt.value[1] <= mother1.Pt() && thnConfigAxisPt.value.back() >= mother1.Pt()) {
            if (mesonConfig.enableTanThetadPhi) {
              float dTheta = photon1.Theta() - photon3.Theta();
              float dPhi = photon1.Phi() - photon3.Phi();
              if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
                registry.fill(HIST("hSparseCalibBack"), mother1.M(), mother1.E() / 2., cent);
              }
            } else {
              registry.fill(HIST("hSparseCalibBack"), mother1.M(), mother1.E() / 2., cent);
            }
          }
        }
      }
      if (iCellIDPhoton2 >= 0) {
        if (std::fabs((photon2.E() - photon3.E()) / (photon2.E() + photon3.E()) < cfgMaxAsymmetry)) { // only use symmetric decays
          ROOT::Math::PtEtaPhiMVector mother2 = photon2 + photon3;
          float openingAngle2 = std::acos(photon2.Vect().Dot(photon3.Vect()) / (photon2.P() * photon3.P()));

          if (openingAngle2 > mesonConfig.minOpenAngle && thnConfigAxisInvMass.value[1] <= mother2.M() && thnConfigAxisInvMass.value.back() >= mother2.M() && thnConfigAxisPt.value[1] <= mother2.Pt() && thnConfigAxisPt.value.back() >= mother2.Pt()) {
            if (mesonConfig.enableTanThetadPhi) {
              float dTheta = photon2.Theta() - photon3.Theta();
              float dPhi = photon2.Phi() - photon3.Phi();
              if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
                registry.fill(HIST("hSparseCalibBack"), mother2.M(), mother2.E() / 2., cent);
              }
            } else {
              registry.fill(HIST("hSparseCalibBack"), mother2.M(), mother2.E() / 2., cent);
            }
          }
        }
      }
    } // end of loop over third photon
    return;
  }

  /// \brief Calculate background using rotation background method for calib
  /// \param meson mother particle from photon1 & photon2 which defines rotation axis
  /// \param photon1 first photon (EMC) which will be rotated
  /// \param photon2 second photon (PCM) which will be rotated
  /// \param photonsEMC sub table of EMC photons of current event which will be combined with rotated photon2
  /// \param photonsPCM sub table of PCM photons of current event which will be combined with rotated photon1
  /// \param ig1 index of photon1
  /// \param ig2 index of photon2
  /// \param cent current collisions centrality
  template <typename TPhotonEMC, typename TPhotonPCM>
  void rotationBackgroundCalibEMCPCM(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photon1, ROOT::Math::PtEtaPhiMVector photon2, TPhotonEMC const& photonsEMC, TPhotonPCM const& photonsPCM, unsigned int ig1, unsigned int ig2, float cent)
  {
    // we need at least 2 clusters or 2 pcm photons for rotation
    if (photonsEMC.size() < NMinPhotonRotBkgMixed || photonsPCM.size() < NMinPhotonRotBkgMixed) {
      return;
    }
    int iCellIDPhoton1 = 0;
    int iCellIDPhoton2 = 0;

    ROOT::Math::AxisAngle rotationAxis(meson.Vect(), cfgRotAngle.value);
    ROOT::Math::Rotation3D rotationMatrix(rotationAxis);
    photon1 = rotationMatrix * photon1;
    photon2 = rotationMatrix * photon2;

    if (checkEtaPhi1D(photon1.Eta(), RecoDecay::constrainAngle(photon1.Phi())) >= cfgEMCalMapLevelBackground.value) {
      iCellIDPhoton1 = -1;
    }

    if (iCellIDPhoton1 == -1 && iCellIDPhoton2 == -1) {
      return;
    }
    // Combining with EMCal photons from event
    if (photonsEMC.size() >= NMinPhotonRotBkgMixed) {
      for (const auto& photon : photonsEMC) {
        if (photon.globalIndex() == ig1) {
          continue;
        }
        if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(photon))) {
          continue;
        }
        if (checkEtaPhi1D(photon.eta(), RecoDecay::constrainAngle(photon.phi())) >= cfgEMCalMapLevelBackground.value) {
          continue;
        }
        float energyCorrectionFactor = 1.f;
        energyCorrectionFactor = fEMCalCorrectionFactor->Eval(photon.e() > MinEnergy ? photon.e() : MinEnergy);
        ROOT::Math::PtEtaPhiMVector photon3(energyCorrectionFactor * photon.pt(), photon.eta(), photon.phi(), 0.);
        if (iCellIDPhoton2 >= 0) {
          ROOT::Math::PtEtaPhiMVector mother = photon2 + photon3;
          float openingAngle2 = std::acos(photon2.Vect().Dot(photon3.Vect()) / (photon2.P() * photon3.P()));

          if (openingAngle2 > mesonConfig.minOpenAngle && thnConfigAxisInvMass.value[1] <= mother.M() && thnConfigAxisInvMass.value.back() >= mother.M() && thnConfigAxisPt.value[1] <= mother.Pt() && thnConfigAxisPt.value.back() >= mother.Pt()) {
            if (mesonConfig.enableTanThetadPhi) {
              float dTheta = photon2.Theta() - photon3.Theta();
              float dPhi = photon2.Phi() - photon3.Phi();
              if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
                registry.fill(HIST("hSparseCalibBack"), mother.M(), mother.E() / 2., cent);
              }
            } else {
              registry.fill(HIST("hSparseCalibBack"), mother.M(), mother.E() / 2., cent);
            }
          }
        }
      } // end of loop over third photon
    } // check that we have at least 2 clusters

    // Combining with PCM photons from event
    if (photonsPCM.size() >= NMinPhotonRotBkgMixed) {
      for (const auto& photon : photonsPCM) {
        if (photon.globalIndex() == ig2) {
          continue;
        }
        if (!(fPCMPhotonCut.IsSelected<aod::V0Legs>(photon))) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector photon3(photon.pt(), photon.eta(), photon.phi(), 0.);
        if (iCellIDPhoton1 >= 0) {
          ROOT::Math::PtEtaPhiMVector mother = photon1 + photon3;
          float openingAngle2 = std::acos(photon1.Vect().Dot(photon3.Vect()) / (photon1.P() * photon3.P()));

          if (openingAngle2 > mesonConfig.minOpenAngle && thnConfigAxisInvMass.value[1] <= mother.M() && thnConfigAxisInvMass.value.back() >= mother.M() && thnConfigAxisPt.value[1] <= mother.Pt() && thnConfigAxisPt.value.back() >= mother.Pt()) {
            if (mesonConfig.enableTanThetadPhi) {
              float dTheta = photon1.Theta() - photon3.Theta();
              float dPhi = photon1.Phi() - photon3.Phi();
              if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
                registry.fill(HIST("hSparseCalibBack"), mother.M(), mother.E() / 2., cent);
              }
            } else {
              registry.fill(HIST("hSparseCalibBack"), mother.M(), mother.E() / 2., cent);
            }
          }
        }
      } // end of loop over third photon
    } // check that we have at least 2 clusters
    return;
  }

  // EMCal calibration same event
  void processEMCalCalib(Colls const& collisions, EMCalPhotons const& clusters, PCMPhotons const&)
  {
    float energyCorrectionFactor = 1.f;
    int nColl = 1;
    for (const auto& collision : collisions) {
      auto photonsPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());
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
      runNow = collision.runNumber();
      if (runNow != runBefore) {
        initCCDB(collision);
        runBefore = runNow;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&registry, collision);
      registry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // accepted
      registry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // accepted

      if (emccuts.cfgEnableQA) {
        for (const auto& photon : photonsPerCollision) {
          registry.fill(HIST("hEClusterBefore"), photon.e()); // before cuts
          if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(photon))) {
            continue;
          }
          registry.fill(HIST("hEClusterAfter"), photon.e()); // accepted after cuts
        }
      }
      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photonsPerCollision, photonsPerCollision))) {
        if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(g1)) || !(fEMCCut.IsSelected<EMCalPhotons::iterator>(g2))) {
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

        energyCorrectionFactor = fEMCalCorrectionFactor->Eval(g1.e() > MinEnergy ? g1.e() : MinEnergy);
        ROOT::Math::PtEtaPhiMVector v1(energyCorrectionFactor * g1.pt(), g1.eta(), g1.phi(), 0.);
        energyCorrectionFactor = fEMCalCorrectionFactor->Eval(g2.e() > MinEnergy ? g2.e() : MinEnergy);
        ROOT::Math::PtEtaPhiMVector v2(energyCorrectionFactor * g2.pt(), g2.eta(), g2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;

        float dTheta = v1.Theta() - v2.Theta();
        float dPhi = v1.Phi() - v2.Phi();
        float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));
        registry.fill(HIST("hClusterCuts"), 1);
        if (openingAngle <= mesonConfig.minOpenAngle) {
          registry.fill(HIST("hClusterCuts"), 2);
          continue;
        }
        if (cfgDoRotation) {
          if (nColl % cfgDownsampling.value == 0) {
            rotationBackgroundCalib<EMCalPhotons>(vMeson, v1, v2, photonsPerCollision, g1.globalIndex(), g2.globalIndex(), collision);
          }
        }
        if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
          registry.fill(HIST("hClusterCuts"), 3);
          continue;
        }
        if (thnConfigAxisPt.value[1] > vMeson.Pt() || thnConfigAxisPt.value.back() < vMeson.Pt()) {
          registry.fill(HIST("hClusterCuts"), 4);
          continue;
        }
        if (mesonConfig.cfgEnableQA) {
          registry.fill(HIST("hInvMassPt"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("hTanThetaPhi"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
          registry.fill(HIST("hAlphaPt"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
          registry.fill(HIST("hClusterCuts"), 5);
          continue;
        }
        if (std::fabs((v1.E() - v2.E()) / (v1.E() + v2.E()) < cfgMaxAsymmetry)) { // only use symmetric decays
          registry.fill(HIST("hClusterCuts"), 6);
          registry.fill(HIST("hSparseCalibSE"), vMeson.M(), vMeson.E() / 2., getCentrality(collision));
        }
      }
      if (cfgDoRotation) {
        if (nColl % cfgDownsampling.value == 0) {
          nColl = 1; // reset counter
        } else {
          nColl++;
        }
      }
    }
  }
  PROCESS_SWITCH(CalibTaskEmc, processEMCalCalib, "Process EMCal calibration same event", true);

  // EMCal calibration
  void processEMCalPCMCalib(Colls const& collisions, EMCalPhotons const& clusters, PCMPhotons const& photons, aod::V0Legs const&)
  {
    float energyCorrectionFactor = 1.f;
    int nColl = 1;
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
      runNow = collision.runNumber();
      if (runNow != runBefore) {
        initCCDB(collision);
        runBefore = runNow;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&registry, collision);
      registry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // accepted
      registry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // accepted

      auto photonsEMCPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());
      auto photonsPCMPerCollision = photons.sliceBy(perCollisionPCM, collision.globalIndex());

      if (emccuts.cfgEnableQA) {
        for (const auto& photon : photonsEMCPerCollision) {
          registry.fill(HIST("hEClusterBefore"), photon.e()); // before cuts
          if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(photon))) {
            continue;
          }
          registry.fill(HIST("hEClusterAfter"), photon.e()); // accepted after cuts
        }
      }
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photonsEMCPerCollision, photonsPCMPerCollision))) {
        if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(g1)) || !(fPCMPhotonCut.IsSelected<aod::V0Legs>(g2))) {
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
        registry.fill(HIST("hClusterCuts"), 1);
        if (openingAngle <= mesonConfig.minOpenAngle) {
          registry.fill(HIST("hClusterCuts"), 2);
          continue;
        }
        if (cfgDoRotation) {
          if (nColl % cfgDownsampling.value == 0) {
            rotationBackgroundCalibEMCPCM<EMCalPhotons, PCMPhotons>(vMeson, v1, v2, photonsEMCPerCollision, photonsPCMPerCollision, g1.globalIndex(), g2.globalIndex(), cent);
          }
        }
        if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
          registry.fill(HIST("hClusterCuts"), 3);
          continue;
        }
        if (thnConfigAxisPt.value[1] > vMeson.Pt() || thnConfigAxisPt.value.back() < vMeson.Pt()) {
          registry.fill(HIST("hClusterCuts"), 4);
          continue;
        }
        if (mesonConfig.cfgEnableQA) {
          registry.fill(HIST("hInvMassPt"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("hTanThetaPhi"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
          registry.fill(HIST("hAlphaPt"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
          registry.fill(HIST("hClusterCuts"), 5);
          continue;
        }
        registry.fill(HIST("hSparseCalibSE"), vMeson.M(), vMeson.E() / 2., getCentrality(collision));
      }
      if (cfgDoRotation) {
        if (nColl % cfgDownsampling.value == 0) {
          nColl = 1; // reset counter
        } else {
          nColl++;
        }
      }
    }
  }
  PROCESS_SWITCH(CalibTaskEmc, processEMCalPCMCalib, "Process EMCal calibration using PCM-EMC same event", true);

  // EMCal calibration mixed event
  void processEMCalCalibMixed(Colls const&, EMCalPhotons const&, PCMPhotons const&)
  {
    float energyCorrectionFactor = 1.f;

    SameKindPair<Colls, EMCalPhotons, BinningType> pair{binningOnPositions, mixingConfig.cfgMixingDepth, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

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
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(clusters1, clusters2))) {
        if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(g1)) || !(fEMCCut.IsSelected<EMCalPhotons::iterator>(g2))) {
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
        energyCorrectionFactor = fEMCalCorrectionFactor->Eval(g1.e() > MinEnergy ? g1.e() : MinEnergy);
        ROOT::Math::PtEtaPhiMVector v1(energyCorrectionFactor * g1.pt(), g1.eta(), g1.phi(), 0.);
        energyCorrectionFactor = fEMCalCorrectionFactor->Eval(g2.e() > MinEnergy ? g2.e() : MinEnergy);
        ROOT::Math::PtEtaPhiMVector v2(energyCorrectionFactor * g2.pt(), g2.eta(), g2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;

        float dTheta = v1.Theta() - v2.Theta();
        float dPhi = v1.Phi() - v2.Phi();
        float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));

        registry.fill(HIST("hClusterCutsMixed"), 1);
        if (openingAngle <= mesonConfig.minOpenAngle) {
          registry.fill(HIST("hClusterCutsMixed"), 2);
          continue;
        }
        if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
          registry.fill(HIST("hClusterCutsMixed"), 3);
          continue;
        }
        if (thnConfigAxisPt.value[1] > vMeson.Pt() || thnConfigAxisPt.value.back() < vMeson.Pt()) {
          registry.fill(HIST("hClusterCutsMixed"), 4);
          continue;
        }
        if (mesonConfig.cfgEnableQA) {
          registry.fill(HIST("hInvMassPtMixed"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("hTanThetaPhiMixed"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
          registry.fill(HIST("hAlphaPtMixed"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
          registry.fill(HIST("hClusterCutsMixed"), 5);
          continue;
        }
        if (std::fabs((v1.E() - v2.E()) / (v1.E() + v2.E()) < cfgMaxAsymmetry)) { // only use symmetric decays
          registry.fill(HIST("hClusterCutsMixed"), 6);
          registry.fill(HIST("hSparseCalibBack"), vMeson.M(), vMeson.E() / 2., getCentrality(c1));
        }
      }
    }
  }
  PROCESS_SWITCH(CalibTaskEmc, processEMCalCalibMixed, "Process EMCal calibration mixed event", false);

  // EMCal calibration
  void processEMCalPCMCalibMixed(Colls const&, EMCalPhotons const&, PCMPhotons const&, aod::V0Legs const&)
  {
    float energyCorrectionFactor = 1.f;

    LOG(info) << "Beginning of processEMCalPCMCalibMixed";

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
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photonEMC, photonPCM))) {
        if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(g1)) || !(fPCMPhotonCut.IsSelected<aod::V0Legs>(g2))) {
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

        registry.fill(HIST("hClusterCutsMixed"), 1);
        if (openingAngle <= mesonConfig.minOpenAngle) {
          registry.fill(HIST("hClusterCutsMixed"), 2);
          continue;
        }
        if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
          registry.fill(HIST("hClusterCutsMixed"), 3);
          continue;
        }
        if (thnConfigAxisPt.value[1] > vMeson.Pt() || thnConfigAxisPt.value.back() < vMeson.Pt()) {
          registry.fill(HIST("hClusterCutsMixed"), 4);
          continue;
        }
        if (mesonConfig.cfgEnableQA) {
          registry.fill(HIST("hInvMassPtMixed"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("hTanThetaPhiMixed"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
          registry.fill(HIST("hAlphaPtMixed"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
          registry.fill(HIST("hClusterCutsMixed"), 5);
          continue;
        }
        registry.fill(HIST("hSparseCalibBack"), vMeson.M(), vMeson.E() / 2., getCentrality(c1));
      }
    }
  }
  PROCESS_SWITCH(CalibTaskEmc, processEMCalPCMCalibMixed, "Process EMCal calibration with PCM-EMC mixed event", false);

}; // End struct CalibTaskEmc

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<CalibTaskEmc>(cfgc)};
}
