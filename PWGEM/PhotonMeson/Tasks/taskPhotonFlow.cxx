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

/// \file taskPhotonFlow.cxx
/// \brief Analysis task for photon flow with PCM or EMCal
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/PhotonMeson/Core/EMBitFlags.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/RecoDecay.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <EMCALBase/Geometry.h>
#include <EMCALBase/GeometryBase.h>
#include <EMCALCalib/BadChannelMap.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Concepts.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <string>
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
  FT0C = 2,
  TPCPos = 3,
  TPCNeg = 4,
  TPCTot = 5,
  FV0A = 6
};

enum CentralityEstimator {
  None = 0,
  CFT0A = 1,
  CFT0C = 2,
  CFT0M = 3,
  NCentralityEstimators = 4
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

struct TaskPhotonFlow {
  // configurable for flow
  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<int> qvecDetector{"qvecDetector", 0, "Detector for Q vector estimation (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5, FV0A: 6)"};
  Configurable<int> qvecSubADetector{"qvecSubADetector", 3, "Sub A Detector for Q vector estimation for resolution (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5, FV0A: 6)"};
  Configurable<int> qvecSubBDetector{"qvecSubBDetector", 4, "Sub B Detector for Q vector estimation for resolution (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5, FV0A: 6)"};
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int> cfgEMCalMapLevelSameEvent{"cfgEMCalMapLevelSameEvent", 4, "Different levels of correction for the same event, the smaller number includes the level of the higher number (4: none, 3: only inside EMCal, 2: remove edges, 1: exclude bad channels)"};
  Configurable<int> cfgDistanceToEdge{"cfgDistanceToEdge", 1, "Distance to edge in cells required for rotated cluster to be accepted"};
  Configurable<float> cfgMaxQVector{"cfgMaxQVector", 20.f, "Maximum allowed absolute QVector value."};

  // configurable axis
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0., 20.}, "pT axis for the neutral meson"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {20, 0., 100.}, "centrality axis for the current event"};
  ConfigurableAxis thnConfigAxisCosDeltaPhi{"thnConfigAxisCosDeltaPhi", {8, -1., 1.}, "cos(delta phi) axis for the current event"};
  ConfigurableAxis thnConfigAxisM02{"thnConfigAxisM02", {200, 0., 5.}, "M02 axis for the EMCal cluster"};
  ConfigurableAxis thnConfigAxisKappa{"thnConfigAxisKappa", {200, -10., 10.}, "Kappa axis for the EMCal cluster"};
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
    Configurable<bool> useEMCal{"useEMCal", false, "flag to use EMCal clusters"};
    Configurable<bool> useDCal{"useDCal", false, "flag to use DCal clusters"};
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
    std::string prefix = "correctionConfig";
    Configurable<std::string> cfgSpresoPath{"cfgSpresoPath", "Users/m/mhemmer/EM/Flow/Resolution", "Path to SP resolution file"};
    Configurable<bool> cfgApplySPresolution{"cfgApplySPresolution", false, "Apply resolution correction"};
  } correctionConfig;

  SliceCache cache;
  EventPlaneHelper epHelper;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb{};
  int runNow = 0;
  int runBefore = -1;

  static constexpr float MaxPhiEMCal = 3.9f;                                 // exatly the middle between EMCal and DCal
  static constexpr uint16_t MaxPhiEMCalUint = static_cast<uint16_t>(39000u); // exatly the middle between EMCal and DCal but as uint16_t that is used for storing phi values in derived data

  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::MinClusters, aod::NonLinEmcClusters>;
  using PCMPhotons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds, aod::NonLinV0s>;
  using Colls = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMEventsQvec_001, aod::EmEmcalObjects>;

  Partition<EMCalPhotons> emcalPhotons = aod::mincluster::storedPhi < MaxPhiEMCalUint;
  Partition<EMCalPhotons> dcalPhotons = aod::mincluster::storedPhi >= MaxPhiEMCalUint;

  static constexpr std::size_t NQVecEntries = 6;

  PresliceOptional<aod::EMCEMEventIds> perCollisionEMC = o2::aod::emccluster::pmeventId;
  PresliceOptional<aod::V0KFEMEventIds> perCollisionPCM = aod::v0photonkf::pmeventId;
  PresliceOptional<MinMTracks> perEMCClusterMT = o2::aod::mintm::minClusterId;
  PresliceOptional<MinMSTracks> perEMCClusterMS = o2::aod::mintm::minClusterId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  o2::emcal::Geometry* emcalGeom = nullptr;
  TH1D* h1SPResolution = nullptr;
  // Constants for eta and phi ranges for the look up table
  static constexpr double EtaMin = -0.75, etaMax = 0.75;
  static constexpr int NBinsEta = 150; // 150 bins for eta

  static constexpr double PhiMin = 1.35, phiMax = 5.75;
  static constexpr int NBinsPhi = 440; // (440 bins = 0.01 step size covering most regions)

  std::array<int8_t, NBinsEta * NBinsPhi> lookupTable1D{};
  float epsilon = 1.e-8;

  // To access the 1D array
  static inline int getIndex(int iEta, int iPhi)
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

    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality (%)"};
    const AxisSpec thnAxisCosDeltaPhi{thnConfigAxisCosDeltaPhi, std::format("cos({}(#varphi - #Psi_{{sub}}))", harmonic.value)};
    const AxisSpec thnAxisM02{thnConfigAxisM02, "M_{02}"};
    const AxisSpec thnAxisKappa{thnConfigAxisKappa, "#Kappa"};
    const AxisSpec thAxisClusterEnergy{thnConfigAxisPt, "#it{E} (GeV)"};
    const AxisSpec thAxisEnergyCalib{thnConfigAxisEnergyCalib, "#it{E}_{clus} (GeV)"};
    const AxisSpec thAxisEnergy{1000, 0., 100., "#it{E}_{clus} (GeV)"};
    const AxisSpec thAxisEta{320, -0.8, 0.8, "#eta"};
    const AxisSpec thAxisPhi{500, 0, 2 * 3.14159, "phi"};

    if (emccuts.cfgEnableQA.value) {
      registry.add("clusterQA/hEClusterBefore", "Histo for cluster energy before cuts", HistType::kTH1D, {thAxisClusterEnergy});
      registry.add("clusterQA/hEClusterAfter", "Histo for cluster energy after cuts", HistType::kTH1D, {thAxisClusterEnergy});
      registry.add("clusterQA/hClusterEtaPhiBefore", "hClusterEtaPhiBefore", HistType::kTH2D, {thAxisPhi, thAxisEta});
      registry.add("clusterQA/hClusterEtaPhiAfter", "hClusterEtaPhiAfter", HistType::kTH2D, {thAxisPhi, thAxisEta});
    }

    if (doprocessEMC) {
      registry.add("EMCal/hSparsePhotonFlow", "<v_n> vs M_{02} vs p_T vs cent", HistType::kTProfile3D, {thnAxisM02, thnAxisPt, thnAxisCent});
      registry.add("EMCal/hSparsePhoton", "M_{02} vs p_T vs cent", HistType::kTH3D, {thnAxisM02, thnAxisPt, thnAxisCent});
    }
    if (doprocessPCM) {
      registry.add("PCM/hSparsePhotonFlow", "<v_n> vs M_{02} vs p_T vs cent", HistType::kTProfile3D, {thnAxisM02, thnAxisPt, thnAxisCent});
      registry.add("PCM/hSparsePhoton", "M_{02} vs p_T vs cent", HistType::kTH3D, {thnAxisM02, thnAxisPt, thnAxisCent});
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    LOG(info) << "thnConfigAxisPt.value[1] = " << thnConfigAxisPt.value[1] << " thnConfigAxisPt.value.back() = " << thnConfigAxisPt.value.back();

  }; // end init

  /// \brief Check whether a photon (by its phi) falls in the EMCal or DCal acceptance
  /// \param phi azimuthal angle of the photon
  static bool isEMCalRegion(float phi)
  {
    return phi < MaxPhiEMCal;
  }

  /// Change radians to degree
  /// \param angle in radians
  /// \return angle in degree
  static constexpr float getAngleDegree(float angle)
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
    } else if ((emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::EMCAL_THIRD) ||
               (emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::DCAL_EXT)) {
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

  template <o2::soa::is_iterator TCollision>
  bool isCellMasked(int cellID, TCollision const& collision)
  {
    bool masked = false;
    auto maskStatus = collision.badChannelMap().getChannelStatus(cellID);
    masked = (maskStatus != o2::emcal::BadChannelMap::MaskType_t::GOOD_CELL);
    return masked;
  }

  template <o2::soa::is_iterator TLeg>
  float getV0Kappa(TLeg const& pos, TLeg const& ele)
  {
    float kappa = (std::fabs(pos.tpcNSigmaEl()) + std::fabs(ele.tpcNSigmaEl())) / 2.f + (2.f * pos.tpcNSigmaEl() + ele.tpcNSigmaEl());
    return kappa;
  }

  template <o2::soa::is_iterator TCollision>
  void initCCDB(TCollision const& collision)
  {
    // Load EMCal geometry
    emcalGeom = o2::emcal::Geometry::GetInstanceFromRunNumber(collision.runNumber());
    lookupTable1D.fill(-1);
    double binWidthEta = (etaMax - EtaMin) / NBinsEta;
    double binWidthPhi = (phiMax - PhiMin) / NBinsPhi;

    if (cfgEMCalMapLevelSameEvent >= static_cast<int>(MapLevel::kAll)) {
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
            } else if (isCellMasked(cellID, collision)) {
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

  /// \brief check if standard event cuts + FT0 occupancy + centrality + QVec good is
  /// \param collision collision that will be checked
  /// \return true if collision survives all checks, otherwise false
  template <o2::soa::is_iterator TCollision>
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

  void processEMC(Colls const& collisions, EMCalPhotons const& clusters, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds)
  {
    if (clusters.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    EMBitFlags emcFlags(clusters.size());
    fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds);

    for (const auto& collision : collisions) {
      float cent = getCentrality(collision);
      if (!isFullEventSelected(collision, true)) {
        continue;
      }
      runNow = collision.runNumber();
      if (runNow != runBefore) {
        initCCDB(collision);
        runBefore = runNow;
      }
      auto photonsPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());

      for (const auto& photon : photonsPerCollision) {
        if (emccuts.cfgEnableQA.value) {
          registry.fill(HIST("clusterQA/hEClusterBefore"), photon.corrE());                  // before cuts
          registry.fill(HIST("clusterQA/hClusterEtaPhiBefore"), photon.phi(), photon.eta()); // before cuts
        }
        if (!(emcFlags.test(photon.globalIndex()))) {
          continue;
        }
        if ((emccuts.useEMCal.value && !emccuts.useDCal.value && (!isEMCalRegion(photon.phi()))) || (!emccuts.useEMCal.value && emccuts.useDCal.value && (isEMCalRegion(photon.phi())))) {
          continue;
        }
        if (cfgDistanceToEdge.value > 0 && (checkEtaPhi1D(photon.eta(), RecoDecay::constrainAngle(photon.phi())) >= cfgEMCalMapLevelSameEvent.value)) {
          continue;
        }
        if (emccuts.cfgEnableQA.value) {
          registry.fill(HIST("clusterQA/hEClusterAfter"), photon.corrE());                  // accepted after cuts
          registry.fill(HIST("clusterQA/hClusterEtaPhiAfter"), photon.phi(), photon.eta()); // after cuts
        }

        auto [xQVec, yQVec] = getQvec(collision, qvecDetector);

        float phiCand = photon.phi();

        float cosNPhi = std::cos(harmonic * phiCand);
        float sinNPhi = std::sin(harmonic * phiCand);
        float scalprodCand = cosNPhi * xQVec + sinNPhi * yQVec;

        if (correctionConfig.cfgApplySPresolution.value) {
          scalprodCand = scalprodCand / h1SPResolution->GetBinContent(h1SPResolution->FindBin(cent + epsilon));
        }
        registry.fill(HIST("EMCal/hSparsePhotonFlow"), photon.m02(), photon.corrPt(), cent, scalprodCand);
        registry.fill(HIST("EMCal/hSparsePhoton"), photon.m02(), photon.corrPt(), cent);
      } // end of loop over single cluster
    } // end of loop over collisions
  } // processEMC
  PROCESS_SWITCH(TaskPhotonFlow, processEMC, "Process single EMCal clusters as function of M02", false);

  void processPCM(Colls const& collisions, PCMPhotons const& photons, aod::V0Legs const& legs)
  {
    if (photons.size() <= 0 || legs.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    auto ele = legs.begin();
    auto pos = legs.begin();
    EMBitFlags v0flags(photons.size());
    fV0PhotonCut.AreSelectedRunning<decltype(photons), aod::V0Legs>(v0flags, photons, &registry);
    for (const auto& collision : collisions) {

      if (!isFullEventSelected(collision, true)) {
        continue;
      }

      float cent = getCentrality(collision);

      runNow = collision.runNumber();
      if (runNow != runBefore) {
        initCCDB(collision);
        runBefore = runNow;
      }

      auto photonsPerCollision = photons.sliceBy(perCollisionPCM, collision.globalIndex());
      for (const auto& photon : photonsPerCollision) {
        if (!(v0flags.test(photon.globalIndex()))) {
          continue;
        }
        ele.setCursor(photon.negTrackId());
        pos.setCursor(photon.posTrackId());

        auto [xQVec, yQVec] = getQvec(collision, qvecDetector);

        float phiCand = photon.phi();

        float cosNPhi = std::cos(harmonic * phiCand);
        float sinNPhi = std::sin(harmonic * phiCand);
        float scalprodCand = cosNPhi * xQVec + sinNPhi * yQVec;

        if (correctionConfig.cfgApplySPresolution.value) {
          scalprodCand = scalprodCand / h1SPResolution->GetBinContent(h1SPResolution->FindBin(cent + epsilon));
        }

        float kappa = getV0Kappa(pos, ele);
        registry.fill(HIST("PCM/hSparsePhotonFlow"), kappa, photon.corrPt(), cent, scalprodCand);
        registry.fill(HIST("PCM/hSparsePhoton"), kappa, photon.corrPt(), cent);
      }
    } // end of loop over collisions
  }
  PROCESS_SWITCH(TaskPhotonFlow, processPCM, "Process PCM Pi0 candidates", true);

}; // End struct TaskPhotonFlow

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<TaskPhotonFlow>(context)};
}
