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

#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
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
#include <cstdint>
#include <memory>
#include <numbers>
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
  TPCTot
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
  Configurable<int> qvecDetector{"qvecDetector", 0, "Detector for Q vector estimation (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5)"};
  Configurable<int> qvecSubADetector{"qvecSubADetector", 3, "Sub A Detector for Q vector estimation for resolution (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5)"};
  Configurable<int> qvecSubBDetector{"qvecSubBDetector", 4, "Sub B Detector for Q vector estimation for resolution (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5)"};
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3)"};
  Configurable<bool> saveEpResoHisto{"saveEpResoHisto", false, "Flag to save event plane resolution histogram"};
  Configurable<bool> saveSPResoHist{"saveSPResoHist", false, "Flag to save scalar product resolution histogram"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> cfgDoRotation{"cfgDoRotation", true, "Flag to enable rotation background method"};
  Configurable<int> cfgDownsampling{"cfgDownsampling", 1, "Calculate rotation background only for every <value> collision"};
  Configurable<int> cfgEMCalMapLevelBackground{"cfgEMCalMapLevelBackground", 2, "Different levels of correction for the background, the smaller number includes the level of the higher number (4: none, 3: only inside EMCal, 2: exclude bad channels, 1: remove edges)"};
  Configurable<int> cfgEMCalMapLevelSameEvent{"cfgEMCalMapLevelSameEvent", 1, "Different levels of correction for the same event, the smaller number includes the level of the higher number (4: none, 3: only inside EMCal, 2: exclude bad channels, 1: remove edges)"};
  Configurable<float> cfgRotAngle{"cfgRotAngle", std::move(const_cast<float&>(o2::constants::math::PIHalf)), "Angle used for the rotation method"};
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
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -1, "min. track occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -1, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<float> cfgMinCent{"cfgMinCent", 0, "min. centrality (%)"};
    Configurable<float> cfgMaxCent{"cfgMaxCent", 90, "max. centrality (%)"};
    Configurable<bool> onlyKeepWeightedEvents{"onlyKeepWeightedEvents", false, "flag to keep only weighted events (for JJ MCs) and remove all MB events (with weight = 1)"};
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
    Configurable<float> cfgEMCEoverp{"cfgEMCEoverp", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
    Configurable<bool> cfgEMCUseExoticCut{"cfgEMCUseExoticCut", true, "FLag to use the EMCal exotic cluster cut"};
    Configurable<bool> cfgEMCUseTM{"cfgEMCUseTM", false, "flag to use EMCal track matching cut or not"};
    Configurable<bool> cfgEnableQA{"cfgEnableQA", false, "flag to turn QA plots on/off"};
  } emccuts;

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

  Filter clusterFilter = aod::skimmedcluster::time >= emccuts.cfgEMCminTime && aod::skimmedcluster::time <= emccuts.cfgEMCmaxTime && aod::skimmedcluster::m02 >= emccuts.cfgEMCminM02 && aod::skimmedcluster::m02 <= emccuts.cfgEMCmaxM02 && aod::skimmedcluster::e >= emccuts.cfgEMCminE;
  Filter collisionFilter = (nabs(aod::collision::posZ) <= eventcuts.cfgZvtxMax) && (aod::evsel::ft0cOccupancyInTimeRange <= eventcuts.cfgFT0COccupancyMax) && (aod::evsel::ft0cOccupancyInTimeRange >= eventcuts.cfgFT0COccupancyMin);
  using FilteredEMCalPhotons = soa::Filtered<soa::Join<aod::EMCEMEventIds, aod::SkimEMCClusters>>;
  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::SkimEMCClusters>;
  using FilteredCollsWithQvecs = soa::Filtered<soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>>;
  using CollsWithQvecs = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
  using Colls = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent>;

  Preslice<EMCalPhotons> perCollisionEMC = aod::emccluster::emeventId;

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
  static constexpr int64_t NMinPhotonRotBkg = 3;

  // Usage when cfgEnableNonLin is enabled
  std::unique_ptr<TF1> fEMCalCorrectionFactor; // ("fEMCalCorrectionFactor","(1 + [0]/x + [1]/x^2) / (1 + [2]/x)", 0.3, 100.);

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
    const float a = emccuts.cfgEMCTMEta->at(0);
    const float b = emccuts.cfgEMCTMEta->at(1);
    const float c = emccuts.cfgEMCTMEta->at(2);

    const float d = emccuts.cfgEMCTMPhi->at(0);
    const float e = emccuts.cfgEMCTMPhi->at(1);
    const float f = emccuts.cfgEMCTMPhi->at(2);
    LOGF(info, "EMCal track matching parameters : a = %f, b = %f, c = %f, d = %f, e = %f, f = %f", a, b, c, d, e, f);
    fEMCCut.SetTrackMatchingEta([a, b, c](float pT) { return a + std::pow(pT + b, c); });
    fEMCCut.SetTrackMatchingPhi([d, e, f](float pT) { return d + std::pow(pT + e, f); });
    fEMCCut.SetMinEoverP(emccuts.cfgEMCEoverp);

    fEMCCut.SetMinE(emccuts.cfgEMCminE);
    fEMCCut.SetMinNCell(emccuts.cfgEMCminNCell);
    fEMCCut.SetM02Range(emccuts.cfgEMCminM02, emccuts.cfgEMCmaxM02);
    fEMCCut.SetTimeRange(emccuts.cfgEMCminTime, emccuts.cfgEMCmaxTime);
    fEMCCut.SetUseExoticCut(emccuts.cfgEMCUseExoticCut);
    fEMCCut.SetClusterizer(emccuts.clusterDefinition);
  }

  void init(InitContext&)
  {
    if (harmonic != kElliptic && harmonic != kTriangluar) {
      LOG(info) << "Harmonic was set to " << harmonic << " but can only be 2 or 3!";
    }

    defineEMEventCut();
    defineEMCCut();
    fEMCCut.SetUseTM(emccuts.cfgEMCUseTM); // disables TM
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&registry);

    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality (%)"};
    const AxisSpec thnAxisCosNPhi{thnConfigAxisCosNPhi, Form("cos(%d#varphi)", harmonic.value)};
    const AxisSpec thnAxisCosDeltaPhi{thnConfigAxisCosDeltaPhi, Form("cos(%d(#varphi - #Psi_{sub}))", harmonic.value)};
    const AxisSpec thnAxisM02{thnConfigAxisM02, "M_{02}"};
    const AxisSpec thAxisTanThetaPhi{mesonConfig.thConfigAxisTanThetaPhi, "atan(#Delta#theta/#Delta#varphi)"};
    const AxisSpec thAxisClusterEnergy{thnConfigAxisPt, "#it{E} (GeV)"};
    const AxisSpec thAxisEnergyCalib{thnConfigAxisEnergyCalib, "#it{E}_{clus} (GeV)"};
    const AxisSpec thAxisAlpha{100, -1., +1, "#alpha"};
    const AxisSpec thAxisMult{1000, 0., +1000, "#it{N}_{ch}"};
    const AxisSpec thAxisEnergy{1000, 0., 100., "#it{E}_{clus} (GeV)"};
    const AxisSpec thAxisTime{1500, -600, 900, "#it{t}_{cl} (ns)"};
    const AxisSpec thAxisEta{320, -0.8, 0.8, "#eta"};
    const AxisSpec thAxisPhi{500, 0, 2 * 3.14159, "phi"};
    const AxisSpec thAxisNCell{17664, 0.5, +17664.5, "#it{N}_{cell}"};
    const AxisSpec thAxisPsi{360 / harmonic.value, -(1. / static_cast<float>(harmonic.value)) * std::numbers::pi_v<float>, (1. / static_cast<float>(harmonic.value)) * std::numbers::pi_v<float>, Form("#Psi_{%d}", harmonic.value)};
    const AxisSpec thAxisCN{8, 0.5, 8.5, "#it{c}_{n}"};
    const AxisSpec thAxisSN{8, 0.5, 8.5, "#it{s}_{n}"};
    const AxisSpec thAxisCPUTime{1000, 0, 10000, "#it{t} (#mus)"};
    const AxisSpec thAxisAzimuth{360, -o2::constants::math::PI, o2::constants::math::PI, "#it{#varphi} (rad)"};

    const AxisSpec thnAxisMixingVtx{mixingConfig.cfgVtxBins, "#it{z} (cm)"};
    const AxisSpec thnAxisMixingCent{mixingConfig.cfgCentBins, "Centrality (%)"};
    const AxisSpec thnAxisMixingEP{mixingConfig.cfgEPBins, Form("cos(%d#varphi)", harmonic.value)};

    registry.add("hSparsePi0Flow", "<v_n> vs m_{inv} vs p_T vs cent for same event", HistType::kTProfile3D, {thnAxisInvMass, thnAxisPt, thnAxisCent});
    registry.add("hSparsePi0", "m_{inv} vs p_T vs cent for same event", HistType::kTH3D, {thnAxisInvMass, thnAxisPt, thnAxisCent});

    registry.add("hSparseBkgMixFlow", "<v_n> vs m_{inv} vs p_T vs cent for mixed event", HistType::kTProfile3D, {thnAxisInvMass, thnAxisPt, thnAxisCent});
    registry.add("hSparseBkgMix", "m_{inv} vs p_T vs cent for mixed event", HistType::kTH3D, {thnAxisInvMass, thnAxisPt, thnAxisCent});

    if (cfgDoRotation.value) {
      registry.add("hSparseBkgRotFlow", "<v_n> vs m_{inv} vs p_T vs cent for rotation background", HistType::kTProfile3D, {thnAxisInvMass, thnAxisPt, thnAxisCent});
      registry.add("hSparseBkgRot", "m_{inv} vs p_T vs cent for rotation background", HistType::kTH3D, {thnAxisInvMass, thnAxisPt, thnAxisCent});
    }

    registry.add("h3DMixingCount", "THn Event Mixing QA", HistType::kTH3D, {thnAxisMixingVtx, thnAxisMixingCent, thnAxisMixingEP});
    if (cfgDoPlaneQA.value) {
      registry.add("hSparsePi0FlowPlane", "THn for SP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCent, thnAxisCosDeltaPhi});
    }
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

    if (saveSPResoHist.value) {
      registry.add("spReso/hSpResoFT0cFT0a", "hSpResoFT0cFT0a; centrality; Q_{FT0c} #bullet Q_{FT0a}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0cTPCpos", "hSpResoFT0cTPCpos; centrality; Q_{FT0c} #bullet Q_{TPCpos}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0cTPCneg", "hSpResoFT0cTPCneg; centrality; Q_{FT0c} #bullet Q_{TPCneg}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0cTPCtot", "hSpResoFT0cTPCtot; centrality; Q_{FT0c} #bullet Q_{TPCtot}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0aTPCpos", "hSpResoFT0aTPCpos; centrality; Q_{FT0a} #bullet Q_{TPCpos}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0aTPCneg", "hSpResoFT0aTPCneg; centrality; Q_{FT0a} #bullet Q_{TPCneg}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0aTPCtot", "hSpResoFT0aTPCtot; centrality; Q_{FT0m} #bullet Q_{TPCtot}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0mTPCpos", "hSpResoFT0mTPCpos; centrality; Q_{FT0m} #bullet Q_{TPCpos}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0mTPCneg", "hSpResoFT0mTPCneg; centrality; Q_{FT0m} #bullet Q_{TPCneg}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0mTPCtot", "hSpResoFT0mTPCtot; centrality; Q_{FT0m} #bullet Q_{TPCtot}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoTPCposTPCneg", "hSpResoTPCposTPCneg; centrality; Q_{TPCpos} #bullet Q_{TPCneg}", HistType::kTProfile, {thnAxisCent});
    }

    if (saveEpResoHisto.value) {
      registry.add("hEventPlaneAngleFT0M", "hEventPlaneAngleFT0M", HistType::kTH2D, {thnAxisCent, thAxisPsi});
      registry.add("hEventPlaneAngleTPCpos", "hEventPlaneAngleTPCpos", HistType::kTH2D, {thnAxisCent, thAxisPsi});
      registry.add("hEventPlaneAngleTPCneg", "hEventPlaneAngleTPCneg", HistType::kTH2D, {thnAxisCent, thAxisPsi});
      registry.add("epReso/hEpResoFT0cFT0a", "hEpResoFT0cFT0a; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0cTPCpos", "hEpResoFT0cTPCpos; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0cTPCneg", "hEpResoFT0cTPCneg; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0cTPCtot", "hEpResoFT0cTPCtot; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0aTPCpos", "hEpResoFT0aTPCpos; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0aTPCneg", "hEpResoFT0aTPCneg; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0aTPCtot", "hEpResoFT0aTPCtot; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0mTPCpos", "hEpResoFT0mTPCpos; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0mTPCneg", "hEpResoFT0mTPCneg; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0mTPCtot", "hEpResoFT0mTPCtot; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoTPCposTPCneg", "hEpResoTPCposTPCneg; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpCosCoefficientsFT0c", "hEpCosCoefficientsFT0c; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsFT0c", "hEpSinCoefficientsFT0c; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("epReso/hEpCosCoefficientsFT0a", "hEpCosCoefficientsFT0a; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsFT0a", "hEpSinCoefficientsFT0a; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("epReso/hEpCosCoefficientsFT0m", "hEpCosCoefficientsFT0m; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsFT0m", "hEpSinCoefficientsFT0m; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("epReso/hEpCosCoefficientsTPCpos", "hEpCosCoefficientsTPCpos; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsTPCpos", "hEpSinCoefficientsTPCpos; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("epReso/hEpCosCoefficientsTPCneg", "hEpCosCoefficientsTPCneg; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsTPCneg", "hEpSinCoefficientsTPCneg; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("epReso/hEpCosCoefficientsTPCTots", "hEpCosCoefficientsTPCTots; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsTPCTots", "hEpSinCoefficientsTPCTots; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("QVector/hQVecMeanRVsPhiFT0a", "hQVecMeanRVsPhiFT0a; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
      registry.add("QVector/hQVecMeanRVsPhiFT0c", "hQVecMeanRVsPhiFT0c; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
      registry.add("QVector/hQVecMeanRVsPhiFT0m", "hQVecMeanRVsPhiFT0m; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
      registry.add("QVector/hQVecMeanRVsPhiTPCpos", "hQVecMeanRVsPhiTPCpos; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
      registry.add("QVector/hQVecMeanRVsPhiTPCneg", "hQVecMeanRVsPhiTPCneg; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
      registry.add("QVector/hQVecMeanRVsPhiTPCTot", "hQVecMeanRVsPhiTPCTot; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
    }

    if (emccuts.cfgEnableQA.value) {
      registry.add("clusterQA/hEClusterBefore", "Histo for cluster energy before cuts", HistType::kTH1D, {thAxisClusterEnergy});
      registry.add("clusterQA/hEClusterAfter", "Histo for cluster energy after cuts", HistType::kTH1D, {thAxisClusterEnergy});
      registry.add("clusterQA/hClusterEtaPhiBefore", "hClusterEtaPhiBefore", HistType::kTH2D, {thAxisPhi, thAxisEta});
      registry.add("clusterQA/hClusterEtaPhiAfter", "hClusterEtaPhiAfter", HistType::kTH2D, {thAxisPhi, thAxisEta});
      if (cfgDoRotation.value) {
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
  float getAngleDegree(float angle)
  {
    return angle * 180.f * std::numbers::inv_pi_v<float>;
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

  /// Get all used Q vector
  /// \param collision is the collision with the Q vector information
  template <typename TCollision>
  std::vector<float> getAllQvec(TCollision const& collision)
  {
    // Retrieve the Q vectors using the helper function for each detector
    auto [xQVecMain, yQVecMain] = getQvec(collision, qvecDetector);
    auto [xQVecSubA, yQVecSubA] = getQvec(collision, qvecSubADetector);
    auto [xQVecSubB, yQVecSubB] = getQvec(collision, qvecSubBDetector);

    return {xQVecMain, yQVecMain, xQVecSubA, yQVecSubA, xQVecSubB, yQVecSubB};
  }

  /// Get the Q vector
  /// \param collision is the collision with the Q vector information
  template <typename TCollision>
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
  bool isQvecGood(std::vector<float> const& QVecs)
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

  template <typename TCollision>
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
  template <typename TPhotons, typename TCollision>
  void rotationBackground(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photon1, ROOT::Math::PtEtaPhiMVector photon2, TPhotons const& photons_coll, unsigned int ig1, unsigned int ig2, TCollision const& collision)
  {
    // if less than 3 clusters are present skip event since we need at least 3 clusters
    if (photons_coll.size() < NMinPhotonRotBkg) {
      return;
    }

    auto [xQVec, yQVec] = getQvec(collision, qvecDetector);
    float cent = getCentrality(collision);
    int iCellIDPhoton1 = 0;
    int iCellIDPhoton2 = 0;

    ROOT::Math::AxisAngle rotationAxis(meson.Vect(), cfgRotAngle.value);
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
            if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
              registry.fill(HIST("hSparseBkgRotFlow"), mother1.M(), mother1.Pt(), cent, scalprodCand1);
              registry.fill(HIST("hSparseBkgRot"), mother1.M(), mother1.Pt(), cent);
            }
          } else {
            registry.fill(HIST("hSparseBkgRotFlow"), mother1.M(), mother1.Pt(), cent, scalprodCand1);
            registry.fill(HIST("hSparseBkgRot"), mother1.M(), mother1.Pt(), cent);
          }
        }
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
            if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
              registry.fill(HIST("hSparseBkgRotFlow"), mother2.M(), mother2.Pt(), cent, scalprodCand2);
              registry.fill(HIST("hSparseBkgRot"), mother2.M(), mother2.Pt(), cent);
            }
          } else {
            registry.fill(HIST("hSparseBkgRotFlow"), mother2.M(), mother2.Pt(), cent, scalprodCand2);
            registry.fill(HIST("hSparseBkgRot"), mother2.M(), mother2.Pt(), cent);
          }
        }
      }
    } // end of loop over third photon
    return;
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
      if (correctionConfig.cfgEnableNonLin.value) {
        energyCorrectionFactor = fEMCalCorrectionFactor->Eval(photon.e() > MinEnergy ? photon.e() : MinEnergy);
      }
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

  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param meson are the selected candidates
  template <const int histType, typename TCollision>
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

  // Pi0 from EMCal
  void processEMCal(CollsWithQvecs const& collisions, EMCalPhotons const& clusters)
  {
    int nColl = 1;
    float energyCorrectionFactor = 1.f;
    if (cfgDoReverseScaling.value) {
      energyCorrectionFactor = 1.0505f;
    }
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

      if (emccuts.cfgEnableQA.value) {
        for (const auto& photon : photonsPerCollision) {
          registry.fill(HIST("clusterQA/hEClusterBefore"), photon.e());                      // before cuts
          registry.fill(HIST("clusterQA/hClusterEtaPhiBefore"), photon.phi(), photon.eta()); // before cuts
          if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(photon))) {
            continue;
          }
          if (cfgDistanceToEdge.value && (checkEtaPhi1D(photon.eta(), RecoDecay::constrainAngle(photon.phi())) >= cfgEMCalMapLevelSameEvent.value)) {
            continue;
          }
          registry.fill(HIST("clusterQA/hEClusterAfter"), photon.e());                      // accepted after cuts
          registry.fill(HIST("clusterQA/hClusterEtaPhiAfter"), photon.phi(), photon.eta()); // after cuts
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
        registry.fill(HIST("hClusterCuts"), 1);
        if (openingAngle <= mesonConfig.minOpenAngle) {
          registry.fill(HIST("hClusterCuts"), 2);
          continue;
        }
        if (cfgDoRotation) {
          if (nColl % cfgDownsampling.value == 0) {
            rotationBackground<EMCalPhotons>(vMeson, v1, v2, photonsPerCollision, g1.globalIndex(), g2.globalIndex(), collision);
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
        if (mesonConfig.cfgEnableQA.value) {
          registry.fill(HIST("mesonQA/hInvMassPt"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("mesonQA/hTanThetaPhi"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
          registry.fill(HIST("mesonQA/hAlphaPt"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
          registry.fill(HIST("hClusterCuts"), 5);
          continue;
        }
        registry.fill(HIST("hClusterCuts"), 6);
        runFlowAnalysis<0>(collision, vMeson);
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
  PROCESS_SWITCH(TaskPi0FlowEMC, processEMCal, "Process EMCal Pi0 candidates", true);

  // Pi0 from EMCal
  void processEMCalMixed(FilteredCollsWithQvecs const& collisions, FilteredEMCalPhotons const& clusters)
  {
    float energyCorrectionFactor = 1.f;
    if (cfgDoReverseScaling.value) {
      energyCorrectionFactor = 1.0505f;
    }
    auto getClustersSize =
      [&clusters, this](FilteredCollsWithQvecs::iterator const& col) {
        auto associatedClusters = clusters.sliceByCached(emccluster::emeventId, col.globalIndex(), this->cache); // it's cached, so slicing/grouping happens only once
        return associatedClusters.size();
      };

    using BinningType = FlexibleBinningPolicy<std::tuple<decltype(getClustersSize)>, aod::collision::PosZ, aod::cent::CentFT0C, emevent::EP2FT0M<emevent::Q2xFT0M, emevent::Q2yFT0M>>;
    BinningType binningWithLambda{{getClustersSize}, {mixingConfig.cfgVtxBins, mixingConfig.cfgCentBins, mixingConfig.cfgEPBins}, true};

    auto clustersTuple = std::make_tuple(clusters);
    SameKindPair<FilteredCollsWithQvecs, FilteredEMCalPhotons, BinningType> pair{binningWithLambda, mixingConfig.cfgMixingDepth, -1, collisions, clustersTuple, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

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
        if (mesonConfig.cfgEnableQA.value) {
          registry.fill(HIST("mesonQA/hInvMassPtMixed"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("mesonQA/hTanThetaPhiMixed"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
          registry.fill(HIST("mesonQA/hAlphaPtMixed"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        if (mesonConfig.enableTanThetadPhi && mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(std::atan(dTheta / dPhi)))) {
          registry.fill(HIST("hClusterCutsMixed"), 5);
          continue;
        }
        registry.fill(HIST("hClusterCutsMixed"), 6);
        runFlowAnalysis<2>(c1, vMeson);
      }
    }
  }
  PROCESS_SWITCH(TaskPi0FlowEMC, processEMCalMixed, "Process EMCal Pi0 mixed event candidates", false);

  // Resolution
  void processResolution(CollsWithQvecs::iterator const& collision)
  {
    o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&registry, collision);
    if (!(fEMEventCut.IsSelected(collision))) {
      // no selection on the centrality is applied on purpose to allow for the resolution study in post-processing
      return;
    }
    if (!(eventcuts.cfgFT0COccupancyMin <= collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax)) {
      return;
    }
    float cent = getCentrality(collision);
    if (cent < eventcuts.cfgMinCent || cent > eventcuts.cfgMaxCent) {
      // event selection
      return;
    }
    if (!isQvecGood(getAllQvec(collision))) {
      // selection based on QVector
      return;
    }
    o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&registry, collision);
    registry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // accepted
    registry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // accepted

    float centrality = getCentrality(collision); // centrality not updated in the rejection mask function
    float xQVecFT0a = -999.f;
    float yQVecFT0a = -999.f;
    float xQVecFT0c = -999.f;
    float yQVecFT0c = -999.f;
    float xQVecFT0m = -999.f;
    float yQVecFT0m = -999.f;
    float xQVecBPos = -999.f;
    float yQVecBPos = -999.f;
    float xQVecBNeg = -999.f;
    float yQVecBNeg = -999.f;
    float xQVecBTot = -999.f;
    float yQVecBTot = -999.f;
    if (harmonic == kElliptic) {
      xQVecFT0a = collision.q2xft0a();
      yQVecFT0a = collision.q2yft0a();
      xQVecFT0c = collision.q2xft0c();
      yQVecFT0c = collision.q2yft0c();
      xQVecFT0m = collision.q2xft0m();
      yQVecFT0m = collision.q2yft0m();
      xQVecBPos = collision.q2xbpos();
      yQVecBPos = collision.q2ybpos();
      xQVecBNeg = collision.q2xbneg();
      yQVecBNeg = collision.q2ybneg();
      xQVecBTot = collision.q2xbtot();
      yQVecBTot = collision.q2ybtot();
    } else if (harmonic == kTriangluar) {
      xQVecFT0a = collision.q3xft0a();
      yQVecFT0a = collision.q3yft0a();
      xQVecFT0c = collision.q3xft0c();
      yQVecFT0c = collision.q3yft0c();
      xQVecFT0m = collision.q3xft0m();
      yQVecFT0m = collision.q3yft0m();
      xQVecBPos = collision.q3xbpos();
      yQVecBPos = collision.q3ybpos();
      xQVecBNeg = collision.q3xbneg();
      yQVecBNeg = collision.q3ybneg();
      xQVecBTot = collision.q3xbtot();
      yQVecBTot = collision.q3ybtot();
    }

    if (saveSPResoHist) {
      registry.fill(HIST("spReso/hSpResoFT0cFT0a"), centrality, xQVecFT0c * xQVecFT0a + yQVecFT0c * yQVecFT0a);
      registry.fill(HIST("spReso/hSpResoFT0cTPCpos"), centrality, xQVecFT0c * xQVecBPos + yQVecFT0c * yQVecBPos);
      registry.fill(HIST("spReso/hSpResoFT0cTPCneg"), centrality, xQVecFT0c * xQVecBNeg + yQVecFT0c * yQVecBNeg);
      registry.fill(HIST("spReso/hSpResoFT0cTPCtot"), centrality, xQVecFT0c * xQVecBTot + yQVecFT0c * yQVecBTot);
      registry.fill(HIST("spReso/hSpResoFT0aTPCpos"), centrality, xQVecFT0a * xQVecBPos + yQVecFT0a * yQVecBPos);
      registry.fill(HIST("spReso/hSpResoFT0aTPCneg"), centrality, xQVecFT0a * xQVecBNeg + yQVecFT0a * yQVecBNeg);
      registry.fill(HIST("spReso/hSpResoFT0aTPCtot"), centrality, xQVecFT0a * xQVecBTot + yQVecFT0a * yQVecBTot);
      registry.fill(HIST("spReso/hSpResoFT0mTPCpos"), centrality, xQVecFT0m * xQVecBPos + yQVecFT0m * yQVecBPos);
      registry.fill(HIST("spReso/hSpResoFT0mTPCneg"), centrality, xQVecFT0m * xQVecBNeg + yQVecFT0m * yQVecBNeg);
      registry.fill(HIST("spReso/hSpResoFT0mTPCtot"), centrality, xQVecFT0m * xQVecBTot + yQVecFT0m * yQVecBTot);
      registry.fill(HIST("spReso/hSpResoTPCposTPCneg"), centrality, xQVecBPos * xQVecBNeg + yQVecBPos * yQVecBNeg);
    }

    if (saveEpResoHisto) {
      float epFT0a = epHelper.GetEventPlane(xQVecFT0a, yQVecFT0a, harmonic);
      float epFT0c = epHelper.GetEventPlane(xQVecFT0c, yQVecFT0c, harmonic);
      float epFT0m = epHelper.GetEventPlane(xQVecFT0m, yQVecFT0m, harmonic);
      float epBPoss = epHelper.GetEventPlane(xQVecBPos, yQVecBPos, harmonic);
      float epBNegs = epHelper.GetEventPlane(xQVecBNeg, yQVecBNeg, harmonic);
      float epBTots = epHelper.GetEventPlane(xQVecBTot, yQVecBTot, harmonic);

      registry.fill(HIST("hEventPlaneAngleFT0M"), centrality, epFT0m);
      registry.fill(HIST("hEventPlaneAngleTPCpos"), centrality, epBPoss);
      registry.fill(HIST("hEventPlaneAngleTPCneg"), centrality, epBNegs);

      registry.fill(HIST("QVector/hQVecMeanRVsPhiFT0a"), centrality, std::atan2(yQVecFT0a, xQVecFT0a), std::hypot(xQVecFT0a, yQVecFT0a));
      registry.fill(HIST("QVector/hQVecMeanRVsPhiFT0c"), centrality, std::atan2(yQVecFT0c, xQVecFT0c), std::hypot(xQVecFT0c, yQVecFT0c));
      registry.fill(HIST("QVector/hQVecMeanRVsPhiFT0m"), centrality, std::atan2(yQVecFT0m, xQVecFT0m), std::hypot(xQVecFT0m, yQVecFT0m));
      registry.fill(HIST("QVector/hQVecMeanRVsPhiTPCpos"), centrality, std::atan2(yQVecBPos, xQVecBPos), std::hypot(xQVecBPos, yQVecBPos));
      registry.fill(HIST("QVector/hQVecMeanRVsPhiTPCneg"), centrality, std::atan2(yQVecBNeg, xQVecBNeg), std::hypot(xQVecBNeg, yQVecBNeg));
      registry.fill(HIST("QVector/hQVecMeanRVsPhiTPCTot"), centrality, std::atan2(yQVecBTot, xQVecBTot), std::hypot(xQVecBTot, yQVecBTot));

      registry.fill(HIST("epReso/hEpResoFT0cFT0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epFT0a)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBTots)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBTots)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBTots)));
      registry.fill(HIST("epReso/hEpResoTPCposTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epBPoss, epBNegs)));
      for (int n = 1; n <= kOctagonal; n++) {
        registry.fill(HIST("epReso/hEpCosCoefficientsFT0c"), centrality, n, std::cos(n * harmonic * epFT0c));
        registry.fill(HIST("epReso/hEpSinCoefficientsFT0c"), centrality, n, std::sin(n * harmonic * epFT0c));
        registry.fill(HIST("epReso/hEpCosCoefficientsFT0a"), centrality, n, std::cos(n * harmonic * epFT0a));
        registry.fill(HIST("epReso/hEpSinCoefficientsFT0a"), centrality, n, std::sin(n * harmonic * epFT0a));
        registry.fill(HIST("epReso/hEpCosCoefficientsFT0m"), centrality, n, std::cos(n * harmonic * epFT0m));
        registry.fill(HIST("epReso/hEpSinCoefficientsFT0m"), centrality, n, std::sin(n * harmonic * epFT0m));
        registry.fill(HIST("epReso/hEpCosCoefficientsTPCpos"), centrality, n, std::cos(n * harmonic * epBPoss));
        registry.fill(HIST("epReso/hEpSinCoefficientsTPCpos"), centrality, n, std::sin(n * harmonic * epBPoss));
        registry.fill(HIST("epReso/hEpCosCoefficientsTPCneg"), centrality, n, std::cos(n * harmonic * epBNegs));
        registry.fill(HIST("epReso/hEpSinCoefficientsTPCneg"), centrality, n, std::sin(n * harmonic * epBNegs));
        registry.fill(HIST("epReso/hEpCosCoefficientsTPCTots"), centrality, n, std::cos(n * harmonic * epBTots));
        registry.fill(HIST("epReso/hEpSinCoefficientsTPCTots"), centrality, n, std::sin(n * harmonic * epBTots));
      }
    }
  }
  PROCESS_SWITCH(TaskPi0FlowEMC, processResolution, "Process resolution", false);

  // EMCal calibration
  void processEMCalCalib(Colls const& collisions, EMCalPhotons const& clusters)
  {
    float energyCorrectionFactor = 1.f;
    if (!correctionConfig.doEMCalCalib) {
      return;
    }
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

      if (emccuts.cfgEnableQA.value) {
        for (const auto& photon : photonsPerCollision) {
          registry.fill(HIST("clusterQA/hEClusterBefore"), photon.e()); // before cuts
          if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(photon))) {
            continue;
          }
          registry.fill(HIST("clusterQA/hEClusterAfter"), photon.e()); // accepted after cuts
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
        if (mesonConfig.cfgEnableQA.value) {
          registry.fill(HIST("mesonQA/hInvMassPt"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("mesonQA/hTanThetaPhi"), vMeson.M(), getAngleDegree(std::atan(dTheta / dPhi)));
          registry.fill(HIST("mesonQA/hAlphaPt"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
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
  PROCESS_SWITCH(TaskPi0FlowEMC, processEMCalCalib, "Process EMCal calibration mixed event", false);

  // Pi0 from EMCal
  void processM02(CollsWithQvecs const& collisions, EMCalPhotons const& clusters)
  {
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

      for (const auto& photon : photonsPerCollision) {
        if (emccuts.cfgEnableQA.value) {
          registry.fill(HIST("clusterQA/hEClusterBefore"), photon.e());                      // before cuts
          registry.fill(HIST("clusterQA/hClusterEtaPhiBefore"), photon.phi(), photon.eta()); // before cuts
        }
        if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(photon))) {
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

}; // End struct TaskPi0FlowEMC

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TaskPi0FlowEMC>(cfgc)};
}
