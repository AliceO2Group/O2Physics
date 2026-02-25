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

/// \file photonResoTask.cxx
/// \brief Analysis task for to obtain photon resolution histograms in MC
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Core/EMBitFlags.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TF1.h>
#include <TH1.h>
#include <TPDGCode.h>

#include <cmath>
#include <memory>
#include <string>
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

enum class MapLevel {
  kGood = 1,
  kNoBad = 2,
  kInEMC = 3,
  kAll = 4
};

struct PhotonResoTask {
  static constexpr float MinEnergy = 0.7f;

  o2::framework::Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  o2::framework::Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  o2::framework::Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  o2::framework::Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};

  // configurable axis
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {400, 0.0, 0.8}, "invariant mass axis for the neutral meson"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0., 20.}, "pT axis for the neutral meson"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {20, 0., 100.}, "centrality axis for the current event"};

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
    Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3)"};
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
    Configurable<bool> cfgEnableQA{"cfgEnableQA", false, "flag to turn QA plots on/off"};
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

  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::MinClusters, aod::EMEMCClusterMCLabels>;
  using PcmPhotons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;

  using PcmMcLegs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;

  using Colls = soa::Join<aod::EMEvents_004, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMMCEventLabels>;

  using McColls = o2::soa::Join<o2::aod::EMMCEvents, o2::aod::BinnedGenPts>;
  using McParticles = EMMCParticles;

  PresliceOptional<EMCalPhotons> perCollisionEMC = o2::aod::emccluster::emphotoneventId;
  PresliceOptional<PcmPhotons> perCollisionPCM = aod::v0photonkf::emphotoneventId;
  PresliceOptional<MinMTracks> perEMCClusterMT = o2::aod::mintm::minClusterId;
  PresliceOptional<MinMSTracks> perEMCClusterMS = o2::aod::mintm::minClusterId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber{-1};
  float dBz{0.f};

  // Usage when cfgEnableNonLin is enabled
  std::unique_ptr<TF1> fEMCalCorrectionFactor; // ("fEMCalCorrectionFactor","(1 + [0]/x + [1]/x^2) / (1 + [2]/x)", 0.3, 100.);
  float energyCorrectionFactor = 1.f;

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
    mRunNumber = 0;
    dBz = 0;

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    defineEMEventCut();
    defineEMCCut();
    fEMCCut.addQAHistograms(&registry);
    definePCMCut();
    fV0PhotonCut.addQAHistograms(&registry);
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&registry);

    const AxisSpec thnAxisPtGen{thnConfigAxisPt, "#it{p}_{T,Gen} (GeV/#it{c})"};
    const AxisSpec thnAxisPtRec{thnConfigAxisPt, "#it{p}_{T,Rec} (GeV/#it{c})"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality (%)"};
    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};

    registry.add("EMCal/hPhotonReso", "EMCal photon rec pT vs true pT vs cent", HistType::kTH3D, {thnAxisPtRec, thnAxisPtGen, thnAxisCent});
    registry.add("EMCal/hConvPhotonReso", "EMCal conversion photon rec pT vs true pT vs cent ", HistType::kTH3D, {thnAxisPtRec, thnAxisPtGen, thnAxisCent});

    registry.add("EMCal/hPi0Reso", "EMCal pi0 rec pT vs true pT vs min vs cent ", HistType::kTHnSparseF, {thnAxisPtRec, thnAxisPtGen, thnConfigAxisInvMass, thnAxisCent});
    registry.add("EMCal/hEtaReso", "EMCal eta rec pT vs true pT vs min vs cent ", HistType::kTHnSparseF, {thnAxisPtRec, thnAxisPtGen, thnConfigAxisInvMass, thnAxisCent});

    registry.add("PCM/hPhotonReso", "PCM  photon rec pT vs true pT vs ", HistType::kTH3D, {thnAxisPtRec, thnAxisPtGen, thnAxisCent});

    auto hMesonCuts = registry.add<TH1>("hMesonCuts", "hMesonCuts;;Counts", kTH1D, {{6, 0.5, 6.5}}, false);
    hMesonCuts->GetXaxis()->SetBinLabel(1, "in");
    hMesonCuts->GetXaxis()->SetBinLabel(2, "opening angle");
    hMesonCuts->GetXaxis()->SetBinLabel(3, "#it{M}_{#gamma#gamma}");
    hMesonCuts->GetXaxis()->SetBinLabel(4, "#it{p}_{T}");
    hMesonCuts->GetXaxis()->SetBinLabel(5, "conversion cut");
    hMesonCuts->GetXaxis()->SetBinLabel(6, "out");
    if (mesonConfig.cfgEnableQA.value) {
      registry.add("mesonQA/hInvMassPt", "Histo for inv pair mass vs pt", HistType::kTH2D, {thnAxisInvMass, thnAxisPtRec});
    }

    fEMCalCorrectionFactor = std::make_unique<TF1>("fEMCalCorrectionFactor", "(1 + [0]/x + [1]/x^2) / (1 + [2]/x)", 0.3, 100.);
    fEMCalCorrectionFactor->SetParameters(-5.33426e-01, 1.40144e-02, -5.24434e-01);
  }; // end init

  template <o2::soa::is_iterator TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }

    auto run3GrpTimestamp = collision.timestamp();
    o2::parameters::GRPObject* grpo = nullptr;
    o2::parameters::GRPMagField* grpmag = nullptr;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3GrpTimestamp);
    if (grpo) {
      // Fetch magnetic field from ccdb for current collision
      dBz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3GrpTimestamp << " with magnetic field of " << dBz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3GrpTimestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3GrpTimestamp;
      }
      // Fetch magnetic field from ccdb for current collision
      dBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3GrpTimestamp << " with magnetic field of " << dBz << " kZG";
    }
    fV0PhotonCut.SetD_Bz(dBz);
    mRunNumber = collision.runNumber();
  }

  /// Get the centrality
  /// \param collision is the collision with the centrality information
  template <o2::soa::is_iterator TCollision>
  float getCentrality(TCollision const& collision)
  {
    float cent = -999.;
    switch (eventcuts.centEstimator) {
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
    if (fillHisto) {
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&registry, collision);
      registry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // accepted
      registry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // accepted
    }
    return true;
  }

  // PCM-EMCal same event
  void processPcmEmcal(Colls const& collisions, EMCalPhotons const& clusters, PcmPhotons const& photons, PcmMcLegs const& legs, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds, EMMCParticles const& mcParticles)
  {
    if (clusters.size() <= 0 && photons.size() < 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    EMBitFlags emcFlags(clusters.size());
    if (clusters.size() > 0) {
      fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds, &registry);
    }
    EMBitFlags v0flags(photons.size());
    if (photons.size() > 0) {
      fV0PhotonCut.AreSelectedRunning<decltype(photons), PcmMcLegs>(v0flags, photons, &registry);
    }

    // create iterators for photon mc particles
    auto mcPhoton1 = mcParticles.begin();
    auto mcPhoton2 = mcParticles.begin();

    // leg iterators for PCM
    auto pos1 = legs.begin();
    auto ele1 = legs.begin();

    // MC leg iterators for PCM
    auto pos1mc = mcParticles.begin();
    auto ele1mc = mcParticles.begin();

    for (const auto& collision : collisions) {
      initCCDB(collision);
      isFullEventSelected(collision, true);

      float cent = getCentrality(collision);

      auto photonsEMCPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());
      auto photonsPCMPerCollision = photons.sliceBy(perCollisionPCM, collision.globalIndex());

      for (const auto& photonEMC : photonsEMCPerCollision) {
        if (!(emcFlags.test(photonEMC.globalIndex()))) {
          continue;
        }
        if (photonEMC.emmcparticleIds().size() <= 0) {
          // this is a cluster with just noise, skip
          continue;
        }
        // we only want to look at the largest contribution
        mcPhoton1.setCursor(photonEMC.emmcparticleIds()[0]);

        if (std::abs(mcPhoton1.pdgCode()) == PDG_t::kGamma) {
          registry.fill(HIST("EMCal/hPhotonReso"), photonEMC.pt(), mcPhoton1.pt(), cent);
        } else if (std::abs(mcPhoton1.pdgCode()) == PDG_t::kElectron) {
          if (!o2::aod::pwgem::photonmeson::utils::mcutil::isMotherPDG(mcPhoton1, PDG_t::kGamma)) {
            continue;
          }
          registry.fill(HIST("EMCal/hConvPhotonReso"), photonEMC.pt(), mcPhoton1.pt(), cent);
        }
      }

      for (const auto& photonPCM : photonsPCMPerCollision) {
        if (!(v0flags.test(photonPCM.globalIndex()))) {
          continue;
        }

        pos1.setCursor(photonPCM.posTrackId());
        ele1.setCursor(photonPCM.negTrackId());

        pos1mc.setCursor(pos1.emmcparticleId());
        ele1mc.setCursor(ele1.emmcparticleId());

        const auto photonid1 = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcParticles);

        if (photonid1 < 0) {
          continue;
        }

        float trueConvRadius = std::hypot(ele1mc.vx(), ele1mc.vy());

        mcPhoton1.setCursor(photonid1);

        if (!fV0PhotonCut.IsConversionPointInAcceptance(mcPhoton1, trueConvRadius)) {
          continue;
        }

        registry.fill(HIST("PCM/hPhotonReso"), photonPCM.pt(), mcPhoton1.pt(), cent);

      } // end of loop over pcm photons

      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photonsEMCPerCollision, photonsEMCPerCollision))) {
        if (!(emcFlags.test(g1.globalIndex())) || !(emcFlags.test(g2.globalIndex()))) {
          continue;
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
        }

        if (g1.emmcparticleIds().size() <= 0 || g2.emmcparticleIds().size() <= 0) {
          // there is a cluster which is just noise, skip
          continue;
        }

        mcPhoton1.setCursor(g1.emmcparticleIds()[0]);
        mcPhoton2.setCursor(g2.emmcparticleIds()[0]);

        int photonid1 = -1, photonid2 = -1, pi0id = -1, etaid = -1;
        photonid1 = o2::aod::pwgem::photonmeson::utils::mcutil::FindMotherInChain(mcPhoton1, mcParticles, std::vector<int>{111, 221});
        photonid2 = o2::aod::pwgem::photonmeson::utils::mcutil::FindMotherInChain(mcPhoton2, mcParticles, std::vector<int>{111, 221});

        if (photonid1 < 0 || photonid2 < 0) {
          continue;
        }
        mcPhoton1.setCursor(photonid1);
        mcPhoton2.setCursor(photonid2);

        pi0id = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(mcPhoton1, mcPhoton2, 22, 22, 111, mcParticles);
        etaid = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(mcPhoton1, mcPhoton2, 22, 22, 221, mcParticles);

        if (pi0id >= 0) {
          const auto pi0mc = mcParticles.iteratorAt(pi0id);
          registry.fill(HIST("EMCal/hPi0Reso"), vMeson.Pt(), pi0mc.pt(), vMeson.M(), cent);
        }
        if (etaid >= 0) {
          const auto etamc = mcParticles.iteratorAt(etaid);
          registry.fill(HIST("EMCal/hEtaReso"), vMeson.Pt(), etamc.pt(), vMeson.M(), cent);
        }
      }
    }
  }
  PROCESS_SWITCH(PhotonResoTask, processPcmEmcal, "Process for pcm and emcal photons", true);

}; // End struct PhotonResoTask

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PhotonResoTask>(cfgc)};
}
