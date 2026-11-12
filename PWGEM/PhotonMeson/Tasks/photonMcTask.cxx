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

/// \file photonMcTask.cxx
/// \brief Analysis task for to obtain photon acceptance, efficiency and purity for photons from PCM and/or EMCal
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Core/EMBitFlags.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Concepts.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <cmath>
#include <string>
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
  CFT0C = 2,
  CFT0M = 3,
  NCentralityEstimators = 4
};

enum class MapLevel {
  kGood = 1,
  kNoBad = 2,
  kInEMC = 3,
  kAll = 4
};

// all values here have a little bit of headroom
constexpr float EmcalMaxEta = 0.75f; // it is actuall 0.7
constexpr float DcalMinEta = 0.2f;   // it is actuall 0.23

constexpr float EmcalLowerPhi = 1.4f - 0.01f;
constexpr float EmcalUpperPhi = 3.26f + 0.01f;
constexpr float DcalLowerPhi = 4.54f + 0.01f;
constexpr float DcalUpperPhosGap = 5.58f + 0.01f;
constexpr float DcalUpperPhi = 5.70f + 0.01f;

struct PhotonMcTask {
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.9f, "max. eta for generated particles"};

  // configurable axis
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {400, 0., 20.}, "pT axis for the neutral meson"};
  ConfigurableAxis thnConfigAxisERelative{"thnConfigAxisERelative", {600, -1., 5.}, "(E rec - E true) / E true axis"};
  ConfigurableAxis thnConfigAxisPRelative{"thnConfigAxisPRelative", {300, -1., 2.}, "(P rec - P true) / P true axis"};
  ConfigurableAxis thnConfigAxisEtaRelative{"thnConfigAxisEtaRelative", {300, -1., 2.}, "(eta rec - eta true) / eta true axis"};
  ConfigurableAxis thnConfigAxisPhiRelative{"thnConfigAxisPhiRelative", {300, -1., 2.}, "(phi rec - phi true) / phi true axis"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {20, 0., 100.}, "centrality axis for the current event"};
  ConfigurableAxis thnConfigAxisMult{"thnConfigAxisMult", {60, 0., 60000.}, "multiplicity axis for the current event"};
  Configurable<bool> useCent{"useCent", false, "flag to enable usage of centrality instead of multiplicity as axis."};

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
    Configurable<std::vector<float>> cfgEMCTMEta{"cfgEMCTMEta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal primary track matching"};
    Configurable<std::vector<float>> cfgEMCTMPhi{"cfgEMCTMPhi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal primary track matching"};
    Configurable<std::vector<float>> emcSecTMEta{"emcSecTMEta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal secondary track matching"};
    Configurable<std::vector<float>> emcSecTMPhi{"emcSecTMPhi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal secondary track matching"};
    Configurable<float> cfgEMCEoverp{"cfgEMCEoverp", 1.75, "Minimum cluster energy over track momentum for EMCal primary track matching"};
    Configurable<bool> cfgEMCUseExoticCut{"cfgEMCUseExoticCut", true, "Flag to use the EMCal exotic cluster cut"};
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

  SliceCache cache;

  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::MinClusters, aod::EMEMCClusterMCLabels_001>;
  using PcmPhotons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;

  using PcmMcLegs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;

  using Colls = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMMCEventLabels, aod::EmMagFields>;

  using McColls = o2::soa::Join<o2::aod::EMMCEvents, o2::aod::BinnedGenPts>;
  using McParticles = EMMCParticles;

  Partition<EMMCParticles> mcPhotons = aod::mcparticle::pdgCode == static_cast<int>(PDG_t::kGamma);
  PresliceUnsorted<EMMCParticles> perMcEvent = aod::emmcparticle::emmceventId;

  PresliceOptional<EMCalPhotons> perCollisionEMC = o2::aod::emccluster::pmeventId;
  PresliceOptional<PcmPhotons> perCollisionPCM = aod::v0photonkf::pmeventId;
  PresliceOptional<MinMTracks> perEMCClusterMT = o2::aod::mintm::minClusterId;
  PresliceOptional<MinMSTracks> perEMCClusterMS = o2::aod::mintm::minClusterId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb{};
  int mRunNumber{-1};
  float dBz{0.f};

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
    fEMEventCut.SetFT0COccupancyRange(eventcuts.cfgFT0COccupancyMin, eventcuts.cfgFT0COccupancyMax);
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

    defineEMEventCut();
    defineEMCCut();
    fEMCCut.addQAHistograms(&registry);
    definePCMCut();
    fV0PhotonCut.addQAHistograms(&registry);
    fEMEventCut.addEventHistograms(&registry);

    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPtGen{thnConfigAxisPt, "#it{p}_{T,Gen} (GeV/#it{c})"};
    const AxisSpec thnAxisPtRec{thnConfigAxisPt, "#it{p}_{T,Rec} (GeV/#it{c})"};
    const AxisSpec thnAxisPGen{thnConfigAxisPt, "#it{p}_{Gen} (GeV/#it{c})"};
    const AxisSpec thnAxisPRec{thnConfigAxisPt, "#it{p}_{Rec} (GeV/#it{c})"};
    const AxisSpec thnAxisPRelative{thnConfigAxisPRelative, "#it{p}_{Rec} - #it{p}_{Gen} / #it{p}_{Gen}"};
    const AxisSpec thnAxisEGen{thnConfigAxisPt, "#it{E}_{Gen} (GeV)"};
    const AxisSpec thnAxisERec{thnConfigAxisPt, "#it{E}_{Rec} (GeV)"};
    const AxisSpec thnAxisERelative{thnConfigAxisERelative, "#it{E}_{Rec} - #it{E}_{Gen} / #it{E}_{Gen}"};

    const AxisSpec thnAxisEtaRelative{thnConfigAxisEtaRelative, "#it{#eta}_{Rec} - #it{#eta}_{Gen} / #it{#eta}_{Gen}"};
    const AxisSpec thnAxisPhiRelative{thnConfigAxisPhiRelative, "#it{#varphi}_{Rec} - #it{#varphi}_{Gen} / #it{#varphi}_{Gen}"};

    const AxisSpec thnAxisEtaGen{280, -0.7, 0.7, "#it{#eta}_{Gen}"};
    const AxisSpec thnAxisPhiGen{360, 0., o2::constants::math::TwoPI, "#it{#varphi}_{Gen} (rad)"};

    AxisSpec thnAxisCentOrMult{1, 0., 1., "Centrality/Multiplicity"}; // placeholder, overwritten in init
    if (useCent.value) {
      // PbPb: use centrality
      thnAxisCentOrMult = {thnConfigAxisCent, "Centrality (%)"};
    } else {
      // pp: use multiplicity
      thnAxisCentOrMult = {thnConfigAxisMult, "FT0C Multiplicity"};
    }

    registry.add("hGenerated", "generated photons pT vs cent/mult", HistType::kTH2D, {thnAxisPt, thnAxisCentOrMult});

    if (doprocessEmcal || doprocessPcmEmcal) {
      registry.add("EMCal/hAccepted", "EMCal accepted photon pt vs cent/mult", HistType::kTH2D, {thnAxisPt, thnAxisCentOrMult});
      registry.add("EMCal/hMeasured", "EMCal measured photon candidate pt vs cent/mult ", HistType::kTH2D, {thnAxisPt, thnAxisCentOrMult});
      registry.add("EMCal/hMeasuredTruePhoton", "EMCal measured and validated photon pt vs cent/mult ", HistType::kTH2D, {thnAxisPt, thnAxisCentOrMult});
    }

    if (doprocessPcm || doprocessPcmEmcal) {
      registry.add("PCM/hAccepted", "PCM accepted photon pt vs cent/mult", HistType::kTH2D, {thnAxisPt, thnAxisCentOrMult});
      registry.add("PCM/hMeasured", "PCM measured photon candidate pt vs cent/mult ", HistType::kTH2D, {thnAxisPt, thnAxisCentOrMult});
      registry.add("PCM/hMeasuredTruePhoton", "PCM measured and validated photon pt vs cent/mult ", HistType::kTH2D, {thnAxisPt, thnAxisCentOrMult});
    }

  }; // end init

  template <o2::soa::is_iterator TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }
    mRunNumber = collision.runNumber();

    dBz = collision.grpMagField().getNominalL3Field();
    LOG(info) << "Retrieved GRP for timestamp " << collision.timestamp() << " with magnetic field of " << dBz << " kZG";
    fV0PhotonCut.SetD_Bz(dBz);
  }

  template <o2::soa::is_iterator TCollision>
  float getCentralityOrMultiplicity(TCollision const& collision)
  {
    if (useCent.value) {
      return getCentrality(collision);
    }
    // pp: use raw FT0C multiplicity
    return collision.multFT0C();
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
      fEMEventCut.fillEventInfo<0>(&registry, collision);
    }
    if (!(fEMEventCut.IsSelected(collision))) {
      // general event selection
      return false;
    }
    float centOrMult = getCentralityOrMultiplicity(collision);
    if (useCent && (centOrMult < eventcuts.cfgMinCent || centOrMult > eventcuts.cfgMaxCent)) {
      // event selection
      return false;
    }

    if (collision.emmceventId() < 0) {
      return false;
    }

    if (fillHisto) {
      fEMEventCut.fillEventInfo<1>(&registry, collision);
      fEMEventCut.fillAccepted(&registry);
    }

    return true;
  }

  /// \brief check if particle is within the EMCal acceptance
  /// \param particle particle will be checked
  /// \return true if particle is within the EMCal, false otherwise
  template <o2::soa::is_iterator TParticle>
  bool inEMCalAcceptance(TParticle const& particle)
  {
    // general EMCal eta acceptance
    if (std::abs(particle.eta()) > EmcalMaxEta) {
      return false;
    }

    // before EMCal or after DCal
    if (particle.phi() < EmcalLowerPhi || particle.phi() > DcalUpperPhi) {
      return false;
    }

    // between EMCal and DCal
    if (EmcalUpperPhi < particle.phi() && particle.phi() < DcalLowerPhi) {
      return false;
    }

    // PHOS GAP in DCal
    if (particle.phi() >= DcalLowerPhi && particle.phi() <= DcalUpperPhosGap) {
      if (std::abs(particle.eta()) < DcalMinEta) {
        return false;
      }
    }

    return true;
  }

  // TODO:
  // - efficiency: number of found photon candidates / number of photons within accpetance
  // - acceptance: number of photons in acceptance / number of photons generated
  // - photon Purity: number of actual photons / number of photon candidates
  // - How to handle EMCal cluster in this regard? Sometimes both conversion leptons can end up in the same cluster.
  // What is a good check when to start treating them as a photon?
  // When they both deposite at least X% of their momentum into the cluster?
  // - electron "purity": number of actual electrons/ number of photon candidates

  // PCM
  void processPcm(Colls const& collisions, PcmPhotons const& photons, PcmMcLegs const& legs, EMMCParticles const& mcParticles)
  {
    EMBitFlags v0flags(photons.size());
    if (photons.size() > 0) {
      fV0PhotonCut.AreSelectedRunning<decltype(photons), PcmMcLegs>(v0flags, photons, &registry);
    }

    // leg iterators for PCM
    auto pos1 = legs.begin();
    auto ele1 = legs.begin();

    // MC leg iterators for PCM
    auto pos1mc = mcParticles.begin();
    auto ele1mc = mcParticles.begin();

    for (const auto& collision : collisions) {

      initCCDB(collision);
      if (!isFullEventSelected(collision, true)) {
        continue;
      }

      auto centOrMult = getCentralityOrMultiplicity(collision);
      auto photonsThisEvent = mcPhotons->sliceByCachedUnsorted(aod::emmcparticle::emmceventId, collision.emmceventId(), cache);

      for (const auto& truePhoton : photonsThisEvent) {
        if (std::abs(truePhoton.eta()) <= cfgMaxEta) {
          registry.fill(HIST("hGenerated"), truePhoton.pt(), centOrMult);
          // trivial case I guess for PCM
          registry.fill(HIST("PCM/hAccepted"), truePhoton.pt(), centOrMult);
        }
      }

      auto photonsPCMPerCollision = photons.sliceBy(perCollisionPCM, collision.globalIndex());

      for (const auto& photon : photonsPCMPerCollision) {
        if (!(v0flags.test(photon.globalIndex()))) {
          continue;
        }
        registry.fill(HIST("PCM/hMeasured"), photon.pt(), centOrMult);

        pos1.setCursor(photon.posTrackId());
        ele1.setCursor(photon.negTrackId());

        pos1mc.setCursor(pos1.emmcparticleId());
        ele1mc.setCursor(ele1.emmcparticleId());

        const auto photonid1 = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(pos1mc, ele1mc, PDG_t::kPositron, PDG_t::kElectron, PDG_t::kGamma, mcParticles);

        if (photonid1 < 0) {
          continue;
        }
        registry.fill(HIST("PCM/hMeasuredTruePhoton"), photon.pt(), centOrMult);
      }
    }
  }
  PROCESS_SWITCH(PhotonMcTask, processPcm, "Process for pcm photons", false)

  // EMCal
  void processEmcal(Colls const& collisions, EMCalPhotons const& clusters, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds, EMMCParticles const& mcParticles)
  {
    EMBitFlags emcFlags(clusters.size());
    if (clusters.size() > 0) {
      fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds, &registry);
    }

    auto mcPhoton = mcParticles.begin();

    for (const auto& collision : collisions) {

      initCCDB(collision);
      if (!isFullEventSelected(collision, true)) {
        continue;
      }

      auto centOrMult = getCentralityOrMultiplicity(collision);
      auto photonsThisEvent = mcPhotons->sliceByCachedUnsorted(aod::emmcparticle::emmceventId, collision.emmceventId(), cache);

      for (const auto& truePhoton : photonsThisEvent) {
        if (std::abs(truePhoton.eta()) <= cfgMaxEta) {
          registry.fill(HIST("hGenerated"), truePhoton.pt(), centOrMult);

          if (inEMCalAcceptance(truePhoton)) {
            registry.fill(HIST("EMCal/hAccepted"), truePhoton.pt(), centOrMult);
          }
        }
      }

      auto photonsEMCPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());

      for (const auto& photonEMC : photonsEMCPerCollision) {
        if (!(emcFlags.test(photonEMC.globalIndex()))) {
          continue;
        }

        registry.fill(HIST("EMCal/hMeasured"), photonEMC.pt(), centOrMult);

        if (photonEMC.emmcparticleIds().empty()) {
          // this is a cluster with just noise, skip
          continue;
        }
        // we only want to look at the largest contribution
        mcPhoton.setCursor(photonEMC.emmcparticleIds()[0]);

        // if the largest contribution is photon, just fill the photon histograms
        if (std::abs(mcPhoton.pdgCode()) != PDG_t::kGamma) {
          continue;
        }

        registry.fill(HIST("EMCal/hMeasuredTruePhoton"), photonEMC.pt(), centOrMult);
      }
    }
  }
  PROCESS_SWITCH(PhotonMcTask, processEmcal, "Process for emcal photons", false)

  // PCM-EMCal same event
  void processPcmEmcal(Colls const& collisions, EMCalPhotons const& clusters, PcmPhotons const& photons, PcmMcLegs const& legs, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds, EMMCParticles const& mcParticles)
  {
    EMBitFlags emcFlags(clusters.size());
    if (clusters.size() > 0) {
      fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds, &registry);
    }
    EMBitFlags v0flags(photons.size());
    if (photons.size() > 0) {
      fV0PhotonCut.AreSelectedRunning<decltype(photons), PcmMcLegs>(v0flags, photons, &registry);
    }

    // create iterators for photon mc particles
    auto mcPhoton = mcParticles.begin();

    // leg iterators for PCM
    auto pos1 = legs.begin();
    auto ele1 = legs.begin();

    // MC leg iterators for PCM
    auto pos1mc = mcParticles.begin();
    auto ele1mc = mcParticles.begin();

    for (const auto& collision : collisions) {
      initCCDB(collision);
      if (!isFullEventSelected(collision, true)) {
        continue;
      }

      float centOrMult = getCentralityOrMultiplicity(collision);

      auto photonsThisEvent = mcPhotons->sliceByCachedUnsorted(aod::emmcparticle::emmceventId, collision.emmceventId(), cache);

      for (const auto& truePhoton : photonsThisEvent) {
        if (std::abs(truePhoton.eta()) <= cfgMaxEta) {
          registry.fill(HIST("hGenerated"), truePhoton.pt(), centOrMult);
          registry.fill(HIST("PCM/hAccepted"), truePhoton.pt(), centOrMult);
        }
        if (inEMCalAcceptance(truePhoton)) {
          registry.fill(HIST("EMCal/hAccepted"), truePhoton.pt(), centOrMult);
        }
      }

      auto photonsEMCPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());
      auto photonsPCMPerCollision = photons.sliceBy(perCollisionPCM, collision.globalIndex());

      for (const auto& photonEMC : photonsEMCPerCollision) {
        if (!(emcFlags.test(photonEMC.globalIndex()))) {
          continue;
        }

        registry.fill(HIST("EMCal/hMeasured"), photonEMC.pt(), centOrMult);

        if (photonEMC.emmcparticleIds().empty()) {
          // this is a cluster with just noise, skip
          continue;
        }
        // we only want to look at the largest contribution
        mcPhoton.setCursor(photonEMC.emmcparticleIds()[0]);

        // if the largest contribution is photon, just fill the photon histograms
        if (std::abs(mcPhoton.pdgCode()) != PDG_t::kGamma) {
          continue;
        }

        registry.fill(HIST("EMCal/hMeasuredTruePhoton"), photonEMC.pt(), centOrMult);
      }

      for (const auto& photon : photonsPCMPerCollision) {
        if (!(v0flags.test(photon.globalIndex()))) {
          continue;
        }
        registry.fill(HIST("PCM/hMeasured"), photon.pt(), centOrMult);

        pos1.setCursor(photon.posTrackId());
        ele1.setCursor(photon.negTrackId());

        pos1mc.setCursor(pos1.emmcparticleId());
        ele1mc.setCursor(ele1.emmcparticleId());

        const auto photonid1 = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(pos1mc, ele1mc, PDG_t::kPositron, PDG_t::kElectron, PDG_t::kGamma, mcParticles);

        if (photonid1 < 0) {
          continue;
        }
        registry.fill(HIST("PCM/hMeasuredTruePhoton"), photon.pt(), centOrMult);
      }
    }
  }
  PROCESS_SWITCH(PhotonMcTask, processPcmEmcal, "Process for pcm and emcal photons", true);

}; // End struct PhotonMcTask

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<PhotonMcTask>(context)};
}
