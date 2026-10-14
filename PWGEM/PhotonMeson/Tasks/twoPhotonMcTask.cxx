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

/// \file twoPhotonMcTask.cxx
/// \brief Analysis task to obtain the invariant mass distributions for different MC truth classes
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Core/EMBitFlags.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

#include "Common/DataModel/Centrality.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Concepts.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;
using namespace o2::aod::pwgem::photonmeson::utils::mcutil;

static constexpr int kMaxChainLength = 16;
using Chain = std::array<int, kMaxChainLength>; // filled with -1 after the last ancestor

/// \brief MC truth of a single photon candidate, ordered from best to worst such that std::max gives the worse flag of a pair
enum class PhotonTruth : uint8_t {
  kTruePhoton = 0, // PCM: both legs from the same photon, EMC: leading contributor is a photon (or e+-, converted in front of the EMCal)
  kNonPhoton,      // PCM: both legs have the same direct mother which is not a photon (K0S, Lambda, pi0 Dalitz), EMC: leading contributor is not a photon
  kFake,           // PCM: leg without MC label or legs without common direct mother, EMC: cluster without MC label
  kNTruth
};

/// \brief Class of the closest common ancestor of the two photon candidates of a pair
enum class PairBucket : uint8_t {
  kNone = 0,         // no common ancestor
  kSameParticle,     // common ancestor is one of the two candidates themselves (EMC conversion split, Dalitz)
  kSignalPi0,        // pi0 is the direct mother of both
  kSignalEta,        // eta is the direct mother of both
  kCrossPi0,         // pi0 further up in the chain: the pi0 has no hadronic decays, so these are pi0 decays with an EM interaction in between (conversion, bremsstrahlung, shower in front of the EMCal)
  kCrossEta,         // eta further up in the chain: hadronic decays (eta -> 3 pi0, photons of different pi0) and eta decays with an EM interaction in between
  kCrossOmega,       // omega
  kCrossK0S,         // K0S
  kCrossEtaPrime,    // eta'
  kCrossOtherHadron, // any other hadron (including the rare Sigma0)
  kPartonString,     // quark, gluon, diquark or string/cluster: jet-like
  kNBuckets
};

/// \brief What the leading contributor of the EMC cluster is in a PCM-EMC same particle pair
enum class SameParticleOrigin : uint8_t {
  kLegElectron = 0, // the e+- that was reconstructed as one of the V0 legs, should have been removed by track matching
  kBremsFromLeg,    // a photon radiated by one of the V0 leg e+-, cannot be removed by track matching
  kOtherElectron,   // any other e+-
  kOther,           // anything else
  kNOrigins
};

/// \brief Why the V0 leg did not veto the EMC cluster, from the tracks matched to the cluster
enum class LegMatchStatus : uint8_t {
  kNoMatchStored = 0,  // no matched track stored for this cluster at all
  kNoCompatibleMatch,  // matched tracks stored, but none with the momentum of the leg
  kOutsideWindow,      // leg found, but outside the track matching window of the analysis
  kInsideWindowEoverP, // leg found inside the window, but E/p above the cut (primaries only)
  kInsideWindow,       // leg found inside the window: the cluster should have been vetoed
  kNStatus
};

enum CentralityEstimator {
  None = 0,
  CFT0A = 1,
  CFT0C = 2,
  CFT0M = 3,
  NCentralityEstimators = 4
};

struct TwoPhotonMcTask {
  // configurable axis
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {400, 0.0, 0.8}, "invariant mass axis for the neutral meson"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0., 20.}, "pT axis for the neutral meson"};
  ConfigurableAxis thnConfigAxisERelative{"thnConfigAxisERelative", {600, -1., 5.}, "(E rec - E true) / E true axis"};
  ConfigurableAxis thnConfigAxisPRelative{"thnConfigAxisPRelative", {300, -1., 2.}, "(P rec - P true) / P true axis"};
  ConfigurableAxis thnConfigAxisEtaRelative{"thnConfigAxisEtaRelative", {300, -1., 2.}, "(eta rec - eta true) / eta true axis"};
  ConfigurableAxis thnConfigAxisPhiRelative{"thnConfigAxisPhiRelative", {300, -1., 2.}, "(phi rec - phi true) / phi true axis"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {10, 0., 100.}, "centrality axis for the current event"};
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

  // same mixing as taskPi0FlowEMC, but without event plane bins: the derived MC has no Q vectors
  struct : ConfigurableGroup {
    std::string prefix = "mixingConfig";
    ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis cfgCentBins{"cfgCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "Mixing bins - centrality"};
    Configurable<int> cfgMixingDepth{"cfgMixingDepth", 2, "Mixing depth"};
    Configurable<bool> cfgCorrectMixedVtxEta{"cfgCorrectMixedVtxEta", true, "EMC-EMC: re-project the eta of cluster 2 onto the vertex of collision 1 in the mixed event"};
    Configurable<float> cfgEmcalEffRadius{"cfgEmcalEffRadius", 430.f, "effective EMCal radius (cm) used for the mixed event vertex correction"};
  } mixingConfig;

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  SliceCache cache;

  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::MinClusters, aod::EMEMCClusterMCLabels_001>;
  using PcmPhotons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;

  using PcmMcLegs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;

  using Colls = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMMCEventLabels, aod::EmMagFields>;

  using McColls = o2::soa::Join<o2::aod::EMMCEvents, o2::aod::BinnedGenPts>;
  using McParticles = EMMCParticles;

  PresliceOptional<EMCalPhotons> perCollisionEMC = o2::aod::emccluster::pmeventId;
  PresliceOptional<PcmPhotons> perCollisionPCM = aod::v0photonkf::pmeventId;
  PresliceOptional<MinMTracks> perEMCClusterMT = o2::aod::mintm::minClusterId;
  PresliceOptional<MinMSTracks> perEMCClusterMS = o2::aod::mintm::minClusterId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int mRunNumber{-1};
  float dBz{0.f};

  static constexpr int PhotonPCM = 0;
  static constexpr int PhotonEMC = 1;
  static constexpr int PairPCMPCM = 0;
  static constexpr int PairEMCEMC = 1;
  static constexpr int PairPCMEMC = 2;

  // ancestor chains and truth flags, filled once per DF, index = globalIndex of the V0 photon / cluster
  std::vector<Chain> mChainsPCM;
  std::vector<Chain> mChainsEMC;
  std::vector<PhotonTruth> mTruthPCM;
  std::vector<PhotonTruth> mTruthEMC;

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

    defineEMEventCut();
    defineEMCCut();
    fEMCCut.addQAHistograms(&registry);
    definePCMCut();
    fV0PhotonCut.addQAHistograms(&registry);
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&registry);

    const AxisSpec thnAxisPtGen{thnConfigAxisPt, "#it{p}_{T,Gen} (GeV/#it{c})"};
    const AxisSpec thnAxisPtRec{thnConfigAxisPt, "#it{p}_{T,Rec} (GeV/#it{c})"};
    const AxisSpec thnAxisPGen{thnConfigAxisPt, "#it{p}_{Gen} (GeV/#it{c})"};
    const AxisSpec thnAxisPRec{thnConfigAxisPt, "#it{p}_{Rec} (GeV/#it{c})"};
    const AxisSpec thnAxisPRelative{thnConfigAxisPRelative, "#it{p}_{Rec} - #it{p}_{Gen} / #it{p}_{Gen}"};
    const AxisSpec thnAxisEGen{thnConfigAxisPt, "#it{E}_{Gen} (GeV)"};
    const AxisSpec thnAxisERec{thnConfigAxisPt, "#it{E}_{Rec} (GeV)"};
    const AxisSpec thnAxisERelative{thnConfigAxisERelative, "#it{E}_{Rec} - #it{E}_{Gen} / #it{E}_{Gen}"};
    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};

    const AxisSpec thnAxisEtaRelative{thnConfigAxisEtaRelative, "#it{#eta}_{Rec} - #it{#eta}_{Gen} / #it{#eta}_{Gen}"};
    const AxisSpec thnAxisPhiRelative{thnConfigAxisPhiRelative, "#it{#varphi}_{Rec} - #it{#varphi}_{Gen} / #it{#varphi}_{Gen}"};

    const AxisSpec thnAxisEtaGen{280, -0.7, 0.7, "#it{#eta}_{Gen}"};
    const AxisSpec thnAxisPhiGen{360, 0., o2::constants::math::TwoPI, "#it{#varphi}_{Gen} (rad)"};

    const AxisSpec thnAxisPairTruth{static_cast<int>(PhotonTruth::kNTruth), -0.5, static_cast<int>(PhotonTruth::kNTruth) - 0.5, "pair truth (worse photon)"};
    const AxisSpec thnAxisPairBucket{static_cast<int>(PairBucket::kNBuckets), -0.5, static_cast<int>(PairBucket::kNBuckets) - 0.5, "closest common ancestor"};

    AxisSpec thnAxisCentOrMult{1, 0., 1., "Centrality/Multiplicity"}; // placeholder, overwritten in init
    if (useCent.value) {
      // PbPb: use centrality
      thnAxisCentOrMult = {thnConfigAxisCent, "Centrality (%)"};
    } else {
      // pp: use multiplicity
      thnAxisCentOrMult = {thnConfigAxisMult, "FT0C Multiplicity"};
    }

    auto hPairAncestor = registry.add<THn>("PCM/hPairAncestor", "minv vs pT vs cent vs pair truth vs closest common ancestor", HistType::kTHnF, {thnAxisInvMass, thnAxisPtRec, thnAxisCentOrMult, thnAxisPairTruth, thnAxisPairBucket});
    static constexpr std::array<std::string_view, static_cast<int>(PhotonTruth::kNTruth)> TruthLabels = {"true #gamma", "non #gamma", "fake"};
    for (std::size_t i = 0; i < TruthLabels.size(); ++i) {
      hPairAncestor->GetAxis(3)->SetBinLabel(static_cast<int>(i) + 1, TruthLabels[i].data());
    }
    static constexpr std::array<std::string_view, static_cast<int>(PairBucket::kNBuckets)> BucketLabels = {"none", "same particle", "#pi^{0} direct", "#eta direct", "#pi^{0} cross", "#eta cross", "#omega", "K^{0}_{S}", "#eta'", "other hadron", "parton/string"};
    for (std::size_t i = 0; i < BucketLabels.size(); ++i) {
      hPairAncestor->GetAxis(4)->SetBinLabel(static_cast<int>(i) + 1, BucketLabels[i].data());
    }
    // mixed event, same axes as hPairAncestor without the bucket (pairs from different events share no ancestor)
    auto hMixed = registry.add<THn>("PCM/hMixed", "minv vs pT vs cent vs pair truth, mixed event", HistType::kTHnF, {thnAxisInvMass, thnAxisPtRec, thnAxisCentOrMult, thnAxisPairTruth});
    for (std::size_t i = 0; i < TruthLabels.size(); ++i) {
      hMixed->GetAxis(3)->SetBinLabel(static_cast<int>(i) + 1, TruthLabels[i].data());
    }
    registry.addClone("PCM/", "EMC/");
    registry.addClone("PCM/", "PCMEMC/");

    // ancestor chain QA
    registry.add("QA/PCM/hChainLength", "length of the ancestor chain;chain length;counts", HistType::kTH1D, {{kMaxChainLength, 0.5, kMaxChainLength + 0.5}});
    auto hChainEnd = registry.add<TH1>("QA/PCM/hChainEnd", "reason the ancestor chain ended;;counts", HistType::kTH1D, {{4, -0.5, 3.5}});
    hChainEnd->GetXaxis()->SetBinLabel(1, "parton-like");
    hChainEnd->GetXaxis()->SetBinLabel(2, "beam nucleus");
    hChainEnd->GetXaxis()->SetBinLabel(3, "no mothers");
    hChainEnd->GetXaxis()->SetBinLabel(4, "truncated");
    auto hChainEndKind = registry.add<TH1>("QA/PCM/hChainEndKind", "kind of the last chain element;;counts", HistType::kTH1D, {{static_cast<int>(AncestorKind::NClasses), -0.5, static_cast<int>(AncestorKind::NClasses) - 0.5}});
    static constexpr std::array<std::string_view, static_cast<int>(AncestorKind::NClasses)> KindLabels = {"none", "string/parton", "beam/nucleus", "#pi^{0}", "#eta", "other meson", "baryon", "#gamma", "e^{#pm}", "other"};
    for (std::size_t i = 0; i < KindLabels.size(); ++i) {
      hChainEndKind->GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, KindLabels[i].data());
    }
    registry.add("QA/PCM/hChainEndPdg", "|PDG| of the last chain element;|PDG| (last bin: #geq 9999, incl. nuclei);counts", HistType::kTH1D, {{10000, -0.5, 9999.5}});
    registry.addClone("QA/PCM/", "QA/EMC/");

    // origin of PCM-EMC same particle pairs and why track matching did not remove them
    if (doprocessPcmEmcal.value) {
      const AxisSpec axisOrigin{static_cast<int>(SameParticleOrigin::kNOrigins), -0.5, static_cast<int>(SameParticleOrigin::kNOrigins) - 0.5, "leading contributor of the cluster"};
      const AxisSpec axisStatus{static_cast<int>(LegMatchStatus::kNStatus), -0.5, static_cast<int>(LegMatchStatus::kNStatus) - 0.5, "V0 leg among the matched tracks"};
      const AxisSpec axisPtLeg{40, 0., 20., "#it{p}_{T,leg} (GeV/#it{c})"};
      const AxisSpec axisDeltaEta{100, -0.1, 0.1, "#Delta#eta (cluster - track)"};
      const AxisSpec axisDeltaPhi{100, -0.1, 0.1, "#Delta#varphi (cluster - track)"};
      const AxisSpec axisDp{200, -1., 1., "(#it{p}_{match} - #it{p}_{leg}) / #it{p}_{leg}"};
      static constexpr std::array<std::string_view, static_cast<int>(SameParticleOrigin::kNOrigins)> OriginLabels = {"leg e^{#pm}", "brems #gamma of leg", "other e^{#pm}", "other"};
      static constexpr std::array<std::string_view, static_cast<int>(LegMatchStatus::kNStatus)> StatusLabels = {"no match stored", "no match with leg #it{p}", "outside window", "in window, #it{E}/#it{p} #geq cut", "in window"};

      auto hOrigin = registry.add<TH3>("PCMEMC/SameParticle/hOrigin", "same particle pairs: leading contributor of the cluster", HistType::kTH3F, {thnAxisInvMass, thnAxisPtRec, axisOrigin});
      for (std::size_t i = 0; i < OriginLabels.size(); ++i) {
        hOrigin->GetZaxis()->SetBinLabel(static_cast<int>(i) + 1, OriginLabels[i].data());
      }
      registry.add("PCMEMC/SameParticle/hLegEoverP", "cluster from a V0 leg: #it{E}_{cls}/#it{p}_{leg}", HistType::kTH2F, {{200, 0., 4., "#it{E}_{cls}/#it{p}_{leg}"}, axisPtLeg});
      for (const auto& [type, title] : {std::pair<std::string_view, std::string_view>{"Prim", "primary"}, {"Sec", "secondary"}}) {
        auto hStatus = registry.add<TH2>(std::format("PCMEMC/SameParticle/hLeg{}Status", type.data()).data(), std::format("cluster from a V0 leg: leg among the matched {} tracks", title.data()).data(), HistType::kTH2F, {axisStatus, axisPtLeg});
        for (std::size_t i = 0; i < StatusLabels.size(); ++i) {
          hStatus->GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, StatusLabels[i].data());
        }
        registry.add(std::format("PCMEMC/SameParticle/hLeg{}Dp", type.data()).data(), std::format("cluster from a V0 leg: closest {} match in #it{{p}}", title.data()).data(), HistType::kTH1F, {axisDp});
        registry.add(std::format("PCMEMC/SameParticle/hLeg{}Residual", type.data()).data(), std::format("cluster from a V0 leg: residual of the {} match", title.data()).data(), HistType::kTH3F, {axisDeltaEta, axisDeltaPhi, axisPtLeg});
      }
    }

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
    float centOrMult = getCentralityOrMultiplicity(collision);
    if (useCent && (centOrMult < eventcuts.cfgMinCent || centOrMult > eventcuts.cfgMaxCent)) {
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

  /// \brief Reset the per DF cache of one photon type
  /// \param chains ancestor chains, set to -1 (not {}, which gives 0, a valid index)
  /// \param truth truth flags, set to kFake
  /// \param size size of the full photon table
  static void resetCache(std::vector<Chain>& chains, std::vector<PhotonTruth>& truth, const std::size_t size)
  {
    Chain empty;
    empty.fill(-1);
    chains.assign(size, empty);
    truth.assign(size, PhotonTruth::kFake);
  }

  /// \brief Fill the ancestor chain QA
  /// \param mcIter MC particle iterator, has to sit where getAncestorChain left it (last chain element or the nucleus that ended the chain)
  /// \param length length of the chain as returned by getAncestorChain
  template <int photonType, o2::soa::is_iterator TMCIter>
  void fillChainQA(TMCIter const& mcIter, const int length)
  {
    static constexpr std::array<std::string_view, 2> HistLength = {"QA/PCM/hChainLength", "QA/EMC/hChainLength"};
    static constexpr std::array<std::string_view, 2> HistEnd = {"QA/PCM/hChainEnd", "QA/EMC/hChainEnd"};
    static constexpr std::array<std::string_view, 2> HistEndKind = {"QA/PCM/hChainEndKind", "QA/EMC/hChainEndKind"};
    static constexpr std::array<std::string_view, 2> HistEndPdg = {"QA/PCM/hChainEndPdg", "QA/EMC/hChainEndPdg"};
    if (length <= 0) {
      return;
    }
    static constexpr int MaxPdgBin = 9999; // nuclei and other large codes go to the last bin
    const int lastPdg = mcIter.pdgCode();
    int endReason = 3; // chain is full without reaching a parton, a nucleus or a particle without mothers
    if (isPartonLike(lastPdg)) {
      endReason = 0;
    } else if (isNucleus(lastPdg)) {
      endReason = 1;
    } else if (!mcIter.has_mothers()) {
      endReason = 2;
    }
    registry.fill(HIST(HistLength[photonType]), length);
    registry.fill(HIST(HistEnd[photonType]), endReason);
    registry.fill(HIST(HistEndKind[photonType]), static_cast<int>(ancestorKindOfPdg(lastPdg)));
    registry.fill(HIST(HistEndPdg[photonType]), std::min(std::abs(lastPdg), MaxPdgBin));
  }

  /// \brief Fill ancestor chains and truth flags of all EMC clusters of the DF
  /// \param clusters full (unfiltered) cluster table
  /// \param mcIter iterator of the full MC particle table, WILL BE CHANGED
  /// \param doChainQA fill the ancestor chain QA (off for the mixed event, so chains are not counted twice)
  template <o2::soa::is_table TClusters, o2::soa::is_iterator TMCIter>
  void fillEmcCache(TClusters const& clusters, TMCIter& mcIter, const bool doChainQA = true)
  {
    for (const auto& cluster : clusters) {
      const auto mcIds = cluster.emmcparticleIds();
      if (mcIds.empty() || mcIds[0] < 0) {
        continue; // noise cluster without MC label: kFake, chain stays -1
      }
      auto& chain = mChainsEMC[cluster.globalIndex()];
      const int length = getAncestorChain(mcIter, mcIds[0], chain); // leading contributor
      if (doChainQA) {
        fillChainQA<PhotonEMC>(mcIter, length);
      }
      if (length <= 0) {
        continue; // leading contributor is a nucleus: kFake, chain stays -1
      }

      // e+- whose mother is a photon (conversion in front of the EMCal): start the chain at the photon,
      // otherwise the direct mother check of the pair classification fails for converted photons
      mcIter.setCursor(chain[0]);
      if (length > 1 && std::abs(mcIter.pdgCode()) == PDG_t::kElectron) {
        mcIter.setCursor(chain[1]);
        if (mcIter.pdgCode() == PDG_t::kGamma) {
          std::shift_left(chain.begin(), chain.end(), 1);
          chain.back() = -1;
        }
      }
      mcIter.setCursor(chain[0]);
      mTruthEMC[cluster.globalIndex()] = mcIter.pdgCode() == PDG_t::kGamma ? PhotonTruth::kTruePhoton : PhotonTruth::kNonPhoton;
    }
  }

  /// \brief Fill ancestor chains and truth flags of all V0 photons of the DF
  /// \param photons full (unfiltered) V0 photon table
  /// \param legs full V0 leg table with MC labels
  /// \param mcParticles full MC particle table the leg labels point into
  /// \param doChainQA fill the ancestor chain QA (off for the mixed event, so chains are not counted twice)
  template <o2::soa::is_table TPhotons, o2::soa::is_table TLegs, o2::soa::is_table TMCParticles>
  void fillPcmCache(TPhotons const& photons, TLegs const& legs, TMCParticles const& mcParticles, const bool doChainQA = true)
  {
    auto mcIter = mcParticles.begin();
    auto posLeg = legs.begin();
    auto negLeg = legs.begin();
    auto posMc = mcParticles.begin();
    auto negMc = mcParticles.begin();

    Chain chainPos{};
    Chain chainNeg{};

    for (const auto& photon : photons) {
      posLeg.setCursor(photon.posTrackId());
      negLeg.setCursor(photon.negTrackId());
      if (posLeg.emmcparticleId() < 0 || negLeg.emmcparticleId() < 0) {
        continue; // leg without MC label: kFake, chain stays -1
      }
      posMc.setCursor(posLeg.emmcparticleId());
      negMc.setCursor(negLeg.emmcparticleId());

      auto& chain = mChainsPCM[photon.globalIndex()];
      auto& truth = mTruthPCM[photon.globalIndex()];
      const int photonId = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(posMc, negMc, PDG_t::kPositron, PDG_t::kElectron, PDG_t::kGamma, mcParticles);
      if (photonId >= 0) {
        truth = PhotonTruth::kTruePhoton;
        const int length = getAncestorChain(mcIter, photonId, chain);
        if (doChainQA) {
          fillChainQA<PhotonPCM>(mcIter, length);
        }
        continue;
      }

      // not a true conversion: K0S/Lambda V0, Dalitz pair or combinatorial fake
      const int lengthPos = getAncestorChain(mcIter, posLeg.emmcparticleId(), chainPos);
      if (doChainQA) {
        fillChainQA<PhotonPCM>(mcIter, lengthPos); // the kept chain is a tail of this one, so it ends the same way
      }
      getAncestorChain(mcIter, negLeg.emmcparticleId(), chainNeg);
      const int ancestor = getFirstCommonAncestor(chainPos, chainNeg);
      if (ancestor < 0) {
        continue; // no common ancestor: kFake, chain stays -1
      }
      const bool isDirectMother = chainPos[1] == ancestor && chainNeg[1] == ancestor;
      truth = isDirectMother ? PhotonTruth::kNonPhoton : PhotonTruth::kFake;
      const auto start = std::ranges::find(chainPos, ancestor);
      std::ranges::fill(std::copy(start, chainPos.end(), chain.begin()), chain.end(), -1);
    }
  }

  /// \brief Classify a pair by the closest common ancestor of the two candidates
  /// \param chain1 ancestor chain of the first candidate
  /// \param chain2 ancestor chain of the second candidate
  /// \param mcIter iterator of the full MC particle table, WILL BE CHANGED
  template <o2::soa::is_iterator TMCIter>
  PairBucket classifyPair(Chain const& chain1, Chain const& chain2, TMCIter& mcIter)
  {
    const int ancestor = getFirstCommonAncestor(chain1, chain2);
    if (ancestor < 0) {
      return PairBucket::kNone;
    }
    if (ancestor == chain1[0] || ancestor == chain2[0]) {
      return PairBucket::kSameParticle;
    }
    mcIter.setCursor(ancestor);
    const int absPdg = std::abs(mcIter.pdgCode());
    const bool isDirectMother = chain1[1] == ancestor && chain2[1] == ancestor;
    if (isDirectMother && absPdg == PDG_t::kPi0) {
      return PairBucket::kSignalPi0;
    }
    if (isDirectMother && absPdg == o2::constants::physics::Pdg::kEta) {
      return PairBucket::kSignalEta;
    }
    switch (absPdg) {
      case PDG_t::kPi0:
        return PairBucket::kCrossPi0;
      case o2::constants::physics::Pdg::kEta:
        return PairBucket::kCrossEta;
      case o2::constants::physics::Pdg::kOmega:
        return PairBucket::kCrossOmega;
      case PDG_t::kK0Short:
        return PairBucket::kCrossK0S;
      case o2::constants::physics::Pdg::kEtaPrime:
        return PairBucket::kCrossEtaPrime;
      default:
        break;
    }
    if (isPartonLike(absPdg)) {
      return PairBucket::kPartonString;
    }
    return PairBucket::kCrossOtherHadron;
  }

  /// \brief Look for the V0 leg among the tracks matched to an EMC cluster and classify why it did not veto the cluster.
  ///        The matched track tables store no track index, so the leg is identified by its momentum.
  /// \tparam isSecondary matches are secondaries (V0 legs, no E/p condition) instead of primaries
  /// \param matches tracks matched to this cluster
  /// \param eCluster cluster energy
  /// \param pLeg momentum of the V0 leg
  /// \param ptLeg transverse momentum of the V0 leg
  template <bool isSecondary, o2::soa::is_table TMatches>
  void fillLegMatch(TMatches const& matches, const float eCluster, const float pLeg, const float ptLeg)
  {
    static constexpr std::array<std::string_view, 2> HistStatus = {"PCMEMC/SameParticle/hLegPrimStatus", "PCMEMC/SameParticle/hLegSecStatus"};
    static constexpr std::array<std::string_view, 2> HistDp = {"PCMEMC/SameParticle/hLegPrimDp", "PCMEMC/SameParticle/hLegSecDp"};
    static constexpr std::array<std::string_view, 2> HistResidual = {"PCMEMC/SameParticle/hLegPrimResidual", "PCMEMC/SameParticle/hLegSecResidual"};
    static constexpr float MaxRelDp = 0.1f; // max |p_match - p_leg| / p_leg to identify the match as the leg

    LegMatchStatus status = LegMatchStatus::kNoMatchStored;
    if (matches.size() > 0) {
      float bestDp = 1e9f;
      float dEta = 0.f, dPhi = 0.f, trackPt = 0.f, trackP = 0.f;
      for (const auto& match : matches) {
        const float dp = (match.trackP() - pLeg) / pLeg;
        if (std::fabs(dp) < std::fabs(bestDp)) {
          bestDp = dp;
          dEta = match.deltaEta();
          dPhi = match.deltaPhi();
          trackPt = match.trackPt();
          trackP = match.trackP();
        }
      }
      registry.fill(HIST(HistDp[isSecondary]), bestDp);
      if (std::fabs(bestDp) > MaxRelDp) {
        status = LegMatchStatus::kNoCompatibleMatch;
      } else {
        registry.fill(HIST(HistResidual[isSecondary]), dEta, dPhi, ptLeg);
        // same windows as EMCPhotonCut::checkTrackMatching
        bool isInWindow = false;
        if constexpr (isSecondary) {
          isInWindow = std::fabs(dEta) <= fEMCCut.GetSecTrackMatchingEta(trackPt) && std::fabs(dPhi) <= fEMCCut.GetSecTrackMatchingPhi(trackPt);
        } else {
          isInWindow = std::fabs(dEta) <= fEMCCut.GetTrackMatchingEta(trackPt) && std::fabs(dPhi) <= fEMCCut.GetTrackMatchingPhi(trackPt);
        }
        if (!isInWindow) {
          status = LegMatchStatus::kOutsideWindow;
        } else if (!isSecondary && eCluster / trackP >= emccuts.cfgEMCEoverp) {
          status = LegMatchStatus::kInsideWindowEoverP;
        } else {
          status = LegMatchStatus::kInsideWindow;
        }
      }
    }
    registry.fill(HIST(HistStatus[isSecondary]), static_cast<int>(status), ptLeg);
  }

  /// \brief Find out where a PCM-EMC same particle pair comes from: what made the EMC cluster and why track matching did not remove it
  /// \param pcm V0 photon
  /// \param emc EMC cluster
  /// \param legs full V0 leg table with MC labels
  /// \param matchedPrims full table of primary tracks matched to clusters
  /// \param matchedSeconds full table of secondary tracks (V0 legs) matched to clusters
  /// \param mcIter iterator of the full MC particle table, WILL BE CHANGED
  template <o2::soa::is_iterator TPCM, o2::soa::is_iterator TEMC, o2::soa::is_table TLegs, o2::soa::is_table TMatchedPrims, o2::soa::is_table TMatchedSecs, o2::soa::is_iterator TMCIter>
  void diagnoseSameParticle(TPCM const& pcm, TEMC const& emc, TLegs const& legs, TMatchedPrims const& matchedPrims, TMatchedSecs const& matchedSeconds, TMCIter& mcIter)
  {
    const auto mcIds = emc.emmcparticleIds();
    if (mcIds.empty() || mcIds[0] < 0) {
      return;
    }
    const int lead = mcIds[0]; // leading contributor, before the conversion shift of the chain
    const auto posLeg = legs.iteratorAt(pcm.posTrackId());
    const auto negLeg = legs.iteratorAt(pcm.negTrackId());
    const std::array<int, 2> legMcIds = {posLeg.emmcparticleId(), negLeg.emmcparticleId()};

    SameParticleOrigin origin = SameParticleOrigin::kOther;
    int legIndex = -1;
    if (lead == legMcIds[0] || lead == legMcIds[1]) {
      origin = SameParticleOrigin::kLegElectron;
      legIndex = lead == legMcIds[0] ? 0 : 1;
    } else {
      mcIter.setCursor(lead);
      const int absPdg = std::abs(mcIter.pdgCode());
      const int motherId = mcIter.has_mothers() ? mcIter.mothersIds()[0] : -1;
      if (absPdg == PDG_t::kGamma && motherId >= 0 && (motherId == legMcIds[0] || motherId == legMcIds[1])) {
        origin = SameParticleOrigin::kBremsFromLeg;
      } else if (absPdg == PDG_t::kElectron) {
        origin = SameParticleOrigin::kOtherElectron;
      }
    }

    ROOT::Math::PtEtaPhiMVector v1(pcm.pt(), pcm.eta(), pcm.phi(), 0.);
    ROOT::Math::PtEtaPhiMVector v2(emc.pt(), emc.eta(), emc.phi(), 0.);
    ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;
    registry.fill(HIST("PCMEMC/SameParticle/hOrigin"), vMeson.M(), vMeson.Pt(), static_cast<int>(origin));
    if (origin != SameParticleOrigin::kLegElectron) {
      return;
    }

    const auto& leg = legIndex == 0 ? posLeg : negLeg;
    const float pLeg = std::hypot(leg.px(), leg.py(), leg.pz());
    const float ptLeg = std::hypot(leg.px(), leg.py());
    if (pLeg <= 0.f) {
      return;
    }
    registry.fill(HIST("PCMEMC/SameParticle/hLegEoverP"), emc.e() / pLeg, ptLeg);
    fillLegMatch<false>(matchedPrims.sliceBy(perEMCClusterMT, emc.globalIndex()), emc.e(), pLeg, ptLeg);
    fillLegMatch<true>(matchedSeconds.sliceBy(perEMCClusterMS, emc.globalIndex()), emc.e(), pLeg, ptLeg);
  }

  /// \brief Apply the pair cuts and fill the pair ancestor histogram
  /// \return bucket of the pair, std::nullopt if the pair was rejected
  template <int pairType, o2::soa::is_iterator TPhoton1, o2::soa::is_iterator TPhoton2, o2::soa::is_iterator TMCIter>
  std::optional<PairBucket> fillPair(TPhoton1 const& g1, TPhoton2 const& g2, Chain const& chain1, Chain const& chain2, const PhotonTruth truth1, const PhotonTruth truth2, const float centOrMult, TMCIter& mcIter)
  {
    static constexpr std::array<std::string_view, 3> HistPairAncestor = {"PCM/hPairAncestor", "EMC/hPairAncestor", "PCMEMC/hPairAncestor"};

    ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
    ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
    ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;
    float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));

    registry.fill(HIST("hMesonCuts"), 1);
    if (openingAngle <= mesonConfig.minOpenAngle) {
      registry.fill(HIST("hMesonCuts"), 2);
      return std::nullopt;
    }
    registry.fill(HIST("hMesonCuts"), 6);

    const PhotonTruth pairTruth = std::max(truth1, truth2);
    const PairBucket bucket = classifyPair(chain1, chain2, mcIter);
    registry.fill(HIST(HistPairAncestor[pairType]), vMeson.M(), vMeson.Pt(), centOrMult, static_cast<int>(pairTruth), static_cast<int>(bucket));
    if (mesonConfig.cfgEnableQA.value) {
      registry.fill(HIST("mesonQA/hInvMassPt"), vMeson.M(), vMeson.Pt());
    }
    return bucket;
  }

  /// \brief eta a cluster would have if its own vertex vzOld is swapped for vzNew, same as in taskPi0FlowEMC
  /// \param etaOld old eta value
  /// \param vzOld old primary vertex z position
  /// \param vzNew new primary vertex z position
  /// \param emcalR radius of the EMCal
  static float correctEtaForVertexShift(float etaOld, float vzOld, float vzNew, float emcalR)
  {
    return std::asinh(std::sinh(etaOld) + (vzOld - vzNew) / emcalR);
  }

  /// \brief Apply the pair cuts to a mixed event pair and fill the mixed event histogram
  /// \param eta2 eta of the second photon, possibly re-projected onto the vertex of the first collision
  template <int pairType, o2::soa::is_iterator TPhoton1, o2::soa::is_iterator TPhoton2>
  void fillMixedPair(TPhoton1 const& g1, TPhoton2 const& g2, const float eta2, const PhotonTruth truth1, const PhotonTruth truth2, const float centOrMult)
  {
    static constexpr std::array<std::string_view, 3> HistMixed = {"PCM/hMixed", "EMC/hMixed", "PCMEMC/hMixed"};

    ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
    ROOT::Math::PtEtaPhiMVector v2(g2.pt(), eta2, g2.phi(), 0.);
    ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;
    float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));
    if (openingAngle <= mesonConfig.minOpenAngle) {
      return;
    }
    registry.fill(HIST(HistMixed[pairType]), vMeson.M(), vMeson.Pt(), centOrMult, static_cast<int>(std::max(truth1, truth2)));
  }

  // PCM-PCM same event
  void processPcm(Colls const& collisions, PcmPhotons const& photons, PcmMcLegs const& legs, McParticles const& mcParticles)
  {
    if (photons.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    // none of the photon tables is filtered, so its size is the size of the full table the globalIndex points into
    resetCache(mChainsPCM, mTruthPCM, photons.size());
    fillPcmCache(photons, legs, mcParticles);

    EMBitFlags v0flags(photons.size());
    fV0PhotonCut.AreSelectedRunning<decltype(photons), PcmMcLegs>(v0flags, photons, &registry);

    auto mcIter = mcParticles.begin();
    for (const auto& collision : collisions) {
      initCCDB(collision);
      if (!isFullEventSelected(collision, true)) {
        continue;
      }
      const float centOrMult = getCentralityOrMultiplicity(collision);
      auto photonsPerCollision = photons.sliceBy(perCollisionPCM, collision.globalIndex());

      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photonsPerCollision, photonsPerCollision))) {
        if (!(v0flags.test(g1.globalIndex())) || !(v0flags.test(g2.globalIndex()))) {
          continue;
        }
        fillPair<PairPCMPCM>(g1, g2, mChainsPCM[g1.globalIndex()], mChainsPCM[g2.globalIndex()], mTruthPCM[g1.globalIndex()], mTruthPCM[g2.globalIndex()], centOrMult, mcIter);
      }
    }
  }
  PROCESS_SWITCH(TwoPhotonMcTask, processPcm, "Process for pcm photons", false);

  // EMCal-EMCal same event
  void processEmcal(Colls const& collisions, EMCalPhotons const& clusters, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds, McParticles const& mcParticles)
  {
    if (clusters.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    auto mcIter = mcParticles.begin();
    // none of the photon tables is filtered, so its size is the size of the full table the globalIndex points into
    resetCache(mChainsEMC, mTruthEMC, clusters.size());
    fillEmcCache(clusters, mcIter);

    EMBitFlags emcFlags(clusters.size());
    fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds, &registry);

    for (const auto& collision : collisions) {
      initCCDB(collision);
      if (!isFullEventSelected(collision, true)) {
        continue;
      }
      const float centOrMult = getCentralityOrMultiplicity(collision);
      auto photonsPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());

      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photonsPerCollision, photonsPerCollision))) {
        if (!(emcFlags.test(g1.globalIndex())) || !(emcFlags.test(g2.globalIndex()))) {
          continue;
        }
        fillPair<PairEMCEMC>(g1, g2, mChainsEMC[g1.globalIndex()], mChainsEMC[g2.globalIndex()], mTruthEMC[g1.globalIndex()], mTruthEMC[g2.globalIndex()], centOrMult, mcIter);
      }
    }
  }
  PROCESS_SWITCH(TwoPhotonMcTask, processEmcal, "Process for emcal photons", false);

  // PCM-EMCal same event
  void processPcmEmcal(Colls const& collisions, EMCalPhotons const& clusters, PcmPhotons const& photons, PcmMcLegs const& legs, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds, McParticles const& mcParticles)
  {
    if (clusters.size() <= 0 || photons.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    auto mcIter = mcParticles.begin();
    // none of the photon tables is filtered, so its size is the size of the full table the globalIndex points into
    resetCache(mChainsPCM, mTruthPCM, photons.size());
    resetCache(mChainsEMC, mTruthEMC, clusters.size());
    fillPcmCache(photons, legs, mcParticles);
    fillEmcCache(clusters, mcIter);

    EMBitFlags emcFlags(clusters.size());
    fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds, &registry);
    EMBitFlags v0flags(photons.size());
    fV0PhotonCut.AreSelectedRunning<decltype(photons), PcmMcLegs>(v0flags, photons, &registry);

    for (const auto& collision : collisions) {
      initCCDB(collision);
      if (!isFullEventSelected(collision, true)) {
        continue;
      }
      const float centOrMult = getCentralityOrMultiplicity(collision);
      auto photonsEMCPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());
      auto photonsPCMPerCollision = photons.sliceBy(perCollisionPCM, collision.globalIndex());

      // PCM - EMC loop, two different tables: all combinations
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photonsPCMPerCollision, photonsEMCPerCollision))) {
        if (!(v0flags.test(g1.globalIndex())) || !(emcFlags.test(g2.globalIndex()))) {
          continue;
        }
        const auto bucket = fillPair<PairPCMEMC>(g1, g2, mChainsPCM[g1.globalIndex()], mChainsEMC[g2.globalIndex()], mTruthPCM[g1.globalIndex()], mTruthEMC[g2.globalIndex()], centOrMult, mcIter);
        if (bucket == PairBucket::kSameParticle) {
          diagnoseSameParticle(g1, g2, legs, matchedPrims, matchedSeconds, mcIter);
        }
      }
    }
  }
  PROCESS_SWITCH(TwoPhotonMcTask, processPcmEmcal, "Process for pcm and emcal photons", true);

  // PCM-PCM mixed event
  void processPcmMixed(Colls const& collisions, PcmPhotons const& photons, PcmMcLegs const& legs, McParticles const& mcParticles)
  {
    if (photons.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    // the mixed pairs are built inside the DF, so the per DF cache is valid for them
    resetCache(mChainsPCM, mTruthPCM, photons.size());
    fillPcmCache(photons, legs, mcParticles, false);

    EMBitFlags v0flags(photons.size());
    fV0PhotonCut.AreSelectedRunning<decltype(photons), PcmMcLegs>(v0flags, photons);

    BinningType binningMixedEvent{{mixingConfig.cfgVtxBins, mixingConfig.cfgCentBins}};
    auto photonTuple = std::make_tuple(photons);
    SameKindPair<Colls, PcmPhotons, BinningType> pair{binningMixedEvent, mixingConfig.cfgMixingDepth, -1, collisions, photonTuple, &cache};

    for (const auto& [c1, photons1, c2, photons2] : pair) {
      if (!isFullEventSelected(c1) || !isFullEventSelected(c2)) {
        continue;
      }
      const float centOrMult = getCentralityOrMultiplicity(c1);
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1, photons2))) {
        if (!(v0flags.test(g1.globalIndex())) || !(v0flags.test(g2.globalIndex()))) {
          continue;
        }
        fillMixedPair<PairPCMPCM>(g1, g2, g2.eta(), mTruthPCM[g1.globalIndex()], mTruthPCM[g2.globalIndex()], centOrMult);
      }
    }
  }
  PROCESS_SWITCH(TwoPhotonMcTask, processPcmMixed, "Process for pcm photons, mixed event", false);

  // EMCal-EMCal mixed event
  void processEmcalMixed(Colls const& collisions, EMCalPhotons const& clusters, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds, McParticles const& mcParticles)
  {
    if (clusters.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    auto mcIter = mcParticles.begin();
    // the mixed pairs are built inside the DF, so the per DF cache is valid for them
    resetCache(mChainsEMC, mTruthEMC, clusters.size());
    fillEmcCache(clusters, mcIter, false);

    EMBitFlags emcFlags(clusters.size());
    fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds);

    BinningType binningMixedEvent{{mixingConfig.cfgVtxBins, mixingConfig.cfgCentBins}};
    auto clusterTuple = std::make_tuple(clusters);
    SameKindPair<Colls, EMCalPhotons, BinningType> pair{binningMixedEvent, mixingConfig.cfgMixingDepth, -1, collisions, clusterTuple, &cache};

    for (const auto& [c1, clusters1, c2, clusters2] : pair) {
      if (!isFullEventSelected(c1) || !isFullEventSelected(c2)) {
        continue;
      }
      const float centOrMult = getCentralityOrMultiplicity(c1);
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(clusters1, clusters2))) {
        if (!(emcFlags.test(g1.globalIndex())) || !(emcFlags.test(g2.globalIndex()))) {
          continue;
        }
        // cluster 2 seen from the vertex of collision 1, as in taskPi0FlowEMC
        const float eta2 = mixingConfig.cfgCorrectMixedVtxEta.value ? correctEtaForVertexShift(g2.eta(), c2.posZ(), c1.posZ(), mixingConfig.cfgEmcalEffRadius.value) : g2.eta();
        fillMixedPair<PairEMCEMC>(g1, g2, eta2, mTruthEMC[g1.globalIndex()], mTruthEMC[g2.globalIndex()], centOrMult);
      }
    }
  }
  PROCESS_SWITCH(TwoPhotonMcTask, processEmcalMixed, "Process for emcal photons, mixed event", false);

  // PCM-EMCal mixed event
  void processPcmEmcalMixed(Colls const& collisions, EMCalPhotons const& clusters, PcmPhotons const& photons, PcmMcLegs const& legs, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds, McParticles const& mcParticles)
  {
    if (clusters.size() <= 0 || photons.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    auto mcIter = mcParticles.begin();
    // the mixed pairs are built inside the DF, so the per DF cache is valid for them
    resetCache(mChainsPCM, mTruthPCM, photons.size());
    resetCache(mChainsEMC, mTruthEMC, clusters.size());
    fillPcmCache(photons, legs, mcParticles, false);
    fillEmcCache(clusters, mcIter, false);

    EMBitFlags emcFlags(clusters.size());
    fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds);
    EMBitFlags v0flags(photons.size());
    fV0PhotonCut.AreSelectedRunning<decltype(photons), PcmMcLegs>(v0flags, photons);

    BinningType binningMixedEvent{{mixingConfig.cfgVtxBins, mixingConfig.cfgCentBins}};
    auto associatedTables = std::make_tuple(photons, clusters);
    Pair<Colls, PcmPhotons, EMCalPhotons, BinningType> pair{binningMixedEvent, mixingConfig.cfgMixingDepth, -1, collisions, associatedTables, &cache};

    for (const auto& [c1, photonsPCM, c2, photonsEMC] : pair) {
      if (!isFullEventSelected(c1) || !isFullEventSelected(c2)) {
        continue;
      }
      const float centOrMult = getCentralityOrMultiplicity(c1);
      for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photonsPCM, photonsEMC))) {
        if (!(v0flags.test(g1.globalIndex())) || !(emcFlags.test(g2.globalIndex()))) {
          continue;
        }
        fillMixedPair<PairPCMEMC>(g1, g2, g2.eta(), mTruthPCM[g1.globalIndex()], mTruthEMC[g2.globalIndex()], centOrMult);
      }
    }
  }
  PROCESS_SWITCH(TwoPhotonMcTask, processPcmEmcalMixed, "Process for pcm and emcal photons, mixed event", false);

}; // End struct TwoPhotonMcTask

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<TwoPhotonMcTask>(context)};
}
