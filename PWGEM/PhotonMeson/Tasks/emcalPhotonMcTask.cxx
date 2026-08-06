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

/// \file emcalPhotonMcTask.cxx
/// \brief Analysis task for to analyse photon clusters and differences to conversion clusters in MC
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/PhotonMeson/Core/EMBitFlags.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/DataModel/ConversionMl.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

#include "Common/Core/RecoDecay.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
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

#include <Math/Vector3D.h> // IWYU pragma: keep (do not replace with Math/Vector3Dfwd.h)
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TPDGCode.h>
#include <TTree.h>

#include <sys/types.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;
using namespace o2::constants::physics;
using namespace o2::aod::pwgem::photonmeson::utils::mcutil;

constexpr float MinAmpThreshold = 0.2f; // Minimum cluster amplitude threshold to count as significant
constexpr float MaxAmpDiff = 0.1f;      // Maximum cluster amplitude difference to leading cluster contribution to count as significant

enum CentralityEstimator {
  None = 0,
  CFT0A,
  CFT0C,
  CFT0M,
  NCentralityEstimators
};

enum class TruthClass {
  Conversion = 0, // 0: true e+/e- pair from the same conversion

  PhotonPairSamePi0, // 1: two photon clusters, same Pi0
  PhotonPairDiffPi0, // 2: two photon clusters, different Pi0s
  PhotonPairOnePi0,  // 3: two photon clusters, only one from a Pi0

  PhotonElectronSamePi0, // 4: photon + electron cluster, same Pi0
  PhotonElectronDiffPi0, // 5: photon + electron cluster, different Pi0s
  PhotonElectronOnePi0,  // 6: photon + electron cluster, only one from a Pi0
  BSPhotonElectron,      // 7: photon + electron cluster, from Bremsstrahlung

  ElectronPairSamePi0, // 8: e+e cluster pair, same Pi0 (conversion and/or Dalitz)
  ElectronPairDiffPi0, // 9: e+e cluster pair, different Pi0s
  ElectronPairOnePi0,  // 10: e+e cluster pair, only one from a Pi0

  SplitPhotonCluster,   // 11: one photon producing two clusters
  SplitLeptonCluster,   // 12: one lepton producing two clusters
  PhotonBSPhotonPair,   // 13: photon + photon from Bremsstrahlung
  ElectronBSPhotonPair, // 14: one cluster from Bremsstrahlung and one electron cluster except case BSPhotonElectron
  BSPhotonPair,         // 15: both photons from Bremsstrahlung

  Background, // 16: else / uncorrelated

  NClasses
};

enum class ClusterTruthClass {
  Photon = 0,    // cluster is a true, unconverted photon
  Conversion,    // cluster is a leg of a true conversion
  Pi0Conversion, // ...and that conversion's photon itself came from pi0/eta (excludes same-gamma-pair case)
  Pi0,           // cluster comes from the pi0/eta decay chain (not via the conversion-leg path above)
  Background,    // cluster is flagged as both meson-chain AND conversion -- ambiguous/contradictory
  NClasses
};

enum class TagDecision {
  NotTagged = 0,
  Tagged,
  NTags
};

struct ClusterMcInfo {
  bool isLepton = false;
  bool isPhoton = false;
  bool isFromConv = false;
  bool isMergedConv = false;
  bool isFromPi0 = false;
  bool isFromBremsstrahlung = false;
  int convMotherId = -1;
  int photonId = -1;
  float purity = 0;
};

template <o2::soa::is_iterator TGroup, o2::soa::is_iterator TIter, o2::soa::is_table McParticles>
ClusterMcInfo classifyCluster(const TGroup& g, TIter& mcCluster, TIter& mcClusterLooper, TIter& mcClusterLooper2, McParticles const& mcParticles)
{
  ClusterMcInfo info;
  mcCluster.setCursor(g.emmcparticleIds()[0]);
  info.isFromBremsstrahlung = isFromBremsstrahlung(mcCluster, mcClusterLooper); // particle has to be a photon and it has to have a e+ or e- as mother!
  float leadingAmplitude = g.amplitude()[0];
  info.purity = leadingAmplitude;
  if (std::abs(mcCluster.pdgCode()) == PDG_t::kElectron) {
    info.isLepton = true;
    info.convMotherId = getMotherIndexFromChain(mcCluster, mcClusterLooper, PDG_t::kGamma);
    info.isFromConv = info.convMotherId >= 0;

    if (mcCluster.mothersIds().size() > 0 && info.isFromConv) {
      for (size_t i = 1; i < g.emmcparticleIds().size(); ++i) {
        mcClusterLooper.setCursor(g.emmcparticleIds()[i]);
        if (std::abs(mcClusterLooper.pdgCode()) == PDG_t::kElectron && mcClusterLooper.pdgCode() == -1 * mcCluster.pdgCode()) {
          int32_t otherConvMotherId = getMotherIndexFromChain(mcClusterLooper, mcClusterLooper2, PDG_t::kGamma);
          if (otherConvMotherId == info.convMotherId) {
            if (g.amplitude()[i] >= leadingAmplitude - MaxAmpDiff && g.amplitude()[i] > MinAmpThreshold) {
              info.isMergedConv = true;
            }
            break;
          }
        }
      }
    }
  }
  if (std::abs(mcCluster.pdgCode()) == PDG_t::kGamma) {
    info.isPhoton = true;
  }

  info.photonId = o2::aod::pwgem::photonmeson::utils::mcutil::FindMotherInChain(mcCluster, mcParticles, std::vector<int>{PDG_t::kPi0, Pdg::kEta, Pdg::kOmega, Pdg::kEtaPrime});
  info.isFromPi0 = info.photonId >= 0;

  return info;
}

struct EmcalPhotonMcTask {
  static constexpr float EMCALRadius = 440.f;
  static constexpr float PhiVUndefined = -999.f;
  static constexpr float Epsilon = 1.e-6f;

  static constexpr std::array<const char*, static_cast<size_t>(TruthClass::NClasses)> kTruthClassNames = {
    "Conversion", "PhotonPairSamePi0", "PhotonPairDiffPi0", "PhotonPairOnePi0",
    "PhotonElectronSamePi0", "PhotonElectronDiffPi0", "PhotonElectronOnePi0", "BSPhotonElectron",
    "ElectronPairSamePi0", "ElectronPairDiffPi0", "ElectronPairOnePi0",
    "SplitPhotonCluster", "SplitLeptonCluster", "Background"};

  Produces<aod::ConvTagCandidates_001> convTagCandidates;

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<bool> writeTable{"writeTable", true, "write table for ML."};
  Configurable<std::vector<int>> classPrescale{"classPrescale", {1, 1, 700, 25, 1, 350, 15, 1, 1, 35, 2, 1, 1, 1, 1, 1, 1000}, "prescale factor per TruthClass, indexed 0..10 matching the enum order"};
  Configurable<uint32_t> bkgPrescaleSeed{"bkgPrescaleSeed", 42, "seed for the background-prescale RNG"};

  // configurable axis
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {400, 0.0, 0.8}, "invariant mass axis for the neutral meson"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {400, 0., 20.}, "pT axis for the neutral meson"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {20, 0., 100.}, "centrality axis for the current event"};
  ConfigurableAxis thnConfigAxisMult{"thnConfigAxisMult", {60, 0., 60000.}, "multiplicity axis for the current event"};
  ConfigurableAxis thnConfigAxisDeltaEta{"thnConfigAxisDeltaEta", {100, -1, 1}, "delta eta axis"};
  ConfigurableAxis thnConfigAxisDeltaPhi{"thnConfigAxisDeltaPhi", {100, -1, 1}, "delta phi axis"};
  Configurable<bool> useCent{"useCent", false, "flag to enable usage of centrality instead of multiplicity as axis."};

  struct : ConfigurableGroup {
    std::string prefix = "conversiontagging";
    Configurable<float> maxMinv{"maxMinv", 0.1f, "maximum invariant mass for tagging conversions."};
    Configurable<float> maxDeltaEta{"maxDeltaEta", 0.05f, "maximum delta eta between two clusters for tagging them as conversions."};
    Configurable<float> maxDeltaPhi{"maxDeltaPhi", 0.1f, "maximum delta phi between two clusters for tagging them as conversions."};
    Configurable<float> minrConv{"minrConv", 370.f, "minimum conversion Radius of two clusters for tagging them as conversions."};
    Configurable<float> maxrConv{"maxrConv", 430.f, "maximum conversion Radius of two clusters for tagging them as conversions."};
  } conversiontaggingcuts;

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

  struct : ConfigurableGroup {
    std::string prefix = "mesonConfig";
    Configurable<bool> cfgEnableQA{"cfgEnableQA", false, "flag to turn QA plots on/off"};
  } mesonConfig;

  SliceCache cache;

  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::MinClusters, aod::EMEMCClusterMCLabels_001>;

  using Colls = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMMCEventLabels>;

  using McColls = o2::soa::Join<o2::aod::EMMCEvents, o2::aod::BinnedGenPts>;
  using McParticles = EMMCParticles;

  PresliceOptional<aod::EMCEMEventIds> perCollisionEMC = o2::aod::emccluster::pmeventId;
  PresliceOptional<MinMTracks> perEMCClusterMT = o2::aod::mintm::minClusterId;
  PresliceOptional<MinMSTracks> perEMCClusterMS = o2::aod::mintm::minClusterId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int8_t bTruthLabel{};

  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  int mRunNumber{0};
  float dBz{0.f};

  std::mt19937 mRandGen;
  std::uniform_int_distribution<int> mPrescaleDist;

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
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&registry);

    const AxisSpec thnAxisERec{thnConfigAxisPt, "#it{E}_{Rec} (GeV)"};
    const AxisSpec thnAxisPtRec{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};

    const AxisSpec thnAxisrConvRec{100, 0, 500, "#it{R}_{rec}"};
    const AxisSpec thnAxisrConvGen{100, 0, 500, "#it{R}_{gen}"};

    const AxisSpec thnAxisDeltaEta{thnConfigAxisDeltaEta, "#Delta#it{eta}"};
    const AxisSpec thnAxisDeltaPhi{thnConfigAxisDeltaPhi, "#Delta#it{#varphi} (rad)"};

    const AxisSpec thnAxisTagging{static_cast<int>(TagDecision::NTags), -0.5, static_cast<double>(TagDecision::NTags) - 0.5, ""};
    const AxisSpec thnAxisClasses{static_cast<int>(TruthClass::NClasses), -0.5, static_cast<double>(TruthClass::NClasses) - 0.5, ""};

    AxisSpec thnAxisCentOrMult{1, 0., 1., "Centrality/Multiplicity"}; // placeholder, overwritten in init
    if (useCent.value) {
      // PbPb: use centrality
      thnAxisCentOrMult = {thnConfigAxisCent, "Centrality (%)"};
    } else {
      // pp: use multiplicity
      thnAxisCentOrMult = {thnConfigAxisMult, "FT0C Multiplicity"};
    }

    auto hTruthLabel = registry.add<TH2>("hTruthLabel", "Truth label distribution;;Counts", HistType::kTH2D, {{static_cast<int>(TruthClass::NClasses), -0.5, static_cast<double>(TruthClass::NClasses) - 0.5}, thnAxisPtRec});

    // set bin labels once at init, so histogram is human-readable without decoding the enum
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::Conversion) + 1, "Conversion");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::PhotonPairSamePi0) + 1, "PhotonPairSamePi0");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::PhotonPairDiffPi0) + 1, "PhotonPairDiffPi0");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::PhotonPairOnePi0) + 1, "PhotonPairOnePi0");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::PhotonElectronSamePi0) + 1, "PhotonElectronSamePi0");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::PhotonElectronDiffPi0) + 1, "PhotonElectronDiffPi0");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::PhotonElectronOnePi0) + 1, "PhotonElectronOnePi0");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::BSPhotonElectron) + 1, "BSPhotonElectron");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::ElectronPairSamePi0) + 1, "ElectronPairSamePi0");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::ElectronPairDiffPi0) + 1, "ElectronPairDiffPi0");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::ElectronPairOnePi0) + 1, "ElectronPairOnePi0");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::SplitPhotonCluster) + 1, "SplitPhotonCluster");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::SplitLeptonCluster) + 1, "SplitLeptonCluster");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::PhotonBSPhotonPair) + 1, "PhotonBSPhotonPair");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::ElectronBSPhotonPair) + 1, "ElectronBSPhotonPair");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::BSPhotonPair) + 1, "BSPhotonPair");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::Background) + 1, "Background");

    auto hPi0BothResolvedLost = registry.add<TH1>("EMCal/hPi0BothResolvedLost", "Confusion matrix for conversion tagging", HistType::kTH1D, {{2, -0.5, 1.5}});
    hPi0BothResolvedLost->GetXaxis()->SetBinLabel(1, "both resolved (denominator)");
    hPi0BothResolvedLost->GetXaxis()->SetBinLabel(2, "lost to tagging");

    auto hConfusionMatrixConversionTagging = registry.add<TH2>("EMCal/ConfusionMatrixConversionTagging", "Confusion matrix for conversion tagging", HistType::kTH2D, {thnAxisTagging, thnAxisClasses});
    hConfusionMatrixConversionTagging->GetXaxis()->SetBinLabel(1, "not tagged");
    hConfusionMatrixConversionTagging->GetXaxis()->SetBinLabel(2, "tagged");

    hConfusionMatrixConversionTagging->GetYaxis()->SetBinLabel(1, "true conv.");
    hConfusionMatrixConversionTagging->GetYaxis()->SetBinLabel(2, "true #pi^{0}");
    hConfusionMatrixConversionTagging->GetYaxis()->SetBinLabel(3, "true #pi^{0} conv.");
    hConfusionMatrixConversionTagging->GetYaxis()->SetBinLabel(4, "background");
    hConfusionMatrixConversionTagging->GetYaxis()->SetBinLabel(5, "#gamma");

    mRandGen.seed(bkgPrescaleSeed.value);

    if (classPrescale.value.size() != kTruthClassNames.size()) {
      LOG(fatal) << "classPrescale has " << classPrescale.value.size() << " entries, "
                 << "but TruthClass has " << kTruthClassNames.size() << " members -- update classPrescale!";
    }

    LOG(info) << "=== classPrescale configuration ===";
    for (size_t i = 0; i < kTruthClassNames.size(); ++i) {
      LOG(info) << "  [" << i << "] " << kTruthClassNames[i] << " -> prescale = " << classPrescale.value[i];
    }
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
    if (!skipGRPOquery) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3GrpTimestamp);
    }
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
    mRunNumber = collision.runNumber();
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
  /// \param fillHisto flag to enable filling of histograms
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

  // EMCal same event
  void processEmcal(Colls const& collisions, McColls const&, EMCalPhotons const& clusters, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds, EMMCParticles const& mcParticles)
  {
    if (clusters.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    if (collisions.size() <= 0) {
      return;
    }
    std::vector<bool> wasMeassured(mcParticles.size(), false);
    EMBitFlags emcFlagsFromTrueMeson(clusters.size());
    EMBitFlags emcFlagsFromTrueMesonSameGamma(clusters.size());
    EMBitFlags emcFlagsFromTrueConversion(clusters.size());
    EMBitFlags emcFlagsTagging(clusters.size());
    EMBitFlags emcFlags(clusters.size());
    if (clusters.size() > 0) {
      fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds, &registry);
    }
    auto mcCluster1 = mcParticles.begin();       // iterator of MC particle of highest contributor for first cluster
    auto mcCluster2 = mcParticles.begin();       // iterator of MC particle of highest contributor for second cluster
    auto mcClusterLooper = mcParticles.begin();  // iterator of MC particle of other contributor for both clusters
    auto mcClusterLooper2 = mcParticles.begin(); // iterator of MC particle of other contributor for both clusters
    auto mcPhoton1 = mcParticles.begin();        // iterator of MC particle of highest contributor for first cluster that is photon
    auto mcPhoton2 = mcParticles.begin();        // iterator of MC particle of highest contributor for second cluster that is photon
    auto mcMother = mcParticles.begin();         // iterator of MC particle of combined mother of both clusters

    for (const auto& collision : collisions) {
      initCCDB(collision);
      isFullEventSelected(collision, true);

      float centOrMult = getCentralityOrMultiplicity(collision);

      auto photonsEMCPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());

      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photonsEMCPerCollision, photonsEMCPerCollision))) {
        if (!(emcFlags.test(g1.globalIndex())) || !(emcFlags.test(g2.globalIndex()))) {
          continue;
        }

        if (g1.emmcparticleIds().empty() || g2.emmcparticleIds().empty()) {
          // there is a cluster which is just noise, skip
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;

        // z-axis unit vector (beam/field direction)
        ROOT::Math::XYZVector p1 = v1.Vect();
        ROOT::Math::XYZVector p2 = v2.Vect();
        ROOT::Math::XYZVector pSum = p1 + p2;
        ROOT::Math::XYZVector zAxis(0, 0, 1);

        ROOT::Math::XYZVector u = pSum.Unit();
        ROOT::Math::XYZVector nDecay = p1.Cross(p2); // normal to the plane containing p1, p2
        ROOT::Math::XYZVector nRef = u.Cross(zAxis); // normal to the plane containing u and z

        float phiV = PhiVUndefined; // sentinel for degenerate geometry
        if (nDecay.R() > Epsilon && nRef.R() > Epsilon) {
          float cosPhiV = static_cast<float>(nDecay.Unit().Dot(nRef.Unit()));
          cosPhiV = std::clamp(cosPhiV, -1.f, 1.f);
          phiV = std::acos(cosPhiV);
        }

        const float deltaPhi = RecoDecay::constrainAngle(v1.Phi() - v2.Phi(), -o2::constants::math::PI);
        const float deltaEta = v1.Eta() - v2.Eta();
        const float eT1 = v1.Et();
        const float eT2 = v2.Et();
        const float harmonicET = (eT1 * eT2) / (eT1 + eT2);

        // calculate the conversion radius in cm. The formula expects radii in m, B field in T and energies in GeV
        const float rConv = 100.f * (EMCALRadius / 100.f - std::fabs(deltaPhi) / (0.15f * (std::fabs(dBz) / 10.f)) * harmonicET);

        // tag conversion pairs
        if (conversiontaggingcuts.maxMinv >= vMeson.M() && conversiontaggingcuts.maxDeltaEta >= std::fabs(deltaEta) && conversiontaggingcuts.minrConv <= rConv && rConv <= conversiontaggingcuts.maxrConv && conversiontaggingcuts.maxDeltaPhi >= std::fabs(deltaPhi)) {
          emcFlagsTagging.set(g1.globalIndex());
          emcFlagsTagging.set(g2.globalIndex());
        }

        // set MC particle cursors to the largest cluster contributor
        mcCluster1.setCursor(g1.emmcparticleIds()[0]);
        mcCluster2.setCursor(g2.emmcparticleIds()[0]);

        bool areFromSamePi0 = false;
        bool areConversionLegs = false;
        bool areSplitPhotonCluster = false;
        bool areSplitLeptonCluster = false;
        bool areBSPhotonElectron = false;

        auto c1 = classifyCluster(g1, mcCluster1, mcClusterLooper, mcClusterLooper2, mcParticles);
        auto c2 = classifyCluster(g2, mcCluster2, mcClusterLooper, mcClusterLooper2, mcParticles);

        // split-cluster check MUST run first and take priority over everything else --
        // if both clusters share the same dominant MC particle, this is one physical
        // shower reconstructed as two clusters, not a genuine pair of anything.
        const bool isSameDominantParticle = (g1.emmcparticleIds()[0] == g2.emmcparticleIds()[0]);
        if (isSameDominantParticle) {
          if (c1.isPhoton) {
            areSplitPhotonCluster = true;
          } else if (c1.isLepton) {
            areSplitLeptonCluster = true;
          }
        }

        const bool isAnyBSPhoton = c1.isFromBremsstrahlung || c2.isFromBremsstrahlung;
        const bool areBSPhotons = c1.isFromBremsstrahlung && c2.isFromBremsstrahlung;

        // if they are not a split cluster check for proper conversion pair
        if (!isSameDominantParticle && c1.isFromConv && c2.isFromConv && c1.convMotherId == c2.convMotherId) {
          emcFlagsFromTrueConversion.set(g1.globalIndex());
          emcFlagsFromTrueConversion.set(g2.globalIndex());
          areConversionLegs = true;
        }

        // if they are not a split cluster check for neutral meson connection
        if (!isSameDominantParticle && c1.isFromPi0 && c2.isFromPi0) {
          mcPhoton1.setCursor(c1.photonId);
          mcPhoton2.setCursor(c2.photonId);
          mcMother.setCursor(mcPhoton1.mothersIds()[0]);
          if (mcMother.producedByGenerator()) {
            if (c1.photonId == c2.photonId) {
              // bremsstrahlung: one side is a photon born from the other side's lepton lineage
              const bool photonIsBS = (c1.isPhoton && c1.isFromBremsstrahlung) || (c2.isPhoton && c2.isFromBremsstrahlung);
              if (photonIsBS && ((c1.isLepton && c2.isPhoton) || (c2.isLepton && c1.isPhoton))) {
                areBSPhotonElectron = true;
              } else {
                areFromSamePi0 = true;
                emcFlagsFromTrueMesonSameGamma.set(g1.globalIndex());
                emcFlagsFromTrueMesonSameGamma.set(g2.globalIndex());
              }
            } else if (mcPhoton1.mothersIds()[0] == mcPhoton2.mothersIds()[0]) {
              areFromSamePi0 = true;
              emcFlagsFromTrueMeson.set(g1.globalIndex());
              emcFlagsFromTrueMeson.set(g2.globalIndex());
            }
          }
        }

        bTruthLabel = static_cast<int8_t>(TruthClass::Background);
        if (areSplitPhotonCluster) {
          bTruthLabel = static_cast<int8_t>(TruthClass::SplitPhotonCluster);
        } else if (areSplitLeptonCluster) {
          bTruthLabel = static_cast<int8_t>(TruthClass::SplitLeptonCluster);
        } else if (areConversionLegs) {
          bTruthLabel = static_cast<int8_t>(TruthClass::Conversion);
        } else if (areBSPhotonElectron) {
          bTruthLabel = static_cast<int8_t>(TruthClass::BSPhotonElectron);
        } else if (areBSPhotons && (c1.isFromPi0 || c2.isFromPi0)) {
          bTruthLabel = static_cast<int8_t>(TruthClass::BSPhotonPair);
        } else if (isAnyBSPhoton && (c1.isFromPi0 || c2.isFromPi0) && ((c1.isPhoton && c2.isLepton) || (c2.isPhoton && c1.isLepton))) {
          bTruthLabel = static_cast<int8_t>(TruthClass::ElectronBSPhotonPair);
        } else if (isAnyBSPhoton && (c1.isFromPi0 || c2.isFromPi0) && (c1.isPhoton && c2.isPhoton)) {
          bTruthLabel = static_cast<int8_t>(TruthClass::PhotonBSPhotonPair);
        } else if (areFromSamePi0) {
          if ((c1.isLepton && c2.isPhoton) || (c2.isLepton && c1.isPhoton)) {
            bTruthLabel = static_cast<int8_t>(TruthClass::PhotonElectronSamePi0);
          } else if (c1.isPhoton && c2.isPhoton) {
            bTruthLabel = static_cast<int8_t>(TruthClass::PhotonPairSamePi0);
          } else if (c1.isLepton && c2.isLepton) {
            bTruthLabel = static_cast<int8_t>(TruthClass::ElectronPairSamePi0);
          }
        } else if (c1.isFromPi0 && c2.isFromPi0) {
          if ((c1.isLepton && c2.isPhoton) || (c2.isLepton && c1.isPhoton)) {
            bTruthLabel = static_cast<int8_t>(TruthClass::PhotonElectronDiffPi0);
          } else if (c1.isPhoton && c2.isPhoton) {
            bTruthLabel = static_cast<int8_t>(TruthClass::PhotonPairDiffPi0);
          } else if (c1.isLepton && c2.isLepton) {
            bTruthLabel = static_cast<int8_t>(TruthClass::ElectronPairDiffPi0);
          }
        } else if ((c1.isFromPi0 && !c2.isFromPi0) || (!c1.isFromPi0 && c2.isFromPi0)) {
          if ((c1.isLepton && c2.isPhoton) || (c2.isLepton && c1.isPhoton)) {
            bTruthLabel = static_cast<int8_t>(TruthClass::PhotonElectronOnePi0);
          } else if (c1.isPhoton && c2.isPhoton) {
            bTruthLabel = static_cast<int8_t>(TruthClass::PhotonPairOnePi0);
          } else if (c1.isLepton && c2.isLepton) {
            bTruthLabel = static_cast<int8_t>(TruthClass::ElectronPairOnePi0);
          }
        }

        registry.fill(HIST("hTruthLabel"), bTruthLabel, vMeson.Pt());

        // final tree values plus filling
        const int prescale = classPrescale.value[static_cast<uint>(bTruthLabel)];
        const bool keepThisRow = (prescale <= 1) || (std::uniform_int_distribution<int>(0, prescale - 1)(mRandGen) == 0);
        if (writeTable.value && keepThisRow) {
          convTagCandidates(collision.globalIndex(), vMeson.M(), harmonicET, deltaEta, deltaPhi, phiV, g1.e(), g2.e(), g1.m02(), g2.m02(), g1.time(), g2.time(), g1.nCells(), g2.nCells(), c1.purity, c2.purity, bTruthLabel, centOrMult);
        }
      } // pair loop
    } // collision loop

    std::vector<bool> photonSeen(mcParticles.size(), false);   // this decay photon has >=1 resolved cluster
    std::vector<bool> photonTagged(mcParticles.size(), false); // >=1 of those clusters got conversion-tagged
    auto collision = collisions.begin();
    for (const auto& cluster : clusters) {
      if (!(emcFlags.test(cluster.globalIndex()))) {
        continue;
      }
      if (cluster.emmcparticleIds().empty()) {
        // there is a cluster which is just noise, skip
        continue;
      }
      if (cluster.pmeventId() > collision.globalIndex()) {
        collision.setCursor(cluster.pmeventId());
      }

      mcCluster1.setCursor(cluster.emmcparticleIds()[0]);
      int photonid1 = o2::aod::pwgem::photonmeson::utils::mcutil::FindMotherInChain(mcCluster1, mcParticles, std::vector<int>{PDG_t::kPi0, Pdg::kEta, Pdg::kOmega, Pdg::kEtaPrime});
      int motherId = -1;
      if (photonid1 >= 0) {
        mcPhoton1.setCursor(photonid1);
        motherId = mcPhoton1.mothersIds()[0];

        photonSeen[static_cast<size_t>(photonid1)] = true;
        const bool isTagged = !emcFlagsTagging.test(cluster.globalIndex());
        if (isTagged) {
          photonTagged[static_cast<size_t>(photonid1)] = true;
        }
      }

      const bool alreadyCountedForMeson = (motherId > -1 && wasMeassured[static_cast<size_t>(motherId)]);

      if (mcCluster1.pdgCode() == PDG_t::kGamma) {
        registry.fill(HIST("EMCal/ConfusionMatrixConversionTagging"), emcFlagsTagging.test(cluster.globalIndex()) ? 0 : 1, static_cast<float>(ClusterTruthClass::Photon));
      }
      if (!emcFlagsFromTrueConversion.test(cluster.globalIndex())) {
        registry.fill(HIST("EMCal/ConfusionMatrixConversionTagging"), emcFlagsTagging.test(cluster.globalIndex()) ? 0 : 1, static_cast<float>(ClusterTruthClass::Conversion));
        if (!emcFlagsFromTrueMesonSameGamma.test(cluster.globalIndex()) && !alreadyCountedForMeson) {
          registry.fill(HIST("EMCal/ConfusionMatrixConversionTagging"), emcFlagsTagging.test(cluster.globalIndex()) ? 0 : 1, static_cast<float>(ClusterTruthClass::Pi0Conversion));
        }
      }
      if (!emcFlagsFromTrueMeson.test(cluster.globalIndex()) && !alreadyCountedForMeson) {
        registry.fill(HIST("EMCal/ConfusionMatrixConversionTagging"), emcFlagsTagging.test(cluster.globalIndex()) ? 0 : 1, static_cast<float>(ClusterTruthClass::Pi0));
      }
      if (emcFlagsFromTrueMeson.test(cluster.globalIndex()) && emcFlagsFromTrueConversion.test(cluster.globalIndex())) {
        registry.fill(HIST("EMCal/ConfusionMatrixConversionTagging"), emcFlagsTagging.test(cluster.globalIndex()) ? 0 : 1, static_cast<float>(ClusterTruthClass::Background));
      }
      if (motherId > -1 && !wasMeassured[static_cast<size_t>(motherId)]) {
        wasMeassured[static_cast<size_t>(motherId)] = true;
      }
    } // end of loop over cluster

    for (const auto& mcPart : mcParticles) {
      if (mcPart.pdgCode() != PDG_t::kPi0 && mcPart.pdgCode() != Pdg::kEta) {
        continue;
      }
      if (!mcPart.producedByGenerator()) {
        continue;
      }
      if (mcPart.daughtersIds().size() != 2) { // o2-linter: disable=magic-number (self explanatory and does not need extra named variable)
        continue;
      }

      const int d0 = mcPart.daughtersIds()[0];
      const int d1 = mcPart.daughtersIds()[1];
      mcPhoton1.setCursor(d0);
      mcPhoton2.setCursor(d1);
      if (mcPhoton1.pdgCode() != PDG_t::kGamma || mcPhoton2.pdgCode() != PDG_t::kGamma) {
        continue; // skip anything that is not pi0/eta -> γγ.
      }

      const bool bothResolved = photonSeen[static_cast<size_t>(d0)] && photonSeen[static_cast<size_t>(d1)];
      if (!bothResolved) {
        continue; // not in the denominator -- ceiling case, not attributable to the cut
      }

      const bool lost = photonTagged[static_cast<size_t>(d0)] || photonTagged[static_cast<size_t>(d1)];
      registry.fill(HIST("EMCal/hPi0BothResolvedLost"), 0.0); // "denominator" bin, once per bothResolved pi0
      if (lost) {
        registry.fill(HIST("EMCal/hPi0BothResolvedLost"), 1.0); // "lost" bin
      }
    }
  }
  PROCESS_SWITCH(EmcalPhotonMcTask, processEmcal, "Process for pcm and emcal photons", true);

}; // End struct EmcalPhotonMcTask

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<EmcalPhotonMcTask>(context)};
}
