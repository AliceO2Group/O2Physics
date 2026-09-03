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
#include "PWGEM/PhotonMeson/Core/EMCConversionCandidate.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/EmMlResponseEMCConversion.h"
#include "PWGEM/PhotonMeson/DataModel/ConversionMl.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Utils/ParticleOrigin.h"

#include "Common/Core/RecoDecay.h"
#include "Tools/ML/MlResponse.h"

#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
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
#include <TMCProcess.h>
#include <TPDGCode.h>
#include <TTree.h>

#include <sys/types.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <unordered_map>
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
using namespace o2::analysis::em;

constexpr float MinAmpThreshold = 0.2f; // Minimum cluster amplitude threshold to count as significant
constexpr float MinAmpFraction = 0.6f;  // Minimum fraction of particle energy that it needs to deposit in cluster to count as significant
enum CentralityEstimator {
  None = 0,
  CFT0A,
  CFT0C,
  CFT0M,
  NCentralityEstimators
};

enum class TruthClass {
  Conversion,                // two electron legs, same conversion vertex
  GammaGammaSamePi0,         // both clusters are DIRECT daughter photons of the same generator-level meson
  GammaGammaAnnihilation,    // two photons from the same e+e- annihilation vertex
  BSPhotonElectron,          // a bremsstrahlung photon paired with the specific lepton that radiated it
  PhotonComptonElectronPair, // a Compton-scattered photon paired with its own recoil electron
  ElectronPairSamePi0,       // two lepton clusters, directly from same meson
  CrossConvertedSiblings,    // two lepton clusters, same meson cross-converted siblings
  DalitzDecaySiblings,       // two lepton clusters, same meson from Dalitz
  SplitPhotonCluster,        // one physical photon shower reconstructed as two clusters
  SplitLeptonCluster,        // one physical lepton shower reconstructed as two clusters
  IndirGammaGammaSamePi0,    // both clusters are indirect daughter photons of the same generator-level meson
  SameMesonIndirect,         // both clusters trace to the SAME generator-level meson, but at least one path passes through extra generations (e.g. a further conversion/BS/scatter) before reaching the cluster -- NOT a clean two-body relationship
  Background,                // no common meson ancestor at all -- genuinely uncorrelated combinatorics
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
  bool isLepton = false;                           // is cluster from a lepton
  bool isPhoton = false;                           // is cluster from a photon
  bool isMergedConv = false;                       // is cluster from a merged conversion
  LeptonOrigin leptonOrigin = LeptonOrigin::Other; // origin of lepton in case isLepton == true
  PhotonOrigin photonOrigin = PhotonOrigin::Other; // origin of photon in case isPhoton == true
  int photonMotherId = -1;                         // for leptons: the photon this lepton traces to
  int mesonId = -1;                                // for Decay photons: the pi0/eta/omega/etaprime id | for leptons from decay photons: the pi0/eta/omega/etaprime id
  int hardPartonId = -1;                           // for Direct photons: the quark/gluon id
  int decayPhotonId = -1;                          // unified: the photon (self, if isPhoton; or photonMotherId, if isLepton) whose immediate mother is a meson. -1 if not applicable.
  float purity = 0.f;                              // fraction of energy the main particle contributed to the cluster
  float radius = 0.f;                              // radius in xy from where the main contributor to the cluster originated
};

template <o2::soa::is_iterator TCluster, o2::soa::is_iterator TIter, o2::soa::is_table McParticles>
ClusterMcInfo classifyCluster(const TCluster& g, TIter& mcCluster, TIter& mcClusterLooper, TIter& mcClusterLooper2, McParticles const& mcParticles)
{
  static const std::array<int, 4> kMesonPdgs{PDG_t::kPi0, Pdg::kEta, Pdg::kOmega, Pdg::kEtaPrime};
  ClusterMcInfo info;
  mcCluster.setCursor(g.emmcparticleIds()[0]);
  info.radius = std::hypot(mcCluster.vx(), mcCluster.vy());
  float leadingAmplitude = g.amplitude()[0];
  info.purity = leadingAmplitude;
  info.mesonId = o2::aod::pwgem::photonmeson::utils::mcutil::GetMesonInChain(mcCluster, mcParticles, kMesonPdgs);
  if (std::abs(mcCluster.pdgCode()) == PDG_t::kElectron) {
    info.isLepton = true;
    info.leptonOrigin = getLeptonOriginType(mcCluster, mcClusterLooper2, kMesonPdgs);
    mcClusterLooper.setCursor(g.emmcparticleIds()[0]);
    if (info.leptonOrigin == LeptonOrigin::Conversion) {
      if (!mcCluster.has_mothers()) [[unlikely]] {
        // conersion with no mother does not make any sense
        info.leptonOrigin = LeptonOrigin::Other;
      } else {
        info.photonMotherId = mcCluster.mothersIds()[0];
        // since this is a real conversion where both daughters exist, check if the other conversion leg entered this cluster as well
        for (size_t i = 1; i < g.emmcparticleIds().size(); ++i) {
          mcClusterLooper.setCursor(g.emmcparticleIds()[i]);
          if (std::abs(mcClusterLooper.pdgCode()) == PDG_t::kElectron && mcClusterLooper.pdgCode() == -1 * mcCluster.pdgCode()) {
            mcClusterLooper2.setCursor(mcClusterLooper.globalIndex());
            int32_t otherConvMotherId = getMotherIndexFromChain(mcClusterLooper2, PDG_t::kGamma);
            if (otherConvMotherId == info.photonMotherId) {
              float energyFraction = (g.amplitude()[i] * g.e()) / mcClusterLooper.e();
              if (energyFraction >= MinAmpFraction && g.amplitude()[i] > MinAmpThreshold) {
                info.isMergedConv = true;
              }
              // we found the sibling leg no need to search further
              break;
            } // if (otherConvMotherId == info.photonMotherId)
          } // if (std::abs(mcClusterLooper.pdgCode()) == PDG_t::kElectron && mcClusterLooper.pdgCode() == -1 * mcCluster.pdgCode())
        } // end of loop over other cluster contributions
      } // particle has mothers
    } // if(info.leptonOrigin == LeptonOrigin::Conversion)
    if (info.photonMotherId >= 0) {
      mcClusterLooper.setCursor(info.photonMotherId);
      info.photonOrigin = getPhotonOriginType(mcClusterLooper, mcClusterLooper2, kMesonPdgs, info.hardPartonId);
      if (info.photonOrigin == PhotonOrigin::Decay) {
        info.decayPhotonId = info.photonMotherId;
        info.mesonId = o2::aod::pwgem::photonmeson::utils::mcutil::GetMesonInChain(mcClusterLooper, mcParticles, kMesonPdgs);
      }
    }
  }
  if (std::abs(mcCluster.pdgCode()) == PDG_t::kGamma) {
    info.isPhoton = true;
    info.photonOrigin = getPhotonOriginType(mcCluster, mcClusterLooper, kMesonPdgs, info.hardPartonId);
    if (info.photonOrigin == PhotonOrigin::Decay) {
      info.mesonId = o2::aod::pwgem::photonmeson::utils::mcutil::GetMesonInChain(mcCluster, mcParticles, kMesonPdgs);
    }
  }
  return info;
}

struct EmcalPhotonMcTask {
  static constexpr float EMCALRadius = 440.f;
  static constexpr float PhiVUndefined = -999.f;
  static constexpr float Epsilon = 1.e-6f;

  static constexpr std::array<std::array<double, 2>, 1> defaultCutsMl{{{0.0, 0.25}}};

  static constexpr std::array<const char*, static_cast<size_t>(TruthClass::NClasses)> kTruthClassNames = {
    "Conversion", "GammaGammaSamePi0", "GammaGammaAnnihilation", "BSPhotonElectron",
    "PhotonComptonElectronPair", "ElectronPairSamePi0", "CrossConvertedSiblings", "DalitzDecaySiblings",
    "SplitPhotonCluster", "SplitLeptonCluster", "IndirGammaGammaSamePi0", "SameMesonIndirect",
    "Background"};

  Produces<aod::ConvTagCandidates_001> convTagCandidates;

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> writeTable{"writeTable", true, "write table for ML."};
  Configurable<std::vector<int>> classPrescale{"classPrescale",
                                               {
                                                 1,    // Conversion
                                                 1,    // GammaGammaSamePi0
                                                 1,    // GammaGammaAnnihilation
                                                 1,    // BSPhotonElectron
                                                 1,    // PhotonComptonElectronPair
                                                 1,    // ElectronPairSamePi0
                                                 1,    // CrossConvertedSiblings
                                                 1,    // DalitzDecaySiblings
                                                 1,    // SplitPhotonCluster
                                                 1,    // SplitLeptonCluster
                                                 1,    // IndirGammaGammaSamePi0
                                                 1,    // SameMesonIndirect
                                                 5000, // Background
                                               },
                                               "prescale factor per TruthClass, indexed 0..12 matching the enum order"};

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

  struct : ConfigurableGroup {
    std::string prefix = "mlConfig";
    Configurable<bool> useMlTagging{"useMlTagging", false, "use ML score instead of box cut for conversion tagging"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "load ML model from CCDB"};
    Configurable<std::vector<std::string>> mlInputFeatures{
      "mlInputFeatures",
      {"minv", "deltaEta", "deltaR", "phiv", "rConv", "totE", "e2", "e1", "deltaPhi"},
      "input feature names -- content and order must match the Python training FEATURES list"};
    Configurable<std::string> mlModelPathLocal{"mlModelPathLocal", "/data/mhemmer/O2ML/code/conversion_tagging_bdt_conversion_splits_brems.onnx", "local ONNX model path"};
    Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"Users/m/mhemmer/EM/ML/"}, "Paths of models on CCDB"};
    Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"conversion_tagging_bdt_conversion_splits_brems.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
    Configurable<float> mlThreshold{"mlThreshold", 0.5f, "positive-class score threshold for tagging"};
    Configurable<LabeledArray<double>> cutsMl{"cutsMl", {defaultCutsMl[0].data(), 1, 2, {"pT bin 0"}, {
                                                                                                        "score photon pairs",
                                                                                                        "score conversion pairs",
                                                                                                      }},
                                              "ML selections per pT bin"};
  } mlConfig;

  SliceCache cache;

  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::MinClusters, aod::EMEMCClusterMCLabels_001>;

  using Colls = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMMCEventLabels, aod::EmMagFields>;

  using McColls = o2::soa::Join<o2::aod::EMMCEvents, o2::aod::BinnedGenPts>;
  using McParticles = EMMCParticles;

  PresliceOptional<aod::EMCEMEventIds> perCollisionEMC = o2::aod::emccluster::pmeventId;
  PresliceOptional<MinMTracks> perEMCClusterMT = o2::aod::mintm::minClusterId;
  PresliceOptional<MinMSTracks> perEMCClusterMS = o2::aod::mintm::minClusterId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int8_t bTruthLabel{};

  o2::ccdb::CcdbApi ccdbApi;
  int mRunNumber{0};
  float dBz{0.f};

  std::mt19937 mRandGen;
  std::uniform_int_distribution<int> mPrescaleDist;

  o2::analysis::em::emcconv::EmMlResponseEMCConversion<float> mMlResponse;

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

    defineEMEventCut();
    defineEMCCut();
    fEMCCut.addQAHistograms(&registry);
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&registry);

    const AxisSpec thnAxisERec{thnConfigAxisPt, "#it{E}_{Rec} (GeV)"};
    const AxisSpec thnAxisPtRec{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};

    const AxisSpec thnAxisrConvRec{1000, 0, 500, "#it{R}_{rec}"};
    const AxisSpec thnAxisrConvGen{1000, 0, 500, "#it{R}_{gen}"};

    const AxisSpec thnAxisDeltaEta{thnConfigAxisDeltaEta, "#Delta#it{eta}"};
    const AxisSpec thnAxisDeltaPhi{thnConfigAxisDeltaPhi, "#Delta#it{#varphi} (rad)"};

    const AxisSpec thnAxisTagging{static_cast<int>(TagDecision::NTags), -0.5, static_cast<double>(TagDecision::NTags) - 0.5, ""};
    const AxisSpec thnAxisClasses{static_cast<int>(TruthClass::NClasses), -0.5, static_cast<double>(TruthClass::NClasses) - 0.5, ""};

    const AxisSpec thnAxisM02{250, 0., 2.5, "#it{M}_{02}"};

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
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::GammaGammaSamePi0) + 1, "GammaGammaSamePi0");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::GammaGammaAnnihilation) + 1, "GammaGammaAnnihilation");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::BSPhotonElectron) + 1, "BSPhotonElectron");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::PhotonComptonElectronPair) + 1, "PhotonComptonElectronPair");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::ElectronPairSamePi0) + 1, "ElectronPairSamePi0");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::CrossConvertedSiblings) + 1, "CrossConvertedSiblings");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::DalitzDecaySiblings) + 1, "DalitzDecaySiblings");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::SplitPhotonCluster) + 1, "SplitPhotonCluster");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::SplitLeptonCluster) + 1, "SplitLeptonCluster");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::IndirGammaGammaSamePi0) + 1, "IndirGammaGammaSamePi0");
    hTruthLabel->GetXaxis()->SetBinLabel(static_cast<int>(TruthClass::SameMesonIndirect) + 1, "SameMesonIndirect");
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

    registry.add<TH1>("hBSRadius", "Radius of BS photons;;counts", HistType::kTH1D, {thnAxisrConvGen});
    registry.add<TH1>("hLeptonRadius", "Radius of leptons;;counts", HistType::kTH1D, {thnAxisrConvGen});
    registry.add<TH1>("hConvLeptonRadius", "Radius of leptons from conversions;;counts", HistType::kTH1D, {thnAxisrConvGen});

    registry.add<TH1>("Photon/M02", "M02 distribution;;counts", HistType::kTH1D, {thnAxisM02});
    auto hPhotonProcess = registry.add<TH1>("Photon/hProcess", "Production process type", HistType::kTH1D, {{kMaxMCProcess, -0.5, kMaxMCProcess - 0.5}});
    for (int i = 0; i < kMaxMCProcess; ++i) {
      hPhotonProcess->GetXaxis()->SetBinLabel(i + 1, TMCProcessName[i]);
    }
    registry.addClone("Photon/", "Electron/");
    registry.addClone("Photon/", "Positron/");
    registry.addClone("Photon/", "BSPhoton/");
    registry.addClone("Photon/", "MergedConv/");
    registry.addClone("Photon/", "ConvElectron/");
    registry.addClone("Photon/", "ConvPositron/");
    registry.addClone("Photon/", "Other/");
    registry.addClone("Photon/", "Lepton/");

    auto hClusterType = registry.add<TH2>("hClusterType", "Truth label distribution;;Counts", HistType::kTH2D, {{8, -0.5, 7.5}, thnAxisPtRec});
    hClusterType->GetXaxis()->SetBinLabel(1, "Photon");
    hClusterType->GetXaxis()->SetBinLabel(2, "Electron");
    hClusterType->GetXaxis()->SetBinLabel(3, "Positron");
    hClusterType->GetXaxis()->SetBinLabel(4, "BSPhoton");
    hClusterType->GetXaxis()->SetBinLabel(5, "MergedConv");
    hClusterType->GetXaxis()->SetBinLabel(6, "Conv electron");
    hClusterType->GetXaxis()->SetBinLabel(7, "Conv positron");
    hClusterType->GetXaxis()->SetBinLabel(8, "Other");

    auto hBSLeptonFate = registry.add<TH1>("hBSLeptonFate", "Fate of the Bremsstrahlungsphotons mother lepton;;Counts", HistType::kTH1D, {{3, -0.5, 2.5}});
    hBSLeptonFate->GetXaxis()->SetBinLabel(1, "dominant");
    hBSLeptonFate->GetXaxis()->SetBinLabel(2, "non dominant");
    hBSLeptonFate->GetXaxis()->SetBinLabel(3, "absent");

    if (mlConfig.useMlTagging.value) {
      registry.add<TH1>("hMlScore", "BDT score;;Counts", HistType::kTH1D, {{100, -10, 10}});
    }

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
    mRunNumber = collision.runNumber();

    auto timestamp = collision.timestamp();
    // Fetch magnetic field from ccdb for current collision
    dBz = collision.grpMagField().getNominalL3Field();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << dBz << " kZG";

    if (mlConfig.useMlTagging.value) {
      // single bin, full range -- see earlier discussion: pT/energy-binned
      // thresholds are a straightforward future extension of this same
      // machinery if ever needed, not used right now
      std::vector<double> binsLimits = {0., 1000.};
      std::vector<int> cutDir = {
        static_cast<int>(o2::cuts_ml::CutDirection::CutNot),     // class 0 (negative-class prob) -- no cut
        static_cast<int>(o2::cuts_ml::CutDirection::CutSmaller), // class 1 (positive-class prob) -- reject if score < threshold
      };

      mMlResponse.configure(binsLimits, mlConfig.cutsMl, cutDir, /*nClasses=*/2);
      mMlResponse.cacheInputFeaturesIndices(mlConfig.mlInputFeatures.value);
      if (mlConfig.loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        mMlResponse.setModelPathsCCDB(mlConfig.onnxFileNames, ccdbApi, mlConfig.modelPathsCCDB.value, timestamp);
      } else {
        mMlResponse.setModelPathsLocal({mlConfig.mlModelPathLocal.value});
      }
      mMlResponse.init();

      LOG(info) << "ML conversion tagging enabled -- model: " << mlConfig.mlModelPathLocal.value
                << ", threshold: " << mlConfig.mlThreshold.value;
    }
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
    EMBitFlags emcFlagsMlTagging(clusters.size());
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
      if (!isFullEventSelected(collision, true)) {
        continue;
      }

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

        if (mlConfig.useMlTagging.value) {
          o2::analysis::em::EMCConversionCandidate candidate{
            .mMinv = static_cast<float>(vMeson.M()), .mDeltaEta = deltaEta, .mDeltaR = std::hypot(deltaEta, deltaPhi), .mPhiv = phiV, .mRConv = rConv, .mTotE = (g2.e() + g1.e()), .mE2 = g2.e(), .mE1 = g1.e(), .mDeltaPhi = deltaPhi};
          std::vector<float> mlInput = mMlResponse.getInputFeatures(candidate);
          std::vector<float> mlOutput;
          bool isTagged = mMlResponse.isSelectedMl(mlInput, 0.f, mlOutput);
          if (isTagged) {
            emcFlagsMlTagging.set(g1.globalIndex());
            emcFlagsMlTagging.set(g2.globalIndex());
          }
          registry.fill(HIST("hMlScore"), mlOutput[1]); // positive-class score, always, tagged or not
        }

        // set MC particle cursors to the largest cluster contributor
        mcCluster1.setCursor(g1.emmcparticleIds()[0]);
        mcCluster2.setCursor(g2.emmcparticleIds()[0]);

        auto c1 = classifyCluster(g1, mcCluster1, mcClusterLooper, mcClusterLooper2, mcParticles);
        auto c2 = classifyCluster(g2, mcCluster2, mcClusterLooper, mcClusterLooper2, mcParticles);

        const bool areFromSamePi0 = c1.mesonId >= 0 && c1.mesonId == c2.mesonId;
        bool areConversionLegs = false;
        bool areSplitPhotonCluster = false;
        bool areSplitLeptonCluster = false;

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

        // if they are not a split cluster check for proper conversion pair
        if (!isSameDominantParticle && c1.leptonOrigin == LeptonOrigin::Conversion && c2.leptonOrigin == LeptonOrigin::Conversion && c1.photonMotherId == c2.photonMotherId) {
          emcFlagsFromTrueConversion.set(g1.globalIndex());
          emcFlagsFromTrueConversion.set(g2.globalIndex());
          areConversionLegs = true;
        }

        // if they are not a split cluster check for neutral meson connection
        if (!isSameDominantParticle && c1.decayPhotonId >= 0 && c2.decayPhotonId >= 0) {
          mcPhoton1.setCursor(c1.decayPhotonId);
          mcPhoton2.setCursor(c2.decayPhotonId);
          mcMother.setCursor(mcPhoton1.mothersIds()[0]);
          if (mcMother.producedByGenerator()) {
            if (c1.mesonId == c2.mesonId) {
              // bremsstrahlung: one side is a photon born from the other side's lepton lineage
              const bool photonIsBS = (c1.isPhoton && c1.photonOrigin == PhotonOrigin::Bremsstrahlung) || (c2.isPhoton && c2.photonOrigin == PhotonOrigin::Bremsstrahlung);
              if (photonIsBS && ((c1.isLepton && c2.isPhoton) || (c2.isLepton && c1.isPhoton))) {
                // nothing
              } else {
                emcFlagsFromTrueMesonSameGamma.set(g1.globalIndex());
                emcFlagsFromTrueMesonSameGamma.set(g2.globalIndex());
              }
            } else if (mcPhoton1.mothersIds()[0] == mcPhoton2.mothersIds()[0]) {
              emcFlagsFromTrueMeson.set(g1.globalIndex());
              emcFlagsFromTrueMeson.set(g2.globalIndex());
            }
          }
        }

        bTruthLabel = static_cast<int8_t>(TruthClass::Background);
        if (!mcCluster1.has_mothers() || !mcCluster2.has_mothers()) {
          registry.fill(HIST("hTruthLabel"), bTruthLabel, vMeson.Pt());

          // final tree values plus filling
          const int prescale = classPrescale.value[static_cast<uint>(bTruthLabel)];
          const bool keepThisRow = (prescale <= 1) || (std::uniform_int_distribution<int>(0, prescale - 1)(mRandGen) == 0);
          if (writeTable.value && keepThisRow) {
            convTagCandidates(collision.globalIndex(), vMeson.M(), harmonicET, deltaEta, deltaPhi, phiV, g1.e(), g2.e(), g1.m02(), g2.m02(), g1.time(), g2.time(), g1.nCells(), g2.nCells(), c1.purity, c2.purity, bTruthLabel, centOrMult);
          }
          continue;
        }
        if (areSplitPhotonCluster) {
          bTruthLabel = static_cast<int8_t>(TruthClass::SplitPhotonCluster);
        } else if (areSplitLeptonCluster) {
          bTruthLabel = static_cast<int8_t>(TruthClass::SplitLeptonCluster);
        } else if (areConversionLegs) {
          bTruthLabel = static_cast<int8_t>(TruthClass::Conversion);
        } else if (c1.photonOrigin == PhotonOrigin::Annihilation && c2.photonOrigin == PhotonOrigin::Annihilation && mcCluster1.mothersIds()[0] == mcCluster2.mothersIds()[0]) {
          bTruthLabel = static_cast<int8_t>(TruthClass::GammaGammaAnnihilation);
        } else if ((c1.photonOrigin == PhotonOrigin::Bremsstrahlung && mcCluster1.mothersIds()[0] == mcCluster2.globalIndex()) || (c2.photonOrigin == PhotonOrigin::Bremsstrahlung && mcCluster2.mothersIds()[0] == mcCluster1.globalIndex())) {
          bTruthLabel = static_cast<int8_t>(TruthClass::BSPhotonElectron);
        } else if ((c1.leptonOrigin == LeptonOrigin::Compton && mcCluster1.mothersIds()[0] == mcCluster2.globalIndex()) || (c2.leptonOrigin == LeptonOrigin::Compton && mcCluster2.mothersIds()[0] == mcCluster1.globalIndex())) {
          bTruthLabel = static_cast<int8_t>(TruthClass::PhotonComptonElectronPair);
        } else if (c1.leptonOrigin == LeptonOrigin::DirectMesonDecay && c2.leptonOrigin == LeptonOrigin::DirectMesonDecay && mcCluster1.mothersIds()[0] == mcCluster2.mothersIds()[0]) { // both cluster are leptons that come from the same meson decay
          mcMother.setCursor(mcCluster1.mothersIds()[0]);
          if (mcMother.daughtersIds().size() == 2) {
            bTruthLabel = static_cast<int8_t>(TruthClass::ElectronPairSamePi0);
          } else if (mcMother.daughtersIds().size() == 3) {
            bTruthLabel = static_cast<int8_t>(TruthClass::DalitzDecaySiblings);
          }
        } else if (c1.leptonOrigin == LeptonOrigin::Conversion && c2.leptonOrigin == LeptonOrigin::Conversion && mcCluster1.mothersIds()[0] != mcCluster2.mothersIds()[0] && areFromSamePi0) { // both cluster are leptons that come from different conversions that come from the same meson
          bTruthLabel = static_cast<int8_t>(TruthClass::CrossConvertedSiblings);
        } else if (c1.photonOrigin == PhotonOrigin::Decay && c2.photonOrigin == PhotonOrigin::Decay && areFromSamePi0) { // both clusters are photons from decay from the same meson
          bTruthLabel = static_cast<int8_t>(TruthClass::GammaGammaSamePi0);
        } else if (c1.isPhoton && c2.isPhoton && areFromSamePi0) {
          bTruthLabel = static_cast<int8_t>(TruthClass::IndirGammaGammaSamePi0);
        } else if (areFromSamePi0) { // both cluster do not fit into one of the categories above, but they share a common meson ancestry
          bTruthLabel = static_cast<int8_t>(TruthClass::SameMesonIndirect);
        }
        registry.fill(HIST("hTruthLabel"), bTruthLabel, vMeson.Pt());

        // final tree values plus filling
        const int prescale = classPrescale.value[static_cast<uint>(bTruthLabel)];
        const bool keepThisRow = (prescale <= 1) || (std::uniform_int_distribution<int>(0, prescale - 1)(mRandGen) == 0);
        if (writeTable.value && keepThisRow) {
          convTagCandidates(collision.globalIndex(), vMeson.M(), harmonicET, deltaEta, deltaPhi, phiV, g1.e(), g2.e(), g1.m02(), g2.m02(), g1.time(), g2.time(), g1.nCells(), g2.nCells(), c1.purity, c2.purity, bTruthLabel, centOrMult);
        }
      } // pair loop

      // key: MC particle global index -> list of (cluster global index, contributor rank)
      std::unordered_map<int, std::vector<std::pair<int64_t, size_t>>> particleToClusterContributions;
      for (const auto& cluster : photonsEMCPerCollision) {
        if (!emcFlags.test(cluster.globalIndex())) {
          continue;
        }
        const auto& ids = cluster.emmcparticleIds();
        for (size_t i = 0; i < ids.size(); ++i) {
          particleToClusterContributions[ids[i]].emplace_back(cluster.globalIndex(), i);
        }
      } // cluster loop

      for (const auto& cluster : photonsEMCPerCollision) {
        if (!emcFlags.test(cluster.globalIndex())) {
          continue;
        }
        auto c1 = classifyCluster(cluster, mcCluster1, mcClusterLooper, mcClusterLooper2, mcParticles);

        // NEW: bremsstrahlung sibling-fate check, mirrors the conversion one above
        if (c1.photonOrigin == PhotonOrigin::Bremsstrahlung) {
          // mcCluster1 is currently sitting on the bremsstrahlung photon itself
          // (classifyCluster leaves it there for the photon branch) -- its
          // immediate mother is the radiating lepton we want to look up.
          if (mcCluster1.has_mothers()) {
            const int radiatingLeptonId = mcCluster1.mothersIds()[0];

            auto it = particleToClusterContributions.find(radiatingLeptonId);
            if (it == particleToClusterContributions.end()) {
              registry.fill(HIST("hBSLeptonFate"), 2); // radiating lepton absent from any cluster
            } else {
              bool isDominantSomewhere = false;
              for (const auto& [clusterId, rank] : it->second) {
                if (clusterId == cluster.globalIndex()) {
                  continue; // skip itself (shouldn't normally match, but same safety as before)
                }
                if (rank == 0) {
                  isDominantSomewhere = true;
                  break;
                }
              }
              registry.fill(HIST("hBSLeptonFate"), isDominantSomewhere ? 0 : 1); // 0=dominant elsewhere, 1=leakage-only
            }
          }
        }
      } // cluster loop
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
      if (!isFullEventSelected(collision, false)) {
        continue;
      }

      auto clusterMcInfo = classifyCluster(cluster, mcCluster1, mcClusterLooper, mcClusterLooper2, mcParticles);
      if (clusterMcInfo.photonOrigin == PhotonOrigin::Bremsstrahlung) {
        registry.fill(HIST("hBSRadius"), clusterMcInfo.radius);
        registry.fill(HIST("BSPhoton/M02"), cluster.m02());
        registry.fill(HIST("hClusterType"), 3, cluster.e());
      } else if (clusterMcInfo.isMergedConv) {
        registry.fill(HIST("MergedConv/M02"), cluster.m02());
        registry.fill(HIST("hClusterType"), 4, cluster.e());
      } else if (clusterMcInfo.leptonOrigin == LeptonOrigin::Conversion) {
        registry.fill(HIST("hConvLeptonRadius"), clusterMcInfo.radius);
        if (mcCluster1.pdgCode() == PDG_t::kElectron) {
          registry.fill(HIST("ConvElectron/M02"), cluster.m02());
          registry.fill(HIST("hClusterType"), 5, cluster.e());
        } else if (mcCluster1.pdgCode() == PDG_t::kPositron) {
          registry.fill(HIST("ConvPositron/M02"), cluster.m02());
          registry.fill(HIST("hClusterType"), 6, cluster.e());
        }
      } else if (clusterMcInfo.isLepton) {
        if (mcCluster1.pdgCode() == PDG_t::kElectron) {
          registry.fill(HIST("Electron/M02"), cluster.m02());
          registry.fill(HIST("hClusterType"), 1, cluster.e());
        } else if (mcCluster1.pdgCode() == PDG_t::kPositron) {
          registry.fill(HIST("Positron/M02"), cluster.m02());
          registry.fill(HIST("hClusterType"), 2, cluster.e());
        }
        registry.fill(HIST("hLeptonRadius"), clusterMcInfo.radius);
      } else if (clusterMcInfo.isPhoton) {
        registry.fill(HIST("hClusterType"), 0, cluster.e());
        registry.fill(HIST("Photon/M02"), cluster.m02());
      } else {
        registry.fill(HIST("hClusterType"), 7, cluster.e());
        registry.fill(HIST("Other/M02"), cluster.m02());
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
        registry.fill(HIST("Photon/hProcess"), mcCluster1.getProcess());
      }
      if (std::abs(mcCluster1.pdgCode()) == PDG_t::kElectron) {
        registry.fill(HIST("Lepton/hProcess"), mcCluster1.getProcess());
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
    } // end of loop over mc particles
  }
  PROCESS_SWITCH(EmcalPhotonMcTask, processEmcal, "Process for pcm and emcal photons", true);

}; // End struct EmcalPhotonMcTask

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<EmcalPhotonMcTask>(context)};
}
