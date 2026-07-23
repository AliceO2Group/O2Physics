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
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

#include "Common/Core/RecoDecay.h"

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

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TF1.h>
#include <TH2.h>
#include <TPDGCode.h>
#include <TTree.h>

#include <cmath>
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

enum CentralityEstimator {
  None = 0,
  CFT0A,
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

enum class TruthClass {
  Conversion = 0,
  Pi0,
  Pi0Conversion,
  Background,
  Photon,
  NClasses
};

enum class TagDecision {
  NotTagged = 0,
  Tagged,
  NTags
};

struct EmcalPhotonMcTask {
  static constexpr float EMCALRadius = 440.f;

  o2::framework::Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  o2::framework::Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  o2::framework::Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  o2::framework::Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  o2::framework::Configurable<bool> writeTree{"writeTree", true, "write tree for ML."};

  // configurable axis
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {400, 0.0, 0.8}, "invariant mass axis for the neutral meson"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {400, 0., 20.}, "pT axis for the neutral meson"};
  ConfigurableAxis thnConfigAxisERelative{"thnConfigAxisERelative", {600, -1., 5.}, "(E rec - E true) / E true axis"};
  ConfigurableAxis thnConfigAxisPRelative{"thnConfigAxisPRelative", {300, -1., 2.}, "(P rec - P true) / P true axis"};
  ConfigurableAxis thnConfigAxisEtaRelative{"thnConfigAxisEtaRelative", {300, -1., 2.}, "(eta rec - eta true) / eta true axis"};
  ConfigurableAxis thnConfigAxisPhiRelative{"thnConfigAxisPhiRelative", {300, -1., 2.}, "(phi rec - phi true) / phi true axis"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {20, 0., 100.}, "centrality axis for the current event"};
  ConfigurableAxis thnConfigAxisMult{"thnConfigAxisMult", {60, 0., 60000.}, "multiplicity axis for the current event"};
  ConfigurableAxis thnConfigAxisDeltaEta{"thnConfigAxisDeltaEta", {100, -1, 1}, "delta eta axis"};
  ConfigurableAxis thnConfigAxisDeltaPhi{"thnConfigAxisDeltaPhi", {100, -1, 1}, "delta phi axis"};
  Configurable<bool> useCent{"useCent", 0, "flag to enable usage of centrality instead of multiplicity as axis."};

  struct : ConfigurableGroup {
    std::string prefix = "conversiontagging";
    Configurable<float> maxMinv{"maxMinv", 0.1f, "maximum invariant mass for tagging conversions."};
    Configurable<float> maxDeltaEta{"maxDeltaEta", 0.05f, "maximum delta eta between two clusters for tagging them as conversions."};
    Configurable<float> maxDeltaPhi{"maxDeltaPhi", 0.1f, "maximum delta phi between two clusters for tagging them as conversions."};
    Configurable<float> minRConv{"minRConv", 370.f, "minimum conversion Radius of two clusters for tagging them as conversions."};
    Configurable<float> maxRConv{"maxRConv", 430.f, "maximum conversion Radius of two clusters for tagging them as conversions."};
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

  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::MinClusters, aod::EMEMCClusterMCLabels>;

  using Colls = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMMCEventLabels>;

  using McColls = o2::soa::Join<o2::aod::EMMCEvents, o2::aod::BinnedGenPts>;
  using McParticles = EMMCParticles;

  PresliceOptional<EMCalPhotons> perCollisionEMC = o2::aod::emccluster::pmeventId;
  PresliceOptional<MinMTracks> perEMCClusterMT = o2::aod::mintm::minClusterId;
  PresliceOptional<MinMSTracks> perEMCClusterMS = o2::aod::mintm::minClusterId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  OutputObj<TTree> fConvTagTree{"convTagTree", OutputObjHandlingPolicy::AnalysisObject};

  float bMinv{}, bRConv{}, bDeltaEta{}, bDeltaPhi{}, bE1{}, bE2{}, bAsymmetry{}, bCentOrMult{}, bM021{}, bM022{}, bTime1{}, bTime2{}, bNCell1{}, bNCell2{}, bOpeningAngle{};
  int bTruthLabel{}, bCollisionIndex{};

  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb{};
  int mRunNumber{0};
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

    const AxisSpec thnAxisRConvRec{100, 0, 500, "#it{R}_{rec}"};
    const AxisSpec thnAxisRConvGen{100, 0, 500, "#it{R}_{gen}"};

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

    registry.add("EMCal/TrueConversion/hMassReco", "minv vs Rconv vs E1 vs E2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisRConvRec, thnAxisERec, thnAxisERec, thnAxisCentOrMult});
    registry.add("EMCal/TrueConversion/hRconvReco", "Rconv vs true RConv vs E1 vs E2", HistType::kTHnSparseF, {thnAxisRConvRec, thnAxisRConvGen, thnAxisERec, thnAxisERec, thnAxisCentOrMult});
    registry.add("EMCal/TrueConversion/hDeltaEtaDeltaPhi", "EMCal Delta Eta vs Delta vs E1 vs E2", HistType::kTHnSparseF, {thnAxisDeltaEta, thnAxisDeltaPhi, thnAxisERec, thnAxisERec, thnAxisCentOrMult});
    registry.add("EMCal/TrueConversion/hEnergyReco", "energy vs cent", HistType::kTH2D, {thnAxisERec, thnAxisCentOrMult});

    registry.addClone("EMCal/TrueConversion/", "EMCal/TrueConversionFromPi0/");
    registry.addClone("EMCal/TrueConversion/", "EMCal/TrueClusterFromPi0/");
    registry.addClone("EMCal/TrueConversion/", "EMCal/Background/");

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

    if (writeTree.value) {
      fConvTagTree.setObject(new TTree("convTagTree", "flattened conversion tagging features"));
      fConvTagTree->Branch("collisionIndex", &bCollisionIndex);
      fConvTagTree->Branch("minv", &bMinv);
      fConvTagTree->Branch("rConv", &bRConv);
      fConvTagTree->Branch("deltaEta", &bDeltaEta);
      fConvTagTree->Branch("deltaPhi", &bDeltaPhi);
      fConvTagTree->Branch("e1", &bE1);
      fConvTagTree->Branch("e2", &bE2);
      fConvTagTree->Branch("asymmetry", &bAsymmetry);
      fConvTagTree->Branch("centOrMult", &bCentOrMult);
      fConvTagTree->Branch("m021", &bM021);
      fConvTagTree->Branch("m022", &bM022);
      fConvTagTree->Branch("time1", &bTime1);
      fConvTagTree->Branch("time2", &bTime2);
      fConvTagTree->Branch("ncell1", &bNCell1);
      fConvTagTree->Branch("ncell2", &bNCell2);
      fConvTagTree->Branch("openingAngle", &bOpeningAngle);
      fConvTagTree->Branch("truthLabel", &bTruthLabel);
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
    std::vector<bool> wasMeassured(mcParticles.size(), false);
    EMBitFlags emcFlagsFromTrueMeson(clusters.size());
    EMBitFlags emcFlagsFromTrueMesonSameGamma(clusters.size());
    EMBitFlags emcFlagsFromTrueConversion(clusters.size());
    EMBitFlags emcFlagsTagging(clusters.size());
    EMBitFlags emcFlags(clusters.size());
    if (clusters.size() > 0) {
      fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds, &registry);
    }
    auto mcCluster1 = mcParticles.begin();
    auto mcCluster2 = mcParticles.begin();
    auto mcPhoton1 = mcParticles.begin();
    auto mcPhoton2 = mcParticles.begin();
    auto mcMother = mcParticles.begin();

    for (const auto& collision : collisions) {
      bCollisionIndex = collision.globalIndex(); // or emmceventId() if you want to split by MC event specifically
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

        const float deltaPhi = RecoDecay::constrainAngle(v1.Phi() - v2.Phi(), -o2::constants::math::PIHalf);
        const float deltaEta = v1.Eta() - v2.Eta();
        const float eT1 = v1.Et();
        const float eT2 = v2.Et();

        const float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));

        // calculate the conversion radius in cm. The formula expects radii in m, B field in T and energies in GeV
        const float RConv = 100.f * (EMCALRadius / 100.f - std::fabs(deltaPhi) / (0.15f * (std::fabs(dBz) / 10.f)) * (eT1 * eT2) / (eT1 + eT2));

        // tree values
        bMinv = vMeson.M();
        bRConv = RConv;
        bDeltaEta = deltaEta;
        bDeltaPhi = deltaPhi;
        bE1 = g1.e();
        bE2 = g2.e();
        bAsymmetry = std::fabs(g1.e() - g2.e()) / (g1.e() + g2.e());
        bCentOrMult = centOrMult;
        bM021 = g1.m02();
        bM022 = g2.m02();
        bTime1 = g1.time();
        bTime2 = g2.time();
        bNCell1 = g1.nCells();
        bNCell2 = g2.nCells();
        bOpeningAngle = openingAngle;

        // tag conversion pairs
        if (conversiontaggingcuts.maxMinv >= vMeson.M() && conversiontaggingcuts.maxDeltaEta >= std::fabs(deltaEta) && conversiontaggingcuts.minRConv <= RConv && RConv <= conversiontaggingcuts.maxRConv && conversiontaggingcuts.maxDeltaPhi >= std::fabs(deltaPhi)) {
          emcFlagsTagging.set(g1.globalIndex());
          emcFlagsTagging.set(g2.globalIndex());
        }

        // set MC particle cursors to the largest cluster contributor
        mcCluster1.setCursor(g1.emmcparticleIds()[0]);
        mcCluster2.setCursor(g2.emmcparticleIds()[0]);

        bool areLeptons = false;        // both clusters are e+/e-
        bool areConversionLegs = false; // both clusters are e+ + e- from one photon
        bool arePhotonFromPi0 = false;  // both cluster come from the same gamma which is from pi0 or eta
        bool areClusterFromPi0 = false; // both cluster come from the decay chain of a pi0 or eta excluding the same gamma case above

        if (std::abs(mcCluster1.pdgCode()) == kElectron && std::abs(mcCluster2.pdgCode()) == kElectron) {
          areLeptons = true;
        }
        if (o2::aod::pwgem::photonmeson::utils::mcutil::isMotherPDG(mcCluster1, kGamma) && o2::aod::pwgem::photonmeson::utils::mcutil::isMotherPDG(mcCluster2, kGamma)) {
          // both cluster come from a photon
          // reset the cursors to the legs!
          mcCluster1.setCursor(g1.emmcparticleIds()[0]);
          mcCluster2.setCursor(g2.emmcparticleIds()[0]);
          if (mcCluster1.mothersIds().size() > 0 && mcCluster2.mothersIds().size() > 0 && mcCluster1.mothersIds()[0] == mcCluster2.mothersIds()[0]) {
            // both cluster come from the same photon
            emcFlagsFromTrueConversion.set(g1.globalIndex());
            emcFlagsFromTrueConversion.set(g2.globalIndex());
            areConversionLegs = true;
          }
        }
        // reset the cursors to the legs!
        mcCluster1.setCursor(g1.emmcparticleIds()[0]);
        mcCluster2.setCursor(g2.emmcparticleIds()[0]);

        // check if the clusters are coming from a chain product of a pi0 or eta decay and then get the first level daughter of the pi0 or eta meson
        int photonid1 = o2::aod::pwgem::photonmeson::utils::mcutil::FindMotherInChain(mcCluster1, mcParticles, std::vector<int>{PDG_t::kPi0, Pdg::kEta});
        int photonid2 = o2::aod::pwgem::photonmeson::utils::mcutil::FindMotherInChain(mcCluster2, mcParticles, std::vector<int>{PDG_t::kPi0, Pdg::kEta});

        if (photonid1 >= 0 && photonid2 >= 0) {
          // set the cursors to the photonId
          mcPhoton1.setCursor(photonid1);
          mcPhoton2.setCursor(photonid2);
          mcMother.setCursor(mcPhoton1.mothersIds()[0]);
          if (mcMother.producedByGenerator()) {
            if (photonid1 == photonid2) {
              // the two clusters came from the same photon which came from a pi0 or eta. Most likely: pi0/eta -> γ + X -> e+e- + X
              emcFlagsFromTrueMesonSameGamma.set(g1.globalIndex());
              emcFlagsFromTrueMesonSameGamma.set(g2.globalIndex());
              arePhotonFromPi0 = true;
            } else if (mcPhoton1.mothersIds()[0] == mcPhoton2.mothersIds()[0]) {
              emcFlagsFromTrueMeson.set(g1.globalIndex());
              emcFlagsFromTrueMeson.set(g2.globalIndex());
              areClusterFromPi0 = true;
            }
          }
        }

        // do we have two conversions
        if (areLeptons && areConversionLegs) {
          registry.fill(HIST("EMCal/TrueConversion/hMassReco"), vMeson.M(), RConv, g1.e(), g2.e(), centOrMult);
          registry.fill(HIST("EMCal/TrueConversion/hRconvReco"), RConv, std::hypot(mcCluster1.vx(), mcCluster1.vy()), g1.e(), g2.e(), centOrMult);
          registry.fill(HIST("EMCal/TrueConversion/hDeltaEtaDeltaPhi"), deltaEta, deltaPhi, g1.e(), g2.e(), centOrMult);
          if (arePhotonFromPi0) {
            registry.fill(HIST("EMCal/TrueConversionFromPi0/hMassReco"), vMeson.M(), RConv, g1.e(), g2.e(), centOrMult);
            registry.fill(HIST("EMCal/TrueConversionFromPi0/hRconvReco"), RConv, std::hypot(mcCluster1.vx(), mcCluster1.vy()), g1.e(), g2.e(), centOrMult);
            registry.fill(HIST("EMCal/TrueConversionFromPi0/hDeltaEtaDeltaPhi"), deltaEta, deltaPhi, g1.e(), g2.e(), centOrMult);
          }
        }
        if (areClusterFromPi0) {
          registry.fill(HIST("EMCal/TrueClusterFromPi0/hMassReco"), vMeson.M(), RConv, g1.e(), g2.e(), centOrMult);
          registry.fill(HIST("EMCal/TrueClusterFromPi0/hRconvReco"), RConv, std::hypot(mcCluster1.vx(), mcCluster1.vy()), g1.e(), g2.e(), centOrMult);
          registry.fill(HIST("EMCal/TrueClusterFromPi0/hDeltaEtaDeltaPhi"), deltaEta, deltaPhi, g1.e(), g2.e(), centOrMult);
        }
        if (!areClusterFromPi0 && !areConversionLegs) {
          registry.fill(HIST("EMCal/Background/hMassReco"), vMeson.M(), RConv, g1.e(), g2.e(), centOrMult);
          registry.fill(HIST("EMCal/Background/hRconvReco"), RConv, std::hypot(mcCluster1.vx(), mcCluster1.vy()), g1.e(), g2.e(), centOrMult);
          registry.fill(HIST("EMCal/Background/hDeltaEtaDeltaPhi"), deltaEta, deltaPhi, g1.e(), g2.e(), centOrMult);
        }

        // final tree values plus filling
        bTruthLabel = (areLeptons && areConversionLegs) ? 0 : (areClusterFromPi0 ? 1 : 2);
        if (writeTree.value) {
          fConvTagTree->Fill();
        }
      } // pair loop
    } // collision loop

    if (collisions.size() <= 0) {
      return;
    }
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
      float centOrMult = getCentralityOrMultiplicity(collision);

      mcCluster1.setCursor(cluster.emmcparticleIds()[0]);
      int photonid1 = o2::aod::pwgem::photonmeson::utils::mcutil::FindMotherInChain(mcCluster1, mcParticles, std::vector<int>{PDG_t::kPi0, Pdg::kEta});
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
        registry.fill(HIST("EMCal/ConfusionMatrixConversionTagging"), emcFlagsTagging.test(cluster.globalIndex()) ? 0 : 1, static_cast<float>(TruthClass::Photon));
      }
      if (!emcFlagsFromTrueConversion.test(cluster.globalIndex())) {
        registry.fill(HIST("EMCal/ConfusionMatrixConversionTagging"), emcFlagsTagging.test(cluster.globalIndex()) ? 0 : 1, static_cast<float>(TruthClass::Conversion));
        registry.fill(HIST("EMCal/TrueConversion/hEnergyReco"), cluster.e(), centOrMult);
        if (!emcFlagsFromTrueMesonSameGamma.test(cluster.globalIndex()) && !alreadyCountedForMeson) {
          registry.fill(HIST("EMCal/ConfusionMatrixConversionTagging"), emcFlagsTagging.test(cluster.globalIndex()) ? 0 : 1, static_cast<float>(TruthClass::Pi0Conversion));
          registry.fill(HIST("EMCal/TrueConversionFromPi0/hEnergyReco"), cluster.e(), centOrMult);
        }
      }
      if (!emcFlagsFromTrueMeson.test(cluster.globalIndex()) && !alreadyCountedForMeson) {
        registry.fill(HIST("EMCal/ConfusionMatrixConversionTagging"), emcFlagsTagging.test(cluster.globalIndex()) ? 0 : 1, static_cast<float>(TruthClass::Pi0));
        registry.fill(HIST("EMCal/TrueClusterFromPi0/hEnergyReco"), cluster.e(), centOrMult);
      }
      if (emcFlagsFromTrueMeson.test(cluster.globalIndex()) && emcFlagsFromTrueConversion.test(cluster.globalIndex())) {
        registry.fill(HIST("EMCal/ConfusionMatrixConversionTagging"), emcFlagsTagging.test(cluster.globalIndex()) ? 0 : 1, static_cast<float>(TruthClass::Background));
        registry.fill(HIST("EMCal/Background/hEnergyReco"), cluster.e(), centOrMult);
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
      if (mcPart.daughtersIds().size() != 2) {
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
