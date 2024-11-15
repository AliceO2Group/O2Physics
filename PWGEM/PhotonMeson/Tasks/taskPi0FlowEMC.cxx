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
///
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include <numbers>
#include <iterator>
#include <array>
#include <string>
#include <map>
#include <algorithm>
#include <vector>
#include <tuple>
#include <utility>

#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/LorentzRotation.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Qvectors.h"

#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsEMCAL/Constants.h"
#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"

#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/emcalHistoDefinitions.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/NMHistograms.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photonmeson::photonpair;
using namespace o2::aod::pwgem::photon;
using namespace o2::aod::pwgem::dilepton::utils;

enum QvecEstimator { FT0M = 0,
                     FT0A = 1,
                     FT0C,
                     TPCPos,
                     TPCNeg,
                     TPCTot };

enum CentralityEstimator { None = 0,
                           CFT0A = 1,
                           CFT0C,
                           CFT0M,
                           NCentralityEstimators
};

struct EMfTaskPi0Flow {
  // configurable for flow
  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<int> qvecDetector{"qvecDetector", 0, "Detector for Q vector estimation (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5)"};
  Configurable<int> qvecSubADetector{"qvecSubADetector", 3, "Sub A Detector for Q vector estimation for resolution (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5)"};
  Configurable<int> qvecSubBDetector{"qvecSubBDetector", 4, "Sub B Detector for Q vector estimation for resolution (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5)"};
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3)"};
  Configurable<bool> saveEpResoHisto{"saveEpResoHisto", false, "Flag to save event plane resolution histogram"};
  Configurable<bool> saveSPResoHist{"saveSPResoHist", false, "Flag to save scalar product resolution histogram"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // configurable axis
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {200, 0.0, 0.4}, ""};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0., 20.}, ""};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {20, 0., 100.}, ""};
  ConfigurableAxis thnConfigAxisCosNPhi{"thnConfigAxisCosNPhi", {100, -1., 1.}, ""};
  ConfigurableAxis thnConfigAxisCosDeltaPhi{"thnConfigAxisCosDeltaPhi", {100, -1., 1.}, ""};
  ConfigurableAxis thnConfigAxisScalarProd{"thnConfigAxisScalarProd", {100, -5., 5.}, ""};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
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
    Configurable<int> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -1, "min. FT0C occupancy"};
    Configurable<int> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<float> cfgMinCent{"cfgMinCent", 0, "min. centrality (%)"};
    Configurable<float> cfgMaxCent{"cfgMaxCent", 90, "max. centrality (%)"};
    Configurable<bool> onlyKeepWeightedEvents{"onlyKeepWeightedEvents", false, "flag to keep only weighted events (for JJ MCs) and remove all MB events (with weight = 1)"};
    Configurable<bool> enableQA{"enableQA", false, "flag to turn QA plots on/off"};
  } eventcuts;

  EMCPhotonCut fEMCCut;
  struct : ConfigurableGroup {
    std::string prefix = "emccut_group";
    Configurable<float> EMC_minTime{"EMC_minTime", -25., "Minimum cluster time for EMCal time cut"};
    Configurable<float> EMC_maxTime{"EMC_maxTime", +30., "Maximum cluster time for EMCal time cut"};
    Configurable<float> EMC_minM02{"EMC_minM02", 0.1, "Minimum M02 for EMCal M02 cut"};
    Configurable<float> EMC_maxM02{"EMC_maxM02", 0.7, "Maximum M02 for EMCal M02 cut"};
    Configurable<float> EMC_minE{"EMC_minE", 0.7, "Minimum cluster energy for EMCal energy cut"};
    Configurable<int> EMC_minNCell{"EMC_minNCell", 1, "Minimum number of cells per cluster for EMCal NCell cut"};
    Configurable<std::vector<float>> EMC_TM_Eta{"EMC_TM_Eta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> EMC_TM_Phi{"EMC_TM_Phi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<float> EMC_Eoverp{"EMC_Eoverp", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
    Configurable<bool> EMC_UseExoticCut{"EMC_UseExoticCut", true, "FLag to use the EMCal exotic cluster cut"};
    Configurable<bool> EMC_UseTM{"EMC_UseTM", false, "flag to use EMCal track matching cut or not"};
    Configurable<bool> enableQA{"enableQA", false, "flag to turn QA plots on/off"};
  } emccuts;

  struct : ConfigurableGroup {
    std::string prefix = "meson";
    Configurable<float> minOpenAngle{"minOpenAngle", 0.0202, "apply min opening angle. Default value one EMCal cell"};
    Configurable<float> minTanThetadPhi{"minTanThetadPhi", 4., "apply min opening angle in delta theta delta phi to cut on late conversion"};
    Configurable<float> maxEnergyAsymmetry{"maxEnergyAsymmetry", 1., "apply max energy asymmetry for meson candidate"};
    Configurable<bool> enableQA{"enableQA", false, "flag to turn QA plots on/off"};
    ConfigurableAxis thConfigAxisTanThetaPhi{"thConfigAxisTanThetaPhi", {180, -90.f, 90.f}, ""};
  } mesonConfig;

  struct : ConfigurableGroup {
    std::string prefix = "event-mixing";
    ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "Mixing bins - centrality"};
    ConfigurableAxis ConfEPBins{"ConfEPBins", {8, -M_PI / 2, +M_PI / 2}, "Mixing bins - event plane angle"};
    ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, 0, 100, 500, 1000, 2000}, "Mixing bins - occupancy"};
    Configurable<int> ConfMixingDepth{"ConfMixingDepth", 2, "Mixing depth"};
  } mixingConfig;

  SliceCache cache;
  EventPlaneHelper epHelper;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  Filter clusterFilter = aod::skimmedcluster::time >= emccuts.EMC_minTime && aod::skimmedcluster::time <= emccuts.EMC_maxTime && aod::skimmedcluster::m02 >= emccuts.EMC_minM02 && aod::skimmedcluster::m02 <= emccuts.EMC_maxM02 && skimmedcluster::e >= emccuts.EMC_minE;
  Filter collisionFilter = aod::evsel::sel8 && nabs(aod::collision::posZ) <= eventcuts.cfgZvtxMax && aod::evsel::trackOccupancyInTimeRange <= eventcuts.cfgTrackOccupancyMax && aod::evsel::trackOccupancyInTimeRange >= eventcuts.cfgTrackOccupancyMin && aod::evsel::ft0cOccupancyInTimeRange <= eventcuts.cfgFT0COccupancyMax && aod::evsel::ft0cOccupancyInTimeRange >= eventcuts.cfgFT0COccupancyMin;
  using FilteredEMCalPhotons = soa::Filtered<soa::Join<aod::EMCEMEventIds, aod::SkimEMCClusters>>;
  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::SkimEMCClusters>;
  using FilteredCollsWithQvecs = soa::Filtered<soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>>;
  using CollsWithQvecs = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;

  Preslice<EMCalPhotons> perCollision_emc = aod::emccluster::emeventId;
  Preslice<FilteredEMCalPhotons> perCollision_emc_filtered = aod::emccluster::emeventId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void DefineEMEventCut()
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

  void DefineEMCCut()
  {
    fEMCCut = EMCPhotonCut("fEMCCut", "fEMCCut");
    const float a = emccuts.EMC_TM_Eta->at(0);
    const float b = emccuts.EMC_TM_Eta->at(1);
    const float c = emccuts.EMC_TM_Eta->at(2);

    const float d = emccuts.EMC_TM_Phi->at(0);
    const float e = emccuts.EMC_TM_Phi->at(1);
    const float f = emccuts.EMC_TM_Phi->at(2);
    LOGF(info, "EMCal track matching parameters : a = %f, b = %f, c = %f, d = %f, e = %f, f = %f", a, b, c, d, e, f);
    fEMCCut.SetTrackMatchingEta([a, b, c](float pT) { return a + pow(pT + b, c); });
    fEMCCut.SetTrackMatchingPhi([d, e, f](float pT) { return d + pow(pT + e, f); });
    fEMCCut.SetMinEoverP(emccuts.EMC_Eoverp);

    fEMCCut.SetMinE(emccuts.EMC_minE);
    fEMCCut.SetMinNCell(emccuts.EMC_minNCell);
    fEMCCut.SetM02Range(emccuts.EMC_minM02, emccuts.EMC_maxM02);
    fEMCCut.SetTimeRange(emccuts.EMC_minTime, emccuts.EMC_maxTime);
    fEMCCut.SetUseExoticCut(emccuts.EMC_UseExoticCut);
  }

  void init(InitContext&)
  {
    if (harmonic != 2 && harmonic != 3) {
      LOG(info) << "Harmonic was set to " << harmonic << " but can only be 2 or 3!";
    }

    DefineEMEventCut();
    DefineEMCCut();
    fEMCCut.SetUseTM(emccuts.EMC_UseTM); // disables TM
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&registry);

    // Load EMCal geometry
    o2::emcal::Geometry::GetInstanceFromRunNumber(300000);

    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality (%)"};
    const AxisSpec thnAxisCosNPhi{thnConfigAxisCosNPhi, Form("cos(%d#varphi)", harmonic.value)};
    const AxisSpec thnAxisCosDeltaPhi{thnConfigAxisCosDeltaPhi, Form("cos(%d(#varphi - #Psi_{sub}))", harmonic.value)};
    const AxisSpec thnAxisScalarProd{thnConfigAxisScalarProd, "SP"};
    const AxisSpec thAxisTanThetaPhi{mesonConfig.thConfigAxisTanThetaPhi, "atan(#Delta#theta/#Delta#varphi)"};
    const AxisSpec thAxisClusterEnergy{thnConfigAxisPt, "#it{E} (GeV)"};
    const AxisSpec thAxisAlpha{100, -1., +1, "#alpha"};
    const AxisSpec thAxisMult{1000, 0., +1000, "#it{N}_{ch}"};
    const AxisSpec thAxisEnergy{1000, 0., 100., "#it{E}_{clus} (GeV)"};
    const AxisSpec thAxisTime{1500, -600, 900, "#it{t}_{cl} (ns)"};
    const AxisSpec thAxisEta{160, -0.8, 0.8, "#eta"};
    const AxisSpec thAxisPhi{72, 0, 2 * 3.14159, "phi"};
    const AxisSpec thAxisNCell{17664, 0.5, +17664.5, "#it{N}_{cell}"};

    const AxisSpec thAxisPsi{360 / harmonic, 0.f, 2. / harmonic * M_PI, Form("#Psi_{%d}", harmonic.value)};

    registry.add("hSparsePi0Flow", "THn for SP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCent, thnAxisScalarProd});
    registry.add("hSparseBkgFlow", "THn for SP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCent, thnAxisScalarProd});
    auto hClusterCuts = registry.add<TH1>("hClusterCuts", "hClusterCuts;;Counts", kTH1D, {{6, 0.5, 6.5}}, false);
    hClusterCuts->GetXaxis()->SetBinLabel(1, "in");
    hClusterCuts->GetXaxis()->SetBinLabel(2, "opening angle");
    hClusterCuts->GetXaxis()->SetBinLabel(3, "#it{M}_{#gamma#gamma}");
    hClusterCuts->GetXaxis()->SetBinLabel(4, "#it{p}_{T}");
    hClusterCuts->GetXaxis()->SetBinLabel(5, "conversion cut");
    hClusterCuts->GetXaxis()->SetBinLabel(6, "out");

    if (saveSPResoHist) {
      registry.add("spReso/hSpResoFT0cFT0a", "hSpResoFT0cFT0a; centrality; Q_{FT0c} #bullet Q_{FT0a}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("spReso/hSpResoFT0cTPCpos", "hSpResoFT0cTPCpos; centrality; Q_{FT0c} #bullet Q_{TPCpos}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("spReso/hSpResoFT0cTPCneg", "hSpResoFT0cTPCneg; centrality; Q_{FT0c} #bullet Q_{TPCneg}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("spReso/hSpResoFT0cTPCtot", "hSpResoFT0cTPCtot; centrality; Q_{FT0c} #bullet Q_{TPCtot}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("spReso/hSpResoFT0aTPCpos", "hSpResoFT0aTPCpos; centrality; Q_{FT0a} #bullet Q_{TPCpos}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("spReso/hSpResoFT0aTPCneg", "hSpResoFT0aTPCneg; centrality; Q_{FT0a} #bullet Q_{TPCneg}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("spReso/hSpResoFT0aTPCtot", "hSpResoFT0aTPCtot; centrality; Q_{FT0m} #bullet Q_{TPCtot}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("spReso/hSpResoFT0mTPCpos", "hSpResoFT0mTPCpos; centrality; Q_{FT0m} #bullet Q_{TPCpos}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("spReso/hSpResoFT0mTPCneg", "hSpResoFT0mTPCneg; centrality; Q_{FT0m} #bullet Q_{TPCneg}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("spReso/hSpResoFT0mTPCtot", "hSpResoFT0mTPCtot; centrality; Q_{FT0m} #bullet Q_{TPCtot}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("spReso/hSpResoTPCposTPCneg", "hSpResoTPCposTPCneg; centrality; Q_{TPCpos} #bullet Q_{TPCneg}", {HistType::kTProfile, {thnAxisCent}});
    }

    if (saveEpResoHisto) {
      registry.add("hEventPlaneAngleFT0M", "hEventPlaneAngleFT0M", HistType::kTH2D, {thnAxisCent, thAxisPsi});
      registry.add("hEventPlaneAngleTPCpos", "hEventPlaneAngleTPCpos", HistType::kTH2D, {thnAxisCent, thAxisPsi});
      registry.add("hEventPlaneAngleTPCneg", "hEventPlaneAngleTPCneg", HistType::kTH2D, {thnAxisCent, thAxisPsi});
      registry.add("epReso/hEpResoFT0cFT0a", "hEpResoFT0cFT0a; centrality; #Delta#Psi_{sub}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("epReso/hEpResoFT0cTPCpos", "hEpResoFT0cTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("epReso/hEpResoFT0cTPCneg", "hEpResoFT0cTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("epReso/hEpResoFT0cTPCtot", "hEpResoFT0cTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("epReso/hEpResoFT0aTPCpos", "hEpResoFT0aTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("epReso/hEpResoFT0aTPCneg", "hEpResoFT0aTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("epReso/hEpResoFT0aTPCtot", "hEpResoFT0aTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("epReso/hEpResoFT0mTPCpos", "hEpResoFT0mTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("epReso/hEpResoFT0mTPCneg", "hEpResoFT0mTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("epReso/hEpResoFT0mTPCtot", "hEpResoFT0mTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTProfile, {thnAxisCent}});
      registry.add("epReso/hEpResoTPCposTPCneg", "hEpResoTPCposTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTProfile, {thnAxisCent}});
    }
    if (eventcuts.enableQA) {
      auto hCollisionEMCCheck = registry.add<TH1>("hCollisionEMCCheck", "collision counter;;Counts", kTH1D, {{7, 0.5, 7.5}}, false);
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(1, "all");
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(2, "EMC MB Readout");
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(3, "has clusters");
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(4, "EMC MB Readout & has clusters");
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(5, "EMC MB Readout but no clusters");
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(6, "No EMC MB Readout but has clusters");
      hCollisionEMCCheck->GetXaxis()->SetBinLabel(7, "No EMC MB Readout and no clusters");
      registry.add("LED/hMult", "multiplicity in LED events", HistType::kTH1D, {thAxisMult});
      registry.add("LED/hClusterEtaPhi", "hClusterEtaPhi", HistType::kTH2D, {thAxisPhi, thAxisEta});
      registry.add("LED/clusterTimeVsE", "Cluster time vs energy", HistType::kTH2D, {thAxisTime, thAxisEnergy});
      registry.add("LED/hNCell", "hNCell", HistType::kTH1D, {thAxisNCell});
    }

    if (emccuts.enableQA) {
      registry.add("hEClusterBefore", "Histo for cluster energy before cuts", HistType::kTH1D, {thAxisClusterEnergy});
      registry.add("hEClusterAfter", "Histo for cluster energy after cuts", HistType::kTH1D, {thAxisClusterEnergy});
    }

    if (mesonConfig.enableQA) {
      registry.add("hInvMassPt", "Histo for inv pair mass vs pt", HistType::kTH2D, {thnAxisInvMass, thnAxisPt});
      registry.add("hTanThetaPhi", "Histo for identification of conversion cluster", HistType::kTH2D, {thnAxisInvMass, thAxisTanThetaPhi});
      registry.add("hAlphaPt", "Histo of meson asymmetry vs pT", HistType::kTH2D, {thAxisAlpha, thnAxisPt});
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    LOG(info) << "thnConfigAxisInvMass.value[1] = " << thnConfigAxisInvMass.value[1] << " thnConfigAxisInvMass.value.back() = " << thnConfigAxisInvMass.value.back();
    LOG(info) << "thnConfigAxisPt.value[1] = " << thnConfigAxisPt.value[1] << " thnConfigAxisPt.value.back() = " << thnConfigAxisPt.value.back();
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
    if (std::abs(deltaPsi) > constants::math::PI / harmonic) {
      if (deltaPsi > 0.)
        deltaPsi -= constants::math::TwoPI / harmonic;
      else
        deltaPsi += constants::math::TwoPI / harmonic;
    }
    return deltaPsi;
  }

  /// Fill THnSparse
  /// \param mass is the invariant mass of the candidate
  /// \param pt is the transverse momentum of the candidate
  /// \param cent is the centrality of the collision
  /// \param sp is the scalar product
  void fillThn(float& mass,
               float& pt,
               float& cent,
               float& sp)
  {
    registry.fill(HIST("hSparsePi0Flow"), mass, pt, cent, sp);
  }

  /// Get the centrality
  /// \param collision is the collision with the centrality information
  float getCentrality(CollsWithQvecs::iterator const& collision)
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
  std::vector<float> getAllQvec(CollsWithQvecs::iterator const& collision)
  {
    // Retrieve the Q vectors using the helper function for each detector
    auto [xQVecMain, yQVecMain] = getQvec(collision, qvecDetector);
    auto [xQVecSubA, yQVecSubA] = getQvec(collision, qvecSubADetector);
    auto [xQVecSubB, yQVecSubB] = getQvec(collision, qvecSubBDetector);

    return {xQVecMain, yQVecMain, xQVecSubA, yQVecSubA, xQVecSubB, yQVecSubB};
  }

  /// Get the Q vector
  /// \param collision is the collision with the Q vector information
  std::pair<float, float> getQvec(CollsWithQvecs::iterator const& collision, int detector)
  {
    float xQVec = -999.f;
    float yQVec = -999.f;

    switch (detector) {
      case QvecEstimator::FT0M:
        if (harmonic == 2) {
          xQVec = collision.q2xft0m();
          yQVec = collision.q2yft0m();
        } else if (harmonic == 3) {
          xQVec = collision.q3xft0m();
          yQVec = collision.q3yft0m();
        }
        break;
      case QvecEstimator::FT0A:
        if (harmonic == 2) {
          xQVec = collision.q2xft0a();
          yQVec = collision.q2yft0a();
        } else if (harmonic == 3) {
          xQVec = collision.q3xft0a();
          yQVec = collision.q3yft0a();
        }
        break;
      case QvecEstimator::FT0C:
        if (harmonic == 2) {
          xQVec = collision.q2xft0c();
          yQVec = collision.q2yft0c();
        } else if (harmonic == 3) {
          xQVec = collision.q3xft0c();
          yQVec = collision.q3yft0c();
        }
        break;
      case QvecEstimator::TPCPos:
        if (harmonic == 2) {
          xQVec = collision.q2xbpos();
          yQVec = collision.q2ybpos();
        } else if (harmonic == 3) {
          xQVec = collision.q3xbpos();
          yQVec = collision.q3ybpos();
        }
        break;
      case QvecEstimator::TPCNeg:
        if (harmonic == 2) {
          xQVec = collision.q2xbneg();
          yQVec = collision.q2ybneg();
        } else if (harmonic == 3) {
          xQVec = collision.q3xbneg();
          yQVec = collision.q3ybneg();
        }
        break;
      case QvecEstimator::TPCTot:
        if (harmonic == 2) {
          xQVec = collision.q2xbtot();
          yQVec = collision.q2ybtot();
        } else if (harmonic == 3) {
          xQVec = collision.q3xbtot();
          yQVec = collision.q3ybtot();
        }
        break;
      default:
        LOG(warning) << "Q vector estimator not valid. Falling back to FT0M";
        if (harmonic == 2) {
          xQVec = collision.q2xft0m();
          yQVec = collision.q2yft0m();
        } else if (harmonic == 3) {
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
    for (auto& QVec : QVecs) {
      if (std::fabs(QVec) > 20.f) {
        isgood = false;
        break;
      }
    }
    return isgood;
  }

  /// \brief Calculate background using rotation background method
  template <typename TPhotons>
  void RotationBackground(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photon1, ROOT::Math::PtEtaPhiMVector photon2, TPhotons const& photons_coll, unsigned int ig1, unsigned int ig2, CollsWithQvecs::iterator const& collision)
  {
    // if less than 3 clusters are present skip event since we need at least 3 clusters
    if (photons_coll.size() < 3) {
      return;
    }

    auto [xQVec, yQVec] = getQvec(collision, qvecDetector);
    float cent = getCentrality(collision);

    const float rotationAngle = M_PI / 2.0; // rotaion angle 90 degree
    ROOT::Math::AxisAngle rotationAxis(meson.Vect(), rotationAngle);
    ROOT::Math::Rotation3D rotationMatrix(rotationAxis);
    photon1 = rotationMatrix * photon1;
    photon2 = rotationMatrix * photon2;

    for (auto& photon : photons_coll) {
      if (photon.globalIndex() == ig1 || photon.globalIndex() == ig2) {
        // only combine rotated photons with other photons
        continue;
      }
      if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(photon))) {
        continue;
      }

      ROOT::Math::PtEtaPhiMVector photon3(photon.pt(), photon.eta(), photon.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector mother1 = photon1 + photon3;
      ROOT::Math::PtEtaPhiMVector mother2 = photon2 + photon3;

      float openingAngle1 = std::acos(photon1.Vect().Dot(photon3.Vect()) / (photon1.P() * photon3.P()));
      float openingAngle2 = std::acos(photon2.Vect().Dot(photon3.Vect()) / (photon2.P() * photon3.P()));

      int iCellID_photon1 = 0;
      int iCellID_photon2 = 0;

      float cosNPhi1 = std::cos(harmonic * mother1.Phi());
      float sinNPhi1 = std::sin(harmonic * mother1.Phi());
      float scalprodCand1 = cosNPhi1 * xQVec + sinNPhi1 * yQVec;

      float cosNPhi2 = std::cos(harmonic * mother2.Phi());
      float sinNPhi2 = std::sin(harmonic * mother2.Phi());
      float scalprodCand2 = cosNPhi2 * xQVec + sinNPhi2 * yQVec;

      try {
        iCellID_photon1 = o2::emcal::Geometry::GetInstance()->GetAbsCellIdFromEtaPhi(photon1.Eta(), photon1.Phi());
      } catch (o2::emcal::InvalidPositionException& e) {
        iCellID_photon1 = -1;
      }
      try {
        iCellID_photon2 = o2::emcal::Geometry::GetInstance()->GetAbsCellIdFromEtaPhi(photon2.Eta(), photon2.Phi());
      } catch (o2::emcal::InvalidPositionException& e) {
        iCellID_photon2 = -1;
      }

      if (openingAngle1 > mesonConfig.minOpenAngle && iCellID_photon1 > 0 && thnConfigAxisInvMass.value[1] <= mother1.M() && thnConfigAxisInvMass.value.back() >= mother1.M() && thnConfigAxisPt.value[1] > mother1.Pt() && thnConfigAxisPt.value.back() < mother1.Pt()) {
        registry.fill(HIST("hSparseBkgFlow"), mother1.M(), mother1.Pt(), cent, scalprodCand1);
      }
      if (openingAngle2 > mesonConfig.minOpenAngle && iCellID_photon2 > 0 && thnConfigAxisInvMass.value[1] <= mother2.M() && thnConfigAxisInvMass.value.back() >= mother2.M() && thnConfigAxisPt.value[1] > mother2.Pt() && thnConfigAxisPt.value.back() < mother2.Pt()) {
        registry.fill(HIST("hSparseBkgFlow"), mother2.M(), mother2.Pt(), cent, scalprodCand2);
      }
    }
  }

  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param meson are the selected candidates
  void runFlowAnalysis(CollsWithQvecs::iterator const& collision, ROOT::Math::PtEtaPhiMVector const& meson)
  {
    auto [xQVec, yQVec] = getQvec(collision, qvecDetector);
    float cent = getCentrality(collision);

    float massCand = meson.M();
    float ptCand = meson.Pt();
    float phiCand = meson.Phi();

    float cosNPhi = std::cos(harmonic * phiCand);
    float sinNPhi = std::sin(harmonic * phiCand);
    float scalprodCand = cosNPhi * xQVec + sinNPhi * yQVec;

    fillThn(massCand, ptCand, cent, scalprodCand);
    return;
  }

  // Pi0 from EMCal
  void processEMCal(CollsWithQvecs const& collisions, EMCalPhotons const& clusters)
  {
    for (auto& collision : collisions) {
      auto photons_per_collision = clusters.sliceBy(perCollision_emc, collision.globalIndex());

      if (eventcuts.enableQA) {
        // TODO: check EMCal NCells in collisions without EMC Readout
        registry.fill(HIST("hCollisionEMCCheck"), 1.); // all
        if (collision.alias_bit(kTVXinEMC) == true) {
          registry.fill(HIST("hCollisionEMCCheck"), 2.); // has EMC read out
          if (photons_per_collision.size() > 0) {
            registry.fill(HIST("hCollisionEMCCheck"), 3.); // has EMC cluster
            registry.fill(HIST("hCollisionEMCCheck"), 4.); // has EMC read out and clusters
          } else {
            registry.fill(HIST("hCollisionEMCCheck"), 5.); // has EMC read out but no clusters
          }
        } else {
          if (photons_per_collision.size() > 0) {
            registry.fill(HIST("hCollisionEMCCheck"), 3.); // has EMC cluster
            registry.fill(HIST("hCollisionEMCCheck"), 6.); // has no EMC read out and clusters
            registry.fill(HIST("LED/hMult"), collision.multFT0C());
            for (auto& photon : photons_per_collision) {
              registry.fill(HIST("LED/hClusterEtaPhi"), photon.phi(), photon.eta());
              registry.fill(HIST("LED/clusterTimeVsE"), photon.time(), photon.e());
              registry.fill(HIST("LED/hNCell"), photon.nCells());
            }
          } else {
            registry.fill(HIST("hCollisionEMCCheck"), 7.); // has no EMC read out and no clusters
          }
        }
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&registry, collision);
      if (!(fEMEventCut.IsSelected(collision))) {
        // general event selection
        continue;
        continue;
      }
      if (!(eventcuts.cfgTrackOccupancyMin <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < eventcuts.cfgTrackOccupancyMax)) {
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
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&registry, collision);
      registry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // accepted
      registry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // accepted

      if (emccuts.enableQA) {
        for (auto& photon : photons_per_collision) {
          registry.fill(HIST("hEClusterBefore"), photon.e()); // before cuts
          if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(photon))) {
            continue;
          }
          registry.fill(HIST("hEClusterAfter"), photon.e()); // accepted after cuts
        }
      }
      for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons_per_collision, photons_per_collision))) {
        if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(g1)) || !(fEMCCut.IsSelected<EMCalPhotons::iterator>(g2))) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector vMeson = v1 + v2;

        RotationBackground<EMCalPhotons>(vMeson, v1, v2, photons_per_collision, g1.globalIndex(), g2.globalIndex(), collision);

        float dTheta = v1.Theta() - v2.Theta();
        float dPhi = v1.Phi() - v2.Phi();
        float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));
        registry.fill(HIST("hClusterCuts"), 1);
        if (openingAngle <= mesonConfig.minOpenAngle) {
          registry.fill(HIST("hClusterCuts"), 2);
          continue;
        }
        if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
          registry.fill(HIST("hClusterCuts"), 3);
          continue;
        }
        if (thnConfigAxisPt.value[1] > vMeson.Pt() || thnConfigAxisPt.value.back() < vMeson.Pt()) {
          registry.fill(HIST("hClusterCuts"), 4);
          continue;
        }
        if (mesonConfig.enableQA) {
          registry.fill(HIST("hInvMassPt"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("hTanThetaPhi"), vMeson.M(), getAngleDegree(atan(dTheta / dPhi)));
          registry.fill(HIST("hAlphaPt"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        if (mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(atan(dTheta / dPhi)))) {
          registry.fill(HIST("hClusterCuts"), 5);
          continue;
        }
        registry.fill(HIST("hClusterCuts"), 6);
        runFlowAnalysis(collision, vMeson);
      }
    }
  }
  PROCESS_SWITCH(EMfTaskPi0Flow, processEMCal, "Process EMCal Pi0 candidates", true);

  // Pi0 from EMCal
  void processEMCalMixed(FilteredCollsWithQvecs const& collisions, FilteredEMCalPhotons const& clusters)
  {
    auto getClustersSize =
      [&clusters, this](FilteredCollsWithQvecs::iterator const& col) {
        auto associatedClusters = clusters.sliceByCached(emccluster::emeventId, col.globalIndex(), this->cache); // it's cached, so slicing/grouping happens only once
        return associatedClusters.size();
      };

    using BinningType = FlexibleBinningPolicy<std::tuple<decltype(getClustersSize)>, aod::collision::PosZ, aod::cent::CentFT0C, emevent::EP2FT0M<emevent::Q2xFT0M, emevent::Q2yFT0M>>;
    BinningType binningWithLambda{{getClustersSize}, {mixingConfig.ConfVtxBins, mixingConfig.ConfCentBins, mixingConfig.ConfEPBins}, true};

    auto clustersTuple = std::make_tuple(clusters);
    SameKindPair<FilteredCollsWithQvecs, FilteredEMCalPhotons, BinningType> pair{binningWithLambda, mixingConfig.ConfMixingDepth, -1, collisions, clustersTuple, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

    for (auto& [c1, clusters1, c2, clusters2] : pair) {
      if (!(fEMEventCut.IsSelected(c1)) || !(fEMEventCut.IsSelected(c2))) {
        // general event selection
        continue;
      }
      if (!(eventcuts.cfgTrackOccupancyMin <= c1.trackOccupancyInTimeRange() && c1.trackOccupancyInTimeRange() < eventcuts.cfgTrackOccupancyMax) || !(eventcuts.cfgTrackOccupancyMin <= c2.trackOccupancyInTimeRange() && c2.trackOccupancyInTimeRange() < eventcuts.cfgTrackOccupancyMax)) {
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
      for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(clusters1, clusters2))) {
        if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(g1)) || !(fEMCCut.IsSelected<EMCalPhotons::iterator>(g2))) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
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
        if (thnConfigAxisInvMass.value[1] > vMeson.M() || thnConfigAxisInvMass.value.back() < vMeson.M()) {
          registry.fill(HIST("hClusterCuts"), 3);
          continue;
        }
        if (thnConfigAxisPt.value[1] > vMeson.Pt() || thnConfigAxisPt.value.back() < vMeson.Pt()) {
          registry.fill(HIST("hClusterCuts"), 4);
          continue;
        }
        if (mesonConfig.enableQA) {
          registry.fill(HIST("hInvMassPt"), vMeson.M(), vMeson.Pt());
          registry.fill(HIST("hTanThetaPhi"), vMeson.M(), getAngleDegree(atan(dTheta / dPhi)));
          registry.fill(HIST("hAlphaPt"), (v1.E() - v2.E()) / (v1.E() + v2.E()), vMeson.Pt());
        }
        if (mesonConfig.minTanThetadPhi > std::fabs(getAngleDegree(atan(dTheta / dPhi)))) {
          registry.fill(HIST("hClusterCuts"), 5);
          continue;
        }
        registry.fill(HIST("hClusterCuts"), 6);
        runFlowAnalysis(c1, vMeson);
      }
    }
  }
  PROCESS_SWITCH(EMfTaskPi0Flow, processEMCalMixed, "Process EMCal Pi0 mixed event candidates", false);

  // Resolution
  void processResolution(CollsWithQvecs::iterator const& collision)
  {
    o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&registry, collision);
    if (!(fEMEventCut.IsSelected(collision))) {
      // no selection on the centrality is applied on purpose to allow for the resolution study in post-processing
      return;
    }
    if (!(eventcuts.cfgTrackOccupancyMin <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < eventcuts.cfgTrackOccupancyMax)) {
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
    if (harmonic == 2) {
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
    } else if (harmonic == 3) {
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
    }
  }
  PROCESS_SWITCH(EMfTaskPi0Flow, processResolution, "Process resolution", false);

}; // End struct EMfTaskPi0Flow

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<EMfTaskPi0Flow>(cfgc)};
}
