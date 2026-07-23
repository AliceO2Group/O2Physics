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

/// \file flowGfwNonflow.cxx
/// \brief Task to analyse scaled non-flow subtraction of flow cumulants
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch
/// \since 06/07/2026

#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/FlowPtContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWConfig.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TF1.h>
#include <TH1.h>
#include <TH3.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TString.h>

#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::analysis::genericframework;

constexpr float DefaultMagneticFieldCut = 99999.f;
static constexpr std::array<std::array<double, 2>, 4> LongArrayDouble = {{{{-0.8, -0.5}}, {{0.5, 0.8}}, {{-2, -2}}, {{-2, -2}}}};

struct FlowGfwNonflow {
  Configurable<int> cfgNbootstrap{"cfgNbootstrap", 10, "Number of subsamples"};
  Configurable<int> cfgMpar{"cfgMpar", 4, "Highest order of pt-pt correlations"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A, 4:NTPV, 5:NGlobal, 6:MFT"};
  Configurable<bool> cfgUseNch{"cfgUseNch", false, "Do correlations as function of Nch"};
  Configurable<int> cfgUseNchCorrection{"cfgUseNchCorrection", 1, "Use correction for Nch; 0: Use size of tracks table, 1: Use efficiency-corrected Nch values, 2: Use uncorrected Nch values"};
  Configurable<bool> cfgRunByRun{"cfgRunByRun", false, "Use run-by-run NUA"};
  Configurable<bool> cfgFillQA{"cfgFillQA", false, "Fill QA histograms"};
  Configurable<bool> cfgUseCentralMoments{"cfgUseCentralMoments", true, "Use central moments in vn-pt calculations"};
  Configurable<bool> cfgUseMultiplicityFlowWeights{"cfgUseMultiplicityFlowWeights", true, "Enable or disable the use of multiplicity-based event weighting"};
  struct : ConfigurableGroup {
    Configurable<std::string> cfgEfficiencyPath{"cfgEfficiencyPath", "", "CCDB path to efficiency object"};
    Configurable<bool> cfgUse2DEfficiency{"cfgUse2DEfficiency", false, "Toggle the use of 2D (pt, centrality) efficiency versus centrality integrated efficiency"};
    Configurable<std::string> cfgAcceptancePath{"cfgAcceptancePath", "", "CCDB path to acceptance object"};
  } cfgCorrections;
  struct : ConfigurableGroup {
    Configurable<std::pair<float, float>> cfgEta{"cfgEta", {-0.8, 0.8}, "eta cut"};
    Configurable<std::pair<float, float>> cfgEtaNch{"cfgEtaNch", {-0.5, 0.5}, "eta cut for nch selection"};
    Configurable<std::pair<float, float>> cfgEtaPtPt{"cfgEtaPtPt", {-0.5, 0.5}, "eta for pt-pt correlation"};
    Configurable<std::pair<float, float>> cfgPtCut{"cfgPtCut", {0.2, 5.0}, "minimum and maximum pt (GeV/c)"};
  } cfgKinematics;
  Configurable<LabeledArray<double>> cfgPtPtGaps{"cfgPtPtGaps", {LongArrayDouble.front().data(), 4, 2, {"subevent 1", "subevent 2", "subevent 3", "subevent 4"}, {"etamin", "etamax"}}, "{etamin,etamax} for all ptpt-subevents"};
  struct : ConfigurableGroup {
    Configurable<float> cfgDCAxyNSigma{"cfgDCAxyNSigma", 7, "Cut on number of sigma deviations from expected DCA in the transverse direction"};
    Configurable<std::string> cfgDCAxyPtDep{"cfgDCAxyPtDep", "(0.0105 + 0.0350/(x^1.1))", "Functional form of pt-dependent 7 sigma DCAxy cut"};
    Configurable<std::string> cfgDCAzPtDep{"cfgDCAzPtDep", "(0.0105 + 0.0350/(x^1.1))", "Functional form of pt-dependent DCAz cut"};
    Configurable<float> cfgDCAz{"cfgDCAz", 2, "Cut on DCA in the longitudinal direction (cm)"};
    Configurable<float> cfgNTPCCls{"cfgNTPCCls", 50, "Cut on number of TPC clusters found"};
    Configurable<float> cfgNTPCXrows{"cfgNTPCXrows", 70, "Cut on number of TPC crossed rows"};
    Configurable<float> cfgMinNITSCls{"cfgMinNITSCls", 5, "Cut on minimum number of ITS clusters found"};
    Configurable<float> cfgChi2PrITSCls{"cfgChi2PrITSCls", 36, "Cut on chi^2 per ITS clusters found"};
    Configurable<float> cfgChi2PrTPCCls{"cfgChi2PrTPCCls", 2.5, "Cut on chi^2 per TPC clusters found"};
    Configurable<bool> cfgTPCSectorCut{"cfgTPCSectorCut", false, "Cut on pt-phi distribution"};
  } cfgTrackCuts;
  struct : ConfigurableGroup {
    Configurable<bool> cfgNoSameBunchPileupCut{"cfgNoSameBunchPileupCut", true, "NoSameBunchPileupCut"};
    Configurable<bool> cfgIsGoodZvtxFT0vsPV{"cfgIsGoodZvtxFT0vsPV", true, "IsGoodZvtxFT0vsPV"};
    Configurable<bool> cfgIsGoodITSLayersAll{"cfgIsGoodITSLayersAll", true, "IsGoodITSLayersAll"};
    Configurable<bool> cfgNoCollInTimeRangeStandard{"cfgNoCollInTimeRangeStandard", true, "NoCollInTimeRangeStandard"};
    Configurable<bool> cfgNoCollInRofStandard{"cfgNoCollInRofStandard", true, "NoCollInRofStandard"};
    Configurable<bool> cfgNoHighMultCollInPrevRof{"cfgNoHighMultCollInPrevRof", true, "NoHighMultCollInPrevRof"};
    Configurable<bool> cfgNoITSROFrameBorder{"cfgNoITSROFrameBorder", true, "NoITSROFrameBorder"};
    Configurable<bool> cfgNoTimeFrameBorder{"cfgNoTimeFrameBorder", true, "NoTimeFrameBorder"};
    Configurable<bool> cfgTVXinTRD{"cfgTVXinTRD", true, "TVXinTRD - Use TVXinTRD (reject TRD triggered events)"};
    Configurable<bool> cfgIsVertexITSTPC{"cfgIsVertexITSTPC", true, "IsVertexITSTPC - Selects collisions with at least one ITS-TPC track"};
  } cfgEventCutFlags;
  struct : ConfigurableGroup {
    Configurable<int> cfgOccupancySelection{"cfgOccupancySelection", 2000, "Max occupancy selection, -999 to disable"};
    Configurable<bool> cfgDoOccupancySel{"cfgDoOccupancySel", true, "Bool for event selection on detector occupancy"};
    Configurable<float> cfgMagField{"cfgMagField", 99999, "Configurable magnetic field; default CCDB will be queried"};
    Configurable<bool> cfgMultCut{"cfgMultCut", false, "Use additional event cut on mult correlations"};
    Configurable<float> cfgVtxZ{"cfgVtxZ", 10, "vertex cut (cm)"};
  } cfgEventSelection;
  struct : ConfigurableGroup {
    Configurable<std::vector<double>> cfgMultGlobalCutPars{"cfgMultGlobalCutPars", std::vector<double>{2272.16, -76.6932, 1.01204, -0.00631545, 1.59868e-05, 136.336, -4.97006, 0.121199, -0.0015921, 7.66197e-06}, "Global vs FT0C multiplicity cut parameter values"};
    Configurable<std::vector<double>> cfgMultPVCutPars{"cfgMultPVCutPars", std::vector<double>{3074.43, -106.192, 1.46176, -0.00968364, 2.61923e-05, 182.128, -7.43492, 0.193901, -0.00256715, 1.22594e-05}, "PV vs FT0C multiplicity cut parameter values"};
    Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars", std::vector<double>{-0.223013, 0.715849, 0.664242, 0.0829653, -0.000503733, 1.21185e-06}, "Global vs PV multiplicity cut parameter values"};
    Configurable<std::string> cfgMultCorrHighCutFunction{"cfgMultCorrHighCutFunction", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut"};
    Configurable<std::string> cfgMultCorrLowCutFunction{"cfgMultCorrLowCutFunction", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut"};
    Configurable<std::string> cfgMultGlobalPVCorrCutFunction{"cfgMultGlobalPVCorrCutFunction", "[0] + [1]*x + 3*([2] + [3]*x + [4]*x*x + [5]*x*x*x)", "Functional for global vs pv multiplicity correlation cut"};
    Configurable<std::string> cfgMultGlobalASideCorrCutFunction{"cfgMultGlobalASideCorrCutFunction", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + [10]*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", "Functional for global vs V0A multiplicity low correlation cut"};
    Configurable<std::vector<double>> cfgMultGlobalV0ACutPars{"cfgMultGlobalV0ACutPars", std::vector<double>{567.785, 172.715, 0.77888, -0.00693466, 1.40564e-05, 679.853, 66.8068, -0.444332, 0.00115002, -4.92064e-07}, "Global vs FV0A multiplicity cut parameter values"};
    Configurable<std::vector<double>> cfgMultGlobalT0ACutPars{"cfgMultGlobalT0ACutPars", std::vector<double>{241.618, 61.8402, 0.348049, -0.00306078, 6.20357e-06, 315.235, 29.1491, -0.188639, 0.00044528, -9.08912e-08}, "Global vs FT0A multiplicity cut parameter values"};
    Configurable<float> cfgGlobalV0ALowSigma{"cfgGlobalV0ALowSigma", -3, "Number of sigma deviations below expected value in global vs V0A correlation"};
    Configurable<float> cfgGlobalV0AHighSigma{"cfgGlobalV0AHighSigma", 4, "Number of sigma deviations above expected value in global vs V0A correlation"};
    Configurable<float> cfgGlobalT0ALowSigma{"cfgGlobalT0ALowSigma", -3., "Number of sigma deviations below expected value in global vs T0A correlation"};
    Configurable<float> cfgGlobalT0AHighSigma{"cfgGlobalT0AHighSigma", 4, "Number of sigma deviations above expected value in global vs T0A correlation"};
  } cfgMultCorrCuts;

  Configurable<GFWBinningCuts> cfgGFWBinning{"cfgGFWBinning", {40, 16, 72, 300, 0, 3000, 0.2, 10.0, 0.2, 3.0, {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90}}, "Configuration for binning"};
  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN", "refP", "refFull"}, {-0.8, 0.4, -0.8}, {-0.4, 0.8, 0.8}, {0, 0, 0}, {1, 1, 1}}, "Configurations for GFW regions"};

  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}", "refFull {2 -2}", "refFull {2 2 -2 -2}"}, {"ChGap22", "ChGap32", "ChGap42", "ChFull22", "ChFull24"}, {0, 0, 0, 0, 0}, {1, 1, 1, 0, 0}}, "Configurations for each correlation to calculate"};

  struct GFWMemberCache {
    std::vector<double> ptbinning = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};
    float ptpoilow = 0.2;
    float ptpoiup = 10.0;
    float ptreflow = 0.2;
    float ptrefup = 3.0;
    float ptlow = 0.2;
    float ptup = 10.0;
    int etabins = 16;
    float etalow = -0.8;
    float etaup = 0.8;
    int vtxZbins = 40;
    int phibins = 72;
    float philow = 0.0;
    float phiup = o2::constants::math::TwoPI;
    int nchbins = 300;
    float nchlow = 0;
    float nchup = 3000;
    std::vector<double> centbinning = std::vector<double>(90);
    GFWRegions regions;
    GFWCorrConfigs configs;
    std::vector<std::pair<double, double>> etagapsPtPt;
    std::vector<double> multGlobalCorrCutPars;
    std::vector<double> multPVCorrCutPars;
    std::vector<double> multGlobalPVCorrCutPars;
    std::vector<double> multGlobalV0ACutPars;
    std::vector<double> multGlobalT0ACutPars;
  } gfwMemberCache;

  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb{};

  struct Config {
    TH1* mEfficiency = nullptr;
    std::vector<GFWWeights*> mAcceptance;
    bool correctionsLoaded = false;
  } correctionsConfig;

  struct MultipletConfig {
    std::size_t representativeConfigIndex = 0;
    int order = 0;
    std::shared_ptr<TProfile> profile;
  };

  std::map<std::string, MultipletConfig> effectiveMultiplets;

  // Outputs
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<FlowContainer> fFCgen{FlowContainer("FlowContainer_gen")};
  OutputObj<FlowPtContainer> fFCpt{FlowPtContainer("FlowPtContainer")};
  HistogramRegistry registry{"registry"};

  // Enums
  enum CentEstimators {
    CentFT0C = 0,
    CentFT0CVariant1,
    CentFT0M,
    CentFV0A,
    CentNTPV,
    CentNGlobal,
    CentMFT
  };
  std::map<int, std::string> centNamesMap = {{CentFT0C, "FT0C"}, {CentFT0CVariant1, "FT0C variant1"}, {CentFT0M, "FT0M"}, {CentFV0A, "FV0A"}, {CentNTPV, "NTPV"}, {CentNGlobal, "NGlobal"}, {CentMFT, "MFT"}};
  enum EventSelFlags {
    FilteredEvent = 1,
    Sel8,
    Occupancy,
    TVXinTRD,
    NoSameBunchPileup,
    IsGoodZvtxFT0vsPV,
    NoCollInTimeRangeStandard,
    NoCollInRofStandard,
    NoHighMultCollInPrevRof,
    NoTimeFrameBorder,
    NoITSROFrameBorder,
    IsVertexITSTPC,
    IsGoodITSLayersAll,
    MultCuts,
    TrackCent
  };
  struct EventCut {
    bool enabled;
    int histBin;
    int flag; // just store the enum
  };
  std::vector<EventCut> eventcutflags;
  enum ParticleIDs {
    ChargedID = 0,
    PionID,
    KaonID,
    ProtonID,
    SpeciesCount
  };

  // Generic Framework
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs{};
  TRandom3* fRndm = new TRandom3(0);
  TAxis* fPtAxis = nullptr;
  int lastRun = -1;
  std::vector<std::string> multipletKeys{};

  // Track selection - DCA functions
  TF1* fPtDepDCAxy = nullptr;
  TF1* fPtDepDCAz = nullptr;

  // Event selection cuts - multiplicity correlation
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultPVGlobalCutHigh = nullptr;
  TF1* fMultGlobalV0ACutLow = nullptr;
  TF1* fMultGlobalV0ACutHigh = nullptr;
  TF1* fMultGlobalT0ACutLow = nullptr;
  TF1* fMultGlobalT0ACutHigh = nullptr;

  // Track selection - pt-phi cuts
  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;

  void init(InitContext const&)
  {
    gfwMemberCache.regions.SetNames(cfgRegions->GetNames());
    gfwMemberCache.regions.SetEtaMin(cfgRegions->GetEtaMin());
    gfwMemberCache.regions.SetEtaMax(cfgRegions->GetEtaMax());
    gfwMemberCache.regions.SetpTDifs(cfgRegions->GetpTDifs());
    gfwMemberCache.regions.SetBitmasks(cfgRegions->GetBitmasks());
    gfwMemberCache.configs.SetCorrs(cfgCorrConfig->GetCorrs());
    gfwMemberCache.configs.SetHeads(cfgCorrConfig->GetHeads());
    gfwMemberCache.configs.SetpTDifs(cfgCorrConfig->GetpTDifs());
    gfwMemberCache.configs.SetpTCorrMasks(cfgCorrConfig->GetpTCorrMasks());
    gfwMemberCache.regions.Print();
    gfwMemberCache.configs.Print();
    gfwMemberCache.ptbinning = cfgGFWBinning->GetPtBinning();
    gfwMemberCache.ptpoilow = cfgGFWBinning->GetPtPOImin();
    gfwMemberCache.ptpoiup = cfgGFWBinning->GetPtPOImax();
    gfwMemberCache.ptreflow = cfgGFWBinning->GetPtRefMin();
    gfwMemberCache.ptrefup = cfgGFWBinning->GetPtRefMax();
    gfwMemberCache.ptlow = cfgKinematics.cfgPtCut->first;
    gfwMemberCache.ptup = cfgKinematics.cfgPtCut->second;
    gfwMemberCache.etabins = cfgGFWBinning->GetEtaBins();
    gfwMemberCache.vtxZbins = cfgGFWBinning->GetVtxZbins();
    gfwMemberCache.phibins = cfgGFWBinning->GetPhiBins();
    gfwMemberCache.philow = 0.0f;
    gfwMemberCache.phiup = o2::constants::math::TwoPI;
    gfwMemberCache.nchbins = cfgGFWBinning->GetNchBins();
    gfwMemberCache.nchlow = cfgGFWBinning->GetNchMin();
    gfwMemberCache.nchup = cfgGFWBinning->GetNchMax();
    gfwMemberCache.centbinning = cfgGFWBinning->GetCentBinning();
    cfgGFWBinning->Print();
    gfwMemberCache.multGlobalCorrCutPars = cfgMultCorrCuts.cfgMultGlobalCutPars;
    gfwMemberCache.multPVCorrCutPars = cfgMultCorrCuts.cfgMultPVCutPars;
    gfwMemberCache.multGlobalPVCorrCutPars = cfgMultCorrCuts.cfgMultGlobalPVCutPars;
    gfwMemberCache.multGlobalV0ACutPars = cfgMultCorrCuts.cfgMultGlobalV0ACutPars;
    gfwMemberCache.multGlobalT0ACutPars = cfgMultCorrCuts.cfgMultGlobalT0ACutPars;

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    // Setup event cuts
    eventcutflags.push_back({cfgEventCutFlags.cfgNoSameBunchPileupCut, NoSameBunchPileup, o2::aod::evsel::kNoSameBunchPileup});
    eventcutflags.push_back({cfgEventCutFlags.cfgIsGoodZvtxFT0vsPV, IsGoodZvtxFT0vsPV, o2::aod::evsel::kIsGoodZvtxFT0vsPV});
    eventcutflags.push_back({cfgEventCutFlags.cfgNoCollInTimeRangeStandard, NoCollInTimeRangeStandard, o2::aod::evsel::kNoCollInTimeRangeStandard});
    eventcutflags.push_back({cfgEventCutFlags.cfgNoCollInRofStandard, NoCollInRofStandard, o2::aod::evsel::kNoCollInRofStandard});
    eventcutflags.push_back({cfgEventCutFlags.cfgNoHighMultCollInPrevRof, NoHighMultCollInPrevRof, o2::aod::evsel::kNoHighMultCollInPrevRof});
    eventcutflags.push_back({cfgEventCutFlags.cfgNoTimeFrameBorder, NoTimeFrameBorder, o2::aod::evsel::kNoTimeFrameBorder});
    eventcutflags.push_back({cfgEventCutFlags.cfgNoITSROFrameBorder, NoITSROFrameBorder, o2::aod::evsel::kNoITSROFrameBorder});
    eventcutflags.push_back({cfgEventCutFlags.cfgIsVertexITSTPC, IsVertexITSTPC, o2::aod::evsel::kIsVertexITSTPC});
    eventcutflags.push_back({cfgEventCutFlags.cfgIsGoodITSLayersAll, IsGoodITSLayersAll, o2::aod::evsel::kIsGoodITSLayersAll});
    for (const auto& cut : eventcutflags) {
      LOGF(info, "Flag %d is %senabled", cut.histBin, (cut.enabled) ? "" : "not ");
    }

    AxisSpec phiAxis = {gfwMemberCache.phibins, gfwMemberCache.philow, gfwMemberCache.phiup, "#phi"};
    AxisSpec phiModAxis = {100, 0, constants::math::PI / 9, "fmod(#varphi,#pi/9)"};
    AxisSpec etaAxis = {gfwMemberCache.etabins, cfgKinematics.cfgEta->first, cfgKinematics.cfgEta->second, "#eta"};
    AxisSpec vtxAxis = {gfwMemberCache.vtxZbins, -cfgEventSelection.cfgVtxZ, cfgEventSelection.cfgVtxZ, "Vtx_{z} (cm)"};
    AxisSpec ptAxis = {gfwMemberCache.ptbinning, "#it{p}_{T} GeV/#it{c}"};
    std::string sCentralityEstimator = centNamesMap[cfgCentEstimator] + " centrality (%)";
    AxisSpec centAxis = {gfwMemberCache.centbinning, sCentralityEstimator.c_str()};
    std::vector<double> nchbinning;
    const double nchskip = (gfwMemberCache.nchup - gfwMemberCache.nchlow) / gfwMemberCache.nchbins;
    for (int i = 0; i <= gfwMemberCache.nchbins; ++i) {
      nchbinning.push_back(nchskip * i + gfwMemberCache.nchlow + 0.5);
    }
    AxisSpec nchAxis = {nchbinning, "N_{ch}"};
    AxisSpec bAxis = {200, 0, 20, "#it{b}"};
    AxisSpec t0cAxis = {1000, 0, 50000, "N_{ch} (T0C)"};
    AxisSpec t0aAxis = {1800, 0, 180000, "N_{ch} (T0A)"};
    AxisSpec v0aAxis = {1800, 0, 180000, "N_{ch} (V0A)"};
    AxisSpec multpvAxis = {3500, 0, 3500, "N_{ch} (PV)"};
    AxisSpec occAxis = {500, 0, 5000, "occupancy"};
    AxisSpec multAxis = (cfgUseNch) ? nchAxis : centAxis;
    AxisSpec dcaZAXis = {200, -2, 2, "DCA_{z} (cm)"};
    AxisSpec dcaXYAXis = {200, -1, 1, "DCA_{xy} (cm)"};

    if (cfgFillQA) {
      registry.add("trackQA/before/phi_eta_vtxZ", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registry.add("trackQA/before/pt_dcaXY_dcaZ", "", {HistType::kTH3D, {ptAxis, dcaXYAXis, dcaZAXis}});
      registry.add("trackQA/before/pt_phi", "", {HistType::kTH2D, {ptAxis, phiModAxis}});
      registry.add("trackQA/before/chi2prTPCcls", "#chi^{2}/cluster for the TPC track segment", {HistType::kTH1D, {{100, 0., 5.}}});
      registry.add("trackQA/before/chi2prITScls", "#chi^{2}/cluster for the ITS track", {HistType::kTH1D, {{100, 0., 50.}}});
      registry.add("trackQA/before/nTPCClusters", "Number of found TPC clusters", {HistType::kTH1D, {{100, 40, 180}}});
      registry.add("trackQA/before/nITSClusters", "Number of found ITS clusters", {HistType::kTH1D, {{100, 0, 20}}});
      registry.add("trackQA/before/nTPCCrossedRows", "Number of crossed TPC Rows", {HistType::kTH1D, {{100, 40, 180}}});
      registry.addClone("trackQA/before/", "trackQA/after/");
      registry.add("trackQA/after/pt_ref", "; #it{p}_{T}; Counts", {HistType::kTH1D, {{100, gfwMemberCache.ptreflow, gfwMemberCache.ptrefup}}});
      registry.add("trackQA/after/pt_poi", "; #it{p}_{T}; Counts", {HistType::kTH1D, {{100, gfwMemberCache.ptpoilow, gfwMemberCache.ptpoiup}}});
      registry.add("trackQA/after/Nch_corrected", "; N_{ch}; Counts", {HistType::kTH1D, {nchAxis}});
      registry.add("trackQA/after/Nch_uncorrected", "; N_{ch}; Counts", {HistType::kTH1D, {nchAxis}});
      registry.add("trackQA/after/etaNch", "; #eta; Counts", {HistType::kTH1D, {etaAxis}});
      registry.add("trackQA/after/etaPtPt", "; #eta; Counts", {HistType::kTH1D, {etaAxis}});

      registry.add("eventQA/before/globalTracks_centT0C", "; FT0C centrality (%); N_{global}", {HistType::kTH2D, {centAxis, nchAxis}});
      registry.add("eventQA/before/PVTracks_centT0C", "; FT0C centrality (%); N_{PV}", {HistType::kTH2D, {centAxis, multpvAxis}});
      registry.add("eventQA/before/globalTracks_PVTracks", "; N_{PV}; N_{global}", {HistType::kTH2D, {multpvAxis, nchAxis}});
      registry.add("eventQA/before/globalTracks_multT0A", "; multT0A; N_{global}", {HistType::kTH2D, {t0aAxis, nchAxis}});
      registry.add("eventQA/before/globalTracks_multV0A", "; multV0A; N_{global}", {HistType::kTH2D, {v0aAxis, nchAxis}});
      registry.add("eventQA/before/multV0A_multT0A", "; multV0A; multT0A", {HistType::kTH2D, {t0aAxis, v0aAxis}});
      registry.add("eventQA/before/multT0C_centT0C", "; multT0C; FT0C centrality (%)", {HistType::kTH2D, {centAxis, t0cAxis}});
      registry.add("eventQA/before/occ_mult_cent", "; occupancy; N_{ch}; centrality (%)", {HistType::kTH3D, {occAxis, nchAxis, centAxis}});
    }
    registry.add("eventQA/before/centrality", "; centrality (%); Counts", {HistType::kTH1D, {centAxis}});
    registry.add("eventQA/before/multiplicity", "; N_{ch}; Counts", {HistType::kTH1D, {nchAxis}});
    registry.addClone("eventQA/before/", "eventQA/after/");
    registry.add("eventQA/eventSel", "Number of Events;; Counts", {HistType::kTH1D, {{15, 0.5, 15.5}}});
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(FilteredEvent, "Filtered event");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(Sel8, "sel8");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(Occupancy, "occupancy");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(TVXinTRD, "TVXinTRD");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(NoSameBunchPileup, "NoSameBunchPileup");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(IsGoodZvtxFT0vsPV, "IsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(NoCollInTimeRangeStandard, "NoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(NoCollInRofStandard, "NoCollInRofStandard");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(NoHighMultCollInPrevRof, "NoHighMultCollInPrevRof");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(NoTimeFrameBorder, "NoTimeFrameBorder");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(NoITSROFrameBorder, "NoITSROFrameBorder");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(IsVertexITSTPC, "IsVertexITSTPC");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(IsGoodITSLayersAll, "IsGoodITSLayersAll");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(MultCuts, "after Mult cuts");
    registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(TrackCent, "has track + within cent");

    const int ptbins = static_cast<int>(gfwMemberCache.ptbinning.size() - 1);
    fPtAxis = new TAxis(ptbins, gfwMemberCache.ptbinning.data());

    if (gfwMemberCache.regions.GetSize() < 0) {
      LOGF(error, "Configuration contains vectors of different size - check the GFWRegions configurable");
    }
    for (auto i(0); i < gfwMemberCache.regions.GetSize(); ++i) {
      fGFW->AddRegion(gfwMemberCache.regions.GetNames()[i], gfwMemberCache.regions.GetEtaMin()[i], gfwMemberCache.regions.GetEtaMax()[i], (gfwMemberCache.regions.GetpTDifs()[i] != 0) ? ptbins + 1 : 1, gfwMemberCache.regions.GetBitmasks()[i]);
    }
    for (auto i = 0; i < gfwMemberCache.configs.GetSize(); ++i) {
      corrconfigs.push_back(fGFW->GetCorrelatorConfig(gfwMemberCache.configs.GetCorrs()[i], gfwMemberCache.configs.GetHeads()[i], gfwMemberCache.configs.GetpTDifs()[i] != 0));
    }
    if (corrconfigs.empty()) {
      LOGF(error, "Configuration contains vectors of different size - check the GFWCorrConfig configurable");
    }
    fGFW->CreateRegions();
    auto oba = new TObjArray();
    addConfigObjectsToObjArray(oba, corrconfigs);
    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fPtAxis);
    fFC->Initialize(oba, multAxis, cfgNbootstrap);
    fFCgen->SetName("FlowContainer_gen");
    fFCgen->SetXAxis(fPtAxis);
    fFCgen->Initialize(oba, multAxis, cfgNbootstrap);
    delete oba;

    // Identify unique multiplet setups and populate with identifying keys
    identifyUniqueMultipletKeys(multipletKeys, corrconfigs, centAxis); // For now keep centrality axis hardcoded
    LOGF(info, "Unique multiplet keys for current configuration setup");
    for (const auto& key : multipletKeys) {
      LOGF(info, key);
    }

    const auto& ptPtGaps = cfgPtPtGaps->getData();
    for (uint32_t i = 0; i < ptPtGaps.rows; ++i) {
      const auto etaMin = ptPtGaps(i, 0);
      const auto etaMax = ptPtGaps(i, 1);

      if (etaMin < -1. || etaMax < -1.) {
        continue;
      }

      gfwMemberCache.etagapsPtPt.emplace_back(etaMin, etaMax);
    }
    for (const auto& [etamin, etamax] : gfwMemberCache.etagapsPtPt) {
      LOGF(info, "pt-pt subevent: {%.1f,%.1f}", etamin, etamax);
    }

    fFCpt->setUseCentralMoments(cfgUseCentralMoments);
    fFCpt->setUseGapMethod(true);
    fFCpt->initialise(multAxis, cfgMpar, gfwMemberCache.configs, cfgNbootstrap);
    fFCpt->initialiseSubevent(multAxis, cfgMpar, gfwMemberCache.etagapsPtPt.size(), cfgNbootstrap);

    fPtDepDCAxy = new TF1("ptDepDCAxy", Form("[0]*%s", cfgTrackCuts.cfgDCAxyPtDep->c_str()), 0.001, 100);
    fPtDepDCAxy->SetParameter(0, cfgTrackCuts.cfgDCAxyNSigma / 7.);
    LOGF(info, "DCAxy pt-dependence function: %s", Form("[0]*%s", cfgTrackCuts.cfgDCAxyPtDep->c_str()));
    if (!cfgTrackCuts.cfgDCAzPtDep.value.empty()) {
      fPtDepDCAz = new TF1("ptDepDCAz", Form("%s", cfgTrackCuts.cfgDCAzPtDep->c_str()), 0.001, 100);
      LOGF(info, "DCAz pt-dependence function: %s", Form("%s", cfgTrackCuts.cfgDCAzPtDep->c_str()));
    }

    // Multiplicity correlation cuts
    if (cfgEventSelection.cfgMultCut) {
      fMultPVCutLow = new TF1("fMultPVCutLow", cfgMultCorrCuts.cfgMultCorrLowCutFunction->c_str(), 0, 100);
      fMultPVCutLow->SetParameters(gfwMemberCache.multPVCorrCutPars.data());
      fMultPVCutHigh = new TF1("fMultPVCutHigh", cfgMultCorrCuts.cfgMultCorrHighCutFunction->c_str(), 0, 100);
      fMultPVCutHigh->SetParameters(gfwMemberCache.multPVCorrCutPars.data());
      fMultCutLow = new TF1("fMultCutLow", cfgMultCorrCuts.cfgMultCorrLowCutFunction->c_str(), 0, 100);
      fMultCutLow->SetParameters(gfwMemberCache.multGlobalCorrCutPars.data());
      fMultCutHigh = new TF1("fMultCutHigh", cfgMultCorrCuts.cfgMultCorrHighCutFunction->c_str(), 0, 100);
      fMultCutHigh->SetParameters(gfwMemberCache.multGlobalCorrCutPars.data());
      fMultPVGlobalCutHigh = new TF1("fMultPVGlobalCutHigh", cfgMultCorrCuts.cfgMultGlobalPVCorrCutFunction->c_str(), 0, nchbinning.back());
      fMultPVGlobalCutHigh->SetParameters(gfwMemberCache.multGlobalPVCorrCutPars.data());

      LOGF(info, "Global V0A function: %s in range 0-%g", cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), v0aAxis.binEdges.back());
      fMultGlobalV0ACutLow = new TF1("fMultGlobalV0ACutLow", cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, v0aAxis.binEdges.back());
      for (std::size_t i = 0; i < gfwMemberCache.multGlobalV0ACutPars.size(); ++i) {
        fMultGlobalV0ACutLow->SetParameter(i, gfwMemberCache.multGlobalV0ACutPars[i]);
      }
      fMultGlobalV0ACutLow->SetParameter(gfwMemberCache.multGlobalV0ACutPars.size(), cfgMultCorrCuts.cfgGlobalV0ALowSigma);
      for (int i = 0; i < fMultGlobalV0ACutLow->GetNpar(); ++i) {
        LOGF(info, "fMultGlobalV0ACutLow par %d = %g", i, fMultGlobalV0ACutLow->GetParameter(i));
      }

      fMultGlobalV0ACutHigh = new TF1("fMultGlobalV0ACutHigh", cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, v0aAxis.binEdges.back());
      for (std::size_t i = 0; i < gfwMemberCache.multGlobalV0ACutPars.size(); ++i) {
        fMultGlobalV0ACutHigh->SetParameter(i, gfwMemberCache.multGlobalV0ACutPars[i]);
      }
      fMultGlobalV0ACutHigh->SetParameter(gfwMemberCache.multGlobalV0ACutPars.size(), cfgMultCorrCuts.cfgGlobalV0AHighSigma);
      for (int i = 0; i < fMultGlobalV0ACutHigh->GetNpar(); ++i) {
        LOGF(info, "fMultGlobalV0ACutHigh par %d = %g", i, fMultGlobalV0ACutHigh->GetParameter(i));
      }

      LOGF(info, "Global T0A function: %s", cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str());
      fMultGlobalT0ACutLow = new TF1("fMultGlobalT0ACutLow", cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, t0aAxis.binEdges.back());
      for (std::size_t i = 0; i < gfwMemberCache.multGlobalT0ACutPars.size(); ++i) {
        fMultGlobalT0ACutLow->SetParameter(i, gfwMemberCache.multGlobalT0ACutPars[i]);
      }
      fMultGlobalT0ACutLow->SetParameter(gfwMemberCache.multGlobalT0ACutPars.size(), cfgMultCorrCuts.cfgGlobalT0ALowSigma);
      for (int i = 0; i < fMultGlobalT0ACutLow->GetNpar(); ++i) {
        LOGF(info, "fMultGlobalT0ACutLow par %d = %g", i, fMultGlobalT0ACutLow->GetParameter(i));
      }

      fMultGlobalT0ACutHigh = new TF1("fMultGlobalT0ACutHigh", cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, t0aAxis.binEdges.back());
      for (std::size_t i = 0; i < gfwMemberCache.multGlobalT0ACutPars.size(); ++i) {
        fMultGlobalT0ACutHigh->SetParameter(i, gfwMemberCache.multGlobalT0ACutPars[i]);
      }
      fMultGlobalT0ACutHigh->SetParameter(gfwMemberCache.multGlobalT0ACutPars.size(), cfgMultCorrCuts.cfgGlobalT0AHighSigma);
      for (int i = 0; i < fMultGlobalT0ACutHigh->GetNpar(); ++i) {
        LOGF(info, "fMultGlobalT0ACutHigh par %d = %g", i, fMultGlobalT0ACutHigh->GetParameter(i));
      }
    }
  }

  static constexpr std::array<std::string_view, 2> FillTimeName = {"before/", "after/"};

  void addConfigObjectsToObjArray(TObjArray* oba, const std::vector<GFW::CorrConfig>& configs)
  {
    for (auto it = configs.begin(); it != configs.end(); ++it) {
      if (it->pTDif) {
        std::string suffix = "_ptDiff";
        for (auto i = 0; i < fPtAxis->GetNbins(); ++i) {
          std::string index = Form("_pt_%i", i + 1);
          oba->Add(new TNamed(it->Head + index, it->Head + suffix));
        }
      } else {
        oba->Add(new TNamed(it->Head.c_str(), it->Head.c_str()));
      }
    }
  }

  void identifyUniqueMultipletKeys(std::vector<std::string>& keys, const std::vector<GFW::CorrConfig>& configs, AxisSpec xAxis)
  {
    std::set<std::string> uniqueKeys(keys.begin(), keys.end());
    const auto& regionNames = gfwMemberCache.regions.GetNames();

    for (std::size_t configIndex = 0; configIndex < configs.size(); ++configIndex) {
      const auto& config = configs[configIndex];
      // TODO: Support pT-differential and explicit-overlap configurations.
      if (config.pTDif) {
        continue;
      }

      if (config.Regs.size() != config.Hars.size()) {
        LOGF(error, "Cannot construct multiplet key for %s: %zu region groups but %zu harmonic groups", config.Head.c_str(), config.Regs.size(), config.Hars.size());
        continue;
      }

      std::string key;
      int totalOrder = 0;
      bool valid = true;

      for (std::size_t group = 0; group < config.Regs.size(); ++group) {
        const auto& regionIndices = config.Regs[group];
        const auto& harmonics = config.Hars[group];

        if (regionIndices.size() != 1 || harmonics.empty()) {
          LOGF(error, "Cannot construct simple multiplet key for %s, group %zu: expected one region and at least one harmonic", config.Head.c_str(), group);
          valid = false;
          break;
        }

        const int regionIndex = regionIndices.front();
        if (regionIndex < 0 ||
            static_cast<std::size_t>(regionIndex) >= regionNames.size()) {
          LOGF(error, "Invalid region index %d in configuration %s", regionIndex, config.Head.c_str());
          valid = false;
          break;
        }

        const int regionOrder = static_cast<int>(harmonics.size());
        totalOrder += regionOrder;

        if (!key.empty()) {
          key += "_";
        }

        key += regionNames[regionIndex];
        key += "_";
        key += std::to_string(regionOrder);
      }

      if (!valid || key.empty()) {
        continue;
      }

      key += "_";
      key += std::to_string(totalOrder);

      if (uniqueKeys.insert(key).second) {
        keys.push_back(std::move(key));
        MultipletConfig currentMultipletConfig{configIndex, totalOrder, registry.add<TProfile>(Form("%s", keys.back().c_str()), Form("; centrality ; N_{%d}", totalOrder), {HistType::kTProfile, {xAxis}})};
        effectiveMultiplets[keys.back()] = currentMultipletConfig;
      }
    }
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    // static o2::parameters::GRPObject* grpo = nullptr;
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      // grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  void loadCorrections(aod::BCsWithTimestamps::iterator const& bc)
  {
    uint64_t timestamp = bc.timestamp();
    if (!cfgRunByRun && correctionsConfig.correctionsLoaded) {
      return;
    }
    if (!cfgCorrections.cfgAcceptancePath.value.empty()) {
      std::string runstr = (cfgRunByRun) ? "RunByRun/" : "";
      correctionsConfig.mAcceptance.clear();
      correctionsConfig.mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgCorrections.cfgAcceptancePath.value + runstr, timestamp));
    }
    // Run-by-run efficiencies are not supported at the moment
    if (correctionsConfig.correctionsLoaded) {
      return;
    }
    if (!cfgCorrections.cfgEfficiencyPath.value.empty()) {
      if (cfgCorrections.cfgUse2DEfficiency) {
        correctionsConfig.mEfficiency = ccdb->getForTimeStamp<TH2D>(cfgCorrections.cfgEfficiencyPath, timestamp);
      } else {
        correctionsConfig.mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgCorrections.cfgEfficiencyPath, timestamp);
      }
      if (correctionsConfig.mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram from %s", cfgCorrections.cfgEfficiencyPath.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s", cfgCorrections.cfgEfficiencyPath.value.c_str());
    }
    correctionsConfig.correctionsLoaded = true;
  }

  template <typename TTrack>
  double getAcceptance(const TTrack& track, const double& vtxz)
  { // 0 ref, 1 ch, 2 pi, 3 ka, 4 pr
    double wacc = 1;
    if (!correctionsConfig.mAcceptance.empty()) {
      wacc = correctionsConfig.mAcceptance[0]->getNUA(track.phi(), track.eta(), vtxz);
    }
    return wacc;
  }

  template <typename TTrack>
  double getEfficiency(const TTrack& track, const float& centrality)
  { //-1 ref, 0 ch, 1 pi, 2 ka, 3 pr, 4 k0, 5 lambda
    double eff = 1.;
    if (!correctionsConfig.mEfficiency) {
      return eff;
    }
    if (cfgCorrections.cfgUse2DEfficiency) {
      eff = dynamic_cast<TH2D*>(correctionsConfig.mEfficiency)->GetBinContent(dynamic_cast<TH2D*>(correctionsConfig.mEfficiency)->FindBin(track.pt(), centrality));
    } else {
      eff = dynamic_cast<TH1D*>(correctionsConfig.mEfficiency)->GetBinContent(dynamic_cast<TH1D*>(correctionsConfig.mEfficiency)->FindBin(track.pt()));
    }
    if (eff == 0) {
      return -1.;
    }
    return 1. / eff;
  }

  template <typename TCollision>
  bool eventSelected(const TCollision& collision, const int multTrk, const float& centrality)
  {
    // Cut on trigger alias
    if (cfgEventCutFlags.cfgTVXinTRD) {
      if (collision.alias_bit(TVXinTRD)) {
        // TRD triggered
        // "CMTVX-B-NOPF-TRD,minbias_TVX"
        return false;
      }
      registry.fill(HIST("eventQA/eventSel"), TVXinTRD);
    }
    // Cut on event selection flags
    for (const auto& cut : eventcutflags) {
      if (!cut.enabled) {
        continue;
      }
      if (!collision.selection_bit(cut.flag)) {
        return false;
      }
      registry.fill(HIST("eventQA/eventSel"), cut.histBin);
    }
    // Cut on vertex
    if (!selectVertex(collision)) {
      return false;
    }
    // Cut on multiplicity correlations - data driven
    if (cfgEventSelection.cfgMultCut) {
      if (!selectMultiplicityCorrelation(collision, multTrk, centrality)) {
        return false;
      }
    }
    return true;
  }

  template <typename TCollision>
  bool selectVertex(const TCollision& collision)
  {
    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = std::sqrt(collision.covZZ());
      float minZRes = 0.25;
      int minNContrib = 20;
      if (zRes > minZRes && collision.numContrib() < minNContrib) {
        vtxz = -999;
      }
    }
    if (vtxz > cfgEventSelection.cfgVtxZ || vtxz < -cfgEventSelection.cfgVtxZ) {
      return false;
    }

    return true;
  }

  template <typename TCollision>
  bool selectMultiplicityCorrelation(const TCollision& collision, const int multTrk, const float& centrality)
  {
    auto multNTracksPV = collision.multNTracksPV();
    if (multNTracksPV < fMultPVCutLow->Eval(centrality)) {
      return false;
    }
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality)) {
      return false;
    }
    if (multTrk < fMultCutLow->Eval(centrality)) {
      return false;
    }
    if (multTrk > fMultCutHigh->Eval(centrality)) {
      return false;
    }
    if (multTrk > fMultPVGlobalCutHigh->Eval(collision.multNTracksPV())) {
      return false;
    }

    if (!(cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->empty()) && static_cast<double>(collision.multFV0A()) < fMultGlobalV0ACutLow->Eval(multTrk)) {
      return false;
    }
    if (!(cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->empty()) && static_cast<double>(collision.multFV0A()) > fMultGlobalV0ACutHigh->Eval(multTrk)) {
      return false;
    }
    if (!(cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->empty()) && static_cast<double>(collision.multFT0A()) < fMultGlobalT0ACutLow->Eval(multTrk)) {
      return false;
    }
    if (!(cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->empty()) && static_cast<double>(collision.multFT0A()) > fMultGlobalT0ACutHigh->Eval(multTrk)) {
      return false;
    }
    registry.fill(HIST("eventQA/eventSel"), MultCuts);
    return true;
  }

  template <typename TTrack>
  bool trackSelected(const TTrack& track, const int field)
  {
    if (cfgTrackCuts.cfgTPCSectorCut) {
      double phimodn = track.phi();
      if (field < 0) { // for negative polarity field
        phimodn = o2::constants::math::TwoPI - phimodn;
      }
      if (track.sign() < 0) { // for negative charge
        phimodn = o2::constants::math::TwoPI - phimodn;
      }
      if (phimodn < 0) {
        LOGF(warning, "phi < 0: %g", phimodn);
      }

      phimodn += o2::constants::math::PI / 18.0; // to center gap in the middle
      phimodn = fmod(phimodn, o2::constants::math::PI / 9.0);
      if (cfgFillQA) {
        registry.fill(HIST("trackQA/before/pt_phi"), track.pt(), phimodn);
      }
      if (phimodn < fPhiCutHigh->Eval(track.pt()) && phimodn > fPhiCutLow->Eval(track.pt())) {
        return false; // reject track
      }
      if (cfgFillQA) {
        registry.fill(HIST("trackQA/after/pt_phi"), track.pt(), phimodn);
      }
    }
    if (cfgTrackCuts.cfgDCAxyNSigma && (std::fabs(track.dcaXY()) > fPtDepDCAxy->Eval(track.pt()))) {
      return false;
    }
    if (!cfgTrackCuts.cfgDCAzPtDep.value.empty() && std::fabs(track.dcaZ() > fPtDepDCAz->Eval(track.pt()))) {
      return false;
    }
    return ((track.tpcNClsCrossedRows() >= cfgTrackCuts.cfgNTPCXrows) && (track.tpcNClsFound() >= cfgTrackCuts.cfgNTPCCls) && (track.itsNCls() >= cfgTrackCuts.cfgMinNITSCls));
  }

  template <typename TTrack>
  bool nchSelected(const TTrack& track)
  {
    // Renormalise to default cut
    const float defaultNsigma = 7;
    if (cfgTrackCuts.cfgDCAxyNSigma && (std::fabs(track.dcaXY()) > defaultNsigma / cfgTrackCuts.cfgDCAxyNSigma * fPtDepDCAxy->Eval(track.pt()))) {
      return false;
    }
    if (!cfgTrackCuts.cfgDCAzPtDep.value.empty() && std::fabs(track.dcaZ() > fPtDepDCAz->Eval(track.pt()))) {
      return false;
    }
    int tpcNClsCrossedRowsDefault = 70;
    int tpcNClsFoundDefault = 50;
    int itsNclsDefault = 5;
    return ((track.tpcNClsCrossedRows() >= tpcNClsCrossedRowsDefault) && (track.tpcNClsFound() >= tpcNClsFoundDefault) && (track.itsNCls() >= itsNclsDefault));
  }

  enum DataType {
    Reco,
    Gen
  };

  template <DataType dt>
  void fillOutputContainers(const float& centmult, const double& rndm, const float& multiplicity)
  {
    fFCpt->calculateCorrelations();
    fFCpt->calculateSubeventCorrelations();
    fFCpt->fillPtProfiles(centmult, rndm);
    fFCpt->fillSubeventPtProfiles(centmult, rndm);
    fFCpt->fillCMProfiles(centmult, rndm);
    fFCpt->fillCMSubeventProfiles(centmult, rndm);

    for (const auto& [key, multiplet] : effectiveMultiplets) {
      const auto& config = corrconfigs.at(multiplet.representativeConfigIndex);
      const double dnx = fGFW->Calculate(config, 0, true).real();

      if (dnx > 0.) {
        const double multPower = std::pow(multiplicity, multiplet.order - 1);
        if (multPower > 0.) {
          multiplet.profile->Fill(centmult, 1. / multPower, dnx);
        }
      }
    }

    for (std::size_t l_ind = 0; l_ind < corrconfigs.size(); ++l_ind) {
      if (!corrconfigs.at(l_ind).pTDif) {
        auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), 0, true).real();
        if (dnx == 0) {
          continue;
        }
        auto val = fGFW->Calculate(corrconfigs.at(l_ind), 0, false).real() / dnx;
        if (std::abs(val) < 1) {
          fFC->FillProfile(corrconfigs.at(l_ind).Head.c_str(), centmult, val, cfgUseMultiplicityFlowWeights ? dnx : 1.0, rndm);
          fFCpt->fillVnPtProfiles(centmult, val, dnx, rndm, gfwMemberCache.configs.GetpTCorrMasks()[l_ind]);
        }
        continue;
      }
      for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
        auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), i - 1, true).real();
        if (dnx == 0) {
          continue;
        }
        auto val = fGFW->Calculate(corrconfigs.at(l_ind), i - 1, false).real() / dnx;
        if (std::abs(val) < 1) {
          fFC->FillProfile(Form("%s_pt_%i", corrconfigs.at(l_ind).Head.c_str(), i), centmult, val, cfgUseMultiplicityFlowWeights ? dnx : 1.0, rndm);
        }
      }
    }
  }

  template <DataType dt, typename TTrack>
  inline void fillPtSums(const TTrack& track, const float& centrality)
  {
    double weff = (dt == Gen) ? 1. : getEfficiency(track, centrality);
    if (weff < 0) {
      return;
    }
    // Fill the nominal sums
    if (track.eta() > cfgKinematics.cfgEtaPtPt->first && track.eta() < cfgKinematics.cfgEtaPtPt->second) {
      fFCpt->fill(weff, track.pt());
    }
    // Fill the subevent sums
    std::size_t index = 0;
    for (const auto& [etamin, etamax] : gfwMemberCache.etagapsPtPt) {
      if (etamin < track.eta() && track.eta() < etamax) {
        fFCpt->fillSub(weff, track.pt(), index);
      }
      ++index;
    }
  }

  template <DataType dt, typename TTrack>
  inline void fillGFW(const TTrack& track, const float& centrality, const double& vtxz)
  {
    bool withinPtRef = (track.pt() > gfwMemberCache.ptreflow && track.pt() < gfwMemberCache.ptrefup);
    bool withinPtPOI = (track.pt() > gfwMemberCache.ptpoilow && track.pt() < gfwMemberCache.ptpoiup);
    if (!withinPtPOI && !withinPtRef) {
      return;
    }
    double weff = (dt == Gen) ? 1. : getEfficiency(track, centrality);
    if (weff < 0) {
      return;
    }
    double wacc = (dt == Gen) ? 1. : getAcceptance(track, vtxz);
    if (withinPtRef) {
      fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 1);
    }
    if (withinPtPOI) {
      fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 2);
    }
    if (withinPtRef && withinPtPOI) {
      fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 4);
    }
  }

  template <typename TCollision>
  float getCentrality(const TCollision& collision)
  {
    switch (cfgCentEstimator) {
      case CentFT0C:
        return collision.centFT0C();
      case CentFT0CVariant1:
        return collision.centFT0CVariant1();
      case CentFT0M:
        return collision.centFT0M();
      case CentFV0A:
        return collision.centFV0A();
      case CentNTPV:
        return collision.centNTPV();
      case CentNGlobal:
        return collision.centNGlobal();
      case CentMFT:
        return collision.centMFT();
      default:
        return collision.centFT0C();
    }
  }

  enum QAFillTime {
    Before,
    After
  };

  template <DataType dt, QAFillTime ft, typename TTrack>
  inline void fillTrackQA(const TTrack& track, const float vtxz)
  {
    if constexpr (dt == Gen) {
      registry.fill(HIST("MCGen/") + HIST(FillTimeName[ft]) + HIST("phi_eta_vtxZ_gen"), track.phi(), track.eta(), vtxz);
      registry.fill(HIST("MCGen/") + HIST(FillTimeName[ft]) + HIST("pt_gen"), track.pt());
    } else {
      double wacc = getAcceptance(track, vtxz);
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("phi_eta_vtxZ"), track.phi(), track.eta(), vtxz, (ft == After) ? wacc : 1.0);
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("pt_dcaXY_dcaZ"), track.pt(), track.dcaXY(), track.dcaZ());

      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("chi2prTPCcls"), track.tpcChi2NCl());
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("chi2prITScls"), track.itsChi2NCl());
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("nTPCClusters"), track.tpcNClsFound());
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("nITSClusters"), track.itsNCls());
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("nTPCCrossedRows"), track.tpcNClsCrossedRows());

      if (ft == After) {
        registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("pt_ref"), track.pt());
        registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("pt_poi"), track.pt());

        if (track.eta() > cfgKinematics.cfgEtaNch->first && track.eta() < cfgKinematics.cfgEtaNch->second) {
          registry.fill(HIST("trackQA/after/etaNch"), track.eta());
        }
        if (track.eta() > cfgKinematics.cfgEtaPtPt->first && track.eta() < cfgKinematics.cfgEtaPtPt->second && cfgFillQA) {
          registry.fill(HIST("trackQA/after/etaPtPt"), track.eta());
        }
      }
    }
  }

  template <QAFillTime ft, typename TCollision, typename TTracks>
  inline void fillEventQA(const TCollision& collision, const TTracks& tracks)
  {
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_centT0C"), collision.centFT0C(), tracks.size());
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_multT0A"), collision.multFT0A(), tracks.size());
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_multV0A"), collision.multFV0A(), tracks.size());
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
  }

  struct AcceptedTracks {
    float total = 0;
    unsigned int totaluncorr = 0;
  };

  template <DataType dt, typename TCollision, typename TTracks>
  void processCollision(const TCollision& collision, const TTracks& tracks, const float& centrality, const float& field)
  {
    if (tracks.size() < 1) {
      return;
    }
    if (dt != Gen && (centrality < gfwMemberCache.centbinning.front() || centrality > gfwMemberCache.centbinning.back())) {
      return;
    }
    if (dt != Gen) {
      registry.fill(HIST("eventQA/eventSel"), TrackCent);
    }
    float vtxz = collision.posZ();

    fGFW->Clear();
    fFCpt->clearVector();

    float lRandom = fRndm->Rndm();

    // process tracks
    AcceptedTracks acceptedTracks;
    for (const auto& track : tracks) {
      processTrack(track, vtxz, field, centrality, acceptedTracks);
    }
    if (dt != Gen && cfgFillQA) {
      registry.fill(HIST("trackQA/after/Nch_corrected"), acceptedTracks.total);
      registry.fill(HIST("trackQA/after/Nch_uncorrected"), acceptedTracks.totaluncorr);
    }

    float multiplicity = 0.f;
    switch (cfgUseNchCorrection) {
      case 0:
        multiplicity = tracks.size();
        break;
      case 1:
        multiplicity = acceptedTracks.total;
        break;
      case 2:
        multiplicity = acceptedTracks.totaluncorr;
        break;
      default:
        multiplicity = tracks.size();
        break;
    }

    fillOutputContainers<dt>((cfgUseNch) ? multiplicity : centrality, lRandom, multiplicity);
  }

  template <typename TTrack>
  inline void processTrack(const TTrack& track, const float& vtxz, const float& field, const float& centrality, AcceptedTracks& acceptedTracks)
  {
    if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TTrack::all_columns>) {
      if (track.mcParticleId() < 0 || !(track.has_mcParticle())) {
        return;
      }
      auto mcParticle = track.mcParticle();
      if (!mcParticle.isPhysicalPrimary()) {
        return;
      }
      if (cfgFillQA) {
        fillTrackQA<Reco, Before>(track, vtxz);
      }
      if (mcParticle.eta() < gfwMemberCache.etalow || mcParticle.eta() > gfwMemberCache.etaup || mcParticle.pt() < gfwMemberCache.ptlow || mcParticle.pt() > gfwMemberCache.ptup) {
        return;
      }
      // Select tracks with nominal cuts always
      if (!nchSelected(track)) {
        return;
      }
      double weffCh = getEfficiency(track, centrality);
      if (track.eta() > cfgKinematics.cfgEtaNch->first && track.eta() < cfgKinematics.cfgEtaNch->second) {
        if (weffCh > 0) {
          acceptedTracks.total += (cfgUseNchCorrection) ? weffCh : 1.0;
        }
        ++acceptedTracks.totaluncorr;
      }

      if (!trackSelected(track, field)) {
        return;
      }

      fillPtSums<Reco>(track, centrality);
      fillGFW<Reco>(mcParticle, centrality, vtxz);
      if (cfgFillQA) {
        fillTrackQA<Reco, After>(track, vtxz);
      }

    } else if constexpr (framework::has_type_v<aod::mcparticle::McCollisionId, typename TTrack::all_columns>) {
      if (!track.isPhysicalPrimary()) {
        return;
      }
      if (cfgFillQA) {
        fillTrackQA<Gen, Before>(track, vtxz);
      }
      if (track.eta() < gfwMemberCache.etalow || track.eta() > gfwMemberCache.etaup || track.pt() < gfwMemberCache.ptlow || track.pt() > gfwMemberCache.ptup) {
        return;
      }

      if (track.eta() > cfgKinematics.cfgEtaNch->first && track.eta() < cfgKinematics.cfgEtaNch->second) {
        ++acceptedTracks.total;
        ++acceptedTracks.totaluncorr;
      }

      fillPtSums<Gen>(track, centrality);
      fillGFW<Gen>(track, centrality, vtxz);
      if (cfgFillQA) {
        fillTrackQA<Gen, After>(track, vtxz);
      }

    } else {
      if (cfgFillQA) {
        fillTrackQA<Reco, Before>(track, vtxz);
      }
      // Select tracks with nominal cuts always
      if (!nchSelected(track)) {
        return;
      }
      double weffCh = getEfficiency(track, centrality);
      if (track.eta() > cfgKinematics.cfgEtaNch->first && track.eta() < cfgKinematics.cfgEtaNch->second) {
        if (weffCh > 0) {
          acceptedTracks.total += (cfgUseNchCorrection) ? weffCh : 1.0;
        }
        ++acceptedTracks.totaluncorr;
      }
      if (!trackSelected(track, field)) {
        return;
      }

      fillPtSums<Reco>(track, centrality);
      fillGFW<Reco>(track, centrality, vtxz);

      if (cfgFillQA) {
        fillTrackQA<Reco, After>(track, vtxz);
      }
    }
  }

  o2::framework::expressions::Filter collisionFilter = nabs(aod::collision::posZ) < cfgEventSelection.cfgVtxZ;
  o2::framework::expressions::Filter trackFilter = (aod::track::eta > cfgKinematics.cfgEta->first) && (aod::track::eta < cfgKinematics.cfgEta->second) && (aod::track::pt > cfgKinematics.cfgPtCut->first) && (aod::track::pt < cfgKinematics.cfgPtCut->second) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true)) && (aod::track::itsChi2NCl < cfgTrackCuts.cfgChi2PrITSCls) && (aod::track::tpcChi2NCl < cfgTrackCuts.cfgChi2PrTPCCls) && nabs(aod::track::dcaZ) < cfgTrackCuts.cfgDCAz;

  using GFWCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs>>;
  using GFWMCCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs, aod::McCollisionLabels>;
  using GFWTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;
  using GFWMCTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>>;

  SliceCache cache;
  Partition<GFWTracks> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<GFWTracks> negTracks = aod::track::signed1Pt < 0.0f;

  double massKaPlus = o2::constants::physics::MassKPlus;
  double massLambda = o2::constants::physics::MassLambda;
  double massK0Short = o2::constants::physics::MassK0Short;

  void processData(GFWCollisions::iterator const& collision, aod::BCsWithTimestamps const&, GFWTracks const& tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    if (run != lastRun) {
      lastRun = run;
      LOGF(info, "run = %d", run);
      if (cfgRunByRun) {
        loadCorrections(bc);
      }
    }
    if (!cfgRunByRun) {
      loadCorrections(bc);
    }
    registry.fill(HIST("eventQA/eventSel"), 0.5);
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("eventQA/eventSel"), 1.5);

    float centrality = getCentrality(collision);
    if (cfgEventSelection.cfgDoOccupancySel) {
      int occupancy = collision.trackOccupancyInTimeRange();
      if (cfgFillQA) {
        registry.fill(HIST("eventQA/before/occ_mult_cent"), occupancy, tracks.size(), centrality);
      }
      if (occupancy < 0 || occupancy > cfgEventSelection.cfgOccupancySelection) {
        return;
      }
      if (cfgFillQA) {
        registry.fill(HIST("eventQA/after/occ_mult_cent"), occupancy, tracks.size(), centrality);
      }
    }
    registry.fill(HIST("eventQA/eventSel"), 2.5);

    if (cfgFillQA) {
      fillEventQA<Before>(collision, tracks);
    }
    registry.fill(HIST("eventQA/before/centrality"), centrality);
    registry.fill(HIST("eventQA/before/multiplicity"), tracks.size());
    if (!eventSelected(collision, tracks.size(), centrality)) {
      return;
    }
    if (cfgFillQA) {
      fillEventQA<After>(collision, tracks);
    }
    registry.fill(HIST("eventQA/after/centrality"), centrality);
    registry.fill(HIST("eventQA/after/multiplicity"), tracks.size());
    // Get magnetic field polarity
    auto field = (cfgEventSelection.cfgMagField == DefaultMagneticFieldCut) ? getMagneticField(bc.timestamp()) : static_cast<int>(cfgEventSelection.cfgMagField);

    processCollision<Reco>(collision, tracks, centrality, field);
  }
  PROCESS_SWITCH(FlowGfwNonflow, processData, "Process analysis for non-derived data", true);

  void processMCReco(GFWCollisions::iterator const& collision, aod::BCsWithTimestamps const&, GFWMCTracks const& tracks, aod::McParticles const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    if (run != lastRun) {
      lastRun = run;
      LOGF(info, "run = %d", run);
      if (cfgRunByRun) {
        loadCorrections(bc);
      }
    }
    if (!cfgRunByRun) {
      loadCorrections(bc);
    }

    registry.fill(HIST("eventQA/eventSel"), 0.5);

    if (!collision.sel8()) {
      return;
    }

    registry.fill(HIST("eventQA/eventSel"), 1.5);

    const auto centrality = getCentrality(collision);

    if (cfgEventSelection.cfgDoOccupancySel) {
      int occupancy = collision.trackOccupancyInTimeRange();
      registry.fill(HIST("eventQA/before/occ_mult_cent"), occupancy, tracks.size(), centrality);
      if (occupancy < 0 || occupancy > cfgEventSelection.cfgOccupancySelection) {
        return;
      }
      registry.fill(HIST("eventQA/after/occ_mult_cent"), occupancy, tracks.size(), centrality);
    }
    registry.fill(HIST("eventQA/eventSel"), 2.5);

    if (cfgFillQA) {
      fillEventQA<Before>(collision, tracks);
    }
    if (!eventSelected(collision, tracks.size(), centrality)) {
      return;
    }
    if (cfgFillQA) {
      fillEventQA<After>(collision, tracks);
    }
    loadCorrections(bc);
    auto field = (cfgEventSelection.cfgMagField == DefaultMagneticFieldCut) ? getMagneticField(bc.timestamp()) : static_cast<int>(cfgEventSelection.cfgMagField);
    processCollision<Reco>(collision, tracks, centrality, field);
  }
  PROCESS_SWITCH(FlowGfwNonflow, processMCReco, "Process analysis for MC reconstructed events", false);

  o2::framework::expressions::Filter mcCollFilter = nabs(aod::mccollision::posZ) < cfgEventSelection.cfgVtxZ;
  void processMCGen(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs>> const& collisions, aod::McParticles const& particles)
  {
    if (collisions.size() != 1) {
      return;
    }
    float centrality = -1;
    for (const auto& collision : collisions) {
      centrality = getCentrality(collision);
    }
    processCollision<Gen>(mcCollision, particles, centrality, -999);
  }
  PROCESS_SWITCH(FlowGfwNonflow, processMCGen, "Process analysis for MC generated events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGfwNonflow>(cfgc),
  };
}
