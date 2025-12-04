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

/// \file flowGfwLightIons.cxx
/// \brief Dedicated GFW task to analyse angular correlations in light-ion collision systems
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch

#include "FlowContainer.h"
#include "FlowPtContainer.h"
#include "GFW.h"
#include "GFWConfig.h"
#include "GFWCumulant.h"
#include "GFWPowerArray.h"
#include "GFWWeights.h"
#include "GFWWeightsList.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>

#include <TF1.h>
#include <TPDGCode.h>
#include <TProfile.h>
#include <TRandom3.h>

#include <algorithm>
#include <complex>
#include <ctime>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

namespace o2::analysis::gfw
{
std::vector<double> ptbinning = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};
float ptpoilow = 0.2, ptpoiup = 10.0;
float ptreflow = 0.2, ptrefup = 3.0;
float ptlow = 0.2, ptup = 10.0;
int etabins = 16;
float etalow = -0.8, etaup = 0.8;
int vtxZbins = 40;
float vtxZlow = -10.0, vtxZup = 10.0;
int phibins = 72;
float philow = 0.0;
float phiup = o2::constants::math::TwoPI;
int nchbins = 300;
float nchlow = 0;
float nchup = 3000;
std::vector<double> centbinning(90);
int nBootstrap = 10;
GFWRegions regions;
GFWCorrConfigs configs;
std::vector<double> multGlobalCorrCutPars;
std::vector<double> multPVCorrCutPars;
std::vector<double> multGlobalPVCorrCutPars;
std::vector<double> multGlobalV0ACutPars;
std::vector<double> multGlobalT0ACutPars;
std::vector<int> firstRunsOfFill;
} // namespace o2::analysis::gfw

struct FlowGfwLightIons {

  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgMpar, int, 4, "Highest order of pt-pt correlations")
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Do correlations as function of Nch")
  O2_DEFINE_CONFIGURABLE(cfgFillWeights, bool, false, "Fill NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgRunByRun, bool, false, "Fill histograms on a run-by-run basis")
  O2_DEFINE_CONFIGURABLE(cfgFillFlowRunByRun, bool, false, "Fill flow profile run-by-run (only for v22)")
  O2_DEFINE_CONFIGURABLE(cfgTimeDependent, bool, false, "Fill output as function of time (for contamination studies)")
  O2_DEFINE_CONFIGURABLE(cfgFirstRunsOfFill, std::vector<int>, {}, "First runs of a fill for time dependent analysis")
  O2_DEFINE_CONFIGURABLE(cfgFillQA, bool, false, "Fill QA histograms")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgUseCentralMoments, bool, true, "Use central moments in vn-pt calculations")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgDCAxyNSigma, float, 7, "Cut on number of sigma deviations from expected DCA in the transverse direction");
  O2_DEFINE_CONFIGURABLE(cfgDCAxy, std::string, "(0.0026+0.005/(x^1.01))", "Functional form of pt-dependent DCAxy cut");
  O2_DEFINE_CONFIGURABLE(cfgDCAz, float, 2, "Cut on DCA in the longitudinal direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgNTPCCls, float, 50, "Cut on number of TPC clusters found");
  O2_DEFINE_CONFIGURABLE(cfgNTPCXrows, float, 70, "Cut on number of TPC crossed rows");
  O2_DEFINE_CONFIGURABLE(cfgMinNITSCls, float, 5, "Cut on minimum number of ITS clusters found");
  O2_DEFINE_CONFIGURABLE(cfgChi2PrITSCls, float, 36, "Cut on chi^2 per ITS clusters found");
  O2_DEFINE_CONFIGURABLE(cfgChi2PrTPCCls, float, 2.5, "Cut on chi^2 per TPC clusters found");
  O2_DEFINE_CONFIGURABLE(cfgPtmin, float, 0.2, "minimum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgPtmax, float, 10, "maximum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgEta, float, 0.8, "eta cut");
  O2_DEFINE_CONFIGURABLE(cfgEtaPtPt, float, 0.4, "eta cut for pt-pt correlations used in subevent vn-pt");
  O2_DEFINE_CONFIGURABLE(cfgEtaPtPtGap, float, 0.4, "eta gap for subevent pt-pt correlations");
  O2_DEFINE_CONFIGURABLE(cfgEtaPtPtFull, float, 0.8, "eta cut for pure pt-pt correlations");
  O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10, "vertex cut (cm)");
  O2_DEFINE_CONFIGURABLE(cfgOccupancySelection, int, 2000, "Max occupancy selection, -999 to disable");
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgNoSameBunchPileupCut, bool, true, "kNoSameBunchPileupCut");
    O2_DEFINE_CONFIGURABLE(cfgIsGoodZvtxFT0vsPV, bool, true, "kIsGoodZvtxFT0vsPV");
    O2_DEFINE_CONFIGURABLE(cfgIsGoodITSLayersAll, bool, true, "kIsGoodITSLayersAll");
    O2_DEFINE_CONFIGURABLE(cfgNoCollInTimeRangeStandard, bool, true, "kNoCollInTimeRangeStandard");
    O2_DEFINE_CONFIGURABLE(cfgNoCollInRofStandard, bool, true, "kNoCollInRofStandard");
    O2_DEFINE_CONFIGURABLE(cfgNoHighMultCollInPrevRof, bool, true, "kNoHighMultCollInPrevRof");
    O2_DEFINE_CONFIGURABLE(cfgNoITSROFrameBorder, bool, true, "kNoITSROFrameBorder");
    O2_DEFINE_CONFIGURABLE(cfgNoTimeFrameBorder, bool, true, "kNoTimeFrameBorder");
    O2_DEFINE_CONFIGURABLE(cfgTVXinTRD, bool, true, "kTVXinTRD - Use kTVXinTRD (reject TRD triggered events)");
    O2_DEFINE_CONFIGURABLE(cfgIsVertexITSTPC, bool, true, "kIsVertexITSTPC - Selects collisions with at least one ITS-TPC track");
  } cfgEventCutFlags;
  O2_DEFINE_CONFIGURABLE(cfgDoOccupancySel, bool, true, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgCentralityFactor, double, 1., "Correction factor for testing centrality robustness");
  O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, true, "Use additional event cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field; default CCDB will be queried");
  O2_DEFINE_CONFIGURABLE(cfgFixedMultMin, int, 1, "Minimum for fixed nch range");
  O2_DEFINE_CONFIGURABLE(cfgFixedMultMax, int, 3000, "Maximum for fixed nch range");
  O2_DEFINE_CONFIGURABLE(cfgUseMultiplicityFlowWeights, bool, true, "Enable or disable the use of multiplicity-based event weighting");
  O2_DEFINE_CONFIGURABLE(cfgConsistentEventFlag, int, 0, "Flag to select consistent events - 0: off, 1: v2{2} gap calculable, 2: v2{4} full calculable, 4: v2{4} gap calculable, 8: v2{4} 3sub calculable");
  O2_DEFINE_CONFIGURABLE(cfgUseDensityDependentCorrection, bool, false, "Use density dependent efficiency correction based on Run 2 measurements");
  Configurable<std::vector<double>> cfgTrackDensityP0{"cfgTrackDensityP0", std::vector<double>{0.7217476707, 0.7384792571, 0.7542625668, 0.7640680200, 0.7701951667, 0.7755299053, 0.7805901710, 0.7849446786, 0.7957356586, 0.8113039262, 0.8211968966, 0.8280558878, 0.8329342135}, "parameter 0 for track density efficiency correction"};
  Configurable<std::vector<double>> cfgTrackDensityP1{"cfgTrackDensityP1", std::vector<double>{-2.169488e-05, -2.191913e-05, -2.295484e-05, -2.556538e-05, -2.754463e-05, -2.816832e-05, -2.846502e-05, -2.843857e-05, -2.705974e-05, -2.477018e-05, -2.321730e-05, -2.203315e-05, -2.109474e-05}, "parameter 1 for track density efficiency correction"};
  struct : ConfigurableGroup {
    Configurable<std::vector<double>> cfgMultGlobalCutPars{"cfgMultGlobalCutPars", std::vector<double>{2272.16, -76.6932, 1.01204, -0.00631545, 1.59868e-05, 136.336, -4.97006, 0.121199, -0.0015921, 7.66197e-06}, "Global vs FT0C multiplicity cut parameter values"};
    Configurable<std::vector<double>> cfgMultPVCutPars{"cfgMultPVCutPars", std::vector<double>{3074.43, -106.192, 1.46176, -0.00968364, 2.61923e-05, 182.128, -7.43492, 0.193901, -0.00256715, 1.22594e-05}, "PV vs FT0C multiplicity cut parameter values"};
    Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars", std::vector<double>{-0.223013, 0.715849, 0.664242, 0.0829653, -0.000503733, 1.21185e-06}, "Global vs PV multiplicity cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultCorrHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultCorrLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalPVCorrCutFunction, std::string, "[0] + [1]*x + 3*([2] + [3]*x + [4]*x*x + [5]*x*x*x)", "Functional for global vs pv multiplicity correlation cut");
  } cfgMultCorrCuts;
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalASideCorrCutFunction, std::string, "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + [10]*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", "Functional for global vs V0A multiplicity low correlation cut");
    Configurable<std::vector<double>> cfgMultGlobalV0ACutPars{"cfgMultGlobalV0ACutPars", std::vector<double>{567.785, 172.715, 0.77888, -0.00693466, 1.40564e-05, 679.853, 66.8068, -0.444332, 0.00115002, -4.92064e-07}, "Global vs FV0A multiplicity cut parameter values"};
    Configurable<std::vector<double>> cfgMultGlobalT0ACutPars{"cfgMultGlobalT0ACutPars", std::vector<double>{241.618, 61.8402, 0.348049, -0.00306078, 6.20357e-06, 315.235, 29.1491, -0.188639, 0.00044528, -9.08912e-08}, "Global vs FT0A multiplicity cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgGlobalV0ALowSigma, float, -3, "Number of sigma deviations below expected value in global vs V0A correlation");
    O2_DEFINE_CONFIGURABLE(cfgGlobalV0AHighSigma, float, 4, "Number of sigma deviations above expected value in global vs V0A correlation");
    O2_DEFINE_CONFIGURABLE(cfgGlobalT0ALowSigma, float, -3., "Number of sigma deviations below expected value in global vs T0A correlation");
    O2_DEFINE_CONFIGURABLE(cfgGlobalT0AHighSigma, float, 4, "Number of sigma deviations above expected value in global vs T0A correlation");
  } cfgGlobalAsideCorrCuts;

  Configurable<GFWBinningCuts> cfgGFWBinning{"cfgGFWBinning", {40, 16, 72, 300, 0, 3000, 0.2, 10.0, 0.2, 3.0, {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90}}, "Configuration for binning"};
  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN", "refP", "refFull"}, {-0.8, 0.4, -0.8}, {-0.4, 0.8, 0.8}, {0, 0, 0}, {1, 1, 1}}, "Configurations for GFW regions"};

  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}", "refFull {2 -2}", "refFull {2 2 -2 -2}"}, {"ChGap22", "ChGap32", "ChGap42", "ChFull22", "ChFull24"}, {0, 0, 0, 0, 0}, {15, 1, 1, 0, 0}}, "Configurations for each correlation to calculate"};

  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  struct Config {
    TH1D* mEfficiency = nullptr;
    GFWWeights* mAcceptance;
    bool correctionsLoaded = false;
  } cfg;

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<FlowPtContainer> fFCpt{FlowPtContainer("FlowPtContainer")};
  OutputObj<FlowPtContainer> fFCptFull{FlowPtContainer("FlowPtContainerFull")};
  OutputObj<FlowContainer> fFCgen{FlowContainer("FlowContainer_gen")};
  OutputObj<FlowPtContainer> fFCptgen{FlowPtContainer("FlowPtContainer_gen")};
  OutputObj<FlowPtContainer> fFCptgenFull{FlowPtContainer("FlowPtContainer_gen")};
  HistogramRegistry registry{"registry"};

  // QA outputs
  std::map<int, std::vector<std::shared_ptr<TH1>>> th1sList;
  std::map<int, std::vector<std::shared_ptr<TProfile>>> tpfsList;
  std::map<int, std::vector<std::shared_ptr<TH3>>> th3sList;
  enum OutputTH1Names {
    hPhi = 0,
    hEta,
    hVtxZ,
    hMult,
    hCent,
    hEventSel,
    kCount_TH1Names
  };
  enum OutputTProfileNames {
    pfCorr22 = 0,
    kCount_TProfileNames
  };
  // NUA outputs
  enum OutputTH3Names {
    hNUAref = 0,
    kCount_TH3Names
  };
  enum CentEstimators {
    kCentFT0C = 0,
    kCentFT0CVariant1,
    kCentFT0M,
    kCentFV0A,
    kCentNTPV,
    kCentNGlobal,
    kCentMFT
  };
  std::map<int, std::string> centNamesMap = {{kCentFT0C, "FT0C"}, {kCentFT0CVariant1, "FT0C variant1"}, {kCentFT0M, "FT0M"}, {kCentFV0A, "FV0A"}, {kCentNTPV, "NTPV"}, {kCentNGlobal, "NGlobal"}, {kCentMFT, "MFT"}};

  enum EventSelFlags {
    kFilteredEvent = 1,
    kSel8,
    kOccupancy,
    kTVXinTRD,
    kNoSameBunchPileup,
    kIsGoodZvtxFT0vsPV,
    kNoCollInTimeRangeStandard,
    kNoCollInRofStandard,
    kNoHighMultCollInPrevRof,
    kNoTimeFrameBorder,
    kNoITSROFrameBorder,
    kIsVertexITSTPC,
    kIsGoodITSLayersAll,
    kMultCuts,
    kTrackCent
  };

  // Define global variables
  // Generic Framework
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;

  TRandom3* fRndm = new TRandom3(0);
  TAxis* fSecondAxis;
  int lastRun = -1;
  std::vector<int>::iterator firstRunOfCurrentFill;
  std::vector<int> runNumbers;

  // Density dependent eff correction
  std::vector<TF1*> funcEff;
  TH1D* hFindPtBin;
  TF1* funcV2;
  TF1* funcV3;
  TF1* funcV4;
  struct DensityCorr {
    double psi2Est;
    double psi3Est;
    double psi4Est;
    double v2;
    double v3;
    double v4;
    int density;
    DensityCorr() : psi2Est(0.), psi3Est(0.), psi4Est(0.), v2(0.), v3(0.), v4(0.), density(0) {}
  };

  // region indices for consistency flag
  int posRegionIndex = -1;
  int negRegionIndex = -1;
  int fullRegionIndex = -1;
  int midRegionIndex = -1;

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultPVGlobalCutHigh = nullptr;
  TF1* fMultGlobalV0ACutLow = nullptr;
  TF1* fMultGlobalV0ACutHigh = nullptr;
  TF1* fMultGlobalT0ACutLow = nullptr;
  TF1* fMultGlobalT0ACutHigh = nullptr;

  TF1* fPtDepDCAxy = nullptr;

  o2::framework::expressions::Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZ;
  o2::framework::expressions::Filter trackFilter = nabs(aod::track::eta) < cfgEta && aod::track::pt > cfgPtmin&& aod::track::pt < cfgPtmax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::itsChi2NCl < cfgChi2PrITSCls) && (aod::track::tpcChi2NCl < cfgChi2PrTPCCls) && nabs(aod::track::dcaZ) < cfgDCAz;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  o2::framework::expressions::Filter mcCollFilter = nabs(aod::mccollision::posZ) < cfgVtxZ;
  o2::framework::expressions::Filter mcParticlesFilter = (aod::mcparticle::eta > o2::analysis::gfw::etalow && aod::mcparticle::eta < o2::analysis::gfw::etaup && aod::mcparticle::pt > o2::analysis::gfw::ptlow && aod::mcparticle::pt < o2::analysis::gfw::ptup);

  using GFWTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;

  void init(InitContext const&)
  {
    LOGF(info, "flowGfwLightIons::init()");
    o2::analysis::gfw::regions.SetNames(cfgRegions->GetNames());
    o2::analysis::gfw::regions.SetEtaMin(cfgRegions->GetEtaMin());
    o2::analysis::gfw::regions.SetEtaMax(cfgRegions->GetEtaMax());
    o2::analysis::gfw::regions.SetpTDifs(cfgRegions->GetpTDifs());
    o2::analysis::gfw::regions.SetBitmasks(cfgRegions->GetBitmasks());
    o2::analysis::gfw::configs.SetCorrs(cfgCorrConfig->GetCorrs());
    o2::analysis::gfw::configs.SetHeads(cfgCorrConfig->GetHeads());
    o2::analysis::gfw::configs.SetpTDifs(cfgCorrConfig->GetpTDifs());
    o2::analysis::gfw::configs.SetpTCorrMasks(cfgCorrConfig->GetpTCorrMasks());
    o2::analysis::gfw::regions.Print();
    o2::analysis::gfw::configs.Print();
    o2::analysis::gfw::ptbinning = cfgGFWBinning->GetPtBinning();
    o2::analysis::gfw::ptpoilow = cfgGFWBinning->GetPtPOImin();
    o2::analysis::gfw::ptpoiup = cfgGFWBinning->GetPtPOImax();
    o2::analysis::gfw::ptreflow = cfgGFWBinning->GetPtRefMin();
    o2::analysis::gfw::ptrefup = cfgGFWBinning->GetPtRefMax();
    o2::analysis::gfw::ptlow = cfgPtmin;
    o2::analysis::gfw::ptup = cfgPtmax;
    o2::analysis::gfw::etabins = cfgGFWBinning->GetEtaBins();
    o2::analysis::gfw::vtxZbins = cfgGFWBinning->GetVtxZbins();
    o2::analysis::gfw::phibins = cfgGFWBinning->GetPhiBins();
    o2::analysis::gfw::philow = 0.0f;
    o2::analysis::gfw::phiup = o2::constants::math::TwoPI;
    o2::analysis::gfw::nchbins = cfgGFWBinning->GetNchBins();
    o2::analysis::gfw::nchlow = cfgGFWBinning->GetNchMin();
    o2::analysis::gfw::nchup = cfgGFWBinning->GetNchMax();
    o2::analysis::gfw::centbinning = cfgGFWBinning->GetCentBinning();
    cfgGFWBinning->Print();
    o2::analysis::gfw::multGlobalCorrCutPars = cfgMultCorrCuts.cfgMultGlobalCutPars;
    o2::analysis::gfw::multPVCorrCutPars = cfgMultCorrCuts.cfgMultPVCutPars;
    o2::analysis::gfw::multGlobalPVCorrCutPars = cfgMultCorrCuts.cfgMultGlobalPVCutPars;
    o2::analysis::gfw::multGlobalV0ACutPars = cfgGlobalAsideCorrCuts.cfgMultGlobalV0ACutPars;
    o2::analysis::gfw::multGlobalT0ACutPars = cfgGlobalAsideCorrCuts.cfgMultGlobalT0ACutPars;
    o2::analysis::gfw::firstRunsOfFill = cfgFirstRunsOfFill;
    if (cfgTimeDependent && !std::is_sorted(o2::analysis::gfw::firstRunsOfFill.begin(), o2::analysis::gfw::firstRunsOfFill.end())) {
      std::sort(o2::analysis::gfw::firstRunsOfFill.begin(), o2::analysis::gfw::firstRunsOfFill.end());
    }
    firstRunOfCurrentFill = o2::analysis::gfw::firstRunsOfFill.begin();

    AxisSpec phiAxis = {o2::analysis::gfw::phibins, o2::analysis::gfw::philow, o2::analysis::gfw::phiup, "#phi"};
    AxisSpec etaAxis = {o2::analysis::gfw::etabins, -cfgEta, cfgEta, "#eta"};
    AxisSpec vtxAxis = {o2::analysis::gfw::vtxZbins, -cfgVtxZ, cfgVtxZ, "Vtx_{z} (cm)"};
    AxisSpec ptAxis = {o2::analysis::gfw::ptbinning, "#it{p}_{T} GeV/#it{c}"};
    std::string sCentralityEstimator = centNamesMap[cfgCentEstimator] + " centrality (%)";
    AxisSpec centAxis = {o2::analysis::gfw::centbinning, sCentralityEstimator.c_str()};
    std::vector<double> nchbinning;
    int nchskip = (o2::analysis::gfw::nchup - o2::analysis::gfw::nchlow) / o2::analysis::gfw::nchbins;
    for (int i = 0; i <= o2::analysis::gfw::nchbins; ++i) {
      nchbinning.push_back(nchskip * i + o2::analysis::gfw::nchlow + 0.5);
    }
    AxisSpec nchAxis = {nchbinning, "N_{ch}"};
    std::vector<double> bbinning(201);
    std::generate(bbinning.begin(), bbinning.end(), [n = -0.1, step = 0.1]() mutable {
      n += step;
      return n;
    });
    AxisSpec bAxis = {bbinning, "#it{b}"};
    AxisSpec t0cAxis = {1000, 0, 10000, "N_{ch} (T0C)"};
    AxisSpec t0aAxis = {300, 0, 30000, "N_{ch} (T0A)"};
    AxisSpec v0aAxis = {800, 0, 80000, "N_{ch} (V0A)"};
    AxisSpec multpvAxis = {600, 0, 600, "N_{ch} (PV)"};
    AxisSpec dcaZAXis = {200, -2, 2, "DCA_{z} (cm)"};
    AxisSpec dcaXYAXis = {200, -0.5, 0.5, "DCA_{xy} (cm)"};
    std::vector<double> timebinning(289);
    std::generate(timebinning.begin(), timebinning.end(), [n = -24 / 288., step = 24 / 288.]() mutable {
      n += step;
      return n;
    });
    AxisSpec timeAxis = {timebinning, "time (hrs)"};

    AxisSpec multAxis = (cfgTimeDependent) ? timeAxis : (doprocessOnTheFly && !cfgUseNch) ? bAxis
                                                      : (cfgUseNch)                       ? nchAxis
                                                                                          : centAxis;
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    int ptbins = o2::analysis::gfw::ptbinning.size() - 1;
    fSecondAxis = (cfgTimeDependent) ? new TAxis(timeAxis.binEdges.size() - 1, &(timeAxis.binEdges[0])) : new TAxis(ptbins, &o2::analysis::gfw::ptbinning[0]);

    if (doprocessMCGen || doprocessOnTheFly) {
      registry.add("MCGen/trackQA/nch_pt", "#it{p}_{T} vs multiplicity; N_{ch}; #it{p}_{T}", {HistType::kTH2D, {nchAxis, ptAxis}});
      registry.add("MCGen/trackQA/phi_eta_vtxZ", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registry.add("MCGen/trackQA/pt_ref", "Reference #it{p}_{T}; #it{p}_{T}; Counts", {HistType::kTH1D, {{100, o2::analysis::gfw::ptreflow, o2::analysis::gfw::ptrefup}}});
      registry.add("MCGen/trackQA/pt_poi", "POI #it{p}_{T}; #it{p}_{T}; Counts", {HistType::kTH1D, {{100, o2::analysis::gfw::ptpoilow, o2::analysis::gfw::ptpoiup}}});
      if (doprocessOnTheFly)
        registry.add("MCGen/impactParameter", "", {HistType::kTH2D, {{bAxis, nchAxis}}});

      registry.add("MCGen/eventQA/multiplicity", "", {HistType::kTH1D, {nchAxis}});
      if (doprocessMCGen)
        registry.add("MCGen/eventQA/centrality", "", {HistType::kTH1D, {centAxis}});
    }
    if (doprocessMCReco || doprocessData) {
      registry.add("trackQA/before/phi_eta_vtxZ", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registry.add("trackQA/before/pt_dcaXY_dcaZ", "", {HistType::kTH3D, {ptAxis, dcaXYAXis, dcaZAXis}});
      registry.add("trackQA/before/nch_pt", "#it{p}_{T} vs multiplicity; N_{ch}; #it{p}_{T}", {HistType::kTH2D, {nchAxis, ptAxis}});
      registry.add("trackQA/before/chi2prTPCcls", "#chi^{2}/cluster for the TPC track segment; #chi^{2}/TPC cluster", {HistType::kTH1D, {{100, 0., 5.}}});
      registry.add("trackQA/before/chi2prITScls", "#chi^{2}/cluster for the ITS track; #chi^{2}/ITS cluster", {HistType::kTH1D, {{100, 0., 50.}}});
      registry.add("trackQA/before/nTPCClusters", "Number of found TPC clusters; TPC N_{cls}; Counts", {HistType::kTH1D, {{100, 40, 180}}});
      registry.add("trackQA/before/nITSClusters", "Number of found ITS clusters; ITS N_{cls}; Counts", {HistType::kTH1D, {{100, 0, 20}}});
      registry.add("trackQA/before/nTPCCrossedRows", "Number of crossed TPC Rows; TPC X-rows; Counts", {HistType::kTH1D, {{100, 40, 180}}});

      registry.addClone("trackQA/before/", "trackQA/after/");
      registry.add("trackQA/after/pt_ref", "", {HistType::kTH1D, {{100, o2::analysis::gfw::ptreflow, o2::analysis::gfw::ptrefup}}});
      registry.add("trackQA/after/pt_poi", "", {HistType::kTH1D, {{100, o2::analysis::gfw::ptpoilow, o2::analysis::gfw::ptpoiup}}});

      registry.add("eventQA/before/multiplicity", "", {HistType::kTH1D, {nchAxis}});
      if (cfgTimeDependent) {
        registry.add("eventQA/before/multiplicity_time", "Multiplicity vs time; time (hrs); N_{ch}", {HistType::kTH2D, {timeAxis, nchAxis}});
        registry.add("eventQA/before/multT0C_time", "T0C Multiplicity vs time; time (hrs); N_{ch} (T0C)", {HistType::kTH2D, {timeAxis, t0cAxis}});
        registry.add("eventQA/before/multT0A_time", "T0A Multiplicity vs time; time (hrs); N_{ch} (T0A)", {HistType::kTH2D, {timeAxis, t0aAxis}});
        registry.add("eventQA/before/multV0A_time", "V0A Multiplicity vs time; time (hrs); N_{ch} (V0A)", {HistType::kTH2D, {timeAxis, v0aAxis}});
        registry.add("eventQA/before/multPV_time", "PV Multiplicity vs time; time (hrs); N_{ch} (PV)", {HistType::kTH2D, {timeAxis, multpvAxis}});
      }
      registry.add("eventQA/before/globalTracks_PVTracks", "", {HistType::kTH2D, {multpvAxis, nchAxis}});
      registry.add("eventQA/before/globalTracks_multT0A", "", {HistType::kTH2D, {t0aAxis, nchAxis}});
      registry.add("eventQA/before/globalTracks_multV0A", "", {HistType::kTH2D, {v0aAxis, nchAxis}});
      registry.add("eventQA/before/multV0A_multT0A", "", {HistType::kTH2D, {t0aAxis, v0aAxis}});

      if (doprocessData || doprocessMCReco) {
        registry.add("eventQA/before/centrality", "", {HistType::kTH1D, {centAxis}});
        registry.add("eventQA/before/globalTracks_centT0C", "", {HistType::kTH2D, {centAxis, nchAxis}});
        registry.add("eventQA/before/PVTracks_centT0C", "", {HistType::kTH2D, {centAxis, multpvAxis}});
        registry.add("eventQA/before/multT0C_centT0C", "", {HistType::kTH2D, {centAxis, t0cAxis}});

        registry.add("eventQA/before/centT0M_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
        registry.add("eventQA/before/centV0A_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
        registry.add("eventQA/before/centGlobal_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
        registry.add("eventQA/before/centNTPV_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
        registry.add("eventQA/before/centMFT_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
      }

      registry.addClone("eventQA/before/", "eventQA/after/");
      registry.add("eventQA/eventSel", "Number of Events;; Counts", {HistType::kTH1D, {{15, 0.5, 15.5}}});
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kFilteredEvent, "Filtered event");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kSel8, "sel8");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kOccupancy, "occupancy");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kTVXinTRD, "kTVXinTRD");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoSameBunchPileup, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kIsGoodZvtxFT0vsPV, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoCollInTimeRangeStandard, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoCollInRofStandard, "kNoCollInRofStandard");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoHighMultCollInPrevRof, "kNoHighMultCollInPrevRof");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoTimeFrameBorder, "kNoTimeFrameBorder");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoITSROFrameBorder, "kNoITSROFrameBorder");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kIsVertexITSTPC, "kIsVertexITSTPC");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kIsGoodITSLayersAll, "kIsGoodITSLayersAll");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kMultCuts, "after Mult cuts");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kTrackCent, "has track + within cent");
      LOGF(info, "eventsel N bins = %d", registry.get<TH1>(HIST("eventQA/eventSel"))->GetNbinsX());
      if (!cfgRunByRun && cfgFillWeights) {
        registry.add<TH3>("phi_eta_vtxz_ref", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      }
    }

    if (o2::analysis::gfw::regions.GetSize() < 0)
      LOGF(error, "Configuration contains vectors of different size - check the GFWRegions configurable");
    for (auto i(0); i < o2::analysis::gfw::regions.GetSize(); ++i) {
      fGFW->AddRegion(o2::analysis::gfw::regions.GetNames()[i], o2::analysis::gfw::regions.GetEtaMin()[i], o2::analysis::gfw::regions.GetEtaMax()[i], (o2::analysis::gfw::regions.GetpTDifs()[i]) ? ptbins + 1 : 1, o2::analysis::gfw::regions.GetBitmasks()[i]);
    }
    for (auto i = 0; i < o2::analysis::gfw::configs.GetSize(); ++i) {
      corrconfigs.push_back(fGFW->GetCorrelatorConfig(o2::analysis::gfw::configs.GetCorrs()[i], o2::analysis::gfw::configs.GetHeads()[i], o2::analysis::gfw::configs.GetpTDifs()[i]));
    }
    if (corrconfigs.empty())
      LOGF(error, "Configuration contains vectors of different size - check the GFWCorrConfig configurable");
    fGFW->CreateRegions();
    TObjArray* oba = new TObjArray();
    addConfigObjectsToObjArray(oba, corrconfigs);
    if (doprocessData || doprocessMCReco) {
      fFC->SetName("FlowContainer");
      fFC->SetXAxis(fSecondAxis);
      fFC->Initialize(oba, multAxis, cfgNbootstrap);
    }
    if (doprocessMCGen || doprocessOnTheFly) {
      fFCgen->SetName("FlowContainer_gen");
      fFCgen->SetXAxis(fSecondAxis);
      fFCgen->Initialize(oba, multAxis, cfgNbootstrap);
    }
    delete oba;
    fFCpt->setUseCentralMoments(cfgUseCentralMoments);
    fFCpt->setUseGapMethod(true);
    fFCpt->initialise(multAxis, cfgMpar, o2::analysis::gfw::configs, cfgNbootstrap);
    fFCptFull->setUseCentralMoments(cfgUseCentralMoments);
    fFCptFull->setUseGapMethod(true);
    fFCptFull->initialise(multAxis, cfgMpar, o2::analysis::gfw::configs, cfgNbootstrap);
    fFCptFull->initialiseSubevent(multAxis, cfgMpar, cfgNbootstrap);
    fFCptgen->setUseCentralMoments(cfgUseCentralMoments);
    fFCptgen->setUseGapMethod(true);
    fFCptgen->initialise(multAxis, cfgMpar, o2::analysis::gfw::configs, cfgNbootstrap);
    fFCptgenFull->setUseCentralMoments(cfgUseCentralMoments);
    fFCptgenFull->setUseGapMethod(true);
    fFCptgenFull->initialise(multAxis, cfgMpar, o2::analysis::gfw::configs, cfgNbootstrap);
    fFCptgenFull->initialiseSubevent(multAxis, cfgMpar, cfgNbootstrap);

    fPtDepDCAxy = new TF1("ptDepDCAxy", Form("[0]*%s", cfgDCAxy->c_str()), 0.001, 100);
    fPtDepDCAxy->SetParameter(0, cfgDCAxyNSigma);
    LOGF(info, "DCAxy pt-dependence function: %s", Form("[0]*%s", cfgDCAxy->c_str()));
    if (cfgUseAdditionalEventCut) {
      fMultPVCutLow = new TF1("fMultPVCutLow", cfgMultCorrCuts.cfgMultCorrLowCutFunction->c_str(), 0, 100);
      fMultPVCutLow->SetParameters(&(o2::analysis::gfw::multPVCorrCutPars[0]));
      fMultPVCutHigh = new TF1("fMultPVCutHigh", cfgMultCorrCuts.cfgMultCorrHighCutFunction->c_str(), 0, 100);
      fMultPVCutHigh->SetParameters(&(o2::analysis::gfw::multPVCorrCutPars[0]));
      fMultCutLow = new TF1("fMultCutLow", cfgMultCorrCuts.cfgMultCorrLowCutFunction->c_str(), 0, 100);
      fMultCutLow->SetParameters(&(o2::analysis::gfw::multGlobalCorrCutPars[0]));
      fMultCutHigh = new TF1("fMultCutHigh", cfgMultCorrCuts.cfgMultCorrHighCutFunction->c_str(), 0, 100);
      fMultCutHigh->SetParameters(&(o2::analysis::gfw::multGlobalCorrCutPars[0]));
      fMultPVGlobalCutHigh = new TF1("fMultPVGlobalCutHigh", cfgMultCorrCuts.cfgMultGlobalPVCorrCutFunction->c_str(), 0, nchbinning.back());
      fMultPVGlobalCutHigh->SetParameters(&(o2::analysis::gfw::multGlobalPVCorrCutPars[0]));

      LOGF(info, "Global V0A function: %s in range 0-%g", cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), v0aAxis.binEdges.back());
      fMultGlobalV0ACutLow = new TF1("fMultGlobalV0ACutLow", cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, v0aAxis.binEdges.back());
      for (std::size_t i = 0; i < o2::analysis::gfw::multGlobalV0ACutPars.size(); ++i)
        fMultGlobalV0ACutLow->SetParameter(i, o2::analysis::gfw::multGlobalV0ACutPars[i]);
      fMultGlobalV0ACutLow->SetParameter(o2::analysis::gfw::multGlobalV0ACutPars.size(), cfgGlobalAsideCorrCuts.cfgGlobalV0ALowSigma);
      for (int i = 0; i < fMultGlobalV0ACutLow->GetNpar(); ++i)
        LOGF(info, "fMultGlobalV0ACutLow par %d = %g", i, fMultGlobalV0ACutLow->GetParameter(i));

      fMultGlobalV0ACutHigh = new TF1("fMultGlobalV0ACutHigh", cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, v0aAxis.binEdges.back());
      for (std::size_t i = 0; i < o2::analysis::gfw::multGlobalV0ACutPars.size(); ++i)
        fMultGlobalV0ACutHigh->SetParameter(i, o2::analysis::gfw::multGlobalV0ACutPars[i]);
      fMultGlobalV0ACutHigh->SetParameter(o2::analysis::gfw::multGlobalV0ACutPars.size(), cfgGlobalAsideCorrCuts.cfgGlobalV0AHighSigma);
      for (int i = 0; i < fMultGlobalV0ACutHigh->GetNpar(); ++i)
        LOGF(info, "fMultGlobalV0ACutHigh par %d = %g", i, fMultGlobalV0ACutHigh->GetParameter(i));

      LOGF(info, "Global T0A function: %s", cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str());
      fMultGlobalT0ACutLow = new TF1("fMultGlobalT0ACutLow", cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, t0aAxis.binEdges.back());
      for (std::size_t i = 0; i < o2::analysis::gfw::multGlobalT0ACutPars.size(); ++i)
        fMultGlobalT0ACutLow->SetParameter(i, o2::analysis::gfw::multGlobalT0ACutPars[i]);
      fMultGlobalT0ACutLow->SetParameter(o2::analysis::gfw::multGlobalT0ACutPars.size(), cfgGlobalAsideCorrCuts.cfgGlobalT0ALowSigma);
      for (int i = 0; i < fMultGlobalT0ACutLow->GetNpar(); ++i)
        LOGF(info, "fMultGlobalT0ACutLow par %d = %g", i, fMultGlobalT0ACutLow->GetParameter(i));

      fMultGlobalT0ACutHigh = new TF1("fMultGlobalT0ACutHigh", cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, t0aAxis.binEdges.back());
      for (std::size_t i = 0; i < o2::analysis::gfw::multGlobalT0ACutPars.size(); ++i)
        fMultGlobalT0ACutHigh->SetParameter(i, o2::analysis::gfw::multGlobalT0ACutPars[i]);
      fMultGlobalT0ACutHigh->SetParameter(o2::analysis::gfw::multGlobalT0ACutPars.size(), cfgGlobalAsideCorrCuts.cfgGlobalT0AHighSigma);
      for (int i = 0; i < fMultGlobalT0ACutHigh->GetNpar(); ++i)
        LOGF(info, "fMultGlobalT0ACutHigh par %d = %g", i, fMultGlobalT0ACutHigh->GetParameter(i));
    }
    if (cfgUseDensityDependentCorrection) {
      std::vector<double> pTEffBins = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0};
      hFindPtBin = new TH1D("hFindPtBin", "hFindPtBin", pTEffBins.size() - 1, &pTEffBins[0]);
      funcEff.resize(pTEffBins.size() - 1);
      // LHC24g3 Eff
      std::vector<double> f1p0 = cfgTrackDensityP0;
      std::vector<double> f1p1 = cfgTrackDensityP1;
      for (uint ifunc = 0; ifunc < pTEffBins.size() - 1; ifunc++) {
        funcEff[ifunc] = new TF1(Form("funcEff%i", ifunc), "[0]+[1]*x", 0, 3000);
        funcEff[ifunc]->SetParameters(f1p0[ifunc], f1p1[ifunc]);
      }
      funcV2 = new TF1("funcV2", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV2->SetParameters(0.0186111, 0.00351907, -4.38264e-05, 1.35383e-07, -3.96266e-10);
      funcV3 = new TF1("funcV3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV3->SetParameters(0.0174056, 0.000703329, -1.45044e-05, 1.91991e-07, -1.62137e-09);
      funcV4 = new TF1("funcV4", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV4->SetParameters(0.008845, 0.000259668, -3.24435e-06, 4.54837e-08, -6.01825e-10);
    }
    if (cfgConsistentEventFlag) {
      posRegionIndex = [&]() {
        auto begin = cfgRegions->GetNames().begin();
        auto end = cfgRegions->GetNames().end();
        auto it = std::find(begin, end, "refP");
        return (it != end) ? std::distance(begin, it) : -1;
      }();
      negRegionIndex = [&]() {
        auto begin = cfgRegions->GetNames().begin();
        auto end = cfgRegions->GetNames().end();
        auto it = std::find(begin, end, "refN");
        return (it != end) ? std::distance(begin, it) : -1;
      }();
      fullRegionIndex = [&]() {
        auto begin = cfgRegions->GetNames().begin();
        auto end = cfgRegions->GetNames().end();
        auto it = std::find(begin, end, "refFull");
        return (it != end) ? std::distance(begin, it) : -1;
      }();
      midRegionIndex = [&]() {
        auto begin = cfgRegions->GetNames().begin();
        auto end = cfgRegions->GetNames().end();
        auto it = std::find(begin, end, "refMid");
        return (it != end) ? std::distance(begin, it) : -1;
      }();
    }
  }

  static constexpr std::string_view FillTimeName[] = {"before/", "after/"};

  enum QAFillTime {
    kBefore,
    kAfter
  };

  void addConfigObjectsToObjArray(TObjArray* oba, const std::vector<GFW::CorrConfig>& configs)
  {
    for (auto it = configs.begin(); it != configs.end(); ++it) {
      if (it->pTDif) {
        std::string suffix = "_ptDiff";
        for (auto i = 0; i < fSecondAxis->GetNbins(); ++i) {
          std::string index = Form("_pt_%i", i + 1);
          oba->Add(new TNamed(it->Head.c_str() + index, it->Head.c_str() + suffix));
        }
      } else {
        oba->Add(new TNamed(it->Head.c_str(), it->Head.c_str()));
      }
    }
  }

  int getMagneticField(uint64_t timestamp)
  {
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
    if (!cfgRunByRun && cfg.correctionsLoaded)
      return;
    if (!cfgAcceptance.value.empty()) {
      std::string runstr = (cfgRunByRun) ? "RunByRun/" : "";
      cfg.mAcceptance = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + runstr, timestamp);
    }
    if (!cfgEfficiency.value.empty()) {
      cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgEfficiency, timestamp);
      if (cfg.mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)cfg.mEfficiency);
    }
    cfg.correctionsLoaded = true;
  }

  template <typename TTrack>
  double getAcceptance(TTrack track, const double& vtxz)
  {
    double wacc = 1;
    if (cfg.mAcceptance)
      wacc = cfg.mAcceptance->getNUA(track.phi(), track.eta(), vtxz);
    return wacc;
  }

  template <typename TTrack>
  double getEfficiency(TTrack track)
  {
    double eff = 1.;
    if (cfg.mEfficiency)
      eff = cfg.mEfficiency->GetBinContent(cfg.mEfficiency->FindBin(track.pt()));
    if (eff == 0)
      return -1.;
    else
      return 1. / eff;
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int& multTrk, const float& centrality, const int& run)
  {
    if (cfgEventCutFlags.cfgTVXinTRD) {
      if (collision.alias_bit(kTVXinTRD)) {
        // TRD triggered
        // "CMTVX-B-NOPF-TRD,minbias_TVX"
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kTVXinTRD);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kTVXinTRD);
    }
    if (cfgEventCutFlags.cfgNoSameBunchPileupCut) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kNoSameBunchPileup);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kNoSameBunchPileup);
    }
    if (cfgEventCutFlags.cfgIsGoodZvtxFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kIsGoodZvtxFT0vsPV);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kIsGoodZvtxFT0vsPV);
    }
    if (cfgEventCutFlags.cfgNoCollInTimeRangeStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kNoCollInTimeRangeStandard);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kNoCollInTimeRangeStandard);
    }

    if (cfgEventCutFlags.cfgNoCollInRofStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kNoCollInRofStandard);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kNoCollInRofStandard);
    }

    if (cfgEventCutFlags.cfgNoHighMultCollInPrevRof) {
      if (!collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kNoHighMultCollInPrevRof);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kNoHighMultCollInPrevRof);
    }

    if (cfgEventCutFlags.cfgIsVertexITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        // selects collisions with at least one ITS-TPC track, and thus rejects vertices built from ITS-only tracks
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kIsVertexITSTPC);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kIsVertexITSTPC);
    }

    if (cfgEventCutFlags.cfgIsGoodITSLayersAll) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kIsGoodITSLayersAll);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kIsGoodITSLayersAll);
    }

    if (cfgEventCutFlags.cfgNoTimeFrameBorder) {
      if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kNoTimeFrameBorder);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kNoTimeFrameBorder);
    }

    if (cfgEventCutFlags.cfgNoITSROFrameBorder) {
      if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kNoITSROFrameBorder);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kNoITSROFrameBorder);
    }

    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = std::sqrt(collision.covZZ());
      float minZRes = 0.25;
      int minNContrib = 20;
      if (zRes > minZRes && collision.numContrib() < minNContrib)
        vtxz = -999;
    }
    auto multNTracksPV = collision.multNTracksPV();

    if (vtxz > o2::analysis::gfw::vtxZup || vtxz < o2::analysis::gfw::vtxZlow)
      return 0;

    if (cfgMultCut) {
      if (multNTracksPV < fMultPVCutLow->Eval(centrality))
        return 0;
      if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
        return 0;
      if (multTrk < fMultCutLow->Eval(centrality))
        return 0;
      if (multTrk > fMultCutHigh->Eval(centrality))
        return 0;
      if (multTrk > fMultPVGlobalCutHigh->Eval(collision.multNTracksPV()))
        return 0;

      if (!(cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->empty()) && static_cast<double>(collision.multFV0A()) < fMultGlobalV0ACutLow->Eval(multTrk))
        return 0;
      if (!(cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->empty()) && static_cast<double>(collision.multFV0A()) > fMultGlobalV0ACutHigh->Eval(multTrk))
        return 0;
      if (!(cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->empty()) && static_cast<double>(collision.multFT0A()) < fMultGlobalT0ACutLow->Eval(multTrk))
        return 0;
      if (!(cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->empty()) && static_cast<double>(collision.multFT0A()) > fMultGlobalT0ACutHigh->Eval(multTrk))
        return 0;
      registry.fill(HIST("eventQA/eventSel"), kMultCuts);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kMultCuts);
    }
    return 1;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    if (cfgDCAxyNSigma && (std::fabs(track.dcaXY()) > fPtDepDCAxy->Eval(track.pt())))
      return false;
    return ((track.tpcNClsCrossedRows() >= cfgNTPCXrows) && (track.tpcNClsFound() >= cfgNTPCCls) && (track.itsNCls() >= cfgMinNITSCls));
  }

  enum DataType {
    kReco,
    kGen
  };

  template <typename TTrack>
  void fillWeights(const TTrack track, const double vtxz, const int& run)
  {
    if (cfgRunByRun)
      th3sList[run][hNUAref]->Fill(track.phi(), track.eta(), vtxz);
    else
      registry.fill(HIST("phi_eta_vtxz_ref"), track.phi(), track.eta(), vtxz);
    return;
  }

  void createRunByRunHistograms(const int& run)
  {
    AxisSpec phiAxis = {o2::analysis::gfw::phibins, o2::analysis::gfw::philow, o2::analysis::gfw::phiup, "#phi"};
    AxisSpec etaAxis = {o2::analysis::gfw::etabins, -cfgEta, cfgEta, "#eta"};
    AxisSpec vtxAxis = {o2::analysis::gfw::vtxZbins, -cfgVtxZ, cfgVtxZ, "Vtx_{z} (cm)"};
    AxisSpec nchAxis = {o2::analysis::gfw::nchbins, o2::analysis::gfw::nchlow, o2::analysis::gfw::nchup, "N_{ch}"};
    AxisSpec centAxis = {o2::analysis::gfw::centbinning, "Centrality (%)"};
    std::vector<std::shared_ptr<TH1>> histos(kCount_TH1Names);
    histos[hPhi] = registry.add<TH1>(Form("%d/phi", run), "", {HistType::kTH1D, {phiAxis}});
    histos[hEta] = registry.add<TH1>(Form("%d/eta", run), "", {HistType::kTH1D, {etaAxis}});
    histos[hVtxZ] = registry.add<TH1>(Form("%d/vtxz", run), "", {HistType::kTH1D, {vtxAxis}});
    histos[hMult] = registry.add<TH1>(Form("%d/mult", run), "", {HistType::kTH1D, {nchAxis}});
    histos[hCent] = registry.add<TH1>(Form("%d/cent", run), "", {HistType::kTH1D, {centAxis}});
    if (cfgFillFlowRunByRun) {
      std::vector<std::shared_ptr<TProfile>> profiles(kCount_TProfileNames);
      profiles[pfCorr22] = registry.add<TProfile>(Form("%d/corr22", run), "", {HistType::kTProfile, {(cfgUseNch) ? nchAxis : centAxis}});
      tpfsList.insert(std::make_pair(run, profiles));
    }
    histos[hEventSel] = registry.add<TH1>(Form("%d/eventSel", run), "Number of Events;; Counts", {HistType::kTH1D, {{15, 0.5, 15.5}}});
    histos[hEventSel]->GetXaxis()->SetBinLabel(kFilteredEvent, "Filtered event");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kSel8, "sel8");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kOccupancy, "occupancy");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kTVXinTRD, "kTVXinTRD");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kNoSameBunchPileup, "kNoSameBunchPileup");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kIsGoodZvtxFT0vsPV, "kIsGoodZvtxFT0vsPV");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kNoCollInTimeRangeStandard, "kNoCollInTimeRangeStandard");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kNoCollInRofStandard, "kNoCollInRofStandard");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kNoHighMultCollInPrevRof, "kNoHighMultCollInPrevRof");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kNoTimeFrameBorder, "kNoTimeFrameBorder");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kNoITSROFrameBorder, "kNoITSROFrameBorder");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kIsVertexITSTPC, "kIsVertexITSTPC");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kIsGoodITSLayersAll, "kIsGoodITSLayersAll");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kMultCuts, "after Mult cuts");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kTrackCent, "has track + within cent");

    th1sList.insert(std::make_pair(run, histos));
    std::vector<std::shared_ptr<TH3>> histos3d(kCount_TH3Names);
    histos3d[hNUAref] = registry.add<TH3>(Form("%d/phi_eta_vtxz_ref", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
    th3sList.insert(std::make_pair(run, histos3d));
    return;
  }

  template <DataType dt>
  void fillOutputContainers(const float& centmult, const double& rndm, const int& run = 0)
  {
    if (dt == kGen) {
      fFCptgen->calculateCorrelations();
      fFCptgenFull->calculateCorrelations();
      fFCptgenFull->calculateSubeventCorrelations();
    } else {
      fFCpt->calculateCorrelations();
      fFCptFull->calculateCorrelations();
      fFCptFull->calculateSubeventCorrelations();
    }
    if (dt == kGen) {
      fFCptgen->fillPtProfiles(centmult, rndm);
      fFCptgen->fillCMProfiles(centmult, rndm);
      fFCptgenFull->fillPtProfiles(centmult, rndm);
      fFCptgenFull->fillCMProfiles(centmult, rndm);
      fFCptgenFull->fillSubeventPtProfiles(centmult, rndm);
      fFCptgenFull->fillCMSubeventProfiles(centmult, rndm);
    } else {
      fFCpt->fillPtProfiles(centmult, rndm);
      fFCpt->fillCMProfiles(centmult, rndm);
      fFCptFull->fillPtProfiles(centmult, rndm);
      fFCptFull->fillSubeventPtProfiles(centmult, rndm);
      fFCptFull->fillCMProfiles(centmult, rndm);
      fFCptFull->fillCMSubeventProfiles(centmult, rndm);
    }
    for (uint l_ind = 0; l_ind < corrconfigs.size(); ++l_ind) {
      if (!corrconfigs.at(l_ind).pTDif) {
        uint8_t vnptmask = o2::analysis::gfw::configs.GetpTCorrMasks()[l_ind];
        auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), 0, kTRUE).real();
        if (dnx == 0) {
          (dt == kGen) ? fFCptgen->skipVnPtProfiles(vnptmask) : fFCpt->skipVnPtProfiles(vnptmask);
          continue;
        }
        auto val = fGFW->Calculate(corrconfigs.at(l_ind), 0, kFALSE).real() / dnx;
        if (std::abs(val) < 1) {
          (dt == kGen) ? fFCgen->FillProfile(corrconfigs.at(l_ind).Head.c_str(), centmult, val, (cfgUseMultiplicityFlowWeights) ? dnx : 1.0, rndm) : fFC->FillProfile(corrconfigs.at(l_ind).Head.c_str(), centmult, val, (cfgUseMultiplicityFlowWeights) ? dnx : 1.0, rndm);
          (dt == kGen) ? fFCptgen->fillVnPtProfiles(centmult, val, (cfgUseMultiplicityFlowWeights) ? dnx : 1.0, rndm, vnptmask) : fFCpt->fillVnPtProfiles(centmult, val, (cfgUseMultiplicityFlowWeights) ? dnx : 1.0, rndm, vnptmask);
          if (cfgRunByRun && cfgFillFlowRunByRun && dt != kGen && l_ind == 0) {
            tpfsList[run][pfCorr22]->Fill(centmult, val, (cfgUseMultiplicityFlowWeights) ? dnx : 1.0);
          }
        } else {
          (dt == kGen) ? fFCptgen->skipVnPtProfiles(vnptmask) : fFCpt->skipVnPtProfiles(vnptmask);
        }
        continue;
      }
      for (int i = 1; i <= fSecondAxis->GetNbins(); i++) {
        auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), i - 1, kTRUE).real();
        if (dnx == 0)
          continue;
        auto val = fGFW->Calculate(corrconfigs.at(l_ind), i - 1, kFALSE).real() / dnx;
        if (std::abs(val) < 1)
          (dt == kGen) ? fFCgen->FillProfile(Form("%s_pt_%i", corrconfigs.at(l_ind).Head.c_str(), i), centmult, val, (cfgUseMultiplicityFlowWeights) ? dnx : 1.0, rndm) : fFC->FillProfile(Form("%s_pt_%i", corrconfigs.at(l_ind).Head.c_str(), i), centmult, val, (cfgUseMultiplicityFlowWeights) ? dnx : 1.0, rndm);
      }
    }
    return;
  }

  struct XAxis {
    float centrality;
    int64_t multiplicity;
    double time;
  };

  struct AcceptedTracks {
    int nPos;
    int nNeg;
    int nFull;
    int nMid;
  };

  template <DataType dt, typename TCollision, typename TTracks>
  void processCollision(TCollision collision, TTracks tracks, const XAxis& xaxis, const int& run)
  {
    if (tracks.size() < 1)
      return;
    if (dt != kGen && xaxis.centrality >= 0 && (xaxis.centrality < o2::analysis::gfw::centbinning.front() || xaxis.centrality > o2::analysis::gfw::centbinning.back()))
      return;
    if (xaxis.multiplicity < cfgFixedMultMin || xaxis.multiplicity > cfgFixedMultMax)
      return;
    if (dt != kGen) {
      registry.fill(HIST("eventQA/eventSel"), kTrackCent);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kTrackCent);
    }
    if (xaxis.centrality >= 0)
      registry.fill(HIST("eventQA/after/centrality"), xaxis.centrality);
    registry.fill(HIST("eventQA/after/multiplicity"), xaxis.multiplicity);
    float vtxz = collision.posZ();
    if (dt != kGen && cfgRunByRun) {
      th1sList[run][hVtxZ]->Fill(vtxz);
      th1sList[run][hMult]->Fill(xaxis.multiplicity);
      th1sList[run][hCent]->Fill(xaxis.centrality);
    }
    fGFW->Clear();
    if (dt == kGen) {
      fFCptgen->clearVector();
      fFCptgenFull->clearVector();
    } else {
      fFCpt->clearVector();
      fFCptFull->clearVector();
    }

    float lRandom = fRndm->Rndm();

    // be cautious, this only works for Pb-Pb
    // esimate the Event plane and vn for this event
    DensityCorr densitycorrections;
    if (cfgUseDensityDependentCorrection) {
      double psi2Est = 0, psi3Est = 0, psi4Est = 0;
      double v2 = 0, v3 = 0, v4 = 0;
      double q2x = 0, q2y = 0;
      double q3x = 0, q3y = 0;
      double q4x = 0, q4y = 0;
      for (const auto& track : tracks) {
        bool withinPtRef = (o2::analysis::gfw::ptreflow < track.pt()) && (track.pt() < o2::analysis::gfw::ptrefup); // within RF pT rang
        if (withinPtRef) {
          q2x += std::cos(2 * track.phi());
          q2y += std::sin(2 * track.phi());
          q3x += std::cos(3 * track.phi());
          q3y += std::sin(3 * track.phi());
          q4x += std::cos(4 * track.phi());
          q4y += std::sin(4 * track.phi());
        }
      }
      psi2Est = std::atan2(q2y, q2x) / 2.;
      psi3Est = std::atan2(q3y, q3x) / 3.;
      psi4Est = std::atan2(q4y, q4x) / 4.;
      v2 = funcV2->Eval(xaxis.centrality);
      v3 = funcV3->Eval(xaxis.centrality);
      v4 = funcV4->Eval(xaxis.centrality);
      densitycorrections.psi2Est = psi2Est;
      densitycorrections.psi3Est = psi3Est;
      densitycorrections.psi4Est = psi4Est;
      densitycorrections.v2 = v2;
      densitycorrections.v3 = v3;
      densitycorrections.v4 = v4;
      densitycorrections.density = tracks.size();
    }
    AcceptedTracks acceptedTracks{0, 0, 0, 0};
    for (const auto& track : tracks) {
      processTrack(track, vtxz, xaxis.multiplicity, run, densitycorrections, acceptedTracks);
      if (cfgConsistentEventFlag & 1)
        if (!acceptedTracks.nPos || !acceptedTracks.nNeg)
          return;
      if (cfgConsistentEventFlag & 2)
        if (acceptedTracks.nFull < 4) // o2-linter: disable=magic-number (at least four tracks in full acceptance)
          return;
      if (cfgConsistentEventFlag & 4)
        if (acceptedTracks.nPos < 2 || acceptedTracks.nNeg < 2) // o2-linter: disable=magic-number (at least two tracks in each subevent)
          return;
      if (cfgConsistentEventFlag & 8)
        if (acceptedTracks.nPos < 2 || acceptedTracks.nMid < 2 || acceptedTracks.nNeg < 2) // o2-linter: disable=magic-number (at least two tracks in all three subevents)
          return;
    }
    if (!cfgFillWeights)
      fillOutputContainers<dt>((cfgTimeDependent) ? xaxis.time : (cfgUseNch) ? xaxis.multiplicity
                                                                             : xaxis.centrality,
                               lRandom, run);
  }

  bool isStable(int pdg)
  {
    if (std::abs(pdg) == PDG_t::kPiPlus)
      return true;
    if (std::abs(pdg) == PDG_t::kKPlus)
      return true;
    if (std::abs(pdg) == PDG_t::kProton)
      return true;
    if (std::abs(pdg) == PDG_t::kElectron)
      return true;
    if (std::abs(pdg) == PDG_t::kMuonMinus)
      return true;
    return false;
  }

  template <typename TTrack>
  void fillAcceptedTracks(TTrack track, AcceptedTracks& acceptedTracks)
  {
    if (posRegionIndex >= 0 && track.eta() > o2::analysis::gfw::regions.GetEtaMin()[posRegionIndex] && track.eta() < o2::analysis::gfw::regions.GetEtaMax()[posRegionIndex])
      ++acceptedTracks.nPos;
    if (negRegionIndex >= 0 && track.eta() > o2::analysis::gfw::regions.GetEtaMin()[negRegionIndex] && track.eta() < o2::analysis::gfw::regions.GetEtaMax()[negRegionIndex])
      ++acceptedTracks.nNeg;
    if (fullRegionIndex >= 0 && track.eta() > o2::analysis::gfw::regions.GetEtaMin()[fullRegionIndex] && track.eta() < o2::analysis::gfw::regions.GetEtaMax()[fullRegionIndex])
      ++acceptedTracks.nFull;
    if (midRegionIndex >= 0 && track.eta() > o2::analysis::gfw::regions.GetEtaMin()[midRegionIndex] && track.eta() < o2::analysis::gfw::regions.GetEtaMax()[midRegionIndex])
      ++acceptedTracks.nMid;
  }

  template <typename TTrack>
  inline void processTrack(TTrack const& track, const float& vtxz, const int& multiplicity, const int& run, DensityCorr densitycorrections, AcceptedTracks& acceptedTracks)
  {
    if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TTrack::all_columns>) {
      if (track.mcParticleId() < 0 || !(track.has_mcParticle()))
        return;

      auto mcParticle = track.mcParticle();
      if (!mcParticle.isPhysicalPrimary())
        return;
      if (!isStable(mcParticle.pdgCode()))
        return;
      if (cfgFillQA) {
        fillTrackQA<kReco, kBefore>(track, vtxz);
        registry.fill(HIST("trackQA/before/nch_pt"), multiplicity, track.pt());
      }
      if (!trackSelected(track))
        return;

      if (cfgFillWeights) {
        fillWeights(track, vtxz, run);
      } else {
        fillPtSums<kReco>(track);
        fillGFW<kReco>(track, vtxz, densitycorrections);
        fillAcceptedTracks(track, acceptedTracks);
      }

      if (cfgFillQA) {
        fillTrackQA<kReco, kAfter>(track, vtxz);
        registry.fill(HIST("trackQA/after/nch_pt"), multiplicity, track.pt());
        if (cfgRunByRun) {
          th1sList[run][hPhi]->Fill(track.phi());
          th1sList[run][hEta]->Fill(track.eta());
        }
      }

    } else if constexpr (framework::has_type_v<aod::mcparticle::McCollisionId, typename TTrack::all_columns>) {
      if (!track.isPhysicalPrimary() || !isStable(track.pdgCode()))
        return;

      fillPtSums<kGen>(track);
      fillGFW<kGen>(track, vtxz, densitycorrections);
      fillAcceptedTracks(track, acceptedTracks);
      if (cfgFillQA) {
        fillTrackQA<kGen, kAfter>(track, vtxz);
        registry.fill(HIST("MCGen/trackQA/nch_pt"), multiplicity, track.pt());
      }
    } else {
      if (cfgFillQA) {
        fillTrackQA<kReco, kBefore>(track, vtxz);
        registry.fill(HIST("trackQA/before/nch_pt"), multiplicity, track.pt());
      }
      if (!trackSelected(track))
        return;

      if (cfgFillWeights) {
        fillWeights(track, vtxz, run);
      } else {
        fillPtSums<kReco>(track);
        fillGFW<kReco>(track, vtxz, densitycorrections);
        fillAcceptedTracks(track, acceptedTracks);
      }
      if (cfgFillQA) {
        fillTrackQA<kReco, kAfter>(track, vtxz);
        registry.fill(HIST("trackQA/after/nch_pt"), multiplicity, track.pt());
        if (cfgRunByRun) {
          th1sList[run][hPhi]->Fill(track.phi());
          th1sList[run][hEta]->Fill(track.eta());
          th3sList[run][hNUAref]->Fill(track.phi(), track.eta(), vtxz, getAcceptance(track, vtxz));
        }
      }
    }
    return;
  }

  template <DataType dt, typename TTrack>
  inline void fillGFW(TTrack track, const double& vtxz, DensityCorr densitycorrections)
  {
    bool withinPtRef = (track.pt() > o2::analysis::gfw::ptreflow && track.pt() < o2::analysis::gfw::ptrefup);
    bool withinPtPOI = (track.pt() > o2::analysis::gfw::ptpoilow && track.pt() < o2::analysis::gfw::ptpoiup);
    if (!withinPtPOI && !withinPtRef)
      return;
    double weff = (dt == kGen) ? 1. : getEfficiency(track);
    if (weff < 0)
      return;
    if (cfgUseDensityDependentCorrection && withinPtRef && dt != kGen) {
      double fphi = densitycorrections.v2 * std::cos(2 * (track.phi() - densitycorrections.psi2Est)) + densitycorrections.v3 * std::cos(3 * (track.phi() - densitycorrections.psi3Est)) + densitycorrections.v4 * std::cos(4 * (track.phi() - densitycorrections.psi4Est));
      fphi = (1 + 2 * fphi);
      int pTBinForEff = hFindPtBin->FindBin(track.pt());
      if (pTBinForEff >= 1 && pTBinForEff <= hFindPtBin->GetNbinsX()) {
        float wEPeff = funcEff[pTBinForEff - 1]->Eval(fphi * densitycorrections.density);
        if (wEPeff > 0.) {
          wEPeff = 1. / wEPeff;
          weff *= wEPeff;
        }
      }
    }
    double wacc = (dt == kGen) ? 1. : getAcceptance(track, vtxz);
    if (withinPtRef)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 1);
    if (withinPtPOI)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 2);
    if (withinPtRef && withinPtPOI)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 4);
    return;
  }

  template <DataType dt, typename TTrack>
  inline void fillPtSums(TTrack track)
  {
    if (track.pt() < o2::analysis::gfw::ptreflow || track.pt() > o2::analysis::gfw::ptrefup)
      return;

    double weff = (dt == kGen) ? 1. : getEfficiency(track);
    if (weff < 0)
      return;

    // Fill pt-pt correlations used in vn-pt correlations (gapped)
    if (std::abs(track.eta()) < cfgEtaPtPt) {
      if (dt == kGen) {
        fFCptgen->fill(1., track.pt());
      } else {
        fFCpt->fill(weff, track.pt());
      }
    }

    // Fill pt-pt correlations for entire eta range
    if (std::abs(track.eta()) < cfgEtaPtPtFull)
      (dt == kGen) ? fFCptgenFull->fill(1., track.pt()) : fFCptFull->fill(weff, track.pt());

    // Fill pt-pt correlations in subevents
    if (track.eta() < -cfgEtaPtPtGap && track.eta() > -cfgEtaPtPtFull)
      (dt == kGen) ? fFCptgenFull->fillSub1(weff, track.pt()) : fFCptFull->fillSub1(weff, track.pt());
    if (track.eta() > cfgEtaPtPtGap && track.eta() < cfgEtaPtPtFull)
      (dt == kGen) ? fFCptgenFull->fillSub2(weff, track.pt()) : fFCptFull->fillSub2(weff, track.pt());
    return;
  }

  template <DataType dt, QAFillTime ft, typename TTrack>
  inline void fillTrackQA(TTrack track, const float vtxz)
  {
    if constexpr (dt == kGen) {
      registry.fill(HIST("MCGen/trackQA/phi_eta_vtxZ"), track.phi(), track.eta(), vtxz);
      registry.fill(HIST("MCGen/trackQA/pt_ref"), track.pt());
      registry.fill(HIST("MCGen/trackQA/pt_poi"), track.pt());
    } else {
      double wacc = getAcceptance(track, vtxz);
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("phi_eta_vtxZ"), track.phi(), track.eta(), vtxz, (ft == kAfter) ? wacc : 1.0);
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("pt_dcaXY_dcaZ"), track.pt(), track.dcaXY(), track.dcaZ());

      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("chi2prTPCcls"), track.tpcChi2NCl());
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("chi2prITScls"), track.itsChi2NCl());
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("nTPCClusters"), track.tpcNClsFound());
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("nITSClusters"), track.itsNCls());
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("nTPCCrossedRows"), track.tpcNClsCrossedRows());

      if (ft == kAfter) {
        registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("pt_ref"), track.pt());
        registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("pt_poi"), track.pt());
      }
    }
  }

  template <typename TCollision>
  float getCentrality(TCollision collision)
  {
    switch (cfgCentEstimator) {
      case kCentFT0C:
        return cfgCentralityFactor * collision.centFT0C();
      case kCentFT0CVariant1:
        return cfgCentralityFactor * collision.centFT0CVariant1();
      case kCentFT0M:
        return cfgCentralityFactor * collision.centFT0M();
      case kCentFV0A:
        return cfgCentralityFactor * collision.centFV0A();
      case kCentNTPV:
        return cfgCentralityFactor * collision.centNTPV();
      case kCentNGlobal:
        return cfgCentralityFactor * collision.centNGlobal();
      case kCentMFT:
        return cfgCentralityFactor * collision.centMFT();
      default:
        return cfgCentralityFactor * collision.centFT0C();
    }
  }

  template <QAFillTime ft, typename TCollision>
  inline void fillEventQA(TCollision collision, XAxis xaxis)
  {
    if constexpr (framework::has_type_v<aod::cent::CentFT0C, typename TCollision::all_columns>) {
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_centT0C"), collision.centFT0C(), xaxis.multiplicity);
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("centT0M_centT0C"), collision.centFT0C(), collision.centFT0M());
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("centV0A_centT0C"), collision.centFT0C(), collision.centFV0A());
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("centGlobal_centT0C"), collision.centFT0C(), collision.centNGlobal());
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("centNTPV_centT0C"), collision.centFT0C(), collision.centNTPV());
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("centMFT_centT0C"), collision.centFT0C(), collision.centMFT());
    }
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_PVTracks"), collision.multNTracksPV(), xaxis.multiplicity);
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_multT0A"), collision.multFT0A(), xaxis.multiplicity);
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_multV0A"), collision.multFV0A(), xaxis.multiplicity);
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
    if (cfgTimeDependent) {
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multiplicity_time"), xaxis.time, xaxis.multiplicity);
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multT0C_time"), xaxis.time, collision.multFT0C());
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multT0A_time"), xaxis.time, collision.multFT0A());
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multV0A_time"), xaxis.time, collision.multFV0A());
      registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multPV_time"), xaxis.time, collision.multNTracksPV());
    }
    return;
  }

  double getTimeSinceStartOfFill(uint64_t timestamp, int firstRun)
  {
    auto runDuration = ccdb->getRunDuration(firstRun);
    uint64_t tsSOF = runDuration.first;
    uint64_t diff = timestamp - tsSOF;
    return static_cast<double>(diff) / 3600000.0;
  }

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs>>::iterator const& collision, aod::BCsWithTimestamps const&, GFWTracks const& tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    if (run != lastRun) {
      lastRun = run;
      LOGF(info, "run = %d", run);
      if (cfgRunByRun) {
        if (std::find(runNumbers.begin(), runNumbers.end(), run) == runNumbers.end()) {
          LOGF(info, "Creating histograms for run %d", run);
          createRunByRunHistograms(run);
          runNumbers.push_back(run);
        } else {
          LOGF(info, "run %d already in runNumbers", run);
        }
        if (!cfgFillWeights)
          loadCorrections(bc);
      }
    }
    if (!cfgFillWeights && !cfgRunByRun)
      loadCorrections(bc);
    registry.fill(HIST("eventQA/eventSel"), kFilteredEvent);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(kFilteredEvent);
    if (!collision.sel8())
      return;
    registry.fill(HIST("eventQA/eventSel"), kSel8);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(kSel8);
    if (cfgDoOccupancySel) {
      int occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < 0 || occupancy > cfgOccupancySelection)
        return;
    }
    registry.fill(HIST("eventQA/eventSel"), kOccupancy);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(kOccupancy);

    const XAxis xaxis{getCentrality(collision), tracks.size(), (cfgTimeDependent) ? getTimeSinceStartOfFill(bc.timestamp(), *firstRunOfCurrentFill) : -1.0};
    if (cfgTimeDependent && run == *firstRunOfCurrentFill && firstRunOfCurrentFill != o2::analysis::gfw::firstRunsOfFill.end() - 1)
      ++firstRunOfCurrentFill;

    if (cfgFillQA)
      fillEventQA<kBefore>(collision, xaxis);
    registry.fill(HIST("eventQA/before/centrality"), xaxis.centrality);
    registry.fill(HIST("eventQA/before/multiplicity"), xaxis.multiplicity);
    if (cfgUseAdditionalEventCut && !eventSelected(collision, xaxis.multiplicity, xaxis.centrality, run))
      return;
    if (cfgFillQA)
      fillEventQA<kAfter>(collision, xaxis);
    processCollision<kReco>(collision, tracks, xaxis, run);
  }
  PROCESS_SWITCH(FlowGfwLightIons, processData, "Process analysis for non-derived data", true);

  void processMCReco(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs>>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>> const& tracks, aod::McParticles const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    if (run != lastRun) {
      lastRun = run;
      if (cfgRunByRun)
        createRunByRunHistograms(run);
    }
    if (!collision.sel8())
      return;

    const XAxis xaxis{getCentrality(collision), tracks.size(), (cfgTimeDependent) ? getTimeSinceStartOfFill(bc.timestamp(), *firstRunOfCurrentFill) : -1.0};
    if (cfgTimeDependent && run == *firstRunOfCurrentFill && firstRunOfCurrentFill != o2::analysis::gfw::firstRunsOfFill.end() - 1)
      ++firstRunOfCurrentFill;

    if (cfgFillQA)
      fillEventQA<kBefore>(collision, xaxis);
    registry.fill(HIST("eventQA/before/centrality"), xaxis.centrality);
    registry.fill(HIST("eventQA/before/multiplicity"), xaxis.multiplicity);
    if (cfgUseAdditionalEventCut && !eventSelected(collision, xaxis.multiplicity, xaxis.centrality, run))
      return;
    if (cfgFillQA)
      fillEventQA<kAfter>(collision, xaxis);
    if (!cfgFillWeights)
      loadCorrections(bc);
    processCollision<kReco>(collision, tracks, xaxis, run);
  }
  PROCESS_SWITCH(FlowGfwLightIons, processMCReco, "Process analysis for MC reconstructed events", false);

  void processMCGen(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs>> const& collisions, aod::McParticles const& particles, GFWTracks const& tracks)
  {
    if (collisions.size() != 1)
      return;
    float centrality = -1;
    for (const auto& collision : collisions) {
      centrality = getCentrality(collision);
    }

    std::vector<int> numberOfTracks;
    for (auto const& collision : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      numberOfTracks.emplace_back(groupedTracks.size());
    }

    const XAxis xaxis{centrality, numberOfTracks[0], -1.0};
    int run = 0;
    processCollision<kGen>(mcCollision, particles, xaxis, run);
  }
  PROCESS_SWITCH(FlowGfwLightIons, processMCGen, "Process analysis for MC generated events", false);

  void processOnTheFly(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    int run = 0;
    registry.fill(HIST("MCGen/impactParameter"), mcCollision.impactParameter(), mcParticles.size());
    const XAxis xaxis{mcCollision.impactParameter(), mcParticles.size(), -1.0};
    processCollision<kGen>(mcCollision, mcParticles, xaxis, run);
  }
  PROCESS_SWITCH(FlowGfwLightIons, processOnTheFly, "Process analysis for MC on-the-fly generated events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGfwLightIons>(cfgc),
  };
}
