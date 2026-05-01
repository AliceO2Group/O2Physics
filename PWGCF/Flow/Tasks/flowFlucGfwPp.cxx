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

/// \file flowFlucGfwPp.cxx
/// \brief GFW task for Event Shape Engineering studies in pp collisions
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch, Wenya Wu, TUM, wenya.wu@cern.ch

#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/FlowPtContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWConfig.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EseTable.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>

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
using namespace o2::analysis::genericframework;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

namespace o2::analysis::gfwflowflucpp
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
} // namespace o2::analysis::gfwflowflucpp

struct FlowFlucGfwPp {
  static constexpr int kInvalidQnBin = -999;
  static constexpr float kInvalidQnSeparator = -999.f;

  static constexpr int kRequireBothEtaSides = 1;
  static constexpr int kRequireFullFourParticleTracks = 2;
  static constexpr int kRequireTwoTracksInBothEtaSides = 4;
  static constexpr int kRequireTwoTracksInThreeEtaRegions = 8;

  static constexpr int kMinTracksForFourParticleCorrelation = 4;
  static constexpr int kMinTracksPerEtaSideForGapCorrelation = 2;
  static constexpr int kMinTracksPerEtaRegionForThreeSubevents = 2;

  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgIsMC, bool, false, "Is MC event")
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FV0A, 4:NTPV, 5:NGlobals, 6:MFT")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Do correlations as function of Nch")
  O2_DEFINE_CONFIGURABLE(cfgFillWeights, bool, false, "Fill NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgRunByRun, bool, false, "Fill histograms on a run-by-run basis")
  O2_DEFINE_CONFIGURABLE(cfgFillFlowRunByRun, bool, false, "Fill flow profile run-by-run (only for v22)")
  O2_DEFINE_CONFIGURABLE(cfgTimeDependent, bool, false, "Fill output as function of time (for contamination studies)")
  O2_DEFINE_CONFIGURABLE(cfgFirstRunsOfFill, std::vector<int>, {}, "First runs of a fill for time dependent analysis")
  O2_DEFINE_CONFIGURABLE(cfgFillQA, bool, false, "Fill QA histograms")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
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
  O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10, "vertex cut (cm)");
  O2_DEFINE_CONFIGURABLE(cfgNoSameBunchPileupCut, bool, false, "kNoSameBunchPileupCut");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodZvtxFT0vsPV, bool, false, "kIsGoodZvtxFT0vsPV");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodITSLayersAll, bool, false, "kIsGoodITSLayersAll");
  O2_DEFINE_CONFIGURABLE(cfgNoCollInTimeRangeStandard, bool, false, "kNoCollInTimeRangeStandard");
  O2_DEFINE_CONFIGURABLE(cfgDoOccupancySel, bool, false, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgOccupancySelection, int, 5000, "Max occupancy selection, -999 to disable");
  O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, false, "Use additional event cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(cfgTVXinTRD, bool, false, "Use kTVXinTRD (reject TRD triggered events)");
  O2_DEFINE_CONFIGURABLE(cfgIsVertexITSTPC, bool, true, "Selects collisions with at least one ITS-TPC track");
  O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field; default CCDB will be queried");
  O2_DEFINE_CONFIGURABLE(cfgFixedMultMin, int, 1, "Minimum for fixed nch range");
  O2_DEFINE_CONFIGURABLE(cfgUseMultiplicityFlowWeights, bool, true, "Enable or disable the use of multiplicity-based event weighting");
  O2_DEFINE_CONFIGURABLE(cfgFixedMultMax, int, 3000, "Maximum for fixed nch range");
  O2_DEFINE_CONFIGURABLE(cfgConsistentEventFlag, int, 0, "Flag to select consistent events - 0: off, 1: v2{2} gap calculable, 2: v2{4} full calculable, 4: v2{4} gap calculable, 8: v2{4} 3sub calculable");
  Configurable<std::vector<double>> cfgMultGlobalCutPars{"cfgMultGlobalCutPars", std::vector<double>{2272.16, -76.6932, 1.01204, -0.00631545, 1.59868e-05, 136.336, -4.97006, 0.121199, -0.0015921, 7.66197e-06}, "Global vs FT0C multiplicity cut parameter values"};
  Configurable<std::vector<double>> cfgMultPVCutPars{"cfgMultPVCutPars", std::vector<double>{3074.43, -106.192, 1.46176, -0.00968364, 2.61923e-05, 182.128, -7.43492, 0.193901, -0.00256715, 1.22594e-05}, "PV vs FT0C multiplicity cut parameter values"};
  Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars", std::vector<double>{-0.223013, 0.715849, 0.664242, 0.0829653, -0.000503733, 1.21185e-06}, "Global vs PV multiplicity cut parameter values"};
  O2_DEFINE_CONFIGURABLE(cfgMultCorrHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
  O2_DEFINE_CONFIGURABLE(cfgMultCorrLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
  O2_DEFINE_CONFIGURABLE(cfgMultGlobalPVCorrCutFunction, std::string, "[0] + [1]*x + 3*([2] + [3]*x + [4]*x*x + [5]*x*x*x)", "Functional for global vs pv multiplicity correlation cut");
  O2_DEFINE_CONFIGURABLE(cfgNumQnBins, int, 10, "Number of qn bins");
  O2_DEFINE_CONFIGURABLE(cfgCentMax, float, 100, "Maximum of centrality or multiplicity");
  O2_DEFINE_CONFIGURABLE(cfgEvtSelCent, bool, true, "Choose event selector as centrality(true) or multicplity(false)");
  O2_DEFINE_CONFIGURABLE(cfgUseNegativeEtaHalfForq2, bool, true, "If true, use -eta half for q2 selection; otherwise use +eta half");
  Configurable<std::vector<float>> qnBinSeparator{"qnBinSeparator", std::vector<float>{-999.f, -999.f, -999.f}, "Qn bin separator"};

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalASideCorrCutFunction, std::string, "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + [10]*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", "Functional for global vs V0A multiplicity low correlation cut");
    Configurable<std::vector<double>> cfgMultGlobalV0ACutPars{"cfgMultGlobalV0ACutPars", std::vector<double>{567.785, 172.715, 0.77888, -0.00693466, 1.40564e-05, 679.853, 66.8068, -0.444332, 0.00115002, -4.92064e-07}, "Global vs FV0A multiplicity cut parameter values"};
    Configurable<std::vector<double>> cfgMultGlobalT0ACutPars{"cfgMultGlobalT0ACutPars", std::vector<double>{241.618, 61.8402, 0.348049, -0.00306078, 6.20357e-06, 315.235, 29.1491, -0.188639, 0.00044528, -9.08912e-08}, "Global vs FT0A multiplicity cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgGlobalV0ALowSigma, float, -3, "Number of sigma deviations below expected value in global vs V0A correlation");
    O2_DEFINE_CONFIGURABLE(cfgGlobalV0AHighSigma, float, 4, "Number of sigma deviations above expected value in global vs V0A correlation");
    O2_DEFINE_CONFIGURABLE(cfgGlobalT0ALowSigma, float, -3., "Number of sigma deviations below expected value in global vs T0A correlation");
    O2_DEFINE_CONFIGURABLE(cfgGlobalT0AHighSigma, float, 4, "Number of sigma deviations above expected value in global vs T0A correlation");
  } cfgGlobalAsideCorrCuts;

  Configurable<GFWBinningCuts> cfgGfwBinning{"cfgGfwBinning",
                                             {40, 16, 72, 300, 0, 3000, 0.2, 10.0, 0.2, 3.0, {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90}},
                                             "Configuration for binning"};

  Configurable<GFWRegions> cfgRegions{"cfgRegions",
                                      {
                                        {"refN", "refP", "refFull"},
                                        {-0.8, 0.4, -0.8},
                                        {-0.4, 0.8, 0.8},
                                        {0, 0, 0}, // pT bins
                                        {1, 1, 1}  // bitmask
                                      },
                                      "Configurations for GFW regions"};

  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig",
                                             {{"refN {2 -2}",
                                               "refN {2 2 -2 -2}",
                                               "refN {2 2 2 -2 -2 -2}",
                                               "refN {2 2 2 2 -2 -2 -2 -2}",
                                               "refP {2 -2}",
                                               "refP {2 2 -2 -2}",
                                               "refP {2 2 2 -2 -2 -2}",
                                               "refP {2 2 2 2 -2 -2 -2 -2}",
                                               "refN {2} refP {-2}",
                                               "refN {2 2} refP {-2 -2}",
                                               "refFull {2 -2}",
                                               "refFull {2 2 -2 -2}",
                                               "refFull {2 2 2 -2 -2 -2}",
                                               "refFull {2 2 2 2 -2 -2 -2 -2}"},
                                              {"ChNeg22",
                                               "ChNeg24",
                                               "ChNeg26",
                                               "ChNeg28",
                                               "ChPos22",
                                               "ChPos24",
                                               "ChPos26",
                                               "ChPos28",
                                               "ChGap22",
                                               "ChGap24",
                                               "ChFull22",
                                               "ChFull24",
                                               "ChFull26",
                                               "ChFull28"},
                                              {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                              {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
                                             "Configurations for pp ESE v2 cumulants"};

  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  struct Config {
    TH1D* mEfficiency = nullptr;
    GFWWeights* mAcceptance;
    bool correctionsLoaded = false;
  } cfg;

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<FlowContainer> fFCgen{FlowContainer("FlowContainer_gen")};
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
  enum EventSelFlags {
    kFilteredEvent = 1,
    kSel8,
    kOccupancy,
    kTVXTRD,
    kNoSamebunchPU,
    kZVtxFT0PV,
    kNoCollTRStd,
    kVtxITSTPC,
    kGoodITSLayers,
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

  std::string getShapeSel() const
  {
    return "ese";
  }

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
  o2::framework::expressions::Filter trackFilter = nabs(aod::track::eta) < cfgEta && aod::track::pt > cfgPtmin&& aod::track::pt < cfgPtmax && (aod::track::itsChi2NCl < cfgChi2PrITSCls) && (aod::track::tpcChi2NCl < cfgChi2PrTPCCls) && nabs(aod::track::dcaZ) < cfgDCAz;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  o2::framework::expressions::Filter mcCollFilter = nabs(aod::mccollision::posZ) < cfgVtxZ;
  o2::framework::expressions::Filter mcParticlesFilter = (aod::mcparticle::eta > o2::analysis::gfwflowflucpp::etalow && aod::mcparticle::eta < o2::analysis::gfwflowflucpp::etaup && aod::mcparticle::pt > o2::analysis::gfwflowflucpp::ptlow && aod::mcparticle::pt < o2::analysis::gfwflowflucpp::ptup);

  using GFWTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;

  void init(InitContext const&)
  {
    LOGF(info, "FlowFlucGfwPp::init()");
    o2::analysis::gfwflowflucpp::regions.SetNames(cfgRegions->GetNames());
    o2::analysis::gfwflowflucpp::regions.SetEtaMin(cfgRegions->GetEtaMin());
    o2::analysis::gfwflowflucpp::regions.SetEtaMax(cfgRegions->GetEtaMax());
    o2::analysis::gfwflowflucpp::regions.SetpTDifs(cfgRegions->GetpTDifs());
    o2::analysis::gfwflowflucpp::regions.SetBitmasks(cfgRegions->GetBitmasks());
    o2::analysis::gfwflowflucpp::configs.SetCorrs(cfgCorrConfig->GetCorrs());
    o2::analysis::gfwflowflucpp::configs.SetHeads(cfgCorrConfig->GetHeads());
    o2::analysis::gfwflowflucpp::configs.SetpTDifs(cfgCorrConfig->GetpTDifs());
    o2::analysis::gfwflowflucpp::configs.SetpTCorrMasks(cfgCorrConfig->GetpTCorrMasks());

    o2::analysis::gfwflowflucpp::regions.Print();
    o2::analysis::gfwflowflucpp::configs.Print();

    o2::analysis::gfwflowflucpp::ptbinning = cfgGfwBinning->GetPtBinning();
    // o2::analysis::gfwflowflucpp::ptpoilow = cfgGfwBinning->GetPtPOImin();
    // o2::analysis::gfwflowflucpp::ptpoiup = cfgGfwBinning->GetPtPOImax();
    o2::analysis::gfwflowflucpp::ptreflow = cfgGfwBinning->GetPtRefMin();
    o2::analysis::gfwflowflucpp::ptrefup = cfgGfwBinning->GetPtRefMax();
    o2::analysis::gfwflowflucpp::ptlow = cfgPtmin;
    o2::analysis::gfwflowflucpp::ptup = cfgPtmax;
    o2::analysis::gfwflowflucpp::etabins = cfgGfwBinning->GetEtaBins();
    o2::analysis::gfwflowflucpp::vtxZbins = cfgGfwBinning->GetVtxZbins();
    o2::analysis::gfwflowflucpp::phibins = cfgGfwBinning->GetPhiBins();
    o2::analysis::gfwflowflucpp::philow = 0.0f;
    o2::analysis::gfwflowflucpp::phiup = o2::constants::math::TwoPI;
    o2::analysis::gfwflowflucpp::nchbins = cfgGfwBinning->GetNchBins();
    o2::analysis::gfwflowflucpp::nchlow = cfgGfwBinning->GetNchMin();
    o2::analysis::gfwflowflucpp::nchup = cfgGfwBinning->GetNchMax();
    o2::analysis::gfwflowflucpp::centbinning = cfgGfwBinning->GetCentBinning();
    cfgGfwBinning->Print();

    o2::analysis::gfwflowflucpp::multGlobalCorrCutPars = cfgMultGlobalCutPars;
    o2::analysis::gfwflowflucpp::multPVCorrCutPars = cfgMultPVCutPars;
    o2::analysis::gfwflowflucpp::multGlobalPVCorrCutPars = cfgMultGlobalPVCutPars;
    o2::analysis::gfwflowflucpp::multGlobalV0ACutPars = cfgGlobalAsideCorrCuts.cfgMultGlobalV0ACutPars;
    o2::analysis::gfwflowflucpp::multGlobalT0ACutPars = cfgGlobalAsideCorrCuts.cfgMultGlobalT0ACutPars;
    o2::analysis::gfwflowflucpp::firstRunsOfFill = cfgFirstRunsOfFill;
    if (cfgTimeDependent && !std::is_sorted(o2::analysis::gfwflowflucpp::firstRunsOfFill.begin(), o2::analysis::gfwflowflucpp::firstRunsOfFill.end())) {
      std::sort(o2::analysis::gfwflowflucpp::firstRunsOfFill.begin(), o2::analysis::gfwflowflucpp::firstRunsOfFill.end());
    }
    firstRunOfCurrentFill = o2::analysis::gfwflowflucpp::firstRunsOfFill.begin();

    AxisSpec phiAxis = {o2::analysis::gfwflowflucpp::phibins, o2::analysis::gfwflowflucpp::philow, o2::analysis::gfwflowflucpp::phiup, "#phi"};
    AxisSpec etaAxis = {o2::analysis::gfwflowflucpp::etabins, -cfgEta, cfgEta, "#eta"};
    AxisSpec vtxAxis = {o2::analysis::gfwflowflucpp::vtxZbins, -cfgVtxZ, cfgVtxZ, "Vtx_{z} (cm)"};
    AxisSpec ptAxis = {o2::analysis::gfwflowflucpp::ptbinning, "#it{p}_{T} GeV/#it{c}"};

    std::string sCentralityEstimator;
    std::map<int, std::string> centEstimatorMap = {
      {kCentFT0C, "FT0C"},
      {kCentFT0CVariant1, "FT0C variant 1"},
      {kCentFT0M, "FT0M"},
      {kCentFV0A, "FV0A"},
      {kCentNTPV, "NTPV"},
      {kCentNGlobal, "NGlobals"},
      {kCentMFT, "MFT"}};
    sCentralityEstimator = centEstimatorMap.at(cfgCentEstimator);
    sCentralityEstimator += " centrality (%)";
    AxisSpec centAxis = {o2::analysis::gfwflowflucpp::centbinning, sCentralityEstimator.c_str()};
    std::vector<double> nchbinning;
    int nchskip = (o2::analysis::gfwflowflucpp::nchup - o2::analysis::gfwflowflucpp::nchlow) / o2::analysis::gfwflowflucpp::nchbins;
    for (int i = 0; i <= o2::analysis::gfwflowflucpp::nchbins; ++i) {
      nchbinning.push_back(nchskip * i + o2::analysis::gfwflowflucpp::nchlow + 0.5);
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

    AxisSpec multAxis = (cfgTimeDependent) ? timeAxis : (cfgUseNch) ? nchAxis
                                                                    : centAxis;

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    int ptbins = o2::analysis::gfwflowflucpp::ptbinning.size() - 1;
    fSecondAxis = (cfgTimeDependent) ? new TAxis(timeAxis.binEdges.size() - 1, &(timeAxis.binEdges[0])) : new TAxis(ptbins, &o2::analysis::gfwflowflucpp::ptbinning[0]);

    if (doprocessq2) {
      registry.add("mq2/eventcounter", "", HistType::kTH1F, {{10, 0, 10}});
      registry.add("mq2/h2_cent_q2_etapos", ";Centrality;#it{q}_{2}^{#eta pos};", HistType::kTH2D, {{100, 0, 100}, {600, 0, 6}});
      registry.add("mq2/h2_cent_q2_etaneg", ";Centrality;#it{q}_{2}^{#eta neg};", HistType::kTH2D, {{100, 0, 100}, {600, 0, 6}});
      registry.add("mq2/h2_mult_q2_etapos", ";Multiplicity;#it{q}_{2}^{#eta pos};", HistType::kTH2D, {{150, 0, 150}, {600, 0, 6}});
      registry.add("mq2/h2_mult_q2_etaneg", ";Multiplicity;#it{q}_{2}^{#eta neg};", HistType::kTH2D, {{150, 0, 150}, {600, 0, 6}});
    }

    if (doprocessData) {
      registry.add("trackQA/before/phi_eta_vtxZ", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registry.add("trackQA/before/pt_dcaXY_dcaZ", "", {HistType::kTH3D, {ptAxis, dcaXYAXis, dcaZAXis}});
      registry.add("trackQA/before/nch_pt", "#it{p}_{T} vs multiplicity; N_{ch}; #it{p}_{T}", {HistType::kTH2D, {nchAxis, ptAxis}});
      registry.add("trackQA/before/chi2prTPCcls", "#chi^{2}/cluster for the TPC track segment; #chi^{2}/TPC cluster", {HistType::kTH1D, {{100, 0., 5.}}});
      registry.add("trackQA/before/chi2prITScls", "#chi^{2}/cluster for the ITS track; #chi^{2}/ITS cluster", {HistType::kTH1D, {{100, 0., 50.}}});
      registry.add("trackQA/before/nTPCClusters", "Number of found TPC clusters; TPC N_{cls}; Counts", {HistType::kTH1D, {{100, 40, 180}}});
      registry.add("trackQA/before/nITSClusters", "Number of found ITS clusters; ITS N_{cls}; Counts", {HistType::kTH1D, {{100, 0, 20}}});
      registry.add("trackQA/before/nTPCCrossedRows", "Number of crossed TPC Rows; TPC X-rows; Counts", {HistType::kTH1D, {{100, 40, 180}}});

      registry.addClone("trackQA/before/", "trackQA/after/");
      registry.add("trackQA/after/pt_ref", "", {HistType::kTH1D, {{100, o2::analysis::gfwflowflucpp::ptreflow, o2::analysis::gfwflowflucpp::ptrefup}}});
      registry.add("trackQA/after/pt_poi", "", {HistType::kTH1D, {{100, o2::analysis::gfwflowflucpp::ptpoilow, o2::analysis::gfwflowflucpp::ptpoiup}}});

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

      if (doprocessData) {
        registry.add("eventQA/before/centrality", "", {HistType::kTH1D, {centAxis}});
        registry.add("eventQA/before/globalTracks_centT0C", "", {HistType::kTH2D, {centAxis, nchAxis}});
        registry.add("eventQA/before/PVTracks_centT0C", "", {HistType::kTH2D, {centAxis, multpvAxis}});
        registry.add("eventQA/before/multT0C_centT0C", "", {HistType::kTH2D, {centAxis, t0cAxis}});

        registry.add("eventQA/before/centT0M_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
        registry.add("eventQA/before/centV0A_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
        registry.add("eventQA/before/centGlobal_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
        registry.add("eventQA/before/centNTPV_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
        registry.add("eventQA/before/centMFT_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});

        if (cfgIsMC) {
          registry.add("MCGen/trackQA/phi_eta_vtxZ", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
          registry.add("MCGen/trackQA/nch_pt", "#it{p}_{T} vs multiplicity; N_{ch}; #it{p}_{T}", {HistType::kTH2D, {nchAxis, ptAxis}});
          registry.add("MCGen/trackQA/pt_ref", "", {HistType::kTH1D, {{100, o2::analysis::gfwflowflucpp::ptreflow, o2::analysis::gfwflowflucpp::ptrefup}}});
          registry.add("MCGen/trackQA/pt_poi", "", {HistType::kTH1D, {{100, o2::analysis::gfwflowflucpp::ptpoilow, o2::analysis::gfwflowflucpp::ptpoiup}}});
        }
      }

      registry.addClone("eventQA/before/", "eventQA/after/");
      registry.add("eventQA/eventSel", "Number of Events;; Counts", {HistType::kTH1D, {{11, 0.5, 11.5}}});
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kFilteredEvent, "Filtered event");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kSel8, "sel8");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kOccupancy, "occupancy");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kTVXTRD, "kTVXinTRD");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoSamebunchPU, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kZVtxFT0PV, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoCollTRStd, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kVtxITSTPC, "kIsVertexITSTPC");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kGoodITSLayers, "kIsGoodITSLayersAll");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kMultCuts, "after Mult cuts");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kTrackCent, "has track + within cent");
      if (!cfgRunByRun && cfgFillWeights) {
        registry.add<TH3>("phi_eta_vtxz_ref", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      }
    }

    if (o2::analysis::gfwflowflucpp::regions.GetSize() < 0)
      LOGF(error, "Configuration contains vectors of different size - check the GFWRegions configurable");
    for (auto i(0); i < o2::analysis::gfwflowflucpp::regions.GetSize(); ++i) {
      fGFW->AddRegion(o2::analysis::gfwflowflucpp::regions.GetNames()[i], o2::analysis::gfwflowflucpp::regions.GetEtaMin()[i], o2::analysis::gfwflowflucpp::regions.GetEtaMax()[i], (o2::analysis::gfwflowflucpp::regions.GetpTDifs()[i]) ? ptbins + 1 : 1, o2::analysis::gfwflowflucpp::regions.GetBitmasks()[i]);
    }
    for (auto i = 0; i < o2::analysis::gfwflowflucpp::configs.GetSize(); ++i) {
      corrconfigs.push_back(fGFW->GetCorrelatorConfig(o2::analysis::gfwflowflucpp::configs.GetCorrs()[i], o2::analysis::gfwflowflucpp::configs.GetHeads()[i], o2::analysis::gfwflowflucpp::configs.GetpTDifs()[i]));
    }
    if (corrconfigs.empty())
      LOGF(error, "Configuration contains vectors of different size - check the GFWCorrConfig configurable");
    fGFW->CreateRegions();
    TObjArray* oba = new TObjArray();
    addConfigObjectsToObjArray(oba, corrconfigs);
    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fSecondAxis);
    fFC->Initialize(oba, multAxis, cfgNbootstrap);

    fFCgen->SetName("FlowContainer_gen");
    fFCgen->SetXAxis(fSecondAxis);
    fFCgen->Initialize(oba, multAxis, cfgNbootstrap);

    delete oba;

    fPtDepDCAxy = new TF1("ptDepDCAxy", Form("[0]*%s", cfgDCAxy->c_str()), 0.001, 100);
    fPtDepDCAxy->SetParameter(0, cfgDCAxyNSigma);
    LOGF(info, "DCAxy pt-dependence function: %s", Form("[0]*%s", cfgDCAxy->c_str()));
    if (cfgUseAdditionalEventCut) {
      fMultPVCutLow = new TF1("fMultPVCutLow", cfgMultCorrLowCutFunction->c_str(), 0, 100);
      fMultPVCutLow->SetParameters(&(o2::analysis::gfwflowflucpp::multPVCorrCutPars[0]));
      fMultPVCutHigh = new TF1("fMultPVCutHigh", cfgMultCorrHighCutFunction->c_str(), 0, 100);
      fMultPVCutHigh->SetParameters(&(o2::analysis::gfwflowflucpp::multPVCorrCutPars[0]));
      fMultCutLow = new TF1("fMultCutLow", cfgMultCorrLowCutFunction->c_str(), 0, 100);
      fMultCutLow->SetParameters(&(o2::analysis::gfwflowflucpp::multGlobalCorrCutPars[0]));
      fMultCutHigh = new TF1("fMultCutHigh", cfgMultCorrHighCutFunction->c_str(), 0, 100);
      fMultCutHigh->SetParameters(&(o2::analysis::gfwflowflucpp::multGlobalCorrCutPars[0]));
      fMultPVGlobalCutHigh = new TF1("fMultPVGlobalCutHigh", cfgMultGlobalPVCorrCutFunction->c_str(), 0, nchbinning.back());
      fMultPVGlobalCutHigh->SetParameters(&(o2::analysis::gfwflowflucpp::multGlobalPVCorrCutPars[0]));

      LOGF(info, "Global V0A function: %s in range 0-%g", cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), v0aAxis.binEdges.back());
      fMultGlobalV0ACutLow = new TF1("fMultGlobalV0ACutLow", cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, v0aAxis.binEdges.back());
      for (std::size_t i = 0; i < o2::analysis::gfwflowflucpp::multGlobalV0ACutPars.size(); ++i)
        fMultGlobalV0ACutLow->SetParameter(i, o2::analysis::gfwflowflucpp::multGlobalV0ACutPars[i]);
      fMultGlobalV0ACutLow->SetParameter(o2::analysis::gfwflowflucpp::multGlobalV0ACutPars.size(), cfgGlobalAsideCorrCuts.cfgGlobalV0ALowSigma);
      for (int i = 0; i < fMultGlobalV0ACutLow->GetNpar(); ++i)
        LOGF(info, "fMultGlobalV0ACutLow par %d = %g", i, fMultGlobalV0ACutLow->GetParameter(i));

      fMultGlobalV0ACutHigh = new TF1("fMultGlobalV0ACutHigh", cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, v0aAxis.binEdges.back());
      for (std::size_t i = 0; i < o2::analysis::gfwflowflucpp::multGlobalV0ACutPars.size(); ++i)
        fMultGlobalV0ACutHigh->SetParameter(i, o2::analysis::gfwflowflucpp::multGlobalV0ACutPars[i]);
      fMultGlobalV0ACutHigh->SetParameter(o2::analysis::gfwflowflucpp::multGlobalV0ACutPars.size(), cfgGlobalAsideCorrCuts.cfgGlobalV0AHighSigma);
      for (int i = 0; i < fMultGlobalV0ACutHigh->GetNpar(); ++i)
        LOGF(info, "fMultGlobalV0ACutHigh par %d = %g", i, fMultGlobalV0ACutHigh->GetParameter(i));

      LOGF(info, "Global T0A function: %s", cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str());
      fMultGlobalT0ACutLow = new TF1("fMultGlobalT0ACutLow", cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, t0aAxis.binEdges.back());
      for (std::size_t i = 0; i < o2::analysis::gfwflowflucpp::multGlobalT0ACutPars.size(); ++i)
        fMultGlobalT0ACutLow->SetParameter(i, o2::analysis::gfwflowflucpp::multGlobalT0ACutPars[i]);
      fMultGlobalT0ACutLow->SetParameter(o2::analysis::gfwflowflucpp::multGlobalT0ACutPars.size(), cfgGlobalAsideCorrCuts.cfgGlobalT0ALowSigma);
      for (int i = 0; i < fMultGlobalT0ACutLow->GetNpar(); ++i)
        LOGF(info, "fMultGlobalT0ACutLow par %d = %g", i, fMultGlobalT0ACutLow->GetParameter(i));

      fMultGlobalT0ACutHigh = new TF1("fMultGlobalT0ACutHigh", cfgGlobalAsideCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, t0aAxis.binEdges.back());
      for (std::size_t i = 0; i < o2::analysis::gfwflowflucpp::multGlobalT0ACutPars.size(); ++i)
        fMultGlobalT0ACutHigh->SetParameter(i, o2::analysis::gfwflowflucpp::multGlobalT0ACutPars[i]);
      fMultGlobalT0ACutHigh->SetParameter(o2::analysis::gfwflowflucpp::multGlobalT0ACutPars.size(), cfgGlobalAsideCorrCuts.cfgGlobalT0AHighSigma);
      for (int i = 0; i < fMultGlobalT0ACutHigh->GetNpar(); ++i)
        LOGF(info, "fMultGlobalT0ACutHigh par %d = %g", i, fMultGlobalT0ACutHigh->GetParameter(i));
    }

    if (cfgConsistentEventFlag) {
      const auto& names = cfgRegions->GetNames();
      auto findRegionIndex = [&](const std::string& name) {
        auto it = std::find(names.begin(), names.end(), name);
        return (it != names.end()) ? std::distance(names.begin(), it) : -1;
      };
      posRegionIndex = findRegionIndex("refP");
      negRegionIndex = findRegionIndex("refN");
      fullRegionIndex = findRegionIndex("refFull");
      midRegionIndex = findRegionIndex("refMid");
    }
  }

  static constexpr std::string_view FillTimeName[] = {"before/", "after/"};

  enum QAFillTime {
    kBefore,
    kAfter
  };

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
    if (cfgTVXinTRD) {
      if (collision.alias_bit(kTVXinTRD)) {
        // TRD triggered
        // "CMTVX-B-NOPF-TRD,minbias_TVX"
        return 0;
      }
      if (cfgFillQA) {
        registry.fill(HIST("eventQA/eventSel"), kTVXTRD);
      }
      if (cfgRunByRun && run != -1)
        th1sList[run][hEventSel]->Fill(kTVXTRD);
    }
    if (cfgNoSameBunchPileupCut) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return 0;
      }
      if (cfgFillQA) {
        registry.fill(HIST("eventQA/eventSel"), kNoSamebunchPU);
      }
      if (cfgRunByRun && run != -1)
        th1sList[run][hEventSel]->Fill(kNoSamebunchPU);
    }
    if (cfgIsGoodZvtxFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return 0;
      }
      if (cfgFillQA) {
        registry.fill(HIST("eventQA/eventSel"), kZVtxFT0PV);
      }
      if (cfgRunByRun && run != -1)
        th1sList[run][hEventSel]->Fill(kZVtxFT0PV);
    }
    if (cfgNoCollInTimeRangeStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      if (cfgFillQA) {
        registry.fill(HIST("eventQA/eventSel"), kNoCollTRStd);
      }
      if (cfgRunByRun && run != -1)
        th1sList[run][hEventSel]->Fill(kNoCollTRStd);
    }

    if (cfgIsVertexITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        // selects collisions with at least one ITS-TPC track, and thus rejects vertices built from ITS-only tracks
        return 0;
      }
      if (cfgFillQA) {
        registry.fill(HIST("eventQA/eventSel"), kVtxITSTPC);
      }
      if (cfgRunByRun && run != -1)
        th1sList[run][hEventSel]->Fill(kVtxITSTPC);
    }

    if (cfgIsGoodITSLayersAll) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        return 0;
      }
      if (cfgFillQA) {
        registry.fill(HIST("eventQA/eventSel"), kGoodITSLayers);
      }
      if (cfgRunByRun && run != -1)
        th1sList[run][hEventSel]->Fill(kGoodITSLayers);
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

    if (vtxz > o2::analysis::gfwflowflucpp::vtxZup || vtxz < o2::analysis::gfwflowflucpp::vtxZlow)
      return 0;

    if (cfgMultCut && cfgUseAdditionalEventCut) {
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
      if (cfgRunByRun && run != -1)
        th1sList[run][hEventSel]->Fill(kMultCuts);
    }
    if (cfgFillQA) {
      registry.fill(HIST("eventQA/eventSel"), kMultCuts);
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
    AxisSpec phiAxis = {o2::analysis::gfwflowflucpp::phibins, o2::analysis::gfwflowflucpp::philow, o2::analysis::gfwflowflucpp::phiup, "#phi"};
    AxisSpec etaAxis = {o2::analysis::gfwflowflucpp::etabins, -cfgEta, cfgEta, "#eta"};
    AxisSpec vtxAxis = {o2::analysis::gfwflowflucpp::vtxZbins, -cfgVtxZ, cfgVtxZ, "Vtx_{z} (cm)"};
    AxisSpec nchAxis = {o2::analysis::gfwflowflucpp::nchbins, o2::analysis::gfwflowflucpp::nchlow, o2::analysis::gfwflowflucpp::nchup, "N_{ch}"};
    AxisSpec centAxis = {o2::analysis::gfwflowflucpp::centbinning, "Centrality (%)"};
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
    histos[hEventSel] = registry.add<TH1>(Form("%d/eventSel", run), "Number of Events;; Counts", {HistType::kTH1D, {{11, 0.5, 11.5}}});
    histos[hEventSel]->GetXaxis()->SetBinLabel(kFilteredEvent, "Filtered event");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kSel8, "sel8");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kOccupancy, "occupancy");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kTVXTRD, "kTVXinTRD");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kNoSamebunchPU, "kNoSameBunchPileup");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kZVtxFT0PV, "kIsGoodZvtxFT0vsPV");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kNoCollTRStd, "kNoCollInTimeRangeStandard");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kVtxITSTPC, "kIsVertexITSTPC");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kGoodITSLayers, "kIsGoodITSLayersAll");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kMultCuts, "after Mult cuts");
    histos[hEventSel]->GetXaxis()->SetBinLabel(kTrackCent, "has track + within cent");
    th1sList.insert(std::make_pair(run, histos));
    std::vector<std::shared_ptr<TH3>> histos3d(kCount_TH3Names);
    histos3d[hNUAref] = registry.add<TH3>(Form("%d/phi_eta_vtxz_ref", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
    th3sList.insert(std::make_pair(run, histos3d));
    return;
  }

  void addConfigObjectsToObjArray(TObjArray* oba, const std::vector<GFW::CorrConfig>& configs)
  {
    const auto shapeSel = getShapeSel();
    for (auto it = configs.begin(); it != configs.end(); ++it) {
      for (int jese = 0; jese < cfgNumQnBins; ++jese) {
        std::string name = Form("%s_%d_%s", shapeSel.c_str(), jese, it->Head.c_str());
        std::string title = it->Head + std::string("_ese");
        oba->Add(new TNamed(name.c_str(), title.c_str()));
      }
    }
  }

  template <DataType dt>
  void fillOutputContainers(const float& centmult, const double& rndm, const int& qPtmp, const int& run = 0)
  {
    const auto shapeSel = getShapeSel();

    for (uint l_ind = 0; l_ind < corrconfigs.size(); ++l_ind) {
      auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), 0, kTRUE).real();
      if (dnx == 0)
        continue;

      auto val = fGFW->Calculate(corrconfigs.at(l_ind), 0, kFALSE).real() / dnx;
      if (std::abs(val) >= 1.)
        continue;

      std::string profName = Form("%s_%d_%s", shapeSel.c_str(), qPtmp, corrconfigs.at(l_ind).Head.c_str());

      if constexpr (dt == kGen) {
        fFCgen->FillProfile(profName.c_str(),
                            centmult,
                            val,
                            cfgUseMultiplicityFlowWeights ? dnx : 1.0,
                            rndm);
      } else {
        fFC->FillProfile(profName.c_str(),
                         centmult,
                         val,
                         cfgUseMultiplicityFlowWeights ? dnx : 1.0,
                         rndm);

        if (cfgRunByRun && cfgFillFlowRunByRun && l_ind == 0) {
          tpfsList[run][pfCorr22]->Fill(centmult,
                                        val,
                                        cfgUseMultiplicityFlowWeights ? dnx : 1.0);
        }
      }
    }
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

  template <typename T>
  float computeqnVec(T const& col, bool useNegativeEtaHalf)
  {
    if (col.qvecTPCposReVec().empty() || col.qvecTPCposImVec().empty() ||
        col.qvecTPCnegReVec().empty() || col.qvecTPCnegImVec().empty()) {
      return -1.f;
    }

    if (col.nTrkTPCpos() <= 0 || col.nTrkTPCneg() <= 0)
      return -1.f;

    const auto qvecPos =
      std::sqrt(col.qvecTPCposReVec()[0] * col.qvecTPCposReVec()[0] +
                col.qvecTPCposImVec()[0] * col.qvecTPCposImVec()[0]) *
      std::sqrt(col.nTrkTPCpos());

    const auto qvecNeg =
      std::sqrt(col.qvecTPCnegReVec()[0] * col.qvecTPCnegReVec()[0] +
                col.qvecTPCnegImVec()[0] * col.qvecTPCnegImVec()[0]) *
      std::sqrt(col.nTrkTPCneg());

    return useNegativeEtaHalf ? qvecNeg : qvecPos;
  }

  /// \return the 1-d qn-vector separator to 2-d
  std::vector<std::vector<float>> getQnBinSeparator2D(std::vector<float> flat, const int numQnBins = 10)
  {
    size_t nBins = numQnBins + 1;

    if (flat.empty() || flat.size() % nBins != 0) {
      LOGP(error, "ConfQnBinSeparator size = {} is not divisible by {}",
           flat.size(), nBins);
      return {{kInvalidQnSeparator, kInvalidQnSeparator}};
    }

    size_t nCent = flat.size() / nBins;
    std::vector<std::vector<float>> res(nCent, std::vector<float>(nBins));

    for (size_t i = 0; i < nCent; ++i) {
      for (size_t j = 0; j < nBins; ++j) {
        res[i][j] = flat[i * nBins + j];
      }
    }
    return res;
  }

  /// Get the bin number of qn-vector(FT0C) of an event
  /// \param centBinWidth centrality bin width, example: per 1%, per 10% ...
  /// \return bin number of qn-vector of the event
  // add a param : bool doFillHisto ?
  int myqnBin(float centrality, float centMax, float qn, std::vector<float> qnBinSprt, const int numQnBins, float centBinWidth = 1.f)
  {
    auto twoDSeparator = getQnBinSeparator2D(qnBinSprt, numQnBins);
    if (twoDSeparator.empty() || twoDSeparator[0][0] == kInvalidQnSeparator) {
      LOGP(warning, "ConfQnBinSeparator not set, using default fallback!");
      return kInvalidQnBin; // safe fallback
    }

    int qnBin = kInvalidQnBin;
    int mycentBin = static_cast<int>(centrality / centBinWidth);
    if (mycentBin >= static_cast<int>(centMax / centBinWidth))
      return qnBin;

    if (mycentBin > static_cast<int>(twoDSeparator.size()) - 1)
      return qnBin;

    for (int iqn(0); iqn < static_cast<int>(twoDSeparator[mycentBin].size()) - 1; ++iqn) {
      if (qn > twoDSeparator[mycentBin][iqn] && qn <= twoDSeparator[mycentBin][iqn + 1]) {
        qnBin = iqn;
        break;
      } else {
        continue;
      }
    }

    return qnBin;
  }

  template <DataType dt, typename TCollision, typename TTracks>
  void processCollision(TCollision collision, TTracks tracks, const XAxis& xaxis, const int& run, const int& qPtmp)
  {
    if (tracks.size() < 1)
      return;
    if (dt != kGen && xaxis.centrality >= 0 &&
        (xaxis.centrality < o2::analysis::gfwflowflucpp::centbinning.front() ||
         xaxis.centrality > o2::analysis::gfwflowflucpp::centbinning.back()))
      return;
    if (xaxis.multiplicity < cfgFixedMultMin || xaxis.multiplicity > cfgFixedMultMax)
      return;

    if (dt != kGen) {
      if (cfgFillQA) {
        registry.fill(HIST("eventQA/eventSel"), kTrackCent);
      }
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kTrackCent);
    }

    if (cfgFillQA && xaxis.centrality >= 0)
      registry.fill(HIST("eventQA/after/centrality"), xaxis.centrality);
    if (cfgFillQA)
      registry.fill(HIST("eventQA/after/multiplicity"), xaxis.multiplicity);

    float vtxz = collision.posZ();
    if (dt != kGen && cfgRunByRun) {
      th1sList[run][hVtxZ]->Fill(vtxz);
      th1sList[run][hMult]->Fill(xaxis.multiplicity);
      th1sList[run][hCent]->Fill(xaxis.centrality);
    }

    fGFW->Clear();
    float lRandom = fRndm->Rndm();

    AcceptedTracks acceptedTracks{0, 0, 0, 0};

    for (const auto& track : tracks) {
      processTrack(track, vtxz, xaxis.multiplicity, run, acceptedTracks);
    }

    if (cfgConsistentEventFlag & kRequireBothEtaSides)
      if (!acceptedTracks.nPos || !acceptedTracks.nNeg)
        return;
    if (cfgConsistentEventFlag & kRequireFullFourParticleTracks)
      if (acceptedTracks.nFull < kMinTracksForFourParticleCorrelation)
        return;
    if (cfgConsistentEventFlag & kRequireTwoTracksInBothEtaSides)
      if (acceptedTracks.nPos < kMinTracksPerEtaSideForGapCorrelation ||
          acceptedTracks.nNeg < kMinTracksPerEtaSideForGapCorrelation)
        return;
    if (cfgConsistentEventFlag & kRequireTwoTracksInThreeEtaRegions)
      if (acceptedTracks.nPos < kMinTracksPerEtaRegionForThreeSubevents ||
          acceptedTracks.nMid < kMinTracksPerEtaRegionForThreeSubevents ||
          acceptedTracks.nNeg < kMinTracksPerEtaRegionForThreeSubevents)
        return;

    fillOutputContainers<dt>(cfgUseNch ? static_cast<float>(xaxis.multiplicity) : xaxis.centrality,
                             lRandom, qPtmp, run);
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
    if (posRegionIndex >= 0 && track.eta() > o2::analysis::gfwflowflucpp::regions.GetEtaMin()[posRegionIndex] && track.eta() < o2::analysis::gfwflowflucpp::regions.GetEtaMax()[posRegionIndex])
      ++acceptedTracks.nPos;
    if (negRegionIndex >= 0 && track.eta() > o2::analysis::gfwflowflucpp::regions.GetEtaMin()[negRegionIndex] && track.eta() < o2::analysis::gfwflowflucpp::regions.GetEtaMax()[negRegionIndex])
      ++acceptedTracks.nNeg;
    if (fullRegionIndex >= 0 && track.eta() > o2::analysis::gfwflowflucpp::regions.GetEtaMin()[fullRegionIndex] && track.eta() < o2::analysis::gfwflowflucpp::regions.GetEtaMax()[fullRegionIndex])
      ++acceptedTracks.nFull;
    if (midRegionIndex >= 0 && track.eta() > o2::analysis::gfwflowflucpp::regions.GetEtaMin()[midRegionIndex] && track.eta() < o2::analysis::gfwflowflucpp::regions.GetEtaMax()[midRegionIndex])
      ++acceptedTracks.nMid;
  }

  template <typename TTrack>
  inline void processTrack(TTrack const& track, const float& vtxz, const int& multiplicity, const int& run, AcceptedTracks& acceptedTracks)
  {
    if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TTrack::all_columns>) {
      if (track.mcParticleId() < 0 || !(track.has_mcParticle()))
        return;

      auto mcParticle = track.mcParticle();
      if (!mcParticle.isPhysicalPrimary())
        return;
      if (!isStable(mcParticle.pdgCode()))
        return;
      if (cfgFillQA && cfgIsMC) {
        fillTrackQA<kReco, kBefore>(track, vtxz);
        registry.fill(HIST("trackQA/before/nch_pt"), multiplicity, track.pt());
      }
      if (!trackSelected(track))
        return;

      if (cfgFillWeights) {
        fillWeights(track, vtxz, run);
      } else {
        fillGFW<kReco>(track, vtxz);
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

      fillGFW<kGen>(track, vtxz);
      fillAcceptedTracks(track, acceptedTracks);
      if (cfgFillQA && cfgIsMC) {
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
        fillGFW<kReco>(track, vtxz);
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
  inline void fillGFW(TTrack track, const double& vtxz)
  {
    bool withinPtRef = (track.pt() > o2::analysis::gfwflowflucpp::ptreflow && track.pt() < o2::analysis::gfwflowflucpp::ptrefup);
    bool withinPtPOI = (track.pt() > o2::analysis::gfwflowflucpp::ptpoilow && track.pt() < o2::analysis::gfwflowflucpp::ptpoiup);
    if (!withinPtPOI && !withinPtRef)
      return;
    double weff = (dt == kGen) ? 1. : getEfficiency(track);
    if (weff < 0)
      return;

    double wacc = (dt == kGen) ? 1. : getAcceptance(track, vtxz);
    if (withinPtRef)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 1);
    if (withinPtPOI)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 2);
    if (withinPtRef && withinPtPOI)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 4);
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
        return collision.centFT0C();
      case kCentFT0CVariant1:
        return collision.centFT0CVariant1();
      case kCentFT0M:
        return collision.centFT0M();
      case kCentFV0A:
        return collision.centFV0A();
      case kCentNTPV:
        return collision.centNTPV();
      case kCentNGlobal:
        return collision.centNGlobal();
      case kCentMFT:
        return collision.centMFT();
      default:
        return collision.centFT0C();
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

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                                           aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms,
                                           aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals,
                                           aod::CentMFTs, aod::Qvectors,
                                           aod::QvectorTPCposVecs, aod::QvectorTPCnegVecs>>::iterator const& collision,
                   aod::BCsWithTimestamps const&, GFWTracks const& tracks)
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

    const XAxis xaxis{
      getCentrality(collision),
      tracks.size(),
      (cfgTimeDependent) ? getTimeSinceStartOfFill(bc.timestamp(), *firstRunOfCurrentFill) : -1.0};

    if (cfgTimeDependent && run == *firstRunOfCurrentFill &&
        firstRunOfCurrentFill != o2::analysis::gfwflowflucpp::firstRunsOfFill.end() - 1)
      ++firstRunOfCurrentFill;

    if (cfgFillQA)
      fillEventQA<kBefore>(collision, xaxis);

    registry.fill(HIST("eventQA/before/centrality"), xaxis.centrality);
    registry.fill(HIST("eventQA/before/multiplicity"), xaxis.multiplicity);

    if (!eventSelected(collision, xaxis.multiplicity, xaxis.centrality, run))
      return;

    if (cfgFillQA)
      fillEventQA<kAfter>(collision, xaxis);

    float qn = computeqnVec(collision, cfgUseNegativeEtaHalfForq2);
    if (qn < 0)
      return;

    int qPtmp = myqnBin(cfgEvtSelCent ? xaxis.centrality : xaxis.multiplicity,
                        cfgCentMax, qn, qnBinSeparator, cfgNumQnBins);
    if (qPtmp < 0)
      return;

    processCollision<kReco>(collision, tracks, xaxis, run, qPtmp);
  }
  PROCESS_SWITCH(FlowFlucGfwPp, processData, "Process analysis for non-derived data", false);

  void processq2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs, aod::Qvectors, aod::QvectorTPCposVecs, aod::QvectorTPCnegVecs>>::iterator const& collision, aod::BCsWithTimestamps const&, GFWTracks const& tracks)
  {
    float count{0.5};
    registry.fill(HIST("mq2/eventcounter"), count++);
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("mq2/eventcounter"), count++);
    if (cfgDoOccupancySel) {
      int occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < 0 || occupancy > cfgOccupancySelection) {
        return;
      }
    }
    registry.fill(HIST("mq2/eventcounter"), count++);

    const XAxis xaxis{getCentrality(collision), tracks.size(), -1.0};
    if (!eventSelected(collision, xaxis.multiplicity, xaxis.centrality, -1))
      return;

    const auto centr = xaxis.centrality;
    const auto multi = xaxis.multiplicity;
    const auto qvecPos = computeqnVec(collision, false);
    const auto qvecNeg = computeqnVec(collision, true);

    if (!std::isfinite(qvecPos) || !std::isfinite(qvecNeg) || qvecPos < 0 || qvecNeg < 0) {
      return;
    }

    registry.fill(HIST("mq2/eventcounter"), count++);
    registry.fill(HIST("mq2/h2_cent_q2_etapos"), centr, qvecPos);
    registry.fill(HIST("mq2/h2_cent_q2_etaneg"), centr, qvecNeg);
    registry.fill(HIST("mq2/h2_mult_q2_etapos"), multi, qvecPos);
    registry.fill(HIST("mq2/h2_mult_q2_etaneg"), multi, qvecNeg);
  }
  PROCESS_SWITCH(FlowFlucGfwPp, processq2, "Process analysis for filling q-vectors", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowFlucGfwPp>(cfgc),
  };
}
