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

/// \file flowGenericFramework.cxx
/// \brief Task to analyse angular and transverse momentum correlations with GFW
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
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
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
} // namespace o2::analysis::gfw

struct FlowGenericFramework {

  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgMpar, int, 8, "Highest order of pt-pt correlations")
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Do correlations as function of Nch")
  O2_DEFINE_CONFIGURABLE(cfgFillWeights, bool, false, "Fill NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgRunByRun, bool, false, "Fill histograms on a run-by-run basis")
  O2_DEFINE_CONFIGURABLE(cfgFillQA, bool, false, "Fill QA histograms")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgUseCentralMoments, bool, true, "Use central moments in vn-pt calculations")
  O2_DEFINE_CONFIGURABLE(cfgUsePID, bool, true, "Enable PID information")
  O2_DEFINE_CONFIGURABLE(cfgUseGapMethod, bool, false, "Use gap method in vn-pt calculations")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgDCAxyNSigma, float, 7, "Cut on number of sigma deviations from expected DCA in the transverse direction");
  O2_DEFINE_CONFIGURABLE(cfgDCAz, float, 2, "Cut on DCA in the longitudinal direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgNTPCCls, float, 70, "Cut on number of TPC clusters found");
  O2_DEFINE_CONFIGURABLE(cfgNTPCXrows, float, 70, "Cut on number of TPC crossed rows");
  O2_DEFINE_CONFIGURABLE(cfgMinNITSCls, float, 5, "Cut on minimum number of ITS clusters found");
  O2_DEFINE_CONFIGURABLE(cfgChi2PrITSCls, float, 36, "Cut on chi^2 per ITS clusters found");
  O2_DEFINE_CONFIGURABLE(cfgChi2PrTPCCls, float, 2.5, "Cut on chi^2 per TPC clusters found");
  O2_DEFINE_CONFIGURABLE(cfgPtmin, float, 0.2, "minimum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgPtmax, float, 10, "maximum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgEta, float, 0.8, "eta cut");
  O2_DEFINE_CONFIGURABLE(cfgEtaPtPt, float, 0.4, "eta cut for pt-pt correlations");
  O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10, "vertex cut (cm)");
  O2_DEFINE_CONFIGURABLE(cfgOccupancySelection, int, 2000, "Max occupancy selection, -999 to disable");
  O2_DEFINE_CONFIGURABLE(cfgNoSameBunchPileupCut, bool, true, "kNoSameBunchPileupCut");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodZvtxFT0vsPV, bool, true, "kIsGoodZvtxFT0vsPV");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodITSLayersAll, bool, true, "kIsGoodITSLayersAll");
  O2_DEFINE_CONFIGURABLE(cfgNoCollInTimeRangeStandard, bool, true, "kNoCollInTimeRangeStandard");
  O2_DEFINE_CONFIGURABLE(cfgDoOccupancySel, bool, true, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, true, "Use additional event cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(cfgTVXinTRD, bool, true, "Use kTVXinTRD (reject TRD triggered events)");
  O2_DEFINE_CONFIGURABLE(cfgIsVertexITSTPC, bool, true, "Selects collisions with at least one ITS-TPC track");
  O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field; default CCDB will be queried");
  O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.5, "pt cut on TOF for PID");
  O2_DEFINE_CONFIGURABLE(cfgUseDensityDependentCorrection, bool, false, "Use density dependent efficiency correction based on Run 2 measurements");
  Configurable<std::vector<double>> cfgTrackDensityP0{"cfgTrackDensityP0", std::vector<double>{0.7217476707, 0.7384792571, 0.7542625668, 0.7640680200, 0.7701951667, 0.7755299053, 0.7805901710, 0.7849446786, 0.7957356586, 0.8113039262, 0.8211968966, 0.8280558878, 0.8329342135}, "parameter 0 for track density efficiency correction"};
  Configurable<std::vector<double>> cfgTrackDensityP1{"cfgTrackDensityP1", std::vector<double>{-2.169488e-05, -2.191913e-05, -2.295484e-05, -2.556538e-05, -2.754463e-05, -2.816832e-05, -2.846502e-05, -2.843857e-05, -2.705974e-05, -2.477018e-05, -2.321730e-05, -2.203315e-05, -2.109474e-05}, "parameter 1 for track density efficiency correction"};
  Configurable<std::vector<double>> cfgMultGlobalCutPars{"cfgMultGlobalCutPars", std::vector<double>{2272.16, -76.6932, 1.01204, -0.00631545, 1.59868e-05, 136.336, -4.97006, 0.121199, -0.0015921, 7.66197e-06}, "Global multiplicity cut parameter values"};
  Configurable<std::vector<double>> cfgMultPVCutPars{"cfgMultPVCutPars", std::vector<double>{3074.43, -106.192, 1.46176, -0.00968364, 2.61923e-05, 182.128, -7.43492, 0.193901, -0.00256715, 1.22594e-05}, "PV multiplicity cut parameter values"};
  O2_DEFINE_CONFIGURABLE(cfgMultCorrHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
  O2_DEFINE_CONFIGURABLE(cfgMultCorrLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");

  Configurable<GFWBinningCuts> cfgGFWBinning{"cfgGFWBinning", {40, 16, 72, 300, 0, 3000, 0.2, 10.0, 0.2, 3.0, {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90}}, "Configuration for binning"};
  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN", "refP", "refFull"}, {-0.8, 0.4, -0.8}, {-0.4, 0.8, 0.8}, {0, 0, 0}, {1, 1, 1}}, "Configurations for GFW regions"};

  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}", "refFull {2 -2}", "refFull {2 2 -2 -2}"}, {"ChGap22", "ChGap32", "ChGap42", "ChFull22", "ChFull24"}, {0, 0, 0, 0, 0}, {15, 1, 1, 0, 0}}, "Configurations for each correlation to calculate"};

  // #include "PWGCF/TwoParticleCorrelations/TableProducer/Productions/skimmingconf_20221115.cxx" // NOLINT
  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  struct Config {
    TH1D* mEfficiency = nullptr;
    std::vector<GFWWeights*> mAcceptance;
    bool correctionsLoaded = false;
  } cfg;

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<FlowPtContainer> fFCpt{FlowPtContainer("FlowPtContainer")};
  OutputObj<FlowContainer> fFCgen{FlowContainer("FlowContainer_gen")};
  HistogramRegistry registry{"registry"};

  // QA outputs
  std::map<int, std::vector<std::shared_ptr<TH1>>> th1sList;
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
  // NUA outputs
  enum OutputTH3Names {
    hNUAref = 0,
    hNUAch,
    hNUApi,
    hNUAka,
    hNUApr,
    kCount_TH3Names
  };
  enum CentEstimators {
    kCentFT0C = 0,
    kCentFT0CVariant1,
    kCentFT0M,
    kCentFV0A,
    kCentNTPV
  };

  // Define global variables
  // Generic Framework
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;

  TRandom3* fRndm = new TRandom3(0);
  TAxis* fPtAxis;
  int lastRun = -1;
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

  // Event selection cuts - Alex
  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(InitContext const&)
  {
    LOGF(info, "flowGenericFramework::init()");
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
    o2::analysis::gfw::multGlobalCorrCutPars = cfgMultGlobalCutPars;
    o2::analysis::gfw::multPVCorrCutPars = cfgMultPVCutPars;

    AxisSpec phiAxis = {o2::analysis::gfw::phibins, o2::analysis::gfw::philow, o2::analysis::gfw::phiup, "#phi"};
    AxisSpec etaAxis = {o2::analysis::gfw::etabins, -cfgEta, cfgEta, "#eta"};
    AxisSpec vtxAxis = {o2::analysis::gfw::vtxZbins, -cfgVtxZ, cfgVtxZ, "Vtx_{z} (cm)"};
    AxisSpec ptAxis = {o2::analysis::gfw::ptbinning, "#it{p}_{T} GeV/#it{c}"};
    std::string sCentralityEstimator;
    switch (cfgCentEstimator) {
      case kCentFT0C:
        sCentralityEstimator = "FT0C";
        break;
      case kCentFT0CVariant1:
        sCentralityEstimator = "FT0C variant 1";
        break;
      case kCentFT0M:
        sCentralityEstimator = "FT0M";
        break;
      case kCentFV0A:
        sCentralityEstimator = "FV0A";
        break;
      case kCentNTPV:
        sCentralityEstimator = "NTPV";
        break;
      default:
        sCentralityEstimator = "FT0C";
    }
    sCentralityEstimator += " centrality (%)";
    AxisSpec centAxis = {o2::analysis::gfw::centbinning, sCentralityEstimator.c_str()};
    std::vector<double> nchbinning;
    int nchskip = (o2::analysis::gfw::nchup - o2::analysis::gfw::nchlow) / o2::analysis::gfw::nchbins;
    for (int i = 0; i <= o2::analysis::gfw::nchbins; ++i) {
      nchbinning.push_back(nchskip * i + o2::analysis::gfw::nchlow + 0.5);
    }
    AxisSpec nchAxis = {nchbinning, "N_{ch}"};
    AxisSpec bAxis = {200, 0, 20, "#it{b}"};
    AxisSpec t0cAxis = {70, 0, 70000, "N_{ch} (T0C)"};
    AxisSpec t0aAxis = {200, 0, 200, "N_{ch}"};
    AxisSpec multpvAxis = {4000, 0, 4000, "N_{ch} (PV)"};
    AxisSpec multAxis = (doprocessOnTheFly && !cfgUseNch) ? bAxis : (cfgUseNch) ? nchAxis
                                                                                : centAxis;
    AxisSpec dcaZAXis = {200, -2, 2, "DCA_{z} (cm)"};
    AxisSpec dcaXYAXis = {200, -1, 1, "DCA_{xy} (cm)"};
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    int ptbins = o2::analysis::gfw::ptbinning.size() - 1;
    fPtAxis = new TAxis(ptbins, &o2::analysis::gfw::ptbinning[0]);

    if (doprocessMCGen || doprocessOnTheFly) {
      registry.add("MCGen/before/pt_gen", "", {HistType::kTH1D, {ptAxis}});
      registry.add("MCGen/before/phi_eta_vtxZ_gen", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registry.addClone("MCGen/before/", "MCGen/after/");
      if (doprocessOnTheFly)
        registry.add("MCGen/impactParameter", "", {HistType::kTH2D, {{bAxis, nchAxis}}});
    }
    if (doprocessMCReco || doprocessData || doprocessRun2) {
      registry.add("trackQA/before/phi_eta_vtxZ", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registry.add("trackQA/before/pt_dcaXY_dcaZ", "", {HistType::kTH3D, {ptAxis, dcaXYAXis, dcaZAXis}});
      registry.add("trackQA/before/chi2prTPCcls", "#chi^{2}/cluster for the TPC track segment", {HistType::kTH1D, {{100, 0., 5.}}});
      registry.add("trackQA/before/chi2prITScls", "#chi^{2}/cluster for the ITS track", {HistType::kTH1D, {{100, 0., 50.}}});
      registry.add("trackQA/before/nTPCClusters", "Number of found TPC clusters", {HistType::kTH1D, {{100, 40, 180}}});
      registry.add("trackQA/before/nITSClusters", "Number of found ITS clusters", {HistType::kTH1D, {{100, 0, 20}}});
      registry.add("trackQA/before/nTPCCrossedRows", "Number of crossed TPC Rows", {HistType::kTH1D, {{100, 40, 180}}});

      registry.addClone("trackQA/before/", "trackQA/after/");
      registry.add("trackQA/after/pt_ref", "", {HistType::kTH1D, {{100, o2::analysis::gfw::ptreflow, o2::analysis::gfw::ptrefup}}});
      registry.add("trackQA/after/pt_poi", "", {HistType::kTH1D, {{100, o2::analysis::gfw::ptpoilow, o2::analysis::gfw::ptpoiup}}});

      registry.add("eventQA/before/centrality", "", {HistType::kTH1D, {centAxis}});
      registry.add("eventQA/before/multiplicity", "", {HistType::kTH1D, {nchAxis}});
      registry.add("eventQA/before/globalTracks_centT0C", "", {HistType::kTH2D, {centAxis, nchAxis}});
      registry.add("eventQA/before/PVTracks_centT0C", "", {HistType::kTH2D, {centAxis, multpvAxis}});
      registry.add("eventQA/before/globalTracks_PVTracks", "", {HistType::kTH2D, {multpvAxis, nchAxis}});
      registry.add("eventQA/before/globalTracks_multT0A", "", {HistType::kTH2D, {t0aAxis, nchAxis}});
      registry.add("eventQA/before/globalTracks_multV0A", "", {HistType::kTH2D, {t0aAxis, nchAxis}});
      registry.add("eventQA/before/multV0A_multT0A", "", {HistType::kTH2D, {t0aAxis, t0aAxis}});
      registry.add("eventQA/before/multT0C_centT0C", "", {HistType::kTH2D, {centAxis, t0cAxis}});
      registry.addClone("eventQA/before/", "eventQA/after/");
      registry.add("eventQA/eventSel", "Number of Events;; Counts", {HistType::kTH1D, {{11, 0, 11}}});
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(1, "Filtered event");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(2, "sel8");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(3, "occupancy");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(4, "kTVXinTRD");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(5, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(6, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(7, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(8, "kIsVertexITSTPC");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(9, "kIsGoodITSLayersAll");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(10, "after Mult cuts");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(11, "has track + within cent");

      if (!cfgRunByRun) {
        if (cfgUsePID) {
          registry.add<TH3>("phi_eta_vtxz_ref", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
          registry.add<TH3>("phi_eta_vtxz_ch", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
          registry.add<TH3>("phi_eta_vtxz_pi", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
          registry.add<TH3>("phi_eta_vtxz_ka", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
          registry.add<TH3>("phi_eta_vtxz_pr", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
        } else {
          registry.add<TH3>("phi_eta_vtxz_ref", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
        }
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
    if (doprocessData || doprocessRun2 || doprocessMCReco) {
      fFC->SetName("FlowContainer");
      fFC->SetXAxis(fPtAxis);
      fFC->Initialize(oba, multAxis, cfgNbootstrap);
    }
    if (doprocessMCGen || doprocessOnTheFly) {
      fFCgen->SetName("FlowContainer_gen");
      fFCgen->SetXAxis(fPtAxis);
      fFCgen->Initialize(oba, multAxis, cfgNbootstrap);
    }
    delete oba;
    fFCpt->setUseCentralMoments(cfgUseCentralMoments);
    fFCpt->setUseGapMethod(cfgUseGapMethod);
    fFCpt->initialise(multAxis, cfgMpar, o2::analysis::gfw::configs, cfgNbootstrap);
    // Event selection - Alex
    if (cfgUseAdditionalEventCut) {
      fMultPVCutLow = new TF1("fMultPVCutLow", cfgMultCorrLowCutFunction->c_str(), 0, 100);
      fMultPVCutLow->SetParameters(&(o2::analysis::gfw::multPVCorrCutPars[0]));
      fMultPVCutHigh = new TF1("fMultPVCutHigh", cfgMultCorrHighCutFunction->c_str(), 0, 100);
      fMultPVCutHigh->SetParameters(&(o2::analysis::gfw::multPVCorrCutPars[0]));
      fMultCutLow = new TF1("fMultCutLow", cfgMultCorrLowCutFunction->c_str(), 0, 100);
      fMultCutLow->SetParameters(&(o2::analysis::gfw::multGlobalCorrCutPars[0]));
      fMultCutHigh = new TF1("fMultCutHigh", cfgMultCorrHighCutFunction->c_str(), 0, 100);
      fMultCutHigh->SetParameters(&(o2::analysis::gfw::multGlobalCorrCutPars[0]));
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
        for (auto i = 0; i < fPtAxis->GetNbins(); ++i) {
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
    if (!cfgRunByRun && cfg.correctionsLoaded)
      return;
    if (!cfgAcceptance.value.empty()) {
      std::string runstr = (cfgRunByRun) ? "RunByRun/" : "";
      cfg.mAcceptance.clear();
      if (cfgUsePID) {
        cfg.mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + runstr + "ref/", timestamp));
        cfg.mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + runstr + "ch/", timestamp));
        cfg.mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + runstr + "pi/", timestamp));
        cfg.mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + runstr + "ka/", timestamp));
        cfg.mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + runstr + "pr/", timestamp));
      } else {
        cfg.mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + runstr, timestamp));
      }
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
  double getAcceptance(TTrack track, const double& vtxz, int index)
  { // 0 ref, 1 ch, 2 pi, 3 ka, 4 pr
    double wacc = 1;
    if (!cfg.mAcceptance.empty())
      wacc = cfg.mAcceptance[index]->getNUA(track.phi(), track.eta(), vtxz);
    return wacc;
  }

  template <typename TTrack>
  double getEfficiency(TTrack track)
  { //-1 ref, 0 ch, 1 pi, 2 ka, 3 pr
    double eff = 1.;
    if (cfg.mEfficiency)
      eff = cfg.mEfficiency->GetBinContent(cfg.mEfficiency->FindBin(track.pt()));
    if (eff == 0)
      return -1.;
    else
      return 1. / eff;
  }

  template <typename TTrack>
  int getNsigmaPID(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaCombined = {std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi()), std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr())};
    int pid = -1;
    float nsigma = 3.0;

    // Choose which nSigma to use
    std::array<float, 3> nSigmaToUse = (track.pt() > cfgTofPtCut && track.hasTOF()) ? nSigmaCombined : nSigmaTPC;
    if (track.pt() >= cfgTofPtCut && !track.hasTOF())
      return -1;

    // Select particle with the lowest nsigma
    const int nspecies = 3;
    for (int i = 0; i < nspecies; ++i) {
      if (std::abs(nSigmaToUse[i]) < nsigma) {
        pid = i;
        nsigma = std::abs(nSigmaToUse[i]);
      }
    }
    return pid + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
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
      registry.fill(HIST("eventQA/eventSel"), 3.5);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(3.5);
    }

    if (cfgNoSameBunchPileupCut) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), 4.5);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(4.5);
    }
    if (cfgIsGoodZvtxFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), 5.5);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(5.5);
    }
    if (cfgNoCollInTimeRangeStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), 6.5);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(6.5);
    }

    if (cfgIsVertexITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        // selects collisions with at least one ITS-TPC track, and thus rejects vertices built from ITS-only tracks
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), 7.5);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(7.5);
    }

    if (cfgIsGoodITSLayersAll) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), 8.5);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(8.5);
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
    // auto multV0A = collision.multFV0A();
    // auto multT0A = collision.multFT0A();
    // auto multT0C = collision.multFT0C();
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
      registry.fill(HIST("eventQA/eventSel"), 9.5);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(9.5);
    }
    return 1;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    if (cfgDCAxyNSigma && (std::fabs(track.dcaXY()) > cfgDCAxyNSigma / 7. * (0.0105f + 0.0035f / track.pt())))
      return false;
    return ((track.tpcNClsCrossedRows() >= cfgNTPCXrows) && (track.tpcNClsFound() >= cfgNTPCCls) && (track.itsNCls() >= cfgMinNITSCls));
  }
  enum DataType {
    kReco,
    kGen
  };

  template <typename TTrack>
  void fillWeights(const TTrack track, const double vtxz, const int& pid_index, const int& run)
  {
    if (cfgUsePID) {
      double ptpidmins[] = {o2::analysis::gfw::ptpoilow, o2::analysis::gfw::ptpoilow, 0.3, 0.5};                  // min pt for ch, pi, ka, pr
      double ptpidmaxs[] = {o2::analysis::gfw::ptpoiup, o2::analysis::gfw::ptpoiup, 6.0, 6.0};                    // max pt for ch, pi, ka, pr
      bool withinPtPOI = (ptpidmins[pid_index] < track.pt()) && (track.pt() < ptpidmaxs[pid_index]);              // within POI pT range
      bool withinPtRef = (o2::analysis::gfw::ptreflow < track.pt()) && (track.pt() < o2::analysis::gfw::ptrefup); // within RF pT range
      if (cfgRunByRun) {
        if (withinPtRef && !pid_index)
          th3sList[run][hNUAref]->Fill(track.phi(), track.eta(), vtxz); // pt-subset of charged particles for ref flow
        if (withinPtPOI)
          th3sList[run][hNUAch + pid_index]->Fill(track.phi(), track.eta(), vtxz); // charged and id'ed particle weights
      } else {
        if (withinPtRef && !pid_index)
          registry.fill(HIST("phi_eta_vtxz_ref"), track.phi(), track.eta(), vtxz); // pt-subset of charged particles for ref flow
        if (withinPtPOI) {
          switch (pid_index) {
            case 0:
              registry.fill(HIST("phi_eta_vtxz_ch"), track.phi(), track.eta(), vtxz); // charged particle weights
              break;
            case 1:
              registry.fill(HIST("phi_eta_vtxz_pi"), track.phi(), track.eta(), vtxz); // pion weights
              break;
            case 2:
              registry.fill(HIST("phi_eta_vtxz_ka"), track.phi(), track.eta(), vtxz); // kaon weights
              break;
            case 3:
              registry.fill(HIST("phi_eta_vtxz_pr"), track.phi(), track.eta(), vtxz); // proton weights
              break;
          }
        }
      }
    } else {
      if (cfgRunByRun)
        th3sList[run][hNUAref]->Fill(track.phi(), track.eta(), vtxz);
      else
        registry.fill(HIST("phi_eta_vtxz_ref"), track.phi(), track.eta(), vtxz);
    }
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
    histos[hEventSel] = registry.add<TH1>(Form("%d/eventSel", run), "Number of Events;; Counts", {HistType::kTH1D, {{11, 0, 11}}});
    histos[hEventSel]->GetXaxis()->SetBinLabel(1, "Filtered event");
    histos[hEventSel]->GetXaxis()->SetBinLabel(2, "sel8");
    histos[hEventSel]->GetXaxis()->SetBinLabel(3, "occupancy");
    histos[hEventSel]->GetXaxis()->SetBinLabel(4, "kTVXinTRD");
    histos[hEventSel]->GetXaxis()->SetBinLabel(5, "kNoSameBunchPileup");
    histos[hEventSel]->GetXaxis()->SetBinLabel(6, "kIsGoodZvtxFT0vsPV");
    histos[hEventSel]->GetXaxis()->SetBinLabel(7, "kNoCollInTimeRangeStandard");
    histos[hEventSel]->GetXaxis()->SetBinLabel(8, "kIsVertexITSTPC");
    histos[hEventSel]->GetXaxis()->SetBinLabel(9, "kIsGoodITSLayersAll");
    histos[hEventSel]->GetXaxis()->SetBinLabel(10, "after Mult cuts");
    histos[hEventSel]->GetXaxis()->SetBinLabel(11, "has track + within cent");
    th1sList.insert(std::make_pair(run, histos));
    std::vector<std::shared_ptr<TH3>> histos3d(kCount_TH3Names);
    if (cfgUsePID) {
      histos3d[hNUAref] = registry.add<TH3>(Form("%d/phi_eta_vtxz_ref", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      histos3d[hNUAch] = registry.add<TH3>(Form("%d/phi_eta_vtxz_ch", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      histos3d[hNUApi] = registry.add<TH3>(Form("%d/phi_eta_vtxz_pi", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      histos3d[hNUAka] = registry.add<TH3>(Form("%d/phi_eta_vtxz_ka", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      histos3d[hNUApr] = registry.add<TH3>(Form("%d/phi_eta_vtxz_pr", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
    } else {
      histos3d[hNUAref] = registry.add<TH3>(Form("%d/phi_eta_vtxz_ref", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
    }
    th3sList.insert(std::make_pair(run, histos3d));
    return;
  }

  template <DataType dt>
  void fillOutputContainers(const float& centmult, const double& rndm)
  {
    fFCpt->calculateCorrelations();
    fFCpt->fillPtProfiles(centmult, rndm);
    fFCpt->fillCMProfiles(centmult, rndm);
    if (!cfgUseGapMethod)
      fFCpt->fillVnPtStdProfiles(centmult, rndm);
    for (uint l_ind = 0; l_ind < corrconfigs.size(); ++l_ind) {
      if (!corrconfigs.at(l_ind).pTDif) {
        auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), 0, kTRUE).real();
        if (dnx == 0)
          continue;
        auto val = fGFW->Calculate(corrconfigs.at(l_ind), 0, kFALSE).real() / dnx;
        if (std::abs(val) < 1) {
          (dt == kGen) ? fFCgen->FillProfile(corrconfigs.at(l_ind).Head.c_str(), centmult, val, dnx, rndm) : fFC->FillProfile(corrconfigs.at(l_ind).Head.c_str(), centmult, val, dnx, rndm);
          if (cfgUseGapMethod)
            fFCpt->fillVnPtProfiles(centmult, val, dnx, rndm, o2::analysis::gfw::configs.GetpTCorrMasks()[l_ind]);
        }
        continue;
      }
      for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
        auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), i - 1, kTRUE).real();
        if (dnx == 0)
          continue;
        auto val = fGFW->Calculate(corrconfigs.at(l_ind), i - 1, kFALSE).real() / dnx;
        if (std::abs(val) < 1)
          (dt == kGen) ? fFCgen->FillProfile(Form("%s_pt_%i", corrconfigs.at(l_ind).Head.c_str(), i), centmult, val, dnx, rndm) : fFC->FillProfile(Form("%s_pt_%i", corrconfigs.at(l_ind).Head.c_str(), i), centmult, val, dnx, rndm);
      }
    }
    return;
  }

  template <DataType dt, typename TCollision, typename TTracks>
  void processCollision(TCollision collision, TTracks tracks, const float& centrality, const int& run)
  {
    if (tracks.size() < 1)
      return;
    if (dt != kGen && (centrality < o2::analysis::gfw::centbinning.front() || centrality > o2::analysis::gfw::centbinning.back()))
      return;
    if (dt != kGen) {
      registry.fill(HIST("eventQA/eventSel"), 10.5);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(10.5);
    }
    float vtxz = collision.posZ();
    if (dt != kGen && cfgRunByRun) {
      th1sList[run][hVtxZ]->Fill(vtxz);
      th1sList[run][hMult]->Fill(tracks.size());
      th1sList[run][hCent]->Fill(centrality);
    }
    fGFW->Clear();
    fFCpt->clearVector();
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
      v2 = funcV2->Eval(centrality);
      v3 = funcV3->Eval(centrality);
      v4 = funcV4->Eval(centrality);
      densitycorrections.psi2Est = psi2Est;
      densitycorrections.psi3Est = psi3Est;
      densitycorrections.psi4Est = psi4Est;
      densitycorrections.v2 = v2;
      densitycorrections.v3 = v3;
      densitycorrections.v4 = v4;
      densitycorrections.density = tracks.size();
    }

    for (const auto& track : tracks) {
      processTrack(track, vtxz, run, densitycorrections);
    }
    if (!cfgFillWeights)
      fillOutputContainers<dt>((cfgUseNch) ? tracks.size() : centrality, lRandom);
  }

  template <typename TTrack>
  inline void processTrack(TTrack const& track, const float& vtxz, const int& run, DensityCorr densitycorrections)
  {
    if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TTrack::all_columns>) {
      if (track.mcParticleId() < 0 || !(track.has_mcParticle()))
        return;

      auto mcParticle = track.mcParticle();
      if (!mcParticle.isPhysicalPrimary())
        return;
      if (cfgFillQA)
        fillTrackQA<kReco, kBefore>(track, vtxz);

      if (mcParticle.eta() < o2::analysis::gfw::etalow || mcParticle.eta() > o2::analysis::gfw::etaup || mcParticle.pt() < o2::analysis::gfw::ptlow || mcParticle.pt() > o2::analysis::gfw::ptup)
        return;

      if (!trackSelected(track))
        return;

      int pidIndex = 0;
      if (cfgUsePID) {
        if (std::abs(mcParticle.pdgCode()) == kPiPlus)
          pidIndex = 1;
        if (std::abs(mcParticle.pdgCode()) == kKPlus)
          pidIndex = 2;
        if (std::abs(mcParticle.pdgCode()) == kProton)
          pidIndex = 3;
      }

      if (cfgFillWeights) {
        fillWeights(mcParticle, vtxz, 0, run);
      } else {
        fillPtSums<kReco>(track, vtxz);
        fillGFW<kReco>(mcParticle, vtxz, pidIndex, densitycorrections);
      }

      if (cfgFillQA) {
        fillTrackQA<kReco, kAfter>(track, vtxz);
        if (cfgRunByRun) {
          th1sList[run][hPhi]->Fill(track.phi());
          th1sList[run][hEta]->Fill(track.eta());
        }
      }

    } else if constexpr (framework::has_type_v<aod::mcparticle::McCollisionId, typename TTrack::all_columns>) {
      if (!track.isPhysicalPrimary())
        return;
      if (cfgFillQA)
        fillTrackQA<kGen, kBefore>(track, vtxz);

      if (track.eta() < o2::analysis::gfw::etalow || track.eta() > o2::analysis::gfw::etaup || track.pt() < o2::analysis::gfw::ptlow || track.pt() > o2::analysis::gfw::ptup)
        return;

      int pidIndex = 0;
      if (cfgUsePID) {
        if (std::abs(track.pdgCode()) == kPiPlus)
          pidIndex = 1;
        if (std::abs(track.pdgCode()) == kKPlus)
          pidIndex = 2;
        if (std::abs(track.pdgCode()) == kProton)
          pidIndex = 3;
      }

      fillPtSums<kGen>(track, vtxz);
      fillGFW<kGen>(track, vtxz, pidIndex, densitycorrections);

      if (cfgFillQA)
        fillTrackQA<kGen, kAfter>(track, vtxz);
    } else {
      if (cfgFillQA)
        fillTrackQA<kReco, kBefore>(track, vtxz);

      if (!trackSelected(track))
        return;

      int pidIndex = 0;
      if (cfgUsePID) {
        // pid_index = getBayesPIDIndex(track);
        pidIndex = getNsigmaPID(track);
      }
      if (cfgFillWeights) {
        fillWeights(track, vtxz, pidIndex, run);
      } else {
        fillPtSums<kReco>(track, vtxz);
        fillGFW<kReco>(track, vtxz, pidIndex, densitycorrections);
      }
      if (cfgFillQA) {
        fillTrackQA<kReco, kAfter>(track, vtxz);
        if (cfgRunByRun) {
          th1sList[run][hPhi]->Fill(track.phi());
          th1sList[run][hEta]->Fill(track.eta());
        }
      }
    }
  }

  template <DataType dt, typename TTrack>
  inline void fillPtSums(TTrack track, const double& vtxz)
  {
    double wacc = (dt == kGen) ? 1. : getAcceptance(track, vtxz, 0);
    double weff = (dt == kGen) ? 1. : getEfficiency(track);
    if (weff < 0)
      return;
    if (std::abs(track.eta()) < cfgEtaPtPt) {
      fFCpt->fill(weff, track.pt());
    }
    if (!cfgUseGapMethod) {
      std::complex<double> q2p = {weff * wacc * std::cos(2 * track.phi()), weff * wacc * std::sin(2 * track.phi())};
      std::complex<double> q2n = {weff * wacc * std::cos(-2 * track.phi()), weff * wacc * std::sin(-2 * track.phi())};
      fFCpt->fillArray(q2p, q2n, weff * track.pt(), weff);
      fFCpt->fillArray(weff * wacc, weff * wacc, weff, weff);
    }
  }

  template <DataType dt, typename TTrack>
  inline void fillGFW(TTrack track, const double& vtxz, int pid_index, DensityCorr densitycorrections)
  {
    if (cfgUsePID) { // Analysing POI flow with id'ed particles
      double ptmins[] = {o2::analysis::gfw::ptpoilow, o2::analysis::gfw::ptpoilow, 0.3, 0.5};
      double ptmaxs[] = {o2::analysis::gfw::ptpoiup, o2::analysis::gfw::ptpoiup, 6.0, 6.0};
      bool withinPtRef = (track.pt() > o2::analysis::gfw::ptreflow && track.pt() < o2::analysis::gfw::ptrefup);
      bool withinPtPOI = (track.pt() > ptmins[pid_index] && track.pt() < ptmaxs[pid_index]);
      bool withinPtNch = (track.pt() > ptmins[0] && track.pt() < ptmaxs[0]);
      if (!withinPtPOI && !withinPtRef)
        return;
      double waccRef = (dt == kGen) ? 1. : getAcceptance(track, vtxz, 0);
      double waccPOI = (dt == kGen) ? 1. : withinPtPOI ? getAcceptance(track, vtxz, pid_index + 1)
                                                       : getAcceptance(track, vtxz, 0); //
      if (withinPtRef && withinPtPOI && pid_index)
        waccRef = waccPOI; // if particle is both (then it's overlap), override ref with POI
      if (withinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), waccRef, 1);
      if (withinPtPOI && pid_index)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), waccPOI, (1 << (pid_index + 1)));
      if (withinPtNch)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), waccPOI, 2);
      if (withinPtPOI && withinPtRef && pid_index)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), waccPOI, (1 << (pid_index + 5)));
      if (withinPtNch && withinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), waccPOI, 32);
    } else { // Analysing only integrated flow
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
      double wacc = (dt == kGen) ? 1. : getAcceptance(track, vtxz, 0);
      if (withinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 1);
      if (withinPtPOI)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 2);
      if (withinPtRef && withinPtPOI)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 4);
    }
    return;
  }

  template <DataType dt, QAFillTime ft, typename TTrack>
  inline void fillTrackQA(TTrack track, const float vtxz)
  {
    if constexpr (dt == kGen) {
      registry.fill(HIST("MCGen/") + HIST(FillTimeName[ft]) + HIST("phi_eta_vtxZ_gen"), track.phi(), track.eta(), vtxz);
      registry.fill(HIST("MCGen/") + HIST(FillTimeName[ft]) + HIST("pt_gen"), track.pt());
    } else {
      double wacc = getAcceptance(track, vtxz, 0);
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

  template <QAFillTime ft, typename CollisionObject, typename TracksObject>
  inline void fillEventQA(CollisionObject collision, TracksObject tracks)
  {
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_centT0C"), collision.centFT0C(), tracks.size());
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_multT0A"), collision.multFT0A(), tracks.size());
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_multV0A"), collision.multFV0A(), tracks.size());
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
    registry.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
    return;
  }

  o2::framework::expressions::Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZ;
  o2::framework::expressions::Filter trackFilter = nabs(aod::track::eta) < cfgEta && aod::track::pt > cfgPtmin&& aod::track::pt < cfgPtmax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::itsChi2NCl < cfgChi2PrITSCls) && (aod::track::tpcChi2NCl < cfgChi2PrTPCCls) && nabs(aod::track::dcaZ) < cfgDCAz;

  using GFWTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFKa, aod::pidTPCKa, aod::pidTOFPr, aod::pidTPCPr>>;

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs>>::iterator const& collision, aod::BCsWithTimestamps const&, GFWTracks const& tracks)
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
    registry.fill(HIST("eventQA/eventSel"), 0.5);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(0.5);
    if (!collision.sel8())
      return;
    registry.fill(HIST("eventQA/eventSel"), 1.5);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(1.5);
    if (cfgDoOccupancySel) {
      int occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < 0 || occupancy > cfgOccupancySelection)
        return;
    }
    registry.fill(HIST("eventQA/eventSel"), 2.5);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(2.5);
    float centrality;
    switch (cfgCentEstimator) {
      case kCentFT0C:
        centrality = collision.centFT0C();
        break;
      case kCentFT0CVariant1:
        centrality = collision.centFT0CVariant1();
        break;
      case kCentFT0M:
        centrality = collision.centFT0M();
        break;
      case kCentFV0A:
        centrality = collision.centFV0A();
        break;
      case kCentNTPV:
        centrality = collision.centNTPV();
        break;
      default:
        centrality = collision.centFT0C();
    }
    if (cfgFillQA)
      fillEventQA<kBefore>(collision, tracks);
    registry.fill(HIST("eventQA/before/centrality"), centrality);
    registry.fill(HIST("eventQA/before/multiplicity"), tracks.size());
    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), centrality, run))
      return;
    if (cfgFillQA)
      fillEventQA<kAfter>(collision, tracks);
    registry.fill(HIST("eventQA/after/centrality"), centrality);
    registry.fill(HIST("eventQA/after/multiplicity"), tracks.size());
    processCollision<kReco>(collision, tracks, centrality, run);
  }
  PROCESS_SWITCH(FlowGenericFramework, processData, "Process analysis for non-derived data", true);

  void processMCReco(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs>>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>> const& tracks, aod::McParticles const&)
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
    const auto centrality = collision.centFT0C();
    if (cfgFillQA)
      fillEventQA<kBefore>(collision, tracks);
    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), centrality, run))
      return;
    if (cfgFillQA)
      fillEventQA<kAfter>(collision, tracks);

    if (!cfgFillWeights)
      loadCorrections(bc);
    processCollision<kReco>(collision, tracks, centrality, run);
  }
  PROCESS_SWITCH(FlowGenericFramework, processMCReco, "Process analysis for MC reconstructed events", false);

  o2::framework::expressions::Filter mcCollFilter = nabs(aod::mccollision::posZ) < cfgVtxZ;
  void processMCGen(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs>> const& collisions, aod::McParticles const& particles)
  {
    if (collisions.size() != 1)
      return;
    float centrality = -1;
    for (const auto& collision : collisions) {
      centrality = collision.centFT0C();
    }
    int run = 0;
    processCollision<kGen>(mcCollision, particles, centrality, run);
  }
  PROCESS_SWITCH(FlowGenericFramework, processMCGen, "Process analysis for MC generated events", false);

  void processOnTheFly(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    int run = 0;
    registry.fill(HIST("MCGen/impactParameter"), mcCollision.impactParameter(), mcParticles.size());
    processCollision<kGen>(mcCollision, mcParticles, mcCollision.impactParameter(), run);
  }
  PROCESS_SWITCH(FlowGenericFramework, processOnTheFly, "Process analysis for MC on-the-fly generated events", false);

  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentRun2V0Ms>>::iterator const& collision, aod::BCsWithTimestamps const&, GFWTracks const& tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    if (run != lastRun) {
      lastRun = run;
      if (cfgRunByRun)
        createRunByRunHistograms(run);
    }
    if (!collision.sel7())
      return;
    const auto centrality = collision.centRun2V0M();
    if (!cfgFillWeights)
      loadCorrections(bc);
    processCollision<kReco>(collision, tracks, centrality, run);
  }
  PROCESS_SWITCH(FlowGenericFramework, processRun2, "Process analysis for Run 2 converted data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGenericFramework>(cfgc),
  };
}
