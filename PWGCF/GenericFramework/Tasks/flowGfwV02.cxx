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
//
/// \file flowGfwV02.cxx
/// \brief Skeleton copy of flowGfwLightIons with empty function bodies
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch

#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWConfig.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"
#include "PWGCF/JCorran/DataModel/JCatalyst.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <TF1.h>
#include <TH1.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TString.h>

#include <boost/algorithm/string/find.hpp>
#include <sys/types.h>

#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace analysis::genericframework;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, (DEFAULT), (HELP)}; // NOLINT(bugprone-macro-parentheses)
static constexpr std::array<std::array<float, 3>, 20> LongArrayFloat = {{{{1.1, 2.1, 3.1}}, {{1.2, 2.2, 3.2}}, {{1.3, 2.3, 3.3}}, {{-1.1, -2.1, -3.1}}, {{-1.2, -2.2, -3.2}}, {{-1.3, -2.3, -3.3}}, {{1.1, 1.1, 1.1}}, {{1.2, 1.2, 1.2}}, {{1.3, 1.3, 1.3}}, {{-1.1, -1.1, -1.1}}, {{-1.2, -1.2, -1.2}}, {{-1.3, -1.3, -1.3}}, {{1.1, 1.1, 1.1}}, {{1.2, 1.2, 1.2}}, {{1.3, 1.3, 1.3}}, {{-1.1, -1.1, -1.1}}, {{-1.2, -1.2, -1.2}}, {{-1.3, -1.3, -1.3}}, {{1.1, 1.1, 1.1}}, {{1.2, 1.2, 1.2}}}};

template <typename T>
auto projectMatrix(Array2D<T> const& mat, std::array<float, 6>& array1, std::array<float, 6>& array2, std::array<float, 6>& array3)
{
  for (auto j = 0; j < static_cast<int>(mat.rows); ++j) {
    array1[j] = mat(j, 0);
    array2[j] = mat(j, 1);
    array3[j] = mat(j, 2);
  }
}

struct FlowGfwV02 {

  struct GFWMemberCache {
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
    float nchup = 300;
    std::vector<double> centbinning{90};
    int nBootstrap = 10;
    GFWRegions regions;
    GFWCorrConfigs configs;
    std::vector<double> multGlobalCorrCutPars;
    std::vector<double> multPVCorrCutPars;
    std::vector<double> multGlobalPVCorrCutPars;
  } gfwMemberCache;

  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgMpar, int, 4, "Highest order of pt-pt correlations")
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
  O2_DEFINE_CONFIGURABLE(cfgCentralityFactor, double, 1., "Correction factor for testing centrality robustness");
  O2_DEFINE_CONFIGURABLE(cfgFillQA, bool, false, "Fill QA histograms")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgPIDEfficiency, bool, false, "Use PID efficiency for efficiency calculation")
  O2_DEFINE_CONFIGURABLE(cfgFixedMultMin, int, 1, "Minimum for fixed nch range");
  O2_DEFINE_CONFIGURABLE(cfgFixedMultMax, int, 3000, "Maximum for fixed nch range");
  O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.5f, "Minimum pt to use TOF N-sigma")
  O2_DEFINE_CONFIGURABLE(cfgUseItsPID, bool, true, "Use ITS PID for particle identification")
  O2_DEFINE_CONFIGURABLE(cfgGetNsigmaQA, bool, true, "Get QA histograms for selection of pions, kaons, and protons")
  O2_DEFINE_CONFIGURABLE(cfgGetdEdx, bool, true, "Get dEdx histograms for pions, kaons, and protons")
  O2_DEFINE_CONFIGURABLE(cfgUseMultiplicityFlowWeights, bool, true, "Enable or disable the use of multiplicity-based event weighting");
  O2_DEFINE_CONFIGURABLE(cfgUseMultiplicityFracWeights, bool, false, "Enable or disable the use of multiplicity-based event weighting for pT fraction");
  O2_DEFINE_CONFIGURABLE(cfgNormalizeByCharged, bool, true, "Enable or disable the normalization by charged particles");
  O2_DEFINE_CONFIGURABLE(cfgConsistentEventFlag, int, 15, "Flag for consistent event selection");
  O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, true, "Use additional event cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(cfgUseV0, bool, false, "Use V0 analysis");

  // Event selection cuts
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

  // Event selection cuts
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgEtaSubAMin, float, -0.8, "Minimum eta for subevent A");
    O2_DEFINE_CONFIGURABLE(cfgEtaSubAMax, float, -0.5, "Maximum eta for subevent A");
    O2_DEFINE_CONFIGURABLE(cfgEtaSubBMin, float, 0.5, "Minimum eta for subevent B");
    O2_DEFINE_CONFIGURABLE(cfgEtaSubBMax, float, 0.8, "Maximum eta for subevent B");
    O2_DEFINE_CONFIGURABLE(cfgEtaSubCMin, float, -0.4, "Minimum eta for subevent C");
    O2_DEFINE_CONFIGURABLE(cfgEtaSubCMax, float, 0.4, "Maximum eta for subevent C");
  } cfgSubeventCuts;

  struct : ConfigurableGroup {
    Configurable<std::vector<double>> cfgMultGlobalCutPars{"cfgMultGlobalCutPars", std::vector<double>{2272.16, -76.6932, 1.01204, -0.00631545, 1.59868e-05, 136.336, -4.97006, 0.121199, -0.0015921, 7.66197e-06}, "Global vs FT0C multiplicity cut parameter values"};
    Configurable<std::vector<double>> cfgMultPVCutPars{"cfgMultPVCutPars", std::vector<double>{3074.43, -106.192, 1.46176, -0.00968364, 2.61923e-05, 182.128, -7.43492, 0.193901, -0.00256715, 1.22594e-05}, "PV vs FT0C multiplicity cut parameter values"};
    Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars", std::vector<double>{-0.223013, 0.715849, 0.664242, 0.0829653, -0.000503733, 1.21185e-06}, "Global vs PV multiplicity cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultCorrHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultCorrLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalPVCorrCutFunction, std::string, "[0] + [1]*x + 3*([2] + [3]*x + [4]*x*x + [5]*x*x*x)", "Functional for global vs pv multiplicity correlation cut");
  } cfgMultCorrCuts;

  // Track selection cuts
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgDCAxyNSigma, float, 7, "Cut on number of sigma deviations from expected DCA in the transverse direction");
    O2_DEFINE_CONFIGURABLE(cfgDCAxy, std::string, "(0.0026+0.005/(x^1.01))", "Functional form of pt-dependent DCAxy cut");
    O2_DEFINE_CONFIGURABLE(cfgDCAz, float, 2, "Cut on DCA in the longitudinal direction (cm)");
    O2_DEFINE_CONFIGURABLE(cfgNTPCCls, float, 50, "Cut on number of TPC clusters found");
    O2_DEFINE_CONFIGURABLE(cfgNTPCXrows, float, 70, "Cut on number of TPC crossed rows");
    O2_DEFINE_CONFIGURABLE(cfgMinNITSCls, float, 5, "Cut on minimum number of ITS clusters found");
    O2_DEFINE_CONFIGURABLE(cfgChi2PrITSCls, float, 36, "Cut on chi^2 per ITS clusters found");
    O2_DEFINE_CONFIGURABLE(cfgChi2PrTPCCls, float, 2.5, "Cut on chi^2 per TPC clusters found");
    O2_DEFINE_CONFIGURABLE(cfgPtMin, float, 0.2, "minimum pt (GeV/c)");
    O2_DEFINE_CONFIGURABLE(cfgPtMax, float, 10, "maximum pt (GeV/c)");
    O2_DEFINE_CONFIGURABLE(cfgEtaMax, float, 0.8, "eta cut");
    O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10, "vertex cut (cm)");
    O2_DEFINE_CONFIGURABLE(cfgOccupancySelection, int, 2000, "Max occupancy selection, -999 to disable");
  } cfgTrackCuts;

  struct : ConfigurableGroup {
    ConfigurableAxis cfgAxisPt{"cfgAxisPt", {VARIABLE_WIDTH, 0.2f, 0.5f, 1.f, 1.5f, 2.f, 3.f, 4.f, 6.f, 10.f}, "pt axis for histograms"};
    ConfigurableAxis cfgAxisNsigmaTPC{"cfgAxisNsigmaTPC", {80, -5.f, 5.f}, "nsigmaTPC axis"};
    ConfigurableAxis cfgAxisNsigmaTOF{"cfgAxisNsigmaTOF", {80, -5.f, 5.f}, "nsigmaTOF axis"};
    ConfigurableAxis cfgAxisNsigmaITS{"cfgAxisNsigmaITS", {80, -5.f, 5.f}, "nsigmaITS axis"};
    ConfigurableAxis cfgAxisTpcSignal{"cfgAxisTpcSignal", {250, 0.f, 250.f}, "dEdx axis for TPC"};
    ConfigurableAxis cfgAxisSigma{"cfgAxisSigma", {200, 0, 20}, "sigma axis for TPC"};
  } cfgPIDHistograms;

  // GFW binning
  Configurable<GFWBinningCuts> cfgGFWBinning{"cfgGFWBinning", {40, 16, 72, 300, 0, 3000, 0.2, 10.0, 0.2, 5.0, {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90}}, "Configuration for binning"};
  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN", "refP", "refFull", "refMid", "piP", "kaP", "prP"}, {-0.8, 0.5, -0.8, -0.4, 0.5, 0.5, 0.5}, {-0.5, 0.8, 0.8, 0.4, 0.8, 0.8, 0.8}, {0, 0, 0, 0, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1}}, "Configurations for GFW regions"};

  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refP {2} refN {-2}", "piP {2} refN {-2}", "kaP {2} refN {-2}", "prP {2} refN {-2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}"}, {"ChGap22", "PiGap22", "KaGap22", "PrGap22", "ChFull22", "nchCh", "nchPi", "nchKa", "nchPr", "v02ptCh", "v02ptPi", "v02ptKa", "v02ptPr"}, {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {15, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, "Configurations for each correlation to calculate"};
  Configurable<LabeledArray<float>> nSigmas{"nSigmas", {LongArrayFloat.front().data(), 6, 3, {"UpCut_pi", "UpCut_ka", "UpCut_pr", "LowCut_pi", "LowCut_ka", "LowCut_pr"}, {"TPC", "TOF", "ITS"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};

  struct : ConfigurableGroup {
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.0f, "Maximum primary vertex cut applied for the events."};
    Configurable<int> cfgMultMin{"cfgMultMin", 10, "Minimum number of particles required for the event to have."};
  } cfgEventCuts;

  // // Filters to be applied to the received data.
  // // The analysis assumes the data has been subjected to a QA of its selection,
  // // and thus only the final distributions of the data for analysis are saved.
  o2::framework::expressions::Filter collFilter = (nabs(aod::collision::posZ) < cfgEventCuts.cfgZvtxMax);
  o2::framework::expressions::Filter trackFilter = (aod::track::pt > cfgTrackCuts.cfgPtMin) && (aod::track::pt < cfgTrackCuts.cfgPtMax) && nabs(aod::track::eta) < cfgTrackCuts.cfgEtaMax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true)) && (aod::track::itsChi2NCl < cfgTrackCuts.cfgChi2PrITSCls) && (aod::track::tpcChi2NCl < cfgTrackCuts.cfgChi2PrTPCCls) && nabs(aod::track::dcaZ) < cfgTrackCuts.cfgDCAz;

  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb{};

  struct Config {
    std::array<TH1D*, 4> mEfficiency{nullptr, nullptr, nullptr, nullptr};
    GFWWeights* mAcceptance = nullptr;
    bool correctionsLoaded = false;
  } cfg{};

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  HistogramRegistry registry{"registry"};

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

  std::unique_ptr<GFW> fGFW{std::make_unique<GFW>()};
  std::unique_ptr<TRandom3> fRndm{std::make_unique<TRandom3>(0)};
  std::unique_ptr<TAxis> fSecondAxis{nullptr};
  std::vector<GFW::CorrConfig> corrconfigs;
  int lastRun = -1;

  // region indices for consistency flag
  int posRegionIndex = -1;
  int negRegionIndex = -1;
  int fullRegionIndex = -1;
  int midRegionIndex = -1;
  // PID

  struct PIDState {
    o2::aod::ITSResponse itsResponse;
    std::array<float, 6> tofNsigmaCut{};
    std::array<float, 6> itsNsigmaCut{};
    std::array<float, 6> tpcNsigmaCut{};
    std::array<std::unique_ptr<TH1D>, 4> hPtMid{};
    std::array<std::unique_ptr<TH1D>, 4> hPtForward{};
    std::array<std::unique_ptr<TH1D>, 4> hPtBackward{};
  };
  PIDState pidStates;

  // Event selection cuts - Alex
  std::unique_ptr<TF1> fMultPVCutLow;
  std::unique_ptr<TF1> fMultPVCutHigh;
  std::unique_ptr<TF1> fMultCutLow;
  std::unique_ptr<TF1> fMultCutHigh;
  std::unique_ptr<TF1> fMultPVGlobalCutHigh;

  std::unique_ptr<TF1> fPtDepDCAxy;

  using GFWTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;

  enum PIDIndex {
    PidCharged = 0,
    PidPions,
    PidKaons,
    PidProtons,
    PidTotal
  };
  enum PiKpArrayIndex {
    IndPionUp = 0,
    IndKaonUp,
    IndProtonUp,
    IndPionLow,
    IndKaonLow,
    IndProtonLow
  };

  void init(InitContext const&)
  {

    projectMatrix(nSigmas->getData(), pidStates.tpcNsigmaCut, pidStates.tofNsigmaCut, pidStates.itsNsigmaCut);

    if (cfgGetNsigmaQA) {
      if (cfgUseItsPID) {
        registry.add("QA_PID/before/TofItsNsigma", "", {HistType::kTHnSparseD, {cfgPIDHistograms.cfgAxisNsigmaITS, cfgPIDHistograms.cfgAxisNsigmaTOF, cfgPIDHistograms.cfgAxisPt}});
      } else {
        registry.add("QA_PID/before/TofTpcNsigma_pions", "", {HistType::kTHnSparseD, {cfgPIDHistograms.cfgAxisNsigmaTPC, cfgPIDHistograms.cfgAxisNsigmaTOF, cfgPIDHistograms.cfgAxisPt}});
        registry.add("QA_PID/before/TofTpcNsigma_kaons", "", {HistType::kTHnSparseD, {cfgPIDHistograms.cfgAxisNsigmaTPC, cfgPIDHistograms.cfgAxisNsigmaTOF, cfgPIDHistograms.cfgAxisPt}});
        registry.add("QA_PID/before/TofTpcNsigma_protons", "", {HistType::kTHnSparseD, {cfgPIDHistograms.cfgAxisNsigmaTPC, cfgPIDHistograms.cfgAxisNsigmaTOF, cfgPIDHistograms.cfgAxisPt}});
      }

      registry.add("QA_PID/before/TpcdEdx_ptwise_pions", "", {HistType::kTHnSparseD, {{cfgPIDHistograms.cfgAxisPt, cfgPIDHistograms.cfgAxisTpcSignal, cfgPIDHistograms.cfgAxisNsigmaTOF}}});
      registry.add("QA_PID/before/ExpTpcdEdx_ptwise_pions", "", {HistType::kTHnSparseD, {{cfgPIDHistograms.cfgAxisPt, cfgPIDHistograms.cfgAxisTpcSignal, cfgPIDHistograms.cfgAxisNsigmaTOF}}});
      registry.add("QA_PID/before/ExpSigma_ptwise_pions", "", {HistType::kTHnSparseD, {{cfgPIDHistograms.cfgAxisPt, cfgPIDHistograms.cfgAxisSigma, cfgPIDHistograms.cfgAxisNsigmaTOF}}});
      registry.add("QA_PID/before/TpcdEdx_ptwise_kaons", "", {HistType::kTHnSparseD, {{cfgPIDHistograms.cfgAxisPt, cfgPIDHistograms.cfgAxisTpcSignal, cfgPIDHistograms.cfgAxisNsigmaTOF}}});
      registry.add("QA_PID/before/ExpTpcdEdx_ptwise_kaons", "", {HistType::kTHnSparseD, {{cfgPIDHistograms.cfgAxisPt, cfgPIDHistograms.cfgAxisTpcSignal, cfgPIDHistograms.cfgAxisNsigmaTOF}}});
      registry.add("QA_PID/before/ExpSigma_ptwise_kaons", "", {HistType::kTHnSparseD, {{cfgPIDHistograms.cfgAxisPt, cfgPIDHistograms.cfgAxisSigma, cfgPIDHistograms.cfgAxisNsigmaTOF}}});
      registry.add("QA_PID/before/TpcdEdx_ptwise_protons", "", {HistType::kTHnSparseD, {{cfgPIDHistograms.cfgAxisPt, cfgPIDHistograms.cfgAxisTpcSignal, cfgPIDHistograms.cfgAxisNsigmaTOF}}});
      registry.add("QA_PID/before/ExpTpcdEdx_ptwise_protons", "", {HistType::kTHnSparseD, {{cfgPIDHistograms.cfgAxisPt, cfgPIDHistograms.cfgAxisTpcSignal, cfgPIDHistograms.cfgAxisNsigmaTOF}}});
      registry.add("QA_PID/before/ExpSigma_ptwise_protons", "", {HistType::kTHnSparseD, {{cfgPIDHistograms.cfgAxisPt, cfgPIDHistograms.cfgAxisSigma, cfgPIDHistograms.cfgAxisNsigmaTOF}}});
      registry.addClone("QA_PID/before/", "QA_PID/after/");
    }

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
    gfwMemberCache.ptlow = cfgTrackCuts.cfgPtMin;
    gfwMemberCache.ptup = cfgTrackCuts.cfgPtMax;
    gfwMemberCache.etabins = cfgGFWBinning->GetEtaBins();
    gfwMemberCache.vtxZbins = cfgGFWBinning->GetVtxZbins();
    gfwMemberCache.phibins = cfgGFWBinning->GetPhiBins();
    gfwMemberCache.philow = 0.0f;
    gfwMemberCache.phiup = o2::constants::math::TwoPI;
    gfwMemberCache.nchbins = cfgGFWBinning->GetNchBins();
    gfwMemberCache.nchlow = cfgGFWBinning->GetNchMin();
    gfwMemberCache.nchup = cfgGFWBinning->GetNchMax();
    gfwMemberCache.centbinning = cfgGFWBinning->GetCentBinning();
    gfwMemberCache.nBootstrap = cfgNbootstrap;
    cfgGFWBinning->Print();
    gfwMemberCache.multGlobalCorrCutPars = cfgMultCorrCuts.cfgMultGlobalCutPars;
    gfwMemberCache.multPVCorrCutPars = cfgMultCorrCuts.cfgMultPVCutPars;
    gfwMemberCache.multGlobalPVCorrCutPars = cfgMultCorrCuts.cfgMultGlobalPVCutPars;

    // Initialise pt spectra histograms for different particles
    pidStates.hPtMid[PidCharged] = std::make_unique<TH1D>("hPtMid_charged", "hPtMid_charged", gfwMemberCache.ptbinning.size() - 1, gfwMemberCache.ptbinning.data());
    pidStates.hPtMid[PidPions] = std::make_unique<TH1D>("hPtMid_pions", "hPtMid_pions", gfwMemberCache.ptbinning.size() - 1, gfwMemberCache.ptbinning.data());
    pidStates.hPtMid[PidKaons] = std::make_unique<TH1D>("hPtMid_kaons", "hPtMid_kaons", gfwMemberCache.ptbinning.size() - 1, gfwMemberCache.ptbinning.data());
    pidStates.hPtMid[PidProtons] = std::make_unique<TH1D>("hPtMid_protons", "hPtMid_protons", gfwMemberCache.ptbinning.size() - 1, gfwMemberCache.ptbinning.data());
    pidStates.hPtMid[PidCharged]->SetDirectory(nullptr);
    pidStates.hPtMid[PidPions]->SetDirectory(nullptr);
    pidStates.hPtMid[PidKaons]->SetDirectory(nullptr);
    pidStates.hPtMid[PidProtons]->SetDirectory(nullptr);

    pidStates.hPtForward[PidCharged] = std::make_unique<TH1D>("hPtForward_charged", "hPtForward_charged", gfwMemberCache.ptbinning.size() - 1, gfwMemberCache.ptbinning.data());
    pidStates.hPtForward[PidPions] = std::make_unique<TH1D>("hPtForward_pions", "hPtForward_pions", gfwMemberCache.ptbinning.size() - 1, gfwMemberCache.ptbinning.data());
    pidStates.hPtForward[PidKaons] = std::make_unique<TH1D>("hPtForward_kaons", "hPtForward_kaons", gfwMemberCache.ptbinning.size() - 1, gfwMemberCache.ptbinning.data());
    pidStates.hPtForward[PidProtons] = std::make_unique<TH1D>("hPtForward_protons", "hPtForward_protons", gfwMemberCache.ptbinning.size() - 1, gfwMemberCache.ptbinning.data());
    pidStates.hPtForward[PidCharged]->SetDirectory(nullptr);
    pidStates.hPtForward[PidPions]->SetDirectory(nullptr);
    pidStates.hPtForward[PidKaons]->SetDirectory(nullptr);
    pidStates.hPtForward[PidProtons]->SetDirectory(nullptr);

    pidStates.hPtBackward[PidCharged] = std::make_unique<TH1D>("hPtBackward_charged", "hPtBackward_charged", gfwMemberCache.ptbinning.size() - 1, gfwMemberCache.ptbinning.data());
    pidStates.hPtBackward[PidPions] = std::make_unique<TH1D>("hPtBackward_pions", "hPtBackward_pions", gfwMemberCache.ptbinning.size() - 1, gfwMemberCache.ptbinning.data());
    pidStates.hPtBackward[PidKaons] = std::make_unique<TH1D>("hPtBackward_kaons", "hPtBackward_kaons", gfwMemberCache.ptbinning.size() - 1, gfwMemberCache.ptbinning.data());
    pidStates.hPtBackward[PidProtons] = std::make_unique<TH1D>("hPtBackward_protons", "hPtBackward_protons", gfwMemberCache.ptbinning.size() - 1, gfwMemberCache.ptbinning.data());
    pidStates.hPtBackward[PidCharged]->SetDirectory(nullptr);
    pidStates.hPtBackward[PidPions]->SetDirectory(nullptr);
    pidStates.hPtBackward[PidKaons]->SetDirectory(nullptr);
    pidStates.hPtBackward[PidProtons]->SetDirectory(nullptr);

    AxisSpec phiAxis = {gfwMemberCache.phibins, gfwMemberCache.philow, gfwMemberCache.phiup, "#phi"};
    AxisSpec etaAxis = {gfwMemberCache.etabins, -cfgTrackCuts.cfgEtaMax, cfgTrackCuts.cfgEtaMax, "#eta"};
    AxisSpec vtxAxis = {gfwMemberCache.vtxZbins, -cfgEventCuts.cfgZvtxMax, cfgEventCuts.cfgZvtxMax, "Vtx_{z} (cm)"};
    AxisSpec ptAxis = {gfwMemberCache.ptbinning, "#it{p}_{T} GeV/#it{c}"};

    std::string sCentralityEstimator = "FT0C centrality (%)";
    AxisSpec centAxis = {gfwMemberCache.centbinning, sCentralityEstimator.c_str()};

    std::vector<double> nchbinning;
    const double nchskip = (gfwMemberCache.nchup - gfwMemberCache.nchlow) / static_cast<double>(gfwMemberCache.nchbins);
    for (int i = 0; i <= gfwMemberCache.nchbins; ++i) {
      nchbinning.push_back(nchskip * i + gfwMemberCache.nchlow + 0.5);
    }
    AxisSpec nchAxis = {nchbinning, "N_{ch}"};
    AxisSpec t0cAxis = {1000, 0, 10000, "N_{ch} (T0C)"};
    AxisSpec t0aAxis = {300, 0, 30000, "N_{ch} (T0A)"};
    AxisSpec v0aAxis = {800, 0, 80000, "N_{ch} (V0A)"};
    AxisSpec multpvAxis = {600, 0, 600, "N_{ch} (PV)"};
    AxisSpec dcaZAxis = {200, -2, 2, "DCA_{z} (cm)"};
    AxisSpec dcaXYAxis = {200, -0.5, 0.5, "DCA_{xy} (cm)"};
    AxisSpec bsAxis = {gfwMemberCache.nBootstrap, -0.5, gfwMemberCache.nBootstrap - 0.5, "Bootstrap Index"};

    registry.add("v02pt", "", {HistType::kTProfile3D, {ptAxis, centAxis, nchAxis}});
    registry.add("nchMid", "", {HistType::kTProfile3D, {ptAxis, centAxis, nchAxis}});
    registry.add("v02centmult", "", {HistType::kTProfile2D, {centAxis, nchAxis}});

    registry.add("analysis/charged/v0AB", "", {HistType::kTProfile3D, {bsAxis, ptAxis, centAxis}});
    registry.add("analysis/charged/v0BA", "", {HistType::kTProfile3D, {bsAxis, ptAxis, centAxis}});
    registry.add("analysis/charged/nchA", "", {HistType::kTProfile3D, {bsAxis, ptAxis, centAxis}});
    registry.add("analysis/charged/nchB", "", {HistType::kTProfile3D, {bsAxis, ptAxis, centAxis}});
    registry.add("analysis/charged/ptA", "", {HistType::kTProfile3D, {bsAxis, centAxis, nchAxis}});
    registry.add("analysis/charged/ptB", "", {HistType::kTProfile3D, {bsAxis, centAxis, nchAxis}});
    registry.add("analysis/charged/ptAB", "", {HistType::kTProfile3D, {bsAxis, centAxis, nchAxis}});

    registry.addClone("analysis/charged/", "analysis/pion/");
    registry.addClone("analysis/charged/", "analysis/kaon/");
    registry.addClone("analysis/charged/", "analysis/proton/");

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    int ptbins = gfwMemberCache.ptbinning.size() - 1;
    fSecondAxis = std::make_unique<TAxis>(ptbins, gfwMemberCache.ptbinning.data());

    // QA histograms
    registry.add("trackQA/before/phi_eta_vtxZ", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
    registry.add("trackQA/before/pt_dcaXY_dcaZ", "", {HistType::kTH3D, {ptAxis, dcaXYAxis, dcaZAxis}});
    registry.add("trackQA/before/nch_pt", "#it{p}_{T} vs multiplicity; N_{ch}; #it{p}_{T}", {HistType::kTH2D, {nchAxis, ptAxis}});
    registry.add("trackQA/after/pt_ref", "", {HistType::kTH1D, {{100, gfwMemberCache.ptreflow, gfwMemberCache.ptrefup}}});
    registry.add("trackQA/after/pt_poi", "", {HistType::kTH1D, {{100, gfwMemberCache.ptpoilow, gfwMemberCache.ptpoiup}}});
    registry.add("trackQA/before/chi2prTPCcls", "#chi^{2}/cluster for the TPC track segment; #chi^{2}/TPC cluster", {HistType::kTH1D, {{100, 0., 5.}}});
    registry.add("trackQA/before/chi2prITScls", "#chi^{2}/cluster for the ITS track; #chi^{2}/ITS cluster", {HistType::kTH1D, {{100, 0., 50.}}});
    registry.add("trackQA/before/nTPCClusters", "Number of found TPC clusters; TPC N_{cls}; Counts", {HistType::kTH1D, {{100, 40, 180}}});
    registry.add("trackQA/before/nITSClusters", "Number of found ITS clusters; ITS N_{cls}; Counts", {HistType::kTH1D, {{100, 0, 20}}});
    registry.add("trackQA/before/nTPCCrossedRows", "Number of crossed TPC Rows; TPC X-rows; Counts", {HistType::kTH1D, {{100, 40, 180}}});
    registry.addClone("trackQA/before/", "trackQA/after/");

    registry.add("eventQA/before/globalTracks_centT0C", "", {HistType::kTH2D, {centAxis, nchAxis}});
    registry.add("eventQA/before/PVTracks_centT0C", "", {HistType::kTH2D, {centAxis, multpvAxis}});
    registry.add("eventQA/before/multT0C_centT0C", "", {HistType::kTH2D, {centAxis, t0cAxis}});

    registry.add("eventQA/before/centT0M_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
    registry.add("eventQA/before/centV0A_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
    registry.add("eventQA/before/centGlobal_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
    registry.add("eventQA/before/centNTPV_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});
    registry.add("eventQA/before/centMFT_centT0C", "", {HistType::kTH2D, {centAxis, centAxis}});

    registry.add("eventQA/before/globalTracks_PVTracks", "", {HistType::kTH2D, {multpvAxis, nchAxis}});
    registry.add("eventQA/before/globalTracks_multT0A", "", {HistType::kTH2D, {t0aAxis, nchAxis}});
    registry.add("eventQA/before/globalTracks_multV0A", "", {HistType::kTH2D, {v0aAxis, nchAxis}});
    registry.add("eventQA/before/multV0A_multT0A", "", {HistType::kTH2D, {t0aAxis, v0aAxis}});

    registry.add("eventQA/before/multiplicity", "", {HistType::kTH1D, {nchAxis}});
    registry.add("eventQA/before/centrality", "", {HistType::kTH1D, {centAxis}});
    registry.addClone("eventQA/before/", "eventQA/after/");

    // Event selection histograms
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

    if (gfwMemberCache.regions.GetSize() < 0)
      LOGF(error, "Configuration contains vectors of different size - check the GFWRegions configurable");
    for (auto i(0); i < gfwMemberCache.regions.GetSize(); ++i) {
      fGFW->AddRegion(gfwMemberCache.regions.GetNames()[i], gfwMemberCache.regions.GetEtaMin()[i], gfwMemberCache.regions.GetEtaMax()[i], (gfwMemberCache.regions.GetpTDifs()[i] != 0) ? ptbins + 1 : 1, gfwMemberCache.regions.GetBitmasks()[i]);
    }
    for (auto i = 0; i < gfwMemberCache.configs.GetSize(); ++i) {
      corrconfigs.push_back(fGFW->GetCorrelatorConfig(gfwMemberCache.configs.GetCorrs()[i], gfwMemberCache.configs.GetHeads()[i], gfwMemberCache.configs.GetpTDifs()[i] != 0));
    }
    if (corrconfigs.empty())
      LOGF(error, "Configuration contains vectors of different size - check the GFWCorrConfig configurable");
    fGFW->CreateRegions();
    auto* oba = new TObjArray();
    oba->SetOwner(kTRUE);
    addConfigObjectsToObjArray(oba, corrconfigs);
    LOGF(info, "Number of correlators: %d", oba->GetEntries());
    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fSecondAxis.get());
    fFC->Initialize(oba, centAxis, cfgNbootstrap);
    delete oba;

    if (cfgConsistentEventFlag != 0) {
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

    if (cfgUseAdditionalEventCut) {
      fMultPVCutLow = std::make_unique<TF1>("fMultPVCutLow", cfgMultCorrCuts.cfgMultCorrLowCutFunction->c_str(), 0, 100);
      fMultPVCutLow->SetParameters(gfwMemberCache.multPVCorrCutPars.data());
      fMultPVCutHigh = std::make_unique<TF1>("fMultPVCutHigh", cfgMultCorrCuts.cfgMultCorrHighCutFunction->c_str(), 0, 100);
      fMultPVCutHigh->SetParameters(gfwMemberCache.multPVCorrCutPars.data());
      fMultCutLow = std::make_unique<TF1>("fMultCutLow", cfgMultCorrCuts.cfgMultCorrLowCutFunction->c_str(), 0, 100);
      fMultCutLow->SetParameters(gfwMemberCache.multGlobalCorrCutPars.data());
      fMultCutHigh = std::make_unique<TF1>("fMultCutHigh", cfgMultCorrCuts.cfgMultCorrHighCutFunction->c_str(), 0, 100);
      fMultCutHigh->SetParameters(gfwMemberCache.multGlobalCorrCutPars.data());
      fMultPVGlobalCutHigh = std::make_unique<TF1>("fMultPVGlobalCutHigh", cfgMultCorrCuts.cfgMultGlobalPVCorrCutFunction->c_str(), 0, nchbinning.back());
      fMultPVGlobalCutHigh->SetParameters(gfwMemberCache.multGlobalPVCorrCutPars.data());
    }
    // Set DCAxy cut
    fPtDepDCAxy = std::make_unique<TF1>("ptDepDCAxy", Form("[0]*%s", cfgTrackCuts.cfgDCAxy->c_str()), 0.001, 100);
    fPtDepDCAxy->SetParameter(0, cfgTrackCuts.cfgDCAxyNSigma);
  }

  static constexpr std::array<std::string_view, 2> FillTimeName = {"before/", "after/"};
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
          oba->Add(new TNamed(it->Head + index, it->Head + suffix));
        }
      } else {
        oba->Add(new TNamed(it->Head.c_str(), it->Head.c_str()));
      }
    }
  }

  template <typename TTrack>
  int getNsigmaPID(TTrack const& track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {pidStates.itsResponse.nSigmaITS<o2::track::PID::Pion>(track), pidStates.itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), pidStates.itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = -1; // -1 = not identified, 1 = pion, 2 = kaon, 3 = proton

    std::array<float, 3> nSigmaToUse = cfgUseItsPID ? nSigmaITS : nSigmaTPC;                                 // Choose which nSigma to use: TPC or ITS
    std::array<float, 6> detectorNsigmaCut = cfgUseItsPID ? pidStates.itsNsigmaCut : pidStates.tpcNsigmaCut; // Choose which nSigma to use: TPC or ITS

    bool isPion = false, isKaon = false, isProton = false;
    bool isDetectedPion = nSigmaToUse[IndPionUp] < detectorNsigmaCut[IndPionUp] && nSigmaToUse[IndPionUp] > detectorNsigmaCut[IndPionLow];
    bool isDetectedKaon = nSigmaToUse[IndKaonUp] < detectorNsigmaCut[IndKaonUp] && nSigmaToUse[IndKaonUp] > detectorNsigmaCut[IndKaonLow];
    bool isDetectedProton = nSigmaToUse[IndProtonUp] < detectorNsigmaCut[IndProtonUp] && nSigmaToUse[IndProtonUp] > detectorNsigmaCut[IndProtonLow];

    bool isTofPion = nSigmaTOF[IndPionUp] < pidStates.tofNsigmaCut[IndPionUp] && nSigmaTOF[IndPionUp] > pidStates.tofNsigmaCut[IndPionLow];
    bool isTofKaon = nSigmaTOF[IndKaonUp] < pidStates.tofNsigmaCut[IndKaonUp] && nSigmaTOF[IndKaonUp] > pidStates.tofNsigmaCut[IndKaonLow];
    bool isTofProton = nSigmaTOF[IndProtonUp] < pidStates.tofNsigmaCut[IndProtonUp] && nSigmaTOF[IndProtonUp] > pidStates.tofNsigmaCut[IndProtonLow];

    if (track.pt() > cfgTofPtCut && !track.hasTOF()) {
      return -1;
    }
    if (track.pt() > cfgTofPtCut && track.hasTOF()) {
      isPion = isTofPion && isDetectedPion;
      isKaon = isTofKaon && isDetectedKaon;
      isProton = isTofProton && isDetectedProton;
    } else {
      isPion = isDetectedPion;
      isKaon = isDetectedKaon;
      isProton = isDetectedProton;
    }

    if ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton)) {
      return -1; // more than one particle satisfy the criteria
    }

    if (isPion) {
      pid = PidPions;
    } else if (isKaon) {
      pid = PidKaons;
    } else if (isProton) {
      pid = PidProtons;
    } else {
      return -1; // no particle satisfies the criteria
    }

    return pid; // -1 = not identified, 1 = pion, 2 = kaon, 3 = proton
  }

  void loadCorrections(aod::BCsWithTimestamps::iterator const& bc)
  {
    uint64_t timestamp = bc.timestamp();
    if (cfg.correctionsLoaded)
      return;
    if (!cfgAcceptance.value.empty()) {
      cfg.mAcceptance = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value, timestamp);
    }
    if (!cfgEfficiency.value.empty()) {
      cfg.mEfficiency[PidCharged] = ccdb->getForTimeStamp<TH1D>(cfgEfficiency, timestamp);
      if (cfg.mEfficiency[PidCharged] == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), static_cast<void*>(cfg.mEfficiency[PidCharged]));
    }
    if (cfgPIDEfficiency) {
      const std::array<std::string, 4> pidStrings = {"ch", "pi", "ka", "pr"};
      for (int i = 1; i < PidTotal; i++) {

        cfg.mEfficiency[i] = ccdb->getForTimeStamp<TH1D>(cfgEfficiency.value + pidStrings[i], timestamp);
        if (cfg.mEfficiency[i] == nullptr) {
          LOGF(fatal, "Could not load PID efficiency histogram from %s", (cfgEfficiency.value + pidStrings[i]).c_str());
        }
        LOGF(info, "Loaded PID efficiency histogram from %s (%p)", (cfgEfficiency.value + pidStrings[i]).c_str(), static_cast<void*>(cfg.mEfficiency[i]));
      }
    }
    cfg.correctionsLoaded = true;
  }

  void loadCorrections(int runnumber)
  {
    if (cfg.correctionsLoaded)
      return;
    if (!cfgAcceptance.value.empty()) {
      cfg.mAcceptance = ccdb->getForRun<GFWWeights>(cfgAcceptance.value, runnumber);
    }
    if (!cfgEfficiency.value.empty()) {
      cfg.mEfficiency[PidCharged] = ccdb->getForRun<TH1D>(cfgEfficiency.value, runnumber);
    }
    cfg.correctionsLoaded = true;
  }

  template <typename TTrack>
  double getJTrackAcceptance(TTrack const& track)
  {
    double wacc = 1;
    if constexpr (requires { track.weightNUA(); })
      wacc = 1. / track.weightNUA();
    return wacc;
  }

  template <typename TTrack>
  double getJTrackEfficiency(TTrack const& track)
  {
    double eff = 1.;
    if constexpr (requires { track.weightEff(); })
      eff = track.weightEff();
    return eff;
  }

  template <typename TTrack>
  double getAcceptance(TTrack const& track, const double& vtxz)
  {
    double wacc = 1;
    if (cfg.mAcceptance)
      wacc = cfg.mAcceptance->getNUA(track.phi(), track.eta(), vtxz);
    return wacc;
  }

  template <typename TTrack>
  double getEfficiency(TTrack const& track, const int& pid = PidCharged)
  {
    double eff = 1.;
    if (cfg.mEfficiency[pid])
      eff = cfg.mEfficiency[pid]->GetBinContent(cfg.mEfficiency[pid]->FindBin(track.pt()));
    if (eff == 0)
      return -1.;
    return 1. / eff;
  }

  template <typename TCollision>
  bool eventSelected(TCollision const& collision, const int& multTrk, const float& centrality)
  {
    if (cfgEventCutFlags.cfgTVXinTRD) {
      if (collision.alias_bit(kTVXinTRD)) {
        // TRD triggered
        // "CMTVX-B-NOPF-TRD,minbias_TVX"
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kTVXinTRD);
    }
    if (cfgEventCutFlags.cfgNoSameBunchPileupCut) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kNoSameBunchPileup);
    }
    if (cfgEventCutFlags.cfgIsGoodZvtxFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kIsGoodZvtxFT0vsPV);
    }
    if (cfgEventCutFlags.cfgNoCollInTimeRangeStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kNoCollInTimeRangeStandard);
    }

    if (cfgEventCutFlags.cfgNoCollInRofStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kNoCollInRofStandard);
    }

    if (cfgEventCutFlags.cfgNoHighMultCollInPrevRof) {
      if (!collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kNoHighMultCollInPrevRof);
    }

    if (cfgEventCutFlags.cfgIsVertexITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        // selects collisions with at least one ITS-TPC track, and thus rejects vertices built from ITS-only tracks
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kIsVertexITSTPC);
    }

    if (cfgEventCutFlags.cfgIsGoodITSLayersAll) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kIsGoodITSLayersAll);
    }

    if (cfgEventCutFlags.cfgNoTimeFrameBorder) {
      if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kNoTimeFrameBorder);
    }

    if (cfgEventCutFlags.cfgNoITSROFrameBorder) {
      if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), kNoITSROFrameBorder);
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

    if (vtxz > gfwMemberCache.vtxZup || vtxz < gfwMemberCache.vtxZlow)
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
      registry.fill(HIST("eventQA/eventSel"), kMultCuts);
    }
    return 1;
  }

  // Define the data type
  enum DataType {
    kReco,
    kGen
  };

  int getPIDIndex(const std::string& corrconfig)
  {
    if (!boost::ifind_first(corrconfig, "pi").empty())
      return PidPions;
    if (!boost::ifind_first(corrconfig, "ka").empty())
      return PidKaons;
    if (!boost::ifind_first(corrconfig, "pr").empty())
      return PidProtons;
    return PidCharged;
  }

  template <DataType dt>
  void fillOutputContainers(const float& centmult, const int& multiplicity, const double& rndm, const int& /*run*/ = 0)
  {
    for (uint l_ind = 0; l_ind < corrconfigs.size(); ++l_ind) {
      if (!corrconfigs.at(l_ind).pTDif) {
        auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), 0, kTRUE).real();
        if (dnx == 0)
          continue;
        auto val = fGFW->Calculate(corrconfigs.at(l_ind), 0, kFALSE).real() / dnx;

        if (std::abs(val) < 1) {
          fFC->FillProfile(corrconfigs.at(l_ind).Head.c_str(), centmult, val, (cfgUseMultiplicityFlowWeights) ? dnx : 1.0, rndm);
        }
        continue;
      }

      // Fill pt profiles for different particles
      int pidInd = getPIDIndex(corrconfigs.at(l_ind).Head);

      auto dnx = fGFW->Calculate(corrconfigs.at(0), 0, kTRUE).real();
      if (dnx == 0)
        continue;
      auto val = fGFW->Calculate(corrconfigs.at(0), 0, kFALSE).real() / dnx;
      double ebyeWeight = (cfgUseMultiplicityFlowWeights) ? dnx : 1.0;
      for (int i = 1; i <= fSecondAxis->GetNbins(); i++) {
        if (corrconfigs.at(l_ind).Head.find("nch") != std::string::npos) {
          ebyeWeight = 1.0;
          val = 1.0;
        }
        if (cfgUseMultiplicityFracWeights && pidStates.hPtMid[PidCharged]->Integral() > 0) {
          ebyeWeight *= pidStates.hPtMid[PidCharged]->Integral();
        }
        double ptFraction = 0;
        double threshold = 1.01;
        int normIndex = (cfgNormalizeByCharged) ? PidCharged : pidInd; // Configured to normalize by charged particles or the selected particle
        if (pidStates.hPtMid[normIndex]->Integral() > 0) {
          ptFraction = pidStates.hPtMid[pidInd]->GetBinContent(i) / pidStates.hPtMid[normIndex]->Integral();
          if (std::abs(val) < threshold)
            fFC->FillProfile(Form("%s_pt_%i", corrconfigs.at(l_ind).Head.c_str(), i), centmult, val * ptFraction, ebyeWeight, rndm);
        }
      }
    }

    if (cfgUseV0) {
      double v0corrAB = 0;
      double v0corrBA = 0;
      double ptMeanForward = pidStates.hPtForward[PidCharged]->GetMean();
      double ptMeanBackward = pidStates.hPtBackward[PidCharged]->GetMean();
      double ptFractionForward = 0.;
      double ptFractionBackward = 0.;
      int bootstrap = fRndm->Integer(gfwMemberCache.nBootstrap);
      for (int pid = 0; pid < PidTotal; pid++) {
        int normIndex = (cfgNormalizeByCharged) ? PidCharged : pid;
        if (!(pidStates.hPtForward[normIndex]->Integral() > 0) || !(pidStates.hPtBackward[normIndex]->Integral() > 0)) {
          continue; // Forward or backward pT distribution is not defined
        }
        for (int i = 1; i <= fSecondAxis->GetNbins(); i++) {
          ptFractionForward = pidStates.hPtForward[pid]->GetBinContent(i) / pidStates.hPtForward[normIndex]->Integral();
          ptFractionBackward = pidStates.hPtBackward[pid]->GetBinContent(i) / pidStates.hPtBackward[normIndex]->Integral();
          v0corrAB = ptFractionForward * ptMeanBackward;
          v0corrBA = ptFractionBackward * ptMeanForward;
          switch (pid) {
            case PidCharged:
              registry.fill(HIST("analysis/charged/v0AB"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, v0corrAB);
              registry.fill(HIST("analysis/charged/v0BA"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, v0corrBA);
              registry.fill(HIST("analysis/charged/nchA"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, ptFractionForward);
              registry.fill(HIST("analysis/charged/nchB"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, ptFractionBackward);
              break;
            case PidPions:
              registry.fill(HIST("analysis/pion/v0AB"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, v0corrAB);
              registry.fill(HIST("analysis/pion/v0BA"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, v0corrBA);
              registry.fill(HIST("analysis/pion/nchA"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, ptFractionForward);
              registry.fill(HIST("analysis/pion/nchB"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, ptFractionBackward);
              break;
            case PidKaons:
              registry.fill(HIST("analysis/kaon/v0AB"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, v0corrAB);
              registry.fill(HIST("analysis/kaon/v0BA"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, v0corrBA);
              registry.fill(HIST("analysis/kaon/nchA"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, ptFractionForward);
              registry.fill(HIST("analysis/kaon/nchB"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, ptFractionBackward);
              break;
            case PidProtons:
              registry.fill(HIST("analysis/proton/v0AB"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, v0corrAB);
              registry.fill(HIST("analysis/proton/v0BA"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, v0corrBA);
              registry.fill(HIST("analysis/proton/nchA"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, ptFractionForward);
              registry.fill(HIST("analysis/proton/nchB"), bootstrap, fSecondAxis->GetBinCenter(i), centmult, ptFractionBackward);
              break;
            default:
              break;
          }
        }
        switch (pid) {
          case PidCharged:
            registry.fill(HIST("analysis/charged/ptA"), bootstrap, centmult, multiplicity, ptMeanForward);
            registry.fill(HIST("analysis/charged/ptB"), bootstrap, centmult, multiplicity, ptMeanBackward);
            registry.fill(HIST("analysis/charged/ptAB"), bootstrap, centmult, multiplicity, ptMeanForward * ptMeanBackward);
            break;
          case PidPions:
            registry.fill(HIST("analysis/pion/ptA"), bootstrap, centmult, multiplicity, ptMeanForward);
            registry.fill(HIST("analysis/pion/ptB"), bootstrap, centmult, multiplicity, ptMeanBackward);
            registry.fill(HIST("analysis/pion/ptAB"), bootstrap, centmult, multiplicity, ptMeanForward * ptMeanBackward);
            break;
          case PidKaons:
            registry.fill(HIST("analysis/kaon/ptA"), bootstrap, centmult, multiplicity, ptMeanForward);
            registry.fill(HIST("analysis/kaon/ptB"), bootstrap, centmult, multiplicity, ptMeanBackward);
            registry.fill(HIST("analysis/kaon/ptAB"), bootstrap, centmult, multiplicity, ptMeanForward * ptMeanBackward);
            break;
          case PidProtons:
            registry.fill(HIST("analysis/proton/ptA"), bootstrap, centmult, multiplicity, ptMeanForward);
            registry.fill(HIST("analysis/proton/ptB"), bootstrap, centmult, multiplicity, ptMeanBackward);
            registry.fill(HIST("analysis/proton/ptAB"), bootstrap, centmult, multiplicity, ptMeanForward * ptMeanBackward);
            break;
          default:
            break;
        }
      }
    }

    // Fill the profiles for each pT bin
    auto dnx = fGFW->Calculate(corrconfigs.at(0), 0, kTRUE).real();
    if (dnx == 0)
      return;
    auto val = fGFW->Calculate(corrconfigs.at(0), 0, kFALSE).real() / dnx;
    for (int i = 1; i <= fSecondAxis->GetNbins(); i++) {
      double ptFraction = 0;
      if (pidStates.hPtMid[PidCharged]->Integral() > 0) {
        ptFraction = pidStates.hPtMid[PidCharged]->GetBinContent(i) / pidStates.hPtMid[PidCharged]->Integral();
        double threshold = 1.01;
        if (std::abs(val) < threshold)
          registry.fill(HIST("v02pt"), fSecondAxis->GetBinCenter(i), centmult, multiplicity, val * ptFraction, (cfgUseMultiplicityFlowWeights) ? dnx : 1.0);
        registry.fill(HIST("nchMid"), fSecondAxis->GetBinCenter(i), centmult, multiplicity, ptFraction);
      }
    }
    registry.fill(HIST("v02centmult"), centmult, multiplicity, val);
  }

  struct XAxis {
    float centrality;
    int64_t multiplicity;
  };

  struct AcceptedTracks {
    int nPos;
    int nNeg;
    int nFull;
    int nMid;
  };

  template <DataType dt, typename TCollision, typename TTracks>
  void processCollision(TCollision const& collision, TTracks const& tracks, const XAxis& xaxis, const int& run)
  {
    float vtxz = collision.posZ();
    if (tracks.size() < 1)
      return;
    if (xaxis.centrality >= 0 && (xaxis.centrality < gfwMemberCache.centbinning.front() || xaxis.centrality > gfwMemberCache.centbinning.back()))
      return;
    if (xaxis.multiplicity < cfgFixedMultMin || xaxis.multiplicity > cfgFixedMultMax)
      return;
    fGFW->Clear();
    pidStates.hPtMid[PidCharged]->Reset();
    pidStates.hPtMid[PidPions]->Reset();
    pidStates.hPtMid[PidKaons]->Reset();
    pidStates.hPtMid[PidProtons]->Reset();
    pidStates.hPtBackward[PidCharged]->Reset();
    pidStates.hPtBackward[PidPions]->Reset();
    pidStates.hPtBackward[PidKaons]->Reset();
    pidStates.hPtBackward[PidProtons]->Reset();
    pidStates.hPtForward[PidCharged]->Reset();
    pidStates.hPtForward[PidPions]->Reset();
    pidStates.hPtForward[PidKaons]->Reset();
    pidStates.hPtForward[PidProtons]->Reset();

    float lRandom = fRndm->Rndm();

    // Loop over tracks and check if they are accepted
    AcceptedTracks acceptedTracks{.nPos = 0, .nNeg = 0, .nFull = 0, .nMid = 0};
    for (const auto& track : tracks) {
      processTrack(track, vtxz, xaxis.multiplicity, run, acceptedTracks);
      if (track.eta() > cfgSubeventCuts.cfgEtaSubCMin && track.eta() < cfgSubeventCuts.cfgEtaSubCMax)
        pidStates.hPtMid[PidCharged]->Fill(track.pt(), getEfficiency(track, PidCharged));
      if (track.eta() > cfgSubeventCuts.cfgEtaSubAMin && track.eta() < cfgSubeventCuts.cfgEtaSubAMax) // add mean pT
        pidStates.hPtBackward[PidCharged]->Fill(track.pt(), getEfficiency(track, PidCharged));
      if (track.eta() > cfgSubeventCuts.cfgEtaSubBMin && track.eta() < cfgSubeventCuts.cfgEtaSubBMax) // add mean pT
        pidStates.hPtForward[PidCharged]->Fill(track.pt(), getEfficiency(track, PidCharged));
      // If PID is identified, fill pt spectrum for the corresponding particle
      int pidInd = getNsigmaPID(track);
      if (pidInd != -1 && track.eta() > cfgSubeventCuts.cfgEtaSubCMin && track.eta() < cfgSubeventCuts.cfgEtaSubCMax) {
        if (cfgPIDEfficiency)
          pidStates.hPtMid[pidInd]->Fill(track.pt(), getEfficiency(track, pidInd));
        else
          pidStates.hPtMid[pidInd]->Fill(track.pt(), getEfficiency(track, PidCharged)); // Default to charged particles if PID efficiency is not used
      }
      if (pidInd != -1 && track.eta() > cfgSubeventCuts.cfgEtaSubAMin && track.eta() < cfgSubeventCuts.cfgEtaSubAMax) {
        if (cfgPIDEfficiency)
          pidStates.hPtBackward[pidInd]->Fill(track.pt(), getEfficiency(track, pidInd));
        else
          pidStates.hPtBackward[pidInd]->Fill(track.pt(), getEfficiency(track, PidCharged)); // Default to charged particles if PID efficiency is not used
      }
      if (pidInd != -1 && track.eta() > cfgSubeventCuts.cfgEtaSubBMin && track.eta() < cfgSubeventCuts.cfgEtaSubBMax) {
        if (cfgPIDEfficiency)
          pidStates.hPtForward[pidInd]->Fill(track.pt(), getEfficiency(track, pidInd));
        else
          pidStates.hPtForward[pidInd]->Fill(track.pt(), getEfficiency(track, PidCharged)); // Default to charged particles if PID efficiency is not used
      }
    }
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
    // Fill output containers
    fillOutputContainers<dt>(xaxis.centrality, xaxis.multiplicity, lRandom, run);
  }

  template <typename TTrack>
  void fillAcceptedTracks(TTrack const& track, AcceptedTracks& acceptedTracks)
  {
    if (posRegionIndex >= 0 && track.eta() > gfwMemberCache.regions.GetEtaMin()[posRegionIndex] && track.eta() < gfwMemberCache.regions.GetEtaMax()[posRegionIndex])
      ++acceptedTracks.nPos;
    if (negRegionIndex >= 0 && track.eta() > gfwMemberCache.regions.GetEtaMin()[negRegionIndex] && track.eta() < gfwMemberCache.regions.GetEtaMax()[negRegionIndex])
      ++acceptedTracks.nNeg;
    if (fullRegionIndex >= 0 && track.eta() > gfwMemberCache.regions.GetEtaMin()[fullRegionIndex] && track.eta() < gfwMemberCache.regions.GetEtaMax()[fullRegionIndex])
      ++acceptedTracks.nFull;
    if (midRegionIndex >= 0 && track.eta() > gfwMemberCache.regions.GetEtaMin()[midRegionIndex] && track.eta() < gfwMemberCache.regions.GetEtaMax()[midRegionIndex])
      ++acceptedTracks.nMid;
  }

  template <typename TTrack>
  bool trackSelected(TTrack const& track)
  {
    if (cfgTrackCuts.cfgDCAxyNSigma && (std::fabs(track.dcaXY()) > fPtDepDCAxy->Eval(track.pt())))
      return false;
    return ((track.tpcNClsCrossedRows() >= cfgTrackCuts.cfgNTPCXrows) && (track.tpcNClsFound() >= cfgTrackCuts.cfgNTPCCls) && (track.itsNCls() >= cfgTrackCuts.cfgMinNITSCls));
  }

  template <typename TCollision>
  float getCentrality(TCollision const& collision)
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

  template <typename TTrack>
  inline void processTrack(TTrack const& track, const float& vtxz, const int& multiplicity, const int& /*run*/, AcceptedTracks& acceptedTracks)
  {

    if (cfgFillQA) {
      fillTrackQA<kBefore>(track, vtxz);
      registry.fill(HIST("trackQA/before/nch_pt"), multiplicity, track.pt());
    }

    if (cfgGetNsigmaQA)
      fillPidQA<kBefore>(track, getNsigmaPID(track));

    if (!trackSelected(track))
      return;

    fillGFW<kReco>(track, vtxz);               // Fill GFW
    fillAcceptedTracks(track, acceptedTracks); // Fill accepted tracks
    if (cfgFillQA) {
      fillTrackQA<kAfter>(track, vtxz);
      registry.fill(HIST("trackQA/after/nch_pt"), multiplicity, track.pt());
    }

    if (cfgGetNsigmaQA)
      fillPidQA<kAfter>(track, getNsigmaPID(track));
  }

  template <DataType dt, typename TTrack>
  inline void fillGFW(TTrack const& track, const double& vtxz)
  {
    int pidInd = getNsigmaPID(track);

    bool withinPtRef = (track.pt() > gfwMemberCache.ptreflow && track.pt() < gfwMemberCache.ptrefup);
    bool withinPtPOI = (track.pt() > gfwMemberCache.ptpoilow && track.pt() < gfwMemberCache.ptpoiup);

    if (!withinPtPOI && !withinPtRef)
      return;
    double weff = getEfficiency(track, PidCharged);
    if (weff < 0)
      return;

    double wacc = getAcceptance(track, vtxz);

    // Fill cumulants for different particles
    // ***Need to add proper weights for each particle!***
    if (withinPtRef)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 1);
    if (withinPtPOI && pidInd == PidPions)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, PidPions + 1);
    if (withinPtPOI && pidInd == PidKaons)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, PidKaons + 1);
    if (withinPtPOI && pidInd == PidProtons)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, PidProtons + 1);
  }

  template <QAFillTime ft, typename TTrack>
  inline void fillPidQA(TTrack const& track, const int& pid)
  {
    // Fill Nsigma QA
    if (!cfgUseItsPID) {
      if (ft == kBefore) {
        registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("TofTpcNsigma_pions"), track.tpcNSigmaPi(), track.tofNSigmaPi(), track.pt());
        if (cfgGetdEdx) {
          double tpcExpSignalPi = track.tpcSignal() - (track.tpcNSigmaPi() * track.tpcExpSigmaPi());

          registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("TpcdEdx_ptwise_pions"), track.pt(), track.tpcSignal(), track.tofNSigmaPi());
          registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("ExpTpcdEdx_ptwise_pions"), track.pt(), tpcExpSignalPi, track.tofNSigmaPi());
          registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("ExpSigma_ptwise_pions"), track.pt(), track.tpcExpSigmaPi(), track.tofNSigmaPi());
        }
        registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("TofTpcNsigma_kaons"), track.tpcNSigmaKa(), track.tofNSigmaKa(), track.pt());
        if (cfgGetdEdx) {
          double tpcExpSignalKa = track.tpcSignal() - (track.tpcNSigmaKa() * track.tpcExpSigmaKa());

          registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("TpcdEdx_ptwise_kaons"), track.pt(), track.tpcSignal(), track.tofNSigmaKa());
          registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("ExpTpcdEdx_ptwise_kaons"), track.pt(), tpcExpSignalKa, track.tofNSigmaKa());
          registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("ExpSigma_ptwise_kaons"), track.pt(), track.tpcExpSigmaKa(), track.tofNSigmaKa());
        }
        registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("TofTpcNsigma_protons"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
        if (cfgGetdEdx) {
          double tpcExpSignalPr = track.tpcSignal() - (track.tpcNSigmaPr() * track.tpcExpSigmaPr());

          registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("TpcdEdx_ptwise_protons"), track.pt(), track.tpcSignal(), track.tofNSigmaPr());
          registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("ExpTpcdEdx_ptwise_protons"), track.pt(), tpcExpSignalPr, track.tofNSigmaPr());
          registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("ExpSigma_ptwise_protons"), track.pt(), track.tpcExpSigmaPr(), track.tofNSigmaPr());
        }
      } else if (ft == kAfter) {
        if (pid == PidPions) {
          registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("TofTpcNsigma_pions"), track.tpcNSigmaPi(), track.tofNSigmaPi(), track.pt());
          if (cfgGetdEdx) {
            double tpcExpSignalPi = track.tpcSignal() - (track.tpcNSigmaPi() * track.tpcExpSigmaPi());

            registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("TpcdEdx_ptwise_pions"), track.pt(), track.tpcSignal(), track.tofNSigmaPi());
            registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("ExpTpcdEdx_ptwise_pions"), track.pt(), tpcExpSignalPi, track.tofNSigmaPi());
            registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("ExpSigma_ptwise_pions"), track.pt(), track.tpcExpSigmaPi(), track.tofNSigmaPi());
          }
        }
        if (pid == PidKaons) {
          registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("TofTpcNsigma_kaons"), track.tpcNSigmaKa(), track.tofNSigmaKa(), track.pt());
          if (cfgGetdEdx) {
            double tpcExpSignalKa = track.tpcSignal() - (track.tpcNSigmaKa() * track.tpcExpSigmaKa());

            registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("TpcdEdx_ptwise_kaons"), track.pt(), track.tpcSignal(), track.tofNSigmaKa());
            registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("ExpTpcdEdx_ptwise_kaons"), track.pt(), tpcExpSignalKa, track.tofNSigmaKa());
            registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("ExpSigma_ptwise_kaons"), track.pt(), track.tpcExpSigmaKa(), track.tofNSigmaKa());
          }
        }
        if (pid == PidProtons) {
          registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("TofTpcNsigma_protons"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
          if (cfgGetdEdx) {
            double tpcExpSignalPr = track.tpcSignal() - (track.tpcNSigmaPr() * track.tpcExpSigmaPr());
            registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("TpcdEdx_ptwise_protons"), track.pt(), track.tpcSignal(), track.tofNSigmaPr());
            registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("ExpTpcdEdx_ptwise_protons"), track.pt(), tpcExpSignalPr, track.tofNSigmaPr());
            registry.fill(HIST("QA_PID/") + HIST(FillTimeName[ft]) + HIST("ExpSigma_ptwise_protons"), track.pt(), track.tpcExpSigmaPr(), track.tofNSigmaPr());
          }
        }
      }
    }
  }

  template <QAFillTime ft, typename TTrack>
  inline void fillTrackQA(TTrack const& track, const float vtxz)
  {
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
  template <QAFillTime ft, typename TCollision>
  inline void fillEventQA(TCollision const& collision, XAxis xaxis)
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
  }

  double getTimeSinceStartOfFill(uint64_t, int) { return 0.0; }

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs>>::iterator const& collision, aod::BCsWithTimestamps const&, GFWTracks const& tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    if (run != lastRun) {
      lastRun = run;
      LOGF(info, "run = %d", run);
    }
    loadCorrections(bc);

    registry.fill(HIST("eventQA/eventSel"), kFilteredEvent);
    if (!collision.sel8())
      return;
    registry.fill(HIST("eventQA/eventSel"), kSel8);
    registry.fill(HIST("eventQA/eventSel"), kOccupancy); // Add occupancy selection later

    const XAxis xaxis{.centrality = getCentrality(collision), .multiplicity = tracks.size()};
    if (cfgFillQA) {
      fillEventQA<kBefore>(collision, xaxis);
      registry.fill(HIST("eventQA/before/centrality"), xaxis.centrality);
      registry.fill(HIST("eventQA/before/multiplicity"), xaxis.multiplicity);
    }
    if (cfgUseAdditionalEventCut && !eventSelected(collision, xaxis.multiplicity, xaxis.centrality))
      return;
    if (cfgFillQA)
      fillEventQA<kAfter>(collision, xaxis);

    registry.fill(HIST("eventQA/after/centrality"), xaxis.centrality);
    registry.fill(HIST("eventQA/after/multiplicity"), xaxis.multiplicity);
    processCollision<kReco>(collision, tracks, xaxis, run);
  }
  PROCESS_SWITCH(FlowGfwV02, processData, "Process analysis for non-derived data", true);

  void processCFDerived(aod::CFCollision const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    int run = collision.runNumber();
    if (run != lastRun) {
      lastRun = run;
      LOGF(info, "run = %d", run);
    }
    loadCorrections(run);
    const XAxis xaxis{.centrality = collision.multiplicity(), .multiplicity = tracks.size()};

    registry.fill(HIST("eventQA/after/centrality"), xaxis.centrality);
    registry.fill(HIST("eventQA/after/multiplicity"), xaxis.multiplicity);

    // processCollision<kReco>(collision, tracks, xaxis, run);
  }
  PROCESS_SWITCH(FlowGfwV02, processCFDerived, "Process analysis for CF derived data", false);
  void processCFDerivedCorrected(aod::CFCollision const& collision, soa::Filtered<soa::Join<aod::CFTracks, aod::JWeights>> const& tracks)
  {
    int run = collision.runNumber();
    if (run != lastRun) {
      lastRun = run;
      LOGF(info, "run = %d", run);
    }
    const XAxis xaxis{.centrality = collision.multiplicity(), .multiplicity = tracks.size()};
    registry.fill(HIST("eventQA/after/centrality"), xaxis.centrality);
    registry.fill(HIST("eventQA/after/multiplicity"), xaxis.multiplicity);
    // processCollision<kReco>(collision, tracks, xaxis, run);
  }
  PROCESS_SWITCH(FlowGfwV02, processCFDerivedCorrected, "Process analysis for CF derived data with corrections", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGfwV02>(cfgc),
  };
}
