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

/// \file flowPbpbPikp.cxx
/// \brief PID flow using the generic framework
/// \author Preet Bhanjan Pati <preet.bhanjan.pati@cern.ch>

#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWConfig.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
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
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <TF1.h>
#include <TH1.h>
#include <TH3.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TString.h>

#include <sys/types.h>

#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, (DEFAULT), (HELP)}; // NOLINT(bugprone-macro-parentheses)

using namespace o2::analysis::genericframework;

template <typename T, typename P>
auto readMatrix(Array2D<T> const& mat, P& array)
{
  for (auto i = 0; i < static_cast<int>(mat.rows); ++i) {
    for (auto j = 0; j < static_cast<int>(mat.cols); ++j) {
      array[i][j] = mat(i, j);
    }
  }

  return;
}

static constexpr std::array<std::array<float, 3>, 20> LongArrayFloat = {{{{1.1, 2.1, 3.1}}, {{1.2, 2.2, 3.2}}, {{1.3, 2.3, 3.3}}, {{-1.1, -2.1, -3.1}}, {{-1.2, -2.2, -3.2}}, {{-1.3, -2.3, -3.3}}, {{1.1, 1.1, 1.1}}, {{1.2, 1.2, 1.2}}, {{1.3, 1.3, 1.3}}, {{-1.1, -1.1, -1.1}}, {{-1.2, -1.2, -1.2}}, {{-1.3, -1.3, -1.3}}, {{1.1, 1.1, 1.1}}, {{1.2, 1.2, 1.2}}, {{1.3, 1.3, 1.3}}, {{-1.1, -1.1, -1.1}}, {{-1.2, -1.2, -1.2}}, {{-1.3, -1.3, -1.3}}, {{1.1, 1.1, 1.1}}, {{1.2, 1.2, 1.2}}}};
static constexpr std::array<std::array<int, 3>, 20> LongArrayInt = {{{{1, 2, 3}}, {{1, 2, 3}}, {{1, 2, 3}}, {{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}, {{1, 1, 1}}, {{1, 1, 1}}, {{1, 1, 1}}, {{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}, {{1, 1, 1}}, {{1, 1, 1}}, {{1, 1, 1}}, {{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}, {{1, 1, 1}}, {{1, 1, 1}}}};

struct FlowPbpbPikp {
  o2::aod::ITSResponse itsResponse;
  Service<ccdb::BasicCCDBManager> ccdb{};

  struct : ConfigurableGroup {
    std::string prefix = "GeneralConfigs";
    O2_DEFINE_CONFIGURABLE(cfgQABasic, bool, true, "Enable QA histograms for event and track selection")
    O2_DEFINE_CONFIGURABLE(cfgRefPtMin, float, 0.2f, "Minimal pT for ref tracks")
    O2_DEFINE_CONFIGURABLE(cfgRefPtMax, float, 3.0f, "Maximal pT for ref tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
    O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
    O2_DEFINE_CONFIGURABLE(cfgStrictTrackCounter, bool, false, "Strict track counter for multiplicity correlation cut, counts only tracks that pass all cuts and are used in the correlation")
  } GeneralConfigs;

  struct : ConfigurableGroup {
    std::string prefix = "Weights";
    O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
    O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
    O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, true, "Fill and output NUA weights")
    O2_DEFINE_CONFIGURABLE(cfgOutputRunByRun, bool, true, "Fill and output NUA weights run by run")
    O2_DEFINE_CONFIGURABLE(cfgUseWeightPhiEtaVtxz, bool, false, "Use Phi, Eta, VertexZ dependent NUA weights")
    O2_DEFINE_CONFIGURABLE(cfgUseWeightPhiPtCent, bool, false, "Use Phi, Pt, Centrality dependent NUA weights")
    O2_DEFINE_CONFIGURABLE(cfgUseWeightPhiEtaPt, bool, true, "Use Phi, Eta, Pt dependent NUA weights")
  } Weights;

  struct : ConfigurableGroup {
    std::string prefix = "EvSels";
    O2_DEFINE_CONFIGURABLE(cfgGlobalCutVertex, float, 10.0f, "Accepted z-vertex range")
    O2_DEFINE_CONFIGURABLE(cfgCutOccupancyMax, int, 2000, "Maximum occupancy cut")
    O2_DEFINE_CONFIGURABLE(cfgCutOccupancyMin, int, 0, "Minimum occupancy cut")
    O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, true, "Use Nch for analysis instead of centrality")
    O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FV0A")
    O2_DEFINE_CONFIGURABLE(cfgCentTableUnavailable, bool, false, "if a dataset does not provide centrality information")
    Configurable<LabeledArray<int>> cfgUseEventCuts{"cfgUseEventCuts", {LongArrayInt.front().data(), 14, 1, {"Filtered Events", "Sel8", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "kNoCollInTimeRangeStandard", "kIsGoodITSLayersAll", "kIsGoodITSLayer0123", "kNoCollInRofStandard", "kNoHighMultCollInPrevRof", "Occupancy", "Multcorrelation", "T0AV0ACut"}, {"EvCuts"}}, "Labeled array (int) for various cuts on resonances"};
  } EvSels;

  struct : ConfigurableGroup {
    std::string prefix = "TrkSels";
    O2_DEFINE_CONFIGURABLE(cfgGlobalPtMin, float, 0.2f, "minimum accepted track pT")
    O2_DEFINE_CONFIGURABLE(cfgGlobalPtMax, float, 10.0f, "maximum accepted track pT")
    O2_DEFINE_CONFIGURABLE(cfgGlobalEta, float, 0.8f, "Eta cut")
    O2_DEFINE_CONFIGURABLE(cfgGlobalChi2prTPCcls, float, 2.5f, "Maximum chi2 per TPC clusters")
    O2_DEFINE_CONFIGURABLE(cfgGlobalDCAz, float, 2.0f, "Maximum DCA to vertex z")
    Configurable<std::vector<std::string>> cfgTrackCutsDCAxy{"cfgTrackCutsDCAxy", std::vector<std::string>{"[0]*(0.0026+0.005/(x^1.01))", "[0]*(0.0026+0.005/(x^1.01))"}, "Functional form of pt-dependent DCAxy cut"};
    Configurable<std::vector<std::string>> cfgTrackCutsDCAz{"cfgTrackCutsDCAz", std::vector<std::string>{"", ""}, "Functional form of pt-dependent DCAz cut"};
    Configurable<LabeledArray<float>> cfgTrackCuts{"cfgTrackCuts", {LongArrayFloat.front().data(), 10, 2, {"PtMin", "PtMax", "EtaCut", "Chi2PrTpcCls", "TpcCluster", "TpcCrossedRows", "ItsCluster", "DCAz", "DCAzNsigma", "DCAxyNsigma"}, {"General", "SelNch"}}, "Labeled array for track selections"};
  } TrkSels;

  struct : ConfigurableGroup {
    std::string prefix = "PIDConfigs";
    O2_DEFINE_CONFIGURABLE(cfgQANsigma, bool, false, "Get QA histograms for PID Nsigma")
    O2_DEFINE_CONFIGURABLE(cfgQAdEdx, bool, false, "Get QA histograms for PID (TpcdEdx, ExpTpc, ExpSigma) vs TOF Nsigma")
    O2_DEFINE_CONFIGURABLE(cfgQABetadEdx, bool, false, "Get QA histograms for TpcdEdx vs TofBeta")
    O2_DEFINE_CONFIGURABLE(cfgPIDTofPtCut, float, 0.4f, "Minimum pt to use TOF N-sigma")
    O2_DEFINE_CONFIGURABLE(cfgPIDUseIts, bool, false, "Use ITS PID for particle identification")
    O2_DEFINE_CONFIGURABLE(cfgPIDUseRejection, bool, true, "True: use exclusion exclusion criteria for PID determination, false: don't use exclusion")
    O2_DEFINE_CONFIGURABLE(cfgPIDNclsPIDCut, int, 70, "Number of clusters used for PID selection")
    O2_DEFINE_CONFIGURABLE(cfgPIDUseNclsPID, bool, false, "Switch for using cfgPIDNclsPIDCut")
    Configurable<LabeledArray<float>> cfgPIDnSigmas{"cfgPIDnSigmas", {LongArrayFloat.front().data(), 6, 3, {"UpCut_pi", "UpCut_ka", "UpCut_pr", "LowCut_pi", "LowCut_ka", "LowCut_pr"}, {"TPC", "TOF", "ITS"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};
  } PIDConfigs;

  struct : ConfigurableGroup {
    std::string prefix = "FuncParas";
    O2_DEFINE_CONFIGURABLE(cfgMultCentHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 10.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultCentLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultT0CCutEnabled, bool, false, "Enable Global multiplicity vs T0C centrality cut")
    Configurable<std::vector<double>> cfgMultT0CCutPars{"cfgMultT0CCutPars", std::vector<double>{143.04, -4.58368, 0.0766055, -0.000727796, 2.86153e-06, 23.3108, -0.36304, 0.00437706, -4.717e-05, 1.98332e-07}, "Global multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultPVT0CCutEnabled, bool, false, "Enable PV multiplicity vs T0C centrality cut")
    Configurable<std::vector<double>> cfgMultPVT0CCutPars{"cfgMultPVT0CCutPars", std::vector<double>{195.357, -6.15194, 0.101313, -0.000955828, 3.74793e-06, 30.0326, -0.43322, 0.00476265, -5.11206e-05, 2.13613e-07}, "PV multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultMultPVHighCutFunction, std::string, "[0]+[1]*x + 5.*([2]+[3]*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultPVLowCutFunction, std::string, "[0]+[1]*x - 5.*([2]+[3]*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalPVCutEnabled, bool, false, "Enable global multiplicity vs PV multiplicity cut")
    Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars", std::vector<double>{-0.140809, 0.734344, 2.77495, 0.0165935}, "PV multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0AHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 4.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0ALowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0ACutEnabled, bool, false, "Enable global multiplicity vs V0A multiplicity cut")
    Configurable<std::vector<double>> cfgMultMultV0ACutPars{"cfgMultMultV0ACutPars", std::vector<double>{534.893, 184.344, 0.423539, -0.00331436, 5.34622e-06, 871.239, 53.3735, -0.203528, 0.000122758, 5.41027e-07}, "Global multiplicity vs V0A multiplicity cut parameter values"};
    std::vector<double> multT0CCutPars;
    std::vector<double> multPVT0CCutPars;
    std::vector<double> multGlobalPVCutPars;
    std::vector<double> multMultV0ACutPars;
    std::unique_ptr<TF1> fMultPVT0CCutLow = nullptr;
    std::unique_ptr<TF1> fMultPVT0CCutHigh = nullptr;
    std::unique_ptr<TF1> fMultT0CCutLow = nullptr;
    std::unique_ptr<TF1> fMultT0CCutHigh = nullptr;
    std::unique_ptr<TF1> fMultGlobalPVCutLow = nullptr;
    std::unique_ptr<TF1> fMultGlobalPVCutHigh = nullptr;
    std::unique_ptr<TF1> fMultMultV0ACutLow = nullptr;
    std::unique_ptr<TF1> fMultMultV0ACutHigh = nullptr;
    std::unique_ptr<TF1> fT0AV0AMean = nullptr;
    std::unique_ptr<TF1> fT0AV0ASigma = nullptr;
    std::unique_ptr<TF1> fPtDepDCAxy = nullptr;
    std::unique_ptr<TF1> fPtDepDCAxyForNch = nullptr;
    std::unique_ptr<TF1> fPtDepDCAz = nullptr;
    std::unique_ptr<TF1> fPtDepDCAzForNch = nullptr;
    O2_DEFINE_CONFIGURABLE(cfgV0AT0Acut, int, 5, "V0AT0A cut")
  } FuncParas;

  struct : ConfigurableGroup {
    std::string prefix = "TrkDensCorr";
    O2_DEFINE_CONFIGURABLE(cfgTrackDensityCorrUse, bool, false, "Use track density efficiency correction")
    Configurable<std::vector<double>> cfgTrackDensityP0{"cfgTrackDensityP0", std::vector<double>{0.7217476707, 0.7384792571, 0.7542625668, 0.7640680200, 0.7701951667, 0.7755299053, 0.7805901710, 0.7849446786, 0.7957356586, 0.8113039262, 0.8211968966, 0.8280558878, 0.8329342135}, "parameter 0 for track density efficiency correction"};
    Configurable<std::vector<double>> cfgTrackDensityP1{"cfgTrackDensityP1", std::vector<double>{-2.169488e-05, -2.191913e-05, -2.295484e-05, -2.556538e-05, -2.754463e-05, -2.816832e-05, -2.846502e-05, -2.843857e-05, -2.705974e-05, -2.477018e-05, -2.321730e-05, -2.203315e-05, -2.109474e-05}, "parameter 1 for track density efficiency correction"};
  } TrkDensCorr;

  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN08", "refP08", "full"}, {-0.8, 0.4, -0.8}, {-0.4, 0.8, 0.8}, {0, 0, 0}, {1, 1, 1}}, "Configurations for GFW regions (Name, minEta, maxEta, pT-diff, BitMask)"};
  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"full {2 -2}", "full {2 -2}", "full {2 -2}"}, {"ChFull22", "PiFull22", "KaFull22"}, {0, 0, 0}, {0, 0, 0}}, "Configurations for each correlation to calculate"};

  // Axes for general QA
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {16, -0.8, 0.8}, "eta axis for histograms"};
  ConfigurableAxis axisPtForQA{"axisPtForQA", {98, 0.2, 10.0}, "pt axis for QA histograms"};
  ConfigurableAxis axisMultForQA{"axisMultForQA", {100, 0, 1000}, "multiplicity axis for QA histograms"};
  ConfigurableAxis axisAllParticles{"axisAllParticles", {4, -0.5, 3.5}, "particle axis for general QA (ch, pi, ka, pr)"};

  // Axes for correlations
  ConfigurableAxis axisMultForCorr{"axisMultForCorr", {VARIABLE_WIDTH, 0, 5, 10, 20, 50, 100, 150, 200, 300, 500}, "multiplicity axis for correlation histograms"};
  ConfigurableAxis axisCentForCorr{"axisCentForCorr", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for correlation histograms"};
  ConfigurableAxis axisPtForCorr{"axisPtForCorr", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt trigger axis for correlation histograms"};

  // Axes for PID QA
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisNsigmaITS{"axisNsigmaITS", {80, -5, 5}, "nsigmaITS axis"};
  ConfigurableAxis axisTpcSignal{"axisTpcSignal", {250, 0, 250}, "dEdx axis for TPC"};
  ConfigurableAxis axisTofBeta{"axisTofBeta", {200, 0, 2}, "axis for TOF beta"};
  ConfigurableAxis axisSigma{"axisSigma", {200, 0, 20}, "sigma axis for TPC"};
  ConfigurableAxis axisPIDParticles{"axisPIDParticles", {3, 0.5, 3.5}, "particle axis for PID QA (pi, ka, pr)"};

  struct GFWMemberCache {
    GFWRegions regions;
    GFWCorrConfigs configs;
  } GfwMemberCache;

  Filter collisionFilter = (nabs(aod::collision::posZ) < EvSels.cfgGlobalCutVertex);
  Filter trackFilter = (nabs(aod::track::eta) < TrkSels.cfgGlobalEta) && (aod::track::pt > TrkSels.cfgGlobalPtMin) && (aod::track::pt < TrkSels.cfgGlobalPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == static_cast<uint8_t>(1))) && (aod::track::tpcChi2NCl < TrkSels.cfgGlobalChi2prTPCcls) && (nabs(aod::track::dcaZ) < TrkSels.cfgGlobalDCAz);

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;
  using AodTracksWithoutBayes = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using AodTracksMcLabels = soa::Filtered<soa::Join<aod::Tracks, aod::McTrackLabels, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;

  // OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<FlowContainer> fFC{"FlowContainer"};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::map<int, std::vector<std::shared_ptr<TH3>>> th3sList;

  // Generic framework intializations
  std::unique_ptr<GFW> fGFW{std::make_unique<GFW>()};
  std::vector<GFW::CorrConfig> corrConfigs;
  TAxis* fPtAxis = nullptr;

  // Initializer for random number
  TRandom3* fRandom = new TRandom3(0);

  // Enum definitions
  enum OutputSpecies {
    UseCharged = 0,
    UsePions,
    UseKaons,
    UseProtons,
    UseRef,
    NumCountTypes
  };

  enum PiKpArrayIndex {
    IndexPionUp = 0,
    IndexKaonUp,
    IndexProtonUp,
    IndexPionLow,
    IndexKaonLow,
    IndexProtonLow
  };

  enum EventCutTypes {
    FilteredEvents = 0,
    AfterSel8,
    UseNoTimeFrameBorder,
    UseNoITSROFrameBorder,
    UseNoSameBunchPileup,
    UseGoodZvtxFT0vsPV,
    UseNoCollInTimeRangeStandard,
    UseGoodITSLayersAll,
    UseGoodITSLayer0123,
    UseNoCollInRofStandard,
    UseNoHighMultCollInPrevRof,
    UseOccupancy,
    UseMultCorrCut,
    UseT0AV0ACut,
    NumEventCuts
  };

  enum TrackCuts {
    TrkCutPtMin = 0,
    TrkCutPtMax,
    TrkCutEta,
    TrkCutChi2PrTpcCls,
    TrkCutTpcCluster,
    TrkCutTpcCrossedRows,
    TrkCutItsCluster,
    TrkCutDCAz,
    TrkCutDCAzNsigma,
    TrkCutDCAxyNsigma
  };

  enum FillTrackCuts {
    FillFilteredTrks = 1,
    FillCutPt,
    FillCutEta,
    FillCutChi2PrTpcCls,
    FillCutTpcCluster,
    FillCutTpcCrossedRows,
    FillCutItsCluster,
    FillCutDCAz,
    FillCutDCAxy
  };

  enum TrackCutGroup {
    UseGenTrkCuts = 0,
    UseNchSelCuts,
    NumTrackCutTypes
  };
  enum CentEstimators {
    UseCentFT0C = 0,
    UseCentFT0CVariant1,
    UseCentFT0M,
    UseCentFV0A,
    NumCentEstimators
  };
  enum DetectorType {
    UseTPC = 0,
    UseTOF,
    UseITS
  };
  enum EventCutType {
    EvCut1 = 0,
    NumEvCutTypes = 1
  };
  enum Stage {
    Before = 0,
    After
  };

  // Vectors for labelled arrays
  std::array<std::array<float, 2>, 10> trackCuts{};
  std::array<std::array<int, 1>, 14> eventCuts{};
  std::array<std::array<float, 3>, 6> nSigmaVals{};

  // For run by run analysis
  int lastRunNumer = -1;
  std::vector<int> runNumbers;
  std::vector<GFWWeights*> mAcceptance;
  bool correctionsLoaded = false;

  TH1D* mEfficiency = nullptr;

  // Local track density correction - Copy from flowTask.cxx
  std::vector<TF1*> funcEff;
  TH1D* hFindPtBin = nullptr;
  TF1* funcV2 = nullptr;
  TF1* funcV3 = nullptr;
  TF1* funcV4 = nullptr;

  void init(InitContext const&)
  {
    readMatrix(TrkSels.cfgTrackCuts->getData(), trackCuts);
    readMatrix(EvSels.cfgUseEventCuts->getData(), eventCuts);
    readMatrix(PIDConfigs.cfgPIDnSigmas->getData(), nSigmaVals);

    if (EvSels.cfgCentTableUnavailable && !EvSels.cfgUseNch) {
      LOGF(fatal, "Centrality table is unavailable, cannot select collisions by centrality");
    }

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    LOGF(info, "flowGenericFramework::init()");
    GfwMemberCache.regions.SetNames(cfgRegions->GetNames());
    GfwMemberCache.regions.SetEtaMin(cfgRegions->GetEtaMin());
    GfwMemberCache.regions.SetEtaMax(cfgRegions->GetEtaMax());
    GfwMemberCache.regions.SetpTDifs(cfgRegions->GetpTDifs());
    GfwMemberCache.regions.SetBitmasks(cfgRegions->GetBitmasks());
    GfwMemberCache.configs.SetCorrs(cfgCorrConfig->GetCorrs());
    GfwMemberCache.configs.SetHeads(cfgCorrConfig->GetHeads());
    GfwMemberCache.configs.SetpTDifs(cfgCorrConfig->GetpTDifs());
    GfwMemberCache.configs.SetpTCorrMasks(cfgCorrConfig->GetpTCorrMasks());
    GfwMemberCache.regions.Print();
    GfwMemberCache.configs.Print();

    // Axis definitions
    std::string hCentTitle = "Centrality distribution, Estimator " + std::to_string(EvSels.cfgCentEstimator);
    const AxisSpec axisCentForQA{100, 0, 100, "Centrality (%)"};
    AxisSpec axisEventClass = {(EvSels.cfgUseNch) ? axisMultForCorr : axisCentForCorr};

    histos.add("hEventCount", "Number of Events;; Count", {HistType::kTH1D, {{14, -0.5, 13.5}}});
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(FilteredEvents + 1, "Filtered events");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(AfterSel8 + 1, "After sel8");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoSameBunchPileup + 1, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseGoodITSLayer0123 + 1, "kIsGoodITSLayer0123");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoCollInRofStandard + 1, "kNoCollInRofStandard");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoHighMultCollInPrevRof + 1, "kNoHighMultCollInPrevRof");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseOccupancy + 1, "Occupancy Cut");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseMultCorrCut + 1, "MultCorrelation Cut");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseT0AV0ACut + 1, "T0AV0A cut");

    histos.add("hTrackCount", "Number of Tracks;; Count", {HistType::kTH1D, {{9, -0.5, 8.5}}});
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(FillFilteredTrks, "Filtered track");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(FillCutPt, "p_{T} Cut");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(FillCutEta, "#eta Cut");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(FillCutChi2PrTpcCls, "#chi^2/nTpcCls Cut");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(FillCutTpcCluster, "TPC Clusters");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(FillCutTpcCrossedRows, "TPC Crossed Rows");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(FillCutItsCluster, "ITS Clusters");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(FillCutDCAz, "DCAz Cut");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(FillCutDCAxy, "DCAxy Cut");

    if (doprocessGFW) {
      if (Weights.cfgOutputNUAWeights && !Weights.cfgOutputRunByRun) {
        histos.add<TH3>("NUA/hPhiEtaVtxz_ref", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
        histos.add<TH3>("NUA/hPhiEtaVtxz_ch", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
        histos.add<TH3>("NUA/hPhiEtaVtxz_pi", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
        histos.add<TH3>("NUA/hPhiEtaVtxz_ka", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
        histos.add<TH3>("NUA/hPhiEtaVtxz_pr", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});

        histos.add<TH3>("NUA/hPhiPtCent_ref", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPtForCorr, axisCentForCorr}});
        histos.add<TH3>("NUA/hPhiPtCent_ch", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPtForCorr, axisCentForCorr}});
        histos.add<TH3>("NUA/hPhiPtCent_pi", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPtForCorr, axisCentForCorr}});
        histos.add<TH3>("NUA/hPhiPtCent_ka", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPtForCorr, axisCentForCorr}});
        histos.add<TH3>("NUA/hPhiPtCent_pr", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPtForCorr, axisCentForCorr}});

        histos.add<TH3>("NUA/hPhiEtaPt_ref", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPtForCorr}});
        histos.add<TH3>("NUA/hPhiEtaPt_ch", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPtForCorr}});
        histos.add<TH3>("NUA/hPhiEtaPt_pi", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPtForCorr}});
        histos.add<TH3>("NUA/hPhiEtaPt_ka", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPtForCorr}});
        histos.add<TH3>("NUA/hPhiEtaPt_pr", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPtForCorr}});
      } // end of output weight histograms condition
      if (!Weights.cfgAcceptance.value.empty()) {
        if (Weights.cfgUseWeightPhiEtaVtxz) {
          histos.add<TH3>("PhiCorrected/hPhiEtaVtxz_pi_corrd", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
          histos.add<TH3>("PhiCorrected/hPhiEtaVtxz_ka_corrd", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
          histos.add<TH3>("PhiCorrected/hPhiEtaVtxz_pr_corrd", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
        }
        if (Weights.cfgUseWeightPhiPtCent) {
          histos.add<TH3>("PhiCorrected/hPhiPtCent_pi_corrd", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPtForCorr, axisCentForCorr}});
          histos.add<TH3>("PhiCorrected/hPhiPtCent_ka_corrd", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPtForCorr, axisCentForCorr}});
          histos.add<TH3>("PhiCorrected/hPhiPtCent_pr_corrd", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPtForCorr, axisCentForCorr}});
        }
        if (Weights.cfgUseWeightPhiEtaPt) {
          histos.add<TH3>("PhiCorrected/hPhiEtaPt_pi_corrd", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPtForCorr}});
          histos.add<TH3>("PhiCorrected/hPhiEtaPt_ka_corrd", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPtForCorr}});
          histos.add<TH3>("PhiCorrected/hPhiEtaPt_pr_corrd", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPtForCorr}});
        }
      } // end of acceptance corrected phi histograms condition

      if (GeneralConfigs.cfgQABasic) {
        histos.add("hPhi", "#phi", {HistType::kTH2D, {axisPhi, axisAllParticles}});
        histos.add("hPhiWeighted", "#phi Corrected", {HistType::kTH2D, {axisPhi, axisAllParticles}});
        histos.add("hEta", "#eta", {HistType::kTH2D, {axisEta, axisAllParticles}});
        histos.add("hEtaCorrected", "#eta Corrected", {HistType::kTH2D, {axisEta, axisAllParticles}});
        histos.add("hPt", "p_{T}", {HistType::kTHnSparseF, {axisPtForQA, axisAllParticles, axisEventClass}});
        histos.add("hPtCorrected", "p_{T} Corrected", {HistType::kTHnSparseF, {axisPtForQA, axisAllParticles, axisEventClass}});
      }

      histos.add("hNch", "", {HistType::kTH1D, {axisMultForQA}});
      histos.add("hNchCorrected", "", {HistType::kTH1D, {axisMultForQA}});
      histos.add("hCent", hCentTitle.c_str(), {HistType::kTH1D, {axisCentForQA}});
      histos.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});

      o2::framework::AxisSpec axis = axisPtForCorr;
      int nPtBins = axis.binEdges.size() - 1;
      double* ptBins = axis.binEdges.data();
      fPtAxis = new TAxis(nPtBins, ptBins);

      // Iniitializing object for GFW
      auto* oba = new TObjArray();

      // Defining the regions
      for (auto i(0); i < GfwMemberCache.regions.GetSize(); i++) {
        fGFW->AddRegion(GfwMemberCache.regions.GetNames().at(i), GfwMemberCache.regions.GetEtaMin().at(i), GfwMemberCache.regions.GetEtaMax().at(i), (GfwMemberCache.regions.GetpTDifs().at(i) != 0) ? nPtBins + 1 : 1, GfwMemberCache.regions.GetBitmasks().at(i));
      }

      // Defining the correlators
      std::string suffix = "_ptDiff";
      for (auto i = 0; i < GfwMemberCache.configs.GetSize(); i++) {
        corrConfigs.push_back(fGFW->GetCorrelatorConfig(GfwMemberCache.configs.GetCorrs().at(i), GfwMemberCache.configs.GetHeads().at(i), GfwMemberCache.configs.GetpTDifs().at(i)));
        if (GfwMemberCache.configs.GetpTDifs().at(i) != 0) {
          for (auto j = 1; j <= nPtBins; j++) {
            std::string index = Form("_pt_%i", j);
            oba->Add(new TNamed(GfwMemberCache.configs.GetHeads().at(i) + index, GfwMemberCache.configs.GetHeads().at(i) + suffix));
            std::string corrName = GfwMemberCache.configs.GetHeads().at(i) + index;
            LOGF(info, Form("Correlator to be calculated: %s", corrName.c_str()));
          }
        } else {
          oba->Add(new TNamed(GfwMemberCache.configs.GetHeads().at(i), GfwMemberCache.configs.GetHeads().at(i)));
          LOGF(info, Form("Correlator to be calculated: %s", GfwMemberCache.configs.GetHeads().at(i).c_str()));
        }
      }
      if (corrConfigs.empty()) {
        LOGF(error, "Configuration contains vectors of different size - check the GFWCorrConfig configurable");
      }
      fGFW->CreateRegions();

      fFC.setObject(new FlowContainer("FlowContainer"));
      fFC->SetXAxis(fPtAxis);
      fFC->Initialize(oba, axisEventClass, GeneralConfigs.cfgNbootstrap);

      delete oba;
    } // end of doprocessGFW

    if (doprocessGFW || doprocessPIDQA || doprocessMCPID) {
      if (PIDConfigs.cfgQANsigma) {
        if (!PIDConfigs.cfgPIDUseIts) { // TPC-TOF PID QA hists
          histos.add("TofTpcNsigma_beforeCut", "", {HistType::kTHnSparseD, {{axisNsigmaTPC, axisNsigmaTOF, axisEventClass, axisPtForCorr, axisPIDParticles}}});
          histos.add("TofTpcNsigma_afterCut", "", {HistType::kTHnSparseD, {{axisNsigmaTPC, axisNsigmaTOF, axisEventClass, axisPtForCorr, axisPIDParticles}}});
        }
        if (PIDConfigs.cfgPIDUseIts) { // ITS-TOF PID QA hists
          histos.add("TofItsNsigma_beforeCut", "", {HistType::kTHnSparseD, {{axisNsigmaITS, axisNsigmaTOF, axisEventClass, axisPtForCorr, axisPIDParticles}}});
          histos.add("TofItsNsigma_afterCut", "", {HistType::kTHnSparseD, {{axisNsigmaITS, axisNsigmaTOF, axisEventClass, axisPtForCorr, axisPIDParticles}}});
        }
      } // end of PID nsigma QA hists
    }

    if (doprocessPIDQA) {
      if (PIDConfigs.cfgQAdEdx) {
        histos.add("TpcdEdx_beforeCut", "", {HistType::kTHnSparseD, {{axisTpcSignal, axisNsigmaTOF, axisEventClass, axisPtForCorr, axisPIDParticles}}});
        histos.add("ExpTpcdEdx_beforeCut", "", {HistType::kTHnSparseD, {{axisTpcSignal, axisNsigmaTOF, axisEventClass, axisPtForCorr, axisPIDParticles}}});
        histos.add("ExpSigma_beforeCut", "", {HistType::kTHnSparseD, {{axisSigma, axisNsigmaTOF, axisEventClass, axisPtForCorr, axisPIDParticles}}});

        histos.add("TpcdEdx_afterCut", "", {HistType::kTHnSparseD, {{axisTpcSignal, axisNsigmaTOF, axisEventClass, axisPtForCorr, axisPIDParticles}}});
        histos.add("ExpTpcdEdx_afterCut", "", {HistType::kTHnSparseD, {{axisTpcSignal, axisNsigmaTOF, axisEventClass, axisPtForCorr, axisPIDParticles}}});
        histos.add("ExpSigma_afterCut", "", {HistType::kTHnSparseD, {{axisSigma, axisNsigmaTOF, axisEventClass, axisPtForCorr, axisPIDParticles}}});
      } // end of PID dEdx QA hists

      if (PIDConfigs.cfgQABetadEdx) {
        histos.add("TpcdEdx_TofBeta_beforeCut", "", {HistType::kTHnSparseD, {{axisTpcSignal, axisTofBeta, axisEventClass, axisPtForCorr, axisPIDParticles}}});
        histos.add("ExpTpcdEdx_TofBeta_beforeCut", "", {HistType::kTHnSparseD, {{axisTpcSignal, axisTofBeta, axisEventClass, axisPtForCorr, axisPIDParticles}}});
        histos.add("ExpSigma_TofBeta_beforeCut", "", {HistType::kTHnSparseD, {{axisSigma, axisTofBeta, axisEventClass, axisPtForCorr, axisPIDParticles}}});

        histos.add("TpcdEdx_TofBeta_afterCut", "", {HistType::kTHnSparseD, {{axisTpcSignal, axisTofBeta, axisEventClass, axisPtForCorr, axisPIDParticles}}});
        histos.add("ExpTpcdEdx_TofBeta_afterCut", "", {HistType::kTHnSparseD, {{axisTpcSignal, axisTofBeta, axisEventClass, axisPtForCorr, axisPIDParticles}}});
        histos.add("ExpSigma_TofBeta_afterCut", "", {HistType::kTHnSparseD, {{axisSigma, axisTofBeta, axisEventClass, axisPtForCorr, axisPIDParticles}}});
      }
    }

    if (doprocessMCPID) {
      histos.add("hNTruePositives", "hNTruePositives", {HistType::kTHnSparseD, {axisEventClass, axisPtForCorr, axisPIDParticles}}); // Fraction of particles that are correctly identified as the given PID particle
      histos.add("hNAllPositives", "hNAllPositives", {HistType::kTHnSparseD, {axisEventClass, axisPtForCorr, axisPIDParticles}});   // Fraction of particles identified as the given PID particle
      histos.add("hNAllTrues", "hNAllTrues", {HistType::kTHnSparseD, {axisEventClass, axisPtForCorr, axisPIDParticles}});           // All true PID particles from MC

      histos.add("TpcdEdx_AllTrues_physprim", "", {HistType::kTHnSparseD, {{axisTpcSignal, axisTofBeta, axisPtForCorr, axisPIDParticles}}});
      histos.add("TpcdEdx_AllPositives_physprim", "", {HistType::kTHnSparseD, {{axisTpcSignal, axisTofBeta, axisPtForCorr, axisPIDParticles}}});
    }

    // Multiplicity correlation cuts
    if (eventCuts[UseMultCorrCut][EvCut1] != 0) {
      FuncParas.multT0CCutPars = FuncParas.cfgMultT0CCutPars;
      FuncParas.multPVT0CCutPars = FuncParas.cfgMultPVT0CCutPars;
      FuncParas.multGlobalPVCutPars = FuncParas.cfgMultGlobalPVCutPars;
      FuncParas.multMultV0ACutPars = FuncParas.cfgMultMultV0ACutPars;

      FuncParas.fMultPVT0CCutLow = std::make_unique<TF1>("fMultPVT0CCutLow", FuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      FuncParas.fMultPVT0CCutLow->SetParameters(FuncParas.multPVT0CCutPars.data());
      FuncParas.fMultPVT0CCutHigh = std::make_unique<TF1>("fMultPVT0CCutHigh", FuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      FuncParas.fMultPVT0CCutHigh->SetParameters(FuncParas.multPVT0CCutPars.data());

      FuncParas.fMultT0CCutLow = std::make_unique<TF1>("fMultT0CCutLow", FuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      FuncParas.fMultT0CCutLow->SetParameters(FuncParas.multT0CCutPars.data());
      FuncParas.fMultT0CCutHigh = std::make_unique<TF1>("fMultT0CCutHigh", FuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      FuncParas.fMultT0CCutHigh->SetParameters(FuncParas.multT0CCutPars.data());

      FuncParas.fMultGlobalPVCutLow = std::make_unique<TF1>("fMultGlobalPVCutLow", FuncParas.cfgMultMultPVLowCutFunction->c_str(), 0, 4000);
      FuncParas.fMultGlobalPVCutLow->SetParameters(FuncParas.multGlobalPVCutPars.data());
      FuncParas.fMultGlobalPVCutHigh = std::make_unique<TF1>("fMultGlobalPVCutHigh", FuncParas.cfgMultMultPVHighCutFunction->c_str(), 0, 4000);
      FuncParas.fMultGlobalPVCutHigh->SetParameters(FuncParas.multGlobalPVCutPars.data());

      FuncParas.fMultMultV0ACutLow = std::make_unique<TF1>("fMultMultV0ACutLow", FuncParas.cfgMultMultV0ALowCutFunction->c_str(), 0, 4000);
      FuncParas.fMultMultV0ACutLow->SetParameters(FuncParas.multMultV0ACutPars.data());
      FuncParas.fMultMultV0ACutHigh = std::make_unique<TF1>("fMultMultV0ACutHigh", FuncParas.cfgMultMultV0AHighCutFunction->c_str(), 0, 4000);
      FuncParas.fMultMultV0ACutHigh->SetParameters(FuncParas.multMultV0ACutPars.data());
    }
    if (eventCuts[UseT0AV0ACut][EvCut1] != 0) {
      FuncParas.fT0AV0AMean = std::make_unique<TF1>("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      FuncParas.fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      FuncParas.fT0AV0ASigma = std::make_unique<TF1>("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      FuncParas.fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

    if (!TrkSels.cfgTrackCutsDCAxy.value[UseGenTrkCuts].empty()) {
      FuncParas.fPtDepDCAxy = std::make_unique<TF1>("ptDepDCAxy", TrkSels.cfgTrackCutsDCAxy.value[UseGenTrkCuts].c_str(), 0.001, 1000);
      FuncParas.fPtDepDCAxy->SetParameter(0, trackCuts[TrkCutDCAxyNsigma][UseGenTrkCuts]);
      LOGF(info, "DCAxy pt-dependence function: %0.1f * %s", trackCuts[TrkCutDCAxyNsigma][UseGenTrkCuts], TrkSels.cfgTrackCutsDCAxy.value[UseGenTrkCuts].c_str());
    }
    if (!TrkSels.cfgTrackCutsDCAxy.value[UseNchSelCuts].empty()) {
      FuncParas.fPtDepDCAxyForNch = std::make_unique<TF1>("ptDepDCAxyForNch", TrkSels.cfgTrackCutsDCAxy.value[UseNchSelCuts].c_str(), 0.001, 1000);
      FuncParas.fPtDepDCAxyForNch->SetParameter(0, trackCuts[TrkCutDCAxyNsigma][UseNchSelCuts]);
      LOGF(info, "DCAxy pt-dependence function for Nch: %0.1f * %s", trackCuts[TrkCutDCAxyNsigma][UseNchSelCuts], TrkSels.cfgTrackCutsDCAxy.value[UseNchSelCuts].c_str());
    }

    if (!TrkSels.cfgTrackCutsDCAz.value[UseGenTrkCuts].empty()) {
      FuncParas.fPtDepDCAz = std::make_unique<TF1>("ptDepDCAz", TrkSels.cfgTrackCutsDCAz.value[UseGenTrkCuts].c_str(), 0.001, 1000);
      FuncParas.fPtDepDCAz->SetParameter(0, trackCuts[TrkCutDCAzNsigma][UseGenTrkCuts]);
      LOGF(info, "DCAz pt-dependence function: %0.1f * %s", trackCuts[TrkCutDCAzNsigma][UseGenTrkCuts], TrkSels.cfgTrackCutsDCAz.value[UseGenTrkCuts].c_str());
    }
    if (!TrkSels.cfgTrackCutsDCAz.value[UseNchSelCuts].empty()) {
      FuncParas.fPtDepDCAzForNch = std::make_unique<TF1>("ptDepDCAzForNch", TrkSels.cfgTrackCutsDCAz.value[UseNchSelCuts].c_str(), 0.001, 1000);
      FuncParas.fPtDepDCAzForNch->SetParameter(0, trackCuts[TrkCutDCAzNsigma][UseNchSelCuts]);
      LOGF(info, "DCAz pt-dependence function for Nch: %0.1f * %s", trackCuts[TrkCutDCAzNsigma][UseNchSelCuts], TrkSels.cfgTrackCutsDCAz.value[UseNchSelCuts].c_str());
    }

    if (TrkDensCorr.cfgTrackDensityCorrUse) {
      std::vector<double> pTEffBins = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0};
      hFindPtBin = new TH1D("hFindPtBin", "hFindPtBin", pTEffBins.size() - 1, pTEffBins.data());
      funcEff.resize(pTEffBins.size() - 1);
      // LHC24g3 Eff
      std::vector<double> f1p0 = TrkDensCorr.cfgTrackDensityP0;
      std::vector<double> f1p1 = TrkDensCorr.cfgTrackDensityP1;
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

  } // End of init()

  template <typename TCollision>
  float getCentrality(TCollision const& collision)
  {
    float cent = 0.0;
    switch (EvSels.cfgCentEstimator) {
      case UseCentFT0C:
        cent = collision.centFT0C();
        break;
      case UseCentFT0CVariant1:
        cent = collision.centFT0CVariant1();
        break;
      case UseCentFT0M:
        cent = collision.centFT0M();
        break;
      case UseCentFV0A:
        cent = collision.centFV0A();
        break;
      default:
        cent = collision.centFT0C();
    }
    return cent;
  }

  template <typename TTrack>
  bool trackSelected(TTrack const& track)
  {
    histos.fill(HIST("hTrackCount"), FillFilteredTrks - 1); // Filtered tracks

    if ((track.pt() <= trackCuts[TrkCutPtMin][UseGenTrkCuts]) || (track.pt() >= trackCuts[TrkCutPtMax][UseGenTrkCuts])) {
      return false;
    }
    histos.fill(HIST("hTrackCount"), FillCutPt - 1); // After pT cut

    if (std::abs(track.eta()) > trackCuts[TrkCutEta][UseGenTrkCuts]) {
      return false;
    }
    histos.fill(HIST("hTrackCount"), FillCutEta - 1); // After eta cut

    if (track.tpcChi2NCl() > trackCuts[TrkCutChi2PrTpcCls][UseGenTrkCuts]) {
      return false;
    }
    histos.fill(HIST("hTrackCount"), FillCutChi2PrTpcCls - 1); // After chi2 per TPC NCls cut

    if (track.tpcNClsFound() < trackCuts[TrkCutTpcCluster][UseGenTrkCuts]) {
      return false;
    }
    histos.fill(HIST("hTrackCount"), FillCutTpcCluster - 1); // After TPC Clusters cut

    if (track.tpcNClsCrossedRows() < trackCuts[TrkCutTpcCrossedRows][UseGenTrkCuts]) {
      return false;
    }
    histos.fill(HIST("hTrackCount"), FillCutTpcCrossedRows - 1); // After TPC Crossed Rows cut

    if (track.itsNCls() < trackCuts[TrkCutItsCluster][UseGenTrkCuts]) {
      return false;
    }
    histos.fill(HIST("hTrackCount"), FillCutItsCluster - 1); // After ITS Clusters cut

    if (!TrkSels.cfgTrackCutsDCAz.value[UseGenTrkCuts].empty()) {
      if (std::fabs(track.dcaZ()) > FuncParas.fPtDepDCAz->Eval(track.pt())) {
        return false;
      }
    } else {
      if (std::fabs(track.dcaZ()) > trackCuts[TrkCutDCAz][UseGenTrkCuts]) {
        return false;
      }
    }
    histos.fill(HIST("hTrackCount"), FillCutDCAz - 1); // After DCAz cut

    if (!TrkSels.cfgTrackCutsDCAxy.value[UseGenTrkCuts].empty() && std::fabs(track.dcaXY()) > FuncParas.fPtDepDCAxy->Eval(track.pt())) {
      return false;
    }
    histos.fill(HIST("hTrackCount"), FillCutDCAxy - 1); // After DCAxy cut

    return true;
  }

  template <typename TTrack>
  bool trackSelectedForNch(TTrack const& track)
  {
    if (!TrkSels.cfgTrackCutsDCAxy.value[UseNchSelCuts].empty() && (std::fabs(track.dcaXY()) > FuncParas.fPtDepDCAxyForNch->Eval(track.pt()))) {
      return false;
    }
    if (!TrkSels.cfgTrackCutsDCAz.value[UseNchSelCuts].empty()) {
      if (std::fabs(track.dcaZ()) > FuncParas.fPtDepDCAzForNch->Eval(track.pt())) {
        return false;
      }
    } else {
      if (std::fabs(track.dcaZ()) > trackCuts[TrkCutDCAz][UseNchSelCuts]) {
        return false;
      }
    }
    return ((track.pt() > trackCuts[TrkCutPtMin][UseNchSelCuts]) && (track.pt() < trackCuts[TrkCutPtMax][UseNchSelCuts]) && (std::abs(track.eta()) < trackCuts[TrkCutEta][UseNchSelCuts]) && (track.tpcChi2NCl() < trackCuts[TrkCutChi2PrTpcCls][UseNchSelCuts]) && (track.tpcNClsFound() >= trackCuts[TrkCutTpcCluster][UseNchSelCuts]) && (track.tpcNClsCrossedRows() >= trackCuts[TrkCutTpcCrossedRows][UseNchSelCuts]) && (track.itsNCls() >= trackCuts[TrkCutItsCluster][UseNchSelCuts]));
  }

  template <typename TTrack>
  int getNsigmaPID(const TTrack& track)
  {
    if (PIDConfigs.cfgPIDUseNclsPID && track.tpcNClsPID() < PIDConfigs.cfgPIDNclsPIDCut) {
      return -1;
    }

    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {itsResponse.nSigmaITS<o2::track::PID::Pion>(track), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = -1; // -1 = not identified, 1 = pion, 2 = kaon, 3 = proton

    std::array<float, 3> nSigmaToUse = PIDConfigs.cfgPIDUseIts ? nSigmaITS : nSigmaTPC; // Choose which nSigma to use: TPC or ITS
    int useIndexDetector = UseTPC;
    useIndexDetector = PIDConfigs.cfgPIDUseIts ? UseITS : UseTPC; // Choose which nSigma to use: TPC or ITS

    bool isPion = false;
    bool isKaon = false;
    bool isProton = false;
    bool isDetectedPion = false;
    bool isDetectedKaon = false;
    bool isDetectedProton = false;
    bool isTofPion = false;
    bool isTofKaon = false;
    bool isTofProton = false;

    isDetectedPion = nSigmaToUse[IndexPionUp] < nSigmaVals[IndexPionUp][useIndexDetector] && nSigmaToUse[IndexPionUp] > nSigmaVals[IndexPionLow][useIndexDetector];
    isDetectedKaon = nSigmaToUse[IndexKaonUp] < nSigmaVals[IndexKaonUp][useIndexDetector] && nSigmaToUse[IndexKaonUp] > nSigmaVals[IndexKaonLow][useIndexDetector];
    isDetectedProton = nSigmaToUse[IndexProtonUp] < nSigmaVals[IndexProtonUp][useIndexDetector] && nSigmaToUse[IndexProtonUp] > nSigmaVals[IndexProtonLow][useIndexDetector];

    isTofPion = nSigmaTOF[IndexPionUp] < nSigmaVals[IndexPionUp][UseTOF] && nSigmaTOF[IndexPionUp] > nSigmaVals[IndexPionLow][UseTOF];
    isTofKaon = nSigmaTOF[IndexKaonUp] < nSigmaVals[IndexKaonUp][UseTOF] && nSigmaTOF[IndexKaonUp] > nSigmaVals[IndexKaonLow][UseTOF];
    isTofProton = nSigmaTOF[IndexProtonUp] < nSigmaVals[IndexProtonUp][UseTOF] && nSigmaTOF[IndexProtonUp] > nSigmaVals[IndexProtonLow][UseTOF];

    if (track.pt() > PIDConfigs.cfgPIDTofPtCut && !track.hasTOF()) {
      return -1;
    }
    if (track.pt() > PIDConfigs.cfgPIDTofPtCut && track.hasTOF()) {
      isPion = isTofPion && isDetectedPion;
      isKaon = isTofKaon && isDetectedKaon;
      isProton = isTofProton && isDetectedProton;
    } else {
      isPion = isDetectedPion;
      isKaon = isDetectedKaon;
      isProton = isDetectedProton;
    }

    if (PIDConfigs.cfgPIDUseRejection && ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton))) {
      return -1; // more than one particle satisfy the criteria
    }

    if (isPion) {
      pid = UsePions;
    } else if (isKaon) {
      pid = UseKaons;
    } else if (isProton) {
      pid = UseProtons;
    } else {
      return -1; // no particle satisfies the criteria
    }

    return pid; // -1 = not identified, 1 = pion, 2 = kaon, 3 = proton
  }

  void loadCorrections(uint64_t timestamp)
  {
    if (correctionsLoaded) {
      return;
    }
    // Loading efficiency histograms
    if (!Weights.cfgEfficiency.value.empty()) {
      mEfficiency = ccdb->getForTimeStamp<TH1D>(Weights.cfgEfficiency, timestamp);
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", Weights.cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", Weights.cfgEfficiency.value.c_str(), static_cast<void*>(mEfficiency));
    }
    // Loading acceptance histograms
    if (!Weights.cfgAcceptance.value.empty()) {
      mAcceptance.clear();
      mAcceptance.resize(NumCountTypes);

      // Charged particles
      mAcceptance[UseCharged] = ccdb->getForTimeStamp<GFWWeights>(Weights.cfgAcceptance.value + "_ch", timestamp);
      if (mAcceptance[UseCharged]) {
        LOGF(info, "Loaded acceptance weights from %s_ch (%p)", Weights.cfgAcceptance.value.c_str(), static_cast<void*>(mAcceptance[UseCharged]));
      } else {
        LOGF(fatal, "Could not load acceptance weights from %s_ch (%p)", Weights.cfgAcceptance.value.c_str(), static_cast<void*>(mAcceptance[UseCharged]));
      }
      // Pions
      mAcceptance[UsePions] = ccdb->getForTimeStamp<GFWWeights>(Weights.cfgAcceptance.value + "_pi", timestamp);
      if (mAcceptance[UsePions]) {
        LOGF(info, "Loaded acceptance weights from %s_pi (%p)", Weights.cfgAcceptance.value.c_str(), static_cast<void*>(mAcceptance[UsePions]));
      } else {
        LOGF(fatal, "Could not load acceptance weights from %s_pi (%p)", Weights.cfgAcceptance.value.c_str(), static_cast<void*>(mAcceptance[UsePions]));
      }
      // Kaons
      mAcceptance[UseKaons] = ccdb->getForTimeStamp<GFWWeights>(Weights.cfgAcceptance.value + "_ka", timestamp);
      if (mAcceptance[UseKaons]) {
        LOGF(info, "Loaded acceptance weights from %s_ka (%p)", Weights.cfgAcceptance.value.c_str(), static_cast<void*>(mAcceptance[UseKaons]));
      } else {
        LOGF(fatal, "Could not load acceptance weights from %s_ka (%p)", Weights.cfgAcceptance.value.c_str(), static_cast<void*>(mAcceptance[UseKaons]));
      }
      // Protons
      mAcceptance[UseProtons] = ccdb->getForTimeStamp<GFWWeights>(Weights.cfgAcceptance.value + "_pr", timestamp);
      if (mAcceptance[UseProtons]) {
        LOGF(info, "Loaded acceptance weights from %s_pr (%p)", Weights.cfgAcceptance.value.c_str(), static_cast<void*>(mAcceptance[UseProtons]));
      } else {
        LOGF(fatal, "Could not load acceptance weights from %s_pr (%p)", Weights.cfgAcceptance.value.c_str(), static_cast<void*>(mAcceptance[UseProtons]));
      }
      // Reference charged particles
      mAcceptance[UseRef] = ccdb->getForTimeStamp<GFWWeights>(Weights.cfgAcceptance.value + "_ref", timestamp);
      if (mAcceptance[UseRef]) {
        LOGF(info, "Loaded acceptance weights from %s_ref (%p)", Weights.cfgAcceptance.value.c_str(), static_cast<void*>(mAcceptance[UseRef]));
      } else {
        LOGF(fatal, "Could not load acceptance weights from %s_ref (%p)", Weights.cfgAcceptance.value.c_str(), static_cast<void*>(mAcceptance[UseRef]));
      }
    }
    correctionsLoaded = true;
  }

  bool getEfficiencyCorrection(float& weight_nue, float pt)
  {
    float effNch = 1.;
    if (mEfficiency) {
      int ptBin = mEfficiency->FindBin(pt);
      effNch = mEfficiency->GetBinContent(ptBin);
    } else {
      effNch = 1.0;
    }
    if (effNch == 0) {
      return false;
    }
    weight_nue = 1. / effNch;
    return true;
  }

  template <typename TTracks>
  void nchCounter(const TTracks& tracks, double& multiplicity) // function to count the number of tracks in the event and fill the histogram
  {
    double nTracksCorrected = 0;
    float weight = 1.0f;
    for (auto const& track : tracks) {
      if (!trackSelectedForNch(track)) {
        continue;
      }

      if (!getEfficiencyCorrection(weight, track.pt())) {
        continue;
      }
      nTracksCorrected += weight;
    }
    multiplicity = nTracksCorrected;
  }

  template <typename TTrack>
  void fillTrackQA(TTrack const& track, int pid, double eventClass, Stage stage, float weff, float wacc) // function to fill the QA after Nsigma selection
  {
    if (stage == Before && pid == UseCharged) {
      histos.fill(HIST("hPhi"), RecoDecay::constrainAngle(track.phi(), 0.0), UseCharged);
      histos.fill(HIST("hPhiWeighted"), RecoDecay::constrainAngle(track.phi(), 0.0), UseCharged, wacc);
      histos.fill(HIST("hEta"), track.eta(), UseCharged);
      histos.fill(HIST("hEtaCorrected"), track.eta(), UseCharged, weff);
      histos.fill(HIST("hPt"), track.pt(), UseCharged, eventClass);
      histos.fill(HIST("hPtCorrected"), track.pt(), UseCharged, eventClass, weff);

    } else if (stage == After && pid > UseCharged) {
      histos.fill(HIST("hPhi"), RecoDecay::constrainAngle(track.phi(), 0.0), pid);
      histos.fill(HIST("hPhiWeighted"), RecoDecay::constrainAngle(track.phi(), 0.0), pid, wacc);
      histos.fill(HIST("hEta"), track.eta(), pid);
      histos.fill(HIST("hEtaCorrected"), track.eta(), pid, weff);
      histos.fill(HIST("hPt"), track.pt(), pid, eventClass);
      histos.fill(HIST("hPtCorrected"), track.pt(), pid, eventClass, weff);
    }
  }

  template <typename TTrack>
  void fillPIDQA(TTrack const& track, int pid, double eventClass, Stage stage)
  {
    const bool useITS = PIDConfigs.cfgPIDUseIts;
    double detectorNSigma = 0.0;
    double tpcExpSigma = 0.0;
    double tofNSigma = 0.0;
    double tpcNSigma = 0.0;
    switch (pid) {
      case UsePions: // For Pions
        tofNSigma = track.tofNSigmaPi();
        tpcExpSigma = track.tpcExpSigmaPi();
        tpcNSigma = track.tpcNSigmaPi();
        detectorNSigma = useITS ? itsResponse.nSigmaITS<o2::track::PID::Pion>(track) : tpcNSigma;
        break;
      case UseKaons: // For Kaons
        tofNSigma = track.tofNSigmaKa();
        tpcExpSigma = track.tpcExpSigmaKa();
        tpcNSigma = track.tpcNSigmaKa();
        detectorNSigma = useITS ? itsResponse.nSigmaITS<o2::track::PID::Kaon>(track) : tpcNSigma;
        break;
      case UseProtons: // For Protons
        tofNSigma = track.tofNSigmaPr();
        tpcExpSigma = track.tpcExpSigmaPr();
        tpcNSigma = track.tpcNSigmaPr();
        detectorNSigma = useITS ? itsResponse.nSigmaITS<o2::track::PID::Proton>(track) : tpcNSigma;
        break;
      default:
        return;
    }

    if (PIDConfigs.cfgQANsigma) {
      if (stage == Before) { // This is filled before deciding the particle type, so this segment is filled four times
        if (!useITS) {
          histos.fill(HIST("TofTpcNsigma_beforeCut"), detectorNSigma, tofNSigma, eventClass, track.pt(), pid);
        } else {
          histos.fill(HIST("TofItsNsigma_beforeCut"), detectorNSigma, tofNSigma, eventClass, track.pt(), pid);
        }
      } else if (stage == After) { // This is filled after deciding the particle type, so this segment is filled once
        if (!useITS) {
          histos.fill(HIST("TofTpcNsigma_afterCut"), detectorNSigma, tofNSigma, eventClass, track.pt(), pid);
        } else {
          histos.fill(HIST("TofItsNsigma_afterCut"), detectorNSigma, tofNSigma, eventClass, track.pt(), pid);
        }
      }
    } // If switched on cfgQANsigma

    if (PIDConfigs.cfgQAdEdx) {
      double tpcExpSignal = track.tpcSignal() - (tpcNSigma * tpcExpSigma);

      if (stage == Before) {
        histos.fill(HIST("TpcdEdx_beforeCut"), track.tpcSignal(), tofNSigma, eventClass, track.pt(), pid);
        histos.fill(HIST("ExpTpcdEdx_beforeCut"), tpcExpSignal, tofNSigma, eventClass, track.pt(), pid);
        histos.fill(HIST("ExpSigma_beforeCut"), tpcExpSigma, tofNSigma, eventClass, track.pt(), pid);
      } else if (stage == After) {
        histos.fill(HIST("TpcdEdx_afterCut"), track.tpcSignal(), tofNSigma, eventClass, track.pt(), pid);
        histos.fill(HIST("ExpTpcdEdx_afterCut"), tpcExpSignal, tofNSigma, eventClass, track.pt(), pid);
        histos.fill(HIST("ExpSigma_afterCut"), tpcExpSigma, tofNSigma, eventClass, track.pt(), pid);
      }
    } // If switched on cfgQAdEdx

    if (PIDConfigs.cfgQABetadEdx) {
      double tpcExpSignal = track.tpcSignal() - (tpcNSigma * tpcExpSigma);

      if (stage == Before) {
        histos.fill(HIST("TpcdEdx_TofBeta_beforeCut"), track.tpcSignal(), track.beta(), eventClass, track.pt(), pid);
        histos.fill(HIST("ExpTpcdEdx_TofBeta_beforeCut"), tpcExpSignal, track.beta(), eventClass, track.pt(), pid);
        histos.fill(HIST("ExpSigma_TofBeta_beforeCut"), tpcExpSigma, track.beta(), eventClass, track.pt(), pid);
      } else if (stage == After) {
        histos.fill(HIST("TpcdEdx_TofBeta_afterCut"), track.tpcSignal(), track.beta(), eventClass, track.pt(), pid);
        histos.fill(HIST("ExpTpcdEdx_TofBeta_afterCut"), tpcExpSignal, track.beta(), eventClass, track.pt(), pid);
        histos.fill(HIST("ExpSigma_TofBeta_afterCut"), tpcExpSigma, track.beta(), eventClass, track.pt(), pid);
      }
    } // If switched on cfgQABetadEdx
  }

  void fillFC(const GFW::CorrConfig& corrconf, const double& eventClass, const double& rndm)
  {
    auto fillSingleContainer = [&](const GFW::CorrConfig& cc, const double& ec, const int& ib, const double& rd, const std::string& cn) {
      auto dnx = fGFW->Calculate(cc, ib, kTRUE).real();
      if (dnx == 0) {
        return;
      }
      auto val = fGFW->Calculate(cc, ib, kFALSE).real() / dnx;
      if (std::abs(val) < 1) {
        fFC->FillProfile(cn.c_str(), ec, val, dnx, rd);
      }
    };

    // For pT-integrated
    if (!corrconf.pTDif) {
      fillSingleContainer(corrconf, eventClass, 0, rndm, corrconf.Head);
    } else { // For pT-differential
      for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
        std::string name = corrconf.Head + "_pt_" + std::to_string(i);
        fillSingleContainer(corrconf, eventClass, i, rndm, name);
      }
    }
  }

  void createRunByRunHistos(int runNumber)
  {
    if (Weights.cfgOutputNUAWeights) {
      std::vector<std::shared_ptr<TH3>> tH3s(NumCountTypes);
      tH3s[UseCharged] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_ch", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      tH3s[UsePions] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_pi", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      tH3s[UseKaons] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_ka", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      tH3s[UseProtons] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_pr", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      tH3s[UseRef] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_ref", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      th3sList.insert(std::make_pair(runNumber, tH3s));
    }
  }

  template <typename TTrack>
  double getAcceptance(const TTrack& track, float cent, float vtxz, int index)
  { // ch = 0, pi = 1, ka = 2, pr = 3, ref = 4
    if (index < 0 || index >= NumCountTypes) {
      return 1;
    }
    double wacc = 1;

    if ((Weights.cfgUseWeightPhiEtaVtxz && Weights.cfgUseWeightPhiPtCent) || (Weights.cfgUseWeightPhiEtaPt && Weights.cfgUseWeightPhiPtCent) || (Weights.cfgUseWeightPhiEtaVtxz && Weights.cfgUseWeightPhiEtaPt)) {
      LOGF(fatal, "Only one of the three weight options can be used at a time");
    }

    if (!mAcceptance.empty() && correctionsLoaded) {
      if (!mAcceptance[index]) {
        LOGF(fatal, "Acceptance weights not loaded for index %d", index);
        return 1;
      }
      if (Weights.cfgUseWeightPhiEtaVtxz) {
        wacc = mAcceptance[index]->getNUA(track.phi(), track.eta(), vtxz);
      }
      if (Weights.cfgUseWeightPhiPtCent) {
        wacc = mAcceptance[index]->getNUA(track.phi(), track.pt(), cent);
      }
      if (Weights.cfgUseWeightPhiEtaPt) {
        wacc = mAcceptance[index]->getNUA(track.phi(), track.eta(), track.pt());
      }
    }
    return wacc;
  }

  template <typename TTrack>
  void fillWeights(const TTrack& track, float cent, float vtxz, const int& pid_index, const int& run)
  {
    double pt = track.pt();
    bool withinPtPOI = (GeneralConfigs.cfgCutPtPOIMin < pt) && (pt < GeneralConfigs.cfgCutPtPOIMax); // within POI pT range
    bool withinPtRef = (GeneralConfigs.cfgRefPtMin < pt) && (pt < GeneralConfigs.cfgRefPtMax);       // within RF pT range

    if (Weights.cfgOutputRunByRun) {
      if (withinPtRef && (pid_index < UseCharged)) {
        th3sList[run][UseRef]->Fill(track.phi(), track.eta(), vtxz); // pt-subset of charged particles for ref flow
      }
      if (withinPtPOI) {
        th3sList[run][UseCharged]->Fill(track.phi(), track.eta(), vtxz); // charged particle weights
        if (pid_index > UseCharged) {
          th3sList[run][pid_index]->Fill(track.phi(), track.eta(), vtxz); // id'ed particle weights
        }
      }
    } else { // Fill NUA histograms for the entire run
      // Fill NUA histograms for reference particles only when they are not identified as PID particle
      // This will only be used for pT-dependent flow calculation for PID particles
      if (withinPtRef && (pid_index < UseCharged)) {
        histos.fill(HIST("NUA/hPhiEtaVtxz_ref"), track.phi(), track.eta(), vtxz); // pt-subset of charged particles for ref flow
        histos.fill(HIST("NUA/hPhiPtCent_ref"), track.phi(), track.pt(), cent);
        histos.fill(HIST("NUA/hPhiEtaPt_ref"), track.phi(), track.eta(), track.pt());
      }

      if (withinPtPOI) {
        histos.fill(HIST("NUA/hPhiEtaVtxz_ch"), track.phi(), track.eta(), vtxz); // charged particle weights
        histos.fill(HIST("NUA/hPhiPtCent_ch"), track.phi(), track.pt(), cent);
        histos.fill(HIST("NUA/hPhiEtaPt_ch"), track.phi(), track.eta(), track.pt());

        switch (pid_index) {
          case UsePions:
            histos.fill(HIST("NUA/hPhiEtaVtxz_pi"), track.phi(), track.eta(), vtxz); // pion weights
            histos.fill(HIST("NUA/hPhiPtCent_pi"), track.phi(), track.pt(), cent);
            histos.fill(HIST("NUA/hPhiEtaPt_pi"), track.phi(), track.eta(), track.pt());
            break;
          case UseKaons:
            histos.fill(HIST("NUA/hPhiEtaVtxz_ka"), track.phi(), track.eta(), vtxz); // kaon weights
            histos.fill(HIST("NUA/hPhiPtCent_ka"), track.phi(), track.pt(), cent);
            histos.fill(HIST("NUA/hPhiEtaPt_ka"), track.phi(), track.eta(), track.pt());
            break;
          case UseProtons:
            histos.fill(HIST("NUA/hPhiEtaVtxz_pr"), track.phi(), track.eta(), vtxz); // proton weights
            histos.fill(HIST("NUA/hPhiPtCent_pr"), track.phi(), track.pt(), cent);
            histos.fill(HIST("NUA/hPhiEtaPt_pr"), track.phi(), track.eta(), track.pt());
            break;
          default:
            break;
        }
      }
    }
  }

  template <typename TCollision>
  bool selectionEvent(TCollision const& collision, const int mult, const double cent, const bool fillCounter)
  {
    if (fillCounter) {
      histos.fill(HIST("hEventCount"), FilteredEvents);
    }
    if (!collision.sel8()) {
      return false;
    }
    if (fillCounter) {
      histos.fill(HIST("hEventCount"), AfterSel8);
    }

    if (eventCuts[UseNoTimeFrameBorder][EvCut1] && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillCounter && eventCuts[UseNoTimeFrameBorder][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseNoTimeFrameBorder);
    }

    if (eventCuts[UseNoITSROFrameBorder][EvCut1] && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillCounter && eventCuts[UseNoITSROFrameBorder][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseNoITSROFrameBorder);
    }

    if (eventCuts[UseNoSameBunchPileup][EvCut1] && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return false;
    }
    if (fillCounter && eventCuts[UseNoSameBunchPileup][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseNoSameBunchPileup);
    }

    if (eventCuts[UseGoodZvtxFT0vsPV][EvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return false;
    }
    if (fillCounter && eventCuts[UseGoodZvtxFT0vsPV][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseGoodZvtxFT0vsPV);
    }

    if (eventCuts[UseNoCollInTimeRangeStandard][EvCut1] && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return false;
    }

    if (fillCounter && eventCuts[UseNoCollInTimeRangeStandard][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseNoCollInTimeRangeStandard);
    }

    if (eventCuts[UseGoodITSLayersAll][EvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return false;
    }

    if (fillCounter && eventCuts[UseGoodITSLayersAll][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseGoodITSLayersAll);
    }

    if (eventCuts[UseGoodITSLayer0123][EvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      return false;
    }
    if (fillCounter && eventCuts[UseGoodITSLayer0123][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseGoodITSLayer0123);
    }

    if (eventCuts[UseNoCollInRofStandard][EvCut1] && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return false;
    }

    if (fillCounter && eventCuts[UseNoCollInRofStandard][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseNoCollInRofStandard);
    }

    if (eventCuts[UseNoHighMultCollInPrevRof][EvCut1] && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return false;
    }
    if (fillCounter && eventCuts[UseNoHighMultCollInPrevRof][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseNoHighMultCollInPrevRof);
    }

    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();

    if (eventCuts[UseOccupancy][EvCut1] && (occupancy < EvSels.cfgCutOccupancyMin || occupancy > EvSels.cfgCutOccupancyMax)) {
      return false;
    }
    if (fillCounter && eventCuts[UseOccupancy][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseOccupancy);
    }

    if (eventCuts[UseMultCorrCut][EvCut1]) {
      if (FuncParas.cfgMultPVT0CCutEnabled) {
        if (multNTracksPV < FuncParas.fMultPVT0CCutLow->Eval(cent)) {
          return false;
        }
        if (multNTracksPV > FuncParas.fMultPVT0CCutHigh->Eval(cent)) {
          return false;
        }
      }
      if (FuncParas.cfgMultT0CCutEnabled) {
        if (mult < FuncParas.fMultT0CCutLow->Eval(cent)) {
          return false;
        }
        if (mult > FuncParas.fMultT0CCutHigh->Eval(cent)) {
          return false;
        }
      }
      if (FuncParas.cfgMultGlobalPVCutEnabled) {
        if (mult < FuncParas.fMultGlobalPVCutLow->Eval(multNTracksPV)) {
          return false;
        }
        if (mult > FuncParas.fMultGlobalPVCutHigh->Eval(multNTracksPV)) {
          return false;
        }
      }
      if (FuncParas.cfgMultMultV0ACutEnabled) {
        if (collision.multFV0A() < FuncParas.fMultMultV0ACutLow->Eval(mult)) {
          return false;
        }
        if (collision.multFV0A() > FuncParas.fMultMultV0ACutHigh->Eval(mult)) {
          return false;
        }
      }
    }

    if (fillCounter && eventCuts[UseMultCorrCut][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseMultCorrCut);
    }

    // V0A T0A 5 sigma cut
    if (eventCuts[UseT0AV0ACut][EvCut1] && (std::fabs(collision.multFV0A() - FuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > FuncParas.cfgV0AT0Acut * FuncParas.fT0AV0ASigma->Eval(collision.multFT0A()))) {
      return false;
    }
    if (fillCounter && eventCuts[UseT0AV0ACut][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseT0AV0ACut);
    }

    return true;
  }

  void processGFW(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracksWithoutBayes const& tracks)
  {
    if (tracks.size() < 1) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float cent = EvSels.cfgCentTableUnavailable ? -1.f : getCentrality(collision);

    if (!selectionEvent(collision, tracks.size(), cent, true)) {
      return;
    }

    int runNumber = bc.runNumber();
    if (Weights.cfgOutputRunByRun && runNumber != lastRunNumer) {
      lastRunNumer = runNumber;
      if (std::find(runNumbers.begin(), runNumbers.end(), runNumber) == runNumbers.end()) {
        // if run number is not in the preconfigured list, create new output histograms for this run
        createRunByRunHistos(runNumber);
        runNumbers.push_back(runNumber);
      }
    }

    fGFW->Clear();

    loadCorrections(bc.timestamp()); // load corrections for each event
    if (!EvSels.cfgCentTableUnavailable) {
      histos.fill(HIST("hCent"), cent);
    }
    histos.fill(HIST("hNch"), tracks.size());
    histos.fill(HIST("hVtxZ"), collision.posZ());

    auto multiplicity = static_cast<double>(tracks.size());
    if (GeneralConfigs.cfgStrictTrackCounter) {
      nchCounter(tracks, multiplicity);
    }
    histos.fill(HIST("hNchCorrected"), multiplicity);

    // Track loop for calculating the Qn angles
    double psi2Est = 0, psi3Est = 0, psi4Est = 0;
    float wEPeff = 1;
    double v2 = 0, v3 = 0, v4 = 0;
    // be cautious, this only works for Pb-Pb
    // esimate the Qn angles and vn for this event
    if (TrkDensCorr.cfgTrackDensityCorrUse) {
      double q2x = 0, q2y = 0;
      double q3x = 0, q3y = 0;
      double q4x = 0, q4y = 0;
      for (const auto& track : tracks) {
        bool withinPtRef = (GeneralConfigs.cfgRefPtMin < track.pt()) && (track.pt() < GeneralConfigs.cfgRefPtMax); // within RF pT range
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

      v2 = funcV2->Eval(cent);
      v3 = funcV3->Eval(cent);
      v4 = funcV4->Eval(cent);
    }

    double fSampleIndex = fRandom->Rndm();
    auto eventClass = (EvSels.cfgUseNch) ? multiplicity : cent;

    // Actual track loop
    for (auto const& track : tracks) {
      if (!trackSelected(track)) {
        continue;
      }

      // Defining the pT regions for POI and REF
      bool withinPtPOI = (GeneralConfigs.cfgCutPtPOIMin < track.pt()) && (track.pt() < GeneralConfigs.cfgCutPtPOIMax); // within POI pT range
      bool withinPtRef = (GeneralConfigs.cfgRefPtMin < track.pt()) && (track.pt() < GeneralConfigs.cfgRefPtMax);       // within RF pT range

      if (!withinPtPOI && !withinPtRef) {
        continue;
      }

      // Initializing weff for each track
      float weff = 1.0;
      float weffPIDPOI = 1.0;
      if (!Weights.cfgEfficiency.value.empty()) {
        getEfficiencyCorrection(weff, track.pt());
      }

      // Get acceptance weights for reference and charged particles
      double waccRef = (!mAcceptance.empty()) ? getAcceptance(track, cent, collision.posZ(), UseRef) : 1.0;
      double waccChargedPOI = (!mAcceptance.empty()) ? getAcceptance(track, cent, collision.posZ(), UseCharged) : 1.0;

      // Filling QA histograms
      if (GeneralConfigs.cfgQABasic) {
        fillTrackQA(track, UseCharged, eventClass, Before, weff, waccChargedPOI); // Fills track properties for all charged particles
      }

      if (PIDConfigs.cfgQANsigma) {
        fillPIDQA(track, UsePions, eventClass, Before);   // Fills Nsigma histograms for pions before PID selection
        fillPIDQA(track, UseKaons, eventClass, Before);   // Fills Nsigma histograms for kaons before PID selection
        fillPIDQA(track, UseProtons, eventClass, Before); // Fills Nsigma histograms for protons before PID selection
      }

      // Select track as PID particle
      int pidIndex = getNsigmaPID(track);

      double waccPIDPOI = (!mAcceptance.empty()) ? getAcceptance(track, cent, collision.posZ(), pidIndex) : 1.0;
      if (withinPtRef && withinPtPOI && (pidIndex > UseCharged)) {
        waccRef = waccPIDPOI; // if particle is both (then it's overlap), override ref with POI
      }

      // Filling QA histograms for PID particles after PID selection
      if (pidIndex > UseCharged) {
        if (GeneralConfigs.cfgQABasic) {
          fillTrackQA(track, pidIndex, eventClass, After, weffPIDPOI, waccPIDPOI); // Fills track properties and Nsigma histograms
        }
        if (PIDConfigs.cfgQANsigma) {
          fillPIDQA(track, pidIndex, eventClass, After);
        }
      }

      // Fill histograms for NUA weights
      if (Weights.cfgOutputNUAWeights) {
        fillWeights(track, cent, collision.posZ(), pidIndex, runNumber);
      }

      // Track density correction
      if (TrkDensCorr.cfgTrackDensityCorrUse && withinPtRef) {
        double fphi = v2 * std::cos(2 * (track.phi() - psi2Est)) + v3 * std::cos(3 * (track.phi() - psi3Est)) + v4 * std::cos(4 * (track.phi() - psi4Est));
        fphi = (1 + 2 * fphi);
        int pTBinForEff = hFindPtBin->FindBin(track.pt());
        if (pTBinForEff >= 1 && pTBinForEff <= hFindPtBin->GetNbinsX()) {
          wEPeff = funcEff[pTBinForEff - 1]->Eval(fphi * tracks.size());
          if (wEPeff > 0.) {
            wEPeff = 1. / wEPeff;
            weff *= wEPeff;
          }
        }
      } // end of track density correction loop

      auto ptBin = fPtAxis->FindBin(track.pt()) - 1;

      // Filling GFW
      if (withinPtRef) {
        fGFW->Fill(track.eta(), ptBin, track.phi(), waccRef * weff, 1);
      }
      if (withinPtPOI) {
        fGFW->Fill(track.eta(), ptBin, track.phi(), waccChargedPOI * weff, 128);
      }
      if (withinPtPOI && withinPtRef) {
        fGFW->Fill(track.eta(), ptBin, track.phi(), waccChargedPOI * weff, 256);
      }

      if (pidIndex > UseCharged) {
        if (withinPtPOI) {
          fGFW->Fill(track.eta(), ptBin, track.phi(), waccPIDPOI * weff, 1 << (pidIndex));
        }
        if (withinPtPOI && withinPtRef) {
          fGFW->Fill(track.eta(), ptBin, track.phi(), waccPIDPOI * weff, 1 << (pidIndex + 3));
        }
      }
    } // track loop ends

    // Filling FlowContainer
    for (auto const& corrconfig : corrConfigs) {
      fillFC(corrconfig, eventClass, fSampleIndex);
    }
  } // end of process
  PROCESS_SWITCH(FlowPbpbPikp, processGFW, "Process for GFW flow analysis", true);

  void processPIDQA(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracksWithoutBayes const& tracks)
  {
    if (tracks.size() < 1) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float cent = EvSels.cfgCentTableUnavailable ? -1.f : getCentrality(collision);

    if (!selectionEvent(collision, tracks.size(), cent, true)) {
      return;
    }

    loadCorrections(bc.timestamp()); // load corrections for each event

    auto multiplicity = static_cast<double>(tracks.size());
    if (GeneralConfigs.cfgStrictTrackCounter) {
      nchCounter(tracks, multiplicity);
    }

    auto eventClass = (EvSels.cfgUseNch) ? multiplicity : cent;

    for (auto const& track : tracks) {
      if (!trackSelected(track)) {
        continue;
      }

      fillPIDQA(track, UsePions, eventClass, Before);   // Fills Nsigma histograms for pions before PID selection
      fillPIDQA(track, UseKaons, eventClass, Before);   // Fills Nsigma histograms for kaons before PID selection
      fillPIDQA(track, UseProtons, eventClass, Before); // Fills Nsigma histograms for protons before PID selection

      // Select track as PID particle
      int pidIndex = getNsigmaPID(track);

      if (pidIndex > UseCharged) {
        fillPIDQA(track, pidIndex, eventClass, After);
      }

    } // track loop ends
  } // end of processPIDQA

  PROCESS_SWITCH(FlowPbpbPikp, processPIDQA, "Process for PID related QA production", true);

  void processMCPID(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracksMcLabels const& tracksmc, aod::McParticles const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float cent = EvSels.cfgCentTableUnavailable ? -1.f : getCentrality(collision);

    if (!selectionEvent(collision, tracksmc.size(), cent, true)) {
      return;
    }

    loadCorrections(bc.timestamp()); // load corrections for each event

    auto multiplicity = static_cast<double>(tracksmc.size());
    if (GeneralConfigs.cfgStrictTrackCounter) {
      nchCounter(tracksmc, multiplicity);
    }

    auto eventClass = (EvSels.cfgUseNch) ? multiplicity : cent;

    // loop over all tracks
    for (auto const& track : tracksmc) {

      if (!trackSelected(track)) {
        continue;
      }

      int pdgIndex = -1;
      if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kPiPlus) {
        pdgIndex = UsePions;
      } else if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kKPlus) {
        pdgIndex = UseKaons;
      } else if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kProton) {
        pdgIndex = UseProtons;
      } else {
        pdgIndex = -1;
      }

      // Fill Counts for MC true particles
      if (pdgIndex > UseCharged) {
        histos.fill(HIST("hNAllTrues"), eventClass, track.pt(), pdgIndex);
      }

      int pidIndex = getNsigmaPID(track);
      // Fill Counts for selection through Nsigma cuts
      if (pidIndex > UseCharged) {
        histos.fill(HIST("hNAllPositives"), eventClass, track.pt(), pidIndex);

        // Fill Counts for selected particles if it matches their PDG particle
        if (pidIndex == pdgIndex) {
          histos.fill(HIST("hNTruePositives"), eventClass, track.pt(), pidIndex);
        }
      }

      if (track.mcParticle().isPhysicalPrimary()) {
        histos.fill(HIST("TpcdEdx_AllPositives_physprim"), track.tpcSignal(), track.beta(), track.pt(), pidIndex);
        histos.fill(HIST("TpcdEdx_AllTrues_physprim"), track.tpcSignal(), track.beta(), track.pt(), pdgIndex);
      }

    } // end of tracks MC loop
  }
  PROCESS_SWITCH(FlowPbpbPikp, processMCPID, "Process for MC PID QA production", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FlowPbpbPikp>(cfgc)};
}
