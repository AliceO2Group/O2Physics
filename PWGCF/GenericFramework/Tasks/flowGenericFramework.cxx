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

#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/FlowPtContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWConfig.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

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
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <ios>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::analysis::genericframework;

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
GFWCorrConfigs configsV02;
GFWCorrConfigs configsV0;
std::vector<std::pair<double, double>> etagapsPtPt;
std::vector<double> multGlobalCorrCutPars;
std::vector<double> multPVCorrCutPars;
std::vector<double> multGlobalPVCorrCutPars;
std::vector<double> multGlobalV0ACutPars;
std::vector<double> multGlobalT0ACutPars;
} // namespace o2::analysis::gfw

template <typename T>
auto projectMatrix(Array2D<T> const& mat, std::array<float, 6>& array1, std::array<float, 6>& array2, std::array<float, 6>& array3)
{
  for (auto j = 0; j < static_cast<int>(mat.cols); ++j) {
    array1[j] = mat(0, j);
    array2[j] = mat(1, j);
    array3[j] = mat(2, j);
  }
  return;
}
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

static constexpr float LongArrayFloat[3][20] = {{1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {2.1, 2.2, 2.3, -2.1, -2.2, -2.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {3.1, 3.2, 3.3, -3.1, -3.2, -3.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}};
static constexpr int LongArrayInt[3][20] = {{1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1}, {2, 2, 2, -2, -2, -2, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1}, {3, 3, 3, -3, -3, -3, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1}};
static constexpr double LongArrayDouble[4][2] = {{-0.8, -0.5}, {0.5, 0.8}, {-2, -2}, {-2, -2}};

struct FlowGenericFramework {
  O2_DEFINE_CONFIGURABLE(cfgUseNewNpt, bool, true, "Use alternative filling for fractions")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgMpar, int, 8, "Highest order of pt-pt correlations")
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A, 4:NTPV, 5:NGlobal, 6:MFT")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Do correlations as function of Nch")
  O2_DEFINE_CONFIGURABLE(cfgUseNchCorrection, int, 1, "Use correction for Nch; 0: Use size of tracks table, 1: Use efficiency-corrected Nch values, 2: Use uncorrected Nch values");
  O2_DEFINE_CONFIGURABLE(cfgFillWeights, bool, false, "Fill NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgRunByRun, bool, false, "Fill histograms on a run-by-run basis")
  O2_DEFINE_CONFIGURABLE(cfgFillQA, bool, false, "Fill QA histograms")
  O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, false, "Use additional event cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(cfgUseCentralMoments, bool, true, "Use central moments in vn-pt calculations")
  O2_DEFINE_CONFIGURABLE(cfgUsePID, bool, true, "Enable PID information")
  O2_DEFINE_CONFIGURABLE(cfgUseGapMethod, bool, false, "Use gap method in vn-pt calculations")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgUsePIDEfficiencies, bool, false, "Use species dependent efficiencies")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgPtmin, float, 0.2, "minimum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgPtmax, float, 10, "maximum pt (GeV/c)");
  Configurable<std::pair<float, float>> cfgEta{"cfgEta", {-0.8, 0.8}, "eta cut"};
  Configurable<std::pair<float, float>> cfgEtaNch{"cfgEtaNch", {-0.5, 0.5}, "eta cut for nch selection"};
  Configurable<std::pair<float, float>> cfgEtaPtPt{"cfgEtaPtPt", {-0.5, 0.5}, "eta for pt-pt correlation"};
  Configurable<std::pair<float, float>> cfgEtaV0Daughters{"cfgEtaV0Daughters", {-0.5, 0.5}, "eta cut on V0 daughter particles"};
  Configurable<LabeledArray<double>> cfgPtPtGaps{"cfgPtPtGaps", {LongArrayDouble[0], 4, 2, {"subevent 1", "subevent 2", "subevent 3", "subevent 4"}, {"etamin", "etamax"}}, "{etamin,etamax} for all ptpt-subevents"};
  O2_DEFINE_CONFIGURABLE(cfgUsePIDTotal, bool, false, "use fraction of PID total");
  O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10, "vertex cut (cm)");
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgDCAxyNSigma, float, 7, "Cut on number of sigma deviations from expected DCA in the transverse direction");
    O2_DEFINE_CONFIGURABLE(cfgDCAz, float, 2, "Cut on DCA in the longitudinal direction (cm)");
    O2_DEFINE_CONFIGURABLE(cfgNTPCCls, float, 50, "Cut on number of TPC clusters found");
    O2_DEFINE_CONFIGURABLE(cfgNTPCXrows, float, 70, "Cut on number of TPC crossed rows");
    O2_DEFINE_CONFIGURABLE(cfgMinNITSCls, float, 5, "Cut on minimum number of ITS clusters found");
    O2_DEFINE_CONFIGURABLE(cfgChi2PrITSCls, float, 36, "Cut on chi^2 per ITS clusters found");
    O2_DEFINE_CONFIGURABLE(cfgChi2PrTPCCls, float, 2.5, "Cut on chi^2 per TPC clusters found");
    O2_DEFINE_CONFIGURABLE(cfgTPCSectorCut, bool, false, "Cut on pt-phi distribution");
  } cfgTrackCuts;
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
  O2_DEFINE_CONFIGURABLE(cfgOccupancySelection, int, 2000, "Max occupancy selection, -999 to disable");
  O2_DEFINE_CONFIGURABLE(cfgDoOccupancySel, bool, true, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field; default CCDB will be queried");
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
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalASideCorrCutFunction, std::string, "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + [10]*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", "Functional for global vs V0A multiplicity low correlation cut");
    Configurable<std::vector<double>> cfgMultGlobalV0ACutPars{"cfgMultGlobalV0ACutPars", std::vector<double>{567.785, 172.715, 0.77888, -0.00693466, 1.40564e-05, 679.853, 66.8068, -0.444332, 0.00115002, -4.92064e-07}, "Global vs FV0A multiplicity cut parameter values"};
    Configurable<std::vector<double>> cfgMultGlobalT0ACutPars{"cfgMultGlobalT0ACutPars", std::vector<double>{241.618, 61.8402, 0.348049, -0.00306078, 6.20357e-06, 315.235, 29.1491, -0.188639, 0.00044528, -9.08912e-08}, "Global vs FT0A multiplicity cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgGlobalV0ALowSigma, float, -3, "Number of sigma deviations below expected value in global vs V0A correlation");
    O2_DEFINE_CONFIGURABLE(cfgGlobalV0AHighSigma, float, 4, "Number of sigma deviations above expected value in global vs V0A correlation");
    O2_DEFINE_CONFIGURABLE(cfgGlobalT0ALowSigma, float, -3., "Number of sigma deviations below expected value in global vs T0A correlation");
    O2_DEFINE_CONFIGURABLE(cfgGlobalT0AHighSigma, float, 4, "Number of sigma deviations above expected value in global vs T0A correlation");
  } cfgMultCorrCuts;
  struct : ConfigurableGroup {
    Configurable<LabeledArray<float>> nSigmas{"nSigmas", {LongArrayFloat[0], 3, 6, {"TPC", "TOF", "ITS"}, {"pos_pi", "pos_ka", "pos_pr", "neg_pi", "neg_ka", "neg_pr"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};
    Configurable<LabeledArray<float>> resonanceCuts{"resonanceCuts", {LongArrayFloat[0], 3, 11, {"K0", "Lambda", "Phi"}, {"cos_PAs", "massMin", "massMax", "PosTrackPt", "NegTrackPt", "DCAPosToPVMin", "DCANegToPVMin", "Lifetime", "RadiusMin", "RadiusMax", "Rapidity"}}, "Labeled array (float) for various cuts on resonances"};
    Configurable<LabeledArray<int>> resonanceSwitches{"resonanceSwitches", {LongArrayInt[0], 3, 6, {"K0", "Lambda", "Phi"}, {"UseParticle", "UseCosPA", "NMassBins", "DCABetDaug", "UseProperLifetime", "UseV0Radius"}}, "Labeled array (int) for various cuts on resonances"};
    O2_DEFINE_CONFIGURABLE(cfgUseLsPhi, bool, true, "Use LikeSign for Phi v2")
    O2_DEFINE_CONFIGURABLE(cfgUseOnlyTPC, bool, true, "Use only TPC PID for daughter selection")
    O2_DEFINE_CONFIGURABLE(cfgFakeKaonCut, float, 0.1f, "Maximum difference in measured momentum and TPC inner ring momentum of particle")
    O2_DEFINE_CONFIGURABLE(cfgUseAsymmetricPID, bool, false, "Use asymmetric PID cuts")
    O2_DEFINE_CONFIGURABLE(cfgTPCNsigmaCut, float, 3.0f, "TPC N-sigma cut for pions, kaons, protons")
    O2_DEFINE_CONFIGURABLE(cfgUseStrictPID, bool, true, "Use strict PID cuts for TPC")
    O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.5, "pt cut on TOF for PID");
    O2_DEFINE_CONFIGURABLE(cfgUseItsPID, bool, true, "Use ITS PID for particle identification")
    O2_DEFINE_CONFIGURABLE(cfgK0SignalMin, float, 0.48, "Minimum cut on K0 mT signal");
    O2_DEFINE_CONFIGURABLE(cfgK0SignalMax, float, 0.51, "Maximum cut on K0 mT signal");
    O2_DEFINE_CONFIGURABLE(cfgLambdaSignalMin, float, 1.1, "Minimum cut on Lambda mT signal");
    O2_DEFINE_CONFIGURABLE(cfgLambdaSignalMax, float, 1.3, "Maximum cut on Lambda mT signal");
    O2_DEFINE_CONFIGURABLE(cfgK0SideBand1Min, float, 0.44, "Minimum cut on K0 side band 1");
    O2_DEFINE_CONFIGURABLE(cfgK0SideBand1Max, float, 0.47, "Maximum cut on K0 side band 1");
    O2_DEFINE_CONFIGURABLE(cfgK0SideBand2Min, float, 0.52, "Minimum cut on K0 side band 2");
    O2_DEFINE_CONFIGURABLE(cfgK0SideBand2Max, float, 0.56, "Maximum cut on K0 side band 2");
    O2_DEFINE_CONFIGURABLE(cfgLambdaSideBand1Min, float, 1.0, "Minimum cut on Lambda side band 1");
    O2_DEFINE_CONFIGURABLE(cfgLambdaSideBand1Max, float, 1.05, "Maximum cut on Lambda side band 1");
    O2_DEFINE_CONFIGURABLE(cfgLambdaSideBand2Min, float, 1.4, "Minimum cut on Lambda side band 2");
    O2_DEFINE_CONFIGURABLE(cfgLambdaSideBand2Max, float, 1.6, "Maximum cut on Lambda side band 2");
  } cfgPIDCuts;

  Configurable<GFWBinningCuts> cfgGFWBinning{"cfgGFWBinning", {40, 16, 72, 300, 0, 3000, 0.2, 10.0, 0.2, 3.0, {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90}}, "Configuration for binning"};
  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN", "refP", "refFull"}, {-0.8, 0.4, -0.8}, {-0.4, 0.8, 0.8}, {0, 0, 0}, {1, 1, 1}}, "Configurations for GFW regions"};

  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}", "refFull {2 -2}", "refFull {2 2 -2 -2}"}, {"ChGap22", "ChGap32", "ChGap42", "ChFull22", "ChFull24"}, {0, 0, 0, 0, 0}, {15, 1, 1, 0, 0}}, "Configurations for each correlation to calculate"};
  Configurable<GFWCorrConfigs> cfgCorrConfigV02{"cfgCorrConfigV02", {{"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}"}, {"ChGap22", "ChGap32", "ChGap42"}, {1, 1, 1}, {0, 0, 0}}, "Configurations for each radial flow correlation to calculate"};
  Configurable<GFWCorrConfigs> cfgCorrConfigV0{"cfgCorrConfigV0", {{"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}"}, {"ChGap22", "ChGap32", "ChGap42"}, {1, 1, 1}, {0, 0, 0}}, "Configurations for each radial flow correlation to calculate"};

  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};

  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  o2::aod::ITSResponse itsResponse;

  struct Config {
    TH1D* mEfficiency = nullptr;
    std::vector<TH1D*> mPIDEfficiencies;
    std::vector<GFWWeights*> mAcceptance;
    bool correctionsLoaded = false;
  } cfg;

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<FlowPtContainer> fFCpt{FlowPtContainer("FlowPtContainer")};
  OutputObj<FlowContainer> fFCgen{FlowContainer("FlowContainer_gen")};
  HistogramRegistry registry{"registry"};
  HistogramRegistry registryQA{"registryQA"};

  std::array<std::array<float, 11>, 3> resoCutVals;
  std::array<std::array<int, 6>, 3> resoSwitchVals;
  std::array<float, 6> tofNsigmaCut;
  std::array<float, 6> itsNsigmaCut;
  std::array<float, 6> tpcNsigmaCut;

  // QA outputs
  std::map<int, std::vector<std::shared_ptr<TH1>>> th1sList;
  std::map<int, std::vector<std::shared_ptr<TH3>>> th3sList;
  std::vector<std::shared_ptr<TH1>> histosNpt;
  std::vector<std::shared_ptr<TH1>> histosResoNpt;
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
    hPtPhiMult,
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
  struct EventCut {
    bool enabled;
    int histBin;
    int flag; // just store the enum
  };
  std::vector<EventCut> eventcutflags;
  enum Particles {
    PIONS,
    KAONS,
    PROTONS
  };
  enum ParticleIDs {
    CHARGEDID = 0,
    PIONID,
    KAONID,
    PROTONID,
    SPECIESCOUNT
  };
  enum ResoIDs {
    K0SIDEBAND1 = 0,
    K0SIGNAL,
    K0SIDEBAND2,
    LAMBDASIDEBAND1,
    LAMBDASIGNAL,
    LAMBDASIDEBAND2,
    RESOCOUNT
  };
  enum OutputSpecies {
    K0 = 0,
    LAMBDA = 1,
    PHI = 2,
    ANLAMBDA = 3,
    REF = 4,
    kCount_OutputSpecies
  };
  enum ParticleCuts {
    kCosPA = 0,
    kMassMin,
    kMassMax,
    kPosTrackPt,
    kNegTrackPt,
    kDCAPosToPVMin,
    kDCANegToPVMin,
    kLifeTime,
    kRadiusMin,
    kRadiusMax,
    kRapidity
  };
  enum ParticleSwitches {
    kUseParticle = 0,
    kUseCosPA,
    kMassBins,
    kDCABetDaug,
    kUseProperLifetime,
    kUseV0Radius
  };
  enum V0Selection {
    kFillCandidate = 1,
    kFillDaughterPt,
    kFillMassCut,
    kFillRapidityCut,
    kFillDCAtoPV,
    kFillDCAxDaughters,
    kFillV0Radius,
    kFillCosPA,
    kFillProperLifetime,
    kFillDaughterTrackSelection
  };

  // Define global variables
  // Generic Framework
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;

  std::vector<GFW::CorrConfig> corrconfigsV02;
  std::vector<GFW::CorrConfig> corrconfigsV0;

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
    LOGF(info, "FlowGenericFramework::init()");
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
    o2::analysis::gfw::configsV02.SetCorrs(cfgCorrConfigV02->GetCorrs());
    o2::analysis::gfw::configsV02.SetHeads(cfgCorrConfigV02->GetHeads());
    o2::analysis::gfw::configsV02.SetpTDifs(cfgCorrConfigV02->GetpTDifs());
    o2::analysis::gfw::configsV02.SetpTCorrMasks(cfgCorrConfigV02->GetpTCorrMasks());
    o2::analysis::gfw::configsV02.Print();
    o2::analysis::gfw::configsV0.SetCorrs(cfgCorrConfigV0->GetCorrs());
    o2::analysis::gfw::configsV0.SetHeads(cfgCorrConfigV0->GetHeads());
    o2::analysis::gfw::configsV0.SetpTDifs(cfgCorrConfigV0->GetpTDifs());
    o2::analysis::gfw::configsV0.SetpTCorrMasks(cfgCorrConfigV0->GetpTCorrMasks());
    o2::analysis::gfw::configsV0.Print();
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
    LOGF(info, "Eta cuts: Filter [%.1f,%.1f] | Nch [%.1f,%.1f] | Pt-Pt [%.1f, %.1f] | V0 daughters [%.1f, %.1f]", cfgEta->first, cfgEta->second, cfgEtaNch->first, cfgEtaNch->second, cfgEtaPtPt->first, cfgEtaPtPt->second, cfgEtaV0Daughters->first, cfgEtaV0Daughters->second);
    o2::analysis::gfw::multGlobalCorrCutPars = cfgMultCorrCuts.cfgMultGlobalCutPars;
    o2::analysis::gfw::multPVCorrCutPars = cfgMultCorrCuts.cfgMultPVCutPars;
    o2::analysis::gfw::multGlobalPVCorrCutPars = cfgMultCorrCuts.cfgMultGlobalPVCutPars;
    o2::analysis::gfw::multGlobalV0ACutPars = cfgMultCorrCuts.cfgMultGlobalV0ACutPars;
    o2::analysis::gfw::multGlobalT0ACutPars = cfgMultCorrCuts.cfgMultGlobalT0ACutPars;

    projectMatrix(cfgPIDCuts.nSigmas->getData(), tpcNsigmaCut, tofNsigmaCut, itsNsigmaCut);
    readMatrix(cfgPIDCuts.resonanceCuts->getData(), resoCutVals);
    readMatrix(cfgPIDCuts.resonanceSwitches->getData(), resoSwitchVals);
    printResoCuts();

    for (int i = 0; i < 4; ++i) { // o2-linter: disable=magic-number (maximum of 4 subevents)
      if (cfgPtPtGaps->getData()[i][0] < -1. || cfgPtPtGaps->getData()[i][1] < -1.)
        continue;
      o2::analysis::gfw::etagapsPtPt.push_back(std::make_pair(cfgPtPtGaps->getData()[i][0], cfgPtPtGaps->getData()[i][1]));
    }

    for (const auto& [etamin, etamax] : o2::analysis::gfw::etagapsPtPt) {
      LOGF(info, "pt-pt subevent: {%.1f,%.1f}", etamin, etamax);
    }

    // Setup event cuts
    eventcutflags.push_back({cfgEventCutFlags.cfgNoSameBunchPileupCut, kNoSameBunchPileup, o2::aod::evsel::kNoSameBunchPileup});
    eventcutflags.push_back({cfgEventCutFlags.cfgIsGoodZvtxFT0vsPV, kIsGoodZvtxFT0vsPV, o2::aod::evsel::kIsGoodZvtxFT0vsPV});
    eventcutflags.push_back({cfgEventCutFlags.cfgNoCollInTimeRangeStandard, kNoCollInTimeRangeStandard, o2::aod::evsel::kNoCollInTimeRangeStandard});
    eventcutflags.push_back({cfgEventCutFlags.cfgNoCollInRofStandard, kNoCollInRofStandard, o2::aod::evsel::kNoCollInRofStandard});
    eventcutflags.push_back({cfgEventCutFlags.cfgNoHighMultCollInPrevRof, kNoHighMultCollInPrevRof, o2::aod::evsel::kNoHighMultCollInPrevRof});
    eventcutflags.push_back({cfgEventCutFlags.cfgNoTimeFrameBorder, kNoTimeFrameBorder, o2::aod::evsel::kNoTimeFrameBorder});
    eventcutflags.push_back({cfgEventCutFlags.cfgNoITSROFrameBorder, kNoITSROFrameBorder, o2::aod::evsel::kNoITSROFrameBorder});
    eventcutflags.push_back({cfgEventCutFlags.cfgIsVertexITSTPC, kIsVertexITSTPC, o2::aod::evsel::kIsVertexITSTPC});
    eventcutflags.push_back({cfgEventCutFlags.cfgIsGoodITSLayersAll, kIsGoodITSLayersAll, o2::aod::evsel::kIsGoodITSLayersAll});
    for (const auto& cut : eventcutflags) {
      LOGF(info, "Flag %d is %senabled", cut.histBin, (cut.enabled) ? "" : "not ");
    }

    AxisSpec phiAxis = {o2::analysis::gfw::phibins, o2::analysis::gfw::philow, o2::analysis::gfw::phiup, "#phi"};
    AxisSpec phiModAxis = {100, 0, constants::math::PI / 9, "fmod(#varphi,#pi/9)"};
    AxisSpec etaAxis = {o2::analysis::gfw::etabins, cfgEta->first, cfgEta->second, "#eta"};
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
    AxisSpec bAxis = {200, 0, 20, "#it{b}"};
    AxisSpec t0cAxis = {1000, 0, 50000, "N_{ch} (T0C)"};
    AxisSpec t0aAxis = {1800, 0, 180000, "N_{ch} (T0A)"};
    AxisSpec v0aAxis = {1800, 0, 180000, "N_{ch} (V0A)"};
    AxisSpec multpvAxis = {3500, 0, 3500, "N_{ch} (PV)"};
    AxisSpec occAxis = {500, 0, 5000, "occupancy"};
    AxisSpec multAxis = (doprocessOnTheFly && !cfgUseNch) ? bAxis : (cfgUseNch) ? nchAxis
                                                                                : centAxis;
    AxisSpec dcaZAXis = {200, -2, 2, "DCA_{z} (cm)"};
    AxisSpec dcaXYAXis = {200, -1, 1, "DCA_{xy} (cm)"};
    AxisSpec singleCount = {1, 0, 1};
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    int ptbins = o2::analysis::gfw::ptbinning.size() - 1;
    fPtAxis = new TAxis(ptbins, &o2::analysis::gfw::ptbinning[0]);

    if (doprocessMCGen || doprocessOnTheFly) {
      registryQA.add("MCGen/before/pt_gen", "", {HistType::kTH1D, {ptAxis}});
      registryQA.add("MCGen/before/phi_eta_vtxZ_gen", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registryQA.addClone("MCGen/before/", "MCGen/after/");
      if (doprocessOnTheFly)
        registryQA.add("MCGen/impactParameter", "", {HistType::kTH2D, {{bAxis, nchAxis}}});
    }
    if (doprocessMCReco || doprocessData || doprocessRun2) {
      registryQA.add("trackQA/before/phi_eta_vtxZ", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registryQA.add("trackQA/before/pt_dcaXY_dcaZ", "", {HistType::kTH3D, {ptAxis, dcaXYAXis, dcaZAXis}});
      registryQA.add("trackQA/before/pt_phi", "", {HistType::kTH2D, {ptAxis, phiModAxis}});
      registryQA.add("trackQA/before/chi2prTPCcls", "#chi^{2}/cluster for the TPC track segment", {HistType::kTH1D, {{100, 0., 5.}}});
      registryQA.add("trackQA/before/chi2prITScls", "#chi^{2}/cluster for the ITS track", {HistType::kTH1D, {{100, 0., 50.}}});
      registryQA.add("trackQA/before/nTPCClusters", "Number of found TPC clusters", {HistType::kTH1D, {{100, 40, 180}}});
      registryQA.add("trackQA/before/nITSClusters", "Number of found ITS clusters", {HistType::kTH1D, {{100, 0, 20}}});
      registryQA.add("trackQA/before/nTPCCrossedRows", "Number of crossed TPC Rows", {HistType::kTH1D, {{100, 40, 180}}});

      registryQA.addClone("trackQA/before/", "trackQA/after/");
      registryQA.add("trackQA/after/pt_ref", "; #it{p}_{T}; Counts", {HistType::kTH1D, {{100, o2::analysis::gfw::ptreflow, o2::analysis::gfw::ptrefup}}});
      registryQA.add("trackQA/after/pt_poi", "; #it{p}_{T}; Counts", {HistType::kTH1D, {{100, o2::analysis::gfw::ptpoilow, o2::analysis::gfw::ptpoiup}}});
      registryQA.add("trackQA/after/Nch_corrected", "; N_{ch}; Counts", {HistType::kTH1D, {nchAxis}});
      registryQA.add("trackQA/after/Nch_uncorrected", "; N_{ch}; Counts", {HistType::kTH1D, {nchAxis}});
      registryQA.add("trackQA/after/etaNch", "; #eta; Counts", {HistType::kTH1D, {etaAxis}});
      registryQA.add("trackQA/after/etaPtPt", "; #eta; Counts", {HistType::kTH1D, {etaAxis}});
      registryQA.add("trackQA/after/etaV0Daughters", "; #eta; Counts", {HistType::kTH1D, {etaAxis}});

      LOGF(info, "Using alternative filling for pt fractions: %d", static_cast<int>(cfgUseNewNpt));
      histosNpt.resize(SPECIESCOUNT);
      histosNpt[CHARGEDID] = registry.add<TH1>("nptCh", "; #it{p}_{T} (GeV/#it{c}; Count)", {HistType::kTH1D, {ptAxis}});
      histosNpt[PIONID] = registry.add<TH1>("nptPi", "; #it{p}_{T} (GeV/#it{c}; Count)", {HistType::kTH1D, {ptAxis}});
      histosNpt[KAONID] = registry.add<TH1>("nptKa", "; #it{p}_{T} (GeV/#it{c}; Count)", {HistType::kTH1D, {ptAxis}});
      histosNpt[PROTONID] = registry.add<TH1>("nptPr", "; #it{p}_{T} (GeV/#it{c}; Count)", {HistType::kTH1D, {ptAxis}});
      histosResoNpt.resize(RESOCOUNT);
      histosResoNpt[K0SIDEBAND1] = registry.add<TH1>("nptK0SB1", "; #it{p}_{T} (GeV/#it{c}; Count", {HistType::kTH1D, {ptAxis}});
      histosResoNpt[K0SIGNAL] = registry.add<TH1>("nptK0Sig", "; #it{p}_{T} (GeV/#it{c}; Count", {HistType::kTH1D, {ptAxis}});
      histosResoNpt[K0SIDEBAND2] = registry.add<TH1>("nptK0SB2", "; #it{p}_{T} (GeV/#it{c}; Count", {HistType::kTH1D, {ptAxis}});
      histosResoNpt[LAMBDASIDEBAND1] = registry.add<TH1>("nptLambdaSB1", "; #it{p}_{T} (GeV/#it{c}; Count", {HistType::kTH1D, {ptAxis}});
      histosResoNpt[LAMBDASIGNAL] = registry.add<TH1>("nptLambdaSig", "; #it{p}_{T} (GeV/#it{c}; Count", {HistType::kTH1D, {ptAxis}});
      histosResoNpt[LAMBDASIDEBAND2] = registry.add<TH1>("nptLambdaSB2", "; #it{p}_{T} (GeV/#it{c}; Count", {HistType::kTH1D, {ptAxis}});

      registryQA.add("eventQA/before/centrality", "; centrality (%); Counts", {HistType::kTH1D, {centAxis}});
      registryQA.add("eventQA/before/multiplicity", "; N_{ch}; Counts", {HistType::kTH1D, {nchAxis}});
      registryQA.add("eventQA/before/globalTracks_centT0C", "; FT0C centrality (%); N_{global}", {HistType::kTH2D, {centAxis, nchAxis}});
      registryQA.add("eventQA/before/PVTracks_centT0C", "; FT0C centrality (%); N_{PV}", {HistType::kTH2D, {centAxis, multpvAxis}});
      registryQA.add("eventQA/before/globalTracks_PVTracks", "; N_{PV}; N_{global}", {HistType::kTH2D, {multpvAxis, nchAxis}});
      registryQA.add("eventQA/before/globalTracks_multT0A", "; multT0A; N_{global}", {HistType::kTH2D, {t0aAxis, nchAxis}});
      registryQA.add("eventQA/before/globalTracks_multV0A", "; multV0A; N_{global}", {HistType::kTH2D, {v0aAxis, nchAxis}});
      registryQA.add("eventQA/before/multV0A_multT0A", "; multV0A; multT0A", {HistType::kTH2D, {t0aAxis, v0aAxis}});
      registryQA.add("eventQA/before/multT0C_centT0C", "; multT0C; FT0C centrality (%)", {HistType::kTH2D, {centAxis, t0cAxis}});
      registryQA.add("eventQA/before/occ_mult_cent", "; occupancy; N_{ch}; centrality (%)", {HistType::kTH3D, {occAxis, nchAxis, centAxis}});
      registryQA.addClone("eventQA/before/", "eventQA/after/");
      registryQA.add("eventQA/eventSel", "Number of Events;; Counts", {HistType::kTH1D, {{15, 0.5, 15.5}}});
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kFilteredEvent, "Filtered event");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kSel8, "sel8");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kOccupancy, "occupancy");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kTVXinTRD, "kTVXinTRD");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoSameBunchPileup, "kNoSameBunchPileup");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kIsGoodZvtxFT0vsPV, "kIsGoodZvtxFT0vsPV");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoCollInTimeRangeStandard, "kNoCollInTimeRangeStandard");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoCollInRofStandard, "kNoCollInRofStandard");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoHighMultCollInPrevRof, "kNoHighMultCollInPrevRof");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoTimeFrameBorder, "kNoTimeFrameBorder");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kNoITSROFrameBorder, "kNoITSROFrameBorder");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kIsVertexITSTPC, "kIsVertexITSTPC");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kIsGoodITSLayersAll, "kIsGoodITSLayersAll");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kMultCuts, "after Mult cuts");
      registryQA.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(kTrackCent, "has track + within cent");

      registry.add("npt_ch", "; #it{p}_{T} (GeV/#it{c}; ; centrality (%); fraction)", {HistType::kTProfile2D, {ptAxis, centAxis}});
      registry.add("npt_pi", "; #it{p}_{T} (GeV/#it{c}; ; centrality (%); fraction)", {HistType::kTProfile2D, {ptAxis, centAxis}});
      registry.add("npt_ka", "; #it{p}_{T} (GeV/#it{c}; ; centrality (%); fraction)", {HistType::kTProfile2D, {ptAxis, centAxis}});
      registry.add("npt_pr", "; #it{p}_{T} (GeV/#it{c}; ; centrality (%); fraction)", {HistType::kTProfile2D, {ptAxis, centAxis}});
      registry.add("npt_K0_sig", "; #it{p}_{T} (GeV/#it{c}; ; centrality (%); fraction)", {HistType::kTProfile2D, {ptAxis, centAxis}});
      registry.add("npt_K0_sb1", "; #it{p}_{T} (GeV/#it{c}; ; centrality (%); fraction)", {HistType::kTProfile2D, {ptAxis, centAxis}});
      registry.add("npt_K0_sb2", "; #it{p}_{T} (GeV/#it{c}; ; centrality (%); fraction)", {HistType::kTProfile2D, {ptAxis, centAxis}});
      registry.add("npt_Lambda_sig", "; #it{p}_{T} (GeV/#it{c}; ; centrality (%); fraction)", {HistType::kTProfile2D, {ptAxis, centAxis}});
      registry.add("npt_Lambda_sb1", "; #it{p}_{T} (GeV/#it{c}; ; centrality (%); fraction)", {HistType::kTProfile2D, {ptAxis, centAxis}});
      registry.add("npt_Lambda_sb2", "; #it{p}_{T} (GeV/#it{c}; ; centrality (%); fraction)", {HistType::kTProfile2D, {ptAxis, centAxis}});

      if (!cfgRunByRun) {
        if (cfgUsePID) {
          registryQA.add<TH3>("phi_eta_vtxz_ref", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
          registryQA.add<TH3>("phi_eta_vtxz_ch", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
          registryQA.add<TH3>("phi_eta_vtxz_pi", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
          registryQA.add<TH3>("phi_eta_vtxz_ka", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
          registryQA.add<TH3>("phi_eta_vtxz_pr", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
        } else {
          registryQA.add<TH3>("phi_eta_vtxz_ref", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
        }
      }

      AxisSpec axisK0Mass = {resoSwitchVals[K0][kMassBins], resoCutVals[K0][kMassMin], resoCutVals[K0][kMassMax]};
      AxisSpec axisLambdaMass = {resoSwitchVals[LAMBDA][kMassBins], resoCutVals[LAMBDA][kMassMin], resoCutVals[LAMBDA][kMassMax]};

      // QA histograms for V0s
      if (resoSwitchVals[K0][kUseParticle]) {
        registryQA.add("K0/PiPlusTPC_K0", "", {HistType::kTH2D, {{ptAxis, axisNsigmaTPC}}});
        registryQA.add("K0/PiMinusTPC_K0", "", {HistType::kTH2D, {{ptAxis, axisNsigmaTPC}}});
        registryQA.add("K0/PiPlusTOF_K0", "", {HistType::kTH2D, {{ptAxis, axisNsigmaTOF}}});
        registryQA.add("K0/PiMinusTOF_K0", "", {HistType::kTH2D, {{ptAxis, axisNsigmaTOF}}});
        registryQA.add("K0/hK0Phi", "", {HistType::kTH1D, {phiAxis}});
        registryQA.add("K0/hK0Eta", "", {HistType::kTH1D, {etaAxis}});
        registryQA.add("K0/hK0Mass_sparse", "", {HistType::kTHnSparseF, {{axisK0Mass, ptAxis, nchAxis}}});
        registryQA.add("K0/hK0s", "", {HistType::kTH1D, {singleCount}});
        registryQA.add("K0/hK0s_corrected", "", {HistType::kTH1D, {singleCount}});

        registryQA.add("K0/hK0Count", "Number of K0;; Count", {HistType::kTH1D, {{10, 0.5, 10.5}}});
        registryQA.get<TH1>(HIST("K0/hK0Count"))->GetXaxis()->SetBinLabel(kFillCandidate, "K0 candidates");
        registryQA.get<TH1>(HIST("K0/hK0Count"))->GetXaxis()->SetBinLabel(kFillDaughterPt, "Daughter pt");
        registryQA.get<TH1>(HIST("K0/hK0Count"))->GetXaxis()->SetBinLabel(kFillMassCut, "Mass cut");
        registryQA.get<TH1>(HIST("K0/hK0Count"))->GetXaxis()->SetBinLabel(kFillRapidityCut, "Rapidity cut");
        registryQA.get<TH1>(HIST("K0/hK0Count"))->GetXaxis()->SetBinLabel(kFillDCAtoPV, "DCA to PV");
        registryQA.get<TH1>(HIST("K0/hK0Count"))->GetXaxis()->SetBinLabel(kFillDCAxDaughters, "DCA between daughters");
        registryQA.get<TH1>(HIST("K0/hK0Count"))->GetXaxis()->SetBinLabel(kFillV0Radius, "V0radius");
        registryQA.get<TH1>(HIST("K0/hK0Count"))->GetXaxis()->SetBinLabel(kFillCosPA, "CosPA");
        registryQA.get<TH1>(HIST("K0/hK0Count"))->GetXaxis()->SetBinLabel(kFillProperLifetime, "Proper lifetime");
        registryQA.get<TH1>(HIST("K0/hK0Count"))->GetXaxis()->SetBinLabel(kFillDaughterTrackSelection, "Daughter track selection");
      }

      if (resoSwitchVals[LAMBDA][kUseParticle]) {
        registryQA.add("Lambda/PrPlusTPC_L", "", {HistType::kTH2D, {{ptAxis, axisNsigmaTPC}}});
        registryQA.add("Lambda/PiMinusTPC_L", "", {HistType::kTH2D, {{ptAxis, axisNsigmaTPC}}});
        registryQA.add("Lambda/PrPlusTOF_L", "", {HistType::kTH2D, {{ptAxis, axisNsigmaTOF}}});
        registryQA.add("Lambda/PiMinusTOF_L", "", {HistType::kTH2D, {{ptAxis, axisNsigmaTOF}}});
        registryQA.add("Lambda/hLambdaPhi", "", {HistType::kTH1D, {phiAxis}});
        registryQA.add("Lambda/hLambdaEta", "", {HistType::kTH1D, {etaAxis}});
        registryQA.add("Lambda/hLambdaMass_sparse", "", {HistType::kTHnSparseF, {{axisLambdaMass, ptAxis, nchAxis}}});
        registryQA.add("Lambda/PiPlusTPC_AL", "", {HistType::kTH2D, {{ptAxis, axisNsigmaTPC}}});
        registryQA.add("Lambda/PrMinusTPC_AL", "", {HistType::kTH2D, {{ptAxis, axisNsigmaTPC}}});
        registryQA.add("Lambda/PiPlusTOF_AL", "", {HistType::kTH2D, {{ptAxis, axisNsigmaTOF}}});
        registryQA.add("Lambda/PrMinusTOF_AL", "", {HistType::kTH2D, {{ptAxis, axisNsigmaTOF}}});
        registryQA.add("Lambda/hAntiLambdaPhi", "", {HistType::kTH1D, {phiAxis}});
        registryQA.add("Lambda/hAntiLambdaEta", "", {HistType::kTH1D, {etaAxis}});
        registryQA.add("Lambda/hAntiLambdaMass_sparse", "", {HistType::kTHnSparseF, {{axisLambdaMass, ptAxis, nchAxis}}});
        registryQA.add("Lambda/hLambdas", "", {HistType::kTH1D, {singleCount}});
        registryQA.add("Lambda/hLambdas_corrected", "", {HistType::kTH1D, {singleCount}});

        registryQA.add("Lambda/hLambdaCount", "Number of Lambda;; Count", {HistType::kTH1D, {{10, 0.5, 10.5}}});
        registryQA.get<TH1>(HIST("Lambda/hLambdaCount"))->GetXaxis()->SetBinLabel(kFillCandidate, "Lambda candidates");
        registryQA.get<TH1>(HIST("Lambda/hLambdaCount"))->GetXaxis()->SetBinLabel(kFillDaughterPt, "Daughter pt");
        registryQA.get<TH1>(HIST("Lambda/hLambdaCount"))->GetXaxis()->SetBinLabel(kFillMassCut, "Mass cut");
        registryQA.get<TH1>(HIST("Lambda/hLambdaCount"))->GetXaxis()->SetBinLabel(kFillRapidityCut, "Rapidity cut");
        registryQA.get<TH1>(HIST("Lambda/hLambdaCount"))->GetXaxis()->SetBinLabel(kFillDCAtoPV, "DCA to PV");
        registryQA.get<TH1>(HIST("Lambda/hLambdaCount"))->GetXaxis()->SetBinLabel(kFillDCAxDaughters, "DCA between daughters");
        registryQA.get<TH1>(HIST("Lambda/hLambdaCount"))->GetXaxis()->SetBinLabel(kFillV0Radius, "V0radius");
        registryQA.get<TH1>(HIST("Lambda/hLambdaCount"))->GetXaxis()->SetBinLabel(kFillCosPA, "CosPA");
        registryQA.get<TH1>(HIST("Lambda/hLambdaCount"))->GetXaxis()->SetBinLabel(kFillProperLifetime, "Proper lifetime");
        registryQA.get<TH1>(HIST("Lambda/hLambdaCount"))->GetXaxis()->SetBinLabel(kFillDaughterTrackSelection, "Daughter track selection");
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

    // Radial flow configs
    for (auto i = 0; i < o2::analysis::gfw::configsV02.GetSize(); ++i) {
      corrconfigsV02.push_back(fGFW->GetCorrelatorConfig(o2::analysis::gfw::configsV02.GetCorrs()[i], o2::analysis::gfw::configsV02.GetHeads()[i], o2::analysis::gfw::configsV02.GetpTDifs()[i]));
    }
    if (corrconfigsV02.empty())
      LOGF(error, "Radial (V02) configuration contains vectors of different size - check the GFWCorrConfig configurable");
    for (auto i = 0; i < o2::analysis::gfw::configsV0.GetSize(); ++i) {
      corrconfigsV0.push_back(fGFW->GetCorrelatorConfig(o2::analysis::gfw::configsV0.GetCorrs()[i], o2::analysis::gfw::configsV0.GetHeads()[i], o2::analysis::gfw::configsV0.GetpTDifs()[i]));
    }
    if (corrconfigsV0.empty())
      LOGF(error, "Radial (V0) configuration contains vectors of different size - check the GFWCorrConfig configurable");

    fGFW->CreateRegions();
    TObjArray* oba = new TObjArray();
    addConfigObjectsToObjArray(oba, corrconfigs);
    addConfigObjectsToObjArray(oba, corrconfigsV02);
    addConfigObjectsToObjArray(oba, corrconfigsV0);

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
    fFCpt->initialiseSubevent(multAxis, cfgMpar, o2::analysis::gfw::etagapsPtPt.size(), cfgNbootstrap);

    // Multiplicity correlation cuts
    if (cfgMultCut) {
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

      LOGF(info, "Global V0A function: %s in range 0-%g", cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), v0aAxis.binEdges.back());
      fMultGlobalV0ACutLow = new TF1("fMultGlobalV0ACutLow", cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, v0aAxis.binEdges.back());
      for (std::size_t i = 0; i < o2::analysis::gfw::multGlobalV0ACutPars.size(); ++i)
        fMultGlobalV0ACutLow->SetParameter(i, o2::analysis::gfw::multGlobalV0ACutPars[i]);
      fMultGlobalV0ACutLow->SetParameter(o2::analysis::gfw::multGlobalV0ACutPars.size(), cfgMultCorrCuts.cfgGlobalV0ALowSigma);
      for (int i = 0; i < fMultGlobalV0ACutLow->GetNpar(); ++i)
        LOGF(info, "fMultGlobalV0ACutLow par %d = %g", i, fMultGlobalV0ACutLow->GetParameter(i));

      fMultGlobalV0ACutHigh = new TF1("fMultGlobalV0ACutHigh", cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, v0aAxis.binEdges.back());
      for (std::size_t i = 0; i < o2::analysis::gfw::multGlobalV0ACutPars.size(); ++i)
        fMultGlobalV0ACutHigh->SetParameter(i, o2::analysis::gfw::multGlobalV0ACutPars[i]);
      fMultGlobalV0ACutHigh->SetParameter(o2::analysis::gfw::multGlobalV0ACutPars.size(), cfgMultCorrCuts.cfgGlobalV0AHighSigma);
      for (int i = 0; i < fMultGlobalV0ACutHigh->GetNpar(); ++i)
        LOGF(info, "fMultGlobalV0ACutHigh par %d = %g", i, fMultGlobalV0ACutHigh->GetParameter(i));

      LOGF(info, "Global T0A function: %s", cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str());
      fMultGlobalT0ACutLow = new TF1("fMultGlobalT0ACutLow", cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, t0aAxis.binEdges.back());
      for (std::size_t i = 0; i < o2::analysis::gfw::multGlobalT0ACutPars.size(); ++i)
        fMultGlobalT0ACutLow->SetParameter(i, o2::analysis::gfw::multGlobalT0ACutPars[i]);
      fMultGlobalT0ACutLow->SetParameter(o2::analysis::gfw::multGlobalT0ACutPars.size(), cfgMultCorrCuts.cfgGlobalT0ALowSigma);
      for (int i = 0; i < fMultGlobalT0ACutLow->GetNpar(); ++i)
        LOGF(info, "fMultGlobalT0ACutLow par %d = %g", i, fMultGlobalT0ACutLow->GetParameter(i));

      fMultGlobalT0ACutHigh = new TF1("fMultGlobalT0ACutHigh", cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->c_str(), 0, t0aAxis.binEdges.back());
      for (std::size_t i = 0; i < o2::analysis::gfw::multGlobalT0ACutPars.size(); ++i)
        fMultGlobalT0ACutHigh->SetParameter(i, o2::analysis::gfw::multGlobalT0ACutPars[i]);
      fMultGlobalT0ACutHigh->SetParameter(o2::analysis::gfw::multGlobalT0ACutPars.size(), cfgMultCorrCuts.cfgGlobalT0AHighSigma);
      for (int i = 0; i < fMultGlobalT0ACutHigh->GetNpar(); ++i)
        LOGF(info, "fMultGlobalT0ACutHigh par %d = %g", i, fMultGlobalT0ACutHigh->GetParameter(i));
    }

    if (cfgTrackCuts.cfgTPCSectorCut) {
      fPhiCutLow = new TF1("fPhiCutLow", "0.06/x+pi/18.0-0.06", 0, 100);
      fPhiCutHigh = new TF1("fPhiCutHigh", "0.1/x+pi/18.0+0.06", 0, 100);
    }
    // Density dependent corrections
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

  void printResoCuts()
  {
    auto printTable = [](const auto& lbl, const auto& valuesMatrix, const std::string& title) {
      LOGF(info, "===== %s =====", title.c_str());

      size_t nRows = lbl.labels_rows.size();
      size_t nCols = lbl.labels_cols.size();

      // Compute column widths: max of header or value string length
      std::vector<size_t> colWidths(nCols);
      for (size_t c = 0; c < nCols; ++c) {
        size_t maxValLen = 0;
        for (size_t r = 0; r < nRows; ++r) {
          std::ostringstream oss;
          oss << std::fixed << std::setprecision(2) << valuesMatrix[r][c];
          maxValLen = std::max(maxValLen, oss.str().length());
        }
        colWidths[c] = std::max(lbl.labels_cols[c].size(), maxValLen);
      }

      // Determine row label width
      size_t rowLabelWidth = 0;
      for (const auto& rowName : lbl.labels_rows)
        rowLabelWidth = std::max(rowLabelWidth, rowName.size());
      rowLabelWidth += 2; // for ": "

      // Print header
      std::ostringstream header;
      header << std::setw(rowLabelWidth) << " ";
      for (size_t c = 0; c < nCols; ++c)
        header << std::setw(colWidths[c]) << lbl.labels_cols[c] << "  ";
      LOGF(info, "%s", header.str().c_str());

      // Print rows
      for (size_t r = 0; r < nRows; ++r) {
        std::ostringstream line;
        line << std::setw(rowLabelWidth) << (lbl.labels_rows[r] + ":");
        for (size_t c = 0; c < nCols; ++c) {
          line << std::setw(colWidths[c]) << std::fixed << std::setprecision(2) << valuesMatrix[r][c] << "  ";
        }
        LOGF(info, "%s", line.str().c_str());
      }
    };

    // ----- nSigma PID -----
    // Map arrays into a 2D vector
    std::vector<std::vector<float>> nSigmaVals = {
      std::vector<float>(tpcNsigmaCut.begin(), tpcNsigmaCut.end()),
      std::vector<float>(tofNsigmaCut.begin(), tofNsigmaCut.end()),
      std::vector<float>(itsNsigmaCut.begin(), itsNsigmaCut.end())};
    printTable(cfgPIDCuts.nSigmas.value, nSigmaVals, "nSigma PID Cuts");

    // ----- Resonance Cuts -----
    std::vector<std::vector<float>> resoCutsVals(resoCutVals.size());
    for (size_t r = 0; r < resoCutVals.size(); ++r)
      resoCutsVals[r] = std::vector<float>(resoCutVals[r].begin(), resoCutVals[r].end());
    printTable(cfgPIDCuts.resonanceCuts.value, resoCutsVals, "Resonance Cuts");

    // ----- Resonance Switches -----
    std::vector<std::vector<float>> resoSwitchValsF(resoSwitchVals.size());
    for (size_t r = 0; r < resoSwitchVals.size(); ++r)
      resoSwitchValsF[r] = std::vector<float>(resoSwitchVals[r].begin(), resoSwitchVals[r].end());
    printTable(cfgPIDCuts.resonanceSwitches.value, resoSwitchValsF, "Resonance Switches");
  }
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
      if (!cfgUsePIDEfficiencies) {
        cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgEfficiency, timestamp);
        if (cfg.mEfficiency == nullptr) {
          LOGF(fatal, "Could not load efficiency histogram from %s", cfgEfficiency.value.c_str());
        }
        LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)cfg.mEfficiency);
      } else {
        std::vector<std::string> species = {"ch", "pi", "ka", "pr"};
        for (const auto& sp : species) {
          cfg.mPIDEfficiencies.push_back(ccdb->getForTimeStamp<TH1D>(cfgEfficiency.value + "/" + sp, timestamp));
          if (cfg.mPIDEfficiencies.back() == nullptr)
            LOGF(fatal, "Could not load PID efficiency histograms from %s", cfgEfficiency.value + "/" + sp);
          LOGF(info, "Loaded PID efficiency histogram from %s", cfgEfficiency.value + "/" + sp);
        }
      }
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
  double getEfficiency(TTrack track, int pidIndex = 0)
  { //-1 ref, 0 ch, 1 pi, 2 ka, 3 pr
    double eff = 1.;
    if (!cfgUsePIDEfficiencies) {
      if (cfg.mEfficiency)
        eff = cfg.mEfficiency->GetBinContent(cfg.mEfficiency->FindBin(track.pt()));
    } else {
      eff = cfg.mPIDEfficiencies[pidIndex]->GetBinContent(cfg.mPIDEfficiencies[pidIndex]->FindBin(track.pt()));
    }
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
    std::array<float, 3> nSigmaToUse = (track.pt() > cfgPIDCuts.cfgTofPtCut && track.hasTOF()) ? nSigmaCombined : nSigmaTPC;
    if (track.pt() > cfgPIDCuts.cfgTofPtCut && !track.hasTOF())
      return 0;

    const int numSpecies = 3;
    int pidCount = 0;
    // Select particle with the lowest nsigma
    for (int i = 0; i < numSpecies; ++i) {
      if (std::abs(nSigmaToUse[i]) < nsigma) {
        if (pidCount > 0 && cfgPIDCuts.cfgUseStrictPID)
          return 0; // more than one particle with low nsigma

        pidCount++;
        pid = i;
        if (!cfgPIDCuts.cfgUseStrictPID)
          nsigma = std::abs(nSigmaToUse[i]);
      }
    }

    return pid + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
  }

  template <typename TTrack>
  int getNsigmaPIDAssymmetric(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {itsResponse.nSigmaITS<o2::track::PID::Pion>(track), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = -1;

    std::array<float, 3> nSigmaToUse = cfgPIDCuts.cfgUseItsPID ? nSigmaITS : nSigmaTPC;             // Choose which nSigma to use: TPC or ITS
    std::array<float, 6> detectorNsigmaCut = cfgPIDCuts.cfgUseItsPID ? itsNsigmaCut : tpcNsigmaCut; // Choose which nSigma to use: TPC or ITS

    bool isPion, isKaon, isProton;
    bool isDetectedPion = nSigmaToUse[0] < detectorNsigmaCut[0] && nSigmaToUse[0] > detectorNsigmaCut[0 + 3];
    bool isDetectedKaon = nSigmaToUse[1] < detectorNsigmaCut[1] && nSigmaToUse[1] > detectorNsigmaCut[1 + 3];
    bool isDetectedProton = nSigmaToUse[2] < detectorNsigmaCut[2] && nSigmaToUse[2] > detectorNsigmaCut[2 + 3];

    bool isTofPion = nSigmaTOF[0] < tofNsigmaCut[0] && nSigmaTOF[0] > tofNsigmaCut[0 + 3];
    bool isTofKaon = nSigmaTOF[1] < tofNsigmaCut[1] && nSigmaTOF[1] > tofNsigmaCut[1 + 3];
    bool isTofProton = nSigmaTOF[2] < tofNsigmaCut[2] && nSigmaTOF[2] > tofNsigmaCut[2 + 3];

    if (track.pt() > cfgPIDCuts.cfgTofPtCut && !track.hasTOF()) {
      return 0;
    } else if (track.pt() > cfgPIDCuts.cfgTofPtCut && track.hasTOF()) {
      isPion = isTofPion && isDetectedPion;
      isKaon = isTofKaon && isDetectedKaon;
      isProton = isTofProton && isDetectedProton;
    } else {
      isPion = isDetectedPion;
      isKaon = isDetectedKaon;
      isProton = isDetectedProton;
    }

    if ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton)) {
      return 0; // more than one particle satisfy the criteria
    }

    if (isPion) {
      pid = PIONS;
    } else if (isKaon) {
      pid = KAONS;
    } else if (isProton) {
      pid = PROTONS;
    } else {
      return 0; // no particle satisfies the criteria
    }

    return pid + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int multTrk, const float& centrality, const int run)
  {
    // Cut on trigger alias
    if (cfgEventCutFlags.cfgTVXinTRD) {
      if (collision.alias_bit(kTVXinTRD)) {
        // TRD triggered
        // "CMTVX-B-NOPF-TRD,minbias_TVX"
        return false;
      }
      registryQA.fill(HIST("eventQA/eventSel"), kTVXinTRD);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kTVXinTRD);
    }
    // Cut on event selection flags
    for (const auto& cut : eventcutflags) {
      if (!cut.enabled)
        continue;
      if (!collision.selection_bit(cut.flag))
        return false;
      registryQA.fill(HIST("eventQA/eventSel"), cut.histBin);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(cut.histBin);
    }
    // Cut on vertex
    if (!selectVertex(collision))
      return false;
    // Cut on multiplicity correlations - data driven
    if (cfgMultCut) {
      if (!selectMultiplicityCorrelation(collision, multTrk, centrality, run))
        return false;
    }
    return true;
  }

  template <typename TCollision>
  bool selectVertex(TCollision collision)
  {
    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = std::sqrt(collision.covZZ());
      float minZRes = 0.25;
      int minNContrib = 20;
      if (zRes > minZRes && collision.numContrib() < minNContrib)
        vtxz = -999;
    }
    if (vtxz > o2::analysis::gfw::vtxZup || vtxz < o2::analysis::gfw::vtxZlow)
      return false;
    else
      return true;
  }

  template <typename TCollision>
  bool selectMultiplicityCorrelation(TCollision collision, const int multTrk, const float& centrality, const int run)
  {
    auto multNTracksPV = collision.multNTracksPV();
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return false;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return false;
    if (multTrk < fMultCutLow->Eval(centrality))
      return false;
    if (multTrk > fMultCutHigh->Eval(centrality))
      return false;
    if (multTrk > fMultPVGlobalCutHigh->Eval(collision.multNTracksPV()))
      return false;

    if (!(cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->empty()) && static_cast<double>(collision.multFV0A()) < fMultGlobalV0ACutLow->Eval(multTrk))
      return false;
    if (!(cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->empty()) && static_cast<double>(collision.multFV0A()) > fMultGlobalV0ACutHigh->Eval(multTrk))
      return false;
    if (!(cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->empty()) && static_cast<double>(collision.multFT0A()) < fMultGlobalT0ACutLow->Eval(multTrk))
      return false;
    if (!(cfgMultCorrCuts.cfgMultGlobalASideCorrCutFunction->empty()) && static_cast<double>(collision.multFT0A()) > fMultGlobalT0ACutHigh->Eval(multTrk))
      return false;
    registryQA.fill(HIST("eventQA/eventSel"), kMultCuts);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(kMultCuts);
    return true;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track, const int field)
  {
    if (cfgTrackCuts.cfgTPCSectorCut) {
      double phimodn = track.phi();
      if (field < 0) // for negative polarity field
        phimodn = o2::constants::math::TwoPI - phimodn;
      if (track.sign() < 0) // for negative charge
        phimodn = o2::constants::math::TwoPI - phimodn;
      if (phimodn < 0)
        LOGF(warning, "phi < 0: %g", phimodn);

      phimodn += o2::constants::math::PI / 18.0; // to center gap in the middle
      phimodn = fmod(phimodn, o2::constants::math::PI / 9.0);
      if (cfgFillQA)
        registryQA.fill(HIST("trackQA/before/pt_phi"), track.pt(), phimodn);
      if (phimodn < fPhiCutHigh->Eval(track.pt()) && phimodn > fPhiCutLow->Eval(track.pt()))
        return false; // reject track
      if (cfgFillQA)
        registryQA.fill(HIST("trackQA/after/pt_phi"), track.pt(), phimodn);
    }
    if (cfgTrackCuts.cfgDCAxyNSigma && (std::fabs(track.dcaXY()) > cfgTrackCuts.cfgDCAxyNSigma / 7. * (0.0105f + 0.0350f / std::pow(track.pt(), 1.1))))
      return false;
    return ((track.tpcNClsCrossedRows() >= cfgTrackCuts.cfgNTPCXrows) && (track.tpcNClsFound() >= cfgTrackCuts.cfgNTPCCls) && (track.itsNCls() >= cfgTrackCuts.cfgMinNITSCls));
  }

  template <typename TTrack>
  bool nchSelected(TTrack track)
  {
    if (std::fabs(track.dcaXY()) > (0.0105f + 0.035f / std::pow(track.pt(), 1.1)))
      return false;
    return ((track.tpcNClsCrossedRows() >= 70) && (track.tpcNClsFound() >= 50) && (track.itsNCls() >= 5)); // o2-linter: disable=magic-number (hard coded default cuts)
  }

  enum DataType {
    kReco,
    kGen
  };

  template <typename TTrack>
  void fillWeights(const TTrack track, const double vtxz, const int pid_index, const int run)
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
          registryQA.fill(HIST("phi_eta_vtxz_ref"), track.phi(), track.eta(), vtxz); // pt-subset of charged particles for ref flow
        if (withinPtPOI) {
          switch (pid_index) {
            case 0:
              registryQA.fill(HIST("phi_eta_vtxz_ch"), track.phi(), track.eta(), vtxz); // charged particle weights
              break;
            case 1:
              registryQA.fill(HIST("phi_eta_vtxz_pi"), track.phi(), track.eta(), vtxz); // pion weights
              break;
            case 2:
              registryQA.fill(HIST("phi_eta_vtxz_ka"), track.phi(), track.eta(), vtxz); // kaon weights
              break;
            case 3:
              registryQA.fill(HIST("phi_eta_vtxz_pr"), track.phi(), track.eta(), vtxz); // proton weights
              break;
          }
        }
      }
    } else {
      if (cfgRunByRun)
        th3sList[run][hNUAref]->Fill(track.phi(), track.eta(), vtxz);
      else
        registryQA.fill(HIST("phi_eta_vtxz_ref"), track.phi(), track.eta(), vtxz);
    }
    return;
  }

  void createRunByRunHistograms(const int run)
  {
    AxisSpec phiAxis = {o2::analysis::gfw::phibins, o2::analysis::gfw::philow, o2::analysis::gfw::phiup, "#phi"};
    AxisSpec phiModAxis = {100, 0, constants::math::PI / 9, "fmod(#varphi,#pi/9)"};
    AxisSpec etaAxis = {o2::analysis::gfw::etabins, cfgEta->first, cfgEta->second, "#eta"};
    AxisSpec vtxAxis = {o2::analysis::gfw::vtxZbins, -cfgVtxZ, cfgVtxZ, "Vtx_{z} (cm)"};
    AxisSpec nchAxis = {o2::analysis::gfw::nchbins, o2::analysis::gfw::nchlow, o2::analysis::gfw::nchup, "N_{ch}"};
    AxisSpec centAxis = {o2::analysis::gfw::centbinning, "Centrality (%)"};
    AxisSpec ptAxis = {o2::analysis::gfw::ptbinning, "#it{p}_{T} GeV/#it{c}"};
    std::vector<std::shared_ptr<TH1>> histos(kCount_TH1Names);
    histos[hPhi] = registryQA.add<TH1>(Form("%d/phi", run), "", {HistType::kTH1D, {phiAxis}});
    histos[hEta] = registryQA.add<TH1>(Form("%d/eta", run), "", {HistType::kTH1D, {etaAxis}});
    histos[hVtxZ] = registryQA.add<TH1>(Form("%d/vtxz", run), "", {HistType::kTH1D, {vtxAxis}});
    histos[hMult] = registryQA.add<TH1>(Form("%d/mult", run), "", {HistType::kTH1D, {nchAxis}});
    histos[hCent] = registryQA.add<TH1>(Form("%d/cent", run), "", {HistType::kTH1D, {centAxis}});
    histos[hEventSel] = registryQA.add<TH1>(Form("%d/eventSel", run), "Number of Events;; Counts", {HistType::kTH1D, {{11, 0, 11}}});
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
      histos3d[hNUAref] = registryQA.add<TH3>(Form("%d/phi_eta_vtxz_ref", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      histos3d[hNUAch] = registryQA.add<TH3>(Form("%d/phi_eta_vtxz_ch", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      histos3d[hNUApi] = registryQA.add<TH3>(Form("%d/phi_eta_vtxz_pi", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      histos3d[hNUAka] = registryQA.add<TH3>(Form("%d/phi_eta_vtxz_ka", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      histos3d[hNUApr] = registryQA.add<TH3>(Form("%d/phi_eta_vtxz_pr", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
    } else {
      histos3d[hNUAref] = registryQA.add<TH3>(Form("%d/phi_eta_vtxz_ref", run), "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
    }
    histos3d[hPtPhiMult] = registryQA.add<TH3>(Form("%d/pt_phi_mult", run), "", {HistType::kTH3D, {ptAxis, phiModAxis, (cfgUseNch) ? nchAxis : centAxis}});
    th3sList.insert(std::make_pair(run, histos3d));
    return;
  }

  struct AcceptedTracks {
    explicit AcceptedTracks(std::size_t nptbins)
      : pidtotal{0, 0, 0},
        nch(nptbins, 0.f),
        npi(nptbins, 0.f),
        nka(nptbins, 0.f),
        npr(nptbins, 0.f)
    {
    }

    float total = 0;
    unsigned int totaluncorr = 0;

    std::vector<float> pidtotal;
    std::vector<double> nch;
    std::vector<double> npi;
    std::vector<double> nka;
    std::vector<double> npr;
  };

  template <DataType dt>
  void fillOutputContainers(const float& centmult, const double& rndm, AcceptedTracks acceptedtracks)
  {
    fFCpt->calculateCorrelations();
    fFCpt->calculateSubeventCorrelations();
    fFCpt->fillPtProfiles(centmult, rndm);
    fFCpt->fillSubeventPtProfiles(centmult, rndm);
    fFCpt->fillCMProfiles(centmult, rndm);
    fFCpt->fillCMSubeventProfiles(centmult, rndm);
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
          if (cfgUseGapMethod) {
            fFCpt->fillVnPtProfiles(centmult, val, dnx, rndm, o2::analysis::gfw::configs.GetpTCorrMasks()[l_ind]);
          }
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

    double chtotal = (cfgUseNchCorrection) ? acceptedtracks.total : acceptedtracks.totaluncorr;
    // calculate fractions
    std::vector<std::vector<double>> inputs = {acceptedtracks.nch, acceptedtracks.npi, acceptedtracks.nka, acceptedtracks.npr};
    std::vector<std::vector<double>> fractions;
    fractions.reserve(inputs.size());
    int pidcounter = 0;
    for (auto& vec : inputs) { // o2-linter: disable=const-ref-in-for-loop (modified through transform)
      fractions.emplace_back();
      fractions.back().reserve(vec.size());

      double total = chtotal;
      if (cfgUsePIDTotal)
        total = (pidcounter) ? acceptedtracks.pidtotal[pidcounter - 1] : chtotal;

      if (total == 0.) {
        ++pidcounter;
        continue;
      }
      std::transform(vec.begin(), vec.end(),
                     std::back_inserter(fractions.back()),
                     [&](double x) { return x / total; });
      ++pidcounter;
    }

    if (cfgUseNewNpt) {
      for (int i = 1; i <= fPtAxis->GetNbins(); ++i) {
        registry.fill(HIST("npt_ch"), fPtAxis->GetBinCenter(i), centmult, histosNpt[CHARGEDID]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral());
        registry.fill(HIST("npt_pi"), fPtAxis->GetBinCenter(i), centmult, histosNpt[PIONID]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral());
        registry.fill(HIST("npt_ka"), fPtAxis->GetBinCenter(i), centmult, histosNpt[KAONID]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral());
        registry.fill(HIST("npt_pr"), fPtAxis->GetBinCenter(i), centmult, histosNpt[PROTONID]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral());
      }
    } else {
      for (std::size_t i = 0; i < fractions[0].size(); ++i)
        registry.fill(HIST("npt_ch"), fPtAxis->GetBinCenter(i + 1), centmult, fractions[0][i]);
      for (std::size_t i = 0; i < fractions[1].size(); ++i)
        registry.fill(HIST("npt_pi"), fPtAxis->GetBinCenter(i + 1), centmult, fractions[1][i]);
      for (std::size_t i = 0; i < fractions[2].size(); ++i)
        registry.fill(HIST("npt_ka"), fPtAxis->GetBinCenter(i + 1), centmult, fractions[2][i]);
      for (std::size_t i = 0; i < fractions[3].size(); ++i)
        registry.fill(HIST("npt_pr"), fPtAxis->GetBinCenter(i + 1), centmult, fractions[3][i]);
    }
    if (corrconfigsV02.size() < SPECIESCOUNT) //
      return;

    for (uint l_ind = 0; l_ind < SPECIESCOUNT; ++l_ind) {
      for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
        auto dnx = fGFW->Calculate(corrconfigsV02.at(l_ind), i - 1, kTRUE).real();
        if (dnx == 0)
          continue;
        auto val = fGFW->Calculate(corrconfigsV02.at(l_ind), i - 1, kFALSE).real() / dnx;
        if (std::abs(val) < 1)
          (dt == kGen) ? fFCgen->FillProfile(Form("%s_pt_%i", corrconfigsV02.at(l_ind).Head.c_str(), i), centmult, val * ((cfgUseNewNpt) ? histosNpt[l_ind]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral() : fractions[l_ind][i - 1]), dnx, rndm) : fFC->FillProfile(Form("%s_pt_%i", corrconfigsV02.at(l_ind).Head.c_str(), i), centmult, val * ((cfgUseNewNpt) ? histosNpt[l_ind]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral() : fractions[l_ind][i - 1]), dnx, rndm);
      }
    }

    if (corrconfigsV0.size() < SPECIESCOUNT)
      return;

    double mpt = 0;
    if (cfgEtaPtPt->first * cfgEtaPtPt->second >= 0) {
      if (fFCpt->corrDen[1] == 0.)
        return;
      mpt = fFCpt->corrNum[1] / fFCpt->corrDen[1];
    } else {
      if (fFCpt->corrDenSub[0][1] == 0. || fFCpt->corrDenSub[1][1] == 0.)
        return;
      double mpt_sub1 = fFCpt->corrNumSub[0][1] / fFCpt->corrDenSub[0][1];
      double mpt_sub2 = fFCpt->corrNumSub[1][1] / fFCpt->corrDenSub[1][1];
      mpt = 0.5 * (mpt_sub1 + mpt_sub2);
    }
    if (std::isnan(mpt))
      return;
    for (uint l_ind = 0; l_ind < SPECIESCOUNT; ++l_ind) {
      for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
        (dt == kGen) ? fFCgen->FillProfile(Form("%s_pt_%i", corrconfigsV0.at(l_ind).Head.c_str(), i), centmult, mpt * ((cfgUseNewNpt) ? histosNpt[l_ind]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral() : fractions[l_ind][i - 1]), 1., rndm) : fFC->FillProfile(Form("%s_pt_%i", corrconfigsV0.at(l_ind).Head.c_str(), i), centmult, mpt * ((cfgUseNewNpt) ? histosNpt[l_ind]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral() : fractions[l_ind][i - 1]), 1., rndm);
      }
    }
    return;
  }

  template <DataType dt, typename TCollision, typename TTracks, typename TV0s>
  void processCollision(TCollision collision, TTracks tracks, TV0s v0s, const float& centrality, const int field, const int run)
  {
    if (tracks.size() < 1)
      return;
    if (dt != kGen && (centrality < o2::analysis::gfw::centbinning.front() || centrality > o2::analysis::gfw::centbinning.back()))
      return;
    if (dt != kGen) {
      registryQA.fill(HIST("eventQA/eventSel"), kTrackCent);
      if (cfgRunByRun)
        th1sList[run][hEventSel]->Fill(kTrackCent);
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
    // process tracks
    AcceptedTracks acceptedTracks(o2::analysis::gfw::ptbinning.size() - 1);
    // Reset fraction histograms per event
    for (const auto& h : histosNpt)
      h->Reset("ICESM");
    for (const auto& track : tracks) {
      processTrack(track, vtxz, field, run, densitycorrections, acceptedTracks);
    }
    registryQA.fill(HIST("trackQA/after/Nch_corrected"), acceptedTracks.total);
    registryQA.fill(HIST("trackQA/after/Nch_uncorrected"), acceptedTracks.totaluncorr);

    int multiplicity = 0;
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

    if (cfgFillWeights)
      return;

    fillOutputContainers<dt>((cfgUseNch) ? multiplicity : centrality, lRandom, acceptedTracks);
    // Reset fraction histograms per event
    for (const auto& h : histosResoNpt)
      h->Reset("ICESM");
    std::vector<std::vector<float>> nptResonances(6, std::vector<float>(o2::analysis::gfw::ptbinning.size()));
    // Process V0s
    for (const auto& v0 : v0s) {
      if (resoSwitchVals[K0][kUseParticle]) {
        double weff = 1;
        if (selectK0(collision, v0, centrality, weff)) {
          int ptBinIndex = fPtAxis->FindBin(v0.pt()) - 1;
          if (!(ptBinIndex < 0 || ptBinIndex >= static_cast<int>(o2::analysis::gfw::ptbinning.size()))) {
            if (v0.mK0Short() > cfgPIDCuts.cfgK0SideBand1Min && v0.mK0Short() < cfgPIDCuts.cfgK0SideBand1Max) {
              nptResonances[0][ptBinIndex] += (cfgUseNchCorrection) ? weff : 1.0;
              histosResoNpt[K0SIDEBAND1]->Fill(v0.pt(), (cfgUseNchCorrection) ? weff : 1.0);
            }
            if (v0.mK0Short() > cfgPIDCuts.cfgK0SignalMin && v0.mK0Short() < cfgPIDCuts.cfgK0SignalMax) {
              nptResonances[1][ptBinIndex] += (cfgUseNchCorrection) ? weff : 1.0;
              histosResoNpt[K0SIGNAL]->Fill(v0.pt(), (cfgUseNchCorrection) ? weff : 1.0);
            }
            if (v0.mK0Short() > cfgPIDCuts.cfgK0SideBand2Min && v0.mK0Short() < cfgPIDCuts.cfgK0SideBand2Max) {
              nptResonances[2][ptBinIndex] += (cfgUseNchCorrection) ? weff : 1.0;
              histosResoNpt[K0SIDEBAND2]->Fill(v0.pt(), (cfgUseNchCorrection) ? weff : 1.0);
            }
          }
        }
      }
      // Add lambdabar
      if (resoSwitchVals[LAMBDA][kUseParticle]) {
        double weff = 1.;
        if (selectLambda(collision, v0, centrality, weff)) {
          int ptBinIndex = fPtAxis->FindBin(v0.pt()) - 1;
          if (!(ptBinIndex < 0 || ptBinIndex >= static_cast<int>(o2::analysis::gfw::ptbinning.size()))) {
            if (v0.mLambda() > cfgPIDCuts.cfgLambdaSideBand1Min && v0.mLambda() < cfgPIDCuts.cfgLambdaSideBand1Max) {
              nptResonances[3][ptBinIndex] += (cfgUseNchCorrection) ? weff : 1.0;
              histosResoNpt[LAMBDASIDEBAND1]->Fill(v0.pt(), (cfgUseNchCorrection) ? weff : 1.0);
            }
            if (v0.mLambda() > cfgPIDCuts.cfgLambdaSignalMin && v0.mLambda() < cfgPIDCuts.cfgLambdaSignalMax) {
              nptResonances[4][ptBinIndex] += (cfgUseNchCorrection) ? weff : 1.0;
              histosResoNpt[LAMBDASIGNAL]->Fill(v0.pt(), (cfgUseNchCorrection) ? weff : 1.0);
            }
            if (v0.mLambda() > cfgPIDCuts.cfgLambdaSideBand2Min && v0.mLambda() < cfgPIDCuts.cfgLambdaSideBand2Max) {
              nptResonances[5][ptBinIndex] += (cfgUseNchCorrection) ? weff : 1.0;
              histosResoNpt[LAMBDASIDEBAND2]->Fill(v0.pt(), (cfgUseNchCorrection) ? weff : 1.0);
            }
          }
        }
      }
    }
    double chtotal = (cfgUseNchCorrection) ? acceptedTracks.total : acceptedTracks.totaluncorr;
    // calculate fractions
    std::vector<std::vector<float>> fractions_resonances = nptResonances;
    int pidcounter = 0;
    for (auto& vec : fractions_resonances) { // o2-linter: disable=const-ref-in-for-loop (modified through transform)
      double total = chtotal;
      if (cfgUsePIDTotal)
        total = (pidcounter) ? std::accumulate(vec.begin(), vec.end(), 0.f) : chtotal;

      if (total == 0.) {
        ++pidcounter;
        continue;
      }
      std::transform(vec.begin(), vec.end(), vec.begin(),
                     [&](float x) { return x / total; });
      ++pidcounter;
    }

    if (cfgUseNewNpt) {
      for (int i = 1; i <= fPtAxis->GetNbins(); ++i) {
        if (histosNpt[CHARGEDID]->Integral() <= 0)
          continue;
        registry.fill(HIST("npt_K0_sb1"), fPtAxis->GetBinCenter(i), centrality, histosResoNpt[K0SIDEBAND1]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral());
        registry.fill(HIST("npt_K0_sig"), fPtAxis->GetBinCenter(i), centrality, histosResoNpt[K0SIGNAL]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral());
        registry.fill(HIST("npt_K0_sb2"), fPtAxis->GetBinCenter(i), centrality, histosResoNpt[K0SIDEBAND2]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral());
        registry.fill(HIST("npt_Lambda_sb1"), fPtAxis->GetBinCenter(i), centrality, histosResoNpt[LAMBDASIDEBAND1]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral());
        registry.fill(HIST("npt_Lambda_sig"), fPtAxis->GetBinCenter(i), centrality, histosResoNpt[LAMBDASIGNAL]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral());
        registry.fill(HIST("npt_Lambda_sb2"), fPtAxis->GetBinCenter(i), centrality, histosResoNpt[LAMBDASIDEBAND2]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral());
      }
    } else {
      for (std::size_t i = 0; i < fractions_resonances[0].size(); ++i)
        registry.fill(HIST("npt_K0_sb1"), fPtAxis->GetBinCenter(i + 1), centrality, fractions_resonances[0][i]);
      for (std::size_t i = 0; i < fractions_resonances[2].size(); ++i)
        registry.fill(HIST("npt_K0_sb2"), fPtAxis->GetBinCenter(i + 1), centrality, fractions_resonances[2][i]);
      for (std::size_t i = 0; i < fractions_resonances[1].size(); ++i)
        registry.fill(HIST("npt_K0_sig"), fPtAxis->GetBinCenter(i + 1), centrality, fractions_resonances[1][i]);
      for (std::size_t i = 0; i < fractions_resonances[3].size(); ++i)
        registry.fill(HIST("npt_Lambda_sb1"), fPtAxis->GetBinCenter(i + 1), centrality, fractions_resonances[3][i]);
      for (std::size_t i = 0; i < fractions_resonances[5].size(); ++i)
        registry.fill(HIST("npt_Lambda_sb2"), fPtAxis->GetBinCenter(i + 1), centrality, fractions_resonances[5][i]);
      for (std::size_t i = 0; i < fractions_resonances[4].size(); ++i)
        registry.fill(HIST("npt_Lambda_sig"), fPtAxis->GetBinCenter(i + 1), centrality, fractions_resonances[4][i]);
    }
    for (uint l_ind = 4; l_ind < corrconfigsV02.size(); ++l_ind) {
      if (histosNpt[CHARGEDID]->Integral() <= 0 && cfgUseNewNpt)
        continue;
      for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
        auto dnx = fGFW->Calculate(corrconfigsV02.at(l_ind), i - 1, kTRUE).real();
        if (dnx == 0)
          continue;
        auto val = fGFW->Calculate(corrconfigsV02.at(l_ind), i - 1, kFALSE).real() / dnx;
        if (std::abs(val) < 1)
          (dt == kGen) ? fFCgen->FillProfile(Form("%s_pt_%i", corrconfigsV02.at(l_ind).Head.c_str(), i), centrality, val * ((cfgUseNewNpt) ? histosResoNpt[l_ind - 4]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral() : fractions_resonances[l_ind - 4][i - 1]), dnx, lRandom) : fFC->FillProfile(Form("%s_pt_%i", corrconfigsV02.at(l_ind).Head.c_str(), i), centrality, val * ((cfgUseNewNpt) ? histosResoNpt[l_ind - 4]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral() : fractions_resonances[l_ind - 4][i - 1]), dnx, lRandom);
      }
    }

    if (fFCpt->corrDenSub[0][1] == 0. || fFCpt->corrDenSub[1][1] == 0.)
      return;

    double mpt_sub1 = fFCpt->corrNumSub[0][1] / fFCpt->corrDenSub[0][1];
    double mpt_sub2 = fFCpt->corrNumSub[1][1] / fFCpt->corrDenSub[1][1];
    double mpt = 0.5 * (mpt_sub1 + mpt_sub2);
    if (std::isnan(mpt))
      return;

    for (uint l_ind = 4; l_ind < corrconfigsV0.size(); ++l_ind) {
      if (histosNpt[CHARGEDID]->Integral() <= 0 && cfgUseNewNpt)
        continue;
      for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
        (dt == kGen) ? fFCgen->FillProfile(Form("%s_pt_%i", corrconfigsV0.at(l_ind).Head.c_str(), i), centrality, mpt * ((cfgUseNewNpt) ? histosResoNpt[l_ind - 4]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral() : fractions_resonances[l_ind - 4][i - 1]), 1.0, lRandom) : fFC->FillProfile(Form("%s_pt_%i", corrconfigsV0.at(l_ind).Head.c_str(), i), centrality, mpt * ((cfgUseNewNpt) ? histosResoNpt[l_ind - 4]->GetBinContent(i) / histosNpt[CHARGEDID]->Integral() : fractions_resonances[l_ind - 4][i - 1]), 1.0, lRandom);
      }
    }
  }

  template <typename TTrack>
  inline void processTrack(TTrack const& track, const float& vtxz, const int field, const int run, DensityCorr densitycorrections, AcceptedTracks& acceptedTracks)
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

      // Select tracks with nominal cuts always
      if (!nchSelected(track))
        return;

      double weffCh = getEfficiency(track, 0);
      if (track.eta() > cfgEtaNch->first && track.eta() < cfgEtaNch->second) {
        if (weffCh > 0)
          acceptedTracks.total += (cfgUseNchCorrection) ? weffCh : 1.0;
        ++acceptedTracks.totaluncorr;
      }

      if (!trackSelected(track, field))
        return;

      int pidIndex = 0;
      if (cfgUsePID) {
        if (std::abs(mcParticle.pdgCode()) == kPiPlus)
          pidIndex = PIONID;
        if (std::abs(mcParticle.pdgCode()) == kKPlus)
          pidIndex = KAONID;
        if (std::abs(mcParticle.pdgCode()) == kProton)
          pidIndex = PROTONID;
      }

      if (track.eta() > cfgEtaNch->first && track.eta() < cfgEtaNch->second) {
        double weff = getEfficiency(track, pidIndex);

        if (pidIndex && weff > 0)
          acceptedTracks.pidtotal[pidIndex - 1] += weff;

        int ptBinIndex = fPtAxis->FindBin(track.pt()) - 1;
        if (!(ptBinIndex < 0 || ptBinIndex >= static_cast<int>(o2::analysis::gfw::ptbinning.size()))) {
          if (weffCh > 0) {
            acceptedTracks.nch[ptBinIndex] += (cfgUseNchCorrection) ? weffCh : 1.0;
            histosNpt[CHARGEDID]->Fill(track.pt(), (cfgUseNchCorrection) ? weffCh : 1.0);
          }
          if (pidIndex == PIONID && weff > 0) {
            acceptedTracks.npi[ptBinIndex] += (cfgUseNchCorrection) ? weff : 1.0;
            histosNpt[PIONID]->Fill(track.pt(), (cfgUseNchCorrection) ? weff : 1.0);
          }
          if (pidIndex == KAONID && weff > 0) {
            acceptedTracks.nka[ptBinIndex] += (cfgUseNchCorrection) ? weff : 1.0;
            histosNpt[KAONID]->Fill(track.pt(), (cfgUseNchCorrection) ? weff : 1.0);
          }
          if (pidIndex == PROTONID && weff > 0) {
            acceptedTracks.npr[ptBinIndex] += (cfgUseNchCorrection) ? weff : 1.0;
            histosNpt[PROTONID]->Fill(track.pt(), (cfgUseNchCorrection) ? weff : 1.0);
          }
        }
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

      if (track.eta() > cfgEtaNch->first && track.eta() < cfgEtaNch->second) {
        ++acceptedTracks.total;
        ++acceptedTracks.totaluncorr;

        if (pidIndex)
          acceptedTracks.pidtotal[pidIndex - 1] += 1;
        int ptBinIndex = fPtAxis->FindBin(track.pt()) - 1;

        if (!(ptBinIndex < 0 || ptBinIndex >= static_cast<int>(o2::analysis::gfw::ptbinning.size()))) {
          acceptedTracks.nch[ptBinIndex] += 1.0;
          histosNpt[CHARGEDID]->Fill(track.pt());
          if (pidIndex == PIONID) {
            acceptedTracks.npi[ptBinIndex] += 1.0;
            histosNpt[PIONID]->Fill(track.pt());
          }
          if (pidIndex == KAONID) {
            acceptedTracks.nka[ptBinIndex] += 1.0;
            histosNpt[KAONID]->Fill(track.pt());
          }
          if (pidIndex == PROTONID) {
            acceptedTracks.npr[ptBinIndex] += 1.0;
            histosNpt[PROTONID]->Fill(track.pt());
          }
        }
      }

      fillPtSums<kGen>(track, vtxz);
      fillGFW<kGen>(track, vtxz, pidIndex, densitycorrections);

      if (cfgFillQA)
        fillTrackQA<kGen, kAfter>(track, vtxz);
    } else {
      if (cfgFillQA)
        fillTrackQA<kReco, kBefore>(track, vtxz);
      // Select tracks with nominal cuts always
      if (!nchSelected(track))
        return;

      double weffCh = getEfficiency(track, 0);
      if (track.eta() > cfgEtaNch->first && track.eta() < cfgEtaNch->second) {
        if (weffCh > 0)
          acceptedTracks.total += (cfgUseNchCorrection) ? weffCh : 1.0;
        ++acceptedTracks.totaluncorr;
      }

      if (!trackSelected(track, field))
        return;
      // int pidIndex = 0;
      // if (cfgUsePID) Need PID for v02
      int pidIndex = getNsigmaPID(track);

      if (track.eta() > cfgEtaNch->first && track.eta() < cfgEtaNch->second) {
        double weff = getEfficiency(track, pidIndex);
        if (pidIndex && weff > 0)
          acceptedTracks.pidtotal[pidIndex - 1] += weff;

        int ptBinIndex = fPtAxis->FindBin(track.pt()) - 1;

        if (!(ptBinIndex < 0 || ptBinIndex >= static_cast<int>(o2::analysis::gfw::ptbinning.size()))) {
          if (weffCh > 0) {
            acceptedTracks.nch[ptBinIndex] += (cfgUseNchCorrection) ? weffCh : 1.0;
            histosNpt[CHARGEDID]->Fill(track.pt(), (cfgUseNchCorrection) ? weffCh : 1.0);
          }
          if (pidIndex == PIONID && weff > 0) {
            acceptedTracks.npi[ptBinIndex] += (cfgUseNchCorrection) ? weff : 1.0;
            histosNpt[PIONID]->Fill(track.pt(), (cfgUseNchCorrection) ? weff : 1.0);
          }
          if (pidIndex == KAONID && weff > 0) {
            acceptedTracks.nka[ptBinIndex] += (cfgUseNchCorrection) ? weff : 1.0;
            histosNpt[KAONID]->Fill(track.pt(), (cfgUseNchCorrection) ? weff : 1.0);
          }
          if (pidIndex == PROTONID && weff > 0) {
            acceptedTracks.npr[ptBinIndex] += (cfgUseNchCorrection) ? weff : 1.0;
            histosNpt[PROTONID]->Fill(track.pt(), (cfgUseNchCorrection) ? weff : 1.0);
          }
        }
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

  template <typename TTrack>
  bool selectionV0Daughter(TTrack const& track, int pid)
  {
    if (!(track.itsNCls() > cfgTrackCuts.cfgMinNITSCls))
      return 0;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < cfgTrackCuts.cfgNTPCCls)
      return false;
    if (!(track.tpcNClsCrossedRows() > cfgTrackCuts.cfgNTPCXrows))
      return 0;

    if (cfgPIDCuts.cfgUseOnlyTPC) {
      if (pid == PIONS && std::abs(track.tpcNSigmaPi()) > cfgPIDCuts.cfgTPCNsigmaCut)
        return false;
      if (pid == KAONS && std::abs(track.tpcNSigmaKa()) > cfgPIDCuts.cfgTPCNsigmaCut)
        return false;
      if (pid == PROTONS && std::abs(track.tpcNSigmaPr()) > cfgPIDCuts.cfgTPCNsigmaCut)
        return false;
    } else {
      int partIndex = cfgPIDCuts.cfgUseAsymmetricPID ? getNsigmaPIDAssymmetric(track) : getNsigmaPID(track);
      int pidIndex = partIndex - 1; // 0 = pion, 1 = kaon, 2 = proton
      if (pidIndex != pid)
        return false;
    }

    // Eta cuts on daughter particles to remove self-correlations with correlated observables
    if (track.eta() < cfgEtaV0Daughters->first || track.eta() > cfgEtaV0Daughters->second)
      return false;
    registryQA.fill(HIST("trackQA/after/etaV0Daughters"), track.eta());
    return true;
  }

  template <typename TCollision, typename TV0>
  bool selectK0(TCollision const& collision, TV0 const& v0, const double& centrality, double& weff)
  {

    double massK0s = v0.mK0Short();

    auto postrack = v0.template posTrack_as<GFWTracks>();
    auto negtrack = v0.template negTrack_as<GFWTracks>();

    registryQA.fill(HIST("K0/hK0Count"), kFillCandidate);
    if (postrack.pt() < resoCutVals[K0][kPosTrackPt] || negtrack.pt() < resoCutVals[K0][kNegTrackPt])
      return false;
    registryQA.fill(HIST("K0/hK0Count"), kFillDaughterPt);
    if (massK0s < resoCutVals[K0][kMassMin] && massK0s > resoCutVals[K0][kMassMax])
      return false;
    registryQA.fill(HIST("K0/hK0Count"), kFillMassCut);
    // Rapidity correction
    if (v0.yK0Short() > resoCutVals[K0][kRapidity])
      return false;
    registryQA.fill(HIST("K0/hK0Count"), kFillRapidityCut);
    // DCA cuts for K0short
    if (std::abs(v0.dcapostopv()) < resoCutVals[K0][kDCAPosToPVMin] || std::abs(v0.dcanegtopv()) < resoCutVals[K0][kDCANegToPVMin])
      return false;
    registryQA.fill(HIST("K0/hK0Count"), kFillDCAtoPV);
    if (std::abs(v0.dcaV0daughters()) > resoSwitchVals[K0][kDCABetDaug])
      return false;
    registryQA.fill(HIST("K0/hK0Count"), kFillDCAxDaughters);
    // v0 radius cuts
    if (resoSwitchVals[K0][kUseV0Radius] && (v0.v0radius() < resoCutVals[K0][kRadiusMin] || v0.v0radius() > resoCutVals[K0][kRadiusMax]))
      return false;
    registryQA.fill(HIST("K0/hK0Count"), kFillV0Radius);
    // cosine pointing angle cuts
    if (v0.v0cosPA() < resoCutVals[K0][kCosPA])
      return false;
    registryQA.fill(HIST("K0/hK0Count"), kFillCosPA);
    // Proper lifetime
    if (resoSwitchVals[K0][kUseProperLifetime] && v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massK0Short > resoCutVals[K0][kLifeTime])
      return false;
    registryQA.fill(HIST("K0/hK0Count"), kFillProperLifetime);
    if (!selectionV0Daughter(postrack, PIONS) || !selectionV0Daughter(negtrack, PIONS))
      return false;
    registryQA.fill(HIST("K0/hK0Count"), kFillDaughterTrackSelection);

    registryQA.fill(HIST("K0/hK0Mass_sparse"), massK0s, v0.pt(), centrality);
    registryQA.fill(HIST("K0/hK0Phi"), v0.phi());
    registryQA.fill(HIST("K0/hK0Eta"), v0.eta());
    registryQA.fill(HIST("K0/PiPlusTPC_K0"), postrack.pt(), postrack.tpcNSigmaKa());
    registryQA.fill(HIST("K0/PiPlusTOF_K0"), postrack.pt(), postrack.tofNSigmaKa());
    registryQA.fill(HIST("K0/PiMinusTPC_K0"), negtrack.pt(), negtrack.tpcNSigmaKa());
    registryQA.fill(HIST("K0/PiMinusTOF_K0"), negtrack.pt(), negtrack.tofNSigmaKa());

    registryQA.fill(HIST("K0/hK0s"), 0.5, 1);
    if (cfgUsePIDEfficiencies) {
      double weffDaughter1 = getEfficiency(postrack, 1);
      double weffDaughter2 = getEfficiency(negtrack, 1);
      weff = weffDaughter1 * weffDaughter2;
      if (weff > 0)
        registryQA.fill(HIST("K0/hK0s_corrected"), 0.5, weff);
    }

    return true;
  }

  template <typename TCollision, typename TV0>
  bool selectLambda(TCollision const& collision, TV0 const& v0, const double& centrality, double& weff)
  {
    bool isL = false;  // Is lambda candidate
    bool isAL = false; // Is anti-lambda candidate

    double mlambda = v0.mLambda();
    double mantilambda = v0.mAntiLambda();

    // separate the positive and negative V0 daughters
    auto postrack = v0.template posTrack_as<GFWTracks>();
    auto negtrack = v0.template negTrack_as<GFWTracks>();

    registryQA.fill(HIST("Lambda/hLambdaCount"), kFillCandidate);
    if (postrack.pt() < resoCutVals[LAMBDA][kPosTrackPt] || negtrack.pt() < resoCutVals[LAMBDA][kNegTrackPt])
      return false;

    registryQA.fill(HIST("Lambda/hLambdaCount"), kFillDaughterPt);
    if (mlambda > resoCutVals[LAMBDA][kMassMin] && mlambda < resoCutVals[LAMBDA][kMassMax])
      isL = true;
    if (mantilambda > resoCutVals[LAMBDA][kMassMin] && mantilambda < resoCutVals[LAMBDA][kMassMax])
      isAL = true;

    if (!isL && !isAL) {
      return false;
    }
    registryQA.fill(HIST("Lambda/hLambdaCount"), kFillMassCut);

    // Rapidity correction
    if (v0.yLambda() > resoCutVals[LAMBDA][kRapidity])
      return false;
    registryQA.fill(HIST("Lambda/hLambdaCount"), kFillRapidityCut);
    // DCA cuts for lambda and antilambda
    if (isL) {
      if (std::abs(v0.dcapostopv()) < resoCutVals[LAMBDA][kDCAPosToPVMin] || std::abs(v0.dcanegtopv()) < resoCutVals[LAMBDA][kDCANegToPVMin])
        return false;
    }
    if (isAL) {
      if (std::abs(v0.dcapostopv()) < resoCutVals[LAMBDA][kDCANegToPVMin] || std::abs(v0.dcanegtopv()) < resoCutVals[LAMBDA][kDCAPosToPVMin])
        return false;
    }
    registryQA.fill(HIST("Lambda/hLambdaCount"), kFillDCAtoPV);
    if (std::abs(v0.dcaV0daughters()) > resoSwitchVals[LAMBDA][kDCABetDaug])
      return false;
    registryQA.fill(HIST("Lambda/hLambdaCount"), kFillDCAxDaughters);
    // v0 radius cuts
    if (resoSwitchVals[LAMBDA][kUseV0Radius] && (v0.v0radius() < resoCutVals[LAMBDA][kRadiusMin] || v0.v0radius() > resoCutVals[LAMBDA][kRadiusMax]))
      return false;
    registryQA.fill(HIST("Lambda/hLambdaCount"), kFillV0Radius);
    // cosine pointing angle cuts
    if (v0.v0cosPA() < resoCutVals[LAMBDA][kCosPA])
      return false;
    registryQA.fill(HIST("Lambda/hLambdaCount"), kFillCosPA);
    // Proper lifetime
    if (resoSwitchVals[LAMBDA][kUseProperLifetime] && v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda > resoCutVals[LAMBDA][kLifeTime])
      return false;
    registryQA.fill(HIST("Lambda/hLambdaCount"), kFillProperLifetime);
    if (isL) {
      if (!selectionV0Daughter(postrack, PROTONS) || !selectionV0Daughter(negtrack, PIONS))
        return false;
    }
    if (isAL) {
      if (!selectionV0Daughter(postrack, PIONS) || !selectionV0Daughter(negtrack, PROTONS))
        return false;
    }
    registryQA.fill(HIST("Lambda/hLambdaCount"), kFillDaughterTrackSelection);

    if (isL) {
      registryQA.fill(HIST("Lambda/hLambdaMass_sparse"), mlambda, v0.pt(), centrality);
      registryQA.fill(HIST("Lambda/hLambdaPhi"), v0.phi());
      registryQA.fill(HIST("Lambda/hLambdaEta"), v0.eta());
      registryQA.fill(HIST("Lambda/PrPlusTPC_L"), postrack.pt(), postrack.tpcNSigmaKa());
      registryQA.fill(HIST("Lambda/PrPlusTOF_L"), postrack.pt(), postrack.tofNSigmaKa());
      registryQA.fill(HIST("Lambda/PiMinusTPC_L"), negtrack.pt(), negtrack.tpcNSigmaKa());
      registryQA.fill(HIST("Lambda/PiMinusTOF_L"), negtrack.pt(), negtrack.tofNSigmaKa());

      registryQA.fill(HIST("Lambda/hLambdas"), 0.5, 1);
      if (cfgUsePIDEfficiencies) {
        double weffDaughter1 = getEfficiency(postrack, 3);
        double weffDaughter2 = getEfficiency(negtrack, 1);
        weff = weffDaughter1 * weffDaughter2;
        if (weff > 0)
          registryQA.fill(HIST("Lambda/hLambdas_corrected"), 0.5, weff);
      }
    }
    if (isAL) {
      registryQA.fill(HIST("Lambda/hAntiLambdaMass_sparse"), mantilambda, v0.pt(), centrality);
      registryQA.fill(HIST("Lambda/hAntiLambdaPhi"), v0.phi());
      registryQA.fill(HIST("Lambda/hAntiLambdaEta"), v0.eta());
      registryQA.fill(HIST("Lambda/PiPlusTPC_AL"), postrack.pt(), postrack.tpcNSigmaKa());
      registryQA.fill(HIST("Lambda/PiPlusTOF_AL"), postrack.pt(), postrack.tofNSigmaKa());
      registryQA.fill(HIST("Lambda/PrMinusTPC_AL"), negtrack.pt(), negtrack.tpcNSigmaKa());
      registryQA.fill(HIST("Lambda/PrMinusTOF_AL"), negtrack.pt(), negtrack.tofNSigmaKa());

      registryQA.fill(HIST("Lambda/hLambdas"), 0.5, 1);
      if (cfgUsePIDEfficiencies) {
        double weffDaughter1 = getEfficiency(postrack, 1);
        double weffDaughter2 = getEfficiency(negtrack, 3);
        weff = weffDaughter1 * weffDaughter2;
        if (weff > 0)
          registryQA.fill(HIST("Lambda/hLambdas_corrected"), 0.5, weff);
      }
    }
    return true;
  }

  template <DataType dt, typename TTrack>
  inline void fillPtSums(TTrack track, const double& vtxz)
  {
    double wacc = (dt == kGen) ? 1. : getAcceptance(track, vtxz, 0);
    double weff = (dt == kGen) ? 1. : getEfficiency(track);
    if (weff < 0)
      return;

    // Fill the nominal sums
    if (track.eta() > cfgEtaPtPt->first && track.eta() < cfgEtaPtPt->second)
      fFCpt->fill(weff, track.pt());

    // Fill the subevent sums
    std::size_t index = 0;
    for (const auto& [etamin, etamax] : o2::analysis::gfw::etagapsPtPt) {
      if (etamin < track.eta() && track.eta() < etamax) {
        fFCpt->fillSub(weff, track.pt(), index);
      }
      ++index;
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
      double weff = (dt == kGen) ? 1. : getEfficiency(track, 0);
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
      registryQA.fill(HIST("MCGen/") + HIST(FillTimeName[ft]) + HIST("phi_eta_vtxZ_gen"), track.phi(), track.eta(), vtxz);
      registryQA.fill(HIST("MCGen/") + HIST(FillTimeName[ft]) + HIST("pt_gen"), track.pt());
    } else {
      double wacc = getAcceptance(track, vtxz, 0);
      registryQA.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("phi_eta_vtxZ"), track.phi(), track.eta(), vtxz, (ft == kAfter) ? wacc : 1.0);
      registryQA.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("pt_dcaXY_dcaZ"), track.pt(), track.dcaXY(), track.dcaZ());

      registryQA.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("chi2prTPCcls"), track.tpcChi2NCl());
      registryQA.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("chi2prITScls"), track.itsChi2NCl());
      registryQA.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("nTPCClusters"), track.tpcNClsFound());
      registryQA.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("nITSClusters"), track.itsNCls());
      registryQA.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("nTPCCrossedRows"), track.tpcNClsCrossedRows());

      if (ft == kAfter) {
        registryQA.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("pt_ref"), track.pt());
        registryQA.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("pt_poi"), track.pt());

        if (track.eta() > cfgEtaNch->first && track.eta() < cfgEtaNch->second)
          registryQA.fill(HIST("trackQA/after/etaNch"), track.eta());
        if (track.eta() > cfgEtaPtPt->first && track.eta() < cfgEtaPtPt->second)
          registryQA.fill(HIST("trackQA/after/etaPtPt"), track.eta());
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

  template <QAFillTime ft, typename CollisionObject, typename TracksObject>
  inline void fillEventQA(CollisionObject collision, TracksObject tracks)
  {
    registryQA.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_centT0C"), collision.centFT0C(), tracks.size());
    registryQA.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
    registryQA.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
    registryQA.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_multT0A"), collision.multFT0A(), tracks.size());
    registryQA.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("globalTracks_multV0A"), collision.multFV0A(), tracks.size());
    registryQA.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
    registryQA.fill(HIST("eventQA/") + HIST(FillTimeName[ft]) + HIST("multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
    return;
  }

  o2::framework::expressions::Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZ;
  o2::framework::expressions::Filter trackFilter = (aod::track::eta > cfgEta->first) && (aod::track::eta < cfgEta->second) && (aod::track::pt > cfgPtmin) && (aod::track::pt < cfgPtmax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true)) && (aod::track::itsChi2NCl < cfgTrackCuts.cfgChi2PrITSCls) && (aod::track::tpcChi2NCl < cfgTrackCuts.cfgChi2PrTPCCls) && nabs(aod::track::dcaZ) < cfgTrackCuts.cfgDCAz;

  using GFWCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs>>;
  // using GFWTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFKa, aod::pidTPCKa, aod::pidTOFPr, aod::pidTPCPr>>;
  using GFWTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;

  SliceCache cache;
  Partition<GFWTracks> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<GFWTracks> negTracks = aod::track::signed1Pt < 0.0f;

  double massKaPlus = o2::constants::physics::MassKPlus;
  double massLambda = o2::constants::physics::MassLambda;
  double massK0Short = o2::constants::physics::MassK0Short;

  void processData(GFWCollisions::iterator const& collision, aod::BCsWithTimestamps const&, GFWTracks const& tracks, aod::V0Datas const& v0s)
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
    registryQA.fill(HIST("eventQA/eventSel"), 0.5);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(0.5);
    if (!collision.sel8())
      return;
    registryQA.fill(HIST("eventQA/eventSel"), 1.5);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(1.5);

    float centrality = getCentrality(collision);
    if (cfgDoOccupancySel) {
      int occupancy = collision.trackOccupancyInTimeRange();
      registryQA.fill(HIST("eventQA/before/occ_mult_cent"), occupancy, tracks.size(), centrality);
      if (occupancy < 0 || occupancy > cfgOccupancySelection)
        return;
      registryQA.fill(HIST("eventQA/after/occ_mult_cent"), occupancy, tracks.size(), centrality);
    }
    registryQA.fill(HIST("eventQA/eventSel"), 2.5);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(2.5);
    if (cfgFillQA)
      fillEventQA<kBefore>(collision, tracks);
    registryQA.fill(HIST("eventQA/before/centrality"), centrality);
    registryQA.fill(HIST("eventQA/before/multiplicity"), tracks.size());
    if (!eventSelected(collision, tracks.size(), centrality, run))
      return;
    if (cfgFillQA)
      fillEventQA<kAfter>(collision, tracks);
    registryQA.fill(HIST("eventQA/after/centrality"), centrality);
    registryQA.fill(HIST("eventQA/after/multiplicity"), tracks.size());
    // Get magnetic field polarity
    auto field = (cfgMagField == 99999) ? getMagneticField(bc.timestamp()) : cfgMagField; // o2-linter: disable=magic-number (hard coded default cut)

    processCollision<kReco>(collision, tracks, v0s, centrality, field, run);
  }
  PROCESS_SWITCH(FlowGenericFramework, processData, "Process analysis for non-derived data", true);

  void processMCReco(GFWCollisions::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>> const& tracks, aod::McParticles const&, aod::V0Datas const& v0s)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    if (run != lastRun) {
      lastRun = run;
      if (cfgRunByRun)
        createRunByRunHistograms(run);
    }

    registryQA.fill(HIST("eventQA/eventSel"), 0.5);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(0.5);

    if (!collision.sel8())
      return;

    registryQA.fill(HIST("eventQA/eventSel"), 1.5);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(1.5);

    const auto centrality = getCentrality(collision);

    if (cfgDoOccupancySel) {
      int occupancy = collision.trackOccupancyInTimeRange();
      registryQA.fill(HIST("eventQA/before/occ_mult_cent"), occupancy, tracks.size(), centrality);
      if (occupancy < 0 || occupancy > cfgOccupancySelection)
        return;
      registryQA.fill(HIST("eventQA/after/occ_mult_cent"), occupancy, tracks.size(), centrality);
    }
    registryQA.fill(HIST("eventQA/eventSel"), 2.5);
    if (cfgRunByRun)
      th1sList[run][hEventSel]->Fill(2.5);

    if (cfgFillQA)
      fillEventQA<kBefore>(collision, tracks);
    if (!eventSelected(collision, tracks.size(), centrality, run))
      return;
    if (cfgFillQA)
      fillEventQA<kAfter>(collision, tracks);

    if (!cfgFillWeights)
      loadCorrections(bc);

    auto field = (cfgMagField == 99999) ? getMagneticField(bc.timestamp()) : cfgMagField; // o2-linter: disable=magic-number (hard coded default cut)
    processCollision<kReco>(collision, tracks, v0s, centrality, field, run);
  }
  PROCESS_SWITCH(FlowGenericFramework, processMCReco, "Process analysis for MC reconstructed events", false);

  o2::framework::expressions::Filter mcCollFilter = nabs(aod::mccollision::posZ) < cfgVtxZ;
  void processMCGen(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs>> const& collisions, aod::McParticles const& particles, aod::V0Datas const& v0s)
  {
    if (collisions.size() != 1)
      return;
    float centrality = -1;
    for (const auto& collision : collisions) {
      centrality = getCentrality(collision);
    }
    int run = 0;
    processCollision<kGen>(mcCollision, particles, v0s, centrality, -999, run);
  }
  PROCESS_SWITCH(FlowGenericFramework, processMCGen, "Process analysis for MC generated events", false);

  void processOnTheFly(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, aod::McParticles const& mcParticles, aod::V0Datas const& v0s)
  {
    int run = 0;
    registryQA.fill(HIST("MCGen/impactParameter"), mcCollision.impactParameter(), mcParticles.size());
    processCollision<kGen>(mcCollision, mcParticles, v0s, mcCollision.impactParameter(), -999, run);
  }
  PROCESS_SWITCH(FlowGenericFramework, processOnTheFly, "Process analysis for MC on-the-fly generated events", false);

  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentRun2V0Ms>>::iterator const& collision, aod::BCsWithTimestamps const&, GFWTracks const& tracks, aod::V0Datas const& v0s)
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

    auto field = (cfgMagField == 99999) ? getMagneticField(bc.timestamp()) : cfgMagField; // o2-linter: disable=magic-number (hard coded default cut)
    processCollision<kReco>(collision, tracks, v0s, centrality, field, run);
  }
  PROCESS_SWITCH(FlowGenericFramework, processRun2, "Process analysis for Run 2 converted data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGenericFramework>(cfgc),
  };
}
