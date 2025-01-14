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
///
/// \brief V0 task for production of invariant mass plots for Optimised Topological Cuts Analysis
/// \author Nikolaos Karatzenis (nikolaos.karatzenis@cern.ch)
/// \author Roman Lietava (roman.lietava@cern.ch)

/*Description
This task creates <=20 histograms (for each of the 5 different V0 topological cuts, namely cosPointingAngle,
DCA[between]V0daughters, v0radius,DCA-positive[daughter]to-primary-vertex and DCA-negative[daughter]to-primary-vertex)
that are filled with the V0 invariant mass under the K0, Lambda and Antilambda mass assumption
(so 20cutsx5parametersx3particles=300 mass invariant plots).It also produces plots of the topological parameters themselves.
The cuts are passed as configurable strings for convenience.
This analysis includes two processes, one for Real Data and one for MC Data switchable at the end of the code, only run one at a time*/

#include <memory>
#include <vector>
#include <string>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonUtils/StringUtils.h"

// namespaces to be used for the plot names and topological cuts that will be given by a configurable string
namespace cuthistoskzerosh
{
std::shared_ptr<TH1> cospaCut[20];
static std::vector<std::string> cospacuts;
std::shared_ptr<TH1> dcaCut[20];
static std::vector<std::string> dcacuts;
std::shared_ptr<TH1> v0radiusCut[20];
static std::vector<std::string> v0radiuscuts;
std::shared_ptr<TH1> dcapostopCut[20];
static std::vector<std::string> dcapostopvcuts;
std::shared_ptr<TH1> dcanegtopCut[20];
static std::vector<std::string> dcanegtopvcuts;
} // namespace cuthistoskzerosh
namespace cuthistoslambda
{
std::shared_ptr<TH1> cospaCut[20];
static std::vector<std::string> cospacuts;
std::shared_ptr<TH1> dcaCut[20];
static std::vector<std::string> dcacuts;
std::shared_ptr<TH1> v0radiusCut[20];
static std::vector<std::string> v0radiuscuts;
std::shared_ptr<TH1> dcapostopCut[20];
static std::vector<std::string> dcapostopvcuts;
std::shared_ptr<TH1> dcanegtopCut[20];
static std::vector<std::string> dcanegtopvcuts;
} // namespace cuthistoslambda
namespace cuthistosantilambda
{
std::shared_ptr<TH1> cospaCut[20];
static std::vector<std::string> cospacuts;
std::shared_ptr<TH1> dcaCut[20];
static std::vector<std::string> dcacuts;
std::shared_ptr<TH1> v0radiusCut[20];
static std::vector<std::string> v0radiuscuts;
std::shared_ptr<TH1> dcapostopCut[20];
static std::vector<std::string> dcapostopvcuts;
std::shared_ptr<TH1> dcanegtopCut[20];
static std::vector<std::string> dcanegtopvcuts;
} // namespace cuthistosantilambda
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct v0topologicalcuts {
  // Histogram Registry includes different V0 Parameteres for all V0s and individual MC-V0s with MC-matching
  HistogramRegistry rV0Parameters_MC_V0match{"V0Parameters_MC_V0Match", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rV0Parameters_MC_K0Smatch{"V0Parameters_MC_K0SMatch", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rV0Parameters_MC_Lambdamatch{"V0Parameters_MC_LambdaMatch", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rV0Parameters_MC_AntiLambdamatch{"V0Parameters_MC_AntiLambdaMatch", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rV0Parameters_Data{"rV0Parameters_Data", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // kzero cut Histogram Registry with MC-matching, each will include 20 histograms for 20 different cuts
  HistogramRegistry rKzeroShort_cospaCut{"KzeroShort_cospaCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShort_dcaCut{"KzeroShort_dcaCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShort_v0radiusCut{"KzeroShort_v0radiusCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShort_dcapostopCut{"KzeroShort_dcapostopvCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShort_dcanegtopCut{"KzeroShort_dcanegtopvCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // lambdas cut histograms with MC-matching (same as in Kzeros above)
  HistogramRegistry rLambda_cospaCut{"Lambda_cospaCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda_dcaCut{"Lambda_dcaCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda_v0radiusCut{"Lambda_v0radiusCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda_dcapostopCut{"Lambda_dcapostopvCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda_dcanegtopCut{"Lambda_dcanegtopvCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // antilambdas cut histograms with MC-matching (same as in Lambdas an Kzeros above)
  HistogramRegistry rAntiLambda_cospaCut{"AntiLambda_cospaCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntiLambda_dcaCut{"AntiLambda_dcaCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntiLambda_v0radiusCut{"AntiLambda_v0radiusCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntiLambda_dcapostopCut{"AntiLambda_dcapostopvCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntiLambda_dcanegtopCut{"AntiLambda_dcanegtopvCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurables for Cuts
  Configurable<float> cutZVertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> nSigmaTPCPion{"nSigmaTPCPion", 4, "nSigmaTPCPion"};
  Configurable<float> nSigmaTPCProton{"nSigmaTOCProton", 4, "nSigmaTPCProton"};
  Configurable<float> compv0masscut{"compv0masscut", 0.01, "CompetitiveV0masscut (GeV)"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};

  // Configurable strings for Kzero cuts
  Configurable<std::string> kzeroshsetting_cospacuts_string{"kzerosetting_cospacuts", {"0_98,0_981,0_982,0_983,0_984,0_985,0_986,0_987,0_988,0_989,0_99,0_991,0_992,0_993,0_994,0_995,0_996,0_997,0_998,0_999"}, "Kzero cosPA Cut Values"};
  Configurable<std::string> kzeroshsetting_dcacuts_string{"kzerosetting_dcacuts", {"0_3,0_285,0_27,0_255,0_24,0_225,0_21,0_195,0_18,0_165,0_15,0_135,0_12,0_105,0_09,0_075,0_06,0_045,0_03,0_015"}, "Kzero DCA Cut Values"};
  Configurable<std::string> kzeroshsetting_v0radius_string{"kzerosetting_v0radiuscuts", {"0_5,0_51,0_52,0_53,0_54,0_55,0_56,0_57,0_58,0_59,0_6,0_61,0_62,0_63,0_64,0_65,0_66,0_67,0_68,0_69"}, "Kzero V0Radius Cut Values"};
  Configurable<std::string> kzeroshsetting_dcapostopv_string{"kzerosetting_dcapostopvcuts", {"0_0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Kzero DCA Pos to PV Cut Values"};
  Configurable<std::string> kzeroshsetting_dcanegtopv_string{"kzerosetting_dcanegtopvcuts", {"0_0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "KzeroDCA Neg to PV Cut Values"};

  // Configurable strings for Lambdacuts
  Configurable<std::string> lambdasetting_cospacuts_string{"lambdasetting_cospacuts", {"0_98,0_981,0_982,0_983,0_984,0_985,0_986,0_987,0_988,0_989,0_99,0_991,0_992,0_993,0_994"}, "Lambda cosPA Cut Values"};
  Configurable<std::string> lambdasetting_dcacuts_string{"lambdasetting_dcacuts", {"0_3,0_285,0_27,0_255,0_24,0_225,0_21,0_195,0_18,0_165,0_15,0_135,0_12,0_105,0_09,0_075,0_06,0_045,0_03,0_015"}, "Lambda DCA Cut Values"};
  Configurable<std::string> lambdasetting_v0radius_string{"lambdasetting_v0radiuscuts", {"0_5,0_51,0_52,0_53,0_54,0_55,0_56,0_57,0_58,0_59,0_6,0_61,0_62,0_63,0_64,0_65,0_66,0_67,0_68,0_69"}, "Lambda V0Radius Cut Values"};
  Configurable<std::string> lambdasetting_dcapostopv_string{"lambdasetting_dcapostopvcuts", {"0_0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Lambda DCA Pos to PV Cut Values"};
  Configurable<std::string> lambdasetting_dcanegtopv_string{"lambdasetting_dcanegtopvcuts", {"0_0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Lambda DCA Neg to PV Cut Values"};

  // Configurable strings for AntiLambdacuts
  Configurable<std::string> antilambdasetting_cospacuts_string{"antilambdasetting_cospacuts", {"0_98,0_981,0_982,0_983,0_984,0_985,0_986,0_987,0_988,0_989,0_99,0_991,0_992,0_993,0_994,0_995,0_996,0_997,0_998,0_999"}, "Antilambda cosPA Cut Values"};
  Configurable<std::string> antilambdasetting_dcacuts_string{"antilambdasetting_dcacuts", {"0_3,0_285,0_27,0_255,0_24,0_225,0_21,0_195,0_18,0_165,0_15,0_135,0_12,0_105,0_09,0_075,0_06,0_045,0_03,0_015"}, "Antilambda DCA Cut Values"};
  Configurable<std::string> antilambdasetting_v0radius_string{"antilambdasetting_v0radiuscuts", {"0_5,0_51,0_52,0_53,0_54,0_55,0_56,0_57,0_58,0_59,0_6,0_61,0_62,0_63,0_64,0_65,0_66,0_67,0_68,0_69"}, "Antilambda V0Radius Cut Values"};
  Configurable<std::string> antilambdasetting_dcapostopv_string{"antilambdasetting_dcapostopvcuts", {"0_0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Antilambda DCA Pos to PV Cut Values"};
  Configurable<std::string> antilambdasetting_dcanegtopv_string{"antilambdasetting_dcanegtopvcuts", {"0_0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Antilambda DCA Neg to PV Cut Values"};

  void init(InitContext const&)
  {
    // kzero filling namespace with configurable strings

    // setting strings from configurable strings in order to manipulate them
    // getting the  cut values for the names of the plots for the five topological cuts
    cuthistoskzerosh::cospacuts = o2::utils::Str::tokenize(kzeroshsetting_cospacuts_string, ',');
    cuthistoskzerosh::dcacuts = o2::utils::Str::tokenize(kzeroshsetting_dcacuts_string, ',');
    cuthistoskzerosh::v0radiuscuts = o2::utils::Str::tokenize(kzeroshsetting_v0radius_string, ',');
    cuthistoskzerosh::dcapostopvcuts = o2::utils::Str::tokenize(kzeroshsetting_dcapostopv_string, ',');
    cuthistoskzerosh::dcanegtopvcuts = o2::utils::Str::tokenize(kzeroshsetting_dcanegtopv_string, ',');

    // lambda filling namespace with configurable strings (same as in Kzeros above)
    cuthistoslambda::cospacuts = o2::utils::Str::tokenize(lambdasetting_cospacuts_string, ',');
    cuthistoslambda::dcacuts = o2::utils::Str::tokenize(lambdasetting_dcacuts_string, ',');
    cuthistoslambda::v0radiuscuts = o2::utils::Str::tokenize(lambdasetting_v0radius_string, ',');
    cuthistoslambda::dcapostopvcuts = o2::utils::Str::tokenize(lambdasetting_dcapostopv_string, ',');
    cuthistoslambda::dcanegtopvcuts = o2::utils::Str::tokenize(lambdasetting_dcanegtopv_string, ',');

    // antilambda filling namespace with configurable strings (same as in Lambdas and Kzeros above)
    cuthistosantilambda::cospacuts = o2::utils::Str::tokenize(antilambdasetting_cospacuts_string, ',');
    cuthistosantilambda::dcacuts = o2::utils::Str::tokenize(antilambdasetting_dcacuts_string, ',');
    cuthistosantilambda::v0radiuscuts = o2::utils::Str::tokenize(antilambdasetting_v0radius_string, ',');
    cuthistosantilambda::dcapostopvcuts = o2::utils::Str::tokenize(antilambdasetting_dcapostopv_string, ',');
    cuthistosantilambda::dcanegtopvcuts = o2::utils::Str::tokenize(antilambdasetting_dcanegtopv_string, ',');

    // Axes for the three invariant mass plots
    AxisSpec K0ShortMassAxis = {nBins, 0.45f, 0.55f, "#it{M} #pi^{+}#pi^{-} [GeV/#it{c}^{2}]"};
    AxisSpec LambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M} p^{+}#pi^{-} [GeV/#it{c}^{2}]"};
    AxisSpec AntiLambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M} p^{-}#pi^{+} [GeV/#it{c}^{2}]"};
    AxisSpec ptAxis = {nBins, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};

    // adding the invariant mass histograms to their Registries using the namespace for kzeros, lambdas and antilambdas
    for (uint32_t i = 0; i < cuthistoskzerosh::cospacuts.size(); i++) {
      cuthistoskzerosh::cospaCut[i] = rKzeroShort_cospaCut.add<TH1>(fmt::format("hKzerocospaCut_{}", cuthistoskzerosh::cospacuts[i]).data(), fmt::format("hKzerocospaCut_{}", cuthistoskzerosh::cospacuts[i]).data(), {HistType::kTH1D, {{K0ShortMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoskzerosh::dcacuts.size(); i++) {
      cuthistoskzerosh::dcaCut[i] = rKzeroShort_dcaCut.add<TH1>(fmt::format("hKzerodcaCut_{}", cuthistoskzerosh::dcacuts[i]).data(), fmt::format("hKzerodcaCut_{}", cuthistoskzerosh::dcacuts[i]).data(), {HistType::kTH1D, {{K0ShortMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoskzerosh::v0radiuscuts.size(); i++) {
      cuthistoskzerosh::v0radiusCut[i] = rKzeroShort_v0radiusCut.add<TH1>(fmt::format("hKzerov0radiusCut_{}", cuthistoskzerosh::v0radiuscuts[i]).data(), fmt::format("hKzerov0radiusCut_{}", cuthistoskzerosh::v0radiuscuts[i]).data(), {HistType::kTH1D, {{K0ShortMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoskzerosh::dcapostopvcuts.size(); i++) {
      cuthistoskzerosh::dcapostopCut[i] = rKzeroShort_dcapostopCut.add<TH1>(fmt::format("hKzerodcapostopCut_{}", cuthistoskzerosh::dcapostopvcuts[i]).data(), fmt::format("hKzerodcapostopCut_{}", cuthistoskzerosh::dcapostopvcuts[i]).data(), {HistType::kTH1D, {{K0ShortMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoskzerosh::dcanegtopvcuts.size(); i++) {
      cuthistoskzerosh::dcanegtopCut[i] = rKzeroShort_dcanegtopCut.add<TH1>(fmt::format("hKzerodcanegtopCut_{}", cuthistoskzerosh::dcanegtopvcuts[i]).data(), fmt::format("hKzerodcanegtopCut_{}", cuthistoskzerosh::dcanegtopvcuts[i]).data(), {HistType::kTH1D, {{K0ShortMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoslambda::cospacuts.size(); i++) {
      cuthistoslambda::cospaCut[i] = rLambda_cospaCut.add<TH1>(fmt::format("hLambdacospaCut_{}", cuthistoslambda::cospacuts[i]).data(), fmt::format("hLambdacospaCut_{}", cuthistoslambda::cospacuts[i]).data(), {HistType::kTH1D, {{LambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoslambda::dcacuts.size(); i++) {
      cuthistoslambda::dcaCut[i] = rLambda_dcaCut.add<TH1>(fmt::format("hLambdadcaCut_{}", cuthistoslambda::dcacuts[i]).data(), fmt::format("hLambdadcaCut_{}", cuthistoslambda::dcacuts[i]).data(), {HistType::kTH1D, {{LambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoslambda::v0radiuscuts.size(); i++) {
      cuthistoslambda::v0radiusCut[i] = rLambda_v0radiusCut.add<TH1>(fmt::format("hLambdav0radiusCut_{}", cuthistoslambda::v0radiuscuts[i]).data(), fmt::format("hLambdav0radiusCut_{}", cuthistoslambda::v0radiuscuts[i]).data(), {HistType::kTH1D, {{LambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoslambda::dcapostopvcuts.size(); i++) {
      cuthistoslambda::dcapostopCut[i] = rLambda_dcapostopCut.add<TH1>(fmt::format("hLambdadcapostopCut_{}", cuthistoslambda::dcapostopvcuts[i]).data(), fmt::format("hLambdadcapostopCut_{}", cuthistoslambda::dcapostopvcuts[i]).data(), {HistType::kTH1D, {{LambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoslambda::dcanegtopvcuts.size(); i++) {
      cuthistoslambda::dcanegtopCut[i] = rLambda_dcanegtopCut.add<TH1>(fmt::format("hLambdadcanegtopCut_{}", cuthistoslambda::dcanegtopvcuts[i]).data(), fmt::format("hLambdadcanegtopCut_{}", cuthistoslambda::dcanegtopvcuts[i]).data(), {HistType::kTH1D, {{LambdaMassAxis}}});
    }

    for (uint32_t i = 0; i < cuthistosantilambda::cospacuts.size(); i++) {
      cuthistosantilambda::cospaCut[i] = rAntiLambda_cospaCut.add<TH1>(fmt::format("hAntiLambdacospaCut_{}", cuthistosantilambda::cospacuts[i]).data(), fmt::format("hAntiLambdacospaCut_{}", cuthistosantilambda::cospacuts[i]).data(), {HistType::kTH1D, {{AntiLambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistosantilambda::dcacuts.size(); i++) {
      cuthistosantilambda::dcaCut[i] = rAntiLambda_dcaCut.add<TH1>(fmt::format("hAntiLambdadcaCut_{}", cuthistosantilambda::dcacuts[i]).data(), fmt::format("hAntiLambdadcaCut_{}", cuthistosantilambda::dcacuts[i]).data(), {HistType::kTH1D, {{AntiLambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistosantilambda::v0radiuscuts.size(); i++) {
      cuthistosantilambda::v0radiusCut[i] = rAntiLambda_v0radiusCut.add<TH1>(fmt::format("hAntiLambdav0radiusCut_{}", cuthistosantilambda::v0radiuscuts[i]).data(), fmt::format("hAntiLambdav0radiusCut_{}", cuthistosantilambda::v0radiuscuts[i]).data(), {HistType::kTH1D, {{AntiLambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistosantilambda::dcapostopvcuts.size(); i++) {
      cuthistosantilambda::dcapostopCut[i] = rAntiLambda_dcapostopCut.add<TH1>(fmt::format("hAntiLambdadcapostopCut_{}", cuthistosantilambda::dcapostopvcuts[i]).data(), fmt::format("hAntiLambdadcapostopCut_{}", cuthistosantilambda::dcapostopvcuts[i]).data(), {HistType::kTH1D, {{AntiLambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistosantilambda::dcanegtopvcuts.size(); i++) {
      cuthistosantilambda::dcanegtopCut[i] = rAntiLambda_dcanegtopCut.add<TH1>(fmt::format("hAntiLambdadcanegtopCut_{}", cuthistosantilambda::dcanegtopvcuts[i]).data(), fmt::format("hAntiLambdadcanegtopCut_{}", cuthistosantilambda::dcanegtopvcuts[i]).data(), {HistType::kTH1D, {{AntiLambdaMassAxis}}});
    }

    // K0s topological cut histograms added and MC-matched
    rV0Parameters_MC_V0match.add("hDCAV0Daughters_V0_Match", "hDCAV0Daughters_No_Match", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0Parameters_MC_V0match.add("hV0CosPA_V0_Match", "hV0CosPA_No_Match", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0Parameters_MC_V0match.add("hV0Radius_V0_Match", "hV0Radius_No_Match", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});
    rV0Parameters_MC_V0match.add("hV0Radius_V0_Match_Full", "hV0Radius_No_Match_Full", {HistType::kTH1F, {{nBins, 0.0f, 40.0f}}});
    rV0Parameters_MC_V0match.add("hDCAPostoPV_V0_Match", "hDCAPostoPV_No_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0Parameters_MC_V0match.add("hDCANegtoPV_V0_Match", "hDCANegtoPV_No_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0Parameters_MC_V0match.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // kzero match
    rV0Parameters_MC_K0Smatch.add("hDCAV0Daughters_KzeroMC_Match", "hDCAV0Daughters_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0Parameters_MC_K0Smatch.add("hV0CosPA_KzeroMC_Match", "hV0CosPA_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0Parameters_MC_K0Smatch.add("hV0Radius_KzeroMC_Match", "hV0Radius_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.2f, 5.0f}}});
    rV0Parameters_MC_K0Smatch.add("hDCAPostoPV_KzeroMC_Match", "hDCAPostoPV_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0Parameters_MC_K0Smatch.add("hDCANegtoPV_KzeroMC_Match", "hDCANegtoPV_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});

    // lambda match
    rV0Parameters_MC_Lambdamatch.add("hDCAV0Daughters_LambdaMC_Match", "hDCAV0Daughters_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0Parameters_MC_Lambdamatch.add("hV0CosPA_LambdaMC_Match", "hV0CosPA_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0Parameters_MC_Lambdamatch.add("hV0Radius_LambdaMC_Match", "hV0Radius_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.2f, 5.0f}}});
    rV0Parameters_MC_Lambdamatch.add("hDCAPostoPV_LambdaMC_Match", "hDCAPostoPV_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0Parameters_MC_Lambdamatch.add("hDCANegtoPV_LambdaMC_Match", "hDCANegtoPV_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});

    // antilambda match
    rV0Parameters_MC_AntiLambdamatch.add("hDCAV0Daughters_AntiLambdaMC_Match", "hDCAV0Daughters_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0Parameters_MC_AntiLambdamatch.add("hV0CosPA_AntiLambdaMC_Match", "hV0CosPA_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0Parameters_MC_AntiLambdamatch.add("hV0Radius_AntiLambdaMC_Match", "hV0Radius_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.2f, 5.0f}}});
    rV0Parameters_MC_AntiLambdamatch.add("hDCAPostoPV_AntiLambdaMC_Match", "hDCAPostoPV_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0Parameters_MC_AntiLambdamatch.add("hDCANegtoPV_AntiLambdaMC_Match", "hDCANegtoPV_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});

    // V0s Data
    rV0Parameters_Data.add("hDCAV0Daughters_V0_Data", "hDCAV0Daughters_V0_Data", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0Parameters_Data.add("hV0CosPA_V0_Data", "hV0CosPA_V0_Data", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0Parameters_Data.add("hV0Radius_V0_Data", "hV0Radius_V0_Data", {HistType::kTH1F, {{nBins, 0.2f, 5.0f}}});
    rV0Parameters_Data.add("hV0Radius_Full_V0_Data", "hV0Radius_Full_V0_Data", {HistType::kTH1F, {{nBins, 0.2f, 40.0f}}});
    rV0Parameters_Data.add("hDCAPostoPV_V0_Data", "hDCAPostoPV_V0_Data", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0Parameters_Data.add("hDCANegtoPV_V0_Data", "hDCANegtoPV_V0_Data", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0Parameters_Data.add("hMassK0ShortNoCuts_V0_Data", "hMassK0ShortNoCuts_V0_Data", {HistType::kTH1F, {{K0ShortMassAxis}}});
    rV0Parameters_Data.add("hMassLambdaNoCuts_V0_Data", "hMassLambdaNoCuts_V0_Data", {HistType::kTH1F, {{LambdaMassAxis}}});
    rV0Parameters_Data.add("hMassAntilambdaNoCuts_V0_Data", "hMassAntilambdaNoCuts_V0_Data", {HistType::kTH1F, {{AntiLambdaMassAxis}}});
    rV0Parameters_Data.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rV0Parameters_Data.add("hMassK0ShortAllAfterCutsandTPCPID", "hMassK0ShortAllAfterCutsandTPCPID", {HistType::kTH1F, {K0ShortMassAxis}});
    rV0Parameters_Data.add("hMassK0ShortAllAfterCutTPCPIDCompV0MassCut", "hMassK0ShortAllAfterCutTPCPIDCompV0MassCut", {HistType::kTH1F, {K0ShortMassAxis}});
    rV0Parameters_Data.add("hNSigmaPosPionFromK0s", "hNSigmaPosPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rV0Parameters_Data.add("hNSigmaNegPionFromK0s", "hNSigmaNegPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rV0Parameters_Data.add("hMassLambdaAllAfterCutsandTPCPID", "hMassLambdaAllAfterCutsandTPCPID", {HistType::kTH1F, {LambdaMassAxis}});
    rV0Parameters_Data.add("hNSigmaPosProtonFromLambda", "hNSigmaPosProtonFromLambda", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rV0Parameters_Data.add("hNSigmaNegPionFromLambda", "hNSigmaNegPionFromLambda", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rV0Parameters_Data.add("hMassAntilambdaAllAfterCutsandTPCPID", "hMassAntilambdaAllAfterCutsandTPCPID", {HistType::kTH1F, {LambdaMassAxis}});
    rV0Parameters_Data.add("hNSigmaNegProtonFromAntilambda", "hNSigmaNegProtonFromAntilambda", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rV0Parameters_Data.add("hNSigmaPosPionFromAntiambda", "hNSigmaPosPionFromAntilambda", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < cutzvertex);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;

  // This is the Process for the MC reconstructed Data
  void RecMCprocess(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                    soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s,
                    DaughterTracks const&, // no need to define a variable for tracks, if we don't access them directly
                    aod::McParticles const&)
  {
    const auto& mLambdaPDG = 1.115683;
    const auto& mK0shPDG = 0.497611;
    for (const auto& v0 : V0s) {
      if (std::abs(v0.posTrack_as<DaughterTracks>().eta()) < etadau && std::abs(v0.negTrack_as<DaughterTracks>().eta()) < etadau) { // daughters pseudorapidity cut
        // filling histograms with V0 values
        rV0Parameters_MC_V0match.fill(HIST("hDCAV0Daughters_V0_Match"), v0.dcaV0daughters());
        rV0Parameters_MC_V0match.fill(HIST("hV0CosPA_V0_Match"), v0.v0cosPA());
        rV0Parameters_MC_V0match.fill(HIST("hV0Radius_V0_Match"), v0.v0radius());
        rV0Parameters_MC_V0match.fill(HIST("hV0Radius_V0_Match_Full"), v0.v0radius());
        rV0Parameters_MC_V0match.fill(HIST("hDCAPostoPV_V0_Match"), v0.dcapostopv());
        rV0Parameters_MC_V0match.fill(HIST("hDCANegtoPV_V0_Match"), v0.dcanegtopv());
        rV0Parameters_MC_V0match.fill(HIST("hVertexZRec"), collision.posZ());

        // Checking that the V0 is a true K0s/Lambdas/Antilambdas and then filling the parameter histograms and the invariant mass plots for different cuts (which are taken from namespace)
        if (v0.has_mcParticle()) {
          auto v0mcParticle = v0.mcParticle();
          if (v0mcParticle.pdgCode() == 310) {                                                                                    // kzero matched
            if (std::abs(v0.mLambda() - mLambdaPDG) > compv0masscut && std::abs(v0.mAntiLambda() - mLambdaPDG) > compv0masscut) { // Kzero competitive v0 mass cut (cut out Lambdas and Anti-Lambdas)
              rV0Parameters_MC_K0Smatch.fill(HIST("hDCAV0Daughters_KzeroMC_Match"), v0.dcaV0daughters());
              rV0Parameters_MC_K0Smatch.fill(HIST("hV0CosPA_KzeroMC_Match"), v0.v0cosPA());
              rV0Parameters_MC_K0Smatch.fill(HIST("hV0Radius_KzeroMC_Match"), v0.v0radius());
              rV0Parameters_MC_K0Smatch.fill(HIST("hDCAPostoPV_KzeroMC_Match"), std::abs(v0.dcapostopv()));
              rV0Parameters_MC_K0Smatch.fill(HIST("hDCANegtoPV_KzeroMC_Match"), std::abs(v0.dcanegtopv()));

              for (uint32_t j = 0; j < cuthistoskzerosh::cospacuts.size(); j++) {
                std::string cospacut = cuthistoskzerosh::cospacuts[j]; // Get the current cut value from the namespace
                size_t pos = cospacut.find("_");                       // find the "_" which needs to change to a "." for it to be a number
                cospacut[pos] = '.';                                   // change the "_" into an "."
                const float cospacutvalue = std::stod(cospacut);       // make the string into a float value
                if (v0.v0cosPA() > cospacutvalue) {                    // enforce the cut value
                  cuthistoskzerosh::cospaCut[j]->Fill(v0.mK0Short());  // fill the corresponding histo from the namespace with the invariant mass (of a Kzero here)
                }
              }
              for (uint32_t j = 0; j < cuthistoskzerosh::dcacuts.size(); j++) {
                std::string dcacut = cuthistoskzerosh::dcacuts[j];
                size_t pos = dcacut.find("_");
                dcacut[pos] = '.';
                const float dcacutvalue = std::stod(dcacut);
                if (v0.dcaV0daughters() < dcacutvalue) {
                  cuthistoskzerosh::dcaCut[j]->Fill(v0.mK0Short());
                }
              }
              for (uint32_t j = 0; j < cuthistoskzerosh::v0radiuscuts.size(); j++) {
                std::string v0radiuscut = cuthistoskzerosh::v0radiuscuts[j];
                size_t pos = v0radiuscut.find("_");
                v0radiuscut[pos] = '.';
                const float v0radiuscutvalue = std::stod(v0radiuscut);
                if (v0.v0radius() > v0radiuscutvalue) {
                  cuthistoskzerosh::v0radiusCut[j]->Fill(v0.mK0Short());
                }
              }
              for (uint32_t j = 0; j < cuthistoskzerosh::dcapostopvcuts.size(); j++) {
                std::string dcapostopcut = cuthistoskzerosh::dcapostopvcuts[j];
                size_t pos = dcapostopcut.find("_");
                dcapostopcut[pos] = '.';
                const float dcapostopcutvalue = std::stod(dcapostopcut);
                if (std::abs(v0.dcapostopv()) > dcapostopcutvalue) {
                  cuthistoskzerosh::dcapostopCut[j]->Fill(v0.mK0Short());
                }
              }
              for (uint32_t j = 0; j < cuthistoskzerosh::dcanegtopvcuts.size(); j++) {
                std::string dcanegtopcut = cuthistoskzerosh::dcanegtopvcuts[j];
                size_t pos = dcanegtopcut.find("_");
                dcanegtopcut[pos] = '.';
                const float dcanegtopcutvalue = std::stod(dcanegtopcut);
                if (std::abs(v0.dcanegtopv()) > dcanegtopcutvalue) {
                  cuthistoskzerosh::dcanegtopCut[j]->Fill(v0.mK0Short());
                }
              }
            }
          }
          if (v0mcParticle.pdgCode() == 3122) {                       // lambda matched
            if (std::abs(v0.mK0Short() - mK0shPDG) > compv0masscut) { // antilambda competitive v0 mass cut (cut out Kaons)
              rV0Parameters_MC_Lambdamatch.fill(HIST("hDCAV0Daughters_LambdaMC_Match"), v0.dcaV0daughters());
              rV0Parameters_MC_Lambdamatch.fill(HIST("hV0CosPA_LambdaMC_Match"), v0.v0cosPA());
              rV0Parameters_MC_Lambdamatch.fill(HIST("hV0Radius_LambdaMC_Match"), v0.v0radius());
              rV0Parameters_MC_Lambdamatch.fill(HIST("hDCAPostoPV_LambdaMC_Match"), std::abs(v0.dcapostopv()));
              rV0Parameters_MC_Lambdamatch.fill(HIST("hDCANegtoPV_LambdaMC_Match"), std::abs(v0.dcanegtopv()));

              // for explanation look at the first Kzero  plot above
              for (uint32_t j = 0; j < cuthistoslambda::cospacuts.size(); j++) {
                std::string cospacutlambda = cuthistoslambda::cospacuts[j];
                size_t pos = cospacutlambda.find("_");
                cospacutlambda[pos] = '.';
                const float cospacutlambdavalue = std::stod(cospacutlambda);
                if (v0.v0cosPA() > cospacutlambdavalue) {
                  cuthistoslambda::cospaCut[j]->Fill(v0.mLambda());
                }
              }
              for (uint32_t j = 0; j < cuthistoslambda::dcacuts.size(); j++) {
                std::string dcacutlambda = cuthistoslambda::dcacuts[j];
                size_t pos = dcacutlambda.find("_");
                dcacutlambda[pos] = '.';
                const float dcacutlambdavalue = std::stod(dcacutlambda);
                if (v0.dcaV0daughters() < dcacutlambdavalue) {
                  cuthistoslambda::dcaCut[j]->Fill(v0.mLambda());
                }
              }
              for (uint32_t j = 0; j < cuthistoslambda::v0radiuscuts.size(); j++) {
                std::string v0radiuscutlambda = cuthistoslambda::v0radiuscuts[j];
                size_t pos = v0radiuscutlambda.find("_");
                v0radiuscutlambda[pos] = '.';
                const float v0radiuscutlambdavalue = std::stod(v0radiuscutlambda);
                if (v0.v0radius() > v0radiuscutlambdavalue) {
                  cuthistoslambda::v0radiusCut[j]->Fill(v0.mLambda());
                }
              }
              for (uint32_t j = 0; j < cuthistoslambda::dcapostopvcuts.size(); j++) {
                std::string dcapostopcutlambda = cuthistoslambda::dcapostopvcuts[j];
                size_t pos = dcapostopcutlambda.find("_");
                dcapostopcutlambda[pos] = '.';
                const float dcapostopcutlambdavalue = std::stod(dcapostopcutlambda);
                if (std::abs(v0.dcapostopv()) > dcapostopcutlambdavalue) {
                  cuthistoslambda::dcapostopCut[j]->Fill(v0.mLambda());
                }
              }
              for (uint32_t j = 0; j < cuthistoslambda::dcanegtopvcuts.size(); j++) {
                std::string dcanegtopcutlambda = cuthistoslambda::dcanegtopvcuts[j];
                size_t pos = dcanegtopcutlambda.find("_");
                dcanegtopcutlambda[pos] = '.';
                const float dcanegtopcutlambdavalue = std::stod(dcanegtopcutlambda);
                if (std::abs(v0.dcanegtopv()) > dcanegtopcutlambdavalue) {
                  cuthistoslambda::dcanegtopCut[j]->Fill(v0.mLambda());
                }
              }
            }
          }
          if (v0mcParticle.pdgCode() == -3122) {                      // antilambda matched
            if (std::abs(v0.mK0Short() - mK0shPDG) > compv0masscut) { // antilambda competitive v0 mass cut (cut out Kaons)
              rV0Parameters_MC_AntiLambdamatch.fill(HIST("hDCAV0Daughters_AntiLambdaMC_Match"), v0.dcaV0daughters());
              rV0Parameters_MC_AntiLambdamatch.fill(HIST("hV0CosPA_AntiLambdaMC_Match"), v0.v0cosPA());
              rV0Parameters_MC_AntiLambdamatch.fill(HIST("hV0Radius_AntiLambdaMC_Match"), v0.v0radius());
              rV0Parameters_MC_AntiLambdamatch.fill(HIST("hDCAPostoPV_AntiLambdaMC_Match"), std::abs(v0.dcapostopv()));
              rV0Parameters_MC_AntiLambdamatch.fill(HIST("hDCANegtoPV_AntiLambdaMC_Match"), std::abs(v0.dcanegtopv()));
              // for explanation look at the first Kzero  plot above
              for (uint32_t j = 0; j < cuthistosantilambda::cospacuts.size(); j++) {
                std::string cospacutantilambda = cuthistosantilambda::cospacuts[j];
                size_t pos = cospacutantilambda.find("_");
                cospacutantilambda[pos] = '.';
                const float cospacutantilambdavalue = std::stod(cospacutantilambda);
                if (v0.v0cosPA() > cospacutantilambdavalue) {
                  cuthistosantilambda::cospaCut[j]->Fill(v0.mAntiLambda());
                }
              }
              for (uint32_t j = 0; j < cuthistosantilambda::dcacuts.size(); j++) {
                std::string dcacutantilambda = cuthistosantilambda::dcacuts[j];
                size_t pos = dcacutantilambda.find("_");
                dcacutantilambda[pos] = '.';
                const float dcacutantilambdavalue = std::stod(dcacutantilambda);
                if (v0.dcaV0daughters() < dcacutantilambdavalue) {
                  cuthistosantilambda::dcaCut[j]->Fill(v0.mAntiLambda());
                }
              }
              for (uint32_t j = 0; j < cuthistosantilambda::v0radiuscuts.size(); j++) {
                std::string v0radiusantilambda = cuthistosantilambda::v0radiuscuts[j];
                size_t pos = v0radiusantilambda.find("_");
                v0radiusantilambda[pos] = '.';
                const float v0radiuscutantilambdavalue = std::stod(v0radiusantilambda);
                if (v0.v0radius() > v0radiuscutantilambdavalue) {
                  cuthistosantilambda::v0radiusCut[j]->Fill(v0.mAntiLambda());
                }
              }
              for (uint32_t j = 0; j < cuthistosantilambda::dcapostopvcuts.size(); j++) {
                std::string dcapostopantilambda = cuthistosantilambda::dcapostopvcuts[j];
                size_t pos = dcapostopantilambda.find("_");
                dcapostopantilambda[pos] = '.';
                const float dcapostopcutantilambdavalue = std::stod(dcapostopantilambda);
                if (std::abs(v0.dcapostopv()) > dcapostopcutantilambdavalue) {
                  cuthistosantilambda::dcapostopCut[j]->Fill(v0.mAntiLambda());
                }
              }
              for (uint32_t j = 0; j < cuthistosantilambda::dcanegtopvcuts.size(); j++) {
                std::string dcanegtopantilambda = cuthistosantilambda::dcanegtopvcuts[j];
                size_t pos = dcanegtopantilambda.find("_");
                dcanegtopantilambda[pos] = '.';
                const float dcanegtopcutantilambdavalue = std::stod(dcanegtopantilambda);
                if (std::abs(v0.dcanegtopv()) > dcanegtopcutantilambdavalue) {
                  cuthistosantilambda::dcanegtopCut[j]->Fill(v0.mAntiLambda());
                }
              }
            }
          }
        }
      }
    }
  }
  // This is the process for Real Data
  void Dataprocess(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                   aod::V0Datas const& V0s,
                   DaughterTracks const&)
  {
    const auto& mLambdaPDG = 1.115683;
    const auto& mK0shPDG = 0.497611;
    // filling histograms with the different V0 parameters
    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<DaughterTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<DaughterTracks>();
      rV0Parameters_Data.fill(HIST("hMassK0ShortNoCuts_V0_Data"), v0.mK0Short());
      rV0Parameters_Data.fill(HIST("hMassLambdaNoCuts_V0_Data"), v0.mLambda());
      rV0Parameters_Data.fill(HIST("hMassAntilambdaNoCuts_V0_Data"), v0.mAntiLambda());
      rV0Parameters_Data.fill(HIST("hDCAV0Daughters_V0_Data"), v0.dcaV0daughters());
      rV0Parameters_Data.fill(HIST("hV0CosPA_V0_Data"), v0.v0cosPA());
      rV0Parameters_Data.fill(HIST("hV0Radius_V0_Data"), v0.v0radius());
      rV0Parameters_Data.fill(HIST("hV0Radius_Full_V0_Data"), v0.v0radius());
      rV0Parameters_Data.fill(HIST("hDCAPostoPV_V0_Data"), std::abs(v0.dcapostopv()));
      rV0Parameters_Data.fill(HIST("hDCANegtoPV_V0_Data"), std::abs(v0.dcanegtopv()));
      rV0Parameters_Data.fill(HIST("hVertexZRec"), collision.posZ());

      if (std::abs(v0.posTrack_as<DaughterTracks>().eta()) < etadau && std::abs(v0.negTrack_as<DaughterTracks>().eta()) < etadau) { // daughters pseudorapidity cut
        if (std::abs(posDaughterTrack.tpcNSigmaPi()) < nSigmaTPCPion && std::abs(negDaughterTrack.tpcNSigmaPi()) < nSigmaTPCPion) { // TPC PID on daughter pions
          rV0Parameters_Data.fill(HIST("hMassK0ShortAllAfterCutsandTPCPID"), v0.mK0Short());
          if (std::abs(v0.mLambda() - mLambdaPDG) > compv0masscut && std::abs(v0.mAntiLambda() - mLambdaPDG) > compv0masscut) { // antilambda competitive v0 mass cut (cut out Lambdas and Anti-Lambdas)
            rV0Parameters_Data.fill(HIST("hMassK0ShortAllAfterCutTPCPIDCompV0MassCut"), v0.mK0Short());
            // Filling the five Kzero invariant mass plots for different cuts (which are taken from namespace), for full explanation see the first kzero cut filling in the MC process
            for (uint32_t j = 0; j < cuthistoskzerosh::cospacuts.size(); j++) {
              std::string cospacut = cuthistoskzerosh::cospacuts[j];
              size_t pos = cospacut.find("_");
              cospacut[pos] = '.';
              const float cospacutvalue = std::stod(cospacut);
              if (v0.v0cosPA() > cospacutvalue) {
                cuthistoskzerosh::cospaCut[j]->Fill(v0.mK0Short());
              }
            }
            for (uint32_t j = 0; j < cuthistoskzerosh::dcacuts.size(); j++) {
              std::string dcacut = cuthistoskzerosh::dcacuts[j];
              size_t pos = dcacut.find("_");
              dcacut[pos] = '.';
              const float dcacutvalue = std::stod(dcacut);
              if (v0.dcaV0daughters() < dcacutvalue) {
                cuthistoskzerosh::dcaCut[j]->Fill(v0.mK0Short());
              }
            }
            for (uint32_t j = 0; j < cuthistoskzerosh::v0radiuscuts.size(); j++) {
              std::string v0radiuscut = cuthistoskzerosh::v0radiuscuts[j];
              size_t pos = v0radiuscut.find("_");
              v0radiuscut[pos] = '.';
              const float v0radiuscutvalue = std::stod(v0radiuscut);
              if (v0.v0radius() > v0radiuscutvalue) {
                cuthistoskzerosh::v0radiusCut[j]->Fill(v0.mK0Short());
              }
            }
            for (uint32_t j = 0; j < cuthistoskzerosh::dcapostopvcuts.size(); j++) {
              std::string dcapostopcut = cuthistoskzerosh::dcapostopvcuts[j];
              size_t pos = dcapostopcut.find("_");
              dcapostopcut[pos] = '.';
              const float dcapostopcutvalue = std::stod(dcapostopcut);
              if (std::abs(v0.dcapostopv()) > dcapostopcutvalue) {
                cuthistoskzerosh::dcapostopCut[j]->Fill(v0.mK0Short());
              }
            }
            for (uint32_t j = 0; j < cuthistoskzerosh::dcanegtopvcuts.size(); j++) {
              std::string dcanegtopcut = cuthistoskzerosh::dcanegtopvcuts[j];
              size_t pos = dcanegtopcut.find("_");
              dcanegtopcut[pos] = '.';
              const float dcanegtopcutvalue = std::stod(dcanegtopcut);
              if (std::abs(v0.dcanegtopv()) > dcanegtopcutvalue) {
                cuthistoskzerosh::dcanegtopCut[j]->Fill(v0.mK0Short());
              }
            }
          }
        }
        if (std::abs(posDaughterTrack.tpcNSigmaPr()) < nSigmaTPCProton && std::abs(negDaughterTrack.tpcNSigmaPi()) < nSigmaTPCPion) { // TPC PID on daughter pion and proton
          rV0Parameters_Data.fill(HIST("hMassLambdaAllAfterCutsandTPCPID"), v0.mLambda());
          if (std::abs(v0.mK0Short() - mK0shPDG) > compv0masscut) { // lambda competitive v0 mass cut (cut out Kaons)
            // Filling the five Lambda invariant mass plots for different cuts (which are taken from namespace), same as with Kzeros above,for full explanation see the first kzero cut filling in the MC process
            for (uint32_t j = 0; j < cuthistoslambda::cospacuts.size(); j++) {
              std::string cospacutlambda = cuthistoslambda::cospacuts[j];
              size_t pos = cospacutlambda.find("_");
              cospacutlambda[pos] = '.';
              const float cospacutlambdavalue = std::stod(cospacutlambda);
              if (v0.v0cosPA() > cospacutlambdavalue) {
                cuthistoslambda::cospaCut[j]->Fill(v0.mLambda());
              }
            }
            for (uint32_t j = 0; j < cuthistoslambda::dcacuts.size(); j++) {
              std::string dcacutlambda = cuthistoslambda::dcacuts[j];
              size_t pos = dcacutlambda.find("_");
              dcacutlambda[pos] = '.';
              const float dcacutlambdavalue = std::stod(dcacutlambda);
              if (v0.dcaV0daughters() < dcacutlambdavalue) {
                cuthistoslambda::dcaCut[j]->Fill(v0.mLambda());
              }
            }
            for (uint32_t j = 0; j < cuthistoslambda::v0radiuscuts.size(); j++) {
              std::string v0radiuscutlambda = cuthistoslambda::v0radiuscuts[j];
              size_t pos = v0radiuscutlambda.find("_");
              v0radiuscutlambda[pos] = '.';
              const float v0radiuscutlambdavalue = std::stod(v0radiuscutlambda);
              if (v0.v0radius() > v0radiuscutlambdavalue) {
                cuthistoslambda::v0radiusCut[j]->Fill(v0.mLambda());
              }
            }
            for (uint32_t j = 0; j < cuthistoslambda::dcanegtopvcuts.size(); j++) {
              std::string dcapostopcutlambda = cuthistoslambda::dcapostopvcuts[j];
              size_t pos = dcapostopcutlambda.find("_");
              dcapostopcutlambda[pos] = '.';
              const float dcapostopcutlambdavalue = std::stod(dcapostopcutlambda);
              if (std::abs(v0.dcapostopv()) > dcapostopcutlambdavalue) {
                cuthistoslambda::dcapostopCut[j]->Fill(v0.mLambda());
              }
            }
            for (uint32_t j = 0; j < cuthistoslambda::dcanegtopvcuts.size(); j++) {
              std::string dcanegtopcutlambda = cuthistoslambda::dcanegtopvcuts[j];
              size_t pos = dcanegtopcutlambda.find("_");
              dcanegtopcutlambda[pos] = '.';
              const float dcanegtopcutlambdavalue = std::stod(dcanegtopcutlambda);
              if (std::abs(v0.dcanegtopv()) > dcanegtopcutlambdavalue) {
                cuthistoslambda::dcanegtopCut[j]->Fill(v0.mLambda());
              }
            }
          }
        }
        if (std::abs(posDaughterTrack.tpcNSigmaPr()) < nSigmaTPCProton && std::abs(negDaughterTrack.tpcNSigmaPi()) < nSigmaTPCPion) { // TPC PID on daughter pion abd proton
          rV0Parameters_Data.fill(HIST("hMassAntilambdaAllAfterCutsandTPCPID"), v0.mAntiLambda());
          if (std::abs(v0.mK0Short() - mK0shPDG) > compv0masscut) { // antilambda competitive v0 mass cut (cut out Kaons)
            // Filling the five Anti-Lambda invariant mass plots for different cuts (which are taken from namespace), same as with Kzeros and Lambdas above,for full explanation see the first kzero cut filling in the MC process
            for (uint32_t j = 0; j < cuthistosantilambda::cospacuts.size(); j++) {
              std::string cospacutantilambda = cuthistosantilambda::cospacuts[j];
              size_t pos = cospacutantilambda.find("_");
              cospacutantilambda[pos] = '.';
              const float cospacutantilambdavalue = std::stod(cospacutantilambda);
              if (v0.v0cosPA() > cospacutantilambdavalue) {
                cuthistosantilambda::cospaCut[j]->Fill(v0.mAntiLambda());
              }
            }
            for (uint32_t j = 0; j < cuthistosantilambda::dcacuts.size(); j++) {
              std::string dcacutantilambda = cuthistosantilambda::dcacuts[j];
              size_t pos = dcacutantilambda.find("_");
              dcacutantilambda[pos] = '.';
              const float dcacutantilambdavalue = std::stod(dcacutantilambda);
              if (v0.dcaV0daughters() < dcacutantilambdavalue) {
                cuthistosantilambda::dcaCut[j]->Fill(v0.mAntiLambda());
              }
            }
            for (uint32_t j = 0; j < cuthistosantilambda::v0radiuscuts.size(); j++) {
              std::string v0radiusantilambda = cuthistosantilambda::v0radiuscuts[j];
              size_t pos = v0radiusantilambda.find("_");
              v0radiusantilambda[pos] = '.';
              const float v0radiuscutantilambdavalue = std::stod(v0radiusantilambda);
              if (v0.v0radius() > v0radiuscutantilambdavalue) {
                cuthistosantilambda::v0radiusCut[j]->Fill(v0.mAntiLambda());
              }
            }
            for (uint32_t j = 0; j < cuthistosantilambda::dcapostopvcuts.size(); j++) {
              std::string dcapostopantilambda = cuthistosantilambda::dcapostopvcuts[j];
              size_t pos = dcapostopantilambda.find("_");
              dcapostopantilambda[pos] = '.';
              const float dcapostopcutantilambdavalue = std::stod(dcapostopantilambda);
              if (std::abs(v0.dcapostopv()) > dcapostopcutantilambdavalue) {
                cuthistosantilambda::dcapostopCut[j]->Fill(v0.mAntiLambda());
              }
            }
            for (uint32_t j = 0; j < cuthistosantilambda::dcanegtopvcuts.size(); j++) {
              std::string dcanegtopantilambda = cuthistosantilambda::dcanegtopvcuts[j];
              size_t pos = dcanegtopantilambda.find("_");
              dcanegtopantilambda[pos] = '.';
              const float dcanegtopcutantilambdavalue = std::stod(dcanegtopantilambda);
              if (std::abs(v0.dcanegtopv()) > dcanegtopcutantilambdavalue) {
                cuthistosantilambda::dcanegtopCut[j]->Fill(v0.mAntiLambda());
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(v0topologicalcuts, RecMCprocess, "Process Run 3 MC:Reconstructed", true);
  PROCESS_SWITCH(v0topologicalcuts, Dataprocess, "Process Run 3 Data,", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0topologicalcuts>(cfgc)};
}
