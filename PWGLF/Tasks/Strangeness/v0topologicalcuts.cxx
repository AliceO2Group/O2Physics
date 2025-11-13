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
/// \file v0topologicalcuts.cxx
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

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include "CommonUtils/StringUtils.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <memory>
#include <string>
#include <vector>

// namespaces to be used for the plot names and topological cuts that will be given by a configurable string
namespace cuthistoskzerosh
{
std::shared_ptr<TH1> cosPACut[20];
static std::vector<std::string> cosPAcuts;
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
std::shared_ptr<TH1> cosPACut[20];
static std::vector<std::string> cosPAcuts;
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
std::shared_ptr<TH1> cosPACut[20];
static std::vector<std::string> cosPAcuts;
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
  HistogramRegistry rV0ParametersMCV0match{"V0ParametersMCV0Match", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rV0ParametersMCK0Smatch{"V0ParametersMCK0SMatch", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rV0ParametersMCLambdamatch{"V0ParametersMCLambdaMatch", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rV0ParametersMCAntiLambdamatch{"V0ParametersMCAntiLambdaMatch", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rV0ParametersData{"rV0ParametersData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // kzero cut Histogram Registry with MC-matching, each will include 20 histograms for 20 different cuts
  HistogramRegistry rKzeroShortCosPACut{"KzeroShortCosPACuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShortDCACut{"KzeroShortDCACuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShortV0radiusCut{"KzeroShortV0radiusCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShortDCApostopCut{"KzeroShortDCApostopvCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShortDCAnegtopCut{"KzeroShortDCAnegtopvCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // lambdas cut histograms with MC-matching (same as in Kzeros above)
  HistogramRegistry rLambdaCosPACut{"LambdaCosPACuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambdaDCACut{"LambdaDCACuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambdaV0radiusCut{"LambdaV0radiusCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambdaDCApostopCut{"LambdaDCApostopvCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambdaDCAnegtopCut{"LambdaDCAnegtopvCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // antilambdas cut histograms with MC-matching (same as in Lambdas an Kzeros above)
  HistogramRegistry rAntiLambdaCosPACut{"AntiLambdaCosPACuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntiLambdaDCACut{"AntiLambdaDCACuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntiLambdaV0radiusCut{"AntiLambdaV0radiusCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntiLambdaDCApostopCut{"AntiLambdaDCApostopvCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntiLambdaDCAnegtopCut{"AntiLambdaDCAnegtopvCuts", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurables for Cuts
  Configurable<float> cutZVertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> nSigmaTPCPion{"nSigmaTPCPion", 4, "nSigmaTPCPion"};
  Configurable<float> nSigmaTPCProton{"nSigmaTPCProton", 4, "nSigmaTPCProton"};
  Configurable<float> compv0masscut{"compv0masscut", 0.01, "CompetitiveV0masscut (GeV)"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};

  // Configurable strings for Kzero cuts
  Configurable<std::string> kzeroshSettingCosPAcutsString{"kzeroshSettingCosPAcutsString", {"0_98,0_981,0_982,0_983,0_984,0_985,0_986,0_987,0_988,0_989,0_99,0_991,0_992,0_993,0_994,0_995,0_996,0_997,0_998,0_999"}, "Kzero cosPA Cut Values"};
  Configurable<std::string> kzeroshSettingDCAcutsString{"kzeroshSettingDCAcutsString", {"0_3,0_285,0_27,0_255,0_24,0_225,0_21,0_195,0_18,0_165,0_15,0_135,0_12,0_105,0_09,0_075,0_06,0_045,0_03,0_015"}, "Kzero DCA Cut Values"};
  Configurable<std::string> kzeroshSettingV0radiusString{"kzeroshSettingV0radiusString", {"0_5,0_51,0_52,0_53,0_54,0_55,0_56,0_57,0_58,0_59,0_6,0_61,0_62,0_63,0_64,0_65,0_66,0_67,0_68,0_69"}, "Kzero V0Radius Cut Values"};
  Configurable<std::string> kzeroshSettingDCApostopvString{"kzeroshSettingDCApostopvString", {"0_0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Kzero DCA Pos to PV Cut Values"};
  Configurable<std::string> kzeroshSettingDCAnegtopvString{"kzeroshSettingDCAnegtopvString", {"0_0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "KzeroDCA Neg to PV Cut Values"};

  // Configurable strings for Lambdacuts
  Configurable<std::string> lambdaSettingCosPAcutsString{"lambdaSettingCosPAcutsString", {"0_98,0_981,0_982,0_983,0_984,0_985,0_986,0_987,0_988,0_989,0_99,0_991,0_992,0_993,0_994"}, "Lambda cosPA Cut Values"};
  Configurable<std::string> lambdaSettingDCAcutsString{"lambdaSettingDCAcutsString", {"0_3,0_285,0_27,0_255,0_24,0_225,0_21,0_195,0_18,0_165,0_15,0_135,0_12,0_105,0_09,0_075,0_06,0_045,0_03,0_015"}, "Lambda DCA Cut Values"};
  Configurable<std::string> lambdaSettingV0radiusString{"lambdaSettingV0radiusString", {"0_5,0_51,0_52,0_53,0_54,0_55,0_56,0_57,0_58,0_59,0_6,0_61,0_62,0_63,0_64,0_65,0_66,0_67,0_68,0_69"}, "Lambda V0Radius Cut Values"};
  Configurable<std::string> lambdaSettingDCApostopvString{"lambdaSettingDCApostopvString", {"0_0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Lambda DCA Pos to PV Cut Values"};
  Configurable<std::string> lambdaSettingDCAnegtopvString{"lambdaSettingDCAnegtopvString", {"0_0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Lambda DCA Neg to PV Cut Values"};

  // Configurable strings for AntiLambdacuts
  Configurable<std::string> antilambdaSettingCosPAcutsString{"antilambdaSettingCosPAcutsString", {"0_98,0_981,0_982,0_983,0_984,0_985,0_986,0_987,0_988,0_989,0_99,0_991,0_992,0_993,0_994,0_995,0_996,0_997,0_998,0_999"}, "Antilambda cosPA Cut Values"};
  Configurable<std::string> antilambdaSettingDCAcutsString{"antilambdaSettingDCAcutsString", {"0_3,0_285,0_27,0_255,0_24,0_225,0_21,0_195,0_18,0_165,0_15,0_135,0_12,0_105,0_09,0_075,0_06,0_045,0_03,0_015"}, "Antilambda DCA Cut Values"};
  Configurable<std::string> antilambdaSettingV0radiusString{"antilambdaSettingV0radiusString", {"0_5,0_51,0_52,0_53,0_54,0_55,0_56,0_57,0_58,0_59,0_6,0_61,0_62,0_63,0_64,0_65,0_66,0_67,0_68,0_69"}, "Antilambda V0Radius Cut Values"};
  Configurable<std::string> antilambdaSettingDCApostopvString{"antilambdaSettingDCApostopvString", {"0_0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Antilambda DCA Pos to PV Cut Values"};
  Configurable<std::string> antilambdaSettingDCAnegtopvString{"antilambdaSettingDCAnegtopvString", {"0_0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Antilambda DCA Neg to PV Cut Values"};

  void init(InitContext const&)
  {
    // kzero filling namespace with configurable strings

    // Setting strings from configurable strings in order to manipulate them
    // getting the  cut values for the names of the plots for the five topological cuts
    cuthistoskzerosh::cosPAcuts = o2::utils::Str::tokenize(kzeroshSettingCosPAcutsString, ',');
    cuthistoskzerosh::dcacuts = o2::utils::Str::tokenize(kzeroshSettingDCAcutsString, ',');
    cuthistoskzerosh::v0radiuscuts = o2::utils::Str::tokenize(kzeroshSettingV0radiusString, ',');
    cuthistoskzerosh::dcapostopvcuts = o2::utils::Str::tokenize(kzeroshSettingDCApostopvString, ',');
    cuthistoskzerosh::dcanegtopvcuts = o2::utils::Str::tokenize(kzeroshSettingDCAnegtopvString, ',');

    // lambda filling namespace with configurable strings (same as in Kzeros above)
    cuthistoslambda::cosPAcuts = o2::utils::Str::tokenize(lambdaSettingCosPAcutsString, ',');
    cuthistoslambda::dcacuts = o2::utils::Str::tokenize(lambdaSettingDCAcutsString, ',');
    cuthistoslambda::v0radiuscuts = o2::utils::Str::tokenize(lambdaSettingV0radiusString, ',');
    cuthistoslambda::dcapostopvcuts = o2::utils::Str::tokenize(lambdaSettingDCApostopvString, ',');
    cuthistoslambda::dcanegtopvcuts = o2::utils::Str::tokenize(lambdaSettingDCAnegtopvString, ',');

    // antilambda filling namespace with configurable strings (same as in Lambdas and Kzeros above)
    cuthistosantilambda::cosPAcuts = o2::utils::Str::tokenize(antilambdaSettingCosPAcutsString, ',');
    cuthistosantilambda::dcacuts = o2::utils::Str::tokenize(antilambdaSettingDCAcutsString, ',');
    cuthistosantilambda::v0radiuscuts = o2::utils::Str::tokenize(antilambdaSettingV0radiusString, ',');
    cuthistosantilambda::dcapostopvcuts = o2::utils::Str::tokenize(antilambdaSettingDCApostopvString, ',');
    cuthistosantilambda::dcanegtopvcuts = o2::utils::Str::tokenize(antilambdaSettingDCAnegtopvString, ',');

    // Axes for the three invariant mass plots
    AxisSpec k0ShortMassAxis = {nBins, 0.45f, 0.55f, "#it{M} #pi^{+}#pi^{-} [GeV/#it{c}^{2}]"};
    AxisSpec lambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M} p^{+}#pi^{-} [GeV/#it{c}^{2}]"};
    AxisSpec antiLambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M} p^{-}#pi^{+} [GeV/#it{c}^{2}]"};
    AxisSpec ptAxis = {nBins, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};

    // adding the invariant mass histograms to their Registries using the namespace for kzeros, lambdas and antilambdas
    for (uint32_t i = 0; i < cuthistoskzerosh::cosPAcuts.size(); i++) {
      cuthistoskzerosh::cosPACut[i] = rKzeroShortCosPACut.add<TH1>(fmt::format("hKzerocosPACut_{}", cuthistoskzerosh::cosPAcuts[i]).data(), fmt::format("hKzerocosPACut_{}", cuthistoskzerosh::cosPAcuts[i]).data(), {HistType::kTH1D, {{k0ShortMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoskzerosh::dcacuts.size(); i++) {
      cuthistoskzerosh::dcaCut[i] = rKzeroShortDCACut.add<TH1>(fmt::format("hKzerodcaCut_{}", cuthistoskzerosh::dcacuts[i]).data(), fmt::format("hKzerodcaCut_{}", cuthistoskzerosh::dcacuts[i]).data(), {HistType::kTH1D, {{k0ShortMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoskzerosh::v0radiuscuts.size(); i++) {
      cuthistoskzerosh::v0radiusCut[i] = rKzeroShortV0radiusCut.add<TH1>(fmt::format("hKzerov0radiusCut_{}", cuthistoskzerosh::v0radiuscuts[i]).data(), fmt::format("hKzerov0radiusCut_{}", cuthistoskzerosh::v0radiuscuts[i]).data(), {HistType::kTH1D, {{k0ShortMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoskzerosh::dcapostopvcuts.size(); i++) {
      cuthistoskzerosh::dcapostopCut[i] = rKzeroShortDCApostopCut.add<TH1>(fmt::format("hKzerodcapostopCut_{}", cuthistoskzerosh::dcapostopvcuts[i]).data(), fmt::format("hKzerodcapostopCut_{}", cuthistoskzerosh::dcapostopvcuts[i]).data(), {HistType::kTH1D, {{k0ShortMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoskzerosh::dcanegtopvcuts.size(); i++) {
      cuthistoskzerosh::dcanegtopCut[i] = rKzeroShortDCAnegtopCut.add<TH1>(fmt::format("hKzerodcanegtopCut_{}", cuthistoskzerosh::dcanegtopvcuts[i]).data(), fmt::format("hKzerodcanegtopCut_{}", cuthistoskzerosh::dcanegtopvcuts[i]).data(), {HistType::kTH1D, {{k0ShortMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoslambda::cosPAcuts.size(); i++) {
      cuthistoslambda::cosPACut[i] = rLambdaCosPACut.add<TH1>(fmt::format("hLambdacosPACut_{}", cuthistoslambda::cosPAcuts[i]).data(), fmt::format("hLambdacosPACut_{}", cuthistoslambda::cosPAcuts[i]).data(), {HistType::kTH1D, {{lambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoslambda::dcacuts.size(); i++) {
      cuthistoslambda::dcaCut[i] = rLambdaDCACut.add<TH1>(fmt::format("hLambdadcaCut_{}", cuthistoslambda::dcacuts[i]).data(), fmt::format("hLambdadcaCut_{}", cuthistoslambda::dcacuts[i]).data(), {HistType::kTH1D, {{lambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoslambda::v0radiuscuts.size(); i++) {
      cuthistoslambda::v0radiusCut[i] = rLambdaV0radiusCut.add<TH1>(fmt::format("hLambdav0radiusCut_{}", cuthistoslambda::v0radiuscuts[i]).data(), fmt::format("hLambdav0radiusCut_{}", cuthistoslambda::v0radiuscuts[i]).data(), {HistType::kTH1D, {{lambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoslambda::dcapostopvcuts.size(); i++) {
      cuthistoslambda::dcapostopCut[i] = rLambdaDCApostopCut.add<TH1>(fmt::format("hLambdadcapostopCut_{}", cuthistoslambda::dcapostopvcuts[i]).data(), fmt::format("hLambdadcapostopCut_{}", cuthistoslambda::dcapostopvcuts[i]).data(), {HistType::kTH1D, {{lambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistoslambda::dcanegtopvcuts.size(); i++) {
      cuthistoslambda::dcanegtopCut[i] = rLambdaDCAnegtopCut.add<TH1>(fmt::format("hLambdadcanegtopCut_{}", cuthistoslambda::dcanegtopvcuts[i]).data(), fmt::format("hLambdadcanegtopCut_{}", cuthistoslambda::dcanegtopvcuts[i]).data(), {HistType::kTH1D, {{lambdaMassAxis}}});
    }

    for (uint32_t i = 0; i < cuthistosantilambda::cosPAcuts.size(); i++) {
      cuthistosantilambda::cosPACut[i] = rAntiLambdaCosPACut.add<TH1>(fmt::format("hAntiLambdacosPACut_{}", cuthistosantilambda::cosPAcuts[i]).data(), fmt::format("hAntiLambdacosPACut_{}", cuthistosantilambda::cosPAcuts[i]).data(), {HistType::kTH1D, {{antiLambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistosantilambda::dcacuts.size(); i++) {
      cuthistosantilambda::dcaCut[i] = rAntiLambdaDCACut.add<TH1>(fmt::format("hAntiLambdadcaCut_{}", cuthistosantilambda::dcacuts[i]).data(), fmt::format("hAntiLambdadcaCut_{}", cuthistosantilambda::dcacuts[i]).data(), {HistType::kTH1D, {{antiLambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistosantilambda::v0radiuscuts.size(); i++) {
      cuthistosantilambda::v0radiusCut[i] = rAntiLambdaV0radiusCut.add<TH1>(fmt::format("hAntiLambdav0radiusCut_{}", cuthistosantilambda::v0radiuscuts[i]).data(), fmt::format("hAntiLambdav0radiusCut_{}", cuthistosantilambda::v0radiuscuts[i]).data(), {HistType::kTH1D, {{antiLambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistosantilambda::dcapostopvcuts.size(); i++) {
      cuthistosantilambda::dcapostopCut[i] = rAntiLambdaDCApostopCut.add<TH1>(fmt::format("hAntiLambdadcapostopCut_{}", cuthistosantilambda::dcapostopvcuts[i]).data(), fmt::format("hAntiLambdadcapostopCut_{}", cuthistosantilambda::dcapostopvcuts[i]).data(), {HistType::kTH1D, {{antiLambdaMassAxis}}});
    }
    for (uint32_t i = 0; i < cuthistosantilambda::dcanegtopvcuts.size(); i++) {
      cuthistosantilambda::dcanegtopCut[i] = rAntiLambdaDCAnegtopCut.add<TH1>(fmt::format("hAntiLambdadcanegtopCut_{}", cuthistosantilambda::dcanegtopvcuts[i]).data(), fmt::format("hAntiLambdadcanegtopCut_{}", cuthistosantilambda::dcanegtopvcuts[i]).data(), {HistType::kTH1D, {{antiLambdaMassAxis}}});
    }

    // K0s topological cut histograms added and MC-matched
    rV0ParametersMCV0match.add("hDCAV0Daughters_V0_Match", "hDCAV0Daughters_No_Match", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0ParametersMCV0match.add("hV0CosPA_V0_Match", "hV0CosPA_No_Match", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0ParametersMCV0match.add("hV0Radius_V0_Match", "hV0Radius_No_Match", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});
    rV0ParametersMCV0match.add("hV0Radius_V0_Match_Full", "hV0Radius_No_Match_Full", {HistType::kTH1F, {{nBins, 0.0f, 40.0f}}});
    rV0ParametersMCV0match.add("hDCAPostoPV_V0_Match", "hDCAPostoPV_No_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0ParametersMCV0match.add("hDCANegtoPV_V0_Match", "hDCANegtoPV_No_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0ParametersMCV0match.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // kzero match
    rV0ParametersMCK0Smatch.add("hDCAV0Daughters_KzeroMC_Match", "hDCAV0Daughters_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0ParametersMCK0Smatch.add("hV0CosPA_KzeroMC_Match", "hV0CosPA_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0ParametersMCK0Smatch.add("hV0Radius_KzeroMC_Match", "hV0Radius_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.2f, 5.0f}}});
    rV0ParametersMCK0Smatch.add("hDCAPostoPV_KzeroMC_Match", "hDCAPostoPV_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0ParametersMCK0Smatch.add("hDCANegtoPV_KzeroMC_Match", "hDCANegtoPV_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});

    // lambda match
    rV0ParametersMCLambdamatch.add("hDCAV0Daughters_LambdaMC_Match", "hDCAV0Daughters_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0ParametersMCLambdamatch.add("hV0CosPA_LambdaMC_Match", "hV0CosPA_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0ParametersMCLambdamatch.add("hV0Radius_LambdaMC_Match", "hV0Radius_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.2f, 5.0f}}});
    rV0ParametersMCLambdamatch.add("hDCAPostoPV_LambdaMC_Match", "hDCAPostoPV_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0ParametersMCLambdamatch.add("hDCANegtoPV_LambdaMC_Match", "hDCANegtoPV_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});

    // antilambda match
    rV0ParametersMCAntiLambdamatch.add("hDCAV0Daughters_AntiLambdaMC_Match", "hDCAV0Daughters_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0ParametersMCAntiLambdamatch.add("hV0CosPA_AntiLambdaMC_Match", "hV0CosPA_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0ParametersMCAntiLambdamatch.add("hV0Radius_AntiLambdaMC_Match", "hV0Radius_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.2f, 5.0f}}});
    rV0ParametersMCAntiLambdamatch.add("hDCAPostoPV_AntiLambdaMC_Match", "hDCAPostoPV_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0ParametersMCAntiLambdamatch.add("hDCANegtoPV_AntiLambdaMC_Match", "hDCANegtoPV_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});

    // V0s Data
    rV0ParametersData.add("hDCAV0Daughters_V0Data", "hDCAV0Daughters_V0Data", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0ParametersData.add("hV0CosPA_V0Data", "hV0CosPA_V0Data", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0ParametersData.add("hV0Radius_V0Data", "hV0Radius_V0Data", {HistType::kTH1F, {{nBins, 0.2f, 5.0f}}});
    rV0ParametersData.add("hV0Radius_Full_V0Data", "hV0Radius_Full_V0Data", {HistType::kTH1F, {{nBins, 0.2f, 40.0f}}});
    rV0ParametersData.add("hDCAPostoPV_V0Data", "hDCAPostoPV_V0Data", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0ParametersData.add("hDCANegtoPV_V0Data", "hDCANegtoPV_V0Data", {HistType::kTH1F, {{nBins, 0.0f, 2.0f}}});
    rV0ParametersData.add("hMassK0ShortNoCuts_V0Data", "hMassK0ShortNoCuts_V0Data", {HistType::kTH1F, {{k0ShortMassAxis}}});
    rV0ParametersData.add("hMassLambdaNoCuts_V0Data", "hMassLambdaNoCuts_V0Data", {HistType::kTH1F, {{lambdaMassAxis}}});
    rV0ParametersData.add("hMassAntilambdaNoCuts_V0Data", "hMassAntilambdaNoCuts_V0Data", {HistType::kTH1F, {{antiLambdaMassAxis}}});
    rV0ParametersData.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    // K0short Data
    rV0ParametersData.add("hMassK0ShortAfterEtaCut", "hMassK0ShortAfterEtaCut", {HistType::kTH1F, {k0ShortMassAxis}});
    rV0ParametersData.add("hMassK0ShortAfterCompmassCut", "hMassK0ShortAfterCompmassCut", {HistType::kTH1F, {k0ShortMassAxis}});
    rV0ParametersData.add("hMassK0ShortAfterAllCuts", "hMassK0ShortAfterAllCuts", {HistType::kTH1F, {k0ShortMassAxis}});
    rV0ParametersData.add("hNSigmaPosPiFromK0s", "hNSigmaPosPiFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rV0ParametersData.add("hNSigmaNegPiFromK0s", "hNSigmaNegPiFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rV0ParametersData.add("hK0shEtaPosDau", "hK0shEtaPosDau", {HistType::kTH1F, {{nBins, -1.2f, 1.2f}}});
    rV0ParametersData.add("hK0shEtaNegDau", "hK0shEtaNegDau", {HistType::kTH1F, {{nBins, -1.2f, 1.2f}}});
    rV0ParametersData.add("hK0shEtaPosDauAfterCut", "hK0shEtaPosDauAfterCut", {HistType::kTH1F, {{nBins, -1.2f, 1.2f}}});
    rV0ParametersData.add("hK0shEtaNegDauAfterCut", "hK0shEtaNegDauAfterCut", {HistType::kTH1F, {{nBins, -1.2f, 1.2f}}});

    // Lambda Data
    rV0ParametersData.add("hMassLambdaAfterEtaCut", "hMassLambdaAfterEtaCut", {HistType::kTH1F, {lambdaMassAxis}});
    rV0ParametersData.add("hMassLambdaAfterCompmassCut", "hMassLambdaAfterCompmassCut", {HistType::kTH1F, {lambdaMassAxis}});
    rV0ParametersData.add("hMassLambdaAfterAllCuts", "hMassLambdaAfterAllCuts", {HistType::kTH1F, {lambdaMassAxis}});
    rV0ParametersData.add("hNSigmaPosProtonFromLambda", "hNSigmaPosProtonFromLambda", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rV0ParametersData.add("hNSigmaNegPionFromLambda", "hNSigmaNegPionFromLambda", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    // Antilambda Data
    rV0ParametersData.add("hMassAntiLambdaAfterEtaCut", "hMassAntiLambdaAfterEtaCut", {HistType::kTH1F, {antiLambdaMassAxis}});
    rV0ParametersData.add("hMassAntiLambdaAfterCompmassCut", "hMassAntiLambdaAfterCompmassCut", {HistType::kTH1F, {antiLambdaMassAxis}});
    rV0ParametersData.add("hMassAntiLambdaAfterAllCuts", "hMassAntiLambdaAfterAllCuts", {HistType::kTH1F, {antiLambdaMassAxis}});
    rV0ParametersData.add("hNSigmaNegProtonFromAntilambda", "hNSigmaNegProtonFromAntilambda", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rV0ParametersData.add("hNSigmaPosPionFromAntilambda", "hNSigmaPosPionFromAntilambda", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < cutZVertex);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutZVertex);

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;

  // This is the Process for the MC reconstructed Data
  void recMCProcess(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                    soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s,
                    DaughterTracks const&, // no need to define a variable for tracks, if we don't access them directly
                    aod::McParticles const&)
  {
    const auto& mLambdaPDG = 1.115683;
    const auto& mK0shPDG = 0.497611;
    for (const auto& v0 : V0s) {
      if (std::abs(v0.posTrack_as<DaughterTracks>().eta()) < etadau && std::abs(v0.negTrack_as<DaughterTracks>().eta()) < etadau) { // daughters pseudorapidity cut
        // filling histograms with V0 values
        rV0ParametersMCV0match.fill(HIST("hDCAV0Daughters_V0_Match"), v0.dcaV0daughters());
        rV0ParametersMCV0match.fill(HIST("hV0CosPA_V0_Match"), v0.v0cosPA());
        rV0ParametersMCV0match.fill(HIST("hV0Radius_V0_Match"), v0.v0radius());
        rV0ParametersMCV0match.fill(HIST("hV0Radius_V0_Match_Full"), v0.v0radius());
        rV0ParametersMCV0match.fill(HIST("hDCAPostoPV_V0_Match"), v0.dcapostopv());
        rV0ParametersMCV0match.fill(HIST("hDCANegtoPV_V0_Match"), v0.dcanegtopv());
        rV0ParametersMCV0match.fill(HIST("hVertexZRec"), collision.posZ());

        // Checking that the V0 is a true K0s/Lambdas/Antilambdas and then filling the parameter histograms and the invariant mass plots for different cuts (which are taken from namespace)
        if (v0.has_mcParticle()) {
          auto v0mcParticle = v0.mcParticle();
          if (v0mcParticle.pdgCode() == 310) {                                                                                    // kzero matched
            if (std::abs(v0.mLambda() - mLambdaPDG) > compv0masscut && std::abs(v0.mAntiLambda() - mLambdaPDG) > compv0masscut) { // Kzero competitive v0 mass cut (cut out Lambdas and Anti-Lambdas)
              rV0ParametersMCK0Smatch.fill(HIST("hDCAV0Daughters_KzeroMC_Match"), v0.dcaV0daughters());
              rV0ParametersMCK0Smatch.fill(HIST("hV0CosPA_KzeroMC_Match"), v0.v0cosPA());
              rV0ParametersMCK0Smatch.fill(HIST("hV0Radius_KzeroMC_Match"), v0.v0radius());
              rV0ParametersMCK0Smatch.fill(HIST("hDCAPostoPV_KzeroMC_Match"), std::abs(v0.dcapostopv()));
              rV0ParametersMCK0Smatch.fill(HIST("hDCANegtoPV_KzeroMC_Match"), std::abs(v0.dcanegtopv()));

              for (uint32_t j = 0; j < cuthistoskzerosh::cosPAcuts.size(); j++) {
                std::string cosPAcut = cuthistoskzerosh::cosPAcuts[j]; // Get the current cut value from the namespace
                size_t pos = cosPAcut.find("_");                       // find the "_" which needs to change to a "." for it to be a number
                cosPAcut[pos] = '.';                                   // change the "_" into an "."
                const float cosPAcutvalue = std::stod(cosPAcut);       // make the string into a float value
                if (v0.v0cosPA() > cosPAcutvalue) {                    // enforce the cut value
                  cuthistoskzerosh::cosPACut[j]->Fill(v0.mK0Short());  // fill the corresponding histo from the namespace with the invariant mass (of a Kzero here)
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
              rV0ParametersMCLambdamatch.fill(HIST("hDCAV0Daughters_LambdaMC_Match"), v0.dcaV0daughters());
              rV0ParametersMCLambdamatch.fill(HIST("hV0CosPA_LambdaMC_Match"), v0.v0cosPA());
              rV0ParametersMCLambdamatch.fill(HIST("hV0Radius_LambdaMC_Match"), v0.v0radius());
              rV0ParametersMCLambdamatch.fill(HIST("hDCAPostoPV_LambdaMC_Match"), std::abs(v0.dcapostopv()));
              rV0ParametersMCLambdamatch.fill(HIST("hDCANegtoPV_LambdaMC_Match"), std::abs(v0.dcanegtopv()));

              // for explanation look at the first Kzero  plot above
              for (uint32_t j = 0; j < cuthistoslambda::cosPAcuts.size(); j++) {
                std::string cosPAcutlambda = cuthistoslambda::cosPAcuts[j];
                size_t pos = cosPAcutlambda.find("_");
                cosPAcutlambda[pos] = '.';
                const float cosPAcutlambdavalue = std::stod(cosPAcutlambda);
                if (v0.v0cosPA() > cosPAcutlambdavalue) {
                  cuthistoslambda::cosPACut[j]->Fill(v0.mLambda());
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
              rV0ParametersMCAntiLambdamatch.fill(HIST("hDCAV0Daughters_AntiLambdaMC_Match"), v0.dcaV0daughters());
              rV0ParametersMCAntiLambdamatch.fill(HIST("hV0CosPA_AntiLambdaMC_Match"), v0.v0cosPA());
              rV0ParametersMCAntiLambdamatch.fill(HIST("hV0Radius_AntiLambdaMC_Match"), v0.v0radius());
              rV0ParametersMCAntiLambdamatch.fill(HIST("hDCAPostoPV_AntiLambdaMC_Match"), std::abs(v0.dcapostopv()));
              rV0ParametersMCAntiLambdamatch.fill(HIST("hDCANegtoPV_AntiLambdaMC_Match"), std::abs(v0.dcanegtopv()));
              // for explanation look at the first Kzero  plot above
              for (uint32_t j = 0; j < cuthistosantilambda::cosPAcuts.size(); j++) {
                std::string cosPAcutantilambda = cuthistosantilambda::cosPAcuts[j];
                size_t pos = cosPAcutantilambda.find("_");
                cosPAcutantilambda[pos] = '.';
                const float cosPAcutantilambdavalue = std::stod(cosPAcutantilambda);
                if (v0.v0cosPA() > cosPAcutantilambdavalue) {
                  cuthistosantilambda::cosPACut[j]->Fill(v0.mAntiLambda());
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
  void dataProcess(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                   aod::V0Datas const& V0s,
                   DaughterTracks const&)
  {
    // filling histograms with the different V0 parameters
    const auto& mLambdaPDG = 1.115683;
    const auto& mK0shPDG = 0.497611;
    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<DaughterTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<DaughterTracks>();
      rV0ParametersData.fill(HIST("hMassK0ShortNoCuts_V0Data"), v0.mK0Short());
      rV0ParametersData.fill(HIST("hMassLambdaNoCuts_V0Data"), v0.mLambda());
      rV0ParametersData.fill(HIST("hMassAntilambdaNoCuts_V0Data"), v0.mAntiLambda());
      rV0ParametersData.fill(HIST("hDCAV0Daughters_V0Data"), v0.dcaV0daughters());
      rV0ParametersData.fill(HIST("hV0CosPA_V0Data"), v0.v0cosPA());
      rV0ParametersData.fill(HIST("hV0Radius_V0Data"), v0.v0radius());
      rV0ParametersData.fill(HIST("hV0Radius_Full_V0Data"), v0.v0radius());
      rV0ParametersData.fill(HIST("hDCAPostoPV_V0Data"), std::abs(v0.dcapostopv()));
      rV0ParametersData.fill(HIST("hDCANegtoPV_V0Data"), std::abs(v0.dcanegtopv()));
      rV0ParametersData.fill(HIST("hVertexZRec"), collision.posZ());
      rV0ParametersData.fill(HIST("hK0shEtaPosDau"), v0.posTrack_as<DaughterTracks>().eta());
      rV0ParametersData.fill(HIST("hK0shEtaNegDau"), v0.negTrack_as<DaughterTracks>().eta());

      if (std::abs(v0.posTrack_as<DaughterTracks>().eta()) < etadau && std::abs(v0.negTrack_as<DaughterTracks>().eta()) < etadau) { // daughters pseudorapidity cut
        rV0ParametersData.fill(HIST("hMassK0ShortAfterEtaCut"), v0.mK0Short());
        rV0ParametersData.fill(HIST("hMassLambdaAfterEtaCut"), v0.mLambda());
        rV0ParametersData.fill(HIST("hMassAntiLambdaAfterEtaCut"), v0.mAntiLambda());
        rV0ParametersData.fill(HIST("hK0shEtaPosDauAfterCut"), v0.posTrack_as<DaughterTracks>().eta());
        rV0ParametersData.fill(HIST("hK0shEtaNegDauAfterCut"), v0.negTrack_as<DaughterTracks>().eta());
        if (std::abs(v0.mLambda() - mLambdaPDG) > compv0masscut && std::abs(v0.mAntiLambda() - mLambdaPDG) > compv0masscut) { // antilambda competitive v0 mass cut (cut out Lambdas and Anti-Lambdas)
          rV0ParametersData.fill(HIST("hMassK0ShortAfterCompmassCut"), v0.mK0Short());
          if (std::abs(posDaughterTrack.tpcNSigmaPi()) < nSigmaTPCPion && std::abs(negDaughterTrack.tpcNSigmaPi()) < nSigmaTPCPion) { // TPC PID on daughter pions
            rV0ParametersData.fill(HIST("hMassK0ShortAfterAllCuts"), v0.mK0Short());
            rV0ParametersData.fill(HIST("hNSigmaPosPiFromK0s"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tpcInnerParam());
            rV0ParametersData.fill(HIST("hNSigmaNegPiFromK0s"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
            // Filling the five Kzero invariant mass plots for different cuts (which are taken from namespace), for full explanation see the first kzero cut filling in the MC process
            for (uint32_t j = 0; j < cuthistoskzerosh::cosPAcuts.size(); j++) {
              std::string cosPAcut = cuthistoskzerosh::cosPAcuts[j];
              size_t pos = cosPAcut.find("_");
              cosPAcut[pos] = '.';
              const float cosPAcutvalue = std::stod(cosPAcut);
              if (v0.v0cosPA() > cosPAcutvalue) {
                cuthistoskzerosh::cosPACut[j]->Fill(v0.mK0Short());
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
        if (std::abs(v0.mK0Short() - mK0shPDG) > compv0masscut) { // lambda competitive v0 mass cut (cut out Kaons)
          rV0ParametersData.fill(HIST("hMassLambdaAfterCompmassCut"), v0.mLambda());
          rV0ParametersData.fill(HIST("hMassAntiLambdaAfterCompmassCut"), v0.mAntiLambda());
          if (std::abs(posDaughterTrack.tpcNSigmaPr()) < nSigmaTPCProton && std::abs(negDaughterTrack.tpcNSigmaPi()) < nSigmaTPCPion) { // TPC PID on daughter pion and proton for Lambda
            rV0ParametersData.fill(HIST("hMassLambdaAfterAllCuts"), v0.mLambda());
            rV0ParametersData.fill(HIST("hNSigmaPosProtonFromLambda"), posDaughterTrack.tpcNSigmaPr(), posDaughterTrack.tpcInnerParam());
            rV0ParametersData.fill(HIST("hNSigmaNegPionFromLambda"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
            // Filling the five Lambda invariant mass plots for different cuts (which are taken from namespace), same as with Kzeros above,for full explanation see the first kzero cut filling in the MC process
            for (uint32_t j = 0; j < cuthistoslambda::cosPAcuts.size(); j++) {
              std::string cosPAcutlambda = cuthistoslambda::cosPAcuts[j];
              size_t pos = cosPAcutlambda.find("_");
              cosPAcutlambda[pos] = '.';
              const float cosPAcutlambdavalue = std::stod(cosPAcutlambda);
              if (v0.v0cosPA() > cosPAcutlambdavalue) {
                cuthistoslambda::cosPACut[j]->Fill(v0.mLambda());
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
          // Filling the five Anti-Lambda invariant mass plots for different cuts (which are taken from namespace), same as with Kzeros and Lambdas above,for full explanation see the first kzero cut filling in the MC process
          if (std::abs(negDaughterTrack.tpcNSigmaPr()) < nSigmaTPCProton && std::abs(posDaughterTrack.tpcNSigmaPi()) < nSigmaTPCPion) { // TPC PID on daughter pion and proton for AntiLambda
            rV0ParametersData.fill(HIST("hMassAntiLambdaAfterAllCuts"), v0.mAntiLambda());
            rV0ParametersData.fill(HIST("hNSigmaPosPionFromAntilambda"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tpcInnerParam());
            rV0ParametersData.fill(HIST("hNSigmaNegProtonFromAntilambda"), negDaughterTrack.tpcNSigmaPr(), negDaughterTrack.tpcInnerParam());
            for (uint32_t j = 0; j < cuthistosantilambda::cosPAcuts.size(); j++) {
              std::string cosPAcutantilambda = cuthistosantilambda::cosPAcuts[j];
              size_t pos = cosPAcutantilambda.find("_");
              cosPAcutantilambda[pos] = '.';
              const float cosPAcutantilambdavalue = std::stod(cosPAcutantilambda);
              if (v0.v0cosPA() > cosPAcutantilambdavalue) {
                cuthistosantilambda::cosPACut[j]->Fill(v0.mAntiLambda());
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
  PROCESS_SWITCH(v0topologicalcuts, recMCProcess, "Process Run 3 MC:Reconstructed", true);
  PROCESS_SWITCH(v0topologicalcuts, dataProcess, "Process Run 3 Data,", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0topologicalcuts>(cfgc)};
}
