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
This task creates 20 histograms (for each of the 5 different V0 topological cuts, namely cosPointingAngle,
DCA[between]V0daughters, v0radius,DCA-positive[daughter]to-primary-vertex and DCA-negative[daughter]to-primary-vertex)
that are filled with the V0 invariant mass under the K0, Lambda and Antilambda mass assumption
(so 20cutsx5parametersx3particles=300 mass invariant plots).It also produces plots of the topological parameters themselves.
The cuts are passed as configurable strings for convenience.
This analysis includes two processes, one for Real Data and one for MC Data switchable at the end of the code, only run one at a time*/

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"

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

  // Configurable strings for Kzero cuts
  Configurable<std::string> kzeroshsetting_cospacuts_string{"kzerosetting_cospacuts", {"0_98,0_981,0_982,0_983,0_984,0_985,0_986,0_987,0_988,0_989,0_99,0_991,0_992,0_993,0_994,0_995,0_996,0_997,0_998,0_999"}, "Kzero cosPA Cut Values"};
  Configurable<std::string> kzeroshsetting_dcacuts_string{"kzerosetting_dcacuts", {"0_3,0_285,0_27,0_255,0_24,0_225,0_21,0_195,0_18,0_165,0_15,0_135,0_12,0_105,0_09,0_075,0_06,0_045,0_03,0_015"}, "Kzero DCA Cut Values"};
  Configurable<std::string> kzeroshsetting_v0radius_string{"kzerosetting_v0radiuscuts", {"0_5,0_51,0_52,0_53,0_54,0_55,0_56,0_57,0_58,0_59,0_6,0_61,0_62,0_63,0_64,0_65,0_66,0_67,0_68,0_69"}, "Kzero V0Radius Cut Values"};
  Configurable<std::string> kzeroshsetting_dcapostopv_string{"kzerosetting_dcapostopvcuts", {"0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Kzero DCA Pos to PV Cut Values"};
  Configurable<std::string> kzeroshsetting_dcanegtopv_string{"kzerosetting_dcanegtopvcuts", {"0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "KzeroDCA Neg to PV Cut Values"};

  // Configurable strings for Lambdacuts
  Configurable<std::string> lambdasetting_cospacuts_string{"lambdasetting_cospacuts", {"0_98,0_981,0_982,0_983,0_984,0_985,0_986,0_987,0_988,0_989,0_99,0_991,0_992,0_993,0_994,0_995,0_996,0_997,0_998,0_999"}, "Lambda cosPA Cut Values"};
  Configurable<std::string> lambdasetting_dcacuts_string{"lambdasetting_dcacuts", {"0_3,0_285,0_27,0_255,0_24,0_225,0_21,0_195,0_18,0_165,0_15,0_135,0_12,0_105,0_09,0_075,0_06,0_045,0_03,0_015"}, "Lambda DCA Cut Values"};
  Configurable<std::string> lambdasetting_v0radius_string{"lambdasetting_v0radiuscuts", {"0_5,0_51,0_52,0_53,0_54,0_55,0_56,0_57,0_58,0_59,0_6,0_61,0_62,0_63,0_64,0_65,0_66,0_67,0_68,0_69"}, "Lambda V0Radius Cut Values"};
  Configurable<std::string> lambdasetting_dcapostopv_string{"lambdasetting_dcapostopvcuts", {"0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Lambda DCA Pos to PV Cut Values"};
  Configurable<std::string> lambdasetting_dcanegtopv_string{"lambdasetting_dcanegtopvcuts", {"0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Lambda DCA Neg to PV Cut Values"};

  // Configurable strings for AntiLambdacuts
  Configurable<std::string> antilambdasetting_cospacuts_string{"antilambdasetting_cospacuts", {"0_98,0_981,0_982,0_983,0_984,0_985,0_986,0_987,0_988,0_989,0_99,0_991,0_992,0_993,0_994,0_995,0_996,0_997,0_998,0_999"}, "Antilambda cosPA Cut Values"};
  Configurable<std::string> antilambdasetting_dcacuts_string{"antilambdasetting_dcacuts", {"0_3,0_285,0_27,0_255,0_24,0_225,0_21,0_195,0_18,0_165,0_15,0_135,0_12,0_105,0_09,0_075,0_06,0_045,0_03,0_015"}, "Antilambda DCA Cut Values"};
  Configurable<std::string> antilambdasetting_v0radius_string{"antilambdasetting_v0radiuscuts", {"0_5,0_51,0_52,0_53,0_54,0_55,0_56,0_57,0_58,0_59,0_6,0_61,0_62,0_63,0_64,0_65,0_66,0_67,0_68,0_69"}, "Antilambda V0Radius Cut Values"};
  Configurable<std::string> antilambdasetting_dcapostopv_string{"antilambdasetting_dcapostopvcuts", {"0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Antilambda DCA Pos to PV Cut Values"};
  Configurable<std::string> antilambdasetting_dcanegtopv_string{"antilambdasetting_dcanegtopvcuts", {"0,0_01,0_02,0_03,0_04,0_05,0_06,0_07,0_08,0_09,0_1,0_11,0_12,0_13,0_14,0_15,0_16,0_17,0_18,0_19"}, "Antilambda DCA Neg to PV Cut Values"};

  void init(InitContext const&)
  {
    // kzero filling namespace with configurable strings

    // setting strings from configurable strings in order to manipulate them
    size_t commapos = 0;
    std::string token1;
    std::string kzeroshsetting_cospacuts = kzeroshsetting_cospacuts_string;
    std::string kzeroshsetting_dcacuts = kzeroshsetting_dcacuts_string;
    std::string kzeroshsetting_v0radiuscuts = kzeroshsetting_v0radius_string;
    std::string kzeroshsetting_dcapostopvcuts = kzeroshsetting_dcapostopv_string;
    std::string kzeroshsetting_dcanegtopvcuts = kzeroshsetting_dcanegtopv_string;

    // getting the  cut values for the names of the plots for the five topological cuts
    for (int i = 0; i < 20; i++) {
      commapos = kzeroshsetting_cospacuts.find(",");         // find comma that separates the values in the string
      token1 = kzeroshsetting_cospacuts.substr(0, commapos); // store the substring (first individual value)
      cuthistoskzerosh::cospacuts.push_back(token1);         //  fill the namespace with the value
      kzeroshsetting_cospacuts.erase(0, commapos + 1);       // erase the value from the set string so it moves to the next
    }
    for (int i = 0; i < 20; i++) {
      commapos = kzeroshsetting_dcacuts.find(",");
      token1 = kzeroshsetting_dcacuts.substr(0, commapos);
      cuthistoskzerosh::dcacuts.push_back(token1);
      kzeroshsetting_dcacuts.erase(0, commapos + 1);
    }
    for (int i = 0; i < 20; i++) {
      commapos = kzeroshsetting_v0radiuscuts.find(",");
      token1 = kzeroshsetting_v0radiuscuts.substr(0, commapos);
      cuthistoskzerosh::v0radiuscuts.push_back(token1);
      kzeroshsetting_v0radiuscuts.erase(0, commapos + 1);
    }
    for (int i = 0; i < 20; i++) {
      commapos = kzeroshsetting_dcapostopvcuts.find(",");
      token1 = kzeroshsetting_dcapostopvcuts.substr(0, commapos);
      cuthistoskzerosh::dcapostopvcuts.push_back(token1);
      kzeroshsetting_dcapostopvcuts.erase(0, commapos + 1);
    }
    for (int i = 0; i < 20; i++) {
      commapos = kzeroshsetting_dcanegtopvcuts.find(",");
      token1 = kzeroshsetting_dcanegtopvcuts.substr(0, commapos);
      cuthistoskzerosh::dcanegtopvcuts.push_back(token1);
      kzeroshsetting_dcanegtopvcuts.erase(0, commapos + 1);
    }

    // lambda filling namespace with configurable strings (same as in Kzeros above)
    std::string lambdasetting_cospacuts = lambdasetting_cospacuts_string;
    std::string lambdasetting_dcacuts = lambdasetting_dcacuts_string;
    std::string lambdasetting_v0radiuscuts = lambdasetting_v0radius_string;
    std::string lambdasetting_dcapostopvcuts = lambdasetting_dcapostopv_string;
    std::string lambdasetting_dcanegtopvcuts = lambdasetting_dcanegtopv_string;

    for (int i = 0; i < 20; i++) {
      commapos = lambdasetting_cospacuts.find(",");
      token1 = lambdasetting_cospacuts.substr(0, commapos);
      cuthistoslambda::cospacuts.push_back(token1);
      lambdasetting_cospacuts.erase(0, commapos + 1);
    }
    for (int i = 0; i < 20; i++) {
      commapos = lambdasetting_dcacuts.find(",");
      token1 = lambdasetting_dcacuts.substr(0, commapos);
      cuthistoslambda::dcacuts.push_back(token1);
      lambdasetting_dcacuts.erase(0, commapos + 1);
    }
    for (int i = 0; i < 20; i++) {
      commapos = lambdasetting_v0radiuscuts.find(",");
      token1 = lambdasetting_v0radiuscuts.substr(0, commapos);
      cuthistoslambda::v0radiuscuts.push_back(token1);
      lambdasetting_v0radiuscuts.erase(0, commapos + 1);
    }
    for (int i = 0; i < 20; i++) {
      commapos = lambdasetting_dcapostopvcuts.find(",");
      token1 = lambdasetting_dcapostopvcuts.substr(0, commapos);
      cuthistoslambda::dcapostopvcuts.push_back(token1);
      lambdasetting_dcapostopvcuts.erase(0, commapos + 1);
    }
    for (int i = 0; i < 20; i++) {
      commapos = lambdasetting_dcanegtopvcuts.find(",");
      token1 = lambdasetting_dcanegtopvcuts.substr(0, commapos);
      cuthistoslambda::dcanegtopvcuts.push_back(token1);
      lambdasetting_dcanegtopvcuts.erase(0, commapos + 1);
    }
    // antilambda filling namespace with configurable strings (same as in Lambdas and Kzeros above)
    std::string antilambdasetting_cospacuts = antilambdasetting_cospacuts_string;
    std::string antilambdasetting_dcacuts = antilambdasetting_dcacuts_string;
    std::string antilambdasetting_v0radiuscuts = antilambdasetting_v0radius_string;
    std::string antilambdasetting_dcapostopvcuts = antilambdasetting_dcapostopv_string;
    std::string antilambdasetting_dcanegtopvcuts = antilambdasetting_dcanegtopv_string;

    for (int i = 0; i < 20; i++) {
      commapos = antilambdasetting_cospacuts.find(",");
      token1 = antilambdasetting_cospacuts.substr(0, commapos);
      cuthistosantilambda::cospacuts.push_back(token1);
      antilambdasetting_cospacuts.erase(0, commapos + 1);
    }
    for (int i = 0; i < 20; i++) {
      commapos = antilambdasetting_dcacuts.find(",");
      token1 = antilambdasetting_dcacuts.substr(0, commapos);
      cuthistosantilambda::dcacuts.push_back(token1);
      antilambdasetting_dcacuts.erase(0, commapos + 1);
    }
    for (int i = 0; i < 20; i++) {
      commapos = antilambdasetting_v0radiuscuts.find(",");
      token1 = antilambdasetting_v0radiuscuts.substr(0, commapos);
      cuthistosantilambda::v0radiuscuts.push_back(token1);
      antilambdasetting_v0radiuscuts.erase(0, commapos + 1);
    }
    for (int i = 0; i < 20; i++) {
      commapos = antilambdasetting_dcapostopvcuts.find(",");
      token1 = antilambdasetting_dcapostopvcuts.substr(0, commapos);
      cuthistosantilambda::dcapostopvcuts.push_back(token1);
      antilambdasetting_dcapostopvcuts.erase(0, commapos + 1);
    }
    for (int i = 0; i < 20; i++) {
      commapos = antilambdasetting_dcanegtopvcuts.find(",");
      token1 = antilambdasetting_dcanegtopvcuts.substr(0, commapos);
      cuthistosantilambda::dcanegtopvcuts.push_back(token1);
      antilambdasetting_dcanegtopvcuts.erase(0, commapos + 1);
    }

    // Axes for the three invariant mass plots
    AxisSpec K0ShortMassAxis = {nBins, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec LambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec AntiLambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};

    // adding the invariant mass histograms to their Registries using the namespace for kzeros, lambdas and antilambdas
    for (int i = 0; i < 20; i++) {
      cuthistoskzerosh::cospaCut[i] = rKzeroShort_cospaCut.add<TH1>(fmt::format("hKzerocospaCut_{}", cuthistoskzerosh::cospacuts[i]).data(), fmt::format("hKzerocospaCut_{}", cuthistoskzerosh::cospacuts[i]).data(), {HistType::kTH1D, {{K0ShortMassAxis}}});
      cuthistoskzerosh::dcaCut[i] = rKzeroShort_dcaCut.add<TH1>(fmt::format("hKzerodcaCut_{}", cuthistoskzerosh::dcacuts[i]).data(), fmt::format("hKzerodcaCut_{}", cuthistoskzerosh::dcacuts[i]).data(), {HistType::kTH1D, {{K0ShortMassAxis}}});
      cuthistoskzerosh::v0radiusCut[i] = rKzeroShort_v0radiusCut.add<TH1>(fmt::format("hKzerov0radiusCut_{}", cuthistoskzerosh::v0radiuscuts[i]).data(), fmt::format("hKzerov0radiusCuts_{}", cuthistoskzerosh::v0radiuscuts[i]).data(), {HistType::kTH1D, {{K0ShortMassAxis}}});
      cuthistoskzerosh::dcapostopCut[i] = rKzeroShort_dcapostopCut.add<TH1>(fmt::format("hKzerodcapostopCut_{}", cuthistoskzerosh::dcapostopvcuts[i]).data(), fmt::format("hKzerodcapostopCut_{}", cuthistoskzerosh::dcapostopvcuts[i]).data(), {HistType::kTH1D, {{K0ShortMassAxis}}});
      cuthistoskzerosh::dcanegtopCut[i] = rKzeroShort_dcanegtopCut.add<TH1>(fmt::format("hKzerodcanegtopCut_{}", cuthistoskzerosh::dcanegtopvcuts[i]).data(), fmt::format("hKzerodcanegtopCut_{}", cuthistoskzerosh::dcanegtopvcuts[i]).data(), {HistType::kTH1D, {{K0ShortMassAxis}}});
    }

    for (int i = 0; i < 20; i++) {
      cuthistoslambda::cospaCut[i] = rLambda_cospaCut.add<TH1>(fmt::format("hLambdacospaCut_{}", cuthistoslambda::cospacuts[i]).data(), fmt::format("hLambdacospaCut_{}", cuthistoslambda::cospacuts[i]).data(), {HistType::kTH1D, {{LambdaMassAxis}}});
      cuthistoslambda::dcaCut[i] = rLambda_dcaCut.add<TH1>(fmt::format("hLambdadcaCut_{}", cuthistoslambda::dcacuts[i]).data(), fmt::format("hLambdadcaCut_{}", cuthistoslambda::dcacuts[i]).data(), {HistType::kTH1D, {{LambdaMassAxis}}});
      cuthistoslambda::v0radiusCut[i] = rLambda_v0radiusCut.add<TH1>(fmt::format("hLambdav0radiusCut_{}", cuthistoslambda::v0radiuscuts[i]).data(), fmt::format("hLambdav0radiusCuts_{}", cuthistoslambda::v0radiuscuts[i]).data(), {HistType::kTH1D, {{LambdaMassAxis}}});
      cuthistoslambda::dcapostopCut[i] = rLambda_dcapostopCut.add<TH1>(fmt::format("hLambdadcapostopCut_{}", cuthistoslambda::dcapostopvcuts[i]).data(), fmt::format("hLambdadcapostopCut_{}", cuthistoslambda::dcapostopvcuts[i]).data(), {HistType::kTH1D, {{LambdaMassAxis}}});
      cuthistoslambda::dcanegtopCut[i] = rLambda_dcanegtopCut.add<TH1>(fmt::format("hLambdadcanegtopCut_{}", cuthistoslambda::dcanegtopvcuts[i]).data(), fmt::format("hLambdadcanegtopCut_{}", cuthistoslambda::dcanegtopvcuts[i]).data(), {HistType::kTH1D, {{LambdaMassAxis}}});
    }

    for (int i = 0; i < 20; i++) {
      cuthistosantilambda::cospaCut[i] = rAntiLambda_cospaCut.add<TH1>(fmt::format("hAntiLambdacospaCut_{}", cuthistosantilambda::cospacuts[i]).data(), fmt::format("hAntiLambdacospaCut_{}", cuthistosantilambda::cospacuts[i]).data(), {HistType::kTH1D, {{AntiLambdaMassAxis}}});
      cuthistosantilambda::dcaCut[i] = rAntiLambda_dcaCut.add<TH1>(fmt::format("hAntiLambdadcaCut_{}", cuthistosantilambda::dcacuts[i]).data(), fmt::format("hAntiLambdadcaCut_{}", cuthistosantilambda::dcacuts[i]).data(), {HistType::kTH1D, {{AntiLambdaMassAxis}}});
      cuthistosantilambda::v0radiusCut[i] = rAntiLambda_v0radiusCut.add<TH1>(fmt::format("hAntiLambdav0radiusCut_{}", cuthistosantilambda::v0radiuscuts[i]).data(), fmt::format("hAntiLambdav0radiusCuts_{}", cuthistosantilambda::v0radiuscuts[i]).data(), {HistType::kTH1D, {{AntiLambdaMassAxis}}});
      cuthistosantilambda::dcapostopCut[i] = rAntiLambda_dcapostopCut.add<TH1>(fmt::format("hAntiLambdadcapostopCut_{}", cuthistosantilambda::dcapostopvcuts[i]).data(), fmt::format("hAntiLambdadcapostopCut_{}", cuthistosantilambda::dcapostopvcuts[i]).data(), {HistType::kTH1D, {{AntiLambdaMassAxis}}});
      cuthistosantilambda::dcanegtopCut[i] = rAntiLambda_dcanegtopCut.add<TH1>(fmt::format("hAntiLambdadcanegtopCut_{}", cuthistosantilambda::dcanegtopvcuts[i]).data(), fmt::format("hAntiLambdadcanegtopCut_{}", cuthistosantilambda::dcanegtopvcuts[i]).data(), {HistType::kTH1D, {{AntiLambdaMassAxis}}});
    }

    // K0s topological/PID cut histograms added and MC-matched
    rV0Parameters_MC_V0match.add("hDCAV0Daughters_V0_Match", "hDCAV0Daughters_No_Match", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0Parameters_MC_V0match.add("hV0CosPA_V0_Match", "hV0CosPA_No_Match", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0Parameters_MC_V0match.add("hV0Radius_V0_Match", "hV0Radius_No_Match", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});
    rV0Parameters_MC_V0match.add("hV0Radius_V0_Match_Full", "hV0Radius_No_Match_Full", {HistType::kTH1F, {{nBins, 0.0f, 40.0f}}});
    rV0Parameters_MC_V0match.add("hDCAPostoPV_V0_Match", "hDCAPostoPV_No_Match", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});
    rV0Parameters_MC_V0match.add("hDCANegtoPV_V0_Match", "hDCANegtoPV_No_Match", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});

    // kzero match
    rV0Parameters_MC_K0Smatch.add("hDCAV0Daughters_KzeroMC_Match", "hDCAV0Daughters_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0Parameters_MC_K0Smatch.add("hV0CosPA_KzeroMC_Match", "hV0CosPA_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0Parameters_MC_K0Smatch.add("hV0Radius_KzeroMC_Match", "hV0Radius_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.2f, 5.0f}}});
    rV0Parameters_MC_K0Smatch.add("hDCAPostoPV_KzeroMC_Match", "hDCAPostoPV_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});
    rV0Parameters_MC_K0Smatch.add("hDCANegtoPV_KzeroMC_Match", "hDCANegtoPV_KzeroMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});

    // lambda match
    rV0Parameters_MC_Lambdamatch.add("hDCAV0Daughters_LambdaMC_Match", "hDCAV0Daughters_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0Parameters_MC_Lambdamatch.add("hV0CosPA_LambdaMC_Match", "hV0CosPA_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0Parameters_MC_Lambdamatch.add("hV0Radius_LambdaMC_Match", "hV0Radius_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.2f, 5.0f}}});
    rV0Parameters_MC_Lambdamatch.add("hDCAPostoPV_LambdaMC_Match", "hDCAPostoPV_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});
    rV0Parameters_MC_Lambdamatch.add("hDCANegtoPV_LambdaMC_Match", "hDCANegtoPV_LambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});

    // antilambda match
    rV0Parameters_MC_AntiLambdamatch.add("hDCAV0Daughters_AntiLambdaMC_Match", "hDCAV0Daughters_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0Parameters_MC_AntiLambdamatch.add("hV0CosPA_AntiLambdaMC_Match", "hV0CosPA_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0Parameters_MC_AntiLambdamatch.add("hV0Radius_AntiLambdaMC_Match", "hV0Radius_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.2f, 5.0f}}});
    rV0Parameters_MC_AntiLambdamatch.add("hDCAPostoPV_AntiLambdaMC_Match", "hDCANegtoPV_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});
    rV0Parameters_MC_AntiLambdamatch.add("hDCANegtoPV_AntiLambdaMC_Match", "hDCANegtoPV_AntiLambdaMC_Match", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});

    // V0s Data
    rV0Parameters_Data.add("hDCAV0Daughters_V0_Data", "hDCAV0Daughters_V0_Data", {HistType::kTH1F, {{nBins, 0.0f, 1.2f}}});
    rV0Parameters_Data.add("hV0CosPA_V0_Data", "hV0CosPA_V0_Data", {HistType::kTH1F, {{nBins, 0.95f, 1.f}}});
    rV0Parameters_Data.add("hV0Radius_V0_Data", "hV0Radius_V0_Data", {HistType::kTH1F, {{nBins, 0.2f, 5.0f}}});
    rV0Parameters_Data.add("hDCAPostoPV_V0_Data", "hDCAPostoPV_V0_Data", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});
    rV0Parameters_Data.add("hDCANegtoPV_V0_Data", "hDCANegtoPV_V0_Data", {HistType::kTH1F, {{nBins, 0.0f, 5.0f}}});
    rV0Parameters_Data.add("hMassK0ShortNoCuts_V0_Data", "hMassK0ShortNoCuts_V0_Data", {HistType::kTH1F, {{K0ShortMassAxis}}});
    rV0Parameters_Data.add("hMassLambdaNoCuts_V0_Data", "hMassLambdaNoCuts_V0_Data", {HistType::kTH1F, {{LambdaMassAxis}}});
    rV0Parameters_Data.add("hMassAntilambdaNoCuts_V0_Data", "hMassAntilambdaNoCuts_V0_Data", {HistType::kTH1F, {{AntiLambdaMassAxis}}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::McTrackLabels>;

  // This is the Process for the MC reconstructed Data
  void RecMCprocess(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const&,
                    soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s,
                    DaughterTracks const&, // no need to define a variable for tracks, if we don't access them directly
                    aod::McParticles const&)
  {
    for (const auto& v0 : V0s) {

      // filling histograms with V0 values
      rV0Parameters_MC_V0match.fill(HIST("hDCAV0Daughters_V0_Match"), v0.dcaV0daughters());
      rV0Parameters_MC_V0match.fill(HIST("hV0CosPA_V0_Match"), v0.v0cosPA());
      rV0Parameters_MC_V0match.fill(HIST("hV0Radius_V0_Match"), v0.v0radius());
      rV0Parameters_MC_V0match.fill(HIST("hV0Radius_V0_Match_Full"), v0.v0radius());
      rV0Parameters_MC_V0match.fill(HIST("hDCAPostoPV_V0_Match"), v0.dcapostopv());
      rV0Parameters_MC_V0match.fill(HIST("hDCANegtoPV_V0_Match"), v0.dcanegtopv());

      // Checking that the V0 is a true K0s/Lambdas/Antilambdas and then filling the parameter histograms and the invariant mass plots for different cuts (which are taken from namespace)
      if (v0.has_mcParticle()) {
        auto v0mcParticle = v0.mcParticle();
        if (v0mcParticle.pdgCode() == 310) { // kzero matched

          rV0Parameters_MC_K0Smatch.fill(HIST("hDCAV0Daughters_KzeroMC_Match"), v0.dcaV0daughters());
          rV0Parameters_MC_K0Smatch.fill(HIST("hV0CosPA_KzeroMC_Match"), v0.v0cosPA());
          rV0Parameters_MC_K0Smatch.fill(HIST("hV0Radius_KzeroMC_Match"), v0.v0radius());
          rV0Parameters_MC_K0Smatch.fill(HIST("hDCAPostoPV_KzeroMC_Match"), v0.dcapostopv());
          rV0Parameters_MC_K0Smatch.fill(HIST("hDCANegtoPV_KzeroMC_Match"), v0.dcanegtopv());

          for (int j = 0; j < 20; j++) {
            std::string cospacut = cuthistoskzerosh::cospacuts[j]; // Get the current cut value from the namespace
            size_t pos = cospacut.find("_");                       // find the "_" which needs to change to a "." for it to be a number
            cospacut[pos] = '.';                                   // change the "_" into an "."
            const float cospacutvalue = std::stod(cospacut);       // make the string into a float value
            if (v0.v0cosPA() > cospacutvalue) {                    // enforce the cut value
              cuthistoskzerosh::cospaCut[j]->Fill(v0.mK0Short());  // fill the corresponding histo from the namespace with the invariant mass (of a Kzero here)
            }
          }
          for (int j = 0; j < 20; j++) {
            std::string dcacut = cuthistoskzerosh::dcacuts[j];
            size_t pos = dcacut.find("_");
            dcacut[pos] = '.';
            const float dcacutvalue = std::stod(dcacut);
            if (v0.dcaV0daughters() < dcacutvalue) {
              cuthistoskzerosh::dcaCut[j]->Fill(v0.mK0Short());
            }
          }
          for (int j = 0; j < 20; j++) {
            std::string v0radiuscut = cuthistoskzerosh::v0radiuscuts[j];
            size_t pos = v0radiuscut.find("_");
            v0radiuscut[pos] = '.';
            const float v0radiuscutvalue = std::stod(v0radiuscut);
            if (v0.v0radius() > v0radiuscutvalue) {
              cuthistoskzerosh::v0radiusCut[j]->Fill(v0.mK0Short());
            }
          }
          for (int j = 0; j < 20; j++) {
            std::string dcapostopcut = cuthistoskzerosh::dcapostopvcuts[j];
            size_t pos = dcapostopcut.find("_");
            dcapostopcut[pos] = '.';
            const float dcapostopcutvalue = std::stod(dcapostopcut);
            if (v0.dcapostopv() > dcapostopcutvalue) {
              cuthistoskzerosh::dcapostopCut[j]->Fill(v0.mK0Short());
            }
          }
          for (int j = 0; j < 20; j++) {
            std::string dcanegtopcut = cuthistoskzerosh::dcanegtopvcuts[j];
            size_t pos = dcanegtopcut.find("_");
            dcanegtopcut[pos] = '.';
            const float dcanegtopcutvalue = std::stod(dcanegtopcut);
            if (v0.dcanegtopv() > dcanegtopcutvalue) {
              cuthistoskzerosh::dcanegtopCut[j]->Fill(v0.mK0Short());
            }
          }
        }
        if (v0mcParticle.pdgCode() == 3122) { // lambda matched
          rV0Parameters_MC_Lambdamatch.fill(HIST("hDCAV0Daughters_LambdaMC_Match"), v0.dcaV0daughters());
          rV0Parameters_MC_Lambdamatch.fill(HIST("hV0CosPA_LambdaMC_Match"), v0.v0cosPA());
          rV0Parameters_MC_Lambdamatch.fill(HIST("hV0Radius_LambdaMC_Match"), v0.v0radius());
          rV0Parameters_MC_Lambdamatch.fill(HIST("hDCAPostoPV_LambdaMC_Match"), v0.dcapostopv());
          rV0Parameters_MC_Lambdamatch.fill(HIST("hDCANegtoPV_LambdaMC_Match"), v0.dcanegtopv());

          // for explanation look at the first Kzero  plot above
          for (int j = 0; j < 20; j++) {
            std::string cospacutlambda = cuthistoslambda::cospacuts[j];
            size_t pos = cospacutlambda.find("_");
            cospacutlambda[pos] = '.';
            const float cospacutlambdavalue = std::stod(cospacutlambda);
            if (v0.v0cosPA() > cospacutlambdavalue) {
              cuthistoslambda::cospaCut[j]->Fill(v0.mLambda());
            }
          }
          for (int j = 0; j < 20; j++) {
            std::string dcacutlambda = cuthistoslambda::dcacuts[j];
            size_t pos = dcacutlambda.find("_");
            dcacutlambda[pos] = '.';
            const float dcacutlambdavalue = std::stod(dcacutlambda);
            if (v0.dcaV0daughters() < dcacutlambdavalue) {
              cuthistoslambda::dcaCut[j]->Fill(v0.mLambda());
            }
          }
          for (int j = 0; j < 20; j++) {
            std::string v0radiuscutlambda = cuthistoslambda::v0radiuscuts[j];
            size_t pos = v0radiuscutlambda.find("_");
            v0radiuscutlambda[pos] = '.';
            const float v0radiuscutlambdavalue = std::stod(v0radiuscutlambda);
            if (v0.v0radius() > v0radiuscutlambdavalue) {
              cuthistoslambda::v0radiusCut[j]->Fill(v0.mLambda());
            }
          }
          for (int j = 0; j < 20; j++) {
            std::string dcapostopcutlambda = cuthistoslambda::dcapostopvcuts[j];
            size_t pos = dcapostopcutlambda.find("_");
            dcapostopcutlambda[pos] = '.';
            const float dcapostopcutlambdavalue = std::stod(dcapostopcutlambda);
            if (v0.dcapostopv() > dcapostopcutlambdavalue) {
              cuthistoslambda::dcapostopCut[j]->Fill(v0.mLambda());
            }
          }
          for (int j = 0; j < 20; j++) {
            std::string dcanegtopcutlambda = cuthistoslambda::dcanegtopvcuts[j];
            size_t pos = dcanegtopcutlambda.find("_");
            dcanegtopcutlambda[pos] = '.';
            const float dcanegtopcutlambdavalue = std::stod(dcanegtopcutlambda);
            if (v0.dcanegtopv() > dcanegtopcutlambdavalue) {
              cuthistoslambda::dcanegtopCut[j]->Fill(v0.mLambda());
            }
          }
        }
        if (v0mcParticle.pdgCode() == -3122) { // antilambda matched
          rV0Parameters_MC_AntiLambdamatch.fill(HIST("hDCAV0Daughters_AntiLambdaMC_Match"), v0.dcaV0daughters());
          rV0Parameters_MC_AntiLambdamatch.fill(HIST("hV0CosPA_AntiLambdaMC_Match"), v0.v0cosPA());
          rV0Parameters_MC_AntiLambdamatch.fill(HIST("hV0Radius_AntiLambdaMC_Match"), v0.v0radius());
          rV0Parameters_MC_AntiLambdamatch.fill(HIST("hDCAPostoPV_AntiLambdaMC_Match"), v0.dcapostopv());
          rV0Parameters_MC_AntiLambdamatch.fill(HIST("hDCANegtoPV_AntiLambdaMC_Match"), v0.dcanegtopv());
          // for explanation look at the first Kzero  plot above
          for (int j = 0; j < 20; j++) {
            std::string cospacutantilambda = cuthistosantilambda::cospacuts[j];
            size_t pos = cospacutantilambda.find("_");
            cospacutantilambda[pos] = '.';
            const float cospacutantilambdavalue = std::stod(cospacutantilambda);
            if (v0.v0cosPA() > cospacutantilambdavalue) {
              cuthistosantilambda::cospaCut[j]->Fill(v0.mAntiLambda());
            }
          }
          for (int j = 0; j < 20; j++) {
            std::string dcacutantilambda = cuthistosantilambda::dcacuts[j];
            size_t pos = dcacutantilambda.find("_");
            dcacutantilambda[pos] = '.';
            const float dcacutantilambdavalue = std::stod(dcacutantilambda);
            if (v0.dcaV0daughters() < dcacutantilambdavalue) {
              cuthistosantilambda::dcaCut[j]->Fill(v0.mAntiLambda());
            }
          }
          for (int j = 0; j < 20; j++) {
            std::string v0radiusantilambda = cuthistosantilambda::v0radiuscuts[j];
            size_t pos = v0radiusantilambda.find("_");
            v0radiusantilambda[pos] = '.';
            const float v0radiuscutantilambdavalue = std::stod(v0radiusantilambda);
            if (v0.v0radius() > v0radiuscutantilambdavalue) {
              cuthistosantilambda::v0radiusCut[j]->Fill(v0.mAntiLambda());
            }
          }
          for (int j = 0; j < 20; j++) {
            std::string dcapostopantilambda = cuthistosantilambda::dcapostopvcuts[j];
            size_t pos = dcapostopantilambda.find("_");
            dcapostopantilambda[pos] = '.';
            const float dcapostopcutantilambdavalue = std::stod(dcapostopantilambda);
            if (v0.dcapostopv() > dcapostopcutantilambdavalue) {
              cuthistosantilambda::dcapostopCut[j]->Fill(v0.mAntiLambda());
            }
          }
          for (int j = 0; j < 20; j++) {
            std::string dcanegtopantilambda = cuthistosantilambda::dcanegtopvcuts[j];
            size_t pos = dcanegtopantilambda.find("_");
            dcanegtopantilambda[pos] = '.';
            const float dcanegtopcutantilambdavalue = std::stod(dcanegtopantilambda);
            if (v0.dcanegtopv() > dcanegtopcutantilambdavalue) {
              cuthistosantilambda::dcanegtopCut[j]->Fill(v0.mAntiLambda());
            }
          }
        }
      }
    }
  }
  // This is the process for Real Data
  void Dataprocess(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const&,
                   aod::V0Datas const& V0s)
  {
    // filling histograms with the different V0 parameters
    for (const auto& v0 : V0s) {
      rV0Parameters_Data.fill(HIST("hMassK0ShortNoCuts_V0_Data"), v0.mK0Short());
      rV0Parameters_Data.fill(HIST("hMassLambdaNoCuts_V0_Data"), v0.mLambda());
      rV0Parameters_Data.fill(HIST("hMassAntiLambdaNoCuts_V0_Data"), v0.mAntiLambda());
      rV0Parameters_Data.fill(HIST("hDCAV0Daughters_V0_Data"), v0.dcaV0daughters());
      rV0Parameters_Data.fill(HIST("hV0CosPA_V0_Data"), v0.v0cosPA());
      rV0Parameters_Data.fill(HIST("hV0Radius_V0_Data"), v0.v0radius());
      rV0Parameters_Data.fill(HIST("hV0Radius_Full_V0_Data"), v0.v0radius());
      rV0Parameters_Data.fill(HIST("hDCAPostoPV_V0_Data"), v0.dcapostopv());
      rV0Parameters_Data.fill(HIST("hDCANegtoPV_V0_Data"), v0.dcanegtopv());

      // Filling the five Kzero invariant mass plots for different cuts (which are taken from namespace), for full explanation see the first kzero cut filling in the MC process
      for (int j = 0; j < 20; j++) {
        std::string cospacut = cuthistoskzerosh::cospacuts[j];
        size_t pos = cospacut.find("_");
        cospacut[pos] = '.';
        const float cospacutvalue = std::stod(cospacut);
        if (v0.v0cosPA() > cospacutvalue) {
          cuthistoskzerosh::cospaCut[j]->Fill(v0.mK0Short());
        }
      }
      for (int j = 0; j < 20; j++) {
        std::string dcacut = cuthistoskzerosh::dcacuts[j];
        size_t pos = dcacut.find("_");
        dcacut[pos] = '.';
        const float dcacutvalue = std::stod(dcacut);
        if (v0.dcaV0daughters() < dcacutvalue) {
          cuthistoskzerosh::dcaCut[j]->Fill(v0.mK0Short());
        }
      }
      for (int j = 0; j < 20; j++) {
        std::string v0radiuscut = cuthistoskzerosh::v0radiuscuts[j];
        size_t pos = v0radiuscut.find("_");
        v0radiuscut[pos] = '.';
        const float v0radiuscutvalue = std::stod(v0radiuscut);
        if (v0.v0radius() > v0radiuscutvalue) {
          cuthistoskzerosh::v0radiusCut[j]->Fill(v0.mK0Short());
        }
      }
      for (int j = 0; j < 20; j++) {
        std::string dcapostopcut = cuthistoskzerosh::dcapostopvcuts[j];
        size_t pos = dcapostopcut.find("_");
        dcapostopcut[pos] = '.';
        const float dcapostopcutvalue = std::stod(dcapostopcut);
        if (v0.dcapostopv() > dcapostopcutvalue) {
          cuthistoskzerosh::dcapostopCut[j]->Fill(v0.mK0Short());
        }
      }
      for (int j = 0; j < 20; j++) {
        std::string dcanegtopcut = cuthistoskzerosh::dcanegtopvcuts[j];
        size_t pos = dcanegtopcut.find("_");
        dcanegtopcut[pos] = '.';
        const float dcanegtopcutvalue = std::stod(dcanegtopcut);
        if (v0.dcanegtopv() > dcanegtopcutvalue) {
          cuthistoskzerosh::dcanegtopCut[j]->Fill(v0.mK0Short());
        }
      }
      // Filling the five Lambda invariant mass plots for different cuts (which are taken from namespace), same as with Kzeros above,for full explanation see the first kzero cut filling in the MC process
      for (int j = 0; j < 20; j++) {
        std::string cospacutlambda = cuthistoslambda::cospacuts[j];
        size_t pos = cospacutlambda.find("_");
        cospacutlambda[pos] = '.';
        const float cospacutlambdavalue = std::stod(cospacutlambda);
        if (v0.v0cosPA() > cospacutlambdavalue) {
          cuthistoslambda::cospaCut[j]->Fill(v0.mLambda());
        }
      }
      for (int j = 0; j < 20; j++) {
        std::string dcacutlambda = cuthistoslambda::dcacuts[j];
        size_t pos = dcacutlambda.find("_");
        dcacutlambda[pos] = '.';
        const float dcacutlambdavalue = std::stod(dcacutlambda);
        if (v0.dcaV0daughters() < dcacutlambdavalue) {
          cuthistoslambda::dcaCut[j]->Fill(v0.mLambda());
        }
      }
      for (int j = 0; j < 20; j++) {
        std::string v0radiuscutlambda = cuthistoslambda::v0radiuscuts[j];
        size_t pos = v0radiuscutlambda.find("_");
        v0radiuscutlambda[pos] = '.';
        const float v0radiuscutlambdavalue = std::stod(v0radiuscutlambda);
        if (v0.v0radius() > v0radiuscutlambdavalue) {
          cuthistoslambda::v0radiusCut[j]->Fill(v0.mLambda());
        }
      }
      for (int j = 0; j < 20; j++) {
        std::string dcapostopcutlambda = cuthistoslambda::dcapostopvcuts[j];
        size_t pos = dcapostopcutlambda.find("_");
        dcapostopcutlambda[pos] = '.';
        const float dcapostopcutlambdavalue = std::stod(dcapostopcutlambda);
        if (v0.dcapostopv() > dcapostopcutlambdavalue) {
          cuthistoslambda::dcapostopCut[j]->Fill(v0.mLambda());
        }
      }
      for (int j = 0; j < 20; j++) {
        std::string dcanegtopcutlambda = cuthistoslambda::dcanegtopvcuts[j];
        size_t pos = dcanegtopcutlambda.find("_");
        dcanegtopcutlambda[pos] = '.';
        const float dcanegtopcutlambdavalue = std::stod(dcanegtopcutlambda);
        if (v0.dcanegtopv() > dcanegtopcutlambdavalue) {
          cuthistoslambda::dcanegtopCut[j]->Fill(v0.mLambda());
        }
      }
      // Filling the five Lambda invariant mass plots for different cuts (which are taken from namespace), same as with Kzeros and Lambdas above,for full explanation see the first kzero cut filling in the MC process
      for (int j = 0; j < 20; j++) {
        std::string cospacutantilambda = cuthistosantilambda::cospacuts[j];
        size_t pos = cospacutantilambda.find("_");
        cospacutantilambda[pos] = '.';
        const float cospacutantilambdavalue = std::stod(cospacutantilambda);
        if (v0.v0cosPA() > cospacutantilambdavalue) {
          cuthistosantilambda::cospaCut[j]->Fill(v0.mAntiLambda());
        }
      }
      for (int j = 0; j < 20; j++) {
        std::string dcacutantilambda = cuthistosantilambda::dcacuts[j];
        size_t pos = dcacutantilambda.find("_");
        dcacutantilambda[pos] = '.';
        const float dcacutantilambdavalue = std::stod(dcacutantilambda);
        if (v0.dcaV0daughters() < dcacutantilambdavalue) {
          cuthistosantilambda::dcaCut[j]->Fill(v0.mAntiLambda());
        }
      }
      for (int j = 0; j < 20; j++) {
        std::string v0radiusantilambda = cuthistosantilambda::v0radiuscuts[j];
        size_t pos = v0radiusantilambda.find("_");
        v0radiusantilambda[pos] = '.';
        const float v0radiuscutantilambdavalue = std::stod(v0radiusantilambda);
        if (v0.v0radius() > v0radiuscutantilambdavalue) {
          cuthistosantilambda::v0radiusCut[j]->Fill(v0.mAntiLambda());
        }
      }
      for (int j = 0; j < 20; j++) {
        std::string dcapostopantilambda = cuthistosantilambda::dcapostopvcuts[j];
        size_t pos = dcapostopantilambda.find("_");
        dcapostopantilambda[pos] = '.';
        const float dcapostopcutantilambdavalue = std::stod(dcapostopantilambda);
        if (v0.dcapostopv() > dcapostopcutantilambdavalue) {
          cuthistosantilambda::dcapostopCut[j]->Fill(v0.mAntiLambda());
        }
      }
      for (int j = 0; j < 20; j++) {
        std::string dcanegtopantilambda = cuthistosantilambda::dcanegtopvcuts[j];
        size_t pos = dcanegtopantilambda.find("_");
        dcanegtopantilambda[pos] = '.';
        const float dcanegtopcutantilambdavalue = std::stod(dcanegtopantilambda);
        if (v0.dcanegtopv() > dcanegtopcutantilambdavalue) {
          cuthistosantilambda::dcanegtopCut[j]->Fill(v0.mAntiLambda());
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
