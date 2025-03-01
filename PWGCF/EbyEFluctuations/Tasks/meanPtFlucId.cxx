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

/// \file meanPtFlucId.cxx
/// \brief Calculate EbyE <pt> fluctuations with cumulant method.
///        For charged particles and identified particles.
///        For RUN-3
///
/// \author Tanu Gahlaut <tanu.gahlaut@cern.ch>

#include <utility>
#include <vector>
#include <TPDGCode.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace std;

struct MeanPtFlucId {
  Configurable<int> nPBins{"nPBins", 300, ""};
  Configurable<int> nPartBins{"nPartBins", 250, ""};
  Configurable<int> nCentBins{"nCentBins", 101, ""};
  Configurable<int> nRapBins{"nRapBins", 100, ""};
  Configurable<int> nPhiBins{"nPhiBins", 100, ""};
  Configurable<float> cfgCutPtMax{"cfgCutPtMax", 3.0, "maximum pT"};
  Configurable<float> cfgCutPtMin{"cfgCutPtMin", 0.15, "minimum pT"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut"};
  Configurable<float> cfgCutRap{"cfgCutRap", 0.5, "Rapidity Cut"};
  Configurable<float> cfgCutDcaXY{"cfgCutDcaXY", 0.12, "DCAxy cut"};
  Configurable<float> cfgCutDcaZ{"cfgCutDcaZ", 0.3, "DCAz cut"};
  Configurable<float> cfgCutPosZ{"cfgCutPosZ", 10.0, "cut for vertex Z"};
  Configurable<float> cfgGammaCut{"cfgGammaCut", 0.003, "Gamma inv Mass Cut for electron-positron rejection"};
  Configurable<float> cfgCutNSig2{"cfgCutNSig2", 2.0, "nSigma cut (2)"};
  Configurable<float> cfgCutNSig3{"cfgCutNSig3", 3.0, "nSigma cut (3)"};
  Configurable<float> cfgCutPiPtMin{"cfgCutPiPtMin", 0.2, "Minimum pion p_{T} cut"};
  Configurable<float> cfgCutKaPtMin{"cfgCutKaPtMin", 0.3, "Minimum kaon p_{T} cut"};
  Configurable<float> cfgCutPrPtMin{"cfgCutPrPtMin", 0.5, "Minimum proton p_{T} cut"};
  Configurable<float> cfgCutPiThrsldP{"cfgCutPiThrsldP", 0.6, "Threshold p cut pion"};
  Configurable<float> cfgCutKaThrsldP{"cfgCutKaThrsldP", 0.6, "Threshold p cut kaon"};
  Configurable<float> cfgCutPrThrsldP{"cfgCutPrThrsldP", 1.0, "Threshold p cut proton "};
  Configurable<float> cfgCutPiP1{"cfgCutPiP1", 0.5, "pion p cut-1"};
  Configurable<float> cfgCutPiP2{"cfgCutPiP2", 0.6, "pion p cut-2"};
  Configurable<float> cfgCutKaP1{"cfgCutKaP1", 0.4, "kaon p cut-1"};
  Configurable<float> cfgCutKaP2{"cfgCutKaP2", 0.6, "kaon p cut-2"};
  Configurable<float> cfgCutKaP3{"cfgCutKaP3", 1.2, "kaon p cut-3"};
  Configurable<float> cfgCutPrP1{"cfgCutPrP1", 0.9, "proton p cut-1"};
  Configurable<float> cfgCutPrP2{"cfgCutPrP2", 1.0, "proton p cut-2"};
  Configurable<bool> cfgRun3{"cfgRun3", true, ""};
  Configurable<bool> cfgRun2{"cfgRun2", false, ""};
  Configurable<bool> cfgCorrection{"cfgCorrection", true, "Efficiency Correction"};
  Configurable<bool> cfgCorrectionPID{"cfgCorrectionPID", true, "ID particles Efficiency Correction"};
  Configurable<bool> cfgCorrectionPtRap{"cfgCorrectionPtRap", false, "Efficiency Correction for pT and eta"};
  Configurable<bool> cfgCorrectionPtRapPID{"cfgCorrectionPtRapPID", false, "ID particles Efficiency Correction for pT and eta"};
  Configurable<bool> cfgPidCut{"cfgPidCut", false, ""};
  Configurable<bool> cfgPDGCodeOnly{"cfgPDGCodeOnly", true, ""};
  Configurable<bool> cfgMCReco{"cfgMCReco", false, ""};
  Configurable<bool> cfgMCTruth{"cfgMCTruth", false, ""};
  Configurable<bool> cfgPosZ{"cfgPosZ", true, "Position Z"};
  Configurable<bool> cfgSel7{"cfgSel7", true, "Run2 Sel7 trigger"};
  Configurable<bool> cfgkINT7{"cfgkINT7", true, "Run2 MB trigger"};
  Configurable<bool> cfgSel8{"cfgSel8", true, "Sel8 trigger"};
  Configurable<bool> cfgNoSameBunchPileup{"cfgNoSameBunchPileup", true, "kNoSameBunchPileup"};
  Configurable<bool> cfgIsVertexITSTPC{"cfgIsVertexITSTPC", true, "kIsVertexITSTPC"};
  Configurable<bool> cfgIsGoodZvtxFT0vsPV{"cfgIsGoodZvtxFT0vsPV", true, "kIsGoodZvtxFT0vsPV"};
  Configurable<bool> cfgTVXinTRD{"cfgTVXinTRD", true, "cfgTVXinTRD"};
  Configurable<bool> cfgNoCollInTimeRangeStandard{"cfgNoCollInTimeRangeStandard", true, "cfgNoCollInTimeRangeStandard"};
  Configurable<bool> cfgRejTrk{"cfgRejTrk", true, "Rejected Tracks"};
  Configurable<bool> cfgInvMass{"cfgInvMass", true, "electron Inv Mass cut selection"};
  Configurable<bool> cfgSelOR{"cfgSelOR", true, "Low OR High momentum "};
  Configurable<bool> cfgSelAND{"cfgSelAND", false, "Low AND High momentum"};
  Configurable<bool> cfgSelLow{"cfgSelLow", true, "PID selection cut for Low momentum"};
  Configurable<bool> cfgSelHigh{"cfgSelHigh", true, "PID selection cut for High momentum"};
  ConfigurableAxis multTPCBins{"multTPCBins", {150, 0, 150}, "TPC Multiplicity bins"};
  ConfigurableAxis multFT0MBins{"multFT0MBins", {1000, 0, 5000}, "Forward Multiplicity bins"};
  ConfigurableAxis multFT0MMCBins{"multFT0MMCBins", {250, 0, 250}, "Forward Multiplicity bins"};
  ConfigurableAxis dcaXYBins{"dcaXYBins", {100, -0.15, 0.15}, "dcaXY bins"};
  ConfigurableAxis dcaZBins{"dcaZBins", {100, -1.2, 1.2}, "dcaZ bins"};
  ConfigurableAxis qNBins{"qNBins", {1000, 0., 100.}, "nth moments bins"};
  ConfigurableAxis tpNBins{"tpNBins", {300, 0., 3000.}, ""};
  ConfigurableAxis tpDBins{"tpDBins", {100, 0., 2000.}, ""};
  Configurable<std::vector<double>> ptBins{"ptBins", {0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00}, "p_{T} bins"};
  Configurable<std::vector<double>> etaBins{"etaBins", {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8}, "#eta bins"};
  Configurable<std::vector<double>> rapBins{"rapBins", {-0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6}, "#rap bins"};
  Configurable<std::vector<float>> effValuesCh{"effValuesCh", {0, 0.429014, 0.487349, 0.491862, 0.487173, 0.493464, 0.502531, 0.510066, 0.517214, 0.524902, 0.529725, 0.537065, 0.542265, 0.546103, 0.549713, 0.555139, 0.55158, 0.562156, 0.563038, 0.568055, 0.570847, 0.580461, 0.580406, 0.585776, 0.587068, 0.598144, 0.590378, 0.609363, 0.607307, 0.604931, 0.6011, 0.593467, 0.61525, 0.61393, 0.61495, 0.610359, 0.622616}, "effeciency values for Charged Particles"};
  Configurable<std::vector<float>> effPtValuesPi{"effPtValuesPi", {0, 0.410663, 0.480289, 0.494895, 0.487076, 0.489786, 0.49886, 0.493927, 0.39043, 0.243861, 0.238888, 0.229684, 0.232042, 0.236374, 0.240662, 0.243322, 0.244936, 0.247454, 0.250458, 0.251617, 0.255598, 0.258227, 0.262528, 0.266772, 0.272183, 0.279049, 0.279705, 0.283223, 0.285635, 0.287154, 0.288375, 0.291491, 0.294697, 0.295954, 0.298417, 0.304913, 0.31268}, "effeciency values for Pions"};
  Configurable<std::vector<float>> effPtRapValuesPi{"effPtRapValuesPi", {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5515, 0.5337, 0.5204, 0.5118, 0.5072, 0.5195, 0.4326, 0.2341, 0.2669, 0.2604, 0.2484, 0.2547, 0.2539, 0.2608, 0.2621, 0.2653, 0.2653, 0.2783, 0.2765, 0.2858, 0.2932, 0.2862, 0.3041, 0.3038, 0.3081, 0.3024, 0.3153, 0.2991, 0.3121, 0.3094, 0.3203, 0.3320, 0.3225, 0.3355, 0.3491, 0.3485, 0.0000, 0.5483, 0.5306, 0.5231, 0.5126, 0.5090, 0.5168, 0.5040, 0.2427, 0.2606, 0.2519, 0.2405, 0.2453, 0.2465, 0.2479, 0.2513, 0.2540, 0.2560, 0.2579, 0.2625, 0.2597, 0.2627, 0.2795, 0.2765, 0.2871, 0.2798, 0.2882, 0.2981, 0.2993, 0.2950, 0.2918, 0.2912, 0.3173, 0.3076, 0.3088, 0.3235, 0.3261, 0.0000, 0.5353, 0.5284, 0.5192, 0.5120, 0.5136, 0.5209, 0.5247, 0.2871, 0.2642, 0.2601, 0.2534, 0.2576, 0.2600, 0.2653, 0.2696, 0.2700, 0.2753, 0.2899, 0.2821, 0.2855, 0.2894, 0.2952, 0.3073, 0.3100, 0.3200, 0.3173, 0.3219, 0.3235, 0.3241, 0.3328, 0.3322, 0.3361, 0.3375, 0.3344, 0.3449, 0.3659, 0.0000, 0.4807, 0.5206, 0.5190, 0.5131, 0.5132, 0.5204, 0.5263, 0.3672, 0.2951, 0.2831, 0.2745, 0.2762, 0.2819, 0.2864, 0.2893, 0.2866, 0.2926, 0.2957, 0.2996, 0.3020, 0.3081, 0.3163, 0.3141, 0.3305, 0.3314, 0.3359, 0.3342, 0.3425, 0.3287, 0.3405, 0.3284, 0.3474, 0.3738, 0.3428, 0.3608, 0.3580, 0.0000, 0.4313, 0.5049, 0.5164, 0.5139, 0.5121, 0.5180, 0.5238, 0.4100, 0.2872, 0.2788, 0.2665, 0.2682, 0.2722, 0.2775, 0.2754, 0.2784, 0.2823, 0.2814, 0.2872, 0.2874, 0.2828, 0.2891, 0.2997, 0.2997, 0.3147, 0.3147, 0.3067, 0.3095, 0.3204, 0.3334, 0.3118, 0.3263, 0.3210, 0.3130, 0.3477, 0.3453, 0.0000, 0.4125, 0.4795, 0.5132, 0.5083, 0.5124, 0.5182, 0.5244, 0.4443, 0.2584, 0.2502, 0.2405, 0.2397, 0.2446, 0.2512, 0.2589, 0.2565, 0.2583, 0.2605, 0.2628, 0.2681, 0.2657, 0.2772, 0.2782, 0.2802, 0.2908, 0.2900, 0.2886, 0.2959, 0.3061, 0.2919, 0.3010, 0.3056, 0.3061, 0.3054, 0.2973, 0.3227, 0.0000, 0.4045, 0.4708, 0.5112, 0.5070, 0.5084, 0.5167, 0.5201, 0.4732, 0.2244, 0.2256, 0.2162, 0.2244, 0.2247, 0.2376, 0.2323, 0.2374, 0.2460, 0.2430, 0.2437, 0.2511, 0.2524, 0.2519, 0.2660, 0.2606, 0.2837, 0.2833, 0.2877, 0.2872, 0.2945, 0.2936, 0.3049, 0.3080, 0.2987, 0.3073, 0.3041, 0.3108, 0.0000, 0.4038, 0.4932, 0.5118, 0.4980, 0.4996, 0.5104, 0.5156, 0.4926, 0.2091, 0.2117, 0.1980, 0.2032, 0.2099, 0.2100, 0.2183, 0.2140, 0.2196, 0.2166, 0.2175, 0.2266, 0.2251, 0.2312, 0.2320, 0.2411, 0.2451, 0.2397, 0.2608, 0.2515, 0.2619, 0.2578, 0.2603, 0.2557, 0.2651, 0.2695, 0.2740, 0.2689, 0.0000, 0.3453, 0.4426, 0.4595, 0.4485, 0.4580, 0.4717, 0.4742, 0.4735, 0.1814, 0.1763, 0.1679, 0.1683, 0.1728, 0.1761, 0.1760, 0.1788, 0.1781, 0.1797, 0.1828, 0.1802, 0.1869, 0.1880, 0.1859, 0.1923, 0.1920, 0.1992, 0.2081, 0.2101, 0.2082, 0.2137, 0.2171, 0.2091, 0.2110, 0.2063, 0.2134, 0.2060, 0.0000, 0.2718, 0.3691, 0.4004, 0.3955, 0.4086, 0.4176, 0.4272, 0.4323, 0.1217, 0.1172, 0.1120, 0.1133, 0.1177, 0.1198, 0.1249, 0.1227, 0.1244, 0.1278, 0.1202, 0.1236, 0.1245, 0.1315, 0.1290, 0.1343, 0.1417, 0.1452, 0.1454, 0.1354, 0.1481, 0.1474, 0.1445, 0.1551, 0.1389, 0.1463, 0.1488, 0.1391, 0.0000, 0.2701, 0.3734, 0.4018, 0.4004, 0.4118, 0.4228, 0.4314, 0.4347, 0.1340, 0.1342, 0.1294, 0.1252, 0.1288, 0.1308, 0.1342, 0.1330, 0.1393, 0.1404, 0.1359, 0.1420, 0.1377, 0.1478, 0.1472, 0.1513, 0.1506, 0.1593, 0.1519, 0.1551, 0.1589, 0.1640, 0.1595, 0.1608, 0.1605, 0.1592, 0.1734, 0.1677, 0.0000, 0.3251, 0.4435, 0.4686, 0.4533, 0.4644, 0.4754, 0.4806, 0.4789, 0.1828, 0.1822, 0.1735, 0.1722, 0.1753, 0.1780, 0.1814, 0.1853, 0.1798, 0.1908, 0.1844, 0.1901, 0.1910, 0.1975, 0.1890, 0.2014, 0.2058, 0.2071, 0.1975, 0.2117, 0.2097, 0.2101, 0.2176, 0.2159, 0.2182, 0.2200, 0.2242, 0.2088, 0.0000, 0.3560, 0.4871, 0.5035, 0.4921, 0.4955, 0.5032, 0.5127, 0.4853, 0.1938, 0.1908, 0.1845, 0.1879, 0.1945, 0.1946, 0.1996, 0.2004, 0.2016, 0.2026, 0.2010, 0.2053, 0.2072, 0.2175, 0.2206, 0.2153, 0.2357, 0.2307, 0.2407, 0.2388, 0.2332, 0.2447, 0.2377, 0.2434, 0.2373, 0.2499, 0.2482, 0.2567, 0.0000, 0.3632, 0.4656, 0.5044, 0.4982, 0.5003, 0.5099, 0.5196, 0.4650, 0.2335, 0.2316, 0.2204, 0.2243, 0.2285, 0.2333, 0.2384, 0.2369, 0.2404, 0.2394, 0.2435, 0.2452, 0.2501, 0.2549, 0.2584, 0.2574, 0.2716, 0.2784, 0.2863, 0.2815, 0.2886, 0.2961, 0.2897, 0.2822, 0.2868, 0.3006, 0.2948, 0.3187, 0.0000, 0.3911, 0.4746, 0.5049, 0.4997, 0.5024, 0.5111, 0.5190, 0.4462, 0.2897, 0.2828, 0.2753, 0.2772, 0.2820, 0.2854, 0.2963, 0.2941, 0.2953, 0.2974, 0.2978, 0.3012, 0.3070, 0.3044, 0.3207, 0.3238, 0.3322, 0.3292, 0.3315, 0.3457, 0.3413, 0.3252, 0.3420, 0.3471, 0.3459, 0.3501, 0.3417, 0.3833, 0.0000, 0.4216, 0.4934, 0.5105, 0.5060, 0.5009, 0.5082, 0.5107, 0.4247, 0.3192, 0.3097, 0.2968, 0.3011, 0.3039, 0.3098, 0.3115, 0.3138, 0.3145, 0.3173, 0.3238, 0.3320, 0.3319, 0.3304, 0.3348, 0.3478, 0.3522, 0.3428, 0.3499, 0.3588, 0.3683, 0.3532, 0.3796, 0.3497, 0.3672, 0.3789, 0.3828, 0.3958, 0.0000, 0.4697, 0.5062, 0.5111, 0.5027, 0.5016, 0.5113, 0.5137, 0.3797, 0.3255, 0.3218, 0.3066, 0.3100, 0.3154, 0.3256, 0.3241, 0.3282, 0.3313, 0.3314, 0.3430, 0.3401, 0.3518, 0.3527, 0.3581, 0.3683, 0.3702, 0.3795, 0.3675, 0.3780, 0.3752, 0.3810, 0.3896, 0.3945, 0.4056, 0.3943, 0.4078, 0.4469, 0.0000, 0.5251, 0.5210, 0.5129, 0.5006, 0.5013, 0.5083, 0.5143, 0.3072, 0.2999, 0.2948, 0.2847, 0.2867, 0.2997, 0.2996, 0.3032, 0.3062, 0.3110, 0.3171, 0.3177, 0.3219, 0.3322, 0.3327, 0.3382, 0.3499, 0.3549, 0.3451, 0.3620, 0.3693, 0.3554, 0.3625, 0.3679, 0.3682, 0.3714, 0.4006, 0.3892, 0.4137, 0.0000, 0.5409, 0.5208, 0.5121, 0.5037, 0.4999, 0.5079, 0.4914, 0.2571, 0.2818, 0.2732, 0.2670, 0.2650, 0.2723, 0.2755, 0.2748, 0.2813, 0.2831, 0.2864, 0.2907, 0.2959, 0.2949, 0.2998, 0.3017, 0.3166, 0.3149, 0.3134, 0.3309, 0.3309, 0.3257, 0.3308, 0.3304, 0.3290, 0.3346, 0.3499, 0.3636, 0.3523, 0.0000, 0.5406, 0.5250, 0.5103, 0.5001, 0.4989, 0.5082, 0.4259, 0.2494, 0.2792, 0.2714, 0.2662, 0.2694, 0.2742, 0.2797, 0.2787, 0.2884, 0.2851, 0.2887, 0.2947, 0.3003, 0.3046, 0.3032, 0.3169, 0.3152, 0.3199, 0.3331, 0.3213, 0.3302, 0.3299, 0.3293, 0.3335, 0.3508, 0.3530, 0.3372, 0.3434, 0.3631, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, "pT rap effeciency values for Pions"};
  Configurable<std::vector<float>> effPtValuesKa{"effPtValuesKa", {0, 0, 0, 0.328845, 0.379771, 0.390088, 0.403074, 0.35504, 0.256438, 0.131726, 0.13796, 0.140295, 0.147229, 0.156968, 0.162245, 0.171312, 0.175851, 0.185823, 0.188763, 0.193965, 0.192999, 0.191121, 0.195547, 0.210082, 0.217502, 0.232456, 0.245035, 0.254051, 0.268206, 0.274664, 0.290428, 0.294979, 0.304817, 0.324206, 0.342578, 0.36466, 0.394134}, "effeciency values for Kaons"};
  Configurable<std::vector<float>> effPtRapValuesKa{"effPtRapValuesKa", {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1710, 0.4333, 0.4561, 0.4428, 0.1393, 0.1635, 0.1836, 0.1883, 0.1719, 0.1795, 0.1780, 0.1907, 0.2059, 0.2126, 0.2202, 0.2160, 0.2130, 0.2434, 0.2642, 0.2601, 0.2988, 0.2967, 0.3087, 0.3210, 0.3659, 0.3808, 0.3560, 0.4011, 0.4004, 0.4507, 0.5358, 0.5475, 0.4919, 0.5736, 0.0000, 0.0000, 0.0000, 0.4064, 0.4340, 0.4523, 0.4884, 0.2246, 0.1275, 0.1524, 0.1638, 0.1678, 0.1899, 0.2055, 0.2134, 0.2186, 0.2245, 0.2655, 0.2590, 0.2420, 0.2367, 0.2698, 0.2533, 0.2606, 0.2782, 0.3164, 0.3342, 0.3403, 0.3238, 0.3297, 0.3572, 0.4031, 0.3991, 0.4431, 0.4791, 0.5160, 0.5619, 0.0000, 0.0000, 0.0000, 0.4011, 0.4186, 0.4499, 0.4929, 0.3641, 0.1376, 0.1716, 0.1813, 0.1748, 0.1816, 0.1987, 0.2031, 0.2157, 0.2189, 0.2383, 0.2336, 0.2463, 0.2388, 0.2608, 0.2532, 0.2906, 0.2875, 0.2858, 0.3218, 0.3383, 0.3852, 0.4082, 0.4257, 0.3886, 0.4276, 0.4136, 0.5055, 0.5401, 0.5559, 0.0000, 0.0000, 0.0000, 0.4018, 0.4292, 0.4693, 0.4800, 0.4941, 0.1452, 0.1656, 0.1789, 0.1754, 0.1932, 0.1992, 0.2124, 0.2334, 0.2529, 0.2683, 0.2569, 0.2610, 0.2560, 0.2817, 0.2783, 0.2701, 0.2927, 0.3348, 0.3479, 0.3598, 0.3752, 0.4069, 0.3838, 0.4251, 0.4384, 0.4489, 0.4720, 0.5712, 0.5559, 0.0000, 0.0000, 0.0000, 0.4044, 0.4158, 0.4518, 0.4736, 0.5041, 0.2438, 0.1757, 0.1872, 0.1927, 0.2077, 0.2163, 0.2303, 0.2333, 0.2466, 0.2500, 0.2546, 0.2598, 0.2597, 0.2549, 0.2765, 0.2824, 0.2953, 0.3066, 0.3252, 0.3230, 0.3371, 0.3696, 0.3785, 0.4025, 0.3631, 0.3894, 0.4227, 0.4968, 0.5269, 0.0000, 0.0000, 0.0000, 0.3772, 0.4052, 0.4427, 0.4735, 0.4973, 0.3424, 0.1813, 0.1817, 0.1755, 0.1913, 0.1964, 0.1931, 0.2108, 0.2139, 0.2267, 0.2449, 0.2446, 0.2306, 0.2203, 0.2281, 0.2395, 0.2564, 0.2729, 0.2793, 0.2971, 0.3204, 0.3486, 0.3403, 0.3598, 0.3336, 0.4031, 0.3827, 0.4153, 0.4825, 0.0000, 0.0000, 0.0000, 0.3692, 0.4236, 0.4238, 0.4535, 0.4874, 0.4023, 0.1581, 0.1600, 0.1752, 0.1779, 0.1921, 0.2002, 0.2115, 0.1938, 0.2353, 0.2190, 0.2317, 0.2499, 0.2108, 0.2228, 0.2372, 0.2598, 0.2897, 0.2968, 0.3133, 0.3346, 0.3145, 0.3446, 0.3797, 0.3506, 0.3971, 0.3812, 0.4328, 0.4627, 0.0000, 0.0000, 0.0000, 0.3574, 0.4078, 0.4161, 0.4555, 0.4859, 0.4549, 0.1234, 0.1394, 0.1447, 0.1518, 0.1617, 0.1720, 0.1870, 0.1899, 0.1990, 0.2060, 0.2290, 0.2154, 0.1940, 0.2001, 0.2104, 0.2297, 0.2476, 0.2475, 0.3002, 0.2997, 0.2965, 0.3335, 0.2995, 0.3265, 0.3485, 0.3586, 0.3593, 0.3895, 0.0000, 0.0000, 0.0000, 0.3515, 0.3866, 0.3904, 0.4175, 0.4436, 0.4482, 0.1154, 0.1325, 0.1284, 0.1375, 0.1461, 0.1515, 0.1668, 0.1598, 0.1779, 0.1718, 0.1851, 0.1923, 0.1751, 0.1677, 0.1788, 0.1905, 0.2016, 0.2277, 0.2293, 0.2357, 0.2365, 0.2547, 0.2547, 0.2904, 0.3060, 0.2770, 0.3024, 0.3044, 0.0000, 0.0000, 0.0000, 0.2892, 0.3182, 0.3325, 0.3723, 0.3948, 0.4164, 0.0746, 0.0849, 0.0802, 0.0831, 0.0977, 0.0977, 0.1076, 0.1055, 0.1078, 0.1292, 0.1338, 0.1203, 0.1220, 0.1222, 0.1269, 0.1520, 0.1471, 0.1481, 0.1551, 0.1632, 0.1747, 0.1881, 0.1929, 0.1952, 0.1882, 0.2317, 0.2041, 0.2184, 0.0000, 0.0000, 0.0000, 0.1412, 0.1523, 0.1380, 0.1249, 0.1134, 0.1096, 0.0206, 0.0215, 0.0217, 0.0234, 0.0245, 0.0235, 0.0264, 0.0260, 0.0259, 0.0306, 0.0301, 0.0329, 0.0278, 0.0293, 0.0291, 0.0285, 0.0338, 0.0367, 0.0327, 0.0347, 0.0379, 0.0398, 0.0394, 0.0415, 0.0432, 0.0459, 0.0455, 0.0448, 0.0000, 0.0000, 0.0000, 0.3380, 0.3800, 0.3929, 0.4257, 0.4659, 0.4770, 0.1180, 0.1269, 0.1332, 0.1331, 0.1554, 0.1509, 0.1550, 0.1511, 0.1614, 0.1707, 0.1756, 0.1879, 0.1745, 0.1814, 0.1935, 0.2082, 0.2145, 0.2267, 0.2253, 0.2399, 0.2584, 0.2321, 0.2688, 0.2835, 0.2905, 0.2703, 0.2646, 0.3305, 0.0000, 0.0000, 0.0000, 0.3511, 0.3996, 0.4104, 0.4351, 0.4728, 0.4420, 0.1167, 0.1296, 0.1344, 0.1395, 0.1582, 0.1586, 0.1630, 0.1761, 0.1778, 0.1758, 0.2058, 0.2141, 0.1884, 0.1912, 0.2007, 0.1999, 0.2372, 0.2419, 0.2498, 0.2777, 0.2696, 0.2919, 0.2636, 0.2878, 0.2997, 0.3078, 0.3120, 0.3634, 0.0000, 0.0000, 0.0000, 0.3515, 0.4008, 0.4128, 0.4487, 0.4714, 0.3919, 0.1661, 0.1741, 0.1862, 0.1905, 0.2126, 0.2116, 0.2105, 0.2297, 0.2320, 0.2338, 0.2458, 0.2438, 0.2333, 0.2475, 0.2702, 0.2490, 0.2799, 0.3052, 0.3312, 0.3008, 0.3384, 0.3727, 0.3729, 0.3636, 0.3633, 0.4090, 0.4434, 0.4768, 0.0000, 0.0000, 0.0000, 0.3819, 0.4277, 0.4347, 0.4577, 0.4789, 0.3477, 0.1967, 0.2060, 0.2049, 0.2209, 0.2227, 0.2315, 0.2301, 0.2438, 0.2524, 0.2502, 0.2623, 0.2576, 0.2564, 0.2699, 0.2855, 0.3143, 0.2979, 0.3263, 0.3232, 0.3717, 0.3753, 0.4165, 0.4346, 0.4080, 0.4615, 0.4266, 0.4392, 0.5890, 0.0000, 0.0000, 0.0000, 0.3925, 0.4204, 0.4355, 0.4655, 0.5137, 0.2522, 0.1977, 0.2036, 0.2264, 0.2299, 0.2466, 0.2628, 0.2639, 0.2634, 0.2815, 0.3062, 0.3054, 0.2843, 0.2735, 0.2903, 0.3024, 0.3271, 0.3474, 0.3701, 0.3710, 0.3882, 0.4187, 0.4192, 0.4064, 0.4700, 0.4275, 0.4634, 0.5384, 0.5804, 0.0000, 0.0000, 0.0000, 0.4013, 0.4256, 0.4442, 0.4596, 0.4860, 0.1637, 0.1880, 0.1854, 0.2001, 0.2260, 0.2395, 0.2359, 0.2605, 0.2761, 0.2912, 0.3005, 0.3060, 0.2783, 0.2967, 0.3172, 0.3477, 0.3610, 0.3755, 0.3853, 0.4036, 0.4232, 0.4345, 0.4714, 0.4609, 0.4807, 0.4760, 0.5307, 0.5920, 0.6667, 0.0000, 0.0000, 0.0000, 0.4080, 0.4225, 0.4587, 0.4873, 0.3585, 0.1490, 0.1847, 0.1978, 0.2024, 0.1907, 0.2021, 0.2301, 0.2312, 0.2396, 0.2437, 0.2615, 0.2626, 0.2670, 0.2763, 0.2692, 0.3377, 0.2918, 0.3532, 0.3584, 0.4009, 0.4231, 0.3914, 0.4495, 0.4513, 0.4265, 0.5075, 0.5419, 0.5530, 0.5991, 0.0000, 0.0000, 0.0000, 0.4028, 0.4344, 0.4443, 0.4706, 0.2189, 0.1332, 0.1607, 0.1688, 0.1826, 0.1941, 0.2104, 0.2330, 0.2464, 0.2556, 0.2657, 0.2593, 0.2694, 0.2570, 0.2585, 0.2601, 0.3015, 0.3106, 0.3096, 0.3325, 0.3396, 0.3926, 0.3547, 0.3922, 0.4017, 0.4041, 0.4452, 0.5183, 0.5656, 0.6020, 0.0000, 0.0000, 0.0000, 0.1668, 0.4392, 0.4493, 0.4341, 0.1520, 0.1720, 0.2046, 0.1954, 0.1838, 0.1763, 0.1911, 0.1864, 0.2099, 0.2275, 0.2336, 0.2244, 0.2312, 0.2482, 0.2544, 0.2762, 0.2956, 0.3206, 0.3364, 0.3618, 0.3787, 0.3705, 0.4065, 0.4193, 0.4130, 0.4612, 0.5291, 0.5680, 0.5802, 0.5690, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, "pT rap effeciency values for Kaons"};
  Configurable<std::vector<float>> effPtValuesPr{"effPtValuesPr", {0, 0, 0, 0, 0, 0, 0, 0.413799, 0.443597, 0.478144, 0.505512, 0.514127, 0.523279, 0.506121, 0.481129, 0.436858, 0.365426, 0.244158, 0.246106, 0.249133, 0.251019, 0.256516, 0.263027, 0.258241, 0.260814, 0.270519, 0.272534, 0.271853, 0.274523, 0.279029, 0.279756, 0.285479, 0.292531, 0.292348, 0.294704, 0.295867, 0.289818}, "effeciency values for Kaons"};
  Configurable<std::vector<float>> effPtRapValuesPr{"effPtRapValuesPr", {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0339, 0.3198, 0.6225, 0.6961, 0.6865, 0.3870, 0.3603, 0.3747, 0.3765, 0.3812, 0.3860, 0.3793, 0.3897, 0.3713, 0.3408, 0.3204, 0.3505, 0.3313, 0.3189, 0.3179, 0.3714, 0.3450, 0.3659, 0.3303, 0.3727, 0.3353, 0.4124, 0.3755, 0.3680, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.2453, 0.5687, 0.6297, 0.6445, 0.6864, 0.7439, 0.6205, 0.3602, 0.3300, 0.3139, 0.3066, 0.2902, 0.2910, 0.3169, 0.3290, 0.3596, 0.3345, 0.3379, 0.3808, 0.3658, 0.3311, 0.3829, 0.4252, 0.3592, 0.3778, 0.4173, 0.4083, 0.3905, 0.4169, 0.3470, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5466, 0.6134, 0.6293, 0.6451, 0.6817, 0.7348, 0.6971, 0.4628, 0.3037, 0.3039, 0.3051, 0.3492, 0.3391, 0.3190, 0.3521, 0.3582, 0.3402, 0.3278, 0.3477, 0.3936, 0.3571, 0.3591, 0.3739, 0.3902, 0.3732, 0.4123, 0.3741, 0.3137, 0.3740, 0.3497, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5707, 0.6136, 0.6606, 0.6479, 0.6613, 0.6752, 0.6961, 0.7150, 0.3865, 0.3336, 0.3168, 0.3175, 0.3289, 0.3310, 0.3365, 0.3517, 0.3521, 0.3274, 0.3496, 0.3850, 0.3615, 0.3780, 0.3802, 0.4004, 0.4011, 0.4142, 0.4057, 0.4173, 0.4073, 0.4222, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5873, 0.6086, 0.6380, 0.6571, 0.6536, 0.7018, 0.7207, 0.7215, 0.5860, 0.3364, 0.3081, 0.3494, 0.3398, 0.3491, 0.3680, 0.3567, 0.4204, 0.3898, 0.3957, 0.4100, 0.3749, 0.4093, 0.3917, 0.4212, 0.3910, 0.4062, 0.4019, 0.4378, 0.3997, 0.3769, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5695, 0.5995, 0.6392, 0.6646, 0.6341, 0.6736, 0.7025, 0.7104, 0.7144, 0.4082, 0.3392, 0.3319, 0.3259, 0.3929, 0.3720, 0.3617, 0.3716, 0.3394, 0.3454, 0.3520, 0.3409, 0.3872, 0.3775, 0.3134, 0.3928, 0.3635, 0.3540, 0.4066, 0.3520, 0.3580, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5524, 0.6210, 0.6435, 0.6720, 0.6714, 0.6821, 0.7348, 0.7403, 0.7531, 0.5534, 0.3513, 0.3381, 0.3319, 0.3384, 0.3426, 0.3254, 0.3360, 0.3503, 0.3264, 0.3306, 0.3642, 0.3545, 0.3478, 0.3527, 0.3625, 0.3678, 0.4036, 0.3569, 0.3321, 0.3761, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5472, 0.5793, 0.6230, 0.6680, 0.6685, 0.6886, 0.7368, 0.7534, 0.7767, 0.6561, 0.2933, 0.2798, 0.2943, 0.2940, 0.3034, 0.3104, 0.3133, 0.3168, 0.3209, 0.3018, 0.3135, 0.3137, 0.3315, 0.3580, 0.3384, 0.3454, 0.3307, 0.3182, 0.3801, 0.3145, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5479, 0.5660, 0.5973, 0.6225, 0.6437, 0.6655, 0.6920, 0.7241, 0.7637, 0.7145, 0.2647, 0.2727, 0.2858, 0.2508, 0.2523, 0.2647, 0.2592, 0.2438, 0.3017, 0.2916, 0.2582, 0.2705, 0.2876, 0.2960, 0.2958, 0.2982, 0.3014, 0.2978, 0.3077, 0.3234, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.4692, 0.5078, 0.5095, 0.5695, 0.6257, 0.5869, 0.6149, 0.6597, 0.7159, 0.6730, 0.1877, 0.1686, 0.1809, 0.1793, 0.1636, 0.1880, 0.1822, 0.1980, 0.2000, 0.1931, 0.1960, 0.1895, 0.1921, 0.2056, 0.2194, 0.2048, 0.1942, 0.1944, 0.2267, 0.1756, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1458, 0.1090, 0.1011, 0.0956, 0.0961, 0.0920, 0.0934, 0.1006, 0.1015, 0.1028, 0.0299, 0.0291, 0.0298, 0.0290, 0.0285, 0.0326, 0.0318, 0.0316, 0.0335, 0.0307, 0.0370, 0.0346, 0.0373, 0.0330, 0.0404, 0.0409, 0.0379, 0.0393, 0.0384, 0.0347, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5459, 0.5794, 0.5870, 0.6118, 0.6806, 0.6603, 0.6685, 0.6983, 0.7300, 0.7229, 0.2701, 0.2476, 0.2332, 0.2516, 0.2832, 0.2878, 0.2401, 0.2429, 0.2737, 0.2909, 0.2821, 0.2739, 0.2935, 0.2763, 0.2905, 0.2901, 0.2624, 0.2743, 0.2765, 0.2768, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5492, 0.5688, 0.6025, 0.6438, 0.6553, 0.6753, 0.6955, 0.7401, 0.7443, 0.6152, 0.2739, 0.2595, 0.2627, 0.2448, 0.2755, 0.2661, 0.2868, 0.2728, 0.2766, 0.2822, 0.2971, 0.2699, 0.2906, 0.2847, 0.2991, 0.2922, 0.2971, 0.3115, 0.3190, 0.2911, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5386, 0.5562, 0.6265, 0.6461, 0.6657, 0.6662, 0.6788, 0.7231, 0.7452, 0.5331, 0.3690, 0.3522, 0.3639, 0.3608, 0.3767, 0.3659, 0.3698, 0.3622, 0.3835, 0.3945, 0.3727, 0.3756, 0.3696, 0.3793, 0.3974, 0.3849, 0.4076, 0.4034, 0.3830, 0.3997, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5549, 0.5906, 0.6211, 0.6459, 0.6337, 0.6527, 0.7139, 0.7135, 0.6978, 0.4328, 0.3874, 0.4117, 0.3946, 0.4086, 0.3969, 0.4171, 0.4177, 0.4252, 0.4013, 0.3983, 0.3879, 0.3580, 0.4246, 0.3776, 0.4004, 0.4148, 0.4322, 0.4258, 0.4097, 0.4332, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5714, 0.6016, 0.6489, 0.6461, 0.6212, 0.6422, 0.6781, 0.6867, 0.5264, 0.3990, 0.3788, 0.3841, 0.3841, 0.4045, 0.4228, 0.4481, 0.4041, 0.4295, 0.4772, 0.4250, 0.4547, 0.4662, 0.4123, 0.4281, 0.3978, 0.4676, 0.4783, 0.4357, 0.4261, 0.5097, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5979, 0.5975, 0.6258, 0.6617, 0.6368, 0.6823, 0.7199, 0.7093, 0.3680, 0.3414, 0.3514, 0.3390, 0.3770, 0.3479, 0.3471, 0.3840, 0.3779, 0.3988, 0.4329, 0.4114, 0.4216, 0.3994, 0.4046, 0.4322, 0.4386, 0.4179, 0.4636, 0.4501, 0.4376, 0.4547, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5787, 0.6295, 0.6553, 0.6404, 0.6299, 0.6723, 0.6928, 0.4525, 0.3031, 0.3146, 0.3160, 0.3328, 0.3610, 0.3694, 0.3414, 0.3894, 0.3549, 0.3514, 0.3854, 0.3937, 0.4015, 0.3991, 0.3783, 0.3765, 0.4228, 0.4194, 0.3932, 0.4179, 0.4323, 0.3799, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.2257, 0.5829, 0.6392, 0.6577, 0.6744, 0.7086, 0.6263, 0.3869, 0.3645, 0.3169, 0.2885, 0.2944, 0.2982, 0.2959, 0.3214, 0.3271, 0.3333, 0.3841, 0.3707, 0.4022, 0.3970, 0.3747, 0.3930, 0.3812, 0.3923, 0.3872, 0.4083, 0.4625, 0.4091, 0.4308, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0343, 0.3174, 0.6091, 0.6870, 0.6767, 0.4265, 0.4035, 0.3921, 0.4054, 0.3786, 0.3852, 0.3952, 0.3551, 0.3799, 0.3574, 0.3478, 0.3342, 0.3269, 0.3239, 0.3679, 0.3209, 0.3090, 0.3767, 0.3253, 0.3707, 0.3915, 0.4005, 0.4359, 0.4303, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, "pT rap effeciency values for Protons"};

  using MyAllTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                                aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                                aod::pidTOFbeta, aod::pidTOFmass>;
  using MyRun3Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::CentFT0Cs, aod::CentFT0Ms>;
  using MyRun3MCCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::CentFT0Cs, aod::CentFT0Ms, aod::McCollisionLabels>;
  using MyMCTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                               aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                               aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                               aod::pidTOFbeta, aod::pidTOFmass, aod::McTrackLabels>;

  Service<o2::framework::O2DatabasePDG> pdg;

  HistogramRegistry hist{"hist", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    const AxisSpec axisEvents{10, 0, 10, "Counts"};
    const AxisSpec axisEta{etaBins, "#eta"};
    const AxisSpec axisPhi{nPhiBins, 0., +7., "#phi (rad)"};
    const AxisSpec axisY{rapBins, "y"};
    const AxisSpec axisPt{ptBins, "p_{T} (GeV/c)"};
    const AxisSpec axisP{nPBins, 0., 3., "p (GeV/c)"};
    const AxisSpec axisInnerParam{nPBins, 0., 3., "p_{InnerParam } (GeV/c)"};
    const AxisSpec axisPart{nPartBins, 0., 18., " "};
    const AxisSpec axisQn{qNBins, ""};
    const AxisSpec axisTpN{tpNBins, "(Q_{1}^{2} - Q_{2})"};
    const AxisSpec axisTpD{tpDBins, "N_{pairs}"};
    const AxisSpec axisDeno{100, 1., 2.0, "#frac{1}{#sqrt{1 - #frac{1}{N}}}"};
    const AxisSpec axisMeanPt{100, 0., 3., "M(p_{T}) (GeV/c)"};
    const AxisSpec axisMult{100, 0, 100, "N_{ch}"};
    const AxisSpec axisMultTPC{multTPCBins, "N_{TPC} "};
    const AxisSpec axisMultFT0M{multFT0MBins, "N_{FT0M}"};
    const AxisSpec axisMultFT0MMC{multFT0MMCBins, "N_{FT0M}"};
    const AxisSpec axisCentFT0C{nCentBins, 0, 101, "FT0C (%)"};
    const AxisSpec axisVtxZ{80, -20., 20., "V_{Z} (cm)"};
    const AxisSpec axisDCAz{dcaZBins, "DCA_{Z} (cm)"};
    const AxisSpec axisDCAxy{dcaXYBins, "DCA_{XY} (cm)"};
    const AxisSpec axisTPCNsigma{500, -5., 5., "n #sigma_{TPC}"};
    const AxisSpec axisTOFNsigma{500, -5., 5., "n #sigma_{TOF}"};
    const AxisSpec axisTPCSignal{100, 20., 500., "#frac{dE}{dx}"};
    const AxisSpec axisTOFSignal{200, 0.2, 1.2, "TOF #beta"};
    const AxisSpec axisChi2{40, 0., 40., "Chi2"};
    const AxisSpec axisCrossedTPC{300, 0, 300, "Crossed TPC"};
    const AxisSpec axisM2{100, 0., 1.4, "#it{m}^{2} (GeV/#it{c}^{2})^{2}"};
    const AxisSpec axisMass{1000, 0., 0.1, "M_{inv} (GeV/#it{c}^2)"};

    HistogramConfigSpec qNHist({HistType::kTHnSparseD, {axisMultTPC, axisQn, axisMultFT0M}});
    HistogramConfigSpec partHist({HistType::kTHnSparseD, {axisMultTPC, axisPart, axisMultFT0M}});
    HistogramConfigSpec denoHist({HistType::kTHnSparseD, {axisMultTPC, axisDeno, axisMultFT0M}});
    HistogramConfigSpec qNMCHist({HistType::kTHnSparseD, {axisMultTPC, axisQn, axisMultFT0M}});
    HistogramConfigSpec partMCHist({HistType::kTHnSparseD, {axisMultTPC, axisPart, axisMultFT0M}});
    HistogramConfigSpec denoMCHist({HistType::kTHnSparseD, {axisMultTPC, axisDeno, axisMultFT0M}});
    HistogramConfigSpec tofNSigmaHist({HistType::kTH2D, {axisP, axisTOFNsigma}});
    HistogramConfigSpec tofSignalHist({HistType::kTH2D, {axisP, axisTOFSignal}});
    HistogramConfigSpec tpcNSigmaHist({HistType::kTH2D, {axisP, axisTPCNsigma}});
    HistogramConfigSpec tpcSignalHist({HistType::kTH2D, {axisP, axisTPCSignal}});
    HistogramConfigSpec tpcTofHist({HistType::kTH2D, {axisTPCNsigma, axisTOFNsigma}});
    HistogramConfigSpec pvsM2Hist({HistType::kTH2D, {axisM2, axisP}});

    HistogramConfigSpec tofNSigmaHist1({HistType::kTH2D, {axisInnerParam, axisTOFNsigma}});
    HistogramConfigSpec tofSignalHist1({HistType::kTH2D, {axisInnerParam, axisTOFSignal}});
    HistogramConfigSpec tpcNSigmaHist1({HistType::kTH2D, {axisInnerParam, axisTPCNsigma}});
    HistogramConfigSpec tpcSignalHist1({HistType::kTH2D, {axisInnerParam, axisTPCSignal}});
    HistogramConfigSpec tpcTofHist1({HistType::kTH2D, {axisTPCNsigma, axisTOFNsigma}});
    HistogramConfigSpec pvsM2Hist1({HistType::kTH2D, {axisM2, axisInnerParam}});

    // QA Plots:
    hist.add("QA/before/h_Counts", "Counts", kTH1D, {axisEvents});
    hist.add("QA/before/h_VtxZ", "V_{Z}", kTH1D, {axisVtxZ});
    hist.add("QA/before/h_Pt", "p_{T}", kTH1D, {axisPt});
    hist.add("QA/before/h_Eta", "#eta ", kTH1D, {axisEta});
    hist.add("QA/before/h_Phi", "#phi ", kTH1D, {axisPhi});
    hist.add("QA/before/h_DcaZ", "DCA_{Z}", kTH1D, {axisDCAz});
    hist.add("QA/before/h_DcaXY", "DCA_{XY}", kTH1D, {axisDCAxy});
    hist.add("QA/before/h2_DcaZ", "DCA_{Z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/before/h2_DcaXY", "DCA_{XY}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/before/h_NTPC", "N_{TPC}", kTH1D, {axisMultTPC});
    hist.add("QA/before/h_NFT0M", "FT0M Multiplicity", kTH1D, {axisMultFT0M});
    hist.add("QA/before/h_Cent", "FT0C (%)", kTH1D, {axisCentFT0C});
    hist.add("QA/before/h_CentM", "FT0M (%)", kTH1D, {axisCentFT0C});

    hist.add("QA/before/h2_TPCSignal", "TPC Signal", tpcSignalHist);
    hist.add("QA/before/h2_TOFSignal", "TOF Signal", tofSignalHist);
    hist.add("QA/before/h2_pvsm2", "p vs m^{2}", pvsM2Hist);

    hist.add("QA/before/innerParam/h2_TPCSignal", "TPC Signal", tpcSignalHist1);
    hist.add("QA/before/innerParam/h2_TOFSignal", "TOF Signal", tofSignalHist1);
    hist.add("QA/before/innerParam/h2_pvsm2", "p vs m^{2}", pvsM2Hist1);

    hist.addClone("QA/before/", "QA/after/");

    hist.add("QA/after/h_TPCChi2perCluster", "TPC #Chi^{2}/Cluster", kTH1D, {axisChi2});
    hist.add("QA/after/h_ITSChi2perCluster", "ITS #Chi^{2}/Cluster", kTH1D, {axisChi2});
    hist.add("QA/after/h_crossedTPC", "Crossed TPC", kTH1D, {axisCrossedTPC});
    hist.add("QA/after/h_NFT0C", "FT0C Multiplicity", kTH1D, {axisMultFT0M});
    hist.add("QA/after/h_invMass_gamma", "Inv Mass of #gamma", kTH1D, {axisMass});
    hist.add("QA/after/h_counts_evSelCuts", "Event selection cuts", kTH1D, {axisEvents});
    hist.add("QA/after/h_VtxZReco", "Simulated Vertex Z", kTH1D, {axisVtxZ});

    hist.add("QA/after/h2_PvsPinner", "p_{InnerParam} vs p", kTH2D, {{axisP}, {axisInnerParam}});
    hist.add("QA/after/h2_Pt_Eta", "p_{T} vs #eta ", kTH2D, {{axisEta}, {axisPt}});
    hist.add("QA/after/h2_NTPC_Cent", "N_{TPC} vs FT0C(%)", kTH2D, {{axisCentFT0C}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NFT0M", "N_{TPC} vs N_{FT0M}", kTH2D, {{axisMultFT0M}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NFT0C", "N_{TPC} vs N_{FT0C}", kTH2D, {{axisMultFT0M}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NCh", "N_{ch} vs N_{TPC}", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NPi", "N_{Pi} vs N_{TPC}", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NKa", "N_{Ka} vs N_{TPC}", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NPr", "N_{Pr} vs N_{TPC}", kTH2D, {{axisMultTPC}, {axisMultTPC}});

    hist.add("QA/after/p_NTPC_NFT0M", "N_{TPC} vs N_{FT0M} (Profile)", kTProfile, {axisMultFT0M});
    hist.add("QA/after/p_NTPC_NFT0C", "N_{TPC} vs N_{FT0C} (Profile)", kTProfile, {axisMultFT0M});
    hist.add("QA/after/p_NTPC_Cent", "N_{TPC} vs FT0C(%) (Profile)", kTProfile, {axisCentFT0C});

    hist.add("QA/after/h_Pt_weighted", "weighted pT distribution", kTH1D, {axisPt});
    hist.add("QA/after/h2_Pt_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/after/h2_pt_nch", "Truth", kTH2D, {{axisMult}, {axisPt}});
    hist.add("QA/after/h2_pt_nch_prof", "Truth", kTProfile, {axisMult});

    hist.add("QA/Pion/before/h2_TPCNsigma", "n #sigma_{TPC}", tpcNSigmaHist);
    hist.add("QA/Pion/before/h2_TOFNsigma", "n #sigma_{TOF}", tofNSigmaHist);
    hist.add("QA/Pion/before/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", tpcTofHist);

    hist.add("QA/Pion/h_Rap", "y ", kTH1D, {axisY});
    hist.add("QA/Pion/h_RapTruth", "y ", kTH1D, {axisY});
    hist.add("QA/Pion/h_Eta", "Pseudorapidity ", kTH1D, {axisEta});
    hist.add("QA/Pion/h_EtaTruth", "Pseudorapidity (Reco Truth) ", kTH1D, {axisEta});
    hist.add("QA/Pion/h_Phi", "Azimuthal Distribution ", kTH1D, {axisPhi});
    hist.add("QA/Pion/h_DcaZ", "DCA_{z}", kTH1D, {axisDCAz});
    hist.add("QA/Pion/h_DcaXY", "DCA_{xy}", kTH1D, {axisDCAxy});
    hist.add("QA/Pion/h_Pt", "p_{T} ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtPos", "p_{T} (positive) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtNeg", "p_{T} (negative) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtTruth", "p_{T} ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtPosTruth", "p_{T} (positive) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtNegTruth", "p_{T} (negative) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_Pt_weighted", "weighted pT distribution", kTH1D, {axisPt});

    hist.add("QA/Pion/h2_Pt_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtPos_rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtNeg_rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtPosTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtNegTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_Pt_Eta", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPt}});
    hist.add("QA/Pion/h2_DcaZ", "DCA_{z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/Pion/h2_DcaXY", "DCA_{xy}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/Pion/h2_Pt_Rap_weighted", "p_{T} vs #eta weighted", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_Pt_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_PtPos_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_PtNeg_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_PtTruth_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_PtPosTruth_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_PtNegTruth_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_pt_nch", "Truth", kTH2D, {{axisMult}, {axisPt}});
    hist.add("QA/Pion/h2_pt_nch_prof", "Truth", kTProfile, {axisMult});

    hist.add("QA/Pion/h2_TPCNsigma", "n #sigma_{TPC}", tpcNSigmaHist);
    hist.add("QA/Pion/h2_TPCNsigma_El", "n #sigma_{TPC, El}", tpcNSigmaHist);
    hist.add("QA/Pion/h2_TOFNsigma_El", "n #sigma_{TOF, El}", tofNSigmaHist);
    hist.add("QA/Pion/h2_TOFNsigma", "n #sigma_{TOF}", tofNSigmaHist);
    hist.add("QA/Pion/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", tpcTofHist);
    hist.add("QA/Pion/h2_TPCSignal", "TPC Signal ", tpcSignalHist);
    hist.add("QA/Pion/h2_TOFSignal", "TOF Signal", tofSignalHist);
    hist.add("QA/Pion/h2_pvsm2", "p vs m^{2}", pvsM2Hist);

    hist.add("QA/Pion/innerParam/before/h2_TPCNsigma", "n #sigma_{TPC}", tpcNSigmaHist1);
    hist.add("QA/Pion/innerParam/before/h2_TOFNsigma", "n #sigma_{TOF}", tofNSigmaHist1);
    hist.add("QA/Pion/innerParam/before/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", tpcTofHist1);
    hist.add("QA/Pion/innerParam/h2_TPCNsigma", "n #sigma_{TPC}", tpcNSigmaHist1);
    hist.add("QA/Pion/innerParam/h2_TPCNsigma_El", "n #sigma_{TPC, El}", tpcNSigmaHist1);
    hist.add("QA/Pion/innerParam/h2_TOFNsigma_El", "n #sigma_{TOF, El}", tofNSigmaHist1);
    hist.add("QA/Pion/innerParam/h2_TOFNsigma", "n #sigma_{TOF}", tofNSigmaHist1);
    hist.add("QA/Pion/innerParam/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", tpcTofHist1);
    hist.add("QA/Pion/innerParam/h2_TPCSignal", "TPC Signal ", tpcSignalHist1);
    hist.add("QA/Pion/innerParam/h2_TOFSignal", "TOF Signal", tofSignalHist1);
    hist.add("QA/Pion/innerParam/h2_pvsm2", "p vs m^{2}", pvsM2Hist1);

    hist.addClone("QA/Pion/", "QA/Kaon/");
    hist.addClone("QA/Pion/", "QA/Proton/");

    // Analysis Plots:
    hist.add("Analysis/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Analysis/Charged/h_Q1", "Q1", qNHist);
    hist.add("Analysis/Charged/h_Q2", "Q2", qNHist);
    hist.add("Analysis/Charged/h_Q3", "Q3", qNHist);
    hist.add("Analysis/Charged/h_Q4", "Q4", qNHist);
    hist.add("Analysis/Charged/h_mean_pT", " <p_{T}> ", kTH1D, {axisMeanPt});
    hist.add("Analysis/Charged/p_mean_pT_Mult_var", " <p_{T}> ", kTProfile, {axisMultTPC});
    hist.add("Analysis/Charged/p_CheckNCh", " 1/denominator vs N_{TPC} ", kTProfile, {axisMultTPC});
    hist.add("Analysis/Charged/h_CheckNCh", " 1/denominator vs N_{TPC} ", denoHist);
    hist.add("Analysis/Charged/h_Q1_var", "Q1 vs N_{TPC}", qNHist);
    hist.add("Analysis/Charged/h_N_var", "N vs N_{TPC}", kTHnSparseD, {axisMultTPC, axisMult, axisMultFT0M});
    hist.add("Analysis/Charged/h_twopart_nume_Mult_var", "twopart numerator", kTHnSparseD, {axisMultTPC, axisTpN, axisMultFT0M});
    hist.add("Analysis/Charged/h_twopart_deno_Mult_var", "twopart denominator", kTHnSparseD, {axisMultTPC, axisTpD, axisMultFT0M});
    hist.add("Analysis/Charged/h_mean_pT_Mult_var", " <p_{T}> vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_mean_pT_Mult_skew", " <p_{T}> vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_mean_pT_Mult_kurto", " <p_{T}> vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_twopart_Mult_var", "Twopart vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_twopart_Mult_skew", "Twopart vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_twopart_Mult_kurto", "Twopart vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_threepart_Mult_skew", "Threepart vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_threepart_Mult_kurto", "Threepart vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_fourpart_Mult_kurto", "Fourpart vs N_{TPC} ", partHist);

    hist.addClone("Analysis/Charged/", "Analysis/Pion/");
    hist.addClone("Analysis/Charged/", "Analysis/Kaon/");
    hist.addClone("Analysis/Charged/", "Analysis/Proton/");

    // MC Generated
    hist.add("Gen/h_Counts", "Counts", kTH1D, {axisEvents});
    hist.add("Gen/h_VtxZ", "Vertex Z ", kTH1D, {axisVtxZ});
    hist.add("Gen/h_VtxZ_b", "Vertex Z ", kTH1D, {axisVtxZ});
    hist.add("Gen/h_NTPC", "Mid rapidity Multiplicity", kTH1D, {axisMultTPC});
    hist.add("Gen/h_NFT0C", "Forward Multiplicity", kTH1D, {axisMultFT0MMC});
    hist.add("Gen/h2_NTPC_NFT0C", "N_{TPC} vs N_{FT0C}", kTH2D, {{axisMultFT0MMC}, {axisMultTPC}});
    hist.add("Gen/h2_NTPC_NFT0M", "N_{TPC} vs N_{FT0M} Reco", kTH2D, {{axisMultFT0M}, {axisMultTPC}});

    hist.add("Gen/h_NSim", "Truth Multiplicity TPC", kTH1D, {axisMultTPC});
    hist.add("Gen/h2_NTPC_NSim", "Reco vs Truth Multiplicty TPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/h2_NChSim_NSim", "Truth Multiplicty NCh vs NTPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/h2_NTPC_NChSim", "Truth Multiplicty NCh vs Reco NTPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/h2_NTPC_NPiSim", "Truth Multiplicty NPi vs Reco NTPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/h2_NTPC_NKaSim", "Truth Multiplicty NKa vs Reco NTPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/h2_NTPC_NPrSim", "Truth Multiplicty NPr vs Reco NTPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/h2_NFT0C_NFT0CSim", "Reco vs Truth Multplicity FT0C", kTH2D, {{axisMultFT0MMC}, {axisMultFT0M}});

    hist.add("Gen/Charged/h_EtaTruth", "#eta ", kTH1D, {axisEta});
    hist.add("Gen/Charged/h_PhiTruth", "#phi ", kTH1D, {axisPhi});
    hist.add("Gen/Charged/h_PtTruth", "p_{T} ", kTH1D, {axisPt});
    hist.add("Gen/Charged/h2_PtTruth_Eta", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPt}});
    hist.add("Gen/Charged/h2_PtTruth_NFT0M", "p_{T} in Multiplicity Classes", kTH2D, {{axisPt}, {axisMultFT0M}});

    hist.add("Gen/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Gen/Charged/h2_pt_nch", "Truth", kTH2D, {{axisMult}, {axisPt}});
    hist.add("Gen/Charged/h2_pt_nch_prof", "Truth", kTProfile, {axisMult});
    hist.add("Gen/Charged/h_mean_pT", " <p_{T}> ", kTH1D, {axisMeanPt});

    hist.add("Gen/Charged/h_Q1", "Q1", qNMCHist);
    hist.add("Gen/Charged/h_Q2", "Q2", qNMCHist);
    hist.add("Gen/Charged/h_Q3", "Q3", qNMCHist);
    hist.add("Gen/Charged/h_Q4", "Q4", qNMCHist);
    hist.add("Gen/Charged/h_Q1_var", "Q1 vs N_{TPC}", qNMCHist);
    hist.add("Gen/Charged/h_N_var", "N vs N_{TPC}", kTHnSparseD, {axisMultTPC, axisMult, axisMultFT0M});
    hist.add("Gen/Charged/h_twopart_nume_Mult_var", "twopart numerator", kTHnSparseD, {axisMultTPC, axisTpN, axisMultFT0M});
    hist.add("Gen/Charged/h_twopart_deno_Mult_var", "twopart denominator", kTHnSparseD, {axisMultTPC, axisTpD, axisMultFT0M});

    hist.add("Gen/Charged/p_mean_pT_Mult_var", " <p_{T}> ", kTProfile, {axisMultTPC});
    hist.add("Gen/Charged/p_CheckNCh", " 1/denominator vs N_{TPC} ", kTProfile, {axisMultTPC});
    hist.add("Gen/Charged/h_CheckNCh", " 1/denominator vs N_{TPC} ", denoMCHist);
    hist.add("Gen/Charged/h_mean_pT_Mult_var", " <p_{T}> vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_mean_pT_Mult_skew", " <p_{T}> vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_mean_pT_Mult_kurto", " <p_{T}> vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_twopart_Mult_var", "Twopart vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_twopart_Mult_skew", "Twopart vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_twopart_Mult_kurto", "Twopart vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_threepart_Mult_skew", "Threepart vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_threepart_Mult_kurto", "Threepart vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_fourpart_Mult_kurto", "Fourpart vs N_{TPC} ", partMCHist);

    hist.addClone("Gen/Charged/", "Gen/Pion/");

    hist.add("Gen/Pion/h_RapTruth", "y", kTH1D, {axisY});
    hist.add("Gen/Pion/h2_PtTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("Gen/Pion/h2_PtPosTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("Gen/Pion/h2_PtNegTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("Gen/Pion/h_PtPosTruth", "p_{T} (Positive)", kTH1D, {axisPt});
    hist.add("Gen/Pion/h_PtNegTruth", "p_{T} (negative)", kTH1D, {axisPt});
    hist.add("Gen/Pion/h2_PtPosTruth_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("Gen/Pion/h2_PtNegTruth_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});

    hist.addClone("Gen/Pion/", "Gen/Kaon/");
    hist.addClone("Gen/Pion/", "Gen/Proton/");
  }

  enum Mode {
    QA_Pion = 0,
    QA_Kaon,
    QA_Proton,
    Analysis_Charged,
    Analysis_Pion,
    Analysis_Kaon,
    Analysis_Proton,
    Gen_Charged,
    Gen_Pion,
    Gen_Kaon,
    Gen_Proton
  };

  static constexpr std::string_view Dire[] = {
    "QA/Pion/",
    "QA/Kaon/",
    "QA/Proton/",
    "Analysis/Charged/",
    "Analysis/Pion/",
    "Analysis/Kaon/",
    "Analysis/Proton/",
    "Gen/Charged/",
    "Gen/Pion/",
    "Gen/Kaon/",
    "Gen/Proton/"};

  // Event selection cuts:
  template <typename T>
  bool selRun3Col(T const& col)
  {
    hist.fill(HIST("QA/after/h_counts_evSelCuts"), 0);

    if (cfgPosZ) {
      if (std::abs(col.posZ()) > cfgCutPosZ) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), 1);
    }

    if (cfgSel8) {
      if (!col.sel8()) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), 2);
    }
    if (cfgNoSameBunchPileup) {
      if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), 4);
    }

    if (cfgIsVertexITSTPC) {
      if (!col.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), 5);
    }

    return true;
  }

  // Track selection cuts:
  template <typename T>
  bool selTrack(T const& track)
  {
    if (!track.isGlobalTrack())
      return false;

    if (track.pt() < cfgCutPtMin)
      return false;

    if (track.pt() > cfgCutPtMax)
      return false;

    if (track.sign() == 0)
      return false;

    if (std::fabs(track.dcaZ()) > cfgCutDcaZ)
      return false;

    if (std::fabs(track.dcaXY()) > cfgCutDcaXY)
      return false;

    return true;
  }

  // Cuts to reject the tracks
  template <typename T>
  bool rejectTracks(T const& track)
  {
    if (((track.tpcNSigmaEl()) > -3. &&
         (track.tpcNSigmaEl()) < 5.) &&
        (std::fabs(track.tpcNSigmaPi()) > 3 &&
         std::fabs(track.tpcNSigmaKa()) > 3 &&
         std::fabs(track.tpcNSigmaPr()) > 3)) {
      return true;
    }

    return false;
  }

  template <typename T>
  bool selElectrons(T const& track)
  {
    if (std::fabs(track.tpcNSigmaEl()) < cfgCutNSig3) {
      return true;
    }

    return false;
  }

  // PID selction cuts for Low momentum Pions
  template <typename T>
  bool selLowPi(T const& track)
  {
    if (track.pt() >= cfgCutPiPtMin &&
        track.p() <= cfgCutPiThrsldP &&
        std::abs(track.rapidity(MassPiPlus)) < cfgCutRap) {
      if (!track.hasTOF() &&
          std::fabs(track.tpcNSigmaPi()) < cfgCutNSig2) {
        return true;
      }

      if (track.hasTOF() &&
          std::fabs(track.tpcNSigmaPi()) < cfgCutNSig2 &&
          std::fabs(track.tofNSigmaPi()) < cfgCutNSig3) {
        return true;
      }
    }
    return false;
  }

  // PID selction cuts for Low momentum Kaons
  template <typename T>
  bool selLowKa(T const& track)
  {
    if (track.pt() >= cfgCutKaPtMin &&
        track.p() <= cfgCutKaThrsldP &&
        std::abs(track.rapidity(MassKPlus)) < cfgCutRap) {
      if (!track.hasTOF() &&
          std::fabs(track.tpcNSigmaKa()) < cfgCutNSig2) {
        return true;
      }

      if (track.hasTOF() &&
          std::fabs(track.tpcNSigmaKa()) < cfgCutNSig2 &&
          std::fabs(track.tofNSigmaKa()) < cfgCutNSig3) {
        return true;
      }
    }

    return false;
  }

  // PID selction cuts for Low momentum Protons
  template <typename T>
  bool selLowPr(T const& track)
  {
    if (track.pt() >= cfgCutPrPtMin &&
        track.p() <= cfgCutPrThrsldP &&
        std::abs(track.rapidity(MassProton)) < cfgCutRap) {
      if (!track.hasTOF() &&
          std::fabs(track.tpcNSigmaPr()) < cfgCutNSig2) {
        return true;
      }

      if (track.hasTOF() &&
          std::fabs(track.tpcNSigmaPr()) < cfgCutNSig2 &&
          std::fabs(track.tofNSigmaPr()) < cfgCutNSig3) {
        return true;
      }
    }

    return false;
  }

  // PID selction cuts for High momentum Protons
  template <typename T>
  bool selHighPi(T const& track)
  {
    if (track.hasTOF() &&
        track.p() > cfgCutPiThrsldP &&
        std::fabs(track.tpcNSigmaPi()) < cfgCutNSig3 &&
        std::fabs(track.tofNSigmaPi()) < cfgCutNSig3) {

      if (std::abs(track.rapidity(MassPiPlus)) < cfgCutRap) {
        return true;
      }
    }

    return false;
  }

  // PID selction cuts for High momentum Kaons
  template <typename T>
  bool selHighKa(T const& track)
  {
    if (track.hasTOF() &&
        track.p() > cfgCutKaThrsldP &&
        std::fabs(track.tpcNSigmaKa()) < cfgCutNSig3 &&
        ((std::fabs(track.tofNSigmaKa()) < cfgCutNSig3 && track.p() <= cfgCutKaP3) ||
         (std::fabs(track.tofNSigmaKa()) < cfgCutNSig2 && track.p() > cfgCutKaP3))) {

      if (std::abs(track.rapidity(MassKPlus)) < cfgCutRap) {
        return true;
      }
    }

    return false;
  }

  // PID selction cuts for High momentum Protons
  template <typename T>
  bool selHighPr(T const& track)
  {
    if (track.hasTOF() &&
        track.p() > cfgCutPrThrsldP &&
        std::fabs(track.tpcNSigmaPr()) < cfgCutNSig3 &&
        std::fabs(track.tofNSigmaPr()) < cfgCutNSig3) {

      if (std::abs(track.rapidity(MassProton)) < cfgCutRap) {
        return true;
      }
    }

    return false;
  }

  // To find the pT bin
  int findBin(float pT, const std::vector<double>& bins)
  {
    for (size_t i = 0; i < bins.size() - 1; ++i) {
      if (pT >= bins[i] && pT < bins[i + 1]) {
        return i;
      }
    }
    return -1;
  }

  // Find bin index for both pT and eta
  std::pair<int, int> find2DBin(float pT, float rap, const std::vector<double>& ptBins, const std::vector<double>& rapBins)
  {
    int ptBin = -1, rapBin = -1;

    // Find pT bin
    for (size_t i = 0; i < ptBins.size() - 1; ++i) {
      if (pT >= ptBins[i] && pT < ptBins[i + 1]) {
        ptBin = i + 1; // ROOT bins start from 1
        break;
      }
    }

    // Find eta bin
    for (size_t j = 0; j < rapBins.size() - 1; ++j) {
      if (rap >= rapBins[j] && rap < rapBins[j + 1]) {
        rapBin = j + 1;
        break;
      }
    }

    return {ptBin, rapBin};
  }

  // Fill hist before selection cuts:
  template <typename T, typename U>
  void fillBeforeQAHistos(T const& col, U const& tracks)
  {
    for (const auto& track : tracks) {
      hist.fill(HIST("QA/before/h_Eta"), track.eta());
      hist.fill(HIST("QA/before/h_Phi"), track.phi());
      hist.fill(HIST("QA/before/h_Pt"), track.pt());
      hist.fill(HIST("QA/before/h_DcaXY"), track.dcaXY());
      hist.fill(HIST("QA/before/h_DcaZ"), track.dcaZ());
      hist.fill(HIST("QA/before/h2_DcaXY"), track.pt(), track.dcaXY());
      hist.fill(HIST("QA/before/h2_DcaZ"), track.pt(), track.dcaZ());
    }
    hist.fill(HIST("QA/before/h_VtxZ"), col.posZ());
    hist.fill(HIST("QA/before/h_Counts"), 2);
    hist.fill(HIST("QA/before/h_NTPC"), col.multNTracksHasTPC());
    hist.fill(HIST("QA/before/h_Cent"), col.centFT0C());
    hist.fill(HIST("QA/after/h_CentM"), col.centFT0M());
    hist.fill(HIST("QA/before/h_NFT0M"), col.multFT0M());
  }

  // Fill hist after selection cuts:
  template <typename T>
  void fillAfterQAHistos(T const& col)
  {
    hist.fill(HIST("QA/after/h_VtxZ"), col.posZ());
    hist.fill(HIST("QA/after/h_Counts"), 2);
    hist.fill(HIST("QA/after/h_NTPC"), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/h_Cent"), col.centFT0C());
    hist.fill(HIST("QA/after/h_CentM"), col.centFT0M());
    hist.fill(HIST("QA/after/h_NFT0M"), col.multFT0M());
    hist.fill(HIST("QA/after/h_NFT0C"), col.multFT0C());
    hist.fill(HIST("QA/after/h2_NTPC_NFT0M"), col.multFT0M(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/h2_NTPC_NFT0C"), col.multFT0C(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/h2_NTPC_Cent"), col.centFT0C(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/p_NTPC_Cent"), col.centFT0C(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/p_NTPC_NFT0M"), col.multFT0M(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/p_NTPC_NFT0C"), col.multFT0C(), col.multNTracksHasTPC());
  }

  // Fill Charged particles QA:
  template <typename T>
  void fillChargedQAHistos(T const& track, int nFT0M)
  {
    hist.fill(HIST("QA/after/h_Eta"), track.eta());
    hist.fill(HIST("QA/after/h_Phi"), track.phi());
    hist.fill(HIST("QA/after/h_Pt"), track.pt());
    hist.fill(HIST("QA/after/h2_Pt_NFT0M"), track.pt(), nFT0M);
    hist.fill(HIST("QA/after/h2_PvsPinner"), track.p(), track.tpcInnerParam());
    hist.fill(HIST("QA/after/h2_Pt_Eta"), track.eta(), track.pt());
    hist.fill(HIST("QA/after/h_DcaZ"), track.dcaZ());
    hist.fill(HIST("QA/after/h_DcaXY"), track.dcaXY());
    hist.fill(HIST("QA/after/h2_DcaXY"), track.pt(), track.dcaXY());
    hist.fill(HIST("QA/after/h2_DcaZ"), track.pt(), track.dcaZ());

    hist.fill(HIST("QA/after/h_TPCChi2perCluster"), track.tpcChi2NCl());
    hist.fill(HIST("QA/after/h_ITSChi2perCluster"), track.itsChi2NCl());
    hist.fill(HIST("QA/after/h_crossedTPC"), track.tpcNClsCrossedRows());
  }

  // Fill before PID cut QA hist:
  template <typename T>
  void fillBeforePIDQAHistos(T const& track)
  {
    hist.fill(HIST("QA/before/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/before/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/before/h2_pvsm2"), track.mass() * track.mass(), track.p());

    hist.fill(HIST("QA/Pion/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaPi());
    hist.fill(HIST("QA/Pion/before/h2_TOFNsigma"), track.p(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Pion/before/h2_TpcTofNsigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Proton/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaPr());
    hist.fill(HIST("QA/Proton/before/h2_TOFNsigma"), track.p(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Proton/before/h2_TpcTofNsigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Kaon/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaKa());
    hist.fill(HIST("QA/Kaon/before/h2_TOFNsigma"), track.p(), track.tofNSigmaKa());
    hist.fill(HIST("QA/Kaon/before/h2_TpcTofNsigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());

    hist.fill(HIST("QA/before/innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST("QA/before/innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST("QA/before/innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());

    hist.fill(HIST("QA/Pion/innerParam/before/h2_TPCNsigma"), track.tpcInnerParam(), track.tpcNSigmaPi());
    hist.fill(HIST("QA/Pion/innerParam/before/h2_TOFNsigma"), track.tpcInnerParam(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Pion/innerParam/before/h2_TpcTofNsigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Proton/innerParam/before/h2_TPCNsigma"), track.tpcInnerParam(), track.tpcNSigmaPr());
    hist.fill(HIST("QA/Proton/innerParam/before/h2_TOFNsigma"), track.tpcInnerParam(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Proton/innerParam/before/h2_TpcTofNsigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Kaon/innerParam/before/h2_TPCNsigma"), track.tpcInnerParam(), track.tpcNSigmaKa());
    hist.fill(HIST("QA/Kaon/innerParam/before/h2_TOFNsigma"), track.tpcInnerParam(), track.tofNSigmaKa());
    hist.fill(HIST("QA/Kaon/innerParam/before/h2_TpcTofNsigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());
  }

  // Moments Calculation:
  void moments(double pt, double weight, double& Q1, double& Q2, double& Q3, double& Q4)
  {
    Q1 += pt * weight;
    Q2 += pt * pt * weight;
    Q3 += pt * pt * pt * weight;
    Q4 += pt * pt * pt * pt * weight;
  }

  double getCorrectedWeight(const std::vector<double>& ptBins, const std::vector<double>& rapBins,
                            const std::vector<float>& effPtValues, const std::vector<float>& effPtRapValues,
                            double pt, double rap, bool cfgCorrectionPtRap, bool cfgCorrection)
  {
    double weight = 1.0;

    if (cfgCorrectionPtRap) {
      auto [ptBin, rapBin] = find2DBin(pt, rap, ptBins, rapBins);

      if (ptBin != -1 && rapBin != -1) {
        int numPtBins = ptBins.size() - 1; // Number of pt bins
        double efficiency = effPtRapValues[rapBin * numPtBins + ptBin];

        if (efficiency > 0) {
          weight = 1.0 / efficiency;
        }
      }
    } else if (cfgCorrection) {
      int binIndex = findBin(pt, ptBins);

      if (binIndex != -1) {
        double efficiency = effPtValues[binIndex];
        if (efficiency > 0) {
          weight = 1.0 / efficiency;
        }
      }
    }

    return weight;
  }

  // Fill after PID cut QA hist:
  template <int Mode, typename T>
  void fillIdParticleQAHistos(T const& track, const std::vector<double>& ptBins, const std::vector<double>& rapBins, const std::vector<float>& effPtValues, const std::vector<float>& effPtRapValues, double rap, double nSigmaTPC, double nSigmaTOF, int nFT0M, int& N, double& Q1, double& Q2, double& Q3, double& Q4)
  {
    double pt = track.pt();
    double weight = getCorrectedWeight(ptBins, rapBins, effPtValues, effPtRapValues, pt, rap, cfgCorrectionPtRapPID, cfgCorrectionPID);

    N += weight;
    moments(pt, weight, Q1, Q2, Q3, Q4);

    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_Rap_weighted"), rap, pt, weight);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Pt_weighted"), pt, weight);

    hist.fill(HIST(Dire[Mode]) + HIST("h_Pt"), track.pt());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_NFT0M"), track.pt(), nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_Eta"), track.eta(), track.pt());
    if (track.sign() > 0) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtPos"), track.pt());
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtPos_NFT0M"), track.pt(), nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtPos_rap"), rap, track.pt());
    }
    if (track.sign() < 0) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtNeg"), track.pt());
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtNeg_NFT0M"), track.pt(), nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtNeg_rap"), rap, track.pt());
    }

    hist.fill(HIST(Dire[Mode]) + HIST("h_Eta"), track.eta());
    hist.fill(HIST(Dire[Mode]) + HIST("h_Phi"), track.phi());
    hist.fill(HIST(Dire[Mode]) + HIST("h_Rap"), rap);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_Rap"), rap, track.pt());
    hist.fill(HIST(Dire[Mode]) + HIST("h_DcaZ"), track.dcaZ());
    hist.fill(HIST(Dire[Mode]) + HIST("h_DcaXY"), track.dcaXY());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_DcaZ"), track.pt(), track.dcaZ());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_DcaXY"), track.pt(), track.dcaXY());

    hist.fill(HIST(Dire[Mode]) + HIST("h2_TPCNsigma_El"), track.p(), track.tpcNSigmaEl());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TOFNsigma_El"), track.p(), track.tofNSigmaEl());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TPCNsigma"), track.p(), nSigmaTPC);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TOFNsigma"), track.p(), nSigmaTOF);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TpcTofNsigma"), nSigmaTPC, nSigmaTOF);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_pvsm2"), track.mass() * track.mass(), track.p());
    hist.fill(HIST("QA/after/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/after/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/after/h2_pvsm2"), track.mass() * track.mass(), track.p());

    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TPCNsigma_El"), track.tpcInnerParam(), track.tpcNSigmaEl());
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TOFNsigma_El"), track.tpcInnerParam(), track.tofNSigmaEl());
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TPCNsigma"), track.tpcInnerParam(), nSigmaTPC);
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TOFNsigma"), track.tpcInnerParam(), nSigmaTOF);
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TpcTofNsigma"), nSigmaTPC, nSigmaTOF);
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());
    hist.fill(HIST("QA/after/innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST("QA/after/innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST("QA/after/innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());
  }

  template <int Mode>
  void fillPtMCHist(double pt, double eta, double rap, int nFT0M, int pid, int pdgCodePos, int pdgCodeNeg)
  {
    hist.fill(HIST(Dire[Mode]) + HIST("h_PtTruth"), pt);
    hist.fill(HIST(Dire[Mode]) + HIST("h_EtaTruth"), eta);
    hist.fill(HIST(Dire[Mode]) + HIST("h_RapTruth"), rap);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_PtTruth_NFT0M"), pt, nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_PtTruth_Rap"), rap, pt);

    if (pid == pdgCodePos) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtPosTruth"), pt);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtPosTruth_NFT0M"), pt, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtPosTruth_Rap"), rap, pt);
    }
    if (pid == pdgCodeNeg) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtNegTruth"), pt);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtNegTruth_NFT0M"), pt, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtNegTruth_Rap"), rap, pt);
    }
  }

  template <int Mode>
  void fillAnalysisHistos(int nTPC, int nFT0M, int N, double Q1, double Q2, double Q3, double Q4)
  {
    if (N == 0) {
      return;
    }
    double twopart1 = ((Q1 * Q1) - Q2);
    double threepart1 = ((Q1 * Q1 * Q1) - (3 * Q2 * Q1) + 2 * Q3);
    double fourpart1 = ((Q1 * Q1 * Q1 * Q1) - (6 * Q2 * Q1 * Q1) + (3 * Q2 * Q2) + (8 * Q3 * Q1) - 6 * Q4);

    hist.fill(HIST(Dire[Mode]) + HIST("h_Mult"), N);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q1"), nTPC, Q1, nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q2"), nTPC, Q2, nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q3"), nTPC, Q3, nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q4"), nTPC, Q4, nFT0M);

    if (N > 1) {
      double meanPt = Q1 / static_cast<double>(N);
      double nPair = (static_cast<double>(N) * (static_cast<double>(N) - 1));
      double twopart = twopart1 / nPair;
      double checkNDenoVar = (1 / std::sqrt(1 - (1 / static_cast<double>(N))));
      hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT"), meanPt);
      hist.fill(HIST(Dire[Mode]) + HIST("p_mean_pT_Mult_var"), nTPC, meanPt);

      hist.fill(HIST(Dire[Mode]) + HIST("h_Q1_var"), nTPC, Q1, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h_N_var"), nTPC, N, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_nume_Mult_var"), nTPC, twopart1, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_deno_Mult_var"), nTPC, nPair, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT_Mult_var"), nTPC, meanPt, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_Mult_var"), nTPC, twopart, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("p_CheckNCh"), nTPC, checkNDenoVar);
      hist.fill(HIST(Dire[Mode]) + HIST("h_CheckNCh"), nTPC, checkNDenoVar, nFT0M);

      if (N > 2) {
        double nTriplet = (static_cast<double>(N) * (static_cast<double>(N) - 1) * (static_cast<double>(N) - 2));
        double threepart = threepart1 / nTriplet;
        hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT_Mult_skew"), nTPC, meanPt, nFT0M);
        hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_Mult_skew"), nTPC, twopart, nFT0M);
        hist.fill(HIST(Dire[Mode]) + HIST("h_threepart_Mult_skew"), nTPC, threepart, nFT0M);

        if (N > 3) {
          double nQuad = (static_cast<double>(N) * (static_cast<double>(N) - 1) * (static_cast<double>(N) - 2) * (static_cast<double>(N) - 3));
          double fourpart = fourpart1 / nQuad;
          hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT_Mult_kurto"), nTPC, meanPt, nFT0M);
          hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_Mult_kurto"), nTPC, twopart, nFT0M);
          hist.fill(HIST(Dire[Mode]) + HIST("h_threepart_Mult_kurto"), nTPC, threepart, nFT0M);
          hist.fill(HIST(Dire[Mode]) + HIST("h_fourpart_Mult_kurto"), nTPC, fourpart, nFT0M);
        }
      }
    }
  }

  template <bool DataFlag, bool RecoFlag, typename T, typename U>
  void fillHistos(T const& col, U const& tracks)
  {
    int nCh = 0, nTPC = 0, nFT0M = 0, nFT0C = 0;

    int nPi = 0, nKa = 0, nPr = 0;
    double ptCh = 0, q1Ch = 0, q2Ch = 0, q3Ch = 0, q4Ch = 0;
    double ptPi = 0, q1Pi = 0, q2Pi = 0, q3Pi = 0, q4Pi = 0;
    double ptPr = 0, q1Pr = 0, q2Pr = 0, q3Pr = 0, q4Pr = 0;
    double ptKa = 0, q1Ka = 0, q2Ka = 0, q3Ka = 0, q4Ka = 0;

    int nChSim = 0, nSim = 0, nFT0CSim = 0;
    int nPiSim = 0, nKaSim = 0, nPrSim = 0;
    double eta = 0, etaSim = 0, rapSim = 0;
    double ptChSim = 0, q1ChSim = 0, q2ChSim = 0, q3ChSim = 0, q4ChSim = 0;
    double ptPiSim = 0, q1PiSim = 0, q2PiSim = 0, q3PiSim = 0, q4PiSim = 0;
    double ptPrSim = 0, q1PrSim = 0, q2PrSim = 0, q3PrSim = 0, q4PrSim = 0;
    double ptKaSim = 0, q1KaSim = 0, q2KaSim = 0, q3KaSim = 0, q4KaSim = 0;

    array<float, 3> p1, p2;
    double invMassGamma = 0.0;

    for (const auto& [trkEl, trkPos] : soa::combinations(soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      if (trkEl.index() == trkPos.index())
        continue;

      if (!selTrack(trkEl) || !selTrack(trkPos))
        continue;

      if (!selElectrons(trkEl) || !selElectrons(trkPos))
        continue;

      p1 = std::array{trkEl.px(), trkEl.py(), trkEl.pz()};
      p2 = std::array{trkPos.px(), trkPos.py(), trkPos.pz()};

      invMassGamma = RecoDecay::m(std::array{p1, p2}, std::array{MassElectron, MassElectron});
      hist.fill(HIST("QA/after/h_invMass_gamma"), invMassGamma);
    }

    fillAfterQAHistos(col);

    if constexpr (DataFlag) {
      nTPC = col.multNTracksHasTPC();
      nFT0M = col.multFT0M();
      nFT0C = col.multFT0C();

      for (const auto& track : tracks) {
        if (!selTrack(track)) {
          continue;
        }

        double nSigmaTPCPi = track.tpcNSigmaPi();
        double nSigmaTPCKa = track.tpcNSigmaKa();
        double nSigmaTPCPr = track.tpcNSigmaPr();
        double nSigmaTOFPi = track.tofNSigmaPi();
        double nSigmaTOFKa = track.tofNSigmaKa();
        double nSigmaTOFPr = track.tofNSigmaPr();
        double rapPi = track.rapidity(MassPiPlus);
        double rapKa = track.rapidity(MassKPlus);
        double rapPr = track.rapidity(MassProton);

        if (std::fabs(track.eta()) < 0.8) {
          ptCh = track.pt();
          double weight = getCorrectedWeight(ptBins, {}, effValuesCh, {}, ptCh, 0.0, false, cfgCorrection);
          nCh += weight;
          hist.fill(HIST("QA/after/h_Pt_weighted"), ptCh, weight);
          moments(ptCh, weight, q1Ch, q2Ch, q3Ch, q4Ch);

          fillChargedQAHistos(track, nFT0M);
        }

        fillBeforePIDQAHistos(track);

        if (rejectTracks(track)) {
          return;
        }

        if (cfgInvMass == true && invMassGamma < cfgGammaCut) {
          continue;
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowPi(track) == cfgSelLow || selHighPi(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Pion>(track, ptBins, rapBins, effPtValuesPi, effPtRapValuesPi, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, nPi, q1Pi, q2Pi, q3Pi, q4Pi);
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowPi(track) == cfgSelLow && selHighPi(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Pion>(track, ptBins, rapBins, effPtValuesPi, effPtRapValuesPi, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, nPi, q1Pi, q2Pi, q3Pi, q4Pi);
          }
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowKa(track) == cfgSelLow || selHighKa(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Kaon>(track, ptBins, rapBins, effPtValuesKa, effPtRapValuesKa, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, nKa, q1Ka, q2Ka, q3Ka, q4Ka);
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowKa(track) == cfgSelLow && selHighKa(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Kaon>(track, ptBins, rapBins, effPtValuesKa, effPtRapValuesKa, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, nKa, q1Ka, q2Ka, q3Ka, q4Ka);
          }
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowPr(track) == cfgSelLow || selHighPr(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Proton>(track, ptBins, rapBins, effPtValuesPr, effPtRapValuesPr, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, nPr, q1Pr, q2Pr, q3Pr, q4Pr);
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowPr(track) == cfgSelLow && selHighPr(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Proton>(track, ptBins, rapBins, effPtValuesPr, effPtRapValuesPr, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, nPr, q1Pr, q2Pr, q3Pr, q4Pr);
          }
        }
      }
    } else if constexpr (RecoFlag) {
      if (!col.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip...");
        return;
      }
      nTPC = col.multNTracksHasTPC();
      nFT0M = col.multFT0M();
      nFT0C = col.multFT0C();
      hist.fill(HIST("Gen/h_VtxZ"), col.mcCollision().posZ());

      for (const auto& track : tracks) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "No MC Particle for this track, skip...");
          continue;
        }
        auto mcPart = track.mcParticle();
        int pid = mcPart.pdgCode();
        if (!mcPart.isPhysicalPrimary()) {
          continue;
        }

        if (std::abs(track.eta()) < 0.8) {
          nTPC++;
        }

        double nSigmaTPCPi = track.tpcNSigmaPi();
        double nSigmaTPCKa = track.tpcNSigmaKa();
        double nSigmaTPCPr = track.tpcNSigmaPr();
        double nSigmaTOFPi = track.tofNSigmaPi();
        double nSigmaTOFKa = track.tofNSigmaKa();
        double nSigmaTOFPr = track.tofNSigmaPr();
        double rapPi = track.rapidity(MassPiPlus);
        double rapKa = track.rapidity(MassKPlus);
        double rapPr = track.rapidity(MassProton);

        //______________________________Reconstructed Level____________________________________________________//

        if (selTrack(track)) {

          if (std::fabs(track.eta()) < 0.8) {
            ptCh = track.pt();
            double weight = getCorrectedWeight(ptBins, {}, effValuesCh, {}, ptCh, 0.0, false, cfgCorrection);
            hist.fill(HIST("QA/after/h_Pt_weighted"), ptCh, weight);
            nCh += weight;
            moments(ptCh, weight, q1Ch, q2Ch, q3Ch, q4Ch);
            fillChargedQAHistos(track, nFT0M);
          }
          fillBeforePIDQAHistos(track);

          if (cfgRejTrk == true && rejectTracks(track)) {
            return;
          }

          if (cfgInvMass == true && invMassGamma < cfgGammaCut) {
            continue;
          }

          eta = track.eta();
          if (cfgPDGCodeOnly == true) {
            if (std::abs(pid) == kPiPlus && std::abs(rapPi) < 0.5 && track.pt() >= cfgCutPiPtMin) {
              ptPi = track.pt();
              fillIdParticleQAHistos<QA_Pion>(track, ptBins, rapBins, effPtValuesPi, effPtRapValuesPi, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, nPi, q1Pi, q2Pi, q3Pi, q4Pi);
              fillPtMCHist<QA_Pion>(ptPi, eta, rapPi, nFT0M, pid, kPiPlus, kPiMinus);
            }

            if (std::abs(pid) == kKPlus && std::abs(rapKa) < 0.5 && track.pt() >= cfgCutKaPtMin) {
              ptKa = track.pt();
              fillIdParticleQAHistos<QA_Kaon>(track, ptBins, rapBins, effPtValuesKa, effPtRapValuesKa, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, nKa, q1Ka, q2Ka, q3Ka, q4Ka);
              fillPtMCHist<QA_Kaon>(ptKa, eta, rapKa, nFT0M, pid, kKPlus, kKMinus);
            }

            if (std::abs(pid) == kProton && std::abs(rapPr) < 0.5 && track.pt() >= cfgCutPrPtMin) {
              ptPr = track.pt();
              fillIdParticleQAHistos<QA_Proton>(track, ptBins, rapBins, effPtValuesPr, effPtRapValuesPr, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, nPr, q1Pr, q2Pr, q3Pr, q4Pr);
              fillPtMCHist<QA_Proton>(ptPr, eta, rapPr, nFT0M, pid, kProton, kProtonBar);
            }
          }

          if (cfgPidCut == true) {
            if (cfgSelOR == true && cfgSelAND == false) {
              if (selLowPi(track) == cfgSelLow || selHighPi(track) == cfgSelHigh) {
                ptPi = track.pt();
                fillIdParticleQAHistos<QA_Pion>(track, ptBins, rapBins, effPtValuesPi, effPtRapValuesPi, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, nPi, q1Pi, q2Pi, q3Pi, q4Pi);
                if (std::abs(pid) == kPiPlus) {
                  fillPtMCHist<QA_Pion>(ptPi, eta, rapPi, nFT0M, pid, kPiPlus, kPiMinus);
                }
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if (selLowPi(track) == cfgSelLow && selHighPi(track) == cfgSelHigh) {
                ptPi = track.pt();
                fillIdParticleQAHistos<QA_Pion>(track, ptBins, rapBins, effPtValuesPi, effPtRapValuesPi, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, nPi, q1Pi, q2Pi, q3Pi, q4Pi);
                if (std::abs(pid) == kPiPlus) {
                  fillPtMCHist<QA_Pion>(ptPi, eta, rapPi, nFT0M, pid, kPiPlus, kPiMinus);
                }
              }
            }

            if (cfgSelOR == true && cfgSelAND == false) {
              if (selLowKa(track) == cfgSelLow || selHighKa(track) == cfgSelHigh) {
                ptKa = track.pt();
                fillIdParticleQAHistos<QA_Kaon>(track, ptBins, rapBins, effPtValuesKa, effPtRapValuesKa, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, nKa, q1Ka, q2Ka, q3Ka, q4Ka);
                if (std::abs(pid) == kKPlus) {
                  fillPtMCHist<QA_Kaon>(ptKa, eta, rapKa, nFT0M, pid, kKPlus, kKMinus);
                }
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if (selLowKa(track) == cfgSelLow && selHighKa(track) == cfgSelHigh) {
                ptKa = track.pt();
                fillIdParticleQAHistos<QA_Kaon>(track, ptBins, rapBins, effPtValuesKa, effPtRapValuesKa, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, nKa, q1Ka, q2Ka, q3Ka, q4Ka);
                if (std::abs(pid) == kKPlus) {
                  fillPtMCHist<QA_Kaon>(ptKa, eta, rapKa, nFT0M, pid, kKPlus, kKMinus);
                }
              }
            }

            if (cfgSelOR == true && cfgSelAND == false) {
              if (selLowPr(track) == cfgSelLow || selHighPr(track) == cfgSelHigh) {
                ptPr = track.pt();
                fillIdParticleQAHistos<QA_Proton>(track, ptBins, rapBins, effPtValuesPr, effPtRapValuesPr, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, nPr, q1Pr, q2Pr, q3Pr, q4Pr);
                if (std::abs(pid) == kProton) {
                  fillPtMCHist<QA_Proton>(ptPr, eta, rapPr, nFT0M, pid, kProton, kProtonBar);
                }
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if (selLowPr(track) == cfgSelLow && selHighPr(track) == cfgSelHigh) {
                ptPr = track.pt();
                fillIdParticleQAHistos<QA_Proton>(track, ptBins, rapBins, effPtValuesPr, effPtRapValuesPr, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, nPr, q1Pr, q2Pr, q3Pr, q4Pr);
                if (std::abs(pid) == kProton) {
                  fillPtMCHist<QA_Proton>(ptPr, eta, rapPr, nFT0M, pid, kProton, kProtonBar);
                }
              }
            }
          }
        }

        //___________________________________Truth Level____________________________________________________//
        auto charge = 0.;
        auto* pd = pdg->GetParticle(pid);
        if (pd != nullptr) {
          charge = pd->Charge();
        }
        if (std::fabs(charge) < 1e-3) {
          continue;
        }
        if (std::abs(pid) != kElectron && std::abs(pid) != kMuonMinus && std::abs(pid) != kPiPlus && std::abs(pid) != kKPlus && std::abs(pid) != kProton) {
          continue;
        }

        if (std::fabs(mcPart.eta()) < 0.8) {
          nSim++;
        }

        if (mcPart.eta() > -3.3 || mcPart.eta() < -2.1) {
          nFT0CSim++;
        }

        if (mcPart.pt() > cfgCutPtMin && mcPart.pt() < cfgCutPtMax) {

          if (std::abs(mcPart.eta()) < 0.8) {
            nChSim++;
            ptChSim = mcPart.pt();
            moments(ptChSim, 1.0, q1ChSim, q2ChSim, q3ChSim, q4ChSim);
            hist.fill(HIST("Gen/Charged/h_PtTruth"), mcPart.pt());
            hist.fill(HIST("Gen/Charged/h2_PtTruth_NFT0M"), mcPart.pt(), nFT0M);
            hist.fill(HIST("Gen/Charged/h2_PtTruth_Eta"), mcPart.eta(), mcPart.pt());
            hist.fill(HIST("Gen/Charged/h_EtaTruth"), mcPart.eta());
            hist.fill(HIST("Gen/Charged/h_PhiTruth"), mcPart.phi());
          }

          if (std::abs(mcPart.y()) > cfgCutRap) {
            continue;
          }

          if (std::abs(pid) == kPiPlus && mcPart.pt() >= cfgCutPiPtMin) {
            etaSim = mcPart.eta();
            rapSim = mcPart.y();

            if (cfgSelOR == true && cfgSelAND == false) {
              if (mcPart.p() <= cfgCutPiThrsldP || mcPart.p() > cfgCutPiThrsldP) {
                nPiSim++;
                ptPiSim = mcPart.pt();
                moments(ptPiSim, 1.0, q1PiSim, q2PiSim, q3PiSim, q4PiSim);
                fillPtMCHist<Gen_Pion>(ptPiSim, etaSim, rapSim, nFT0M, pid, kPiPlus, kPiMinus);
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutPiThrsldP) && (cfgSelHigh == true && mcPart.p() > cfgCutPiThrsldP)) {
                nPiSim++;
                ptPiSim = mcPart.pt();
                moments(ptPiSim, 1.0, q1PiSim, q2PiSim, q3PiSim, q4PiSim);
                fillPtMCHist<Gen_Pion>(ptPiSim, etaSim, rapSim, nFT0M, pid, kPiPlus, kPiMinus);
              }
            }
            hist.fill(HIST("Gen/Pion/h_PhiTruth"), mcPart.phi());
          }

          if (std::abs(pid) == kKPlus && mcPart.pt() >= cfgCutKaPtMin) {
            if (cfgSelOR == true && cfgSelAND == false) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutPiThrsldP) || (cfgSelHigh == true && mcPart.p() > cfgCutPiThrsldP)) {
                nKaSim++;
                ptKaSim = mcPart.pt();
                moments(ptKaSim, 1.0, q1KaSim, q2KaSim, q3KaSim, q4KaSim);
                fillPtMCHist<Gen_Kaon>(ptKaSim, etaSim, rapSim, nFT0M, pid, kKPlus, kKMinus);
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutKaThrsldP) && (cfgSelHigh == true && mcPart.p() > cfgCutKaThrsldP)) {
                nKaSim++;
                ptKaSim = mcPart.pt();
                moments(ptKaSim, 1.0, q1KaSim, q2KaSim, q3KaSim, q4KaSim);
                fillPtMCHist<Gen_Kaon>(ptKaSim, etaSim, rapSim, nFT0M, pid, kKPlus, kKMinus);
              }
            }
            hist.fill(HIST("Gen/Kaon/h_PhiTruth"), mcPart.phi());
          }

          if (std::abs(pid) == kProton && mcPart.pt() >= cfgCutPrPtMin) {
            if (cfgSelOR == true && cfgSelAND == false) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutPrThrsldP) || (cfgSelHigh == true && mcPart.p() > cfgCutPrThrsldP)) {
                nPrSim++;
                ptPrSim = mcPart.pt();
                moments(ptPrSim, 1.0, q1PrSim, q2PrSim, q3PrSim, q4PrSim);
                fillPtMCHist<Gen_Proton>(ptPrSim, etaSim, rapSim, nFT0M, pid, kProton, kProtonBar);
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutPrThrsldP) && (cfgSelHigh == true && mcPart.p() > cfgCutPrThrsldP)) {
                nPrSim++;
                ptPrSim = mcPart.pt();
                moments(ptPrSim, 1.0, q1PrSim, q2PrSim, q3PrSim, q4PrSim);
                fillPtMCHist<Gen_Proton>(ptPrSim, etaSim, rapSim, nFT0M, pid, kProton, kProtonBar);
              }
            }
            hist.fill(HIST("Gen/Proton/h_PhiTruth"), mcPart.phi());
          }
        }
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////

      for (const auto& track : tracks) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "No MC Particle for this track, skip...");
          continue;
        }
        auto mcPart = track.mcParticle();
        int pid = mcPart.pdgCode();
        if (!mcPart.isPhysicalPrimary()) {
          continue;
        }

        double rapPi = track.rapidity(MassPiPlus);
        double rapKa = track.rapidity(MassKPlus);
        double rapPr = track.rapidity(MassProton);

        if (selTrack(track)) {
          if (std::abs(track.eta()) < 0.8) {
            double pt = track.pt();
            double weight = getCorrectedWeight(ptBins, {}, effValuesCh, {}, pt, 0.0, false, cfgCorrection);
            hist.fill(HIST("QA/after/h2_pt_nch"), nCh, pt, weight);
            hist.fill(HIST("QA/after/h2_pt_nch_prof"), nCh, pt, weight);
          }

          if (selLowPi(track) == cfgSelLow || selHighPi(track) == cfgSelHigh) {
            ptPi = track.pt();
            double weight = getCorrectedWeight(ptBins, rapBins, effPtValuesPi, effPtRapValuesPi, ptPi, rapPi, cfgCorrectionPtRapPID, cfgCorrectionPID);
            hist.fill(HIST("QA/Pion/h2_pt_nch"), nPi, ptPi, weight);
            hist.fill(HIST("QA/Pion/h2_pt_nch_prof"), nPi, ptPi, weight);
          }

          if (selLowKa(track) == cfgSelLow || selHighKa(track) == cfgSelHigh) {
            ptKa = track.pt();
            double weight = getCorrectedWeight(ptBins, rapBins, effPtValuesKa, effPtRapValuesKa, ptKa, rapKa, cfgCorrectionPtRapPID, cfgCorrectionPID);
            hist.fill(HIST("QA/Kaon/h2_pt_nch"), nKa, ptKa, weight);
            hist.fill(HIST("QA/Kaon/h2_pt_nch_prof"), nKa, ptKa, weight);
          }

          if (selLowPr(track) == cfgSelLow || selHighPr(track) == cfgSelHigh) {
            ptPr = track.pt();
            double weight = getCorrectedWeight(ptBins, rapBins, effPtValuesPr, effPtRapValuesPr, ptPr, rapPr, cfgCorrectionPtRapPID, cfgCorrectionPID);
            hist.fill(HIST("QA/Proton/h2_pt_nch"), nPr, ptPr, weight);
            hist.fill(HIST("QA/Proton/h2_pt_nch_prof"), nPr, ptPr, weight);
          }
        }

        auto charge = 0.;
        auto* pd = pdg->GetParticle(pid);
        if (pd != nullptr) {
          charge = pd->Charge();
        }
        if (std::fabs(charge) < 1e-3) {
          continue;
        }
        if (std::abs(pid) != kElectron && std::abs(pid) != kMuonMinus && std::abs(pid) != kPiPlus && std::abs(pid) != kKPlus && std::abs(pid) != kProton) {
          continue;
        }
        if (mcPart.pt() > cfgCutPtMin && mcPart.pt() < cfgCutPtMax) {
          if (std::abs(mcPart.eta()) < 0.8) {
            double pt = mcPart.pt();
            hist.fill(HIST("Gen/Charged/h2_pt_nch"), nChSim, pt);
            hist.fill(HIST("Gen/Charged/h2_pt_nch_prof"), nChSim, pt);
          }

          if (std::abs(mcPart.y()) > 0.5)
            continue;

          if (std::abs(pid) == kPiPlus && mcPart.pt() >= cfgCutPiPtMin) {
            hist.fill(HIST("Gen/Pion/h2_pt_nch"), nPiSim, mcPart.pt());
            hist.fill(HIST("Gen/Pion/h2_pt_nch_prof"), nPiSim, mcPart.pt());
          }
          if (std::abs(pid) == kKPlus && mcPart.pt() >= cfgCutKaPtMin) {
            hist.fill(HIST("Gen/Kaon/h2_pt_nch"), nKaSim, mcPart.pt());
            hist.fill(HIST("Gen/Kaon/h2_pt_nch_prof"), nKaSim, mcPart.pt());
          }
          if (std::abs(pid) == kProton && mcPart.pt() >= cfgCutPrPtMin) {
            hist.fill(HIST("Gen/Proton/h2_pt_nch"), nPrSim, mcPart.pt());
            hist.fill(HIST("Gen/Proton/h2_pt_nch_prof"), nPrSim, mcPart.pt());
          }
        }
      }

      hist.fill(HIST("Gen/h_Counts"), 2);
      hist.fill(HIST("QA/after/h_VtxZReco"), col.posZ());

      if (nSim > 0)
        hist.fill(HIST("Gen/h_NSim"), nSim);

      if (nSim > 0 && nChSim > 0)
        hist.fill(HIST("Gen/h2_NChSim_NSim"), nSim, nChSim);

      if (nSim > 0 && nTPC > 0)
        hist.fill(HIST("Gen/h2_NTPC_NSim"), nSim, nTPC);

      if (nChSim > 0 && nTPC > 0)
        hist.fill(HIST("Gen/h2_NTPC_NChSim"), nTPC, nChSim);

      if (nPiSim > 0 && nTPC > 0)
        hist.fill(HIST("Gen/h2_NTPC_NPiSim"), nTPC, nPiSim);

      if (nKaSim > 0 && nTPC > 0)
        hist.fill(HIST("Gen/h2_NTPC_NKaSim"), nTPC, nKaSim);

      if (nPrSim > 0 && nTPC > 0)
        hist.fill(HIST("Gen/h2_NTPC_NPrSim"), nTPC, nPrSim);

      hist.fill(HIST("Gen/h_NTPC"), nTPC);
      hist.fill(HIST("Gen/h_NFT0C"), nFT0CSim);
      hist.fill(HIST("Gen/h2_NTPC_NFT0C"), nFT0CSim, nTPC);
      hist.fill(HIST("Gen/h2_NTPC_NFT0M"), nFT0M, nTPC);

      if (nFT0C != 0 && nFT0CSim != 0)
        hist.fill(HIST("Gen/h2_NFT0C_NFT0CSim"), nFT0CSim, nFT0C);

      fillAnalysisHistos<Gen_Charged>(nTPC, nFT0M, nChSim, q1ChSim, q2ChSim, q3ChSim, q4ChSim);
      fillAnalysisHistos<Gen_Pion>(nTPC, nFT0M, nPiSim, q1PiSim, q2PiSim, q3PiSim, q4PiSim);
      fillAnalysisHistos<Gen_Kaon>(nTPC, nFT0M, nKaSim, q1KaSim, q2KaSim, q3KaSim, q4KaSim);
      fillAnalysisHistos<Gen_Proton>(nTPC, nFT0M, nPrSim, q1PrSim, q2PrSim, q3PrSim, q4PrSim);
    }

    if (nTPC > 0 && nCh > 0)
      hist.fill(HIST("QA/after/h2_NTPC_NCh"), nTPC, nCh);

    if (nPi > 0 && nTPC > 0)
      hist.fill(HIST("QA/after/h2_NTPC_NPi"), nTPC, nPi);

    if (nKa > 0 && nTPC > 0)
      hist.fill(HIST("QA/after/h2_NTPC_NKa"), nTPC, nKa);

    if (nPr > 0 && nTPC > 0)
      hist.fill(HIST("QA/after/h2_NTPC_NPr"), nTPC, nPr);

    fillAnalysisHistos<Analysis_Charged>(nTPC, nFT0M, nCh, q1Ch, q2Ch, q3Ch, q4Ch);
    fillAnalysisHistos<Analysis_Pion>(nTPC, nFT0M, nPi, q1Pi, q2Pi, q3Pi, q4Pi);
    fillAnalysisHistos<Analysis_Kaon>(nTPC, nFT0M, nKa, q1Ka, q2Ka, q3Ka, q4Ka);
    fillAnalysisHistos<Analysis_Proton>(nTPC, nFT0M, nPr, q1Pr, q2Pr, q3Pr, q4Pr);
  }

  void processRun3(MyRun3Collisions::iterator const& col, MyAllTracks const& tracks)
  {
    // Before Collision and Track Cuts:
    fillBeforeQAHistos(col, tracks);

    // After Collision and Track Cuts:
    if (selRun3Col(col)) {
      fillHistos<true, false>(col, tracks);
    }
  }
  PROCESS_SWITCH(MeanPtFlucId, processRun3, "Process for Run-3", false);

  void processMCRecoSimRun3(MyRun3MCCollisions::iterator const& col, aod::McCollisions const&, MyMCTracks const& tracks, aod::McParticles const&)
  {
    // Before Collision and Track Cuts:
    fillBeforeQAHistos(col, tracks);

    hist.fill(HIST("Gen/h_VtxZ_b"), col.mcCollision().posZ());

    // After Collision and Track Cuts:
    if (selRun3Col(col)) {
      fillHistos<false, true>(col, tracks);
    }
  }
  PROCESS_SWITCH(MeanPtFlucId, processMCRecoSimRun3, "process MC Reconstructed & Truth Run-3", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MeanPtFlucId>(cfgc)};
}
