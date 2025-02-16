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
  Configurable<std::vector<float>> effPtValuesPi{"effPtValuesPi", {0, 0.408075, 0.473332, 0.48221, 0.469699, 0.472676, 0.482403, 0.478351, 0.38468, 0.249696, 0.244316, 0.235498, 0.236493, 0.241719, 0.245363, 0.248324, 0.251595, 0.254327, 0.257727, 0.260208, 0.263414, 0.267699, 0.270322, 0.275128, 0.280835, 0.284328, 0.288791, 0.294786, 0.292418, 0.299766, 0.299413, 0.301257, 0.305466, 0.304929, 0.316837, 0.317915, 0.316018}, "effeciency values for Pions"};
  Configurable<std::vector<float>> effPtRapValuesPi{"effPtRapValuesPi", {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5520, 0.5331, 0.5201, 0.5140, 0.5135, 0.5271, 0.4378, 0.2341, 0.2668, 0.2606, 0.2484, 0.2553, 0.2533, 0.2606, 0.2621, 0.2645, 0.2651, 0.2785, 0.2767, 0.2867, 0.2927, 0.2862, 0.3046, 0.3050, 0.3078, 0.3025, 0.3156, 0.2991, 0.3137, 0.3093, 0.3210, 0.3300, 0.3207, 0.3356, 0.3495, 0.3475, 0.0000, 0.5487, 0.5299, 0.5228, 0.5152, 0.5165, 0.5260, 0.5123, 0.2426, 0.2604, 0.2523, 0.2404, 0.2459, 0.2472, 0.2482, 0.2519, 0.2541, 0.2563, 0.2577, 0.2620, 0.2598, 0.2636, 0.2791, 0.2770, 0.2866, 0.2798, 0.2887, 0.2969, 0.3003, 0.2947, 0.2927, 0.2917, 0.3183, 0.3067, 0.3099, 0.3234, 0.3285, 0.0000, 0.5350, 0.5276, 0.5192, 0.5150, 0.5209, 0.5296, 0.5339, 0.2891, 0.2645, 0.2591, 0.2535, 0.2577, 0.2597, 0.2652, 0.2697, 0.2695, 0.2754, 0.2892, 0.2815, 0.2849, 0.2885, 0.2946, 0.3080, 0.3094, 0.3208, 0.3174, 0.3211, 0.3239, 0.3231, 0.3332, 0.3319, 0.3381, 0.3369, 0.3325, 0.3442, 0.3663, 0.0000, 0.4808, 0.5197, 0.5194, 0.5159, 0.5217, 0.5304, 0.5366, 0.3708, 0.2952, 0.2832, 0.2743, 0.2762, 0.2824, 0.2870, 0.2895, 0.2874, 0.2922, 0.2953, 0.2996, 0.3024, 0.3084, 0.3159, 0.3138, 0.3305, 0.3309, 0.3366, 0.3340, 0.3404, 0.3292, 0.3403, 0.3279, 0.3486, 0.3729, 0.3447, 0.3603, 0.3566, 0.0000, 0.4303, 0.5043, 0.5166, 0.5174, 0.5203, 0.5282, 0.5341, 0.4151, 0.2868, 0.2790, 0.2662, 0.2682, 0.2718, 0.2771, 0.2761, 0.2790, 0.2824, 0.2820, 0.2873, 0.2869, 0.2835, 0.2888, 0.3004, 0.3002, 0.3143, 0.3159, 0.3056, 0.3083, 0.3197, 0.3315, 0.3146, 0.3280, 0.3207, 0.3138, 0.3467, 0.3454, 0.0000, 0.4125, 0.4794, 0.5124, 0.5108, 0.5186, 0.5265, 0.5342, 0.4518, 0.2585, 0.2503, 0.2408, 0.2399, 0.2449, 0.2506, 0.2588, 0.2568, 0.2575, 0.2597, 0.2636, 0.2693, 0.2664, 0.2763, 0.2773, 0.2808, 0.2902, 0.2899, 0.2872, 0.2968, 0.3056, 0.2918, 0.3007, 0.3048, 0.3077, 0.3041, 0.2981, 0.3222, 0.0000, 0.4042, 0.4703, 0.5100, 0.5087, 0.5131, 0.5238, 0.5280, 0.4796, 0.2243, 0.2254, 0.2157, 0.2245, 0.2249, 0.2374, 0.2325, 0.2372, 0.2463, 0.2433, 0.2444, 0.2506, 0.2525, 0.2526, 0.2661, 0.2605, 0.2831, 0.2841, 0.2878, 0.2874, 0.2964, 0.2937, 0.3050, 0.3100, 0.2976, 0.3073, 0.3083, 0.3122, 0.0000, 0.4036, 0.4928, 0.5121, 0.4988, 0.5041, 0.5162, 0.5238, 0.5000, 0.2092, 0.2117, 0.1977, 0.2032, 0.2100, 0.2097, 0.2180, 0.2140, 0.2196, 0.2170, 0.2174, 0.2269, 0.2259, 0.2318, 0.2318, 0.2415, 0.2447, 0.2391, 0.2604, 0.2516, 0.2619, 0.2569, 0.2600, 0.2573, 0.2637, 0.2717, 0.2740, 0.2681, 0.0000, 0.3452, 0.4427, 0.4595, 0.4499, 0.4631, 0.4782, 0.4818, 0.4794, 0.1810, 0.1763, 0.1678, 0.1683, 0.1725, 0.1760, 0.1765, 0.1788, 0.1778, 0.1792, 0.1834, 0.1799, 0.1872, 0.1882, 0.1866, 0.1915, 0.1917, 0.1996, 0.2100, 0.2097, 0.2089, 0.2145, 0.2168, 0.2099, 0.2127, 0.2073, 0.2149, 0.2060, 0.0000, 0.2717, 0.3686, 0.4011, 0.3970, 0.4118, 0.4221, 0.4316, 0.4357, 0.1215, 0.1172, 0.1120, 0.1131, 0.1175, 0.1202, 0.1243, 0.1230, 0.1247, 0.1276, 0.1210, 0.1237, 0.1244, 0.1311, 0.1279, 0.1349, 0.1413, 0.1449, 0.1452, 0.1355, 0.1478, 0.1471, 0.1443, 0.1558, 0.1395, 0.1451, 0.1489, 0.1403, 0.0000, 0.2701, 0.3727, 0.4013, 0.4013, 0.4164, 0.4267, 0.4355, 0.4404, 0.1342, 0.1344, 0.1293, 0.1253, 0.1289, 0.1308, 0.1340, 0.1330, 0.1395, 0.1402, 0.1358, 0.1421, 0.1385, 0.1481, 0.1470, 0.1513, 0.1497, 0.1598, 0.1516, 0.1558, 0.1605, 0.1635, 0.1601, 0.1610, 0.1599, 0.1601, 0.1727, 0.1671, 0.0000, 0.3255, 0.4427, 0.4688, 0.4548, 0.4695, 0.4819, 0.4880, 0.4857, 0.1827, 0.1819, 0.1737, 0.1721, 0.1753, 0.1782, 0.1812, 0.1846, 0.1798, 0.1907, 0.1844, 0.1904, 0.1916, 0.1968, 0.1888, 0.2013, 0.2075, 0.2080, 0.1968, 0.2124, 0.2114, 0.2100, 0.2169, 0.2154, 0.2192, 0.2189, 0.2242, 0.2097, 0.0000, 0.3559, 0.4877, 0.5037, 0.4938, 0.4998, 0.5092, 0.5187, 0.4925, 0.1932, 0.1911, 0.1847, 0.1877, 0.1948, 0.1948, 0.1996, 0.2005, 0.2017, 0.2029, 0.2009, 0.2049, 0.2076, 0.2177, 0.2211, 0.2158, 0.2359, 0.2317, 0.2399, 0.2390, 0.2345, 0.2436, 0.2378, 0.2438, 0.2376, 0.2495, 0.2477, 0.2544, 0.0000, 0.3632, 0.4649, 0.5039, 0.4997, 0.5064, 0.5164, 0.5270, 0.4722, 0.2336, 0.2317, 0.2199, 0.2248, 0.2287, 0.2334, 0.2378, 0.2367, 0.2405, 0.2390, 0.2434, 0.2463, 0.2497, 0.2564, 0.2580, 0.2573, 0.2702, 0.2775, 0.2865, 0.2808, 0.2899, 0.2966, 0.2892, 0.2826, 0.2860, 0.2987, 0.2915, 0.3188, 0.0000, 0.3912, 0.4740, 0.5051, 0.5038, 0.5093, 0.5205, 0.5300, 0.4531, 0.2901, 0.2830, 0.2754, 0.2772, 0.2822, 0.2859, 0.2954, 0.2940, 0.2955, 0.2977, 0.2976, 0.3012, 0.3065, 0.3059, 0.3208, 0.3244, 0.3326, 0.3296, 0.3297, 0.3454, 0.3402, 0.3265, 0.3414, 0.3474, 0.3438, 0.3496, 0.3428, 0.3895, 0.0000, 0.4220, 0.4939, 0.5100, 0.5090, 0.5090, 0.5183, 0.5215, 0.4323, 0.3189, 0.3098, 0.2970, 0.3009, 0.3036, 0.3102, 0.3106, 0.3132, 0.3145, 0.3170, 0.3242, 0.3323, 0.3320, 0.3308, 0.3345, 0.3471, 0.3511, 0.3423, 0.3487, 0.3591, 0.3691, 0.3523, 0.3801, 0.3492, 0.3678, 0.3793, 0.3818, 0.3993, 0.0000, 0.4704, 0.5064, 0.5109, 0.5062, 0.5106, 0.5219, 0.5250, 0.3847, 0.3256, 0.3213, 0.3061, 0.3100, 0.3152, 0.3253, 0.3235, 0.3294, 0.3316, 0.3300, 0.3438, 0.3404, 0.3513, 0.3530, 0.3572, 0.3688, 0.3709, 0.3792, 0.3676, 0.3775, 0.3744, 0.3786, 0.3887, 0.3953, 0.4057, 0.3939, 0.4099, 0.4465, 0.0000, 0.5253, 0.5212, 0.5124, 0.5030, 0.5088, 0.5185, 0.5247, 0.3098, 0.2995, 0.2949, 0.2843, 0.2866, 0.2991, 0.3002, 0.3030, 0.3063, 0.3106, 0.3176, 0.3178, 0.3224, 0.3323, 0.3324, 0.3369, 0.3503, 0.3554, 0.3458, 0.3615, 0.3689, 0.3563, 0.3619, 0.3656, 0.3689, 0.3717, 0.3991, 0.3880, 0.4147, 0.0000, 0.5410, 0.5214, 0.5112, 0.5068, 0.5073, 0.5173, 0.5006, 0.2571, 0.2815, 0.2729, 0.2671, 0.2652, 0.2720, 0.2758, 0.2746, 0.2820, 0.2832, 0.2870, 0.2906, 0.2963, 0.2961, 0.2997, 0.3028, 0.3165, 0.3137, 0.3139, 0.3317, 0.3318, 0.3241, 0.3317, 0.3324, 0.3305, 0.3340, 0.3501, 0.3671, 0.3516, 0.0000, 0.5411, 0.5255, 0.5099, 0.5021, 0.5054, 0.5161, 0.4317, 0.2497, 0.2792, 0.2713, 0.2659, 0.2699, 0.2748, 0.2805, 0.2782, 0.2890, 0.2846, 0.2897, 0.2942, 0.3005, 0.3039, 0.3034, 0.3167, 0.3165, 0.3181, 0.3328, 0.3217, 0.3308, 0.3305, 0.3281, 0.3323, 0.3493, 0.3511, 0.3402, 0.3450, 0.3665, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, "pT rap effeciency values for Pions"};
  Configurable<std::vector<float>> effPtValuesKa{"effPtValuesKa", {0, 0, 0, 0.312144, 0.369847, 0.38878, 0.413275, 0.393619, 0.315429, 0.1375, 0.146659, 0.147163, 0.155197, 0.163588, 0.168412, 0.177936, 0.17782, 0.186872, 0.190744, 0.199436, 0.197739, 0.192307, 0.198484, 0.19927, 0.218019, 0.221942, 0.237642, 0.235765, 0.249873, 0.251034, 0.259014, 0.268821, 0.275786, 0.280998, 0.29936, 0.304559, 0.312684}, "pT eta effeciency values for Kaons"};
  Configurable<std::vector<float>> effPtRapValuesKa{"effPtRapValuesKa", {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1711, 0.4322, 0.4572, 0.4474, 0.1392, 0.1648, 0.1836, 0.1891, 0.1732, 0.1794, 0.1778, 0.1914, 0.2058, 0.2122, 0.2190, 0.2180, 0.2156, 0.2448, 0.2679, 0.2626, 0.2984, 0.2953, 0.3137, 0.3228, 0.3696, 0.3786, 0.3517, 0.4040, 0.3901, 0.4604, 0.5460, 0.5458, 0.4860, 0.5755, 0.0000, 0.0000, 0.0000, 0.4031, 0.4326, 0.4559, 0.4946, 0.2264, 0.1289, 0.1516, 0.1633, 0.1679, 0.1906, 0.2048, 0.2132, 0.2183, 0.2247, 0.2651, 0.2572, 0.2406, 0.2365, 0.2701, 0.2525, 0.2643, 0.2806, 0.3142, 0.3317, 0.3364, 0.3298, 0.3212, 0.3564, 0.3873, 0.3982, 0.4500, 0.4643, 0.5220, 0.6104, 0.0000, 0.0000, 0.0000, 0.3993, 0.4172, 0.4474, 0.4974, 0.3636, 0.1374, 0.1723, 0.1809, 0.1743, 0.1815, 0.1985, 0.2034, 0.2150, 0.2187, 0.2393, 0.2346, 0.2479, 0.2409, 0.2599, 0.2530, 0.2897, 0.2897, 0.2922, 0.3214, 0.3364, 0.3855, 0.4070, 0.4153, 0.3828, 0.4269, 0.4225, 0.5101, 0.5487, 0.5465, 0.0000, 0.0000, 0.0000, 0.4002, 0.4271, 0.4676, 0.4789, 0.4990, 0.1447, 0.1659, 0.1784, 0.1734, 0.1930, 0.2012, 0.2124, 0.2340, 0.2522, 0.2673, 0.2594, 0.2609, 0.2565, 0.2824, 0.2783, 0.2727, 0.2953, 0.3379, 0.3473, 0.3603, 0.3743, 0.4059, 0.3903, 0.4152, 0.4254, 0.4361, 0.4638, 0.5754, 0.5599, 0.0000, 0.0000, 0.0000, 0.4020, 0.4132, 0.4547, 0.4768, 0.5111, 0.2441, 0.1760, 0.1872, 0.1925, 0.2073, 0.2169, 0.2318, 0.2329, 0.2453, 0.2519, 0.2538, 0.2599, 0.2565, 0.2559, 0.2720, 0.2844, 0.2968, 0.3048, 0.3230, 0.3233, 0.3402, 0.3664, 0.3795, 0.4017, 0.3661, 0.3992, 0.4338, 0.4848, 0.5206, 0.0000, 0.0000, 0.0000, 0.3770, 0.4073, 0.4435, 0.4825, 0.5040, 0.3451, 0.1814, 0.1806, 0.1758, 0.1899, 0.1964, 0.1949, 0.2105, 0.2144, 0.2261, 0.2461, 0.2451, 0.2293, 0.2176, 0.2297, 0.2395, 0.2615, 0.2716, 0.2832, 0.2936, 0.3194, 0.3466, 0.3317, 0.3598, 0.3355, 0.4052, 0.3900, 0.4260, 0.4737, 0.0000, 0.0000, 0.0000, 0.3688, 0.4227, 0.4245, 0.4550, 0.4903, 0.4048, 0.1576, 0.1597, 0.1740, 0.1776, 0.1911, 0.2010, 0.2130, 0.1963, 0.2356, 0.2175, 0.2320, 0.2517, 0.2136, 0.2220, 0.2380, 0.2635, 0.2890, 0.2978, 0.3152, 0.3344, 0.3124, 0.3439, 0.3770, 0.3441, 0.3936, 0.3892, 0.4360, 0.4348, 0.0000, 0.0000, 0.0000, 0.3585, 0.4101, 0.4185, 0.4567, 0.4896, 0.4543, 0.1239, 0.1382, 0.1451, 0.1523, 0.1614, 0.1722, 0.1873, 0.1906, 0.2008, 0.2057, 0.2276, 0.2153, 0.1922, 0.2014, 0.2141, 0.2305, 0.2502, 0.2477, 0.2940, 0.2962, 0.2987, 0.3237, 0.3023, 0.3161, 0.3414, 0.3568, 0.3466, 0.3992, 0.0000, 0.0000, 0.0000, 0.3535, 0.3858, 0.3940, 0.4202, 0.4460, 0.4512, 0.1152, 0.1319, 0.1280, 0.1366, 0.1450, 0.1519, 0.1673, 0.1597, 0.1797, 0.1727, 0.1865, 0.1909, 0.1746, 0.1674, 0.1788, 0.1923, 0.2043, 0.2303, 0.2264, 0.2317, 0.2345, 0.2547, 0.2597, 0.2948, 0.3002, 0.2719, 0.3049, 0.2998, 0.0000, 0.0000, 0.0000, 0.2898, 0.3194, 0.3320, 0.3707, 0.3963, 0.4181, 0.0748, 0.0840, 0.0799, 0.0835, 0.0968, 0.0985, 0.1080, 0.1054, 0.1080, 0.1286, 0.1336, 0.1204, 0.1228, 0.1232, 0.1261, 0.1495, 0.1487, 0.1535, 0.1564, 0.1631, 0.1768, 0.1956, 0.1909, 0.1957, 0.1873, 0.2248, 0.2199, 0.2283, 0.0000, 0.0000, 0.0000, 0.1413, 0.1521, 0.1386, 0.1263, 0.1148, 0.1107, 0.0205, 0.0212, 0.0217, 0.0234, 0.0245, 0.0235, 0.0265, 0.0263, 0.0262, 0.0308, 0.0299, 0.0331, 0.0280, 0.0296, 0.0291, 0.0286, 0.0337, 0.0370, 0.0326, 0.0355, 0.0382, 0.0408, 0.0384, 0.0408, 0.0433, 0.0463, 0.0443, 0.0464, 0.0000, 0.0000, 0.0000, 0.3353, 0.3807, 0.3944, 0.4307, 0.4725, 0.4790, 0.1181, 0.1269, 0.1314, 0.1350, 0.1546, 0.1514, 0.1554, 0.1507, 0.1604, 0.1682, 0.1766, 0.1862, 0.1731, 0.1818, 0.1945, 0.2096, 0.2152, 0.2253, 0.2221, 0.2372, 0.2559, 0.2340, 0.2753, 0.2901, 0.2818, 0.2684, 0.2700, 0.3042, 0.0000, 0.0000, 0.0000, 0.3494, 0.4016, 0.4101, 0.4368, 0.4755, 0.4433, 0.1166, 0.1291, 0.1342, 0.1404, 0.1585, 0.1586, 0.1631, 0.1759, 0.1788, 0.1769, 0.2053, 0.2112, 0.1885, 0.1921, 0.2021, 0.2006, 0.2387, 0.2434, 0.2547, 0.2705, 0.2690, 0.2966, 0.2698, 0.2877, 0.3023, 0.3092, 0.3123, 0.3627, 0.0000, 0.0000, 0.0000, 0.3515, 0.4035, 0.4132, 0.4501, 0.4755, 0.3929, 0.1655, 0.1738, 0.1867, 0.1915, 0.2132, 0.2122, 0.2105, 0.2303, 0.2302, 0.2343, 0.2467, 0.2437, 0.2364, 0.2472, 0.2688, 0.2485, 0.2821, 0.3039, 0.3322, 0.3020, 0.3388, 0.3728, 0.3620, 0.3690, 0.3572, 0.4131, 0.4549, 0.4871, 0.0000, 0.0000, 0.0000, 0.3803, 0.4282, 0.4350, 0.4591, 0.4858, 0.3540, 0.1968, 0.2057, 0.2057, 0.2208, 0.2226, 0.2314, 0.2284, 0.2425, 0.2542, 0.2500, 0.2618, 0.2562, 0.2512, 0.2674, 0.2827, 0.3082, 0.2932, 0.3296, 0.3233, 0.3689, 0.3794, 0.4202, 0.4337, 0.4032, 0.4574, 0.4255, 0.4487, 0.6078, 0.0000, 0.0000, 0.0000, 0.3940, 0.4184, 0.4363, 0.4688, 0.5184, 0.2523, 0.1979, 0.2038, 0.2260, 0.2290, 0.2466, 0.2631, 0.2637, 0.2640, 0.2802, 0.3055, 0.3082, 0.2831, 0.2731, 0.2934, 0.2987, 0.3306, 0.3460, 0.3697, 0.3769, 0.3792, 0.4081, 0.4194, 0.4079, 0.4766, 0.4226, 0.4606, 0.5413, 0.5926, 0.0000, 0.0000, 0.0000, 0.4037, 0.4264, 0.4461, 0.4638, 0.4895, 0.1648, 0.1881, 0.1862, 0.1997, 0.2273, 0.2403, 0.2347, 0.2588, 0.2769, 0.2900, 0.3019, 0.3046, 0.2762, 0.2977, 0.3163, 0.3498, 0.3600, 0.3686, 0.3807, 0.4021, 0.4237, 0.4305, 0.4928, 0.4549, 0.4803, 0.4745, 0.5351, 0.6006, 0.6772, 0.0000, 0.0000, 0.0000, 0.4067, 0.4231, 0.4615, 0.4870, 0.3601, 0.1493, 0.1847, 0.1981, 0.2003, 0.1911, 0.2023, 0.2315, 0.2318, 0.2385, 0.2464, 0.2636, 0.2626, 0.2690, 0.2777, 0.2707, 0.3373, 0.2916, 0.3545, 0.3536, 0.4019, 0.4201, 0.3934, 0.4489, 0.4565, 0.4178, 0.5096, 0.5585, 0.5518, 0.6223, 0.0000, 0.0000, 0.0000, 0.4070, 0.4335, 0.4495, 0.4745, 0.2224, 0.1328, 0.1608, 0.1702, 0.1827, 0.1939, 0.2113, 0.2331, 0.2482, 0.2555, 0.2676, 0.2594, 0.2693, 0.2592, 0.2608, 0.2601, 0.3022, 0.3137, 0.3099, 0.3338, 0.3411, 0.3938, 0.3477, 0.3983, 0.3982, 0.3937, 0.4450, 0.4986, 0.5582, 0.5724, 0.0000, 0.0000, 0.0000, 0.1654, 0.4371, 0.4530, 0.4372, 0.1515, 0.1719, 0.2036, 0.1941, 0.1823, 0.1778, 0.1908, 0.1862, 0.2106, 0.2292, 0.2369, 0.2251, 0.2314, 0.2476, 0.2559, 0.2732, 0.2931, 0.3197, 0.3373, 0.3648, 0.3785, 0.3693, 0.4063, 0.4223, 0.4156, 0.4548, 0.5356, 0.5663, 0.5616, 0.5979, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, "pT rap effeciency values for Kaons"};
  Configurable<std::vector<float>> effPtValuesPr{"effPtValuesPr", {0, 0, 0, 0, 0, 0, 0, 0.394712, 0.425251, 0.458426, 0.489121, 0.509505, 0.516103, 0.517117, 0.491584, 0.450721, 0.379836, 0.253402, 0.257575, 0.261382, 0.260373, 0.269008, 0.266811, 0.265011, 0.272768, 0.269553, 0.276003, 0.279878, 0.284216, 0.276346, 0.293437, 0.294727, 0.281017, 0.287609, 0.292402, 0.28614, 0.307208}, "effeciency values for Kaons"};
  Configurable<std::vector<float>> effPtRapValuesPr{"effPtRapValuesPr", {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0345, 0.3209, 0.6345, 0.7076, 0.6968, 0.3865, 0.3562, 0.3741, 0.3760, 0.3825, 0.3930, 0.3762, 0.3869, 0.3708, 0.3414, 0.3246, 0.3439, 0.3257, 0.3233, 0.3185, 0.3693, 0.3426, 0.3601, 0.3260, 0.3715, 0.3373, 0.4214, 0.3766, 0.3692, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.2474, 0.5755, 0.6398, 0.6618, 0.6986, 0.7487, 0.6339, 0.3570, 0.3286, 0.3126, 0.3066, 0.2918, 0.2906, 0.3151, 0.3293, 0.3585, 0.3321, 0.3348, 0.3798, 0.3705, 0.3326, 0.3756, 0.4221, 0.3602, 0.3791, 0.4167, 0.4121, 0.3948, 0.4191, 0.3476, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5490, 0.6266, 0.6415, 0.6566, 0.6976, 0.7392, 0.7155, 0.4702, 0.3043, 0.3034, 0.3037, 0.3458, 0.3393, 0.3188, 0.3545, 0.3541, 0.3383, 0.3254, 0.3477, 0.3959, 0.3584, 0.3581, 0.3751, 0.3855, 0.3655, 0.4074, 0.3771, 0.3188, 0.3816, 0.3437, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5701, 0.6201, 0.6612, 0.6533, 0.6684, 0.6916, 0.7126, 0.7208, 0.3838, 0.3316, 0.3169, 0.3159, 0.3281, 0.3361, 0.3345, 0.3542, 0.3535, 0.3339, 0.3527, 0.3803, 0.3599, 0.3745, 0.3810, 0.4013, 0.4036, 0.4118, 0.4016, 0.4161, 0.4036, 0.4198, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5860, 0.6117, 0.6459, 0.6697, 0.6637, 0.7082, 0.7381, 0.7367, 0.5968, 0.3374, 0.3051, 0.3492, 0.3423, 0.3498, 0.3690, 0.3573, 0.4151, 0.3900, 0.3929, 0.4092, 0.3775, 0.4070, 0.3861, 0.4141, 0.3865, 0.4036, 0.3960, 0.4424, 0.3997, 0.3739, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5673, 0.6094, 0.6451, 0.6705, 0.6427, 0.6821, 0.7109, 0.7315, 0.7319, 0.4056, 0.3433, 0.3316, 0.3208, 0.3931, 0.3721, 0.3669, 0.3688, 0.3418, 0.3440, 0.3496, 0.3414, 0.3877, 0.3804, 0.3131, 0.3899, 0.3649, 0.3495, 0.4038, 0.3501, 0.3558, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5497, 0.6281, 0.6455, 0.6752, 0.6818, 0.6935, 0.7480, 0.7510, 0.7712, 0.5599, 0.3500, 0.3410, 0.3289, 0.3403, 0.3424, 0.3240, 0.3404, 0.3540, 0.3268, 0.3331, 0.3636, 0.3601, 0.3469, 0.3593, 0.3637, 0.3669, 0.4036, 0.3625, 0.3350, 0.3712, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5518, 0.5855, 0.6229, 0.6740, 0.6720, 0.6969, 0.7505, 0.7625, 0.7747, 0.6585, 0.2959, 0.2823, 0.2930, 0.2922, 0.3070, 0.3106, 0.3171, 0.3162, 0.3181, 0.3020, 0.3095, 0.3160, 0.3324, 0.3606, 0.3356, 0.3447, 0.3284, 0.3208, 0.3748, 0.3126, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5482, 0.5687, 0.5996, 0.6375, 0.6556, 0.6755, 0.6978, 0.7299, 0.7756, 0.7260, 0.2638, 0.2749, 0.2854, 0.2523, 0.2516, 0.2650, 0.2604, 0.2451, 0.3007, 0.2926, 0.2598, 0.2695, 0.2873, 0.2964, 0.2977, 0.2946, 0.3014, 0.2974, 0.3122, 0.3164, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.4730, 0.5139, 0.5155, 0.5724, 0.6327, 0.6003, 0.6216, 0.6769, 0.7379, 0.6791, 0.1884, 0.1699, 0.1829, 0.1780, 0.1624, 0.1892, 0.1804, 0.1976, 0.1994, 0.1923, 0.1981, 0.1878, 0.1947, 0.2073, 0.2188, 0.2078, 0.1938, 0.1953, 0.2260, 0.1745, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1458, 0.1096, 0.1024, 0.0974, 0.0973, 0.0919, 0.0947, 0.1018, 0.1024, 0.1042, 0.0298, 0.0291, 0.0299, 0.0289, 0.0281, 0.0325, 0.0319, 0.0316, 0.0329, 0.0310, 0.0377, 0.0349, 0.0370, 0.0328, 0.0397, 0.0414, 0.0379, 0.0395, 0.0388, 0.0340, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5397, 0.5801, 0.5897, 0.6112, 0.6864, 0.6721, 0.6859, 0.7086, 0.7534, 0.7322, 0.2720, 0.2496, 0.2315, 0.2494, 0.2812, 0.2877, 0.2433, 0.2465, 0.2750, 0.2883, 0.2848, 0.2793, 0.2909, 0.2771, 0.2888, 0.2923, 0.2572, 0.2723, 0.2758, 0.2731, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5495, 0.5804, 0.6113, 0.6502, 0.6647, 0.6813, 0.7079, 0.7474, 0.7657, 0.6306, 0.2714, 0.2598, 0.2601, 0.2452, 0.2792, 0.2720, 0.2845, 0.2724, 0.2774, 0.2813, 0.2949, 0.2674, 0.2935, 0.2860, 0.3022, 0.2930, 0.2968, 0.3176, 0.3171, 0.2917, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5416, 0.5600, 0.6354, 0.6564, 0.6746, 0.6809, 0.6946, 0.7371, 0.7607, 0.5364, 0.3708, 0.3458, 0.3653, 0.3599, 0.3769, 0.3640, 0.3667, 0.3653, 0.3829, 0.3929, 0.3713, 0.3772, 0.3674, 0.3810, 0.4024, 0.3838, 0.4130, 0.4039, 0.3830, 0.4046, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5512, 0.5967, 0.6277, 0.6578, 0.6424, 0.6641, 0.7317, 0.7343, 0.7142, 0.4351, 0.3927, 0.4098, 0.3937, 0.4046, 0.3986, 0.4117, 0.4140, 0.4294, 0.4020, 0.4036, 0.3858, 0.3582, 0.4257, 0.3741, 0.3968, 0.4209, 0.4387, 0.4221, 0.4059, 0.4374, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5739, 0.6044, 0.6561, 0.6594, 0.6280, 0.6595, 0.6962, 0.7008, 0.5354, 0.3991, 0.3788, 0.3769, 0.3824, 0.4046, 0.4187, 0.4477, 0.4040, 0.4276, 0.4801, 0.4211, 0.4550, 0.4688, 0.4087, 0.4302, 0.4011, 0.4649, 0.4776, 0.4433, 0.4286, 0.5134, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.6076, 0.6062, 0.6302, 0.6758, 0.6457, 0.6936, 0.7342, 0.7175, 0.3698, 0.3439, 0.3499, 0.3405, 0.3785, 0.3427, 0.3453, 0.3833, 0.3716, 0.3945, 0.4290, 0.4172, 0.4177, 0.3970, 0.4051, 0.4350, 0.4395, 0.4208, 0.4605, 0.4489, 0.4315, 0.4499, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5830, 0.6400, 0.6661, 0.6551, 0.6424, 0.6899, 0.7067, 0.4596, 0.2994, 0.3129, 0.3194, 0.3341, 0.3655, 0.3699, 0.3448, 0.3879, 0.3518, 0.3476, 0.3879, 0.3961, 0.3959, 0.4016, 0.3773, 0.3704, 0.4236, 0.4153, 0.3890, 0.4157, 0.4405, 0.3765, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.2265, 0.5939, 0.6490, 0.6726, 0.6857, 0.7228, 0.6359, 0.3889, 0.3641, 0.3183, 0.2852, 0.2965, 0.2959, 0.2933, 0.3224, 0.3247, 0.3324, 0.3821, 0.3729, 0.3990, 0.3941, 0.3753, 0.3909, 0.3909, 0.3953, 0.3842, 0.4036, 0.4685, 0.4010, 0.4241, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0364, 0.3236, 0.6224, 0.6982, 0.6916, 0.4223, 0.3984, 0.3900, 0.4060, 0.3757, 0.3848, 0.3947, 0.3534, 0.3863, 0.3527, 0.3452, 0.3373, 0.3268, 0.3212, 0.3690, 0.3220, 0.3078, 0.3798, 0.3162, 0.3735, 0.3842, 0.4035, 0.4351, 0.4374, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, "pT rap effeciency values for Protons"};

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
    hist.add("QA/after/h_vtxZReco", "Simulated Vertex Z", kTH1D, {axisVtxZ});

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

    hist.add("QA/Pion/h_rap", "y ", kTH1D, {axisY});
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

    hist.add("QA/Pion/h2_Pt_rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtPos_rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtNeg_rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtPosTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtNegTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_Pt_Eta", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPt}});
    hist.add("QA/Pion/h2_DcaZ", "DCA_{z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/Pion/h2_DcaXY", "DCA_{xy}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/Pion/h2_Pt_rap_weighted", "p_{T} vs #eta weighted", kTH2D, {{axisY}, {axisPt}});
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
    hist.add("Gen/h_vtxZ", "Vertex Z ", kTH1D, {axisVtxZ});
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
          std::fabs(track.tpcNSigmaPi()) < cfgCutNSig3 &&
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
          std::fabs(track.tpcNSigmaKa()) < cfgCutNSig3 &&
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
          std::fabs(track.tpcNSigmaPr()) < cfgCutNSig3 &&
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

    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_rap_weighted"), rap, pt, weight);
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
    hist.fill(HIST(Dire[Mode]) + HIST("h_rap"), rap);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_rap"), rap, track.pt());
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
            hist.fill(HIST("Gen/Charged/h2_PtTruth_Eta"), mcPart.pt(), mcPart.eta());
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
      hist.fill(HIST("QA/after/h_vtxZReco"), col.posZ());
      hist.fill(HIST("Gen/h_vtxZ"), col.mcCollision().posZ());

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
